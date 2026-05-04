"""
Unit tests for ImageProperties.

This module verifies EPL+shear image-property solving in
``ler.image_properties.ImageProperties`` (njit and lenstronomy backends), and
the underlying caustic / source-sampling / lens-equation njit routines used by
the analytical solver. Lens-parameter draws use
``LensGalaxyParameterDistribution.epl_shear_sl_parameters_rvs``; the
fixed-lens validation tests mirror ``examples/image_properties/``.

Test Coverage:
--------------
- Class initialization sanity checks (attributes, Da interpolators, njit solver)
- Wrong ``lens_model_list`` for the njit backend raises ``ValueError``
- Custom astropy cosmology propagates into Da interpolators
- Image-property output sanity for both backends
  (``image_properties_epl_shear_njit`` vs ``image_properties_epl_shear_lenstronomy``)
  on the same ``epl_shear_sl_parameters_rvs`` draw
  (image positions are NOT compared row-wise: source positions are sampled
  independently inside each backend)
- CAUSTIC VALIDATION: ler ``caustic_points_epl_shear`` vs lenstronomy
  ``caustics_epl_shear`` (boundary distance + convex-hull area)
- SOURCE SAMPLING VALIDATION: ``sample_source_from_double_caustic`` returns
  finite source positions, all inside the LER caustic polygon
- LENS EQUATION SOLVING VALIDATION: ler ``image_position_analytical_njit`` vs
  lenstronomy ``LensEquationSolver.image_position_from_source`` on the SAME
  source positions (image counts, positions, magnifications)
- ``recover_redundant_parameters`` and ``produce_effective_params`` on one njit
  dict (shared setup; keys from ``tests_utils``).
- Optional ``slow`` test: njit vs lenstronomy backend wall time on the same
  lens batch.
"""

import copy

import numpy as np
import pytest
from astropy.cosmology import LambdaCDM
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
from lenstronomy.LensModel.Solver.epl_shear_solver import (
    caustics_epl_shear as caustics_lenstronomy,
)

from tests_utils import (
    CommonTestUtils,
    median_call_time,
    EXPECTED_EFFECTIVE_KEYS,
    EXPECTED_IMAGE_KEYS,
    EXPECTED_REDUNDANT_KEYS,
    EXPECTED_SOURCE_POS_KEYS,
)
from ler.image_properties import (
    ImageProperties,
    caustic_points_epl_shear,
    sample_source_from_double_caustic,
    image_position_analytical_njit,
)
from ler.image_properties.cross_section_njit import phi_q2_ellipticity
from ler.lens_galaxy_population import LensGalaxyParameterDistribution


# ---------------------------------------------------------------------------
# Test configuration
# ---------------------------------------------------------------------------

# Number of lens-parameter samples drawn per backend test
N_SAMPLES = 10

# Light proposal count for cross-section-based EPL+shear sampling (fast tests)
N_PROP_FAST = 50

DEFAULT_LENS_CONFIG = dict(
    npool=6,
    z_min=0.0,
    z_max=10.0,
    create_new_interpolator=False,
    verbose=False,
    buffer_size=200,
)


# ---------------------------------------------------------------------------
# Small helpers to keep tests readable
# ---------------------------------------------------------------------------

def _make_config(interpolator_dir, **overrides):
    # Start from defaults and allow per-test overrides.
    config = DEFAULT_LENS_CONFIG.copy()
    config["directory"] = interpolator_dir
    config.update(overrides)
    return config


@pytest.fixture(scope="module")
def lens_epl_joint(interpolator_dir):
    """
    Module-scoped LensGalaxyParameterDistribution (with ImageProperties mixed in).

    Initialization is expensive (Da interpolators + njit EPL+shear solver compile);
    one instance is reused across all tests in this module. The cross-section-based
    sampler is reconfigured to a small ``n_prop`` for speed.
    """
    cfg = _make_config(interpolator_dir)
    lens = LensGalaxyParameterDistribution(lens_type="epl_shear_galaxy", **cfg)
    lens.lens_functions["cross_section_based_sampler"] = "importance_sampler_partial"
    lens.lens_functions_params["cross_section_based_sampler"]["n_prop"] = N_PROP_FAST
    lens.cross_section_based_sampler = lens._initialization_cross_section_sampler()
    return lens


# Reference fixed lens used by the caustic / sampling / lens-equation tests
# (matches examples/image_properties/lens_equation_EPL_SHEAR.ipynb).
_FIXED_LENS = dict(
    theta_E=1.0,
    q=0.6,
    phi=np.pi / 4,
    gamma=1.8,
    gamma1=-0.05,
    gamma2=-0.05,
    zl=0.8,
    zs=2.0,
)


def _kwargs_lens_lenstronomy(p):
    """Build a lenstronomy ``kwargs_lens`` for an EPL+shear lens dict."""
    e1, e2 = phi_q2_ellipticity(p["phi"], p["q"])
    return [
        dict(
            theta_E=p["theta_E"], e1=float(e1), e2=float(e2),
            gamma=p["gamma"], center_x=0.0, center_y=0.0,
        ),
        dict(gamma1=p["gamma1"], gamma2=p["gamma2"], ra_0=0.0, dec_0=0.0),
    ]


class TestImageProperties(CommonTestUtils):
    """Tests for ImageProperties initialization, EPL+shear image solving and helpers."""

    def test_init(self, lens_epl_joint):
        """
        Tests
        -----
        - LensGalaxyParameterDistribution + ImageProperties initializes.
        - lens_type, lens_model_list, npool, n_max_images, cosmo are set correctly.
        - Da interpolators are present (set by OpticalDepth before ImageProperties).
        - njit EPL+shear solver ``epl_solver`` is constructed.
        """
        lens = lens_epl_joint
        assert lens.lens_type == "epl_shear_galaxy", f"lens_type: expected 'epl_shear_galaxy', got {lens.lens_type}"
        assert lens.lens_model_list == ["EPL_NUMBA", "SHEAR"], \
            f"lens_model_list: expected ['EPL_NUMBA', 'SHEAR'], got {lens.lens_model_list}"
        assert lens.npool == DEFAULT_LENS_CONFIG["npool"], \
            f"npool: expected {DEFAULT_LENS_CONFIG['npool']}, got {lens.npool}"
        assert lens.n_max_images >= 2, f"n_max_images must be >= 2, got {lens.n_max_images}"
        assert lens.cosmo is not None, "cosmo must not be None after initialization"
        assert hasattr(lens, "angular_diameter_distance"), "angular_diameter_distance interpolator is missing"
        assert hasattr(lens, "angular_diameter_distance_z1z2"), "angular_diameter_distance_z1z2 interpolator is missing"
        assert lens.epl_solver is not None, "epl_solver must not be None after initialization"

    def test_standalone_image_properties_init(self, interpolator_dir):
        """
        Tests
        -----
        - Standalone ImageProperties (no OpticalDepth) builds Da interpolators
          and stores ``directory`` / ``z_min`` / ``z_max``.
        - ``epl_solver`` is constructed for the njit backend.
        - Wrong ``lens_model_list`` for the njit backend raises ``ValueError``.
        """
        ip = ImageProperties(
            npool=DEFAULT_LENS_CONFIG["npool"],
            directory=interpolator_dir,
            z_min=DEFAULT_LENS_CONFIG["z_min"],
            z_max=DEFAULT_LENS_CONFIG["z_max"],
            create_new_interpolator=DEFAULT_LENS_CONFIG["create_new_interpolator"],
            multiprocessing_verbose=False,
        )
        assert ip.directory == interpolator_dir, \
            f"directory: expected {interpolator_dir}, got {ip.directory}"
        assert ip.z_min == DEFAULT_LENS_CONFIG["z_min"], \
            f"z_min: expected {DEFAULT_LENS_CONFIG['z_min']}, got {ip.z_min}"
        assert ip.z_max == DEFAULT_LENS_CONFIG["z_max"], \
            f"z_max: expected {DEFAULT_LENS_CONFIG['z_max']}, got {ip.z_max}"
        assert hasattr(ip, "angular_diameter_distance"), "angular_diameter_distance interpolator is missing"
        assert hasattr(ip, "angular_diameter_distance_z1z2"), "angular_diameter_distance_z1z2 interpolator is missing"
        assert ip.epl_solver is not None, "epl_solver must not be None after initialization"

        with pytest.raises(ValueError, match="lens_model_list"):
            ImageProperties(
                lens_model_list=["EPL_NUMBA"],
                image_properties_function="image_properties_epl_shear_njit",
                directory=interpolator_dir,
                create_new_interpolator=DEFAULT_LENS_CONFIG["create_new_interpolator"],
                multiprocessing_verbose=False,
            )

    def test_custom_cosmology_propagates_to_da(self, interpolator_dir):
        """
        Tests
        -----
        - User-supplied astropy cosmology is stored on ``ImageProperties.cosmo``.
        - The Da interpolators agree with the same cosmology evaluated via astropy.
        """
        cosmo = LambdaCDM(H0=67.4, Om0=0.315, Ode0=0.685, Tcmb0=2.725)
        ip = ImageProperties(
            npool=DEFAULT_LENS_CONFIG["npool"],
            cosmology=cosmo,
            directory=interpolator_dir,
            z_min=DEFAULT_LENS_CONFIG["z_min"],
            z_max=DEFAULT_LENS_CONFIG["z_max"],
            # Force a fresh interpolator so it uses the supplied cosmology.
            create_new_interpolator=True,
            multiprocessing_verbose=False,
        )
        assert ip.cosmo is cosmo, \
            f"cosmo: expected the same astropy cosmology object, got {type(ip.cosmo)}"

        zs = np.array([0.5, 1.0, 2.0])
        zl = np.array([0.2, 0.5, 1.0])
        Da_zs = ip.angular_diameter_distance.function(zs)
        Da_zl_zs = ip.angular_diameter_distance_z1z2.function(zl, zs)
        np.testing.assert_allclose(
            Da_zs, cosmo.angular_diameter_distance(zs).value, rtol=2e-3,
        )
        np.testing.assert_allclose(
            Da_zl_zs,
            cosmo.angular_diameter_distance_z1z2(zl, zs).value,
            rtol=2e-3,
        )

    def test_image_properties_njit_vs_lenstronomy_output_sanity(self, lens_epl_joint):
        """
        Tests
        -----
        - Both ``image_properties_epl_shear_njit`` and
          ``image_properties_epl_shear_lenstronomy`` run on the SAME
          ``epl_shear_sl_parameters_rvs`` draw and return the expected keys.
        - Image-property arrays have shape ``(N_SAMPLES, n_max_images)``.
        - Each event has at least 2 images, finite image data in the first
          ``n_images`` slots, and time delays sorted ascendingly.

        Notes
        -----
        Source positions are sampled independently inside each backend, so the
        per-image positions / magnifications / time delays are NOT compared
        row-wise across backends.
        """
        lens = lens_epl_joint
        # Single base draw: each backend gets its own copy because
        # the methods may delete `theta_E` from the input dict.
        lp = lens.epl_shear_sl_parameters_rvs(size=N_SAMPLES)

        ip_ls = ImageProperties(
            npool=DEFAULT_LENS_CONFIG["npool"],
            image_properties_function="image_properties_epl_shear_lenstronomy",
            multiprocessing_verbose=False,
        )

        res_njit = lens.image_properties_epl_shear_njit(dict(lp))
        res_ls = ip_ls.image_properties_epl_shear_lenstronomy(dict(lp))

        for label, res, n_max in [
            ("njit", res_njit, lens.n_max_images),
            ("lenstronomy", res_ls, ip_ls.n_max_images),
        ]:
            shape2 = (N_SAMPLES, n_max)
            self._assert_param_dict_valid(
                res,
                EXPECTED_IMAGE_KEYS,
                size=shape2,
                nan_to_num=True,
            )
            self._assert_param_dict_valid(
                res,
                EXPECTED_SOURCE_POS_KEYS,
                size=N_SAMPLES,
            )
            assert res["x0_image_positions"].shape == (N_SAMPLES, n_max), \
                f"{label}: x0_image_positions shape expected ({N_SAMPLES}, {n_max}), got {res['x0_image_positions'].shape}"
            assert res["n_images"].shape == (N_SAMPLES,), \
                f"{label}: n_images shape expected ({N_SAMPLES},), got {res['n_images'].shape}"

            for i in range(N_SAMPLES):
                n_img = int(res["n_images"][i])
                assert n_img >= 2, f"{label} row {i}: expected n_images>=2, got {n_img}"
                td = res["time_delays"][i, :n_img]
                mu = res["magnifications"][i, :n_img]
                x0 = res["x0_image_positions"][i, :n_img]
                x1 = res["x1_image_positions"][i, :n_img]
                self._assert_array_valid(td, name=f"{label} time_delays[{i}]", size=n_img)
                self._assert_array_valid(mu, name=f"{label} magnifications[{i}]", size=n_img)
                self._assert_array_valid(x0, name=f"{label} x0[{i}]", size=n_img)
                self._assert_array_valid(x1, name=f"{label} x1[{i}]", size=n_img)
                assert np.all(np.diff(td) >= 0), (
                    f"{label} row {i}: time delays not sorted ascending"
                )

    def test_caustic_njit_vs_lenstronomy(self):
        """
        Tests
        -----
        - LER ``caustic_points_epl_shear`` (njit) and lenstronomy
          ``caustics_epl_shear`` produce equivalent double-caustic boundaries:
          - both polylines are finite and have the expected shape
          - each LER point lies within ``boundary_threshold`` of the lenstronomy
            curve and vice versa (max nearest-neighbor distance)
          - convex-hull areas agree to within 5%
        """
        from scipy.spatial import ConvexHull
        from scipy.spatial.distance import cdist

        p = _FIXED_LENS
        num_th = 500

        pts_ler = caustic_points_epl_shear(
            p["theta_E"], p["q"], p["phi"], p["gamma"], p["gamma1"], p["gamma2"],
            num_th=num_th, maginf=-100.0,
        )
        pts_ls = caustics_lenstronomy(
            _kwargs_lens_lenstronomy(p),
            return_which="double", maginf=-100.0,
            num_th=num_th, sourceplane=True,
        )

        assert pts_ler.shape == (2, num_th), \
            f"LER caustic shape expected (2, {num_th}), got {pts_ler.shape}"
        assert np.isfinite(pts_ler).all(), "LER caustic points contain non-finite values"
        assert pts_ls.shape[0] == 2 and pts_ls.shape[1] >= 3, \
            f"lenstronomy caustic shape expected (2, >=3), got {pts_ls.shape}"
        assert np.isfinite(pts_ls).all(), "lenstronomy caustic points contain non-finite values"

        ler_T = np.column_stack((pts_ler[0], pts_ler[1]))
        ls_T = np.column_stack((pts_ls[0], pts_ls[1]))
        max_d_ler = float(np.max(np.min(cdist(ler_T, ls_T), axis=1)))
        max_d_ls = float(np.max(np.min(cdist(ls_T, ler_T), axis=1)))

        boundary_threshold = 1e-3  # arcsec-scale (theta_E = 1)
        assert max_d_ler < boundary_threshold, (
            f"max LER->LS boundary distance {max_d_ler:.2e} > {boundary_threshold:.2e}"
        )
        assert max_d_ls < boundary_threshold, (
            f"max LS->LER boundary distance {max_d_ls:.2e} > {boundary_threshold:.2e}"
        )

        area_ler = float(ConvexHull(ler_T).volume)
        area_ls = float(ConvexHull(ls_T).volume)
        np.testing.assert_allclose(area_ler / area_ls, 1.0, atol=0.05)

    def test_source_sampling_inside_caustic(self):
        """
        Tests
        -----
        - ``sample_source_from_double_caustic`` returns finite source coordinates.
        - At least 95% of sampled sources lie inside the LER caustic polygon
          (the small slack covers points on or just outside the discretized boundary).
        """
        from matplotlib.path import Path

        p = _FIXED_LENS
        size = 100

        pts_ler = caustic_points_epl_shear(
            p["theta_E"], p["q"], p["phi"], p["gamma"], p["gamma1"], p["gamma2"],
            num_th=500, maginf=-100.0,
        )

        np.random.seed(0)
        xs = np.empty(size)
        ys = np.empty(size)
        for i in range(size):
            xs[i], ys[i] = sample_source_from_double_caustic(
                p["theta_E"], p["q"], p["phi"], p["gamma"], p["gamma1"], p["gamma2"],
                num_th=500, maginf=-100.0,
            )
        self._assert_array_valid(xs, name="xs", size=size, finite=True)
        self._assert_array_valid(ys, name="ys", size=size, finite=True)

        polygon = Path(np.column_stack((pts_ler[0], pts_ler[1])))
        # Small radius accounts for boundary-resolution noise.
        inside = polygon.contains_points(np.column_stack((xs, ys)), radius=1e-4)
        frac = float(np.mean(inside))
        assert frac >= 0.95, f"only {100 * frac:.1f}% of sources inside caustic"

    def test_lens_equation_njit_vs_lenstronomy(self):
        """
        Tests
        -----
        - For a fixed lens and source positions sampled inside its double caustic,
          ler ``image_position_analytical_njit`` and lenstronomy
          ``LensEquationSolver.image_position_from_source`` find the SAME images:
          - same image count (>= 2) for at least 90% of sources
          - matched image positions agree within 0.02 (theta_E units)
          - matched magnifications agree within 5% on average
        """
        p = _FIXED_LENS
        size = 20

        cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        ls_model = LensModel(
            lens_model_list=["EPL_NUMBA", "SHEAR"],
            z_lens=p["zl"], z_source=p["zs"], cosmo=cosmo,
        )
        ls_solver = LensEquationSolver(ls_model)
        kwargs_lens = _kwargs_lens_lenstronomy(p)

        np.random.seed(0)
        diff_pos = []
        diff_mu_rel = []
        n_match = 0
        n_total = 0
        for _ in range(size):
            x_src, y_src = sample_source_from_double_caustic(
                p["theta_E"], p["q"], p["phi"], p["gamma"], p["gamma1"], p["gamma2"],
                num_th=500, maginf=-100.0,
            )
            if not (np.isfinite(x_src) and np.isfinite(y_src)):
                continue

            x_ler, y_ler, _, mu_ler, _, n_ler = image_position_analytical_njit(
                x_src=x_src, y_src=y_src,
                q=p["q"], phi=p["phi"], gamma=p["gamma"],
                gamma1=p["gamma1"], gamma2=p["gamma2"], theta_E=p["theta_E"],
            )
            x_ls, y_ls = ls_solver.image_position_from_source(
                sourcePos_x=x_src, sourcePos_y=y_src,
                kwargs_lens=kwargs_lens,
                solver="analytical", magnification_limit=0.01,
                arrival_time_sort=True,
            )
            n_ls = len(x_ls)
            n_total += 1

            if n_ler == n_ls and n_ler >= 2:
                n_match += 1
                # Greedy nearest-neighbor pairing (LER -> lenstronomy).
                for i in range(n_ler):
                    j = int(np.argmin(
                        np.hypot(x_ler[i] - x_ls, y_ler[i] - y_ls)
                    ))
                    diff_pos.append(float(np.hypot(x_ler[i] - x_ls[j], y_ler[i] - y_ls[j])))
                    mu_j = float(ls_model.magnification(x_ls[j], y_ls[j], kwargs_lens))
                    diff_mu_rel.append(abs((mu_ler[i] - mu_j) / mu_j))

        assert n_total > 0, "no valid sources sampled from the caustic"
        assert n_match / n_total >= 0.9, (
            f"image-count agreement only {n_match}/{n_total}"
        )
        assert np.max(diff_pos) < 0.02, (
            f"max image-position discrepancy {np.max(diff_pos):.3e} > 0.02"
        )
        assert np.mean(diff_mu_rel) < 0.05, (
            f"mean magnification discrepancy {100 * np.mean(diff_mu_rel):.2f}% > 5%"
        )

    def test_recover_redundant_and_effective_parameters(self, lens_epl_joint):
        """
        Tests
        -----
        - Shared njit image dict: GW / distance fields added as in the rates pipeline.
        - After removing ``EXPECTED_REDUNDANT_KEYS``, ``recover_redundant_parameters``
          restores them (mass sources and ``n_images`` cross-checked).
        - ``produce_effective_params`` adds ``EXPECTED_EFFECTIVE_KEYS`` with shape
          ``(N_SAMPLES, n_max_images)``; first-image identities when ``time_delays==0``
          and minimum image type.
        """
        lens = lens_epl_joint
        lp = lens.epl_shear_sl_parameters_rvs(size=N_SAMPLES)
        base = lens.image_properties_epl_shear_njit(dict(lp))

        m1_src = 30.0
        m2_src = 25.0
        base["mass_1"] = m1_src * (1.0 + base["zs"])
        base["mass_2"] = m2_src * (1.0 + base["zs"])
        base["luminosity_distance"] = lens.luminosity_distance(base["zs"])
        base["geocent_time"] = np.full(N_SAMPLES, 1238166018.0)
        base["ra"] = np.full(N_SAMPLES, 0.5)
        base["dec"] = np.full(N_SAMPLES, 0.1)
        base["phase"] = np.full(N_SAMPLES, 1.0)

        # Strip redundant keys, then recover (same list as LeR when strips are applied).
        for_recovery = copy.deepcopy(base)
        for k in EXPECTED_REDUNDANT_KEYS:
            for_recovery.pop(k, None)
        recovered = lens.recover_redundant_parameters(for_recovery)
        self._assert_param_dict_valid(
            recovered,
            EXPECTED_REDUNDANT_KEYS,
            size=N_SAMPLES,
        )
        np.testing.assert_allclose(recovered["mass_1_source"], m1_src)
        np.testing.assert_allclose(recovered["mass_2_source"], m2_src)
        n_finite = np.sum(~np.isnan(recovered["x0_image_positions"]), axis=1)
        np.testing.assert_array_equal(recovered["n_images"], n_finite.astype(int))

        # Effective parameters mutate the dict; use a copy of the fully prepared base.
        eff_input = copy.deepcopy(base)
        out = lens.produce_effective_params(eff_input)
        nmax = lens.n_max_images
        self._assert_param_dict_valid(
            out,
            EXPECTED_EFFECTIVE_KEYS,
            size=(N_SAMPLES, nmax),
            nan_to_num=True,
        )
        for i in range(N_SAMPLES):
            if out["time_delays"][i, 0] == 0.0 and out["image_type"][i, 0] == 1.0:
                assert out["effective_geocent_time"][i, 0] == base["geocent_time"][i], (
                    f"row {i}: effective_geocent_time mismatch: "
                    f"{out['effective_geocent_time'][i, 0]} != {base['geocent_time'][i]}"
                )
                assert out["effective_phase"][i, 0] == base["phase"][i], (
                    f"row {i}: effective_phase mismatch: "
                    f"{out['effective_phase'][i, 0]} != {base['phase'][i]}"
                )

    @pytest.mark.slow
    def test_njit_vs_lenstronomy_speed(self, lens_epl_joint, interpolator_dir):
        """
        Tests
        -----
        - Coarse wall-time comparison on the same lens-parameter batch:
          ``image_properties_epl_shear_njit`` should be faster than
          ``image_properties_epl_shear_lenstronomy``.

        Notes
        -----
        Not a strict benchmark; multiprocessing and tqdm overhead can dominate
        for very small batches, so we use a moderate batch and a small speedup
        threshold to avoid noise-induced flakes.
        """
        lens = lens_epl_joint
        size = 200
        repeats = 3
        min_speedup = 1.5

        ip_ls = ImageProperties(
            npool=DEFAULT_LENS_CONFIG["npool"],
            image_properties_function="image_properties_epl_shear_lenstronomy",
            multiprocessing_verbose=False,
        )

        # Warm up njit on a small batch to amortize compilation.
        lp_warm = lens.epl_shear_sl_parameters_rvs(size=4)
        _ = lens.image_properties_epl_shear_njit(dict(lp_warm))

        lp = lens.epl_shear_sl_parameters_rvs(size=size)
        t_njit = median_call_time(
            lambda: lens.image_properties_epl_shear_njit(dict(lp)),
            repeats=repeats,
        )
        t_ls = median_call_time(
            lambda: ip_ls.image_properties_epl_shear_lenstronomy(dict(lp)),
            repeats=repeats,
        )

        assert t_njit > 0.0 and t_ls > 0.0, f"invalid timings: t_njit={t_njit}, t_ls={t_ls}"
        speedup = t_ls / t_njit
        assert np.isfinite(speedup), \
            f"invalid speedup: t_njit={t_njit:.4f}s, t_ls={t_ls:.4f}s"
        assert speedup >= min_speedup, (
            f"Expected njit to be faster than lenstronomy, got speedup={speedup:.2f}x "
            f"(njit={t_njit:.4f}s, lenstronomy={t_ls:.4f}s)."
        )
