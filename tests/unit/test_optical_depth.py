"""
Unit tests for OpticalDepth.

This module verifies strong-lensing optical depth and lens-parameter sampling
for the supported lens models in ``ler.lens_galaxy_population.OpticalDepth``.

Test Coverage:
--------------
- Class initialization for supported lens models (EPL+shear / SIE / SIS)
- Optical depth evaluation:
  - ``optical_depth(zs)`` returns finite non-negative values
  - optical depth is close (order-of-magnitude) to the SIS analytic baseline
- Lens-redshift distribution for strongly-lensed events:
  - ``lens_redshift_sl.pdf`` and ``lens_redshift_sl.rvs`` are callable and finite
- Lens-parameter prior objects:
  - for each lens-parameter sampler, ``rvs`` and ``pdf`` are callable and finite
  - for samplers used in numerical integrals, a ``function`` method is also validated
- Cross-section evaluation:
  - cross-section values are finite and non-negative
  - Optional ``slow`` test (EPL+shear): njit vs lenstronomy cross-section accuracy and timing.
"""

import numpy as np
import time
import pytest
from tests_utils import CommonTestUtils
from ler.lens_galaxy_population import OpticalDepth


# ---------------------------------------------------------------------------
# Test configuration
# ---------------------------------------------------------------------------

# Small sample size for fast tests
N_SAMPLES = 10

# Lens models covered by the basic OpticalDepth sanity test.
LENS_TYPES = ["epl_shear_galaxy", "sie_galaxy", "sis_galaxy"]

DEFAULT_CONFIG = dict(
    npool=6,
    z_min=0.0,
    z_max=10.0,
    create_new_interpolator=False,
    verbose=False,
)

# ---------------------------------------------------------------------------
# Small helpers to keep tests readable
# ---------------------------------------------------------------------------

def _make_config(interpolator_dir, **overrides):
    # Start from defaults and allow per-test overrides.
    config = DEFAULT_CONFIG.copy()
    config["directory"] = interpolator_dir
    config.update(overrides)
    return config


@pytest.fixture(scope="module", params=LENS_TYPES)
def optical_depth_obj(interpolator_dir, request):
    """
    Module-scoped fixture: create one OpticalDepth instance per lens type.

    This keeps unit tests fast by reusing already-initialized cached
    interpolators within this module.
    """
    lens_type = request.param
    config = _make_config(interpolator_dir)
    od = OpticalDepth(lens_type=lens_type, **config)
    return lens_type, od


@pytest.fixture(scope="module")
def optical_depth_epl_shear(interpolator_dir):
    """
    Module-scoped fixture: EPL+shear OpticalDepth instance.

    Reused by the lenstronomy comparison test to avoid repeated initialization.
    """
    config = _make_config(interpolator_dir)
    return OpticalDepth(lens_type="epl_shear_galaxy", **config)


class TestOpticalDepth(CommonTestUtils):
    """Tests for OpticalDepth initialization and basic callables."""

    def test_optical_depth_vs_sis_analytic_baseline(self, optical_depth_obj):
        """
        Tests
        -----
        - ``optical_depth(zs)`` returns finite, non-negative values.
        - Optical depth is within a loose order-of-magnitude envelope relative
          to the SIS analytic optical-depth baseline.

        Notes
        -----
        - Very small redshift values are avoided here because spline interpolation
          near boundaries can introduce edge artifacts when the baseline tends to 0.
        """
        lens_type, od = optical_depth_obj

        # Optical depth sanity check: compare with SIS analytic baseline.
        size = 10
        zs = np.linspace(0.1, od.z_max, size)
        tau = od.optical_depth(zs)
        self._assert_array_valid(
            tau, name=f"tau_{lens_type}", size=len(zs), positive=True
        )

        tau_baseline = od.optical_depth_sis_analytic(zs)
        # Loose envelope: allow model differences but reject gross failures.
        assert np.all(tau / tau_baseline < 5.0) and np.all(tau / tau_baseline > 0.2), \
            f"tau/tau_baseline out of [0.2, 5.0]: min={float(np.min(tau/tau_baseline)):.3f}, max={float(np.max(tau/tau_baseline)):.3f}"

    def test_optical_depth_output_value_test(self, optical_depth_obj):
        """
        Tests
        -----
        - Key priors/samplers (pdf/rvs/function) are callable and produce finite
          outputs with expected ranges/sign constraints.
        - Cross-section callables return finite, non-negative values for a set of
          deterministic lens-parameter points.
        """
        lens_type, od = optical_depth_obj

        size = 10
        zs = np.linspace(0.1, od.z_max, size)

        # Strongly-lensed lens-redshift distribution: pdf and rvs are defined.
        # Use a safe deterministic support point: 0 < zl < zs elementwise.
        zl = 0.5 * zs
        self._assert_valid_object(od.lens_redshift_sl, x_array=zl, y_array=zs, name=f"lens_redshift_sl_{lens_type}", size=size, function=False, pdf=True, rvs=True, lo=0.001, hi=od.z_max, positive=True, finite=True)

        # Parameter ranges used for cross-section sanity checks.
        min_sigma = od.lens_functions_params['cross_section_based_sampler']['sigma_min'] 
        max_sigma = od.lens_functions_params['cross_section_based_sampler']['sigma_max']
        sigma = np.linspace(min_sigma, max_sigma, size)

        min_q = od.lens_functions_params['cross_section_based_sampler']['q_min']
        max_q = od.lens_functions_params['cross_section_based_sampler']['q_max']
        q = np.linspace(min_q, max_q, size)

        min_phi = od.lens_functions_params['cross_section_based_sampler']['phi_min']
        max_phi = od.lens_functions_params['cross_section_based_sampler']['phi_max']
        phi = np.linspace(min_phi, max_phi, size)

        min_gamma = od.lens_functions_params['cross_section_based_sampler']['gamma_min']
        max_gamma = od.lens_functions_params['cross_section_based_sampler']['gamma_max']
        gamma = np.linspace(min_gamma, max_gamma, size)

        min_shear = od.lens_functions_params['cross_section_based_sampler']['shear_min']
        max_shear = od.lens_functions_params['cross_section_based_sampler']['shear_max']
        shear = np.linspace(min_shear, max_shear, size)

        # The numerical lens-redshift sampler requires velocity-dispersion bounds.
        assert 'sigma_min' in od.lens_priors_params['velocity_dispersion'], f"sigma_min not found in od.lens_priors_params['velocity_dispersion']; this is required for lens_redshift_strongly_lensed_numerical"
        assert 'sigma_max' in od.lens_priors_params['velocity_dispersion'], f"sigma_max not found in od.lens_priors_params['velocity_dispersion']; this is required for lens_redshift_strongly_lensed_numerical"

        # Velocity dispersion sampler may be conditioned on zl (some models use zl).
        if self._arg_count(od.velocity_dispersion.pdf) == 2:
            self._assert_valid_object(od.velocity_dispersion, x_array=sigma, y_array=zl, name=f"velocity_dispersion_function_{lens_type}", size=size, function=True, pdf=True, rvs=True, lo=min_sigma, hi=max_sigma, positive=True, finite=True)
        elif self._arg_count(od.velocity_dispersion.pdf) == 1:
            self._assert_valid_object(od.velocity_dispersion, x_array=sigma, name=f"velocity_dispersion_function_{lens_type}", size=size, function=True, pdf=True, rvs=True, lo=min_sigma, hi=max_sigma, positive=True, finite=True)
        else:
            raise ValueError(f"velocity_dispersion function has {self._arg_count(od.velocity_dispersion.pdf)} arguments, expected 1 or 2")

        if lens_type in ["epl_shear_galaxy", "sie_galaxy"]:
            # Axis ratio may be conditioned on sigma.
            if self._arg_count(od.axis_ratio.pdf) == 2:
                self._assert_valid_object(od.axis_ratio, x_array=q, y_array=sigma, name=f"axis_ratio_function_{lens_type}", size=size, function=False, pdf=True, rvs=True, lo=min_q, hi=max_q, positive=True, finite=True)
            elif self._arg_count(od.axis_ratio.pdf) == 1:
                self._assert_valid_object(od.axis_ratio, x_array=q, name=f"axis_ratio_function_{lens_type}", size=size, function=False, pdf=True, rvs=True, lo=min_q, hi=max_q, positive=True, finite=True)
            else:
                raise ValueError(f"axis_ratio function has {self._arg_count(od.axis_ratio.pdf)} arguments, expected 1 or 2")

            # Orientation is unconditioned and bounded by [phi_min, phi_max].
            self._assert_valid_object(od.axis_rotation_angle, x_array=phi, name=f"axis_rotation_angle_function_{lens_type}", size=size, function=False, pdf=True, rvs=True, lo=min_phi, hi=max_phi, positive=True, finite=True)

            if lens_type == "epl_shear_galaxy":
                # EPL slope is positive over the configured support.
                self._assert_valid_object(od.density_profile_slope, x_array=gamma, name=f"density_profile_slope_function_{lens_type}", size=size, function=False, pdf=True, rvs=True, lo=min_gamma, hi=max_gamma, positive=True, finite=True)

                # Shear components can be negative; only finite/range checks apply.
                self._assert_valid_object(od.external_shear1, x_array=shear, name=f"external_shear1_function_{lens_type}", size=size, function=False, pdf=True, rvs=True, lo=min_shear, hi=max_shear, positive=False, finite=True)

                self._assert_valid_object(od.external_shear2, x_array=shear, name=f"external_shear2_function_{lens_type}", size=size, function=False, pdf=True, rvs=True, lo=min_shear, hi=max_shear, positive=False, finite=True)


        # Cross-section is a physical area: must be non-negative.
        if self._arg_count(od.cross_section) == 3:
            cs = od.cross_section(zs, zl, sigma)
            self._assert_array_valid(cs, name=f"cross_section_function_{lens_type}", size=size, positive=True, finite=True, lo=0.0, hi=None)
        elif self._arg_count(od.cross_section) == 4:
            cs = od.cross_section(zs, zl, sigma, q)
            self._assert_array_valid(cs, name=f"cross_section_function_{lens_type}", size=size, positive=True, finite=True, lo=0.0, hi=None)
        elif self._arg_count(od.cross_section) == 8:
            cs = od.cross_section(zs, zl, sigma, q, phi, gamma, shear, shear)
            self._assert_array_valid(cs, name=f"cross_section_function_{lens_type}", size=size, positive=True, finite=True, lo=0.0, hi=None)
        else:
            raise ValueError(f"cross_section function has {self._arg_count(od.cross_section)} arguments, expected 3, 4, or 8")


        # Optional: if an EPL+shear cross-section interpolator is enabled, check it.
        if od.lens_priors_params["lens_redshift_sl"] == "lens_redshift_strongly_lensed_numerical":
            if od.lens_priors_params['lens_redshift_sl']['lens_redshift_strongly_lensed_numerical']['cross_section_epl_shear_interpolation']:
                cs_interpolated = od.cross_section_epl_shear_interpolation(
                    zs=zs, zl=zl, sigma=sigma, q=q, phi=phi, gamma=gamma, gamma1=shear, gamma2=shear)
                self._assert_array_valid(cs_interpolated, name=f"cross_section_interpolated_function_{lens_type}", size=size, positive=True, finite=True, lo=0.0, hi=None)

                # Interpolated cross-section should be broadly consistent with direct njit call.
                percent_err = np.abs(cs_interpolated - cs) / np.abs(cs) * 100
                assert np.mean(percent_err) < 5.0, \
                    f"interpolated cross-section mean percent error {np.mean(percent_err):.2f}% >= 5%"
                assert np.max(percent_err) < 500.0, \
                    f"interpolated cross-section max percent error {np.max(percent_err):.2f}% >= 500%"


    @pytest.mark.slow
    @pytest.mark.parametrize("lens_type", ["epl_shear_galaxy"])
    def test_cross_section_comparison_with_lenstronomy(self, optical_depth_epl_shear, lens_type):
        """
        Tests
        -----
        - The njit EPL+shear cross-section implementation is consistent with a
          direct lenstronomy numerical evaluation for the same lens parameters.
        - Basic performance sanity: njit implementation should be faster than
          the pure-python/lenstronomy numerical loop.

        Notes
        -----
        - This is not a strict benchmark. We only assert large, robust timing
          inequalities (njit faster than lenstronomy), and avoid assuming
          multiprocessing is always faster than serial due to process overhead.
        """

        od = optical_depth_epl_shear
        rng = np.random.default_rng(1)

        def generate_parameters(size):
            zs = np.linspace(0.01, od.z_max, size)
            zl = rng.uniform(0.001, zs)

            # parameters for testing
            if self._arg_count(od.velocity_dispersion.rvs) == 2:
                sigma = od.velocity_dispersion.rvs(size, zl)
            elif self._arg_count(od.velocity_dispersion.rvs) == 1:
                sigma = od.velocity_dispersion.rvs(size)
            else:
                raise ValueError(f"velocity_dispersion function has {self._arg_count(od.velocity_dispersion.rvs)} arguments, expected 1 or 2")

            if self._arg_count(od.axis_ratio.rvs) == 2:
                q = od.axis_ratio.rvs(size, sigma)
            elif self._arg_count(od.axis_ratio.rvs) == 1:
                q = od.axis_ratio.rvs(size)
            else:
                raise ValueError(f"axis_ratio function has {self._arg_count(od.axis_ratio.rvs)} arguments, expected 1 or 2")

            phi = od.axis_rotation_angle.rvs(size)
            gamma = od.density_profile_slope.rvs(size)
            gamma1 = od.external_shear1.rvs(size)
            gamma2 = od.external_shear2.rvs(size)

            return zs, zl, sigma, q, phi, gamma, gamma1, gamma2

        # --------------------------
        # ler's njit cross_section
        # --------------------------

        # warm up the njit cross_section
        time_start = time.time()
        zs, zl, sigma, q, phi, gamma, gamma1, gamma2 = generate_parameters(2)
        time_end = time.time()
        njit_compile_time_lens_parameters = time_end - time_start
        assert njit_compile_time_lens_parameters < 30.0, f"njit compile time for lens parameters is too long: {njit_compile_time_lens_parameters} seconds"

        time_start = time.time()
        cs_njit = od.cross_section(zs, zl, sigma, q, phi, gamma, gamma1, gamma2)
        time_end = time.time()
        njit_compile_time_cross_section = time_end - time_start
        assert njit_compile_time_cross_section < 30.0, f"njit compile time for cross_section is too long: {njit_compile_time_cross_section} seconds"

        # Keep this moderate: unit tests should be quick but still reduce noise.
        size = 2000
        zs, zl, sigma, q, phi, gamma, gamma1, gamma2 = generate_parameters(size)
        time_start = time.time()
        cs_njit = od.cross_section(zs, zl, sigma, q, phi, gamma, gamma1, gamma2)
        time_end = time.time()
        time_cs_njit = time_end - time_start

        # --------------------------
        # vs lenstronomy cross_section with multiprocessing
        # --------------------------
        theta_E = od.compute_einstein_radii(sigma, zl, zs)
        time_start = time.time()
        cs_lenstronomy = od.cross_section_epl_shear_numerical_mp(
            theta_E,
            gamma,
            gamma1,
            gamma2,
            q,
            phi,
        )
        time_end = time.time()
        time_cs_lenstronomy_mp = time_end - time_start

        self._assert_array_valid(cs_lenstronomy, name=f"cross_section_lenstronomy_function_{lens_type}", size=size, positive=True, finite=True, lo=0.0, hi=None)

        # assert lenstronomy cross_section values are close to od.cross_section . mean percent error < 5%
        percent_err = np.abs(cs_lenstronomy - cs_njit) / np.abs(cs_njit) * 100
        assert np.mean(percent_err) < 5.0, \
            f"lenstronomy vs njit cross-section mean percent error {np.mean(percent_err):.2f}% >= 5%"
        assert np.max(percent_err) < 500.0, \
            f"lenstronomy vs njit cross-section max percent error {np.max(percent_err):.2f}% >= 500%"

        # --------------------------
        # vs lenstronomy cross_section without multiprocessing
        # --------------------------
        time_start = time.time()
        cs_lenstronomy_serial = od.cross_section_epl_shear_numerical(
            zs, zl, sigma, q, phi, gamma, gamma1, gamma2
        )
        time_end = time.time()
        time_cs_lenstronomy = time_end - time_start

        # Shapes and basic positivity (physical area).
        self._assert_array_valid(
            cs_lenstronomy_serial,
            name=f"cross_section_lenstronomy_serial_function_{lens_type}",
            size=size,
            positive=True,
            finite=True,
            lo=0.0,
            hi=None,
        )

        # compare the runtime
        assert time_cs_lenstronomy_mp > time_cs_njit, \
            f"expected lenstronomy (mp) to be slower than njit, but got lenstronomy_mp={time_cs_lenstronomy_mp:.4f}s vs njit={time_cs_njit:.4f}s"
        assert time_cs_lenstronomy > time_cs_njit, \
            f"expected lenstronomy (serial) to be slower than njit, but got lenstronomy={time_cs_lenstronomy:.4f}s vs njit={time_cs_njit:.4f}s"
        # Do not assert ordering between mp and serial: depends on overhead / cores.
