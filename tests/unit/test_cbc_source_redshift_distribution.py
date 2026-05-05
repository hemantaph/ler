"""
Unit tests for CBCSourceRedshiftDistribution class.

This module validates the redshift distribution sampling functionality
of the CBCSourceRedshiftDistribution class, including class initialization,
merger rate density evaluation, and source redshift sampling.

Test Coverage:
--------------
- Class initialization for supported event types (BBH/BNS/NSBH)
- Merger rate density evaluation (array input) and detector-frame consistency:
  R_det(z) = R_src(z)/(1+z) * dVc/dz
- Source redshift sampling: output size, finiteness, and range within [z_min, z_max]
- Normalization constant `normalization_pdf_z` is finite and positive
- Cosmological interpolators (luminosity distance and dVc/dz) agree with `astropy`
  for a custom cosmology
- Custom merger-rate-density parameters for:
  - ler-available merger rate density models (string-based)
  - user-provided callables (with identifier metadata)
- Sampler sanity check: consecutive draws differ (non-deterministic sampler)
- Input validation: invalid `event_type` and invalid merger-rate-density model string raise `ValueError`
- Optional ``slow`` test: no-JIT subprocess vs in-process njit redshift sampling.
"""

import os
import subprocess
import sys
import numpy as np
import pytest
from astropy.cosmology import LambdaCDM
from tests_utils import CommonTestUtils, clamp_npool_for_numba, median_call_time
from ler.gw_source_population import CBCSourceRedshiftDistribution

# ---------------------------------------------------------------------------
# Test configuration
# ---------------------------------------------------------------------------

# Small sample size to keep tests fast
N_SAMPLES = 10

# Default configuration – uses existing cached interpolators for speed
DEFAULT_CONFIG = dict(
    npool=clamp_npool_for_numba(6),
    z_min=0.0,
    z_max=10.0,
    create_new_interpolator=False,
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


def _make_cbc(interpolator_dir, event_type="BBH", **kwargs):
    # Centralize object creation so tests stay concise.
    return CBCSourceRedshiftDistribution(
        event_type=event_type,
        **_make_config(interpolator_dir, **kwargs),
    )


class TestCBCSourceRedshiftDistribution(CommonTestUtils):
    """Tests for CBCSourceRedshiftDistribution initialization and sampling."""

    @pytest.mark.parametrize("event_type", ["BBH", "BNS", "NSBH"])
    def test_init(self, interpolator_dir, event_type):
        """
        Tests
        -----
        - CBCSourceRedshiftDistribution initializes without error for supported event types.
        - Key attributes (z_min, z_max, event_type, cosmo) are set correctly.
        - normalization_pdf_z is a finite positive float.
        """
        # Create a sampler for the requested event type.
        config = _make_config(interpolator_dir)
        cbc = CBCSourceRedshiftDistribution(event_type=event_type, **config)

        # Validate basic attributes and normalization.
        assert cbc.z_min == config["z_min"], f"z_min: expected {config['z_min']}, got {cbc.z_min}"
        assert cbc.z_max == config["z_max"], f"z_max: expected {config['z_max']}, got {cbc.z_max}"
        assert cbc.event_type == event_type, f"event_type: expected {event_type}, got {cbc.event_type}"
        assert cbc.cosmo is not None, "cosmo must not be None after initialization"
        assert np.isfinite(cbc.normalization_pdf_z), f"normalization_pdf_z is not finite: {cbc.normalization_pdf_z}"
        assert cbc.normalization_pdf_z > 0, f"normalization_pdf_z must be positive, got {cbc.normalization_pdf_z}"

    def test_merger_rate_density_and_source_sampling(self, interpolator_dir):
        """
        Tests
        -----
        - merger_rate_density returns finite positive values for array input.
        - merger_rate_density_detector_frame returns finite positive values.
        - Detector-frame relation holds: R_det(z) = R_src(z)/(1+z) * dVc/dz.
        - Sampled redshifts have correct size (N_SAMPLES).
        - All samples are within [z_min, z_max].
        - All values are finite.
        """
        # Create a BBH sampler for deterministic evaluation points.
        config = _make_config(interpolator_dir)
        cbc = CBCSourceRedshiftDistribution(event_type="BBH", **config)

        # Check source-frame and detector-frame merger-rate definitions.
        z_arr = np.array([0.1, 0.5, 1.0, 2.0])
        rate_src = cbc.merger_rate_density(z_arr)
        self._assert_array_valid(
            rate_src,
            name="merger_rate_density",
            size=len(z_arr),
            positive=True,
        )

        rate_det = cbc.merger_rate_density_detector_frame(z_arr)
        self._assert_array_valid(
            rate_det,
            name="merger_rate_density_detector_frame",
            size=len(z_arr),
            positive=True,
        )

        dvc = cbc.differential_comoving_volume(z_arr)
        expected_det = rate_src / (1.0 + z_arr) * dvc
        np.testing.assert_allclose(rate_det, expected_det, rtol=1e-3)

        # Draw redshift samples and validate range + finiteness.
        zs = cbc.merger_rate_density_based_source_redshift(N_SAMPLES)
        self._assert_array_valid(
            zs,
            name="zs",
            size=N_SAMPLES,
            lo=config["z_min"],
            hi=config["z_max"],
        )

    def test_cosmological_functions_interpolator(self, interpolator_dir):
        """
        Tests
        -----
        - luminosity_distance interpolator matches astropy for a custom cosmology.
        - differential_comoving_volume interpolator matches astropy for a custom cosmology.
        - Both return finite positive values.
        """
        # Custom cosmology should propagate into interpolators.
        config = _make_config(interpolator_dir)
        cosmo = LambdaCDM(H0=67.4, Om0=0.315, Ode0=0.685, Tcmb0=2.725)
        cbc = CBCSourceRedshiftDistribution(
            event_type="BBH",
            cosmology=cosmo,
            **config,
        )

        z_arr = np.array([0.1, 0.5, 1.0])
        dl = cbc.luminosity_distance(z_arr)
        self._assert_array_valid(dl, name="luminosity_distance", size=3, positive=True)
        dl_astropy = cosmo.luminosity_distance(z_arr).value
        np.testing.assert_allclose(dl, dl_astropy, rtol=2e-3)

        dvc = cbc.differential_comoving_volume(z_arr)
        self._assert_array_valid(dvc, name="dVc_dz", size=3, positive=True)
        dvc_astropy = cosmo.differential_comoving_volume(z_arr).value * 4 * np.pi
        np.testing.assert_allclose(dvc, dvc_astropy, rtol=2e-3)

    def test_custom_merger_rate_density_parameters(self, interpolator_dir):
        """
        Tests
        -----
        - Custom parameters work for (i) ler-available merger rate density functions.
        - Custom parameters work for (ii) user-provided merger rate density callables.
        """
        # This test is intentionally split into two parts:
        # (i) builtin ler function with user-provided parameters,
        # (ii) user-provided callable with its identifier metadata.
        config = _make_config(interpolator_dir)

        # ------------------------------------------------------------------
        # Part (i): ler-available function + user-provided parameter override
        # ------------------------------------------------------------------
        cbc_builtin_default = CBCSourceRedshiftDistribution(
            event_type="BBH",
            merger_rate_density="merger_rate_density_bbh_oguri2018",
            **config,
        )
        default_R0 = float(cbc_builtin_default.merger_rate_density_param["R0"])

        cbc_builtin_custom = CBCSourceRedshiftDistribution(
            event_type="BBH",
            merger_rate_density="merger_rate_density_bbh_oguri2018",
            merger_rate_density_param=dict(R0=2.0 * default_R0),
            **config,
        )

        z0 = 0.5
        r_default = np.asarray(cbc_builtin_default.merger_rate_density(z0)).item()
        r_custom = np.asarray(cbc_builtin_custom.merger_rate_density(z0)).item()
        assert np.isfinite(r_default) and r_default > 0.0, f"default R0 rate not finite/positive: {r_default}"
        assert np.isfinite(r_custom) and r_custom > 0.0, f"custom R0 rate not finite/positive: {r_custom}"
        # For this model, R0 is a multiplicative normalization.
        np.testing.assert_allclose(r_custom / r_default, 2.0, rtol=1e-2)

        # ------------------------------------------------------------------
        # Part (ii): user-provided callable + identifier metadata
        # ------------------------------------------------------------------
        def custom_rate(z):
            z = np.asarray(z)
            return 1.0e-9 * (1.0 + 0.4 * z + 0.1 * z**2)

        cbc = CBCSourceRedshiftDistribution(
            event_type="BBH",
            merger_rate_density=custom_rate,
            merger_rate_density_param=dict(
                param_name="merger_rate_density",
                function_type="custom_rate_poly",
            ),
            **config,
        )

        z_arr = np.array([0.2, 0.7, 1.4, 2.1])
        rate_src = cbc.merger_rate_density(z_arr)
        expected = custom_rate(z_arr)
        self._assert_array_valid(rate_src, name="custom_merger_rate", size=len(z_arr), positive=True)
        np.testing.assert_allclose(rate_src, expected, rtol=1e-2)

    def test_sampling_reproducibility(self, interpolator_dir):
        """
        Tests
        -----
        - Consecutive sampler calls produce different draws (sanity check).

        Notes
        -----
        The sampler may use internal RNG state that is not fully controlled by
        `np.random.seed()`, so we avoid asserting strict equality under reseeding.
        """
        # Sanity check: consecutive draws should differ.
        cbc = _make_cbc(interpolator_dir, event_type="BBH")

        np.random.seed(42)
        zs1 = cbc.merger_rate_density_based_source_redshift(20)
        zs2 = cbc.merger_rate_density_based_source_redshift(20)
        assert np.any(zs1 != zs2), "consecutive samples should differ"

    def test_invalid_event_type_raises(self, interpolator_dir):
        config = _make_config(interpolator_dir)
        with pytest.raises(ValueError):
            _ = CBCSourceRedshiftDistribution(event_type="NOT_A_REAL_TYPE", **config)

    def test_invalid_merger_rate_density_string_raises(self, interpolator_dir):
        config = _make_config(interpolator_dir)
        with pytest.raises(ValueError):
            _ = CBCSourceRedshiftDistribution(
                event_type="BBH",
                merger_rate_density="definitely_not_a_model_name",
                **config,
            )

    @pytest.mark.slow
    def test_njit_speed(self, interpolator_dir):
        """
        Tests
        -----
        - Compare speed of the sampler with numba JIT disabled (no-JIT baseline)
          against njit-backed source-redshift sampling (post-compilation).
        - The no-JIT baseline is executed in a subprocess so that
          NUMBA_DISABLE_JIT=1 is in effect before any numba import.
        - The njit path is warmed up before timing, and the timed section is
          repeated to reduce noise.
        """
        sample_size = 100000
        repeats = 3
        # Require a small speed-up to account for timing noise.
        # (Perf tests are inherently noisy across machines/loads.)
        min_speedup = 1.01

        # ------------------------------------------------------------------
        # Test 1: no-JIT baseline — run in a subprocess with NUMBA_DISABLE_JIT=1
        # so the env var is set before numba is first imported.
        # ------------------------------------------------------------------
        script_no_jit = "\n".join(
            [
                "import time",
                "import numpy as np",
                "from ler.gw_source_population import CBCSourceRedshiftDistribution",
                "cbc = CBCSourceRedshiftDistribution(",
                "    npool=1, event_type='BBH',",
                "    z_min=0.0, z_max=10.0,",
                "    create_new_interpolator=False,",
                f"    directory={repr(interpolator_dir)},",
                ")",
                # one warm-up call (still no-JIT) to reduce one-off overheads
                "_ = cbc.merger_rate_density_based_source_redshift(2)",
                "times = []",
                f"sample_size = {sample_size}",
                f"repeats = {repeats}",
                "for _ in range(repeats):",
                "    t0 = time.perf_counter()",
                "    _ = cbc.merger_rate_density_based_source_redshift(sample_size)",
                "    times.append(time.perf_counter() - t0)",
                "import sys; print(float(np.median(np.asarray(times))), file=sys.stderr)",
            ]
        )
        env_no_jit = os.environ.copy()
        env_no_jit.update(dict(NUMBA_DISABLE_JIT="1"))
        proc = subprocess.run(
            [sys.executable, "-c", script_no_jit],
            capture_output=True,
            text=True,
            env=env_no_jit,
        )
        assert proc.returncode == 0, (
            f"No-JIT subprocess failed:\nSTDOUT: {proc.stdout}\nSTDERR: {proc.stderr}"
        )
        # timing is written to stderr to avoid mixing with CBCSourceRedshiftDistribution banner output
        t_no_jit = float(proc.stderr.strip().splitlines()[-1])

        # ------------------------------------------------------------------
        # Test 2: njit-backed sampler (in-process, JIT enabled).
        # ------------------------------------------------------------------
        config = DEFAULT_CONFIG.copy()
        config["directory"] = interpolator_dir
        config.pop("npool", None)

        # Use npool=1 to avoid multiprocessing overhead dominating the comparison.
        cbc_njit = CBCSourceRedshiftDistribution(npool=1, event_type="BBH", **config)
        # Warm up (compile) njit path before timing.
        _ = cbc_njit.merger_rate_density_based_source_redshift(2)

        t_njit = median_call_time(
            lambda: cbc_njit.merger_rate_density_based_source_redshift(sample_size),
            repeats=repeats,
        )

        # Both paths must produce valid positive timings.
        assert t_no_jit > 0.0 and t_njit > 0.0, f"invalid timings: t_no_jit={t_no_jit}, t_njit={t_njit}"
        speed_ratio = t_no_jit / t_njit

        assert np.isfinite(speed_ratio), (
            f"invalid speed ratio: t_no_jit={t_no_jit:.6f}s, t_njit={t_njit:.6f}s"
        )

        assert speed_ratio >= min_speedup, (
            f"Expected njit to be faster after warm-up, but got speedup={speed_ratio:.2f}x "
            f"(no-JIT={t_no_jit:.6f}s, njit={t_njit:.6f}s)."
        )
