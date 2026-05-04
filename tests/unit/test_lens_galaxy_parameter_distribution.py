"""
Unit tests for LensGalaxyParameterDistribution.

This module verifies lens-galaxy parameter sampling and lensed-source redshift
sampling for ``ler.lens_galaxy_population.LensGalaxyParameterDistribution``.

Test Coverage:
--------------
- Class initialization sanity checks (attributes, cached interpolators)
- Strongly-lensed source redshift sampler:
  - ``strongly_lensed_source_redshift`` returns finite samples within [z_min, z_max]
  - ``get_attribute=True`` returns a ``FunctionConditioning`` object
- Cross-section-based EPL+shear lens-parameter sampling (``epl_shear_sl_parameters_rvs``):
  - ``test_cross_section_based_sampler_variants``: ``importance_sampler_partial``
    on default runs (slow tests omitted); uses ``n_prop`` from
    ``N_PROP_CROSS_SECTION`` for a fast CI-friendly check.
  - ``test_cross_section_based_sampler_variants_heavy`` *(slow)*: the three
    slower samplers — ``rejection_sampler_full``, ``importance_sampler_full``,
    ``rejection_sampler_partial`` — same output assertions as the fast test.
  - ``test_njit_speed`` *(slow)*: no-JIT subprocess vs in-process njit for
    ``importance_sampler_partial`` with ``npool=1`` and ``npool=6``.
- Intrinsic lens-parameter sampling (no strong-lensing weighting):
  - ``sample_all_routine_epl_shear_intrinsic`` returns finite samples with expected keys.
"""

import numpy as np
import pytest
from tests_utils import CommonTestUtils, median_call_time
from ler.lens_galaxy_population import LensGalaxyParameterDistribution
from ler.utils import FunctionConditioning


# ---------------------------------------------------------------------------
# Test configuration
# ---------------------------------------------------------------------------

N_SAMPLES = 10

DEFAULT_CONFIG = dict(
    npool=6,
    z_min=0.0,
    z_max=10.0,
    create_new_interpolator=False,
    verbose=False,
    buffer_size=200,
)

# Fast sampler: exercised on every default run (slow-marked tests deselected).
CROSS_SECTION_SAMPLER_FAST = "importance_sampler_partial"

# Heavier cross-section samplers: only under ``@pytest.mark.slow`` (slow rejection /
# importance paths with larger proposal counts).
CROSS_SECTION_SAMPLERS_SLOW_ONLY = [
    "rejection_sampler_full",
    "importance_sampler_full",
    "rejection_sampler_partial",
]

# Proposal counts per sampler (interface/sanity coverage, not production efficiency).
N_PROP_CROSS_SECTION = dict(
    rejection_sampler_full=500,
    importance_sampler_full=100,
    rejection_sampler_partial=100,
    importance_sampler_partial=50,
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


def _run_cross_section_sampler_variant(lens, sampler_type):
    """
    Reconfigure ``lens`` for the given cross-section sampler, sample EPL+shear
    lens parameters, and return the parameter dict.

    Mutates ``lens.lens_functions``, ``lens.lens_functions_params``, and
    ``lens.cross_section_based_sampler`` in place (same pattern as the tests).
    """
    lens.lens_functions["cross_section_based_sampler"] = sampler_type
    lens.lens_functions_params["cross_section_based_sampler"]["n_prop"] = (
        N_PROP_CROSS_SECTION[sampler_type]
    )
    lens.cross_section_based_sampler = lens._initialization_cross_section_sampler()
    return lens.epl_shear_sl_parameters_rvs(size=N_SAMPLES)


@pytest.fixture(scope="module")
def lens_epl_shear(interpolator_dir):
    """
    Module-scoped LensGalaxyParameterDistribution instance for EPL+shear tests.

    Initialization can be expensive due to cached interpolator setup; this fixture
    reuses the same instance for multiple tests in this module.
    """
    cfg = _make_config(interpolator_dir, npool=6)
    return LensGalaxyParameterDistribution(lens_type="epl_shear_galaxy", **cfg)


@pytest.fixture(scope="module")
def lens_epl_shear_npool1(interpolator_dir):
    """Same as `lens_epl_shear` but with `npool=1` (speed/warmup comparisons)."""
    cfg = _make_config(interpolator_dir, npool=1)
    return LensGalaxyParameterDistribution(lens_type="epl_shear_galaxy", **cfg)


class TestLensGalaxyParameterDistribution(CommonTestUtils):
    """Tests for LensGalaxyParameterDistribution initialization and sampling."""

    def _assert_cross_section_sampler_output(self, lens, params):
        """
        Assert dict shape, keys, finiteness, and basic supports for outputs of
        ``epl_shear_sl_parameters_rvs`` after a cross-section-based sampler run.

        Shared by ``test_cross_section_based_sampler_variants`` and
        ``test_cross_section_based_sampler_variants_heavy``.
        """
        assert isinstance(params, dict), f"epl_shear_sl_parameters_rvs: expected dict, got {type(params)}"

        expected_keys = ["zl", "zs", "sigma", "theta_E", "q", "phi", "gamma", "gamma1", "gamma2"]
        for k in expected_keys:
            assert k in params, f"missing key: {k}"
            self._assert_array_valid(params[k], name=k, size=N_SAMPLES, finite=True)

        assert np.all(params["zs"] >= lens.z_min) and np.all(params["zs"] <= lens.z_max), \
            f"zs out of [z_min={lens.z_min}, z_max={lens.z_max}]: min={params['zs'].min():.3f}, max={params['zs'].max():.3f}"
        assert np.all(params["zl"] >= lens.z_min) and np.all(params["zl"] <= params["zs"]), \
            f"zl out of [z_min={lens.z_min}, zs]: min={params['zl'].min():.3f}"
        assert np.all(params["sigma"] > 0.0), f"sigma must be positive, min={params['sigma'].min():.3f}"
        assert np.all(params["q"] > 0.0) and np.all(params["q"] <= 1.0), \
            f"q must be in (0, 1]: min={params['q'].min():.3f}, max={params['q'].max():.3f}"
        assert np.all(params["theta_E"] >= 0.0), f"theta_E must be >= 0, min={params['theta_E'].min():.3f}"

    def test_init(self, lens_epl_shear):
        """
        Tests
        -----
        - LensGalaxyParameterDistribution initializes without error.
        - Key attributes are set and finite:
          - z_min, z_max, lens_type, npool
          - normalization_pdf_z_lensed is finite and positive
        - Key samplers are present as FunctionConditioning objects:
          - zs_sl (strongly-lensed source redshift)
          - lens_redshift_sl (strongly-lensed lens redshift, conditioned on zs)
        """
        lens = lens_epl_shear
        assert lens.lens_type == "epl_shear_galaxy", f"lens_type: expected 'epl_shear_galaxy', got {lens.lens_type}"
        assert lens.z_min == DEFAULT_CONFIG["z_min"], f"z_min: expected {DEFAULT_CONFIG['z_min']}, got {lens.z_min}"
        assert lens.z_max == DEFAULT_CONFIG["z_max"], f"z_max: expected {DEFAULT_CONFIG['z_max']}, got {lens.z_max}"
        assert lens.npool == 6, f"npool: expected 6, got {lens.npool}"

        assert np.isfinite(lens.normalization_pdf_z_lensed), f"normalization_pdf_z_lensed is not finite: {lens.normalization_pdf_z_lensed}"
        assert lens.normalization_pdf_z_lensed > 0.0, f"normalization_pdf_z_lensed must be positive, got {lens.normalization_pdf_z_lensed}"

        assert isinstance(lens.zs_sl, FunctionConditioning), \
            f"zs_sl: expected FunctionConditioning, got {type(lens.zs_sl)}"
        assert isinstance(lens.lens_redshift_sl, FunctionConditioning), \
            f"lens_redshift_sl: expected FunctionConditioning, got {type(lens.lens_redshift_sl)}"

    def test_strongly_lensed_source_redshift(self, lens_epl_shear):
        """
        Tests
        -----
        - strongly_lensed_source_redshift(size) returns finite redshifts in [z_min, z_max].
        - strongly_lensed_source_redshift(get_attribute=True) returns FunctionConditioning.
        """
        lens = lens_epl_shear

        zs_obj = lens.strongly_lensed_source_redshift(size=10, get_attribute=True)
        assert isinstance(zs_obj, FunctionConditioning), \
            f"get_attribute=True: expected FunctionConditioning, got {type(zs_obj)}"

        zs = lens.strongly_lensed_source_redshift(size=N_SAMPLES)
        self._assert_array_valid(
            zs,
            name="strongly_lensed_source_redshift",
            size=N_SAMPLES,
            finite=True,
            lo=lens.z_min,
            hi=lens.z_max,
        )

    def test_cross_section_based_sampler_variants(self, lens_epl_shear):
        """
        Tests
        -----
        - ``importance_sampler_partial`` (fast path for CI): after reconfiguring
          via ``_run_cross_section_sampler_variant``, ``epl_shear_sl_parameters_rvs``
          returns arrays with expected keys and size; values are finite and in
          basic physical ranges. Uses ``N_PROP_CROSS_SECTION['importance_sampler_partial']``.

        Notes
        -----
        Heavier samplers are covered in ``test_cross_section_based_sampler_variants_heavy``
        (``@pytest.mark.slow``).
        """
        lens = lens_epl_shear
        params = _run_cross_section_sampler_variant(lens, CROSS_SECTION_SAMPLER_FAST)
        self._assert_cross_section_sampler_output(lens, params)

    @pytest.mark.slow
    @pytest.mark.parametrize("sampler_type", CROSS_SECTION_SAMPLERS_SLOW_ONLY)
    def test_cross_section_based_sampler_variants_heavy(self, lens_epl_shear, sampler_type):
        """
        Tests
        -----
        - For each slow cross-section-based sampler type — parametrized over
          ``rejection_sampler_full``, ``importance_sampler_full``, and
          ``rejection_sampler_partial`` — ``epl_shear_sl_parameters_rvs`` output
          passes the same checks as ``test_cross_section_based_sampler_variants``
          (via ``_assert_cross_section_sampler_output``). Proposal counts use
          ``N_PROP_CROSS_SECTION[sampler_type]``.

        Excluded from default runs by the ``slow`` marker (see
        ``tests/conftest.py``).
        """
        lens = lens_epl_shear
        params = _run_cross_section_sampler_variant(lens, sampler_type)
        self._assert_cross_section_sampler_output(lens, params)

    def test_intrinsic_lens_parameter_sampling(self, lens_epl_shear):
        """
        Tests
        -----
        - sample_all_routine_epl_shear_intrinsic returns expected lens-parameter keys.
        - Outputs are finite and satisfy basic support constraints.
        """
        lens = lens_epl_shear
        params = lens.sample_all_routine_epl_shear_intrinsic(size=N_SAMPLES)
        assert isinstance(params, dict), f"sample_all_routine_epl_shear_intrinsic: expected dict, got {type(params)}"

        expected_keys = ["zl", "zs", "sigma", "q", "phi", "gamma", "gamma1", "gamma2"]
        for k in expected_keys:
            assert k in params, f"missing key: {k}"
            self._assert_array_valid(params[k], name=k, size=N_SAMPLES, finite=True)

        assert np.all(params["zs"] >= lens.z_min) and np.all(params["zs"] <= lens.z_max), \
            f"zs out of [z_min={lens.z_min}, z_max={lens.z_max}]: min={params['zs'].min():.3f}, max={params['zs'].max():.3f}"
        assert np.all(params["zl"] >= lens.z_min) and np.all(params["zl"] <= lens.z_max), \
            f"zl out of [z_min={lens.z_min}, z_max={lens.z_max}]: min={params['zl'].min():.3f}, max={params['zl'].max():.3f}"
        assert np.all(params["sigma"] > 0.0), f"sigma must be positive, min={params['sigma'].min():.3f}"
        assert np.all(params["q"] > 0.0) and np.all(params["q"] <= 1.0), \
            f"q must be in (0, 1]: min={params['q'].min():.3f}, max={params['q'].max():.3f}"

    @pytest.mark.slow
    def test_njit_speed(self, interpolator_dir, lens_epl_shear_npool1, lens_epl_shear):
        """
        Tests
        -----
        - Compare wall time of cross-section-based ``epl_shear_sl_parameters_rvs``
          with ``importance_sampler_partial`` (same sampler as the default
          fast test ``test_cross_section_based_sampler_variants``) in three modes:
          i) njit-backed sampler with npool=6
          ii) njit-backed sampler with npool=1
          iii) no-JIT baseline (NUMBA_DISABLE_JIT=1) with npool=1
        - The no-JIT baseline is executed in a subprocess so that
          NUMBA_DISABLE_JIT=1 is in effect before any numba import.
        - njit paths are warmed up before timing, and the timed section is
          repeated to reduce noise.

        Notes
        -----
        Coarse speed sanity check (not a strict benchmark). Slower cross-section
        modes are correctness-tested in ``test_cross_section_based_sampler_variants_heavy``,
        not timed here.
        """
        sample_size = 16
        repeats = 3
        min_speedup = 1.01

        # ------------------------------------------------------------------
        # Test 1: no-JIT baseline (subprocess)
        # ------------------------------------------------------------------
        import os
        import subprocess
        import sys

        script_no_jit = "\n".join(
            [
                "import time",
                "import numpy as np",
                "from ler.lens_galaxy_population import LensGalaxyParameterDistribution",
                "lens = LensGalaxyParameterDistribution(",
                "    npool=1, lens_type='epl_shear_galaxy',",
                "    z_min=0.0, z_max=10.0,",
                "    create_new_interpolator=False,",
                f"    directory={repr(interpolator_dir)},",
                "    buffer_size=200,",
                ")",
                "lens.lens_functions['cross_section_based_sampler'] = 'importance_sampler_partial'",
                "lens.cross_section_based_sampler = lens._initialization_cross_section_sampler()",
                # one warm-up call to reduce one-off overheads
                "_ = lens.epl_shear_sl_parameters_rvs(10)",
                "times = []",
                f"sample_size = {sample_size}",
                f"repeats = {repeats}",
                "for _ in range(repeats):",
                "    t0 = time.perf_counter()",
                "    _ = lens.epl_shear_sl_parameters_rvs(sample_size)",
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
        t_no_jit = float(proc.stderr.strip().splitlines()[-1])

        # ------------------------------------------------------------------
        # Test 2/3: njit-backed sampler (in-process)
        # ------------------------------------------------------------------
        def _time_sampler(lens):
            lens.lens_functions["cross_section_based_sampler"] = "importance_sampler_partial"
            lens.cross_section_based_sampler = lens._initialization_cross_section_sampler()
            _ = lens.epl_shear_sl_parameters_rvs(10)
            return median_call_time(
                lambda: lens.epl_shear_sl_parameters_rvs(sample_size),
                repeats=repeats,
            )

        t_njit_1 = _time_sampler(lens_epl_shear_npool1)
        t_njit_6 = _time_sampler(lens_epl_shear)

        assert t_no_jit > 0.0 and t_njit_1 > 0.0 and t_njit_6 > 0.0, \
            f"invalid timings: t_no_jit={t_no_jit}, t_njit_1={t_njit_1}, t_njit_6={t_njit_6}"

        speedup_1 = t_no_jit / t_njit_1
        speedup_6 = t_no_jit / t_njit_6

        assert np.isfinite(speedup_1) and np.isfinite(speedup_6), \
            f"invalid speedups: speedup_1={speedup_1}, speedup_6={speedup_6}"
        assert speedup_1 >= min_speedup, (
            f"Expected njit (npool=1) to be faster after warm-up, but got speedup={speedup_1:.2f}x "
            f"(no-JIT={t_no_jit:.6f}s, njit npool=1={t_njit_1:.6f}s)."
        )
        assert speedup_6 >= min_speedup, (
            f"Expected njit (npool=6) to be faster after warm-up, but got speedup={speedup_6:.2f}x "
            f"(no-JIT={t_no_jit:.6f}s, njit npool=6={t_njit_6:.6f}s)."
        )
