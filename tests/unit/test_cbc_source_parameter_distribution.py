"""
Unit tests for CBCSourceParameterDistribution.

This module verifies compact-binary source-parameter sampling across event
types, spin configurations, custom prior overrides, and invalid input
handling for CBCSourceParameterDistribution.

Test Coverage:
--------------
- Output structure and finite sampled values for BBH, BNS, and NSBH
- Absence of spin and precession keys when ``spin_zero=True``
- Presence of precessing-spin keys for BBH when ``spin_precession=True``
- Custom ``gw_priors`` callables for redshift, mass, angular, time, and spin parameters
- Custom ``gw_priors_params`` overrides for built-in prior samplers
- Both supported custom mass pathways: ``mass_ratio`` and ``mass_2_source``
- Sampling reproducibility for custom ``gw_priors`` under explicit reseeding
- Fixed-parameter overrides via ``param`` in ``sample_gw_parameters``
- Invalid ``event_type`` rejection via ``ValueError``
- ``gw_parameters_rvs_njit`` output sanity (structure and valid ranges) for BBH with ``spin_zero=True``
- Optional ``slow`` test: ``gw_parameters_rvs_njit`` no-JIT subprocess vs in-process JIT (``npool=1`` and clamped parallel ``npool``, typically ``6``).
"""

import numpy as np
import pytest
from ler.gw_source_population import CBCSourceParameterDistribution
from tests_utils import (
    CommonTestUtils,
    EXPECTED_SPIN_PRECESSING_KEYS,
    EXPECTED_SPIN_ZERO_FORBIDDEN_BBH_KEYS,
    clamp_npool_for_numba,
    median_call_time,
)

# Matches ``desired=6`` on large hosts; clamps for Numba on low-vCPU CI runners.
NPOOL_PARALLEL = clamp_npool_for_numba(6)

# ---------------------------------------------------------------------------
# Test configuration
# ---------------------------------------------------------------------------

# Small sample size for fast tests
N_SAMPLES = 30

# Shared configuration
DEFAULT_CONFIG = dict(
    npool=NPOOL_PARALLEL,
    z_min=0.0,
    z_max=10.0,
    create_new_interpolator=False,
)

# Minimum expected keys in all GW parameter samples
BASE_KEYS = [
    "zs", "mass_1_source", "mass_2_source",
    "luminosity_distance", "mass_1", "mass_2",
    "geocent_time", "ra", "dec", "phase", "psi", "theta_jn",
]


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
    return CBCSourceParameterDistribution(
        event_type=event_type,
        **_make_config(interpolator_dir, **kwargs),
    )


class TestCBCSourceParameterDistribution(CommonTestUtils):
    """Tests for CBC source parameter sampling and customization."""

    @pytest.mark.parametrize(
        "event_type, spin_zero, forbidden_spin_keys",
        [
            ("BBH", True, EXPECTED_SPIN_ZERO_FORBIDDEN_BBH_KEYS),
            ("BNS", True, EXPECTED_SPIN_PRECESSING_KEYS),
            ("NSBH", True, EXPECTED_SPIN_PRECESSING_KEYS),
        ],
    )
    def test_output_sanity_spin_zero(self, interpolator_dir, event_type, spin_zero, forbidden_spin_keys):
        """
        Tests
        -----
        - Initialization succeeds for supported event types.
        - ``gw_parameters_rvs`` returns expected base keys with requested size.
        - Spin/precession keys are absent when ``spin_zero=True``.
        """
        # Create a sampler with spins disabled for the chosen event type.
        cbc = _make_cbc(interpolator_dir, event_type=event_type, spin_zero=spin_zero)

        # Draw a small sample of parameters.
        params = cbc.gw_parameters_rvs(size=N_SAMPLES)

        # Validate output structure and basic invariants.
        assert isinstance(params, dict), f"gw_parameters_rvs: expected dict, got {type(params)}"
        self._assert_param_dict_valid(params, expected_keys=BASE_KEYS, size=N_SAMPLES)
        assert cbc.event_type == event_type, f"event_type: expected {event_type}, got {cbc.event_type}"
        assert cbc.spin_zero is True, f"spin_zero: expected True, got {cbc.spin_zero}"

        # Ensure spin-related keys are not present when spins are disabled.
        for spin_key in forbidden_spin_keys:
            assert spin_key not in params, f"spin key '{spin_key}' should be absent"


    def test_bbh_precessing_spin_output_sanity(self, interpolator_dir):
        """
        Tests
        -----
        - Precessing-spin setup initializes for BBH.
        - All precession-specific keys are present in ``gw_parameters_rvs`` output.
        - All output arrays have requested sample size and finite values.
        """
        # Create a BBH sampler with precessing spins enabled.
        cbc = _make_cbc(
            interpolator_dir,
            event_type="BBH",
            spin_zero=False,
            spin_precession=True,
        )

        # Draw a small sample.
        params = cbc.gw_parameters_rvs(size=N_SAMPLES)

        # Validate that all expected keys exist and values are finite.
        self._assert_param_dict_valid(
            params,
            expected_keys=BASE_KEYS + EXPECTED_SPIN_PRECESSING_KEYS,
            size=N_SAMPLES,
        )


    def test_custom_all_gw_priors_sanity(self, interpolator_dir):
        """
        Tests
        -----
        - User-defined callables are accepted for supported ``gw_priors`` hooks.
        - Both mass pathways are covered:
          1) BBH with custom ``mass_ratio``
          2) BNS with custom ``mass_2_source``
        - Returned outputs are structurally sane (keys, sizes, finite values).
        """
        # Use the shared config, but keep it local to this test.
        config = _make_config(interpolator_dir)

        # Define lightweight custom priors (pure numpy RNG) so we can validate
        # user-overrides without relying on internal prior choices.
        def custom_zs(size, *args, **kwargs):
            return np.random.uniform(0.01, 2.0, size)

        def custom_mass_1(size, *args, **kwargs):
            return np.random.uniform(20.0, 40.0, size)

        def custom_mass_ratio(size, *args, **kwargs):
            return np.random.uniform(0.4, 1.0, size)

        def custom_mass_2(size, *args, **kwargs):
            return np.random.uniform(1.1, 2.2, size)

        def custom_geocent_time(size, *args, **kwargs):
            return np.random.uniform(1238166018, 1269702018, size)

        def custom_ra(size, *args, **kwargs):
            return np.random.uniform(0.0, 2.0 * np.pi, size)

        def custom_dec(size, *args, **kwargs):
            return np.random.uniform(-np.pi / 2.0, np.pi / 2.0, size)

        def custom_phase(size, *args, **kwargs):
            return np.random.uniform(0.0, 2.0 * np.pi, size)

        def custom_psi(size, *args, **kwargs):
            return np.random.uniform(0.0, np.pi, size)

        def custom_theta_jn(size, *args, **kwargs):
            return np.random.uniform(0.0, np.pi, size)

        def custom_a_1(size, *args, **kwargs):
            return np.random.uniform(0.0, 0.9, size)

        def custom_a_2(size, *args, **kwargs):
            return np.random.uniform(0.0, 0.9, size)

        def custom_tilt_1(size, *args, **kwargs):
            return np.random.uniform(0.0, np.pi, size)

        def custom_tilt_2(size, *args, **kwargs):
            return np.random.uniform(0.0, np.pi, size)

        def custom_phi_12(size, *args, **kwargs):
            return np.random.uniform(0.0, 2.0 * np.pi, size)

        def custom_phi_jl(size, *args, **kwargs):
            return np.random.uniform(0.0, 2.0 * np.pi, size)

        # Path 1: BBH with custom mass_ratio (m2 derived from m1 * q).
        cbc_bbh = CBCSourceParameterDistribution(
            event_type="BBH",
            spin_zero=False,
            spin_precession=True,
            directory=config["directory"],
            gw_priors=dict(
                zs=custom_zs,
                mass_1_source=custom_mass_1,
                mass_ratio=custom_mass_ratio,
                geocent_time=custom_geocent_time,
                ra=custom_ra,
                dec=custom_dec,
                phase=custom_phase,
                psi=custom_psi,
                theta_jn=custom_theta_jn,
                a_1=custom_a_1,
                a_2=custom_a_2,
                tilt_1=custom_tilt_1,
                tilt_2=custom_tilt_2,
                phi_12=custom_phi_12,
                phi_jl=custom_phi_jl,
            ),
            z_min=config["z_min"],
            z_max=config["z_max"],
            create_new_interpolator=config["create_new_interpolator"],
        )

        # Draw samples and validate structure.
        params_bbh = cbc_bbh.gw_parameters_rvs(size=N_SAMPLES)

        self._assert_param_dict_valid(
            params_bbh,
            expected_keys=BASE_KEYS + EXPECTED_SPIN_PRECESSING_KEYS,
            size=N_SAMPLES,
        )
        assert "mass_ratio" in params_bbh, "mass_ratio key missing from BBH custom prior output"
        self._assert_array_valid(params_bbh["mass_ratio"], name="mass_ratio", size=N_SAMPLES)
        assert np.all((params_bbh["mass_1_source"] >= 20.0) & (params_bbh["mass_1_source"] <= 40.0)), \
            f"mass_1_source out of [20, 40]: min={params_bbh['mass_1_source'].min():.2f}, max={params_bbh['mass_1_source'].max():.2f}"
        assert np.all((params_bbh["mass_ratio"] >= 0.4) & (params_bbh["mass_ratio"] <= 1.0)), \
            f"mass_ratio out of [0.4, 1.0]: min={params_bbh['mass_ratio'].min():.2f}, max={params_bbh['mass_ratio'].max():.2f}"

        # Path 2: BNS with custom mass_2_source (explicit secondary-mass prior).
        cbc_bns = CBCSourceParameterDistribution(
            event_type="BNS",
            spin_zero=True,
            directory=config["directory"],
            gw_priors=dict(
                zs=custom_zs,
                mass_1_source=custom_mass_1,
                mass_2_source=custom_mass_2,
                geocent_time=custom_geocent_time,
                ra=custom_ra,
                dec=custom_dec,
                phase=custom_phase,
                psi=custom_psi,
                theta_jn=custom_theta_jn,
            ),
            z_min=config["z_min"],
            z_max=config["z_max"],
            create_new_interpolator=config["create_new_interpolator"],
        )

        # Draw samples and validate structure.
        params_bns = cbc_bns.gw_parameters_rvs(size=N_SAMPLES)

        self._assert_param_dict_valid(params_bns, expected_keys=BASE_KEYS, size=N_SAMPLES)
        assert "mass_ratio" not in params_bns, "mass_ratio should be absent for BNS with explicit mass_2_source prior"
        assert np.all((params_bns["mass_1_source"] >= 20.0) & (params_bns["mass_1_source"] <= 40.0)), \
            f"mass_1_source out of [20, 40]: min={params_bns['mass_1_source'].min():.2f}, max={params_bns['mass_1_source'].max():.2f}"
        assert np.all((params_bns["mass_2_source"] >= 1.1) & (params_bns["mass_2_source"] <= 2.2)), \
            f"mass_2_source out of [1.1, 2.2]: min={params_bns['mass_2_source'].min():.2f}, max={params_bns['mass_2_source'].max():.2f}"

    def test_custom_prior_parameters(self, interpolator_dir):
        """
        Tests
        -----
        - Custom parameters work for (i) ler-available prior samplers (string-based).
        - User-provided callables work for (ii) custom prior samplers.

        Notes
        -----
        This test intentionally exercises "parameter overrides" (`gw_priors_params`)
        across *all* sampled parameters, not only masses.
        """
        config = _make_config(interpolator_dir)

        # ------------------------------------------------------------------
        # Part (i): ler-available sampler + user-provided parameter override
        # ------------------------------------------------------------------
        zs_lo, zs_hi = 0.2, 0.25
        m1_lo, m1_hi = 1.85, 1.95
        m2_lo, m2_hi = 1.15, 1.25
        t_lo, t_hi = 1238166018.0, 1238166020.0
        cbc_builtin = CBCSourceParameterDistribution(
            event_type="BNS",
            spin_zero=True,
            directory=config["directory"],
            gw_priors=dict(
                mass_1_source="uniform",
                mass_2_source="uniform",
                geocent_time="uniform",
                ra="uniform",
                dec="uniform",
                phase="uniform",
                psi="uniform",
                theta_jn="uniform",
            ),
            gw_priors_params=dict(
                mass_1_source=dict(x_min=m1_lo, x_max=m1_hi),
                mass_2_source=dict(x_min=m2_lo, x_max=m2_hi),
                geocent_time=dict(x_min=t_lo, x_max=t_hi),
                # angular priors use the default support; we still run through
                # `gw_priors_params` to ensure overrides are accepted broadly.
                ra=dict(x_min=0.0, x_max=2.0 * np.pi),
                dec=dict(x_min=-np.pi / 2.0, x_max=np.pi / 2.0),
                phase=dict(x_min=0.0, x_max=2.0 * np.pi),
                psi=dict(x_min=0.0, x_max=np.pi),
                theta_jn=dict(x_min=0.0, x_max=np.pi),
            ),
            z_min=config["z_min"],
            z_max=config["z_max"],
            create_new_interpolator=config["create_new_interpolator"],
        )
        params_builtin = cbc_builtin.gw_parameters_rvs(size=N_SAMPLES)
        self._assert_param_dict_valid(params_builtin, expected_keys=BASE_KEYS, size=N_SAMPLES)
        assert "mass_ratio" not in params_builtin, "mass_ratio should be absent for BNS with uniform mass prior"
        # `zs` sampler is internally tied to the merger-rate-density based sampler;
        # here we only check it respects the configured z-range.
        assert np.all((params_builtin["zs"] >= config["z_min"]) & (params_builtin["zs"] <= config["z_max"])), \
            f"zs out of [{config['z_min']}, {config['z_max']}]: min={params_builtin['zs'].min():.3f}, max={params_builtin['zs'].max():.3f}"
        assert np.all((params_builtin["mass_1_source"] >= m1_lo) & (params_builtin["mass_1_source"] <= m1_hi)), \
            f"mass_1_source out of [{m1_lo}, {m1_hi}]: min={params_builtin['mass_1_source'].min():.3f}, max={params_builtin['mass_1_source'].max():.3f}"
        assert np.all((params_builtin["mass_2_source"] >= m2_lo) & (params_builtin["mass_2_source"] <= m2_hi)), \
            f"mass_2_source out of [{m2_lo}, {m2_hi}]: min={params_builtin['mass_2_source'].min():.3f}, max={params_builtin['mass_2_source'].max():.3f}"
        assert np.all((params_builtin["geocent_time"] >= t_lo) & (params_builtin["geocent_time"] <= t_hi)), \
            f"geocent_time out of [{t_lo}, {t_hi}]: min={params_builtin['geocent_time'].min():.3f}, max={params_builtin['geocent_time'].max():.3f}"
        assert np.all((params_builtin["ra"] >= 0.0) & (params_builtin["ra"] <= 2.0 * np.pi)), \
            f"ra out of [0, 2pi]: min={params_builtin['ra'].min():.3f}, max={params_builtin['ra'].max():.3f}"
        assert np.all((params_builtin["dec"] >= -np.pi / 2.0) & (params_builtin["dec"] <= np.pi / 2.0)), \
            f"dec out of [-pi/2, pi/2]: min={params_builtin['dec'].min():.3f}, max={params_builtin['dec'].max():.3f}"
        assert np.all((params_builtin["phase"] >= 0.0) & (params_builtin["phase"] <= 2.0 * np.pi)), \
            f"phase out of [0, 2pi]: min={params_builtin['phase'].min():.3f}, max={params_builtin['phase'].max():.3f}"
        assert np.all((params_builtin["psi"] >= 0.0) & (params_builtin["psi"] <= np.pi)), \
            f"psi out of [0, pi]: min={params_builtin['psi'].min():.3f}, max={params_builtin['psi'].max():.3f}"
        assert np.all((params_builtin["theta_jn"] >= 0.0) & (params_builtin["theta_jn"] <= np.pi)), \
            f"theta_jn out of [0, pi]: min={params_builtin['theta_jn'].min():.3f}, max={params_builtin['theta_jn'].max():.3f}"

        # ------------------------------------------------------------------
        # Part (ii): user-provided callable
        # ------------------------------------------------------------------
        def custom_zs(size, *args, **kwargs):
            return np.random.uniform(zs_lo, zs_hi, size)

        def custom_mass_1(size, *args, **kwargs):
            return np.random.uniform(10.0, 12.0, size)

        def custom_mass_2(size, *args, **kwargs):
            return np.random.uniform(1.2, 1.3, size)

        def custom_geocent_time(size, *args, **kwargs):
            return np.random.uniform(1238166025.0, 1238166026.0, size)

        def custom_ra(size, *args, **kwargs):
            return np.random.uniform(0.1, 0.2, size)

        def custom_dec(size, *args, **kwargs):
            return np.random.uniform(-0.2, -0.1, size)

        def custom_phase(size, *args, **kwargs):
            return np.random.uniform(0.2, 0.3, size)

        def custom_psi(size, *args, **kwargs):
            return np.random.uniform(0.3, 0.4, size)

        def custom_theta_jn(size, *args, **kwargs):
            return np.random.uniform(0.4, 0.5, size)

        cbc_custom = CBCSourceParameterDistribution(
            event_type="BNS",
            spin_zero=True,
            directory=config["directory"],
            gw_priors=dict(
                zs=custom_zs,
                mass_1_source=custom_mass_1,
                mass_2_source=custom_mass_2,
                geocent_time=custom_geocent_time,
                ra=custom_ra,
                dec=custom_dec,
                phase=custom_phase,
                psi=custom_psi,
                theta_jn=custom_theta_jn,
            ),
            z_min=config["z_min"],
            z_max=config["z_max"],
            create_new_interpolator=config["create_new_interpolator"],
        )
        params_custom = cbc_custom.gw_parameters_rvs(size=N_SAMPLES)
        self._assert_param_dict_valid(params_custom, expected_keys=BASE_KEYS, size=N_SAMPLES)
        assert "mass_ratio" not in params_custom, "mass_ratio should be absent for BNS with custom mass_2_source prior"
        assert np.all((params_custom["zs"] >= zs_lo) & (params_custom["zs"] <= zs_hi)), \
            f"zs out of [{zs_lo}, {zs_hi}]: min={params_custom['zs'].min():.3f}, max={params_custom['zs'].max():.3f}"
        assert np.all((params_custom["mass_1_source"] >= 10.0) & (params_custom["mass_1_source"] <= 12.0)), \
            f"mass_1_source out of [10, 12]: min={params_custom['mass_1_source'].min():.2f}, max={params_custom['mass_1_source'].max():.2f}"
        assert np.all((params_custom["mass_2_source"] >= 1.2) & (params_custom["mass_2_source"] <= 1.3)), \
            f"mass_2_source out of [1.2, 1.3]: min={params_custom['mass_2_source'].min():.3f}, max={params_custom['mass_2_source'].max():.3f}"
        assert np.all((params_custom["geocent_time"] >= 1238166025.0) & (params_custom["geocent_time"] <= 1238166026.0)), \
            f"geocent_time out of expected range: min={params_custom['geocent_time'].min():.2f}, max={params_custom['geocent_time'].max():.2f}"
        assert np.all((params_custom["ra"] >= 0.1) & (params_custom["ra"] <= 0.2)), \
            f"ra out of [0.1, 0.2]: min={params_custom['ra'].min():.3f}, max={params_custom['ra'].max():.3f}"
        assert np.all((params_custom["dec"] >= -0.2) & (params_custom["dec"] <= -0.1)), \
            f"dec out of [-0.2, -0.1]: min={params_custom['dec'].min():.3f}, max={params_custom['dec'].max():.3f}"
        assert np.all((params_custom["phase"] >= 0.2) & (params_custom["phase"] <= 0.3)), \
            f"phase out of [0.2, 0.3]: min={params_custom['phase'].min():.3f}, max={params_custom['phase'].max():.3f}"
        assert np.all((params_custom["psi"] >= 0.3) & (params_custom["psi"] <= 0.4)), \
            f"psi out of [0.3, 0.4]: min={params_custom['psi'].min():.3f}, max={params_custom['psi'].max():.3f}"
        assert np.all((params_custom["theta_jn"] >= 0.4) & (params_custom["theta_jn"] <= 0.5)), \
            f"theta_jn out of [0.4, 0.5]: min={params_custom['theta_jn'].min():.3f}, max={params_custom['theta_jn'].max():.3f}"

    def test_sampling_reproducibility_with_reseeding(self, interpolator_dir):
        """
        Tests
        -----
        - Re-seeding NumPy with the same seed reproduces identical samples
          for custom ``gw_priors``.
        """
        # Use a local config and custom priors that depend only on NumPy RNG.
        config = _make_config(interpolator_dir)

        def custom_zs(size, *args, **kwargs):
            return np.random.uniform(0.01, 2.0, size)

        def custom_mass_1(size, *args, **kwargs):
            return np.random.uniform(20.0, 40.0, size)

        def custom_mass_ratio(size, *args, **kwargs):
            return np.random.uniform(0.4, 1.0, size)

        def custom_geocent_time(size, *args, **kwargs):
            return np.random.uniform(1238166018, 1269702018, size)

        def custom_ra(size, *args, **kwargs):
            return np.random.uniform(0.0, 2.0 * np.pi, size)

        def custom_dec(size, *args, **kwargs):
            return np.random.uniform(-np.pi / 2.0, np.pi / 2.0, size)

        def custom_phase(size, *args, **kwargs):
            return np.random.uniform(0.0, 2.0 * np.pi, size)

        def custom_psi(size, *args, **kwargs):
            return np.random.uniform(0.0, np.pi, size)

        def custom_theta_jn(size, *args, **kwargs):
            return np.random.uniform(0.0, np.pi, size)

        cbc = CBCSourceParameterDistribution(
            event_type="BBH",
            spin_zero=True,
            directory=config["directory"],
            gw_priors=dict(
                zs=custom_zs,
                mass_1_source=custom_mass_1,
                mass_ratio=custom_mass_ratio,
                geocent_time=custom_geocent_time,
                ra=custom_ra,
                dec=custom_dec,
                phase=custom_phase,
                psi=custom_psi,
                theta_jn=custom_theta_jn,
            ),
            z_min=config["z_min"],
            z_max=config["z_max"],
            create_new_interpolator=config["create_new_interpolator"],
        )

        # Draw once, then reseed and draw again; outputs must match exactly.
        np.random.seed(2026)
        params_1 = cbc.gw_parameters_rvs(size=N_SAMPLES)

        np.random.seed(2026)
        params_2 = cbc.gw_parameters_rvs(size=N_SAMPLES)

        assert params_1.keys() == params_2.keys(), \
            f"keys differ between reseeded draws: {params_1.keys()} vs {params_2.keys()}"
        for key in params_1:
            np.testing.assert_array_equal(params_1[key], params_2[key])

    def test_sample_gw_parameters_fixed_param_override(self, interpolator_dir):
        """
        Tests
        -----
        - ``sample_gw_parameters`` honors fixed-parameter overrides via ``param``.

        Notes
        -----
        Internally, ``sample_gw_parameters`` forwards ``param`` to the configured
        sampling routine as keyword arguments, so values provided in ``param``
        should be returned unchanged.
        """
        config = _make_config(interpolator_dir)
        cbc = CBCSourceParameterDistribution(
            event_type="BNS",
            spin_zero=True,
            directory=config["directory"],
            z_min=config["z_min"],
            z_max=config["z_max"],
            create_new_interpolator=config["create_new_interpolator"],
        )

        size = 20
        zs0 = 0.123
        m10 = 1.9
        m20 = 1.2
        ra0 = 0.42
        dec0 = -0.17
        phase0 = 0.9
        psi0 = 0.3
        theta0 = 1.1
        t0 = 1238166019.0

        fixed = dict(
            zs=np.full(size, zs0),
            mass_1_source=np.full(size, m10),
            mass_2_source=np.full(size, m20),
            geocent_time=np.full(size, t0),
            ra=np.full(size, ra0),
            dec=np.full(size, dec0),
            phase=np.full(size, phase0),
            psi=np.full(size, psi0),
            theta_jn=np.full(size, theta0),
        )

        params = cbc.sample_gw_parameters(size=size, param=fixed)
        self._assert_param_dict_valid(params, expected_keys=BASE_KEYS, size=size)
        assert "mass_ratio" not in params, "mass_ratio should be absent for BNS with explicit mass_2_source prior"

        for k, v in fixed.items():
            np.testing.assert_array_equal(params[k], v)

        # Derived quantities should still be consistent with the fixed inputs.
        np.testing.assert_allclose(params["mass_1"], m10 * (1.0 + zs0), rtol=1e-12)
        np.testing.assert_allclose(params["mass_2"], m20 * (1.0 + zs0), rtol=1e-12)
        np.testing.assert_allclose(
            params["luminosity_distance"],
            cbc.luminosity_distance(np.full(size, zs0)),
            rtol=1e-10,
        )

    def test_invalid_event_type_raises(self, interpolator_dir):
        """
        Tests
        -----
        - Unknown ``event_type`` raises ``ValueError``.
        """
        # Ensure invalid event_type fails fast with a ValueError.
        config = _make_config(interpolator_dir)
        with pytest.raises(ValueError, match="event_type"):
            CBCSourceParameterDistribution(
                event_type="UNKNOWN",
                spin_zero=True,
                directory=config["directory"],
                z_min=config["z_min"],
                z_max=config["z_max"],
                create_new_interpolator=config["create_new_interpolator"],
            )

    def test_gw_parameters_rvs_njit_output_sanity(self, interpolator_dir):
        """
        Tests
        -----
        - ``gw_parameters_rvs_njit`` returns expected base keys and valid
          finite samples for BBH with ``spin_zero=True`` (no spin keys).
        """
        cbc = _make_cbc(
            interpolator_dir,
            event_type="BBH",
            spin_zero=True,
        )
        params = cbc.gw_parameters_rvs_njit(size=N_SAMPLES)
        assert isinstance(params, dict), f"gw_parameters_rvs_njit: expected dict, got {type(params)}"
        self._assert_param_dict_valid(params, expected_keys=BASE_KEYS, size=N_SAMPLES)
        for k in EXPECTED_SPIN_PRECESSING_KEYS:
            assert k not in params, f"spin key '{k}' should be absent when spin_zero=True"

    @pytest.mark.slow
    def test_njit_speed_gw_parameters_rvs(self, interpolator_dir):
        """
        Tests
        -----
        - Compare wall time of ``gw_parameters_rvs_njit`` in three modes:
          i) njit with parallel ``npool`` (``clamp_npool_for_numba(6)``, in-process JIT)
          ii) njit with ``npool=1`` (in-process, JIT enabled)
          iii) no-JIT baseline: subprocess with ``NUMBA_DISABLE_JIT=1`` and ``npool=1``
        - The no-JIT baseline runs in a subprocess so the env var is set before
          numba is first imported.
        - The njit paths are warmed up before timing; the timed section is
          repeated to reduce noise.

        Notes
        -----
        Coarse speed sanity check (not a strict benchmark). This does not compare
        ``gw_parameters_rvs`` to ``gw_parameters_rvs_njit``; correctness of the
        njit entry point is covered in ``test_gw_parameters_rvs_njit_output_sanity``.
        """
        import os
        import subprocess
        import sys

        sample_size = 50000
        repeats = 3
        min_speedup = 1.01

        # ------------------------------------------------------------------
        # Test 1: no-JIT baseline (subprocess)
        # ------------------------------------------------------------------
        script_no_jit = "\n".join(
            [
                "import time",
                "import numpy as np",
                "from ler.gw_source_population import CBCSourceParameterDistribution",
                "cbc = CBCSourceParameterDistribution(",
                "    npool=1, event_type='BBH', spin_zero=True,",
                "    z_min=0.0, z_max=10.0,",
                "    create_new_interpolator=False,",
                f"    directory={repr(interpolator_dir)},",
                ")",
                "_ = cbc.gw_parameters_rvs_njit(2)",
                "times = []",
                f"sample_size = {sample_size}",
                f"repeats = {repeats}",
                "for _ in range(repeats):",
                "    t0 = time.perf_counter()",
                "    _ = cbc.gw_parameters_rvs_njit(sample_size)",
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
        # Test 2/3: njit-backed sampler (in-process), npool=1 and parallel npool
        # ------------------------------------------------------------------
        def _time_gw_parameters_rvs_njit(npool):
            cfg = _make_config(interpolator_dir, npool=npool)
            cbc = CBCSourceParameterDistribution(
                event_type="BBH",
                spin_zero=True,
                **cfg,
            )
            _ = cbc.gw_parameters_rvs_njit(2)
            return median_call_time(
                lambda: cbc.gw_parameters_rvs_njit(sample_size),
                repeats=repeats,
            )

        t_njit_1 = _time_gw_parameters_rvs_njit(1)
        t_njit_hi = _time_gw_parameters_rvs_njit(NPOOL_PARALLEL)

        assert t_no_jit > 0.0 and t_njit_1 > 0.0 and t_njit_hi > 0.0, (
            f"invalid timings: t_no_jit={t_no_jit}, t_njit_1={t_njit_1}, "
            f"t_njit_hi={t_njit_hi} (npool={NPOOL_PARALLEL})"
        )

        speedup_1 = t_no_jit / t_njit_1
        speedup_hi = t_no_jit / t_njit_hi

        assert np.isfinite(speedup_1) and np.isfinite(speedup_hi), (
            f"invalid speedups: speedup_1={speedup_1}, speedup_hi={speedup_hi}"
        )
        assert speedup_1 >= min_speedup, (
            f"Expected njit (npool=1) to be faster after warm-up, but got speedup={speedup_1:.2f}x "
            f"(no-JIT={t_no_jit:.6f}s, njit npool=1={t_njit_1:.6f}s)."
        )
        assert speedup_hi >= min_speedup, (
            f"Expected njit (npool={NPOOL_PARALLEL}) to be faster after warm-up, "
            f"but got speedup={speedup_hi:.2f}x "
            f"(no-JIT={t_no_jit:.6f}s, njit npool={NPOOL_PARALLEL} wall={t_njit_hi:.6f}s)."
        )
