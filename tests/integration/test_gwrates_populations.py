"""
Integration tests for ``GWRATES`` (population sampling and BBH hyperparameter stress).

Test Coverage:
--------------
- ``test_gwrates_population_rates_and_outputs``: BBH / BNS / NSBH;
  ``gw_cbc_statistics`` at ``SIZE``; parameter dict keys and finite arrays;
  ``gw_rate`` above a loose floor and **below** ``normalization_pdf_z``
  (intrinsic merger-rate normalization from the CBC model; ``gw_rate`` is that
  rate times the simulated detectable fraction).
- ``test_gwrates_hyper_uncertainty_rates``: BBH only; vary one of ``H0``,
  merger ``R0``, or ``mass_1_source`` ``alpha_1`` between two values; compares
  rates so the endpoint with larger astrophysical rate gives a larger number
  (``normalization_pdf_z`` for the ``H0`` case, ``gw_rate`` for ``R0`` /
  ``alpha_1``).

They use large Monte Carlo samples (``SIZE``, ``BATCH_SIZE`` below) and real
``gwsnr`` detection, and run in default CI (``.github/workflows/tests.yml``).

"""

import os

import numpy as np
import pytest
from astropy.cosmology import LambdaCDM

from ler.rates import GWRATES
from tests_utils import (
    CommonTestUtils,
    EXPECTED_UNLENSED_GW_KEYS_NO_SPIN,
    EXPECTED_UNLENSED_GW_KEYS_PRECESSING_SPIN,
    clamp_npool_for_numba,
)

# Prefer six workers where Numba permits; clamps on four-vCPU GitHub-hosted runners.
NPOOL = clamp_npool_for_numba(6)

SIZE = 20000
BATCH_SIZE = 10000

Z_MIN = 0.0
Z_MAX = 10.0

COSMO_OM0 = 0.315
COSMO_ODE0 = 0.685
COSMO_TCMB = 2.725

# Planck 2018-like mean with a ±0.42 km/s/Mpc bracket for the H0-only check.
H0_LOW = 67.66 - 0.42
H0_HIGH = 67.66 + 0.42

# Merger-rate amplitude R0 (Gpc^-3 yr^-1): low vs high vs BBH default ~19e-9.
R0_LOW = 14.0 * 1e-9
R0_HIGH = 26.0 * 1e-9

# Broken-power-law slope below break for ``mass_1_source`` (BBH only in these tests).
ALPHA_1_LOW = -0.05
ALPHA_1_HIGH = 2.86


def _make_cosmo(h0):
    """Flat ``LambdaCDM`` with fixed densities; only ``H0`` varies."""
    return LambdaCDM(
        H0=h0,
        Om0=COSMO_OM0,
        Ode0=COSMO_ODE0,
        Tcmb0=COSMO_TCMB,
    )


class TestGWRATESPopulationsUncertainty(CommonTestUtils):
    """Integration tests for ``GWRATES`` with real ``gwsnr`` detection."""

    @pytest.mark.parametrize(
        "event_type, rate_low, expected_gw_keys",
        [
            ("BBH", 1.0e2, EXPECTED_UNLENSED_GW_KEYS_PRECESSING_SPIN),
            ("BNS", 1e5, EXPECTED_UNLENSED_GW_KEYS_NO_SPIN),
            ("NSBH", 1e4, EXPECTED_UNLENSED_GW_KEYS_NO_SPIN),
        ],
    )
    def test_gwrates_population_rates_and_outputs(
        self,
        interpolator_directory,
        ler_directory,
        event_type,
        rate_low,
        expected_gw_keys,
    ):
        """
        Tests
        -----
        - ``gw_cbc_statistics`` returns dict-shaped output with ``SIZE`` samples;
          required keys present; arrays finite (via ``_assert_param_dict_valid``).
        - ``normalization_pdf_z`` is the intrinsic merger-rate scale (yr^-1) for
          the CBC + cosmology setup (total rate before the detection cut).
        - ``gw_rate`` equals ``normalization_pdf_z * (n_det / n_total)`` in
          default boolean-detection mode, so ``gw_rate < normalization_pdf_z``;
          also ``gw_rate > rate_low`` (parametrized loose floor).
        """
        ler_dir = os.path.join(ler_directory, f"pop_{event_type}")
        os.makedirs(ler_dir, exist_ok=True)

        if event_type == "BBH":
            ler = GWRATES(
                npool=NPOOL,
                z_min=Z_MIN,
                z_max=Z_MAX,
                event_type=event_type,
                interpolator_directory=interpolator_directory,
                ler_directory=ler_dir,
                spin_zero=False,
                spin_precession=True,
                waveform_approximant="IMRPhenomXPHM",
                snr_recalculation=True,
                snr_recalculation_range=[6, 14],
                snr_recalculation_waveform_approximant="IMRPhenomXPHM",
            )
        elif event_type == "BNS":
            ler = GWRATES(
                npool=NPOOL,
                z_min=Z_MIN,
                z_max=Z_MAX,
                event_type=event_type,
                interpolator_directory=interpolator_directory,
                ler_directory=ler_dir,
                snr_method="interpolation_no_spins",
                ifos=["ET", "CE"],
                snr_type="optimal_snr",
                spin_zero=True,
                spin_precession=False,
                waveform_approximant="IMRPhenomD",
            )
        elif event_type == "NSBH":
            ler = GWRATES(
                npool=NPOOL,
                z_min=Z_MIN,
                z_max=Z_MAX,
                event_type=event_type,
                interpolator_directory=interpolator_directory,
                ler_directory=ler_dir,
                snr_method="interpolation_no_spins",
                ifos=["ET", "CE"],
                snr_type="optimal_snr",
                spin_zero=True,
                spin_precession=False,
                waveform_approximant="IMRPhenomD",
            )
        else:
            raise ValueError(f"unknown event_type: {event_type}")

        param = ler.gw_cbc_statistics(
            size=SIZE,
            batch_size=BATCH_SIZE,
            resume=False,
            save_batch=False,
            output_jsonfile=False,
        )

        self._assert_unlensed_cbc_outputs(
            param, SIZE, expected_gw_keys=expected_gw_keys,
        )

        rate, _ = ler.gw_rate(gw_param=param, output_jsonfile=False)
        cap = ler.normalization_pdf_z

        print(
            f"\n{ler.event_type}: gw_rate={rate} yr^-1, "
            f"floor={rate_low} yr^-1, norm={cap} yr^-1\n"
        )

        assert rate > rate_low, (
            f"{ler.event_type}: rate {rate} yr^-1 below floor {rate_low} yr^-1"
        )
        assert rate < cap, (
            f"{ler.event_type}: rate {rate} yr^-1 above norm {cap} yr^-1"
        )

    @pytest.mark.parametrize(
        "uncertainty_type, val_a, val_b",
        [
            # Intrinsic coalescence rate in the volume integral drops when H0 increases.
            ("H0", H0_HIGH, H0_LOW),
            # Higher R0 raises the detected rate (same network cut).
            ("R0", R0_LOW, R0_HIGH),
            # Higher alpha_1 de-weights high mass_1; rate drops.
            ("alpha_1", ALPHA_1_HIGH, ALPHA_1_LOW),
        ],
    )
    def test_gwrates_hyper_uncertainty_rates(
        self,
        interpolator_directory,
        ler_directory,
        uncertainty_type,
        val_a,
        val_b,
    ):
        """
        Tests
        -----
        - Two BBH ``GWRATES`` instances differ in exactly one hyperparameter
          (``H0``, merger ``R0``, or ``mass_1_source`` ``alpha_1``); outputs
          isolated under ``ler_directory/hyper_*``.
        - Parametrization orders ``val_a`` / ``val_b`` so the larger physical
          rate is at ``val_b`` (see comments on each parametrized row).
        - For ``H0`` only: compare ``normalization_pdf_z`` (intrinsic rate).
        - For ``R0`` and ``alpha_1``: compare ``gw_rate`` after
          ``gw_cbc_statistics`` (detection-weighted rate).
        - Assert ``rate_b > rate_a`` and a positive half-width between them.
        """
        ler_a = os.path.join(ler_directory, f"hyper_{uncertainty_type}_a")
        ler_b = os.path.join(ler_directory, f"hyper_{uncertainty_type}_b")
        os.makedirs(ler_a, exist_ok=True)
        os.makedirs(ler_b, exist_ok=True)

        if uncertainty_type == "H0":
            g_a = GWRATES(
                npool=NPOOL,
                z_min=Z_MIN,
                z_max=Z_MAX,
                event_type="BBH",
                interpolator_directory=interpolator_directory,
                ler_directory=ler_a,
                cosmology=_make_cosmo(val_a),
            )
            g_b = GWRATES(
                npool=NPOOL,
                z_min=Z_MIN,
                z_max=Z_MAX,
                event_type="BBH",
                interpolator_directory=interpolator_directory,
                ler_directory=ler_b,
                cosmology=_make_cosmo(val_b),
            )
        elif uncertainty_type == "R0":
            g_a = GWRATES(
                npool=NPOOL,
                z_min=Z_MIN,
                z_max=Z_MAX,
                event_type="BBH",
                interpolator_directory=interpolator_directory,
                ler_directory=ler_a,
                gw_functions_params={"merger_rate_density": {"R0": val_a}},
            )
            g_b = GWRATES(
                npool=NPOOL,
                z_min=Z_MIN,
                z_max=Z_MAX,
                event_type="BBH",
                interpolator_directory=interpolator_directory,
                ler_directory=ler_b,
                gw_functions_params={"merger_rate_density": {"R0": val_b}},
            )
        elif uncertainty_type == "alpha_1":
            g_a = GWRATES(
                npool=NPOOL,
                z_min=Z_MIN,
                z_max=Z_MAX,
                event_type="BBH",
                interpolator_directory=interpolator_directory,
                ler_directory=ler_a,
                gw_priors_params={"mass_1_source": {"alpha_1": val_a}},
            )
            g_b = GWRATES(
                npool=NPOOL,
                z_min=Z_MIN,
                z_max=Z_MAX,
                event_type="BBH",
                interpolator_directory=interpolator_directory,
                ler_directory=ler_b,
                gw_priors_params={"mass_1_source": {"alpha_1": val_b}},
            )
        else:
            raise ValueError(f"unknown uncertainty_type: {uncertainty_type}")

        param_a = g_a.gw_cbc_statistics(
            size=SIZE,
            batch_size=BATCH_SIZE,
            resume=False,
            save_batch=False,
            output_jsonfile=False,
        )
        param_b = g_b.gw_cbc_statistics(
            size=SIZE,
            batch_size=BATCH_SIZE,
            resume=False,
            save_batch=False,
            output_jsonfile=False,
        )

        if uncertainty_type == "H0":
            rate_a = float(g_a.normalization_pdf_z)
            rate_b = float(g_b.normalization_pdf_z)
        else:
            rate_a, _ = g_a.gw_rate(gw_param=param_a, output_jsonfile=False)
            rate_b, _ = g_b.gw_rate(gw_param=param_b, output_jsonfile=False)

        print(
            f"\n{uncertainty_type}: val_a={val_a}, val_b={val_b} -> "
            f"rate_a={rate_a}, rate_b={rate_b}\n"
        )

        assert np.isfinite(rate_a), f"{uncertainty_type}: rate_a not finite"
        assert np.isfinite(rate_b), f"{uncertainty_type}: rate_b not finite"
        assert rate_b > rate_a, (
            f"{uncertainty_type}: expected rate_b > rate_a, got {rate_b} vs {rate_a}"
        )

        half = 0.5 * abs(rate_b - rate_a)
        assert np.isfinite(half) and half > 0, (
            f"{uncertainty_type}: rate bracket half-width invalid"
        )
