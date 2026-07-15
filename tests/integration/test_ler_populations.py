"""
Integration tests for ``LeR`` (lensed / unlensed populations and lens-galaxy ``phistar``).

Test Coverage:
--------------
- ``test_ler_population_rates_and_outputs``: ``lens_type`` in
  ``epl_shear_galaxy``, ``sie_galaxy``, ``sis_galaxy``; large-sample
  ``unlensed_cbc_statistics`` and ``lensed_cbc_statistics``; dict keys from
  ``tests_utils`` (including image columns and ``pdet_net``); ``unlensed_rate``,
  ``lensed_rate``, and ratio ``unlensed / lensed`` with loose floors vs
  ``normalization_pdf_z`` / ``normalization_pdf_z_lensed``.
- ``test_ler_hyper_uncertainty_rates``: two ``LeR`` instances with
  ``lens_type='sie_galaxy'`` differ only in ``phistar`` (velocity-dispersion
  prior, comoving number density normalization). Asserts higher ``phistar``
  gives a larger *detection-weighted* ``lensed_rate``, and
  ``unlensed_rate / lensed_rate`` is larger for the low-``phistar`` run.

Runs in default CI (``.github/workflows/tests.yml``).

"""

import os

import numpy as np
import pytest

from ler.rates import LeR
from tests_utils import (
    CommonTestUtils,
    EXPECTED_IMAGE_KEYS,
    EXPECTED_LENSED_PARAM_KEYS_ALIGNED_SPIN,
    EXPECTED_SOURCE_POS_KEYS,
    EXPECTED_UNLENSED_GW_KEYS_ALIGNED_SPIN,
    clamp_npool_for_numba,
)

# Monte Carlo sample size (same for unlensed and lensed draws).
SIZE = 20000
# Batch size for ``*_cbc_statistics`` I/O and memory.
BATCH_SIZE = 10000

# Worker count for multiprocessing inside LeR.
NPOOL = clamp_npool_for_numba(6)

# Source redshift range for CBC sampling and optical-depth integrals.
Z_MIN = 0.0
Z_MAX = 10.0

# Loose floors (yr^-1) so rates are not near zero from noise; upper bounds use
# ``normalization_pdf_z`` / ``normalization_pdf_z_lensed`` from the instance.
RATE_LOW_UNLENSED = 1.0e2
RATE_LOW_LENSED = 1.0e-2

# ``phistar`` in ``velocity_dispersion_ewoud`` (h^3 Mpc^-3): normalization of the
# lens galaxy comoving number density in that prior. Bracket two literature
# values so optical depth and lensed rates differ between runs.
PHISTAR_LOW = 0.00274  # Choi et al. (2005)
PHISTAR_HIGH = 2.099e-2  # Bernardi et al. (2010)


class TestLeRPopulations(CommonTestUtils):
    """Integration tests for ``LeR`` using real ``gwsnr`` detection."""

    @pytest.mark.parametrize(
        "lens_type",
        ["epl_shear_galaxy", "sie_galaxy", "sis_galaxy"],
    )
    def test_ler_population_rates_and_outputs(
        self,
        interpolator_directory,
        ler_directory,
        lens_type,
    ):
        """
        Tests
        -----
        - Build ``LeR`` with the parametrized ``lens_type``, ``event_type='BBH'``,
          and default BBH GW settings; write under ``ler_directory/pop_<lens_type>/``.
        - Run ``unlensed_cbc_statistics`` and ``lensed_cbc_statistics`` at ``SIZE``.
        - Check dict keys and array shapes (unlensed GW + ``pdet_net``; lensed 1-D,
          2-D image block, source positions, 2-D ``pdet_net``).
        - Compute ``unlensed_rate`` and ``lensed_rate`` from those dicts; require
          finite rates, above loose floors, below the corresponding pdf
          normalizations, and ``unlensed / lensed`` in a sensible range.
        """
        # Isolated output directory per lens model so runs do not clobber each other.
        ler_dir = os.path.join(ler_directory, f"pop_{lens_type}")
        os.makedirs(ler_dir, exist_ok=True)

        # Construct LeR with bundled interpolators and BBH defaults (aligned spins).
        ler = LeR(
            npool=NPOOL,
            z_min=Z_MIN,
            z_max=Z_MAX,
            event_type="BBH",
            lens_type=lens_type,
            interpolator_directory=interpolator_directory,
            ler_directory=ler_dir,
            create_new_interpolator=False,
        )

        # Draw intrinsic + detected CBC parameters (unlensed pipeline).
        unlensed = ler.unlensed_cbc_statistics(
            size=SIZE,
            batch_size=BATCH_SIZE,
            resume=False,
            save_batch=False,
            output_jsonfile=False,
        )
        # Draw lens + GW parameters for strongly lensed systems.
        lensed = ler.lensed_cbc_statistics(
            size=SIZE,
            batch_size=BATCH_SIZE,
            resume=False,
            save_batch=False,
            output_jsonfile=False,
        )

        # Unlensed: expected GW columns for aligned-spin BBH + scalar ``pdet_net``.
        self._assert_unlensed_cbc_outputs(
            unlensed,
            SIZE,
            expected_gw_keys=EXPECTED_UNLENSED_GW_KEYS_ALIGNED_SPIN,
        )
        # Lensed: 1-D keys including spin magnitudes (same BBH setup).
        self._assert_param_dict_valid(
            lensed,
            expected_keys=EXPECTED_LENSED_PARAM_KEYS_ALIGNED_SPIN,
            size=SIZE,
        )
        # Lensed: per-image arrays are 2-D with NaN padding where there is no image.
        shape2 = (SIZE, ler.n_max_images)
        self._assert_param_dict_valid(
            lensed,
            expected_keys=EXPECTED_IMAGE_KEYS,
            size=shape2,
            nan_to_num=True,
        )
        # Lensed: source position in the lens plane (1-D per event).
        self._assert_param_dict_valid(
            lensed,
            expected_keys=EXPECTED_SOURCE_POS_KEYS,
            size=SIZE,
        )
        # Lensed: network detection probability per image row (NaNs allowed on padding).
        self._assert_lensed_pdet_net_2d(lensed, shape2)

        # Detection-weighted rates from the Monte Carlo sample (boolean ``pdet`` cut).
        r_unl, _ = ler.unlensed_rate(
            unlensed_param=unlensed, output_jsonfile=False,
        )
        r_lens, _ = ler.lensed_rate(
            lensed_param=lensed, output_jsonfile=False,
        )
        # Typical expectation: many more unlensed detections than lensed per year.
        ratio = r_unl / r_lens

        print(
            f"\n{lens_type}: unlensed={r_unl} yr^-1, lensed={r_lens} yr^-1, "
            f"ratio={ratio}\n"
        )

        # Upper limits from the merger-rate × selection normalization constants.
        rate_max_unl = ler.normalization_pdf_z
        rate_max_lens = ler.normalization_pdf_z_lensed

        # Unlensed rate must be physical and strictly below the pdf normalization.
        assert np.isfinite(r_unl) and r_unl > RATE_LOW_UNLENSED and r_unl < rate_max_unl, (
            f"{lens_type}: unlensed_rate"
        )
        # Lensed rate same check against lensed normalization.
        assert np.isfinite(r_lens) and r_lens > RATE_LOW_LENSED and r_lens < rate_max_lens, (
            f"{lens_type}: lensed_rate"
        )
        # Ratio should be finite and above 1 (unlensed dominates); cap avoids bogus blow-ups.
        assert np.isfinite(ratio) and ratio > 1 and ratio < 1e4, f"{lens_type}: rate ratio"

    def test_ler_hyper_uncertainty_rates(
        self,
        interpolator_directory,
        ler_directory,
    ):
        """
        Tests
        -----
        - Build two ``LeR`` instances that only differ in ``phistar`` for the
          velocity-dispersion prior (``velocity_dispersion_ewoud``). Both use
          ``lens_type='sie_galaxy'`` and default BBH settings otherwise.
        - For each instance, run ``lensed_cbc_statistics`` at ``SIZE`` and
          ``lensed_rate``; compare high vs low ``phistar``.
        - Run ``unlensed_cbc_statistics`` once on ``ler_lo`` and
          ``unlensed_rate`` on that dict.  Form
          ``ratio_* = unlensed_rate / lensed_rate_*`` using the same numerator
          for both ``phistar`` values so only the lensed denominator changes.

        Expectations:
            - ``lensed_rate`` increases with ``phistar`` (more lens galaxies →
              higher detection-weighted lensed rate for this setup).
            - ``unlensed_rate / lensed_rate`` is smaller for high ``phistar``
              (lensed channel picks up relatively more events).
        """

        def velocity_dispersion_prior(phistar):
            """Keyword args for ``lens_priors_params['velocity_dispersion']`` (SIE path)."""
            return {
                "param_name": "velocity_dispersion",
                "sampler_type": "velocity_dispersion_ewoud",
                "sigma_min": 100.0,
                "sigma_max": 400.0,
                "alpha": 0.94,
                "beta": 1.85,
                "phistar": float(phistar),
                "sigmastar": 113.78,
            }

        def make_ler_for_phistar(subdir, phistar):
            """Create LeR under ``ler_directory/<subdir>/`` with the given ``phistar`` only."""
            # Subdirectory keeps JSON and scratch separate for the two hyperparameter values.
            sub = os.path.join(ler_directory, subdir)
            os.makedirs(sub, exist_ok=True)
            return LeR(
                npool=NPOOL,
                z_min=Z_MIN,
                z_max=Z_MAX,
                event_type="BBH",
                lens_type="sie_galaxy",
                interpolator_directory=interpolator_directory,
                ler_directory=sub,
                create_new_interpolator=False,
                lens_priors_params={
                    "velocity_dispersion": velocity_dispersion_prior(phistar),
                },
            )

        # Two populations: low vs high lens galaxy density normalization.
        ler_lo = make_ler_for_phistar("phistar_lo", PHISTAR_LOW)
        ler_hi = make_ler_for_phistar("phistar_hi", PHISTAR_HIGH)

        # Monte Carlo lensed sample for the low-``phistar`` universe.
        lensed_lo = ler_lo.lensed_cbc_statistics(
            size=SIZE,
            batch_size=BATCH_SIZE,
            resume=False,
            save_batch=False,
            output_jsonfile=False,
        )
        # Monte Carlo lensed sample for the high-``phistar`` universe.
        lensed_hi = ler_hi.lensed_cbc_statistics(
            size=SIZE,
            batch_size=BATCH_SIZE,
            resume=False,
            save_batch=False,
            output_jsonfile=False,
        )

        # Unlensed sample from the low-``phistar`` instance (reference CBC + network).
        unlensed = ler_lo.unlensed_cbc_statistics(
            size=SIZE,
            batch_size=BATCH_SIZE,
            resume=False,
            save_batch=False,
            output_jsonfile=False,
        )

        # Convert lensed MC counts to yr^-1 using ``normalization_pdf_z_lensed``.
        rate_lo, _ = ler_lo.lensed_rate(
            lensed_param=lensed_lo, output_jsonfile=False,
        )
        rate_hi, _ = ler_hi.lensed_rate(
            lensed_param=lensed_hi, output_jsonfile=False,
        )

        # Same unlensed sample for ratio numerator (only ``ler_lo`` cosmology/CBC path used).
        rate_unl, _ = ler_lo.unlensed_rate(
            unlensed_param=unlensed, output_jsonfile=False,
        )

        # How much the unlensed channel dominates vs lensed for each ``phistar``.
        ratio_lo = rate_unl / rate_lo
        ratio_hi = rate_unl / rate_hi

        print(
            f"\nphistar: low={PHISTAR_LOW} high={PHISTAR_HIGH} -> "
            f"lensed_rate: {rate_lo} vs {rate_hi} yr^-1\n"
            f"ratio: {ratio_lo} vs {ratio_hi}\n"
        )

        # Sanity: rates are physical and below the lensed pdf normalization.
        # Do not apply RATE_LOW_LENSED here: the deliberately low ``phistar``
        # population can produce a valid rate below that generic stochastic floor.
        rate_max_lens_lo = ler_lo.normalization_pdf_z_lensed
        rate_max_lens_hi = ler_hi.normalization_pdf_z_lensed

        assert (
            np.isfinite(rate_lo)
            and rate_lo > 0.0
            and rate_lo < rate_max_lens_lo
        ), f"phistar={PHISTAR_LOW}: rate_lo not finite or not in expected range"
        assert (
            np.isfinite(rate_hi)
            and rate_hi > 0.0
            and rate_hi < rate_max_lens_hi
        ), f"phistar={PHISTAR_HIGH}: rate_hi not finite or not in expected range"
        # More lens galaxies should increase the detection-weighted lensed rate here.
        assert rate_hi > rate_lo, (
            f"phistar={PHISTAR_HIGH}: expected higher rate for larger phistar, got {rate_hi} vs {rate_lo}"
        )

        # Higher ``phistar`` boosts lensed rate more than unlensed, so this ratio drops.
        assert ratio_hi < ratio_lo, (
            f"phistar={PHISTAR_HIGH}: expected higher ratio for smaller phistar, got {ratio_hi} vs {ratio_lo}"
        )
