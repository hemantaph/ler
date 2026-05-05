"""
Integration test: selection function ``xi(Lambda)`` with ``GWRATES``.

This follows ``docs/examples/selection_function.ipynb``: for hyperparameters
``Lambda`` of the primary-mass ``broken_powerlaw_plus_2peaks`` model, the
selection function is estimated as the Monte Carlo average

    xi(Lambda) = <p_det(theta) >_{theta ~ p(theta | Lambda)} ,

using the network ``pdet_net`` column from ``gw_cbc_statistics``.

Test Coverage:
--------------
- ``test_gwrates_selection_function_mean_pdet``: two fixed hyperparameter
  points; swap ``GWRATES.mass_1_source`` to a callable that draws from
  ``broken_powerlaw_plus_2peaks_rvs``; run ``gw_cbc_statistics``; assert
  ``mean(pdet_net)`` is finite and in ``(0, 1]`` for each point, and the two
  values differ measurably.

Runs in default CI (``.github/workflows/tests.yml``).

"""

import os

import numpy as np
import pytest

from ler.rates import GWRATES
from ler.gw_source_population import broken_powerlaw_plus_2peaks_rvs

from tests_utils import CommonTestUtils, EXPECTED_ALIGNED_MODE_SPIN_KEYS, clamp_npool_for_numba

# ---------------------------------------------------------------------------
# Test configuration
# ---------------------------------------------------------------------------

# Smaller than the full demonstration scale (100k) so the integration run stays bounded;
# large enough for a stable mean ``pdet_net``.
SIZE = 8000
BATCH_SIZE = 4000

NPOOL = clamp_npool_for_numba(6)

Z_MIN = 0.0
Z_MAX = 10.0

def make_mass_1_source_sampler(hyp):
    """Return ``mass_1_source(size)`` drawn from ``broken_powerlaw_plus_2peaks_rvs``."""

    def sampler(size):
        return broken_powerlaw_plus_2peaks_rvs(
            size=size,
            lam_0=hyp["lam_0"],
            lam_1=hyp["lam_1"],
            mpp_1=hyp["mpp_1"],
            sigpp_1=hyp["sigpp_1"],
            mpp_2=hyp["mpp_2"],
            sigpp_2=hyp["sigpp_2"],
            mlow_1=hyp["mlow_1"],
            delta_m_1=hyp["delta_m_1"],
            break_mass=hyp["break_mass"],
            alpha_1=hyp["alpha_1"],
            alpha_2=hyp["alpha_2"],
            mmax=300.0,
        )

    return sampler

class TestGWRATESSelectionFunction(CommonTestUtils):
    """``GWRATES`` + ``gwsnr``: selection function as mean ``pdet_net``."""

    def test_gwrates_selection_function_mean_pdet(
        self,
        interpolator_directory,
        ler_directory,
    ):
        """
        Tests
        -----
        - Build ``GWRATES`` (BBH, precessing, same network setup as other
          integration tests).
        - For each of two hyperparameter dicts ``Lambda``: set
          ``ler.mass_1_source`` to a sampler from ``broken_powerlaw_plus_2peaks_rvs``;
          run ``gw_cbc_statistics``; take ``xi = mean(pdet_net)``.
        - ``xi`` must be finite and in ``(0, 1]`` (detection probabilities).
        - The two hyperparameter rows must give measurably different ``xi`` so
          the test checks that ``Lambda`` actually moves the mean ``pdet_net``.
        """
        # Reference hyperparameters: defaults documented on ``broken_powerlaw_plus_2peaks_rvs``.
        hyper_ref = {
            "lam_0": 0.361,
            "lam_1": 0.586,
            "mpp_1": 9.764,
            "sigpp_1": 0.649,
            "mpp_2": 32.763,
            "sigpp_2": 3.918,
            "mlow_1": 5.059,
            "delta_m_1": 4.321,
            "break_mass": 35.622,
            "alpha_1": 1.728,
            "alpha_2": 4.512,
        }
        # Second point: emphasize high-mass power-law tail (different typical m1 / SNR / pdet).
        hyper_stiff = {
            "lam_0": 0.05,
            "lam_1": 0.05,
            "mpp_1": 12.0,
            "sigpp_1": 2.0,
            "mpp_2": 45.0,
            "sigpp_2": 5.0,
            "mlow_1": 8.0,
            "delta_m_1": 5.0,
            "break_mass": 55.0,
            "alpha_1": 0.5,
            "alpha_2": 9.0,
        }

        ler_dir = os.path.join(ler_directory, "gwrates_selection_function")
        os.makedirs(ler_dir, exist_ok=True)

        # Single instance; only ``mass_1_source`` changes between hyperparameter points.
        ler = GWRATES(
            npool=NPOOL,
            z_min=Z_MIN,
            z_max=Z_MAX,
            event_type="BBH",
            interpolator_directory=interpolator_directory,
            ler_directory=ler_dir,
            verbose=False,
        )

        selection_xi = []

        for hyp in (hyper_ref, hyper_stiff):
            # Inject primary-mass prior for this ``Lambda``.
            ler.mass_1_source = make_mass_1_source_sampler(hyp)

            # Joint CBC draw + ``gwsnr`` ``pdet_net`` for each event.
            gw_param = ler.gw_cbc_statistics(
                size=SIZE,
                batch_size=BATCH_SIZE,
                resume=False,
                save_batch=False,
                output_jsonfile=False,
            )

            # Same columns as other BBH precessing integration runs.
            self._assert_unlensed_cbc_outputs(
                gw_param,
                SIZE,
                expected_gw_keys=EXPECTED_ALIGNED_MODE_SPIN_KEYS,
            )

            # xi(Lambda) = average detection probability under p(theta | Lambda).
            pdet_net = gw_param["pdet_net"]
            xi = float(np.mean(pdet_net))
            selection_xi.append(xi)

        xi_ref, xi_stiff = selection_xi

        # Strict interior (0, 1]: Monte Carlo mean of Bernoulli-ish weights.
        assert np.isfinite(xi_ref) and 0.0 < xi_ref <= 1.0, f"xi(ref): {xi_ref}"
        assert np.isfinite(xi_stiff) and 0.0 < xi_stiff <= 1.0, f"xi(stiff): {xi_stiff}"

        # Mass model changed; we only require the two averages are not bitwise identical.
        # Require a clear separation in mean pdet (two distinct mass models).
        diff_xi = abs(xi_ref - xi_stiff)
        assert diff_xi > 5e-5, (
            f"expected measurably different xi for different Lambda, got {xi_ref} vs {xi_stiff} (|Δ|={diff_xi})"
        )
