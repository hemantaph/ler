"""
Unit tests for the LeR class.

Test Coverage:
--------------
- ``test_init``: core attributes, default ``json_file_names``, ``ler_directory``
  exists on disk, inheritance from ``LensGalaxyParameterDistribution``,
  custom ``pdet_finder`` stored, custom astropy cosmology stored.
- ``test_unlensed_cbc_statistics``: output keys and finite values; ``batch_size``
  and ``resume`` (first run 20 samples, second run resumes from batch 3 → 40 total,
  file content matches returned dict); ``output_jsonfile`` modes (``True``:
  default name, content verified; ``False``: no file written, verified).
- ``test_lensed_cbc_statistics``: output keys and finite values; ``batch_size``
  and ``resume`` (first run 10 samples, second run resumes → 20 total, file
  content matches returned dict); ``output_jsonfile`` modes; image columns and
  ``pdet_net`` shapes; optional ``include_redundant_parameters`` /
  ``include_effective_parameters`` key checks.
- ``test_lensed_cbc_statistics_redundant_and_effective``: small run with both
  flags true; asserts :data:`EXPECTED_REDUNDANT_KEYS` and
  :data:`EXPECTED_EFFECTIVE_KEYS`.
- ``test_unlensed_rate``: returns ``(rate, unlensed_param_detectable)``; rate
  equals ``normalization_pdf_z * detectable / total``; ``pdet_type=
  'probability_distribution'`` sums pdet_net directly; ``pdet_threshold``
  boundary cases: 0.0 (all detectable), 2.0 (none, rate=0); ``unlensed_param=None``
  loads from the default JSON file on disk; invalid ``pdet_type`` raises
  ``ValueError``.
- ``test_lensed_rate``: returns ``(rate, lensed_param_detectable)``; rate equals
  ``normalization_pdf_z_lensed * detectable / total``; ``pdet_threshold=2.0``
  (none detectable, rate=0); ``lensed_param=None`` loads from the default JSON
  file on disk.
- ``test_utilities``: ``rate_function`` math for both ``param_type='unlensed'``
  and ``param_type='lensed'``; ``_load_param`` returns a copy.
- ``test_selecting_n_unlensed_detectable_events``: ``stopping_criteria=None``
  accumulates batches until detectable events exceed ``size``; meta-data file has
  the expected keys and the final rate matches the return value; ``trim_to_size=False``
  gives final_size >= size; ``resume=True`` + ``trim_to_size=True`` randomly
  subsamples to ``new_size`` (every value present in the original collection);
  ``stopping_criteria=dict(...)`` stops when cumulative rate has converged and
  collected size exceeds ``size``.
- ``test_selecting_n_lensed_detectable_events``: same three test cases as the
  unlensed variant; uses smaller ``batch_size`` and ``size`` because lensed
  sampling is slower; ``pdet_threshold=[0.5, 0.5]`` (default 2-image condition).

"""

import os
import numpy as np
import pytest
from astropy.cosmology import LambdaCDM
from ler.rates import LeR
from ler.lens_galaxy_population import LensGalaxyParameterDistribution
from ler.utils import get_param_from_json
from tests_utils import (
    CommonTestUtils,
    EXPECTED_EFFECTIVE_KEYS,
    EXPECTED_IMAGE_KEYS,
    EXPECTED_SOURCE_POS_KEYS,
    EXPECTED_LENSED_PARAM_KEYS,
    EXPECTED_REDUNDANT_KEYS,
    EXPECTED_UNLENSED_GW_KEYS_ALIGNED_SPIN,
    EXPECTED_UNLENSED_PARAM_KEYS,
)


# ---------------------------------------------------------------------------
# Test configuration
# ---------------------------------------------------------------------------

# Default keys created by LeR.__init__ for parameter storage.
DEFAULT_JSON_KEYS = [
    "ler_params",
    "unlensed_param",
    "unlensed_param_detectable",
    "lensed_param",
    "lensed_param_detectable",
]

N_SAMPLES = 20

DEFAULT_CONFIG = dict(
    npool=1,
    z_min=0.0,
    z_max=10.0,
    event_type="BBH",
    lens_type="epl_shear_galaxy",
    create_new_interpolator=False,
    verbose=False,
)


def mock_pdet_finder(gw_param_dict):
    """Lightweight mock: detection probability = 1 for every event.

    Uses the first value in the dict to infer size, since per-image calls from
    lensed sampling do not include 'zs'.
    """
    size = len(list(gw_param_dict.values())[0])
    return dict(pdet_net=np.ones(size))


def _make_ler(interpolator_directory, ler_directory, **overrides):
    # Always provide a mock pdet_finder so gwsnr is not initialized.
    cfg = DEFAULT_CONFIG.copy()
    cfg.update(overrides)
    return LeR(
        pdet_finder=mock_pdet_finder,
        interpolator_directory=interpolator_directory,
        ler_directory=ler_directory,
        **cfg,
    )


@pytest.fixture(scope="module")
def ler_instance(interpolator_directory, ler_directory):
    """
    Module-scoped LeR instance reused across all tests in this file.

    Initialization is expensive (builds lens, optical depth, image properties,
    and CBC parameter samplers); a single instance is sufficient for unit-level
    checks.
    """
    return _make_ler(interpolator_directory, ler_directory)


@pytest.fixture(scope="module")
def unlensed_param(ler_instance):
    """
    Module-scoped unlensed parameter dict shared across unlensed_rate tests.

    Sampled once and reused to avoid redundant calls to unlensed_cbc_statistics.
    """
    return ler_instance.unlensed_cbc_statistics(
        size=N_SAMPLES, batch_size=N_SAMPLES, resume=False, output_jsonfile=False,
    )


@pytest.fixture(scope="module")
def lensed_param(ler_instance):
    """
    Module-scoped lensed parameter dict shared across lensed_rate tests.

    Sampled once and reused to avoid redundant calls to lensed_cbc_statistics.
    """
    return ler_instance.lensed_cbc_statistics(
        size=N_SAMPLES, batch_size=N_SAMPLES, resume=False, output_jsonfile=False,
    )


class TestLeR(CommonTestUtils):
    """Unit tests for the LeR class."""

    # -----------------------------------------------------------------------
    # Initialization
    # -----------------------------------------------------------------------

    def test_init(self, ler_instance, ler_directory, interpolator_directory):
        """
        Tests
        -----
        - Core attributes are set correctly: ``z_min``, ``z_max``,
          ``event_type``, ``lens_type``, ``cosmo``, ``npool``,
          ``ler_directory``.
        - Default ``json_file_names`` contains the expected keys.
        - ``ler_directory`` exists on disk.
        - LeR inherits from ``LensGalaxyParameterDistribution``.
        - A user-supplied ``pdet_finder`` is stored on ``self.pdet_finder``
          without triggering gwsnr initialization.
        """
        ler = ler_instance
        assert ler.z_min == DEFAULT_CONFIG["z_min"], \
            f"z_min: expected {DEFAULT_CONFIG['z_min']}, got {ler.z_min}"
        assert ler.z_max == DEFAULT_CONFIG["z_max"], \
            f"z_max: expected {DEFAULT_CONFIG['z_max']}, got {ler.z_max}"
        assert ler.event_type == DEFAULT_CONFIG["event_type"], \
            f"event_type: expected {DEFAULT_CONFIG['event_type']}, got {ler.event_type}"
        assert ler.lens_type == DEFAULT_CONFIG["lens_type"], \
            f"lens_type: expected {DEFAULT_CONFIG['lens_type']}, got {ler.lens_type}"
        assert ler.npool == DEFAULT_CONFIG["npool"], \
            f"npool: expected {DEFAULT_CONFIG['npool']}, got {ler.npool}"
        assert ler.cosmo is not None, "cosmo must not be None after initialization"
        assert ler.ler_directory == ler_directory, \
            f"ler_directory: expected {ler_directory}, got {ler.ler_directory}"
        assert os.path.isdir(ler.ler_directory), \
            f"ler_directory '{ler.ler_directory}' does not exist on disk"

        # json_file_names maps logical names → JSON filenames used for all disk I/O;
        # if any key is missing, the corresponding save/load call will silently fail.
        for key in DEFAULT_JSON_KEYS:
            assert key in ler.json_file_names, f"missing default json key: {key}"

        # LeR inherits the full lensed pipeline via LensGalaxyParameterDistribution
        assert isinstance(ler, LensGalaxyParameterDistribution), \
            "LeR must inherit from LensGalaxyParameterDistribution"

        # mock_pdet_finder simply returns 1 for every event, so no gwsnr initialization is triggered
        assert ler.pdet_finder is mock_pdet_finder, \
            "pdet_finder was not stored on self.pdet_finder after initialization"

    # -----------------------------------------------------------------------
    # unlensed_cbc_statistics
    # -----------------------------------------------------------------------

    def test_unlensed_cbc_statistics(self, ler_instance, ler_directory):
        """
        Tests
        -----
        - Output dict has all expected keys; arrays are finite with correct size.
        - ``pdet_net`` is in [0, 1] (mock returns 1 for every event).

        - ``batch_size`` and ``resume`` (uses a custom filename for this sub-test):
            - ``size=20, batch_size=10, resume=False`` → 2 batches, 20 samples saved to file.
            - ``size=40, batch_size=10, resume=True`` → resumes from batch 3,
              adds 2 more batches → 40 total samples; file content matches returned dict.

        - ``output_jsonfile`` modes:
            - ``True``  → file saved under ``self.json_file_names["unlensed_param"]``
              inside ``ler_directory``; file exists and content matches the returned dict.
            - ``False`` → no file written; return value is still a valid dict and
              the file does not exist on disk.
        """
        # --- basic output ---
        unlensed_param = ler_instance.unlensed_cbc_statistics(
            size=N_SAMPLES, batch_size=N_SAMPLES, resume=False, output_jsonfile=False,
        )
        self._assert_unlensed_cbc_outputs(
            unlensed_param,
            N_SAMPLES,
            expected_gw_keys=EXPECTED_UNLENSED_GW_KEYS_ALIGNED_SPIN,
        )

        # --- batch_size and resume ---
        fname_resume = "unit_unlensed_resume.json"

        # first run: 2 batches of 10 → 20 samples written to file
        unlensed_20 = ler_instance.unlensed_cbc_statistics(
            size=20, batch_size=10, resume=False, output_jsonfile=fname_resume,
        )
        assert len(unlensed_20["zs"]) == 20, \
            f"first run (resume=False, size=20): expected 20 samples, got {len(unlensed_20['zs'])}"

        # second run: resume=True, size=40 → batch_handler loads the 20 existing
        # samples and computes that batches 3 and 4 are still needed → 40 total
        unlensed_40 = ler_instance.unlensed_cbc_statistics(
            size=40, batch_size=10, resume=True, output_jsonfile=fname_resume,
        )
        assert len(unlensed_40["zs"]) == 40, \
            f"second run (resume=True, size=40): expected 40 samples, got {len(unlensed_40['zs'])}"

        path = os.path.join(ler_directory, fname_resume)
        assert os.path.isfile(path), \
            f"resume: file '{fname_resume}' not found in ler_directory after resume run"
        unlensed_from_file = get_param_from_json(path)
        self._assert_param_dicts_equal(unlensed_from_file, unlensed_40, label="resume")

        # --- output_jsonfile=True: file stored under the default name ---
        unlensed_param = ler_instance.unlensed_cbc_statistics(
            size=10, batch_size=10, resume=False, output_jsonfile=True,
        )
        default_fname = ler_instance.json_file_names["unlensed_param"]
        path = os.path.join(ler_directory, default_fname)
        assert os.path.isfile(path), \
            f"output_jsonfile=True: file '{default_fname}' not found in ler_directory"
        unlensed_from_file = get_param_from_json(path)
        self._assert_param_dicts_equal(
            unlensed_from_file, unlensed_param, label="output_jsonfile=True",
        )

        # --- output_jsonfile=False: no file written, but return value is valid ---
        os.remove(path)
        unlensed_no_file = ler_instance.unlensed_cbc_statistics(
            size=10, batch_size=10, resume=False, output_jsonfile=False,
        )
        self._assert_unlensed_cbc_outputs(
            unlensed_no_file,
            10,
            expected_gw_keys=EXPECTED_UNLENSED_GW_KEYS_ALIGNED_SPIN,
        )
        assert not os.path.isfile(path), \
            "output_jsonfile=False: file should not exist"

    # -----------------------------------------------------------------------
    # lensed_cbc_statistics
    # -----------------------------------------------------------------------

    def test_lensed_cbc_statistics(self, ler_instance, ler_directory):
        """
        Tests
        -----
        - Output dict has all expected 1-D lens keys with finite values and correct size.
        - Image keys: 2-D block with shape ``(N_SAMPLES, n_max_images)`` and
          ``nan_to_num`` in :meth:`_assert_param_dict_valid`; source positions 1-D;
        - ``pdet_net`` 2-D bounds checked via :meth:`_assert_lensed_pdet_net_2d`.

        - ``batch_size`` and ``resume`` (uses a custom filename for this sub-test):
            - ``size=10, batch_size=5, resume=False`` → 2 batches, 10 samples saved to file.
            - ``size=20, batch_size=5, resume=True`` → resumes, adds 2 more batches → 20
              total samples; file content matches returned dict.

        - ``output_jsonfile`` modes:
            - ``True``  → file saved under ``self.json_file_names["lensed_param"]``
              inside ``ler_directory``; file exists and content matches the returned dict.
            - ``False`` → no file written; return value is still a valid dict and
              the file does not exist on disk.
        """
        # --- basic output ---
        lensed_param = ler_instance.lensed_cbc_statistics(
            size=N_SAMPLES, batch_size=N_SAMPLES, resume=False, output_jsonfile=False,
        )
        nmax = ler_instance.n_max_images
        shape2 = (N_SAMPLES, nmax)
        self._assert_param_dict_valid(
            lensed_param, expected_keys=EXPECTED_LENSED_PARAM_KEYS, size=N_SAMPLES,
        )
        self._assert_param_dict_valid(
            lensed_param,
            expected_keys=EXPECTED_IMAGE_KEYS,
            size=shape2,
            nan_to_num=True,
        )
        self._assert_param_dict_valid(
            lensed_param,
            expected_keys=EXPECTED_SOURCE_POS_KEYS,
            size=N_SAMPLES,
        )
        self._assert_lensed_pdet_net_2d(lensed_param, shape2)

        # --- batch_size and resume ---
        fname_resume = "unit_lensed_resume.json"

        # first run: 2 batches of 5 → 10 samples written to file
        lensed_10 = ler_instance.lensed_cbc_statistics(
            size=10, batch_size=5, resume=False, output_jsonfile=fname_resume,
        )
        assert len(lensed_10["zs"]) == 10, \
            f"first run (resume=False, size=10): expected 10 samples, got {len(lensed_10['zs'])}"

        # second run: resume=True, size=20 → 20 total
        lensed_20 = ler_instance.lensed_cbc_statistics(
            size=20, batch_size=5, resume=True, output_jsonfile=fname_resume,
        )
        assert len(lensed_20["zs"]) == 20, \
            f"second run (resume=True, size=20): expected 20 samples, got {len(lensed_20['zs'])}"

        path = os.path.join(ler_directory, fname_resume)
        assert os.path.isfile(path), \
            f"resume: file '{fname_resume}' not found in ler_directory after resume run"
        lensed_from_file = get_param_from_json(path)
        self._assert_param_dicts_equal(lensed_from_file, lensed_20, label="resume")

        # --- output_jsonfile=True: file stored under the default name ---
        lensed_param = ler_instance.lensed_cbc_statistics(
            size=5, batch_size=5, resume=False, output_jsonfile=True,
        )
        default_fname = ler_instance.json_file_names["lensed_param"]
        path = os.path.join(ler_directory, default_fname)
        assert os.path.isfile(path), \
            f"output_jsonfile=True: file '{default_fname}' not found in ler_directory"
        lensed_from_file = get_param_from_json(path)
        self._assert_param_dicts_equal(
            lensed_from_file, lensed_param, label="output_jsonfile=True",
        )

        # --- output_jsonfile=False: no file written, but return value is valid ---
        os.remove(path)
        lensed_no_file = ler_instance.lensed_cbc_statistics(
            size=5, batch_size=5, resume=False, output_jsonfile=False,
        )
        shape5 = (5, nmax)
        self._assert_param_dict_valid(
            lensed_no_file, expected_keys=EXPECTED_LENSED_PARAM_KEYS, size=5,
        )
        self._assert_param_dict_valid(
            lensed_no_file,
            expected_keys=EXPECTED_IMAGE_KEYS,
            size=shape5,
            nan_to_num=True,
        )
        self._assert_param_dict_valid(
            lensed_no_file,
            expected_keys=EXPECTED_SOURCE_POS_KEYS,
            size=5,
        )
        self._assert_lensed_pdet_net_2d(lensed_no_file, shape5)
        assert not os.path.isfile(path), \
            "output_jsonfile=False: file should not exist"

    def test_lensed_cbc_statistics_redundant_and_effective(
        self,
        interpolator_directory,
        ler_directory,
    ):
        """
        Tests
        -----
        - With ``include_redundant_parameters=True`` and
          ``include_effective_parameters=True``, ``lensed_cbc_statistics`` includes
          :data:`EXPECTED_REDUNDANT_KEYS` (1-D) and :data:`EXPECTED_EFFECTIVE_KEYS`
          (2-D with NaN padding), in addition to the default lensed / image / ``pdet_net``
          columns validated elsewhere.
        """
        sub = os.path.join(ler_directory, "unit_lensed_redundant_effective")
        os.makedirs(sub, exist_ok=True)
        ler = _make_ler(
            interpolator_directory,
            sub,
            include_redundant_parameters=True,
            include_effective_parameters=True,
        )
        n_small = 5
        lp = ler.lensed_cbc_statistics(
            size=n_small,
            batch_size=n_small,
            resume=False,
            output_jsonfile=False,
        )
        nmax = ler.n_max_images
        shape_small = (n_small, nmax)
        self._assert_param_dict_valid(
            lp, expected_keys=EXPECTED_LENSED_PARAM_KEYS, size=n_small,
        )
        self._assert_param_dict_valid(
            lp,
            expected_keys=EXPECTED_IMAGE_KEYS,
            size=shape_small,
            nan_to_num=True,
        )
        self._assert_param_dict_valid(
            lp,
            expected_keys=EXPECTED_SOURCE_POS_KEYS,
            size=n_small,
        )
        self._assert_lensed_pdet_net_2d(lp, shape_small)
        self._assert_param_dict_valid(
            lp, expected_keys=EXPECTED_REDUNDANT_KEYS, size=n_small,
        )
        self._assert_param_dict_valid(
            lp,
            expected_keys=EXPECTED_EFFECTIVE_KEYS,
            size=shape_small,
            nan_to_num=True,
        )

    # -----------------------------------------------------------------------
    # unlensed_rate
    # -----------------------------------------------------------------------

    def test_unlensed_rate(self, ler_instance, unlensed_param):
        """
        Tests
        -----
        - Returns ``(rate, unlensed_param_detectable)``; rate is a finite positive float.
        - ``pdet_type='boolean'``: rate = ``normalization_pdf_z * detectable / total``
          where detectable counts events with pdet_net >= pdet_threshold (default 0.5).
          With mock pdet=1, all N_SAMPLES events are detectable.
        - ``pdet_type='probability_distribution'``: rate uses sum(pdet_net) directly.
        - ``pdet_threshold=0.0``: all events pass; ``pdet_threshold=2.0``: none pass,
          rate == 0.
        - ``unlensed_param=None``: loads parameters from the default JSON file on disk.
        - Invalid ``pdet_type`` raises ``ValueError``.
        """
        # --- boolean mode ---
        rate, detectable_param = ler_instance.unlensed_rate(
            unlensed_param=unlensed_param, pdet_type="boolean", output_jsonfile=False,
        )
        assert isinstance(rate, float) and np.isfinite(rate) and rate > 0, \
            f"unlensed_rate (boolean): expected finite positive float, got {rate}"
        # boolean mode: count events where pdet_net >= pdet_threshold (default 0.5)
        detectable = float(np.sum(np.asarray(unlensed_param["pdet_net"]) >= 0.5))
        # Monte Carlo rate formula: R = normalization_pdf_z * (detectable / total)
        np.testing.assert_allclose(
            rate,
            ler_instance.normalization_pdf_z * detectable / N_SAMPLES,
            rtol=1e-10,
        )
        # with mock pdet=1 and threshold=0.5, all N_SAMPLES events are detectable
        self._assert_unlensed_cbc_outputs(
            detectable_param,
            N_SAMPLES,
            expected_gw_keys=EXPECTED_UNLENSED_GW_KEYS_ALIGNED_SPIN,
        )

        # --- probability_distribution mode ---
        rate, _ = ler_instance.unlensed_rate(
            unlensed_param=unlensed_param, pdet_type="probability_distribution",
            output_jsonfile=False,
        )
        assert isinstance(rate, float) and np.isfinite(rate) and rate > 0, \
            f"unlensed_rate (probability_distribution): expected finite positive float, got {rate}"
        # probability_distribution mode: sum pdet_net directly instead of thresholding
        detectable = float(np.sum(np.asarray(unlensed_param["pdet_net"])))
        np.testing.assert_allclose(
            rate,
            ler_instance.normalization_pdf_z * detectable / N_SAMPLES,
            rtol=1e-10,
        )

        # --- pdet_threshold: boundary cases ---
        # threshold=0.0 → all events pass (pdet_net=1 >= 0.0)
        rate_all, params_all = ler_instance.unlensed_rate(
            unlensed_param=unlensed_param, pdet_threshold=0.0, pdet_type="boolean",
            output_jsonfile=False,
        )
        assert isinstance(rate_all, float) and np.isfinite(rate_all) and rate_all > 0, \
            f"pdet_threshold=0.0: expected finite positive float rate, got {rate_all}"
        assert len(params_all["zs"]) == N_SAMPLES, \
            f"pdet_threshold=0.0: expected all {N_SAMPLES} events detectable, got {len(params_all['zs'])}"

        # threshold=2.0 → no event has pdet_net >= 2.0 → zero detectable, rate = 0
        rate_none, params_none = ler_instance.unlensed_rate(
            unlensed_param=unlensed_param, pdet_threshold=2.0, pdet_type="boolean",
            output_jsonfile=False,
        )
        assert len(params_none["zs"]) == 0, \
            f"pdet_threshold=2.0: expected 0 detectable events, got {len(params_none['zs'])}"
        assert rate_none == 0.0, \
            f"pdet_threshold=2.0: expected rate=0, got {rate_none}"

        # --- unlensed_param=None: load from file ---
        # write the default file first so unlensed_rate can find it
        ler_instance.unlensed_cbc_statistics(
            size=N_SAMPLES, batch_size=N_SAMPLES, resume=False, output_jsonfile=True,
        )
        fname = "unit_unlensed_detectable.json"
        rate, loaded_param = ler_instance.unlensed_rate(
            unlensed_param=None, output_jsonfile=fname,
        )
        assert isinstance(rate, float) and np.isfinite(rate) and rate > 0, \
            f"unlensed_rate (unlensed_param=None): expected finite positive float, got {rate}"
        self._assert_unlensed_cbc_outputs(
            loaded_param,
            N_SAMPLES,
            expected_gw_keys=EXPECTED_UNLENSED_GW_KEYS_ALIGNED_SPIN,
        )
        loaded_param_from_file = get_param_from_json(
            os.path.join(ler_instance.ler_directory, fname)
        )
        self._assert_param_dicts_equal(
            loaded_param_from_file, loaded_param, label="unlensed_rate (param=None)",
        )

        # --- invalid pdet_type raises ValueError ---
        with pytest.raises(ValueError):
            ler_instance.unlensed_rate(
                unlensed_param=unlensed_param, pdet_type="not_valid",
                output_jsonfile=False,
            )

    # -----------------------------------------------------------------------
    # lensed_rate
    # -----------------------------------------------------------------------

    def test_lensed_rate(self, ler_instance, lensed_param):
        """
        Tests
        -----
        - Returns ``(rate, lensed_param_detectable)``; rate is a finite positive float.
        - ``pdet_type='boolean'`` (default): rate equals
          ``normalization_pdf_z_lensed * len(detectable) / total``.
        - ``pdet_threshold=2.0`` (per image): no event has pdet >= 2.0, so rate = 0.
        - ``lensed_param=None``: loads parameters from the default JSON file on disk.
        """
        # --- boolean mode ---
        rate, detectable_param = ler_instance.lensed_rate(
            lensed_param=lensed_param, pdet_type="boolean", output_jsonfile=False,
        )
        assert isinstance(rate, float) and np.isfinite(rate) and rate > 0, \
            f"lensed_rate (boolean): expected finite positive float, got {rate}"
        # boolean mode: rate = normalization_pdf_z_lensed * n_detectable / total
        # with multi-image threshold [0.5, 0.5], detectable = len(detectable_param)
        detectable_size = len(detectable_param["zs"])
        np.testing.assert_allclose(
            rate,
            ler_instance.normalization_pdf_z_lensed * detectable_size / N_SAMPLES,
            rtol=1e-10,
        )
        self._assert_param_dict_valid(
            detectable_param, expected_keys=EXPECTED_LENSED_PARAM_KEYS, size=detectable_size,
        )

        # --- pdet_threshold=2.0 → no event has pdet >= 2.0 → rate = 0 ---
        rate_none, params_none = ler_instance.lensed_rate(
            lensed_param=lensed_param, pdet_threshold=2.0, num_img=1,
            pdet_type="boolean", output_jsonfile=False,
        )
        assert len(params_none["zs"]) == 0, \
            f"pdet_threshold=2.0: expected 0 detectable events, got {len(params_none['zs'])}"
        assert rate_none == 0.0, \
            f"pdet_threshold=2.0: expected rate=0, got {rate_none}"

        # --- lensed_param=None: load from file ---
        ler_instance.lensed_cbc_statistics(
            size=N_SAMPLES, batch_size=N_SAMPLES, resume=False, output_jsonfile=True,
        )
        fname = "unit_lensed_detectable.json"
        rate, loaded_param = ler_instance.lensed_rate(
            lensed_param=None, output_jsonfile=fname,
        )
        assert isinstance(rate, float) and np.isfinite(rate) and rate > 0, \
            f"lensed_rate (lensed_param=None): expected finite positive float, got {rate}"
        loaded_size = len(loaded_param["zs"])
        self._assert_param_dict_valid(
            loaded_param, expected_keys=EXPECTED_LENSED_PARAM_KEYS, size=loaded_size,
        )
        loaded_param_from_file = get_param_from_json(
            os.path.join(ler_instance.ler_directory, fname)
        )
        self._assert_param_dicts_equal(
            loaded_param_from_file, loaded_param, label="lensed_rate (param=None)",
        )

    # -----------------------------------------------------------------------
    # rate_function and _load_param
    # -----------------------------------------------------------------------

    def test_utilities(self, ler_instance):
        """
        Tests
        -----
        - ``rate_function`` returns the correct value for both ``param_type='unlensed'``
          (``normalization_pdf_z * detectable / total``) and ``param_type='lensed'``
          (``normalization_pdf_z_lensed * detectable / total``); outputs are finite
          positive floats.
        - ``_load_param`` returns a copy (not the same dict object) when given a dict,
          so downstream mutations do not affect the caller's dict; tested for both
          ``param_type='unlensed'`` and ``param_type='lensed'``.
        """
        # --- rate_function math ---
        detectable, total = 100, 1000

        rate_u = ler_instance.rate_function(detectable, total, param_type="unlensed", verbose=False)
        rate_l = ler_instance.rate_function(detectable, total, param_type="lensed", verbose=False)

        assert isinstance(rate_u, float) and np.isfinite(rate_u) and rate_u > 0, \
            f"rate_function (unlensed): expected finite positive float, got {rate_u}"
        assert isinstance(rate_l, float) and np.isfinite(rate_l) and rate_l > 0, \
            f"rate_function (lensed): expected finite positive float, got {rate_l}"
        # Monte Carlo rate formula: R = normalization * (detectable / total)
        np.testing.assert_allclose(
            rate_u, ler_instance.normalization_pdf_z * detectable / total, rtol=1e-12,
            err_msg="rate_function (unlensed): rate does not match normalization_pdf_z * det / total",
        )
        np.testing.assert_allclose(
            rate_l, ler_instance.normalization_pdf_z_lensed * detectable / total, rtol=1e-12,
            err_msg="rate_function (lensed): rate does not match normalization_pdf_z_lensed * det / total",
        )

        # --- _load_param returns a copy ---
        original = dict(
            zs=np.array([1.0, 2.0]),
            pdet_net=np.array([0.6, 0.9]),
        )
        for param_type in ("unlensed", "lensed"):
            loaded = ler_instance._load_param(original, param_type=param_type)
            # must be a different object so downstream mutations don't affect the caller
            assert loaded is not original, \
                f"_load_param (param_type='{param_type}'): returned the same dict object, not a copy"
            assert set(loaded.keys()) == set(original.keys()), \
                f"_load_param (param_type='{param_type}'): key mismatch"
            np.testing.assert_array_equal(loaded["zs"], original["zs"])
            np.testing.assert_array_equal(loaded["pdet_net"], original["pdet_net"])

    # -----------------------------------------------------------------------
    # selecting_n_unlensed_detectable_events
    # -----------------------------------------------------------------------

    def test_selecting_n_unlensed_detectable_events(self, ler_instance):
        """
        Tests
        -----
        - ``stopping_criteria=None``: accumulates batches until detectable events
          exceed ``size``; meta-data file has keys ``events_total``,
          ``detectable_events``, ``total_rate``; final rate matches the return value.
        - ``trim_to_size=False``: returned dict has final_size >= size; stored JSON
          file matches the returned dict.
        - ``resume=True`` + ``trim_to_size=True``: no new sampling (existing file
          already satisfies size); result is randomly trimmed to ``new_size`` and
          every value is a member of the original collection.
        - ``stopping_criteria=dict(...)``: stops when the cumulative rate has
          converged (relative difference of the last 4 batches < 0.5%) and
          collected size exceeds ``size``.
        """
        size = 500
        batch_size = 10000
        output_jsonfile = "unit_n_unlensed_detectable.json"
        meta_data_file = "unit_meta_unlensed.json"

        # --- stopping_criteria=None, trim_to_size=False ---
        rate, param = ler_instance.selecting_n_unlensed_detectable_events(
            size=size,
            batch_size=batch_size,
            stopping_criteria=None,
            pdet_threshold=0.5,
            resume=False,
            trim_to_size=False,
            output_jsonfile=output_jsonfile,
            meta_data_file=meta_data_file,
        )

        meta_data = get_param_from_json(
            os.path.join(ler_instance.ler_directory, meta_data_file)
        )
        expected_meta_keys = ("events_total", "detectable_events", "total_rate", "batch_rate")
        assert isinstance(meta_data, dict) and len(meta_data) == len(expected_meta_keys), \
            f"meta_data: expected dict with {len(expected_meta_keys)} keys, got {len(meta_data)}"
        for key in expected_meta_keys:
            assert key in meta_data, f"meta_data: missing key '{key}'"
        assert isinstance(meta_data["total_rate"][-1], float) \
            and np.isfinite(meta_data["total_rate"][-1]) \
            and meta_data["total_rate"][-1] > 0, \
            f"meta_data: expected finite positive float total_rate, got {meta_data['total_rate'][-1]}"
        assert isinstance(meta_data["batch_rate"][-1], float) \
            and np.isfinite(meta_data["batch_rate"][-1]) \
            and meta_data["batch_rate"][-1] > 0, \
            f"meta_data: expected finite positive float batch_rate, got {meta_data['batch_rate'][-1]}"
        assert meta_data["total_rate"][-1] == rate, \
            f"meta_data: final total_rate {meta_data['total_rate'][-1]} != returned rate {rate}"

        # trim_to_size=False → all collected events returned (final_size >= size)
        final_size = meta_data["detectable_events"][-1]
        self._assert_unlensed_cbc_outputs(
            param,
            final_size,
            expected_gw_keys=EXPECTED_UNLENSED_GW_KEYS_ALIGNED_SPIN,
        )
        dict_from_file = get_param_from_json(
            os.path.join(ler_instance.ler_directory, output_jsonfile)
        )
        self._assert_param_dicts_equal(
            dict_from_file, param, label="selecting_n_unlensed (no_trim)",
        )

        # --- resume=True, trim_to_size=True ---
        # existing file already has final_size >= new_size, so no new sampling occurs;
        # _trim_results_to_size uses np.random.choice → new_param is a random subset
        new_size = 600
        _, new_param = ler_instance.selecting_n_unlensed_detectable_events(
            size=new_size,
            batch_size=batch_size,
            stopping_criteria=None,
            pdet_threshold=0.5,
            resume=True,
            trim_to_size=True,
            output_jsonfile=output_jsonfile,
            meta_data_file=meta_data_file,
        )
        self._assert_unlensed_cbc_outputs(
            new_param,
            new_size,
            expected_gw_keys=EXPECTED_UNLENSED_GW_KEYS_ALIGNED_SPIN,
        )
        for key in EXPECTED_UNLENSED_PARAM_KEYS:
            assert np.isin(new_param[key], param[key]).all(), \
                f"selecting_n_unlensed (resume+trim): '{key}' contains values not in original param"

        # --- stopping_criteria: stops when rate has converged ---
        stopping_criteria = dict(relative_diff_percentage=0.5, number_of_last_batches_to_check=4)
        sc_output = "unit_n_unlensed_sc.json"
        sc_meta = "unit_meta_unlensed_sc.json"
        _, _ = ler_instance.selecting_n_unlensed_detectable_events(
            size=new_size,
            batch_size=batch_size,
            stopping_criteria=stopping_criteria,
            pdet_threshold=0.5,
            resume=False,
            trim_to_size=False,
            output_jsonfile=sc_output,
            meta_data_file=sc_meta,
        )
        sc_meta_data = get_param_from_json(
            os.path.join(ler_instance.ler_directory, sc_meta)
        )
        # last 4 cumulative rates must differ by < 0.5%
        idx_last4 = [-4, -3, -2, -1]
        rates_last4 = np.array(sc_meta_data["total_rate"])[idx_last4]
        pct_diff = np.abs((rates_last4 - rates_last4[-1]) / rates_last4[-1]) * 100
        assert np.all(pct_diff < 0.5), \
            "selecting_n_unlensed (stopping_criteria): last 4 cumulative rates must differ by < 0.5%"
        assert sc_meta_data["detectable_events"][-1] > new_size, \
            f"selecting_n_unlensed (stopping_criteria): expected collected_size > {new_size}, " \
            f"got {sc_meta_data['detectable_events'][-1]}"

    # -----------------------------------------------------------------------
    # selecting_n_lensed_detectable_events
    # -----------------------------------------------------------------------

    def test_selecting_n_lensed_detectable_events(self, ler_instance):
        """
        Tests
        -----
        - ``stopping_criteria=None``: accumulates batches until detectable events
          exceed ``size``; meta-data file has keys ``events_total``,
          ``detectable_events``, ``total_rate``, ``batch_rate``; final rate matches the return value.
        - ``trim_to_size=False``: returned dict has final_size >= size; stored JSON
          file matches the returned dict.
        - ``resume=True`` + ``trim_to_size=True``: no new sampling (existing file
          already satisfies size); result is randomly trimmed to ``new_size`` and
          every value is a member of the original collection.
        - ``stopping_criteria=dict(...)``: stops when the cumulative rate has
          converged (relative difference of the last 4 batches < 2%) and
          collected size exceeds ``size``.

        Note: ``batch_size`` and ``size`` are kept small because each lensed batch
        involves lens equation solving, which is significantly slower than unlensed.
        """
        size = 30
        batch_size = 500
        output_jsonfile = "unit_n_lensed_detectable.json"
        meta_data_file = "unit_meta_lensed.json"

        # --- stopping_criteria=None, trim_to_size=False ---
        rate, param = ler_instance.selecting_n_lensed_detectable_events(
            size=size,
            batch_size=batch_size,
            stopping_criteria=None,
            pdet_threshold=[0.5, 0.5],
            num_img=[1, 1],
            resume=False,
            trim_to_size=False,
            output_jsonfile=output_jsonfile,
            meta_data_file=meta_data_file,
        )

        meta_data = get_param_from_json(
            os.path.join(ler_instance.ler_directory, meta_data_file)
        )
        expected_meta_keys = ("events_total", "detectable_events", "total_rate", "batch_rate")
        assert isinstance(meta_data, dict) and len(meta_data) == len(expected_meta_keys), \
            f"meta_data: expected dict with {len(expected_meta_keys)} keys, got {len(meta_data)}"
        for key in expected_meta_keys:
            assert key in meta_data, f"meta_data: missing key '{key}'"
        assert isinstance(meta_data["total_rate"][-1], float) \
            and np.isfinite(meta_data["total_rate"][-1]) \
            and meta_data["total_rate"][-1] > 0, \
            f"meta_data: expected finite positive float total_rate, got {meta_data['total_rate'][-1]}"
        assert isinstance(meta_data["batch_rate"][-1], float) \
            and np.isfinite(meta_data["batch_rate"][-1]) \
            and meta_data["batch_rate"][-1] > 0, \
            f"meta_data: expected finite positive float batch_rate, got {meta_data['batch_rate'][-1]}"
        assert meta_data["total_rate"][-1] == rate, \
            f"meta_data: final total_rate {meta_data['total_rate'][-1]} != returned rate {rate}"

        # trim_to_size=False → all collected events returned (final_size >= size)
        final_size = meta_data["detectable_events"][-1]
        self._assert_param_dict_valid(
            param, expected_keys=EXPECTED_LENSED_PARAM_KEYS, size=final_size,
        )
        dict_from_file = get_param_from_json(
            os.path.join(ler_instance.ler_directory, output_jsonfile)
        )
        self._assert_param_dicts_equal(
            dict_from_file, param, label="selecting_n_lensed (no_trim)",
        )

        # --- resume=True, trim_to_size=True ---
        # existing file already has final_size >= new_size, so no new sampling occurs;
        # _trim_results_to_size uses np.random.choice → new_param is a random subset
        new_size = 60
        _, new_param = ler_instance.selecting_n_lensed_detectable_events(
            size=new_size,
            batch_size=batch_size,
            stopping_criteria=None,
            pdet_threshold=[0.5, 0.5],
            num_img=[1, 1],
            resume=True,
            trim_to_size=True,
            output_jsonfile=output_jsonfile,
            meta_data_file=meta_data_file,
        )
        self._assert_param_dict_valid(
            new_param, expected_keys=EXPECTED_LENSED_PARAM_KEYS, size=new_size,
        )
        for key in EXPECTED_LENSED_PARAM_KEYS:
            assert np.isin(new_param[key], param[key]).all(), \
                f"selecting_n_lensed (resume+trim): '{key}' contains values not in original param"

        # --- stopping_criteria: stops when rate has converged ---
        # default lensed stopping_criteria uses 2% threshold (rate converges more slowly
        # than unlensed because lensed sampling is more stochastic)
        stopping_criteria = dict(relative_diff_percentage=2.0, number_of_last_batches_to_check=4)
        sc_output = "unit_n_lensed_sc.json"
        sc_meta = "unit_meta_lensed_sc.json"
        _, _ = ler_instance.selecting_n_lensed_detectable_events(
            size=new_size,
            batch_size=batch_size,
            stopping_criteria=stopping_criteria,
            pdet_threshold=[0.5, 0.5],
            num_img=[1, 1],
            resume=False,
            trim_to_size=False,
            output_jsonfile=sc_output,
            meta_data_file=sc_meta,
        )
        sc_meta_data = get_param_from_json(
            os.path.join(ler_instance.ler_directory, sc_meta)
        )
        # last 4 cumulative rates must differ by < 2%
        idx_last4 = [-4, -3, -2, -1]
        rates_last4 = np.array(sc_meta_data["total_rate"])[idx_last4]
        pct_diff = np.abs((rates_last4 - rates_last4[-1]) / rates_last4[-1]) * 100
        assert np.all(pct_diff < 2.0), \
            "selecting_n_lensed (stopping_criteria): last 4 cumulative rates must differ by < 2%"
        assert sc_meta_data["detectable_events"][-1] > new_size, \
            f"selecting_n_lensed (stopping_criteria): expected collected_size > {new_size}, " \
            f"got {sc_meta_data['detectable_events'][-1]}"
