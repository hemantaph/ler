"""
Unit tests for the GWRATES class.

Test Coverage:
--------------
- ``test_init``: core attributes, default ``json_file_names``, ``ler_directory``
  exists on disk, inheritance from ``CBCSourceParameterDistribution``,
  custom ``pdet_finder`` stored, custom astropy cosmology stored.
- ``test_gw_cbc_statistics``: output keys and finite values; ``batch_size`` and
  ``resume`` (first run 20 samples, second run resumes from batch 3 → 40 total,
  file content matches returned dict); ``output_jsonfile`` modes (``True``:
  default name, content verified; ``False``: no file written, verified).
- ``test_gw_rate``: returns ``(rate, gw_param_detectable)``; rate equals
  ``normalization_pdf_z * detectable / total``; ``pdet_type='probability_distribution'``
  sums pdet_net directly; ``pdet_threshold`` boundary cases: 0.0 (all detectable),
  2.0 (none, rate=0); ``gw_param=None`` loads from the default JSON file on disk;
  invalid ``pdet_type`` raises ``ValueError``.
- ``test_selecting_n_detectable_events``: collects at least ``size`` detectable
  events when ``stopping_criteria=None``; ``trim_to_size=True`` trims randomly to
  exactly ``size`` and values are a subset of the original; ``resume=True`` from
  an already-complete file skips resampling; ``stopping_criteria`` stops when the
  rate converges and the collected size exceeds ``size``.
- ``test_utilities``: ``rate_function`` math; ``_load_param`` returns a copy.

"""

import os
import numpy as np
import pytest
from astropy.cosmology import LambdaCDM
from tests_utils import (
    CommonTestUtils,
    EXPECTED_UNLENSED_GW_KEYS_ALIGNED_SPIN,
    EXPECTED_UNLENSED_PARAM_KEYS,
    clamp_npool_for_numba,
)

from ler.rates import GWRATES
from ler.gw_source_population import CBCSourceParameterDistribution
from ler.utils import get_param_from_json


# ---------------------------------------------------------------------------
# Test configuration
# ---------------------------------------------------------------------------

# Default keys created by GWRATES.__init__ for parameter storage.
DEFAULT_JSON_KEYS = [
    "gwrates_params",
    "gw_param",
    "gw_param_detectable",
]

# Minimum keys expected in gw_cbc_statistics output (see ``tests_utils``).
EXPECTED_PARAM_KEYS = EXPECTED_UNLENSED_PARAM_KEYS

N_SAMPLES = 20

DEFAULT_CONFIG = dict(
    # Preferred high npool on big machines (6); clamp so Numba succeeds on slim CI VMs.
    npool=clamp_npool_for_numba(6),
    z_min=0.0,
    z_max=10.0,
    event_type="BBH",
    create_new_interpolator=False,
    verbose=False,
)

def mock_pdet_finder(gw_param_dict):
    """Lightweight mock: detection probability = 1 for every event."""
    size = len(gw_param_dict["zs"])
    return dict(pdet_net=np.ones(size))


def _make_gwrates(interpolator_directory, ler_directory, **overrides):
    # Always provide a mock pdet_finder so gwsnr is not initialized.
    cfg = DEFAULT_CONFIG.copy()
    cfg.update(overrides)
    return GWRATES(
        pdet_finder=mock_pdet_finder,
        interpolator_directory=interpolator_directory,
        ler_directory=ler_directory,
        **cfg,
    )


@pytest.fixture(scope="module")
def gwrates_instance(interpolator_directory, ler_directory):
    """
    Module-scoped GWRATES instance reused across all tests in this file.

    Uses a mock pdet_finder so gwsnr is never initialized.
    """
    return _make_gwrates(interpolator_directory, ler_directory)


@pytest.fixture(scope="module")
def gw_param(gwrates_instance):
    """
    Module-scoped sample dict shared across gw_rate tests.

    Sampled once and reused to avoid redundant calls to gw_cbc_statistics.
    """
    return gwrates_instance.gw_cbc_statistics(
        size=N_SAMPLES, batch_size=N_SAMPLES, resume=False, output_jsonfile=False,
    )


class TestGWRATES(CommonTestUtils):
    """Unit tests for the GWRATES class."""

    # -----------------------------------------------------------------------
    # Initialization
    # -----------------------------------------------------------------------

    def test_init(self, gwrates_instance, ler_directory, interpolator_directory):
        """
        Tests
        -----
        - Core attributes are set correctly: ``z_min``, ``z_max``,
          ``event_type``, ``cosmo``, ``npool``, ``ler_directory``.
        - Default ``json_file_names`` contains the expected keys
          (``gwrates_params``, ``gw_param``, ``gw_param_detectable``).
        - ``ler_directory`` exists on disk.
        - GWRATES inherits from ``CBCSourceParameterDistribution``.
        - A user-supplied ``pdet_finder`` is stored on ``self.pdet_finder``
          without triggering gwsnr initialization.
        - A user-supplied astropy cosmology is stored as ``self.cosmo``.
        """
        gwrates = gwrates_instance
        assert gwrates.z_min == DEFAULT_CONFIG["z_min"], f"z_min: expected {DEFAULT_CONFIG['z_min']}, got {gwrates.z_min}"
        assert gwrates.z_max == DEFAULT_CONFIG["z_max"], f"z_max: expected {DEFAULT_CONFIG['z_max']}, got {gwrates.z_max}"
        assert gwrates.event_type == DEFAULT_CONFIG["event_type"], f"event_type: expected {DEFAULT_CONFIG['event_type']}, got {gwrates.event_type}"
        assert gwrates.npool == DEFAULT_CONFIG["npool"], f"npool: expected {DEFAULT_CONFIG['npool']}, got {gwrates.npool}"
        assert gwrates.cosmo is not None, "cosmo must not be None after initialization"
        assert gwrates.ler_directory == ler_directory, f"ler_directory: expected {ler_directory}, got {gwrates.ler_directory}"
        assert os.path.isdir(gwrates.ler_directory), f"ler_directory '{gwrates.ler_directory}' does not exist on disk"

        # json_file_names maps logical names → JSON filenames used for all disk I/O;
        # if any key is missing, the corresponding save/load call will silently fail.
        for key in DEFAULT_JSON_KEYS:
            assert key in gwrates.json_file_names, f"missing default json key: {key}"

        # GWRATES inherits CBC source-param samplers via CBCSourceParameterDistribution
        assert isinstance(gwrates, CBCSourceParameterDistribution), \
            "GWRATES must inherit from CBCSourceParameterDistribution"

        # mock_pdet_finder simply returns 1 for every event, so no gwsnr initialization is triggered
        assert gwrates.pdet_finder is mock_pdet_finder, \
            "pdet_finder was not stored on self.pdet_finder after initialization"

        # custom astropy cosmology must be stored as-is (not copied or replaced)
        cosmo = LambdaCDM(H0=67.4, Om0=0.315, Ode0=0.685, Tcmb0=2.725)
        gwrates_custom = _make_gwrates(interpolator_directory, ler_directory, cosmology=cosmo)
        assert gwrates_custom.cosmo is cosmo, \
            "custom cosmology was not stored as self.cosmo (identity check failed)"

    # -----------------------------------------------------------------------
    # gw_cbc_statistics
    # -----------------------------------------------------------------------

    def test_gw_cbc_statistics(self, gwrates_instance, ler_directory):
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
            - ``True``  → file saved under ``self.json_file_names["gw_param"]`` inside
              ``ler_directory``; file exists and content matches the returned dict.
            - ``False`` → no file written; return value is still a valid dict and
              the file does not exist on disk.
        """
        # --- basic output ---
        gw_param = gwrates_instance.gw_cbc_statistics(
            size=N_SAMPLES, batch_size=N_SAMPLES, resume=False, output_jsonfile=False,
        )
        self._assert_unlensed_cbc_outputs(
            gw_param,
            N_SAMPLES,
            expected_gw_keys=EXPECTED_UNLENSED_GW_KEYS_ALIGNED_SPIN,
        )

        # --- batch_size and resume ---
        # this also test the custom file name mode
        fname_resume = "unit_gw_resume.json"

        # first run: 2 batches of 10 → 20 samples written to file
        gw_param_20 = gwrates_instance.gw_cbc_statistics(
            size=20, batch_size=10, resume=False, output_jsonfile=fname_resume,
        )
        assert len(gw_param_20["zs"]) == 20, \
            f"first run (resume=False, size=20): expected 20 samples, got {len(gw_param_20['zs'])}"

        # second run: resume=True, size=40 → batch_handler loads the 20 existing
        # samples and computes that batches 3 and 4 are still needed → 40 total
        gw_param_40 = gwrates_instance.gw_cbc_statistics(
            size=40, batch_size=10, resume=True, output_jsonfile=fname_resume,
        )
        assert len(gw_param_40["zs"]) == 40, \
            f"second run (resume=True, size=40): expected 40 samples, got {len(gw_param_40['zs'])}"

        # test save file is the same as the output dict
        path = os.path.join(ler_directory, fname_resume)
        assert os.path.isfile(path), \
            f"resume: file '{fname_resume}' not found in ler_directory after resume run"
        gw_param_from_file = get_param_from_json(path)
        self._assert_param_dicts_equal(gw_param_from_file, gw_param_40, label="resume")

        # --- output_jsonfile=True: file stored under the default name ---
        gw_param = gwrates_instance.gw_cbc_statistics(
            size=10, batch_size=10, resume=False, output_jsonfile=True,
        )
        default_fname = gwrates_instance.json_file_names["gw_param"]
        path = os.path.join(ler_directory, default_fname)
        assert os.path.isfile(path), \
            f"output_jsonfile=True: file '{default_fname}' not found in ler_directory"
        # verify that the file content (dict) matches the returned dict
        gw_param_from_file = get_param_from_json(path)
        self._assert_param_dicts_equal(gw_param_from_file, gw_param, label="output_jsonfile=True")

        # --- output_jsonfile=False: no file written, but return value is valid ---
        os.remove(path)
        gw_param_no_file = gwrates_instance.gw_cbc_statistics(
            size=10, batch_size=10, resume=False, output_jsonfile=False,
        )
        self._assert_unlensed_cbc_outputs(
            gw_param_no_file,
            10,
            expected_gw_keys=EXPECTED_UNLENSED_GW_KEYS_ALIGNED_SPIN,
        )
        # verify that the file does not exist
        assert not os.path.isfile(path), \
            "output_jsonfile=False: file should not exist"

    # -----------------------------------------------------------------------
    # gw_rate
    # -----------------------------------------------------------------------

    def test_gw_rate(self, gwrates_instance, gw_param):
        """
        Tests
        -----
        - Returns ``(rate, gw_param_detectable)``; rate is a finite positive float.
        - ``pdet_type='boolean'``: rate = ``normalization_pdf_z * detectable / total``
          where detectable counts events with pdet_net >= pdet_threshold (default 0.5).
          With mock pdet=1, all N_SAMPLES events are detectable.
        - ``pdet_type='probability_distribution'``: rate uses sum(pdet_net) instead
          of a hard threshold.
        - ``pdet_threshold=0.0``: all events pass; ``pdet_threshold=2.0``: none pass,
          rate == 0.
        - ``gw_param=None``: loads parameters from the default JSON file on disk.
        - Invalid ``pdet_type`` raises ``ValueError``.
        """
        # --- boolean mode ---
        rate, detectable_param = gwrates_instance.gw_rate(
            gw_param=gw_param, pdet_type="boolean", output_jsonfile=False,
        )
        assert isinstance(rate, float) and np.isfinite(rate) and rate > 0, \
            f"gw_rate (boolean): expected finite positive float, got {rate}"
        # boolean mode: count events where pdet_net >= pdet_threshold (default 0.5)
        detectable = float(np.sum(np.asarray(gw_param["pdet_net"]) >= 0.5))
        # Monte Carlo rate formula: R = normalization_pdf_z * (detectable / total)
        np.testing.assert_allclose(
            rate,
            gwrates_instance.normalization_pdf_z * detectable / N_SAMPLES,
            rtol=1e-10,
        )
        # with mock pdet=1 and threshold=0.5, all N_SAMPLES events are detectable
        self._assert_unlensed_cbc_outputs(
            detectable_param,
            N_SAMPLES,
            expected_gw_keys=EXPECTED_UNLENSED_GW_KEYS_ALIGNED_SPIN,
        )

        # --- probability_distribution mode ---
        rate, _ = gwrates_instance.gw_rate(
            gw_param=gw_param, pdet_type="probability_distribution", output_jsonfile=False,
        )
        assert isinstance(rate, float) and np.isfinite(rate) and rate > 0, \
            f"gw_rate (probability_distribution): expected finite positive float, got {rate}"
        # probability_distribution mode: sum pdet_net directly instead of thresholding
        detectable = float(np.sum(np.asarray(gw_param["pdet_net"])))
        np.testing.assert_allclose(
            rate,
            gwrates_instance.normalization_pdf_z * detectable / N_SAMPLES,
            rtol=1e-10,
        )

        # --- pdet_threshold: boundary cases ---
        # threshold=0.0 → all events pass (pdet_net=1 >= 0.0)
        rate_all, params_all = gwrates_instance.gw_rate(
            gw_param=gw_param, pdet_threshold=0.0, pdet_type="boolean",
            output_jsonfile=False,
        )
        assert isinstance(rate_all, float) and np.isfinite(rate_all) and rate_all > 0, \
            f"pdet_threshold=0.0: expected finite positive float rate, got {rate_all}"
        assert len(params_all["zs"]) == N_SAMPLES, \
            f"pdet_threshold=0.0: expected all {N_SAMPLES} events detectable, got {len(params_all['zs'])}"

        # threshold=2.0 → no event has pdet_net >= 2.0 → zero detectable, rate = 0
        rate_none, params_none = gwrates_instance.gw_rate(
            gw_param=gw_param, pdet_threshold=2.0, pdet_type="boolean",
            output_jsonfile=False,
        )
        assert len(params_none["zs"]) == 0, \
            f"pdet_threshold=2.0: expected 0 detectable events, got {len(params_none['zs'])}"
        assert rate_none == 0.0, \
            f"pdet_threshold=2.0: expected rate=0, got {rate_none}"

        # --- gw_param=None: load from file ---
        # write the default file first so gw_rate can find it
        # save the detectable param in the file
        # check output file is the same as the returned dict
        fname = 'unit_gw_param.json'
        gwrates_instance.gw_cbc_statistics(
            size=N_SAMPLES, batch_size=N_SAMPLES, resume=False, output_jsonfile=True,
        )
        rate, loaded_param = gwrates_instance.gw_rate(
            gw_param=None, output_jsonfile=fname,
        )
        assert isinstance(rate, float) and np.isfinite(rate) and rate > 0, \
            f"gw_rate (gw_param=None): expected finite positive float, got {rate}"
        self._assert_unlensed_cbc_outputs(
            loaded_param,
            N_SAMPLES,
            expected_gw_keys=EXPECTED_UNLENSED_GW_KEYS_ALIGNED_SPIN,
        )
        # check output file is the same as the returned dict
        loaded_param_from_file = get_param_from_json(os.path.join(gwrates_instance.ler_directory, fname))
        self._assert_param_dicts_equal(loaded_param_from_file, loaded_param, label="gw_rate (gw_param=None)")

        # --- invalid pdet_type raises ValueError ---
        with pytest.raises(ValueError, match="pdet_type"):
            gwrates_instance.gw_rate(
                gw_param=gw_param, pdet_type="not_valid", output_jsonfile=False,
            )

    # -----------------------------------------------------------------------
    # selecting_n_gw_detectable_events
    # -----------------------------------------------------------------------

    def test_selecting_n_detectable_events(self, gwrates_instance):
        """
        Tests
        -----
        - ``stopping_criteria=None``: accumulates batches until the number of
          detectable events exceeds ``size``; meta-data file has the expected
          keys and the reported rate matches the function return value.
        - ``trim_to_size=False``: returns all collected events (final_size >= size);
          stored JSON file matches the returned dict.
        - ``resume=True`` with the same output file (already has final_size events):
          no new sampling; ``trim_to_size=True`` trims randomly to ``new_size``
          and every value in the result is a member of the original collection.
        - ``stopping_criteria=dict(...)``: stops when the cumulative rate has
          converged (relative difference of the last 4 batches < 0.5%) and
          collected size exceeds ``size``.
        """

        # trim_to_size=False.
        # stopping_criteria=None: accumulate until >= size detectable events are found
        size = 500
        batch_size = 10000
        output_jsonfile = "unit_n_detectable.json"
        meta_data_file = "unit_meta_gw.json"
        rate, param = gwrates_instance.selecting_n_gw_detectable_events(
            size=size, 
            batch_size=batch_size, 
            stopping_criteria=None,
            pdet_threshold=0.5, 
            resume=False, 
            trim_to_size=False,
            output_jsonfile=output_jsonfile,
            meta_data_file=meta_data_file,
        )

        # get the meta file 
        # check the keys ['events_total', 'detectable_events', 'total_rate']
        # the last total_rate should be the same as the rate returned by the function
        # rate should be finite and positive
        meta_data = get_param_from_json(os.path.join(gwrates_instance.ler_directory, meta_data_file))
        assert isinstance(meta_data, dict) and len(meta_data) == 3, \
            f"meta_data: expected dict with 3 keys, got {len(meta_data)}"
        assert 'events_total' in meta_data and 'detectable_events' in meta_data and 'total_rate' in meta_data, \
            f"meta_data: expected keys ['events_total', 'detectable_events', 'total_rate'], got {list(meta_data.keys())}"
        assert isinstance(meta_data['total_rate'][-1], float) and np.isfinite(meta_data['total_rate'][-1]) and meta_data['total_rate'][-1] > 0, \
            f"meta_data: expected finite positive float total_rate, got {meta_data['total_rate'][-1]}"
        assert meta_data['total_rate'][-1] == rate, \
            f"meta_data: expected total_rate to be the same as the rate returned by the function, got {meta_data['total_rate'][-1]} and {rate}"

        # check the param dict
        # since trim_to_size=False, the param dict should have all the sample events with final_size > size.
        final_size = meta_data['detectable_events'][-1]
        self._assert_unlensed_cbc_outputs(
            param,
            final_size,
            expected_gw_keys=EXPECTED_UNLENSED_GW_KEYS_ALIGNED_SPIN,
        )

        # stored dict in the json file should be the same as the param returned by the function
        dict_from_file = get_param_from_json(os.path.join(gwrates_instance.ler_directory, output_jsonfile))
        self._assert_param_dicts_equal(dict_from_file, param, label="selecting_n_detectable_events")

        # resume from the same file which already has final_size events (>= new_size);
        # no new sampling occurs; trim_to_size=True randomly subsamples to new_size.
        new_size = 600 
        batch_size = 10000
        new_rate, new_param = gwrates_instance.selecting_n_gw_detectable_events(
            size=new_size, 
            batch_size=batch_size, 
            stopping_criteria=None,
            pdet_threshold=0.5, 
            resume=True, 
            trim_to_size=True,
            output_jsonfile=output_jsonfile,
            meta_data_file=meta_data_file,
        )
        # check the param dict
        self._assert_unlensed_cbc_outputs(
            new_param,
            new_size,
            expected_gw_keys=EXPECTED_UNLENSED_GW_KEYS_ALIGNED_SPIN,
        )
        # resumed from the same file: new_param is a random subset of param;
        # every value in new_param must appear in param (trim_to_size uses np.random.choice)
        for key in EXPECTED_PARAM_KEYS:
            assert np.isin(new_param[key], param[key]).all(), \
                f"selecting_n_detectable_events (resume+trim): '{key}' contains values not in original param"

        # use stopping_criteria=dict(relative_diff_percentage=0.5, number_of_last_batches_to_check=4)
        # check meta data file for the last total_rate. the last 4 values should relative difference less than 0.5%
        # collected size should be more than new_size.
        stopping_criteria = dict(relative_diff_percentage=0.5, number_of_last_batches_to_check=4)
        new_rate, new_param = gwrates_instance.selecting_n_gw_detectable_events(
            size=new_size, 
            batch_size=batch_size, 
            stopping_criteria=stopping_criteria,
            pdet_threshold=0.5, 
            resume=False, 
            trim_to_size=False,
            output_jsonfile="unit_n_detectable_resume_stopping.json",
            meta_data_file="unit_meta_gw_resume_stopping.json",
        )
        meta_data = get_param_from_json(os.path.join(gwrates_instance.ler_directory, "unit_meta_gw_resume_stopping.json"))
        idx_converged = [-4, -3, -2, -1]
        rates_converged = np.array(meta_data['total_rate'])[idx_converged]
        percentage_diff = (
                    np.abs((rates_converged - rates_converged[-1]) / rates_converged[-1])
                    * 100
                )
        assert np.all(percentage_diff < 0.5), \
            f"selecting_n_detectable_events: relative difference of total rate for the last 4 cumulative batches should be less than 0.5% for stopping criteria to be met"

        num_collected = meta_data['detectable_events'][-1]
        assert num_collected > new_size, \
            f"selecting_n_detectable_events: collected size should be more than new_size, got {num_collected} and {new_size}"

    # -----------------------------------------------------------------------
    # rate_function and _load_param
    # -----------------------------------------------------------------------

    def test_utilities(self, gwrates_instance):
        """
        Tests
        -----
        - ``rate_function`` returns ``normalization_pdf_z * detectable / total``;
          output is a finite positive float.
        - ``_load_param`` returns a copy (not the same dict object) when given a
          dict, so downstream mutations do not affect the caller's dict.
        """
        # --- rate_function math ---
        detectable, total = 100, 1000
        rate = gwrates_instance.rate_function(detectable, total, verbose=False)
        assert isinstance(rate, float) and np.isfinite(rate) and rate > 0, \
            f"rate_function: expected finite positive float, got {rate}"
        # Monte Carlo rate formula: R = normalization_pdf_z * (detectable / total)
        np.testing.assert_allclose(
            rate, gwrates_instance.normalization_pdf_z * detectable / total, rtol=1e-12,
        )

        # --- _load_param returns a copy ---
        original = dict(
            zs=np.array([1.0, 2.0]),
            pdet_net=np.array([0.6, 0.9]),
        )
        loaded = gwrates_instance._load_param(original, param_type="gw")

        # must be a different object so downstream mutations don't affect the caller
        assert loaded is not original, "_load_param returned the same dict object, not a copy"
        assert set(loaded.keys()) == set(original.keys()), \
            f"_load_param: key mismatch: expected {set(original.keys())}, got {set(loaded.keys())}"
        np.testing.assert_array_equal(loaded["zs"], original["zs"])
        np.testing.assert_array_equal(loaded["pdet_net"], original["pdet_net"])
