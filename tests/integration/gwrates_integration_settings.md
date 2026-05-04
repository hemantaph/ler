# GWRATES integration settings

This note records the detector, waveform, and sample-size assumptions used by
the integration tests. It is mainly for readers who want to understand why the
rate floors in `test_gwrates_populations.py` are deliberately loose.

## How to run the related tests

Run the integration tests that use these settings:

```bash
python -m pytest tests/integration/test_gwrates_populations.py tests/integration/test_gwrates_selection_function.py --tb=short
```

Run all active integration tests:

```bash
python -m pytest tests/integration/ --tb=short
```

## What this file explains

- Which detector network `GWRATES` uses when `psds=None` and `ifos=None`.
- Which `GWRATES` to `GWSNR` defaults matter for these tests.
- Which event-type-specific overrides are used in `test_gwrates_populations.py`.
- Why the integration tests use loose rate floors.
- Which Monte Carlo sizes are used by each integration test.

## Default network

When `psds=None` and `ifos=None`, `gwsnr` asks Bilby for its default PSDs:

| IFO | PSD file |
|-----|----------|
| L1 | `aLIGO_O4_high_asd.txt` |
| H1 | `aLIGO_O4_high_asd.txt` |
| V1 | `AdV_asd.txt` |

So the default `GWRATES` integration setup is an O4-style LIGO high range plus
Virgo design/AdV three-detector network: `L1`, `H1`, `V1`.

The population-test rate floors are loose sanity checks for that setup:

| Event type | Floor |
|------------|-------|
| BBH | `1e2 yr^-1` |
| BNS | `1e5 yr^-1` |
| NSBH | `1e4 yr^-1` |

## Default GWRATES to GWSNR arguments

These defaults mirror `ler/rates/gwrates.py::_gwsnr_initialization` unless an
individual test overrides them.

| Setting | Default |
|---------|---------|
| `snr_method` | `interpolation_aligned_spins` |
| `snr_type` | `optimal_snr` |
| `pdet_kwargs` | `snr_th=10`, `snr_th_net=10`, `pdet_type='boolean'`, `distribution_type='noncentral_chi2'`, `include_optimal_snr=False`, `include_observed_snr=False` |
| `mtot_min` / `mtot_max` | `1.0` / `500.0` |
| `ratio_min` / `ratio_max` | `0.1` / `1.0` |
| `spin_max` | `0.99` |
| `mtot_resolution` / `ratio_resolution` / `spin_resolution` | `200` / `20` / `10` |
| `batch_size_interpolation` | `1_000_000` |
| `interpolator_dir` | From `GWRATES(interpolator_directory=...)` |
| `create_new_interpolator` | `False` |
| `waveform_approximant` | `IMRPhenomD` |
| `frequency_domain_source_model` | `lal_binary_black_hole` |
| `sampling_frequency` / `minimum_frequency` | `2048 Hz` / `20 Hz` |
| `reference_frequency` | `None` |
| `duration_max` / `duration_min` / `fixed_duration` | `None` |
| `mtot_cut` | `False` |
| `psds` / `ifos` | `None` / `None` |
| `noise_realization` | `None` |
| `ann_path_dict` | `None` |
| `snr_recalculation` | `False` |
| `snr_recalculation_range` | `[6, 14]` |
| `snr_recalculation_waveform_approximant` | `IMRPhenomXPHM` |

Tests construct `GWRATES` with `verbose=False`, `gwsnr_verbose=False`, and
`multiprocessing_verbose=False` where applicable.

## Population-test overrides

`tests/integration/test_gwrates_populations.py` keeps the default BBH-oriented
path for BBH and uses a faster, non-degenerate detector setup for compact
NSBH/BNS draws.

| Event type | Extra kwargs | Reason |
|------------|--------------|--------|
| BBH | `spin_precession=True`, `waveform_approximant='IMRPhenomXPHM'`, `snr_recalculation=True`, `snr_recalculation_range=[6, 14]` | Exercises precessing BBH with SNR refinement. |
| NSBH | `snr_method='interpolation_no_spins'`, `ifos=['ET', 'CE']`, `snr_type='optimal_snr'`, `spin_zero=True`, `waveform_approximant='IMRPhenomD'` | Keeps the NSBH detection probability away from a near-zero default-grid corner. |
| BNS | Same pattern as NSBH. | Fast BNS path without the BBH spin grid. |

`gw_rate` uses the default `pdet_type='boolean'` in these tests.

## Selection function

`tests/integration/test_gwrates_selection_function.py` uses BBH precession and
SNR recalculation, then injects two fixed `broken_powerlaw_plus_2peaks_rvs`
hyperparameter dictionaries through `GWRATES.mass_1_source`. It checks that
`xi(Lambda) = mean(pdet_net)` is finite, lies in `(0, 1]`, and changes between
the two hyperparameter points.

## Monte Carlo sizes

| Test file | Size | Batch size | Default CI? |
|-----------|------|------------|-------------|
| `test_gwrates_populations.py` | `20000` | `10000` | Yes |
| `test_ler_populations.py` | `20000` | `10000` | Yes |
| `test_gwrates_selection_function.py` | `8000` | `4000` | Yes |
| `test_ler_speed_benchmark.py` | `100` unlensed, `10` lensed | Same as size | No; marked `slow` |

The default pytest run deselects `slow` tests, so the speed benchmark only runs
with `--run-slow` or `-m slow`.
