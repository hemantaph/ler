# Integration tests

These tests exercise end-to-end workflows: populations, detection
probabilities, rates, selection functions, and the optional LeR speed
benchmark.

## How to run

Run active integration tests:

```bash
python -m pytest tests/integration/ --tb=short
```

Include slow integration tests:

```bash
python -m pytest tests/integration/ --run-slow --tb=short
```

Run only slow integration tests:

```bash
python -m pytest tests/integration/ -m slow --tb=short
```

## What is tested

| File | What it checks | Default run? |
|------|----------------|--------------|
| `test_gwrates_populations.py` | `GWRATES` BBH/BNS/NSBH `gw_cbc_statistics`; rates versus `normalization_pdf_z`; BBH hyperparameter checks for `H0`, `R0`, and `alpha_1`. | Yes |
| `test_gwrates_selection_function.py` | `GWRATES` selection function `xi(Lambda) = mean(pdet_net)` with custom `mass_1_source` from `broken_powerlaw_plus_2peaks_rvs`. | Yes |
| `test_ler_populations.py` | `LeR` EPL/SIE/SIS lens types; unlensed/lensed statistics and rates; `phistar` bracket in the velocity-dispersion prior. | Yes |
| `test_ler_speed_benchmark.py` | JIT `npool=6` versus no-JIT subprocesses; mean of 3 passes through unlensed `size=100` and lensed `size=10`. | No; marked `slow` |

## Shared fixtures

Integration tests use the same fixtures from `tests/conftest.py` as the unit
tests:

- `interpolator_directory`: extracted bundled interpolator data.
- `ler_directory`: temporary output directory for generated JSON files and
  benchmark cache directories.

## Detector and sample-size notes

`test_gwrates_populations.py`, `test_ler_populations.py`, and
`test_gwrates_selection_function.py` use real `gwsnr` detection probabilities.
See `tests/integration/gwrates_integration_settings.md` for detector defaults,
per-event-type overrides, and Monte Carlo sizes.

## CI

Default CI keeps `slow` tests deselected, but the population and
selection-function integration tests are included.
