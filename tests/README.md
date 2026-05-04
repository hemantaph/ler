# ler test suite

This directory contains the unit and integration tests for `ler`.

## How to run

Run the default test suite:

```bash
python -m pytest tests/ --tb=short
```

Run only unit tests:

```bash
python -m pytest tests/unit/ --tb=short
```

Run only integration tests:

```bash
python -m pytest tests/integration/ --tb=short
```

Run slow tests as well:

```bash
python -m pytest tests/ --run-slow --tb=short
```

Run only tests marked `slow`:

```bash
python -m pytest tests/ -m slow --tb=short
```

## What is tested

| Path | Purpose | Default run? |
|------|---------|--------------|
| `tests/unit/` | Component checks for source populations, lens populations, image properties, `GWRATES`, `LeR`, and shared utilities. | Yes |
| `tests/integration/` | End-to-end population, rate, and selection-function checks using real `gwsnr` where applicable. | Yes |
| `@pytest.mark.slow` tests | Wall-clock/JIT comparisons, heavy samplers, and the subprocess LeR speed benchmark. | No; opt in |

The main coverage areas are:

- Input/output shape and finite-value checks for priors, samplers, interpolators, and rate outputs.
- Custom priors and custom functions where the public API supports them.
- Unlensed rate calculation through `GWRATES`.
- Lensed and unlensed rate calculation through `LeR`.
- Lens-galaxy optical depth, cross sections, source sampling, image positions, magnifications, and time delays.
- Integration-level population and selection-function behavior with `gwsnr`.

## Integration tests

| File | What it checks |
|------|----------------|
| `test_gwrates_populations.py` | BBH/BNS/NSBH population rates and BBH hyperparameter checks. |
| `test_ler_populations.py` | `LeR` lens-type sweep for EPL/SIE/SIS plus `phistar` uncertainty. |
| `test_gwrates_selection_function.py` | Selection function `xi(Lambda) = mean(pdet_net)` for two mass-model hyperparameter points. |
| `test_ler_speed_benchmark.py` | `slow`: JIT `npool=6` versus `NUMBA_DISABLE_JIT=1`, using subprocesses and small benchmark sizes. |

See `tests/integration/README.md` for the integration test map and
`tests/integration/gwrates_integration_settings.md` for detector defaults,
sample sizes, and `gwsnr` settings.

## Slow tests

The default pytest collection deselects `@pytest.mark.slow` in
`tests/conftest.py`. Slow tests include:

- `test_njit_speed`
- `test_njit_speed_gw_parameters_rvs`
- `test_cross_section_comparison_with_lenstronomy`
- `test_cross_section_based_sampler_variants_heavy`
- `test_njit_vs_lenstronomy_speed`
- `tests/integration/test_ler_speed_benchmark.py`

Slow tests are mostly timing-sensitive and are meaningful only after relevant
Numba caches are warm.

## Cache and output directories

`tests/conftest.py` provides shared fixtures:

- `interpolator_directory`: extracts bundled `interpolator_json.zip` into
  `tests/.cache/interpolators/` once per test session. Override with
  `LER_TEST_INTERPOLATOR_CACHE_DIR` if needed.
- `ler_directory`: temporary output directory for JSON files written during
  tests.

## CI

GitHub Actions (`.github/workflows/tests.yml`) uses the default test selection:
unit tests and active integration tests run, while `slow` tests are deselected.

Sphinx deployment is separate: `.github/workflows/documentation.yml`.
