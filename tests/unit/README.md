# Unit Tests

Component-level tests covering source-population sampling, lensing statistics,
unlensed (`GWRATES`) and lensed+unlensed (`LeR`) rate calculation, and image
properties. Each test module targets one class and is largely independent.

## How to run

Run all active unit tests:

```bash
python -m pytest tests/unit/ --tb=short
```

Include slow unit tests:

```bash
python -m pytest tests/unit/ --run-slow --tb=short
```

Run only slow unit tests:

```bash
python -m pytest tests/unit/ -m slow --tb=short
```

Show the slowest unit tests:

```bash
python -m pytest tests/unit/ --durations=25
```

## What is tested

Tests are executed in dependency order (light sampling modules first,
the full pipeline last) enforced by `conftest.py`.

### `test_cbc_source_redshift_distribution.py`

Validates `CBCSourceRedshiftDistribution`, which draws source redshifts from
a user-specified merger-rate density.

| Test | What is checked |
|------|-----------------|
| `test_init` | Correct `z_min`/`z_max`, merger-rate attribute, and `normalization_pdf_z` for BBH, BNS, NSBH. |
| `test_merger_rate_density_and_source_sampling` | Sampled redshifts are finite, positive, and within the specified range. |
| `test_cosmological_functions_interpolator` | Luminosity- and angular-diameter-distance interpolators return finite, positive arrays. |
| `test_custom_merger_rate_density_parameters` | Custom MRD supplied as a string name or callable produces valid samples. |
| `test_sampling_reproducibility` | Same random seed → identical output arrays. |
| `test_invalid_*_raises` | `ValueError` on unknown event type or unknown MRD function name. |
| `test_njit_speed` *(slow)* | Numba-accelerated sampler is faster than the no-JIT baseline after warm-up. |

---

### `test_cbc_source_parameter_distribution.py`

Validates `CBCSourceParameterDistribution`, which draws the full set of GW
source parameters (masses, spins, sky position, redshift, luminosity distance).

| Test | What is checked |
|------|-----------------|
| `test_output_sanity_spin_zero` | All expected parameter keys present; arrays are finite and within physical bounds; spin keys absent when `spin_zero=True`. Parametrized over BBH, BNS, NSBH. |
| `test_bbh_precessing_spin_output_sanity` | Precessing-spin keys (`tilt_*`, `phi_*`) present and within bounds. |
| `test_custom_all_gw_priors_sanity` | All default priors replaced with custom lambdas; output still has correct shape and finite values. |
| `test_custom_prior_parameters` | Narrowing the mass-ratio prior bounds is reflected in the sampled `mass_2_source`. |
| `test_sampling_reproducibility_with_reseeding` | Same seed → identical arrays. |
| `test_sample_gw_parameters_fixed_param_override` | A pre-drawn `zs` array passed via `param=` is preserved while other parameters are resampled. |
| `test_invalid_event_type_raises` | `ValueError` on unknown `event_type`. |
| `test_gw_parameters_rvs_njit_output_sanity` | Numba-backed sampler returns the same keys and valid values as the Python path. |
| `test_njit_speed_gw_parameters_rvs` *(slow)* | Numba path outperforms the no-JIT baseline. |

---

### `test_gwrates.py`

Validates `GWRATES`, the unlensed GW detection-rate calculator.

| Test | What is checked |
|------|-----------------|
| `test_init` | Core attributes, `json_file_names`, class inheritance, custom `pdet_finder` and cosmology stored correctly. |
| `test_gw_cbc_statistics` | Output keys and finite values; batched sampling with resume (20 → 40 events); `output_jsonfile` on/off. |
| `test_gw_rate` | Rate formula (`normalization_pdf_z × detectable / total`) for `boolean` and `probability_distribution` detection modes; boundary `pdet_threshold` values; loading parameters from file; invalid `pdet_type` raises `ValueError`. |
| `test_selecting_n_detectable_events` | Accumulate exactly *N* detectable events across batches; resume + trim to a random subset; convergence via `stopping_criteria`. |
| `test_utilities` | `rate_function` arithmetic; `_load_param` returns an independent copy. |

---

### `test_optical_depth.py`

Validates `OpticalDepth`, which computes the strong-lensing optical depth and
draws lens parameters.

| Test | What is checked |
|------|-----------------|
| `test_optical_depth_vs_sis_analytic_baseline` | Computed optical depth is within ±50% of the SIS analytic formula at several source redshifts. |
| `test_optical_depth_output_value_test` | Lens-redshift PDF/RVS, source-redshift RVS, velocity-dispersion, axis-ratio, and cross-section interfaces all return valid outputs for all supported lens types. |
| `test_cross_section_comparison_with_lenstronomy` *(slow)* | njit cross-section agrees with lenstronomy to within 1% for EPL+shear. |

---

### `test_lens_galaxy_parameter_distribution.py`

Validates `LensGalaxyParameterDistribution`, which samples joint lens + source
parameters conditioned on strong lensing.

| Test | What is checked |
|------|-----------------|
| `test_init` | `z_min`, `z_max`, `lens_type`, `npool`, cosmology, and inheritance from `OpticalDepth`. |
| `test_strongly_lensed_source_redshift` | Conditioned source redshifts are finite, positive, and within bounds. |
| `test_cross_section_based_sampler_variants` | Cross-section EPL+shear sampling with `importance_sampler_partial` only (`n_prop` from module-level `N_PROP_CROSS_SECTION`). Default / CI run. |
| `test_cross_section_based_sampler_variants_heavy` *(slow)* | Same assertions for `rejection_sampler_full`, `importance_sampler_full`, and `rejection_sampler_partial` (larger `n_prop` for those modes). |
| `test_intrinsic_lens_parameter_sampling` | Intrinsic EPL+shear sampler returns the expected parameter keys with valid arrays. |
| `test_njit_speed` *(slow)* | Numba path outperforms no-JIT baseline; `npool=6` outperforms `npool=1`. |

---

### `test_image_properties.py`

Validates `ImageProperties`, which solves the EPL+shear lens equation and
computes image magnifications, time delays, and related observables.

| Test | What is checked |
|------|-----------------|
| `test_init` / `test_standalone_image_properties_init` | Attributes set correctly; class works independently of `LensGalaxyParameterDistribution`. |
| `test_custom_cosmology_propagates_to_da` | A custom astropy cosmology is forwarded to the angular-diameter-distance interpolator. |
| `test_image_properties_njit_vs_lenstronomy_output_sanity` | njit and lenstronomy backends agree on image count, magnification signs, and time-delay ordering. |
| `test_caustic_njit_vs_lenstronomy` | Caustic boundary points from both backends are finite and agree within 1%. |
| `test_source_sampling_inside_caustic` | Sampled source positions lie strictly inside the caustic. |
| `test_lens_equation_njit_vs_lenstronomy` | For fixed source positions, both backends find the same number of images and agree on positions to < 1%. |
| `test_recover_redundant_parameters` | `theta_E`, `n_images`, source-frame masses, and luminosity distance are correctly recovered from the compact lensed-param dict. |
| `test_produce_effective_params` | Effective luminosity distance = `d_L / sqrt(|μ|)` and effective geocentric time = `t + Δt` for the first image. |
| `test_njit_vs_lenstronomy_speed` *(slow)* | njit solver is faster than the lenstronomy loop. |

---

### `test_ler.py`

Validates `LeR`, the end-to-end pipeline combining unlensed and lensed GW event
statistics and rate estimation.

| Test | What is checked |
|------|-----------------|
| `test_init` | All core attributes, `json_file_names`, directory created on disk, class inheritance, custom `pdet_finder` and cosmology. |
| `test_unlensed_cbc_statistics` | Output keys and finite values; batch + resume (20 → 40 events); `output_jsonfile` on/off. |
| `test_lensed_cbc_statistics` | Same as unlensed, plus: `pdet_net` is 2-D `(n_events, n_max_images)` with NaN padding for missing images. |
| `test_unlensed_rate` | Rate formula for `boolean` and `probability_distribution` modes; boundary `pdet_threshold` values; loading from file; invalid `pdet_type` raises `ValueError`. |
| `test_lensed_rate` | Rate formula with per-image threshold `[0.5, 0.5]`; impossible threshold gives rate = 0; loading from file. |
| `test_utilities` | `rate_function` arithmetic for both `'unlensed'` and `'lensed'` modes; `_load_param` returns an independent copy. |
| `test_selecting_n_unlensed_detectable_events` | Accumulate *N* detectable unlensed events; resume + trim; convergence via `stopping_criteria`. |
| `test_selecting_n_lensed_detectable_events` | Same for lensed events (smaller batch size; 2% convergence threshold because lensed rates converge more slowly). |

---

## `slow` marker (wall-clock / heavy)

`@pytest.mark.slow` tests compare wall-clock times (e.g. Numba vs no-JIT, or
multi-process vs single-process) and/or run slower sampling variants that are
not needed on every CI run. The default pytest collection deselects them
(see `tests/conftest.py`); they are only meaningful after the Numba JIT cache
is warm.

Files that define `slow` tests (in addition to any mentioned in the module
tables above):

| File | `slow` tests |
|------|----------------|
| `test_cbc_source_redshift_distribution.py` | `test_njit_speed` |
| `test_cbc_source_parameter_distribution.py` | `test_njit_speed_gw_parameters_rvs` |
| `test_optical_depth.py` | `test_cross_section_comparison_with_lenstronomy` |
| `test_lens_galaxy_parameter_distribution.py` | `test_cross_section_based_sampler_variants_heavy`, `test_njit_speed` |
| `test_image_properties.py` | `test_njit_vs_lenstronomy_speed` |

---

## Shared infrastructure

**`tests/conftest.py`** provides two session-scoped fixtures used by all tests:

- `interpolator_directory` — extracts the bundled `interpolator_json.zip` once
  per session into a persistent cache directory, avoiding repeated interpolator
  construction.
- `ler_directory` — a temporary directory for JSON output files written during
  tests.

**`tests/tests_utils.py`** provides `CommonTestUtils`, a base class inherited by every
`Test*` class, with helpers for validating arrays (`_assert_array_valid`),
parameter dictionaries (`_assert_param_dict_valid`, `_assert_param_dicts_equal`),
`FunctionConditioning` objects (`_assert_valid_object`), and the module-level
timing helper `median_call_time` for wall-clock / `slow`-style tests.

## Integration tests

Full-pipeline integration checks live under `tests/integration/` (see `tests/integration/README.md` and `tests/README.md`). Default CI keeps `slow` tests deselected while population and selection pipelines stay on.
