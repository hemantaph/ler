# ler Test Suite

Comprehensive test suite for validating ler functionality across unit and integration scenarios. The test suite ensures reliability, accuracy, and performance of gravitational wave, lensed and unlensed, events simulations and rate calculations. 

## Test Organization

### Unit Tests (`unit/`)
Component-level tests validating individual ler components and backends:

- `gw_source_population`: Tests for source population generation and sampling.
- `lens_galaxy_population`: Tests for optical depth calculation and lens galaxy population generation and sampling.
- `image_properties`: Tests for lensing magnification, time delay, and image properties calculations.
- `detection_probability` with `LeR` and `GWRATES`: Tests for detection probability calculations using gwsnr and other methods.

what to be tested? input output test. interpolator generation for priors and other functions. custom priors and custom functions.

### Integration Tests (`integration/`)
End-to-end workflow tests for astrophysical applications:
- `BBH source and EPL+Shear lens simulation and rates`: Simulate BBH population with EPL+Shear lens model, compute lensed and unlensed rates, and validate results against expected values.
- `BNS source and EPL+Shear lens simulation and rates`: Simulate BNS population with EPL+Shear lens model, compute lensed and unlensed rates, and validate results against expected values.
- `BBH source and SIE lens simulation and rates`: Simulate BBH population with SIE lens model, compute lensed and unlensed rates, and validate results against expected values.
- `BBH source and SIS lens simulation and rates`: Simulate BBH population with SIS lens model, compute lensed and unlensed rates, and validate results against expected values.