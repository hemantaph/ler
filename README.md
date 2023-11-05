# LeR
[![DOI](https://zenodo.org/badge/626733473.svg)](https://zenodo.org/badge/latestdoi/626733473) [![PyPI version](https://badge.fury.io/py/ler.svg)](https://badge.fury.io/py/ler) [![DOCS](https://readthedocs.org/projects/ler/badge/?version=latest)](https://ler.readthedocs.io/en/latest/)

`LeR` is a statistical-based python package whose core function is designed for the computation of detectable rates pertaining to both lensed and unlensed gravitational wave (GW) events. This calculation intricately hinges upon the interplay of various components within the package, which can be categorized into three primary segments: 1. Sampling the properties of compact-binary sources, 2. Sampling characteristics of lens galaxies, and 3. Solving the lens equation to derive image attributes of the source. The holistic functionality of the package is built upon leveraging array operations and linear algebra from the `numpy` library, interpolation from `scipy`, and the `multiprocessing` capability inherent to Python. This design optimizes both speed and functionality while upholding user-friendliness. The architecture of the "LeR" API is deliberately organized such that each distinct functionality holds its own significance in scientific research. Simultaneously, these functionalities seamlessly integrate and can be employed collectively based on specific research requirements. Key features of `LeR` and its dependencies can be summarized as follows,

- Detectable merger rates: 
    * The calculation depends not only on simulated event properties but also on GW detector detectability. We compute optimal signal-to-noise ratios (SNRs) for simulated events, which can be computationally intensive. `LeR` mitigates this using [`gwsnr`](https://github.com/hemantaph/gwsnr) for efficient  and rapid SNR calculation. `gwsnr` enables rate calculation for current and future detectors with customizable sensitivities.
    * The merger rates of both the simulated unlensed and lensed events can be computed and subsequently compared. 
- Sampling GW sources:
    * The distribution of the source's red-shift is based on the merger rate density of compact binaries, including BBHs, BNSs etc. The code is meticulously structured to facilitate straightforward incorporation of future updates or additional distributions of such sources by users.
    * The sampling of intrinsic and extrinsic parameters of gravitational wave sources is conducted employing the prior distributions encompassed within the `gwcosmo` and `bilby` Python packages. Prior to parameterizing the rate calculation, users retain the flexibility to manually substitute any relevant parameters as needed.
- Sampling of lens galaxies:
    * The Lens distribution follows [(Oguri et al. 2018](https://arxiv.org/abs/1807.02584). It depends on the sampled source red-shifts and also on the optical depth.
    * `LeR` employs the Elliptical Power Law model with the external shear (EPL+Shear) model for sampling other galaxy features, which is available in the `Lenstronomy` package.
    * Rejection sampling is applied on the above samples on condition that whether the event is strongly lensed or not.
- Generation of image properties:
    * Source position is sampled from the caustic in the source plane.
    * Sampled lens properties and source position is fed in `Lenstronomy` to generate properties of the images.
    * Properties like magnification and time delay are essential as it modifies the source signal strength, changing the SNR and detection ability.
    * `LeR` can handle both super-threshold and sub-threshold events in picking detectable events and rate computation.

The `LeR` software has been developed to cater to the requirements of both the LIGO scientific collaboration and research scholars engaged in astrophysics studies. It is currently used in generating detectable lensing events and GW lensing rates with the available information on current and future detectors. The results will predict the feasibility of various detectors for detecting and studying such lensing events. Statistics generated from `LeR` will be used in event validation of the ongoing effort to detect lensed gravitational waves. Lastly, `LeR` was designed with upgradability in mind to include additional statistics as required by the related research.

# Installation

Follow the installation instruction at [ler.readthedoc](https://ler.readthedocs.io/en/latest/installation.html)




