# LeR
[![DOI](https://zenodo.org/badge/626733473.svg)](https://zenodo.org/badge/latestdoi/626733473) [![PyPI version](https://badge.fury.io/py/ler.svg)](https://badge.fury.io/py/ler) [![DOCS](https://readthedocs.org/projects/ler/badge/?version=latest)](https://ler.readthedocs.io/en/latest/)

<p align="center">
  <img src="lerlogo.png" alt="Your Logo" width="30%">
</p>

## Installation

```
pip install ler
```

## About

<p align="center">
  <img src="sl.png" alt="Your Logo" width="100%">
</p>

<span style="font-family:cambria; color:maroon;">LeR</span> is a statistical-based python package whose core function is designed for the computation of detectable rates pertaining to both lensed and unlensed gravitational wave (GW) events. This calculation intricately hinges upon the interplay of various components within the package, which can be categorized into four primary segments: *1.* Sampling the properties of compact-binary sources, *2.* Sampling the characteristics of lens galaxies, *3.* Solving the lens equations to derive the properties of the resultant images, and finally, *4.* Computing the rates of detectable GW events. The holistic functionality of the package is built upon leveraging array operations and linear algebra from the *numpy* library, complemented by interpolation methods from *scipy* and Python’s native *multiprocessing* capabilities. The software's efficiency is notably enhanced by the *numba* library's Just-In-Time (*njit*, n: no-python mode) compilation, which dynamically converts Python and NumPy code into machine code, thereby optimizing performance in numerical computations involving extensive loops and array operations. Predominantly, <span style="font-family:cambria; color:maroon;">LeR</span> employs the inverse transform sampling method, utilizing interpolated inverse cumulative distribution functions (CDFs). This approach eschews the need for rejection sampling, which often requires numerous loops and the handling of relatively tedious probability density functions (PDFs) of the sampled parameters. The overall design of the package optimizes both speed and functionality while upholding user-friendliness. The architecture of the <span style="font-family:cambria; color:maroon;">LeR</span> API is deliberately organized such that each distinct functionality holds its own significance in scientific research. Simultaneously, these functionalities seamlessly integrate and can be employed collectively based on specific research requirements. Key features of <span style="font-family:cambria; color:maroon;">LeR</span> and its dependencies can be summarized as follows:

- Sampling Gravitational Wave (GW) Source Properties:
    * For the unlensed events, The sampling distribution $(\,R_m^U(z_s)\,)$ for the source's redshift $(\,z_s\,)$ is derived from the merger rate density of compact binaries, which, in turn, is based on the star formation rate. The code is meticulously designed to enable the straightforward integration of future updates or user-specified distributions of these sources.
    * The process of sampling both intrinsic and extrinsic parameters of GW sources, represented by $\theta_i$, utilizes the prior distributions $(P(\theta_i))$ available within the *gwcosmo* and *bilby* Python packages. Before parameterizing the rate calculation, the software allows users the flexibility to manually alter any pertinent parameters as required. The said feature is employed when one doesn't specifically have a prior distribution for a particular parameter.

- Sampling of lens galaxies attributes and source red-shifts:
    * For lensed case, the source redshift $(z_s)$ is sampled under the strong lensing condition $(\text{SL})$, based on the precomputed probability of strong lensing with source at $z_s$ $(\text{optical depth: }P(z_s|\text{SL})\,)$. This probability can be recalculated for specified configurations of lens galaxies, leveraging *multiprocessing* and *njit* functionalities for enhanced efficiency.
    * The package utilizes the Elliptical Power Law with external shear (EPL+Shear) model for galaxy parameters $(\theta_L)$ sampling, following [Wierda et. al 2021](https://arxiv.org/abs/2106.06303). Rejection sampling is applied on the above samples on condition that whether the event is strongly lensed or not, $P(\theta_L|z_s,\text{SL})$.

- Generation of image properties:
    * Source position $(\beta)$ is sampled from the caustic in the source plane.
    * Sampled lens properties and source position is fed in *Lenstronomy* to generate properties of the images. This is the slowest part of the entire simulation, which <span style="font-family:cambria; color:maroon;">LeR</span> tackles through parallelization with multiprocessing.
    * Image properties like magnification $(\mu_i)$ and time delay ($\Delta t_i$) modifies the original source signal strength, changing the signal-to-noise ratio SNR and our ability to detect.

- Calculation of Detectable Merger Rates Per Year:
    * The calculation of rates necessitates integration over simulated events that meet specific detection criteria. This process includes computing SNRs $(\rho)$ for each event or its lensed images, followed by an assessment against a predetermined threshold(s) $(\rho_{th})$.
    * SNR calculations are optimized using [*gwsnr*](https://github.com/hemantaph/gwsnr), leveraging interpolation and multiprocessing for accuracy and speed.
    * Simulated events and rate results, along with input configurations, are systematically archived for easy access and future analysis. Additionally, all interpolators used in the process are preserved for future applications.

The <span style="font-family:cambria; color:maroon;">LeR</span> software has been developed to cater to the requirements of both the LIGO-Virgo-KAGRA Scientific Collaboration and research scholars engaged in astrophysics studies. It is currently used in generating detectable lensing events and GW lensing rates with the available information on current and future detectors. The results will predict the capacity of various detectors to detect and study such lensing events. Statistics generated from <span style="font-family:cambria; color:maroon;">LeR</span> will be used in event validation for the ongoing effort to detect lensed GWs. Lastly, <span style="font-family:cambria; color:maroon;">LeR</span> was designed with upgradability in mind to include additional statistics as required by the related research.

# Documentation

The `ler` package documentation is available at [ReadTheDocs](https://ler.readthedocs.io/en/latest/installation.html).




