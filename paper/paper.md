---
title: 'LeR: A Python package for generating gravitational waves lensing rates and related statistics'
tags:
  - Python
  - astrophysics
  - statistics
  - gravitational waves
  - LIGO
authors:
  - name: Phurailatpam Hemantakumar
    orcid: 0000-0000-0000-0000
    # equal-contrib: true
    affiliation: "1 , 2"
  - name: Otto A. HANNUKSELA 
    # equal-contrib: False 
    affiliation: "1 , 2"
  - name: Harsh Narola
    # equal-contrib: False
    affiliation: "3 , 2"
  - name: Anupreeta More
    # equal-contrib: False
    affiliation: "5 , 6"
  - name: Ng Chung Yin (Leo)
    # equal-contrib: False 
    affiliation: "1 , 2"
  - name: Justin Janquart
    # equal-contrib: False 
    affiliation: "3 , 2, 4"
  - name: Chris Van Den Broeck
    # equal-contrib: False 
    affiliation: "3 , 2, 4"
affiliations:
  - name: The Chinese University of Hong Kong, Hong Kong
    index: 1
  - name: LIGO scientific collaboration
    index: 2
  - name: Department of Physics, Utrecht University, Princetonplein 1, 3584 CC Utrecht, The Netherlands
    index: 3
  - name: Nikhef – National Institute for Subatomic Physics, Science Park, 1098 XG Amsterdam, The Netherlands
    index: 4
  - name:  The Inter-University Centre for Astronomy and Astrophysics (IUCAA), Post Bag 4, Ganeshkhind, Pune 411007, India
    index: 5
  - name: Kavli Institute for the Physics and Mathematics of the Universe (IPMU), 5-1-5 Kashiwanoha, Kashiwa-shi, Chiba 277-8583, Japan
    index: 6
date: 17 August 2023
bibliography: paper.bib
---

# Summary

Gravitational waves (GWs) are ripples in the fabric of space and time caused by the acceleration of unsymmetrically distributed mass/masses. Observable GWs are created especially during the violent cosmic events of merging compact binaries, such as 'binary black holes' (BBH) and 'binary neutron stars' (BNS). The gravitational waves emitted by these events are often distorted or magnified by the gravitational fields of massive objects such as galaxies or galaxy clusters, a phenomenon known as gravitational lensing. Profound comprehension of gravitational lensing's impact on GW signals is imperative to the accurate interpretation of these signals and the extraction of astrophysical insights therein. In this field of physics, statistical modelling of GWs lensing can provide valuable insights into the properties of the lensing objects and the sources of gravitational waves. Such statistics require accurate and efficient means to calculate the detectable lensing rates, which depend on up-to-date modelling and implementation of lens and source properties and its distribution. The outcomes of these computational analyses not only contribute to generating dependable forecasts but also play a pivotal role in validating forthcoming lensing events [ref](https://arxiv.org/abs/2306.03827). 

Obtaining precise outcomes in statistical analyses of this nature necessitates the utilization of large-scale sampling, often numbering in the millions. However, this process is computationally demanding. The `LeR` framework employs innovative techniques to parallelize the workflow efficiently. Additionally, it simplifies calculations through the utilization of interpolation and linear algebra. A further challenge entailed the seamless integration of disparate statistical components in a modular manner that facilitates adaptability, upgradability, and extendability. The LeR framework effectively addresses these complexities, offering user-friendly features that ensure seamless compatibility with other associated software packages.

# Equations

$\textbf{Detectable Unlensed rates:}$

\begin{equation*}
\begin{split}
R_U = \int & dz_s R_m^U(z_s)\left\{\Theta[\rho(z_s,\theta)-\rho_{th}] P(\theta) d\theta \right\}
\end{split}
\end{equation*}

* $z_s$: source red-shift, $R_m^U(z_s)$: source frame merger rate density in the co-moving volume at $z_s$, $\theta$: source parameters, $P$: probability distribution, $\rho$: SNR, $\rho_{th}$: SNR threshold, $\Theta$: function to select detectable events.

$\textbf{Detectable Lensed rates:}$

\begin{equation*}
\begin{split}
R_L = \int & dz_s R_m^L(z_s) \,\mathcal{O}_{images}(z_s,\theta,\mu_i,\Delta t_i, \rho_{th}) \, \\ 
& \, P(\theta) P(\theta_L, z_L|\text{SL},z_s) P(\beta|\text{SL}) d\theta d\beta dz_L d\theta_L dz_s 
\end{split}
\end{equation*}

* $R_m^L(z_s)$: strongly lensed (optical depth applied) source frame merger rate density in the co-moving volume at $z_s$, $\theta_L$: lens parameters, $\beta$: image properties, $\mu$: image magnification, $\Delta t$: image time delay, $\mathcal{O}$: function to select detectable lensed events, $\text{SL}$: strong lensing condition.

# Statement of need

`LeR` is a statistical-based python package whose core function is designed for the computation of detectable rates pertaining to both lensed and unlensed gravitational wave (GW) events. This calculation intricately hinges upon the interplay of various components within the package, which can be categorized into three primary segments: 1. Sampling the properties of compact-binary sources, 2. Sampling characteristics of lens galaxies, and 3. Solving the lens equation to derive image attributes of the source. The holistic functionality of the package is built upon leveraging array operations and linear algebra from the numpy library, interpolation from scipy, and the `multiprocessing` capability inherent to Python. This design optimizes both speed and functionality while upholding user-friendliness. The architecture of the "LeR" API is deliberately organized such that each distinct functionality holds its own significance in scientific research. Simultaneously, these functionalities seamlessly integrate and can be employed collectively based on specific research requirements. Key features of `LeR` and its dependencies can be summarized as follows,

- Detectable merger rates: 
    * The calculation depends not only on simulated event properties but also on GW detector detectability. We compute optimal signal-to-noise ratios (SNRs) for simulated events, which can be computationally intensive. `LeR` mitigates this using [`gwsnr`](https://github.com/hemantaph/gwsnr) for efficient  and rapid SNR calculation. `gwsnr` enables rate calculation for current and future detectors with customizable sensitivities.
    * The merger rates of both the simulated unlensed and lensed events can be computed and subsequently compared. 
- Sampling GW sources:
    * The distribution of the source's red-shift is based on the merger rate density of compact binaries, including binary black hole (BBH), binary neutron star (BNS) etc. The code is meticulously structured to facilitate straightforward incorporation of future updates or additional distributions of such sources by users.
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

# Acknowledgements
The author acknowledges the essential contributions that have facilitated the completion of this research. Gratitude is extended to the academic advisors for their invaluable guidance and unwavering support. The author appreciates the collaborative efforts and insightful discussions with fellow research colleagues, enriching the quality of the study. Special recognition is directed towards the Department of Physics, The Chinese University of Hong Kong, for providing the Postgraduate Studentship that supported this research endeavor. The author also acknowledges the participation of individuals who contributed to the empirical foundation of this work. 

# References