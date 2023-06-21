---
title: 'LeR: A Python package for generating gravitational waves lensing statistics'
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
    # equal-contrib: true 
    affiliation: "1 , 2"
affiliations:
 - name: The Chinese University of Hong Kong, Hong Kong
   index: 1
 - name: LIGO scientific collaboration
   index: 2
date: 21 June 2023
bibliography: paper.bib
---

# Summary

Gravitational waves (GWs) are ripples in the fabric of space and time caused by acceleration of unsymmetrically distributed mass/masses. Observable GWs are created especially during the violent events of merging compact binaries, such as 'binary black-holes' (BBH) and 'binary neutron stars' (BNS). The gravitational waves emitted by these events are often distorted or magnified by the gravitational fields of massive objects such as galaxies or galaxy clusters, a phenomenon known as gravitational lensing. Understanding the effects of gravitational lensing on GW signals is crucial for accurately interpreting these signals and extracting astrophysical information from them. In this field of physics, statistical modeling of GWs lensing can provide valuable insights into the properties of the lensing objects and the sources of gravitational waves. Such statistics requires accurate and efficient means to calculate the detectable lensing rates which in turn depends on upto-date modeling and implementation of lens and source properties and its distribution. These computational results will not only help in producing reliable predictions but helps in event validation of future lensing events [cite](https://arxiv.org/abs/2306.03827).

# Statement of need

`LeR` is a statistical based python package whose core function is to calculate detectable rates of both lensing and unlensed GW events. This calculation very much dependent on the other functionality of the package, which can be subdivided into three parts; 1. Sampling of compact binary source properties, 2. Sampling of lens galaxy characteristics and 3. Solving the lens equation to get image properties of the source. The package as a whole relies on `numpy` array operation, `scipy` interpolation and `multiprocessing` functionality of python to increase speed and functionality without compromising on the ease-of-use. The API of `LeR` is structure such that each functionality mentioned stands on this own right for scientific research but also can be used together as needed. Keys features of `LeR` and its dependencies can be summarized as follows,

- Detectable merger rates: 
    * Calculation not only relies on the properties of simulated events but also on detectability provided by the condition of the GW detectors. For this, `LeR` relies on `gwsnr` for the calculation of optimal signl-to-noise ratio (SNR). Due to prowess of `gwsnr`, rate calulation can be done both for present and future detectors with customizable sensitivities. 
    * Merger rates of both the simulated unlensed and lensed events can be calculated and compared. 
- Sampling GW sources:
    * Distribution of source's red-shift is based on the merger rate density of compact binaries, which can be BBH [cite](https://arxiv.org/abs/2306.03827), BNS [cite](https://arxiv.org/abs/2306.03827), primodial black holes (PBHs) [cite](https://arxiv.org/abs/2306.03827) etc. The code is designed to accomodate easy updates or additions of such distribution by the users in the future. 
    * Sampling of BBH masses is done using `gwcosmo` follwing the powerlaw+peak model. Other related properties are sampled from available priors of `bilby`. Each of them can me manually replaced by the user before feeding in for rate computation.
- Sampling of lens galaxies:
    * Lens distribution of follows [cite](https://arxiv.org/abs/2306.03827). It depends on the sampled source red-shifts and also on the optical depth [cite](https://arxiv.org/abs/2306.03827).
    * `LeR` employs Elliptical Power Law model with external shear (EPL+Shear) model for sampling other features of the galaxy, which is available in the `Lenstronomy` package.
    * Rejection sampling is applied on the above samples on condition that whether event is strongly lensed or not.
- Generation of image properties:
    * Source position is sampled from the caustic in the source plane.
    * Sampled lens' properties and source position is fed in `Lenstronomy` to generate properties of the images.
    * Properties like magnification and time-delay is important as it modifies the source signal strength which in turns changes the SNR and detect ability.
    * `LeR` can handle both super-therhold and sub-threshold events in picking detectable events and rate computation.

`LeR` was written to used by both LIGO scientific collaboration and research students for related works in astrophysics. It is currently use in generating detectable lensing events and GW lensing rates with the available information on current and future detectors. The results will predicts the feasibility of various detectors on the detection of such lensing events. Statistics generated from `LeR` will be use in event validation of the ongoing effort to detected lensed gravitational waves. Lastly, `LeR` was design with upgradability in mind to include additional statistics as required by the related research. 

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

# Acknowledgements

We acknowledge NG Chung Yin (Leo) for bug reports.

# References