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
  - name: Harsh Narola
    # equal-contrib: False
    affiliation: "3 , 2"
  - name: Otto A. HANNUKSELA 
    # equal-contrib: False 
    affiliation: "1 , 2"
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
date: 21 June 2023
bibliography: paper.bib
---

# Summary

Gravitational waves (GWs) are ripples in the fabric of space and time caused by the acceleration of unsymmetrically distributed mass/masses. Observable GWs are created especially during the violent events of merging compact binaries, such as 'binary black holes' (BBH) and 'binary neutron stars' (BNS). The gravitational waves emitted by these events are often distorted or magnified by the gravitational fields of massive objects such as galaxies or galaxy clusters, a phenomenon known as gravitational lensing. Understanding the effects of gravitational lensing on GW signals is crucial for accurately interpreting these signals and extracting astrophysical information from them. In this field of physics, statistical modelling of GWs lensing can provide valuable insights into the properties of the lensing objects and the sources of gravitational waves. Such statistics require accurate and efficient means to calculate the detectable lensing rates, which depend on up-to-date modelling and implementation of lens and source properties and its distribution. These computational results will not only help in producing reliable predictions but helps in event validation of future lensing events [cite](https://arxiv.org/abs/2306.03827). 

Producing accurate results in such statistics requires sampling in the order of millions and it's computationally expensive to do so.  LeR employs ingenious methods in parallelizing the workflow and also simplifying the calculation through interpolation and linear algebra. Besides, it was also a challenge to seamlessly integrate different parts of the statistics in a modular fashion that allows upgradability and extendability. LeR not only tackles these issues but also has easy-to-use features to work flawlessly with other related packages.  

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

`LeR` is a statistical-based python package whose core function is to calculate detectable rates of both lensing and unlensed GW events. This calculation is very much dependent on the other functionality of the package, which can be subdivided into three parts; 1. Sampling of compact-binary source properties, 2. Sampling of lens galaxy characteristics and 3. Solving the lens equation to get image properties of the source. The package as a whole relies on `numpy` array operation and linear algebra, `scipy` interpolation and `multiprocessing` functionality of python to increase speed and functionality without compromising on the ease of use. The API of `LeR` is structured such that each functionality mentioned stands in its own right for scientific research but also can be used together as needed. Key features of `LeR` and its dependencies can be summarized as follows,

- Detectable merger rates: 
    * Calculation not only relies on the properties of simulated events but also on detectability provided by the condition of the GW detectors. For this, the optimal signal-to-noise ratio (SNR) is calculated for each of the simulated events and it can be computationally expensive. This is mitigated because `LeR` relies on [`gwsnr`]{https://github.com/hemantaph/gwsnr/tree/main} for efficient and rapid calculation of SNRs. Due to the prowess of `gwsnr`, rate calculation can also be done both for present and future detectors with customizable sensitivities. 
    * Merger rates of both the simulated unlensed and lensed events can be calculated and compared. 
- Sampling GW sources:
    * Distribution of the source's red-shift is based on the merger rate density of compact binaries, which can be BBH, BNS, primordial black holes (PBHs) etc. The code is designed to accommodate easy updates or additions of such distribution by the users in the future. 
    * Sampling of BBH masses is done using `gwcosmo` following the powerlaw+peak model. Other related properties are sampled from available priors of `bilby`. The user can manually replace any before feeding the parameters in for rate computation.
- Sampling of lens galaxies:
    * Lens distribution follows [(Oguri et al. 2018](https://arxiv.org/abs/1807.02584). It depends on the sampled source red-shifts and also on the optical depth [cite](https://arxiv.org/abs/2306.03827).
    * `LeR` employs the Elliptical Power Law model with the external shear (EPL+Shear) model for sampling other galaxy features, which is available in the `Lenstronomy` package.
    * Rejection sampling is applied on the above samples on condition that whether the event is strongly lensed or not.
- Generation of image properties:
    * Source position is sampled from the caustic in the source plane.
    * Sampled lens properties and source position is fed in `Lenstronomy` to generate properties of the images.
    * Properties like magnification and time delay are essential as it modifies the source signal strength, changing the SNR and detection ability.
    * `LeR` can handle both super-threshold and sub-threshold events in picking detectable events and rate computation.

`LeR` was written to be used by both LIGO scientific collaboration and research students for related works in astrophysics. It is currently used in generating detectable lensing events and GW lensing rates with the available information on current and future detectors. The results will predict the feasibility of various detectors for detecting and studying such lensing events. Statistics generated from `LeR` will be used in event validation of the ongoing effort to detect lensed gravitational waves. Lastly, `LeR` was designed with upgradability in mind to include additional statistics as required by the related research. 

# Acknowledgements
I acknowledge support from the Department of Physics, The Chinese University of Hong Kong with the Postgraduate Studentship. 

# References