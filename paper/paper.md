---
title: '$ler$ : LVK (LIGO-Virgo-KAGRA collaboration) event (compact-binary mergers) rate calculator and simulator'
tags:
  - Python
  - astrophysics
  - statistics
  - gravitational waves
  - LIGO
authors:
  - name: Hemantakumar Phurailatpam
    orcid: 0000-0002-0471-3724
    # equal-contrib: True
    affiliation: "1"
  - name: Anupreeta More
    # equal-contrib: False
    orcid: 0000-0001-7714-7076
    # equal-contrib: False
    affiliation: "2,3"
  - name: Harsh Narola
    orcid: 0000-0001-9161-7919
    # equal-contrib: False
    affiliation: "4,5"
  - name: Ng Chung Yin (Leo)
    # equal-contrib: False
    affiliation: "1"
  - name: Justin Janquart
    orcid: 0000-0003-2888-7152
    # equal-contrib: False
    affiliation: "4,5,6,7"
  - name: Chris Van Den Broeck
    orcid: 0000-0001-6800-4006
    # equal-contrib: False
    affiliation: "4,5"
  - name: Otto Akseli Hannuksela
    orcid: 0000-0002-3887-7137
    # equal-contrib: False
    affiliation: "1"
  - name: Neha Singh
    orcid: 0000-0002-1135-3456
    # equal-contrib: False
    affiliation: "8"
  - name: David Keitel
    orcid: 0000-0002-2824-626X
    # equal-contrib: False
    affiliation: "8"
affiliations:
  - name: Department of Physics, The Chinese University of Hong Kong, Shatin, New Territories, Hong Kong
    index: 1
  - name: The Inter-University Centre for Astronomy and Astrophysics (IUCAA), Post Bag 4, Ganeshkhind, Pune 411007, India
    index: 2
  - name: Kavli Institute for the Physics and Mathematics of the Universe (IPMU), 5-1-5 Kashiwanoha, Kashiwa-shi, Chiba 277-8583, Japan
    index: 3
  - name: Department of Physics, Utrecht University, Princetonplein 1, 3584 CC Utrecht, The Netherlands
    index: 4
  - name: Nikhef – National Institute for Subatomic Physics, Science Park, 1098 XG Amsterdam, The Netherlands
    index: 5
  - name: Center for Cosmology, Particle Physics and Phenomenology - CP3, Universit\'e Catholique de Louvain, Louvain-La-Neuve, B-1348, Belgium
    index: 6
  - name: Royal Observatory of Belgium, Avenue Circulaire, 3, 1180 Uccle, Belgium
    index: 7
  - name: Departament de Física, Universitat de les Illes Balears, IAC3-IEEC, Crta. Valldemossa km 7.5, E-07122 Palma, Spain
    index: 8
date: 6 May 2026
bibliography: paper.bib
---

# Summary

Gravitational waves (GWs) are ripples in spacetime caused by accelerating massive objects, primarily produced by violent cosmic mergers of compact binaries such as binary black holes (BBHs), binary neutron stars (BNSs), and neutron star-black hole pairs (NSBHs). Since 2015, the LIGO-Virgo-KAGRA (LVK) collaboration [@LIGO:2015] [@VIRGO:2015] [@KAGRA:2021] has routinely detected such GW transients during observing runs [@GWTC4:2026] [@GWTC3:2023] [@GWTC2:2021] [@GWTC1:2019]. Most of these signals reach detectors unobstructed, forming a baseline unlensed population, but a small fraction can encounter massive intervening objects like galaxies or galaxy clusters. These objects act as lenses, splitting the signals into multiple images that are magnified, time-shifted, and phase-shifted through strong gravitational lensing [@Abbott:2021] [@Wierda:2021].

Accurate modeling of unlensed and lensed gravitational-wave populations is essential for astrophysical interpretation, inferring merger rate histories, and validating future lensing events [@Abbott:2021] [@Janquart:2023] [@ligolensing:2023] [@Wierda:2021]. The $ler$ software is a Python framework that provides a unified pipeline to model these populations statistically and calculate their detectable rates. Producing precise results requires computationally demanding tasks such as large-scale sampling of source and lens attributes, solving lens equations, and simulating detector responses. $ler$ overcomes these traditional bottlenecks through optimization and parallelization, reducing computational time significantly and achieving up to a thousandfold speedup [@lertestgithub:2026] over conventional methods. Source code, validation examples, and comprehensive documentation are available through the project repository [@lergithub:2026], supporting examples [@lertestgithub:2026], and documentation site [@lerdoc:2026], respectively.

# Statement of need

$ler$ is a statistics-based Python tool designed for the joint simulation and computation of detectable rates for both unlensed and lensed gravitational-wave events. It is intended for LVK analyses and for astrophysics researchers working on compact-binary populations and strong lensing. Its core functionality integrates the sampling of compact-binary source properties and lens-galaxy characteristics, solving lens equations to derive image properties, and computing detectable event rates. This functionality leverages *numpy* [@numpy:2022] for array operations and linear algebra alongside *scipy* [@scipy:2020] interpolation methods and Python multiprocessing. Efficiency is improved using *numba* [@numba:2022] for just-in-time compilation and multithreading, and by using inverse-transform sampling and importance sampling where rejection sampling is inefficient. $ler$ can use *lenstronomy* [@Birrer:2021] for lensing calculations and *gwsnr* [@Phurailatpam:2025] for efficient detection-probability calculations $P_{\rm det}$.

$ler$ produces parameter distributions for intrinsic populations and for populations selected by detection criteria, including strong-lensing selection and image-level detectability. This supports validation studies for lensing candidates [@Janquart:2023] and forecasting for present and future detector networks [@Ng:2025]. The design can be adapted to multi-messenger applications by changing the detection criteria, including studies of strongly lensed gamma-ray bursts [@More:2025]. $ler$ also supports integration with existing inference tools by providing $P_{\rm det}$ estimates that can be used to construct selection functions [@lertestgithub:2026] for hierarchical population inference frameworks [@Thrane:2019] with *bilby* software [@Ashton:2019].

# State of the field

Joint population studies of unlensed and strongly lensed gravitational waves require end-to-end simulations that combine source and lens-galaxy models with detector selection effects. Common inference packages such as *bilby* [@Ashton:2019] and cosmology and rate tools such as *gwcosmo* [@gwcosmo:2022] provide general-purpose likelihood evaluation, sampling, and population inference utilities. However, they do not provide a single, configuration-driven pipeline to generate unlensed and strongly lensed populations jointly at scale, nor to propagate lensing image properties into detectability. Lens-modelling tools such as lenstronomy [@Birrer:2021] implement fast analytical routines for lensing calculations but are not designed as high-throughput population-level simulators.

Previous strong-lensing population studies often relied on pipelines tightly coupled to specific lens models or detection criteria [@Haris:2018] [@Wierda:2021]. Many studies report rate estimates based on limited sample sizes without showing convergence diagnostics or uncertainty estimates. They also commonly neglect uncertainty in population-model hyperparameters, which should propagate into the detectable population and rate predictions. $ler$ addresses these gaps by providing an integrated, high-throughput workflow linking population draws, lensing cross-sections, image properties, and detector selection in one reproducible pipeline.

# Software design

$ler$ is designed as an end-to-end population simulator with a clear separation between the scientific steps and the computational machinery. The pipeline is organized around two rate calculations. Unlensed rates are computed by drawing compact-binary populations from astrophysical priors, evaluating detectability, and integrating over the detected events. Lensed rates extend this flow by sampling lens-galaxy and source-redshift parameters under a strong-lensing condition, solving the lens equation to generate image properties, applying detectability criteria at the image level, and integrating over the detectable lensed events. This separation keeps the rate logic explicit while allowing the source model, lens model, and detection model to be replaced without rewriting the full pipeline. The mathematical foundations for these implementations are detailed in the Mathematics section.

A core design choice is modularity with stable interfaces between modules that handle source populations, lens populations, image-property calculations, and rate estimation. Users can select built-in models such as BBH, BNS, or NSBH populations and lens models such as SIS (singular isothermal sphere), SIE (singular isothermal ellipsoid), or EPL+Shear (elliptical power-law with external shear), or provide their own distributions and detection criteria. Detector selection is handled through a backend interface so that signal-to-noise and detection probability calculations can be delegated to *gwsnr* or replaced by alternative implementations. This architecture supports extension and testing while keeping the scientific assumptions visible in configuration rather than embedded in ad hoc scripts.

Performance is obtained by compiling the repeatedly evaluated parts of the pipeline with *Numba* `njit` and parallelizing the main Monte Carlo loops. $ler$ uses *lenstronomy* [@Birrer:2021] as an optional backend, but it provides an in-house analytical EPL+Shear solver and caustic and cross-section calculator that are compiled with *Numba* and run in parallel, following the same mathematical formulation as *lenstronomy* and Tessore (2015) [@Tessore:2015]. This leads to two practical trade-offs. First, the rate estimates are Monte Carlo approximations, so accuracy is controlled by sample size rather than closed-form expressions, and $ler$ reports convergence statistics from repeated batches. Second, just-in-time compilation adds overhead on first use, which is about 2 s for the first unlensed batch and about 25 s for the first lensed batch that includes image-property calculations. After compilation, the compiled kernels reduce per-sample cost and make large simulations feasible. In benchmarks that disable just-in-time compilation and parallelism, these choices provide speed gains of about $180\times$ for unlensed populations and $5189\times$ for lensed populations. On a standard six-core processor, a batch of $10^5$ samples takes about 300 ms for unlensed calculations and about 30 s for lensed calculations, and repeated batches converge to a stable rate estimate.

# Mathematics

The mathematical workflow in $ler$ can be summarised as Monte Carlo estimates of detectable event rates for unlensed and strongly lensed populations, obtained by averaging detection probabilities over simulated events drawn from the relevant population models. For unlensed events, the detector-frame rate of detectable events is

$$
\frac{\Delta N^{\rm obs}_{\rm U}}{\Delta t} = {\cal N}_{\rm U} \bigg\langle P({\rm obs} \mid \vec{\theta}) \bigg\rangle_{\vec{\theta} \sim P(\vec{\theta})}
$$

where ${\cal N}_{\rm U}$ is the total intrinsic merger rate in the detector frame, $\vec{\theta}$ denotes the source parameters drawn from $P(\vec{\theta})$, and $P({\rm obs}\mid \vec{\theta})$ is the detection probability.

For strongly lensed events, the detector-frame rate of detectable lensed events is

$$
\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t} = {\cal N}_{\rm L} \bigg\langle P({\rm obs}\mid \vec{\theta}_{\rm U}, \vec{\theta}_{\rm L}, \vec{\beta}, {\rm SL}) \bigg\rangle_{\substack{ \vec{\theta}_{\rm U},\vec{\theta}_{\rm L} \sim P(\vec{\theta}_{\rm U},\vec{\theta}_{\rm L} \mid z_L, z_s, {\rm SL}) \\ \vec{\beta} \sim P(\vec{\beta} \mid z_s, \vec{\theta}_{\rm L}, {\rm SL}) }} \, ,
$$

where ${\cal N}_{\rm L}$ is the total intrinsic merger rate in the detector frame for the lensed population, $\vec{\theta}_{\rm U}$ and $\vec{\theta}_{\rm L}$ denote the source and lens parameters, and $\vec{\beta}$ is the source position in the source plane. The distributions $P(\vec{\theta}_{\rm U},\vec{\theta}_{\rm L}\mid z_L,z_s,{\rm SL})$ and $P(\vec{\beta}\mid z_s,\vec{\theta}_{\rm L},{\rm SL})$ define the hierarchical sampling under the strong-lensing condition ${\rm SL}$, and $P({\rm obs}\mid \vec{\theta}_{\rm U},\vec{\theta}_{\rm L},\vec{\beta},{\rm SL})$ is evaluated using the image properties implied by $(\vec{\theta}_{\rm L},\vec{\beta})$ and the chosen detectability criteria.

# Research impact statement

$ler$ has established a measurable impact within the gravitational-wave community through its integration into research workflows and its growing user base. The package has recorded more than 3000 downloads from PyPI [@lerpypi:2026] and supports a diverse community ranging from students to faculty members at multiple institutions. Development is driven by contributions from researchers across various organizations and has undergone peer review within the LIGO-Virgo-KAGRA Scientific Collaboration [@lerlvkpnpreview:2024]. Active community engagement is maintained through bug reports and feature requests on the GitHub issue tracker and through collaboration communication channels such as LIGO Mattermost.

The software has facilitated several scientific publications related to gravitational-wave lensing statistics. It has been used to analyze the distribution of time delays and magnification ratios for the validation of lensing candidates [@Janquart:2023] and to forecast the detection of sub-threshold lensed events with future observatories [@Ng:2025]. Additional research applications include the generation of lensed-event statistics and the calculation of posterior odds [@Hannuksela:2026], as well as forecasting electromagnetic counterparts for strongly lensed gamma-ray bursts [@More:2025]. These implementations demonstrate that the tool provides a robust framework for high-level population studies and astrophysical forecasting.


# AI usage disclosure

No generative AI tools were used to develop the software, write this manuscript, or prepare the accompanying materials.


# Acknowledgements

The authors thank their academic advisors for guidance and support. Hemantakumar Phurailatpam acknowledges the Department of Physics at The Chinese University of Hong Kong for the Postgraduate Studentship that facilitated this research. Hemantakumar Phurailatpam and Otto A. Hannuksela acknowledge support from the Research Grants Council of Hong Kong (Project Nos. CUHK 14304622 and 14307923), the start-up grant from The Chinese University of Hong Kong, and the Direct Grant for Research from the Research Committee of The Chinese University of Hong Kong. The authors also acknowledge support from the Netherlands Organisation for Scientific Research (NWO). N. Singh and D. Keitel are supported by Universitat de les Illes Balears (UIB); the Spanish Agencia Estatal de Investigación grants CNS2022-135440, PID2022-138626NB-I00, RED2022-134204-E, RED2022-134411-T, funded by MICIU/AEI/10.13039/501100011033, the European Union NextGenerationEU/PRTR, and the ERDF/EU; and the Comunitat Autònoma de les Illes Balears through the Direcció General de Recerca, Innovació I Transformació Digital with funds from the Tourist Stay Tax Law (PDR2020/11 - ITS2017-006) as well as through the Conselleria d’Economia, Hisenda i Innovació with grant numbers SINCO2022/6719 (European Union NextGenerationEU/PRTR-C17.I1) and SINCO2022/18146 (co-financed by the European Union and FEDER Operational Program 2021-2027 of the Balearic Islands). The authors acknowledge the computational resources provided by the LIGO Laboratory, supported by National Science Foundation Grants No. PHY-0757058 and No. PHY-0823459.

# References
