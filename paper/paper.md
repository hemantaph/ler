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
  - name:  The Inter-University Centre for Astronomy and Astrophysics (IUCAA), Post Bag 4, Ganeshkhind, Pune 411007, India
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
date: 19 June 2024
bibliography: paper.bib
---

# Summary

<!-- Summary should include: A description of the high-level functionality and purpose of the software for a diverse, non-specialist audience. -->

<!-- domain knowledge 1 -->
<!-- Gravitational waves (GWs) are ripples in the fabric of space and time caused by the acceleration of unevenly distributed mass or masses. Observable GWs are created especially during the violent cosmic events of merging compact binaries, such as 'binary black holes' (BBHs), 'binary neutron stars' (BNSs), and 'neutron star and black hole pair' (NSBHs). -->
<!-- domain knowledge 2 -->
<!-- As these signals travel across the universe, the vast majority reach our detectors unobstructed, forming the baseline (unlensed) GW population.
However, a fraction (<0.1%) of these signals encounter massive objects (lenses) like galaxies or galaxy clusters along their path, which can split them into multiple images, that are magnified (de-magnified), time and phase shifted, through strong gravitational lensing [Wierda2021] [@Abbott2021]. -->
<!-- narrower problem. What's the current state. Purpose of the work -->
<!-- A profound comprehension and accurate modelling of both the unlensed population and lens population and the lensing effects is imperative for accurate astrophysical interpretation, inferring merger rate histories, and validating forthcoming lensing events [@Janquart2023] [@ligolensing2023] [@Abbott2021] [@Wierda2021]. -->

Gravitational waves are spacetime ripples caused by accelerating massive objects, primarily observable from the violent cosmic mergers of compact binaries like 'binary black holes' (BBHs), 'binary neutron stars' (BNSs), and 'neutron star and black hole pairs' (NSBHs). While most of these signals reach our detectors unobstructed to form the baseline unlensed population, a small fraction encounter massive intervening objects such as galaxies or galaxy clusters. These objects act as lenses, splitting the signals into multiple images that are magnified, time-shifted, and phase-shifted through strong gravitational lensing [@Abbott2021] [@Wierda2021].

<!-- what does software do? high level functionality. -->
<!-- $ler$ is a Python framework designed for the comprehensive statistical modeling of GW populations, providing a unified pipeline to calculate detectable rates for both unlensed and lensed events. Obtaining precise outcomes in these statistical analyses necessitates large-scale sampling of source properties, lens galaxy attributes, and lensing-induced image properties calculation (lens equation solves), along with the computation of detectable rates based on simulated detector responses, often numbering in the millions. This process is traditionally highly demanding. The $ler$ framework addresses this by employing ... techniques to optimize the workflow and computation efficiency required for handling large-scale statistical analyses, which leads to significant reductions in computational time (1000x speed up than traditional methods). Detailed description, source code, and examples are available in $ler$ [documentation](https://ler.hemantaph.com). -->

Accurate modeling of unlensed and lensed gravitational wave populations is essential for astrophysical interpretation, inferring merger rate histories, and validating future lensing events [@Abbott2021] [@Janquart2023] [@ligolensing2023] [@Wierda2021]. The $ler$ software is a Python framework that provides a unified pipeline to model these populations statistically and calculate their detectable rates. Generating precise outcomes requires computationally demanding tasks, such as large-scale sampling of source and lens attributes, solving lens equations, and simulating detector responses. $ler$ overcomes these traditional bottlenecks by utilizing advanced interpolation and parallelization techniques. This optimized workflow reduces computational time significantly, achieving a thousandfold speedup over conventional methods. Source code, examples, and comprehensive details are available in the $ler$ [documentation](https://ler.hemantaph.com).

# Statement of need

<!-- Statement of need should include: A section that clearly illustrates the research purpose of the software and places it in the context of related work. This should clearly state what problems the software is designed to solve, who the target audience is, and its relation to other work. -->

<!-- what is ler? Why is it needed? What does it do? Who is it for? How does it relate to other work? -->
<!-- $ler$ is a statistics-based Python package specifically designed for efficient simulation and computing detectable rates of both lensed and unlensed GW events, catering to the requirements of the LIGO-Virgo-KAGRA Scientific Collaboration [@LIGO2015] [@VIRGO2015] [@KAGRA2021] and astrophysics research scholars. -->

<!-- The core functionality of $ler$ seemless integration of various components which include sampling the properties of compact-binary sources, lens galaxies characteristics, solving lens equations to derive properties of resultant images, and computing detectable GW rates. This comprehensive functionality builds on the leveraging of array operations and linear algebra from the *numpy* [@numpy] library, enhanced by interpolation methods from *scipy* [@scipy] and Python’s *multiprocessing* capabilities. Efficiency is further boosted by the *numba* [@numba] library's Just-In-Time (*njit*) compilation and *prange* multithreading, optimizing extensive numerical computations and employing the inverse transform sampling method and importance sampling to replace more cumbersome rejection sampling. $ler$ lenverages the *lenstronomy* [@Birrer2021] package for lensing calculations and the *gwsnr* package for efficient probability of detection (Pdet), ensuring accurate and rapid computations of lensing effects and detectability assessments. -->

$ler$ is a statistics-based Python tool designed for the efficient simulation and computation of detectable rates for both unlensed and lensed gravitational wave events, catering directly to the LIGO-Virgo-KAGRA Scientific Collaboration [@KAGRA2021; @LIGO2015; @VIRGO2015] and astrophysics researchers. Its core functionality seamlessly integrates the sampling of compact-binary source properties and lens galaxy characteristics, solving lens equations to derive image properties, and computing detectable event rates. This functionality leverages *numpy* [@numpy] for array operations and linear algebra alongside *scipy* [@scipy] interpolation methods and Python multiprocessing. Efficiency is further enhanced by *numba* [@numba] just-in-time compilation and multithreading, employing inverse transform sampling and importance sampling over cumbersome rejection techniques. Furthermore, $ler$ integrates *lenstronomy* [@Birrer2021] for lensing calculations and *gwsnr* [@Phurailatpam2025] for efficient probability of detection ($P_{\rm det}$) computations.

<!-- Moreover, the modular design of $ler$ not only optimizes speed and functionality but also ensures adaptability and upgradability, supporting the integration of additional statistics as research evolves. Currently, $ler$ is an important tool in generating simulated GW events—both lensed and unlensed—and provides astrophysically accurate distributions of event-related parameters for both detectable and non-detectable events. This functionality aids in event validation and enhances the forecasting of detection capabilities across various GW detectors to study such events [paper1] [paper2], which also include inference of EM-counterpart inference (with just defining a new detection criteria for GRB detectors) [paper3]. The architecture of the $ler$ API can facilitates seamless compatibility with other software packages, allowing researchers to integrate and utilize its functionalities based on specific scientific requirements, e.g. rapid population generation with Pdet to study selection effects, which can integarte to hierarchical inference of population level properties with *bilby* [@Ashton2019]. -->

Beyond computational speed, the modular architecture of $ler$ ensures adaptability and upgradability, supporting the integration of additional statistics as research evolves. It generates astrophysically motivated parameter distributions for both detectable and non-detectable simulated populations. This capability is important for validating lensing candidates from the observed events, and forecasting future detections across various detectors. The flexible design also extends to multi-messenger astronomy, allowing users to adapt detection criteria to study electromagnetic counterparts like strongly lensed gamma-ray bursts [@more2025]. Finally, the $ler$ API allows for straightforward integration with other widely used ecosystem tools. For example, its rapid population generation with $P_{\rm det}$ can be incorporated directly into packages like bilby [@Ashton2019] to account for selection effects during the hierarchical inference of population-level properties.

# State of the field

<!-- State of the field should include: A description of how this software compares to other commonly-used packages in the research area. If related tools exist, provide a clear “build vs. contribute” justification explaining your unique scholarly contribution and why existing alternatives are insufficient. -->

<!-- Combine population studies of GWs and strongly lensed GWs require end-to-end simulations that combine astrophysical priors for compact-binary sources and lens galaxies with image generation (magnifications and time delays) and detector selection. While general-purpose GW inference tools provide waveform likelihood evaluation and generic sampling utilities, they do not provide an integrated, configuration-driven workflow for generating *strongly lensed* GW populations at scale, propagating lensing-induced image properties into detectability, and producing rate predictions and associated parameter distributions suitable for quick forecasting and validation [bilby] [gwcosmo] [lenstronomy]. Moreover, statistics and framework implemented in the previous literature [paper1] [paper2] [paper3] for population studies and forcast are often specific to particular lens models, source populations, or detection criteria, and are not easily adaptable to updated astrophysical models and updated lensing statistics (e.g. implementation lens model specific optical depth calculation, source-lens based configuration based cross-section calculation and cross-section based strong lensing selection with importance sampling). Furthermore, rate computation requires numerical integration over large simulated populations, where number of events usually 1 million in previous literature, without clear reason for the choice of number of events, and it doesn't gaurantee the convergence of the results. -->

Joint population studies of standard unlensed and strongly lensed gravitational waves require end-to-end simulations that combine astrophysical priors for compact-binary sources and lens galaxies with image generation and detector selection metrics. While general-purpose GW inference tools provide waveform likelihood evaluation and generic sampling utilities, they do not provide an integrated, configuration-driven workflow for generating strongly lensed GW populations at scale, propagating lensing-induced image properties into detectability, and producing rate predictions and associated parameter distributions suitable for quick forecasting and validation [@Ashton2019] [@gwcosmo] [@Birrer2021]. Furthermore, statistical frameworks in previous literature are often rigid and tailored to specific lens models, source populations, or detection criteria [@Haris2018]. This makes them difficult to adapt for updated astrophysical models or advanced statistics like cross-section-based strong lensing selection. Previous studies also typically limit numerical rate integration to around one million events without clear mathematical justification for convergence.

<!-- $ler$ addresses this gap. It is designed, following statistics from previous work [paper] [paper2] (but improved upon), to fill this gap by providing a modular user friendly flexible pipeline specialized for population statistics (jointly for lensed and unlensed events) in the GW context. With it's cached interpolators, just in time compilation and parallelization, allows quick simulation and rate comuration for large populations (e.g., 10 million events) and convergence of rate integrals in a reasonable time frame (5 min for unlensed and 1 hour for lensed events on 6 CPU cores), which is crucial for robust forecasting studies and building statistics for validation of lensing events. The modular design of $ler$ not only optimizes speed and functionality but also ensures adaptability and upgradability, supporting the integration of additional statistics as research evolves. -->

$ler$ addresses these gaps by providing a flexible pipeline specialized for the joint population statistics of lensed and unlensed events with the upto date statistics and and population models. By utilizing cached interpolators, just-in-time compilation, and parallelization, $ler$ enables the rapid simulation and rate computation of massive populations. It successfully processes ten million events and achieves rate integral convergence efficiently, requiring approximately five minutes for unlensed events and one hour for lensed events on a standard six-core processor. This performance is vital for robust forecasting and building the extensive statistics required to validate actual lensing events.

# Software design

<!-- Software design should include: An explanation of the trade-offs you weighed, the design/architecture you chose, and why it matters for your research application. This should demonstrate meaningful design thinking beyond a superficial code structure description. -->
Design philosophy of $ler$ is centered around modularity, efficiency, and adaptability. To understand this, we can break down $ler$ functionality into 2 main components (also refer to Mathematics section): 
(1) Unlnesd rates: requiring integration over simulated events that meet specific detection criteria. simulated events depends on module that handles sampling of GW intrinsic and extrinsic parameters.
(2) Lensed rates: requires integration over simulated events that meet strongl lensing and specific detection criteria for lensed events. simulated events depends on module that handles sampling of lens galaxy attributes and source red-shifts under strong lensing condition, generation of image properties through lens equation solving, and selection of detectable lensed events based on image properties and detection criteria. 

<!-- ler philosophy in action -->
<!-- modularity -->
The architecture of the $ler$ API is deliberately organized such that each distinct functionality holds its own significance in scientific research. Simultaneously, these functionalities seamlessly integrate and can be employed collectively to accommodate diverse scientific objectives. 

Modular pipeline specialized for population statistics in the GW context, allowing easy plug and play with user provided settings (e.g., source related: 'BBH', 'BNS', 'NSBH', lens related: 'EPL+Shear', 'SIE', 'SIS') and models (e.g., source related: 'merger rate density', 'binary source masses' etc., lens related: 'velocity dispersion', 'axis-atio', etc.)

Modularity is a key design principle, where each component is developed as an independent module, allowing for easy maintenance, testing, and future enhancements. This includes modules that deals with GW source properties, lens galaxy attributes, image properties, and rate calculations. SNR computation and Pdet calculation is handle seperately by `gwsnr` backend. 

<!-- efficiency -->
Most functions in $ler$ is a class object, which encompases cosmology related (e.g., luminosity distance, comoving volume [@astropy]), sampling related (e.g., sampling of source and lens properties), selection effect related (e.g., cross-section calculation). This class object leverages njitted interpolation of any smooth fuctions (single or multi-dimensional) or PDFs or CDFs (for inverse transform sampling) to optimize the workflow. Furthermore, for sampling strongly lensed events (with cross-section based strong lensing selection with importance sampling), parallelization is implemented through *prange* to further speed up the process. The major bottleneck of the entire workflow is lens equation solving, which is handled by *lenstronomy* [@Birrer2021], but $ler$ tackles this through parallelization with multiprocessing. SNR calculations are optimized using [*gwsnr*](https://gwsnr.readthedocs.io/en/latest/) python package backend, which again uses njitted interpolation and multiprocessing. This combination of techniques allows $ler$ to achieve significant speed-ups (1000x) compared to traditional methods, enabling the generation of large simulated populations (e.g., 10 million events) and convergence of rate integrals in a reasonable time frame (5 min for unlensed and 1 hour for lensed events on 6 CPU cores), which is crucial for robust forecasting studies and building statistics for validation of lensing events.
<!-- adaptability -->

The modular design also facilitates the integration of additional statistics and models as research evolves, ensuring that $ler$ remains adaptable to new scientific insights and requirements. $ler$ inherently supports customization, allowing users to input their own distributions for source and lens properties, as well as custom detection criteria, making it a versatile tool for a wide range of astrophysical studies.

# Mathematics

The entire workflow of $ler$ hinges on the computation of detectable rates, which involves numerical integration (monte-carlo method) over simulated events that meet specific detection criteria; so you can get distributions of intrinsic populations, selected population with strong lensing condition, and selected population with just unlensed event detection criteria or strong lensing plus lensed event detection criteria. The simulation of events depends on sampling of source and lens properties, generation of image properties through lens equation solving, and selection of detectable events based on image properties and detection criteria.

The mathematical formulation of these rates is as follows:

The annual rate of detectable/observable (unlensed) GW events $\frac{\Delta N^{\rm obs}_{\rm U}}{\Delta t}$, reads,

$$
\begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm U}}{\Delta t}
&= {\cal N}_{\rm U} \bigg\langle P({\rm obs} \mid  \vec{\theta}) \bigg\rangle_{\vec{\theta} \sim P(\vec{\theta})} \,,
\end{split}
\end{equation}
$$

where ${\cal N}_{\rm U}$ is the total intrinsic merger rate in the detector-frame, $\vec{\theta}$ is the vector of GW intrinsic and extrinsic parameters, $P(\vec{\theta})$ is the joint distribution of source parameters from where sampling is done, and $P({\rm obs} \mid  \vec{\theta})$ is the probability of detection for a source with parameters $\vec{\theta}$.

The annual rate of detectable/observable (lensed) GW events $\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}$, reads,

$$
\begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}
&=
{\cal N}_{\rm L}
\bigg\langle
P({\rm obs}\mid \vec{\theta}_{\rm U}, \vec{\theta}_{\rm L}, \vec{\beta}, {\rm SL})
\bigg\rangle_{\substack{
\vec{\theta}_{\rm U},\vec{\theta}_{\rm L} \sim P(\vec{\theta}_{\rm U},\vec{\theta}_{\rm L} \mid z_L, z_s, {\rm SL}) \\
\vec{\beta} \sim P(\vec{\beta} \mid z_s, \vec{\theta}_{\rm L}, {\rm SL})
}} \, ,
\end{split}
\end{equation}
$$

where ${\cal N}_{\rm L}$ is the total intrinsic merger rate in the detector-frame of the lensed sources, $\vec{\theta}_{\rm L}$ are the vectors lens galaxy related parameters, $P(\vec{\theta}_{\rm U},\vec{\theta}_{\rm L} \mid z_L, z_s, {\rm SL})$ is the joint distribution of GW and lens parameters under strong lensing condition from where sampling is done hierarchichally based on optical depth and cross-section, $\vec{\beta}$ is the vector of source position in the source plane, $P(\vec{\beta} \mid z_s, \vec{\theta}_{\rm L}, {\rm SL})$ is the distribution of source position in the source plane under strong lensing condition from where sampling is done and this results in image properties, and $P({\rm obs}\mid \vec{\theta}_{\rm U}, \vec{\theta}_{\rm L}, \vec{\beta}, {\rm SL})$ is the probability of detection for a lensed source with parameters $\vec{\theta}_{\rm U}$, $\vec{\theta}_{\rm L}$, $\vec{\beta}$ under strong lensing condition.

# Research impact statement

<!-- Research impact statement should include: Evidence of realized impact (publications, external use, integrations) or credible near-term significance (benchmarks, reproducible materials, community-readiness signals). The evidence should be compelling and specific, not aspirational. -->

- substantial research impact
- user base and community engagement
- the package involves contributions from multiple researchers across different institutions
- peer reviewed within the LIGO-Virgo-KAGRA Scientific Collaboration
- reporting bugs in ligo mattermost channel and github issues

- the package has been used in multiple research papers and cited, including 
  - [@Janquart2023] : used in time delays vs magnification ratio distribution of lensed events and their implications for lensing event candidate validation. 
  - [@Ng2025] : used in forecasting sub threshold lensed event detection with current and future GW detectors.
  - [@more2025] : used in forecasting EM-counterpart (GRBs) inference for lensed events with current and future GRB detectors.
  - [@hannuksela2025] : used in generation of various statistics of lensed events and their implications for lensing events.

- ler has been downloaded over 3000 times from PyPi over the past few years (showing active user base) and has been cited in multiple research papers (showing the impact of the package in the research community).
- users include graduate and undergraduate students, postdoctoral researchers, and faculty members in the field of astrophysics and gravitational wave research, indicating a broad user base and community engagement.


# AI usage disclosure

<!-- AI usage disclosure should include: Transparent disclosure of any use of generative AI in the software creation, documentation, or paper authoring. If no AI tools were used, state this explicitly. If AI tools were used, describe how they were used and how the quality and correctness of AI-generated content was verified. -->

No generative AI tools were used in the development of this software, the writing of this manuscript, or the preparation of supporting materials.


# Acknowledgements

The authors express their sincere appreciation for the significant contributions that have been instrumental in completing this research. Special thanks are extended to the academic advisors for their invaluable guidance and steadfast support. Acknowledgement is given to the Department of Physics, The Chinese University of Hong Kong, for the Postgraduate Studentship that facilitated this research. Hemantakumar Phurailatpam and Otto A. Hannuksela acknowledge support by grants from the Research Grants Council of Hong Kong (Project No. CUHK 14304622 and 14307923), the start-up grant from the Chinese University of Hong Kong, and the Direct Grant for Research from the Research Committee of The Chinese University of Hong Kong. Further gratitude is extended to the Netherlands Organisation for Scientific Research (NWO) for their support. N. Singh and D. Keitel are supported by Universitat de les Illes Balears (UIB); the Spanish Agencia Estatal de Investigación grants CNS2022-135440, PID2022-138626NB-I00, RED2022-134204-E, RED2022-134411-T, funded by MICIU/AEI/10.13039/501100011033, the European Union NextGenerationEU/PRTR, and the ERDF/EU; and the Comunitat Autònoma de les Illes Balears through the Direcció General de Recerca, Innovació I Transformació Digital with funds from the Tourist Stay Tax Law (PDR2020/11 - ITS2017-006) as well as through the Conselleria d’Economia, Hisenda i Innovació with grant numbers SINCO2022/6719 (European Union NextGenerationEU/PRTR-C17.I1) and SINCO2022/18146 (co-financed by the European Union and FEDER Operational Program 2021-2027 of the Balearic Islands). The authors also recognize the contributions of individuals who added empirical depth to this work. Appreciation is conveyed for the computational resources provided by the LIGO Laboratory, supported by National Science Foundation Grants No. PHY-0757058 and No. PHY-0823459.

# References
