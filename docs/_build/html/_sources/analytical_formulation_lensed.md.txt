# Analytical Formulation for Lensed Gravitational Wave Event Rates

Written by [Phurailatpam Hemantakumar](https://hemantaph.com). Last updated: 2026-02-06.

## Overview

This document extends the framework in [Analytical Formulation for Gravitational Wave Event Rates](https://ler.hemantaph.com/analytical_formulation_unlensed.html) to incorporate strong gravitational lensing. The goal is to compute the rate of detectable lensed gravitational-wave (GW) events from compact-binary coalescences (CBCs) by combining the intrinsic merger population with a lens population model, multi-image lensing geometry, and a detection criterion applied to the lensed images. The formulation introduces the optical depth for strong lensing, conditional distributions for source and lens properties given strong lensing, the multi-image caustic cross-section for an EPL+Shear lens model, and the probability of detecting at least two lensed images in a detector network.


## Table of Contents

- [Introduction](#introduction)  
- [Parameter-Marginalized Event Rate](#parameter-marginalized-event-rate)
- [Decomposing the Joint Source–Lens Parameter Distribution](#decomposing-the-joint-sourcelens-parameter-distribution)
- [Source Redshift Distributions of Lensed Events](#source-redshift-distributions-of-lensed-events)
- [Optical Depth](#optical-depth)
- [Lens Redshift Distribution of Lensed Events](#lens-redshift-distribution-of-lensed-events)
- [Multi-Image Caustic Cross-Section](#multi-image-caustic-cross-section)
- [Lens Parameter Distributions of Lensed Events](#lens-parameter-distributions-of-lensed-events)
- [Source Position Distribution and Image Properties](#source-position-distribution-and-image-properties)
- [Detection Probability of Lensed Events](#detection-probability-of-lensed-events)

## Introduction

The annual rate of detectable strongly lensed GW events, $\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}$, is the expected number of observed lensed CBC mergers per year in a given detector network. It is obtained by scaling the total intrinsic (unlensed) merger rate in the detector frame, $\frac{\Delta N_{\rm U}}{\Delta t}$, by the joint probability that an intrinsic merger is strongly lensed and detected, denoted $P({\rm SL,obs})$:

$$ 
\begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}
= \frac{\Delta N_{\rm U}}{\Delta t}\, P({\rm SL,obs}) \, .
\end{split}
\end{equation} 
$$


Using the product rule, $P({\rm SL,obs})=P({\rm SL})\,P({\rm obs}\mid{\rm SL})$, where $P({\rm SL})$ is the probability that an intrinsic merger is strongly lensed and $P({\rm obs}\mid{\rm SL})$ is the conditional probability of detection given strong lensing. This yields

$$ 
\begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}
= \frac{\Delta N_{\rm U}}{\Delta t}\, P({\rm SL})\, P({\rm obs}\mid{\rm SL})
= {\cal N}_{\rm L}\, P({\rm obs}\mid{\rm SL}) \, ,
\end{split}
\end{equation} 
$$

where ${\cal N}_{\rm L}\equiv \frac{\Delta N_{\rm U}}{\Delta t}\,P({\rm SL})$ is the total intrinsic rate of strongly lensed mergers, irrespective of detectability.


## Parameter-Marginalized Event Rate

The observed lensed event rate is obtained by averaging the detection probability over the unlensed GW parameters $\vec{\theta}_{\rm U}$, the lens parameters $\vec{\theta}_{\rm L}$, and the source position $\vec{\beta}$ within the lens caustic, conditioned on the occurrence of strong lensing. The unlensed GW parameters follow the definitions provided in the unlensed framework, with intrinsic and extrinsic components described in [Analytical Formulation for Gravitational Wave Event Rates](https://ler.hemantaph.com/analytical_formulation_unlensed.html).

For the Elliptical Power Law with external shear (EPL+Shear) lens model, the lens-parameter vector is

$$ 
\begin{equation}
\begin{split}
\vec{\theta}_{\rm L}\in \{z_L,\sigma,q,\phi_{\rm rot},\gamma,\gamma_1,\gamma_2\} \, ,
\end{split}
\end{equation} 
$$

where $z_L$ is the lens redshift, $\sigma$ is the velocity dispersion, $q$ is the projected axis ratio, $\phi_{\rm rot}$ is the major-axis orientation, $\gamma$ is the EPL logarithmic density slope, and $(\gamma_1,\gamma_2)$ are the external shear components. The source position on the source plane is $\vec{\beta}\in\{\beta_x,\beta_y\}$, defined relative to the lens center. The adopted priors are summarized in [Lens Parameter Priors](#lens-parameter-priors).

Given a joint prior distribution $P(\vec{\theta}_{\rm U},\vec{\theta}_{\rm L},\vec{\beta}\mid{\rm SL})$ conditioned on strong lensing,  and a detection probability model $P({\rm obs}\mid \vec{\theta}_{\rm U},\vec{\theta}_{\rm L},\vec{\beta},{\rm SL})$ for a specific source-lens configuration, the population-averaged detection probability becomes

$$ 
\begin{equation}
\begin{split}
P({\rm obs}\mid{\rm SL})
= \int P({\rm obs}\mid \vec{\theta}_{\rm U},\vec{\theta}_{\rm L},\vec{\beta},{\rm SL})\,
P(\vec{\theta}_{\rm U},\vec{\theta}_{\rm L},\vec{\beta}\mid{\rm SL})\,
d\vec{\theta}_{\rm U}\,d\vec{\theta}_{\rm L}\,d\vec{\beta} \, .
\end{split}
\end{equation} 
$$


To evaluate this integral with [Monte Carlo integration](https://en.wikipedia.org/wiki/Monte_Carlo_integration), the unlensed GW parameters are split into the source redshift $z_s$ and the remaining set $\vec{\theta}^{*}_{\rm U}$,

$$ 
\begin{equation}
\begin{split}
\vec{\theta}^{*}_{\rm U}\in \{m_1,m_2,a_1,a_2,\theta_1,\theta_2,\phi_{12},\phi_{\rm JL},\iota,\phi,\psi,t_c, {\rm RA},{\rm Dec}\} \, ,
\end{split}
\end{equation} 
$$

The lens parameters are split into the lens redshift $z_L$ and the structural parameters $\vec{\theta}^{*}_{\rm L}=\{\sigma,q,\phi_{\rm rot},\gamma,\gamma_1,\gamma_2\}$. Under the assumption of isotropic sky locations for sources and lenses, the strong-lensing-conditioned distribution of $(z_s,z_L,\vec{\theta}^{*}_{\rm L})$ is independent of $({\rm RA},{\rm Dec})$, so that

$$ 
\begin{equation}
\begin{split}
P(z_s,z_L,\vec{\theta}^{*}_{\rm L}\mid {\rm SL})
\equiv P(z_s,z_L,\vec{\theta}^{*}_{\rm L}\mid {\rm RA},{\rm Dec},{\rm SL}) \, .
\end{split}
\end{equation} 
$$


With these decompositions, the lensed detectable rate can be written as the expectation value

$$ 
\begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}
\simeq {\cal N}_{\rm L}
\bigg\langle
P({\rm obs}\mid \vec{\theta}^{*}_{\rm U},z_s,z_L,\vec{\theta}^{*}_{\rm L},\vec{\beta},{\rm SL})
\bigg\rangle_{\substack{
\vec{\theta}^{*}_{\rm U}\sim P(\vec{\theta}^{*}_{\rm U}) \\
(z_s,z_L,\vec{\theta}^{*}_{\rm L})\sim P(z_s,z_L,\vec{\theta}^{*}_{\rm L}\mid {\rm SL}) \\
\vec{\beta}\sim P(\vec{\beta}\mid z_s,z_L,\vec{\theta}^{*}_{\rm L},{\rm SL})
}} \, .
\end{split}
\end{equation} 
$$

The expectation value is evaluated by hierarchical sampling from the listed distributions. The sky location $({\rm RA},{\rm Dec})$ is drawn once from $P({\rm RA},{\rm Dec})$ and is then held fixed while sampling $(z_s,z_L,\vec{\theta}^{*}_{\rm L})$ and $\vec{\beta}$ for the same event.


## Decomposing the Joint Source–Lens Parameter Distribution

Hierarchical sampling of the lensing configuration is enabled by the chain-rule decomposition, i.e.

$$ 
\begin{equation}
\begin{split}
P(z_s,z_L,\vec{\theta}^{*}_{\rm L}\mid{\rm SL})
= P(z_s\mid{\rm SL})\,P(z_L\mid z_s,{\rm SL})\,P(\vec{\theta}^{*}_{\rm L}\mid z_L,z_s,{\rm SL}) \, .
\end{split}
\end{equation} 
$$

The lens-parameter factor can be expressed using Bayes’ theorem as

$$ 
\begin{equation}
\begin{split}
P(\vec{\theta}^{*}_{\rm L}\mid z_L,z_s,{\rm SL})
&= \frac{P({\rm SL}\mid z_L,z_s,\vec{\theta}^{*}_{\rm L})\,P(\vec{\theta}^{*}_{\rm L}\mid z_L,z_s)}{P({\rm SL}\mid z_L,z_s)} \,,\\ 
&\propto \sigma_{\rm SL}(z_L,z_s,\vec{\theta}^{*}_{\rm L})\,P(\vec{\theta}^{*}_{\rm L}\mid z_L,z_s) \, ,
\end{split}
\end{equation} 
$$


where $\sigma_{\rm SL}$ is the strong-lensing angular cross-section ([multi-image caustic cross-section](#multi-image-caustic-cross-section) in this documentation) and $P(\vec{\theta}^{*}_{\rm L}\mid z_L,z_s)$ is the prior distribution of lens parameters at fixed redshifts. This factorization motivates drawing $z_s$ first, then $z_L$, followed by sampling $\vec{\theta}^{*}_{\rm L}$ with importance sampling to account for the cross-section weighting.


## Source Redshift Distributions of Lensed Events

The redshift distribution of lensed sources differs from the intrinsic distribution because the probability of strong lensing increases with source distance. By Bayes’ theorem,

$$ 
\begin{equation}
\begin{split}
P(z_s\mid{\rm SL})
= \frac{P({\rm SL}\mid z_s)\,P(z_s)}{P({\rm SL})}
\propto P({\rm SL}\mid z_s)\,P(z_s) \, .
\end{split}
\end{equation} 
$$

Using the intrinsic unlensed redshift distribution,

$$ 
\begin{equation}
\begin{split}
P(z_s)
= \frac{1}{{\cal N}_{\rm U}}\,
\frac{R_{\rm U}(z_s)}{1+z_s}\,
\frac{dV_c}{dz_s} \, ,
\end{split}
\end{equation} 
$$

the lensed-source redshift distribution becomes

$$ 
\begin{equation}
\begin{split}
P(z_s\mid{\rm SL})
= \frac{1}{{\cal N}_{\rm L}}\,
P({\rm SL}\mid z_s)\,
\frac{R_{\rm U}(z_s)}{1+z_s}\,
\frac{dV_c}{dz_s} \, ,
\end{split}
\end{equation} 
$$

where ${\cal N}_{\rm U}$ is the total intrinsic merger rate per unit detector-frame time and ${\cal N}_{\rm L}$ is the total intrinsic strongly lensed merger rate per unit detector-frame time. The normalization satisfies

$$ 
\begin{equation}
\begin{split}
{\cal N}_{\rm L}
= {\cal N}_{\rm U}\,P({\rm SL})
= \int P({\rm SL}\mid z_s)\,
\frac{R_{\rm U}(z_s)}{1+z_s}\,
\frac{dV_c}{dz_s}\,dz_s \, ,
\end{split}
\end{equation} 
$$

so that $P({\rm SL})={\cal N}_{\rm L}/{\cal N}_{\rm U}$, which is the fraction of intrinsically occurring mergers that are strongly lensed in the population model. In this work, $P({\rm SL}\mid z_s)$ is computed for the adopted EPL+Shear lens population as detailed in the subsequent section. Once this probability is determined and the conditional distribution $P(z_s\mid {\rm SL})$ is established, source redshifts are generated using [inverse transform sampling](https://en.wikipedia.org/wiki/Inverse_transform_sampling).


## Optical Depth

In general radiative-transfer language, the optical depth $\tau(z_s)$ is the line-of-sight integral of an effective opacity and can be interpreted as the expected number of interactions (absorption, scattering, etc.) experienced by a signal as it propagates through a medium. By direct analogy in strong gravitational lensing, $\tau(z_s)$ is the expected number of strong-lensing encounters for a source at redshift $z_s$ under a randomly distributed lens population, where an encounter corresponds to a lens–source configuration in which the source lies inside the multi-image caustic and multiple images form. The optical depth is dimensionless and set by the lens number density and the strong-lensing cross-section. In the low-optical-depth regime $\tau(z_s)\ll 1$ relevant for galaxy-scale lensing, the probability that the source is strongly lensed satisfies $P({\rm SL}\mid z_s)\simeq \tau(z_s)$.

A convenient way to derive this relationship is to treat strong-lensing encounters as a Poisson process along the line of sight. If $\lambda$ denotes the expected number of encounters, the probability of observing exactly $k$ encounters is given by

$$ 
\begin{equation}
\begin{split}
P(k) = \frac{\lambda^k e^{-\lambda}}{k!} \, .
\end{split}
\end{equation}
$$

For strong lensing, the Poisson mean $\lambda$ corresponds to the integrated encounter expectation out to $z_s$, obtained by summing the angular target size of each potential lens over the lens population along the path. This gives

$$ 
\begin{equation}
\begin{split}
\lambda
&= \int ({\rm target\ size}) \times ({\rm number\ density}) \times ({\rm path\ length}) \\
&= \int_{0}^{z_s} \int_{\vec{\theta}^{*}_{\rm L}}
\frac{\sigma_{\rm SL}(z_L, z_s, \vec{\theta}^{*}_{\rm L})}{4\pi}\,
\frac{d^2N(z_L, \vec{\theta}^{*}_{\rm L})}{dV_c\,d\vec{\theta}^{*}_{\rm L}}\,
\frac{dV_c}{dz_L}\,
d\vec{\theta}^{*}_{\rm L}\,dz_L \\
&\equiv \tau(z_s) \, .
\end{split}
\end{equation}
$$

Geometrically, $\tau(z_s)$ represents the fraction of the sky covered by the effective strong-lensing cross-sections of all lenses along the line of sight, calculated under the single-lens approximation.

The probability of encountering exactly one strong lens is therefore

$$ 
\begin{equation}
\begin{split}
P({\rm SL}\mid z_s)
= P(k=1\mid \lambda=\tau(z_s))
= \tau(z_s)\,e^{-\tau(z_s)} \, .
\end{split}
\end{equation}
$$

Since $\tau(z_s)\ll 1$ for galaxy populations on cosmological volumes, this reduces to $P({\rm SL}\mid z_s)\simeq \tau(z_s)$. A rough estimate using `ler` shows that typically $\tau(z_s) \lesssim 0.004$, with the relative difference $\frac{P({\rm SL} \mid z_s) - \tau(z_s)}{P({\rm SL} \mid z_s)} \times 100 \lesssim 0.4\%$.

In practice, the lens population is often specified by the velocity-dispersion function $\phi(\sigma,z_L)=\frac{d^2N(z_L,\sigma)}{dV_c\,d\sigma}$. Separating the remaining lens (caustic) shape parameters $\vec{\theta}^{**}_{\rm L}$ from $\sigma$, the strong-lensing probability can be written as

$$ 
\begin{equation}
\begin{split}
P({\rm SL}\mid z_s)
&\simeq \int_{0}^{z_s}\int_{\sigma,\vec{\theta}^{**}_{\rm L}}
\frac{\sigma_{\rm SL}(z_L, z_s, \sigma, \vec{\theta}^{**}_{\rm L})}{4\pi}\,
\phi(\sigma,z_L)\\
&\quad\quad\quad\quad \times P(\vec{\theta}^{**}_{\rm L}\mid\sigma)\,
\frac{dV_c}{dz_L}\,
d\sigma\,d\vec{\theta}^{**}_{\rm L}\,dz_L \, .
\end{split}
\end{equation}
$$

It is convenient to introduce the differential optical depth with respect to lens redshift,

$$ 
\begin{equation}
\begin{split}
P({\rm SL}\mid z_s)
\simeq \int_{0}^{z_s} \frac{d\tau}{dz_L}\,dz_L \, ,
\end{split}
\end{equation}
$$

where

$$ 
\begin{equation}
\begin{split}
\frac{d\tau}{dz_L}
= \int_{\sigma,\vec{\theta}^{**}_{\rm L}}
\frac{\sigma_{\rm SL}(z_L, z_s, \sigma, \vec{\theta}^{**}_{\rm L})}{4\pi}\,
\phi(\sigma,z_L)\,
P(\vec{\theta}^{**}_{\rm L}\mid\sigma)\,
\frac{dV_c}{dz_L}\,
d\sigma\,d\vec{\theta}^{**}_{\rm L} \, .
\end{split}
\end{equation}
$$

The function $d\tau/dz_L$ is used both to evaluate $\tau(z_s)$ and to construct the lens-redshift distribution $P(z_L\mid z_s,{\rm SL})$ in the next section.

For the EPL+Shear model, the shape-parameter set is $\vec{\theta}^{**}_{\rm L}=\{q,\phi_{\rm rot},\gamma,\gamma_1,\gamma_2\}$, and the differential optical depth becomes

$$ 
\begin{equation}
\begin{split}
\frac{d\tau}{dz_L}
&= \int_{\sigma,q,\phi_{\rm rot},\gamma,\gamma_1,\gamma_2} \frac{\sigma_{\rm SL}(z_L, \sigma, q,\phi_{\rm rot},\gamma,\gamma_1,\gamma_2)}{4\pi} \phi(\sigma,z_L) \\
& \quad \times P(\gamma_1, \gamma_2) \, P(\gamma) \, P(\phi_{\rm rot}) \, P(q \mid \sigma)  \frac{dV_c}{dz_L} \, d\sigma \, dq \, d\phi_{\rm rot} \, d\gamma \, d\gamma_1 \, d\gamma_2\, .
\end{split}
\end{equation} 
$$

Numerically, this integral is approximated using Monte Carlo integration. To simplify the sampling process, we employ a uniform proposal distribution for the velocity dispersion, $P_o(\sigma) = 1/\Delta\sigma$, defined over the range of interest $\Delta\sigma$. Rewriting the integral as an expectation value, the differential optical depth reads

$$ 
\begin{equation}
\begin{split}
\frac{d\tau}{dz_L}
= \Delta\sigma
\left\langle
\frac{\sigma^{\rm EPL}_{\rm SL}}{4\pi}\,
\phi(\sigma,z_L)\,
\frac{dV_c}{dz_L}
\right\rangle_{\substack{
\sigma\sim P_o(\sigma) \\
q\sim P(q\mid\sigma) \\
\phi_{\rm rot}\sim P(\phi_{\rm rot}) \\
\gamma\sim P(\gamma) \\
\gamma_1,\gamma_2\sim P(\gamma_1,\gamma_2)
}} \, ,
\end{split}
\end{equation} 
$$

where $\sigma^{\rm EPL}_{\rm SL}$ denotes the multi-image caustic cross-section for the EPL+Shear model, and each parameter is sampled hierarchically from its respective prior.

## Lens Redshift Distribution of Lensed Events

For fixed $z_s$, the lens-redshift distribution conditioned on the occurance of strong lensing is directly proportional to the differential optical depth, and is given by

$$ 
\begin{equation}
\begin{split}
P(z_L\mid z_s,{\rm SL})
= \frac{1}{\tau(z_s)}\,\frac{d\tau}{dz_L} \, ,
\end{split}
\end{equation} 
$$

where the normalization factor is the total optical depth $\tau(z_s)=\int_{0}^{z_s}(d\tau/dz_L)\,dz_L$. In `ler`, the differential optical depth $d\tau/dz_L$ is pre-computed on a grid spanning the relevant redshift space ($z_L, z_s$). This grid enables efficient interpolation for inverse-transform sampling of $z_L$ and allows for rapid evaluation of the strong-lensing probability $P({\rm SL}\mid z_s)$.

For comparison, the intrinsic redshift distribution of the lens population, without conditioning on strong lensing and independent of any specific background source, is obtained by marginalizing the lens number density over the lens parameters

$$ 
\begin{equation}
\begin{split}
P(z_L)
= \frac{1}{{\cal N}_{z_L}}
\int
\frac{d^2N(z_L,\vec{\theta}^{*}_{\rm L})}{dV_c\, d\vec{\theta}^{*}_{\rm L}}\,
\frac{dV_c}{dz_L}\,
d\vec{\theta}^{*}_{\rm L} \, ,
\end{split}
\end{equation} 
$$

where ${\cal N}_{z_L}$ is a normalization constant. 

The impact of the optical depth on the source redshift distributions of the lensed events is illustrated in Figure 1.

<div align="center">
  <img src="_static/lensed_redshifts.png" alt="Merger rate density and PDF of redshift for BBH mergers" width="600"/>
</div>

<div id="fig1"></div>

>**Figure 1:** Redshift distributions (blue, left axis) for the intrinsic source population $P(z=z_s)$, the lensed-source population $P(z=z_s\mid{\rm SL})$, and the conditional lens distribution $P(z=z_L\mid z_s,{\rm SL})$, shown together with the optical depth $\tau(z=z_s)$ (orange, right axis), which quantifies the strong-lensing probability for a source at redshift $z_s$ in the low-optical-depth regime. The intrinsic $P(z_s)$ follows the unlensed merger-rate model, while $P(z_s\mid{\rm SL})$ is obtained by reweighting the intrinsic population by the lensing probability, $P(z_s\mid{\rm SL})\propto \tau(z_s),P(z_s)$, implemented via rejection sampling using $\tau(z_s)$. For each lensed source, the lens redshift is drawn from $0<z_L<z_s$ with a density proportional to $\frac{d\tau}{dz_L}$, so that $P(z_L\mid z_s,{\rm SL})$ traces where along the line of sight the lenses that produce strong lensing are most likely to occur. Since $\tau(z_s)$ increases with redshift, strong lensing preferentially selects more distant sources and shifts $P(z_s\mid{\rm SL})$ to higher redshift relative to $P(z_s)$, while $P(z_L\mid z_s,{\rm SL})$ describes the redshift distribution of the lenses responsible for that lensed-source population.


## Multi-Image Caustic Cross-Section

The strong-lensing condition is satisfied when a source lies within the caustic region of the lens, resulting in multiple images. The angular area of this region is the multi-image caustic cross-section $\sigma_{\rm SL}$ (or $\sigma^{\rm EPL}_{\rm SL}$ for the EPL+Shear model). The probability of strong lensing for a specific configuration is the ratio of this cross-section to the total sky area

$$ 
\begin{equation}
\begin{split}
P({\rm SL}\mid z_L, z_s, \vec{\theta}^{*}_{\rm L})
= \frac{\sigma_{\rm SL}}{4\pi} \, .
\end{split}
\end{equation} 
$$

Direct evaluation of $\sigma_{\rm SL}$ by tracing caustics with [lenstronomy](https://lenstronomy.readthedocs.io) is computationally expensive. In `ler`, $\sigma^{\rm EPL}_{\rm SL}$ is therefore precomputed on the grid of lens (caustic) shape parameters $\vec{\theta}^{**}_{\rm L}\in \{q,\phi_{\rm rot},\gamma,\gamma_1,\gamma_2\}$ at unit Einstein radius $\theta_{\rm E}=1$, and interpolated during sampling. For arbitrary $(\sigma,z_L,z_s)$, the Einstein radius is

$$ 
\begin{equation}
\begin{split}
\theta_{\rm E}
= 4\pi\left(\frac{\sigma}{c}\right)^2
\frac{D_{LS}(z_L,z_s)}{D_S(z_s)} \, ,
\end{split}
\end{equation} 
$$

where $D_{LS}$ and $D_S$ are angular-diameter distances.

The cross-section is then rescaled from the $\theta_{\rm E}=1$ values using the calibrated relation

$$ 
\begin{equation}
\begin{split}
\sigma^{\rm EPL}_{\rm SL}(\theta_{\rm E}, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2)
&= \sigma^{\rm EPL}_{\rm SL}(1, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2) \\
&\quad\times \left({\rm intercept} + {\rm slope}\,\pi\,\theta_{\rm E}^2\right) \, ,
\end{split}
\end{equation}
$$

where ${\rm intercept}$ and ${\rm slope}$ are precomputed constants. This interpolation and rescaling scheme accelerates cross-section evaluation and, consequently, the computation of $\frac{d\tau}{dz_L}$.

Going forward, I simply use the term 'cross-section' to refer to $\sigma_{\rm SL}$ or $\sigma^{\rm EPL}_{\rm SL}$.

## Lens Parameter Distributions of Lensed Events

For a fixed redshift pair $(z_L,z_s)$, the conditional distribution of lens parameters $P(\vec{\theta}^{*}_{\rm L}\mid z_L,z_s,{\rm SL})$ is determined by weighting the intrinsic lens population with the strong-lensing cross-section $\sigma^{\rm EPL}_{\rm SL}(z_L,z_s,\vec{\theta}^{*}_{\rm L})$. For the EPL+Shear model, this target distribution is expressed as

$$
\begin{equation}
\begin{split}
&P(\sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2 \mid z_L, z_s, {\rm SL})\\
&\quad \propto P({\rm SL}\mid z_L, z_s, \sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2) \, P(\gamma_1, \gamma_2) \, P(\gamma) \, P(\phi_{\rm rot}) \, P(q \mid \sigma) \, P(\sigma \mid z_L) \\
&\quad \propto \frac{\sigma^{\rm EPL}_{\rm SL}}{4\pi} \, P(\gamma_1, \gamma_2) \, P(\gamma) \, P(\phi_{\rm rot}) \, P(q \mid \sigma) \, P(\sigma \mid z_L) \, .
\end{split}
\end{equation}
$$



Direct rejection sampling from the intrinsic distribution $P(\sigma\mid z_L)$ is computationally inefficient because the lensing cross-section $\sigma^{\rm EPL}_{\rm SL}$ scales roughly as $\sigma^4$, heavily favoring high-mass lenses that are exponentially rare in the intrinsic population. Furthermore, the lack of a strict analytical maximum for the cross-section complicates the definition of a rejection envelope.

To address these challenges, `ler` employs [importance sampling](https://en.wikipedia.org/wiki/Importance_sampling). A proposal distribution $P_o(\sigma)$, typically chosen to be uniform over $[\sigma_{\min},\sigma_{\max}]$, is introduced to rewrite the target density. By defining the intrinsic distribution $P(\sigma \mid z_L)$ as the normalized velocity dispersion function $\phi(\sigma, z_L)/\int \phi(\sigma', z_L) d\sigma'$, the target distribution becomes

$$
\begin{equation}
\begin{split}
&P(\sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2 \mid z_L, z_s, {\rm SL})\\
&\quad \propto \left[\frac{ \sigma^{\rm EPL}_{\rm SL}\, P(\sigma \mid z_L)}{4\pi\,P_o(\sigma)}\right] P(\gamma_1, \gamma_2) \, P(\gamma) \, P(\phi_{\rm rot}) \, P(q \mid \sigma) \, P_o(\sigma) \\
&\quad = \frac{1}{{\cal W}(z_L)} \left[\frac{ \sigma^{\rm EPL}_{\rm SL}\, \phi(\sigma, z_L)}{P_o(\sigma)}\right] P(\gamma_1, \gamma_2) \, P(\gamma) \, P(\phi_{\rm rot}) \, P(q \mid \sigma) \, P_o(\sigma) \, .
\end{split}
\end{equation}
$$

In this formulation, constant factors including $4\pi$ and the normalization of the intrinsic velocity dispersion function are absorbed into the overall normalization constant ${\cal W}(z_L)$. The term in the square brackets represents the unnormalized importance weight required to convert samples from the proposal distribution into samples from the cross-section-weighted target distribution.

Computationally, this is implemented by generating a candidate batch of size $N_{\rm prop}$ (defaulting to 200 in `ler`) for each $(z_L, z_s)$ pair. The velocity dispersion $\sigma_i$ is drawn from the proposal $P_o(\sigma)$, while the remaining shape parameters are sampled directly from their intrinsic priors. Each candidate in the batch is assigned a normalized importance weight

$$ 
\begin{equation}
\begin{split}
w_i(\sigma_i, q_i, \phi_{{\rm rot},i}, \gamma_i, \gamma_{1,i}, \gamma_{2,i};z_L,z_s)=\frac{1}{{\cal W}(z_L)}\,\frac{\sigma^{\rm EPL}_{{\rm SL}, i}\,\phi(\sigma_i, z_L)}{P_o(\sigma_i)} \, ,
\end{split}
\end{equation} 
$$

where the batch normalization constant is defined as the sum of the unnormalized weights

$$ 
\begin{equation}
\begin{split}
{\cal W}(z_L)=\sum_{i=1}^{N_{\rm prop}}\frac{\sigma^{\rm EPL}_{{\rm SL}, i}\,\phi(\sigma_i, z_L)}{P_o(\sigma_i)} \, .
\end{split}
\end{equation} 
$$

A final lens configuration is selected from this batch by resampling with probabilities proportional to $w_i$, ensuring the final parameter set accurately reflects the properties of the strongly lensed population.


## Source Position Distribution and Image Properties

Given $(z_s,z_L,\vec{\theta}^{*}_{\rm L})$, the source position $\vec{\beta}$ is sampled uniformly within the multi-image caustic so that each draw produces multiple images. The lens equation is solved with [`lenstronomy`](https://lenstronomy.readthedocs.io) to obtain the image positions $\vec{\theta}_i\in \{\theta_{x,i},\theta_{y,i}\}$, magnifications $\mu_i$, and arrival-time delays $t_i$. These quantities are initially computed at $\theta_{\rm E}=1$ and then rescaled to physical units using the system’s $\theta_{\rm E}$.

Each image is assigned a Morse index $n_i\in\{0,1/2,1\}$ corresponding to minima (Type I), saddles (Type II), and maxima (Type III), respectively. Type III images are typically highly demagnified and are neglected in the detectability criterion adopted here. Lensing modifies the observed GW parameters for each image $i$ as follows


$$ 
\begin{equation}
\begin{split}
d_L &\rightarrow \frac{d_L(z_s)}{\sqrt{\lvert\mu_i\rvert}}\,, \\
\phi_c &\rightarrow \phi_c - n_i\pi\,, \\
t_c &\rightarrow t_c + t_i \, , \\
{\rm RA} &\rightarrow {\rm RA} + \frac{\theta_{x,i}-\beta_x}{\cos({\rm Dec})}\,, \\
{\rm Dec} &\rightarrow {\rm Dec} + (\theta_{y,i}-\beta_y)\,, \\
\end{split}
\end{equation} 
$$


while the other extrinsic and intrinsic parameters are unchanged from the unlensed values, including the detector-frame component masses $\left(m_1 \, (1+z_s),m_2 \, (1+z_s)\right)$. The sky-location offsets associated with $\vec{\theta}_i-\vec{\beta}$ are typically negligible for current GW detector networks and have no impact on SNR. Detectability is therefore evaluated by computing the SNR for each image using the modified GW pararameters $\vec{\theta}_{{\rm GW},i}$, as described in the next section.

> **Note:** The lens-centre coordinates in the sky is given by $\left\{ {\rm RA} - \frac{\beta_x}{\cos({\rm Dec})}, {\rm Dec} - \beta_y \right\}$.


## Detection Probability of Lensed Events

For a specific lensed configuration, the detection probability is determined by the detectability of its images. The signal-to-noise ratio (SNR), $\rho(\vec{\theta}_{{\rm GW},i})$, is computed for each image using the modified GW parameters. An event is considered detected if at least two images meet the detection criterion. The probability is therefore a binary condition

$$
\begin{equation}
\begin{split}
&P({\rm obs} \mid \vec{\theta}^*_{\rm U}, z_s, \vec{\theta_{\rm L}}, \vec{\beta}, {\rm SL}) = \Theta\left[\left(\sum_i P_{\rm det}(\vec{\theta}_{{\rm GW},i}, \rho_{\rm th}) \right) - 2\right] \, ,
\end{split}
\end{equation}
$$

where $\Theta$ is the Heaviside step function and the sum is over all images $i$. In the simplest step-function model, the individual image detection probability is defined as

$$
\begin{equation}
\begin{split}
P_{\rm det} (\vec{\theta}, \rho_{\rm th}) = 
\Theta[\rho(\vec{\theta}_{{\rm GW},i}) - \rho_{\rm th}] \, ,
\end{split}
\end{equation}
$$

and in `ler`, it is estimated using the [gwsnr](https://gwsnr.hemantaph.com) package. 

## Complete Expression for Lensed Event Rate

Summarizing the derivation steps, the annual rate of detectable strongly lensed GW events is given by

$$
\begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}
&=
{\cal N}_{\rm L}
\bigg\langle
P({\rm obs}\mid \vec{\theta}^{*}_{\rm U}, z_s, \vec{\theta}_{\rm L}, \vec{\beta}, {\rm SL})
\bigg\rangle_{\substack{
\vec{\theta}^{*}_{\rm U}\sim P(\vec{\theta}^{*}_{\rm U}) \\
z_s \sim P(z_s \mid {\rm SL}) \\
z_L \sim P(z_L \mid z_s, {\rm SL}) \\
\vec{\theta}^{*}_{\rm L} \sim P(\vec{\theta}^{*}_{\rm L} \mid z_L, z_s, {\rm SL}) \\
\vec{\beta} \sim P(\vec{\beta} \mid z_s, \vec{\theta}_{\rm L}, {\rm SL})
}} \, .
\end{split}
\end{equation}
$$

The computational procedure for evaluating the rate via Monte Carlo integration is summarized as follows:

1. **Intrinsic GW Parameters:** Sample the unlensed GW parameters $\vec{\theta}^{*}_{\rm U}$ from their intrinsic priors $P(\vec{\theta}^{*}_{\rm U})$.

2. **Lensed Source Redshift:** Sample the source redshift $z_s$ via inverse-transform sampling from the lensed distribution $P(z_s \mid {\rm SL})$, which is derived from the intrinsic redshift distribution $P(z_s)$ and the precomputed optical depth $\tau(z_s)$.

3. **Lens Redshift:** Conditioned on the sampled $z_s$, sample the lens redshift $z_L$ via inverse-transform sampling from the distribution $P(z_L \mid z_s, {\rm SL})$, which is determined numerically from the differential optical depth $d\tau/dz_L$.

4. **Lens Properties:** Conditioned on $(z_s, z_L)$, generate a candidate batch of lens parameters $\{ \vec{\theta}^{*}_{{\rm L},i} \}_{i=1}^{N_{\rm prop}}$. Draw shape parameters $\vec{\theta}^{**}_{{\rm L},i}$ from their intrinsic priors $P(\vec{\theta}^{**}_{\rm L} \mid z_L, z_s)$ and velocity dispersions $\sigma_i$ from a uniform proposal $P_o(\sigma)$. Select a single final configuration from this batch using importance sampling with normalized weights $w_i = \frac{1}{{\cal W}(z_L)}\frac{\sigma^{\rm EPL}_{{\rm SL},i}\,\phi(\sigma_i\mid z_L)}{P_o(\sigma_i)}$, where ${\cal W}(z_L)$ ensures $\sum_i w_i=1$.

5. **Source Position:** Sample the source position $\vec{\beta}$ uniformly within the multi-image caustic region defined by the selected lens parameters, calculated at a reference Einstein radius $\theta_{\rm E}=1$.

6. **Image Properties:** Solve the lens equation to obtain the normalized image positions $\vec{\theta}_i$, magnifications $\mu_i$, and time delays $t_i$. Rescale the geometric parameters $\vec{\beta}$, $\vec{\theta}_i$, and $t_i$ to physical units using the physical Einstein radius $\theta_{\rm E}(\sigma, z_L, z_s)$.

7. **Observed GW Parameters:** Construct the observed parameters $\vec{\theta}_{{\rm GW},i}$ for each image by applying lensing transformations (magnification, time delay, Morse phase shift and source position offsets) to the unlensed parameters $\vec{\theta}^{*}_{\rm U}$.

8. **Detection Criteria:** Compute the SNR $\rho(\vec{\theta}_{{\rm GW},i})$ for each image. An event is classified as detected if at least two images satisfy the detection threshold, i.e., $P({\rm obs}\mid \vec{\theta}^{*}_{\rm U}, z_s, \vec{\theta}_{\rm L}, \vec{\beta}, {\rm SL})=1$ if $\sum_i \rho(\vec{\theta}_{{\rm GW},i})\geq 2$ otherwise $0$.

9. **Rate Estimation:** Estimate the final observable rate by averaging the detection probability $P({\rm obs}\mid \vec{\theta}^{*}_{\rm U}, z_s, \vec{\theta}_{\rm L}, \vec{\beta}, {\rm SL})$ over a large ensemble (typically $\gtrsim 10^6$) of parameter sets $\{\vec{\theta}^{*}_{\rm U}, z_s, \vec{\theta}_{\rm L}, \vec{\beta} \}$ and scaling the result by the total lensed rate ${\cal N}_{\rm L}$.

The specific parameter priors and their functional forms are detailed below.


## Simulation Results

This section reports results (related to strong lensing) generated with the default configuration of the ler package, compare it with unlensed events. Although ler supports alternative population prescriptions, including user-defined models, all numbers and figures shown here use the standard settings and assumptions summarized below.

### Simulation Settings

All calculations assume a flat $\Lambda{\rm CDM}$ cosmology with $H_0 = 70\,{\rm km}\,{\rm s}^{-1}\,{\rm Mpc}^{-1}$, $\Omega_m = 0.3$, and $\Omega_\Lambda = 0.7$.

Detection probabilities are evaluated through the [`gwsnr`](https://gwsnr.hemantaph.com) backend. Waveforms are generated with the IMRPhenomXPHM approximant, and detector responses are computed with a sampling frequency of $2048\,{\rm Hz}$. The lower frequency cutoff is set to $f_{\rm low}=20\,{\rm Hz}$ for the LIGO–Virgo–KAGRA network at O4 design sensitivity and $f_{\rm low}=10\,{\rm Hz}$ for the third-generation (3G) network consisting of Einstein Telescope (ET) and Cosmic Explorer (CE). The detection threshold is set to $\rho_{\rm th}=10$ for both single detectors and networks. The detector duty cycle is assumed to be 100%.

### Lens and source Parameter Priors

GW Source Parameter Priors are kept the same as in the section [GW Source Parameter Priors](https://ler.hemantaph.com/analytical_formulation_unlensed.html#gw-source-parameter-priors) of [Analytical Formulation for Gravitational Wave Event Rates](https://ler.hemantaph.com/analytical_formulation_unlensed.html). The lens prior distributions and parameter ranges used in the event-rate calculations are summarized in Table 1, Table 2 and Table 3. These choices follows common conventions such as GW population forllowing GWTC-3–motivated models for source and SDSS and SLACS observations of galaxy populations.

#### Table 1: Redshifts

The source and lens redshifts are drawn from distributions conditioned on the occurrence of strong lensing. The optical depth $\tau(z_s)$ and differential optical depth $d\tau/dz_L$ are pre-computed based on the adopted EPL+Shear lens population and cosmology.

| Parameter | Unit | Prior <br>Distribution | Range <br>[Min, Max] | Description |
| :--- | :--- | :--- | :--- | :--- |
| $z_s$ | - | $P(z_s\mid {\rm SL})$ <br>$\propto \tau(z_s)\,P(z_s)$ | [0, 10] | Source redshift <br>conditioned on ${\rm SL}$ |
| $z_L$ | - | $P(z_L\mid z_s,{\rm SL})$ <br>$\propto \frac{d\tau}{dz_L}$ | [0, $z_s$] | Lens redshift <br>conditioned on $z_s$ and ${\rm SL}$ |

#### Table 2: Lens parameter priors

Lens parameters priors, intrinsic or proposed, are listed below. Shape parameters assume a local flat coordinate system where the x-axis and the y-axis aligns with sky coordinates RA and Dec, respectively. 

| Parameter | Unit | Prior <br>Distribution | Range <br>[Min, Max] | Description |
| :--- | :--- | :--- | :--- | :--- |
| $\sigma$ | km/s | **Proposed prior** <br>Uniform <br><br>**Physical prior** <br>$P(\sigma\mid z_L)\propto$ <br>$\phi(\sigma, z_L)$ | [100, 400] | Velocity dispersion. <br><br>**Proposed prior :** <br> Used for importance sampling <br>and computing $\frac{d\tau}{dz_L}$. <br><br>**Physical prior :** <br>Based on the velocity <br>dispersion function $\phi(\sigma, z_L)$ <br>(Oguri et al. 2018; SDSS). <br>It is used to calculate <br>importance weights and <br>$d\tau/dz_L$, rather than for <br>direct sampling. |
| $q$ | - | $P(q\mid \sigma)$ <br>Rayleigh <br>distribution | [0.2, 1.0] | Projected axis ratio, <br>conditioned on $\sigma$ <br>(Collett et al. 2015; SDSS). |
| $\phi_{\rm rot}$ | rad | Uniform | [0, $\pi$] | Lens orientation angle, <br>measured counter-clockwise <br>from the x-axis. |
| $\gamma$ | - | Normal | - | Density profile slope <br>(Sonnenfeld et al. 2024; SLACS). <br>Mean: 2.0, Std: 0.1 |
| $\gamma_1, \gamma_2$ | - | Normal | - | External shear <br>(Collett et al. 2015; SDSS).<br> x and y components. <br>Mean: 0.0, Std: 0.05 |

#### Table 3: Source position prior

Source position $\beta$ are defined in a local flat coordinate system with the x-axis aligned with RA and the y-axis aligned with Dec.

| Parameter | Unit | Prior Distribution | Description |
| :--- | :--- | :--- | :--- |
| $\beta_x$ | $\theta_{\rm E}$ | Uniform within<br> multi-image caustic | x-component <br>of source position from lens center |
| $\beta_y$ | $\theta_{\rm E}$ | Uniform within<br> multi-image caustic | y-component <br>of source position from lens center |

### Plot Comparison for Lensed+Detectable, Lensed and Intrinsic Populations

Selection effects in lensed gravitational-wave observations are illustrated by comparing the intrinsic BBH population to the subsets that satisfies the strong-lensing condition and the detection criteria. The former subset is referred to as the `lensed population`, while the latter is referred to as the `lensed+detectable population`. 

## Rate Estimates and Comparisons

Estimated detectable annual merger rates for BBH (Pop I–II)  populations are listed below for 

| Detector <br>Configuration | Unlensed <br>Merger Rate (${\rm yr}^{-1}$) | Lensed <br>Merger Rate (${\rm yr}^{-1}$) | Ratio <br>Unlensed:Lensed |
| :--- | :--- | :--- | :--- |
| [L1, H1, V1] (O4) | 299.56 | 0.10 | 2996:1 |