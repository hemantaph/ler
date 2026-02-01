# Analytical Formulation for Lensed Gravitational Wave Event Rates

Written by [Phurailatpam Hemantakumar](https://hemantaph.com).

## Overview

This documentation extends the [Analytical Formulation for Gravitational Wave Event Rates](https://ler.hemantaph.com/analytical_formulation_unlensed.html) to account for strong gravitational lensing effects. We present a comprehensive framework for computing the rates of detectable lensed gravitational-wave events from compact-binary mergers. The formulation addresses key aspects of lensed systems: the optical depth describing lensing probability, the distributions of source and lens properties conditioned on strong lensing, the multi-image lensing cross-sections, and the detection prospects for lensed signals. This approach enables accurate predictions of observable lensed merger rates in current and future gravitational-wave detector networks.


## Table of Contents

1. [Introduction](#introduction)
2. [Parameter-Marginalized Event Rate](#parameter-marginalized-event-rate)
3. [Decomposing the Joint Source–Lens Parameter Distribution](#decomposing-the-joint-source–lens-parameter-distribution)
4. [Source Redshift Distributions of Lensed Events](#source-redshift-distributions-of-lensed-events)
5. [Optical Depth](#optical-depth)
6. [Lens Redshift Distribution of Lensed Events](#lens-redshift-distribution-of-lensed-events)
7. [Multi-Image Lensing Cross-Section](#multi-image-lensing-cross-section)
8. [Lens Parameter Distributions of Lensed Events](#lens-parameter-distributions-of-lensed-events)
9. [Source Position Distribution and Image Properties](#source-position-distribution-and-image-properties)
10. [Detection Probability of Lensed Events](#detection-probability-of-lensed-events)
11. [Complete Expression for Lensed Event Rate](#complete-expression-for-lensed-event-rate)
12. [Lens Parameter Priors](#lens-parameter-priors)

## Introduction

The annual rate of detectable strongly lensed gravitational-wave (GW) events, $\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}$, represents the expected number of observed lensed compact-binary mergers per year for a specific detector network. This rate is derived by scaling the total intrinsic merger rate, $\frac{\Delta N_{\rm U}}{\Delta t}$, by the joint probability that a merger is both strongly lensed and detected, denoted as $P({\rm SL, obs})$. This relationship is expressed as

$$ \begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}
= \frac{\Delta N_{\rm U}}{\Delta t}\, P({\rm SL, obs}) \, .
\end{split}
\end{equation} $$

The joint probability is decomposed using the product rule into $P({\rm SL})\,P({\rm obs}\mid{\rm SL})$, where $P({\rm SL})$ is the probability of strong lensing occurring within the intrinsic population and $P({\rm obs}\mid{\rm SL})$ is the conditional probability of detecting an event given that it is strongly lensed. The rate equation then becomes

$$ \begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}
= \frac{\Delta N_{\rm U}}{\Delta t}\, P({\rm SL})\, P({\rm obs}\mid{\rm SL})
= {\cal N}_{\rm L}\, P({\rm obs}\mid{\rm SL}) \, ,
\end{split}
\end{equation} $$

where ${\cal N}_{\rm L} \equiv \frac{\Delta N_{\rm U}}{\Delta t}\, P({\rm SL})$ is defined as the total rate of strongly lensed mergers regardless of their detectability.


## Parameter-Marginalized Event Rate

The conditional detection probability $P({\rm obs}\mid{\rm SL})$ requires marginalization over the unlensed GW parameters $\vec{\theta}_{\rm U}$, the lens parameters $\vec{\theta_{\rm L}}$, and the source position $\vec{\beta}$. The GW parameters $\vec{\theta}_{\rm U}$ follow the definitions provided in the [Analytical Formulation for Gravitational Wave Event Rates](https://ler.hemantaph.com/analytical_formulation_unlensed.html). For the Elliptical Power Law with external shear (EPL+Shear) model used in this work, the lens parameter set is defined as

$$ \begin{equation}
\begin{split}
\vec{\theta_{\rm L}}\in \{z_L,\sigma,q,\phi_{\rm rot},\gamma,\gamma_1,\gamma_2\} \, ,
\end{split}
\end{equation} $$

which includes the lens redshift $z_L$, velocity dispersion $\sigma$, projected axis ratio $q$, axis rotation angle $\phi_{\rm rot}$, density-profile slope $\gamma$, and external shear components $(\gamma_1,\gamma_2)$. The source position $\vec{\beta}\in\{\beta_x,\beta_y\}$ represents the angular offset of the source from the lens center in the source plane.

The observable lensed event rate is calculated by integrating the detection probability over the parameter space weighted by the joint distribution conditioned on strong lensing. It reads

$$ \begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}
&= {\cal N}_{\rm L}
\int_{\vec{\theta_{\rm U}}}\int_{\vec{\theta_{\rm L}}}\int_{\vec{\beta}} P({\rm obs} \mid \vec{\theta_{\rm U}}, \vec{\theta_{\rm L}}, \vec{\beta}, {\rm SL}) \\
&\quad\quad\quad\quad\quad\quad\quad\quad \times P(\vec{\theta_{\rm U}}, \vec{\theta_{\rm L}}, \vec{\beta} \mid {\rm SL}) d\vec{\theta_{\rm U}}\, d\vec{\theta_{\rm L}}\, d\vec{\beta} \, .
\end{split}
\end{equation} $$

To facilitate sampling, the GW parameters are split into the redshift $z_s$ and the remaining intrinsic parameters $\vec{\theta}^{*}_{\rm U}$. Similarly, the lens parameters are split into the redshift $z_L$ and the structural parameters $\vec{\theta}^{*}_{\rm L}$. Assuming $\vec{\theta}^{*}_{\rm U}$ is independent of the lensing configuration once $z_s$ is fixed, the integral is evaluated via Monte Carlo sampling as

$$
\begin{equation}
\begin{split}
&\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}\\
&\quad\quad= {\cal N}_{\rm L}
\int_{\vec{\theta}_{\rm U}}\int_{\vec{\theta_{\rm L}}}\int_{\vec{\beta}}
 P({\rm obs} \mid \vec{\theta}^*_{\rm U}, z_s, z_L,\vec{\theta}^{*}_{\rm L}, \vec{\beta}, {\rm SL}) \\
&\quad\quad\quad\quad\quad\quad\quad\quad \times P(\vec{\theta}^*_{\rm U})\, P(z_s, z_L,\vec{\theta}^{*}_{\rm L} \mid {\rm SL}) \\
&\quad\quad\quad\quad\quad\quad\quad\quad \times P(\vec{\beta} \mid z_s, z_L,\vec{\theta}^{*}_{\rm L}, {\rm SL}) \, d\vec{\theta}_{\rm U}\, d\vec{\theta_{\rm L}}\, d\vec{\beta} 
\\
&\quad\quad= {\cal N}_{\rm L}
\left\langle
P({\rm obs}\mid \vec{\theta}^{*}_{\rm U},z_s,z_L,\vec{\theta}^{*}_{\rm L},\vec{\beta},{\rm SL})
\right\rangle_{\substack{
\vec{\theta}^{*}_{\rm U}\sim P(\vec{\theta}^{*}_{\rm U}) \\
(z_s,z_L,\vec{\theta}^{*}_{\rm L})\sim P(z_s,z_L,\vec{\theta}^{*}_{\rm L}\mid{\rm SL}) \\
\vec{\beta}\sim P(\vec{\beta}\mid z_s,z_L,\vec{\theta}^{*}_{\rm L},{\rm SL})
}} \, ,
\end{split}
\end{equation}
$$

where the brackets denote the average over samples drawn from the specified conditional distributions.

## Decomposing the Joint Source–Lens Parameter Distribution

The joint distribution $P(z_s, z_L, \vec{\theta}^{*}_{\rm L}\mid{\rm SL})$ is decomposed using the chain rule to enable a hierarchical sampling strategy. It reads

$$ \begin{equation}
\begin{split}
P(z_s, z_L, \vec{\theta}^{*}_{\rm L}\mid{\rm SL})
&= P(z_s\mid{\rm SL})\,P(z_L\mid z_s,{\rm SL})\,P(\vec{\theta}^{*}_{\rm L}\mid z_L, z_s,{\rm SL}) \, .
\end{split}
\end{equation} $$

The first term $P(z_s\mid{\rm SL})$ represents the redshift distribution of lensed sources. The second term $P(z_L\mid z_s,{\rm SL})$ is the redshift distribution of lenses for a specific source redshift. The third term $P(\vec{\theta}^{*}_{\rm L}\mid z_L, z_s,{\rm SL})$ represents the distribution of lens structural parameters conditioned on the occurrence of strong lensing. Using Bayes' theorem, this final term can be further decomposed as

$$ \begin{equation}
\begin{split}
P(\vec{\theta}^{*}_{\rm L}\mid z_L, z_s,{\rm SL})
&= \frac{P({\rm SL}\mid z_L, z_s, \vec{\theta}^{*}_{\rm L}) \, P(\vec{\theta}^{*}_{\rm L} \mid z_L, z_s)}{P({\rm SL}\mid z_L, z_s)} \, .
\end{split}
\end{equation} $$

For EPL+Shear lenses, it is proportional to the [strong-lensing cross-section](#multi-image-lensing-cross-section) and the intrinsic priors

$$ \begin{equation}
\begin{split}
P(\vec{\theta}^{*}_{\rm L}\mid z_L, z_s,{\rm SL})
&\propto P({\rm SL}\mid z_L, z_s, \vec{\theta}^{*}_{\rm L}) \, P(\vec{\theta}^{*}_{\rm L} \mid z_L, z_s) \\
&\propto \sigma^{\rm MI}_{\rm EPL} \, P(\gamma_1, \gamma_2) \, P(\gamma) \, P(\phi_{\rm rot}) \, P(q \mid \sigma) \, P(\sigma \mid z_L) \, .
\end{split}
\end{equation} $$

This decomposition allows for a three-step sampling process: drawing the source redshift, then the lens redshift, and finally the lens properties using importance sampling to account for the cross-section weight. More details on constructing each factor are provided in the subsequent sections.


## Source Redshift Distributions of Lensed Events

The redshift distribution of strongly lensed sources differs from the intrinsic distribution because the probability of lensing increases with redshift. Using Bayes' theorem, the conditional distribution is

$$ \begin{equation}
\begin{split}
P(z_s\mid{\rm SL})
= \frac{P({\rm SL}\mid z_s)\,P(z_s)}{P({\rm SL})} \, ,
\end{split}
\end{equation} $$

where $P({\rm SL}\mid z_s)$ is the strong-lensing optical depth. Substituting the intrinsic redshift distribution for the unlensed sources, the expression becomes

$$ 
\begin{equation}
\begin{split}
P(z_s)
= \frac{1}{{\cal N}_{\rm U}}
\frac{R_{\rm U}(z_s)}{1+z_s}\,
\frac{dV_c}{dz_s} \, ,
\end{split}
\end{equation} $$
the lensed-source redshift distribution becomes
$$ \begin{equation}
\begin{split}
P(z_s\mid{\rm SL})
&= \frac{P({\rm SL}\mid z_s)}{P({\rm SL})}
\left[
\frac{1}{{\cal N}_{\rm U}}
\frac{R_{\rm U}(z_s)}{1+z_s}\,
\frac{dV_c}{dz_s}
\right] \\
&= \frac{1}{{\cal N}_{\rm L}}\,
P({\rm SL}\mid z_s)\,
\frac{R_{\rm U}(z_s)}{1+z_s}\,
\frac{dV_c}{dz_s} \, .
\end{split}
\end{equation} 
$$

Here $R_{\rm U}(z_s)$ is the intrinsic merger-rate density in the source frame, $dV_c/dz_s$ is the comoving volume element, and $(1+z_s)^{-1}$ accounts for cosmological time dilation when converting from source-frame to detector-frame time. The constant ${\cal N}_{\rm U}$ is the total intrinsic (unlensed) merger rate per unit detector-frame time.

Furthermore, the lensed-rate normalization, representing the total intrinsic strongly lensed merger rate per unit detector-frame time, follows

$$
\begin{equation}
\begin{split}
{\cal N}_{\rm L}
&= {\cal N}_{\rm U}\,P({\rm SL}) \\
&= {\cal N}_{\rm U}\int P({\rm SL}\mid z_s)\,P(z_s)\,dz_s \\
&= \int P({\rm SL}\mid z_s)\,
\frac{R_{\rm U}(z_s)}{1+z_s}\,
\frac{dV_c}{dz_s}\,dz_s \, .
\end{split}
\end{equation}
$$

The overall strong-lensing probability is therefore $P({\rm SL})={\cal N}_{\rm L}/{\cal N}_{\rm U}$, which is the fraction of intrinsically occurring mergers that are strongly lensed in the population model. In this work, the optical depth $P({\rm SL}\mid z_s)$ is computed for the adopted EPL+Shear lens population and is detailed in the next section.


## Optical Depth

The strong-lensing optical depth $P({\rm SL}\mid z_s)$ is the probability that a GW source at redshift $z_s$ is strongly lensed by an intervening galaxy population, resulting in multiple images. Geometrically, it represents the fraction of the sky covered by the effective strong-lensing cross-sections of the intervening galaxy population.  Mathematically, it is obtained by integrating the angular cross-section $\sigma_{\rm SL}$ over the lens parameters $\vec{\theta}^{*}_{\rm L}$ and redshift $z_L$

$$ 
\begin{equation}
\begin{split}
P({\rm SL}\mid z_s)
&= \int_{0}^{z_s}
\left[
\int_{\vec{\theta}^{*}_{\rm L}}
\frac{\sigma_{\rm SL}(z_L,z_s,\vec{\theta}^{*}_{\rm L})}{4\pi}\,
\frac{d^2N(z_L,\vec{\theta}^{*}_{\rm L})}{dV_c\,d\vec{\theta}^{*}_{\rm L}}\,
\frac{dV_c}{dz_L}\,
d\vec{\theta}^{*}_{\rm L}
\right] dz_L \, .
\end{split}
\end{equation}
$$

To simplify the calculation, we define the effective lensing kernel $\Phi_{\rm SL}(z_L,z_s)$ as the integrand marginalized over the lens properties $\vec{\theta}^{*}_{\rm L}$

$$ 
\begin{equation}
\begin{split}
\Phi_{\rm SL}(z_L,z_s)
&= \int_{\vec{\theta}^{*}_{\rm L}}
\frac{\sigma_{\rm SL}(z_L,z_s,\vec{\theta}^{*}_{\rm L})}{4\pi}\,
\frac{d^2N(z_L,\vec{\theta}^{*}_{\rm L})}{dV_c\,d\vec{\theta}^{*}_{\rm L}}\,
\frac{dV_c}{dz_L}\,
d\vec{\theta}^{*}_{\rm L} \, .
\end{split}
\end{equation} 
$$

The total optical depth is then simply the integral of this kernel along the line of sight, $P({\rm SL}\mid z_s) = \int_{0}^{z_s} \Phi_{\rm SL}(z_L,z_s)\,dz_L$.

For the EPL+Shear model, the lens population is described by the velocity-dispersion function $\phi(\sigma,z_L)$ and conditional priors for the shape parameters. The kernel expansion becomes

$$ 
\begin{equation}
\begin{split}
\Phi_{\rm SL}(z_L, z_s)
&= \int_{\sigma,q,\phi_{\rm rot},\gamma,\gamma_1,\gamma_2} P({\rm SL}|z_s,z_L,\sigma,q,\phi_{\rm rot},\gamma,\gamma_1,\gamma_2) \\
& \quad \times P(\gamma_1, \gamma_2) \, P(\gamma) \, P(\phi_{\rm rot}) \, P(q \mid \sigma) \\
& \quad \times \frac{d^2N(z_L, \sigma)}{dV_c \, d\sigma} \, \frac{dV_c}{dz_L} \, d\sigma \, dq \, d\phi_{\rm rot} \, d\gamma \, d\gamma_1 \, d\gamma_2\, .
\end{split}
\end{equation} 
$$

Numerically, this integral is approximated using Monte Carlo sampling. By drawing the velocity dispersion $\sigma$ from a uniform proposal distribution $P_o(\sigma)=1/\Delta\sigma$, the kernel is estimated as

$$ 
\begin{equation}
\begin{split}
\Phi_{\rm SL}(z_L,z_s)
= \Delta\sigma
\left\langle
\frac{\sigma^{\rm MI}_{\rm EPL}}{4\pi}\,
\phi(\sigma,z_L)\,
\frac{dV_c}{dz_L}
\right\rangle_{\substack{
\sigma\sim P_o(\sigma),\\
q\sim P(q\mid\sigma),\\
\phi_{\rm rot}\sim P(\phi_{\rm rot}),\\
\gamma\sim P(\gamma),\\
(\gamma_1,\gamma_2)\sim P(\gamma_1,\gamma_2)
}} \, ,
\end{split}
\end{equation} 
$$

where $\sigma^{\rm MI}_{\rm EPL}$ is the multi-image caustic cross-section. This kernel $\Phi_{\rm SL}(z_L,z_s)$ is central to the analysis, as it determines both the total optical depth and the probability distribution of lens redshifts $P(z_L\mid z_s,{\rm SL})$.


## Lens Redshift Distribution of Lensed Events

The distribution of lens redshifts for a fixed source redshift $z_s$ is directly proportional to the effective lensing kernel derived in the optical depth calculation. It is given by

$$ \begin{equation}
\begin{split}
P(z_L\mid z_s,{\rm SL})
= \frac{\Phi_{\rm SL}(z_L,z_s)}{{\cal N}_\Phi(z_s)} \, ,
\end{split}
\end{equation} $$

where ${\cal N}_\Phi(z_s)$ normalizes the distribution. In the [`ler`](https://ler.hemantaph.com/) package, this distribution is precomputed on a grid to allow for efficient inverse-transform sampling of $z_L$ given a sampled $z_s$.

For comparison, the intrinsic lens redshift distribution without conditioning on strong lensing is obtained by marginalizing over lens parameters,

$$ \begin{equation}
\begin{split}
P(z_L)
&= \frac{1}{{\cal N}_{z_L}}
\int_{\vec{\theta}^{*}_{\rm L}}
\frac{d^2N(z_L,\vec{\theta}^{*}_{\rm L})}{dV_c\, d\vec{\theta}^{*}_{\rm L}}\,
\frac{dV_c}{dz_L}\,
d\vec{\theta}^{*}_{\rm L} \\
&= \frac{\Delta\sigma}{{\cal N}_{z_L}}
\left\langle
\phi(\sigma,z_L)\,
\frac{dV_c}{dz_L}
\right\rangle_{\substack{
\sigma\sim P_o(\sigma),\\
q\sim P(q\mid\sigma),\\
\phi_{\rm rot}\sim P(\phi_{\rm rot}),\\
\gamma\sim P(\gamma),\\
(\gamma_1,\gamma_2)\sim P(\gamma_1,\gamma_2)
}} \, ,
\end{split}
\end{equation} 
$$

with normalization ${\cal N}_{z_L}$.


## Multi-Image Lensing Cross-Section

The strong-lensing condition is satisfied when a source lies within the caustic region of the lens, resulting in multiple images. The angular area of this region is the multi-image cross-section $\sigma^{\rm MI}_{\rm EPL}$. The probability of strong lensing for a specific configuration is the ratio of this cross-section to the total sky area

$$ \begin{equation}
\begin{split}
P({\rm SL}\mid z_L, z_s, \vec{\theta}^{*}_{\rm L})
= \frac{\sigma^{\rm MI}_{\rm EPL}}{4\pi} \, .
\end{split}
\end{equation} $$

Direct evaluation of $\sigma^{\rm MI}_{\rm EPL}$ by tracing caustics with [`lenstronomy`](https://lenstronomy.readthedocs.io) is computationally expensive. In [`ler`](https://ler.hemantaph.com/), $\sigma^{\rm MI}_{\rm EPL}$ is therefore precomputed on a grid in $(q,\phi_{\rm rot},\gamma,\gamma_1,\gamma_2)$ at unit Einstein radius $\theta_{\rm E}=1$, and interpolated during sampling. For arbitrary $(\sigma,z_L,z_s)$, the Einstein radius is

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
\sigma^{\rm MI}_{\rm EPL}(\theta_{\rm E}, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2)
&= \sigma^{\rm MI}_{\rm EPL}(1, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2) \\
&\quad\times \left({\rm intercept} + {\rm slope}\,\pi\,\theta_{\rm E}^2\right) \, ,
\end{split}
\end{equation}
$$

where ${\rm intercept}$ and ${\rm slope}$ are precomputed constants. This interpolation and rescaling scheme accelerates cross-section evaluation and, consequently, the computation of $\Phi_{\rm SL}(z_L,z_s)$ and the optical depth $P({\rm SL}\mid z_s)$.

## Lens Parameter Distributions of Lensed Events

For a fixed $(z_L,z_s)$, the distribution of lens parameters for strongly lensed events is obtained by weighting the intrinsic lens population by the multi-image cross-section. For the EPL+Shear model, it reads

$$
\begin{equation}
\begin{split}
&P(\sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2 \mid z_L, z_s, {\rm SL})\\
&\quad \propto \sigma^{\rm MI}_{\rm EPL} \, P(\gamma_1, \gamma_2) \, P(\gamma) \, P(\phi_{\rm rot}) \, P(q \mid \sigma) \, P(\sigma \mid z_L) \, .
\end{split}
\end{equation}
$$

Direct rejection sampling from the intrinsic distribution $P(\sigma\mid z_L)$ is highly inefficient because the lensing weight $\sigma^{\rm MI}_{\rm EPL}$ scales strongly with $\sigma$, favoring high-mass lenses that are rare in the intrinsic population. Furthermore, the lack of an analytical maximum for $\sigma^{\rm MI}_{\rm EPL}$ complicates standard rejection techniques.

To address this, [`ler`](https://ler.hemantaph.com/) employs importance sampling. A proposal distribution $P_o(\sigma)$, which is typically uniform over $[\sigma_{\min},\sigma_{\max}]$, is introduced to rewrite the target distribution as

$$
\begin{equation}
\begin{split}
&P(\sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2 \mid z_L, z_s, {\rm SL})\\
&\quad \propto \left[\frac{ \sigma^{\rm MI}_{\rm EPL}\, P(\sigma \mid z_L)}{P_o(\sigma)}\right] P(\gamma_1, \gamma_2) \, P(\gamma) \, P(\phi_{\rm rot}) \, P(q \mid \sigma) \, P_o(\sigma) \, .
\end{split}
\end{equation}
$$

In this scheme, $\sigma$ is drawn from $P_o(\sigma)$ and the other parameters from their intrinsic priors. Each sample set is then assigned an importance weight

$$ 
\begin{equation}
\begin{split}
w(\sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2;z_L,z_s)=\frac{\sigma^{\rm MI}_{\rm EPL}\,P(\sigma\mid z_L)}{P_o(\sigma)} \, .
\end{split}
\end{equation} 
$$

These weights are applied during the final event-rate integration to accurately represent the lensed population.



## Source Position Distribution and Image Properties



The source position $\vec{\beta}$ is sampled uniformly within the multi-image caustic region defined by the lens parameters. This ensures that every sampled configuration produces multiple images (typically two, three or four, ignoring the demagnified central image).

The lens equation is solved using [`lenstronomy`](https://lenstronomy.readthedocs.io) to determine image positions $\vec{\theta}_i$, magnifications $\mu_i$, and arrival-time delays $t_i$. These quantities are initially computed at $\theta_{\rm E}=1$ and then rescaled to physical units using the actual Einstein radius of the system.

Images are classified by their Morse index $n_i$ into Type-I (minima, $n_i=0$) and Type-II (saddle, $n_i=1/2$) images. Type-III (maxima) images are generally too faint to be detected and are excluded. Lensing modifies the observed GW parameters for each image $i$ as follows

$$
\begin{equation}
\begin{split}
\vec{\theta}_{{\rm GW},i} &\in \Bigg\{
\frac{d_L(z_s)}{\sqrt{\lvert \mu_i\rvert}},\;
\phi_c - n_i\pi,\;
t_c + t_i,\; \\
&\quad\quad {\rm RA} + \frac{\theta_{x,i}-\beta_x}{\cos({\rm Dec})},\;
{\rm Dec} + (\theta_{y,i}-\beta_y),\;\ldots
\Bigg\}
\end{split} \, ,
\end{equation}
$$

while detector-frame component masses and spins are unchanged from the unlensed values. The sky-location offsets associated with $\vec{\theta}_i-\vec{\beta}$ are typically negligible for current GW detector networks and have no impact on SNR. Detectability is therefore evaluated by computing the SNR for each image using $\vec{\theta}_{{\rm GW},i}$, as described in the next section.

## Detection Probability of Lensed Events

For a specific lensed configuration, the detection probability is determined by the detectability of its images. The signal-to-noise ratio (SNR), $\rho(\vec{\theta}_{{\rm GW},i})$, is computed for each image using the modified GW parameters. An event is considered detected if at least two images meet the detection criterion. The probability is therefore a binary condition

$$
\begin{equation}
\begin{split}
&P({\rm obs} \mid \vec{\theta}^*_{\rm U}, z_s, \vec{\theta_{\rm L}}, \vec{\beta}, {\rm SL})\\
&\quad\quad = \Theta\left[\left(\sum_i P_{\rm det}(\vec{\theta}_{{\rm GW},i}, \rho_{\rm th}) \right) - 2\right] \, ,
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

and in `ler`, it is estimated using the [`gwsnr`](https://gwsnr.hemantaph.com) package. 


## Complete Expression for Lensed Event Rate

Combining all components, the final expression for the annual rate of detectable strongly lensed GW events is

$$ 
\begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}
=
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

This equation represents the average detection probability over the population of strongly lensed events, weighted by their respective occurrence probabilities.

The corresponding parameter priors and their adopted forms are summarized below.

## Lens Parameter Priors

The priors used for the redshift and lens parameter distributions are summarized below.

### Table 1: Redshifts conditioned on strong lensing

| Parameter | Unit | Prior Distribution | Range [Min, Max] | Description |
| :--- | :--- | :--- | :--- | :--- |
| $z_s$ | - | $\propto \frac{P({\rm SL} \mid z_s) R_{\rm U}(z_s)}{1+z_s}\,\frac{dV_c}{dz_s}$ | [0.0, 10.0] | Source redshift weighted by optical depth |
| $z_L$ | - | $\frac{\Phi(z_L,z_s)}{{\cal N}_\Phi(z_s)}$ | [0.0, $z_s$] | Lens redshift weighted by effective kernel |

### Table 2: Lens intrinsic priors used with cross-section weighting

| Parameter | Unit | Prior Distribution | Range [Min, Max] | Description |
| :--- | :--- | :--- | :--- | :--- |
| $\sigma$ | km/s | $\frac{d^2N(z_L, \sigma)}{dV_c \, d\sigma}$ | [100, 400] | Velocity dispersion (Oguri et al. 2018) |
| $q$ | - | Rayleigh distribution | [0.2, 1.0] | Axis ratio (Collett et al. 2015) |
| $\phi_{\rm rot}$ | rad | Uniform | [0, $\pi$] | Lens orientation |
| $\gamma$ | - | Normal | - | Density profile slope (Mean: 2.0, Std: 0.1) |
| $\gamma_1, \gamma_2$ | - | Normal | - | External shear (Mean: 0.0, Std: 0.05) |

### Table 3: Source position prior



| Parameter | Unit | Prior Distribution | Description |
| :--- | :--- | :--- | :--- |
| $\beta_x,$ | $\theta_{\rm E}$ | Uniform within<br> multi-image caustic | x-component of source position from lens center |
| $\beta_y$ | $\theta_{\rm E}$ | Uniform within<br> multi-image caustic | y-component of source position from lens center |