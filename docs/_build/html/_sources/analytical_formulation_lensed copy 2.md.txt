# Analytical Formulation for Lensed Gravitational Wave Event Rates

## Table of Contents

1. [Introduction](#introduction)
2. [Parameter-Marginalized Event Rate](#parameter-marginalized-event-rate)

## Introduction

The annual rate of detectable strongly lensed gravitational-wave (GW) events, $\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}$, is the expected number of observed lensed compact-binary mergers per year for a given detector or detector network. This rate is obtained by scaling the total intrinsic merger rate, $\frac{\Delta N_{\rm U}}{\Delta t}$, by the probability that an intrinsically occurring merger is both strongly lensed and detectable. Denoting this joint probability by $P({\rm SL, obs})$, one may write

$$ 
\begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}
= \frac{\Delta N_{\rm U}}{\Delta t}\, P({\rm SL, obs}) \, .
\end{split}
\end{equation} 
$$

Using the product rule, the joint probability can be decomposed as $P({\rm SL, obs}) = P({\rm SL})\,P({\rm obs}\mid{\rm SL})$, where $P({\rm SL})$ is the strong-lensing probability in the intrinsic population and $P({\rm obs}\mid{\rm SL})$ is the conditional detection probability for an event known to be strongly lensed. Substitution gives
$$ \begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}
= \frac{\Delta N_{\rm U}}{\Delta t}\, P({\rm SL})\, P({\rm obs}\mid{\rm SL})
= \frac{\Delta N_{\rm L}}{\Delta t}\, P({\rm obs}\mid{\rm SL}) \, .
\end{split}
\end{equation} $$

Here $\frac{\Delta N_{\rm L}}{\Delta t}$ is the total rate of strongly lensed mergers irrespective of detectability. For notational convenience, it is denoted by ${\cal N}_{\rm L}$ and defined by
$$ \begin{equation}
\begin{split}
\frac{\Delta N_{\rm L}}{\Delta t}
= \frac{\Delta N_{\rm U}}{\Delta t}\, P({\rm SL})
\equiv {\cal N}_{\rm L} \, .
\end{split}
\end{equation} $$

## Parameter-Marginalized Event Rate

The conditional detection probability $P({\rm obs}\mid{\rm SL})$ is obtained by marginalizing over the unlensed GW source parameters $\vec{\theta}_{\rm U}$, the lens parameters $\vec{\theta_{\rm L}}$, and the source position on the lens-source plane, $\vec{\beta}$. The definition of $\vec{\theta}_{\rm U}$ follows the unlensed-event formulation in [Analytical Formulation for Gravitational Wave Event Rates](https://ler.hemantaph.com/analytical_formulation_unlensed.html). For an Elliptical Power Law lens with external shear (EPL+Shear), the lens-parameter vector is

$$
\begin{equation}
\begin{split}
\vec{\theta_{\rm L}}\in \{z_L,\sigma,q,\phi_{\rm rot},\gamma,\gamma_1,\gamma_2\} \,,
\end{split}
\end{equation}
$$

where $z_L$ is the lens redshift, $\sigma$ is the velocity dispersion, $q$ is the projected axis ratio, $\phi_{\rm rot}$ is the lens orientation, $\gamma$ is the density-profile slope, and $(\gamma_1,\gamma_2)$ are the external shear components. The source position is written as $\vec{\beta}\in\{\beta_x,\beta_y\}$ and denotes the angular offset of the source from the lens center on the source plane. The sampling strategy for lens-paramters and source positons are dicussed in the following sections. The priors adopted for these parameters are summarized in [Table 2](#) and [Table 3](#).

With this notation, the annual rate of detectable strongly lensed events can be expressed as a parameter-marginalized expectation over strongly lensed configurations,

$$
\begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}
&= \frac{\Delta N_{\rm L}}{\Delta t}
\int_{\vec{\theta_{\rm U}}}\int_{\vec{\theta_{\rm L}}}\int_{\vec{\beta}} P({\rm obs} \mid \vec{\theta_{\rm U}}, \vec{\theta_{\rm L}}, \vec{\beta}, {\rm SL}) \\
&\quad\quad\quad\quad\quad\quad\quad\quad \times P(\vec{\theta_{\rm U}}, \vec{\theta_{\rm L}}, \vec{\beta} \mid {\rm SL}) d\vec{\theta_{\rm U}}\, d\vec{\theta_{\rm L}}\, d\vec{\beta} \,,
\end{split}
\end{equation}
$$

where $P({\rm obs}\mid \vec{\theta}_{\rm U},\vec{\theta_{\rm L}},\vec{\beta},{\rm SL})$ is the detection probability for a specific source–lens configuration, and $P(\vec{\theta}_{\rm U},\vec{\theta_{\rm L}},\vec{\beta}\mid{\rm SL})$ is the joint distribution of parameters conditioned on the event being strongly lensed.

It is convenient to separate the source redshift from the remaining intrinsic GW parameters by writing $\vec{\theta}_{\rm U}=(z_s,\vec{\theta}^{*}_{\rm U})$, where $\vec{\theta}^{*}_{\rm U}$ denotes the unlensed source parameters excluding $z_s$. Similarly, we write $\vec{\theta_{\rm L}}=(z_L,\vec{\theta}^{*}_{\rm L})$, where $\vec{\theta}^{*}_{\rm L}$ denotes the lens parameters excluding $z_L$.

Assuming that $\vec{\theta}^{*}_{\rm U}$ is independent of the detailed lensing configuration once $z_s$ is fixed and the event is conditioned to be strongly lensed, the joint distribution factorizes as

$$ 
\begin{equation}
\begin{split}
P(\vec{\theta}_{\rm U},\vec{\theta_{\rm L}},\vec{\beta}\mid{\rm SL})
&= P(\vec{\theta}^{*}_{\rm U})\,
P(\vec{\beta}\mid z_s, z_L,\vec{\theta}^{*}_{\rm L},{\rm SL})\\
&\quad\times
P(z_s, z_L,\vec{\theta}^{*}_{\rm L}\mid{\rm SL}) \, .
\end{split}
\end{equation} 
$$

Substituting this factorization into the rate expression yields

$$
\begin{equation}
\begin{split}
&\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}\\
&\quad\quad= \frac{\Delta N_{\rm L}}{\Delta t}
\int_{\vec{\theta}_{\rm U}}\int_{\vec{\theta_{\rm L}}}\int_{\vec{\beta}}
 P({\rm obs} \mid \vec{\theta}^*_{\rm U}, z_s, z_L,\vec{\theta}^{*}_{\rm L}, \vec{\beta}, {\rm SL}) \\
&\quad\quad\quad\quad\quad\quad\quad\quad \times P(\vec{\theta}^*_{\rm U})\, P(\vec{\beta} \mid z_s, z_L,\vec{\theta}^{*}_{\rm L}, {\rm SL}) \\
&\quad\quad\quad\quad\quad\quad\quad\quad \times P(z_s, z_L,\vec{\theta}^{*}_{\rm L} \mid {\rm SL}) \, d\vec{\theta}_{\rm U}\, d\vec{\theta_{\rm L}}\, d\vec{\beta} \,,
\end{split}
\end{equation}
$$

In practice, the multidimensional integral is evaluated with Monte Carlo sampling over the strongly lensed population, so that the observable lensed rate can be written as

$$ 
\begin{equation}
\begin{split}
&\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}\\
&\quad\quad = {\cal N}_{\rm L}
\left\langle
P({\rm obs}\mid \vec{\theta}^{*}_{\rm U},z_s,z_L,\vec{\theta}^{*}_{\rm L},\vec{\beta},{\rm SL})
\right\rangle_{\substack{
\vec{\theta}^{*}_{\rm U}\sim P(\vec{\theta}^{*}_{\rm U}) \\
(z_s,z_L,\vec{\theta}^{*}_{\rm L})\sim P(z_s,z_L,\vec{\theta}^{*}_{\rm L}\mid{\rm SL}) \\
\vec{\beta}\sim P(\vec{\beta}\mid z_s,z_L,\vec{\theta}^{*}_{\rm L},{\rm SL})
}} \,,
\end{split}
\end{equation} 
$$

where $\langle\cdot\rangle$ denotes an average over samples drawn from the indicated distributions.

## Decomposing the Joint Source–Lens Parameter Distribution

The joint distribution of source redshift and lens parameters conditioned on strong lensing can then be written using the chain rule as

$$
\begin{equation}
\begin{split}
P(z_s, z_L, \vec{\theta}^{*}_{\rm L}\mid{\rm SL})
&= P(z_s\mid{\rm SL})\,P(z_L\mid z_s,{\rm SL})\\
&\quad\times P(\vec{\theta}^{*}_{\rm L}\mid z_L, z_s,{\rm SL}) \,.
\end{split}
\end{equation} 
$$

This factorization is useful because each term has a direct physical meaning and a corresponding sampling strategy. The first factor, $P(z_s\mid{\rm SL})$, is the lensed-source redshift distribution obtained by weighting the intrinsic source redshift distribution $P(z_s)$ by the strong-lensing optical depth $P({\rm SL}\mid z_s)$. The second factor, $P(z_L\mid z_s,{\rm SL})$, describes the lens-redshift distribution for a given lensed source redshift and is determined by the effective lensing kernel $\Phi_{\rm SL}(z_L,z_s)$ (see next section). 
The third factor, $P(\vec{\theta}^{*}_{\rm L}\mid z_L,z_s,{\rm SL})$, is the conditional distribution of lens structural parameters at fixed $(z_L,z_s)$ under the strong-lensing condition, and it is proportional to the strong-lensing cross-section for that configuration times the intrinsic lens-parameter priors. For EPL+Shear lenses, the conditional distribution of $\vec{\theta}^{*}_{\rm L}$ can be expressed as

$$
\begin{equation}
\begin{split}
&P(\vec{\theta}^{*}_{\rm L}\mid z_L, z_s,{\rm SL})\\
&\quad \propto P({\rm SL}\mid z_L, z_s, \vec{\theta}^{*}_{\rm L}) \\
&\quad\quad\times P(\vec{\theta_{\rm L}^*} \mid z_L, z_s) \\
&\quad \propto P({\rm SL} \mid \sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2, z_L, z_s) \\
&\quad\quad\times P(\gamma_1, \gamma_2) \, P(\gamma) \, P(\phi_{\rm rot}) \, P(q \mid \sigma) \, P(\sigma \mid z_L) \, ,
\end{split}
\end{equation}
$$

where $P({\rm SL}\mid z_L, z_s, \vec{\theta}^{*}_{\rm L})$ is set by the lensing cross-section (more details on cross-section in the following section) at fixed $(z_L,z_s,\vec{\theta}^{*}_{\rm L})$. Normalization is implicit in the proportionality and is handled numerically in the Monte Carlo procedure.

With this decomposition, sampling from $P(z_s, z_L, \vec{\theta}^{*}_{\rm L}\mid{\rm SL})$ proceeds naturally in three stages. One first draws $z_s$ from $P(z_s\mid{\rm SL})$, then draws $z_L$ from $P(z_L\mid z_s,{\rm SL})$ at that $z_s$, and finally draws $\vec{\theta}^{*}_{\rm L}$ from its intrinsic priors while accounting for the strong-lensing cross-section through weighting or importance sampling at the chosen $(z_s,z_L)$. The construction of each factor is detailed in the sections below.


## Source Redshift Distributions of Lensed Events

The redshift distribution of strongly lensed sources, $P(z_s\mid{\rm SL})$, is not the same as the intrinsic source redshift distribution $P(z_s)$ because the strong-lensing probability varies with redshift. By Bayes’ theorem,

$$ 
\begin{equation}
\begin{split}
P(z_s\mid{\rm SL})
= \frac{P({\rm SL}\mid z_s)\,P(z_s)}{P({\rm SL})} \, ,
\end{split}
\end{equation}
$$

where $P({\rm SL}\mid z_s)$ is the strong-lensing optical depth for a source at redshift $z_s$, and $P({\rm SL})$ is the population-averaged strong-lensing probability.

Using the intrinsic redshift distribution implied by the unlensed-event rate density (see [Analytical Formulation for Gravitational Wave Event Rates](https://ler.hemantaph.com/analytical_formulation_unlensed.html#redshift-distribution-and-intrinsic-merger-rates)),

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

Here $R_{\rm U}(z_s)$ is the intrinsic merger-rate density in the source frame, $dV_c/dz_s$ is the comoving volume element, and $(1+z_s)^{-1}$ accounts for cosmological time dilation when converting from source-frame to detector-frame time. The normalization constant ${\cal N}_{\rm U}$ is the total intrinsic (unlensed) merger rate per unit detector-frame time, while ${\cal N}_{\rm L}$ is the total intrinsic strongly lensed merger rate per unit detector-frame time.

The lensed-rate normalization follows directly from the requirement $\int P(z_s\mid{\rm SL})\,dz_s=1$:

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

The overall strong-lensing probability is therefore $P({\rm SL})={\cal N}_{\rm L}/{\cal N}_{\rm U}$, which is the fraction of intrinsically occurring mergers that are strongly lensed in the population model. In this work, the optical depth $P({\rm SL}\mid z_s)$ is computed for the adopted EPL+Shear lens population and is described in the next section.

## Optical Depth

The strong-lensing optical depth $P({\rm SL}\mid z_s)$ is the probability that a GW source at redshift $z_s$ is strongly lensed by an intervening galaxy population, resulting in multiple images. It can be written as the sky-fraction covered by the strong-lensing angular cross-sections of the lens population between the observer and the source, and it reads

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
\right] dz_L \, ,
\end{split}
\end{equation}
$$

where $\sigma_{\rm SL}(z_L,z_s,\vec{\theta}^{*}_{\rm L})$ is the angular cross-section for a single lens with parameters $(z_L,\vec{\theta}^{*}_{\rm L})$ and it is detailed in the following section.

Defining the effective lensing kernel
$$ \begin{equation}
\begin{split}
\Phi_{\rm SL}(z_L,z_s)
&= \int_{\vec{\theta}^{*}_{\rm L}}
\frac{\sigma_{\rm SL}(z_L,z_s,\vec{\theta}^{*}_{\rm L})}{4\pi}\,
\frac{d^2N(z_L,\vec{\theta}^{*}_{\rm L})}{dV_c\,d\vec{\theta}^{*}_{\rm L}}\,
\frac{dV_c}{dz_L}\,
d\vec{\theta}^{*}_{\rm L} \, ,
\end{split}
\end{equation} $$
the optical depth becomes
$$ \begin{equation}
\begin{split}
P({\rm SL}\mid z_s) = \int_{0}^{z_s} \Phi_{\rm SL}(z_L,z_s)\,dz_L \, .
\end{split}
\end{equation} $$

For EPL+Shear lenses with $\vec{\theta}^{*}_{\rm L}=\{\sigma,q,\phi_{\rm rot},\gamma,\gamma_1,\gamma_2\}$, we use the velocity-dispersion function $\phi(\sigma,z_L)=d^2N(z_L,\sigma)/(dV_c\,d\sigma)$ and priors for the remaining parameters, giving

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

The integration is numerically evaluated using Monte Carlo sampling over the lens parameter space. While lens parameters are sampled from their intrinsic ppopulation, $\sigma$ is drawn from a proposal distribution $P_o(\sigma)=\frac{1}{\sigma_{\max}-\sigma_{\min}}=\frac{1}{\Delta \sigma}$ for convenience. Thus, the effective lensing kernel can be expressed as

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

where $\sigma^{\rm MI}_{\rm EPL}$ is the multi-image caustic cross-section for EPL+Shear lenses (see following section).

The kernel $\Phi_{\rm SL}(z_L,z_s)$, which also represents the effective lensing cross-section for lenses at redshift $z_L$, marginalized over the lens parameters $\vec{\theta_{\rm L}^*}$ and normalized by the total solid angle of the sky, is used both to evaluate $P({\rm SL}\mid z_s)$ and to sample $z_L$ from $P(z_L\mid z_s,{\rm SL})$ as described in the next section.


## Lens Redshift Distribution of Lensed Events

For sources at fixed redshift $z_s$ that are conditioned to be strongly lensed, the lens-redshift distribution is proportional to the effective lensing kernel $\Phi_{\rm SL}(z_L,z_s)$, giving

$$ 
\begin{equation}
\begin{split}
P(z_L\mid z_s,{\rm SL})
&= \frac{1}{{\cal N}_\Phi(z_s)}\, \Phi_{\rm SL}(z_L,z_s) \, ,
\end{split}
\end{equation} 
$$

where ${\cal N}_\Phi(z_s)=\int_{0}^{z_s}\Phi_{\rm SL}(z_L,z_s)\,dz_L$ is the normalization constant at fixed $z_s$.

In [`ler`](https://ler.hemantaph.com/), $\Phi_{\rm SL}(z_L,z_s)$ and the corresponding $P(z_L\mid z_s,{\rm SL})$ are precomputed on a $(z_L,z_s)$ grid for the adopted EPL+Shear lens population. Values at arbitrary $(z_L,z_s)$ are obtained by two-dimensional interpolation, and sampling at fixed $z_s$ is performed by inverse-transform sampling using the cumulative distribution function of $P(z_L\mid z_s,{\rm SL})$.

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
}} \, .
\end{split}
\end{equation} 
$$

with normalization ${\cal N}_{z_L}$.


## Multi-Image Lensing Cross-Section

In the lensed-event sampling, the strong-lensing condition is enforced by requiring that the lens-source configuration produces multiple images. Configurations that do not satisfy the multi-imaging condition are rejected and resampled.

For EPL+Shear lenses, multi-imaging occurs when the lens and source are aligned and source position lies inside the relevant caustic region on the source plane. The corresponding angular area defines the strong-lensing cross-section, denoted $\sigma^{\rm MI}_{\rm EPL}$. For a given lens configuration $(z_L,z_s,\vec{\theta}^{*}_{\rm L})$, the probability that a randomly oriented line of sight yields strong lensing is the cross-section normalized by the full sky,

$$ 
\begin{equation}
\begin{split}
P({\rm SL}\mid z_L, z_s, \vec{\theta}^{*}_{\rm L})
= \frac{\sigma^{\rm MI}_{\rm EPL}}{4\pi} \, .
\end{split}
\end{equation}
$$

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

At fixed $(z_L,z_s)$, the conditional distribution of EPL+Shear lens parameters under the strong-lensing condition is obtained by weighting the intrinsic lens population by the multi-image cross-section, and it reads

$$
\begin{equation}
\begin{split}
&P(\sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2 \mid z_L, z_s, {\rm SL})\\
&\quad \propto P({\rm SL} \mid \sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2, z_L, z_s) \\
&\quad\quad \times P(\gamma_1, \gamma_2) \, P(\gamma) \, P(\phi_{\rm rot}) \, P(q \mid \sigma) \, P(\sigma \mid z_L) \\
&\quad \propto \sigma^{\rm MI}_{\rm EPL} \, P(\gamma_1, \gamma_2) \, P(\gamma) \, P(\phi_{\rm rot}) \, P(q \mid \sigma) \, P(\sigma \mid z_L) \,.
\end{split}
\end{equation}
$$

Direct rejection sampling from $P(\sigma\mid z_L)$ is inefficient because the strong-lensing weight $\sigma^{\rm MI}_{\rm EPL}$ favours high-$\sigma$ lenses, while $P(\sigma\mid z_L)$ typically favours low-$\sigma$ lenses, leading to a low acceptance rate. Moreover, maximum value of $\sigma^{\rm DC}_{\rm EPL}$ is not known analytically for EPL+Shear lenses, making it difficult to use the standard rejection sampling method. In [`ler`](https://ler.hemantaph.com/), sampling is therefore performed with importance sampling by introducing a proposal distribution $P_o(\sigma)$ for the velocity dispersion. Writing

$$
\begin{equation}
\begin{split}
&P(\sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2 \mid z_L, z_s, {\rm SL})\\
&\quad \propto \left[\frac{ \sigma^{\rm MI}_{\rm EPL}\, P(\sigma \mid z_L)}{P_o(\sigma)}\right] \\
&\quad\quad \times P(\gamma_1, \gamma_2) \, P(\gamma) \, P(\phi_{\rm rot}) \, P(q \mid \sigma) \, P_o(\sigma) \,,
\end{split}
\end{equation}
$$

one draws $\sigma\sim P_o(\sigma)$ and samples $(q,\phi_{\rm rot},\gamma,\gamma_1,\gamma_2)$ from their intrinsic priors, then assigns each configuration the importance weight

$$ 
\begin{equation}
\begin{split}
w(\sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2;z_L,z_s)=\frac{\sigma^{\rm MI}_{\rm EPL}\,P(\sigma\mid z_L)}{P_o(\sigma)} \, .
\end{split}
\end{equation} $$
A convenient choice is a uniform proposal on $[\sigma_{\min},\sigma_{\max}]$, for which $P_o(\sigma)=1/(\sigma_{\max}-\sigma_{\min})$. The weighted samples are then used to estimate expectations over $P(\vec{\theta}^{*}_{\rm L}\mid z_L,z_s,{\rm SL})$ in the event-rate calculation.


## Source Position Distribution and Image Properties

Source-plane coordinates are expressed in units of the Einstein radius by setting $\theta_{\rm E}=1$ when sampling $\vec{\beta}$ and solving the lens equation. All angular and time-delay quantities are subsequently rescaled to physical units using the actual $\theta_{\rm E}(\sigma,z_L,z_s)$ of each lens configuration.

For a given set of lens-shape parameters $(q,\phi_{\rm rot},\gamma,\gamma_1,\gamma_2)$, the source position $\vec{\beta}\in\{\beta_x,\beta_y\}$ is drawn uniformly over the multi-imaging caustic region on the source plane. This guarantees that the configuration produces multiple images, typically two or four in the EPL+Shear model, with the central demagnified image ignored in this work.

For each sampled configuration, the lens equation is solved with [`lenstronomy`](https://lenstronomy.readthedocs.io) to obtain image positions $\vec{\theta}_i\in\{\theta_{x,i},\theta_{y,i}\}$, magnifications $\mu_i$, and arrival-time delays $t_i$ for each image $i$. The dimensionless quantities computed at $\theta_{\rm E}=1$ are rescaled to the physical angular and temporal units using the corresponding $\theta_{\rm E}$.

Image parity is determined from the Jacobian of the lens mapping and summarized by the Morse index $n_i$, with Type-I (minimum) images having $n_i=0$, Type-II (saddle) images having $n_i=1/2$, and Type-III (maximum) images having $n_i=1$. Only Type-I and Type-II images are retained. Images are ordered by arrival time, $t_1<t_2<t_3<t_4$, and up to four images are considered.

Lensing modifies the inferred GW parameters for each image through magnification, phase shifts, and arrival-time delays. The image-specific parameter set is updated as

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
\end{split} \notag
\end{equation}
$$

while detector-frame component masses and spins are unchanged from the unlensed values. The sky-location offsets associated with $\vec{\theta}_i-\vec{\beta}$ are typically negligible for current GW detector networks and have no impact on SNR. Detectability is therefore evaluated by computing the SNR for each image using $\vec{\theta}_{{\rm GW},i}$, as described in the next section.


## Detection Probability of Lensed Events

For a given lensed configuration specified by $(\vec{\theta}^{*}_{\rm U}, z_s, \vec{\theta}_{\rm L}, \vec{\beta})$ under the strong-lensing condition, the detection probability $P({\rm obs}\mid \vec{\theta}^{*}_{\rm U}, z_s, \vec{\theta}_{\rm L}, \vec{\beta}, {\rm SL})$ is determined by the detectability of the individual lensed images. For each image $i$, the signal-to-noise ratio $\rho(\vec{\theta}_{{\rm GW},i})$ is evaluated using the lensed GW parameters $\vec{\theta}_{{\rm GW},i}$, which incorporate magnification, arrival-time delay, and phase shift.

An event is counted as observable when at least two images satisfy the detection criterion, i.e. having SNR above a chosen threshold $\rho_{\rm th}$. Thus, the detection probability is given by

$$
\begin{equation}
\begin{split}
&P({\rm obs} \mid \vec{\theta}^*_{\rm U}, z_s, \vec{\theta_{\rm L}}, \vec{\beta}, {\rm SL}) \\
&\quad = \left\{ 
  \begin{array}{ c l }
    1 & \sum_i^{\rm images} \Theta[\rho(\vec{\theta}_{{\rm GW},i}) - \rho_{th}]\ge 2 \\
    0 & {\rm otherwise}
  \end{array}
\right. \,,
\end{split}
\end{equation}
$$

where $\Theta$ is the Heaviside step function, giving 1 when its argument is positive and 0 otherwise.


## Complete Expression for Lensed Event Rate

The annual rate of detectable strongly lensed GW events is written as the total intrinsic strong-lensing rate ${\cal N}_{\rm L}$ multiplied by the mean detection probability over the strongly lensed population:
$$ \begin{equation}
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
\end{equation} $$

The distributions appearing in the average define the Monte Carlo sampling procedure. The corresponding parameter priors and their adopted forms are summarized below.

## Lens Parameter Priors

Redshifts conditioned on strong lensing:

| Parameter | Unit | Prior <br>Distribution | Range <br>[Min, Max] | Description |
| :--- | :--- | :--- | :--- | :--- |
| $z_s$ | - | $\propto \frac{P({\rm SL} \mid z_s) R_{\rm U}(z_s)}{1+z_s}\,\frac{dV_c}{dz_s}$<br> Intrinsic merger rate<br> weighted by SL<br> optical depth | [0.0, 10.0] | Source redshift |
| $z_L$ | - | $\frac{\Phi(z_L,z_s)}{{\cal N}_\Phi(z_s)}$<br> Integrated number<br> density weighted<br> by SL cross-section | [0.0, $z_s$] | Lens redshift |

Lens intrinsic priors used with cross-section weighting:

| Parameter | Unit | Prior <br>Distribution | Range <br>[Min, Max] | Description |
| :--- | :--- | :--- | :--- | :--- |
| $\sigma$ | km/s | $\frac{d^2N(z_L, \sigma)}{dV_c \, d\sigma}$<br> number density | [100, 400] | Velocity dispersion<br> Oguri et al. 2018 |
| $q$ | - | Rayleigh<br> distribution | [0.2, 1.0] | Axis ratio<br> Collett et al. 2015 |
| $\phi_{\rm rot}$ | rad | Uniform | [0, $\pi$] | Lens orientation |
| $\gamma$ | - | Normal | - | Density profile slope<br> (Mean: 2.0, Std: 0.1)<br> Sonnenfeld et al. 2024 |
| $\gamma_1, \gamma_2$ | - | Normal | - | External shear components<br> (Mean: 0.0, Std: 0.05)<br> Collet et al. 2015 |

The intrinsic lens-parameter priors above are combined with the multi-image cross-section through importance weighting to obtain samples from $P(\vec{\theta}^{*}_{\rm L}\mid z_L,z_s,{\rm SL})$, as described in the previous section.
