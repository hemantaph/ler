# Analytical Formulation for Lensed Gravitational Wave Event Rates

## Table of Contents

1. [Introduction](#introduction)
2. [Parameter-Marginalized Event Rate](#parameter-marginalized-event-rate)

## Introduction

The annual rate of detectable strongly lensed gravitational-wave (GW) events, $\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}$, gives the expected number of observed lensed mergers per year for a given detector or detector network. It is obtained by weighting the total intrinsic compact-binary merger rate, $\frac{\Delta N_{\rm U}}{\Delta t}$, by the joint probability that a merger is strongly lensed and detected, $P({\rm SL, obs})$, such that

$$
\begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}
= \frac{\Delta N_{\rm U}}{\Delta t} \, P({\rm SL, obs}) \,,
\end{split}
\end{equation}
$$

where $P({\rm SL, obs})$ can be written as the product of the strong-lensing probability $P({\rm SL})$ and the conditional detection probability given strong lensing, $P({\rm obs \mid SL})$. It follows that

$$
\begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}
= \frac{\Delta N_{\rm U}}{\Delta t} \, P({\rm SL}) \, P({\rm obs \mid SL})
= \frac{\Delta N_{\rm L}}{\Delta t} \, P({\rm obs \mid SL}) \,,
\end{split}
\end{equation}
$$

with $\frac{\Delta N_{\rm L}}{\Delta t}=\frac{\Delta N_{\rm U}}{\Delta t}\,P({\rm SL})={\cal N}_{\rm L}$ denoting the total rate of strongly lensed mergers, irrespective of detectability.

## Parameter-Marginalized Event Rate

The conditional probability $P({\rm obs \mid SL})$ is evaluated by marginalizing over the relevant GW source parameters $\vec{\theta}$, the lens parameters $\vec{\theta_{\rm L}}$, and the source position relative to the lens $\beta$. For $\vec{\theta}$ see example [Analytical Formulation for Gravitational Wave Event Rates](https://ler.hemantaph.com/analytical_formulation_unlensed.html). For Elliptical Power Law (EPL) $+$ external shear lens models, these sets of parameters $\vec{\theta_{\rm L}}$ are defined as [Analytical Formulation for Gravitational Wave Event Rates](file:///Users/phurailatpamhemantakumar/phd/mypackages/ler/docs/_build/html/analytical_formulation_unlensed.html)

$$
\vec{\theta_{\rm L}}\in \{z_L, \sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2\},
$$

where $z_L$ is the lens redshift, $\sigma$ is the velocity dispersion, $q$ is the projected axis ratio, $\phi _{\rm rot}$ is the lens orientation, and $\gamma$ is the denisty profile slope, with $(\gamma_1,\gamma_2)$ describing external shear. The prior distributions used in this analysis are summarized in [Table 1](#). 

The source position $\vec{\beta} \in\{\beta_x,\beta_y\}$ represents the angular offset between the source and the lens center in the source plane. More about the $\beta$ sampling given in [section](#).

Making the parameter dependence explicit, the observed lensed event rate can be written as

$$
\begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}
= \frac{\Delta N_{\rm L}}{\Delta t}
\int_{\vec{\theta}}\int_{\vec{\theta_{\rm L}}}\int_{\vec{\beta}}
& P({\rm obs} \mid \vec{\theta}, \vec{\theta_{\rm L}}, \vec{\beta}, {\rm SL})\,
P(\vec{\theta}, \vec{\theta_{\rm L}}, \vec{\beta} \mid {\rm SL}) d\vec{\theta}\, d\vec{\theta_{\rm L}}\, d\vec{\beta} \,,
\end{split}
\end{equation}
$$

where $P({\rm obs} \mid \vec{\theta}, \vec{\theta_{\rm L}}, \vec{\beta}, {\rm SL})$ is the detection probability for a specific source–lens configuration, and $P(\vec{\theta}, \vec{\theta_{\rm L}}, \vec{\beta} \mid {\rm SL})$ is the joint distribution of parameters conditioned on strong lensing. Assuming that the intrinsic GW parameters $\vec{\theta^*}$ (* indicating without $z_s$) are independent of the detailed lensing configuration (and the strong lensing condition) once $z_s$ is specified and the event is known to be lensed, this joint distribution can be decomposed using the chain rule as

$$
\begin{equation}
\begin{split}
P(\vec{\theta}, \vec{\theta_{\rm L}}, \vec{\beta} \mid {\rm SL}) &= P(\vec{\theta^*}) \, P(\vec{\beta} \mid z_s, \vec{\theta_{\rm L}}, {\rm SL}) \, P(z_s, \vec{\theta_{\rm L}} \mid {\rm SL}) \,.
\end{split}
\end{equation}
$$

Numerically, the integral is evaluated using Monte Carlo sampling over the full parameter space. Substituting this expression back into the event-rate equation gives

$$
\begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}
&= \frac{\Delta N_{\rm L}}{\Delta t}
\int_{\vec{\theta^*}}\int_{\vec{\theta_{\rm L}}}\int_{\vec{\beta}}
 P({\rm obs} \mid \vec{\theta^*}, z_s, \vec{\theta_{\rm L}}, \vec{\beta}, {\rm SL})\,
P(\vec{\theta^*}) \, P(\vec{\beta} \mid z_s, \vec{\theta_{\rm L}}, {\rm SL}) \\
&\quad\quad\quad\quad\quad\quad\quad\quad \times P(z_s, \vec{\theta_{\rm L}} \mid {\rm SL}) \, d\vec{\theta^*}\, d\vec{\theta_{\rm L}}\, d\vec{\beta} \\
&= {\cal N}_{\rm L} \bigg\langle P({\rm obs} \mid \vec{\theta^*}, z_s, \vec{\theta_{\rm L}}, \vec{\beta}, {\rm SL}) \bigg\rangle_{\substack{
\vec{\theta^*}\sim P(\vec{\theta^*}) \\
\vec{\theta_{\rm L}}, z_s \sim P(z_s, \vec{\theta_{\rm L}} \mid {\rm SL}) \\
\vec{\beta} \sim P(\vec{\beta} \mid z_s, \vec{\theta_{\rm L}}, {\rm SL})
}} \,,
\end{split}
\end{equation}
$$

where $\langle \cdots \rangle$ denotes the average over samples drawn from the specified distributions.


## Decomposing the Joint Source–Lens Parameter Distribution

Let's use $\vec{\theta_{\rm L}^*}\in \{\sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2\}$, the lens related parameter without $z_L$. To have an idea of how to sample source redshift and lens paramters conditioned on strong lensing, $P(z_s, z_L, \vec{\theta_{\rm L}^*} \mid {\rm SL})$ needs to be further expanded using chain rule and Bayes' rule as

$$
\begin{equation}
\begin{split}
&P(z_s, z_L, \vec{\theta_{\rm L}^*} \mid {\rm SL})
\\
&\quad = P(\vec{\theta_{\rm L}^*} \mid z_L, z_s, {\rm SL}) \, P(z_L \mid z_s, {\rm SL}) \, P(z_s \mid {\rm SL}) \\
&\quad \propto P({\rm SL} \mid \vec{\theta_{\rm L}^*}, z_L, z_s) \\
&\quad\quad\times P(\vec{\theta_{\rm L}^*} \mid z_L, z_s) \\
&\quad\quad\times P(z_L \mid z_s, {\rm SL}) \\
&\quad\quad\times P({\rm SL} \mid z_s) \, P(z_s) \,, \\
&\quad \propto P({\rm SL} \mid \vec{\theta_{\rm L}^*}, z_L, z_s) \\
&\quad\quad\times P(\gamma_1, \gamma_2) \, P(\gamma) \, P(\phi_{\rm rot}) \, P(q \mid \sigma) \, P(\sigma \mid z_L) \\
&\quad\quad\times P(z_L \mid z_s, {\rm SL}) \\
&\quad\quad\times P({\rm SL} \mid z_s) \, P(z_s) \, ,
\end{split}
\end{equation}
$$

with the constant of proportionality $P({\rm SL} \mid z_L, z_s)\, P({\rm SL})$ absorbed into the normalization.

where $P({\rm SL} \mid \vec{\theta_{\rm L}^*}, z_L, z_s)$ is the cross-section giving the strong-lensing probability for a given source–lens configuration, while $P({\rm SL} \mid z_s)$ in the last term is the strong-lensing optical depth. Substituting this expression back into the event-rate equation gives

where each factor corresponds to a prior or conditional distribution used for sampling and numerical integration.

The sampling procedures should go like this: $z_s$ is first drawn from the it's intrinsic distribution weighted by the strong-lensing optical depth; then $z_L$ is drawn from the effective lensing kernel at that $z_s$; next the lens parameters are drawn from their intrinsic distributions weighted by the strong-lensing cross-section at that $(z_s,z_L)$. Each of these steps is described in the following sections.

## Source Redshift Distributions of Lensed Events

The redshift distribution of lensed sources, $P(z_s \mid {\rm SL})$, differs from the intrinsic source redshift distribution $P(z_s)$ because the probability of strong lensing depends on redshift. Using Bayes' theorem, the conditional probability can be dissociated as

$$
\begin{equation}
\begin{split}
P(z_s \mid {\rm SL})
= \frac{P({\rm SL} \mid z_s)\,P(z_s)}{P({\rm SL})} \,,
\end{split}
\end{equation}
$$

where $P({\rm SL} \mid z_s)$ is the strong-lensing optical depth. Substituting the [intrinsic redshift distribution](https://ler.hemantaph.com/analytical_formulation_unlensed.html#redshift-distribution-and-intrinsic-merger-rates) used for unlensed events gives

$$
\begin{equation}
\begin{split}
P(z_s \mid {\rm SL})
&= \frac{P({\rm SL} \mid z_s)}{P({\rm SL})} \left\{ \frac{1}{{\cal N}_{\rm U}} \frac{R_{\rm U}(z_s)}{1+z_s} \frac{dV_c}{dz_s} \right\}\\
&= \frac{1}{{\cal N}_{\rm L}}\,P({\rm SL} \mid z_s)\,\frac{R_{\rm U}(z_s)}{1+z_s}\,\frac{dV_c}{dz_s} \,,
\end{split}
\end{equation}
$$

where $R_{\rm U}(z_s)$ is the intrinsic merger rate density in the source frame, 
$\frac{dV_c}{dz_s}$ is the comoving volume element, $1/(1+z_s)$ is the cosmological time-dilation factor, and ${\cal N}_{\rm U}$ is the normalization constant for the unlensed population also representing the total unlensed intrinsic merger rate per year in the detector frame.

The normalization constant ${\cal N}_{\rm L}$ is given by

$$
\begin{equation}
\begin{split}
{\cal N}_{\rm L}
&= {\cal N}_{\rm U} \, P({\rm SL}) \,, \\
&= {\cal N}_{\rm U} \int_{z_s} P({\rm SL} \mid z_s)\, P(z_s) \, dz_s \,, \\
&= {\cal N}_{\rm U} \int_{z_s} P({\rm SL} \mid z_s) \left\{ \frac{1}{{\cal N}_{\rm U}} \frac{R_{\rm U}(z_s)}{1+z_s}\,\frac{dV_c}{dz_s} \right\} dz_s \,, \\
&= \int_{z_s} P({\rm SL} \mid z_s)\, \frac{R_{\rm U}(z_s)}{1+z_s}\, \frac{dV_c}{dz_s} \, dz_s \,.
\end{split}
\end{equation}
$$

The constant ${\cal N}_{\rm L}$ is the total rate of strongly lensed mergers per year, which I denoted earlier as $\frac{\Delta N_{\rm L}}{\Delta t}$. $P({\rm SL})$ is the overall strong-lensing probability, given by $P({\rm SL})={\cal N}_{\rm L}/{\cal N}_{\rm U}$ which is the fraction of lensed events in the total population.

In this work, the optical depth $P({\rm SL} \mid z_s)$ is computed for the adopted EPL+Shear lens model and it is detailed in the following section.

<!-- ## Strong-Lensing Optical Depth

The strong-lensing optical depth $P({\rm SL} \mid z_s)$ quantifies the probability that a source at redshift $z_s$ is strongly lensed by intervening galaxies. It is computed by integrating the effective lensing kernel $\Phi_{\rm EPL}(z_L,z_s)$ over all possible lens redshifts $z_L$ as -->

## Optical Depth

The optical depth for strong gravitational lensing (SL) of gravitational waves, $P({\rm SL}|z_s)$, gives the probability that a GW source at redshift $z_s$ is strongly lensed by intervening galaxies along the line of sight. Strong lensing occurs when a lens galaxy at redshift $z_L<z_s$ produces multiple images because the source position lies inside the double-image caustic region of the lens on the source plane.

Conceptually, the optical depth measures the fraction of the sky that is effectively covered by the strong-lensing angular cross-sections of all potential lenses between the observer and the source. With $\sigma_{\rm SL}(z_L,z_s)$, the cummulative angular cross-section of lenses at redshift $z_L$ for producing strong lensing of a source at $z_s$, and the all-sky given by $4\pi$, and $dN(z_L)$ represents the differential number of lenses in the redshift interval $dz_L$, it can be written as an integral over lens redshift as

$$
\begin{equation}
\begin{split}
P({\rm SL}|z_s)
&= \int_{0}^{z_s} \frac{\sigma_{\rm SL}(z_L,z_s)}{4\pi}\, dN(z_L) \\
&= \int_{0}^{z_s} P({\rm SL}|z_s,z_L)\, dN(z_L) \\
&= \int_{0}^{z_s} \left\{ \int_{\vec{\theta_{\rm L}^*}} P({\rm SL}|z_s,z_L,\vec{\theta_{\rm L}^*})\,\frac{d^2N(z_L,\vec{\theta_{\rm L}^*})}{dV_c\, d\vec{\theta_{\rm L}^*}}\,\frac{dV_c}{dz_L}\, d\vec{\theta_{\rm L}^*} \right\} dz_L \,,
\end{split}
\end{equation}
$$

<!-- where $\sigma_{\rm SL}(z_L,z_s)$ is the angular cross-section of a single lens at redshift $z_L$ for producing strong lensing of a source at $z_s$. The factor $1/(4\pi)$ normalizes this angular area by the total sky solid angle, so that $P({\rm SL}|z_s,z_L)=\sigma_{\rm SL}(z_L,z_s)/(4\pi)$ is the probability that a source at $z_s$ is strongly lensed by that lens. The term $dN(z_L)$ represents the differential number of lenses in the redshift interval $dz_L$. -->

The term $\frac{d^2N(z_L,\vec{\theta_{\rm L}^*})}{dV_c\, d\vec{\theta_{\rm L}^*}}$ is the comoving number density of lenses per unit comoving volume and per unit interval in $\vec{\theta_{\rm L}^*}$, and $\frac{dV_c}{dz_L}$ is the comoving volume element.

It is useful to define the effective lensing kernel $\Phi(z_L,z_s)$, which collects the contribution of all lenses at redshift $z_L$ after marginalizing over lens properties. The optical depth then becomes

$$
\begin{equation}
\begin{split}
P({\rm SL}|z_s)
= \int_{0}^{z_s} \Phi(z_L,z_s)\,dz_L \,,
\end{split}
\end{equation}
$$

with

$$
\begin{equation}
\begin{split}
\Phi(z_L,z_s)
= \int_{\vec{\theta_{\rm L}^*}} P({\rm SL}|z_s,z_L,\vec{\theta_{\rm L}^*})\,\frac{d^2N(z_L,\vec{\theta_{\rm L}^*})}{dV_c\, d\vec{\theta_{\rm L}^*}}\,\frac{dV_c}{dz_L}\, d\vec{\theta_{\rm L}^*} \,,
\end{split}
\end{equation}
$$

$\Phi(z_L,z_s)$ which also represents the effective lensing cross-section for lenses at redshift $z_L$, marginalized over the lens parameters $\vec{\theta_{\rm L}^*}$ and normalized by the total solid angle of the sky, can be use for sampling lens redshifts of lensed events as described in the next section.

## Lens Redshift Distribution of Lensed Events

In `ler`, $\Phi(z_L,z_s)$ is numerically computed, in the grid of $z_L$ and $z_s$, for the adopted EPL+Shear lens model by integrating over the lens parameter space using Monte Carlo sampling.

$$
\begin{equation}
\begin{split}
\Phi(z_L,z_s)
\approx \frac{1}{N_{\rm samp}} \sum_{i=1}^{N_{\rm samp}} P({\rm SL}|z_s,z_L,\vec{\theta_{\rm L,i}^*})\,\frac{d^2N(z_L,\vec{\theta_{\rm L,i}^*})}{dV_c\, d\vec{\theta_{\rm L}^*}}\,\frac{dV_c}{dz_L} \,,
\end{split}
\end{equation}
$$


The conditional lens-redshift distribution for strongly lensed sources, $P(z_L|z_s,{\rm SL})$, follows directly from $\Phi(z_L,z_s)$ and can be written as

$$
\begin{equation}
\begin{split}
P(z_L|z_s,{\rm SL})
= \frac{\Phi(z_L,z_s)}{\int_{0}^{z_s}\Phi(z_L,z_s)\,dz_L} \,,
\end{split}
\end{equation}
$$

which is normalized over $z_L\in[0,z_s]$ for each fixed source redshift $z_s$.


---

The conditional lens redshift distribution $P(z_L \mid z_s,{\rm SL})$ scales with the effective lensing kernel $\Phi_{\rm EPL}(z_L,z_s)$, which combines the lens comoving number density with the lensing cross-section. For EPL$+$shear lenses, it can be written as

$$
\begin{equation}
\begin{split}
P(z_L \mid z_s,{\rm SL})
= \frac{1}{{\cal N}_\Phi(z_s)}\,\Phi_{\rm EPL}(z_L,z_s) \,,
\end{split}
\end{equation}
$$

where ${\cal N}_\Phi(z_s)$ is the normalization factor for each fixed source redshift $z_s$.

## Lens Parameter Distribution and Sampling Approximation

The lens parameters are denoted by $\theta_{\rm L}\in\{\sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2\}$, where $\sigma$ is the velocity dispersion, $q$ is the projected axis ratio, $\phi_{\rm rot}$ is the lens orientation, and $(\gamma,\gamma_1,\gamma_2)$ describe external shear. The conditional distribution $P(\theta_{\rm L} \mid z_L,z_s,{\rm SL})$ can be constructed using Bayes' theorem as

$$
\begin{equation}
\begin{split}
P(\theta_{\rm L} \mid z_L, z_s, {\rm SL})
\propto P({\rm SL} \mid z_L, z_s, \theta_{\rm L})\,P(\theta_{\rm L} \mid z_L, z_s) \,,
\end{split}
\end{equation}
$$

where $P(\theta_{\rm L} \mid z_L, z_s)$ is the intrinsic prior over lens parameters. For EPL lenses, the strong-lensing probability scales with the double-image caustic cross-section $\sigma^{\rm DC}_{\rm EPL}$, so that

$$
\begin{equation}
\begin{split}
P({\rm SL} \mid z_L, z_s, \theta_{\rm L}) \propto \sigma^{\rm DC}_{\rm EPL} \,.
\end{split}
\end{equation}
$$

Direct evaluation of $\sigma^{\rm DC}_{\rm EPL}$ for each proposed lens is computationally expensive. To reduce cost, `ler` uses the analytic SIE cross-section $\sigma_{\rm SIE}$ as a proxy when performing rejection sampling for the core parameters $(\sigma,q)$. The remaining lens parameters are sampled from distributions conditioned on strong lensing. Under this approximation, the conditional lens-parameter distribution can be written as

$$
\begin{equation}
\begin{split}
P(\theta_{\rm L} \mid z_L, z_s, {\rm SL})
= P(\sigma,q \mid z_L, z_s, {\rm SL})\,
P(\phi_{\rm rot} \mid {\rm SL})\,
P(\gamma \mid {\rm SL})\,
P(\gamma_1,\gamma_2 \mid {\rm SL}) \,,
\end{split}
\end{equation}
$$

with the proxy weight

$$
\begin{equation}
\begin{split}
P({\rm SL} \mid \sigma,q,z_L,z_s) \propto \sigma_{\rm SIE} \,.
\end{split}
\end{equation}
$$

## Source Position Distribution and Image Properties

The source position in the source plane is denoted by $\beta\in\{\beta_x,\beta_y\}$, representing the angular offset between the source and the lens center. Under the strong-lensing condition, `ler` samples $\beta$ uniformly within the region enclosed by the tangential caustic of the lens model. For each sampled configuration, the lens equation is solved using [`lenstronomy`](https://lenstronomy.readthedocs.io), yielding the image positions $\theta_i$, arrival-time delays $t_i$, and magnifications $\mu_i$ for each image.

Image types are identified using the trace and determinant of the Jacobian matrix of the lens mapping, which determines the Morse index $n_i$. The images are classified as Type-I (minima, $n_i=0$), Type-II (saddle points, $n_i=1/2$), and Type-III (maxima, $n_i=1$). In this study, central demagnified Type-III images are neglected. Images are ordered by arrival time such that $t_1<t_2<t_3<t_4$, and configurations with up to four images are considered.

Strong lensing modifies the observed GW signal for each image. The coalescence time becomes $t_c+t_i$, the effective luminosity distance becomes $d_L/\sqrt{ \mid \mu_i \mid }$, and the coalescence phase is shifted to $\phi_c-n_i\pi$. The subscript $i$ denotes the $i^{\rm th}$ image.
