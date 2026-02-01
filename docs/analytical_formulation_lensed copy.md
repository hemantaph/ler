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

The conditional probability $P({\rm obs \mid SL})$ is evaluated by marginalizing over the relevant GW (unlensed) source parameters $\vec{\theta_{\rm U}}$, the lens parameters $\vec{\theta_{\rm L}}$, and the source position relative to the lens $\beta$. For $\vec{\theta_{\rm U}}$ see example [Analytical Formulation for Gravitational Wave Event Rates](https://ler.hemantaph.com/analytical_formulation_unlensed.html). For Elliptical Power Law with external shear lens models (EPL+Shear), these sets of parameters $\vec{\theta_{\rm L}}$ are defined as

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
&= \frac{\Delta N_{\rm L}}{\Delta t}
\int_{\vec{\theta_{\rm U}}}\int_{\vec{\theta_{\rm L}}}\int_{\vec{\beta}} P({\rm obs} \mid \vec{\theta_{\rm U}}, \vec{\theta_{\rm L}}, \vec{\beta}, {\rm SL}) \\
&\quad\quad\quad\quad\quad\quad\quad\quad \times P(\vec{\theta_{\rm U}}, \vec{\theta_{\rm L}}, \vec{\beta} \mid {\rm SL}) d\vec{\theta_{\rm U}}\, d\vec{\theta_{\rm L}}\, d\vec{\beta} \,,
\end{split}
\end{equation}
$$

where $P({\rm obs} \mid \vec{\theta_{\rm U}}, \vec{\theta_{\rm L}}, \vec{\beta}, {\rm SL})$ is the detection probability for a specific source–lens configuration, and $P(\vec{\theta_{\rm U}}, \vec{\theta_{\rm L}}, \vec{\beta} \mid {\rm SL})$ is the joint distribution of parameters conditioned on strong lensing. Assuming that the intrinsic GW parameters $\vec{\theta}^*_{\rm U}$ (* indicating without $z_s$) are independent of the detailed lensing configuration (and the strong lensing condition) once $z_s$ is specified and the event is known to be lensed, this joint distribution can be decomposed using the chain rule as

$$
\begin{equation}
\begin{split}
P(\vec{\theta_{\rm U}}, \vec{\theta_{\rm L}}, \vec{\beta} \mid {\rm SL}) &= P(\vec{\theta}^*_{\rm U}) \, P(\vec{\beta} \mid z_s, \vec{\theta_{\rm L}}, {\rm SL}) \, P(z_s, \vec{\theta_{\rm L}} \mid {\rm SL}) \,.
\end{split}
\end{equation}
$$

Numerically, the integral is evaluated using Monte Carlo sampling over the full parameter space. Substituting this expression back into the event-rate equation gives

$$
\begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}
&= \frac{\Delta N_{\rm L}}{\Delta t}
\int_{\vec{\theta}^*_{\rm U}}\int_{\vec{\theta_{\rm L}}}\int_{\vec{\beta}}
 P({\rm obs} \mid \vec{\theta}^*_{\rm U}, z_s, \vec{\theta_{\rm L}}, \vec{\beta}, {\rm SL})\,
P(\vec{\theta}^*_{\rm U}) \, P(\vec{\beta} \mid z_s, \vec{\theta_{\rm L}}, {\rm SL}) \\
&\quad\quad\quad\quad\quad\quad\quad\quad \times P(z_s, \vec{\theta_{\rm L}} \mid {\rm SL}) \, d\vec{\theta}^*_{\rm U}\, d\vec{\theta_{\rm L}}\, d\vec{\beta} \\
&= {\cal N}_{\rm L} \bigg\langle P({\rm obs} \mid \vec{\theta}^*_{\rm U}, z_s, \vec{\theta_{\rm L}}, \vec{\beta}, {\rm SL}) \bigg\rangle_{\substack{
\vec{\theta}^*_{\rm U}\sim P(\vec{\theta}^*_{\rm U}) \\
\vec{\theta_{\rm L}}, z_s \sim P(z_s, \vec{\theta_{\rm L}} \mid {\rm SL}) \\
\vec{\beta} \sim P(\vec{\beta} \mid z_s, \vec{\theta_{\rm L}}, {\rm SL})
}} \,,
\end{split}
\end{equation}
$$

where $\langle \cdots \rangle$ denotes the average over samples drawn from the specified distributions.


## Decomposing the Joint Source–Lens Parameter Distribution

For EPL+Shear lenses, let's use $\vec{\theta_{\rm L}^*}\in \{\sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2\}$, which is the the lens related parameter without $z_L$. To have an idea of how to sample source redshift and lens paramters conditioned on strong lensing, $P(z_s, z_L, \vec{\theta_{\rm L}^*} \mid {\rm SL})$ needs to be further expanded using chain rule and Bayes' rule as

$$
\begin{equation}
\begin{split}
&P(z_s, z_L, \vec{\theta_{\rm L}^*} \mid {\rm SL})
\\
&\quad = P(\vec{\theta_{\rm L}^*} \mid z_L, z_s, {\rm SL}) \, P(z_L \mid z_s, {\rm SL}) \, P(z_s \mid {\rm SL}) \\
&\quad \propto P({\rm SL} \mid \vec{\theta_{\rm L}^*}, z_L, z_s) \\
&\quad\quad\times P(\vec{\theta_{\rm L}^*} \mid z_L, z_s) \\
&\quad\quad\times P(z_L \mid z_s, {\rm SL}) \\
&\quad\quad\times P({\rm SL} \mid z_s) \, P(z_s) \,,
\end{split}
\end{equation}
$$

and for EPL+Shear lenses, it becomes

$$
\begin{equation}
\begin{split}
&P(z_s, z_L, \sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2 \mid {\rm SL})
\\
&\quad \propto P({\rm SL} \mid \sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2, z_L, z_s) \\
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

<!-- ## Stro
ng-Lensing Optical Depth

The strong-lensing optical depth $P({\rm SL} \mid z_s)$ quantifies the probability that a source at redshift $z_s$ is strongly lensed by intervening galaxies. It is computed by integrating the effective lensing kernel $\Phi_{\rm EPL}(z_L,z_s)$ over all possible lens redshifts $z_L$ as -->

## Optical Depth

The optical depth for strong gravitational lensing (SL) of gravitational waves, $P({\rm SL}|z_s)$, gives the probability that a GW source at redshift $z_s$ is strongly lensed by intervening galaxies along the line of sight. Strong lensing occurs when a lens galaxy at redshift $z_L<z_s$ produces multiple images because the source position lies inside the [double-image caustic](#) region (which also include the diamond caustic) of the lens on the source plane, resulting in multiple detectable images.

Conceptually, the optical depth measures the fraction of the sky that is effectively covered by the strong-lensing angular cross-sections of all potential lenses between the observer and the source. With $\sigma_{\rm SL}(z_L,z_s)$, the cummulative angular cross-section of lenses at redshift $z_L$ for producing strong lensing of a source at $z_s$, and the all-sky given by $4\pi$, and $dN(z_L)$ representing the differential number of lenses in the redshift interval $dz_L$, optical-depth can be written as an integral over lens redshift as

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

It is useful to define the effective lensing kernel $\Phi_{\rm SL}(z_L,z_s)$, which collects the contribution of all lenses at redshift $z_L$ after marginalizing over lens properties. The optical depth then becomes

$$
\begin{equation}
\begin{split}
P({\rm SL}|z_s)
= \int_{0}^{z_s} \Phi_{\rm SL}(z_L,z_s)\,dz_L \,,
\end{split}
\end{equation}
$$

with

$$
\begin{equation}
\begin{split}
\Phi_{\rm SL}(z_L,z_s)
= \int_{\vec{\theta_{\rm L}^*}} P({\rm SL}|z_s,z_L,\vec{\theta_{\rm L}^*})\,\frac{d^2N(z_L,\vec{\theta_{\rm L}^*})}{dV_c\, d\vec{\theta_{\rm L}^*}}\,\frac{dV_c}{dz_L}\, d\vec{\theta_{\rm L}^*} \,,
\end{split}
\end{equation}
$$

and with for EPL+Shear lenses, it reads,

$$
\begin{equation}
\begin{split}
\Phi_{\rm SL}(z_L, z_s)
&= \int_{\sigma,q,\phi_{\rm rot},\gamma,\gamma_1,\gamma_2} P({\rm SL}|z_s,z_L,\sigma,q,\phi_{\rm rot},\gamma,\gamma_1,\gamma_2) \\
& \quad \times P(\gamma_1, \gamma_2) \, P(\gamma) \, P(\phi_{\rm rot}) \, P(q \mid \sigma) \\
& \quad \times \frac{d^2N(z_L, \sigma)}{dV_c \, d\sigma} \, \frac{dV_c}{dz_L} \, d\sigma \, dq \, d\phi_{\rm rot} \, d\gamma \, d\gamma_1 \, d\gamma_2 \\
&= \Delta \sigma \left\langle 
\frac{\sigma^{\rm DC}_{\rm EPL}}{4\pi} \, \phi(\sigma, z_L) \, \frac{dV_c}{dz_L} 
\right\rangle_{\substack{
q \in P(q|\sigma), \\
\phi_{\rm rot} \in P(\phi_{\rm rot}), \\
\gamma \in P(\gamma), \\
\gamma_1, \gamma_2 \in P(\gamma_1, \gamma_2), \\
\sigma \in P_o(\sigma)
}} \,,
\end{split}
\end{equation}
$$

where $\phi(\sigma, z_L)$ is the velocity-dispersion function giving the comoving number density of lenses per unit velocity-dispersion interval, and $\sigma^{\rm DC}_{\rm EPL}$ is the double-image caustic cross-section for EPL+Shear lenses. The last line shows how $\Phi_{\rm SL}(z_L,z_s)$ is numerically computed using Monte Carlo sampling over the lens parameter space, with $\Delta \sigma$ arrising from using uniform sampling in $\sigma$ with $P_o(\sigma)=\frac{1}{\sigma_{\max}-\sigma_{\min}}=\frac{1}{\Delta \sigma}$.

$\Phi_{\rm SL}(z_L,z_s)$ which also represents the effective lensing cross-section for lenses at redshift $z_L$, marginalized over the lens parameters $\vec{\theta_{\rm L}^*}$ and normalized by the total solid angle of the sky, and it can be use for sampling lens redshifts of lensed events as described in the next section.

## Lens Redshift Distribution of Lensed Events

The conditional lens-redshift distribution for strongly lensed sources, $P(z_L|z_s,{\rm SL})$, follows directly from $\Phi_{\rm SL}(z_L,z_s)$ and can be written as

$$
\begin{equation}
\begin{split}
P(z_L|z_s,{\rm SL})
&= \frac{\Phi_{\rm SL}(z_L,z_s)}{\int_{0}^{z_s}\Phi_{\rm SL}(z_L,z_s)\,dz_L} \\
&= \frac{\Phi_{\rm SL}(z_L,z_s)}{{\cal N}_\Phi(z_s)} \,,
\end{split}
\end{equation}
$$

with normalization factor ${\cal N}_\Phi(z_s)$ at each $z_s$.

In `ler`, $\Phi_{\rm SL}(z_L,z_s)$ and $P(z_L|z_s,{\rm SL})$ are numerically precomputed, in the grid of $z_L$ and $z_s$, and for the adopted EPL+Shear lens model. 2D interpolation is then used to evaluate these functions at arbitrary $(z_L,z_s)$. For sampling at fixed $z_s$, the inverse-transform sampling method is employed using the cumulative distribution function constructed from $P(z_L|z_s,{\rm SL})$.

If we want the intrinsic lens redshift distribution without the strong-lensing condition, it can be obtained by marginalizing over the lens parameters as

$$
\begin{equation}
\begin{split}
P(z_L) = \frac{1}{{\cal N}_{z_L}} \int_{\vec{\theta_{\rm L}^*}} \frac{d^2N(z_L,\vec{\theta_{\rm L}^*})}{dV_c\, d\vec{\theta_{\rm L}^*}}\,\frac{dV_c}{dz_L}\, d\vec{\theta_{\rm L}^*}   \,,
\end{split}
\end{equation}
$$

where ${\cal N}_{z_L}$ is the normalization constant.

For EPL+Shear lenses, it becomes

$$
\begin{equation}
\begin{split}
P(z_L)
&= \frac{1}{{\cal N}_{z_L}} \Delta \sigma \left\langle 
\phi(\sigma, z_L) \, \frac{dV_c}{dz_L} 
\right\rangle_{\substack{
q \in P(q|\sigma), \\
\phi_{\rm rot} \in P(\phi_{\rm rot}), \\
\gamma \in P(\gamma), \\
\gamma_1, \gamma_2 \in P(\gamma_1, \gamma_2), \\
\sigma \in P_o(\sigma)
}}  \,.
\end{split}
\end{equation}
$$

## Double-Image Caustic as Lensing Cross-Section

The strong-lensing condition applied in our $\vec{\theta_{\rm L}^*}$ and $z_s$ sampling, already implied that lens and source are aligned such that multiple images (2 or more) are formed. If not, the configuration is discarded during sampling and a new one is drawn.

For EPL+Shear lenses, strong lensing producing multiple images occurs when the source position lies within the double-image caustic region on the source plane. The angular cross-section for this region, denoted as $\sigma^{\rm DC}_{\rm EPL}$, serves as the effective lensing cross-section for strong lensing. Thus, 

$$
\begin{equation}
\begin{split}
P({\rm SL} \mid z_L, z_s, \vec{\theta_{\rm L}^*}) = \frac{\sigma^{\rm DC}_{\rm EPL}}{4\pi} \,.
\end{split}
\end{equation}
$$

Typically, $\sigma^{\rm DC}_{\rm EPL}$ is computed numerically by tracing the caustic curves on the source plane using [`lenstronomy`](https://lenstronomy.readthedocs.io). This is computationally expensive and time consuming. However, in `ler`, we implement multi-dimensiional interpolation of precomputed $\sigma^{\rm DC}_{\rm EPL}$ values over the lens parameter space ($q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2$) with unit Einstein radius $\theta_{\rm E}=1$, and scaling with the $\theta_{\rm E}$ for arbitrary lens parameters ($\sigma, z_L, z_s$), with $\theta_{\rm E}$ given by

$$
\begin{equation}
\begin{split}
\theta_{\rm E} = 4\pi \left( \frac{\sigma}{c} \right)^2 \frac{D_{LS}(z_L,z_s)}{D_S(z_s)} \,,
\end{split}
\end{equation}
$$

where $D_{LS}$ and $D_S$ are the angular diameter distances between lens and source, and observer and source, respectively. 

<!-- area_sis = np.pi * theta_E**2
cs = cs_unit * (csunit_to_cs_intercept + csunit_to_cs_slope * area_sis) -->

The scaling relation is given by
$$
\begin{equation}
\begin{split}
&\sigma^{\rm DC}_{\rm EPL}(\theta_{\rm E}(\sigma, z_L, z_s), q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2) \\
&\quad\quad= \sigma^{\rm DC}_{\rm EPL}(\theta_{\rm E}=1, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2) \\
&\quad\quad\quad\times \left({\rm intercept} + {\rm slope} \times \pi \, \theta_{\rm E}^2 \right) \,,
\end{split}
\end{equation}
$$

where `intercept` and `slope` are precomputed constants. This approach significantly speeds up the computation of $\sigma^{\rm DC}_{\rm EPL}$ during sampling and also enables speedy calculation of the effective lensing kernel $\Phi(z_L,z_s)$ and thus the optical depth $P({\rm SL}|z_s)$. 

## Lens Parameter Distributions of Lensed Events

The conditional lens-parameter distribution for strongly lensed sources, $P(\vec{\theta_{\rm L}^*} \mid z_L, z_s, {\rm SL})$ at fixed $z_L$ and $z_s$, for EPL+Shear lenses, can be expressed as

$$
\begin{equation}
\begin{split}
&P(\sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2 \mid z_L, z_s, {\rm SL})\\
&\quad \propto P({\rm SL} \mid \sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2, z_L, z_s) \\
&\quad\quad \times P(\gamma_1, \gamma_2) \, P(\gamma) \, P(\phi_{\rm rot}) \, P(q \mid \sigma) \, P(\sigma \mid z_L) \\
&\quad \propto \sigma^{\rm DC}_{\rm EPL} \, P(\gamma_1, \gamma_2) \, P(\gamma) \, P(\phi_{\rm rot}) \, P(q \mid \sigma) \, P(\sigma \mid z_L) \,.
\end{split}
\end{equation}
$$

Rejection sampling is in-efficient as for $P(\sigma \mid z_L)$, as most of the drawn configurations will not satisfy the strong-lensing condition. Moreover, maximum value of $\sigma^{\rm DC}_{\rm EPL}$ is not known analytically for EPL+Shear lenses, making it difficult to use the standard rejection sampling method. To overcome these challenges, `ler` employs an importance sampling approach. Rewriting the sampling equation as

$$
\begin{equation}
\begin{split}
&P(\sigma, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2 \mid z_L, z_s, {\rm SL})\\
&\quad \propto\frac{ \sigma^{\rm DC}_{\rm EPL}\, P(\sigma \mid z_L)}{P_o(\sigma)} \, P(\gamma_1, \gamma_2) \, P(\gamma) \, P(\phi_{\rm rot}) \, P(q \mid \sigma) \, P_o(\sigma) \,,
\end{split}
\end{equation}
$$

where $P_o(\sigma)$ is a convenient proposal distribution for $\sigma$ that allows efficient sampling such as uniform distribution $\left( P_o(\sigma)=\frac{1}{\sigma_{\max}-\sigma_{\min}}\right)$. 

We use $\frac{ \sigma^{\rm DC}_{\rm EPL}\, P(\sigma \mid z_L)}{P_o(\sigma)}$ as the importance weight for each sampled configuration. The lens parameters are then drawn from their intrinsic distributions, with $\sigma$ drawn from the proposal distribution $P_o(\sigma)$.

## Source Position Distribution and Image Properties

We set $\theta_{\rm E}$ as the unit of angular scale on the source plane for sampling source positions and lens equation solving. The image properties are later rescaled to physical units using the actual $\theta_{\rm E}$ for each lens configuration.

Thus, the source position $\vec{\beta}$ is sampled uniformly within the double-image caustic region defined by the lens parameters ($\theta_{\rm E}=1, q, \phi_{\rm rot}, \gamma, \gamma_1, \gamma_2$), resulting in multiple images (2, 3 or 4 etc.) for each sampled configuration.

For each sampled configuration, the lens equation is solved using [`lenstronomy`](https://lenstronomy.readthedocs.io), yielding the image positions $\vec{\theta_i}\in\{\theta_{x,i}, \theta_{y,i}\}$, arrival-time delays $t_i$, and magnifications $\mu_i$ for each image, with the subscript $i$ denotes the $i^{\rm th}$ image.. $\vec{\beta}$, $\theta_i$, and $t_i$ are then re-scaled using the actual $\theta_{\rm E}$.

Image types are identified using the trace and determinant of the Jacobian matrix of the lens mapping, which determines the Morse index $n_i$. The images are classified as Type-I (minima, $n_i=0$), Type-II (saddle points, $n_i=1/2$), and Type-III (maxima, $n_i=1$). In this study, central demagnified Type-III images are neglected. Images are ordered by arrival time such that $t_1<t_2<t_3<t_4$, and configurations with up to four images are considered.

Strong lensing modifies the observed GW signal for each image. 
<!-- The coalescence time becomes $t_c+t_i$, the effective luminosity distance becomes $d_L/\sqrt{ \mid \mu_i \mid }$, and the coalescence phase is shifted to $\phi_c-n_i\pi$.  -->
So the GW parameters for each image are updated to

$$
\begin{equation}
\begin{split}
\vec{\theta_{{\rm GW},i}} &\in \Bigg\{
\frac{d_L(z_s)}{\sqrt{\lvert \mu_i\rvert}},\;
\phi_c - n_i\pi,\;
t_c + t_i,\; \\
&\quad\quad {\rm RA} + \frac{\theta_{x,i}-\beta_x}{\cos({\rm Dec})},\;
{\rm Dec} + (\theta_{y,i}-\beta_y),\;\ldots
\Bigg\}
\end{split} \notag
\end{equation}
$$

the detector-frame masses and spins remain unchanged from the unlensed values. Even though the sky-location changes (with image positions) due to lensing, the effect is negligible for GW detectors given their poor angular resolution, and also this will not affect the SNR calculation significantly. SNR for each image and the detectability are then evaluated using these updated GW parameters as described in the next section.

## Detection Probability of Lensed Events

The detection probability for a specific lensed source–lens configuration, $P({\rm obs} \mid \vec{\theta}^*_{\rm U}, z_s, \vec{\theta_{\rm L}}, \vec{\beta}, {\rm SL})$, is determined by evaluating the detectability of each image produced by the lensing. For each image $i$, the signal-to-noise ratio (SNR) is computed using the updated GW parameters $\vec{\theta}_{{\rm GW},i}$ that account for lensing effects. An event is considered detectable if at least two of its images meets the detection criteria. Thus, the detection probability can be expressed as
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

where $\Theta$ is the Heaviside step function, and $\rho_{th}$ is the SNR detection threshold. The summation runs over all images produced by the lensing configuration. If two or more images have SNRs exceeding the threshold, the lensed event is deemed detectable, and the detection probability is set to 1; otherwise, it is set to 0.

## Complete Expression for Lensed Event Rate

Combining all the components discussed above, the final expression for the annual rate of detectable strongly lensed GW events is given by

$$
\begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm L}}{\Delta t}
&= {\cal N}_{\rm L} \bigg\langle P({\rm obs} \mid \vec{\theta}^*_{\rm U}, z_s, \vec{\theta_{\rm L}}, \vec{\beta}, {\rm SL}) \bigg\rangle_{\substack{
\vec{\theta}^*_{\rm U}\sim P(\vec{\theta}^*_{\rm U}) \\
z_s \sim P(z_s \mid {\rm SL}) \\
z_L \sim P(z_L \mid z_s, {\rm SL}) \\
\vec{\theta_{\rm L}^*} \sim P(\vec{\theta_{\rm L}^*} \mid z_L, z_s, {\rm SL}) \\
\vec{\beta} \sim P(\vec{\beta} \mid z_s, \vec{\theta_{\rm L}}, {\rm SL})
}} \,,
\end{split}
\end{equation}
$$


where the lens related parameter prior and description are listed in the following table.

## Lens Parameter Priors

Redshifts with strong-lensing condition:

| Parameter | Unit | Prior <br>Distribution | Range <br>[Min, Max] | Description |
| :--- | :--- | :--- | :--- | :--- |
| $z_s$ | - | $\propto \frac{P({\rm SL} \mid z_s) R_{\rm U}(z_s)}{1+z_s}\,\frac{dV_c}{dz_s}$<br> Intrinsic merger rate<br> weighted by SL<br> optical depth | [0.0, 10.0] | Source redshift |
| $z_L$ | - | $\frac{\Phi(z_L,z_s)}{{\cal N}_\Phi(z_s)}$<br> Integrated number<br> density weighted<br> by SL cross-section | [0.0, $z_s$] | Lens redshift |

Lens parameters (intrinsic distributions):

| Parameter | Unit | Prior <br>Distribution | Range <br>[Min, Max] | Description |
| :--- | :--- | :--- | :--- | :--- |
| $\sigma$ | km/s | $\frac{d^2N(z_L, \sigma)}{dV_c \, d\sigma}$<br> number density | [100, 400] | Velocity dispersion<br> Oguri et al. 2018 |
| $q$ | - | Rayleigh<br> distribution | [0.2, 1.0] | Axis ratio<br> Collett et al. 2015 |
| $\phi_{\rm rot}$ | rad | Uniform | [0, $\pi$] | Lens orientation |
| $\gamma$ | - | Normal | - | Density profile slope<br> (Mean: 2.0, Std: 0.1)<br> Sonnenfeld et al. 2024 |
| $\gamma_1, \gamma_2$ | - | Normal | - | External shear components<br> (Mean: 0.0, Std: 0.05)<br> Collet et al. 2015 |

These intrinsic distributions are used for sampling lens parameters of lensed events, and then weighted by the strong-lensing cross-section as described in the previous sections.