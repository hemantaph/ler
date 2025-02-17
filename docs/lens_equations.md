# Lensed event sampling (Analytical formulation)

## Initial setup

Note: I will interchangeably use terms like observable events and detectable events.

Define all the parameters involved.

- Source parameters: $\theta \in \{$ $m_1$ (mass of the heavier one), $m_2$ (mass of the lighter one), $\iota$ (inclination-angle), $\phi$ (phase-of-coalescence), $\psi$ (polarization-angle), $ra$ (right-ascension), $dec$ (declination) $\}$ and $z_s$ : red-shift of the source.
- Lens parameters: $\theta_L \in \{$ $\sigma$ (velocity-dispersion), $q$ (axis-ratio), $\psi$ (axis-rotation), $\gamma$ (spectral-index), $[\gamma_1,\gamma_2]$ (external-shear), $[e_1,e_2]$ (ellipticity), $\beta$ (source position) $\}$
- $z_L$ : red-shift of the galaxy lens
- image properties: $\{$ $\mu_i$ (magnification), $dt_i$ (time-delay), $n_i$ (morse-phase) $\}$. There is subscript $i$ because there can be multiple images.
- $SL$: Strong lensing condition. It means the original gravitational wave signal is lensed and split into multiple images $(\geq 2)$.
- Distance representation of lens object and source object will be in terms of redshifts ($z_l$ and $z_s$).

Montecarlo integration: Given a function $f(\theta)$ and normalized probability distribution function PDF $P(\theta)$, then,

$$
\begin{equation}
\begin{split}
\int f(\theta) P(\theta) d\theta &\approx \frac{1}{N} \sum_{i=1}^{N} f(\theta_i) \nonumber \\
&\approx \left< f(\theta) \right>_{\theta \in P(\theta)}
\end{split}
\end{equation}
$$

where $\theta_i$ are samples drawn from $P(\theta)$.

## General formulation of rate equation.

What we need to find is the number of lensed gravitational events that are detectable by the detector in a given time interval (1 year). Forget about lensing and detection for a moment. If we want to know the total number of compact binary coalescence (CBC) events then we can write it as,

$$
\begin{equation}
\begin{split}
\text{Rate of CBC events at redshift }z_s &= \frac{\text{Number CBC events at redshift }z_s}{\text{Time}} \nonumber \\
&= \text{Merger rate density at }z_s\times \text{Volume at }z_s \\
\end{split}
\end{equation}
$$

$$
\begin{equation}
\begin{split}
&\text{Rate of CBC events in the entire observable universe} \nonumber \\
&\;\;\;= \sum_{z_{min}}^{z_{max}} \text{Rate of CBC events at redshift }z_s \nonumber \\
\end{split}
\end{equation}
$$

Similarly, we can write the rate of lensed detectable events as,

$$
\begin{equation}
\begin{split}
&\text{Rate of detectable lensed events} \nonumber \\
&\;\;\;= \sum_{z_{min}}^{z_{max}} \text{observable rate density at }z_s\text{ which are also lensed}\times \text{Volume at }z_s
\end{split}
\end{equation}
$$

Merger rate density at $z_s$ is known. 

Let's try to formulate it with proper mathematics.

Given $d N^L_{obs}(z_s)$ is the number of lensed GW detectable events from sources at red-shift $z_s$ in a spherical shell of thickness $d z_s$, then, the rate observing lensed GW events is given by,

$$
\begin{equation}
\begin{split}
\mathcal{R}_L &= \int_{z_{min}}^{z_{max}} \frac{d N^L_{obs}(z_s)}{d t} \\
&= \int_{z_{min}}^{z_{max}} \frac{d^2 N^L_{obs}(z_s)}{d t \;d\mathcal{V}_c}  \frac{d\mathcal{V}_c}{d\Omega dz_s} d\Omega dz_s \\
&= \int_{z_{min}}^{z_{max}} \frac{d^2 N^L_{obs}(z_s)}{\left(d t/(1+z_s)\right) \;d\mathcal{V}_c}\frac{1}{(1+z_s)}\frac{d\mathcal{V}_c}{d\Omega dz_s} d\Omega dz_s \\
&= \int_{z_{min}}^{z_{max}} \frac{d^2 N^L_{obs}(z_s)}{d\tau \;d\mathcal{V}_c}\frac{1}{(1+z_s)} \frac{d\mathcal{V}_c}{d\Omega dz_s} d\Omega dz_s \\
&= \int_{z_{min}}^{z_{max}} \frac{d^2 N^L_{obs}(z_s)}{d\tau \;d\mathcal{V}_c}\frac{1}{(1+z_s)} \frac{d\mathcal{V}_c}{d\Omega dz_s} \Bigg\{ \int_{all\; sky} d\Omega \Bigg\}dz_s \\
&= \int_{z_{min}}^{z_{max}} \frac{d^2 N^L_{obs}(z_s)}{d\tau \;d\mathcal{V}_c}\frac{1}{(1+z_s)} \Bigg(\frac{d\mathcal{V}_c}{d\Omega dz_s} 4\pi \Bigg)dz_s \\
&= \int_{z_{min}}^{z_{max}} \frac{d^2 N^L_{obs}(z_s)}{d\tau \;d\mathcal{V}_c}\frac{1}{(1+z_s)} \frac{dV_c}{dz_s} dz_s \\
\end{split} \tag{1}
\end{equation}
$$

Note: Reader should not confuse between lensed events and lensed images. Each lensed event can have multiple lensed images. These images can be detected as seperate GW signals. In case of un-lensed events, there is only one GW signal for each event.

Given $d N^L_{obs}(z_s)$ is the number of lensed GW detectable events from sources at red-shift $z_s$ in a spherical shell of thickness $d z_s$, then, the rate observing lensed GW events is given by,

$$
\begin{equation}
\begin{split}
\mathcal{R}_L &= \int_{z_{min}}^{z_{max}} \frac{d N^L_{obs}(z_s)}{d t} \\
&= \int_{z_{min}}^{z_{max}} \frac{d^2 N^L_{obs}(z_s)}{d t \;d\mathcal{V}_c}  \frac{d\mathcal{V}_c}{d\Omega dz_s} d\Omega dz_s \\
&= \int_{z_{min}}^{z_{max}} \frac{d^2 N^L_{obs}(z_s)}{\left(d t/(1+z_s)\right) \;d\mathcal{V}_c}\frac{1}{(1+z_s)}\frac{d\mathcal{V}_c}{d\Omega dz_s} d\Omega dz_s \\
&= \int_{z_{min}}^{z_{max}} \frac{d^2 N^L_{obs}(z_s)}{d\tau \;d\mathcal{V}_c}\frac{1}{(1+z_s)} \frac{d\mathcal{V}_c}{d\Omega dz_s} d\Omega dz_s \\
&= \int_{z_{min}}^{z_{max}} \frac{d^2 N^L_{obs}(z_s)}{d\tau \;d\mathcal{V}_c}\frac{1}{(1+z_s)} \frac{d\mathcal{V}_c}{d\Omega dz_s} \Bigg\{ \int_{all\; sky} d\Omega \Bigg\}dz_s \\
&= \int_{z_{min}}^{z_{max}} \frac{d^2 N^L_{obs}(z_s)}{d\tau \;d\mathcal{V}_c}\frac{1}{(1+z_s)} \Bigg(\frac{d\mathcal{V}_c}{d\Omega dz_s} 4\pi \Bigg)dz_s \\
&= \int_{z_{min}}^{z_{max}} \frac{d^2 N^L_{obs}(z_s)}{d\tau \;d\mathcal{V}_c}\frac{1}{(1+z_s)} \frac{dV_c}{dz_s} dz_s \\
\end{split} \tag{1}
\end{equation}
$$

Note: Reader should not confuse between lensed events and lensed images. Each lensed event can have multiple lensed images. These images can be detected as seperate GW signals. In case of un-lensed events, there is only one GW signal for each event.

Note that if $d\tau=dt/(1+z_s)$ is considered the proper time and it can be converted from the time at detector frame $dt$ using the time-dilation factor $(1+z_s)$. Consequently, $\frac{d^2 N^L_{obs}(z_s)}{d t \;d\mathcal{V}_c}$ and $\frac{d^2 N^L_{obs}(z_s)}{d \tau \;d\mathcal{V}_c}$ are the observed merger rate density at detector-frame and source-frame respectively. We want to use the $R^L_{obs}(z_s)=\frac{d^2 N^L_{obs}(z_s)}{d \tau \;d\mathcal{V}_c}$ for our analysis as most observational papers and the output of theoretical predictions are in the source-frame. $\frac{dV_c}{dz_s}dz_s$ is considered a spherical-shell volume element in co-moving coordinates at red-shift $z_s$. So, the rate expression simplifies to integrating (density) $\times$ (time dilation effect) over the shell volume element. For more information on the volume element refer to this [page]() of the documentation.

Note: $\frac{dV_c}{dz_s}$ is the differential co-moving volume at red-shift $z_s$ and you can get the value by using `astropy` cosmology package for a given cosmology.

Now we have,

$$
\begin{equation}
\begin{split}
\mathcal{R}_L &= \int_{z_{min}}^{z_{max}} \frac{d^2 N^L_{obs}(z_s)}{d \tau\;d\mathcal{V}_c} \frac{1}{(1+z_s)} \frac{dV_c}{dz_s} dz_s \\
&= \int_{z_{min}}^{z_{max}} R^L_{obs}(z_s) \frac{1}{(1+z_s)} \frac{dV_c}{dz_s} dz_s
\end{split} \tag{2}
\end{equation}
$$

We want to re-write it in terms intrinsic merger rate distribution $R(z_s)$. $R(z_s)$ is the merger rate density distribution regardless of whether the source is detectable or not.

$$
\begin{equation}
\begin{split}
\mathcal{R}_L &= \int_{z_{min}}^{z_{max}} R^L_{obs}(z_s) \frac{1}{(1+z_s)} \frac{dV_c}{dz_s} dz_s \\
&= \int_{z_{min}}^{z_{max}} R(z_s) P(obs, SL|z_s) \frac{1}{(1+z_s)} \frac{dV_c}{dz_s} dz_s \\
&= \int_{z_{min}}^{z_{max}} R(z_s) P(obs|z_s, SL) P(SL|z_s) \frac{1}{(1+z_s)} \frac{dV_c}{dz_s} dz_s \\
\end{split} \tag{3}
\end{equation}
$$

$$
\begin{equation}
\begin{split}
\mathcal{R}_L &= \mathcal{N}^L \int_{z_{min}}^{z_{max}} P(z_s| SL) P(obs|z_s, SL) dz_s \\
\end{split} \tag{4}
\end{equation}
$$

$P(obs|z_s, SL)$ is the probability of observing a lensed GW event given that the source is at red-shift $z_s$ and the signal is strongly lensed. Here, it is considered marginalized over all the gravitational wave parameters $\theta$ (without $z_s$), lens related parameters $\theta_L$, lens redshifts $z_L$ and possible source positions $\beta$, i.e., $P(obs|z_s,SL) = \int_{\theta, \theta_L, z_l, \beta} d\theta\, d\theta_L\, dz_l\, d\beta\, P(obs|z_s, \theta, \theta_L, z_l, \beta, SL)$.

$P(SL|z_s)$ is the probability of strong lensing given that the source is at red-shift $z_s$. It is also called the optical depth of strong lensing. It is discussed in detail in the next section.

In Eqn(4) the merger rate density becomes the prior distribution $P(z_s|SL)$ for the lensed source redshift, and it is normalized over the red-shift range. So the normalizing factor is given by,

$$
\begin{equation}
\begin{split}
\mathcal{N}^L = \int_{z_{min}}^{z_{max}} R(z_s) P(SL|z_s) \frac{1}{(1+z_s)} \frac{dV_c}{dz_s} dz_s 
\end{split} \tag{5}
\end{equation}
$$

Therefore,

$$
\begin{equation}
\begin{split}
P(z_s|SL) = \frac{R(z_s)}{\mathcal{N}^L} P(SL|z_s) \frac{1}{(1+z_s)} \frac{dV_c}{dz_s}
\end{split} \tag{6}
\end{equation}
$$

This is done because we want to sample the red-shift values from this prior distribution, which later will be used in the monte-carlo integration to estimate the rate. 

Below I have included the Bayesian way of arriving at $P(z_s|SL)$ expression.

$$
\begin{equation}
\begin{split}
P(z_s|SL) &= \frac{P(SL|z_s) P(z_s)}{P(SL)} \\
&= \frac{P(SL|z_s) P(z_s)}{\int_{z_{min}}^{z_{max}} P(SL|z_s) P(z_s) dz_s} \\
&= \frac{P(SL|z_s) P(z_s)}{\mathcal{N}^L/\mathcal{N}^U} \\
&= \frac{R(z_s)}{\mathcal{N}^L} P(SL|z_s) \frac{1}{(1+z_s)} \frac{dV_c}{dz_s}
\end{split} \tag{7}
\end{equation}
$$

In the above expression, I have used the redshift prior $P(z_s)$ of unlensed events as given in Eqn(6) of [Analytical formulation of the un-lensed event sampling](https://ler.readthedocs.io/en/latest/GW_equations.html#Gravitational-wave-events-sampling-(Analytical-formulation)). 

$$
\begin{equation}
\begin{split}
P(z_s) = \frac{R(z_s)}{\mathcal{N}^U} \frac{1}{(1+z_s)} \frac{dV_c}{dz_s}
\end{split} \tag{8.1}
\end{equation}
$$

Let's re-write the integrand in rate expression in Eqn(4) so that it also become dependent on  lens parameter $\theta_L$, lens redshift $z_l$, source position $\beta$, and othet GW parameters $\theta$.

$$
\begin{equation}
\begin{split}
\mathcal{R}_L &= \mathcal{N}^L \int_{z_{min}}^{z_{max}} P(z_s| SL)\, P(obs|z_s, SL)\, dz_s \\ \nonumber
\text{where,}& \\
P(obs|z_s, SL) &= \int_{\theta, \theta_L, z_l, \beta} P(obs|z_s, \theta, \theta_L, z_l, \beta, SL) \\
&\;\; P(\theta, \theta_L, z_l, \beta|z_s, SL) d\theta\, d\theta_L\, dz_l\, d\beta \\
\end{split}
\end{equation}
$$

Let's break down the integrand and reassemble it, so that it is easier from the sampling perspective.

$$
\begin{equation}
\begin{split}
P(\theta, \theta_L, z_l, \beta|z_s, SL) P(z_s| SL) &= P(\theta) P(\theta_L, z_l, \beta|z_s, SL) P(z_s| SL) \\
&= P(\theta) P(\beta|\theta_L, z_l, z_s, SL) P(\theta_L, z_l|z_s, SL)P(z_s| SL) \\
&= P(\theta) P(\beta|\theta_L, z_l, z_s, SL) P(\theta_L, z_l, z_s|SL) \\
\end{split} \tag{8.2}
\end{equation}
$$

$$
\begin{equation}
\begin{split}
P(\theta_L, z_l, z_s|SL) &= P(SL|\theta_L, z_l, z_s) P(\theta_L, z_l, z_s) \\
\end{split} \tag{8.3}
\end{equation}
$$

$P(SL|\theta_L, z_l, z_s)$ is the cross-section of the lensing event. In `ler`, it is assumed to the double-caustic area and solving it with `lenstronomy` is expensive. The parameters can be sample from intrinsic distribution $P(\theta_L, z_l, z_s)$ and reject or select it based on the lensing cross-section. This is ineffecient and slow as most of the samples will be rejected due to small lensing cross-section. Thus, we re-strucutre the mathematics to have a more efficient sampling.

Consider lens parameters $\theta_L\in \{\sigma, q, \psi, \gamma, [\gamma_1,\gamma_2]\}$ and the lensing cross-section, $\phi_{\text{SIE}} = P(SL|q,\sigma,z_l,z_s)$. Let's first seperate the co-related parameters from the uncorrelated one from Eqn(8.2).

$$
\begin{equation}
\begin{split}
P(\sigma, q, \psi, \gamma, \gamma_1, \gamma_2, z_l, z_s|SL) &= P(\psi|SL) P(\gamma|SL) P(\gamma_1, \gamma_2|SL) P(\sigma, q, z_l, z_s|SL) \\
P(\sigma, q, z_l, z_s|SL) &=  P(\sigma, q | z_l, z_s, SL) P( z_l| z_s, SL) P( z_s| SL)\\
P(\sigma, q | z_l, z_s, SL) &\propto P(SL| \sigma, q, z_l, z_s) P(\sigma, q | z_l, z_s) = \phi_{\text{SIE}} P(q| \sigma, z_l, z_s) P(\sigma | z_l, z_s) \\
P( z_s| SL) &\propto P(SL| z_s) P(z_s)\\
\end{split} \tag{8.3}
\end{equation}
$$

The sampling is done in the following order,

* sample $z_s$ from astrophysical prior $P(z_s)$ and reject or select it based on the optical-depth $P(SL|z_s)$.

* given $z_s$, sample $z_l$ from $P(z_l|z_s, SL)$. $P(z_l|z_s, SL)$ derivation is discussed in the next section.

* given $z_s$ and $z_l$, sample $\sigma$ from an astrophysical distribution $P(\sigma|z_l, z_s)$ and then sample $q$ from $P(q| \sigma, z_l, z_s)$ Collett et al. 2015. The sampled $\sigma$ $q$ pair is then accepted or rejected based on the lensing cross-section $\phi_{\text{SIE}}$.

* Other lensing paramter distributions are pre-computed using equation (8.3), under strong lensing condition. So the sampling of $\psi$, $\gamma$, $[\gamma_1,\gamma_2]$ is done directly from the pre-computed distribution. This is possible because they are not correlated with the $\sigma$, $q$, $z_l$ and $z_s$.

$$
\begin{equation}
\begin{split}
\mathcal{R}_L &= \mathcal{N}^L \int P(z_s| SL)\, P(obs|z_s, \theta, \theta_L, z_l, \beta, SL) \\
&\;\; P(\theta)\, P(\theta_L|z_l, z_s, SL)\, P(z_l|z_s)\, P(\beta|z_s,z_l, \theta_L, SL) dz_s d\theta d\theta_L dz_l d\beta \\
\end{split} \tag{9}
\end{equation}
$$

$$
\begin{equation}
\begin{split}
\mathcal{R}_L &=  \mathcal{N}^L \int \frac{R(z_s)}{\mathcal{N}^L} P(SL|z_s) \frac{1}{(1+z_s)} \frac{dV_c}{dz_s}\, P(obs|z_s, \theta, \theta_L, z_l, \beta, SL)  \\
&\;\; P(\theta)\, P(\theta_L|z_l, z_s, SL)\, P(z_l|z_s, SL)\, P(\beta|z_s,z_l, \theta_L, SL) dz_s d\theta d\theta_L dz_l d\beta \\
\end{split} \tag{10}
\end{equation}
$$

Eqn(10) is the complete expression for the rate of lensed detectable events.

In arriving at Eqn(9) I have used the following relations, 

$$
\begin{equation}
\begin{split}
P(z_s|SL)\, P(obs|z_s, SL) &= \int_{z_l}\int_{\beta}\int_{\theta}\int_{\theta_L} P(z_s, z_l, \beta, \theta, \theta_L|SL)\, P(obs|z_s, \theta, \theta_L, z_l, \beta, SL) d\theta d\theta_L dz_l d\beta \nonumber \\
P(z_s, z_l, \beta, \theta, \theta_L|SL) &= P(\theta) P(z_s, z_l, \beta, \theta_L|SL) \\
P(z_s, z_l, \beta, \theta_L|SL) &= P(\beta|z_s,z_l, \theta_L, SL) P(z_s, z_l, \theta_L|SL) \\
P(z_s, z_l, \theta_L|SL) &= P(z_l, \theta_L|z_s, SL) P(z_s|SL) \\
P(z_l, \theta_L|z_s, SL) &= P(\theta_L|z_s, z_l, SL) P(z_l|z_s, SL) \\
&= P(\theta_L|z_l, z_s, SL) P(z_l|z_s, SL) \\
\end{split}
\end{equation}
$$

Therefore,

$$
\begin{equation}
\begin{split}
&P(z_s|SL)\, P(obs|z_s, SL) = \nonumber \\
&\;\;\; \int_{z_l}\int_{\beta}\int_{\theta}\int_{\theta_L} P(obs|z_s, \theta, \theta_L, z_l, \beta, SL) \nonumber \\
&\;\;\; P(z_s|SL) P(\theta) P(\theta_L|z_l, z_s, SL)\, P(z_l|z_s, SL)\, P(\beta|z_s,z_l, \theta_L, SL) d\theta d\theta_L dz_l d\beta 
\end{split}
\end{equation}
$$

## Optical depth

**Strong lensing probability (where does it come from?)**

### General formulation

Optical depth $P(SL|z_s)$ is defined as the probability of a source being strongly lensed. One can also think of it as cross-section of the lensing, i.e..

$$
\begin{equation}
\begin{split}
&\text{Optical-depth of a source at redshift }z_s\\ 
&\;\;=\frac{\text{Sum of small portions of the sky where lensing can happen}}{\text{Entire sky}}
\end{split} \tag{9}
\end{equation}
$$

The numerator corresponds to relevant lensing area of all possible lensing galaxies present in the universe between the observer and the source. For a given source at red-shift $z_s$, the probability of strong lensing is given by $P(SL|z_s)$. And let there be a distribution of lenses $dN(z_l)$ at red-shift $z_l$ in a spherical shell of thickness $d z_l$. Then the probability of strong lensing is given by,

$$
\begin{equation}
\begin{split}
P(SL|z_s) &=  \frac{1}{4\pi} \int^{z_s}_{z_l=0} P(SL|z_s, z_l) dN(z_l) \\
\end{split} \tag{10}
\end{equation}
$$

Let's try to write the integrand in terms of number density. Consequently, we need to consider the integrand's dependence on lens parameters $\theta_L$ and lens redshift $z_l$. Now the Eqn(10) becomes,

$$
\begin{equation}
\begin{split}
P(SL|z_s) &=  \frac{1}{4\pi} \int^{z_s}_{z_l=0} \int_{\theta_L} P(SL|z_s, z_l, \theta_L) \frac{d^2N(z_l, \theta_L)}{dV_c d\theta_L} \frac{dV_c}{dz_l} dz_l d\theta_L \\
\end{split} \tag{11}
\end{equation}
$$

The number density of lenses $d^3N(z_l)/dV_c dz_ld\theta_L$ is the number of lenses per unit volume per unit red-shift per unit lens parameter. The lensing probability $P(SL|z_s, z_l, \theta_L)$ is the probability of strong lensing given that the source is at red-shift $z_s$, the lens is at red-shift $z_l$ and the lens parameters are $\theta_L$. The integral is over all possible lens red-shifts and lens parameters. The factor $1/4\pi$ is the solid angle of the entire sky.

### Optical depth for Singular isothermal sphere SIS lens

This is in reference to [Haris et al. 2018](https://arxiv.org/abs/1807.07062). 

Take $P(SL|z_s, z_l, \theta_L)$ as $\phi_{\text{SIS}}=\pi \theta_E^2$ ($\theta_E$: Einstein radius) and $\theta_L\in \{\sigma\}$, where $\sigma$ is the velocity dispersion of the lens. The probability of strong lensing for a SIS lens is given by,

$$
\begin{equation}
\begin{split}
P(SL|z_s) &=  \frac{1}{4\pi} \int^{z_s}_{z_l=0} \int_{\sigma} P(SL|z_s, z_l, \sigma) \frac{d^2N(z_l, \sigma)}{dV_c d\sigma} \frac{dV_c}{dz_l} dz_l d\sigma \\
&=  \frac{1}{4\pi} \int^{z_s}_{z_l=0} \int_{\sigma} \phi_{\text{SIS}} \frac{d^2N(z_l, \sigma)}{dV_c d\sigma} \frac{dV_c}{dz_l} dz_l d\sigma \\
\end{split} \tag{12a}
\end{equation}
$$

Haris, following Choi et al. 2008 for early type galaxy, has considered the number density of the lens, $\frac{d^2N(z_l, \sigma)}{dV_c d\sigma}=\,<n>_{\sigma\in P(\sigma)} P(\sigma)$ where $P(\sigma)$ is the PDF of velocity dispersion and $<n>_{\sigma\in P(\sigma)}=n_o=8\times 10^{-3} h^3 Mpc^{-3}$ is the average number density of the lens which constant over the red-shift range of local universe. But I will consider a general case, like in Oguri et al. 2018 for all-type galaxy, where $\frac{d^2N(z_l, \sigma)}{dV_c d\sigma}=\phi(\sigma,z_l)=\phi_\text{loc}(\sigma)\frac{\phi_\text{hyd}(\sigma,z_l)}{\phi_\text{hyd}(\sigma,0)}$, and has both function dependence on both red-shift and velocity dispersion. The optical depth becomes,

$$
\begin{equation}
\begin{split}
P(SL|z_s) &=  \int \frac{\pi \theta_E^2}{4 \pi} \phi(\sigma,z_l) \frac{dV_c}{dz_l} dz_l d\sigma = \int_0^{z_s} \Phi_{SIS}(z_l| z_s) dz_l \\
\end{split} \tag{12b}
\end{equation}
$$

Considering $P_o(z_l)$ and $P_o(\sigma)$ are the uniform distribution over the respective parameters, the integral can be numerically calculated as calculated as,

$$
\begin{equation}
\begin{split}
P(SL|z_s) &=  \int \frac{\theta_E^2}{4} \phi(\sigma,z_l) \frac{dV_c}{dz_l} \frac{P_o(z_l)}{P_o(z_l)} \frac{P_o(\sigma)}{P_o(\sigma)} dz_l d\sigma  \\ \nonumber
&=  (z_{l,max}-z_{l,min})(\sigma_{max}-\sigma_{min})\left< \frac{\theta_E^2}{4} \phi(\sigma,z_l) \frac{dV_c}{dz_l}\right>_{\sigma\in P_o(\sigma)\, z_l\in P_o(z_l)} \\ \nonumber
&=  \Delta z_l \Delta \sigma \left< \frac{\theta_E^2}{4} \phi(\sigma,z_l) \frac{dV_c}{dz_l}\right>_{\sigma\in P_o(\sigma)\, z_l\in P_o(z_l)}
\end{split} \tag{12c}
\end{equation}
$$

Or we can first compute $\Phi_{SIS}(z_l| z_s)$ and then integrate over $z_l$ to get $P(\text{SL}|z_s)$. `ler` prefer the second approach, as it allows us to compute $P(z_l|z_s, \text{SL})$ as follows,

$$
\begin{equation}
\begin{split}
P(z_l|z_s, SL) &= \frac{\Phi_{SIS}(z_l| z_s)}{\kappa(z_s)}= \frac{\Delta \sigma}{\kappa(z_s)} \left< \frac{\theta_E^2}{4} \phi(\sigma,z_l) \frac{dV_c}{dz_l}\right>_{\sigma\in P_o(\sigma)}
\end{split} \tag{12d}
\end{equation}
$$

Where $\kappa(z_s)$ is the normalization factor given the source redshift $z_s$.

### Optical depth for SIE lens

Take $P(SL|z_s, z_l, \theta_L)$ as $\phi_{\text{SIE}}=\phi_{\text{SIS}} \phi_{\text{CUT}}^{\text{SIE}}(q)$ and $\theta_L\in \{\sigma, q\}$, where $\sigma$ is the velocity dispersion and $q$ is the axis-ratio. $\phi^{SIE}_{CUT}(q)$ is derive from the expression given in [Fei Xu et al. 2022](https://iopscience.iop.org/article/10.3847/1538-4357/ac58f8). Let's start with Eqn(12a).

$$
\begin{equation}
\begin{split}
P(SL|z_s) &=  \frac{1}{4\pi} \int^{z_s}_{z_l=0} \int_{\sigma} P(SL|z_s, z_l, \sigma) \frac{d^2N(z_l, \sigma)}{dV_c d\sigma} \frac{dV_c}{dz_l} dz_l d\sigma \\
&=  \frac{1}{4\pi} \int^{z_s}_{z_l=0} \int_{\sigma} \int_{q} P(SL|z_s, z_l, \sigma, q) P(q) \frac{d^2N(z_l, \sigma)}{dV_c d\sigma} \frac{dV_c}{dz_l} dz_l d\sigma dq \\
&=  \frac{1}{4\pi} \int^{z_s}_{z_l=0} \int_{\sigma} \int_{q} \phi_{\text{SIE}} P(q) \phi_(\sigma, z_l) \frac{dV_c}{dz_l} dz_l d\sigma dq \\
&= \frac{\Delta z_l \Delta \sigma}{4\pi} \left< \phi_{\text{SIE}}\; \phi_(\sigma, z_l)\; \frac{dV_c}{dz_l}\right>_{\sigma\in P_o(\sigma)\, z_l\in P_o(z_l)\, q\in P(q)}
\end{split} \tag{12a}
\end{equation}
$$

### Optical depth for EPL+Shear lens

Take $P(SL|z_s, z_l, \theta_L)$ as $\phi_{\text{EPL}}$ and $\theta_L \in$ {$\sigma$ (velocity-dispersion), $q$ (axis-ratio), $\psi$ (axis-rotation), $\gamma$ (mass density spectral-index), $[\gamma_1,\gamma_2]$ (external-shear), $[e_1,e_2]$ (ellipticity)}.

$$
\begin{equation}
\begin{split}
P(SL|z_s) &= \frac{\Delta z_l \Delta \sigma}{4\pi} \left< \phi_{\text{EPL}}\; \phi_(\sigma, z_l)\; \frac{dV_c}{dz_l}\right>_{\sigma\in P_o(\sigma)\, z_l\in P_o(z_l)\, q\in P(q) \, \psi\in P(\psi) \, \gamma\in P(\gamma) \, \gamma_1,\gamma_2\in P(\gamma_1,\gamma_2) \, e_1,e_2\in P(e_1,e_2)}
\end{split} \tag{12a}
\end{equation}
$$


where P_o(\sigma) and P_o(z_l) are the uniform distribution, and $P(q)$, $P(\psi)$, $P(\gamma)$, $P(\gamma_1,\gamma_2)$, $P(e_1,e_2)$ are the astrophysical PDF of the respective lens parameters.

## The order of sampling and rate calculation steps in LeR are listed below.

1. Sample $z_s$ from $P(z_s)$. And apply rejection sampling with optical depth, $P(SL|z_s)$. Other GW parameters are sampled separately, $P(\theta)$.
2. $z_l$ from $P(z_l|z_s)$.
3. $\sigma$ from $P(\sigma|z_l)$.
4. $q$ from $P(q|\sigma)$.
5. Calculation of Einstein radius and application of lensing condition to the sampled lens parameters, $P(SL|z_s, z_l, \sigma, q) \propto \theta_E^2\,\phi^{SIE}_{CUT}$.
6. Other lens parameters ($e_1$, $e_2$, $\gamma_1$, $\gamma_2$, $\gamma$) are sampled independent of the SL condition, $P(e_1,e_2,\gamma_1,\gamma_2,\gamma)$. But, this will be rejection sampled later along with the image position.
7. Draw image position, $\beta$, from within the caustic boundary and solve the lens equation to get image positions. Accept it if it results in 2 or more images, otherwise resample $\beta$. Note that the source postion and image positions are wrt to the center of the lens galaxy (thin lens approximation). But in the sky map the source position will be sampled ra and dec.
8. Sometimes (once in 100-200 thousand), the strong lensing condition cannot be satisfied. For these particular events, resample lens parameters and draw image positions, i.e. repeat steps 2-7.
9. Calculate the magnification, $\mu_i$, time-delay, $dt_i$ and morse phase, $n_i$ for each of the lensed event.
10. Modify the luminosity distance, $D_l$ to $D_l^{eff}=D_l/\sqrt{|\mu_i|}$, and geocent_time to $t_{eff}=t_{geocent}+dt_i$. 
11. Calculate SNR with [gwsnr](https://gwsnr.readthedocs.io/en/latest/)
12. Apply the SNR threshold and check whether the event is detectable or not.
13. Calculate rate of lensed events, $\mathcal{R}_L$ using equation 13.
