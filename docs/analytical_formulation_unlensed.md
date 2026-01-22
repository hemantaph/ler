# Analytical Formulation for Gravitational Wave Event Rates

The annual rate of detectable unlensed gravitational wave (GW) events, $\frac{\Delta N^{\rm obs}_{\rm U}}{\Delta t}$, can be determined by integrating the intrinsic rate of compact binary mergers across cosmic time and volume, weighted by the probability of detecting these events. This relationship can be expressed by first considering the total intrinsic merger rate per year, $\frac{\Delta N_{\rm U}}{\Delta t}$, and the average detection probability, $P({\rm obs})$, and it reads,

$$
\begin{equation}
\frac{\Delta N^{\rm obs}_{\rm U}}{\Delta t} = \frac{\Delta N_{\rm U}}{\Delta t} \times P({\rm obs}) \tag{1}
\end{equation}
$$

The detection probability $P({\rm obs})$ is effectively marginalized over all relevant GW source parameters $\vec{\theta}$, i.e.\,

$$\vec{\theta}\in \{z_s, m_1, m_2, a_1, a_2, \theta_1, \theta_2, \phi_{12}, \phi_{JL}, \theta_{\rm jn}, \phi, \psi, {\rm RA, Dec}, t_c\}$$

Refer to [Table 1](#table-1-prior-distributions-for-gw-source-parameters-used-in-the-unlensed-event-rate-calculations) for the complete list of GW parameters and their prior distributions used in this analysis.

To account explicitly for the parameters dependence, the observed event rate can be written as an integral over $\vec{\theta}$, as,

$$
\begin{equation}
\frac{\Delta N^{\rm obs}_{\rm U}}{\Delta t} = \frac{\Delta N_{\rm U}}{\Delta t} \int_{\vec{\theta}} P({\rm obs}|\vec{\theta}) \, P(\vec{\theta}) \, d\vec{\theta} \tag{2}
\end{equation}
$$

where $P({\rm obs}|\vec{\theta})$ is the conditional probability of detecting a source given its parameters $\vec{\theta}$, and $P(\vec{\theta})$ is the joint prior distribution of the source parameters. $P(\vec{\theta})$ can be expanded to it's constituent parts, starting with the redshift distribution of sources, and it reads,

$$
\begin{equation}
\begin{split}
P(\vec{\theta}) &= P(z_s) \, P(m_1, m_2) \, P(a_1) \, P(a_2) \, P(\theta_1, \theta_2) \, P(\phi_{12}) \, P(\phi_{JL}) \\
& \quad \times P(\theta_{\rm jn}) \, P(\phi) \, P(\psi) \, P({\rm RA}) \, P({\rm Dec}) \, P(t_c)
\end{split}
\end{equation}
$$

`ler` uses redshift ($z_s$) instead of luminosity distance ($D_L$) as one of the source parameters, and we assume that it is uncorrelated with the other parameters. Then, it's probability density distribution $P(z_s)$ is the normalized redshift distribution of sources and it propotionally relates to the differential intrinsic merger rate at detector frame, i.e.,

$$
\begin{equation}
\begin{split}
P(z_s) &\propto \frac{d^2 N}{dt \, dV_c} \frac{dV_c}{dz_s} \\
&\propto \frac{1}{(1+z_s)}\frac{d^2 N}{d\tau \, dV_c} \frac{dV_c}{dz_s} \\
& \propto \frac{1}{(1+z_s)} R_{\rm U}(z_s) \frac{dV_c}{dz_s}
\end{split} \tag{3}
\end{equation}
$$

where $R_{\rm U}(z_s)$ denotes the intrinsic merger rate in the source frame. Note that the intrinsic merger rate is expressed in terms of source-frame proper time $\tau$ (with $d\tau = dt/(1+z_s)$) and comoving volume $V_c$. The factor $1/(1+z_s)$ accounts for time dilation due to cosmic expansion. We use $R_{\rm U}(z_s)$ because most observational studies and theoretical predictions report merger rates in the source frame. The term $\frac{dV_c}{dz_s}dz_s$ represents a spherical-shell volume element in comoving coordinates at redshift $z_s$, with integration over $z_s$ spanning the entire cosmic history.



$$
P(z_s) = \frac{1}{{\cal N}_{\rm U}} \frac{R_{\rm U}(z_s)}{1+z_s} \frac{dV_c}{dz_s}
$$

where $\frac{dV_c}{dz_s}$ is the differential comoving volume element. The factor $1/(1+z_s)$ accounts for cosmic time dilation, converting the rate from the source's proper time to the observer frame. ${\cal N}_{\rm U}$ is the normalization factor representing the total intrinsic merger rate per year in the detector frame, and it reads,

$$
{\cal N}_{\rm U} = \int_{z_s} \frac{R_{\rm U}(z_s)}{1+z_s} \frac{dV_c}{dz_s} \, dz_s = \frac{\Delta N_{\rm U}}{\Delta t}
$$

A more complete formulation incorporates the set of relevant GW source parameters $\theta$, which typically includes component masses $(m_1, m_2)$, luminosity distance $D_L$ (determined by $z_s$), inclination angle $\theta_{\rm jn}$, orbital phase $\phi$, polarization angle $\psi$, sky location (RA, Dec), and coalescence time $t_c$. The observed event rate then becomes,

$$
\frac{\Delta N^{\rm obs}_{\rm U}}{\Delta t} = {\cal N}_{\rm U} \int_{z_s} \int_{\theta} P({\rm obs}|z_s, \theta) \, P(z_s, \theta) \, dz_s \, d\theta
$$

To simplify the calculation and make use of the established redshift-dependent merger rate $R_{\rm U}(z_s)$, we assume that the redshift distribution is independent of the other parameters $\theta$, i.e., $P(z_s, \theta) = P(z_s)\, P(\theta)$. Using this factorization, the integral becomes,

$$
\frac{\Delta N^{\rm obs}_{\rm U}}{\Delta t} = {\cal N}_{\rm U} \int_{z_s} \int_{\theta} P({\rm obs}|z_s, \theta)\, P(z_s)\, P(\theta) dz_s d\theta = {\cal N}_{\rm U} \bigg\langle P({\rm obs}|z_s, \theta) \bigg\rangle_{z_s \in P(z_s), \theta \in P(\theta)}
$$

The final expression indicates a numerical Monte Carlo integration, where $z_s$ is sampled from $P(z_s)$ and $\theta$ is sampled from the prior distributions specified in the table below.

The conditional detection probability $P({\rm obs}|z_s, \theta)$ is defined using a threshold on the optimal signal-to-noise ratio (SNR), $\rho$, of the GW signal, and it reads,

$$
P({\rm obs}|z_s, \theta) = \begin{cases}
1, & \text{if } \rho(z_s, \theta) > \rho_{\rm th} \\
0, & \text{otherwise}
\end{cases}
$$

where $\rho(z_s, \theta)$ is the optimal SNR, which depends on redshift and source parameters. For this analysis, we adopt $\rho_{\rm th} = 8$ as the detection threshold, consistent with standard detectability criteria in the Gaussian noise regime. The SNR $\rho$ is efficiently computed using the `gwsnr` Python package.

### Table 1: Prior distributions for GW source parameters used in the unlensed event rate calculations.

| Parameter | Unit | Prior | Min | Max | Description |
|-----------|------|-------|-----|-----|-------------|
| $m_{1,2}$ | $M_\odot$ | PowerLaw+Peak (GWTC-3: $\alpha$=3.78, $\mu_g$=32.27, $\sigma_g$=3.88, $\lambda_p$=0.03, $\delta_m$=4.8, $\beta$=0.81) | 4.98 | 112.5 | Component masses of the binary |
| $a_{1,2}$ | dimensionless | uniform | 0 | 0.99 | Dimensionless spin magnitudes of the components |
| $\theta_{1,2}$ | dimensionless | sine | 0 | $\pi$ | Tilt angles of the spin vectors |
| $\phi_{12}$ | rad | uniform | 0 | $2\pi$ | Azimuthal angle between the spin vectors |
| $\phi_{JL}$ | rad | uniform | 0 | $2\pi$ | Azimuthal angle between total angular momentum and orbital angular momentum |
| RA | rad | uniform | 0 | $2\pi$ | Right ascension sky coordinate |
| Dec | rad | cos | $-\pi/2$ | $\pi/2$ | Declination sky coordinate |
| $\theta_{jn}$ | rad | sin | 0 | $\pi$ | Inclination angle relative to line of sight |
| $\psi$ | rad | uniform | 0 | $\pi$ | Polarization angle of the GW signal |
| $\phi$ | rad | uniform | 0 | $2\pi$ | Orbital phase at coalescence |


