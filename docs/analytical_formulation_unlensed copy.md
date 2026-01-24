# Analytical Formulation for Gravitational Wave Event Rates

## Table of Contents

1. [Introduction](#introduction)
2. [Parameter-Marginalized Event Rate](#parameter-marginalized-event-rate)
3. [Redshift Distribution and Intrinsic Merger Rates](#redshift-distribution-and-intrinsic-merger-rates)
4. [Detection Criterion and SNR Modeling](#detection-criterion-and-snr-modeling)
5. [GW Source Parameter Priors](#gw-source-parameter-priors)

## Introduction

The annual rate of detectable (unlensed) gravitational-wave (GW) events, $\frac{\Delta N^{\rm obs}_{\rm U}}{\Delta t}$, gives the expected number of observed compact-binary coalescences (BBH, BNS, and NSBH) per year for a given detector or detector network. It is obtained by combining the total intrinsic merger rate in the detector frame, $\frac{\Delta N_{\rm U}}{\Delta t}$, with the population-averaged probability of detection, $P({\rm obs})$, such that

$$
\begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm U}}{\Delta t}
= \frac{\Delta N_{\rm U}}{\Delta t} \times P({\rm obs}) \,.
\end{split}
\end{equation}
$$

## Parameter-Marginalized Event Rate

The detection probability $P({\rm obs})$ is obtained by marginalizing over the relevant GW source parameters $\vec{\theta}$, i.e.,

$$
\vec{\theta}\in \{z_s, m_1, m_2, a_1, a_2, \theta_1, \theta_2, \phi_{12}, \phi_{JL}, \iota, \phi, \psi, {\rm RA}, {\rm Dec}, t_c\},
$$

where $z_s$ is the source redshift, $(m_1,m_2)$ are component masses, $(a_1,a_2)$ are dimensionless spin magnitudes, $(\theta_1,\theta_2,\phi_{12},\phi_{JL})$ describe spin orientations, $\iota$ is the inclination angle, $(\phi,\psi)$ are the coalescence phase and polarization angle, $({\rm RA},{\rm Dec})$ specify sky position, and $t_c$ is the coalescence time. The prior distributions used in this analysis are summarized in [Table 1](#gw-source-parameter-priors).

Making the parameter dependence explicit, the observed event rate can be written as

$$
\begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm U}}{\Delta t}
= \frac{\Delta N_{\rm U}}{\Delta t} \int_{\vec{\theta}} P({\rm obs}|\vec{\theta}) \, P(\vec{\theta}) \, d\vec{\theta} \,,
\end{split}
\end{equation}
$$

where $P({\rm obs}|\vec{\theta})$ is the conditional probability of detection for a source with parameters $\vec{\theta}$, and $P(\vec{\theta})$ is the joint prior distribution. In most applications, the joint prior is constructed from standard assumptions such as isotropic sky location and orientation, along with independent or dependent priors on intrinsic parameters.

## Redshift Distribution and Intrinsic Merger Rates

The `ler` package uses redshift $z_s$ to represent source distance rather than luminosity distance $D_L$, and assumes that $z_s$ is uncorrelated with the other source parameters. The redshift probability density $P(z_s)$ is defined as the normalized distribution of sources over cosmic history and is proportional to the detector frame intrinsic merger rate density $\frac{d^2 N}{dt \, dV_c}$ and the comoving volume element $\frac{dV_c}{dz_s}$. Equivalently, using the source-frame intrinsic merger rate density $R_{\rm U}(z_s)=\frac{d^2 N}{d\tau , dV_c}$, it can be written as

$$
\begin{equation}
\begin{split}
P(z_s) \propto \frac{1}{(1+z_s)} R_{\rm U}(z_s) \frac{dV_c}{dz_s} \,,
\end{split}
\end{equation}
$$

where $R_{\rm U}(z_s)$ is expressed per unit source-frame proper time $\tau$ and per unit comoving volume, with $d\tau = dt/(1+z_s)$. The factor $1/(1+z_s)$ accounts for cosmological time dilation between the source-frame time $\tau$ and the detector-frame time $t$. The term $\frac{dV_c}{dz_s}dz_s$ represents the comoving shell volume element at redshift $z_s$, and the integration over $z_s$ in the event-rate calculation is carried out over the full redshift range of interest.

Normalizing the redshift distribution introduces the constant ${\cal N}_{\rm U}$, which is equal to the total intrinsic merger rate per year in the detector frame. The normalized form becomes

$$
\begin{equation}
\begin{split}
P(z_s) = \frac{1}{{\cal N}_{\rm U}} \frac{R_{\rm U}(z_s)}{1+z_s} \frac{dV_c}{dz_s} \,,
\end{split}
\end{equation}
$$

with

$$
\begin{equation}
\begin{split}
{\cal N}_{\rm U} = \int_{z_s} \frac{R_{\rm U}(z_s)}{1+z_s} \frac{dV_c}{dz_s} \, dz_s \,.
\end{split}
\end{equation}
$$

Using this normalization, the observed event rate can be expressed as the intrinsic detector-frame rate multiplied by the expectation value of the detection probability over the prior distribution,

$$
\begin{equation}
\begin{split}
\frac{\Delta N^{\rm obs}_{\rm U}}{\Delta t}
&= {\cal N}_{\rm U} \int_{\vec{\theta}} P({\rm obs}|\vec{\theta}) \, P(\vec{\theta}) \, d\vec{\theta} \,, \\
&= {\cal N}_{\rm U} \bigg\langle P({\rm obs}| \vec{\theta}) \bigg\rangle_{\vec{\theta} \sim P(\vec{\theta})} \,,
\end{split}
\end{equation}
$$

which is evaluated numerically using Monte Carlo integration by drawing samples from $P(\vec{\theta})$.

## Detection Criterion and SNR Modeling

The conditional detection probability $P({\rm obs}|\vec{\theta})$ is determined by a detection threshold $\rho_{\rm th}$ on the observed signal-to-noise ratio (SNR) $\rho_{\rm obs}$. In the simplest step-function model, the detection probability is

$$
\begin{equation}
\begin{split}
P({\rm obs}| \vec{\theta}) = P_{\rm det} (\vec{\theta}, \rho_{\rm th}) = 
\begin{cases}
1, & \text{if } \rho_{\rm obs}(\vec{\theta}) > \rho_{\rm th} \\
0, & \text{otherwise},
\end{cases}
\end{split}
\end{equation}
$$

where $P_{\rm det}$ denotes the probability of detection under the chosen detection criterion. In `ler`, $P_{\rm det}$ is estimated using the [`gwsnr`](https://gwsnr.hemantaph.com) package. The observed SNR $\rho_{\rm obs}$ is modeled either as a Gaussian random variate centered at $\rho_{\rm opt}$ (or $\rho_{\rm opt,net}$ for a detector network) with unit variance (Fishbach et al. 2020; Abbott et al. 2019), or as a non-central $\chi$ distribution (`scipy.stats.ncx2`) with non-centrality parameter $\lambda = \rho_{\rm opt}$ (or $\rho_{\rm opt,net}$) and two degrees of freedom for a single detector, extended to $2N$ for a network of $N$ detectors (Essick 2023). Refer to the `gwsnr` documentation for more details on SNR modeling and detection probability calculations.

## GW Source Parameter Priors

The prior distributions and parameter ranges used in the unlensed event rate calculation are summarized in Table 1. These choices follow standard assumptions used in GW population analyses and are consistent with GWTC-3 motivated population models.

### Table 1: Prior distributions for GW source parameters used in the unlensed event rate calculations.

| Parameter | Unit | Prior <br>Distribution | Range <br>[Min, Max] | Description |
| :--- | :--- | :--- | :--- | :--- |
| $m_{1,2}$ | $M_\odot$ | PowerLaw<br>+Peak | [4.98, 112.5] | Component masses ($\alpha$=3.78, $\mu_g$=32.27, $\sigma_g$=3.88, $\lambda_p$=0.03, $\delta_m$=4.8, $\beta$=0.81) |
| $a_{1,2}$ | - | Uniform | [0, 0.99] | Dimensionless spin magnitudes |
| $\theta_{1,2}$ | rad | Sine | [0, $\pi$] | Spin vector tilt angles |
| $\phi_{12}$ | rad | Uniform | [0, $2\pi$] | Azimuthal angle between spin vectors |
| $\phi_{JL}$ | rad | Uniform | [0, $2\pi$] | Angle between total and orbital angular momentum |
| RA | rad | Uniform | [0, $2\pi$] | Right ascension |
| Dec | rad | Cosine | [$-\frac{\pi}{2}$, $\frac{\pi}{2}$] | Declination |
| $\iota$ | rad | Sine | [0, $\pi$] | Inclination relative to line of sight |
| $\psi$ | rad | Uniform | [0, $\pi$] | GW signal polarization angle |
| $\phi$ | rad | Uniform | [0, $2\pi$] | Coalescence orbital phase |
