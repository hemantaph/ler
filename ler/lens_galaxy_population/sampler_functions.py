# -*- coding: utf-8 -*-
"""
Module for lens galaxy parameter sampling functions.

This module provides probability density functions (PDFs) and random variable
samplers (RVS) for lens galaxy parameters including redshift, velocity dispersion,
and axis ratio. It also includes rejection and importance sampling algorithms
for sampling lens parameters weighted by gravitational lensing cross sections.

Key Components: \n
- Lens redshift samplers (SIS model from Haris et al. 2018) \n
- Velocity dispersion samplers (generalized gamma distribution) \n
- Axis ratio samplers (Rayleigh and Padilla-Strauss distributions) \n
- Rejection and importance sampling for cross-section weighted parameters \n

Copyright (C) 2026 Hemantakumar Phurailatpam. Distributed under MIT License.
"""

import numpy as np
from numba import njit, prange
from scipy.interpolate import CubicSpline
from astropy.cosmology import LambdaCDM
from ..utils import is_njitted
from multiprocessing import Pool
from tqdm import tqdm
from ..utils import (
    save_pickle,
    load_pickle,
    inverse_transform_sampler,
    redshift_optimal_spacing,
)
import os


def available_sampler_list():
    """
    Return list of available lens parameter samplers.

    Returns
    -------
    sampler_list : ``list``
        List of available sampler function names.

    Examples
    --------
    >>> samplers = available_sampler_list()
    >>> print(samplers)
    ['lens_redshift_strongly_lensed_sis_haris', 'velocity_dispersion_gengamma', 'axis_ratio_rayleigh', 'axis_ratio_padilla_strauss']
    """
    return [
        "lens_redshift_strongly_lensed_sis_haris_pdf",
        "lens_redshift_strongly_lensed_sis_haris_rvs",
        "velocity_dispersion_ewoud_denisty_function",
        "velocity_dispersion_bernardi_denisty_function",
        "velocity_dispersion_gengamma_density_function",
        "velocity_dispersion_gengamma_pdf",
        "velocity_dispersion_gengamma_rvs",
        "axis_ratio_rayleigh_rvs",
        "axis_ratio_rayleigh_pdf",
        "axis_ratio_padilla_strauss_rvs",
        "axis_ratio_padilla_strauss_pdf",
        "bounded_normal_sample",
        "rejection_sampler",
        "importance_sampler",
        "importance_sampler_mp",
    ]


# ---------------------------------
# Lens redshift sampler functions
# ---------------------------------
@njit
def lens_redshift_strongly_lensed_sis_haris_pdf(
    zl,
    zs,
    cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0),
):
    """
    Compute lens redshift PDF for SIS model (Haris et al. 2018).

    Computes the probability density function for lens redshift between
    zl=0 and zl=zs using the analytical form from Haris et al. (2018)
    equation A7, based on the SIS (Singular Isothermal Sphere) lens model.

    Parameters
    ----------
    zl : ``float``
        Redshift of the lens galaxy.
    zs : ``float``
        Redshift of the source.
    cosmo : ``astropy.cosmology``
        Cosmology object for distance calculations. \n
        default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)

    Returns
    -------
    pdf : ``float``
        Probability density at the given lens redshift.

    Examples
    --------
    >>> from astropy.cosmology import LambdaCDM
    >>> cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)
    >>> pdf = lens_redshift_strongly_lensed_sis_haris_pdf(zl=0.5, zs=1.0, cosmo=cosmo)
    >>> print(f"PDF at zl=0.5: {pdf:.4f}")
    PDF at zl=0.5: 1.8750
    """
    Dc_zl = cosmo.comoving_distance(zl).value
    Dc_zs = cosmo.comoving_distance(zs).value
    x = Dc_zl / Dc_zs
    return 30 * x**2 * (1 - x) ** 2


def lens_redshift_strongly_lensed_sis_haris_rvs(
    size,
    zs,
    z_min=0.001,
    z_max=10.0,
    cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0),
):
    """
    Sample lens redshifts for SIS model (Haris et al. 2018).

    Uses inverse transform sampling with the analytical CDF of the Haris et al.
    (2018) lens redshift distribution to efficiently generate samples.

    Parameters
    ----------
    size : ``int``
        Number of samples to draw.
    zs : ``float``
        Redshift of the source.
    z_min : ``float``
        Minimum redshift for interpolation grid. \n
        default: 0.001
    z_max : ``float``
        Maximum redshift for interpolation grid. \n
        default: 10.0
    cosmo : ``astropy.cosmology``
        Cosmology object for distance calculations. \n
        default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)

    Returns
    -------
    zl : ``numpy.ndarray``
        Array of sampled lens redshifts with shape (size,).

    Examples
    --------
    >>> from astropy.cosmology import LambdaCDM
    >>> cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)
    >>> zl_samples = lens_redshift_strongly_lensed_sis_haris_rvs(
    ...     size=1000,
    ...     zs=1.5,
    ...     z_min=0.001,
    ...     z_max=10.0,
    ...     cosmo=cosmo,
    ... )
    >>> print(f"Mean lens redshift: {zl_samples.mean():.2f}")
    >>> print(f"Redshift range: [{zl_samples.min():.3f}, {zl_samples.max():.3f}]")
    """
    # Create comoving distance to redshift mapping
    zs_array = redshift_optimal_spacing(z_min, z_max, 500)
    Dc_array = cosmo.comoving_distance(zs_array).value
    inverse_spline_Dc = CubicSpline(Dc_array, zs_array)

    # Create CDF values using analytical form: CDF = 6x^5 - 15x^4 + 10x^3
    x_array = np.linspace(0.0, 1.0, 500)
    cdf_values = 6 * x_array**5 - 15 * x_array**4 + 10 * x_array**3

    # Inverse transform sampling
    r = inverse_transform_sampler(size, cdf_values, x_array)

    zs_Dc = cosmo.comoving_distance(zs).value
    zl_Dc = zs_Dc * r
    return inverse_spline_Dc(zl_Dc)


# ---------------------------------------
# Velocity dispersion sampler functions
# ---------------------------------------
@njit
def _gamma(x):
    """
    Compute the gamma function using the Lanczos approximation.

    Parameters
    ----------
    x : ``float``
        Input value.

    Returns
    -------
    result : ``float``
        Gamma function value at x.
    """
    g = 7
    p = np.array(
        [
            0.99999999999980993,
            676.5203681218851,
            -1259.1392167224028,
            771.32342877765313,
            -176.61502916214059,
            12.507343278686905,
            -0.13857109526572012,
            9.9843695780195716e-6,
            1.5056327351493116e-7,
        ]
    )

    if x < 0.5:
        return np.pi / (np.sin(np.pi * x) * _gamma(1 - x))
    else:
        x -= 1
        y = p[0]
        for i in range(1, g + 2):
            y += p[i] / (x + i)
        t = x + g + 0.5
        return np.sqrt(2 * np.pi) * t ** (x + 0.5) * np.exp(-t) * y


@njit
def _cvdf_fit(log_vd, redshift):
    """
    Compute the cumulative velocity dispersion function fit.

    This uses the fit coefficients from Bernardi et al. (2010) for the
    velocity dispersion function at the local universe.

    Parameters
    ----------
    log_vd : ``float``
        Log10 of the velocity dispersion (km/s).
    redshift : ``float``
        Redshift of the lens galaxy.

    Returns
    -------
    result : ``float``
        Cumulative velocity dispersion function value.
    """
    this_vars = np.array(
        [
            [7.39149763, 5.72940031, -1.12055245],
            [-6.86339338, -5.27327109, 1.10411386],
            [2.85208259, 1.25569600, -0.28663846],
            [0.06703215, -0.04868317, 0.00764841],
        ]
    )
    coeffs = [
        this_vars[i][0] + this_vars[i][1] * redshift + this_vars[i][2] * redshift**2
        for i in range(4)
    ]
    mstar = log_vd - coeffs[3]
    return coeffs[0] + coeffs[1] * mstar + coeffs[2] * mstar**2 - np.exp(mstar)


@njit
def _cvdf_derivative(log_vd, redshift, dx):
    """
    Compute the numerical derivative of the CVDF fit function.

    Parameters
    ----------
    log_vd : ``float``
        Log10 of the velocity dispersion (km/s).
    redshift : ``float``
        Redshift of the lens galaxy.
    dx : ``float``
        Step size for numerical differentiation.

    Returns
    -------
    derivative : ``float``
        Numerical derivative value.
    """
    return (
        0.5 * (_cvdf_fit(log_vd + dx, redshift) - _cvdf_fit(log_vd - dx, redshift)) / dx
    )


@njit
def _pdf_phi_z_ratio(sigma, z):
    """
    Compute the PDF ratio of velocity dispersion function at redshift z to z=0.

    This ratio is used in the derivation of the velocity dispersion function
    at redshift z following Oguri et al. (2018b). This lacks the scaling factor.

    Parameters
    ----------
    sigma : ``float``
        Velocity dispersion (km/s).
    z : ``float``
        Redshift of the lens galaxy.

    Returns
    -------
    ratio : ``float``
        Ratio of phi(sigma, z) / phi(sigma, 0).
    """
    log_vd = np.log10(sigma)
    phi_sim_z = 10 ** _cvdf_fit(log_vd, z) / sigma * _cvdf_derivative(log_vd, z, 1e-8)
    phi_sim_0 = 10 ** _cvdf_fit(log_vd, 0) / sigma * _cvdf_derivative(log_vd, 0, 1e-8)

    return phi_sim_z / phi_sim_0


@njit
def velocity_dispersion_ewoud_denisty_function(
    sigma, z, alpha=0.94, beta=1.85, phistar=2.099e-2, sigmastar=113.78
):
    """
    Calculate the lens galaxy velocity dispersion function at redshift z (Oguri et al. (2018b) + Wempe et al. (2022)).

    Parameters
    ----------
    sigma : ``numpy.ndarray``
        Velocity dispersion of the lens galaxy (km/s).
    z : ``float``
        Redshift of the lens galaxy.
    alpha : ``float``
        Shape parameter of the velocity dispersion function. \n
        default: 0.94
    beta : ``float``
        Slope parameter of the velocity dispersion function. \n
        default: 1.85
    phistar : ``float``
        Normalization of the velocity dispersion function (Mpc^-3). \n
        default: 2.099e-2
    sigmastar : ``float``
        Characteristic velocity dispersion (km/s). \n
        default: 113.78

    Returns
    -------
    result : ``numpy.ndarray``
        Velocity dispersion function values. \n
        Negative values are clipped to 0.
    """
    result = _pdf_phi_z_ratio(sigma, z) * velocity_dispersion_bernardi_denisty_function(
        sigma=sigma, alpha=alpha, beta=beta, phistar=phistar, sigmastar=sigmastar
    )
    result[result < 0.0] = 0.0
    return result


@njit
def velocity_dispersion_bernardi_denisty_function(
    sigma, alpha, beta, phistar, sigmastar
):
    """
    Calculate the local universe velocity dispersion function.

    This implements the velocity dispersion function from Bernardi et al. (2010).

    Parameters
    ----------
    sigma : ``numpy.ndarray``
        Velocity dispersion of the lens galaxy (km/s).
    alpha : ``float``
        Shape parameter (alpha/beta is used in the gamma function). \n
        For Oguri et al. (2018b): alpha=0.94 \n
        For Choi et al. (2008): alpha=2.32/2.67
    beta : ``float``
        Slope parameter of the velocity dispersion function. \n
        For Oguri et al. (2018b): beta=1.85 \n
        For Choi et al. (2008): beta=2.67
    phistar : ``float``
        Normalization of the velocity dispersion function (Mpc^-3). \n
        For Oguri et al. (2018b): phistar=2.099e-2*(h/0.7)^3 \n
        For Choi et al. (2008): phistar=8.0e-3*h^3
    sigmastar : ``float``
        Characteristic velocity dispersion (km/s). \n
        For Oguri et al. (2018b): sigmastar=113.78 \n
        For Choi et al. (2008): sigmastar=161.0

    Returns
    -------
    philoc : ``numpy.ndarray``
        Local velocity dispersion function values.
    """
    philoc = (
        phistar
        * (sigma / sigmastar) ** alpha
        * np.exp(-((sigma / sigmastar) ** beta))
        * beta
        / _gamma(alpha / beta)
        / sigma
    )
    return philoc


def velocity_dispersion_gengamma_density_function(
    sigma,
    alpha=0.94,
    beta=1.85,
    phistar=2.099e-2,
    sigmastar=113.78,
    **kwargs,
):
    """
    Compute unnormalized velocity dispersion function using generalized gamma.

    Computes the galaxy velocity dispersion function (VDF) using the generalized
    gamma distribution formulation from Choi et al. (2007).

    Parameters
    ----------
    sigma : ``float`` or ``numpy.ndarray``
        Velocity dispersion in km/s.
    alpha : ``float``
        Power-law index governing the low-velocity slope. \n
        default: 0.94
    beta : ``float``
        Exponential parameter for high-velocity cutoff sharpness. \n
        default: 1.85
    phistar : ``float``
        Normalization constant (comoving number density, Mpc^-3). \n
        default: 8.0e-3
    sigmastar : ``float``
        Characteristic velocity scale in km/s. \n
        default: 113.78

    Returns
    -------
    density : ``float`` or ``numpy.ndarray``
        Unnormalized velocity dispersion function value.

    Examples
    --------
    >>> import numpy as np
    >>> sigma = np.array([150.0, 200.0, 250.0])
    >>> density = velocity_dispersion_gengamma_density_function(sigma)
    >>> print(f"Density at sigma=200 km/s: {density[1]:.6f}")
    """
    from scipy.stats import gengamma

    density = phistar * gengamma.pdf(
        sigma / sigmastar,
        a=alpha / beta,
        c=beta,
    )

    return density


def velocity_dispersion_gengamma_pdf(
    sigma,
    sigma_min=100.0,
    sigma_max=400.0,
    alpha=0.94,
    beta=1.85,
    sigmastar=113.78,
):
    """
    Compute normalized velocity dispersion PDF using generalized gamma.

    Computes the probability density function for velocity dispersion using
    the generalized gamma distribution, normalized over the specified range.

    Parameters
    ----------
    sigma : ``float`` or ``numpy.ndarray``
        Velocity dispersion in km/s.
    sigma_min : ``float``
        Minimum velocity dispersion for normalization (km/s). \n
        default: 100.0
    sigma_max : ``float``
        Maximum velocity dispersion for normalization (km/s). \n
        default: 400.0
    alpha : ``float``
        Power-law index governing the low-velocity slope. \n
        default: 0.94
    beta : ``float``
        Exponential parameter for high-velocity cutoff sharpness. \n
        default: 1.85
    sigmastar : ``float``
        Characteristic velocity scale in km/s. \n
        default: 113.78

    Returns
    -------
    pdf : ``float`` or ``numpy.ndarray``
        Normalized probability density at the given velocity dispersion.

    Examples
    --------
    >>> pdf = velocity_dispersion_gengamma_pdf(
    ...     sigma=200.0,
    ...     sigma_min=100.0,
    ...     sigma_max=400.0,
    ... )
    >>> print(f"PDF at sigma=200 km/s: {pdf:.6f}")
    """
    # Compute normalization constant
    sigma_array = np.linspace(sigma_min, sigma_max, 500)
    density = velocity_dispersion_gengamma_density_function(
        sigma=sigma_array,
        alpha=alpha,
        beta=beta,
        sigmastar=sigmastar,
    )
    norm_const = np.trapz(density, sigma_array)

    pdf = (
        velocity_dispersion_gengamma_density_function(
            sigma=sigma,
            alpha=alpha,
            beta=beta,
            sigmastar=sigmastar,
        )
        / norm_const
    )

    return pdf


def _truncated_gengamma_rvs(size, a, c, loc, scale, lower_bound, upper_bound):
    """
    Sample from truncated generalized gamma distribution via rejection.

    Parameters
    ----------
    size : ``int``
        Number of samples to draw.
    a : ``float``
        Shape parameter of gengamma distribution.
    c : ``float``
        Shape parameter of gengamma distribution.
    loc : ``float``
        Location parameter of gengamma distribution.
    scale : ``float``
        Scale parameter of gengamma distribution.
    lower_bound : ``float``
        Lower bound of the truncated distribution.
    upper_bound : ``float``
        Upper bound of the truncated distribution.

    Returns
    -------
    rvs : ``numpy.ndarray``
        Truncated gengamma random samples with shape (size,).
    """
    from scipy.stats import gengamma

    total_samples = []
    while len(total_samples) < size:
        batch = gengamma.rvs(a, c, loc=loc, scale=scale, size=size * 2)
        valid = batch[(batch >= lower_bound) & (batch <= upper_bound)]
        total_samples.extend(valid)

    return np.array(total_samples[:size])


def velocity_dispersion_gengamma_rvs(
    size,
    sigma_min=100.0,
    sigma_max=400.0,
    alpha=0.94,
    beta=1.85,
    sigmastar=113.78,
):
    """
    Sample velocity dispersions from generalized gamma distribution.

    Uses truncated generalized gamma sampling via rejection to generate
    velocity dispersion samples within the specified bounds.

    Parameters
    ----------
    size : ``int``
        Number of samples to draw.
    sigma_min : ``float``
        Minimum velocity dispersion (km/s). \n
        default: 100.0
    sigma_max : ``float``
        Maximum velocity dispersion (km/s). \n
        default: 400.0
    alpha : ``float``
        Power-law index governing the low-velocity slope. \n
        default: 0.94
    beta : ``float``
        Exponential parameter for high-velocity cutoff sharpness. \n
        default: 1.85
    sigmastar : ``float``
        Characteristic velocity scale in km/s. \n
        default: 113.78

    Returns
    -------
    sigma : ``numpy.ndarray``
        Sampled velocity dispersions in km/s with shape (size,).

    Examples
    --------
    >>> sigma_samples = velocity_dispersion_gengamma(
    ...     size=1000,
    ...     sigma_min=100.0,
    ...     sigma_max=400.0,
    ... )
    >>> print(f"Mean sigma: {sigma_samples.mean():.2f} km/s")
    >>> print(f"Range: [{sigma_samples.min():.1f}, {sigma_samples.max():.1f}] km/s")
    """
    sigma_samples = _truncated_gengamma_rvs(
        size=size,
        a=alpha / beta,
        c=beta,
        loc=0,
        scale=sigmastar,
        lower_bound=sigma_min,
        upper_bound=sigma_max,
    )
    return sigma_samples


# ------------------------------
# Axis ratio sampler functions
# ------------------------------
@njit
def axis_ratio_rayleigh_rvs(size, sigma, q_min=0.2, q_max=1.0):
    """
    Sample axis ratios from velocity-dependent Rayleigh distribution.

    Generates axis ratio samples using the Rayleigh distribution with
    scale parameter dependent on velocity dispersion, as described in
    Wierda et al. (2021) Appendix C.

    Parameters
    ----------
    size : ``int``
        Number of samples to draw.
    sigma : ``numpy.ndarray``
        Velocity dispersions in km/s with shape (size,).
    q_min : ``float``
        Minimum allowed axis ratio. \n
        default: 0.2
    q_max : ``float``
        Maximum allowed axis ratio. \n
        default: 1.0

    Returns
    -------
    q : ``numpy.ndarray``
        Sampled axis ratios with shape (size,).

    Examples
    --------
    >>> import numpy as np
    >>> sigma = np.random.uniform(100, 300, 1000)
    >>> q_samples = axis_ratio_rayleigh_rvs(size=1000, sigma=sigma, q_min=0.2, q_max=1.0)
    >>> print(f"Mean axis ratio: {q_samples.mean():.2f}")
    >>> print(f"Range: [{q_samples.min():.2f}, {q_samples.max():.2f}]")
    """
    a = sigma / 161.0
    q = np.ones(size)
    idx = np.arange(size)
    size_ = size

    while size_ != 0:
        # Scale parameter for Rayleigh distribution (Wierda et al. 2021)
        s = 0.38 - 0.09177 * a[idx]
        s[s <= 0] = 0.0001
        u = np.random.uniform(0, 1, size=size_)
        b = s * np.sqrt(-2 * np.log(u))  # Inverse CDF of Rayleigh
        q_ = 1.0 - b

        # Select samples within bounds
        idx2 = (q_ >= q_min) & (q_ <= q_max)
        q[idx[idx2]] = q_[idx2]

        # Track remaining samples outside bounds
        idx = idx[(q_ <= q_min) | (q_ >= q_max)]
        size_ = len(idx)

    return q


@njit(parallel=True)
def axis_ratio_rayleigh_pdf(q, sigma, q_min=0.2, q_max=1.0):
    """
    Compute truncated Rayleigh PDF for axis ratio.

    Computes the probability density function for axis ratio using the
    truncated Rayleigh distribution with velocity-dependent scale parameter
    (Wierda et al. 2021 equation C16).

    Parameters
    ----------
    q : ``numpy.ndarray``
        Axis ratios at which to evaluate PDF.
    sigma : ``numpy.ndarray``
        Velocity dispersions in km/s (same shape as q).
    q_min : ``float``
        Minimum axis ratio for truncation. \n
        default: 0.2
    q_max : ``float``
        Maximum axis ratio for truncation. \n
        default: 1.0

    Returns
    -------
    pdf : ``numpy.ndarray``
        Probability density values with same shape as q.

    Examples
    --------
    >>> import numpy as np
    >>> q = np.array([0.5, 0.7, 0.9])
    >>> sigma = np.array([150.0, 200.0, 250.0])
    >>> pdf = axis_ratio_rayleigh_pdf(q, sigma)
    >>> print(f"PDF values: {pdf}")
    """
    out = np.zeros_like(q)
    b_lo = 1.0 - q_max
    b_hi = 1.0 - q_min

    for i in prange(q.size):
        s = 0.38 - 0.09177 * (sigma[i] / 161.0)
        if s <= 0.0:
            s = 1e-4

        qi = q[i]
        if qi < q_min or qi > q_max:
            out[i] = 0.0
            continue

        b = 1.0 - qi
        # Base (untruncated) Rayleigh density
        base = (b / (s * s)) * np.exp(-0.5 * (b * b) / (s * s))

        # Truncation normalization
        Z = np.exp(-0.5 * (b_lo * b_lo) / (s * s)) - np.exp(
            -0.5 * (b_hi * b_hi) / (s * s)
        )
        if Z <= 0.0:
            out[i] = 0.0
        else:
            out[i] = base / Z

    return out


@njit
def _axis_ratio_padilla_strauss_data():
    """
    Return axis ratio PDF data points from Padilla & Strauss (2008).

    Returns
    -------
    q_array : ``numpy.ndarray``
        Axis ratio data points.
    pdf_array : ``numpy.ndarray``
        PDF values at each axis ratio.
    """
    q_array = np.array(
        [
            0.04903276402927845,
            0.09210526315789469,
            0.13596491228070173,
            0.20789473684210524,
            0.2899703729522482,
            0.3230132450331126,
            0.35350877192982455,
            0.37946148483792264,
            0.4219298245614036,
            0.4689525967235971,
            0.5075026141512723,
            0.5226472638550018,
            0.5640350877192983,
            0.6096491228070177,
            0.6500000000000001,
            0.6864848379226213,
            0.7377192982456142,
            0.7787295224817011,
            0.8007581038689441,
            0.822786685256187,
            0.8668438480306729,
            0.8973684210526317,
            0.9254385964912283,
        ]
    )
    pdf = np.array(
        [
            0.04185262687135349,
            0.06114520695141845,
            0.096997499638376,
            0.1932510900336828,
            0.39547914337673706,
            0.49569751276216234,
            0.6154609137685201,
            0.7182049959882812,
            0.920153741243567,
            1.1573982157399754,
            1.3353263628106684,
            1.413149656448315,
            1.5790713532948977,
            1.7280185150744938,
            1.8132994441344819,
            1.8365803753840484,
            1.8178662203211204,
            1.748929843583365,
            1.688182592496342,
            1.6274353414093188,
            1.4948487090314488,
            1.402785526832393,
            1.321844068356993,
        ]
    )
    return q_array, pdf


@njit
def axis_ratio_padilla_strauss_rvs(size):
    """
    Sample axis ratios from Padilla & Strauss (2008) distribution.

    Uses inverse transform sampling with the empirical PDF from
    Padilla & Strauss (2008) for early-type galaxy axis ratios.

    Parameters
    ----------
    size : ``int``
        Number of samples to draw.

    Returns
    -------
    q : ``numpy.ndarray``
        Sampled axis ratios with shape (size,).

    Examples
    --------
    >>> q_samples = axis_ratio_padilla_strauss_rvs(size=1000)
    >>> print(f"Mean axis ratio: {q_samples.mean():.2f}")
    >>> print(f"Range: [{q_samples.min():.2f}, {q_samples.max():.2f}]")
    """
    q_array, pdf = _axis_ratio_padilla_strauss_data()

    # Compute CDF and normalize
    cdf_values = np.cumsum(pdf)
    cdf_values = cdf_values / cdf_values[-1]

    return inverse_transform_sampler(size, cdf_values, q_array)


@njit
def axis_ratio_padilla_strauss_pdf(q):
    """
    Compute axis ratio PDF from Padilla & Strauss (2008).

    Evaluates the probability density function for axis ratio using
    cubic spline interpolation of the Padilla & Strauss (2008) data.

    Parameters
    ----------
    q : ``numpy.ndarray``
        Axis ratios at which to evaluate PDF.

    Returns
    -------
    pdf : ``numpy.ndarray``
        Probability density values with same shape as q.

    Examples
    --------
    >>> import numpy as np
    >>> q = np.array([0.3, 0.5, 0.7, 0.9])
    >>> pdf = axis_ratio_padilla_strauss_pdf(q)
    >>> print(f"PDF at q=0.5: {pdf[1]:.4f}")
    """
    q_array, pdf = _axis_ratio_padilla_strauss_data()
    spline = CubicSpline(q_array, pdf, extrapolate=True)
    return spline(q)


@njit
def bounded_normal_sample(size, mean, std, low, high):
    """
    Sample from truncated normal distribution via rejection.

    Generates samples from a normal distribution with specified mean and
    standard deviation, rejecting samples outside the specified bounds.

    Parameters
    ----------
    size : ``int``
        Number of samples to draw.
    mean : ``float``
        Mean of the normal distribution.
    std : ``float``
        Standard deviation of the normal distribution.
    low : ``float``
        Lower bound for samples.
    high : ``float``
        Upper bound for samples.

    Returns
    -------
    samples : ``numpy.ndarray``
        Bounded normal samples with shape (size,).

    Examples
    --------
    >>> samples = bounded_normal_sample(size=1000, mean=2.0, std=0.2, low=1.5, high=2.5)
    >>> print(f"Mean: {samples.mean():.2f}, Std: {samples.std():.2f}")
    >>> print(f"Range: [{samples.min():.2f}, {samples.max():.2f}]")
    """
    samples = np.empty(size)
    for i in range(size):
        while True:
            sample = np.random.normal(mean, std)
            if low <= sample <= high:
                break
        samples[i] = sample

    return samples


# -------------------------------------------------
# Rejection sampling of strongly lensed parameters
# -------------------------------------------------
def rejection_sampler(
    zs,
    zl,
    sigma_max,
    sigma_rvs,
    q_rvs,
    phi_rvs,
    gamma_rvs,
    shear_rvs,
    cross_section,
    safety_factor=1.2,
):
    """
    Core rejection sampling algorithm for lens parameters.

    Parameters
    ----------
    zs : ``numpy.ndarray``
        Source redshifts.
    zl : ``numpy.ndarray``
        Lens redshifts.
    sigma_max : ``float``
        Maximum velocity dispersion (km/s) for computing upper bound.
    sigma_rvs : ``callable``
        Function to sample velocity dispersion: sigma_rvs(n, zl) -> array.
    q_rvs : ``callable``
        Function to sample axis ratio: q_rvs(n, sigma) -> array.
    phi_rvs : ``callable``
        Function to sample orientation angle: phi_rvs(n) -> array.
    gamma_rvs : ``callable``
        Function to sample power-law index: gamma_rvs(n) -> array.
    shear_rvs : ``callable``
        Function to sample external shear: shear_rvs(n) -> (gamma1, gamma2).
    cross_section : ``callable``
        Function to compute lensing cross section.
    safety_factor : ``float``
        Multiplicative safety factor for the upper bound. \n
        default: 1.2

    Returns
    -------
    sigma_array : ``numpy.ndarray``
        Sampled velocity dispersions (km/s).
    q_array : ``numpy.ndarray``
        Sampled axis ratios.
    phi_array : ``numpy.ndarray``
        Sampled orientation angles (rad).
    gamma_array : ``numpy.ndarray``
        Sampled power-law indices.
    gamma1_array : ``numpy.ndarray``
        Sampled external shear component 1.
    gamma2_array : ``numpy.ndarray``
        Sampled external shear component 2.
    """
    n_samples = zl.size

    sigma_array = np.zeros(n_samples)
    q_array = np.zeros(n_samples)
    phi_array = np.zeros(n_samples)
    gamma_array = np.zeros(n_samples)
    gamma1_array = np.zeros(n_samples)
    gamma2_array = np.zeros(n_samples)
    idx_remaining = np.arange(n_samples)

    # Compute maximum cross section for rejection bound
    cs_max = (
        cross_section(
            zs=zs,
            zl=zl,
            sigma=sigma_max * np.ones(n_samples),
            q=0.9 * np.ones(n_samples),
            phi=np.zeros(n_samples),
            gamma=2.645 * np.ones(n_samples),
            gamma1=np.zeros(n_samples),
            gamma2=np.zeros(n_samples),
        )
        * safety_factor
    )

    while len(idx_remaining) > 0:
        n_remaining = len(idx_remaining)
        sigma_samples = sigma_rvs(n_remaining, zl[idx_remaining])
        q_samples = q_rvs(n_remaining, sigma_samples)
        phi_samples = phi_rvs(n_remaining)
        gamma_samples = gamma_rvs(n_remaining)
        gamma1_samples, gamma2_samples = shear_rvs(n_remaining)

        cs = cross_section(
            zs=zs[idx_remaining],
            zl=zl[idx_remaining],
            sigma=sigma_samples,
            q=q_samples,
            phi=phi_samples,
            gamma=gamma_samples,
            gamma1=gamma1_samples,
            gamma2=gamma2_samples,
        )

        accept = np.random.random(n_remaining) < (cs / cs_max[idx_remaining])

        accepted_indices = idx_remaining[accept]
        sigma_array[accepted_indices] = sigma_samples[accept]
        q_array[accepted_indices] = q_samples[accept]
        phi_array[accepted_indices] = phi_samples[accept]
        gamma_array[accepted_indices] = gamma_samples[accept]
        gamma1_array[accepted_indices] = gamma1_samples[accept]
        gamma2_array[accepted_indices] = gamma2_samples[accept]

        idx_remaining = idx_remaining[~accept]

    return sigma_array, q_array, phi_array, gamma_array, gamma1_array, gamma2_array


def create_rejection_sampler(
    sigma_max,
    sigma_rvs,
    q_rvs,
    phi_rvs,
    gamma_rvs,
    shear_rvs,
    cross_section,
    safety_factor=1.2,
    use_njit_sampler=True,
):
    """
    Create a rejection sampler for cross-section weighted lens parameters.

    Returns a callable that samples lens parameters using rejection sampling,
    weighting by the gravitational lensing cross section. Optionally uses
    Numba JIT compilation for improved performance.

    Parameters
    ----------
    sigma_max : ``float``
        Maximum velocity dispersion (km/s) for computing upper bound.
    sigma_rvs : ``callable``
        Function to sample velocity dispersion: sigma_rvs(n, zl) -> array.
    q_rvs : ``callable``
        Function to sample axis ratio: q_rvs(n, sigma) -> array.
    phi_rvs : ``callable``
        Function to sample orientation angle: phi_rvs(n) -> array.
    gamma_rvs : ``callable``
        Function to sample power-law index: gamma_rvs(n) -> array.
    shear_rvs : ``callable``
        Function to sample external shear: shear_rvs(n) -> (gamma1, gamma2).
    cross_section : ``callable``
        Function to compute lensing cross section.
    safety_factor : ``float``
        Multiplicative safety factor for the upper bound. \n
        default: 1.2
    use_njit_sampler : ``bool``
        If True, uses Numba JIT compilation for faster execution. \n
        default: True

    Returns
    -------
    rejection_sampler_wrapper : ``callable``
        Function with signature (zs, zl) -> (sigma, q, phi, gamma, gamma1, gamma2).

    Examples
    --------
    >>> import numpy as np
    >>> from numba import njit
    >>> @njit
    ... def sigma_rvs(n, zl):
    ...     return 100 + 200 * np.random.random(n)
    >>> @njit
    ... def q_rvs(n, sigma):
    ...     return 0.5 + 0.5 * np.random.random(n)
    >>> @njit
    ... def phi_rvs(n):
    ...     return np.pi * np.random.random(n)
    >>> @njit
    ... def gamma_rvs(n):
    ...     return 2.0 + 0.2 * np.random.randn(n)
    >>> @njit
    ... def shear_rvs(n):
    ...     return 0.05 * np.random.randn(n), 0.05 * np.random.randn(n)
    >>> @njit
    ... def cross_section(zs, zl, sigma, q, phi, gamma, gamma1, gamma2):
    ...     return sigma**4
    >>> sampler = create_rejection_sampler(
    ...     sigma_max=400.0,
    ...     sigma_rvs=sigma_rvs,
    ...     q_rvs=q_rvs,
    ...     phi_rvs=phi_rvs,
    ...     gamma_rvs=gamma_rvs,
    ...     shear_rvs=shear_rvs,
    ...     cross_section=cross_section,
    ... )
    >>> zs = np.array([1.0, 1.5, 2.0])
    >>> zl = np.array([0.3, 0.5, 0.7])
    >>> sigma, q, phi, gamma, gamma1, gamma2 = sampler(zs, zl)
    """
    if use_njit_sampler:
        _base_sampler = njit(rejection_sampler)
        print(
            "Faster, njitted and rejection sampling based lens parameter sampler will be used."
        )

        @njit
        def rejection_sampler_wrapper(zs, zl):
            return _base_sampler(
                zs=zs,
                zl=zl,
                sigma_max=sigma_max,
                sigma_rvs=sigma_rvs,
                q_rvs=q_rvs,
                phi_rvs=phi_rvs,
                gamma_rvs=gamma_rvs,
                shear_rvs=shear_rvs,
                cross_section=cross_section,
                safety_factor=safety_factor,
            )

    else:
        print(
            "Slower, non-njit and rejection sampling based lens parameter sampler will be used."
        )

        def rejection_sampler_wrapper(zs, zl):
            return rejection_sampler(
                zs=zs,
                zl=zl,
                sigma_max=sigma_max,
                sigma_rvs=sigma_rvs,
                q_rvs=q_rvs,
                phi_rvs=phi_rvs,
                gamma_rvs=gamma_rvs,
                shear_rvs=shear_rvs,
                cross_section=cross_section,
                safety_factor=safety_factor,
            )

    return rejection_sampler_wrapper


# --------------------------------------------------
# Importance sampling of strongly lensed parameters
# --------------------------------------------------
@njit
def _sigma_proposal_uniform(n, sigma_min, sigma_max):
    """
    Draw uniform samples for velocity dispersion proposal.

    Parameters
    ----------
    n : ``int``
        Number of samples to draw.
    sigma_min : ``float``
        Minimum velocity dispersion (km/s).
    sigma_max : ``float``
        Maximum velocity dispersion (km/s).

    Returns
    -------
    sigma : ``numpy.ndarray``
        Uniform samples in [sigma_min, sigma_max].
    """
    return sigma_min + (sigma_max - sigma_min) * np.random.random(n)


@njit
def _weighted_choice_1d(weights):
    """
    Draw an index with probability proportional to weights.

    Numba-safe replacement for np.random.choice(n, p=weights).

    Parameters
    ----------
    weights : ``numpy.ndarray``
        Non-negative weights (need not be normalized).

    Returns
    -------
    idx : ``int``
        Randomly selected index.
    """
    total = 0.0
    for i in range(weights.size):
        total += weights[i]

    if not (total > 0.0):
        return np.random.randint(weights.size)

    u = np.random.random() * total
    c = 0.0
    for i in range(weights.size):
        c += weights[i]
        if u <= c:
            return i
    return weights.size - 1


def importance_sampler(
    zs,
    zl,
    sigma_min,
    sigma_max,
    q_rvs,
    phi_rvs,
    gamma_rvs,
    shear_rvs,
    number_density,
    cross_section,
    n_prop,
):
    """
    Core importance sampling algorithm for lens parameters.

    This function samples lens galaxy parameters weighted by their lensing
    cross sections using importance sampling with a uniform proposal distribution
    for velocity dispersion.

    Algorithm
    ---------
    For each lens-source pair (zl_i, zs_i):

    1. **Draw proposal samples** (n_prop samples per lens):
       - sigma_k ~ Uniform(sigma_min, sigma_max)  with proposal density q(sigma) = 1/(sigma_max - sigma_min)
       - q_k ~ p(q|sigma_k)  (axis ratio conditioned on velocity dispersion)
       - φ_k ~ p(φ)  (orientation angle prior)
       - gamma_k ~ p(gamma)  (power-law index prior)
       - (gamma1_k, gamma2_k) ~ p(gamma1, gamma2)  (external shear prior)

    2. **Compute lensing cross sections**:
       - cs_k = CrossSection(zs_i, zl_i, sigma_k, q_k, φ_k, gamma_k, gamma1_k, gamma2_k)
       - Normalize: cs_k ← cs_k / Σ_k cs_k

    3. **Compute importance weights**:
       - The target distribution is: p(θ|zl) ∝ p(sigma|zl) x CrossSection(θ)
       - The proposal distribution is: q(θ) ∝ Uniform(sigma) x p(q|sigma) x p(φ) x p(gamma) x p(shear)
       - Importance weight: w_k = cs_k x [p(sigma_k|zl) / q(sigma)]
       - Normalize: w_k ← w_k / Σ_k w_k

    4. **Resample**:
       - Draw one sample index from {1, ..., n_prop} with probabilities {w_1, ..., w_n_prop}
       - Return the corresponding parameter values as the posterior sample

    The algorithm produces samples from the posterior distribution of lens
    parameters weighted by their contribution to the strong lensing cross section.

    Parameters
    ----------
    zs : ``numpy.ndarray``
        Source redshifts.
    zl : ``numpy.ndarray``
        Lens redshifts.
    sigma_min : ``float``
        Minimum velocity dispersion (km/s) for uniform proposal.
    sigma_max : ``float``
        Maximum velocity dispersion (km/s) for uniform proposal.
    q_rvs : ``callable``
        Function to sample axis ratio: q_rvs(n, sigma) -> array.
    phi_rvs : ``callable``
        Function to sample orientation angle: phi_rvs(n) -> array.
    gamma_rvs : ``callable``
        Function to sample power-law index: gamma_rvs(n) -> array.
    shear_rvs : ``callable``
        Function to sample external shear: shear_rvs(n) -> (gamma1, gamma2).
    number_density : ``callable``
        Number density or velocity dispersion function: number_density(sigma, zl) -> array.
    cross_section : ``callable``
        Function to compute lensing cross section.
    n_prop : ``int``
        Number of proposal samples per lens.

    Returns
    -------
    sigma_post : ``numpy.ndarray``
        Sampled velocity dispersions (km/s).
    q_post : ``numpy.ndarray``
        Sampled axis ratios.
    phi_post : ``numpy.ndarray``
        Sampled orientation angles (rad).
    gamma_post : ``numpy.ndarray``
        Sampled power-law indices.
    gamma1_post : ``numpy.ndarray``
        Sampled external shear component 1.
    gamma2_post : ``numpy.ndarray``
        Sampled external shear component 2.
    """
    n_samples = zl.size

    sigma_post = np.zeros(n_samples)
    q_post = np.zeros(n_samples)
    phi_post = np.zeros(n_samples)
    gamma_post = np.zeros(n_samples)
    gamma1_post = np.zeros(n_samples)
    gamma2_post = np.zeros(n_samples)

    p0 = 1.0 / (sigma_max - sigma_min)

    for i in prange(n_samples):  # for each (zl, zs) pair
        # Draw proposals from uniform distribution
        sigma_prop = _sigma_proposal_uniform(n_prop, sigma_min, sigma_max)

        # Draw other parameters from their priors
        q_prop = q_rvs(n_prop, sigma_prop)
        phi_prop = phi_rvs(n_prop)
        gamma_prop = gamma_rvs(n_prop)
        gamma1_prop, gamma2_prop = shear_rvs(n_prop)

        # Compute cross sections
        zs_arr = zs[i] * np.ones(n_prop)
        zl_arr = zl[i] * np.ones(n_prop)

        cs = cross_section(
            zs_arr,
            zl_arr,
            sigma_prop,
            q_prop,
            phi_prop,
            gamma_prop,
            gamma1_prop,
            gamma2_prop,
        )

        # Compute importance weights
        sigma_function = number_density(sigma_prop, zl_arr)

        w = cs * (sigma_function / p0)
        w = np.where(w > 0.0, w, 0.0)
        w_sum = np.sum(w)

        # Normalize weights
        if w_sum > 0.0:
            w = w / w_sum
        else:
            w = np.ones(n_prop) / n_prop

        # Draw posterior sample via weighted choice
        idx = _weighted_choice_1d(w)

        sigma_post[i] = sigma_prop[idx]
        q_post[i] = q_prop[idx]
        phi_post[i] = phi_prop[idx]
        gamma_post[i] = gamma_prop[idx]
        gamma1_post[i] = gamma1_prop[idx]
        gamma2_post[i] = gamma2_prop[idx]

    return sigma_post, q_post, phi_post, gamma_post, gamma1_post, gamma2_post


# --------------------------------------
# Importance sampling (multiprocessing)
# --------------------------------------
def _importance_sampler_worker(params):
    """
    Worker function for multiprocessing importance sampling.

    Parameters
    ----------
    params : ``tuple``
        Packed parameters: (zs_i, zl_i, worker_idx).

    Returns
    -------
    worker_idx : ``int``
        Worker index for result ordering.
    result : ``tuple``
        (sigma, q, phi, gamma, gamma1, gamma2) for this sample.
    """
    zs_i, zl_i, worker_idx = params

    # Load shared data from pickle file
    shared_data = load_pickle("importance_sampler_shared.pkl")
    sigma_min = shared_data["sigma_min"]
    sigma_max = shared_data["sigma_max"]
    q_rvs = shared_data["q_rvs"]
    phi_rvs = shared_data["phi_rvs"]
    gamma_rvs = shared_data["gamma_rvs"]
    shear_rvs = shared_data["shear_rvs"]
    number_density = shared_data["number_density"]
    cross_section = shared_data["cross_section"]
    n_prop = shared_data["n_prop"]

    p0 = 1.0 / (sigma_max - sigma_min)

    # Draw proposals
    sigma_prop = np.random.uniform(sigma_min, sigma_max, n_prop)

    # Draw other parameters from their priors
    q_prop = q_rvs(n_prop, sigma_prop)
    phi_prop = phi_rvs(n_prop)
    gamma_prop = gamma_rvs(n_prop)
    gamma1_prop, gamma2_prop = shear_rvs(n_prop)

    # Compute cross sections
    zs_arr = zs_i * np.ones(n_prop)
    zl_arr = zl_i * np.ones(n_prop)

    cs = cross_section(
        zs_arr,
        zl_arr,
        sigma_prop,
        q_prop,
        phi_prop,
        gamma_prop,
        gamma1_prop,
        gamma2_prop,
    )/ (4.0*np.pi)

    # Compute importance weights
    sigma_function = number_density(sigma_prop, zl_arr)

    w = cs * (sigma_function / p0)
    w = np.where(w > 0.0, w, 0.0)
    w_sum = np.sum(w)

    # Normalize weights
    if w_sum > 0.0:
        w = w / w_sum
    else:
        w = np.ones(n_prop) / n_prop

    # Draw posterior sample via weighted choice
    idx = np.random.choice(n_prop, p=w)

    result = (
        sigma_prop[idx],
        q_prop[idx],
        phi_prop[idx],
        gamma_prop[idx],
        gamma1_prop[idx],
        gamma2_prop[idx],
    )

    return worker_idx, result


def importance_sampler_mp(
    zs,
    zl,
    sigma_min,
    sigma_max,
    q_rvs,
    phi_rvs,
    gamma_rvs,
    shear_rvs,
    number_density,
    cross_section,
    n_prop,
    npool=4,
):
    """
    Multiprocessing version of importance sampling for lens parameters.

    Parameters
    ----------
    zs : ``numpy.ndarray``
        Source redshifts.
    zl : ``numpy.ndarray``
        Lens redshifts.
    sigma_min : ``float``
        Minimum velocity dispersion (km/s) for uniform proposal.
    sigma_max : ``float``
        Maximum velocity dispersion (km/s) for uniform proposal.
    q_rvs : ``callable``
        Function to sample axis ratio: q_rvs(n, sigma) -> array.
    phi_rvs : ``callable``
        Function to sample orientation angle: phi_rvs(n) -> array.
    gamma_rvs : ``callable``
        Function to sample power-law index: gamma_rvs(n) -> array.
    shear_rvs : ``callable``
        Function to sample external shear: shear_rvs(n) -> (gamma1, gamma2).
    number_density : ``callable``
        Number density or velocity dispersion function: number_density(sigma, zl) -> array.
    cross_section : ``callable``
        Function to compute lensing cross section.
    n_prop : ``int``
        Number of proposal samples per lens.
    npool : ``int``
        Number of parallel processes to use. \n
        default: 4

    Returns
    -------
    sigma_post : ``numpy.ndarray``
        Sampled velocity dispersions (km/s).
    q_post : ``numpy.ndarray``
        Sampled axis ratios.
    phi_post : ``numpy.ndarray``
        Sampled orientation angles (rad).
    gamma_post : ``numpy.ndarray``
        Sampled power-law indices.
    gamma1_post : ``numpy.ndarray``
        Sampled external shear component 1.
    gamma2_post : ``numpy.ndarray``
        Sampled external shear component 2.
    """
    n_samples = zl.size

    # Save shared data to pickle file for workers
    shared_data = {
        "sigma_min": sigma_min,
        "sigma_max": sigma_max,
        "q_rvs": q_rvs,
        "phi_rvs": phi_rvs,
        "gamma_rvs": gamma_rvs,
        "shear_rvs": shear_rvs,
        "number_density": number_density,
        "cross_section": cross_section,
        "n_prop": n_prop,
    }
    save_pickle("importance_sampler_shared.pkl", shared_data)

    # Prepare input parameters for workers
    input_params = [(zs[i], zl[i], i) for i in range(n_samples)]

    # Initialize output arrays
    sigma_post = np.zeros(n_samples)
    q_post = np.zeros(n_samples)
    phi_post = np.zeros(n_samples)
    gamma_post = np.zeros(n_samples)
    gamma1_post = np.zeros(n_samples)
    gamma2_post = np.zeros(n_samples)

    # Run multiprocessing with progress bar
    with Pool(processes=npool) as pool:
        for worker_idx, result in tqdm(
            pool.imap_unordered(_importance_sampler_worker, input_params),
            total=n_samples,
            ncols=100,
            desc="Importance sampling",
        ):
            sigma_post[worker_idx] = result[0]
            q_post[worker_idx] = result[1]
            phi_post[worker_idx] = result[2]
            gamma_post[worker_idx] = result[3]
            gamma1_post[worker_idx] = result[4]
            gamma2_post[worker_idx] = result[5]

    # Cleanup pickle file
    os.remove("importance_sampler_shared.pkl")

    return sigma_post, q_post, phi_post, gamma_post, gamma1_post, gamma2_post


def create_importance_sampler(
    sigma_min,
    sigma_max,
    q_rvs,
    phi_rvs,
    gamma_rvs,
    shear_rvs,
    number_density,
    cross_section,
    n_prop,
    use_njit_sampler=True,
    npool=4,
):
    """
    Create an importance sampler for cross-section weighted lens parameters.

    Returns a callable that samples lens parameters using importance sampling
    with uniform proposal distribution, optionally JIT-compiled for improved
    performance.

    Parameters
    ----------
    sigma_min : ``float``
        Minimum velocity dispersion (km/s) for uniform proposal.
    sigma_max : ``float``
        Maximum velocity dispersion (km/s) for uniform proposal.
    q_rvs : ``callable``
        Function to sample axis ratio: q_rvs(n, sigma) -> array.
    phi_rvs : ``callable``
        Function to sample orientation angle: phi_rvs(n) -> array.
    gamma_rvs : ``callable``
        Function to sample power-law index: gamma_rvs(n) -> array.
    shear_rvs : ``callable``
        Function to sample external shear: shear_rvs(n) -> (gamma1, gamma2).
    number_density : ``callable``
        Number density or velocity dispersion function: number_density(sigma, zl) -> array.
    cross_section : ``callable``
        Function to compute lensing cross section.
    n_prop : ``int``
        Number of proposal samples per lens.
    use_njit_sampler : ``bool``
        If True, uses Numba JIT compilation for faster execution. \n
        default: True
    npool : ``int``
        Number of parallel processes (only used when use_njit_sampler=False). \n
        default: 4

    Returns
    -------
    importance_sampler_wrapper : ``callable``
        Function with signature (zs, zl) -> (sigma, q, phi, gamma, gamma1, gamma2).

    Examples
    --------
    >>> import numpy as np
    >>> from numba import njit
    >>> @njit
    ... def q_rvs(n, sigma):
    ...     return 0.5 + 0.5 * np.random.random(n)
    >>> @njit
    ... def phi_rvs(n):
    ...     return np.pi * np.random.random(n)
    >>> @njit
    ... def gamma_rvs(n):
    ...     return 2.0 + 0.2 * np.random.randn(n)
    >>> @njit
    ... def shear_rvs(n):
    ...     return 0.05 * np.random.randn(n), 0.05 * np.random.randn(n)
    >>> @njit
    ... def number_density(sigma, zl):
    ...     return np.ones_like(sigma)
    >>> @njit
    ... def cross_section(zs, zl, sigma, q, phi, gamma, gamma1, gamma2):
    ...     return sigma**4
    >>> sampler = create_importance_sampler(
    ...     sigma_min=100.0,
    ...     sigma_max=400.0,
    ...     q_rvs=q_rvs,
    ...     phi_rvs=phi_rvs,
    ...     gamma_rvs=gamma_rvs,
    ...     shear_rvs=shear_rvs,
    ...     number_density=number_density,
    ...     cross_section=cross_section,
    ...     n_prop=100,
    ... )
    >>> zs = np.array([1.0, 1.5, 2.0])
    >>> zl = np.array([0.3, 0.5, 0.7])
    >>> sigma, q, phi, gamma, gamma1, gamma2 = sampler(zs, zl)
    """
    if use_njit_sampler:
        print(
            "Faster, njitted and importance sampling based lens parameter sampler will be used."
        )
        _base_sampler = njit(parallel=True)(importance_sampler)

        @njit(parallel=True)
        def importance_sampler_wrapper(zs, zl):
            return _base_sampler(
                zs=zs,
                zl=zl,
                sigma_min=sigma_min,
                sigma_max=sigma_max,
                q_rvs=q_rvs,
                phi_rvs=phi_rvs,
                gamma_rvs=gamma_rvs,
                shear_rvs=shear_rvs,
                number_density=number_density,
                cross_section=cross_section,
                n_prop=n_prop,
            )

    else:
        print(
            "Slower, non-njit and importance sampling based lens parameter sampler will be used."
        )

        def importance_sampler_wrapper(zs, zl):
            return importance_sampler_mp(
                zs=zs,
                zl=zl,
                sigma_min=sigma_min,
                sigma_max=sigma_max,
                q_rvs=q_rvs,
                phi_rvs=phi_rvs,
                gamma_rvs=gamma_rvs,
                shear_rvs=shear_rvs,
                number_density=number_density,
                cross_section=cross_section,
                n_prop=n_prop,
                npool=npool,
            )

    return importance_sampler_wrapper


def _njit_checks(
    sigma_rvs_,
    q_rvs_,
    phi_rvs,
    gamma_rvs,
    shear_rvs,
    sigma_pdf_,
    number_density_,
    cross_section_,
):
    """
    Check and wrap sampler functions for JIT compatibility.

    Parameters
    ----------
    sigma_rvs_ : ``callable``
        Function to sample velocity dispersion.
    q_rvs_ : ``callable``
        Function to sample axis ratio.
    phi_rvs : ``callable``
        Function to sample orientation angle.
    gamma_rvs : ``callable``
        Function to sample power-law index.
    shear_rvs : ``callable``
        Function to sample external shear.
    sigma_pdf_ : ``callable``
        PDF of velocity dispersion.
    number_density_ : ``callable``
        Number density function of lens galaxies.
    cross_section_ : ``callable``
        Function to compute lensing cross section.

    Returns
    -------
    use_njit_sampler : ``bool``
        True if all functions are JIT compiled, False otherwise.
    dict_ : ``dict``
        Dictionary containing wrapped/compiled versions of all functions.
    """
    # Wrap cross_section function based on argument count
    if cross_section_.__code__.co_argcount == 4:
        cross_section_function = njit(
            lambda zs, zl, sigma, q, phi, gamma, gamma1, gamma2: cross_section_(
                zs, zl, sigma, q
            )
        )
    elif cross_section_.__code__.co_argcount == 3:
        cross_section_function = njit(
            lambda zs, zl, sigma, q, phi, gamma, gamma1, gamma2: cross_section_(
                zs, zl, sigma
            )
        )
    else:
        cross_section_function = cross_section_

    # Wrap samplers and PDFs based on argument count
    if sigma_rvs_.__code__.co_argcount == 1:
        if is_njitted(sigma_rvs_):
            sigma_rvs = njit(lambda size, zl: sigma_rvs_(size))
            sigma_pdf = njit(lambda sigma, zl: sigma_pdf_(sigma))
            number_density = njit(lambda sigma, zl: number_density_(sigma))
        else:
            sigma_rvs = lambda size, zl: sigma_rvs_(size)
            sigma_pdf = lambda sigma, zl: sigma_pdf_(sigma)
            number_density = lambda sigma, zl: number_density_(sigma)
    else:
        sigma_rvs = sigma_rvs_
        sigma_pdf = sigma_pdf_
        number_density = number_density_

    if q_rvs_.__code__.co_argcount == 1:
        if is_njitted(q_rvs_):
            q_rvs = njit(lambda size, sigma: q_rvs_(size))
        else:
            q_rvs = lambda size, sigma: q_rvs_(size)
    else:
        q_rvs = q_rvs_

    # Build dictionary of wrapped functions
    dict_ = {
        "sigma_rvs": sigma_rvs,
        "q_rvs": q_rvs,
        "phi_rvs": phi_rvs,
        "gamma_rvs": gamma_rvs,
        "shear_rvs": shear_rvs,
        "sigma_pdf": sigma_pdf,
        "number_density": number_density,
        "cross_section_function": cross_section_function,
    }

    # Check if all functions are JIT compiled
    use_njit_sampler = True
    for key, value in dict_.items():
        if not is_njitted(value):
            print(f"Warning: {key} is not njitted.")
            use_njit_sampler = False

    return use_njit_sampler, dict_
