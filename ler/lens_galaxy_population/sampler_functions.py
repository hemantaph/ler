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
from numba import njit, prange, set_num_threads
from scipy.interpolate import CubicSpline
from astropy.cosmology import LambdaCDM
from ..utils import is_njitted
from multiprocessing import Pool
from tqdm import tqdm
from ..utils import (
    inverse_transform_sampler,
    generate_mixed_grid,
)

# ---------------------------------
# Lens redshift sampler functions
# ---------------------------------
# NOTE: cache=True intentionally NOT set. ``cosmo.comoving_distance`` is
# an astropy method, not a Numba-compatible call, so this function only
# runs in object-mode fallback and is not cacheable.
@njit(fastmath=True)
def lens_redshift_strongly_lensed_sis_analytical_pdf(
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
    >>> pdf = lens_redshift_strongly_lensed_sis_analytical_pdf(zl=0.5, zs=1.0, cosmo=cosmo)
    >>> print(f"PDF at zl=0.5: {pdf:.4f}")
    PDF at zl=0.5: 1.8750
    """
    Dc_zl = cosmo.comoving_distance(zl).value
    Dc_zs = cosmo.comoving_distance(zs).value
    x = Dc_zl / Dc_zs
    return 30 * x**2 * (1 - x) ** 2


def lens_redshift_strongly_lensed_sis_analytical_rvs(
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
    >>> zl_samples = lens_redshift_strongly_lensed_sis_analytical_rvs(
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
    zs_array = generate_mixed_grid(z_min, z_max, 500)
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
@njit(cache=True, fastmath=True)
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


@njit(cache=True, fastmath=True)
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


@njit(cache=True, fastmath=True)
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


@njit(cache=True, fastmath=True)
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


@njit(cache=True, fastmath=True)
def velocity_dispersion_ewoud_denisty_function(
    sigma, z, alpha=0.94, beta=1.85, phistar=2.099e-2, sigmastar=113.78
):
    """
    Calculate the lens galaxy velocity dispersion function at redshift z (Oguri et al. (2018b) + Wempe et al. (2022)).

    Notes
    -----
    The function name contains the misspelling ``denisty`` for backward
    compatibility with existing serialized configuration dictionaries.

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


@njit(cache=True, fastmath=True)
def velocity_dispersion_bernardi_denisty_function(
    sigma, alpha, beta, phistar, sigmastar
):
    """
    Calculate the local universe velocity dispersion function.

    This implements the velocity dispersion function from Bernardi et al. (2010).

    .. math::

        \\phi(\\sigma) = \\phi_* \\left(\\frac{\\sigma}{\\sigma_*}\\right)^\\alpha
        \\exp\\left[-\\left(\\frac{\\sigma}{\\sigma_*}\\right)^\\beta\\right]
        \\frac{\\beta}{\\Gamma(\\alpha/\\beta)\\sigma}.

    Notes
    -----
    The function name contains the misspelling ``denisty`` for backward
    compatibility with existing serialized configuration dictionaries.

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


def gengamma_function(
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
    >>> density = gengamma_function(sigma)
    >>> print(f"Density at sigma=200 km/s: {density[1]:.6f}")
    """
    from scipy.stats import gengamma

    density = phistar * gengamma.pdf(
        sigma / sigmastar,
        a=alpha / beta,
        c=beta,
    )

    return density


def gengamma_pdf(
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
    >>> pdf = gengamma_pdf(
    ...     sigma=200.0,
    ...     sigma_min=100.0,
    ...     sigma_max=400.0,
    ... )
    >>> print(f"PDF at sigma=200 km/s: {pdf:.6f}")
    """
    # Compute normalization constant
    sigma_array = np.linspace(sigma_min, sigma_max, 500)
    density = gengamma_function(
        sigma=sigma_array,
        alpha=alpha,
        beta=beta,
        sigmastar=sigmastar,
    )
    integrate = np.trapezoid if hasattr(np, "trapezoid") else np.trapz
    norm_const = integrate(density, sigma_array)

    pdf = (
        gengamma_function(
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


def gengamma_rvs(
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
    >>> sigma_samples = gengamma_rvs(
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
@njit(cache=True, fastmath=True)
def rayleigh_rvs(size, sigma, q_min=0.2, q_max=1.0):
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
    >>> q_samples = rayleigh_rvs(size=1000, sigma=sigma, q_min=0.2, q_max=1.0)
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


@njit(parallel=True, fastmath=True)
def rayleigh_pdf(q, sigma, q_min=0.2, q_max=1.0):
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
    >>> pdf = rayleigh_pdf(q, sigma)
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


@njit(cache=True, fastmath=True)
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


@njit(cache=True, fastmath=True)
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


# NOTE: cache=True intentionally NOT set. This function calls
# ``scipy.interpolate.CubicSpline``, which is not Numba-compatible and
# falls back to object mode, so the compiled version cannot be cached.
@njit(fastmath=True)
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


# @njit
# def bounded_normal_sample(size, mean, std, low, high):
#     """
#     Sample from truncated normal distribution via rejection.

#     Generates samples from a normal distribution with specified mean and
#     standard deviation, rejecting samples outside the specified bounds.

#     Parameters
#     ----------
#     size : ``int``
#         Number of samples to draw.
#     mean : ``float``
#         Mean of the normal distribution.
#     std : ``float``
#         Standard deviation of the normal distribution.
#     low : ``float``
#         Lower bound for samples.
#     high : ``float``
#         Upper bound for samples.

#     Returns
#     -------
#     samples : ``numpy.ndarray``
#         Bounded normal samples with shape (size,).

#     Examples
#     --------
#     >>> samples = bounded_normal_sample(size=1000, mean=2.0, std=0.2, low=1.5, high=2.5)
#     >>> print(f"Mean: {samples.mean():.2f}, Std: {samples.std():.2f}")
#     >>> print(f"Range: [{samples.min():.2f}, {samples.max():.2f}]")
#     """
#     samples = np.empty(size)
#     for i in range(size):
#         while True:
#             sample = np.random.normal(mean, std)
#             if low <= sample <= high:
#                 break
#         samples[i] = sample

#     return samples


# -------------------------------------------------
# Cross-section sampler helpers
# -------------------------------------------------
def _njit_checks(
    zs_pdf=None,
    zs_rvs=None,
    zl_pdf=None,
    zl_rvs=None,
    sigma_pdf=None,
    sigma_rvs=None,
    q_pdf=None,
    q_rvs=None,
    phi_pdf=None,
    phi_rvs=None,
    gamma_pdf=None,
    gamma_rvs=None,
    shear1_pdf=None,
    shear1_rvs=None,
    shear2_pdf=None,
    shear2_rvs=None,
    cross_section_function=None,
    number_density_function=None,
    dVcdz_function=None,
    create_njit_sampler=False,
):
    """Normalize callable signatures used by the lens redshift integrators."""

    def _arg_count(function):
        if function is None:
            return 0
        if hasattr(function, "__code__"):
            return function.__code__.co_argcount
        if hasattr(function, "py_func"):
            return function.py_func.__code__.co_argcount
        return 0

    def _one_arg_pdf(pdf):
        return pdf

    def _two_arg_pdf(pdf):
        if pdf is None or _arg_count(pdf) != 1:
            return pdf

        def wrapped(x, y=None):
            return pdf(x)

        if is_njitted(pdf) and create_njit_sampler:
            wrapped = njit(wrapped)
        return wrapped

    def _two_arg_rvs(rvs):
        if rvs is None or _arg_count(rvs) != 1:
            return rvs

        def wrapped(size, y=None):
            return rvs(size)

        if is_njitted(rvs) and create_njit_sampler:
            wrapped = njit(wrapped)
        return wrapped

    def _cross_section(cross_section):
        if cross_section is None:
            return None
        arg_count = _arg_count(cross_section)
        if arg_count == 3:

            def wrapped(zs, zl, sigma, q, phi, gamma, gamma1, gamma2):
                return cross_section(zs, zl, sigma)

        elif arg_count == 4:

            def wrapped(zs, zl, sigma, q, phi, gamma, gamma1, gamma2):
                return cross_section(zs, zl, sigma, q)

        else:
            return cross_section

        if is_njitted(cross_section) and create_njit_sampler:
            wrapped = njit(wrapped)
        return wrapped

    dict_ = dict(
        zs_pdf=_one_arg_pdf(zs_pdf),
        zs_rvs=zs_rvs,
        zl_pdf=_two_arg_pdf(zl_pdf),
        zl_rvs=_two_arg_rvs(zl_rvs),
        sigma_pdf=_two_arg_pdf(sigma_pdf),
        sigma_rvs=_two_arg_rvs(sigma_rvs),
        q_pdf=_two_arg_pdf(q_pdf),
        q_rvs=_two_arg_rvs(q_rvs),
        phi_pdf=_one_arg_pdf(phi_pdf),
        phi_rvs=phi_rvs,
        gamma_pdf=gamma_pdf,
        gamma_rvs=gamma_rvs,
        shear1_pdf=shear1_pdf,
        shear1_rvs=shear1_rvs,
        shear2_pdf=shear2_pdf,
        shear2_rvs=shear2_rvs,
        cross_section_function=_cross_section(cross_section_function),
        number_density_function=_two_arg_pdf(number_density_function),
        dVcdz_function=_one_arg_pdf(dVcdz_function),
    )

    use_njit_sampler = False
    if create_njit_sampler:
        use_njit_sampler = True
        for key, value in dict_.items():
            if value is not None and not is_njitted(value):
                print(f"Warning: {key} is not njitted.")
                use_njit_sampler = False

    return use_njit_sampler, dict_


# -------------------------------------------------
# Simple uniform proposal samplers
# -------------------------------------------------
# def _wrap_cross_section(cross_section):
#     """Return an 8-argument cross-section function."""
#     if hasattr(cross_section, "__code__"):
#         arg_count = cross_section.__code__.co_argcount
#     elif hasattr(cross_section, "py_func"):
#         arg_count = cross_section.py_func.__code__.co_argcount
#     else:
#         arg_count = 8

#     if arg_count == 3:
#         return lambda zs, zl, sigma, q, phi, gamma, gamma1, gamma2: cross_section(
#             zs, zl, sigma
#         )
#     if arg_count == 4:
#         return lambda zs, zl, sigma, q, phi, gamma, gamma1, gamma2: cross_section(
#             zs, zl, sigma, q
#         )
#     return cross_section

# -------------------------------------------------
# Main rejection sampler function
# -------------------------------------------------
@njit(cache=True, fastmath=True)
def _uniform_pdf(x_min, x_max):
    if x_max > x_min:
        return 1.0 / (x_max - x_min)
    return 1.0


@njit(cache=True, fastmath=True)
def _chunk_bounds(i, size, npool):
    size_i = size // npool
    start = i * size_i
    end = start + size_i
    if i == npool - 1:
        end += size % npool
    return start, end


@njit(cache=True, fastmath=True)
def _weighted_sample_index(w, n_prop):
    w_sum = np.sum(w)
    if w_sum <= 0.0:
        idx = int(np.random.random() * n_prop)
        if idx >= n_prop:
            idx = n_prop - 1
        return idx

    cdf = 0.0
    r = np.random.random() * w_sum
    for i in range(n_prop):
        cdf += w[i]
        if r <= cdf:
            return i

    return n_prop - 1


# @njit(parallel=True, fastmath=True)
def rejection_sampler_full(
    size,
    zs_pdf,
    number_density,
    q_pdf,
    phi_pdf,
    gamma_pdf,
    shear1_pdf,
    shear2_pdf,
    cross_section,
    dVdz,
    ranges,
    threshold_factor=1e-4,
    n_prop=1000,
    lens_type="epl_shear_galaxy",
    npool=4,
):
    """
    Joint rejection sample over (zs, zl, sigma, nuisance).

    Intrinsic target (before cross section) uses ``p(zs) * dV/dz_l * n(sigma,zl) * ...``;
    ``z_l`` is proposed uniform in ``[z_min, z_s)`` — no separate ``p(zl|zs)``
    factor (see ``importance_sampler_full`` docstring).
    """
    zs_post = np.zeros(size)
    zl_post = np.zeros(size)
    sigma_post = np.zeros(size)
    q_post = np.ones(size)
    phi_post = np.zeros(size)
    gamma_post = np.ones(size) * 2.0
    gamma1_post = np.zeros(size)
    gamma2_post = np.zeros(size)

    z_min, z_max, sigma_min, sigma_max, q_min, q_max, phi_min, phi_max, gamma_min, gamma_max, shear_min, shear_max = ranges

    proposal_pdf = 1.0
    proposal_pdf *= _uniform_pdf(z_min, z_max)
    proposal_pdf *= _uniform_pdf(sigma_min, sigma_max)

    if lens_type in ["epl_shear_galaxy"]:

        proposal_pdf *= _uniform_pdf(q_min, q_max)
        proposal_pdf *= _uniform_pdf(phi_min, phi_max)
        proposal_pdf *= _uniform_pdf(gamma_min, gamma_max)
        proposal_pdf *= _uniform_pdf(shear_min, shear_max)
        proposal_pdf *= _uniform_pdf(shear_min, shear_max)

        for i in prange(npool):

            start, end = _chunk_bounds(i, size, npool)

            size_i = end - start

            zs_collected = np.zeros(size_i)
            zl_collected = np.zeros(size_i)
            sigma_collected = np.zeros(size_i)
            q_collected = np.zeros(size_i)
            phi_collected = np.zeros(size_i)
            gamma_collected = np.zeros(size_i)
            gamma1_collected = np.zeros(size_i)
            gamma2_collected = np.zeros(size_i)

            count = 0
            w_max = 0.0

            while count < size_i:

                zs_prop = z_min + (z_max - z_min) * np.random.random(n_prop)
                zl_prop = z_min + (zs_prop - z_min) * np.random.random(n_prop)

                proposal_pdf_i = proposal_pdf / (zs_prop - z_min)

                sigma_prop = sigma_min + (sigma_max - sigma_min) * np.random.random(n_prop)
                q_prop = q_min + (q_max - q_min) * np.random.random(n_prop)
                phi_prop = phi_min + (phi_max - phi_min) * np.random.random(n_prop)
                gamma_prop = gamma_min + (gamma_max - gamma_min) * np.random.random(n_prop)
                gamma1_prop = shear_min + (shear_max - shear_min) * np.random.random(n_prop)
                gamma2_prop = shear_min + (shear_max - shear_min) * np.random.random(n_prop)

                intrinsic_pdf = (
                    zs_pdf(zs_prop)
                    * dVdz(zl_prop)
                    * number_density(sigma_prop, zl_prop)
                    * q_pdf(q_prop, sigma_prop)
                    * phi_pdf(phi_prop)
                    * gamma_pdf(gamma_prop)
                    * shear1_pdf(gamma1_prop)
                    * shear2_pdf(gamma2_prop)
                )

                valid = (
                    np.isfinite(proposal_pdf_i)
                    & (proposal_pdf_i > 0.0)
                    & np.isfinite(intrinsic_pdf)
                    & (intrinsic_pdf > 0.0)
                )
                intrinsic_by_proposal = np.zeros(n_prop)
                intrinsic_by_proposal[valid] = intrinsic_pdf[valid] / proposal_pdf_i[valid]

                if threshold_factor <= 0.0:
                    cs = cross_section(
                        zs_prop,
                        zl_prop,
                        sigma_prop,
                        q_prop,
                        phi_prop,
                        gamma_prop,
                        gamma1_prop,
                        gamma2_prop,
                    )
                else:
                    w_filter = sigma_prop**4 * intrinsic_by_proposal
                    max_filter = np.max(w_filter)

                    cs = np.zeros(n_prop)
                    if max_filter > 0.0:
                        threshold = max_filter * threshold_factor
                        idx_filter = w_filter > threshold
                        if np.any(idx_filter):
                            cs[idx_filter] = cross_section(
                                zs_prop[idx_filter],
                                zl_prop[idx_filter],
                                sigma_prop[idx_filter],
                                q_prop[idx_filter],
                                phi_prop[idx_filter],
                                gamma_prop[idx_filter],
                                gamma1_prop[idx_filter],
                                gamma2_prop[idx_filter],
                            )

                w = np.zeros(n_prop)
                idx_weight = valid & np.isfinite(cs) & (cs > 0.0)
                w[idx_weight] = cs[idx_weight] * intrinsic_by_proposal[idx_weight] / (4.0 * np.pi)

                batch_max = np.max(w)
                if batch_max <= 0.0:
                    continue

                if batch_max > w_max:
                    if w_max > 0.0:
                        keep_prob = w_max / batch_max
                        new_count = 0

                        for j in range(count):
                            if np.random.random() < keep_prob:
                                zs_collected[new_count] = zs_collected[j]
                                zl_collected[new_count] = zl_collected[j]
                                sigma_collected[new_count] = sigma_collected[j]
                                q_collected[new_count] = q_collected[j]
                                phi_collected[new_count] = phi_collected[j]
                                gamma_collected[new_count] = gamma_collected[j]
                                gamma1_collected[new_count] = gamma1_collected[j]
                                gamma2_collected[new_count] = gamma2_collected[j]
                                new_count += 1

                        count = new_count

                    w_max = batch_max

                for j in range(n_prop):
                    if count >= size_i:
                        break

                    if np.random.random() < (w[j] / w_max):
                        zs_collected[count] = zs_prop[j]
                        zl_collected[count] = zl_prop[j]
                        sigma_collected[count] = sigma_prop[j]
                        q_collected[count] = q_prop[j]
                        phi_collected[count] = phi_prop[j]
                        gamma_collected[count] = gamma_prop[j]
                        gamma1_collected[count] = gamma1_prop[j]
                        gamma2_collected[count] = gamma2_prop[j]
                        count += 1

            zs_post[start:end] = zs_collected
            zl_post[start:end] = zl_collected
            sigma_post[start:end] = sigma_collected
            q_post[start:end] = q_collected
            phi_post[start:end] = phi_collected
            gamma_post[start:end] = gamma_collected
            gamma1_post[start:end] = gamma1_collected
            gamma2_post[start:end] = gamma2_collected

    if lens_type in ["sie_galaxy"]:

        proposal_pdf *= _uniform_pdf(q_min, q_max)
        proposal_pdf *= _uniform_pdf(phi_min, phi_max)

        gamma_prop = np.ones(n_prop) * 2.0
        gamma1_prop = np.zeros(n_prop)
        gamma2_prop = np.zeros(n_prop)

        for i in prange(npool):

            start, end = _chunk_bounds(i, size, npool)

            size_i = end - start

            zs_collected = np.zeros(size_i)
            zl_collected = np.zeros(size_i)
            sigma_collected = np.zeros(size_i)
            q_collected = np.zeros(size_i)
            phi_collected = np.zeros(size_i)

            count = 0
            w_max = 0.0

            while count < size_i:

                zs_prop = z_min + (z_max - z_min) * np.random.random(n_prop)
                zl_prop = z_min + (zs_prop - z_min) * np.random.random(n_prop)

                proposal_pdf_i = proposal_pdf / (zs_prop - z_min)

                sigma_prop = sigma_min + (sigma_max - sigma_min) * np.random.random(n_prop)
                q_prop = q_min + (q_max - q_min) * np.random.random(n_prop)
                phi_prop = phi_min + (phi_max - phi_min) * np.random.random(n_prop)

                intrinsic_pdf = (
                    zs_pdf(zs_prop)
                    * dVdz(zl_prop)
                    * number_density(sigma_prop, zl_prop)
                    * q_pdf(q_prop, sigma_prop)
                    * phi_pdf(phi_prop)
                )

                valid = (
                    np.isfinite(proposal_pdf_i)
                    & (proposal_pdf_i > 0.0)
                    & np.isfinite(intrinsic_pdf)
                    & (intrinsic_pdf > 0.0)
                )
                intrinsic_by_proposal = np.zeros(n_prop)
                intrinsic_by_proposal[valid] = intrinsic_pdf[valid] / proposal_pdf_i[valid]

                if threshold_factor <= 0.0:
                    cs = cross_section(
                        zs_prop,
                        zl_prop,
                        sigma_prop,
                        q_prop,
                        phi_prop,
                        gamma_prop,
                        gamma1_prop,
                        gamma2_prop,
                    )
                else:
                    w_filter = sigma_prop**4 * intrinsic_by_proposal
                    max_filter = np.max(w_filter)

                    cs = np.zeros(n_prop)
                    if max_filter > 0.0:
                        threshold = max_filter * threshold_factor
                        idx_filter = w_filter > threshold
                        if np.any(idx_filter):
                            cs[idx_filter] = cross_section(
                                zs_prop[idx_filter],
                                zl_prop[idx_filter],
                                sigma_prop[idx_filter],
                                q_prop[idx_filter],
                                phi_prop[idx_filter],
                                gamma_prop[idx_filter],
                                gamma1_prop[idx_filter],
                                gamma2_prop[idx_filter],
                            )

                w = np.zeros(n_prop)
                idx_weight = valid & np.isfinite(cs) & (cs > 0.0)
                w[idx_weight] = cs[idx_weight] * intrinsic_by_proposal[idx_weight] / (4.0 * np.pi)

                batch_max = np.max(w)
                if batch_max <= 0.0:
                    continue

                if batch_max > w_max:
                    if w_max > 0.0:
                        keep_prob = w_max / batch_max
                        new_count = 0

                        for j in range(count):
                            if np.random.random() < keep_prob:
                                zs_collected[new_count] = zs_collected[j]
                                zl_collected[new_count] = zl_collected[j]
                                sigma_collected[new_count] = sigma_collected[j]
                                q_collected[new_count] = q_collected[j]
                                phi_collected[new_count] = phi_collected[j]
                                new_count += 1

                        count = new_count

                    w_max = batch_max

                for j in range(n_prop):
                    if count >= size_i:
                        break

                    if np.random.random() < (w[j] / w_max):
                        zs_collected[count] = zs_prop[j]
                        zl_collected[count] = zl_prop[j]
                        sigma_collected[count] = sigma_prop[j]
                        q_collected[count] = q_prop[j]
                        phi_collected[count] = phi_prop[j]
                        count += 1

            zs_post[start:end] = zs_collected
            zl_post[start:end] = zl_collected
            sigma_post[start:end] = sigma_collected
            q_post[start:end] = q_collected
            phi_post[start:end] = phi_collected

    if lens_type in ["sis_galaxy"]:

        q_prop = np.ones(n_prop)
        phi_prop = np.zeros(n_prop)
        gamma_prop = np.ones(n_prop) * 2.0
        gamma1_prop = np.zeros(n_prop)
        gamma2_prop = np.zeros(n_prop)

        for i in prange(npool):

            start, end = _chunk_bounds(i, size, npool)

            size_i = end - start

            zs_collected = np.zeros(size_i)
            zl_collected = np.zeros(size_i)
            sigma_collected = np.zeros(size_i)

            count = 0
            w_max = 0.0

            while count < size_i:

                zs_prop = z_min + (z_max - z_min) * np.random.random(n_prop)
                zl_prop = z_min + (zs_prop - z_min) * np.random.random(n_prop)

                proposal_pdf_i = proposal_pdf / (zs_prop - z_min)

                sigma_prop = sigma_min + (sigma_max - sigma_min) * np.random.random(n_prop)

                intrinsic_pdf = (
                    zs_pdf(zs_prop)
                    * dVdz(zl_prop)
                    * number_density(sigma_prop, zl_prop)
                )

                valid = (
                    np.isfinite(proposal_pdf_i)
                    & (proposal_pdf_i > 0.0)
                    & np.isfinite(intrinsic_pdf)
                    & (intrinsic_pdf > 0.0)
                )
                intrinsic_by_proposal = np.zeros(n_prop)
                intrinsic_by_proposal[valid] = intrinsic_pdf[valid] / proposal_pdf_i[valid]

                if threshold_factor <= 0.0:
                    cs = cross_section(
                        zs_prop,
                        zl_prop,
                        sigma_prop,
                        q_prop,
                        phi_prop,
                        gamma_prop,
                        gamma1_prop,
                        gamma2_prop,
                    )
                else:
                    w_filter = sigma_prop**4 * intrinsic_by_proposal
                    max_filter = np.max(w_filter)

                    cs = np.zeros(n_prop)
                    if max_filter > 0.0:
                        threshold = max_filter * threshold_factor
                        idx_filter = w_filter > threshold
                        if np.any(idx_filter):
                            cs[idx_filter] = cross_section(
                                zs_prop[idx_filter],
                                zl_prop[idx_filter],
                                sigma_prop[idx_filter],
                                q_prop[idx_filter],
                                phi_prop[idx_filter],
                                gamma_prop[idx_filter],
                                gamma1_prop[idx_filter],
                                gamma2_prop[idx_filter],
                            )

                w = np.zeros(n_prop)
                idx_weight = valid & np.isfinite(cs) & (cs > 0.0)
                w[idx_weight] = cs[idx_weight] * intrinsic_by_proposal[idx_weight] / (4.0 * np.pi)

                batch_max = np.max(w)
                if batch_max <= 0.0:
                    continue

                if batch_max > w_max:
                    if w_max > 0.0:
                        keep_prob = w_max / batch_max
                        new_count = 0

                        for j in range(count):
                            if np.random.random() < keep_prob:
                                zs_collected[new_count] = zs_collected[j]
                                zl_collected[new_count] = zl_collected[j]
                                sigma_collected[new_count] = sigma_collected[j]
                                new_count += 1

                        count = new_count

                    w_max = batch_max

                for j in range(n_prop):
                    if count >= size_i:
                        break

                    if np.random.random() < (w[j] / w_max):
                        zs_collected[count] = zs_prop[j]
                        zl_collected[count] = zl_prop[j]
                        sigma_collected[count] = sigma_prop[j]
                        count += 1

            zs_post[start:end] = zs_collected
            zl_post[start:end] = zl_collected
            sigma_post[start:end] = sigma_collected

    return zs_post, zl_post, sigma_post, q_post, phi_post, gamma_post, gamma1_post, gamma2_post

# @njit(parallel=True, fastmath=True)
def rejection_sampler_partial(
    size,
    zs,
    zl,
    number_density,
    q_pdf,
    phi_pdf,
    gamma_pdf,
    shear1_pdf,
    shear2_pdf,
    cross_section,
    ranges,
    threshold_factor=1e-4,
    n_prop=1000,
    lens_type="epl_shear_galaxy",
    npool=4,
):
    sigma_post = np.zeros(size)
    q_post = np.ones(size)
    phi_post = np.zeros(size)
    gamma_post = np.ones(size) * 2.0
    gamma1_post = np.zeros(size)
    gamma2_post = np.zeros(size)

    z_min, z_max, sigma_min, sigma_max, q_min, q_max, phi_min, phi_max, gamma_min, gamma_max, shear_min, shear_max = ranges

    proposal_pdf = 1.0
    proposal_pdf *= _uniform_pdf(sigma_min, sigma_max)

    if lens_type in ["epl_shear_galaxy"]:

        proposal_pdf *= _uniform_pdf(q_min, q_max)
        proposal_pdf *= _uniform_pdf(phi_min, phi_max)
        proposal_pdf *= _uniform_pdf(gamma_min, gamma_max)
        proposal_pdf *= _uniform_pdf(shear_min, shear_max)
        proposal_pdf *= _uniform_pdf(shear_min, shear_max)

        for i in prange(npool):

            start, end = _chunk_bounds(i, size, npool)

            for k in range(start, end):

                zs_i = np.ones(n_prop) * zs[k]
                zl_i = np.ones(n_prop) * zl[k]

                sigma_sample = 0.0
                q_sample = 1.0
                phi_sample = 0.0
                gamma_sample = 2.0
                gamma1_sample = 0.0
                gamma2_sample = 0.0

                count = 0
                w_max = 0.0

                while count < 1:

                    sigma_prop = sigma_min + (sigma_max - sigma_min) * np.random.random(n_prop)
                    q_prop = q_min + (q_max - q_min) * np.random.random(n_prop)
                    phi_prop = phi_min + (phi_max - phi_min) * np.random.random(n_prop)
                    gamma_prop = gamma_min + (gamma_max - gamma_min) * np.random.random(n_prop)
                    gamma1_prop = shear_min + (shear_max - shear_min) * np.random.random(n_prop)
                    gamma2_prop = shear_min + (shear_max - shear_min) * np.random.random(n_prop)

                    intrinsic_pdf = (
                        number_density(sigma_prop, zl_i)
                        * q_pdf(q_prop, sigma_prop)
                        * phi_pdf(phi_prop)
                        * gamma_pdf(gamma_prop)
                        * shear1_pdf(gamma1_prop)
                        * shear2_pdf(gamma2_prop)
                    )

                    valid = (
                        np.isfinite(proposal_pdf)
                        & (proposal_pdf > 0.0)
                        & np.isfinite(intrinsic_pdf)
                        & (intrinsic_pdf > 0.0)
                    )
                    intrinsic_by_proposal = np.zeros(n_prop)
                    intrinsic_by_proposal[valid] = intrinsic_pdf[valid] / proposal_pdf

                    if threshold_factor <= 0.0:
                        cs = cross_section(
                            zs_i,
                            zl_i,
                            sigma_prop,
                            q_prop,
                            phi_prop,
                            gamma_prop,
                            gamma1_prop,
                            gamma2_prop,
                        )
                    else:
                        w_filter = sigma_prop**4 * intrinsic_by_proposal
                        max_filter = np.max(w_filter)

                        cs = np.zeros(n_prop)
                        if max_filter > 0.0:
                            threshold = max_filter * threshold_factor
                            idx_filter = w_filter > threshold
                            if np.any(idx_filter):
                                cs[idx_filter] = cross_section(
                                    zs_i[idx_filter],
                                    zl_i[idx_filter],
                                    sigma_prop[idx_filter],
                                    q_prop[idx_filter],
                                    phi_prop[idx_filter],
                                    gamma_prop[idx_filter],
                                    gamma1_prop[idx_filter],
                                    gamma2_prop[idx_filter],
                                )

                    w = np.zeros(n_prop)
                    idx_weight = valid & np.isfinite(cs) & (cs > 0.0)
                    w[idx_weight] = cs[idx_weight] * intrinsic_by_proposal[idx_weight] / (4.0 * np.pi)

                    batch_max = np.max(w)
                    if batch_max <= 0.0:
                        continue

                    if batch_max > w_max:
                        if w_max > 0.0 and count == 1:
                            if np.random.random() >= (w_max / batch_max):
                                count = 0

                        w_max = batch_max

                    for j in range(n_prop):
                        if count >= 1:
                            break

                        if np.random.random() < (w[j] / w_max):
                            sigma_sample = sigma_prop[j]
                            q_sample = q_prop[j]
                            phi_sample = phi_prop[j]
                            gamma_sample = gamma_prop[j]
                            gamma1_sample = gamma1_prop[j]
                            gamma2_sample = gamma2_prop[j]
                            count = 1

                sigma_post[k] = sigma_sample
                q_post[k] = q_sample
                phi_post[k] = phi_sample
                gamma_post[k] = gamma_sample
                gamma1_post[k] = gamma1_sample
                gamma2_post[k] = gamma2_sample

    if lens_type in ["sie_galaxy"]:

        proposal_pdf *= _uniform_pdf(q_min, q_max)
        proposal_pdf *= _uniform_pdf(phi_min, phi_max)

        gamma_prop = np.ones(n_prop) * 2.0
        gamma1_prop = np.zeros(n_prop)
        gamma2_prop = np.zeros(n_prop)

        for i in prange(npool):

            start, end = _chunk_bounds(i, size, npool)

            for k in range(start, end):

                zs_i = np.ones(n_prop) * zs[k]
                zl_i = np.ones(n_prop) * zl[k]

                sigma_sample = 0.0
                q_sample = 1.0
                phi_sample = 0.0

                count = 0
                w_max = 0.0

                while count < 1:

                    sigma_prop = sigma_min + (sigma_max - sigma_min) * np.random.random(n_prop)
                    q_prop = q_min + (q_max - q_min) * np.random.random(n_prop)
                    phi_prop = phi_min + (phi_max - phi_min) * np.random.random(n_prop)

                    intrinsic_pdf = (
                        number_density(sigma_prop, zl_i)
                        * q_pdf(q_prop, sigma_prop)
                        * phi_pdf(phi_prop)
                    )

                    valid = (
                        np.isfinite(proposal_pdf)
                        & (proposal_pdf > 0.0)
                        & np.isfinite(intrinsic_pdf)
                        & (intrinsic_pdf > 0.0)
                    )
                    intrinsic_by_proposal = np.zeros(n_prop)
                    intrinsic_by_proposal[valid] = intrinsic_pdf[valid] / proposal_pdf

                    if threshold_factor <= 0.0:
                        cs = cross_section(
                            zs_i,
                            zl_i,
                            sigma_prop,
                            q_prop,
                            phi_prop,
                            gamma_prop,
                            gamma1_prop,
                            gamma2_prop,
                        )
                    else:
                        w_filter = sigma_prop**4 * intrinsic_by_proposal
                        max_filter = np.max(w_filter)

                        cs = np.zeros(n_prop)
                        if max_filter > 0.0:
                            threshold = max_filter * threshold_factor
                            idx_filter = w_filter > threshold
                            if np.any(idx_filter):
                                cs[idx_filter] = cross_section(
                                    zs_i[idx_filter],
                                    zl_i[idx_filter],
                                    sigma_prop[idx_filter],
                                    q_prop[idx_filter],
                                    phi_prop[idx_filter],
                                    gamma_prop[idx_filter],
                                    gamma1_prop[idx_filter],
                                    gamma2_prop[idx_filter],
                                )

                    w = np.zeros(n_prop)
                    idx_weight = valid & np.isfinite(cs) & (cs > 0.0)
                    w[idx_weight] = cs[idx_weight] * intrinsic_by_proposal[idx_weight] / (4.0 * np.pi)

                    batch_max = np.max(w)
                    if batch_max <= 0.0:
                        continue

                    if batch_max > w_max:
                        if w_max > 0.0 and count == 1:
                            if np.random.random() >= (w_max / batch_max):
                                count = 0

                        w_max = batch_max

                    for j in range(n_prop):
                        if count >= 1:
                            break

                        if np.random.random() < (w[j] / w_max):
                            sigma_sample = sigma_prop[j]
                            q_sample = q_prop[j]
                            phi_sample = phi_prop[j]
                            count = 1

                sigma_post[k] = sigma_sample
                q_post[k] = q_sample
                phi_post[k] = phi_sample

    if lens_type in ["sis_galaxy"]:

        q_prop = np.ones(n_prop)
        phi_prop = np.zeros(n_prop)
        gamma_prop = np.ones(n_prop) * 2.0
        gamma1_prop = np.zeros(n_prop)
        gamma2_prop = np.zeros(n_prop)

        for i in prange(npool):

            start, end = _chunk_bounds(i, size, npool)

            for k in range(start, end):

                zs_i = np.ones(n_prop) * zs[k]
                zl_i = np.ones(n_prop) * zl[k]

                sigma_sample = 0.0

                count = 0
                w_max = 0.0

                while count < 1:

                    sigma_prop = sigma_min + (sigma_max - sigma_min) * np.random.random(n_prop)

                    intrinsic_pdf = number_density(sigma_prop, zl_i)

                    valid = (
                        np.isfinite(proposal_pdf)
                        & (proposal_pdf > 0.0)
                        & np.isfinite(intrinsic_pdf)
                        & (intrinsic_pdf > 0.0)
                    )
                    intrinsic_by_proposal = np.zeros(n_prop)
                    intrinsic_by_proposal[valid] = intrinsic_pdf[valid] / proposal_pdf

                    w_filter = sigma_prop**4 * intrinsic_by_proposal
                    max_filter = np.max(w_filter)

                    cs = np.zeros(n_prop)
                    if max_filter > 0.0:
                        threshold = max_filter * threshold_factor
                        idx_filter = w_filter > threshold
                        if np.any(idx_filter):
                            cs[idx_filter] = cross_section(
                                zs_i[idx_filter],
                                zl_i[idx_filter],
                                sigma_prop[idx_filter],
                                q_prop[idx_filter],
                                phi_prop[idx_filter],
                                gamma_prop[idx_filter],
                                gamma1_prop[idx_filter],
                                gamma2_prop[idx_filter],
                            )

                    w = np.zeros(n_prop)
                    idx_weight = valid & np.isfinite(cs) & (cs > 0.0)
                    w[idx_weight] = cs[idx_weight] * intrinsic_by_proposal[idx_weight] / (4.0 * np.pi)

                    batch_max = np.max(w)
                    if batch_max <= 0.0:
                        continue

                    if batch_max > w_max:
                        if w_max > 0.0 and count == 1:
                            if np.random.random() >= (w_max / batch_max):
                                count = 0

                        w_max = batch_max

                    for j in range(n_prop):
                        if count >= 1:
                            break

                        if np.random.random() < (w[j] / w_max):
                            sigma_sample = sigma_prop[j]
                            count = 1

                sigma_post[k] = sigma_sample

    return sigma_post, q_post, phi_post, gamma_post, gamma1_post, gamma2_post

# -------------------------------------------------
# Main importance sampler function
# -------------------------------------------------
# @njit(parallel=True, fastmath=True)
def importance_sampler_partial(
    size,
    zs,
    zl,
    number_density,
    q_pdf,
    phi_pdf,
    gamma_pdf,
    shear1_pdf,
    shear2_pdf,
    cross_section,
    ranges,
    threshold_factor=1e-4,
    n_prop=100, 
    lens_type="epl_shear_galaxy", 
):
    """
    Importance/resampling for lens nuisance parameters at fixed (zs, zl).

    Redshifts are fixed (typically drawn from ``zs_sl`` and ``lens_redshift_sl``).
    Conditional weights target ``p(sigma,lambda|zs,zl,SL)`` up to a constant:
    proportional to ``n(sigma,zl)``, axis-ratio / shear / slope PDFs, and the
    strong-lensing cross section. Terms that depend only on ``(zs,zl)`` ---
    e.g. ``P(zs|SL)``, ``P(zl|zs,SL)``, ``dV/dz_l`` --- are omitted (constant
    for fixed redshifts and identical scale across sigma proposals).
    """

    sigma_post = np.zeros(size)
    q_post = np.ones(size)
    phi_post = np.zeros(size)
    gamma_post = np.ones(size) * 2.0
    gamma1_post = np.zeros(size)
    gamma2_post = np.zeros(size)

    z_min, z_max, sigma_min, sigma_max, q_min, q_max, phi_min, phi_max, gamma_min, gamma_max, shear_min, shear_max = ranges

    # wi = (cs/4*pi) * (number_density(sigma, zl) * shape PDFs) / proposal_pdf

    proposal_pdf = np.ones(n_prop)
    # sigma
    proposal_pdf *= _uniform_pdf(sigma_min, sigma_max)

    if lens_type in ["epl_shear_galaxy"]:
        # q
        proposal_pdf *= _uniform_pdf(q_min, q_max)
        # phi
        proposal_pdf *= _uniform_pdf(phi_min, phi_max)
        # gamma
        proposal_pdf *= _uniform_pdf(gamma_min, gamma_max)
        # shear1
        proposal_pdf *= _uniform_pdf(shear_min, shear_max)
        # shear2
        proposal_pdf *= _uniform_pdf(shear_min, shear_max)

        for i in prange(size):  # for each (zl, zs) pair

            zs_i = np.ones(n_prop) * zs[i]
            zl_i = np.ones(n_prop) * zl[i]
            
            # Draw proposals from uniform distribution
            sigma_prop = sigma_min + (sigma_max - sigma_min) * np.random.random(n_prop)
            q_prop = q_min + (q_max - q_min) * np.random.random(n_prop)
            phi_prop = phi_min + (phi_max - phi_min) * np.random.random(n_prop)
            gamma_prop = gamma_min + (gamma_max - gamma_min) * np.random.random(n_prop)
            gamma1_prop = shear_min + (shear_max - shear_min) * np.random.random(n_prop)
            gamma2_prop = shear_min + (shear_max - shear_min) * np.random.random(n_prop)

            # Target factors that vary with proposals (n(sigma,zl) and shape PDFs)
            intrinsic_pdf = (
                number_density(sigma_prop, zl_i)
                * q_pdf(q_prop, sigma_prop)
                * phi_pdf(phi_prop)
                * gamma_pdf(gamma_prop)
                * shear1_pdf(gamma1_prop)
                * shear2_pdf(gamma2_prop)
            )

            # Early filtering: only compute cross-section for proposals with significant intrinsic PDF
            # This can reduce cross-section calls by 50-90% with negligible effect on ESS

            # filter out proposals with based on sis cross-section (\propto sigma_prop**4) and intrinsic pdf, which are cheap to compute
            valid = (
                np.isfinite(proposal_pdf)
                & (proposal_pdf > 0.0)
                & np.isfinite(intrinsic_pdf)
                & (intrinsic_pdf > 0.0)
            )
            intrinsic_by_proposal = np.zeros(n_prop)
            intrinsic_by_proposal[valid] = intrinsic_pdf[valid] / proposal_pdf[valid]

            if threshold_factor <= 0.0:
                cs = cross_section(
                    zs_i,
                    zl_i,
                    sigma_prop,
                    q_prop,
                    phi_prop,
                    gamma_prop,
                    gamma1_prop,
                    gamma2_prop,
                )
            else:
                w_filter = sigma_prop**4 * intrinsic_by_proposal
                # w_filter[(~np.isfinite(w_filter)) | (w_filter < 0.0)] = 0.0
                max_filter = np.max(w_filter)

                cs = np.zeros(n_prop)
                if max_filter > 0.0:
                    threshold = max_filter * threshold_factor
                    idx_filter = w_filter > threshold

                    if np.any(idx_filter):
                        cs[idx_filter] = cross_section(
                            zs_i[idx_filter],
                            zl_i[idx_filter],
                            sigma_prop[idx_filter],
                            q_prop[idx_filter],
                            phi_prop[idx_filter],
                            gamma_prop[idx_filter],
                            gamma1_prop[idx_filter],
                            gamma2_prop[idx_filter],
                        )

            w = np.zeros(n_prop)
            idx_weight = valid & np.isfinite(cs) & (cs > 0.0)

            w[idx_weight] = (
                cs[idx_weight]
                * intrinsic_by_proposal[idx_weight]
                / (4.0 * np.pi)
            )

            idx = _weighted_sample_index(w, n_prop)

            sigma_post[i] = sigma_prop[idx]
            q_post[i] = q_prop[idx]
            phi_post[i] = phi_prop[idx]
            gamma_post[i] = gamma_prop[idx]
            gamma1_post[i] = gamma1_prop[idx]
            gamma2_post[i] = gamma2_prop[idx]

    if lens_type in ["sie_galaxy"]:

        # q
        proposal_pdf *= _uniform_pdf(q_min, q_max)
        # phi
        proposal_pdf *= _uniform_pdf(phi_min, phi_max)

        gamma_prop = np.ones(n_prop) * 2.0
        gamma1_prop = np.zeros(n_prop)
        gamma2_prop = np.zeros(n_prop)

        for i in prange(size):  # for each (zl, zs) pair

            zs_i = np.ones(n_prop) * zs[i]
            zl_i = np.ones(n_prop) * zl[i]
            
            # Draw proposals from uniform distribution
            sigma_prop = sigma_min + (sigma_max - sigma_min) * np.random.random(n_prop)
            q_prop = q_min + (q_max - q_min) * np.random.random(n_prop)
            phi_prop = phi_min + (phi_max - phi_min) * np.random.random(n_prop)

            intrinsic_pdf = (
                number_density(sigma_prop, zl_i)
                * q_pdf(q_prop, sigma_prop)
                * phi_pdf(phi_prop)
            )

            valid = (
                np.isfinite(proposal_pdf)
                & (proposal_pdf > 0.0)
                & np.isfinite(intrinsic_pdf)
                & (intrinsic_pdf > 0.0)
            )
            intrinsic_by_proposal = np.zeros(n_prop)
            intrinsic_by_proposal[valid] = intrinsic_pdf[valid] / proposal_pdf[valid]

            if threshold_factor <= 0.0:
                cs = cross_section(
                    zs_i,
                    zl_i,
                    sigma_prop,
                    q_prop,
                    phi_prop,
                    gamma_prop,
                    gamma1_prop,
                    gamma2_prop,
                )
            else:
                w_filter = sigma_prop**4 * intrinsic_by_proposal
                max_filter = np.max(w_filter)

                cs = np.zeros(n_prop)
                if max_filter > 0.0:
                    threshold = max_filter * threshold_factor
                    idx_filter = w_filter > threshold

                    if np.any(idx_filter):
                        cs[idx_filter] = cross_section(
                            zs_i[idx_filter],
                            zl_i[idx_filter],
                            sigma_prop[idx_filter],
                            q_prop[idx_filter],
                            phi_prop[idx_filter],
                            gamma_prop[idx_filter],
                            gamma1_prop[idx_filter],
                            gamma2_prop[idx_filter],
                        )

            # Compute importance weights
            w = np.zeros(n_prop)
            idx_weight = valid & np.isfinite(cs) & (cs > 0.0)
            w[idx_weight] = cs[idx_weight] * intrinsic_by_proposal[idx_weight] / (4.0 * np.pi)
            # Draw posterior sample via weighted choice
            idx = _weighted_sample_index(w, n_prop)

            sigma_post[i] = sigma_prop[idx]
            q_post[i] = q_prop[idx]
            phi_post[i] = phi_prop[idx]

    if lens_type in ["sis_galaxy"]:

        q_prop = np.ones(n_prop)
        phi_prop = np.zeros(n_prop)
        gamma_prop = np.ones(n_prop) * 2.0
        gamma1_prop = np.zeros(n_prop)
        gamma2_prop = np.zeros(n_prop)

        for i in prange(size):  # for each (zl, zs) pair

            zs_i = np.ones(n_prop) * zs[i]
            zl_i = np.ones(n_prop) * zl[i]
            
            # Draw proposals from uniform distribution
            sigma_prop = sigma_min + (sigma_max - sigma_min) * np.random.random(n_prop)

            intrinsic_pdf = number_density(sigma_prop, zl_i)

            valid = (
                np.isfinite(proposal_pdf)
                & (proposal_pdf > 0.0)
                & np.isfinite(intrinsic_pdf)
                & (intrinsic_pdf > 0.0)
            )
            intrinsic_by_proposal = np.zeros(n_prop)
            intrinsic_by_proposal[valid] = intrinsic_pdf[valid] / proposal_pdf[valid]

            if threshold_factor <= 0.0:
                cs = cross_section(
                    zs_i,
                    zl_i,
                    sigma_prop,
                    q_prop,
                    phi_prop,
                    gamma_prop,
                    gamma1_prop,
                    gamma2_prop,
                )
            else:
                w_filter = sigma_prop**4 * intrinsic_by_proposal
                max_filter = np.max(w_filter)

                cs = np.zeros(n_prop)
                if max_filter > 0.0:
                    threshold = max_filter * threshold_factor
                    idx_filter = w_filter > threshold

                    if np.any(idx_filter):
                        cs[idx_filter] = cross_section(
                            zs_i[idx_filter],
                            zl_i[idx_filter],
                            sigma_prop[idx_filter],
                            q_prop[idx_filter],
                            phi_prop[idx_filter],
                            gamma_prop[idx_filter],
                            gamma1_prop[idx_filter],
                            gamma2_prop[idx_filter],
                        )

            # Compute importance weights
            w = np.zeros(n_prop)
            idx_weight = valid & np.isfinite(cs) & (cs > 0.0)
            w[idx_weight] = cs[idx_weight] * intrinsic_by_proposal[idx_weight] / (4.0 * np.pi)
            # Draw posterior sample via weighted choice
            idx = _weighted_sample_index(w, n_prop)

            sigma_post[i] = sigma_prop[idx]

    return sigma_post, q_post, phi_post, gamma_post, gamma1_post, gamma2_post

# @njit(parallel=True, fastmath=True)
def importance_sampler_full(
    size,
    zs_pdf,
    number_density,
    q_pdf,
    phi_pdf,
    gamma_pdf,
    shear1_pdf,
    shear2_pdf,
    cross_section,
    dVdz,
    ranges,
    threshold_factor=1e-4,
    n_prop=50, 
    lens_type="epl_shear_galaxy", 
):
    """
    Full joint importance sample over (zs, zl, sigma, nuisance).

    Target for lens+source intrinsic factors (before cross section) is taken as

        p(zs) * (dV_c/dz_l) * n(sigma, zl) * (axis ratio, shear, slope PDFs),

    i.e. **no separate** ``p(zl|zs)`` factor: ``z_l`` is proposed uniformly in
    ``[z_min, z_s)`` and the comoving volume and Schechter-style ``n(sigma,zl)``
    carry the lens-redshift abundance. This matches
    ``P(SL|...) p(zs) n(sigma,zl|...) dV/dz_l ...`` up to overall normalisation.
    """

    zs_post = np.zeros(size)
    zl_post = np.zeros(size)
    sigma_post = np.zeros(size)
    q_post = np.ones(size)
    phi_post = np.zeros(size)
    gamma_post = np.ones(size) * 2.0
    gamma1_post = np.zeros(size)
    gamma2_post = np.zeros(size)

    z_min, z_max, sigma_min, sigma_max, q_min, q_max, phi_min, phi_max, gamma_min, gamma_max, shear_min, shear_max = ranges

    # wi = (cs/4*pi) * ( zs_pdf(zs) * dVdz(zl) * number_density(sigma, zl) * ... ) / proposal_pdf

    proposal_pdf = np.ones(n_prop)
    # zs
    proposal_pdf *= _uniform_pdf(z_min, z_max)
    # sigma
    proposal_pdf *= _uniform_pdf(sigma_min, sigma_max)

    if lens_type in ["epl_shear_galaxy"]:
        # q
        proposal_pdf *= _uniform_pdf(q_min, q_max)
        # phi
        proposal_pdf *= _uniform_pdf(phi_min, phi_max)
        # gamma
        proposal_pdf *= _uniform_pdf(gamma_min, gamma_max)
        # shear1
        proposal_pdf *= _uniform_pdf(shear_min, shear_max)
        # shear2
        proposal_pdf *= _uniform_pdf(shear_min, shear_max)

        for i in prange(size):  # for each (zl, zs) pair

            # Draw proposals from uniform distribution
            zs_prop = z_min + (z_max - z_min) * np.random.random(n_prop)

            # zl
            zl_prop = z_min + (zs_prop - z_min) * np.random.random(n_prop)
            proposal_pdf_i = proposal_pdf / (zs_prop - z_min)

            sigma_prop = sigma_min + (sigma_max - sigma_min) * np.random.random(n_prop)
            q_prop = q_min + (q_max - q_min) * np.random.random(n_prop)
            phi_prop = phi_min + (phi_max - phi_min) * np.random.random(n_prop)
            gamma_prop = gamma_min + (gamma_max - gamma_min) * np.random.random(n_prop)
            gamma1_prop = shear_min + (shear_max - shear_min) * np.random.random(n_prop)
            gamma2_prop = shear_min + (shear_max - shear_min) * np.random.random(n_prop)

            # intrinsic pdfs
            intrinsic_pdf = (zs_pdf(zs_prop) * dVdz(zl_prop) * number_density(sigma_prop, zl_prop) * q_pdf(q_prop, sigma_prop) * phi_pdf(phi_prop) * gamma_pdf(gamma_prop) * shear1_pdf(gamma1_prop) * shear2_pdf(gamma2_prop))

            valid = (
                np.isfinite(proposal_pdf_i)
                & (proposal_pdf_i > 0.0)
                & np.isfinite(intrinsic_pdf)
                & (intrinsic_pdf > 0.0)
            )
            intrinsic_by_proposal = np.zeros(n_prop)
            intrinsic_by_proposal[valid] = intrinsic_pdf[valid] / proposal_pdf_i[valid]

            if threshold_factor <= 0.0:
                cs = cross_section(
                    zs_prop,
                    zl_prop,
                    sigma_prop,
                    q_prop,
                    phi_prop,
                    gamma_prop,
                    gamma1_prop,
                    gamma2_prop,
                )
            else:
                w_filter = sigma_prop**4 * intrinsic_by_proposal
                max_filter = np.max(w_filter)

                cs = np.zeros(n_prop)
                if max_filter > 0.0:
                    threshold = max_filter * threshold_factor
                    idx_filter = w_filter > threshold
                    if np.any(idx_filter):
                        cs[idx_filter] = cross_section(
                            zs_prop[idx_filter],
                            zl_prop[idx_filter],
                            sigma_prop[idx_filter],
                            q_prop[idx_filter],
                            phi_prop[idx_filter],
                            gamma_prop[idx_filter],
                            gamma1_prop[idx_filter],
                            gamma2_prop[idx_filter],
                        )

            # Compute importance weights
            w = np.zeros(n_prop)
            idx_weight = valid & np.isfinite(cs) & (cs > 0.0)
            w[idx_weight] = cs[idx_weight] * intrinsic_by_proposal[idx_weight] / (4.0 * np.pi)
            # Draw posterior sample via weighted choice
            idx = _weighted_sample_index(w, n_prop)

            zs_post[i] = zs_prop[idx]
            zl_post[i] = zl_prop[idx]
            sigma_post[i] = sigma_prop[idx]
            q_post[i] = q_prop[idx]
            phi_post[i] = phi_prop[idx]
            gamma_post[i] = gamma_prop[idx]
            gamma1_post[i] = gamma1_prop[idx]
            gamma2_post[i] = gamma2_prop[idx]

    if lens_type in ["sie_galaxy"]:

        # q
        proposal_pdf *= _uniform_pdf(q_min, q_max)
        # phi
        proposal_pdf *= _uniform_pdf(phi_min, phi_max)

        gamma_prop = np.ones(n_prop) * 2.0
        gamma1_prop = np.zeros(n_prop)
        gamma2_prop = np.zeros(n_prop)

        for i in prange(size):  # for each (zl, zs) pair

            # Draw proposals from uniform distribution
            zs_prop = z_min + (z_max - z_min) * np.random.random(n_prop)

            # zl
            zl_prop = z_min + (zs_prop - z_min) * np.random.random(n_prop)
            proposal_pdf_i = proposal_pdf / (zs_prop - z_min)

            sigma_prop = sigma_min + (sigma_max - sigma_min) * np.random.random(n_prop)
            q_prop = q_min + (q_max - q_min) * np.random.random(n_prop)
            phi_prop = phi_min + (phi_max - phi_min) * np.random.random(n_prop)

            # intrinsic pdfs
            intrinsic_pdf = (zs_pdf(zs_prop) * dVdz(zl_prop) * number_density(sigma_prop, zl_prop) * q_pdf(q_prop, sigma_prop) * phi_pdf(phi_prop) )

            valid = (
                np.isfinite(proposal_pdf_i)
                & (proposal_pdf_i > 0.0)
                & np.isfinite(intrinsic_pdf)
                & (intrinsic_pdf > 0.0)
            )
            intrinsic_by_proposal = np.zeros(n_prop)
            intrinsic_by_proposal[valid] = intrinsic_pdf[valid] / proposal_pdf_i[valid]

            if threshold_factor <= 0.0:
                cs = cross_section(
                    zs_prop,
                    zl_prop,
                    sigma_prop,
                    q_prop,
                    phi_prop,
                    gamma_prop,
                    gamma1_prop,
                    gamma2_prop,
                )
            else:
                w_filter = sigma_prop**4 * intrinsic_by_proposal
                max_filter = np.max(w_filter)

                cs = np.zeros(n_prop)
                if max_filter > 0.0:
                    threshold = max_filter * threshold_factor
                    idx_filter = w_filter > threshold
                    if np.any(idx_filter):
                        cs[idx_filter] = cross_section(
                            zs_prop[idx_filter],
                            zl_prop[idx_filter],
                            sigma_prop[idx_filter],
                            q_prop[idx_filter],
                            phi_prop[idx_filter],
                            gamma_prop[idx_filter],
                            gamma1_prop[idx_filter],
                            gamma2_prop[idx_filter],
                        )

            # Compute importance weights
            w = np.zeros(n_prop)
            idx_weight = valid & np.isfinite(cs) & (cs > 0.0)
            w[idx_weight] = cs[idx_weight] * intrinsic_by_proposal[idx_weight] / (4.0 * np.pi)
            # Draw posterior sample via weighted choice
            idx = _weighted_sample_index(w, n_prop)

            zs_post[i] = zs_prop[idx]
            zl_post[i] = zl_prop[idx]
            sigma_post[i] = sigma_prop[idx]
            q_post[i] = q_prop[idx]
            phi_post[i] = phi_prop[idx]

    if lens_type in ["sis_galaxy"]:

        q_prop = np.ones(n_prop)
        phi_prop = np.zeros(n_prop)
        gamma_prop = np.ones(n_prop) * 2.0
        gamma1_prop = np.zeros(n_prop)
        gamma2_prop = np.zeros(n_prop)

        for i in prange(size):  # for each (zl, zs) pair

            # Draw proposals from uniform distribution
            zs_prop = z_min + (z_max - z_min) * np.random.random(n_prop)

            # zl
            zl_prop = z_min + (zs_prop - z_min) * np.random.random(n_prop)
            proposal_pdf_i = proposal_pdf / (zs_prop - z_min)

            sigma_prop = sigma_min + (sigma_max - sigma_min) * np.random.random(n_prop)

            # intrinsic pdfs
            intrinsic_pdf = (zs_pdf(zs_prop) * dVdz(zl_prop) * number_density(sigma_prop, zl_prop) )

            valid = (
                np.isfinite(proposal_pdf_i)
                & (proposal_pdf_i > 0.0)
                & np.isfinite(intrinsic_pdf)
                & (intrinsic_pdf > 0.0)
            )
            intrinsic_by_proposal = np.zeros(n_prop)
            intrinsic_by_proposal[valid] = intrinsic_pdf[valid] / proposal_pdf_i[valid]

            if threshold_factor <= 0.0:
                cs = cross_section(
                    zs_prop,
                    zl_prop,
                    sigma_prop,
                    q_prop,
                    phi_prop,
                    gamma_prop,
                    gamma1_prop,
                    gamma2_prop,
                )
            else:
                w_filter = sigma_prop**4 * intrinsic_by_proposal
                max_filter = np.max(w_filter)

                cs = np.zeros(n_prop)
                if max_filter > 0.0:
                    threshold = max_filter * threshold_factor
                    idx_filter = w_filter > threshold
                    if np.any(idx_filter):
                        cs[idx_filter] = cross_section(
                            zs_prop[idx_filter],
                            zl_prop[idx_filter],
                            sigma_prop[idx_filter],
                            q_prop[idx_filter],
                            phi_prop[idx_filter],
                            gamma_prop[idx_filter],
                            gamma1_prop[idx_filter],
                            gamma2_prop[idx_filter],
                        )

            # Compute importance weights
            w = np.zeros(n_prop)
            idx_weight = valid & np.isfinite(cs) & (cs > 0.0)
            w[idx_weight] = cs[idx_weight] * intrinsic_by_proposal[idx_weight] / (4.0 * np.pi)
            # Draw posterior sample via weighted choice
            idx = _weighted_sample_index(w, n_prop)

            zs_post[i] = zs_prop[idx]
            zl_post[i] = zl_prop[idx]
            sigma_post[i] = sigma_prop[idx]

    return zs_post, zl_post, sigma_post, q_post, phi_post, gamma_post, gamma1_post, gamma2_post

def _range_dict(
    z_min=0.0,
    z_max=10.0,
    sigma_min=100.0,
    sigma_max=400.0,
    q_min=0.2,
    q_max=1.0,
    phi_min=0.0,
    phi_max=2 * np.pi,
    gamma_min=1.5,
    gamma_max=2.5,
    shear_min=-0.2,
    shear_max=0.2,
):
    return (
        z_min,
        z_max,
        sigma_min,
        sigma_max,
        q_min,
        q_max,
        phi_min,
        phi_max,
        gamma_min,
        gamma_max,
        shear_min,
        shear_max,
    )

importance_sampler_partial_njit = njit(parallel=True, fastmath=True)(importance_sampler_partial)
importance_sampler_full_njit = njit(parallel=True, fastmath=True)(importance_sampler_full)
rejection_sampler_partial_njit = njit(parallel=True, fastmath=True)(rejection_sampler_partial)
rejection_sampler_full_njit = njit(parallel=True, fastmath=True)(rejection_sampler_full)

def create_sampler(
    number_density,
    q_pdf,
    phi_pdf,
    gamma_pdf,
    shear1_pdf,
    shear2_pdf,
    cross_section,
    dVdz,
    zs_pdf=None,
    use_njit_sampler=True,
    sampler_type="importance_sampler_partial",
    threshold_factor=1e-4,
    n_prop=50, 
    lens_type="epl_shear_galaxy",
    npool=4, 
    **range_kwargs
):
    """
    Build a closure that runs one of the cross-section-weighted lens samplers.

    ``use_njit_sampler`` selects the implementation, not only a print message:

    - ``True`` (default): Numba ``njit`` kernels. ``cross_section`` and other
      callables must be Numba-compatible (e.g. ``@njit`` cross sections).
    - ``False``: interpreted Python samplers, which accept ordinary callables
      (e.g. numerical / SciPy-backed cross sections).
    """
    set_num_threads(npool) # Centralize Numba thread setting

    # min max
    ranges = _range_dict(**range_kwargs)

    if sampler_type == "importance_sampler_partial":

        _imp_partial = (
            importance_sampler_partial_njit
            if use_njit_sampler
            else importance_sampler_partial
        )

        def sampler_wrapper(size, zs, zl):
            return _imp_partial(
                size=size,
                zs=zs,
                zl=zl,
                threshold_factor=threshold_factor,
                n_prop=n_prop,
                lens_type=lens_type,
                number_density=number_density,
                q_pdf=q_pdf,
                phi_pdf=phi_pdf,
                gamma_pdf=gamma_pdf,
                shear1_pdf=shear1_pdf,
                shear2_pdf=shear2_pdf,
                cross_section=cross_section,
                ranges=ranges,
            )

        if use_njit_sampler:
            print(
                "Faster, njitted and importance sampling based lens parameter sampler will be used."
            )
        else: 
            print(
                "Slower, non-njit and importance sampling based lens parameter sampler will be used."
            )

    elif sampler_type == "importance_sampler_full":

        _imp_full = (
            importance_sampler_full_njit
            if use_njit_sampler
            else importance_sampler_full
        )

        def sampler_wrapper(size):
            return _imp_full(
                size=size,
                zs_pdf=zs_pdf,
                threshold_factor=threshold_factor,
                n_prop=n_prop,
                lens_type=lens_type,
                number_density=number_density,
                q_pdf=q_pdf,
                phi_pdf=phi_pdf,
                gamma_pdf=gamma_pdf,
                shear1_pdf=shear1_pdf,
                shear2_pdf=shear2_pdf,
                cross_section=cross_section,
                dVdz=dVdz,
                ranges=ranges,
            )

        if use_njit_sampler:
            print(
                "Faster, njitted and importance sampling based lens parameter sampler will be used."
            )
        else:
            print(
                "Slower, non-njit and importance sampling based lens parameter sampler will be used."
            )

    elif sampler_type == "rejection_sampler_partial":

        _rej_partial = (
            rejection_sampler_partial_njit
            if use_njit_sampler
            else rejection_sampler_partial
        )

        def sampler_wrapper(size, zs, zl):
            return _rej_partial(
                size=size,
                zs=zs,
                zl=zl,
                threshold_factor=threshold_factor,
                n_prop=n_prop,
                lens_type=lens_type,
                npool=npool,
                number_density=number_density,
                q_pdf=q_pdf,
                phi_pdf=phi_pdf,
                gamma_pdf=gamma_pdf,
                shear1_pdf=shear1_pdf,
                shear2_pdf=shear2_pdf,
                cross_section=cross_section,
                ranges=ranges,
            )

        if use_njit_sampler:
            print(
                "Faster, njitted and rejection sampling based lens parameter sampler will be used."
            )
        else:
            print(
                "Slower, non-njit and rejection sampling based lens parameter sampler will be used."
            )

    elif sampler_type == "rejection_sampler_full":

        _rej_full = (
            rejection_sampler_full_njit
            if use_njit_sampler
            else rejection_sampler_full
        )

        def sampler_wrapper(size):
            return _rej_full(
                size=size,
                zs_pdf=zs_pdf,
                threshold_factor=threshold_factor,
                n_prop=n_prop,
                lens_type=lens_type,
                npool=npool,
                number_density=number_density,
                q_pdf=q_pdf,
                phi_pdf=phi_pdf,
                gamma_pdf=gamma_pdf,
                shear1_pdf=shear1_pdf,
                shear2_pdf=shear2_pdf,
                cross_section=cross_section,
                dVdz=dVdz,
                ranges=ranges,
            )

        if use_njit_sampler:
            print(
                "Faster, njitted and rejection sampling based lens parameter sampler will be used."
            )
        else:
            print(
                "Slower, non-njit and rejection sampling based lens parameter sampler will be used."
            )

    else:
        raise ValueError(
            f"Invalid cross_section_based_sampler: {sampler_type}. Available options are 'rejection_sampler_full', 'importance_sampler_full', 'rejection_sampler_partial' and 'importance_sampler_partial'."
        )

    return sampler_wrapper
