# -*- coding: utf-8 -*-
"""
Broken Power-Law with Two Gaussian Peaks Mass Distribution Model.

This module implements a primary mass distribution model combining a broken
power-law component with two Gaussian peaks, as detailed in the accompanying
paper. The model is parameterized with 12 parameters for the primary mass and
3 parameters for the mass ratio distribution.

Key Features:
- Broken power-law distribution with configurable break mass and power-law indices
- Two Gaussian peaks for modeling specific mass populations
- Smoothing taper function for physical support enforcement
- Efficient numerical integration using log-space methods
- Full support for both PDF evaluation and random sample generation
- Mass ratio distribution with power-law form

Usage:
    Basic workflow for sampling binary masses:

    >>> from ler.gw_source_population.broken_powerlaw_plus_2peaks import broken_powerlaw_plus_2peaks_rvs
    >>> m1, m2 = broken_powerlaw_plus_2peaks_rvs(size=1000)

Copyright (C) 2026 Author Name. Distributed under MIT License.
"""

import numpy as np
from numba import njit, prange
from .prior_functions import _smoothing_S, powerlaw_with_smoothing, truncated_normal_pdf
from ..utils import inverse_transform_sampler, cumulative_trapezoid

# GWTC-4 : # lam_0=0.361, lam_1=0.586, mpp_1=9.764, sigpp_1=0.649, mpp_2=32.763, sigpp_2=3.918, mlow_1=5.059, delta_m_1=4.321, break_mass=35.622, alpha_1=1.728, alpha_2=4.512, mmax=300.0, beta=1.171, mlow_2=3.551, delta_m_2=4.910, normalization_size=500

# ---------------------------------
# Primary Mass Distribution
# ---------------------------------

@njit(fastmath=True)
def _broken_powerlaw(
    m1, break_mass=35.622, alpha_1=1.728, alpha_2=4.512, mlow_1=5.059, mmax=300.0
):
    """
    Compute unnormalized broken power-law mass distribution (equation B18).

    Evaluates the unnormalized broken power-law following the form:
    $(1/N)(m_1 / m_{break})^{-α_1}$ for $m_{low,1} ≤ m_1 < m_{break}$
    $(1/N)(m_1 / m_{break})^{-α_2}$ for $m_{break} ≤ m_1 < m_{max}$
    Support is enforced by returning 0 outside the valid interval.

    Parameters
    ----------
    m1 : ``numpy.ndarray``
        Primary masses (solar masses).
    break_mass : ``float``
        Mass at which the power-law index transitions.
        default: 35.622
    alpha_1 : ``float``
        Power-law index below the break mass.
        default: 1.728
    alpha_2 : ``float``
        Power-law index above the break mass.
        default: 4.512
    mlow_1 : ``float``
        Minimum primary mass (solar masses).
        default: 5.059
    mmax : ``float``
        Maximum primary mass (solar masses).
        default: 300.0

    Returns
    -------
    result : ``numpy.ndarray``
        Unnormalized probability density at each mass value.
    """

    result = np.zeros_like(m1)
    N = break_mass * (
        ((1 - (mlow_1 / break_mass) ** (1 - alpha_1)) / (1 - alpha_1))
        + (((mmax / break_mass) ** (1 - alpha_2) - 1) / (1 - alpha_2))
    )

    low_mask = (m1 >= mlow_1) & (m1 < break_mass)
    high_mask = (m1 >= break_mass) & (m1 < mmax)
    result[low_mask] = (m1[low_mask] / break_mass) ** (-alpha_1) / N
    result[high_mask] = (m1[high_mask] / break_mass) ** (-alpha_2) / N

    return result

# @njit
# def left_truncated_normal_pdf(m1, mu, sigma, mlow):
#     """
#     Compute left-truncated normal probability density function (equation B20).

#     Evaluates the left-truncated normal distribution $N_{lt}(m_1 | μ, σ, m_{low})$,
#     which is a Gaussian distribution with support only above a minimum value.

#     Parameters
#     ----------
#     m1 : ``numpy.ndarray``
#         Primary masses (solar masses).
#     mu : ``float``
#         Mean of the Gaussian distribution.
#     sigma : ``float``
#         Standard deviation of the Gaussian distribution.
#     mlow : ``float``
#         Minimum mass (left truncation point).

#     Returns
#     -------
#     pdf : ``numpy.ndarray``
#         Probability density values, 0 for $m_1 < m_{low}$.
#     """
#     a = (mlow - mu) / sigma
#     tail_norm = 0.5 * (1.0 - math.erf(a / math.sqrt(2.0)))
#     z = (m1 - mu) / sigma
#     pdf = np.exp(-0.5 * z * z) / (sigma * math.sqrt(2.0 * math.pi) * tail_norm)
#     return np.where(m1 >= mlow, pdf, 0.0)

@njit(fastmath=True)
def _broken_powerlaw_plus_2peaks_unnormalized(
    m1,
    lam_0=0.361,
    lam_1=0.586,
    mpp_1=9.764,
    sigpp_1=0.649,
    mpp_2=32.763,
    sigpp_2=3.918,
    mlow_1=5.059,
    delta_m_1=4.321,
    break_mass=35.622,
    alpha_1=1.728,
    alpha_2=4.512,
    mmax=300.0,
):
    """
    Compute unnormalized mixture of broken power-law and two Gaussian peaks (equation B14).

    Evaluates the unnormalized primary-mass distribution prior to normalization,
    combining a broken power-law with two left-truncated Gaussian peaks. The mixture
    is tapered at low masses for smooth support enforcement.

    The mixture model is:
    $π(m_1 | Λ) ∝ [λ_0 p_{BP}(m_1) + λ_1 N_{lt}(m_1|μ_1,σ_1,m_{low,1}) +$
    $(1-λ_0-λ_1) N_{lt}(m_1|μ_2,σ_2,m_{low,1})] × S(m_1 | m_{low,1}, δm_1)$

    Parameters
    ----------
    m1 : ``numpy.ndarray``
        Primary masses (solar masses).
    lam_0, lam_1, mpp_1, sigpp_1, mpp_2, sigpp_2, mlow_1, delta_m_1, break_mass, alpha_1, alpha_2, mmax : ``float``
        Distribution parameters.
        +------------+---------+-----------------------------------------------------------+
        | param      | default | description                                               |
        +============+=========+===========================================================+
        | lam_0      | 0.361   | Fraction of broken power-law component (0 <= lam_0 <= 1)  |
        +------------+---------+-----------------------------------------------------------+
        | lam_1      | 0.586   | Fraction of first Gaussian peak (0 <= lam_1 <= 1)         |
        +------------+---------+-----------------------------------------------------------+
        | mpp_1      | 9.764   | Mean of first Gaussian peak (solar masses)                |
        +------------+---------+-----------------------------------------------------------+
        | sigpp_1    | 0.649   | Standard deviation of first Gaussian peak (solar masses)  |
        +------------+---------+-----------------------------------------------------------+
        | mpp_2      | 32.763  | Mean of second Gaussian peak (solar masses)               |
        +------------+---------+-----------------------------------------------------------+
        | sigpp_2    | 3.918   | Standard deviation of second Gaussian peak (solar masses) |
        +------------+---------+-----------------------------------------------------------+
        | mlow_1     | 5.059   | Minimum primary mass (solar masses)                       |
        +------------+---------+-----------------------------------------------------------+
        | delta_m_1  | 4.321   | Taper width (solar masses)                                |
        +------------+---------+-----------------------------------------------------------+
        | break_mass | 35.622  | Power-law break mass (solar masses)                       |
        +------------+---------+-----------------------------------------------------------+
        | alpha_1    | 1.728   | Power-law index below break mass                          |
        +------------+---------+-----------------------------------------------------------+
        | alpha_2    | 4.512   | Power-law index above break mass                          |
        +------------+---------+-----------------------------------------------------------+
        | mmax       | 300.0   | Maximum primary mass (solar masses)                       |
        +------------+---------+-----------------------------------------------------------+

    Returns
    -------
    unnormalized : ``numpy.ndarray``
        Unnormalized probability density at each mass value.
    """

    bp_term = lam_0 * _broken_powerlaw(
        m1,
        break_mass=break_mass,
        alpha_1=alpha_1,
        alpha_2=alpha_2,
        mlow_1=mlow_1,
        mmax=mmax,
    )
    peak1_term = lam_1 * truncated_normal_pdf(m1, mpp_1, sigpp_1, mlow_1)
    peak2_term = (1.0 - lam_0 - lam_1) * truncated_normal_pdf(
        m1, mpp_2, sigpp_2, mlow_1
    )

    mixture = bp_term + peak1_term + peak2_term

    return mixture * _smoothing_S(m1, mlow_1, delta_m_1)


@njit(fastmath=True)
def broken_powerlaw_plus_2peaks_pdf(
    m1,
    lam_0=0.361,
    lam_1=0.586,
    mpp_1=9.764,
    sigpp_1=0.649,
    mpp_2=32.763,
    sigpp_2=3.918,
    mlow_1=5.059,
    delta_m_1=4.321,
    break_mass=35.622,
    alpha_1=1.728,
    alpha_2=4.512,
    mmax=300.0,
    normalization_size=500,
):
    """
    Evaluate fully normalized primary-mass probability density (equation B20).

    Computes the normalized probability density $π(m_1 | Λ)$ by integrating
    the unnormalized mixture in log-space for numerical stability and consistency.

    Parameters
    ----------
    m1 : ``numpy.ndarray``
        Primary masses (solar masses).
    lam_0, lam_1, mpp_1, sigpp_1, mpp_2, sigpp_2, mlow_1, delta_m_1, break_mass, alpha_1, alpha_2, mmax : ``float``
        Distribution parameters.
        +------------+---------+-----------------------------------------------------------+
        | param      | default | description                                               |
        +============+=========+===========================================================+
        | lam_0      | 0.361   | Fraction of broken power-law component (0 <= lam_0 <= 1)  |
        +------------+---------+-----------------------------------------------------------+
        | lam_1      | 0.586   | Fraction of first Gaussian peak (0 <= lam_1 <= 1)         |
        +------------+---------+-----------------------------------------------------------+
        | mpp_1      | 9.764   | Mean of first Gaussian peak (solar masses)                |
        +------------+---------+-----------------------------------------------------------+
        | sigpp_1    | 0.649   | Standard deviation of first Gaussian peak (solar masses)  |
        +------------+---------+-----------------------------------------------------------+
        | mpp_2      | 32.763  | Mean of second Gaussian peak (solar masses)               |
        +------------+---------+-----------------------------------------------------------+
        | sigpp_2    | 3.918   | Standard deviation of second Gaussian peak (solar masses) |
        +------------+---------+-----------------------------------------------------------+
        | mlow_1     | 5.059   | Minimum primary mass (solar masses)                       |
        +------------+---------+-----------------------------------------------------------+
        | delta_m_1  | 4.321   | Taper width (solar masses)                                |
        +------------+---------+-----------------------------------------------------------+
        | break_mass | 35.622  | Power-law break mass (solar masses)                       |
        +------------+---------+-----------------------------------------------------------+
        | alpha_1    | 1.728   | Power-law index below break mass                          |
        +------------+---------+-----------------------------------------------------------+
        | alpha_2    | 4.512   | Power-law index above break mass                          |
        +------------+---------+-----------------------------------------------------------+
        | mmax       | 300.0   | Maximum primary mass (solar masses)                       |
        +------------+---------+-----------------------------------------------------------+
    normalization_size : ``int``
        Number of grid points for numerical integration.
        default: 500

    Returns
    -------
    pdf : ``numpy.ndarray``
        Normalized probability density at each mass value.

    Examples
    --------
    >>> from ler.gw_source_population.broken_powerlaw_plus_2peaks import broken_powerlaw_plus_2peaks_pdf
    >>> import numpy as np
    >>> m1 = np.array([10.0, 20.0, 50.0])
    >>> pdf_values = broken_powerlaw_plus_2peaks_pdf(m1)
    >>> print(pdf_values)  # doctest: +SKIP
    """

    # compute normalization constant for p(m1 | Λ) by integrating the unnormalized _broken_powerlaw_plus_2peaks_unnormalized over m1 from mlow_1 to mmax

    m1_grid = np.geomspace(mlow_1, mmax, normalization_size)
    pdf_unnormalized = _broken_powerlaw_plus_2peaks_unnormalized(m1_grid, lam_0, lam_1, mpp_1, sigpp_1, mpp_2, sigpp_2, mlow_1, delta_m_1, break_mass, alpha_1, alpha_2, mmax)
    _, _, m1_norm = cumulative_trapezoid(x=m1_grid, y=pdf_unnormalized, initial=0.0)

    return _broken_powerlaw_plus_2peaks_unnormalized(m1, lam_0, lam_1, mpp_1, sigpp_1, mpp_2, sigpp_2, mlow_1, delta_m_1, break_mass, alpha_1, alpha_2, mmax) / m1_norm


@njit(fastmath=True)
def broken_powerlaw_plus_2peaks_function(
    m1,
    lam_0=0.361,
    lam_1=0.586,
    mpp_1=9.764,
    sigpp_1=0.649,
    mpp_2=32.763,
    sigpp_2=3.918,
    mlow_1=5.059,
    delta_m_1=4.321,
    break_mass=35.622,
    alpha_1=1.728,
    alpha_2=4.512,
    mmax=300.0,
    normalization_size=500,
    R0=16.158, 
    kappa=3.166, 
    z_eval=0.2
):
    """
    Compute event rate as function of primary mass (equation B20).

    Combines the normalized primary-mass PDF with a redshift-dependent rate,
    modeling the merger rate evolution $R(z) = R_0 (1+z)^κ$.

    Parameters
    ----------
    m1 : ``numpy.ndarray``
        Primary masses (solar masses).
    lam_0, lam_1, mpp_1, sigpp_1, mpp_2, sigpp_2, mlow_1, delta_m_1, break_mass, alpha_1, alpha_2, mmax : ``float``
        Distribution parameters.
        +------------+---------+-----------------------------------------------------------+
        | param      | default | description                                               |
        +============+=========+===========================================================+
        | lam_0      | 0.361   | Fraction of broken power-law component (0 <= lam_0 <= 1)  |
        +------------+---------+-----------------------------------------------------------+
        | lam_1      | 0.586   | Fraction of first Gaussian peak (0 <= lam_1 <= 1)         |
        +------------+---------+-----------------------------------------------------------+
        | mpp_1      | 9.764   | Mean of first Gaussian peak (solar masses)                |
        +------------+---------+-----------------------------------------------------------+
        | sigpp_1    | 0.649   | Standard deviation of first Gaussian peak (solar masses)  |
        +------------+---------+-----------------------------------------------------------+
        | mpp_2      | 32.763  | Mean of second Gaussian peak (solar masses)               |
        +------------+---------+-----------------------------------------------------------+
        | sigpp_2    | 3.918   | Standard deviation of second Gaussian peak (solar masses) |
        +------------+---------+-----------------------------------------------------------+
        | mlow_1     | 5.059   | Minimum primary mass (solar masses)                       |
        +------------+---------+-----------------------------------------------------------+
        | delta_m_1  | 4.321   | Taper width (solar masses)                                |
        +------------+---------+-----------------------------------------------------------+
        | break_mass | 35.622  | Power-law break mass (solar masses)                       |
        +------------+---------+-----------------------------------------------------------+
        | alpha_1    | 1.728   | Power-law index below break mass                          |
        +------------+---------+-----------------------------------------------------------+
        | alpha_2    | 4.512   | Power-law index above break mass                          |
        +------------+---------+-----------------------------------------------------------+
        | mmax       | 300.0   | Maximum primary mass (solar masses)                       |
        +------------+---------+-----------------------------------------------------------+
    normalization_size : ``int``
        Number of grid points for numerical integration.
        default: 500
    R0 : ``float``
        Local merger rate (Gpc$^{-3}$ yr$^{-1}$).
        default: 16.158
    kappa : ``float``
        Power-law index for redshift evolution.
        default: 3.166
    z_eval : ``float``
        Redshift at which to evaluate the rate.
        default: 0.2

    Returns
    -------
    rate : ``numpy.ndarray``
        Differential merger rate $dR/dm_1$ at each mass value.

    Examples
    --------
    >>> from ler.gw_source_population.broken_powerlaw_plus_2peaks import broken_powerlaw_plus_2peaks_function
    >>> import numpy as np
    >>> m1 = np.array([10.0, 20.0, 50.0])
    >>> rate_values = broken_powerlaw_plus_2peaks_function(m1)
    >>> print(rate_values)  # doctest: +SKIP
    """

    # get pdf value
    pdf_m1 = broken_powerlaw_plus_2peaks_pdf(m1, lam_0, lam_1, mpp_1, sigpp_1, mpp_2, sigpp_2, mlow_1, delta_m_1, break_mass, alpha_1, alpha_2, mmax, normalization_size)

    # rate
    R_z = R0 * (1.0 + z_eval) ** kappa

    return R_z * pdf_m1


@njit(fastmath=True)
def broken_powerlaw_plus_2peaks_cdf(
    size,
    lam_0=0.361,
    lam_1=0.586,
    mpp_1=9.764,
    sigpp_1=0.649,
    mpp_2=32.763,
    sigpp_2=3.918,
    mlow_1=5.059,
    delta_m_1=4.321,
    break_mass=35.622,
    alpha_1=1.728,
    alpha_2=4.512,
    mmax=300.0,
):
    """
    Compute cumulative distribution function (CDF) for primary masses.

    Evaluates the CDF using log-space numerical integration for consistency
    with the PDF normalization, enabling inverse transform sampling.

    Parameters
    ----------
    size : ``int``
        Number of grid points for CDF computation.
    lam_0, lam_1, mpp_1, sigpp_1, mpp_2, sigpp_2, mlow_1, delta_m_1, break_mass, alpha_1, alpha_2, mmax : ``float``
        Distribution parameters.
        +------------+---------+-----------------------------------------------------------+
        | param      | default | description                                               |
        +============+=========+===========================================================+
        | lam_0      | 0.361   | Fraction of broken power-law component (0 <= lam_0 <= 1)  |
        +------------+---------+-----------------------------------------------------------+
        | lam_1      | 0.586   | Fraction of first Gaussian peak (0 <= lam_1 <= 1)         |
        +------------+---------+-----------------------------------------------------------+
        | mpp_1      | 9.764   | Mean of first Gaussian peak (solar masses)                |
        +------------+---------+-----------------------------------------------------------+
        | sigpp_1    | 0.649   | Standard deviation of first Gaussian peak (solar masses)  |
        +------------+---------+-----------------------------------------------------------+
        | mpp_2      | 32.763  | Mean of second Gaussian peak (solar masses)               |
        +------------+---------+-----------------------------------------------------------+
        | sigpp_2    | 3.918   | Standard deviation of second Gaussian peak (solar masses) |
        +------------+---------+-----------------------------------------------------------+
        | mlow_1     | 5.059   | Minimum primary mass (solar masses)                       |
        +------------+---------+-----------------------------------------------------------+
        | delta_m_1  | 4.321   | Taper width (solar masses)                                |
        +------------+---------+-----------------------------------------------------------+
        | break_mass | 35.622  | Power-law break mass (solar masses)                       |
        +------------+---------+-----------------------------------------------------------+
        | alpha_1    | 1.728   | Power-law index below break mass                          |
        +------------+---------+-----------------------------------------------------------+
        | alpha_2    | 4.512   | Power-law index above break mass                          |
        +------------+---------+-----------------------------------------------------------+
        | mmax       | 300.0   | Maximum primary mass (solar masses)                       |
        +------------+---------+-----------------------------------------------------------+

    Returns
    -------
    cdf_values : ``numpy.ndarray``
        Normalized CDF values at each grid point.
    m_grid : ``numpy.ndarray``
        Mass grid points (solar masses).

    Examples
    --------
    >>> from ler.gw_source_population.broken_powerlaw_plus_2peaks import broken_powerlaw_plus_2peaks_cdf
    >>> cdf_values, m_grid = broken_powerlaw_plus_2peaks_cdf(size=100)
    >>> print(cdf_values[-1])  # CDF should approach 1
    1.0
    """

    m_try = np.geomspace(mlow_1, mmax, size)
    pdf_unnormalized = _broken_powerlaw_plus_2peaks_unnormalized(m_try, lam_0, lam_1, mpp_1, sigpp_1, mpp_2, sigpp_2, mlow_1, delta_m_1, break_mass, alpha_1, alpha_2, mmax)

    cdf_values, m_try, _ = cumulative_trapezoid(
        x=m_try, y=pdf_unnormalized, initial=0.0
    )

    return cdf_values, m_try


@njit(fastmath=True)
def broken_powerlaw_plus_2peaks_rvs(
    size,
    lam_0=0.361,
    lam_1=0.586,
    mpp_1=9.764,
    sigpp_1=0.649,
    mpp_2=32.763,
    sigpp_2=3.918,
    mlow_1=5.059,
    delta_m_1=4.321,
    break_mass=35.622,
    alpha_1=1.728,
    alpha_2=4.512,
    mmax=300.0,
    normalization_size=500,
):
    """
    Generate random samples from the primary-mass distribution.

    Draws random samples using inverse transform sampling applied to the CDF,
    which is computed using numerically stable log-space integration.

    Parameters
    ----------
    size : ``int``
        Number of samples to generate.
    lam_0, lam_1, mpp_1, sigpp_1, mpp_2, sigpp_2, mlow_1, delta_m_1, break_mass, alpha_1, alpha_2, mmax : ``float``
        Distribution parameters.
        +------------+---------+-----------------------------------------------------------+
        | param      | default | description                                               |
        +============+=========+===========================================================+
        | lam_0      | 0.361   | Fraction of broken power-law component (0 <= lam_0 <= 1)  |
        +------------+---------+-----------------------------------------------------------+
        | lam_1      | 0.586   | Fraction of first Gaussian peak (0 <= lam_1 <= 1)         |
        +------------+---------+-----------------------------------------------------------+
        | mpp_1      | 9.764   | Mean of first Gaussian peak (solar masses)                |
        +------------+---------+-----------------------------------------------------------+
        | sigpp_1    | 0.649   | Standard deviation of first Gaussian peak (solar masses)  |
        +------------+---------+-----------------------------------------------------------+
        | mpp_2      | 32.763  | Mean of second Gaussian peak (solar masses)               |
        +------------+---------+-----------------------------------------------------------+
        | sigpp_2    | 3.918   | Standard deviation of second Gaussian peak (solar masses) |
        +------------+---------+-----------------------------------------------------------+
        | mlow_1     | 5.059   | Minimum primary mass (solar masses)                       |
        +------------+---------+-----------------------------------------------------------+
        | delta_m_1  | 4.321   | Taper width (solar masses)                                |
        +------------+---------+-----------------------------------------------------------+
        | break_mass | 35.622  | Power-law break mass (solar masses)                       |
        +------------+---------+-----------------------------------------------------------+
        | alpha_1    | 1.728   | Power-law index below break mass                          |
        +------------+---------+-----------------------------------------------------------+
        | alpha_2    | 4.512   | Power-law index above break mass                          |
        +------------+---------+-----------------------------------------------------------+
        | mmax       | 300.0   | Maximum primary mass (solar masses)                       |
        +------------+---------+-----------------------------------------------------------+
    normalization_size : ``int``
        Number of grid points for CDF computation.
        default: 500

    Returns
    -------
    samples : ``numpy.ndarray``
        Primary masses (solar masses) randomly sampled from the distribution.

    Examples
    --------
    >>> from ler.gw_source_population.broken_powerlaw_plus_2peaks import broken_powerlaw_plus_2peaks_rvs
    >>> m1_samples = broken_powerlaw_plus_2peaks_rvs(size=1000)
    >>> print(m1_samples.shape)
    (1000,)
    """

    cdf_values, x = broken_powerlaw_plus_2peaks_cdf(normalization_size, lam_0, lam_1, mpp_1, sigpp_1, mpp_2, sigpp_2, mlow_1, delta_m_1, break_mass, alpha_1, alpha_2, mmax)

    samples = inverse_transform_sampler(
        size=size,
        cdf=cdf_values,
        x=x,
    )

    return samples

# ---------------------------------
# Mass Ratio Distribution
# ---------------------------------
@njit(fastmath=True)
def mass_ratio_powerlaw_with_smoothing_scalar(
    q, 
    m1, 
    beta=1.171, 
    mlow_2=3.551, 
    delta_m_2=4.910, 
    normalization_size=500
):
    """
    Compute normalized conditional mass ratio density for a scalar secondary mass.

    Evaluates the normalized PDF $p(q | m_1, Λ)$ for a single primary mass value,
    integrating the unnormalized density to determine the normalization constant.

    Parameters
    ----------
    q : ``float``
        Mass ratio ($m_2 / m_1$), bounded 0 < q ≤ 1.
    m1 : ``float``
        Primary mass (solar masses).
    beta : ``float``
        Power-law index for mass ratio distribution.
        default: 1.171
    mlow_2 : ``float``
        Minimum secondary mass (solar masses).
        default: 3.551
    delta_m_2 : ``float``
        Taper width for secondary mass (solar masses).
        default: 4.910
    normalization_size : ``int``
        Number of grid points for numerical integration.
        default: 500

    Returns
    -------
    pdf_value : ``float``
        Normalized probability density at the given (q, m1) pair.
    """

    # compute normalization constant for p(q | m1, Λ) by integrating over q from q_min to 1
    q_min = mlow_2 / m1
    q_grid = np.linspace(q_min, 1.0, normalization_size)
    q_pdf_unnormalized = powerlaw_with_smoothing(
        q_grid, q_grid * m1, mlow_2, beta, delta_m_2
    )

    pdf = (q**beta) * _smoothing_S(np.array([m1 * q]), mlow_2, delta_m_2)[0]
    _, _, norm = cumulative_trapezoid(
        x=q_grid, y=q_pdf_unnormalized, initial=0.0
    )  # total integral for normalization

    return pdf / norm


@njit(fastmath=True)
def mass_ratio_powerlaw_with_smoothing_pdf(
    q, 
    m1, 
    beta=1.171, 
    mlow_2=3.551, 
    delta_m_2=4.910, 
    normalization_size=500
):
    """
    Evaluate fully normalized conditional mass ratio density (equation B22).

    Computes the normalized conditional PDF $p(q | m_1, Λ)$ for mass ratio distributions.
    Supports both scalar and array inputs through element-wise evaluation.

    Parameters
    ----------
    q : ``numpy.ndarray``
        Mass ratios ($m_2 / m_1$), bounded 0 < q ≤ 1. Can be scalar or array.
    m1 : ``numpy.ndarray``
        Primary masses (solar masses). Must have same shape as q.
    beta : ``float``
        Power-law index for mass ratio distribution.
        default: 1.171
    mlow_2 : ``float``
        Minimum secondary mass (solar masses).
        default: 3.551
    delta_m_2 : ``float``
        Taper width for secondary mass (solar masses).
        default: 4.910
    normalization_size : ``int``
        Number of grid points for numerical integration.
        default: 500

    Returns
    -------
    pdf : ``numpy.ndarray`` or ``float``
        Normalized probability density at each (q, m1) pair. Returns scalar if input is scalar.

    Examples
    --------
    >>> from ler.gw_source_population.broken_powerlaw_plus_2peaks import mass_ratio_powerlaw_with_smoothing_pdf
    >>> import numpy as np
    >>> q = np.array([0.5, 0.7, 0.9])
    >>> m1 = np.array([20.0, 30.0, 40.0])
    >>> pdf_values = mass_ratio_powerlaw_with_smoothing_pdf(q, m1)
    >>> print(pdf_values)  # doctest: +SKIP
    """
    size = q.shape[0]
    result = np.zeros(size, dtype=float)

    for i in prange(size):
        result[i] = mass_ratio_powerlaw_with_smoothing_scalar(q[i], m1[i], beta, mlow_2, delta_m_2, normalization_size)

    return result


@njit(fastmath=True)
def mass_ratio_powerlaw_with_smoothing_cdf(size, m1, beta=1.171, mlow_2=3.551, delta_m_2=4.910):
    """
    Compute cumulative distribution function for mass ratio distribution.

    Evaluates the CDF $P(q ≤ q_0 | m_1, Λ)$ using spline-based numerical integration,
    enabling efficient inverse transform sampling for a single primary mass value.

    Parameters
    ----------
    size : ``int``
        Number of grid points for CDF computation.
    m1 : ``float``
        Primary mass (solar masses).
    beta : ``float``
        Power-law index for mass ratio distribution.
        default: 1.171
    mlow_2 : ``float``
        Minimum secondary mass (solar masses).
        default: 3.551
    delta_m_2 : ``float``
        Taper width for secondary mass (solar masses).
        default: 4.910

    Returns
    -------
    cdf_values : ``numpy.ndarray``
        Normalized CDF values at each grid point.
    q_grid : ``numpy.ndarray``
        Mass ratio grid points (0 < q ≤ 1).

    Examples
    --------
    >>> from ler.gw_source_population.broken_powerlaw_plus_2peaks import mass_ratio_powerlaw_with_smoothing_cdf
    >>> cdf_values, q_grid = mass_ratio_powerlaw_with_smoothing_cdf(size=100, m1=30.0, beta=1.171, mlow_2=3.551, delta_m_2=4.910)
    >>> print(cdf_values[-1])  # CDF should approach 1
    1.0
    """
    # beta, mlow_2, delta_m_2 = params

    q_try = np.linspace(mlow_2 / m1, 1.0, size)
    pdf_unnormalized = powerlaw_with_smoothing(
        q_try, q_try * m1, mlow_2, beta, delta_m_2
    )

    cdf_values, q_try, _ = cumulative_trapezoid(
        x=q_try, y=pdf_unnormalized, initial=0.0
    )

    return cdf_values, q_try


@njit(parallel=True, fastmath=True)
def mass_ratio_powerlaw_with_smoothing_rvs(m1, beta=1.171, mlow_2=3.551, delta_m_2=4.910, normalization_size=500):
    """
    Generate random samples from the mass ratio distribution.

    Draws random mass ratios using inverse transform sampling applied to the CDF
    for each primary mass value independently. This function uses parallel processing
    to accelerate sampling for large arrays.

    Parameters
    ----------
    m1 : ``numpy.ndarray``
        Primary masses (solar masses).
    beta : ``float``
        Power-law index for mass ratio distribution.
        default: 1.171
    mlow_2 : ``float``
        Minimum secondary mass (solar masses).
        default: 3.551
    delta_m_2 : ``float``
        Taper width for secondary mass (solar masses).
        default: 4.910
    normalization_size : ``int``
        Number of grid points for CDF computation.
        default: 500

    Returns
    -------
    q_samples : ``numpy.ndarray``
        Mass ratio samples ($m_2 / m_1$), bounded 0 < q ≤ 1.

    Examples
    --------
    >>> from ler.gw_source_population.broken_powerlaw_plus_2peaks import mass_ratio_powerlaw_with_smoothing_rvs
    >>> import numpy as np
    >>> m1 = np.array([20.0, 30.0, 40.0])
    >>> q_samples = mass_ratio_powerlaw_with_smoothing_rvs(m1)
    >>> print(q_samples.shape)
    (3,)
    """
    size = m1.shape[0]
    samples = np.empty(size, dtype=float)

    for i in prange(size):
        # Guard: if m1 is too small for a valid q range, return q=1
        q_min = mlow_2 / m1[i]
        if q_min >= 1.0:
            samples[i] = 1.0
            continue

        cdf_values, x = mass_ratio_powerlaw_with_smoothing_cdf(normalization_size, m1[i], beta, mlow_2, delta_m_2)

        samples[i] = inverse_transform_sampler(
            size=1,
            cdf=cdf_values,
            x=x,
        )[0]

    return samples
# ---------------------------------