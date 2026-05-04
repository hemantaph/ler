# -*- coding: utf-8 -*-
"""
Module for JIT-compiled functions used in GW source population simulations.

This module provides Numba JIT-compiled functions for efficient sampling and
computation of merger rate densities, star formation rates, and mass distributions for various gravitational wave source populations including BBH, BNS, and NSBH.

Key Features: \n
- Merger rate density functions for PopI/II, PopIII, and Primordial BBH \n
- Star formation rate models (Madau & Dickinson 2014, Madau & Fragos 2017) \n
- Mass distribution samplers (power-law Gaussian, broken power-law, bimodal) \n
- JIT-compiled for high performance with Numba \n

Copyright (C) 2026 Hemanta Ph. Distributed under MIT License.
"""

import numpy as np
import math
from numba import njit
from scipy.interpolate import CubicSpline
from astropy.cosmology import LambdaCDM

from ..utils import inverse_transform_sampler, cumulative_trapezoid, is_njitted

# ------------------------------
# Merger rate density functions
# ------------------------------
@njit(fastmath=True)
def merger_rate_density_bbh_oguri2018_function(zs, R0=19 * 1e-9, b2=1.6, b3=2.1, b4=30):
    """
    Compute the merger rate density for PopI/II BBH.

    Reference: Oguri et al. (2018). The source-frame rate density is

    .. math::

        R(z) = R_0 \\frac{(b_4 + 1) e^{b_2 z}}{b_4 + e^{b_3 z}}.

    Parameters
    ----------
    zs : ``float`` or ``numpy.ndarray``
        Source redshifts.
    R0 : ``float``
        Local merger rate density at low redshift (Mpc^-3 yr^-1). \n
        default: 19e-9 (GWTC-4)
    b2 : ``float``
        Fitting parameter. \n
        default: 1.6
    b3 : ``float``
        Fitting parameter. \n
        default: 2.1
    b4 : ``float``
        Fitting parameter. \n
        default: 30

    Returns
    -------
    rate_density : ``float`` or ``numpy.ndarray``
        Merger rate density.

    Examples
    --------
    >>> import numpy as np
    >>> from ler.gw_source_population import merger_rate_density_bbh_oguri2018_function
    >>> rate_density = merger_rate_density_bbh_oguri2018_function(zs=np.array([0.1]))
    """
    return R0 * (b4 + 1) * np.exp(b2 * zs) / (b4 + np.exp(b3 * zs))

@njit(fastmath=True)
def merger_rate_density_bbh_popIII_ken2022_function(zs, R0=19.2 * 1e-9, aIII=0.66, bIII=0.3, zIII=11.6):
    """
    Compute the unnormalized merger rate density for PopIII BBH.

    Reference: Ng et al. (2022). The model is

    .. math::

        R(z) = R_0 \\frac{e^{a_{\\rm III}(z-z_{\\rm III})}}
        {b_{\\rm III} + a_{\\rm III} e^{(a_{\\rm III}+b_{\\rm III})(z-z_{\\rm III})}}.

    Parameters
    ----------
    zs : ``float`` or ``numpy.ndarray``
        Source redshifts.
    R0 : ``float``
        Normalization constant. \n
        default: 19.2e-9
    aIII : ``float``
        Fitting parameter. \n
        default: 0.66
    bIII : ``float``
        Fitting parameter. \n
        default: 0.3
    zIII : ``float``
        Characteristic redshift. \n
        default: 11.6

    Returns
    -------
    rate_density : ``float`` or ``numpy.ndarray``
        Merger rate density.

    Examples
    --------
    >>> import numpy as np
    >>> from ler.gw_source_population import merger_rate_density_bbh_popIII_ken2022_function
    >>> rate_density = merger_rate_density_bbh_popIII_ken2022_function(zs=np.array([0.1]))
    """
    return (
        R0
        * np.exp(aIII * (zs - zIII))
        / (bIII + aIII * np.exp((aIII + bIII) * (zs - zIII)))
    )

@njit(fastmath=True)
def merger_rate_density_madau_dickinson2014_function(zs, R0=19 * 1e-9, a=0.015, b=2.7, c=2.9, d=5.6):
    """
    Compute the merger rate density for BBH using Madau & Dickinson (2014) model.

    The shape follows

    .. math::

        \\psi(z) = a \\frac{(1+z)^b}{1 + \\left((1+z)/c\\right)^d},
        \\qquad
        R(z) = R_0 \\frac{\\psi(z)}{\\psi(0)}.

    Reference: Eq. 15 of https://arxiv.org/pdf/1403.0007

    Parameters
    ----------
    zs : ``float`` or ``numpy.ndarray``
        Source redshifts.
    R0 : ``float``
        Local merger rate density (Mpc^-3 yr^-1). \n
        default: 19e-9
    a : ``float``
        Normalization parameter. \n
        default: 0.015
    b : ``float``
        Low-redshift power-law slope. \n
        default: 2.7
    c : ``float``
        Turnover redshift parameter. \n
        default: 2.9
    d : ``float``
        High-redshift power-law slope. \n
        default: 5.6

    Returns
    -------
    rate_density : ``float`` or ``numpy.ndarray``
        Merger rate density (Mpc^-3 yr^-1).

    Examples
    --------
    >>> import numpy as np
    >>> from ler.gw_source_population import merger_rate_density_madau_dickinson2014_function
    >>> rate_density = merger_rate_density_madau_dickinson2014_function(zs=np.array([0.1]))
    """

    def density_helper(zs):
        return sfr_madau_dickinson2014(
            zs=zs,
            a=a,
            b=b,
            c=c,
            d=d,
        )

    density_zs = R0 * density_helper(zs)/ density_helper(np.array([0.]))[0]

    return density_zs

@njit(fastmath=True)
def merger_rate_density_madau_dickinson_belczynski_ng_function(zs, R0=19 * 1e-9, alpha_F=2.57, beta_F=5.83, c_F=3.36):
    """
    Compute BBH merger rate density following Ng et al. (2021).

    This model uses a Madau-Dickinson-like functional form to fit the 
    merger rate density of field BHs, accounting for time delays and 
    metallicity effects. Coefficients from Madau & Dickinson (2014) are translated as: B-> alpha_F, D-> beta_F, C-> c_F.

    .. math::

        \\psi_F(z) = \\frac{(1+z)^{\\alpha_F}}
        {1 + \\left((1+z)/c_F\\right)^{\\beta_F}},
        \\qquad
        R(z) = R_0 \\frac{\\psi_F(z)}{\\psi_F(0)}.

    Parameters
    ----------
    zs : ``float`` or ``numpy.ndarray``
        Source redshifts.
    R0 : ``float``
        Local merger rate density (Mpc^-3 yr^-1). \n
        default: 19e-9
    alpha_F : ``float``
        Low-redshift power-law slope. \n
        default: 2.57
    beta_F : ``float``
        High-redshift power-law slope. \n
        default: 5.83
    c_F : ``float``
        Turnover redshift parameter. \n
        default: 3.36

    Returns
    -------
    rate_density : ``float`` or ``numpy.ndarray``
        Merger rate density (Mpc^-3 yr^-1).

    Examples
    --------
    >>> import numpy as np
    >>> from ler.gw_source_population import merger_rate_density_madau_dickinson_belczynski_ng_function
    >>> rate_density = merger_rate_density_madau_dickinson_belczynski_ng_function(zs=np.array([0.1]))
    """

    def density_helper(zs):
        return sfr_madau_dickinson2014(
            zs=zs,
            a=1.0,
            b=alpha_F,
            c=c_F,
            d=beta_F,
        )
    density_zs = R0 * density_helper(zs)/ density_helper(np.array([0.]))[0]

    return density_zs

def merger_rate_density_bbh_primordial_ken2022_function(zs, cosmology=None, R0=0.044 * 1e-9, t0=13.786885302009708):
    """
    Compute the merger rate density for Primordial BBH.

    Reference: Ng et al. (2022). The rate is parameterized by cosmic age:

    .. math::

        R(z) = R_0 \\left(\\frac{t(z)}{t_0}\\right)^{-34/37}.

    Parameters
    ----------
    zs : ``float`` or ``numpy.ndarray``
        Source redshifts.
    cosmology : ``astropy.cosmology`` or ``None``
        Cosmology object for age calculations. \n
        default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)
    R0 : ``float``
        Normalization constant. \n
        default: 0.044e-9
    t0 : ``float``
        Present age of the Universe (Gyr). \n
        default: 13.786885302009708

    Returns
    -------
    rate_density : ``float`` or ``numpy.ndarray``
        Merger rate density.

    Examples
    --------
    >>> import numpy as np
    >>> from ler.gw_source_population import merger_rate_density_bbh_primordial_ken2022_function
    >>> rate_density = merger_rate_density_bbh_primordial_ken2022_function(zs=np.array([0.1]))
    """
    if cosmology is None:
        cosmology = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)

    rate_density = R0 * (cosmology.age(z=zs).value / t0) ** (-34 / 37)
    return rate_density


def sfr_madau_fragos2017_with_bbh_td(zs, R0=19 * 1e-9):
    """
    Compute the merger rate density for BBH. This is computed from star formation rate, Madau & Fragos (2017), with an additional time delay. This function relies on pre-generated data points. 

    Parameters
    ----------
    zs : ``float`` or ``numpy.ndarray``
        Source redshifts.
    R0 : ``float``
        Local merger rate density (Mpc^-3 yr^-1). \n
        default: 19e-9

    Returns
    -------
    SFR : ``float`` or ``numpy.ndarray``
        Star formation rate (Mpc^-3 yr^-1).
    """
    rm = np.array([1.00304765, 1.00370075, 1.00449545, 1.00546251, 1.00663937, 1.00807168, 1.00981505, 1.01193727, 1.01046483, 1.01359803, 1.01741386, 1.02206193, 1.02772495, 1.03462601, 1.04303746, 1.05329142, 1.07093106, 1.08624215, 1.10489848, 1.12760683,1.15519183, 1.18858451, 1.22878158, 1.27676494, 1.33727882, 1.40335222, 1.47956936, 1.56759515, 1.6711375 , 1.79690371, 1.95410462, 2.15201042, 2.36151109, 2.66742932, 3.04354598, 3.49048755, 3.98122536, 4.42347511, 4.61710896, 4.30190679, 3.50890876, 2.37699066, 1.41830834, 0.77944771,0.40667706, 0.20463758, 0.09975143, 0.04745116])
    zs_ = np.geomspace(0.001, 10, 48)

    spline = CubicSpline(zs_, rm, extrapolate=True)
    SFR = spline(zs)*R0
    return SFR

def sfr_madau_dickinson2014_with_bbh_td(zs, R0=19 * 1e-9):
    """
    Compute the merger rate density for BBH. This is computed from star formation rate, Madau & Dickinson (2014), with an additional time delay. This function relies on pre-generated data points. 

    Parameters
    ----------
    zs : ``float`` or ``numpy.ndarray``
        Source redshifts.
    R0 : ``float``
        Local merger rate density (Mpc^-3 yr^-1). \n
        default: 19e-9

    Returns
    -------
    SFR : ``float`` or ``numpy.ndarray``
        Star formation rate (Mpc^-3 yr^-1).
    """
    rm = np.array([1.00292325, 1.0035494 , 1.0043112 , 1.00523807, 1.0063658 , 1.00773798, 1.00940767, 1.01143948, 1.00997839, 1.01297699, 1.01662662, 1.02106895, 1.02647649, 1.03305927, 1.04107277, 1.05082743, 1.06802831, 1.08255749, 1.10022606, 1.12169013, 1.14772134, 1.1792097 , 1.21715406, 1.26263694, 1.32051095, 1.38462461, 1.45997648, 1.54851567, 1.65349288, 1.78046645, 1.93811129, 2.1354612 , 2.34086287, 2.63664802, 2.98892341, 3.38353439, 3.76990612, 4.03489696, 4.00806904, 3.56766897, 2.86966689, 2.01282062, 1.29696347, 0.78913584, 0.46166281, 0.26226345, 0.14509118, 0.07854392])
    zs_ = np.geomspace(0.001, 10, 48)

    spline = CubicSpline(zs_, rm, extrapolate=True)
    SFR = spline(zs)*R0
    return SFR

def sfr_madau_fragos2017_with_bns_td(zs, R0=89 * 1e-9):
    """
    Compute the merger rate density for BNS. This is computed from star formation rate, Madau & Fragos (2017), with an additional time delay. This function relies on pre-generated data points. 

    Parameters
    ----------
    zs : ``float`` or ``numpy.ndarray``
        Source redshifts.
    R0 : ``float``
        Local merger rate density (Mpc^-3 yr^-1). \n
        default: 89e-9

    Returns
    -------
    SFR : ``float`` or ``numpy.ndarray``
        Star formation rate (Mpc^-3 yr^-1).
    """
    rm = np.array([1.00309364, 1.00375139, 1.00455175, 1.00552568, 1.00671091, 1.00815339, 1.00990912, 1.01204635, 1.00757017, 1.01071962, 1.01455507, 1.01922677, 1.02491815, 1.03185311, 1.04030479, 1.05060602, 1.06970166, 1.08508957, 1.10382838, 1.12661829, 1.15427005, 1.18768774, 1.22781836, 1.27555711, 1.31791484, 1.38209039, 1.4555543 , 1.5397332 , 1.63806934, 1.75685668, 1.90448546, 2.08862044, 2.34440211, 2.63899295, 2.99729389, 3.41567274, 3.86324106, 4.24545603, 4.37018218, 4.00555831, 3.10525751, 2.06354992, 1.20906304, 0.65233811, 0.33356891, 0.16397688, 0.08024945, 0.036953])
    zs_ = np.geomspace(0.001, 10, 48)

    spline = CubicSpline(zs_, rm, extrapolate=True)
    SFR = spline(zs)*R0
    return SFR

def sfr_madau_dickinson2014_with_bns_td(zs, R0=89 * 1e-9):
    """
    Compute the merger rate density for BNS. This is computed from star formation rate, Madau & Dickinson (2014), with an additional time delay. This function relies on pre-generated data points. 

    Parameters
    ----------
    zs : ``float`` or ``numpy.ndarray``
        Source redshifts.
    R0 : ``float``
        Local merger rate density (Mpc^-3 yr^-1). \n
        default: 89e-9

    Returns
    -------
    SFR : ``float`` or ``numpy.ndarray``
        Star formation rate (Mpc^-3 yr^-1).
    """
    rm = np.array([1.0029945 , 1.00362259, 1.00438674, 1.00531645, 1.00644763, 1.00782396, 1.00949865, 1.01153645, 1.00240992, 1.00539605, 1.00903013, 1.01345293, 1.01883579, 1.02538714, 1.03336026, 1.04306247, 1.05841698, 1.07283625, 1.09035966, 1.11162909, 1.1373949 , 1.16851509, 1.20594024, 1.25068092, 1.3085267 , 1.37111306, 1.44421094, 1.52948237, 1.62985636, 1.75058453, 1.90010572, 2.0870216 , 2.33573104, 2.6218286 , 2.96031682, 3.3343522 , 3.69149889, 3.92099769, 3.86227814, 3.40811745, 2.59314381, 1.79588097, 1.14260538, 0.686002  , 0.3954134 , 0.22083291, 0.11548455, 0.06064368])
    zs_ = np.geomspace(0.001, 10, 48)

    spline = CubicSpline(zs_, rm, extrapolate=True)
    SFR = spline(zs)*R0
    return SFR

# ------------------------------
# Star formation rate functions
# ------------------------------
@njit(fastmath=True)
def sfr_madau_fragos2017(zs, a=0.01, b=2.6, c=3.2, d=6.2):
    """
    Compute star formation rate using Madau & Fragos (2017) model.

    Reference: https://arxiv.org/pdf/1606.07887.pdf

    Parameters
    ----------
    zs : ``float`` or ``numpy.ndarray``
        Source redshifts.
    a : ``float``
        Normalization parameter. \n
        default: 0.01
    b : ``float``
        Low-redshift power-law slope. \n
        default: 2.6
    c : ``float``
        Turnover redshift parameter. \n
        default: 3.2
    d : ``float``
        High-redshift power-law slope. \n
        default: 6.2

    Returns
    -------
    SFR : ``float`` or ``numpy.ndarray``
        Star formation rate (Msun yr^-1 Mpc^-3).
    """
    return a * (1+zs)**b / (1 + ((1+zs)/c)**d)

@njit(fastmath=True)
def sfr_madau_dickinson2014(zs, a=0.015, b=2.7, c=2.9, d=5.6):
    """
    Compute star formation rate using Madau & Dickinson (2014) model.

    Reference: Eq. 15 of https://arxiv.org/pdf/1403.0007

    Parameters
    ----------
    zs : ``float`` or ``numpy.ndarray``
        Source redshifts.
    a : ``float``
        Normalization parameter. \n
        default: 0.015
    b : ``float``
        Low-redshift power-law slope. \n
        default: 2.7
    c : ``float``
        Turnover redshift parameter. \n
        default: 2.9
    d : ``float``
        High-redshift power-law slope. \n
        default: 5.6

    Returns
    -------
    SFR : ``float`` or ``numpy.ndarray``
        Star formation rate (Msun yr^-1 Mpc^-3).

    Examples
    --------
    >>> from ler.gw_source_population import sfr_madau_dickinson2014
    >>> sfr = sfr_madau_dickinson2014(zs=np.array([0.1]))
    """
    return a * (1 + zs) ** b / (1 + ((1 + zs) / c) ** d)

# ------------------------------
# Binary mass functions
# ------------------------------
@njit(fastmath=True)
def _lognormal_psi(m, Mc, sigma):
    """
    Compute lognormal mass function (equation 1, Ng et al. 2022).

    Parameters
    ----------
    m : ``numpy.ndarray``
        Mass values (solar masses).
    Mc : ``float``
        Characteristic mass scale (solar masses).
    sigma : ``float``
        Width parameter of the lognormal distribution.

    Returns
    -------
    psi : ``numpy.ndarray``
        Mass function values.
    """
    return np.exp(-np.log(m / Mc) ** 2 / (2 * sigma**2)) / (
        np.sqrt(2 * np.pi) * sigma * m
    )

@njit(fastmath=True)
def ng2022_lognormal_joint_pdf(m1, m2, Mc, sigma):
    """
    Compute joint probability density for lognormal 2D mass distribution.

    Parameters
    ----------
    m1 : ``numpy.ndarray``
        Primary masses (solar masses).
    m2 : ``numpy.ndarray``
        Secondary masses (solar masses).
    Mc : ``float``
        Characteristic mass scale (solar masses).
    sigma : ``float``
        Width parameter of the lognormal distribution.

    Returns
    -------
    pdf : ``numpy.ndarray``
        Joint probability density values.
    """
    return (
        (m1 + m2) ** (36 / 37)
        * (m1 * m2) ** (32 / 37)
        * _lognormal_psi(m1, Mc, sigma)
        * _lognormal_psi(m2, Mc, sigma)
    )

@njit(fastmath=True)
def binary_masses_BBH_popIII_lognormal_rvs(size, m_min=1.0, m_max=100.0, Mc=20.0, sigma=0.3, chunk_size=10000):
    """
    Generate random samples of binary masses for PopIII BBH (lognormal model).

    Draws samples from a 2D lognormal distribution in mass space using rejection
    sampling. Reference: Ng et al. (2022).

    Parameters
    ----------
    size : ``int``
        Number of binary systems to sample.
    m_min : ``float``
        Minimum mass (solar masses). 
        default: 1.0
    m_max : ``float``
        Maximum mass (solar masses). 
        default: 100.0
    Mc : ``float``
        Characteristic mass scale (solar masses). 
        default: 20.0
    sigma : ``float``
        Width parameter of the lognormal distribution. 
        default: 0.3
    chunk_size : ``int``
        Number of samples per rejection sampling chunk. 
        default: 10000

    Returns
    -------
    m1 : ``numpy.ndarray``
        Primary mass samples (solar masses).
    m2 : ``numpy.ndarray``
        Secondary mass samples (solar masses).

    Examples
    --------
    >>> from ler.gw_source_population.prior_functions import binary_masses_BBH_popIII_lognormal_rvs
    >>> m1, m2 = binary_masses_BBH_popIII_lognormal_rvs(size=1000)
    >>> print(m1.shape, m2.shape)
    (1000,) (1000,)
    """
    # rejection sampling initialization
    m1 = np.random.uniform(m_min, m_max, chunk_size)
    m2 = np.random.uniform(m_min, m_max, chunk_size)
    z = ng2022_lognormal_joint_pdf(m1, m2, Mc, sigma)
    zmax = np.max(z)

    # rejection sampling in chunks
    m1_sample = np.zeros(size)
    m2_sample = np.zeros(size)
    old_num = 0
    while True:
        m1_try = np.random.uniform(m_min, m_max, size=chunk_size)
        m2_try = np.random.uniform(m_min, m_max, size=chunk_size)

        z_try = np.random.uniform(0, zmax, size=chunk_size)
        zmax = max(zmax, np.max(z_try))
        idx = z_try < ng2022_lognormal_joint_pdf(m1_try, m2_try, Mc, sigma)
        new_num = old_num + np.sum(idx)
        if new_num >= size:
            m1_sample[old_num:size] = m1_try[idx][: size - old_num]
            m2_sample[old_num:size] = m2_try[idx][: size - old_num]
            break
        else:
            m1_sample[old_num:new_num] = m1_try[idx]
            m2_sample[old_num:new_num] = m2_try[idx]
            old_num = new_num

    # swap masses to ensure m1 >= m2
    idx = m1_sample < m2_sample
    m1_sample[idx], m2_sample[idx] = m2_sample[idx], m1_sample[idx]
    return m1_sample, m2_sample

@njit(fastmath=True)
def binary_masses_BBH_primordial_lognormal_rvs(size, m_min=1.0, m_max=100.0, Mc=20.0, sigma=0.3, chunk_size=10000):
    """
    Generate random samples of binary masses for Primordial BBH (lognormal model).

    Draws samples from a 2D lognormal distribution in mass space using rejection
    sampling. Reference: Ng et al. (2022).

    Parameters
    ----------
    size : ``int``
        Number of binary systems to sample.
    m_min : ``float``
        Minimum mass (solar masses). 
        default: 1.0
    m_max : ``float``
        Maximum mass (solar masses). 
        default: 100.0
    Mc : ``float``
        Characteristic mass scale (solar masses). 
        default: 20.0
    sigma : ``float``
        Width parameter of the lognormal distribution. 
        default: 0.3
    chunk_size : ``int``
        Number of samples per rejection sampling chunk. 
        default: 10000

    Returns
    -------
    m1 : ``numpy.ndarray``
        Primary mass samples (solar masses).
    m2 : ``numpy.ndarray``
        Secondary mass samples (solar masses).

    Examples
    --------
    >>> from ler.gw_source_population.prior_functions import binary_masses_BBH_primordial_lognormal_rvs
    >>> m1, m2 = binary_masses_BBH_primordial_lognormal_rvs(size=1000)
    >>> print(m1.shape, m2.shape)
    (1000,) (1000,)
    """
    # rejection sampling initialization
    m1 = np.random.uniform(m_min, m_max, chunk_size)
    m2 = np.random.uniform(m_min, m_max, chunk_size)
    z = ng2022_lognormal_joint_pdf(m1, m2, Mc, sigma)
    zmax = np.max(z)

    # rejection sampling in chunks
    m1_sample = np.zeros(size)
    m2_sample = np.zeros(size)
    old_num = 0
    while True:
        m1_try = np.random.uniform(m_min, m_max, size=chunk_size)
        m2_try = np.random.uniform(m_min, m_max, size=chunk_size)

        z_try = np.random.uniform(0, zmax, size=chunk_size)
        zmax = max(zmax, np.max(z_try))
        idx = z_try < ng2022_lognormal_joint_pdf(m1_try, m2_try, Mc, sigma)
        new_num = old_num + np.sum(idx)
        if new_num >= size:
            m1_sample[old_num:size] = m1_try[idx][: size - old_num]
            m2_sample[old_num:size] = m2_try[idx][: size - old_num]
            break
        else:
            m1_sample[old_num:new_num] = m1_try[idx]
            m2_sample[old_num:new_num] = m2_try[idx]
            old_num = new_num

    # swap masses to ensure m1 >= m2
    idx = m1_sample < m2_sample
    m1_sample[idx], m2_sample[idx] = m2_sample[idx], m1_sample[idx]
    return m1_sample, m2_sample

@njit(fastmath=True)
def _erf(x):
    """
    Compute the error function using Abramowitz & Stegun approximation.

    Parameters
    ----------
    x : ``float`` or ``numpy.ndarray``
        Input value(s).

    Returns
    -------
    result : ``float`` or ``numpy.ndarray``
        Error function value(s).
    """
    # constants for A&S formula 7.1.26
    p = 0.3275911
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429

    sign = np.sign(x)
    x = np.abs(x)

    t = 1.0 / (1.0 + p * x)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * np.exp(-x * x)

    return sign * y

@njit(fastmath=True)
def _compute_normalization_factor(mu, sigma, mmin, mmax):
    """
    Compute normalization factor for truncated Gaussian.

    Parameters
    ----------
    mu : ``float``
        Mean of the Gaussian.
    sigma : ``float``
        Standard deviation of the Gaussian.
    mmin : ``float``
        Minimum mass bound.
    mmax : ``float``
        Maximum mass bound.

    Returns
    -------
    N : ``float``
        Normalization factor.
    """
    part1 = (mmax - mu) / (np.sqrt(2) * sigma)
    part2 = (mmin - mu) / (np.sqrt(2) * sigma)
    N = np.sqrt(2 * np.pi) * sigma * (0.5 * (_erf(part1) - _erf(part2)))
    return N

@njit(fastmath=True)
def _bimodal_unnormalized(m, w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1.0, mmax=2.3):
    """
    Compute unnormalized bimodal Gaussian density for BNS mass distribution.

    Parameters
    ----------
    m : ``numpy.ndarray``
        Mass values (solar masses).
    w : ``float``
        Weight of the left (low-mass) peak. 
        default: 0.643
    muL : ``float``
        Mean of the left peak (solar masses). 
        default: 1.352
    sigmaL : ``float``
        Standard deviation of the left peak (solar masses). 
        default: 0.08
    muR : ``float``
        Mean of the right peak (solar masses). 
        default: 1.88
    sigmaR : ``float``
        Standard deviation of the right peak (solar masses). 
        default: 0.3
    mmin : ``float``
        Minimum mass (solar masses). 
        default: 1.0
    mmax : ``float``
        Maximum mass (solar masses). 
        default: 2.3

    Returns
    -------
    density : ``numpy.ndarray``
        Unnormalized probability density values.
    """
    # left peak (before truncation)
    pdf_unnormL = np.exp(-((m - muL) ** 2) / (2 * sigmaL**2))
    # right peak (before truncation)
    pdf_unnormR = np.exp(-((m - muR) ** 2) / (2 * sigmaR**2))
    # mixture before normalization
    density = w * pdf_unnormL + (1 - w) * pdf_unnormR
    return density
    
@njit(fastmath=True)
def bimodal_pdf(m, w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1.0, mmax=2.3):
    """
    Evaluate fully normalized bimodal Gaussian PDF for BNS mass distribution.

    Computes the normalized probability density function combining two truncated
    Gaussian peaks. Reference: Will M. Farr et al. 2020.

    Parameters
    ----------
    m : ``numpy.ndarray``
        Mass values (solar masses).
    w : ``float``
        Weight of the left (low-mass) peak. 
        default: 0.643
    muL : ``float``
        Mean of the left peak (solar masses). 
        default: 1.352
    sigmaL : ``float``
        Standard deviation of the left peak (solar masses). 
        default: 0.08
    muR : ``float``
        Mean of the right peak (solar masses). 
        default: 1.88
    sigmaR : ``float``
        Standard deviation of the right peak (solar masses). 
        default: 0.3
    mmin : ``float``
        Minimum mass (solar masses). 
        default: 1.0
    mmax : ``float``
        Maximum mass (solar masses). 
        default: 2.3

    Returns
    -------
    pdf : ``numpy.ndarray``
        Normalized probability density values.

    Examples
    --------
    >>> from ler.gw_source_population.prior_functions import bimodal_pdf
    >>> import numpy as np
    >>> m = np.array([1.2, 1.4, 1.8])
    >>> pdf_values = bimodal_pdf(m)
    >>> print(pdf_values)  # doctest: +SKIP
    """
    # left peak
    pdf_unnormL = np.exp(-((m - muL) ** 2) / (2 * sigmaL**2))
    normL = _compute_normalization_factor(muL, sigmaL, mmin, mmax)
    # right peak
    pdf_unnormR = np.exp(-((m - muR) ** 2) / (2 * sigmaR**2))
    normR = _compute_normalization_factor(muR, sigmaR, mmin, mmax)
    # total pdf
    pdf = w * pdf_unnormL / normL + (1 - w) * pdf_unnormR / normR

    return pdf

@njit(fastmath=True)
def bimodal_cdf(size, w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1.0, mmax=2.3):
    """
    Compute cumulative distribution function for BNS bimodal mass distribution.

    Parameters
    ----------
    size : ``int``
        Resolution of the mass grid.
    w : ``float``
        Weight of the left (low-mass) peak. 
        default: 0.643
    muL : ``float``
        Mean of the left peak (solar masses). 
        default: 1.352
    sigmaL : ``float``
        Standard deviation of the left peak (solar masses). 
        default: 0.08
    muR : ``float``
        Mean of the right peak (solar masses). 
        default: 1.88
    sigmaR : ``float``
        Standard deviation of the right peak (solar masses). 
        default: 0.3
    mmin : ``float``
        Minimum mass (solar masses). 
        default: 1.0
    mmax : ``float``
        Maximum mass (solar masses). 
        default: 2.3

    Returns
    -------
    cdf_values : ``numpy.ndarray``
        Cumulative distribution function values.
    mass_grid : ``numpy.ndarray``
        Mass grid corresponding to CDF values.
    """
    mass_grid = np.linspace(mmin, mmax, size)
    pdf_values = bimodal_pdf(mass_grid, w, muL, sigmaL, muR, sigmaR, mmin, mmax)
    cdf_values, mass_grid, _ = cumulative_trapezoid(
        x=mass_grid, y=pdf_values, initial=0.0
    )
    
    return cdf_values, mass_grid

@njit(fastmath=True)
def binary_masses_BNS_bimodal_rvs(size, w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1.0, mmax=2.3, resolution=500):
    """
    Generate random samples of binary masses for BNS (bimodal Gaussian model).

    Uses inverse transform sampling to draw samples from a bimodal Gaussian
    distribution of neutron star masses.

    Parameters
    ----------
    size : ``int``
        Number of binary systems to sample.
    w : ``float``
        Weight of the left (low-mass) peak. 
        default: 0.643
    muL : ``float``
        Mean of the left peak (solar masses). 
        default: 1.352
    sigmaL : ``float``
        Standard deviation of the left peak (solar masses). 
        default: 0.08
    muR : ``float``
        Mean of the right peak (solar masses). 
        default: 1.88
    sigmaR : ``float``
        Standard deviation of the right peak (solar masses). 
        default: 0.3
    mmin : ``float``
        Minimum mass (solar masses). 
        default: 1.0
    mmax : ``float``
        Maximum mass (solar masses). 
        default: 2.3
    resolution : ``int``
        Resolution of the mass grid for CDF computation.
        default: 500

    Returns
    -------
    mass_samples : ``numpy.ndarray``
        Array of mass samples (solar masses).
    """
    cdf_values, mass_grid = bimodal_cdf(resolution, w, muL, sigmaL, muR, sigmaR, mmin, mmax)
    
    return inverse_transform_sampler(size=size, cdf=cdf_values, x=mass_grid)

# @njit
# def _inverse_transform_sampler_m1m2(size, cdf_values, x):
#     """
#     Sample m1 and m2 using inverse transform sampling for BNS.

#     This is a helper function for the BNS Alsing mass distribution function.

#     Parameters
#     ----------
#     size : ``int``
#         Number of samples to draw.
#     cdf_values : ``numpy.ndarray``
#         Cumulative distribution function values.
#     x : ``numpy.ndarray``
#         Mass values corresponding to the CDF.

#     Returns
#     -------
#     m1 : ``numpy.ndarray``
#         Primary mass samples (Msun).
#     m2 : ``numpy.ndarray``
#         Secondary mass samples (Msun).

#     Examples
#     --------
#     >>> from ler.gw_source_population import _inverse_transform_sampler_m1m2
#     >>> m1, m2 = _inverse_transform_sampler_m1m2(size=1000, cdf_values=cdf_values, x=mass_arr)
#     """
#     m1 = inverse_transform_sampler(size=size, cdf=cdf_values, x=x)
#     m2 = inverse_transform_sampler(size=size, cdf=cdf_values, x=x)
#     # swap m1 and m2 if m1 < m2
#     idx = m1 < m2
#     m1[idx], m2[idx] = m2[idx], m1[idx]

#     return m1, m2

# ------------------------------
# broken_powerlaw
# ------------------------------

@njit(fastmath=True)
def _broken_powerlaw_unnormalized(m, mminbh=26., mmaxbh=125., alpha_1=6.75, alpha_2=0., b=0.5, delta_m=5.):
    """
    Compute unnormalized PDF for broken power-law distribution.

    Parameters
    ----------
    m : ``numpy.ndarray``
        Mass values.
    mminbh : ``float``
        Minimum BH mass.
    mmaxbh : ``float``
        Maximum BH mass.
    alpha_1 : ``float``
        Power-law index below break.
    alpha_2 : ``float``
        Power-law index above break.
    b : ``float``
        Break location parameter.
    delta_m : ``float``
        Smoothing width.

    Returns
    -------
    pdf_unnormalized : ``numpy.ndarray``
        Unnormalized PDF values.
    """
    mbreak = mminbh + b * (mmaxbh - mminbh)
    idx_1 = (m > mminbh) & (m < mbreak)
    idx_2 = (m >= mbreak) & (m < mmaxbh)

    pdf_unnormalized = np.zeros_like(m)
    pdf_unnormalized[idx_1] = powerlaw_with_smoothing(m[idx_1], m[idx_1], mminbh, -alpha_1, delta_m)
    norm_1 = pdf_unnormalized[idx_1][np.sum(idx_1)-1]
    pdf_unnormalized[idx_2] = powerlaw_with_smoothing(m[idx_2], m[idx_2], mminbh, -alpha_2, delta_m)
    norm_2 = pdf_unnormalized[idx_2][0]
    pdf_unnormalized[idx_2] = pdf_unnormalized[idx_2] * (norm_1 / norm_2)

    return pdf_unnormalized

def broken_powerlaw_pdf(m, mminbh=26., mmaxbh=125., alpha_1=6.75, alpha_2=0., b=0.5, delta_m=5.):
    """
    Compute normalized PDF for broken power-law distribution.

    Parameters
    ----------
    m : ``numpy.ndarray``
        Mass values.
    mminbh : ``float``
        Minimum BH mass.
    mmaxbh : ``float``
        Maximum BH mass.
    alpha_1 : ``float``
        Power-law index below break.
    alpha_2 : ``float``
        Power-law index above break.
    b : ``float``
        Break location parameter.
    delta_m : ``float``
        Smoothing width.

    Returns
    -------
    pdf : ``numpy.ndarray``
        Normalized PDF values.
    """
    # compute unnormalized PDF on a fine grid for numerical integration
    m_try = np.geomspace(mminbh, mmaxbh, 1000)
    pdf_unnormalized = _broken_powerlaw_unnormalized(m_try, mminbh=mminbh, mmaxbh=mmaxbh, alpha_1=alpha_1, alpha_2=alpha_2, b=b, delta_m=delta_m)
    _, _, norm = cumulative_trapezoid(y=pdf_unnormalized, x=m_try, initial=0.0)

    pdf = _broken_powerlaw_unnormalized(m, mminbh=mminbh, mmaxbh=mmaxbh, alpha_1=alpha_1, alpha_2=alpha_2, b=b, delta_m=delta_m) / norm

    return pdf

@njit(fastmath=True)
def _gaussian_pdf(m, mu_g=32.27, sigma_g=3.88):
    """
    Compute Gaussian distribution.

    Parameters
    ----------
    m : ``numpy.ndarray``
        Mass values.
    mu_g : ``float``
        Mean of the Gaussian.
    sigma_g : ``float``
        Standard deviation.

    Returns
    -------
    pdf : ``numpy.ndarray``
        Gaussian PDF values.
    """

    normalization = 1.0 / (sigma_g * np.sqrt(2 * np.pi))
    exponent = -0.5 * ((m - mu_g) / sigma_g) ** 2
    pdf = normalization * np.exp(exponent)
    return pdf

# -----------------------
# Spin Tilt Angles
# -----------------------
def gaussian_plus_isotropic_pdf(x, mu_t=0.426, sigma_t=1.222, zeta=0.652):
    """
    1D marginal PDF for a single cosine tilt x = cos(theta).

    p(x) = zeta * TruncNorm[-1,1](x | mu_t, sigma_t) + (1-zeta)/2
    """

    return (
        zeta * truncated_normal_pdf(x, mu_t, sigma_t, -1.0, 1.0)
        + (1.0 - zeta) * 0.5
    )

def gaussian_plus_isotropic_joint_pdf(x1, x2, mu_t=0.426, sigma_t=1.222, zeta=0.652):
    """
    2D joint PDF for (x1, x2) = (cos(theta1), cos(theta2)).

    p(x1, x2) =
        zeta * TN(x1) * TN(x2)
        + (1-zeta) / 4
    on [-1,1]^2.
    """

    g1 = truncated_normal_pdf(x1, mu_t, sigma_t, -1.0, 1.0)
    g2 = truncated_normal_pdf(x2, mu_t, sigma_t, -1.0, 1.0)

    pdf = zeta * g1 * g2 + (1.0 - zeta) * 0.25
    return pdf

# -----------------------
# common helper functions
# -----------------------
@njit(fastmath=True)
def powerlaw_pdf(x, alpha=-7.7, x_min=1.0, x_max=2.5):
    """
        Compute normalized power-law distribution.
        p(x) ∝ x^{-alpha}, x in [x_min, x_max]

        Parameters
        ----------
        x : ``numpy.ndarray``
            Input values.
        alpha : ``float``
            Power-law spectral index.
        x_min : ``float``
            Minimum value.
        x_max : ``float``
            Maximum value.
            Maximum mass.

        Returns
        -------
        pdf : ``numpy.ndarray``
            Normalized power-law PDF.
    """
    normalization = (x_max ** (-alpha + 1)) / (-alpha + 1) - (x_min ** (-alpha + 1)) / (-alpha + 1)
    pdf = x ** (-alpha) / normalization
    return pdf

@njit(fastmath=True)
def powerlaw_rvs(size, alpha, x_min, x_max):
    """
    Inverse transform sampling for a power-law distribution.
    
    p(x) ∝ x^{-alpha}, x in [x_min, x_max]

    Parameters
    ----------
    size : ``int``
        Number of samples to generate.
    alpha : ``float``
        Power-law index (alpha).
    x_min : ``float``
        Minimum value (lower bound).
    x_max : ``float``
        Maximum value (upper bound).

    Returns
    -------
    x : ``numpy.ndarray``
        Array of sampled values.
    """

    u = np.random.uniform(0, 1, size)

    if alpha == 1.0:
        # Special case α=1
        x = x_min * (x_max / x_min) ** u
    elif alpha == 0.0:
        # Special case α=0 (uniform distribution)
        x = x_min + (x_max - x_min) * u
    else:
        pow1 = 1.0 - alpha
        x_min_pow = x_min**pow1
        x_max_pow = x_max**pow1
        x = (u * (x_max_pow - x_min_pow) + x_min_pow) ** (1.0 / pow1)

    return x

@njit()
def truncated_normal_pdf(x, mu, sigma, x_min, x_max=np.inf):
    """
        Compute left-truncated or left-and-right-truncated normal probability density function.

        Evaluates the truncated normal distribution $N_{[x_min, x_max]}(x | μ, σ)$,
        which is a Gaussian distribution with support only between a minimum and maximum value.
        If x_max is not provided (or set to np.inf), it defaults to left-only truncation.

        Parameters
        ----------
        x : ``numpy.ndarray``
            Input values.
        mu : ``float``
            Mean of the Gaussian distribution.
        sigma : ``float``
            Standard deviation of the Gaussian distribution.
        x_min : ``float``
            Minimum value (left truncation point).
        x_max : ``float``, optional
            Maximum value (right truncation point). Default is np.inf (no right truncation).

        Returns
        -------
        pdf : ``numpy.ndarray``
            Probability density values, 0 for $x < x_min$ or $x > x_max$.
    """
    # 1. Compute CDF at the lower bound
    a = (x_min - mu) / sigma
    cdf_a = 0.5 * (1.0 + math.erf(a / math.sqrt(2.0)))
    
    # 2. Compute CDF at the upper bound
    if x_max == np.inf:
        cdf_b = 1.0
    else:
        b = (x_max - mu) / sigma
        cdf_b = 0.5 * (1.0 + math.erf(b / math.sqrt(2.0)))
        
    # 3. Calculate normalization factor
    norm = cdf_b - cdf_a
    
    # Safety guard against numerical precision issues yielding 0 normalization
    if norm <= 0.0:
        return np.zeros_like(x)
        
    # 4. Evaluate Gaussian
    z = (x - mu) / sigma
    pdf = np.exp(-0.5 * z * z) / (sigma * math.sqrt(2.0 * math.pi) * norm)
    
    # 5. Apply bounds
    mask = (x >= x_min) & (x <= x_max)
    return np.where(mask, pdf, 0.0)

@njit()
def _standard_normal_ppf(p):
    """
    Approximate inverse CDF (quantile) of the standard normal distribution.

    Uses the Acklam rational approximation, accurate enough for sampling.
    """
    if p <= 0.0:
        return -np.inf
    if p >= 1.0:
        return np.inf

    plow = 0.02425
    phigh = 1.0 - plow

    a1 = -39.69683028665376
    a2 = 220.9460984245205
    a3 = -275.9285104469687
    a4 = 138.3577518672690
    a5 = -30.66479806614716
    a6 = 2.506628277459239

    b1 = -54.47609879822406
    b2 = 161.5858368580409
    b3 = -155.6989798598866
    b4 = 66.80131188771972
    b5 = -13.28068155288572

    c1 = -0.007784894002430293
    c2 = -0.3223964580411365
    c3 = -2.400758277161838
    c4 = -2.549732539343734
    c5 = 4.374664141464968
    c6 = 2.938163982698783

    d1 = 0.007784695709041462
    d2 = 0.3224671290700398
    d3 = 2.445134137142996
    d4 = 3.754408661907416

    if p < plow:
        q = math.sqrt(-2.0 * math.log(p))
        return (
            (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6)
            / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0)
        )

    if p > phigh:
        q = math.sqrt(-2.0 * math.log(1.0 - p))
        return -(
            (((((c1 * q + c2) * q + c3) * q + c4) * q + c5) * q + c6)
            / ((((d1 * q + d2) * q + d3) * q + d4) * q + 1.0)
        )

    q = p - 0.5
    r = q * q
    return (
        (((((a1 * r + a2) * r + a3) * r + a4) * r + a5) * r + a6) * q
        / (((((b1 * r + b2) * r + b3) * r + b4) * r + b5) * r + 1.0)
    )

@njit()
def _truncated_normal_ppf(u, mu, sigma, x_min, x_max=np.inf):
    """
    Inverse CDF (quantile) for a truncated normal distribution.
    """
    a = (x_min - mu) / sigma
    cdf_a = 0.5 * (1.0 + math.erf(a / math.sqrt(2.0)))

    if x_max == np.inf:
        cdf_b = 1.0
    else:
        b = (x_max - mu) / sigma
        cdf_b = 0.5 * (1.0 + math.erf(b / math.sqrt(2.0)))

    norm = cdf_b - cdf_a
    if norm <= 0.0:
        return x_min

    if u <= 0.0:
        return x_min
    if u >= 1.0:
        return x_max

    p = cdf_a + u * norm
    z = _standard_normal_ppf(p)
    x = mu + sigma * z

    if x < x_min:
        return x_min
    if x > x_max:
        return x_max
    return x

@njit()
def truncated_normal_rvs(size, mu, sigma, x_min, x_max=np.inf):
    """
    Draw samples from a truncated normal distribution using analytical inverse CDF.
    """
    u = np.random.uniform(0.0, 1.0, size)
    samples = np.zeros(size)
    for i in range(size):
        samples[i] = _truncated_normal_ppf(u[i], mu, sigma, x_min, x_max)

    return samples

# -----------------------
@njit()
def _smoothing_S(m, mmin, delta_m, threshold=709.0):
    """
    Compute low-mass smoothing function to avoid sharp cutoffs.

    Parameters
    ----------
    m : ``numpy.ndarray``
        Mass values.
    mmin : ``float``
        Minimum mass.
    delta_m : ``float``
        Smoothing width.
    threshold : ``float``
        Maximum exponent to avoid overflow. \n
        default: 709.0

    Returns
    -------
    s : ``numpy.ndarray``
        Smoothing function values.
    """
    s = np.zeros_like(m)

    # region where smoothing is not needed: m >= mmin + delta_m
    idx_2 = m >= mmin + delta_m
    s[idx_2] = 1.0

    # region where smoothing is applied: mmin <= m < mmin + delta_m
    idx_1 = (m >= mmin) & (m < mmin + delta_m)
    mprime = m[idx_1] - mmin
    exponent = delta_m / mprime + delta_m / (mprime - delta_m)

    # safe exponentiation only where exponent is below threshold
    safe_idx = exponent <= threshold
    s_vals = np.zeros_like(mprime)
    s_vals[safe_idx] = 1.0 / (np.exp(exponent[safe_idx]) + 1.0)

    s[idx_1] = s_vals

    return s

@njit(fastmath=True)
def powerlaw_with_smoothing(q, m, mmin, beta, delta_m):
    """
    Compute power-law distribution with low-mass smoothing.

    """

    return q ** (beta) * _smoothing_S(m, mmin, delta_m)

def _njit_checks(
    zs_rvs_,
    m1_rvs_,
    q_rvs_,
    m2_rvs_,
    tc_rvs_,
    ra_rvs_,
    dec_rvs_,
    phase_rvs_,
    psi_rvs_,
    theta_jn_rvs_,
    a_1_rvs_=None,
    a_2_rvs_=None,
    tilt_1_rvs_=None,
    tilt_2_rvs_=None,
    phi_12_rvs_=None,
    phi_jl_rvs_=None,
    spin_zero=False,
    spin_precession=False,
):
    """
    Check and wrap sampler functions for JIT compatibility.

    """

    if m1_rvs_.__code__.co_argcount == 1:
        if is_njitted(m1_rvs_):
            @njit
            def m1_rvs(size, zs=None):
                return m1_rvs_(size)
        else:
            def m1_rvs(size, zs=None):
                return m1_rvs_(size)
    else:
        m1_rvs = m1_rvs_

    if m2_rvs_ is None:
        if q_rvs_.__code__.co_argcount == 1:
            if is_njitted(q_rvs_):
                @njit
                def q_rvs(size, m1=None):
                    return q_rvs_(size)
            else:
                def q_rvs(size, m1=None):
                    return q_rvs_(size)
        else:
            q_rvs = q_rvs_

    else:
        if m2_rvs_.__code__.co_argcount == 1:
            if is_njitted(m2_rvs_):
                @njit
                def m2_rvs(size, zs=None):
                    return m2_rvs_(size)
            else:
                def m2_rvs(size, zs=None):
                    return m2_rvs_(size)
        else:
            m2_rvs = m2_rvs_

    if not spin_zero:
        if a_2_rvs_.__code__.co_argcount == 1:
            if is_njitted(a_2_rvs_):
                @njit
                def a_2_rvs(size, a_1=None):
                    return a_2_rvs_(size)
            else:
                def a_2_rvs(size, a_1=None):
                    return a_2_rvs_(size)
        else:
            a_2_rvs = a_2_rvs_

        if spin_precession:

            if tilt_2_rvs_.__code__.co_argcount == 1:
                if is_njitted(tilt_2_rvs_):
                    @njit
                    def tilt_2_rvs(size, tilt_1=None):
                        return tilt_2_rvs_(size)
                else:
                    def tilt_2_rvs(size, tilt_1=None):
                        return tilt_2_rvs_(size)
            else:
                tilt_2_rvs = tilt_2_rvs_

    # Build dictionary of wrapped functions
    dict_ = {
        "zs_rvs": zs_rvs_,
        "m1_rvs": m1_rvs,
        "q_rvs": q_rvs if m2_rvs_ is None else None,
        "m2_rvs": m2_rvs if m2_rvs_ is not None else None,
        "tc_rvs": tc_rvs_,
        "ra_rvs": ra_rvs_,
        "dec_rvs": dec_rvs_,
        "phase_rvs": phase_rvs_,
        "psi_rvs": psi_rvs_,
        "theta_jn_rvs": theta_jn_rvs_,
        "a_1_rvs": a_1_rvs_ if not spin_zero else None,
        "a_2_rvs": a_2_rvs if not spin_zero else None,
        "tilt_1_rvs": tilt_1_rvs_ if (not spin_zero and spin_precession) else None,
        "tilt_2_rvs": tilt_2_rvs if (not spin_zero and spin_precession) else None,
        "phi_12_rvs": phi_12_rvs_ if (not spin_zero and spin_precession) else None,
        "phi_jl_rvs": phi_jl_rvs_ if (not spin_zero and spin_precession) else None,
    }

    # Check if all functions are JIT compiled
    use_njit_sampler = True
    for key, value in dict_.items():
        if not is_njitted(value) and value is not None:
            print(f"Warning: {key} is not njitted.")
            use_njit_sampler = False

    return use_njit_sampler, dict_

def create_gw_parameters_sampler(
    zs_rvs,
    m1_rvs,
    q_rvs,
    m2_rvs,
    tc_rvs,
    ra_rvs,
    dec_rvs,
    phase_rvs,
    psi_rvs,
    theta_jn_rvs,
    a_1_rvs,
    a_2_rvs,
    tilt_1_rvs,
    tilt_2_rvs,
    phi_12_rvs,
    phi_jl_rvs,
    use_njit_sampler=True,
    spin_zero=False,
    spin_precession=False,
):
    
    if m2_rvs is None:

        if spin_zero:
            def sampler(size):
                zs = zs_rvs(size)
                m1 = m1_rvs(size, zs)
                q = q_rvs(size, m1)
                m2 = q * m1
                tc = tc_rvs(size)
                ra = ra_rvs(size)
                dec = dec_rvs(size)
                phase = phase_rvs(size)
                psi = psi_rvs(size)
                theta_jn = theta_jn_rvs(size)

                # swap m1 and m2 to ensure m1 >= m2
                m1, m2 = np.where(m1 > m2, m1, m2), np.where(m1 > m2, m2, m1)

                return zs, m1, m2, tc, ra, dec, phase, psi, theta_jn

        else:
            
            if spin_precession:
                def sampler(size):
                    zs = zs_rvs(size)
                    m1 = m1_rvs(size, zs)
                    q = q_rvs(size, m1)
                    m2 = q * m1
                    tc = tc_rvs(size)
                    ra = ra_rvs(size)
                    dec = dec_rvs(size)
                    phase = phase_rvs(size)
                    psi = psi_rvs(size)
                    theta_jn = theta_jn_rvs(size)
                    a_1 = a_1_rvs(size)
                    a_2 = a_2_rvs(size)
                    tilt_1 = tilt_1_rvs(size)
                    tilt_2 = tilt_2_rvs(size, tilt_1)
                    phi_12 = phi_12_rvs(size)
                    phi_jl = phi_jl_rvs(size)

                    # swap m1 and m2 to ensure m1 >= m2
                    m1, m2 = np.where(m1 > m2, m1, m2), np.where(m1 > m2, m2, m1)

                    return zs, m1, m2, tc, ra, dec, phase, psi, theta_jn, a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl
            else:
                def sampler(size):
                    zs = zs_rvs(size)
                    m1 = m1_rvs(size, zs)
                    q = q_rvs(size, m1)
                    m2 = q * m1
                    tc = tc_rvs(size)
                    ra = ra_rvs(size)
                    dec = dec_rvs(size)
                    phase = phase_rvs(size)
                    psi = psi_rvs(size)
                    theta_jn = theta_jn_rvs(size)
                    a_1 = a_1_rvs(size)
                    a_2 = a_2_rvs(size)

                    # swap m1 and m2 to ensure m1 >= m2
                    m1, m2 = np.where(m1 > m2, m1, m2), np.where(m1 > m2, m2, m1)

                    return zs, m1, m2, tc, ra, dec, phase, psi, theta_jn, a_1, a_2

    else: # if m2_rvs is not None
        if spin_zero:
            def sampler(size):
                zs = zs_rvs(size)
                m1 = m1_rvs(size, zs)
                m2 = m2_rvs(size, m1)
                tc = tc_rvs(size)
                ra = ra_rvs(size)
                dec = dec_rvs(size)
                phase = phase_rvs(size)
                psi = psi_rvs(size)
                theta_jn = theta_jn_rvs(size)

                # swap m1 and m2 to ensure m1 >= m2
                m1, m2 = np.where(m1 > m2, m1, m2), np.where(m1 > m2, m2, m1)

                return zs, m1, m2, tc, ra, dec, phase, psi, theta_jn

        else:
            if spin_precession:
                def sampler(size):
                    zs = zs_rvs(size)
                    m1 = m1_rvs(size, zs)
                    m2 = m2_rvs(size, m1)
                    tc = tc_rvs(size)
                    ra = ra_rvs(size)
                    dec = dec_rvs(size)
                    phase = phase_rvs(size)
                    psi = psi_rvs(size)
                    theta_jn = theta_jn_rvs(size)
                    a_1 = a_1_rvs(size, zs)
                    a_2 = a_2_rvs(size, a_1)
                    tilt_1 = tilt_1_rvs(size, zs)
                    tilt_2 = tilt_2_rvs(size, tilt_1)
                    phi_12 = phi_12_rvs(size, zs)
                    phi_jl = phi_jl_rvs(size, zs)

                    # swap m1 and m2 to ensure m1 >= m2
                    m1, m2 = np.where(m1 > m2, m1, m2), np.where(m1 > m2, m2, m1)

                    return zs, m1, m2, tc, ra, dec, phase, psi, theta_jn, a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl
            else:
                def sampler(size):
                    zs = zs_rvs(size)
                    m1 = m1_rvs(size, zs)
                    m2 = m2_rvs(size, m1)
                    tc = tc_rvs(size)
                    ra = ra_rvs(size)
                    dec = dec_rvs(size)
                    phase = phase_rvs(size)
                    psi = psi_rvs(size)
                    theta_jn = theta_jn_rvs(size)
                    a_1 = a_1_rvs(size, zs)
                    a_2 = a_2_rvs(size, a_1)

                    # swap m1 and m2 to ensure m1 >= m2
                    m1, m2 = np.where(m1 > m2, m1, m2), np.where(m1 > m2, m2, m1)

                    return zs, m1, m2, tc, ra, dec, phase, psi, theta_jn, a_1, a_2
        

    if use_njit_sampler:
        return njit(sampler)
    else:
        return sampler
