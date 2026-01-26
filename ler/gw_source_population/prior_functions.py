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
from numba import njit
from scipy.interpolate import CubicSpline
from astropy.cosmology import LambdaCDM

from ..utils import inverse_transform_sampler, sample_from_powerlaw_distribution

# ------------------------------
# Merger rate density functions
# ------------------------------
@njit
def merger_rate_density_bbh_oguri2018_function(zs, R0=19 * 1e-9, b2=1.6, b3=2.1, b4=30):
    """
    Compute the merger rate density for PopI/II BBH.

    Reference: Oguri et al. (2018). The output is in detector frame and is
    unnormalized.

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
    >>> from ler.gw_source_population import merger_rate_density_bbh_oguri2018
    >>> rate_density = merger_rate_density_bbh_oguri2018(zs=np.array([0.1]))
    """
    return R0 * (b4 + 1) * np.exp(b2 * zs) / (b4 + np.exp(b3 * zs))

@njit
def merger_rate_density_bbh_popIII_ken2022_function(zs, n0=19.2 * 1e-9, aIII=0.66, bIII=0.3, zIII=11.6):
    """
    Compute the unnormalized merger rate density for PopIII BBH.

    Reference: Ng et al. (2022). The output is in detector frame and is
    unnormalized.

    Parameters
    ----------
    zs : ``float`` or ``numpy.ndarray``
        Source redshifts.
    n0 : ``float``
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
    >>> from ler.gw_source_population import merger_rate_density_bbh_popIII_ken2022
    >>> rate_density = merger_rate_density_bbh_popIII_ken2022(zs=np.array([0.1]))
    """
    return (
        n0
        * np.exp(aIII * (zs - zIII))
        / (bIII + aIII * np.exp((aIII + bIII) * (zs - zIII)))
    )

@njit
def merger_rate_density_madau_dickinson2014_function(zs, R0=19 * 1e-9, a=0.015, b=2.7, c=2.9, d=5.6):
    """
    Compute the merger rate density for BBH using Madau & Dickinson (2014) model.

    Reference: Eqn. 15 of https://arxiv.org/pdf/1403.0007

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
    >>> from ler.gw_source_population import merger_rate_density_madau_dickinson2014
    >>> rate_density = merger_rate_density_madau_dickinson2014(zs=np.array([0.1]))
    """

    density_helper = lambda zs: sfr_madau_dickinson2014(  
        zs=zs, 
        a=a, 
        b=b, 
        c=c,
        d=d,
    )

    density_zs = R0 * density_helper(zs)/ density_helper(np.array([0.]))[0]

    return density_zs

@njit
def merger_rate_density_madau_dickinson_belczynski_ng_function(zs, R0=19 * 1e-9, alpha_F=2.57, beta_F=5.83, c_F=3.36):
    """
    Compute BBH merger rate density following Ng et al. (2021).

    This model uses a Madau-Dickinson-like functional form to fit the 
    merger rate density of field BHs, accounting for time delays and 
    metallicity effects. Coefficients from Madau & Dickinson (2014) are translated as: B-> alpha_F, D-> beta_F, C-> c_F.

    density(zs) ∝ (1 + zs) ** alpha_F / (1 + ((1 + zs) / c_F) ** beta_F)

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
    >>> from ler.gw_source_population import merger_rate_density_madau_dickinson_belczynski_ng
    >>> rate_density = merger_rate_density_madau_dickinson_belczynski_ng(zs=np.array([0.1]))
    """

    density_helper = lambda zs: sfr_madau_dickinson2014(  
        zs=zs, 
        a=1.0, 
        b=alpha_F, 
        c=c_F,
        d=beta_F,
    )
    density_zs = R0 * density_helper(zs)/ density_helper(np.array([0.]))[0]

    return density_zs

def merger_rate_density_bbh_primordial_ken2022_function(zs, cosmology=None, n0=0.044 * 1e-9, t0=13.786885302009708):
    """
    Compute the merger rate density for Primordial BBH.

    Reference: Ng et al. (2022). The output is in detector frame and is
    unnormalized.

    Parameters
    ----------
    zs : ``float`` or ``numpy.ndarray``
        Source redshifts.
    cosmology : ``astropy.cosmology`` or ``None``
        Cosmology object for age calculations. \n
        default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)
    n0 : ``float``
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
    >>> from ler.gw_source_population import merger_rate_density_bbh_primordial_ken2022
    >>> rate_density = merger_rate_density_bbh_primordial_ken2022(zs=np.array([0.1]))
    """
    if cosmology is None:
        cosmology = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)

    rate_density = n0 * (cosmology.age(z=zs).value / t0) ** (-34 / 37)
    return rate_density


def sfr_madau_fragos2017_with_bbh_td(zs, R0=19 * 1e-9):
    """
    Compute the merger rate density for BBH. This is computed from star formation rate, Madau & Fragos (2017), with an additional time delay. This function is relies on pre-generated data points. 

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
    Compute the merger rate density for BBH. This is computed from star formation rate, Madau & Dickinson (2014), with an additional time delay. This function is relies on pre-generated data points. 

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
    Compute the merger rate density for BNS. This is computed from star formation rate, Madau & Fragos (2017), with an additional time delay. This function is relies on pre-generated data points. 

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
    Compute the merger rate density for BNS. This is computed from star formation rate, Madau & Dickinson (2014), with an additional time delay. This function is relies on pre-generated data points. 

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
@njit
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

@njit
def sfr_madau_dickinson2014(zs, a=0.015, b=2.7, c=2.9, d=5.6):
    """
    Compute star formation rate using Madau & Dickinson (2014) model.

    Reference: Eqn. 15 of https://arxiv.org/pdf/1403.0007

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
@njit
def binary_masses_BBH_popIII_lognormal_rvs(size, m_min=1.0, m_max=100.0, Mc=20.0, sigma=0.3, chunk_size=10000):
    """
    Sample from a lognormal distribution in 2D mass space.

    Reference: Ng et al. (2022). This is a helper function for PopIII BBH
    and primordial BBH merger rate density distribution functions.

    Parameters
    ----------
    size : ``int``
        Number of samples to draw.
    m_min : ``float``
        Minimum mass (Msun). \n
        default: 1.0
    m_max : ``float``
        Maximum mass (Msun). \n
        default: 100.0
    Mc : ``float``
        Characteristic mass scale (Msun). \n
        default: 20.0
    sigma : ``float``
        Width of the distribution. \n
        default: 0.3
    chunk_size : ``int``
        Number of samples per rejection sampling chunk. \n
        default: 10000

    Returns
    -------
    m1_sample : ``numpy.ndarray``
        Primary mass samples (Msun).
    m2_sample : ``numpy.ndarray``
        Secondary mass samples (Msun).

    Examples
    --------
    >>> from ler.gw_source_population import binary_masses_BBH_popIII_lognormal
    >>> m1, m2 = binary_masses_BBH_popIII_lognormal(size=1000)
    """
    # mass function (Eqn. 1 of Ng et al. 2022)
    psi = lambda m: np.exp(-np.log(m / Mc) ** 2 / (2 * sigma**2)) / (
        np.sqrt(2 * np.pi) * sigma * m
    )
    # probability density function (Eqn. 4 of Ng et al. 2022)
    pdf = (
        lambda m1, m2: (m1 + m2) ** (36 / 37)
        * (m1 * m2) ** (32 / 37)
        * psi(m1)
        * psi(m2)
    )

    # rejection sampling initialization
    m1 = np.random.uniform(m_min, m_max, chunk_size)
    m2 = np.random.uniform(m_min, m_max, chunk_size)
    z = pdf(m1, m2)
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
        idx = z_try < pdf(m1_try, m2_try)
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

@njit
def binary_masses_BBH_primordial_lognormal_rvs(size, m_min=1.0, m_max=100.0, Mc=20.0, sigma=0.3, chunk_size=10000):
    """
    Sample from a lognormal distribution in 2D mass space.

    Based on Eqn. 1 and 4 of Ng et al. 2022 for primordial black holes.

    Parameters
    ----------
    size : ``int``
        Number of samples to draw.
    m_min : ``float``
        Minimum mass (Msun). \n
        default: 1.0
    m_max : ``float``
        Maximum mass (Msun). \n
        default: 100.0
    Mc : ``float``
        Characteristic mass scale (Msun). \n
        default: 20.0
    sigma : ``float``
        Width of the distribution. \n
        default: 0.3
    chunk_size : ``int``
        Number of samples per rejection sampling chunk. \n
        default: 10000

    Returns
    -------
    m1_sample : ``numpy.ndarray``
        Primary mass samples (Msun).
    m2_sample : ``numpy.ndarray``
        Secondary mass samples (Msun).

    Examples
    --------
    >>> from ler.gw_source_population import binary_masses_BBH_primordial_lognormal
    >>> m1, m2 = binary_masses_BBH_primordial_lognormal(size=1000)
    """
    # mass function (Eqn. 1 of Ng et al. 2022)
    psi = lambda m: np.exp(-np.log(m / Mc) ** 2 / (2 * sigma**2)) / (  
        np.sqrt(2 * np.pi) * sigma * m
    )
    # probability density function (Eqn. 4 of Ng et al. 2022)
    pdf = (  
        lambda m1, m2: (m1 + m2) ** (36 / 37)
        * (m1 * m2) ** (32 / 37)
        * psi(m1)
        * psi(m2)
    )

    # rejection sampling initialization
    m1 = np.random.uniform(m_min, m_max, chunk_size)
    m2 = np.random.uniform(m_min, m_max, chunk_size)
    z = pdf(m1, m2)
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
        idx = z_try < pdf(m1_try, m2_try)
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

@njit
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

@njit
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

@njit
def _bns_bimodal_pdf(m, w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1.0, mmax=2.3):
    """
    Compute the bimodal Gaussian PDF for BNS mass distribution.

    Parameters
    ----------
    m : ``float`` or ``numpy.ndarray``
        Mass values (Msun).
    w : ``float``
        Weight of the left (low-mass) peak. \n
        default: 0.643
    muL : ``float``
        Mean of the left peak (Msun). \n
        default: 1.352
    sigmaL : ``float``
        Standard deviation of the left peak (Msun). \n
        default: 0.08
    muR : ``float``
        Mean of the right peak (Msun). \n
        default: 1.88
    sigmaR : ``float``
        Standard deviation of the right peak (Msun). \n
        default: 0.3
    mmin : ``float``
        Minimum mass (Msun). \n
        default: 1.0
    mmax : ``float``
        Maximum mass (Msun). \n
        default: 2.3

    Returns
    -------
    pdf : ``float`` or ``numpy.ndarray``
        Probability density values.
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

@njit
def binary_masses_BNS_bimodal_rvs(size, w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1.0, mmax=2.3, resolution=500):
    """
    Sample BNS masses from bimodal Gaussian distribution.

    Based on Will M. Farr et al. 2020 Eqn. 6 for neutron star mass
    distribution combining two Gaussian peaks.

    Parameters
    ----------
    size : ``int``
        Number of samples to draw.
    w : ``float``
        Weight of the left (low-mass) peak. \n
        default: 0.643
    muL : ``float``
        Mean of the left peak (Msun). \n
        default: 1.352
    sigmaL : ``float``
        Standard deviation of the left peak (Msun). \n
        default: 0.08
    muR : ``float``
        Mean of the right peak (Msun). \n
        default: 1.88
    sigmaR : ``float``
        Standard deviation of the right peak (Msun). \n
        default: 0.3
    mmin : ``float``
        Minimum mass (Msun). \n
        default: 1.0
    mmax : ``float``
        Maximum mass (Msun). \n
        default: 2.3
    resolution : ``int``
        Number of points to use for the CDF. \n
        default: 500
    
    Returns
    -------
    m1 : ``numpy.ndarray``
        Primary mass samples (Msun).
    m2 : ``numpy.ndarray``
        Secondary mass samples (Msun).
    
    """
    mass = np.linspace(mmin, mmax, resolution)
    pdf = _bns_bimodal_pdf(mass, w, muL, sigmaL, muR, sigmaR, mmin, mmax)
    cdf = np.cumsum(pdf) / np.sum(pdf)
    
    return _inverse_transform_sampler_m1m2(size, cdf, mass)

@njit
def _inverse_transform_sampler_m1m2(size, cdf_values, x):
    """
    Sample m1 and m2 using inverse transform sampling for BNS.

    This is a helper function for the BNS Alsing mass distribution function.

    Parameters
    ----------
    size : ``int``
        Number of samples to draw.
    cdf_values : ``numpy.ndarray``
        Cumulative distribution function values.
    x : ``numpy.ndarray``
        Mass values corresponding to the CDF.

    Returns
    -------
    m1 : ``numpy.ndarray``
        Primary mass samples (Msun).
    m2 : ``numpy.ndarray``
        Secondary mass samples (Msun).

    Examples
    --------
    >>> from ler.gw_source_population import _inverse_transform_sampler_m1m2
    >>> m1, m2 = _inverse_transform_sampler_m1m2(size=1000, cdf_values=cdf_values, x=mass_arr)
    """
    m1 = inverse_transform_sampler(size, cdf_values, x)
    m2 = inverse_transform_sampler(size, cdf_values, x)
    # swap m1 and m2 if m1 < m2
    idx = m1 < m2
    m1[idx], m2[idx] = m2[idx], m1[idx]

    return m1, m2

@njit
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

@njit
def _powerlaw_with_smoothing(m, mmin, alpha, delta_m):
    """
    Compute power-law distribution with low-mass smoothing.

    Parameters
    ----------
    m : ``numpy.ndarray``
        Mass values.
    mmin : ``float``
        Minimum mass.
    alpha : ``float``
        Power-law spectral index.
    delta_m : ``float``
        Smoothing width.

    Returns
    -------
    pdf : ``numpy.ndarray``
        Power-law PDF with smoothing.
    """
    s = _smoothing_S(m, mmin, delta_m)
    return m ** (-alpha) * s

@njit
def _broken_powerlaw_cdf(size=1000, mminbh=26., mmaxbh=125., alpha_1=6.75, alpha_2=0.0, b=0.5, delta_m=5.):
    """
    Compute CDF for broken power-law mass distribution.

    Parameters
    ----------
    size : ``int``
        Number of grid points. \n
        default: 1000
    mminbh : ``float``
        Minimum BH mass (Msun). \n
        default: 26.0
    mmaxbh : ``float``
        Maximum BH mass (Msun). \n
        default: 125.0
    alpha_1 : ``float``
        Power-law index below the break. \n
        default: 6.75
    alpha_2 : ``float``
        Power-law index above the break. \n
        default: 0.0
    b : ``float``
        Break location parameter (0-1). \n
        default: 0.5
    delta_m : ``float``
        Smoothing width (Msun). \n
        default: 5.0

    Returns
    -------
    cdf_values : ``numpy.ndarray``
        Normalized CDF values.
    """
    m_try = np.linspace(mminbh, mmaxbh, size)
    pdf_unnormalized = _broken_powerlaw_unnormalized(m_try, mminbh=mminbh, mmaxbh=mmaxbh, alpha_1=alpha_1, alpha_2=alpha_2, b=b, delta_m=delta_m)
    cdf_values = _cumulative_trapezoid(y=pdf_unnormalized, x=m_try, dx=1.0, initial=0.0)
    normalization = cdf_values[size-1]
    cdf_values /= normalization

    return cdf_values

@njit
def _sample_broken_powerlaw(size=1000, mminbh=26., mmaxbh=125., alpha_1=6.75, alpha_2=0., b=0.5, delta_m=5., normalization_size=1000):
    """
    Generate samples from the broken power-law mass distribution.

    Parameters
    ----------
    size : ``int``
        Number of samples to draw. \n
        default: 1000
    mminbh : ``float``
        Minimum BH mass (Msun). \n
        default: 26.0
    mmaxbh : ``float``
        Maximum BH mass (Msun). \n
        default: 125.0
    alpha_1 : ``float``
        Power-law index below the break. \n
        default: 6.75
    alpha_2 : ``float``
        Power-law index above the break. \n
        default: 0.0
    b : ``float``
        Break location parameter (0-1). \n
        default: 0.5
    delta_m : ``float``
        Smoothing width (Msun). \n
        default: 5.0
    normalization_size : ``int``
        Grid size for CDF computation. \n
        default: 1000

    Returns
    -------
    samples : ``numpy.ndarray``
        Mass samples (Msun).
    """
    cdf_values = _broken_powerlaw_cdf(size=normalization_size, mminbh=mminbh, mmaxbh=mmaxbh, alpha_1=alpha_1, alpha_2=alpha_2, b=b, delta_m=delta_m)
    cdf_values = cdf_values.astype(np.float64)

    x = np.linspace(mminbh, mmaxbh, normalization_size)
    idx = np.argwhere(cdf_values > 0)[0][0]
    cdf_values = cdf_values[idx:]
    x = x[idx:]
    samples = inverse_transform_sampler(size, cdf_values, x)

    return samples

@njit
def binary_masses_NSBH_broken_powerlaw_rvs(size=1000, mminbh=26., mmaxbh=125., alpha_1=6.75, alpha_2=0., b=0.5, delta_m=5., mminns=1.0, mmaxns=3.0, alphans=0.0, normalization_size=1000):
    """
    Generate NSBH mass samples from broken power-law (BH) and power-law (NS).

    Parameters
    ----------
    size : ``int``
        Number of samples to draw. \n
        default: 1000
    mminbh : ``float``
        Minimum BH mass (Msun). \n
        default: 26.0
    mmaxbh : ``float``
        Maximum BH mass (Msun). \n
        default: 125.0
    alpha_1 : ``float``
        BH power-law index below break. \n
        default: 6.75
    alpha_2 : ``float``
        BH power-law index above break. \n
        default: 0.0
    b : ``float``
        Break location parameter (0-1). \n
        default: 0.5
    delta_m : ``float``
        Smoothing width (Msun). \n
        default: 5.0
    mminns : ``float``
        Minimum NS mass (Msun). \n
        default: 1.0
    mmaxns : ``float``
        Maximum NS mass (Msun). \n
        default: 3.0
    alphans : ``float``
        NS power-law index. \n
        default: 0.0
    normalization_size : ``int``
        Grid size for CDF computation. \n
        default: 1000

    Returns
    -------
    m1_samples : ``numpy.ndarray``
        BH mass samples (Msun).
    m2_samples : ``numpy.ndarray``
        NS mass samples (Msun).
    """
    m1_samples = _sample_broken_powerlaw(size=size, mminbh=mminbh, mmaxbh=mmaxbh, alpha_1=alpha_1, alpha_2=alpha_2, b=b, delta_m=delta_m, normalization_size=normalization_size)
    m2_samples = sample_from_powerlaw_distribution(size, alphans, mminns, mmaxns)

    return m1_samples, m2_samples

@njit
def _broken_powerlaw_pdf(m, mminbh=26., mmaxbh=125., alpha_1=6.75, alpha_2=0., b=0.5, delta_m=5., normalization_size=1000):
    """
    Compute the normalized PDF for broken power-law mass distribution.

    Parameters
    ----------
    m : ``numpy.ndarray``
        Mass values to evaluate (Msun).
    mminbh : ``float``
        Minimum BH mass (Msun). \n
        default: 26.0
    mmaxbh : ``float``
        Maximum BH mass (Msun). \n
        default: 125.0
    alpha_1 : ``float``
        Power-law index below the break. \n
        default: 6.75
    alpha_2 : ``float``
        Power-law index above the break. \n
        default: 0.0
    b : ``float``
        Break location parameter (0-1). \n
        default: 0.5
    delta_m : ``float``
        Smoothing width (Msun). \n
        default: 5.0
    normalization_size : ``int``
        Grid size for normalization. \n
        default: 1000

    Returns
    -------
    pdf : ``numpy.ndarray``
        Normalized probability density values.
    """
    m_try = np.linspace(mminbh, mmaxbh, normalization_size)
    pdf_unnormalized = _broken_powerlaw_unnormalized(m_try, mminbh=mminbh, mmaxbh=mmaxbh, alpha_1=alpha_1, alpha_2=alpha_2, b=b, delta_m=delta_m)
    normalization = np.trapz(pdf_unnormalized, m_try)

    pdf_unnormalized = _broken_powerlaw_unnormalized(m, mminbh=mminbh, mmaxbh=mmaxbh, alpha_1=alpha_1, alpha_2=alpha_2, b=b, delta_m=delta_m)
    pdf = pdf_unnormalized / normalization

    return pdf

@njit
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
    pdf_unnormalized[idx_1] = _powerlaw_with_smoothing(m[idx_1], mminbh, alpha_1, delta_m)
    norm_1 = pdf_unnormalized[idx_1][np.sum(idx_1)-1]
    pdf_unnormalized[idx_2] = _powerlaw_with_smoothing(m[idx_2], mminbh, alpha_2, delta_m)
    norm_2 = pdf_unnormalized[idx_2][0]
    pdf_unnormalized[idx_2] = pdf_unnormalized[idx_2] * (norm_1 / norm_2)

    return pdf_unnormalized

@njit
def _powerlaw_B(m, alpha, mminbh, mmaxbh):
    """
    Compute normalized power-law distribution.

    Parameters
    ----------
    m : ``numpy.ndarray``
        Mass values.
    alpha : ``float``
        Power-law spectral index.
    mminbh : ``float``
        Minimum mass.
    mmaxbh : ``float``
        Maximum mass.

    Returns
    -------
    pdf : ``numpy.ndarray``
        Normalized power-law PDF.
    """
    normalization = (mmaxbh ** (-alpha + 1)) / (-alpha + 1) - (mminbh ** (-alpha + 1)) / (-alpha + 1)
    pdf = m ** (-alpha) / normalization
    return pdf

@njit
def _gaussian_G(m, mu_g, sigma_g):
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

@njit
def _powerlaw_gaussian_pdf(m, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m, normalization_size=1000):
    """
    Compute the normalized PDF for power-law + Gaussian mass model.

    Parameters
    ----------
    m : ``numpy.ndarray``
        Mass values to evaluate (Msun).
    mminbh : ``float``
        Minimum BH mass (Msun).
    mmaxbh : ``float``
        Maximum BH mass (Msun).
    alpha : ``float``
        Power-law spectral index.
    mu_g : ``float``
        Mean of the Gaussian peak (Msun).
    sigma_g : ``float``
        Standard deviation of the Gaussian peak (Msun).
    lambda_peak : ``float``
        Fraction of sources in the Gaussian component (0-1).
    delta_m : ``float``
        Low-mass smoothing width (Msun).
    normalization_size : ``int``
        Grid size for normalization. \n
        default: 1000

    Returns
    -------
    pdf : ``numpy.ndarray``
        Normalized probability density values.
    """
    m_try = np.linspace(mminbh, mmaxbh, normalization_size)
    pdf_unnormalized = _powerlaw_gaussian_unnormalized(
        m_try, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m
    )
    normalization = np.trapz(pdf_unnormalized, m_try)

    pdf_unnormalized = _powerlaw_gaussian_unnormalized(
        m, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m
    )
    pdf = pdf_unnormalized / normalization

    return pdf

@njit
def _powerlaw_gaussian_cdf(size, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m):
    """
    Compute CDF for power-law + Gaussian mass model.

    Parameters
    ----------
    size : ``int``
        Number of grid points.
    mminbh : ``float``
        Minimum BH mass.
    mmaxbh : ``float``
        Maximum BH mass.
    alpha : ``float``
        Power-law spectral index.
    mu_g : ``float``
        Mean of the Gaussian peak.
    sigma_g : ``float``
        Standard deviation of the Gaussian peak.
    lambda_peak : ``float``
        Fraction in Gaussian component.
    delta_m : ``float``
        Smoothing width.

    Returns
    -------
    cdf_values : ``numpy.ndarray``
        Normalized CDF values.
    """
    m_try = np.linspace(mminbh, mmaxbh, size)
    pdf_unnormalized = _powerlaw_gaussian_unnormalized(
        m_try, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m
    )
    cdf_values = _cumulative_trapezoid(y=pdf_unnormalized, x=m_try, dx=1.0, initial=0.0)
    normalization = cdf_values[size-1]
    cdf_values /= normalization

    return cdf_values

@njit
def _sample_powerlaw_gaussian(size, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m, normalization_size=1000):
    """
    Generate samples from the power-law + Gaussian mass model.

    Parameters
    ----------
    size : ``int``
        Number of samples to draw.
    mminbh : ``float``
        Minimum BH mass (Msun).
    mmaxbh : ``float``
        Maximum BH mass (Msun).
    alpha : ``float``
        Power-law spectral index.
    mu_g : ``float``
        Mean of the Gaussian peak (Msun).
    sigma_g : ``float``
        Standard deviation of the Gaussian peak (Msun).
    lambda_peak : ``float``
        Fraction in Gaussian component (0-1).
    delta_m : ``float``
        Low-mass smoothing width (Msun).
    normalization_size : ``int``
        Grid size for CDF computation. \n
        default: 1000

    Returns
    -------
    samples : ``numpy.ndarray``
        Mass samples (Msun).
    """
    cdf_values = _powerlaw_gaussian_cdf(normalization_size, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m)

    x = np.linspace(mminbh, mmaxbh, normalization_size)
    idx = np.argwhere(cdf_values > 0)[0][0]
    cdf_values = cdf_values[idx:]
    x = x[idx:]
    samples = inverse_transform_sampler(size, cdf_values, x)

    return samples

@njit
def binary_masses_BBH_powerlaw_gaussian_rvs(size, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m, beta, normalization_size=1000):
    """
    Generate BBH mass samples from power-law + Gaussian model with mass ratio.

    Parameters
    ----------
    size : ``int``
        Number of samples to draw.
    mminbh : ``float``
        Minimum BH mass (Msun).
    mmaxbh : ``float``
        Maximum BH mass (Msun).
    alpha : ``float``
        Power-law spectral index for m1.
    mu_g : ``float``
        Mean of the Gaussian peak (Msun).
    sigma_g : ``float``
        Standard deviation of the Gaussian peak (Msun).
    lambda_peak : ``float``
        Fraction in Gaussian component (0-1).
    delta_m : ``float``
        Low-mass smoothing width (Msun).
    beta : ``float``
        Power-law index for mass ratio distribution.
    normalization_size : ``int``
        Grid size for CDF computation. \n
        default: 1000

    Returns
    -------
    m1 : ``numpy.ndarray``
        Primary mass samples (Msun).
    m2 : ``numpy.ndarray``
        Secondary mass samples (Msun).
    """
    m1 = _sample_powerlaw_gaussian(size, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m, normalization_size)

    q = np.zeros(size)
    for i in range(size):
        q[i] = _sample_mass_ratio(m1[i], mminbh, beta, delta_m)

    m2 = m1 * q

    return m1, m2

@njit
def _sample_mass_ratio(m1, mminbh, beta, delta_m):
    """
    Sample mass ratio using rejection sampling with smoothing.

    Parameters
    ----------
    m1 : ``float``
        Primary mass.
    mminbh : ``float``
        Minimum BH mass.
    beta : ``float``
        Power-law index for mass ratio.
    delta_m : ``float``
        Smoothing width.

    Returns
    -------
    q : ``float``
        Sampled mass ratio.
    """
    qmin = mminbh / m1
    pow_beta = beta + 1.0

    while True:
        u_q = np.random.rand()
        q = (u_q * (1.0**pow_beta - qmin**pow_beta) + qmin**pow_beta)**(1.0 / pow_beta)
        m2 = m1 * q
        s_m2 = _smoothing_S(np.array([m2]), mminbh, delta_m)[0]
        u_smooth = np.random.rand()
        if u_smooth < s_m2:
            break

    return q

@njit
def _powerlaw_gaussian_unnormalized(m, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m):
    """
    Compute unnormalized PDF for power-law + Gaussian model.

    Parameters
    ----------
    m : ``numpy.ndarray``
        Mass values.
    mminbh : ``float``
        Minimum BH mass.
    mmaxbh : ``float``
        Maximum BH mass.
    alpha : ``float``
        Power-law spectral index.
    mu_g : ``float``
        Mean of Gaussian peak.
    sigma_g : ``float``
        Standard deviation of Gaussian peak.
    lambda_peak : ``float``
        Fraction in Gaussian component.
    delta_m : ``float``
        Smoothing width.

    Returns
    -------
    pdf_unnormalized : ``numpy.ndarray``
        Unnormalized PDF values.
    """
    pdf_unnormalized = ((1-lambda_peak)*_powerlaw_B(m, alpha, mminbh, mmaxbh) + (lambda_peak * _gaussian_G(m, mu_g, sigma_g)))* _smoothing_S(m, mminbh, delta_m)

    return pdf_unnormalized

def available_prior_list():
    """
    Returns a list of available priors.
    """
    return [
        'merger_rate_density_bbh_oguri2018_function',
        'merger_rate_density_bbh_popIII_ken2022_function',
        'merger_rate_density_madau_dickinson2014_function',
        'merger_rate_density_madau_dickinson_belczynski_ng_function',
        'merger_rate_density_bbh_primordial_ken2022_function',
        'sfr_madau_fragos2017_with_bbh_td',
        'sfr_madau_dickinson2014_with_bbh_td',
        'sfr_madau_fragos2017_with_bns_td',
        'sfr_madau_dickinson2014_with_bns_td',
        'sfr_madau_fragos2017',
        'sfr_madau_dickinson2014',
        'binary_masses_BBH_popIII_lognormal_rvs',
        'binary_masses_BBH_primordial_lognormal_rvs',
        'binary_masses_BNS_bimodal_rvs',
        'binary_masses_NSBH_broken_powerlaw_rvs',
        'binary_masses_BBH_powerlaw_gaussian_rvs',
    ]

# ------------------------
# Other utility functions
# ------------------------
@njit
def _cumulative_trapezoid(y, x=None, dx=1.0, initial=0.0):
    """
    Compute the cumulative integral using the trapezoidal rule.

    Parameters
    ----------
    y : ``numpy.ndarray``
        Function values to integrate.
    x : ``numpy.ndarray`` or ``None``
        x-coordinates. If None, uses evenly spaced points with spacing dx. \n
        default: None
    dx : ``float``
        Spacing between x-coordinates if x is None. \n
        default: 1.0
    initial : ``float``
        Initial value for the cumulative sum. \n
        default: 0.0

    Returns
    -------
    cumsum : ``numpy.ndarray``
        Cumulative integral values.
    """
    if x is None:
        x = np.arange(len(y)) * dx

    cumsum = np.zeros_like(y)
    cumsum[0] = initial
    for i in range(1, len(y)):
        cumsum[i] = cumsum[i - 1] + (y[i - 1] + y[i]) * (x[i] - x[i - 1]) / 2.0

    return cumsum



        

