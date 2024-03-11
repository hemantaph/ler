# -*- coding: utf-8 -*-
"""
This module contains various functions use for simulating GW source population.
"""

import numpy as np
from numba import njit, jit
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

from ler.utils import inverse_transform_sampler

# import pickle
# # call the interpolator
# try:
#     with open('./mp_interpolator/zs_inv_cdf.pickle', 'rb') as input_:
#         print('Loading the interpolator...')
#         zs_inv_cdf = pickle.load(input_)
# except:
#     print("can't Loading the interpolator...")
#     pass

@njit
def sample_source_redshift(size, zs_inv_cdf=None):

    u = np.random.uniform(0, 1, size=size)
    x = zs_inv_cdf[0]  # cdf values
    y = zs_inv_cdf[1]  # redshift values
    return np.interp(u, x, y)

@njit
def merger_rate_density_bbh_popI_II_oguri2018(
        zs, R0=23.9 * 1e-9, b2=1.6, b3=2.0, b4=30,
    ):
    """
    Function to compute the merger rate density (PopI/PopII). Reference: Oguri et al. (2018). The output is in detector frame and is unnormalized.

    Parameters
    ----------
    zs : `float` or `numpy.ndarray` (nD array of floats)
        Source redshifts
    R0 : `float`
        local merger rate density at low redshift
        default: 23.9*1e-9 Mpc^-3 yr^-1
    b2 : `float`
        Fitting paramters
        default: 1.6
    b3 : `float`
        Fitting paramters
        default: 2.0
    b4 : `float`
        Fitting paramters
        default: 30

    Returns
    ----------
    rate_density : `float` or `numpy.ndarray` (nD array of floats)
        merger rate density

    Examples
    ----------
    >>> from ler.gw_source_population import merger_rate_density_bbh_popI_II_oguri2018
    >>> rate_density = merger_rate_density_bbh_popI_II_oguri2018(zs=0.1)
    """
    
    # rate_density
    return R0 * (b4 + 1) * np.exp(b2 * zs) / (b4 + np.exp(b3 * zs))

@njit
def merger_rate_density_bbh_popIII_ken2022(zs, n0=19.2 * 1e-9, aIII=0.66, bIII=0.3, zIII=11.6,
    ):
    """
    Function to compute the unnormalized merger rate density (PopIII). Reference: Ng et al. 2022. The output is in detector frame and is unnormalized.

    Parameters
    ----------
    zs : `float` or `numpy.ndarray` (nD array of floats)
        Source redshifts
    n0 : `float`
        normalization constant
        default: 19.2*1e-9
    aIII : `float`
        Fitting paramters
        default: 0.66
    bIII : `float`
        Fitting paramters
        default: 0.3
    zIII : `float`
        Fitting paramters
        default: 11.6

    Returns
    ----------
    rate_density : `float` or `numpy.ndarray` (nD array of floats)
        merger rate density

    Examples
    ----------
    >>> from ler.gw_source_population import merger_rate_density_bbh_popIII_ken2022
    >>> rate_density = merger_rate_density_bbh_popIII_ken2022(zs=0.1)
    """

    # rate density
    return (
        n0
        * np.exp(aIII * (zs - zIII))
        / (bIII + aIII * np.exp((aIII + bIII) * (zs - zIII)))
    )

@njit
def star_formation_rate_madau_dickinson2014(
        zs, af=2.7, bf=5.6, cf=2.9,
    ):
    """
    Function to compute star formation rate as given in Eqn. 15 Madau & Dickinson (2014). The output is in detector frame and is unnormalized.

    Parameters
    ----------
    zs : `float` or `numpy.ndarray` (nD array of floats)
        Source redshifts
    af : `float`
        Fitting paramters
        default: 2.7
    bf : `float`
        Fitting paramters
        default: 5.6
    cf : `float`
        Fitting paramters
        default: 2.9
        
    Returns
    ----------
    rate_density : `float` or `numpy.ndarray` (nD array of floats)
        merger rate density

    Examples
    ----------
    >>> from ler.gw_source_population import star_formation_rate_madau_dickinson2014
    >>> rate_density = star_formation_rate_madau_dickinson2014(zs=0.1)
    """

    # rate density
    return 0.015 * (1 + zs) ** af / (1 + ((1 + zs) / cf) ** bf)


@jit
def merger_rate_density_bbh_primordial_ken2022(
        zs, cosmology=cosmo, n0=0.044 * 1e-9, t0=13.786885302009708
    ):
        """
        Function to compute the merger rate density (Primordial). Reference: Ng et al. 2022. The output is in detector frame and is unnormalized.

        Parameters
        ----------
        zs : `float`
            Source redshifts
        n0 : `float`
            normalization constant
            default: 0.044*1e-9
        t0 : `float`
            Present age of the Universe in Gyr
            default: 13.786885302009708
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. param = dict(t0=13.786885302009708)

        Returns
        ----------
        rate_density : `float`
            merger rate density

        Examples
        ----------
        """

        # rate density
        rate_density = n0 * (cosmology.age(z=zs).value / t0) ** (-34 / 37)

        return rate_density

@njit
def lognormal_distribution_2D(size, m_min=1.0, m_max=100.0, Mc=20.0, sigma=0.3, chunk_size=10000):
    """
    Function to sample from a lognormal distribution in 2D space. Reference: Ng et al. 2022. This a helper function for popIII BBH and primordial BBH merger rate density distribution functions.

    Parameters
    ----------
    size : `int`
        Number of samples to draw
    m_min : `float`
        Minimum mass
        default: 1.0
    m_max : `float`
        Maximum mass
        default: 100.0
    Mc : `float`
        Mass scale
        default: 20.0
    sigma : `float`
        width of the distribution
        default: 0.3
    chunk_size : `int`
        Number of samples to draw in each chunk
        default: 10000

    Returns
    ----------
    m1_sample : `numpy.ndarray` (1D array of floats)
        Mass of the primary
    m2_sample : `numpy.ndarray` (1D array of floats)
        Mass of the secondary

    Examples
    ----------
    >>> from ler.gw_source_population import lognormal_distribution_2D
    >>> m1_sample, m2_sample = lognormal_distribution_2D(size=1000)
    """
       
    # mass function. Eqn. 1 of Ng et al. 2022
    psi = lambda m: np.exp(-np.log(m / Mc) ** 2 / (2 * sigma**2)) / (
        np.sqrt(2 * np.pi) * sigma * m
    )
    # probability density function
    # Eqn. 4 of Ng et al. 2022
    pdf = (
        lambda m1, m2: (m1 + m2) ** (36 / 37)
        * (m1 * m2) ** (32 / 37)
        * psi(m1)
        * psi(m2)
    )

    # rejection sampling
    m1 = np.random.uniform(m_min, m_max, chunk_size)
    m2 = np.random.uniform(m_min, m_max, chunk_size)
    z = pdf(m1, m2)
    zmax = np.max(z)

    # Rejection sample in chunks
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
            
    # swap the mass if m1 < m2
    idx = m1_sample < m2_sample
    m1_sample[idx], m2_sample[idx] = m2_sample[idx], m1_sample[idx]
    return m1_sample, m2_sample

@njit
def inverse_transform_sampler_m1m2(size, inv_cdf, x):
    """
    Function to sample from a distribution using inverse transform sampling. This is a helper function BNS Alsing mass distribution function.

    Parameters
    ----------
    size : `int`
        Number of samples to draw
    inv_cdf : `numpy.ndarray` (1D array of floats)
        Inverse cumulative distribution function
    x : `numpy.ndarray` (1D array of floats)
        array of mass values for which the inverse cumulative distribution function is computed

    Returns
    ----------
    m1 : `numpy.ndarray` (1D array of floats)
        Mass of the primary
    m2 : `numpy.ndarray` (1D array of floats)
        Mass of the secondary

    Examples
    ----------
    >>> from ler.gw_source_population import inverse_transform_sampler_m1m2
    >>> m1, m2 = inverse_transform_sampler_m1m2(size=1000, inv_cdf=inv_cdf, x=x)
    """
    
    m1 = inverse_transform_sampler(size, inv_cdf, x)
    m2 = inverse_transform_sampler(size, inv_cdf, x)
    # swap m1 and m2 if m1 < m2
    idx = m1 < m2
    m1[idx], m2[idx] = m2[idx], m1[idx]
    
    return m1, m2
