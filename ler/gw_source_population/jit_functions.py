# -*- coding: utf-8 -*-
"""
This module contains various functions use for simulating GW source population.
"""

import numpy as np
from numba import njit, jit
from scipy.interpolate import CubicSpline
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
from scipy.integrate import quad
from scipy.optimize import fsolve

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
        zs, R0=23.9 * 1e-9, b2=1.6, b3=2.1, b4=30,
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
        default: 2.1
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

# @njit
def sfr_madau_fragos2017_with_bbh_td(zs, R0=23.9 * 1e-9):
    """
    """

    rm = np.array([1.00304765, 1.00370075, 1.00449545, 1.00546251, 1.00663937, 1.00807168, 1.00981505, 1.01193727, 1.01046483, 1.01359803, 1.01741386, 1.02206193, 1.02772495, 1.03462601, 1.04303746, 1.05329142, 1.07093106, 1.08624215, 1.10489848, 1.12760683,1.15519183, 1.18858451, 1.22878158, 1.27676494, 1.33727882, 1.40335222, 1.47956936, 1.56759515, 1.6711375 , 1.79690371, 1.95410462, 2.15201042, 2.36151109, 2.66742932, 3.04354598, 3.49048755, 3.98122536, 4.42347511, 4.61710896, 4.30190679, 3.50890876, 2.37699066, 1.41830834, 0.77944771,0.40667706, 0.20463758, 0.09975143, 0.04745116])
    zs_ = np.geomspace(0.001, 10, 48)

    spline = CubicSpline(zs_, rm, extrapolate=True)
    SFR = spline(zs)*R0 # in Mpc^-3 yr^-1
    return SFR

# @njit
def sfr_madau_dickinson2014_with_bbh_td(zs, R0=23.9 * 1e-9):
    """
    """

    rm = np.array([1.00292325, 1.0035494 , 1.0043112 , 1.00523807, 1.0063658 , 1.00773798, 1.00940767, 1.01143948, 1.00997839, 1.01297699, 1.01662662, 1.02106895, 1.02647649, 1.03305927, 1.04107277, 1.05082743, 1.06802831, 1.08255749, 1.10022606, 1.12169013, 1.14772134, 1.1792097 , 1.21715406, 1.26263694, 1.32051095, 1.38462461, 1.45997648, 1.54851567, 1.65349288, 1.78046645, 1.93811129, 2.1354612 , 2.34086287, 2.63664802, 2.98892341, 3.38353439, 3.76990612, 4.03489696, 4.00806904, 3.56766897, 2.86966689, 2.01282062, 1.29696347, 0.78913584, 0.46166281, 0.26226345, 0.14509118, 0.07854392])
    zs_ = np.geomspace(0.001, 10, 48)

    spline = CubicSpline(zs_, rm, extrapolate=True)
    SFR = spline(zs)*R0 # in Mpc^-3 yr^-1
    return SFR

# @njit
def sfr_madau_fragos2017_with_bns_td(zs, R0=105.5 * 1e-9):
    """
    """

    rm = np.array([1.00309364, 1.00375139, 1.00455175, 1.00552568, 1.00671091, 1.00815339, 1.00990912, 1.01204635, 1.00757017, 1.01071962, 1.01455507, 1.01922677, 1.02491815, 1.03185311, 1.04030479, 1.05060602, 1.06970166, 1.08508957, 1.10382838, 1.12661829, 1.15427005, 1.18768774, 1.22781836, 1.27555711, 1.31791484, 1.38209039, 1.4555543 , 1.5397332 , 1.63806934, 1.75685668, 1.90448546, 2.08862044, 2.34440211, 2.63899295, 2.99729389, 3.41567274, 3.86324106, 4.24545603, 4.37018218, 4.00555831, 3.10525751, 2.06354992, 1.20906304, 0.65233811, 0.33356891, 0.16397688, 0.08024945, 0.036953])
    zs_ = np.geomspace(0.001, 10, 48)

    spline = CubicSpline(zs_, rm, extrapolate=True)
    SFR = spline(zs)*R0 # in Mpc^-3 yr^-1
    return SFR

# @njit
def sfr_madau_dickinson2014_with_bns_td(zs, R0=105.5 * 1e-9):
    """
    """

    rm = np.array([1.0029945 , 1.00362259, 1.00438674, 1.00531645, 1.00644763, 1.00782396, 1.00949865, 1.01153645, 1.00240992, 1.00539605, 1.00903013, 1.01345293, 1.01883579, 1.02538714, 1.03336026, 1.04306247, 1.05841698, 1.07283625, 1.09035966, 1.11162909, 1.1373949 , 1.16851509, 1.20594024, 1.25068092, 1.3085267 , 1.37111306, 1.44421094, 1.52948237, 1.62985636, 1.75058453, 1.90010572, 2.0870216 , 2.33573104, 2.6218286 , 2.96031682, 3.3343522 , 3.69149889, 3.92099769, 3.86227814, 3.40811745, 2.59314381, 1.79588097, 1.14260538, 0.686002  , 0.3954134 , 0.22083291, 0.11548455, 0.06064368])
    zs_ = np.geomspace(0.001, 10, 48)

    spline = CubicSpline(zs_, rm, extrapolate=True)
    SFR = spline(zs)*R0 # in Mpc^-3 yr^-1
    return SFR


@njit
def sfr_madau_fragos2017(zs, a=0.01, b=2.6, c=3.2, d=6.2):
    """
    https://arxiv.org/pdf/1606.07887.pdf
    """

    return a * (1+zs)**b / (1 + ((1+zs)/c)**d) # [Msun yr-1 Mpc-3]

@njit
def sfr_madau_dickinson2014(
        zs, a=0.015, b=2.7, c=2.9, d=5.6,
    ):
    """
    Function to compute star formation rate as given in Eqn. 15 Madau & Dickinson (2014). The output is in detector frame and is unnormalized. https://arxiv.org/pdf/1403.0007

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
    >>> from ler.gw_source_population import sfr_madau_dickinson2014
    >>> rate_density = sfr_madau_dickinson2014(zs=0.1)
    """

    # rate density
    return a * (1 + zs) ** b / (1 + ((1 + zs) / c) ** d) # [Msun yr-1 Mpc-3]


# @jit
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

@njit
def erf(x):
    # Constants for the approximation
    p = 0.3275911
    a1 = 0.254829592
    a2 = -0.284496736
    a3 = 1.421413741
    a4 = -1.453152027
    a5 = 1.061405429
    
    # Save the sign of x
    sign = np.sign(x)
    x = abs(x)
    
    # A&S formula 7.1.26 given in Handbook of Mathematical Functions
    t = 1.0 / (1.0 + p * x)
    y = 1.0 - (((((a5 * t + a4) * t) + a3) * t + a2) * t + a1) * t * np.exp(-x * x)
    
    return sign * y

@njit
def compute_normalization_factor(mu, sigma, mmin, mmax):
    part1 = (mmax - mu) / (np.sqrt(2) * sigma)
    part2 = (mmin - mu) / (np.sqrt(2) * sigma)
    N = np.sqrt(2 * np.pi) * sigma * (0.5 * (erf(part1) - erf(part2)))
    return N

@njit
def bns_bimodal_pdf(m, w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1.0, mmax=2.3,):

    # left peak
    pdf_unnormL = np.exp(-((m - muL) ** 2) / (2 * sigmaL**2))
    normL = compute_normalization_factor(muL, sigmaL, mmin, mmax)  # normalization constant
    # right peak
    pdf_unnormR = np.exp(-((m - muR) ** 2) / (2 * sigmaR**2))
    normR = compute_normalization_factor(muR, sigmaR, mmin, mmax)
    # total pdf
    pdf = w * pdf_unnormL / normL + (1 - w) * pdf_unnormR / normR

    return pdf


# def bns_bimodal_pdf_scipy(m, w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1.0, mmax=2.3,):

#     mass_arr = np.linspace(mmin, mmax, 1000)
#     # left and right peak
#     pdf_unnormL = lambda m: np.exp(-((m - muL) ** 2) / (2 * sigmaL**2))
#     normL = quad(pdf_unnormL, mmin, mmax)[0]  # normalization constant
#     pdf_unnormR = lambda m: np.exp(-((m - muR) ** 2) / (2 * sigmaR**2))
#     normR = quad(pdf_unnormR, mmin, mmax)[0]  # normalization constant
#     # total pdf
#     pdf = w * pdf_unnormL(m) / normL + (1 - w) * pdf_unnormR(m) / normR

#     return pdf
