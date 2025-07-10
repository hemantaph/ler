:py:mod:`ler.gw_source_population.jit_functions`
================================================

.. py:module:: ler.gw_source_population.jit_functions

.. autoapi-nested-parse::

   This module contains various functions use for simulating GW source population.

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.gw_source_population.jit_functions.cumulative_trapezoid
   ler.gw_source_population.jit_functions.merger_rate_density_bbh_popI_II_oguri2018
   ler.gw_source_population.jit_functions.merger_rate_density_bbh_popIII_ken2022
   ler.gw_source_population.jit_functions.sfr_madau_fragos2017_with_bbh_td
   ler.gw_source_population.jit_functions.sfr_madau_dickinson2014_with_bbh_td
   ler.gw_source_population.jit_functions.sfr_madau_fragos2017_with_bns_td
   ler.gw_source_population.jit_functions.sfr_madau_dickinson2014_with_bns_td
   ler.gw_source_population.jit_functions.sfr_madau_fragos2017
   ler.gw_source_population.jit_functions.sfr_madau_dickinson2014
   ler.gw_source_population.jit_functions.merger_rate_density_bbh_primordial_ken2022
   ler.gw_source_population.jit_functions.lognormal_distribution_2D
   ler.gw_source_population.jit_functions.inverse_transform_sampler_m1m2
   ler.gw_source_population.jit_functions.powerlaw_with_smoothing
   ler.gw_source_population.jit_functions.inverse_transform_sampler
   ler.gw_source_population.jit_functions.sample_broken_powerlaw
   ler.gw_source_population.jit_functions.sample_broken_powerlaw_nsbh_masses
   ler.gw_source_population.jit_functions.broken_powerlaw_pdf
   ler.gw_source_population.jit_functions.broken_powerlaw_unormalized
   ler.gw_source_population.jit_functions.powerlaw_B
   ler.gw_source_population.jit_functions.gaussian_G
   ler.gw_source_population.jit_functions.powerlaw_gaussian_pdf
   ler.gw_source_population.jit_functions.powerlaw_gaussian_cdf
   ler.gw_source_population.jit_functions.sample_powerlaw_gaussian
   ler.gw_source_population.jit_functions.sample_powerlaw_gaussian_source_bbh_masses
   ler.gw_source_population.jit_functions.powerlaw_gaussian_unnormalized



.. py:function:: cumulative_trapezoid(y, x=None, dx=1.0, initial=0.0)

   
   Compute the cumulative integral of a function using the trapezoidal rule.
















   ..
       !! processed by numpydoc !!

.. py:function:: merger_rate_density_bbh_popI_II_oguri2018(zs, R0=23.9 * 1e-09, b2=1.6, b3=2.1, b4=30)

   
   Function to compute the merger rate density (PopI/PopII). Reference: Oguri et al. (2018). The output is in detector frame and is unnormalized.


   :Parameters:

       **zs** : `float` or `numpy.ndarray` (nD array of floats)
           Source redshifts

       **R0** : `float`
           local merger rate density at low redshift
           default: 23.9*1e-9 Mpc^-3 yr^-1

       **b2** : `float`
           Fitting paramters
           default: 1.6

       **b3** : `float`
           Fitting paramters
           default: 2.1

       **b4** : `float`
           Fitting paramters
           default: 30

   :Returns:

       **rate_density** : `float` or `numpy.ndarray` (nD array of floats)
           merger rate density










   .. rubric:: Examples

   >>> from ler.gw_source_population import merger_rate_density_bbh_popI_II_oguri2018
   >>> rate_density = merger_rate_density_bbh_popI_II_oguri2018(zs=0.1)



   ..
       !! processed by numpydoc !!

.. py:function:: merger_rate_density_bbh_popIII_ken2022(zs, n0=19.2 * 1e-09, aIII=0.66, bIII=0.3, zIII=11.6)

   
   Function to compute the unnormalized merger rate density (PopIII). Reference: Ng et al. 2022. The output is in detector frame and is unnormalized.


   :Parameters:

       **zs** : `float` or `numpy.ndarray` (nD array of floats)
           Source redshifts

       **n0** : `float`
           normalization constant
           default: 19.2*1e-9

       **aIII** : `float`
           Fitting paramters
           default: 0.66

       **bIII** : `float`
           Fitting paramters
           default: 0.3

       **zIII** : `float`
           Fitting paramters
           default: 11.6

   :Returns:

       **rate_density** : `float` or `numpy.ndarray` (nD array of floats)
           merger rate density










   .. rubric:: Examples

   >>> from ler.gw_source_population import merger_rate_density_bbh_popIII_ken2022
   >>> rate_density = merger_rate_density_bbh_popIII_ken2022(zs=0.1)



   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_fragos2017_with_bbh_td(zs, R0=23.9 * 1e-09)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_dickinson2014_with_bbh_td(zs, R0=23.9 * 1e-09)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_fragos2017_with_bns_td(zs, R0=105.5 * 1e-09)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_dickinson2014_with_bns_td(zs, R0=105.5 * 1e-09)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_fragos2017(zs, a=0.01, b=2.6, c=3.2, d=6.2)

   
   https://arxiv.org/pdf/1606.07887.pdf
















   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_dickinson2014(zs, a=0.015, b=2.7, c=2.9, d=5.6)

   
   Function to compute star formation rate as given in Eqn. 15 Madau & Dickinson (2014). The output is in detector frame and is unnormalized. https://arxiv.org/pdf/1403.0007


   :Parameters:

       **zs** : `float` or `numpy.ndarray` (nD array of floats)
           Source redshifts

       **af** : `float`
           Fitting paramters
           default: 2.7

       **bf** : `float`
           Fitting paramters
           default: 5.6

       **cf** : `float`
           Fitting paramters
           default: 2.9

   :Returns:

       **rate_density** : `float` or `numpy.ndarray` (nD array of floats)
           merger rate density










   .. rubric:: Examples

   >>> from ler.gw_source_population import sfr_madau_dickinson2014
   >>> rate_density = sfr_madau_dickinson2014(zs=0.1)



   ..
       !! processed by numpydoc !!

.. py:function:: merger_rate_density_bbh_primordial_ken2022(zs, cosmology=cosmo, n0=0.044 * 1e-09, t0=13.786885302009708)

   
   Function to compute the merger rate density (Primordial). Reference: Ng et al. 2022. The output is in detector frame and is unnormalized.


   :Parameters:

       **zs** : `float`
           Source redshifts

       **n0** : `float`
           normalization constant
           default: 0.044*1e-9

       **t0** : `float`
           Present age of the Universe in Gyr
           default: 13.786885302009708

       **param** : `dict`
           Allows to pass in above parameters as dict.
           e.g. param = dict(t0=13.786885302009708)

   :Returns:

       **rate_density** : `float`
           merger rate density













   ..
       !! processed by numpydoc !!

.. py:function:: lognormal_distribution_2D(size, m_min=1.0, m_max=100.0, Mc=20.0, sigma=0.3, chunk_size=10000)

   
   Function to sample from a lognormal distribution in 2D space. Reference: Ng et al. 2022. This a helper function for popIII BBH and primordial BBH merger rate density distribution functions.


   :Parameters:

       **size** : `int`
           Number of samples to draw

       **m_min** : `float`
           Minimum mass
           default: 1.0

       **m_max** : `float`
           Maximum mass
           default: 100.0

       **Mc** : `float`
           Mass scale
           default: 20.0

       **sigma** : `float`
           width of the distribution
           default: 0.3

       **chunk_size** : `int`
           Number of samples to draw in each chunk
           default: 10000

   :Returns:

       **m1_sample** : `numpy.ndarray` (1D array of floats)
           Mass of the primary

       **m2_sample** : `numpy.ndarray` (1D array of floats)
           Mass of the secondary










   .. rubric:: Examples

   >>> from ler.gw_source_population import lognormal_distribution_2D
   >>> m1_sample, m2_sample = lognormal_distribution_2D(size=1000)



   ..
       !! processed by numpydoc !!

.. py:function:: inverse_transform_sampler_m1m2(size, inv_cdf, x)

   
   Function to sample from a distribution using inverse transform sampling. This is a helper function BNS Alsing mass distribution function.


   :Parameters:

       **size** : `int`
           Number of samples to draw

       **inv_cdf** : `numpy.ndarray` (1D array of floats)
           Inverse cumulative distribution function

       **x** : `numpy.ndarray` (1D array of floats)
           array of mass values for which the inverse cumulative distribution function is computed

   :Returns:

       **m1** : `numpy.ndarray` (1D array of floats)
           Mass of the primary

       **m2** : `numpy.ndarray` (1D array of floats)
           Mass of the secondary










   .. rubric:: Examples

   >>> from ler.gw_source_population import inverse_transform_sampler_m1m2
   >>> m1, m2 = inverse_transform_sampler_m1m2(size=1000, inv_cdf=inv_cdf, x=x)



   ..
       !! processed by numpydoc !!

.. py:function:: powerlaw_with_smoothing(m, mmin, alpha, delta_m)

   
   Power law with smoothing applied.
















   ..
       !! processed by numpydoc !!

.. py:function:: inverse_transform_sampler(size, cdf, x)

   
   Function to sample from the inverse transform method.
















   ..
       !! processed by numpydoc !!

.. py:function:: sample_broken_powerlaw(size=1000, mminbh=26.0, mmaxbh=125.0, alpha_1=6.75, alpha_2=0.0, b=0.5, delta_m=5.0, normalization_size=1000)

   
   Generates samples from the broken powerlaw distribution.
















   ..
       !! processed by numpydoc !!

.. py:function:: sample_broken_powerlaw_nsbh_masses(size=1000, mminbh=26.0, mmaxbh=125.0, alpha_1=6.75, alpha_2=0.0, b=0.5, delta_m=5.0, mminns=1.0, mmaxns=3.0, alphans=0.0, normalization_size=1000)

   
   Generates samples from the broken powerlaw distribution for NSBH masses.
















   ..
       !! processed by numpydoc !!

.. py:function:: broken_powerlaw_pdf(m, mminbh=26.0, mmaxbh=125.0, alpha_1=6.75, alpha_2=0.0, b=0.5, delta_m=5.0, normalization_size=1000)

   
   Generates samples using a Numba-jitted loop for high performance.
















   ..
       !! processed by numpydoc !!

.. py:function:: broken_powerlaw_unormalized(m, mminbh=26.0, mmaxbh=125.0, alpha_1=6.75, alpha_2=0.0, b=0.5, delta_m=5.0)

   
   Probability density function for the broken powerlaw model.
















   ..
       !! processed by numpydoc !!

.. py:function:: powerlaw_B(m, alpha, mminbh, mmaxbh)

   
   normalised power-law distribution with spectral index -alpha and cut-off mmaxbh
















   ..
       !! processed by numpydoc !!

.. py:function:: gaussian_G(m, mu_g, sigma_g)

   
   Gaussian distribution with mean mu_g and standard deviation sigma_g.
















   ..
       !! processed by numpydoc !!

.. py:function:: powerlaw_gaussian_pdf(m, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m, normalization_size=1000)

   
   Calculate the PDF for the power-law Gaussian model.
















   ..
       !! processed by numpydoc !!

.. py:function:: powerlaw_gaussian_cdf(size, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m)

   
   Sample from the power-law Gaussian model.
















   ..
       !! processed by numpydoc !!

.. py:function:: sample_powerlaw_gaussian(size, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m, normalization_size=1000)

   
   Sample from the power-law Gaussian model.
















   ..
       !! processed by numpydoc !!

.. py:function:: sample_powerlaw_gaussian_source_bbh_masses(size, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m, beta, normalization_size=1000)

   
   Sample from the power-law Gaussian model for source masses.
















   ..
       !! processed by numpydoc !!

.. py:function:: powerlaw_gaussian_unnormalized(m, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m)

   
   Calculate the unnormalized PDF for the power-law Gaussian model.
















   ..
       !! processed by numpydoc !!

