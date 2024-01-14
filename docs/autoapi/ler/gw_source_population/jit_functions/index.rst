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

   ler.gw_source_population.jit_functions.merger_rate_density_bbh_popI_II_oguri2018
   ler.gw_source_population.jit_functions.merger_rate_density_bbh_popIII_ken2022
   ler.gw_source_population.jit_functions.star_formation_rate_madau_dickinson2014
   ler.gw_source_population.jit_functions.merger_rate_density_bbh_primordial_ken2022
   ler.gw_source_population.jit_functions.lognormal_distribution_2D
   ler.gw_source_population.jit_functions.inverse_transform_sampler_m1m2



.. py:function:: merger_rate_density_bbh_popI_II_oguri2018(zs, R0=23.9 * 1e-09, b2=1.6, b3=2.0, b4=30)

   
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
           default: 2.0

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

.. py:function:: star_formation_rate_madau_dickinson2014(zs, af=2.7, bf=5.6, cf=2.9)

   
   Function to compute star formation rate as given in Eqn. 15 Madau & Dickinson (2014). The output is in detector frame and is unnormalized.


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

   >>> from ler.gw_source_population import star_formation_rate_madau_dickinson2014
   >>> rate_density = star_formation_rate_madau_dickinson2014(zs=0.1)



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

