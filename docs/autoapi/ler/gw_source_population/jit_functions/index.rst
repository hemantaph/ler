:py:mod:`ler.gw_source_population.jit_functions`
================================================

.. py:module:: ler.gw_source_population.jit_functions

.. autoapi-nested-parse::

   Module for JIT-compiled functions used in GW source population simulations.

   This module provides Numba JIT-compiled functions for efficient sampling and
   computation of merger rate densities, star formation rates, and mass distributions for various gravitational wave source populations including BBH, BNS, and NSBH.

   Key Features:

   - Merger rate density functions for PopI/II, PopIII, and Primordial BBH

   - Star formation rate models (Madau & Dickinson 2014, Madau & Fragos 2017)

   - Mass distribution samplers (power-law Gaussian, broken power-law, bimodal)

   - JIT-compiled for high performance with Numba


   Copyright (C) 2026 Hemanta Ph. Distributed under MIT License.

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

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
   ler.gw_source_population.jit_functions.bns_bimodal_pdf
   ler.gw_source_population.jit_functions.sample_broken_powerlaw
   ler.gw_source_population.jit_functions.sample_broken_powerlaw_nsbh_masses
   ler.gw_source_population.jit_functions.broken_powerlaw_pdf
   ler.gw_source_population.jit_functions.powerlaw_gaussian_pdf
   ler.gw_source_population.jit_functions.sample_powerlaw_gaussian
   ler.gw_source_population.jit_functions.sample_powerlaw_gaussian_source_bbh_masses



.. py:function:: merger_rate_density_bbh_popI_II_oguri2018(zs, R0=23.9 * 1e-09, b2=1.6, b3=2.1, b4=30)

   
   Compute the merger rate density for PopI/II BBH.

   Reference: Oguri et al. (2018). The output is in detector frame and is
   unnormalized.

   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **R0** : ``float``
           Local merger rate density at low redshift (Mpc^-3 yr^-1).

           default: 23.9e-9

       **b2** : ``float``
           Fitting parameter.

           default: 1.6

       **b3** : ``float``
           Fitting parameter.

           default: 2.1

       **b4** : ``float``
           Fitting parameter.

           default: 30

   :Returns:

       **rate_density** : ``float`` or ``numpy.ndarray``
           Merger rate density.










   .. rubric:: Examples

   >>> from ler.gw_source_population import merger_rate_density_bbh_popI_II_oguri2018
   >>> rate_density = merger_rate_density_bbh_popI_II_oguri2018(zs=0.1)



   ..
       !! processed by numpydoc !!

.. py:function:: merger_rate_density_bbh_popIII_ken2022(zs, n0=19.2 * 1e-09, aIII=0.66, bIII=0.3, zIII=11.6)

   
   Compute the unnormalized merger rate density for PopIII BBH.

   Reference: Ng et al. (2022). The output is in detector frame and is
   unnormalized.

   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **n0** : ``float``
           Normalization constant.

           default: 19.2e-9

       **aIII** : ``float``
           Fitting parameter.

           default: 0.66

       **bIII** : ``float``
           Fitting parameter.

           default: 0.3

       **zIII** : ``float``
           Characteristic redshift.

           default: 11.6

   :Returns:

       **rate_density** : ``float`` or ``numpy.ndarray``
           Merger rate density.










   .. rubric:: Examples

   >>> from ler.gw_source_population import merger_rate_density_bbh_popIII_ken2022
   >>> rate_density = merger_rate_density_bbh_popIII_ken2022(zs=0.1)



   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_fragos2017_with_bbh_td(zs, R0=23.9 * 1e-09)

   
   Compute star formation rate with BBH time delay using Madau & Fragos (2017).


   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **R0** : ``float``
           Local merger rate density (Mpc^-3 yr^-1).

           default: 23.9e-9

   :Returns:

       **SFR** : ``float`` or ``numpy.ndarray``
           Star formation rate (Mpc^-3 yr^-1).













   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_dickinson2014_with_bbh_td(zs, R0=23.9 * 1e-09)

   
   Compute star formation rate with BBH time delay using Madau & Dickinson (2014).


   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **R0** : ``float``
           Local merger rate density (Mpc^-3 yr^-1).

           default: 23.9e-9

   :Returns:

       **SFR** : ``float`` or ``numpy.ndarray``
           Star formation rate (Mpc^-3 yr^-1).













   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_fragos2017_with_bns_td(zs, R0=105.5 * 1e-09)

   
   Compute star formation rate with BNS time delay using Madau & Fragos (2017).


   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **R0** : ``float``
           Local merger rate density (Mpc^-3 yr^-1).

           default: 105.5e-9

   :Returns:

       **SFR** : ``float`` or ``numpy.ndarray``
           Star formation rate (Mpc^-3 yr^-1).













   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_dickinson2014_with_bns_td(zs, R0=105.5 * 1e-09)

   
   Compute star formation rate with BNS time delay using Madau & Dickinson (2014).


   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **R0** : ``float``
           Local merger rate density (Mpc^-3 yr^-1).

           default: 105.5e-9

   :Returns:

       **SFR** : ``float`` or ``numpy.ndarray``
           Star formation rate (Mpc^-3 yr^-1).













   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_fragos2017(zs, a=0.01, b=2.6, c=3.2, d=6.2)

   
   Compute star formation rate using Madau & Fragos (2017) model.

   Reference: https://arxiv.org/pdf/1606.07887.pdf

   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **a** : ``float``
           Normalization parameter.

           default: 0.01

       **b** : ``float``
           Low-redshift power-law slope.

           default: 2.6

       **c** : ``float``
           Turnover redshift parameter.

           default: 3.2

       **d** : ``float``
           High-redshift power-law slope.

           default: 6.2

   :Returns:

       **SFR** : ``float`` or ``numpy.ndarray``
           Star formation rate (Msun yr^-1 Mpc^-3).













   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_dickinson2014(zs, a=0.015, b=2.7, c=2.9, d=5.6)

   
   Compute star formation rate using Madau & Dickinson (2014) model.

   Reference: Eqn. 15 of https://arxiv.org/pdf/1403.0007

   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **a** : ``float``
           Normalization parameter.

           default: 0.015

       **b** : ``float``
           Low-redshift power-law slope.

           default: 2.7

       **c** : ``float``
           Turnover redshift parameter.

           default: 2.9

       **d** : ``float``
           High-redshift power-law slope.

           default: 5.6

   :Returns:

       **SFR** : ``float`` or ``numpy.ndarray``
           Star formation rate (Msun yr^-1 Mpc^-3).










   .. rubric:: Examples

   >>> from ler.gw_source_population import sfr_madau_dickinson2014
   >>> sfr = sfr_madau_dickinson2014(zs=0.1)



   ..
       !! processed by numpydoc !!

.. py:function:: merger_rate_density_bbh_primordial_ken2022(zs, cosmology=None, n0=0.044 * 1e-09, t0=13.786885302009708)

   
   Compute the merger rate density for Primordial BBH.

   Reference: Ng et al. (2022). The output is in detector frame and is
   unnormalized.

   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **cosmology** : ``astropy.cosmology`` or ``None``
           Cosmology object for age calculations.

           default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

       **n0** : ``float``
           Normalization constant.

           default: 0.044e-9

       **t0** : ``float``
           Present age of the Universe (Gyr).

           default: 13.786885302009708

   :Returns:

       **rate_density** : ``float`` or ``numpy.ndarray``
           Merger rate density.










   .. rubric:: Examples

   >>> from ler.gw_source_population import merger_rate_density_bbh_primordial_ken2022
   >>> rate_density = merger_rate_density_bbh_primordial_ken2022(zs=0.1)



   ..
       !! processed by numpydoc !!

.. py:function:: lognormal_distribution_2D(size, m_min=1.0, m_max=100.0, Mc=20.0, sigma=0.3, chunk_size=10000)

   
   Sample from a lognormal distribution in 2D mass space.

   Reference: Ng et al. (2022). This is a helper function for PopIII BBH
   and primordial BBH merger rate density distribution functions.

   :Parameters:

       **size** : ``int``
           Number of samples to draw.

       **m_min** : ``float``
           Minimum mass (Msun).

           default: 1.0

       **m_max** : ``float``
           Maximum mass (Msun).

           default: 100.0

       **Mc** : ``float``
           Characteristic mass scale (Msun).

           default: 20.0

       **sigma** : ``float``
           Width of the distribution.

           default: 0.3

       **chunk_size** : ``int``
           Number of samples per rejection sampling chunk.

           default: 10000

   :Returns:

       **m1_sample** : ``numpy.ndarray``
           Primary mass samples (Msun).

       **m2_sample** : ``numpy.ndarray``
           Secondary mass samples (Msun).










   .. rubric:: Examples

   >>> from ler.gw_source_population import lognormal_distribution_2D
   >>> m1, m2 = lognormal_distribution_2D(size=1000)



   ..
       !! processed by numpydoc !!

.. py:function:: inverse_transform_sampler_m1m2(size, inv_cdf, x)

   
   Sample m1 and m2 using inverse transform sampling for BNS.

   This is a helper function for the BNS Alsing mass distribution function.

   :Parameters:

       **size** : ``int``
           Number of samples to draw.

       **inv_cdf** : ``numpy.ndarray``
           Cumulative distribution function values.

       **x** : ``numpy.ndarray``
           Mass values corresponding to the CDF.

   :Returns:

       **m1** : ``numpy.ndarray``
           Primary mass samples (Msun).

       **m2** : ``numpy.ndarray``
           Secondary mass samples (Msun).










   .. rubric:: Examples

   >>> from ler.gw_source_population import inverse_transform_sampler_m1m2
   >>> m1, m2 = inverse_transform_sampler_m1m2(size=1000, inv_cdf=cdf, x=mass_arr)



   ..
       !! processed by numpydoc !!

.. py:function:: bns_bimodal_pdf(m, w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1.0, mmax=2.3)

   
   Compute the bimodal Gaussian PDF for BNS mass distribution.


   :Parameters:

       **m** : ``float`` or ``numpy.ndarray``
           Mass values (Msun).

       **w** : ``float``
           Weight of the left (low-mass) peak.

           default: 0.643

       **muL** : ``float``
           Mean of the left peak (Msun).

           default: 1.352

       **sigmaL** : ``float``
           Standard deviation of the left peak (Msun).

           default: 0.08

       **muR** : ``float``
           Mean of the right peak (Msun).

           default: 1.88

       **sigmaR** : ``float``
           Standard deviation of the right peak (Msun).

           default: 0.3

       **mmin** : ``float``
           Minimum mass (Msun).

           default: 1.0

       **mmax** : ``float``
           Maximum mass (Msun).

           default: 2.3

   :Returns:

       **pdf** : ``float`` or ``numpy.ndarray``
           Probability density values.













   ..
       !! processed by numpydoc !!

.. py:function:: sample_broken_powerlaw(size=1000, mminbh=26.0, mmaxbh=125.0, alpha_1=6.75, alpha_2=0.0, b=0.5, delta_m=5.0, normalization_size=1000)

   
   Generate samples from the broken power-law mass distribution.


   :Parameters:

       **size** : ``int``
           Number of samples to draw.

           default: 1000

       **mminbh** : ``float``
           Minimum BH mass (Msun).

           default: 26.0

       **mmaxbh** : ``float``
           Maximum BH mass (Msun).

           default: 125.0

       **alpha_1** : ``float``
           Power-law index below the break.

           default: 6.75

       **alpha_2** : ``float``
           Power-law index above the break.

           default: 0.0

       **b** : ``float``
           Break location parameter (0-1).

           default: 0.5

       **delta_m** : ``float``
           Smoothing width (Msun).

           default: 5.0

       **normalization_size** : ``int``
           Grid size for CDF computation.

           default: 1000

   :Returns:

       **samples** : ``numpy.ndarray``
           Mass samples (Msun).













   ..
       !! processed by numpydoc !!

.. py:function:: sample_broken_powerlaw_nsbh_masses(size=1000, mminbh=26.0, mmaxbh=125.0, alpha_1=6.75, alpha_2=0.0, b=0.5, delta_m=5.0, mminns=1.0, mmaxns=3.0, alphans=0.0, normalization_size=1000)

   
   Generate NSBH mass samples from broken power-law (BH) and power-law (NS).


   :Parameters:

       **size** : ``int``
           Number of samples to draw.

           default: 1000

       **mminbh** : ``float``
           Minimum BH mass (Msun).

           default: 26.0

       **mmaxbh** : ``float``
           Maximum BH mass (Msun).

           default: 125.0

       **alpha_1** : ``float``
           BH power-law index below break.

           default: 6.75

       **alpha_2** : ``float``
           BH power-law index above break.

           default: 0.0

       **b** : ``float``
           Break location parameter (0-1).

           default: 0.5

       **delta_m** : ``float``
           Smoothing width (Msun).

           default: 5.0

       **mminns** : ``float``
           Minimum NS mass (Msun).

           default: 1.0

       **mmaxns** : ``float``
           Maximum NS mass (Msun).

           default: 3.0

       **alphans** : ``float``
           NS power-law index.

           default: 0.0

       **normalization_size** : ``int``
           Grid size for CDF computation.

           default: 1000

   :Returns:

       **m1_samples** : ``numpy.ndarray``
           BH mass samples (Msun).

       **m2_samples** : ``numpy.ndarray``
           NS mass samples (Msun).













   ..
       !! processed by numpydoc !!

.. py:function:: broken_powerlaw_pdf(m, mminbh=26.0, mmaxbh=125.0, alpha_1=6.75, alpha_2=0.0, b=0.5, delta_m=5.0, normalization_size=1000)

   
   Compute the normalized PDF for broken power-law mass distribution.


   :Parameters:

       **m** : ``numpy.ndarray``
           Mass values to evaluate (Msun).

       **mminbh** : ``float``
           Minimum BH mass (Msun).

           default: 26.0

       **mmaxbh** : ``float``
           Maximum BH mass (Msun).

           default: 125.0

       **alpha_1** : ``float``
           Power-law index below the break.

           default: 6.75

       **alpha_2** : ``float``
           Power-law index above the break.

           default: 0.0

       **b** : ``float``
           Break location parameter (0-1).

           default: 0.5

       **delta_m** : ``float``
           Smoothing width (Msun).

           default: 5.0

       **normalization_size** : ``int``
           Grid size for normalization.

           default: 1000

   :Returns:

       **pdf** : ``numpy.ndarray``
           Normalized probability density values.













   ..
       !! processed by numpydoc !!

.. py:function:: powerlaw_gaussian_pdf(m, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m, normalization_size=1000)

   
   Compute the normalized PDF for power-law + Gaussian mass model.


   :Parameters:

       **m** : ``numpy.ndarray``
           Mass values to evaluate (Msun).

       **mminbh** : ``float``
           Minimum BH mass (Msun).

       **mmaxbh** : ``float``
           Maximum BH mass (Msun).

       **alpha** : ``float``
           Power-law spectral index.

       **mu_g** : ``float``
           Mean of the Gaussian peak (Msun).

       **sigma_g** : ``float``
           Standard deviation of the Gaussian peak (Msun).

       **lambda_peak** : ``float``
           Fraction of sources in the Gaussian component (0-1).

       **delta_m** : ``float``
           Low-mass smoothing width (Msun).

       **normalization_size** : ``int``
           Grid size for normalization.

           default: 1000

   :Returns:

       **pdf** : ``numpy.ndarray``
           Normalized probability density values.













   ..
       !! processed by numpydoc !!

.. py:function:: sample_powerlaw_gaussian(size, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m, normalization_size=1000)

   
   Generate samples from the power-law + Gaussian mass model.


   :Parameters:

       **size** : ``int``
           Number of samples to draw.

       **mminbh** : ``float``
           Minimum BH mass (Msun).

       **mmaxbh** : ``float``
           Maximum BH mass (Msun).

       **alpha** : ``float``
           Power-law spectral index.

       **mu_g** : ``float``
           Mean of the Gaussian peak (Msun).

       **sigma_g** : ``float``
           Standard deviation of the Gaussian peak (Msun).

       **lambda_peak** : ``float``
           Fraction in Gaussian component (0-1).

       **delta_m** : ``float``
           Low-mass smoothing width (Msun).

       **normalization_size** : ``int``
           Grid size for CDF computation.

           default: 1000

   :Returns:

       **samples** : ``numpy.ndarray``
           Mass samples (Msun).













   ..
       !! processed by numpydoc !!

.. py:function:: sample_powerlaw_gaussian_source_bbh_masses(size, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m, beta, normalization_size=1000)

   
   Generate BBH mass samples from power-law + Gaussian model with mass ratio.


   :Parameters:

       **size** : ``int``
           Number of samples to draw.

       **mminbh** : ``float``
           Minimum BH mass (Msun).

       **mmaxbh** : ``float``
           Maximum BH mass (Msun).

       **alpha** : ``float``
           Power-law spectral index for m1.

       **mu_g** : ``float``
           Mean of the Gaussian peak (Msun).

       **sigma_g** : ``float``
           Standard deviation of the Gaussian peak (Msun).

       **lambda_peak** : ``float``
           Fraction in Gaussian component (0-1).

       **delta_m** : ``float``
           Low-mass smoothing width (Msun).

       **beta** : ``float``
           Power-law index for mass ratio distribution.

       **normalization_size** : ``int``
           Grid size for CDF computation.

           default: 1000

   :Returns:

       **m1** : ``numpy.ndarray``
           Primary mass samples (Msun).

       **m2** : ``numpy.ndarray``
           Secondary mass samples (Msun).













   ..
       !! processed by numpydoc !!

