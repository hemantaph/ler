:py:mod:`ler.gw_source_population.prior_functions`
==================================================

.. py:module:: ler.gw_source_population.prior_functions

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

   ler.gw_source_population.prior_functions.merger_rate_density_bbh_oguri2018_function
   ler.gw_source_population.prior_functions.merger_rate_density_bbh_popIII_ken2022_function
   ler.gw_source_population.prior_functions.merger_rate_density_madau_dickinson2014_function
   ler.gw_source_population.prior_functions.merger_rate_density_madau_dickinson_belczynski_ng_function
   ler.gw_source_population.prior_functions.merger_rate_density_bbh_primordial_ken2022_function
   ler.gw_source_population.prior_functions.sfr_madau_fragos2017_with_bbh_td
   ler.gw_source_population.prior_functions.sfr_madau_dickinson2014_with_bbh_td
   ler.gw_source_population.prior_functions.sfr_madau_fragos2017_with_bns_td
   ler.gw_source_population.prior_functions.sfr_madau_dickinson2014_with_bns_td
   ler.gw_source_population.prior_functions.sfr_madau_fragos2017
   ler.gw_source_population.prior_functions.sfr_madau_dickinson2014
   ler.gw_source_population.prior_functions.binary_masses_BBH_popIII_lognormal_rvs
   ler.gw_source_population.prior_functions.binary_masses_BBH_primordial_lognormal_rvs
   ler.gw_source_population.prior_functions.binary_masses_BNS_bimodal_rvs
   ler.gw_source_population.prior_functions.binary_masses_NSBH_broken_powerlaw_rvs
   ler.gw_source_population.prior_functions.binary_masses_BBH_powerlaw_gaussian_rvs
   ler.gw_source_population.prior_functions.available_prior_list



.. py:function:: merger_rate_density_bbh_oguri2018_function(zs, R0=19 * 1e-09, b2=1.6, b3=2.1, b4=30)

   
   Compute the merger rate density for PopI/II BBH.

   Reference: Oguri et al. (2018). The output is in detector frame and is
   unnormalized.

   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **R0** : ``float``
           Local merger rate density at low redshift (Mpc^-3 yr^-1).

           default: 19e-9 (GWTC-4)

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

   >>> from ler.gw_source_population import merger_rate_density_bbh_oguri2018
   >>> rate_density = merger_rate_density_bbh_oguri2018(zs=np.array([0.1]))



   ..
       !! processed by numpydoc !!

.. py:function:: merger_rate_density_bbh_popIII_ken2022_function(zs, n0=19.2 * 1e-09, aIII=0.66, bIII=0.3, zIII=11.6)

   
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
   >>> rate_density = merger_rate_density_bbh_popIII_ken2022(zs=np.array([0.1]))



   ..
       !! processed by numpydoc !!

.. py:function:: merger_rate_density_madau_dickinson2014_function(zs, R0=19 * 1e-09, a=0.015, b=2.7, c=2.9, d=5.6)

   
   Compute the merger rate density for BBH using Madau & Dickinson (2014) model.

   Reference: Eqn. 15 of https://arxiv.org/pdf/1403.0007

   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **R0** : ``float``
           Local merger rate density (Mpc^-3 yr^-1).

           default: 19e-9

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

       **rate_density** : ``float`` or ``numpy.ndarray``
           Merger rate density (Mpc^-3 yr^-1).










   .. rubric:: Examples

   >>> from ler.gw_source_population import merger_rate_density_madau_dickinson2014
   >>> rate_density = merger_rate_density_madau_dickinson2014(zs=np.array([0.1]))



   ..
       !! processed by numpydoc !!

.. py:function:: merger_rate_density_madau_dickinson_belczynski_ng_function(zs, R0=19 * 1e-09, alpha_F=2.57, beta_F=5.83, c_F=3.36)

   
   Compute BBH merger rate density following Ng et al. (2021).

   This model uses a Madau-Dickinson-like functional form to fit the
   merger rate density of field BHs, accounting for time delays and
   metallicity effects. Coefficients from Madau & Dickinson (2014) are translated as: B-> alpha_F, D-> beta_F, C-> c_F.

   density(zs) ∝ (1 + zs) ** alpha_F / (1 + ((1 + zs) / c_F) ** beta_F)

   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **R0** : ``float``
           Local merger rate density (Mpc^-3 yr^-1).

           default: 19e-9

       **alpha_F** : ``float``
           Low-redshift power-law slope.

           default: 2.57

       **beta_F** : ``float``
           High-redshift power-law slope.

           default: 5.83

       **c_F** : ``float``
           Turnover redshift parameter.

           default: 3.36

   :Returns:

       **rate_density** : ``float`` or ``numpy.ndarray``
           Merger rate density (Mpc^-3 yr^-1).










   .. rubric:: Examples

   >>> from ler.gw_source_population import merger_rate_density_madau_dickinson_belczynski_ng
   >>> rate_density = merger_rate_density_madau_dickinson_belczynski_ng(zs=np.array([0.1]))



   ..
       !! processed by numpydoc !!

.. py:function:: merger_rate_density_bbh_primordial_ken2022_function(zs, cosmology=None, n0=0.044 * 1e-09, t0=13.786885302009708)

   
   Compute the merger rate density for Primordial BBH.

   Reference: Ng et al. (2022). The output is in detector frame and is
   unnormalized.

   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **cosmology** : ``astropy.cosmology`` or ``None``
           Cosmology object for age calculations.

           default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)

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
   >>> rate_density = merger_rate_density_bbh_primordial_ken2022(zs=np.array([0.1]))



   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_fragos2017_with_bbh_td(zs, R0=19 * 1e-09)

   
   Compute the merger rate density for BBH. This is computed from star formation rate, Madau & Fragos (2017), with an additional time delay. This function is relies on pre-generated data points.


   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **R0** : ``float``
           Local merger rate density (Mpc^-3 yr^-1).

           default: 19e-9

   :Returns:

       **SFR** : ``float`` or ``numpy.ndarray``
           Star formation rate (Mpc^-3 yr^-1).













   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_dickinson2014_with_bbh_td(zs, R0=19 * 1e-09)

   
   Compute the merger rate density for BBH. This is computed from star formation rate, Madau & Dickinson (2014), with an additional time delay. This function is relies on pre-generated data points.


   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **R0** : ``float``
           Local merger rate density (Mpc^-3 yr^-1).

           default: 19e-9

   :Returns:

       **SFR** : ``float`` or ``numpy.ndarray``
           Star formation rate (Mpc^-3 yr^-1).













   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_fragos2017_with_bns_td(zs, R0=89 * 1e-09)

   
   Compute the merger rate density for BNS. This is computed from star formation rate, Madau & Fragos (2017), with an additional time delay. This function is relies on pre-generated data points.


   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **R0** : ``float``
           Local merger rate density (Mpc^-3 yr^-1).

           default: 89e-9

   :Returns:

       **SFR** : ``float`` or ``numpy.ndarray``
           Star formation rate (Mpc^-3 yr^-1).













   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_dickinson2014_with_bns_td(zs, R0=89 * 1e-09)

   
   Compute the merger rate density for BNS. This is computed from star formation rate, Madau & Dickinson (2014), with an additional time delay. This function is relies on pre-generated data points.


   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **R0** : ``float``
           Local merger rate density (Mpc^-3 yr^-1).

           default: 89e-9

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
   >>> sfr = sfr_madau_dickinson2014(zs=np.array([0.1]))



   ..
       !! processed by numpydoc !!

.. py:function:: binary_masses_BBH_popIII_lognormal_rvs(size, m_min=1.0, m_max=100.0, Mc=20.0, sigma=0.3, chunk_size=10000)

   
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

   >>> from ler.gw_source_population import binary_masses_BBH_popIII_lognormal
   >>> m1, m2 = binary_masses_BBH_popIII_lognormal(size=1000)



   ..
       !! processed by numpydoc !!

.. py:function:: binary_masses_BBH_primordial_lognormal_rvs(size, m_min=1.0, m_max=100.0, Mc=20.0, sigma=0.3, chunk_size=10000)

   
   Sample from a lognormal distribution in 2D mass space.

   Based on Eqn. 1 and 4 of Ng et al. 2022 for primordial black holes.

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

   >>> from ler.gw_source_population import binary_masses_BBH_primordial_lognormal
   >>> m1, m2 = binary_masses_BBH_primordial_lognormal(size=1000)



   ..
       !! processed by numpydoc !!

.. py:function:: binary_masses_BNS_bimodal_rvs(size, w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1.0, mmax=2.3, resolution=500)

   
   Sample BNS masses from bimodal Gaussian distribution.

   Based on Will M. Farr et al. 2020 Eqn. 6 for neutron star mass
   distribution combining two Gaussian peaks.

   :Parameters:

       **size** : ``int``
           Number of samples to draw.

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

       **resolution** : ``int``
           Number of points to use for the CDF.

           default: 500

   :Returns:

       **m1** : ``numpy.ndarray``
           Primary mass samples (Msun).

       **m2** : ``numpy.ndarray``
           Secondary mass samples (Msun).













   ..
       !! processed by numpydoc !!

.. py:function:: binary_masses_NSBH_broken_powerlaw_rvs(size=1000, mminbh=26.0, mmaxbh=125.0, alpha_1=6.75, alpha_2=0.0, b=0.5, delta_m=5.0, mminns=1.0, mmaxns=3.0, alphans=0.0, normalization_size=1000)

   
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

.. py:function:: binary_masses_BBH_powerlaw_gaussian_rvs(size, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m, beta, normalization_size=1000)

   
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

.. py:function:: available_prior_list()

   
   Returns a list of available priors.
















   ..
       !! processed by numpydoc !!

