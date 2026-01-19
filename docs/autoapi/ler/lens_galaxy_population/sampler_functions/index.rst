:py:mod:`ler.lens_galaxy_population.sampler_functions`
======================================================

.. py:module:: ler.lens_galaxy_population.sampler_functions

.. autoapi-nested-parse::

   Module for lens galaxy parameter sampling functions.

   This module provides probability density functions (PDFs) and random variable
   samplers (RVS) for lens galaxy parameters including redshift, velocity dispersion,
   and axis ratio. It also includes rejection and importance sampling algorithms
   for sampling lens parameters weighted by gravitational lensing cross sections.

   Key Components:

   - Lens redshift samplers (SIS model from Haris et al. 2018)

   - Velocity dispersion samplers (generalized gamma distribution)

   - Axis ratio samplers (Rayleigh and Padilla-Strauss distributions)

   - Rejection and importance sampling for cross-section weighted parameters


   Copyright (C) 2026 Hemantakumar Phurailatpam. Distributed under MIT License.

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.sampler_functions.available_sampler_list
   ler.lens_galaxy_population.sampler_functions.lens_redshift_strongly_lensed_sis_haris_pdf
   ler.lens_galaxy_population.sampler_functions.lens_redshift_strongly_lensed_sis_haris_rvs
   ler.lens_galaxy_population.sampler_functions.velocity_dispersion_ewoud_denisty_function
   ler.lens_galaxy_population.sampler_functions.velocity_dispersion_bernardi_denisty_function
   ler.lens_galaxy_population.sampler_functions.velocity_dispersion_gengamma_density_function
   ler.lens_galaxy_population.sampler_functions.velocity_dispersion_gengamma_pdf
   ler.lens_galaxy_population.sampler_functions.velocity_dispersion_gengamma_rvs
   ler.lens_galaxy_population.sampler_functions.axis_ratio_rayleigh_rvs
   ler.lens_galaxy_population.sampler_functions.axis_ratio_rayleigh_pdf
   ler.lens_galaxy_population.sampler_functions.axis_ratio_padilla_strauss_rvs
   ler.lens_galaxy_population.sampler_functions.axis_ratio_padilla_strauss_pdf
   ler.lens_galaxy_population.sampler_functions.bounded_normal_sample
   ler.lens_galaxy_population.sampler_functions.rejection_sampler
   ler.lens_galaxy_population.sampler_functions.create_rejection_sampler
   ler.lens_galaxy_population.sampler_functions.importance_sampler
   ler.lens_galaxy_population.sampler_functions.importance_sampler_mp
   ler.lens_galaxy_population.sampler_functions.create_importance_sampler



.. py:function:: available_sampler_list()

   
   Return list of available lens parameter samplers.



   :Returns:

       **sampler_list** : ``list``
           List of available sampler function names.










   .. rubric:: Examples

   >>> samplers = available_sampler_list()
   >>> print(samplers)
   ['lens_redshift_strongly_lensed_sis_haris', 'velocity_dispersion_gengamma', 'axis_ratio_rayleigh', 'axis_ratio_padilla_strauss']



   ..
       !! processed by numpydoc !!

.. py:function:: lens_redshift_strongly_lensed_sis_haris_pdf(zl, zs, cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7))

   
   Compute lens redshift PDF for SIS model (Haris et al. 2018).

   Computes the probability density function for lens redshift between
   zl=0 and zl=zs using the analytical form from Haris et al. (2018)
   equation A7, based on the SIS (Singular Isothermal Sphere) lens model.

   :Parameters:

       **zl** : ``float``
           Redshift of the lens galaxy.

       **zs** : ``float``
           Redshift of the source.

       **cosmo** : ``astropy.cosmology``
           Cosmology object for distance calculations.

           default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

   :Returns:

       **pdf** : ``float``
           Probability density at the given lens redshift.










   .. rubric:: Examples

   >>> from astropy.cosmology import LambdaCDM
   >>> cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
   >>> pdf = lens_redshift_strongly_lensed_sis_haris_pdf(zl=0.5, zs=1.0, cosmo=cosmo)
   >>> print(f"PDF at zl=0.5: {pdf:.4f}")
   PDF at zl=0.5: 1.8750



   ..
       !! processed by numpydoc !!

.. py:function:: lens_redshift_strongly_lensed_sis_haris_rvs(size, zs, z_min=0.001, z_max=10.0, cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7))

   
   Sample lens redshifts for SIS model (Haris et al. 2018).

   Uses inverse transform sampling with the analytical CDF of the Haris et al.
   (2018) lens redshift distribution to efficiently generate samples.

   :Parameters:

       **size** : ``int``
           Number of samples to draw.

       **zs** : ``float``
           Redshift of the source.

       **z_min** : ``float``
           Minimum redshift for interpolation grid.

           default: 0.001

       **z_max** : ``float``
           Maximum redshift for interpolation grid.

           default: 10.0

       **cosmo** : ``astropy.cosmology``
           Cosmology object for distance calculations.

           default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

   :Returns:

       **zl** : ``numpy.ndarray``
           Array of sampled lens redshifts with shape (size,).










   .. rubric:: Examples

   >>> from astropy.cosmology import LambdaCDM
   >>> cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
   >>> zl_samples = lens_redshift_strongly_lensed_sis_haris_rvs(
   ...     size=1000,
   ...     zs=1.5,
   ...     z_min=0.001,
   ...     z_max=10.0,
   ...     cosmo=cosmo,
   ... )
   >>> print(f"Mean lens redshift: {zl_samples.mean():.2f}")
   >>> print(f"Redshift range: [{zl_samples.min():.3f}, {zl_samples.max():.3f}]")



   ..
       !! processed by numpydoc !!

.. py:function:: velocity_dispersion_ewoud_denisty_function(sigma, z, alpha=0.94, beta=1.85, phistar=0.02099, sigmastar=113.78)

   
   Calculate the lens galaxy velocity dispersion function at redshift z (Oguri et al. (2018b) + Wempe et al. (2022)).


   :Parameters:

       **sigma** : ``numpy.ndarray``
           Velocity dispersion of the lens galaxy (km/s).

       **z** : ``float``
           Redshift of the lens galaxy.

       **alpha** : ``float``
           Shape parameter of the velocity dispersion function.

           default: 0.94

       **beta** : ``float``
           Slope parameter of the velocity dispersion function.

           default: 1.85

       **phistar** : ``float``
           Normalization of the velocity dispersion function (Mpc^-3).

           default: 2.099e-2

       **sigmastar** : ``float``
           Characteristic velocity dispersion (km/s).

           default: 113.78

   :Returns:

       **result** : ``numpy.ndarray``
           Velocity dispersion function values.

           Negative values are clipped to 0.













   ..
       !! processed by numpydoc !!

.. py:function:: velocity_dispersion_bernardi_denisty_function(sigma, alpha, beta, phistar, sigmastar)

   
   Calculate the local universe velocity dispersion function.

   This implements the velocity dispersion function from Bernardi et al. (2010).

   :Parameters:

       **sigma** : ``numpy.ndarray``
           Velocity dispersion of the lens galaxy (km/s).

       **alpha** : ``float``
           Shape parameter (alpha/beta is used in the gamma function).

           For Oguri et al. (2018b): alpha=0.94

           For Choi et al. (2008): alpha=2.32/2.67

       **beta** : ``float``
           Slope parameter of the velocity dispersion function.

           For Oguri et al. (2018b): beta=1.85

           For Choi et al. (2008): beta=2.67

       **phistar** : ``float``
           Normalization of the velocity dispersion function (Mpc^-3).

           For Oguri et al. (2018b): phistar=2.099e-2*(h/0.7)^3

           For Choi et al. (2008): phistar=8.0e-3*h^3

       **sigmastar** : ``float``
           Characteristic velocity dispersion (km/s).

           For Oguri et al. (2018b): sigmastar=113.78

           For Choi et al. (2008): sigmastar=161.0

   :Returns:

       **philoc** : ``numpy.ndarray``
           Local velocity dispersion function values.













   ..
       !! processed by numpydoc !!

.. py:function:: velocity_dispersion_gengamma_density_function(sigma, alpha=0.94, beta=1.85, phistar=0.02099, sigmastar=113.78, **kwargs)

   
   Compute unnormalized velocity dispersion function using generalized gamma.

   Computes the galaxy velocity dispersion function (VDF) using the generalized
   gamma distribution formulation from Choi et al. (2007).

   :Parameters:

       **sigma** : ``float`` or ``numpy.ndarray``
           Velocity dispersion in km/s.

       **alpha** : ``float``
           Power-law index governing the low-velocity slope.

           default: 0.94

       **beta** : ``float``
           Exponential parameter for high-velocity cutoff sharpness.

           default: 1.85

       **phistar** : ``float``
           Normalization constant (comoving number density, Mpc^-3).

           default: 8.0e-3

       **sigmastar** : ``float``
           Characteristic velocity scale in km/s.

           default: 113.78

   :Returns:

       **density** : ``float`` or ``numpy.ndarray``
           Unnormalized velocity dispersion function value.










   .. rubric:: Examples

   >>> import numpy as np
   >>> sigma = np.array([150.0, 200.0, 250.0])
   >>> density = velocity_dispersion_gengamma_density_function(sigma)
   >>> print(f"Density at sigma=200 km/s: {density[1]:.6f}")



   ..
       !! processed by numpydoc !!

.. py:function:: velocity_dispersion_gengamma_pdf(sigma, sigma_min=100.0, sigma_max=400.0, alpha=0.94, beta=1.85, sigmastar=113.78)

   
   Compute normalized velocity dispersion PDF using generalized gamma.

   Computes the probability density function for velocity dispersion using
   the generalized gamma distribution, normalized over the specified range.

   :Parameters:

       **sigma** : ``float`` or ``numpy.ndarray``
           Velocity dispersion in km/s.

       **sigma_min** : ``float``
           Minimum velocity dispersion for normalization (km/s).

           default: 100.0

       **sigma_max** : ``float``
           Maximum velocity dispersion for normalization (km/s).

           default: 400.0

       **alpha** : ``float``
           Power-law index governing the low-velocity slope.

           default: 0.94

       **beta** : ``float``
           Exponential parameter for high-velocity cutoff sharpness.

           default: 1.85

       **sigmastar** : ``float``
           Characteristic velocity scale in km/s.

           default: 113.78

   :Returns:

       **pdf** : ``float`` or ``numpy.ndarray``
           Normalized probability density at the given velocity dispersion.










   .. rubric:: Examples

   >>> pdf = velocity_dispersion_gengamma_pdf(
   ...     sigma=200.0,
   ...     sigma_min=100.0,
   ...     sigma_max=400.0,
   ... )
   >>> print(f"PDF at sigma=200 km/s: {pdf:.6f}")



   ..
       !! processed by numpydoc !!

.. py:function:: velocity_dispersion_gengamma_rvs(size, sigma_min=100.0, sigma_max=400.0, alpha=0.94, beta=1.85, sigmastar=113.78)

   
   Sample velocity dispersions from generalized gamma distribution.

   Uses truncated generalized gamma sampling via rejection to generate
   velocity dispersion samples within the specified bounds.

   :Parameters:

       **size** : ``int``
           Number of samples to draw.

       **sigma_min** : ``float``
           Minimum velocity dispersion (km/s).

           default: 100.0

       **sigma_max** : ``float``
           Maximum velocity dispersion (km/s).

           default: 400.0

       **alpha** : ``float``
           Power-law index governing the low-velocity slope.

           default: 0.94

       **beta** : ``float``
           Exponential parameter for high-velocity cutoff sharpness.

           default: 1.85

       **sigmastar** : ``float``
           Characteristic velocity scale in km/s.

           default: 113.78

   :Returns:

       **sigma** : ``numpy.ndarray``
           Sampled velocity dispersions in km/s with shape (size,).










   .. rubric:: Examples

   >>> sigma_samples = velocity_dispersion_gengamma(
   ...     size=1000,
   ...     sigma_min=100.0,
   ...     sigma_max=400.0,
   ... )
   >>> print(f"Mean sigma: {sigma_samples.mean():.2f} km/s")
   >>> print(f"Range: [{sigma_samples.min():.1f}, {sigma_samples.max():.1f}] km/s")



   ..
       !! processed by numpydoc !!

.. py:function:: axis_ratio_rayleigh_rvs(size, sigma, q_min=0.2, q_max=1.0)

   
   Sample axis ratios from velocity-dependent Rayleigh distribution.

   Generates axis ratio samples using the Rayleigh distribution with
   scale parameter dependent on velocity dispersion, as described in
   Wierda et al. (2021) Appendix C.

   :Parameters:

       **size** : ``int``
           Number of samples to draw.

       **sigma** : ``numpy.ndarray``
           Velocity dispersions in km/s with shape (size,).

       **q_min** : ``float``
           Minimum allowed axis ratio.

           default: 0.2

       **q_max** : ``float``
           Maximum allowed axis ratio.

           default: 1.0

   :Returns:

       **q** : ``numpy.ndarray``
           Sampled axis ratios with shape (size,).










   .. rubric:: Examples

   >>> import numpy as np
   >>> sigma = np.random.uniform(100, 300, 1000)
   >>> q_samples = axis_ratio_rayleigh_rvs(size=1000, sigma=sigma, q_min=0.2, q_max=1.0)
   >>> print(f"Mean axis ratio: {q_samples.mean():.2f}")
   >>> print(f"Range: [{q_samples.min():.2f}, {q_samples.max():.2f}]")



   ..
       !! processed by numpydoc !!

.. py:function:: axis_ratio_rayleigh_pdf(q, sigma, q_min=0.2, q_max=1.0)

   
   Compute truncated Rayleigh PDF for axis ratio.

   Computes the probability density function for axis ratio using the
   truncated Rayleigh distribution with velocity-dependent scale parameter
   (Wierda et al. 2021 equation C16).

   :Parameters:

       **q** : ``numpy.ndarray``
           Axis ratios at which to evaluate PDF.

       **sigma** : ``numpy.ndarray``
           Velocity dispersions in km/s (same shape as q).

       **q_min** : ``float``
           Minimum axis ratio for truncation.

           default: 0.2

       **q_max** : ``float``
           Maximum axis ratio for truncation.

           default: 1.0

   :Returns:

       **pdf** : ``numpy.ndarray``
           Probability density values with same shape as q.










   .. rubric:: Examples

   >>> import numpy as np
   >>> q = np.array([0.5, 0.7, 0.9])
   >>> sigma = np.array([150.0, 200.0, 250.0])
   >>> pdf = axis_ratio_rayleigh_pdf(q, sigma)
   >>> print(f"PDF values: {pdf}")



   ..
       !! processed by numpydoc !!

.. py:function:: axis_ratio_padilla_strauss_rvs(size)

   
   Sample axis ratios from Padilla & Strauss (2008) distribution.

   Uses inverse transform sampling with the empirical PDF from
   Padilla & Strauss (2008) for early-type galaxy axis ratios.

   :Parameters:

       **size** : ``int``
           Number of samples to draw.

   :Returns:

       **q** : ``numpy.ndarray``
           Sampled axis ratios with shape (size,).










   .. rubric:: Examples

   >>> q_samples = axis_ratio_padilla_strauss_rvs(size=1000)
   >>> print(f"Mean axis ratio: {q_samples.mean():.2f}")
   >>> print(f"Range: [{q_samples.min():.2f}, {q_samples.max():.2f}]")



   ..
       !! processed by numpydoc !!

.. py:function:: axis_ratio_padilla_strauss_pdf(q)

   
   Compute axis ratio PDF from Padilla & Strauss (2008).

   Evaluates the probability density function for axis ratio using
   cubic spline interpolation of the Padilla & Strauss (2008) data.

   :Parameters:

       **q** : ``numpy.ndarray``
           Axis ratios at which to evaluate PDF.

   :Returns:

       **pdf** : ``numpy.ndarray``
           Probability density values with same shape as q.










   .. rubric:: Examples

   >>> import numpy as np
   >>> q = np.array([0.3, 0.5, 0.7, 0.9])
   >>> pdf = axis_ratio_padilla_strauss_pdf(q)
   >>> print(f"PDF at q=0.5: {pdf[1]:.4f}")



   ..
       !! processed by numpydoc !!

.. py:function:: bounded_normal_sample(size, mean, std, low, high)

   
   Sample from truncated normal distribution via rejection.

   Generates samples from a normal distribution with specified mean and
   standard deviation, rejecting samples outside the specified bounds.

   :Parameters:

       **size** : ``int``
           Number of samples to draw.

       **mean** : ``float``
           Mean of the normal distribution.

       **std** : ``float``
           Standard deviation of the normal distribution.

       **low** : ``float``
           Lower bound for samples.

       **high** : ``float``
           Upper bound for samples.

   :Returns:

       **samples** : ``numpy.ndarray``
           Bounded normal samples with shape (size,).










   .. rubric:: Examples

   >>> samples = bounded_normal_sample(size=1000, mean=2.0, std=0.2, low=1.5, high=2.5)
   >>> print(f"Mean: {samples.mean():.2f}, Std: {samples.std():.2f}")
   >>> print(f"Range: [{samples.min():.2f}, {samples.max():.2f}]")



   ..
       !! processed by numpydoc !!

.. py:function:: rejection_sampler(zs, zl, sigma_max, sigma_rvs, q_rvs, phi_rvs, gamma_rvs, shear_rvs, cross_section, safety_factor=1.2)

   
   Core rejection sampling algorithm for lens parameters.


   :Parameters:

       **zs** : ``numpy.ndarray``
           Source redshifts.

       **zl** : ``numpy.ndarray``
           Lens redshifts.

       **sigma_max** : ``float``
           Maximum velocity dispersion (km/s) for computing upper bound.

       **sigma_rvs** : ``callable``
           Function to sample velocity dispersion: sigma_rvs(n, zl) -> array.

       **q_rvs** : ``callable``
           Function to sample axis ratio: q_rvs(n, sigma) -> array.

       **phi_rvs** : ``callable``
           Function to sample orientation angle: phi_rvs(n) -> array.

       **gamma_rvs** : ``callable``
           Function to sample power-law index: gamma_rvs(n) -> array.

       **shear_rvs** : ``callable``
           Function to sample external shear: shear_rvs(n) -> (gamma1, gamma2).

       **cross_section** : ``callable``
           Function to compute lensing cross section.

       **safety_factor** : ``float``
           Multiplicative safety factor for the upper bound.

           default: 1.2

   :Returns:

       **sigma_array** : ``numpy.ndarray``
           Sampled velocity dispersions (km/s).

       **q_array** : ``numpy.ndarray``
           Sampled axis ratios.

       **phi_array** : ``numpy.ndarray``
           Sampled orientation angles (rad).

       **gamma_array** : ``numpy.ndarray``
           Sampled power-law indices.

       **gamma1_array** : ``numpy.ndarray``
           Sampled external shear component 1.

       **gamma2_array** : ``numpy.ndarray``
           Sampled external shear component 2.













   ..
       !! processed by numpydoc !!

.. py:function:: create_rejection_sampler(sigma_max, sigma_rvs, q_rvs, phi_rvs, gamma_rvs, shear_rvs, cross_section, safety_factor=1.2, use_njit_sampler=True)

   
   Create a rejection sampler for cross-section weighted lens parameters.

   Returns a callable that samples lens parameters using rejection sampling,
   weighting by the gravitational lensing cross section. Optionally uses
   Numba JIT compilation for improved performance.

   :Parameters:

       **sigma_max** : ``float``
           Maximum velocity dispersion (km/s) for computing upper bound.

       **sigma_rvs** : ``callable``
           Function to sample velocity dispersion: sigma_rvs(n, zl) -> array.

       **q_rvs** : ``callable``
           Function to sample axis ratio: q_rvs(n, sigma) -> array.

       **phi_rvs** : ``callable``
           Function to sample orientation angle: phi_rvs(n) -> array.

       **gamma_rvs** : ``callable``
           Function to sample power-law index: gamma_rvs(n) -> array.

       **shear_rvs** : ``callable``
           Function to sample external shear: shear_rvs(n) -> (gamma1, gamma2).

       **cross_section** : ``callable``
           Function to compute lensing cross section.

       **safety_factor** : ``float``
           Multiplicative safety factor for the upper bound.

           default: 1.2

       **use_njit_sampler** : ``bool``
           If True, uses Numba JIT compilation for faster execution.

           default: True

   :Returns:

       **rejection_sampler_wrapper** : ``callable``
           Function with signature (zs, zl) -> (sigma, q, phi, gamma, gamma1, gamma2).










   .. rubric:: Examples

   >>> import numpy as np
   >>> from numba import njit
   >>> @njit
   ... def sigma_rvs(n, zl):
   ...     return 100 + 200 * np.random.random(n)
   >>> @njit
   ... def q_rvs(n, sigma):
   ...     return 0.5 + 0.5 * np.random.random(n)
   >>> @njit
   ... def phi_rvs(n):
   ...     return np.pi * np.random.random(n)
   >>> @njit
   ... def gamma_rvs(n):
   ...     return 2.0 + 0.2 * np.random.randn(n)
   >>> @njit
   ... def shear_rvs(n):
   ...     return 0.05 * np.random.randn(n), 0.05 * np.random.randn(n)
   >>> @njit
   ... def cross_section(zs, zl, sigma, q, phi, gamma, gamma1, gamma2):
   ...     return sigma**4
   >>> sampler = create_rejection_sampler(
   ...     sigma_max=400.0,
   ...     sigma_rvs=sigma_rvs,
   ...     q_rvs=q_rvs,
   ...     phi_rvs=phi_rvs,
   ...     gamma_rvs=gamma_rvs,
   ...     shear_rvs=shear_rvs,
   ...     cross_section=cross_section,
   ... )
   >>> zs = np.array([1.0, 1.5, 2.0])
   >>> zl = np.array([0.3, 0.5, 0.7])
   >>> sigma, q, phi, gamma, gamma1, gamma2 = sampler(zs, zl)



   ..
       !! processed by numpydoc !!

.. py:function:: importance_sampler(zs, zl, sigma_min, sigma_max, q_rvs, phi_rvs, gamma_rvs, shear_rvs, sigma_pdf, cross_section, n_prop)

   
   Core importance sampling algorithm for lens parameters.


   :Parameters:

       **zs** : ``numpy.ndarray``
           Source redshifts.

       **zl** : ``numpy.ndarray``
           Lens redshifts.

       **sigma_min** : ``float``
           Minimum velocity dispersion (km/s) for uniform proposal.

       **sigma_max** : ``float``
           Maximum velocity dispersion (km/s) for uniform proposal.

       **q_rvs** : ``callable``
           Function to sample axis ratio: q_rvs(n, sigma) -> array.

       **phi_rvs** : ``callable``
           Function to sample orientation angle: phi_rvs(n) -> array.

       **gamma_rvs** : ``callable``
           Function to sample power-law index: gamma_rvs(n) -> array.

       **shear_rvs** : ``callable``
           Function to sample external shear: shear_rvs(n) -> (gamma1, gamma2).

       **sigma_pdf** : ``callable``
           PDF of velocity dispersion: sigma_pdf(sigma, zl) -> array.

       **cross_section** : ``callable``
           Function to compute lensing cross section.

       **n_prop** : ``int``
           Number of proposal samples per lens.

   :Returns:

       **sigma_post** : ``numpy.ndarray``
           Sampled velocity dispersions (km/s).

       **q_post** : ``numpy.ndarray``
           Sampled axis ratios.

       **phi_post** : ``numpy.ndarray``
           Sampled orientation angles (rad).

       **gamma_post** : ``numpy.ndarray``
           Sampled power-law indices.

       **gamma1_post** : ``numpy.ndarray``
           Sampled external shear component 1.

       **gamma2_post** : ``numpy.ndarray``
           Sampled external shear component 2.













   ..
       !! processed by numpydoc !!

.. py:function:: importance_sampler_mp(zs, zl, sigma_min, sigma_max, q_rvs, phi_rvs, gamma_rvs, shear_rvs, sigma_pdf, cross_section, n_prop, npool=4)

   
   Multiprocessing version of importance sampling for lens parameters.


   :Parameters:

       **zs** : ``numpy.ndarray``
           Source redshifts.

       **zl** : ``numpy.ndarray``
           Lens redshifts.

       **sigma_min** : ``float``
           Minimum velocity dispersion (km/s) for uniform proposal.

       **sigma_max** : ``float``
           Maximum velocity dispersion (km/s) for uniform proposal.

       **q_rvs** : ``callable``
           Function to sample axis ratio: q_rvs(n, sigma) -> array.

       **phi_rvs** : ``callable``
           Function to sample orientation angle: phi_rvs(n) -> array.

       **gamma_rvs** : ``callable``
           Function to sample power-law index: gamma_rvs(n) -> array.

       **shear_rvs** : ``callable``
           Function to sample external shear: shear_rvs(n) -> (gamma1, gamma2).

       **sigma_pdf** : ``callable``
           PDF of velocity dispersion: sigma_pdf(sigma, zl) -> array.

       **cross_section** : ``callable``
           Function to compute lensing cross section.

       **n_prop** : ``int``
           Number of proposal samples per lens.

       **npool** : ``int``
           Number of parallel processes to use.

           default: 4

   :Returns:

       **sigma_post** : ``numpy.ndarray``
           Sampled velocity dispersions (km/s).

       **q_post** : ``numpy.ndarray``
           Sampled axis ratios.

       **phi_post** : ``numpy.ndarray``
           Sampled orientation angles (rad).

       **gamma_post** : ``numpy.ndarray``
           Sampled power-law indices.

       **gamma1_post** : ``numpy.ndarray``
           Sampled external shear component 1.

       **gamma2_post** : ``numpy.ndarray``
           Sampled external shear component 2.













   ..
       !! processed by numpydoc !!

.. py:function:: create_importance_sampler(sigma_min, sigma_max, q_rvs, phi_rvs, gamma_rvs, shear_rvs, sigma_pdf, cross_section, n_prop, use_njit_sampler=True, npool=4)

   
   Create an importance sampler for cross-section weighted lens parameters.

   Returns a callable that samples lens parameters using importance sampling
   with uniform proposal distribution, optionally JIT-compiled for improved
   performance.

   :Parameters:

       **sigma_min** : ``float``
           Minimum velocity dispersion (km/s) for uniform proposal.

       **sigma_max** : ``float``
           Maximum velocity dispersion (km/s) for uniform proposal.

       **q_rvs** : ``callable``
           Function to sample axis ratio: q_rvs(n, sigma) -> array.

       **phi_rvs** : ``callable``
           Function to sample orientation angle: phi_rvs(n) -> array.

       **gamma_rvs** : ``callable``
           Function to sample power-law index: gamma_rvs(n) -> array.

       **shear_rvs** : ``callable``
           Function to sample external shear: shear_rvs(n) -> (gamma1, gamma2).

       **sigma_pdf** : ``callable``
           PDF of velocity dispersion: sigma_pdf(sigma, zl) -> array.

       **cross_section** : ``callable``
           Function to compute lensing cross section.

       **n_prop** : ``int``
           Number of proposal samples per lens.

       **use_njit_sampler** : ``bool``
           If True, uses Numba JIT compilation for faster execution.

           default: True

       **npool** : ``int``
           Number of parallel processes (only used when use_njit_sampler=False).

           default: 4

   :Returns:

       **importance_sampler_wrapper** : ``callable``
           Function with signature (zs, zl) -> (sigma, q, phi, gamma, gamma1, gamma2).










   .. rubric:: Examples

   >>> import numpy as np
   >>> from numba import njit
   >>> @njit
   ... def q_rvs(n, sigma):
   ...     return 0.5 + 0.5 * np.random.random(n)
   >>> @njit
   ... def phi_rvs(n):
   ...     return np.pi * np.random.random(n)
   >>> @njit
   ... def gamma_rvs(n):
   ...     return 2.0 + 0.2 * np.random.randn(n)
   >>> @njit
   ... def shear_rvs(n):
   ...     return 0.05 * np.random.randn(n), 0.05 * np.random.randn(n)
   >>> @njit
   ... def sigma_pdf(sigma, zl):
   ...     return np.ones_like(sigma)
   >>> @njit
   ... def cross_section(zs, zl, sigma, q, phi, gamma, gamma1, gamma2):
   ...     return sigma**4
   >>> sampler = create_importance_sampler(
   ...     sigma_min=100.0,
   ...     sigma_max=400.0,
   ...     q_rvs=q_rvs,
   ...     phi_rvs=phi_rvs,
   ...     gamma_rvs=gamma_rvs,
   ...     shear_rvs=shear_rvs,
   ...     sigma_pdf=sigma_pdf,
   ...     cross_section=cross_section,
   ...     n_prop=100,
   ... )
   >>> zs = np.array([1.0, 1.5, 2.0])
   >>> zl = np.array([0.3, 0.5, 0.7])
   >>> sigma, q, phi, gamma, gamma1, gamma2 = sampler(zs, zl)



   ..
       !! processed by numpydoc !!

