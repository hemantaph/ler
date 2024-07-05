:orphan:

:py:mod:`ler.lens_galaxy_population.jit_functions`
==================================================

.. py:module:: ler.lens_galaxy_population.jit_functions


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.jit_functions.axis_ratio_SIS
   ler.lens_galaxy_population.jit_functions.phi
   ler.lens_galaxy_population.jit_functions.phi_loc_bernardi
   ler.lens_galaxy_population.jit_functions.phi_cut_SIE
   ler.lens_galaxy_population.jit_functions.axis_ratio_rayleigh
   ler.lens_galaxy_population.jit_functions.velocity_dispersion_z_dependent
   ler.lens_galaxy_population.jit_functions.lens_redshift_SDSS_catalogue
   ler.lens_galaxy_population.jit_functions.bounded_normal_sample



.. py:function:: axis_ratio_SIS(sigma)

   
   Function to sample axis ratio from the SIS distribution with given velocity dispersion.


   :Parameters:

       **sigma** : `float: array`
           velocity dispersion of the lens galaxy

   :Returns:

       **q** : `float: array`
           axis ratio of the lens galaxy













   ..
       !! processed by numpydoc !!

.. py:function:: phi(s, z, cosmology_h=0.7)

   
   Function to calculate the lens galaxy velocity dispersion function at redshift z.


   :Parameters:

       **s** : `float: array`
           velocity dispersion of the lens galaxy

       **z** : `float: array`
           redshift of the lens galaxy

       **cosmology_h** : `float`
           Hubble constant

   :Returns:

       **result** : `float: array`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: phi_loc_bernardi(sigma, alpha=0.94, beta=1.85, phistar=0.02099, sigmastar=113.78, cosmology_h=0.7)

   
   Function to calculate the local universe velocity dispersion function. Bernardi et al. (2010).


   :Parameters:

       **sigma** : `float: array`
           velocity dispersion of the lens galaxy

       **alpha, beta, phistar, sigmastar** : `float`
           parameters of the velocity dispersion function

       **cosmology_h** : `float`
           Hubble constant with respect to 100 km/s/Mpc

   :Returns:

       **philoc_** : `float: array`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: phi_cut_SIE(q)

   
   Function to calculate cross-section scaling factor for the SIE lens galaxy from SIS lens galaxy.


   :Parameters:

       **q** : `float: array`
           axis ratio of the lens galaxy

   :Returns:

       **result** : `float: array`
           scaling factor













   ..
       !! processed by numpydoc !!

.. py:function:: axis_ratio_rayleigh(sigma, q_min=0.2, q_max=1.0)

   
   Function to sample axis ratio from rayleigh distribution with given velocity dispersion.


   :Parameters:

       **sigma** : `float: array`
           velocity dispersion of the lens galaxy

   :Returns:

       **q** : `float: array`
           axis ratio of the lens galaxy













   ..
       !! processed by numpydoc !!

.. py:function:: velocity_dispersion_z_dependent(size, zl, zl_list, vd_inv_cdf)

   
   Function to sample velocity dispersion from the interpolator


   :Parameters:

       **size: int**
           Number of samples to draw

       **zl: `numpy.ndarray` (1D array of float of size=size)**
           Redshift of the lens galaxy

   :Returns:

       samples: numpy.ndarray
           Samples of velocity dispersion













   ..
       !! processed by numpydoc !!

.. py:function:: lens_redshift_SDSS_catalogue(zs, splineDc, splineDcInv, u, cdf)

   
   Function to sample lens redshift from the SDSS catalogue.


   :Parameters:

       **zs: `numpy.ndarray` (1D array of float of size=size)**
           Redshift of the source galaxy

       **splineDc: `list`**
           List of spline coefficients for the comoving distance and redshifts

       **splineDcInv: `list`**
           List of spline coefficients for the inverse of comoving distance and redshifts

       **u: `numpy.ndarray` (1D array of float of size=size)**
           corresponding x values wrt to the cdf values
           e.g. u = np.linspace(0, 1, 500)

       **cdf: `numpy.ndarray` (1D array of float of size=size)**
           Cumulative distribution function of the lens redshift distribution between 0 and 1

   :Returns:

       zl: `numpy.ndarray` (1D array of float of size=size)
           Redshift of the lens galaxy corresponding to the zs













   ..
       !! processed by numpydoc !!

.. py:function:: bounded_normal_sample(size, mean, std, low, high)

   
   Function to sample from a normal distribution with bounds.


   :Parameters:

       **mean: `float`**
           Mean of the normal distribution

       **std: `float`**
           Standard deviation of the normal distribution

       **low: `float`**
           Lower bound

       **high: `float`**
           Upper bound














   ..
       !! processed by numpydoc !!

