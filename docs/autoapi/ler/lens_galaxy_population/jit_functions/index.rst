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
   ler.lens_galaxy_population.jit_functions.axis_ratio_rayleigh_rvs
   ler.lens_galaxy_population.jit_functions.velocity_dispersion_z_dependent
   ler.lens_galaxy_population.jit_functions.lens_redshift_sis_haris
   ler.lens_galaxy_population.jit_functions.bounded_normal_sample
   ler.lens_galaxy_population.jit_functions.phi_q2_ellipticity_hemanta



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

.. py:function:: phi(s, z, alpha=0.94, beta=1.85, phistar=0.02099, sigmastar=113.78)

   
   Function to calculate the lens galaxy velocity dispersion function at redshift z.
   For Oguri et al. (2018b) model: alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78


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

.. py:function:: phi_loc_bernardi(sigma, alpha, beta, phistar, sigmastar)

   
   Function to calculate the local universe velocity dispersion function. Bernardi et al. (2010).
   For Oguri et al. (2018b) model: alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78
   For Choi et al. (2008) model: alpha = 2.32 / 2.67, beta = 2.67, phistar = 8.0e-3*self.cosmo.h**3, sigmastar = 161.0


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

.. py:function:: axis_ratio_rayleigh_rvs(sigma, q_min=0.2, q_max=1.0)

   
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

.. py:function:: lens_redshift_sis_haris(zs, splineDc, splineDcInv, u, cdf)

   
   Function to sample lens redshift from the SDSS catalogue. Haris et al. (2018) cdf = (10 * u**3 - 15 * u**4 + 6 * u**5)


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

.. py:function:: phi_q2_ellipticity_hemanta(phi, q)

   
   Function to convert phi and q to ellipticity e1 and e2.


   :Parameters:

       **phi** : `float: array`
           angle of the major axis in radians

       **q** : `float: array`
           axis ratio

   :Returns:

       **e1** : `float: array`
           ellipticity component 1

       **e2** : `float: array`
           ..













   ..
       !! processed by numpydoc !!

