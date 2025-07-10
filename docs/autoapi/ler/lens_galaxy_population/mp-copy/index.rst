:orphan:

:py:mod:`ler.lens_galaxy_population.mp copy`
============================================

.. py:module:: ler.lens_galaxy_population.mp copy


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.mp copy.optical_depth_sis4_mp
   ler.lens_galaxy_population.mp copy.optical_depth_sis3_mp
   ler.lens_galaxy_population.mp copy.optical_depth_sis1_mp
   ler.lens_galaxy_population.mp copy.optical_depth_sis2_mp
   ler.lens_galaxy_population.mp copy.optical_depth_sie2_mp
   ler.lens_galaxy_population.mp copy.optical_depth_epl_shear1_mp
   ler.lens_galaxy_population.mp copy.optical_depth_epl_shear2_mp



.. py:function:: optical_depth_sis4_mp(zs_array, splineDa, splinedVcdz, size=1000000, min_max=np.array([[50.0, 420.0]]), sigma_args=np.array([2.32, 2.67, 0.0027439999999999995, 161.0]))

   
   Function to calculate the optical depth for SIS lens with velocity dispersion distribution independent on redshift. With montecarlo integration.
















   ..
       !! processed by numpydoc !!

.. py:function:: optical_depth_sis3_mp(params)

   
   Function to calculate the optical depth for SIS lens with velocity dispersion distribution independent on redshift. With montecarlo integration.
















   ..
       !! processed by numpydoc !!

.. py:function:: optical_depth_sis1_mp(params)

   
   Function to calculate the optical depth for SIS lens with velocity dispersion distribution independent on redshift.
















   ..
       !! processed by numpydoc !!

.. py:function:: optical_depth_sis2_mp(params)

   
   Function to calculate the optical depth for SIS lens with velocity dispersion distribution depending on redshift.


   :Parameters:

       **params** : `list`
           list of parameters
           params[0] = zs (source redshift, float)
           params[1] = no (number density of lens galaxies, float)
           params[2] = vd_inv_cdf (velocity dispersion inverse cdf coefficients and redshift list, list). This vd_inv_cdf(s) of each redshift.
           params[3] = splineVcdz (differential comoving volume spline interpolator coefficients, list)
           params[4] = splineDa (angular diameter distance spline interpolator coefficients and redshift list, list)
           params[5] = idx (index to keep track of the operation, int)
           params[6] = zl_list (list of lens redshifts, list). This use for choosing the right vd_inv_cdf(s) for each lens redshifts.














   ..
       !! processed by numpydoc !!

.. py:function:: optical_depth_sie2_mp(params)

   
   Function to calculate the optical depth for SIE lens with velocity dispersion distribution depending on redshift.


   :Parameters:

       **params** : `list`
           list of parameters
           params[0] = zs (source redshift, float)
           params[1] = no (number density of lens galaxies, float)
           params[2] = vd_inv_cdf (velocity dispersion inverse cdf coefficients and redshift list, list). This vd_inv_cdf(s) of each redshift.
           params[3] = splineVcdz (differential comoving volume spline interpolator coefficients, list)
           params[4] = splineDa (angular diameter distance spline interpolator coefficients and redshift list, list)
           params[5] = idx (index to keep track of the operation, int)
           params[6] = zl_list (list of lens redshifts, list). This use for choosing the right vd_inv_cdf(s) for each lens redshifts.














   ..
       !! processed by numpydoc !!

.. py:function:: optical_depth_epl_shear1_mp(params)

   
   Function to calculate the optical depth for SIE lens with velocity dispersion distribution depending on redshift.
















   ..
       !! processed by numpydoc !!

.. py:function:: optical_depth_epl_shear2_mp(params)

   
   Function to calculate the optical depth for EPL+Shear lens with velocity dispersion distribution depending on redshift.


   :Parameters:

       **params** : `list`
           list of parameters
           params[0] = zs (source redshift, float)
           params[1] = no (number density of lens galaxies, float)
           params[2] = vd_inv_cdf (velocity dispersion inverse cdf coefficients and redshift list, list). This vd_inv_cdf(s) of each redshift.
           params[3] = splineVcdz (differential comoving volume spline interpolator coefficients, list)
           params[4] = splineDa (angular diameter distance spline interpolator coefficients and redshift list, list)
           params[5] = idx (index to keep track of the operation, int)
           params[6] = zl_list (list of lens redshifts, list). This use for choosing the right vd_inv_cdf(s) for each lens redshifts.














   ..
       !! processed by numpydoc !!

