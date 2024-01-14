:orphan:

:py:mod:`ler.lens_galaxy_population.mp`
=======================================

.. py:module:: ler.lens_galaxy_population.mp


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.mp.optical_depth_sie2_mp



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

