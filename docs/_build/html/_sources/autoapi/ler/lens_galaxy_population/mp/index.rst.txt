:py:mod:`ler.lens_galaxy_population.mp`
=======================================

.. py:module:: ler.lens_galaxy_population.mp

.. autoapi-nested-parse::

   Multiprocessing helper functions for lens galaxy population calculations.

   This module provides parallelized functions for computing optical depth
   and cross-sections for strong gravitational lensing. These functions are
   designed to be used with Python's multiprocessing module for efficient
   Monte Carlo integration over lens parameters.

   Key Features:

   - JIT-compiled parallel function for lens redshift sampling

   - Multiprocessing-compatible wrapper functions for optical depth

   - Cross-section calculation using EPL+Shear lens models


   Copyright (C) 2026 Phurailatpam Hemantakumar. Distributed under MIT License.

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.mp.lens_redshift_strongly_lensed_njit
   ler.lens_galaxy_population.mp.lens_redshift_strongly_lensed_mp
   ler.lens_galaxy_population.mp.cross_section_unit_mp
   ler.lens_galaxy_population.mp.cross_section_mp



Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.mp.CS_UNIT_SLOPE
   ler.lens_galaxy_population.mp.CS_UNIT_INTERCEPT


.. py:data:: CS_UNIT_SLOPE
   :value: '0.31830988618379075'

   

.. py:data:: CS_UNIT_INTERCEPT
   :value: '-3.2311742677852644e-27'

   

.. py:function:: lens_redshift_strongly_lensed_njit(zs_array, zl_scaled, sigma_min, sigma_max, q_rvs, phi_rvs, gamma_rvs, shear_rvs, number_density, cross_section, dVcdz_function, integration_size)

   
   JIT-compiled parallel computation of lens redshift optical depth.

   Computes the differential optical depth for strong lensing as a function
   of source and lens redshifts using Monte Carlo integration with parallel
   execution over source redshifts.

   :Parameters:

       **zs_array** : ``numpy.ndarray``
           1D array of source redshifts.

       **zl_scaled** : ``numpy.ndarray``
           2D array of scaled lens redshifts (zl/zs).

       **sigma_min** : ``float``
           Minimum velocity dispersion (km/s).

       **sigma_max** : ``float``
           Maximum velocity dispersion (km/s).

       **q_rvs** : ``callable``
           Function to sample axis ratios given size and sigma.

       **phi_rvs** : ``callable``
           Function to sample axis rotation angles.

       **gamma_rvs** : ``callable``
           Function to sample density profile slopes.

       **shear_rvs** : ``callable``
           Function to sample external shear components (gamma1, gamma2).

       **number_density** : ``callable``
           Function to compute velocity dispersion number density.

       **cross_section** : ``callable``
           Function to compute lensing cross-section.

       **dVcdz_function** : ``callable``
           Function to compute differential comoving volume.

       **integration_size** : ``int``
           Number of Monte Carlo samples per (zs, zl) pair.

   :Returns:

       **result_array** : ``numpy.ndarray``
           2D array of optical depth values with shape (len(zs_array), len(zl_scaled[0])).













   ..
       !! processed by numpydoc !!

.. py:function:: lens_redshift_strongly_lensed_mp(params)

   
   Multiprocessing worker for lens redshift optical depth calculation.

   Computes the differential optical depth for a single source redshift
   across multiple lens redshifts. Designed to be called via multiprocessing
   Pool.map() for parallel computation.

   :Parameters:

       **params** : ``tuple``
           Packed parameters tuple containing:

           - params[0]: Source redshift (float)

           - params[1]: Scaled lens redshift array (1D array)

           - params[2]: Worker index (int)

   :Returns:

       **idx** : ``int``
           Worker index for result ordering.

       **result_array** : ``numpy.ndarray``
           1D array of optical depth values for each lens redshift.













   ..
       !! processed by numpydoc !!

.. py:function:: cross_section_unit_mp(params)

   
   Multiprocessing worker for unit Einstein radius cross-section.

   Computes the lensing cross-section for a lens with unit Einstein radius
   (theta_E = 1). Used for building cross-section interpolation grids.

   :Parameters:

       **params** : ``tuple``
           Packed parameters (e1, e2, gamma, gamma1, gamma2, idx).

   :Returns:

       **idx** : ``int``
           Worker index for result ordering.

       **area** : ``float``
           Cross-section area in square arcseconds.













   ..
       !! processed by numpydoc !!

.. py:function:: cross_section_mp(params)

   
   Multiprocessing worker for cross-section calculation.

   Computes the lensing cross-section for given lens parameters.
   Designed to be called via multiprocessing Pool.map().

   :Parameters:

       **params** : ``tuple``
           Packed parameters (theta_E, e1, e2, gamma, gamma1, gamma2, idx).

   :Returns:

       **idx** : ``int``
           Worker index for result ordering.

       **area** : ``float``
           Cross-section area in square arcseconds.













   ..
       !! processed by numpydoc !!

