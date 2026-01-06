:py:mod:`ler.lens_galaxy_population.cross_section_interpolator`
===============================================================

.. py:module:: ler.lens_galaxy_population.cross_section_interpolator

.. autoapi-nested-parse::

   Cross Section Interpolator Module
   ==================================

   This module provides highly optimized Numba-compiled functions for interpolating
   gravitational lensing cross sections in a 5-dimensional parameter space.

   The module uses cubic B-spline interpolation on prefiltered coefficient grids to
   achieve efficient and accurate interpolation of cross section values based on:
     - Ellipticity components (e1, e2) derived from axis ratio q and position angle phi
     - Density profile slope (gamma)
     - External shear components (gamma1, gamma2)

   Key Features:
   -------------
   - JIT-compiled with Numba for high performance
   - Parallel processing support for batch evaluations
   - Cubic B-spline interpolation matching scipy's map_coordinates behavior
   - Automatic scaling by Einstein radius and affine calibration

   Author: Hemanta Kumar Phurailatpam

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.cross_section_interpolator.cross_section
   ler.lens_galaxy_population.cross_section_interpolator.make_cross_section_reinit



Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.cross_section_interpolator.C_LIGHT


.. py:data:: C_LIGHT
   :value: '299792.458'

   

.. py:function:: cross_section(zs, zl, sigma, q, phi, gamma, gamma1, gamma2, e1_grid, e2_grid, gamma_grid, gamma1_grid, gamma2_grid, cs_spline_coeff_grid, Da_instance, csunit_to_cs_slope, csunit_to_cs_intercept)

   
   theta_E, gamma, gamma1, gamma2, q, phi: 1D arrays (same length)
   grids: 1D
   cs_spline_coeff: spline-filtered coefficients (same shape as cs_spline_coeff)
   returns: 1D array
















   ..
       !! processed by numpydoc !!

.. py:function:: make_cross_section_reinit(e1_grid, e2_grid, gamma_grid, gamma1_grid, gamma2_grid, cs_spline_coeff_grid, Da_instance, csunit_to_cs_slope=0.31830988618379075, csunit_to_cs_intercept=-3.2311742677852644e-27)

   
















   ..
       !! processed by numpydoc !!

