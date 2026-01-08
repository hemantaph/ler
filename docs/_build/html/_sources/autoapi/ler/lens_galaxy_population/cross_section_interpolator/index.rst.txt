:py:mod:`ler.lens_galaxy_population.cross_section_interpolator`
===============================================================

.. py:module:: ler.lens_galaxy_population.cross_section_interpolator

.. autoapi-nested-parse::

   Module for gravitational lensing cross section interpolation.

   This module provides highly optimized Numba-compiled functions for interpolating
   gravitational lensing cross sections in a 5-dimensional parameter space using
   cubic B-spline interpolation on prefiltered coefficient grids.

   The interpolation is performed based on:

   - Ellipticity components (e1, e2) derived from axis ratio q and position angle phi

   - Density profile slope (gamma)

   - External shear components (gamma1, gamma2)


   Key Features:

   - JIT-compiled with Numba for high performance

   - Parallel processing support for batch evaluations

   - Cubic B-spline interpolation matching scipy's map_coordinates behavior

   - Automatic scaling by Einstein radius and affine calibration


   Usage:
       Basic workflow example:

       >>> from ler.lens_galaxy_population.cross_section_interpolator import make_cross_section_reinit
       >>> cs_func = make_cross_section_reinit(e1_grid, e2_grid, gamma_grid, ...)
       >>> cross_sections = cs_func(zs, zl, sigma, q, phi, gamma, gamma1, gamma2)

   Copyright (C) 2024 Hemanta Kumar Phurailatpam. Distributed under MIT License.

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.cross_section_interpolator.make_cross_section_reinit



Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.cross_section_interpolator.C_LIGHT


.. py:data:: C_LIGHT
   :value: '299792.458'

   

.. py:function:: make_cross_section_reinit(e1_grid, e2_grid, gamma_grid, gamma1_grid, gamma2_grid, cs_grid, Da_instance, csunit_to_cs_slope=0.31830988618379075, csunit_to_cs_intercept=-3.2311742677852644e-27)

   
   Factory function to create a JIT-compiled cross section calculator.

   This function precomputes B-spline coefficients and creates a closure
   that captures the grid parameters, returning a fast Numba-compiled
   function for computing cross sections.

   :Parameters:

       **e1_grid** : ``numpy.ndarray``
           Grid values for ellipticity component e1, shape (n_e1,).

       **e2_grid** : ``numpy.ndarray``
           Grid values for ellipticity component e2, shape (n_e2,).

       **gamma_grid** : ``numpy.ndarray``
           Grid values for density slope gamma, shape (n_g,).

       **gamma1_grid** : ``numpy.ndarray``
           Grid values for shear component gamma1, shape (n_g1,).

       **gamma2_grid** : ``numpy.ndarray``
           Grid values for shear component gamma2, shape (n_g2,).

       **cs_grid** : ``numpy.ndarray``
           Raw cross section grid data (before spline filtering),

           shape (n_e1, n_e2, n_g, n_g1, n_g2).

       **Da_instance** : ``callable``
           Angular diameter distance function.

           Signature: ``Da_instance(z) -> distance``

       **csunit_to_cs_slope** : ``float``
           Slope for affine calibration from unit cross section.

           default: 0.31830988618379075

       **csunit_to_cs_intercept** : ``float``
           Intercept for affine calibration from unit cross section.

           default: -3.2311742677852644e-27

   :Returns:

       **cross_section_reinit** : ``callable``
           JIT-compiled function with signature:

           ``cross_section_reinit(zs, zl, sigma, q, phi, gamma, gamma1, gamma2)``

           Returns cross sections as ``numpy.ndarray`` of shape (N,).










   .. rubric:: Examples

   >>> from ler.lens_galaxy_population.cross_section_interpolator import make_cross_section_reinit
   >>> cs_func = make_cross_section_reinit(
   ...     e1_grid, e2_grid, gamma_grid, gamma1_grid, gamma2_grid,
   ...     cs_grid, Da_instance
   ... )
   >>> cross_sections = cs_func(zs, zl, sigma, q, phi, gamma, gamma1, gamma2)



   ..
       !! processed by numpydoc !!

