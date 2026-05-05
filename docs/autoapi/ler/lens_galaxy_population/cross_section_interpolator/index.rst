:py:mod:`ler.lens_galaxy_population.cross_section_interpolator`
===============================================================

.. py:module:: ler.lens_galaxy_population.cross_section_interpolator

.. autoapi-nested-parse::

   Module for interpolating gravitational-lensing cross sections.

   This module provides fast Numba-compiled interpolation in a 5D parameter
   space for lensing cross sections. The interpolated dimensions are ellipticity
   (``e1``, ``e2``), density slope (``gamma``), and external shear
   (``gamma1``, ``gamma2``). Inputs from lens parameters are mapped to this grid,
   then rescaled by Einstein-radius geometry and affine calibration.

   Key Features:

   - 5D local interpolation using precomputed tensor basis weights

   - Zero memory allocation inside execution loops for maximum throughput

   - Numba JIT compilation with fused parallel loops

   - Direct conversion from (q, phi) to (e1, e2) before interpolation


   Copyright (C) 2024 Hemantakumar Phurailatpam. Distributed under MIT License.

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.cross_section_interpolator.make_cross_section_area_reinit



Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.cross_section_interpolator.C_LIGHT


.. py:data:: C_LIGHT
   :value: '299792.458'

   

.. py:function:: make_cross_section_area_reinit(e1_grid, e2_grid, gamma_grid, gamma1_grid, gamma2_grid, cs_unit_grid, Da_instance, csunit_to_cs_slope=0.31830988618379075, csunit_to_cs_intercept=-3.2311742677852644e-27)


