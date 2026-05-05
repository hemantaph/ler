:py:mod:`ler.image_properties.cross_section_njit copy 2`
========================================================

.. py:module:: ler.image_properties.cross_section_njit copy 2

.. autoapi-nested-parse::

   Module for analytical caustic computation and cross-section evaluation
   of EPL (Elliptical Power-Law) + external shear lens models.

   Provides numba-accelerated routines for computing critical curves, caustic
   boundaries, polygon areas, and lensing cross sections. All core functions
   are decorated with ``@njit`` for high performance.

   Usage:
       Compute the double-caustic boundary for a single lens:

       >>> from ler.image_properties.cross_section_njit import caustics_epl_shear
       >>> pts = caustics_epl_shear(q=0.8, phi=0.0, gamma=2.0, gamma1=0.03, gamma2=-0.01)

   Copyright (C) 2026 Phurailatpam Hemantakumar. Distributed under MIT License.

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.image_properties.cross_section_njit copy 2.omega_scalar
   ler.image_properties.cross_section_njit copy 2.solve_quad_scalar_real
   ler.image_properties.cross_section_njit copy 2.partition_polar
   ler.image_properties.cross_section_njit copy 2.quicksort_polar
   ler.image_properties.cross_section_njit copy 2.caustic_area_epl_shear
   ler.image_properties.cross_section_njit copy 2.caustic_points_epl_shear
   ler.image_properties.cross_section_njit copy 2.make_cross_section_area_reinit
   ler.image_properties.cross_section_njit copy 2.pol_to_ell
   ler.image_properties.cross_section_njit copy 2.select_physical_root
   ler.image_properties.cross_section_njit copy 2.cdot
   ler.image_properties.cross_section_njit copy 2.cross_section_epl_shear_unit
   ler.image_properties.cross_section_njit copy 2.phi_q2_ellipticity
   ler.image_properties.cross_section_njit copy 2.ellipticity2phi_q



Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.image_properties.cross_section_njit copy 2.C_LIGHT
   ler.image_properties.cross_section_njit copy 2.PI
   ler.image_properties.cross_section_njit copy 2.TWO_PI


.. py:data:: C_LIGHT
   :value: '299792458.0'

   

.. py:data:: PI

   

.. py:data:: TWO_PI

   

.. py:function:: omega_scalar(phi, t, q, niter_max=200, tol=1e-16)


.. py:function:: solve_quad_scalar_real(a, b, c)


.. py:function:: partition_polar(th, r, low, high)


.. py:function:: quicksort_polar(th, r, low, high)


.. py:function:: caustic_area_epl_shear(q, t, gamma1, gamma2, b, num_th, cos_th, sin_th, cos_2th, sin_2th, rcut, thcut, r_main, th_main, r2, maginf=-100.0)

   
   Calculates only the area, skipping rotation and coordinate array construction.
















   ..
       !! processed by numpydoc !!

.. py:function:: caustic_points_epl_shear(theta_E, q, phi, gamma, gamma1, gamma2, num_th=500, maginf=-100.0)

   
   Calculates the 2D coordinates of the double caustic for a SINGLE lens.
   Accepts scalar float values.
















   ..
       !! processed by numpydoc !!

.. py:function:: make_cross_section_area_reinit(Da_instance, num_th=500, maginf=-100.0)

   
   Drop-in replacement with:
   - trig arrays closed over (faster)
   - safe memory layout
   - optional chunking for huge batches (no thread-ID hack)
















   ..
       !! processed by numpydoc !!

.. py:function:: pol_to_ell(r, theta, q)


.. py:function:: select_physical_root(u0, u1, prefer_second)


.. py:function:: cdot(a, b)

   
   Compute the real-valued dot product of two complex numbers.
   Used in epl_shear_njit.py.

   Equivalent to ``Re(a) * Re(b) + Im(a) * Im(b)``.

   :Parameters:

       **a** : ``complex``
           First complex number.

       **b** : ``complex``
           Second complex number.

   :Returns:

       **result** : ``float``
           Real-valued dot product.










   .. rubric:: Examples

   >>> cdot(1+2j, 3+4j)
   11.0



   ..
       !! processed by numpydoc !!

.. py:function:: cross_section_epl_shear_unit(e1, e2, gamma, gamma1, gamma2, num_th=500, maginf=-100.0)

   
   Compute double-caustic cross sections for batched lens parameters.
















   ..
       !! processed by numpydoc !!

.. py:function:: phi_q2_ellipticity(phi, q)

   
   Convert lens orientation and axis ratio to ellipticity components.


   :Parameters:

       **phi** : ``float``
           Position angle of the lens major axis in radians.

       **q** : ``float``
           Axis ratio (minor/major), where ``0 < q <= 1``.

   :Returns:

       **e1** : ``float``
           First ellipticity component.

       **e2** : ``float``
           Second ellipticity component.










   .. rubric:: Examples

   >>> e1, e2 = phi_q2_ellipticity(phi=0.25, q=0.8)



   ..
       !! processed by numpydoc !!

.. py:function:: ellipticity2phi_q(e1, e2)

   
   Convert complex ellipticity moduli to orientation angle and axis ratio.


   :Parameters:

       **e1** : ``float``
           Eccentricity in x-direction.

       **e2** : ``float``
           Eccentricity in xy-direction.

   :Returns:

       **phi** : ``float``
           Orientation angle in radians.

       **q** : ``float``
           Axis ratio (minor/major).













   ..
       !! processed by numpydoc !!

