:py:mod:`ler.image_properties.cross_section_njit test1`
=======================================================

.. py:module:: ler.image_properties.cross_section_njit test1

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

   ler.image_properties.cross_section_njit test1.phi_q2_ellipticity
   ler.image_properties.cross_section_njit test1.ellipticity2phi_q
   ler.image_properties.cross_section_njit test1.pol_to_ell
   ler.image_properties.cross_section_njit test1.omega
   ler.image_properties.cross_section_njit test1.omega_scalar
   ler.image_properties.cross_section_njit test1.cdot
   ler.image_properties.cross_section_njit test1.pol_to_cart
   ler.image_properties.cross_section_njit test1.cart_to_pol
   ler.image_properties.cross_section_njit test1.caustic_points_epl_shear
   ler.image_properties.cross_section_njit test1.caustic_area_epl_shear
   ler.image_properties.cross_section_njit test1.polygon_area
   ler.image_properties.cross_section_njit test1.make_cross_section_area_reinit
   ler.image_properties.cross_section_njit test1.cross_section_epl_shear_unit



Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.image_properties.cross_section_njit test1.C_LIGHT
   ler.image_properties.cross_section_njit test1.PI
   ler.image_properties.cross_section_njit test1.TWO_PI


.. py:data:: C_LIGHT
   :value: '299792458.0'

   

.. py:data:: PI

   

.. py:data:: TWO_PI

   

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
















   ..
       !! processed by numpydoc !!

.. py:function:: pol_to_ell(r, theta, q)

   
   Convert polar coordinates to elliptical coordinates.
















   ..
       !! processed by numpydoc !!

.. py:function:: omega(phi, t, q, niter_max=200, tol=1e-16)

   
   Evaluate the complex angular function Omega for the EPL profile.

   This series expansion converges geometrically with ratio
   ``f = (1 - q)/(1 + q)``. The ``fastmath`` flag provides ~4x speedup
   due to the reduction nature of the summation.

   :Parameters:

       **phi** : ``numpy.ndarray``
           Azimuthal angles in radians.

       **t** : ``float``
           EPL slope exponent (``t = gamma - 1``).

       **q** : ``float``
           Axis ratio.

       **niter_max** : ``int``
           Maximum number of series terms.

           default: 200

       **tol** : ``float``
           Convergence tolerance.

           default: 1e-16

   :Returns:

       **omegas** : ``numpy.ndarray``
           Complex Omega values at each angle (same object as input buffer).










   .. rubric:: Examples

   >>> import numpy as np
   >>> phi = np.linspace(0, 2 * np.pi, 100)
   >>> result = np.empty_like(phi, dtype=np.complex128)
   >>> omega(phi, t=1.0, q=0.8, omegas=result)



   ..
       !! processed by numpydoc !!

.. py:function:: omega_scalar(phi, t, q, niter_max=200, tol=1e-16)

   
   Scalar version of ``omega`` that avoids temporary array allocation.
















   ..
       !! processed by numpydoc !!

.. py:function:: cdot(a, b)

   
   Compute the real-valued dot product of two complex numbers.

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

.. py:function:: pol_to_cart(r, th)

   
   Convert polar coordinates to Cartesian coordinates.


   :Parameters:

       **r** : ``float`` or ``numpy.ndarray``
           Radial coordinate.

       **th** : ``float`` or ``numpy.ndarray``
           Polar angle in radians.

   :Returns:

       **x** : ``float`` or ``numpy.ndarray``
           Cartesian x-coordinate.

       **y** : ``float`` or ``numpy.ndarray``
           Cartesian y-coordinate.










   .. rubric:: Examples

   >>> x, y = pol_to_cart(1.0, np.pi / 4)



   ..
       !! processed by numpydoc !!

.. py:function:: cart_to_pol(x, y)

   
   Convert Cartesian coordinates to polar coordinates.

   The returned angle is wrapped to ``[0, 2π)``.

   :Parameters:

       **x** : ``float`` or ``numpy.ndarray``
           Cartesian x-coordinate.

       **y** : ``float`` or ``numpy.ndarray``
           Cartesian y-coordinate.

   :Returns:

       **r** : ``float`` or ``numpy.ndarray``
           Radial coordinate.

       **theta** : ``float`` or ``numpy.ndarray``
           Polar angle in radians, wrapped to ``[0, 2π)``.










   .. rubric:: Examples

   >>> r, theta = cart_to_pol(1.0, 1.0)



   ..
       !! processed by numpydoc !!

.. py:function:: caustic_points_epl_shear(theta_E, q, phi, gamma, gamma1, gamma2, num_th=500, maginf=-100.0)

   
   Calculates the 2D coordinates of the double caustic for a SINGLE lens.
   Accepts scalar float values.
















   ..
       !! processed by numpydoc !!

.. py:function:: caustic_area_epl_shear(q, phi, gamma, gamma1, gamma2, theta_E, theta, cos_th, sin_th, cos_2th, sin_2th, maginf=-100.0)

   
   Calculates the 2D coordinates of the double caustic for a SINGLE lens.
   Accepts scalar float values.
















   ..
       !! processed by numpydoc !!

.. py:function:: polygon_area(xv, yv)

   
   Compute the area of a simple polygon using the Shoelace formula.
















   ..
       !! processed by numpydoc !!

.. py:function:: make_cross_section_area_reinit(Da_instance, num_th=500, maginf=-100.0)

   
   Drop-in replacement with:
   - trig arrays closed over (faster)
   - safe memory layout
   - optional chunking for huge batches (no thread-ID hack)
















   ..
       !! processed by numpydoc !!

.. py:function:: cross_section_epl_shear_unit(e1, e2, gamma, gamma1, gamma2, num_th=500, maginf=-100.0)

   
   Compute double-caustic cross sections for batched lens parameters.
















   ..
       !! processed by numpydoc !!

