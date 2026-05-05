:py:mod:`ler.image_properties.cross_section_njit copy`
======================================================

.. py:module:: ler.image_properties.cross_section_njit copy

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

   ler.image_properties.cross_section_njit copy.phi_q2_ellipticity
   ler.image_properties.cross_section_njit copy.ellipticity2phi_q
   ler.image_properties.cross_section_njit copy.pol_to_ell
   ler.image_properties.cross_section_njit copy.omega_scalar
   ler.image_properties.cross_section_njit copy.cdot
   ler.image_properties.cross_section_njit copy.solve_quad_scalar
   ler.image_properties.cross_section_njit copy.alpha_epl_shear_scalar
   ler.image_properties.cross_section_njit copy.select_physical_root
   ler.image_properties.cross_section_njit copy.caustics_epl_shear
   ler.image_properties.cross_section_njit copy.polygon_area
   ler.image_properties.cross_section_njit copy.make_cross_section_area_reinit
   ler.image_properties.cross_section_njit copy.cross_section_epl_shear_unit



Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.image_properties.cross_section_njit copy.C_LIGHT
   ler.image_properties.cross_section_njit copy.PI
   ler.image_properties.cross_section_njit copy.TWO_PI


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

.. py:function:: pol_to_ell(r, theta, q)

   
   Convert polar coordinates to elliptical coordinates.
   Use in epl_shear_njit.py


   :Parameters:

       **r** : ``float`` or ``numpy.ndarray``
           Radial coordinate.

       **theta** : ``float`` or ``numpy.ndarray``
           Polar angle in radians.

       **q** : ``float``
           Axis ratio.

   :Returns:

       **rell** : ``float`` or ``numpy.ndarray``
           Elliptical radial coordinate.

       **phi** : ``float`` or ``numpy.ndarray``
           Elliptical angle in radians.













   ..
       !! processed by numpydoc !!

.. py:function:: omega_scalar(phi, t, q, niter_max=200, tol=1e-16)

   
   Return Re(Omega), Im(Omega) for the EPL angular function.

   Uses a purely real recurrence to avoid complex arithmetic in the hot path.















   ..
       !! processed by numpydoc !!

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

.. py:function:: solve_quad_scalar(a, b, c)


.. py:function:: alpha_epl_shear_scalar(x, y, b, q, t, gamma1, gamma2, omega_real, omega_imag)

   
   Calculates total deflection components (real, imag) for a single point.
















   ..
       !! processed by numpydoc !!

.. py:function:: select_physical_root(u0, u1, prefer_second)

   
   Select a physically meaningful inverse-radius root.

   Priority:
   1) Root with small imaginary part and positive real part.
   2) If both valid, keep the analytically expected branch.
   3) If neither valid, fall back to preferred branch real part if positive,
      otherwise to the alternate branch if positive, else a small floor.















   ..
       !! processed by numpydoc !!

.. py:function:: caustics_epl_shear(q, phi, t, gamma1, gamma2, b, theta, maginf=-100.0)


.. py:function:: polygon_area(xv, yv)

   
   Compute the area of a simple polygon using the Shoelace formula.


   :Parameters:

       **xv** : ``numpy.ndarray``
           Polygon vertex x-coordinates.

       **yv** : ``numpy.ndarray``
           Polygon vertex y-coordinates.

   :Returns:

       **area** : ``float``
           The enclosed geometric area.










   .. rubric:: Examples

   >>> area = polygon_area(np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0]))



   ..
       !! processed by numpydoc !!

.. py:function:: make_cross_section_area_reinit(Da_instance, num_th=500, maginf=-100.0)

   
   Create a JIT-compiled cross-section evaluator for batched systems.
















   ..
       !! processed by numpydoc !!

.. py:function:: cross_section_epl_shear_unit(e1, e2, gamma, gamma1, gamma2, num_th=500, maginf=-100.0)

   
   Compute double-caustic cross sections for batched lens parameters.
















   ..
       !! processed by numpydoc !!

