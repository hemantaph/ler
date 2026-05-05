:py:mod:`ler.image_properties.cross_section_njit`
=================================================

.. py:module:: ler.image_properties.cross_section_njit

.. autoapi-nested-parse::

   Module for analytical caustic computation and cross-section evaluation.

   Provides numba-accelerated routines for computing critical curves, caustic
   boundaries, polygon areas, and lensing cross sections for EPL (Elliptical
   Power-Law) + external shear lens models. All core functions are decorated
   with ``@njit`` for high performance.

   Usage:
       Basic workflow example:

       >>> from ler.image_properties.cross_section_njit import caustic_points_epl_shear
       >>> pts = caustic_points_epl_shear(theta_E=1.0, q=0.8, phi=0.0, gamma=2.0, gamma1=0.03, gamma2=-0.01)

   Copyright (C) 2026 Phurailatpam Hemantakumar. Distributed under MIT License.

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.image_properties.cross_section_njit.phi_q2_ellipticity
   ler.image_properties.cross_section_njit.ellipticity2phi_q
   ler.image_properties.cross_section_njit.omega
   ler.image_properties.cross_section_njit.cdot
   ler.image_properties.cross_section_njit.pol_to_cart
   ler.image_properties.cross_section_njit.cart_to_pol
   ler.image_properties.cross_section_njit.caustic_points_epl_shear
   ler.image_properties.cross_section_njit.caustic_area_epl_shear
   ler.image_properties.cross_section_njit.polygon_area
   ler.image_properties.cross_section_njit.make_cross_section_area_reinit
   ler.image_properties.cross_section_njit.cross_section_epl_shear_unit



Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.image_properties.cross_section_njit.C_LIGHT
   ler.image_properties.cross_section_njit.PI
   ler.image_properties.cross_section_njit.TWO_PI
   ler.image_properties.cross_section_njit.EPS


.. py:data:: C_LIGHT
   :value: '299792458.0'

   

.. py:data:: PI

   

.. py:data:: TWO_PI

   

.. py:data:: EPS
   :value: '1e-16'

   

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
           First ellipticity component.

       **e2** : ``float``
           Second ellipticity component.

   :Returns:

       **phi** : ``float``
           Orientation angle in radians.

       **q** : ``float``
           Axis ratio (minor/major).










   .. rubric:: Examples

   >>> phi, q = ellipticity2phi_q(0.1, 0.05)



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
           Complex Omega values at each angle.










   .. rubric:: Examples

   >>> import numpy as np
   >>> phi = np.linspace(0, 2 * np.pi, 100)
   >>> omegas = omega(phi, t=1.0, q=0.8)



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


   :Parameters:

       **theta_E** : ``float``
           Einstein radius.

       **q** : ``float``
           Axis ratio.

       **phi** : ``float``
           Position angle of the lens major axis in radians.

       **gamma** : ``float``
           EPL slope exponent.

       **gamma1** : ``float``
           First shear component.

       **gamma2** : ``float``
           Second shear component.

       **num_th** : ``int``
           Number of angular samples.

           default: 500

       **maginf** : ``float``
           Magnification cut threshold.

           default: -100.0

   :Returns:

       **rotated** : ``numpy.ndarray``
           Shape ``(2, num_th)`` Cartesian coordinates of the double caustic.










   .. rubric:: Examples

   >>> pts = caustic_points_epl_shear(theta_E=1.0, q=0.8, phi=0.0, gamma=2.0, gamma1=0.03, gamma2=-0.01)



   ..
       !! processed by numpydoc !!

.. py:function:: caustic_area_epl_shear(q, phi, gamma, gamma1, gamma2, theta_E, theta, cos_th, sin_th, cos_2th, sin_2th, maginf=-100.0)

   
   Calculates the area of the double caustic for a SINGLE lens.
   Accepts scalar float values.


   :Parameters:

       **q** : ``float``
           Axis ratio.

       **phi** : ``float``
           Position angle of the lens major axis in radians.

       **gamma** : ``float``
           EPL slope exponent.

       **gamma1** : ``float``
           First shear component.

       **gamma2** : ``float``
           Second shear component.

       **theta_E** : ``float``
           Einstein radius.

       **theta** : ``numpy.ndarray``
           Array of angles.

       **cos_th** : ``numpy.ndarray``
           Array of cosine of angles.

       **sin_th** : ``numpy.ndarray``
           Array of sine of angles.

       **cos_2th** : ``numpy.ndarray``
           Array of cosine of 2*angles.

       **sin_2th** : ``numpy.ndarray``
           Array of sine of 2*angles.

       **maginf** : ``float``
           Magnification cut threshold.

           default: -100.0

   :Returns:

       **area** : ``float``
           Area of the caustic.










   .. rubric:: Examples

   >>> import numpy as np
   >>> theta = np.linspace(0, 2*np.pi, 500)
   >>> area = caustic_area_epl_shear(0.8, 0.0, 2.0, 0.03, -0.01, 1.0, theta, np.cos(theta), np.sin(theta), np.cos(2*theta), np.sin(2*theta))



   ..
       !! processed by numpydoc !!

.. py:function:: polygon_area(xv, yv)

   
   Compute the area of a simple polygon using the Shoelace formula.


   :Parameters:

       **xv** : ``numpy.ndarray``
           x-coordinates of the polygon vertices.

       **yv** : ``numpy.ndarray``
           y-coordinates of the polygon vertices.

   :Returns:

       **area** : ``float``
           Area of the polygon.










   .. rubric:: Examples

   >>> import numpy as np
   >>> area = polygon_area(np.array([0., 1., 0.]), np.array([0., 0., 1.]))



   ..
       !! processed by numpydoc !!

.. py:function:: make_cross_section_area_reinit(Da_instance, num_th=500, maginf=-100.0)

   
   Make a jitted function to compute double-caustic cross sections.

   Drop-in replacement with trig arrays closed over (faster),
   safe memory layout, and optional chunking for huge batches.

   :Parameters:

       **Da_instance** : ``callable``
           Angular diameter distance function.

       **num_th** : ``int``
           Number of points to sample.

           default: 500

       **maginf** : ``float``
           Magnification cut threshold.

           default: -100.0

   :Returns:

       **cross_section_area** : ``callable``
           Jitted function to compute cross section areas.










   .. rubric:: Examples

   >>> cross_section_fn = make_cross_section_area_reinit(Da_instance)



   ..
       !! processed by numpydoc !!

.. py:function:: cross_section_epl_shear_unit(e1, e2, gamma, gamma1, gamma2, num_th=500, maginf=-100.0)

   
   Compute double-caustic cross sections for batched lens parameters.


   :Parameters:

       **e1** : ``numpy.ndarray``
           First ellipticity component array.

       **e2** : ``numpy.ndarray``
           Second ellipticity component array.

       **gamma** : ``numpy.ndarray``
           EPL slope exponent array.

       **gamma1** : ``numpy.ndarray``
           First shear component array.

       **gamma2** : ``numpy.ndarray``
           Second shear component array.

       **num_th** : ``int``
           Number of points to sample.

           default: 500

       **maginf** : ``float``
           Magnification cut threshold.

           default: -100.0

   :Returns:

       **cross_section_area** : ``callable``
           Jitted function to compute cross section areas.










   .. rubric:: Examples

   >>> import numpy as np
   >>> calculate_area = cross_section_epl_shear_unit(np.array([0.1]), np.array([0.05]), np.array([2.0]), np.array([0.03]), np.array([-0.01]))



   ..
       !! processed by numpydoc !!

