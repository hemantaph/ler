:py:mod:`ler.image_properties.epl_shear_njit_old`
=================================================

.. py:module:: ler.image_properties.epl_shear_njit_old

.. autoapi-nested-parse::

   Module for semi-analytical image finding in EPL + external shear lens models.

   Provides a numba-accelerated solver that locates multiple lensed images,
   computes their magnifications, arrival times, and image types for an
   elliptical power-law (EPL) mass profile with external shear.

   Usage:
       Solve for image positions of a single source:

       >>> from ler.image_properties.epl_shear_njit import image_position_analytical_njit
       >>> x_img, y_img, tau, mu, itype, n = image_position_analytical_njit(
       ...     x=0.1, y=0.05, q=0.8, phi=0.3, gamma=2.0,
       ...     gamma1=0.03, gamma2=-0.01
       ... )

   Copyright (C) 2026 Phurailatpam Hemantakumar. Distributed under MIT License.

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.image_properties.epl_shear_njit_old.lensing_diagnostics_scalar
   ler.image_properties.epl_shear_njit_old.fermat_potential_scalar
   ler.image_properties.epl_shear_njit_old.image_position_analytical_njit
   ler.image_properties.epl_shear_njit_old.create_epl_shear_solver



Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.image_properties.epl_shear_njit_old.EPS
   ler.image_properties.epl_shear_njit_old.MAX_ROOTS
   ler.image_properties.epl_shear_njit_old.MAX_IMGS
   ler.image_properties.epl_shear_njit_old.C_LIGHT


.. py:data:: EPS
   :value: '1e-14'

   

.. py:data:: MAX_ROOTS
   :value: '64'

   

.. py:data:: MAX_IMGS
   :value: '5'

   

.. py:data:: C_LIGHT
   :value: '299792458.0'

   

.. py:function:: lensing_diagnostics_scalar(x, y, theta_E, gamma, gamma1, gamma2, q, phi, center_x=0.0, center_y=0.0)

   
   Compute all Hessian-derived local diagnostics at one image position.


   :Parameters:

       **x** : ``float``
           Image-plane x-coordinate.

       **y** : ``float``
           Image-plane y-coordinate.

       **theta_E** : ``float``
           Einstein radius.

       **gamma** : ``float``
           EPL power-law slope.

       **gamma1** : ``float``
           External shear component 1.

       **gamma2** : ``float``
           External shear component 2.

       **q** : ``float``
           Axis ratio.

       **phi** : ``float``
           Lens position angle in radians.

       **center_x** : ``float``
           Lens center x-coordinate.

       **center_y** : ``float``
           Lens center y-coordinate.

   :Returns:

       **f_xx** : ``float``
           Hessian component d²ψ/dx².

       **f_xy** : ``float``
           Hessian component d²ψ/dxdy.

       **f_yx** : ``float``
           Hessian component d²ψ/dydx.

       **f_yy** : ``float``
           Hessian component d²ψ/dy².

       **detA** : ``float``
           Jacobian determinant.

       **traceA** : ``float``
           Jacobian trace.

       **mu** : ``float``
           Signed magnification.

       **image_type** : ``int``
           Image classification.

           1 = Type I (minimum), 2 = Type II (saddle),

           3 = Type III (maximum), 0 = undefined.










   .. rubric:: Examples

   >>> f_xx, f_xy, f_yx, f_yy, detA, traceA, mu, itype = lensing_diagnostics_scalar(
   ...     x=0.5, y=0.3, theta_E=1.0, gamma=2.0,
   ...     gamma1=0.03, gamma2=-0.01, q=0.8, phi=0.3
   ... )



   ..
       !! processed by numpydoc !!

.. py:function:: fermat_potential_scalar(x_image, y_image, x_source, y_source, theta_E, gamma, gamma1, gamma2, q, phi, center_x=0.0, center_y=0.0, ra_0=0.0, dec_0=0.0)

   
   Compute the Fermat potential (geometric delay minus lensing potential).


   :Parameters:

       **x_image** : ``float``
           Image-plane x-coordinate.

       **y_image** : ``float``
           Image-plane y-coordinate.

       **x_source** : ``float``
           Source-plane x-coordinate.

       **y_source** : ``float``
           Source-plane y-coordinate.

       **theta_E** : ``float``
           Einstein radius.

       **gamma** : ``float``
           EPL power-law slope.

       **gamma1** : ``float``
           External shear component 1.

       **gamma2** : ``float``
           External shear component 2.

       **q** : ``float``
           Axis ratio.

       **phi** : ``float``
           Lens position angle in radians.

       **center_x** : ``float``
           Lens center x-coordinate.

       **center_y** : ``float``
           Lens center y-coordinate.

       **ra_0** : ``float``
           Shear center x-coordinate.

       **dec_0** : ``float``
           Shear center y-coordinate.

   :Returns:

       **tau** : ``float``
           Fermat potential value.










   .. rubric:: Examples

   >>> tau = fermat_potential_scalar(
   ...     x_image=0.5, y_image=0.3, x_source=0.1, y_source=0.05,
   ...     theta_E=1.0, gamma=2.0, gamma1=0.03, gamma2=-0.01, q=0.8, phi=0.3
   ... )



   ..
       !! processed by numpydoc !!

.. py:function:: image_position_analytical_njit(x, y, q, phi, gamma, gamma1, gamma2, theta_E=1.0, alpha_scaling=1.0, magnification_limit=0.01, Nmeas=400, Nmeas_extra=80)

   
   Standalone EPL + shear analytical image finder.

   Locates lensed images for a given source position, computes their
   magnifications, Fermat potential (arrival time proxy), and image
   types. Results are sorted by ascending arrival time and filtered
   by a minimum magnification threshold.

   :Parameters:

       **x** : ``float``
           Source-plane x-coordinate.

       **y** : ``float``
           Source-plane y-coordinate.

       **q** : ``float``
           Lens axis ratio.

       **phi** : ``float``
           Lens position angle in radians.

       **gamma** : ``float``
           EPL power-law slope (lenstronomy convention).

           The Tessore exponent is ``t = gamma - 1``.

       **gamma1** : ``float``
           External shear component 1.

       **gamma2** : ``float``
           External shear component 2.

       **theta_E** : ``float``
           Einstein radius.

           default: 1.0

       **alpha_scaling** : ``float``
           Deflection scaling factor applied to theta_E.

           default: 1.0

       **magnification_limit** : ``float``
           Minimum |mu| threshold; images below this are discarded.

           default: 0.01

       **Nmeas** : ``int``
           Angular root-finding grid size.

           default: 400

       **Nmeas_extra** : ``int``
           Extra refinement points near the source angle.

           default: 80

   :Returns:

       **x_img** : ``numpy.ndarray``
           Image x-positions sorted by arrival time.

       **y_img** : ``numpy.ndarray``
           Image y-positions sorted by arrival time.

       **arrival_time** : ``numpy.ndarray``
           Fermat potential values in increasing order.

       **magnification** : ``numpy.ndarray``
           Signed magnifications.

       **image_type** : ``numpy.ndarray``
           Image classification codes (int64).

           1 = Type I (minimum), 2 = Type II (saddle),

           3 = Type III (maximum), 0 = undefined.

       **nimg** : ``int``
           Number of images returned.










   .. rubric:: Examples

   >>> x_img, y_img, tau, mu, itype, n = image_position_analytical_njit(
   ...     x=0.1, y=0.05, q=0.8, phi=0.3, gamma=2.0,
   ...     gamma1=0.03, gamma2=-0.01
   ... )



   ..
       !! processed by numpydoc !!

.. py:function:: create_epl_shear_solver(arrival_time_sort=True, max_img=4, num_th=500, maginf=-100.0, max_tries=100, alpha_scaling=1.0, magnification_limit=0.01, Nmeas=400, Nmeas_extra=80)

   
   Create a parallel EPL + shear solver for batched lens systems.

   Returns a JIT-compiled function that, for each system in a batch,
   samples a source from the double caustic and solves for image
   positions, magnifications, time delays, and image types.

   :Parameters:

       **arrival_time_sort** : ``bool``
           Whether to sort images by arrival time.

           default: True

       **max_img** : ``int``
           Maximum number of images to store per system.

           default: 4

       **num_th** : ``int``
           Angular samples for caustic construction.

           default: 500

       **maginf** : ``float``
           Magnification cutoff for caustic boundary.

           default: -100.0

       **max_tries** : ``int``
           Maximum rejection-sampling attempts per source.

           default: 100

       **alpha_scaling** : ``float``
           Deflection scaling factor.

           default: 1.0

       **magnification_limit** : ``float``
           Minimum |mu| threshold for image retention.

           default: 0.01

       **Nmeas** : ``int``
           Angular root-finding grid size.

           default: 400

       **Nmeas_extra** : ``int``
           Extra refinement points.

           default: 80

   :Returns:

       **solve_epl_shear_multithreaded** : ``callable``
           Parallel solver function with signature

           ``(theta_E, D_dt, q, phi, gamma, gamma1, gamma2)``

           returning a tuple of result arrays.










   .. rubric:: Examples

   >>> solver = create_epl_shear_solver()
   >>> results = solver(theta_E, D_dt, q, phi, gamma, gamma1, gamma2)



   ..
       !! processed by numpydoc !!

