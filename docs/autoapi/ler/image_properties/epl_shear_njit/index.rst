:py:mod:`ler.image_properties.epl_shear_njit`
=============================================

.. py:module:: ler.image_properties.epl_shear_njit

.. autoapi-nested-parse::

   Numba-JIT-compiled EPL + external shear lensing solver.

   Provides analytical image-position finding, magnification computation,
   Fermat-potential evaluation, and a parallel batch solver for
   elliptical power-law (EPL) lens models with external shear.

   The module uses the 1-D lens-equation approach, reducing the 2-D vector
   equation to a scalar root-finding problem parameterised by the
   image-plane angle. All computationally intensive functions are compiled
   with ``@njit`` (Numba) for performance.

   The top-level entry points are:

   - :func:`image_position_analytical_njit` — single-system solver

   - :func:`create_epl_shear_solver` — factory for a parallel batch solver


   Copyright (C) 2026 Author Name. Distributed under MIT License.

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.image_properties.epl_shear_njit.pol_to_ell
   ler.image_properties.epl_shear_njit.omega_scalar
   ler.image_properties.epl_shear_njit.lensing_diagnostics_scalar
   ler.image_properties.epl_shear_njit.fermat_potential_scalar
   ler.image_properties.epl_shear_njit.image_position_analytical_njit
   ler.image_properties.epl_shear_njit.create_epl_shear_solver



Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.image_properties.epl_shear_njit.EPS
   ler.image_properties.epl_shear_njit.MAX_ROOTS
   ler.image_properties.epl_shear_njit.MAX_IMGS
   ler.image_properties.epl_shear_njit.C_LIGHT


.. py:data:: EPS
   :value: '1e-16'

   

.. py:data:: MAX_ROOTS
   :value: '16'

   

.. py:data:: MAX_IMGS
   :value: '5'

   

.. py:data:: C_LIGHT
   :value: '299792458.0'

   

.. py:function:: pol_to_ell(r, theta, q)

   
   Convert polar coordinates to elliptical coordinates.


   :Parameters:

       **r** : ``float``
           Polar radial coordinate.

       **theta** : ``float``
           Polar angle in radians.

       **q** : ``float``
           Axis ratio of the ellipse (0 < q ≤ 1).

   :Returns:

       **rell** : ``float``
           Elliptical radial coordinate,
           ``rell = r * sqrt(q^2*cos^2(theta) + sin^2(theta))``.

       **phi** : ``float``
           Elliptical angle in radians,
           ``phi = arctan2(sin(theta), cos(theta)*q)``.










   .. rubric:: Examples

   >>> import numpy as np
   >>> from ler.image_properties.epl_shear_njit import pol_to_ell
   >>> rell, phi = pol_to_ell(1.0, np.pi / 4, 0.7)



   ..
       !! processed by numpydoc !!

.. py:function:: omega_scalar(phi, t, q, niter_max=200, tol=1e-16)

   
   Scalar series expansion of the EPL deflection kernel Omega.

   Evaluates the convergent series for the EPL Omega function at a
   single elliptical angle, avoiding temporary array allocation.
   The series converges geometrically with ratio
   ``f = (1-q)/(1+q)``.

   :Parameters:

       **phi** : ``float``
           Elliptical angle in radians.

       **t** : ``float``
           EPL slope parameter (``t = gamma - 1``).

       **q** : ``float``
           Axis ratio (0 < q ≤ 1).

       **niter_max** : ``int``
           Maximum number of series iterations.
           default: 200

       **tol** : ``float``
           Convergence tolerance for series truncation.
           default: 1e-16

   :Returns:

       **omega_sum** : ``complex``
           Complex EPL deflection kernel value at angle ``phi``.










   .. rubric:: Examples

   >>> import numpy as np
   >>> from ler.image_properties.epl_shear_njit import omega_scalar
   >>> omega = omega_scalar(np.pi / 3, t=1.0, q=0.8)



   ..
       !! processed by numpydoc !!

.. py:function:: lensing_diagnostics_scalar(z, b, t, gamma1, gamma2, q, phi, Omega)

   
   Compute magnification and image type at one image position.

   Evaluates the total Hessian of the lensing potential (EPL +
   external shear), then derives the Jacobian determinant and trace
   to classify the image and compute its signed magnification.
   Inputs ``z``, ``b``, ``t``, ``q``, ``phi``, and ``Omega`` are all
   scalars, not arrays.

   Image type convention:
   - 1: Type I  (minimum of the Fermat potential)
   - 2: Type II (saddle point)
   - 3: Type III (maximum of the Fermat potential)
   - 0: undefined / degenerate (on a critical curve)

   :Parameters:

       **z** : ``complex``
           Image position in the axis-aligned frame
           (``z = exp(-i*phi) * (x + i*y)``).

       **b** : ``float``
           EPL scale parameter (``b = theta_E * sqrt(q)``).

       **t** : ``float``
           EPL slope parameter (``t = gamma - 1``).

       **gamma1** : ``float``
           External shear component 1.

       **gamma2** : ``float``
           External shear component 2.

       **q** : ``float``
           Axis ratio (0 < q ≤ 1).

       **phi** : ``float``
           Lens position angle in radians.

       **Omega** : ``complex``
           EPL deflection kernel at the image elliptical angle.

   :Returns:

       **mu** : ``float``
           Signed magnification ``1/det(A)``; ``np.inf`` on a critical curve.

       **image_type** : ``int``
           Image-type code (0, 1, 2, or 3).










   .. rubric:: Examples

   >>> import numpy as np
   >>> from ler.image_properties.epl_shear_njit import (
   ...     lensing_diagnostics_scalar, omega_scalar
   ... )
   >>> z = np.exp(-1j * 0.0) * (0.8 + 1j * 0.3)
   >>> phi_ell = np.angle(z.real * 0.8 + 1j * z.imag)
   >>> Omega = omega_scalar(phi_ell, t=1.0, q=0.8)
   >>> mu, itype = lensing_diagnostics_scalar(
   ...     z, b=0.894, t=1.0, gamma1=0.0, gamma2=0.0, q=0.8, phi=0.0, Omega=Omega
   ... )



   ..
       !! processed by numpydoc !!

.. py:function:: fermat_potential_scalar(z, x, y, x_source, y_source, b, t, gamma1, gamma2, q, phi, Omega)

   
   Compute the Fermat potential at one image position.

   Returns the geometric minus gravitational time-delay contribution:
   ``tau = 0.5*|theta - beta|^2 - psi_EPL(theta) - psi_shear(theta)``

   :Parameters:

       **z** : ``complex``
           Image position in the axis-aligned frame
           (``z = exp(-i*phi) * (x + i*y)``).

       **x** : ``float``
           Image x-coordinate in sky frame.

       **y** : ``float``
           Image y-coordinate in sky frame.

       **x_source** : ``float``
           Source x-coordinate in sky frame.

       **y_source** : ``float``
           Source y-coordinate in sky frame.

       **b** : ``float``
           EPL scale parameter (``b = theta_E * sqrt(q)``).

       **t** : ``float``
           EPL slope parameter (``t = gamma - 1``).

       **gamma1** : ``float``
           External shear component 1.

       **gamma2** : ``float``
           External shear component 2.

       **q** : ``float``
           Axis ratio (0 < q ≤ 1).

       **phi** : ``float``
           Lens position angle in radians.

       **Omega** : ``complex``
           EPL deflection kernel at the image elliptical angle.

   :Returns:

       **tau** : ``float``
           Fermat potential (dimensionless; proportional to arrival-time delay).













   ..
       !! processed by numpydoc !!

.. py:function:: image_position_analytical_njit(x_src, y_src, q, phi, gamma, gamma1, gamma2, theta_E=1.0, alpha_scaling=1.0, magnification_limit=0.01, Nmeas=400, Nmeas_extra=80)

   
   Standalone EPL + external shear analytical image finder.

   Locates all lensed images for a given source position, computes
   signed magnifications, Fermat potentials (arrival-time proxies),
   and image types. Results are sorted by ascending arrival time and
   filtered by a minimum magnification threshold.

   :Parameters:

       **x_src** : ``float``
           Source x-coordinate (normalised to Einstein radius).

       **y_src** : ``float``
           Source y-coordinate (normalised to Einstein radius).

       **q** : ``float``
           Lens axis ratio (0 < q ≤ 1).

       **phi** : ``float``
           Lens position angle in radians.

       **gamma** : ``float``
           EPL power-law slope (``gamma = 2`` for isothermal).

       **gamma1** : ``float``
           External shear component 1.

       **gamma2** : ``float``
           External shear component 2.

       **theta_E** : ``float``
           Einstein radius.
           default: 1.0

       **alpha_scaling** : ``float``
           Deflection scaling applied as
           ``theta_E_eff = theta_E * alpha_scaling^(1/(gamma-1))``.
           default: 1.0

       **magnification_limit** : ``float``
           Minimum ``|mu|`` for an image to be retained.
           default: 0.01

       **Nmeas** : ``int``
           Number of uniformly-spaced angle samples for root finding.
           default: 400

       **Nmeas_extra** : ``int``
           Number of extra refinement samples near the source angle.
           default: 80

   :Returns:

       **x_img** : ``numpy.ndarray``
           Image x-coordinates.

       **y_img** : ``numpy.ndarray``
           Image y-coordinates.

       **fermat_pot** : ``numpy.ndarray``
           Fermat potentials at each image.

       **magnification** : ``numpy.ndarray``
           Signed magnifications at each image.

       **image_type** : ``numpy.ndarray``
           Image-type codes (1 = min, 2 = saddle, 3 = max, 0 = degenerate).

       **nimg** : ``int``
           Number of valid images retained.










   .. rubric:: Examples

   >>> from ler.image_properties.epl_shear_njit import image_position_analytical_njit
   >>> x_img, y_img, fermat_pot, mu, itype, nimg = image_position_analytical_njit(
   ...     x_src=0.1, y_src=0.05,
   ...     q=0.8, phi=0.3, gamma=2.0,
   ...     gamma1=0.05, gamma2=0.02,
   ... )
   >>> print(nimg)



   ..
       !! processed by numpydoc !!

.. py:function:: create_epl_shear_solver(arrival_time_sort=True, max_img=4, num_th=500, maginf=-100.0, alpha_scaling=1.0, magnification_limit=0.01, Nmeas=400, Nmeas_extra=80)

   
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

