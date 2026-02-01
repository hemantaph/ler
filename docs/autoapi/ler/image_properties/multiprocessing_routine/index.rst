:py:mod:`ler.image_properties.multiprocessing_routine`
======================================================

.. py:module:: ler.image_properties.multiprocessing_routine

.. autoapi-nested-parse::

   Module for solving lens equations using multiprocessing.

   This sub-module contains functions to solve the lens equation for a given set
   of lens parameters. The lens equation is solved using the analytical solver in
   lenstronomy. These functions are used in the multiprocessing routine within
   the ImageProperties class.

   Usage:
       Basic workflow example:

       >>> from ler.image_properties.multiprocessing_routine import solve_lens_equation
       >>> import numpy as np
       >>> from multiprocessing import Pool
       >>> lens_params = np.array([2, 0.02, -0.01, 1.9, 0.1, 0.09, 0.25, 0.94, 1e-6, 0, 'EPL_NUMBA', 'SHEAR'], dtype=object)
       >>> result = solve_lens_equation(lens_params)

   Copyright (C) 2026 Phurailatpam Hemanta Kumar. Distributed under MIT License.

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.image_properties.multiprocessing_routine.solve_lens_equation



Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.image_properties.multiprocessing_routine.MAX_RETRIES
   ler.image_properties.multiprocessing_routine.MIN_MAGNIFICATION


.. py:data:: MAX_RETRIES
   :value: '100'

   

.. py:data:: MIN_MAGNIFICATION
   :value: '0.01'

   

.. py:function:: solve_lens_equation(lens_parameters)

   
   Solve the lens equation to find image properties.

   Uses the analytical solver from lenstronomy to find image positions,
   magnifications, time delays, and hessian properties for strongly
   lensed sources. Source positions are sampled from within the caustic
   region to ensure multiple imaging.

   :Parameters:

       **lens_parameters** : ``numpy.ndarray``
           Array of lens configuration parameters with the following structure:

           - [0]: n_min_images - minimum number of images required

           - [1]: e1 - ellipticity component 1

           - [2]: e2 - ellipticity component 2

           - [3]: gamma - power-law slope of mass density

           - [4]: gamma1 - external shear component 1

           - [5]: gamma2 - external shear component 2

           - [6]: zl - lens redshift

           - [7]: zs - source redshift

           - [8]: einstein_radius - Einstein radius (units: radians)

           - [9]: iteration - iteration index for tracking

           - [10:]: lens_model_list - lens model names (e.g., 'EPL_NUMBA', 'SHEAR')

   :Returns:

       **x_source** : ``float``
           Source x-position (units: radians).

       **y_source** : ``float``
           Source y-position (units: radians).

       **x0_image_position** : ``numpy.ndarray``
           Image x-positions (units: radians).

       **x1_image_position** : ``numpy.ndarray``
           Image y-positions (units: radians).

       **magnifications** : ``numpy.ndarray``
           Magnification factors for each image.

       **time_delays** : ``numpy.ndarray``
           Time delays for each image (units: seconds).

       **nImages** : ``int``
           Number of images formed.

       **determinant** : ``numpy.ndarray``
           Determinant of the lensing Jacobian for each image.

       **trace** : ``numpy.ndarray``
           Trace of the lensing Jacobian for each image.

       **iteration** : ``int``
           Iteration index passed through for tracking.










   .. rubric:: Examples

   >>> from ler.image_properties.multiprocessing_routine import solve_lens_equation
   >>> import numpy as np
   >>> from multiprocessing import Pool
   >>> lens_parameters1 = np.array([2, 0.024, -0.016, 1.89, 0.10, 0.09, 0.25, 0.94, 2.5e-06, 0, 'EPL_NUMBA', 'SHEAR'], dtype=object)
   >>> lens_parameters2 = np.array([2, -0.040, -0.014, 2.00, 0.08, -0.01, 1.09, 2.55, 1.0e-06, 1, 'EPL_NUMBA', 'SHEAR'], dtype=object)
   >>> input_arguments = np.vstack((lens_parameters1, lens_parameters2))
   >>> with Pool(2) as p:
   ...     result = p.map(solve_lens_equation, input_arguments)
   >>> print(f"Number of images: {result[0][6]}")



   ..
       !! processed by numpydoc !!

