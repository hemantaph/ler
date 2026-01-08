:py:mod:`ler.image_properties`
==============================

.. py:module:: ler.image_properties


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   image_properties/index.rst
   multiprocessing_routine/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   ler.image_properties.ImageProperties



Functions
~~~~~~~~~

.. autoapisummary::

   ler.image_properties.solve_lens_equation
   ler.image_properties.phi_q2_ellipticity_hemanta
   ler.image_properties.solve_lens_equation



Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.image_properties.cosmo
   ler.image_properties.MAX_RETRIES
   ler.image_properties.MIN_MAGNIFICATION


.. py:data:: cosmo

   

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

           - [8]: einstein_radius - Einstein radius (units: arcsec)

           - [9]: iteration - iteration index for tracking

           - [10:]: lens_model_list - lens model names (e.g., 'EPL_NUMBA', 'SHEAR')

   :Returns:

       **x_source** : ``float``
           Source x-position (units: arcsec).

       **y_source** : ``float``
           Source y-position (units: arcsec).

       **x0_image_position** : ``numpy.ndarray``
           Image x-positions (units: arcsec).

       **x1_image_position** : ``numpy.ndarray``
           Image y-positions (units: arcsec).

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

.. py:function:: phi_q2_ellipticity_hemanta(phi, q)

   
   Function to convert phi and q to ellipticity e1 and e2.


   :Parameters:

       **phi** : `float: array`
           angle of the major axis in radians

       **q** : `float: array`
           axis ratio

   :Returns:

       **e1** : `float: array`
           ellipticity component 1

       **e2** : `float: array`
           ..













   ..
       !! processed by numpydoc !!

.. py:class:: ImageProperties(npool=4, n_min_images=2, n_max_images=4, time_window=365 * 24 * 3600 * 20, lens_model_list=['EPL_NUMBA', 'SHEAR'], cosmology=None, spin_zero=True, spin_precession=False)


   
   Class to compute image properties of strongly lensed gravitational wave events.

   This class solves the lens equation to find image positions, magnifications,
   time delays, and image types (morse phase) for strongly lensed sources. It uses
   multiprocessing for efficient computation of large samples.

   Key Features:

   - Solves lens equations using multiprocessing for efficiency

   - Computes image positions, magnifications, and time delays

   - Classifies image types using morse phase

   - Calculates detection probabilities for lensed images

   :Parameters:

       **npool** : ``int``
           Number of processes for multiprocessing.

           default: 4

       **n_min_images** : ``int``
           Minimum number of images required for a valid lensing event.

           default: 2

       **n_max_images** : ``int``
           Maximum number of images to consider per event.

           default: 4

       **time_window** : ``float``
           Time window for lensed events (units: seconds).

           default: 365*24*3600*20 (20 years)

       **lens_model_list** : ``list``
           List of lens models to use.

           default: ['EPL_NUMBA', 'SHEAR']

       **cosmology** : ``astropy.cosmology`` or ``None``
           Cosmology for distance calculations.

           If None, uses default LambdaCDM.

           default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

       **spin_zero** : ``bool``
           If True, spin parameters are set to zero (no spin sampling).

           default: False

       **spin_precession** : ``bool``
           If True (and spin_zero=False), sample precessing spin parameters.

           If False (and spin_zero=False), sample aligned/anti-aligned spins.

           default: False











   .. rubric:: Examples

   Basic usage:

   >>> from ler.image_properties import ImageProperties
   >>> ip = ImageProperties()
   >>> lens_parameters = dict(
   ...     zs=np.array([2.0]),
   ...     zl=np.array([0.5]),
   ...     gamma1=np.array([0.0]),
   ...     gamma2=np.array([0.0]),
   ...     phi=np.array([0.0]),
   ...     q=np.array([0.8]),
   ...     gamma=np.array([2.0]),
   ...     theta_E=np.array([1.0])
   ... )
   >>> result = ip.image_properties(lens_parameters)
   >>> print(result.keys())

   Instance Methods
   ----------
   ImageProperties has the following methods:

   +-----------------------------------------------------+------------------------------------------------+
   | Method                                              | Description                                    |
   +=====================================================+================================================+
   | :meth:`~image_properties`                           | Compute image properties for lensed events     |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~get_lensed_snrs`                            | Compute detection probability for lensed images|
   +-----------------------------------------------------+------------------------------------------------+

   Instance Attributes
   ----------
   ImageProperties has the following attributes:

   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
   | Attribute                                           | Type                      | Unit     | Description                                    |
   +=====================================================+===========================+==========+================================================+
   | :attr:`~npool`                                      | ``int``                   |          | Number of multiprocessing workers              |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
   | :attr:`~n_min_images`                               | ``int``                   |          | Minimum number of images required              |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
   | :attr:`~n_max_images`                               | ``int``                   |          | Maximum number of images per event             |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
   | :attr:`~time_window`                                | ``float``                 | s        | Time window for lensed events                  |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
   | :attr:`~lens_model_list`                            | ``list``                  |          | List of lens models                            |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
   | :attr:`~cosmo`                                      | ``astropy.cosmology``     |          | Cosmology for calculations                     |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
   | :attr:`~spin_zero`                                  | ``bool``                  |          | Flag for zero spin assumption                  |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
   | :attr:`~spin_precession`                            | ``bool``                  |          | Flag for spin precession                       |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:property:: npool

      
      Number of multiprocessing workers.



      :Returns:

          **npool** : ``int``
              Number of processes for multiprocessing.

              default: 4













      ..
          !! processed by numpydoc !!

   .. py:property:: n_min_images

      
      Minimum number of images required for a valid lensing event.



      :Returns:

          **n_min_images** : ``int``
              Minimum number of images required.

              default: 2













      ..
          !! processed by numpydoc !!

   .. py:property:: n_max_images

      
      Maximum number of images per event.



      :Returns:

          **n_max_images** : ``int``
              Maximum number of images to consider per event.

              default: 4













      ..
          !! processed by numpydoc !!

   .. py:property:: lens_model_list

      
      List of lens models to use.



      :Returns:

          **lens_model_list** : ``list``
              List of lens model names.

              default: ['EPL_NUMBA', 'SHEAR']













      ..
          !! processed by numpydoc !!

   .. py:property:: spin_zero

      
      Flag for zero spin assumption.



      :Returns:

          **spin_zero** : ``bool``
              Whether to assume zero spin for compact objects.

              If True, spin parameters are set to zero (no spin sampling).

              If False, spin parameters are sampled.

              default: False













      ..
          !! processed by numpydoc !!

   .. py:property:: spin_precession

      
      Flag for spin precession.



      :Returns:

          **spin_precession** : ``bool``
              Whether to include spin precession effects.

              If True (and spin_zero=False), sample precessing spin parameters.

              If False (and spin_zero=False), sample aligned/anti-aligned spins.

              default: False













      ..
          !! processed by numpydoc !!

   .. py:property:: time_window

      
      Time window for lensed events.



      :Returns:

          **time_window** : ``float``
              Time window for lensed events (units: s).

              default: 365*24*3600*20 (20 years)













      ..
          !! processed by numpydoc !!

   .. py:property:: cosmo

      
      Astropy cosmology object for calculations.



      :Returns:

          **cosmo** : ``astropy.cosmology``
              Cosmology used for distance calculations.

              default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)













      ..
          !! processed by numpydoc !!

   .. py:method:: image_properties(lens_parameters)

      
      Compute image properties for strongly lensed events.

      Solves the lens equation using multiprocessing to find image positions,
      magnifications, time delays, and image types for each lensing event.

      :Parameters:

          **lens_parameters** : ``dict``
              Dictionary containing lens and source parameters with keys:

              - 'zs': source redshift (array)

              - 'zl': lens redshift (array)

              - 'gamma1': external shear component 1 (array)

              - 'gamma2': external shear component 2 (array)

              - 'phi': position angle of lens ellipticity (array)

              - 'q': axis ratio of lens (array)

              - 'gamma': power-law slope of mass density (array)

              - 'theta_E': Einstein radius in radians (array)

      :Returns:

          **lens_parameters** : ``dict``
              Updated dictionary with additional image properties:

              - 'x0_image_positions': x-coordinates of images (shape: size x n_max_images)

              - 'x1_image_positions': y-coordinates of images (shape: size x n_max_images)

              - 'magnifications': magnification factors (shape: size x n_max_images)

              - 'time_delays': time delays relative to first image (shape: size x n_max_images, units: s)

              - 'image_type': morse phase classification (1=minimum, 2=saddle, 3=maximum)

              - 'n_images': number of images per event (array)

              - 'x_source': source x-position (array)

              - 'y_source': source y-position (array)













      ..
          !! processed by numpydoc !!

   .. py:method:: get_lensed_snrs(lensed_param, pdet_calculator, list_of_detectors=None)

      
      Compute detection probability for each lensed image.

      Calculates the effective luminosity distance, geocent time, and phase
      for each image accounting for magnification and morse phase, then
      computes detection probabilities using the provided calculator.

      :Parameters:

          **lensed_param** : ``dict``
              Dictionary containing lensed source and image parameters with keys:

              - 'mass_1', 'mass_2': detector-frame masses (array)

              - 'luminosity_distance' or 'effective_luminosity_distance': distance (array)

              - 'geocent_time' or 'effective_geocent_time': GPS time (array)

              - 'phase' or 'effective_phase': coalescence phase (array)

              - 'theta_jn', 'psi', 'ra', 'dec': orientation and position (arrays)

              - 'magnifications': image magnifications (shape: size x n_max_images)

              - 'time_delays': image time delays (shape: size x n_max_images)

              - 'image_type': morse phase type (shape: size x n_max_images)

          **pdet_calculator** : ``callable``
              Function that computes detection probability given GW parameters.

          **list_of_detectors** : ``list`` or ``None``
              List of detector names (e.g., ['H1', 'L1', 'V1']) for per-detector results.

              default: None

      :Returns:

          **result_dict** : ``dict``
              Dictionary containing:

              - 'pdet_net': network detection probability (shape: size x n_max_images)

              - Individual detector probabilities if list_of_detectors provided

          **lensed_param** : ``dict``
              Updated dictionary with effective parameters:

              - 'effective_luminosity_distance': magnification-corrected distance

              - 'effective_geocent_time': time-delay-corrected GPS time

              - 'effective_phase': morse-phase-corrected coalescence phase













      ..
          !! processed by numpydoc !!


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

           - [8]: einstein_radius - Einstein radius (units: arcsec)

           - [9]: iteration - iteration index for tracking

           - [10:]: lens_model_list - lens model names (e.g., 'EPL_NUMBA', 'SHEAR')

   :Returns:

       **x_source** : ``float``
           Source x-position (units: arcsec).

       **y_source** : ``float``
           Source y-position (units: arcsec).

       **x0_image_position** : ``numpy.ndarray``
           Image x-positions (units: arcsec).

       **x1_image_position** : ``numpy.ndarray``
           Image y-positions (units: arcsec).

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

