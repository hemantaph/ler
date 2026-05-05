:py:mod:`ler.image_properties`
==============================

.. py:module:: ler.image_properties


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   cross_section_njit/index.rst
   cross_section_njit copy/index.rst
   cross_section_njit copy 2/index.rst
   cross_section_njit test1/index.rst
   epl_shear_njit/index.rst
   epl_shear_njit copy/index.rst
   epl_shear_njit copy 2/index.rst
   epl_shear_njit_old/index.rst
   image_properties/index.rst
   multiprocessing_routine_epl_shear/index.rst
   sample_caustic_points_njit/index.rst


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
   ler.image_properties.phi_q2_ellipticity
   ler.image_properties.solve_lens_equation
   ler.image_properties.sample_source_from_double_caustic
   ler.image_properties.cdot
   ler.image_properties.pol_to_ell
   ler.image_properties.omega_scalar
   ler.image_properties.lensing_diagnostics_scalar
   ler.image_properties.fermat_potential_scalar
   ler.image_properties.image_position_analytical_njit
   ler.image_properties.create_epl_shear_solver
   ler.image_properties.phi_q2_ellipticity
   ler.image_properties.ellipticity2phi_q
   ler.image_properties.omega
   ler.image_properties.cdot
   ler.image_properties.pol_to_cart
   ler.image_properties.cart_to_pol
   ler.image_properties.caustic_points_epl_shear
   ler.image_properties.caustic_area_epl_shear
   ler.image_properties.polygon_area
   ler.image_properties.make_cross_section_area_reinit
   ler.image_properties.cross_section_epl_shear_unit
   ler.image_properties.caustic_points_epl_shear
   ler.image_properties.sample_source_from_double_caustic



Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.image_properties.Mpc_to_m
   ler.image_properties.MAX_RETRIES
   ler.image_properties.MIN_MAGNIFICATION
   ler.image_properties.EPS
   ler.image_properties.MAX_ROOTS
   ler.image_properties.MAX_IMGS
   ler.image_properties.C_LIGHT
   ler.image_properties.C_LIGHT
   ler.image_properties.PI
   ler.image_properties.TWO_PI
   ler.image_properties.EPS


.. py:function:: solve_lens_equation(lens_parameters)

   
   Solve the lens equation to find image properties.

   Uses the analytical solver from lenstronomy to find image positions,
   magnifications, time delays, and hessian properties for strongly
   lensed sources. Source positions are sampled from within the caustic
   region to ensure multiple imaging.

   :Parameters:

       **lens_parameters** : ``numpy.ndarray``
           Array of lens configuration parameters with the following structure:

           - [0]: e1 - ellipticity component 1

           - [1]: e2 - ellipticity component 2

           - [2]: gamma - power-law slope of mass density

           - [3]: gamma1 - external shear component 1

           - [4]: gamma2 - external shear component 2

           - [5]: zl - lens redshift

           - [6]: zs - source redshift

           - [7]: einstein_radius - Einstein radius (units: radians)

           - [8]: iteration - iteration index for tracking

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

   >>> from ler.image_properties.multiprocessing_routine import solve_lens_equation, _init_worker_multiprocessing
   >>> import numpy as np
   >>> from multiprocessing import Pool
   >>> lens_parameters1 = np.array([0.024, -0.016, 1.89, 0.10, 0.09, 0.25, 0.94, 2.5e-06, 0])
   >>> lens_parameters2 = np.array([-0.040, -0.014, 2.00, 0.08, -0.01, 1.09, 2.55, 1.0e-06, 1])
   >>> input_arguments = np.vstack((lens_parameters1, lens_parameters2))
   >>> with Pool(
   ...     processes=2, # Number of worker processes
   ...     initializer=_init_worker_multiprocessing, # common
   ...     initargs=(
   ...         2, # n_min_images
   ...         ['EPL_NUMBA', 'SHEAR'], # lensModelList
   ...     ),
   ... ) as pool:
   ...     result = pool.map(solve_lens_equation, input_arguments)



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

.. py:data:: Mpc_to_m
   :value: '3.085677581491367e+22'

   

.. py:class:: ImageProperties(npool=4, n_min_images=2, n_max_images=4, lens_model_list=['EPL_NUMBA', 'SHEAR'], image_properties_function='image_properties_epl_shear_njit', image_properties_function_params=None, cosmology=None, time_window=365.0 * 24.0 * 3600.0 * 1.0, spin_zero=True, spin_precession=False, pdet_finder=None, include_effective_parameters=False, multiprocessing_verbose=True, include_redundant_parameters=False)


   
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
           Time window for lensed events from min(geocent_time) (units: seconds).

           default: 365*24*3600*2 (2 years)

       **include_effective_parameters** : ``bool``
           Whether to include effective parameters (effective_phase, effective_ra, effective_dec) in the output.

           default: True

       **lens_model_list** : ``list``
           List of lens models to use.

           default: ['EPL_NUMBA', 'SHEAR']

       **cosmology** : ``astropy.cosmology`` or ``None``
           Cosmology for distance calculations.

           If None, uses default LambdaCDM.

           default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)

       **spin_zero** : ``bool``
           If True, spin parameters are set to zero (no spin sampling).

           default: False

       **spin_precession** : ``bool``
           If True (and spin_zero=False), sample precessing spin parameters.

           If False (and spin_zero=False), sample aligned/anti-aligned spins.

           default: False

       **multiprocessing_verbose** : ``bool``
           If True, shows a progress bar for multiprocessing tasks.

           default: True

       **include_redundant_parameters** : ``bool``
           If True, removes redundant parameters (e.g., theta_E, n_images, mass_1, mass_2, luminosity_distance) from output to save memory.











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
   | :meth:`~image_properties_epl_shear`                 | Compute image properties for lensed events     |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~get_lensed_snrs`                            | Compute detection probability for lensed images|
   +-----------------------------------------------------+------------------------------------------------+

   Instance Attributes
   ----------
   ImageProperties has the following attributes:

   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------------------------+
   | Attribute                                           | Type                      | Unit     | Description                                                      |
   +=====================================================+===========================+==========+==================================================================+
   | :attr:`~npool`                                      | ``int``                   |          | Number of multiprocessing workers                                |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------------------------+
   | :attr:`~multiprocessing_verbose`                    | ``bool``                  |          | If True, shows a progress bar for multiprocessing tasks          |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------------------------+
   | :attr:`~n_min_images`                               | ``int``                   |          | Minimum number of images required                                |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------------------------+
   | :attr:`~n_max_images`                               | ``int``                   |          | Maximum number of images per event                               |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------------------------+
   | :attr:`~time_window`                                | ``float``                 | s        | Time window for lensed events                                    |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------------------------+
   | :attr:`~include_effective_parameters`               | ``bool``                  |          | To include effective parameters in output                        |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------------------------+
   | :attr:`~include_redundant_parameters`               | ``bool``                  |          | If True, removes redundant parameters from output to save memory |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------------------------+
   | :attr:`~lens_model_list`                            | ``list``                  |          | List of lens models                                              |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------------------------+
   | :attr:`~cosmo`                                      | ``astropy.cosmology``     |          | Cosmology for calculations                                       |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------------------------+
   | :attr:`~spin_zero`                                  | ``bool``                  |          | Flag for zero spin assumption                                    |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------------------------+
   | :attr:`~spin_precession`                            | ``bool``                  |          | Flag for spin precession                                         |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------------------------+
   | :attr:`~pdet_finder`                                | ``callable``              |          | Probability of detection calculator                              |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------------------------+
   | :attr:`~pdet_finder_output_keys`                    | ``list``                  |          | Keys for probability of detection outputs                        |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------------------------+



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

              default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)













      ..
          !! processed by numpydoc !!

   .. py:property:: pdet_finder

      
      Detection probability finder function.



      :Returns:

          **pdet_finder** : ``callable``
              Function that calculates detection probability for GW events.

              The function signature should be:

              ``pdet_finder(gw_param_dict) -> dict`` with key 'pdet_net'.













      ..
          !! processed by numpydoc !!

   .. py:property:: pdet_finder_output_keys

      
      Output keys from the detection probability finder function.



      :Returns:

          **pdet_finder_output_keys** : ``list``
              List of keys returned by the pdet_finder function.

              default: None













      ..
          !! processed by numpydoc !!

   .. py:property:: include_effective_parameters

      
      Flag to include effective parameters in output.



      :Returns:

          **include_effective_parameters** : ``bool``
              Whether to include effective parameters in the output of get_lensed_snrs.

              default: False













      ..
          !! processed by numpydoc !!

   .. py:property:: available_image_properties_functions

      
      Dictionary of available functions for computing image properties.



      :Returns:

          **available_image_properties_functions** : ``dict``
              Dictionary with function names and default parameters.













      ..
          !! processed by numpydoc !!

   .. py:attribute:: image_properties_function

      

   .. py:attribute:: multiprocessing_verbose
      :value: 'True'

      

   .. py:attribute:: include_redundant_parameters
      :value: 'False'

      

   .. py:attribute:: image_properties_function_params

      

   .. py:method:: image_properties_epl_shear_njit(lens_parameters)

      
      Compute image properties for strongly lensed events. This use functions similar to lenstronomy but rewritten in numba njit for speed.

      Solves the lens equation using multiprocessing to find image positions,
      magnifications, time delays, and image types for each lensing event.

      :Parameters:

          **lens_parameters** : ``dict``
              Dictionary containing lens and source parameters shown in the table:

              +------------------------------+-----------+-------------------------------------------------------+
              | Parameter                    | Units     | Description                                           |
              +==============================+===========+=======================================================+
              | zl                           |           | redshift of the lens                                  |
              +------------------------------+-----------+-------------------------------------------------------+
              | zs                           |           | redshift of the source                                |
              +------------------------------+-----------+-------------------------------------------------------+
              | sigma                        | km s^-1   | velocity dispersion                                   |
              +------------------------------+-----------+-------------------------------------------------------+
              | q                            |           | axis ratio                                            |
              +------------------------------+-----------+-------------------------------------------------------+
              | theta_E                      | radian    | Einstein radius                                       |
              +------------------------------+-----------+-------------------------------------------------------+
              | phi                          | rad       | axis rotation angle. counter-clockwise from the       |
              |                              |           | positive x-axis (RA-like axis) to the major axis of   |
              |                              |           | the projected mass distribution.                      |
              +------------------------------+-----------+-------------------------------------------------------+
              | gamma                        |           | density profile slope of EPL galaxy                   |
              +------------------------------+-----------+-------------------------------------------------------+
              | gamma1                       |           | external shear component in the x-direction           |
              |                              |           | (RA-like axis)                                        |
              +------------------------------+-----------+-------------------------------------------------------+
              | gamma2                       |           | external shear component in the y-direction           |
              |                              |           | (Dec-like axis)                                       |
              +------------------------------+-----------+-------------------------------------------------------+

      :Returns:

          **lens_parameters** : ``dict``
              Updated dictionary with additional image properties with the following description:

              +------------------------------+-----------+-------------------------------------------------------+
              | x0_image_positions           | radian    | x-coordinate (RA-like axis) of the images             |
              +------------------------------+-----------+-------------------------------------------------------+
              | x1_image_positions           | radian    | y-coordinate (Dec-like axis) of the images            |
              +------------------------------+-----------+-------------------------------------------------------+
              | magnifications               |           | magnifications                                        |
              +------------------------------+-----------+-------------------------------------------------------+
              | time_delays                  |           | time delays                                           |
              +------------------------------+-----------+-------------------------------------------------------+
              | image_type                   |           | image type                                            |
              +------------------------------+-----------+-------------------------------------------------------+
              | n_images                     |           | number of images                                      |
              +------------------------------+-----------+-------------------------------------------------------+
              | x_source                     | radian    | x-coordinate (RA-like axis) of the source             |
              +------------------------------+-----------+-------------------------------------------------------+
              | y_source                     | radian    | y-coordinate (Dec-like axis) of the source            |
              +------------------------------+-----------+-------------------------------------------------------+













      ..
          !! processed by numpydoc !!

   .. py:method:: image_properties_epl_shear_lenstronomy(lens_parameters)

      
      Compute image properties for strongly lensed events. This use functions from lenstronomy.

      Solves the lens equation using multiprocessing to find image positions,
      magnifications, time delays, and image types for each lensing event.

      :Parameters:

          **lens_parameters** : ``dict``
              Dictionary containing lens and source parameters shown in the table:

              +------------------------------+-----------+-------------------------------------------------------+
              | Parameter                    | Units     | Description                                           |
              +==============================+===========+=======================================================+
              | zl                           |           | redshift of the lens                                  |
              +------------------------------+-----------+-------------------------------------------------------+
              | zs                           |           | redshift of the source                                |
              +------------------------------+-----------+-------------------------------------------------------+
              | sigma                        | km s^-1   | velocity dispersion                                   |
              +------------------------------+-----------+-------------------------------------------------------+
              | q                            |           | axis ratio                                            |
              +------------------------------+-----------+-------------------------------------------------------+
              | theta_E                      | radian    | Einstein radius                                       |
              +------------------------------+-----------+-------------------------------------------------------+
              | phi                          | rad       | axis rotation angle. counter-clockwise from the       |
              |                              |           | positive x-axis (RA-like axis) to the major axis of   |
              |                              |           | the projected mass distribution.                      |
              +------------------------------+-----------+-------------------------------------------------------+
              | gamma                        |           | density profile slope of EPL galaxy                   |
              +------------------------------+-----------+-------------------------------------------------------+
              | gamma1                       |           | external shear component in the x-direction           |
              |                              |           | (RA-like axis)                                        |
              +------------------------------+-----------+-------------------------------------------------------+
              | gamma2                       |           | external shear component in the y-direction           |
              |                              |           | (Dec-like axis)                                       |
              +------------------------------+-----------+-------------------------------------------------------+

      :Returns:

          **lens_parameters** : ``dict``
              Updated dictionary with additional image properties with the following description:

              +------------------------------+-----------+-------------------------------------------------------+
              | x0_image_positions           | radian    | x-coordinate (RA-like axis) of the images             |
              +------------------------------+-----------+-------------------------------------------------------+
              | x1_image_positions           | radian    | y-coordinate (Dec-like axis) of the images            |
              +------------------------------+-----------+-------------------------------------------------------+
              | magnifications               |           | magnifications                                        |
              +------------------------------+-----------+-------------------------------------------------------+
              | time_delays                  |           | time delays                                           |
              +------------------------------+-----------+-------------------------------------------------------+
              | image_type                   |           | image type                                            |
              +------------------------------+-----------+-------------------------------------------------------+
              | n_images                     |           | number of images                                      |
              +------------------------------+-----------+-------------------------------------------------------+
              | x_source                     | radian    | x-coordinate (RA-like axis) of the source             |
              +------------------------------+-----------+-------------------------------------------------------+
              | y_source                     | radian    | y-coordinate (Dec-like axis) of the source            |
              +------------------------------+-----------+-------------------------------------------------------+













      ..
          !! processed by numpydoc !!

   .. py:method:: get_lensed_snrs(lensed_param, pdet_finder=None, include_effective_parameters=False)

      
      Compute detection probability for each lensed image.

      Calculates the effective luminosity distance, geocent time, and phase
      for each image accounting for magnification and morse phase, then
      computes detection probabilities using the provided calculator.

      :Parameters:

          **lensed_param** : ``dict``
              Dictionary containing lensed source and image parameters given below:

              +------------------------------+-----------+-------------------------------------------------------+
              | Parameter                    | Units     | Description                                           |
              +==============================+===========+=======================================================+
              | geocent_time                 | s         | geocent time                                          |
              +------------------------------+-----------+-------------------------------------------------------+
              | ra                           | rad       | right ascension                                       |
              +------------------------------+-----------+-------------------------------------------------------+
              | dec                          | rad       | declination                                           |
              +------------------------------+-----------+-------------------------------------------------------+
              | phase                        | rad       | phase of GW at reference freq                         |
              +------------------------------+-----------+-------------------------------------------------------+
              | psi                          | rad       | polarization angle                                    |
              +------------------------------+-----------+-------------------------------------------------------+
              | theta_jn                     | rad       | inclination angle                                     |
              +------------------------------+-----------+-------------------------------------------------------+
              | a_1                          |           | spin of the primary compact binary                    |
              +------------------------------+-----------+-------------------------------------------------------+
              | a_2                          |           | spin of the secondary compact binary                  |
              +------------------------------+-----------+-------------------------------------------------------+
              | luminosity_distance          | Mpc       | luminosity distance of the source                     |
              +------------------------------+-----------+-------------------------------------------------------+
              | mass_1                       | Msun      | mass of the primary compact binary (detector frame)   |
              +------------------------------+-----------+-------------------------------------------------------+
              | mass_2                       | Msun      | mass of the secondary compact binary (detector frame) |
              +------------------------------+-----------+-------------------------------------------------------+
              | x0_image_positions           | radian    | x-coordinate (RA-like axis) of the images             |
              +------------------------------+-----------+-------------------------------------------------------+
              | x1_image_positions           | radian    | y-coordinate (Dec-like axis) of the images            |
              +------------------------------+-----------+-------------------------------------------------------+
              | magnifications               |           | magnifications                                        |
              +------------------------------+-----------+-------------------------------------------------------+
              | time_delays                  |           | time delays                                           |
              +------------------------------+-----------+-------------------------------------------------------+
              | image_type                   |           | image type                                            |
              +------------------------------+-----------+-------------------------------------------------------+

          **pdet_finder** : ``callable``
              Function that computes detection probability given GW parameters.

          **include_effective_parameters** : ``bool``
              If True, includes effective parameters in output lensed_param.

      :Returns:

          **result_dict** : ``dict``
              Dictionary containing:

              - 'pdet_net': network detection probability (shape: size x n_max_images)

              - Individual detector probabilities if pdet_finder outputs them

          **lensed_param** : ``dict``
              Updated dictionary with effective parameters shown below:

              +----------------------------------+-----------+------------------------------------------------+
              | Parameter                        | Units     | Description                                    |
              +==================================+===========+================================================+
              | effective_luminosity_distance    | Mpc       | magnification-corrected distance               |
              |                                  |           | luminosity_distance / sqrt(|magnifications_i|) |
              +----------------------------------+-----------+------------------------------------------------|
              | effective_geocent_time           | s         | time-delay-corrected GPS time                  |
              |                                  |           | geocent_time + time_delays_i                   |
              +----------------------------------+-----------+------------------------------------------------|
              | effective_phase                  | rad       | morse-phase-corrected phase                    |
              |                                  |           | phi - morse_phase_i                            |
              +----------------------------------+-----------+------------------------------------------------+
              | effective_ra                     | rad       | RA of the image                                |
              |                                  |           | ra + (x0_image_positions_i - x_source)/cos(dec)|
              +----------------------------------+-----------+------------------------------------------------+
              | effective_dec                    | rad       | Dec of the image                               |
              |                                  |           | dec + (x1_image_positions_i - y_source)        |
              +----------------------------------+-----------+------------------------------------------------+













      ..
          !! processed by numpydoc !!

   .. py:method:: recover_redundant_parameters(lensed_param)

      
      Recover redundant parameters in lensed_param, i.e. theta_E, n_images, mass_1, mass_2, luminosity_distance.
















      ..
          !! processed by numpydoc !!

   .. py:method:: produce_effective_params(lensed_param)

      
      Produce effective parameters for each lensed image.

      Calculates the effective luminosity distance, geocent time, phase,
      RA, and Dec for each image accounting for magnification and morse phase.

      :Parameters:

          **lensed_param** : ``dict``
              Dictionary containing lensed source and image parameters.

      :Returns:

          **lensed_param** : ``dict``
              Updated dictionary with effective parameters shown below:

              +----------------------------------+-----------+------------------------------------------------+
              | Parameter                        | Units     | Description                                    |
              +==================================+===========+================================================+
              | effective_luminosity_distance    | Mpc       | magnification-corrected distance               |
              |                                  |           | luminosity_distance / sqrt(|magnifications_i|) |
              +----------------------------------+-----------+------------------------------------------------|
              | effective_geocent_time           | s         | time-delay-corrected GPS time                  |
              |                                  |           | geocent_time + time_delays_i                   |
              +----------------------------------+-----------+------------------------------------------------|
              | effective_phase                  | rad       | morse-phase-corrected phase                    |
              |                                  |           | phi - morse_phase_i                            |
              +----------------------------------+-----------+------------------------------------------------+
              | effective_ra                     | rad       | RA of the image                                |
              |                                  |           | ra + (x0_image_positions_i - x_source)/cos(dec)|
              +----------------------------------+-----------+------------------------------------------------+
              | effective_dec                    | rad       | Dec of the image                               |
              |                                  |           | dec + (x1_image_positions_i - y_source)        |
              +----------------------------------+-----------+------------------------------------------------+













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

           - [0]: e1 - ellipticity component 1

           - [1]: e2 - ellipticity component 2

           - [2]: gamma - power-law slope of mass density

           - [3]: gamma1 - external shear component 1

           - [4]: gamma2 - external shear component 2

           - [5]: zl - lens redshift

           - [6]: zs - source redshift

           - [7]: einstein_radius - Einstein radius (units: radians)

           - [8]: iteration - iteration index for tracking

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

   >>> from ler.image_properties.multiprocessing_routine import solve_lens_equation, _init_worker_multiprocessing
   >>> import numpy as np
   >>> from multiprocessing import Pool
   >>> lens_parameters1 = np.array([0.024, -0.016, 1.89, 0.10, 0.09, 0.25, 0.94, 2.5e-06, 0])
   >>> lens_parameters2 = np.array([-0.040, -0.014, 2.00, 0.08, -0.01, 1.09, 2.55, 1.0e-06, 1])
   >>> input_arguments = np.vstack((lens_parameters1, lens_parameters2))
   >>> with Pool(
   ...     processes=2, # Number of worker processes
   ...     initializer=_init_worker_multiprocessing, # common
   ...     initargs=(
   ...         2, # n_min_images
   ...         ['EPL_NUMBA', 'SHEAR'], # lensModelList
   ...     ),
   ... ) as pool:
   ...     result = pool.map(solve_lens_equation, input_arguments)



   ..
       !! processed by numpydoc !!

.. py:function:: sample_source_from_double_caustic(theta_E, q, phi, gamma, gamma1, gamma2, num_th=500, maginf=-100.0)

   
   Sample a single source position from the double caustic region.

   Precomputes the caustic boundary once (the expensive part) and
   draws one uniform sample via exact fan-triangulation from the origin.

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
           Number of angular samples for caustic boundary.

           default: 500

       **maginf** : ``float``
           Magnification cut threshold.

           default: -100.0

   :Returns:

       **xs** : ``float``
           Sampled source x-coordinate (``NaN`` if caustic is invalid).

       **ys** : ``float``
           Sampled source y-coordinate (``NaN`` if caustic is invalid).










   .. rubric:: Examples

   >>> xs, ys = sample_source_from_double_caustic(
   ...     theta_E=1.0, q=0.8, phi=0.0, gamma=2.0, gamma1=0.03, gamma2=-0.01
   ... )



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

.. py:function:: sample_source_from_double_caustic(theta_E, q, phi, gamma, gamma1, gamma2, num_th=500, maginf=-100.0)

   
   Sample a single source position from the double caustic region.

   Precomputes the caustic boundary once (the expensive part) and
   draws one uniform sample via exact fan-triangulation from the origin.

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
           Number of angular samples for caustic boundary.

           default: 500

       **maginf** : ``float``
           Magnification cut threshold.

           default: -100.0

   :Returns:

       **xs** : ``float``
           Sampled source x-coordinate (``NaN`` if caustic is invalid).

       **ys** : ``float``
           Sampled source y-coordinate (``NaN`` if caustic is invalid).










   .. rubric:: Examples

   >>> xs, ys = sample_source_from_double_caustic(
   ...     theta_E=1.0, q=0.8, phi=0.0, gamma=2.0, gamma1=0.03, gamma2=-0.01
   ... )



   ..
       !! processed by numpydoc !!

