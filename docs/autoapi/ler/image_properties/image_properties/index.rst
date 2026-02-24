:py:mod:`ler.image_properties.image_properties`
===============================================

.. py:module:: ler.image_properties.image_properties

.. autoapi-nested-parse::

   Module for computing image properties of strongly lensed gravitational wave events.

   This module contains the ImageProperties class, which computes image positions,
   magnifications, time delays, and other lensing-related quantities for gravitational
   wave sources that are strongly lensed by intervening galaxies. The class uses
   multiprocessing to efficiently solve lens equations for large samples.

   Usage:
       Basic workflow example:

       >>> from ler.image_properties import ImageProperties
       >>> ip = ImageProperties()
       >>> lens_parameters = dict(zs=np.array([2.0]), zl=np.array([0.5]), ...)
       >>> result = ip.image_properties(lens_parameters)

   Copyright (C) 2026 Phurailatpam Hemanta Kumar. Distributed under MIT License.

   ..
       !! processed by numpydoc !!


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   ler.image_properties.image_properties.ImageProperties




Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.image_properties.image_properties.cosmo


.. py:data:: cosmo

   

.. py:class:: ImageProperties(npool=4, n_min_images=2, n_max_images=4, lens_model_list=['EPL_NUMBA', 'SHEAR'], cosmology=None, time_window=365 * 24 * 3600 * 2, spin_zero=True, spin_precession=False, pdet_finder=None, include_effective_parameters=False, multiprocessing_verbose=True, include_redundant_parameters=False)


   
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
   | :meth:`~image_properties`                           | Compute image properties for lensed events     |
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
   | :attr:`~include_redundant_parameters`                | ``bool``                  |          | If True, removes redundant parameters from output to save memory |
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

   .. py:attribute:: multiprocessing_verbose
      :value: 'True'

      

   .. py:attribute:: include_redundant_parameters
      :value: 'False'

      

   .. py:method:: image_properties(lens_parameters)

      
      Compute image properties for strongly lensed events.

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
              | luminosity_distance          | Mpc       | luminosity distance of the source                      |
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

              +----------------------------------+-----------+------------------------------------------------|
              | Parameter                        | Units     | Description
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

              +----------------------------------+-----------+------------------------------------------------|
              | Parameter                        | Units     | Description
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


