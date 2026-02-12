:py:mod:`ler.lens_galaxy_population.lens_galaxy_parameter_distribution`
=======================================================================

.. py:module:: ler.lens_galaxy_population.lens_galaxy_parameter_distribution

.. autoapi-nested-parse::

   Module for lens galaxy parameter distributions.

   This module contains the ``LensGalaxyParameterDistribution`` class which samples
   lens galaxy parameters and source parameters conditioned on the source being
   strongly lensed. It handles the distribution of lens parameters such as velocity
   dispersion, axis ratio, rotation angle, shear, and density profile slope.

   Inheritance hierarchy:

   - :class:`~ler.gw_source_population.CBCSourceParameterDistribution`

   - :class:`~ler.image_properties.ImageProperties`

   - :class:`~ler.lens_galaxy_population.OpticalDepth`


   Usage:
       Basic workflow example:

       >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
       >>> lens = LensGalaxyParameterDistribution()
       >>> lensed_params = lens.sample_lens_parameters(size=1000)
       >>> print(lensed_params.keys())

   Copyright (C) 2026 Phurailatpam Hemantakumar. Distributed under MIT License.

   ..
       !! processed by numpydoc !!


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.lens_galaxy_parameter_distribution.LensGalaxyParameterDistribution




.. py:class:: LensGalaxyParameterDistribution(npool=4, z_min=0.0, z_max=10.0, cosmology=None, event_type='BBH', lens_type='epl_shear_galaxy', lens_functions=None, lens_functions_params=None, lens_param_samplers=None, lens_param_samplers_params=None, directory='./interpolator_json', create_new_interpolator=False, buffer_size=1000, **kwargs)


   Bases: :py:obj:`ler.gw_source_population.CBCSourceParameterDistribution`, :py:obj:`ler.image_properties.ImageProperties`, :py:obj:`ler.lens_galaxy_population.optical_depth.OpticalDepth`

   
   Sample lens galaxy parameters conditioned on strong lensing.

   This class handles the distribution of lens galaxy parameters such as velocity
   dispersion, axis ratio, axis rotation angle, shear, and density profile slope.
   It samples source parameters conditioned on the source being strongly lensed
   using cross-section based rejection or importance sampling.

   Key Features:

   - Samples lens parameters using EPL+shear galaxy model

   - Supports rejection and importance sampling based on cross-section

   - Computes optical depth weighted source redshift distributions

   - Integrates with GW source population and image property calculations

   :Parameters:

       **npool** : ``int``
           Number of processors to use for parallel sampling.

           default: 4

       **z_min** : ``float``
           Minimum redshift for source and lens populations.

           default: 0.0

       **z_max** : ``float``
           Maximum redshift for source and lens populations.

           default: 10.0

       **cosmology** : ``astropy.cosmology`` or ``None``
           Cosmology object for distance calculations.

           default: None (uses FlatLambdaCDM with H0=70, Om0=0.3)

       **event_type** : ``str``
           Type of compact binary coalescence event.

           Options:

           - 'BBH': Binary black hole

           - 'BNS': Binary neutron star

           - 'NSBH': Neutron star-black hole

           default: 'BBH'

       **lens_type** : ``str``
           Type of lens galaxy model to use.

           default: 'epl_shear_galaxy'

       **lens_functions** : ``dict`` or ``None``
           Dictionary specifying lens-related functions.

           default: None (uses defaults from OpticalDepth)

       **lens_functions_params** : ``dict`` or ``None``
           Parameters for lens functions.

           default: None

       **lens_param_samplers** : ``dict`` or ``None``
           Dictionary specifying lens parameter sampling functions.

           default: None (uses defaults from OpticalDepth)

       **lens_param_samplers_params** : ``dict`` or ``None``
           Parameters for lens parameter samplers.

           default: None

       **directory** : ``str``
           Directory for storing interpolator files.

           default: './interpolator_json'

       **create_new_interpolator** : ``bool``
           If True, recreates interpolators even if files exist.

           default: False

       **buffer_size** : ``int``
           Buffer size for batch sampling of lens parameters.

           default: 1000

       **\*\*kwargs** : ``dict``
           Additional keyword arguments passed to parent classes:

           :class:`~ler.gw_source_population.CBCSourceParameterDistribution`,

           :class:`~ler.image_properties.ImageProperties`,

           :class:`~ler.lens_galaxy_population.OpticalDepth`.











   .. rubric:: Examples

   Basic usage:

   >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
   >>> lens = LensGalaxyParameterDistribution()
   >>> lensed_params = lens.sample_lens_parameters(size=1000)
   >>> print(lensed_params.keys())

   Instance Methods
   ----------
   LensGalaxyParameterDistribution has the following methods:

   +-----------------------------------------------------+------------------------------------------------+
   | Method                                              | Description                                    |
   +=====================================================+================================================+
   | :meth:`~sample_lens_parameters`                     | Sample lens and source parameters              |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~sample_all_routine_epl_shear_intrinsic`    | Sample EPL+shear lens parameters from intrinsic |
   |                                                     | distributions                                  |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~sample_all_routine_epl_shear_sl`            | Sample EPL+shear lens parameters with strong   |
   |                                                     | lensing condition                              |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~strongly_lensed_source_redshift`           | Sample source redshifts with lensing condition |
   +-----------------------------------------------------+------------------------------------------------+

   Instance Attributes
   ----------
   LensGalaxyParameterDistribution has the following attributes:

   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | Attribute                                      | Type                 | Unit  | Description                                    |
   +================================================+======================+=======+================================================+
   | :attr:`~npool`                                 | ``int``              |       | Number of processors for parallel computation  |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | :attr:`~z_min`                                 | ``float``            |       | Minimum redshift                               |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | :attr:`~z_max`                                 | ``float``            |       | Maximum redshift                               |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | :attr:`~cosmo`                                 | ``astropy.cosmology``|       | Cosmology object for calculations              |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | :attr:`~event_type`                            | ``str``              |       | Type of CBC event (BBH, BNS, NSBH)             |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | :attr:`~directory`                             | ``str``              |       | Path to interpolator storage directory         |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | :attr:`~lens_param_samplers`                   | ``dict``             |       | Dictionary of lens parameter sampler names     |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | :attr:`~lens_param_samplers_params`            | ``dict``             |       | Parameters for lens parameter samplers         |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | :attr:`~lens_functions`                        | ``dict``             |       | Dictionary of lens function names              |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | :attr:`~normalization_pdf_z_lensed`            | ``float``            |       | Normalization constant for lensed source z pdf |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:property:: source_redshift_sl

      
      Function to sample source redshifts conditioned on strong lensing.



      :Returns:

          **source_redshift_sl** : ``ler.functions.FunctionConditioning``
              Function for sampling source redshifts conditioned on strong lensing.













      ..
          !! processed by numpydoc !!

   .. py:property:: normalization_pdf_z_lensed

      
      Normalization constant for the lensed source redshift pdf.

      This constant is used to normalize the probability distribution

      of source redshifts conditioned on strong lensing. It is computed

      by integrating the merger rate density times optical depth.


      :Returns:

          **normalization_pdf_z_lensed** : ``float``
              Normalization constant for lensed redshift distribution.













      ..
          !! processed by numpydoc !!

   .. py:property:: lens_param_samplers

      
      Dictionary of lens parameter sampler function names.



      :Returns:

          **lens_param_samplers** : ``dict``
              Dictionary mapping parameter names to sampler function names.

              Keys include: 'source_redshift_sl', 'lens_redshift',

              'velocity_dispersion', 'axis_ratio', 'axis_rotation_angle',

              'external_shear', 'density_profile_slope'.













      ..
          !! processed by numpydoc !!

   .. py:property:: lens_param_samplers_params

      
      Dictionary of parameters for lens parameter samplers.



      :Returns:

          **lens_param_samplers_params** : ``dict``
              Dictionary with sampler parameters.

              Each key corresponds to a sampler in lens_param_samplers.













      ..
          !! processed by numpydoc !!

   .. py:property:: lens_functions

      
      Dictionary of lens-related function names.



      :Returns:

          **lens_functions** : ``dict``
              Dictionary mapping function types to function names.

              Keys include: 'param_sampler_type', 'cross_section_based_sampler',

              'optical_depth', 'cross_section'.













      ..
          !! processed by numpydoc !!

   .. py:attribute:: event_type
      :value: "'BBH'"

      
      ``str``

      Type of event to generate.

      e.g. 'BBH', 'BNS', 'NSBH'















      ..
          !! processed by numpydoc !!

   .. py:attribute:: sample_lens_parameters_routine

      

   .. py:attribute:: cross_section_based_sampler

      

   .. py:method:: sample_lens_parameters(size=1000)

      
      Sample lens galaxy and source parameters conditioned on strong lensing.

      This method samples both lens galaxy parameters (velocity dispersion, axis
      ratio, shear, etc.) and gravitational wave source parameters, with the
      source redshift distribution weighted by strong lensing optical depth.

      :Parameters:

          **size** : ``int``
              Number of lens-source parameter sets to sample.

              default: 1000

      :Returns:

          **lens_parameters** : ``dict``
              Dictionary containing sampled lens and source parameters.

              The included parameters and their units are as follows (for default settings):

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
              | mass_1_source                | Msun      | mass of the primary compact binary (source frame)     |
              +------------------------------+-----------+-------------------------------------------------------+
              | mass_2_source                | Msun      | mass of the secondary compact binary (source frame)   |
              +------------------------------+-----------+-------------------------------------------------------+
              | mass_1                       | Msun      | mass of the primary compact binary (detector frame)   |
              +------------------------------+-----------+-------------------------------------------------------+
              | mass_2                       | Msun      | mass of the secondary compact binary (detector frame) |
              +------------------------------+-----------+-------------------------------------------------------+










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> params = lens.sample_lens_parameters(size=1000)
      >>> print(params.keys())



      ..
          !! processed by numpydoc !!

   .. py:method:: sample_all_routine_epl_shear_sl(size=1000)

      
      Sample EPL+shear galaxy lens parameters with strong lensing condition.


      :Parameters:

          **size** : ``int``
              Number of lens parameters to sample.

              default: 1000

      :Returns:

          **lens_parameters** : ``dict``
              Dictionary of sampled lens parameters.

              The included parameters and their units are as follows (for default settings):

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













      ..
          !! processed by numpydoc !!

   .. py:method:: strongly_lensed_source_redshift(size, get_attribute=False, **kwargs)

      
      Sample source redshifts conditioned on strong lensing.


      :Parameters:

          **size** : ``int``
              Number of samples to generate.

              default: 1000

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters

      :Returns:

          **redshifts** : ``numpy.ndarray``
              Array of source redshifts conditioned on strong lensing.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> zs = lens.strongly_lensed_source_redshift(size=1000)
      >>> print(f"strongly lensed source redshift: {zs.mean():.2f}")



      ..
          !! processed by numpydoc !!

   .. py:method:: sample_all_routine_epl_shear_intrinsic(size=1000)

      
      Sample EPL+shear galaxy lens parameters from intrinsic distributions.

      Samples lens parameters from their intrinsic distributions without
      applying strong lensing cross-section weighting.

      :Parameters:

          **size** : ``int``
              Number of lens parameters to sample.

              default: 1000

      :Returns:

          **lens_parameters** : ``dict``
              Dictionary of sampled lens parameters.

              The included parameters and their units are as follows (for default settings):

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













      ..
          !! processed by numpydoc !!


