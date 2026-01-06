:py:mod:`ler.lens_galaxy_population.lens_galaxy_parameter_distribution`
=======================================================================

.. py:module:: ler.lens_galaxy_population.lens_galaxy_parameter_distribution

.. autoapi-nested-parse::

   This module contains the ``LensGalaxyParameterDistribution`` class.

   The class is used to sample lens galaxy parameters and source parameters conditioned on
   the source being strongly lensed.
   It inherits from the ``ImageProperties`` class to calculate image properties (magnification,
   time delays, source position, image position, morse phase).
   It also inherits from ``CBCSourceParameterDistribution``.

   Copyright (c) 2026 Phurailatpam Hemantakumar
   License: MIT

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

   
   Class to sample lens galaxy parameters and source parameters conditioned on the source being strongly lensed.

   This class deals with the distribution of lens galaxy parameters, such as velocity dispersion,
   axis ratio, axis rotation angle, shear, and density profile slope. It also handles the
   sampling of source parameters conditioned on the source being strongly lensed.

   :Parameters:

       **npool** : int, optional
           Number of processors to use.
           Default is 4.

       **z_min** : float, optional
           Minimum redshift.
           Default is 0.0.

       **z_max** : float, optional
           Maximum redshift.
           Default is 10.0.

       **cosmology** : astropy.cosmology, optional
           Cosmology to use.
           Default is None, which falls back to ``astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)``.

       **event_type** : str, optional
           Type of event to generate. e.g. 'BBH', 'BNS', 'NSBH'.
           Default is 'BBH'.

       **lens_type** : str, optional
           Type of lens galaxy to generate.
           Default is 'epl_shear_galaxy'.

       **lens_functions** : dict, optional
           Dictionary of lens functions.

       **lens_functions_params** : dict, optional
           Dictionary of parameters for lens functions.

       **lens_param_samplers** : dict, optional
           Dictionary of lens parameter samplers.

       **lens_param_samplers_params** : dict, optional
           Dictionary of parameters for lens parameter samplers.

       **directory** : str, optional
           Directory to store the interpolators.
           Default is './interpolator_json'.

       **create_new_interpolator** : bool, optional
           If True, creates a new interpolator.
           Default is False.

       **buffer_size** : int, optional
           Buffer size for sampling lens parameters.
           Default is 1000.

       **\*\*kwargs**
           Keyword arguments to pass to the parent classes.











   .. rubric:: Examples

   >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
   >>> lens = LensGalaxyParameterDistribution()
   >>> lensed_params = lens.sample_lens_parameters(size=1000)
   >>> print(lensed_params.keys())

   :Attributes:

       **npool** : int
           Number of processors to use.

       **z_min** : float
           Minimum redshift.

       **z_max** : float
           Maximum redshift.

       **cosmo** : astropy.cosmology
           Cosmology object.

       **event_type** : str
           Type of event to generate.

       **directory** : str
           Directory to store the interpolators.

       **create_new_interpolator** : dict
           Dictionary to check if new interpolator is created.

       **lens_param_samplers** : dict
           Dictionary of lens parameter samplers.

       **lens_param_samplers_params** : dict
           Dictionary of lens parameter sampler parameters.

       **lens_functions** : dict
           Dictionary of lens functions.

       **normalization_pdf_z_lensed** : float
           Normalization constant of the pdf p(z) for lensed events.


   ..
       !! processed by numpydoc !!
   .. py:attribute:: cbc_pop
      :value: 'None'

      
      Inherited class for sampling source parameters.
















      ..
          !! processed by numpydoc !!

      :type: :class:`~ler.gw_source_population.CBCSourceParameterDistribution`

   .. py:attribute:: z_min
      :value: 'None'

      
      Minimum redshift.
















      ..
          !! processed by numpydoc !!

      :type: float

   .. py:attribute:: z_max
      :value: 'None'

      
      Maximum redshift.
















      ..
          !! processed by numpydoc !!

      :type: float

   .. py:attribute:: m_min
      :value: 'None'

      
      Minimum mass in detector frame.
















      ..
          !! processed by numpydoc !!

      :type: float

   .. py:attribute:: m_max
      :value: 'None'

      
      Maximum mass in detector frame.
















      ..
          !! processed by numpydoc !!

      :type: float

   .. py:attribute:: normalization_pdf_z
      :value: 'None'

      
      Normalization constant of the pdf p(z).
















      ..
          !! processed by numpydoc !!

      :type: float

   .. py:attribute:: event_type
      :value: "'BBH'"

      
      ``str``

      Type of event to generate.

      e.g. 'BBH', 'BNS', 'NSBH'















      ..
          !! processed by numpydoc !!

   .. py:attribute:: sample_source_redshift_sl

      

   .. py:attribute:: sample_lens_parameters_routine

      

   .. py:attribute:: cross_section_based_sampler

      

   .. py:attribute:: normalization_pdf_z_lensed

      
      Normalization constant of the pdf p(z) for lensed events.
















      ..
          !! processed by numpydoc !!

      :type: float

   .. py:method:: class_initialization_lens(npool, z_min, z_max, cosmology, lens_type, lens_functions, lens_functions_params, lens_param_samplers, lens_param_samplers_params, directory, create_new_interpolator, params)

      
      Initialize the LensGalaxyParameterDistribution class.


      :Parameters:

          **npool** : int
              Number of processors to use for sampling.

          **z_min** : float
              Minimum redshift of the lens galaxy.

          **z_max** : float
              Maximum redshift of the lens galaxy.

          **cosmology** : astropy.cosmology
              Cosmology object.

          **lens_type** : str
              Type of the lens galaxy.

          **lens_functions** : dict
              Dictionary with the lens related functions.

          **lens_functions_params** : dict
              Dictionary with the parameters for the lens related functions.

          **lens_param_samplers** : dict
              Dictionary with the priors for the sampler.

          **lens_param_samplers_params** : dict
              Dictionary with the parameters for the priors of the sampler.

          **directory** : str
              Directory where the interpolators are saved.

          **create_new_interpolator** : bool
              If True, creates a new interpolator.

          **params** : dict
              Additional parameters for the ``CBCSourceParameterDistribution`` and ``ImageProperties`` classes.














      ..
          !! processed by numpydoc !!

   .. py:method:: sample_lens_parameters(size=1000)

      
      Sample lens galaxy parameters along with the source parameters, conditioned on the source being strongly lensed.


      :Parameters:

          **size** : int, optional
              Number of lens parameters to sample.
              Default is 1000.

      :Returns:

          **lens_parameters** : dict
              Dictionary of sampled lens parameters and source parameters.
              Keys include ``zl``, ``zs``, ``sigma``, ``q``, ``theta_E``, ``phi``, ``e1``, ``e2``,
              ``gamma1``, ``gamma2``, ``gamma``, ``geocent_time``, ``phase``, ``psi``, ``theta_jn``,
              ``luminosity_distance``, ``mass_1_source``, ``mass_2_source``, ``ra``, ``dec``.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> od = LensGalaxyParameterDistribution(lens_param_samplers=dict(velocity_dispersion="velocity_dispersion_ewoud"))
      >>> print(od.sample_lens_parameters(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: sample_all_routine_epl_shear_sl(size=1000)

      
      Sample galaxy lens parameters. EPL shear cross section is used for rejection sampling.


      :Parameters:

          **size** : int, optional
              Number of lens parameters to sample.
              Default is 1000.

      :Returns:

          **lens_parameters** : dict
              Dictionary of lens parameters and source parameters (lens conditions applied).
              Keys include ``zl``, ``zs``, ``sigma``, ``q``, ``theta_E``, ``phi``, ``e1``, ``e2``,
              ``gamma1``, ``gamma2``, ``gamma``.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.sample_all_routine_epl_shear_sl(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: strongly_lensed_source_redshifts(size=1000)

      
      Sample source redshifts, conditioned on the source being strongly lensed.


      :Parameters:

          **size** : int, optional
              Number of lens parameters to sample.
              Default is 1000.

      :Returns:

          **redshifts** : numpy.ndarray
              Source redshifts conditioned on the source being strongly lensed.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.strongly_lensed_source_redshifts(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: sample_all_routine_epl_shear_intrinsic(size=1000)

      
      Sample galaxy lens parameters. EPL shear cross section is used for rejection sampling.


      :Parameters:

          **size** : int, optional
              Number of lens parameters to sample.
              Default is 1000.

      :Returns:

          **lens_parameters** : dict
              Dictionary of lens parameters and source parameters (lens conditions applied).
              Keys include ``zl``, ``zs``, ``sigma``, ``q``, ``theta_E``, ``phi``, ``e1``, ``e2``,
              ``gamma1``, ``gamma2``, ``gamma``.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.sample_all_routine_epl_shear_intrinsic(size=1000)



      ..
          !! processed by numpydoc !!


