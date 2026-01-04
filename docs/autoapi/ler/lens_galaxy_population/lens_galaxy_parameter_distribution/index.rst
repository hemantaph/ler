:py:mod:`ler.lens_galaxy_population.lens_galaxy_parameter_distribution`
=======================================================================

.. py:module:: ler.lens_galaxy_population.lens_galaxy_parameter_distribution

.. autoapi-nested-parse::

   This module contains the LensGalaxyPopulation class, which is used to sample lens galaxy parameters, source parameters conditioned on the source being strongly lensed.

   The class inherits from the ImageProperties class, which is used calculate image properties (magnification, timedelays, source position, image position, morse phase).

   Either the class takes in initialized CBCSourceParameterDistribution class as input or inherits the CBCSourceParameterDistribution class with default params (if no input)

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

   
   Class to sample lens galaxy parameters, source parameters conditioned on the source being strongly lensed, and image properties


   :Parameters:

       **npool** : `int`
           number of processors to use

       **z_min** : `float`
           minimum redshift

       **z_max** : `float`
           maximum redshift

       **cosmology** : `astropy.cosmology`
           Cosmology to use
           default: None/astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)

       **event_type** : `str`
           Type of event to generate.
           e.g. 'BBH', 'BNS', 'NSBH'
           default: 'BBH'

       **lens_type** : `str`
           Type of lens galaxy to generate.
           default: 'epl_shear_galaxy'

       **lens_functions, lens_param_samplers, lens_param_samplers_params** : `dict`, `dict`, `dict`
           dictionary of lens functions, priors, and priors parameters
           Check for default/available lens functions, priors and corresponding input parameters by running,

           >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
           >>> lens = LensGalaxyParameterDistribution()
           >>> print(lens.lens_functions)
           >>> print(lens.lens_param_samplers)
           >>> print(lens.lens_param_samplers_params)

       **directory** : `str`
           directory to store the interpolators
           default: './interpolator_json'

       **\*\*kwargs**
           keyword arguments to pass to the parent classes











   .. rubric:: Examples

   >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
   >>> lens = LensGalaxyParameterDistribution()
   >>> lensed_params = lens.sample_lens_parameters(size=1000)
   >>> lensed_params.keys()

   Instance Attributes
   ----------
   LensGalaxyPopulation class has the following instance attributes:

   +-------------------------------------+----------------------------------+
   | Atrributes                          | Type                             |
   +=====================================+==================================+
   |:attr:`~npool`                       | `int`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~z_min`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~z_max`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~cosmo`                       | `astropy.cosmology`              |
   +-------------------------------------+----------------------------------+
   |:attr:`~event_type`                  | `str`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~directory`                   | `str`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~create_new_interpolator`     | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~lens_param_samplers`         | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~lens_param_samplers_params`  | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~lens_sampler_names`          | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~lens_functions`              | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~normalization_pdf_z_lensed`  | `float`                          |
   +-------------------------------------+----------------------------------+

   Instance Methods
   ----------
   LensGalaxyPopulation class has the following instance methods:

   +-------------------------------------+----------------------------------+
   | Methods                             | Type                             |
   +=====================================+==================================+
   |:meth:`~sample_lens_parameters`      | Function to call the specific    |
   |                                     | galaxy lens parameters sampler   |
   |                                     | routine.                         |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_all_routine_sie`      | Function to sample galaxy lens   |
   |                                     | parameters along with the source |
   |                                     | parameters.                      |
   +-------------------------------------+----------------------------------+
   |:meth:`~strongly_lensed_source_redshifts`                               |
   +-------------------------------------+----------------------------------+
   |                                     | Function to sample source        |
   |                                     | redshifts conditioned on the     |
   |                                     | source being strongly lensed     |
   +-------------------------------------+----------------------------------+
   |:meth:`~source_parameters`           | Function to sample gw source     |
   |                                     | parameters                       |
   +-------------------------------------+----------------------------------+
   |:meth:`~lens_redshift_SDSS_catalogue`| Function to sample lens          |
   |                                     | redshifts, conditioned on the    |
   |                                     | lens being strongly lensed       |
   +-------------------------------------+----------------------------------+
   |:meth:`~axis_rotation_angle_uniform` | Function to sample the axis      |
   |                                     | rotation angle of the elliptical |
   |                                     | lens galaxy from a uniform       |
   |                                     | distribution                     |
   +-------------------------------------+----------------------------------+
   |:meth:`~shear_norm`                  | Function to sample the           |
   |                                     | elliptical lens galaxy shear     |
   |                                     | from a normal distribution       |
   +-------------------------------------+----------------------------------+
   |:meth:`~density_profile_slope_normal`                             |
   +-------------------------------------+----------------------------------+
   |                                     | Function to sample the lens      |
   |                                     | galaxy spectral index of the     |
   |                                     | mass density profile from a      |
   |                                     | normal distribution              |
   +-------------------------------------+----------------------------------+
   |:meth:`~compute_einstein_radii`      | Function to compute the Einstein |
   |                                     | radii of the lens galaxies       |
   +-------------------------------------+----------------------------------+
   |:meth:`~rejection_sampling_with_cross_section_sis`  | Function to conduct rejection    |
   |                                     | sampling wrt einstein radius     |
   +-------------------------------------+----------------------------------+
   |:meth:`~rejection_sampling_with_cross_section_sie`  | Function to conduct rejection    |
   |                                     | sampling wrt cross_section       |
   +-------------------------------------+----------------------------------+
   |:attr:`~rejection_sample_sl`         | Function to conduct rejection    |
   |                                     | sampling with the given rejection|
   |                                     | sampling function                |
   +-------------------------------------+----------------------------------+
   |:attr:`~sample_source_redshift_sl`   | Function to sample source        |
   |                                     | redshifts conditioned on the     |
   |                                     | source being strongly lensed     |
   +-------------------------------------+----------------------------------+
   |:attr:`~sample_lens_redshift`        | Function to sample lens          |
   |                                     | redshifts, conditioned on the    |
   |                                     | lens being strongly lensed       |
   +-------------------------------------+----------------------------------+
   |:attr:`~sample_axis_rotation_angle`  | Function to sample the axis      |
   |                                     | rotation angle of the elliptical |
   |                                     | lens galaxy from a uniform       |
   |                                     | distribution                     |
   +-------------------------------------+----------------------------------+
   |:attr:`~sample_shear`                | Function to sample the           |
   |                                     | elliptical lens galaxy shear     |
   |                                     | from a normal distribution       |
   +-------------------------------------+----------------------------------+
   |:attr:`~sample_density_profile_slope`                             |
   +-------------------------------------+----------------------------------+
   |                                     | Function to sample the lens      |
   |                                     | galaxy spectral index of the     |
   |                                     | mass density profile from a      |
   |                                     | normal distribution              |
   +-------------------------------------+----------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:attribute:: cbc_pop

      
      :class:`~CBCSourceParameterDistribution` class

      This is an already initialized class that contains a function (CBCSourceParameterDistribution.sample_gw_parameters) that actually samples the source parameters.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: z_min

      
      `float`

      minimum redshift















      ..
          !! processed by numpydoc !!

   .. py:attribute:: z_max

      
      `float`

      maximum redshift















      ..
          !! processed by numpydoc !!

   .. py:attribute:: m_min

      
      `float`

      minimum mass in detector frame















      ..
          !! processed by numpydoc !!

   .. py:attribute:: m_max

      
      `float`

      maximum mass in detector frame















      ..
          !! processed by numpydoc !!

   .. py:attribute:: normalization_pdf_z

      
      `float`

      normalization constant of the pdf p(z)















      ..
          !! processed by numpydoc !!

   .. py:method:: class_initialization_lens(npool, z_min, z_max, cosmology, lens_type, lens_functions, lens_functions_params, lens_param_samplers, lens_param_samplers_params, directory, create_new_interpolator, params)

      
      Initialize the LensGalaxyParameterDistribution class.


      :Parameters:

          **npool** : `int`
              number of processors to use for sampling

          **z_min** : `float`
              minimum redshift of the lens galaxy

          **z_max** : `float`
              maximum redshift of the lens galaxy

          **cosmology** : `astropy.cosmology`
              cosmology object

          **lens_type** : `str`
              type of the lens galaxy

          **lens_functions** : `dict`
              dictionary with the lens related functions

          **lens_functions_params** : `dict`
              dictionary with the parameters for the lens related functions

          **lens_param_samplers** : `dict`
              dictionary with the priors for the sampler

          **lens_param_samplers_params** : `dict`
              dictionary with the parameters for the priors of the sampler

          **directory** : `str`
              directory where the interpolators are saved

          **create_new_interpolator** : `bool`
              if True, creates a new interpolator

          **params** : `dict`
              additional parameters for the CBCSourceParameterDistribution and ImageProperties classes














      ..
          !! processed by numpydoc !!

   .. py:method:: sample_lens_parameters(size=1000)

      
      Function to sample galaxy lens parameters along with the source parameters, conditioned on the source being strongly lensed.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **lens_parameters** : `dict`
              dictionary of sampled lens parameters and source parameters.

              keys:

              zl: lens redshifts

              zs: source redshifts, lensed condition applied

              sigma: velocity dispersions

              q: axis ratios

              theta_E: Einstein radii

              phi: axis rotation angle

              e1: ellipticity component 1

              e2: ellipticity component 2

              gamma1: shear component 1

              gamma2: shear component 2

              gamma: density profile slope distribution

              geocent_time: time of arrival of the unlensed signal

              phase: phase of the unlensed signal

              psi: polarization angle of the unlensed signal

              theta_jn: inclination angle of the unlensed signal

              luminosity_distance: luminosity distance of the source

              mass_1_source: mass 1 (larger) of the source

              mass_2_source: mass 2 (smaller) of the source

              ra: right ascension of the source

              dec: declination of the source










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> od = LensGalaxyParameterDistribution(lens_param_samplers=dict(velocity_dispersion="velocity_dispersion_ewoud"))
      >>> print(od.sample_lens_parameters(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: sample_all_routine_sie_sl(size=1000)

      
      Function to sample galaxy lens parameters. SIE cross section is used for rejection sampling.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **lens_parameters_input** : `dict`
              dictionary of lens parameters to sample

      :Returns:

          **lens_parameters** : `dict`
              dictionary of lens parameters and source parameters (lens conditions applied):

              zl: lens redshifts

              zs: source redshifts, lensed condition applied

              sigma: velocity dispersions

              q: axis ratios

              theta_E: Einstein radii

              phi: axis rotation angle

              e1: ellipticity component 1

              e2: ellipticity component 2

              gamma1: shear component 1

              gamma2: shear component 2

              gamma: density profile slope distribution










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.sample_all_routine_sie(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: sample_all_routine_epl_shear_sl(size=1000)

      
      Function to sample galaxy lens parameters along. EPL shear cross section is used for rejection sampling.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **lens_parameters_input** : `dict`
              dictionary of lens parameters to sample

      :Returns:

          **lens_parameters** : `dict`
              dictionary of lens parameters and source parameters (lens conditions applied):

              zl: lens redshifts

              zs: source redshifts, lensed condition applied

              sigma: velocity dispersions

              q: axis ratios

              theta_E: Einstein radii

              phi: axis rotation angle

              e1: ellipticity component 1

              e2: ellipticity component 2

              gamma1: shear component 1

              gamma2: shear component 2

              gamma: density profile slope distribution










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.sample_all_routine_sie(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: sampling_routine_sis_nsl(zl, zs, size=1000)

      
      Function to sample SIS lens related parameters.


      :Parameters:

          **zl** : `float`
              lens redshifts

          **zs** : `float`
              source redshifts

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **lens_parameters** : `dict`
              dictionary of sampled lens parameters.
              keys: sigma, theta_E













      ..
          !! processed by numpydoc !!

   .. py:method:: sampling_routine_sie_nsl(zl, zs, size=1000)

      
      Function to sample SIE lens related parameters.


      :Parameters:

          **zl** : `float`
              lens redshifts

          **zs** : `float`
              source redshifts

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **lens_parameters** : `dict`
              dictionary of sampled lens parameters.
              keys: sigma, q, phi













      ..
          !! processed by numpydoc !!

   .. py:method:: sampling_routine_epl_shear_nsl(zl, zs, size=1000)

      
      Function to sample EPL and shear related parameters.


      :Parameters:

          **zl** : `float`
              lens redshifts

          **zs** : `float`
              source redshifts

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **lens_parameters** : `dict`
              dictionary of sampled lens parameters.
              keys: sigma, q, phi, gamma, gamma1, gamma2













      ..
          !! processed by numpydoc !!

   .. py:method:: strongly_lensed_source_redshifts(size=1000)

      
      Function to sample source redshifts, conditioned on the source being strongly lensed.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **redshifts** : `float`
              source redshifts conditioned on the source being strongly lensed










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.strongly_lensed_source_redshifts(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: rejection_sampling_with_cross_section_sis(param_dict, cross_section_max=0.0)

      
      Function to conduct rejection sampling wrt einstein radius


      :Parameters:

          **param_dict** : `dict`
              dictionary of lens parameters and source parameters

      :Returns:

          **lens_params** : `dict`
              dictionary of lens parameters after rejection sampling













      ..
          !! processed by numpydoc !!

   .. py:method:: rejection_sampling_with_cross_section_sie_feixu(param_dict, cross_section_max=0.0)

      
      Function to conduct rejection sampling wrt cross_section


      :Parameters:

          **param_dict** : `dict`
              dictionary of lens parameters and source parameters

      :Returns:

          **lens_params** : `dict`
              dictionary of lens parameters after rejection sampling













      ..
          !! processed by numpydoc !!

   .. py:method:: rejection_sampling_with_cross_section(param_dict, cross_section_max=0.0)

      
      Function to conduct rejection sampling wrt cross_section of EPL+Shear lens


      :Parameters:

          **param_dict** : `dict`
              dictionary of lens parameters and source parameters

      :Returns:

          **lens_params** : `dict`
              dictionary of lens parameters after rejection sampling













      ..
          !! processed by numpydoc !!

   .. py:method:: rejection_sampling_with_cross_section(param_dict, cross_section_max=0.0)

      
      Function to conduct rejection sampling wrt cross_section, multiprocessing


      :Parameters:

          **param_dict** : `dict`
              dictionary of lens parameters and source parameters

      :Returns:

          **lens_params** : `dict`
              dictionary of lens parameters after rejection sampling













      ..
          !! processed by numpydoc !!


