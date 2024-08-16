:orphan:

:py:mod:`ler.rates`
===================

.. py:module:: ler.rates


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   gwrates/index.rst
   ler/index.rst
   ler copy/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   ler.rates.LensGalaxyParameterDistribution
   ler.rates.LeR
   ler.rates.CBCSourceParameterDistribution
   ler.rates.GWRATES



Functions
~~~~~~~~~

.. autoapisummary::

   ler.rates.load_json
   ler.rates.append_json
   ler.rates.get_param_from_json
   ler.rates.batch_handler
   ler.rates.add_dict_values
   ler.rates.load_json
   ler.rates.append_json
   ler.rates.get_param_from_json
   ler.rates.batch_handler



.. py:class:: LensGalaxyParameterDistribution(npool=4, z_min=0.0, z_max=10.0, cosmology=None, event_type='BBH', lens_type='epl_galaxy', lens_functions=None, lens_priors=None, lens_priors_params=None, directory='./interpolator_pickle', create_new_interpolator=False, **kwargs)


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
           default: 'epl_galaxy'

       **lens_functions, lens_priors, lens_priors_params** : `dict`, `dict`, `dict`
           dictionary of lens functions, priors, and priors parameters
           Check for default/available lens functions, priors and corresponding input parameters by running,

           >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
           >>> lens = LensGalaxyParameterDistribution()
           >>> print(lens.lens_functions)
           >>> print(lens.lens_priors)
           >>> print(lens.lens_priors_params)

       **directory** : `str`
           directory to store the interpolators
           default: './interpolator_pickle'

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
   |:meth:`~sample_all_routine`          | Function to sample galaxy lens   |
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
   |:meth:`~mass_density_spectral_index_normal`                             |
   +-------------------------------------+----------------------------------+
   |                                     | Function to sample the lens      |
   |                                     | galaxy spectral index of the     |
   |                                     | mass density profile from a      |
   |                                     | normal distribution              |
   +-------------------------------------+----------------------------------+
   |:meth:`~compute_einstein_radii`      | Function to compute the Einstein |
   |                                     | radii of the lens galaxies       |
   +-------------------------------------+----------------------------------+
   |:meth:`~rjs_with_cross_section_SIE`  | Function to conduct rejection    |
   |                                     | sampling wrt einstein radius     |
   +-------------------------------------+----------------------------------+
   |:meth:`~rjs_with_cross_section_SIE`  | Function to conduct rejection    |
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
   |:attr:`~sample_mass_density_spectral_index`                             |
   +-------------------------------------+----------------------------------+
   |                                     | Function to sample the lens      |
   |                                     | galaxy spectral index of the     |
   |                                     | mass density profile from a      |
   |                                     | normal distribution              |
   +-------------------------------------+----------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:property:: sample_source_redshift_sl

      
      Function to sample source redshifts conditioned on the source being strongly lensed


      :Parameters:

          **size** : `int`
              number samples to draw

      :Returns:

          **zs** : `numpy.ndarray` (1D array of floats)
              source redshifts conditioned on the source being strongly lensed










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.sample_source_redshift_sl(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:property:: sample_source_parameters

      
      Function to sample source parameters conditioned on the source being strongly lensed


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **source_parameters** : `dict`
              dictionary of source parameters conditioned on the source being strongly lensed













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_lens_redshift

      
      Function to sample lens redshifts, conditioned on the lens being strongly lensed


      :Parameters:

          **zs** : `numpy.ndarray` (1D array of floats)
              source redshifts

      :Returns:

          **zl** : `numpy.ndarray` (1D array of floats)
              lens redshifts corresponding to the source redshifts










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> zs = lens.sample_source_redshift_sl(size=1000)
      >>> lens.sample_lens_redshift(zs=zs)



      ..
          !! processed by numpydoc !!

   .. py:property:: sample_axis_rotation_angle

      
      Function to sample the axis rotation angle of the elliptical lens galaxy from a uniform distribution


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **phi** : `float`
              axis rotation angle of the elliptical lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.sample_axis_rotation_angle(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:property:: sample_shear

      
      Function to sample the elliptical lens galaxy shear from a normal distribution


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **gamma_1** : `float`
              shear component in the x-direction

          **gamma_2** : `float`
              shear component in the y-direction










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> gamma_1, gamma_2 = lens.shear_norm(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:property:: sample_mass_density_spectral_index

      
      Function to sample the lens galaxy spectral index of the mass density profile from a normal distribution


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **gamma** : `float`
              spectral index of the density profile










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.mass_density_spectral_index_normal(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:property:: sample_source_parameters

      
      Function to sample source parameters conditioned on the source being strongly lensed


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **source_parameters** : `dict`
              dictionary of source parameters conditioned on the source being strongly lensed













      ..
          !! processed by numpydoc !!

   .. py:property:: available_lens_prior_list_and_its_params

      
      Dictionary with list all the available priors and it's corresponding parameters. This is an immutable instance attribute.
















      ..
          !! processed by numpydoc !!

   .. py:property:: available_lens_functions

      
      Dictionary with list all the available lens functions. This is an immutable instance attribute.
















      ..
          !! processed by numpydoc !!

   .. py:attribute:: cbc_pop

      
      :class:`~CompactBinaryPopulation` class

      This is an already initialized class that contains a function (CompactBinaryPopulation.sample_gw_parameters) that actually samples the source parameters.















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

   .. py:method:: class_initialization_lens(params=None)

      
      Function to initialize the parent classes


      :Parameters:

          **params** : `dict`
              dictionary of parameters to initialize the parent classes














      ..
          !! processed by numpydoc !!

   .. py:method:: lens_priors_categorization(lens_type, lens_priors=None, lens_priors_params=None, lens_functions=None)

      
      Function to categorize the lens priors/samplers


      :Parameters:

          **lens_type** : `str`
              lens type
              e.g. 'epl_galaxy' for elliptical power-law galaxy

          **lens_priors** : `dict`
              dictionary of priors

          **lens_priors_params** : `dict`
              dictionary of priors parameters

          **lens_functions** : `dict`
              dictionary of lens functions

      :Returns:

          **lens_priors_** : `dict`
              dictionary of priors

          **lens_priors_params_** : `dict`
              dictionary of priors parameters

          **lens_sampler_names_** : `dict`
              dictionary of sampler names

          **lens_functions_** : `dict`
              dictionary of lens functions













      ..
          !! processed by numpydoc !!

   .. py:method:: sample_lens_parameters(size=1000, lens_parameters_input=None)

      
      Function to call the specific galaxy lens parameters sampler routine.
















      ..
          !! processed by numpydoc !!

   .. py:method:: sample_all_routine(size=1000, lens_parameters_input=None)

      
      Function to sample galaxy lens parameters along with the source parameters.


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

              gamma: spectral index of the mass density distribution

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
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.sample_all_routine(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: strongly_lensed_source_redshifts(size=1000)

      
      Function to sample source redshifts and other parameters, conditioned on the source being strongly lensed.


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

   .. py:method:: source_parameters(size, get_attribute=False, param=None)

      
      Function to sample gw source parameters


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **param** : `dict`
              Allows to pass in parameters as dict.
              param =

      :Returns:

          **source_parameters** : `dict`
              Dictionary of source parameters
              source_parameters.keys() = ['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'inclination', 'polarization_angle', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl']










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.source_parameters(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: lens_redshift_SDSS_catalogue(zs, get_attribute=False, param=None)

      
      Function to sample lens redshifts, conditioned on the lens being strongly lensed


      :Parameters:

          **zs** : `float`
              source redshifts

          **get_attribute** : `bool`
              If True, returns a function that can be called with zs as input

      :Returns:

          **zl** : `float`
              lens redshifts










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.lens_redshift_SDSS_catalogue(zs=1.0)



      ..
          !! processed by numpydoc !!

   .. py:method:: axis_rotation_angle_uniform(size=1000, phi_min=0.0, phi_max=2 * np.pi, get_attribute=False, param=None)

      
      Function to sample the axis rotation angle of the elliptical lens galaxy from a uniform distribution.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **phi_min** : `float`
              minimum axis rotation angle of the elliptical lens galaxy

          **phi_max** : `float`
              maximum axis rotation angle of the elliptical lens galaxy

          **get_attribute** : `bool`
              If True, returns a function that can be called with size as input

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(phi_min=0.0, phi_max=2 * np.pi)

      :Returns:

          **phi** : `float`
              axis rotation angle of the elliptical lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.axis_rotation_angle_uniform(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: shear_norm(size, scale=0.05, get_attribute=False, param=None)

      
      Function to sample the elliptical lens galaxy shear from a normal distribution


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **scale** : `float`
              standard deviation of the normal distribution

          **get_attribute** : `bool`
              If True, returns a function that can be called with size as input

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(scale=0.05)

      :Returns:

          **gamma_1** : `float`
              shear component in the x-direction

          **gamma_2** : `float`
              shear component in the y-direction










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> gamma_1, gamma_2 = lens.shear_norm(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: mass_density_spectral_index_normal(size=1000, mean=2.0, std=0.2, get_attribute=False, param=None)

      
      Function to sample the lens galaxy spectral index of the mass density profile from a normal distribution


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **mean** : `float`
              mean of the normal distribution

          **std** : `float`
              standard deviation of the normal distribution

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(mean=2.0, std=0.2)

      :Returns:

          **gamma** : `float`
              spectral index of the density profile










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.mass_density_spectral_index_normal(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: compute_einstein_radii(sigma, zl, zs)

      
      Function to compute the Einstein radii of the lens galaxies


      :Parameters:

          **sigma** : `float`
              velocity dispersion of the lens galaxy

          **zl** : `float`
              lens redshifts

          **zs** : `float`
              source redshifts

      :Returns:

          **theta_E** : `float`
              Einstein radii of the lens galaxies










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> sigma = 200.0
      >>> zl = 0.5
      >>> zs = 1.0
      >>> lens.compute_einstein_radii(sigma, zl, zs)



      ..
          !! processed by numpydoc !!

   .. py:method:: rjs_with_cross_section_SIS(param_dict)

      
      Function to conduct rejection sampling wrt einstein radius


      :Parameters:

          **param_dict** : `dict`
              dictionary of lens parameters and source parameters

      :Returns:

          **lens_params** : `dict`
              dictionary of lens parameters after rejection sampling













      ..
          !! processed by numpydoc !!

   .. py:method:: rjs_with_cross_section_SIE(param_dict)

      
      Function to conduct rejection sampling wrt cross_section


      :Parameters:

          **param_dict** : `dict`
              dictionary of lens parameters and source parameters

      :Returns:

          **lens_params** : `dict`
              dictionary of lens parameters after rejection sampling













      ..
          !! processed by numpydoc !!


.. py:function:: load_json(file_name)

   
   Load a json file.


   :Parameters:

       **file_name** : `str`
           json file name for storing the parameters.

   :Returns:

       **param** : `dict`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: append_json(file_name, new_dictionary, old_dictionary=None, replace=False)

   
   Append (values with corresponding keys) and update a json file with a dictionary. There are four options:

   1. If old_dictionary is provided, the values of the new dictionary will be appended to the old dictionary and save in the 'file_name' json file.
   2. If replace is True, replace the json file (with the 'file_name') content with the new_dictionary.
   3. If the file does not exist, create a new one with the new_dictionary.
   4. If none of the above, append the new dictionary to the content of the json file.

   :Parameters:

       **file_name** : `str`
           json file name for storing the parameters.

       **new_dictionary** : `dict`
           dictionary to be appended to the json file.

       **old_dictionary** : `dict`, optional
           If provided the values of the new dictionary will be appended to the old dictionary and save in the 'file_name' json file.
           Default is None.

       **replace** : `bool`, optional
           If True, replace the json file with the dictionary. Default is False.














   ..
       !! processed by numpydoc !!

.. py:function:: get_param_from_json(json_file)

   
   Function to get the parameters from json file.


   :Parameters:

       **json_file** : `str`
           json file name for storing the parameters.

   :Returns:

       **param** : `dict`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: batch_handler(size, batch_size, sampling_routine, output_jsonfile, save_batch=True, resume=False)

   
   Function to run the sampling in batches.


   :Parameters:

       **size** : `int`
           number of samples.

       **batch_size** : `int`
           batch size.

       **sampling_routine** : `function`
           function to sample the parameters.
           e.g. unlensed_sampling_routine() or lensed_sampling_routine()

       **output_jsonfile** : `str`
           name of the json file to store the parameters.

       **resume** : `bool`
           if True, it will resume the sampling from the last batch.
           default resume = False.














   ..
       !! processed by numpydoc !!

.. py:function:: add_dict_values(dict1, dict2)

   
   Adds the values of two dictionaries together.


   :Parameters:

       **dict1** : `dict`
           dictionary to be added.

       **dict2** : `dict`
           dictionary to be added.

   :Returns:

       **dict1** : `dict`
           dictionary with added values.













   ..
       !! processed by numpydoc !!

.. py:class:: LeR(npool=int(4), z_min=0.0, z_max=10.0, event_type='BBH', size=100000, batch_size=50000, cosmology=None, snr_finder=None, pdet_finder=None, list_of_detectors=None, json_file_names=None, interpolator_directory='./interpolator_pickle', ler_directory='./ler_data', verbose=True, **kwargs)


   Bases: :py:obj:`ler.lens_galaxy_population.LensGalaxyParameterDistribution`

   
   Class to sample of lensed and unlensed events and calculate it's rates. Please note that parameters of the simulated events are stored in json file but not as an attribute of the class. This saves RAM memory.


   :Parameters:

       **npool** : `int`
           number of cores to use.
           default npool = 4.

       **z_min** : `float`
           minimum redshift.
           default z_min = 0.
           for popI_II, popIII, primordial, BNS z_min = 0., 5., 5., 0. respectively.

       **z_max** : `float`
           maximum redshift.
           default z_max = 10.
           for popI_II, popIII, primordial, BNS z_max = 10., 40., 40., 5. respectively.

       **event_type** : `str`
           type of event to generate.
           default event_type = 'BBH'. Other options are 'BNS', 'NSBH'.

       **size** : `int`
           number of samples for sampling.
           default size = 100000. To get stable rates, size should be large (>=1e6).

       **batch_size** : `int`
           batch size for SNR calculation.
           default batch_size = 50000.
           reduce the batch size if you are getting memory error.
           recommended batch_size = 200000, if size = 1000000.

       **cosmology** : `astropy.cosmology`
           cosmology to use for the calculation.
           default cosmology = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7).

       **snr_finder** : `str` or `function`
           default snr_finder = 'gwsnr'.
           if None, the SNR will be calculated using the gwsnr package.
           if custom snr finder function is provided, the SNR will be calculated using a custom function. The custom function should follow the following signature:
           def snr_finder(gw_param_dict):
               ...
               return optimal_snr_dict
           where optimal_snr_dict.keys = ['optimal_snr_net']. Refer to `gwsnr` package's GWSNR.snr attribute for more details.

       **pdet_finder** : `function`
           default pdet_finder = None.
           The rate calculation uses either the pdet_finder or the snr_finder to calculate the detectable events. The custom pdet finder function should follow the following signature:
           def pdet_finder(gw_param_dict):
               ...
               return pdet_net_dict
           where pdet_net_dict.keys = ['pdet_net']. For example uses, refer to [GRB pdet example](https://ler.readthedocs.io/en/latest/examples/rates/grb%20detection%20rate.html).

       **list_of_detectors** : `list`
           list of detectors.
           default list_of_detectors = ['H1', 'L1', 'V1']. This is used for lensed SNR calculation wrt to the detectors. Provide 'None' if you only need net SNR/Pdet. Refer to ImageProperties.get_lensed_snrs for more details.

       **json_file_names: `dict`**
           names of the json files to strore the necessary parameters.
           default json_file_names = {'ler_params': 'LeR_params.json', 'unlensed_param': 'unlensed_param.json', 'unlensed_param_detectable': 'unlensed_param_detectable.json'}.

       **interpolator_directory** : `str`
           directory to store the interpolators.
           default interpolator_directory = './interpolator_pickle'. This is used for storing the various interpolators related to `ler` and `gwsnr` package.

       **ler_directory** : `str`
           directory to store the parameters.
           default ler_directory = './ler_data'. This is used for storing the parameters of the simulated events.

       **verbose** : `bool`
           default verbose = True.
           if True, the function will print all chosen parameters.
           Choose False to prevent anything from printing.

       **kwargs** : `keyword arguments`
           Note : kwargs takes input for initializing the :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution`, :class:`~ler.gw_source_population.CBCSourceParameterDistribution`, :class:`~ler.gw_source_population.CBCSourceRedshiftDistribution` and :class:`~ler.image_properties.ImageProperties` classes. If snr_finder='gwsnr', then kwargs also takes input for initializing the :class:`~gwsnr.GWSNR` class. Please refer to the respective classes for more details.











   .. rubric:: Examples

   >>> from ler.rates import LeR
   >>> ler = LeR()
   >>> unlensed_params = ler.unlensed_cbc_statistics();
   >>> ler.unlensed_rate();
   >>> lensed_params = ler.lensed_cbc_statistics();
   >>> ler.lensed_rate();
   >>> ler.rate_ratio();

   Instance Attributes
   ----------
   LeR class has the following attributes,

   +-------------------------------------+----------------------------------+
   | Atrributes                          | Type                             |
   +=====================================+==================================+
   |:attr:`~npool`                       | `int`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~z_min`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~z_max`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~event_type`                  | `str`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~cosmo`                       | `astropy.cosmology`              |
   +-------------------------------------+----------------------------------+
   |:attr:`~size`                        | `int`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~batch_size`                  | `int`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~json_file_names`             | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~interpolator_directory`      | `str`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~ler_directory`               | `str`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~gwsnr`                       | `bool`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~gw_param_sampler_dict`       | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~snr_calculator_dict`         | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~list_of_detectors`           | `list`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~unlensed_param`              | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~unlensed_param_detectable`   | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~lensed_param`                | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~lensed_param_detectable`     | `dict`                           |
   +-------------------------------------+----------------------------------+

   Instance Methods
   ----------
   LeR class has the following methods,

   +-------------------------------------+----------------------------------+
   | Methods                             | Description                      |
   +=====================================+==================================+
   |:meth:`~class_initialization`        | Function to initialize the       |
   |                                     | parent classes                   |
   +-------------------------------------+----------------------------------+
   |:meth:`~gwsnr_intialization`         | Function to initialize the       |
   |                                     | gwsnr class                      |
   +-------------------------------------+----------------------------------+
   |:meth:`~snr`                         | Function to get the snr with the |
   |                                     | given parameters.                |
   +-------------------------------------+----------------------------------+
   |:meth:`~snr_bilby`                   | Function to get the snr with the |
   |                                     | given parameters using inner-    |
   |                                     | product method.                  |
   +-------------------------------------+----------------------------------+
   |:meth:`~pdet`                        | Function to get the pdet with    |
   |                                     | the given parameters.            |
   +-------------------------------------+----------------------------------+
   |:meth:`~store_ler_params`            | Function to store the all the    |
   |                                     | necessary parameters.            |
   +-------------------------------------+----------------------------------+
   |:meth:`~unlensed_cbc_statistics`     | Function to generate unlensed    |
   |                                     | GW source parameters in batches. |
   +-------------------------------------+----------------------------------+
   |:meth:`~unlensed_sampling_routine`   | Function to generate unlensed    |
   |                                     | GW source parameters. It stores  |
   |                                     | the parameters of the generated  |
   |                                     | events in a json file.           |
   +-------------------------------------+----------------------------------+
   |:meth:`~unlensed_rate`               | Function to calculate the        |
   |                                     | unlensed rate. It also stores    |
   |                                     | the parameters of the detectable |
   |                                     | unlesed events in a json file.   |
   +-------------------------------------+----------------------------------+
   |:meth:`~lensed_cbc_statistics`       | Function to generate lensed      |
   |                                     | GW source parameters.            |
   +-------------------------------------+----------------------------------+
   |:meth:`~lensed_sampling_routine`     | Function to generate lensed      |
   |                                     | GW source parameters. It stores  |
   |                                     | the parameters of the generated  |
   |                                     | events in a json file.           |
   +-------------------------------------+----------------------------------+
   |:meth:`~lensed_rate`                 | Function to calculate the        |
   |                                     | lensed rate. It also stores the  |
   |                                     | parameters of the detectable     |
   |                                     | lensed events in a json file.    |
   +-------------------------------------+----------------------------------+
   |:meth:`~rate_ratio`                  | Function to calculate the rate   |
   |                                     | ratio between lensed and         |
   |                                     | unlensed events.                 |
   +-------------------------------------+----------------------------------+
   |:meth:`~rate_comparision_with_rate_calculation                          |
   +-------------------------------------+----------------------------------+
   |                                     | Function to calculate rates for  |
   |                                     | unleesed and lensed events and   |
   |                                     | compare it with the rate. It also|
   |                                     | stores the parameters of the     |
   |                                     | detectable events in a json file.|
   +-------------------------------------+----------------------------------+
   |:meth:`~selecting_n_unlensed_detectable_events`                         |
   +-------------------------------------+----------------------------------+
   |                                     | Function to select n unlensed    |
   |                                     | detectable events. It stores the |
   |                                     | parameters of the detectable     |
   |                                     | unlesed events in a json file.   |
   +-------------------------------------+----------------------------------+
   |:meth:`~selecting_n_lensed_detectable_events`                           |
   +-------------------------------------+----------------------------------+
   |                                     | Function to select n lensed      |
   |                                     | detectable events. It stores the |
   |                                     | parameters of the detectable     |
   |                                     | lensed events in a json file.    |
   +-------------------------------------+----------------------------------+

   Note: `LeR` class also inherits all the instances from the :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution` class. Please refer to the :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution` class for more details.



   ..
       !! processed by numpydoc !!
   .. py:property:: snr

      
      Function to get the snr with the given parameters.


      :Parameters:

          **gw_param_dict** : `dict`
              dictionary of GW source parameters.
              mass_1 : `numpy.ndarray` or `float`
                  mass_1 of the compact binary (detector frame) (Msun).
              mass_2 : `numpy.ndarray` or `float`
                  mass_2 of the compact binary (detector frame) (Msun).
              luminosity_distance : `numpy.ndarray` or `float`
                  luminosity distance of the source (Mpc).
              theta_jn : `numpy.ndarray` or `float`
                  inclination angle of the source (rad).
              psi : `numpy.ndarray` or `float`
                  polarization angle of the source (rad).
              phase : `numpy.ndarray` or `float`
                  phase of GW at reference frequency  (rad).
              geocent_time : `numpy.ndarray` or `float`
                  GPS time of coalescence (s).
              ra : `numpy.ndarray` or `float`
                  right ascension of the source (rad).
              dec : `numpy.ndarray` or `float`
                  declination of the source (rad).
              a_1 : `numpy.ndarray` or `float`
                  dimensionless spin magnitude of the more massive object.
              a_2 : `numpy.ndarray` or `float`
                  dimensionless spin magnitude of the less massive object.
              tilt_1 : `numpy.ndarray` or `float`
                  tilt angle of the more massive object spin.
              tilt_2 : `numpy.ndarray` or `float`
                  tilt angle of the less massive object spin.
              phi_12 : `numpy.ndarray` or `float`
                  azimuthal angle between the two spin vectors.
              phi_jl : `numpy.ndarray` or `float`
                  azimuthal angle between total angular momentum and the orbital angular momentum.

      :Returns:

          **optimal_snr_list** : `list`
              e.g. [optimal_snr_net, 'L1', 'H1', 'V1']
              optimal_snr_net : `numpy.ndarray` or `float`
                  optimal snr of the network.
              'H1' : `numpy.ndarray` or `float`
                  optimal snr of H1.
              'L1' : `numpy.ndarray` or `float`
                  optimal snr of L1.
              'V1' : `numpy.ndarray` or `float`
                  optimal snr of V1.













      ..
          !! processed by numpydoc !!

   .. py:property:: unlensed_param

      
      Function to get data from the json file self.json_file_names["unlensed_param"].



      :Returns:

          **unlensed_param** : `dict`
              dictionary of unlensed GW source parameters.













      ..
          !! processed by numpydoc !!

   .. py:property:: unlensed_param_detectable

      
      Function to get data from the json file self.json_file_names["unlensed_param_detectable"].



      :Returns:

          **unlensed_param_detectable** : `dict`
              dictionary of unlensed GW source parameters.













      ..
          !! processed by numpydoc !!

   .. py:property:: lensed_param

      
      Function to get data from the json file self.json_file_names["lensed_param"].



      :Returns:

          **lensed_param** : `dict`
              dictionary of lensed GW source parameters.













      ..
          !! processed by numpydoc !!

   .. py:property:: lensed_param_detectable

      
      Function to get data from the json file self.json_file_names["lensed_param_detectable"].



      :Returns:

          **lensed_param_detectable** : `dict`
              dictionary of lensed GW source parameters.













      ..
          !! processed by numpydoc !!

   .. py:attribute:: npool

      
      ``int``

      Number of logical cores to use.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: z_min

      
      ``float``

      Minimum redshift of the source population















      ..
          !! processed by numpydoc !!

   .. py:attribute:: z_max

      
      ``float``

      Maximum redshift of the source population















      ..
          !! processed by numpydoc !!

   .. py:attribute:: event_type

      
      ``str``

      Type of event to generate.

      e.g. 'BBH', 'BNS', 'NSBH'















      ..
          !! processed by numpydoc !!

   .. py:attribute:: cosmo

      
      ``astropy.cosmology``

      Cosmology to use for the calculation.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: size

      
      ``int``

      Number of samples for sampling.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: batch_size

      
      ``int``

      Batch size for sampling.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: json_file_names

      
      ``dict``

      Names of the json files to store the necessary parameters.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: interpolator_directory

      
      ``str``

      Directory to store the interpolators.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: ler_directory

      
      ``str``

      Directory to store the parameters.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: gwsnr

      
      ``bool``

      If True, the SNR will be calculated using the gwsnr package.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: gw_param_sampler_dict

      
      ``dict``

      Dictionary of parameters to initialize the ``CBCSourceParameterDistribution`` class.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: snr_calculator_dict

      
      ``dict``

      Dictionary of parameters to initialize the ``GWSNR`` class.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: list_of_detectors

      
      ``list``

      List of detectors.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: unlensed_param

      
      ``dict``

      Dictionary of unlensed GW source parameters. The included parameters and their units are as follows (for default settings):

      +--------------------+--------------+--------------------------------------+
      | Parameter          | Units        | Description                          |
      +====================+==============+======================================+
      | zs                 |              | redshift of the source               |
      +--------------------+--------------+--------------------------------------+
      | geocent_time       | s            | GPS time of coalescence              |
      +--------------------+--------------+--------------------------------------+
      | ra                 | rad          | right ascension                      |
      +--------------------+--------------+--------------------------------------+
      | dec                | rad          | declination                          |
      +--------------------+--------------+--------------------------------------+
      | phase              | rad          | phase of GW at reference frequency   |
      +--------------------+--------------+--------------------------------------+
      | psi                | rad          | polarization angle                   |
      +--------------------+--------------+--------------------------------------+
      | theta_jn           | rad          | inclination angle                    |
      +--------------------+--------------+--------------------------------------+
      | luminosity_distance| Mpc          | luminosity distance                  |
      +--------------------+--------------+--------------------------------------+
      | mass_1_source      | Msun         | mass_1 of the compact binary         |
      |                    |              | (source frame)                       |
      +--------------------+--------------+--------------------------------------+
      | mass_2_source      | Msun         | mass_2 of the compact binary         |
      |                    |              | (source frame)                       |
      +--------------------+--------------+--------------------------------------+
      | mass_1             | Msun         | mass_1 of the compact binary         |
      |                    |              | (detector frame)                     |
      +--------------------+--------------+--------------------------------------+
      | mass_2             | Msun         | mass_2 of the compact binary         |
      |                    |              | (detector frame)                     |
      +--------------------+--------------+--------------------------------------+
      | L1                 |              | optimal snr of L1                    |
      +--------------------+--------------+--------------------------------------+
      | H1                 |              | optimal snr of H1                    |
      +--------------------+--------------+--------------------------------------+
      | V1                 |              | optimal snr of V1                    |
      +--------------------+--------------+--------------------------------------+
      | optimal_snr_net    |              | optimal snr of the network           |
      +--------------------+--------------+--------------------------------------+















      ..
          !! processed by numpydoc !!

   .. py:attribute:: unlensed_param_detectable

      
      ``dict``

      Dictionary of detectable unlensed GW source parameters. It includes the same parameters as the :attr:`~unlensed_param` attribute.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: lensed_param

      
      ``dict``

      Dictionary of lens parameters, images parameters and lensed GW source parameters. The included parameters and their units are as follows (for default settings):

      +------------------------------+-----------+-------------------------------+
      | Parameter                    | Units     | Description                   |
      +==============================+===========+===============================+
      | zl                           |           | redshift of the lens          |
      +------------------------------+-----------+-------------------------------+
      | zs                           |           | redshift of the source        |
      +------------------------------+-----------+-------------------------------+
      | sigma                        |km s^-1    | velocity dispersion           |
      +------------------------------+-----------+-------------------------------+
      | q                            |           | axis ratio                    |
      +------------------------------+-----------+-------------------------------+
      | theta_E                      | arcsec    | Einstein radius               |
      +------------------------------+-----------+-------------------------------+
      | phi                          | rad       | axis rotation angle           |
      +------------------------------+-----------+-------------------------------+
      | e1                           |           | ellipticity component 1       |
      +------------------------------+-----------+-------------------------------+
      | e2                           |           | ellipticity component 2       |
      +------------------------------+-----------+-------------------------------+
      | gamma1                       |           | shear component 1             |
      +------------------------------+-----------+-------------------------------+
      | gamma2                       |           | shear component 2             |
      +------------------------------+-----------+-------------------------------+
      | gamma                        |           | shear                         |
      +------------------------------+-----------+-------------------------------+
      | ra                           | rad       | right ascension               |
      +------------------------------+-----------+-------------------------------+
      | dec                          | rad       | declination                   |
      +------------------------------+-----------+-------------------------------+
      | phase                        | rad       | phase of GW at reference freq |
      +------------------------------+-----------+-------------------------------+
      | psi                          | rad       | polarization angle            |
      +------------------------------+-----------+-------------------------------+
      | theta_jn                     | rad       | inclination angle             |
      +------------------------------+-----------+-------------------------------+
      | mass_1_source                | Msun      | mass_1 of the compact binary  |
      |                              |           | (source frame)                |
      +------------------------------+-----------+-------------------------------+
      | mass_2_source                | Msun      | mass_2 of the compact binary  |
      |                              |           | (source frame)                |
      +------------------------------+-----------+-------------------------------+
      | mass_1                       | Msun      | mass_1 of the compact binary  |
      |                              |           | (detector frame)              |
      +------------------------------+-----------+-------------------------------+
      | mass_2                       | Msun      | mass_2 of the compact binary  |
      |                              |           | (detector frame)              |
      +------------------------------+-----------+-------------------------------+
      | x0_image_positions           |           | x0 image positions            |
      +------------------------------+-----------+-------------------------------+
      | x1_image_positions           |           | x1 image positions            |
      +------------------------------+-----------+-------------------------------+
      | magnifications               |           | magnifications                |
      +------------------------------+-----------+-------------------------------+
      | time_delays                  |           | time delays                   |
      +------------------------------+-----------+-------------------------------+
      | image_type                   |           | image type                    |
      +------------------------------+-----------+-------------------------------+
      | n_images                     |           | number of images              |
      +------------------------------+-----------+-------------------------------+
      | effective_luminosity_distance| Mpc       | effective luminosity distance |
      +------------------------------+-----------+-------------------------------+
      | effective_geocent_time       | s         | effective GPS time of coalesc |
      +------------------------------+-----------+-------------------------------+
      | L1                           |           | optimal snr of L1             |
      +------------------------------+-----------+-------------------------------+
      | H1                           |           | optimal snr of H1             |
      +------------------------------+-----------+-------------------------------+
      | V1                           |           | optimal snr of V1             |
      +------------------------------+-----------+-------------------------------+
      | optimal_snr_net              |           | optimal snr of the network    |
      +------------------------------+-----------+-------------------------------+















      ..
          !! processed by numpydoc !!

   .. py:attribute:: lensed_param_detectable

      
      ``dict``

      Dictionary of detectable lensed GW source parameters.















      ..
          !! processed by numpydoc !!

   .. py:method:: print_all_params()

      
      Function to print all the parameters.
















      ..
          !! processed by numpydoc !!

   .. py:method:: class_initialization(params=None)

      
      Function to initialize the parent classes.


      :Parameters:

          **params** : `dict`
              dictionary of parameters to initialize the parent classes














      ..
          !! processed by numpydoc !!

   .. py:method:: gwsnr_intialization(params=None)

      
      Function to initialize the GWSNR class from the `gwsnr` package.


      :Parameters:

          **params** : `dict`
              dictionary of parameters to initialize the gwsnr class














      ..
          !! processed by numpydoc !!

   .. py:method:: store_ler_params(output_jsonfile='ler_params.json')

      
      Function to store the all the necessary parameters. This is useful for reproducing the results. All the parameters stored are in string format to make it json compatible.


      :Parameters:

          **output_jsonfile** : `str`
              name of the json file to store the parameters














      ..
          !! processed by numpydoc !!

   .. py:method:: unlensed_cbc_statistics(size=None, resume=False, save_batch=False, output_jsonfile=None)

      
      Function to generate unlensed GW source parameters. This function calls the unlensed_sampling_routine function to generate the parameters in batches. The generated parameters are stored in a json file; and if save_batch=True, it keeps updating the file in batches.


      :Parameters:

          **size** : `int`
              number of samples.
              default size = 100000.

          **resume** : `bool`
              resume = False (default) or True.
              if True, the function will resume from the last batch.

          **save_batch** : `bool`
              if True, the function will save the parameters in batches. if False, the function will save all the parameters at the end of sampling. save_batch=False is faster.

          **output_jsonfile** : `str`
              json file name for storing the parameters.
              default output_jsonfile = 'unlensed_params.json'. Note that this file will be stored in the self.ler_directory.

      :Returns:

          **unlensed_param** : `dict`
              dictionary of unlensed GW source parameters. Refer to :attr:`~unlensed_param` for details.










      .. rubric:: Examples

      >>> from ler.rates import LeR
      >>> ler = LeR()
      >>> unlensed_param = ler.unlensed_cbc_statistics()



      ..
          !! processed by numpydoc !!

   .. py:method:: unlensed_sampling_routine(size, output_jsonfile, resume=False, save_batch=True)

      
      Function to generate unlensed GW source parameters. This function also stores the parameters in json file in the current batch if save_batch=True.


      :Parameters:

          **size** : `int`
              number of samples.
              default size = 100000.

          **output_jsonfile** : `str`
              json file name for storing the parameters.
              default output_jsonfile = 'unlensed_params.json'. Note that this file will be stored in the self.ler_directory.

          **resume** : `bool`
              resume = False (default) or True.
              if True, it appends the new samples to the existing json file.

          **save_batch** : `bool`
              if True, the function will save the parameters in batches. if False, the function will save all the parameters at the end of sampling. save_batch=False is faster.

      :Returns:

          **unlensed_param** : `dict`
              dictionary of unlensed GW source parameters. Refer to :attr:`~unlensed_param` for details.













      ..
          !! processed by numpydoc !!

   .. py:method:: unlensed_rate(unlensed_param=None, snr_threshold=8.0, pdet_threshold=0.5, output_jsonfile=None, detectability_condition='step_function', snr_recalculation=False, snr_threshold_recalculation=5.5)

      
      Function to calculate the unlensed rate. This function also stores the parameters of the detectable events in json file. There are two conditions for detectability: 'step_function' and 'pdet'.

      1. 'step_function': If two images have SNR>8.0, then the event is detectable. This is a step function. This is with the assumption that SNR function is provided and not None.
      2. 'pdet':
          i) If self.pdet is None and self.snr is not None, then it will calculate the pdet from the snr. There is no hard cut for this pdet and can have value ranging from 0 to 1 near the threshold.
          ii) If self.pdet is not None, then it will use the generated pdet.

      :Parameters:

          **unlensed_param** : `dict` or `str`
              dictionary of GW source parameters or json file name.
              default unlensed_param = 'unlensed_params.json'.

          **snr_threshold** : `float`
              threshold for detection signal to noise ratio.
              e.g. snr_threshold = 8.

          **pdet_threshold** : `float`
              threshold for detection probability.
              e.g. pdet_threshold = 0.5.

          **output_jsonfile** : `str`
              json file name for storing the parameters of the detectable events.
              default output_jsonfile = 'unlensed_params_detectable.json'.

          **detectability_condition** : `str`
              detectability condition.
              default detectability_condition = 'step_function'.
              other options are 'pdet'.

          **snr_recalculation** : `bool`
              if True, the SNR of centain events (snr>snr_threshold_recalculation)will be recalculate with 'inner-product' method. This is useful when the snr is calculated with 'ann' method of `gwsnr`.
              default snr_recalculation = False.

          **snr_threshold_recalculation** : `float`
              threshold for recalculation of detection signal to noise ratio.

      :Returns:

          **total_rate** : `float`
              total unlensed rate (Mpc^-3 yr^-1).

          **unlensed_param** : `dict`
              dictionary of unlensed GW source parameters of the detectable events. Refer to :attr:`~unlensed_param` for details.










      .. rubric:: Examples

      >>> from ler.rates import LeR
      >>> ler = LeR()
      >>> ler.unlensed_cbc_statistics();
      >>> total_rate, unlensed_param_detectable = ler.unlensed_rate()



      ..
          !! processed by numpydoc !!

   .. py:method:: rate_function(detectable_size, total_size, param_type='unlensed', verbose=True)

      
      General helper function to calculate the rate for unlensed and lensed events.


      :Parameters:

          **detectable_size** : `int`
              number of detectable events.

          **total_size** : `int`
              total number of events.

          **param_type** : `str`
              type of parameters.

      :Returns:

          **rate** : `float`
              rate of the events.










      .. rubric:: Examples

      >>> from ler.rates import LeR
      >>> ler = LeR()
      >>> rate = ler.rate_function(detectable_size=100, total_size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: lensed_cbc_statistics(size=None, save_batch=False, resume=False, output_jsonfile=None)

      
      Function to generate lensed GW source parameters. This function calls the lensed_sampling_routine function to generate the parameters in batches. The generated parameters are stored in a json file; and if save_batch=True, it keeps updating the file in batches.


      :Parameters:

          **size** : `int`
              number of samples.
              default size = 100000.

          **save_batch** : `bool`
              if True, the function will save the parameters in batches. if False, the function will save all the parameters at the end of sampling. save_batch=False is faster.

          **resume** : `bool`
              resume = False (default) or True.
              if True, the function will resume from the last batch.

          **output_jsonfile** : `str`
              json file name for storing the parameters.
              default output_jsonfile = 'lensed_params.json'.

      :Returns:

          **lensed_param** : `dict`
              dictionary of lensed GW source parameters. Refer to :attr:`~lensed_param` for details.










      .. rubric:: Examples

      >>> from ler.rates import LeR
      >>> ler = LeR()
      >>> lensed_param = ler.lensed_cbc_statistics()



      ..
          !! processed by numpydoc !!

   .. py:method:: lensed_sampling_routine(size, output_jsonfile, save_batch=True, resume=False)

      
      Function to generate lensed GW source parameters. This function also stores the parameters in json file in the current batch if save_batch=True.


      :Parameters:

          **size** : `int`
              number of samples.
              default size = 100000.

          **output_jsonfile** : `str`
              json file name for storing the parameters.
              default output_jsonfile = 'lensed_params.json'. Note that this file will be stored in the self.ler_directory.

          **save_batch** : `bool`
              if True, the function will save the parameters in batches. if False, the function will save all the parameters at the end of sampling. save_batch=False is faster.

          **resume** : `bool`
              resume = False (default) or True.
              if True, it appends the new samples to the existing json file.

      :Returns:

          **lensed_param** : `dict`
              dictionary of lensed GW source parameters. Refer to :attr:`~lensed_param` for details.













      ..
          !! processed by numpydoc !!

   .. py:method:: lensed_rate(lensed_param=None, snr_threshold=[8.0, 8.0], pdet_threshold=0.5, num_img=[1, 1], output_jsonfile=None, nan_to_num=True, detectability_condition='step_function', snr_recalculation=False, snr_threshold_recalculation=[5.5, 5.5])

      
      Function to calculate the lensed rate. This function also stores the parameters of the detectable events in json file. There are two conditions for detectability: 'step_function' and 'pdet'.

      1. 'step_function': If two images have SNR>8.0, then the event is detectable. This is a step function. This is with the assumption that SNR function is provided and not None.
      2. 'pdet':
          i) If self.pdet is None and self.snr is not None, then it will calculate the pdet from the snr. There is no hard cut for this pdet and can have value ranging from 0 to 1 near the threshold.
          ii) If self.pdet is not None, then it will use the generated pdet.

      :Parameters:

          **lensed_param** : `dict` or `str`
              dictionary of GW source parameters or json file name.
              default lensed_param = 'lensed_params.json'.

          **snr_threshold** : `float`
              threshold for detection signal to noise ratio. This is use when self.snr is provided.
              default snr_threshold = [8.0,8.0].

          **pdet_threshold** : `float`
              threshold for detection probability. This is use when self.pdet is provided.
              default pdet_threshold = 0.5.

          **num_img** : `int`
              number of images corresponding to the snr_threshold.
              default num_img = [1,1]. Together with snr_threshold = [8.0,8.0], it means that two images with snr>8.0. Same condition can also be represented by snr_threshold = 8.0 and num_img = 2.

          **output_jsonfile** : `str`
              json file name for storing the parameters of the detectable events.
              default output_jsonfile = 'lensed_params_detectable.json'.

          **nan_to_num** : `bool`
              if True, nan values will be converted to 0.
              default nan_to_num = True.

          **detectability_condition** : `str`
              detectability condition.
              default detectability_condition = 'step_function'.
              other options are 'pdet'.

          **snr_recalculation** : `bool`
              if True, the SNR of centain events (snr>snr_threshold_recalculation)will be recalculate with 'inner-product' method. This is useful when the snr is calculated with 'ann' method of `gwsnr`.
              default snr_recalculation = False.

          **snr_threshold_recalculation** : `float`
              threshold for recalculation of detection signal to noise ratio.
              default snr_threshold_recalculation = [5.5,5.5].

      :Returns:

          **total_rate** : `float`
              total lensed rate (Mpc^-3 yr^-1).

          **lensed_param** : `dict`
              dictionary of lensed GW source parameters of the detectable events. Refer to :attr:`~lensed_param` for details.










      .. rubric:: Examples

      >>> from ler.rates import LeR
      >>> ler = LeR()
      >>> ler.lensed_cbc_statistics();
      >>> total_rate, lensed_param_detectable = ler.lensed_rate()



      ..
          !! processed by numpydoc !!

   .. py:method:: rate_comparision_with_rate_calculation(unlensed_param=None, snr_threshold_unlensed=8.0, output_jsonfile_unlensed=None, lensed_param=None, snr_threshold_lensed=[8.0, 8.0], num_img=[1, 1], output_jsonfile_lensed=None, nan_to_num=True, detectability_condition='step_function')

      
      Function to calculate the unlensed and lensed rate and compare by computing the ratio. This function also stores the parameters of the detectable events in json file. If you use this function, you do not need to call the functions unlensed_rate and lensed_rate separately.


      :Parameters:

          **unlensed_param** : `dict` or `str`
              dictionary of GW source parameters or json file name.
              default unlensed_param = 'unlensed_params.json'.

          **snr_threshold_unlensed** : `float`
              threshold for detection signal to noise ratio.
              e.g. snr_threshold_unlensed = 8.

          **output_jsonfile_unlensed** : `str`
              json file name for storing the parameters of the detectable events.
              default output_jsonfile_unlensed = 'unlensed_params_detectable.json'.

          **lensed_param** : `dict` or `str`
              dictionary of GW source parameters or json file name.
              default lensed_param = 'lensed_params.json'.

          **snr_threshold_lensed** : `float`
              threshold for detection signal to noise ratio.
              default snr_threshold_lensed = [8.0,8.0].

          **num_img** : `int`
              number of images.
              default num_img = [1,1]. Together with snr_threshold = [8.0,8.0], it means that two images with snr>8.0. Same condition can also be represented by snr_threshold = 8.0 and num_img = 2.

          **output_jsonfile_lensed** : `str`
              json file name for storing the parameters of the detectable events.
              default output_jsonfile_lensed = 'lensed_params_detectable.json'.

          **nan_to_num** : `bool`
              if True, nan values will be converted to 0.
              default nan_to_num = True.

          **detectability_condition** : `str`
              detectability condition.
              default detectability_condition = 'step_function'.
              other options are 'pdet'.

      :Returns:

          **rate_ratio** : `float`
              rate ratio.

          **unlensed_param** : `dict`
              dictionary of unlensed GW source parameters of the detectable events. Refer to :attr:`~unlensed_param` for details.

          **lensed_param** : `dict`
              dictionary of lensed GW source parameters of the detectable events. Refer to :attr:`~lensed_param` for details.










      .. rubric:: Examples

      >>> from ler.rates import LeR
      >>> ler = LeR()
      >>> ler.unlensed_cbc_statistics();
      >>> ler.lensed_cbc_statistics();
      >>> rate_ratio, unlensed_param, lensed_param = ler.rate_comparision_with_rate_calculation()



      ..
          !! processed by numpydoc !!

   .. py:method:: rate_ratio()

      
      Function to calculate and display unlensed and lensed merger rate ratio. It will get the unlensed_rate and lensed_rate from files corresponding to the names included in self.json_file_ler_param.



      :Returns:

          **rate_ratio** : `float`
              rate ratio.










      .. rubric:: Examples

      >>> from ler.rates import LeR
      >>> ler = LeR()
      >>> ler.unlensed_cbc_statistics();
      >>> ler.lensed_cbc_statistics();
      >>> ler.unlensed_rate();
      >>> ler.lensed_rate();
      >>> ler.rate_ratio()



      ..
          !! processed by numpydoc !!

   .. py:method:: selecting_n_unlensed_detectable_events(size=100, batch_size=None, snr_threshold=8.0, pdet_threshold=0.5, resume=False, output_jsonfile='n_unlensed_param_detectable.json', meta_data_file='meta_unlensed.json', detectability_condition='step_function', trim_to_size=True, snr_recalculation=False, snr_threshold_recalculation=5.5)

      
      Function to generate n unlensed detectable events. This fuction samples the unlensed parameters and save only the detectable events in json file. It also records metadata in the JSON file, which includes the total number of events and the cumulative rate of events. This functionality is particularly useful for generating a fixed or large number of detectable events until the event rates stabilize.


      :Parameters:

          **size** : `int`
              number of samples to be selected.
              default size = 100.

          **batch_size** : `int`
              batch size for sampling.
              default batch_size = 50000.

          **snr_threshold** : `float`
              threshold for detection signal to noise ratio.
              e.g. snr_threshold = 8.

          **pdet_threshold** : `float`
              threshold for detection probability.
              default pdet_threshold = 0.5.

          **resume** : `bool`
              resume = False (default) or True.
              if True, the function will resume from the last batch.

          **output_jsonfile** : `str`
              json file name for storing the parameters of the detectable events.
              default output_jsonfile = 'n_unlensed_param_detectable.json'.

          **meta_data_file** : `str`
              json file name for storing the metadata.
              default meta_data_file = 'meta_unlensed.json'.

          **detectability_condition** : `str`
              detectability condition.
              default detectability_condition = 'step_function'.
              other options are 'pdet'.

          **trim_to_size** : `bool`
              if True, the final result will be trimmed to size.
              default trim_to_size = True.

          **snr_recalculation** : `bool`
              if True, the SNR of centain events (snr>snr_threshold_recalculation)will be recalculate with 'inner-product' method. This is useful when the snr is calculated with 'ann' method of `gwsnr`.
              default snr_recalculation = False.

          **snr_threshold_recalculation** : `float`
              threshold for recalculation of detection signal to noise ratio.
              default snr_threshold_recalculation = 5.5.

      :Returns:

          **param_final** : `dict`
              dictionary of unlensed GW source parameters of the detectable events. Refer to :attr:`~unlensed_param` for details.










      .. rubric:: Examples

      >>> from ler.rates import LeR
      >>> ler = LeR()
      >>> unlensed_param = ler.selecting_n_unlensed_detectable_events(size=100)



      ..
          !! processed by numpydoc !!

   .. py:method:: selecting_n_lensed_detectable_events(size=100, batch_size=None, snr_threshold=[8.0, 8.0], pdet_threshold=0.5, num_img=[1, 1], resume=False, detectability_condition='step_function', output_jsonfile='n_lensed_params_detectable.json', meta_data_file='meta_lensed.json', trim_to_size=True, nan_to_num=False, snr_recalculation=False, snr_threshold_recalculation=[5.5, 5.5])

      
      Function to generate n lensed detectable events. This fuction only samples the lensed parameters and save only the detectable events in json file. It also records metadata in the JSON file, which includes the total number of events and the cumulative rate of events. This functionality is particularly useful for generating a fixed or large number of detectable events until the event rates stabilize.


      :Parameters:

          **size** : `int`
              number of samples.
              default size = 100.

          **batch_size** : `int`
              batch size for sampling.
              default batch_size = 50000.

          **snr_threshold** : `float`
              threshold for detection signal to noise ratio.
              default snr_threshold = [8.0,8.0].

          **pdet_threshold** : `float`
              threshold for detection probability.
              default pdet_threshold = 0.5.

          **num_img** : `int`
              number of images.
              default num_img = [1,1]. Together with snr_threshold = [8.0,8.0], it means that two images with snr>8.0. Same condition can also be represented by snr_threshold = 8.0 and num_img = 2.

          **resume** : `bool`
              resume = False (default) or True.
              if True, it appends the new samples to the existing json file.

          **detectability_condition** : `str`
              detectability condition.
              default detectability_condition = 'step_function'.
              other options are 'pdet'.

          **output_jsonfile** : `str`
              json file name for storing the parameters of the detectable events.
              default output_jsonfile = 'n_lensed_params_detectable.json'.

          **meta_data_file** : `str`
              json file name for storing the metadata.
              default meta_data_file = 'meta_lensed.json'.

          **trim_to_size** : `bool`
              if True, the final result will be trimmed to size.
              default trim_to_size = True.

          **nan_to_num** : `bool`
              if True, nan values will be converted to 0.
              default nan_to_num = False.

          **snr_recalculation** : `bool`
              if True, the SNR of centain events (snr>snr_threshold_recalculation)will be recalculate with 'inner-product' method. This is useful when the snr is calculated with 'ann' method of `gwsnr`.
              default snr_recalculation = False.

          **snr_threshold_recalculation** : `float`
              threshold for recalculation of detection signal to noise ratio.
              default snr_threshold_recalculation = [5.5,5.5].

      :Returns:

          **param_final** : `dict`
              dictionary of lensed GW source parameters of the detectable events. Refer to :attr:`~lensed_param` for details.










      .. rubric:: Examples

      >>> from ler.rates import LeR
      >>> ler = LeR()
      >>> lensed_param = ler.selecting_n_lensed_detectable_events(size=100)



      ..
          !! processed by numpydoc !!


.. py:class:: CBCSourceParameterDistribution(z_min=0.0, z_max=10.0, event_type='BBH', source_priors=None, source_priors_params=None, cosmology=None, spin_zero=True, spin_precession=False, directory='./interpolator_pickle', create_new_interpolator=False)


   Bases: :py:obj:`ler.gw_source_population.cbc_source_redshift_distribution.CBCSourceRedshiftDistribution`

   
   Class to generate a population of compact binaries. It helps sample all the intrinsic and extrinsic parameters of compact binaries. This daughter class inherits from :class:`~ler.ler.CBCSourceRedshiftDistribution` class.


   :Parameters:

       **z_min** : `float`
           Minimum redshift of the source population
           default: 0.001

       **z_max** : `float`
           Maximum redshift of the source population
           default: 10.

       **event_type** : `str`
           Type of event to generate.
           e.g. 'BBH', 'BNS', 'NSBH'

       **source_priors, source_priors_params** : `dict`, `dict`
           Dictionary of prior sampler functions and its input parameters.
           Check for available priors and corresponding input parameters by running,
           >>> from ler.gw_source_population import CBCSourceParameterDistribution
           >>> cbc = CompactBinaryPopulation()
           >>> cbc.available_gw_prior_list_and_its_params()
           # To check the current chosen priors and its parameters, run,
           >>> print("default priors=",cbc.gw_param_samplers)
           >>> print("default priors's parameters=",cbc.gw_param_samplers_params)

       **cosmology** : `astropy.cosmology`
           Cosmology to use
           default: None/astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)

       **spin_zero** : `bool`
           If True, spin parameters are completely ignore in the sampling.
           default: True

       **spin_precession** : `bool`
           If spin_zero=True and spin_precession=True, spin parameters are sampled for precessing binaries.
           if spin_zero=True and spin_precession=False, spin parameters are sampled for aligned/anti-aligned spin binaries.
           default: False

       **directory** : `str`
           Directory to store the interpolator pickle files
           default: './interpolator_pickle'

       **create_new_interpolator** : `dict`
           Dictionary of boolean values and resolution to create new interpolator.
           default: dict(redshift_distribution=dict(create_new=False, resolution=500), z_to_luminosity_distance=dict(create_new=False, resolution=500), differential_comoving_volume=dict(create_new=False, resolution=500))











   .. rubric:: Examples

   >>> from ler.gw_source_population import CBCSourceParameterDistribution
   >>> cbc = CBCSourceParameterDistribution()
   >>> params = cbc.sample_gw_parameters(size=1000)
   >>> print("sampled parameters=",list(params.keys()))

   Instance Attributes
   ----------
   CompactBinaryPopulation has the following instance attributes:

   +-------------------------------------+----------------------------------+
   | Atrributes                          | Type                             |
   +=====================================+==================================+
   |:attr:`~z_min`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~z_max`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~event_type`                  | `str`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~source_priors`               | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~source_priors_params`        | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~cosmo`                       | `astropy.cosmology`              |
   +-------------------------------------+----------------------------------+
   |:attr:`~spin_zero`                   | `bool`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~spin_precession`             | `bool`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~directory`                   | `str`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~create_new_interpolator`     | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~available_gw_prior_list_and_its_params`                            |
   +-------------------------------------+----------------------------------+
   |                                     | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~gw_param_samplers`           | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~gw_param_samplers_params`    | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~sampler_names`               | `dict`                           |
   +-------------------------------------+----------------------------------+

   Instance Methods
   ----------
   CompactBinaryPopulation has the following instance methods:

   +-------------------------------------+----------------------------------+
   | Methods                             | Type                             |
   +=====================================+==================================+
   |:meth:`~source_priors_categorization`                                   |
   +-------------------------------------+----------------------------------+
   |                                     | Function to categorize the event |
   |                                     | priors and its parameters        |
   +-------------------------------------+----------------------------------+
   |:meth:`~lookup_table_luminosity_distance`                               |
   |                                     | Function to create a lookup      |
   |                                     | table for converting redshift    |
   |                                     | to luminosity distance           |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_gw_parameters`        | Function to sample all the       |
   |                                     | intrinsic and extrinsic          |
   |                                     | parameters of compact binaries   |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_source_frame_masses`  | Function to sample source mass1  |
   |                                     | and mass2                        |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_geocent_time`         | Function to sample geocent time  |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_zs`                   | Function to sample source        |
   |                                     | redshift                         |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_ra`                   | Function to sample right         |
   |                                     | ascension (sky position)         |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_dec`                  | Function to sample declination   |
   |                                     | (sky position)                   |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_phase`                | Function to sample coalescence   |
   |                                     | phase                            |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_psi`                  | Function to sample polarization  |
   |                                     | angle                            |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_theta_jn`             | Function to sample inclination   |
   |                                     | angle                            |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_a1`                   | Function to sample spin1         |
   |                                     | magnitude                        |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_a2`                   | Function to sample spin2         |
   |                                     | magnitude                        |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_tilt_1`               | Function to sample tilt1 angle   |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_tilt_2`               | Function to sample tilt2 angle   |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_phi_12`               | Function to sample phi12 angle   |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_phi_jl`               | Function to sample phi_jl angle  |
   +-------------------------------------+----------------------------------+
   |:meth:`~binary_masses_BBH_popI_II_powerlaw_gaussian`                    |
   +-------------------------------------+----------------------------------+
   |                                     | Function to sample source mass1  |
   |                                     | and mass2 with PowerLaw+PEAK     |
   |                                     | model                            |
   +-------------------------------------+----------------------------------+
   |:meth:`~binary_masses_BBH_popIII_lognormal`                             |
   +-------------------------------------+----------------------------------+
   |                                     | Function to sample source mass1  |
   |                                     | and mass2 with popIII orgin from |
   |                                     | lognormal distribution. Refer to |
   |                                     | Ng et al. 2022. Eqn. 1 and 4     |
   +-------------------------------------+----------------------------------+
   |:meth:`~binary_masses_BBH_primordial_lognormal`                         |
   +-------------------------------------+----------------------------------+
   |                                     | Function to sample source mass1  |
   |                                     | and mass2 with primordial orgin  |
   |                                     | from lognormal distribution.     |
   |                                     | Refer to Ng et al. 2022. Eqn. 1  |
   |                                     | and 4                            |
   +-------------------------------------+----------------------------------+
   |:meth:`~binary_masses_BNS_gwcosmo`                                      |
   +-------------------------------------+----------------------------------+
   |                                     | Function to sample source mass1  |
   |                                     | and mass2 from powerlaw          |
   |                                     | distribution.                    |
   +-------------------------------------+----------------------------------+
   |:meth:`~binary_masses_BNS_bimodal`   | Function to sample source mass1  |
   |                                     | and mass2 from bimodal           |
   |                                     | distribution. Refer to           |
   |                                     | Will M. Farr et al. 2020 Eqn. 6  |
   +-------------------------------------+----------------------------------+
   |:meth:`~constant_values_n_size`      | Function to return array of      |
   |                                     | constant values of size n        |
   +-------------------------------------+----------------------------------+
   |:meth:`~sampler_uniform`             | Function to sample from uniform  |
   |                                     | distribution                     |
   +-------------------------------------+----------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:property:: available_gw_prior_list_and_its_params

      
      Dictionary with list all the available priors and it's corresponding parameters. This is an immutable instance attribute.













      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CompactBinaryPopulation()
      >>> priors = cbc.available_gw_prior_list_and_its_params
      >>> priors.keys()  # type of priors
      dict_keys(['merger_rate_density', 'source_frame_masses', 'spin', 'geocent_time', 'ra', 'phase', 'psi', 'theta_jn'])
      >>> priors['source_frame_masses'].keys()  # type of source_frame_masses priors
      dict_keys(['binary_masses_BBH_popI_II_powerlaw_gaussian', 'binary_masses_BBH_popIII_lognormal', 'binary_masses_BBH_primordial_lognormal', 'binary_masses_BNS_gwcosmo', 'binary_masses_BNS_bimodal'])
      >>> priors['source_frame_masses']['binary_masses_BBH_popI_II_powerlaw_gaussian'].keys()  # parameters of binary_masses_BBH_popI_II_powerlaw_gaussian
      dict_keys(['mminbh', 'mmaxbh', 'alpha', 'mu_g', 'sigma_g', 'lambda_peak', 'delta_m', 'beta'])



      ..
          !! processed by numpydoc !!

   .. py:property:: sample_zs

      
      Function to sample redshifts with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **zs** : `numpy.ndarray` (1D array of floats)
              Array of redshifts













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_source_frame_masses

      
      Function to sample source frame masses (mass1_source, mass2_source) with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **mass_1_source** : `numpy.ndarray` (1D array of floats)
              Array of mass1 in source frame

          **mass_2_source** : `numpy.ndarray` (1D array of floats)
              Array of mass2 in source frame













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_geocent_time

      
      Function to sample geocent time with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **geocent_time** : `numpy.ndarray` (1D array of floats)
              Array of geocent_time or time of coalescence













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_ra

      
      Function to sample right ascension of sky position with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **ra** : `numpy.ndarray` (1D array of floats)
              Array of right ascension of sky position













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_dec

      
      Function to sample declination of sky position with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **dec** : `numpy.ndarray` (1D array of floats)
              Array of declination of sky position













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_phase

      
      Function to sample coalescence phase with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **phase** : `numpy.ndarray` (1D array of floats)
              Array of coalescence phase













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_psi

      
      Function to sample polarization angle with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **psi** : `numpy.ndarray` (1D array of floats)
              Array of polarization angle













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_theta_jn

      
      Function to sample theta_jn with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **theta_jn** : `numpy.ndarray` (1D array of floats)
              Array of theta_jn













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_a_1

      
      Function to sample spin magnitude of the compact binaries (body1) with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **a_1** : `numpy.ndarray` (1D array of floats)
              Array of spin magnitude of the compact binaries (body1)













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_a_2

      
      Function to sample spin magnitude of the compact binaries (body2) with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **a_2** : `numpy.ndarray` (1D array of floats)
              Array of spin magnitude of the compact binaries (body2)













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_tilt_1

      
      Function to sample tilt angle of the compact binaries (body1) with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **tilt_1** : `numpy.ndarray` (1D array of floats)
              Array of tilt angle of the compact binaries (body1)













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_tilt_2

      
      Function to sample tilt angle of the compact binaries (body2) with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **tilt_2** : `numpy.ndarray` (1D array of floats)
              Array of tilt angle of the compact binaries (body2)













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_phi_12

      
      Function to sample azimuthal angle between the two spins with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **phi_12** : `numpy.ndarray` (1D array of floats)
              Array of azimuthal angle between the two spins













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_phi_jl

      
      Function to sample azimuthal angle between the total angular momentum and the orbital angular momentum with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **phi_jl** : `numpy.ndarray` (1D array of floats)
              Array of azimuthal angle between the total angular momentum and the orbital angular momentum













      ..
          !! processed by numpydoc !!

   .. py:attribute:: z_min

      
      ``float``

      Minimum redshift of the source population















      ..
          !! processed by numpydoc !!

   .. py:attribute:: z_max

      
      ``float``

      Maximum redshift of the source population















      ..
          !! processed by numpydoc !!

   .. py:attribute:: event_type

      
      ``str``

      Type of event to generate.

      e.g. 'BBH', 'BNS', 'NSBH'















      ..
          !! processed by numpydoc !!

   .. py:attribute:: source_priors

      
      ``dict``

      Dictionary of prior sampler functions.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: source_priors_params

      
      ``dict``

      Dictionary of prior sampler functions' input parameters.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: cosmo

      
      ``astropy.cosmology``

      Cosmology to use.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: spin_zero

      
      ``bool``

      If True, spin prior is set to zero.















      ..
          !! processed by numpydoc !!

   .. py:method:: lookup_table_luminosity_distance(z_min, z_max, directory)

      
      Function to create a lookup table for the differential comoving volume
      and luminosity distance wrt redshift.


      :Parameters:

          **z_min** : `float`
              Minimum redshift of the source population

          **z_max** : `float`
              Maximum redshift of the source population












      :Attributes:

          **z_to_luminosity_distance** : `scipy.interpolate.interpolate`
              Function to convert redshift to luminosity distance

          **differential_comoving_volume** : `scipy.interpolate.interpolate`
              Function to calculate the differential comoving volume


      ..
          !! processed by numpydoc !!

   .. py:method:: sample_gw_parameters(size=1000, param=None)

      
      Function to sample BBH/BNS/NSBH intrinsic and extrinsics parameters.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **gw_parameters** : `dict`
              Dictionary of sampled parameters
              gw_parameters.keys() = ['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'theta_jn', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl']










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> params = cbc.sample_gw_parameters(size=1000)
      >>> print("sampled parameters=",list(params.keys()))



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BBH_popI_II_powerlaw_gaussian(size, mminbh=4.98, mmaxbh=112.5, alpha=3.78, mu_g=32.27, sigma_g=3.88, lambda_peak=0.03, delta_m=4.8, beta=0.81, get_attribute=False, param=None)

      
      Function to sample source mass1 and mass2 with PowerLaw+PEAK model


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **mminbh** : `float`
              Minimum mass of the black hole (Msun)
              default: 4.98

          **mmaxbh** : `float`
              Maximum mass of the black hole (Msun)
              default: 86.22

          **alpha** : `float`
              Spectral index for the powerlaw of the primary mass distribution
              default: 2.63

          **mu_g** : `float`
              Mean of the Gaussian component in the primary mass distribution
              default: 33.07

          **sigma_g** : `float`
              Width of the Gaussian component in the primary mass distribution
              default: 5.69

          **lambda_peak** : `float`
              Fraction of the model in the Gaussian component
              default: 0.10

          **delta_m** : `float`
              Range of mass tapering on the lower end of the mass distribution
              default: 4.82

          **beta** : `float`
              Spectral index for the powerlaw of the mass ratio distribution

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(mminbh=4.98, mmaxbh=86.22, alpha=2.63, mu_g=33.07, sigma_g=5.69, lambda_peak=0.10, delta_m=4.82, beta=1.26)

      :Returns:

          **mass_1_source** : `numpy.ndarray` (1D array of floats)
              Array of mass1 in source frame (Msun)

          **mass_2_source** : `numpy.ndarray` (1D array of floats)
              Array of mass2 in source frame (Msun)










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> m1_src, m2_src = cbc.binary_masses_BBH_popI_II_powerlaw_gaussian(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BBH_popIII_lognormal(size, m_min=5.0, m_max=150.0, Mc=30.0, sigma=0.3, chunk_size=10000, get_attribute=False, param=None)

      
      Function to sample source mass1 and mass2 with pop III origin. Refer to Eqn. 1 and 4 of Ng et al. 2022


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **m_min** : `float`
              Minimum mass of the black hole (popIII) (Msun)
              default: 10.

          **m_max** : `float`
              Maximum mass of the black hole (popIII) (Msun)
              default: 100.

          **Mc** : `float`
              Mass scale; the distribution is centered around Mc
              default: 30.0

          **sigma** : `float`
              Width of the distribution
              default: 0.3

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(m_min=10., m_max=100., Mc=30.0, sigma=0.3)

      :Returns:

          **mass_1_source** : `numpy.ndarray` (1D array of floats)
              Array of mass1 in source frame (Msun)

          **mass_2_source** : `numpy.ndarray` (1D array of floats)
              Array of mass2 in source frame (Msun)










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> m1_src, m2_src = cbc.binary_masses_BBH_popIII_lognormal(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BBH_primordial_lognormal(size, m_min=1.0, m_max=100.0, Mc=20.0, sigma=0.3, chunk_size=10000, get_attribute=False, param=None)

      
      Function to sample source mass1 and mass2 with primordial origin. Refer to Eqn. 1 and 4 of Ng et al. 2022


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **m_min** : `float`
              Minimum mass of the black hole (primordial) (Msun)
              default: 10.

          **m_max** : `float`
              Maximum mass of the black hole (primordial) (Msun)
              default: 100.

          **Mc, sigma** : `float`
              Fitting parameters
              default: Mc=30.0, sigma=0.3

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(m_min=10., m_max=100., Mc=30.0, sigma=0.3)

      :Returns:

          **mass_1_source** : `numpy.ndarray` (1D array of floats)
              Array of mass1 in source frame (Msun)

          **mass_2_source** : `numpy.ndarray` (1D array of floats)
              Array of mass2 in source frame (Msun)













      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BNS_gwcosmo(size, mminns=1.0, mmaxns=3.0, alphans=0.0, get_attribute=False, param=None)

      
      Function to calculate source mass1 and mass2 of BNS from powerlaw distribution (gwcosmo)


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **mminns** : `float`
              Minimum mass of the BNS (Msun)
              default: 1.0

          **mmaxns** : `float`
              Maximum mass of the BNS (Msun)
              default: 3.0

          **alphans** : `float`
              Power law index
              default: 0.0

      :Returns:

          **mass_1_source** : `numpy.ndarray` (1D array of floats)
              Array of mass1 in source frame (Msun)

          **mass_2_source** : `numpy.ndarray` (1D array of floats)
              Array of mass2 in source frame (Msun)










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> m1_src, m2_src = cbc.binary_masses_BNS_gwcosmo(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_NSBH_broken_powerlaw(size, mminbh=26, mmaxbh=125, alpha_1=6.75, alpha_2=6.75, b=0.5, delta_m=5, mminns=1.0, mmaxns=3.0, alphans=0.0, get_attribute=False, param=None)

      
      Function to calculate source mass1 and mass2 of NSBH from powerlaw distribution (gwcosmo). Parameters are mminbh=26,mmaxbh=125,alpha_1=6.75,alpha_2=6.75,b=0.5,delta_m=5,mminns=1.0,mmaxns=3.0,alphans=0.0.


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **mminbh** : `float`
              Minimum mass of the black hole (Msun)
              default: 26

          **mmaxbh** : `float`
              Maximum mass of the black hole (Msun)
              default: 125

          **alpha_1** : `float`
              Power law index for the primary mass distribution
              default: 6.75

          **alpha_2** : `float`
              Power law index for the secondary mass distribution
              default: 6.75

          **b** : `float`
              Break point of the power law
              default: 0.5

          **delta_m** : `float`
              Range of mass tapering on
              default: 5

          **mminns** : `float`
              Minimum mass of the neutron star (Msun)
              default: 1.0

          **mmaxns** : `float`
              Maximum mass of the neutron star (Msun)
              default: 3.0

          **alphans** : `float`
              Power law index for the neutron star mass distribution
              default: 0.0

          **get_attribute** : `bool`
              If True, return a sampler function with size as the only input where parameters are fixed to the given values.

          **param** : `dict`
              Allows to pass in above parameters as dict.

      :Returns:

          **mass_1_source** : `numpy.ndarray` (1D array of floats)
              Array of mass1 in source frame (Msun)

          **mass_2_source** : `numpy.ndarray` (1D array of floats)
              Array of mass2 in source frame (Msun)










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> m1_src, m2_src = cbc.binary_masses_NSBH_broken_powerlaw(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_uniform(size, m_min=1.0, m_max=3.0, get_attribute=False, param=None)

      
      Function to sample source mass1 and mass2 from uniform distribution.


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **m_min** : `float`
              Minimum mass of the BNS
              default: 1.0

          **m_max** : `float`
              Maximum mass of the BNS
              default: 3.0

          **get_attribute** : `bool`
              If True, return a sampler function with size as the only input where parameters are fixed to the given values.

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(m_min=1.0, m_max=3.0)

      :Returns:

          **mass_1_source** : `numpy.ndarray` (1D array of floats)
              Array of mass1 in source frame (Msun)

          **mass_2_source** : `numpy.ndarray` (1D array of floats)
              Array of mass2 in source frame (Msun)










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> m1_src, m2_src = cbc.binary_masses_uniform(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BNS_bimodal(size, w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1.0, mmax=2.3, resolution=500, create_new=False, get_attribute=False, param=None)

      
      Function to sample source mass1 and mass2 from bimodal distribution. Refer to Will M. Farr et al. 2020 Eqn. 6, https://arxiv.org/pdf/2005.00032.pdf .


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **w** : `float`
              Weight of the left peak
              default: 0.643

          **muL** : `float`
              Mean of the left peak
              default: 1.352

          **sigmaL** : `float`
              Width of the left peak
              default: 0.08

          **muR** : `float`
              Mean of the right peak
              default: 1.88

          **sigmaR** : `float`
              Width of the right peak
              default: 0.3

          **mmin** : `float`
              Minimum mass of the BNS
              default: 1.0

          **mmax** : `float`
              Maximum mass of the BNS
              default: 2.3

          **resolution** : `int`
              Number of points to sample
              default: 500

          **create_new** : `bool`
              If True, create new interpolator
              default: False

          **get_attribute** : `bool`
              If True, return a sampler function with size as the only input where parameters are fixed to the given values.

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1.0, mmax=2.3, resolution=500)

      :Returns:

          **mass_1_source** : `numpy.ndarray` (1D array of floats)
              Array of mass1 in source frame (Msun)

          **mass_2_source** : `numpy.ndarray` (1D array of floats)
              Array of mass2 in source frame (Msun)










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> m1_src, m2_src = cbc.binary_masses_BNS_bimodal(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: constant_values_n_size(size=100, value=0.0, get_attribute=False, param=None)

      
      Function to sample constant values of size n.


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **value** : `float`
              Constant value
              default: 0.0

          **get_attribute** : `bool`
              If True, return the njitted sampler function with size as the only input where parameters are fixed to the given values.

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(value=0.0)

      :Returns:

          **values** : `numpy.ndarray` (1D array of floats)
              Array of constant values










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> value = cbc.constant_values_n_size(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: sampler_uniform(size, min_=0, max_=np.pi, get_attribute=False, param=None)

      
      Function to sample values from uniform distribution.


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **start_time** : `float`
              Start time of the uniform distribution
              default: 1238166018

          **end_time** : `float`
              End time of the uniform distribution
              default: 1238166018 + 31536000

          **get_attribute** : `bool`
              If True, return the njitted sampler function with size as the only input where parameters are fixed to the given values.

          **param** : `dict`
              Allows to pass in above parameters as dict.

      :Returns:

          **values** : `numpy.ndarray` (1D array of floats)
              Array of uniformly distributed values in the range of [min_, max_]










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> value = cbc.sampler_uniform(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: sampler_cosine(size, get_attribute=False, param=None)

      
      Function to sample from sine distribution at the limit of [-np.pi/2, np.pi/2]


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **get_attribute** : `bool`
              If True, return the njitted sampler function with size as the only input where parameters are fixed to the given values.

          **param** : None
              This parameter is not used. It is only here to make the function signature consistent with other samplers.

      :Returns:

          **sine** : `numpy.ndarray` (1D array of floats)
              Array of values in the range of [-np.pi/2, np.pi/2]













      ..
          !! processed by numpydoc !!

   .. py:method:: sampler_sine(size, get_attribute=False, param=None)

      
      Function to sample from sine distribution at the limit of [0, np.pi]


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **get_attribute** : `bool`
              If True, return the njitted sampler function with size as the only input where parameters are fixed to the given values.

          **param** : None
              This parameter is not used. It is only here to make the function signature consistent with other samplers.

      :Returns:

          **sine** : `numpy.ndarray` (1D array of floats)
              Array of values in the range of [0, np.pi]













      ..
          !! processed by numpydoc !!

   .. py:method:: source_priors_categorization(event_type, source_priors, event_prior_params)

      
      Function to categorize the event priors and its parameters.


      :Parameters:

          **event_type** : `str`
              Type of event to generate.
              e.g. 'BBH', 'BNS', 'BBH_popIII', 'BBH_primordial', 'NSBH'

          **source_priors** : `dict`
              Dictionary of prior sampler functions for each parameter

          **event_prior_params** : `dict`
              Dictionary of sampler parameters for each GW parameter

      :Returns:

          **source_priors_** : `dict`
              Dictionary of prior sampler functions for each parameter

          **event_prior_params_** : `dict`
              Dictionary of sampler parameters for each parameter

          **sampler_names_** : `dict`
              Dictionary of sampler names with description










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> source_priors, event_prior_params, sampler_names = cbc.source_priors_categorization(event_type='BBH', source_priors=None, event_prior_params=None)
      >>> print(source_priors.keys())
      >>> print(event_prior_params.keys())
      >>> print(sampler_names.keys())



      ..
          !! processed by numpydoc !!


.. py:function:: load_json(file_name)

   
   Load a json file.


   :Parameters:

       **file_name** : `str`
           json file name for storing the parameters.

   :Returns:

       **param** : `dict`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: append_json(file_name, new_dictionary, old_dictionary=None, replace=False)

   
   Append (values with corresponding keys) and update a json file with a dictionary. There are four options:

   1. If old_dictionary is provided, the values of the new dictionary will be appended to the old dictionary and save in the 'file_name' json file.
   2. If replace is True, replace the json file (with the 'file_name') content with the new_dictionary.
   3. If the file does not exist, create a new one with the new_dictionary.
   4. If none of the above, append the new dictionary to the content of the json file.

   :Parameters:

       **file_name** : `str`
           json file name for storing the parameters.

       **new_dictionary** : `dict`
           dictionary to be appended to the json file.

       **old_dictionary** : `dict`, optional
           If provided the values of the new dictionary will be appended to the old dictionary and save in the 'file_name' json file.
           Default is None.

       **replace** : `bool`, optional
           If True, replace the json file with the dictionary. Default is False.














   ..
       !! processed by numpydoc !!

.. py:function:: get_param_from_json(json_file)

   
   Function to get the parameters from json file.


   :Parameters:

       **json_file** : `str`
           json file name for storing the parameters.

   :Returns:

       **param** : `dict`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: batch_handler(size, batch_size, sampling_routine, output_jsonfile, save_batch=True, resume=False)

   
   Function to run the sampling in batches.


   :Parameters:

       **size** : `int`
           number of samples.

       **batch_size** : `int`
           batch size.

       **sampling_routine** : `function`
           function to sample the parameters.
           e.g. unlensed_sampling_routine() or lensed_sampling_routine()

       **output_jsonfile** : `str`
           name of the json file to store the parameters.

       **resume** : `bool`
           if True, it will resume the sampling from the last batch.
           default resume = False.














   ..
       !! processed by numpydoc !!

.. py:class:: GWRATES(npool=int(4), z_min=0.0, z_max=10.0, event_type='BBH', size=100000, batch_size=50000, cosmology=None, snr_finder=None, pdet_finder=None, list_of_detectors=None, json_file_names=None, interpolator_directory='./interpolator_pickle', ler_directory='./ler_data', verbose=True, **kwargs)


   Bases: :py:obj:`ler.gw_source_population.CBCSourceParameterDistribution`

   
   Class to calculate both the rates of lensed and gw events. Please note that parameters of the simulated events are stored in json file but not as an attribute of the class. This saves RAM memory.


   :Parameters:

       **npool** : `int`
           number of cores to use.
           default npool = 4.

       **z_min** : `float`
           minimum redshift.
           default z_min = 0.
           for popI_II, popIII, primordial, BNS z_min = 0., 5., 5., 0. respectively.

       **z_max** : `float`
           maximum redshift.
           default z_max = 10.
           for popI_II, popIII, primordial, BNS z_max = 10., 40., 40., 5. respectively.

       **event_type** : `str`
           type of event to generate.
           default event_type = 'BBH'. Other options are 'BNS', 'NSBH'.

       **size** : `int`
           number of samples for sampling.
           default size = 100000. To get stable rates, size should be large (>=1e6).

       **batch_size** : `int`
           batch size for SNR calculation.
           default batch_size = 50000.
           reduce the batch size if you are getting memory error.
           recommended batch_size = 200000, if size = 1000000.

       **cosmology** : `astropy.cosmology`
           cosmology to use for the calculation.
           default cosmology = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7).

       **snr_finder** : `str` or `function`
           default snr_finder = 'gwsnr'.
           if None, the SNR will be calculated using the gwsnr package.
           if custom snr finder function is provided, the SNR will be calculated using a custom function. The custom function should follow the following signature:
           def snr_finder(gw_param_dict):
               ...
               return optimal_snr_dict
           where optimal_snr_dict.keys = ['optimal_snr_net']. Refer to `gwsnr` package's GWSNR.snr attribute for more details.

       **pdet_finder** : `function`
           default pdet_finder = None.
           The rate calculation uses either the pdet_finder or the snr_finder to calculate the detectable events. The custom pdet finder function should follow the following signature:
           def pdet_finder(gw_param_dict):
               ...
               return pdet_net_dict
           where pdet_net_dict.keys = ['pdet_net']. For example uses, refer to [GRB pdet example](https://ler.readthedocs.io/en/latest/examples/rates/grb%20detection%20rate.html).

       **json_file_names: `dict`**
           names of the json files to strore the necessary parameters.
           default json_file_names = {'gwrates_param':'gwrates_params.json', 'gw_param': 'gw_param.json', 'gw_param_detectable': 'gw_param_detectable.json'}.

       **interpolator_directory** : `str`
           directory to store the interpolators.
           default interpolator_directory = './interpolator_pickle'. This is used for storing the various interpolators related to `ler` and `gwsnr` package.

       **ler_directory** : `str`
           directory to store the parameters.
           default ler_directory = './ler_data'. This is used for storing the parameters of the simulated events.

       **verbose** : `bool`
           default verbose = True.
           if True, the function will print all chosen parameters.
           Choose False to prevent anything from printing.

       **kwargs** : `keyword arguments`
           Note : kwargs takes input for initializing the :class:`~ler.gw_source_population.CBCSourceParameterDistribution` and :class:`~ler.gw_source_population.CBCSourceRedshiftDistribution` classes. If snr_finder='gwsnr', then kwargs also takes input for initializing the :class:`~gwsnr.GWSNR` class. Please refer to the respective classes for more details.











   .. rubric:: Examples

   >>> from ler.rates import GWRATES
   >>> ler = GWRATES()
   >>> ler.gw_cbc_statistics();
   >>> ler.gw_rate();

   Instance Attributes
   ----------
   LeR class has the following attributes,

   +-------------------------------------+----------------------------------+
   | Atrributes                          | Type                             |
   +=====================================+==================================+
   |:attr:`~npool`                       | `int`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~z_min`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~z_max`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~event_type`                  | `str`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~cosmo`                       | `astropy.cosmology`              |
   +-------------------------------------+----------------------------------+
   |:attr:`~size`                        | `int`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~batch_size`                  | `int`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~json_file_names`             | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~interpolator_directory`      | `str`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~ler_directory`               | `str`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~gwsnr`                       | `bool`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~gw_param_sampler_dict`       | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~snr_calculator_dict`         | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~gw_param`                    | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~gw_param_detectable`         | `dict`                           |
   +-------------------------------------+----------------------------------+

   Instance Methods
   ----------
   LeR class has the following methods,

   +-------------------------------------+----------------------------------+
   | Methods                             | Description                      |
   +=====================================+==================================+
   |:meth:`~class_initialization`        | Function to initialize the       |
   |                                     | parent classes                   |
   +-------------------------------------+----------------------------------+
   |:meth:`~gwsnr_intialization`         | Function to initialize the       |
   |                                     | gwsnr class                      |
   +-------------------------------------+----------------------------------+
   |:meth:`~snr`                         | Function to get the snr with the |
   |                                     | given parameters.                |
   +-------------------------------------+----------------------------------+
   |:meth:`~store_gwrates_params`        | Function to store the all the    |
   |                                     | necessary parameters.            |
   +-------------------------------------+----------------------------------+
   |:meth:`~gw_cbc_statistics`           | Function to generate gw          |
   |                                     | GW source parameters.            |
   +-------------------------------------+----------------------------------+
   |:meth:`~gw_sampling_routine`         | Function to generate gw          |
   |                                     | GW source parameters.            |
   +-------------------------------------+----------------------------------+
   |:meth:`~gw_rate`                     | Function to calculate the        |
   |                                     | gw rate.                         |
   +-------------------------------------+----------------------------------+
   |:meth:`~selecting_n_gw_detectable_events`                               |
   +-------------------------------------+----------------------------------+
   |                                     | Function to select n gw    |
   |                                     | detectable events.               |
   +-------------------------------------+----------------------------------+
   |:meth:`~gw_param_plot`               | Function to plot the             |
   |                                     | distribution of the GW source    |
   |                                     | parameters.                      |
   +-------------------------------------+----------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:property:: snr

      
      Function to get the snr with the given parameters.


      :Parameters:

          **gw_param_dict** : `dict`
              dictionary of GW source parameters.
              mass_1 : `numpy.ndarray` or `float`
                  mass_1 of the compact binary (detector frame) (Msun).
              mass_2 : `numpy.ndarray` or `float`
                  mass_2 of the compact binary (detector frame) (Msun).
              luminosity_distance : `numpy.ndarray` or `float`
                  luminosity distance of the source (Mpc).
              theta_jn : `numpy.ndarray` or `float`
                  inclination angle of the source (rad).
              psi : `numpy.ndarray` or `float`
                  polarization angle of the source (rad).
              phase : `numpy.ndarray` or `float`
                  phase of GW at reference frequency  (rad).
              geocent_time : `numpy.ndarray` or `float`
                  GPS time of coalescence (s).
              ra : `numpy.ndarray` or `float`
                  right ascension of the source (rad).
              dec : `numpy.ndarray` or `float`
                  declination of the source (rad).
              a_1 : `numpy.ndarray` or `float`
                  dimensionless spin magnitude of the more massive object.
              a_2 : `numpy.ndarray` or `float`
                  dimensionless spin magnitude of the less massive object.
              tilt_1 : `numpy.ndarray` or `float`
                  tilt angle of the more massive object spin.
              tilt_2 : `numpy.ndarray` or `float`
                  tilt angle of the less massive object spin.
              phi_12 : `numpy.ndarray` or `float`
                  azimuthal angle between the two spin vectors.
              phi_jl : `numpy.ndarray` or `float`
                  azimuthal angle between total angular momentum and the orbital angular momentum.

      :Returns:

          **optimal_snr_list** : `list`
              e.g. [optimal_snr_net, 'L1', 'H1', 'V1']
              optimal_snr_net : `numpy.ndarray` or `float`
                  optimal snr of the network.
              'H1' : `numpy.ndarray` or `float`
                  optimal snr of H1.
              'L1' : `numpy.ndarray` or `float`
                  optimal snr of L1.
              'V1' : `numpy.ndarray` or `float`
                  optimal snr of V1.













      ..
          !! processed by numpydoc !!

   .. py:property:: gw_param

      
      Function to get data from the json file self.json_file_names["gw_param"].



      :Returns:

          **gw_param** : `dict`
              dictionary of gw GW source parameters.













      ..
          !! processed by numpydoc !!

   .. py:property:: gw_param_detectable

      
      Function to get data from the json file self.json_file_names["gw_param_detectable"].



      :Returns:

          **gw_param_detectable** : `dict`
              dictionary of gw GW source parameters.













      ..
          !! processed by numpydoc !!

   .. py:attribute:: z_min

      
      ``float``

      Minimum redshift of the source population















      ..
          !! processed by numpydoc !!

   .. py:attribute:: z_max

      
      ``float``

      Maximum redshift of the source population















      ..
          !! processed by numpydoc !!

   .. py:attribute:: event_type

      
      ``str``

      Type of event to generate.

      e.g. 'BBH', 'BNS', 'NSBH'















      ..
          !! processed by numpydoc !!

   .. py:attribute:: cosmo

      
      ``astropy.cosmology``

      Cosmology to use for the calculation.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: size

      
      ``int``

      Number of samples for sampling.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: batch_size

      
      ``int``

      Batch size for sampling.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: json_file_names

      
      ``dict``

      Names of the json files to store the necessary parameters.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: interpolator_directory

      
      ``str``

      Directory to store the interpolators.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: ler_directory

      
      ``str``

      Directory to store the parameters.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: gwsnr

      
      ``bool``

      If True, the SNR will be calculated using the gwsnr package.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: gw_param_sampler_dict

      
      ``dict``

      Dictionary of parameters to initialize the ``CBCSourceParameterDistribution`` class.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: snr_calculator_dict

      
      ``dict``

      Dictionary of parameters to initialize the ``GWSNR`` class.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: gw_param

      
      ``dict``

      Dictionary of GW source parameters. The included parameters and their units are as follows (for default settings):

      +--------------------+--------------+--------------------------------------+
      | Parameter          | Units        | Description                          |
      +====================+==============+======================================+
      | zs                 |              | redshift of the source               |
      +--------------------+--------------+--------------------------------------+
      | geocent_time       | s            | GPS time of coalescence              |
      +--------------------+--------------+--------------------------------------+
      | ra                 | rad          | right ascension                      |
      +--------------------+--------------+--------------------------------------+
      | dec                | rad          | declination                          |
      +--------------------+--------------+--------------------------------------+
      | phase              | rad          | phase of GW at reference frequency   |
      +--------------------+--------------+--------------------------------------+
      | psi                | rad          | polarization angle                   |
      +--------------------+--------------+--------------------------------------+
      | theta_jn           | rad          | inclination angle                    |
      +--------------------+--------------+--------------------------------------+
      | luminosity_distance| Mpc          | luminosity distance                  |
      +--------------------+--------------+--------------------------------------+
      | mass_1_source      | Msun         | mass_1 of the compact binary         |
      |                    |              | (source frame)                       |
      +--------------------+--------------+--------------------------------------+
      | mass_2_source      | Msun         | mass_2 of the compact binary         |
      |                    |              | (source frame)                       |
      +--------------------+--------------+--------------------------------------+
      | mass_1             | Msun         | mass_1 of the compact binary         |
      |                    |              | (detector frame)                     |
      +--------------------+--------------+--------------------------------------+
      | mass_2             | Msun         | mass_2 of the compact binary         |
      |                    |              | (detector frame)                     |
      +--------------------+--------------+--------------------------------------+
      | L1                 |              | optimal snr of L1                    |
      +--------------------+--------------+--------------------------------------+
      | H1                 |              | optimal snr of H1                    |
      +--------------------+--------------+--------------------------------------+
      | V1                 |              | optimal snr of V1                    |
      +--------------------+--------------+--------------------------------------+
      | optimal_snr_net    |              | optimal snr of the network           |
      +--------------------+--------------+--------------------------------------+















      ..
          !! processed by numpydoc !!

   .. py:attribute:: gw_param_detectable

      
      ``dict``

      Dictionary of detectable GW source parameters. It includes the same parameters as the :attr:`~gw_param` attribute.















      ..
          !! processed by numpydoc !!

   .. py:method:: print_all_params()

      
      Function to print all the parameters.
















      ..
          !! processed by numpydoc !!

   .. py:method:: class_initialization(params=None)

      
      Function to initialize the parent classes.


      :Parameters:

          **params** : `dict`
              dictionary of parameters to initialize the parent classes














      ..
          !! processed by numpydoc !!

   .. py:method:: gwsnr_intialization(params=None)

      
      Function to initialize the GWSNR class from the `gwsnr` package.


      :Parameters:

          **params** : `dict`
              dictionary of parameters to initialize the gwsnr class














      ..
          !! processed by numpydoc !!

   .. py:method:: store_gwrates_params(output_jsonfile='gwrates_params.json')

      
      Function to store the all the necessary parameters. This is useful for reproducing the results. All the parameters stored are in string format to make it json compatible.


      :Parameters:

          **output_jsonfile** : `str`
              name of the json file to store the parameters














      ..
          !! processed by numpydoc !!

   .. py:method:: gw_cbc_statistics(size=None, resume=False, save_batch=False, output_jsonfile=None)

      
      Function to generate gw GW source parameters. This function calls the gw_sampling_routine function to generate the parameters in batches. The generated parameters are stored in a json file; and if save_batch=True, it keeps updating the file in batches.


      :Parameters:

          **size** : `int`
              number of samples.
              default size = 100000.

          **resume** : `bool`
              resume = False (default) or True.
              if True, the function will resume from the last batch.

          **save_batch** : `bool`
              if True, the function will save the parameters in batches. if False, the function will save all the parameters at the end of sampling. save_batch=False is faster.

          **output_jsonfile** : `str`
              json file name for storing the parameters.
              default output_jsonfile = 'gw_params.json'. Note that this file will be stored in the self.ler_directory.

      :Returns:

          **gw_param** : `dict`
              dictionary of gw GW source parameters. Refer to :attr:`~gw_param` for details.










      .. rubric:: Examples

      >>> from ler.rates import GWRATES
      >>> ler = GWRATES()
      >>> param = ler.gw_cbc_statistics()



      ..
          !! processed by numpydoc !!

   .. py:method:: gw_sampling_routine(size, output_jsonfile, resume=False, save_batch=True)

      
      Function to generate GW source parameters. This function also stores the parameters in json file in the current batch if save_batch=True.


      :Parameters:

          **size** : `int`
              number of samples.
              default size = 100000.

          **output_jsonfile** : `str`
              json file name for storing the parameters.
              default output_jsonfile = 'gw_params.json'. Note that this file will be stored in the self.ler_directory.

          **resume** : `bool`
              resume = False (default) or True.
              if True, it appends the new samples to the existing json file.

          **save_batch** : `bool`
              if True, the function will save the parameters in batches. if False, the function will save all the parameters at the end of sampling. save_batch=False is faster.

      :Returns:

          **gw_param** : `dict`
              dictionary of gw GW source parameters. Refer to :attr:`~gw_param` for details.













      ..
          !! processed by numpydoc !!

   .. py:method:: gw_rate(gw_param=None, snr_threshold=8.0, pdet_threshold=0.5, output_jsonfile=None, detectability_condition='step_function', snr_recalculation=False, snr_threshold_recalculation=5.5)

      
      Function to calculate the GW rate. This function also stores the parameters of the detectable events in json file. There are two conditions for detectability: 'step_function' and 'pdet'.

      1. 'step_function': If two images have SNR>8.0, then the event is detectable. This is a step function. This is with the assumption that SNR function is provided and not None.
      2. 'pdet':
          i) If self.pdet is None and self.snr is not None, then it will calculate the pdet from the snr. There is no hard cut for this pdet and can have value ranging from 0 to 1 near the threshold.
          ii) If self.pdet is not None, then it will use the generated pdet.

      :Parameters:

          **gw_param** : `dict` or `str`
              dictionary of GW source parameters or json file name.
              default gw_param = self.json_file_names["gw_param"]

          **snr_threshold** : `float`
              threshold for detection signal to noise ratio.
              e.g. snr_threshold = 8.

          **pdet_threshold** : `float`
              threshold for detection probability.
              e.g. pdet_threshold = 0.5.

          **output_jsonfile** : `str`
              json file name for storing the parameters of the detectable events.
              default output_jsonfile = 'gw_params_detectable.json'.

          **detectability_condition** : `str`
              detectability condition.
              default detectability_condition = 'step_function'.
              other options are 'pdet'.

          **snr_recalculation** : `bool`
              if True, the SNR of centain events (snr>snr_threshold_recalculation)will be recalculate with 'inner-product' method. This is useful when the snr is calculated with 'ann' method.
              default snr_recalculation = False.

          **threshold_snr_recalculation** : `float`
              threshold for recalculation of detection signal to noise ratio.

      :Returns:

          **total_rate** : `float`
              total gw rate (Mpc^-3 yr^-1).

          **gw_param** : `dict`
              dictionary of gw GW source parameters of the detectable events. Refer to :attr:`~gw_param` for details.










      .. rubric:: Examples

      >>> from ler.rates import GWRATES
      >>> ler = GWRATES()
      >>> ler.gw_cbc_statistics();
      >>> total_rate, gw_param = ler.gw_rate()



      ..
          !! processed by numpydoc !!

   .. py:method:: rate_function(detectable_size, total_size, verbose=True)

      
      General helper function to calculate the rate for GW events.


      :Parameters:

          **detectable_size** : `int`
              number of detectable events.

          **total_size** : `int`
              total number of events.

          **param_type** : `str`
              type of parameters.

      :Returns:

          **rate** : `float`
              rate of the events.










      .. rubric:: Examples

      >>> from ler.rates import LeR
      >>> ler = LeR()
      >>> rate = ler.rate_function(detectable_size=100, total_size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: selecting_n_gw_detectable_events(size=100, batch_size=None, snr_threshold=8.0, pdet_threshold=0.5, resume=False, output_jsonfile='gw_params_n_detectable.json', meta_data_file='meta_gw.json', detectability_condition='step_function', trim_to_size=True, snr_recalculation=False, snr_threshold_recalculation=5.5)

      
      Function to generate n GW detectable events. This fuction samples the GW parameters and save only the detectable events in json file. It also records metadata in the JSON file, which includes the total number of events and the cumulative rate of events. This functionality is particularly useful for generating a fixed or large number of detectable events until the event rates stabilize.


      :Parameters:

          **size** : `int`
              number of samples to be selected.
              default size = 100.

          **batch_size** : `int`
              batch size for sampling.
              default batch_size = 50000.

          **snr_threshold** : `float`
              threshold for detection signal to noise ratio.
              e.g. snr_threshold = 8.

          **pdet_threshold** : `float`
              threshold for detection probability.
              default pdet_threshold = 0.5.

          **resume** : `bool`
              resume = False (default) or True.
              if True, the function will resume from the last batch.

          **output_jsonfile** : `str`
              json file name for storing the parameters of the detectable events.
              default output_jsonfile = 'n_gw_param_detectable.json'.

          **meta_data_file** : `str`
              json file name for storing the metadata.
              default meta_data_file = 'meta_gw.json'.

          **detectability_condition** : `str`
              detectability condition.
              default detectability_condition = 'step_function'.
              other options are 'pdet'.

          **trim_to_size** : `bool`
              if True, the final result will be trimmed to size.
              default trim_to_size = True.

          **snr_recalculation** : `bool`
              if True, the SNR of centain events (snr>snr_threshold_recalculation)will be recalculate with 'inner-product' method. This is useful when the snr is calculated with 'ann' method of `gwsnr`.
              default snr_recalculation = False.

          **snr_threshold_recalculation** : `float`
              threshold for recalculation of detection signal to noise ratio.
              default snr_threshold_recalculation = 5.5.

      :Returns:

          **param_final** : `dict`
              dictionary of gw GW source parameters of the detectable events. Refer to :attr:`~gw_param` for details.










      .. rubric:: Examples

      >>> from ler.rates import LeR
      >>> ler = LeR()
      >>> gw_param = ler.selecting_n_gw_detectable_events(size=100)



      ..
          !! processed by numpydoc !!


