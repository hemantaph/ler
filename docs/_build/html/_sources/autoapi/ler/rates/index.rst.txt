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
   ler.rates.remove_file
   ler.rates.load_json
   ler.rates.append_json
   ler.rates.get_param_from_json
   ler.rates.batch_handler
   ler.rates.remove_file



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

.. py:function:: batch_handler(size, batch_size, sampling_routine, output_jsonfile, save_batch=True, resume=False, param_name='parameters')

   
   Function to run the sampling in batches.


   :Parameters:

       **size** : `int`
           number of samples.

       **batch_size** : `int`
           batch size.

       **sampling_routine** : `function`
           sampling function. It should have 'size' as input and return a dictionary.

       **output_jsonfile** : `str`
           json file name for storing the parameters.

       **save_batch** : `bool`, optional
           if True, save sampled parameters in each iteration. Default is True.

       **resume** : `bool`, optional
           if True, resume sampling from the last batch. Default is False.

       **param_name** : `str`, optional
           name of the parameter. Default is 'parameters'.

   :Returns:

       **dict_buffer** : `dict`
           dictionary of parameters.













   ..
       !! processed by numpydoc !!

.. py:function:: remove_file(file_name)

   
   Remove a file.
















   ..
       !! processed by numpydoc !!

.. py:class:: LeR(npool=int(4), z_min=0.0, z_max=10.0, event_type='BBH', lens_type='epl_shear_galaxy', cosmology=None, pdet_finder=None, json_file_names=None, interpolator_directory='./interpolator_json', create_new_interpolator=False, ler_directory='./ler_data', verbose=True, **kwargs)


   Bases: :py:obj:`ler.lens_galaxy_population.LensGalaxyParameterDistribution`

   
   Class to sample lensed and unlensed GW events and calculate their detection rates.

   This class provides functionality for sampling gravitational wave source parameters,
   detection probabilities, and computing detection rates for both lensed and unlensed
   compact binary coalescence events.
   Parameters of simulated events are stored in JSON files (not as class attributes)
   to conserve RAM memory.

   Key Features:

   - Sampling of unlensed and lensed CBC event parameters

   - Detection probability calculation using ``gwsnr`` package or custom functions

   - Rate calculation for detectable events

   - Batch processing for memory efficiency

   - JSON-based parameter storage for reproducibility

   :Parameters:

       **npool** : ``int``
           Number of cores to use for parallel processing.

           default: 4

       **z_min** : ``float``
           Minimum redshift of the source population.

           default: 0.0

       **z_max** : ``float``
           Maximum redshift of the source population.

           default: 10.0

       **event_type** : ``str``
           Type of event to generate. source_priors and source_priors_params will be set accordingly.

           Options:

           - 'BBH': Binary Black Hole

           - 'BNS': Binary Neutron Star

           - 'NSBH': Neutron Star-Black Hole

           default: 'BBH'

       **lens_type** : ``str``
           Type of lens model to use. lens_functions, lens_functions_params, lens_param_samplers and lens_param_samplers_params will be set accordingly.

           Options:

           - 'epl_shear_galaxy': Exponential Power Law Shear Galaxy

           - 'sie_galaxy': Singular Isothermal Ellipsoid Galaxy

           - 'sis_galaxy': Singular Isothermal Sphere Galaxy

           default: 'epl_shear_galaxy'

       **cosmology** : ``astropy.cosmology``
           Cosmology to use for the calculation.

           default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)

       **pdet_finder** : ``function`` or ``None``
           Custom detection probability finder function.

           If None, uses gwsnr's pdet calculator.

           The function should follow the signature:

           ``def pdet_finder(gw_param_dict): return pdet_net_dict``

           where pdet_net_dict.keys = ['pdet_net'].

           default: None

       **json_file_names** : ``dict``
           Names of the JSON files to store the necessary parameters.

           default: dict(
               ler_params="ler_params.json",
               unlensed_param="unlensed_param.json",
               unlensed_param_detectable="unlensed_param_detectable.json",
               lensed_param="lensed_param.json",
               lensed_param_detectable="lensed_param_detectable.json"
           )

       **interpolator_directory** : ``str``
           Directory to store the interpolators.

           default: './interpolator_json'

       **create_new_interpolator** : ``bool`` or ``dict``
           Whether to create new interpolators. Look at :meth:`~ler.ler_rates.LER.create_new_interpolator` for details.

           Options:

           - True: Create all interpolators anew

           - False: Load existing interpolators if available

           - dict: Specify which interpolators to create new

           default: False

       **ler_directory** : ``str``
           Directory to store the output parameters.

           default: './ler_data'

       **verbose** : ``bool``
           If True, print all chosen parameters during initialization.

           default: True

       **\*\*kwargs** : ``dict``
           Additional keyword arguments passed to parent classes:

           :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution`,

           :class:`~ler.gw_source_population.CBCSourceParameterDistribution`,

           :class:`~ler.image_properties.ImageProperties`, and

           :class:`~gwsnr.GWSNR` (if snr_finder='gwsnr').









   .. rubric:: Notes

   - ``LeR`` class inherits from :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution`.

     Refer to that class for additional inherited attributes and methods.

   - Parameters are stored in JSON files for memory efficiency and reproducibility.

   - For stable rate estimates, use size >= 1e6 samples.


   .. rubric:: Examples

   Basic usage:

   >>> from ler import LeR
   >>> ler = LeR()
   >>> unlensed_params = ler.unlensed_cbc_statistics()
   >>> ler.unlensed_rate()
   >>> lensed_params = ler.lensed_cbc_statistics()
   >>> ler.lensed_rate()
   >>> ler.rate_ratio()

   Instance Methods
   ----------
   LeR class has the following methods:

   +-----------------------------------------------------+------------------------------------------------+
   | Method                                              | Description                                    |
   +=====================================================+================================================+
   | :meth:`~unlensed_cbc_statistics`                    | Generate unlensed GW source parameters         |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~unlensed_sampling_routine`                  | Generate unlensed parameters with batching     |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~unlensed_rate`                              | Calculate the unlensed detection rate          |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~lensed_cbc_statistics`                      | Generate lensed GW source parameters           |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~lensed_sampling_routine`                    | Generate lensed parameters with batching       |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~lensed_rate`                                | Calculate the lensed detection rate            |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~rate_function`                              | General helper for rate calculation            |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~rate_ratio`                                 | Calculate lensed/unlensed rate ratio           |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~rate_comparison_with_rate_calculation`      | Calculate and compare lensed/unlensed rates    |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~selecting_n_unlensed_detectable_events`     | Select n unlensed detectable events            |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~selecting_n_lensed_detectable_events`       | Select n lensed detectable events              |
   +-----------------------------------------------------+------------------------------------------------+

   Instance Attributes
   ----------
   LeR class has the following attributes:

   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | Attribute                                      | Type             | Unit  | Description                                    |
   +================================================+==================+=======+================================================+
   | :meth:`~npool`                                 | ``int``          |       | Number of parallel processing cores            |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~z_min`                                 | ``float``        |       | Minimum source redshift                        |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~z_max`                                 | ``float``        |       | Maximum source redshift                        |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~event_type`                            | ``str``          |       | Type of CBC event (BBH, BNS, NSBH)             |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~lens_type`                             | ``str``          |       | Type of lens galaxy model                      |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~cosmo`                                 | ``Cosmology``    |       | Astropy cosmology object                       |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~json_file_names`                       | ``dict``         |       | JSON file names for parameter storage          |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~interpolator_directory`                | ``str``          |       | Directory for interpolator files               |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~ler_directory`                         | ``str``          |       | Directory for output parameter files           |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~pdet_finder`                           | ``callable``     |       | Detection probability finder function          |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~ler_args`                              | ``dict``         |       | All LeR initialization arguments               |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~create_new_interpolator`               | ``dict``         |       | Interpolator creation settings                 |
   +------------------------------------------------+------------------+-------+------------------------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:property:: npool

      
      Number of parallel processing cores.



      :Returns:

          **npool** : ``int``
              Number of logical cores to use for multiprocessing.

              default: 4













      ..
          !! processed by numpydoc !!

   .. py:property:: z_min

      
      Minimum redshift of the source population.



      :Returns:

          **z_min** : ``float``
              Minimum source redshift for sampling.

              default: 0.0













      ..
          !! processed by numpydoc !!

   .. py:property:: z_max

      
      Maximum redshift of the source population.



      :Returns:

          **z_max** : ``float``
              Maximum source redshift for sampling.

              default: 10.0













      ..
          !! processed by numpydoc !!

   .. py:property:: event_type

      
      Type of compact binary coalescence event.



      :Returns:

          **event_type** : ``str``
              Type of CBC event.

              Options:

              - 'BBH': Binary Black Hole

              - 'BNS': Binary Neutron Star

              - 'NSBH': Neutron Star-Black Hole

              default: 'BBH'













      ..
          !! processed by numpydoc !!

   .. py:property:: lens_type

      
      Type of lens galaxy model.



      :Returns:

          **lens_type** : ``str``
              Type of lens model.

              Options:

              - 'epl_shear_galaxy': Elliptical Power Law with external shear

              - 'sie_galaxy': Singular Isothermal Ellipsoid

              - 'sis_galaxy': Singular Isothermal Sphere

              default: 'epl_shear_galaxy'













      ..
          !! processed by numpydoc !!

   .. py:property:: cosmo

      
      Astropy cosmology object for distance calculations.



      :Returns:

          **cosmo** : ``astropy.cosmology``
              Cosmology used for luminosity distance and comoving volume calculations.

              default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)













      ..
          !! processed by numpydoc !!

   .. py:property:: json_file_names

      
      Dictionary of JSON file names for parameter storage.



      :Returns:

          **json_file_names** : ``dict``
              Dictionary with keys:

              - 'ler_params': LeR initialization parameters

              - 'unlensed_param': Unlensed event parameters

              - 'unlensed_param_detectable': Detectable unlensed events

              - 'lensed_param': Lensed event parameters

              - 'lensed_param_detectable': Detectable lensed events













      ..
          !! processed by numpydoc !!

   .. py:property:: interpolator_directory

      
      Directory path for interpolator JSON files.



      :Returns:

          **interpolator_directory** : ``str``
              Path to directory containing interpolator data files.

              default: './interpolator_json'













      ..
          !! processed by numpydoc !!

   .. py:property:: ler_directory

      
      Directory path for LeR output files.



      :Returns:

          **ler_directory** : ``str``
              Path to directory for storing output parameter files.

              default: './ler_data'













      ..
          !! processed by numpydoc !!

   .. py:property:: create_new_interpolator

      
      Configuration dictionary for interpolator creation settings.



      :Returns:

          **create_new_interpolator** : ``dict``
              Dictionary specifying which interpolators to create.

              Each key is an interpolator name, and values are dicts with:

              - 'create_new': bool - Whether to create new interpolator

              - 'resolution': int or list - Grid resolution for interpolation

              Special key 'gwsnr' is a bool for GWSNR interpolator creation.
              Default: dict(
                  merger_rate_density = {'create_new': False, 'resolution': 500},
                  redshift_distribution = {'create_new': False, 'resolution': 500},
                  luminosity_distance = {'create_new': False, 'resolution': 500},
                  differential_comoving_volume = {'create_new': False, 'resolution': 500},
                  source_frame_masses = {'create_new': False, 'resolution': 500},
                  geocent_time = {'create_new': False, 'resolution': 500},
                  ra = {'create_new': False, 'resolution': 500},
                  dec = {'create_new': False, 'resolution': 500},
                  phase = {'create_new': False, 'resolution': 500},
                  psi = {'create_new': False, 'resolution': 500},
                  theta_jn = {'create_new': False, 'resolution': 500},
                  a_1 = {'create_new': False, 'resolution': 500},
                  a_2 = {'create_new': False, 'resolution': 500},
                  tilt_1 = {'create_new': False, 'resolution': 500},
                  tilt_2 = {'create_new': False, 'resolution': 500},
                  phi_12 = {'create_new': False, 'resolution': 500},
                  phi_jl = {'create_new': False, 'resolution': 500},
                  velocity_dispersion = {'create_new': False, 'resolution': 500, 'zl_resolution': 48},
                  axis_ratio = {'create_new': False, 'resolution': 500, 'sigma_resolution': 48},
                  lens_redshift = {'create_new': False, 'resolution': 48, 'zl_resolution': 48},
                  lens_redshift_intrinsic = {'create_new': False, 'resolution': 500},
                  optical_depth = {'create_new': False, 'resolution': 48},
                  comoving_distance = {'create_new': False, 'resolution': 500},
                  angular_diameter_distance = {'create_new': False, 'resolution': 500},
                  angular_diameter_distance_z1z2 = {'create_new': False, 'resolution': 500},
                  density_profile_slope = {'create_new': False, 'resolution': 100},
                  lens_parameters_kde_sl = {'create_new': False, 'resolution': 5000},
                  cross_section = {'create_new': False, 'resolution': [25, 25, 45, 15, 15]},
                  gwsnr = False,
              )













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

   .. py:property:: ler_args

      
      Dictionary of all LeR initialization arguments.



      :Returns:

          **ler_args** : ``dict``
              Dictionary containing all parameters used to initialize LeR and

              its parent classes, useful for reproducibility.













      ..
          !! processed by numpydoc !!

   .. py:method:: unlensed_cbc_statistics(size=100000, batch_size=50000, resume=True, save_batch=False, output_jsonfile=None)

      
      Generate unlensed GW source parameters.

      This function calls the unlensed_sampling_routine function to generate
      the parameters in batches. The generated parameters are stored in a JSON
      file; and if save_batch=True, it keeps updating the file in batches.

      :Parameters:

          **size** : ``int``
              Number of samples to generate.

              default: 100000

          **batch_size** : ``int``
              Batch size for sampling.

              default: 50000

          **resume** : ``bool``
              If True, the function will resume from the last batch.

              default: True

          **save_batch** : ``bool``
              If True, saves parameters in batches during sampling.

              If False, saves all parameters at the end (faster).

              default: False

          **output_jsonfile** : ``str``
              JSON file name for storing the parameters.

              default: None (uses self.json_file_names["unlensed_param"])

      :Returns:

          **unlensed_param** : ``dict``
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
              | a_1                |              | spin_1 of the compact binary         |
              +--------------------+--------------+--------------------------------------+
              | a_2                |              | spin_2 of the compact binary         |
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
              | pdet_L1            |              | pdet of L1                           |
              +--------------------+--------------+--------------------------------------+
              | pdet_H1            |              | pdet of H1                           |
              +--------------------+--------------+--------------------------------------+
              | pdet_V1            |              | pdet of V1                           |
              +--------------------+--------------+--------------------------------------+
              | pdet_net           |              | pdet of the network                  |
              +--------------------+--------------+--------------------------------------+










      .. rubric:: Examples

      >>> from ler import LeR
      >>> ler = LeR()
      >>> unlensed_param = ler.unlensed_cbc_statistics()



      ..
          !! processed by numpydoc !!

   .. py:method:: unlensed_sampling_routine(size, output_jsonfile, resume=True, save_batch=True)

      
      Generate unlensed GW source parameters for a single batch.

      This is the core sampling routine called by unlensed_cbc_statistics.
      It samples GW source parameters and calculates detection probabilities.

      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **output_jsonfile** : ``str``
              JSON file name for storing the parameters.

          **resume** : ``bool``
              If True, appends new samples to existing JSON file.

              default: True

          **save_batch** : ``bool``
              If True, saves parameters in batches during sampling.

              default: True

      :Returns:

          **unlensed_param** : ``dict``
              Dictionary of unlensed GW source parameters.













      ..
          !! processed by numpydoc !!

   .. py:method:: unlensed_rate(unlensed_param=None, pdet_threshold=0.5, pdet_type='boolean', output_jsonfile=None)

      
      Function to calculate the unlensed rate.

      This function calculates the detection rate for unlensed events and stores
      the parameters of the detectable events in a JSON file.

      :Parameters:

          **unlensed_param** : ``dict`` or ``str``
              Dictionary of GW source parameters or JSON file name.

              default: None (uses self.json_file_names["unlensed_param"])

          **pdet_threshold** : ``float``
              Threshold for detection probability.

              default: 0.5

          **pdet_type** : ``str``
              Detectability condition type.

              Options:

              - 'boolean': Binary detection based on pdet_threshold

              - 'probability_distribution': Uses pdet values directly

              default: 'boolean'

          **output_jsonfile** : ``str``
              JSON file name for storing the parameters of the detectable events.

              default: None (uses self.json_file_names["unlensed_param_detectable"])

      :Returns:

          **total_rate** : ``float``
              Total unlensed rate (yr^-1).

          **unlensed_param** : ``dict``
              Dictionary of unlensed GW source parameters of the detectable events.










      .. rubric:: Examples

      >>> from ler import LeR
      >>> ler = LeR()
      >>> ler.unlensed_cbc_statistics()
      >>> total_rate, unlensed_param_detectable = ler.unlensed_rate()



      ..
          !! processed by numpydoc !!

   .. py:method:: rate_function(detectable_size, total_size, param_type='unlensed', verbose=True)

      
      Calculate the detection rate for unlensed or lensed events.

      This is a general helper function that computes the rate based on
      Monte Carlo integration using the ratio of detectable to total events.

      :Parameters:

          **detectable_size** : ``int`` or ``float``
              Number of detectable events (or sum of pdet values).

          **total_size** : ``int``
              Total number of simulated events.

          **param_type** : ``str``
              Type of parameters.

              Options:

              - 'unlensed': Use unlensed normalization

              - 'lensed': Use lensed normalization

              default: 'unlensed'

          **verbose** : ``bool``
              If True, print rate information.

              default: True

      :Returns:

          **rate** : ``float``
              Event rate (yr^-1).










      .. rubric:: Examples

      >>> from ler import LeR
      >>> ler = LeR()
      >>> rate = ler.rate_function(detectable_size=100, total_size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: lensed_cbc_statistics(size=100000, batch_size=50000, save_batch=False, resume=True, output_jsonfile=None)

      
      Generate lensed GW source parameters.

      This function calls the lensed_sampling_routine function to generate
      the parameters in batches. The generated parameters are stored in a JSON
      file; and if save_batch=True, it keeps updating the file in batches.

      :Parameters:

          **size** : ``int``
              Number of samples to generate.

              default: 100000

          **batch_size** : ``int``
              Batch size for sampling.

              default: 50000

          **save_batch** : ``bool``
              If True, saves parameters in batches during sampling.

              If False, saves all parameters at the end (faster).

              default: True

          **resume** : ``bool``
              If True, the function will resume from the last batch.

              default: True

          **output_jsonfile** : ``str``
              JSON file name for storing the parameters.

              default: None (uses self.json_file_names["lensed_param"])

      :Returns:

          **lensed_param** : ``dict``
              Dictionary of lensed GW source parameters. The included parameters and their units are as follows (for default settings):

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
              | effective_luminosity_distance| Mpc       | effective luminosity distance of the images           |
              |                              |           | luminosity_distance / sqrt(|magnifications_i|)        |
              +------------------------------+-----------+-------------------------------------------------------+
              | effective_geocent_time       | s         | effective GPS time of coalescence of the images       |
              |                              |           | geocent_time + time_delays_i                          |
              +------------------------------+-----------+-------------------------------------------------------+
              | effective_phase              | rad       | morse-phase-corrected phase                           |
              |                              |           | phi - morse_phase_i                                   |
              +------------------------------+-----------+-------------------------------------------------------+
              | effective_ra                 | rad       | RA of the image                                       |
              |                              |           | ra + (x0_image_positions_i - x_source)/cos(dec)       |
              +------------------------------+-----------+-------------------------------------------------------+
              | effective_dec                | rad       | Dec of the image                                      |
              |                              |           | dec + (x1_image_positions_i - y_source)               |
              +------------------------------+-----------+-------------------------------------------------------+
              | effective_geocent_time       | s         | effective GPS time of coalescence of the images       |
              +------------------------------+-----------+-------------------------------------------------------+
              | pdet_L1                      |           | detection probability of L1                           |
              +------------------------------+-----------+-------------------------------------------------------+
              | pdet_H1                      |           | detection probability of H1                           |
              +------------------------------+-----------+-------------------------------------------------------+
              | pdet_V1                      |           | detection probability of V1                           |
              +------------------------------+-----------+-------------------------------------------------------+
              | pdet_net                     |           | detection probability of the network                  |
              +------------------------------+-----------+-------------------------------------------------------+










      .. rubric:: Examples

      >>> from ler import LeR
      >>> ler = LeR()
      >>> lensed_param = ler.lensed_cbc_statistics()



      ..
          !! processed by numpydoc !!

   .. py:method:: lensed_sampling_routine(size, output_jsonfile, save_batch=True, resume=True)

      
      Generate lensed GW source parameters for a single batch.

      This is the core sampling routine called by lensed_cbc_statistics.
      It samples lens parameters, calculates image properties, and computes
      detection probabilities for the images of lensed events.

      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **output_jsonfile** : ``str``
              JSON file name for storing the parameters.

          **save_batch** : ``bool``
              If True, saves parameters in batches during sampling.

              default: True

          **resume** : ``bool``
              If True, appends new samples to existing JSON file.

              default: True

      :Returns:

          **lensed_param** : ``dict``
              Dictionary of lensed GW source parameters.













      ..
          !! processed by numpydoc !!

   .. py:method:: lensed_rate(lensed_param=None, pdet_threshold=[0.5, 0.5], num_img=[1, 1], output_jsonfile=None, nan_to_num=True, pdet_type='boolean')

      
      Function to calculate the lensed rate.

      This function calculates the detection rate for lensed events and stores
      the parameters of the detectable events in a JSON file.

      :Parameters:

          **lensed_param** : ``dict`` or ``str``
              Dictionary of lensed GW source parameters or JSON file name.

              default: None (uses self.json_file_names["lensed_param"])

          **pdet_threshold** : ``float`` or ``list``
              Threshold for detection probability.

              default: [0.5, 0.5]

          **num_img** : ``int`` or ``list``
              Number of images corresponding to the pdet_threshold.

              Together with pdet_threshold = [0.5, 0.5], it means that two images with pdet > 0.5.

              Same condition can also be represented by pdet_threshold = 0.5 and num_img = 2.

              default: [1, 1]

          **output_jsonfile** : ``str``
              JSON file name for storing the parameters of the detectable events.

              default: None (uses self.json_file_names["lensed_param_detectable"])

          **nan_to_num** : ``bool``
              If True, NaN values will be converted to 0.

              default: True

          **pdet_type** : ``str``
              Detectability condition type.

              Options:

              - 'boolean': Binary detection based on pdet_threshold

              - 'probability_distribution': Uses pdet values directly

              default: 'boolean'

      :Returns:

          **total_rate** : ``float``
              Total lensed rate (yr^-1).

          **lensed_param** : ``dict``
              Dictionary of lensed GW source parameters of the detectable events. The included parameters and their units are as follows (for default settings):

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
              | effective_luminosity_distance| Mpc       | effective luminosity distance of the images           |
              |                              |           | luminosity_distance / sqrt(|magnifications_i|)        |
              +------------------------------+-----------+-------------------------------------------------------+
              | effective_geocent_time       | s         | effective GPS time of coalescence of the images       |
              |                              |           | geocent_time + time_delays_i                          |
              +------------------------------+-----------+-------------------------------------------------------+
              | effective_phase              | rad       | morse-phase-corrected phase                           |
              |                              |           | phi - morse_phase_i                                   |
              +------------------------------+-----------+-------------------------------------------------------+
              | effective_ra                 | rad       | RA of the image                                       |
              |                              |           | ra + (x0_image_positions_i - x_source)/cos(dec)       |
              +------------------------------+-----------+-------------------------------------------------------+
              | effective_dec                | rad       | Dec of the image                                      |
              |                              |           | dec + (x1_image_positions_i - y_source)               |
              +------------------------------+-----------+-------------------------------------------------------+
              | pdet_L1                      |           | detection probability of L1                           |
              +------------------------------+-----------+-------------------------------------------------------+
              | pdet_H1                      |           | detection probability of H1                           |
              +------------------------------+-----------+-------------------------------------------------------+
              | pdet_V1                      |           | detection probability of V1                           |
              +------------------------------+-----------+-------------------------------------------------------+
              | pdet_net                     |           | detection probability of the network                  |
              +------------------------------+-----------+-------------------------------------------------------+










      .. rubric:: Examples

      >>> from ler import LeR
      >>> ler = LeR()
      >>> ler.lensed_cbc_statistics()
      >>> total_rate, lensed_param_detectable = ler.lensed_rate()



      ..
          !! processed by numpydoc !!

   .. py:method:: rate_comparison_with_rate_calculation(unlensed_param=None, pdet_threshold_unlensed=0.5, output_jsonfile_unlensed=None, lensed_param=None, pdet_threshold_lensed=[0.5, 0.5], num_img_lensed=[1, 1], output_jsonfile_lensed=None, nan_to_num=True, pdet_type='boolean')

      
      Calculate and compare unlensed and lensed detection rates.

      This function calculates both unlensed and lensed rates and computes
      their ratio. It stores the parameters of the detectable events in JSON
      files. Using this function eliminates the need to call unlensed_rate
      and lensed_rate separately.

      :Parameters:

          **unlensed_param** : ``dict`` or ``str``
              Dictionary of GW source parameters or JSON file name.

              default: None (uses self.json_file_names["unlensed_param"])

          **pdet_threshold_unlensed** : ``float``
              Detection probability threshold for unlensed events.

              default: 0.5

          **output_jsonfile_unlensed** : ``str``
              JSON file name for storing detectable unlensed parameters.

              default: None

          **lensed_param** : ``dict`` or ``str``
              Dictionary of lensed GW source parameters or JSON file name.

              default: None (uses self.json_file_names["lensed_param"])

          **pdet_threshold_lensed** : ``float`` or ``list``
              Detection probability threshold for lensed events.

              default: [0.5, 0.5]

          **num_img_lensed** : ``list``
              Number of images for lensed events.

              default: [1, 1]

          **output_jsonfile_lensed** : ``str``
              JSON file name for storing detectable lensed parameters.

              default: None

          **nan_to_num** : ``bool``
              If True, NaN values will be converted to 0.

              default: True

          **pdet_type** : ``str``
              Detectability condition type.

              Options:

              - 'boolean': Binary detection based on pdet_threshold

              - 'probability_distribution': Uses pdet values directly

              default: 'boolean'

      :Returns:

          **rate_ratio** : ``float``
              Ratio of unlensed rate to lensed rate.

          **unlensed_param_detectable** : ``dict``
              Dictionary of detectable unlensed GW source parameters.

          **lensed_param_detectable** : ``dict``
              Dictionary of detectable lensed GW source parameters.










      .. rubric:: Examples

      >>> from ler import LeR
      >>> ler = LeR()
      >>> ler.unlensed_cbc_statistics()
      >>> ler.lensed_cbc_statistics()
      >>> rate_ratio, unlensed_param, lensed_param = ler.rate_comparison_with_rate_calculation()



      ..
          !! processed by numpydoc !!

   .. py:method:: rate_ratio()

      
      Calculate and display the unlensed to lensed merger rate ratio.

      This function retrieves the unlensed_rate and lensed_rate from the
      JSON file specified in self.json_file_names["ler_params"] and computes
      their ratio.


      :Returns:

          **rate_ratio** : ``float``
              Ratio of unlensed rate to lensed rate.










      .. rubric:: Examples

      >>> from ler import LeR
      >>> ler = LeR()
      >>> ler.unlensed_cbc_statistics()
      >>> ler.lensed_cbc_statistics()
      >>> ler.unlensed_rate()
      >>> ler.lensed_rate()
      >>> ler.rate_ratio()



      ..
          !! processed by numpydoc !!

   .. py:method:: selecting_n_unlensed_detectable_events(size=100, batch_size=50000, stopping_criteria=dict(relative_diff_percentage=0.5, number_of_last_batches_to_check=4), pdet_threshold=0.5, resume=True, output_jsonfile='n_unlensed_param_detectable.json', meta_data_file='meta_unlensed.json', pdet_type='boolean', trim_to_size=False)

      
      Generate a target number of detectable unlensed events by sampling in batches, with the option to stop once the cumulative rate has stabilized.

      This function samples unlensed parameters and saves only the detectable
      events in a JSON file. It also records metadata including the total
      number of events and the cumulative rate.

      :Parameters:

          **size** : ``int``
              Target number of detectable samples to collect.

              default: 100

          **batch_size** : ``int``
              Batch size for sampling.

              default: 50000

          **stopping_criteria** : ``dict`` or ``None``
              Criteria for stopping sample collection (but will not stop until n>size).

              Keys:

              - 'relative_diff_percentage': Maximum relative difference in rate (float)

              - 'number_of_last_batches_to_check': Number of batches for comparison (int)

              If None, stops when detectable events exceed size.

              default: dict(relative_diff_percentage=0.5, number_of_last_batches_to_check=4)

          **pdet_threshold** : ``float``
              Detection probability threshold.

              default: 0.5

          **resume** : ``bool``
              If True, resumes from last saved batch.

              default: True

          **output_jsonfile** : ``str``
              JSON file name for storing detectable parameters.

              default: 'n_unlensed_param_detectable.json'

          **meta_data_file** : ``str``
              JSON file name for storing metadata.

              default: 'meta_unlensed.json'

          **pdet_type** : ``str``
              Detectability condition type.

              Options:

              - 'boolean': Binary detection based on pdet_threshold

              - 'probability_distribution': Uses pdet values directly

              default: 'boolean'

          **trim_to_size** : ``bool``
              If True, trims final result to exactly size events.

              default: False

      :Returns:

          **param_final** : ``dict``
              Dictionary of unlensed GW source parameters of detectable events.










      .. rubric:: Examples

      >>> from ler import LeR
      >>> ler = LeR()
      >>> unlensed_param = ler.selecting_n_unlensed_detectable_events(size=100)



      ..
          !! processed by numpydoc !!

   .. py:method:: selecting_n_lensed_detectable_events(size=100, stopping_criteria=dict(relative_diff_percentage=2, number_of_last_batches_to_check=4), batch_size=50000, pdet_threshold=[0.5, 0.5], num_img=[1, 1], resume=True, pdet_type='boolean', output_jsonfile='n_lensed_params_detectable.json', meta_data_file='meta_lensed.json', trim_to_size=False, nan_to_num=True)

      
      Generate a target number of detectable lensed events by sampling in batches, with the option to stop once the cumulative rate has stabilized.

      This function samples lensed parameters and saves only the detectable
      events in a JSON file. It also records metadata including the total
      number of events and the cumulative rate.

      :Parameters:

          **size** : ``int``
              Target number of detectable samples to collect.

              default: 100

          **stopping_criteria** : ``dict`` or ``None``
              Criteria for stopping sample collection (but will not stop until n>size).

              Keys:

              - 'relative_diff_percentage': Maximum relative difference in rate (float)

              - 'number_of_last_batches_to_check': Number of batches for comparison (int)

              If None, stops when detectable events exceed size.

              default: dict(relative_diff_percentage=2, number_of_last_batches_to_check=4)

          **batch_size** : ``int``
              Batch size for sampling.

              default: 50000

          **pdet_threshold** : ``float`` or ``list``
              Detection probability threshold.

              default: [0.5, 0.5]

          **num_img** : ``list``
              Number of images corresponding to each pdet_threshold.

              default: [1, 1]

          **resume** : ``bool``
              If True, resumes from last saved batch.

              default: True

          **pdet_type** : ``str``
              Detectability condition type.

              Options:

              - 'boolean': Binary detection based on pdet_threshold

              - 'probability_distribution': Uses pdet values directly

              default: 'boolean'

          **output_jsonfile** : ``str``
              JSON file name for storing detectable parameters.

              default: 'n_lensed_params_detectable.json'

          **meta_data_file** : ``str``
              JSON file name for storing metadata.

              default: 'meta_lensed.json'

          **trim_to_size** : ``bool``
              If True, trims final result to exactly size events.

              default: False

          **nan_to_num** : ``bool``
              If True, NaN values will be converted to 0.

              default: False

      :Returns:

          **param_final** : ``dict``
              Dictionary of lensed GW source parameters of detectable events.










      .. rubric:: Examples

      >>> from ler import LeR
      >>> ler = LeR()
      >>> lensed_param = ler.selecting_n_lensed_detectable_events(size=100)



      ..
          !! processed by numpydoc !!


.. py:class:: CBCSourceParameterDistribution(z_min=0.0, z_max=10.0, event_type='BBH', source_priors=None, source_priors_params=None, cosmology=None, spin_zero=False, spin_precession=False, directory='./interpolator_json', create_new_interpolator=False)


   Bases: :py:obj:`ler.gw_source_population.cbc_source_redshift_distribution.CBCSourceRedshiftDistribution`

   
   Class for sampling compact binary coalescence source parameters.

   This class generates complete sets of intrinsic and extrinsic gravitational
   wave parameters for compact binary sources including masses, spins, sky
   positions, and orbital parameters. It supports BBH, BNS, NSBH, and primordial
   black hole populations with configurable prior distributions.

   Key Features:

   - Multiple mass distribution models (PowerLaw+Gaussian, lognormal, bimodal)

   - Configurable spin priors (zero, aligned, precessing)

   - Isotropic sky position and orientation sampling

   - Built-in support for population III and primordial black holes

   :Parameters:

       **z_min** : ``float``
           Minimum redshift of the source population.

           default: 0.0

       **z_max** : ``float``
           Maximum redshift of the source population.

           default: 10.0

       **event_type** : ``str``
           Type of compact binary event to generate.

           Options:

           - 'BBH': Binary black hole (Population I/II)

           - 'BNS': Binary neutron star

           - 'NSBH': Neutron star-black hole

           - 'BBH_popIII': Population III binary black hole

           - 'BBH_primordial': Primordial binary black hole

           default: 'BBH'

       **source_priors** : ``dict`` or ``None``
           Dictionary of prior sampler functions for each parameter.

           If None, uses default priors based on event_type.

           default: None

       **source_priors_params** : ``dict`` or ``None``
           Dictionary of parameters for each prior sampler function.

           If None, uses default parameters based on event_type.

           default: None

       **cosmology** : ``astropy.cosmology`` or ``None``
           Cosmology to use for distance calculations.

           default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)

       **spin_zero** : ``bool``
           If True, spin parameters are set to zero (no spin sampling).

           default: False

       **spin_precession** : ``bool``
           If True (and spin_zero=False), sample precessing spin parameters.

           If False (and spin_zero=False), sample aligned/anti-aligned spins.

           default: False

       **directory** : ``str``
           Directory to store interpolator JSON files.

           default: './interpolator_json'

       **create_new_interpolator** : ``dict`` or ``bool``
           Configuration for creating new interpolators.

           If bool, applies to all interpolators.

           default: False











   .. rubric:: Examples

   >>> from ler.gw_source_population import CBCSourceParameterDistribution
   >>> cbc = CBCSourceParameterDistribution(event_type='BBH')
   >>> params = cbc.sample_gw_parameters(size=1000)
   >>> print(list(params.keys()))

   Instance Methods
   ----------
   CBCSourceParameterDistribution has the following methods:

   +-----------------------------------------------------+------------------------------------------------+
   | Method                                              | Description                                    |
   +=====================================================+================================================+
   | :meth:`~sample_gw_parameters`                       | Sample all GW parameters for compact binaries  |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~binary_masses_BBH_powerlaw_gaussian`| Sample BBH masses with PowerLaw+PEAK model     |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~binary_masses_BBH_popIII_lognormal`         | Sample pop III BBH masses from lognormal       |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~binary_masses_BBH_primordial_lognormal`     | Sample primordial BBH masses from lognormal    |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~binary_masses_NSBH_broken_powerlaw`         | Sample NSBH masses from broken powerlaw        |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~binary_masses_uniform`                      | Sample masses from uniform distribution        |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~binary_masses_BNS_bimodal`                  | Sample BNS masses from bimodal distribution    |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~constant_values_n_size`                     | Return array of constant values                |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~sampler_uniform`                            | Sample from uniform distribution               |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~sampler_cosine`                             | Sample from cosine distribution                |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~sampler_sine`                               | Sample from sine distribution                  |
   +-----------------------------------------------------+------------------------------------------------+

   Instance Attributes
   ----------
   CBCSourceParameterDistribution has the following attributes:

   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | Attribute                                      | Type                   | Unit  | Description                                    |
   +================================================+========================+=======+================================================+
   | :attr:`~z_min`                                 | ``float``              |       | Minimum redshift of source population          |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~z_max`                                 | ``float``              |       | Maximum redshift of source population          |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~cosmo`                                 | ``astropy.cosmology``  |       | Cosmology for distance calculations            |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~spin_zero`                             | ``bool``               |       | Whether to ignore spin parameters              |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~spin_precession`                       | ``bool``               |       | Whether to use precessing spins                |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~directory`                             | ``str``                |       | Directory for interpolator files               |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~gw_param_samplers`                     | ``dict``               |       | Dictionary of parameter sampler functions      |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~gw_param_samplers_params`              | ``dict``               |       | Dictionary of sampler function parameters      |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~available_gw_prior`                    | ``dict``               |       | Available prior distributions                  |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~source_frame_masses`                   | ``callable``           |       | Sampler for source frame masses                |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~zs`                                    | ``callable``           |       | Sampler for source redshift                    |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~geocent_time`                          | ``callable``           |       | Sampler for geocentric time                    |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~ra`                                    | ``callable``           |       | Sampler for right ascension                    |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~dec`                                   | ``callable``           |       | Sampler for declination                        |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~phase`                                 | ``callable``           |       | Sampler for coalescence phase                  |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~psi`                                   | ``callable``           |       | Sampler for polarization angle                 |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~theta_jn`                              | ``callable``           |       | Sampler for inclination angle                  |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~a_1`                                   | ``callable``           |       | Sampler for spin1 magnitude                    |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~a_2`                                   | ``callable``           |       | Sampler for spin2 magnitude                    |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~tilt_1`                                | ``callable``           |       | Sampler for tilt1 angle                        |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~tilt_2`                                | ``callable``           |       | Sampler for tilt2 angle                        |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~phi_12`                                | ``callable``           |       | Sampler for phi_12 angle                       |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~phi_jl`                                | ``callable``           |       | Sampler for phi_jl angle                       |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:property:: zs

      
      Class object (of FunctionConditioning) for source redshift, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the redshift distribution

      - `pdf`: returns the probability density function of the redshift distribution

      - `function`: returns the redshift distribution function.


      :Returns:

          **zs** : ``numpy.ndarray``
              Array of redshift values.













      ..
          !! processed by numpydoc !!

   .. py:property:: source_frame_masses

      
      Class object (of FunctionConditioning) for source frame masses, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the density profile slope distribution


      :Returns:

          **mass_1_source** : ``numpy.ndarray``
              Array of mass_1_source values in solar masses.

          **mass_2_source** : ``numpy.ndarray``
              Array of mass_2_source values in solar masses.










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc_source_param_dist = CBCSourceParameterDistribution()
      >>> cbc_source_param_dist.source_frame_masses(size=10)



      ..
          !! processed by numpydoc !!

   .. py:property:: geocent_time

      
      Class object (of FunctionConditioning) for geocentric time, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the geocentric time distribution

      - `pdf`: returns the probability density function of the geocentric time distribution

      - `function`: returns the geocentric time distribution function.


      :Returns:

          **geocent_time** : ``numpy.ndarray``
              Array of geocentric time values.













      ..
          !! processed by numpydoc !!

   .. py:property:: ra

      
      Class object (of FunctionConditioning) for right ascension, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the right ascension distribution

      - `pdf`: returns the probability density function of the right ascension distribution

      - `function`: returns the right ascension distribution function.


      :Returns:

          **ra** : ``numpy.ndarray``
              Array of right ascension values.













      ..
          !! processed by numpydoc !!

   .. py:property:: dec

      
      Class object (of FunctionConditioning) for declination, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the declination distribution

      - `pdf`: returns the probability density function of the declination distribution

      - `function`: returns the declination distribution function.


      :Returns:

          **dec** : ``numpy.ndarray``
              Array of declination values.













      ..
          !! processed by numpydoc !!

   .. py:property:: phase

      
      Class object (of FunctionConditioning) for coalescence phase, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the coalescence phase distribution

      - `pdf`: returns the probability density function of the coalescence phase distribution

      - `function`: returns the coalescence phase distribution function.


      :Returns:

          **phase** : ``numpy.ndarray``
              Array of coalescence phase values.













      ..
          !! processed by numpydoc !!

   .. py:property:: psi

      
      Class object (of FunctionConditioning) for polarization angle, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the polarization angle distribution

      - `pdf`: returns the probability density function of the polarization angle distribution

      - `function`: returns the polarization angle distribution function.


      :Returns:

          **geocent_time** : ``numpy.ndarray``
              Array of polarization angle values.













      ..
          !! processed by numpydoc !!

   .. py:property:: theta_jn

      
      Class object (of FunctionConditioning) for inclination angle, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the inclination angle distribution

      - `pdf`: returns the probability density function of the inclination angle distribution

      - `function`: returns the inclination angle distribution function.


      :Returns:

          **theta_jn** : ``numpy.ndarray``
              Array of inclination angle values, i.e. the angle between the line of sight and the orbital angular momentum (rad).













      ..
          !! processed by numpydoc !!

   .. py:property:: a_1

      
      Class object (of FunctionConditioning) for spin1 magnitude, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the spin1 magnitude distribution

      - `pdf`: returns the probability density function of the spin1 magnitude distribution

      - `function`: returns the spin1 magnitude distribution function.


      :Returns:

          **a_1** : ``numpy.ndarray``
              Array of spin magnitude values for the primary body.













      ..
          !! processed by numpydoc !!

   .. py:property:: a_2

      
      Class object (of FunctionConditioning) for spin2 magnitude, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the spin2 magnitude distribution

      - `pdf`: returns the probability density function of the spin2 magnitude distribution

      - `function`: returns the spin2 magnitude distribution function.


      :Returns:

          **a_2** : ``numpy.ndarray``
              Array of spin magnitude values for the secondary body.













      ..
          !! processed by numpydoc !!

   .. py:property:: tilt_1

      
      Class object (of FunctionConditioning) for tilt1 angle, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the tilt1 angle distribution

      - `pdf`: returns the probability density function of the tilt1 angle distribution

      - `function`: returns the tilt1 angle distribution function.


      :Returns:

          **tilt_1** : ``numpy.ndarray``
              Array of the spin tilt angle of the primary body, i.e. the angle between the spin vector and the orbital angular momentum for the primary body (rad).













      ..
          !! processed by numpydoc !!

   .. py:property:: tilt_2

      
      Class object (of FunctionConditioning) for tilt2 angle, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the tilt2 angle distribution

      - `pdf`: returns the probability density function of the tilt2 angle distribution

      - `function`: returns the tilt2 angle distribution function.


      :Returns:

          **tilt_2** : ``numpy.ndarray``
              Array of the spin tilt angle of the secondary body, i.e. the angle between the spin vector and the orbital angular momentum for the secondary body (rad).













      ..
          !! processed by numpydoc !!

   .. py:property:: phi_12

      
      Class object (of FunctionConditioning) for phi_12 angle, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the phi_12 angle distribution

      - `pdf`: returns the probability density function of the phi_12 angle distribution

      - `function`: returns the phi_12 angle distribution function.


      :Returns:

          **phi_12** : ``numpy.ndarray``
              Array of the spin tilt angle between the two spins, i.e., angle between the projections of the two spins onto the orbital plane (rad).













      ..
          !! processed by numpydoc !!

   .. py:property:: phi_jl

      
      Class object (of FunctionConditioning) for phi_jl angle, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the phi_jl angle distribution

      - `pdf`: returns the probability density function of the phi_jl angle distribution

      - `function`: returns the phi_jl angle distribution function.


      :Returns:

          **phi_jl** : ``numpy.ndarray``
              Array of the angle values between the orientation of the total angular momentum around the orbital angular momentum (rad).













      ..
          !! processed by numpydoc !!

   .. py:property:: available_gw_prior

      
      Dictionary of all available prior distributions and their parameters.

      This is a dynamically generated dictionary containing available samplers
      for each GW parameter type and their default parameter values.


      :Returns:

          **available_gw_prior** : ``dict``
              Nested dictionary organized by parameter type (e.g., 'source_frame_masses',

              'geocent_time', etc.) with sampler names and default parameters.













      ..
          !! processed by numpydoc !!

   .. py:attribute:: z_min
      :value: 'None'

      
      ``float``

      Minimum redshift of the source population















      ..
          !! processed by numpydoc !!

   .. py:attribute:: z_max
      :value: 'None'

      
      ``float``

      Maximum redshift of the source population















      ..
          !! processed by numpydoc !!

   .. py:attribute:: event_type
      :value: 'None'

      
      ``str``

      Type of event to generate.

      e.g. 'BBH', 'BNS', 'NSBH'















      ..
          !! processed by numpydoc !!

   .. py:attribute:: source_priors
      :value: 'None'

      
      ``dict``

      Dictionary of prior sampler functions.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: source_priors_params
      :value: 'None'

      
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
      :value: 'None'

      
      ``bool``

      If True, spin prior is set to zero.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: spin_precession
      :value: 'False'

      

   .. py:attribute:: directory
      :value: "'./interpolator_json'"

      
      Directory path for storing interpolator JSON files.



      :Returns:

          **directory** : ``str``
              Path to the interpolator storage directory.

              default: './interpolator_json'













      ..
          !! processed by numpydoc !!

   .. py:method:: sample_gw_parameters(size=1000, param=None)

      
      Sample all gravitational wave parameters for compact binaries.

      Generates a complete set of intrinsic and extrinsic parameters including
      masses, redshift, luminosity distance, sky position, orientation, and
      optionally spin parameters.

      :Parameters:

          **size** : ``int``
              Number of samples to draw.

              default: 1000

          **param** : ``dict`` or ``None``
              Dictionary of fixed parameter values.

              Parameters in this dict will not be sampled.

              default: None

      :Returns:

          **gw_parameters** : ``dict``
              Dictionary of sampled GW parameters. The included parameters and their units are as follows (for default settings):

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
              | a_1                |              | spin_1 of the compact binary         |
              +--------------------+--------------+--------------------------------------+
              | a_2                |              | spin_2 of the compact binary         |
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










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> params = cbc.sample_gw_parameters(size=1000)
      >>> print(list(params.keys()))



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BBH_powerlaw_gaussian(size, get_attribute=False, **kwargs)

      
      Sample source masses with PowerLaw+PEAK model for Population I/II BBH.

      Implements the mass distribution model from LIGO-Virgo population analyses
      combining a power-law with a Gaussian peak component.

      :Parameters:

          **size** : ``int``
              Number of samples to draw.

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Model parameters:

              - mminbh: Minimum BH mass (Msun), default: 4.98

              - mmaxbh: Maximum BH mass (Msun), default: 112.5

              - alpha: Power-law spectral index, default: 3.78

              - mu_g: Gaussian peak mean (Msun), default: 32.27

              - sigma_g: Gaussian peak width (Msun), default: 3.88

              - lambda_peak: Fraction in Gaussian component, default: 0.03

              - delta_m: Low-mass tapering range (Msun), default: 4.8

              - beta: Mass ratio power-law index, default: 0.81

      :Returns:

          **mass_1_source** : ``numpy.ndarray``
              Array of primary masses in source frame (Msun).

          **mass_2_source** : ``numpy.ndarray``
              Array of secondary masses in source frame (Msun).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> m1_src, m2_src = cbc.binary_masses_BBH_powerlaw_gaussian(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BBH_popIII_lognormal(size, get_attribute=False, **kwargs)

      
      Sample source masses for Population III BBH from lognormal distribution.

      Based on Eqn. 1 and 4 of Ng et al. 2022 for Population III black holes.

      :Parameters:

          **size** : ``int``
              Number of samples to draw.

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Model parameters:

              - m_min: Minimum BH mass (Msun), default: 5.0

              - m_max: Maximum BH mass (Msun), default: 150.0

              - Mc: Central mass scale (Msun), default: 30.0

              - sigma: Distribution width, default: 0.3

      :Returns:

          **mass_1_source** : ``numpy.ndarray``
              Array of primary masses in source frame (Msun).

          **mass_2_source** : ``numpy.ndarray``
              Array of secondary masses in source frame (Msun).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution(event_type='BBH_popIII')
      >>> m1_src, m2_src = cbc.binary_masses_BBH_popIII_lognormal(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BBH_primordial_lognormal(size, get_attribute=False, **kwargs)

      
      Sample source masses for primordial BBH from lognormal distribution.

      Based on Eqn. 1 and 4 of Ng et al. 2022 for primordial black holes.

      :Parameters:

          **size** : ``int``
              Number of samples to draw.

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Model parameters:

              - m_min: Minimum BH mass (Msun), default: 1.0

              - m_max: Maximum BH mass (Msun), default: 100.0

              - Mc: Central mass scale (Msun), default: 20.0

              - sigma: Distribution width, default: 0.3

      :Returns:

          **mass_1_source** : ``numpy.ndarray``
              Array of primary masses in source frame (Msun).

          **mass_2_source** : ``numpy.ndarray``
              Array of secondary masses in source frame (Msun).













      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_NSBH_broken_powerlaw(size, get_attribute=False, **kwargs)

      
      Sample source masses for NSBH from broken power-law distribution.

      Uses gwcosmo-style broken power-law for black hole mass and power-law
      for neutron star mass.

      :Parameters:

          **size** : ``int``
              Number of samples to draw.

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Model parameters:

              - mminbh: Minimum BH mass (Msun), default: 26

              - mmaxbh: Maximum BH mass (Msun), default: 125

              - alpha_1: Primary power-law index, default: 6.75

              - alpha_2: Secondary power-law index, default: 6.75

              - b: Break point, default: 0.5

              - delta_m: Tapering range (Msun), default: 5

              - mminns: Minimum NS mass (Msun), default: 1.0

              - mmaxns: Maximum NS mass (Msun), default: 3.0

              - alphans: NS mass power-law index, default: 0.0

      :Returns:

          **mass_1_source** : ``numpy.ndarray``
              Array of BH masses in source frame (Msun).

          **mass_2_source** : ``numpy.ndarray``
              Array of NS masses in source frame (Msun).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution(event_type='NSBH')
      >>> m1_src, m2_src = cbc.binary_masses_NSBH_broken_powerlaw(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_uniform(size, get_attribute=False, **kwargs)

      
      Sample source masses from uniform distribution.


      :Parameters:

          **size** : ``int``
              Number of samples to draw.

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Model parameters:

              - m_min: Minimum mass (Msun), default: 1.0

              - m_max: Maximum mass (Msun), default: 3.0

      :Returns:

          **mass_1_source** : ``numpy.ndarray``
              Array of primary masses in source frame (Msun).

          **mass_2_source** : ``numpy.ndarray``
              Array of secondary masses in source frame (Msun).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> m1_src, m2_src = cbc.binary_masses_uniform(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BNS_bimodal(size, get_attribute=False, **kwargs)

      
      Sample BNS masses from bimodal Gaussian distribution.

      Based on Will M. Farr et al. 2020 Eqn. 6 for neutron star mass
      distribution combining two Gaussian peaks.

      :Parameters:

          **size** : ``int``
              Number of samples to draw.

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Model parameters:

              - w: Weight of left peak, default: 0.643

              - muL: Mean of left peak (Msun), default: 1.352

              - sigmaL: Width of left peak (Msun), default: 0.08

              - muR: Mean of right peak (Msun), default: 1.88

              - sigmaR: Width of right peak (Msun), default: 0.3

              - mmin: Minimum mass (Msun), default: 1.0

              - mmax: Maximum mass (Msun), default: 2.3

      :Returns:

          **mass_1_source** : ``numpy.ndarray``
              Array of primary masses in source frame (Msun).

          **mass_2_source** : ``numpy.ndarray``
              Array of secondary masses in source frame (Msun).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution(event_type='BNS')
      >>> m1_src, m2_src = cbc.binary_masses_BNS_bimodal(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: constant_values_n_size(size=100, get_attribute=False, **kwargs)

      
      Return array of constant values.


      :Parameters:

          **size** : ``int``
              Number of values to return.

              default: 100

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Model parameters:

              - value: Constant value to return, default: 0.0

      :Returns:

          **values** : ``numpy.ndarray``
              Array of constant values.













      ..
          !! processed by numpydoc !!

   .. py:method:: sampler_uniform(size, get_attribute=False, **kwargs)

      
      Sample values from uniform distribution.


      :Parameters:

          **size** : ``int``
              Number of samples to draw.

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Model parameters:

              - xmin: Minimum value, default: 0.0

              - xmax: Maximum value, default: 1.0

      :Returns:

          **values** : ``numpy.ndarray``
              Array of uniformly distributed values in range [xmin, xmax].













      ..
          !! processed by numpydoc !!

   .. py:method:: sampler_cosine(size, get_attribute=False, **kwargs)

      
      Sample from cosine distribution for declination.

      Samples values in range [-pi/2, pi/2] following a cosine distribution,
      appropriate for isotropic sky position declination.

      :Parameters:

          **size** : ``int``
              Number of samples to draw.

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

      :Returns:

          **values** : ``numpy.ndarray``
              Array of values in range [-pi/2, pi/2] (rad).













      ..
          !! processed by numpydoc !!

   .. py:method:: sampler_sine(size, get_attribute=False, **kwargs)

      
      Sample from sine distribution for inclination angles.

      Samples values in range [0, pi] following a sine distribution,
      appropriate for isotropic orientation angles.

      :Parameters:

          **size** : ``int``
              Number of samples to draw.

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

      :Returns:

          **values** : ``numpy.ndarray``
              Array of values in range [0, pi] (rad).













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

.. py:function:: batch_handler(size, batch_size, sampling_routine, output_jsonfile, save_batch=True, resume=False, param_name='parameters')

   
   Function to run the sampling in batches.


   :Parameters:

       **size** : `int`
           number of samples.

       **batch_size** : `int`
           batch size.

       **sampling_routine** : `function`
           sampling function. It should have 'size' as input and return a dictionary.

       **output_jsonfile** : `str`
           json file name for storing the parameters.

       **save_batch** : `bool`, optional
           if True, save sampled parameters in each iteration. Default is True.

       **resume** : `bool`, optional
           if True, resume sampling from the last batch. Default is False.

       **param_name** : `str`, optional
           name of the parameter. Default is 'parameters'.

   :Returns:

       **dict_buffer** : `dict`
           dictionary of parameters.













   ..
       !! processed by numpydoc !!

.. py:function:: remove_file(file_name)

   
   Remove a file.
















   ..
       !! processed by numpydoc !!

.. py:class:: GWRATES(npool=int(4), z_min=0.0, z_max=10.0, event_type='BBH', cosmology=None, pdet_finder=None, json_file_names=None, interpolator_directory='./interpolator_json', create_new_interpolator=False, ler_directory='./ler_data', verbose=True, **kwargs)


   Bases: :py:obj:`ler.gw_source_population.CBCSourceParameterDistribution`

   
   Class to sample GW events and calculate their detection rates.

   This class provides functionality for sampling gravitational wave source
   parameters, detection probabilities, and computing detection rates for
   compact binary coalescence events. Parameters of simulated events are
   stored in JSON files (not as class attributes) to conserve RAM memory.

   Key Features:

   - Sampling of GW event parameters

   - Detection probability calculation using ``gwsnr`` package or custom functions

   - Rate calculation for detectable events

   - Batch processing for memory efficiency

   - JSON-based parameter storage for reproducibility

   :Parameters:

       **npool** : ``int``
           Number of cores to use for parallel processing.

           default: 4

       **z_min** : ``float``
           Minimum redshift of the source population.

           default: 0.0

       **z_max** : ``float``
           Maximum redshift of the source population.

           default: 10.0

       **event_type** : ``str``
           Type of event to generate.

           Options:

           - 'BBH': Binary Black Hole

           - 'BNS': Binary Neutron Star

           - 'NSBH': Neutron Star-Black Hole

           default: 'BBH'

       **cosmology** : ``astropy.cosmology`` or ``None``
           Cosmology to use for the calculation.

           default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)

       **pdet_finder** : ``callable`` or ``None``
           Custom detection probability finder function.

           If None, uses gwsnr's pdet calculator.

           The function should follow the signature:

           ``def pdet_finder(gw_param_dict): return pdet_net_dict``

           where pdet_net_dict.keys = ['pdet_net'].

           default: None

       **json_file_names** : ``dict`` or ``None``
           Names of the JSON files to store the necessary parameters.

           default: dict(gwrates_params="gwrates_params.json", gw_param="gw_param.json", gw_param_detectable="gw_param_detectable.json")

       **interpolator_directory** : ``str``
           Directory to store the interpolators.

           default: './interpolator_json'

       **create_new_interpolator** : ``bool`` or ``dict``
           Whether to create new interpolators.

           default: False

       **ler_directory** : ``str``
           Directory to store the output parameters.

           default: './ler_data'

       **verbose** : ``bool``
           If True, print all chosen parameters during initialization.

           default: True

       **\*\*kwargs** : ``dict``
           Additional keyword arguments passed to parent classes:

           :class:`~ler.gw_source_population.CBCSourceParameterDistribution` and

           :class:`~gwsnr.GWSNR` (if pdet_finder is not provided).









   .. rubric:: Notes

   - ``GWRATES`` class inherits from :class:`~ler.gw_source_population.CBCSourceParameterDistribution`.

     Refer to that class for additional inherited attributes and methods.

   - Parameters are stored in JSON files for memory efficiency and reproducibility.

   - For stable rate estimates, use size >= 1e6 samples.


   .. rubric:: Examples

   Basic usage:

   >>> from ler.rates import GWRATES
   >>> gwrates = GWRATES()
   >>> gw_params = gwrates.gw_cbc_statistics()
   >>> gwrates.gw_rate()

   Instance Methods
   ----------
   GWRATES class has the following methods:

   +-----------------------------------------------------+------------------------------------------------+
   | Method                                              | Description                                    |
   +=====================================================+================================================+
   | :meth:`~gw_cbc_statistics`                          | Generate GW source parameters                  |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~gw_rate`                                    | Calculate the GW detection rate                |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~selecting_n_gw_detectable_events`           | Select n GW detectable events                  |
   +-----------------------------------------------------+------------------------------------------------+

   Instance Attributes
   ----------
   GWRATES class has the following attributes:

   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | Attribute                                      | Type             | Unit  | Description                                    |
   +================================================+==================+=======+================================================+
   | :meth:`~npool`                                 | ``int``          |       | Number of parallel processing cores            |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~z_min`                                 | ``float``        |       | Minimum source redshift                        |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~z_max`                                 | ``float``        |       | Maximum source redshift                        |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~event_type`                            | ``str``          |       | Type of CBC event                              |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~cosmo`                                 | ``astropy.cosmology`` |  | Cosmology for calculations                     |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~json_file_names`                       | ``dict``         |       | JSON file names for parameter storage          |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~interpolator_directory`                | ``str``          |       | Directory for interpolator files               |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~ler_directory`                         | ``str``          |       | Directory for output parameter files           |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~list_of_detectors`                     | ``list``         |       | List of detector names                         |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~pdet_finder`                           | ``callable``     |       | Detection probability finder function          |
   +------------------------------------------------+------------------+-------+------------------------------------------------+
   | :meth:`~gwrates_args`                          | ``dict``         |       | All GWRATES initialization arguments           |
   +------------------------------------------------+------------------+-------+------------------------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:property:: npool

      
      Number of parallel processing cores.



      :Returns:

          **npool** : ``int``
              Number of logical cores to use for multiprocessing.

              default: 4













      ..
          !! processed by numpydoc !!

   .. py:property:: z_min

      
      Minimum redshift of the source population.



      :Returns:

          **z_min** : ``float``
              Minimum source redshift for sampling.

              default: 0.0













      ..
          !! processed by numpydoc !!

   .. py:property:: z_max

      
      Maximum redshift of the source population.



      :Returns:

          **z_max** : ``float``
              Maximum source redshift for sampling.

              default: 10.0













      ..
          !! processed by numpydoc !!

   .. py:property:: event_type

      
      Type of compact binary coalescence event.



      :Returns:

          **event_type** : ``str``
              Type of CBC event.

              Options:

              - 'BBH': Binary Black Hole

              - 'BNS': Binary Neutron Star

              - 'NSBH': Neutron Star-Black Hole

              default: 'BBH'













      ..
          !! processed by numpydoc !!

   .. py:property:: cosmo

      
      Astropy cosmology object for distance calculations.



      :Returns:

          **cosmo** : ``astropy.cosmology``
              Cosmology used for luminosity distance and comoving volume calculations.

              default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)













      ..
          !! processed by numpydoc !!

   .. py:property:: json_file_names

      
      Dictionary of JSON file names for parameter storage.



      :Returns:

          **json_file_names** : ``dict``
              Dictionary with keys:

              - 'gwrates_params': GWRATES initialization parameters

              - 'gw_param': GW event parameters

              - 'gw_param_detectable': Detectable GW events













      ..
          !! processed by numpydoc !!

   .. py:property:: interpolator_directory

      
      Directory path for interpolator JSON files.



      :Returns:

          **interpolator_directory** : ``str``
              Path to directory containing interpolator data files.

              default: './interpolator_json'













      ..
          !! processed by numpydoc !!

   .. py:property:: ler_directory

      
      Directory path for LeR output files.



      :Returns:

          **ler_directory** : ``str``
              Path to directory for storing output parameter files.

              default: './ler_data'













      ..
          !! processed by numpydoc !!

   .. py:property:: create_new_interpolator

      
      Configuration dictionary for interpolator creation settings.



      :Returns:

          **create_new_interpolator** : ``dict``
              Dictionary specifying which interpolators to create.

              Each key is an interpolator name, and values are dicts with:

              - 'create_new': bool - Whether to create new interpolator

              - 'resolution': int or list - Grid resolution for interpolation

              Special key 'gwsnr' is a bool for GWSNR interpolator creation.
              Default: dict(
                  merger_rate_density = {'create_new': False, 'resolution': 500},
                  redshift_distribution = {'create_new': False, 'resolution': 500},
                  luminosity_distance = {'create_new': False, 'resolution': 500},
                  differential_comoving_volume = {'create_new': False, 'resolution': 500},
                  source_frame_masses = {'create_new': False, 'resolution': 500},
                  geocent_time = {'create_new': False, 'resolution': 500},
                  ra = {'create_new': False, 'resolution': 500},
                  dec = {'create_new': False, 'resolution': 500},
                  phase = {'create_new': False, 'resolution': 500},
                  psi = {'create_new': False, 'resolution': 500},
                  theta_jn = {'create_new': False, 'resolution': 500},
                  a_1 = {'create_new': False, 'resolution': 500},
                  a_2 = {'create_new': False, 'resolution': 500},
                  tilt_1 = {'create_new': False, 'resolution': 500},
                  tilt_2 = {'create_new': False, 'resolution': 500},
                  phi_12 = {'create_new': False, 'resolution': 500},
                  phi_jl = {'create_new': False, 'resolution': 500},
              )













      ..
          !! processed by numpydoc !!

   .. py:property:: list_of_detectors

      
      List of gravitational wave detector names.



      :Returns:

          **list_of_detectors** : ``list``
              List of detector identifiers used for pdet calculations.

              Typically set from gwsnr.detector_list during initialization.













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

   .. py:property:: gwrates_args

      
      Dictionary of all GWRATES initialization arguments.



      :Returns:

          **gwrates_args** : ``dict``
              Dictionary containing all parameters used to initialize GWRATES and

              its parent classes, useful for reproducibility.













      ..
          !! processed by numpydoc !!

   .. py:method:: gw_cbc_statistics(size=100000, batch_size=50000, resume=True, save_batch=False, output_jsonfile=None)

      
      Generate GW source parameters with detection probabilities.

      This function calls the _gw_sampling_routine function to generate
      the parameters in batches. The generated parameters are stored in
      a JSON file; if save_batch=True, it updates the file after each batch.

      :Parameters:

          **size** : ``int``
              Number of samples to generate.

              default: 100000

          **batch_size** : ``int``
              Batch size for sampling.

              default: 50000

          **resume** : ``bool``
              If True, resume from the last batch.

              default: True

          **save_batch** : ``bool``
              If True, save parameters after each batch.

              If False, save all parameters at the end (faster).

              default: False

          **output_jsonfile** : ``str`` or ``None``
              JSON file name for storing the parameters.

              default: 'gw_param.json' (stored in ler_directory)

      :Returns:

          **gw_param** : ``dict``
              dictionary of GW source parameters. The included parameters and their units are as follows (for default settings):

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
              | a_1                |              | spin of the primary compact binary         |
              +--------------------+--------------+--------------------------------------+
              | a_2                |              | spin_2 of the compact binary         |
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
              | pdet_L1            |              | pdet of L1                           |
              +--------------------+--------------+--------------------------------------+
              | pdet_H1            |              | pdet of H1                           |
              +--------------------+--------------+--------------------------------------+
              | pdet_V1            |              | pdet of V1                           |
              +--------------------+--------------+--------------------------------------+
              | pdet_net           |              | pdet of the network                  |
              +--------------------+--------------+--------------------------------------+










      .. rubric:: Examples

      >>> from ler.rates import GWRATES
      >>> gwrates = GWRATES()
      >>> param = gwrates.gw_cbc_statistics()



      ..
          !! processed by numpydoc !!

   .. py:method:: gw_rate(gw_param=None, pdet_threshold=0.5, pdet_type='boolean', output_jsonfile=None)

      
      Calculate the GW detection rate.

      This function calculates the detection rate and stores the parameters
      of detectable events in a JSON file.

      :Parameters:

          **gw_param** : ``dict`` or ``str`` or ``None``
              Dictionary of GW source parameters or JSON file name.

              default: None (uses self.json_file_names["gw_param"])

          **pdet_threshold** : ``float``
              Threshold for detection probability.

              default: 0.5

          **pdet_type** : ``str``
              Detectability condition type.

              Options:

              - 'boolean': Binary detection based on pdet_threshold

              - 'probability_distribution': Uses pdet values directly

              default: 'boolean'

          **output_jsonfile** : ``str`` or ``None``
              JSON file name for storing detectable event parameters.

              default: 'gw_param_detectable.json'

      :Returns:

          **total_rate** : ``float``
              Total GW detection rate (yr^-1).

          **gw_param** : ``dict``
              dictionary of GW source parameters of the detectable events. The included parameters and their units are as follows (for default settings):

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
              | a_1                |              | spin of the primary compact binary         |
              +--------------------+--------------+--------------------------------------+
              | a_2                |              | spin_2 of the compact binary         |
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
              | pdet_L1            |              | pdet of L1                           |
              +--------------------+--------------+--------------------------------------+
              | pdet_H1            |              | pdet of H1                           |
              +--------------------+--------------+--------------------------------------+
              | pdet_V1            |              | pdet of V1                           |
              +--------------------+--------------+--------------------------------------+
              | pdet_net           |              | pdet of the network                  |
              +--------------------+--------------+--------------------------------------+










      .. rubric:: Examples

      >>> from ler.rates import GWRATES
      >>> gwrates = GWRATES()
      >>> gwrates.gw_cbc_statistics()
      >>> total_rate, gw_param = gwrates.gw_rate()



      ..
          !! processed by numpydoc !!

   .. py:method:: rate_function(detectable_size, total_size, verbose=True)

      
      Helper function to calculate the detection rate via Monte Carlo integration.


      :Parameters:

          **detectable_size** : ``int`` or ``float``
              Number of detectable events.

          **total_size** : ``int``
              Total number of simulated events.

          **verbose** : ``bool``
              If True, print rate information.

              default: True

      :Returns:

          **rate** : ``float``
              Detection rate (yr^-1).













      ..
          !! processed by numpydoc !!

   .. py:method:: selecting_n_gw_detectable_events(size=100, batch_size=50000, stopping_criteria=dict(relative_diff_percentage=0.5, number_of_last_batches_to_check=4), pdet_threshold=0.5, resume=True, output_jsonfile='gw_params_n_detectable.json', meta_data_file='meta_gw.json', pdet_type='boolean', trim_to_size=False)

      
      Generate a target number of detectable GW events by iterative batch sampling.

      This function samples GW parameters in batches and saves only the
      detectable events to a JSON file. It optionally stops when the
      cumulative rate has stabilized based on the stopping criteria.

      :Parameters:

          **size** : ``int``
              Target number of detectable samples to collect.

              default: 100

          **batch_size** : ``int``
              Batch size for sampling.

              default: 50000

          **stopping_criteria** : ``dict`` or ``None``
              Criteria for stopping sample collection (but will not stop until n>size).

              Keys:

              - 'relative_diff_percentage': Maximum relative difference in rate (float)

              - 'number_of_last_batches_to_check': Number of batches for comparison (int)

              If None, stops when detectable events exceed size.

              default: dict(relative_diff_percentage=0.5, number_of_last_batches_to_check=4)

          **pdet_threshold** : ``float``
              Threshold for detection probability.

              default: 0.5

          **resume** : ``bool``
              If True, resume from the last batch.

              default: True

          **output_jsonfile** : ``str``
              JSON file name for storing detectable event parameters.

              default: 'gw_params_n_detectable.json'

          **meta_data_file** : ``str``
              JSON file name for storing metadata.

              default: 'meta_gw.json'

          **pdet_type** : ``str``
              Detectability condition type.

              Options:

              - 'boolean': Binary detection based on pdet_threshold

              - 'probability_distribution': Uses pdet values directly

              default: 'boolean'

          **trim_to_size** : ``bool``
              If True, trim final result to exactly the specified size.

              default: False

      :Returns:

          **param_final** : ``dict``
              dictionary of GW source parameters of the detectable events. The included parameters and their units are as follows (for default settings):

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
              | a_1                |              | spin of the primary compact binary   |
              +--------------------+--------------+--------------------------------------+
              | a_2                |              | spin of the secondary compact binary |
              +--------------------+--------------+--------------------------------------+
              | luminosity_distance| Mpc          | luminosity distance                  |
              +--------------------+--------------+--------------------------------------+
              | mass_1_source      | Msun         | mass of the primary compact binary   |
              |                    |              | (source frame)                       |
              +--------------------+--------------+--------------------------------------+
              | mass_2_source      | Msun         | mass of the secondary compact binary |
              |                    |              | (source frame)                       |
              +--------------------+--------------+--------------------------------------+
              | mass_1             | Msun         | mass_1 of the compact binary         |
              |                    |              | (detector frame)                     |
              +--------------------+--------------+--------------------------------------+
              | mass_2             | Msun         | mass_2 of the compact binary         |
              |                    |              | (detector frame)                     |
              +--------------------+--------------+--------------------------------------+
              | pdet_L1            |              | pdet of L1                           |
              +--------------------+--------------+--------------------------------------+
              | pdet_H1            |              | pdet of H1                           |
              +--------------------+--------------+--------------------------------------+
              | pdet_V1            |              | pdet of V1                           |
              +--------------------+--------------+--------------------------------------+
              | pdet_net           |              | pdet of the network                  |
              +--------------------+--------------+--------------------------------------+










      .. rubric:: Examples

      >>> from ler.rates import GWRATES
      >>> gwrates = GWRATES()
      >>> gw_param = gwrates.selecting_n_gw_detectable_events(size=100)



      ..
          !! processed by numpydoc !!


