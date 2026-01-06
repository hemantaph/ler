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
   ler.rates.load_json
   ler.rates.append_json
   ler.rates.get_param_from_json
   ler.rates.batch_handler



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

           default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

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
   | :meth:`~list_of_detectors`                     | ``list``         |       | List of detector names                         |
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

              default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)













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

   .. py:method:: lensed_cbc_statistics(size=100000, batch_size=50000, save_batch=True, resume=True, output_jsonfile=None)

      
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
              | theta_E                      | arcsec    | Einstein radius                                       |
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
              | mass_1_source                | Msun      | mass of the primary compact binary (source frame)     |
              +------------------------------+-----------+-------------------------------------------------------+
              | mass_2_source                | Msun      | mass of the secondary compact binary (source frame)   |
              +------------------------------+-----------+-------------------------------------------------------+
              | mass_1                       | Msun      | mass of the primary compact binary (detector frame)   |
              +------------------------------+-----------+-------------------------------------------------------+
              | mass_2                       | Msun      | mass of the secondary compact binary (detector frame) |
              +------------------------------+-----------+-------------------------------------------------------+
              | x0_image_positions           | arcsec    | x-coordinate (RA-like axis) of the images             |
              +------------------------------+-----------+-------------------------------------------------------+
              | x1_image_positions           | arcsec    | y-coordinate (Dec-like axis) of the images            |
              +------------------------------+-----------+-------------------------------------------------------+
              | magnifications               |           | magnifications                                        |
              +------------------------------+-----------+-------------------------------------------------------+
              | time_delays                  |           | time delays                                           |
              +------------------------------+-----------+-------------------------------------------------------+
              | image_type                   |           | image type                                            |
              +------------------------------+-----------+-------------------------------------------------------+
              | n_images                     |           | number of images                                      |
              +------------------------------+-----------+-------------------------------------------------------+
              | x_source                     | arcsec    | x-coordinate (RA-like axis) of the source             |
              +------------------------------+-----------+-------------------------------------------------------+
              | y_source                     | arcsec    | y-coordinate (Dec-like axis) of the source            |
              +------------------------------+-----------+-------------------------------------------------------+
              | effective_luminosity_distance| Mpc       | effective luminosity distance of the images           |
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
              | theta_E                      | arcsec    | Einstein radius                                       |
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
              | mass_1_source                | Msun      | mass of the primary compact binary (source frame)     |
              +------------------------------+-----------+-------------------------------------------------------+
              | mass_2_source                | Msun      | mass of the secondary compact binary (source frame)   |
              +------------------------------+-----------+-------------------------------------------------------+
              | mass_1                       | Msun      | mass of the primary compact binary (detector frame)   |
              +------------------------------+-----------+-------------------------------------------------------+
              | mass_2                       | Msun      | mass of the secondary compact binary (detector frame) |
              +------------------------------+-----------+-------------------------------------------------------+
              | x0_image_positions           | arcsec    | x-coordinate (RA-like axis) of the images             |
              +------------------------------+-----------+-------------------------------------------------------+
              | x1_image_positions           | arcsec    | y-coordinate (Dec-like axis) of the images            |
              +------------------------------+-----------+-------------------------------------------------------+
              | magnifications               |           | magnifications                                        |
              +------------------------------+-----------+-------------------------------------------------------+
              | time_delays                  |           | time delays                                           |
              +------------------------------+-----------+-------------------------------------------------------+
              | image_type                   |           | image type                                            |
              +------------------------------+-----------+-------------------------------------------------------+
              | n_images                     |           | number of images                                      |
              +------------------------------+-----------+-------------------------------------------------------+
              | x_source                     | arcsec    | x-coordinate (RA-like axis) of the source             |
              +------------------------------+-----------+-------------------------------------------------------+
              | y_source                     | arcsec    | y-coordinate (Dec-like axis) of the source            |
              +------------------------------+-----------+-------------------------------------------------------+
              | effective_luminosity_distance| Mpc       | effective luminosity distance of the images           |
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

          **unlensed_param** : ``dict``
              Dictionary of detectable unlensed GW source parameters.

          **lensed_param** : ``dict``
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

   .. py:method:: selecting_n_lensed_detectable_events(size=100, stopping_criteria=dict(relative_diff_percentage=2, number_of_last_batches_to_check=4), batch_size=50000, pdet_threshold=[0.5, 0.5], num_img=[1, 1], resume=True, pdet_type='boolean', output_jsonfile='n_lensed_params_detectable.json', meta_data_file='meta_lensed.json', trim_to_size=False, nan_to_num=False)

      
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

           default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

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
   | :meth:`~binary_masses_BBH_popI_II_powerlaw_gaussian`| Sample BBH masses with PowerLaw+PEAK model     |
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

   .. py:method:: binary_masses_BBH_popI_II_powerlaw_gaussian(size, get_attribute=False, **kwargs)

      
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
      >>> m1_src, m2_src = cbc.binary_masses_BBH_popI_II_powerlaw_gaussian(size=1000)



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

.. py:class:: GWRATES(npool=int(4), z_min=0.0, z_max=10.0, event_type='BBH', size=100000, batch_size=50000, cosmology=None, snr_finder=None, pdet_finder=None, list_of_detectors=None, json_file_names=None, interpolator_directory='./interpolator_json', create_new_interpolator=False, ler_directory='./ler_data', verbose=True, **kwargs)


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
           where optimal_snr_dict.keys = ['snr_net']. Refer to `gwsnr` package's GWSNR.snr attribute for more details.

       **pdet_finder** : `function`
           default pdet_finder = None.
           The rate calculation uses either the pdet_finder or the snr_finder to calculate the detectable events. The custom pdet finder function should follow the following signature:
           def pdet_finder(gw_param_dict):
               ...
               return pdet_net_dict
           where pdet_net_dict.keys = ['pdet_net']. For example uses, refer to [GRB pdet example](https://ler.readthedocs.io/en/latest/examples/rates/grb%20detection%20rate.html).

       **json_file_names: `dict`**
           names of the json files to strore the necessary parameters.
           default json_file_names = {'gwrates_params':'gwrates_params.json', 'gw_param': 'gw_param.json', 'gw_param_detectable': 'gw_param_detectable.json'}.

       **interpolator_directory** : `str`
           directory to store the interpolators.
           default interpolator_directory = './interpolator_json'. This is used for storing the various interpolators related to `ler` and `gwsnr` package.

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
   LeR class has the following attributes:

   +-------------------------+----------------------+
   | Atrributes                          | Type                             |
   +=====================================+==================================+
   |:attr:`~npool`                       | `int`                            |
   +-------------------------+----------------------+
   |:attr:`~z_min`                       | `float`                          |
   +-------------------------+----------------------+
   |:attr:`~z_max`                       | `float`                          |
   +-------------------------+----------------------+
   |:attr:`~event_type`                  | `str`                            |
   +-------------------------+----------------------+
   |:attr:`~cosmo`                       | `astropy.cosmology`              |
   +-------------------------+----------------------+
   |:attr:`~size`                        | `int`                            |
   +-------------------------+----------------------+
   |:attr:`~batch_size`                  | `int`                            |
   +-------------------------+----------------------+
   |:attr:`~json_file_names`             | `dict`                           |
   +-------------------------+----------------------+
   |:attr:`~interpolator_directory`      | `str`                            |
   +-------------------------+----------------------+
   |:attr:`~ler_directory`               | `str`                            |
   +-------------------------+----------------------+
   |:attr:`~gwsnr`                       | `bool`                           |
   +-------------------------+----------------------+
   |:attr:`~gw_param_sampler_dict`       | `dict`                           |
   +-------------------------+----------------------+
   |:attr:`~snr_calculator_dict`         | `dict`                           |
   +-------------------------+----------------------+
   |:attr:`~gw_param`                    | `dict`                           |
   +-------------------------+----------------------+
   |:attr:`~gw_param_detectable`         | `dict`                           |
   +-------------------------+----------------------+

   Instance Methods
   ----------
   LeR class has the following methods:

   +-------------------------+----------------------+
   | Methods                             | Description                      |
   +=====================================+==================================+
   |:meth:`~class_initialization`        | Function to initialize the       |
   |                                     | parent classes                   |
   +-------------------------+----------------------+
   |:meth:`~gwsnr_initialization`         | Function to initialize the       |
   |                                     | gwsnr class                      |
   +-------------------------+----------------------+
   |:meth:`~snr`                         | Function to get the snr with the |
   |                                     | given parameters.                |
   +-------------------------+----------------------+
   |:meth:`~store_gwrates_params`        | Function to store the all the    |
   |                                     | necessary parameters.            |
   +-------------------------+----------------------+
   |:meth:`~gw_cbc_statistics`           | Function to generate gw          |
   |                                     | GW source parameters.            |
   +-------------------------+----------------------+
   |:meth:`~gw_sampling_routine`         | Function to generate gw          |
   |                                     | GW source parameters.            |
   +-------------------------+----------------------+
   |:meth:`~gw_rate`                     | Function to calculate the        |
   |                                     | gw rate.                         |
   +-------------------------+----------------------+
   |:meth:`~selecting_n_gw_detectable_events`                               |
   +-------------------------+----------------------+
   |                                     | Function to select n gw    |
   |                                     | detectable events.               |
   +-------------------------+----------------------+
   |:meth:`~gw_param_plot`               | Function to plot the             |
   |                                     | distribution of the GW source    |
   |                                     | parameters.                      |
   +-------------------------+----------------------+



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

   .. py:attribute:: cosmo
      :value: 'None'

      
      ``astropy.cosmology``

      Cosmology to use for the calculation.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: size
      :value: 'None'

      
      ``int``

      Number of samples for sampling.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: batch_size
      :value: 'None'

      
      ``int``

      Batch size for sampling.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: json_file_names
      :value: 'None'

      
      ``dict``

      Names of the json files to store the necessary parameters.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: interpolator_directory
      :value: 'None'

      
      ``str``

      Directory to store the interpolators.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: ler_directory
      :value: 'None'

      
      ``str``

      Directory to store the parameters.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: gwsnr
      :value: 'None'

      
      ``bool``

      If True, the SNR will be calculated using the gwsnr package.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: gw_param_sampler_dict
      :value: 'None'

      
      ``dict``

      Dictionary of parameters to initialize the ``CBCSourceParameterDistribution`` class.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: snr_calculator_dict
      :value: 'None'

      
      ``dict``

      Dictionary of parameters to initialize the ``GWSNR`` class.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: npool
      :value: '4'

      
      Number of processors for multiprocessing.



      :Returns:

          **npool** : ``int``
              Number of parallel processes to use.

              default: 4













      ..
          !! processed by numpydoc !!

   .. py:method:: print_all_params_ler()

      
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

   .. py:method:: gwsnr_initialization(params=None)

      
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

   .. py:method:: gw_rate(gw_param=None, snr_threshold=10.0, pdet_threshold=0.5, output_jsonfile=None, detectability_condition='step_function', snr_recalculation=False, snr_threshold_recalculation=[4, 20])

      
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

          **snr_threshold_recalculation** : `list`
              lower and upper threshold for recalculation of detection signal to noise ratio.
              default snr_threshold_recalculation = [4, 20].

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

   .. py:method:: selecting_n_gw_detectable_events(size=100, batch_size=None, stopping_criteria=dict(relative_diff_percentage=0.5, number_of_last_batches_to_check=4), snr_threshold=10.0, pdet_threshold=0.5, resume=False, output_jsonfile='gw_params_n_detectable.json', meta_data_file='meta_gw.json', detectability_condition='step_function', trim_to_size=False, snr_recalculation=False, snr_threshold_recalculation=[4, 12])

      
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

          **snr_threshold_recalculation** : `list`
              lower and upper threshold for recalculation of detection signal to noise ratio.
              default snr_threshold_recalculation = [4, 12].

      :Returns:

          **param_final** : `dict`
              dictionary of gw GW source parameters of the detectable events. Refer to :attr:`~gw_param` for details.










      .. rubric:: Examples

      >>> from ler.rates import LeR
      >>> ler = LeR()
      >>> gw_param = ler.selecting_n_gw_detectable_events(size=100)



      ..
          !! processed by numpydoc !!


