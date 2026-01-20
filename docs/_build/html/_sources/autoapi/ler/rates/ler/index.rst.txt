:py:mod:`ler.rates.ler`
=======================

.. py:module:: ler.rates.ler

.. autoapi-nested-parse::

   Module for calculating detection rates of gravitational wave events.

   This module contains the main ``LeR`` class for calculating the rates of
   detectable gravitational wave events, both lensed and unlensed. The class
   inherits from :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution`
   for source and lens parameters sampling, and utilizes image property calculations.

   The inheritance hierarchy is as follows:

   - :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution`

     - :class:`~ler.lens_galaxy_population.OpticalDepth`

     - :class:`~ler.image_properties.ImageProperties`

     - :class:`~ler.gw_source_population.CBCSourceParameterDistribution`

     - :class:`~ler.gw_source_population.CBCSourceRedshiftDistribution`

   - Uses the ``gwsnr`` package for pdet calculation.

   Usage:
       Basic workflow for rate calculation:

       >>> from ler.rates import LeR
       >>> ler = LeR()
       >>> unlensed_params = ler.unlensed_cbc_statistics()
       >>> ler.unlensed_rate()
       >>> lensed_params = ler.lensed_cbc_statistics()
       >>> ler.lensed_rate()
       >>> ler.rate_ratio()

   Copyright (C) 2026 Phurailatpam Hemantakumar. Distributed under MIT License.

   ..
       !! processed by numpydoc !!


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   ler.rates.ler.LeR




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


