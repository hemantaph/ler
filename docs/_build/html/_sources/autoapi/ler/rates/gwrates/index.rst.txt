:py:mod:`ler.rates.gwrates`
===========================

.. py:module:: ler.rates.gwrates

.. autoapi-nested-parse::

   Module for calculating detection rates of gravitational wave events.

   This module contains the main ``GWRATES`` class for calculating the rates of
   detectable gravitational wave events. The class inherits from
   :class:`~ler.gw_source_population.CBCSourceParameterDistribution` for source
   parameters sampling and utilizes the ``gwsnr`` package for detection probability
   calculation.

   Inheritance hierarchy:

   - :class:`~ler.gw_source_population.CBCSourceParameterDistribution`

     - :class:`~ler.gw_source_population.CBCSourceRedshiftDistribution`

   - Uses the ``gwsnr`` package for pdet calculation.

   Usage:
       Basic workflow for rate calculation:

       >>> from ler.rates import GWRATES
       >>> gwrates = GWRATES()
       >>> gw_params = gwrates.gw_cbc_statistics()
       >>> gwrates.gw_rate()

   Copyright (C) 2026 Phurailatpam Hemantakumar. Distributed under MIT License.

   ..
       !! processed by numpydoc !!


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   ler.rates.gwrates.GWRATES




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


