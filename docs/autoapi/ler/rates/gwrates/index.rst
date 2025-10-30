:py:mod:`ler.rates.gwrates`
===========================

.. py:module:: ler.rates.gwrates

.. autoapi-nested-parse::

   This module contains the main class for calculating the rates of detectable gravitational waves events. The class inherits the :class:`~ler.gw_source_population.CBCSourceParameterDistribution` class for source parameters sampling and uses `gwsnr` package for SNR calculation.

   ..
       !! processed by numpydoc !!


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   ler.rates.gwrates.GWRATES




.. py:class:: GWRATES(npool=int(4), z_min=0.0, z_max=10.0, event_type='BBH', size=100000, batch_size=50000, cosmology=None, snr_finder=None, pdet_finder=None, list_of_detectors=None, json_file_names=None, interpolator_directory='./interpolator_pickle', create_new_interpolator=False, ler_directory='./ler_data', verbose=True, **kwargs)


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
   LeR class has the following attributes:

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
   LeR class has the following methods:

   +-------------------------------------+----------------------------------+
   | Methods                             | Description                      |
   +=====================================+==================================+
   |:meth:`~class_initialization`        | Function to initialize the       |
   |                                     | parent classes                   |
   +-------------------------------------+----------------------------------+
   |:meth:`~gwsnr_initialization`         | Function to initialize the       |
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

   .. py:method:: gw_rate(gw_param=None, snr_threshold=8.0, pdet_threshold=0.5, output_jsonfile=None, detectability_condition='step_function', snr_recalculation=False, snr_threshold_recalculation=[4, 20])

      
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

   .. py:method:: selecting_n_gw_detectable_events(size=100, batch_size=None, stopping_criteria=dict(relative_diff_percentage=0.5, number_of_last_batches_to_check=4), snr_threshold=8.0, pdet_threshold=0.5, resume=False, output_jsonfile='gw_params_n_detectable.json', meta_data_file='meta_gw.json', detectability_condition='step_function', trim_to_size=False, snr_recalculation=False, snr_threshold_recalculation=[4, 12])

      
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


