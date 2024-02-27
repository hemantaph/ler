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




.. py:class:: GWRATES(npool=int(4), z_min=0.0, z_max=10.0, event_type='BBH', size=100000, batch_size=25000, cosmology=None, snr_finder='gwsnr', json_file_names=None, interpolator_directory='./interpolator_pickle', ler_directory='./ler_data', verbose=True, **kwargs)


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
           for popI_II, popIII, primordial, BNS z_max = 10., 40., 40., 2. respectively.

       **size** : `int`
           number of samples for sampling.
           default size = 100000.

       **batch_size** : `int`
           batch size for SNR calculation.
           default batch_size = 25000.
           reduce the batch size if you are getting memory error.
           recommended batch_size = 50000, if size = 1000000.

       **snr_finder** : `str`
           default snr_finder = 'gwsnr'.
           if 'gwsnr', the SNR will be calculated using the gwsnr package.
           if 'custom', the SNR will be calculated using a custom function.
           The custom function should have input and output as given in GWSNR.snr method.

       **json_file_names: `dict`**
           names of the json files to strore the necessary parameters.
           default json_file_names = {'ler_param': 'LeR_params.json', 'gw_param': 'gw_param.json', 'gw_param_detectable': 'gw_param_detectable.json'}.

       **kwargs** : `keyword arguments`
           Note : kwargs takes input for initializing the :class:`~ler.gw_source_population.CBCSourceParameterDistribution`, :meth:`~gwsnr_intialization`.











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
   |:attr:`~directory`                   | `str`                            |
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

      Names of the json files to strore the necessary parameters.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: directory

      
      ``str``

      Directory to store the interpolators.















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

   .. py:method:: print_all_params()

      
      Function to print all the parameters.
















      ..
          !! processed by numpydoc !!

   .. py:method:: class_initialization(params=None)

      
      Function to initialize the parent classes. List of relevant initialized instances,

      1. self.sample_source_redshift
      2. self.sample_gw_parameters
      3. self.normalization_pdf_z

      :Parameters:

          **params** : `dict`
              dictionary of parameters to initialize the parent classes














      ..
          !! processed by numpydoc !!

   .. py:method:: gwsnr_intialization(params=None)

      
      Function to initialize the gwsnr class


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

   .. py:method:: gw_cbc_statistics(size=None, resume=False, save_batch=True, output_jsonfile=None)

      
      Function to generate gw GW source parameters. This function also stores the parameters in json file.


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
              default output_jsonfile = 'gw_params.json'.

      :Returns:

          **gw_param** : `dict`
              dictionary of gw GW source parameters.
              gw_param.keys() = ['zs', 'geocent_time', 'ra', 'dec', 'phase', 'psi', 'theta_jn', 'luminosity_distance', 'mass_1_source', 'mass_2_source', 'mass_1', 'mass_2', 'optimal_snr_net', 'L1', 'H1', 'V1']










      .. rubric:: Examples

      >>> from ler.rates import GWRATES
      >>> ler = GWRATES()
      >>> param = ler.gw_cbc_statistics()



      ..
          !! processed by numpydoc !!

   .. py:method:: gw_sampling_routine(size, output_jsonfile, resume=False, save_batch=True)

      
      Function to generate gw GW source parameters. This function also stores the parameters in json file.


      :Parameters:

          **size** : `int`
              number of samples.
              default size = 100000.

          **resume** : `bool`
              resume = False (default) or True.
              if True, the function will resume from the last batch.

          **output_jsonfile** : `str`
              json file name for storing the parameters.
              default output_jsonfile = 'gw_params.json'.

      :Returns:

          **gw_param** : `dict`
              dictionary of gw GW source parameters.
              gw_param.keys() = ['zs', 'geocent_time', 'ra', 'dec', 'phase', 'psi', 'theta_jn', 'luminosity_distance', 'mass_1_source', 'mass_2_source', 'mass_1', 'mass_2', 'optimal_snr_net', 'L1', 'H1', 'V1']













      ..
          !! processed by numpydoc !!

   .. py:method:: gw_rate(gw_param=None, snr_threshold=8.0, output_jsonfile=None, detectability_condition='step_function', snr_recalculation=False, threshold_snr_recalculation=6.0)

      
      Function to calculate the gw rate. This function also stores the parameters of the detectable events in json file.


      :Parameters:

          **gw_param** : `dict` or `str`
              dictionary of GW source parameters or json file name.
              default gw_param = self.json_file_names["gw_param"]

          **snr_threshold** : `float`
              threshold for detection signal to noise ratio.
              e.g. snr_threshold = 8.

          **output_jsonfile** : `str`
              json file name for storing the parameters of the detectable events.
              default output_jsonfile = 'gw_params_detectable.json'.

          **detectability_condition** : `str`
              detectability condition.
              default detectability_condition = 'step_function'.
              other options are 'pdet'.

          **snr_recalculation** : `bool`
              if True, the SNR of centain events (snr>threshold_snr_recalculation)will be recalculate with 'inner product'. This is useful when the snr is calculated with 'ann' method.
              default snr_recalculation = False.

          **threshold_snr_recalculation** : `float`
              threshold for recalculation of detection signal to noise ratio.

      :Returns:

          **total_rate** : `float`
              total gw rate (Mpc^-3 yr^-1).

          **gw_param** : `dict`
              dictionary of gw GW source parameters of the detectable events.
              gw_param.keys() = ['zs', 'geocent_time', 'ra', 'dec', 'phase', 'psi', 'theta_jn', 'luminosity_distance', 'mass_1_source', 'mass_2_source', 'mass_1', 'mass_2', 'optimal_snr_net', 'L1', 'H1', 'V1']










      .. rubric:: Examples

      >>> from ler.rates import GWRATES
      >>> ler = GWRATES()
      >>> total_rate, gw_param = ler.gw_rate()



      ..
          !! processed by numpydoc !!

   .. py:method:: selecting_n_gw_detectable_events(size=100, batch_size=None, snr_threshold=8.0, resume=False, output_jsonfile='gw_params_n_detectable.json', meta_data_file='meta_gw.json', trim_to_size=True)

      
      Function to select n gw detectable events.


      :Parameters:

          **size** : `int`
              number of samples to be selected.
              default size = 100.

          **snr_threshold** : `float`
              threshold for detection signal to noise ratio.
              e.g. snr_threshold = 8.

          **resume** : `bool`
              if True, it will resume the sampling from the last batch.
              default resume = False.

          **output_jsonfile** : `str`
              json file name for storing the parameters.
              default output_jsonfile = 'gw_params_detectable.json'.

      :Returns:

          **param_final** : `dict`
              dictionary of gw GW source parameters of the detectable events.
              param_final.keys() = ['zs', 'geocent_time', 'ra', 'dec', 'phase', 'psi', 'theta_jn', 'luminosity_distance', 'mass_1_source', 'mass_2_source', 'mass_1', 'mass_2', 'optimal_snr_net', 'L1', 'H1', 'V1']










      .. rubric:: Examples

      >>> from ler.rates import GWRATES
      >>> ler = GWRATES()
      >>> param_final = ler.selecting_n_gw_detectable_events(size=500)



      ..
          !! processed by numpydoc !!


