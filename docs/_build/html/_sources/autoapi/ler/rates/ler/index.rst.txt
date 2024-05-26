:py:mod:`ler.rates.ler`
=======================

.. py:module:: ler.rates.ler

.. autoapi-nested-parse::

   This module contains the main class for calculating the rates of detectable gravitational waves events. The class inherits the :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution` class for source parameters and lens parameters sampling. It also finds the image properties. :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution` inherits the :class:`~ler.gw_source_population.CBCSourceParameterDistribution`, :class:`~ler.image_properties.ImageProperties` and uses `gwsnr` package for SNR calculation.

   ..
       !! processed by numpydoc !!


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   ler.rates.ler.LeR




.. py:class:: LeR(npool=int(4), z_min=0.0, z_max=10.0, event_type='BBH', size=100000, batch_size=50000, cosmology=None, snr_finder=None, pdet_finder=None, list_of_detectors=None, json_file_names=None, interpolator_directory='./interpolator_pickle', ler_directory='./ler_data', verbose=True, **kwargs)


   Bases: :py:obj:`ler.lens_galaxy_population.LensGalaxyParameterDistribution`

   
   Class to calculate both the rates of lensed and unlensed events. Please note that parameters of the simulated events are stored in json file but not as an attribute of the class. This saves RAM memory.


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
           default batch_size = 50000.
           reduce the batch size if you are getting memory error.
           recommended batch_size = 50000, if size = 1000000.

       **snr_finder** : `str`
           default snr_finder = 'gwsnr'.
           if None, the SNR will be calculated using the gwsnr package.
           if 'custom', the SNR will be calculated using a custom function.
           The custom function should have input and output as given in GWSNR.snr method.

       **json_file_names: `dict`**
           names of the json files to strore the necessary parameters.
           default json_file_names = {'ler_params': 'LeR_params.json', 'unlensed_param': 'unlensed_param.json', 'unlensed_param_detectable': 'unlensed_param_detectable.json'}.

       **kwargs** : `keyword arguments`
           Note : kwargs takes input for initializing the :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution`, :meth:`~gwsnr_intialization`.











   .. rubric:: Examples

   >>> from ler.rates import LeR
   >>> ler = LeR()
   >>> unlensed_params = ler.unlensed_cbc_statistics();
   >>> ler.unlensed_rate();

   Instance Attributes
   ----------
   LeR class has the following attributes,

   +-------------------------------------+----------------------------------+
   | Atrributes                          | Type                             |
   +=====================================+==================================+
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
   |:meth:`~store_ler_params`            | Function to store the all the    |
   |                                     | necessary parameters.            |
   +-------------------------------------+----------------------------------+
   |:meth:`~unlensed_cbc_statistics`     | Function to generate unlensed    |
   |                                     | GW source parameters.            |
   +-------------------------------------+----------------------------------+
   |:meth:`~unlensed_sampling_routine`   | Function to generate unlensed    |
   |                                     | GW source parameters.            |
   +-------------------------------------+----------------------------------+
   |:meth:`~unlensed_rate`               | Function to calculate the        |
   |                                     | unlensed rate.                   |
   +-------------------------------------+----------------------------------+
   |:meth:`~selecting_n_unlensed_detectable_events`                         |
   +-------------------------------------+----------------------------------+
   |                                     | Function to select n unlensed    |
   |                                     | detectable events.               |
   +-------------------------------------+----------------------------------+
   |:meth:`~lensed_cbc_statistics`       | Function to generate lensed      |
   |                                     | GW source parameters.            |
   +-------------------------------------+----------------------------------+
   |:meth:`~lensed_sampling_routine`     | Function to generate lensed      |
   |                                     | GW source parameters.            |
   +-------------------------------------+----------------------------------+
   |:meth:`~lensed_rate`                 | Function to calculate the        |
   |                                     | lensed rate.                     |
   +-------------------------------------+----------------------------------+
   |:meth:`~rate_ratio`                  | Function to calculate the rate   |
   |                                     | ratio.                           |
   +-------------------------------------+----------------------------------+
   |:meth:`~rate_comparision_with_rate_calculation                          |
   +-------------------------------------+----------------------------------+
   |                                     | Function to compare the rates    |
   |                                     | calculated using LeR between     |
   |                                     | unlensed and lensed events.      |
   +-------------------------------------+----------------------------------+
   |:meth:`~param_plot`                  | Function to plot the             |
   |                                     | distribution of various          |
   |                                     | parameters.                      |
   +-------------------------------------+----------------------------------+
   |:meth:`~relative_mu_dt_lensed`       | Function to calculate the        |
   |                                     | relative magnification and       |
   |                                     | relative time-delay of lensed    |
   |                                     | events.                          |
   +-------------------------------------+----------------------------------+
   |:meth:`~relative_mu_dt_unlensed`     | Function to calculate the        |
   |                                     | relative magnification and       |
   |                                     | relative time-delay of unlensed  |
   |                                     | events.                          |
   +-------------------------------------+----------------------------------+
   |:meth:`~ mu_vs_dt_plot`              | Function to plot the             |
   |                                     | relative magnification vs        |
   |                                     | relative time-delay.             |
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

   .. py:attribute:: gw_param_sampler_dict

      
      ``dict``

      Dictionary of parameters to initialize the ``CBCSourceParameterDistribution`` class.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: lens_param_sampler_dict

      
      ``dict``

      Dictionary of parameters to initialize the ``LensGalaxyParameterDistribution`` class.















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

   .. py:method:: print_all_params()

      
      Function to print all the parameters.
















      ..
          !! processed by numpydoc !!

   .. py:method:: class_initialization(params=None)

      
      Function to initialize the parent classes. List of relevant initialized instances,

      1. self.sample_source_redshift
      2. self.sample_unlensed_parameters
      3. self.normalization_pdf_z
      4. self.sample_lens_parameters
      5. self.normalization_pdf_z_lensed
      6. self.image_properties
      7. self.get_lensed_snrs

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

   .. py:method:: store_ler_params(output_jsonfile='ler_params.json')

      
      Function to store the all the necessary parameters. This is useful for reproducing the results. All the parameters stored are in string format to make it json compatible.


      :Parameters:

          **output_jsonfile** : `str`
              name of the json file to store the parameters














      ..
          !! processed by numpydoc !!

   .. py:method:: unlensed_cbc_statistics(size=None, resume=False, save_batch=False, output_jsonfile=None)

      
      Function to generate unlensed GW source parameters. This function also stores the parameters in json file.


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
              default output_jsonfile = 'unlensed_params.json'.

      :Returns:

          **unlensed_param** : `dict`
              dictionary of unlensed GW source parameters.
              unlensed_param.keys() = ['zs', 'geocent_time', 'ra', 'dec', 'phase', 'psi', 'theta_jn', 'luminosity_distance', 'mass_1_source', 'mass_2_source', 'mass_1', 'mass_2', 'optimal_snr_net', 'L1', 'H1', 'V1']










      .. rubric:: Examples

      >>> from ler.rates import LeR
      >>> ler = LeR()
      >>> unlensed_param = ler.unlensed_cbc_statistics()



      ..
          !! processed by numpydoc !!

   .. py:method:: unlensed_sampling_routine(size, output_jsonfile, resume=False, save_batch=True)

      
      Function to generate unlensed GW source parameters. This function also stores the parameters in json file.


      :Parameters:

          **size** : `int`
              number of samples.
              default size = 100000.

          **resume** : `bool`
              resume = False (default) or True.
              if True, the function will resume from the last batch.

          **output_jsonfile** : `str`
              json file name for storing the parameters.
              default output_jsonfile = 'unlensed_params.json'.

      :Returns:

          **unlensed_param** : `dict`
              dictionary of unlensed GW source parameters.
              unlensed_param.keys() = ['zs', 'geocent_time', 'ra', 'dec', 'phase', 'psi', 'theta_jn', 'luminosity_distance', 'mass_1_source', 'mass_2_source', 'mass_1', 'mass_2', 'optimal_snr_net', 'L1', 'H1', 'V1']













      ..
          !! processed by numpydoc !!

   .. py:method:: unlensed_rate(unlensed_param=None, snr_threshold=8.0, output_jsonfile=None, detectability_condition='step_function', snr_recalculation=False, snr_threshold_recalculation=5.5)

      
      Function to calculate the unlensed rate. This function also stores the parameters of the detectable events in json file.


      :Parameters:

          **unlensed_param** : `dict` or `str`
              dictionary of GW source parameters or json file name.
              default unlensed_param = 'unlensed_params.json'.

          **snr_threshold** : `float`
              threshold for detection signal to noise ratio.
              e.g. snr_threshold = 8.

          **output_jsonfile** : `str`
              json file name for storing the parameters of the detectable events.
              default output_jsonfile = 'unlensed_params_detectable.json'.

          **detectability_condition** : `str`
              detectability condition.
              default detectability_condition = 'step_function'.
              other options are 'pdet'.

          **snr_recalculation** : `bool`
              if True, the SNR of centain events (snr>snr_threshold_recalculation)will be recalculate with 'inner product'. This is useful when the snr is calculated with 'ann' method.
              default snr_recalculation = False.

          **snr_threshold_recalculation** : `float`
              threshold for recalculation of detection signal to noise ratio.

      :Returns:

          **total_rate** : `float`
              total unlensed rate (Mpc^-3 yr^-1).

          **unlensed_param** : `dict`
              dictionary of unlensed GW source parameters of the detectable events.
              unlensed_param.keys() = ['zs', 'geocent_time', 'ra', 'dec', 'phase', 'psi', 'theta_jn', 'luminosity_distance', 'mass_1_source', 'mass_2_source', 'mass_1', 'mass_2', 'optimal_snr_net', 'L1', 'H1', 'V1']










      .. rubric:: Examples

      >>> from ler.rates import LeR
      >>> ler = LeR()
      >>> ler.unlensed_cbc_statistics();
      >>> total_rate, unlensed_param_detectable = ler.unlensed_rate()



      ..
          !! processed by numpydoc !!

   .. py:method:: lensed_cbc_statistics(size=None, save_batch=False, resume=False, output_jsonfile=None)

      
      Function to generate lensed GW source parameters. This function also stores the parameters in json file.


      :Parameters:

          **size** : `int`
              number of samples.
              default size = 100000.

          **resume** : `bool`
              resume = False (default) or True.
              if True, the function will resume from the last batch.

          **output_jsonfile** : `str`
              json file name for storing the parameters.
              default output_jsonfile = 'lensed_params.json'.

      :Returns:

          **lensed_param** : `dict`
              dictionary of lensed GW source parameters.
              lensed_param.keys() =










      .. rubric:: Examples

      >>> from ler.rates import LeR
      >>> ler = LeR()
      >>> lensed_param = ler.lensed_cbc_statistics()



      ..
          !! processed by numpydoc !!

   .. py:method:: lensed_sampling_routine(size, output_jsonfile, save_batch=True, resume=False)

      
      Function to generate lensed GW source parameters. This function also stores the parameters in json file.


      :Parameters:

          **size** : `int`
              number of samples.
              default size = 100000.

          **resume** : `bool`
              resume = False (default) or True.
              if True, the function will resume from the last batch.

          **output_jsonfile** : `str`
              json file name for storing the parameters.
              default output_jsonfile = 'lensed_params.json'.

      :Returns:

          **lensed_param** : `dict`
              dictionary of lensed GW source parameters.
              lensed_param.keys() =













      ..
          !! processed by numpydoc !!

   .. py:method:: lensed_rate(lensed_param=None, snr_threshold=[8.0, 8.0], num_img=[1, 1], output_jsonfile=None, nan_to_num=True, detectability_condition='step_function', snr_recalculation=False, snr_threshold_recalculation=[5.5, 5.5])

      
      Function to calculate the lensed rate. This function also stores the parameters of the detectable events in json file.


      :Parameters:

          **lensed_param** : `dict` or `str`
              dictionary of GW source parameters or json file name.
              default lensed_param = 'lensed_params.json'.

          **snr_threshold** : `float`
              threshold for detection signal to noise ratio.
              default snr_threshold = [8.0,8.0].

          **num_img** : `int`
              number of images.
              default num_img = [1,1].

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
              if True, the SNR of centain events (snr>snr_threshold_recalculation)will be recalculate with 'inner product'. This is useful when the snr is calculated with 'ann' method.
              default snr_recalculation = False.

          **snr_threshold_recalculation** : `float`
              threshold for recalculation of detection signal to noise ratio.

      :Returns:

          **total_rate** : `float`
              total lensed rate (Mpc^-3 yr^-1).

          **lensed_param** : `dict`
              dictionary of lensed GW source parameters of the detectable events.
              lensed_param.keys() =










      .. rubric:: Examples

      >>> from ler.rates import LeR
      >>> ler = LeR()
      >>> ler.lensed_cbc_statistics();
      >>> total_rate, lensed_param_detectable = ler.lensed_rate()



      ..
          !! processed by numpydoc !!

   .. py:method:: rate_ratio()

      
      Function to calculate and display unlensed and lensed merger rate ratio.
      It will get the unlensed_rate and lensed_rate from self.json_file_ler_param



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

   .. py:method:: rate_comparision_with_rate_calculation(unlensed_param=None, snr_threshold_unlensed=8.0, output_jsonfile_unlensed=None, lensed_param=None, snr_threshold_lensed=[8.0, 8.0], num_img=[1, 1], output_jsonfile_lensed=None, nan_to_num=True, detectability_condition='step_function')

      
      Function to calculate the unlensed and lensed rate and compare by computing the ratio. This function also stores the parameters of the detectable events in json file.


      :Parameters:

          **unlensed_param** : `dict` or `str`
              dictionary of GW source parameters or json file name.
              default unlensed_param = 'unlensed_params.json'.

          **snr_threshold_unlensed** : `float`
              threshold for detection signal to noise ratio.
              e.g. snr_threshold = 8.

          **output_jsonfile_unlensed** : `str`
              json file name for storing the parameters of the detectable events.
              default output_jsonfile = 'unlensed_params_detectable.json'.

          **lensed_param** : `dict` or `str`
              dictionary of GW source parameters or json file name.
              default lensed_param = 'lensed_params.json'.

          **snr_threshold_lensed** : `float`
              threshold for detection signal to noise ratio.
              e.g. snr_threshold = 8.

          **output_jsonfile_lensed** : `str`
              json file name for storing the parameters of the detectable events.
              default output_jsonfile = 'lensed_params_detectable.json'.

          **detectability_condition** : `str`
              detectability condition.
              default detectability_condition = 'step_function'.
              other options are 'pdet'.

      :Returns:

          **rate_ratio** : `float`
              rate ratio.

          **unlensed_param** : `dict`
              dictionary of unlensed GW source parameters of the detectable events.

          **lensed_param** : `dict`
              dictionary of lensed GW source parameters of the detectable events.










      .. rubric:: Examples

      >>> from ler.rates import LeR
      >>> ler = LeR()
      >>> ler.unlensed_cbc_statistics();
      >>> ler.lensed_cbc_statistics();
      >>> rate_ratio, unlensed_param, lensed_param = ler.rate_comparision_with_rate_calculation()



      ..
          !! processed by numpydoc !!

   .. py:method:: selecting_n_unlensed_detectable_events(size=100, batch_size=None, snr_threshold=8.0, resume=False, output_jsonfile='n_unlensed_param_detectable.json', meta_data_file='meta_unlensed.json', trim_to_size=True, snr_recalculation=False, snr_threshold_recalculation=5.5)

      
      Function to select n unlensed detectable events.


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
              default output_jsonfile = 'n_unlensed_params_detectable.json'.

      :Returns:

          **param_final** : `dict`
              dictionary of unlensed GW source parameters of the detectable events.
              param_final.keys() = ['zs', 'geocent_time', 'ra', 'dec', 'phase', 'psi', 'theta_jn', 'luminosity_distance', 'mass_1_source', 'mass_2_source', 'mass_1', 'mass_2', 'optimal_snr_net', 'L1', 'H1', 'V1']










      .. rubric:: Examples

      >>> from ler.rates import LeR
      >>> ler = LeR()
      >>> unlensed_param_final = ler.selecting_n_unlensed_detectable_events(size=500)



      ..
          !! processed by numpydoc !!

   .. py:method:: selecting_n_lensed_detectable_events(size=100, batch_size=None, snr_threshold=[8.0, 8.0], num_img=[1, 1], resume=False, detectability_condition='step_function', output_jsonfile='n_lensed_params_detectable.json', meta_data_file='meta_lensed.json', trim_to_size=True, nan_to_num=False, snr_recalculation=False, snr_threshold_recalculation=[5.5, 5.5])

      
      Function to select n lensed detectable events.


      :Parameters:

          **size** : `int`
              number of samples to be selected.
              default size = 100.

          **snr_threshold** : `float`
              threshold for detection signal to noise ratio.
              e.g. snr_threshold = 8.

          **num_img** : `int`
              number of images.
              default num_img = 2.

          **resume** : `bool`
              if True, it will resume the sampling from the last batch.
              default resume = False.

          **detectability_condition** : `str`
              detectability condition.
              default detectability_condition = 'step_function'.
              other options are 'pdet'.

          **output_jsonfile** : `str`
              json file name for storing the parameters.

      :Returns:

          **param_final** : `dict`
              dictionary of lensed GW source parameters of the detectable events.
              param_final.keys() =










      .. rubric:: Examples

      >>> from ler.rates import LeR
      >>> ler = LeR()
      >>> lensed_param_final = ler.selecting_n_lensed_detectable_events(size=500)



      ..
          !! processed by numpydoc !!


