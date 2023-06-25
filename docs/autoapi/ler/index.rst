:orphan:

:py:mod:`ler`
=============

.. py:module:: ler


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   helperroutines/index.rst
   lens_galaxy_population/index.rst
   ler/index.rst
   multiprocessing_routine/index.rst
   source_population/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   ler.LeR
   ler.LensGalaxyPopulation
   ler.SourceGalaxyPopulationModel
   ler.CompactBinaryPopulation



Functions
~~~~~~~~~

.. autoapisummary::

   ler.solve_lens_equation1
   ler.solve_lens_equation2
   ler.add_dictionaries_together
   ler.rejection_sample



.. py:class:: LeR(nsamples=100000, npool=int(4), z_min=0.0, z_max=10.0, batch_size=25000, snr_finder='gwsnr', json_file_ler_param='./LeR_params.json', **kwargs)

   
   Class to calculate both the rates of lensed and unlensed events.


   :Parameters:

       **nsamples** : `int`
           number of samples for sampling.
           default nsamples = 100000.

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

       **batch_size** : `int`
           batch size for SNR calculation.
           default batch_size = 25000.
           reduce the batch size if you are getting memory error.

       **snr_finder** : `str`
           default snr_finder = 'gwsnr'.
           if 'gwsnr', the SNR will be calculated using the gwsnr package.
           if 'custom', the SNR will be calculated using a custom function.

       **kwargs** : `keyword arguments`
           Note : kwargs takes input for initializing the :class:`~ler.CompactBinaryPopulation`, :class:`LensGalaxyPopulation`, :meth:`~gwsnr_intialization`.











   .. rubric:: Examples

   - class initialization
   - ``ler`` needs `gwsnr <https://github.com/hemantaph/gwsnr/>`_.
   - generation of ``gwsnr`` snr interpolator will take time at the first initialization. The interpolator will be stored in the working dir.
   - ``m_min``, ``m_max`` were used for initializing the ``CompactBinaryPopulation`` class. ``waveform_approximant`` was used for initializing the ``snr_calculator`` (``gwsnr``) class. ``min_lensed_images`` was used for initializing the ``LensGalaxyPopulation`` class.

   >>> from ler import LeR
   >>> ler_ = LeR(nsamples=100000, npool=int(4), z_min=0., z_max=10., batch_size=25000, snr_finder='gwsnr', m_min=4.59, m_max=86.22, waveform_approximant='IMRPhenomD', min_lensed_images=2)
   Given: IMR waveform
   psds not given. Choosing bilby's default psds
   getting stored interpolator...
   In case if you need regeneration of interpolator of the given gwsnr param, please delete this file, ./interpolator_pickle/halfSNR_dict_0.pickle

   Instance Attributes
   ----------
   LeR class has the following attributes,

   +-------------------------------------+----------------------------------+
   | Atrributes                          | Type                             |
   +=====================================+==================================+
   |:attr:`~gw_param`                    |`dict`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~gw_param_detectable`         |`dict`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~lensed_param`                |`dict`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~lensed_param_detectable`     |`dict`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~gw_param_sampler_dict`       |`dict`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~lensed_param_sampler_dict`   |`dict`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~snr_calculator_dict`         |`dict`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~z_to_Dc`                     |`scipy.interpolate.interp1d`      |
   +-------------------------------------+----------------------------------+
   |:attr:`~Dc_to_z`                     |`scipy.interpolate.interp1d`      |
   +-------------------------------------+----------------------------------+
   |:attr:`~z_to_luminosity_distance`    |`scipy.interpolate.interp1d`      |
   +-------------------------------------+----------------------------------+
   |:attr:`~differential_comoving_volume`|`scipy.interpolate.interp1d`      |
   +-------------------------------------+----------------------------------+
   |:attr:`~compact_binary_pop`          |`CompactBinaryPopulation class`   |
   +-------------------------------------+----------------------------------+
   |:attr:`~lens_galaxy_pop`             |`LensGalaxyPopulation class`      |
   +-------------------------------------+----------------------------------+
   | :attr:`~snr`                        |``gwsnr`` `package`               |
   +-------------------------------------+----------------------------------+

   Instance Methods
   ----------
   LeR class has the following method(s),

   +------------------------------------+-------------------------------------+
   | Method(s)                          | Description                         |
   +====================================+=====================================+
   |:meth:`~gwsnr_intialization`        |Function for initializing the        |
   |                                    |``gwsnr`` package.                   |
   +------------------------------------+-------------------------------------+
   |:meth:`~create_lookup_tables`       |To creating lookup tables for fast   |
   |                                    |calculation for the following        |
   |                                    |conversion operations,               |
   |                                    |redshift to co-moving distance.      |
   |                                    |co-moving distance to redshift.      |
   |                                    |redshift to luminosity distance.     |
   +------------------------------------+-------------------------------------+
   |:meth:`~unlensed_cbc_statistics`    |Function to generate unlensed GW     |
   |                                    |source parameters.                   |
   +------------------------------------+-------------------------------------+
   |:meth:`~unlensed_rate`              |Function to calculate unlensed       |
   |                                    |merger rate.                         |
   +------------------------------------+-------------------------------------+
   |:meth:`~lensed_cbc_statistics`      |Function to generate lensed GW       |
   |                                    |source parameters.                   |
   +------------------------------------+-------------------------------------+
   |:meth:`~lensed_rate`                |Function to calculate lensed         |
   |                                    |merger rate.                         |
   +------------------------------------+-------------------------------------+
   |:meth:`~batch_handler`              |Function to handle the batch size.   |
   +------------------------------------+-------------------------------------+
   |:meth:`~store_ler_params`           |Fuction to store the parameters of   |
   |                                    |the LER model.                       |
   +------------------------------------+-------------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:property:: gw_param

      
      ``bool``, ``dict``

      gw_param is a dictionary of unlensed parameters (source parameters)

      it will be populated when unlened_cbc_statistics() is called

      if unavailable, the unlensed parameters will be sampled when unlensed_rate() is called

      gw_param.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time']















      ..
          !! processed by numpydoc !!

   .. py:property:: gw_param_detectable

      
      ``bool``, ``dict``

      gw_param_detectable is a dictionary of unlensed parameters (source parameters)

      it will be populated when unlened_cbc_statistics() is called

      if unavailable, the unlensed parameters will be sampled when unlensed_rate() is called

      gw_param_detectable.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time']















      ..
          !! processed by numpydoc !!

   .. py:property:: lensed_param

      
      ``bool``, ``dict``

      lensed_param is a dictionary of lensed parameters

      it will be populated when lensed_cbc_statistics() is called

      if unavailable, the lensed parameters will be sampled when lensed_rate() is called

      lensed_param.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time', 'lensed_images']















      ..
          !! processed by numpydoc !!

   .. py:property:: lensed_param_detectable

      
      ``bool``, ``dict``

      lensed_param_detectable is a dictionary of lensed parameters

      it will be populated when lensed_cbc_statistics() is called

      if unavailable, the lensed parameters will be sampled when lensed_rate() is called

      lensed_param_detectable.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time', 'lensed_images']















      ..
          !! processed by numpydoc !!

   .. py:attribute:: gw_param_sampler_dict

      
      ``dict``

      dictionary of params for initializing ``CompactBinaryPopulation`` class

      this will be used for GW unlensed parameters sampling

      gw_param_sampler_dict.keys() = ['nsamples', 'm_min', 'm_max', 'z_min', 'z_max', 'event_type', 'src_model_params']















      ..
          !! processed by numpydoc !!

   .. py:attribute:: lensed_param_sampler_dict

      
      ``dict``

      dictionary of params for initializing ``LensGalaxyPopulation`` class

      this will be used for GW lensed parameters sampling

      lensed_param_sampler_dict.keys() = ['nsamples', 'min_lensed_images', 'max_lensed_images', 'lensModelList']















      ..
          !! processed by numpydoc !!

   .. py:attribute:: snr_calculator_dict

      
      ``dict``

      dictionary of params for initializing ``snr_calculator`` (``gwsnr``) class

      this will be used for SNR calculation

      snr_calculator_dict.keys() = ['mtot_min', 'mtot_max', 'nsamples_mtot', 'nsamples_mass_ratio', 'sampling_frequency', 'waveform_approximant', 'minimum_frequency', 'snr_type', 'waveform_inspiral_must_be_above_fmin', 'psds', 'psd_file', 'ifos']















      ..
          !! processed by numpydoc !!

   .. py:attribute:: z_to_Dc

      
      ``scipy.interpolate.interp1d``

      redshift to co-moving distance.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: Dc_to_z

      
      ``scipy.interpolate.interp1d``

      co-moving distance to redshift.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: z_to_luminosity_distance

      
      ``scipy.interpolate.interp1d``

      redshift to luminosity distance.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: differential_comoving_volume

      
      ``scipy.interpolate.interp1d``

      differential comoving volume.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: compact_binary_pop

      
      ``CompactBinaryPopulation class``

      class for sampling GW parameters.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: lens_galaxy_pop

      
      ``LensGalaxyPopulation class``

      class for sampling lensed GW parameters.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: snr

      
      ``gwsnr package``

      class for calculating SNR.















      ..
          !! processed by numpydoc !!

   .. py:method:: class_initialization()

      
      Function for initializing the ``CompactBinaryPopulation`` and ``LensGalaxyPopulation`` classes.
















      ..
          !! processed by numpydoc !!

   .. py:method:: store_ler_params(json_file='./LeR_params.json')

      
      Fuction to store the parameters of the LER model. This is useful for reproducing the results.
















      ..
          !! processed by numpydoc !!

   .. py:method:: gwsnr_intialization(kwargs_dict)

      
      Function for initializing the `gwsnr <https://github.com/hemantaph/gwsnr/>`_ package.


      :Parameters:

          **kwargs_dict** : 'dict'
              keyword arguments for the initialization of the `gwsnr` package.
              kwargs_dict.keys()

              ``nsamples_mtot`` : `int`
                  nsamples_mtot = 200 (recommended for accurate results)
              ``nsamples_mass_ratio`` : `int`
                  nsamples_mass_ratio = 500 (recommended for accurate results)
              ``sampling_frequency`` : `float`
                  sampling_frequency = 4096. (recommended for accurate results)
              ``waveform_approximant`` : `str`
                  waveform_approximant = "IMRPhenomD" (for BBH) or "TaylorF2" (for BNS)
                  if you want to use other approximants, please set ``snr_type`` = 'inner_product'
              ``minimum_frequency`` : `float`
                  minimum_frequency = 20. (for O3 and O4 runs) or 10. (for 3G detectors)
              ``snr_type`` : `str`
                  snr_type = 'interpolation' (for fast results) or 'inner_product' (for bilby like results)
              ``waveform_inspiral_must_be_above_fmin`` : `bool`
                  False if dont want minimum frequency cut-off as higher mass BBH can merger below that frequency.
              ``psds`` : `bool` or `dict` or `str` (txt file)
                  e.g. For O4 design sensitivity

                      psds = {'L1':'aLIGOaLIGODesignSensitivityT1800044',

                      'H1':'aLIGOaLIGODesignSensitivityT1800044',

                      'V1':'AdvVirgo'}
              ``psd_file`` : `bool`, `list`
                  psd_file = False (if ASD) or True (if PSD file)
                  psd_file = [False,True] if psds[0] is a asd and psds[1] is a psd
              ``ifos`` : `list`
                  interferometer object name list
                  ifos = ['L1', 'H1', 'V1'] (for O4 design sensitivity)

      :Returns:

          **snr_** : `the gwsnr object`
              gwsnr object is used to calculate the SNR and pdet (probability of detection)













      ..
          !! processed by numpydoc !!

   .. py:method:: create_lookup_tables(z_min, z_max)

      
      To creating lookup tables for fast calculation for the following conversion operations,

      #. redshift to co-moving distance.
      #. co-moving distance to redshift.
      #. redshift to luminosity distance.

      :Parameters:

          **z_min** : `float`
              minimum redshift.
              for popI_II, popIII, primordial, BNS z_min = 0., 5., 5., 0. respectively.

          **z_max** : `float`
              maximum redshift.
              for popI_II, popIII, primordial, BNS z_max = 10., 40., 40., 2. respectively.












      :Attributes:

          **z_to_Dc** : `scipy.interpolate.interp1d`
              redshift to co-moving distance.

          **Dc_to_z** : `scipy.interpolate.interp1d`
              co-moving distance to redshift.

          **z_to_luminosity_distance** : `scipy.interpolate.interp1d`
              redshift to luminosity distance.

          **differential_comoving_volume** : `scipy.interpolate.interp1d`
              differential comoving volume.


      ..
          !! processed by numpydoc !!

   .. py:method:: batch_handler(nsamples, sampling_routine, json_file, resume=False)

      
      Function to handle the batch size.


      :Parameters:

          **nsamples** : `int`
              number of samples.

          **sampling_routine** : `function`
              function to sample the parameters.
              e.g. unlensed_sampling_routine() or lensed_sampling_routine()

          **json_file** : `str`
              name of the json file to store the parameters.

          **resume** : `bool`
              if True, it will resume the sampling from the last batch.
              default resume = False.














      ..
          !! processed by numpydoc !!

   .. py:method:: unlensed_sampling_routine(nsamples, file_name, resume=False)

      
      Function to generate unlensed GW source parameters.


      :Parameters:

          **nsamples** : `int`
              number of samples.
              default nsamples = 100000.

          **file_name** : `str`
              name of the json file to store the parameters.

          **resume** : `bool`
              if True, it will resume the sampling from the last batch.
              default resume = False.














      ..
          !! processed by numpydoc !!

   .. py:method:: unlensed_cbc_statistics(nsamples=None, resume=False, json_file='./gw_params.json', **kwargs)

      
      Function to generate unlensed GW source parameters.


      :Parameters:

          **nsamples** : `int`
              number of samples.
              default nsamples = 100000.

          **resume** : `bool`
              resume = False (default) or True.
              if True, the function will resume from the last batch.

          **json_file** : `str`
              json file name for storing the parameters.
              default json_file = './gw_params.json'.

          **kwargs** : `dict`
              key word arguments for initializing the ``CompactBinaryPopulation`` class.

              This initialization is either done at the time of class initialization or at the time of calling this function.

              Following parameters can be provided,

              ``m_min`` : `float`
                  minimum mass of the compact binary (single).
              ``m_max`` : `float`
                  maximum mass of the compact binary (single).
              ``event_type`` : `str`
                  event_type = 'popI_II' or `popIII` or `primordial`.
              ``src_model_params`` : `dict`
                  src_model_params = {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82,

                  'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08,

                  'mu_g': 33.07, 'sigma_g': 5.69}}

      :Returns:

          **unlensed_gw_params** : `dict`
              dictionary of unlensed GW source parameters.
              unlensed_gw_params.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time']













      ..
          !! processed by numpydoc !!

   .. py:method:: unlensed_rate(gw_param='./gw_params.json', snr_threshold=8.0, jsonfile='./gw_params_detectable.json')

      
      Function to calculate unlensed merger rate.

      .. math::
          R_U = \mathcal{N}^U\int dz_s R_o^U(z_s)\bigg\{\Theta[\rho(z_s,\theta)-\rho_{th}] P(\theta) d\theta \bigg\}

      - where :math:`\mathcal{N}^U` is the normalization factor of the unlensed merger rate distribution wrt redshift.

      :Parameters:

          **gw_param** : `dict` or `str` for json file name.
              dictionary of unlensed GW source parameters.
              default gw_param = './gw_params.json'.

          **snr_threshold** : `float`
              SNR threshold for detection.
              default snr_threshold = 8.

          **jsonfile** : `str`
              json file name for storing the detectable parameters.
              default jsonfile = './gw_params_detectable.json'.

      :Returns:

          **unlensed_rate** : (`float`,`float`)
              unlensed merger rate in a year
              unlensed_rate[0] = total unlensed rate with step function
              unlensed_rate[1] = total unlensed rate with pdet function

          **gw_param_detectable** : `dict`
              dictionary of detectable unlensed GW source parameters.
              gw_param_detectable.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time']













      ..
          !! processed by numpydoc !!

   .. py:method:: lensed_sampling_routine(nsamples, file_name, resume=False)

      
      Function to generate lensed GW source parameters, lens galaxy parameters and image paramters.


      :Parameters:

          **nsamples** : `int`
              number of samples.

          **file_name** : `str`
              name of the json file to store the parameters.

          **resume** : `bool`
              if True, it will resume the sampling from the last batch.
              default resume = False.














      ..
          !! processed by numpydoc !!

   .. py:method:: lensed_cbc_statistics(nsamples=None, resume=False, json_file='./lensed_params.json', **kwargs)

      
      Function to generate lensed GW source parameters, lens galaxy parameters and image paramters.


      :Parameters:

          **nsamples** : `int`
              number of samples.
              default nsamples = 100000.

          **resume** : `bool`
              resume = False (default) or True.
              if True, the function will resume from the last batch.

          **json_file** : `str`
              json file name for storing the parameters.
              default json_file = './lensed_params.json'.

          **kwargs** : `dict`
              key word arguments for initializing the ``LensGalaxyPopulation`` class.

              This initialization is either done at the time of class initialization or at the time of calling this function.

              Following parameters can be provided,

              ``min_lensed_images`` : `int`
                  minimum number of lensed images.
              ``max_lensed_images`` : `int`
                  maximum number of lensed images.
              ``lensModelList`` : `list`
                  list of lens models.
                  e.g. lensModelList = ['EPL_NUMBA', 'SHEAR'].

      :Returns:

          **lensed_param** : `dict`
              dictionary of lensed GW source parameters, lens galaxy parameters and image paramters.
              lensed_param.keys() = ['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2', 'Dl',
              'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source',
              'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'n_images',
              'x0_image_positions', 'x1_image_positions', 'magnifications', 'time_delays', 'traces',
              'determinants', 'image_type', 'weights', 'opt_snr_net', 'L1', 'H1', 'V1']













      ..
          !! processed by numpydoc !!

   .. py:method:: lensed_rate(lensed_param='./lensed_params.json', snr_threshold=8.0, num_img=2, jsonfile='./lensed_params_detectable.json', none_as_nan=True)

      
      Function to calculate lensed merger rate.

      .. math::
          R_L = \mathcal{N}^L\int dz_s R_o^L(z_s)\bigg\{\Theta[\rho(z_s,\theta)-\rho_{th}] P(\theta) d\theta \bigg\}

      - where :math:`\mathcal{N}^L` is the normalization factor of the lensed merger rate distribution wrt redshift.

      :Parameters:

          **lensed_param** : `dict` or `str` for the json file name.
              dictionary of lensed GW source parameters, lens galaxy parameters and image paramters.
              lensed_param.keys() = ['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2', 'Dl',
              'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source',
              'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'n_images',
              'x0_image_positions', 'x1_image_positions', 'magnifications', 'time_delays', 'traces',
              'determinants', 'image_type', 'weights', 'opt_snr_net', 'L1', 'H1', 'V1']

          **snr_threshold** : `float`
              threshold for detection signal to noise ratio.
              e.g. snr_threshold = 8.

          **num_img** : `int`
              number of images.
              e.g. num_img = 2.

          **jsonfile** : `str`
              json file name for storing the parameters.
              default jsonfile = './lensed_params_detectable.json'.

          **none_as_nan** : `bool`
              if True, replace None with np.nan in the lensed_param dictionary.
              default none_as_nan = True.

      :Returns:

          **lensed_rate** : `float`
              lensed merger rate in a year.
              lensed_rate[0] = total lensed rate with step function
              lensed_rate[1] = total lensed rate with pdet function

          **detectable_lensed_param** : `dict`
              dictionary of detectable lensed GW source parameters, lens galaxy parameters and image paramters.
              detectable_lensed_param.keys() = ['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2',
              'Dl', 'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source',
              'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'n_images',
              'x0_image_positions', 'x1_image_positions', 'magnifications', 'time_delays', 'traces',
              'determinants', 'image_type', 'weights', 'opt_snr_net', 'L1', 'H1', 'V1']













      ..
          !! processed by numpydoc !!

   .. py:method:: rate_comparision(snr_threshold_unlensed=8.0, unlened_param='./gw_params.json', snr_threshold_lensed=8.0, num_img=2, lensed_param='./lensed_params.json', jsonfile_unlensed='./gw_params_detectable.json', jsonfile_lensed='./lensed_params_detectable.json')

      
      Function to calculate unlensed and lensed merger rate and their ratio.


      :Parameters:

          **snr_threshold_unlensed** : `float`
              threshold for detection signal to noise ratio for unlensed case.
              e.g. snr_threshold_unlensed = 8.

          **unlened_param** : `dict`
              dictionary of unlensed GW source parameters.
              unlened_param.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time']

          **snr_threshold_lensed** : `float`
              threshold for detection signal to noise ratio for lensed case.
              e.g. snr_threshold_lensed = 8.

          **num_img** : `int`
              number of images crossing the threshold.
              e.g. num_img = 2.

          **lensed_param** : `dict`
              dictionary of lensed GW source parameters, lens galaxy parameters and image paramters.
              lensed_param.keys() = ['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2', 'Dl',
              'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source',
              'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'n_images',
              'x0_image_positions', 'x1_image_positions', 'magnifications', 'time_delays', 'traces',
              'determinants', 'image_type', 'weights', 'opt_snr_net', 'L1', 'H1', 'V1']

          **jsonfile_unlensed** : `str`
              json file name for storing the parameters for unlensed detectable case.
              default jsonfile_unlensed = './gw_params_detectable.json'.

          **jsonfile_lensed** : `str`
              json file name for storing the parameters for lensed detectable case.
              default jsonfile_lensed = './lensed_params_detectable.json'.

      :Returns:

          **unlensed_rate** : (`float`,`float`)
              unlensed merger rate in a year
              unlensed_rate[0] = total unlensed rate with step function
              unlensed_rate[1] = total unlensed rate with pdet function

          **lensed_rate** : (`float`,`float`)
              lensed merger rate in a year
              lensed_rate[0] = total lensed rate with step function
              lensed_rate[1] = total lensed rate with pdet function

          **rate_ratio** : (`float`,`float`)
              unlensed/lensed rate ratio
              rate_ratio[0] = total unlensed/lensed rate ratio with step function
              rate_ratio[1] = total unlensed/lensed rate ratio with pdet function













      ..
          !! processed by numpydoc !!

   .. py:method:: selecting_n_lensed_detectable_events(nsamples=100, snr_threshold=8.0, num_img=2, resume=False, json_file='./lensed_params_detectable.json')

      
      Function to select n lensed detectable events.


      :Parameters:

          **nsamples** : `int`
              number of samples to be selected.
              default size = 100.

          **snr_threshold** : `float`
              threshold for detection signal to noise ratio.
              e.g. snr_threshold = 8. or [8.,6.]

          **num_img** : `int`
              number of images crossing the threshold.
              e.g. num_img = 2 or [1,1]

          **resume** : `bool`
              if True, it will resume the sampling from the last batch.
              default resume = False.

          **json_file** : `str`
              json file name for storing the parameters.
              default json_file = './lensed_params_detectable.json'.

      :Returns:

          **param_final** : `dict`
              dictionary of lensed GW source parameters, lens galaxy parameters and image paramters.
              param_final.keys() = ['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2',
              'Dl', 'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source',
              'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'n_images',
              'x0_image_positions', 'x1_image_positions', 'magnifications', 'time_delays', 'image_type',
              'weights', 'opt_snr_net', 'L1', 'H1', 'V1']













      ..
          !! processed by numpydoc !!

   .. py:method:: relative_mu_dt_lensed(lensed_param, snr_threshold=[8.0, 8.0])

      
      Function to classify the lensed images wrt to the morse phase difference.


      :Parameters:

          **lensed_param** : `dict`
              dictionary of lensed GW source parameters, lens galaxy parameters and image paramters.
              lensed_param.keys() = ['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2', 'Dl',
              'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source',
              'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'n_images',
              'x0_image_positions', 'x1_image_positions', 'magnifications', 'time_delays', 'traces',
              'determinants', 'image_type', 'weights', 'opt_snr_net', 'L1', 'H1', 'V1']

          **snr_threshold** : `float`
              threshold for detection signal to noise ratio.
              e.g. snr_threshold = [8.,8.] or [8.,6.] for subthreshold

      :Returns:

          **mu_rel0** : `float.array`
              relative magnification for 0 degree phase difference.

          **dt_rel0** : `float.array`
              relative time delay for 0 degree phase difference.

          **mu_rel90** : `float.array`
              relative magnification for 90 degree phase difference.

          **dt_rel90** : `float.array`
              relative time delay for 90 degree phase difference.













      ..
          !! processed by numpydoc !!

   .. py:method:: mu_vs_dt_plot(x_array, y_array, savefig=False, ax=None, colors='blue', linestyles='-', origin='upper', alpha=0.6, extent=[0.01, 500.0, 0.01, 100.0], contour_levels=[0.1, 0.4, 0.68, 0.95])

      
      Function to generate 2D KDE and plot the relative magnification vs time delay difference for lensed samples.


      :Parameters:

          **x_array** : `float.array`
              x array.

          **y_array** : `float.array`
              y array.

          **xlabel** : `str`
              x label.

          **ylabel** : `str`
              y label.

          **title** : `str`
              title.

          **savefig** : `bool`
              if True, it will save the figure.
              default savefig = False.

          **ax** : `matplotlib.axes`
              matplotlib axes.
              default ax = None.

          **colors** : `str`
              color of the plot.
              default colors = 'blue'.

          **linestyles** : `str`
              linestyle of the plot.
              default linestyles = '-'.

          **origin** : `str`
              origin of the plot.
              default origin = 'upper'.

          **alpha** : `float`
              alpha of the plot.
              default alpha = 0.6.

          **extent** : `list`
              extent of the plot.
              default extent = [1e-2,5e2,1e-2,1e2].

          **contour_levels** : `list`
              contour levels of the plot.
              default contour_levels = [0.10,0.40,0.68,0.95] which corresponds to 1,2,3,4 sigma.

      :Returns:

          None
              ..













      ..
          !! processed by numpydoc !!

   .. py:method:: selecting_n_unlensed_detectable_events(nsamples=100, snr_threshold=8.0, resume=False, json_file='./gw_params_detectable.json')

      
      Function to select n unlensed detectable events.


      :Parameters:

          **nsamples** : `int`
              number of samples to be selected.
              default size = 100.

          **snr_threshold** : `float`
              threshold for detection signal to noise ratio.
              e.g. snr_threshold = 8.

          **resume** : `bool`
              if True, it will resume the sampling from the last batch.
              default resume = False.

          **json_file** : `str`
              json file name for storing the parameters.
              default json_file = './gw_params_detectable.json'.

      :Returns:

          **param_final** : `dict`
              dictionary of unlensed GW source parameters.
              param_final.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time']













      ..
          !! processed by numpydoc !!

   .. py:method:: relative_mu_dt_unlensed(param, size=100)

      
      Function to generate relative magnification vs time delay difference for unlensed samples.


      :Parameters:

          **param** : `dict`
              dictionary of unlensed GW source parameters.
              unlensed_param.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time']

      :Returns:

          **dmu** : `float.array`
              relative magnification.

          **dt** : `float.array`
              relative time delay.













      ..
          !! processed by numpydoc !!


.. py:class:: LensGalaxyPopulation(CompactBinaryPopulation_=False)

   
   Class to sample lens galaxy parameters


   :Parameters:

       **CompactBinaryPopulation_** : CompactBinaryPopulation class
           This is an already initialized class that contains a function (CompactBinaryPopulation.sample_gw_parameters) that actually samples the source parameters.

           :class:`~ler.source_population.CompactBinaryPopulation`











   .. rubric:: Examples

   >>> from ler.lens_galaxy_population import LensGalaxyPopulation
   >>> lens_pop = LensGalaxyPopulation()
   >>> # list all the methods of the class
   >>> print([method for method in dir(lens_pop) if method.startswith('__') is False])
   ['Dc_to_z', 'angular_diameter_distance', 'angular_diameter_distance_z1z2', 'cbc_pop', 'compute_einstein_radii', 'create_lookup_table', 'differential_comoving_volume', 'get_image_properties', 'get_lensed_snrs', 'lens_redshift_sampler_helper_function', 'm_max', 'm_min', 'normalization_pdf_z', 'rejection_sample_lensing_probability', 'sample_axis_ratio_angle_phi', 'sample_galaxy_shear', 'sample_gamma', 'sample_lens_parameters', 'sample_lens_parameters_routine', 'sample_lens_redshifts', 'sample_strongly_lensed_source_parameters', 'sample_velocity_dispersion_axis_ratio', 'strong_lensing_optical_depth', 'z_max', 'z_min', 'z_to_Dc', 'z_to_luminosity_distance']
   >>> # sample lens parameters
   >>> lens_parameters = lens_pop.sample_lens_parameters(size=1000)
   >>> lens_parameters.keys()
   dict_keys(['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2', 'Dl', 'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a2', 'tilt1', 'tilt2', 'phi12', 'phi_jl'])
   >>> # get image properties
   >>> lens_parameters = lens_pop.get_image_properties(lens_parameters, n_min_images=2, n_max_images=4, lensModelList=['EPL_NUMBA', 'SHEAR'], npool=4)
   solving lens equations...
   100%|█████████████████████████████████████████████████████████| 1000/1000 [00:00<00:00, 1258.38it/s]
   >>> lens_parameters.keys()
   dict_keys(['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2', 'Dl', 'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a2', 'tilt1', 'tilt2', 'phi12', 'phi_jl', 'n_images', 'x0_image_positions', 'x1_image_positions', 'magnifications', 'time_delays', 'image_type', 'weights'])
   >>> # get lensed SNRs
   >>> from gwsnr import GWSNR
   >>> snr = GWSNR()
   Given: IMR waveform
   psds not given. Choosing bilby's default psds
   given psds:  {'L1': 'aLIGO_O4_high_asd.txt', 'H1': 'aLIGO_O4_high_asd.txt', 'V1': 'AdV_asd.txt'}
   Interpolator will be generated for L1 detector at ./interpolator_pickle/L1/halfSNR_dict_0.pickle
   Interpolator will be generated for H1 detector at ./interpolator_pickle/H1/halfSNR_dict_0.pickle
   Interpolator will be generated for V1 detector at ./interpolator_pickle/V1/halfSNR_dict_0.pickle
   Generating interpolator for ['L1', 'H1', 'V1'] detectors
   interpolation for each mass_ratios: 100%|███████████████████████████| 50/50 [00:23<00:00,  2.10it/s]
   interpolator generated
   >>> lens_snrs = lens_pop.get_lensed_snrs(snr, lens_parameters, n_max_images=4)
   >>> lens_snrs.keys()
   dict_keys(['opt_snr_net', 'L1', 'H1', 'V1'])

   Instance Attributes
   ----------
   LensGalaxyPopulation class has the following instance attributes:

   +-------------------------------------+----------------------------------+
   | Atrributes                          | Type                             |
   +=====================================+==================================+
   |:attr:`~cbc_pop`                     | CompactBinaryPopulation class    |
   +-------------------------------------+----------------------------------+
   |:attr:`~z_min`                       | float                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~z_max`                       | float                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~m_min`                       | float                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~m_max`                       | float                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~normalization_pdf_z`         | float                            |
   +-------------------------------------+----------------------------------+

   Instance Methods
   ----------
   LensGalaxyPopulation class has the following instance methods:

   +-------------------------------------+----------------------------------+
   | Methods                             | Type                             |
   +=====================================+==================================+
   |:meth:`~create_lookup_table`         | Function to create a lookup      |
   |                                     | table for the differential       |
   |                                     | comoving volume and luminosity   |
   |                                     | distance wrt redshift            |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_lens_parameters`      | Function to sample lens galaxy   |
   |                                     | parameters                       |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_lens_parameters_routine`                                 |
   +-------------------------------------+----------------------------------+
   |                                     | Function to sample lens galaxy   |
   |                                     | parameters                       |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_strongly_lensed_source_parameters`                       |
   +-------------------------------------+----------------------------------+
   |                                     | Function to sample source        |
   |                                     | parameters conditioned on the    |
   |                                     | source being strongly lensed     |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_lens_redshifts`       | Function to sample lens redshifts|
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_velocity_dispersion_axis_ratio`                          |
   +-------------------------------------+----------------------------------+
   |                                     | Function to sample velocity      |
   |                                     | dispersion and axis ratio of the |
   |                                     | lens galaxy                      |
   +-------------------------------------+----------------------------------+
   |:meth:`~compute_einstein_radii`      | Function to compute the Einstein |
   |                                     | radii of the lens galaxies       |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_axis_ratio_angle_phi` | Function to sample the axis      |
   |                                     | rotation angle of the elliptical |
   |                                     | lens galaxy                      |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_galaxy_shear`         | Function to sample the lens      |
   |                                     | galaxy shear                     |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_gamma`                | Function to sample the lens      |
   |                                     | galaxy spectral index of the     |
   |                                     | density profile                  |
   +-------------------------------------+----------------------------------+
   |:meth:`~rejection_sample_lensing_probability`                           |
   +-------------------------------------+----------------------------------+
   |                                     | Function to conduct rejection    |
   |                                     | sampling wrt einstein radius     |
   +-------------------------------------+----------------------------------+
   |:meth:`~strong_lensing_optical_depth`| Function to compute the strong   |
   |                                     | lensing optical depth            |
   +-------------------------------------+----------------------------------+
   |:meth:`~get_image_properties`        | Function to get the image        |
   |                                     | properties e.g. image positions, |
   |                                     | magnifications, time delays, etc.|
   +-------------------------------------+----------------------------------+
   |:meth:`~get_lensed_snrs`             | Function to get the lensed SNRs  |
   +-------------------------------------+----------------------------------+



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

   .. py:method:: create_lookup_table(z_min, z_max)

      
      Functions to create lookup tables
      1. Redshift to co-moving distance.
      2. Co-moving distance to redshift.
      3. Redshift to luminosity distance
      4. Redshift to angular diameter distance.
      5. Lens redshift sampler helper function.
      6. Redshift to differential comoving volume.


      :Parameters:

          **z_min** : `float`
              minimum redshift

          **z_max** : `float`
              maximum redshift














      ..
          !! processed by numpydoc !!

   .. py:method:: sample_lens_parameters(size=1000, lens_parameters_input={}, verbose=False)

      
      Function to sample galaxy lens parameters


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **lens_parameters_input** : `dict`
              dictionary of lens parameters to sample

      :Returns:

          **lens_parameters** : `dict`
              dictionary of lens parameters and source parameters (lens conditions applied)
              e.g. dictionary keys:

              lensing related=>['zl':redshift of lens, 'zs': redshift of source, 'sigma':velocity dispersion, 'q':axis ratios, 'e1':ellipticity, 'e2':ellipticity, 'gamma1':external-shear, 'gamma2':external-shear, 'Dl':angular diameter distance of lens, 'Ds':angular diameter distance of source, 'Dls':angular diameter distance between lens and source, 'theta_E': einstein radius in radian, 'gamma':spectral index of mass density distribution]

              source related=>['mass_1': mass in detector frame (mass1>mass2), 'mass_2': mass in detector frame, 'mass_1_source':mass in source frame, 'mass_2_source':mass source frame, 'luminosity_distance': luminosity distance, 'iota': inclination angle, 'psi': polarization angle, 'phase': coalesence phase, 'geocent_time': coalensence GPS time at geocenter, 'ra': right ascension, 'dec': declination, 'a_1': spin magnitude of the more massive black hole, 'a2': spin magnitude of the less massive black hole, 'tilt_1': tilt angle of the more massive black hole, 'tilt_2': tilt angle of the less massive black hole, 'phi_12': azimuthal angle between the two spins, 'phi_jl': azimuthal angle between the total angular momentum and the orbital angular momentum]













      ..
          !! processed by numpydoc !!

   .. py:method:: sample_lens_parameters_routine(size=1000, lens_parameters_input={})

      
      Function to sample galaxy lens parameters


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **lens_parameters_input** : `dict`
              dictionary of lens parameters to sample

      :Returns:

          **lens_parameters** : `dict`
              dictionary of lens parameters and source parameters (lens conditions applied)
              e.g. dictionary keys:

              lensing related=>['zl':redshift of lens, 'zs': redshift of source, 'sigma':velocity dispersion, 'q':axis ratios, 'e1':ellipticity, 'e2':ellipticity, 'gamma1':external-shear, 'gamma2':external-shear, 'Dl':angular diameter distance of lens, 'Ds':angular diameter distance of source, 'Dls':angular diameter distance between lens and source, 'theta_E': einstein radius in radian, 'gamma':spectral index of mass density distribution]

              source related=>['mass_1': mass in detector frame (mass1>mass2), 'mass_2': mass in detector frame, 'mass_1_source':mass in source frame, 'mass_2_source':mass source frame, 'luminosity_distance': luminosity distance, 'iota': inclination angle, 'psi': polarization angle, 'phase': coalesence phase, 'geocent_time': coalensence GPS time at geocenter, 'ra': right ascension, 'dec': declination, 'a_1': spin magnitude of the more massive black hole, 'a2': spin magnitude of the less massive black hole, 'tilt_1': tilt angle of the more massive black hole, 'tilt_2': tilt angle of the less massive black hole, 'phi_12': azimuthal angle between the two spins, 'phi_jl': azimuthal angle between the total angular momentum and the orbital angular momentum]













      ..
          !! processed by numpydoc !!

   .. py:method:: sample_strongly_lensed_source_parameters(size=1000)

      
      Function to sample source redshifts and other parameters, conditioned on the source being strongly lensed.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **gw_param_strongly_lensed** : `dict`
              dictionary of source parameters. `zs` is sampled considering the merger rate density at source frame, comoving volume and strong lensing optical depth.

              e.g. gw_param_strongly_lensed.keys() = ['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a2', 'tilt1', 'tilt2', 'phi12', 'phi_jl']













      ..
          !! processed by numpydoc !!

   .. py:method:: sample_lens_redshifts(zs)

      
      Function to sample lens redshifts, conditioned on the lens being strongly lensed
      Input parameters:
          zs : source redshifts
      Output parameters:
          zl : lens redshifts
















      ..
          !! processed by numpydoc !!

   .. py:method:: sample_velocity_dispersion_axis_ratio(zs)

      
      Function to sample velocity dispersion and axis ratio of the lens galaxy


      :Parameters:

          **zs** : `float`
              source redshifts

      :Returns:

          **sigma** : `float`
              velocity dispersion of the lens galaxy

          **q** : `float`
              axis ratio of the lens galaxy













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













      ..
          !! processed by numpydoc !!

   .. py:method:: sample_axis_ratio_angle_phi(size=1000)

      
      Function to sample the axis rotation angle of the elliptical lens galaxy


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **phi** : `float`
              axis rotation angle of the elliptical lens galaxy













      ..
          !! processed by numpydoc !!

   .. py:method:: sample_galaxy_shear(size)

      
      Function to sample the lens galaxy shear


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **gamma_1** : `float`
              shear component in the x-direction

          **gamma_2** : `float`
              shear component in the y-direction













      ..
          !! processed by numpydoc !!

   .. py:method:: sample_gamma(size=1000)

      
      Function to sample the lens galaxy spectral index of the density profile


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **gamma** : `float`
              spectral index of the density profile













      ..
          !! processed by numpydoc !!

   .. py:method:: rejection_sample_lensing_probability(theta_E)

      
      Function to conduct rejection sampling wrt einstein radius


      :Parameters:

          **theta_E** : `float`
              Einstein radii of the lens galaxies

      :Returns:

          **idx** : `bool`
              boolean array of size len(theta_E) indicating whether the sample is accepted or not













      ..
          !! processed by numpydoc !!

   .. py:method:: strong_lensing_optical_depth(zs)

      
      Function to compute the strong lensing optical depth


      :Parameters:

          **zs** : `float`
              source redshifts

      :Returns:

          **tau** : `float`
              strong lensing optical depth













      ..
          !! processed by numpydoc !!

   .. py:method:: get_image_properties(lens_parameters, n_min_images=int(2), n_max_images=int(4), lensModelList=['EPL_NUMBA', 'SHEAR'], npool=4)

      
      Function to get the image properties e.g. image positions, magnifications, time delays, etc.


      :Parameters:

          **lens_parameters** : `dict`
              dictionary of lens parameters
              e.g. lens_parameters.keys() = ['zs', 'zl', 'gamma1', 'gamma2', 'e1', 'e2', 'gamma', 'theta_E']

          **n_min_images** : `int`
              minimum number of images to consider
              default: 2

          **n_max_images** : `int`
              maximum number of images to consider
              default: 4

          **lensModelList** : `list`
              list of lens models
              default: ['EPL_NUMBA', 'SHEAR']

          **npool** : `int`
              number of processes to use
              default: 4

      :Returns:

          **lens_parameters** : `dict`
              dictionary of lens parameters and image properties
              e.g. lens_parameters contains the following keys:

              lens related=>['zs': source redshift, 'zl': lens redshift, 'gamma1': shear component in the x-direction, 'gamma2': shear component in the y-direction, 'e1': ellipticity component in the x-direction, 'e2': ellipticity component in the y-direction, 'gamma': spectral index of the mass density distribution, 'theta_E': einstein radius in radian]

              source related=>['mass_1': mass in detector frame (mass1>mass2), 'mass_2': mass in detector frame, 'mass_1_source':mass in source frame, 'mass_2_source':mass source frame, 'luminosity_distance': luminosity distance, 'iota': inclination angle, 'psi': polarization angle, 'phase': coalesence phase, 'geocent_time': coalensence GPS time at geocenter, 'ra': right ascension, 'dec': declination, 'a_1': spin magnitude of the more massive black hole, 'a2': spin magnitude of the less massive black hole, 'tilt_1': tilt angle of the more massive black hole, 'tilt_2': tilt angle of the less massive black hole, 'phi_12': azimuthal angle between the two spins, 'phi_jl': azimuthal angle between the total angular momentum and the orbital angular momentum]

              image related=>['x_source': source position in the x-direction, 'y_source': source position in the y-direction, 'x0_image_position': image position in the x-direction, 'x1_image_position': image position in the y-direction, 'magnifications': magnifications, 'time_delays': time delays, 'n_images': number of images formed, 'determinant': determinants, 'trace': traces, 'iteration': to keep track of the iteration number, 'weights': weights for the caustic considered]













      ..
          !! processed by numpydoc !!

   .. py:method:: get_lensed_snrs(snr_calculator, lensed_param, n_max_images=4)

      
      Function to calculate the signal to noise ratio for each image in each event.


      :Parameters:

          **snr_calculator** : `class`
              snr_calculator class
              this is an already initialized class that contains a function (snr_calculator.snr) that actually calculates snr with the given gw_params.

              Luminosity distance and time delay are modified to be effective luminosity distance and effective time delay, respectively, for each image using the magnifications and time delays.

          **lensed_param** : `dict`
              dictionary containing the both already lensed source paramters and image parameters.
              e.g. lensed_param.keys() = ['mass_1', 'mass_2', 'zs', 'luminosity_distance', 'iota', 'psi', 'phi', 'ra', 'dec', 'geocent_time', 'phase', 'a_1', 'a2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl', 'magnifications', 'time_delays']

          **n_max_images** : `int`
              maximum number of images to consider
              default: 4

      :Returns:

          **snrs** : `dict`
              signal to noise ratio for each image in each event.
              (dictionary containing 'H1', 'L1', ..., and 'opt_snr_net', which is the network snr, for each image as an array with dimensions (number_of_lensed_events,n_max_images) )













      ..
          !! processed by numpydoc !!


.. py:class:: SourceGalaxyPopulationModel(z_min=0.0, z_max=10.0, event_type='popI_II', merger_rate_density_param=None)

   
   Class to generate a population of source galaxies.
   This class is inherited by :class:`~ler.ler.CompactBinaryPopulation` class.


   :Parameters:

       **z_min** : `float`
           Minimum redshift of the source population
           default: 0.

       **z_max** : `float`
           Maximum redshift of the source population
           default: 10.

       **event_type** : `str`
           Type of event to generate
           e.g. 'popI_II', 'BNS', 'popIII', 'primordial', 'popI_II_Madau_Dickinson'
           default: 'popI_II'











   .. rubric:: Examples

   >>> from ler import SourceGalaxyPopulationModel
   >>> pop = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, event_type = "popI_II")
   >>> zs = pop.sample_source_redshifts(size=1000)
   >>> zs
   array([0.0001, 0.0001, 0.0001, ..., 9.9999, 9.9999, 9.9999])

   Instance Attributes
   ----------
   SourceGalaxyPopulationModel has the following instance attributes:

   +-------------------------------------+----------------------------------+
   | Atrributes                          | Type                             |
   +=====================================+==================================+
   |:attr:`~z_min`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~z_max`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~normalization_pdf_z`         | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~z_to_luminosity_distance`    | `scipy.interpolate.interpolate`  |
   +-------------------------------------+----------------------------------+
   |:attr:`~differential_comoving_volume`| `scipy.interpolate.interpolate`  |
   +-------------------------------------+----------------------------------+

   Instance Methods
   ----------
   SourceGalaxyPopulationModel has the following instance methods:

   +-------------------------------------+----------------------------------+
   | Methods                             | Type                             |
   +=====================================+==================================+
   |:meth:`~create_lookup_table`         | Function to create a lookup      |
   |                                     | table for the differential       |
   |                                     | comoving volume and luminosity   |
   |                                     | distance wrt redshift            |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_source_redshifts`     | Function to sample source        |
   |                                     | redshifts from the source        |
   |                                     | galaxy population model          |
   +-------------------------------------+----------------------------------+
   |:meth:`~merger_rate_density_popI_II` | Function to compute the merger   |
   |                                     | rate density (PopI/PopII)        |
   +-------------------------------------+----------------------------------+
   |:meth:`~merger_rate_density_popI_II_Madau_Dickinson`                    |
   +-------------------------------------+----------------------------------+
   |                                     | Function to compute the          |
   |                                     | merger rate density (PopI/PopII) |
   |                                     | from Madau & Dickinson (2014)    |
   +-------------------------------------+----------------------------------+
   |:meth:`~merger_rate_density_popIII`  | Function to compute the merger   |
   |                                     | rate density (PopIII)            |
   +-------------------------------------+----------------------------------+
   |:meth:`~merger_rate_density_primordial`                                 |
   +-------------------------------------+----------------------------------+
   |                                     | Function to compute the merger   |
   |                                     | rate density (Primordial)        |
   +-------------------------------------+----------------------------------+



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

   .. py:attribute:: normalization_pdf_z

      
      ``float``

      Normalization constant of the pdf p(z)















      ..
          !! processed by numpydoc !!

   .. py:attribute:: z_to_luminosity_distance

      
      ``scipy.interpolate.interpolate``

      Function to convert redshift to luminosity distance















      ..
          !! processed by numpydoc !!

   .. py:attribute:: differential_comoving_volume

      
      ``scipy.interpolate.interpolate``

      Function to calculate the differential comoving volume















      ..
          !! processed by numpydoc !!

   .. py:method:: create_lookup_table(z_min, z_max)

      
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

   .. py:method:: sample_source_redshifts(size=1000, z_min=0.0, z_max=10.0)

      
      Function to sample source redshifts from the source galaxy population
      model


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **z_min** : `float`
              Minimum redshift of the source population

          **z_max** : `float`
              Maximum redshift of the source population

      :Returns:

          **zs** : `array`
              Array of sampled redshifts










      .. rubric:: Examples

      >>> from ler import SourceGalaxyPopulationModel
      >>> pop = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, event_type = "popI_II")
      >>> zs = pop.sample_source_redshifts(size=1000)
      >>> zs
      array([0.0001, 0.0001, 0.0001, ..., 9.9999, 9.9999, 9.9999])



      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_popI_II(zs, R0=23.9 * 1e-09, b2=1.6, b3=2.0, b4=30)

      
      Function to compute the merger rate density (PopI/PopII)


      :Parameters:

          **zs** : `float`
              Source redshifts

          **R0** : `float`
              Normalization constant
              default: 23.9*1e-9 Mpc^-3 yr^-1

          **b2** : `float`
              Fitting paramters
              default: 1.6

          **b3** : `float`
              Fitting paramters
              default: 2.0

          **b4** : `float`
              Fitting paramters
              default: 30

      :Returns:

          **rate_density** : `float`
              merger rate density










      .. rubric:: Examples

      >>> from ler import SourceGalaxyPopulationModel
      >>> pop = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, event_type = "popI_II")
      >>> rate_density = pop.merger_rate_density_popI_II(zs=0.1)
      >>> rate_density
      2.7848018586883885e-08



      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_popI_II_Madau_Dickinson(zs, af=2.7, bf=5.6, cf=1.9)

      
      Function to compute the unormalized merger rate density (PopI/PopII) from Madau & Dickinson (2014)


      :Parameters:

          **zs** : `float`
              Source redshifts

          **af** : `float`
              Fitting paramters
              default: 2.7

          **bf** : `float`
              Fitting paramters
              default: 5.6

          **cf** : `float`
              Fitting paramters
              default: 1.9

      :Returns:

          **rate_density** : `float`
              merger rate density










      .. rubric:: Examples

      >>> from ler import SourceGalaxyPopulationModel
      >>> pop = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, event_type = "popI_II_Madau_Dickinson")
      >>> rate_density = pop.merger_rate_density_popI_II_Madau_Dickinson(zs=0.1)
      >>> rate_density
      1.2355851838964846



      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_popIII(zs, aIII=0.66, bIII=0.3, zIII=11.6)

      
      Function to compute the unnormalized merger rate density (PopIII)


      :Parameters:

          **zs** : `float`
              Source redshifts

          **aIII** : `float`
              Fitting paramters
              default: 0.66

          **bIII** : `float`
              Fitting paramters
              default: 0.3

          **zIII** : `float`
              Fitting paramters
              default: 11.6

      :Returns:

          **rate_density** : `float`
              merger rate density










      .. rubric:: Examples

      >>> from ler import SourceGalaxyPopulationModel
      >>> pop = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, event_type = "popIII")
      >>> rate_density = pop.merger_rate_density_popIII(zs=0.1)
      >>> rate_density
      0.00010000000000000002



      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_primordial(zs, t0=13.786885302009708)

      
      Function to compute the merger rate density (Primordial)


      :Parameters:

          **zs** : `float`
              Source redshifts

          **t0** : `float`
              Present ge of the Universe in Gyr
              default: 13.786885302009708

      :Returns:

          **rate_density** : `float`
              merger rate density










      .. rubric:: Examples

      >>> from ler import SourceGalaxyPopulationModel
      >>> pop = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, event_type = "primordial")
      >>> rate_density = pop.merger_rate_density_primordial(zs=0.1)
      >>> rate_density
      0.00010000000000000002



      ..
          !! processed by numpydoc !!


.. py:class:: CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type='popI_II', merger_rate_density_param=None, src_model_params=None)

   Bases: :py:obj:`SourceGalaxyPopulationModel`

   
   Class to generate a population of compact binaries. Inherits from :class:`~ler.ler.SourceGalaxyPopulationModel` class.


   :Parameters:

       **z_min** : `float`
           Minimum redshift of the source population

       **z_max** : `float`
           Maximum redshift of the source population

       **m_min** : `float`
           Minimum mass of the BBHs

       **m_max** : `float`
           Maximum mass of the BBHs

       **event_type** : `str`
           Type of event to generate.
           e.g. 'popI_II', 'BNS', 'popIII', 'primordial', 'popI_II_Madau_Dickinson'

       **src_model_params** : `dict`
           Dictionary of model parameters.
           e.g. for popI_II: {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, 'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, 'mu_g': 33.07, 'sigma_g': 5.69}











   .. rubric:: Examples

   >>> from ler import CompactBinaryPopulation
   >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popI_II")
   >>> gw_parameters = pop.sample_gw_parameters(nsamples=1000)
   >>> gw_parameters.keys()
   dict_keys(['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl'])

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
   |:attr:`~m_min`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~m_max`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~event_type`                  | `str`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~src_model_params`                  | `dict`                           |
   +-------------------------------------+----------------------------------+

   Instance Methods
   ----------
   CompactBinaryPopulation has the following instance methods:

   +-------------------------------------+----------------------------------+
   | Methods                             | Type                             |
   +=====================================+==================================+
   |:meth:`~sample_gw_parameters`        | Function for sampling GW         |
   |                                     | parameters from the source       |
   |                                     | galaxy population model          |
   +-------------------------------------+----------------------------------+
   |:meth:`~binary_masses_popI_II`       | Function to calculate source     |
   |                                     | mass1 and mass2 with             |
   |                                     | PowerLaw+PEAK model              |
   +-------------------------------------+----------------------------------+
   |:meth:`~binary_masses_popIII`        | Function to calculate source     |
   |                                     | mass1 and mass2 with pop III     |
   |                                     | origin                           |
   +-------------------------------------+----------------------------------+
   |:meth:`~binary_masses_primordial`    | Function to calculate source     |
   |                                     | mass1 and mass2 for primordial   |
   |                                     | BBHs                             |
   +-------------------------------------+----------------------------------+
   |:meth:`~binary_masses_BNS`           | Function to calculate source     |
   |                                     | mass1 and mass2 of BNS           |
   +-------------------------------------+----------------------------------+
   |:meth:`~mass_ratio`                  | Function to calculate mass ratio |
   +-------------------------------------+----------------------------------+



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

   .. py:attribute:: m_min

      
      ``float``

      Minimum mass of the BBHs















      ..
          !! processed by numpydoc !!

   .. py:attribute:: m_max

      
      ``float``

      Maximum mass of the BBHs















      ..
          !! processed by numpydoc !!

   .. py:attribute:: event_type

      
      ``str``

      Type of event to generate.

      e.g. 'popI_II', 'BNS', 'popIII', 'primordial', 'popI_II_Madau_Dickinson'















      ..
          !! processed by numpydoc !!

   .. py:attribute:: src_model_params

      
      ``dict``

      Dictionary of model parameters.

      e.g. for popI_II: {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, 'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, 'mu_g': 33.07, 'sigma_g': 5.69}

      for popI_II_Madau_Dickinson: {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, 'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, 'mu_g': 33.07, 'sigma_g': 5.69}

      for popIII: None

      for primordial: {'Mc':30.,'sigma':0.3,'beta':1.1}

      for BNS: None












      .. rubric:: Examples

      >>> from ler import CompactBinaryPopulation
      >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popI_II")
      >>> method_list = [method for method in dir(pop) if method.startswith('__') is False]
      >>> print(method_list)
      ['create_lookup_table', 'differential_comoving_volume', 'merger_rate_density', 'merger_rate_density_popIII', 'merger_rate_density_popI_II', 'merger_rate_density_popI_II_Madau_Dickinson', 'merger_rate_density_primordial', 'normalization_pdf_z', 'sample_source_redshifts', 'z_max', 'z_min', 'z_to_luminosity_distance']



      ..
          !! processed by numpydoc !!

   .. py:method:: sample_gw_parameters(nsamples=1000, verbose=False, **kwargs)

      
      Function to sample BBH parameters from the source galaxy population
      model


      :Parameters:

          **nsamples** : `int`
              Number of samples to draw

          **kwargs** : `dict`
              Keyword arguments to pass in parameter values
              e.g. zs = np.array([0.1,0.2,0.3])

      :Returns:

          **gw_parameters** : `dict`
              Dictionary of sampled parameters
              gw_parameters.keys() = ['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl']










      .. rubric:: Examples

      >>> from ler import CompactBinaryPopulation
      >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popI_II")
      >>> gw_parameters = pop.sample_gw_parameters(nsamples=1000)
      >>> gw_parameters.keys()
      dict_keys(['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl'])



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_popI_II(size, alpha=3.63, beta=1.26, delta_m=4.82, mmin=4.59, mmax=86.22, lambda_peak=0.08, mu_g=33.07, sigma_g=5.69)

      
      Function to calculate source mass1 and mass2 with PowerLaw+PEAK model


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **src_model_params** : `dict`
              Dictionary of model parameters
              e.g. {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, 'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, 'mu_g': 33.07, 'sigma_g': 5.69}

      :Returns:

          **mass_1_source** : `array`
              Array of mass1 in source frame

          **mass_2_source** : `array`
              Array of mass2 in source frame










      .. rubric:: Examples

      >>> from ler import CompactBinaryPopulation
      >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popI_II")
      >>> src_model_params = {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, 'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, 'mu_g': 33.07, 'sigma_g': 5.69}
      >>> mass_1_source, mass_2_source = pop.binary_masses_popI_II(size=1000, src_model_params=src_model_params)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_popIII(size)

      
      Function to calculate source mass1 and mass2 with pop III origin


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **src_model_params** : `dict`
              Dictionary of model parameters

      :Returns:

          **mass_1_source** : `array`
              Array of mass1 in source frame

          **mass_2_source** : `array`
              Array of mass2 in source frame










      .. rubric:: Examples

      >>> from ler import CompactBinaryPopulation
      >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popIII")
      >>> mass_1_source, mass_2_source = pop.binary_masses_popIII(size=1000, src_model_params=None)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_primordial(size, Mc=30.0, sigma=0.3, beta=1.1)

      
      Function to calculate source mass1 and mass2 for primordial BBHs


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **src_model_params** : `dict`
              Dictionary of model parameters
              e.g. {'Mc':30.,'sigma':0.3,'beta':1.1}

      :Returns:

          **mass_1_source** : `array`
              Array of mass1 in source frame

          **mass_2_source** : `array`
              Array of mass2 in source frame










      .. rubric:: Examples

      >>> from ler import CompactBinaryPopulation
      >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "primordial")
      >>> src_model_params = {'Mc':30.,'sigma':0.3,'beta':1.1}
      >>> mass_1_source, mass_2_source = pop.binary_masses_primordial(size=1000, src_model_params=src_model_params)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BNS(size)

      
      Function to calculate source mass1 and mass2 of BNS


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **src_model_params** : `dict`
              Dictionary of model parameters

      :Returns:

          **mass_1_source** : `array`
              Array of mass1 in source frame

          **mass_2_source** : `array`
              Array of mass2 in source frame










      .. rubric:: Examples

      >>> from ler import CompactBinaryPopulation
      >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=1.0, m_max=3.0, event_type = "BNS")
      >>> mass_1_source, mass_2_source = pop.binary_masses_BNS(size=1000, src_model_params=None)



      ..
          !! processed by numpydoc !!

   .. py:method:: mass_ratio(size, beta=1.1)

      
      Function to calculate mass ratio with power law distribution


      :Parameters:

          **size** : `int`
              Number of samples

          **beta** : `float`
              Power law index

      :Returns:

          **q** : `array`
              Array of mass ratio










      .. rubric:: Examples

      >>> from ler import CompactBinaryPopulation
      >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=1.0, m_max=3.0, event_type = "BNS")
      >>> q = pop.mass_ratio(size=1000, beta=1.1)



      ..
          !! processed by numpydoc !!


.. py:function:: solve_lens_equation1(lens_parameters)

   
   Function to solve the lens equation (min_image = 2)


   :Parameters:

       **lens_parameters** : `list`
           a list of parameters
           lens_parameters[0] = min_images : minimum number of images
           lens_parameters[1] = e1 : ellipticity
           lens_parameters[2] = e2 : ellipticity
           lens_parameters[3] = gamma : power-law index
           lens_parameters[4] = gamma1 : shear
           lens_parameters[5] = gamma2 : shear
           lens_parameters[6] = zl : redshift of the lens
           lens_parameters[7] = zs : redshift of the source
           lens_parameters[8] = einstein_radius : Einstein radius
           lens_parameters[9] = iteration : iteration number
           lens_parameters[10:] = lens_model_list : numpy array of lens models

   :Returns:

       **x_source** : `float`
           x position of the source in the source plane

       **y_source** : `float`
           y position of the source in the source plane

       **x0_image_position** : `float`
           x position of the images in the source plane

       **x1_image_position** : `float`
           y position of the images in the source plane

       **magnifications** : `float`
           magnification of the images

       **time_delays** : `float`
           time-delay of the images

       **nImages** : `int`
           number of images

       **determinant** : `float`
           determinant of the hessian matrix

       **trace** : `float`
           trace of the hessian matrix

       **iteration** : `int`
           iteration number

       **weights** : `float`
           weights for the caustic










   .. rubric:: Examples

   >>> from ler.multiprocessing_routine import solve_lens_equation1
   >>> import numpy as np
   >>> from multiprocessing import Pool
   >>> # lens parameters input contains 12 parameters [e1, e2, gamma, gamma1, gamma2, zl, zs, einstein_radius, iteration, lens_model_list]
   >>> lens_parameters1 = np.array([2, 0.024069457093642648, -0.016002190961948142, 1.8945414936459974, 0.10117465203892329, 0.09600089396968613, 0.2503743800068136, 0.9418211055453296, 2.5055790287104725e-06, 0, 'EPL_NUMBA', 'SHEAR'], dtype=object)
   >>> lens_parameters2 = np.array([2, -0.04030088581646998, -0.01419438113690042, 2.0068239327017, 0.08482718989370612, -0.015393332086560785, 1.0952303138971118, 2.5534097159384417, 1.0125570159563301e-06, 1, 'EPL_NUMBA', 'SHEAR'], dtype=object)
   >>> input_arguments = np.vstack((lens_parameters1, lens_parameters2))
   >>> # solve the lens equation for each set of lens parameters
   >>> with Pool(2) as p:
   ...     result = p.map(solve_lens_equation1, input_arguments)
   >>> # result is a list of tuples
   >>> # each tuple contains the output parameters of the function
   >>> # each output parameter contains x_source, y_source, x0_image_position, x1_image_position, magnifications, time_delays, nImages, determinant, trace, iteration, weights
   >>> print(f"magnification of images with lens parameters 'lens_parameters1' is {result[0][6]}")
   magnification of images with lens parameters 'lens_parameters1' is [ 2.18973765 -1.27542831]



   ..
       !! processed by numpydoc !!

.. py:function:: solve_lens_equation2(lens_parameters)

   
   Function to solve the lens equation (min_image > 2)


   :Parameters:

       **lens_parameters** : `list`
           a list of parameters
           lens_parameters[0] = min_images : minimum number of images
           lens_parameters[1] = e1 : ellipticity
           lens_parameters[2] = e2 : ellipticity
           lens_parameters[3] = gamma : power-law index
           lens_parameters[4] = gamma1 : shear
           lens_parameters[5] = gamma2 : shear
           lens_parameters[6] = zl : redshift of the lens
           lens_parameters[7] = zs : redshift of the source
           lens_parameters[8] = einstein_radius : Einstein radius
           lens_parameters[9] = iteration : iteration number
           lens_parameters[10:] = lens_model_list : numpy array of lens models

   :Returns:

       **x_source** : `float`
           x position of the source in the source plane

       **y_source** : `float`
           y position of the source in the source plane

       **x0_image_position** : `float`
           x position of the images in the source plane

       **x1_image_position** : `float`
           y position of the images in the source plane

       **magnifications** : `float`
           magnification of the images

       **time_delays** : `float`
           time-delay of the images

       **nImages** : `int`
           number of images

       **determinant** : `float`
           determinant of the hessian matrix

       **trace** : `float`
           trace of the hessian matrix

       **iteration** : `int`
           iteration number

       **weights** : `float`
           weights for the caustic










   .. rubric:: Examples

   >>> from ler.multiprocessing_routine import solve_lens_equation2
   >>> import numpy as np
   >>> from multiprocessing import Pool
   >>> # lens parameters input contains 12 parameters [e1, e2, gamma, gamma1, gamma2, zl, zs, einstein_radius, iteration, lens_model_list]
   >>> lens_parameters1 = np.array([3, 0.024069457093642648, -0.016002190961948142, 1.8945414936459974, 0.10117465203892329, 0.09600089396968613, 0.2503743800068136, 0.9418211055453296, 2.5055790287104725e-06, 0, 'EPL_NUMBA', 'SHEAR'], dtype=object)
   >>> lens_parameters2 = np.array([3, -0.04030088581646998, -0.01419438113690042, 2.0068239327017, 0.08482718989370612, -0.015393332086560785, 1.0952303138971118, 2.5534097159384417, 1.0125570159563301e-06, 1, 'EPL_NUMBA', 'SHEAR'], dtype=object)
   >>> input_arguments = np.vstack((lens_parameters1, lens_parameters2))
   >>> # solve the lens equation for each set of lens parameters
   >>> with Pool(2) as p:
   ...     result = p.map(solve_lens_equation2, input_arguments)
   >>> # result is a list of tuples
   >>> # each tuple contains the output parameters of the function
   >>> # each output parameter contains x_source, y_source, x0_image_position, x1_image_position, magnifications, time_delays, nImages, determinant, trace, iteration, weights
   >>> print(f"magnification of images with lens parameters 'lens_parameters1' is {result[0][6]}")
   magnification of images with lens parameters 'lens_parameters1' is [ 2.18973765 -1.27542831]



   ..
       !! processed by numpydoc !!

.. py:function:: add_dictionaries_together(dictionary1, dictionary2)

   
   Adds two dictionaries with the same keys together.
















   ..
       !! processed by numpydoc !!

.. py:function:: rejection_sample(pdf, xmin, xmax, size=100)

   
   Helper function for rejection sampling from a pdf with maximum and minimum arguments.
   Input parameters:
       pdf: the pdf to sample from
       xmin: the minimum argument of the pdf
       xmax: the maximum argument of the pdf
       size: the number of samples to draw
   Output:
       samples: the samples drawn from the pdf
















   ..
       !! processed by numpydoc !!

