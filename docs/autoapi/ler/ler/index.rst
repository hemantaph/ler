:py:mod:`ler.ler`
=================

.. py:module:: ler.ler

.. autoapi-nested-parse::

   This module contains the main class for calculating the rates of lensed and unlensed events.

   ..
       !! processed by numpydoc !!


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   ler.ler.LeR




.. py:class:: LeR(nsamples=100000, npool=int(4), z_min=0.0, z_max=10.0, batch_size=25000, snr_finder='gwsnr', **kwargs)

   
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

      gw_param_sampler_dict.keys() = ['nsamples', 'm_min', 'm_max', 'z_min', 'z_max', 'event_type', 'model_pars']















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

   .. py:method:: store_ler_params()

      
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
              ``model_pars`` : `dict`
                  model_pars = {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82,

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


