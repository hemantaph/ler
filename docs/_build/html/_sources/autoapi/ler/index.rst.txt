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
   ler/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   ler.LeR
   ler.SourceGalaxyPopulationModel
   ler.CompactBinaryPopulation



Functions
~~~~~~~~~

.. autoapisummary::

   ler.solve_lens_equation1
   ler.solve_lens_equation2
   ler.add_dictionaries_together
   ler.rejection_sample



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

       +---------------------------------------+--------------------------------------------------+
       | Atrributes                            | Type                                             |
       +=======================================+==================================================+
       | :attr:`~gw_param`                     | ``dict``                                         |
       +---------------------------------------+--------------------------------------------------+
       | :attr:`~gw_param_detectable`          | ``dict``                                         |
       +---------------------------------------+--------------------------------------------------+
       | :attr:`~lensed_param`                 | ``dict``                                         |
       +---------------------------------------+--------------------------------------------------+
       | :attr:`~lensed_param_detectable`      | ``dict``                                         |
       +---------------------------------------+--------------------------------------------------+
       | :attr:`~gw_param_sampler_dict`        | ``dict``                                         |
       +---------------------------------------+--------------------------------------------------+
       | :attr:`~lensed_param_sampler_dict`    | ``dict``                                         |
       +---------------------------------------+--------------------------------------------------+
       | :attr:`~snr_calculator_dict`          | ``dict``                                         |
       +---------------------------------------+--------------------------------------------------+
       | :attr:`~z_to_Dc`                      | ``scipy.interpolate.interp1d``                   |
       +---------------------------------------+--------------------------------------------------+
       | :attr:`~Dc_to_z`                      | ``scipy.interpolate.interp1d``                   |
       +---------------------------------------+--------------------------------------------------+
       | :attr:`~z_to_luminosity_distance`     | ``scipy.interpolate.interp1d``                   |
       +---------------------------------------+--------------------------------------------------+
       | :attr:`~differential_comoving_volume` | ``scipy.interpolate.interp1d``                   |
       +---------------------------------------+--------------------------------------------------+
       | :attr:`~compact_binary_pop`           | ``CompactBinaryPopulation class``                |
       +---------------------------------------+--------------------------------------------------+
       | :attr:`~lens_galaxy_pop`              | ``LensGalaxyPopulation class``                   |
       +---------------------------------------+--------------------------------------------------+
       | :attr:`~snr`                          | ``gwsnr package``                                |
       +---------------------------------------+--------------------------------------------------+

   Instance Methods
   ----------
   LeR class has the following method(s),

       +---------------------------------------+--------------------------------------------------+
       | Method(s)                             | Description                                      |
       +=======================================+==================================================+
       | :meth:`~gwsnr_intialization`          | Function for initializing the ``gwsnr`` package. |
       +---------------------------------------+--------------------------------------------------+
       | :meth:`~create_lookup_tables`         | To creating lookup tables for fast calculation   |
       |                                       | for the following conversion operations,         |
       |                                       | redshift to co-moving distance.                  |
       |                                       | co-moving distance to redshift.                  |
       |                                       | redshift to luminosity distance.                 |
       +---------------------------------------+--------------------------------------------------+
       | :meth:`~unlensed_cbc_statistics`      | Function to generate unlensed GW source          |
       |                                       | parameters.                                      |
       +---------------------------------------+--------------------------------------------------+
       | :meth:`~unlensed_rate`                | Function to calculate unlensed merger rate.      |
       +---------------------------------------+--------------------------------------------------+
       | :meth:`~lensed_cbc_statistics`        | Function to generate lensed GW source            |
       |                                       | parameters.                                      |
       +---------------------------------------+--------------------------------------------------+
       | :meth:`~lensed_rate`                  | Function to calculate lensed merger rate.        |
       +---------------------------------------+--------------------------------------------------+



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

   .. py:attribute:: gw_param_detectable

      
      ``dict``

      gw_param_detectable is a dictionary of unlensed parameters for detectable sources

      it will be populated when unlensed_rate() is called

      gw_param_detectable.keys() = ['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'opt_snr_net', 'L1', 'H1', 'V1', 'pdet_net']















      ..
          !! processed by numpydoc !!

   .. py:attribute:: lensed_param

      
      ``bool``, ``dict``

      lensed_param is a dictionary of lensed parameters (for lens galaxy, source properties and image properties)

      it will be populated when lensed_cbc_statistics() is called

      if unavailable, then lensed parameters will be sampled when lensed_rate() is called

      lensed_param.keys() = ['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2', 'Dl', 'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'n_images', 'x0_image_positions', 'x1_image_positions', 'magnifications', 'time_delays', 'traces', 'determinants', 'image_type', 'weights', 'opt_snr_net', 'L1', 'H1', 'V1', 'pdet_net']















      ..
          !! processed by numpydoc !!

   .. py:attribute:: lensed_param_detectable

      
      ``dict``

      lensed_param_detectable is a dictionary of lensed parameters for detectable sources

      it will be populated when lensed_rate() is called

      lensed_param_detectable.keys() = ['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2', 'Dl', 'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'n_images', 'x0_image_positions', 'x1_image_positions', 'magnifications', 'time_delays', 'traces', 'determinants', 'image_type', 'weights', 'opt_snr_net', 'L1', 'H1', 'V1', 'pdet_net']















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

   .. py:method:: unlensed_cbc_statistics(nsamples=False, jsonfile=False, **kwargs)

      
      Function to generate unlensed GW source parameters.


      :Parameters:

          **nsamples** : `int`
              number of samples.

          **jsonfile** : `bool`, `str`
              if `str`, store all gravitational waves source parameters in jsonfile

              (for all sources, detected and undetected).

          **kwargs** : `dict`
              if new paramteres are provided, it will be used for sampling source parameters.
              Following parameters can be provided,

              ``m_min`` : `float`
                  minimum mass of the compact binary (single).
              ``m_max`` : `float`
                  maximum mass of the compact binary (single).
              ``event_type`` : `str`
                  event_type = 'popI_II' etc.
              ``model_pars`` : `dict`
                  model_pars = {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82,

                  'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08,

                  'mu_g': 33.07, 'sigma_g': 5.69}}

      :Returns:

          **unlensed_gw_params** : `dict`
              dictionary of unlensed GW source parameters.
              unlensed_gw_params.keys() = ['m1', 'm2', 'z', 'Dc', 'Dl', 'Dl_obs', 'snr', 'pdet', 'event_type']













      ..
          !! processed by numpydoc !!

   .. py:method:: unlensed_rate(size=False, snr_threshold=8.0, jsonfile=True)

      
      Function to calculate unlensed merger rate.

      .. math::
          R_U = \mathcal{N}^U\int dz_s R_o^U(z_s)\bigg\{\Theta[\rho(z_s,\theta)-\rho_{th}] P(\theta) d\theta \bigg\}

      - where :math:`\mathcal{N}^U` is the normalization factor of the unlensed merger rate distribution wrt redshift.

      :Parameters:

          **size** : `int`
              number of samples.

          **snr_threshold** : `float`
              threshold for detection signal to noise ratio.
              e.g. snr_threshold = 8.

          **jsonfile** : `bool`
              if True, store all gravitational waves source parameters in json file

              (for detected sources).

      :Returns:

          **unlensed_rate** : `float`
              unlensed merger rate in a year.













      ..
          !! processed by numpydoc !!

   .. py:method:: lensed_cbc_statistics(nsamples=False, jsonfile=True, **kwargs)

      
      function to generate lensed GW source parameters, lens parameters and image parameters
      Intput Parameters:
          nsamples: number of samples
          snr_threshold: threshold for detection signal to noise ratio
          jsonfile: if True, store lensed GW source parameters, lens parameters and image parameters in json file
                      (both for detected and undetected sources)
          **kwargs: if new parameters are provided, it will be used for sampling
      Output Parameters:
          lensed_param: `dict`ionary of lensed GW source parameters, lens parameters and image parameters
















      ..
          !! processed by numpydoc !!

   .. py:method:: lensed_rate(size=False, snr_threshold=8.0, num_img=2, jsonfile=True, none_as_nan=True)

      
      Function to calculate detectable lensed merger rate
      Intput Parameters:
          size (int): number of samples
          snr_threshold (float/array): threshold for detection signal to noise ratio
          num_img (int/array): number of images
                              e.g. For Sub-thershold events, snr_threshold=[8.,6.], num_img=[1,1]
                              The event will contain 1 image with snr>8 and 1 image with snr>6
          jsonfile (bool): if True, store all gravitational waves source parameters in json file
          none_as_nan (bool): if True,  no value is kept as np.nan
                              if False, no value is kept as 0.
      Output Parameters:
          lensed_rate (float): lensed merger rate in yr^-1
















      ..
          !! processed by numpydoc !!

   .. py:method:: selecting_n_lensed_detectable_events_from_dict(snr_threshold=8.0, num_img=2, none_as_nan=True, lenstype='I')

      
      Function to select n lensed detectable events from self.lensed_param
      Input Parameters:
          snr_threshold (float/array): threshold for detection signal to noise ratio
          num_img (int/array): number of images
                              e.g. For Sub-thershold events, snr_threshold=[8.,6.], num_img=[1,1]
                              The event will contain 1 image with snr>8 and 1 image with snr>6
          none_as_nan (bool): if True,  no value is kept as np.nan
                              if False, no value is kept as 0.
          lenstype (str): lens type, 'I' or 'II'
      Output Parameters:
          lensed_param (dict): `dict`ionary of lensed parameters
















      ..
          !! processed by numpydoc !!

   .. py:method:: selecting_n_lensed_detectable_events_with_sampling(snr_threshold=8.0, num_img=2, none_as_nan=True, size=100, lenstype='I', min_img=2, max_img=4)

      
      Function to select n lensed detectable events with sampling
















      ..
          !! processed by numpydoc !!

   .. py:method:: selecting_n_unlensed_detectable_events_from_dict(snr_threshold=8.0)

      
      Function to select n lensed detectable events from self.gw_param
















      ..
          !! processed by numpydoc !!

   .. py:method:: selecting_n_unlensed_detectable_events_with_sampling(snr_threshold=8.0, size=100)

      
      Function to select n lensed detectable events with sampling
















      ..
          !! processed by numpydoc !!

   .. py:method:: selecting_n_detectable_events(snr_threshold=8.0, num_img=2, none_as_nan=True, lensed=True, jsonfile=True, lenstype='I', new=False, size=100, min_img=2, max_img=4, batch_size=10000, **kwargs)

      
      Function to select n detectable events
      Input parameters:
          snr_threshold (float): the threshold for the SNR
          num_img (int): the number of images
          none_as_nan (bool): if True, then replace None with np.nan
          lensed (bool): if True, then select lensed events
          jsonfile (bool): if True, then save the dictionary as a json file
          lenstype (str): the lens type
                          e.g. 'I' for image type I, 'II' for image type II, 'III' for image type III, 'any' for any type
          new (bool): if True, then sample new events
          size (int): the number of events to be sampled
          min_img (int): the minimum number of images
          max_img (int): the maximum number of images
          batch_size (int): the number of events to be sampled in each iteration
      Output parameters:
          param (dict): the dictionary containing the parameters of the selected events
















      ..
          !! processed by numpydoc !!

   .. py:method:: rate_comparision(size=False, snr_threshold=8.0, num_img=2, jsonfile=True, none_as_nan=True)

      
      Function to compare the detectable lensed merger rate with the unlensed merger rate
      Intput Parameters:
          size (int): number of samples
          snr_threshold (float/array): threshold for detection signal to noise ratio
          num_img (int/array): number of images
                              e.g. For Sub-thershold events, snr_threshold=[8.,6.], num_img=[1,1]
                              The event will contain 1 image with snr>8 and 1 image with snr>6
          jsonfile (bool): if True, store all gravitational waves source parameters in json file
          none_as_nan (bool): if True,  no value is kept as np.nan
                              if False, no value is kept as 0.
      Output Parameters:
          unlened_rate (float): unlensed merger rate in yr^-1
          lensed_rate (float): lensed merger rate in yr^-1
          rate_ratio (float): lensed/unlensed merger rate ratio
















      ..
          !! processed by numpydoc !!


.. py:class:: SourceGalaxyPopulationModel(z_min=0.0, z_max=10.0, event_type='popI_II')

   
   Class to create lookup tables for redshifts and distances, see :func:`~ler.ler.SourceGalaxyPopulationModel.__init__` method.
















   ..
       !! processed by numpydoc !!
   .. py:method:: create_lookup_table(z_min, z_max)

      
      Functions to create lookup tables for redshifts and comoving volume
      Input parameters:
          z_min (float): minimum redshift of the source population
          z_max (float): maximum redshift of the source population
      Output parameters:
          None
















      ..
          !! processed by numpydoc !!

   .. py:method:: sample_source_redshifts(size=1000, z_min=0.0, z_max=12.0)

      
      Function to sample source redshifts from the source galaxy population model
      Input parameters:
          size (int): number of source redshifts to sample
          z_min (float): minimum redshift of the source population
          z_max (float): maximum redshift of the source population
      Output parameters:
          zs (array): array of source redshifts
















      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_popI_II(zs, R0=23.9 * 1e-09, b2=1.6, b3=2.0, b4=30)

      
      Function to compute the merger rate density (PopI/PopII)
      Input parameters:
          zs (float/array): source redshifts
          R0 (float)      : normalization constant [default: 23.9*1e-9 Mpc^-3 yr^-1]
          b2 (float)      : fitting paramters [default: 1.6]
          b3 (float)      : fitting paramters [default: 2.0]
          b4 (float)      : fitting paramters [default: 30]
      Output parameters:
          rate_density (float/array): merger rate density (Mpc^-3 yr^-1)
















      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_popI_II_Madau_Dickinson(zs, af=2.7, bf=5.6, cf=1.9)

      
      Function to compute the merger rate density (PopI/PopII)
      :param zs:
      :type zs: `float`

      :Returns:

          **rate_density** : `float`
              ..













      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_popIII(zs, aIII=0.66, bIII=0.3, zIII=11.6)

      
      Function to compute the merger rate density (PopIII)
      :param zs:
      :type zs: `float`

      :Returns:

          **rate_density** : `float`
              ..













      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_primordial(zs, t0=13.786885302009708)

      
      Function to compute the merger rate density
















      ..
          !! processed by numpydoc !!


.. py:class:: CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type='popI_II', model_pars=False)

   Bases: :py:obj:`SourceGalaxyPopulationModel`

   
   Class to create lookup tables for redshifts and distances, see :func:`~ler.ler.SourceGalaxyPopulationModel.__init__` method.
















   ..
       !! processed by numpydoc !!
   .. py:method:: sample_gw_parameters(nsamples=1000, **kwargs)

      
      Function to sample BBH parameters from the source galaxy population model
      Input parameters:
          nsamples (int)  : number of BBHs to sample
          **kwargs        : keyword arguments
                          e.g. zs = np.array([0.1,0.2,0.3])
                              model_pars = {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82,                                    'mmin': m_min, 'mmax': m_max, 'lambda_peak': 0.08,                                    'mu_g': 33.07, 'sigma_g': 5.69}
      Output parameters:
          gw_parameters (dict): dictionary of GW parameters
                              e.g. gw_parameters.keys() = ['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source',                                     'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec']
















      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_popI_II(size, model_pars)

      
      Function to calculate source mass1 and mass2 with PowerLaw+PEAK model
















      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_popIII(size, model_pars)

      
      Function to calculate source mass1 and mass2 with pop III origin
















      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_primordial(size, model_pars={'Mc': 30.0, 'sigma': 0.3, 'beta': 1.1})

      
      Function to calculate source mass1 and mass2 with PowerLaw+PEAK model
















      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BNS(size, model_pars)

      
      Function to calculate source mass1 and mass2 with pop III origin
















      ..
          !! processed by numpydoc !!


.. py:function:: solve_lens_equation1(lens_parameters)

   
   Function to solve the lens equation (min_image = 2)
   Input parameters:
       lens_parameters : a list of parameters
                       lens_parameters[0] = e1 : ellipticity
                       lens_parameters[1] = e2 : ellipticity
                       lens_parameters[2] = gamma : power-law index
                       lens_parameters[3] = gamma1 : shear
                       lens_parameters[4] = gamma2 : shear
                       lens_parameters[5] = zl : redshift of the lens
                       lens_parameters[6] = zs : redshift of the source
                       lens_parameters[7] = einstein_radius : Einstein radius
                       lens_parameters[8] = iteration : iteration number
                       lens_parameters[9:] = lens_model_list : list of lens models
   Output parameters:
       x_source : x position of the source in the source plane
       y_source : y position of the source in the source plane
       eta : polar coordinate of the source in the source plane
       phi : polar coordinate of the source in the source plane
       x0_image_position : x position of the images in the source plane
       x1_image_position : y position of the images in the source plane
       magnifications : magnification of the images
       time_delays : time-delay of the images
       nImages : number of images
       determinant : determinant of the hessian matrix
       trace : trace of the hessian matrix
       iteration : iteration number
       weights : weights for the caustic
















   ..
       !! processed by numpydoc !!

.. py:function:: solve_lens_equation2(lens_parameters)

   
   Function to solve the lens equation (min_image > 2)
   Input parameters:
       lens_parameters : a list of parameters
                       lens_parameters[0] = e1 : ellipticity
                       lens_parameters[1] = e2 : ellipticity
                       lens_parameters[2] = gamma : power-law index
                       lens_parameters[3] = gamma1 : shear
                       lens_parameters[4] = gamma2 : shear
                       lens_parameters[5] = zl : redshift of the lens
                       lens_parameters[6] = zs : redshift of the source
                       lens_parameters[7] = einstein_radius : Einstein radius
                       lens_parameters[8] = iteration : iteration number
                       lens_parameters[9:] = lens_model_list : list of lens models
   Output parameters:
       x_source : x position of the source in the source plane
       y_source : y position of the source in the source plane
       eta : polar coordinate of the source in the source plane
       phi : polar coordinate of the source in the source plane
       x0_image_position : x position of the images in the source plane
       x1_image_position : y position of the images in the source plane
       magnifications : magnification of the images
       time_delays : time-delay of the images
       nImages : number of images
       determinant : determinant of the hessian matrix
       trace : trace of the hessian matrix
       iteration : iteration number
       weights : weights for the caustic
















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

