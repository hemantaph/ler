:orphan:

:py:mod:`ler.lens_galaxy_population`
====================================

.. py:module:: ler.lens_galaxy_population


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   lens_galaxy_population/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.CompactBinaryPopulation
   ler.lens_galaxy_population.ImageProperties
   ler.lens_galaxy_population.LensGalaxyPopulation



Functions
~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.add_dictionaries_together
   ler.lens_galaxy_population.trim_dictionary



.. py:class:: CompactBinaryPopulation(z_min=0.0001, z_max=10, event_type='BBH', event_priors=None, event_priors_params=None, cosmology=None)


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

   .. py:method:: event_priors_categorization(event_type, event_priors, event_prior_params)

      
      Function to sample BBH parameters from the source galaxy population
      model


      :Parameters:

          **event_type** : `str`
              Type of event to generate.
              e.g. 'BBH', 'BNS', 'BBH_popIII', 'BBH_primordial', 'NSBH'

          **event_priors** : `dict`
              Dictionary of prior sampler functions for each parameter

          **event_prior_params** : `dict`
              Dictionary of sampler parameters for each parameter

      :Returns:

          **event_priors_** : `dict`
              Dictionary of prior sampler functions for each parameter

          **event_prior_params_** : `dict`
              Dictionary of sampler parameters for each parameter













      ..
          !! processed by numpydoc !!

   .. py:method:: sample_gw_parameters(nsamples=1000, **kwargs)

      
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

   .. py:method:: binary_masses_BBH_popI_II_powerlaw_gaussian(size, mminbh=4.98, mmaxbh=112.5, alpha=3.78, mu_g=32.27, sigma_g=3.88, lambda_peak=0.03, delta_m=4.8, beta=0.81, param=None)

      
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

   .. py:method:: binary_masses_BBH_popIII_gwcosmo(size, Mc=30.0, sigma=0.3, beta=1.1)

      
      Function to calculate source mass1 and mass2 with pop III origin


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **Mc, sigma, beta** : `float`
              Fitting parameters
              default: Mc=30.0, sigma=0.3, beta=1.1

      :Returns:

          **mass_1_source** : `array`
              Array of mass1 in source frame

          **mass_2_source** : `array`
              Array of mass2 in source frame










      .. rubric:: Examples

      >>> from ler import CompactBinaryPopulation
      >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popIII")
      >>> mass_1_source, mass_2_source = pop.binary_masses_popIII(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BBH_primordial_lognormal(size, Mc=30.0, sigma=0.3, beta=1.1)

      
      Function to calculate source mass1 and mass2 for primordial BBHs


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **Mc, sigma, beta** : `float`
              Fitting parameters
              default: Mc=30.0, sigma=0.3, beta=1.1

      :Returns:

          **mass_1_source** : `array`
              Array of mass1 in source frame

          **mass_2_source** : `array`
              Array of mass2 in source frame










      .. rubric:: Examples

      >>> from ler import CompactBinaryPopulation
      >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "primordial")
      >>> mass_1_source, mass_2_source = pop.binary_masses_primordial(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BNS_popI_II_gwcosmo(size, mminns=1.0, mmaxns=3.0, alphans=0.0)

      
      Function to calculate source mass1 and mass2 of BNS (gwcosmo)


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **mminns** : `float`
              Minimum mass of the BNS
              default: 1.0

          **mmaxns** : `float`
              Maximum mass of the BNS
              default: 3.0

          **alphans** : `float`
              Power law index
              default: 0.0

      :Returns:

          **mass_1_source** : `array`
              Array of mass1 in source frame

          **mass_2_source** : `array`
              Array of mass2 in source frame













      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BNS_popI_II_Alsing(size, param=dict(w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1.0, mmax=2.3))

      
      Function to calculate source mass1 and mass2 of BNS (Alsing)


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **w, muL, sigmaL, muR, sigmaR** : `float`
              Fitting parameters
              default: w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3

          **mmin** : `float`
              Minimum mass of the BNS
              default: 1.0

          **mmax** : `float`
              Maximum mass of the BNS
              default: 3.0

      :Returns:

          **mass_1_source** : `array`
              Array of mass1 in source frame

          **mass_2_source** : `array`
              Array of mass2 in source frame













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

   .. py:method:: binary_spin_BBH(size)

      
      Function to calculate spin parameters with PowerLaw+PEAK model


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **a_1** : `array`
              Array of spin1

          **a_2** : `array`
              Array of spin2

          **tilt_1** : `array`
              Array of tilt1

          **tilt_2** : `array`
              Array of tilt2

          **phi_12** : `array`
              Array of phi12

          **phi_jl** : `array`
              Array of phi_jl










      .. rubric:: Examples

      >>> from ler import CompactBinaryPopulation
      >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popI_II")
      >>> a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl = pop.binary_spin_BBH(size=1000)



      ..
          !! processed by numpydoc !!


.. py:class:: ImageProperties(npool=4, n_min_images=2, n_max_images=4, lens_model_list=['EPL_NUMBA', 'SHEAR'])


   
   Class to find the image properties of a lensed event. Image properties include image positions, magnifications, time delays, etc.


   :Parameters:

       **npool** : `int`
           number of processes to use
           default: 4

       **n_min_images** : `int`
           minimum number of images to consider
           default: 2

       **n_max_images** : `int`
           maximum number of images to consider
           default: 4

       **lens_model_list** : `list`
           list of lens models
           default: ['EPL_NUMBA', 'SHEAR']














   ..
       !! processed by numpydoc !!
   .. py:method:: image_properties(lens_parameters)

      
      Function to get the image properties e.g. image positions, magnifications, time delays, etc.


      :Parameters:

          **lens_parameters** : `dict`
              dictionary of lens parameters
              e.g. lens_parameters.keys() = ['zs', 'zl', 'gamma1', 'gamma2', 'e1', 'e2', 'gamma', 'theta_E']

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


.. py:function:: add_dictionaries_together(dictionary1, dictionary2)

   
   Adds two dictionaries with the same keys together.
















   ..
       !! processed by numpydoc !!

.. py:function:: trim_dictionary(dictionary, size)

   
   Filters an event dictionary to only contain the size.
















   ..
       !! processed by numpydoc !!

.. py:class:: LensGalaxyPopulation(CompactBinaryPopulation_=False)


   Bases: :py:obj:`ler.image_properties.ImageProperties`

   
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
              Einstein radii of the lens galaxies in radian













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

   .. py:method:: strong_lensing_optical_depth_SIE(zs)

      
      Function to compute the strong lensing optical depth SIE


      :Parameters:

          **zs** : `float`
              source redshifts

      :Returns:

          **tau** : `float`
              strong lensing optical depth













      ..
          !! processed by numpydoc !!

   .. py:method:: strong_lensing_optical_depth_SIS(zs)

      
      Function to compute the strong lensing optical depth (SIS)


      :Parameters:

          **zs** : `float`
              source redshifts

      :Returns:

          **tau** : `float`
              strong lensing optical depth













      ..
          !! processed by numpydoc !!


