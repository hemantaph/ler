:py:mod:`lens_galaxy_population-checkpoint`
===========================================

.. py:module:: lens_galaxy_population-checkpoint

.. autoapi-nested-parse::

   This module contains the LensGalaxyPopulation class, which is used to sample lens galaxy parameters, source parameters conditioned on the source being strongly lensed, image properties, and lensed SNRs.

   The class inherits from the CompactBinaryPopulation class, which is used to sample source parameters.

   ..
       !! processed by numpydoc !!


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   lens_galaxy_population-checkpoint.LensGalaxyPopulation




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


