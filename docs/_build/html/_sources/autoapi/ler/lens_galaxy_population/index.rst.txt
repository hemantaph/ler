:orphan:

:py:mod:`ler.lens_galaxy_population`
====================================

.. py:module:: ler.lens_galaxy_population


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   lens_galaxy_parameter_distribution/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.CBCSourceParameterDistribution
   ler.lens_galaxy_population.OpticalDepth
   ler.lens_galaxy_population.ImageProperties
   ler.lens_galaxy_population.LensGalaxyParameterDistribution
   ler.lens_galaxy_population.OpticalDepth



Functions
~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.add_dictionaries_together
   ler.lens_galaxy_population.trim_dictionary
   ler.lens_galaxy_population.interpolator_pickle_path
   ler.lens_galaxy_population.phi_cut_SIE
   ler.lens_galaxy_population.phi_q2_ellipticity_hemanta
   ler.lens_galaxy_population.cubic_spline_interpolator
   ler.lens_galaxy_population.inverse_transform_sampler
   ler.lens_galaxy_population.cubic_spline_interpolator2d_array
   ler.lens_galaxy_population.save_pickle
   ler.lens_galaxy_population.interpolator_pickle_path
   ler.lens_galaxy_population.inverse_transform_sampler2d
   ler.lens_galaxy_population.pdf_cubic_spline_interpolator2d_array
   ler.lens_galaxy_population.normal_pdf
   ler.lens_galaxy_population.normal_pdf_2d
   ler.lens_galaxy_population.load_txt_from_module
   ler.lens_galaxy_population.phi_cut_SIE
   ler.lens_galaxy_population.phi
   ler.lens_galaxy_population.phi_loc_bernardi
   ler.lens_galaxy_population.phi_q2_ellipticity_hemanta
   ler.lens_galaxy_population.lens_redshift_SDSS_catalogue_sis
   ler.lens_galaxy_population.phi_cut_SIE
   ler.lens_galaxy_population.phi_q2_ellipticity_hemanta
   ler.lens_galaxy_population.inverse_transform_sampler
   ler.lens_galaxy_population.cubic_spline_interpolator
   ler.lens_galaxy_population.cubic_spline_interpolator2d_array
   ler.lens_galaxy_population.inverse_transform_sampler2d
   ler.lens_galaxy_population.phi_loc_bernardi
   ler.lens_galaxy_population.phi
   ler.lens_galaxy_population.inverse_transform_sampler
   ler.lens_galaxy_population.cubic_spline_interpolator
   ler.lens_galaxy_population.axis_ratio_SIS
   ler.lens_galaxy_population.phi
   ler.lens_galaxy_population.phi_loc_bernardi
   ler.lens_galaxy_population.phi_cut_SIE
   ler.lens_galaxy_population.axis_ratio_rayleigh_rvs
   ler.lens_galaxy_population.velocity_dispersion_z_dependent
   ler.lens_galaxy_population.lens_redshift_SDSS_catalogue_sis
   ler.lens_galaxy_population.bounded_normal_sample
   ler.lens_galaxy_population.phi_q2_ellipticity_hemanta



.. py:class:: CBCSourceParameterDistribution(z_min=0.0, z_max=10.0, event_type='BBH', source_priors=None, source_priors_params=None, cosmology=None, spin_zero=True, spin_precession=False, directory='./interpolator_pickle', create_new_interpolator=False)


   Bases: :py:obj:`ler.gw_source_population.cbc_source_redshift_distribution.CBCSourceRedshiftDistribution`

   
   Class to generate a population of compact binaries. It helps sample all the intrinsic and extrinsic parameters of compact binaries. This daughter class inherits from :class:`~ler.ler.CBCSourceRedshiftDistribution` class.


   :Parameters:

       **z_min** : `float`
           Minimum redshift of the source population
           default: 0.001

       **z_max** : `float`
           Maximum redshift of the source population
           default: 10.

       **event_type** : `str`
           Type of event to generate.
           e.g. 'BBH', 'BNS', 'NSBH'

       **source_priors, source_priors_params** : `dict`, `dict`
           Dictionary of prior sampler functions and its input parameters.
           Check for available priors and corresponding input parameters by running,
           >>> from ler.gw_source_population import CBCSourceParameterDistribution
           >>> cbc = CBCSourceParameterDistribution()
           >>> cbc.available_gw_prior_list_and_its_params()
           # To check the current chosen priors and its parameters, run,
           >>> print("default priors=",cbc.gw_param_samplers)
           >>> print("default priors's parameters=",cbc.gw_param_samplers_params)

       **cosmology** : `astropy.cosmology`
           Cosmology to use
           default: None/astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)

       **spin_zero** : `bool`
           If True, spin parameters are completely ignore in the sampling.
           default: True

       **spin_precession** : `bool`
           If spin_zero=False and spin_precession=True, spin parameters are sampled for precessing binaries.
           if spin_zero=False and spin_precession=False, spin parameters are sampled for aligned/anti-aligned spin binaries.
           default: False

       **directory** : `str`
           Directory to store the interpolator pickle files
           default: './interpolator_pickle'

       **create_new_interpolator** : `dict`
           Dictionary of boolean values and resolution to create new interpolator.
           default: dict(redshift_distribution=dict(create_new=False, resolution=500), z_to_luminosity_distance=dict(create_new=False, resolution=500), differential_comoving_volume=dict(create_new=False, resolution=500))











   .. rubric:: Examples

   >>> from ler.gw_source_population import CBCSourceParameterDistribution
   >>> cbc = CBCSourceParameterDistribution()
   >>> params = cbc.gw_parameters(size=1000)
   >>> print("sampled parameters=",list(params.keys()))

   Instance Attributes
   ----------
   CBCSourceParameterDistribution has the following instance attributes:

   +-------------------------------------+----------------------------------+
   | Atrributes                          | Type                             |
   +=====================================+==================================+
   |:attr:`~z_min`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~z_max`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~event_type`                  | `str`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~source_priors`               | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~source_priors_params`        | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~cosmo`                       | `astropy.cosmology`              |
   +-------------------------------------+----------------------------------+
   |:attr:`~spin_zero`                   | `bool`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~spin_precession`             | `bool`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~directory`                   | `str`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~create_new_interpolator`     | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~available_gw_prior_list_and_its_params`                            |
   +-------------------------------------+----------------------------------+
   |                                     | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~gw_param_samplers`           | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~gw_param_samplers_params`    | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~sampler_names`               | `dict`                           |
   +-------------------------------------+----------------------------------+

   Instance Methods
   ----------
   CBCSourceParameterDistribution has the following instance methods:

   +-------------------------------------+----------------------------------+
   | Methods                             | Type                             |
   +=====================================+==================================+
   |:meth:`~source_priors_categorization`                                   |
   +-------------------------------------+----------------------------------+
   |                                     | Function to categorize the event |
   |                                     | priors and its parameters        |
   +-------------------------------------+----------------------------------+
   |:meth:`~lookup_table_luminosity_distance`                               |
   |                                     | Function to create a lookup      |
   |                                     | table for converting redshift    |
   |                                     | to luminosity distance           |
   +-------------------------------------+----------------------------------+
   |:meth:`~gw_parameters`        | Function to sample all the       |
   |                                     | intrinsic and extrinsic          |
   |                                     | parameters of compact binaries   |
   +-------------------------------------+----------------------------------+
   |:meth:`~source_frame_masses`  | Function to sample source mass1  |
   |                                     | and mass2                        |
   +-------------------------------------+----------------------------------+
   |:meth:`~geocent_time`         | Function to sample geocent time  |
   +-------------------------------------+----------------------------------+
   |:meth:`~zs`                   | Function to sample source        |
   |                                     | redshift                         |
   +-------------------------------------+----------------------------------+
   |:meth:`~ra`                   | Function to sample right         |
   |                                     | ascension (sky position)         |
   +-------------------------------------+----------------------------------+
   |:meth:`~dec`                  | Function to sample declination   |
   |                                     | (sky position)                   |
   +-------------------------------------+----------------------------------+
   |:meth:`~phase`                | Function to sample coalescence   |
   |                                     | phase                            |
   +-------------------------------------+----------------------------------+
   |:meth:`~psi`                  | Function to sample polarization  |
   |                                     | angle                            |
   +-------------------------------------+----------------------------------+
   |:meth:`~theta_jn`             | Function to sample inclination   |
   |                                     | angle                            |
   +-------------------------------------+----------------------------------+
   |:meth:`~a_1`                   | Function to sample spin1         |
   |                                     | magnitude                        |
   +-------------------------------------+----------------------------------+
   |:meth:`~a_2`                   | Function to sample spin2         |
   |                                     | magnitude                        |
   +-------------------------------------+----------------------------------+
   |:meth:`~tilt_1`               | Function to sample tilt1 angle   |
   +-------------------------------------+----------------------------------+
   |:meth:`~tilt_2`               | Function to sample tilt2 angle   |
   +-------------------------------------+----------------------------------+
   |:meth:`~phi_12`               | Function to sample phi12 angle   |
   +-------------------------------------+----------------------------------+
   |:meth:`~phi_jl`               | Function to sample phi_jl angle  |
   +-------------------------------------+----------------------------------+
   |:meth:`~binary_masses_BBH_popI_II_powerlaw_gaussian`                    |
   +-------------------------------------+----------------------------------+
   |                                     | Function to sample source mass1  |
   |                                     | and mass2 with PowerLaw+PEAK     |
   |                                     | model                            |
   +-------------------------------------+----------------------------------+
   |:meth:`~binary_masses_BBH_popIII_lognormal`                             |
   +-------------------------------------+----------------------------------+
   |                                     | Function to sample source mass1  |
   |                                     | and mass2 with popIII orgin from |
   |                                     | lognormal distribution. Refer to |
   |                                     | Ng et al. 2022. Eqn. 1 and 4     |
   +-------------------------------------+----------------------------------+
   |:meth:`~binary_masses_BBH_primordial_lognormal`                         |
   +-------------------------------------+----------------------------------+
   |                                     | Function to sample source mass1  |
   |                                     | and mass2 with primordial orgin  |
   |                                     | from lognormal distribution.     |
   |                                     | Refer to Ng et al. 2022. Eqn. 1  |
   |                                     | and 4                            |
   +-------------------------------------+----------------------------------+
   |:meth:`~binary_masses_BNS_bimodal`   | Function to sample source mass1  |
   |                                     | and mass2 from bimodal           |
   |                                     | distribution. Refer to           |
   |                                     | Will M. Farr et al. 2020 Eqn. 6  |
   +-------------------------------------+----------------------------------+
   |:meth:`~constant_values_n_size`      | Function to return array of      |
   |                                     | constant values of size n        |
   +-------------------------------------+----------------------------------+
   |:meth:`~sampler_uniform`             | Function to sample from uniform  |
   |                                     | distribution                     |
   +-------------------------------------+----------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:property:: source_frame_masses

      
      Function to sample source frame masses (mass1_source, mass2_source) with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **mass_1_source** : `numpy.ndarray` (1D array of floats)
              Array of mass1 in source frame

          **mass_2_source** : `numpy.ndarray` (1D array of floats)
              Array of mass2 in source frame













      ..
          !! processed by numpydoc !!

   .. py:property:: zs

      
      Function to sample source redshift with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **zs** : `numpy.ndarray` (1D array of floats)
              Array of source redshift













      ..
          !! processed by numpydoc !!

   .. py:property:: geocent_time

      
      Function to sample geocent time with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **geocent_time** : `numpy.ndarray` (1D array of floats)
              Array of geocent_time or time of coalescence













      ..
          !! processed by numpydoc !!

   .. py:property:: ra

      
      Function to sample right ascension of sky position with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **ra** : `numpy.ndarray` (1D array of floats)
              Array of right ascension of sky position













      ..
          !! processed by numpydoc !!

   .. py:property:: dec

      
      Function to sample declination of sky position with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **dec** : `numpy.ndarray` (1D array of floats)
              Array of declination of sky position













      ..
          !! processed by numpydoc !!

   .. py:property:: phase

      
      Function to sample coalescence phase with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **phase** : `numpy.ndarray` (1D array of floats)
              Array of coalescence phase













      ..
          !! processed by numpydoc !!

   .. py:property:: psi

      
      Function to sample polarization angle with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **psi** : `numpy.ndarray` (1D array of floats)
              Array of polarization angle













      ..
          !! processed by numpydoc !!

   .. py:property:: theta_jn

      
      Function to sample theta_jn with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **theta_jn** : `numpy.ndarray` (1D array of floats)
              Array of theta_jn













      ..
          !! processed by numpydoc !!

   .. py:property:: a_1

      
      Function to sample spin magnitude of the compact binaries (body1) with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **a_1** : `numpy.ndarray` (1D array of floats)
              Array of spin magnitude of the compact binaries (body1)













      ..
          !! processed by numpydoc !!

   .. py:property:: a_2

      
      Function to sample spin magnitude of the compact binaries (body2) with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **a_2** : `numpy.ndarray` (1D array of floats)
              Array of spin magnitude of the compact binaries (body2)













      ..
          !! processed by numpydoc !!

   .. py:property:: tilt_1

      
      Function to sample tilt angle of the compact binaries (body1) with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **tilt_1** : `numpy.ndarray` (1D array of floats)
              Array of tilt angle of the compact binaries (body1)













      ..
          !! processed by numpydoc !!

   .. py:property:: tilt_2

      
      Function to sample tilt angle of the compact binaries (body2) with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **tilt_2** : `numpy.ndarray` (1D array of floats)
              Array of tilt angle of the compact binaries (body2)













      ..
          !! processed by numpydoc !!

   .. py:property:: phi_12

      
      Function to sample azimuthal angle between the two spins with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **phi_12** : `numpy.ndarray` (1D array of floats)
              Array of azimuthal angle between the two spins













      ..
          !! processed by numpydoc !!

   .. py:property:: phi_jl

      
      Function to sample azimuthal angle between the total angular momentum and the orbital angular momentum with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **phi_jl** : `numpy.ndarray` (1D array of floats)
              Array of azimuthal angle between the total angular momentum and the orbital angular momentum













      ..
          !! processed by numpydoc !!

   .. py:property:: available_gw_prior_list_and_its_params

      
      Dictionary with list all the available priors and it's corresponding parameters. This is an immutable instance attribute.













      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> priors = cbc.available_gw_prior_list_and_its_params
      >>> priors.keys()  # type of priors
      dict_keys(['merger_rate_density', 'source_frame_masses', 'spin', 'geocent_time', 'ra', 'phase', 'psi', 'theta_jn'])
      >>> priors['source_frame_masses'].keys()  # type of source_frame_masses priors
      dict_keys(['binary_masses_BBH_popI_II_powerlaw_gaussian', 'binary_masses_BBH_popIII_lognormal', 'binary_masses_BBH_primordial_lognormal', 'binary_masses_BNS_bimodal'])
      >>> priors['source_frame_masses']['binary_masses_BBH_popI_II_powerlaw_gaussian'].keys()  # parameters of binary_masses_BBH_popI_II_powerlaw_gaussian
      dict_keys(['mminbh', 'mmaxbh', 'alpha', 'mu_g', 'sigma_g', 'lambda_peak', 'delta_m', 'beta'])



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

   .. py:attribute:: source_priors

      
      ``dict``

      Dictionary of prior sampler functions.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: source_priors_params

      
      ``dict``

      Dictionary of prior sampler functions' input parameters.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: cosmo

      
      ``astropy.cosmology``

      Cosmology to use.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: spin_zero

      
      ``bool``

      If True, spin prior is set to zero.















      ..
          !! processed by numpydoc !!

   .. py:method:: setup_decision_dictionary_gw_params(create_new_interpolator)

      
      Method to set up a decision dictionary for interpolator creation.


      :Parameters:

          **create_new_interpolator** : `dict`, `bool`
              If `dict`, dictionary of boolean values and resolution to create new interpolator.
              If `bool`, boolean value to create new interpolator for all quantities.

      :Returns:

          **create_new_interpolator_** : `dict`
              Dictionary of boolean values and resolution to create new interpolator.
              e.g. dict(redshift_distribution=dict(create_new=False, resolution=1000), luminosity_distance=dict(create_new=False, resolution=1000), differential_comoving_volume=dict(create_new=False, resolution=1000))













      ..
          !! processed by numpydoc !!

   .. py:method:: source_priors_categorization(event_type, source_priors, source_prior_params)

      
      Function to categorize the event priors and its parameters.


      :Parameters:

          **event_type** : `str`
              Type of event to generate.
              e.g. 'BBH', 'BNS', 'BBH_popIII', 'BBH_primordial', 'NSBH'

          **source_priors** : `dict`
              Dictionary of prior sampler functions for each parameter

          **source_prior_params** : `dict`
              Dictionary of sampler parameters for each GW parameter

      :Returns:

          **source_priors_** : `dict`
              Dictionary of prior sampler functions for each parameter

          **source_prior_params_** : `dict`
              Dictionary of sampler parameters for each parameter

          **sampler_names_** : `dict`
              Dictionary of sampler names with description










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> source_priors, source_prior_params, sampler_names = cbc.source_priors_categorization(event_type='BBH', source_priors=None, source_prior_params=None)
      >>> print(source_priors.keys())
      >>> print(source_prior_params.keys())
      >>> print(sampler_names.keys())



      ..
          !! processed by numpydoc !!

   .. py:method:: sample_gw_parameters(size=1000, param=None)

      
      Function to sample BBH/BNS/NSBH intrinsic and extrinsics parameters.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **gw_parameters** : `dict`
              Dictionary of sampled parameters
              gw_parameters.keys() = ['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'theta_jn', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl']










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> params = cbc.gw_parameters(size=1000)
      >>> print("sampled parameters=",list(params.keys()))



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BBH_popI_II_powerlaw_gaussian(size, get_attribute=False, **kwargs)

      
      Function to sample source mass1 and mass2 with PowerLaw+PEAK model


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **mminbh** : `float`
              Minimum mass of the black hole (Msun)
              default: 4.98

          **mmaxbh** : `float`
              Maximum mass of the black hole (Msun)
              default: 86.22

          **alpha** : `float`
              Spectral index for the powerlaw of the primary mass distribution
              default: 2.63

          **mu_g** : `float`
              Mean of the Gaussian component in the primary mass distribution
              default: 33.07

          **sigma_g** : `float`
              Width of the Gaussian component in the primary mass distribution
              default: 5.69

          **lambda_peak** : `float`
              Fraction of the model in the Gaussian component
              default: 0.10

          **delta_m** : `float`
              Range of mass tapering on the lower end of the mass distribution
              default: 4.82

          **beta** : `float`
              Spectral index for the powerlaw of the mass ratio distribution

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(mminbh=4.98, mmaxbh=86.22, alpha=2.63, mu_g=33.07, sigma_g=5.69, lambda_peak=0.10, delta_m=4.82, beta=1.26)

      :Returns:

          **mass_1_source** : `numpy.ndarray` (1D array of floats)
              Array of mass1 in source frame (Msun)

          **mass_2_source** : `numpy.ndarray` (1D array of floats)
              Array of mass2 in source frame (Msun)










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> m1_src, m2_src = cbc.binary_masses_BBH_popI_II_powerlaw_gaussian(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BBH_popIII_lognormal(size, get_attribute=False, **kwargs)

      
      Function to sample source mass1 and mass2 with pop III origin. Refer to Eqn. 1 and 4 of Ng et al. 2022


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **m_min** : `float`
              Minimum mass of the black hole (popIII) (Msun)
              default: 10.

          **m_max** : `float`
              Maximum mass of the black hole (popIII) (Msun)
              default: 100.

          **Mc** : `float`
              Mass scale; the distribution is centered around Mc
              default: 30.0

          **sigma** : `float`
              Width of the distribution
              default: 0.3

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(m_min=10., m_max=100., Mc=30.0, sigma=0.3)

      :Returns:

          **mass_1_source** : `numpy.ndarray` (1D array of floats)
              Array of mass1 in source frame (Msun)

          **mass_2_source** : `numpy.ndarray` (1D array of floats)
              Array of mass2 in source frame (Msun)










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> m1_src, m2_src = cbc.binary_masses_BBH_popIII_lognormal(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BBH_primordial_lognormal(size, get_attribute=False, **kwargs)

      
      Function to sample source mass1 and mass2 with primordial origin. Refer to Eqn. 1 and 4 of Ng et al. 2022


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **m_min** : `float`
              Minimum mass of the black hole (primordial) (Msun)
              default: 10.

          **m_max** : `float`
              Maximum mass of the black hole (primordial) (Msun)
              default: 100.

          **Mc, sigma** : `float`
              Fitting parameters
              default: Mc=30.0, sigma=0.3

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(m_min=10., m_max=100., Mc=30.0, sigma=0.3)

      :Returns:

          **mass_1_source** : `numpy.ndarray` (1D array of floats)
              Array of mass1 in source frame (Msun)

          **mass_2_source** : `numpy.ndarray` (1D array of floats)
              Array of mass2 in source frame (Msun)













      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_NSBH_broken_powerlaw(size, get_attribute=False, **kwargs)

      
      Function to calculate source mass1 and mass2 of NSBH from powerlaw distribution (gwcosmo). Parameters are mminbh=26,mmaxbh=125,alpha_1=6.75,alpha_2=6.75,b=0.5,delta_m=5,mminns=1.0,mmaxns=3.0,alphans=0.0.


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **mminbh** : `float`
              Minimum mass of the black hole (Msun)
              default: 26

          **mmaxbh** : `float`
              Maximum mass of the black hole (Msun)
              default: 125

          **alpha_1** : `float`
              Power law index for the primary mass distribution
              default: 6.75

          **alpha_2** : `float`
              Power law index for the secondary mass distribution
              default: 6.75

          **b** : `float`
              Break point of the power law
              default: 0.5

          **delta_m** : `float`
              Range of mass tapering on
              default: 5

          **mminns** : `float`
              Minimum mass of the neutron star (Msun)
              default: 1.0

          **mmaxns** : `float`
              Maximum mass of the neutron star (Msun)
              default: 3.0

          **alphans** : `float`
              Power law index for the neutron star mass distribution
              default: 0.0

          **get_attribute** : `bool`
              If True, return a sampler function with size as the only input where parameters are fixed to the given values.

          **param** : `dict`
              Allows to pass in above parameters as dict.

      :Returns:

          **mass_1_source** : `numpy.ndarray` (1D array of floats)
              Array of mass1 in source frame (Msun)

          **mass_2_source** : `numpy.ndarray` (1D array of floats)
              Array of mass2 in source frame (Msun)










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> m1_src, m2_src = cbc.binary_masses_NSBH_broken_powerlaw(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_uniform(size, get_attribute=False, **kwargs)

      
      Function to sample source mass1 and mass2 from uniform distribution.


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **m_min** : `float`
              Minimum mass of the BNS
              default: 1.0

          **m_max** : `float`
              Maximum mass of the BNS
              default: 3.0

          **get_attribute** : `bool`
              If True, return a sampler function with size as the only input where parameters are fixed to the given values.

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(m_min=1.0, m_max=3.0)

      :Returns:

          **mass_1_source** : `numpy.ndarray` (1D array of floats)
              Array of mass1 in source frame (Msun)

          **mass_2_source** : `numpy.ndarray` (1D array of floats)
              Array of mass2 in source frame (Msun)










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> m1_src, m2_src = cbc.binary_masses_uniform(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BNS_bimodal(size, get_attribute=False, **kwargs)

      
      Function to sample source mass1 and mass2 from bimodal distribution. Refer to Will M. Farr et al. 2020 Eqn. 6, https://arxiv.org/pdf/2005.00032.pdf .


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **w** : `float`
              Weight of the left peak
              default: 0.643

          **muL** : `float`
              Mean of the left peak
              default: 1.352

          **sigmaL** : `float`
              Width of the left peak
              default: 0.08

          **muR** : `float`
              Mean of the right peak
              default: 1.88

          **sigmaR** : `float`
              Width of the right peak
              default: 0.3

          **mmin** : `float`
              Minimum mass of the BNS
              default: 1.0

          **mmax** : `float`
              Maximum mass of the BNS
              default: 2.3

          **resolution** : `int`
              Number of points to sample
              default: 500

          **create_new** : `bool`
              If True, create new interpolator
              default: False

          **get_attribute** : `bool`
              If True, return a sampler function with size as the only input where parameters are fixed to the given values.

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1.0, mmax=2.3, resolution=500)

      :Returns:

          **mass_1_source** : `numpy.ndarray` (1D array of floats)
              Array of mass1 in source frame (Msun)

          **mass_2_source** : `numpy.ndarray` (1D array of floats)
              Array of mass2 in source frame (Msun)










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> m1_src, m2_src = cbc.binary_masses_BNS_bimodal(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: constant_values_n_size(size=100, get_attribute=False, **kwargs)

      
      Function to sample constant values of size n.


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **get_attribute** : `bool`
              If True, return the njitted sampler function with size as the only input where parameters are fixed to the given values.

          **kwargs** : `keyword arguments`
              Additional parameters to pass to the function

      :Returns:

          **values** : `numpy.ndarray` (1D array of floats)
              Array of constant values










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> value = cbc.constant_values_n_size(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: sampler_uniform(size, get_attribute=False, **kwargs)

      
      Function to sample values from uniform distribution.


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **get_attribute** : `bool`
              If True, return the njitted sampler function with size as the only input where parameters are fixed to the given values.

          **param** : `dict`
              Allows to pass in above parameters as dict.

      :Returns:

          **values** : `numpy.ndarray` (1D array of floats)
              Array of uniformly distributed values in the range of [min_, max_]










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> value = cbc.sampler_uniform(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: sampler_cosine(size, get_attribute=False, **kwargs)

      
      Function to sample from sine distribution at the limit of [-np.pi/2, np.pi/2]


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **get_attribute** : `bool`
              If True, return the njitted sampler function with size as the only input where parameters are fixed to the given values.

          **param** : None
              This parameter is not used. It is only here to make the function signature consistent with other samplers.

      :Returns:

          **sine** : `numpy.ndarray` (1D array of floats)
              Array of values in the range of [-np.pi/2, np.pi/2]













      ..
          !! processed by numpydoc !!

   .. py:method:: sampler_sine(size, get_attribute=False, **kwargs)

      
      Function to sample from sine distribution at the limit of [0, np.pi]


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **get_attribute** : `bool`
              If True, return the njitted sampler function with size as the only input where parameters are fixed to the given values.

          **param** : None
              This parameter is not used. It is only here to make the function signature consistent with other samplers.

      :Returns:

          **sine** : `numpy.ndarray` (1D array of floats)
              Array of values in the range of [0, np.pi]













      ..
          !! processed by numpydoc !!


.. py:class:: OpticalDepth(npool=4, z_min=0.0, z_max=10.0, cosmology=None, lens_type='epl_shear_galaxy', lens_functions=None, lens_functions_params=None, lens_param_samplers=None, lens_param_samplers_params=None, directory='./interpolator_pickle', create_new_interpolator=False, verbose=False)


   
   Class to calculate the optical depth, velocity dispersion and axis-ratio of a lens galaxy population.


   :Parameters:

       **npool** : int, optional
           Number of processors to use for multiprocessing (default is 4).

       **z_min** : float, optional
           Minimum redshift of the lens galaxy population (default is 0.0).

       **z_max** : float, optional
           Maximum redshift of the lens galaxy population (default is 10.0).

       **cosmology** : astropy.cosmology, optional
           Cosmology object to use (default is FlatLambdaCDM with H0=70, Om0=0.3, Ode0=0.7).

       **lens_type** : str, optional
           Type of the lens galaxy. Must be one of ['sie_galaxy', 'epl_shear_galaxy', 'sis_galaxy'] (default is 'epl_shear_galaxy').

       **lens_functions** : dict, optional
           Dictionary with lens-related functions.

       **lens_functions_params** : dict, optional
           Dictionary with parameters for the lens-related functions.

       **lens_param_samplers** : dict, optional
           Dictionary of sampler functions for velocity dispersion and axis-ratio.

       **lens_param_samplers_params** : dict, optional
           Dictionary with parameters for the priors of the samplers.

       **directory** : str, optional
           Directory where the interpolators are saved (default is './interpolator_pickle').
           If True, creates a new interpolator (default is False).

       **verbose** : bool, optional
           If True, prints additional information during initialization (default is False).





   :Raises:

       ValueError
           If `lens_type` is not in ['sie_galaxy', 'epl_shear_galaxy', 'sis_galaxy'].









   ..
       !! processed by numpydoc !!
   .. py:property:: optical_depth

      
      Function to compute the strong lensing optical depth.


      :Parameters:

          **zs** : `numpy.ndarray` (1D array of floats)
              source redshifts

      :Returns:

          **tau** : `numpy.ndarray` (1D array of floats)
              strong lensing optical depth










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> print(od.optical_depth(np.array([0.1,0.2,0.3])))



      ..
          !! processed by numpydoc !!

   .. py:property:: velocity_dispersion

      
      Class object to sample velocity dispersion. `zl` is required only if velocity dispersion sampler is redshift dependent.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **zl** : `float`
              redshift of the lens galaxy

      :Returns:

          **sigma** : `numpy.ndarray` (1D array of floats)
              velocity dispersion of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> print(od.velocity_dispersion(size=10))



      ..
          !! processed by numpydoc !!

   .. py:property:: axis_ratio

      
      Function to sample axis ratio from rayleigh distribution with given velocity dispersion.


      :Parameters:

          **sigma** : `numpy.ndarray` (1D array of floats)
              velocity dispersion of the lens galaxy

      :Returns:

          **q** : `numpy.ndarray` (1D array of floats)
              axis ratio of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> print(od.axis_ratio(sigma=200.))



      ..
          !! processed by numpydoc !!

   .. py:property:: available_lens_prior_list_and_its_params

      
      Dictionary with list all the available priors and it's corresponding parameters. This is an immutable instance attribute.
















      ..
          !! processed by numpydoc !!

   .. py:property:: available_lens_functions_and_its_params

      
      Dictionary with list all the available lens functions. This is an immutable instance attribute.
















      ..
          !! processed by numpydoc !!

   .. py:method:: default_lens_samplers_and_functions(lens_type)

      
      Function to categorize the lens priors/samplers


      :Parameters:

          **lens_type** : `str`
              lens type
              e.g. 'epl_shear_galaxy' for elliptical power-law galaxy

      :Returns:

          **lens_priors_** : `dict`
              dictionary of priors

          **lens_priors_params_** : `dict`
              dictionary of priors parameters

          **lens_sampler_names_** : `dict`
              dictionary of sampler names

          **lens_functions_** : `dict`
              dictionary of lens functions













      ..
          !! processed by numpydoc !!

   .. py:method:: initialize_decision_dictionary(create_new_interpolator)

      
      Function to initialize decision dictionary for creating interpolator


      :Parameters:

          **create_new_interpolator** : `dict` or `bool`
              dictionary to create new interpolator for velocity dispersion and optical depth.














      ..
          !! processed by numpydoc !!

   .. py:method:: lens_functions_and_sampler_categorization(lens_param_samplers, lens_param_samplers_params, lens_functions, lens_functions_params)

      
      Function to initialize velocity dispersion sampler with it's settings. The reason I am seperating this from lens_priors_categorization is only a specific parameters needs special attention.


      :Parameters:

          **lens_param_samplers** : `str` or `function`
              sampler name or function

          **lens_param_samplers_params** : `dict`
              sampler parameters

          **lens_functions** : `str` or `function`
              lens function name or function

          **lens_functions_params** : `dict`
              lens function parameters














      ..
          !! processed by numpydoc !!

   .. py:method:: axis_ratio_rayleigh(size, sigma, get_attribute=False, **kwargs)

      
      Function to sample axis ratio from rayleigh distribution with given velocity dispersion.


      :Parameters:

          **sigma** : `float: array`
              velocity dispersion of the lens galaxy

          **q_min, q_max** : `float`
              minimum and maximum axis ratio

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample axis ratio

      :Returns:

          **q** : `float: array`
              axis ratio of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(axis_ratio="axis_ratio_rayleigh"))
      >>> print(od.axis_ratio(sigma=200.))



      ..
          !! processed by numpydoc !!

   .. py:method:: axis_ratio_padilla_strauss(size=1000, get_attribute=False, **kwargs)

      
      Function to sample axis ratio using Padilla and Strauss 2008 distribution for axis ratio


      :Parameters:

          **size** : `int`
              sample size

          **q_min, q_max** : `float`
              minimum and maximum axis ratio

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample axis ratio

      :Returns:

          **q** : `float: array`
              axis ratio of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(axis_ratio="axis_ratio_padilla_strauss"))
      >>> print(od.axis_ratio(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: lens_redshift_SDSS_catalogue_numerical(size=1000, zs=None, get_attribute=False, **kwargs)

      
      Function to sample lens redshifts, conditioned on the lens being strongly lensed


      :Parameters:

          **size** : `int`
              sample size

          **zs** : `float`
              source redshifts

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample lens redshifts

      :Returns:

          **zs** : `float: array`
              lens redshifts










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(lens_redshift="lens_redshift_SDSS_catalogue_numerical"))
      >>> print(od.lens_redshift(size=10, zs=1.0))



      ..
          !! processed by numpydoc !!

   .. py:method:: lens_redshift_SDSS_catalogue_hemanta(size=1000, zs=None, get_attribute=False, **kwargs)

      
      Function to sample lens redshifts, conditioned on the lens being strongly lensed


      :Parameters:

          **size** : `int`
              sample size

          **zs** : `float`
              source redshifts

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample lens redshifts

      :Returns:

          **zs** : `float: array`
              lens redshifts










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(lens_redshift="lens_redshift_SDSS_catalogue_numerical"))
      >>> print(od.lens_redshift(size=10, zs=1.0))



      ..
          !! processed by numpydoc !!

   .. py:method:: intrinsic_lens_redshift(size=1000, get_attribute=False, **kwargs)

      
      Function to sample intrinsic lens redshifts, based on the intrinsic velocity dispersion of the lens galaxy.


      :Parameters:

          **size** : `int`
              sample size

          **zs** : `float`
              source redshifts

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample lens redshifts

      :Returns:

          **zs** : `float: array`
              lens redshifts













      ..
          !! processed by numpydoc !!

   .. py:method:: lens_redshift_multiprocessing(zl_scaled2d, zs1d)

      
      Compute the lens redshift distribution using multiprocessing.


      :Parameters:

          **zl_scaled2d** : array_like
              2D array of lens redshifts, scaled by the source redshift.

          **zs1d** : array_like
              1D array of source redshifts.

          **zl_distribution_name** : str
              Name of the lens redshift distribution to compute.

      :Returns:

          **density_array** : array_like
              2D array of the lens redshift distribution.













      ..
          !! processed by numpydoc !!

   .. py:method:: axis_rotation_angle_uniform(size, get_attribute=False, **kwargs)

      
      Function to sample the axis rotation angle of the elliptical lens galaxy from a uniform distribution.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **get_attribute** : `bool`
              if True, returns a function that can be called with size as input

      :Returns:

          **phi** : `numpy.ndarray`
              axis rotation angle of the elliptical lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(axis_rotation_angle="axis_rotation_angle_uniform"))
      >>> print(od.axis_rotation_angle_uniform(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: axis_ratio_uniform(size, get_attribute=False, **kwargs)

      
      Function to sample the axis ratio of the elliptical lens galaxy from a uniform distribution.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **get_attribute** : `bool`
              if True, returns a function that can be called with size as input

      :Returns:

          **q** : `numpy.ndarray`
              axis ratio of the elliptical lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(axis_ratio="axis_ratio_uniform"))
      >>> print(od.axis_ratio_uniform(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: external_shear_normal(size, get_attribute=False, **kwargs)

      
      Function to sample the external shear parameters from a normal distribution.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **get_attribute** : `bool`
              if True, returns a function that can be called with size as input

      :Returns:

          **gamma_1** : `numpy.ndarray`
              shear component in the x-direction

          **gamma_2** : `numpy.ndarray`
              shear component in the y-direction










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(external_shear="external_shear_normal"))
      >>> print(od.external_shear_normal(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: external_shear_numerical_hemanta(size, get_attribute=False, **kwargs)

      
      Function to sample the external shear parameters from a normal distribution.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **get_attribute** : `bool`
              if True, returns a function that can be called with size as input

      :Returns:

          **gamma_1** : `numpy.ndarray`
              shear component in the x-direction

          **gamma_2** : `numpy.ndarray`
              shear component in the y-direction










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(external_shear="external_shear_normal"))
      >>> print(od.external_shear_normal(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: density_profile_slope_normal(size, get_attribute=False, **kwargs)

      
      Function to sample the lens galaxy density profile slope with normal distribution.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample velocity dispersion

          **\*\*kwargs** : `dict`
              additional parameters to be passed to the function,
              e.g. `mean` and `std` for the normal distribution

      :Returns:

          **slope** : `float`
              density profile slope of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(density_profile_slope="density_profile_slope_normal"))
      >>> print(od.density_profile_slope_normal(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: density_profile_slope_numerical_hemanta(size, get_attribute=False, **kwargs)

      
      Function to sample the lens galaxy density profile slope with normal distribution.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample velocity dispersion

          **\*\*kwargs** : `dict`
              additional parameters to be passed to the function,
              e.g. `mean` and `std` for the normal distribution

      :Returns:

          **slope** : `float`
              density profile slope of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(density_profile_slope="density_profile_slope_normal"))
      >>> print(od.density_profile_slope_normal(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: lens_redshift_SDSS_catalogue_sis(size, zs, get_attribute=False, **kwargs)

      
      Function to sample lens redshifts, conditioned on the lens being strongly lensed


      :Parameters:

          **zs** : `float`
              source redshifts

          **get_attribute** : `bool`
              If True, returns a function that can be called with zs as input

      :Returns:

          **zl** : `float`
              lens redshifts










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.lens_redshift_SDSS_catalogue_sis(zs=1.0)



      ..
          !! processed by numpydoc !!

   .. py:method:: velocity_dispersion_gengamma(size, get_attribute=False, **kwargs)

      
      Function to sample velocity dispersion from gengamma distribution


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **a,c** : `float`
              parameters of gengamma distribution
              refer to https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gengamma.html

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample velocity dispersion

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(a=2.32 / 2.67, c=2.67)

      :Returns:

          **sigma** : `numpy.ndarray` (1D array of floats)
              velocity dispersion of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(velocity_dispersion="velocity_dispersion_gengamma"), lens_param_samplers_params=dict(velocity_dispersion=dict(a=2.32 / 2.67, c=2.67)))
      >>> print(od.velocity_dispersion(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: velocity_dispersion_bernardi(size, get_attribute=False, **kwargs)

      
      Function to sample velocity dispersion from Bernardi et al. (2010). This uses inverse transform sampling.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample velocity dispersion

      :Returns:

          **sigma** : `numpy.ndarray` (1D array of floats)
              velocity dispersion of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(velocity_dispersion="velocity_dispersion_bernardi"))
      >>> print(od.velocity_dispersion(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: velocity_dispersion_ewoud(size, zl, get_attribute=False, **kwargs)

      
      Function to sample velocity dispersion (redshift dependent) from Wempe et al. (2022). This uses inverse transform sampling.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **zl** : `float`
              redshift of the lens galaxy

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample velocity dispersion

      :Returns:

          **sigma** : `numpy.ndarray` (1D array of floats)
              velocity dispersion of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(velocity_dispersion="velocity_dispersion_ewoud"))
      >>> print(od.velocity_dispersion(size=10, zl=0.5))



      ..
          !! processed by numpydoc !!

   .. py:method:: cross_section_sis(sigma, zl, zs, **kwargs)

      
      Function to compute the SIS cross-section


      :Parameters:

          **sigma** : `float`
              velocity dispersion of the lens galaxy

          **zl** : `float`
              redshift of the lens galaxy

          **zs** : `float`
              redshift of the source galaxy

      :Returns:

          **cross_section** : `float`
              SIS cross-section










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> print(od.cross_section_sis(sigma=200., zl=0.5, zs=1.0))



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
              Einstein radii of the lens galaxies in radians. Multiply by










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> sigma = 200.0
      >>> zl = 0.5
      >>> zs = 1.0
      >>> lens.compute_einstein_radii(sigma, zl, zs)



      ..
          !! processed by numpydoc !!

   .. py:method:: cross_section_sie_feixu(sigma, zl, zs, q, **kwargs)

      
      Function to compute the SIE cross-section from Fei Xu et al. (2021)


      :Parameters:

          **sigma** : `float`
              velocity dispersion of the lens galaxy

          **zl** : `float`
              redshift of the lens galaxy

          **zs** : `float`
              redshift of the source galaxy

      :Returns:

          **cross_section** : `float`
              SIE cross-section










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> print(od.cross_section_sie_feixu(sigma=200., zl=0.5, zs=1.0, q=1.0))



      ..
          !! processed by numpydoc !!

   .. py:method:: optical_depth_epl_shear_hemanta(zs, get_attribute=False, **kwargs)

      
      Function to compute the strong lensing optical depth (EPL with shear).

      LambdaCDM(H0=70, Om0=0.3, Ode0=0.7) was used to derive the following equation. This is the analytic version of optical depth from z=0 to z=zs.

      :Parameters:

          **zs** : `float`
              source redshifts

      :Returns:

          **tau** : `float`
              strong lensing optical depth










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> print(od.optical_depth_epl_shear_lambdacdm(zs=1.0))



      ..
          !! processed by numpydoc !!

   .. py:method:: optical_depth_sis_haris(zs, get_attribute=False, **kwargs)

      
      Function to compute the strong lensing optical depth (SIS).

      LambdaCDM(H0=70, Om0=0.3, Ode0=0.7) was used to derive the following equation. This is the analytic version of optical depth from z=0 to z=zs.

      :Parameters:

          **zs** : `float`
              source redshifts

      :Returns:

          **tau** : `float`
              strong lensing optical depth










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> print(od.optical_depth_sis_haris(zs=1.0))



      ..
          !! processed by numpydoc !!

   .. py:method:: interpolated_cross_section_function(theta_E, e1, e2, gamma, gamma1, gamma2, get_attribute=False, **kwargs)

      
      Function to compute the cross-section correction factor
















      ..
          !! processed by numpydoc !!

   .. py:method:: create_lookup_table_fuction()

      
      Functions to create lookup tables
      1. Redshift to co-moving distance.
      2. Co-moving distance to redshift.
      3. Redshift to angular diameter distance.
















      ..
          !! processed by numpydoc !!


.. py:class:: ImageProperties(npool=4, z_min=0.0, z_max=10, n_min_images=2, n_max_images=4, geocent_time_min=1126259462.4, geocent_time_max=1126259462.4 + 365 * 24 * 3600 * 20, lens_model_list=['EPL_NUMBA', 'SHEAR'], cosmology=None, spin_zero=True, spin_precession=False, directory='./interpolator_pickle', create_new_interpolator=False)


   
   Class to find the image properties of a lensed event. Image properties include image positions, magnifications, time delays, etc.


   :Parameters:

       **npool** : `int`
           number of processes to use
           default: 4

       **z_min** : `float`
           minimum redshift to consider
           default: 0.0

       **z_max** : `float`
           maximum redshift to consider
           default: 10.0

       **n_min_images** : `int`
           minimum number of images to consider
           default: 2

       **n_max_images** : `int`
           maximum number of images to consider
           default: 4

       **geocent_time_min** : `float`
           minimum geocent time to consider
           default: 1126259462.4 , which is the GPS time of the first GW detection

       **geocent_time_max** : `float`
           maximum geocent time to consider
           default: 1126259462.4+365*24*3600*100 , which is the GPS time of the first GW detection + 100 years. Some time delays can be very large.

       **lens_model_list** : `list`
           list of lens models
           default: ['EPL_NUMBA', 'SHEAR']

       **cosmology** : `astropy.cosmology`
           cosmology
           default: None/astropy.cosmology.LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

       **spin_zero** : `bool`
           whether to assume spin zero or not
           default: True

       **spin_precession** : `bool`
           whether to assume spin precession or not
           default: False

       **directory** : `str`
           directory to save the interpolator pickle files
           default: "./interpolator_pickle"

       **create_new_interpolator** : `dict`
           dictionary to create new interpolator pickle files
           default: dict(luminosity_distance_to_z=dict(create_new=False, resolution=1000))











   .. rubric:: Examples

   >>> from ler.image_properties import ImageProperties
   >>> image_properties = ImageProperties()
   >>> lens_parameters = dict(zs=2.0, zl=0.5, gamma1=0.0, gamma2=0.0, e1=0.0, e2=0.0, gamma=2.0, theta_E=1.0)
   >>> lens_parameters = image_properties.image_properties(lens_parameters)
   >>> print(lens_parameters.keys())

   Instance Attributes
   ----------
   ImageProperties has the following instance attributes:

   +-------------------------------------+----------------------------------+
   | Atrributes                          | Type                             |
   +=====================================+==================================+
   |:attr:`npool`                        | `int`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`z_min`                        | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`z_max`                        | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`n_min_images`                 | `int`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`n_max_images`                 | `int`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`geocent_time_min`             | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`geocent_time_max`             | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`lens_model_list`              | `list`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`cosmo`                        | `astropy.cosmology`              |
   +-------------------------------------+----------------------------------+
   |:attr:`spin_zero`                    | `bool`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`spin_precession`              | `bool`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`directory`                    | `str`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`create_new_interpolator`      | `dict`                           |
   +-------------------------------------+----------------------------------+



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

              source related=>['mass_1': mass in detector frame (mass1>mass2), 'mass_2': mass in detector frame, 'mass_1_source':mass in source frame, 'mass_2_source':mass source frame, 'luminosity_distance': luminosity distance, 'theta_jn': inclination angle, 'psi': polarization angle, 'phase': coalesence phase, 'geocent_time': coalensence GPS time at geocenter, 'ra': right ascension, 'dec': declination, 'a_1': spin magnitude of the more massive black hole, 'a2': spin magnitude of the less massive black hole, 'tilt_1': tilt angle of the more massive black hole, 'tilt_2': tilt angle of the less massive black hole, 'phi_12': azimuthal angle between the two spins, 'phi_jl': azimuthal angle between the total angular momentum and the orbital angular momentum]

              image related=>['x_source': source position in the x-direction, 'y_source': source position in the y-direction, 'x0_image_position': image position in the x-direction, 'x1_image_position': image position in the y-direction, 'magnifications': magnifications, 'time_delays': time delays, 'n_images': number of images formed, 'determinant': determinants, 'trace': traces, 'iteration': to keep track of the iteration number













      ..
          !! processed by numpydoc !!

   .. py:method:: get_lensed_snrs(lensed_param, list_of_detectors=None, snr_calculator=None, pdet_calculator=None)

      
      Function to calculate the signal to noise ratio for each image in each event.


      :Parameters:

          **snr_calculator** : `function`
              snr function, as describe in the :class:`~ler.rates.GWRATES` class.

          **list_of_detectors** : `list`
              list of detectors
              e.g. ['H1', 'L1', 'V1']

          **lensed_param** : `dict`
              dictionary containing the both already lensed source paramters and image parameters.
              e.g. lensed_param.keys() = ['mass_1', 'mass_2', 'zs', 'luminosity_distance', 'theta_jn', 'psi', 'phi', 'ra', 'dec', 'geocent_time', 'phase', 'a_1', 'a2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl', 'magnifications', 'time_delays']

          **n_max_images** : `int`
              maximum number of images to consider
              default: 4

      :Returns:

          **snrs** : `dict`
              signal to noise ratio for each image in each event.
              (dictionary containing 'H1', 'L1', ..., and 'snr_net', which is the network snr, for each image as an array with dimensions (number_of_lensed_events,n_max_images) )













      ..
          !! processed by numpydoc !!


.. py:function:: add_dictionaries_together(dictionary1, dictionary2)

   
   Adds two dictionaries with the same keys together.


   :Parameters:

       **dictionary1** : `dict`
           dictionary to be added.

       **dictionary2** : `dict`
           dictionary to be added.

   :Returns:

       **dictionary** : `dict`
           dictionary with added values.













   ..
       !! processed by numpydoc !!

.. py:function:: trim_dictionary(dictionary, size)

   
   Filters an event dictionary to only contain the size.


   :Parameters:

       **dictionary** : `dict`
           dictionary to be trimmed.

       **size** : `int`
           size to trim the dictionary to.

   :Returns:

       **dictionary** : `dict`
           trimmed dictionary.













   ..
       !! processed by numpydoc !!

.. py:function:: interpolator_pickle_path(param_dict_given, directory, sub_directory, interpolator_name)

   
   Function to create the interpolator pickle file path.


   :Parameters:

       **param_dict_given** : `dict`
           dictionary of parameters.

       **directory** : `str`
           directory to store the interpolator.

       **sub_directory** : `str`
           sub-directory to store the interpolator.

       **interpolator_name** : `str`
           name of the interpolator.

   :Returns:

       **path_inv_cdf** : `str`
           path of the interpolator pickle file.

       **it_exist** : `bool`
           if True, the interpolator exists.













   ..
       !! processed by numpydoc !!

.. py:function:: phi_cut_SIE(q)

   
   Function to calculate cross-section scaling factor for the SIE lens galaxy from SIS lens galaxy.


   :Parameters:

       **q** : `float: array`
           axis ratio of the lens galaxy

   :Returns:

       **result** : `float: array`
           scaling factor













   ..
       !! processed by numpydoc !!

.. py:function:: phi_q2_ellipticity_hemanta(phi, q)

   
   Function to convert phi and q to ellipticity e1 and e2.


   :Parameters:

       **phi** : `float: array`
           angle of the major axis in radians

       **q** : `float: array`
           axis ratio

   :Returns:

       **e1** : `float: array`
           ellipticity component 1

       **e2** : `float: array`
           ..













   ..
       !! processed by numpydoc !!

.. py:class:: LensGalaxyParameterDistribution(npool=4, z_min=0.0, z_max=10.0, cosmology=None, event_type='BBH', lens_type='epl_shear_galaxy', lens_functions=None, lens_functions_params=None, lens_param_samplers=None, lens_param_samplers_params=None, directory='./interpolator_pickle', create_new_interpolator=False, buffer_size=1000, **kwargs)


   Bases: :py:obj:`ler.gw_source_population.CBCSourceParameterDistribution`, :py:obj:`ler.image_properties.ImageProperties`, :py:obj:`ler.lens_galaxy_population.optical_depth.OpticalDepth`

   
   Class to sample lens galaxy parameters, source parameters conditioned on the source being strongly lensed, and image properties


   :Parameters:

       **npool** : `int`
           number of processors to use

       **z_min** : `float`
           minimum redshift

       **z_max** : `float`
           maximum redshift

       **cosmology** : `astropy.cosmology`
           Cosmology to use
           default: None/astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)

       **event_type** : `str`
           Type of event to generate.
           e.g. 'BBH', 'BNS', 'NSBH'
           default: 'BBH'

       **lens_type** : `str`
           Type of lens galaxy to generate.
           default: 'epl_shear_galaxy'

       **lens_functions, lens_priors, lens_priors_params** : `dict`, `dict`, `dict`
           dictionary of lens functions, priors, and priors parameters
           Check for default/available lens functions, priors and corresponding input parameters by running,

           >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
           >>> lens = LensGalaxyParameterDistribution()
           >>> print(lens.lens_functions)
           >>> print(lens.lens_priors)
           >>> print(lens.lens_priors_params)

       **directory** : `str`
           directory to store the interpolators
           default: './interpolator_pickle'

       **\*\*kwargs**
           keyword arguments to pass to the parent classes











   .. rubric:: Examples

   >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
   >>> lens = LensGalaxyParameterDistribution()
   >>> lensed_params = lens.sample_lens_parameters(size=1000)
   >>> lensed_params.keys()

   Instance Attributes
   ----------
   LensGalaxyPopulation class has the following instance attributes:

   +-------------------------------------+----------------------------------+
   | Atrributes                          | Type                             |
   +=====================================+==================================+
   |:attr:`~npool`                       | `int`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~z_min`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~z_max`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~cosmo`                       | `astropy.cosmology`              |
   +-------------------------------------+----------------------------------+
   |:attr:`~event_type`                  | `str`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~directory`                   | `str`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~create_new_interpolator`     | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~lens_param_samplers`         | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~lens_param_samplers_params`  | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~lens_sampler_names`          | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~lens_functions`              | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~normalization_pdf_z_lensed`  | `float`                          |
   +-------------------------------------+----------------------------------+

   Instance Methods
   ----------
   LensGalaxyPopulation class has the following instance methods:

   +-------------------------------------+----------------------------------+
   | Methods                             | Type                             |
   +=====================================+==================================+
   |:meth:`~sample_lens_parameters`      | Function to call the specific    |
   |                                     | galaxy lens parameters sampler   |
   |                                     | routine.                         |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_all_routine_sie`      | Function to sample galaxy lens   |
   |                                     | parameters along with the source |
   |                                     | parameters.                      |
   +-------------------------------------+----------------------------------+
   |:meth:`~strongly_lensed_source_redshifts`                               |
   +-------------------------------------+----------------------------------+
   |                                     | Function to sample source        |
   |                                     | redshifts conditioned on the     |
   |                                     | source being strongly lensed     |
   +-------------------------------------+----------------------------------+
   |:meth:`~source_parameters`           | Function to sample gw source     |
   |                                     | parameters                       |
   +-------------------------------------+----------------------------------+
   |:meth:`~lens_redshift_SDSS_catalogue`| Function to sample lens          |
   |                                     | redshifts, conditioned on the    |
   |                                     | lens being strongly lensed       |
   +-------------------------------------+----------------------------------+
   |:meth:`~axis_rotation_angle_uniform` | Function to sample the axis      |
   |                                     | rotation angle of the elliptical |
   |                                     | lens galaxy from a uniform       |
   |                                     | distribution                     |
   +-------------------------------------+----------------------------------+
   |:meth:`~shear_norm`                  | Function to sample the           |
   |                                     | elliptical lens galaxy shear     |
   |                                     | from a normal distribution       |
   +-------------------------------------+----------------------------------+
   |:meth:`~density_profile_slope_normal`                             |
   +-------------------------------------+----------------------------------+
   |                                     | Function to sample the lens      |
   |                                     | galaxy spectral index of the     |
   |                                     | mass density profile from a      |
   |                                     | normal distribution              |
   +-------------------------------------+----------------------------------+
   |:meth:`~compute_einstein_radii`      | Function to compute the Einstein |
   |                                     | radii of the lens galaxies       |
   +-------------------------------------+----------------------------------+
   |:meth:`~rjs_with_cross_section_sis`  | Function to conduct rejection    |
   |                                     | sampling wrt einstein radius     |
   +-------------------------------------+----------------------------------+
   |:meth:`~rjs_with_cross_section_sie`  | Function to conduct rejection    |
   |                                     | sampling wrt cross_section       |
   +-------------------------------------+----------------------------------+
   |:attr:`~rejection_sample_sl`         | Function to conduct rejection    |
   |                                     | sampling with the given rejection|
   |                                     | sampling function                |
   +-------------------------------------+----------------------------------+
   |:attr:`~sample_source_redshift_sl`   | Function to sample source        |
   |                                     | redshifts conditioned on the     |
   |                                     | source being strongly lensed     |
   +-------------------------------------+----------------------------------+
   |:attr:`~sample_lens_redshift`        | Function to sample lens          |
   |                                     | redshifts, conditioned on the    |
   |                                     | lens being strongly lensed       |
   +-------------------------------------+----------------------------------+
   |:attr:`~sample_axis_rotation_angle`  | Function to sample the axis      |
   |                                     | rotation angle of the elliptical |
   |                                     | lens galaxy from a uniform       |
   |                                     | distribution                     |
   +-------------------------------------+----------------------------------+
   |:attr:`~sample_shear`                | Function to sample the           |
   |                                     | elliptical lens galaxy shear     |
   |                                     | from a normal distribution       |
   +-------------------------------------+----------------------------------+
   |:attr:`~sample_density_profile_slope`                             |
   +-------------------------------------+----------------------------------+
   |                                     | Function to sample the lens      |
   |                                     | galaxy spectral index of the     |
   |                                     | mass density profile from a      |
   |                                     | normal distribution              |
   +-------------------------------------+----------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:attribute:: cbc_pop

      
      :class:`~CBCSourceParameterDistribution` class

      This is an already initialized class that contains a function (CBCSourceParameterDistribution.sample_gw_parameters) that actually samples the source parameters.















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

   .. py:method:: class_initialization_lens(npool, z_min, z_max, cosmology, lens_type, lens_functions, lens_functions_params, lens_param_samplers, lens_param_samplers_params, directory, create_new_interpolator, params)

      
      Initialize the LensGalaxyParameterDistribution class.


      :Parameters:

          **npool** : `int`
              number of processors to use for sampling

          **z_min** : `float`
              minimum redshift of the lens galaxy

          **z_max** : `float`
              maximum redshift of the lens galaxy

          **cosmology** : `astropy.cosmology`
              cosmology object

          **lens_type** : `str`
              type of the lens galaxy

          **lens_functions** : `dict`
              dictionary with the lens related functions

          **lens_functions_params** : `dict`
              dictionary with the parameters for the lens related functions

          **lens_param_samplers** : `dict`
              dictionary with the priors for the sampler

          **lens_param_samplers_params** : `dict`
              dictionary with the parameters for the priors of the sampler

          **directory** : `str`
              directory where the interpolators are saved

          **create_new_interpolator** : `bool`
              if True, creates a new interpolator

          **params** : `dict`
              additional parameters for the CBCSourceParameterDistribution and ImageProperties classes














      ..
          !! processed by numpydoc !!

   .. py:method:: sample_lens_parameters(size=1000)

      
      Function to sample galaxy lens parameters along with the source parameters, conditioned on the source being strongly lensed.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **lens_parameters** : `dict`
              dictionary of sampled lens parameters and source parameters.

              keys:

              zl: lens redshifts

              zs: source redshifts, lensed condition applied

              sigma: velocity dispersions

              q: axis ratios

              theta_E: Einstein radii

              phi: axis rotation angle

              e1: ellipticity component 1

              e2: ellipticity component 2

              gamma1: shear component 1

              gamma2: shear component 2

              gamma: density profile slope distribution

              geocent_time: time of arrival of the unlensed signal

              phase: phase of the unlensed signal

              psi: polarization angle of the unlensed signal

              theta_jn: inclination angle of the unlensed signal

              luminosity_distance: luminosity distance of the source

              mass_1_source: mass 1 (larger) of the source

              mass_2_source: mass 2 (smaller) of the source

              ra: right ascension of the source

              dec: declination of the source










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> od = LensGalaxyParameterDistribution(lens_param_samplers=dict(velocity_dispersion="velocity_dispersion_ewoud"))
      >>> print(od.sample_lens_parameters(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: sample_all_routine_sie_sl(size=1000)

      
      Function to sample galaxy lens parameters. SIE cross section is used for rejection sampling.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **lens_parameters_input** : `dict`
              dictionary of lens parameters to sample

      :Returns:

          **lens_parameters** : `dict`
              dictionary of lens parameters and source parameters (lens conditions applied):

              zl: lens redshifts

              zs: source redshifts, lensed condition applied

              sigma: velocity dispersions

              q: axis ratios

              theta_E: Einstein radii

              phi: axis rotation angle

              e1: ellipticity component 1

              e2: ellipticity component 2

              gamma1: shear component 1

              gamma2: shear component 2

              gamma: density profile slope distribution










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.sample_all_routine_sie(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: sample_all_routine_epl_shear_sl(size=1000)

      
      Function to sample galaxy lens parameters along. EPL shear cross section is used for rejection sampling.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **lens_parameters_input** : `dict`
              dictionary of lens parameters to sample

      :Returns:

          **lens_parameters** : `dict`
              dictionary of lens parameters and source parameters (lens conditions applied):

              zl: lens redshifts

              zs: source redshifts, lensed condition applied

              sigma: velocity dispersions

              q: axis ratios

              theta_E: Einstein radii

              phi: axis rotation angle

              e1: ellipticity component 1

              e2: ellipticity component 2

              gamma1: shear component 1

              gamma2: shear component 2

              gamma: density profile slope distribution










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.sample_all_routine_sie(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: sampling_routine_sis_nsl(zl, zs, size=1000)

      
      Function to sample SIS lens related parameters.


      :Parameters:

          **zl** : `float`
              lens redshifts

          **zs** : `float`
              source redshifts

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **lens_parameters** : `dict`
              dictionary of sampled lens parameters.
              keys: sigma, theta_E













      ..
          !! processed by numpydoc !!

   .. py:method:: sampling_routine_sie_nsl(zl, zs, size=1000)

      
      Function to sample SIE lens related parameters.


      :Parameters:

          **zl** : `float`
              lens redshifts

          **zs** : `float`
              source redshifts

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **lens_parameters** : `dict`
              dictionary of sampled lens parameters.
              keys: sigma, q, phi













      ..
          !! processed by numpydoc !!

   .. py:method:: sampling_routine_epl_shear_nsl(zl, zs, size=1000)

      
      Function to sample EPL and shear related parameters.


      :Parameters:

          **zl** : `float`
              lens redshifts

          **zs** : `float`
              source redshifts

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **lens_parameters** : `dict`
              dictionary of sampled lens parameters.
              keys: sigma, q, phi, gamma, gamma1, gamma2













      ..
          !! processed by numpydoc !!

   .. py:method:: strongly_lensed_source_redshifts(size=1000)

      
      Function to sample source redshifts, conditioned on the source being strongly lensed.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **redshifts** : `float`
              source redshifts conditioned on the source being strongly lensed










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.strongly_lensed_source_redshifts(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: rjs_with_cross_section_sis(param_dict, cross_section_max=0.0)

      
      Function to conduct rejection sampling wrt einstein radius


      :Parameters:

          **param_dict** : `dict`
              dictionary of lens parameters and source parameters

      :Returns:

          **lens_params** : `dict`
              dictionary of lens parameters after rejection sampling













      ..
          !! processed by numpydoc !!

   .. py:method:: rjs_with_cross_section_sie_feixu(param_dict, cross_section_max=0.0)

      
      Function to conduct rejection sampling wrt cross_section


      :Parameters:

          **param_dict** : `dict`
              dictionary of lens parameters and source parameters

      :Returns:

          **lens_params** : `dict`
              dictionary of lens parameters after rejection sampling













      ..
          !! processed by numpydoc !!

   .. py:method:: rjs_with_cross_section(param_dict, cross_section_max=0.0)

      
      Function to conduct rejection sampling wrt cross_section of EPL+Shear lens


      :Parameters:

          **param_dict** : `dict`
              dictionary of lens parameters and source parameters

      :Returns:

          **lens_params** : `dict`
              dictionary of lens parameters after rejection sampling













      ..
          !! processed by numpydoc !!

   .. py:method:: rjs_with_cross_section_mp(param_dict, cross_section_max=0.0)

      
      Function to conduct rejection sampling wrt cross_section, multiprocessing


      :Parameters:

          **param_dict** : `dict`
              dictionary of lens parameters and source parameters

      :Returns:

          **lens_params** : `dict`
              dictionary of lens parameters after rejection sampling













      ..
          !! processed by numpydoc !!


.. py:function:: cubic_spline_interpolator(xnew, coefficients, x)

   
   Function to interpolate using cubic spline.


   :Parameters:

       **xnew** : `numpy.ndarray`
           new x values.

       **coefficients** : `numpy.ndarray`
           coefficients of the cubic spline.

       **x** : `numpy.ndarray`
           x values.

   :Returns:

       **result** : `numpy.ndarray`
           interpolated values.













   ..
       !! processed by numpydoc !!

.. py:function:: inverse_transform_sampler(size, cdf, x)

   
   Function to sample from the inverse transform method.


   :Parameters:

       **size** : `int`
           number of samples.

       **cdf** : `numpy.ndarray`
           cdf values.

       **x** : `numpy.ndarray`
           x values.

   :Returns:

       **samples** : `numpy.ndarray`
           samples from the cdf.













   ..
       !! processed by numpydoc !!

.. py:function:: cubic_spline_interpolator2d_array(xnew_array, ynew_array, coefficients, x, y)

   
   Function to calculate the interpolated value of snr_partialscaled given the mass ratio (ynew) and total mass (xnew). This is based off 2D bicubic spline interpolation.


   :Parameters:

       **xnew_array** : `numpy.ndarray`
           Total mass of the binary.

       **ynew_array** : `numpy.ndarray`
           Mass ratio of the binary.

       **coefficients** : `numpy.ndarray`
           Array of coefficients for the cubic spline interpolation.

       **x** : `numpy.ndarray`
           Array of total mass values for the coefficients.

       **y** : `numpy.ndarray`
           Array of mass ratio values for the coefficients.

   :Returns:

       **result** : `float`
           Interpolated value of snr_partialscaled.













   ..
       !! processed by numpydoc !!

.. py:function:: save_pickle(file_name, param)

   
   Save a dictionary as a pickle file.


   :Parameters:

       **file_name** : `str`
           pickle file name for storing the parameters.

       **param** : `dict`
           dictionary to be saved as a pickle file.














   ..
       !! processed by numpydoc !!

.. py:function:: interpolator_pickle_path(param_dict_given, directory, sub_directory, interpolator_name)

   
   Function to create the interpolator pickle file path.


   :Parameters:

       **param_dict_given** : `dict`
           dictionary of parameters.

       **directory** : `str`
           directory to store the interpolator.

       **sub_directory** : `str`
           sub-directory to store the interpolator.

       **interpolator_name** : `str`
           name of the interpolator.

   :Returns:

       **path_inv_cdf** : `str`
           path of the interpolator pickle file.

       **it_exist** : `bool`
           if True, the interpolator exists.













   ..
       !! processed by numpydoc !!

.. py:function:: inverse_transform_sampler2d(size, conditioned_y, cdf2d, x2d, y1d)

   
   Function to find sampler interpolator coefficients from the conditioned y.


   :Parameters:

       **size: `int`**
           Size of the sample.

       **conditioned_y: `float`**
           Conditioned y value.

       **cdf2d: `numpy.ndarray`**
           2D array of cdf values.

       **x2d: `numpy.ndarray`**
           2D array of x values.

       **y1d: `numpy.ndarray`**
           1D array of y values.

   :Returns:

       samples: `numpy.ndarray`
           Samples of the conditioned y.













   ..
       !! processed by numpydoc !!

.. py:function:: pdf_cubic_spline_interpolator2d_array(xnew_array, ynew_array, norm_array, coefficients, x, y)

   
   Function to calculate the interpolated value of snr_partialscaled given the mass ratio (ynew) and total mass (xnew). This is based off 2D bicubic spline interpolation.


   :Parameters:

       **xnew_array** : `numpy.ndarray`
           Total mass of the binary.

       **ynew_array** : `numpy.ndarray`
           Mass ratio of the binary.

       **coefficients** : `numpy.ndarray`
           Array of coefficients for the cubic spline interpolation.

       **x** : `numpy.ndarray`
           Array of total mass values for the coefficients.

       **y** : `numpy.ndarray`
           Array of mass ratio values for the coefficients.

   :Returns:

       **result** : `float`
           Interpolated value of snr_partialscaled.













   ..
       !! processed by numpydoc !!

.. py:function:: normal_pdf(x, mean=0.0, std=0.05)

   
   Calculate the value of a normal probability density function.


   :Parameters:

       **x** : `float` or `numpy.ndarray`
           The value(s) at which to evaluate the PDF.

       **mean** : `float`, optional
           The mean of the normal distribution. Default is 0.

       **std** : `float`, optional
           The standard deviation of the normal distribution. Default is 0.05.

   :Returns:

       **pdf** : `float` or `numpy.ndarray`
           The probability density function value(s) at x.













   ..
       !! processed by numpydoc !!

.. py:function:: normal_pdf_2d(x, y, mean_x=0.0, std_x=0.05, mean_y=0.0, std_y=0.05)

   
   Calculate the value of a 2D normal probability density function.


   :Parameters:

       **x** : `float`
           The x-coordinate for which the PDF is evaluated.

       **y** : `float`
           The y-coordinate for which the PDF is evaluated.

       **mean_x** : `float`, optional
           The mean of the normal distribution along the x-axis. Default is 0.

       **std_x** : `float`, optional
           The standard deviation of the normal distribution along the x-axis. Default is 0.05.

       **mean_y** : `float`, optional
           The mean of the normal distribution along the y-axis. Default is 0.

       **std_y** : `float`, optional
           The standard deviation of the normal distribution along the y-axis. Default is 0.05.

   :Returns:

       `float`
           The probability density function value at the given x and y coordinates.













   ..
       !! processed by numpydoc !!

.. py:function:: load_txt_from_module(package, directory, filename)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: phi_cut_SIE(q)

   
   Function to calculate cross-section scaling factor for the SIE lens galaxy from SIS lens galaxy.


   :Parameters:

       **q** : `float: array`
           axis ratio of the lens galaxy

   :Returns:

       **result** : `float: array`
           scaling factor













   ..
       !! processed by numpydoc !!

.. py:function:: phi(s, z, alpha=0.94, beta=1.85, phistar=0.02099, sigmastar=113.78)

   
   Function to calculate the lens galaxy velocity dispersion function at redshift z.
   For Oguri et al. (2018b) model: alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78


   :Parameters:

       **s** : `float: array`
           velocity dispersion of the lens galaxy

       **z** : `float: array`
           redshift of the lens galaxy

       **cosmology_h** : `float`
           Hubble constant

   :Returns:

       **result** : `float: array`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: phi_loc_bernardi(sigma, alpha, beta, phistar, sigmastar)

   
   Function to calculate the local universe velocity dispersion function. Bernardi et al. (2010).
   For Oguri et al. (2018b) model: alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78
   For Choi et al. (2008) model: alpha = 2.32 / 2.67, beta = 2.67, phistar = 8.0e-3*self.cosmo.h**3, sigmastar = 161.0


   :Parameters:

       **sigma** : `float: array`
           velocity dispersion of the lens galaxy

       **alpha, beta, phistar, sigmastar** : `float`
           parameters of the velocity dispersion function

       **cosmology_h** : `float`
           Hubble constant with respect to 100 km/s/Mpc

   :Returns:

       **philoc_** : `float: array`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: phi_q2_ellipticity_hemanta(phi, q)

   
   Function to convert phi and q to ellipticity e1 and e2.


   :Parameters:

       **phi** : `float: array`
           angle of the major axis in radians

       **q** : `float: array`
           axis ratio

   :Returns:

       **e1** : `float: array`
           ellipticity component 1

       **e2** : `float: array`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: lens_redshift_SDSS_catalogue_sis(zs, splineDc, splineDcInv, u, cdf)

   
   Function to sample lens redshift from the SDSS catalogue. Haris et al. (2018) cdf = (10 * u**3 - 15 * u**4 + 6 * u**5)


   :Parameters:

       **zs: `numpy.ndarray` (1D array of float of size=size)**
           Redshift of the source galaxy

       **splineDc: `list`**
           List of spline coefficients for the comoving distance and redshifts

       **splineDcInv: `list`**
           List of spline coefficients for the inverse of comoving distance and redshifts

       **u: `numpy.ndarray` (1D array of float of size=size)**
           corresponding x values wrt to the cdf values
           e.g. u = np.linspace(0, 1, 500)

       **cdf: `numpy.ndarray` (1D array of float of size=size)**
           Cumulative distribution function of the lens redshift distribution between 0 and 1

   :Returns:

       zl: `numpy.ndarray` (1D array of float of size=size)
           Redshift of the lens galaxy corresponding to the zs













   ..
       !! processed by numpydoc !!

.. py:class:: OpticalDepth(npool=4, z_min=0.0, z_max=10.0, cosmology=None, lens_type='epl_shear_galaxy', lens_functions=None, lens_functions_params=None, lens_param_samplers=None, lens_param_samplers_params=None, directory='./interpolator_pickle', create_new_interpolator=False, verbose=False)


   
   Class to calculate the optical depth, velocity dispersion and axis-ratio of a lens galaxy population.


   :Parameters:

       **npool** : int, optional
           Number of processors to use for multiprocessing (default is 4).

       **z_min** : float, optional
           Minimum redshift of the lens galaxy population (default is 0.0).

       **z_max** : float, optional
           Maximum redshift of the lens galaxy population (default is 10.0).

       **cosmology** : astropy.cosmology, optional
           Cosmology object to use (default is FlatLambdaCDM with H0=70, Om0=0.3, Ode0=0.7).

       **lens_type** : str, optional
           Type of the lens galaxy. Must be one of ['sie_galaxy', 'epl_shear_galaxy', 'sis_galaxy'] (default is 'epl_shear_galaxy').

       **lens_functions** : dict, optional
           Dictionary with lens-related functions.

       **lens_functions_params** : dict, optional
           Dictionary with parameters for the lens-related functions.

       **lens_param_samplers** : dict, optional
           Dictionary of sampler functions for velocity dispersion and axis-ratio.

       **lens_param_samplers_params** : dict, optional
           Dictionary with parameters for the priors of the samplers.

       **directory** : str, optional
           Directory where the interpolators are saved (default is './interpolator_pickle').
           If True, creates a new interpolator (default is False).

       **verbose** : bool, optional
           If True, prints additional information during initialization (default is False).





   :Raises:

       ValueError
           If `lens_type` is not in ['sie_galaxy', 'epl_shear_galaxy', 'sis_galaxy'].









   ..
       !! processed by numpydoc !!
   .. py:property:: optical_depth

      
      Function to compute the strong lensing optical depth.


      :Parameters:

          **zs** : `numpy.ndarray` (1D array of floats)
              source redshifts

      :Returns:

          **tau** : `numpy.ndarray` (1D array of floats)
              strong lensing optical depth










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> print(od.optical_depth(np.array([0.1,0.2,0.3])))



      ..
          !! processed by numpydoc !!

   .. py:property:: velocity_dispersion

      
      Class object to sample velocity dispersion. `zl` is required only if velocity dispersion sampler is redshift dependent.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **zl** : `float`
              redshift of the lens galaxy

      :Returns:

          **sigma** : `numpy.ndarray` (1D array of floats)
              velocity dispersion of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> print(od.velocity_dispersion(size=10))



      ..
          !! processed by numpydoc !!

   .. py:property:: axis_ratio

      
      Function to sample axis ratio from rayleigh distribution with given velocity dispersion.


      :Parameters:

          **sigma** : `numpy.ndarray` (1D array of floats)
              velocity dispersion of the lens galaxy

      :Returns:

          **q** : `numpy.ndarray` (1D array of floats)
              axis ratio of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> print(od.axis_ratio(sigma=200.))



      ..
          !! processed by numpydoc !!

   .. py:property:: available_lens_prior_list_and_its_params

      
      Dictionary with list all the available priors and it's corresponding parameters. This is an immutable instance attribute.
















      ..
          !! processed by numpydoc !!

   .. py:property:: available_lens_functions_and_its_params

      
      Dictionary with list all the available lens functions. This is an immutable instance attribute.
















      ..
          !! processed by numpydoc !!

   .. py:method:: default_lens_samplers_and_functions(lens_type)

      
      Function to categorize the lens priors/samplers


      :Parameters:

          **lens_type** : `str`
              lens type
              e.g. 'epl_shear_galaxy' for elliptical power-law galaxy

      :Returns:

          **lens_priors_** : `dict`
              dictionary of priors

          **lens_priors_params_** : `dict`
              dictionary of priors parameters

          **lens_sampler_names_** : `dict`
              dictionary of sampler names

          **lens_functions_** : `dict`
              dictionary of lens functions













      ..
          !! processed by numpydoc !!

   .. py:method:: initialize_decision_dictionary(create_new_interpolator)

      
      Function to initialize decision dictionary for creating interpolator


      :Parameters:

          **create_new_interpolator** : `dict` or `bool`
              dictionary to create new interpolator for velocity dispersion and optical depth.














      ..
          !! processed by numpydoc !!

   .. py:method:: lens_functions_and_sampler_categorization(lens_param_samplers, lens_param_samplers_params, lens_functions, lens_functions_params)

      
      Function to initialize velocity dispersion sampler with it's settings. The reason I am seperating this from lens_priors_categorization is only a specific parameters needs special attention.


      :Parameters:

          **lens_param_samplers** : `str` or `function`
              sampler name or function

          **lens_param_samplers_params** : `dict`
              sampler parameters

          **lens_functions** : `str` or `function`
              lens function name or function

          **lens_functions_params** : `dict`
              lens function parameters














      ..
          !! processed by numpydoc !!

   .. py:method:: axis_ratio_rayleigh(size, sigma, get_attribute=False, **kwargs)

      
      Function to sample axis ratio from rayleigh distribution with given velocity dispersion.


      :Parameters:

          **sigma** : `float: array`
              velocity dispersion of the lens galaxy

          **q_min, q_max** : `float`
              minimum and maximum axis ratio

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample axis ratio

      :Returns:

          **q** : `float: array`
              axis ratio of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(axis_ratio="axis_ratio_rayleigh"))
      >>> print(od.axis_ratio(sigma=200.))



      ..
          !! processed by numpydoc !!

   .. py:method:: axis_ratio_padilla_strauss(size=1000, get_attribute=False, **kwargs)

      
      Function to sample axis ratio using Padilla and Strauss 2008 distribution for axis ratio


      :Parameters:

          **size** : `int`
              sample size

          **q_min, q_max** : `float`
              minimum and maximum axis ratio

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample axis ratio

      :Returns:

          **q** : `float: array`
              axis ratio of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(axis_ratio="axis_ratio_padilla_strauss"))
      >>> print(od.axis_ratio(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: lens_redshift_SDSS_catalogue_numerical(size=1000, zs=None, get_attribute=False, **kwargs)

      
      Function to sample lens redshifts, conditioned on the lens being strongly lensed


      :Parameters:

          **size** : `int`
              sample size

          **zs** : `float`
              source redshifts

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample lens redshifts

      :Returns:

          **zs** : `float: array`
              lens redshifts










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(lens_redshift="lens_redshift_SDSS_catalogue_numerical"))
      >>> print(od.lens_redshift(size=10, zs=1.0))



      ..
          !! processed by numpydoc !!

   .. py:method:: lens_redshift_SDSS_catalogue_hemanta(size=1000, zs=None, get_attribute=False, **kwargs)

      
      Function to sample lens redshifts, conditioned on the lens being strongly lensed


      :Parameters:

          **size** : `int`
              sample size

          **zs** : `float`
              source redshifts

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample lens redshifts

      :Returns:

          **zs** : `float: array`
              lens redshifts










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(lens_redshift="lens_redshift_SDSS_catalogue_numerical"))
      >>> print(od.lens_redshift(size=10, zs=1.0))



      ..
          !! processed by numpydoc !!

   .. py:method:: intrinsic_lens_redshift(size=1000, get_attribute=False, **kwargs)

      
      Function to sample intrinsic lens redshifts, based on the intrinsic velocity dispersion of the lens galaxy.


      :Parameters:

          **size** : `int`
              sample size

          **zs** : `float`
              source redshifts

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample lens redshifts

      :Returns:

          **zs** : `float: array`
              lens redshifts













      ..
          !! processed by numpydoc !!

   .. py:method:: lens_redshift_multiprocessing(zl_scaled2d, zs1d)

      
      Compute the lens redshift distribution using multiprocessing.


      :Parameters:

          **zl_scaled2d** : array_like
              2D array of lens redshifts, scaled by the source redshift.

          **zs1d** : array_like
              1D array of source redshifts.

          **zl_distribution_name** : str
              Name of the lens redshift distribution to compute.

      :Returns:

          **density_array** : array_like
              2D array of the lens redshift distribution.













      ..
          !! processed by numpydoc !!

   .. py:method:: axis_rotation_angle_uniform(size, get_attribute=False, **kwargs)

      
      Function to sample the axis rotation angle of the elliptical lens galaxy from a uniform distribution.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **get_attribute** : `bool`
              if True, returns a function that can be called with size as input

      :Returns:

          **phi** : `numpy.ndarray`
              axis rotation angle of the elliptical lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(axis_rotation_angle="axis_rotation_angle_uniform"))
      >>> print(od.axis_rotation_angle_uniform(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: axis_ratio_uniform(size, get_attribute=False, **kwargs)

      
      Function to sample the axis ratio of the elliptical lens galaxy from a uniform distribution.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **get_attribute** : `bool`
              if True, returns a function that can be called with size as input

      :Returns:

          **q** : `numpy.ndarray`
              axis ratio of the elliptical lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(axis_ratio="axis_ratio_uniform"))
      >>> print(od.axis_ratio_uniform(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: external_shear_normal(size, get_attribute=False, **kwargs)

      
      Function to sample the external shear parameters from a normal distribution.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **get_attribute** : `bool`
              if True, returns a function that can be called with size as input

      :Returns:

          **gamma_1** : `numpy.ndarray`
              shear component in the x-direction

          **gamma_2** : `numpy.ndarray`
              shear component in the y-direction










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(external_shear="external_shear_normal"))
      >>> print(od.external_shear_normal(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: external_shear_numerical_hemanta(size, get_attribute=False, **kwargs)

      
      Function to sample the external shear parameters from a normal distribution.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **get_attribute** : `bool`
              if True, returns a function that can be called with size as input

      :Returns:

          **gamma_1** : `numpy.ndarray`
              shear component in the x-direction

          **gamma_2** : `numpy.ndarray`
              shear component in the y-direction










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(external_shear="external_shear_normal"))
      >>> print(od.external_shear_normal(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: density_profile_slope_normal(size, get_attribute=False, **kwargs)

      
      Function to sample the lens galaxy density profile slope with normal distribution.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample velocity dispersion

          **\*\*kwargs** : `dict`
              additional parameters to be passed to the function,
              e.g. `mean` and `std` for the normal distribution

      :Returns:

          **slope** : `float`
              density profile slope of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(density_profile_slope="density_profile_slope_normal"))
      >>> print(od.density_profile_slope_normal(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: density_profile_slope_numerical_hemanta(size, get_attribute=False, **kwargs)

      
      Function to sample the lens galaxy density profile slope with normal distribution.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample velocity dispersion

          **\*\*kwargs** : `dict`
              additional parameters to be passed to the function,
              e.g. `mean` and `std` for the normal distribution

      :Returns:

          **slope** : `float`
              density profile slope of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(density_profile_slope="density_profile_slope_normal"))
      >>> print(od.density_profile_slope_normal(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: lens_redshift_SDSS_catalogue_sis(size, zs, get_attribute=False, **kwargs)

      
      Function to sample lens redshifts, conditioned on the lens being strongly lensed


      :Parameters:

          **zs** : `float`
              source redshifts

          **get_attribute** : `bool`
              If True, returns a function that can be called with zs as input

      :Returns:

          **zl** : `float`
              lens redshifts










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.lens_redshift_SDSS_catalogue_sis(zs=1.0)



      ..
          !! processed by numpydoc !!

   .. py:method:: velocity_dispersion_gengamma(size, get_attribute=False, **kwargs)

      
      Function to sample velocity dispersion from gengamma distribution


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **a,c** : `float`
              parameters of gengamma distribution
              refer to https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gengamma.html

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample velocity dispersion

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(a=2.32 / 2.67, c=2.67)

      :Returns:

          **sigma** : `numpy.ndarray` (1D array of floats)
              velocity dispersion of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(velocity_dispersion="velocity_dispersion_gengamma"), lens_param_samplers_params=dict(velocity_dispersion=dict(a=2.32 / 2.67, c=2.67)))
      >>> print(od.velocity_dispersion(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: velocity_dispersion_bernardi(size, get_attribute=False, **kwargs)

      
      Function to sample velocity dispersion from Bernardi et al. (2010). This uses inverse transform sampling.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample velocity dispersion

      :Returns:

          **sigma** : `numpy.ndarray` (1D array of floats)
              velocity dispersion of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(velocity_dispersion="velocity_dispersion_bernardi"))
      >>> print(od.velocity_dispersion(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: velocity_dispersion_ewoud(size, zl, get_attribute=False, **kwargs)

      
      Function to sample velocity dispersion (redshift dependent) from Wempe et al. (2022). This uses inverse transform sampling.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **zl** : `float`
              redshift of the lens galaxy

          **get_attribute** : `bool`
              if True, returns a function that can be used to sample velocity dispersion

      :Returns:

          **sigma** : `numpy.ndarray` (1D array of floats)
              velocity dispersion of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(velocity_dispersion="velocity_dispersion_ewoud"))
      >>> print(od.velocity_dispersion(size=10, zl=0.5))



      ..
          !! processed by numpydoc !!

   .. py:method:: cross_section_sis(sigma, zl, zs, **kwargs)

      
      Function to compute the SIS cross-section


      :Parameters:

          **sigma** : `float`
              velocity dispersion of the lens galaxy

          **zl** : `float`
              redshift of the lens galaxy

          **zs** : `float`
              redshift of the source galaxy

      :Returns:

          **cross_section** : `float`
              SIS cross-section










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> print(od.cross_section_sis(sigma=200., zl=0.5, zs=1.0))



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
              Einstein radii of the lens galaxies in radians. Multiply by










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> sigma = 200.0
      >>> zl = 0.5
      >>> zs = 1.0
      >>> lens.compute_einstein_radii(sigma, zl, zs)



      ..
          !! processed by numpydoc !!

   .. py:method:: cross_section_sie_feixu(sigma, zl, zs, q, **kwargs)

      
      Function to compute the SIE cross-section from Fei Xu et al. (2021)


      :Parameters:

          **sigma** : `float`
              velocity dispersion of the lens galaxy

          **zl** : `float`
              redshift of the lens galaxy

          **zs** : `float`
              redshift of the source galaxy

      :Returns:

          **cross_section** : `float`
              SIE cross-section










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> print(od.cross_section_sie_feixu(sigma=200., zl=0.5, zs=1.0, q=1.0))



      ..
          !! processed by numpydoc !!

   .. py:method:: optical_depth_epl_shear_hemanta(zs, get_attribute=False, **kwargs)

      
      Function to compute the strong lensing optical depth (EPL with shear).

      LambdaCDM(H0=70, Om0=0.3, Ode0=0.7) was used to derive the following equation. This is the analytic version of optical depth from z=0 to z=zs.

      :Parameters:

          **zs** : `float`
              source redshifts

      :Returns:

          **tau** : `float`
              strong lensing optical depth










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> print(od.optical_depth_epl_shear_lambdacdm(zs=1.0))



      ..
          !! processed by numpydoc !!

   .. py:method:: optical_depth_sis_haris(zs, get_attribute=False, **kwargs)

      
      Function to compute the strong lensing optical depth (SIS).

      LambdaCDM(H0=70, Om0=0.3, Ode0=0.7) was used to derive the following equation. This is the analytic version of optical depth from z=0 to z=zs.

      :Parameters:

          **zs** : `float`
              source redshifts

      :Returns:

          **tau** : `float`
              strong lensing optical depth










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> print(od.optical_depth_sis_haris(zs=1.0))



      ..
          !! processed by numpydoc !!

   .. py:method:: interpolated_cross_section_function(theta_E, e1, e2, gamma, gamma1, gamma2, get_attribute=False, **kwargs)

      
      Function to compute the cross-section correction factor
















      ..
          !! processed by numpydoc !!

   .. py:method:: create_lookup_table_fuction()

      
      Functions to create lookup tables
      1. Redshift to co-moving distance.
      2. Co-moving distance to redshift.
      3. Redshift to angular diameter distance.
















      ..
          !! processed by numpydoc !!


.. py:function:: phi_cut_SIE(q)

   
   Function to calculate cross-section scaling factor for the SIE lens galaxy from SIS lens galaxy.


   :Parameters:

       **q** : `float: array`
           axis ratio of the lens galaxy

   :Returns:

       **result** : `float: array`
           scaling factor













   ..
       !! processed by numpydoc !!

.. py:function:: phi_q2_ellipticity_hemanta(phi, q)

   
   Function to convert phi and q to ellipticity e1 and e2.


   :Parameters:

       **phi** : `float: array`
           angle of the major axis in radians

       **q** : `float: array`
           axis ratio

   :Returns:

       **e1** : `float: array`
           ellipticity component 1

       **e2** : `float: array`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: inverse_transform_sampler(size, cdf, x)

   
   Function to sample from the inverse transform method.


   :Parameters:

       **size** : `int`
           number of samples.

       **cdf** : `numpy.ndarray`
           cdf values.

       **x** : `numpy.ndarray`
           x values.

   :Returns:

       **samples** : `numpy.ndarray`
           samples from the cdf.













   ..
       !! processed by numpydoc !!

.. py:function:: cubic_spline_interpolator(xnew, coefficients, x)

   
   Function to interpolate using cubic spline.


   :Parameters:

       **xnew** : `numpy.ndarray`
           new x values.

       **coefficients** : `numpy.ndarray`
           coefficients of the cubic spline.

       **x** : `numpy.ndarray`
           x values.

   :Returns:

       **result** : `numpy.ndarray`
           interpolated values.













   ..
       !! processed by numpydoc !!

.. py:function:: cubic_spline_interpolator2d_array(xnew_array, ynew_array, coefficients, x, y)

   
   Function to calculate the interpolated value of snr_partialscaled given the mass ratio (ynew) and total mass (xnew). This is based off 2D bicubic spline interpolation.


   :Parameters:

       **xnew_array** : `numpy.ndarray`
           Total mass of the binary.

       **ynew_array** : `numpy.ndarray`
           Mass ratio of the binary.

       **coefficients** : `numpy.ndarray`
           Array of coefficients for the cubic spline interpolation.

       **x** : `numpy.ndarray`
           Array of total mass values for the coefficients.

       **y** : `numpy.ndarray`
           Array of mass ratio values for the coefficients.

   :Returns:

       **result** : `float`
           Interpolated value of snr_partialscaled.













   ..
       !! processed by numpydoc !!

.. py:function:: inverse_transform_sampler2d(size, conditioned_y, cdf2d, x2d, y1d)

   
   Function to find sampler interpolator coefficients from the conditioned y.


   :Parameters:

       **size: `int`**
           Size of the sample.

       **conditioned_y: `float`**
           Conditioned y value.

       **cdf2d: `numpy.ndarray`**
           2D array of cdf values.

       **x2d: `numpy.ndarray`**
           2D array of x values.

       **y1d: `numpy.ndarray`**
           1D array of y values.

   :Returns:

       samples: `numpy.ndarray`
           Samples of the conditioned y.













   ..
       !! processed by numpydoc !!

.. py:function:: phi_loc_bernardi(sigma, alpha, beta, phistar, sigmastar)

   
   Function to calculate the local universe velocity dispersion function. Bernardi et al. (2010).
   For Oguri et al. (2018b) model: alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78
   For Choi et al. (2008) model: alpha = 2.32 / 2.67, beta = 2.67, phistar = 8.0e-3*self.cosmo.h**3, sigmastar = 161.0


   :Parameters:

       **sigma** : `float: array`
           velocity dispersion of the lens galaxy

       **alpha, beta, phistar, sigmastar** : `float`
           parameters of the velocity dispersion function

       **cosmology_h** : `float`
           Hubble constant with respect to 100 km/s/Mpc

   :Returns:

       **philoc_** : `float: array`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: phi(s, z, alpha=0.94, beta=1.85, phistar=0.02099, sigmastar=113.78)

   
   Function to calculate the lens galaxy velocity dispersion function at redshift z.
   For Oguri et al. (2018b) model: alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78


   :Parameters:

       **s** : `float: array`
           velocity dispersion of the lens galaxy

       **z** : `float: array`
           redshift of the lens galaxy

       **cosmology_h** : `float`
           Hubble constant

   :Returns:

       **result** : `float: array`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: inverse_transform_sampler(size, cdf, x)

   
   Function to sample from the inverse transform method.


   :Parameters:

       **size** : `int`
           number of samples.

       **cdf** : `numpy.ndarray`
           cdf values.

       **x** : `numpy.ndarray`
           x values.

   :Returns:

       **samples** : `numpy.ndarray`
           samples from the cdf.













   ..
       !! processed by numpydoc !!

.. py:function:: cubic_spline_interpolator(xnew, coefficients, x)

   
   Function to interpolate using cubic spline.


   :Parameters:

       **xnew** : `numpy.ndarray`
           new x values.

       **coefficients** : `numpy.ndarray`
           coefficients of the cubic spline.

       **x** : `numpy.ndarray`
           x values.

   :Returns:

       **result** : `numpy.ndarray`
           interpolated values.













   ..
       !! processed by numpydoc !!

.. py:function:: axis_ratio_SIS(sigma)

   
   Function to sample axis ratio from the SIS distribution with given velocity dispersion.


   :Parameters:

       **sigma** : `float: array`
           velocity dispersion of the lens galaxy

   :Returns:

       **q** : `float: array`
           axis ratio of the lens galaxy













   ..
       !! processed by numpydoc !!

.. py:function:: phi(s, z, alpha=0.94, beta=1.85, phistar=0.02099, sigmastar=113.78)

   
   Function to calculate the lens galaxy velocity dispersion function at redshift z.
   For Oguri et al. (2018b) model: alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78


   :Parameters:

       **s** : `float: array`
           velocity dispersion of the lens galaxy

       **z** : `float: array`
           redshift of the lens galaxy

       **cosmology_h** : `float`
           Hubble constant

   :Returns:

       **result** : `float: array`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: phi_loc_bernardi(sigma, alpha, beta, phistar, sigmastar)

   
   Function to calculate the local universe velocity dispersion function. Bernardi et al. (2010).
   For Oguri et al. (2018b) model: alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78
   For Choi et al. (2008) model: alpha = 2.32 / 2.67, beta = 2.67, phistar = 8.0e-3*self.cosmo.h**3, sigmastar = 161.0


   :Parameters:

       **sigma** : `float: array`
           velocity dispersion of the lens galaxy

       **alpha, beta, phistar, sigmastar** : `float`
           parameters of the velocity dispersion function

       **cosmology_h** : `float`
           Hubble constant with respect to 100 km/s/Mpc

   :Returns:

       **philoc_** : `float: array`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: phi_cut_SIE(q)

   
   Function to calculate cross-section scaling factor for the SIE lens galaxy from SIS lens galaxy.


   :Parameters:

       **q** : `float: array`
           axis ratio of the lens galaxy

   :Returns:

       **result** : `float: array`
           scaling factor













   ..
       !! processed by numpydoc !!

.. py:function:: axis_ratio_rayleigh_rvs(sigma, q_min=0.2, q_max=1.0)

   
   Function to sample axis ratio from rayleigh distribution with given velocity dispersion.


   :Parameters:

       **sigma** : `float: array`
           velocity dispersion of the lens galaxy

   :Returns:

       **q** : `float: array`
           axis ratio of the lens galaxy













   ..
       !! processed by numpydoc !!

.. py:function:: velocity_dispersion_z_dependent(size, zl, zl_list, vd_inv_cdf)

   
   Function to sample velocity dispersion from the interpolator


   :Parameters:

       **size: int**
           Number of samples to draw

       **zl: `numpy.ndarray` (1D array of float of size=size)**
           Redshift of the lens galaxy

   :Returns:

       samples: numpy.ndarray
           Samples of velocity dispersion













   ..
       !! processed by numpydoc !!

.. py:function:: lens_redshift_SDSS_catalogue_sis(zs, splineDc, splineDcInv, u, cdf)

   
   Function to sample lens redshift from the SDSS catalogue. Haris et al. (2018) cdf = (10 * u**3 - 15 * u**4 + 6 * u**5)


   :Parameters:

       **zs: `numpy.ndarray` (1D array of float of size=size)**
           Redshift of the source galaxy

       **splineDc: `list`**
           List of spline coefficients for the comoving distance and redshifts

       **splineDcInv: `list`**
           List of spline coefficients for the inverse of comoving distance and redshifts

       **u: `numpy.ndarray` (1D array of float of size=size)**
           corresponding x values wrt to the cdf values
           e.g. u = np.linspace(0, 1, 500)

       **cdf: `numpy.ndarray` (1D array of float of size=size)**
           Cumulative distribution function of the lens redshift distribution between 0 and 1

   :Returns:

       zl: `numpy.ndarray` (1D array of float of size=size)
           Redshift of the lens galaxy corresponding to the zs













   ..
       !! processed by numpydoc !!

.. py:function:: bounded_normal_sample(size, mean, std, low, high)

   
   Function to sample from a normal distribution with bounds.


   :Parameters:

       **mean: `float`**
           Mean of the normal distribution

       **std: `float`**
           Standard deviation of the normal distribution

       **low: `float`**
           Lower bound

       **high: `float`**
           Upper bound














   ..
       !! processed by numpydoc !!

.. py:function:: phi_q2_ellipticity_hemanta(phi, q)

   
   Function to convert phi and q to ellipticity e1 and e2.


   :Parameters:

       **phi** : `float: array`
           angle of the major axis in radians

       **q** : `float: array`
           axis ratio

   :Returns:

       **e1** : `float: array`
           ellipticity component 1

       **e2** : `float: array`
           ..













   ..
       !! processed by numpydoc !!

