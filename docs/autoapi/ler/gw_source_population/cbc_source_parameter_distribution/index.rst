:py:mod:`ler.gw_source_population.cbc_source_parameter_distribution`
====================================================================

.. py:module:: ler.gw_source_population.cbc_source_parameter_distribution

.. autoapi-nested-parse::

   This module contains the classes to generate source compact binary's intrinsic and extrinsic gravitational waves parameters.

   ..
       !! processed by numpydoc !!


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   ler.gw_source_population.cbc_source_parameter_distribution.CBCSourceParameterDistribution




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
           >>> from ler.gw_source_population import CompactBinaryPopulation
           >>> cbc = CompactBinaryPopulation()
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
           If spin_zero=True and spin_precession=True, spin parameters are sampled for precessing binaries.
           if spin_zero=True and spin_precession=False, spin parameters are sampled for aligned/anti-aligned spin binaries.
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
   >>> params = cbc.sample_gw_parameters(size=1000)
   >>> print("sampled parameters=",list(params.keys()))

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
   CompactBinaryPopulation has the following instance methods:

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
   |:meth:`~sample_gw_parameters`        | Function to sample all the       |
   |                                     | intrinsic and extrinsic          |
   |                                     | parameters of compact binaries   |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_source_frame_masses`  | Function to sample source mass1  |
   |                                     | and mass2                        |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_geocent_time`         | Function to sample geocent time  |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_zs`                   | Function to sample source        |
   |                                     | redshift                         |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_ra`                   | Function to sample right         |
   |                                     | ascension (sky position)         |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_dec`                  | Function to sample declination   |
   |                                     | (sky position)                   |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_phase`                | Function to sample coalescence   |
   |                                     | phase                            |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_psi`                  | Function to sample polarization  |
   |                                     | angle                            |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_theta_jn`             | Function to sample inclination   |
   |                                     | angle                            |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_a1`                   | Function to sample spin1         |
   |                                     | magnitude                        |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_a2`                   | Function to sample spin2         |
   |                                     | magnitude                        |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_tilt_1`               | Function to sample tilt1 angle   |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_tilt_2`               | Function to sample tilt2 angle   |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_phi_12`               | Function to sample phi12 angle   |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_phi_jl`               | Function to sample phi_jl angle  |
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
   |:meth:`~binary_masses_BNS_gwcosmo`                                      |
   +-------------------------------------+----------------------------------+
   |                                     | Function to sample source mass1  |
   |                                     | and mass2 from powerlaw          |
   |                                     | distribution.                    |
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
   .. py:property:: available_gw_prior_list_and_its_params

      
      Dictionary with list all the available priors and it's corresponding parameters. This is an immutable instance attribute.













      .. rubric:: Examples

      >>> from ler.gw_source_population import CompactBinaryPopulation
      >>> cbc = CompactBinaryPopulation()
      >>> priors = cbc.available_gw_prior_list_and_its_params
      >>> priors.keys()  # type of priors
      dict_keys(['merger_rate_density', 'source_frame_masses', 'spin', 'geocent_time', 'ra', 'phase', 'psi', 'theta_jn'])
      >>> priors['source_frame_masses'].keys()  # type of source_frame_masses priors
      dict_keys(['binary_masses_BBH_popI_II_powerlaw_gaussian', 'binary_masses_BBH_popIII_lognormal', 'binary_masses_BBH_primordial_lognormal', 'binary_masses_BNS_gwcosmo', 'binary_masses_BNS_bimodal'])
      >>> priors['source_frame_masses']['binary_masses_BBH_popI_II_powerlaw_gaussian'].keys()  # parameters of binary_masses_BBH_popI_II_powerlaw_gaussian
      dict_keys(['mminbh', 'mmaxbh', 'alpha', 'mu_g', 'sigma_g', 'lambda_peak', 'delta_m', 'beta'])



      ..
          !! processed by numpydoc !!

   .. py:property:: sample_zs

      
      Function to sample redshifts with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **zs** : `numpy.ndarray` (1D array of floats)
              Array of redshifts













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_source_frame_masses

      
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

   .. py:property:: sample_geocent_time

      
      Function to sample geocent time with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **geocent_time** : `numpy.ndarray` (1D array of floats)
              Array of geocent_time or time of coalescence













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_ra

      
      Function to sample right ascension of sky position with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **ra** : `numpy.ndarray` (1D array of floats)
              Array of right ascension of sky position













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_dec

      
      Function to sample declination of sky position with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **dec** : `numpy.ndarray` (1D array of floats)
              Array of declination of sky position













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_phase

      
      Function to sample coalescence phase with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **phase** : `numpy.ndarray` (1D array of floats)
              Array of coalescence phase













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_psi

      
      Function to sample polarization angle with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **psi** : `numpy.ndarray` (1D array of floats)
              Array of polarization angle













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_theta_jn

      
      Function to sample theta_jn with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **theta_jn** : `numpy.ndarray` (1D array of floats)
              Array of theta_jn













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_a_1

      
      Function to sample spin magnitude of the compact binaries (body1) with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **a_1** : `numpy.ndarray` (1D array of floats)
              Array of spin magnitude of the compact binaries (body1)













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_a_2

      
      Function to sample spin magnitude of the compact binaries (body2) with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **a_2** : `numpy.ndarray` (1D array of floats)
              Array of spin magnitude of the compact binaries (body2)













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_tilt_1

      
      Function to sample tilt angle of the compact binaries (body1) with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **tilt_1** : `numpy.ndarray` (1D array of floats)
              Array of tilt angle of the compact binaries (body1)













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_tilt_2

      
      Function to sample tilt angle of the compact binaries (body2) with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **tilt_2** : `numpy.ndarray` (1D array of floats)
              Array of tilt angle of the compact binaries (body2)













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_phi_12

      
      Function to sample azimuthal angle between the two spins with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **phi_12** : `numpy.ndarray` (1D array of floats)
              Array of azimuthal angle between the two spins













      ..
          !! processed by numpydoc !!

   .. py:property:: sample_phi_jl

      
      Function to sample azimuthal angle between the total angular momentum and the orbital angular momentum with the initialized prior.


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **phi_jl** : `numpy.ndarray` (1D array of floats)
              Array of azimuthal angle between the total angular momentum and the orbital angular momentum













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

   .. py:method:: lookup_table_luminosity_distance(z_min, z_max, directory)

      
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
      >>> params = cbc.sample_gw_parameters(size=1000)
      >>> print("sampled parameters=",list(params.keys()))



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BBH_popI_II_powerlaw_gaussian(size, mminbh=4.98, mmaxbh=112.5, alpha=3.78, mu_g=32.27, sigma_g=3.88, lambda_peak=0.03, delta_m=4.8, beta=0.81, get_attribute=False, param=None)

      
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

   .. py:method:: binary_masses_BBH_popIII_lognormal(size, m_min=5.0, m_max=150.0, Mc=30.0, sigma=0.3, chunk_size=10000, get_attribute=False, param=None)

      
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

   .. py:method:: binary_masses_BBH_primordial_lognormal(size, m_min=1.0, m_max=100.0, Mc=20.0, sigma=0.3, chunk_size=10000, get_attribute=False, param=None)

      
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

   .. py:method:: binary_masses_BNS_gwcosmo(size, mminns=1.0, mmaxns=3.0, alphans=0.0, get_attribute=False, param=None)

      
      Function to calculate source mass1 and mass2 of BNS from powerlaw distribution (gwcosmo)


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **mminns** : `float`
              Minimum mass of the BNS (Msun)
              default: 1.0

          **mmaxns** : `float`
              Maximum mass of the BNS (Msun)
              default: 3.0

          **alphans** : `float`
              Power law index
              default: 0.0

      :Returns:

          **mass_1_source** : `numpy.ndarray` (1D array of floats)
              Array of mass1 in source frame (Msun)

          **mass_2_source** : `numpy.ndarray` (1D array of floats)
              Array of mass2 in source frame (Msun)










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> m1_src, m2_src = cbc.binary_masses_BNS_gwcosmo(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BNS_bimodal(size, w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1.0, mmax=2.3, resolution=500, create_new=False, get_attribute=False, param=None)

      
      Function to sample source mass1 and mass2 from bimodal distribution. Refer to Will M. Farr et al. 2020 Eqn. 6


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

   .. py:method:: constant_values_n_size(size=100, value=0.0, get_attribute=False, param=None)

      
      Function to sample constant values of size n.


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **value** : `float`
              Constant value
              default: 0.0

          **get_attribute** : `bool`
              If True, return the njitted sampler function with size as the only input where parameters are fixed to the given values.

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(value=0.0)

      :Returns:

          **values** : `numpy.ndarray` (1D array of floats)
              Array of constant values










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> value = cbc.constant_values_n_size(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: sampler_uniform(size, min_=0, max_=np.pi, get_attribute=False, param=None)

      
      Function to sample values from uniform distribution.


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **start_time** : `float`
              Start time of the uniform distribution
              default: 1238166018

          **end_time** : `float`
              End time of the uniform distribution
              default: 1238166018 + 31536000

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

   .. py:method:: sampler_cosine(size, get_attribute=False, param=None)

      
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

   .. py:method:: sampler_sine(size, get_attribute=False, param=None)

      
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

   .. py:method:: source_priors_categorization(event_type, source_priors, event_prior_params)

      
      Function to categorize the event priors and its parameters.


      :Parameters:

          **event_type** : `str`
              Type of event to generate.
              e.g. 'BBH', 'BNS', 'BBH_popIII', 'BBH_primordial', 'NSBH'

          **source_priors** : `dict`
              Dictionary of prior sampler functions for each parameter

          **event_prior_params** : `dict`
              Dictionary of sampler parameters for each GW parameter

      :Returns:

          **source_priors_** : `dict`
              Dictionary of prior sampler functions for each parameter

          **event_prior_params_** : `dict`
              Dictionary of sampler parameters for each parameter

          **sampler_names_** : `dict`
              Dictionary of sampler names with description










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> source_priors, event_prior_params, sampler_names = cbc.source_priors_categorization(event_type='BBH', source_priors=None, event_prior_params=None)
      >>> print(source_priors.keys())
      >>> print(event_prior_params.keys())
      >>> print(sampler_names.keys())



      ..
          !! processed by numpydoc !!


