:py:mod:`ler.lens_galaxy_population`
====================================

.. py:module:: ler.lens_galaxy_population


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   cross_section_interpolator/index.rst
   jit_functions/index.rst
   lens_galaxy_parameter_distribution/index.rst
   lens_parameter_sampler/index.rst
   mp/index.rst
   optical_depth/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.CBCSourceParameterDistribution
   ler.lens_galaxy_population.OpticalDepth
   ler.lens_galaxy_population.ImageProperties
   ler.lens_galaxy_population.LensGalaxyParameterDistribution
   ler.lens_galaxy_population.FunctionConditioning
   ler.lens_galaxy_population.OpticalDepth



Functions
~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.is_njitted
   ler.lens_galaxy_population.cubic_spline_interpolator
   ler.lens_galaxy_population.inverse_transform_sampler
   ler.lens_galaxy_population.cubic_spline_interpolator2d_array
   ler.lens_galaxy_population.save_json
   ler.lens_galaxy_population.load_json
   ler.lens_galaxy_population.interpolator_json_path
   ler.lens_galaxy_population.inverse_transform_sampler2d
   ler.lens_galaxy_population.pdf_cubic_spline_interpolator2d_array
   ler.lens_galaxy_population.normal_pdf
   ler.lens_galaxy_population.normal_pdf_2d
   ler.lens_galaxy_population.comoving_distance
   ler.lens_galaxy_population.angular_diameter_distance
   ler.lens_galaxy_population.angular_diameter_distance_z1z2
   ler.lens_galaxy_population.differential_comoving_volume
   ler.lens_galaxy_population.is_njitted
   ler.lens_galaxy_population.redshift_optimal_spacing
   ler.lens_galaxy_population.phi_cut_SIE
   ler.lens_galaxy_population.phi
   ler.lens_galaxy_population.phi_loc_bernardi
   ler.lens_galaxy_population.phi_q2_ellipticity
   ler.lens_galaxy_population.axis_ratio_rayleigh_pdf
   ler.lens_galaxy_population.phi_cut_SIE
   ler.lens_galaxy_population.phi_q2_ellipticity
   ler.lens_galaxy_population.inverse_transform_sampler
   ler.lens_galaxy_population.cubic_spline_interpolator
   ler.lens_galaxy_population.cubic_spline_interpolator2d_array
   ler.lens_galaxy_population.inverse_transform_sampler2d
   ler.lens_galaxy_population.load_json
   ler.lens_galaxy_population.make_cross_section_reinit
   ler.lens_galaxy_population.lens_redshift_strongly_lensed_njit
   ler.lens_galaxy_population.lens_redshift_strongly_lensed_mp
   ler.lens_galaxy_population.cross_section_unit_mp
   ler.lens_galaxy_population.cross_section_mp
   ler.lens_galaxy_population.cross_section
   ler.lens_galaxy_population.inverse_transform_sampler
   ler.lens_galaxy_population.cubic_spline_interpolator
   ler.lens_galaxy_population.axis_ratio_SIS
   ler.lens_galaxy_population.gamma_
   ler.lens_galaxy_population.cvdf_fit
   ler.lens_galaxy_population.my_derivative
   ler.lens_galaxy_population.pdf_phi_z_div_0
   ler.lens_galaxy_population.phi
   ler.lens_galaxy_population.phi_loc_bernardi
   ler.lens_galaxy_population.phi_cut_SIE
   ler.lens_galaxy_population.axis_ratio_rayleigh_rvs
   ler.lens_galaxy_population.axis_ratio_rayleigh_pdf
   ler.lens_galaxy_population.velocity_dispersion_z_dependent
   ler.lens_galaxy_population.bounded_normal_sample
   ler.lens_galaxy_population.phi_q2_ellipticity
   ler.lens_galaxy_population.sample_sigma_zl
   ler.lens_galaxy_population.phi_q2_ellipticity
   ler.lens_galaxy_population.make_cross_section_reinit
   ler.lens_galaxy_population.create_rejection_sampler
   ler.lens_galaxy_population.create_importance_sampler



Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.CS_UNIT_SLOPE
   ler.lens_galaxy_population.CS_UNIT_INTERCEPT
   ler.lens_galaxy_population.C_LIGHT


.. py:class:: CBCSourceParameterDistribution(z_min=0.0, z_max=10.0, event_type='BBH', source_priors=None, source_priors_params=None, cosmology=None, spin_zero=False, spin_precession=False, directory='./interpolator_json', create_new_interpolator=False)


   Bases: :py:obj:`ler.gw_source_population.cbc_source_redshift_distribution.CBCSourceRedshiftDistribution`

   
   Class for sampling compact binary coalescence source parameters.

   This class generates complete sets of intrinsic and extrinsic gravitational
   wave parameters for compact binary sources including masses, spins, sky
   positions, and orbital parameters. It supports BBH, BNS, NSBH, and primordial
   black hole populations with configurable prior distributions.

   Key Features:

   - Multiple mass distribution models (PowerLaw+Gaussian, lognormal, bimodal)

   - Configurable spin priors (zero, aligned, precessing)

   - Isotropic sky position and orientation sampling

   - Built-in support for population III and primordial black holes

   :Parameters:

       **z_min** : ``float``
           Minimum redshift of the source population.

           default: 0.0

       **z_max** : ``float``
           Maximum redshift of the source population.

           default: 10.0

       **event_type** : ``str``
           Type of compact binary event to generate.

           Options:

           - 'BBH': Binary black hole (Population I/II)

           - 'BNS': Binary neutron star

           - 'NSBH': Neutron star-black hole

           - 'BBH_popIII': Population III binary black hole

           - 'BBH_primordial': Primordial binary black hole

           default: 'BBH'

       **source_priors** : ``dict`` or ``None``
           Dictionary of prior sampler functions for each parameter.

           If None, uses default priors based on event_type.

           default: None

       **source_priors_params** : ``dict`` or ``None``
           Dictionary of parameters for each prior sampler function.

           If None, uses default parameters based on event_type.

           default: None

       **cosmology** : ``astropy.cosmology`` or ``None``
           Cosmology to use for distance calculations.

           default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

       **spin_zero** : ``bool``
           If True, spin parameters are set to zero (no spin sampling).

           default: False

       **spin_precession** : ``bool``
           If True (and spin_zero=False), sample precessing spin parameters.

           If False (and spin_zero=False), sample aligned/anti-aligned spins.

           default: False

       **directory** : ``str``
           Directory to store interpolator JSON files.

           default: './interpolator_json'

       **create_new_interpolator** : ``dict`` or ``bool``
           Configuration for creating new interpolators.

           If bool, applies to all interpolators.

           default: False











   .. rubric:: Examples

   >>> from ler.gw_source_population import CBCSourceParameterDistribution
   >>> cbc = CBCSourceParameterDistribution(event_type='BBH')
   >>> params = cbc.sample_gw_parameters(size=1000)
   >>> print(list(params.keys()))

   Instance Methods
   ----------
   CBCSourceParameterDistribution has the following methods:

   +-----------------------------------------------------+------------------------------------------------+
   | Method                                              | Description                                    |
   +=====================================================+================================================+
   | :meth:`~sample_gw_parameters`                       | Sample all GW parameters for compact binaries  |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~binary_masses_BBH_powerlaw_gaussian`| Sample BBH masses with PowerLaw+PEAK model     |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~binary_masses_BBH_popIII_lognormal`         | Sample pop III BBH masses from lognormal       |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~binary_masses_BBH_primordial_lognormal`     | Sample primordial BBH masses from lognormal    |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~binary_masses_NSBH_broken_powerlaw`         | Sample NSBH masses from broken powerlaw        |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~binary_masses_uniform`                      | Sample masses from uniform distribution        |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~binary_masses_BNS_bimodal`                  | Sample BNS masses from bimodal distribution    |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~constant_values_n_size`                     | Return array of constant values                |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~sampler_uniform`                            | Sample from uniform distribution               |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~sampler_cosine`                             | Sample from cosine distribution                |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~sampler_sine`                               | Sample from sine distribution                  |
   +-----------------------------------------------------+------------------------------------------------+

   Instance Attributes
   ----------
   CBCSourceParameterDistribution has the following attributes:

   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | Attribute                                      | Type                   | Unit  | Description                                    |
   +================================================+========================+=======+================================================+
   | :attr:`~z_min`                                 | ``float``              |       | Minimum redshift of source population          |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~z_max`                                 | ``float``              |       | Maximum redshift of source population          |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~cosmo`                                 | ``astropy.cosmology``  |       | Cosmology for distance calculations            |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~spin_zero`                             | ``bool``               |       | Whether to ignore spin parameters              |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~spin_precession`                       | ``bool``               |       | Whether to use precessing spins                |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~directory`                             | ``str``                |       | Directory for interpolator files               |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~gw_param_samplers`                     | ``dict``               |       | Dictionary of parameter sampler functions      |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~gw_param_samplers_params`              | ``dict``               |       | Dictionary of sampler function parameters      |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~available_gw_prior`                    | ``dict``               |       | Available prior distributions                  |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~source_frame_masses`                   | ``callable``           |       | Sampler for source frame masses                |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~zs`                                    | ``callable``           |       | Sampler for source redshift                    |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~geocent_time`                          | ``callable``           |       | Sampler for geocentric time                    |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~ra`                                    | ``callable``           |       | Sampler for right ascension                    |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~dec`                                   | ``callable``           |       | Sampler for declination                        |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~phase`                                 | ``callable``           |       | Sampler for coalescence phase                  |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~psi`                                   | ``callable``           |       | Sampler for polarization angle                 |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~theta_jn`                              | ``callable``           |       | Sampler for inclination angle                  |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~a_1`                                   | ``callable``           |       | Sampler for spin1 magnitude                    |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~a_2`                                   | ``callable``           |       | Sampler for spin2 magnitude                    |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~tilt_1`                                | ``callable``           |       | Sampler for tilt1 angle                        |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~tilt_2`                                | ``callable``           |       | Sampler for tilt2 angle                        |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~phi_12`                                | ``callable``           |       | Sampler for phi_12 angle                       |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+
   | :attr:`~phi_jl`                                | ``callable``           |       | Sampler for phi_jl angle                       |
   +------------------------------------------------+------------------------+-------+------------------------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:property:: zs

      
      Class object (of FunctionConditioning) for source redshift, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the redshift distribution

      - `pdf`: returns the probability density function of the redshift distribution

      - `function`: returns the redshift distribution function.


      :Returns:

          **zs** : ``numpy.ndarray``
              Array of redshift values.













      ..
          !! processed by numpydoc !!

   .. py:property:: source_frame_masses

      
      Class object (of FunctionConditioning) for source frame masses, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the density profile slope distribution


      :Returns:

          **mass_1_source** : ``numpy.ndarray``
              Array of mass_1_source values in solar masses.

          **mass_2_source** : ``numpy.ndarray``
              Array of mass_2_source values in solar masses.










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc_source_param_dist = CBCSourceParameterDistribution()
      >>> cbc_source_param_dist.source_frame_masses(size=10)



      ..
          !! processed by numpydoc !!

   .. py:property:: geocent_time

      
      Class object (of FunctionConditioning) for geocentric time, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the geocentric time distribution

      - `pdf`: returns the probability density function of the geocentric time distribution

      - `function`: returns the geocentric time distribution function.


      :Returns:

          **geocent_time** : ``numpy.ndarray``
              Array of geocentric time values.













      ..
          !! processed by numpydoc !!

   .. py:property:: ra

      
      Class object (of FunctionConditioning) for right ascension, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the right ascension distribution

      - `pdf`: returns the probability density function of the right ascension distribution

      - `function`: returns the right ascension distribution function.


      :Returns:

          **ra** : ``numpy.ndarray``
              Array of right ascension values.













      ..
          !! processed by numpydoc !!

   .. py:property:: dec

      
      Class object (of FunctionConditioning) for declination, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the declination distribution

      - `pdf`: returns the probability density function of the declination distribution

      - `function`: returns the declination distribution function.


      :Returns:

          **dec** : ``numpy.ndarray``
              Array of declination values.













      ..
          !! processed by numpydoc !!

   .. py:property:: phase

      
      Class object (of FunctionConditioning) for coalescence phase, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the coalescence phase distribution

      - `pdf`: returns the probability density function of the coalescence phase distribution

      - `function`: returns the coalescence phase distribution function.


      :Returns:

          **phase** : ``numpy.ndarray``
              Array of coalescence phase values.













      ..
          !! processed by numpydoc !!

   .. py:property:: psi

      
      Class object (of FunctionConditioning) for polarization angle, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the polarization angle distribution

      - `pdf`: returns the probability density function of the polarization angle distribution

      - `function`: returns the polarization angle distribution function.


      :Returns:

          **geocent_time** : ``numpy.ndarray``
              Array of polarization angle values.













      ..
          !! processed by numpydoc !!

   .. py:property:: theta_jn

      
      Class object (of FunctionConditioning) for inclination angle, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the inclination angle distribution

      - `pdf`: returns the probability density function of the inclination angle distribution

      - `function`: returns the inclination angle distribution function.


      :Returns:

          **theta_jn** : ``numpy.ndarray``
              Array of inclination angle values, i.e. the angle between the line of sight and the orbital angular momentum (rad).













      ..
          !! processed by numpydoc !!

   .. py:property:: a_1

      
      Class object (of FunctionConditioning) for spin1 magnitude, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the spin1 magnitude distribution

      - `pdf`: returns the probability density function of the spin1 magnitude distribution

      - `function`: returns the spin1 magnitude distribution function.


      :Returns:

          **a_1** : ``numpy.ndarray``
              Array of spin magnitude values for the primary body.













      ..
          !! processed by numpydoc !!

   .. py:property:: a_2

      
      Class object (of FunctionConditioning) for spin2 magnitude, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the spin2 magnitude distribution

      - `pdf`: returns the probability density function of the spin2 magnitude distribution

      - `function`: returns the spin2 magnitude distribution function.


      :Returns:

          **a_2** : ``numpy.ndarray``
              Array of spin magnitude values for the secondary body.













      ..
          !! processed by numpydoc !!

   .. py:property:: tilt_1

      
      Class object (of FunctionConditioning) for tilt1 angle, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the tilt1 angle distribution

      - `pdf`: returns the probability density function of the tilt1 angle distribution

      - `function`: returns the tilt1 angle distribution function.


      :Returns:

          **tilt_1** : ``numpy.ndarray``
              Array of the spin tilt angle of the primary body, i.e. the angle between the spin vector and the orbital angular momentum for the primary body (rad).













      ..
          !! processed by numpydoc !!

   .. py:property:: tilt_2

      
      Class object (of FunctionConditioning) for tilt2 angle, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the tilt2 angle distribution

      - `pdf`: returns the probability density function of the tilt2 angle distribution

      - `function`: returns the tilt2 angle distribution function.


      :Returns:

          **tilt_2** : ``numpy.ndarray``
              Array of the spin tilt angle of the secondary body, i.e. the angle between the spin vector and the orbital angular momentum for the secondary body (rad).













      ..
          !! processed by numpydoc !!

   .. py:property:: phi_12

      
      Class object (of FunctionConditioning) for phi_12 angle, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the phi_12 angle distribution

      - `pdf`: returns the probability density function of the phi_12 angle distribution

      - `function`: returns the phi_12 angle distribution function.


      :Returns:

          **phi_12** : ``numpy.ndarray``
              Array of the spin tilt angle between the two spins, i.e., angle between the projections of the two spins onto the orbital plane (rad).













      ..
          !! processed by numpydoc !!

   .. py:property:: phi_jl

      
      Class object (of FunctionConditioning) for phi_jl angle, with rvs/sampler as callback. Can also be a user defined callable sampler.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the phi_jl angle distribution

      - `pdf`: returns the probability density function of the phi_jl angle distribution

      - `function`: returns the phi_jl angle distribution function.


      :Returns:

          **phi_jl** : ``numpy.ndarray``
              Array of the angle values between the orientation of the total angular momentum around the orbital angular momentum (rad).













      ..
          !! processed by numpydoc !!

   .. py:property:: available_gw_prior

      
      Dictionary of all available prior distributions and their parameters.

      This is a dynamically generated dictionary containing available samplers
      for each GW parameter type and their default parameter values.


      :Returns:

          **available_gw_prior** : ``dict``
              Nested dictionary organized by parameter type (e.g., 'source_frame_masses',

              'geocent_time', etc.) with sampler names and default parameters.













      ..
          !! processed by numpydoc !!

   .. py:attribute:: z_min
      :value: 'None'

      
      ``float``

      Minimum redshift of the source population















      ..
          !! processed by numpydoc !!

   .. py:attribute:: z_max
      :value: 'None'

      
      ``float``

      Maximum redshift of the source population















      ..
          !! processed by numpydoc !!

   .. py:attribute:: event_type
      :value: 'None'

      
      ``str``

      Type of event to generate.

      e.g. 'BBH', 'BNS', 'NSBH'















      ..
          !! processed by numpydoc !!

   .. py:attribute:: source_priors
      :value: 'None'

      
      ``dict``

      Dictionary of prior sampler functions.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: source_priors_params
      :value: 'None'

      
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
      :value: 'None'

      
      ``bool``

      If True, spin prior is set to zero.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: spin_precession
      :value: 'False'

      

   .. py:attribute:: directory
      :value: "'./interpolator_json'"

      
      Directory path for storing interpolator JSON files.



      :Returns:

          **directory** : ``str``
              Path to the interpolator storage directory.

              default: './interpolator_json'













      ..
          !! processed by numpydoc !!

   .. py:method:: sample_gw_parameters(size=1000, param=None)

      
      Sample all gravitational wave parameters for compact binaries.

      Generates a complete set of intrinsic and extrinsic parameters including
      masses, redshift, luminosity distance, sky position, orientation, and
      optionally spin parameters.

      :Parameters:

          **size** : ``int``
              Number of samples to draw.

              default: 1000

          **param** : ``dict`` or ``None``
              Dictionary of fixed parameter values.

              Parameters in this dict will not be sampled.

              default: None

      :Returns:

          **gw_parameters** : ``dict``
              Dictionary of sampled GW parameters. The included parameters and their units are as follows (for default settings):

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
              | a_1                |              | spin_1 of the compact binary         |
              +--------------------+--------------+--------------------------------------+
              | a_2                |              | spin_2 of the compact binary         |
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










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> params = cbc.sample_gw_parameters(size=1000)
      >>> print(list(params.keys()))



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BBH_powerlaw_gaussian(size, get_attribute=False, **kwargs)

      
      Sample source masses with PowerLaw+PEAK model for Population I/II BBH.

      Implements the mass distribution model from LIGO-Virgo population analyses
      combining a power-law with a Gaussian peak component.

      :Parameters:

          **size** : ``int``
              Number of samples to draw.

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Model parameters:

              - mminbh: Minimum BH mass (Msun), default: 4.98

              - mmaxbh: Maximum BH mass (Msun), default: 112.5

              - alpha: Power-law spectral index, default: 3.78

              - mu_g: Gaussian peak mean (Msun), default: 32.27

              - sigma_g: Gaussian peak width (Msun), default: 3.88

              - lambda_peak: Fraction in Gaussian component, default: 0.03

              - delta_m: Low-mass tapering range (Msun), default: 4.8

              - beta: Mass ratio power-law index, default: 0.81

      :Returns:

          **mass_1_source** : ``numpy.ndarray``
              Array of primary masses in source frame (Msun).

          **mass_2_source** : ``numpy.ndarray``
              Array of secondary masses in source frame (Msun).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> m1_src, m2_src = cbc.binary_masses_BBH_powerlaw_gaussian(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BBH_popIII_lognormal(size, get_attribute=False, **kwargs)

      
      Sample source masses for Population III BBH from lognormal distribution.

      Based on Eqn. 1 and 4 of Ng et al. 2022 for Population III black holes.

      :Parameters:

          **size** : ``int``
              Number of samples to draw.

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Model parameters:

              - m_min: Minimum BH mass (Msun), default: 5.0

              - m_max: Maximum BH mass (Msun), default: 150.0

              - Mc: Central mass scale (Msun), default: 30.0

              - sigma: Distribution width, default: 0.3

      :Returns:

          **mass_1_source** : ``numpy.ndarray``
              Array of primary masses in source frame (Msun).

          **mass_2_source** : ``numpy.ndarray``
              Array of secondary masses in source frame (Msun).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution(event_type='BBH_popIII')
      >>> m1_src, m2_src = cbc.binary_masses_BBH_popIII_lognormal(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BBH_primordial_lognormal(size, get_attribute=False, **kwargs)

      
      Sample source masses for primordial BBH from lognormal distribution.

      Based on Eqn. 1 and 4 of Ng et al. 2022 for primordial black holes.

      :Parameters:

          **size** : ``int``
              Number of samples to draw.

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Model parameters:

              - m_min: Minimum BH mass (Msun), default: 1.0

              - m_max: Maximum BH mass (Msun), default: 100.0

              - Mc: Central mass scale (Msun), default: 20.0

              - sigma: Distribution width, default: 0.3

      :Returns:

          **mass_1_source** : ``numpy.ndarray``
              Array of primary masses in source frame (Msun).

          **mass_2_source** : ``numpy.ndarray``
              Array of secondary masses in source frame (Msun).













      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_NSBH_broken_powerlaw(size, get_attribute=False, **kwargs)

      
      Sample source masses for NSBH from broken power-law distribution.

      Uses gwcosmo-style broken power-law for black hole mass and power-law
      for neutron star mass.

      :Parameters:

          **size** : ``int``
              Number of samples to draw.

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Model parameters:

              - mminbh: Minimum BH mass (Msun), default: 26

              - mmaxbh: Maximum BH mass (Msun), default: 125

              - alpha_1: Primary power-law index, default: 6.75

              - alpha_2: Secondary power-law index, default: 6.75

              - b: Break point, default: 0.5

              - delta_m: Tapering range (Msun), default: 5

              - mminns: Minimum NS mass (Msun), default: 1.0

              - mmaxns: Maximum NS mass (Msun), default: 3.0

              - alphans: NS mass power-law index, default: 0.0

      :Returns:

          **mass_1_source** : ``numpy.ndarray``
              Array of BH masses in source frame (Msun).

          **mass_2_source** : ``numpy.ndarray``
              Array of NS masses in source frame (Msun).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution(event_type='NSBH')
      >>> m1_src, m2_src = cbc.binary_masses_NSBH_broken_powerlaw(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_uniform(size, get_attribute=False, **kwargs)

      
      Sample source masses from uniform distribution.


      :Parameters:

          **size** : ``int``
              Number of samples to draw.

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Model parameters:

              - m_min: Minimum mass (Msun), default: 1.0

              - m_max: Maximum mass (Msun), default: 3.0

      :Returns:

          **mass_1_source** : ``numpy.ndarray``
              Array of primary masses in source frame (Msun).

          **mass_2_source** : ``numpy.ndarray``
              Array of secondary masses in source frame (Msun).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution()
      >>> m1_src, m2_src = cbc.binary_masses_uniform(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BNS_bimodal(size, get_attribute=False, **kwargs)

      
      Sample BNS masses from bimodal Gaussian distribution.

      Based on Will M. Farr et al. 2020 Eqn. 6 for neutron star mass
      distribution combining two Gaussian peaks.

      :Parameters:

          **size** : ``int``
              Number of samples to draw.

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Model parameters:

              - w: Weight of left peak, default: 0.643

              - muL: Mean of left peak (Msun), default: 1.352

              - sigmaL: Width of left peak (Msun), default: 0.08

              - muR: Mean of right peak (Msun), default: 1.88

              - sigmaR: Width of right peak (Msun), default: 0.3

              - mmin: Minimum mass (Msun), default: 1.0

              - mmax: Maximum mass (Msun), default: 2.3

      :Returns:

          **mass_1_source** : ``numpy.ndarray``
              Array of primary masses in source frame (Msun).

          **mass_2_source** : ``numpy.ndarray``
              Array of secondary masses in source frame (Msun).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceParameterDistribution
      >>> cbc = CBCSourceParameterDistribution(event_type='BNS')
      >>> m1_src, m2_src = cbc.binary_masses_BNS_bimodal(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: constant_values_n_size(size=100, get_attribute=False, **kwargs)

      
      Return array of constant values.


      :Parameters:

          **size** : ``int``
              Number of values to return.

              default: 100

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Model parameters:

              - value: Constant value to return, default: 0.0

      :Returns:

          **values** : ``numpy.ndarray``
              Array of constant values.













      ..
          !! processed by numpydoc !!

   .. py:method:: sampler_uniform(size, get_attribute=False, **kwargs)

      
      Sample values from uniform distribution.


      :Parameters:

          **size** : ``int``
              Number of samples to draw.

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Model parameters:

              - xmin: Minimum value, default: 0.0

              - xmax: Maximum value, default: 1.0

      :Returns:

          **values** : ``numpy.ndarray``
              Array of uniformly distributed values in range [xmin, xmax].













      ..
          !! processed by numpydoc !!

   .. py:method:: sampler_cosine(size, get_attribute=False, **kwargs)

      
      Sample from cosine distribution for declination.

      Samples values in range [-pi/2, pi/2] following a cosine distribution,
      appropriate for isotropic sky position declination.

      :Parameters:

          **size** : ``int``
              Number of samples to draw.

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

      :Returns:

          **values** : ``numpy.ndarray``
              Array of values in range [-pi/2, pi/2] (rad).













      ..
          !! processed by numpydoc !!

   .. py:method:: sampler_sine(size, get_attribute=False, **kwargs)

      
      Sample from sine distribution for inclination angles.

      Samples values in range [0, pi] following a sine distribution,
      appropriate for isotropic orientation angles.

      :Parameters:

          **size** : ``int``
              Number of samples to draw.

          **get_attribute** : ``bool``
              If True, return the sampler object instead of samples.

              default: False

      :Returns:

          **values** : ``numpy.ndarray``
              Array of values in range [0, pi] (rad).













      ..
          !! processed by numpydoc !!


.. py:class:: OpticalDepth(npool=4, z_min=0.0, z_max=10.0, cosmology=None, lens_type='epl_shear_galaxy', lens_functions=None, lens_functions_params=None, lens_param_samplers=None, lens_param_samplers_params=None, directory='./interpolator_json', create_new_interpolator=False, verbose=False)


   
   Class for computing optical depth and lens galaxy population parameters.

   This class calculates strong lensing optical depth, velocity dispersion,
   axis ratio, and other parameters for a lens galaxy population. It supports
   SIS, SIE, and EPL + external shear lens models with customizable samplers
   and interpolators for efficient computation.

   Key Features:

   - Multiple lens model support (SIS, SIE, EPL + shear)

   - Configurable velocity dispersion distributions

   - Cached interpolators for fast optical depth computation

   - Flexible parameter sampling with user-defined priors

   :Parameters:

       **npool** : ``int``
           Number of processors for multiprocessing.

           default: 4

       **z_min** : ``float``
           Minimum redshift of the lens galaxy population.

           default: 0.0

       **z_max** : ``float``
           Maximum redshift of the lens galaxy population.

           default: 10.0

       **cosmology** : ``astropy.cosmology`` or ``None``
           Cosmology object for distance calculations.

           default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

       **lens_type** : ``str``
           Type of lens galaxy model.

           Options:

           - 'epl_shear_galaxy': Elliptical power-law with external shear

           - 'sie_galaxy': Singular isothermal ellipsoid

           - 'sis_galaxy': Singular isothermal sphere

           default: 'epl_shear_galaxy'

       **lens_functions** : ``dict`` or ``None``
           Dictionary with lens-related functions.

           default: None (uses defaults for lens_type)

       **lens_functions_params** : ``dict`` or ``None``
           Dictionary with parameters for lens-related functions.

           default: None

       **lens_param_samplers** : ``dict`` or ``None``
           Dictionary of sampler functions for lens parameters.

           default: None (uses defaults for lens_type)

       **lens_param_samplers_params** : ``dict`` or ``None``
           Dictionary with parameters for the samplers.

           default: None

       **directory** : ``str``
           Directory where interpolators are saved.

           default: './interpolator_json'

       **create_new_interpolator** : ``bool`` or ``dict``
           Whether to create new interpolators.

           default: False

       **verbose** : ``bool``
           If True, prints additional information.

           default: False











   .. rubric:: Examples

   Basic usage:

   >>> from ler.lens_galaxy_population import OpticalDepth
   >>> od = OpticalDepth()
   >>> tau = od.optical_depth(zs=np.array([1.0, 2.0]))

   Instance Methods
   ----------
   OpticalDepth has the following instance methods:

   +-----------------------------------------------------+----------------------------------------------------------+
   | Method                                              | Description                                              |
   +=====================================================+==========================================================+
   | :meth:`~axis_ratio_rayleigh`                        | Sample axis ratio from Rayleigh distribution             |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~axis_ratio_padilla_strauss`                 | Sample axis ratio from Padilla & Strauss 2008            |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~axis_ratio_uniform`                         | Sample axis ratio from uniform distribution              |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~axis_rotation_angle_uniform`                | Sample axis rotation angle from uniform distribution     |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~lens_redshift_strongly_lensed_numerical`    | Sample lens redshift for strong lensing                  |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~lens_redshift_strongly_lensed_sis_haris`                    | Sample SIS lens redshift (Haris et al. 2018)             |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~velocity_dispersion_gengamma`               | Sample velocity dispersion from gengamma distribution    |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~velocity_dispersion_bernardi`               | Sample velocity dispersion (Bernardi et al. 2010)        |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~velocity_dispersion_ewoud`                  | Sample redshift-dependent velocity dispersion            |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~external_shear_normal`                      | Sample external shear from normal distribution           |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~density_profile_slope_normal`               | Sample density profile slope from normal distribution    |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~optical_depth_sis_analytic`                    | Compute SIS optical depth (Haris et al. 2018)            |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~cross_section_sis`                          | Compute SIS cross-section                                |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~cross_section_sie_feixu`                    | Compute SIE cross-section (Fei Xu et al. 2021)           |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~cross_section_epl_shear_numerical`          | Compute EPL+shear cross-section numerically              |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~cross_section_epl_shear_interpolation`      | Compute EPL+shear cross-section via interpolation        |
   +-----------------------------------------------------+----------------------------------------------------------+

   Instance Attributes
   ----------
   OpticalDepth has the following instance attributes:

   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | Attribute                                      | Type                         | Unit  | Description                                              |
   +================================================+==============================+=======+==========================================================+
   | :attr:`~npool`                                 | ``int``                      |       | Number of processors for multiprocessing                 |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~z_min`                                 | ``float``                    |       | Minimum redshift                                         |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~z_max`                                 | ``float``                    |       | Maximum redshift                                         |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~cosmo`                                 | ``astropy.cosmology``        |       | Cosmology object                                         |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~lens_type`                             | ``str``                      |       | Type of lens galaxy model                                |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~directory`                             | ``str``                      |       | Directory for interpolator storage                       |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~optical_depth`                         | ``FunctionConditioning``     |       | Optical depth calculator                                 |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~velocity_dispersion`                   | ``FunctionConditioning``     | km/s  | Velocity dispersion sampler                              |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~axis_ratio`                            | ``FunctionConditioning``     |       | Axis ratio sampler                                       |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~axis_rotation_angle`                   | ``FunctionConditioning``     | rad   | Axis rotation angle sampler                              |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~lens_redshift`                         | ``FunctionConditioning``     |       | Lens redshift sampler                                    |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~external_shear`                        | ``FunctionConditioning``     |       | External shear sampler                                   |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~density_profile_slope`                 | ``FunctionConditioning``     |       | Density profile slope sampler                            |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~cross_section`                         | ``callable``                 | rad²  | Cross-section calculator                                 |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~available_lens_samplers`               | ``dict``                     |       | Available lens parameter samplers                        |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~available_lens_functions`              | ``dict``                     |       | Available lens functions                                 |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:property:: lens_type

      
      Type of lens galaxy model.



      :Returns:

          **lens_type** : ``str``
              Lens type ('epl_shear_galaxy', 'sie_galaxy', or 'sis_galaxy').













      ..
          !! processed by numpydoc !!

   .. py:property:: npool

      
      Number of processors for multiprocessing.



      :Returns:

          **npool** : ``int``
              Number of parallel processors.













      ..
          !! processed by numpydoc !!

   .. py:property:: z_min

      
      Minimum redshift of the lens galaxy population.



      :Returns:

          **z_min** : ``float``
              Minimum redshift.













      ..
          !! processed by numpydoc !!

   .. py:property:: z_max

      
      Maximum redshift of the lens galaxy population.



      :Returns:

          **z_max** : ``float``
              Maximum redshift.













      ..
          !! processed by numpydoc !!

   .. py:property:: cosmo

      
      Cosmology object for distance calculations.



      :Returns:

          **cosmo** : ``astropy.cosmology``
              Cosmology object.













      ..
          !! processed by numpydoc !!

   .. py:property:: directory

      
      Directory for interpolator storage.



      :Returns:

          **directory** : ``str``
              Path to interpolator JSON files.













      ..
          !! processed by numpydoc !!

   .. py:property:: velocity_dispersion

      
      Velocity dispersion sampler object.

      Returns a ``FunctionConditioning`` object with methods:

      - ``rvs(size, zl)``: Sample velocity dispersion values

      - ``pdf(sigma, zl)``: Get probability density

      - ``function(sigma, zl)``: Get number density function


      :Returns:

          **velocity_dispersion** : ``FunctionConditioning``
              Sampler object for velocity dispersion (km/s).










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> sigma = od.velocity_dispersion(size=100, zl=np.ones(100)*0.5)



      ..
          !! processed by numpydoc !!

   .. py:property:: axis_ratio

      
      Axis ratio sampler object.

      Returns a ``FunctionConditioning`` object with methods:

      - ``rvs(size, sigma)``: Sample axis ratio values

      - ``pdf(q, sigma)``: Get probability density

      - ``function(q, sigma)``: Get distribution function


      :Returns:

          **axis_ratio** : ``FunctionConditioning``
              Sampler object for axis ratio.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> q = od.axis_ratio(size=100, sigma=np.ones(100)*200.)



      ..
          !! processed by numpydoc !!

   .. py:property:: axis_rotation_angle

      
      Axis rotation angle sampler object.

      Returns a ``FunctionConditioning`` object with methods:

      - ``rvs(size)``: Sample axis rotation angles

      - ``pdf(phi)``: Get probability density


      :Returns:

          **axis_rotation_angle** : ``FunctionConditioning``
              Sampler object for axis rotation angle (rad).










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> phi = od.axis_rotation_angle(size=100)



      ..
          !! processed by numpydoc !!

   .. py:property:: density_profile_slope

      
      Density profile slope sampler object.

      Returns a ``FunctionConditioning`` object with methods:

      - ``rvs(size)``: Sample density profile slope values

      - ``pdf(gamma)``: Get probability density


      :Returns:

          **density_profile_slope** : ``FunctionConditioning``
              Sampler object for density profile slope.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> gamma = od.density_profile_slope(size=100)



      ..
          !! processed by numpydoc !!

   .. py:property:: external_shear

      
      External shear sampler object.

      Returns a ``FunctionConditioning`` object with methods:

      - ``rvs(size)``: Sample shear components (gamma1, gamma2)

      - ``pdf(gamma1, gamma2)``: Get probability density


      :Returns:

          **external_shear** : ``FunctionConditioning``
              Sampler object for external shear.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> gamma1, gamma2 = od.external_shear(size=100)



      ..
          !! processed by numpydoc !!

   .. py:property:: cross_section

      
      Lensing cross-section calculator.

      Returns a callable that computes lensing cross-section for individual

      lensing events. Input parameters depend on lens type:

      - EPL+shear: zs, zl, sigma, q, phi, gamma, gamma1, gamma2

      - SIE: zs, zl, sigma, q

      - SIS: zs, zl, sigma


      :Returns:

          **cross_section** : ``callable``
              Cross-section function (rad² units).










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> cs = od.cross_section(zs=zs, zl=zl, sigma=sigma, ...)



      ..
          !! processed by numpydoc !!

   .. py:property:: lens_redshift

      
      Lens redshift sampler object.

      Returns a ``FunctionConditioning`` object with methods:

      - ``rvs(size, zs)``: Sample lens redshifts given source redshifts

      - ``pdf(zl, zs)``: Get probability density

      - ``function(zl, zs)``: Get effective lensing cross-section


      :Returns:

          **lens_redshift** : ``FunctionConditioning``
              Sampler object for lens redshift.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> zl = od.lens_redshift(size=100, zs=np.ones(100)*2.0)



      ..
          !! processed by numpydoc !!

   .. py:property:: density_profile_slope_sl

      
      Density profile slope sampler object (strong lensing conditioned).

      Returns a ``FunctionConditioning`` object with methods:

      - ``rvs(size)``: Sample density profile slope values

      - ``pdf(gamma)``: Get probability density


      :Returns:

          **density_profile_slope_sl** : ``FunctionConditioning``
              Sampler object for density profile slope (strong lensing).










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> gamma = od.density_profile_slope_sl(size=100)



      ..
          !! processed by numpydoc !!

   .. py:property:: external_shear_sl

      
      External shear sampler object (strong lensing conditioned).

      Returns a ``FunctionConditioning`` object with methods:

      - ``rvs(size)``: Sample shear components (gamma1, gamma2)

      - ``pdf(gamma1, gamma2)``: Get probability density


      :Returns:

          **external_shear_sl** : ``FunctionConditioning``
              Sampler object for external shear (strong lensing).










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> gamma1, gamma2 = od.external_shear_sl(size=100)



      ..
          !! processed by numpydoc !!

   .. py:property:: optical_depth

      
      Strong lensing optical depth calculator.



      :Returns:

          **optical_depth** : ``FunctionConditioning``
              Function object with `.function(zs)` method that returns \n
              optical depth for given source redshifts.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> tau = od.optical_depth.function(np.array([1.0, 2.0]))



      ..
          !! processed by numpydoc !!

   .. py:property:: available_lens_samplers

      
      Dictionary of available lens parameter samplers and their default parameters.



      :Returns:

          **available_lens_samplers** : ``dict``
              Dictionary with sampler names and default parameters.













      ..
          !! processed by numpydoc !!

   .. py:property:: available_lens_functions

      
      Dictionary of available lens functions and their default parameters.



      :Returns:

          **available_lens_functions** : ``dict``
              Dictionary with function names and default parameters.













      ..
          !! processed by numpydoc !!

   .. py:attribute:: comoving_distance

      

   .. py:attribute:: angular_diameter_distance

      

   .. py:attribute:: angular_diameter_distance_z1z2

      

   .. py:attribute:: differential_comoving_volume

      

   .. py:attribute:: lens_redshift_intrinsic

      

   .. py:method:: axis_ratio_rayleigh(size, sigma, get_attribute=False, **kwargs)

      
      Sample axis ratio from Rayleigh distribution conditioned on velocity dispersion.


      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **sigma** : ``numpy.ndarray``
              Velocity dispersion of the lens galaxy (km/s).

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters (q_min, q_max).

      :Returns:

          **q** : ``numpy.ndarray`` or ``FunctionConditioning``
              Axis ratio samples or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(axis_ratio="axis_ratio_rayleigh"))
      >>> q = od.axis_ratio(size=100, sigma=np.ones(100)*200.)



      ..
          !! processed by numpydoc !!

   .. py:method:: axis_ratio_padilla_strauss(size=1000, get_attribute=False, **kwargs)

      
      Sample axis ratio from Padilla & Strauss (2008) distribution.


      :Parameters:

          **size** : ``int``
              Number of samples to generate.

              default: 1000

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters (q_min, q_max).

      :Returns:

          **q** : ``numpy.ndarray`` or ``FunctionConditioning``
              Axis ratio samples or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(axis_ratio="axis_ratio_padilla_strauss"))
      >>> q = od.axis_ratio(size=100)



      ..
          !! processed by numpydoc !!

   .. py:method:: lens_redshift_strongly_lensed_numerical(size=1000, zs=None, get_attribute=False, **kwargs)

      
      Sample lens redshifts conditioned on strong lensing (numerical method).

      This method computes the lens redshift distribution by numerically
      integrating over the velocity dispersion distribution (galaxy density distribution wrt), cross-section and differential comoving volume.

      :Parameters:

          **size** : ``int``
              Number of samples to generate. \n
              default: 1000

          **zs** : ``numpy.ndarray``
              Source redshifts.

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples. \n
              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters.

      :Returns:

          **zl** : ``numpy.ndarray`` or ``FunctionConditioning``
              Lens redshift samples or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> zl = od.lens_redshift(size=100, zs=np.ones(100)*2.0)



      ..
          !! processed by numpydoc !!

   .. py:method:: axis_rotation_angle_uniform(size, get_attribute=False, **kwargs)

      
      Sample axis rotation angle from uniform distribution.


      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters (phi_min, phi_max).

      :Returns:

          **phi** : ``numpy.ndarray`` or ``FunctionConditioning``
              Axis rotation angle samples (rad) or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> phi = od.axis_rotation_angle(size=100)



      ..
          !! processed by numpydoc !!

   .. py:method:: axis_ratio_uniform(size, get_attribute=False, **kwargs)

      
      Sample axis ratio from uniform distribution.


      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters (q_min, q_max).

      :Returns:

          **q** : ``numpy.ndarray`` or ``FunctionConditioning``
              Axis ratio samples or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(axis_ratio="axis_ratio_uniform"))
      >>> q = od.axis_ratio(size=100)



      ..
          !! processed by numpydoc !!

   .. py:method:: external_shear_normal(size, get_attribute=False, **kwargs)

      
      Sample external shear parameters from 2D normal distribution.


      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters (mean, std).

      :Returns:

          **shear** : ``numpy.ndarray`` or ``FunctionConditioning``
              Array of shape (2, size) with gamma1, gamma2 or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> gamma1, gamma2 = od.external_shear(size=100)



      ..
          !! processed by numpydoc !!

   .. py:method:: density_profile_slope_normal(size, get_attribute=False, **kwargs)

      
      Sample density profile slope from normal distribution.


      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters (mean, std).

      :Returns:

          **gamma** : ``numpy.ndarray`` or ``FunctionConditioning``
              Density profile slope samples or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> gamma = od.density_profile_slope(size=100)



      ..
          !! processed by numpydoc !!

   .. py:method:: lens_redshift_strongly_lensed_sis_haris(size, zs, get_attribute=False, **kwargs)

      
      Sample SIS lens redshifts using Haris et al. (2018) distribution.


      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **zs** : ``numpy.ndarray``
              Source redshifts.

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters.

      :Returns:

          **zl** : ``numpy.ndarray`` or ``FunctionConditioning``
              Lens redshift samples or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_type='sis_galaxy')
      >>> zl = od.lens_redshift(size=100, zs=np.ones(100)*2.0)



      ..
          !! processed by numpydoc !!

   .. py:method:: velocity_dispersion_gengamma(size, get_attribute=False, **kwargs)

      
      Sample velocity dispersion from generalized gamma distribution.


      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters (a, c for gengamma distribution).

              See scipy.stats.gengamma for details.

      :Returns:

          **sigma** : ``numpy.ndarray`` or ``FunctionConditioning``
              Velocity dispersion samples (km/s) or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(
      ...     velocity_dispersion="velocity_dispersion_gengamma"))
      >>> sigma = od.velocity_dispersion(size=100)



      ..
          !! processed by numpydoc !!

   .. py:method:: velocity_dispersion_bernardi(size, get_attribute=False, **kwargs)

      
      Sample velocity dispersion from Bernardi et al. (2010) distribution.

      Uses inverse transform sampling on the velocity dispersion function.

      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters.

      :Returns:

          **sigma** : ``numpy.ndarray`` or ``FunctionConditioning``
              Velocity dispersion samples (km/s) or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(
      ...     velocity_dispersion="velocity_dispersion_bernardi"))
      >>> sigma = od.velocity_dispersion(size=100)



      ..
          !! processed by numpydoc !!

   .. py:method:: velocity_dispersion_ewoud(size, zl, get_attribute=False, **kwargs)

      
      Sample redshift-dependent velocity dispersion from Wempe et al. (2022).

      Uses inverse transform sampling with redshift-dependent velocity
      dispersion function.

      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **zl** : ``numpy.ndarray``
              Lens redshifts.

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters.

      :Returns:

          **sigma** : ``numpy.ndarray`` or ``FunctionConditioning``
              Velocity dispersion samples (km/s) or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> sigma = od.velocity_dispersion(size=100, zl=np.ones(100)*0.5)



      ..
          !! processed by numpydoc !!

   .. py:method:: optical_depth_numerical(zs, get_attribute=False, **kwargs)

      
      Helper to compute optical depth numerically by integrating lens redshift.


      :Parameters:

          **zs** : ``numpy.ndarray``
              Source redshifts.

          **get_attribute** : ``bool``
              If True, returns the function object instead of values. \n
              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters.

      :Returns:

          **tau** : ``numpy.ndarray`` or ``FunctionConditioning``
              Optical depth values or function object.













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

   .. py:method:: optical_depth_sis_analytic(zs, get_attribute=False, **kwargs)

      
      Function to compute the strong lensing optical depth (SIS).
      LambdaCDM(H0=70, Om0=0.3, Ode0=0.7) was used to derive the following equation. This is the analytic version of optical depth from z=0 to z=zs.


      :Parameters:

          **zs** : ``numpy.ndarray``
              Source redshifts.

          **get_attribute** : ``bool``
              If True, returns the function object instead of values. \n
              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters.

      :Returns:

          **tau** : ``numpy.ndarray`` or ``FunctionConditioning``
              Optical depth values or function object.













      ..
          !! processed by numpydoc !!

   .. py:method:: cross_section_sis(zs=None, zl=None, sigma=None, get_attribute=False, **kwargs)

      
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
      >>> print(self.cross_section_sis(sigma=200., zl=0.5, zs=1.0))



      ..
          !! processed by numpydoc !!

   .. py:method:: cross_section_sie_feixu(zs=None, zl=None, sigma=None, q=None, get_attribute=False, **kwargs)

      
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
      >>> print(self.cross_section_sie_feixu(sigma=200., zl=0.5, zs=1.0, q=1.0))



      ..
          !! processed by numpydoc !!

   .. py:method:: cross_section_epl_shear_numerical(theta_E, gamma, gamma1, gamma2, q=None, phi=None, e1=None, e2=None, verbose=False, **kwargs)

      
      Function to compute the strong lensing cross-section numerically for EPL + external shear lenses.


      :Parameters:

          **theta_E** : `numpy.ndarray`
              Einstein radii of the lens galaxies in radians

          **gamma** : `numpy.ndarray`
              external shear magnitudes of the lens galaxies

          **gamma1** : `numpy.ndarray`
              external shear component 1 of the lens galaxies

          **gamma2** : `numpy.ndarray`
              external shear component 2 of the lens galaxies

          **q** : `numpy.ndarray`
              axis ratios of the lens galaxies

          **phi** : `numpy.ndarray`
              axis rotation angles of the lens galaxies in radians

          **e1** : `numpy.ndarray`
              ellipticity component 1 of the lens galaxies

          **e2** : `numpy.ndarray`
              ellipticity component 2 of the lens galaxies

      :Returns:

          **cross_section** : `numpy.ndarray`
              strong lensing cross-section of the lens galaxies in square radians













      ..
          !! processed by numpydoc !!

   .. py:method:: create_parameter_grid(size_list=[25, 25, 45, 15, 15])


   .. py:method:: cross_section_epl_shear_interpolation_init(file_path, size_list)


   .. py:method:: cross_section_epl_shear_interpolation(zs, zl, sigma, q, phi, gamma, gamma1, gamma2, get_attribute=False, size_list=[25, 25, 45, 15, 15], **kwargs)

      
      Function to compute the cross-section correction factor
















      ..
          !! processed by numpydoc !!


.. py:class:: ImageProperties(npool=4, n_min_images=2, n_max_images=4, time_window=365 * 24 * 3600 * 20, lens_model_list=['EPL_NUMBA', 'SHEAR'], cosmology=None, spin_zero=True, spin_precession=False)


   
   Class to compute image properties of strongly lensed gravitational wave events.

   This class solves the lens equation to find image positions, magnifications,
   time delays, and image types (morse phase) for strongly lensed sources. It uses
   multiprocessing for efficient computation of large samples.

   Key Features:

   - Solves lens equations using multiprocessing for efficiency

   - Computes image positions, magnifications, and time delays

   - Classifies image types using morse phase

   - Calculates detection probabilities for lensed images

   :Parameters:

       **npool** : ``int``
           Number of processes for multiprocessing.

           default: 4

       **n_min_images** : ``int``
           Minimum number of images required for a valid lensing event.

           default: 2

       **n_max_images** : ``int``
           Maximum number of images to consider per event.

           default: 4

       **time_window** : ``float``
           Time window for lensed events (units: seconds).

           default: 365*24*3600*20 (20 years)

       **lens_model_list** : ``list``
           List of lens models to use.

           default: ['EPL_NUMBA', 'SHEAR']

       **cosmology** : ``astropy.cosmology`` or ``None``
           Cosmology for distance calculations.

           If None, uses default LambdaCDM.

           default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

       **spin_zero** : ``bool``
           If True, spin parameters are set to zero (no spin sampling).

           default: False

       **spin_precession** : ``bool``
           If True (and spin_zero=False), sample precessing spin parameters.

           If False (and spin_zero=False), sample aligned/anti-aligned spins.

           default: False











   .. rubric:: Examples

   Basic usage:

   >>> from ler.image_properties import ImageProperties
   >>> ip = ImageProperties()
   >>> lens_parameters = dict(
   ...     zs=np.array([2.0]),
   ...     zl=np.array([0.5]),
   ...     gamma1=np.array([0.0]),
   ...     gamma2=np.array([0.0]),
   ...     phi=np.array([0.0]),
   ...     q=np.array([0.8]),
   ...     gamma=np.array([2.0]),
   ...     theta_E=np.array([1.0])
   ... )
   >>> result = ip.image_properties(lens_parameters)
   >>> print(result.keys())

   Instance Methods
   ----------
   ImageProperties has the following methods:

   +-----------------------------------------------------+------------------------------------------------+
   | Method                                              | Description                                    |
   +=====================================================+================================================+
   | :meth:`~image_properties`                           | Compute image properties for lensed events     |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~get_lensed_snrs`                            | Compute detection probability for lensed images|
   +-----------------------------------------------------+------------------------------------------------+

   Instance Attributes
   ----------
   ImageProperties has the following attributes:

   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
   | Attribute                                           | Type                      | Unit     | Description                                    |
   +=====================================================+===========================+==========+================================================+
   | :attr:`~npool`                                      | ``int``                   |          | Number of multiprocessing workers              |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
   | :attr:`~n_min_images`                               | ``int``                   |          | Minimum number of images required              |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
   | :attr:`~n_max_images`                               | ``int``                   |          | Maximum number of images per event             |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
   | :attr:`~time_window`                                | ``float``                 | s        | Time window for lensed events                  |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
   | :attr:`~lens_model_list`                            | ``list``                  |          | List of lens models                            |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
   | :attr:`~cosmo`                                      | ``astropy.cosmology``     |          | Cosmology for calculations                     |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
   | :attr:`~spin_zero`                                  | ``bool``                  |          | Flag for zero spin assumption                  |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
   | :attr:`~spin_precession`                            | ``bool``                  |          | Flag for spin precession                       |
   +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:property:: npool

      
      Number of multiprocessing workers.



      :Returns:

          **npool** : ``int``
              Number of processes for multiprocessing.

              default: 4













      ..
          !! processed by numpydoc !!

   .. py:property:: n_min_images

      
      Minimum number of images required for a valid lensing event.



      :Returns:

          **n_min_images** : ``int``
              Minimum number of images required.

              default: 2













      ..
          !! processed by numpydoc !!

   .. py:property:: n_max_images

      
      Maximum number of images per event.



      :Returns:

          **n_max_images** : ``int``
              Maximum number of images to consider per event.

              default: 4













      ..
          !! processed by numpydoc !!

   .. py:property:: lens_model_list

      
      List of lens models to use.



      :Returns:

          **lens_model_list** : ``list``
              List of lens model names.

              default: ['EPL_NUMBA', 'SHEAR']













      ..
          !! processed by numpydoc !!

   .. py:property:: spin_zero

      
      Flag for zero spin assumption.



      :Returns:

          **spin_zero** : ``bool``
              Whether to assume zero spin for compact objects.

              If True, spin parameters are set to zero (no spin sampling).

              If False, spin parameters are sampled.

              default: False













      ..
          !! processed by numpydoc !!

   .. py:property:: spin_precession

      
      Flag for spin precession.



      :Returns:

          **spin_precession** : ``bool``
              Whether to include spin precession effects.

              If True (and spin_zero=False), sample precessing spin parameters.

              If False (and spin_zero=False), sample aligned/anti-aligned spins.

              default: False













      ..
          !! processed by numpydoc !!

   .. py:property:: time_window

      
      Time window for lensed events.



      :Returns:

          **time_window** : ``float``
              Time window for lensed events (units: s).

              default: 365*24*3600*20 (20 years)













      ..
          !! processed by numpydoc !!

   .. py:property:: cosmo

      
      Astropy cosmology object for calculations.



      :Returns:

          **cosmo** : ``astropy.cosmology``
              Cosmology used for distance calculations.

              default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)













      ..
          !! processed by numpydoc !!

   .. py:method:: image_properties(lens_parameters)

      
      Compute image properties for strongly lensed events.

      Solves the lens equation using multiprocessing to find image positions,
      magnifications, time delays, and image types for each lensing event.

      :Parameters:

          **lens_parameters** : ``dict``
              Dictionary containing lens and source parameters with keys:

              - 'zs': source redshift (array)

              - 'zl': lens redshift (array)

              - 'gamma1': external shear component 1 (array)

              - 'gamma2': external shear component 2 (array)

              - 'phi': position angle of lens ellipticity (array)

              - 'q': axis ratio of lens (array)

              - 'gamma': power-law slope of mass density (array)

              - 'theta_E': Einstein radius in radians (array)

      :Returns:

          **lens_parameters** : ``dict``
              Updated dictionary with additional image properties:

              - 'x0_image_positions': x-coordinates of images (shape: size x n_max_images)

              - 'x1_image_positions': y-coordinates of images (shape: size x n_max_images)

              - 'magnifications': magnification factors (shape: size x n_max_images)

              - 'time_delays': time delays relative to first image (shape: size x n_max_images, units: s)

              - 'image_type': morse phase classification (1=minimum, 2=saddle, 3=maximum)

              - 'n_images': number of images per event (array)

              - 'x_source': source x-position (array)

              - 'y_source': source y-position (array)













      ..
          !! processed by numpydoc !!

   .. py:method:: get_lensed_snrs(lensed_param, pdet_calculator, list_of_detectors=None)

      
      Compute detection probability for each lensed image.

      Calculates the effective luminosity distance, geocent time, and phase
      for each image accounting for magnification and morse phase, then
      computes detection probabilities using the provided calculator.

      :Parameters:

          **lensed_param** : ``dict``
              Dictionary containing lensed source and image parameters with keys:

              - 'mass_1', 'mass_2': detector-frame masses (array)

              - 'luminosity_distance' or 'effective_luminosity_distance': distance (array)

              - 'geocent_time' or 'effective_geocent_time': GPS time (array)

              - 'phase' or 'effective_phase': coalescence phase (array)

              - 'theta_jn', 'psi', 'ra', 'dec': orientation and position (arrays)

              - 'magnifications': image magnifications (shape: size x n_max_images)

              - 'time_delays': image time delays (shape: size x n_max_images)

              - 'image_type': morse phase type (shape: size x n_max_images)

          **pdet_calculator** : ``callable``
              Function that computes detection probability given GW parameters.

          **list_of_detectors** : ``list`` or ``None``
              List of detector names (e.g., ['H1', 'L1', 'V1']) for per-detector results.

              default: None

      :Returns:

          **result_dict** : ``dict``
              Dictionary containing:

              - 'pdet_net': network detection probability (shape: size x n_max_images)

              - Individual detector probabilities if list_of_detectors provided

          **lensed_param** : ``dict``
              Updated dictionary with effective parameters:

              - 'effective_luminosity_distance': magnification-corrected distance

              - 'effective_geocent_time': time-delay-corrected GPS time

              - 'effective_phase': morse-phase-corrected coalescence phase













      ..
          !! processed by numpydoc !!


.. py:function:: is_njitted(func)


.. py:class:: LensGalaxyParameterDistribution(npool=4, z_min=0.0, z_max=10.0, cosmology=None, event_type='BBH', lens_type='epl_shear_galaxy', lens_functions=None, lens_functions_params=None, lens_param_samplers=None, lens_param_samplers_params=None, directory='./interpolator_json', create_new_interpolator=False, buffer_size=1000, **kwargs)


   Bases: :py:obj:`ler.gw_source_population.CBCSourceParameterDistribution`, :py:obj:`ler.image_properties.ImageProperties`, :py:obj:`ler.lens_galaxy_population.optical_depth.OpticalDepth`

   
   Sample lens galaxy parameters conditioned on strong lensing.

   This class handles the distribution of lens galaxy parameters such as velocity
   dispersion, axis ratio, axis rotation angle, shear, and density profile slope.
   It samples source parameters conditioned on the source being strongly lensed
   using cross-section based rejection or importance sampling.

   Key Features:

   - Samples lens parameters using EPL+shear galaxy model

   - Supports rejection and importance sampling based on cross-section

   - Computes optical depth weighted source redshift distributions

   - Integrates with GW source population and image property calculations

   :Parameters:

       **npool** : ``int``
           Number of processors to use for parallel sampling.

           default: 4

       **z_min** : ``float``
           Minimum redshift for source and lens populations.

           default: 0.0

       **z_max** : ``float``
           Maximum redshift for source and lens populations.

           default: 10.0

       **cosmology** : ``astropy.cosmology`` or ``None``
           Cosmology object for distance calculations.

           default: None (uses FlatLambdaCDM with H0=70, Om0=0.3)

       **event_type** : ``str``
           Type of compact binary coalescence event.

           Options:

           - 'BBH': Binary black hole

           - 'BNS': Binary neutron star

           - 'NSBH': Neutron star-black hole

           default: 'BBH'

       **lens_type** : ``str``
           Type of lens galaxy model to use.

           default: 'epl_shear_galaxy'

       **lens_functions** : ``dict`` or ``None``
           Dictionary specifying lens-related functions.

           default: None (uses defaults from OpticalDepth)

       **lens_functions_params** : ``dict`` or ``None``
           Parameters for lens functions.

           default: None

       **lens_param_samplers** : ``dict`` or ``None``
           Dictionary specifying lens parameter sampling functions.

           default: None (uses defaults from OpticalDepth)

       **lens_param_samplers_params** : ``dict`` or ``None``
           Parameters for lens parameter samplers.

           default: None

       **directory** : ``str``
           Directory for storing interpolator files.

           default: './interpolator_json'

       **create_new_interpolator** : ``bool``
           If True, recreates interpolators even if files exist.

           default: False

       **buffer_size** : ``int``
           Buffer size for batch sampling of lens parameters.

           default: 1000

       **\*\*kwargs** : ``dict``
           Additional keyword arguments passed to parent classes:

           :class:`~ler.gw_source_population.CBCSourceParameterDistribution`,

           :class:`~ler.image_properties.ImageProperties`,

           :class:`~ler.lens_galaxy_population.OpticalDepth`.











   .. rubric:: Examples

   Basic usage:

   >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
   >>> lens = LensGalaxyParameterDistribution()
   >>> lensed_params = lens.sample_lens_parameters(size=1000)
   >>> print(lensed_params.keys())

   Instance Methods
   ----------
   LensGalaxyParameterDistribution has the following methods:

   +-----------------------------------------------------+------------------------------------------------+
   | Method                                              | Description                                    |
   +=====================================================+================================================+
   | :meth:`~sample_lens_parameters`                     | Sample lens and source parameters              |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~sample_all_routine_epl_shear_intrinsic`    | Sample EPL+shear lens parameters from intrinsic |
   |                                                     | distributions                                  |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~sample_all_routine_epl_shear_sl`            | Sample EPL+shear lens parameters with strong   |
   |                                                     | lensing condition                              |
   +-----------------------------------------------------+------------------------------------------------+
   | :meth:`~strongly_lensed_source_redshifts`           | Sample source redshifts with lensing condition |
   +-----------------------------------------------------+------------------------------------------------+

   Instance Attributes
   ----------
   LensGalaxyParameterDistribution has the following attributes:

   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | Attribute                                      | Type                 | Unit  | Description                                    |
   +================================================+======================+=======+================================================+
   | :attr:`~npool`                                 | ``int``              |       | Number of processors for parallel computation  |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | :attr:`~z_min`                                 | ``float``            |       | Minimum redshift                               |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | :attr:`~z_max`                                 | ``float``            |       | Maximum redshift                               |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | :attr:`~cosmo`                                 | ``astropy.cosmology``|       | Cosmology object for calculations              |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | :attr:`~event_type`                            | ``str``              |       | Type of CBC event (BBH, BNS, NSBH)             |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | :attr:`~directory`                             | ``str``              |       | Path to interpolator storage directory         |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | :attr:`~lens_param_samplers`                   | ``dict``             |       | Dictionary of lens parameter sampler names     |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | :attr:`~lens_param_samplers_params`            | ``dict``             |       | Parameters for lens parameter samplers         |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | :attr:`~lens_functions`                        | ``dict``             |       | Dictionary of lens function names              |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+
   | :attr:`~normalization_pdf_z_lensed`            | ``float``            |       | Normalization constant for lensed source z pdf |
   +------------------------------------------------+----------------------+-------+------------------------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:property:: normalization_pdf_z_lensed

      
      Normalization constant for the lensed source redshift pdf.

      This constant is used to normalize the probability distribution

      of source redshifts conditioned on strong lensing. It is computed

      by integrating the merger rate density times optical depth.


      :Returns:

          **normalization_pdf_z_lensed** : ``float``
              Normalization constant for lensed redshift distribution.













      ..
          !! processed by numpydoc !!

   .. py:property:: lens_param_samplers

      
      Dictionary of lens parameter sampler function names.



      :Returns:

          **lens_param_samplers** : ``dict``
              Dictionary mapping parameter names to sampler function names.

              Keys include: 'source_redshift_sl', 'lens_redshift',

              'velocity_dispersion', 'axis_ratio', 'axis_rotation_angle',

              'external_shear', 'density_profile_slope'.













      ..
          !! processed by numpydoc !!

   .. py:property:: lens_param_samplers_params

      
      Dictionary of parameters for lens parameter samplers.



      :Returns:

          **lens_param_samplers_params** : ``dict``
              Dictionary with sampler parameters.

              Each key corresponds to a sampler in lens_param_samplers.













      ..
          !! processed by numpydoc !!

   .. py:property:: lens_functions

      
      Dictionary of lens-related function names.



      :Returns:

          **lens_functions** : ``dict``
              Dictionary mapping function types to function names.

              Keys include: 'param_sampler_type', 'cross_section_based_sampler',

              'optical_depth', 'cross_section'.













      ..
          !! processed by numpydoc !!

   .. py:attribute:: event_type
      :value: "'BBH'"

      
      ``str``

      Type of event to generate.

      e.g. 'BBH', 'BNS', 'NSBH'















      ..
          !! processed by numpydoc !!

   .. py:attribute:: sample_source_redshift_sl

      

   .. py:attribute:: sample_lens_parameters_routine

      

   .. py:attribute:: cross_section_based_sampler

      

   .. py:method:: sample_lens_parameters(size=1000)

      
      Sample lens galaxy and source parameters conditioned on strong lensing.

      This method samples both lens galaxy parameters (velocity dispersion, axis
      ratio, shear, etc.) and gravitational wave source parameters, with the
      source redshift distribution weighted by strong lensing optical depth.

      :Parameters:

          **size** : ``int``
              Number of lens-source parameter sets to sample.

              default: 1000

      :Returns:

          **lens_parameters** : ``dict``
              Dictionary containing sampled lens and source parameters.

              The included parameters and their units are as follows (for default settings):

              +------------------------------+-----------+-------------------------------------------------------+
              | Parameter                    | Units     | Description                                           |
              +==============================+===========+=======================================================+
              | zl                           |           | redshift of the lens                                  |
              +------------------------------+-----------+-------------------------------------------------------+
              | zs                           |           | redshift of the source                                |
              +------------------------------+-----------+-------------------------------------------------------+
              | sigma                        | km s^-1   | velocity dispersion                                   |
              +------------------------------+-----------+-------------------------------------------------------+
              | q                            |           | axis ratio                                            |
              +------------------------------+-----------+-------------------------------------------------------+
              | theta_E                      | arcsec    | Einstein radius                                       |
              +------------------------------+-----------+-------------------------------------------------------+
              | phi                          | rad       | axis rotation angle. counter-clockwise from the       |
              |                              |           | positive x-axis (RA-like axis) to the major axis of   |
              |                              |           | the projected mass distribution.                      |
              +------------------------------+-----------+-------------------------------------------------------+
              | gamma                        |           | density profile slope of EPL galaxy                   |
              +------------------------------+-----------+-------------------------------------------------------+
              | gamma1                       |           | external shear component in the x-direction           |
              |                              |           | (RA-like axis)                                        |
              +------------------------------+-----------+-------------------------------------------------------+
              | gamma2                       |           | external shear component in the y-direction           |
              |                              |           | (Dec-like axis)                                       |
              +------------------------------+-----------+-------------------------------------------------------+
              | geocent_time                 | s         | geocent time                                          |
              +------------------------------+-----------+-------------------------------------------------------+
              | ra                           | rad       | right ascension                                       |
              +------------------------------+-----------+-------------------------------------------------------+
              | dec                          | rad       | declination                                           |
              +------------------------------+-----------+-------------------------------------------------------+
              | phase                        | rad       | phase of GW at reference freq                         |
              +------------------------------+-----------+-------------------------------------------------------+
              | psi                          | rad       | polarization angle                                    |
              +------------------------------+-----------+-------------------------------------------------------+
              | theta_jn                     | rad       | inclination angle                                     |
              +------------------------------+-----------+-------------------------------------------------------+
              | a_1                          |           | spin of the primary compact binary                    |
              +------------------------------+-----------+-------------------------------------------------------+
              | a_2                          |           | spin of the secondary compact binary                  |
              +------------------------------+-----------+-------------------------------------------------------+
              | mass_1_source                | Msun      | mass of the primary compact binary (source frame)     |
              +------------------------------+-----------+-------------------------------------------------------+
              | mass_2_source                | Msun      | mass of the secondary compact binary (source frame)   |
              +------------------------------+-----------+-------------------------------------------------------+
              | mass_1                       | Msun      | mass of the primary compact binary (detector frame)   |
              +------------------------------+-----------+-------------------------------------------------------+
              | mass_2                       | Msun      | mass of the secondary compact binary (detector frame) |
              +------------------------------+-----------+-------------------------------------------------------+










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> params = lens.sample_lens_parameters(size=1000)
      >>> print(params.keys())



      ..
          !! processed by numpydoc !!

   .. py:method:: sample_all_routine_epl_shear_sl(size=1000)

      
      Sample EPL+shear galaxy lens parameters with strong lensing condition.


      :Parameters:

          **size** : ``int``
              Number of lens parameters to sample.

              default: 1000

      :Returns:

          **lens_parameters** : ``dict``
              Dictionary of sampled lens parameters.

              The included parameters and their units are as follows (for default settings):

              +------------------------------+-----------+-------------------------------------------------------+
              | Parameter                    | Units     | Description                                           |
              +==============================+===========+=======================================================+
              | zl                           |           | redshift of the lens                                  |
              +------------------------------+-----------+-------------------------------------------------------+
              | zs                           |           | redshift of the source                                |
              +------------------------------+-----------+-------------------------------------------------------+
              | sigma                        | km s^-1   | velocity dispersion                                   |
              +------------------------------+-----------+-------------------------------------------------------+
              | q                            |           | axis ratio                                            |
              +------------------------------+-----------+-------------------------------------------------------+
              | theta_E                      | arcsec    | Einstein radius                                       |
              +------------------------------+-----------+-------------------------------------------------------+
              | phi                          | rad       | axis rotation angle. counter-clockwise from the       |
              |                              |           | positive x-axis (RA-like axis) to the major axis of   |
              |                              |           | the projected mass distribution.                      |
              +------------------------------+-----------+-------------------------------------------------------+
              | gamma                        |           | density profile slope of EPL galaxy                   |
              +------------------------------+-----------+-------------------------------------------------------+
              | gamma1                       |           | external shear component in the x-direction           |
              |                              |           | (RA-like axis)                                        |
              +------------------------------+-----------+-------------------------------------------------------+
              | gamma2                       |           | external shear component in the y-direction           |
              |                              |           | (Dec-like axis)                                       |
              +------------------------------+-----------+-------------------------------------------------------+













      ..
          !! processed by numpydoc !!

   .. py:method:: strongly_lensed_source_redshifts(size=1000)

      
      Sample source redshifts conditioned on strong lensing.

      Uses rejection sampling to generate source redshifts from the CBC source
      population weighted by the optical depth, which increases with redshift.

      :Parameters:

          **size** : ``int``
              Number of redshifts to sample.

              default: 1000

      :Returns:

          **redshifts** : ``numpy.ndarray``
              Array of source redshifts conditioned on strong lensing.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> zs = lens.strongly_lensed_source_redshifts(size=1000)
      >>> print(f"Mean source redshift: {zs.mean():.2f}")



      ..
          !! processed by numpydoc !!

   .. py:method:: sample_all_routine_epl_shear_intrinsic(size=1000)

      
      Sample EPL+shear galaxy lens parameters from intrinsic distributions.

      Samples lens parameters from their intrinsic distributions without
      applying strong lensing cross-section weighting.

      :Parameters:

          **size** : ``int``
              Number of lens parameters to sample.

              default: 1000

      :Returns:

          **lens_parameters** : ``dict``
              Dictionary of sampled lens parameters.

              The included parameters and their units are as follows (for default settings):

              +------------------------------+-----------+-------------------------------------------------------+
              | Parameter                    | Units     | Description                                           |
              +==============================+===========+=======================================================+
              | zl                           |           | redshift of the lens                                  |
              +------------------------------+-----------+-------------------------------------------------------+
              | zs                           |           | redshift of the source                                |
              +------------------------------+-----------+-------------------------------------------------------+
              | sigma                        | km s^-1   | velocity dispersion                                   |
              +------------------------------+-----------+-------------------------------------------------------+
              | q                            |           | axis ratio                                            |
              +------------------------------+-----------+-------------------------------------------------------+
              | theta_E                      | arcsec    | Einstein radius                                       |
              +------------------------------+-----------+-------------------------------------------------------+
              | phi                          | rad       | axis rotation angle. counter-clockwise from the       |
              |                              |           | positive x-axis (RA-like axis) to the major axis of   |
              |                              |           | the projected mass distribution.                      |
              +------------------------------+-----------+-------------------------------------------------------+
              | gamma                        |           | density profile slope of EPL galaxy                   |
              +------------------------------+-----------+-------------------------------------------------------+
              | gamma1                       |           | external shear component in the x-direction           |
              |                              |           | (RA-like axis)                                        |
              +------------------------------+-----------+-------------------------------------------------------+
              | gamma2                       |           | external shear component in the y-direction           |
              |                              |           | (Dec-like axis)                                       |
              +------------------------------+-----------+-------------------------------------------------------+













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

.. py:function:: save_json(file_name, param)

   
   Save a dictionary as a json file.


   :Parameters:

       **file_name** : `str`
           json file name for storing the parameters.

       **param** : `dict`
           dictionary to be saved as a json file.














   ..
       !! processed by numpydoc !!

.. py:function:: load_json(file_name)

   
   Load a json file.


   :Parameters:

       **file_name** : `str`
           json file name for storing the parameters.

   :Returns:

       **param** : `dict`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: interpolator_json_path(identifier_dict, directory, sub_directory, interpolator_name)

   
   Function to create the interpolator json file path.


   :Parameters:

       **identifier_dict** : `dict`
           dictionary of identifiers.

       **directory** : `str`
           directory to store the interpolator.

       **sub_directory** : `str`
           sub-directory to store the interpolator.

       **interpolator_name** : `str`
           name of the interpolator.

   :Returns:

       **path_inv_cdf** : `str`
           path of the interpolator json file.

       **it_exist** : `bool`
           if True, the interpolator exists.













   ..
       !! processed by numpydoc !!

.. py:class:: FunctionConditioning(function=None, x_array=None, conditioned_y_array=None, y_array=None, non_zero_function=False, gaussian_kde=False, gaussian_kde_kwargs={}, identifier_dict={}, directory='./interpolator_json', sub_directory='default', name='default', create_new=False, create_function=False, create_function_inverse=False, create_pdf=False, create_rvs=False, multiprocessing_function=False, callback=None)


   .. py:attribute:: info

      

   .. py:attribute:: callback
      :value: 'None'

      

   .. py:method:: __call__(*args)


   .. py:method:: create_decision_function(create_function, create_function_inverse, create_pdf, create_rvs)


   .. py:method:: create_gaussian_kde(x_array, y_array, gaussian_kde_kwargs)


   .. py:method:: create_interpolator(function, x_array, conditioned_y_array, create_function_inverse, create_pdf, create_rvs, multiprocessing_function)


   .. py:method:: create_z_array(x_array, function, conditioned_y_array, create_pdf, create_rvs, multiprocessing_function)


   .. py:method:: cdf_values_generator(x_array, z_array, conditioned_y_array)


   .. py:method:: pdf_norm_const_generator(x_array, function_spline, conditioned_y_array)


   .. py:method:: function_spline_generator(x_array, z_array, conditioned_y_array)



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

.. py:function:: comoving_distance(z=None, z_min=0.001, z_max=10.0, cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7), directory='./interpolator_json', create_new=False, resolution=500, get_attribute=True)


.. py:function:: angular_diameter_distance(z=None, z_min=0.001, z_max=10.0, cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7), directory='./interpolator_json', create_new=False, resolution=500, get_attribute=True)


.. py:function:: angular_diameter_distance_z1z2(z1=None, z2=None, z_min=0.001, z_max=10.0, cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7), directory='./interpolator_json', create_new=False, resolution=500, get_attribute=True)


.. py:function:: differential_comoving_volume(z=None, z_min=0.001, z_max=10.0, cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7), directory='./interpolator_json', create_new=False, resolution=500, get_attribute=True)


.. py:function:: is_njitted(func)


.. py:function:: redshift_optimal_spacing(z_min, z_max, resolution)


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

.. py:function:: phi_q2_ellipticity(phi, q)

   
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

.. py:function:: axis_ratio_rayleigh_pdf(q, sigma, q_min=0.2, q_max=1.0)


.. py:class:: OpticalDepth(npool=4, z_min=0.0, z_max=10.0, cosmology=None, lens_type='epl_shear_galaxy', lens_functions=None, lens_functions_params=None, lens_param_samplers=None, lens_param_samplers_params=None, directory='./interpolator_json', create_new_interpolator=False, verbose=False)


   
   Class for computing optical depth and lens galaxy population parameters.

   This class calculates strong lensing optical depth, velocity dispersion,
   axis ratio, and other parameters for a lens galaxy population. It supports
   SIS, SIE, and EPL + external shear lens models with customizable samplers
   and interpolators for efficient computation.

   Key Features:

   - Multiple lens model support (SIS, SIE, EPL + shear)

   - Configurable velocity dispersion distributions

   - Cached interpolators for fast optical depth computation

   - Flexible parameter sampling with user-defined priors

   :Parameters:

       **npool** : ``int``
           Number of processors for multiprocessing.

           default: 4

       **z_min** : ``float``
           Minimum redshift of the lens galaxy population.

           default: 0.0

       **z_max** : ``float``
           Maximum redshift of the lens galaxy population.

           default: 10.0

       **cosmology** : ``astropy.cosmology`` or ``None``
           Cosmology object for distance calculations.

           default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

       **lens_type** : ``str``
           Type of lens galaxy model.

           Options:

           - 'epl_shear_galaxy': Elliptical power-law with external shear

           - 'sie_galaxy': Singular isothermal ellipsoid

           - 'sis_galaxy': Singular isothermal sphere

           default: 'epl_shear_galaxy'

       **lens_functions** : ``dict`` or ``None``
           Dictionary with lens-related functions.

           default: None (uses defaults for lens_type)

       **lens_functions_params** : ``dict`` or ``None``
           Dictionary with parameters for lens-related functions.

           default: None

       **lens_param_samplers** : ``dict`` or ``None``
           Dictionary of sampler functions for lens parameters.

           default: None (uses defaults for lens_type)

       **lens_param_samplers_params** : ``dict`` or ``None``
           Dictionary with parameters for the samplers.

           default: None

       **directory** : ``str``
           Directory where interpolators are saved.

           default: './interpolator_json'

       **create_new_interpolator** : ``bool`` or ``dict``
           Whether to create new interpolators.

           default: False

       **verbose** : ``bool``
           If True, prints additional information.

           default: False











   .. rubric:: Examples

   Basic usage:

   >>> from ler.lens_galaxy_population import OpticalDepth
   >>> od = OpticalDepth()
   >>> tau = od.optical_depth(zs=np.array([1.0, 2.0]))

   Instance Methods
   ----------
   OpticalDepth has the following instance methods:

   +-----------------------------------------------------+----------------------------------------------------------+
   | Method                                              | Description                                              |
   +=====================================================+==========================================================+
   | :meth:`~axis_ratio_rayleigh`                        | Sample axis ratio from Rayleigh distribution             |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~axis_ratio_padilla_strauss`                 | Sample axis ratio from Padilla & Strauss 2008            |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~axis_ratio_uniform`                         | Sample axis ratio from uniform distribution              |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~axis_rotation_angle_uniform`                | Sample axis rotation angle from uniform distribution     |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~lens_redshift_strongly_lensed_numerical`    | Sample lens redshift for strong lensing                  |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~lens_redshift_strongly_lensed_sis_haris`                    | Sample SIS lens redshift (Haris et al. 2018)             |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~velocity_dispersion_gengamma`               | Sample velocity dispersion from gengamma distribution    |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~velocity_dispersion_bernardi`               | Sample velocity dispersion (Bernardi et al. 2010)        |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~velocity_dispersion_ewoud`                  | Sample redshift-dependent velocity dispersion            |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~external_shear_normal`                      | Sample external shear from normal distribution           |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~density_profile_slope_normal`               | Sample density profile slope from normal distribution    |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~optical_depth_sis_analytic`                    | Compute SIS optical depth (Haris et al. 2018)            |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~cross_section_sis`                          | Compute SIS cross-section                                |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~cross_section_sie_feixu`                    | Compute SIE cross-section (Fei Xu et al. 2021)           |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~cross_section_epl_shear_numerical`          | Compute EPL+shear cross-section numerically              |
   +-----------------------------------------------------+----------------------------------------------------------+
   | :meth:`~cross_section_epl_shear_interpolation`      | Compute EPL+shear cross-section via interpolation        |
   +-----------------------------------------------------+----------------------------------------------------------+

   Instance Attributes
   ----------
   OpticalDepth has the following instance attributes:

   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | Attribute                                      | Type                         | Unit  | Description                                              |
   +================================================+==============================+=======+==========================================================+
   | :attr:`~npool`                                 | ``int``                      |       | Number of processors for multiprocessing                 |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~z_min`                                 | ``float``                    |       | Minimum redshift                                         |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~z_max`                                 | ``float``                    |       | Maximum redshift                                         |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~cosmo`                                 | ``astropy.cosmology``        |       | Cosmology object                                         |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~lens_type`                             | ``str``                      |       | Type of lens galaxy model                                |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~directory`                             | ``str``                      |       | Directory for interpolator storage                       |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~optical_depth`                         | ``FunctionConditioning``     |       | Optical depth calculator                                 |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~velocity_dispersion`                   | ``FunctionConditioning``     | km/s  | Velocity dispersion sampler                              |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~axis_ratio`                            | ``FunctionConditioning``     |       | Axis ratio sampler                                       |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~axis_rotation_angle`                   | ``FunctionConditioning``     | rad   | Axis rotation angle sampler                              |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~lens_redshift`                         | ``FunctionConditioning``     |       | Lens redshift sampler                                    |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~external_shear`                        | ``FunctionConditioning``     |       | External shear sampler                                   |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~density_profile_slope`                 | ``FunctionConditioning``     |       | Density profile slope sampler                            |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~cross_section`                         | ``callable``                 | rad²  | Cross-section calculator                                 |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~available_lens_samplers`               | ``dict``                     |       | Available lens parameter samplers                        |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
   | :attr:`~available_lens_functions`              | ``dict``                     |       | Available lens functions                                 |
   +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:property:: lens_type

      
      Type of lens galaxy model.



      :Returns:

          **lens_type** : ``str``
              Lens type ('epl_shear_galaxy', 'sie_galaxy', or 'sis_galaxy').













      ..
          !! processed by numpydoc !!

   .. py:property:: npool

      
      Number of processors for multiprocessing.



      :Returns:

          **npool** : ``int``
              Number of parallel processors.













      ..
          !! processed by numpydoc !!

   .. py:property:: z_min

      
      Minimum redshift of the lens galaxy population.



      :Returns:

          **z_min** : ``float``
              Minimum redshift.













      ..
          !! processed by numpydoc !!

   .. py:property:: z_max

      
      Maximum redshift of the lens galaxy population.



      :Returns:

          **z_max** : ``float``
              Maximum redshift.













      ..
          !! processed by numpydoc !!

   .. py:property:: cosmo

      
      Cosmology object for distance calculations.



      :Returns:

          **cosmo** : ``astropy.cosmology``
              Cosmology object.













      ..
          !! processed by numpydoc !!

   .. py:property:: directory

      
      Directory for interpolator storage.



      :Returns:

          **directory** : ``str``
              Path to interpolator JSON files.













      ..
          !! processed by numpydoc !!

   .. py:property:: velocity_dispersion

      
      Velocity dispersion sampler object.

      Returns a ``FunctionConditioning`` object with methods:

      - ``rvs(size, zl)``: Sample velocity dispersion values

      - ``pdf(sigma, zl)``: Get probability density

      - ``function(sigma, zl)``: Get number density function


      :Returns:

          **velocity_dispersion** : ``FunctionConditioning``
              Sampler object for velocity dispersion (km/s).










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> sigma = od.velocity_dispersion(size=100, zl=np.ones(100)*0.5)



      ..
          !! processed by numpydoc !!

   .. py:property:: axis_ratio

      
      Axis ratio sampler object.

      Returns a ``FunctionConditioning`` object with methods:

      - ``rvs(size, sigma)``: Sample axis ratio values

      - ``pdf(q, sigma)``: Get probability density

      - ``function(q, sigma)``: Get distribution function


      :Returns:

          **axis_ratio** : ``FunctionConditioning``
              Sampler object for axis ratio.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> q = od.axis_ratio(size=100, sigma=np.ones(100)*200.)



      ..
          !! processed by numpydoc !!

   .. py:property:: axis_rotation_angle

      
      Axis rotation angle sampler object.

      Returns a ``FunctionConditioning`` object with methods:

      - ``rvs(size)``: Sample axis rotation angles

      - ``pdf(phi)``: Get probability density


      :Returns:

          **axis_rotation_angle** : ``FunctionConditioning``
              Sampler object for axis rotation angle (rad).










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> phi = od.axis_rotation_angle(size=100)



      ..
          !! processed by numpydoc !!

   .. py:property:: density_profile_slope

      
      Density profile slope sampler object.

      Returns a ``FunctionConditioning`` object with methods:

      - ``rvs(size)``: Sample density profile slope values

      - ``pdf(gamma)``: Get probability density


      :Returns:

          **density_profile_slope** : ``FunctionConditioning``
              Sampler object for density profile slope.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> gamma = od.density_profile_slope(size=100)



      ..
          !! processed by numpydoc !!

   .. py:property:: external_shear

      
      External shear sampler object.

      Returns a ``FunctionConditioning`` object with methods:

      - ``rvs(size)``: Sample shear components (gamma1, gamma2)

      - ``pdf(gamma1, gamma2)``: Get probability density


      :Returns:

          **external_shear** : ``FunctionConditioning``
              Sampler object for external shear.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> gamma1, gamma2 = od.external_shear(size=100)



      ..
          !! processed by numpydoc !!

   .. py:property:: cross_section

      
      Lensing cross-section calculator.

      Returns a callable that computes lensing cross-section for individual

      lensing events. Input parameters depend on lens type:

      - EPL+shear: zs, zl, sigma, q, phi, gamma, gamma1, gamma2

      - SIE: zs, zl, sigma, q

      - SIS: zs, zl, sigma


      :Returns:

          **cross_section** : ``callable``
              Cross-section function (rad² units).










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> cs = od.cross_section(zs=zs, zl=zl, sigma=sigma, ...)



      ..
          !! processed by numpydoc !!

   .. py:property:: lens_redshift

      
      Lens redshift sampler object.

      Returns a ``FunctionConditioning`` object with methods:

      - ``rvs(size, zs)``: Sample lens redshifts given source redshifts

      - ``pdf(zl, zs)``: Get probability density

      - ``function(zl, zs)``: Get effective lensing cross-section


      :Returns:

          **lens_redshift** : ``FunctionConditioning``
              Sampler object for lens redshift.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> zl = od.lens_redshift(size=100, zs=np.ones(100)*2.0)



      ..
          !! processed by numpydoc !!

   .. py:property:: density_profile_slope_sl

      
      Density profile slope sampler object (strong lensing conditioned).

      Returns a ``FunctionConditioning`` object with methods:

      - ``rvs(size)``: Sample density profile slope values

      - ``pdf(gamma)``: Get probability density


      :Returns:

          **density_profile_slope_sl** : ``FunctionConditioning``
              Sampler object for density profile slope (strong lensing).










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> gamma = od.density_profile_slope_sl(size=100)



      ..
          !! processed by numpydoc !!

   .. py:property:: external_shear_sl

      
      External shear sampler object (strong lensing conditioned).

      Returns a ``FunctionConditioning`` object with methods:

      - ``rvs(size)``: Sample shear components (gamma1, gamma2)

      - ``pdf(gamma1, gamma2)``: Get probability density


      :Returns:

          **external_shear_sl** : ``FunctionConditioning``
              Sampler object for external shear (strong lensing).










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> gamma1, gamma2 = od.external_shear_sl(size=100)



      ..
          !! processed by numpydoc !!

   .. py:property:: optical_depth

      
      Strong lensing optical depth calculator.



      :Returns:

          **optical_depth** : ``FunctionConditioning``
              Function object with `.function(zs)` method that returns \n
              optical depth for given source redshifts.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> tau = od.optical_depth.function(np.array([1.0, 2.0]))



      ..
          !! processed by numpydoc !!

   .. py:property:: available_lens_samplers

      
      Dictionary of available lens parameter samplers and their default parameters.



      :Returns:

          **available_lens_samplers** : ``dict``
              Dictionary with sampler names and default parameters.













      ..
          !! processed by numpydoc !!

   .. py:property:: available_lens_functions

      
      Dictionary of available lens functions and their default parameters.



      :Returns:

          **available_lens_functions** : ``dict``
              Dictionary with function names and default parameters.













      ..
          !! processed by numpydoc !!

   .. py:attribute:: comoving_distance

      

   .. py:attribute:: angular_diameter_distance

      

   .. py:attribute:: angular_diameter_distance_z1z2

      

   .. py:attribute:: differential_comoving_volume

      

   .. py:attribute:: lens_redshift_intrinsic

      

   .. py:method:: axis_ratio_rayleigh(size, sigma, get_attribute=False, **kwargs)

      
      Sample axis ratio from Rayleigh distribution conditioned on velocity dispersion.


      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **sigma** : ``numpy.ndarray``
              Velocity dispersion of the lens galaxy (km/s).

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters (q_min, q_max).

      :Returns:

          **q** : ``numpy.ndarray`` or ``FunctionConditioning``
              Axis ratio samples or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(axis_ratio="axis_ratio_rayleigh"))
      >>> q = od.axis_ratio(size=100, sigma=np.ones(100)*200.)



      ..
          !! processed by numpydoc !!

   .. py:method:: axis_ratio_padilla_strauss(size=1000, get_attribute=False, **kwargs)

      
      Sample axis ratio from Padilla & Strauss (2008) distribution.


      :Parameters:

          **size** : ``int``
              Number of samples to generate.

              default: 1000

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters (q_min, q_max).

      :Returns:

          **q** : ``numpy.ndarray`` or ``FunctionConditioning``
              Axis ratio samples or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(axis_ratio="axis_ratio_padilla_strauss"))
      >>> q = od.axis_ratio(size=100)



      ..
          !! processed by numpydoc !!

   .. py:method:: lens_redshift_strongly_lensed_numerical(size=1000, zs=None, get_attribute=False, **kwargs)

      
      Sample lens redshifts conditioned on strong lensing (numerical method).

      This method computes the lens redshift distribution by numerically
      integrating over the velocity dispersion distribution (galaxy density distribution wrt), cross-section and differential comoving volume.

      :Parameters:

          **size** : ``int``
              Number of samples to generate. \n
              default: 1000

          **zs** : ``numpy.ndarray``
              Source redshifts.

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples. \n
              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters.

      :Returns:

          **zl** : ``numpy.ndarray`` or ``FunctionConditioning``
              Lens redshift samples or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> zl = od.lens_redshift(size=100, zs=np.ones(100)*2.0)



      ..
          !! processed by numpydoc !!

   .. py:method:: axis_rotation_angle_uniform(size, get_attribute=False, **kwargs)

      
      Sample axis rotation angle from uniform distribution.


      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters (phi_min, phi_max).

      :Returns:

          **phi** : ``numpy.ndarray`` or ``FunctionConditioning``
              Axis rotation angle samples (rad) or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> phi = od.axis_rotation_angle(size=100)



      ..
          !! processed by numpydoc !!

   .. py:method:: axis_ratio_uniform(size, get_attribute=False, **kwargs)

      
      Sample axis ratio from uniform distribution.


      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters (q_min, q_max).

      :Returns:

          **q** : ``numpy.ndarray`` or ``FunctionConditioning``
              Axis ratio samples or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(axis_ratio="axis_ratio_uniform"))
      >>> q = od.axis_ratio(size=100)



      ..
          !! processed by numpydoc !!

   .. py:method:: external_shear_normal(size, get_attribute=False, **kwargs)

      
      Sample external shear parameters from 2D normal distribution.


      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters (mean, std).

      :Returns:

          **shear** : ``numpy.ndarray`` or ``FunctionConditioning``
              Array of shape (2, size) with gamma1, gamma2 or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> gamma1, gamma2 = od.external_shear(size=100)



      ..
          !! processed by numpydoc !!

   .. py:method:: density_profile_slope_normal(size, get_attribute=False, **kwargs)

      
      Sample density profile slope from normal distribution.


      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters (mean, std).

      :Returns:

          **gamma** : ``numpy.ndarray`` or ``FunctionConditioning``
              Density profile slope samples or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> gamma = od.density_profile_slope(size=100)



      ..
          !! processed by numpydoc !!

   .. py:method:: lens_redshift_strongly_lensed_sis_haris(size, zs, get_attribute=False, **kwargs)

      
      Sample SIS lens redshifts using Haris et al. (2018) distribution.


      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **zs** : ``numpy.ndarray``
              Source redshifts.

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters.

      :Returns:

          **zl** : ``numpy.ndarray`` or ``FunctionConditioning``
              Lens redshift samples or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_type='sis_galaxy')
      >>> zl = od.lens_redshift(size=100, zs=np.ones(100)*2.0)



      ..
          !! processed by numpydoc !!

   .. py:method:: velocity_dispersion_gengamma(size, get_attribute=False, **kwargs)

      
      Sample velocity dispersion from generalized gamma distribution.


      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters (a, c for gengamma distribution).

              See scipy.stats.gengamma for details.

      :Returns:

          **sigma** : ``numpy.ndarray`` or ``FunctionConditioning``
              Velocity dispersion samples (km/s) or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(
      ...     velocity_dispersion="velocity_dispersion_gengamma"))
      >>> sigma = od.velocity_dispersion(size=100)



      ..
          !! processed by numpydoc !!

   .. py:method:: velocity_dispersion_bernardi(size, get_attribute=False, **kwargs)

      
      Sample velocity dispersion from Bernardi et al. (2010) distribution.

      Uses inverse transform sampling on the velocity dispersion function.

      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters.

      :Returns:

          **sigma** : ``numpy.ndarray`` or ``FunctionConditioning``
              Velocity dispersion samples (km/s) or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth(lens_param_samplers=dict(
      ...     velocity_dispersion="velocity_dispersion_bernardi"))
      >>> sigma = od.velocity_dispersion(size=100)



      ..
          !! processed by numpydoc !!

   .. py:method:: velocity_dispersion_ewoud(size, zl, get_attribute=False, **kwargs)

      
      Sample redshift-dependent velocity dispersion from Wempe et al. (2022).

      Uses inverse transform sampling with redshift-dependent velocity
      dispersion function.

      :Parameters:

          **size** : ``int``
              Number of samples to generate.

          **zl** : ``numpy.ndarray``
              Lens redshifts.

          **get_attribute** : ``bool``
              If True, returns the sampler object instead of samples.

              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters.

      :Returns:

          **sigma** : ``numpy.ndarray`` or ``FunctionConditioning``
              Velocity dispersion samples (km/s) or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> sigma = od.velocity_dispersion(size=100, zl=np.ones(100)*0.5)



      ..
          !! processed by numpydoc !!

   .. py:method:: optical_depth_numerical(zs, get_attribute=False, **kwargs)

      
      Helper to compute optical depth numerically by integrating lens redshift.


      :Parameters:

          **zs** : ``numpy.ndarray``
              Source redshifts.

          **get_attribute** : ``bool``
              If True, returns the function object instead of values. \n
              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters.

      :Returns:

          **tau** : ``numpy.ndarray`` or ``FunctionConditioning``
              Optical depth values or function object.













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

   .. py:method:: optical_depth_sis_analytic(zs, get_attribute=False, **kwargs)

      
      Function to compute the strong lensing optical depth (SIS).
      LambdaCDM(H0=70, Om0=0.3, Ode0=0.7) was used to derive the following equation. This is the analytic version of optical depth from z=0 to z=zs.


      :Parameters:

          **zs** : ``numpy.ndarray``
              Source redshifts.

          **get_attribute** : ``bool``
              If True, returns the function object instead of values. \n
              default: False

          **\*\*kwargs** : ``dict``
              Additional parameters.

      :Returns:

          **tau** : ``numpy.ndarray`` or ``FunctionConditioning``
              Optical depth values or function object.













      ..
          !! processed by numpydoc !!

   .. py:method:: cross_section_sis(zs=None, zl=None, sigma=None, get_attribute=False, **kwargs)

      
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
      >>> print(self.cross_section_sis(sigma=200., zl=0.5, zs=1.0))



      ..
          !! processed by numpydoc !!

   .. py:method:: cross_section_sie_feixu(zs=None, zl=None, sigma=None, q=None, get_attribute=False, **kwargs)

      
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
      >>> print(self.cross_section_sie_feixu(sigma=200., zl=0.5, zs=1.0, q=1.0))



      ..
          !! processed by numpydoc !!

   .. py:method:: cross_section_epl_shear_numerical(theta_E, gamma, gamma1, gamma2, q=None, phi=None, e1=None, e2=None, verbose=False, **kwargs)

      
      Function to compute the strong lensing cross-section numerically for EPL + external shear lenses.


      :Parameters:

          **theta_E** : `numpy.ndarray`
              Einstein radii of the lens galaxies in radians

          **gamma** : `numpy.ndarray`
              external shear magnitudes of the lens galaxies

          **gamma1** : `numpy.ndarray`
              external shear component 1 of the lens galaxies

          **gamma2** : `numpy.ndarray`
              external shear component 2 of the lens galaxies

          **q** : `numpy.ndarray`
              axis ratios of the lens galaxies

          **phi** : `numpy.ndarray`
              axis rotation angles of the lens galaxies in radians

          **e1** : `numpy.ndarray`
              ellipticity component 1 of the lens galaxies

          **e2** : `numpy.ndarray`
              ellipticity component 2 of the lens galaxies

      :Returns:

          **cross_section** : `numpy.ndarray`
              strong lensing cross-section of the lens galaxies in square radians













      ..
          !! processed by numpydoc !!

   .. py:method:: create_parameter_grid(size_list=[25, 25, 45, 15, 15])


   .. py:method:: cross_section_epl_shear_interpolation_init(file_path, size_list)


   .. py:method:: cross_section_epl_shear_interpolation(zs, zl, sigma, q, phi, gamma, gamma1, gamma2, get_attribute=False, size_list=[25, 25, 45, 15, 15], **kwargs)

      
      Function to compute the cross-section correction factor
















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

.. py:function:: phi_q2_ellipticity(phi, q)

   
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

.. py:function:: load_json(file_name)

   
   Load a json file.


   :Parameters:

       **file_name** : `str`
           json file name for storing the parameters.

   :Returns:

       **param** : `dict`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: make_cross_section_reinit(e1_grid, e2_grid, gamma_grid, gamma1_grid, gamma2_grid, cs_grid, Da_instance, csunit_to_cs_slope=0.31830988618379075, csunit_to_cs_intercept=-3.2311742677852644e-27)

   
   Factory function to create a JIT-compiled cross section calculator.

   This function precomputes B-spline coefficients and creates a closure
   that captures the grid parameters, returning a fast Numba-compiled
   function for computing cross sections.

   :Parameters:

       **e1_grid** : ``numpy.ndarray``
           Grid values for ellipticity component e1, shape (n_e1,).

       **e2_grid** : ``numpy.ndarray``
           Grid values for ellipticity component e2, shape (n_e2,).

       **gamma_grid** : ``numpy.ndarray``
           Grid values for density slope gamma, shape (n_g,).

       **gamma1_grid** : ``numpy.ndarray``
           Grid values for shear component gamma1, shape (n_g1,).

       **gamma2_grid** : ``numpy.ndarray``
           Grid values for shear component gamma2, shape (n_g2,).

       **cs_grid** : ``numpy.ndarray``
           Raw cross section grid data (before spline filtering),

           shape (n_e1, n_e2, n_g, n_g1, n_g2).

       **Da_instance** : ``callable``
           Angular diameter distance function.

           Signature: ``Da_instance(z) -> distance``

       **csunit_to_cs_slope** : ``float``
           Slope for affine calibration from unit cross section.

           default: 0.31830988618379075

       **csunit_to_cs_intercept** : ``float``
           Intercept for affine calibration from unit cross section.

           default: -3.2311742677852644e-27

   :Returns:

       **cross_section_reinit** : ``callable``
           JIT-compiled function with signature:

           ``cross_section_reinit(zs, zl, sigma, q, phi, gamma, gamma1, gamma2)``

           Returns cross sections as ``numpy.ndarray`` of shape (N,).










   .. rubric:: Examples

   >>> from ler.lens_galaxy_population.cross_section_interpolator import make_cross_section_reinit
   >>> cs_func = make_cross_section_reinit(
   ...     e1_grid, e2_grid, gamma_grid, gamma1_grid, gamma2_grid,
   ...     cs_grid, Da_instance
   ... )
   >>> cross_sections = cs_func(zs, zl, sigma, q, phi, gamma, gamma1, gamma2)



   ..
       !! processed by numpydoc !!

.. py:data:: CS_UNIT_SLOPE
   :value: '0.31830988618379075'

   

.. py:data:: CS_UNIT_INTERCEPT
   :value: '-3.2311742677852644e-27'

   

.. py:function:: lens_redshift_strongly_lensed_njit(zs_array, zl_scaled, sigma_min, sigma_max, q_rvs, phi_rvs, gamma_rvs, shear_rvs, number_density, cross_section, dVcdz_function, integration_size)

   
   JIT-compiled parallel computation of lens redshift optical depth.

   Computes the differential optical depth for strong lensing as a function
   of source and lens redshifts using Monte Carlo integration with parallel
   execution over source redshifts.

   :Parameters:

       **zs_array** : ``numpy.ndarray``
           1D array of source redshifts.

       **zl_scaled** : ``numpy.ndarray``
           2D array of scaled lens redshifts (zl/zs).

       **sigma_min** : ``float``
           Minimum velocity dispersion (km/s).

       **sigma_max** : ``float``
           Maximum velocity dispersion (km/s).

       **q_rvs** : ``callable``
           Function to sample axis ratios given size and sigma.

       **phi_rvs** : ``callable``
           Function to sample axis rotation angles.

       **gamma_rvs** : ``callable``
           Function to sample density profile slopes.

       **shear_rvs** : ``callable``
           Function to sample external shear components (gamma1, gamma2).

       **number_density** : ``callable``
           Function to compute velocity dispersion number density.

       **cross_section** : ``callable``
           Function to compute lensing cross-section.

       **dVcdz_function** : ``callable``
           Function to compute differential comoving volume.

       **integration_size** : ``int``
           Number of Monte Carlo samples per (zs, zl) pair.

   :Returns:

       **result_array** : ``numpy.ndarray``
           2D array of optical depth values with shape (len(zs_array), len(zl_scaled[0])).













   ..
       !! processed by numpydoc !!

.. py:function:: lens_redshift_strongly_lensed_mp(params)

   
   Multiprocessing worker for lens redshift optical depth calculation.

   Computes the differential optical depth for a single source redshift
   across multiple lens redshifts. Designed to be called via multiprocessing
   Pool.map() for parallel computation.

   :Parameters:

       **params** : ``tuple``
           Packed parameters tuple containing:

           - params[0]: Source redshift (float)

           - params[1]: Scaled lens redshift array (1D array)

           - params[2]: Sigma args (min, max, spline_z, spline_z2d, spline_coeff)

           - params[3]: Axis ratio args (type, min/spline, max/coeff, coeff_2d)

           - params[4]: Angular diameter distance args (x_values, coeffs)

           - params[5]: dVcdz args (x_values, coeffs)

           - params[6]: Worker index (int)

           - params[7]: Cross-section args (type, json_path)

           - params[8]: Axis rotation args (type, min, max or spline)

           - params[9]: Shear args (type, loc, scale)

           - params[10]: Slope args (type, loc/spline, scale/coeff)

           - params[11]: Integration size (int)

   :Returns:

       **idx** : ``int``
           Worker index for result ordering.

       **result_array** : ``numpy.ndarray``
           1D array of optical depth values for each lens redshift.













   ..
       !! processed by numpydoc !!

.. py:function:: cross_section_unit_mp(params)

   
   Multiprocessing worker for unit Einstein radius cross-section.

   Computes the lensing cross-section for a lens with unit Einstein radius
   (theta_E = 1). Used for building cross-section interpolation grids.

   :Parameters:

       **params** : ``tuple``
           Packed parameters (e1, e2, gamma, gamma1, gamma2, idx).

   :Returns:

       **idx** : ``int``
           Worker index for result ordering.

       **area** : ``float``
           Cross-section area in square arcseconds.













   ..
       !! processed by numpydoc !!

.. py:function:: cross_section_mp(params)

   
   Multiprocessing worker for cross-section calculation.

   Computes the lensing cross-section for given lens parameters.
   Designed to be called via multiprocessing Pool.map().

   :Parameters:

       **params** : ``tuple``
           Packed parameters (theta_E, e1, e2, gamma, gamma1, gamma2, idx).

   :Returns:

       **idx** : ``int``
           Worker index for result ordering.

       **area** : ``float``
           Cross-section area in square arcseconds.













   ..
       !! processed by numpydoc !!

.. py:function:: cross_section(theta_E, e1, e2, gamma, gamma1, gamma2)

   
   Compute the strong lensing cross-section for an EPL+Shear lens.

   Uses lenstronomy to compute the caustic structure and returns the
   area enclosed by the double-image (outer) caustic.

   :Parameters:

       **theta_E** : ``float``
           Einstein radius in arcseconds.

       **e1** : ``float``
           First ellipticity component.

       **e2** : ``float``
           Second ellipticity component.

       **gamma** : ``float``
           Power-law density profile slope.

       **gamma1** : ``float``
           First external shear component.

       **gamma2** : ``float``
           Second external shear component.

   :Returns:

       **area** : ``float``
           Cross-section area in square arcseconds.

           Returns 0.0 if caustic computation fails.













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

.. py:function:: gamma_(x)


.. py:function:: cvdf_fit(log_vd, redshift)


.. py:function:: my_derivative(log_vd, redshift, dx)


.. py:function:: pdf_phi_z_div_0(s, z)


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

.. py:function:: axis_ratio_rayleigh_pdf(q, sigma, q_min=0.2, q_max=1.0)


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

.. py:function:: phi_q2_ellipticity(phi, q)

   
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

.. py:function:: sample_sigma_zl(pdf, sigma_min, sigma_max, zl_min, zl_max, zs, chunk_size=10000)


.. py:function:: phi_q2_ellipticity(phi, q)

   
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

.. py:data:: C_LIGHT
   :value: '299792.458'

   

.. py:function:: make_cross_section_reinit(e1_grid, e2_grid, gamma_grid, gamma1_grid, gamma2_grid, cs_grid, Da_instance, csunit_to_cs_slope=0.31830988618379075, csunit_to_cs_intercept=-3.2311742677852644e-27)

   
   Factory function to create a JIT-compiled cross section calculator.

   This function precomputes B-spline coefficients and creates a closure
   that captures the grid parameters, returning a fast Numba-compiled
   function for computing cross sections.

   :Parameters:

       **e1_grid** : ``numpy.ndarray``
           Grid values for ellipticity component e1, shape (n_e1,).

       **e2_grid** : ``numpy.ndarray``
           Grid values for ellipticity component e2, shape (n_e2,).

       **gamma_grid** : ``numpy.ndarray``
           Grid values for density slope gamma, shape (n_g,).

       **gamma1_grid** : ``numpy.ndarray``
           Grid values for shear component gamma1, shape (n_g1,).

       **gamma2_grid** : ``numpy.ndarray``
           Grid values for shear component gamma2, shape (n_g2,).

       **cs_grid** : ``numpy.ndarray``
           Raw cross section grid data (before spline filtering),

           shape (n_e1, n_e2, n_g, n_g1, n_g2).

       **Da_instance** : ``callable``
           Angular diameter distance function.

           Signature: ``Da_instance(z) -> distance``

       **csunit_to_cs_slope** : ``float``
           Slope for affine calibration from unit cross section.

           default: 0.31830988618379075

       **csunit_to_cs_intercept** : ``float``
           Intercept for affine calibration from unit cross section.

           default: -3.2311742677852644e-27

   :Returns:

       **cross_section_reinit** : ``callable``
           JIT-compiled function with signature:

           ``cross_section_reinit(zs, zl, sigma, q, phi, gamma, gamma1, gamma2)``

           Returns cross sections as ``numpy.ndarray`` of shape (N,).










   .. rubric:: Examples

   >>> from ler.lens_galaxy_population.cross_section_interpolator import make_cross_section_reinit
   >>> cs_func = make_cross_section_reinit(
   ...     e1_grid, e2_grid, gamma_grid, gamma1_grid, gamma2_grid,
   ...     cs_grid, Da_instance
   ... )
   >>> cross_sections = cs_func(zs, zl, sigma, q, phi, gamma, gamma1, gamma2)



   ..
       !! processed by numpydoc !!

.. py:function:: create_rejection_sampler(sigma_max, sigma_rvs, q_rvs, phi_rvs, gamma_rvs, shear_rvs, cross_section, safety_factor=1.2, use_njit_sampler=True)

   
   Create a rejection sampler for lens parameters weighted by cross section.

   Returns a callable that samples lens parameters using rejection sampling,
   optionally JIT-compiled for improved performance.

   :Parameters:

       **sigma_max** : ``float``
           Maximum velocity dispersion (km/s) for computing upper bound.

       **sigma_rvs** : ``callable``
           Function to sample velocity dispersion: sigma_rvs(n, zl) -> array.

       **q_rvs** : ``callable``
           Function to sample axis ratio: q_rvs(n, sigma) -> array.

       **phi_rvs** : ``callable``
           Function to sample orientation angle: phi_rvs(n) -> array.

       **gamma_rvs** : ``callable``
           Function to sample power-law index: gamma_rvs(n) -> array.

       **shear_rvs** : ``callable``
           Function to sample external shear: shear_rvs(n) -> (gamma1, gamma2).

       **cross_section** : ``callable``
           Function to compute cross section.

       **safety_factor** : ``float``
           Multiplicative safety factor for the upper bound.

           default: 1.2

       **use_njit_sampler** : ``bool``
           If True, uses Numba JIT compilation for faster execution.

           default: True

   :Returns:

       **rejection_sampler_wrapper** : ``callable``
           Function with signature (zs, zl) -> (sigma, q, phi, gamma, gamma1, gamma2).













   ..
       !! processed by numpydoc !!

.. py:function:: create_importance_sampler(sigma_min, sigma_max, q_rvs, phi_rvs, gamma_rvs, shear_rvs, sigma_pdf, cross_section, n_prop, use_njit_sampler=True)

   
   Create an importance sampler for lens parameters weighted by cross section.

   Returns a callable that samples lens parameters using importance sampling
   with uniform proposal, optionally JIT-compiled for improved performance.

   :Parameters:

       **sigma_min** : ``float``
           Minimum velocity dispersion (km/s) for uniform proposal.

       **sigma_max** : ``float``
           Maximum velocity dispersion (km/s) for uniform proposal.

       **q_rvs** : ``callable``
           Function to sample axis ratio: q_rvs(n, sigma) -> array.

       **phi_rvs** : ``callable``
           Function to sample orientation angle: phi_rvs(n) -> array.

       **gamma_rvs** : ``callable``
           Function to sample power-law index: gamma_rvs(n) -> array.

       **shear_rvs** : ``callable``
           Function to sample external shear: shear_rvs(n) -> (gamma1, gamma2).

       **sigma_pdf** : ``callable``
           PDF of velocity dispersion: sigma_pdf(sigma, zl) -> array.

       **cross_section** : ``callable``
           Function to compute cross section.

       **n_prop** : ``int``
           Number of proposal samples per lens.

       **use_njit_sampler** : ``bool``
           If True, uses Numba JIT compilation for faster execution.

           default: True

   :Returns:

       **importance_sampler_wrapper** : ``callable``
           Function with signature (zs, zl) -> (sigma, q, phi, gamma, gamma1, gamma2).













   ..
       !! processed by numpydoc !!

