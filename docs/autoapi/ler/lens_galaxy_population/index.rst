:py:mod:`ler.lens_galaxy_population`
====================================

.. py:module:: ler.lens_galaxy_population


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   cross_section_interpolator/index.rst
   epl_shear_cross_section/index.rst
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
   ler.lens_galaxy_population.phi_q2_ellipticity_hemanta
   ler.lens_galaxy_population.axis_ratio_rayleigh_pdf
   ler.lens_galaxy_population.cross_section_mp
   ler.lens_galaxy_population.phi_cut_SIE
   ler.lens_galaxy_population.phi_q2_ellipticity_hemanta
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
   ler.lens_galaxy_population.phi_q2_ellipticity_hemanta
   ler.lens_galaxy_population.sample_sigma_zl
   ler.lens_galaxy_population.phi_q2_ellipticity_hemanta
   ler.lens_galaxy_population.cross_section
   ler.lens_galaxy_population.make_cross_section_reinit
   ler.lens_galaxy_population.rejection_sampler
   ler.lens_galaxy_population.create_rejection_sampler
   ler.lens_galaxy_population.sigma_proposal_uniform
   ler.lens_galaxy_population.weighted_choice_1d
   ler.lens_galaxy_population.importance_sampler
   ler.lens_galaxy_population.create_importance_sampler



Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.C_LIGHT
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
   | :meth:`~binary_masses_BBH_popI_II_powerlaw_gaussian`| Sample BBH masses with PowerLaw+PEAK model     |
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

   .. py:method:: binary_masses_BBH_popI_II_powerlaw_gaussian(size, get_attribute=False, **kwargs)

      
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
      >>> m1_src, m2_src = cbc.binary_masses_BBH_popI_II_powerlaw_gaussian(size=1000)



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
           Directory where the interpolators are saved (default is './interpolator_json').
           If True, creates a new interpolator (default is False).

       **verbose** : bool, optional
           If True, prints additional information during initialization (default is False).





   :Raises:

       ValueError
           If `lens_type` is not in ['sie_galaxy', 'epl_shear_galaxy', 'sis_galaxy'].









   ..
       !! processed by numpydoc !!
   .. py:property:: velocity_dispersion

      
      Class object (of FunctionConditioning) for velocity dispersion of lens galaxy, with rvs/sampler as callback. Lens redshift `zl` is required only if velocity dispersion is redshift dependent. The class object contains the following attribute methods:
      - `rvs`: returns random samples from the velocity dispersion distribution
      - `pdf`: returns the probability density function of the velocity dispersion distribution
      - `function`: returns the velocity dispersion distribution function which represents the number density of lens galaxies as a function of velocity dispersion


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **zl** : `float`
              redshift of the lens galaxy. Should be of shape (size,)

      :Returns:

          **sigma** : `numpy.ndarray`
              velocity dispersion of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> import numpy as np
      >>> ler = OpticalDepth()
      >>> size = 10
      >>> zl = np.ones(size)*1.5
      >>> print(ler.velocity_dispersion(size, zl))



      ..
          !! processed by numpydoc !!

   .. py:property:: axis_ratio

      
      Class object (of FunctionConditioning) for axis ratio of lens galaxy, with rvs/sampler as callback. Velocity dispersion `sigma` is required only if axis ratio is velocity dispersion dependent. The class object contains the following attribute methods:
      - `rvs`: returns random samples from the axis ratio distribution
      - `pdf`: returns the probability density function of the axis ratio distribution
      - `function`: returns the un-normalized axis ratio distribution function.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **sigma** : `numpy.ndarray`
              velocity dispersion of the lens galaxy. Should be of shape (size,)

      :Returns:

          **q** : `numpy.ndarray`
              axis ratio of the lens galaxy.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> import numpy as np
      >>> ler = OpticalDepth()
      >>> size = 10
      >>> sigma = np.ones(size)*150.0
      >>> print(ler.axis_ratio(size, sigma))



      ..
          !! processed by numpydoc !!

   .. py:property:: axis_rotation_angle

      
      Class object (of FunctionConditioning) for axis rotation angle of lens galaxy, with rvs/sampler as callback. The class object contains the following attribute methods:
      - `rvs`: returns random samples from the axis rotation angle distribution
      - `pdf`: returns the probability density function of the axis rotation angle distribution
      - `function`: returns the axis rotation angle distribution function.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **theta** : `numpy.ndarray`
              axis rotation angle of the lens galaxy.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> import numpy as np
      >>> ler = OpticalDepth()
      >>> size = 10
      >>> print(ler.axis_rotation_angle(size))



      ..
          !! processed by numpydoc !!

   .. py:property:: density_profile_slope

      
      Class object (of FunctionConditioning) for density profile slope of lens galaxy, with rvs/sampler as callback. The class object contains the following attribute methods:
      - `rvs`: returns random samples from the density profile slope distribution
      - `pdf`: returns the probability density function of the density profile slope distribution
      - `function`: returns the density profile slope distribution function.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **gamma** : `numpy.ndarray`
              density profile slope of the lens galaxy.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> import numpy as np
      >>> ler = OpticalDepth()
      >>> size = 10
      >>> print(ler.density_profile_slope(size))



      ..
          !! processed by numpydoc !!

   .. py:property:: external_shear

      
      Class object (of FunctionConditioning) for external shear of lens galaxy, with rvs/sampler as callback. The class object contains the following attribute methods:
      - `rvs`: returns random samples from the external shear distribution
      - `pdf`: returns the probability density function of the external shear distribution
      - `function`: returns the external shear distribution function.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **gamma1** : `numpy.ndarray`
              external shear of the lens galaxy.

          **gamma2** : `numpy.ndarray`
              external shear of the lens galaxy.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> import numpy as np
      >>> ler = OpticalDepth()
      >>> size = 10
      >>> print(ler.external_shear(size))



      ..
          !! processed by numpydoc !!

   .. py:property:: cross_section

      
      Lensing cross section for individual lensing events. It can be of the following lens type and corresponding input parameters:
      - `epl_shear_galaxy`: `zs`, `zl`, `sigma`, `q`, `phi`, `gamma`, `gamma1`, `gamma2`
      - `sie_galaxy`: `zs`, `zl`, `sigma`, `q`
      - `sis_galaxy`: `zs`, `zl`, `sigma`


      :Parameters:

          **zs** : `numpy.ndarray`
              Redshift of the source

          **zl** : `numpy.ndarray`
              Redshift of the lens

          **sigma** : `numpy.ndarray`
              Angular size of the lens

          **q** : `numpy.ndarray`
              Axis ratio of the lens

          **phi** : `numpy.ndarray`
              Position angle of the lens

          **gamma** : `numpy.ndarray`
              density profile slope of the lens

          **gamma1** : `numpy.ndarray`
              Shear of the lens (x-direction)

          **gamma2** : `numpy.ndarray`
              Shear of the lens (y-direction)

      :Returns:

          **cross_section** : `numpy.ndarray`
              Lensing cross section for individual lensing events













      ..
          !! processed by numpydoc !!

   .. py:property:: lens_redshift

      
      Class object (of FunctionConditioning) for lens redshift, with rvs/sampler as callback. The class object contains the following attribute methods:
      - `rvs`: returns random samples from the lens redshift distribution
      - `pdf`: returns the probability density function of the lens redshift distribution
      - `function`: returns the lens redshift distribution function which represents effective lensing cross-section for lenses at redshift zl,


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **zs** : `numpy.ndarray`
              redshift of the lens galaxy. Should be of shape (size,)

      :Returns:

          **zl** : `numpy.ndarray`
              redshift of the lens galaxy.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> import numpy as np
      >>> ler = OpticalDepth()
      >>> size = 10
      >>> zs = np.ones(size)*1.5
      >>> print(ler.lens_redshift(size, zs))



      ..
          !! processed by numpydoc !!

   .. py:property:: density_profile_slope_sl

      
      Class object (of FunctionConditioning) for density profile slope of lens galaxy (strong lensing condition applied), with rvs/sampler as callback. The class object contains the following attribute methods:
      - `rvs`: returns random samples from the density profile slope distribution
      - `pdf`: returns the probability density function of the density profile slope distribution
      - `function`: returns the density profile slope distribution function.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **gamma** : `numpy.ndarray`
              density profile slope of the lens galaxy.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> import numpy as np
      >>> ler = OpticalDepth()
      >>> size = 10
      >>> print(ler.density_profile_slope_sl(size))



      ..
          !! processed by numpydoc !!

   .. py:property:: external_shear_sl

      
      Class object (of FunctionConditioning) for external shear of lens galaxy (strong lensing condition applied), with rvs/sampler as callback. The class object contains the following attribute methods:
      - `rvs`: returns random samples from the external shear distribution
      - `pdf`: returns the probability density function of the external shear distribution
      - `function`: returns the external shear distribution function.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **gamma1** : `numpy.ndarray`
              external shear of the lens galaxy.

          **gamma2** : `numpy.ndarray`
              external shear of the lens galaxy.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> import numpy as np
      >>> ler = OpticalDepth()
      >>> size = 10
      >>> print(ler.external_shear_sl(size))



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
      >>> print(self.optical_depth(np.array([0.1,0.2,0.3])))



      ..
          !! processed by numpydoc !!

   .. py:property:: available_lens_samplers

      
      Dictionary with list all the available priors and it's corresponding parameters. This is an immutable instance attribute.
















      ..
          !! processed by numpydoc !!

   .. py:property:: available_lens_functions

      
      Dictionary with list all the available lens functions. This is an immutable instance attribute.
















      ..
          !! processed by numpydoc !!

   .. py:attribute:: lens_type
      :value: "'epl_shear_galaxy'"

      

   .. py:attribute:: npool
      :value: '4'

      

   .. py:attribute:: z_min
      :value: '0.0'

      

   .. py:attribute:: z_max
      :value: '10.0'

      

   .. py:attribute:: cosmo

      

   .. py:attribute:: directory
      :value: "'./interpolator_json'"

      

   .. py:attribute:: comoving_distance

      

   .. py:attribute:: angular_diameter_distance

      

   .. py:attribute:: angular_diameter_distance_z1z2

      

   .. py:attribute:: differential_comoving_volume

      

   .. py:attribute:: lens_redshift_intrinsic

      

   .. py:method:: default_lens_samplers_and_functions(lens_type)

      
      Function to categorize the lens priors/samplers


      :Parameters:

          **lens_type** : `str`
              lens type
              e.g. 'epl_shear_galaxy' for elliptical power-law galaxy

      :Returns:

          **lens_param_samplers_** : `dict`
              dictionary of priors

          **lens_param_samplers_params_** : `dict`
              dictionary of priors parameters

          **lens_sampler_names_** : `dict`
              dictionary of sampler names

          **lens_functions_** : `dict`
              dictionary of lens functions













      ..
          !! processed by numpydoc !!

   .. py:method:: initialize_decision_dictionary(create_new_interpolator, lens_type)

      
      Function to initialize decision dictionary for creating interpolator


      :Parameters:

          **create_new_interpolator** : `dict` or `bool`
              dictionary to create new interpolator for velocity dispersion and optical depth.














      ..
          !! processed by numpydoc !!

   .. py:method:: lens_functions_and_sampler_categorization(lens_param_samplers, lens_param_samplers_params, lens_functions, lens_functions_params)

      
      Function to initialize velocity dispersion sampler with it's settings. The reason I am seperating this from lens_param_samplers_categorization is only a specific parameters needs special attention.


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
      >>> print(self.axis_ratio(sigma=200.))



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
      >>> print(self.axis_ratio(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: lens_redshift_strongly_lensed_numerical(size=1000, zs=None, get_attribute=False, **kwargs)

      
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
      >>> od = OpticalDepth(lens_param_samplers=dict(lens_redshift="lens_redshift_strongly_lensed_numerical"))
      >>> print(self.lens_redshift(size=10, zs=1.0))



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
      >>> print(self.axis_rotation_angle_uniform(size=10))



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
      >>> print(self.axis_ratio_uniform(size=10))



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
      >>> print(self.external_shear_normal(size=10))



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
      >>> print(self.density_profile_slope_normal(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: lens_redshift_sis_haris(size, zs, get_attribute=False, **kwargs)

      
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
      >>> lens.lens_redshift_sis_haris(zs=1.0)



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
      >>> print(self.velocity_dispersion(size=10))



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
      >>> print(self.velocity_dispersion(size=10))



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
      >>> print(self.velocity_dispersion(size=10, zl=0.5))



      ..
          !! processed by numpydoc !!

   .. py:method:: optical_depth_numerical(zs, get_attribute=False, **kwargs)


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
      >>> print(self.optical_depth_sis_haris(zs=1.0))



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


.. py:class:: ImageProperties(npool=4, z_min=0.0, z_max=10, n_min_images=2, n_max_images=4, time_window=365 * 24 * 3600 * 20, lens_model_list=['EPL_NUMBA', 'SHEAR'], cosmology=None, spin_zero=True, spin_precession=False, directory='./interpolator_json', create_new_interpolator=False)


   
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
           default: "./interpolator_json"

       **create_new_interpolator** : `dict`
           dictionary to create new interpolator pickle files
           default: dict(luminosity_distance=dict(create_new=False, resolution=1000))











   .. rubric:: Examples

   >>> from ler.image_properties import ImageProperties
   >>> image_properties = ImageProperties()
   >>> lens_parameters = dict(zs=2.0, zl=0.5, gamma1=0.0, gamma2=0.0, e1=0.0, e2=0.0, gamma=2.0, theta_E=1.0)
   >>> lens_parameters = image_properties.image_properties(lens_parameters)
   >>> print(lens_parameters.keys())

   Instance Attributes
   ----------
   ImageProperties has the following instance attributes:

   +-------------------------+----------------------+
   | Atrributes                          | Type                             |
   +=====================================+==================================+
   |:attr:`npool`                        | `int`                            |
   +-------------------------+----------------------+
   |:attr:`z_min`                        | `float`                          |
   +-------------------------+----------------------+
   |:attr:`z_max`                        | `float`                          |
   +-------------------------+----------------------+
   |:attr:`n_min_images`                 | `int`                            |
   +-------------------------+----------------------+
   |:attr:`n_max_images`                 | `int`                            |
   +-------------------------+----------------------+
   |:attr:`geocent_time_min`             | `float`                          |
   +-------------------------+----------------------+
   |:attr:`geocent_time_max`             | `float`                          |
   +-------------------------+----------------------+
   |:attr:`lens_model_list`              | `list`                           |
   +-------------------------+----------------------+
   |:attr:`cosmo`                        | `astropy.cosmology`              |
   +-------------------------+----------------------+
   |:attr:`spin_zero`                    | `bool`                           |
   +-------------------------+----------------------+
   |:attr:`spin_precession`              | `bool`                           |
   +-------------------------+----------------------+
   |:attr:`directory`                    | `str`                            |
   +-------------------------+----------------------+
   |:attr:`create_new_interpolator`      | `dict`                           |
   +-------------------------+----------------------+



   ..
       !! processed by numpydoc !!
   .. py:attribute:: npool
      :value: '4'

      

   .. py:attribute:: n_min_images
      :value: '2'

      

   .. py:attribute:: n_max_images
      :value: '4'

      

   .. py:attribute:: lens_model_list
      :value: "['EPL_NUMBA', 'SHEAR']"

      

   .. py:attribute:: spin_zero
      :value: 'True'

      

   .. py:attribute:: spin_precession
      :value: 'False'

      

   .. py:attribute:: time_window
      :value: '630720000'

      

   .. py:attribute:: cosmo

      

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

              image related=>['x_source': source position in the x-direction, 'y_source': source position in the y-direction, 'x0_image_position': image position in the x-direction, 'x1_image_position': image position in the y-direction, 'magnifications': magnifications, 'time_delays': time delays: number of images formed, 'determinant': determinants, 'trace': traces, 'iteration': to keep track of the iteration number













      ..
          !! processed by numpydoc !!

   .. py:method:: get_lensed_snrs(lensed_param, pdet_calculator, list_of_detectors=None)

      
      Function to calculate the signal to noise ratio for each image in each event.


      :Parameters:

          **list_of_detectors** : `list`
              list of detectors
              e.g. ['H1', 'L1', 'V1']

          **lensed_param** : `dict`
              dictionary containing the both already lensed source paramters and image parameters.
              e.g. lensed_param.keys() = ['mass_1', 'mass_2', 'zs', 'luminosity_distance', 'theta_jn', 'psi', 'phi', 'ra', 'dec', 'geocent_time', 'phase', 'a_1', 'a2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl', 'magnifications', 'time_delays']

      :Returns:

          **snrs** : `dict`
              signal to noise ratio for each image in each event.
              (dictionary containing 'H1', 'L1', ..., and 'snr_net', which is the network snr, for each image as an array with dimensions (number_of_lensed_events,n_max_images) )













      ..
          !! processed by numpydoc !!


.. py:function:: is_njitted(func)


.. py:class:: LensGalaxyParameterDistribution(npool=4, z_min=0.0, z_max=10.0, cosmology=None, event_type='BBH', lens_type='epl_shear_galaxy', lens_functions=None, lens_functions_params=None, lens_param_samplers=None, lens_param_samplers_params=None, directory='./interpolator_json', create_new_interpolator=False, buffer_size=1000, **kwargs)


   Bases: :py:obj:`ler.gw_source_population.CBCSourceParameterDistribution`, :py:obj:`ler.image_properties.ImageProperties`, :py:obj:`ler.lens_galaxy_population.optical_depth.OpticalDepth`

   
   Class to sample lens galaxy parameters and source parameters conditioned on the source being strongly lensed.

   This class deals with the distribution of lens galaxy parameters, such as velocity dispersion,
   axis ratio, axis rotation angle, shear, and density profile slope. It also handles the
   sampling of source parameters conditioned on the source being strongly lensed.

   :Parameters:

       **npool** : int, optional
           Number of processors to use.
           Default is 4.

       **z_min** : float, optional
           Minimum redshift.
           Default is 0.0.

       **z_max** : float, optional
           Maximum redshift.
           Default is 10.0.

       **cosmology** : astropy.cosmology, optional
           Cosmology to use.
           Default is None, which falls back to ``astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)``.

       **event_type** : str, optional
           Type of event to generate. e.g. 'BBH', 'BNS', 'NSBH'.
           Default is 'BBH'.

       **lens_type** : str, optional
           Type of lens galaxy to generate.
           Default is 'epl_shear_galaxy'.

       **lens_functions** : dict, optional
           Dictionary of lens functions.

       **lens_functions_params** : dict, optional
           Dictionary of parameters for lens functions.

       **lens_param_samplers** : dict, optional
           Dictionary of lens parameter samplers.

       **lens_param_samplers_params** : dict, optional
           Dictionary of parameters for lens parameter samplers.

       **directory** : str, optional
           Directory to store the interpolators.
           Default is './interpolator_json'.

       **create_new_interpolator** : bool, optional
           If True, creates a new interpolator.
           Default is False.

       **buffer_size** : int, optional
           Buffer size for sampling lens parameters.
           Default is 1000.

       **\*\*kwargs**
           Keyword arguments to pass to the parent classes.











   .. rubric:: Examples

   >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
   >>> lens = LensGalaxyParameterDistribution()
   >>> lensed_params = lens.sample_lens_parameters(size=1000)
   >>> print(lensed_params.keys())

   :Attributes:

       **npool** : int
           Number of processors to use.

       **z_min** : float
           Minimum redshift.

       **z_max** : float
           Maximum redshift.

       **cosmo** : astropy.cosmology
           Cosmology object.

       **event_type** : str
           Type of event to generate.

       **directory** : str
           Directory to store the interpolators.

       **create_new_interpolator** : dict
           Dictionary to check if new interpolator is created.

       **lens_param_samplers** : dict
           Dictionary of lens parameter samplers.

       **lens_param_samplers_params** : dict
           Dictionary of lens parameter sampler parameters.

       **lens_functions** : dict
           Dictionary of lens functions.

       **normalization_pdf_z_lensed** : float
           Normalization constant of the pdf p(z) for lensed events.


   ..
       !! processed by numpydoc !!
   .. py:attribute:: cbc_pop
      :value: 'None'

      
      Inherited class for sampling source parameters.
















      ..
          !! processed by numpydoc !!

      :type: :class:`~ler.gw_source_population.CBCSourceParameterDistribution`

   .. py:attribute:: z_min
      :value: 'None'

      
      Minimum redshift.
















      ..
          !! processed by numpydoc !!

      :type: float

   .. py:attribute:: z_max
      :value: 'None'

      
      Maximum redshift.
















      ..
          !! processed by numpydoc !!

      :type: float

   .. py:attribute:: m_min
      :value: 'None'

      
      Minimum mass in detector frame.
















      ..
          !! processed by numpydoc !!

      :type: float

   .. py:attribute:: m_max
      :value: 'None'

      
      Maximum mass in detector frame.
















      ..
          !! processed by numpydoc !!

      :type: float

   .. py:attribute:: normalization_pdf_z
      :value: 'None'

      
      Normalization constant of the pdf p(z).
















      ..
          !! processed by numpydoc !!

      :type: float

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

      

   .. py:attribute:: normalization_pdf_z_lensed

      
      Normalization constant of the pdf p(z) for lensed events.
















      ..
          !! processed by numpydoc !!

      :type: float

   .. py:method:: class_initialization_lens(npool, z_min, z_max, cosmology, lens_type, lens_functions, lens_functions_params, lens_param_samplers, lens_param_samplers_params, directory, create_new_interpolator, params)

      
      Initialize the LensGalaxyParameterDistribution class.


      :Parameters:

          **npool** : int
              Number of processors to use for sampling.

          **z_min** : float
              Minimum redshift of the lens galaxy.

          **z_max** : float
              Maximum redshift of the lens galaxy.

          **cosmology** : astropy.cosmology
              Cosmology object.

          **lens_type** : str
              Type of the lens galaxy.

          **lens_functions** : dict
              Dictionary with the lens related functions.

          **lens_functions_params** : dict
              Dictionary with the parameters for the lens related functions.

          **lens_param_samplers** : dict
              Dictionary with the priors for the sampler.

          **lens_param_samplers_params** : dict
              Dictionary with the parameters for the priors of the sampler.

          **directory** : str
              Directory where the interpolators are saved.

          **create_new_interpolator** : bool
              If True, creates a new interpolator.

          **params** : dict
              Additional parameters for the ``CBCSourceParameterDistribution`` and ``ImageProperties`` classes.














      ..
          !! processed by numpydoc !!

   .. py:method:: sample_lens_parameters(size=1000)

      
      Sample lens galaxy parameters along with the source parameters, conditioned on the source being strongly lensed.


      :Parameters:

          **size** : int, optional
              Number of lens parameters to sample.
              Default is 1000.

      :Returns:

          **lens_parameters** : dict
              Dictionary of sampled lens parameters and source parameters.
              Keys include ``zl``, ``zs``, ``sigma``, ``q``, ``theta_E``, ``phi``, ``e1``, ``e2``,
              ``gamma1``, ``gamma2``, ``gamma``, ``geocent_time``, ``phase``, ``psi``, ``theta_jn``,
              ``luminosity_distance``, ``mass_1_source``, ``mass_2_source``, ``ra``, ``dec``.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> od = LensGalaxyParameterDistribution(lens_param_samplers=dict(velocity_dispersion="velocity_dispersion_ewoud"))
      >>> print(od.sample_lens_parameters(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: sample_all_routine_epl_shear_sl(size=1000)

      
      Sample galaxy lens parameters. EPL shear cross section is used for rejection sampling.


      :Parameters:

          **size** : int, optional
              Number of lens parameters to sample.
              Default is 1000.

      :Returns:

          **lens_parameters** : dict
              Dictionary of lens parameters and source parameters (lens conditions applied).
              Keys include ``zl``, ``zs``, ``sigma``, ``q``, ``theta_E``, ``phi``, ``e1``, ``e2``,
              ``gamma1``, ``gamma2``, ``gamma``.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.sample_all_routine_epl_shear_sl(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: strongly_lensed_source_redshifts(size=1000)

      
      Sample source redshifts, conditioned on the source being strongly lensed.


      :Parameters:

          **size** : int, optional
              Number of lens parameters to sample.
              Default is 1000.

      :Returns:

          **redshifts** : numpy.ndarray
              Source redshifts conditioned on the source being strongly lensed.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.strongly_lensed_source_redshifts(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: sample_all_routine_epl_shear_intrinsic(size=1000)

      
      Sample galaxy lens parameters. EPL shear cross section is used for rejection sampling.


      :Parameters:

          **size** : int, optional
              Number of lens parameters to sample.
              Default is 1000.

      :Returns:

          **lens_parameters** : dict
              Dictionary of lens parameters and source parameters (lens conditions applied).
              Keys include ``zl``, ``zs``, ``sigma``, ``q``, ``theta_E``, ``phi``, ``e1``, ``e2``,
              ``gamma1``, ``gamma2``, ``gamma``.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
      >>> lens = LensGalaxyParameterDistribution()
      >>> lens.sample_all_routine_epl_shear_intrinsic(size=1000)



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

.. py:function:: axis_ratio_rayleigh_pdf(q, sigma, q_min=0.2, q_max=1.0)


.. py:function:: cross_section_mp(params)


.. py:class:: OpticalDepth(npool=4, z_min=0.0, z_max=10.0, cosmology=None, lens_type='epl_shear_galaxy', lens_functions=None, lens_functions_params=None, lens_param_samplers=None, lens_param_samplers_params=None, directory='./interpolator_json', create_new_interpolator=False, verbose=False)


   
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
           Directory where the interpolators are saved (default is './interpolator_json').
           If True, creates a new interpolator (default is False).

       **verbose** : bool, optional
           If True, prints additional information during initialization (default is False).





   :Raises:

       ValueError
           If `lens_type` is not in ['sie_galaxy', 'epl_shear_galaxy', 'sis_galaxy'].









   ..
       !! processed by numpydoc !!
   .. py:property:: velocity_dispersion

      
      Class object (of FunctionConditioning) for velocity dispersion of lens galaxy, with rvs/sampler as callback. Lens redshift `zl` is required only if velocity dispersion is redshift dependent. The class object contains the following attribute methods:
      - `rvs`: returns random samples from the velocity dispersion distribution
      - `pdf`: returns the probability density function of the velocity dispersion distribution
      - `function`: returns the velocity dispersion distribution function which represents the number density of lens galaxies as a function of velocity dispersion


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **zl** : `float`
              redshift of the lens galaxy. Should be of shape (size,)

      :Returns:

          **sigma** : `numpy.ndarray`
              velocity dispersion of the lens galaxy










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> import numpy as np
      >>> ler = OpticalDepth()
      >>> size = 10
      >>> zl = np.ones(size)*1.5
      >>> print(ler.velocity_dispersion(size, zl))



      ..
          !! processed by numpydoc !!

   .. py:property:: axis_ratio

      
      Class object (of FunctionConditioning) for axis ratio of lens galaxy, with rvs/sampler as callback. Velocity dispersion `sigma` is required only if axis ratio is velocity dispersion dependent. The class object contains the following attribute methods:
      - `rvs`: returns random samples from the axis ratio distribution
      - `pdf`: returns the probability density function of the axis ratio distribution
      - `function`: returns the un-normalized axis ratio distribution function.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **sigma** : `numpy.ndarray`
              velocity dispersion of the lens galaxy. Should be of shape (size,)

      :Returns:

          **q** : `numpy.ndarray`
              axis ratio of the lens galaxy.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> import numpy as np
      >>> ler = OpticalDepth()
      >>> size = 10
      >>> sigma = np.ones(size)*150.0
      >>> print(ler.axis_ratio(size, sigma))



      ..
          !! processed by numpydoc !!

   .. py:property:: axis_rotation_angle

      
      Class object (of FunctionConditioning) for axis rotation angle of lens galaxy, with rvs/sampler as callback. The class object contains the following attribute methods:
      - `rvs`: returns random samples from the axis rotation angle distribution
      - `pdf`: returns the probability density function of the axis rotation angle distribution
      - `function`: returns the axis rotation angle distribution function.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **theta** : `numpy.ndarray`
              axis rotation angle of the lens galaxy.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> import numpy as np
      >>> ler = OpticalDepth()
      >>> size = 10
      >>> print(ler.axis_rotation_angle(size))



      ..
          !! processed by numpydoc !!

   .. py:property:: density_profile_slope

      
      Class object (of FunctionConditioning) for density profile slope of lens galaxy, with rvs/sampler as callback. The class object contains the following attribute methods:
      - `rvs`: returns random samples from the density profile slope distribution
      - `pdf`: returns the probability density function of the density profile slope distribution
      - `function`: returns the density profile slope distribution function.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **gamma** : `numpy.ndarray`
              density profile slope of the lens galaxy.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> import numpy as np
      >>> ler = OpticalDepth()
      >>> size = 10
      >>> print(ler.density_profile_slope(size))



      ..
          !! processed by numpydoc !!

   .. py:property:: external_shear

      
      Class object (of FunctionConditioning) for external shear of lens galaxy, with rvs/sampler as callback. The class object contains the following attribute methods:
      - `rvs`: returns random samples from the external shear distribution
      - `pdf`: returns the probability density function of the external shear distribution
      - `function`: returns the external shear distribution function.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **gamma1** : `numpy.ndarray`
              external shear of the lens galaxy.

          **gamma2** : `numpy.ndarray`
              external shear of the lens galaxy.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> import numpy as np
      >>> ler = OpticalDepth()
      >>> size = 10
      >>> print(ler.external_shear(size))



      ..
          !! processed by numpydoc !!

   .. py:property:: cross_section

      
      Lensing cross section for individual lensing events. It can be of the following lens type and corresponding input parameters:
      - `epl_shear_galaxy`: `zs`, `zl`, `sigma`, `q`, `phi`, `gamma`, `gamma1`, `gamma2`
      - `sie_galaxy`: `zs`, `zl`, `sigma`, `q`
      - `sis_galaxy`: `zs`, `zl`, `sigma`


      :Parameters:

          **zs** : `numpy.ndarray`
              Redshift of the source

          **zl** : `numpy.ndarray`
              Redshift of the lens

          **sigma** : `numpy.ndarray`
              Angular size of the lens

          **q** : `numpy.ndarray`
              Axis ratio of the lens

          **phi** : `numpy.ndarray`
              Position angle of the lens

          **gamma** : `numpy.ndarray`
              density profile slope of the lens

          **gamma1** : `numpy.ndarray`
              Shear of the lens (x-direction)

          **gamma2** : `numpy.ndarray`
              Shear of the lens (y-direction)

      :Returns:

          **cross_section** : `numpy.ndarray`
              Lensing cross section for individual lensing events













      ..
          !! processed by numpydoc !!

   .. py:property:: lens_redshift

      
      Class object (of FunctionConditioning) for lens redshift, with rvs/sampler as callback. The class object contains the following attribute methods:
      - `rvs`: returns random samples from the lens redshift distribution
      - `pdf`: returns the probability density function of the lens redshift distribution
      - `function`: returns the lens redshift distribution function which represents effective lensing cross-section for lenses at redshift zl,


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

          **zs** : `numpy.ndarray`
              redshift of the lens galaxy. Should be of shape (size,)

      :Returns:

          **zl** : `numpy.ndarray`
              redshift of the lens galaxy.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> import numpy as np
      >>> ler = OpticalDepth()
      >>> size = 10
      >>> zs = np.ones(size)*1.5
      >>> print(ler.lens_redshift(size, zs))



      ..
          !! processed by numpydoc !!

   .. py:property:: density_profile_slope_sl

      
      Class object (of FunctionConditioning) for density profile slope of lens galaxy (strong lensing condition applied), with rvs/sampler as callback. The class object contains the following attribute methods:
      - `rvs`: returns random samples from the density profile slope distribution
      - `pdf`: returns the probability density function of the density profile slope distribution
      - `function`: returns the density profile slope distribution function.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **gamma** : `numpy.ndarray`
              density profile slope of the lens galaxy.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> import numpy as np
      >>> ler = OpticalDepth()
      >>> size = 10
      >>> print(ler.density_profile_slope_sl(size))



      ..
          !! processed by numpydoc !!

   .. py:property:: external_shear_sl

      
      Class object (of FunctionConditioning) for external shear of lens galaxy (strong lensing condition applied), with rvs/sampler as callback. The class object contains the following attribute methods:
      - `rvs`: returns random samples from the external shear distribution
      - `pdf`: returns the probability density function of the external shear distribution
      - `function`: returns the external shear distribution function.


      :Parameters:

          **size** : `int`
              number of lens parameters to sample

      :Returns:

          **gamma1** : `numpy.ndarray`
              external shear of the lens galaxy.

          **gamma2** : `numpy.ndarray`
              external shear of the lens galaxy.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> import numpy as np
      >>> ler = OpticalDepth()
      >>> size = 10
      >>> print(ler.external_shear_sl(size))



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
      >>> print(self.optical_depth(np.array([0.1,0.2,0.3])))



      ..
          !! processed by numpydoc !!

   .. py:property:: available_lens_samplers

      
      Dictionary with list all the available priors and it's corresponding parameters. This is an immutable instance attribute.
















      ..
          !! processed by numpydoc !!

   .. py:property:: available_lens_functions

      
      Dictionary with list all the available lens functions. This is an immutable instance attribute.
















      ..
          !! processed by numpydoc !!

   .. py:attribute:: lens_type
      :value: "'epl_shear_galaxy'"

      

   .. py:attribute:: npool
      :value: '4'

      

   .. py:attribute:: z_min
      :value: '0.0'

      

   .. py:attribute:: z_max
      :value: '10.0'

      

   .. py:attribute:: cosmo

      

   .. py:attribute:: directory
      :value: "'./interpolator_json'"

      

   .. py:attribute:: comoving_distance

      

   .. py:attribute:: angular_diameter_distance

      

   .. py:attribute:: angular_diameter_distance_z1z2

      

   .. py:attribute:: differential_comoving_volume

      

   .. py:attribute:: lens_redshift_intrinsic

      

   .. py:method:: default_lens_samplers_and_functions(lens_type)

      
      Function to categorize the lens priors/samplers


      :Parameters:

          **lens_type** : `str`
              lens type
              e.g. 'epl_shear_galaxy' for elliptical power-law galaxy

      :Returns:

          **lens_param_samplers_** : `dict`
              dictionary of priors

          **lens_param_samplers_params_** : `dict`
              dictionary of priors parameters

          **lens_sampler_names_** : `dict`
              dictionary of sampler names

          **lens_functions_** : `dict`
              dictionary of lens functions













      ..
          !! processed by numpydoc !!

   .. py:method:: initialize_decision_dictionary(create_new_interpolator, lens_type)

      
      Function to initialize decision dictionary for creating interpolator


      :Parameters:

          **create_new_interpolator** : `dict` or `bool`
              dictionary to create new interpolator for velocity dispersion and optical depth.














      ..
          !! processed by numpydoc !!

   .. py:method:: lens_functions_and_sampler_categorization(lens_param_samplers, lens_param_samplers_params, lens_functions, lens_functions_params)

      
      Function to initialize velocity dispersion sampler with it's settings. The reason I am seperating this from lens_param_samplers_categorization is only a specific parameters needs special attention.


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
      >>> print(self.axis_ratio(sigma=200.))



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
      >>> print(self.axis_ratio(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: lens_redshift_strongly_lensed_numerical(size=1000, zs=None, get_attribute=False, **kwargs)

      
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
      >>> od = OpticalDepth(lens_param_samplers=dict(lens_redshift="lens_redshift_strongly_lensed_numerical"))
      >>> print(self.lens_redshift(size=10, zs=1.0))



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
      >>> print(self.axis_rotation_angle_uniform(size=10))



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
      >>> print(self.axis_ratio_uniform(size=10))



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
      >>> print(self.external_shear_normal(size=10))



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
      >>> print(self.density_profile_slope_normal(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: lens_redshift_sis_haris(size, zs, get_attribute=False, **kwargs)

      
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
      >>> lens.lens_redshift_sis_haris(zs=1.0)



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
      >>> print(self.velocity_dispersion(size=10))



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
      >>> print(self.velocity_dispersion(size=10))



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
      >>> print(self.velocity_dispersion(size=10, zl=0.5))



      ..
          !! processed by numpydoc !!

   .. py:method:: optical_depth_numerical(zs, get_attribute=False, **kwargs)


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
      >>> print(self.optical_depth_sis_haris(zs=1.0))



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

.. py:function:: make_cross_section_reinit(e1_grid, e2_grid, gamma_grid, gamma1_grid, gamma2_grid, cs_spline_coeff_grid, Da_instance, csunit_to_cs_slope=0.31830988618379075, csunit_to_cs_intercept=-3.2311742677852644e-27)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: lens_redshift_strongly_lensed_njit(zs_array, zl_scaled, sigma_min, sigma_max, q_rvs, phi_rvs, gamma_rvs, shear_rvs, number_density, cross_section, dVcdz_function, integration_size)


.. py:function:: lens_redshift_strongly_lensed_mp(params)


.. py:function:: cross_section_unit_mp(params)


.. py:function:: cross_section_mp(params)


.. py:function:: cross_section(theta_E, e1, e2, gamma, gamma1, gamma2)


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

.. py:function:: sample_sigma_zl(pdf, sigma_min, sigma_max, zl_min, zl_max, zs, chunk_size=10000)


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

.. py:data:: C_LIGHT
   :value: '299792.458'

   

.. py:function:: cross_section(zs, zl, sigma, q, phi, gamma, gamma1, gamma2, e1_grid, e2_grid, gamma_grid, gamma1_grid, gamma2_grid, cs_spline_coeff_grid, Da_instance, csunit_to_cs_slope, csunit_to_cs_intercept)

   
   theta_E, gamma, gamma1, gamma2, q, phi: 1D arrays (same length)
   grids: 1D
   cs_spline_coeff: spline-filtered coefficients (same shape as cs_spline_coeff)
   returns: 1D array
















   ..
       !! processed by numpydoc !!

.. py:function:: make_cross_section_reinit(e1_grid, e2_grid, gamma_grid, gamma1_grid, gamma2_grid, cs_spline_coeff_grid, Da_instance, csunit_to_cs_slope=0.31830988618379075, csunit_to_cs_intercept=-3.2311742677852644e-27)

   
















   ..
       !! processed by numpydoc !!

.. py:data:: C_LIGHT
   :value: '299792.458'

   

.. py:function:: rejection_sampler(zs, zl, sigma_max, sigma_rvs, q_rvs, phi_rvs, gamma_rvs, shear_rvs, cross_section, saftey_factor=1.2)


.. py:function:: create_rejection_sampler(sigma_max, sigma_rvs, q_rvs, phi_rvs, gamma_rvs, shear_rvs, cross_section, saftey_factor=1.2, use_njit_sampler=True)


.. py:function:: sigma_proposal_uniform(n, sigma_min, sigma_max)


.. py:function:: weighted_choice_1d(weights)

   
   Draw an index with probability proportional to 'weights' (assumed >=0).
   Numba-safe replacement for np.random.choice(n, p=weights).
















   ..
       !! processed by numpydoc !!

.. py:function:: importance_sampler(zs, zl, sigma_min, sigma_max, q_rvs, phi_rvs, gamma_rvs, shear_rvs, sigma_pdf, cross_section, n_prop)


.. py:function:: create_importance_sampler(sigma_min, sigma_max, q_rvs, phi_rvs, gamma_rvs, shear_rvs, sigma_pdf, cross_section, n_prop, use_njit_sampler=True)


