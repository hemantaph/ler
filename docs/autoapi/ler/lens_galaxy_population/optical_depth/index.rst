:py:mod:`ler.lens_galaxy_population.optical_depth`
==================================================

.. py:module:: ler.lens_galaxy_population.optical_depth

.. autoapi-nested-parse::

   Module for optical depth and lens parameter distribution calculations.

   This module provides the ``OpticalDepth`` class for computing strong lensing
   optical depth, cross-section, sampling velocity dispersion, axis ratio, and other lens galaxy
   population parameters. It supports multiple lens models including SIS, SIE,
   and EPL + external shear.

   Key features:

   - Optical depth computation for strong gravitational lensing

   - Velocity dispersion sampling with multiple models

   - Lens redshift distribution sampling

   - Cross-section calculations for various lens models


   Copyright (C) 2024 Hemantakumar Phurailatpam. Distributed under MIT License.

   ..
       !! processed by numpydoc !!


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.optical_depth.OpticalDepth




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
              Additional parameters.

              +-------------------+-------------------+---------------------------------------+
              | Parameter         | Unit              | Description                           |
              +===================+===================+=======================================+
              | sigma_min         | km/s              | minimum velocity dispersion           |
              +-------------------+-------------------+---------------------------------------+
              | sigma_max         | km/s              | maximum velocity dispersion           |
              +-------------------+-------------------+---------------------------------------+
              | alpha             | dimensionless     | Power-law index governing the slope   |
              |                   |                   | of the distribution at low velocities |
              +-------------------+-------------------+---------------------------------------+
              | beta              | dimensionless     | Exponential parameter determining the |
              |                   |                   | sharpness of the high-velocity cutoff |
              +-------------------+-------------------+---------------------------------------+
              | phistar           | h^3 Mpc^-3        | Normalization constant representing   |
              |                   |                   | the comoving number density of galaxy |
              +-------------------+-------------------+---------------------------------------+
              | sigmastar         | km/s              | Characteristic velocity scale marking |
              |                   |                   | the "knee" transition in the VDF      |
              +-------------------+-------------------+---------------------------------------+

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

              +-------------------+-------------------+---------------------------------------+
              | Parameter         | Unit              | Description                           |
              +===================+===================+=======================================+
              | sigma_min         | km/s              | minimum velocity dispersion           |
              +-------------------+-------------------+---------------------------------------+
              | sigma_max         | km/s              | maximum velocity dispersion           |
              +-------------------+-------------------+---------------------------------------+
              | alpha             | dimensionless     | Power-law index governing the slope   |
              |                   |                   | of the distribution at low velocities |
              +-------------------+-------------------+---------------------------------------+
              | beta              | dimensionless     | Exponential parameter determining the |
              |                   |                   | sharpness of the high-velocity cutoff |
              +-------------------+-------------------+---------------------------------------+
              | phistar           | h^3 Mpc^-3        | Normalization constant representing   |
              |                   |                   | the comoving number density of galaxy |
              +-------------------+-------------------+---------------------------------------+
              | sigmastar         | km/s              | Characteristic velocity scale marking |
              |                   |                   | the "knee" transition in the VDF      |
              +-------------------+-------------------+---------------------------------------+

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

              +-------------------+-------------------+---------------------------------------+
              | Parameter         | Unit              | Description                           |
              +===================+===================+=======================================+
              | sigma_min         | km/s              | minimum velocity dispersion           |
              +-------------------+-------------------+---------------------------------------+
              | sigma_max         | km/s              | maximum velocity dispersion           |
              +-------------------+-------------------+---------------------------------------+
              | alpha             | dimensionless     | Power-law index governing the slope   |
              |                   |                   | of the distribution at low velocities |
              +-------------------+-------------------+---------------------------------------+
              | beta              | dimensionless     | Exponential parameter determining the |
              |                   |                   | sharpness of the high-velocity cutoff |
              +-------------------+-------------------+---------------------------------------+
              | phistar           | h^3 Mpc^-3        | Normalization constant representing   |
              |                   |                   | the comoving number density of galaxy |
              +-------------------+-------------------+---------------------------------------+
              | sigmastar         | km/s              | Characteristic velocity scale marking |
              |                   |                   | the "knee" transition in the VDF      |
              +-------------------+-------------------+---------------------------------------+

      :Returns:

          **sigma** : ``numpy.ndarray`` or ``FunctionConditioning``
              Velocity dispersion samples (km/s) or sampler object.










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> sigma = od.velocity_dispersion(size=100, zl=np.ones(100)*0.5)



      ..
          !! processed by numpydoc !!

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

   .. py:method:: cross_section_epl_shear_numerical(zs, zl, sigma, q, phi, gamma, gamma1, gamma2)

      
      Function to compute the strong lensing cross-section numerically for EPL + external shear lenses.


      :Parameters:

          **zs** : `numpy.ndarray`
              redshift of the source galaxies

          **zl** : `numpy.ndarray`
              redshift of the lens galaxies

          **sigma** : `numpy.ndarray`
              velocity dispersion of the lens galaxies

          **q** : `numpy.ndarray`
              axis ratios of the lens galaxies

          **phi** : `numpy.ndarray`
              axis rotation angles of the lens galaxies in radians

          **gamma** : `numpy.ndarray`
              external shear magnitudes of the lens galaxies

          **gamma1** : `numpy.ndarray`
              external shear component 1 of the lens galaxies

          **gamma2** : `numpy.ndarray`
              external shear component 2 of the lens galaxies

      :Returns:

          **cross_section** : `numpy.ndarray`
              strong lensing cross-section of the lens galaxies in square radians













      ..
          !! processed by numpydoc !!

   .. py:method:: cross_section_epl_shear_numerical_mp(theta_E, gamma, gamma1, gamma2, q=None, phi=None, e1=None, e2=None, verbose=False, **kwargs)

      
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

      
      Create a parameter grid for lens galaxies.


      :Parameters:

          **size_list** : list
              List of sizes for each parameter grid.

      :Returns:

          **zl** : numpy.ndarray
              Lens redshifts.

          **sigma** : numpy.ndarray
              Velocity dispersions.

          **q** : numpy.ndarray
              Axis ratios.













      ..
          !! processed by numpydoc !!

   .. py:method:: cross_section_epl_shear_interpolation_init(file_path, size_list)


   .. py:method:: cross_section_epl_shear_interpolation(zs, zl, sigma, q, phi, gamma, gamma1, gamma2, get_attribute=False, size_list=[25, 25, 45, 15, 15], **kwargs)

      
      Function to compute the cross-section correction factor
















      ..
          !! processed by numpydoc !!


