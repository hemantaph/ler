:py:mod:`ler.lens_galaxy_population.optical_depth`
==================================================

.. py:module:: ler.lens_galaxy_population.optical_depth


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.optical_depth.OpticalDepth




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


