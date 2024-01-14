:orphan:

:py:mod:`ler.lens_galaxy_population.optical_depth`
==================================================

.. py:module:: ler.lens_galaxy_population.optical_depth


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.optical_depth.OpticalDepth




.. py:class:: OpticalDepth(npool=4, z_min=0.001, z_max=10.0, optical_depth_function=None, sampler_priors=None, sampler_priors_params=None, cosmology=None, directory='./interpolator_pickle', create_new_interpolator=False)


   
   Class to calculate the optical depth, velocity dispersion and axis-ratio of a lens galaxy population.


   :Parameters:

       **npool** : `int`
           number of processors to use for multiprocessing

       **z_min** : `float`
           minimum redshift of the lens galaxy population

       **z_max** : `float`
           maximum redshift of the lens galaxy population

       **optical_depth_function** : `str` or `callable`
           Function or function name to calculate optical depth.
           Check for default/available optical depth functions by running,
           >>> from ler.lens_galaxy_population import OpticalDepth
           >>> print(OpticalDepth().available_optical_depth_list_and_its_params)

       **sampler_priors, sampler_priors_params** : `dict`, `dict`
           dictionary of sampler functions and it's parameters to sample velocity dispersion and axis-ratio.
           Check for default/available sampler priors and corresponding input parameters by running,
           >>> from ler.lens_galaxy_population import OpticalDepth
           >>> print(OpticalDepth().available_velocity_dispersion_list_and_its_params)
           >>> print(OpticalDepth().available_axis_ratio_list_and_its_params)

       **cosmology** : `astropy.cosmology`
           Cosmology to use
           default: None/astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)

       **directory** : `str`
           directory to store interpolator pickle files
           default: "./interpolator_pickle"

       **create_new_interpolator** : `dict`
           dictionary to create new interpolator for velocity dispersion and optical depth.











   .. rubric:: Examples

   >>> from ler.lens_galaxy_population import OpticalDepth
   >>> od = OpticalDepth()
   >>> print(od.strong_lensing_optical_depth(0.5))



   ..
       !! processed by numpydoc !!
   .. py:property:: strong_lensing_optical_depth

      
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
      >>> print(od.strong_lensing_optical_depth(np.array([0.1,0.2,0.3])))



      ..
          !! processed by numpydoc !!

   .. py:property:: sample_velocity_dispersion

      
      Function to sample velocity dispersion. `zl` is required only if velocity dispersion sampler is redshift dependent.


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
      >>> print(od.sample_velocity_dispersion(size=10))



      ..
          !! processed by numpydoc !!

   .. py:property:: sample_axis_ratio

      
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
      >>> print(od.sample_axis_ratio(sigma=200.))



      ..
          !! processed by numpydoc !!

   .. py:property:: available_velocity_dispersion_list_and_its_params

      
      Function to list all available velocity dispersion sampler and its parameters.
















      ..
          !! processed by numpydoc !!

   .. py:property:: available_axis_ratio_list_and_its_params

      
      Function to list all available axis ratio sampler.
















      ..
          !! processed by numpydoc !!

   .. py:property:: available_optical_depth_list_and_its_params

      
      Function to list all available optical depth sampler.
















      ..
          !! processed by numpydoc !!

   .. py:method:: initialize_velocity_dispersion_sampler(vd_name)

      
      Function to initialize velocity dispersion sampler


      :Parameters:

          **vd_name** : `str`
              name of velocity dispersion sampler














      ..
          !! processed by numpydoc !!

   .. py:method:: initialize_optical_depth_function(tau_name, vd_name)

      
      Function to initialize optical depth function.


      :Parameters:

          **tau_name** : `str`
              name of optical depth function

          **vd_name** : `str`
              name of velocity dispersion sampler














      ..
          !! processed by numpydoc !!

   .. py:method:: axis_ratio_rayleigh(sigma, q_min=0.2, q_max=1.0, get_attribute=False, param=None, **kwargs)

      
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
      >>> od = OpticalDepth(sampler_priors=dict(axis_ratio="axis_ratio_rayleigh"))
      >>> print(od.sample_axis_ratio(sigma=200.))



      ..
          !! processed by numpydoc !!

   .. py:method:: axis_ratio_padilla_strauss(size=1000, q_min=0.2, q_max=1.0, get_attribute=False, param=None, **kwargs)

      
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
      >>> od = OpticalDepth(sampler_priors=dict(axis_ratio="axis_ratio_padilla_strauss"))
      >>> print(od.sample_axis_ratio(size=10))



      ..
          !! processed by numpydoc !!

   .. py:method:: velocity_dispersion_gengamma(size, a=2.32 / 2.67, c=2.67, get_attribute=False, param=None, **kwargs)

      
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
      >>> od = OpticalDepth(sampler_priors=dict(velocity_dispersion="velocity_dispersion_gengamma"), sampler_priors_params=dict(velocity_dispersion=dict(a=2.32 / 2.67, c=2.67)))
      >>> print(od.sample_velocity_dispersion(size=10))



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
      >>> od = OpticalDepth(sampler_priors=dict(velocity_dispersion="velocity_dispersion_bernardi"))
      >>> print(od.sample_velocity_dispersion(size=10))



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
      >>> od = OpticalDepth(sampler_priors=dict(velocity_dispersion="velocity_dispersion_ewoud"))
      >>> print(od.sample_velocity_dispersion(size=10, zl=0.5))



      ..
          !! processed by numpydoc !!

   .. py:method:: cross_section_SIS(sigma, zl, zs)

      
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
      >>> print(od.cross_section_SIS(sigma=200., zl=0.5, zs=1.0))



      ..
          !! processed by numpydoc !!

   .. py:method:: tau_zl_zs(zl, zs)

      
      Function to compute the optical depth for a given lens redshift and source redshift


      :Parameters:

          **zl** : `float`
              redshift of the lens galaxy

          **zs** : `float`
              redshift of the source galaxy

      :Returns:

          **tau** : `float`
              optical depth










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> print(od.tau_zl_zs(zl=0.5, zs=1.0))



      ..
          !! processed by numpydoc !!

   .. py:method:: optical_depth_calculator(zs)

      
      Function to compute the optical depth without multiprocessing. This is the integrated version of tau_zl_zs from z=0 to z=zs.


      :Parameters:

          **zs** : `float`
              source redshifts

      :Returns:

          **tau** : `float`
              optical depth










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> print(od.optical_depth_calculator(zs=1.0))



      ..
          !! processed by numpydoc !!

   .. py:method:: optical_depth_multiprocessing(zs)

      
      Function to compute the optical depth with multiprocessing. This is the integrated version of optical depth from z=0 to z=zs.


      :Parameters:

          **zs** : `float`
              source redshifts

      :Returns:

          **tau** : `float`
              optical depth










      .. rubric:: Examples

      >>> from ler.lens_galaxy_population import OpticalDepth
      >>> od = OpticalDepth()
      >>> print(od.optical_depth_multiprocessing(zs=1.0))



      ..
          !! processed by numpydoc !!

   .. py:method:: optical_depth_SIS_haris(zs)

      
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
      >>> print(od.optical_depth_SIS_haris(zs=1.0))



      ..
          !! processed by numpydoc !!

   .. py:method:: create_lookup_table_fuction(z_max)

      
      Functions to create lookup tables
      1. Redshift to co-moving distance.
      2. Co-moving distance to redshift.
      3. Redshift to angular diameter distance.


      :Parameters:

          **z_max** : `float`
              maximum redshift of the lens galaxy population














      ..
          !! processed by numpydoc !!


