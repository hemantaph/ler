:py:mod:`ler.gw_source_population`
==================================

.. py:module:: ler.gw_source_population


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   cbc_source_parameter_distribution/index.rst
   cbc_source_redshift_distribution/index.rst
   prior_functions/index.rst
   sfr_with_time_delay/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   ler.gw_source_population.FunctionConditioning
   ler.gw_source_population.CBCSourceRedshiftDistribution
   ler.gw_source_population.FunctionConditioning
   ler.gw_source_population.CBCSourceRedshiftDistribution
   ler.gw_source_population.CBCSourceParameterDistribution



Functions
~~~~~~~~~

.. autoapisummary::

   ler.gw_source_population.interpolator_json_path
   ler.gw_source_population.luminosity_distance
   ler.gw_source_population.differential_comoving_volume
   ler.gw_source_population.inverse_transform_sampler
   ler.gw_source_population.sample_from_powerlaw_distribution
   ler.gw_source_population.merger_rate_density_bbh_oguri2018_function
   ler.gw_source_population.merger_rate_density_bbh_popIII_ken2022_function
   ler.gw_source_population.merger_rate_density_madau_dickinson2014_function
   ler.gw_source_population.merger_rate_density_madau_dickinson_belczynski_ng_function
   ler.gw_source_population.merger_rate_density_bbh_primordial_ken2022_function
   ler.gw_source_population.sfr_madau_fragos2017_with_bbh_td
   ler.gw_source_population.sfr_madau_dickinson2014_with_bbh_td
   ler.gw_source_population.sfr_madau_fragos2017_with_bns_td
   ler.gw_source_population.sfr_madau_dickinson2014_with_bns_td
   ler.gw_source_population.sfr_madau_fragos2017
   ler.gw_source_population.sfr_madau_dickinson2014
   ler.gw_source_population.binary_masses_BBH_popIII_lognormal_rvs
   ler.gw_source_population.binary_masses_BBH_primordial_lognormal_rvs
   ler.gw_source_population.binary_masses_BNS_bimodal_rvs
   ler.gw_source_population.binary_masses_NSBH_broken_powerlaw_rvs
   ler.gw_source_population.binary_masses_BBH_powerlaw_gaussian_rvs
   ler.gw_source_population.available_prior_list
   ler.gw_source_population.sfr_madau_fragos2017
   ler.gw_source_population.sfr_with_time_delay_function



Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.gw_source_population.chunk_size


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

.. py:function:: luminosity_distance(z=None, z_min=0.001, z_max=10.0, cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0), directory='./interpolator_json', create_new=False, resolution=500, get_attribute=True)

   
   Function to create a lookup table for the luminosity distance wrt redshift.


   :Parameters:

       **z** : `numpy.ndarray` or `float`
           Source redshifts

       **z_min** : `float`
           Minimum redshift of the source population

       **z_max** : `float`
           Maximum redshift of the source population












   :Attributes:

       **z_to_luminosity_distance** : `ler.utils.FunctionConditioning`
           Object of FunctionConditioning class containing the luminosity distance wrt redshift


   ..
       !! processed by numpydoc !!

.. py:function:: differential_comoving_volume(z=None, z_min=0.001, z_max=10.0, cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0), directory='./interpolator_json', create_new=False, resolution=500, get_attribute=True)


.. py:class:: CBCSourceRedshiftDistribution(npool=4, z_min=0.001, z_max=10.0, event_type='BBH', merger_rate_density=None, merger_rate_density_param=None, cosmology=None, directory='./interpolator_json', create_new_interpolator=False)


   Bases: :py:obj:`object`

   
   Class for generating compact binary coalescence source redshift distributions.

   This class generates source redshift distributions for compact binary
   coalescence events (BBH, BNS, NSBH) using various astrophysical merger rate
   density models. It provides interpolated functions for efficient sampling of
   source redshifts weighted by the merger rate density in the detector frame.

   Key Features:

   - Multiple merger rate density models (PopI/II, PopIII, Primordial)

   - Configurable cosmology for distance calculations

   - Cached interpolators for computational efficiency

   - Support for user-defined merger rate density functions

   :Parameters:

       **npool** : ``int``
           Number of processors to use for multiprocessing.

           default: 4

       **z_min** : ``float``
           Minimum redshift of the source population.

           default: 0.001

       **z_max** : ``float``
           Maximum redshift of the source population.

           default: 10.0

       **event_type** : ``str``
           Type of compact binary event.

           Options:

           - 'BBH': Binary black hole

           - 'BNS': Binary neutron star

           - 'NSBH': Neutron star-black hole

           default: 'BBH'

       **merger_rate_density** : ``str`` or ``callable`` or ``None``
           Merger rate density model to use.

           Options:

           - 'merger_rate_density_bbh_oguri2018': PopI/II BBH (Oguri 2018)

           - 'merger_rate_density_bbh_popIII_ken2022': PopIII BBH (Ng 2022)

           - 'merger_rate_density_bbh_primordial_ken2022': Primordial BBH (Ng 2022)

           - callable: User-defined function f(z) -> rate density

           default: None (uses 'merger_rate_density_bbh_oguri2018')

       **merger_rate_density_param** : ``dict`` or ``None``
           Parameters for the merger rate density function.

           default: None (uses dict(R0=19 * 1e-9, b2=1.6, b3=2.1, b4=30))

       **cosmology** : ``astropy.cosmology`` or ``None``
           Cosmology for distance calculations.

           default: None (uses LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0))

       **directory** : ``str``
           Directory to store interpolator JSON files.

           default: './interpolator_json'

       **create_new_interpolator** : ``dict`` or ``bool``
           Control interpolator creation.

           If ``bool``: Apply to all interpolators.

           If ``dict``: Per-quantity settings with keys 'create_new' and 'resolution'.

           default: False











   .. rubric:: Examples

   Basic usage:

   >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
   >>> cbc = CBCSourceRedshiftDistribution(z_min=0.001, z_max=10)
   >>> zs_samples = cbc.source_redshift(size=1000)
   >>> rate = cbc.merger_rate_density(zs=0.5)

   Instance Methods
   ----------
   CBCSourceRedshiftDistribution has the following methods:

   +-----------------------------------------------------+----------------------------------------------------+
   | Method                                              | Description                                        |
   +=====================================================+====================================================+
   | :meth:`~merger_rate_density_detector_frame`         | Compute merger rate density in detector frame      |
   +-----------------------------------------------------+----------------------------------------------------+
   | :meth:`~merger_rate_density_bbh_oguri2018`          | PopI/II merger rate density (Oguri 2018)           |
   +-----------------------------------------------------+----------------------------------------------------+
   | :meth:`~sfr_madau_dickinson2014`                    | Star formation rate (Madau & Dickinson 2014)       |
   +-----------------------------------------------------+----------------------------------------------------+
   | :meth:`~sfr_with_time_delay`                        | SFR with time delay convolution                    |
   +-----------------------------------------------------+----------------------------------------------------+
   | :meth:`~merger_rate_density_bbh_popIII_ken2022`     | PopIII merger rate density (Ng 2022)               |
   +-----------------------------------------------------+----------------------------------------------------+
   | :meth:`~merger_rate_density_bbh_primordial_ken2022` | Primordial BBH merger rate density (Ng 2022)       |
   +-----------------------------------------------------+----------------------------------------------------+

   Instance Attributes
   ----------
   CBCSourceRedshiftDistribution has the following attributes:

   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | Attribute                                      | Type                      | Unit  | Description                                  |
   +================================================+===========================+=======+==============================================+
   | :attr:`~z_min`                                 | ``float``                 |       | Minimum source redshift                      |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~z_max`                                 | ``float``                 |       | Maximum source redshift                      |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~event_type`                            | ``str``                   |       | Type of CBC event (BBH/BNS/NSBH)             |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~cosmo`                                 | ``astropy.cosmology``     |       | Cosmology for calculations                   |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~directory`                             | ``str``                   |       | Path for storing interpolators               |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~merger_rate_density_param`             | ``dict``                  |       | Merger rate density parameters               |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~normalization_pdf_z`                   | ``float``                 |       | Normalization constant for p(z)              |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~merger_rate_density`                   | ``callable``              |       | Merger rate density function R(z)            |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~available_merger_rate_density_model`   | ``dict``                  |       | Available merger rate density models         |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~source_redshift`                       | ``FunctionConditioning``  |       | Source redshift sampler                      |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~luminosity_distance`                   | ``FunctionConditioning``  |       | Luminosity distance interpolator             |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~differential_comoving_volume`          | ``FunctionConditioning``  |       | dVc/dz interpolator                          |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:property:: npool

      
      Number of processors for multiprocessing.



      :Returns:

          **npool** : ``int``
              Number of parallel processes to use.

              default: 4













      ..
          !! processed by numpydoc !!

   .. py:property:: z_min

      
      Minimum source redshift.



      :Returns:

          **z_min** : ``float``
              Lower bound of the redshift range.

              default: 0.001













      ..
          !! processed by numpydoc !!

   .. py:property:: z_max

      
      Maximum source redshift.



      :Returns:

          **z_max** : ``float``
              Upper bound of the redshift range.

              default: 10.0













      ..
          !! processed by numpydoc !!

   .. py:property:: directory

      
      Directory path for storing interpolator JSON files.



      :Returns:

          **directory** : ``str``
              Path to the interpolator storage directory.

              default: './interpolator_json'













      ..
          !! processed by numpydoc !!

   .. py:property:: event_type

      
      Type of compact binary coalescence event.



      :Returns:

          **event_type** : ``str``
              CBC event type.

              Options:

              - 'BBH': Binary black hole

              - 'BNS': Binary neutron star

              - 'NSBH': Neutron star-black hole

              default: 'BBH'













      ..
          !! processed by numpydoc !!

   .. py:property:: cosmo

      
      Astropy cosmology object for distance calculations.



      :Returns:

          **cosmo** : ``astropy.cosmology``
              Cosmology used for redshift-distance conversions.

              default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)













      ..
          !! processed by numpydoc !!

   .. py:property:: create_new_interpolator

      
      Dictionary controlling interpolator creation settings.



      :Returns:

          **create_new_interpolator** : ``dict``
              Dictionary with that controls the creation of new interpolators.
              Default: {'merger_rate_density': {'create_new': False, 'resolution': 100}, 'luminosity_distance': {'create_new': False, 'resolution': 100}, 'differential_comoving_volume': {'create_new': False, 'resolution': 100}}













      ..
          !! processed by numpydoc !!

   .. py:property:: luminosity_distance

      
      Class object (of FunctionConditioning) for the luminosity distance, with function as callback, which converts redshift to luminosity distance (in Mpc) for the selected cosmology.

      The class object contains the following attribute methods:

      - `function`: returns the luminosity distance distribution function.

      - `function_inverse`: returns the inverse luminosity distance distribution function, which converts luminosity distance (in Mpc) to redshift.


      :Returns:

          **luminosity_distance** : ``numpy.ndarray``
              Array of luminosity distances (in Mpc).













      ..
          !! processed by numpydoc !!

   .. py:property:: differential_comoving_volume

      
      Class object (of FunctionConditioning) for the differential comoving volume function, with function as callback, which returns dVc/dz (in Mpc^3 sr^-1) for the selected cosmology.

      The class object contains the following attribute methods:

      - `function`: returns the differential comoving volume distribution function.


      :Returns:

          **differential_comoving_volume** : ``numpy.ndarray``
              Array of differential comoving volumes (in Mpc^3 sr^-1).













      ..
          !! processed by numpydoc !!

   .. py:property:: merger_rate_density

      
      Source-frame merger rate density object.

      Returns a ``FunctionConditioning`` object with methods:

      - ``function(zs)``: Get merger rate density in source frame

      - ``rvs(size)``: Sample source redshifts in source frame

      - ``pdf(zs)``: Get probability density


      :Returns:

          **merger_rate_density** : ``FunctionConditioning``
              Callable that accepts redshift(s) and returns merger rate density
              in source frame (units: Mpc^-3 yr^-1).













      ..
          !! processed by numpydoc !!

   .. py:property:: source_redshift

      
      Class object (of FunctionConditioning) for the source redshift sampler, with rvs/sampler as callback, which samples source redshifts from p(z) ∝ R(z)/(1+z) dVc/dz , where p(z) is the redshift probability distribution, R(z) is the merger rate density, and dVc/dz is the differential comoving volume.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the source redshift distribution.

      - `pdf`: returns the source redshift probability density function.

      - `function`: returns the source redshift distribution function.


      :Returns:

          **source_redshift** : ``numpy.ndarray``
              Array of source redshifts (detector frame)













      ..
          !! processed by numpydoc !!

   .. py:property:: normalization_pdf_z

      
      Normalization constant for the redshift probability distribution.



      :Returns:

          **normalization_pdf_z** : ``float``
              Integral of the unnormalized p(z) over [z_min, z_max].













      ..
          !! processed by numpydoc !!

   .. py:property:: merger_rate_density_param

      
      Parameters for the merger rate density function.



      :Returns:

          **merger_rate_density_param** : ``dict``
              Dictionary of parameters for the selected merger rate density model.













      ..
          !! processed by numpydoc !!

   .. py:property:: available_merger_rate_density_model

      
      Dictionary of available merger rate density models and default parameters.



      :Returns:

          **available_merger_rate_density_model** : ``dict``
              Dictionary with model names as keys and parameter dicts as values.

              Available models:

              - 'merger_rate_density_bbh_oguri2018'

              - 'sfr_with_time_delay'

              - 'merger_rate_density_bbh_popIII_ken2022'

              - 'merger_rate_density_bbh_primordial_ken2022'













      ..
          !! processed by numpydoc !!

   .. py:attribute:: merger_rate_density_detector_frame

      

   .. py:method:: merger_rate_density_bbh_oguri2018(zs, get_attribute=False, **kwargs)

      
      Compute PopI/II BBH merger rate density (Oguri et al. 2018).

      Returns the source-frame merger rate density following the
      Oguri et al. (2018) prescription for PopI/II stellar populations.

      :Parameters:

          **zs** : ```numpy.ndarray``
              Source redshift(s) at which to evaluate.

          **get_attribute** : ``bool``
              If True, return the FunctionConditioning object.
              default: False

          **\*\*kwargs** : ``dict``
              Override default fitting parameters:
              R0=19e-9, b2=1.6, b3=2.1, b4=30.

      :Returns:

          **rate_density** : ```numpy.ndarray`` or ``FunctionConditioning``
              Merger rate density in source frame (units: Mpc^-3 yr^-1).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
      >>> cbc = CBCSourceRedshiftDistribution(merger_rate_density="merger_rate_density_bbh_oguri2018")
      >>> rate = cbc.merger_rate_density(zs=0.5)



      ..
          !! processed by numpydoc !!

   .. py:method:: sfr_with_time_delay(zs, get_attribute=False, **kwargs)

      
      Compute merger rate density with time delay convolution.

      Convolves the star formation rate with a time delay distribution
      to compute the merger rate density. Uses multiprocessing for
      numerical integration (Borhanian & Sathyaprakash 2024).

      :Parameters:

          **zs** : ```numpy.ndarray``
              Source redshift(s) at which to evaluate.

          **get_attribute** : ``bool``
              If True, return the FunctionConditioning object.
              default: False

          **\*\*kwargs** : ``dict``
              Override default SFR and time delay parameters.

      :Returns:

          **rate_density** : ```numpy.ndarray`` or ``FunctionConditioning``
              Merger rate density (units: Mpc^-3 yr^-1).













      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_madau_dickinson2014(zs, get_attribute=False, **kwargs)

      
      Compute star formation rate following Madau & Dickinson (2014).

      Returns the cosmic star formation rate density as given in
      Equation 15 of Madau & Dickinson (2014).

      :Parameters:

          **zs** : ```numpy.ndarray``
              Source redshift(s) at which to evaluate.

          **get_attribute** : ``bool``
              If True, return the FunctionConditioning object.
              default: False

          **\*\*kwargs** : ``dict``
              Override default fitting parameters: R0=19 * 1e-9, a=0.015, b=2.7, c=2.9, d=5.6.

      :Returns:

          **rate_density** : ```numpy.ndarray`` or ``FunctionConditioning``
              Star formation rate density (units: M_sun yr^-1 Mpc^-3).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
      >>> cbc = CBCSourceRedshiftDistribution(merger_rate_density="merger_rate_density_madau_dickinson2014")
      >>> sfr = cbc.merger_rate_density(zs=2.0)



      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_madau_dickinson_belczynski_ng(zs, get_attribute=False, **kwargs)

      
      Compute BBH merger rate density following Ng et al. (2021).

      This model uses a Madau-Dickinson-like functional form to fit the
      merger rate density of field BHs, accounting for time delays and
      metallicity effects.

      density(zs) ∝ (1 + zs) ** alpha_F / (1 + ((1 + zs) / c_F) ** beta_F)

      :Parameters:

          **zs** : ```numpy.ndarray``
              Source redshift(s) at which to evaluate.

          **get_attribute** : ``bool``
              If True, return the FunctionConditioning object.
              default: False

          **\*\*kwargs** : ``dict``
              Override default fitting parameters: R0, alpha_F, beta_F, c_F.

      :Returns:

          **rate_density** : ```numpy.ndarray`` or ``FunctionConditioning``
              Star formation rate density (units: M_sun yr^-1 Mpc^-3).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
      >>> cbc = CBCSourceRedshiftDistribution(merger_rate_density="merger_rate_density_madau_dickinson_belczynski_ng")
      >>> sfr = cbc.merger_rate_density(zs=2.0)



      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_bbh_popIII_ken2022(zs, get_attribute=False, **kwargs)

      
      Compute PopIII BBH merger rate density (Ng et al. 2022).

      Returns the merger rate density for Population III binary black
      holes following the Ng et al. (2022) prescription.

      :Parameters:

          **zs** : ```numpy.ndarray``
              Source redshift(s) at which to evaluate.

          **get_attribute** : ``bool``
              If True, return the FunctionConditioning object.
              default: False

          **\*\*kwargs** : ``dict``
              Override default fitting parameters:
              R0=19.2e-9, aIII=0.66, bIII=0.3, zIII=11.6.

      :Returns:

          **rate_density** : ```numpy.ndarray`` or ``FunctionConditioning``
              Merger rate density (units: Mpc^-3 yr^-1).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
      >>> cbc = CBCSourceRedshiftDistribution(
      ...     z_min=5, z_max=40,
      ...     merger_rate_density="merger_rate_density_bbh_popIII_ken2022"
      ... )
      >>> rate = cbc.merger_rate_density(zs=10)



      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_bbh_primordial_ken2022(zs, get_attribute=False, **kwargs)

      
      Compute primordial BBH merger rate density (Ng et al. 2022).

      Returns the merger rate density for primordial binary black holes
      following the Ng et al. (2022) prescription.

      :Parameters:

          **zs** : ```numpy.ndarray``
              Source redshift(s) at which to evaluate.

          **get_attribute** : ``bool``
              If True, return the FunctionConditioning object.
              default: False

          **\*\*kwargs** : ``dict``
              Override default fitting parameters:
              R0=0.044e-9, t0=13.786885302009708.

      :Returns:

          **rate_density** : ```numpy.ndarray`` or ``FunctionConditioning``
              Merger rate density (units: Mpc^-3 yr^-1).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
      >>> cbc = CBCSourceRedshiftDistribution(
      ...     z_min=5, z_max=40,
      ...     merger_rate_density="merger_rate_density_bbh_primordial_ken2022"
      ... )
      >>> rate = cbc.merger_rate_density(zs=10)



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



.. py:class:: CBCSourceRedshiftDistribution(npool=4, z_min=0.001, z_max=10.0, event_type='BBH', merger_rate_density=None, merger_rate_density_param=None, cosmology=None, directory='./interpolator_json', create_new_interpolator=False)


   Bases: :py:obj:`object`

   
   Class for generating compact binary coalescence source redshift distributions.

   This class generates source redshift distributions for compact binary
   coalescence events (BBH, BNS, NSBH) using various astrophysical merger rate
   density models. It provides interpolated functions for efficient sampling of
   source redshifts weighted by the merger rate density in the detector frame.

   Key Features:

   - Multiple merger rate density models (PopI/II, PopIII, Primordial)

   - Configurable cosmology for distance calculations

   - Cached interpolators for computational efficiency

   - Support for user-defined merger rate density functions

   :Parameters:

       **npool** : ``int``
           Number of processors to use for multiprocessing.

           default: 4

       **z_min** : ``float``
           Minimum redshift of the source population.

           default: 0.001

       **z_max** : ``float``
           Maximum redshift of the source population.

           default: 10.0

       **event_type** : ``str``
           Type of compact binary event.

           Options:

           - 'BBH': Binary black hole

           - 'BNS': Binary neutron star

           - 'NSBH': Neutron star-black hole

           default: 'BBH'

       **merger_rate_density** : ``str`` or ``callable`` or ``None``
           Merger rate density model to use.

           Options:

           - 'merger_rate_density_bbh_oguri2018': PopI/II BBH (Oguri 2018)

           - 'merger_rate_density_bbh_popIII_ken2022': PopIII BBH (Ng 2022)

           - 'merger_rate_density_bbh_primordial_ken2022': Primordial BBH (Ng 2022)

           - callable: User-defined function f(z) -> rate density

           default: None (uses 'merger_rate_density_bbh_oguri2018')

       **merger_rate_density_param** : ``dict`` or ``None``
           Parameters for the merger rate density function.

           default: None (uses dict(R0=19 * 1e-9, b2=1.6, b3=2.1, b4=30))

       **cosmology** : ``astropy.cosmology`` or ``None``
           Cosmology for distance calculations.

           default: None (uses LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0))

       **directory** : ``str``
           Directory to store interpolator JSON files.

           default: './interpolator_json'

       **create_new_interpolator** : ``dict`` or ``bool``
           Control interpolator creation.

           If ``bool``: Apply to all interpolators.

           If ``dict``: Per-quantity settings with keys 'create_new' and 'resolution'.

           default: False











   .. rubric:: Examples

   Basic usage:

   >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
   >>> cbc = CBCSourceRedshiftDistribution(z_min=0.001, z_max=10)
   >>> zs_samples = cbc.source_redshift(size=1000)
   >>> rate = cbc.merger_rate_density(zs=0.5)

   Instance Methods
   ----------
   CBCSourceRedshiftDistribution has the following methods:

   +-----------------------------------------------------+----------------------------------------------------+
   | Method                                              | Description                                        |
   +=====================================================+====================================================+
   | :meth:`~merger_rate_density_detector_frame`         | Compute merger rate density in detector frame      |
   +-----------------------------------------------------+----------------------------------------------------+
   | :meth:`~merger_rate_density_bbh_oguri2018`          | PopI/II merger rate density (Oguri 2018)           |
   +-----------------------------------------------------+----------------------------------------------------+
   | :meth:`~sfr_madau_dickinson2014`                    | Star formation rate (Madau & Dickinson 2014)       |
   +-----------------------------------------------------+----------------------------------------------------+
   | :meth:`~sfr_with_time_delay`                        | SFR with time delay convolution                    |
   +-----------------------------------------------------+----------------------------------------------------+
   | :meth:`~merger_rate_density_bbh_popIII_ken2022`     | PopIII merger rate density (Ng 2022)               |
   +-----------------------------------------------------+----------------------------------------------------+
   | :meth:`~merger_rate_density_bbh_primordial_ken2022` | Primordial BBH merger rate density (Ng 2022)       |
   +-----------------------------------------------------+----------------------------------------------------+

   Instance Attributes
   ----------
   CBCSourceRedshiftDistribution has the following attributes:

   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | Attribute                                      | Type                      | Unit  | Description                                  |
   +================================================+===========================+=======+==============================================+
   | :attr:`~z_min`                                 | ``float``                 |       | Minimum source redshift                      |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~z_max`                                 | ``float``                 |       | Maximum source redshift                      |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~event_type`                            | ``str``                   |       | Type of CBC event (BBH/BNS/NSBH)             |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~cosmo`                                 | ``astropy.cosmology``     |       | Cosmology for calculations                   |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~directory`                             | ``str``                   |       | Path for storing interpolators               |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~merger_rate_density_param`             | ``dict``                  |       | Merger rate density parameters               |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~normalization_pdf_z`                   | ``float``                 |       | Normalization constant for p(z)              |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~merger_rate_density`                   | ``callable``              |       | Merger rate density function R(z)            |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~available_merger_rate_density_model`   | ``dict``                  |       | Available merger rate density models         |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~source_redshift`                       | ``FunctionConditioning``  |       | Source redshift sampler                      |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~luminosity_distance`                   | ``FunctionConditioning``  |       | Luminosity distance interpolator             |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+
   | :attr:`~differential_comoving_volume`          | ``FunctionConditioning``  |       | dVc/dz interpolator                          |
   +------------------------------------------------+---------------------------+-------+----------------------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:property:: npool

      
      Number of processors for multiprocessing.



      :Returns:

          **npool** : ``int``
              Number of parallel processes to use.

              default: 4













      ..
          !! processed by numpydoc !!

   .. py:property:: z_min

      
      Minimum source redshift.



      :Returns:

          **z_min** : ``float``
              Lower bound of the redshift range.

              default: 0.001













      ..
          !! processed by numpydoc !!

   .. py:property:: z_max

      
      Maximum source redshift.



      :Returns:

          **z_max** : ``float``
              Upper bound of the redshift range.

              default: 10.0













      ..
          !! processed by numpydoc !!

   .. py:property:: directory

      
      Directory path for storing interpolator JSON files.



      :Returns:

          **directory** : ``str``
              Path to the interpolator storage directory.

              default: './interpolator_json'













      ..
          !! processed by numpydoc !!

   .. py:property:: event_type

      
      Type of compact binary coalescence event.



      :Returns:

          **event_type** : ``str``
              CBC event type.

              Options:

              - 'BBH': Binary black hole

              - 'BNS': Binary neutron star

              - 'NSBH': Neutron star-black hole

              default: 'BBH'













      ..
          !! processed by numpydoc !!

   .. py:property:: cosmo

      
      Astropy cosmology object for distance calculations.



      :Returns:

          **cosmo** : ``astropy.cosmology``
              Cosmology used for redshift-distance conversions.

              default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)













      ..
          !! processed by numpydoc !!

   .. py:property:: create_new_interpolator

      
      Dictionary controlling interpolator creation settings.



      :Returns:

          **create_new_interpolator** : ``dict``
              Dictionary with that controls the creation of new interpolators.
              Default: {'merger_rate_density': {'create_new': False, 'resolution': 100}, 'luminosity_distance': {'create_new': False, 'resolution': 100}, 'differential_comoving_volume': {'create_new': False, 'resolution': 100}}













      ..
          !! processed by numpydoc !!

   .. py:property:: luminosity_distance

      
      Class object (of FunctionConditioning) for the luminosity distance, with function as callback, which converts redshift to luminosity distance (in Mpc) for the selected cosmology.

      The class object contains the following attribute methods:

      - `function`: returns the luminosity distance distribution function.

      - `function_inverse`: returns the inverse luminosity distance distribution function, which converts luminosity distance (in Mpc) to redshift.


      :Returns:

          **luminosity_distance** : ``numpy.ndarray``
              Array of luminosity distances (in Mpc).













      ..
          !! processed by numpydoc !!

   .. py:property:: differential_comoving_volume

      
      Class object (of FunctionConditioning) for the differential comoving volume function, with function as callback, which returns dVc/dz (in Mpc^3 sr^-1) for the selected cosmology.

      The class object contains the following attribute methods:

      - `function`: returns the differential comoving volume distribution function.


      :Returns:

          **differential_comoving_volume** : ``numpy.ndarray``
              Array of differential comoving volumes (in Mpc^3 sr^-1).













      ..
          !! processed by numpydoc !!

   .. py:property:: merger_rate_density

      
      Source-frame merger rate density object.

      Returns a ``FunctionConditioning`` object with methods:

      - ``function(zs)``: Get merger rate density in source frame

      - ``rvs(size)``: Sample source redshifts in source frame

      - ``pdf(zs)``: Get probability density


      :Returns:

          **merger_rate_density** : ``FunctionConditioning``
              Callable that accepts redshift(s) and returns merger rate density
              in source frame (units: Mpc^-3 yr^-1).













      ..
          !! processed by numpydoc !!

   .. py:property:: source_redshift

      
      Class object (of FunctionConditioning) for the source redshift sampler, with rvs/sampler as callback, which samples source redshifts from p(z) ∝ R(z)/(1+z) dVc/dz , where p(z) is the redshift probability distribution, R(z) is the merger rate density, and dVc/dz is the differential comoving volume.

      The class object contains the following attribute methods:

      - `rvs`: returns random samples from the source redshift distribution.

      - `pdf`: returns the source redshift probability density function.

      - `function`: returns the source redshift distribution function.


      :Returns:

          **source_redshift** : ``numpy.ndarray``
              Array of source redshifts (detector frame)













      ..
          !! processed by numpydoc !!

   .. py:property:: normalization_pdf_z

      
      Normalization constant for the redshift probability distribution.



      :Returns:

          **normalization_pdf_z** : ``float``
              Integral of the unnormalized p(z) over [z_min, z_max].













      ..
          !! processed by numpydoc !!

   .. py:property:: merger_rate_density_param

      
      Parameters for the merger rate density function.



      :Returns:

          **merger_rate_density_param** : ``dict``
              Dictionary of parameters for the selected merger rate density model.













      ..
          !! processed by numpydoc !!

   .. py:property:: available_merger_rate_density_model

      
      Dictionary of available merger rate density models and default parameters.



      :Returns:

          **available_merger_rate_density_model** : ``dict``
              Dictionary with model names as keys and parameter dicts as values.

              Available models:

              - 'merger_rate_density_bbh_oguri2018'

              - 'sfr_with_time_delay'

              - 'merger_rate_density_bbh_popIII_ken2022'

              - 'merger_rate_density_bbh_primordial_ken2022'













      ..
          !! processed by numpydoc !!

   .. py:attribute:: merger_rate_density_detector_frame

      

   .. py:method:: merger_rate_density_bbh_oguri2018(zs, get_attribute=False, **kwargs)

      
      Compute PopI/II BBH merger rate density (Oguri et al. 2018).

      Returns the source-frame merger rate density following the
      Oguri et al. (2018) prescription for PopI/II stellar populations.

      :Parameters:

          **zs** : ```numpy.ndarray``
              Source redshift(s) at which to evaluate.

          **get_attribute** : ``bool``
              If True, return the FunctionConditioning object.
              default: False

          **\*\*kwargs** : ``dict``
              Override default fitting parameters:
              R0=19e-9, b2=1.6, b3=2.1, b4=30.

      :Returns:

          **rate_density** : ```numpy.ndarray`` or ``FunctionConditioning``
              Merger rate density in source frame (units: Mpc^-3 yr^-1).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
      >>> cbc = CBCSourceRedshiftDistribution(merger_rate_density="merger_rate_density_bbh_oguri2018")
      >>> rate = cbc.merger_rate_density(zs=0.5)



      ..
          !! processed by numpydoc !!

   .. py:method:: sfr_with_time_delay(zs, get_attribute=False, **kwargs)

      
      Compute merger rate density with time delay convolution.

      Convolves the star formation rate with a time delay distribution
      to compute the merger rate density. Uses multiprocessing for
      numerical integration (Borhanian & Sathyaprakash 2024).

      :Parameters:

          **zs** : ```numpy.ndarray``
              Source redshift(s) at which to evaluate.

          **get_attribute** : ``bool``
              If True, return the FunctionConditioning object.
              default: False

          **\*\*kwargs** : ``dict``
              Override default SFR and time delay parameters.

      :Returns:

          **rate_density** : ```numpy.ndarray`` or ``FunctionConditioning``
              Merger rate density (units: Mpc^-3 yr^-1).













      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_madau_dickinson2014(zs, get_attribute=False, **kwargs)

      
      Compute star formation rate following Madau & Dickinson (2014).

      Returns the cosmic star formation rate density as given in
      Equation 15 of Madau & Dickinson (2014).

      :Parameters:

          **zs** : ```numpy.ndarray``
              Source redshift(s) at which to evaluate.

          **get_attribute** : ``bool``
              If True, return the FunctionConditioning object.
              default: False

          **\*\*kwargs** : ``dict``
              Override default fitting parameters: R0=19 * 1e-9, a=0.015, b=2.7, c=2.9, d=5.6.

      :Returns:

          **rate_density** : ```numpy.ndarray`` or ``FunctionConditioning``
              Star formation rate density (units: M_sun yr^-1 Mpc^-3).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
      >>> cbc = CBCSourceRedshiftDistribution(merger_rate_density="merger_rate_density_madau_dickinson2014")
      >>> sfr = cbc.merger_rate_density(zs=2.0)



      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_madau_dickinson_belczynski_ng(zs, get_attribute=False, **kwargs)

      
      Compute BBH merger rate density following Ng et al. (2021).

      This model uses a Madau-Dickinson-like functional form to fit the
      merger rate density of field BHs, accounting for time delays and
      metallicity effects.

      density(zs) ∝ (1 + zs) ** alpha_F / (1 + ((1 + zs) / c_F) ** beta_F)

      :Parameters:

          **zs** : ```numpy.ndarray``
              Source redshift(s) at which to evaluate.

          **get_attribute** : ``bool``
              If True, return the FunctionConditioning object.
              default: False

          **\*\*kwargs** : ``dict``
              Override default fitting parameters: R0, alpha_F, beta_F, c_F.

      :Returns:

          **rate_density** : ```numpy.ndarray`` or ``FunctionConditioning``
              Star formation rate density (units: M_sun yr^-1 Mpc^-3).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
      >>> cbc = CBCSourceRedshiftDistribution(merger_rate_density="merger_rate_density_madau_dickinson_belczynski_ng")
      >>> sfr = cbc.merger_rate_density(zs=2.0)



      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_bbh_popIII_ken2022(zs, get_attribute=False, **kwargs)

      
      Compute PopIII BBH merger rate density (Ng et al. 2022).

      Returns the merger rate density for Population III binary black
      holes following the Ng et al. (2022) prescription.

      :Parameters:

          **zs** : ```numpy.ndarray``
              Source redshift(s) at which to evaluate.

          **get_attribute** : ``bool``
              If True, return the FunctionConditioning object.
              default: False

          **\*\*kwargs** : ``dict``
              Override default fitting parameters:
              R0=19.2e-9, aIII=0.66, bIII=0.3, zIII=11.6.

      :Returns:

          **rate_density** : ```numpy.ndarray`` or ``FunctionConditioning``
              Merger rate density (units: Mpc^-3 yr^-1).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
      >>> cbc = CBCSourceRedshiftDistribution(
      ...     z_min=5, z_max=40,
      ...     merger_rate_density="merger_rate_density_bbh_popIII_ken2022"
      ... )
      >>> rate = cbc.merger_rate_density(zs=10)



      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_bbh_primordial_ken2022(zs, get_attribute=False, **kwargs)

      
      Compute primordial BBH merger rate density (Ng et al. 2022).

      Returns the merger rate density for primordial binary black holes
      following the Ng et al. (2022) prescription.

      :Parameters:

          **zs** : ```numpy.ndarray``
              Source redshift(s) at which to evaluate.

          **get_attribute** : ``bool``
              If True, return the FunctionConditioning object.
              default: False

          **\*\*kwargs** : ``dict``
              Override default fitting parameters:
              R0=0.044e-9, t0=13.786885302009708.

      :Returns:

          **rate_density** : ```numpy.ndarray`` or ``FunctionConditioning``
              Merger rate density (units: Mpc^-3 yr^-1).










      .. rubric:: Examples

      >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
      >>> cbc = CBCSourceRedshiftDistribution(
      ...     z_min=5, z_max=40,
      ...     merger_rate_density="merger_rate_density_bbh_primordial_ken2022"
      ... )
      >>> rate = cbc.merger_rate_density(zs=10)



      ..
          !! processed by numpydoc !!


.. py:data:: chunk_size
   :value: '10000'

   

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

           default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)

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
   | :meth:`~binary_masses_BBH_powerlaw_gaussian`.       | Sample BBH masses with PowerLaw+PEAK model     |
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
              | a_2                |              | spin of the secondary compact binary |
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

.. py:function:: sample_from_powerlaw_distribution(size, alphans, mminns, mmaxns)

   
   Inverse transform sampling for a power-law mass distribution:
   p(m) ∝ m^{-alphans}, m in [mminns, mmaxns]


   :Parameters:

       **size** : int
           Number of samples to generate.

       **alphans** : float
           Power-law index (alpha).

       **mminns** : float
           Minimum neutron star mass (lower bound).

       **mmaxns** : float
           Maximum neutron star mass (upper bound).

       **random_state** : int, np.random.Generator, or None
           Seed or random generator for reproducibility.

   :Returns:

       **m** : ndarray
           Array of sampled neutron star masses.













   ..
       !! processed by numpydoc !!

.. py:function:: merger_rate_density_bbh_oguri2018_function(zs, R0=19 * 1e-09, b2=1.6, b3=2.1, b4=30)

   
   Compute the merger rate density for PopI/II BBH.

   Reference: Oguri et al. (2018). The output is in detector frame and is
   unnormalized.

   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **R0** : ``float``
           Local merger rate density at low redshift (Mpc^-3 yr^-1).

           default: 19e-9 (GWTC-4)

       **b2** : ``float``
           Fitting parameter.

           default: 1.6

       **b3** : ``float``
           Fitting parameter.

           default: 2.1

       **b4** : ``float``
           Fitting parameter.

           default: 30

   :Returns:

       **rate_density** : ``float`` or ``numpy.ndarray``
           Merger rate density.










   .. rubric:: Examples

   >>> from ler.gw_source_population import merger_rate_density_bbh_oguri2018
   >>> rate_density = merger_rate_density_bbh_oguri2018(zs=np.array([0.1]))



   ..
       !! processed by numpydoc !!

.. py:function:: merger_rate_density_bbh_popIII_ken2022_function(zs, n0=19.2 * 1e-09, aIII=0.66, bIII=0.3, zIII=11.6)

   
   Compute the unnormalized merger rate density for PopIII BBH.

   Reference: Ng et al. (2022). The output is in detector frame and is
   unnormalized.

   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **n0** : ``float``
           Normalization constant.

           default: 19.2e-9

       **aIII** : ``float``
           Fitting parameter.

           default: 0.66

       **bIII** : ``float``
           Fitting parameter.

           default: 0.3

       **zIII** : ``float``
           Characteristic redshift.

           default: 11.6

   :Returns:

       **rate_density** : ``float`` or ``numpy.ndarray``
           Merger rate density.










   .. rubric:: Examples

   >>> from ler.gw_source_population import merger_rate_density_bbh_popIII_ken2022
   >>> rate_density = merger_rate_density_bbh_popIII_ken2022(zs=np.array([0.1]))



   ..
       !! processed by numpydoc !!

.. py:function:: merger_rate_density_madau_dickinson2014_function(zs, R0=19 * 1e-09, a=0.015, b=2.7, c=2.9, d=5.6)

   
   Compute the merger rate density for BBH using Madau & Dickinson (2014) model.

   density(zs) ∝ (1 + zs) ** b / (1 + ((1 + zs) / c) ** d)

   Reference: Eqn. 15 of https://arxiv.org/pdf/1403.0007

   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **R0** : ``float``
           Local merger rate density (Mpc^-3 yr^-1).

           default: 19e-9

       **a** : ``float``
           Normalization parameter.

           default: 0.015

       **b** : ``float``
           Low-redshift power-law slope.

           default: 2.7

       **c** : ``float``
           Turnover redshift parameter.

           default: 2.9

       **d** : ``float``
           High-redshift power-law slope.

           default: 5.6

   :Returns:

       **rate_density** : ``float`` or ``numpy.ndarray``
           Merger rate density (Mpc^-3 yr^-1).










   .. rubric:: Examples

   >>> from ler.gw_source_population import merger_rate_density_madau_dickinson2014
   >>> rate_density = merger_rate_density_madau_dickinson2014(zs=np.array([0.1]))



   ..
       !! processed by numpydoc !!

.. py:function:: merger_rate_density_madau_dickinson_belczynski_ng_function(zs, R0=19 * 1e-09, alpha_F=2.57, beta_F=5.83, c_F=3.36)

   
   Compute BBH merger rate density following Ng et al. (2021).

   This model uses a Madau-Dickinson-like functional form to fit the
   merger rate density of field BHs, accounting for time delays and
   metallicity effects. Coefficients from Madau & Dickinson (2014) are translated as: B-> alpha_F, D-> beta_F, C-> c_F.

   density(zs) ∝ (1 + zs) ** alpha_F / (1 + ((1 + zs) / c_F) ** beta_F)

   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **R0** : ``float``
           Local merger rate density (Mpc^-3 yr^-1).

           default: 19e-9

       **alpha_F** : ``float``
           Low-redshift power-law slope.

           default: 2.57

       **beta_F** : ``float``
           High-redshift power-law slope.

           default: 5.83

       **c_F** : ``float``
           Turnover redshift parameter.

           default: 3.36

   :Returns:

       **rate_density** : ``float`` or ``numpy.ndarray``
           Merger rate density (Mpc^-3 yr^-1).










   .. rubric:: Examples

   >>> from ler.gw_source_population import merger_rate_density_madau_dickinson_belczynski_ng
   >>> rate_density = merger_rate_density_madau_dickinson_belczynski_ng(zs=np.array([0.1]))



   ..
       !! processed by numpydoc !!

.. py:function:: merger_rate_density_bbh_primordial_ken2022_function(zs, cosmology=None, n0=0.044 * 1e-09, t0=13.786885302009708)

   
   Compute the merger rate density for Primordial BBH.

   Reference: Ng et al. (2022). The output is in detector frame and is
   unnormalized.

   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **cosmology** : ``astropy.cosmology`` or ``None``
           Cosmology object for age calculations.

           default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)

       **n0** : ``float``
           Normalization constant.

           default: 0.044e-9

       **t0** : ``float``
           Present age of the Universe (Gyr).

           default: 13.786885302009708

   :Returns:

       **rate_density** : ``float`` or ``numpy.ndarray``
           Merger rate density.










   .. rubric:: Examples

   >>> from ler.gw_source_population import merger_rate_density_bbh_primordial_ken2022
   >>> rate_density = merger_rate_density_bbh_primordial_ken2022(zs=np.array([0.1]))



   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_fragos2017_with_bbh_td(zs, R0=19 * 1e-09)

   
   Compute the merger rate density for BBH. This is computed from star formation rate, Madau & Fragos (2017), with an additional time delay. This function is relies on pre-generated data points.


   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **R0** : ``float``
           Local merger rate density (Mpc^-3 yr^-1).

           default: 19e-9

   :Returns:

       **SFR** : ``float`` or ``numpy.ndarray``
           Star formation rate (Mpc^-3 yr^-1).













   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_dickinson2014_with_bbh_td(zs, R0=19 * 1e-09)

   
   Compute the merger rate density for BBH. This is computed from star formation rate, Madau & Dickinson (2014), with an additional time delay. This function is relies on pre-generated data points.


   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **R0** : ``float``
           Local merger rate density (Mpc^-3 yr^-1).

           default: 19e-9

   :Returns:

       **SFR** : ``float`` or ``numpy.ndarray``
           Star formation rate (Mpc^-3 yr^-1).













   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_fragos2017_with_bns_td(zs, R0=89 * 1e-09)

   
   Compute the merger rate density for BNS. This is computed from star formation rate, Madau & Fragos (2017), with an additional time delay. This function is relies on pre-generated data points.


   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **R0** : ``float``
           Local merger rate density (Mpc^-3 yr^-1).

           default: 89e-9

   :Returns:

       **SFR** : ``float`` or ``numpy.ndarray``
           Star formation rate (Mpc^-3 yr^-1).













   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_dickinson2014_with_bns_td(zs, R0=89 * 1e-09)

   
   Compute the merger rate density for BNS. This is computed from star formation rate, Madau & Dickinson (2014), with an additional time delay. This function is relies on pre-generated data points.


   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **R0** : ``float``
           Local merger rate density (Mpc^-3 yr^-1).

           default: 89e-9

   :Returns:

       **SFR** : ``float`` or ``numpy.ndarray``
           Star formation rate (Mpc^-3 yr^-1).













   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_fragos2017(zs, a=0.01, b=2.6, c=3.2, d=6.2)

   
   Compute star formation rate using Madau & Fragos (2017) model.

   Reference: https://arxiv.org/pdf/1606.07887.pdf

   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **a** : ``float``
           Normalization parameter.

           default: 0.01

       **b** : ``float``
           Low-redshift power-law slope.

           default: 2.6

       **c** : ``float``
           Turnover redshift parameter.

           default: 3.2

       **d** : ``float``
           High-redshift power-law slope.

           default: 6.2

   :Returns:

       **SFR** : ``float`` or ``numpy.ndarray``
           Star formation rate (Msun yr^-1 Mpc^-3).













   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_dickinson2014(zs, a=0.015, b=2.7, c=2.9, d=5.6)

   
   Compute star formation rate using Madau & Dickinson (2014) model.

   Reference: Eqn. 15 of https://arxiv.org/pdf/1403.0007

   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **a** : ``float``
           Normalization parameter.

           default: 0.015

       **b** : ``float``
           Low-redshift power-law slope.

           default: 2.7

       **c** : ``float``
           Turnover redshift parameter.

           default: 2.9

       **d** : ``float``
           High-redshift power-law slope.

           default: 5.6

   :Returns:

       **SFR** : ``float`` or ``numpy.ndarray``
           Star formation rate (Msun yr^-1 Mpc^-3).










   .. rubric:: Examples

   >>> from ler.gw_source_population import sfr_madau_dickinson2014
   >>> sfr = sfr_madau_dickinson2014(zs=np.array([0.1]))



   ..
       !! processed by numpydoc !!

.. py:function:: binary_masses_BBH_popIII_lognormal_rvs(size, m_min=1.0, m_max=100.0, Mc=20.0, sigma=0.3, chunk_size=10000)

   
   Sample from a lognormal distribution in 2D mass space.

   Reference: Ng et al. (2022). This is a helper function for PopIII BBH
   and primordial BBH merger rate density distribution functions.

   :Parameters:

       **size** : ``int``
           Number of samples to draw.

       **m_min** : ``float``
           Minimum mass (Msun).

           default: 1.0

       **m_max** : ``float``
           Maximum mass (Msun).

           default: 100.0

       **Mc** : ``float``
           Characteristic mass scale (Msun).

           default: 20.0

       **sigma** : ``float``
           Width of the distribution.

           default: 0.3

       **chunk_size** : ``int``
           Number of samples per rejection sampling chunk.

           default: 10000

   :Returns:

       **m1_sample** : ``numpy.ndarray``
           Primary mass samples (Msun).

       **m2_sample** : ``numpy.ndarray``
           Secondary mass samples (Msun).










   .. rubric:: Examples

   >>> from ler.gw_source_population import binary_masses_BBH_popIII_lognormal
   >>> m1, m2 = binary_masses_BBH_popIII_lognormal(size=1000)



   ..
       !! processed by numpydoc !!

.. py:function:: binary_masses_BBH_primordial_lognormal_rvs(size, m_min=1.0, m_max=100.0, Mc=20.0, sigma=0.3, chunk_size=10000)

   
   Sample from a lognormal distribution in 2D mass space.

   Based on Eqn. 1 and 4 of Ng et al. 2022 for primordial black holes.

   :Parameters:

       **size** : ``int``
           Number of samples to draw.

       **m_min** : ``float``
           Minimum mass (Msun).

           default: 1.0

       **m_max** : ``float``
           Maximum mass (Msun).

           default: 100.0

       **Mc** : ``float``
           Characteristic mass scale (Msun).

           default: 20.0

       **sigma** : ``float``
           Width of the distribution.

           default: 0.3

       **chunk_size** : ``int``
           Number of samples per rejection sampling chunk.

           default: 10000

   :Returns:

       **m1_sample** : ``numpy.ndarray``
           Primary mass samples (Msun).

       **m2_sample** : ``numpy.ndarray``
           Secondary mass samples (Msun).










   .. rubric:: Examples

   >>> from ler.gw_source_population import binary_masses_BBH_primordial_lognormal
   >>> m1, m2 = binary_masses_BBH_primordial_lognormal(size=1000)



   ..
       !! processed by numpydoc !!

.. py:function:: binary_masses_BNS_bimodal_rvs(size, w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1.0, mmax=2.3, resolution=500)

   
   Sample BNS masses from bimodal Gaussian distribution.

   Based on Will M. Farr et al. 2020 Eqn. 6 for neutron star mass
   distribution combining two Gaussian peaks.

   :Parameters:

       **size** : ``int``
           Number of samples to draw.

       **w** : ``float``
           Weight of the left (low-mass) peak.

           default: 0.643

       **muL** : ``float``
           Mean of the left peak (Msun).

           default: 1.352

       **sigmaL** : ``float``
           Standard deviation of the left peak (Msun).

           default: 0.08

       **muR** : ``float``
           Mean of the right peak (Msun).

           default: 1.88

       **sigmaR** : ``float``
           Standard deviation of the right peak (Msun).

           default: 0.3

       **mmin** : ``float``
           Minimum mass (Msun).

           default: 1.0

       **mmax** : ``float``
           Maximum mass (Msun).

           default: 2.3

       **resolution** : ``int``
           Number of points to use for the CDF.

           default: 500

   :Returns:

       **m1** : ``numpy.ndarray``
           Primary mass samples (Msun).

       **m2** : ``numpy.ndarray``
           Secondary mass samples (Msun).













   ..
       !! processed by numpydoc !!

.. py:function:: binary_masses_NSBH_broken_powerlaw_rvs(size=1000, mminbh=26.0, mmaxbh=125.0, alpha_1=6.75, alpha_2=0.0, b=0.5, delta_m=5.0, mminns=1.0, mmaxns=3.0, alphans=0.0, normalization_size=1000)

   
   Generate NSBH mass samples from broken power-law (BH) and power-law (NS).


   :Parameters:

       **size** : ``int``
           Number of samples to draw.

           default: 1000

       **mminbh** : ``float``
           Minimum BH mass (Msun).

           default: 26.0

       **mmaxbh** : ``float``
           Maximum BH mass (Msun).

           default: 125.0

       **alpha_1** : ``float``
           BH power-law index below break.

           default: 6.75

       **alpha_2** : ``float``
           BH power-law index above break.

           default: 0.0

       **b** : ``float``
           Break location parameter (0-1).

           default: 0.5

       **delta_m** : ``float``
           Smoothing width (Msun).

           default: 5.0

       **mminns** : ``float``
           Minimum NS mass (Msun).

           default: 1.0

       **mmaxns** : ``float``
           Maximum NS mass (Msun).

           default: 3.0

       **alphans** : ``float``
           NS power-law index.

           default: 0.0

       **normalization_size** : ``int``
           Grid size for CDF computation.

           default: 1000

   :Returns:

       **m1_samples** : ``numpy.ndarray``
           BH mass samples (Msun).

       **m2_samples** : ``numpy.ndarray``
           NS mass samples (Msun).













   ..
       !! processed by numpydoc !!

.. py:function:: binary_masses_BBH_powerlaw_gaussian_rvs(size, mminbh, mmaxbh, alpha, mu_g, sigma_g, lambda_peak, delta_m, beta, normalization_size=1000)

   
   Generate BBH mass samples from power-law + Gaussian model with mass ratio.


   :Parameters:

       **size** : ``int``
           Number of samples to draw.

       **mminbh** : ``float``
           Minimum BH mass (Msun).

       **mmaxbh** : ``float``
           Maximum BH mass (Msun).

       **alpha** : ``float``
           Power-law spectral index for m1.

       **mu_g** : ``float``
           Mean of the Gaussian peak (Msun).

       **sigma_g** : ``float``
           Standard deviation of the Gaussian peak (Msun).

       **lambda_peak** : ``float``
           Fraction in Gaussian component (0-1).

       **delta_m** : ``float``
           Low-mass smoothing width (Msun).

       **beta** : ``float``
           Power-law index for mass ratio distribution.

       **normalization_size** : ``int``
           Grid size for CDF computation.

           default: 1000

   :Returns:

       **m1** : ``numpy.ndarray``
           Primary mass samples (Msun).

       **m2** : ``numpy.ndarray``
           Secondary mass samples (Msun).













   ..
       !! processed by numpydoc !!

.. py:function:: available_prior_list()

   
   Returns a list of available priors.
















   ..
       !! processed by numpydoc !!

.. py:function:: sfr_madau_fragos2017(zs, a=0.01, b=2.6, c=3.2, d=6.2)

   
   Compute star formation rate using Madau & Fragos (2017) model.

   Reference: https://arxiv.org/pdf/1606.07887.pdf

   :Parameters:

       **zs** : ``float`` or ``numpy.ndarray``
           Source redshifts.

       **a** : ``float``
           Normalization parameter.

           default: 0.01

       **b** : ``float``
           Low-redshift power-law slope.

           default: 2.6

       **c** : ``float``
           Turnover redshift parameter.

           default: 3.2

       **d** : ``float``
           High-redshift power-law slope.

           default: 6.2

   :Returns:

       **SFR** : ``float`` or ``numpy.ndarray``
           Star formation rate (Msun yr^-1 Mpc^-3).













   ..
       !! processed by numpydoc !!

.. py:function:: sfr_with_time_delay_function(input_args)

   
   Compute star formation rate at observed redshift with time delay.

   The star formation rate is time-delayed relative to the observed redshift,
   with a time delay uniformly distributed between td_min and td_max. The
   formation redshift is computed using the cosmological age-redshift relation.

   :Parameters:

       **input_args** : ``list``
           List containing the following elements in order:

           - z (``float``): Observed redshift

           - idx (``int``): Index identifier for the computation

           - td_min (``float``): Minimum time delay (Gyr)

           - td_max (``float``): Maximum time delay (Gyr)

           - H0 (``float``): Hubble constant (km/s/Mpc)

           - Omega_M (``float``): Matter density parameter

           - Omega_Lambda (``float``): Dark energy density parameter

           - a (``float``): Madau-Fragos SFR normalization parameter

           - b (``float``): Madau-Fragos low-z power-law slope

           - c (``float``): Madau-Fragos turnover parameter

           - d (``float``): Madau-Fragos high-z power-law slope

   :Returns:

       **idx** : ``int``
           Index identifier (same as input).

       **result** : ``float``
           Time-averaged star formation rate at observed redshift z.










   .. rubric:: Examples

   >>> from ler.gw_source_population.sfr_with_time_delay import sfr_with_time_delay
   >>> args = [0.5, 0, 0.02, 13.0, 70.0, 0.3, 0.7, 0.01, 2.6, 3.2, 6.2]
   >>> idx, sfr = sfr_with_time_delay(args)



   ..
       !! processed by numpydoc !!

