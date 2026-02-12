:py:mod:`ler.gw_source_population.cbc_source_redshift_distribution`
===================================================================

.. py:module:: ler.gw_source_population.cbc_source_redshift_distribution

.. autoapi-nested-parse::

   Module for compact binary coalescence (CBC) source redshift distribution.

   This module provides the :class:`~CBCSourceRedshiftDistribution` class for
   generating redshift distributions of compact binary sources (BBH, BNS, NSBH)
   based on various merger rate density models. It supports multiple astrophysical
   merger rate density prescriptions including PopI/II, PopIII, and primordial
   black hole models.

   Usage:
       Basic workflow example:

       >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
       >>> cbc = CBCSourceRedshiftDistribution(z_min=0.001, z_max=10)
       >>> zs_samples = cbc.source_redshift(size=1000)

   Copyright (C) 2024 Hemantakumar Phurailatpam. Distributed under MIT License.

   ..
       !! processed by numpydoc !!


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   ler.gw_source_population.cbc_source_redshift_distribution.CBCSourceRedshiftDistribution




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
   | :meth:`~merger_rate_density_bbh_oguri2018`  | PopI/II merger rate density (Oguri 2018)           |
   +-----------------------------------------------------+----------------------------------------------------+
   | :meth:`~sfr_madau_dickinson2014`                    | Star formation rate (Madau & Dickinson 2014)       |
   +-----------------------------------------------------+----------------------------------------------------+
   | :meth:`~sfr_with_time_delay`                                | SFR with time delay convolution                    |
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


