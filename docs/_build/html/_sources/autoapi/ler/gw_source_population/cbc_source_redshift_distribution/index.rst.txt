:py:mod:`ler.gw_source_population.cbc_source_redshift_distribution`
===================================================================

.. py:module:: ler.gw_source_population.cbc_source_redshift_distribution

.. autoapi-nested-parse::

   This module contains the classes to generate source galaxy population model
   and compact binary population model. The compact binary population model
   inherits from the source galaxy population model. The source galaxy population
   model is used to generate the source galaxy population model. The compact binary
   population model is used to generate the compact binary population model.

   ..
       !! processed by numpydoc !!


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   ler.gw_source_population.cbc_source_redshift_distribution.CBCSourceRedshiftDistribution




.. py:class:: CBCSourceRedshiftDistribution(z_min=0.001, z_max=10.0, event_type='BBH', merger_rate_density='merger_rate_density_bbh_popI_II_oguri2018', merger_rate_density_param=None, cosmology=None, directory='./interpolator_pickle', create_new_interpolator=False)


   Bases: :py:obj:`object`

   
   Class to generate a population of source galaxies.
   This class is inherited by :class:`~ler.ler.CompactBinaryPopulation` and :class:`~ler.ler.LensGalaxyParameterDistribution` class.


   :Parameters:

       **z_min** : `float`
           Minimum redshift of the source population
           default: 0.

       **z_max** : `float`
           Maximum redshift of the source population
           default: 10.

       **event_type** : `str`
           Type of event to generate.
           e.g. 'BBH', 'BNS', 'NSBH'

       **cosmology** : `astropy.cosmology`
           Cosmology to use
           default: None/astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)

       **merger_rate_density** : `str` or `function`
           Type of merger rate density function to use
           default: 'merger_rate_density_popI_II_oguri2018'
           for others see instance method in :class:`~ler.ler.merger_rate_density_model_list`

       **merger_rate_density_param** : `dict`
           Dictionary of merger rate density function parameters
           default: None/dict(R0=25 * 1e-9, b2=1.6, b3=2.0, b4=30)

       **directory** : `str`
           Directory to store the interpolator pickle files
           default: './interpolator_pickle'

       **create_new_interpolator** : `dict`
           Dictionary of interpolator creation parameters
           default: None/dict(redshift_distribution=dict(create_new=False, resolution=1000), z_to_luminosity_distance=dict(create_new=False, resolution=1000), differential_comoving_volume=dict(create_new=False, resolution=1000))











   .. rubric:: Examples

   >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
   >>> cbc = CBCSourceRedshiftDistribution(z_min=0.001, z_max=10, merger_rate_density="merger_rate_density_bbh_popI_II_oguri2018")
   >>> cbc.merger_rate_density(zs=0.0001) # local merger rate density at low redshift

   Instance Attributes
   ----------
   SourceGalaxyPopulationModel has the following instance attributes:

   +-------------------------------------+----------------------------------+
   | Atrributes                          | Type                             |
   +=====================================+==================================+
   |:attr:`~z_min`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~z_max`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~event_type`                  | `str`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~cosmo`                       | `astropy.cosmology`              |
   +-------------------------------------+----------------------------------+
   |:attr:`~merger_rate_density_param`   | `dict`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`~normalization_pdf_z`         | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~merger_rate_density_model_list`                                 |
   +-------------------------------------+----------------------------------+
   |                                     | List of available                |
   |                                     | merger rate density functions    |
   |                                     | and its parameters               |
   +-------------------------------------+----------------------------------+
   |:attr:`~sample_source_redshift`      | Function to sample source        |
   |                                     | redshifts                        |
   +-------------------------------------+----------------------------------+

   Instance Methods
   ----------
   SourceGalaxyPopulationModel has the following instance methods:

   +-------------------------------------+----------------------------------+
   | Methods                             | Type                             |
   +=====================================+==================================+
   |:attr:`~merger_rate_density`         | `function`                       |
   +-------------------------------------+----------------------------------+
   |:attr:`~z_to_luminosity_distance`    | `function`                       |
   +-------------------------------------+----------------------------------+
   |:attr:`~differential_comoving_volume`| `function`                       |
   +-------------------------------------+----------------------------------+
   |:meth:`~pdf_z`                       | Function to compute the pdf      |
   |                                     | p(z)                             |
   +-------------------------------------+----------------------------------+
   |:meth:`~merger_rate_density_detector_frame`                                  |
   +-------------------------------------+----------------------------------+
   |                                     | Function to compute the merger   |
   |                                     | rate density (detector frame)    |
   +-------------------------------------+----------------------------------+
   |:meth:`~create_lookup_table`         | Function to create a lookup      |
   |                                     | table for the differential       |
   |                                     | comoving volume and luminosity   |
   |                                     | distance wrt redshift            |
   +-------------------------------------+----------------------------------+
   |:meth:`~merger_rate_density_bbh_popI_II_oguri2018`                      |
   +-------------------------------------+----------------------------------+
   |                                     | Function to compute the merger   |
   |                                     | rate density (PopI/PopII)        |
   |                                     | from Oguri et al. (2018)         |
   +-------------------------------------+----------------------------------+
   |:meth:`~star_formation_rate_madau_dickinson2014`                        |
   +-------------------------------------+----------------------------------+
   |                                     | Function to compute star         |
   |                                     | formation rate as given in       |
   |                                     | Eqn. 15 Madau & Dickinson (2014) |
   +-------------------------------------+----------------------------------+
   |:meth:`~merger_rate_density_bbh_popIII_ken2022`                         |
   +-------------------------------------+----------------------------------+
   |                                     | Function to compute the merger   |
   |                                     | rate density (PopIII)            |
   +-------------------------------------+----------------------------------+
   |:meth:`~merger_rate_density_bbh_primordial_ken2022`                     |
   +-------------------------------------+----------------------------------+
   |                                     | Function to compute the merger   |
   |                                     | rate density (Primordial)        |
   +-------------------------------------+----------------------------------+



   ..
       !! processed by numpydoc !!
   .. py:property:: sample_source_redshift

      
      Function to sample source redshifts (detector frame) between z_min and z_max from the detector galaxy population


      :Parameters:

          **size** : `int`
              Number of samples to generate

      :Returns:

          zs
              Source redshifts













      ..
          !! processed by numpydoc !!

   .. py:property:: merger_rate_density

      
      Function to get the merger rate density function wrt redshift. The output is in detector frame and is unnormalized.


      :Parameters:

          **zs** : `float`
              1D array of floats
              Source redshifts

      :Returns:

          **merger_rate_density** : `float`
              merger rate density in detector frame (Mpc^-3 yr^-1)










      .. rubric:: Examples

      >>> from ler.gw_source_population import SourceGalaxyPopulationModel
      >>> cbc = SourceGalaxyPopulationModel()
      >>> merger_rate_density = cbc.merger_rate_density(zs=0.1)



      ..
          !! processed by numpydoc !!

   .. py:property:: merger_rate_density_model_list

      
      Dictionary of available merger rate density functions and its parameters.
















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

   .. py:attribute:: cosmo

      
      ``astropy.cosmology``

      Cosmology to use for the redshift distribution.

      e.g. Planck18, WMAP9, FlatLambdaCDM(H0=70, Om0=0.3) etc.















      ..
          !! processed by numpydoc !!

   .. py:attribute:: merger_rate_density_param

      
      ``dict``

      Dictionary of merger rate density function input parameters















      ..
          !! processed by numpydoc !!

   .. py:attribute:: c_n_i

      
      ``dict``

      c_n_i stands for 'create new interpolator'. Dictionary of interpolator creation parameters.

      e.g. dict(redshift_distribution=dict(create_new=False, resolution=1000), z_to_luminosity_distance=dict(create_new=False, resolution=1000), differential_comoving_volume=dict(create_new=False, resolution=1000))















      ..
          !! processed by numpydoc !!

   .. py:attribute:: normalization_pdf_z

      
      ``float``

      Normalization constant of the pdf p(z)















      ..
          !! processed by numpydoc !!

   .. py:attribute:: z_to_luminosity_distance

      
      ``scipy.interpolate.interpolate``

      Function to convert redshift to luminosity distance















      ..
          !! processed by numpydoc !!

   .. py:attribute:: differential_comoving_volume

      
      ``scipy.interpolate.interpolate``

      Function to calculate the differential comoving volume















      ..
          !! processed by numpydoc !!

   .. py:method:: pdf_z(zs, param=None)

      
      Function to compute the pdf p(z). The output is in detector frame and is normalized.


      :Parameters:

          **zs** : `float` or `numpy.ndarray` (1D array of floats)
              Source redshifts

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. if the merger_rate_density is merger_rate_density_bbh_popI_II_oguri2018
              param = dict(R0=23.9*1e-9, b2=1.6, b3=2.0, b4=30)

      :Returns:

          **pdf** : `numpy.ndarray`
              1D array of floats
              pdf p(z)










      .. rubric:: Examples

      >>> from ler.gw_source_population import SourceGalaxyPopulationModel
      >>> cbc = SourceGalaxyPopulationModel()
      >>> pdf = cbc.pdf_z(zs=0.1)



      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_detector_frame(zs, param=None)

      
      Function to compute the merger rate density (detector frame). The output is in detector frame and is unnormalized.


      :Parameters:

          **zs** : `float` or `numpy.ndarray` (1D array of floats)
              Source redshifts

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. if the merger_rate_density is merger_rate_density_bbh_popI_II_oguri2018
              param = dict(R0=23.9*1e-9, b2=1.6, b3=2.0, b4=30)

      :Returns:

          **rate_density** : `numpy.ndarray`
              1D array of floats
              merger rate density (detector frame) (Mpc^-3 yr^-1)










      .. rubric:: Examples

      >>> from ler.gw_source_population import SourceGalaxyPopulationModel
      >>> cbc = SourceGalaxyPopulationModel()
      >>> rate_density = cbc.merger_rate_density_detector_frame(zs=0.1)



      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_bbh_popI_II_oguri2018(zs, R0=23.9 * 1e-09, b2=1.6, b3=2.0, b4=30, param=None)

      
      Function to compute the merger rate density (PopI/PopII). Reference: Oguri et al. (2018). The output is in source frame and is unnormalized.


      :Parameters:

          **zs** : `float` or `numpy.ndarray` (nD array of floats)
              Source redshifts

          **R0** : `float`
              local merger rate density at low redshift
              default: 23.9*1e-9 Mpc^-3 yr^-1

          **b2** : `float`
              Fitting paramters
              default: 1.6

          **b3** : `float`
              Fitting paramters
              default: 2.0

          **b4** : `float`
              Fitting paramters
              default: 30

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(R0=23.9*1e-9, b2=1.6, b3=2.0, b4=30)
              default: None

      :Returns:

          **rate_density** : `float` or `numpy.ndarray` (nD array of floats)
              merger rate density










      .. rubric:: Examples

      >>> from ler.gw_source_population import SourceGalaxyPopulationModel
      >>> pop = SourceGalaxyPopulationModel(z_min=0.0, z_max=10, event_type = "BBH", merger_rate_density="merger_rate_density_bbh_popI_II_oguri2018")
      >>> rate_density = pop.merger_rate_density(zs=0.1)



      ..
          !! processed by numpydoc !!

   .. py:method:: star_formation_rate_madau_dickinson2014(zs, af=2.7, bf=5.6, cf=2.9, param=None)

      
      Function to compute star formation rate as given in Eqn. 15 Madau & Dickinson (2014). The output is in detector frame and is unnormalized.


      :Parameters:

          **zs** : `float` or `numpy.ndarray` (nD array of floats)
              Source redshifts

          **af** : `float`
              Fitting paramters
              default: 2.7

          **bf** : `float`
              Fitting paramters
              default: 5.6

          **cf** : `float`
              Fitting paramters
              default: 2.9

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(af=2.7, bf=5.6, cf=2.9)
              default: None

      :Returns:

          **rate_density** : `float` or `numpy.ndarray` (nD array of floats)
              merger rate density










      .. rubric:: Examples

      >>> from ler.gw_source_population import SourceGalaxyPopulationModel
      >>> pop = SourceGalaxyPopulationModel(z_min=5., z_max=40., event_type = "BBH", merger_rate_density="star_formation_rate_madau_dickinson2014")
      >>> rate_density = pop.merger_rate_density(zs=10)



      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_bbh_popIII_ken2022(zs, n0=19.2 * 1e-09, aIII=0.66, bIII=0.3, zIII=11.6, param=None)

      
      Function to compute the unnormalized merger rate density (PopIII). Reference: Ng et al. 2022. The output is in detector frame and is unnormalized.


      :Parameters:

          **zs** : `float` or `numpy.ndarray` (nD array of floats)
              Source redshifts

          **n0** : `float`
              normalization constant
              default: 19.2*1e-9

          **aIII** : `float`
              Fitting paramters
              default: 0.66

          **bIII** : `float`
              Fitting paramters
              default: 0.3

          **zIII** : `float`
              Fitting paramters
              default: 11.6

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(aIII=0.66, bIII=0.3, zIII=11.6)
              default: None

      :Returns:

          **rate_density** : `float` or `numpy.ndarray` (nD array of floats)
              merger rate density










      .. rubric:: Examples

      >>> from ler.gw_source_population import SourceGalaxyPopulationModel
      >>> pop = SourceGalaxyPopulationModel(z_min=5, z_max=40, event_type = "BBH", merger_rate_density="merger_rate_density_popIII_ken2022")
      >>> rate_density = pop.merger_rate_density(zs=10)
      >>> rate_density  # Mpc^-3 yr^-1
      1.5107979464621443e-08



      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_bbh_primordial_ken2022(zs, n0=0.044 * 1e-09, t0=13.786885302009708, param=None)

      
      Function to compute the merger rate density (Primordial). Reference: Ng et al. 2022. The output is in detector frame and is unnormalized.


      :Parameters:

          **zs** : `float`
              Source redshifts

          **n0** : `float`
              normalization constant
              default: 0.044*1e-9

          **t0** : `float`
              Present age of the Universe in Gyr
              default: 13.786885302009708

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(t0=13.786885302009708)

      :Returns:

          **rate_density** : `float`
              merger rate density










      .. rubric:: Examples

      >>> from ler.gw_source_population import SourceGalaxyPopulationModel
      >>> pop = SourceGalaxyPopulationModel(z_min=5, z_max=40, event_type = "BBH", merger_rate_density="merger_rate_density_bbh_primordial_ken2022")
      >>> rate_density = pop.merger_rate_density(zs=10)
      >>> rate_density  # Mpc^-3 yr^-1
      9.78691173794454e-10



      ..
          !! processed by numpydoc !!

   .. py:method:: create_lookup_table(z_min, z_max, directory)

      
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


