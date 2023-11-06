:py:mod:`ler.source_population`
===============================

.. py:module:: ler.source_population

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

   ler.source_population.SourceGalaxyPopulationModel
   ler.source_population.CompactBinaryPopulation




.. py:class:: SourceGalaxyPopulationModel(z_min=0.0, z_max=10.0, event_type='BBH', merger_rate_density='merger_rate_density_bbh_popI_II_oguri', merger_rate_density_param=None)


   
   Class to generate a population of source galaxies.
   This class is inherited by :class:`~ler.ler.CompactBinaryPopulation` class.


   :Parameters:

       **z_min** : `float`
           Minimum redshift of the source population
           default: 0.

       **z_max** : `float`
           Maximum redshift of the source population
           default: 10.

       **event_type** : `str`
           Type of event to generate
           e.g. 'BBH', 'BNS', 'BBH_popIII', 'BBH_primordial', 'NSBH'
           default: 'BBH'

       **merger_rate_density** : `str`
           Type of merger rate density function to use
           default: None/'merger_rate_density_popI_II_Oguri'
           for others see instance method in :class:`~ler.ler.SourceGalaxyPopulationModel`

       **merger_rate_density_param** : `dict`
           Dictionary of merger rate density function parameters
           default: None











   .. rubric:: Examples

   >>> from ler import SourceGalaxyPopulationModel
   >>> pop = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, event_type = "BBH")
   >>> zs = pop.sample_source_redshifts(size=1000)
   >>> zs
   array([0.0001, 0.0001, 0.0001, ..., 9.9999, 9.9999, 9.9999])

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
   |:attr:`~normalization_pdf_z`         | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~z_to_luminosity_distance`    | `scipy.interpolate.interpolate`  |
   +-------------------------------------+----------------------------------+
   |:attr:`~differential_comoving_volume`| `scipy.interpolate.interpolate`  |
   +-------------------------------------+----------------------------------+

   Instance Methods
   ----------
   SourceGalaxyPopulationModel has the following instance methods:

   +-------------------------------------+----------------------------------+
   | Methods                             | Type                             |
   +=====================================+==================================+
   |:meth:`~create_lookup_table`         | Function to create a lookup      |
   |                                     | table for the differential       |
   |                                     | comoving volume and luminosity   |
   |                                     | distance wrt redshift            |
   +-------------------------------------+----------------------------------+
   |:meth:`~sample_source_redshifts`     | Function to sample source        |
   |                                     | redshifts from the source        |
   |                                     | galaxy population model          |
   +-------------------------------------+----------------------------------+
   |:meth:`~merger_rate_density_bbh_popI_II_oguri`                          |
   +-------------------------------------+----------------------------------+
   |                                     | Function to compute the merger   |
   |                                     | rate density (PopI/PopII)        |
   |                                     | from Oguri et al. (2018)         |
   +-------------------------------------+----------------------------------+
   |:meth:`~merger_rate_density_bbh_popI_II_madau_dickinson`                |
   +-------------------------------------+----------------------------------+
   |                                     | Function to compute the          |
   |                                     | merger rate density (PopI/PopII) |
   |                                     | from Madau & Dickinson (2014)    |
   +-------------------------------------+----------------------------------+
   |:meth:`~merger_rate_density_bbh_popIII`                                 |
   +-------------------------------------+----------------------------------+
   |                                     | Function to compute the merger   |
   |                                     | rate density (PopIII)            |
   +-------------------------------------+----------------------------------+
   |:meth:`~merger_rate_density_primordial`                                 |
   +-------------------------------------+----------------------------------+
   |                                     | Function to compute the merger   |
   |                                     | rate density (Primordial)        |
   +-------------------------------------+----------------------------------+



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

      e.g. 'BBH', 'BNS', 'BBH_popIII', 'BBH_primordial', 'NSBH'















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

   .. py:method:: create_lookup_table(z_min, z_max)

      
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

   .. py:method:: sample_source_redshifts(size=1000, z_min=0.0, z_max=10.0, param=None)

      
      Function to sample source redshifts (source frame) from the source galaxy population
      model


      :Parameters:

          **size** : `int`
              Number of samples to draw
              default: 1000

          **z_min** : `float`
              Minimum redshift of the source population
              default: 0.

          **z_max** : `float`
              Maximum redshift of the source population
              default: 10.

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(z_min=0.0, z_max=10.0)
              default: None

      :Returns:

          **zs** : `array`
              Array of sampled redshifts










      .. rubric:: Examples

      >>> from ler import SourceGalaxyPopulationModel
      >>> pop = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, event_type = "BBH")
      >>> zs = pop.sample_source_redshifts(size=1000)
      >>> zs
      array([0.0001, 0.0001, 0.0001, ..., 9.9999, 9.9999, 9.9999])



      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_bbh_popI_II_oguri(zs, R0=23.9 * 1e-09, b2=1.6, b3=2.0, b4=30, param=None)

      
      Function to compute the merger rate density (PopI/PopII). Reference: Oguri et al. (2018)


      :Parameters:

          **zs** : `float`
              Source redshifts

          **R0** : `float`
              Normalization constant
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

          **rate_density** : `float`
              merger rate density










      .. rubric:: Examples

      >>> from ler import SourceGalaxyPopulationModel
      >>> pop = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, event_type = "BBH", merger_rate_density="merger_rate_density_popI_II_Oguri")
      >>> rate_density = pop.merger_rate_density(zs=0.1)
      >>> rate_density
      2.7848018586883885e-08



      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_bbh_popI_II_madau_dickinson(zs, af=2.7, bf=5.6, cf=1.9, param=None)

      
      Function to compute the unormalized merger rate density (PopI/PopII) from Madau & Dickinson (2014)


      :Parameters:

          **zs** : `float`
              Source redshifts

          **af** : `float`
              Fitting paramters
              default: 2.7

          **bf** : `float`
              Fitting paramters
              default: 5.6

          **cf** : `float`
              Fitting paramters
              default: 1.9

          **param** : `dict`
              Allows to pass in above parameters as dict.
              e.g. param = dict(af=2.7, bf=5.6, cf=1.9)
              default: None

      :Returns:

          **rate_density** : `float`
              merger rate density










      .. rubric:: Examples

      >>> from ler import SourceGalaxyPopulationModel
      >>> pop = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, event_type = "BBH", merger_rate_density="merger_rate_density_bbh_popI_II_madau_dickinson")
      >>> rate_density = pop.merger_rate_density(zs=0.1)
      >>> rate_density
      1.2355851838964846



      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_popIII_ken(zs, aIII=0.66, bIII=0.3, zIII=11.6, param=None)

      
      Function to compute the unnormalized merger rate density (PopIII). Reference: Ken K. Y. Ng et al. (2022)


      :Parameters:

          **zs** : `float`
              Source redshifts

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

          **rate_density** : `float`
              merger rate density










      .. rubric:: Examples

      >>> from ler import SourceGalaxyPopulationModel
      >>> pop = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, event_type = "BBH", merger_rate_density="merger_rate_density_popIII_Ken")
      >>> rate_density = pop.merger_rate_density(zs=0.1)
      >>> rate_density
      0.00010000000000000002



      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_primordial_ken(zs, t0=13.786885302009708, param=None)

      
      Function to compute the merger rate density (Primordial). Reference: Ken K. Y. Ng et al. (2022)


      :Parameters:

          **zs** : `float`
              Source redshifts

          **t0** : `float`
              Present ge of the Universe in Gyr
              default: 13.786885302009708

      :Returns:

          **rate_density** : `float`
              merger rate density










      .. rubric:: Examples

      >>> from ler import SourceGalaxyPopulationModel
      >>> pop = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, event_type = "primordial")
      >>> rate_density = pop.merger_rate_density_primordial(zs=0.1)
      >>> rate_density
      0.00010000000000000002



      ..
          !! processed by numpydoc !!


.. py:class:: CompactBinaryPopulation(z_min=0.0001, z_max=10, event_type='BBH', event_priors=None, event_priors_params=None)


   Bases: :py:obj:`SourceGalaxyPopulationModel`

   
   Class to generate a population of compact binaries. Inherits from :class:`~ler.ler.SourceGalaxyPopulationModel` class.


   :Parameters:

       **z_min** : `float`
           Minimum redshift of the source population

       **z_max** : `float`
           Maximum redshift of the source population

       **m_min** : `float`
           Minimum mass of the BBHs

       **m_max** : `float`
           Maximum mass of the BBHs

       **event_type** : `str`
           Type of event to generate.
           e.g. 'popI_II', 'BNS', 'popIII', 'primordial', 'popI_II_Madau_Dickinson'

       **src_model_params** : `dict`
           Dictionary of model parameters.
           e.g. for popI_II: {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, 'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, 'mu_g': 33.07, 'sigma_g': 5.69}











   .. rubric:: Examples

   >>> from ler import CompactBinaryPopulation
   >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popI_II")
   >>> gw_parameters = pop.sample_gw_parameters(nsamples=1000)
   >>> gw_parameters.keys()
   dict_keys(['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl'])

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
   |:attr:`~m_min`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~m_max`                       | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`~event_type`                  | `str`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`~src_model_params`                  | `dict`                           |
   +-------------------------------------+----------------------------------+

   Instance Methods
   ----------
   CompactBinaryPopulation has the following instance methods:

   +-------------------------------------+----------------------------------+
   | Methods                             | Type                             |
   +=====================================+==================================+
   |:meth:`~sample_gw_parameters`        | Function for sampling GW         |
   |                                     | parameters from the source       |
   |                                     | galaxy population model          |
   +-------------------------------------+----------------------------------+
   |:meth:`~binary_masses_popI_II`       | Function to calculate source     |
   |                                     | mass1 and mass2 with             |
   |                                     | PowerLaw+PEAK model              |
   +-------------------------------------+----------------------------------+
   |:meth:`~binary_masses_popIII`        | Function to calculate source     |
   |                                     | mass1 and mass2 with pop III     |
   |                                     | origin                           |
   +-------------------------------------+----------------------------------+
   |:meth:`~binary_masses_primordial`    | Function to calculate source     |
   |                                     | mass1 and mass2 for primordial   |
   |                                     | BBHs                             |
   +-------------------------------------+----------------------------------+
   |:meth:`~binary_masses_BNS`           | Function to calculate source     |
   |                                     | mass1 and mass2 of BNS           |
   +-------------------------------------+----------------------------------+
   |:meth:`~mass_ratio`                  | Function to calculate mass ratio |
   +-------------------------------------+----------------------------------+



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

   .. py:attribute:: m_min

      
      ``float``

      Minimum mass of the BBHs















      ..
          !! processed by numpydoc !!

   .. py:attribute:: m_max

      
      ``float``

      Maximum mass of the BBHs















      ..
          !! processed by numpydoc !!

   .. py:attribute:: event_type

      
      ``str``

      Type of event to generate.

      e.g. 'popI_II', 'BNS', 'popIII', 'primordial', 'popI_II_Madau_Dickinson'















      ..
          !! processed by numpydoc !!

   .. py:method:: event_priors_categorization(event_type, event_priors, event_prior_params)

      
      Function to sample BBH parameters from the source galaxy population
      model


      :Parameters:

          **event_type** : `str`
              Type of event to generate.
              e.g. 'BBH', 'BNS', 'BBH_popIII', 'BBH_primordial', 'NSBH'

          **event_priors** : `dict`
              Dictionary of prior sampler functions for each parameter

          **event_prior_params** : `dict`
              Dictionary of sampler parameters for each parameter

      :Returns:

          **event_priors_** : `dict`
              Dictionary of prior sampler functions for each parameter

          **event_prior_params_** : `dict`
              Dictionary of sampler parameters for each parameter













      ..
          !! processed by numpydoc !!

   .. py:method:: sample_gw_parameters(nsamples=1000, **kwargs)

      
      Function to sample BBH parameters from the source galaxy population
      model


      :Parameters:

          **nsamples** : `int`
              Number of samples to draw

          **kwargs** : `dict`
              Keyword arguments to pass in parameter values
              e.g. zs = np.array([0.1,0.2,0.3])

      :Returns:

          **gw_parameters** : `dict`
              Dictionary of sampled parameters
              gw_parameters.keys() = ['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl']










      .. rubric:: Examples

      >>> from ler import CompactBinaryPopulation
      >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popI_II")
      >>> gw_parameters = pop.sample_gw_parameters(nsamples=1000)
      >>> gw_parameters.keys()
      dict_keys(['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl'])



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BBH_popI_II_powerlaw_gaussian(size, mminbh=4.98, mmaxbh=112.5, alpha=3.78, mu_g=32.27, sigma_g=3.88, lambda_peak=0.03, delta_m=4.8, beta=0.81, param=None)

      
      Function to calculate source mass1 and mass2 with PowerLaw+PEAK model


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **src_model_params** : `dict`
              Dictionary of model parameters
              e.g. {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, 'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, 'mu_g': 33.07, 'sigma_g': 5.69}

      :Returns:

          **mass_1_source** : `array`
              Array of mass1 in source frame

          **mass_2_source** : `array`
              Array of mass2 in source frame










      .. rubric:: Examples

      >>> from ler import CompactBinaryPopulation
      >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popI_II")
      >>> src_model_params = {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, 'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, 'mu_g': 33.07, 'sigma_g': 5.69}
      >>> mass_1_source, mass_2_source = pop.binary_masses_popI_II(size=1000, src_model_params=src_model_params)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BBH_popIII_gwcosmo(size, Mc=30.0, sigma=0.3, beta=1.1)

      
      Function to calculate source mass1 and mass2 with pop III origin


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **Mc, sigma, beta** : `float`
              Fitting parameters
              default: Mc=30.0, sigma=0.3, beta=1.1

      :Returns:

          **mass_1_source** : `array`
              Array of mass1 in source frame

          **mass_2_source** : `array`
              Array of mass2 in source frame










      .. rubric:: Examples

      >>> from ler import CompactBinaryPopulation
      >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popIII")
      >>> mass_1_source, mass_2_source = pop.binary_masses_popIII(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BBH_primordial_lognormal(size, Mc=30.0, sigma=0.3, beta=1.1)

      
      Function to calculate source mass1 and mass2 for primordial BBHs


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **Mc, sigma, beta** : `float`
              Fitting parameters
              default: Mc=30.0, sigma=0.3, beta=1.1

      :Returns:

          **mass_1_source** : `array`
              Array of mass1 in source frame

          **mass_2_source** : `array`
              Array of mass2 in source frame










      .. rubric:: Examples

      >>> from ler import CompactBinaryPopulation
      >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "primordial")
      >>> mass_1_source, mass_2_source = pop.binary_masses_primordial(size=1000)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BNS_popI_II_gwcosmo(size, mminns=1.0, mmaxns=3.0, alphans=0.0)

      
      Function to calculate source mass1 and mass2 of BNS (gwcosmo)


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **mminns** : `float`
              Minimum mass of the BNS
              default: 1.0

          **mmaxns** : `float`
              Maximum mass of the BNS
              default: 3.0

          **alphans** : `float`
              Power law index
              default: 0.0

      :Returns:

          **mass_1_source** : `array`
              Array of mass1 in source frame

          **mass_2_source** : `array`
              Array of mass2 in source frame













      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BNS_popI_II_Alsing(size, param=dict(w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1.0, mmax=2.3))

      
      Function to calculate source mass1 and mass2 of BNS (Alsing)


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **w, muL, sigmaL, muR, sigmaR** : `float`
              Fitting parameters
              default: w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3

          **mmin** : `float`
              Minimum mass of the BNS
              default: 1.0

          **mmax** : `float`
              Maximum mass of the BNS
              default: 3.0

      :Returns:

          **mass_1_source** : `array`
              Array of mass1 in source frame

          **mass_2_source** : `array`
              Array of mass2 in source frame













      ..
          !! processed by numpydoc !!

   .. py:method:: mass_ratio(size, beta=1.1)

      
      Function to calculate mass ratio with power law distribution


      :Parameters:

          **size** : `int`
              Number of samples

          **beta** : `float`
              Power law index

      :Returns:

          **q** : `array`
              Array of mass ratio










      .. rubric:: Examples

      >>> from ler import CompactBinaryPopulation
      >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=1.0, m_max=3.0, event_type = "BNS")
      >>> q = pop.mass_ratio(size=1000, beta=1.1)



      ..
          !! processed by numpydoc !!

   .. py:method:: binary_spin_BBH(size)

      
      Function to calculate spin parameters with PowerLaw+PEAK model


      :Parameters:

          **size** : `int`
              Number of samples to draw

      :Returns:

          **a_1** : `array`
              Array of spin1

          **a_2** : `array`
              Array of spin2

          **tilt_1** : `array`
              Array of tilt1

          **tilt_2** : `array`
              Array of tilt2

          **phi_12** : `array`
              Array of phi12

          **phi_jl** : `array`
              Array of phi_jl










      .. rubric:: Examples

      >>> from ler import CompactBinaryPopulation
      >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popI_II")
      >>> a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl = pop.binary_spin_BBH(size=1000)



      ..
          !! processed by numpydoc !!


