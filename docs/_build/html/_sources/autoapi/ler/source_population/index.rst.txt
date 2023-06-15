:orphan:

:py:mod:`ler.source_population`
===============================

.. py:module:: ler.source_population


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   ler.source_population.SourceGalaxyPopulationModel
   ler.source_population.CompactBinaryPopulation




.. py:class:: SourceGalaxyPopulationModel(z_min=0.0, z_max=10.0, event_type='popI_II')

   
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
           e.g. 'popI_II', 'BNS', 'popIII', 'primordial', 'popI_II_Madau_Dickinson'
           default: 'popI_II'











   .. rubric:: Examples

   >>> from ler import SourceGalaxyPopulationModel
   >>> pop = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, event_type = "popI_II")
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
   |:meth:`~merger_rate_density_popI_II` | Function to compute the merger   |
   |                                     | rate density (PopI/PopII)        |
   +-------------------------------------+----------------------------------+
   |:meth:`~merger_rate_density_popI_II_Madau_Dickinson`                    |
   +-------------------------------------+----------------------------------+
   |                                     | Function to compute the          |
   |                                     | merger rate density (PopI/PopII) |
   |                                     | from Madau & Dickinson (2014)    |
   +-------------------------------------+----------------------------------+
   |:meth:`~merger_rate_density_popIII`  | Function to compute the merger   |
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

   .. py:method:: sample_source_redshifts(size=1000, z_min=0.0, z_max=10.0)

      
      Function to sample source redshifts from the source galaxy population
      model


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **z_min** : `float`
              Minimum redshift of the source population

          **z_max** : `float`
              Maximum redshift of the source population

      :Returns:

          **zs** : `array`
              Array of sampled redshifts













      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_popI_II(zs, R0=23.9 * 1e-09, b2=1.6, b3=2.0, b4=30)

      
      Function to compute the merger rate density (PopI/PopII)


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














      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_popI_II_Madau_Dickinson(zs, af=2.7, bf=5.6, cf=1.9)

      
      Function to compute the merger rate density (PopI/PopII) from Madau & Dickinson (2014)


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

      :Returns:

          **rate_density** : `float`
              ..













      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_popIII(zs, aIII=0.66, bIII=0.3, zIII=11.6)

      
      Function to compute the merger rate density (PopIII)


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

      :Returns:

          **rate_density** : `float`
              ..













      ..
          !! processed by numpydoc !!

   .. py:method:: merger_rate_density_primordial(zs, t0=13.786885302009708)

      
      Function to compute the merger rate density (Primordial)


      :Parameters:

          **zs** : `float`
              Source redshifts

          **t0** : `float`
              Present ge of the Universe in Gyr
              default: 13.786885302009708

      :Returns:

          **rate_density** : `float`
              ..













      ..
          !! processed by numpydoc !!


.. py:class:: CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type='popI_II', model_pars={'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, 'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, 'mu_g': 33.07, 'sigma_g': 5.69})

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

       **model_pars** : `dict`
           Dictionary of model parameters.
           e.g. for popI_II: {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, 'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, 'mu_g': 33.07, 'sigma_g': 5.69}











   .. rubric:: Examples

   >>> from ler import CompactBinaryPopulation
   >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popI_II")
   >>> gw_parameters = pop.sample_gw_parameters(nsamples=1000)
   >>> gw_parameters.keys()
   dict_keys(['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a1', 'a2', 'tilt1', 'tilt2', 'phi12', 'phi_jl'])

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
   |:attr:`~model_pars`                  | `dict`                           |
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

   .. py:attribute:: model_pars

      
      ``dict``

      Dictionary of model parameters.

      e.g. for popI_II: {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, 'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, 'mu_g': 33.07, 'sigma_g': 5.69}

      for popI_II_Madau_Dickinson: {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, 'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, 'mu_g': 33.07, 'sigma_g': 5.69}

      for popIII: None

      for primordial: {'Mc':30.,'sigma':0.3,'beta':1.1}

      for BNS: None















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
              gw_parameters.keys() = ['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a1', 'a2', 'tilt1', 'tilt2', 'phi12', 'phi_jl']













      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_popI_II(size, model_pars)

      
      Function to calculate source mass1 and mass2 with PowerLaw+PEAK model


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **model_pars** : `dict`
              Dictionary of model parameters
              e.g. {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, 'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, 'mu_g': 33.07, 'sigma_g': 5.69}

      :Returns:

          **mass_1_source** : `array`
              Array of mass1 in source frame

          **mass_2_source** : `array`
              Array of mass2 in source frame













      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_popIII(size, model_pars)

      
      Function to calculate source mass1 and mass2 with pop III origin


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **model_pars** : `dict`
              Dictionary of model parameters

      :Returns:

          **mass_1_source** : `array`
              Array of mass1 in source frame

          **mass_2_source** : `array`
              Array of mass2 in source frame













      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_primordial(size, model_pars={'Mc': 30.0, 'sigma': 0.3, 'beta': 1.1})

      
      Function to calculate source mass1 and mass2 for primordial BBHs


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **model_pars** : `dict`
              Dictionary of model parameters
              e.g. {'Mc':30.,'sigma':0.3,'beta':1.1}

      :Returns:

          **mass_1_source** : `array`
              Array of mass1 in source frame

          **mass_2_source** : `array`
              Array of mass2 in source frame













      ..
          !! processed by numpydoc !!

   .. py:method:: binary_masses_BNS(size, model_pars)

      
      Function to calculate source mass1 and mass2 of BNS


      :Parameters:

          **size** : `int`
              Number of samples to draw

          **model_pars** : `dict`
              Dictionary of model parameters

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













      ..
          !! processed by numpydoc !!


