:py:mod:`ler.image_properties.image_properties`
===============================================

.. py:module:: ler.image_properties.image_properties

.. autoapi-nested-parse::

   This module contains the LensGalaxyPopulation class, which is used to sample lens galaxy parameters, source parameters conditioned on the source being strongly lensed, image properties, and lensed SNRs.

   The class inherits from the CBCSourceParameterDistribution class, which is used to sample source parameters.

   ..
       !! processed by numpydoc !!


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   ler.image_properties.image_properties.ImageProperties




.. py:class:: ImageProperties(npool=4, z_min=0.0, z_max=10, n_min_images=2, n_max_images=4, geocent_time_min=1126259462.4, geocent_time_max=1126259462.4 + 365 * 24 * 3600 * 20, lens_model_list=['EPL_NUMBA', 'SHEAR'], cosmology=None, spin_zero=True, spin_precession=False, directory='./interpolator_json', create_new_interpolator=False)


   
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
           default: dict(luminosity_distance_to_z=dict(create_new=False, resolution=1000))











   .. rubric:: Examples

   >>> from ler.image_properties import ImageProperties
   >>> image_properties = ImageProperties()
   >>> lens_parameters = dict(zs=2.0, zl=0.5, gamma1=0.0, gamma2=0.0, e1=0.0, e2=0.0, gamma=2.0, theta_E=1.0)
   >>> lens_parameters = image_properties.image_properties(lens_parameters)
   >>> print(lens_parameters.keys())

   Instance Attributes
   ----------
   ImageProperties has the following instance attributes:

   +-------------------------------------+----------------------------------+
   | Atrributes                          | Type                             |
   +=====================================+==================================+
   |:attr:`npool`                        | `int`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`z_min`                        | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`z_max`                        | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`n_min_images`                 | `int`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`n_max_images`                 | `int`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`geocent_time_min`             | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`geocent_time_max`             | `float`                          |
   +-------------------------------------+----------------------------------+
   |:attr:`lens_model_list`              | `list`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`cosmo`                        | `astropy.cosmology`              |
   +-------------------------------------+----------------------------------+
   |:attr:`spin_zero`                    | `bool`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`spin_precession`              | `bool`                           |
   +-------------------------------------+----------------------------------+
   |:attr:`directory`                    | `str`                            |
   +-------------------------------------+----------------------------------+
   |:attr:`create_new_interpolator`      | `dict`                           |
   +-------------------------------------+----------------------------------+



   ..
       !! processed by numpydoc !!
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

              image related=>['x_source': source position in the x-direction, 'y_source': source position in the y-direction, 'x0_image_position': image position in the x-direction, 'x1_image_position': image position in the y-direction, 'magnifications': magnifications, 'time_delays': time delays, 'n_images': number of images formed, 'determinant': determinants, 'trace': traces, 'iteration': to keep track of the iteration number













      ..
          !! processed by numpydoc !!

   .. py:method:: get_lensed_snrs(lensed_param, list_of_detectors=None, snr_calculator=None, pdet_calculator=None)

      
      Function to calculate the signal to noise ratio for each image in each event.


      :Parameters:

          **snr_calculator** : `function`
              snr function, as describe in the :class:`~ler.rates.GWRATES` class.

          **list_of_detectors** : `list`
              list of detectors
              e.g. ['H1', 'L1', 'V1']

          **lensed_param** : `dict`
              dictionary containing the both already lensed source paramters and image parameters.
              e.g. lensed_param.keys() = ['mass_1', 'mass_2', 'zs', 'luminosity_distance', 'theta_jn', 'psi', 'phi', 'ra', 'dec', 'geocent_time', 'phase', 'a_1', 'a2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl', 'magnifications', 'time_delays']

          **n_max_images** : `int`
              maximum number of images to consider
              default: 4

      :Returns:

          **snrs** : `dict`
              signal to noise ratio for each image in each event.
              (dictionary containing 'H1', 'L1', ..., and 'snr_net', which is the network snr, for each image as an array with dimensions (number_of_lensed_events,n_max_images) )













      ..
          !! processed by numpydoc !!


