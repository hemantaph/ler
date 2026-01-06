:py:mod:`ler.image_properties`
==============================

.. py:module:: ler.image_properties


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   image_properties/index.rst
   multiprocessing_routine/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   ler.image_properties.FunctionConditioning
   ler.image_properties.ImageProperties



Functions
~~~~~~~~~

.. autoapisummary::

   ler.image_properties.solve_lens_equation
   ler.image_properties.phi_q2_ellipticity_hemanta
   ler.image_properties.luminosity_distance
   ler.image_properties.solve_lens_equation



Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.image_properties.cosmo


.. py:data:: cosmo

   

.. py:function:: solve_lens_equation(lens_parameters)

   
   Function to solve the lens equation (min_image = 2)


   :Parameters:

       **lens_parameters** : `list`
           a list of parameters
           lens_parameters[0] = min_images : minimum number of images
           lens_parameters[1] = e1 : ellipticity
           lens_parameters[2] = e2 : ellipticity
           lens_parameters[3] = gamma : power-law index
           lens_parameters[4] = gamma1 : shear
           lens_parameters[5] = gamma2 : shear
           lens_parameters[6] = zl : redshift of the lens
           lens_parameters[7] = zs : redshift of the source
           lens_parameters[8] = einstein_radius : Einstein radius
           lens_parameters[9] = iteration : iteration number
           lens_parameters[10:] = lens_model_list : numpy array of lens models

   :Returns:

       **x_source** : `float`
           x position of the source in the source plane, unit: arcsec

       **y_source** : `float`
           y position of the source in the source plane

       **x0_image_position** : `float`
           x position of the images in the source plane

       **x1_image_position** : `float`
           y position of the images in the source plane

       **magnifications** : `float`
           magnification of the images

       **time_delays** : `float`
           time-delay of the images

       **nImages** : `int`
           number of images

       **determinant** : `float`
           determinant of the hessian matrix

       **trace** : `float`
           trace of the hessian matrix

       **iteration** : `int`
           iteration number










   .. rubric:: Examples

   >>> from ler.image_properties.multiprocessing_routine import solve_lens_equation
   >>> import numpy as np
   >>> from multiprocessing import Pool
   >>> # lens parameters input contains 12 parameters [e1, e2, gamma, gamma1, gamma2, zl, zs, einstein_radius, iteration, lens_model_list]
   >>> lens_parameters1 = np.array([2, 0.024069457093642648, -0.016002190961948142, 1.8945414936459974, 0.10117465203892329, 0.09600089396968613, 0.2503743800068136, 0.9418211055453296, 2.5055790287104725e-06, 0, 'EPL_NUMBA', 'SHEAR'], dtype=object)
   >>> lens_parameters2 = np.array([2, -0.04030088581646998, -0.01419438113690042, 2.0068239327017, 0.08482718989370612, -0.015393332086560785, 1.0952303138971118, 2.5534097159384417, 1.0125570159563301e-06, 1, 'EPL_NUMBA', 'SHEAR'], dtype=object)
   >>> input_arguments = np.vstack((lens_parameters1, lens_parameters2))
   >>> # solve the lens equation for each set of lens parameters
   >>> with Pool(2) as p:
   ...     result = p.map(solve_lens_equation1, input_arguments)
   >>> # result is a list of tuples
   >>> # each tuple contains the output parameters of the function
   >>> # each output parameter contains x_source, y_source, x0_image_position, x1_image_position, magnifications, time_delays, nImages, determinant, trace, iteration
   >>> print(f"magnification of images with lens parameters 'lens_parameters1' is {result[0][6]}")
   magnification of images with lens parameters 'lens_parameters1' is [ 2.18973765 -1.27542831]



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



.. py:function:: luminosity_distance(z=None, z_min=0.001, z_max=10.0, cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7), directory='./interpolator_json', create_new=False, resolution=500, get_attribute=True)

   
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


.. py:function:: solve_lens_equation(lens_parameters)

   
   Function to solve the lens equation (min_image = 2)


   :Parameters:

       **lens_parameters** : `list`
           a list of parameters
           lens_parameters[0] = min_images : minimum number of images
           lens_parameters[1] = e1 : ellipticity
           lens_parameters[2] = e2 : ellipticity
           lens_parameters[3] = gamma : power-law index
           lens_parameters[4] = gamma1 : shear
           lens_parameters[5] = gamma2 : shear
           lens_parameters[6] = zl : redshift of the lens
           lens_parameters[7] = zs : redshift of the source
           lens_parameters[8] = einstein_radius : Einstein radius
           lens_parameters[9] = iteration : iteration number
           lens_parameters[10:] = lens_model_list : numpy array of lens models

   :Returns:

       **x_source** : `float`
           x position of the source in the source plane, unit: arcsec

       **y_source** : `float`
           y position of the source in the source plane

       **x0_image_position** : `float`
           x position of the images in the source plane

       **x1_image_position** : `float`
           y position of the images in the source plane

       **magnifications** : `float`
           magnification of the images

       **time_delays** : `float`
           time-delay of the images

       **nImages** : `int`
           number of images

       **determinant** : `float`
           determinant of the hessian matrix

       **trace** : `float`
           trace of the hessian matrix

       **iteration** : `int`
           iteration number










   .. rubric:: Examples

   >>> from ler.image_properties.multiprocessing_routine import solve_lens_equation
   >>> import numpy as np
   >>> from multiprocessing import Pool
   >>> # lens parameters input contains 12 parameters [e1, e2, gamma, gamma1, gamma2, zl, zs, einstein_radius, iteration, lens_model_list]
   >>> lens_parameters1 = np.array([2, 0.024069457093642648, -0.016002190961948142, 1.8945414936459974, 0.10117465203892329, 0.09600089396968613, 0.2503743800068136, 0.9418211055453296, 2.5055790287104725e-06, 0, 'EPL_NUMBA', 'SHEAR'], dtype=object)
   >>> lens_parameters2 = np.array([2, -0.04030088581646998, -0.01419438113690042, 2.0068239327017, 0.08482718989370612, -0.015393332086560785, 1.0952303138971118, 2.5534097159384417, 1.0125570159563301e-06, 1, 'EPL_NUMBA', 'SHEAR'], dtype=object)
   >>> input_arguments = np.vstack((lens_parameters1, lens_parameters2))
   >>> # solve the lens equation for each set of lens parameters
   >>> with Pool(2) as p:
   ...     result = p.map(solve_lens_equation1, input_arguments)
   >>> # result is a list of tuples
   >>> # each tuple contains the output parameters of the function
   >>> # each output parameter contains x_source, y_source, x0_image_position, x1_image_position, magnifications, time_delays, nImages, determinant, trace, iteration
   >>> print(f"magnification of images with lens parameters 'lens_parameters1' is {result[0][6]}")
   magnification of images with lens parameters 'lens_parameters1' is [ 2.18973765 -1.27542831]



   ..
       !! processed by numpydoc !!

