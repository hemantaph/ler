:orphan:

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

   ler.image_properties.ImageProperties



Functions
~~~~~~~~~

.. autoapisummary::

   ler.image_properties.solve_lens_equation1
   ler.image_properties.solve_lens_equation2
   ler.image_properties.solve_lens_equation1
   ler.image_properties.solve_lens_equation2



.. py:function:: solve_lens_equation1(lens_parameters)

   
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
           x position of the source in the source plane

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

       **weights** : `float`
           weights for the caustic










   .. rubric:: Examples

   >>> from ler.multiprocessing_routine import solve_lens_equation1
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
   >>> # each output parameter contains x_source, y_source, x0_image_position, x1_image_position, magnifications, time_delays, nImages, determinant, trace, iteration, weights
   >>> print(f"magnification of images with lens parameters 'lens_parameters1' is {result[0][6]}")
   magnification of images with lens parameters 'lens_parameters1' is [ 2.18973765 -1.27542831]



   ..
       !! processed by numpydoc !!

.. py:function:: solve_lens_equation2(lens_parameters)

   
   Function to solve the lens equation (min_image > 2)


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
           x position of the source in the source plane

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

       **weights** : `float`
           weights for the caustic










   .. rubric:: Examples

   >>> from ler.multiprocessing_routine import solve_lens_equation2
   >>> import numpy as np
   >>> from multiprocessing import Pool
   >>> # lens parameters input contains 12 parameters [e1, e2, gamma, gamma1, gamma2, zl, zs, einstein_radius, iteration, lens_model_list]
   >>> lens_parameters1 = np.array([3, 0.024069457093642648, -0.016002190961948142, 1.8945414936459974, 0.10117465203892329, 0.09600089396968613, 0.2503743800068136, 0.9418211055453296, 2.5055790287104725e-06, 0, 'EPL_NUMBA', 'SHEAR'], dtype=object)
   >>> lens_parameters2 = np.array([3, -0.04030088581646998, -0.01419438113690042, 2.0068239327017, 0.08482718989370612, -0.015393332086560785, 1.0952303138971118, 2.5534097159384417, 1.0125570159563301e-06, 1, 'EPL_NUMBA', 'SHEAR'], dtype=object)
   >>> input_arguments = np.vstack((lens_parameters1, lens_parameters2))
   >>> # solve the lens equation for each set of lens parameters
   >>> with Pool(2) as p:
   ...     result = p.map(solve_lens_equation2, input_arguments)
   >>> # result is a list of tuples
   >>> # each tuple contains the output parameters of the function
   >>> # each output parameter contains x_source, y_source, x0_image_position, x1_image_position, magnifications, time_delays, nImages, determinant, trace, iteration, weights
   >>> print(f"magnification of images with lens parameters 'lens_parameters1' is {result[0][6]}")
   magnification of images with lens parameters 'lens_parameters1' is [ 2.18973765 -1.27542831]



   ..
       !! processed by numpydoc !!

.. py:class:: ImageProperties(npool=4, n_min_images=2, n_max_images=4, lens_model_list=['EPL_NUMBA', 'SHEAR'])


   
   Class to find the image properties of a lensed event. Image properties include image positions, magnifications, time delays, etc.


   :Parameters:

       **npool** : `int`
           number of processes to use
           default: 4

       **n_min_images** : `int`
           minimum number of images to consider
           default: 2

       **n_max_images** : `int`
           maximum number of images to consider
           default: 4

       **lens_model_list** : `list`
           list of lens models
           default: ['EPL_NUMBA', 'SHEAR']














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

              source related=>['mass_1': mass in detector frame (mass1>mass2), 'mass_2': mass in detector frame, 'mass_1_source':mass in source frame, 'mass_2_source':mass source frame, 'luminosity_distance': luminosity distance, 'iota': inclination angle, 'psi': polarization angle, 'phase': coalesence phase, 'geocent_time': coalensence GPS time at geocenter, 'ra': right ascension, 'dec': declination, 'a_1': spin magnitude of the more massive black hole, 'a2': spin magnitude of the less massive black hole, 'tilt_1': tilt angle of the more massive black hole, 'tilt_2': tilt angle of the less massive black hole, 'phi_12': azimuthal angle between the two spins, 'phi_jl': azimuthal angle between the total angular momentum and the orbital angular momentum]

              image related=>['x_source': source position in the x-direction, 'y_source': source position in the y-direction, 'x0_image_position': image position in the x-direction, 'x1_image_position': image position in the y-direction, 'magnifications': magnifications, 'time_delays': time delays, 'n_images': number of images formed, 'determinant': determinants, 'trace': traces, 'iteration': to keep track of the iteration number, 'weights': weights for the caustic considered]













      ..
          !! processed by numpydoc !!

   .. py:method:: get_lensed_snrs(snr_calculator, lensed_param, n_max_images=4)

      
      Function to calculate the signal to noise ratio for each image in each event.


      :Parameters:

          **snr_calculator** : `class`
              snr_calculator class
              this is an already initialized class that contains a function (snr_calculator.snr) that actually calculates snr with the given gw_params.

              Luminosity distance and time delay are modified to be effective luminosity distance and effective time delay, respectively, for each image using the magnifications and time delays.

          **lensed_param** : `dict`
              dictionary containing the both already lensed source paramters and image parameters.
              e.g. lensed_param.keys() = ['mass_1', 'mass_2', 'zs', 'luminosity_distance', 'iota', 'psi', 'phi', 'ra', 'dec', 'geocent_time', 'phase', 'a_1', 'a2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl', 'magnifications', 'time_delays']

          **n_max_images** : `int`
              maximum number of images to consider
              default: 4

      :Returns:

          **snrs** : `dict`
              signal to noise ratio for each image in each event.
              (dictionary containing 'H1', 'L1', ..., and 'opt_snr_net', which is the network snr, for each image as an array with dimensions (number_of_lensed_events,n_max_images) )













      ..
          !! processed by numpydoc !!


.. py:function:: solve_lens_equation1(lens_parameters)

   
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
           x position of the source in the source plane

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

       **weights** : `float`
           weights for the caustic










   .. rubric:: Examples

   >>> from ler.multiprocessing_routine import solve_lens_equation1
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
   >>> # each output parameter contains x_source, y_source, x0_image_position, x1_image_position, magnifications, time_delays, nImages, determinant, trace, iteration, weights
   >>> print(f"magnification of images with lens parameters 'lens_parameters1' is {result[0][6]}")
   magnification of images with lens parameters 'lens_parameters1' is [ 2.18973765 -1.27542831]



   ..
       !! processed by numpydoc !!

.. py:function:: solve_lens_equation2(lens_parameters)

   
   Function to solve the lens equation (min_image > 2)


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
           x position of the source in the source plane

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

       **weights** : `float`
           weights for the caustic










   .. rubric:: Examples

   >>> from ler.multiprocessing_routine import solve_lens_equation2
   >>> import numpy as np
   >>> from multiprocessing import Pool
   >>> # lens parameters input contains 12 parameters [e1, e2, gamma, gamma1, gamma2, zl, zs, einstein_radius, iteration, lens_model_list]
   >>> lens_parameters1 = np.array([3, 0.024069457093642648, -0.016002190961948142, 1.8945414936459974, 0.10117465203892329, 0.09600089396968613, 0.2503743800068136, 0.9418211055453296, 2.5055790287104725e-06, 0, 'EPL_NUMBA', 'SHEAR'], dtype=object)
   >>> lens_parameters2 = np.array([3, -0.04030088581646998, -0.01419438113690042, 2.0068239327017, 0.08482718989370612, -0.015393332086560785, 1.0952303138971118, 2.5534097159384417, 1.0125570159563301e-06, 1, 'EPL_NUMBA', 'SHEAR'], dtype=object)
   >>> input_arguments = np.vstack((lens_parameters1, lens_parameters2))
   >>> # solve the lens equation for each set of lens parameters
   >>> with Pool(2) as p:
   ...     result = p.map(solve_lens_equation2, input_arguments)
   >>> # result is a list of tuples
   >>> # each tuple contains the output parameters of the function
   >>> # each output parameter contains x_source, y_source, x0_image_position, x1_image_position, magnifications, time_delays, nImages, determinant, trace, iteration, weights
   >>> print(f"magnification of images with lens parameters 'lens_parameters1' is {result[0][6]}")
   magnification of images with lens parameters 'lens_parameters1' is [ 2.18973765 -1.27542831]



   ..
       !! processed by numpydoc !!

