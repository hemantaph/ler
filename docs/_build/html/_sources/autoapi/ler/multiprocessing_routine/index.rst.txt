:orphan:

:py:mod:`ler.multiprocessing_routine`
=====================================

.. py:module:: ler.multiprocessing_routine


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.multiprocessing_routine.solve_lens_equation1
   ler.multiprocessing_routine.solve_lens_equation2



.. py:function:: solve_lens_equation1(lens_parameters)

   
   Function to solve the lens equation (min_image = 2)
   Input parameters:
       lens_parameters : a list of parameters
                       lens_parameters[0] = e1 : ellipticity
                       lens_parameters[1] = e2 : ellipticity
                       lens_parameters[2] = gamma : power-law index
                       lens_parameters[3] = gamma1 : shear
                       lens_parameters[4] = gamma2 : shear
                       lens_parameters[5] = zl : redshift of the lens
                       lens_parameters[6] = zs : redshift of the source
                       lens_parameters[7] = einstein_radius : Einstein radius
                       lens_parameters[8] = iteration : iteration number
                       lens_parameters[9:] = lens_model_list : list of lens models
   Output parameters:
       x_source : x position of the source in the source plane
       y_source : y position of the source in the source plane
       eta : polar coordinate of the source in the source plane
       phi : polar coordinate of the source in the source plane
       x0_image_position : x position of the images in the source plane
       x1_image_position : y position of the images in the source plane
       magnifications : magnification of the images
       time_delays : time-delay of the images
       nImages : number of images
       determinant : determinant of the hessian matrix
       trace : trace of the hessian matrix
       iteration : iteration number
       weights : weights for the caustic
















   ..
       !! processed by numpydoc !!

.. py:function:: solve_lens_equation2(lens_parameters)

   
   Function to solve the lens equation (min_image > 2)
   Input parameters:
       lens_parameters : a list of parameters
                       lens_parameters[0] = e1 : ellipticity
                       lens_parameters[1] = e2 : ellipticity
                       lens_parameters[2] = gamma : power-law index
                       lens_parameters[3] = gamma1 : shear
                       lens_parameters[4] = gamma2 : shear
                       lens_parameters[5] = zl : redshift of the lens
                       lens_parameters[6] = zs : redshift of the source
                       lens_parameters[7] = einstein_radius : Einstein radius
                       lens_parameters[8] = iteration : iteration number
                       lens_parameters[9:] = lens_model_list : list of lens models
   Output parameters:
       x_source : x position of the source in the source plane
       y_source : y position of the source in the source plane
       eta : polar coordinate of the source in the source plane
       phi : polar coordinate of the source in the source plane
       x0_image_position : x position of the images in the source plane
       x1_image_position : y position of the images in the source plane
       magnifications : magnification of the images
       time_delays : time-delay of the images
       nImages : number of images
       determinant : determinant of the hessian matrix
       trace : trace of the hessian matrix
       iteration : iteration number
       weights : weights for the caustic
















   ..
       !! processed by numpydoc !!

