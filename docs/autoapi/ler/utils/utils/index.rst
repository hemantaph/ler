:py:mod:`ler.utils.utils`
=========================

.. py:module:: ler.utils.utils

.. autoapi-nested-parse::

   This module contains helper routines for other modules in the ler package.

   ..
       !! processed by numpydoc !!


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   ler.utils.utils.NumpyEncoder



Functions
~~~~~~~~~

.. autoapisummary::

   ler.utils.utils.is_njitted
   ler.utils.utils.load_pickle
   ler.utils.utils.save_pickle
   ler.utils.utils.load_hdf5
   ler.utils.utils.save_hdf5
   ler.utils.utils.load_json
   ler.utils.utils.save_json
   ler.utils.utils.append_json
   ler.utils.utils.concatenate_dict_values
   ler.utils.utils.get_param_from_json
   ler.utils.utils.load_txt_from_module
   ler.utils.utils.rejection_sample
   ler.utils.utils.rejection_sample2d
   ler.utils.utils.add_dictionaries_together
   ler.utils.utils.trim_dictionary
   ler.utils.utils.create_func_pdf_invcdf
   ler.utils.utils.create_conditioned_pdf_invcdf
   ler.utils.utils.create_func
   ler.utils.utils.create_func_inv
   ler.utils.utils.create_pdf
   ler.utils.utils.create_inv_cdf_array
   ler.utils.utils.create_conditioned_pdf
   ler.utils.utils.create_conditioned_inv_cdf_array
   ler.utils.utils.interpolator_from_json
   ler.utils.utils.interpolator_json_path
   ler.utils.utils.batch_handler
   ler.utils.utils.create_batch_params
   ler.utils.utils.monte_carlo_integration
   ler.utils.utils.cubic_spline_interpolator
   ler.utils.utils.pdf_cubic_spline_interpolator
   ler.utils.utils.pdf_cubic_spline_interpolator2d_array
   ler.utils.utils.cubic_spline_interpolator2d_array
   ler.utils.utils.coefficients_generator_ler
   ler.utils.utils.inverse_transform_sampler2d
   ler.utils.utils.inverse_transform_sampler
   ler.utils.utils.normal_pdf
   ler.utils.utils.normal_pdf_2d
   ler.utils.utils.load_txt_from_module
   ler.utils.utils.cumulative_trapezoid
   ler.utils.utils.sample_from_powerlaw_distribution



.. py:function:: is_njitted(func)


.. py:class:: NumpyEncoder(*, skipkeys=False, ensure_ascii=True, check_circular=True, allow_nan=True, sort_keys=False, indent=None, separators=None, default=None)


   Bases: :py:obj:`json.JSONEncoder`

   
   Class for storing a numpy.ndarray or any nested-list composition as JSON file. This is required for dealing np.nan and np.inf.


   :Parameters:

       **json.JSONEncoder** : `class`
           class for encoding JSON file

   :Returns:

       **json.JSONEncoder.default** : `function`
           function for encoding JSON file













   ..
       !! processed by numpydoc !!
   .. py:method:: default(obj)

      
      function for encoding JSON file
















      ..
          !! processed by numpydoc !!


.. py:function:: load_pickle(file_name)

   
   Load a pickle file.


   :Parameters:

       **file_name** : `str`
           pickle file name for storing the parameters.

   :Returns:

       **param** : `dict`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: save_pickle(file_name, param)

   
   Save a dictionary as a pickle file.


   :Parameters:

       **file_name** : `str`
           pickle file name for storing the parameters.

       **param** : `dict`
           dictionary to be saved as a pickle file.














   ..
       !! processed by numpydoc !!

.. py:function:: load_hdf5(file_name)

   
   Load a hdf5 file.


   :Parameters:

       **file_name** : `str`
           hdf5 file name for storing the parameters.

   :Returns:

       **param** : `dict`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: save_hdf5(file_name, param)

   
   Save a dictionary as a hdf5 file.


   :Parameters:

       **file_name** : `str`
           hdf5 file name for storing the parameters.

       **param** : `dict`
           dictionary to be saved as a hdf5 file.














   ..
       !! processed by numpydoc !!

.. py:function:: load_json(file_name)

   
   Load a json file.


   :Parameters:

       **file_name** : `str`
           json file name for storing the parameters.

   :Returns:

       **param** : `dict`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: save_json(file_name, param)

   
   Save a dictionary as a json file.


   :Parameters:

       **file_name** : `str`
           json file name for storing the parameters.

       **param** : `dict`
           dictionary to be saved as a json file.














   ..
       !! processed by numpydoc !!

.. py:function:: append_json(file_name, new_dictionary, old_dictionary=None, replace=False)

   
   Append (values with corresponding keys) and update a json file with a dictionary. There are four options:

   1. If old_dictionary is provided, the values of the new dictionary will be appended to the old dictionary and save in the 'file_name' json file.
   2. If replace is True, replace the json file (with the 'file_name') content with the new_dictionary.
   3. If the file does not exist, create a new one with the new_dictionary.
   4. If none of the above, append the new dictionary to the content of the json file.

   :Parameters:

       **file_name** : `str`
           json file name for storing the parameters.

       **new_dictionary** : `dict`
           dictionary to be appended to the json file.

       **old_dictionary** : `dict`, optional
           If provided the values of the new dictionary will be appended to the old dictionary and save in the 'file_name' json file.
           Default is None.

       **replace** : `bool`, optional
           If True, replace the json file with the dictionary. Default is False.














   ..
       !! processed by numpydoc !!

.. py:function:: concatenate_dict_values(dict1, dict2)

   
   Adds the values of two dictionaries together.


   :Parameters:

       **dict1** : `dict`
           dictionary to be added.

       **dict2** : `dict`
           dictionary to be added.

   :Returns:

       **dict1** : `dict`
           dictionary with added values.













   ..
       !! processed by numpydoc !!

.. py:function:: get_param_from_json(json_file)

   
   Function to get the parameters from json file.


   :Parameters:

       **json_file** : `str`
           json file name for storing the parameters.

   :Returns:

       **param** : `dict`
           ..













   ..
       !! processed by numpydoc !!

.. py:function:: load_txt_from_module(package, directory, filename)

   
   Function to load a specific dataset from a .txt file within the package


   :Parameters:

       **package** : str
           name of the package

       **directory** : str
           name of the directory within the package

       **filename** : str
           name of the .txt file

   :Returns:

       **data** : `dict`
           Dictionary loaded from the .txt file













   ..
       !! processed by numpydoc !!

.. py:function:: rejection_sample(pdf, xmin, xmax, size=100, chunk_size=10000)

   
   Helper function for rejection sampling from a pdf with maximum and minimum arguments.


   :Parameters:

       **pdf** : `function`
           pdf function.

       **xmin** : `float`
           minimum value of the pdf.

       **xmax** : `float`
           maximum value of the pdf.

       **size** : `int`, optional
           number of samples. Default is 100.

       **chunk_size** : `int`, optional
           chunk size for sampling. Default is 10000.

   :Returns:

       **x_sample** : `numpy.ndarray`
           samples from the pdf.













   ..
       !! processed by numpydoc !!

.. py:function:: rejection_sample2d(pdf, xmin, xmax, ymin, ymax, size=100, chunk_size=10000)

   
   Helper function for rejection sampling from a 2D pdf with maximum and minimum arguments.


   :Parameters:

       **pdf** : `function`
           2D pdf function.

       **xmin** : `float`
           minimum value of the pdf in the x-axis.

       **xmax** : `float`
           maximum value of the pdf in the x-axis.

       **ymin** : `float`
           minimum value of the pdf in the y-axis.

       **ymax** : `float`
           maximum value of the pdf in the y-axis.

       **size** : `int`, optional
           number of samples. Default is 100.

       **chunk_size** : `int`, optional
           chunk size for sampling. Default is 10000.

   :Returns:

       **x_sample** : `numpy.ndarray`
           samples from the pdf in the x-axis.













   ..
       !! processed by numpydoc !!

.. py:function:: add_dictionaries_together(dictionary1, dictionary2)

   
   Adds two dictionaries with the same keys together.


   :Parameters:

       **dictionary1** : `dict`
           dictionary to be added.

       **dictionary2** : `dict`
           dictionary to be added.

   :Returns:

       **dictionary** : `dict`
           dictionary with added values.













   ..
       !! processed by numpydoc !!

.. py:function:: trim_dictionary(dictionary, size)

   
   Filters an event dictionary to only contain the size.


   :Parameters:

       **dictionary** : `dict`
           dictionary to be trimmed.

       **size** : `int`
           size to trim the dictionary to.

   :Returns:

       **dictionary** : `dict`
           trimmed dictionary.













   ..
       !! processed by numpydoc !!

.. py:function:: create_func_pdf_invcdf(x, y, category='function')

   
   Function to create a interpolated function, inverse function or inverse cdf from the input x and y.


   :Parameters:

       **x** : `numpy.ndarray`
           x values. This has to sorted in ascending order.

       **y** : `numpy.ndarray`
           y values. Corresponding to the x values.

       **category** : `str`, optional
           category of the function. Default is "function". Other options are "function_inverse", "pdf" and "inv_cdf".

   :Returns:

       **pdf** : `pdf function`
           interpolated pdf function.

       **inv_pdf** : `function inverse`
           interpolated inverse pdf function.

       **inv_cdf** : `function`
           interpolated inverse cdf.













   ..
       !! processed by numpydoc !!

.. py:function:: create_conditioned_pdf_invcdf(x, conditioned_y, pdf_func, category)

   
   pdf_func is the function to calculate the pdf of x given y
   x is an array and the output of pdf_func is an array
   y is the condition
   we consider parameter plane of x and y


   :Parameters:

       **x** : `numpy.ndarray`
           x values.

       **conditioned_y** : `numpy.ndarray`
           conditioned y values.

       **pdf_func** : `function`
           function to calculate the pdf of x given y.

       **category** : `str`, optional
           category of the function. Default is "function". Other options are "function_inverse", "pdf" and "inv_cdf".














   ..
       !! processed by numpydoc !!

.. py:function:: create_func(x, y)

   
   Function to create a spline interpolated function from the input x and y.


   :Parameters:

       **x** : `numpy.ndarray`
           x values.

       **y** : `numpy.ndarray`
           y values.

   :Returns:

       **c** : `numpy.ndarray`
           spline coefficients.













   ..
       !! processed by numpydoc !!

.. py:function:: create_func_inv(x, y)

   
   Function to create a spline interpolated inverse function from the input x and y.


   :Parameters:

       **x** : `numpy.ndarray`
           x values.

       **y** : `numpy.ndarray`
           y values.

   :Returns:

       **c** : `numpy.ndarray`
           spline coefficients.













   ..
       !! processed by numpydoc !!

.. py:function:: create_pdf(x, y)

   
   Function to create a spline interpolated normalized pdf from the input x and y.


   :Parameters:

       **x** : `numpy.ndarray`
           x values.

       **y** : `numpy.ndarray`
           y values.

   :Returns:

       **c** : `numpy.ndarray`
           spline coefficients.













   ..
       !! processed by numpydoc !!

.. py:function:: create_inv_cdf_array(x, y)

   
   Function to create a spline interpolated inverse cdf from the input x and y.


   :Parameters:

       **x** : `numpy.ndarray`
           x values.

       **y** : `numpy.ndarray`
           y values.

   :Returns:

       **c** : `numpy.ndarray`
           spline coefficients.













   ..
       !! processed by numpydoc !!

.. py:function:: create_conditioned_pdf(x, conditioned_y, pdf_func)

   
   Function to create a conditioned pdf from the input x and y.


   :Parameters:

       **x** : `numpy.ndarray`
           x values.

       **conditioned_y** : `numpy.ndarray`
           conditioned y values.

       **pdf_func** : `function`
           function to calculate the pdf of x given y.

   :Returns:

       **list_** : `list`
           list of pdfs.













   ..
       !! processed by numpydoc !!

.. py:function:: create_conditioned_inv_cdf_array(x, conditioned_y, pdf_func)

   
   Function to create a conditioned inv_cdf from the input x and y.


   :Parameters:

       **x** : `numpy.ndarray`
           x values.

       **conditioned_y** : `numpy.ndarray`
           conditioned y values.

       **pdf_func** : `function`
           function to calculate the pdf of x given y.

   :Returns:

       **list_** : `list`
           list of inv_cdfs.













   ..
       !! processed by numpydoc !!

.. py:function:: interpolator_from_json(identifier_dict, directory, sub_directory, name, x, pdf_func=None, y=None, conditioned_y=None, dimension=1, category='pdf', create_new=False)

   
   Function to decide which interpolator to use.


   :Parameters:

       **identifier_dict** : `dict`
           dictionary of identifiers.

       **directory** : `str`
           directory to store the interpolator.

       **sub_directory** : `str`
           sub-directory to store the interpolator.

       **name** : `str`
           name of the interpolator.

       **x** : `numpy.ndarray`
           x values.

       **pdf_func** : `function`
           function to calculate the pdf of x given y.

       **y** : `numpy.ndarray`
           y values.

       **conditioned_y** : `numpy.ndarray`
           conditioned y values.

       **dimension** : `int`
           dimension of the interpolator. Default is 1.

       **category** : `str`
           category of the function. Default is "pdf".

       **create_new** : `bool`
           if True, create a new interpolator. Default is False.

   :Returns:

       **interpolator** : `function`
           interpolator function.













   ..
       !! processed by numpydoc !!

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

.. py:function:: batch_handler(size, batch_size, sampling_routine, output_jsonfile, save_batch=True, resume=False, param_name='parameters')

   
   Function to run the sampling in batches.


   :Parameters:

       **size** : `int`
           number of samples.

       **batch_size** : `int`
           batch size.

       **sampling_routine** : `function`
           sampling function. It should have 'size' as input and return a dictionary.

       **output_jsonfile** : `str`
           json file name for storing the parameters.

       **save_batch** : `bool`, optional
           if True, save sampled parameters in each iteration. Default is True.

       **resume** : `bool`, optional
           if True, resume sampling from the last batch. Default is False.

       **param_name** : `str`, optional
           name of the parameter. Default is 'parameters'.

   :Returns:

       **dict_buffer** : `dict`
           dictionary of parameters.













   ..
       !! processed by numpydoc !!

.. py:function:: create_batch_params(sampling_routine, frac_batches, dict_buffer, save_batch, output_jsonfile, track_batches, resume=False)

   
   Helper function to batch_handler. It create batch parameters and store in a dictionary.


   :Parameters:

       **sampling_routine** : `function`
           sampling function. It should have 'size' as input and return a dictionary.

       **frac_batches** : `int`
           batch size.

       **dict_buffer** : `dict`
           dictionary of parameters.

       **save_batch** : `bool`
           if True, save sampled parameters in each iteration.

       **output_jsonfile** : `str`
           json file name for storing the parameters.

       **track_batches** : `int`
           track the number of batches.

       **resume** : `bool`, optional
           if True, resume sampling from the last batch. Default is False.

   :Returns:

       **track_batches** : `int`
           track the number of batches.













   ..
       !! processed by numpydoc !!

.. py:function:: monte_carlo_integration(function, uniform_prior, size=10000)

   
   Function to perform Monte Carlo integration.


   :Parameters:

       **function** : `function`
           function to be integrated.

       **prior** : `function`
           prior function.

   :Returns:

       **integral** : `float`
           integral value.













   ..
       !! processed by numpydoc !!

.. py:function:: cubic_spline_interpolator(xnew, coefficients, x)

   
   Function to interpolate using cubic spline.


   :Parameters:

       **xnew** : `numpy.ndarray`
           new x values.

       **coefficients** : `numpy.ndarray`
           coefficients of the cubic spline.

       **x** : `numpy.ndarray`
           x values.

   :Returns:

       **result** : `numpy.ndarray`
           interpolated values.













   ..
       !! processed by numpydoc !!

.. py:function:: pdf_cubic_spline_interpolator(xnew, norm_const, coefficients, x)

   
   Function to interpolate pdf using cubic spline.


   :Parameters:

       **xnew** : `numpy.ndarray`
           new x values.

       **coefficients** : `numpy.ndarray`
           coefficients of the cubic spline.

       **x** : `numpy.ndarray`
           x values.

   :Returns:

       **result** : `numpy.ndarray`
           interpolated values.













   ..
       !! processed by numpydoc !!

.. py:function:: pdf_cubic_spline_interpolator2d_array(xnew_array, ynew_array, norm_array, coefficients, x, y)

   
   Function to calculate the interpolated value of snr_partialscaled given the mass ratio (ynew) and total mass (xnew). This is based off 2D bicubic spline interpolation.


   :Parameters:

       **xnew_array** : `numpy.ndarray`
           Total mass of the binary.

       **ynew_array** : `numpy.ndarray`
           Mass ratio of the binary.

       **coefficients** : `numpy.ndarray`
           Array of coefficients for the cubic spline interpolation.

       **x** : `numpy.ndarray`
           Array of total mass values for the coefficients.

       **y** : `numpy.ndarray`
           Array of mass ratio values for the coefficients.

   :Returns:

       **result** : `float`
           Interpolated value of snr_partialscaled.













   ..
       !! processed by numpydoc !!

.. py:function:: cubic_spline_interpolator2d_array(xnew_array, ynew_array, coefficients, x, y)

   
   Function to calculate the interpolated value of snr_partialscaled given the mass ratio (ynew) and total mass (xnew). This is based off 2D bicubic spline interpolation.


   :Parameters:

       **xnew_array** : `numpy.ndarray`
           Total mass of the binary.

       **ynew_array** : `numpy.ndarray`
           Mass ratio of the binary.

       **coefficients** : `numpy.ndarray`
           Array of coefficients for the cubic spline interpolation.

       **x** : `numpy.ndarray`
           Array of total mass values for the coefficients.

       **y** : `numpy.ndarray`
           Array of mass ratio values for the coefficients.

   :Returns:

       **result** : `float`
           Interpolated value of snr_partialscaled.













   ..
       !! processed by numpydoc !!

.. py:function:: coefficients_generator_ler(y1, y2, y3, y4, z1, z2, z3, z4)

   
   Function to generate the coefficients for the cubic spline interpolation of fn(y)=z.


   :Parameters:

       **y1, y2, y3, y4, z1, z2, z3, z4: `float`**
           Values of y and z for the cubic spline interpolation.

   :Returns:

       coefficients: `numpy.ndarray`
           Coefficients for the cubic spline interpolation.













   ..
       !! processed by numpydoc !!

.. py:function:: inverse_transform_sampler2d(size, conditioned_y, cdf2d, x2d, y1d)

   
   Function to find sampler interpolator coefficients from the conditioned y.


   :Parameters:

       **size: `int`**
           Size of the sample.

       **conditioned_y: `float`**
           Conditioned y value.

       **cdf2d: `numpy.ndarray`**
           2D array of cdf values.

       **x2d: `numpy.ndarray`**
           2D array of x values.

       **y1d: `numpy.ndarray`**
           1D array of y values.

   :Returns:

       samples: `numpy.ndarray`
           Samples of the conditioned y.













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

.. py:function:: normal_pdf(x, mean=0.0, std=0.05)

   
   Calculate the value of a normal probability density function.


   :Parameters:

       **x** : `float` or `numpy.ndarray`
           The value(s) at which to evaluate the PDF.

       **mean** : `float`, optional
           The mean of the normal distribution. Default is 0.

       **std** : `float`, optional
           The standard deviation of the normal distribution. Default is 0.05.

   :Returns:

       **pdf** : `float` or `numpy.ndarray`
           The probability density function value(s) at x.













   ..
       !! processed by numpydoc !!

.. py:function:: normal_pdf_2d(x, y, mean_x=0.0, std_x=0.05, mean_y=0.0, std_y=0.05)

   
   Calculate the value of a 2D normal probability density function.


   :Parameters:

       **x** : `float`
           The x-coordinate for which the PDF is evaluated.

       **y** : `float`
           The y-coordinate for which the PDF is evaluated.

       **mean_x** : `float`, optional
           The mean of the normal distribution along the x-axis. Default is 0.

       **std_x** : `float`, optional
           The standard deviation of the normal distribution along the x-axis. Default is 0.05.

       **mean_y** : `float`, optional
           The mean of the normal distribution along the y-axis. Default is 0.

       **std_y** : `float`, optional
           The standard deviation of the normal distribution along the y-axis. Default is 0.05.

   :Returns:

       `float`
           The probability density function value at the given x and y coordinates.













   ..
       !! processed by numpydoc !!

.. py:function:: load_txt_from_module(package, directory, filename)

   
















   ..
       !! processed by numpydoc !!

.. py:function:: cumulative_trapezoid(y, x=None, dx=1.0, initial=0.0)

   
   Compute the cumulative integral of a function using the trapezoidal rule.
















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

