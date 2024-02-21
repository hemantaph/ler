:orphan:

:py:mod:`ler.utils`
===================

.. py:module:: ler.utils


Submodules
----------
.. toctree::
   :titlesonly:
   :maxdepth: 1

   plots/index.rst
   utils/index.rst


Package Contents
----------------

Classes
~~~~~~~

.. autoapisummary::

   ler.utils.NumpyEncoder



Functions
~~~~~~~~~

.. autoapisummary::

   ler.utils.load_json
   ler.utils.append_json
   ler.utils.add_dict_values
   ler.utils.get_param_from_json
   ler.utils.rejection_sample
   ler.utils.rejection_sample2d
   ler.utils.add_dictionaries_together
   ler.utils.trim_dictionary
   ler.utils.create_func_pdf_invcdf
   ler.utils.create_conditioned_pdf_invcdf
   ler.utils.create_func
   ler.utils.create_func_inv
   ler.utils.create_pdf
   ler.utils.create_inv_cdf_array
   ler.utils.create_conditioned_pdf
   ler.utils.create_conditioned_inv_cdf_array
   ler.utils.interpolator_from_pickle
   ler.utils.interpolator_pickle_path
   ler.utils.interpolator_pdf_conditioned
   ler.utils.interpolator_sampler_conditioned
   ler.utils.cubic_spline_interpolator
   ler.utils.inverse_transform_sampler
   ler.utils.batch_handler
   ler.utils.get_param_from_json
   ler.utils.param_plot
   ler.utils.relative_mu_dt_unlensed
   ler.utils.relative_mu_dt_lensed
   ler.utils.mu_vs_dt_plot



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

.. py:function:: append_json(file_name, new_dictionary, old_dictionary=None, replace=False)

   
   Append and update a json file with a dictionary.


   :Parameters:

       **file_name** : `str`
           json file name for storing the parameters.

       **new_dictionary** : `dict`
           dictionary to be appended to the json file.

       **replace** : `bool`, optional
           If True, replace the json file with the dictionary. Default is False.














   ..
       !! processed by numpydoc !!

.. py:function:: add_dict_values(dict1, dict2)

   
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
           x values.

       **y** : `numpy.ndarray`
           y values.

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

.. py:function:: interpolator_from_pickle(param_dict_given, directory, sub_directory, name, x, pdf_func=None, y=None, conditioned_y=None, dimension=1, category='pdf', create_new=False)

   
   Function to decide which interpolator to use.


   :Parameters:

       **param_dict_given** : `dict`
           dictionary of parameters.

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

.. py:function:: interpolator_pickle_path(param_dict_given, directory, sub_directory, interpolator_name)

   
   Function to create the interpolator pickle file path.


   :Parameters:

       **param_dict_given** : `dict`
           dictionary of parameters.

       **directory** : `str`
           directory to store the interpolator.

       **sub_directory** : `str`
           sub-directory to store the interpolator.

       **interpolator_name** : `str`
           name of the interpolator.

   :Returns:

       **path_inv_cdf** : `str`
           path of the interpolator pickle file.

       **it_exist** : `bool`
           if True, the interpolator exists.













   ..
       !! processed by numpydoc !!

.. py:function:: interpolator_pdf_conditioned(x, conditioned_y, y_array, interpolator_list)

   
   Function to find the pdf interpolator coefficients from the conditioned y.


   :Parameters:

       **x** : `numpy.ndarray`
           x values.

       **conditioned_y** : `float`
           conditioned y value.

       **y_array** : `numpy.ndarray`
           y values.

       **interpolator_list** : `list`
           list of interpolators.

   :Returns:

       **interpolator_list[idx](x)** : `numpy.ndarray`
           samples from the interpolator.













   ..
       !! processed by numpydoc !!

.. py:function:: interpolator_sampler_conditioned(conditioned_y, y_array, interpolator_list, size=1000)

   
   Function to find sampler interpolator coefficients from the conditioned y.


   :Parameters:

       **conditioned_y** : `float`
           conditioned y value.

       **y_array** : `numpy.ndarray`
           y values.

       **interpolator_list** : `list`
           list of interpolators.

       **size** : `int`
           number of samples.

   :Returns:


           ..













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

.. py:function:: batch_handler(size, batch_size, sampling_routine, output_jsonfile, save_batch=True, resume=False)

   
   Function to run the sampling in batches.


   :Parameters:

       **size** : `int`
           number of samples.

       **batch_size** : `int`
           batch size.

       **sampling_routine** : `function`
           function to sample the parameters.
           e.g. unlensed_sampling_routine() or lensed_sampling_routine()

       **output_jsonfile** : `str`
           name of the json file to store the parameters.

       **resume** : `bool`
           if True, it will resume the sampling from the last batch.
           default resume = False.














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

.. py:function:: param_plot(param_name='zs', param_dict='./gw_params.json', plot_label='zs', param_min=None, param_max=None, kde=True, kde_bandwidth=0.2, histogram=True, histogram_bins=30)

   
   Function to plot the distribution of the GW source parameters.


   :Parameters:

       **param_name** : `str`
           name of the parameter to plot.
           default param_name = 'zs'.

       **param_dict** : `dict` or `str`
           dictionary of GW source parameters or json file name.
           default param_dict = './gw_params.json'.

       **param_xlabel** : `str`
           x-axis label.
           default param_xlabel = 'source redshift'.

       **param_ylabel** : `str`
           y-axis label.
           default param_ylabel = 'probability density'.

       **param_min** : `float`
           minimum value of the parameter.
           default param_min = None.

       **param_max** : `float`
           maximum value of the parameter.
           default param_max = None.

       **figsize** : `tuple`
           figure size.
           default figsize = (4, 4).

       **kde** : `bool`
           if True, kde will be plotted.
           default kde = True.

       **kde_bandwidth** : `float`
           bandwidth for kde.
           default kde_bandwidth = 0.2.

       **histogram** : `bool`
           if True, histogram will be plotted.
           default histogram = True.

       **histogram_bins** : `int`
           number of bins for histogram.
           default histogram_bins = 30.











   .. rubric:: Examples

   >>> import matplotlib.pyplot as plt
   >>> from ler.rates import LeR
   >>> ler = LeR()
   >>> param = ler.unlensed_cbc_statistics();
   >>> rate, param_detectable = ler.unlensed_rate()
   >>> plt.figure(figsize=(6, 4))
   >>> ler.param_plot(param_name='zs', param_dict=param, plot_label='all events')
   >>> ler.param_plot(param_name='zs', param_dict=param_detectable, plot_label='detectable events')
   >>> plt.xlabel('source redshift')
   >>> plt.ylabel('probability density')
   >>> plt.title('source redshift distribution')
   >>> plt.grid(alpha=0.5)
   >>> plt.savefig('source_redshift_distribution.png')



   ..
       !! processed by numpydoc !!

.. py:function:: relative_mu_dt_unlensed(param, size=100)

   
   Function to generate relative magnification vs time delay difference for unlensed samples.


   :Parameters:

       **param** : `dict`
           dictionary of unlensed GW source parameters.
           unlensed_param.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time']

   :Returns:

       **dmu** : `float.array`
           relative magnification.

       **dt** : `float.array`
           relative time delay.













   ..
       !! processed by numpydoc !!

.. py:function:: relative_mu_dt_lensed(lensed_param, snr_threshold=[8.0, 8.0])

   
   Function to classify the lensed images wrt to the morse phase difference.


   :Parameters:

       **lensed_param** : `dict`
           dictionary of lensed GW source parameters, lens galaxy parameters and image paramters.
           lensed_param.keys() = ['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2', 'Dl',
           'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source',
           'luminosity_distance', 'theta_jn', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'n_images',
           'x0_image_positions', 'x1_image_positions', 'magnifications', 'time_delays', 'traces',
           'determinants', 'image_type', 'weights', 'optimal_snr_net', 'L1', 'H1', 'V1']

       **snr_threshold** : `float`
           threshold for detection signal to noise ratio.
           e.g. snr_threshold = [8.,8.] or [8.,6.] for subthreshold

   :Returns:

       **mu_rel0** : `float.array`
           relative magnification for 0 degree phase difference.

       **dt_rel0** : `float.array`
           relative time delay for 0 degree phase difference.

       **mu_rel90** : `float.array`
           relative magnification for 90 degree phase difference.

       **dt_rel90** : `float.array`
           relative time delay for 90 degree phase difference.













   ..
       !! processed by numpydoc !!

.. py:function:: mu_vs_dt_plot(x_array, y_array, savefig=False, ax=None, colors='blue', linestyles='-', origin='upper', alpha=0.6, extent=[0.01, 500.0, 0.01, 100.0], contour_levels=[0.1, 0.4, 0.68, 0.95])

   
   Function to generate 2D KDE and plot the relative magnification vs time delay difference for lensed samples.


   :Parameters:

       **x_array** : `float.array`
           x array.

       **y_array** : `float.array`
           y array.

       **xlabel** : `str`
           x label.

       **ylabel** : `str`
           y label.

       **title** : `str`
           title.

       **savefig** : `bool`
           if True, it will save the figure.
           default savefig = False.

       **ax** : `matplotlib.axes`
           matplotlib axes.
           default ax = None.

       **colors** : `str`
           color of the plot.
           default colors = 'blue'.

       **linestyles** : `str`
           linestyle of the plot.
           default linestyles = '-'.

       **origin** : `str`
           origin of the plot.
           default origin = 'upper'.

       **alpha** : `float`
           alpha of the plot.
           default alpha = 0.6.

       **extent** : `list`
           extent of the plot.
           default extent = [1e-2,5e2,1e-2,1e2].

       **contour_levels** : `list`
           contour levels of the plot.
           default contour_levels = [0.10,0.40,0.68,0.95] which corresponds to 1,2,3,4 sigma.

   :Returns:

       None
           ..













   ..
       !! processed by numpydoc !!

