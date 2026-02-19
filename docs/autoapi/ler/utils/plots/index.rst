:py:mod:`ler.utils.plots`
=========================

.. py:module:: ler.utils.plots

.. autoapi-nested-parse::

   This module contains functions to plot the results.

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.utils.plots.param_plot
   ler.utils.plots.relative_mu_dt_unlensed
   ler.utils.plots.relative_mu_dt_lensed
   ler.utils.plots.mu_vs_dt_plot



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
   >>> from ler.utils import param_plot
   >>> from ler.rates import LeR
   >>> ler = LeR(verbose=False)
   >>> param = ler.unlensed_cbc_statistics();
   >>> rate, param_detectable = ler.unlensed_rate()
   >>> plt.figure(figsize=(6, 4))
   >>> param_plot(param_name='zs', param_dict=param, plot_label='all events')
   >>> param_plot(param_name='zs', param_dict=param_detectable, plot_label='detectable events')
   >>> plt.xlabel('source redshift')
   >>> plt.ylabel('probability density')
   >>> plt.title('source redshift distribution')
   >>> plt.grid(alpha=0.5)
   >>> plt.savefig('source_redshift_distribution.png')



   ..
       !! processed by numpydoc !!

.. py:function:: relative_mu_dt_unlensed(param, size=100, randomize=True)

   
   Function to generate relative magnification vs time delay difference for unlensed samples.


   :Parameters:

       **param** : `dict`
           dictionary of unlensed GW source parameters.
           unlensed_param.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time']

       **size** : `int`
           number of samples.
           default size = 100.

       **randomize** : `bool`
           if True, it will randomize the samples.
           default randomize = True.

   :Returns:

       **dmu** : `float.array`
           relative magnification: abs(mu2/mu1) or abs(dl1/dl2)**2.

       **dt** : `float.array`
           relative time delay: abs(t1-t2) in days.













   ..
       !! processed by numpydoc !!

.. py:function:: relative_mu_dt_lensed(lensed_param, pdet_threshold=[0.5, 0.5], classification_type='morse_phase')

   
   Function to classify the lensed images wrt to the morse phase difference.


   :Parameters:

       **lensed_param** : `dict`
           dictionary of lensed GW source parameters.
           lensed_param.keys() = ['x_source', 'y_source', 'x0_image_position', 'x1_image_position', 'magnifications', 'time_delays', 'pdet_net', 'image_type']

       **pdet_threshold** : `list`
           threshold for pdet_net to consider the image as detectable.
           default pdet_threshold = [0.5, 0.5].

       **classification_type** : `str`
           type of classification to be done.
           default classification_type = 'morse_phase'. other options: 'arrival_time'.

   :Returns:

       dict
           dictionary containing the relative magnification and time delay difference for the classified images.
            if classification_type = 'morse_phase':
               {
                   "dt_rel0": np.array of relative time delay difference for images with morse phase difference = 0,
                   "mu_rel0": np.array of relative magnification for images with morse phase difference = 0,
                   "dt_rel90": np.array of relative time delay difference for images with morse phase difference = 90,
                   "mu_rel90": np.array of relative magnification for images with morse phase difference = 90,
               }
            if classification_type = 'arrival_time':
               {
                   "dt_12": np.array of relative time delay difference for image 1 and image 2,
                   "mu_12": np.array of relative magnification for image 1 and image 2,
                   "dt_13": np.array of relative time delay difference for image 1 and image 3,
                   "mu_13": np.array of relative magnification for image 1 and image 3,
                   "dt_14": np.array of relative time delay difference for image 1 and image 4,
                   "mu_14": np.array of relative magnification for image 1 and image 4,
                   "dt_23": np.array of relative time delay difference for image 2 and image 3,
                   "mu_23": np.array of relative magnification for image 2 and image 3,
                   "dt_24": np.array of relative time delay difference for image 2 and image 4,
                   "mu_24": np.array of relative magnification for image 2 and image 4,
                   "dt_34": np.array of relative time delay difference for image 3 and image 4,
                   "mu_34": np.array of relative magnification for image 3 and image 4,
               }













   ..
       !! processed by numpydoc !!

.. py:function:: mu_vs_dt_plot(x_array, y_array, xscale='log10', yscale='log10', alpha=0.6, extent=None, contour_levels=[10, 40, 68, 95], colors=['blue', 'blue', 'blue', 'blue', 'blue'])

   
   Function to generate 2D KDE and plot the relative magnification vs time delay difference for lensed samples.


   :Parameters:

       **x_array** : `float.array`
           x array.

       **y_array** : `float.array`
           y array.

       **xscale** : `str`
           x-axis scale.
           default xscale = 'log10'. other options: 'log', None.

       **yscale** : `str`
           y-axis scale.
           default yscale = 'log10'. other options: 'log', None.

       **alpha** : `float`
           transparency of the contour plot.
           default alpha = 0.6.

       **extent** : `list`
           extent of the plot.
           default extent = None. It will consider the full range of x_array and y_array.

       **contour_levels** : `list`
           levels for contour plot.
           default contour_levels = [10, 40, 68, 95].

       **colors** : `str`
           colors for contour plot.
           default colors = ['blue', 'blue', 'blue', 'blue', 'blue'].











   .. rubric:: Examples

   >>> import numpy as np
   >>> import matplotlib.pyplot as plt
   >>> from ler.utils import param_plot, mu_vs_dt_plot, get_param_from_json, relative_mu_dt_unlensed, relative_mu_dt_lensed
   >>> # get the parameters. For data generation, refer to the 'LeR complete example' in the documentation.
   >>> unlensed_param = get_param_from_json('ler_data/unlensed_param.json')
   >>> unlensed_param_detectable = get_param_from_json('ler_data/unlensed_param_detectable.json')
   >>> lensed_param = get_param_from_json('ler_data/lensed_param.json')
   >>> lensed_param_detectable = get_param_from_json('ler_data/lensed_param_detectable.json')
   >>> # get the relative mu and dt
   >>> ans = relative_mu_dt_lensed(lensed_param_detectable)
   >>> dmu, dt = relative_mu_dt_unlensed(unlensed_param_detectable, size=1000, randomize=True)
   >>> # plot
   >>> plt.figure(figsize=(4, 4))
   >>> mu_vs_dt_plot(ans['dt_rel90'], ans['mu_rel90'], colors='b')
   >>> mu_vs_dt_plot(ans['dt_rel0'], ans['mu_rel0'], colors='g')
   >>> mu_vs_dt_plot(dt, dmu, colors='r')
   >>> # Create proxy artists for legend
   >>> proxy1 = plt.Line2D([0], [0], linestyle='-', color='b', label=r'Lensed ($\Delta \phi=90$)')
   >>> proxy2 = plt.Line2D([0], [0], linestyle='-', color='g', label=r'Lensed ($\Delta \phi=0$)')
   >>> proxy3 = plt.Line2D([0], [0], linestyle='-', color='r', label=r'Unlensed')
   >>> plt.legend(handles=[proxy1, proxy2, proxy3], loc='upper left')
   >>> plt.xlim(-5, 2.5)
   >>> plt.ylim(-2.5, 2.5)
   >>> plt.grid(alpha=0.4)
   >>> plt.show()



   ..
       !! processed by numpydoc !!

