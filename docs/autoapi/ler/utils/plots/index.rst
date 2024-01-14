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

