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

.. py:function:: relative_mu_dt_lensed(lensed_param, snr_threshold=[8.0, 8.0], classification_type='morse_phase')

   
   Function to classify the lensed images wrt to the morse phase difference.
















   ..
       !! processed by numpydoc !!

.. py:function:: mu_vs_dt_plot(x_array, y_array, xscale='log10', yscale='log10', savefig=False, linestyles='-', origin='upper', alpha=0.6, extent=None, contour_levels=[10, 40, 68, 95], colors=['blue', 'blue', 'blue', 'blue', 'blue'])

   
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
   >>> plt.legend(handles=[proxy1, proxy2, proxy3])
   >>> plt.xlim(-5, 2.5)
   >>> plt.ylim(-2.5, 2.5)
   >>> plt.grid(alpha=0.4)
   >>> plt.show()



   ..
       !! processed by numpydoc !!

