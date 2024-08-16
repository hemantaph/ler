# -*- coding: utf-8 -*-
"""
This module contains functions to plot the results.
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from scipy.interpolate import interp1d

from ler.utils.utils import get_param_from_json

def param_plot(
        param_name="zs",
        param_dict="./gw_params.json",
        plot_label="zs",
        param_min=None,
        param_max=None,
        kde=True,
        kde_bandwidth=0.2,
        histogram=True,
        histogram_bins=30,
):
    """
    Function to plot the distribution of the GW source parameters.

    Parameters
    ----------
    param_name : `str`
        name of the parameter to plot.
        default param_name = 'zs'.
    param_dict : `dict` or `str`
        dictionary of GW source parameters or json file name.
        default param_dict = './gw_params.json'.
    param_xlabel : `str`
        x-axis label.
        default param_xlabel = 'source redshift'.
    param_ylabel : `str`
        y-axis label.
        default param_ylabel = 'probability density'.
    param_min : `float`
        minimum value of the parameter.
        default param_min = None.
    param_max : `float`
        maximum value of the parameter.
        default param_max = None.
    figsize : `tuple`
        figure size.
        default figsize = (4, 4).
    kde : `bool`
        if True, kde will be plotted.
        default kde = True.
    kde_bandwidth : `float`
        bandwidth for kde.
        default kde_bandwidth = 0.2.
    histogram : `bool`
        if True, histogram will be plotted.
        default histogram = True.
    histogram_bins : `int`
        number of bins for histogram.
        default histogram_bins = 30.

    Examples
    ----------
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
    """

    # get gw params from json file if not provided
    if type(param_dict) == str:
        print(f"getting gw_params from json file {param_dict}...")
        param_dict = get_param_from_json(param_dict)

    if param_min is None:
        param_min = np.min(param_dict[param_name])
    if param_max is None:
        param_max = np.max(param_dict[param_name])

    # plot the distribution of the parameter
    if histogram:
        plt.hist(
            param_dict[param_name],
            bins=histogram_bins,
            density=True,
            histtype="step",
            label=plot_label,
        )
    if kde:
        kde = gaussian_kde(param_dict[param_name], bw_method=kde_bandwidth)
        x = np.linspace(param_min, param_max, 1000)
        plt.plot(x, kde(x), label=plot_label+" kde")
    plt.legend()

def relative_mu_dt_unlensed(param, size=100, randomize=True):
    """
    Function to generate relative magnification vs time delay difference for unlensed samples.

    Parameters
    ----------
    param : `dict`
        dictionary of unlensed GW source parameters.
        unlensed_param.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time']
    size : `int`
        number of samples.
        default size = 100.
    randomize : `bool`
        if True, it will randomize the samples.
        default randomize = True.

    Returns
    ----------
    dmu : `float.array`
        relative magnification: abs(mu2/mu1) or abs(dl1/dl2)**2.
    dt : `float.array`
        relative time delay: abs(t1-t2) in days.
    """

    t = param["geocent_time"]
    mu = param["luminosity_distance"]
    len_ = len(t)
    # randomize it
    if randomize:
        idx_ = np.random.permutation(len_)
        t = t[idx_]
        mu = mu[idx_]
        
    # Ensure enough unique pairs can be formed
    if size > (len(t) * (len(t) - 1)) // 2:
        raise ValueError(f"size should be less than the number of unique pairs {len(t) * (len(t) - 1) // 2}")

    # Generate unique pairs
    # find idx1 and idx2
    idx1 = np.array([])
    idx2 = np.array([])
    while len(idx1) < size:
        idx1_ = np.random.choice(len_, size=size, replace=True)
        idx2_ = np.random.choice(len_, size=size, replace=True)
        idx1 = np.concatenate((idx1, idx1_))
        idx2 = np.concatenate((idx2, idx2_))
        idx = np.where(idx1 != idx2)[0]
        idx1 = idx1[idx]
        idx2 = idx2[idx]
    idx1 = idx1[:size].astype(int)
    idx2 = idx2[:size].astype(int)

    dt = abs(t[idx1] - t[idx2]) / (60 * 60 * 24)  # in days
    dmu = abs(mu[idx1]/mu[idx2])**2

    return dmu, dt

def relative_mu_dt_lensed(
    lensed_param, 
    snr_threshold=[8.0, 8.0], 
    classification_type='morse_phase'
    ):
    """
    Function to classify the lensed images wrt to the morse phase difference.

    
    """

    # get magnifications, time_delays and snr
    mu = lensed_param["magnifications"]
    dt = lensed_param["time_delays"]
    snr = lensed_param["optimal_snr_net"]
    image_type = lensed_param["image_type"]

    # pair images wrt to image_type
    if classification_type == 'morse_phase':
        dt_rel0 = []
        mu_rel0 = []
        dt_rel90 = []
        mu_rel90 = []
        for i in range(len(image_type)):
            if image_type[i,0]==image_type[i,1]:
                # snr check
                # below will also take care of the nan values
                if snr[i,0]>snr_threshold[0] and snr[i,1]>snr_threshold[1]:
                    dt_rel0.append(abs(dt[i,1]-dt[i,0])/ (60 * 60 * 24))
                    mu_rel0.append(abs(mu[i,1]/mu[i,0]))
            else:
                if snr[i,0]>snr_threshold[0] and snr[i,1]>snr_threshold[1]:
                    dt_rel90.append(abs(dt[i,1]-dt[i,0])/ (60 * 60 * 24))
                    mu_rel90.append(abs(mu[i,1]/mu[i,0]))
            if image_type[i,0]==image_type[i,2]:
                # snr check
                # below will also take care of the nan values
                if snr[i,0]>snr_threshold[0] and snr[i,2]>snr_threshold[1]:
                    dt_rel0.append(abs(dt[i,2]-dt[i,0])/ (60 * 60 * 24))
                    mu_rel0.append(abs(mu[i,2]/mu[i,0]))
            else:
                if snr[i,0]>snr_threshold[0] and snr[i,2]>snr_threshold[1]:
                    dt_rel90.append(abs(dt[i,2]-dt[i,0])/ (60 * 60 * 24))
                    mu_rel90.append(abs(mu[i,2]/mu[i,0]))
            if image_type[i,0]==image_type[i,3]:
                # snr check
                # below will also take care of the nan values
                if snr[i,0]>snr_threshold[0] and snr[i,3]>snr_threshold[1]:
                    dt_rel0.append(abs(dt[i,3]-dt[i,0])/ (60 * 60 * 24))
                    mu_rel0.append(abs(mu[i,3]/mu[i,0]))
            else:
                if snr[i,0]>snr_threshold[0] and snr[i,3]>snr_threshold[1]:
                    dt_rel90.append(abs(dt[i,3]-dt[i,0])/ (60 * 60 * 24))
                    mu_rel90.append(abs(mu[i,3]/mu[i,0]))
            if image_type[i,1]==image_type[i,2]:
                # snr check
                # below will also take care of the nan values
                if snr[i,1]>snr_threshold[0] and snr[i,2]>snr_threshold[1]:
                    dt_rel0.append(abs(dt[i,2]-dt[i,1])/ (60 * 60 * 24))
                    mu_rel0.append(abs(mu[i,2]/mu[i,1]))
            else:
                if snr[i,1]>snr_threshold[0] and snr[i,2]>snr_threshold[1]:
                    dt_rel90.append(abs(dt[i,2]-dt[i,1])/ (60 * 60 * 24))
                    mu_rel90.append(abs(mu[i,2]/mu[i,1]))
            if image_type[i,1]==image_type[i,3]:
                # snr check
                # below will also take care of the nan values
                if snr[i,1]>snr_threshold[0] and snr[i,3]>snr_threshold[1]:
                    dt_rel0.append(abs(dt[i,3]-dt[i,1])/ (60 * 60 * 24))
                    mu_rel0.append(abs(mu[i,3]/mu[i,1]))
            else:
                if snr[i,1]>snr_threshold[0] and snr[i,3]>snr_threshold[1]:
                    dt_rel90.append(abs(dt[i,3]-dt[i,1])/ (60 * 60 * 24))
                    mu_rel90.append(abs(mu[i,3]/mu[i,1]))
            if image_type[i,2]==image_type[i,3]:
                # snr check
                # below will also take care of the nan values
                if snr[i,2]>snr_threshold[0] and snr[i,3]>snr_threshold[1]:
                    dt_rel0.append(abs(dt[i,3]-dt[i,2])/ (60 * 60 * 24))
                    mu_rel0.append(abs(mu[i,3]/mu[i,2]))
            else:
                if snr[i,2]>snr_threshold[0] and snr[i,3]>snr_threshold[1]:
                    dt_rel90.append(abs(dt[i,3]-dt[i,2])/ (60 * 60 * 24))
                    mu_rel90.append(abs(mu[i,3]/mu[i,2]))

        return {
            "dt_rel0": np.array(dt_rel0), "mu_rel0": np.array(mu_rel0),
            "dt_rel90": np.array(dt_rel90), "mu_rel90": np.array(mu_rel90),
        }

    if classification_type == 'arrival_time':
        print('classification_type = arrival_time')
        print('make sure that the images are sorted wrt to arrival time')
        print('direct output from "ler" should be sorted')
        dt_12, dt_13, dt_14, dt_23, dt_24, dt_34 = [], [], [], [], [], []
        mu_12, mu_13, mu_14, mu_23, mu_24, mu_34 = [], [], [], [], [], []

        for i in range(len(image_type)):
            if snr[i,0]>snr_threshold[0] and snr[i,1]>snr_threshold[1]:
                dt_12.append(abs(dt[i,1]-dt[i,0])/ (60 * 60 * 24))
                mu_12.append(abs(mu[i,1]/mu[i,0]))
            if snr[i,0]>snr_threshold[0] and snr[i,2]>snr_threshold[1]:
                dt_13.append(abs(dt[i,2]-dt[i,0])/ (60 * 60 * 24))
                mu_13.append(abs(mu[i,2]/mu[i,0]))
            if snr[i,0]>snr_threshold[0] and snr[i,3]>snr_threshold[1]:
                dt_14.append(abs(dt[i,3]-dt[i,0])/ (60 * 60 * 24))
                mu_14.append(abs(mu[i,3]/mu[i,0]))
            if snr[i,1]>snr_threshold[0] and snr[i,2]>snr_threshold[1]:
                dt_23.append(abs(dt[i,2]-dt[i,1])/ (60 * 60 * 24))
                mu_23.append(abs(mu[i,2]/mu[i,1]))
            if snr[i,1]>snr_threshold[0] and snr[i,3]>snr_threshold[1]:
                dt_24.append(abs(dt[i,3]-dt[i,1])/ (60 * 60 * 24))
                mu_24.append(abs(mu[i,3]/mu[i,1]))
            if snr[i,2]>snr_threshold[0] and snr[i,3]>snr_threshold[1]:
                dt_34.append(abs(dt[i,3]-dt[i,2])/ (60 * 60 * 24))
                mu_34.append(abs(mu[i,3]/mu[i,2]))

        return {
            "dt_12": np.array(dt_12), "mu_12": np.array(mu_12),
            "dt_13": np.array(dt_13), "mu_13": np.array(mu_13),
            "dt_14": np.array(dt_14), "mu_14": np.array(mu_14),
            "dt_23": np.array(dt_23), "mu_23": np.array(mu_23),
            "dt_24": np.array(dt_24), "mu_24": np.array(mu_24),
            "dt_34": np.array(dt_34), "mu_34": np.array(mu_34),
        }

def mu_vs_dt_plot(
    x_array,
    y_array,
    xscale = 'log10',
    yscale = 'log10',
    alpha=0.6,
    extent=None,
    contour_levels=[10, 40, 68, 95],
    colors=['blue', 'blue', 'blue', 'blue', 'blue'],
):
    """
        Function to generate 2D KDE and plot the relative magnification vs time delay difference for lensed samples.

        Parameters
        ----------
        x_array : `float.array`
            x array.
        y_array : `float.array`
            y array.
        xscale : `str`
            x-axis scale.
            default xscale = 'log10'. other options: 'log', None.
        yscale : `str`
            y-axis scale.
            default yscale = 'log10'. other options: 'log', None.
        alpha : `float`
            transparency of the contour plot.
            default alpha = 0.6.
        extent : `list`
            extent of the plot.
            default extent = None. It will consider the full range of x_array and y_array.
        contour_levels : `list`
            levels for contour plot.
            default contour_levels = [10, 40, 68, 95].
        colors : `str`
            colors for contour plot.
            default colors = ['blue', 'blue', 'blue', 'blue', 'blue'].

        Examples
        ----------
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
    """

    x_min = min(x_array)
    x_max = max(x_array)
    y_min = min(y_array)
    y_max = max(y_array)

    # applying cutt-off
    if extent:
        x_min, x_max, y_min, y_max = extent
        x_array = x_array[(x_array >= x_min) & (x_array <= x_max)]
        y_array = y_array[(y_array >= y_min) & (y_array <= y_max)]
        
    # convert to log scale
    if xscale == 'log10':
        x_array = np.log10(x_array)
        x_min = np.log10(x_min)
        x_max = np.log10(x_max)
    if yscale == 'log10':
        y_array = np.log10(y_array)
        y_min = np.log10(y_min)
        y_max = np.log10(y_max)
    if xscale == 'log':
        x_array = np.log(x_array)
        x_min = np.log(x_min)
        x_max = np.log(x_max)
    if yscale == 'log':
        y_array = np.log(y_array)
        y_min = np.log(y_min)
        y_max = np.log(y_max)

    # Perform a kernel density estimation (KDE)
    xy = np.vstack([x_array, y_array])
    kde = gaussian_kde(xy)(xy)

    # Define the levels for contour as percentiles of the density
    levels = np.percentile(kde, [10, 40, 68, 95])

    # Create a grid for contour plot
    xgrid = np.linspace(x_min, x_max, 1000)
    ygrid = np.linspace(y_min, y_max, 1000)
    X1, Y1 = np.meshgrid(xgrid, ygrid)
    Z1 = gaussian_kde(xy)(np.vstack([X1.ravel(), Y1.ravel()])).reshape(X1.shape)

    if isinstance(colors, str):
        colors = [colors] * len(contour_levels)

    plt.contour(X1, Y1, Z1, levels=levels, colors=colors, alpha=alpha)
    