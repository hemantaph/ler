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

def relative_mu_dt_unlensed(param, size=100):
    """
    Function to generate relative magnification vs time delay difference for unlensed samples.

    Parameters
    ----------
    param : `dict`
        dictionary of unlensed GW source parameters.
        unlensed_param.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time']

    Returns
    ----------
    dmu : `float.array`
        relative magnification.
    dt : `float.array`
        relative time delay.

    """

    t = param["geocent_time"]
    mu = param["luminosity_distance"]

    len_ = len(t)
    t_ = []
    mu_ = []
    while len(t_) < size:
        idx1 = np.random.choice(np.arange(0,len_), size, replace=False)
        idx2 = np.random.choice(np.arange(0,len_), size, replace=False)
        t_.append(t[idx2] - t[idx1])
        mu_.append(mu[idx2] / mu[idx1])

    dt = np.abs(np.array(t_)) / (60 * 60 * 24)  # in days
    dmu = np.sqrt(np.abs(np.array(mu_)))

    return (dmu, dt)

def relative_mu_dt_lensed(lensed_param, snr_threshold=[8.0, 8.0]):
    """
    Function to classify the lensed images wrt to the morse phase difference.

    Parameters
    ----------
    lensed_param : `dict`
        dictionary of lensed GW source parameters, lens galaxy parameters and image paramters.
        lensed_param.keys() = ['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2', 'Dl',
        'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source',
        'luminosity_distance', 'theta_jn', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'n_images',
        'x0_image_positions', 'x1_image_positions', 'magnifications', 'time_delays', 'traces',
        'determinants', 'image_type', 'weights', 'optimal_snr_net', 'L1', 'H1', 'V1']
    snr_threshold : `float`
        threshold for detection signal to noise ratio.
        e.g. snr_threshold = [8.,8.] or [8.,6.] for subthreshold

    Returns
    ----------
    mu_rel0 : `float.array`
        relative magnification for 0 degree phase difference.
    dt_rel0 : `float.array`
        relative time delay for 0 degree phase difference.
    mu_rel90 : `float.array`
        relative magnification for 90 degree phase difference.
    dt_rel90 : `float.array`
        relative time delay for 90 degree phase difference.
    """

    # get magnifications, time_delays and snr
    mu = np.nan_to_num(lensed_param["magnifications"])
    dt = np.nan_to_num(lensed_param["time_delays"])
    snr = np.nan_to_num(lensed_param["optimal_snr_net"])

    # for 0 degree phase difference
    # get the index of the image which cross the threshold
    # get snr_threshold sorted first in descending order
    snr_threshold = -np.sort(-np.array(snr_threshold))
    # for type I
    snr1 = -np.sort(-snr[:, [0, 1]], axis=1)
    # for type II
    snr2 = -np.sort(-snr[:, [2, 3]], axis=1)

    # checking for zero values
    # check for threshold condition
    idx1, idx2 = [], []
    for i in range(len(snr)):
        if (
            any(x != 0.0 for x in snr1[i])
            and snr1[i][0] > snr_threshold[0]
            and snr1[i][1] > snr_threshold[1]
        ):
            idx1.append(i)
        if (
            any(x != 0.0 for x in snr2[i])
            and snr2[i][0] > snr_threshold[0]
            and snr2[i][1] > snr_threshold[1]
        ):
            idx2.append(i)

    # combine magnifications and time_delays
    mu_ = np.concatenate((mu[idx1][:, [0, 1]], mu[idx2][:, [2, 3]]), axis=0)
    dt_ = np.concatenate((dt[idx1][:, [0, 1]], dt[idx2][:, [2, 3]]), axis=0) / (
        60 * 60 * 24
    )  # to days

    # relative magnification
    mu_rel0 = np.abs(mu_[:, 1] / mu_[:, 0])
    # relative time delay
    dt_rel0 = np.abs(dt_[:, 1] - dt_[:, 0])

    # for 90 degree phase difference
    # for type I
    snr1 = -np.sort(-snr[:, [0, 2]], axis=1)
    # for type II
    snr2 = -np.sort(-snr[:, [1, 3]], axis=1)

    # checking for zero values
    # check for threshold condition
    idx1, idx2 = [], []
    for i in range(len(snr)):
        if (
            any(x != 0.0 for x in snr1[i])
            and snr1[i][0] > snr_threshold[0]
            and snr1[i][1] > snr_threshold[1]
        ):
            idx1.append(i)
        if (
            any(x != 0.0 for x in snr2[i])
            and snr2[i][0] > snr_threshold[0]
            and snr2[i][1] > snr_threshold[1]
        ):
            idx2.append(i)

    # combine magnifications and time_delays
    mu_ = np.concatenate((mu[idx1][:, [0, 2]], mu[idx2][:, [1, 3]]), axis=0)
    dt_ = np.concatenate((dt[idx1][:, [0, 2]], dt[idx2][:, [1, 3]]), axis=0) / (
        60 * 60 * 24
    )  # in days

    # relative magnification
    mu_rel90 = np.abs(mu_[:, 1] / mu_[:, 0])
    # relative time delay
    dt_rel90 = np.abs(dt_[:, 1] - dt_[:, 0])

    return (mu_rel0, dt_rel0, mu_rel90, dt_rel90)

def mu_vs_dt_plot(
    x_array,
    y_array,
    savefig=False,
    ax=None,
    colors="blue",
    linestyles="-",
    origin="upper",
    alpha=0.6,
    extent=[1e-2, 5e2, 1e-2, 1e2],
    contour_levels=[0.10, 0.40, 0.68, 0.95],
):
    """
    Function to generate 2D KDE and plot the relative magnification vs time delay difference for lensed samples.

    Parameters
    ----------
    x_array : `float.array`
        x array.
    y_array : `float.array`
        y array.
    xlabel : `str`
        x label.
    ylabel : `str`
        y label.
    title : `str`
        title.
    savefig : `bool`
        if True, it will save the figure.
        default savefig = False.
    ax : `matplotlib.axes`
        matplotlib axes.
        default ax = None.
    colors : `str`
        color of the plot.
        default colors = 'blue'.
    linestyles : `str`
        linestyle of the plot.
        default linestyles = '-'.
    origin : `str`
        origin of the plot.
        default origin = 'upper'.
    alpha : `float`
        alpha of the plot.
        default alpha = 0.6.
    extent : `list`
        extent of the plot.
        default extent = [1e-2,5e2,1e-2,1e2].
    contour_levels : `list`
        contour levels of the plot.
        default contour_levels = [0.10,0.40,0.68,0.95] which corresponds to 1,2,3,4 sigma.

    Returns
    ----------
    None

    """
    # applying cutt-off
    idx = (
        (x_array > extent[0])
        & (x_array < extent[1])
        & (y_array > extent[2])
        & (y_array < extent[3])
    )
    x_array = x_array[idx]
    y_array = y_array[idx]

    xu = np.log10(x_array)
    yu = np.log10(y_array)

    xmin = np.log10(1e-2)
    xmax = np.log10(5e2)
    ymin = np.log10(1e-2)
    ymax = np.log10(1e2)

    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])
    values = np.vstack([xu, yu])
    kernel = gaussian_kde(values)
    ff = np.reshape(kernel(positions).T, xx.shape)

    zsort = -np.sort(-ff.flatten())

    cumz = np.cumsum(zsort) / np.sum(zsort)
    spl = interp1d(cumz, zsort, kind="cubic", fill_value="extrapolate")

    levels = []
    for i in contour_levels:
        levels.append(spl(i))
    levels = np.array(levels)[::-1]

    ax.contour(
        np.rot90(ff),
        levels,
        colors=colors,
        linestyles=linestyles,
        origin=origin,
        alpha=alpha,
        extent=np.log10(extent),
    )

    # labels
    ax.xlabel(r"$log_{10}\Delta t$ (days)")
    ax.ylabel(r"$\Delta log_{10}\mu$")
    ax.title(r"relative magnification vs relative time delay")

    # save figure
    if savefig:
        ax.savefig("mu_vs_dt.png", dpi=300, bbox_inches="tight")

    return None
    