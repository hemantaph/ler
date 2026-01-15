# -*- coding: utf-8 -*-
"""
Lens functions for gravitational lensing calculations.

This module provides utility functions for computing lens galaxy velocity 
dispersion functions, ellipticity conversions, and strong lensing cross-sections.
These functions support the lens galaxy population modeling in the ``ler`` package.

The velocity dispersion functions follow the models from Oguri et al. (2018b) 
and Bernardi et al. (2010). Cross-section calculations use lenstronomy for 
EPL+Shear lens models.
"""

import numpy as np
from numba import njit
from lenstronomy.LensModel.Solver.epl_shear_solver import caustics_epl_shear
from shapely.geometry import Polygon


@njit
def phi_cut_SIE(q):
    """
    Calculate cross-section scaling factor for SIE lens galaxy from SIS.

    Computes the ratio of the SIE (Singular Isothermal Ellipsoid) cross-section
    to the SIS (Singular Isothermal Sphere) cross-section for a given axis ratio.

    Parameters
    ----------
    q : ``numpy.ndarray``
        Axis ratio of the lens galaxy (0 < q <= 1).

    Returns
    -------
    result : ``numpy.ndarray``
        Scaling factor (normalized to pi). \n
        For q -> 1 (spherical): returns 1.0 \n
        For q -> 0 (highly elliptical): returns ~0.
    """
    n = len(q)
    result = np.empty(n)
    for i in range(n):
        val = q[i]
        if 0.01 < val < 0.99:
            result[i] = (2 * np.pi * val * np.log(val)) / (val ** 2 - 1)
        elif val < 0.01:
            result[i] = -2 * (np.pi * np.log(val)) * val
        else:
            result[i] = np.pi
    return result/np.pi


@njit
def phi_q2_ellipticity(phi, q):
    """
    Convert position angle and axis ratio to ellipticity components.

    Parameters
    ----------
    phi : ``numpy.ndarray``
        Position angle of the major axis (radians).
    q : ``numpy.ndarray``
        Axis ratio (0 < q <= 1).

    Returns
    -------
    e1 : ``numpy.ndarray``
        First ellipticity component.
    e2 : ``numpy.ndarray``
        Second ellipticity component.
    """
    e_1 = (1.0 - q) / (1.0 + q) * np.cos(2 * phi)
    e_2 = (1.0 - q) / (1.0 + q) * np.sin(2 * phi)
    return e_1, e_2


def cross_section(theta_E, e1, e2, gamma, gamma1, gamma2):
    """
    Compute the strong lensing cross-section for an EPL+Shear lens.

    Uses lenstronomy to compute the caustic structure and returns the 
    area enclosed by the double-image (outer) caustic.

    Parameters
    ----------
    theta_E : ``float``
        Einstein radius (arcseconds).
    e1 : ``float``
        First ellipticity component.
    e2 : ``float``
        Second ellipticity component.
    gamma : ``float``
        Power-law density profile slope.
    gamma1 : ``float``
        First external shear component.
    gamma2 : ``float``
        Second external shear component.

    Returns
    -------
    area : ``float``
        Cross-section area (arcseconds^2). \n
        Returns 0.0 if caustic computation fails.
    """
    kwargs_lens = [
        {
            "theta_E": theta_E,
            "e1": e1,
            "e2": e2,
            "gamma": gamma,
            "center_x": 0.0,
            "center_y": 0.0,
        },
        {
            "gamma1": gamma1,
            "gamma2": gamma2,
            "ra_0": 0,
            "dec_0": 0,
        },
    ]
    
    caustic_double_points = caustics_epl_shear(
        kwargs_lens, return_which="double", maginf=-100
    )
    caustic = np.logical_not(np.isnan(caustic_double_points).any())

    if caustic:
        area = Polygon(caustic_double_points.T).area
    else:
        area = 0.0

    return area
