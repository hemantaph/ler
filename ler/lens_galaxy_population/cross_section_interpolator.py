# -*- coding: utf-8 -*-
"""
Module for gravitational lensing cross section interpolation.

This module provides highly optimized Numba-compiled functions for interpolating
gravitational lensing cross sections in a 5-dimensional parameter space using 
cubic B-spline interpolation on prefiltered coefficient grids.

The interpolation is performed based on: \n
- Ellipticity components (e1, e2) derived from axis ratio q and position angle phi \n
- Density profile slope (gamma) \n
- External shear components (gamma1, gamma2) \n

Key Features: \n
- JIT-compiled with Numba for high performance \n
- Parallel processing support for batch evaluations \n
- Cubic B-spline interpolation matching scipy's map_coordinates behavior \n
- Automatic scaling by Einstein radius and affine calibration \n

Usage:
    Basic workflow example:

    >>> from ler.lens_galaxy_population.cross_section_interpolator import make_cross_section_reinit
    >>> cs_func = make_cross_section_reinit(e1_grid, e2_grid, gamma_grid, ...)
    >>> cross_sections = cs_func(zs, zl, sigma, q, phi, gamma, gamma1, gamma2)

Copyright (C) 2024 Hemanta Kumar Phurailatpam. Distributed under MIT License.
"""

import numpy as np
from numba import njit, prange
from scipy import ndimage

from .lens_functions import phi_q2_ellipticity

C_LIGHT = 299792.458  # km/s


@njit
def _bspline3(t):
    """
    Evaluate the centered cubic B-spline basis function B3(t).

    Parameters
    ----------
    t : ``float``
        Distance from the spline center. \n
        The function has compact support on [-2, 2].

    Returns
    -------
    value : ``float``
        The B-spline basis value at position t. \n
        - Maximum value of 2/3 at t=0 \n
        - Smoothly decays to 0 at |t|=2 \n
        - Returns 0 for |t| >= 2 \n
    """
    at = abs(t)
    if at < 1.0:
        return (4.0 - 6.0 * at * at + 3.0 * at * at * at) / 6.0
    elif at < 2.0:
        u = 2.0 - at
        return (u * u * u) / 6.0
    else:
        return 0.0


@njit
def _clamp_int(i, n):
    """
    Clamp an index to valid array bounds using 'nearest' boundary mode.

    Parameters
    ----------
    i : ``int``
        The index to clamp. Can be negative or exceed array bounds.
    n : ``int``
        The size of the array dimension (valid indices are 0 to n-1).

    Returns
    -------
    clamped_index : ``int``
        The clamped index in the range [0, n-1]. \n
        - If i < 0, returns 0 \n
        - If i > n-1, returns n-1 \n
        - Otherwise returns i unchanged \n
    """
    if i < 0:
        return 0
    if i > n - 1:
        return n - 1
    return i


@njit
def _physical_to_index(x, x0, x1, n):
    """
    Convert a physical coordinate to a fractional grid index.

    Parameters
    ----------
    x : ``float``
        The physical coordinate value to convert.
    x0 : ``float``
        The minimum physical coordinate (maps to index 0).
    x1 : ``float``
        The maximum physical coordinate (maps to index n-1).
    n : ``int``
        The number of grid points in this dimension.

    Returns
    -------
    fractional_index : ``float``
        The fractional index in [0, n-1]. \n
        Non-integer values indicate positions between grid points.
    """
    return (x - x0) * (n - 1) / (x1 - x0)


# @njit(parallel=True)
@njit
def _map_coordinates_5d_cubic_nearest(coeff, ie1, ie2, ig, ig1, ig2):
    """
    Perform 5D cubic B-spline interpolation on prefiltered coefficients.

    This function evaluates a cubic (order-3) B-spline interpolation in 5 dimensions.
    It matches the behavior of ``scipy.ndimage.map_coordinates`` with order=3,
    mode='nearest', and prefilter=False.

    Parameters
    ----------
    coeff : ``numpy.ndarray``
        Prefiltered B-spline coefficients with shape (n1, n2, n3, n4, n5). \n
        Should be the output of ``scipy.ndimage.spline_filter(data, order=3)``.
    ie1 : ``numpy.ndarray``
        Fractional indices for the first dimension (e1), shape (N,).
    ie2 : ``numpy.ndarray``
        Fractional indices for the second dimension (e2), shape (N,).
    ig : ``numpy.ndarray``
        Fractional indices for the third dimension (gamma), shape (N,).
    ig1 : ``numpy.ndarray``
        Fractional indices for the fourth dimension (gamma1), shape (N,).
    ig2 : ``numpy.ndarray``
        Fractional indices for the fifth dimension (gamma2), shape (N,).

    Returns
    -------
    out : ``numpy.ndarray``
        Interpolated values at the specified fractional indices, shape (N,).
    """
    n1, n2, n3, n4, n5 = coeff.shape
    N = ie1.size
    out = np.empty(N, dtype=np.float64)

    for k in prange(N):
        x1 = ie1[k]
        x2 = ie2[k]
        x3 = ig[k]
        x4 = ig1[k]
        x5 = ig2[k]

        b1 = int(np.floor(x1)) - 1
        b2 = int(np.floor(x2)) - 1
        b3 = int(np.floor(x3)) - 1
        b4 = int(np.floor(x4)) - 1
        b5 = int(np.floor(x5)) - 1

        acc = 0.0

        for a1 in range(4):
            j1 = _clamp_int(b1 + a1, n1)
            w1 = _bspline3(x1 - (b1 + a1))
            if w1 == 0.0:
                continue

            for a2 in range(4):
                j2 = _clamp_int(b2 + a2, n2)
                w2 = _bspline3(x2 - (b2 + a2))
                if w2 == 0.0:
                    continue

                for a3 in range(4):
                    j3 = _clamp_int(b3 + a3, n3)
                    w3 = _bspline3(x3 - (b3 + a3))
                    if w3 == 0.0:
                        continue

                    for a4 in range(4):
                        j4 = _clamp_int(b4 + a4, n4)
                        w4 = _bspline3(x4 - (b4 + a4))
                        if w4 == 0.0:
                            continue

                        for a5 in range(4):
                            j5 = _clamp_int(b5 + a5, n5)
                            w5 = _bspline3(x5 - (b5 + a5))
                            if w5 == 0.0:
                                continue

                            acc += (w1 * w2 * w3 * w4 * w5) * coeff[j1, j2, j3, j4, j5]

        out[k] = acc

    return out


@njit
def _cross_section(
    zs,
    zl,
    sigma,
    q,
    phi,
    gamma,
    gamma1,
    gamma2,
    e1_grid,
    e2_grid,
    gamma_grid,
    gamma1_grid,
    gamma2_grid,
    cs_spline_coeff_grid,
    Da_instance,
    csunit_to_cs_slope,
    csunit_to_cs_intercept,
):
    """
    Compute lensing cross sections via 5D B-spline interpolation.

    Parameters
    ----------
    zs : ``numpy.ndarray``
        Source redshifts, shape (N,).
    zl : ``numpy.ndarray``
        Lens redshifts, shape (N,).
    sigma : ``numpy.ndarray``
        Velocity dispersions (units: km/s), shape (N,).
    q : ``numpy.ndarray``
        Axis ratios, shape (N,).
    phi : ``numpy.ndarray``
        Position angles (units: radians), shape (N,).
    gamma : ``numpy.ndarray``
        Density profile slopes, shape (N,).
    gamma1 : ``numpy.ndarray``
        External shear component 1, shape (N,).
    gamma2 : ``numpy.ndarray``
        External shear component 2, shape (N,).
    e1_grid : ``numpy.ndarray``
        Grid values for ellipticity component e1, shape (n_e1,).
    e2_grid : ``numpy.ndarray``
        Grid values for ellipticity component e2, shape (n_e2,).
    gamma_grid : ``numpy.ndarray``
        Grid values for density slope gamma, shape (n_g,).
    gamma1_grid : ``numpy.ndarray``
        Grid values for shear component gamma1, shape (n_g1,).
    gamma2_grid : ``numpy.ndarray``
        Grid values for shear component gamma2, shape (n_g2,).
    cs_spline_coeff_grid : ``numpy.ndarray``
        Prefiltered B-spline coefficients for cross section, \n
        shape (n_e1, n_e2, n_g, n_g1, n_g2).
    Da_instance : ``callable``
        Angular diameter distance function. \n
        Signature: ``Da_instance(z) -> distance``
    csunit_to_cs_slope : ``float``
        Slope for affine calibration from unit cross section.
    csunit_to_cs_intercept : ``float``
        Intercept for affine calibration from unit cross section.

    Returns
    -------
    cs : ``numpy.ndarray``
        Computed cross sections (units: radians^2), shape (N,). \n
        Negative values are clipped to zero.
    """
    # Angular diameter distances
    Ds = Da_instance(zs)
    Dl = Da_instance(zl)
    Dls = (Ds * (1 + zs) - Dl * (1 + zl)) / (1 + zs)

    # Einstein radius
    theta_E = 4.0 * np.pi * (sigma / C_LIGHT) ** 2 * (Dls / Ds)

    # Convert (q, phi) -> (e1, e2)
    e1, e2 = phi_q2_ellipticity(phi, q)

    N = zs.size

    # Compute fractional indices
    ie1 = np.empty(N, dtype=np.float64)
    ie2 = np.empty(N, dtype=np.float64)
    ig = np.empty(N, dtype=np.float64)
    ig1 = np.empty(N, dtype=np.float64)
    ig2 = np.empty(N, dtype=np.float64)

    e1_min, e1_max = e1_grid[0], e1_grid[-1]
    e2_min, e2_max = e2_grid[0], e2_grid[-1]
    g_min, g_max = gamma_grid[0], gamma_grid[-1]
    g1_min, g1_max = gamma1_grid[0], gamma1_grid[-1]
    g2_min, g2_max = gamma2_grid[0], gamma2_grid[-1]

    n_e1 = e1_grid.size
    n_e2 = e2_grid.size
    n_g = gamma_grid.size
    n_g1 = gamma1_grid.size
    n_g2 = gamma2_grid.size

    for i in range(N):
        ie1[i] = _physical_to_index(e1[i], e1_min, e1_max, n_e1)
        ie2[i] = _physical_to_index(e2[i], e2_min, e2_max, n_e2)
        ig[i] = _physical_to_index(gamma[i], g_min, g_max, n_g)
        ig1[i] = _physical_to_index(gamma1[i], g1_min, g1_max, n_g1)
        ig2[i] = _physical_to_index(gamma2[i], g2_min, g2_max, n_g2)

    # Unit cross-section via cubic spline interpolation
    cs_unit = _map_coordinates_5d_cubic_nearest(cs_spline_coeff_grid, ie1, ie2, ig, ig1, ig2)

    # Scale by SIS area and apply affine calibration
    area_sis = np.pi * theta_E**2
    cs = cs_unit * (csunit_to_cs_intercept + csunit_to_cs_slope * area_sis)
    cs = np.maximum(cs, 0.0)

    return cs


def make_cross_section_reinit(
    e1_grid,
    e2_grid,
    gamma_grid,
    gamma1_grid,
    gamma2_grid,
    cs_spline_coeff_grid,
    Da_instance,
    csunit_to_cs_slope=0.31830988618379075,
    csunit_to_cs_intercept=-3.2311742677852644e-27,
):
    """
    Factory function to create a JIT-compiled cross section calculator.

    This function precomputes B-spline coefficients and creates a closure 
    that captures the grid parameters, returning a fast Numba-compiled 
    function for computing cross sections.

    Parameters
    ----------
    e1_grid : ``numpy.ndarray``
        Grid values for ellipticity component e1, shape (n_e1,).
    e2_grid : ``numpy.ndarray``
        Grid values for ellipticity component e2, shape (n_e2,).
    gamma_grid : ``numpy.ndarray``
        Grid values for density slope gamma, shape (n_g,).
    gamma1_grid : ``numpy.ndarray``
        Grid values for shear component gamma1, shape (n_g1,).
    gamma2_grid : ``numpy.ndarray``
        Grid values for shear component gamma2, shape (n_g2,).
    cs_spline_coeff_grid : ``numpy.ndarray``
        Raw cross section grid data (before spline filtering), \n
        shape (n_e1, n_e2, n_g, n_g1, n_g2).
    Da_instance : ``callable``
        Angular diameter distance function. \n
        Signature: ``Da_instance(z) -> distance``
    csunit_to_cs_slope : ``float``
        Slope for affine calibration from unit cross section. \n
        default: 0.31830988618379075
    csunit_to_cs_intercept : ``float``
        Intercept for affine calibration from unit cross section. \n
        default: -3.2311742677852644e-27

    Returns
    -------
    cross_section_reinit : ``callable``
        JIT-compiled function with signature: \n
        ``cross_section_reinit(zs, zl, sigma, q, phi, gamma, gamma1, gamma2)`` \n
        Returns cross sections as ``numpy.ndarray`` of shape (N,).

    Examples
    --------
    >>> from ler.lens_galaxy_population.cross_section_interpolator import make_cross_section_reinit
    >>> cs_func = make_cross_section_reinit(
    ...     e1_grid, e2_grid, gamma_grid, gamma1_grid, gamma2_grid,
    ...     cs_spline_coeff_grid, Da_instance
    ... )
    >>> cross_sections = cs_func(zs, zl, sigma, q, phi, gamma, gamma1, gamma2)
    """
    e1g = np.asarray(e1_grid, dtype=np.float64)
    e2g = np.asarray(e2_grid, dtype=np.float64)
    gg = np.asarray(gamma_grid, dtype=np.float64)
    g1g = np.asarray(gamma1_grid, dtype=np.float64)
    g2g = np.asarray(gamma2_grid, dtype=np.float64)
    cs = np.asarray(cs_spline_coeff_grid, dtype=np.float64)

    # Precompute spline coefficients
    cs_coeff = ndimage.spline_filter(cs, order=3)

    slope_given = float(csunit_to_cs_slope)
    intercept_given = float(csunit_to_cs_intercept)

    Da_function = njit(lambda z: Da_instance(z))

    @njit
    def cross_section_reinit(zs, zl, sigma, q, phi, gamma, gamma1, gamma2):
        return _cross_section(
            zs=zs,
            zl=zl,
            sigma=sigma,
            q=q,
            phi=phi,
            gamma=gamma,
            gamma1=gamma1,
            gamma2=gamma2,
            e1_grid=e1g,
            e2_grid=e2g,
            gamma_grid=gg,
            gamma1_grid=g1g,
            gamma2_grid=g2g,
            cs_spline_coeff_grid=cs_coeff,
            Da_instance=Da_function,
            csunit_to_cs_slope=slope_given,
            csunit_to_cs_intercept=intercept_given,
        )

    return cross_section_reinit
