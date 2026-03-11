# -*- coding: utf-8 -*-
"""
Module for interpolating gravitational-lensing cross sections.

This module provides fast Numba-compiled interpolation in a 5D parameter
space for lensing cross sections. The interpolated dimensions are ellipticity
(``e1``, ``e2``), density slope (``gamma``), and external shear
(``gamma1``, ``gamma2``). Inputs from lens parameters are mapped to this grid,
then rescaled by Einstein-radius geometry and affine calibration.

Key Features: \n
- 5D local interpolation using precomputed tensor basis weights \n
- Zero memory allocation inside execution loops for maximum throughput \n
- Numba JIT compilation with fused parallel loops \n
- Direct conversion from (q, phi) to (e1, e2) before interpolation \n

Copyright (C) 2024 Hemantakumar Phurailatpam. Distributed under MIT License.
"""

import numpy as np
from numba import njit, prange

from .lens_functions import phi_q2_ellipticity

C_LIGHT = 299792.458  # km/s


@njit(cache=True)
def _precompute_1d_coeffs(x_array):
    """
    Precompute grid invariants (spacings, ratios) for cubic Hermite interpolation.
    Runs once during factory initialization.
    """
    n = x_array.shape[0]
    m1_0 = np.zeros(n, dtype=np.float64)
    m1_1 = np.zeros(n, dtype=np.float64)
    m1_2 = np.zeros(n, dtype=np.float64)
    m2_1 = np.zeros(n, dtype=np.float64)
    m2_2 = np.zeros(n, dtype=np.float64)
    m2_3 = np.zeros(n, dtype=np.float64)
    denom = np.zeros(n, dtype=np.float64)

    for i in range(1, n - 2):
        x0, x1, x2, x3 = x_array[i-1], x_array[i], x_array[i+1], x_array[i+2]
        dx10 = x1 - x0
        dx21 = x2 - x1
        dx20 = x2 - x0
        dx32 = x3 - x2
        dx31 = x3 - x1

        denom[i] = dx21

        w10 = dx21 / dx20 if dx20 != 0.0 else 0.0
        w21_left = dx10 / dx20 if dx20 != 0.0 else 0.0
        w21_right = dx32 / dx31 if dx31 != 0.0 else 0.0
        w32 = dx21 / dx31 if dx31 != 0.0 else 0.0

        m1_0[i] = -w10 / dx10 if dx10 != 0.0 else 0.0
        m1_1[i] = (w10 / dx10 if dx10 != 0.0 else 0.0) - (w21_left / dx21 if dx21 != 0.0 else 0.0)
        m1_2[i] = w21_left / dx21 if dx21 != 0.0 else 0.0

        m2_1[i] = -w21_right / dx21 if dx21 != 0.0 else 0.0
        m2_2[i] = (w21_right / dx21 if dx21 != 0.0 else 0.0) - (w32 / dx32 if dx32 != 0.0 else 0.0)
        m2_3[i] = w32 / dx32 if dx32 != 0.0 else 0.0

    return m1_0, m1_1, m1_2, m2_1, m2_2, m2_3, denom


@njit(cache=True, fastmath=True, inline="always")
def _get_1d_weights(x_eval, x_grid, m1_0, m1_1, m1_2, m2_1, m2_2, m2_3, denom_arr):
    """
    Finds the base 4-point index and returns 4 scalar weights (C0, C1, C2, C3)
    representing the tensor basis for interpolation. 
    Uses full Cubic Hermite everywhere, including one-sided finite differences for edges.
    """
    n = x_grid.shape[0]

    i = np.searchsorted(x_grid, x_eval, side="right") - 1
    if i < 1:
        i = 1
    elif i > n - 3:
        i = n - 3

    # -----------------------------------------------------------------
    # Left edge / Extrapolation (Cubic Hermite via Forward Difference)
    # -----------------------------------------------------------------
    if x_eval <= x_grid[1]:
        x0 = x_grid[0]
        x1 = x_grid[1]
        x2 = x_grid[2]
        
        h0 = x1 - x0
        h1 = x2 - x1
        
        t = (x_eval - x0) / h0 if h0 != 0.0 else 0.0
        t2 = t * t
        t3 = t2 * t
        
        h00 = 2.0 * t3 - 3.0 * t2 + 1.0
        h10 = t3 - 2.0 * t2 + t
        h01 = -2.0 * t3 + 3.0 * t2
        h11 = t3 - t2
        
        # Tangent at x0 using 3-point forward difference
        d0 = -(2.0 * h0 + h1) / (h0 * (h0 + h1)) if h0 != 0.0 else 0.0
        d1 = (h0 + h1) / (h0 * h1) if h0 != 0.0 else 0.0
        d2 = -h0 / (h1 * (h0 + h1)) if h1 != 0.0 else 0.0
        
        # Tangent at x1 using 3-point central difference
        e0 = -h1 / (h0 * (h0 + h1)) if h0 != 0.0 else 0.0
        e1 = (h1 - h0) / (h0 * h1) if h0 != 0.0 else 0.0
        e2 = h0 / (h1 * (h0 + h1)) if h1 != 0.0 else 0.0
        
        d_h10 = h10 * h0
        d_h11 = h11 * h0
        
        # Map to base index 0: y0, y1, y2, y3 (C3 is unused)
        C0 = h00 + d_h10 * d0 + d_h11 * e0
        C1 = h01 + d_h10 * d1 + d_h11 * e1
        C2 = d_h10 * d2 + d_h11 * e2
        C3 = 0.0
        
        return 0, C0, C1, C2, C3
        
    # -----------------------------------------------------------------
    # Right edge / Extrapolation (Cubic Hermite via Backward Difference)
    # -----------------------------------------------------------------
    elif x_eval >= x_grid[n - 2]:
        x1 = x_grid[n - 3]
        x2 = x_grid[n - 2]
        x3 = x_grid[n - 1]
        
        h1 = x2 - x1
        h2 = x3 - x2
        
        t = (x_eval - x2) / h2 if h2 != 0.0 else 0.0
        t2 = t * t
        t3 = t2 * t
        
        h00 = 2.0 * t3 - 3.0 * t2 + 1.0
        h10 = t3 - 2.0 * t2 + t
        h01 = -2.0 * t3 + 3.0 * t2
        h11 = t3 - t2
        
        # Tangent at x2 using 3-point central difference
        e1 = -h2 / (h1 * (h1 + h2)) if h1 != 0.0 else 0.0
        e2 = (h2 - h1) / (h1 * h2) if h1 != 0.0 else 0.0
        e3 = h1 / (h2 * (h1 + h2)) if h2 != 0.0 else 0.0
        
        # Tangent at x3 using 3-point backward difference
        d1 = h2 / (h1 * (h1 + h2)) if h1 != 0.0 else 0.0
        d2 = -(h1 + h2) / (h1 * h2) if h1 != 0.0 else 0.0
        d3 = (2.0 * h2 + h1) / (h2 * (h1 + h2)) if h2 != 0.0 else 0.0
        
        d_h10 = h10 * h2
        d_h11 = h11 * h2
        
        # Map to base index n-4: y_{n-4}, y_{n-3}, y_{n-2}, y_{n-1}
        C0 = 0.0
        C1 = d_h10 * e1 + d_h11 * d1
        C2 = h00 + d_h10 * e2 + d_h11 * d2
        C3 = h01 + d_h10 * e3 + d_h11 * d3
        
        return n - 4, C0, C1, C2, C3
        
    # -----------------------------------------------------------------
    # Interior (Standard precomputed Cubic Hermite)
    # -----------------------------------------------------------------
    else:
        x1 = x_grid[i]
        denom = denom_arr[i]
        t = (x_eval - x1) / denom if denom != 0.0 else 0.0
        t2 = t * t
        t3 = t2 * t

        h00 = 2.0 * t3 - 3.0 * t2 + 1.0
        h10 = t3 - 2.0 * t2 + t
        h01 = -2.0 * t3 + 3.0 * t2
        h11 = t3 - t2

        d_h10 = h10 * denom
        d_h11 = h11 * denom

        C0 = d_h10 * m1_0[i]
        C1 = h00 + d_h10 * m1_1[i] + d_h11 * m2_1[i]
        C2 = h01 + d_h10 * m1_2[i] + d_h11 * m2_2[i]
        C3 = d_h11 * m2_3[i]
        
        return i - 1, C0, C1, C2, C3


@njit(cache=True, fastmath=True, parallel=True)
def _cross_section_eval(
    zs, zl, ds_arr, dl_arr, sigma, q, phi, gamma, gamma1, gamma2,
    e1_grid, e2_grid, gamma_grid, gamma1_grid, gamma2_grid, cs_unit_grid,
    e1_c, e2_c, g_c, g1_c, g2_c,
    csunit_to_cs_slope, csunit_to_cs_intercept
):
    N = zs.shape[0]
    cross_section = np.empty(N, dtype=np.float64)
    e1_arr, e2_arr = phi_q2_ellipticity(phi, q)

    for i in prange(N):
        # Read precomputed distances from the arrays
        ds = ds_arr[i]
        dl = dl_arr[i]
        
        # Physical computation
        dsl = (ds * (1.0 + zs[i]) - dl * (1.0 + zl[i])) / (1.0 + zs[i])
        einstein_radius = 4.0 * np.pi * (sigma[i] / C_LIGHT)**2 * (dsl / ds)
        sis_area = np.pi * einstein_radius**2

        # 1D Index & Weight lookups
        idx_e1, we1_0, we1_1, we1_2, we1_3 = _get_1d_weights(e1_arr[i], e1_grid, e1_c[0], e1_c[1], e1_c[2], e1_c[3], e1_c[4], e1_c[5], e1_c[6])
        idx_e2, we2_0, we2_1, we2_2, we2_3 = _get_1d_weights(e2_arr[i], e2_grid, e2_c[0], e2_c[1], e2_c[2], e2_c[3], e2_c[4], e2_c[5], e2_c[6])
        idx_g,  wg_0,  wg_1,  wg_2,  wg_3  = _get_1d_weights(gamma[i], gamma_grid, g_c[0], g_c[1], g_c[2], g_c[3], g_c[4], g_c[5], g_c[6])
        idx_g1, wg1_0, wg1_1, wg1_2, wg1_3 = _get_1d_weights(gamma1[i], gamma1_grid, g1_c[0], g1_c[1], g1_c[2], g1_c[3], g1_c[4], g1_c[5], g1_c[6])
        idx_g2, wg2_0, wg2_1, wg2_2, wg2_3 = _get_1d_weights(gamma2[i], gamma2_grid, g2_c[0], g2_c[1], g2_c[2], g2_c[3], g2_c[4], g2_c[5], g2_c[6])

        w_e1 = (we1_0, we1_1, we1_2, we1_3)
        w_e2 = (we2_0, we2_1, we2_2, we2_3)
        w_g  = (wg_0,  wg_1,  wg_2,  wg_3)
        w_g1 = (wg1_0, wg1_1, wg1_2, wg1_3)
        w_g2 = (wg2_0, wg2_1, wg2_2, wg2_3)

        unit_cs_val = 0.0
        for i1 in range(4):
            we1 = w_e1[i1]
            if we1 == 0.0: continue
            for i2 in range(4):
                we2 = w_e2[i2]
                if we2 == 0.0: continue
                for i3 in range(4):
                    wg = w_g[i3]
                    if wg == 0.0: continue
                    for i4 in range(4):
                        wg1 = w_g1[i4]
                        if wg1 == 0.0: continue
                        for i5 in range(4):
                            wg2 = w_g2[i5]
                            if wg2 == 0.0: continue

                            unit_cs_val += we1 * we2 * wg * wg1 * wg2 * cs_unit_grid[
                                idx_e1 + i1, idx_e2 + i2, idx_g  + i3, idx_g1 + i4, idx_g2 + i5
                            ]

        cs = unit_cs_val * (csunit_to_cs_intercept + csunit_to_cs_slope * sis_area)
        cross_section[i] = cs if cs > 0.0 else 0.0

    return cross_section


def make_cross_section_reinit(
    e1_grid, e2_grid, gamma_grid, gamma1_grid, gamma2_grid,
    cs_unit_grid, Da_instance,
    csunit_to_cs_slope=0.31830988618379075,
    csunit_to_cs_intercept=-3.2311742677852644e-27,
):
    # ... (Keep the grid precomputations the exact same as before) ...
    e1_grid_val = np.asarray(e1_grid, dtype=np.float64)
    e2_grid_val = np.asarray(e2_grid, dtype=np.float64)
    gamma_grid_val = np.asarray(gamma_grid, dtype=np.float64)
    gamma1_grid_val = np.asarray(gamma1_grid, dtype=np.float64)
    gamma2_grid_val = np.asarray(gamma2_grid, dtype=np.float64)
    cs_unit_grid_val = np.asarray(cs_unit_grid, dtype=np.float64)

    e1_c = _precompute_1d_coeffs(e1_grid_val)
    e2_c = _precompute_1d_coeffs(e2_grid_val)
    g_c  = _precompute_1d_coeffs(gamma_grid_val)
    g1_c = _precompute_1d_coeffs(gamma1_grid_val)
    g2_c = _precompute_1d_coeffs(gamma2_grid_val)

    slope = float(csunit_to_cs_slope)
    intercept = float(csunit_to_cs_intercept)

    @njit
    def angular_diameter_distance_numba(z_array):
        # Now explicitly accepts and processes the full NumPy array
        return Da_instance(z_array)

    @njit(cache=True, fastmath=True)
    def cross_section_reinit(zs, zl, sigma, q, phi, gamma, gamma1, gamma2):
        # 1. Precompute distances for the entire arrays using the fast vectorised function
        ds_arr = angular_diameter_distance_numba(zs)
        dl_arr = angular_diameter_distance_numba(zl)
        
        # 2. Pass the arrays into the parallel evaluator
        return _cross_section_eval(
            zs, zl, ds_arr, dl_arr, sigma, q, phi, gamma, gamma1, gamma2,
            e1_grid_val, e2_grid_val, gamma_grid_val, gamma1_grid_val, gamma2_grid_val,
            cs_unit_grid_val,
            e1_c, e2_c, g_c, g1_c, g2_c,
            slope, intercept
        )

    return cross_section_reinit