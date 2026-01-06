"""
Cross Section Interpolator Module
==================================

This module provides highly optimized Numba-compiled functions for interpolating
gravitational lensing cross sections in a 5-dimensional parameter space.

The module uses cubic B-spline interpolation on prefiltered coefficient grids to
achieve efficient and accurate interpolation of cross section values based on:
  - Ellipticity components (e1, e2) derived from axis ratio q and position angle phi
  - Density profile slope (gamma)
  - External shear components (gamma1, gamma2)

Key Features:
-------------
- JIT-compiled with Numba for high performance
- Parallel processing support for batch evaluations
- Cubic B-spline interpolation matching scipy's map_coordinates behavior
- Automatic scaling by Einstein radius and affine calibration

Author: Hemanta Kumar Phurailatpam
"""

import numpy as np
from numba import njit, prange
from scipy import ndimage

# Import Numba-compatible ellipticity conversion function
from ler.lens_galaxy_population.jit_functions import phi_q2_ellipticity_hemanta
C_LIGHT = 299792.458


# ---------------------
# Cubic B-spline basis (order=3)
# ---------------------
@njit(inline="always")
def _bspline3(t):
    """
    Evaluate the centered cubic B-spline basis function B3(t).

    This is the standard cubic B-spline with compact support on [-2, 2].
    It provides C2 continuity (twice continuously differentiable) which is
    essential for smooth interpolation.

    Parameters
    ----------
    t : float
        Distance from the spline center. The function has compact support,
        meaning it returns 0 for |t| >= 2.

    Returns
    -------
    float
        The B-spline basis value at position t.
        - Maximum value of 2/3 occurs at t=0
        - Smoothly decays to 0 at |t|=2
        - Returns 0 for |t| >= 2

    Notes
    -----
    The cubic B-spline is defined piecewise:
    - For |t| < 1: (4 - 6t² + 3t³) / 6
    - For 1 <= |t| < 2: (2-|t|)³ / 6
    - For |t| >= 2: 0
    """
    at = abs(t)  # Use absolute value for symmetric spline
    if at < 1.0:
        # Central region: polynomial approximation
        return (4.0 - 6.0 * at * at + 3.0 * at * at * at) / 6.0
    elif at < 2.0:
        # Outer region: cubic decay
        u = 2.0 - at
        return (u * u * u) / 6.0
    else:
        # Beyond support: zero
        return 0.0


@njit(inline="always")
def _clamp_int(i, n):
    """
    Clamp an index to valid array bounds using 'nearest' boundary mode.

    This implements the 'nearest' boundary condition used in scipy's
    map_coordinates. Out-of-bounds indices are clamped to the nearest
    valid index rather than wrapping or raising an error.

    Parameters
    ----------
    i : int
        The index to clamp. Can be negative or exceed array bounds.
    n : int
        The size of the array dimension (valid indices are 0 to n-1).

    Returns
    -------
    int
        The clamped index in the range [0, n-1]:
        - If i < 0, returns 0
        - If i > n-1, returns n-1
        - Otherwise returns i unchanged
    """
    if i < 0:
        return 0  # Clamp to lower bound
    if i > n - 1:
        return n - 1  # Clamp to upper bound
    return i  # Within bounds, return as-is


@njit(inline="always")
def _physical_to_index(x, x0, x1, n):
    """
    Convert a physical coordinate to a fractional grid index.

    Maps a value from the physical coordinate space [x0, x1] to the
    corresponding fractional index in the discrete grid space [0, n-1].
    This is essential for interpolation as it determines which grid points
    contribute to the interpolated value.

    Parameters
    ----------
    x : float
        The physical coordinate value to convert.
    x0 : float
        The minimum physical coordinate (maps to index 0).
    x1 : float
        The maximum physical coordinate (maps to index n-1).
    n : int
        The number of grid points in this dimension.

    Returns
    -------
    float
        The fractional index in [0, n-1]. Non-integer values indicate
        positions between grid points that require interpolation.

    Examples
    --------
    >>> _physical_to_index(0.5, 0.0, 1.0, 11)
    5.0  # Exactly at grid point 5
    >>> _physical_to_index(0.55, 0.0, 1.0, 11)
    5.5  # Halfway between grid points 5 and 6
    """
    return (x - x0) * (n - 1) / (x1 - x0)


# ---------------------
# 5D cubic B-spline interpolation on prefiltered coefficients
# This matches scipy.ndimage.map_coordinates with:
#   order=3 (cubic), mode='nearest', prefilter=False
# ---------------------
@njit(parallel=True)
def _map_coordinates_5d_cubic_nearest(coeff, ie1, ie2, ig, ig1, ig2):
    """
    Perform 5D cubic B-spline interpolation on prefiltered coefficients.

    This function evaluates a cubic (order-3) B-spline interpolation in 5 dimensions.
    It matches the behavior of scipy.ndimage.map_coordinates with order=3,
    mode='nearest', and prefilter=False. The function is parallelized across
    the input points for optimal performance.

    Parameters
    ----------
    coeff : ndarray, shape (n1, n2, n3, n4, n5)
        Prefiltered B-spline coefficients. This should be the output of
        scipy.ndimage.spline_filter(data, order=3) applied to the original
        grid data.
    ie1, ie2, ig, ig1, ig2 : ndarray, shape (N,)
        Fractional indices for each of the 5 dimensions. These are the
        positions at which to interpolate, expressed in grid index coordinates.
        Values can be non-integer (for positions between grid points) and
        can be outside [0, n-1] (will use 'nearest' boundary mode).

    Returns
    -------
    out : ndarray, shape (N,)
        Interpolated values at the specified fractional indices.

    Algorithm
    ---------
    For each output point:
    1. Determine the 4x4x4x4x4 neighborhood of grid points that contribute
       to the interpolation (cubic B-spline has compact support of width 4)
    2. Compute the B-spline basis weight for each dimension
    3. Accumulate the weighted sum of coefficient values
    4. Skip zero-weight contributions for efficiency

    Notes
    -----
    - The function uses parallel loops (prange) for processing multiple points
    - Zero-weight contributions are skipped to improve performance
    - Boundary handling uses 'nearest' mode (clamping to valid indices)
    - This is significantly faster than scipy for large batches due to JIT compilation

    See Also
    --------
    scipy.ndimage.map_coordinates : The scipy function this mimics
    scipy.ndimage.spline_filter : For computing the prefiltered coefficients
    """
    n1, n2, n3, n4, n5 = coeff.shape
    N = ie1.size
    out = np.empty(N, dtype=np.float64)

    # Parallel loop over all interpolation points
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

# ---------------------
# User-facing jitted function with the signature you requested
# IMPORTANT: Always returns a 1D array (consistent type).
# ---------------------
@njit
def cross_section(
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
    theta_E, gamma, gamma1, gamma2, q, phi: 1D arrays (same length)
    grids: 1D
    cs_spline_coeff: spline-filtered coefficients (same shape as cs_spline_coeff)
    returns: 1D array
    """

    # angular diameter distance
    Ds = Da_instance(zs)
    Dl = Da_instance(zl)
    Dls = (Ds*(1+zs) - Dl*(1+zl))/(1+zs)
    # find theta_E
    theta_E = 4.0 * np.pi * (sigma / C_LIGHT) ** 2 * (Dls / Ds)

    # Convert (q, phi) -> (e1, e2)
    e1, e2 = phi_q2_ellipticity_hemanta(phi, q)

    N = zs.size

    # fractional indices
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

    # unit cross-section via cubic spline on coefficients
    cs_unit = _map_coordinates_5d_cubic_nearest(cs_spline_coeff_grid, ie1, ie2, ig, ig1, ig2)

    # scaling by area_sis and affine calibration
    # out = np.empty(N, dtype=np.float64)
    # for i in range(N):
    #     area_sis = np.pi * theta_E[i] * theta_E[i]
    #     cs = cs_unit[i] * (csunit_to_cs_intercept + csunit_to_cs_slope * area_sis)
    #     if cs < 0.0:
    #         cs = 0.0
    #     out[i] = cs
    area_sis = np.pi * theta_E**2
    cs = cs_unit * (csunit_to_cs_intercept + csunit_to_cs_slope * area_sis)
    idx = np.where(cs < 0.0)
    cs[idx] = 0.0

    return cs


# ---------------------
# Factory: computes coeffs once and returns a reinitialized njit closure
# ---------------------
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
    
    """
    e1g = np.asarray(e1_grid, dtype=np.float64)
    e2g = np.asarray(e2_grid, dtype=np.float64)
    gg = np.asarray(gamma_grid, dtype=np.float64)
    g1g = np.asarray(gamma1_grid, dtype=np.float64)
    g2g = np.asarray(gamma2_grid, dtype=np.float64)
    cs = np.asarray(cs_spline_coeff_grid, dtype=np.float64)

    # Same preprocessing as your class:
    cs_coeff = ndimage.spline_filter(cs, order=3)

    slope_given = float(csunit_to_cs_slope)
    intercept_given = float(csunit_to_cs_intercept)

    _Da_instance = Da_instance
    
    Da_function = njit(lambda z: _Da_instance(z))

    # Close over "given" parameters
    @njit
    def cross_section_reinit(zs, zl, sigma, q, phi, gamma, gamma1, gamma2):
        return cross_section(
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
