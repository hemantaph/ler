"""
Module for sampling source positions from EPL + shear caustic regions.

Provides numba-accelerated routines for exact fan-triangulation sampling
inside arbitrary star-convex polygons, with a convenience wrapper to
sample from the double caustic of an EPL + external shear lens model.

Usage:
    Basic workflow example:

    >>> from ler.image_properties.sample_caustic_points_njit import sample_source_from_double_caustic
    >>> xs, ys = sample_source_from_double_caustic(
    ...     theta_E=1.0, q=0.8, phi=0.0, gamma=2.0, gamma1=0.03, gamma2=-0.01
    ... )

Copyright (C) 2026 Phurailatpam Hemantakumar. Distributed under MIT License.
"""

from .cross_section_njit import caustic_points_epl_shear
import numpy as np
from numba import njit

# ------------------------------------------------------------------
# Fan-triangulation sampler from origin (exact, zero-rejection)
# ------------------------------------------------------------------
@njit(cache=True, fastmath=True)
def _build_fan_cumulative(xv, yv):
    """
    Precompute cumulative triangle areas for fan triangulation from (0,0).
    Assumes the polygon is simple and contains the origin
    (true for EPL+shear caustics).

    Parameters
    ----------
    xv : ``numpy.ndarray``
        x-coordinates of the polygon vertices. \n
    yv : ``numpy.ndarray``
        y-coordinates of the polygon vertices. \n

    Returns
    -------
    cum_area : ``numpy.ndarray``
        Cumulative triangle areas of length ``n + 1``. \n
    total : ``float``
        Total polygon area. \n
    """
    n = xv.shape[0]
    cum_area = np.empty(n + 1, dtype=np.float64)
    cum_area[0] = 0.0
    total = 0.0
    for i in range(n):
        j = (i + 1) % n
        tri_area = 0.5 * np.abs(xv[i] * yv[j] - xv[j] * yv[i])
        total += tri_area
        cum_area[i + 1] = total
    return cum_area, total


@njit(cache=True, fastmath=True)
def _sample_one_fan(xv, yv, cum_area, total_area):
    """
    Sample a single point uniformly inside a fan-triangulated polygon.

    Parameters
    ----------
    xv : ``numpy.ndarray``
        x-coordinates of the polygon vertices. \n
    yv : ``numpy.ndarray``
        y-coordinates of the polygon vertices. \n
    cum_area : ``numpy.ndarray``
        Cumulative triangle areas from :func:`_build_fan_cumulative`. \n
    total_area : ``float``
        Total polygon area. \n

    Returns
    -------
    x : ``float``
        Sampled x-coordinate (``NaN`` if area is zero). \n
    y : ``float``
        Sampled y-coordinate (``NaN`` if area is zero). \n
    """
    if total_area <= 0.0:
        return np.nan, np.nan

    u = np.random.random() * total_area
    # Binary search for the triangle (O(log N) with N=500)
    i = np.searchsorted(cum_area[1:], u)
    j = (i + 1) % xv.shape[0]

    # Uniform sample inside triangle (0, p_i, p_j)
    a = np.random.random()
    b = np.random.random()
    if a + b > 1.0:
        a = 1.0 - a
        b = 1.0 - b

    x = a * xv[i] + b * xv[j]
    y = a * yv[i] + b * yv[j]
    return x, y


# ------------------------------------------------------------------
# Factory: precompute caustic once → ultra-fast batched sampling
# ------------------------------------------------------------------
@njit(cache=True, fastmath=True)
def sample_source_from_double_caustic(
    theta_E,
    q,
    phi,
    gamma,
    gamma1,
    gamma2,
    num_th = 500,
    maginf = -100.0,
):
    """
    Sample a single source position from the double caustic region.

    Precomputes the caustic boundary once (the expensive part) and
    draws one uniform sample via exact fan-triangulation from the origin.

    Parameters
    ----------
    theta_E : ``float``
        Einstein radius. \n
    q : ``float``
        Axis ratio. \n
    phi : ``float``
        Position angle of the lens major axis in radians. \n
    gamma : ``float``
        EPL slope exponent. \n
    gamma1 : ``float``
        First shear component. \n
    gamma2 : ``float``
        Second shear component. \n
    num_th : ``int``
        Number of angular samples for caustic boundary. \n
        default: 500
    maginf : ``float``
        Magnification cut threshold. \n
        default: -100.0

    Returns
    -------
    xs : ``float``
        Sampled source x-coordinate (``NaN`` if caustic is invalid). \n
    ys : ``float``
        Sampled source y-coordinate (``NaN`` if caustic is invalid). \n

    Examples
    --------
    >>> xs, ys = sample_source_from_double_caustic(
    ...     theta_E=1.0, q=0.8, phi=0.0, gamma=2.0, gamma1=0.03, gamma2=-0.01
    ... )
    """

    pts = caustic_points_epl_shear(
        theta_E, q, phi, gamma, gamma1, gamma2, num_th=num_th, maginf=maginf
    )

    cum_area, total_area = _build_fan_cumulative(pts[0], pts[1])

    if total_area <= 0.0:
        return np.nan, np.nan  # invalid caustic, return a dummy sampler that fails

    xs, ys = _sample_one_fan(pts[0], pts[1], cum_area, total_area)

    return xs, ys