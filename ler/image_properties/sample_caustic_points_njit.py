"""
Symmetry-accelerated source-position sampling for EPL + shear caustics.

This version avoids constructing the full caustic boundary for sampling.
It asks cross_section_njit._helper_caustic_epl_shear for only the half-boundary
0 <= theta < pi, rotates only that half-boundary, builds a half-fan cumulative
area table, samples from that half, and then flips the sampled point to the
antipodal half with probability 1/2.

No parallel=True and no cache=True are used.
"""

import numpy as np
from numba import njit

from .cross_section_njit import _helper_caustic_epl_shear

PI = np.pi

# ------------------------------------------------------------------
# Half-symmetry caustic generator
# ------------------------------------------------------------------
@njit(cache=True, fastmath=True)
def caustic_points_epl_shear_half(
    theta_E,
    q,
    phi,
    gamma,
    gamma1,
    gamma2,
    num_th=500,
    maginf=-100.0,
    quad=False,
):
    """
    Return only the rotated half-boundary of the EPL+shear caustic.

    The returned points correspond to the first half of the full angular grid,
    0 <= theta < pi.  The other half is obtained exactly by central inversion,
    p(theta + pi) = -p(theta).  This is intended for area/sampling routines
    that can exploit central symmetry and therefore do not need the explicit
    full boundary.

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
        Number of angular samples. \n
        default: 500
    maginf : ``float``
        Magnification cut threshold. \n
        default: -100.0
    quad : ``bool``
        If True, return the quad (inner) caustic. If False, return the double (outer) caustic. \n
        default: False

    Returns
    -------
    pos_half : numpy.ndarray
        Shape (2, num_th//2).  Rotated source-plane caustic half-boundary.
    """
    h = num_th // 2
    theta_h = np.arange(h) * (2.0 * PI / num_th)

    cos_th = np.cos(theta_h)
    sin_th = np.sin(theta_h)
    cos_2th = np.cos(2.0 * theta_h)
    sin_2th = np.sin(2.0 * theta_h)

    pos_half = _helper_caustic_epl_shear(
        q, phi, gamma, gamma1, gamma2, theta_E,
        theta_h, cos_th, sin_th, cos_2th, sin_2th,
        maginf=maginf, quad=quad, half_only=True,
    )

    # Equivalent to _rotmat(-phi) @ pos_half, but avoids constructing a 2x2 matrix
    # and avoids matrix multiplication for a small 2 x h array.
    c = np.cos(phi)
    s = np.sin(phi)
    rotated = np.empty_like(pos_half)
    for i in range(h):
        x = pos_half[0, i]
        y = pos_half[1, i]
        rotated[0, i] = c * x - s * y
        rotated[1, i] = s * x + c * y

    return rotated

# ------------------------------------------------------------------
# Fan-triangulation sampler from origin (exact, zero-rejection)
# ------------------------------------------------------------------
@njit(cache=True, fastmath=True)
def _build_half_fan_cumulative(xh, yh):
    """
    Build cumulative fan-triangle areas for one centrally symmetric half.

    The half polygon is closed by the antipodal of the first point:
    [p0, p1, ..., p_{h-1}, -p0].  It contains exactly half of the full
    centrally symmetric caustic area.

    Parameters
    ----------
    xh : numpy.ndarray
        X-coordinates of the half-boundary points. \n
    yh : numpy.ndarray
        Y-coordinates of the half-boundary points. \n

    Returns
    -------
    cum_area : numpy.ndarray
        Cumulative areas of the fan-triangles. \n
    total : float
        Total area of the half-boundary polygon. \n
    """
    h = xh.shape[0]
    cum_area = np.empty(h + 1, dtype=np.float64)
    cum_area[0] = 0.0

    total = 0.0
    for i in range(h):
        if i + 1 < h:
            xj = xh[i + 1]
            yj = yh[i + 1]
        else:
            xj = -xh[0]
            yj = -yh[0]

        tri_area = 0.5 * abs(xh[i] * yj - xj * yh[i])
        total += tri_area
        cum_area[i + 1] = total

    return cum_area, total


@njit(cache=True, fastmath=True)
def _sample_one_half_fan(xh, yh, cum_area, half_area):
    """
    Uniformly sample one point from a centrally symmetric polygon using only
    the half-boundary fan table.

    First sample uniformly from one half of the fan triangulation, then choose
    the antipodal copy with probability 1/2.  This is exactly equivalent to
    sampling from the full fan triangulation because every triangle has an
    antipodal triangle with the same area.

    Parameters
    ----------
    xh : numpy.ndarray
        X-coordinates of the half-boundary points. \n
    yh : numpy.ndarray
        Y-coordinates of the half-boundary points. \n
    cum_area : numpy.ndarray
        Cumulative areas of the half-boundary fan triangulation. \n
    half_area : float
        Total area of the half-boundary polygon. The full polygon area is twice this. \n

    Returns
    -------
    x : float
        Sampled x-coordinate. \n
    y : float
        Sampled y-coordinate. \n
    """
    if half_area <= 0.0:
        return np.nan, np.nan

    h = xh.shape[0]

    u = np.random.random() * half_area
    i = np.searchsorted(cum_area[1:], u)

    if i + 1 < h:
        xj = xh[i + 1]
        yj = yh[i + 1]
    else:
        xj = -xh[0]
        yj = -yh[0]

    # Uniform sample in triangle (0, p_i, p_j)
    a = np.random.random()
    b = np.random.random()
    if a + b > 1.0:
        a = 1.0 - a
        b = 1.0 - b

    x = a * xh[i] + b * xj
    y = a * yh[i] + b * yj

    # Choose the antipodal half with probability 1/2.
    if np.random.random() < 0.5:
        x = -x
        y = -y

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
    num_th=500,
    maginf=-100.0,
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

    pts_half = caustic_points_epl_shear_half(
        theta_E, q, phi, gamma, gamma1, gamma2,
        num_th=num_th, maginf=maginf, quad=False,
    )

    cum_area, half_area = _build_half_fan_cumulative(pts_half[0], pts_half[1])
    if half_area <= 0.0:
        return np.nan, np.nan

    return _sample_one_half_fan(pts_half[0], pts_half[1], cum_area, half_area)


# @njit(cache=True, fastmath=True)
# def prepare_half_fan_sampler(
#     theta_E,
#     q,
#     phi,
#     gamma,
#     gamma1,
#     gamma2,
#     num_th=500,
#     maginf=-100.0,
# ):
#     """
#     Precompute the half-boundary and cumulative half-fan table for repeated
#     sampling from the same lens.
#     """
#     pts_half = caustic_points_epl_shear_half(
#         theta_E, q, phi, gamma, gamma1, gamma2,
#         num_th=num_th, maginf=maginf, quad=False,
#     )
#     cum_area, half_area = _build_half_fan_cumulative(pts_half[0], pts_half[1])
#     return pts_half[0], pts_half[1], cum_area, half_area


# @njit(cache=True, fastmath=True)
# def sample_many_from_prepared_half_fan(xh, yh, cum_area, half_area, size):
#     """
#     Draw many source positions from a precomputed half-fan sampler.
#     """
#     xs = np.empty(size)
#     ys = np.empty(size)
#     for i in range(size):
#         xs[i], ys[i] = _sample_one_half_fan(xh, yh, cum_area, half_area)
#     return xs, ys
