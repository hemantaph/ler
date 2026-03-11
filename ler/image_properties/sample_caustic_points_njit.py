"""
Module for sampling source positions from EPL + shear caustic regions.

Provides numba-accelerated routines for uniform rejection sampling inside
arbitrary polygons, with convenience wrappers to sample from the double
caustic of an EPL + external shear lens model.

Usage:
    Draw a single source from the double-caustic region:

    >>> from ler.image_properties.sample_caustic_points_njit import sample_source_from_double_caustic
    >>> beta, pts = sample_source_from_double_caustic(
    ...     q=0.8, phi=0.0, gamma=2.0, gamma1=0.03, gamma2=-0.01
    ... )

Copyright (C) 2026 Phurailatpam Hemantakumar. Distributed under MIT License.
"""

from .cross_section_njit import caustics_epl_shear, polygon_area
import numpy as np
from numba import njit

TWO_PI = 2.0 * np.pi


@njit(cache=True, fastmath=True)
def _poly_bbox(xv, yv):
    """
    Compute the axis-aligned bounding box of a polygon.

    Parameters
    ----------
    xv : ``numpy.ndarray``
        Polygon vertex x-coordinates.
    yv : ``numpy.ndarray``
        Polygon vertex y-coordinates.

    Returns
    -------
    xmin : ``float``
        Minimum x value.
    xmax : ``float``
        Maximum x value.
    ymin : ``float``
        Minimum y value.
    ymax : ``float``
        Maximum y value.
    """
    n = xv.size
    xmin = xv[0]
    xmax = xv[0]
    ymin = yv[0]
    ymax = yv[0]
    for i in range(1, n):
        x = xv[i]
        y = yv[i]
        if x < xmin:
            xmin = x
        elif x > xmax:
            xmax = x
        if y < ymin:
            ymin = y
        elif y > ymax:
            ymax = y
    return xmin, xmax, ymin, ymax


@njit(cache=True, fastmath=True)
def _build_edge_coeffs(xv, yv):
    """
    Build precomputed edge coefficients for point-in-polygon tests.

    Parameters
    ----------
    xv : ``numpy.ndarray``
        Polygon vertex x-coordinates.
    yv : ``numpy.ndarray``
        Polygon vertex y-coordinates.

    Returns
    -------
    xi : ``numpy.ndarray``
        Edge start x-coordinate for each edge.
    yi : ``numpy.ndarray``
        Edge start y-coordinate for each edge.
    yj : ``numpy.ndarray``
        Edge end y-coordinate for each edge.
    slope : ``numpy.ndarray``
        Precomputed slope ``(xj - xi)/(yj - yi)``. 
        Horizontal edges are encoded with 0.0 and skipped by the test.
    """
    n = xv.size
    xi = np.empty(n, dtype=np.float64)
    yi = np.empty(n, dtype=np.float64)
    yj = np.empty(n, dtype=np.float64)
    slope = np.empty(n, dtype=np.float64)

    # edge from i to j=i+1 (last to 0)
    for i in range(n - 1):
        x0 = xv[i]
        y0 = yv[i]
        x1 = xv[i + 1]
        y1 = yv[i + 1]
        xi[i] = x0
        yi[i] = y0
        yj[i] = y1
        dy = y1 - y0
        if abs(dy) < 1e-300:
            slope[i] = 0.0
        else:
            slope[i] = (x1 - x0) / dy

    # last edge
    x0 = xv[n - 1]
    y0 = yv[n - 1]
    x1 = xv[0]
    y1 = yv[0]
    xi[n - 1] = x0
    yi[n - 1] = y0
    yj[n - 1] = y1
    dy = y1 - y0
    if abs(dy) < 1e-300:
        slope[n - 1] = 0.0
    else:
        slope[n - 1] = (x1 - x0) / dy

    return xi, yi, yj, slope


@njit(cache=True)
def _point_in_poly_precomp(x, y, xi, yi, yj, slope):
    """
    Test whether a point is inside a polygon using the even-odd rule.

    Parameters
    ----------
    x : ``float``
        Point x-coordinate.
    y : ``float``
        Point y-coordinate.
    xi : ``numpy.ndarray``
        Edge start x-coordinates.
    yi : ``numpy.ndarray``
        Edge start y-coordinates.
    yj : ``numpy.ndarray``
        Edge end y-coordinates.
    slope : ``numpy.ndarray``
        Precomputed edge slopes.

    Returns
    -------
    inside : ``int``
        1 if the point is inside, otherwise 0.
    """
    inside = False
    n = xi.size
    for k in range(n):
        s = slope[k]
        if s == 0.0:
            continue
        y0 = yi[k]
        y1 = yj[k]
        # (y0 > y) != (y1 > y) is the standard crossing condition
        c0 = y0 > y
        c1 = y1 > y
        if c0 != c1:
            xint = xi[k] + (y - y0) * s
            if x < xint:
                inside = not inside
    return 1 if inside else 0


@njit(cache=True)
def _points_in_poly_precomp(xs, ys, xi, yi, yj, slope):
    """
    Evaluate point-in-polygon for multiple query points.

    Parameters
    ----------
    xs : ``numpy.ndarray``
        Query x-coordinates.
    ys : ``numpy.ndarray``
        Query y-coordinates.
    xi : ``numpy.ndarray``
        Edge start x-coordinates.
    yi : ``numpy.ndarray``
        Edge start y-coordinates.
    yj : ``numpy.ndarray``
        Edge end y-coordinates.
    slope : ``numpy.ndarray``
        Precomputed edge slopes.

    Returns
    -------
    out : ``numpy.ndarray``
        UInt8 mask of inside flags with 1 for inside and 0 for outside.
    """
    n_points = xs.size
    out = np.empty(n_points, dtype=np.uint8)
    for j in range(n_points):
        out[j] = _point_in_poly_precomp(xs[j], ys[j], xi, yi, yj, slope)
    return out


@njit(cache=True)
def sample_uniform_in_polygon(xv, yv, max_tries=100000):
    """
    Draw one uniform random point inside a polygon via rejection sampling.

    Parameters
    ----------
    xv : ``numpy.ndarray``
        Polygon vertex x-coordinates.
    yv : ``numpy.ndarray``
        Polygon vertex y-coordinates.
    max_tries : ``int``
        Maximum number of rejection-sampling attempts. 
        default: 100000

    Returns
    -------
    x : ``float``
        Sampled x-coordinate.
    y : ``float``
        Sampled y-coordinate.
    ok : ``int``
        Sampling flag where 1 indicates success and 0 indicates failure.

    Examples
    --------
    >>> import numpy as np
    >>> xv = np.array([0.0, 1.0, 1.0, 0.0])
    >>> yv = np.array([0.0, 0.0, 1.0, 1.0])
    >>> x, y, ok = sample_uniform_in_polygon(xv, yv)
    """
    xmin, xmax, ymin, ymax = _poly_bbox(xv, yv)
    dx = xmax - xmin
    dy = ymax - ymin
    if dx <= 0.0 or dy <= 0.0:
        return np.nan, np.nan, 0

    xi, yi, yj, slope = _build_edge_coeffs(xv, yv)

    for _ in range(max_tries):
        x = xmin + dx * np.random.random()
        y = ymin + dy * np.random.random()
        if _point_in_poly_precomp(x, y, xi, yi, yj, slope) == 1:
            return x, y, 1

    return np.nan, np.nan, 0


@njit(cache=True)
def sample_many_uniform_in_polygon(xv, yv, n_samples, max_tries_per=100000):
    """
    Draw multiple uniform random points inside a polygon.

    Sampling uses batched rejection from the polygon bounding box.

    Parameters
    ----------
    xv : ``numpy.ndarray``
        Polygon vertex x-coordinates.
    yv : ``numpy.ndarray``
        Polygon vertex y-coordinates.
    n_samples : ``int``
        Number of requested samples.
    max_tries_per : ``int``
        Maximum attempts per requested sample. 
        default: 100000

    Returns
    -------
    x_src : ``numpy.ndarray``
        Sampled x-coordinates with capacity ``n_samples``.
    y_src : ``numpy.ndarray``
        Sampled y-coordinates with capacity ``n_samples``.
    n_ok : ``int``
        Number of valid samples in the leading entries.

    Examples
    --------
    >>> import numpy as np
    >>> xv = np.array([0.0, 1.0, 1.0, 0.0])
    >>> yv = np.array([0.0, 0.0, 1.0, 1.0])
    >>> x_src, y_src, n_ok = sample_many_uniform_in_polygon(xv, yv, 100)
    """
    xmin, xmax, ymin, ymax = _poly_bbox(xv, yv)
    dx = xmax - xmin
    dy = ymax - ymin
    if dx <= 0.0 or dy <= 0.0:
        return np.empty(0), np.empty(0), 0

    area = polygon_area(xv, yv)
    bbox_area = dx * dy
    if area <= 0.0 or bbox_area <= 0.0:
        return np.empty(0), np.empty(0), 0

    eff = area / bbox_area
    if eff <= 0.0:
        return np.empty(0), np.empty(0), 0

    xi, yi, yj, slope = _build_edge_coeffs(xv, yv)

    x_src = np.empty(n_samples, dtype=np.float64)
    y_src = np.empty(n_samples, dtype=np.float64)
    n_ok = 0

    # choose a batch size that amortizes overhead but doesn't explode memory
    batch_size = 65536
    max_total = n_samples * max_tries_per
    total_tries = 0

    while n_ok < n_samples and total_tries < max_total:
        n_need = n_samples - n_ok

        # expected candidates to get n_need accepts ~ n_need/eff
        n_gen = int(np.ceil(n_need / eff * 1.1))
        if n_gen > batch_size:
            n_gen = batch_size
        # also respect max_total
        remain = max_total - total_tries
        if n_gen > remain:
            n_gen = remain
        if n_gen <= 0:
            break

        x_cands = xmin + dx * np.random.random(n_gen)
        y_cands = ymin + dy * np.random.random(n_gen)

        inside = _points_in_poly_precomp(x_cands, y_cands, xi, yi, yj, slope)

        # gather without np.where (cheaper in njit)
        for j in range(n_gen):
            if inside[j] == 1:
                x_src[n_ok] = x_cands[j]
                y_src[n_ok] = y_cands[j]
                n_ok += 1
                if n_ok == n_samples:
                    break

        total_tries += n_gen

    return x_src, y_src, n_ok


@njit(cache=True)
def sample_source_from_double_caustic(
    q,
    phi,
    gamma,
    gamma1,
    gamma2,
    theta_E=1.0,
    num_th=500,
    maginf=-100.0,
    max_tries=10,
):
    """
    Draw one source position uniformly inside the double caustic.

    Parameters
    ----------
    q : ``float``
        Lens axis ratio.
    phi : ``float``
        Lens position angle (rad).
    gamma : ``float``
        EPL slope parameter.
    gamma1 : ``float``
        External shear component 1.
    gamma2 : ``float``
        External shear component 2.
    theta_E : ``float``
        Einstein radius used to construct the caustic. 
        default: 1.0
    num_th : ``int``
        Number of angular samples for the boundary curve. 
        default: 500
    maginf : ``float``
        Magnification cutoff used in the caustic construction. 
        default: -100.0
    max_tries : ``int``
        Maximum rejection-sampling attempts. 
        default: 100000

    Returns
    -------
    beta_x : ``float``
        Sampled source x-coordinate.
    beta_y : ``float``
        Sampled source y-coordinate.
    ok : ``int``
        Sampling status flag where 1 indicates success.

    Examples
    --------
    >>> beta_x, beta_y, ok = sample_source_from_double_caustic(
    ...     q=0.8, phi=0.0, gamma=2.0, gamma1=0.03, gamma2=-0.01
    ... )
    """

    pts = caustics_epl_shear(
        q,
        phi,
        gamma,
        gamma1,
        gamma2,
        theta_E=theta_E,
        num_th=num_th,
        maginf=maginf,
        sourceplane=True,
    )
    return sample_uniform_in_polygon(pts[0], pts[1], max_tries=max_tries), pts


@njit(cache=True)
def sample_many_sources_from_double_caustic(
    q,
    phi,
    gamma,
    gamma1,
    gamma2,
    n_samples,
    theta_E=1.0,
    num_th=500,
    maginf=-100.0,
    max_tries_per=10,
):
    """
    Draw many source positions uniformly inside the double caustic.

    Parameters
    ----------
    q : ``float``
        Lens axis ratio.
    phi : ``float``
        Lens position angle (rad).
    gamma : ``float``
        EPL slope parameter.
    gamma1 : ``float``
        External shear component 1.
    gamma2 : ``float``
        External shear component 2.
    n_samples : ``int``
        Number of source samples requested.
    theta_E : ``float``
        Einstein radius used to construct the caustic. 
        default: 1.0
    num_th : ``int``
        Number of angular samples for the boundary curve. 
        default: 500
    maginf : ``float``
        Magnification cutoff used in the caustic construction. 
        default: -100.0
    max_tries_per : ``int``
        Maximum rejection attempts per requested sample. 
        default: 100000

    Returns
    -------
    beta_x : ``numpy.ndarray``
        Sampled source x-coordinates.
    beta_y : ``numpy.ndarray``
        Sampled source y-coordinates.
    n_ok : ``int``
        Number of valid samples in the leading entries.

    Examples
    --------
    >>> bx, by, n_ok = sample_many_sources_from_double_caustic(
    ...     q=0.8, phi=0.0, gamma=2.0, gamma1=0.03, gamma2=-0.01, n_samples=1000
    ... )
    """

    pts = caustics_epl_shear(
        q,
        phi,
        gamma,
        gamma1,
        gamma2,
        theta_E=theta_E,
        num_th=num_th,
        maginf=maginf,
        sourceplane=True,
    )
    return sample_many_uniform_in_polygon(
        pts[0], pts[1], n_samples, max_tries_per=max_tries_per
    )