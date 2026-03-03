# sample_caustic_points_njit.py
from .cross_section_njit import caustic_double_points_epl_shear, polygon_area
import numpy as np
from numba import njit, prange

TWO_PI = 2.0 * np.pi


# ---------------------------------------------------------------------------
# Geometry helpers
# ---------------------------------------------------------------------------


@njit(cache=True)
def _poly_bbox(xv, yv):
    """Compute polygon bounding box (manual reduction; fast in njit)."""
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


@njit(cache=True)
def _build_edge_coeffs(xv, yv):
    """
    Precompute edge coefficients for the even–odd (ray crossing) point-in-polygon test.

    For each edge (i -> j), store:
      xi, yi, yj, slope = (xj-xi)/(yj-yi)

    Horizontal edges (yj == yi) are marked with slope=0 and will be skipped by the test.
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


@njit(cache=True, inline="always")
def _point_in_poly_precomp(x, y, xi, yi, yj, slope):
    """
    Even–odd rule with precomputed edge coefficients.

    Skips horizontal edges via slope==0 sentinel.
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


@njit(parallel=True, cache=True)
def _points_in_poly_precomp(xs, ys, xi, yi, yj, slope):
    """Parallel PIP for many points, using precomputed edge coefficients."""
    n_points = xs.size
    out = np.empty(n_points, dtype=np.uint8)
    for j in prange(n_points):
        out[j] = _point_in_poly_precomp(xs[j], ys[j], xi, yi, yj, slope)
    return out


# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------


@njit(cache=True)
def sample_uniform_in_polygon(xv, yv, max_tries=100000):
    """
    Sample one point uniformly inside a simple polygon via bbox rejection.
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
    Sample multiple points uniformly inside a polygon via batched rejection.
    Returns arrays (x_src, y_src, n_ok) where only first n_ok entries are valid.
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


# ---------------------------------------------------------------------------
# Public API — sampling from the double caustic
# ---------------------------------------------------------------------------


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
    max_tries=200000,
):
    """
    Sample one source position uniformly inside the double caustic polygon.
    """

    pts = caustic_double_points_epl_shear(
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
    return sample_uniform_in_polygon(pts[0], pts[1], max_tries=max_tries)


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
    max_tries_per=200000,
):
    """
    Sample multiple source positions uniformly inside the double caustic polygon.
    """

    pts = caustic_double_points_epl_shear(
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
