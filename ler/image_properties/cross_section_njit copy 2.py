"""
Module for analytical caustic computation and cross-section evaluation
of EPL (Elliptical Power-Law) + external shear lens models.

Provides numba-accelerated routines for computing critical curves, caustic
boundaries, polygon areas, and lensing cross sections. All core functions
are decorated with ``@njit`` for high performance.

Usage:
    Compute the double-caustic boundary for a single lens:

    >>> from ler.image_properties.cross_section_njit import caustics_epl_shear
    >>> pts = caustics_epl_shear(q=0.8, phi=0.0, gamma=2.0, gamma1=0.03, gamma2=-0.01)

Copyright (C) 2026 Phurailatpam Hemantakumar. Distributed under MIT License.
"""

import numpy as np
from numba import njit, prange

C_LIGHT = 299792458.0  # m/s
PI = np.pi
TWO_PI = 2.0 * PI

# === ORIGINAL PUBLIC API (For backwards compatibility with other modules) ===
@njit(cache=True, fastmath=True, inline="always")
def omega_scalar(phi, t, q, niter_max=200, tol=1e-16):
    f = (1.0 - q) / (1.0 + q)
    if f <= 0.0:
        return np.cos(phi), np.sin(phi)

    niter = min(niter_max, int(np.log(tol) / np.log(f)) + 2)
    if niter < 1: niter = 1

    c1, s1 = np.cos(phi), np.sin(phi)
    c2, s2 = np.cos(2.0 * phi), np.sin(2.0 * phi)

    ur, ui = c1, s1
    fr, fi = -f * c2, -f * s2

    res_r, res_i = 0.0, 0.0

    for n in range(1, niter):
        res_r += ur
        res_i += ui

        an = (2.0 * n - (2.0 - t)) / (2.0 * n + (2.0 - t))
        tmp_r = ur * fr - ui * fi
        tmp_i = ur * fi + ui * fr

        ur, ui = an * tmp_r, an * tmp_i

    res_r += ur
    res_i += ui
    return res_r, res_i

# === INTERNAL FAST LOOP OPTIMIZATION ===
@njit(cache=True, fastmath=True, inline="always")
def _omega_scalar_opt(phi, t, f, niter):
    """Optimized internal omega: f and niter are pre-computed once per system."""
    if f <= 0.0:
        return np.cos(phi), np.sin(phi)

    c1, s1 = np.cos(phi), np.sin(phi)
    c2, s2 = np.cos(2.0 * phi), np.sin(2.0 * phi)

    ur, ui = c1, s1
    fr, fi = -f * c2, -f * s2

    res_r, res_i = 0.0, 0.0

    for n in range(1, niter):
        res_r += ur
        res_i += ui

        an = (2.0 * n - (2.0 - t)) / (2.0 * n + (2.0 - t))
        tmp_r = ur * fr - ui * fi
        tmp_i = ur * fi + ui * fr

        ur, ui = an * tmp_r, an * tmp_i

    res_r += ur
    res_i += ui
    return res_r, res_i

@njit(cache=True, fastmath=True, inline="always")
def solve_quad_scalar_real(a, b, c):
    if b != 0.0:
        if a != 0.0:
            disc = b * b - 4.0 * a * c
            if disc < 0.0:
                disc = 0.0

            sD = np.sqrt(disc)
            sign_b = 1.0 if b > 0.0 else -1.0
            denominator = -b - sign_b * sD

            if np.abs(denominator) > 1e-30:
                x1 = denominator / (2.0 * a)
                x2 = (2.0 * c) / denominator
            else:
                x1 = (-b + sD) / (2.0 * a)
                x2 = x1
        else:
            x1, x2 = -c / b, -c / b
    else:
        if a != 0.0:
            val = -c / a
            if val < 0.0:
                val = 0.0
            x2 = np.sqrt(val)
            x1 = -x2
        else:
            x1, x2 = 0.0, 0.0
    return x1, x2

# === IN-PLACE SORTING (Zero memory allocation) ===
@njit(cache=True, fastmath=True, inline="always")
def partition_polar(th, r, low, high):
    pivot = th[high]
    i = low - 1
    for j in range(low, high):
        if th[j] <= pivot:
            i += 1
            th[i], th[j] = th[j], th[i]
            r[i], r[j] = r[j], r[i]
    th[i + 1], th[high] = th[high], th[i + 1]
    r[i + 1], r[high] = r[high], r[i + 1]
    return i + 1

@njit(cache=True, fastmath=True)
def quicksort_polar(th, r, low, high):
    if low < high:
        pi = partition_polar(th, r, low, high)
        quicksort_polar(th, r, low, pi - 1)
        quicksort_polar(th, r, pi + 1, high)

@njit(cache=True, fastmath=True, inline="always")
def _interp_periodic_inplace(x, xp, fp, period, out):
    n = xp.shape[0]
    for i in range(x.shape[0]):
        xi = x[i] % period
        lo, hi = 0, n - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if xp[mid] < xi:
                lo = mid + 1
            else:
                hi = mid

        if lo == 0:
            x0, x1 = xp[n - 1] - period, xp[0]
            f0, f1 = fp[n - 1], fp[0]
        else:
            x0, x1 = xp[lo - 1], xp[lo]
            f0, f1 = fp[lo - 1], fp[lo]

        dx = x1 - x0
        if abs(dx) < 1e-30:
            out[i] = f0
        else:
            out[i] = f0 + (f1 - f0) * (xi - x0) / dx

# 1. THE HELPER 
@njit(cache=True, fastmath=True, inline="always")
def _helper_caustic_epl_shear(
    q, t, gamma1, gamma2, b,
    num_th, cos_th, sin_th, cos_2th, sin_2th,
    rcut, thcut, r_main, th_main, r2,
    maginf=-100.0
):
    """
    Generates the unrotated radial boundaries of the EPL+shear caustic 
    and writes them into the provided workspace arrays.
    """
    abs_gamma_sq = gamma1**2 + gamma2**2
    one_minus_t = 1.0 - t
    one_minus_t_sq = one_minus_t**2
    two_minus_t = 2.0 - t
    two_by_one_plus_q = 2.0 / (1.0 + q)
    minus_one_by_t = -1.0 / t
    aa = 1.0 - abs_gamma_sq

    f_omega = (1.0 - q) / (1.0 + q)
    niter_omega = 1
    if f_omega > 0.0:
        niter_omega = min(200, int(np.log(1e-16) / np.log(f_omega)) + 2)
        if niter_omega < 1: niter_omega = 1

    if t > 1.0:
        aa_cut = aa - maginf
        use_x2_cut = True
    else:
        aa_cut = aa + maginf
        use_x2_cut = False

    for j in range(num_th):
        c_th, s_th = cos_th[j], sin_th[j]
        c_2th, s_2th = cos_2th[j], sin_2th[j]

        R_val = np.sqrt(q**2 * c_th**2 + s_th**2)
        phi_ell = np.arctan2(s_th, c_th * q)

        omega_real, omega_imag = _omega_scalar_opt(phi_ell, t, f_omega, niter_omega)
        frac_roverR = 1.0 / R_val

        cc = one_minus_t * two_minus_t * (c_th * omega_real + s_th * omega_imag) * R_val * two_by_one_plus_q
        cc -= one_minus_t_sq * (two_by_one_plus_q**2) * (omega_real**2 + omega_imag**2) * (R_val * R_val)

        shear_prefac = one_minus_t * two_by_one_plus_q * R_val
        gammaint_fac_real = -c_2th * two_minus_t / 2.0 + shear_prefac * (c_th * omega_real - s_th * omega_imag)
        gammaint_fac_imag = -s_2th * two_minus_t / 2.0 + shear_prefac * (c_th * omega_imag + s_th * omega_real)

        bb = -two_minus_t - 2.0 * (gamma1 * gammaint_fac_real + gamma2 * gammaint_fac_imag)

        # Main curve
        x1_4, x2_4 = solve_quad_scalar_real(cc, bb, aa)
        u4 = np.abs(x2_4)
        if u4 < 1e-15: u4 = 1e-15
        pow_u4 = u4 ** minus_one_by_t         
        r_4 = b * pow_u4 * frac_roverR
        R4 = b * pow_u4                       

        xcr_4, ycr_4 = r_4 * c_th, r_4 * s_th

        # Cut curve 
        x1_cut, x2_cut = solve_quad_scalar_real(cc, bb, aa_cut)
        u_cut = np.abs(x2_cut if use_x2_cut else x1_cut)
        if u_cut < 1e-15: u_cut = 1e-15
        pow_ucut = u_cut ** minus_one_by_t
        r_cut_val = b * pow_ucut * frac_roverR
        R_cut = b * pow_ucut

        xcr_cut, ycr_cut = r_cut_val * c_th, r_cut_val * s_th

        # Inline Alpha Main
        if R4 == 0.0:
            al_4_r = al_4_i = 0.0
        else:
            prefac = 2.0 / (1.0 + q) * u4 * R4
            al_4_r, al_4_i = prefac * omega_real, prefac * omega_imag
        al_4_r += gamma1 * xcr_4 + gamma2 * ycr_4
        al_4_i += gamma2 * xcr_4 - gamma1 * ycr_4

        # Inline Alpha Cut
        if R_cut == 0.0:
            al_cut_r = al_cut_i = 0.0
        else:
            prefac_cut = 2.0 / (1.0 + q) * u_cut * R_cut
            al_cut_r, al_cut_i = prefac_cut * omega_real, prefac_cut * omega_imag
        al_cut_r += gamma1 * xcr_cut + gamma2 * ycr_cut
        al_cut_i += gamma2 * xcr_cut - gamma1 * ycr_cut

        # Source-plane mapping
        xca_4, yca_4 = xcr_4 - al_4_r, ycr_4 - al_4_i
        xca_cut, yca_cut = xcr_cut - al_cut_r, ycr_cut - al_cut_i

        rcut[j] = np.sqrt(xca_cut**2 + yca_cut**2)
        thcut[j] = np.arctan2(yca_cut, xca_cut) % TWO_PI
        r_main[j] = np.sqrt(xca_4**2 + yca_4**2)
        th_main[j] = np.arctan2(yca_4, xca_4) % TWO_PI

    idx = np.argsort(thcut)
    thcut[:] = thcut[idx]
    rcut[:] = rcut[idx]

    _interp_periodic_inplace(th_main, thcut, rcut, TWO_PI, r2)

# 2. THE AREA CALCULATOR
@njit(cache=True, fastmath=True)
def caustic_area_epl_shear(
    q, t, gamma1, gamma2, b,
    num_th, cos_th, sin_th, cos_2th, sin_2th,
    rcut, thcut, r_main, th_main, r2,
    maginf=-100.0
):
    """
    Calculates only the area, skipping rotation and coordinate array construction.
    """
    _helper_caustic_epl_shear(
        q, t, gamma1, gamma2, b,
        num_th, cos_th, sin_th, cos_2th, sin_2th,
        rcut, thcut, r_main, th_main, r2,
        maginf
    )

    area = 0.0
    r_f_prev = r_main[num_th - 1] if r_main[num_th - 1] > r2[num_th - 1] else r2[num_th - 1]
    x_prev = r_f_prev * np.cos(th_main[num_th - 1])
    y_prev = r_f_prev * np.sin(th_main[num_th - 1])

    for i in range(num_th):
        r_f_curr = r_main[i] if r_main[i] > r2[i] else r2[i]
        x_curr = r_f_curr * np.cos(th_main[i])
        y_curr = r_f_curr * np.sin(th_main[i])

        area += x_prev * y_curr - x_curr * y_prev
        x_prev, y_prev = x_curr, y_curr

    return abs(area) / 2.0

# 3. THE POINT CALCULATOR
@njit(cache=True, fastmath=True)
def caustic_points_epl_shear(
    theta_E, q, phi, gamma, gamma1, gamma2, num_th=500, maginf=-100.0
):
    """
    Calculates the 2D coordinates of the double caustic for a SINGLE lens.
    Accepts scalar float values.
    """
    theta = np.linspace(0, 2 * PI * (1.0 - 1.0 / num_th), num_th)
    
    # Pre-compute trig
    cos_th, sin_th = np.cos(theta), np.sin(theta)
    cos_2th, sin_2th = np.cos(2.0 * theta), np.sin(2.0 * theta)

    # 1D Allocations for a single lens (Fast when not in a prange loop)
    rcut = np.empty(num_th)
    thcut = np.empty(num_th)
    r_main = np.empty(num_th)
    th_main = np.empty(num_th)
    r2 = np.empty(num_th)

    # Shear rotation
    th_gamma = np.arctan2(gamma2, gamma1) / 2.0
    g_mag = np.sqrt(gamma1**2 + gamma2**2)
    th_gamma -= phi
    g1_rot = g_mag * np.cos(2.0 * th_gamma)
    g2_rot = g_mag * np.sin(2.0 * th_gamma)

    b_val = np.sqrt(q) * theta_E
    t_val = gamma - 1.0
    
    # Call the helper
    _helper_caustic_epl_shear(
        q, t_val, g1_rot, g2_rot, b_val,
        num_th, cos_th, sin_th, cos_2th, sin_2th,
        rcut, thcut, r_main, th_main, r2,
        maginf
    )
    
    # Rotate final coordinates back to original frame
    mphi = -phi
    cos_mphi = np.cos(mphi)
    sin_mphi = np.sin(mphi)

    pts_out = np.empty((2, num_th))
    for j in range(num_th):
        r_f = r_main[j] if r_main[j] > r2[j] else r2[j]
        pos_x = r_f * np.cos(th_main[j])
        pos_y = r_f * np.sin(th_main[j])

        pts_out[0, j] = cos_mphi * pos_x + sin_mphi * pos_y
        pts_out[1, j] = -sin_mphi * pos_x + cos_mphi * pos_y

    return pts_out


def make_cross_section_area_reinit(Da_instance, num_th=500, maginf=-100.0):
    """Drop-in replacement with:
       - trig arrays closed over (faster)
       - safe memory layout
       - optional chunking for huge batches (no thread-ID hack)"""
    theta = np.linspace(0, 2 * np.pi * (1.0 - 1.0 / num_th), num_th)
    
    # Closed-over constants (no recomputation)
    cos_th = np.cos(theta)
    sin_th = np.sin(theta)
    cos_2th = np.cos(2.0 * theta)
    sin_2th = np.sin(2.0 * theta)

    @njit(parallel=True, cache=True, fastmath=True)
    def cross_section_area(zs, zl, sigma, q, phi, gamma, gamma1, gamma2):
        size = zs.size
        cs_area = np.empty(size)

        Ds = Da_instance(zs)
        Dl = Da_instance(zl)

        # Original safe layout - still the fastest and simplest
        rcut_work  = np.empty((size, num_th))
        thcut_work = np.empty((size, num_th))
        rmain_work = np.empty((size, num_th))
        thmain_work = np.empty((size, num_th))
        r2_work    = np.empty((size, num_th))

        for i in prange(size):
            Dls = (Ds[i] * (1.0 + zs[i]) - Dl[i] * (1.0 + zl[i])) / (1.0 + zs[i])
            theta_E = 4.0 * PI * (sigma[i] * 1000.0 / C_LIGHT) ** 2 * (Dls / Ds[i])

            th_gamma = np.arctan2(gamma2[i], gamma1[i]) / 2.0
            g_mag = np.sqrt(gamma1[i]**2 + gamma2[i]**2)
            th_gamma -= phi[i]
            g1_rot = g_mag * np.cos(2.0 * th_gamma)
            g2_rot = g_mag * np.sin(2.0 * th_gamma)

            b_val = np.sqrt(q[i]) * theta_E
            t_val = gamma[i] - 1.0

            area = caustic_area_epl_shear(
                q[i], t_val, g1_rot, g2_rot, b_val,
                num_th, cos_th, sin_th, cos_2th, sin_2th,
                rcut_work[i], thcut_work[i], rmain_work[i], thmain_work[i], r2_work[i],
                maginf
            )

            cs_area[i] = area if np.isfinite(area) else 0.0

        return cs_area

    return cross_section_area

@njit(cache=True, fastmath=True)
def pol_to_ell(r, theta, q):
    phi = np.arctan2(np.sin(theta), np.cos(theta) * q)
    rell = r * np.sqrt(q**2 * np.cos(theta) ** 2 + np.sin(theta) ** 2)
    return rell, phi

@njit(cache=True, fastmath=True, inline="always")
def select_physical_root(u0, u1, prefer_second):
    tiny = 1e-15
    imag_tol = 1e-10

    u0r = u0.real
    u1r = u1.real

    tol0 = imag_tol * (1.0 + np.abs(u0r))
    tol1 = imag_tol * (1.0 + np.abs(u1r))

    valid0 = np.isfinite(u0r) and (u0r > tiny) and (np.abs(u0.imag) <= tol0)
    valid1 = np.isfinite(u1r) and (u1r > tiny) and (np.abs(u1.imag) <= tol1)

    if valid0 and valid1:
        return u1r if prefer_second else u0r
    if valid0:
        return u0r
    if valid1:
        return u1r

    if prefer_second:
        if np.isfinite(u1r) and (u1r > tiny):
            return u1r
        if np.isfinite(u0r) and (u0r > tiny):
            return u0r
    else:
        if np.isfinite(u0r) and (u0r > tiny):
            return u0r
        if np.isfinite(u1r) and (u1r > tiny):
            return u1r
    return tiny

# -------------------------------------------
@njit(cache=True, fastmath=True)
def cdot(a, b):
    """
    Compute the real-valued dot product of two complex numbers.
    Used in epl_shear_njit.py.

    Equivalent to ``Re(a) * Re(b) + Im(a) * Im(b)``.

    Parameters
    ----------
    a : ``complex``
        First complex number.
    b : ``complex``
        Second complex number.

    Returns
    -------
    result : ``float``
        Real-valued dot product.

    Examples
    --------
    >>> cdot(1+2j, 3+4j)
    11.0
    """
    return a.real * b.real + a.imag * b.imag

@njit(cache=True)
def _interp_periodic(x, xp, fp, period):
    """
    Perform numba-compatible periodic linear interpolation.

    Equivalent to ``np.interp`` with the ``period`` keyword.

    Parameters
    ----------
    x : ``numpy.ndarray``
        Query points.
    xp : ``numpy.ndarray``
        Known x-values, sorted and within ``[0, period)``.
    fp : ``numpy.ndarray``
        Known function values corresponding to ``xp``.
    period : ``float``
        Period of the function.

    Returns
    -------
    result : ``numpy.ndarray``
        Interpolated values at query points.
    """
    n = xp.shape[0]
    result = np.empty_like(x)
    for i in range(x.shape[0]):
        xi = x[i] % period
        # binary search in sorted xp (assumed sorted and within [0, period))
        lo = 0
        hi = n - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if xp[mid] < xi:
                lo = mid + 1
            else:
                hi = mid
        # lo is the index of the first element >= xi
        if lo == 0:
            # xi is before the first point – wrap around
            x0 = xp[n - 1] - period
            x1 = xp[0]
            f0 = fp[n - 1]
            f1 = fp[0]
        else:
            x0 = xp[lo - 1]
            x1 = xp[lo]
            f0 = fp[lo - 1]
            f1 = fp[lo]
        dx = x1 - x0
        if abs(dx) < 1e-30:
            result[i] = f0
        else:
            result[i] = f0 + (f1 - f0) * (xi - x0) / dx
    return result

@njit(parallel=True, cache=True, fastmath=True)
def cross_section_epl_shear_unit(
    e1, 
    e2,
    gamma, 
    gamma1, 
    gamma2,
    num_th=500,
    maginf=-100.0,
):
    """
    Compute double-caustic cross sections for batched lens parameters.
    """
    theta = np.linspace(0, 2 * np.pi * (1.0 - 1.0 / num_th), num_th)

    phi, q = ellipticity2phi_q(e1, e2)
    
    # Closed-over constants (no recomputation)
    cos_th = np.cos(theta)
    sin_th = np.sin(theta)
    cos_2th = np.cos(2.0 * theta)
    sin_2th = np.sin(2.0 * theta)

    size = gamma.size
    cs_area = np.empty(size)

    # Original safe layout - still the fastest and simplest
    rcut_work  = np.empty((size, num_th))
    thcut_work = np.empty((size, num_th))
    rmain_work = np.empty((size, num_th))
    thmain_work = np.empty((size, num_th))
    r2_work    = np.empty((size, num_th))

    for i in prange(size):

        th_gamma = np.arctan2(gamma2[i], gamma1[i]) / 2.0
        g_mag = np.sqrt(gamma1[i]**2 + gamma2[i]**2)
        th_gamma -= phi[i]
        g1_rot = g_mag * np.cos(2.0 * th_gamma)
        g2_rot = g_mag * np.sin(2.0 * th_gamma)

        b_val = np.sqrt(q[i]) # * theta_E, theta_E=1.0
        t_val = gamma[i] - 1.0

        area = caustic_area_epl_shear(
            q[i], t_val, g1_rot, g2_rot, b_val,
            num_th, cos_th, sin_th, cos_2th, sin_2th,
            rcut_work[i], thcut_work[i], rmain_work[i], thmain_work[i], r2_work[i],
            maginf
        )

        cs_area[i] = area if np.isfinite(area) else 0.0

    return cs_area

@njit(cache=True, fastmath=True)
def phi_q2_ellipticity(phi, q):
    """
    Convert lens orientation and axis ratio to ellipticity components.

    Parameters
    ----------
    phi : ``float``
        Position angle of the lens major axis in radians.
    q : ``float``
        Axis ratio (minor/major), where ``0 < q <= 1``.

    Returns
    -------
    e1 : ``float``
        First ellipticity component.
    e2 : ``float``
        Second ellipticity component.

    Examples
    --------
    >>> e1, e2 = phi_q2_ellipticity(phi=0.25, q=0.8)
    """
    e = (1.0 - q) / (1.0 + q)
    return e * np.cos(2.0 * phi), e * np.sin(2.0 * phi)

@njit(cache=True, fastmath=True)
def ellipticity2phi_q(e1, e2):
    """
    Convert complex ellipticity moduli to orientation angle and axis ratio.

    Parameters
    ----------
    e1 : ``float``
        Eccentricity in x-direction.
    e2 : ``float``
        Eccentricity in xy-direction.

    Returns
    -------
    phi : ``float``
        Orientation angle in radians.
    q : ``float``
        Axis ratio (minor/major).
    """
    phi = np.arctan2(e2, e1) / 2
    c = np.sqrt(e1**2 + e2**2)
    c = np.minimum(c, 0.9999)
    q = (1 - c) / (1 + c)
    return phi, q

@njit(cache=True, fastmath=True)
def _omega_scalar(phi, t, q, niter_max=200, tol=1e-16):
    """
    Scalar version of ``omega`` that avoids temporary array allocation.
    """
    f = (1.0 - q) / (1.0 + q)
    omega_sum = 0.0 + 0.0j
    niter = min(niter_max, int(np.log(tol) / np.log(f)) + 2)
    Omega = np.exp(1j * phi)
    fact = -f * np.exp(2j * phi)

    for n in range(1, niter):
        omega_sum += Omega
        Omega *= (2.0 * n - (2.0 - t)) / (2.0 * n + (2.0 - t)) * fact
    omega_sum += Omega
    return omega_sum

@njit(cache=True)
def _alpha_epl_shear_scalar(x, y, b, q, t=1, gamma1=0, gamma2=0):
    """
    Scalar complex deflection for EPL + external shear.
    Use by epl_shear_njit.py.
    """
    
    zz = x * q + 1j * y
    R = np.abs(zz)
    phi = np.angle(zz)
    Omega = _omega_scalar(phi, t, q)
    alph = (2 * b) / (1 + q) * np.nan_to_num((b / R) ** t * R / b) * Omega

    return (
        alph
        + (gamma1 * x + gamma2 * y)
        + 1j * (gamma2 * x - gamma1 * y)
    )


@njit(cache=True)
def _interp_periodic(x, xp, fp, period):
    """
    Perform numba-compatible periodic linear interpolation.

    Equivalent to ``np.interp`` with the ``period`` keyword.

    Parameters
    ----------
    x : ``numpy.ndarray``
        Query points.
    xp : ``numpy.ndarray``
        Known x-values, sorted and within ``[0, period)``.
    fp : ``numpy.ndarray``
        Known function values corresponding to ``xp``.
    period : ``float``
        Period of the function.

    Returns
    -------
    result : ``numpy.ndarray``
        Interpolated values at query points.
    """
    n = xp.shape[0]
    result = np.empty_like(x)
    for i in range(x.shape[0]):
        xi = x[i] % period
        # binary search in sorted xp (assumed sorted and within [0, period))
        lo = 0
        hi = n - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if xp[mid] < xi:
                lo = mid + 1
            else:
                hi = mid
        # lo is the index of the first element >= xi
        if lo == 0:
            # xi is before the first point – wrap around
            x0 = xp[n - 1] - period
            x1 = xp[0]
            f0 = fp[n - 1]
            f1 = fp[0]
        else:
            x0 = xp[lo - 1]
            x1 = xp[lo]
            f0 = fp[lo - 1]
            f1 = fp[lo]
        dx = x1 - x0
        if abs(dx) < 1e-30:
            result[i] = f0
        else:
            result[i] = f0 + (f1 - f0) * (xi - x0) / dx
    return result

# @njit(cache=True, fastmath=True)
# def cart_to_pol(x, y):
#     """
#     Convert Cartesian coordinates to polar coordinates.

#     The returned angle is wrapped to ``[0, 2π)``.

#     Parameters
#     ----------
#     x : ``float`` or ``numpy.ndarray``
#         Cartesian x-coordinate.
#     y : ``float`` or ``numpy.ndarray``
#         Cartesian y-coordinate.

#     Returns
#     -------
#     r : ``float`` or ``numpy.ndarray``
#         Radial coordinate.
#     theta : ``float`` or ``numpy.ndarray``
#         Polar angle in radians, wrapped to ``[0, 2π)``.

#     Examples
#     --------
#     >>> r, theta = cart_to_pol(1.0, 1.0)
#     """
#     return np.sqrt(x**2 + y**2), np.arctan2(y, x) % (2 * np.pi)

# @njit(cache=True, fastmath=True)
# def _shear_cartesian2polar(gamma1, gamma2):
#     """
#     Convert Cartesian shear components to polar form.

#     Parameters
#     ----------
#     gamma1 : ``float``
#         First Cartesian shear component.
#     gamma2 : ``float``
#         Second Cartesian shear component.

#     Returns
#     -------
#     phi : ``float``
#         Shear position angle in radians.
#     gamma : ``float``
#         Shear magnitude.
#     """
#     phi = np.arctan2(gamma2, gamma1) / 2
#     gamma = np.sqrt(gamma1**2 + gamma2**2)
#     return phi, gamma

# @njit(cache=True, fastmath=True)
# def _shear_polar2cartesian(phi, gamma):
#     """
#     Convert polar shear representation to Cartesian components.

#     Parameters
#     ----------
#     phi : ``float``
#         Shear angle in radians.
#     gamma : ``float``
#         Shear magnitude.

#     Returns
#     -------
#     gamma1 : ``float``
#         First Cartesian shear component.
#     gamma2 : ``float``
#         Second Cartesian shear component.
#     """
#     gamma1 = gamma * np.cos(2 * phi)
#     gamma2 = gamma * np.sin(2 * phi)
#     return gamma1, gamma2

# @njit(cache=True, fastmath=True)
# def _rotmat(th):
#     """
#     Compute a 2D rotation matrix.

#     Parameters
#     ----------
#     th : ``float``
#         Rotation angle in radians.

#     Returns
#     -------
#     M : ``numpy.ndarray``
#         2×2 rotation matrix.
#     """
#     return np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]])

# @njit(cache=True, fastmath=True)
# def omega(phi, t, q, omegas, niter_max=200, tol=1e-16):
#     """
#     Evaluate the complex angular function Omega for the EPL profile.

#     This series expansion converges geometrically with ratio
#     ``f = (1 - q)/(1 + q)``. The ``fastmath`` flag provides ~4x speedup
#     due to the reduction nature of the summation.

#     Parameters
#     ----------
#     phi : ``numpy.ndarray``
#         Azimuthal angles in radians.
#     t : ``float``
#         EPL slope exponent (``t = gamma - 1``).
#     q : ``float``
#         Axis ratio.
#     omegas : ``numpy.ndarray``
#         Pre-allocated complex output buffer with the same shape as ``phi``.
#         Filled in-place.
#     niter_max : ``int``
#         Maximum number of series terms. \n
#         default: 200
#     tol : ``float``
#         Convergence tolerance. \n
#         default: 1e-16

#     Returns
#     -------
#     omegas : ``numpy.ndarray``
#         Complex Omega values at each angle (same object as input buffer).

#     Examples
#     --------
#     >>> import numpy as np
#     >>> phi = np.linspace(0, 2 * np.pi, 100)
#     >>> result = np.empty_like(phi, dtype=np.complex128)
#     >>> omega(phi, t=1.0, q=0.8, omegas=result)
#     """
#     f = (1 - q) / (1 + q)
#     # omegas[:] = 0.0 + 0.0j
#     niter = min(
#         niter_max, int(np.log(tol) / np.log(f)) + 2
#     )  # The absolute value of each summand is always less than f, hence this limit for the number of iterations.
#     Omega = 1 * np.exp(1j * phi)
#     fact = -f * np.exp(2j * phi)
#     for n in range(1, niter):
#         omegas += Omega
#         Omega *= (2 * n - (2 - t)) / (2 * n + (2 - t)) * fact
#     omegas += Omega
#     return omegas

# @njit(cache=True)
# def _solvequadeq(a, b, c):
#     """
#     Solve a quadratic equation with numerically careful root selection.

#     Uses sign-stabilized formulas to avoid loss of significance.
#     See https://en.wikipedia.org/wiki/Loss_of_significance.

#     Parameters
#     ----------
#     a : ``numpy.ndarray``
#         Quadratic coefficient.
#     b : ``numpy.ndarray``
#         Linear coefficient.
#     c : ``numpy.ndarray``
#         Constant coefficient.

#     Returns
#     -------
#     x1 : ``numpy.ndarray``
#         First root.
#     x2 : ``numpy.ndarray``
#         Second root.
#     """
#     sD = (b**2 - 4 * a * c) ** 0.5
#     x1 = (-b - np.sign(b) * sD) / (2 * a)
#     x2 = 2 * c / (-b - np.sign(b) * sD)
#     return np.where(b != 0, np.where(a != 0, x1, -c / b), -((-c / a) ** 0.5)), np.where(
#         b != 0, np.where(a != 0, x2, -c / b + 1e-8), +((-c / a) ** 0.5)
#     )

# @njit(cache=True, fastmath=True)
# def pol_to_cart(r, th):
#     """
#     Convert polar coordinates to Cartesian coordinates.

#     Parameters
#     ----------
#     r : ``float`` or ``numpy.ndarray``
#         Radial coordinate.
#     th : ``float`` or ``numpy.ndarray``
#         Polar angle in radians.

#     Returns
#     -------
#     x : ``float`` or ``numpy.ndarray``
#         Cartesian x-coordinate.
#     y : ``float`` or ``numpy.ndarray``
#         Cartesian y-coordinate.

#     Examples
#     --------
#     >>> x, y = pol_to_cart(1.0, np.pi / 4)
#     """
#     return r * np.cos(th), r * np.sin(th)

# @njit(cache=True)
# def _alpha_epl_shear(x, y, b, q, t, gamma1, gamma2, Omega):
#     """
#     Compute the complex deflection of an EPL + external shear lens.

#     Parameters
#     ----------
#     x : ``float`` or ``numpy.ndarray``
#         Image-plane x-coordinate.
#     y : ``float`` or ``numpy.ndarray``
#         Image-plane y-coordinate.
#     b : ``float``
#         EPL strength parameter (``theta_E * sqrt(q)``).
#     q : ``float``
#         Axis ratio.
#     t : ``float``
#         EPL slope exponent (``t = gamma - 1``).
#     gamma1 : ``float``
#         External shear component 1.
#     gamma2 : ``float``
#         External shear component 2.
#     Omega : ``numpy.ndarray``
#         Precomputed Omega array.

#     Returns
#     -------
#     deflection : ``complex`` or ``numpy.ndarray``
#         Complex deflection angle.
#     """
#     zz = x * q + 1j * y
#     R = np.abs(zz)
#     alph = (2 * b) / (1 + q) * np.nan_to_num((b / R) ** t * R / b) * Omega

#     return (
#         alph
#         + (gamma1 * x + gamma2 * y)
#         + 1j * (gamma2 * x - gamma1 * y)
#     )

## without root selection, for testing/debugging only
# @njit(cache=True, fastmath=True)
# def caustics_epl_shear(q, phi, t, gamma1, gamma2, b, theta, maginf=-100.0):
#     num_th = theta.size
    
#     # Pre-allocate 1D arrays ONLY for outputs we need to sort/interpolate
#     xca_cut = np.empty(num_th)
#     yca_cut = np.empty(num_th)
#     xca_4 = np.empty(num_th)
#     yca_4 = np.empty(num_th)
    
#     abs_gamma_sq = gamma1**2 + gamma2**2
    
#     # FUSED LOOP: Calculate everything element-by-element
#     for j in range(num_th):
#         th_val = theta[j]
#         cos_th = np.cos(th_val)
#         sin_th = np.sin(th_val)
        
#         R_val = np.sqrt(q**2 * cos_th**2 + sin_th**2)
#         phi_ell = np.arctan2(sin_th, cos_th * q)
        
#         omega_val = omega_scalar(phi_ell, t, q)
#         frac_roverR = 1.0 / R_val
        
#         # Assemble quadratic coefficients (cc = a, bb = b, aa = c for the solver)
#         one_minus_t = 1.0 - t
#         two_minus_t = 2.0 - t
#         one_plus_q = 1.0 + q
#         minus_one_by_t = -1.0 / t
#         two_pi = 2.0 * np.pi


#         cc = one_minus_t * two_minus_t * (cos_th * omega_val.real + sin_th * omega_val.imag) / frac_roverR * 2.0 / one_plus_q

#         cc -= one_minus_t**2 * (2.0 / one_plus_q)**2 * (omega_val.real**2 + omega_val.imag**2) / (frac_roverR**2)
        
#         exp_2j_th = np.cos(2*th_val) + 1j * np.sin(2*th_val)
#         exp_1j_th = cos_th + 1j * sin_th
        
#         gammaint_fac = -exp_2j_th * two_minus_t / 2.0 + one_minus_t * exp_1j_th * 2.0 / one_plus_q * omega_val / frac_roverR
        
#         bb = -two_minus_t - 2.0 * (gamma1 * gammaint_fac.real + gamma2 * gammaint_fac.imag)
#         aa = 1.0 - abs_gamma_sq
        
#         # Solve for main (quad) critical curve branch
#         x1_4, x2_4 = solve_quad_scalar(cc, bb, aa)
#         r_4 = b * (x2_4 ** minus_one_by_t) * frac_roverR
#         xcr_4 = r_4.real * cos_th
#         ycr_4 = r_4.real * sin_th
        
#         # Solve secondary (cut) branch
#         if t > 1.0:
#             x1_cut, x2_cut = solve_quad_scalar(cc, bb, aa - maginf)
#             u_cut = x2_cut
#         else:
#             x1_cut, x2_cut = solve_quad_scalar(cc, bb, aa + maginf)
#             u_cut = x1_cut
            
#         r_cut = b * (u_cut ** minus_one_by_t) * frac_roverR
#         xcr_cut = r_cut.real * cos_th
#         ycr_cut = r_cut.real * sin_th
        
#         # Compute deflections and map to source plane
#         al_4 = alpha_epl_shear_scalar(xcr_4, ycr_4, b, q, t, gamma1, gamma2, omega_val)
#         al_cut = alpha_epl_shear_scalar(xcr_cut, ycr_cut, b, q, t, gamma1, gamma2, omega_val)
        
#         xca_4[j] = xcr_4 - al_4.real
#         yca_4[j] = ycr_4 - al_4.imag
#         xca_cut[j] = xcr_cut - al_cut.real
#         yca_cut[j] = ycr_cut - al_cut.imag

#     # Extract radii and map angles strictly to [0, 2*pi]
#     rcut = np.sqrt(xca_cut**2 + yca_cut**2)
#     thcut = np.arctan2(yca_cut, xca_cut) % two_pi
    
#     r_main = np.sqrt(xca_4**2 + yca_4**2)
#     th_main = np.arctan2(yca_4, xca_4) % two_pi
    
#     # Sort cut branch for interpolation
#     sort_idx = np.argsort(thcut)
#     thcut_sorted = thcut[sort_idx]
#     rcut_sorted = rcut[sort_idx]
    
#     # --- Manual Periodic Interpolation for Numba ---
#     # Numba's np.interp doesn't support kwargs, so we manually wrap the endpoints
#     thcut_pad = np.empty(num_th + 2)
#     rcut_pad = np.empty(num_th + 2)
#     thcut_pad[1:-1] = thcut_sorted
#     rcut_pad[1:-1] = rcut_sorted
#     thcut_pad[0] = thcut_sorted[-1] - two_pi
#     rcut_pad[0] = rcut_sorted[-1]
#     thcut_pad[-1] = thcut_sorted[0] + two_pi
#     rcut_pad[-1] = rcut_sorted[0]
    
#     r2 = np.interp(th_main, thcut_pad, rcut_pad)
    
#     # Compose outer double-image boundary
#     r_final = np.fmax(r_main, r2)
    
#     pos_x = r_final * np.cos(th_main)
#     pos_y = r_final * np.sin(th_main)
    
#     # Rotate final coordinates back to original frame
#     mphi = -phi
#     cos_mphi = np.cos(mphi)
#     sin_mphi = np.sin(mphi)
    
#     out = np.empty((2, num_th))
#     out[0, :] = cos_mphi * pos_x + sin_mphi * pos_y
#     out[1, :] = -sin_mphi * pos_x + cos_mphi * pos_y
    
#     return out