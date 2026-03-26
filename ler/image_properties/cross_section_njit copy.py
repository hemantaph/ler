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
def pol_to_ell(r, theta, q):
    """
    Convert polar coordinates to elliptical coordinates.
    Use in epl_shear_njit.py

    Parameters
    ----------
    r : ``float`` or ``numpy.ndarray``
        Radial coordinate.
    theta : ``float`` or ``numpy.ndarray``
        Polar angle in radians.
    q : ``float``
        Axis ratio.

    Returns
    -------
    rell : ``float`` or ``numpy.ndarray``
        Elliptical radial coordinate.
    phi : ``float`` or ``numpy.ndarray``
        Elliptical angle in radians.
    """
    phi = np.arctan2(np.sin(theta), np.cos(theta) * q)
    rell = r * np.sqrt(q**2 * np.cos(theta) ** 2 + np.sin(theta) ** 2)
    return rell, phi

@njit(cache=True, fastmath=True, inline="always")
def omega_scalar(phi, t, q, niter_max=200, tol=1e-16):
    """
    Return Re(Omega), Im(Omega) for the EPL angular function.

    Uses a purely real recurrence to avoid complex arithmetic in the hot path.
    """
    f = (1.0 - q) / (1.0 + q)

    # Circular limit
    if f <= 0.0:
        return np.cos(phi), np.sin(phi)

    niter = min(niter_max, int(np.log(tol) / np.log(f)) + 2)
    if niter < 1:
        niter = 1

    c1 = np.cos(phi)
    s1 = np.sin(phi)
    c2 = np.cos(2.0 * phi)
    s2 = np.sin(2.0 * phi)

    # Omega_0 = exp(i phi)
    ur = c1
    ui = s1

    # fact = -f * exp(2 i phi)
    fr = -f * c2
    fi = -f * s2

    # running sum
    res_r = 0.0
    res_i = 0.0

    for n in range(1, niter):
        res_r += ur
        res_i += ui

        an = (2.0 * n - (2.0 - t)) / (2.0 * n + (2.0 - t))

        tmp_r = ur * fr - ui * fi
        tmp_i = ur * fi + ui * fr

        ur = an * tmp_r
        ui = an * tmp_i

    res_r += ur
    res_i += ui
    return res_r, res_i

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

@njit(cache=True, fastmath=True, inline="always")
def solve_quad_scalar(a, b, c):
    if b != 0.0:
        if a != 0.0:
            disc = b * b - 4.0 * a * c + 0j
            sD = np.sqrt(disc)
            sign_b = 1.0 if b > 0.0 else -1.0
            denominator = -b - sign_b * sD

            # instead of 'if denominator != 0.0:' for numerical safety
            if np.abs(denominator) > 1e-30:
                x1 = denominator / (2.0 * a)
                x2 = (2.0 * c) / denominator
            else:
                x1 = (-b + sD) / (2.0 * a)
                x2 = x1
        else:
            x1 = -c / b + 0j
            x2 = x1
    else:
        if a != 0.0:
            val = np.sqrt(-c / a + 0j)
            x1 = -val
            x2 = val
        else:
            x1 = 0.0 + 0.0j
            x2 = 0.0 + 0.0j
    return x1, x2

@njit(cache=True, fastmath=True, inline="always")
def alpha_epl_shear_scalar(x, y, b, q, t, gamma1, gamma2, omega_real, omega_imag):
    """Calculates total deflection components (real, imag) for a single point."""
    xq = x * q
    R = np.sqrt(xq * xq + y * y)
    prefac = (2* b**t)/(1 + q) * R**(1-t)
    alph_real = prefac * omega_real
    alph_imag = prefac * omega_imag

    return (
        alph_real + (gamma1 * x + gamma2 * y),
        alph_imag + (gamma2 * x - gamma1 * y),
    )

@njit(cache=True, fastmath=True, inline="always")
def select_physical_root(u0, u1, prefer_second):
    """
    Select a physically meaningful inverse-radius root.

    Priority:
    1) Root with small imaginary part and positive real part.
    2) If both valid, keep the analytically expected branch.
    3) If neither valid, fall back to preferred branch real part if positive,
       otherwise to the alternate branch if positive, else a small floor.
    """
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

@njit(cache=True, fastmath=True)
def caustics_epl_shear(q, phi, t, gamma1, gamma2, b, theta, maginf=-100.0):
    # b = np.sqrt(q[i]) * theta_E[i]
    num_th = theta.size
    
    # Pre-allocate 1D arrays ONLY for outputs we need to sort/interpolate
    xca_cut = np.empty(num_th)
    yca_cut = np.empty(num_th)
    xca_4 = np.empty(num_th)
    yca_4 = np.empty(num_th)
    
    abs_gamma_sq = gamma1**2 + gamma2**2
    one_minus_t = 1.0 - t
    one_minus_t_sq = one_minus_t**2
    two_minus_t = 2.0 - t
    two_by_one_plus_q = 2.0 / (1.0 + q)
    minus_one_by_t = -1.0 / t
    aa = 1.0 - abs_gamma_sq
    
    # FUSED LOOP: Calculate everything element-by-element
    for j in range(num_th):
        th_val = theta[j]
        cos_th = np.cos(th_val)
        sin_th = np.sin(th_val)
        
        R_val = np.sqrt(q**2 * cos_th**2 + sin_th**2)
        phi_ell = np.arctan2(sin_th, cos_th * q)
        
        omega_real, omega_imag = omega_scalar(phi_ell, t, q)
        frac_roverR = 1.0 / R_val
        
        # Assemble quadratic coefficients (cc = a, bb = b, aa = c for the solver)
        cc = one_minus_t * two_minus_t * (cos_th * omega_real + sin_th * omega_imag) * R_val * two_by_one_plus_q
        cc -= one_minus_t_sq * (two_by_one_plus_q**2) * (omega_real**2 + omega_imag**2) * (R_val * R_val)

        cos_2th = np.cos(2.0 * th_val)
        sin_2th = np.sin(2.0 * th_val)
        shear_prefac = one_minus_t * two_by_one_plus_q * R_val

        gammaint_fac_real = -cos_2th * two_minus_t / 2.0
        gammaint_fac_real += shear_prefac * (cos_th * omega_real - sin_th * omega_imag)

        gammaint_fac_imag = -sin_2th * two_minus_t / 2.0
        gammaint_fac_imag += shear_prefac * (cos_th * omega_imag + sin_th * omega_real)

        bb = -two_minus_t - 2.0 * (gamma1 * gammaint_fac_real + gamma2 * gammaint_fac_imag)
        
        # Solve for main (quad) critical curve branch
        x1_4, x2_4 = solve_quad_scalar(cc, bb, aa)
        u4 = select_physical_root(x1_4, x2_4, True)
        r_4 = b * (u4 ** minus_one_by_t) * frac_roverR
        xcr_4 = r_4 * cos_th
        ycr_4 = r_4 * sin_th
        
        # Solve secondary (cut) branch
        if t > 1.0:
            x1_cut, x2_cut = solve_quad_scalar(cc, bb, aa - maginf)
            u_cut = select_physical_root(x1_cut, x2_cut, True)
        else:
            x1_cut, x2_cut = solve_quad_scalar(cc, bb, aa + maginf)
            u_cut = select_physical_root(x1_cut, x2_cut, False)
            
        r_cut = b * (u_cut ** minus_one_by_t) * frac_roverR
        xcr_cut = r_cut * cos_th
        ycr_cut = r_cut * sin_th
        
        # Compute deflections and map to source plane
        al_4_real, al_4_imag = alpha_epl_shear_scalar(xcr_4, ycr_4, b, q, t, gamma1, gamma2, omega_real, omega_imag)
        al_cut_real, al_cut_imag = alpha_epl_shear_scalar(xcr_cut, ycr_cut, b, q, t, gamma1, gamma2, omega_real, omega_imag)
        
        # xca_4_ = xcr_4 - al_4_real
        # yca_4_ = ycr_4 - al_4_imag
        # xca_cut_ = xcr_cut - al_cut_real
        # yca_cut_ = ycr_cut - al_cut_imag
        # xca_4[j] = xca_4_
        # yca_4[j] = yca_4_
        # xca_cut[j] = xca_cut_
        # yca_cut[j] = yca_cut_
        xca_4[j] = xcr_4 - al_4_real
        yca_4[j] = ycr_4 - al_4_imag
        xca_cut[j] = xcr_cut - al_cut_real
        yca_cut[j] = ycr_cut - al_cut_imag

    # Extract radii and map angles strictly to [0, 2*pi]
    rcut = np.sqrt(xca_cut**2 + yca_cut**2)
    thcut = np.arctan2(yca_cut, xca_cut) % TWO_PI
    
    r_main = np.sqrt(xca_4**2 + yca_4**2)
    th_main = np.arctan2(yca_4, xca_4) % TWO_PI
    
    # Sort cut branch for interpolation
    sort_idx = np.argsort(thcut)
    thcut_sorted = thcut[sort_idx]
    rcut_sorted = rcut[sort_idx]
    
    # # --- Manual Periodic Interpolation for Numba ---
    # # Numba's np.interp doesn't support kwargs, so we manually wrap the endpoints
    # thcut_pad = np.empty(num_th + 2)
    # rcut_pad = np.empty(num_th + 2)
    # thcut_pad[1:-1] = thcut_sorted
    # rcut_pad[1:-1] = rcut_sorted
    # thcut_pad[0] = thcut_sorted[-1] - TWO_PI
    # rcut_pad[0] = rcut_sorted[-1]
    # thcut_pad[-1] = thcut_sorted[0] + TWO_PI
    # rcut_pad[-1] = rcut_sorted[0]
    # r2 = np.interp(th_main, thcut_pad, rcut_pad)

    # use _interp_periodic function
    r2 = _interp_periodic(th_main, thcut_sorted, rcut_sorted, TWO_PI)
    
    # Compose outer double-image boundary
    r_final = np.fmax(r_main, r2)
    
    pos_x = r_final * np.cos(th_main)
    pos_y = r_final * np.sin(th_main)
    
    # Rotate final coordinates back to original frame
    mphi = -phi
    cos_mphi = np.cos(mphi)
    sin_mphi = np.sin(mphi)
    
    out = np.empty((2, num_th))
    out[0, :] = cos_mphi * pos_x + sin_mphi * pos_y
    out[1, :] = -sin_mphi * pos_x + cos_mphi * pos_y
    
    return out

@njit(cache=True, fastmath=True)
def polygon_area(xv, yv):
    """
    Compute the area of a simple polygon using the Shoelace formula.

    Parameters
    ----------
    xv : ``numpy.ndarray``
        Polygon vertex x-coordinates.
    yv : ``numpy.ndarray``
        Polygon vertex y-coordinates.

    Returns
    -------
    area : ``float``
        The enclosed geometric area.

    Examples
    --------
    >>> area = polygon_area(np.array([0.0, 1.0, 0.0]), np.array([0.0, 0.0, 1.0]))
    """
    area = 0.0
    n = xv.size
    j = n - 1
    for i in range(n):
        area += xv[j] * yv[i] - xv[i] * yv[j]
        j = i
    return abs(area) / 2.0

def make_cross_section_area_reinit(
    Da_instance,
    num_th=500, 
    maginf=-100.0,
):
    """
    Create a JIT-compiled cross-section evaluator for batched systems.
    """
    theta = np.linspace(0, 2 * np.pi * (1.0 - 1.0 / num_th), num_th)
    
    @njit(parallel=True, cache=True, fastmath=True)
    def cross_section_area(zs, zl, sigma, q, phi, gamma, gamma1, gamma2):
        """
        Compute double-caustic cross sections for batched lens parameters.
        """
        size = zs.size

        # Output array
        cs_area = np.empty(size) 

        # Distance and Einstein radius scalar math
        Ds = Da_instance(zs)
        Dl = Da_instance(zl)
        Dls = (Ds * (1.0 + zs) - Dl * (1.0 + zl)) / (1.0 + zs)
        theta_E = 4.0 * PI * (sigma * 1000.0 / C_LIGHT) ** 2 * (Dls / Ds)
        
        # Local polar shear conversion
        th_gamma = np.arctan2(gamma2, gamma1) / 2.0
        g_mag = np.sqrt(gamma1**2 + gamma2**2)
        
        # Rotate shear
        th_gamma -= phi
        g1_rot = g_mag * np.cos(2.0 * th_gamma)
        g2_rot = g_mag * np.sin(2.0 * th_gamma)
        
        b_val = np.sqrt(q) * theta_E
        t_val = gamma - 1.0
        
        # Process everything as scalars inside the parallel loop
        for i in prange(size):
            
            pts = caustics_epl_shear(
                q=q[i], 
                phi=phi[i], 
                t=t_val[i], 
                gamma1=g1_rot[i], 
                gamma2=g2_rot[i], 
                b=b_val[i], 
                theta=theta, 
                maginf=maginf
            )
            
            area = polygon_area(pts[0, :], pts[1, :])
            if not np.isfinite(area):
                area = 0.0
            cs_area[i] = area
            
        return cs_area
        
    return cross_section_area

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
    size = gamma.size
    # e1 and e2 are converted to phi and q
    phi, q = ellipticity2phi_q(e1, e2)

    # Output array only (no pts storage)
    cs_area = np.empty(size) 
    
    # Local polar shear conversion
    th_gamma = np.arctan2(gamma2, gamma1) / 2.0
    g_mag = np.sqrt(gamma1**2 + gamma2**2)
    
    # Rotate shear
    th_gamma -= phi
    g1_rot = g_mag * np.cos(2.0 * th_gamma)
    g2_rot = g_mag * np.sin(2.0 * th_gamma)
    
    b_val = np.sqrt(q) * 1.0
    t_val = gamma - 1.0

    theta = np.linspace(0, 2 * np.pi * (1.0 - 1.0 / num_th), num_th)

    for i in prange(size):

        pts = caustics_epl_shear(
            q=q[i], 
            phi=phi[i], 
            t=t_val[i], 
            gamma1=g1_rot[i], 
            gamma2=g2_rot[i], 
            b=b_val[i], 
            theta=theta, 
            maginf=maginf
        )
        
        area = polygon_area(pts[0, :], pts[1, :])
        if not np.isfinite(area):
            area = 0.0
        cs_area[i] = area

    return np.maximum(cs_area, 0.0)

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