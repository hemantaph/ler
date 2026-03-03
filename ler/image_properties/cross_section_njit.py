import numpy as np
from numba import njit, prange

C_LIGHT = 299792458.0  # m/s


@njit(cache=True)
def phi_q2_ellipticity(phi, q):
    """
    Convert orientation angle and axis ratio to ellipticity components.
    Parameters
    ----------
    phi : ``float``
        Position angle of the lens major axis (rad).
    q : ``float``
        Axis ratio (minor/major), 0 < q <= 1.
    Returns
    -------
    e1 : ``float``
        First ellipticity component.
    e2 : ``float``
        Second ellipticity component.
    """
    e = (1.0 - q) / (1.0 + q)
    return e * np.cos(2.0 * phi), e * np.sin(2.0 * phi)


@njit(cache=True)
def _cdot(a, b):
    """
    Real dot-product for complex numbers: Re(a * conj(b)).
    Parameters
    ----------
    a : ``complex128``
        First complex number.
    b : ``complex128``
        Second complex number.
    Returns
    -------
    result : ``float``
        a.real * b.real + a.imag * b.imag.
    """
    return a.real * b.real + a.imag * b.imag


@njit(cache=True)
def _pol_to_ell_R_phi(theta, q):
    """
    Polar to elliptical coordinate transform (unit radius).
    Reproduces lenstronomy's ``pol_to_ell(r=1, theta, q)``.
    Parameters
    ----------
    theta : ``float``
        Polar angle (rad).
    q : ``float``
        Axis ratio.
    Returns
    -------
    R : ``float``
        Elliptical radial coordinate.
    phi : ``float``
        Elliptical angular coordinate.
    """
    c = np.cos(theta)
    s = np.sin(theta)
    R = np.sqrt((q * q) * (c * c) + (s * s))
    phi = np.arctan2(s, c * q)
    return R, phi


@njit(cache=True)
def _solve_quad_all(cc, bb, aa):
    """
    Solve the quadratic equation cc * u^2 + bb * u + aa = 0.
    Uses a numerically stable formulation to reduce cancellation.
    Parameters
    ----------
    cc : ``numpy.ndarray``
        Quadratic coefficient.
    bb : ``numpy.ndarray``
        Linear coefficient.
    aa : ``numpy.ndarray``
        Constant coefficient.
    Returns
    -------
    u0 : ``numpy.ndarray``
        First root (NaN if discriminant < 0).
    u1 : ``numpy.ndarray``
        Second root (NaN if discriminant < 0).
    """
    d = bb * bb - 4.0 * cc * aa
    u0 = np.full_like(d, np.nan)
    u1 = np.full_like(d, np.nan)
    mask = d >= 0.0
    if not np.any(mask):
        return u0, u1
    s = np.zeros_like(d)
    s[mask] = np.sqrt(d[mask])
    dn = np.zeros_like(d)
    mask_pos = mask & (bb >= 0.0)
    dn[mask_pos] = -bb[mask_pos] - s[mask_pos]
    mask_neg = mask & (bb < 0.0)
    dn[mask_neg] = -bb[mask_neg] + s[mask_neg]
    small_dn = np.abs(dn) < 1e-30
    mask_small = mask & small_dn
    u0[mask_small] = (-bb[mask_small] + s[mask_small]) / (2.0 * cc[mask_small])
    u1[mask_small] = (-bb[mask_small] - s[mask_small]) / (2.0 * cc[mask_small])
    mask_not_small = mask & ~small_dn
    u0[mask_not_small] = dn[mask_not_small] / (2.0 * cc[mask_not_small])
    u1[mask_not_small] = (2.0 * aa[mask_not_small]) / dn[mask_not_small]
    return u0, u1


@njit(cache=True)
def _interp_periodic(theta, xp, fp):
    """
    Periodic linear interpolation on [0, 2 pi).
    Parameters
    ----------
    theta : ``numpy.ndarray``
        Query angles (rad).
    xp : ``numpy.ndarray``
        Sorted sample angles in [0, 2 pi).
    fp : ``numpy.ndarray``
        Function values at ``xp``.
    Returns
    -------
    value : ``numpy.ndarray``
        Interpolated function values at ``theta``.
    """
    twopi = 2.0 * np.pi
    t = np.mod(theta, twopi)
    N = xp.size
    idx = np.searchsorted(xp, t)
    val = np.zeros_like(t)
    mask_low = idx == 0
    mask_high = idx == N
    mask_mid = ~(mask_low | mask_high)
    # mid
    m = mask_mid
    if np.any(m):
        lo = idx[m] - 1
        hi = idx[m]
        x0 = xp[lo]
        x1 = xp[hi]
        y0 = fp[lo]
        y1 = fp[hi]
        w = (t[m] - x0) / (x1 - x0)
        val[m] = y0 + w * (y1 - y0)
    # low
    mask_low_strict = mask_low & (t < xp[0])
    mask_low_eq = mask_low & (t >= xp[0])
    # low_eq
    m = mask_low_eq
    if np.any(m):
        lo = 0
        hi = 1
        x0 = xp[lo]
        x1 = xp[hi]
        y0 = fp[lo]
        y1 = fp[hi]
        w = (t[m] - x0) / (x1 - x0)
        val[m] = y0 + w * (y1 - y0)
    # low_strict
    m = mask_low_strict
    if np.any(m):
        x0 = xp[-1]
        x1 = xp[0] + twopi
        y0 = fp[-1]
        y1 = fp[0]
        w = (t[m] + twopi - x0) / (x1 - x0)
        val[m] = y0 + w * (y1 - y0)
    # high
    m = mask_high
    if np.any(m):
        x0 = xp[-1]
        x1 = xp[0] + twopi
        y0 = fp[-1]
        y1 = fp[0]
        w = (t[m] - x0) / (x1 - x0)
        val[m] = y0 + w * (y1 - y0)
    return val


@njit(cache=True)
def _shear_c2p(g1, g2):
    """
    Convert Cartesian shear components to polar form.
    Parameters
    ----------
    g1 : ``float``
        First shear component.
    g2 : ``float``
        Second shear component.
    Returns
    -------
    phi_gamma : ``float``
        Shear position angle (rad).
    gamma : ``float``
        Shear magnitude.
    """
    phi = 0.5 * np.arctan2(g2, g1)
    g = np.sqrt(g1 * g1 + g2 * g2)
    return phi, g


@njit(cache=True)
def _shear_p2c(phi, g):
    """
    Convert polar shear to Cartesian components.
    Parameters
    ----------
    phi : ``float``
        Shear position angle (rad).
    g : ``float``
        Shear magnitude.
    Returns
    -------
    g1 : ``float``
        First shear component.
    g2 : ``float``
        Second shear component.
    """
    return g * np.cos(2.0 * phi), g * np.sin(2.0 * phi)


@njit(cache=True)
def _omega(phi, t, q):
    """
    Omega series for EPL deflections.
    The series length is adaptively chosen as
    ``ni ~ log(1e-16) / log(|f|)``, clipped to [2, 200].
    Parameters
    ----------
    phi : ``numpy.ndarray``
        Elliptical angles (rad).
    t : ``float``
        Power-law slope parameter, ``t = gamma - 1``.
    q : ``float``
        Axis ratio.
    Returns
    -------
    Om : ``numpy.ndarray``
        Complex Omega values.
    """
    f = (1.0 - q) / (1.0 + q)
    absf = abs(f)
    if absf < 1e-15:
        return np.exp(1j * phi)
    ni = int(np.floor(np.log(1e-16) / np.log(absf + 1e-30)) + 2.0)
    ni = max(2, min(ni, 200))
    Om = np.exp(1j * phi)
    fc = (-f) * np.exp(1j * 2.0 * phi)
    res = np.zeros_like(Om)
    for n in range(1, ni):
        res += Om
        num = 2.0 * n - (2.0 - t)
        den = 2.0 * n + (2.0 - t)
        Om = Om * (num / den) * fc
    return res + Om


@njit(cache=True)
def _alpha_epl_xy(x, y, theta_E, q, gamma):
    """
    EPL deflection angle in the solver frame (major axis along x).
    Parameters
    ----------
    x : ``numpy.ndarray``
        Image x-coordinates (solver frame).
    y : ``numpy.ndarray``
        Image y-coordinates (solver frame).
    theta_E : ``float``
        Einstein radius.
    q : ``float``
        Axis ratio.
    gamma : ``float``
        Power-law slope.
    Returns
    -------
    ax : ``numpy.ndarray``
        Deflection in x.
    ay : ``numpy.ndarray``
        Deflection in y.
    """
    t = gamma - 1.0
    b = theta_E * np.sqrt(q)
    R = np.sqrt((q * x) * (q * x) + y * y)
    pf = (2.0 * b) / (1.0 + q)
    ax = np.zeros_like(x)
    ay = np.zeros_like(y)
    mask = R > 1e-15
    if not np.any(mask):
        return ax, ay
    xm = x[mask]
    ym = y[mask]
    Rm = R[mask]
    fr = (b / Rm) ** t * (Rm / b)
    ph = np.arctan2(ym, xm * q)
    Om = _omega(ph, t, q)
    a = pf * fr * Om
    ax[mask] = np.real(a)
    ay[mask] = np.imag(a)
    return ax, ay


@njit(cache=True)
def _alpha_shear_xy(x, y, g1, g2):
    """
    External shear deflection.
    Parameters
    ----------
    x : ``numpy.ndarray``
        Image x-coordinates.
    y : ``numpy.ndarray``
        Image y-coordinates.
    g1 : ``float``
        Shear component 1 (solver frame).
    g2 : ``float``
        Shear component 2 (solver frame).
    Returns
    -------
    ax : ``numpy.ndarray``
        Shear deflection in x.
    ay : ``numpy.ndarray``
        Shear deflection in y.
    """
    return g1 * x + g2 * y, g2 * x - g1 * y


@njit(cache=True)
def _alpha_total_xy(x, y, theta_E, q, gamma, g1, g2):
    """
    Total deflection (EPL + shear) in the solver frame.
    Parameters
    ----------
    x : ``numpy.ndarray``
        Image x-coordinates.
    y : ``numpy.ndarray``
        Image y-coordinates.
    theta_E : ``float``
        Einstein radius.
    q : ``float``
        Axis ratio.
    gamma : ``float``
        Power-law slope.
    g1 : ``float``
        Shear component 1 (solver frame).
    g2 : ``float``
        Shear component 2 (solver frame).
    Returns
    -------
    ax : ``numpy.ndarray``
        Total deflection in x.
    ay : ``numpy.ndarray``
        Total deflection in y.
    """
    ax, ay = _alpha_epl_xy(x, y, theta_E, q, gamma)
    sx, sy = _alpha_shear_xy(x, y, g1, g2)
    return ax + sx, ay + sy


@njit(cache=True)
def _rot_apply(theta, x, y):
    """
    Apply 2-D rotation matrix to points.
    The rotation convention is ``[[c, s], [-s, c]]``.
    Parameters
    ----------
    theta : ``float``
        Rotation angle (rad).
    x : ``numpy.ndarray``
        x-coordinates.
    y : ``numpy.ndarray``
        y-coordinates.
    Returns
    -------
    x_rot : ``numpy.ndarray``
        Rotated x-coordinates.
    y_rot : ``numpy.ndarray``
        Rotated y-coordinates.
    """
    c = np.cos(theta)
    s = np.sin(theta)
    return c * x + s * y, -s * x + c * y


# ---------------------------------------------------------------------------
# Public — caustic computation
# ---------------------------------------------------------------------------
@njit(cache=True)
def caustic_double_points_epl_shear(
    q,
    phi,
    gamma,
    gamma1_unrot,
    gamma2_unrot,
    theta_E=1.0,
    num_th=500,
    maginf=-100.0,
    sourceplane=True,
):
    """
    Compute the double-imaging boundary (double caustic) for an EPL + shear lens.
    """
    t = gamma - 1.0
    b = theta_E * np.sqrt(q)
    # rotate shear into lens major-axis frame
    theta_gamma, gamma_mag = _shear_c2p(gamma1_unrot, gamma2_unrot)
    theta_gamma -= phi
    g1, g2 = _shear_p2c(theta_gamma, gamma_mag)
    # sample angles
    theta = np.arange(num_th) * (2.0 * np.pi / num_th)
    # vectorized computations
    c = np.cos(theta)
    s = np.sin(theta)
    R = np.sqrt(q**2 * c**2 + s**2)
    phiell = np.arctan2(s, c * q)
    frac_roverR = 1.0 / R
    Om = _omega(phiell, t, q)
    e1t = c + 1j * s
    e2t = np.cos(2.0 * theta) + 1j * np.sin(2.0 * theta)
    # quadratic coefficients
    aa0 = 1.0 - (g1 * g1 + g2 * g2)
    bb0 = -(2.0 - t)
    const_fac = 2.0 / (1.0 + q)
    # cc coefficient
    cc = (1.0 - t) * (2.0 - t) * np.real(e1t * np.conj(Om)) / frac_roverR * const_fac
    cc -= (
        (1.0 - t)
        * (1.0 - t)
        * (const_fac * const_fac)
        * (np.abs(Om) ** 2)
        / (frac_roverR * frac_roverR)
    )
    # shear interaction
    gammaint_fac = (-e2t * (2.0 - t) / 2.0) + (
        1.0 - t
    ) * e1t * const_fac * Om / frac_roverR
    bb = bb0 - 2.0 * np.real((g1 + 1j * g2) * np.conj(gammaint_fac))
    aa = np.full_like(bb, aa0)
    # 4-image critical curve (second root)
    u0, u1 = _solve_quad_all(cc, bb, aa)
    u = u1
    r4 = np.full_like(u, np.nan)
    mask_r4 = np.isfinite(u) & (u > 0.0)
    r4[mask_r4] = b * (u[mask_r4] ** (-1.0 / t)) * frac_roverR[mask_r4]
    # outer (cut) curve
    if t > 1.0:
        u0c, u1c = _solve_quad_all(cc, bb, aa - maginf)
        uc = u1c
    else:
        u0c, u1c = _solve_quad_all(cc, bb, aa + maginf)
        uc = u0c
    rc = np.full_like(uc, np.nan)
    mask_rc = np.isfinite(uc) & (uc > 0.0)
    rc[mask_rc] = b * (uc[mask_rc] ** (-1.0 / t)) * frac_roverR[mask_rc]
    # image-plane points
    xcr4 = r4 * c
    ycr4 = r4 * s
    xcrc = rc * c
    ycrc = rc * s
    # map to source plane if requested
    if sourceplane:
        a4x, a4y = _alpha_total_xy(xcr4, ycr4, theta_E, q, gamma, g1, g2)
        acx, acy = _alpha_total_xy(xcrc, ycrc, theta_E, q, gamma, g1, g2)
        xca4 = xcr4 - a4x
        yca4 = ycr4 - a4y
        xcac = xcrc - acx
        ycac = ycrc - acy
    else:
        xca4 = xcr4
        yca4 = ycr4
        xcac = xcrc
        ycac = ycrc
    # convert both curves to polar (r, theta)
    th4 = np.mod(np.arctan2(yca4, xca4), 2.0 * np.pi)
    rc4 = np.sqrt(xca4 * xca4 + yca4 * yca4)
    thc = np.mod(np.arctan2(ycac, xcac), 2.0 * np.pi)
    rcc = np.sqrt(xcac * xcac + ycac * ycac)
    # sort by angle for interpolation
    idxc = np.argsort(thc)
    thc_s = thc[idxc]
    rcc_s = rcc[idxc]
    idx4 = np.argsort(th4)
    th4_s = th4[idx4]
    rc4_s = rc4[idx4]
    # build double boundary: r = max(r4, r_cut(th))
    r_out = np.maximum(rc4_s, _interp_periodic(th4_s, thc_s, rcc_s))
    # back to Cartesian
    x_out = r_out * np.cos(th4_s)
    y_out = r_out * np.sin(th4_s)
    # rotate back to sky frame
    xs, ys = _rot_apply(-phi, x_out, y_out)
    pts = np.empty((2, num_th), dtype=np.float64)
    pts[0] = xs
    pts[1] = ys
    return pts


# ---------------------------------------------------------------------------
# Public API — caustic area
# ---------------------------------------------------------------------------
@njit(cache=True)
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
    """
    area = 0.0
    n = xv.size
    j = n - 1
    for i in range(n):
        area += xv[j] * yv[i] - xv[i] * yv[j]
        j = i
    return abs(area) / 2.0


@njit(cache=True)
def caustic_double_area(
    q, phi, gamma, gamma1, gamma2, theta_E=1.0, num_th=500, maginf=-100.0
):
    """
    Compute the exact area inside the double caustic of an EPL+Shear lens.
    Parameters
    ----------
    (Matches caustic_double_points_epl_shear signatures)
    Returns
    -------
    area : ``float``
        Area of the double caustic in units of theta_E^2.
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
    return polygon_area(pts[0], pts[1])


# ---------------------------------------------------------------------------
# Public API — cross section
# ---------------------------------------------------------------------------
def make_cross_section_reinit(
    Da_instance,
    num_th=500,
    maginf=-100.0,
):
    """
    Factory function to create a JIT-compiled cross section calculator.
    """

    @njit(parallel=True, cache=True)
    def cross_section_reinit(zs, zl, sigma, q, phi, gamma, gamma1, gamma2):
        """
        Compute the cross section for a given set of lens parameters.
        """
        size = zs.size
        Ds = Da_instance(zs)
        Dl = Da_instance(zl)
        Dls = (Ds * (1 + zs) - Dl * (1 + zl)) / (1 + zs)
        theta_E = 4.0 * np.pi * (sigma * 1000.0 / C_LIGHT) ** 2 * (Dls / Ds)

        cs = np.zeros(size)
        for i in prange(size):

            cs[i] = caustic_double_area(
                q[i],
                phi[i],
                gamma[i],
                gamma1[i],
                gamma2[i],
                theta_E=theta_E[i],
                num_th=num_th,
                maginf=maginf,
            )
        return cs

    return cross_section_reinit
