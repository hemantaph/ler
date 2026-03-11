"""
Module for semi-analytical image finding in EPL + external shear lens models.

Provides a numba-accelerated solver that locates multiple lensed images,
computes their magnifications, arrival times, and image types for an
elliptical power-law (EPL) mass profile with external shear.

Usage:
    Solve for image positions of a single source:

    >>> from ler.image_properties.epl_shear_njit import image_position_analytical_njit
    >>> x_img, y_img, tau, mu, itype, n = image_position_analytical_njit(
    ...     x=0.1, y=0.05, q=0.8, phi=0.3, gamma=2.0,
    ...     gamma1=0.03, gamma2=-0.01
    ... )

Copyright (C) 2026 Phurailatpam Hemantakumar. Distributed under MIT License.
"""

import numpy as np
from numba import njit, prange
from .sample_caustic_points_njit import sample_source_from_double_caustic
from .cross_section_njit import omega_scalar, cdot, pol_to_ell, _alpha_epl_shear_scalar


EPS = 1e-14
MAX_ROOTS = 64
MAX_IMGS = 5
C_LIGHT = 299792458.0  # m/s

@njit(cache=True)
def _image_type_name(code):
    """
    Return the integer image-type code unchanged.

    Convention: 1 = Type I (minimum), 2 = Type II (saddle),
    3 = Type III (maximum), 0 = undefined / degenerate.

    Parameters
    ----------
    code : ``int``
        Image-type code.

    Returns
    -------
    code : ``int``
        Same code, passed through.
    """
    return code


@njit(cache=True)
def _unique_points(x, y, n, tol):
    """
    Remove near-duplicate points in place.
    """
    keep = np.ones(n, dtype=np.uint8)

    for i in range(n):
        if keep[i] == 0:
            continue
        xi = x[i]
        yi = y[i]
        for j in range(i + 1, n):
            if keep[j] == 1:
                if abs(x[j] - xi) < tol and abs(y[j] - yi) < tol:
                    keep[j] = 0

    m = 0
    for i in range(n):
        if keep[i] == 1:
            x[m] = x[i]
            y[m] = y[i]
            m += 1
    return m


@njit(cache=True)
def _insertion_sort5(x, y, t, mu, itype, n):
    """
    Sort first n elements by ascending arrival time t.
    """
    for i in range(1, n):
        x0 = x[i]
        y0 = y[i]
        t0 = t[i]
        mu0 = mu[i]
        it0 = itype[i]

        j = i - 1
        while j >= 0 and t[j] > t0:
            x[j + 1] = x[j]
            y[j + 1] = y[j]
            t[j + 1] = t[j]
            mu[j + 1] = mu[j]
            itype[j + 1] = itype[j]
            j -= 1

        x[j + 1] = x0
        y[j + 1] = y0
        t[j + 1] = t0
        mu[j + 1] = mu0
        itype[j + 1] = it0


@njit(cache=True)
def _compress_keep5(mask, x, y, t, mu, itype, n):
    """
    In-place compression using boolean-like uint8 mask.
    """
    m = 0
    for i in range(n):
        if mask[i] == 1:
            x[m] = x[i]
            y[m] = y[i]
            t[m] = t[i]
            mu[m] = mu[i]
            itype[m] = itype[i]
            m += 1
    return m


@njit(cache=True, fastmath=True)
def _shear_function(x, y, gamma1, gamma2, ra_0=0.0, dec_0=0.0):
    """
    Compute the external shear lensing potential at a point.

    Parameters
    ----------
    x : ``float``
        Image-plane x-coordinate.
    y : ``float``
        Image-plane y-coordinate.
    gamma1 : ``float``
        External shear component 1.
    gamma2 : ``float``
        External shear component 2.
    ra_0 : ``float``
        Shear center x-coordinate.
    dec_0 : ``float``
        Shear center y-coordinate.

    Returns
    -------
    psi_shear : ``float``
        Shear potential value.
    """
    x_ = x - ra_0
    y_ = y - dec_0
    return 0.5 * (gamma1 * x_ * x_ + 2.0 * gamma2 * x_ * y_ - gamma1 * y_ * y_)


@njit(cache=True, fastmath=True)
def _shear_hessian(gamma1, gamma2):
    """
    Compute the Hessian matrix components of the external shear.

    Parameters
    ----------
    gamma1 : ``float``
        External shear component 1.
    gamma2 : ``float``
        External shear component 2.

    Returns
    -------
    f_xx : ``float``
        Second derivative d²ψ/dx².
    f_xy : ``float``
        Cross derivative d²ψ/dxdy.
    f_yx : ``float``
        Cross derivative d²ψ/dydx.
    f_yy : ``float``
        Second derivative d²ψ/dy².
    """
    f_xx = gamma1
    f_yy = -gamma1
    f_xy = gamma2
    return f_xx, f_xy, f_xy, f_yy


@njit(cache=True)
def _param_transform(x, y, theta_E, gamma, q, phi, center_x=0.0, center_y=0.0):
    """
    Convert lenstronomy-like EPL parameters into the Tessore+ convention.

    Parameters
    ----------
    x : ``float``
        Image-plane x-coordinate.
    y : ``float``
        Image-plane y-coordinate.
    theta_E : ``float``
        Einstein radius.
    gamma : ``float``
        EPL power-law slope.
    q : ``float``
        Axis ratio.
    phi : ``float``
        Lens position angle in radians.
    center_x : ``float``
        Lens center x-coordinate.
    center_y : ``float``
        Lens center y-coordinate.

    Returns
    -------
    z : ``complex``
        Rotated complex coordinate.
    b : ``float``
        EPL strength parameter.
    t : ``float``
        EPL slope exponent.
    q : ``float``
        Axis ratio.
    phi : ``float``
        Position angle.
    """
    t = gamma - 1.0
    x_shift = x - center_x
    y_shift = y - center_y
    z = np.exp(-1j * phi) * (x_shift + 1j * y_shift)
    b = theta_E * np.sqrt(q)
    return z, b, t, q, phi


@njit(cache=True)
def _alpha_epl(x, y, b, q, t):
    """
    Compute the complex EPL deflection in the rotated frame.

    Parameters
    ----------
    x : ``float``
        Rotated x-coordinate.
    y : ``float``
        Rotated y-coordinate.
    b : ``float``
        EPL strength parameter.
    q : ``float``
        Axis ratio.
    t : ``float``
        EPL slope exponent.

    Returns
    -------
    alph : ``complex``
        Complex deflection angle.
    """
    zz = x * q + 1j * y
    R = np.abs(zz)
    phi = np.angle(zz)
    Omega = omega_scalar(phi, t, q)

    if R < EPS:
        return 0.0 + 0.0j

    fac = (b / R) ** t * (R / b)
    return (2.0 * b / (1.0 + q)) * fac * Omega


@njit(cache=True)
def _epl_function(x, y, theta_E, gamma, q, phi, center_x=0.0, center_y=0.0):
    """
    Compute the EPL lensing potential at a point.

    Parameters
    ----------
    x : ``float``
        Image-plane x-coordinate.
    y : ``float``
        Image-plane y-coordinate.
    theta_E : ``float``
        Einstein radius.
    gamma : ``float``
        EPL power-law slope.
    q : ``float``
        Axis ratio.
    phi : ``float``
        Lens position angle in radians.
    center_x : ``float``
        Lens center x-coordinate.
    center_y : ``float``
        Lens center y-coordinate.

    Returns
    -------
    psi : ``float``
        EPL potential value.
    """
    z, b, t, q_, _ = _param_transform(x, y, theta_E, gamma, q, phi, center_x, center_y)
    alph = _alpha_epl(z.real, z.imag, b, q_, t)
    return (z.real * alph.real + z.imag * alph.imag) / (2.0 - t)


@njit(cache=True)
def _epl_hessian(x, y, theta_E, gamma, q, phi, center_x=0.0, center_y=0.0):
    """
    Compute the EPL Hessian components at a point.

    Parameters
    ----------
    x : ``float``
        Image-plane x-coordinate.
    y : ``float``
        Image-plane y-coordinate.
    theta_E : ``float``
        Einstein radius.
    gamma : ``float``
        EPL power-law slope.
    q : ``float``
        Axis ratio.
    phi : ``float``
        Lens position angle in radians.
    center_x : ``float``
        Lens center x-coordinate.
    center_y : ``float``
        Lens center y-coordinate.

    Returns
    -------
    f_xx : ``float``
        Second derivative d²ψ/dx².
    f_xy : ``float``
        Cross derivative d²ψ/dxdy.
    f_yx : ``float``
        Cross derivative d²ψ/dydx.
    f_yy : ``float``
        Second derivative d²ψ/dy².
    """
    z, b, t, q_, ang_ell = _param_transform(x, y, theta_E, gamma, q, phi, center_x, center_y)

    ang = np.angle(z)
    zz_ell = z.real * q_ + 1j * z.imag
    R = np.abs(zz_ell)
    phi_ell = np.angle(zz_ell)

    if R > EPS:
        u = (b / R) ** t
        if u > 1.0e10:
            u = 1.0e10
    else:
        u = 1.0e10

    kappa = 0.5 * (2.0 - t)
    Roverr = np.sqrt((np.cos(ang) ** 2) * q_ * q_ + np.sin(ang) ** 2)

    Omega = omega_scalar(phi_ell, t, q_)
    alph = (2.0 / (1.0 + q_)) * Omega

    gamma_shear = (
        -np.exp(2j * (ang + ang_ell)) * kappa
        + (1.0 - t) * np.exp(1j * (ang + 2.0 * ang_ell)) * alph * Roverr
    )

    f_xx = (kappa + gamma_shear.real) * u
    f_yy = (kappa - gamma_shear.real) * u
    f_xy = gamma_shear.imag * u

    return f_xx, f_xy, f_xy, f_yy


@njit(cache=True)
def _potential_scalar(x, y, theta_E, gamma, gamma1, gamma2, q, phi,
                     center_x=0.0, center_y=0.0, ra_0=0.0, dec_0=0.0):
    """
    Compute the total lensing potential (EPL + external shear).

    Parameters
    ----------
    x : ``float``
        Image-plane x-coordinate.
    y : ``float``
        Image-plane y-coordinate.
    theta_E : ``float``
        Einstein radius.
    gamma : ``float``
        EPL power-law slope.
    gamma1 : ``float``
        External shear component 1.
    gamma2 : ``float``
        External shear component 2.
    q : ``float``
        Axis ratio.
    phi : ``float``
        Lens position angle in radians.
    center_x : ``float``
        Lens center x-coordinate.
    center_y : ``float``
        Lens center y-coordinate.
    ra_0 : ``float``
        Shear center x-coordinate.
    dec_0 : ``float``
        Shear center y-coordinate.

    Returns
    -------
    psi : ``float``
        Total potential value.
    """
    return (
        _epl_function(x, y, theta_E, gamma, q, phi, center_x, center_y)
        + _shear_function(x, y, gamma1, gamma2, ra_0, dec_0)
    )


@njit(cache=True)
def _hessian_scalar(x, y, theta_E, gamma, gamma1, gamma2, q, phi,
                   center_x=0.0, center_y=0.0):
    """
    Compute the total Hessian (EPL + external shear).

    Parameters
    ----------
    x : ``float``
        Image-plane x-coordinate.
    y : ``float``
        Image-plane y-coordinate.
    theta_E : ``float``
        Einstein radius.
    gamma : ``float``
        EPL power-law slope.
    gamma1 : ``float``
        External shear component 1.
    gamma2 : ``float``
        External shear component 2.
    q : ``float``
        Axis ratio.
    phi : ``float``
        Lens position angle in radians.
    center_x : ``float``
        Lens center x-coordinate.
    center_y : ``float``
        Lens center y-coordinate.

    Returns
    -------
    f_xx : ``float``
        Total d²ψ/dx².
    f_xy : ``float``
        Total d²ψ/dxdy.
    f_yx : ``float``
        Total d²ψ/dydx.
    f_yy : ``float``
        Total d²ψ/dy².
    """
    f_xx_e, f_xy_e, f_yx_e, f_yy_e = _epl_hessian(
        x, y, theta_E, gamma, q, phi, center_x, center_y
    )
    f_xx_s, f_xy_s, f_yx_s, f_yy_s = _shear_hessian(gamma1, gamma2)

    return (
        f_xx_e + f_xx_s,
        f_xy_e + f_xy_s,
        f_yx_e + f_yx_s,
        f_yy_e + f_yy_s,
    )


@njit(cache=True)
def lensing_diagnostics_scalar(x, y, theta_E, gamma, gamma1, gamma2, q, phi,
                               center_x=0.0, center_y=0.0):
    """
    Compute all Hessian-derived local diagnostics at one image position.

    Parameters
    ----------
    x : ``float``
        Image-plane x-coordinate.
    y : ``float``
        Image-plane y-coordinate.
    theta_E : ``float``
        Einstein radius.
    gamma : ``float``
        EPL power-law slope.
    gamma1 : ``float``
        External shear component 1.
    gamma2 : ``float``
        External shear component 2.
    q : ``float``
        Axis ratio.
    phi : ``float``
        Lens position angle in radians.
    center_x : ``float``
        Lens center x-coordinate.
    center_y : ``float``
        Lens center y-coordinate.

    Returns
    -------
    f_xx : ``float``
        Hessian component d²ψ/dx².
    f_xy : ``float``
        Hessian component d²ψ/dxdy.
    f_yx : ``float``
        Hessian component d²ψ/dydx.
    f_yy : ``float``
        Hessian component d²ψ/dy².
    detA : ``float``
        Jacobian determinant.
    traceA : ``float``
        Jacobian trace.
    mu : ``float``
        Signed magnification.
    image_type : ``int``
        Image classification. \n
        1 = Type I (minimum), 2 = Type II (saddle), \n
        3 = Type III (maximum), 0 = undefined.

    Examples
    --------
    >>> f_xx, f_xy, f_yx, f_yy, detA, traceA, mu, itype = lensing_diagnostics_scalar(
    ...     x=0.5, y=0.3, theta_E=1.0, gamma=2.0,
    ...     gamma1=0.03, gamma2=-0.01, q=0.8, phi=0.3
    ... )
    """
    f_xx, f_xy, f_yx, f_yy = _hessian_scalar(
        x, y, theta_E, gamma, gamma1, gamma2, q, phi, center_x, center_y
    )

    detA = (1.0 - f_xx) * (1.0 - f_yy) - f_xy * f_yx
    traceA = 2.0 - f_xx - f_yy

    if abs(detA) < EPS:
        mu = np.inf
        image_type = 0
    else:
        mu = 1.0 / detA

        if detA < 0.0:
            image_type = 2
        elif traceA > 0.0:
            image_type = 1
        elif traceA < 0.0:
            image_type = 3
        else:
            image_type = 0

    return f_xx, f_xy, f_yx, f_yy, detA, traceA, mu, image_type


@njit(cache=True)
def fermat_potential_scalar(x_image, y_image, x_source, y_source,
                            theta_E, gamma, gamma1, gamma2, q, phi,
                            center_x=0.0, center_y=0.0, ra_0=0.0, dec_0=0.0):
    """
    Compute the Fermat potential (geometric delay minus lensing potential).

    Parameters
    ----------
    x_image : ``float``
        Image-plane x-coordinate.
    y_image : ``float``
        Image-plane y-coordinate.
    x_source : ``float``
        Source-plane x-coordinate.
    y_source : ``float``
        Source-plane y-coordinate.
    theta_E : ``float``
        Einstein radius.
    gamma : ``float``
        EPL power-law slope.
    gamma1 : ``float``
        External shear component 1.
    gamma2 : ``float``
        External shear component 2.
    q : ``float``
        Axis ratio.
    phi : ``float``
        Lens position angle in radians.
    center_x : ``float``
        Lens center x-coordinate.
    center_y : ``float``
        Lens center y-coordinate.
    ra_0 : ``float``
        Shear center x-coordinate.
    dec_0 : ``float``
        Shear center y-coordinate.

    Returns
    -------
    tau : ``float``
        Fermat potential value.

    Examples
    --------
    >>> tau = fermat_potential_scalar(
    ...     x_image=0.5, y_image=0.3, x_source=0.1, y_source=0.05,
    ...     theta_E=1.0, gamma=2.0, gamma1=0.03, gamma2=-0.01, q=0.8, phi=0.3
    ... )
    """
    pot = _potential_scalar(
        x_image, y_image, theta_E, gamma, gamma1, gamma2, q, phi,
        center_x, center_y, ra_0, dec_0
    )
    geom = 0.5 * ((x_image - x_source) ** 2 + (y_image - y_source) ** 2)
    return geom - pot


@njit(cache=True)
def _geomlinspace(a, b, N):
    """
    Generate a hybrid linear-geometric spacing array.

    Parameters
    ----------
    a : ``float``
        Start value.
    b : ``float``
        End value.
    N : ``int``
        Number of geometric points.

    Returns
    -------
    out : ``numpy.ndarray``
        Combined linear + geometric spacing array.
    """
    ratio = (b / a) ** (1.0 / (N - 1))
    delta = a * (ratio - 1.0)

    nlin = int(a / delta)
    if nlin < 1:
        nlin = 1
    if nlin > N:
        nlin = N

    out = np.empty(nlin + N, dtype=np.float64)

    for i in range(nlin):
        out[i] = a * i / nlin

    for i in range(N):
        out[nlin + i] = a * (ratio ** i)

    return out


@njit(cache=True)
def _ps(x, p):
    """
    Compute the signed power function ``|x|^p * sign(x)``.

    Parameters
    ----------
    x : ``float``
        Input value.
    p : ``float``
        Exponent.

    Returns
    -------
    result : ``float``
        Signed power.
    """
    return np.abs(x) ** p * np.sign(x)


@njit(cache=True)
def _ell_to_pol(rell, theta, q):
    """
    Convert elliptical coordinates to polar coordinates.

    Parameters
    ----------
    rell : ``float``
        Elliptical radial coordinate.
    theta : ``float``
        Elliptical angle in radians.
    q : ``float``
        Axis ratio.

    Returns
    -------
    r : ``float``
        Polar radial coordinate.
    phi : ``float``
        Polar angle in radians.
    """
    phi = np.arctan2(np.sin(theta) * q, np.cos(theta))
    r = rell * np.sqrt((np.cos(theta) ** 2) / (q * q) + np.sin(theta) ** 2)
    return r, phi


@njit(cache=True)
def _one_dim_lens_eq_calcs(args, phi):
    """
    Intermediate quantities for the 1-D lens equation.
    """
    b, t, y1, y2, q, gamma1, gamma2 = args
    y = y1 + 1j * y2

    den = 1.0 - gamma1 * gamma1 - gamma2 * gamma2
    if abs(den) < EPS:
        den = EPS if den >= 0.0 else -EPS

    c = np.cos(phi)
    s = np.sin(phi)

    rhat = (
        ((1.0 + gamma1) * c + gamma2 * s)
        + 1j * (gamma2 * c + (1.0 - gamma1) * s)
    ) / den

    thetahat = (
        ((1.0 + gamma1) * s - gamma2 * c)
        + 1j * (gamma2 * s - (1.0 - gamma1) * c)
    ) / den

    frac_Roverrsh, phiell = pol_to_ell(1.0, phi, q)
    Omega = omega_scalar(phiell, t, q)
    const = 2.0 * b / (1.0 + q)

    if abs(t - 1.0) > 1e-4:
        denom = const * cdot(Omega, thetahat)
        if abs(denom) < EPS:
            denom = EPS if denom >= 0.0 else -EPS

        val = -cdot(y, thetahat) / denom
        aval = abs(val)
        if aval < EPS:
            R = 0.0
        else:
            R = b * aval ** (1.0 / (1.0 - t)) * np.sign(val)
    else:
        Omega_ort = 1j * Omega
        x = ((1.0 - gamma1) * c - gamma2 * s) + 1j * (
            -gamma2 * c + (1.0 + gamma1) * s
        )
        denom = cdot(Omega_ort, x)
        if abs(denom) < EPS:
            denom = EPS if denom >= 0.0 else -EPS
        R = cdot(Omega_ort, y) / denom * frac_Roverrsh

    r, theta = _ell_to_pol(R, phiell, q)
    return Omega, const, phiell, q, r, rhat, t, b, thetahat, y


@njit(cache=True)
def _one_dim_lens_eq(phi, args):
    """
    Smooth branch of the 1-D lens equation.
    """
    Omega, const, phiell, q, r, rhat, t, b, thetahat, y = _one_dim_lens_eq_calcs(args, phi)
    rr, _ = _ell_to_pol(1.0, phiell, q)

    ip = cdot(y, rhat) * cdot(Omega, thetahat) - cdot(Omega, rhat) * cdot(y, thetahat)

    eq = (rr * b) ** (2.0 / t - 2.0) * _ps(cdot(y, thetahat) / const, 2.0 / t) * ip * ip
    eq += _ps(ip, 2.0 / t) * (cdot(Omega, thetahat) ** 2)
    return eq


@njit(cache=True)
def _one_dim_lens_eq_unsmooth(phi, args):
    """
    Non-smooth branch of the 1-D lens equation.
    """
    Omega, const, phiell, q, r, rhat, t, b, thetahat, y = _one_dim_lens_eq_calcs(args, phi)
    rr, _ = _ell_to_pol(1.0, phiell, q)

    ip = cdot(y, rhat) * cdot(Omega, thetahat) - cdot(Omega, rhat) * cdot(y, thetahat)

    eq = _ps(rr * b, 1.0 - t) * (cdot(y, thetahat) / const) * (np.abs(ip) ** t)
    eq += ip * (np.abs(cdot(Omega, thetahat)) ** t)
    return eq


@njit(cache=True)
def _one_dim_lens_eq_both(phi, args):
    n = phi.shape[0]
    y = np.empty(n, dtype=np.float64)
    y_ns = np.empty(n, dtype=np.float64)

    for i in range(n):
        y[i] = _one_dim_lens_eq(phi[i], args)
        y_ns[i] = _one_dim_lens_eq_unsmooth(phi[i], args)

    return y, y_ns


@njit(cache=True)
def _min_approx(x1, x2, x3, y1, y2, y3):
    """
    Estimate the location of a local extremum via parabolic interpolation.

    Parameters
    ----------
    x1 : ``float``
        First x-coordinate.
    x2 : ``float``
        Second x-coordinate.
    x3 : ``float``
        Third x-coordinate.
    y1 : ``float``
        Function value at x1.
    y2 : ``float``
        Function value at x2.
    y3 : ``float``
        Function value at x3.

    Returns
    -------
    x_min : ``float``
        Estimated extremum location.
    """
    div = 2.0 * (x3 * (y1 - y2) + x1 * (y2 - y3) + x2 * (-y1 + y3))
    if abs(div) < EPS:
        return x2
    num = x3 * x3 * (y1 - y2) + x1 * x1 * (y2 - y3) + x2 * x2 * (-y1 + y3)
    return num / div


@njit(cache=True)
def _brentq_inline(f, xa, xb, xtol=2e-14, rtol=16 * np.finfo(np.float64).eps, maxiter=100, args=()):
    """
    Numba-compatible Brent root finder.

    Returns ``np.nan`` if the bracket is invalid or numerically unstable.

    Parameters
    ----------
    f : ``callable``
        Scalar function of ``(x, args)``.
    xa : ``float``
        Left bracket.
    xb : ``float``
        Right bracket.
    xtol : ``float``
        Absolute tolerance.
    rtol : ``float``
        Relative tolerance.
    maxiter : ``int``
        Maximum iterations.
    args : ``tuple``
        Extra arguments forwarded to ``f``.

    Returns
    -------
    root : ``float``
        Approximate root, or ``np.nan`` on failure.
    """
    xpre = xa
    xcur = xb
    xblk = 0.0
    fblk = 0.0
    spre = 0.0
    scur = 0.0

    fpre = f(xpre, args)
    fcur = f(xcur, args)

    if np.isnan(fpre) or np.isnan(fcur):
        return np.nan

    if fpre * fcur > 0.0:
        return np.nan
    if fpre == 0.0:
        return xpre
    if fcur == 0.0:
        return xcur

    for _ in range(maxiter):
        if fpre * fcur < 0.0:
            xblk = xpre
            fblk = fpre

        if abs(fblk) < abs(fcur):
            xpre = xcur
            xcur = xblk
            xblk = xpre

            fpre = fcur
            fcur = fblk
            fblk = fpre

        delta = 0.5 * (xtol + rtol * abs(xcur))
        sbis = 0.5 * (xblk - xcur)

        if fcur == 0.0 or abs(sbis) < delta:
            return xcur

        use_bisect = True
        stry = 0.0

        if abs(spre) > delta and abs(fcur) < abs(fpre):
            if xpre == xblk:
                denom = fcur - fpre
                if abs(denom) >= EPS:
                    stry = -fcur * (xcur - xpre) / denom
                    use_bisect = False
            else:
                dpre_den = xpre - xcur
                dblk_den = xblk - xcur
                if abs(dpre_den) >= EPS and abs(dblk_den) >= EPS:
                    dpre = (fpre - fcur) / dpre_den
                    dblk = (fblk - fcur) / dblk_den
                    denom = dblk * dpre * (fblk - fpre)
                    if abs(denom) >= EPS:
                        stry = -fcur * (fblk * dblk - fpre * dpre) / denom
                        use_bisect = False

            if (not use_bisect) and (2.0 * abs(stry) < min(abs(spre), 3.0 * abs(sbis) - delta)):
                spre = scur
                scur = stry
            else:
                spre = sbis
                scur = sbis
        else:
            spre = sbis
            scur = sbis

        xpre = xcur
        fpre = fcur

        if abs(scur) > delta:
            xcur += scur
        else:
            xcur += delta if sbis > 0.0 else -delta

        fcur = f(xcur, args)
        if np.isnan(fcur):
            return np.nan

    return xcur


@njit(cache=True)
def _getr(phi, args):
    """
    Compute radius for a given image-plane angle phi.
    """
    _, _, _, _, r, _, _, _, _, _ = _one_dim_lens_eq_calcs(args, phi)
    return r


@njit(cache=True)
def _getphi(thpl, args):
    """
    Find all angular roots of the 1-D lens equation on the supplied grid.
    """
    y, y_ns = _one_dim_lens_eq_both(thpl, args)
    nphi = thpl.shape[0]

    roots = np.empty(MAX_ROOTS, dtype=np.float64)
    nroots = 0

    for i in range(nphi - 1):
        if y[i + 1] * y[i] <= 0.0:
            r = _brentq_inline(_one_dim_lens_eq, thpl[i], thpl[i + 1], args=args)
            if (not np.isnan(r)) and nroots < MAX_ROOTS:
                roots[nroots] = r % (2.0 * np.pi)
                nroots += 1
        elif y_ns[i + 1] * y_ns[i] <= 0.0:
            r = _brentq_inline(_one_dim_lens_eq_unsmooth, thpl[i], thpl[i + 1], args=args)
            if (not np.isnan(r)) and nroots < MAX_ROOTS:
                roots[nroots] = r % (2.0 * np.pi)
                nroots += 1

    for i in range(1, nphi):
        y1 = y[i - 1]
        y2 = y[i]
        y3 = y[i + 1] if i + 1 < nphi else y[0]

        y1n = y_ns[i - 1]
        y2n = y_ns[i]
        y3n = y_ns[i + 1] if i + 1 < nphi else y_ns[0]

        if (y3 - y2) * (y2 - y1) <= 0.0 or (y3n - y2n) * (y2n - y1n) <= 0.0:
            if y3 * y2 <= 0.0 or y1 * y2 <= 0.0:
                continue
            if i > nphi - 2:
                continue

            x1 = thpl[i - 1]
            x2 = thpl[i]
            x3 = thpl[i + 1]

            xmin = _min_approx(x1, x2, x3, y1, y2, y3)
            xmin_ns = _min_approx(x1, x2, x3, y1n, y2n, y3n)

            ymin = _one_dim_lens_eq(xmin, args)
            ymin_ns = _one_dim_lens_eq_unsmooth(xmin_ns, args)

            if ymin * y2 <= 0.0 and x2 <= xmin <= x3:
                if nroots < MAX_ROOTS:
                    r = _brentq_inline(_one_dim_lens_eq, x2, xmin, args=args)
                    if not np.isnan(r):
                        roots[nroots] = r % (2.0 * np.pi)
                        nroots += 1
                if nroots < MAX_ROOTS:
                    r = _brentq_inline(_one_dim_lens_eq, xmin, x3, args=args)
                    if not np.isnan(r):
                        roots[nroots] = r % (2.0 * np.pi)
                        nroots += 1

            elif ymin * y2 <= 0.0 and x1 <= xmin <= x2:
                if nroots < MAX_ROOTS:
                    r = _brentq_inline(_one_dim_lens_eq, x1, xmin, args=args)
                    if not np.isnan(r):
                        roots[nroots] = r % (2.0 * np.pi)
                        nroots += 1
                if nroots < MAX_ROOTS:
                    r = _brentq_inline(_one_dim_lens_eq, xmin, x2, args=args)
                    if not np.isnan(r):
                        roots[nroots] = r % (2.0 * np.pi)
                        nroots += 1

            elif ymin_ns * y2n <= 0.0 and x2 <= xmin_ns <= x3:
                if nroots < MAX_ROOTS:
                    r = _brentq_inline(_one_dim_lens_eq_unsmooth, x2, xmin_ns, args=args)
                    if not np.isnan(r):
                        roots[nroots] = r % (2.0 * np.pi)
                        nroots += 1
                if nroots < MAX_ROOTS:
                    r = _brentq_inline(_one_dim_lens_eq_unsmooth, xmin_ns, x3, args=args)
                    if not np.isnan(r):
                        roots[nroots] = r % (2.0 * np.pi)
                        nroots += 1

            elif ymin_ns * y2n <= 0.0 and x1 <= xmin_ns <= x2:
                if nroots < MAX_ROOTS:
                    r = _brentq_inline(_one_dim_lens_eq_unsmooth, x1, xmin_ns, args=args)
                    if not np.isnan(r):
                        roots[nroots] = r % (2.0 * np.pi)
                        nroots += 1
                if nroots < MAX_ROOTS:
                    r = _brentq_inline(_one_dim_lens_eq_unsmooth, xmin_ns, x2, args=args)
                    if not np.isnan(r):
                        roots[nroots] = r % (2.0 * np.pi)
                        nroots += 1

    return roots, nroots


@njit(cache=True)
def _solvelenseq_majoraxis(b, t, y1, y2, q, gamma1, gamma2, Nmeas=200, Nmeas_extra=50):
    """
    Solve the lens equation in the major-axis frame.

    Parameters
    ----------
    b : ``float``
        EPL strength parameter.
    t : ``float``
        EPL slope exponent.
    y1 : ``float``
        Source x-coordinate in the rotated frame.
    y2 : ``float``
        Source y-coordinate in the rotated frame.
    q : ``float``
        Axis ratio.
    gamma1 : ``float``
        Rotated shear component 1.
    gamma2 : ``float``
        Rotated shear component 2.
    Nmeas : ``int``
        Number of angular grid points.
    Nmeas_extra : ``int``
        Extra refinement points near the source angle.

    Returns
    -------
    xsol : ``numpy.ndarray``
        Image x-coordinates.
    ysol : ``numpy.ndarray``
        Image y-coordinates.
    nimg : ``int``
        Number of images found.
    """
    p1 = np.arctan2(
        y2 * (1.0 - gamma1) + gamma2 * y1,
        y1 * (1.0 + gamma1) + gamma2 * y2
    )

    geom = _geomlinspace(1e-4, 0.1, Nmeas_extra)
    ngeom = geom.shape[0]
    ntot = Nmeas + 2 * ngeom

    thpl = np.empty(ntot, dtype=np.float64)

    for i in range(Nmeas):
        thpl[i] = np.pi * i / (Nmeas - 1)

    pmod = p1 % np.pi
    for i in range(ngeom):
        thpl[Nmeas + i] = pmod - geom[i]
        thpl[Nmeas + ngeom + i] = pmod + geom[i]

    thpl.sort()

    args = (b, t, y1, y2, q, gamma1, gamma2)
    roots, nroots = _getphi(thpl, args)

    xsol = np.empty(MAX_IMGS, dtype=np.float64)
    ysol = np.empty(MAX_IMGS, dtype=np.float64)
    nimg = 0

    for i in range(nroots):
        for add_pi in range(2):
            th = roots[i] + add_pi * np.pi
            R = _getr(th, args)

            if R <= 0.0:
                continue

            xs = R * np.cos(th)
            ys = R * np.sin(th)

            diff = (-y1 - 1j * y2) + (xs + 1j * ys) - _alpha_epl_shear_scalar(
                xs, ys, b, q, t, gamma1=gamma1, gamma2=gamma2
            )

            if abs(diff) < 1e-8:
                if nimg < MAX_IMGS:
                    xsol[nimg] = xs
                    ysol[nimg] = ys
                    nimg += 1

    nimg = _unique_points(xsol, ysol, nimg, 1e-8)
    return xsol, ysol, nimg


@njit(cache=True)
def _solve_lenseq_pemd(xsrc, ysrc, q, phi, gamma, gamma1, gamma2,
                      theta_E=1.0, Nmeas=400, Nmeas_extra=80):
    """
    Solve the lens equation in the rotated major-axis frame and rotate back.

    Parameters
    ----------
    xsrc : ``float``
        Source x-coordinate.
    ysrc : ``float``
        Source y-coordinate.
    q : ``float``
        Axis ratio.
    phi : ``float``
        Lens position angle in radians.
    gamma : ``float``
        EPL power-law slope.
    gamma1 : ``float``
        External shear component 1.
    gamma2 : ``float``
        External shear component 2.
    theta_E : ``float``
        Einstein radius.
    Nmeas : ``int``
        Angular grid size.
    Nmeas_extra : ``int``
        Extra refinement points.

    Returns
    -------
    xout : ``numpy.ndarray``
        Image x-coordinates in the original frame.
    yout : ``numpy.ndarray``
        Image y-coordinates in the original frame.
    nimg : ``int``
        Number of images found.
    """
    t = gamma - 1.0
    b = theta_E * np.sqrt(q)

    rot = np.exp(-1j * phi)
    rot2 = rot * rot

    p = (xsrc + 1j * ysrc) * rot
    g = (gamma1 + 1j * gamma2) * rot2

    xloc, yloc, nimg = _solvelenseq_majoraxis(
        b, t, p.real, p.imag, q, g.real, g.imag, Nmeas, Nmeas_extra
    )

    xout = np.empty(MAX_IMGS, dtype=np.float64)
    yout = np.empty(MAX_IMGS, dtype=np.float64)

    invrot = np.exp(1j * phi)
    for i in range(nimg):
        z = (xloc[i] + 1j * yloc[i]) * invrot
        xout[i] = z.real
        yout[i] = z.imag

    return xout, yout, nimg


@njit(cache=True)
def image_position_analytical_njit(
    x,
    y,
    q,
    phi,
    gamma,
    gamma1,
    gamma2,
    theta_E=1.0,
    alpha_scaling=1.0,
    magnification_limit=0.01,
    Nmeas=400,
    Nmeas_extra=80,
):
    """
    Standalone EPL + shear analytical image finder.

    Locates lensed images for a given source position, computes their
    magnifications, Fermat potential (arrival time proxy), and image
    types. Results are sorted by ascending arrival time and filtered
    by a minimum magnification threshold.

    Parameters
    ----------
    x : ``float``
        Source-plane x-coordinate.
    y : ``float``
        Source-plane y-coordinate.
    q : ``float``
        Lens axis ratio.
    phi : ``float``
        Lens position angle in radians.
    gamma : ``float``
        EPL power-law slope (lenstronomy convention). \n
        The Tessore exponent is ``t = gamma - 1``.
    gamma1 : ``float``
        External shear component 1.
    gamma2 : ``float``
        External shear component 2.
    theta_E : ``float``
        Einstein radius. \n
        default: 1.0
    alpha_scaling : ``float``
        Deflection scaling factor applied to theta_E. \n
        default: 1.0
    magnification_limit : ``float``
        Minimum |mu| threshold; images below this are discarded. \n
        default: 0.01
    Nmeas : ``int``
        Angular root-finding grid size. \n
        default: 400
    Nmeas_extra : ``int``
        Extra refinement points near the source angle. \n
        default: 80

    Returns
    -------
    x_img : ``numpy.ndarray``
        Image x-positions sorted by arrival time.
    y_img : ``numpy.ndarray``
        Image y-positions sorted by arrival time.
    arrival_time : ``numpy.ndarray``
        Fermat potential values in increasing order.
    magnification : ``numpy.ndarray``
        Signed magnifications.
    image_type : ``numpy.ndarray``
        Image classification codes (int64). \n
        1 = Type I (minimum), 2 = Type II (saddle), \n
        3 = Type III (maximum), 0 = undefined.
    nimg : ``int``
        Number of images returned.

    Examples
    --------
    >>> x_img, y_img, tau, mu, itype, n = image_position_analytical_njit(
    ...     x=0.1, y=0.05, q=0.8, phi=0.3, gamma=2.0,
    ...     gamma1=0.03, gamma2=-0.01
    ... )
    """
    theta_E_eff = theta_E * alpha_scaling ** (1.0 / (gamma - 1.0))

    x_img, y_img, nimg = _solve_lenseq_pemd(
        x, y, q, phi, gamma, gamma1, gamma2,
        theta_E=theta_E_eff,
        Nmeas=Nmeas,
        Nmeas_extra=Nmeas_extra,
    )

    arrival_time = np.empty(MAX_IMGS, dtype=np.float64)
    magnification = np.empty(MAX_IMGS, dtype=np.float64)
    image_type = np.empty(MAX_IMGS, dtype=np.int64)

    for i in range(nimg):
        _, _, _, _, _, _, mu, itype = lensing_diagnostics_scalar(
            x_img[i], y_img[i],
            theta_E_eff, gamma, gamma1, gamma2, q, phi
        )

        magnification[i] = mu
        image_type[i] = itype

        arrival_time[i] = fermat_potential_scalar(
            x_img[i], y_img[i], x, y,
            theta_E_eff, gamma, gamma1, gamma2, q, phi
        )

    _insertion_sort5(x_img, y_img, arrival_time, magnification, image_type, nimg)

    keep = np.zeros(MAX_IMGS, dtype=np.uint8)
    for i in range(nimg):
        if np.isinf(magnification[i]) or abs(magnification[i]) >= magnification_limit:
            keep[i] = 1

    nimg = _compress_keep5(
        keep, x_img, y_img, arrival_time, magnification, image_type, nimg
    )

    return (
        x_img[:nimg],
        y_img[:nimg],
        arrival_time[:nimg],
        magnification[:nimg],
        image_type[:nimg],
        nimg,
    )

def create_epl_shear_solver(
    arrival_time_sort=True,
    max_img=4,
    num_th=500,
    maginf=-100.0,
    max_tries=100,
    alpha_scaling=1.0,
    magnification_limit=0.01,
    Nmeas=400,
    Nmeas_extra=80,
):
    """
    Create a parallel EPL + shear solver for batched lens systems.

    Returns a JIT-compiled function that, for each system in a batch,
    samples a source from the double caustic and solves for image
    positions, magnifications, time delays, and image types.

    Parameters
    ----------
    arrival_time_sort : ``bool``
        Whether to sort images by arrival time. \n
        default: True
    max_img : ``int``
        Maximum number of images to store per system. \n
        default: 4
    num_th : ``int``
        Angular samples for caustic construction. \n
        default: 500
    maginf : ``float``
        Magnification cutoff for caustic boundary. \n
        default: -100.0
    max_tries : ``int``
        Maximum rejection-sampling attempts per source. \n
        default: 100
    alpha_scaling : ``float``
        Deflection scaling factor. \n
        default: 1.0
    magnification_limit : ``float``
        Minimum |mu| threshold for image retention. \n
        default: 0.01
    Nmeas : ``int``
        Angular root-finding grid size. \n
        default: 400
    Nmeas_extra : ``int``
        Extra refinement points. \n
        default: 80

    Returns
    -------
    solve_epl_shear_multithreaded : ``callable``
        Parallel solver function with signature \n
        ``(theta_E, D_dt, q, phi, gamma, gamma1, gamma2)`` \n
        returning a tuple of result arrays.

    Examples
    --------
    >>> solver = create_epl_shear_solver()
    >>> results = solver(theta_E, D_dt, q, phi, gamma, gamma1, gamma2)
    """

    @njit(parallel=True, cache=True)
    def solve_epl_shear_multithreaded(theta_E, D_dt, q, phi, gamma, gamma1, gamma2):
        """
        Solve EPL + shear lens systems in parallel for a batch of parameters.

        Parameters
        ----------
        theta_E : ``numpy.ndarray``
            Einstein radii.
        D_dt : ``numpy.ndarray``
            Time-delay distances in metres.
        q : ``numpy.ndarray``
            Lens axis ratios.
        phi : ``numpy.ndarray``
            Lens position angles in radians.
        gamma : ``numpy.ndarray``
            EPL slopes.
        gamma1 : ``numpy.ndarray``
            External shear component 1.
        gamma2 : ``numpy.ndarray``
            External shear component 2.

        Returns
        -------
        beta_x_arr : ``numpy.ndarray``
            Source x-coordinates.
        beta_y_arr : ``numpy.ndarray``
            Source y-coordinates.
        x_img : ``numpy.ndarray``
            Image x-coordinates (size × max_img).
        y_img : ``numpy.ndarray``
            Image y-coordinates (size × max_img).
        mu_arr : ``numpy.ndarray``
            Signed magnifications (size × max_img).
        tau_arr : ``numpy.ndarray``
            Physical time delays in seconds (size × max_img).
        nimg : ``numpy.ndarray``
            Number of images per system.
        itype : ``numpy.ndarray``
            Image type codes (size × max_img).
        """
        size = theta_E.size

        # result arrays
        nimg = np.zeros(size, dtype=np.int64)
        x_img = np.full((size, max_img), np.nan)
        y_img = np.full((size, max_img), np.nan)
        mu_arr = np.full((size, max_img), np.nan)
        tau_arr = np.full((size, max_img), np.nan)
        itype = np.zeros((size, max_img), dtype=np.int64)
        beta_x_arr = np.full((size,), np.nan)
        beta_y_arr = np.full((size,), np.nan)

        ok_arr = np.zeros(size, dtype=np.uint8)

        # find source position for each lens system and sample caustic points
        for i in prange(size):
            beta, pts = sample_source_from_double_caustic(
                q[i],
                phi[i],
                gamma[i],
                gamma1[i],
                gamma2[i],
                theta_E=1.0,
                num_th=num_th,
                maginf=maginf,
                max_tries=max_tries,
            )
            beta_x, beta_y, ok = beta

            ok_arr[i] = ok
            beta_x_arr[i] = beta_x
            beta_y_arr[i] = beta_y

        # solve for image positions in systems with valid source samples
        idx = np.where(ok_arr == 1)[0]

        for i in prange(idx.size):

            idx_i = idx[i]

            theta_E_val = theta_E[idx_i]

            x, y, arrival_time, mu, image_type, n = image_position_analytical_njit(
                x=beta_x_arr[idx_i], 
                y=beta_y_arr[idx_i], 
                q=q[idx_i], 
                phi=phi[idx_i], 
                gamma=gamma[idx_i], 
                gamma1=gamma1[idx_i], 
                gamma2=gamma2[idx_i],
                theta_E=1.0,
                alpha_scaling=alpha_scaling,
                magnification_limit=magnification_limit,
                Nmeas=Nmeas,
                Nmeas_extra=Nmeas_extra,
            )

            nimg[idx_i] = n
            mu_arr[idx_i, :n] = mu[:n]

            # rescale to physical units
            x_img[idx_i, :n] = x[:n] * theta_E_val
            y_img[idx_i, :n] = y[:n] * theta_E_val
            beta_x_arr[idx_i] = beta_x_arr[idx_i] * theta_E_val
            beta_y_arr[idx_i] = beta_y_arr[idx_i] * theta_E_val

            # rescale time delays
            tau_hat = arrival_time[:n]
            # Fermat potential rescales as theta_E^2
            tau_phys = tau_hat * (theta_E_val * theta_E_val)
            # physical time delays in days: Δt = (D_dt/c) * Δtau
            tau_arr[idx_i, :n] = (D_dt[idx_i] / C_LIGHT) * tau_phys

            # image type
            itype[idx_i, :n] = image_type[:n]

        return (
            beta_x_arr,
            beta_y_arr,
            x_img,
            y_img,
            mu_arr,
            tau_arr,
            nimg,
            itype,
        )

    return solve_epl_shear_multithreaded
