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
def _shear_cartesian2polar(gamma1, gamma2):
    """
    Convert Cartesian shear components to polar form.

    Parameters
    ----------
    gamma1 : ``float``
        First Cartesian shear component.
    gamma2 : ``float``
        Second Cartesian shear component.

    Returns
    -------
    phi : ``float``
        Shear position angle in radians.
    gamma : ``float``
        Shear magnitude.
    """
    phi = np.arctan2(gamma2, gamma1) / 2
    gamma = np.sqrt(gamma1**2 + gamma2**2)
    return phi, gamma


@njit(cache=True, fastmath=True)
def _ellipticity2phi_q(e1, e2):
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
def _shear_polar2cartesian(phi, gamma):
    """
    Convert polar shear representation to Cartesian components.

    Parameters
    ----------
    phi : ``float``
        Shear angle in radians.
    gamma : ``float``
        Shear magnitude.

    Returns
    -------
    gamma1 : ``float``
        First Cartesian shear component.
    gamma2 : ``float``
        Second Cartesian shear component.
    """
    gamma1 = gamma * np.cos(2 * phi)
    gamma2 = gamma * np.sin(2 * phi)
    return gamma1, gamma2

@njit(cache=True, fastmath=True)
def pol_to_ell(r, theta, q):
    """
    Convert polar coordinates to elliptical coordinates.

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

    Examples
    --------
    >>> rell, phi = pol_to_ell(1.0, 0.5, 0.8)
    """
    phi = np.arctan2(np.sin(theta), np.cos(theta) * q)
    rell = r * np.sqrt(q**2 * np.cos(theta) ** 2 + np.sin(theta) ** 2)
    return rell, phi

@njit(cache=True, fastmath=True)
def _rotmat(th):
    """
    Compute a 2D rotation matrix.

    Parameters
    ----------
    th : ``float``
        Rotation angle in radians.

    Returns
    -------
    M : ``numpy.ndarray``
        2×2 rotation matrix.
    """
    return np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]])

@njit(cache=True, fastmath=True)
def omega(phi, t, q, omegas, niter_max=200, tol=1e-16):
    """
    Evaluate the complex angular function Omega for the EPL profile.

    This series expansion converges geometrically with ratio
    ``f = (1 - q)/(1 + q)``. The ``fastmath`` flag provides ~4x speedup
    due to the reduction nature of the summation.

    Parameters
    ----------
    phi : ``numpy.ndarray``
        Azimuthal angles in radians.
    t : ``float``
        EPL slope exponent (``t = gamma - 1``).
    q : ``float``
        Axis ratio.
    omegas : ``numpy.ndarray``
        Pre-allocated complex output buffer with the same shape as ``phi``.
        Filled in-place.
    niter_max : ``int``
        Maximum number of series terms. \n
        default: 200
    tol : ``float``
        Convergence tolerance. \n
        default: 1e-16

    Returns
    -------
    omegas : ``numpy.ndarray``
        Complex Omega values at each angle (same object as input buffer).

    Examples
    --------
    >>> import numpy as np
    >>> phi = np.linspace(0, 2 * np.pi, 100)
    >>> result = np.empty_like(phi, dtype=np.complex128)
    >>> omega(phi, t=1.0, q=0.8, omegas=result)
    """
    f = (1 - q) / (1 + q)
    omegas[:] = 0.0 + 0.0j
    niter = min(
        niter_max, int(np.log(tol) / np.log(f)) + 2
    )  # The absolute value of each summand is always less than f, hence this limit for the number of iterations.
    Omega = 1 * np.exp(1j * phi)
    fact = -f * np.exp(2j * phi)
    for n in range(1, niter):
        omegas += Omega
        Omega *= (2 * n - (2 - t)) / (2 * n + (2 - t)) * fact
    omegas += Omega
    return omegas


@njit(cache=True, fastmath=True)
def omega_scalar(phi, t, q, niter_max=200, tol=1e-16):
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

@njit(cache=True, fastmath=True)
def cdot(a, b):
    """
    Compute the real-valued dot product of two complex numbers.

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
def _solvequadeq(a, b, c):
    """
    Solve a quadratic equation with numerically careful root selection.

    Uses sign-stabilized formulas to avoid loss of significance.
    See https://en.wikipedia.org/wiki/Loss_of_significance.

    Parameters
    ----------
    a : ``numpy.ndarray``
        Quadratic coefficient.
    b : ``numpy.ndarray``
        Linear coefficient.
    c : ``numpy.ndarray``
        Constant coefficient.

    Returns
    -------
    x1 : ``numpy.ndarray``
        First root.
    x2 : ``numpy.ndarray``
        Second root.
    """
    sD = (b**2 - 4 * a * c) ** 0.5
    x1 = (-b - np.sign(b) * sD) / (2 * a)
    x2 = 2 * c / (-b - np.sign(b) * sD)
    return np.where(b != 0, np.where(a != 0, x1, -c / b), -((-c / a) ** 0.5)), np.where(
        b != 0, np.where(a != 0, x2, -c / b + 1e-8), +((-c / a) ** 0.5)
    )

@njit(cache=True, fastmath=True)
def pol_to_cart(r, th):
    """
    Convert polar coordinates to Cartesian coordinates.

    Parameters
    ----------
    r : ``float`` or ``numpy.ndarray``
        Radial coordinate.
    th : ``float`` or ``numpy.ndarray``
        Polar angle in radians.

    Returns
    -------
    x : ``float`` or ``numpy.ndarray``
        Cartesian x-coordinate.
    y : ``float`` or ``numpy.ndarray``
        Cartesian y-coordinate.

    Examples
    --------
    >>> x, y = pol_to_cart(1.0, np.pi / 4)
    """
    return r * np.cos(th), r * np.sin(th)

@njit(cache=True)
def _alpha_epl_shear(x, y, b, q, t, gamma1, gamma2, Omega):
    """
    Compute the complex deflection of an EPL + external shear lens.

    Parameters
    ----------
    x : ``float`` or ``numpy.ndarray``
        Image-plane x-coordinate.
    y : ``float`` or ``numpy.ndarray``
        Image-plane y-coordinate.
    b : ``float``
        EPL strength parameter (``theta_E * sqrt(q)``).
    q : ``float``
        Axis ratio.
    t : ``float``
        EPL slope exponent (``t = gamma - 1``).
    gamma1 : ``float``
        External shear component 1.
    gamma2 : ``float``
        External shear component 2.
    Omega : ``numpy.ndarray``
        Precomputed Omega array.

    Returns
    -------
    deflection : ``complex`` or ``numpy.ndarray``
        Complex deflection angle.
    """
    zz = x * q + 1j * y
    R = np.abs(zz)
    alph = (2 * b) / (1 + q) * np.nan_to_num((b / R) ** t * R / b) * Omega

    return (
        alph
        + (gamma1 * x + gamma2 * y)
        + 1j * (gamma2 * x - gamma1 * y)
    )


@njit(cache=True)
def _alpha_epl_shear_scalar(x, y, b, q, t=1, gamma1=0, gamma2=0):
    """
    Scalar complex deflection for EPL + external shear.
    """
    zz = x * q + 1j * y
    R = np.abs(zz)
    phi = np.angle(zz)
    Omega = omega_scalar(phi, t, q)
    alph = (2 * b) / (1 + q) * np.nan_to_num((b / R) ** t * R / b) * Omega
    return (
        alph
        + (gamma1 * x + gamma2 * y)
        + 1j * (gamma2 * x - gamma1 * y)
    )

@njit(cache=True, fastmath=True)
def cart_to_pol(x, y):
    """
    Convert Cartesian coordinates to polar coordinates.

    The returned angle is wrapped to ``[0, 2π)``.

    Parameters
    ----------
    x : ``float`` or ``numpy.ndarray``
        Cartesian x-coordinate.
    y : ``float`` or ``numpy.ndarray``
        Cartesian y-coordinate.

    Returns
    -------
    r : ``float`` or ``numpy.ndarray``
        Radial coordinate.
    theta : ``float`` or ``numpy.ndarray``
        Polar angle in radians, wrapped to ``[0, 2π)``.

    Examples
    --------
    >>> r, theta = cart_to_pol(1.0, 1.0)
    """
    return np.sqrt(x**2 + y**2), np.arctan2(y, x) % (2 * np.pi)

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


@njit(cache=True)
def caustics_epl_shear(
    q, phi, gamma, gamma1, gamma2, theta_E=1.0, num_th=500, maginf=-100.0, sourceplane=True, return_which="double"
):
    """
    Analytically compute the caustics of an EPL + external shear lens model.

    For ``gamma > 2`` the outer critical curve does not exist. In that case
    the routine finds the curve at a finite magnification ``maginf`` instead
    of the true caustic.

    Parameters
    ----------
    q : ``float``
        Lens axis ratio.
    phi : ``float``
        Lens position angle in radians.
    gamma : ``float``
        EPL power-law slope.
    gamma1 : ``float``
        External shear component 1.
    gamma2 : ``float``
        External shear component 2.
    theta_E : ``float``
        Einstein radius. \n
        default: 1.0
    num_th : ``int``
        Number of azimuthal samples for the boundary curve. \n
        default: 500
    maginf : ``float``
        Magnification cutoff for the outer curve. \n
        default: -100.0
    sourceplane : ``bool``
        If True, map critical curves to the source plane. \n
        default: True
    return_which : ``str``
        Which boundary to return. \n
        Options: \n
        - 'double': Double-image region boundary \n
        - 'quad': Quad-image region boundary \n
        - 'caustic': Inner caustic (astroid) only \n
        - 'cut': Outer cut curve only \n
        default: 'double'

    Returns
    -------
    pts : ``numpy.ndarray``
        2×N array of (x, y) boundary coordinates.

    Examples
    --------
    >>> pts = caustics_epl_shear(
    ...     q=0.8, phi=0.0, gamma=2.0, gamma1=0.03, gamma2=-0.01
    ... )
    """

    # Convert lens slope to internal EPL exponent (t = gamma - 1)
    t = gamma - 1.0

    # Convert external shear to polar form
    theta_ell = phi
    gamma1unr, gamma2unr = gamma1, gamma2
    theta_gamma, gamma_mag = _shear_cartesian2polar(gamma1unr, gamma2unr)

    # Build Einstein-radius normalization and center offset
    b = np.sqrt(q) * theta_E
    cen = np.expand_dims(
        np.array([0.0, 0.0]), 1
    )

    # Rotate shear into lens-aligned coordinates
    theta_gamma -= theta_ell
    gamma1, gamma2 = _shear_polar2cartesian(theta_gamma, gamma_mag)
    M = _rotmat(-theta_ell)

    # Sample azimuthal angles and compute elliptical polar quantities
    theta = np.linspace(0, 2 * np.pi * (1.0 - 1.0 / num_th), num_th)
    r = 1
    R, phi = pol_to_ell(1, theta, q)
    Omega = np.empty_like(phi, dtype=np.complex128)
    omega(phi, t, q, Omega)

    # Assemble quadratic coefficients for inverse-radius solutions
    aa = np.ones_like(theta)
    bb = np.full_like(theta, -(2 - t))
    frac_roverR = r / R
    cc = (
        (1 - t)
        * (2 - t)
        * (cdot(np.exp(1j * theta), Omega))
        / frac_roverR
        * 2
        / (1 + q)
    )
    cc -= (1 - t) ** 2 * (2 / (1 + q)) ** 2 * np.abs(Omega) ** 2 / frac_roverR**2

    # Add external shear contribution to the quadratic system
    gammaint_fac = (
        -np.exp(2j * theta) * (2 - t) / 2
        + (1 - t) * np.exp(1j * theta) * 2 / (1 + q) * Omega / frac_roverR
    )
    gamma = gamma1 + 1j * gamma2
    aa -= np.abs(gamma) ** 2
    bb -= 2 * cdot(gamma, gammaint_fac)

    # Solve for the main (quad) critical curve branch
    usol = np.stack(_solvequadeq(cc, bb, aa)).T
    xcr_4, ycr_4 = pol_to_cart(b * usol[:, 1] ** (-1 / t) * frac_roverR, theta)

    # Solve the secondary branch (cut) with t-dependent magnification convention
    if (
        t > 1
    ):  # If t>1, get the approximate outer caustic instead (where inverse magnification = maginf).
        usol = np.stack(_solvequadeq(cc, bb, aa - maginf)).T
        xcr_cut, ycr_cut = pol_to_cart(b * usol[:, 1] ** (-1 / t) * frac_roverR, theta)
    else:
        usol = np.stack(_solvequadeq(cc, bb, aa + maginf)).T
        xcr_cut, ycr_cut = pol_to_cart(b * usol[:, 0] ** (-1 / t) * frac_roverR, theta)

    # Compute deflection and map critical curves to caustics if requested
    al_cut = _alpha_epl_shear(xcr_cut, ycr_cut, b, q, t, gamma1, gamma2, Omega=Omega)
    al_4 = _alpha_epl_shear(xcr_4, ycr_4, b, q, t, gamma1, gamma2, Omega=Omega)
    if sourceplane:
        xca_cut, yca_cut = xcr_cut - al_cut.real, ycr_cut - al_cut.imag
        xca_4, yca_4 = xcr_4 - al_4.real, ycr_4 - al_4.imag
    else:
        xca_cut, yca_cut = xcr_cut, ycr_cut
        xca_4, yca_4 = xcr_4, ycr_4

    # Return individual branches early when explicitly requested
    if return_which == "caustic":
        return M @ np.stack((xca_4, yca_4)) + cen
    if return_which == "cut":
        return M @ np.stack((xca_cut, yca_cut)) + cen

    # Interpolate cut radius onto caustic angles for boundary composition
    rcut, thcut = cart_to_pol(xca_cut, yca_cut)
    r, th = cart_to_pol(xca_4, yca_4)
    # Sort cut data by angle for interpolation
    sort_idx = np.argsort(thcut)
    thcut_sorted = thcut[sort_idx]
    rcut_sorted = rcut[sort_idx]
    r2 = _interp_periodic(th, thcut_sorted, rcut_sorted, 2.0 * np.pi)

    # Build either the double-image or quad-image sampling boundary
    if return_which == "double":
        r = np.fmax(r, r2)
    else:  # Quad
        r = np.fmin(r, r2)

    # Convert the selected radial boundary back to Cartesian samples
    pos_tosample = np.empty((2, num_th))
    pos_tosample[0], pos_tosample[1] = pol_to_cart(r, th)

    # Return requested region
    return M @ pos_tosample + cen

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


@njit(cache=True, fastmath=True)
def caustic_double_area(
    q, phi, gamma, gamma1, gamma2, theta_E=1.0, num_th=500, maginf=-100.0, sourceplane=True, return_which="double"
):
    """
    Compute the area enclosed by the double caustic of an EPL + shear lens.

    Parameters
    ----------
    q : ``float``
        Lens axis ratio.
    phi : ``float``
        Lens position angle in radians.
    gamma : ``float``
        Power-law slope.
    gamma1 : ``float``
        External shear component 1.
    gamma2 : ``float``
        External shear component 2.
    theta_E : ``float``
        Einstein radius. 
        default: 1.0
    num_th : ``int``
        Number of angular samples. 
        default: 500
    maginf : ``float``
        Magnification cutoff parameter. 
        default: -100.0

    Returns
    -------
    area : ``float``
        Area of the double caustic in units of theta_E^2.

    Examples
    --------
    >>> area = caustic_double_area(0.8, 0.0, 2.0, 0.03, -0.01)
    """
    pts = caustics_epl_shear(
        q, phi, gamma, gamma1, gamma2, theta_E=theta_E, num_th=num_th, maginf=maginf, sourceplane=sourceplane, return_which=return_which
    )
    if np.any(~np.isfinite(pts)):
        return 0.0
    area = polygon_area(pts[0], pts[1])
    if not np.isfinite(area):
        return 0.0
    return area

def make_cross_section_reinit(
    Da_instance,
    num_th=500, 
    maginf=-100.0, 
    sourceplane=True, 
    return_which="double",
):
    """
    Create a JIT-compiled cross-section evaluator for batched systems.

    Parameters
    ----------
    Da_instance : ``callable``
        Angular-diameter-distance function that accepts redshift arrays.
    num_th : ``int``
        Number of angular samples for caustic construction. 
        default: 500
    maginf : ``float``
        Magnification cutoff used by the caustic routine. 
        default: -100.0

    Returns
    -------
    cross_section_reinit : ``callable``
        Parallel numba function with signature 
        ``(zs, zl, sigma, q, phi, gamma, gamma1, gamma2)`` returning 
        the double-caustic cross section array.

    Examples
    --------
    >>> cross_section = make_cross_section_reinit(Da_instance)
    >>> cs = cross_section(zs, zl, sigma, q, phi, gamma, gamma1, gamma2)
    """

    @njit(parallel=True, cache=True)
    def cross_section_reinit(zs, zl, sigma, q, phi, gamma, gamma1, gamma2):
        """
        Compute double-caustic cross sections for batched lens parameters.

        Parameters
        ----------
        zs : ``numpy.ndarray``
            Source redshifts.
        zl : ``numpy.ndarray``
            Lens redshifts.
        sigma : ``numpy.ndarray``
            Velocity dispersions in km/s.
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
        cs : ``numpy.ndarray``
            Cross section values in angular units.
        """
        size = zs.size
        Ds = Da_instance(zs)
        Dl = Da_instance(zl)
        Dls = (Ds * (1 + zs) - Dl * (1 + zl)) / (1 + zs)
        theta_E = 4.0 * np.pi * (sigma * 1000.0 / C_LIGHT) ** 2 * (Dls / Ds)

        cs = np.zeros(size)
        for i in prange(size):

            area = caustic_double_area(
                q[i],
                phi[i],
                gamma[i],
                gamma1[i],
                gamma2[i],
                theta_E=theta_E[i],
                num_th=num_th, 
                maginf=maginf, 
                sourceplane=sourceplane, 
                return_which=return_which,
            )
            if np.isfinite(area):
                cs[i] = area
            else:
                cs[i] = 0.0
        return cs

    return cross_section_reinit


@njit(parallel=True, cache=True)
def cross_section_epl_shear_unit(
    e1, 
    e2, 
    gamma, 
    gamma1, 
    gamma2,
    theta_E=None,
    num_th=500,
    maginf=-100.0,
    sourceplane=True,
    return_which="double",
):
    """
    Compute double-caustic cross sections for batched lens parameters.

    Parameters
    ----------
    e1 : ``numpy.ndarray``
        Lens ellipticity component 1.
    e2 : ``numpy.ndarray``
        Lens ellipticity component 2.
    gamma : ``numpy.ndarray``
        EPL slopes.
    gamma1 : ``numpy.ndarray``
        External shear component 1.
    gamma2 : ``numpy.ndarray``
        External shear component 2.
    theta_E : ``numpy.ndarray`` or ``None``
        Einstein radius in angular units. If None, assumed to be 1.0 for all
    Returns
    -------
    cs : ``numpy.ndarray``
        Cross section values in angular units.
    """
    size = e1.size

    if theta_E is None:
        theta_E_arr = np.ones(size, dtype=np.float64)
    else:
        theta_E_arr = theta_E

    # convert ellipticity to axis ratio and position angle
    phi, q = _ellipticity2phi_q(e1, e2)

    cs = np.zeros(size)
    for i in prange(size):

        cs[i] = caustic_double_area(
            q=q[i],
            phi=phi[i],
            gamma=gamma[i],
            gamma1=gamma1[i],
            gamma2=gamma2[i],
            theta_E=theta_E_arr[i],
            num_th=num_th, 
            maginf=maginf, 
            sourceplane=sourceplane, 
            return_which=return_which,
        )
        # if np.isfinite(area):
        #     cs[i] = area
        # else:
        #     cs[i] = 0.0
    cs[~np.isfinite(cs)] = 0.0
    return np.maximum(cs, 0.0)
