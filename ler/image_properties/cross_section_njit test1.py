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
def ellipticity2phi_q(e1, e2):
    """
    Convert complex ellipticity moduli to orientation angle and axis ratio.
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
    """
    gamma1 = gamma * np.cos(2 * phi)
    gamma2 = gamma * np.sin(2 * phi)
    return gamma1, gamma2

@njit(cache=True, fastmath=True)
def pol_to_ell(r, theta, q):
    """
    Convert polar coordinates to elliptical coordinates.
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
def omega(phi, t, q, niter_max=200, tol=1e-16):
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
    omegas = np.zeros_like(phi, dtype=np.complex128)
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
    # # Clamp discriminant: angles with no real critical curve get a repeated root
    # # instead of NaN from sqrt(negative).
    # disc = np.maximum(b**2 - 4.0 * a * c, 0.0)
    # sD = np.sqrt(disc)
    # denominator = -b - np.sign(b) * sD
    # x1 = np.where(np.abs(denominator) > 1e-30, denominator / (2 * a), -sD / (2.0 * a))
    # x2 = np.where(np.abs(denominator) > 1e-30, 2 * c / denominator, sD / (2.0 * a))
    # return x1, x2

    sD = (b**2 - 4 * a * c) ** 0.5
    x1 = (-b - np.sign(b) * sD) / (2 * a)
    x2 = 2 * c / (-b - np.sign(b) * sD)

    return (
        np.where(
            b != 0, 
            np.where(
                a != 0, 
                x1, 
                -c / b
            ), 
            -((-c / a) ** 0.5)
        ), 
        np.where(
            b != 0, 
            np.where(
                a != 0, 
                x2, -c / b + 1e-8
            ), 
            +((-c / a) ** 0.5)
        )
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

    alpha_EPL = (2 * b) / (1 + q) * (b / R)^t * R / b * Omega
    alpha_Shear = gamma1 * x + gamma2 * y + i * (gamma2 * x - gamma1 * y)
    """
    zz = x * q + 1j * y
    R = np.abs(zz)
    #alph = (2 * b) / (1 + q) * np.nan_to_num((b / R) ** t * R / b) * Omega
    alph = (2* b**t)/(1 + q) * R**(1-t) * Omega

    return (
        alph
        + (gamma1 * x + gamma2 * y)
        + 1j * (gamma2 * x - gamma1 * y)
    )

@njit(cache=True)
def _alpha_epl_shear_scalar(x, y, b, q, t=1, gamma1=0, gamma2=0):
    """
    Scalar complex deflection for EPL + external shear. Same as _alpha_epl_shear but avoids temporary array allocation by calling the scalar omega function.

    alpha_EPL = (2 * b) / (1 + q) * (b / R)^t * R / b * Omega
    alpha_Shear = gamma1 * x + gamma2 * y + i * (gamma2 * x - gamma1 * y)
    """
    zz = x * q + 1j * y
    R = np.abs(zz)
    phi = np.angle(zz)
    Omega = omega_scalar(phi, t, q)
    # alph = (2 * b) / (1 + q) * np.nan_to_num((b / R) ** t * R / b) * Omega
    alph = (2* b**t)/(1 + q) * R**(1-t) * Omega
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
        lo, hi = 0, n - 1
        while lo < hi:
            mid = (lo + hi) // 2
            if xp[mid] < xi:
                lo = mid + 1
            else:
                hi = mid
        # lo is the index of the first element >= xi
        if lo == 0:
            # xi is before the first point – wrap around
            x0, x1 = xp[n - 1] - period, xp[0]
            f0, f1 = fp[n - 1], fp[0]
        else:
            x0, x1 = xp[lo - 1], xp[lo]
            f0, f1 = fp[lo - 1], fp[lo]

        dx = x1 - x0
        if abs(dx) < 1e-30:
            result[i] = f0
        else:
            result[i] = f0 + (f1 - f0) * (xi - x0) / dx
    return result

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
def _helper_caustic_epl_shear(
    q, phi, gamma, gamma1, gamma2, theta_E,
    theta, cos_th, sin_th, cos_2th, sin_2th,
    maginf=-100.0
):
    """
    Generates the unrotated radial boundaries of the EPL+shear caustic 
    and writes them into the provided workspace arrays.
    """

    t = gamma - 1.0
    b = np.sqrt(q) * theta_E

    # rotation M is done outside this helper to avoid redundant computation

    # convert shear to polar form and rotate into lens-aligned coordinates
    theta_gamma = np.arctan2(gamma2, gamma1) / 2
    gamma_mag = np.sqrt(gamma1**2 + gamma2**2)
    theta_gamma -= phi # rotate shear into lens-aligned coordinates
    gamma1_rot = gamma_mag * np.cos(2 * theta_gamma)
    gamma2_rot = gamma_mag * np.sin(2 * theta_gamma)

    # R=np.sqrt(q**2 * x**2 + y**2), with x=cos(theta), y=sin(theta)
    # (x,y) are the coordinates of the unit circle, so this is just the elliptical radius at each angle.
    # theta is the polar angle of the unit circle, so phi_theta_ell is the corresponding angle in the elliptical coordinate system.
    # theta = np.linspace(0, 2 * np.pi * (1.0 - 1.0 / num_th), num_th)
    r = 1.0
    R_theta_ell = r * np.sqrt(q**2 * cos_th
     ** 2 + sin_th ** 2)
    phi_theta_ell = np.arctan2(sin_th, cos_th * q)
    # angular part of the deflection
    # deflection = fn1(R) x fn2(theta), where fn1 is the radial scaling and fn2 is the angular function Omega.
    # no output, just populate the given Omega array
    Omega = omega(phi_theta_ell, t, q) # check with the scalar version also

    # Assemble quadratic coefficients for inverse-radius solutions
    # (K^2 - |\vec{\Gamma}|^2) u^2 - (2K + 2\vec{\Gamma} \cdot \vec{\gamma}_{ext}) u + (1 - |\vec{\gamma}_{ext}|^2) = 0
    # \vec{\Gamma} = gammaint_fac , 
    # aa = 1 - |gamma_ext|^2
    # bb = -2(2-t)/2 - 2 * cdot(gamma, gammaint_fac)
    # cc = K**2 - gammaint_fac**2 , K = (2-t)/2
    # cc = 2*K*A*Re(exp(i*theta)*Omega) - A^2*|Omega|^2, A = (1-t)*2/(1+q)/R

    # constant term 
    aa = np.ones_like(theta)
    gamma_ext = gamma1_rot + 1j * gamma2_rot # external shear, np.abs(gamma_ext)**2 = gamma1_rot**2 + gamma2_rot**2
    aa -= np.abs(gamma_ext) ** 2
    
    # linear term
    bb = np.full_like(theta, -(2 - t))
    frac_roverR = r / R_theta_ell
    # internal shear contribution
    gammaint_fac = (
        -np.exp(2j * theta) * (2 - t) / 2
        + (1 - t) * np.exp(1j * theta) * 2 / (1 + q) * Omega / frac_roverR
    ) # gamma_internal = gammaint_fac * u
    # external shear contribution to the linear term of the quadratic system
    bb -= 2 * cdot(gamma_ext, gammaint_fac)

    # quadratic term
    cc = (
        (1 - t)
        * (2 - t)
        * (cdot(np.exp(1j * theta), Omega))
        / frac_roverR
        * 2
        / (1 + q)
    )
    cc -= (1 - t) ** 2 * (2 / (1 + q)) ** 2 * np.abs(Omega) ** 2 / frac_roverR**2

    # Solve for the main (quad) critical curve branch
    # solve for u = (R/r)^t
    # u = (R/r)^t is physically a positive radial ratio; negative roots are
    # unphysical (occur when cc < 0, aa > 0). np.abs preserves the correct
    # magnitude while keeping the result real.
    usol = np.stack(_solvequadeq(cc, bb, aa)).T
    # if usol[:, 1] < 1e-30:
    #     u4 = usol[:, 0]
    # else:
    #     u4 = usol[:, 1]
    u4 = np.where(usol[:, 1] > 0., usol[:, 1], usol[:, 0])
    xcr_4, ycr_4 = pol_to_cart(
        b * u4 ** (-1 / t) * frac_roverR,
        theta
    )

    # Solve the secondary branch (cut) with t-dependent magnification convention
    if (
        t > 1
    ):  # If t>1, get the approximate outer caustic instead (where inverse magnification = maginf).
        usol = np.stack(_solvequadeq(cc, bb, aa - maginf)).T
        u_cut = np.where(usol[:, 1] > 0., usol[:, 1], usol[:, 0])
        xcr_cut, ycr_cut = pol_to_cart(
            b * u_cut ** (-1 / t) * frac_roverR,
            theta
        )
    else:
        usol = np.stack(_solvequadeq(cc, bb, aa + maginf)).T
        u_cut = np.where(usol[:, 1] > 0., usol[:, 1], usol[:, 0])
        xcr_cut, ycr_cut = pol_to_cart(
            b * u_cut ** (-1 / t) * frac_roverR,
            theta
        )

    # Compute deflection and map critical curves to caustics 
    al_cut = _alpha_epl_shear(xcr_cut, ycr_cut, b, q, t, gamma1_rot, gamma2_rot, Omega=Omega)
    al_4 = _alpha_epl_shear(xcr_4, ycr_4, b, q, t, gamma1_rot, gamma2_rot, Omega=Omega)

    xca_cut, yca_cut = xcr_cut - al_cut.real, ycr_cut - al_cut.imag
    xca_4, yca_4 = xcr_4 - al_4.real, ycr_4 - al_4.imag
    # If not source plane, the caustic points are just the critical curve points, so no deflection mapping.
    # xca_cut, yca_cut = xcr_cut, ycr_cut
    # xca_4, yca_4 = xcr_4, ycr_4

    # # Return individual branches early when explicitly requested
    # if return_which == "caustic":
    #     return M @ np.stack((xca_4, yca_4)) + cen
    # if return_which == "cut":
    #     return M @ np.stack((xca_cut, yca_cut)) + cen

    # Interpolate cut radius onto caustic angles for boundary composition
    rcut, thcut = cart_to_pol(xca_cut, yca_cut)
    r, th = cart_to_pol(xca_4, yca_4)
    # Sort cut data by angle for interpolation
    sort_idx = np.argsort(thcut)
    thcut_sorted = thcut[sort_idx]
    rcut_sorted = rcut[sort_idx]
    r2 = _interp_periodic(th, thcut_sorted, rcut_sorted, 2.0 * np.pi)

    # Build either the double-image or quad-image sampling boundary
    # if return_which == "double":
    r = np.fmax(r, r2)
    # else:  # Quad
    #     r = np.fmin(r, r2)

    # Convert the selected radial boundary back to Cartesian samples
    num_th = theta.shape[0]
    pos_tosample = np.empty((2, num_th))
    pos_tosample[0], pos_tosample[1] = pol_to_cart(r, th)
    return pos_tosample

    # return pol_to_cart(r, th)

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
    
    # Call the helper
    pos_tosample = _helper_caustic_epl_shear(
        q, phi, gamma, gamma1, gamma2, theta_E,
        theta, cos_th, sin_th, cos_2th, sin_2th,
        maginf=maginf
    )

    M = _rotmat(-phi)
    return M @ pos_tosample
    # # Numba-friendly explicit 2D rotation (equivalent to _rotmat(-phi) @ pos_tosample)
    # c = np.cos(phi)
    # s = np.sin(phi)
    # rotated = np.empty_like(pos_tosample)
    # for i in range(pos_tosample.shape[1]):
    #     x = pos_tosample[0, i]
    #     y = pos_tosample[1, i]
    #     rotated[0, i] = c * x - s * y
    #     rotated[1, i] = s * x + c * y
    # return rotated

# 2. THE AREA CALCULATOR
@njit(cache=True, fastmath=True)
def caustic_area_epl_shear(
    q, phi, gamma, gamma1, gamma2, theta_E,
    theta, cos_th, sin_th, cos_2th, sin_2th,
    maginf=-100.0
):
    """
    Calculates the 2D coordinates of the double caustic for a SINGLE lens.
    Accepts scalar float values.
    """
    
    # Call the helper
    # pos_tosample = np.empty((2, num_th))
    pos_tosample = _helper_caustic_epl_shear(
        q, phi, gamma, gamma1, gamma2, theta_E,
        theta, cos_th, sin_th, cos_2th, sin_2th,
        maginf=maginf
    )

    # M = _rotmat(-phi) # no rotation needed for area calculation
    # return area
    
    return polygon_area(pos_tosample[0, :], pos_tosample[1, :])

@njit(cache=True, fastmath=True)
def polygon_area(xv, yv):
    """
    Compute the area of a simple polygon using the Shoelace formula.
    """
    area = 0.0
    n = xv.shape[0]
    j = n - 1
    for i in range(n):
        area += xv[j] * yv[i] - xv[i] * yv[j]
        j = i
    return abs(area) / 2.0

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
        size = zs.shape[0]
        cs_area = np.empty(size)

        Ds = Da_instance(zs)
        Dl = Da_instance(zl)

        for i in prange(size):
            Dls = (Ds[i] * (1.0 + zs[i]) - Dl[i] * (1.0 + zl[i])) / (1.0 + zs[i])
            theta_E = 4.0 * PI * (sigma[i] * 1000.0 / C_LIGHT) ** 2 * (Dls / Ds[i])

            area = caustic_area_epl_shear(
                q[i], phi[i], gamma[i], gamma1[i], gamma2[i], theta_E,
                theta, cos_th, sin_th, cos_2th, sin_2th,
                maginf
            )

            cs_area[i] = area if np.isfinite(area) else 0.0

        return cs_area

    return cross_section_area

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

    @njit(parallel=True, cache=True, fastmath=True)
    def cross_section_area(zs, zl, sigma, q, phi, gamma, gamma1, gamma2):
        size = zs.shape[0]
        cs_area = np.empty(size)

        for i in prange(size):
            theta_E = 1.0

            area = caustic_area_epl_shear(
                q[i], phi[i], gamma[i], gamma1[i], gamma2[i], theta_E,
                theta, cos_th, sin_th, cos_2th, sin_2th,
                maginf
            )

            cs_area[i] = area if np.isfinite(area) else 0.0

        return cs_area

    return cross_section_area
