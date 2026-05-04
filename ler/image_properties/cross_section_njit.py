"""
Module for analytical caustic computation and cross-section evaluation.

Provides numba-accelerated routines for computing critical curves, caustic
boundaries, polygon areas, and lensing cross sections for EPL (Elliptical
Power-Law) + external shear lens models. All core functions are decorated
with ``@njit`` for high performance.

Usage:
    Basic workflow example:

    >>> from ler.image_properties.cross_section_njit import caustic_points_epl_shear
    >>> pts = caustic_points_epl_shear(theta_E=1.0, q=0.8, phi=0.0, gamma=2.0, gamma1=0.03, gamma2=-0.01)

Copyright (C) 2026 Phurailatpam Hemantakumar. Distributed under MIT License.
"""

import numpy as np
from numba import njit, prange

C_LIGHT = 299792458.0  # m/s
PI = np.pi
TWO_PI = 2.0 * PI
EPS = 1e-16

@njit(cache=True, fastmath=True)
def phi_q2_ellipticity(phi, q):
    """
    Convert lens orientation and axis ratio to ellipticity components.

    Parameters
    ----------
    phi : ``float``
        Position angle of the lens major axis in radians. \n
    q : ``float``
        Axis ratio (minor/major), where ``0 < q <= 1``. \n

    Returns
    -------
    e1 : ``float``
        First ellipticity component. \n
    e2 : ``float``
        Second ellipticity component. \n

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
        First Cartesian shear component. \n
    gamma2 : ``float``
        Second Cartesian shear component. \n

    Returns
    -------
    phi : ``float``
        Shear position angle in radians. \n
    gamma : ``float``
        Shear magnitude. \n
    """
    phi = np.arctan2(gamma2, gamma1) / 2
    gamma = np.sqrt(gamma1**2 + gamma2**2)
    return phi, gamma


@njit(cache=True, fastmath=True)
def ellipticity2phi_q(e1, e2):
    """
    Convert complex ellipticity moduli to orientation angle and axis ratio.

    Parameters
    ----------
    e1 : ``float``
        First ellipticity component. \n
    e2 : ``float``
        Second ellipticity component. \n

    Returns
    -------
    phi : ``float``
        Orientation angle in radians. \n
    q : ``float``
        Axis ratio (minor/major). \n

    Examples
    --------
    >>> phi, q = ellipticity2phi_q(0.1, 0.05)
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
        Shear position angle in radians. \n
    gamma : ``float``
        Shear magnitude. \n

    Returns
    -------
    gamma1 : ``float``
        First Cartesian shear component. \n
    gamma2 : ``float``
        Second Cartesian shear component. \n
    """
    gamma1 = gamma * np.cos(2 * phi)
    gamma2 = gamma * np.sin(2 * phi)
    return gamma1, gamma2

@njit(cache=True, fastmath=True)
def _rotmat(th):
    """
    Compute a 2D rotation matrix.

    Parameters
    ----------
    th : ``float``
        Rotation angle in radians. \n

    Returns
    -------
    M : ``numpy.ndarray``
        2×2 rotation matrix. \n
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
        Azimuthal angles in radians. \n
    t : ``float``
        EPL slope exponent (``t = gamma - 1``). \n
    q : ``float``
        Axis ratio. \n
    niter_max : ``int``
        Maximum number of series terms. \n
        default: 200
    tol : ``float``
        Convergence tolerance. \n
        default: 1e-16

    Returns
    -------
    omegas : ``numpy.ndarray``
        Complex Omega values at each angle. \n

    Examples
    --------
    >>> import numpy as np
    >>> phi = np.linspace(0, 2 * np.pi, 100)
    >>> omegas = omega(phi, t=1.0, q=0.8)
    """
    f = (1 - q) / (1 + q)
    if f <= EPS:
        return np.exp(1j * phi)
        
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
def cdot(a, b):
    """
    Compute the real-valued dot product of two complex numbers.

    Equivalent to ``Re(a) * Re(b) + Im(a) * Im(b)``.

    Parameters
    ----------
    a : ``complex``
        First complex number. \n
    b : ``complex``
        Second complex number. \n

    Returns
    -------
    result : ``float``
        Real-valued dot product. \n

    Examples
    --------
    >>> cdot(1+2j, 3+4j)
    11.0
    """
    return a.real * b.real + a.imag * b.imag

@njit(cache=True, fastmath=True)
def _solvequadeq(a, b, c):
    """
    Solve a quadratic equation with numerically careful root selection.

    Uses sign-stabilized formulas to avoid loss of significance.
    See https://en.wikipedia.org/wiki/Loss_of_significance.

    Parameters
    ----------
    a : ``numpy.ndarray``
        Quadratic coefficient. \n
    b : ``numpy.ndarray``
        Linear coefficient. \n
    c : ``numpy.ndarray``
        Constant coefficient. \n

    Returns
    -------
    x1 : ``numpy.ndarray``
        First root. \n
    x2 : ``numpy.ndarray``
        Second root. \n
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
        Radial coordinate. \n
    th : ``float`` or ``numpy.ndarray``
        Polar angle in radians. \n

    Returns
    -------
    x : ``float`` or ``numpy.ndarray``
        Cartesian x-coordinate. \n
    y : ``float`` or ``numpy.ndarray``
        Cartesian y-coordinate. \n

    Examples
    --------
    >>> x, y = pol_to_cart(1.0, np.pi / 4)
    """
    return r * np.cos(th), r * np.sin(th)

@njit(cache=True, fastmath=True)
def _alpha_epl_shear(x, y, b, q, t, gamma1, gamma2, Omega):
    """
    Compute the complex deflection of an EPL + external shear lens.

    alpha_EPL = (2 * b) / (1 + q) * (b / R)^t * R / b * Omega
    alpha_Shear = gamma1 * x + gamma2 * y + i * (gamma2 * x - gamma1 * y)

    Parameters
    ----------
    x : ``numpy.ndarray``
        Cartesian x-coordinate. \n
    y : ``numpy.ndarray``
        Cartesian y-coordinate. \n
    b : ``float``
        EPL scale radius. \n
    q : ``float``
        Axis ratio. \n
    t : ``float``
        EPL slope exponent (``t = gamma - 1``). \n
    gamma1 : ``float``
        First shear component. \n
    gamma2 : ``float``
        Second shear component. \n
    Omega : ``numpy.ndarray``
        Complex angular function Omega. \n

    Returns
    -------
    alpha : ``numpy.ndarray``
        Complex deflection array. \n
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

@njit(cache=True, fastmath=True)
def cart_to_pol(x, y):
    """
    Convert Cartesian coordinates to polar coordinates.

    The returned angle is wrapped to ``[0, 2π)``.

    Parameters
    ----------
    x : ``float`` or ``numpy.ndarray``
        Cartesian x-coordinate. \n
    y : ``float`` or ``numpy.ndarray``
        Cartesian y-coordinate. \n

    Returns
    -------
    r : ``float`` or ``numpy.ndarray``
        Radial coordinate. \n
    theta : ``float`` or ``numpy.ndarray``
        Polar angle in radians, wrapped to ``[0, 2π)``. \n

    Examples
    --------
    >>> r, theta = cart_to_pol(1.0, 1.0)
    """
    return np.sqrt(x**2 + y**2), np.arctan2(y, x) % (2 * np.pi)

@njit(cache=True, fastmath=True)
def _interp_periodic(x, xp, fp, period):
    """
    Perform numba-compatible periodic linear interpolation.

    Equivalent to ``np.interp`` with the ``period`` keyword.

    Parameters
    ----------
    x : ``numpy.ndarray``
        Query points. \n
    xp : ``numpy.ndarray``
        Known x-values, sorted and within ``[0, period)``. \n
    fp : ``numpy.ndarray``
        Known function values corresponding to ``xp``. \n
    period : ``float``
        Period of the function. \n

    Returns
    -------
    result : ``numpy.ndarray``
        Interpolated values at query points. \n
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
        Rotation angle in radians. \n

    Returns
    -------
    M : ``numpy.ndarray``
        2×2 rotation matrix. \n
    """
    return np.array([[np.cos(th), np.sin(th)], [-np.sin(th), np.cos(th)]])


@njit(cache=True, fastmath=True)
def _helper_caustic_epl_shear(
    q, phi, gamma, gamma1, gamma2, theta_E,
    theta, cos_th, sin_th, cos_2th, sin_2th,
    maginf=-100.0, quad=False, half_only=False
):
    """
    Generates the unrotated radial boundaries of the EPL+shear caustic 
    and writes them into the provided workspace arrays.

    Parameters
    ----------
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
    theta_E : ``float``
        Einstein radius. \n
    theta : ``numpy.ndarray``
        Array of angles. \n
    cos_th : ``numpy.ndarray``
        Array of cosine of angles. \n
    sin_th : ``numpy.ndarray``
        Array of sine of angles. \n
    cos_2th : ``numpy.ndarray``
        Array of cosine of 2*angles. \n
    sin_2th : ``numpy.ndarray``
        Array of sine of 2*angles. \n
    maginf : ``float``
        Magnification cut threshold. \n
        default: -100.0

    Returns
    -------
    pos_tosample : ``numpy.ndarray``
        Cartesian coordinates of the radial boundaries. \n
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
    # replace np.exp(2j * theta) with cos_2th + 1j*sin_2th and np.exp(1j * theta) with cos_th + 1j*sin_th 
    gammaint_fac = (
        -(cos_2th + 1j * sin_2th) * (2 - t) / 2
        + (1 - t) * (cos_th + 1j * sin_th) * 2 / (1 + q) * Omega / frac_roverR
    ) # gamma_internal = gammaint_fac * u
    # external shear contribution to the linear term of the quadratic system
    bb -= 2 * cdot(gamma_ext, gammaint_fac)

    # quadratic term
    cc = (
        (1 - t)
        * (2 - t)
        * (cdot(cos_th + 1j * sin_th, Omega))
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
    x1_4, x2_4 = _solvequadeq(cc, bb, aa)
    # if usol[:, 1] < 1e-30:
    #     u4 = usol[:, 0]
    # else:
    #     u4 = usol[:, 1]
    u4 = np.where(x2_4 > 0., x2_4, x1_4)
    xcr_4, ycr_4 = pol_to_cart(
        b * u4 ** (-1 / t) * frac_roverR,
        theta
    )

    # Solve the secondary branch (cut) with t-dependent magnification convention
    if (
        t > 1
    ):  # If t>1, get the approximate outer caustic instead (where inverse magnification = maginf).
        x1_cut, x2_cut = _solvequadeq(cc, bb, aa - maginf)
        u_cut = np.where(x2_cut > 0., x2_cut, x1_cut)
        xcr_cut, ycr_cut = pol_to_cart(
            b * u_cut ** (-1 / t) * frac_roverR,
            theta
        )
    else:
        x1_cut, x2_cut = _solvequadeq(cc, bb, aa + maginf)
        u_cut = np.where(x1_cut > 0., x1_cut, x2_cut)
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
    # Sort cut data by angle for interpolation.
    # When half_only=True the input theta covers only [0, π).  Both caustic
    # branches have inversion symmetry: rcut[i+h] == rcut[i] and
    # thcut[i+h] == (thcut[i] + π) % 2π.  Mirror the half cut-branch to
    # reconstruct full-circle coverage before the periodic interpolation.
    if half_only:
        nh = thcut.shape[0]
        # Compute mirrored angles via Cartesian negation: arctan2(-y,-x) is
        # more numerically stable than (arctan2(y,x) + π) % 2π, especially
        # near 0 and π where modular arithmetic can lose a ULP.
        thcut_mir = np.empty(2 * nh)
        rcut_mir = np.empty(2 * nh)
        for i in range(nh):
            th_i = thcut[i]
            thcut_mir[i] = th_i
            th_shift = th_i + PI
            if th_shift >= TWO_PI:
                th_shift -= TWO_PI
            thcut_mir[i + nh] = th_shift
            rcut_mir[i] = rcut[i]
            rcut_mir[i + nh] = rcut[i]
        sort_idx = np.argsort(thcut_mir)
        thcut_sorted = thcut_mir[sort_idx]
        rcut_sorted = rcut_mir[sort_idx]
    else:
        sort_idx = np.argsort(thcut)
        thcut_sorted = thcut[sort_idx]
        rcut_sorted = rcut[sort_idx]
    r2 = _interp_periodic(th, thcut_sorted, rcut_sorted, 2.0 * np.pi)

    # Build either the double-image or quad-image sampling boundary
    if quad:
        r = np.fmin(r, r2)
    else:
        r = np.fmax(r, r2)

    # Convert the selected radial boundary back to Cartesian samples
    num_th = theta.shape[0]
    pos_tosample = np.empty((2, num_th))
    pos_tosample[0], pos_tosample[1] = pol_to_cart(r, th)
    return pos_tosample

    # return pol_to_cart(r, th)

# 3. THE POINT CALCULATOR
@njit(cache=True, fastmath=True)
def caustic_points_epl_shear(
    theta_E, q, phi, gamma, gamma1, gamma2, num_th=500, maginf=-100.0, quad=False
):
    """
    Calculates the 2D coordinates of the caustic for a SINGLE lens.
    Accepts scalar float values.

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
    rotated : ``numpy.ndarray``
        Shape ``(2, num_th)`` Cartesian coordinates of the caustic. \n

    Examples
    --------
    >>> pts = caustic_points_epl_shear(theta_E=1.0, q=0.8, phi=0.0, gamma=2.0, gamma1=0.03, gamma2=-0.01)
    >>> pts_quad = caustic_points_epl_shear(theta_E=1.0, q=0.8, phi=0.0, gamma=2.0, gamma1=0.03, gamma2=-0.01, quad=True)
    """
    # Exploit exact 2-fold symmetry C(θ+π) = -C(θ): only compute the first
    # half of the boundary (theta in [0, π)) and mirror to reconstruct the
    # second half.  half_only=True makes _helper extend the cut-branch angles
    # to [0, 2π) for correct _interp_periodic coverage even though the input
    # only covers [0, π).  The full (2, num_th) output has exactly the same
    # values as calling the helper with all num_th angles.
    h = num_th // 2
    theta_h = np.arange(h) * (2.0 * PI / num_th)

    cos_th, sin_th = np.cos(theta_h), np.sin(theta_h)
    cos_2th, sin_2th = np.cos(2.0 * theta_h), np.sin(2.0 * theta_h)

    pos_half = _helper_caustic_epl_shear(
        q, phi, gamma, gamma1, gamma2, theta_E,
        theta_h, cos_th, sin_th, cos_2th, sin_2th,
        maginf=maginf, quad=quad, half_only=True
    )

    # Reconstruct full (2, num_th) output: second half is the antipodal of the first
    pos_full = np.empty((2, num_th))
    pos_full[:, :h] = pos_half
    pos_full[:, h:] = -pos_half

    M = _rotmat(-phi)
    return M @ pos_full
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
    Calculates the area of the double caustic for a SINGLE lens.
    Accepts scalar float values.

    Parameters
    ----------
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
    theta_E : ``float``
        Einstein radius. \n
    theta : ``numpy.ndarray``
        Array of angles. \n
    cos_th : ``numpy.ndarray``
        Array of cosine of angles. \n
    sin_th : ``numpy.ndarray``
        Array of sine of angles. \n
    cos_2th : ``numpy.ndarray``
        Array of cosine of 2*angles. \n
    sin_2th : ``numpy.ndarray``
        Array of sine of 2*angles. \n
    maginf : ``float``
        Magnification cut threshold. \n
        default: -100.0

    Returns
    -------
    area : ``float``
        Area of the caustic. \n

    Examples
    --------
    >>> import numpy as np
    >>> num_th = 500
    >>> theta = np.linspace(0, 2*np.pi*(num_th//2 - 1)/num_th, num_th//2)
    >>> area = caustic_area_epl_shear(0.8, 0.0, 2.0, 0.03, -0.01, 1.0, theta, np.cos(theta), np.sin(theta), np.cos(2*theta), np.sin(2*theta))
    """
    
    # theta, cos_th, sin_th, cos_2th, sin_2th must cover [0, π) with num_th//2
    # points — i.e. the first half of the full linspace.  The caller is
    # responsible for pre-building half-size arrays so that no slicing or
    # extra allocation is needed here.  Close the half-polygon via the
    # antipodal of the start point (exact by 2-fold symmetry) and multiply
    # the area by 2.
    h = theta.shape[0]
    pos_half = _helper_caustic_epl_shear(
        q, phi, gamma, gamma1, gamma2, theta_E,
        theta, cos_th, sin_th, cos_2th, sin_2th,
        maginf=maginf, half_only=True
    )
    return half_symmetric_polygon_area(pos_half[0], pos_half[1])


@njit(cache=True, fastmath=True)
def half_symmetric_polygon_area(xh, yh):
    """
    Area of a centrally symmetric closed polygon from half-boundary samples.

    The full boundary is [v_0, ..., v_{h-1}, -v_0, ..., -v_{h-1}],
    so the total area is twice the area of the half polygon closed by -v_0.
    This avoids allocating temporary arrays of length h+1.
    """
    h = xh.shape[0]
    if h < 2:
        return 0.0

    area2_half = 0.0

    x_prev = xh[0]
    y_prev = yh[0]
    for i in range(1, h):
        x_i = xh[i]
        y_i = yh[i]
        area2_half += x_prev * y_i - x_i * y_prev
        x_prev = x_i
        y_prev = y_i

    area2_half += x_prev * (-yh[0]) - (-xh[0]) * y_prev

    if area2_half < 0.0:
        return -area2_half
    return area2_half

@njit(cache=True, fastmath=True)
def polygon_area(xv, yv):
    """
    Compute the area of a simple polygon using the Shoelace formula.

    Parameters
    ----------
    xv : ``numpy.ndarray``
        x-coordinates of the polygon vertices. \n
    yv : ``numpy.ndarray``
        y-coordinates of the polygon vertices. \n

    Returns
    -------
    area : ``float``
        Area of the polygon. \n

    Examples
    --------
    >>> import numpy as np
    >>> area = polygon_area(np.array([0., 1., 0.]), np.array([0., 0., 1.]))
    """
    area = 0.0
    n = xv.shape[0]
    j = n - 1
    for i in range(n):
        area += xv[j] * yv[i] - xv[i] * yv[j]
        j = i
    return abs(area) / 2.0

def make_cross_section_area_reinit(Da_instance, num_th=500, maginf=-100.0):
    """
    Make a jitted function to compute double-caustic cross sections.

    Drop-in replacement with trig arrays closed over (faster),
    safe memory layout, and optional chunking for huge batches.

    Parameters
    ----------
    Da_instance : ``callable``
        Angular diameter distance function. \n
    num_th : ``int``
        Number of points to sample. \n
        default: 500
    maginf : ``float``
        Magnification cut threshold. \n
        default: -100.0

    Returns
    -------
    cross_section_area : ``callable``
        Jitted function to compute cross section areas. \n

    Examples
    --------
    >>> cross_section_fn = make_cross_section_area_reinit(Da_instance)
    """
    # Build only the first half of the theta grid — covers [0, π).
    # This halves memory for all five trig arrays (closed over by the inner
    # function) and removes the [:h] slicing overhead on every call to
    # caustic_area_epl_shear.
    h = num_th // 2
    theta = np.arange(h) * (2.0 * np.pi / num_th)

    # Closed-over constants (no recomputation)
    cos_th = np.cos(theta)
    sin_th = np.sin(theta)
    cos_2th = np.cos(2.0 * theta)
    sin_2th = np.sin(2.0 * theta)

    @njit(parallel=True, cache=False, fastmath=True)
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

@njit(parallel=True, cache=False, fastmath=True)
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

    Parameters
    ----------
    e1 : ``numpy.ndarray``
        First ellipticity component array. \n
    e2 : ``numpy.ndarray``
        Second ellipticity component array. \n
    gamma : ``numpy.ndarray``
        EPL slope exponent array. \n
    gamma1 : ``numpy.ndarray``
        First shear component array. \n
    gamma2 : ``numpy.ndarray``
        Second shear component array. \n
    num_th : ``int``
        Number of points to sample. \n
        default: 500
    maginf : ``float``
        Magnification cut threshold. \n
        default: -100.0

    Returns
    -------
    cross_section_area : ``callable``
        Jitted function to compute cross section areas. \n

    Examples
    --------
    >>> import numpy as np
    >>> calculate_area = cross_section_epl_shear_unit(np.array([0.1]), np.array([0.05]), np.array([2.0]), np.array([0.03]), np.array([-0.01]))
    """
    # Build only the first half of the theta grid — covers [0, π).
    # Halves memory for all five trig arrays and removes [:h] slicing overhead
    # on every caustic_area_epl_shear call inside the serial loop.
    h = num_th // 2
    theta = np.arange(h) * (2.0 * np.pi / num_th)

    phi, q = ellipticity2phi_q(e1, e2)

    cos_th = np.cos(theta)
    sin_th = np.sin(theta)
    cos_2th = np.cos(2.0 * theta)
    sin_2th = np.sin(2.0 * theta)

    size = e1.shape[0]
    cs_area = np.empty(size)
    theta_E = 1.0

    for i in prange(size):

        area = caustic_area_epl_shear(
            q[i], phi[i], gamma[i], gamma1[i], gamma2[i], theta_E,
            theta, cos_th, sin_th, cos_2th, sin_2th,
            maginf
        )

        cs_area[i] = area if np.isfinite(area) else 0.0

    return cs_area
