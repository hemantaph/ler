"""
Numba-JIT-compiled EPL + external shear lensing solver.

Provides analytical image-position finding, magnification computation,
Fermat-potential evaluation, and a parallel batch solver for
elliptical power-law (EPL) lens models with external shear.

The module uses the 1-D lens-equation approach, reducing the 2-D vector
equation to a scalar root-finding problem parameterised by the
image-plane angle. All computationally intensive functions are compiled
with ``@njit`` (Numba) for performance.

The top-level entry points are: \n
- :func:`image_position_analytical_njit` — single-system solver \n
- :func:`create_epl_shear_solver` — factory for a parallel batch solver \n

Copyright (C) 2026 Author Name. Distributed under MIT License.
"""

import numpy as np
from numba import njit, prange
from .sample_caustic_points_njit import sample_source_from_double_caustic
from .cross_section_njit import cdot


EPS = 1e-16
MAX_ROOTS = 16
MAX_IMGS = 5
C_LIGHT = 299792458.0  # m/s

@njit(cache=True, fastmath=True)
def pol_to_ell(r, theta, q):
    """
    Convert polar coordinates to elliptical coordinates.

    Parameters
    ----------
    r : ``float``
        Polar radial coordinate.
    theta : ``float``
        Polar angle in radians.
    q : ``float``
        Axis ratio of the ellipse (0 < q ≤ 1).

    Returns
    -------
    rell : ``float``
        Elliptical radial coordinate,
        ``rell = r * sqrt(q^2*cos^2(theta) + sin^2(theta))``.
    phi : ``float``
        Elliptical angle in radians,
        ``phi = arctan2(sin(theta), cos(theta)*q)``.

    Examples
    --------
    >>> import numpy as np
    >>> from ler.image_properties.epl_shear_njit import pol_to_ell
    >>> rell, phi = pol_to_ell(1.0, np.pi / 4, 0.7)
    """

    phi = np.arctan2(np.sin(theta), np.cos(theta) * q)
    rell = r * np.sqrt(q**2 * np.cos(theta) ** 2 + np.sin(theta) ** 2)
    return rell, phi

@njit(cache=True, fastmath=True)
def omega_scalar(phi, t, q, niter_max=200, tol=1e-16):
    """
    Scalar series expansion of the EPL deflection kernel Omega. 

    Evaluates the convergent series for the EPL Omega function at a
    single elliptical angle, avoiding temporary array allocation.
    The series converges geometrically with ratio
    ``f = (1-q)/(1+q)``.

    Parameters
    ----------
    phi : ``float``
        Elliptical angle in radians.
    t : ``float``
        EPL slope parameter (``t = gamma - 1``).
    q : ``float``
        Axis ratio (0 < q ≤ 1).
    niter_max : ``int``
        Maximum number of series iterations. 
        default: 200
    tol : ``float``
        Convergence tolerance for series truncation. 
        default: 1e-16

    Returns
    -------
    omega_sum : ``complex``
        Complex EPL deflection kernel value at angle ``phi``.

    Examples
    --------
    >>> import numpy as np
    >>> from ler.image_properties.epl_shear_njit import omega_scalar
    >>> omega = omega_scalar(np.pi / 3, t=1.0, q=0.8)
    """
    f = (1.0 - q) / (1.0 + q)

    # Floating-point safeguard: if q is extremely close to 1, f approaches 0.
    # This prevents np.log(f) from evaluating to -inf and crashing the loop counter.
    if f <= EPS:
        return np.exp(1j * phi)

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
def _alpha_epl_shear_scalar(x, y, b, q, t=1, gamma1=0, gamma2=0):
    """
    Scalar complex deflection for EPL + external shear. 

    Same as the array version but avoids temporary allocation by
    calling the scalar :func:`omega_scalar`. 
    ``alpha_EPL   = (2 b^t) / (1 + q) * R^(1-t) * Omega`` 
    ``alpha_Shear = gamma1*x + gamma2*y + i*(gamma2*x - gamma1*y)``

    Parameters
    ----------
    x : ``float``
        Image-plane x-coordinate.
    y : ``float``
        Image-plane y-coordinate.
    b : ``float``
        EPL scale parameter (``b = theta_E * sqrt(q)``).
    q : ``float``
        Axis ratio (0 < q ≤ 1).
    t : ``float``
        EPL slope parameter (``t = gamma - 1``). 
        default: 1
    gamma1 : ``float``
        External shear component 1. 
        default: 0
    gamma2 : ``float``
        External shear component 2. 
        default: 0

    Returns
    -------
    alpha : ``complex``
        Total complex deflection angle (EPL + shear).
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


@njit(cache=True, fastmath=True)
def _unique_points(x, y, n, tol):
    """
    Remove near-duplicate image positions in-place. 

    Iterates over the first ``n`` entries of ``x`` and ``y``,
    marking duplicates as inactive and compressing the arrays.

    Parameters
    ----------
    x : ``numpy.ndarray``
        Image x-coordinates (modified in place).
    y : ``numpy.ndarray``
        Image y-coordinates (modified in place).
    n : ``int``
        Number of active entries to check.
    tol : ``float``
        Distance threshold below which two points are considered identical.

    Returns
    -------
    m : ``int``
        Number of unique points remaining after deduplication.
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


@njit(cache=True, fastmath=True)
def _insertion_sort5(x, y, t, mu, itype, n):
    """
    Sort the first ``n`` image records by ascending arrival time. 

    Operates in-place on all five parallel arrays simultaneously
    using insertion sort.

    Parameters
    ----------
    x : ``numpy.ndarray``
        Image x-coordinates.
    y : ``numpy.ndarray``
        Image y-coordinates.
    t : ``numpy.ndarray``
        Arrival times (Fermat potentials) used as the sort key.
    mu : ``numpy.ndarray``
        Signed magnifications.
    itype : ``numpy.ndarray``
        Image-type codes.
    n : ``int``
        Number of elements to sort.
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


@njit(cache=True, fastmath=True)
def _compress_keep5(mask, x, y, t, mu, itype, n):
    """
    Compress five parallel image arrays using a boolean mask in-place. 

    Retains only entries where ``mask[i] == 1`` and shifts them to
    the front of the arrays.

    Parameters
    ----------
    mask : ``numpy.ndarray``
        Boolean-like ``uint8`` array; 1 = keep, 0 = discard.
    x : ``numpy.ndarray``
        Image x-coordinates.
    y : ``numpy.ndarray``
        Image y-coordinates.
    t : ``numpy.ndarray``
        Arrival times (Fermat potentials).
    mu : ``numpy.ndarray``
        Signed magnifications.
    itype : ``numpy.ndarray``
        Image-type codes.
    n : ``int``
        Number of elements to process.

    Returns
    -------
    m : ``int``
        Number of retained entries after compression.
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


@njit(cache=True, fastmath=True)
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


@njit(cache=True, fastmath=True)
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

# ----------------
# Image Positions
# ----------------
@njit(cache=True, fastmath=True)
def _one_dim_lens_eq_unsmooth(phi, args):
    """
    Non-smooth branch of the 1-D lens equation at a single angle. 

    Used when the smooth branch is numerically unstable due to sign
    ambiguity in the cross-product term ``ip``. The branch uses
    ``|ip|^t`` instead of ``ps(ip, 2/t)`` to avoid fractional-power
    sign issues.

    Parameters
    ----------
    phi : ``float``
        Trial image-plane angle in radians.
    args : ``tuple``
        Packed lens parameters ``(b, t, y1, y2, q, gamma1, gamma2)``.

    Returns
    -------
    eq_notsmooth : ``float``
        Residual of the non-smooth lens equation at angle ``phi``.
    """
    # lens eqn: y = (1-Gamma) r.rhat - lambda(R) Omega(phiell), where lambda(R) = const (R/b)^(1-t), b = theta_E * sqrt(q) 
    # r = R * sqrt(cos^2(phiell)*q^2 + sin^2(phiell))
    # phiell = arctan(sin(phi)*q/cos(phi))
    Omega, const, phiell, q, r, rhat, t, b, thetahat, y = _one_dim_lens_eq_calcs(
        args, phi
    )
    # rr = r * sqrt(cos^2(phiell)*q^2 + sin^2(phiell)); r=1
    rr, thetaa = _ell_to_pol(1, phiell, q)
    ip = cdot(y, rhat) * cdot(Omega, thetahat) - cdot(Omega, rhat) * cdot(y, thetahat)
    eq_notsmooth = _ps(rr * b, 1 - t) * (cdot(y, thetahat) / const) * np.abs(
        ip
    ) ** t + ip * np.abs(cdot(Omega, thetahat)) ** (+t)
    return eq_notsmooth

@njit(cache=True, fastmath=True)
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
    if not abs(div) > EPS:
        return x2
    num = x3 * x3 * (y1 - y2) + x1 * x1 * (y2 - y3) + x2 * x2 * (-y1 + y3)
    return num / div


@njit(cache=True, fastmath=True)
def _brentq_one_dim_lens_residual(
    xa,
    xb,
    args,
    smooth_branch=True,
    xtol=2e-14,
    rtol=16 * np.finfo(np.float64).eps,
    maxiter=100,
):
    """
    Brent root finder for the 1-D EPL+shear lens equation (scalar residual).

    This is the same Brent-q logic as previously used for the 1-D residual,
    but **does not** accept a callable argument so the implementation is
    pickle-cacheable under ``@njit(cache=True)``
    (no first-class Numba dispatchers in the ABI).

    The residual is taken from ``_one_dim_lens_eq`` (smooth branch) or
    ``_one_dim_lens_eq_unsmooth`` (non-smooth branch) depending on
    ``smooth_branch``.

    Returns ``np.nan`` if the bracket is invalid or numerically unstable.

    Parameters
    ----------
    xa : ``float``
        Left bracket boundary.
    xb : ``float``
        Right bracket boundary.
    args : ``tuple``
        Packed lens parameters ``(b, t, y1, y2, q, gamma1, gamma2)``.
    smooth_branch : ``bool``
        If True, evaluate ``_one_dim_lens_eq``; otherwise
        ``_one_dim_lens_eq_unsmooth``.
    xtol : ``float``
        Absolute tolerance.
    rtol : ``float``
        Relative tolerance.
    maxiter : ``int``
        Maximum iterations.

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

    if smooth_branch:
        fpre = _one_dim_lens_eq(xpre, args)
        fcur = _one_dim_lens_eq(xcur, args)
    else:
        fpre = _one_dim_lens_eq_unsmooth(xpre, args)
        fcur = _one_dim_lens_eq_unsmooth(xcur, args)

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

            if (not use_bisect) and (
                2.0 * abs(stry) < min(abs(spre), 3.0 * abs(sbis) - delta)
            ):
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

        if smooth_branch:
            fcur = _one_dim_lens_eq(xcur, args)
        else:
            fcur = _one_dim_lens_eq_unsmooth(xcur, args)

        if np.isnan(fcur):
            return np.nan

    return xcur


@njit(cache=True, fastmath=True)
def _getr(phi, args):
    """
    Compute the image-plane radial distance for a trial angle ``phi``.

    Parameters
    ----------
    phi : ``float``
        Image-plane angle in radians.
    args : ``tuple``
        Packed lens parameters ``(b, t, y1, y2, q, gamma1, gamma2)``.

    Returns
    -------
    r : ``float``
        Estimated image-plane radial distance at angle ``phi``.
    """
    _, _, _, _, r, _, _, _, _, _ = _one_dim_lens_eq_calcs(args, phi)
    return r


@njit(cache=True, fastmath=True)
def _one_dim_lens_eq(phi, args):
    """
    Smooth branch of the 1-D lens equation at a single angle. 

    Call chain:
    ``_solvelenseq_majoraxis -> _getphi -> _one_dim_lens_eq_both -> _one_dim_lens_eq``.

    The residual is: 
    ``(rr*b)^(2/t-2) * ps(y.thetahat/const, 2/t) * ip^2
    + ps(ip, 2/t) * (Omega.thetahat)^2 = 0``

    Parameters
    ----------
    phi : ``float``
        Trial image-plane angle in radians.
    args : ``tuple``
        Packed lens parameters ``(b, t, y1, y2, q, gamma1, gamma2)``.

    Returns
    -------
    eq : ``float``
        Residual of the smooth lens equation at angle ``phi``.
    """
    # lens eqn: y = (1-Gamma) r.rhat - lambda(R) Omega(phiell), where lambda(R) = const (R/b)^(1-t), b = theta_E * sqrt(q) 
    # r = R * sqrt(cos^2(phiell)*q^2 + sin^2(phiell))
    # phiell = arctan(sin(phi)*q/cos(phi))
    Omega, const, phiell, q, r, rhat, t, b, thetahat, y = _one_dim_lens_eq_calcs(
        args, phi
    )
    # rr = r * sqrt(cos^2(phiell)*q^2 + sin^2(phiell)); r=1
    rr, _ = _ell_to_pol(1.0, phiell, q)

    # multiply thetahat to eliminate the unknown image-plane radius in the lens equation, and get expression for lambda(R)
    # multiply rhat to get the expression for the lens equation, get scaler image positions, and use lambda(R)
    # ip = r. Omega . thetahat = (beta . rhat) * (Omega . thetahat) - (Omega . rhat) * (beta . thetahat) 
    ip = cdot(y, rhat) * cdot(Omega, thetahat) - cdot(Omega, rhat) * cdot(y, thetahat)
    # using lambda(R) = const (R/b)^(1-t), 
    # (R/b)^(1-t) = - (beta . thetahat) / (const * Omega . thetahat)
    # probably b is scaled with rr.
    # A*ip^2 + B*ip^(2/t) = 0.
    eq = (rr * b) ** (2 / t - 2) * _ps(
        (cdot(y, thetahat) / const), 2 / t
    ) * ip**2 + _ps(ip, 2 / t) * cdot(Omega, thetahat) ** 2

    return eq

@njit(cache=True, fastmath=True)
def _one_dim_lens_eq_calcs(args, phi):
    """
    Compute intermediate geometric quantities for the 1-D lens equation. 

    Evaluates the shear-corrected unit vectors ``rhat`` and ``thetahat``,
    the EPL deflection kernel ``Omega``, and the estimated image radius
    ``r`` for a trial angle ``phi``. The function does not compute the
    physical image radius; it only computes the geometric conversion
    factors used in both smooth and non-smooth branches of the lens
    equation. The actual image radius is computed separately in each
    branch after eliminating the unknown radius via the dot product
    with ``thetahat``. 

    Call chain: 
    ``_solvelenseq_majoraxis -> _getphi -> _one_dim_lens_eq_both
    -> _one_dim_lens_eq -> _one_dim_lens_eq_calcs``.

    Parameters
    ----------
    args : ``tuple``
        Packed lens parameters ``(b, t, y1, y2, q, gamma1, gamma2)``.
    phi : ``float``
        Trial image-plane angle in radians.

    Returns
    -------
    Omega : ``complex``
        EPL deflection kernel value at the elliptical angle.
    const : ``float``
        EPL amplitude constant ``2*b/(1+q)``.
    phiell : ``float``
        Elliptical angle corresponding to ``phi``.
    q : ``float``
        Axis ratio.
    r : ``float``
        Estimated image-plane radial distance.
    rhat : ``complex``
        Shear-corrected radial unit vector.
    t : ``float``
        EPL slope parameter.
    b : ``float``
        EPL scale parameter.
    thetahat : ``complex``
        Shear-corrected tangential unit vector.
    y : ``complex``
        Source position as a complex number (``y1 + i*y2``).
    """
    b, t, y1, y2, q, gamma1, gamma2 = args
    # Note: y1 and y2 are the source-plane coordinates in the rotated frame, not the original source-plane coordinates.
    y = y1 + 1j * y2

    # rhat' : with shear, rhat : without shear [cos(phi), sin(phi)]
    # rhat' = M^-1 . rhat, where M^-1 = [[M22, -M12], [-M21, M11]] / det(M)
    # M = [[1 + gamma1, gamma2], [gamma2, 1 - gamma1]]
    # det(M) = 1 - gamma1^2 - gamma2^2 
    den = 1.0 - gamma1 * gamma1 - gamma2 * gamma2 
    # if abs(den) < EPS:
    #     den = EPS if den >= 0.0 else -EPS
    c = np.cos(phi)
    s = np.sin(phi)

    rhat = (
        ((1.0 + gamma1) * c + gamma2 * s)
        + 1j * (gamma2 * c + (1.0 - gamma1) * s)
    ) / den

    # construct an unit vector thetahat perpendicular to rhat, with the same shear transformation applied.
    # thetahat is used to eliminate the unknown image-plane radius in the lens equation, by taking the dot product with the source-plane vector y and the deflection vector Omega.
    # thetahat' : with shear, thetahat : without shear [-sin(phi), cos(phi)]
    # but code uses [sin(phi), -cos(phi)]
    thetahat = (
        ((1.0 + gamma1) * s - gamma2 * c)
        + 1j * (gamma2 * s - (1.0 - gamma1) * c)
    ) / den

    # convert the trial image angle phi into elliptical angle
    # frac_Roverrsh = (R/r), r=1
    # R = sqrt(cos^2(phi)*q^2 + sin^2(phi))
    # phiell = arctan(sin(phi)*q/cos(phi))
    frac_Roverrsh, phiell = pol_to_ell(1.0, phi, q)

    # alpha_epl = (2.0*b/(1.0+q)) * (b/frac_Roverrsh)**(t-1) * Omega(phiell)
    Omega = omega_scalar(phiell, t, q)
    const = 2.0 * b / (1.0 + q)

    if abs(t - 1.0) > 1e-4: # avoid numerical instability near isothermal, gamma = 2.0

        # after taking the dot product with thetahat, the lens equation becomes a scalar equation in phi, with the unknown image-plane radius eliminated:
        # R = b * |y . thetahat / (const * Omega . thetahat)|^(1/(1-t)) * sign(y . thetahat / (const * Omega . thetahat))
        val = -cdot(y, thetahat) / (const * cdot(Omega, thetahat))
        # aval = abs(val)
        # if aval < EPS:
        #     R = 0.0
        # else:
        # R = b * abs(val) ** (1.0 / (1.0 - t)) * np.sign(val)
        R = b * _ps(val, 1.0 / (1.0 - t))
    else:
        # isothermal case, gamma = 2.0, t = 1.0, the lens equation becomes:
        # construct an vector Omega_ort perpendicular to Omega, with the same shear transformation applied.
        Omega_ort = 1j * Omega
        x = ((1.0 - gamma1) * c - gamma2 * s) + 1j * (
            -gamma2 * c + (1.0 + gamma1) * s
        )
        denom = cdot(Omega_ort, x)
        if abs(denom) > EPS:
            R = cdot(Omega_ort, y) / denom * frac_Roverrsh
        else:
            # Fallback for SIS (q ≈ 1): cdot(Omega_ort, x) → 0 because
            # Omega_ort ⊥ x when the elliptical and polar angles coincide.
            # Project the isothermal lens equation y = r·den·rhat − const·Omega
            # onto rhat to solve for r directly, then convert to R.
            rhat_sq = cdot(rhat, rhat)
            if abs(rhat_sq) > EPS:
                r_direct = (cdot(y, rhat) + const * cdot(Omega, rhat)) / (den * rhat_sq)
                R = r_direct * frac_Roverrsh
            else:
                R = 0.0

    # r = R * sqrt(cos^2(phiell)*q^2 + sin^2(phiell))
    # phiell = arctan(sin(phiell)*q/cos(phiell))
    r, _ = _ell_to_pol(R, phiell, q)

    return Omega, const, phiell, q, r, rhat, t, b, thetahat, y

@njit(cache=True, fastmath=True)
def _one_dim_lens_eq_both(phi, args):
    """
    Evaluate both branches of the 1-D lens equation on a phi grid. 

    Computes the smooth and non-smooth residuals at every angle in
    ``phi`` and returns them as two parallel arrays. 

    Call chain:
    ``_solvelenseq_majoraxis -> _getphi -> _one_dim_lens_eq_both``.

    Parameters
    ----------
    phi : ``numpy.ndarray``
        Array of trial image-plane angles in radians.
    args : ``tuple``
        Packed lens parameters ``(b, t, y1, y2, q, gamma1, gamma2)``.

    Returns
    -------
    eq : ``numpy.ndarray``
        Smooth-branch residuals at each angle.
    eq_notsmooth : ``numpy.ndarray``
        Non-smooth-branch residuals at each angle.
    """
    n = phi.shape[0]
    eq = np.empty(n, dtype=np.float64)
    eq_notsmooth = np.empty(n, dtype=np.float64)

    for i in range(n):
        # lens eqn: y = (1-Gamma) r.rhat - lambda(R) Omega(phiell), where lambda(R) = const (R/b)^(1-t), b = theta_E * sqrt(q) 
        # r = R * sqrt(cos^2(phiell)*q^2 + sin^2(phiell))
        # phiell = arctan(sin(phi)*q/cos(phi))
        Omega, const, phiell, q, r, rhat, t, b, thetahat, y = _one_dim_lens_eq_calcs(args, phi[i])
        # rr = r * sqrt(cos^2(phiell)*q^2 + sin^2(phiell)); r=1
        rr, thetaa = _ell_to_pol(1.0, phiell, q)

        # multiply thetahat to eliminate the unknown image-plane radius in the lens equation, and get expression for lambda(R)
        # multiply rhat to get the expression for the lens equation, get scaler image positions, and use lambda(R)
        # ip = r. Omega . thetahat = (beta . rhat) * (Omega . thetahat) - (Omega . rhat) * (beta . thetahat) 
        ip = cdot(y, rhat) * cdot(Omega, thetahat) - cdot(Omega, rhat) * cdot(y, thetahat)
        # using lambda(R) = const (R/b)^(1-t), 
        # (R/b)^(1-t) = - (beta . thetahat) / (const * Omega . thetahat)
        # probably b is scaled with rr.
        # A*ip^2 + B*ip^(2/t) = 0.
        eq[i] = (rr * b) ** (2.0 / t - 2.0) * _ps(cdot(y, thetahat) / const, 2.0 / t) * ip * ip + _ps(ip, 2.0 / t) * (cdot(Omega, thetahat) ** 2)

        eq_notsmooth[i] = _ps(rr * b, 1 - t) * (cdot(y, thetahat) / const) * np.abs(
            ip
        ) ** t + ip * np.abs(cdot(Omega, thetahat)) ** (+t)

    return eq, eq_notsmooth

@njit(cache=True, fastmath=True)
def _getphi(thpl, args):
    """
    Find all angular roots of the 1-D lens equation on the supplied grid. 

    Scans both smooth and non-smooth branches for sign changes and
    local extrema,     then refines each bracket with :func:`_brentq_one_dim_lens_residual`.
    Handles the case where a root sits at an
    extremum (double root / tangent crossing) via parabolic
    interpolation. 

    Call chain: ``_solvelenseq_majoraxis -> _getphi``.

    Parameters
    ----------
    thpl : ``numpy.ndarray``
        Sorted array of trial angles in radians.
    args : ``tuple``
        Packed lens parameters ``(b, t, y1, y2, q, gamma1, gamma2)``.

    Returns
    -------
    roots : ``numpy.ndarray``
        Array of root angles (length ``MAX_ROOTS``; only first
        ``nroots`` entries are valid).
    nroots : ``int``
        Number of valid roots found.
    """

    # Evaluate both branches of the lens equation on the grid to find sign changes.
    # y = A * ip^2 + B * ip^(2/t) = 0, where A and B are functions of phi, and ip is a function of phi. The non-smooth branch is the one with the absolute value of ip, which can have a different sign change pattern than the smooth branch.
    y, y_ns = _one_dim_lens_eq_both(thpl, args)
    nphi = thpl.shape[0]

    roots = np.empty(MAX_ROOTS, dtype=np.float64)
    nroots = 0

    for i in range(nphi - 1):
        # Check for sign changes in both branches of the lens equation to identify root brackets.
        if y[i + 1] * y[i] <= 0.0:
            if nroots < MAX_ROOTS:
                roots[nroots] = _brentq_one_dim_lens_residual(
                    thpl[i], thpl[i + 1], args,
                ) % (2 * np.pi)
                nroots += 1
        elif y_ns[i + 1] * y_ns[i] <= 0.0:
            if nroots < MAX_ROOTS:
                roots[nroots] = _brentq_one_dim_lens_residual(
                    thpl[i], thpl[i + 1], args, False,
                ) % (2 * np.pi)
                nroots += 1

    for i in range(1, nphi+1):
        y1 = y[i - 1]
        y2 = y[i % nphi]
        y3 = y[(i + 1) % nphi]

        y1n = y_ns[i - 1]
        y2n = y_ns[i % nphi]
        y3n = y_ns[(i + 1) % nphi]

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
                    roots[nroots] = _brentq_one_dim_lens_residual(x2, xmin, args) % (
                        2.0 * np.pi
                    )
                    nroots += 1
                if nroots < MAX_ROOTS:
                    roots[nroots] = _brentq_one_dim_lens_residual(xmin, x3, args) % (
                        2.0 * np.pi
                    )
                    nroots += 1

            elif ymin * y2 <= 0.0 and x1 <= xmin <= x2:

                if nroots < MAX_ROOTS:
                    roots[nroots] = _brentq_one_dim_lens_residual(x1, xmin, args) % (
                        2.0 * np.pi
                    )
                    nroots += 1
                if nroots < MAX_ROOTS:
                    roots[nroots] = _brentq_one_dim_lens_residual(xmin, x2, args) % (
                        2.0 * np.pi
                    )
                    nroots += 1

            elif ymin_ns * y2n <= 0.0 and x2 <= xmin_ns <= x3:

                if nroots < MAX_ROOTS:
                    roots[nroots] = _brentq_one_dim_lens_residual(
                        x2, xmin_ns, args, False,
                    ) % (2.0 * np.pi)
                    nroots += 1
                if nroots < MAX_ROOTS:
                    roots[nroots] = _brentq_one_dim_lens_residual(
                        xmin_ns, x3, args, False,
                    ) % (2.0 * np.pi)
                    nroots += 1

            elif ymin_ns * y2n <= 0.0 and x1 <= xmin_ns <= x2:

                if nroots < MAX_ROOTS:
                    roots[nroots] = _brentq_one_dim_lens_residual(
                        x1, xmin_ns, args, False,
                    ) % (2.0 * np.pi)
                    nroots += 1
                if nroots < MAX_ROOTS:
                    roots[nroots] = _brentq_one_dim_lens_residual(
                        xmin_ns, x2, args, False,
                    ) % (2.0 * np.pi)
                    nroots += 1

    return roots, nroots


@njit(cache=True, fastmath=True)
def _solvelenseq_majoraxis(b, t, y1, y2, q, gamma1, gamma2, Nmeas=200, Nmeas_extra=50):
    """
    Solve the lens equation in the lens major-axis-aligned frame. 

    Builds an angular grid with extra refinement near the source
    direction, finds all angular roots via :func:`_getphi`, then
    reconstructs 2-D image positions and verifies each candidate
    against the full deflection equation.

    Parameters
    ----------
    b : ``float``
        EPL scale parameter (``b = theta_E * sqrt(q)``).
    t : ``float``
        EPL slope parameter (``t = gamma - 1``).
    y1 : ``float``
        Source x-coordinate in the axis-aligned frame.
    y2 : ``float``
        Source y-coordinate in the axis-aligned frame.
    q : ``float``
        Axis ratio (0 < q ≤ 1).
    gamma1 : ``float``
        External shear component 1 (axis-aligned frame).
    gamma2 : ``float``
        External shear component 2 (axis-aligned frame).
    Nmeas : ``int``
        Number of uniformly-spaced angle samples. 
        default: 200
    Nmeas_extra : ``int``
        Number of extra refinement samples near the source direction. 
        default: 50

    Returns
    -------
    xsol : ``numpy.ndarray``
        Image x-coordinates (length ``MAX_IMGS``).
    ysol : ``numpy.ndarray``
        Image y-coordinates (length ``MAX_IMGS``).
    nimg : ``int``
        Number of valid images found.
    """
    # alpha_shear = [[gamma1, gamma2], [gamma2, -gamma1]] * theta = Gamma * theta
    # alpha_epl = lambda(R) * Omega(phi_ell)
    # lens equation with shear:
    # beta = (1-Gamma) * theta - alpha_epl(theta); (1-Gamma) = M
    # M^-1 * beta = theta - M^-1 * alpha_epl(theta)
    # beta_shear = [[1+gamma1, gamma2], [gamma2, 1-gamma1]] . [y1, y2] = r * exp(i*p1)
    p1 = np.arctan2(
        y2 * (1.0 - gamma1) + gamma2 * y1,
        y1 * (1.0 + gamma1) + gamma2 * y2
    ) # direction angle of the effective source position after shear correction

    # Create a grid of angles for root finding, with extra sampling near the source angle.
    # 0 to 1e-4 for the linear part, then geometric spacing up to 0.1.
    geom = _geomlinspace(1e-4, 0.1, Nmeas_extra)
    ngeom = geom.shape[0]
    # Use the actual geometric grid length (ngeom), which can exceed Nmeas_extra.
    ntot = Nmeas + 2 * ngeom

    thpl = np.empty(ntot, dtype=np.float64)

    for i in range(Nmeas):
        # Uniformly spaced angles from 0 to pi.
        # the equation has π-periodicity due to the quadratic nature of shear + EPL
        thpl[i] = np.pi * i / (Nmeas - 1)

    pmod = p1 % np.pi # mod by pi to get the angle in the range [0, pi)
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

            # alpha_epl = lambda(R) * Omega(phi_ell)
            diff = (-y1 - 1j * y2) + (xs + 1j * ys) - _alpha_epl_shear_scalar(
                xs, ys, b, q, t, gamma1=gamma1, gamma2=gamma2
            )

            # 1e-8 is a tight tolerance to ensure we only keep valid roots, it directly follows the lenstronomy code.
            if abs(diff) < 1e-8:
                if nimg < MAX_IMGS:
                    xsol[nimg] = xs
                    ysol[nimg] = ys
                    nimg += 1

    nimg = _unique_points(xsol, ysol, nimg, 1e-8)
    return xsol, ysol, nimg

@njit(cache=True, fastmath=True)
def _solve_lenseq_pemd(xsrc, ysrc, q, phi, t, gamma1, gamma2,
                      b, Nmeas=400, Nmeas_extra=80):
    """
    Solve the EPL + shear lens equation in the rotated major-axis frame. 

    Rotates source coordinates and shear components into the lens
    major-axis frame, calls :func:`_solvelenseq_majoraxis`, then
    rotates image positions back to the original sky frame.

    Parameters
    ----------
    xsrc : ``float``
        Source x-coordinate in sky frame.
    ysrc : ``float``
        Source y-coordinate in sky frame.
    q : ``float``
        Axis ratio (0 < q ≤ 1).
    phi : ``float``
        Lens position angle in radians.
    t : ``float``
        EPL slope parameter (``t = gamma - 1``).
    gamma1 : ``float``
        External shear component 1.
    gamma2 : ``float``
        External shear component 2.
    b : ``float``
        EPL scale parameter (``b = theta_E * sqrt(q)``).
    Nmeas : ``int``
        Number of angle samples for root finding. 
        default: 400
    Nmeas_extra : ``int``
        Number of extra refinement samples. 
        default: 80

    Returns
    -------
    xout : ``numpy.ndarray``
        Image x-coordinates in sky frame (length ``MAX_IMGS``).
    yout : ``numpy.ndarray``
        Image y-coordinates in sky frame (length ``MAX_IMGS``).
    nimg : ``int``
        Number of valid images found.
    """
    # t = gamma - 1.0
    # b = theta_E * np.sqrt(q)

    # rotate source and shear into the axis-aligned frame
    # p = (xsrc + 1j * ysrc) * np.exp(-1j * phi) # expand
    # g = (gamma1 + 1j * gamma2) * np.exp(-1j * 2*phi) # expand
    cp = np.cos(phi)
    sp = np.sin(phi)
    c2p = cp * cp - sp * sp   # cos(2*phi)
    s2p = 2.0 * cp * sp       # sin(2*phi)
    x_rot =  xsrc * cp + ysrc * sp
    y_rot = -xsrc * sp + ysrc * cp
    gamma1_rot =  gamma1 * c2p + gamma2 * s2p
    gamma2_rot = -gamma1 * s2p + gamma2 * c2p

    xloc, yloc, nimg = _solvelenseq_majoraxis(
        b, # scale parameter for the lens
        t, # power-law slope minus one
        x_rot, # source position rotated into the lens frame
        y_rot, # source position rotated into the lens frame
        q, # axis ratio of the lens
        gamma1_rot, # shear component 1 rotated into the lens frame
        gamma2_rot, # shear component 2 rotated into the lens frame
        Nmeas, # number of angular samples for root finding
        Nmeas_extra # extra angular samples near the source angle for better root finding
    )

    xout = np.empty(MAX_IMGS, dtype=np.float64)
    yout = np.empty(MAX_IMGS, dtype=np.float64)

    for i in range(nimg):
        xout[i] = xloc[i] * cp - yloc[i] * sp
        yout[i] = xloc[i] * sp + yloc[i] * cp

    return xout, yout, nimg
# ---------------

# ---------------
# Magnifications
# ---------------
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
    # Gamma = [[gamma1, gamma2], [gamma2, -gamma1]]
    return f_xx, f_xy, f_xy, f_yy

@njit(cache=True, fastmath=True)
def _epl_hessian(z, b, t, q, phi, Omega):
    """
    Compute the EPL Hessian components at a complex image position. 

    Returns the second derivatives of the EPL lensing potential
    in the axis-aligned frame.

    Parameters
    ----------
    z : ``complex``
        Image position in the axis-aligned frame
        (``z = exp(-i*phi) * (x + i*y)``).
    b : ``float``
        EPL scale parameter (``b = theta_E * sqrt(q)``).
    t : ``float``
        EPL slope parameter (``t = gamma - 1``).
    q : ``float``
        Axis ratio (0 < q ≤ 1).
    phi : ``float``
        Lens position angle in radians.
    Omega : ``complex``
        EPL deflection kernel at the image elliptical angle.

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

    # t = gamma - 1.0
    # z = np.exp(-1j * phi) * (x + 1j * y); rotated complex coordinate
    # b = theta_E * np.sqrt(q)

    # Return the angle of the complex argument.
    # ang is the angle of z, at the axis aligned frame.
    ang = np.angle(z)
    # z in elliptical frame coordinates
    zz_ell = z.real * q + 1j * z.imag
    # R is the elliptical radius
    # R = np.sqrt((z.real * q) ** 2 + z.imag**2)
    # R = r * np.sqrt(cos^2(phi)*q^2 + sin^2(phi)); phi not phiell
    R = np.abs(zz_ell)
    # phi_ell is the angle in the elliptical frame
    # phi_ell = np.angle(zz_ell)

    # x, y are scalars
    # The deflection magnitude factor u = (b/R)^t is computed with safeguards to prevent overflow when R is very small.
    if R > EPS:
        u = (b / R) ** t
        if u > 1.0e10:
            u = 1.0e10
    else:
        u = 1.0e10
    
    # lensing convergence kappa
    #  ∇²ψ = 2 * kappa; kappa = 0.5 * (2 - t) * (b/R)^(t)
    # ψ_xx + ψ_yy = 2 * kappa; kappa = 0.5 * (2 - t) * u
    kappa = 0.5 * (2.0 - t)
    # R/r = sqrt(cos^2(phi)*q^2 + sin^2(phi)); phi not phiell
    Roverr = np.sqrt((np.cos(ang) ** 2) * q * q + np.sin(ang) ** 2)

    # Omega = omega_scalar(phi_ell, t, q)
    # deflection: alpha_EPL = (2 b)/(1+q) * (b/R)^t * (R/b) * Omega = (2 b)/(1+q) * u * Roverr * Omega
    alph = (2.0 / (1.0 + q)) * Omega

    # internal shear 
    gamma_shear = (
        -np.exp(2j * (ang + phi)) * kappa
        + (1.0 - t) * np.exp(1j * (ang + 2.0 * phi)) * alph * Roverr
    )

    f_xx = (kappa + gamma_shear.real) * u
    f_yy = (kappa - gamma_shear.real) * u
    f_xy = gamma_shear.imag * u

    return f_xx, f_xy, f_xy, f_yy

@njit(cache=True, fastmath=True)
def _hessian_scalar(z, b, t, gamma1, gamma2, q, phi, Omega):
    """
    Compute the total lensing Hessian (EPL + external shear) at one point.

    Parameters
    ----------
    z : ``complex``
        Image position in the axis-aligned frame.
    b : ``float``
        EPL scale parameter (``b = theta_E * sqrt(q)``).
    t : ``float``
        EPL slope parameter (``t = gamma - 1``).
    gamma1 : ``float``
        External shear component 1.
    gamma2 : ``float``
        External shear component 2.
    q : ``float``
        Axis ratio (0 < q ≤ 1).
    phi : ``float``
        Lens position angle in radians.
    Omega : ``complex``
        EPL deflection kernel at the image elliptical angle.

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
    # Jacobian matrix
    # A = I - Gamma - H
    # H is the Hessian of the EPL potential
    # H = [[f_xx_e, f_xy_e], [f_yx_e, f_yy_e]]
    f_xx_e, f_xy_e, f_yx_e, f_yy_e = _epl_hessian(
        z, b, t, q, phi, Omega
    )
    # Gamma = [[gamma1, gamma2], [gamma2, -gamma1]]
    f_xx_s, f_xy_s, f_yx_s, f_yy_s = _shear_hessian(gamma1, gamma2)

    return (
        f_xx_e + f_xx_s,
        f_xy_e + f_xy_s,
        f_yx_e + f_yx_s,
        f_yy_e + f_yy_s,
    )

@njit(cache=True, fastmath=True)
def lensing_diagnostics_scalar(z, b, t, gamma1, gamma2, q, phi, Omega):
    """
    Compute magnification and image type at one image position. 

    Evaluates the total Hessian of the lensing potential (EPL +
    external shear), then derives the Jacobian determinant and trace
    to classify the image and compute its signed magnification.
    Inputs ``z``, ``b``, ``t``, ``q``, ``phi``, and ``Omega`` are all
    scalars, not arrays.

    Image type convention: 
    - 1: Type I  (minimum of the Fermat potential) 
    - 2: Type II (saddle point) 
    - 3: Type III (maximum of the Fermat potential) 
    - 0: undefined / degenerate (on a critical curve)

    Parameters
    ----------
    z : ``complex``
        Image position in the axis-aligned frame
        (``z = exp(-i*phi) * (x + i*y)``).
    b : ``float``
        EPL scale parameter (``b = theta_E * sqrt(q)``).
    t : ``float``
        EPL slope parameter (``t = gamma - 1``).
    gamma1 : ``float``
        External shear component 1.
    gamma2 : ``float``
        External shear component 2.
    q : ``float``
        Axis ratio (0 < q ≤ 1).
    phi : ``float``
        Lens position angle in radians.
    Omega : ``complex``
        EPL deflection kernel at the image elliptical angle.

    Returns
    -------
    mu : ``float``
        Signed magnification ``1/det(A)``; ``np.inf`` on a critical curve.
    image_type : ``int``
        Image-type code (0, 1, 2, or 3).

    Examples
    --------
    >>> import numpy as np
    >>> from ler.image_properties.epl_shear_njit import (
    ...     lensing_diagnostics_scalar, omega_scalar
    ... )
    >>> z = np.exp(-1j * 0.0) * (0.8 + 1j * 0.3)
    >>> phi_ell = np.angle(z.real * 0.8 + 1j * z.imag)
    >>> Omega = omega_scalar(phi_ell, t=1.0, q=0.8)
    >>> mu, itype = lensing_diagnostics_scalar(
    ...     z, b=0.894, t=1.0, gamma1=0.0, gamma2=0.0, q=0.8, phi=0.0, Omega=Omega
    ... )
    """
    # Jacobian matrix
    # A = I - Gamma - H
    # hessian_scalar returns the total Hessian (EPL + external shear)
    f_xx, f_xy, f_yx, f_yy = _hessian_scalar(
        z, b, t, gamma1, gamma2, q, phi, Omega
    )
    # determinant and trace of the Jacobian
    detA = (1.0 - f_xx) * (1.0 - f_yy) - f_xy * f_yx
    traceA = 2.0 - f_xx - f_yy

    if abs(detA) < EPS:
        mu = 10000.0 # effectively infinite magnification on the critical curve
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

    # return f_xx, f_xy, f_yx, f_yy, detA, traceA, mu, image_type
    return mu, image_type
# -------------------------------------------------

# ------------------
# Fermat potentials
# ------------------
@njit(cache=True, fastmath=True)
def _shear_function(x, y, gamma1, gamma2):
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

    Returns
    -------
    psi_shear : ``float``
        Shear potential value.
    """

    return 0.5 * (gamma1 * x * x + 2.0 * gamma2 * x * y - gamma1 * y * y)

@njit(cache=True, fastmath=True)
def _epl_function(x, y, b, t, q, phi, Omega):
    """
    Compute the EPL lensing potential at a point in the axis-aligned frame.

    Parameters
    ----------
    x : ``float``
        Image x-coordinate in the axis-aligned frame.
    y : ``float``
        Image y-coordinate in the axis-aligned frame.
    b : ``float``
        EPL scale parameter (``b = theta_E * sqrt(q)``).
    t : ``float``
        EPL slope parameter (``t = gamma - 1``).
    q : ``float``
        Axis ratio (0 < q ≤ 1).
    phi : ``float``
        Lens position angle in radians (unused; retained for interface symmetry).
    Omega : ``complex``
        EPL deflection kernel at the image elliptical angle.

    Returns
    -------
    psi_epl : ``float``
        EPL lensing potential value.
    """
    # t = gamma - 1.0
    # z = np.exp(-1j * phi) * (x + 1j * y)
    # b = theta_E * np.sqrt(q)

    # R is the elliptical radius
    # R = np.sqrt((z.real * q) ** 2 + z.imag**2)
    # R = r * np.sqrt(cos^2(phi)*q^2 + sin^2(phi)); phi not phiell

    # this x and y are in the axis-aligned frame
    R = np.abs(x * q + 1j * y)

    if R < EPS:
        return 0.0

    alph = (2.0 / (1.0 + q)) * (b / R) ** t * R * Omega

    return (x * alph.real + y * alph.imag) / (2.0 - t)

@njit(cache=True, fastmath=True)
def fermat_potential_scalar(z, x, y, x_source, y_source,
                            b, t, gamma1, gamma2, q, phi, Omega):
    """
    Compute the Fermat potential at one image position. 

    Returns the geometric minus gravitational time-delay contribution: 
    ``tau = 0.5*|theta - beta|^2 - psi_EPL(theta) - psi_shear(theta)``

    Parameters
    ----------
    z : ``complex``
        Image position in the axis-aligned frame
        (``z = exp(-i*phi) * (x + i*y)``).
    x : ``float``
        Image x-coordinate in sky frame.
    y : ``float``
        Image y-coordinate in sky frame.
    x_source : ``float``
        Source x-coordinate in sky frame.
    y_source : ``float``
        Source y-coordinate in sky frame.
    b : ``float``
        EPL scale parameter (``b = theta_E * sqrt(q)``).
    t : ``float``
        EPL slope parameter (``t = gamma - 1``).
    gamma1 : ``float``
        External shear component 1.
    gamma2 : ``float``
        External shear component 2.
    q : ``float``
        Axis ratio (0 < q ≤ 1).
    phi : ``float``
        Lens position angle in radians.
    Omega : ``complex``
        EPL deflection kernel at the image elliptical angle.

    Returns
    -------
    tau : ``float``
        Fermat potential (dimensionless; proportional to arrival-time delay).
    """

    # x, y in axis-aligned frame
    x_rot = z.real
    y_rot = z.imag
    # Gravitational time delay 
    # EPL potential 
    tau_grav = (
        _epl_function(x_rot, y_rot, b, t, q, phi, Omega)
        + _shear_function(x, y, gamma1, gamma2)
    )
    # Geometric time delay
    tau_geom = 0.5 * ((x - x_source) ** 2 + (y - y_source) ** 2)

    return tau_geom - tau_grav
# ------------------------------------------------

@njit(cache=True, fastmath=True)
def image_position_analytical_njit(
    x_src,
    y_src,
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
    Standalone EPL + external shear analytical image finder. 

    Locates all lensed images for a given source position, computes
    signed magnifications, Fermat potentials (arrival-time proxies),
    and image types. Results are sorted by ascending arrival time and
    filtered by a minimum magnification threshold.

    Parameters
    ----------
    x_src : ``float``
        Source x-coordinate (normalised to Einstein radius).
    y_src : ``float``
        Source y-coordinate (normalised to Einstein radius).
    q : ``float``
        Lens axis ratio (0 < q ≤ 1).
    phi : ``float``
        Lens position angle in radians.
    gamma : ``float``
        EPL power-law slope (``gamma = 2`` for isothermal).
    gamma1 : ``float``
        External shear component 1.
    gamma2 : ``float``
        External shear component 2.
    theta_E : ``float``
        Einstein radius. 
        default: 1.0
    alpha_scaling : ``float``
        Deflection scaling applied as
        ``theta_E_eff = theta_E * alpha_scaling^(1/(gamma-1))``. 
        default: 1.0
    magnification_limit : ``float``
        Minimum ``|mu|`` for an image to be retained. 
        default: 0.01
    Nmeas : ``int``
        Number of uniformly-spaced angle samples for root finding. 
        default: 400
    Nmeas_extra : ``int``
        Number of extra refinement samples near the source angle. 
        default: 80

    Returns
    -------
    x_img : ``numpy.ndarray``
        Image x-coordinates.
    y_img : ``numpy.ndarray``
        Image y-coordinates.
    fermat_pot : ``numpy.ndarray``
        Fermat potentials at each image.
    magnification : ``numpy.ndarray``
        Signed magnifications at each image.
    image_type : ``numpy.ndarray``
        Image-type codes (1 = min, 2 = saddle, 3 = max, 0 = degenerate).
    nimg : ``int``
        Number of valid images retained.

    Examples
    --------
    >>> from ler.image_properties.epl_shear_njit import image_position_analytical_njit
    >>> x_img, y_img, fermat_pot, mu, itype, nimg = image_position_analytical_njit(
    ...     x_src=0.1, y_src=0.05,
    ...     q=0.8, phi=0.3, gamma=2.0,
    ...     gamma1=0.05, gamma2=0.02,
    ... )
    >>> print(nimg)
    """
    theta_E_eff = theta_E * alpha_scaling ** (1.0 / (gamma - 1.0))

    t = gamma - 1.0
    b = theta_E_eff * np.sqrt(q)

    x_img, y_img, nimg = _solve_lenseq_pemd(
        x_src, y_src, q, phi, t, gamma1, gamma2,
        b=b,
        Nmeas=Nmeas,
        Nmeas_extra=Nmeas_extra,
    )

    fermat_pot = np.empty(MAX_IMGS, dtype=np.float64)
    magnification = np.empty(MAX_IMGS, dtype=np.float64)
    image_type = np.empty(MAX_IMGS, dtype=np.int64)

    # Rotate source position and shear to the axis-aligned (major-axis) frame.
    # fermat_potential_scalar works with image position z = exp(-i*phi)*(x+iy),
    # so the source and shear must also be in the same axis-aligned frame.
    # z_src = np.exp(-1j * phi) * (x_src + 1j * y_src)
    # x_src_rot = z_src.real
    # y_src_rot = z_src.imag
    # g_rot = (gamma1 + 1j * gamma2) * np.exp(-2j * phi)
    # gamma1_rot = g_rot.real
    # gamma2_rot = g_rot.imag

    # precompute rotation factors once (reused for every image)
    cp_phi = np.cos(phi)
    sp_phi = np.sin(phi)

    for i in range(nimg):
        # z is in axis-aligned frame: z = exp(-i*phi)*(x+iy)
        xi = x_img[i]
        yi = y_img[i]
        z_re =  xi * cp_phi + yi * sp_phi
        z_im = -xi * sp_phi + yi * cp_phi
        z = z_re + 1j * z_im
        # angle in the elliptical frame
        phi_ell = np.angle(z_re * q + 1j * z_im)
        # Omega is calculated only once per image, and reused for both magnification and Fermat potential calculations.
        Omega = omega_scalar(phi_ell, t, q)
        mu, itype = lensing_diagnostics_scalar(
            z, b, t, gamma1, gamma2, q, phi, Omega
        )

        magnification[i] = mu
        image_type[i] = itype

        fermat_pot[i] = fermat_potential_scalar(
            z, xi, yi, x_src, y_src,
            b, t, gamma1, gamma2, q, phi, Omega
        )

    # sort by arrival time (Fermat potential)
    _insertion_sort5(x_img, y_img, fermat_pot, magnification, image_type, nimg)

    # inline compress: keep images above magnification_limit (no keep[] array allocation)
    m = 0
    for i in range(nimg):
        if not np.isfinite(magnification[i]) or abs(magnification[i]) >= magnification_limit:
            x_img[m] = x_img[i]
            y_img[m] = y_img[i]
            fermat_pot[m] = fermat_pot[i]
            magnification[m] = magnification[i]
            image_type[m] = image_type[i]
            m += 1
    nimg = m

    return (
        x_img[:nimg],
        y_img[:nimg],
        fermat_pot[:nimg],
        magnification[:nimg],
        image_type[:nimg],
        nimg,
    )

def create_epl_shear_solver(
    arrival_time_sort=True,
    max_img=4,
    num_th=500,
    maginf=-100.0,
    max_source_sample_attempts=128,
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
    max_source_sample_attempts : ``int``
        Maximum draws from ``sample_source_from_double_caustic`` per system
        when the returned ``(beta_x, beta_y)`` are non-finite (invalid
        caustic / geometry). Stops early on the first finite pair. \n
        default: 128
    alpha_scaling : ``float``
        Deflection scaling factor. \n
        default: 1.0
    magnification_limit : ``float``
        Minimum ``abs(mu)`` threshold for image retention. \n
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

    @njit(parallel=True, cache=False, fastmath=True)
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
            Source x-coordinates (``NaN`` if sampling stayed invalid after
            ``max_source_sample_attempts`` tries).
        beta_y_arr : ``numpy.ndarray``
            Source y-coordinates (same convention as ``beta_x_arr``).
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

        # find source position for each lens system and sample caustic points
        for i in prange(size):
            # Retry until finite coordinates or attempts exhausted (NaNs stay NaN).
            beta_x = np.nan
            beta_y = np.nan
            for _ in range(max_source_sample_attempts):
                bx, by = sample_source_from_double_caustic(
                    theta_E=1.0,
                    q=q[i],
                    phi=phi[i],
                    gamma=gamma[i],
                    gamma1=gamma1[i],
                    gamma2=gamma2[i],
                    num_th=num_th,
                    maginf=maginf,
                )
                if np.isfinite(bx) and np.isfinite(by):
                    beta_x = bx
                    beta_y = by
                    break

            beta_x_arr[i] = beta_x
            beta_y_arr[i] = beta_y

        # solve for image positions in systems with valid source samples
        idx = np.where(np.isfinite(beta_x_arr) & np.isfinite(beta_y_arr))[0]

        for i in prange(idx.size):

            idx_i = idx[i]

            theta_E_val = theta_E[idx_i]

            x, y, fermat_pot, mu, image_type, n = image_position_analytical_njit(
                x_src=beta_x_arr[idx_i], 
                y_src=beta_y_arr[idx_i], 
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

            n = min(n, max_img)  # enforce max_img limit

            nimg[idx_i] = n

            if n == 0:
                # no valid images; result arrays already NaN/0
                beta_x_arr[idx_i] = beta_x_arr[idx_i] * theta_E_val
                beta_y_arr[idx_i] = beta_y_arr[idx_i] * theta_E_val
                continue

            mu_arr[idx_i, :n] = mu[:n]

            # rescale to physical units
            x_img[idx_i, :n] = x[:n] * theta_E_val
            y_img[idx_i, :n] = y[:n] * theta_E_val
            beta_x_arr[idx_i] = beta_x_arr[idx_i] * theta_E_val
            beta_y_arr[idx_i] = beta_y_arr[idx_i] * theta_E_val

            # rescale time delays
            tau_hat = fermat_pot[:n]
            # Fermat potential rescales as theta_E^2
            tau_phys = tau_hat * (theta_E_val * theta_E_val)
            # physical time delays in days: Δt = (D_dt/c) * Δtau
            tau_ = (D_dt[idx_i] / C_LIGHT) * tau_phys
            # positive time delays 
            tau_arr[idx_i, :n] = tau_ - np.min(tau_)

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
