import numpy as np
from numba import njit, prange

# Import only the truly-scalar helpers (rotation + shear conversion).
# Do NOT import vectorized alpha functions.
from .cross_section_njit import _shear_c2p, _shear_p2c, _rot_apply, _omega
from .sample_caustic_points_njit import sample_source_from_double_caustic

C_LIGHT = 299792458.0  # m/s


# ---------------------------------------------------------------------------
# Scalar deflections: EPL and shear (NO tiny array allocations)
# ---------------------------------------------------------------------------


@njit(cache=True)
def _alpha_epl_xy_s(x, y, theta_E, q, gamma):
    """
    Scalar EPL deflection in solver frame.
    Matches the scalar version used in your earlier cross-section implementation.
    """
    t = gamma - 1.0
    b = theta_E * np.sqrt(q)

    R = np.sqrt((q * x) * (q * x) + y * y)
    if R <= 1e-15:
        return 0.0, 0.0

    pf = (2.0 * b) / (1.0 + q)
    fr = (b / R) ** t * (R / b)
    ph = np.arctan2(y, x * q)
    Om = _omega(ph, t, q)
    a = pf * fr * Om
    return a.real, a.imag


@njit(cache=True)
def _alpha_shear_xy_s(x, y, g1, g2):
    """Scalar shear deflection."""
    return g1 * x + g2 * y, g2 * x - g1 * y


@njit(cache=True)
def _alpha_total_xy_s(x, y, theta_E, q, gamma, g1, g2):
    """Scalar total deflection."""
    ax, ay = _alpha_epl_xy_s(x, y, theta_E, q, gamma)
    sx, sy = _alpha_shear_xy_s(x, y, g1, g2)
    return ax + sx, ay + sy


# ---------------------------------------------------------------------------
# Private helpers — potential, Fermat (scalar)
# ---------------------------------------------------------------------------


@njit(cache=True)
def _psi_shear_xy(x, y, g1, g2):
    return 0.5 * (g1 * (x * x - y * y) + 2.0 * g2 * x * y)


@njit(cache=True)
def _psi_epl_via_euler_xy(x, y, theta_E, q, gamma):
    """
    EPL potential via Euler relation: psi = (theta · alpha_EPL)/(3-gamma)
    """
    t = gamma - 1.0
    denom = 2.0 - t  # 3-gamma
    ax, ay = _alpha_epl_xy_s(x, y, theta_E, q, gamma)
    return (x * ax + y * ay) / denom


@njit(cache=True)
def _fermat_tau_xy(x, y, bx, by, theta_E, q, gamma, g1, g2):
    geom = 0.5 * ((x - bx) * (x - bx) + (y - by) * (y - by))
    psi = _psi_epl_via_euler_xy(x, y, theta_E, q, gamma) + _psi_shear_xy(x, y, g1, g2)
    return geom - psi


# ---------------------------------------------------------------------------
# Jacobian, image type, Newton refinement (scalar)
# ---------------------------------------------------------------------------


@njit(cache=True)
def _jacobian_A_fd(x, y, theta_E, q, gamma, g1, g2, h):
    """
    Jacobian A = I - d(alpha)/d(theta) via forward finite differences.
    Uses scalar deflections (fast).
    """
    a0x, a0y = _alpha_total_xy_s(x, y, theta_E, q, gamma, g1, g2)

    axx, axy = _alpha_total_xy_s(x + h, y, theta_E, q, gamma, g1, g2)
    ayx, ayy = _alpha_total_xy_s(x, y + h, theta_E, q, gamma, g1, g2)

    dalpha_dx_x = (axx - a0x) / h
    dalpha_dx_y = (axy - a0y) / h
    dalpha_dy_x = (ayx - a0x) / h
    dalpha_dy_y = (ayy - a0y) / h

    A00 = 1.0 - dalpha_dx_x
    A10 = -dalpha_dx_y
    A01 = -dalpha_dy_x
    A11 = 1.0 - dalpha_dy_y
    return A00, A01, A10, A11


@njit(cache=True)
def _image_type_from_A(A00, A01, A10, A11):
    det = A00 * A11 - A01 * A10
    tr = A00 + A11
    if det < 0.0:
        return 2
    if tr > 0.0:
        return 1
    return 3


@njit(cache=True)
def _newton_refine(x0, y0, bx, by, theta_E, q, gamma, g1, g2, h):
    """
    Newton refinement for f(theta)=theta-alpha(theta)-beta=0.
    Scalar throughout.
    """
    x = x0
    y = y0

    for _ in range(30):  # 40 -> 30 usually enough if seeds are reasonable
        ax, ay = _alpha_total_xy_s(x, y, theta_E, q, gamma, g1, g2)
        fx = x - ax - bx
        fy = y - ay - by

        if fx * fx + fy * fy < 1e-20:
            return x, y, 1

        A00, A01, A10, A11 = _jacobian_A_fd(x, y, theta_E, q, gamma, g1, g2, h)
        detA = A00 * A11 - A01 * A10
        if abs(detA) < 1e-18:
            return x, y, 0

        dx = (A11 * fx - A01 * fy) / detA
        dy = (-A10 * fx + A00 * fy) / detA

        x -= dx
        y -= dy

        # optional damping if divergence is an issue (cheap safeguard)
        # if dx*dx + dy*dy > 1.0:
        #     x += 0.5 * dx
        #     y += 0.5 * dy

    ax, ay = _alpha_total_xy_s(x, y, theta_E, q, gamma, g1, g2)
    fx = x - ax - bx
    fy = y - ay - by
    if fx * fx + fy * fy < 1e-12:
        return x, y, 1
    return x, y, 0


# ---------------------------------------------------------------------------
# Public — solver
# ---------------------------------------------------------------------------


@njit(cache=True)
def solve_epl_shear(
    beta_x,
    beta_y,
    theta_E,
    q,
    phi,
    gamma,
    gamma1,
    gamma2,
    arrival_time_sort=True,
    ngrid=160,
    rng=3.5,
):
    """
    Solve lens equation for EPL+shear.
    Keeps your overall algorithm, but removes pathological allocations.
    """
    # rotate source into solver frame
    bx_s, by_s = _rot_apply(phi, beta_x, beta_y)

    # scale to unit Einstein radius
    if theta_E == 0.0:
        # invalid; return empty
        return (
            0,
            np.empty(4),
            np.empty(4),
            np.empty(4),
            np.empty(4),
            np.empty(4),
            np.empty(4, dtype=np.int64),
        )
    bx_unit = bx_s / theta_E
    by_unit = by_s / theta_E
    theta_E_unit = 1.0

    # rotate shear into solver frame
    phi_g, gmag = _shear_c2p(gamma1, gamma2)
    phi_g_rot = phi_g - phi
    g1_s, g2_s = _shear_p2c(phi_g_rot, gmag)

    # grid parameters
    xs = np.linspace(-rng, rng, ngrid)
    dxg = xs[1] - xs[0]
    thr2 = (dxg * 2.5) * (dxg * 2.5)

    # candidate list (avoid np.where; fill arrays directly)
    max_cand = ngrid * ngrid
    cand_x = np.empty(max_cand, dtype=np.float64)
    cand_y = np.empty(max_cand, dtype=np.float64)
    nc = 0

    # seed candidates by scanning the grid (no big temporary arrays)
    for iy in range(ngrid):
        y0 = xs[iy]
        for ix in range(ngrid):
            x0 = xs[ix]
            ax, ay = _alpha_total_xy_s(x0, y0, theta_E_unit, q, gamma, g1_s, g2_s)
            fx = x0 - ax - bx_unit
            fy = y0 - ay - by_unit
            if fx * fx + fy * fy < thr2:
                cand_x[nc] = x0
                cand_y[nc] = y0
                nc += 1

    # refine candidates with Newton and de-duplicate
    max_img = 4
    img_x_unit = np.empty(max_img, dtype=np.float64)
    img_y_unit = np.empty(max_img, dtype=np.float64)
    nimg = 0

    h = 1e-5

    for k in range(nc):
        x0 = cand_x[k]
        y0 = cand_y[k]
        x, y, ok = _newton_refine(
            x0, y0, bx_unit, by_unit, theta_E_unit, q, gamma, g1_s, g2_s, h
        )
        if ok == 0:
            continue

        # de-duplicate
        dup = 0
        for m in range(nimg):
            dxm = img_x_unit[m] - x
            dym = img_y_unit[m] - y
            if dxm * dxm + dym * dym < 1e-8:
                dup = 1
                break
        if dup == 1:
            continue

        if nimg < max_img:
            img_x_unit[nimg] = x
            img_y_unit[nimg] = y
            nimg += 1
            if nimg == max_img:
                # cannot store more anyway; keep going is wasted work
                # (optional: break if you trust max 4 always)
                pass

    # diagnostics and rotate back to sky frame
    x_img = np.empty(max_img, dtype=np.float64)
    y_img = np.empty(max_img, dtype=np.float64)
    mu_arr = np.empty(max_img, dtype=np.float64)
    detA_arr = np.empty(max_img, dtype=np.float64)
    tau_arr = np.empty(max_img, dtype=np.float64)
    itype = np.empty(max_img, dtype=np.int64)

    theta_E_sq = theta_E * theta_E

    for m in range(nimg):
        x_unit = img_x_unit[m]
        y_unit = img_y_unit[m]

        A00, A01, A10, A11 = _jacobian_A_fd(
            x_unit, y_unit, theta_E_unit, q, gamma, g1_s, g2_s, h
        )
        detA = A00 * A11 - A01 * A10

        mu = np.inf
        if abs(detA) > 1e-18:
            mu = 1.0 / detA

        tau_unit = _fermat_tau_xy(
            x_unit, y_unit, bx_unit, by_unit, theta_E_unit, q, gamma, g1_s, g2_s
        )
        tau = tau_unit * theta_E_sq
        tp = _image_type_from_A(A00, A01, A10, A11)

        x_s = x_unit * theta_E
        y_s = y_unit * theta_E
        xk, yk = _rot_apply(-phi, x_s, y_s)

        x_img[m] = xk
        y_img[m] = yk
        mu_arr[m] = mu
        detA_arr[m] = detA
        tau_arr[m] = tau
        itype[m] = tp

    # fill unused slots
    for m in range(nimg, max_img):
        x_img[m] = np.nan
        y_img[m] = np.nan
        mu_arr[m] = np.nan
        detA_arr[m] = np.nan
        tau_arr[m] = np.nan
        itype[m] = 0

    # sort by arrival time (selection sort; small n)
    if arrival_time_sort and nimg > 1:
        for i in range(nimg - 1):
            kmin = i
            tmin = tau_arr[i]
            for j in range(i + 1, nimg):
                if tau_arr[j] < tmin:
                    tmin = tau_arr[j]
                    kmin = j
            if kmin != i:
                tau_arr[i], tau_arr[kmin] = tau_arr[kmin], tau_arr[i]
                x_img[i], x_img[kmin] = x_img[kmin], x_img[i]
                y_img[i], y_img[kmin] = y_img[kmin], y_img[i]
                mu_arr[i], mu_arr[kmin] = mu_arr[kmin], mu_arr[i]
                detA_arr[i], detA_arr[kmin] = detA_arr[kmin], detA_arr[i]
                itype[i], itype[kmin] = itype[kmin], itype[i]

    return nimg, x_img, y_img, mu_arr, detA_arr, tau_arr, itype


def create_epl_shear_solver(
    arrival_time_sort=True,
    ngrid=160,
    rng=3.5,
    max_img=4,
    num_th=500,
    maginf=-100.0,
    max_tries=100000,
):
    """
    Same public API as your current factory, but uses the optimized scalar solve.
    """

    @njit(parallel=True, cache=True)
    def solve_epl_shear_multithreaded(theta_E, D_dt, q, phi, gamma, gamma1, gamma2):
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

        for i in prange(size):

            beta_x, beta_y, ok = sample_source_from_double_caustic(
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

            if ok == 1:

                theta_E_val = theta_E[i]
                n, x, y, mu, dA, tau, tp = solve_epl_shear(
                    beta_x,
                    beta_y,
                    1.0,
                    q[i],
                    phi[i],
                    gamma[i],
                    gamma1[i],
                    gamma2[i],
                    arrival_time_sort=arrival_time_sort,
                    ngrid=ngrid,
                    rng=rng,
                )

                nimg[i] = n
                mu_arr[i, :n] = mu[:n]

                # rescale to physical units
                x_img[i, :n] = x[:n] * theta_E_val
                y_img[i, :n] = y[:n] * theta_E_val
                beta_x_arr[i] = beta_x * theta_E_val
                beta_y_arr[i] = beta_y * theta_E_val

                # rescale time delays
                tau_hat = tau[:n]
                # Fermat potential rescales as theta_E^2
                tau_phys = tau_hat * (theta_E_val * theta_E_val)
                # physical time delays in days: Δt = (D_dt/c) * Δtau
                tau_arr[i, :n] = (D_dt[i] / C_LIGHT) * tau_phys

                # image type
                itype[i, :n] = tp[:n]

            else:
                pass

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
