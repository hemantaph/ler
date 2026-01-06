from numba import njit, prange
import numpy as np
C_LIGHT = 299792.458  # km/s

# ----------------------
# Rejection sampler
# ----------------------
def rejection_sampler(
    zs, 
    zl,
    sigma_max,
    sigma_rvs,
    q_rvs,
    phi_rvs,
    gamma_rvs,
    shear_rvs,
    cross_section,
    saftey_factor = 1.2
):
    size = zl.size

    sigma_array = np.zeros(size)
    q_array = np.zeros(size)
    phi_array = np.zeros(size)
    gamma_array = np.zeros(size)
    gamma1_array = np.zeros(size)
    gamma2_array = np.zeros(size)
    idx_remaining = np.arange(size)

    # Compute the maximum cross section
    # to be used as the upper bound for the rejection sampler
    cs_max = cross_section(
                zs=zs,
                zl=zl,
                sigma=sigma_max*np.ones(size),
                q=0.9*np.ones(size),
                phi=np.zeros(size),
                gamma=2.645*np.ones(size),
                gamma1=np.zeros(size),
                gamma2=np.zeros(size)
              ) * saftey_factor

    while len(idx_remaining) > 0:
        n_remaining = len(idx_remaining)
        sigma_samples = sigma_rvs(n_remaining, zl[idx_remaining])
        q_samples = q_rvs(n_remaining, sigma_samples)
        phi_samples = phi_rvs(n_remaining)
        gamma_samples = gamma_rvs(n_remaining)
        gamma1_samples, gamma2_samples = shear_rvs(n_remaining)

        cs = cross_section(
            zs=zs[idx_remaining],
            zl=zl[idx_remaining],
            sigma=sigma_samples,
            q=q_samples,
            phi=phi_samples,
            gamma=gamma_samples,
            gamma1=gamma1_samples,
            gamma2=gamma2_samples
        )

        accept = np.random.random(n_remaining) < (cs / cs_max[idx_remaining])

        accepted_indices = idx_remaining[accept]
        sigma_array[accepted_indices] = sigma_samples[accept]
        q_array[accepted_indices] = q_samples[accept]
        phi_array[accepted_indices] = phi_samples[accept]
        gamma_array[accepted_indices] = gamma_samples[accept]
        gamma1_array[accepted_indices] = gamma1_samples[accept]
        gamma2_array[accepted_indices] = gamma2_samples[accept]

        idx_remaining = idx_remaining[~accept]

    return sigma_array, q_array, phi_array, gamma_array, gamma1_array, gamma2_array

def create_rejection_sampler(
    sigma_max,
    sigma_rvs,
    q_rvs,
    phi_rvs,
    gamma_rvs,
    shear_rvs,
    cross_section,
    saftey_factor=1.2,
    use_njit_sampler=True
):
    
    if use_njit_sampler:
        _base_sampler = njit(rejection_sampler)
        print("Faster, njitted and rejection sampling based lens parameter sampler will be used.")

        @njit
        def rejection_sampler_wrapper(
            zs, zl,
        ):
            return _base_sampler(
                zs=zs,
                zl=zl,
                sigma_max=sigma_max,
                sigma_rvs=sigma_rvs,
                q_rvs=q_rvs,
                phi_rvs=phi_rvs,
                gamma_rvs=gamma_rvs,
                shear_rvs=shear_rvs,
                cross_section=cross_section,
                saftey_factor=saftey_factor
            )
    else:
        print("Slower, non-njit and rejection sampling based lens parameter sampler will be used.")   
        def rejection_sampler_wrapper(
            zs, zl,
        ):
            return rejection_sampler(
                zs=zs,
                zl=zl,
                sigma_max=sigma_max,
                sigma_rvs=sigma_rvs,
                q_rvs=q_rvs,
                phi_rvs=phi_rvs,
                gamma_rvs=gamma_rvs,
                shear_rvs=shear_rvs,
                cross_section=cross_section,
                saftey_factor=saftey_factor
            )

    return rejection_sampler_wrapper


# ----------------------
# Importance sampler
# ----------------------
@njit
def sigma_proposal_uniform(n, sigma_min, sigma_max):
    return sigma_min + (sigma_max - sigma_min) * np.random.random(n)


@njit
def weighted_choice_1d(weights):
    """
    Draw an index with probability proportional to 'weights' (assumed >=0).
    Numba-safe replacement for np.random.choice(n, p=weights).
    """
    total = 0.0
    for i in range(weights.size):
        total += weights[i]

    # If all weights are zero (or NaN-ish), fall back to uniform.
    if not (total > 0.0):
        return np.random.randint(weights.size)

    u = np.random.random() * total
    c = 0.0
    for i in range(weights.size):
        c += weights[i]
        if u <= c:
            return i
    return weights.size - 1  # numerical fallback


# ----------------------
# Importance sampler
# ----------------------
def importance_sampler(
    zs,
    zl,
    sigma_min,
    sigma_max,
    q_rvs,
    phi_rvs,
    gamma_rvs,
    shear_rvs,
    sigma_pdf,
    cross_section,
    n_prop,
):

    N = zl.size

    sigma_post = np.zeros(N)
    q_post = np.zeros(N)
    phi_post = np.zeros(N)
    gamma_post = np.zeros(N)
    gamma1_post = np.zeros(N)
    gamma2_post = np.zeros(N)

    p0 = 1.0 / (sigma_max - sigma_min)

    for i in prange(N):
        # (1) draw proposals
        sigma_prop = sigma_proposal_uniform(n_prop, sigma_min, sigma_max)

        # draw other parameters from their priors
        q_prop = q_rvs(n_prop, sigma_prop)
        phi_prop = phi_rvs(n_prop)
        gamma_prop = gamma_rvs(n_prop)
        gamma1_prop, gamma2_prop = shear_rvs(n_prop)

        # (2) compute cross sections
        # cross_section expects (zs, zl, sigma, q, phi, gamma, gamma1, gamma2)
        zs_arr = zs[i] * np.ones(n_prop)
        zl_arr = zl[i] * np.ones(n_prop)

        cs = cross_section(
            zs_arr, zl_arr, sigma_prop, q_prop, phi_prop, gamma_prop, gamma1_prop, gamma2_prop
        )

        # Stabilize: normalize cs to avoid overflow in weights
        cs_sum = 0.0
        for k in range(cs.size):
            cs_sum += cs[k]

        if cs_sum > 0.0:
            for k in range(cs.size):
                cs[k] = cs[k] / cs_sum
        else:
            # if cs all zeros, sampling should revert to prior/proposal only
            # which is p_sigma/p0; keep cs = 1 so w ∝ p_sigma/p0
            for k in range(cs.size):
                cs[k] = 1.0

        # (3) compute weights: w = cs * (p_sigma / p0)
        zl_vec = zl[i] * np.ones(n_prop)
        p_sigma = sigma_pdf(sigma_prop, zl_vec)

        w = np.empty(n_prop)
        w_sum = 0.0
        for k in range(n_prop):
            wk = cs[k] * (p_sigma[k] / p0)
            # guard against negatives / NaNs
            if not (wk > 0.0):
                wk = 0.0
            w[k] = wk
            w_sum += wk

        # normalize weights
        if w_sum > 0.0:
            for k in range(n_prop):
                w[k] /= w_sum
        else:
            # fully degenerate case -> uniform
            inv = 1.0 / n_prop
            for k in range(n_prop):
                w[k] = inv

        # (4) draw 1 posterior sample via weighted choice
        idx = weighted_choice_1d(w)

        sigma_post[i] = sigma_prop[idx]
        q_post[i] = q_prop[idx]
        phi_post[i] = phi_prop[idx]
        gamma_post[i] = gamma_prop[idx]
        gamma1_post[i] = gamma1_prop[idx]
        gamma2_post[i] = gamma2_prop[idx]

    return sigma_post, q_post, phi_post, gamma_post, gamma1_post, gamma2_post


def create_importance_sampler(
    sigma_min,
    sigma_max,
    q_rvs,
    phi_rvs,
    gamma_rvs,
    shear_rvs,
    sigma_pdf,
    cross_section,
    n_prop,
    use_njit_sampler=True,
):
    if use_njit_sampler:
        print("Faster, njitted and importance sampling based lens parameter sampler will be used.")
        _base_sampler = njit(parallel=True)(importance_sampler)

        @njit(parallel=True)
        def importance_sampler_wrapper(
            zs, zl,
        ):
            return _base_sampler(
                zs=zs,
                zl=zl,
                sigma_min=sigma_min,
                sigma_max=sigma_max,
                q_rvs=q_rvs,
                phi_rvs=phi_rvs,
                gamma_rvs=gamma_rvs,
                shear_rvs=shear_rvs,
                sigma_pdf=sigma_pdf,
                cross_section=cross_section,
                n_prop=n_prop,
            )
    else:
        print("Slower, non-njit and importance sampling based lens parameter sampler will be used.")
        def importance_sampler_wrapper(
            zs, zl,
        ):
            return importance_sampler(
                zs=zs,
                zl=zl,
                sigma_min=sigma_min,
                sigma_max=sigma_max,
                q_rvs=q_rvs,
                phi_rvs=phi_rvs,
                gamma_rvs=gamma_rvs,
                shear_rvs=shear_rvs,
                sigma_pdf=sigma_pdf,
                cross_section=cross_section,
                n_prop=n_prop,
            )

    return importance_sampler_wrapper
