import numpy as np
from numba import njit, jit
from scipy.stats import rayleigh

# def rjs_with_einstein_radius(self, param_dict):
#     theta_E = param_dict["theta_E"]
#     size = len(theta_E)
#     theta_E_max = np.max(theta_E)  # maximum einstein radius
#     u = np.random.uniform(0, theta_E_max**2, size=size)
#     mask = u < theta_E**2

#     # return the dictionary with the mask applied
#     return {key: val[mask] for key, val in param_dict.items()}

@njit
def rjs_with_cross_section(theta_E, q):
    pass



@njit
def axis_ratio_SIS(sigma):
    return np.ones(len(sigma))

# sampler for the lens redshift
@njit
def lens_redshift_SDSS_catalogue(zs):
    """
    Function to sample lens redshift from the SDSS catalogue.
    """

    size = zs.size
    r = np.random.uniform(0, 1, size=size)
    # sample from inverse of cdf of the lens redshift distribution
    # See the integral of Eq. A7 of https://arxiv.org/pdf/1807.07062.pdf (cdf)
    return 10 * r**3 - 15 * r**4 + 6 * r**5
    
@njit
def gamma_(x):
    # Coefficients for the Lanczos approximation
    g = 7
    p = np.array([0.99999999999980993, 676.5203681218851, -1259.1392167224028,
                  771.32342877765313, -176.61502916214059, 12.507343278686905,
                  -0.13857109526572012, 9.9843695780195716e-6, 1.5056327351493116e-7])

    if x < 0.5:
        # Reflection formula
        return np.pi / (np.sin(np.pi * x) * gamma_(1 - x))
    else:
        x -= 1
        y = p[0]
        for i in range(1, g + 2):
            y += p[i] / (x + i)
        t = x + g + 0.5
        return np.sqrt(2 * np.pi) * t**(x + 0.5) * np.exp(-t) * y
    
@njit
def cvdf_fit(log_vd, redshift):
    this_vars = np.array([
        [7.39149763, 5.72940031, -1.12055245],
        [-6.86339338, -5.27327109, 1.10411386],
        [2.85208259, 1.25569600, -0.28663846],
        [0.06703215, -0.04868317, 0.00764841]])
    coeffs = [this_vars[i][0] + this_vars[i][1] * redshift + this_vars[i][2] * redshift ** 2 for i in range(4)]
    mstar = log_vd - coeffs[3]
    return coeffs[0] + coeffs[1] * mstar + coeffs[2] * mstar ** 2 - np.exp(mstar)

@njit
def my_derivative(log_vd, redshift, dx):
    return 0.5 * (cvdf_fit(log_vd + dx, redshift) - cvdf_fit(log_vd - dx, redshift)) / dx

@njit
def pdf_phi_z_div_0(s, z):
    log_vd = np.log10(s)
    phi_sim_z = 10 ** cvdf_fit(log_vd, z) / s * my_derivative(log_vd, z, 1e-8)
    phi_sim_0 = 10 ** cvdf_fit(log_vd, 0) / s * my_derivative(log_vd, 0, 1e-8)

    return phi_sim_z / phi_sim_0

@njit
def phi(s,z):
    result = s**4*pdf_phi_z_div_0(s,z)*phi_loc_bernardi(s)
    result[result < 0.] = 0.
    return result

# elliptical lens galaxy
@njit
def phi_cut_SIE(q):
    n = len(q)
    result = np.empty(n)
    for i in range(n):
        val = q[i]
        if 0.01 < val < 0.99:
            result[i] = (2 * np.pi * val * np.log(val)) / (val ** 2 - 1)
        elif val < 0.01:
            result[i] = -2 * (np.pi * np.log(val)) * val
        else:
            result[i] = np.pi
    return result/np.pi

@jit
def axis_ratio_rayleigh(sigma, q_min=0.2, param=None):
        """
        Function to sample axis ratio from rayleigh distribution with given velocity dispersion.

        Parameters
        ----------
        sigma : `float: array`
            velocity dispersion of the lens galaxy

        Returns
        -------
        q : `float: array`
            axis ratio of the lens galaxy
        """

        if param:
            q_min = param["q_min"]

        size = len(sigma)
        a = sigma / 161.0
        q = np.ones(size)
        idx = np.arange(size)  # idx tracker
        size_ = size

        while size_ != 0:
            # Draw the axis ratio see Appendix of https://arxiv.org/pdf/1807.07062.pdf
            s = abs(0.38 - 0.09177 * a[idx])
            b = rayleigh.rvs(scale=s, size=size_)
            q_ = 1.0 - b

            # Weed out axis ratios that have axis ratio below q_min
            idx2 = q_ > q_min
            q[idx[idx2]] = q_[idx2]

            # remaining idx from the original array
            # that still not have axis ratio above q_min
            idx = idx[q <= q_min]
            size_ = len(idx)

        return q