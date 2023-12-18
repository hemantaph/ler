import numpy as np
from numba import njit
#from scipy.stats import rayleigh
from ..utils import inverse_transform_sampler, cubic_spline_interpolator

from lenstronomy.LensModel.Solver.epl_shear_solver import caustics_epl_shear
from shapely.geometry import Polygon

# def rjs_with_einstein_radius(self, param_dict):
#     theta_E = param_dict["theta_E"]
#     size = len(theta_E)
#     theta_E_max = np.max(theta_E)  # maximum einstein radius
#     u = np.random.uniform(0, theta_E_max**2, size=size)
#     mask = u < theta_E**2

#     # return the dictionary with the mask applied
#     return {key: val[mask] for key, val in param_dict.items()}



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
def phi(s,z, cosmology_h=0.7):

    result = s**4*pdf_phi_z_div_0(s,z)*phi_loc_bernardi(sigma=s, cosmology_h=cosmology_h)
    # result[result < 0.] = 0.
    return result

@njit
def phi_loc_bernardi(sigma, alpha=0.94, beta=1.85, phistar=2.099e-2, sigmastar=113.78, cosmology_h=0.7):
    phistar = phistar * (cosmology_h / 0.7) ** 3  # Mpc**-3
    philoc_ = phistar*(sigma/sigmastar)**alpha * np.exp(-(sigma/sigmastar)**beta) * beta/gamma_(alpha/beta)/sigma
    return philoc_

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

@njit
def axis_ratio_rayleigh(sigma, q_min=0.2, q_max=1.0):
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

        size = len(sigma)
        a = sigma / 161.0
        q = np.ones(size)
        idx = np.arange(size)  # idx tracker
        size_ = size

        while size_ != 0:
            # Draw the axis ratio see Appendix of https://arxiv.org/pdf/1807.07062.pdf
            s = 0.38 - 0.09177 * a[idx]
            s[s<=0] = 0.0001
            u = np.random.uniform(0, 1, size=size_)
            b = s * np.sqrt(-2 * np.log(u))  # inverse cdf rayleigh distribution
            q_ = 1.0 - b

            # selecting only q that is within the range
            idx2 = (q_ >= q_min) & (q_ <= q_max)
            q[idx[idx2]] = q_[idx2]

            # remaining idx from the original array
            # that still not have axis ratio above q_min
            idx = idx[(q_ <= q_min) | (q_ >= q_max)]

            size_ = len(idx)

        return q

@njit
def velocity_dispersion_z_dependent(size, zl, zl_list, vd_inv_cdf):
    """
    Function to sample velocity dispersion from the interpolator

    Parameters
    ----------
    size: int
        Number of samples to draw
    zl: `numpy.ndarray` (1D array of float of size=size)
        Redshift of the lens galaxy

    Returns
    ----------
    samples: numpy.ndarray
        Samples of velocity dispersion
    """

    index = np.searchsorted(zl_list, zl)
    u = np.random.uniform(0, 1, size)
    samples = np.zeros(size)
        
    for i in range(size):
        cdf, x = vd_inv_cdf[index[i],0], vd_inv_cdf[index[i],1]
        idx = np.searchsorted(cdf, u[i])  # vd cdf
        x1, x0, y1, y0 = cdf[idx], cdf[idx-1], x[idx], x[idx-1]
        samples[i] = y0 + (y1 - y0) * (u[i] - x0) / (x1 - x0)

    return samples

@njit
def lens_redshift_SDSS_catalogue(zs, splineDc, splineDcInv, u, cdf):
    """
    Function to sample lens redshift from the SDSS catalogue.

    Parameters
    ----------
    zs: `numpy.ndarray` (1D array of float of size=size)
        Redshift of the source galaxy
    splineDc: `list`
        List of spline coefficients for the comoving distance and redshifts
    splineDcInv: `list`
        List of spline coefficients for the inverse of comoving distance and redshifts
    u: `numpy.ndarray` (1D array of float of size=size)
        e.g. u = np.linspace(0, 1, 500)
    cdf: `numpy.ndarray` (1D array of float of size=size)
        Cumulative distribution function of the lens redshift distribution between 0 and 1

    Returns
    ----------
    zl: `numpy.ndarray` (1D array of float of size=size)
        Redshift of the lens galaxy corresponding to the zs
    """

    splineDc_coeff = splineDc[0]
    splineDc_z_list = splineDc[1]
    splineDcInv_coeff = splineDcInv[0]
    splineDcInv_z_list = splineDcInv[1]

    size = len(zs)
    r = inverse_transform_sampler(size, cdf, u)
    lens_galaxy_Dc = cubic_spline_interpolator(zs, splineDc_coeff, splineDc_z_list) * r  # corresponding element-wise multiplication between 2 arrays

    return cubic_spline_interpolator(lens_galaxy_Dc, splineDcInv_coeff, splineDcInv_z_list)


def phi_cut_epl_shear(e1, e2, gamma, gamma1, gamma2):

    results = []
    size = len(e1)  # size of the array
    for i in range(size):

        while True:
            kwargs_lens = [
                {
                    "theta_E": 1.,
                    "e1": e1[i],
                    "e2": e2[i],
                    "gamma": gamma[i],
                    "center_x": 0.0,
                    "center_y": 0.0,
                },
                {
                    "gamma1": gamma1[i],
                    "gamma2": gamma2[i],
                    "ra_0": 0,
                    "dec_0": 0,
                },
            ]
            caustic_double_points = caustics_epl_shear(
                kwargs_lens, return_which="double", maginf=-5000
            )
            caustic = np.logical_not(np.isnan(caustic_double_points).any())
            # If there is a nan, caustic=False, draw a new gamma
            if caustic:
                break
            else:
                gamma[i] = np.random.normal(loc=2.0, scale=0.2, size=1)[0]

        caustic_double = Polygon(caustic_double_points.T)
        area = caustic_double.area

        # kwargs_lens = [
        #     {
        #         "theta_E": 1.,
        #         "e1": e1[i],
        #         "e2": e2[i],
        #         "gamma": gamma[i],
        #         "center_x": 0.0,
        #         "center_y": 0.0,
        #     },
        #     {
        #         "gamma1": gamma1[i],
        #         "gamma2": gamma2[i],
        #         "ra_0": 0,
        #         "dec_0": 0,
        #     },
        # ]

        # caustic_double_points = caustics_epl_shear(
        #             kwargs_lens, return_which="double", maginf=-5000
        #         )

        # try:
        #     caustic_double = Polygon(caustic_double_points.T)
        
        #     area = caustic_double.area
        # except:
        #     print("\n",e1[i], e2[i], gamma[i], gamma1[i], gamma2[i])

        results.append(np.round(area/np.pi, 3))

    return np.array(results)

@njit
def bounded_normal_sample(size, mean, std, low, high):
    """
    Function to sample from a normal distribution with bounds.

    Parameters
    ----------
    mean: `float`
        Mean of the normal distribution
    std: `float`
        Standard deviation of the normal distribution
    low: `float`
        Lower bound
    high: `float`
        Upper bound
    """
    samples = np.empty(size)
    for i in range(size):
        while True:
            sample = np.random.normal(mean, std)
            if low <= sample <= high:
                break
        samples[i] = sample
    return samples

    
