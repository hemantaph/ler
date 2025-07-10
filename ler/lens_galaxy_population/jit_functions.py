import numpy as np
from numba import njit
#from scipy.stats import rayleigh
from ..utils import inverse_transform_sampler, cubic_spline_interpolator


@njit
def axis_ratio_SIS(sigma):
    """
    Function to sample axis ratio from the SIS distribution with given velocity dispersion.

    Parameters
    ----------
    sigma : `float: array`
        velocity dispersion of the lens galaxy

    Returns
    -------
    q : `float: array`
        axis ratio of the lens galaxy
    """

    return np.ones(len(sigma))
    
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
    # Coefficients for the fit. Use in the derivation velocity dispersion function (at local universe), Bernardi et al. (2010).
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
    # Derivative of the cvdf_fit function. Use in the derivation velocity dispersion function (at local universe), Bernardi et al. (2010).
    return 0.5 * (cvdf_fit(log_vd + dx, redshift) - cvdf_fit(log_vd - dx, redshift)) / dx

@njit
def pdf_phi_z_div_0(s, z):
    # Derivation of the pdf of velocity dispersion function (at redshift z), Oguri et al. (2018b). This lacks the scaling factor.
    log_vd = np.log10(s)
    phi_sim_z = 10 ** cvdf_fit(log_vd, z) / s * my_derivative(log_vd, z, 1e-8)
    phi_sim_0 = 10 ** cvdf_fit(log_vd, 0) / s * my_derivative(log_vd, 0, 1e-8)

    return phi_sim_z / phi_sim_0

@njit
def phi(s, z, alpha=0.94, beta=1.85, phistar=2.099e-2, sigmastar=113.78):
    """
    Function to calculate the lens galaxy velocity dispersion function at redshift z.
    For Oguri et al. (2018b) model: alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78

    Parameters
    ----------
    s : `float: array`
        velocity dispersion of the lens galaxy
    z : `float: array`
        redshift of the lens galaxy
    cosmology_h : `float`
        Hubble constant

    Returns
    -------
    result : `float: array`
    """
    
    result = pdf_phi_z_div_0(s,z) * phi_loc_bernardi(sigma=s, alpha=alpha, beta=beta, phistar=phistar, sigmastar=sigmastar)
    result[result < 0.] = 0.
    return result

@njit
def phi_loc_bernardi(sigma, alpha, beta, phistar, sigmastar):
    """
    Function to calculate the local universe velocity dispersion function. Bernardi et al. (2010).
    For Oguri et al. (2018b) model: alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78
    For Choi et al. (2008) model: alpha = 2.32 / 2.67, beta = 2.67, phistar = 8.0e-3*self.cosmo.h**3, sigmastar = 161.0

    Parameters
    ----------
    sigma : `float: array`
        velocity dispersion of the lens galaxy
    alpha, beta, phistar, sigmastar : `float`
        parameters of the velocity dispersion function
    cosmology_h : `float`
        Hubble constant with respect to 100 km/s/Mpc

    Returns
    -------
    philoc_ : `float: array`
    """

    # phistar = phistar * (cosmology_h / 0.7) ** 3  # Mpc**-3
    philoc_ = phistar*(sigma/sigmastar)**alpha * np.exp(-(sigma/sigmastar)**beta) * beta/gamma_(alpha/beta)/sigma
    return philoc_

# For elliptical lens galaxy
@njit
def phi_cut_SIE(q):
    """
    Function to calculate cross-section scaling factor for the SIE lens galaxy from SIS lens galaxy.

    Parameters
    ----------
    q : `float: array`
        axis ratio of the lens galaxy

    Returns
    -------
    result : `float: array`
        scaling factor
    """

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
def axis_ratio_rayleigh_rvs(sigma, q_min=0.2, q_max=1.0):
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
def axis_ratio_rayleigh_pdf(q, sigma, q_min=0.2, q_max=1.0):

    s = 0.38 - 5.7e-4 * sigma
    idx = (q >= q_min) & (q <= q_max)
    q_pdf = np.zeros_like(q)
    q_pdf[idx] = (1-q[idx])/s[idx]**2 * np.exp(-0.5*(1-q[idx])**2/s[idx]**2)
    return q_pdf


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
def lens_redshift_SDSS_catalogue_sis(zs, splineDc, splineDcInv, u, cdf):
    """
    Function to sample lens redshift from the SDSS catalogue. Haris et al. (2018) cdf = (10 * u**3 - 15 * u**4 + 6 * u**5)

    Parameters
    ----------
    zs: `numpy.ndarray` (1D array of float of size=size)
        Redshift of the source galaxy
    splineDc: `list`
        List of spline coefficients for the comoving distance and redshifts
    splineDcInv: `list`
        List of spline coefficients for the inverse of comoving distance and redshifts
    u: `numpy.ndarray` (1D array of float of size=size)
        corresponding x values wrt to the cdf values
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


@njit
def phi_q2_ellipticity_hemanta(phi, q):
    """Function to convert phi and q to ellipticity e1 and e2.

    Parameters
    ----------
    phi : `float: array`
        angle of the major axis in radians
    q : `float: array`
        axis ratio

    Returns
    -------
    e1 : `float: array`
        ellipticity component 1
    e2 : `float: array`
    """

    e_1 = (1.0 - q) / (1.0 + q) * np.cos(2 * phi)
    e_2 = (1.0 - q) / (1.0 + q) * np.sin(2 * phi)
    return e_1, e_2

# @njit
# def sample_sigma_zl(pdf, sigma_min, sigma_max, zl_min, zl_max, zs, chunk_size=10000):
#     x_sample = []
#     y_sample = []
    
#     for z in zs:
#         xmin, xmax = sigma_min, sigma_max
#         ymin, ymax = zl_min, z
        
#         # Keep sampling until a valid (sigma, zl) pair is found for this z
#         while True:
#             x_try = np.random.uniform(xmin, xmax, chunk_size)
#             y_try = np.random.uniform(ymin, ymax, chunk_size)
#             pdf_xy_try = pdf(x_try, y_try)
#             zmax = np.max(pdf_xy_try)  # Maximum of the PDF for the current batch

#             # Generate acceptance thresholds
#             z_try = np.random.uniform(0, zmax, chunk_size)

#             # Check which samples are accepted
#             accepted_indices = z_try < pdf_xy_try
#             if np.any(accepted_indices):
#                 # Accept the first valid sample
#                 x_sample.append(x_try[accepted_indices][0])
#                 y_sample.append(y_try[accepted_indices][0])
#                 break  # Exit the loop once the first valid sample is collected

#     return np.array(x_sample), np.array(y_sample)

@njit
def sample_sigma_zl(pdf, sigma_min, sigma_max, zl_min, zl_max, zs, chunk_size=10000):
    x_sample = []
    y_sample = []
    
    for z in zs:
        xmin, xmax = sigma_min, sigma_max
        ymin, ymax = zl_min, zl_max  # Use full range and then filter
        
        valid_sample_found = False
        while not valid_sample_found:
            x_try = np.random.uniform(xmin, xmax, chunk_size)
            y_try = np.random.uniform(ymin, ymax, chunk_size)
            pdf_xy_try = pdf(x_try, y_try)
            zmax = np.max(pdf_xy_try)  # Find maximum PDF value for current batch

            # Generate acceptance thresholds
            z_try = np.random.uniform(0, zmax, chunk_size)

            # Check and accept samples
            accepted = (z_try < pdf_xy_try) & (y_try < z)  # Ensure zl < zs condition
            if np.any(accepted):
                idx = np.argmax(accepted)  # First valid index
                x_sample.append(x_try[idx])
                y_sample.append(y_try[idx])
                valid_sample_found = True

    return np.array(x_sample), np.array(y_sample)