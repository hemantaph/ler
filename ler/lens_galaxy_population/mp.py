import numpy as np
from scipy.integrate import quad
from numba import njit
from lenstronomy.LensModel.Solver.epl_shear_solver import caustics_epl_shear
from shapely.geometry import Polygon

from .jit_functions import phi_cut_SIE, axis_ratio_rayleigh, phi_q2_ellipticity_hemanta
from ..utils import inverse_transform_sampler, cubic_spline_interpolator


def optical_depth_sis_mp(params):

    # integrand
    zs = params[0]
    no = params[1]
    vd_inv_cdf = params[2]
    vd_inv_cdf_coeff = vd_inv_cdf[0]
    vd_list = vd_inv_cdf[1]
    vd_sampler = njit(lambda size_: inverse_transform_sampler(size_, vd_inv_cdf_coeff, vd_list))

    splineDa = params[4]
    splineDa_coeff = splineDa[0]
    splineDa_z_list = splineDa[1]
    Da_zs = cubic_spline_interpolator(np.array([zs]), splineDa_coeff, splineDa_z_list)[0]
    Da = njit(lambda zl_: cubic_spline_interpolator(np.array([zl_]), splineDa_coeff, splineDa_z_list)[0])

    splinedVcdz = params[3]
    dVcdz = njit(lambda zl_: cubic_spline_interpolator(np.array([zl_]), splinedVcdz[0], splinedVcdz[1])[0])

    # integrand
    # @njit makes it slower
    def tau_integrand(zl, zs):

        # velocity dispersion 
        size = 5000
        sigma = vd_sampler(size)

        # einstein radius 
        Dls = (Da_zs*(1+zs) - Da(zl)*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / Da_zs
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        # cross section 
        cross_section_SIS = theta_E ** 2  # np.pi is ommited

        result = cross_section_SIS/4 * no * dVcdz(zl)
        # average
        return np.mean(result)

    return(int(params[5]), quad(tau_integrand, 0, params[0], args=(params[0]))[0])

def optical_depth_sie1_mp(params):

    # integrand
    zs = params[0]
    no = params[1]
    vd_inv_cdf = params[2]
    vd_inv_cdf_coeff = vd_inv_cdf[0]  # coefficients of the inverse cdf
    vd_list = vd_inv_cdf[1]  # original list of velocity dispersions
    vd_sampler = njit(lambda size_: inverse_transform_sampler(size_, vd_inv_cdf_coeff, vd_list))

    splineDa = params[4]
    splineDa_coeff = splineDa[0]  # coefficients of the angular diameter distance spline interpolator
    splineDa_z_list = splineDa[1]  # redshift list of the angular diameter distance spline interpolator
    Da_zs = cubic_spline_interpolator(np.array([zs]), splineDa_coeff, splineDa_z_list)[0]
    Da = njit(lambda zl_: cubic_spline_interpolator(np.array([zl_]), splineDa_coeff, splineDa_z_list)[0])

    splinedVcdz = params[3]
    dVcdz = njit(lambda zl_: cubic_spline_interpolator(np.array([zl_]), splinedVcdz[0], splinedVcdz[1])[0])

    # integrand
    # @njit makes it slower
    def tau_integrand(zl, zs):
        """
        Function to calculate the optical depth for SIE lens with velocity dispersion distribution depending on redshift.

        Parameters
        ----------
        params : `list`
            list of parameters
            params[0] = zs (source redshift, float)
            params[1] = no (number density of lens galaxies, float)
            params[2] = vd_inv_cdf (velocity dispersion inverse cdf coefficients and redshift list, list). This vd_inv_cdf(s) of each redshift.
            params[3] = splineVcdz (differential comoving volume spline interpolator coefficients, list)
            params[4] = splineDa (angular diameter distance spline interpolator coefficients and redshift list, list)
            params[5] = idx (index to keep track of the operation, int)
            params[6] = zl_list (list of lens redshifts, list). This use for choosing the right vd_inv_cdf(s) for each lens redshifts.
        
        """

        # velocity dispersion 
        size = 5000
        sigma = vd_sampler(size)

        # axis ratio 
        q = axis_ratio_rayleigh(sigma)  # if SIS, q=array of 1.0

        # einstein radius 
        Dls = (Da_zs*(1+zs) - Da(zl)*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / Da_zs
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        # cross section 
        cross_section_SIS = theta_E ** 2  # np.pi is ommited

        result = phi_cut_SIE(q) * cross_section_SIS/4 * no * dVcdz(zl)
        # average
        return np.mean(result)

    return(int(params[5]), quad(tau_integrand, 0, params[0], args=(params[0]))[0])

def optical_depth_sie2_mp(params):
    """
    Function to calculate the optical depth for SIE lens with velocity dispersion distribution depending on redshift.

    Parameters
    ----------
    params : `list`
        list of parameters
        params[0] = zs (source redshift, float)
        params[1] = no (number density of lens galaxies, float)
        params[2] = vd_inv_cdf (velocity dispersion inverse cdf coefficients and redshift list, list). This vd_inv_cdf(s) of each redshift.
        params[3] = splineVcdz (differential comoving volume spline interpolator coefficients, list)
        params[4] = splineDa (angular diameter distance spline interpolator coefficients and redshift list, list)
        params[5] = idx (index to keep track of the operation, int)
        params[6] = zl_list (list of lens redshifts, list). This use for choosing the right vd_inv_cdf(s) for each lens redshifts.
    
    """

    # integrand
    zs = params[0]
    no = params[1]
    z_list = params[6]
    vd_inv_cdf = params[2]
    vd_inv_cdf_coeff = vd_inv_cdf[:,0]
    vd_list = vd_inv_cdf[0,1]
    vd_sampler = njit(lambda size_,zl_: inverse_transform_sampler(size_, vd_inv_cdf_coeff[np.searchsorted(z_list, zl_)], vd_list))

    splineDa = params[4]
    splineDa_coeff = splineDa[0]
    splineDa_z_list = splineDa[1]
    Da_zs = cubic_spline_interpolator(np.array([zs]), splineDa_coeff, splineDa_z_list)[0]
    Da = njit(lambda zl_: cubic_spline_interpolator(np.array([zl_]), splineDa_coeff, splineDa_z_list)[0])

    splinedVcdz = params[3]
    dVcdz = njit(lambda zl_: cubic_spline_interpolator(np.array([zl_]), splinedVcdz[0], splinedVcdz[1])[0])

    # @njit makes it slower
    def tau_integrand(zl, zs):

        # velocity dispersion #
        size = 5000
        sigma = vd_sampler(size, zl)

        # axis ratio 
        q = axis_ratio_rayleigh(sigma)  # if SIS, q=array of 1.0

        # einstein radius 
        Dls = (Da_zs*(1+zs) - Da(zl)*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / Da_zs
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        # cross section 
        cross_section_SIS = theta_E ** 2  # np.pi is ommited
        
        result = phi_cut_SIE(q) * cross_section_SIS/4 * no * dVcdz(zl)
        # average
        return np.mean(result)

    return(int(params[5]), quad(tau_integrand, 0, params[0], args=(params[0]))[0])


def optical_depth_sie3_mp(params):
    """
        Function to calculate the optical depth for SIE lens with velocity dispersion distribution depending on redshift.
    """

    # integrand
    zs = params[0]
    no = params[1]
    vd_inv_cdf = params[2]
    vd_inv_cdf_coeff = vd_inv_cdf[0]  # coefficients of the inverse cdf
    vd_list = vd_inv_cdf[1]  # original list of velocity dispersions
    vd_sampler = njit(lambda size_: inverse_transform_sampler(size_, vd_inv_cdf_coeff, vd_list))

    splineDa = params[4]
    splineDa_coeff = splineDa[0]  # coefficients of the angular diameter distance spline interpolator
    splineDa_z_list = splineDa[1]  # redshift list of the angular diameter distance spline interpolator
    Da_zs = cubic_spline_interpolator(np.array([zs]), splineDa_coeff, splineDa_z_list)[0]
    Da = njit(lambda zl_: cubic_spline_interpolator(np.array([zl_]), splineDa_coeff, splineDa_z_list)[0])

    splinedVcdz = params[3]
    dVcdz = njit(lambda zl_: cubic_spline_interpolator(np.array([zl_]), splinedVcdz[0], splinedVcdz[1])[0])

    # @njit makes it slower
    def sample_lens_params(zl, zs, size):

        # velocity dispersion #
        sigma = vd_sampler(size)
        # Sample axis ratios
        q = axis_ratio_rayleigh(sigma)
        # axis rotation angle
        phi = np.random.uniform(0.0, 2 * np.pi, size)
        # Transform the axis ratio and the angle, to ellipticities e1, e2, using lenstronomy
        e1, e2 = phi_q2_ellipticity_hemanta(phi, q)
        # Sample shears
        gamma1, gamma2 = np.random.normal(loc=0, scale=0.05,size=(2,size))
        # Sample the spectral index of the mass density distribution
        gamma = np.random.normal(loc=2.0, scale=0.2, size=size)

        # einstein radius 
        Dls = (Da_zs*(1+zs) - Da(zl)*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / Da_zs
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        return theta_E, e1, e2, gamma1, gamma2, gamma

    def tau_integrand(zl, zs):
        size=5000
        final_size = 0
        area_list = []

        while final_size < size:
            size = size - final_size
            theta_E, e1, e2, gamma1, gamma2, gamma = sample_lens_params(zl, zs, size=size)

            for i in range(size):
                kwargs_lens = [
                    {
                        "theta_E": theta_E[i],
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
                    kwargs_lens, return_which="double", maginf=-100
                )
                caustic = np.logical_not(np.isnan(caustic_double_points).any())

                # If there is a nan, caustic=False, draw a new gamma
                if caustic:
                    area_list.append(Polygon(caustic_double_points.T).area)

            final_size = len(area_list)

        result_EPL = np.array(area_list)/(4*np.pi) * no * dVcdz(zl)

        return np.mean(result_EPL[result_EPL<1.0]) #, np.mean(result_SIS)

    return(int(params[5]), quad(tau_integrand, 0, params[0], args=(params[0]))[0])

def optical_depth_sie4_mp(params):
    """
    Function to calculate the optical depth for EPL+Shear lens with velocity dispersion distribution depending on redshift.

    Parameters
    ----------
    params : `list`
        list of parameters
        params[0] = zs (source redshift, float)
        params[1] = no (number density of lens galaxies, float)
        params[2] = vd_inv_cdf (velocity dispersion inverse cdf coefficients and redshift list, list). This vd_inv_cdf(s) of each redshift.
        params[3] = splineVcdz (differential comoving volume spline interpolator coefficients, list)
        params[4] = splineDa (angular diameter distance spline interpolator coefficients and redshift list, list)
        params[5] = idx (index to keep track of the operation, int)
        params[6] = zl_list (list of lens redshifts, list). This use for choosing the right vd_inv_cdf(s) for each lens redshifts.
    
    """

    # integrand
    zs = params[0]
    no = params[1]
    z_list = params[6]
    vd_inv_cdf = params[2]
    vd_inv_cdf_coeff = vd_inv_cdf[:,0]
    vd_list = vd_inv_cdf[0,1]
    vd_sampler = njit(lambda size_,zl_: inverse_transform_sampler(size_, vd_inv_cdf_coeff[np.searchsorted(z_list, zl_)], vd_list))

    splineDa = params[4]
    splineDa_coeff = splineDa[0]
    splineDa_z_list = splineDa[1]
    Da_zs = cubic_spline_interpolator(np.array([zs]), splineDa_coeff, splineDa_z_list)[0]
    Da = njit(lambda zl_: cubic_spline_interpolator(np.array([zl_]), splineDa_coeff, splineDa_z_list)[0])

    splinedVcdz = params[3]
    dVcdz = njit(lambda zl_: cubic_spline_interpolator(np.array([zl_]), splinedVcdz[0], splinedVcdz[1])[0])

    # @njit makes it slower
    def sample_lens_params(zl, zs, size):

        # velocity dispersion #
        sigma = vd_sampler(size, zl)
        # Sample axis ratios
        q = axis_ratio_rayleigh(sigma)
        # axis rotation angle
        phi = np.random.uniform(0.0, 2 * np.pi, size)
        # Transform the axis ratio and the angle, to ellipticities e1, e2, using lenstronomy
        e1, e2 = phi_q2_ellipticity_hemanta(phi, q)
        # Sample shears
        gamma1, gamma2 = np.random.normal(loc=0, scale=0.05,size=(2,size))
        # Sample the spectral index of the mass density distribution
        gamma = np.random.normal(loc=2.0, scale=0.2, size=size)

        # einstein radius 
        Dls = (Da_zs*(1+zs) - Da(zl)*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / Da_zs
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        return theta_E, e1, e2, gamma1, gamma2, gamma

    def tau_integrand(zl, zs):
        size=5000
        final_size = 0
        area_list = []

        while final_size < size:
            size = size - final_size
            theta_E, e1, e2, gamma1, gamma2, gamma = sample_lens_params(zl, zs, size=size)

            for i in range(size):
                kwargs_lens = [
                    {
                        "theta_E": theta_E[i],
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
                    kwargs_lens, return_which="double", maginf=-100
                )
                caustic = np.logical_not(np.isnan(caustic_double_points).any())

                # If there is a nan, caustic=False, draw a new gamma
                if caustic:
                    area_list.append(Polygon(caustic_double_points.T).area)

            final_size = len(area_list)

        result_EPL = np.array(area_list)/(4*np.pi) * no * dVcdz(zl)

        return np.mean(result_EPL[result_EPL<1.0]) #, np.mean(result_SIS)

    return(int(params[5]), quad(tau_integrand, 0, params[0], args=(params[0]))[0])