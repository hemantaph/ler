import numpy as np
from scipy.integrate import quad

from .jit_functions import phi_cut_SIE, axis_ratio_rayleigh

import pickle
# call the interpolator
with open('./mp_interpolator/differential_comoving_volume.pickle', 'rb') as input_:
    dVcdz = pickle.load(input_)
with open('./mp_interpolator/velocity_dispersion.pickle', 'rb') as input_:
    vd_inv_cdf = pickle.load(input_)
with open('./mp_interpolator/angular_diameter_distance.pickle', 'rb') as input_:
    Da = pickle.load(input_)
with open('./mp_interpolator/redshift_list.pickle', 'rb') as input_:
    zlist = pickle.load(input_)

vd_min = 0.0
vd_max = 600.0

def optical_depth_sis_mp(params):
    # integrand
    def tau_integrand(zl, zs):

        # velocity dispersion 
        u = np.random.uniform(0, 1, size=5000)
        sigma = vd_inv_cdf(u)

        # einstein radius 
        Ds = Da(zs)
        Dls = (Da(zs)*(1+zs) - Da(zl)*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / (Ds)
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        # cross section 
        cross_section_SIS = theta_E ** 2  # np.pi is ommited 

        no = params[1]
        result = cross_section_SIS/4 * no * dVcdz(zl)
        # average
        return np.mean(result)
    
    return quad(tau_integrand, 0, params[0], args=(params[0]))[0]

def optical_depth_sie1_mp(params):
    # integrand
    def tau_integrand(zl, zs):

        # velocity dispersion 
        u = np.random.uniform(0, 1, size=5000)
        sigma = vd_inv_cdf(u)

        # axis ratio 
        q = axis_ratio_rayleigh(sigma)  # if SIS, q=array of 1.0

        # einstein radius 
        Ds = Da(zs)
        Dls = (Da(zs)*(1+zs) - Da(zl)*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / (Ds)
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        # cross section 
        cross_section_SIS = theta_E ** 2  # np.pi is ommited

        no = params[1]
        result = phi_cut_SIE(q) * cross_section_SIS/4 * no * dVcdz(zl)
        # average
        return np.mean(result)

    return quad(tau_integrand, 0, params[0], args=(params[0]))[0]

def optical_depth_sie2_mp(params):
    # integrand
    def tau_integrand(zl, zs):

        # velocity dispersion #
        size = 5000
        u = np.random.uniform(0, 1, size=size)
        idx = np.searchsorted(zlist, zl)
        sigma = np.array(vd_inv_cdf[idx](u))
        sigma = sigma[(sigma>vd_min) & (sigma<vd_max)]

        # axis ratio 
        q = axis_ratio_rayleigh(sigma)  # if SIS, q=array of 1.0

        # einstein radius 
        Ds = Da(zs)
        Dls = (Da(zs)*(1+zs) - Da(zl)*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / (Ds)
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        # cross section 
        cross_section_SIS = theta_E ** 2  # np.pi is ommited

        no = params[1]
        result = phi_cut_SIE(q) * cross_section_SIS/4 * no * dVcdz(zl)
        # average
        return np.mean(result)

    return quad(tau_integrand, 0, params[0], args=(params[0]))[0]


