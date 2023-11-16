import numpy as np
from scipy.integrate import quad
# from multiprocessing import Pool, Manager
from ler.lens_galaxy_population import optical_depth
cosmology_h = optical_depth.cosmology_h
phi_cut_SIE = optical_depth.phi_cut_SIE
sample_axis_ratio = optical_depth.axis_ratio_rayleigh
# od = optical_depth.OpticalDepth()
# tau_integrand = od.tau_integrand

import pickle
with open('./mp_interpolator/differential_comoving_volume.pickle', 'rb') as input_:
    dVcdz = pickle.load(input_)
with open('./mp_interpolator/velocity_dispersion.pickle', 'rb') as input_:
    vd_inv_cdf = pickle.load(input_)
with open('./mp_interpolator/angular_diameter_distance.pickle', 'rb') as input_:
    Da = pickle.load(input_)
with open('./mp_interpolator/redshift_list.pickle', 'rb') as input_:
    zlist = pickle.load(input_)

def optical_depth_sis_mp(params):
    # integrand
    def tau_integrand(zl, zs):
        #######################
        # velocity dispersion #
        #######################
        u = np.random.uniform(0, 1, size=5000)
        sigma = vd_inv_cdf(u)
        ###################
        # einstein radius #
        ###################
        Ds = Da(zs)
        Dls = (Da(zs)*(1+zs) - Da(zl)*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / (Ds)
        )  # Note: km/s for sigma; Dls, Ds are in Mpc
        #################
        # cross section #
        #################
        cross_section_SIS = theta_E ** 2  # np.pi is ommited 
        #############################################
        no = params[1]
        result = cross_section_SIS/4 * no * dVcdz(zl)
        # average
        return np.mean(result)
    
    return quad(tau_integrand, 0, params[0], args=(params[0]))[0]

def optical_depth_sie1_mp(params):
    # integrand
    def tau_integrand(zl, zs):
        #######################
        # velocity dispersion #
        #######################
        u = np.random.uniform(0, 1, size=5000)
        sigma = vd_inv_cdf(u)
        ##############
        # axis ratio #
        ##############
        q = sample_axis_ratio(sigma)  # if SIS, q=array of 1.0
        ###################
        # einstein radius #
        ###################
        Ds = Da(zs)
        Dls = (Da(zs)*(1+zs) - Da(zl)*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / (Ds)
        )  # Note: km/s for sigma; Dls, Ds are in Mpc
        #################
        # cross section #
        #################
        cross_section_SIS = theta_E ** 2  # np.pi is ommited
        #############################################
        no = params[1]
        result = phi_cut_SIE(q) * cross_section_SIS/4 * no * dVcdz(zl)
        # average
        return np.mean(result)

    return quad(tau_integrand, 0, params[0], args=(params[0]))[0]

def optical_depth_sie2_mp(params):
    # integrand
    def tau_integrand(zl, zs):
        #######################
        # velocity dispersion #
        #######################
        u = np.random.uniform(0, 1, size=5000)
        idx = np.searchsorted(zlist, zl)
        sigma = vd_inv_cdf[idx](u)
        ##############
        # axis ratio #
        ##############
        q = sample_axis_ratio(sigma)  # if SIS, q=array of 1.0
        ###################
        # einstein radius #
        ###################
        Ds = Da(zs)
        Dls = (Da(zs)*(1+zs) - Da(zl)*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / (Ds)
        )  # Note: km/s for sigma; Dls, Ds are in Mpc
        #################
        # cross section #
        #################
        cross_section_SIS = theta_E ** 2  # np.pi is ommited
        #############################################
        no = params[1]
        result = phi_cut_SIE(q) * cross_section_SIS/4 * no * dVcdz(zl)
        # average
        return np.mean(result)

    return quad(tau_integrand, 0, params[0], args=(params[0]))[0]


def worker(x, shared_dict):
    # Modify or access shared data
    shared_dict['result'] = x * shared_dict['attribute']
    return shared_dict['result']

#####################
# Pass as Arguments #
#####################
# don't work
# class MyClass:
#     def __init__(self, attribute):
#         self.attribute = attribute

#     def worker(self, x):
#         # Do something with self.attribute and x
#         return x * self.attribute

#     def run_parallel(self, data):
#         with Pool() as pool:
#             # Pass self.attribute as an argument
#             results = pool.map(lambda x: self.worker(x), data)
#         return results

########################
# Use Global Variables #
########################
# don't work
# class_attribute = None

# def worker(x):
#     global class_attribute
#     return x * class_attribute

# class MyClass:
#     def __init__(self, attribute):
#         self.attribute = attribute

#     def run_parallel(self, data):
#         global class_attribute
#         class_attribute = self.attribute
#         with Pool() as pool:
#             results = pool.map(worker, data)
#         return results

#################
# Use a Manager #
#################
# class MyClass:
#     def __init__(self, attribute):
#         self.attribute = attribute

#     def worker(self, x, attribute):
#         return x * attribute

#     def run_parallel(self, data):
#         with Manager() as manager:
#             # Create a managed object
#             attribute = manager.Value('i', self.attribute)
#             with Pool() as pool:
#                 # Pass the managed object as an argument
#                 results = pool.starmap(self.worker, [(x, attribute.value) for x in data])
#             return results

