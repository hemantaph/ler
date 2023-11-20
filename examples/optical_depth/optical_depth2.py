import os
import pickle
import numpy as np
from multiprocessing import Pool
from scipy.integrate import quad
from scipy.stats import gengamma, rayleigh
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
from numba import njit, jit
from ler.utils import  interpolator_from_pickle
#import mp

class OpticalDepth():
    """
    Class to calculate the optical depth of a lens galaxy population.

    Parameters
    ----------
    nsamples_z : int
        Number of samples in redshift for the interpolator
    nsamples_sigma : int
        Number of samples in velocity dispersion for the interpolator
    param_dict_given : dict
        Dictionary of parameters for the interpolator
    directory : str
        Directory to save the interpolator
    optical_depth_method : str
        Method to calculate the optical depth
        e.g. 'optical_depth_SIS_haris', 'optical_depth_SIS_hemata', 'optical_depth_SIE_hemanta'
    vel_disp_method : str
        Method to sample velocity dispersion
        e.g. 'velocity_dispersion_gengamma', 'velocity_dispersion_bernardi', 'velocity_dispersion_ewoud' 
    """

    def __init__(self, 
        nsamples_z=100,
        nsamples_sigma=100,
        functions=dict(
            strong_lensing_condition="rjs_with_cross_section",
            optical_depth="optical_depth_SIS_hemanta",
        ),
        sampler_priors=dict(
            lens_redshift="lens_redshift_SDSS_catalogue",
            velocity_dispersion="velocity_dispersion_bernardi",
            axis_ratio="axis_ratio_rayleigh",
        ),
        sampler_priors_params=dict(
            lens_redshift=dict(z_min=0.001, z_max=10),
            velocity_dispersion=dict(vd_min=60, vd_max=600),
            axis_ratio=dict(q_min=0.01, q_max=1),    
        ),
        cosmology=cosmo,
        directory="./interpolator_pickle",
        create_new_interpolator=dict(velocity_dispersion=False, optical_depth=False, z_to_Dc=False, Dc_to_z=False, angular_diameter_distance=False, differential_comoving_volume=False)
        ):

        self.sampler_priors = sampler_priors
        self.sampler_priors_params = sampler_priors_params
        self.directory = directory
        self.nsamples_z = nsamples_z
        self.nsamples_sigma = nsamples_sigma
        self.cosmo = cosmology
        self.create_new_interpolator = create_new_interpolator

        z_min = sampler_priors_params['lens_redshift']['z_min']
        z_max = sampler_priors_params['lens_redshift']['z_max']
        vd_min = sampler_priors_params['velocity_dispersion']['vd_min']
        vd_max = sampler_priors_params['velocity_dispersion']['vd_max']
        vd_name = sampler_priors["velocity_dispersion"]
        tau_name = functions["optical_depth"]

        self.create_lookup_table_fuction(z_max, create_new=create_new_interpolator);

        self.rejection_sample_sl = getattr(
            self, functions["strong_lensing_condition"]
        )  # SL: Strong Lensing

        #######################
        # velocity dispersion #
        #######################
        # generating inverse cdf interpolator for velocity dispersion
        # if velocity_dispersion_haris, you don't need to generate interpolator
        if vd_name == "velocity_dispersion_haris" or vd_name == "velocity_dispersion_gengamma":
            self.vd_inv_cdf = None
                
        elif vd_name == "velocity_dispersion_bernardi":
            self.vd_inv_cdf = interpolator_from_pickle(
                param_dict_given= dict(z_min=z_min, 
                                       z_max=z_max, 
                                       vd_min=vd_min, 
                                       vd_max=vd_max,
                                       cosmology=cosmo,
                                       name=vd_name
                                       ),
                directory=directory,
                sub_directory=vd_name,
                name=vd_name,
                x = np.linspace(vd_min, vd_max, nsamples_sigma),
                pdf_func= lambda vd_: vd_**4*phi_loc_bernardi(vd_),
                conditioned_y=None, 
                dimension=1,
                category="inv_cdf",
                create_new=create_new_interpolator["velocity_dispersion"],
            )

        elif vd_name == "velocity_dispersion_ewoud":
            self.vd_inv_cdf = interpolator_from_pickle(
                param_dict_given= dict(z_min=z_min, 
                                       z_max=z_max, 
                                       vd_min=vd_min, 
                                       vd_max=vd_max,
                                       cosmology=cosmo,
                                       name=vd_name
                                       ),
                directory=directory,
                sub_directory=vd_name,
                name=vd_name,
                x = np.linspace(vd_min, vd_max, nsamples_sigma),
                pdf_func= lambda vd_, zl_: phi(vd_,zl_)*self.differential_comoving_volume(zl_), 
                conditioned_y=np.linspace(z_min, z_max, self.nsamples_z), 
                dimension=2,
                category="inv_cdf",
                create_new=create_new_interpolator["velocity_dispersion"],
            )

        # this will initialize the interpolator
        self.sample_velocity_dispersion = sampler_priors["velocity_dispersion"]

        ##############
        # for multiprocessing test
        # pickle dump interpolator
        # create dir
        if not os.path.exists("./mp_interpolator"):
            os.makedirs("./mp_interpolator")
            
        with open("./mp_interpolator/differential_comoving_volume.pickle", "wb") as handle:
            pickle.dump(self.differential_comoving_volume, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open("./mp_interpolator/velocity_dispersion.pickle", "wb") as handle:
            pickle.dump(self.vd_inv_cdf, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open("./mp_interpolator/angular_diameter_distance.pickle", "wb") as handle:
            pickle.dump(self.angular_diameter_distance, handle, protocol=pickle.HIGHEST_PROTOCOL)
        with open("./mp_interpolator/redshift_list.pickle", "wb") as handle:
            pickle.dump(np.linspace(z_min, z_max, self.nsamples_z), handle, protocol=pickle.HIGHEST_PROTOCOL)
        ##############

        #################
        # optical depth #
        #################
        # generating function interpolator for velocity dispersion
        if tau_name == "optical_depth_SIS_hemanta" or tau_name == "SIS":
            self.sample_axis_ratio = axis_ratio_SIS
            self.tau_mp_routine = mp.optical_depth_sis_mp

            self.optical_depth = interpolator_from_pickle(
                param_dict_given = dict(z_min=z_min,
                                       z_max=z_max,
                                       vd_min=vd_min,
                                       vd_max=vd_max,
                                       cosmology=cosmo,
                                       name=tau_name
                                       ), 
                directory="./interpolator",
                sub_directory=tau_name, 
                name=tau_name,
                x = np.geomspace(z_min, z_max, 50),
                pdf_func= self.optical_depth_multiprocessing,
                conditioned_y=None,
                dimension=1,
                category="function",
                create_new=create_new_interpolator["optical_depth"],
            )
        elif tau_name=="optical_depth_SIS_haris":
            self.sample_axis_ratio = axis_ratio_SIS
            self.optical_depth = getattr(self, tau_name)

        elif tau_name=="optical_depth_SIE_hemanta":
            # axis-ratio sampler
            if sampler_priors["axis_ratio"]=="axis_ratio_rayleigh":
                self.sample_axis_ratio = axis_ratio_rayleigh
            else:
                self.sample_axis_ratio = sampler_priors["axis_ratio"]
            if vd_name == "velocity_dispersion_ewoud":
                self.tau_mp_routine = mp.optical_depth_sie2_mp
            else:
                self.tau_mp_routine = mp.optical_depth_sie1_mp

            self.optical_depth = interpolator_from_pickle(
                param_dict_given = dict(z_min=z_min, 
                                       z_max=z_max, 
                                       vd_min=vd_min, 
                                       vd_max=vd_max,
                                       cosmology=cosmo,
                                       tau_name=tau_name,
                                       vd_name=vd_name,
                                       ), 
                directory="./interpolator",
                sub_directory=tau_name, 
                name=tau_name,
                x = np.geomspace(z_min, z_max, 50),
                pdf_func= self.optical_depth_multiprocessing, 
                conditioned_y=None, 
                dimension=1,
                category="function",
                create_new=create_new_interpolator["optical_depth"],
            )

    def rjs_with_cross_section(self, param_dict):
        """
        Function to conduct rejection sampling wrt einstein radius

        Parameters
        ----------
        param_dict : `dict`
            dictionary of lens parameters

        Returns
        -------
        lens_params : `dict`
            dictionary of lens parameters after rejection sampling
        """

        theta_E = param_dict["theta_E"]
        size = len(theta_E)
        theta_E_max = np.max(theta_E)  # maximum einstein radius
        u = np.random.uniform(0, theta_E_max**2, size=size)
        mask = u < theta_E**2

        # return the dictionary with the mask applied
        return {key: val[mask] for key, val in param_dict.items()}

    def velocity_dispersion_haris(self, size, a=2.32 / 2.67, c=2.67, param=None, **kwargs):
        """
        Function to sample velocity dispersion from gengamma distribution

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        a,c : `float`
            parameters of gengamma distribution
            refer to https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gengamma.html

        Returns
        -------
        sigma : `array`
            velocity dispersion of the lens galaxy
        """

        if param:
            a = param["a"]
            c = param["c"]

        # sample velocity dispersion from gengamma distribution
        return 161.*gengamma.rvs(a, c, size=size)  # km/s
    
    def velocity_dispersion_bernardi(self, size, **kwargs):
        """
        Function to sample velocity dispersion from the interpolator
        """

        # get the interpolator (inverse cdf) and sample
        return self.vd_inv_cdf(np.random.uniform(0, 1, size=size))

    def velocity_dispersion_ewoud(self, size, zl):
        """
        Function to sample velocity dispersion from the interpolator
        """

        z_max = self.sampler_priors_params['lens_redshift']['z_max']
        z_min = self.sampler_priors_params['lens_redshift']['z_min']
        zlist = np.linspace(z_min, z_max, self.nsamples_z)
        # find the index of z in zlist
        idx = np.searchsorted(zlist, zl)
        # get the interpolator (inverse cdf) and sample
        u = np.random.uniform(0, 1, size=size)

        return self.vd_inv_cdf[idx](u)

    # SIS crossection 
    def cross_section_SIS(self, sigma, zl, zs):

        Ds = self.angular_diameter_distance(zs)
        Dls = self.angular_diameter_distance_z1z2(zl, zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / (Ds)
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        return np.pi * theta_E**2
    
    def tau_integrand(self, zl, zs):
            # size=5000 will take ~ 48s to run, with ewoud vd sampler
            sigma = self.sample_velocity_dispersion(size=5000, zl=zl)
            q = self.sample_axis_ratio(sigma)  # if SIS, q=array of 1.0
            no = 8*1e-3*self.cosmo.h**3
            test = phi_cut_SIE(q)*self.cross_section_SIS(sigma=sigma, zl=zl, zs=zs)/(4*np.pi)*no*self.differential_comoving_volume(zl)
            # average
            return np.mean(test)
    
    def optical_depth_calculator(self, zs):

        zs = np.array([zs]).reshape(-1)
        tau_list = []
        for z in zs:
            tau_list.append(quad(self.tau_integrand, 0, z, args=(z))[0])

        return np.array(tau_list)
    
    def optical_depth_multiprocessing(self, zs):

        zs = np.array([zs]).reshape(-1)
        no = np.ones(len(zs))*8*1e-3*self.cosmo.h**3
        input_params = np.array([zs, no]).T

        # Create a pool of workers and parallelize the integration
        with Pool(processes=4) as pool:
            tau_list = list(pool.imap(self.tau_mp_routine, input_params))

        return np.array(tau_list)
    

    def optical_depth_SIS_haris(self, zs):
        """
        Function to compute the strong lensing optical depth (SIS). \n
        LambdaCDM(H0=70, Om0=0.3, Ode0=0.7) was used to derive the following equation.

        Parameters
        ----------
        zs : `float`
            source redshifts

        Returns
        -------
        tau : `float`
            strong lensing optical depth
        """

        # z to luminosity_distance (luminosity_distance) conversion
        Dc = self.z_to_Dc(zs) * 1e-3  # 1e-3 converts Mpc to Gpc

        return (Dc / 62.2) ** 3  # 62.2 is the critical density in Gpc
    
    def create_lookup_table_fuction(self, z_max, create_new=False):
        """
        Functions to create lookup tables
        1. Redshift to co-moving distance.
        2. Co-moving distance to redshift.
        3. Redshift to angular diameter distance.
        4. Lens redshift sampler helper function.
        """

        z = np.linspace(0.0, z_max, 500)  # red-shift
        Dc = lambda z_: self.cosmo.comoving_distance(z_).value  # co-moving distance in Mpc
        self.z_to_Dc = interpolator_from_pickle(
            param_dict_given= dict(z_min=0.001, z_max=z_max, cosmology=self.cosmo),
            directory="./interpolator",
            sub_directory="z_to_Dc", 
            name="z_to_Dc",
            x = z,
            pdf_func= Dc, 
            conditioned_y=None, 
            dimension=1,
            category="function",
            create_new=create_new["z_to_Dc"],
        )
        self.Dc_to_z = interpolator_from_pickle(
            param_dict_given= dict(z_min=0.001, z_max=z_max, cosmology=self.cosmo), 
            directory="./interpolator",
            sub_directory="Dc_to_z", 
            name="Dc_to_z",
            x = z,
            pdf_func= Dc, 
            conditioned_y=None, 
            dimension=1,
            category="function_inverse",
            create_new=create_new["Dc_to_z"],
        )

        # for angular diameter distance
        Da = lambda z_: self.cosmo.angular_diameter_distance(z_).value
        self.angular_diameter_distance = interpolator_from_pickle(
                    param_dict_given= dict(z_min=0.001, z_max=z_max, cosmology=self.cosmo), 
                    directory="./interpolator",
                    sub_directory="angular_diameter_distance", 
                    name="angular_diameter_distance",
                    x = z,
                    pdf_func= Da, 
                    conditioned_y=None, 
                    dimension=1,
                    category="function",
                    create_new=create_new["angular_diameter_distance"],
                )
        # for angular diameter distance between two redshifts
        self.angular_diameter_distance_z1z2 = lambda zl0, zs0: (self.angular_diameter_distance(zs0)*(1.+zs0) - self.angular_diameter_distance(zl0)*(1.+zl0))/(1.+zs0)

        # get differential co-moving volume interpolator
        self.differential_comoving_volume = interpolator_from_pickle(
            param_dict_given= dict(z_min=0.001, z_max=z_max, cosmology=self.cosmo), 
            directory="./interpolator",
            sub_directory="differential_comoving_volume", 
            name="differential_comoving_volume",
            x = z,
            pdf_func= lambda z_: self.cosmo.differential_comoving_volume(z_).value * 4 * np.pi,
            conditioned_y=None, 
            dimension=1,
            category="function",
            create_new=create_new["differential_comoving_volume"],
        )

    @property
    def sample_lens_redshift(self):
        """
        Function to sample lens redshifts, conditioned on the lens being strongly lensed
        """
        return self._sample_lens_redshift

    @sample_lens_redshift.setter
    def sample_lens_redshift(self, prior):
        try:
            self._sample_lens_redshift = getattr(self, prior)
        except:
            self._sample_lens_redshift = prior

    @property
    def sample_velocity_dispersion(self):
        """
        Function to sample velocity dispersion from gengamma distribution
        """
        return self._sample_velocity_dispersion

    @sample_velocity_dispersion.setter
    def sample_velocity_dispersion(self, prior):
        try:
            self._sample_velocity_dispersion = getattr(self, prior)
        except:
            self._sample_velocity_dispersion = prior

    @property
    def sample_axis_ratio(self):
        """
        Function to sample axis ratio from rayleigh distribution with given velocity dispersion.
        """
        return self._sample_axis_ratio

    @sample_axis_ratio.setter
    def sample_axis_ratio(self, prior):
        try:
            self._sample_axis_ratio = getattr(self, prior)
        except:
            self._sample_axis_ratio = prior


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

@njit
def axis_ratio_SIS(sigma):
    return np.ones(len(sigma))

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
def phi_loc_bernardi(sigma, alpha=0.94, beta=1.85, phistar=2.099e-2, sigmastar=113.78):
    
    phistar = phistar * (cosmology_h / 0.7) ** 3  # Mpc**-3
    philoc_ = phistar*(sigma/sigmastar)**alpha * np.exp(-(sigma/sigmastar)**beta) * beta/gamma_(alpha/beta)/sigma
    return philoc_

@njit
def phi(s,z):
    return s**4*pdf_phi_z_div_0(s,z)*phi_loc_bernardi(s)