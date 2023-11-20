from numba import njit, jit
import os
import pickle
from multiprocessing import Pool
import numpy as np
from scipy.integrate import quad
from scipy.stats import gengamma
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
cosmology_h = 0.7
from ..utils import  interpolator_from_pickle
from .jit_functions import phi_cut_SIE, axis_ratio_rayleigh, axis_ratio_SIS, gamma_, phi

class OpticalDepth():
    """
    Class to calculate the optical depth of a lens galaxy population.

    Parameters
    ----------
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
        z_min=0.001,
        z_max=10.0,
        functions=dict(
            strong_lensing_condition="rjs_with_cross_section",
            optical_depth="optical_depth_SIS_hemanta",
        ),
        sampler_priors=dict(
            velocity_dispersion="velocity_dispersion_bernardi",
            axis_ratio="axis_ratio_rayleigh",
        ),
        sampler_priors_params=dict(
            velocity_dispersion=dict(vd_min=60., vd_max=600.),
            axis_ratio=dict(q_min=0.01, q_max=1.),
        ),
        cosmology=None,
        directory="./interpolator_pickle",
        create_new_interpolator=dict(
            velocity_dispersion=dict(create_new=True, resolution=100), 
            optical_depth=dict(create_new=True, resolution=100), 
            z_to_Dc=dict(create_new=True, resolution=100), 
            Dc_to_z=dict(create_new=True, resolution=100),
            angular_diameter_distance=dict(create_new=True, resolution=100),
            differential_comoving_volume=dict(create_new=True, resolution=100),
        ),
        ):

        self.cosmo = cosmology if cosmology else cosmo
        global cosmology_h
        cosmology_h = self.cosmo.h
        self.sampler_priors = sampler_priors
        self.sampler_priors_params = sampler_priors_params
        self.directory = directory
        self.c_n_i = create_new_interpolator

        self.z_min = z_min
        self.z_max = z_max
        self.vd_min = sampler_priors_params['velocity_dispersion']['vd_min']
        self.vd_max = sampler_priors_params['velocity_dispersion']['vd_max']
        vd_name = sampler_priors["velocity_dispersion"]
        tau_name = functions["optical_depth"]

        self.create_lookup_table_fuction(self.z_max)

        self.rejection_sample_sl = getattr(
            self, functions["strong_lensing_condition"]
        )  # SL: Strong Lensing

        #######################
        # velocity dispersion #
        #######################
        # generating inverse cdf interpolator for velocity dispersion
        param_dict_given_ = dict(z_min=self.z_min, z_max=self.z_max, vd_min=self.vd_min, vd_max=self.vd_max, cosmology=cosmo, name=vd_name)
        sub_directory_ = vd_name
        x_ = np.linspace(self.vd_min, self.vd_max, self.c_n_i["velocity_dispersion"]["resolution"])
        category_ = "inv_cdf"
        create_new_ = self.c_n_i["velocity_dispersion"]["create_new"]

        # if velocity_dispersion_haris, you don't need to generate interpolator
        check = True  # check if the velocity dispersion name is in the library
        if vd_name == "velocity_dispersion_haris" or vd_name == "velocity_dispersion_gengamma":
            # no interpolation needed
            self.vd_inv_cdf = None

        else:      
            if vd_name == "velocity_dispersion_bernardi":
                # setting up input parameters for interpolation
                pdf_func_ = lambda vd_: vd_**4*phi_loc_bernardi(vd_)
                conditioned_y_ = None
                dimension_ = 1

            elif vd_name == "velocity_dispersion_ewoud":
                # setting up input parameters for interpolation
                pdf_func_ = lambda vd_, zl_: phi(vd_,zl_)*self.differential_comoving_volume(zl_)
                conditioned_y_ = np.linspace(self.z_min, self.z_max, 100)
                dimension_ = 2

            else:
                print("velocity dispersion name not in `ler` library.")
                check = False

            # this will initialize the interpolator
            self.vd_inv_cdf = interpolator_from_pickle(
                param_dict_given = param_dict_given_,
                directory=self.directory,
                sub_directory=sub_directory_,
                name=vd_name,
                x = x_,
                pdf_func= pdf_func_,
                conditioned_y=conditioned_y_,
                dimension=dimension_,
                category=category_,
                create_new=create_new_,
            )

        # this will initialize the sampler
        # check sampler_priors["velocity_dispersion"] is a callable function or not
        if check:
            pass
        else:
            if callable(sampler_priors["velocity_dispersion"]):
                pass
            else:
                # raise error and terminate
                raise ValueError("velocity dispersion name not in `ler` library. And sampler_priors['velocity_dispersion'] is not a callable function.")

        self.sample_velocity_dispersion = sampler_priors["velocity_dispersion"]

        ##############
        # for multiprocessing task
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
        # this z list is related to the velocity dispersion vs redshift interpolator
        with open("./mp_interpolator/redshift_list.pickle", "wb") as handle:
            pickle.dump(np.linspace(self.z_min, self.z_max, 100), handle, protocol=pickle.HIGHEST_PROTOCOL)
        ##############

        #################
        # optical depth #
        #################
        # setting up input parameters for interpolation
        param_dict_given_ = dict(z_min=self.z_min, z_max=self.z_max, vd_min=self.vd_min, vd_max=self.vd_max, cosmology=cosmo, tau_name=tau_name, vd_name=vd_name)
        sub_directory_ = tau_name
        if self.z_min==0.0:
            x_ = np.linspace(self.z_min+0.001, self.z_max, self.c_n_i["optical_depth"]["resolution"])
        else:
            x_ = np.geomspace(self.z_min, self.z_max, self.c_n_i["optical_depth"]["resolution"])
        pdf_func_ = self.optical_depth_multiprocessing
        dimension_ = 1
        category_ = "function"
        create_new_ = self.c_n_i["optical_depth"]["create_new"]

        # generating function interpolator for velocity dispersion
        if tau_name=="optical_depth_SIS_haris":
            self.sample_axis_ratio = axis_ratio_SIS
            # no interpolation needed
            self.strong_lensing_optical_depth = getattr(self, tau_name)

        else:
            if tau_name == "optical_depth_SIS_hemanta" or tau_name == "SIS":
                from .mp import optical_depth_sis_mp
                self.sample_axis_ratio = axis_ratio_SIS
                self.tau_mp_routine = optical_depth_sis_mp

            elif tau_name=="optical_depth_SIE_hemanta":
                # axis-ratio sampler
                if sampler_priors["axis_ratio"]=="axis_ratio_rayleigh":
                    self.sample_axis_ratio = axis_ratio_rayleigh
                else:
                    self.sample_axis_ratio = sampler_priors["axis_ratio"]

                if vd_name == "velocity_dispersion_ewoud":
                    from .mp import optical_depth_sie2_mp
                    self.tau_mp_routine = optical_depth_sie2_mp
                else:
                    from .mp import optical_depth_sie1_mp
                    self.tau_mp_routine = optical_depth_sie1_mp

            # this will initialize the interpolator
            self.strong_lensing_optical_depth = interpolator_from_pickle(
                param_dict_given = param_dict_given_,
                directory=self.directory,
                sub_directory=sub_directory_,
                name=tau_name,
                x = x_,
                pdf_func= pdf_func_,
                conditioned_y=None,
                dimension=dimension_,
                category=category_,
                create_new=create_new_,
            )

    def test(self):
        print("cosmology_h: ", cosmology_h)

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

    def velocity_dispersion_ewoud(self, size, zl, init_sigma=[]):
        """
        Function to sample velocity dispersion from the interpolator
        """

        z_max = self.z_max
        z_min = self.z_min
        zlist = np.linspace(z_min, z_max, 100)
        # find the index of z in zlist
        idx = np.searchsorted(zlist, zl)
        # get the interpolator (inverse cdf) and sample
        u = np.random.uniform(0, 1, size=size)
        sigma = self.vd_inv_cdf[idx](u)
        vd_min = self.vd_min
        vd_max = self.vd_max
        # choose sigma that is within the range
        sigma = np.concatenate((init_sigma,sigma[(sigma>vd_min) & (sigma<vd_max)]))
        size_ = size - len(sigma)
        if size_ <= 0:
            return sigma
        else:
            return self.velocity_dispersion_ewoud(size_, zl, sigma)

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
    
    def create_lookup_table_fuction(self, z_max):
        """
        Functions to create lookup tables
        1. Redshift to co-moving distance.
        2. Co-moving distance to redshift.
        3. Redshift to angular diameter distance.
        4. Lens redshift sampler helper function.
        """

        Dc = lambda z_: self.cosmo.comoving_distance(z_).value  # co-moving distance in Mpc
        self.z_to_Dc = interpolator_from_pickle(
            param_dict_given= dict(z_min=0.001, z_max=z_max, cosmology=self.cosmo),
            directory=self.directory,
            sub_directory="z_to_Dc", 
            name="z_to_Dc",
            x = np.linspace(0.001, z_max, self.c_n_i["z_to_Dc"]["resolution"]),
            pdf_func= Dc, 
            conditioned_y=None, 
            dimension=1,
            category="function",
            create_new=self.c_n_i["z_to_Dc"]["create_new"],
        )
        self.Dc_to_z = interpolator_from_pickle(
            param_dict_given= dict(z_min=0.001, z_max=z_max, cosmology=self.cosmo), 
            directory=self.directory,
            sub_directory="Dc_to_z", 
            name="Dc_to_z",
            x = np.linspace(0.001, z_max, self.c_n_i["Dc_to_z"]["resolution"]),
            pdf_func= Dc, 
            conditioned_y=None, 
            dimension=1,
            category="function_inverse",
            create_new=self.c_n_i["Dc_to_z"]["create_new"],
        )

        # for angular diameter distance
        Da = lambda z_: self.cosmo.angular_diameter_distance(z_).value
        self.angular_diameter_distance = interpolator_from_pickle(
                    param_dict_given= dict(z_min=0.001, z_max=z_max, cosmology=self.cosmo), 
                    directory=self.directory,
                    sub_directory="angular_diameter_distance", 
                    name="angular_diameter_distance",
                    x = np.linspace(0.001, z_max, self.c_n_i["angular_diameter_distance"]["resolution"]),
                    pdf_func= Da, 
                    conditioned_y=None, 
                    dimension=1,
                    category="function",
                    create_new=self.c_n_i["angular_diameter_distance"]["create_new"],
                )
        # for angular diameter distance between two redshifts
        self.angular_diameter_distance_z1z2 = lambda zl0, zs0: (self.angular_diameter_distance(zs0)*(1.+zs0) - self.angular_diameter_distance(zl0)*(1.+zl0))/(1.+zs0)

        # get differential co-moving volume interpolator
        self.differential_comoving_volume = interpolator_from_pickle(
            param_dict_given= dict(z_min=0.001, z_max=z_max, cosmology=self.cosmo), 
            directory=self.directory,
            sub_directory="differential_comoving_volume", 
            name="differential_comoving_volume",
            x = np.linspace(0.001, z_max, self.c_n_i["differential_comoving_volume"]["resolution"]),
            pdf_func= lambda z_: self.cosmo.differential_comoving_volume(z_).value * 4 * np.pi,
            conditioned_y=None, 
            dimension=1,
            category="function",
            create_new=self.c_n_i["differential_comoving_volume"]["create_new"],
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


@njit
def phi_loc_bernardi(sigma, alpha=0.94, beta=1.85, phistar=2.099e-2, sigmastar=113.78):
    phistar = phistar * (cosmology_h / 0.7) ** 3  # Mpc**-3
    philoc_ = phistar*(sigma/sigmastar)**alpha * np.exp(-(sigma/sigmastar)**beta) * beta/gamma_(alpha/beta)/sigma
    return philoc_





