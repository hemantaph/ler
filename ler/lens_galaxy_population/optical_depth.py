from numba import njit
from multiprocessing import Pool
from tqdm import tqdm
import numpy as np
from scipy.integrate import quad
from scipy.stats import gengamma
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
from ..utils import  interpolator_from_pickle, cubic_spline_interpolator, inverse_transform_sampler
from .jit_functions import phi_cut_SIE, axis_ratio_rayleigh, axis_ratio_SIS, gamma_, phi, phi_loc_bernardi

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
        npool=4,
        z_min=0.001,
        z_max=10.0,
        functions=dict(
            optical_depth="optical_depth_SIE_hemanta",
        ),
        sampler_priors=dict(
            velocity_dispersion="velocity_dispersion_bernardi",
            axis_ratio="axis_ratio_rayleigh",
        ),
        sampler_priors_params=dict(
            velocity_dispersion=dict(vd_min=0., vd_max=600.),
            axis_ratio=dict(q_min=0.01, q_max=1.),
        ),
        cosmology=None,
        directory="./interpolator_pickle",
        create_new_interpolator=dict(
            velocity_dispersion=dict(create_new=False, resolution=100), 
            optical_depth=dict(create_new=False, resolution=100), 
            z_to_Dc=dict(create_new=False, resolution=500), 
            Dc_to_z=dict(create_new=False, resolution=500),
            angular_diameter_distance=dict(create_new=False, resolution=500),
            differential_comoving_volume=dict(create_new=False, resolution=500),
        ),
        ):

        self.npool = npool
        self.cosmo = cosmology if cosmology else cosmo
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

        #######################
        # velocity dispersion #
        #######################
        # generating inverse cdf interpolator for velocity dispersion
        param_dict_given_ = dict(z_min=self.z_min, z_max=self.z_max, vd_min=self.vd_min, vd_max=self.vd_max, cosmology=cosmo, name=vd_name, resolution=self.c_n_i["velocity_dispersion"]["resolution"])
        sub_directory_ = vd_name
        x_ = np.linspace(self.vd_min, self.vd_max, self.c_n_i["velocity_dispersion"]["resolution"])
        category_ = "inv_cdf"
        create_new_ = self.c_n_i["velocity_dispersion"]["create_new"]
        
        if vd_name == "velocity_dispersion_gengamma":
            # setting up input parameters for interpolation
            pdf_func_ = lambda vd_: gengamma.pdf(vd_/161., a=2.32 / 2.67, c=2.67)   # gengamma pdf
            conditioned_y_ = None
            dimension_ = 1
     
        elif vd_name == "velocity_dispersion_bernardi":
            # setting up input parameters for interpolation
            pdf_func_ = lambda vd_: vd_**4*phi_loc_bernardi(sigma=vd_, cosmology_h=cosmology_h)
            conditioned_y_ = None
            dimension_ = 1

        elif vd_name == "velocity_dispersion_ewoud":
            if self.z_min==0.0:
                self.zl_list = np.linspace(self.z_min+0.001, self.z_max, 100)
            else:
                self.zl_list = np.linspace(self.z_min, self.z_max, 100)
            # setting up input parameters for interpolation
            pdf_func_ = lambda vd_, zl_: phi(vd_, zl_, cosmology_h=cosmology_h)*self.differential_comoving_volume(np.array([zl_]))
            conditioned_y_ = self.zl_list
            dimension_ = 2

        else:
            # this is to allow user to pass in their own sampler
            if callable(sampler_priors["velocity_dispersion"]):
                pass
            else:
                # raise error and terminate
                raise ValueError("velocity dispersion name not in `ler` library. And sampler_priors['velocity_dispersion'] is not a callable function.")
            
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

        self.sample_velocity_dispersion = sampler_priors["velocity_dispersion"]

        #################
        # optical depth #
        #################
        if tau_name=="optical_depth_SIS_haris":
            self.sample_axis_ratio = axis_ratio_SIS
            # no interpolation needed
            optical_depth_setter = getattr(self, tau_name)

        elif callable(tau_name):
            self.sample_axis_ratio = axis_ratio_SIS
            optical_depth_setter = tau_name

        else:
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

            if tau_name == "optical_depth_SIS_hemanta" or tau_name == "SIS":
                from .mp import optical_depth_sis_mp
                self.sample_axis_ratio = axis_ratio_SIS
                self.tau_mp_routine = optical_depth_sis_mp

            elif tau_name=="optical_depth_SIE_hemanta" or tau_name == "SIE":
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
            optical_depth_setter = interpolator_from_pickle(
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
                
        self.strong_lensing_optical_depth = optical_depth_setter


    def velocity_dispersion_gengamma(self, size, a=2.32 / 2.67, c=2.67, get_attribute=False, param=None, **kwargs):
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

        if get_attribute:
            return lambda size: 161.*gengamma.rvs(a, c, size=size)
        else:
            # sample velocity dispersion from gengamma distribution
            return 161.*gengamma.rvs(a, c, size=size)  # km/s
    
    def velocity_dispersion_bernardi(self, size, get_attribute=False, **kwargs):
        """
        Function to sample velocity dispersion from the interpolator
        """

        vd_inv_cdf = self.vd_inv_cdf.copy()
        # get the interpolator (inverse cdf) and sample
        if get_attribute:
            return njit(lambda size: inverse_transform_sampler(size, vd_inv_cdf[0], vd_inv_cdf[1]))
        else:
            return inverse_transform_sampler(size, vd_inv_cdf[0], vd_inv_cdf[1])

    def velocity_dispersion_ewoud(self, size, zl, get_attribute=False):
        """
        Function to sample velocity dispersion from the interpolator
        """

        z_max = self.z_max
        z_min = self.z_min
        zlist = np.linspace(z_min, z_max, 100)
        vd_inv_cdf = self.vd_inv_cdf.copy()

        @njit
        def sampler(size, zl, inv_cdf, zlist):

            idx = np.searchsorted(zlist, zl)
            return inverse_transform_sampler(size, inv_cdf[idx][0], inv_cdf[idx][1])
        
        if get_attribute:
            return njit(lambda size, zl: sampler(size, zl, vd_inv_cdf, zlist))
        else:
            return sampler(size, zl, vd_inv_cdf, zlist)

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
        no = 8*1e-3*self.cosmo.h**3
        vd_inv_cdf = self.vd_inv_cdf
        splinedVcdz = self.splinedVcdz
        splineDa = self.splineDa
        idx = np.arange(len(zs))
        try:
            zl_list = self.zl_list
            input_params = [(zs[i], no, vd_inv_cdf, splinedVcdz, splineDa, idx[i], zl_list) for i in range(len(zs))]
        except:
            input_params = [(zs[i], no, vd_inv_cdf, splinedVcdz, splineDa, idx[i]) for i in range(len(zs))]

        # Create a pool of workers and parallelize the integration
        with Pool(processes=self.npool) as pool:
            result = list(pool.map(self.tau_mp_routine, input_params))

        result = np.array(result)
        tau_list = result[:,1][np.array(result[:,0], dtype=int)]

        return tau_list
    

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
        resolution = self.c_n_i["z_to_Dc"]["resolution"]
        create_new = self.c_n_i["z_to_Dc"]["create_new"]
        splineDc = interpolator_from_pickle(
            param_dict_given= dict(z_min=0.001, z_max=z_max, cosmology=self.cosmo, resolution=resolution),
            directory=self.directory,
            sub_directory="z_to_Dc", 
            name="z_to_Dc",
            x = np.geomspace(0.001, z_max, resolution),
            pdf_func= Dc, 
            conditioned_y=None,
            dimension=1,
            category="function",
            create_new= create_new,
        )
        self.z_to_Dc = njit(lambda z_: cubic_spline_interpolator(z_, splineDc[0], splineDc[1]))

        # for co-moving distance to redshift
        resolution = self.c_n_i["Dc_to_z"]["resolution"]
        create_new = self.c_n_i["Dc_to_z"]["create_new"]
        splineDcInv = interpolator_from_pickle(
            param_dict_given= dict(z_min=0.001, z_max=z_max, cosmology=self.cosmo, resolution=resolution), 
            directory=self.directory,
            sub_directory="Dc_to_z", 
            name="Dc_to_z",
            x = np.geomspace(0.001, z_max, resolution),
            pdf_func= Dc, 
            conditioned_y=None, 
            dimension=1,
            category="function_inverse",
            create_new=create_new,
        )
        self.Dc_to_z = njit(lambda Dc_: cubic_spline_interpolator(Dc_, splineDcInv[0], splineDcInv[1]))

        # for angular diameter distance
        resolution = self.c_n_i["angular_diameter_distance"]["resolution"]
        create_new = self.c_n_i["angular_diameter_distance"]["create_new"]
        Da = lambda z_: self.cosmo.angular_diameter_distance(z_).value
        splineDa = interpolator_from_pickle(
                    param_dict_given= dict(z_min=0.001, z_max=z_max, cosmology=self.cosmo, resolution=resolution), 
                    directory=self.directory,
                    sub_directory="angular_diameter_distance", 
                    name="angular_diameter_distance",
                    x = np.geomspace(0.001, z_max, resolution),
                    pdf_func= Da, 
                    conditioned_y=None, 
                    dimension=1,
                    category="function",
                    create_new=create_new,
                )
        self.splineDa = splineDa
        angular_diameter_distance = njit(lambda z_: cubic_spline_interpolator(z_, splineDa[0], splineDa[1]))
        self.angular_diameter_distance = angular_diameter_distance

        # for angular diameter distance between two redshifts
        self.angular_diameter_distance_z1z2 = njit(lambda zl0, zs0: (angular_diameter_distance(zs0)*(1.+zs0) - angular_diameter_distance(zl0)*(1.+zl0))/(1.+zs0))

        # get differential co-moving volume interpolator
        resolution = self.c_n_i["differential_comoving_volume"]["resolution"]
        create_new = self.c_n_i["differential_comoving_volume"]["create_new"]
        splinedVcdz = interpolator_from_pickle(
            param_dict_given= dict(z_min=0.001, z_max=z_max, cosmology=self.cosmo, resolution=resolution), 
            directory=self.directory,
            sub_directory="differential_comoving_volume", 
            name="differential_comoving_volume",
            x = np.geomspace(0.001, z_max, resolution),
            pdf_func= lambda z_: self.cosmo.differential_comoving_volume(z_).value * 4 * np.pi,
            conditioned_y=None, 
            dimension=1,
            category="function",
            create_new=create_new,
        )
        self.splinedVcdz = splinedVcdz
        self.differential_comoving_volume = njit(lambda z_: cubic_spline_interpolator(z_, splinedVcdz[0], splinedVcdz[1]))

    @property
    def strong_lensing_optical_depth(self):
        """
        Function to compute the strong lensing optical depth.
        """
        return self._strong_lensing_optical_depth

    @strong_lensing_optical_depth.setter
    def strong_lensing_optical_depth(self, input_):
        if callable(input_):
            self._strong_lensing_optical_depth = input_
        else:
            # input can be function or spline interpolator array
            try:
                self._strong_lensing_optical_depth = njit(lambda z_: cubic_spline_interpolator(z_, input_[0], input_[1]))
            except:
                raise ValueError("strong_lensing_optical_depth must be a callable function or spline interpolator array.")

    @property
    def sample_velocity_dispersion(self):
        """
        Function to sample velocity dispersion from gengamma distribution
        """
        return self._sample_velocity_dispersion

    @sample_velocity_dispersion.setter
    def sample_velocity_dispersion(self, prior):
        try:
            self._sample_velocity_dispersion = getattr(self, prior)(size=None, zl=None, get_attribute=True)
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

    @property
    def available_velocity_dispersion_list_and_its_params(self):
        """
        Function to list all available velocity dispersion sampler and its parameters.
        """
        self._available_velocity_dispersion_list_and_its_params = dict(
            velocity_dispersion_haris=dict(a=2.32 / 2.67, c=2.67),
            velocity_dispersion_gengamma=dict(a=2.32 / 2.67, c=2.67),
            velocity_dispersion_bernardi=None,
            velocity_dispersion_ewoud=None,
        )

        return self._available_velocity_dispersion_list_and_its_params
    
    @property
    def available_axis_ratio_list(self):
        """
        Function to list all available axis ratio sampler.
        """
        self._available_axis_ratio_list = dict(
            axis_ratio_rayleigh=axis_ratio_rayleigh,
            axis_ratio_SIS=axis_ratio_SIS,
        )

        return self._available_axis_ratio_list
    
    @property
    def available_optical_depth_list(self):
        """
        Function to list all available optical depth sampler.
        """
        self._available_optical_depth_list = dict(
            optical_depth_SIS_haris=self.optical_depth_SIS_haris,
            optical_depth_SIS_hemanta=self.optical_depth_SIS_hemanta,
            optical_depth_SIE_hemanta=self.optical_depth_SIE_hemanta,
        )

        return self._available_optical_depth_list








