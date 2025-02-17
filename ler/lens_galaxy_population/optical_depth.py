import os
from numba import njit
from multiprocessing import Pool
import numpy as np
from scipy.integrate import quad
from scipy.stats import gengamma, gaussian_kde
from scipy.interpolate import CubicSpline
from astropy.cosmology import LambdaCDM
from tqdm import tqdm
from lenstronomy.LensModel.Solver.epl_shear_solver import caustics_epl_shear
from shapely.geometry import Polygon


from ..utils import  interpolator_from_pickle, cubic_spline_interpolator, inverse_transform_sampler, cubic_spline_interpolator2d_array, save_pickle, interpolator_pickle_path, FunctionConditioning, inverse_transform_sampler2d, pdf_cubic_spline_interpolator2d_array, normal_pdf, normal_pdf_2d

from .jit_functions import phi_cut_SIE, phi, phi_loc_bernardi, phi_q2_ellipticity_hemanta, lens_redshift_SDSS_catalogue_sis, axis_ratio_rayleigh_pdf

class OpticalDepth():
    """
        Class to calculate the optical depth, velocity dispersion and axis-ratio of a lens galaxy population.

        Parameters
        ----------
        npool : int, optional
            Number of processors to use for multiprocessing (default is 4).
        z_min : float, optional
            Minimum redshift of the lens galaxy population (default is 0.0).
        z_max : float, optional
            Maximum redshift of the lens galaxy population (default is 10.0).
        cosmology : astropy.cosmology, optional
            Cosmology object to use (default is FlatLambdaCDM with H0=70, Om0=0.3, Ode0=0.7).
        lens_type : str, optional
            Type of the lens galaxy. Must be one of ['sie_galaxy', 'epl_shear_galaxy', 'sis_galaxy'] (default is 'epl_shear_galaxy').
        lens_functions : dict, optional
            Dictionary with lens-related functions.
        lens_functions_params : dict, optional
            Dictionary with parameters for the lens-related functions.
        lens_param_samplers : dict, optional
            Dictionary of sampler functions for velocity dispersion and axis-ratio.
        lens_param_samplers_params : dict, optional
            Dictionary with parameters for the priors of the samplers.
        directory : str, optional
            Directory where the interpolators are saved (default is './interpolator_pickle').
            If True, creates a new interpolator (default is False).
        verbose : bool, optional
            If True, prints additional information during initialization (default is False).

        Raises
        ------
        ValueError
            If `lens_type` is not in ['sie_galaxy', 'epl_shear_galaxy', 'sis_galaxy'].
    """

    def __init__(self,
        npool=4,
        z_min=0.0,
        z_max=10.0,
        cosmology=None,
        lens_type="epl_shear_galaxy",
        lens_functions= None,
        lens_functions_params=None,
        lens_param_samplers=None,
        lens_param_samplers_params=None,
        directory="./interpolator_pickle",
        create_new_interpolator=False,
        verbose=False,
    ):
        
        print("\nInitializing OpticalDepth class\n")
        if lens_type not in ["sie_galaxy", "epl_shear_galaxy", "sis_galaxy"]:
            raise ValueError("lens_type not in ['sie_galaxy', 'epl_shear_galaxy', 'sis_galaxy']")
        self.lens_type = lens_type
        
        self.npool = npool
        self.z_min = z_min
        self.z_max = z_max
        self.cosmo = cosmology if cosmology else LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        self.directory = directory

        # decision dict fot creating interpolator
        self.initialize_decision_dictionary(create_new_interpolator);

        # dealing with prior functions and categorization
        # considering user input samplers and functions
        self.lens_param_samplers, self.lens_param_samplers_params, self.lens_functions, self.lens_functions_params = self.lens_priors_categorization(lens_type, lens_param_samplers, lens_param_samplers_params, lens_functions, lens_functions_params)

        # considering user input samplers and functions
        self.initialize_all_lens_samplers_and_functions(lens_param_samplers, lens_param_samplers_params, lens_functions, lens_functions_params);

        # setting up astropy cosmology functions
        self.create_lookup_table_fuction();

        print(f"\n{self.lens_functions['cross_section']}\n")
        if self.lens_functions['cross_section'] == "interpolated_cross_section_function":
            self.create_cross_section_function();
            self.cross_section = getattr(self, "interpolated_cross_section_function")
        else:
            self.cross_section = getattr(self, self.lens_functions['cross_section'])
            self.nbrs=None
            self.values=None
            self.cross_section_spline=None
            self.sis_area_array=None

        # lens sampler initialization
        self.velocity_dispersion = self.lens_param_samplers['velocity_dispersion']
        self.axis_ratio = self.lens_param_samplers['axis_ratio']
        self.axis_rotation_angle = self.lens_param_samplers['axis_rotation_angle']
        self.density_profile_slope = self.lens_param_samplers['density_profile_slope']
        self.external_shear = self.lens_param_samplers['external_shear']
        self.lens_redshift = self.lens_param_samplers['lens_redshift']

        # lens function initialization
        self.optical_depth = self.lens_functions['optical_depth']

    ##################
    # Initialization #
    ##################
    def lens_priors_categorization(
        self, lens_type, lens_priors=None, lens_priors_params=None, lens_functions=None, lens_functions_params=None,
    ):
        """
            Function to categorize the lens priors/samplers

            Parameters
            ----------
                lens_type : `str`
                    lens type
                    e.g. 'epl_shear_galaxy' for elliptical power-law galaxy
                lens_priors : `dict`
                    dictionary of priors
                lens_priors_params : `dict`
                    dictionary of priors parameters
                lens_functions : `dict`
                    dictionary of lens functions

            Returns
            -------
            lens_priors_ : `dict`
                dictionary of priors
            lens_priors_params_ : `dict`
                dictionary of priors parameters
            lens_sampler_names_ : `dict`
                dictionary of sampler names
            lens_functions_ : `dict`
                dictionary of lens functions
        """

        if lens_type == "epl_shear_galaxy":
            lens_priors_ = dict(
                source_redshift_sl="strongly_lensed_source_redshifts",
                lens_redshift="lens_redshift_SDSS_catalogue_numerical",
                velocity_dispersion="velocity_dispersion_ewoud",
                axis_ratio="axis_ratio_rayleigh",
                axis_rotation_angle="axis_rotation_angle_uniform",
                external_shear="external_shear_normal",
                density_profile_slope="density_profile_slope_normal",
            )
            lens_priors_params_ = dict(
                source_redshift_sl=None,
                lens_redshift=None,
                velocity_dispersion=dict(sigma_min=10., sigma_max=420., alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78),
                axis_ratio=dict(q_min=0.2, q_max=1.),
                axis_rotation_angle=dict(phi_min=0.0, phi_max=2 * np.pi),
                external_shear=dict(mean=0., std=0.05),
                density_profile_slope=dict(mean=1.99, std=0.149),
            )
            lens_functions_ = dict(
                param_sampler_type="sample_all_routine_sie_sl",
                strong_lensing_condition="rjs_with_cross_section_sie_feixu",
                optical_depth="optical_depth_numerical",
                cross_section="interpolated_cross_section_function",
                galaxy_density_distribution="galaxy_density_distribution_oguri",
            )
            lens_functions_params_ = dict(
                strong_lensing_condition=None,
                optical_depth=None,
                param_sampler_type=dict(fast_sampler=True),
                cross_section=None,
                galaxy_density_distribution=None,
            )
        elif lens_type == "sie_galaxy":
            lens_priors_ = dict(
                source_redshift_sl="strongly_lensed_source_redshifts",
                lens_redshift="lens_redshift_SDSS_sie_numerical",
                velocity_dispersion="velocity_dispersion_ewoud",
                axis_ratio="axis_ratio_rayleigh",
                axis_rotation_angle="axis_rotation_angle_uniform",
                external_shear=njit( lambda size: np.zeros((2,size)) ),
                density_profile_slope=njit(lambda size: np.ones(size)*2.0),
            )
            lens_priors_params_ = dict(
                source_redshift_sl=None,
                lens_redshift=None,
                velocity_dispersion=dict(sigma_min=10., sigma_max=420., alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78),
                axis_ratio=dict(q_min=0.2, q_max=1.),
                axis_rotation_angle=dict(phi_min=0.0, phi_max=2 * np.pi),
                external_shear=None,
                density_profile_slope=None,
            )
            lens_functions_ = dict(
                strong_lensing_condition="rjs_with_cross_section_sie_feixu",
                optical_depth="optical_depth_sie_hemanta",
                param_sampler_type="sample_all_routine_sie_sl",
                cross_section="cross_section_sie_feixu",
            )
            lens_functions_params_ = dict(
                strong_lensing_condition=None,
                optical_depth=dict(interpolated_cross_section=True),
                param_sampler_type=None,
                cross_section=None,
            )
        else:
            raise ValueError("lens_type not recognized")

        # update the priors if input is given
        if lens_priors:
            lens_priors_.update(lens_priors)
        if lens_priors_params:
            lens_priors_params_.update(lens_priors_params)
        if lens_functions:
            lens_functions_.update(lens_functions)
        if lens_functions_params:
            lens_functions_params_.update(lens_functions_params)

        return(lens_priors_, lens_priors_params_, lens_functions_, lens_functions_params_)

    def initialize_decision_dictionary(self, create_new_interpolator):
        """
        Function to initialize decision dictionary for creating interpolator

        Parameters
        ----------
        create_new_interpolator : `dict` or `bool`
            dictionary to create new interpolator for velocity dispersion and optical depth.
        """

        # initialize the interpolator's parameters
        self.create_new_interpolator = dict(
            velocity_dispersion=dict(create_new=False, resolution=500),
            axis_ratio=dict(create_new=False, resolution=500),
            lens_redshift=dict(create_new=False, resolution=50),
            optical_depth=dict(create_new=False, resolution=48),
            comoving_distance=dict(create_new=False, resolution=500),
            angular_diameter_distance=dict(create_new=False, resolution=500),
            differential_comoving_volume=dict(create_new=False, resolution=500),
            cross_section=dict(create_new=False, resolution=1000000),
            density_profile_slope=dict(create_new=False, resolution=100),
            lens_parameters_kde_sl=dict(create_new=False, resolution=5000),
        )
        if isinstance(create_new_interpolator, dict):
            self.create_new_interpolator.update(create_new_interpolator)
        # if create_new_interpolator is True, create new interpolator for all
        elif create_new_interpolator is True:
            for key in self.create_new_interpolator.keys():
                self.create_new_interpolator[key]['create_new'] = True
        
    def initialize_all_lens_samplers_and_functions(self, lens_param_samplers, lens_param_samplers_params, lens_functions, lens_functions_params):
        """
        Function to initialize velocity dispersion sampler with it's settings. The reason I am seperating this from lens_priors_categorization is only a specific parameters needs special attention.

        Parameters
        ----------
        lens_param_samplers : `str` or `function`
            sampler name or function
        lens_param_samplers_params : `dict`
            sampler parameters
        lens_functions : `str` or `function`
            lens function name or function
        lens_functions_params : `dict`
            lens function parameters
        """

        sampler_prior_names = ['velocity_dispersion', 'axis_ratio', 'axis_rotation_angle', 'density_profile_slope', 'external_shear', 'lens_redshift']

        # if there is user input function, update the sampler priors
        for name in sampler_prior_names:
            # the following is already taken care in lens_priors_categorization
            # if (lens_param_samplers is not None) and (name in lens_param_samplers):
            #     self.lens_param_samplers[name] = lens_param_samplers[name]
            # if (lens_param_samplers_params is not None) and (name in lens_param_samplers_params):
            #     if self.lens_param_samplers_params[name] is None: # not a dictionary
            #         self.lens_param_samplers_params[name] = lens_param_samplers_params[name]
            #     else:  # if there is user input lens_param_samplers_params
            #         self.lens_param_samplers_params[name].update(lens_param_samplers_params[name])

            # if there is user input function name, update the sampler priors
            if isinstance(lens_param_samplers[name], str):
                # setting up sample prior params
                sampler_name = lens_param_samplers[name]  # e.g. 'axis_ratio_padilla_strauss'
                dict_ = self.available_lens_prior_list_and_its_params[name]  # e.g. {'axis_ratio_padilla_strauss': {'q_min': 0.2, 'q_max': 1.0}, ....}
                if sampler_name in dict_:
                    param_dict = dict_[sampler_name]  # {'q_min': 0.2, 'q_max': 1.0}
                    if lens_param_samplers_params[name] is None:  # not a dictionary
                        self.lens_param_samplers_params[name] = param_dict
                    else:  # if there is user input lens_param_samplers_params
                        self.lens_param_samplers_params[name] = param_dict.update(lens_param_samplers_params[name])  # user inputs will override the default values
                else:
                    raise ValueError(f"{name} sampler {sampler_name} not available.\n Available {name} samplers and its parameters are: {dict_[name]}")
            elif not callable(lens_param_samplers[name]):
                raise ValueError(f"Given {name} sampler should be either a string name of available sampler or a function")

                
        lens_function_names = ['optical_depth', 'cross_section']

        # if there is user input function, update the sampler priors
        for name in lens_function_names:
            if isinstance(lens_functions[name], str):
                # setting up sample prior params
                function_name = lens_functions[name]  
                dict_ = self.available_lens_functions_and_its_params[name]  
                if function_name in dict_:
                    param_dict = dict_[function_name]  
                    if lens_functions_params[name] is None:  # not a dictionary
                        self.lens_functions_params[name] = param_dict
                    else:  # if there is user input lens_functions_params
                        self.lens_functions_params[name] = param_dict.update(lens_functions_params[name])  # user inputs will override the default values
                else:
                    raise ValueError(f"{name} function {function_name} not available.\n Available {name} functions and its parameters are: {dict_[name]}")
            elif not callable(lens_functions[name]):
                raise ValueError(f"Given {name} function should be either a string name of available function or a function")
                
        # setting up multiprocessor function for lens redshift
        if 'epl_shear_galaxy':
            if self.lens_param_samplers['velocity_dispersion'] == 'velocity_dispersion_ewoud':
                from .mp import lens_redshift_epl_shear2_mp
                self._helper_sl_disribution_mp = lens_redshift_epl_shear2_mp
            else:
                from .mp import lens_redshift_epl_shear1_mp
                self._helper_sl_disribution_mp = lens_redshift_epl_shear1_mp
        elif 'sie_galaxy':
            if self.lens_param_samplers['velocity_dispersion'] == 'velocity_dispersion_ewoud':
                if self.lens_functions_params['cross_section']['interpolated_cross_section_function']:
                    from .mp import lens_redshift_sie4_mp
                    self._helper_sl_disribution_mp = lens_redshift_sie4_mp
                else:
                    from .mp import lens_redshift_sie2_mp
                    self._helper_sl_disribution_mp = lens_redshift_sie2_mp
            else:
                if self.lens_functions_params['cross_section']['interpolated_cross_section_function']:
                    from .mp import lens_redshift_sie3_mp
                    self._helper_sl_disribution_mp = lens_redshift_sie3_mp
                else:
                    from .mp import lens_redshift_sie1_mp
                    self._helper_sl_disribution_mp = lens_redshift_sie1_mp
        elif 'sis_galaxy':
            if self.lens_param_samplers['velocity_dispersion'] == 'velocity_dispersion_ewoud':
                from .mp import lens_redshift_sis2_mp
                self._helper_sl_disribution_mp = lens_redshift_sis2_mp
            else:
                from .mp import lens_redshift_sis1_mp
                self._helper_sl_disribution_mp = lens_redshift_sis1_mp


    #####################
    # Sampler functions #
    #####################
    def axis_ratio_rayleigh(self, size, sigma, get_attribute=False, **kwargs):
        """
        Function to sample axis ratio from rayleigh distribution with given velocity dispersion.

        Parameters
        ----------
        sigma : `float: array`
            velocity dispersion of the lens galaxy
        q_min, q_max : `float`
            minimum and maximum axis ratio
        get_attribute : `bool`
            if True, returns a function that can be used to sample axis ratio

        Returns
        -------
        q : `float: array`
            axis ratio of the lens galaxy

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_param_samplers=dict(axis_ratio="axis_ratio_rayleigh"))
        >>> print(od.axis_ratio(sigma=200.))
        """

        param_dict_given_ = self.lens_param_samplers_params["axis_ratio"].copy()
        param_dict_given_.update(self.lens_param_samplers_params["velocity_dispersion"].copy())
        param_dict_given_['name'] = "axis_ratio_rayleigh"
        param_dict_given_['resolution'] = self.create_new_interpolator["axis_ratio"]["resolution"]
        param_dict_given_.update(kwargs)

        q_array = np.linspace(
            param_dict_given_["q_min"], 
            param_dict_given_["q_max"],
            param_dict_given_["resolution"],
        )
        sigma_array = np.linspace(
            param_dict_given_["sigma_min"],
            param_dict_given_["sigma_max"],
            500,
        )

        q_pdf = lambda q, sigma: axis_ratio_rayleigh_pdf(
            q=q,
            sigma=sigma,
            q_min=param_dict_given_["q_min"],
            q_max=param_dict_given_["q_max"],
        )

        q_object = FunctionConditioning(
            function=q_pdf,
            x_array=q_array,
            conditioned_y_array=sigma_array,
            param_dict_given=param_dict_given_,
            directory=self.directory,
            sub_directory="axis_ratio",
            name="axis_ratio_rayleigh",
            create_new=self.create_new_interpolator["axis_ratio"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='rvs',
        )

        if get_attribute:
            return q_object
        else:
            return q_object(size, sigma)
        
    def axis_ratio_padilla_strauss(self, size=1000, get_attribute=False, **kwargs):
        """
        Function to sample axis ratio using Padilla and Strauss 2008 distribution for axis ratio

        Parameters
        ----------
        size : `int`
            sample size
        q_min, q_max : `float`
            minimum and maximum axis ratio
        get_attribute : `bool`
            if True, returns a function that can be used to sample axis ratio

        Returns
        -------
        q : `float: array`
            axis ratio of the lens galaxy

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_param_samplers=dict(axis_ratio="axis_ratio_padilla_strauss"))
        >>> print(od.axis_ratio(size=10))
        """

        param_dict_given_ = self.lens_param_samplers_params["axis_ratio"].copy()
        param_dict_given_['name'] = "axis_ratio_padilla_strauss"
        param_dict_given_['resolution'] = self.create_new_interpolator["axis_ratio"]["resolution"]
        param_dict_given_.update(kwargs)

        # Using Padilla and Strauss 2008 distribution for axis ratio
        q_array = np.array([0.04903276402927845, 0.09210526315789469, 0.13596491228070173, 0.20789473684210524, 0.2899703729522482, 0.3230132450331126, 0.35350877192982455, 0.37946148483792264, 0.4219298245614036, 0.4689525967235971, 0.5075026141512723, 0.5226472638550018, 0.5640350877192983, 0.6096491228070177, 0.6500000000000001, 0.6864848379226213, 0.7377192982456142, 0.7787295224817011, 0.8007581038689441, 0.822786685256187, 0.8668438480306729, 0.8973684210526317, 0.9254385964912283])
        pdf = np.array([0.04185262687135349, 0.06114520695141845, 0.096997499638376, 0.1932510900336828, 0.39547914337673706, 0.49569751276216234, 0.6154609137685201, 0.7182049959882812, 0.920153741243567, 1.1573982157399754, 1.3353263628106684, 1.413149656448315, 1.5790713532948977, 1.7280185150744938, 1.8132994441344819, 1.8365803753840484, 1.8178662203211204, 1.748929843583365, 1.688182592496342, 1.6274353414093188, 1.4948487090314488, 1.402785526832393, 1.321844068356993])

        q_object = FunctionConditioning(
            function=pdf,
            x_array=q_array,
            conditioned_y_array=None,
            param_dict_given=param_dict_given_,
            directory=self.directory,
            sub_directory="axis_ratio",
            name="axis_ratio_padilla_strauss",
            create_new=self.create_new_interpolator["axis_ratio"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='rvs',
        )

        if get_attribute:
            return q_object
        else:
            return q_object(size)
        
    def lens_redshift_SDSS_catalogue_numerical(self, size=1000, zs=None, get_attribute=False, **kwargs):
        """
        Function to sample lens redshifts, conditioned on the lens being strongly lensed

        Parameters
        ----------
        size : `int`
            sample size
        zs : `float`
            source redshifts
        get_attribute : `bool`
            if True, returns a function that can be used to sample lens redshifts

        Returns
        -------
        zs : `float: array`
            lens redshifts

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_param_samplers=dict(lens_redshift="lens_redshift_SDSS_catalogue_numerical"))
        >>> print(od.lens_redshift(size=10, zs=1.0))
        """
    
        param_dict_given_ = self.lens_param_samplers_params["lens_redshift"].copy() if self.lens_param_samplers_params["lens_redshift"] is not None else {}
        try:
            param_dict_given_.update(self.velocity_dispersion.info)
        except:
            param_dict_given_['z_min'] = self.z_min
            param_dict_given_['z_max'] = self.z_max
            param_dict_given_['cosmology'] = self.cosmo
        param_dict_given_['name'] = "lens_redshift_numerical"+self.lens_type
        param_dict_given_['resolution'] = self.create_new_interpolator["lens_redshift"]["resolution"]
        param_dict_given_.update(kwargs)

        print("Numerically solving the lens redshift distribution...")
        zs_resolution = self.create_new_interpolator["optical_depth"]["resolution"]
        zs_array = np.geomspace(self.z_min+0.001, self.z_max, zs_resolution) if self.z_min==0 else np.geomspace(self.z_min, self.z_max, zs_resolution)

        zl_scaled2d = []
        for i, zs_ in enumerate(zs_array):
            buffer_ = np.linspace(0.0001, zs_-0.0001, param_dict_given_['resolution'])
            zl_scaled2d.append(buffer_/zs_)
        zl_scaled2d = np.array(zl_scaled2d)

        _, it_exist = interpolator_pickle_path(
            param_dict_given=param_dict_given_,
            directory=self.directory,
            sub_directory="lens_redshift",
            interpolator_name=param_dict_given_['name'],
        )

        create_new = self.create_new_interpolator["lens_redshift"]["create_new"]
        if not it_exist or create_new:
            number_density = self.lens_redshift_multiprocessing(zl_scaled2d, zs_array)
            # Adding zero at the first element of each row
            zl_scaled2d = np.hstack((np.zeros((zl_scaled2d.shape[0], 1)), zl_scaled2d))
            number_density = np.hstack((np.zeros((number_density.shape[0], 1)), number_density))
            # Adding one at the last element of each row of zl_scaled2d
            zl_scaled2d = np.hstack((zl_scaled2d, np.ones((zl_scaled2d.shape[0], 1))))
            # Adding zero at the last element of each row of density
            number_density = np.hstack((number_density, np.zeros((number_density.shape[0], 1))))
        else:
            number_density=None

        zl_object = FunctionConditioning(
            function=number_density,
            x_array=zl_scaled2d,
            conditioned_y_array=zs_array,
            param_dict_given=param_dict_given_,
            directory=self.directory,
            sub_directory="lens_redshift",
            name=param_dict_given_['name'],
            create_new=create_new,
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='rvs',
        )

        # un-scaled lens redshift
        cdf_values = zl_object.cdf_values
        x_array = zl_object.x_array
        y_array = zl_object.conditioned_y_array
        function_spline = zl_object.function_spline
        pdf_norm_const = zl_object.pdf_norm_const

        zl_object.function = njit(lambda x, y: cubic_spline_interpolator2d_array(x/y, y, function_spline, x_array, y_array))

        zl_object.pdf = njit(lambda x, y: pdf_cubic_spline_interpolator2d_array(x/y, y, pdf_norm_const, function_spline, x_array, y_array)) 

        zl_object.rvs = njit(lambda size, y: inverse_transform_sampler2d(size, y, cdf_values, x_array, y_array)*y)

        if get_attribute:
            return zl_object
        else:
            return zl_object(size, zs)
        
    def intrinsic_lens_redshift(self, size=1000, get_attribute=False, **kwargs):
        """
        Function to sample lens redshifts, conditioned on the lens being strongly lensed

        Parameters
        ----------
        size : `int`
            sample size
        zs : `float`
            source redshifts
        get_attribute : `bool`
            if True, returns a function that can be used to sample lens redshifts

        Returns
        -------
        zs : `float: array`
            lens redshifts
        """
    
        param_dict_given_ = self.lens_param_samplers_params["lens_redshift"].copy() if self.lens_param_samplers_params["lens_redshift"] is not None else {}
        try:
            param_dict_given_.update(self.velocity_dispersion.info)
        except:
            param_dict_given_['z_min'] = self.z_min
            param_dict_given_['z_max'] = self.z_max
            param_dict_given_['cosmology'] = self.cosmo
        param_dict_given_['name'] = "intrinsic_lens_redshift"
        param_dict_given_['resolution'] = 500

        _, it_exist = interpolator_pickle_path(
            param_dict_given=param_dict_given_,
            directory=self.directory,
            sub_directory="lens_redshift",
            interpolator_name=param_dict_given_['name'],
        )

        create_new = self.create_new_interpolator["lens_redshift"]["create_new"]
        if not it_exist or create_new:
            z_min = 0.001 if self.z_min==0 else self.z_min
            z_max = self.z_max
            zl_array = np.linspace(z_min, z_max, param_dict_given_['resolution']) 
            integrand = lambda sigma, z: self.velocity_dispersion.function(np.array([sigma]), np.array([z]))[0]
            integral = [quad(integrand, 10., 420., args=(z))[0] for z in zl_array]
        else:
            zl_array = None
            integral=None

        zl_object = FunctionConditioning(
            function=integral,
            x_array=zl_array,
            param_dict_given=param_dict_given_,
            directory=self.directory,
            sub_directory="lens_redshift",
            name=param_dict_given_['name'],
            create_new=create_new,
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='rvs',
        )

        if get_attribute:
            return zl_object
        else:
            return zl_object(size)
        
    def lens_redshift_multiprocessing(self, zl_scaled2d, zs1d):
        """
            Compute the lens redshift distribution using multiprocessing.
            
            Parameters
            ----------
            zl_scaled2d : array_like
                2D array of lens redshifts, scaled by the source redshift.
            zs1d : array_like
                1D array of source redshifts.
            zl_distribution_name : str
                Name of the lens redshift distribution to compute.
            
            Returns
            -------
            density_array : array_like
                2D array of the lens redshift distribution.
        """
   
        sigma_args = self.lens_param_samplers_params["velocity_dispersion"]
        sigma_args = [
            sigma_args["sigma_min"],
            sigma_args["sigma_max"],
            sigma_args["alpha"],
            sigma_args["beta"],
            sigma_args["phistar"],
            sigma_args["sigmastar"],
        ]
        # 3. q; q_rvs = self.axis_ratio.rvs
        q_args = [
            np.array(self.axis_ratio.cdf_values),
            np.array(self.axis_ratio.x_array),
            np.array(self.axis_ratio.conditioned_y_array),
        ]

        # 4. Da; Da_function = self.angular_diameter_distance.function
        Da_args = [
            np.array(self.angular_diameter_distance.function_spline),
            np.array(self.angular_diameter_distance.x_array),
        ]

        # 5. dVcdz; dVcdz_function = self.differential_comoving_volume.function
        dVcdz_args = [
            np.array(self.differential_comoving_volume.function_spline),
            np.array(self.differential_comoving_volume.x_array),
        ]

        # 6. idx; index to order the results
        idx = np.arange(len(zs1d))

        # 7. cross section
        # cs_fn = self.cross_section_function
        cs_args = [
            self.nbrs,
            np.array(self.values),
            np.array(self.cross_section_spline),
            np.array(self.sis_area_array),
        ]

        # 8. phi; axis_rotation_angle = self.axis_rotation_angle.rvs
        phi_args = self.lens_param_samplers_params["axis_rotation_angle"]
        phi_args = [
            phi_args["phi_min"],
            phi_args["phi_max"],
        ]

        # 9. shear; external_shear = self.external_shear.rvs
        shear_args = self.lens_param_samplers_params["external_shear"]
        shear_args = [
            shear_args['mean'],
            shear_args['std'],
        ]

        # 10. slope; density_profile_slope = self.density_profile_slope.rvs
        slope_args = self.lens_param_samplers_params["density_profile_slope"]
        slope_args = [
            slope_args['mean'],
            slope_args['std'],
        ]

        input_params = np.array([(zs1d[i], zl_scaled2d[i], sigma_args, q_args, Da_args, dVcdz_args, idx[i], cs_args, phi_args, shear_args, slope_args) for i in range(len(zs1d))], dtype=object)

        print("Computing lens redshift distribution with multiprocessing...")
        # with tqdm
        mp_fn = self._helper_sl_disribution_mp
        density_array = np.zeros_like(zl_scaled2d)
        with Pool(processes=self.npool) as pool:            
            for result in tqdm(
                pool.imap_unordered(mp_fn, input_params),
                total=len(zs1d),
                ncols=100,
                disable=False,
            ):
                # print(result)
                (
                    iter_i,
                    density_,
                ) = result

                density_array[iter_i] = density_
                
        # # without multiprocessing, just for loop with tqdm
        # density_array = np.zeros_like(zl_scaled2d)
        # for i in range(len(zs1d)):
        #     iter_i, density_ = mp_fn(input_params[i])
        #     density_array[iter_i] = density_
        # except:
        #     return input_params

        return density_array

    def axis_rotation_angle_uniform(self, size, get_attribute=False, **kwargs):

        phi_min = self.lens_param_samplers_params["axis_rotation_angle"]["phi_min"]
        phi_max = self.lens_param_samplers_params["axis_rotation_angle"]["phi_max"]
        phi_rvs = njit(lambda size: np.random.uniform(
            low=phi_min,
            high=phi_max,
            size=size,
        ))
        phi_pdf = lambda phi: 1/(phi_max-phi_min)

        phi_object =FunctionConditioning(
            function=None,
            x_array=None,
            create_rvs=phi_rvs,
            create_pdf=phi_pdf,
            callback='rvs',
        )

        if get_attribute:
            return phi_object
        else:
            return phi_rvs(size)
        
    def external_shear_normal(self, size, get_attribute=False, **kwargs):

        mean = self.lens_param_samplers_params["external_shear"]["mean"]
        std = self.lens_param_samplers_params["external_shear"]["std"]
        shear_rvs = njit(lambda size: np.random.normal(
            loc=mean,
            scale=std,
            size=(2,size),
        ))
        shear_pdf = njit(lambda shear1, shear2: normal_pdf_2d(
            x=shear1,
            y=shear2,
            mean_x=mean,
            mean_y=mean,
            std_x=std,
            std_y=std,
        ))

        shear_object =FunctionConditioning(
            function=None,
            x_array=None,
            create_rvs=shear_rvs,
            create_pdf=shear_pdf,
            callback='rvs',
        )

        if get_attribute:
            return shear_object
        else:
            return shear_rvs(size)
        
    def density_profile_slope_normal(self, size, get_attribute=False, **kwargs):

        mean = self.lens_param_samplers_params["density_profile_slope"]["mean"]
        std = self.lens_param_samplers_params["density_profile_slope"]["std"]
        slope_rvs = njit(lambda size: np.random.normal(
            loc=mean,
            scale=std,
            size=size,
        ))
        slope_pdf = njit(lambda slope: normal_pdf(
            x=slope,
            mean=mean,
            std=std,
        ))

        slope_object =FunctionConditioning(
            function=None,
            x_array=None,
            create_rvs=slope_rvs,
            create_pdf=slope_pdf,
            callback='rvs',
        )

        if get_attribute:
            return slope_object
        else:
            return slope_rvs(size)
    
    def lens_redshift_SDSS_catalogue_sis(self, size, zs, get_attribute=False, **kwargs):
        """
        Function to sample lens redshifts, conditioned on the lens being strongly lensed

        Parameters
        ----------
        zs : `float`
            source redshifts
        get_attribute : `bool`
            If True, returns a function that can be called with zs as input

        Returns
        -------
        zl : `float`
            lens redshifts

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> lens.lens_redshift_SDSS_catalogue_sis(zs=1.0)
        """

        # # Old way
        # splineDc = self.comoving_distance.function_spline  # spline coefficients for the comoving distance and redshifts
        # splineDcInv = self.comoving_distance.function_inverse_spline # spline coefficients for the redshifts and comoving distance
        # u = np.linspace(0, 1, 500)
        # cdf = (10 * u**3 - 15 * u**4 + 6 * u**5)  # See the integral of Eq. A7 of https://arxiv.org/pdf/1807.07062.pdf (cdf)
        # zs = np.array([zs]).reshape(-1)

        # New way
        param_dict_given_ = self.lens_param_samplers_params["lens_redshift"].copy() if self.lens_param_samplers_params["lens_redshift"] is not None else {}
        param_dict_given_.update(self.comoving_distance.info)
        param_dict_given_['name'] = "lens_redshift_numerical"+self.lens_type
        param_dict_given_['resolution'] = self.create_new_interpolator["lens_redshift"]["resolution"]
        param_dict_given_.update(kwargs)

        zl_resolution = self.create_new_interpolator["lens_redshift"]["resolution"]
        x_array = np.linspace(0., 1., zl_resolution)
        pdf_ = lambda x: 30* x**2 * (1-x)**2  # Haris et al. 2018 (A7)

        zl_object = FunctionConditioning(
            function=pdf_,
            x_array=x_array,
            param_dict_given=param_dict_given_,
            directory=self.directory,
            sub_directory="lens_redshift",
            name=param_dict_given_['name'],
            create_new=self.create_new_interpolator["lens_redshift"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='rvs',
        )

        cdf_values_zl = zl_object.cdf_values
        function_spline_zl = zl_object.function_spline
        pdf_norm_const = zl_object.pdf_norm_const
        x_array_zl = zl_object.x_array

        spline_Dc = self.comoving_distance.function_spline
        x_array_Dc = self.comoving_distance.x_array
        inverse_spline_Dc = self.comoving_distance.function_inverse_spline
        z_array_Dc = self.comoving_distance.z_array

        @njit
        def zl_rvs(size, zs):
            r = inverse_transform_sampler(size, cdf_values_zl, x_array_zl)
            zs_Dc = cubic_spline_interpolator(zs, spline_Dc, x_array_Dc)
            zl_Dc = zs_Dc * r
            return cubic_spline_interpolator(zl_Dc, inverse_spline_Dc, z_array_Dc)
        
        @njit
        def zl_pdf(zl, zs):
            r = zl / zs
            return cubic_spline_interpolator(r, function_spline_zl, x_array_zl)/pdf_norm_const
        
        @njit
        def zl_function(zl, zs):
            r = zl / zs
            return cubic_spline_interpolator(r, function_spline_zl, x_array_zl)
        
        zl_object.function = zl_function
        zl_object.pdf = zl_pdf
        zl_object.rvs = zl_rvs

        # lens redshifts
        #return self.Dc_to_z(lens_galaxy_Dc)
        if get_attribute:
            return zl_object
        else:
            return zl_object(size, zs)
        
         
    def velocity_dispersion_gengamma(self, size, get_attribute=False, **kwargs):
        """
        Function to sample velocity dispersion from gengamma distribution

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        a,c : `float`
            parameters of gengamma distribution
            refer to https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gengamma.html
        get_attribute : `bool`
            if True, returns a function that can be used to sample velocity dispersion
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. param = dict(a=2.32 / 2.67, c=2.67)

        Returns
        -------
        sigma : `numpy.ndarray` (1D array of floats)
            velocity dispersion of the lens galaxy

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_param_samplers=dict(velocity_dispersion="velocity_dispersion_gengamma"), lens_param_samplers_params=dict(velocity_dispersion=dict(a=2.32 / 2.67, c=2.67)))
        >>> print(od.velocity_dispersion(size=10))
        """

        param_dict_given_ = self.lens_param_samplers_params["velocity_dispersion"].copy()
        param_dict_given_['z_min'] = self.z_min
        param_dict_given_['z_max'] = self.z_max
        param_dict_given_['cosmology'] = self.cosmo
        param_dict_given_['name'] = "velocity_dispersion_gengamma"
        param_dict_given_['resolution'] = self.create_new_interpolator["velocity_dispersion"]["resolution"]
        # update with kwargs
        param_dict_given_.update(kwargs)
        # setting up inputs for the interpolator
        sigma_array = np.linspace(
            param_dict_given_['sigma_min'],
            param_dict_given_['sigma_max'],
            param_dict_given_["resolution"]
        )
        pdf_func_ = lambda sigma_: gengamma.pdf(
            sigma_/param_dict_given_['sigmastar'],
            a= param_dict_given_['alpha']/param_dict_given_['beta'],
            c= param_dict_given_['beta'],
        )   # gengamma pdf

        sigma_object = FunctionConditioning(
            function=pdf_func_,
            x_array=sigma_array,
            param_dict_given=param_dict_given_,
            directory=self.directory,
            sub_directory="velocity_dispersion",
            name=param_dict_given_['name'],
            create_new=self.create_new_interpolator["velocity_dispersion"]['create_new'],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='rvs',
        )
        
        if get_attribute:
            return sigma_object
        else:
            return sigma_object.rvs(size)
    
    def velocity_dispersion_bernardi(self, size, get_attribute=False, **kwargs):
        """
        Function to sample velocity dispersion from Bernardi et al. (2010). This uses inverse transform sampling.

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        get_attribute : `bool`
            if True, returns a function that can be used to sample velocity dispersion

        Returns
        -------
        sigma : `numpy.ndarray` (1D array of floats)
            velocity dispersion of the lens galaxy

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_param_samplers=dict(velocity_dispersion="velocity_dispersion_bernardi"))
        >>> print(od.velocity_dispersion(size=10))
        """

        param_dict_given_ = self.lens_param_samplers_params["velocity_dispersion"].copy()
        param_dict_given_['z_min'] = self.z_min
        param_dict_given_['z_max'] = self.z_max
        param_dict_given_['cosmology'] = self.cosmo
        param_dict_given_['name'] = "velocity_dispersion_bernardi"
        param_dict_given_['resolution'] = self.create_new_interpolator["velocity_dispersion"]["resolution"]
        # update with kwargs
        param_dict_given_.update(kwargs)
        # setting up inputs for the interpolator
        sigma_array = np.linspace(
            param_dict_given_['sigma_min'],
            param_dict_given_['sigma_max'],
            param_dict_given_["resolution"]
        )
        number_density_function = lambda sigma: phi_loc_bernardi(
            sigma, 
            alpha=self.lens_param_samplers_params["velocity_dispersion"]['alpha'],
            beta=self.lens_param_samplers_params["velocity_dispersion"]['beta'],
            phistar=self.lens_param_samplers_params["velocity_dispersion"]['phistar'],
            sigmastar=self.lens_param_samplers_params["velocity_dispersion"]['sigmastar'],
        )

        sigma_object = FunctionConditioning(
            function=number_density_function,
            x_array=sigma_array,
            param_dict_given=param_dict_given_,
            directory=self.directory,
            sub_directory="velocity_dispersion",
            name=param_dict_given_['name'],
            create_new=self.create_new_interpolator["velocity_dispersion"]['create_new'],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='rvs',
        )
        
        if get_attribute:
            return sigma_object
        else:
            return sigma_object.rvs(size)

    def velocity_dispersion_ewoud(self, size, zl, get_attribute=False, **kwargs):
        """
        Function to sample velocity dispersion (redshift dependent) from Wempe et al. (2022). This uses inverse transform sampling.

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        zl : `float`
            redshift of the lens galaxy
        get_attribute : `bool`
            if True, returns a function that can be used to sample velocity dispersion

        Returns
        -------
        sigma : `numpy.ndarray` (1D array of floats)
            velocity dispersion of the lens galaxy

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_param_samplers=dict(velocity_dispersion="velocity_dispersion_ewoud"))
        >>> print(od.velocity_dispersion(size=10, zl=0.5))
        """

        param_dict_given_ = self.lens_param_samplers_params["velocity_dispersion"].copy()
        param_dict_given_['z_min'] = self.z_min
        param_dict_given_['z_max'] = self.z_max
        param_dict_given_['cosmology'] = self.cosmo
        param_dict_given_['name'] = "velocity_dispersion_ewoud"
        param_dict_given_['resolution'] = self.create_new_interpolator["velocity_dispersion"]["resolution"]
        # update with kwargs
        param_dict_given_.update(kwargs)
        # setting up inputs for the interpolator
        sigma_array = np.linspace(
            param_dict_given_['sigma_min'],
            param_dict_given_['sigma_max'],
            param_dict_given_["resolution"]
        )
        number_density_function = lambda sigma, zl: phi(
            sigma, 
            zl,
            alpha=self.lens_param_samplers_params["velocity_dispersion"]['alpha'],
            beta=self.lens_param_samplers_params["velocity_dispersion"]['beta'],
            phistar=self.lens_param_samplers_params["velocity_dispersion"]['phistar'],
            sigmastar=self.lens_param_samplers_params["velocity_dispersion"]['sigmastar'],
        )

        zl_array = np.geomspace(self.z_min+0.001, self.z_max, 500) if self.z_min==0 else np.geomspace(self.z_min, self.z_max, 500)

        sigma_object = FunctionConditioning(
            function=number_density_function,
            x_array=sigma_array,
            conditioned_y_array=zl_array,
            param_dict_given=param_dict_given_,
            directory=self.directory,
            sub_directory="velocity_dispersion",
            name=param_dict_given_['name'],
            create_new=self.create_new_interpolator["velocity_dispersion"]['create_new'],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='rvs',
        )
        
        if get_attribute:
            return sigma_object
        else:
            return sigma_object.rvs(size, zl)

    ##################
    # Lens functions #
    ##################
    def optical_depth_numerical(self, zs, get_attribute=False, **kwargs):

        param_dict_given_ = self.lens_functions_params["optical_depth"] if self.lens_functions_params["optical_depth"] is not None else {}
        param_dict_given_['z_min'] = self.z_min
        param_dict_given_['z_max'] = self.z_max
        param_dict_given_['cosmology'] = self.cosmo
        param_dict_given_['name'] = "optical_depth_numerical"
        param_dict_given_['resolution'] = self.create_new_interpolator["optical_depth"]["resolution"]
        param_dict_given_['sigma_args'] = self.lens_param_samplers_params["velocity_dispersion"]
        param_dict_given_['axis_rotation_angle_args'] = self.lens_param_samplers_params["axis_rotation_angle"]
        param_dict_given_['external_shear_args'] = self.lens_param_samplers_params["external_shear"]
        param_dict_given_['density_profile_slope_args'] = self.lens_param_samplers_params["density_profile_slope"]
        param_dict_given_['axis_ratio_args'] = self.lens_param_samplers_params["axis_ratio"]
        param_dict_given_.update(kwargs)

        z_min = self.z_min if self.z_min>0. else 0.001
        z_max = self.z_max
        resolution = self.create_new_interpolator["optical_depth"]["resolution"]
        zs_array = np.geomspace(
            z_min, 
            z_max, 
            resolution
        )

        def tau(zs):
            # self.lens_redshift.function gives cross-section
            integrand = lambda zl_, zs_: self.lens_redshift.function(np.array([zl_]), np.array([zs_]))[0]
            integral = [quad(integrand, 0.0, z, args=(z))[0] for z in zs] 
            return integral
        
        tau_object = FunctionConditioning(
            function=tau,
            x_array=zs_array,
            conditioned_y_array=None,
            param_dict_given=param_dict_given_,
            directory=self.directory,
            sub_directory="optical_depth",
            name=param_dict_given_['name'],
            create_new=self.create_new_interpolator["optical_depth"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='function',
        )

        if get_attribute:
            return tau_object
        else:
            return tau(zs)
        
    def cross_section_sis(self, sigma, zl, zs, **kwargs):
        """
        Function to compute the SIS cross-section

        Parameters
        ----------
        sigma : `float`
            velocity dispersion of the lens galaxy
        zl : `float`
            redshift of the lens galaxy
        zs : `float`
            redshift of the source galaxy

        Returns
        -------
        cross_section : `float`
            SIS cross-section

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> print(od.cross_section_sis(sigma=200., zl=0.5, zs=1.0))
        """

        theta_E = self.compute_einstein_radii(sigma, zl, zs)
        return np.pi * theta_E**2
    
    def compute_einstein_radii(self, sigma, zl, zs):
        """
        Function to compute the Einstein radii of the lens galaxies

        Parameters
        ----------
        sigma : `float`
            velocity dispersion of the lens galaxy
        zl : `float`
            lens redshifts
        zs : `float`
            source redshifts

        Returns
        -------
        theta_E : `float`
            Einstein radii of the lens galaxies

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> sigma = 200.0
        >>> zl = 0.5
        >>> zs = 1.0
        >>> lens.compute_einstein_radii(sigma, zl, zs)
        """

        # Compute the angular diameter distances
        Ds = self.angular_diameter_distance(zs)
        Dls = self.angular_diameter_distance_z1z2(zl, zs)
        # Compute the Einstein radii
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / (Ds)
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        return theta_E

    def cross_section_sie_feixu(self, sigma, zl, zs, q, **kwargs):
        """
        Function to compute the SIE cross-section from Fei Xu et al. (2021)

        Parameters
        ----------
        sigma : `float`
            velocity dispersion of the lens galaxy
        zl : `float`
            redshift of the lens galaxy
        zs : `float`
            redshift of the source galaxy

        Returns
        -------
        cross_section : `float`
            SIE cross-section

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> print(od.cross_section_sie_feixu(sigma=200., zl=0.5, zs=1.0, q=1.0))
        """

        # phi_cut_SIE already includes the factor of pi
        return phi_cut_SIE(q)*self.cross_section_sis(sigma=sigma, zl=zl, zs=zs)
    
    def optical_depth_epl_shear_hemanta(self, zs, get_attribute=False, **kwargs):
        """
        Function to compute the strong lensing optical depth (EPL with shear). \n
        LambdaCDM(H0=70, Om0=0.3, Ode0=0.7) was used to derive the following equation. This is the analytic version of optical depth from z=0 to z=zs.

        Parameters
        ----------
        zs : `float`
            source redshifts

        Returns
        -------
        tau : `float`
            strong lensing optical depth

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> print(od.optical_depth_epl_shear_lambdacdm(zs=1.0))
        """

        zs_arr = np.array([1.00000000e-03, 1.21648395e-03, 1.47983320e-03, 1.80019333e-03, 2.18990629e-03, 2.66398586e-03, 3.24069604e-03, 3.94225471e-03, 4.79568958e-03, 5.83387940e-03, 7.09682065e-03, 8.63316841e-03, 1.05021108e-02, 1.27756492e-02, 1.55413722e-02, 1.89058298e-02, 2.29986385e-02, 2.79774746e-02, 3.40341488e-02, 4.14019958e-02, 5.03648633e-02, 6.12680478e-02, 7.45315967e-02, 9.06664911e-02, 1.10294331e-01, 1.34171284e-01, 1.63217213e-01, 1.98551120e-01, 2.41534250e-01, 2.93822538e-01, 3.57430402e-01, 4.34808347e-01, 5.28937375e-01, 6.43443826e-01, 7.82739087e-01, 9.52189535e-01, 1.15832329e+00, 1.40908169e+00, 1.71412525e+00, 2.08520586e+00, 2.53661946e+00, 3.08575685e+00, 3.75377368e+00, 4.56640543e+00, 5.55495891e+00, 6.75751836e+00, 8.22041261e+00, 1.00000000e+01])
        # zs and tau data points
        if self.lens_param_samplers["velocity_dispersion"] == "velocity_dispersion_ewoud":
            tau_arr = np.array([6.03160386e-13, 1.10059613e-12, 1.95432798e-12,3.51678900e-12, 6.33417074e-12, 1.14135959e-11,2.05172549e-11, 3.69566725e-11, 6.64560594e-11,1.19567171e-10, 2.15045832e-10, 3.86804139e-10,6.95544002e-10, 1.25051410e-09, 2.24678876e-09,4.03770124e-09, 7.22902834e-09, 1.29770147e-08,2.32749843e-08, 4.17109576e-08, 7.46815960e-08,1.33544615e-07, 2.38334487e-07, 4.24696339e-07,7.55269773e-07, 1.34275891e-06, 2.36388764e-06,4.15418141e-06, 7.25947183e-06, 1.25895566e-05,2.16379755e-05, 3.67860126e-05, 6.16981466e-05,1.01690398e-04, 1.64327884e-04, 2.59006953e-04,3.97439317e-04, 5.91890529e-04, 8.52861969e-04,1.18758691e-03, 1.59535092e-03, 2.07060366e-03,2.59787615e-03, 3.15410459e-03, 3.71873528e-03,4.27359478e-03, 4.80010742e-03, 5.29316594e-03])
        elif self.lens_param_samplers["velocity_dispersion"] == "velocity_dispersion_bernardi":
            tau_arr = np.array([6.00850294e-13, 1.08151426e-12, 1.94596619e-12,3.50298950e-12, 6.30402062e-12, 1.13442972e-11,2.04185899e-11, 3.67329162e-11, 6.62640787e-11,1.19216485e-10, 2.14427873e-10, 3.85652158e-10,6.93134949e-10, 1.24617524e-09, 2.23879004e-09,4.02030613e-09, 7.23573234e-09, 1.29811881e-08,2.32726057e-08, 4.16808742e-08, 7.45702352e-08,1.33237672e-07, 2.37571038e-07, 4.22869164e-07,7.50078130e-07, 1.32708456e-06, 2.33911764e-06,4.10269165e-06, 7.14815719e-06, 1.23628518e-05,2.11810176e-05, 3.58674556e-05, 5.96848327e-05,9.79654101e-05, 1.57538383e-04, 2.47616378e-04,3.79005699e-04, 5.63804261e-04, 8.13766486e-04,1.13922731e-03, 1.55230925e-03, 2.04780867e-03,2.62674567e-03, 3.28454887e-03, 4.00924305e-03,4.79080880e-03, 5.61454580e-03, 6.46916402e-03])
        elif self.lens_param_samplers["velocity_dispersion"] == "velocity_dispersion_choi":
            # print("\nChoi\n")
            tau_arr = np.array([4.08414474e-13, 7.35144978e-13, 1.32295988e-12, 2.38118175e-12, 4.28588534e-12, 7.71120932e-12, 1.38771914e-11, 2.49711498e-11, 4.48291558e-11, 8.06450938e-11, 1.45027675e-10, 2.60849771e-10, 4.69063665e-10, 8.43134806e-10, 1.51479934e-09, 2.71997503e-09, 4.87999147e-09, 8.75163578e-09, 1.56877772e-08, 2.81048883e-08, 5.02784547e-08, 8.98264456e-08, 1.60176352e-07, 2.85127183e-07, 5.08168496e-07, 8.98942430e-07, 1.58429969e-06, 2.77820454e-06, 4.84094964e-06, 8.37289340e-06, 1.43435176e-05, 2.42868938e-05, 4.04335816e-05, 6.63637616e-05, 1.06745981e-04, 1.67759099e-04, 2.56906334e-04, 3.82278437e-04, 5.51906916e-04, 7.72602939e-04, 1.04818537e-03, 1.38241558e-03, 1.77355736e-03, 2.21595390e-03, 2.70449762e-03, 3.23214166e-03, 3.78761428e-03, 4.36316348e-03])
        else:
            raise ValueError("velocity_dispersion should be one of velocity_dispersion_ewoud, velocity_dispersion_bernardi, velocity_dispersion_choi.\n Otherwise use optical_depth_numerical with lens_redshift_SDSS_catalogue_numerical.")

        # setting up input parameters for interpolation
        param_dict_given_ = self.lens_functions_params["optical_depth"] if self.lens_functions_params["optical_depth"] is not None else {}
        param_dict_given_['z_min'] = self.z_min
        param_dict_given_['z_max'] = self.z_max
        param_dict_given_['cosmology'] = self.cosmo
        param_dict_given_['name'] = "optical_depth_epl_shear_hemanta"
        param_dict_given_['resolution'] = self.create_new_interpolator["optical_depth"]["resolution"]
        param_dict_given_['sigma_args'] = self.lens_param_samplers_params["velocity_dispersion"]
        param_dict_given_['sigma_args']['name'] = self.lens_param_samplers["velocity_dispersion"]
        param_dict_given_['axis_rotation_angle_args'] = self.lens_param_samplers_params["axis_rotation_angle"]
        param_dict_given_['external_shear_args'] = self.lens_param_samplers_params["external_shear"]
        param_dict_given_['density_profile_slope_args'] = self.lens_param_samplers_params["density_profile_slope"]
        param_dict_given_['axis_ratio_args'] = self.lens_param_samplers_params["axis_ratio"]
        param_dict_given_.update(kwargs)

        tau_object = FunctionConditioning(
            function=tau_arr,
            x_array=zs_arr,
            conditioned_y_array=None,
            param_dict_given=param_dict_given_,
            directory=self.directory,
            sub_directory="optical_depth",
            name=param_dict_given_['name'],
            create_new=self.create_new_interpolator["optical_depth"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='function',
        )

        if get_attribute:
            return tau_object
        else:
            return tau_object.function(zs)

        
    def optical_depth_sie_hemanta(self, zs, get_attribute=False, **kwargs):
        """
        Function to compute the strong lensing optical depth (EPL with shear). \n
        LambdaCDM(H0=70, Om0=0.3, Ode0=0.7) was used to derive the following equation. This is the analytic version of optical depth from z=0 to z=zs.

        Parameters
        ----------
        zs : `float`
            source redshifts

        Returns
        -------
        tau : `float`
            strong lensing optical depth

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> print(od.optical_depth_epl_shear_lambdacdm(zs=1.0))
        """

        zs_arr = np.array([1.00000000e-03, 1.21648395e-03, 1.47983320e-03, 1.80019333e-03, 2.18990629e-03, 2.66398586e-03, 3.24069604e-03, 3.94225471e-03, 4.79568958e-03, 5.83387940e-03, 7.09682065e-03, 8.63316841e-03, 1.05021108e-02, 1.27756492e-02, 1.55413722e-02, 1.89058298e-02, 2.29986385e-02, 2.79774746e-02, 3.40341488e-02, 4.14019958e-02, 5.03648633e-02, 6.12680478e-02, 7.45315967e-02, 9.06664911e-02, 1.10294331e-01, 1.34171284e-01, 1.63217213e-01, 1.98551120e-01, 2.41534250e-01, 2.93822538e-01, 3.57430402e-01, 4.34808347e-01, 5.28937375e-01, 6.43443826e-01, 7.82739087e-01, 9.52189535e-01, 1.15832329e+00, 1.40908169e+00, 1.71412525e+00, 2.08520586e+00, 2.53661946e+00, 3.08575685e+00, 3.75377368e+00, 4.56640543e+00, 5.55495891e+00, 6.75751836e+00, 8.22041261e+00, 1.00000000e+01])
        # zs and tau data points
        if self.lens_param_samplers["velocity_dispersion"] == "velocity_dispersion_ewoud":
            tau_arr = np.array([4.54036644e-13, 8.17268345e-13, 1.47099928e-12,2.64758548e-12, 4.76546219e-12, 8.57629879e-12,1.54350154e-11, 2.77685962e-11, 4.99750460e-11,8.99132847e-11, 1.61762940e-10, 2.90914905e-10,5.23059702e-10, 9.40436875e-10, 1.69006517e-09,3.03640554e-09, 5.44428663e-09, 9.76944171e-09,1.75267955e-08, 3.14050085e-08, 5.62138789e-08,1.00536493e-07, 1.79457116e-07, 3.19758127e-07,5.70201751e-07, 1.01048160e-06, 1.78413346e-06,3.13567775e-06, 5.47906920e-06, 9.50476063e-06,1.63384307e-05, 2.77745755e-05, 4.64364701e-05,7.65441573e-05, 1.23597343e-04, 1.94901461e-04,2.99138988e-04, 4.45383678e-04, 6.41658926e-04,8.93600436e-04, 1.20324039e-03, 1.56217727e-03,1.95954919e-03, 2.38176031e-03, 2.80927989e-03,3.22777949e-03, 3.62813939e-03, 3.99877537e-03])
        elif self.lens_param_samplers["velocity_dispersion"] == "velocity_dispersion_bernardi":
            tau_arr = np.array([4.54063374e-13, 8.17279328e-13, 1.47114189e-12, 2.64715333e-12, 4.76477438e-12, 8.57443341e-12, 1.54282211e-11, 2.77619047e-11, 4.99557984e-11, 8.98575333e-11, 1.61640635e-10, 2.90702164e-10, 5.22577159e-10, 9.39402651e-10, 1.68784451e-09, 3.03139742e-09, 5.43477142e-09, 9.75051198e-09, 1.74810705e-08, 3.13081553e-08, 5.60075774e-08, 1.00062419e-07, 1.78479781e-07, 3.17642651e-07, 5.65649750e-07, 1.00075359e-06, 1.76376952e-06, 3.09304545e-06, 5.39042972e-06, 9.32256438e-06, 1.59713696e-05, 2.70505750e-05, 4.50540923e-05, 7.39462404e-05, 1.18936516e-04, 1.86949693e-04, 2.86271754e-04, 4.25919631e-04, 6.14957660e-04, 8.60665526e-04, 1.17141003e-03, 1.54454572e-03, 1.98089667e-03, 2.47619125e-03, 3.02253022e-03, 3.61140114e-03, 4.23206664e-03, 4.87465742e-03])
        elif self.lens_param_samplers["velocity_dispersion"] == "velocity_dispersion_choi":
            tau_arr = np.array([3.08366701e-13, 5.54949810e-13, 9.98838724e-13,1.79781610e-12, 3.23592802e-12, 5.82346322e-12,1.04792802e-11, 1.88535240e-11, 3.38875656e-11,6.09717917e-11, 1.09663271e-10, 1.97183649e-10,3.54528312e-10, 6.37217292e-10, 1.14494304e-09,2.05681955e-09, 3.68773233e-09, 6.61558072e-09,1.18637936e-08, 2.12438273e-08, 3.80072216e-08,6.78958282e-08, 1.21115645e-07, 2.15545918e-07,3.84052022e-07, 6.79353353e-07, 1.19734646e-06,2.09966521e-06, 3.65958234e-06, 6.32972798e-06,1.08436671e-05, 1.83628014e-05, 3.05685722e-05,5.01655770e-05, 8.07097438e-05, 1.26851057e-04,1.94231292e-04, 2.88923844e-04, 4.17255033e-04,5.84257923e-04, 7.94755792e-04, 1.04802273e-03,1.34395427e-03, 1.67941176e-03, 2.05010705e-03,2.44910251e-03, 2.86979677e-03, 3.30478671e-03])
        else:
            raise ValueError(
                f"velocity_dispersion must be one of 'velocity_dispersion_gengamma', 'velocity_dispersion_choi', 'velocity_dispersion_gaussian'.")
        
        # setting up input parameters for interpolation
        param_dict_given_ = self.lens_functions_params["optical_depth"] if self.lens_functions_params["optical_depth"] is not None else {}
        param_dict_given_['z_min'] = self.z_min
        param_dict_given_['z_max'] = self.z_max
        param_dict_given_['cosmology'] = self.cosmo
        param_dict_given_['name'] = "optical_depth_sie_hemanta"
        param_dict_given_['resolution'] = self.create_new_interpolator["optical_depth"]["resolution"]
        param_dict_given_['sigma_args'] = self.lens_param_samplers_params["velocity_dispersion"]
        param_dict_given_['sigma_args']['name'] = self.lens_param_samplers["velocity_dispersion"]
        param_dict_given_['axis_rotation_angle_args'] = self.lens_param_samplers_params["axis_rotation_angle"]
        param_dict_given_['external_shear_args'] = self.lens_param_samplers_params["external_shear"]
        param_dict_given_['density_profile_slope_args'] = self.lens_param_samplers_params["density_profile_slope"]
        param_dict_given_['axis_ratio_args'] = self.lens_param_samplers_params["axis_ratio"]
        param_dict_given_.update(kwargs)

        tau_object = FunctionConditioning(
            function=tau_arr,
            x_array=zs_arr,
            conditioned_y_array=None,
            param_dict_given=param_dict_given_,
            directory=self.directory,
            sub_directory="optical_depth",
            name=param_dict_given_['name'],
            create_new=self.create_new_interpolator["optical_depth"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='function',
        )

        if get_attribute:
            return tau_object
        else:
            return tau_object.function(zs)

    def optical_depth_sis_haris(self, zs, get_attribute=False, **kwargs):
        """
        Function to compute the strong lensing optical depth (SIS). \n
        LambdaCDM(H0=70, Om0=0.3, Ode0=0.7) was used to derive the following equation. This is the analytic version of optical depth from z=0 to z=zs.

        Parameters
        ----------
        zs : `float`
            source redshifts

        Returns
        -------
        tau : `float`
            strong lensing optical depth

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> print(od.optical_depth_sis_haris(zs=1.0))
        """

        # z to luminosity_distance (luminosity_distance) conversion
        #Dc = self.z_to_Dc(zs) * 1e-3  # 1e-3 converts Mpc to Gpc
        splineDc = self.splineDc

        @njit
        def tau(zs):
            Dc = cubic_spline_interpolator(zs, splineDc[0], splineDc[1])
            return (Dc * 1e-3 / 62.2) ** 3

        if get_attribute:
            return tau
        else:
            return tau(zs)
    
    #######################################
    # interpolated cross section function #
    #######################################
    def create_data_set(self, size=1000000):
        zs = np.random.uniform(0.01, 10.0, size-2)
        # inlcude z=0.01 and z=10.0
        zs = np.insert(zs, 0, 0.01)
        zs = np.append(zs, 10.0)

        zl = np.zeros(size)
        for i in range(size):
            zl[i] = np.random.uniform(0.001, zs[i]-0.001)

        sigma = np.random.uniform(self.lens_param_samplers_params["velocity_dispersion"]['sigma_min'], self.lens_param_samplers_params["velocity_dispersion"]['sigma_max'], size-2)
        # inlcude sigma_min and sigma_max
        sigma = np.insert(sigma, 0, self.lens_param_samplers_params["velocity_dispersion"]['sigma_min'])
        sigma = np.append(sigma, self.lens_param_samplers_params["velocity_dispersion"]['sigma_max'])

        q = np.random.uniform(self.lens_param_samplers_params["axis_ratio"]['q_min'], self.lens_param_samplers_params["axis_ratio"]['q_max'], size-2)
        # inlcude q_min and q_max
        q = np.insert(q, 0, self.lens_param_samplers_params["axis_ratio"]['q_min'])
        q = np.append(q, self.lens_param_samplers_params["axis_ratio"]['q_max'])

        phi = np.random.uniform(0.0, 2 * np.pi, size-2)
        # inlcude 0 and 2 * np.pi
        phi = np.insert(phi, 0, 0.0)
        phi = np.append(phi, 2 * np.pi)

        e1, e2 = phi_q2_ellipticity_hemanta(phi, q)
        gamma1, gamma2 = np.random.normal(loc=0, scale=0.05, size=(2, size))
        gamma = np.random.normal(loc=2.0, scale=0.2, size=size)
        # Da_zs = od_epl_ewoud.angular_diameter_distance(np.array([zs]))[0]
        # einstein radius 
        Ds = self.angular_diameter_distance.function(zs)
        Dls = self.angular_diameter_distance_z1z2(zl, zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / (Ds)
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        return theta_E, e1, e2, gamma, gamma1, gamma2, q
    
    def cross_section_caustic_area(self, theta_E, e1, e2, gamma, gamma1, gamma2):

        size = theta_E.shape[0]
        area_list = np.zeros(size)
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
                area_list[i] = (Polygon(caustic_double_points.T).area)
            
        return np.array(area_list)
    
    def cross_section_caustic_area_unit(self, e1, e2, gamma, gamma1, gamma2):

        size = e1.shape[0]
        area_list = np.zeros(size)
        for i in range(size):
            kwargs_lens = [
                {
                    "theta_E": 1.0,
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
                area_list[i] = (Polygon(caustic_double_points.T).area)
            
        return np.array(area_list)
    
    def create_cross_section_function(self):
        
        param_dict_given = dict(
            z_min=self.z_min,
            z_max=self.z_max,
            cosmology=self.cosmo,
            resolution=self.create_new_interpolator["cross_section"]["resolution"],
            sigma_min=self.lens_param_samplers_params["velocity_dispersion"]['sigma_min'],
            sigma_max=self.lens_param_samplers_params["velocity_dispersion"]['sigma_max'],
            q_min = self.lens_param_samplers_params["axis_ratio"]['q_min'],
            q_max = self.lens_param_samplers_params["axis_ratio"]['q_max'],
        )

        file_path, it_exist = interpolator_pickle_path(
            param_dict_given=param_dict_given,
            directory=self.directory,
            sub_directory="cross_section_function",
            interpolator_name="cross_section_function",
        )

        if (self.create_new_interpolator["cross_section"]["create_new"] is True) or (it_exist is False):
            
            print(f"Interpolated cross section function will be created at {file_path}")
            print("Creating interpolated cross section function...")

            size_interpolation=1000
            theta_E, e1, e2, gamma, gamma1, gamma2, _ = self.create_data_set(size_interpolation)
            # sort
            idx = np.argsort(theta_E)
            theta_E, e1, e2, gamma, gamma1, gamma2 = theta_E[idx], e1[idx], e2[idx], gamma[idx], gamma1[idx], gamma2[idx]
            values1 = self.cross_section_caustic_area(theta_E, e1, e2, gamma, gamma1, gamma2)
            values2 = self.cross_section_caustic_area_unit(e1, e2, gamma, gamma1, gamma2)

            x = np.pi*theta_E**2

            idx_cs = (values2 > 0.0) & (values1 > 0.0) & (x > 0.0)
            value3 = values1[idx_cs]/values2[idx_cs]
            x = x[idx_cs]
            spline_coeff = CubicSpline(x, value3).c
            
            # K nearest neighbors
            size = self.create_new_interpolator["cross_section"]["resolution"]
            iter_ = np.arange(size)
            _, e1, e2, gamma, gamma1, gamma2, _ = self.create_data_set(size)
            points = np.column_stack((e1, e2, gamma, gamma1, gamma2, iter_))
            from .mp import cross_section_unit_mp
            area_list = np.ones(size)
            with Pool(processes=8) as pool:            
                for result in tqdm(
                    pool.imap_unordered(cross_section_unit_mp, points),
                    total=size,
                    ncols=100,
                    disable=False,
                ):
                    (
                        iter_i,
                        area,
                    ) = result
                    area_list[int(iter_i)] = area

            values = np.array(area_list)
            points = points[:, :-1]

            from sklearn.neighbors import NearestNeighbors

            # Fit nearest neighbor model
            nbrs = NearestNeighbors(n_neighbors=10).fit(points)
            
            #print(mean_squared_error(values1, y_pred))
            
            # save the cross_section_function as pickle
            save_pickle(file_path, [nbrs,values,spline_coeff,x])
            self.nbrs = nbrs
            self.values = values
            self.cross_section_spline = spline_coeff
            self.sis_area_array = x

        else:
            from ..utils import load_pickle
            print(f"Interpolated cross section function loaded from {file_path}")
            load_ = load_pickle(file_path)
            self.nbrs, self.values, self.cross_section_spline, self.sis_area_array = load_

    def interpolated_cross_section_function(self, theta_E, e1, e2, gamma, gamma1, gamma2, get_attribute=False, **kwargs):
        """
        Function to compute the cross-section correction factor
        """
        nbrs = self.nbrs
        values = self.values
        
        points_ = np.array([e1, e2, gamma, gamma1, gamma2]).T
        _, indices = nbrs.kneighbors(points_)
        new_values_nn1 = np.mean(values[indices], axis=1)

        theta_E_correction = cubic_spline_interpolator(np.pi*theta_E**2, self.cross_section_spline, self.sis_area_array)

        return new_values_nn1*theta_E_correction


    ##########################
    # Cosmological functions #
    ##########################
    def create_lookup_table_fuction(self):
        """
        Functions to create lookup tables
        1. Redshift to co-moving distance.
        2. Co-moving distance to redshift.
        3. Redshift to angular diameter distance.

        """

        z_max = self.z_max
        z_min = 0.001 if self.z_min == 0. else self.z_min

        # for co-moving distance
        Dc = self.cosmo.comoving_distance(zs).value  # co-moving distance in Mpc
        resolution = self.create_new_interpolator["comoving_distance"]["resolution"]
        create_new = self.create_new_interpolator["comoving_distance"]["create_new"]
        zs = np.geomspace(z_min, z_max, resolution)
        self.comoving_distance = FunctionConditioning(
            function=Dc,
            x_array=zs,
            conditioned_y_array=None,
            param_dict_given=dict(z_min=z_min, z_max=z_max, cosmology=self.cosmo, resolution=resolution, details="comoving_distance from astropy.cosmology"),
            directory=self.directory,
            sub_directory="comoving_distance",
            name="comoving_distance",
            create_new=create_new,
            create_function_inverse=True,
            create_function=True,
            create_pdf=False,
            create_rvs=False,
            callback='function',
        )
        self.comoving_distance.__doc__ = """
        Redshift to comoving distance conversion.

        Parameters
        ----------
        zs : `numpy.ndarray` or `float`
            Source redshifts

        Returns
        ----------
        comoving_distance : `numpy.ndarray`
            comoving distance in Mpc

        Examples
        ----------
        >>> from ler.len_galaxy_population import OpticalDepth
        >>> ler = OpticalDepth()  # with default LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        >>> comoving_distance = ler.comoving_distance(1.)
        >>> comoving_distance = ler.comoving_distance.function(np.array([1., 2.]))
        >>> redshift = ler.comoving_distance.function_inverse(np.array([100., 200.]))
        """

        # for angular diameter distance
        resolution = self.create_new_interpolator["angular_diameter_distance"]["resolution"]
        create_new = self.create_new_interpolator["angular_diameter_distance"]["create_new"]
        zs = np.geomspace(z_min, z_max, resolution)
        Da = self.cosmo.angular_diameter_distance(zs).value
        self.angular_diameter_distance = FunctionConditioning(
            function=Da,
            x_array=zs,
            conditioned_y_array=None,
            param_dict_given=dict(z_min=z_min, z_max=z_max, cosmology=self.cosmo, resolution=resolution, details="angular_diameter_distance from astropy.cosmology"),
            directory=self.directory,
            sub_directory="angular_diameter_distance",
            name="angular_diameter_distance",
            create_new=create_new,
            create_function_inverse=False,
            create_function=True,
            create_pdf=False,
            create_rvs=False,
            callback='function',
        )
        self.angular_diameter_distance.__doc__ = """
        Redshift to angular diameter distance conversion.

        Parameters
        ----------
        zs : `numpy.ndarray` or `float`
            Source redshifts

        Returns
        ----------
        angular_diameter_distance : `numpy.ndarray`
            angular diameter distance in Mpc

        Examples
        ----------
        >>> from ler.len_galaxy_population import OpticalDepth
        >>> ler = OpticalDepth()  # with default LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        >>> angular_diameter_distance = ler.angular_diameter_distance(1.)
        >>> angular_diameter_distance = ler.angular_diameter_distance.function(np.array([1., 2.]))
        >>> redshift = ler.angular_diameter_distance.function_inverse(np.array([100., 200.]))
        """

        # for angular diameter distance between two redshifts
        angular_diameter_distance_z1z2 = lambda zl0, zs0: (self.angular_diameter_distance.function(zs0)*(1.+zs0) - self.angular_diameter_distance.function(zl0)*(1.+zl0))/(1.+zs0)
        self.angular_diameter_distance_z1z2 = FunctionConditioning(
            function=None,
            x_array=None,
            param_dict_given=dict(z_min=z_min, z_max=z_max, cosmology=self.cosmo, resolution=resolution, details="angular_diameter_distance_z1z2 from astropy.cosmology"),
            create_function=angular_diameter_distance_z1z2,
            callback='function',
        )
        self.angular_diameter_distance_z1z2.__doc__ = """
        Redshift to angular diameter distance conversion.

        Parameters
        ----------
        zl0 : `numpy.ndarray` or `float`
            Lens redshifts
        zs0 : `numpy.ndarray` or `float`
            Source redshifts

        Returns
        ----------
        angular_diameter_distance_z1z2 : `numpy.ndarray`
            angular diameter distance in Mpc

        Examples
        ----------
        >>> from ler.len_galaxy_population import OpticalDepth
        >>> ler = OpticalDepth()  # with default LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        >>> angular_diameter_distance_z1z2 = ler.angular_diameter_distance_z1z2(1., 2.)
        >>> angular_diameter_distance_z1z2 = ler.angular_diameter_distance_z1z2.function(np.array([1., 2.]), np.array([1., 2.]))
        """

        # get differential co-moving volume interpolator
        resolution = self.create_new_interpolator["differential_comoving_volume"]["resolution"]
        create_new = self.create_new_interpolator["differential_comoving_volume"]["create_new"]
        zs = np.geomspace(z_min, z_max, resolution)
        dVcdz = self.cosmo.differential_comoving_volume(zs).value * 4 * np.pi  # volume of shell in Mpc^3
        self.differential_comoving_volume = FunctionConditioning(
            function=dVcdz,
            x_array=zs,
            conditioned_y_array=None,
            param_dict_given=dict(z_min=z_min, z_max=z_max, cosmology=self.cosmo, resolution=resolution, details="differential_comoving_volume from astropy.cosmology"),
            directory=self.directory,
            sub_directory="differential_comoving_volume",
            name="differential_comoving_volume",
            create_new=create_new,
            create_function_inverse=False,
            create_function=True,
            create_pdf=False,
            create_rvs=False,
            callback='function',
        )
        self.differential_comoving_volume.__doc__ = """
        Redshift to differential comoving volume conversion.

        Parameters
        ----------
        zs : `numpy.ndarray` or `float`
            Source redshifts

        Returns
        ----------
        differential_comoving_volume : `numpy.ndarray`
            differential comoving volume in Mpc^3

        Examples
        ----------
        >>> from ler.len_galaxy_population import OpticalDepth
        >>> ler = OpticalDepth()  # with default LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        >>> differential_comoving_volume = ler.differential_comoving_volume(1.)
        >>> differential_comoving_volume = ler.differential_comoving_volume.function(np.array([1., 2.]))
        """

    ########################
    # Attribute properties #
    ########################
    # @property
    # def cross_section(self):
    #     return self._cross_section
    
    # @cross_section.setter
    # def cross_section(self, input_function):
    #     if input_function in self.available_cross_section_list_and_its_params:
    #         print(f"using ler available cross section function : {input_function}")
    #         self._cross_section = getattr(self, input_function)(zs=None, get_attribute=True)
    #     elif callable(input_function):
    #         print(f"using user provided custom cross section function")
    #         self._cross_section = input_function
    #     elif isinstance(input_function, object):
    #         print(f"using user provided custom cross section class/object")
    #         self._cross_section = input_function
    #     else:
    #         raise ValueError("input_function not in available_cross_section_list_and_its_params")


    @property
    def optical_depth(self):
        """
        Function to compute the strong lensing optical depth.

        Parameters
        ----------
        zs : `numpy.ndarray` (1D array of floats)
            source redshifts

        Returns
        -------
        tau : `numpy.ndarray` (1D array of floats)
            strong lensing optical depth

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> print(od.optical_depth(np.array([0.1,0.2,0.3])))
        """

        return self._optical_depth

    @optical_depth.setter
    def optical_depth(self, input_function):
        if input_function in self.available_lens_functions_and_its_params['optical_depth']:
            print(f"using ler available optical depth function : {input_function}")
            self._optical_depth = getattr(self, input_function)(zs=None, get_attribute=True)
        elif callable(input_function):
            print(f"using user provided custom optical depth function")
            self._optical_depth = FunctionConditioning(function=None, x_array=None, create_function=input_function)
        elif isinstance(input_function, object):
            print(f"using user provided custom optical depth class/object")
            self._optical_depth = input_function
        else:
            raise ValueError("input_function not in available_lens_functions_and_its_params['optical_depth']")
        

    @property
    def velocity_dispersion(self):
        """
        Class object to sample velocity dispersion. `zl` is required only if velocity dispersion sampler is redshift dependent.

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        zl : `float`
            redshift of the lens galaxy

        Returns
        -------
        sigma : `numpy.ndarray` (1D array of floats)
            velocity dispersion of the lens galaxy

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> print(od.velocity_dispersion(size=10))
        """

        return self._velocity_dispersion

    @velocity_dispersion.setter
    def velocity_dispersion(self, prior):
        if prior=="velocity_dispersion_choi":
            prior = "velocity_dispersion_bernardi"
    
        if prior in self.available_lens_prior_list_and_its_params["velocity_dispersion"]:
            print(f"using ler available velocity dispersion function : {prior}")
            self._velocity_dispersion = getattr(self, prior)(size=None, zl=None, get_attribute=True)
        elif callable(prior):
            print("using user provided custom velocity_dispersion function")
            self._velocity_dispersion = FunctionConditioning(function=None, x_array=None,create_rvs=prior)
        elif isinstance(prior, object):
            print("using user provided custom velocity_dispersion class/object")
            self._velocity_dispersion = prior
        else:
            raise ValueError(f"velocity_dispersion prior not in available_lens_prior_list_and_its_params['velocity_dispersion']")

    @property
    def axis_ratio(self):
        """
        Function to sample axis ratio from rayleigh distribution with given velocity dispersion.

        Parameters
        ----------
        sigma : `numpy.ndarray` (1D array of floats)
            velocity dispersion of the lens galaxy

        Returns
        -------
        q : `numpy.ndarray` (1D array of floats)
            axis ratio of the lens galaxy

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> print(od.axis_ratio(sigma=200.))
        """

        return self._axis_ratio

    @axis_ratio.setter
    def axis_ratio(self, prior):
        if prior in self.available_lens_prior_list_and_its_params["axis_ratio"]:
            print(f"using ler available axis_ratio function : {prior}")
            self._axis_ratio = getattr(self, prior)(size=None,sigma=None, get_attribute=True)
        elif callable(prior):
            print("using user provided custom axis_ratio function")
            self._axis_ratio = FunctionConditioning(function=None, x_array=None,create_rvs=prior)
        elif isinstance(prior, object):
            print("using user provided custom axis_ratio class/object")
            self._axis_ratio = prior
        else:
            raise ValueError("prior not in available_lens_prior_list_and_its_params['axis_ratio']")
        
    @property
    def lens_redshift(self):
        return self._lens_redshift
    
    @lens_redshift.setter
    def lens_redshift(self, prior):
        if prior in self.available_lens_prior_list_and_its_params["lens_redshift"]:
            print(f"using ler available lens_redshift function : {prior}")
            self._lens_redshift = getattr(self, prior)(size=None, zs=None, get_attribute=True)
        elif callable(prior):
            print("using user provided custom lens_redshift function")
            self._lens_redshift = FunctionConditioning(function=prior, original_function=True)
        elif isinstance(prior, object):
            print("using user provided custom lens_redshift class/object")
            self._lens_redshift = prior
        else:
            raise ValueError("prior not in available_lens_prior_list_and_its_params['lens_redshift']")
    
    @property
    def axis_rotation_angle(self):
        return self._axis_rotation_angle
    
    @axis_rotation_angle.setter
    def axis_rotation_angle(self, prior):
        if prior in self.available_lens_prior_list_and_its_params["axis_rotation_angle"]:
            print(f"using ler available axis_rotation_angle function : {prior}")
            self._axis_rotation_angle = getattr(self, prior)(size=None, get_attribute=True)
        elif callable(prior):
            print("using user provided custom axis_rotation_angle function")
            self._axis_rotation_angle = FunctionConditioning(function=prior, original_function=True)
        elif isinstance(prior, object):
            print("using user provided custom axis_rotation_angle class/object")
            self._axis_rotation_angle = prior
        else:
            raise ValueError("prior not in available_lens_prior_list_and_its_params['axis_rotation_angle']")
        
    @property
    def external_shear(self):
        return self._external_shear
    
    @external_shear.setter
    def external_shear(self, prior):
        if prior in self.available_lens_prior_list_and_its_params["external_shear"]:
            print(f"using ler available external_shear function : {prior}")
            self._external_shear = getattr(self, prior)(size=None, get_attribute=True)
        elif callable(prior):
            print("using user provided custom external_shear function")
            self._external_shear = FunctionConditioning(function=prior, original_function=True)
        elif isinstance(prior, object):
            print("using user provided custom external_shear class/object")
            self._external_shear = prior
        else:
            raise ValueError("prior not in available_lens_prior_list_and_its_params['external_shear']")
        
    @property
    def density_profile_slope(self):
        return self._density_profile_slope
    
    @density_profile_slope.setter
    def density_profile_slope(self, prior):
        if prior in self.available_lens_prior_list_and_its_params["density_profile_slope"]:
            print(f"using ler available density_profile_slope function : {prior}")
            self._density_profile_slope = getattr(self, prior)(size=None, get_attribute=True)
        elif callable(prior):
            print("using user provided custom density_profile_slope function")
            self._density_profile_slope = FunctionConditioning(function=prior, original_function=True)
        elif isinstance(prior, object):
            print("using user provided custom density_profile_slope class/object")
            self._density_profile_slope = prior
        else:
            raise ValueError("prior not in available_lens_prior_list_and_its_params['density_profile_slope']")
        

    ####################################
    # Available samplers and functions #
    ####################################
    @property
    def available_lens_prior_list_and_its_params(self):
        """
        Dictionary with list all the available priors and it's corresponding parameters. This is an immutable instance attribute.
        """
        
        self._available_lens_prior_list_and_its_params = dict(
            source_redshift_sl = dict(
                strongly_lensed_source_redshifts=None
            ),
            lens_redshift=dict(
                lens_redshift_SDSS_catalogue_sis=None,
                lens_redshift_SDSS_catalogue_numerical=None,
            ),
            velocity_dispersion = dict(
                velocity_dispersion_gengamma = dict(
                    sigma_min=10., sigma_max=420., alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78
                ),
                velocity_dispersion_choi = dict(
                    sigma_min = 50., sigma_max = 420., alpha = 2.32, beta = 2.67, phistar = 8.0e-3*self.cosmo.h**3, sigmastar = 161.0
                ),
                velocity_dispersion_bernardi = dict(
                    sigma_min=10., sigma_max=420., alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78
                ),
                velocity_dispersion_ewoud = dict(
                    sigma_min=10., sigma_max=420., alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78
                ),
            ),
            axis_ratio = dict(
                axis_ratio_rayleigh = dict(q_min=0.2, q_max=1.0),
                axis_ratio_padilla_strauss = dict(q_min=0.2, q_max=1.0),
                axis_ratio_SIS = None,
            ),
            axis_rotation_angle = dict(
                axis_rotation_angle_uniform = dict(
                    phi_min=0.0, phi_max=2 * np.pi
                ),
            ),
            external_shear = dict(
                external_shear_normal = dict(mean=0., std=0.05),
            ),
            density_profile_slope = dict(
                density_profile_slope_normal=dict(mean=1.99, std=0.149),
            ),
            source_parameters=dict(
                sample_gw_parameters=None
            ),
        )

        return self._available_lens_prior_list_and_its_params
    
    @property
    def available_lens_functions_and_its_params(self):
        """
        Dictionary with list all the available lens functions. This is an immutable instance attribute.
        """

        self._available_lens_functions_and_its_params = dict(
            strong_lensing_condition=dict(
                rjs_with_cross_section_sie_feixu=None,
                rjs_with_cross_section_sis=None,
                rjs_with_cross_section_epl_shear=None,
            ),
            optical_depth=dict(
                optical_depth_sis_haris=None,
                optical_depth_sie_hemanta=None,
                optical_depth_epl_shear_hemanta=None,
                optical_depth_numerical=None,
            ),
            param_sampler_type=dict(
                sample_all_routine_sie=None,
                sample_all_routine_sis=None,
                sample_all_routine_epl_shear=None,
            ),
            cross_section=dict(
                cross_section_sie_feixu=None,
                cross_section_sis=None,
                cross_section_caustic_area=None,
                interpolated_cross_section_function=None,
            ),
        )

        return self._available_lens_functions_and_its_params









