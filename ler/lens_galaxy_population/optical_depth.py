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


from ..utils import cubic_spline_interpolator, inverse_transform_sampler, cubic_spline_interpolator2d_array, save_pickle, interpolator_pickle_path, FunctionConditioning, inverse_transform_sampler2d, pdf_cubic_spline_interpolator2d_array, normal_pdf, normal_pdf_2d, load_txt_from_module

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

        # default lens functions and samplers
        self.lens_param_samplers, self.lens_param_samplers_params, self.lens_functions, self.lens_functions_params = self.default_lens_samplers_and_functions(lens_type)

        # change default lens functions and samplers if user input is given
        # if lens functions and samplers are given but not parameters, then use default parameters
        # it updates self.lens_param_samplers, self.lens_param_samplers_params, self.lens_functions, self.lens_functions_params
        self.lens_functions_and_sampler_categorization(lens_param_samplers, lens_param_samplers_params, lens_functions, lens_functions_params);

        # setting up astropy cosmology functions
        self.create_lookup_table_fuction();

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
        self.density_profile_slope_sl = self.lens_param_samplers['density_profile_slope_sl']
        self.external_shear_sl = self.lens_param_samplers['external_shear_sl']
        
        # lens function initialization
        self.optical_depth = self.lens_functions['optical_depth']

    ##################
    # Initialization #
    ##################
    def default_lens_samplers_and_functions(self, lens_type):
        """
            Function to categorize the lens priors/samplers

            Parameters
            ----------
                lens_type : `str`
                    lens type
                    e.g. 'epl_shear_galaxy' for elliptical power-law galaxy

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
                lens_redshift="lens_redshift_SDSS_catalogue_hemanta", #"lens_redshift_SDSS_catalogue_numerical",
                velocity_dispersion="velocity_dispersion_ewoud",
                axis_ratio="axis_ratio_rayleigh",
                axis_rotation_angle="axis_rotation_angle_uniform",
                external_shear="external_shear_normal",
                density_profile_slope="density_profile_slope_normal",
                external_shear_sl="external_shear_normal",
                density_profile_slope_sl="density_profile_slope_normal",
            )
            lens_priors_params_ = dict(
                source_redshift_sl=None,
                lens_redshift=None,
                velocity_dispersion=dict(sigma_min=100., sigma_max=400., alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78),
                axis_ratio=dict(q_min=0.2, q_max=1.),
                axis_rotation_angle=dict(phi_min=0.0, phi_max=2 * np.pi),
                external_shear=dict(mean=0., std=0.05),
                density_profile_slope=dict(mean=1.99, std=0.149),
                external_shear_sl=dict(mean=0., std=0.05, name="external_shear_sl"),
                density_profile_slope_sl=dict(mean=2.078, std=0.16, name="density_profile_slope_sl"),
            )
            lens_functions_ = dict(
                param_sampler_type="sample_all_routine_sie_sl",
                strong_lensing_condition="rjs_with_cross_section_sie_feixu",
                optical_depth="optical_depth_numerical",
                cross_section="interpolated_cross_section_function",
            )
            lens_functions_params_ = dict(
                strong_lensing_condition=None,
                optical_depth=None,
                param_sampler_type=dict(fast_sampler=True), # not implemented yet
                cross_section=None,
            )
        elif lens_type == "sie_galaxy":
            lens_priors_ = dict(
                source_redshift_sl="strongly_lensed_source_redshifts",
                lens_redshift="lens_redshift_SDSS_catalogue_hemanta", #"lens_redshift_SDSS_catalogue_numerical",
                velocity_dispersion="velocity_dispersion_ewoud",
                axis_ratio="axis_ratio_rayleigh",
                axis_rotation_angle="axis_rotation_angle_uniform",
                external_shear="external_shear_normal",
                density_profile_slope="density_profile_slope_normal",
                external_shear_sl="external_shear_normal",
                density_profile_slope_sl="density_profile_slope_normal",
            )
            lens_priors_params_ = dict(
                source_redshift_sl=None,
                lens_redshift=None,
                velocity_dispersion=dict(sigma_min=100., sigma_max=400., alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78),
                axis_ratio=dict(q_min=0.2, q_max=1.),
                axis_rotation_angle=dict(phi_min=0.0, phi_max=2 * np.pi),
                external_shear=dict(mean=0., std=0.),
                density_profile_slope=dict(mean=2., std=0.),
                external_shear_sl=dict(mean=0., std=0.),
                density_profile_slope_sl=dict(mean=2., std=0.),
            )
            lens_functions_ = dict(
                strong_lensing_condition="rjs_with_cross_section_sie_feixu",
                optical_depth="optical_depth_numerical",
                param_sampler_type="sample_all_routine_sie_sl",
                cross_section="cross_section_sie_feixu",
            )
            lens_functions_params_ = dict(
                strong_lensing_condition=None,
                optical_depth=dict(interpolated_cross_section=True),
                param_sampler_type=None,
                cross_section=None,
            )
        elif lens_type == "sis_galaxy":
            lens_priors_ = dict(
                source_redshift_sl="strongly_lensed_source_redshifts",
                lens_redshift="lens_redshift_SDSS_catalogue_numerical",
                velocity_dispersion="velocity_dispersion_choi",
                axis_ratio="axis_ratio_uniform",
                axis_rotation_angle="axis_rotation_angle_uniform",
                external_shear="external_shear_normal",
                density_profile_slope="density_profile_slope_normal",
                external_shear_sl="external_shear_normal",
                density_profile_slope_sl="density_profile_slope_normal",
            )
            lens_priors_params_ = dict(
                source_redshift_sl=None,
                lens_redshift=None,
                velocity_dispersion=dict(sigma_min = 10., sigma_max = 420., alpha = 2.32, beta = 2.67, phistar = 8.0e-3*self.cosmo.h**3, sigmastar = 161.0),
                axis_ratio=dict(q_min=1., q_max=1.),
                axis_rotation_angle=dict(phi_min=0.0, phi_max=0.),
                external_shear=dict(mean=0., std=0.),
                density_profile_slope=dict(mean=2., std=0.),
                external_shear_sl=dict(mean=0., std=0.),
                density_profile_slope_sl=dict(mean=2., std=0.),
            )
            lens_functions_ = dict(
                strong_lensing_condition="rjs_with_cross_section_sis",
                optical_depth="optical_depth_numerical",
                param_sampler_type="sample_all_routine_sis",
                cross_section="cross_section_sis",
            )
            lens_functions_params_ = dict(
                strong_lensing_condition=None,
                optical_depth=dict(interpolated_cross_section=True),
                param_sampler_type=None,
                cross_section=None,
            )
        else:
            raise ValueError("lens_type not recognized")

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
        
    def lens_functions_and_sampler_categorization(self, lens_param_samplers, lens_param_samplers_params, lens_functions, lens_functions_params):
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

        # update the priors if input is given
        if lens_param_samplers:
            self.lens_param_samplers.update(lens_param_samplers)
        if lens_param_samplers_params:
            self.lens_param_samplers_params.update(lens_param_samplers_params)
        if lens_functions:
            self.lens_functions.update(lens_functions)
        if lens_functions_params:
            self.lens_functions_params.update(lens_functions_params)

        sampler_prior_names = ['lens_redshift', 'velocity_dispersion', 'axis_ratio', 'axis_rotation_angle', 'external_shear', 'density_profile_slope', 'external_shear_sl', 'density_profile_slope_sl']

        # if there is user input sampler prior params
        # you can't only give lens_param_samplers_params, you have to give lens_param_samplers as well if you want to update the sampler priors
        for name in sampler_prior_names:  # e.g. name='axis_ratio'
            # avoid None
            if (lens_param_samplers is not None) and (name in lens_param_samplers):
                # if there is user input function name, update the sampler priors
                sampler_name = lens_param_samplers[name]  # e.g. 'axis_ratio_rayleigh'
                # check input sampler is string or function
                if isinstance(sampler_name, str):
                    # available lens samplers for name e.g. 'axis_ratio'
                    dict_ = self.available_lens_prior_list_and_its_params[name]  # e.g. {'axis_ratio_padilla_strauss': {'q_min': 0.2, 'q_max': 1.0}, ....}
                    if sampler_name in dict_:  # e.g. 'axis_ratio_padilla_strauss'
                        param_dict = dict_[sampler_name]  # e.g. {'q_min': 0.2, 'q_max': 1.0}
                        if (lens_param_samplers_params is None) or (lens_param_samplers_params[name] is None):  # not a dictionary
                            self.lens_param_samplers_params[name] = param_dict
                        else:  # if there is user input lens_param_samplers_params
                            param_dict.update(lens_param_samplers_params[name])
                            self.lens_param_samplers_params[name].update(param_dict)  # user inputs will override the default values
                    else:
                        raise ValueError(f"{name} sampler {sampler_name} not available.\n Available {name} samplers and its parameters are: {dict_[name]}")
                elif not callable(lens_param_samplers[name]):
                    raise ValueError(f"Given {name} sampler should be either a string name of available sampler or a function")

                
        lens_function_names = ['optical_depth', 'cross_section']

        # if there is user input function, update the sampler priors
        for name in lens_function_names:
            if (lens_functions is not None) and (name in lens_functions):
                function_name = lens_functions[name] 
                if isinstance(function_name, str):
                     # available lens functions for name e.g. 'optical_depth'
                    dict_ = self.available_lens_functions_and_its_params[name]  
                    if function_name in dict_:
                        param_dict = dict_[function_name]  
                        if (lens_functions_params is None) or (lens_functions_params[name] is None):  # not a dictionary
                            self.lens_functions_params[name] = param_dict
                        else:  # if there is user input lens_functions_params
                            param_dict.update(lens_functions_params[name])
                            self.lens_functions_params[name].update(param_dict)   # user inputs will override the default values
                    else:
                        raise ValueError(f"{name} function {function_name} not available.\n Available {name} functions and its parameters are: {dict_[name]}")
                elif not callable(lens_functions[name]):
                    raise ValueError(f"Given {name} function should be either a string name of available function or a function")
                
        # setting up multiprocessor function for lens redshift
        if self.lens_type == 'epl_shear_galaxy':
            if self.lens_param_samplers['velocity_dispersion'] == 'velocity_dispersion_ewoud':
                print('using lens_redshift_epl_shear2_mp')
                from .mp import lens_redshift_epl_shear2_mp
                self._helper_sl_disribution_mp = lens_redshift_epl_shear2_mp
            else:
                print('using lens_redshift_epl_shear1_mp')
                from .mp import lens_redshift_epl_shear1_mp
                self._helper_sl_disribution_mp = lens_redshift_epl_shear1_mp
        elif self.lens_type == 'sie_galaxy':
            if self.lens_param_samplers['velocity_dispersion'] == 'velocity_dispersion_ewoud':
                if self.lens_functions['cross_section']=='interpolated_cross_section_function':
                    print('using lens_redshift_sie4_mp')
                    from .mp import lens_redshift_sie4_mp
                    self._helper_sl_disribution_mp = lens_redshift_sie4_mp
                else:
                    print('using lens_redshift_sie2_mp')
                    from .mp import lens_redshift_sie2_mp
                    self._helper_sl_disribution_mp = lens_redshift_sie2_mp
            else:
                if self.lens_functions['cross_section']=='interpolated_cross_section_function':
                    print('using lens_redshift_sie3_mp')
                    from .mp import lens_redshift_sie3_mp
                    self._helper_sl_disribution_mp = lens_redshift_sie3_mp
                else:
                    print('using lens_redshift_sie1_mp')
                    from .mp import lens_redshift_sie1_mp
                    self._helper_sl_disribution_mp = lens_redshift_sie1_mp
        elif self.lens_type == 'sis_galaxy':
            if self.lens_param_samplers['velocity_dispersion'] == 'velocity_dispersion_ewoud':
                print('using lens_redshift_sis2_mp')
                from .mp import lens_redshift_sis2_mp
                self._helper_sl_disribution_mp = lens_redshift_sis2_mp
            else:
                print('using lens_redshift_sis1_mp')
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

        identifier_dict = {'name': "axis_ratio_rayleigh"}
        identifier_dict['velocity_dispersion'] = self.lens_param_samplers_params['velocity_dispersion']
        if identifier_dict['velocity_dispersion']: # if velocity_dispersion is not None
            identifier_dict['velocity_dispersion']['name'] = str(self.lens_param_samplers['velocity_dispersion'])
        identifier_dict['resolution'] = self.create_new_interpolator["axis_ratio"]["resolution"]
        param_dict = self.available_lens_prior_list_and_its_params["axis_ratio"]["axis_ratio_rayleigh"]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        q_array = np.linspace(
            identifier_dict["q_min"], 
            identifier_dict["q_max"],
            identifier_dict["resolution"],
        )
        sigma_array = np.linspace(
            identifier_dict['velocity_dispersion']["sigma_min"],
            identifier_dict['velocity_dispersion']["sigma_max"],
            500,
        )

        q_pdf = lambda q, sigma: axis_ratio_rayleigh_pdf(
            q=q,
            sigma=sigma,
            q_min=identifier_dict["q_min"],
            q_max=identifier_dict["q_max"],
        )

        q_object = FunctionConditioning(
            function=q_pdf,
            x_array=q_array,
            conditioned_y_array=sigma_array,
            param_dict_given=identifier_dict,
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

        return q_object if get_attribute else q_object(size, sigma)
        
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

        identifier_dict = {'name': "axis_ratio_padilla_strauss"}
        identifier_dict['resolution'] = self.create_new_interpolator["axis_ratio"]["resolution"]
        param_dict = self.available_lens_prior_list_and_its_params["axis_ratio"]["axis_ratio_padilla_strauss"]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        # Using Padilla and Strauss 2008 distribution for axis ratio
        q_array = np.array([0.04903276402927845, 0.09210526315789469, 0.13596491228070173, 0.20789473684210524, 0.2899703729522482, 0.3230132450331126, 0.35350877192982455, 0.37946148483792264, 0.4219298245614036, 0.4689525967235971, 0.5075026141512723, 0.5226472638550018, 0.5640350877192983, 0.6096491228070177, 0.6500000000000001, 0.6864848379226213, 0.7377192982456142, 0.7787295224817011, 0.8007581038689441, 0.822786685256187, 0.8668438480306729, 0.8973684210526317, 0.9254385964912283])
        pdf = np.array([0.04185262687135349, 0.06114520695141845, 0.096997499638376, 0.1932510900336828, 0.39547914337673706, 0.49569751276216234, 0.6154609137685201, 0.7182049959882812, 0.920153741243567, 1.1573982157399754, 1.3353263628106684, 1.413149656448315, 1.5790713532948977, 1.7280185150744938, 1.8132994441344819, 1.8365803753840484, 1.8178662203211204, 1.748929843583365, 1.688182592496342, 1.6274353414093188, 1.4948487090314488, 1.402785526832393, 1.321844068356993])

        q_object = FunctionConditioning(
            function=pdf,
            x_array=q_array,
            conditioned_y_array=None,
            param_dict_given=identifier_dict,
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

        return q_object if get_attribute else q_object(size)
        
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
    
        identifier_dict = {}
        identifier_dict['name'] = "lens_redshift_numerical_"+self.lens_type
        identifier_dict['resolution'] = self.create_new_interpolator["lens_redshift"]["resolution"]
        identifier_dict['z_min'] = self.z_min
        identifier_dict['z_max'] = self.z_max
        identifier_dict['cosmology'] = self.cosmo
        identifier_dict['velocity_dispersion'] = self.lens_param_samplers_params['velocity_dispersion']
        if identifier_dict['velocity_dispersion']: # if velocity_dispersion is not None
            identifier_dict['velocity_dispersion']['name'] = str(self.lens_param_samplers['velocity_dispersion'])
        identifier_dict['axis_ratio'] = self.lens_param_samplers_params['axis_ratio']
        if identifier_dict['axis_ratio']: # if axis_ratio is not None
            identifier_dict['axis_ratio']['name'] = str(self.lens_param_samplers['axis_ratio'])
        identifier_dict['axis_rotation_angle'] = self.lens_param_samplers_params['axis_rotation_angle']
        if identifier_dict['axis_rotation_angle']: # if axis_rotation_angle is not None
            identifier_dict['axis_rotation_angle']['name'] = str(self.lens_param_samplers['axis_rotation_angle'])
        identifier_dict['density_profile_slope'] = self.lens_param_samplers_params['density_profile_slope']
        if identifier_dict['density_profile_slope']: # if density_profile_slope is not None
            identifier_dict['density_profile_slope']['name'] = str(self.lens_param_samplers['density_profile_slope'])
        identifier_dict['external_shear'] = self.lens_param_samplers_params['external_shear']
        if identifier_dict['external_shear']: # if external_shear is not None
            identifier_dict['external_shear']['name'] = str(self.lens_param_samplers['external_shear'])
        
        param_dict = self.available_lens_prior_list_and_its_params["lens_redshift"]["lens_redshift_SDSS_catalogue_numerical"]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        print("Numerically solving the lens redshift distribution...")
        zs_resolution = self.create_new_interpolator["optical_depth"]["resolution"]
        zs_array = np.geomspace(self.z_min+0.001, self.z_max, zs_resolution) if self.z_min==0 else np.geomspace(self.z_min, self.z_max, zs_resolution)

        zl_scaled2d = []
        for i, zs_ in enumerate(zs_array):
            buffer_ = np.linspace(0.0001, zs_-0.0001, identifier_dict['resolution'])
            zl_scaled2d.append(buffer_/zs_)
        zl_scaled2d = np.array(zl_scaled2d)

        _, it_exist = interpolator_pickle_path(
            param_dict_given=identifier_dict,
            directory=self.directory,
            sub_directory="lens_redshift",
            interpolator_name=identifier_dict['name'],
        )

        create_new = self.create_new_interpolator["lens_redshift"]["create_new"]
        if not it_exist or create_new:
            number_density = self.lens_redshift_multiprocessing(zl_scaled2d, zs_array)
            # number density is zero for zl=0 and infinite for zl=zs
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
            param_dict_given=identifier_dict,
            directory=self.directory,
            sub_directory="lens_redshift",
            name=identifier_dict['name'],
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

        return zl_object if get_attribute else zl_object(size, zs)
    
    def lens_redshift_SDSS_catalogue_hemanta(self, size=1000, zs=None, get_attribute=False, **kwargs):
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
    
        # default parameters
        identifier_dict = {'name': 'lens_redshift_numerical_epl_shear_galaxy',
            'resolution': 50,
            'z_min': 0.0,
            'z_max': 10.0,
            'cosmology': LambdaCDM(H0=70, Om0=0.3, Ode0=0.7),
            'velocity_dispersion': {'sigma_min': 100.0,
            'sigma_max': 400.0,
            'alpha': 0.94,
            'beta': 1.85,
            'phistar': 0.02099,
            'sigmastar': 113.78,
            'name': 'velocity_dispersion_ewoud'},
            'axis_ratio': {'q_min': 0.2, 'q_max': 1.0, 'name': 'axis_ratio_rayleigh'},
            'axis_rotation_angle': {'phi_min': 0.0,
            'phi_max': 6.283185307179586,
            'name': 'axis_rotation_angle_uniform'},
            'density_profile_slope': {'mean': 1.99,
            'std': 0.149,
            'name': 'density_profile_slope_normal'},
            'external_shear': {'mean': 0.0, 'std': 0.05, 'name': 'external_shear_normal'}
        }
        
        print("Getting pre-computed lens redshift distribution...")
        zs_resolution = 48
        zs_array = np.geomspace(0.001, 10., zs_resolution)
        zl_scaled2d = []
        for i, zs_ in enumerate(zs_array):
            buffer_ = np.linspace(0.0001, zs_-0.0001, identifier_dict['resolution'])
            zl_scaled2d.append(buffer_/zs_)
        zl_scaled2d = np.array(zl_scaled2d)

        # get number_density
        number_density = load_txt_from_module('ler', 'lens_galaxy_population.lens_param_data', 'number_density_zl_zs.txt')

        # number density is zero for zl=0 and infinite for zl=zs
        # Adding zero at the first element of each row
        zl_scaled2d = np.hstack((np.zeros((zl_scaled2d.shape[0], 1)), zl_scaled2d))
        # number_density = np.hstack((np.zeros((number_density.shape[0], 1)), number_density))
        # Adding one at the last element of each row of zl_scaled2d
        zl_scaled2d = np.hstack((zl_scaled2d, np.ones((zl_scaled2d.shape[0], 1))))
        # Adding zero at the last element of each row of density
       # number_density = np.hstack((number_density, np.zeros((number_density.shape[0], 1))))

        zl_object = FunctionConditioning(
            function=number_density,
            x_array=zl_scaled2d,
            conditioned_y_array=zs_array,
            param_dict_given=identifier_dict,
            directory=self.directory,
            sub_directory="lens_redshift",
            name=identifier_dict['name'],
            create_new=False,
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

        return zl_object if get_attribute else zl_object(size, zs)
        
    def intrinsic_lens_redshift(self, size=1000, get_attribute=False, **kwargs):
        """
        Function to sample intrinsic lens redshifts, based on the intrinsic velocity dispersion of the lens galaxy.

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
    
        identifier_dict = {'name': "intrinsic_lens_redshift"}
        identifier_dict['z_min'] = self.z_min
        identifier_dict['z_max'] = self.z_max
        identifier_dict['cosmology'] = self.cosmo
        identifier_dict['velocity_dispersion'] = self.lens_param_samplers_params['velocity_dispersion']
        if identifier_dict['velocity_dispersion']: # if velocity_dispersion is not None
            identifier_dict['velocity_dispersion']['name'] = str(self.lens_param_samplers['velocity_dispersion'])

        identifier_dict['resolution'] = 500
        # param_dict = self.available_lens_prior_list_and_its_params["lens_redshift"]["intrinsic_lens_redshift"]
        # if param_dict:
        #     param_dict.update(kwargs)
        # else:
        #     param_dict = kwargs
        # identifier_dict.update(param_dict)

        _, it_exist = interpolator_pickle_path(
            param_dict_given=identifier_dict,
            directory=self.directory,
            sub_directory="lens_redshift",
            interpolator_name=identifier_dict['name'],
        )

        create_new = self.create_new_interpolator["lens_redshift"]["create_new"]
        if not it_exist or create_new:
            z_min = 0.001 if self.z_min==0 else self.z_min
            z_max = self.z_max
            zl_array = np.linspace(z_min, z_max, identifier_dict['resolution']) 
            integrand = lambda sigma, z: self.velocity_dispersion.function(np.array([sigma]), np.array([z]))[0]*self.differential_comoving_volume.function(np.array([z]))[0]
            # integral = [quad(integrand, 10., 420., args=(z))[0] for z in zl_array]
            integral = [quad(integrand, identifier_dict['velocity_dispersion']['sigma_min'], identifier_dict['velocity_dispersion']['sigma_max'], args=(z))[0] for z in zl_array]
        else:
            zl_array = None
            integral=None

        zl_object = FunctionConditioning(
            function=integral,
            x_array=zl_array,
            param_dict_given=identifier_dict,
            directory=self.directory,
            sub_directory="lens_redshift",
            name=identifier_dict['name'],
            create_new=create_new,
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='rvs',
        )

        return zl_object if get_attribute else zl_object(size)
        
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
        try:
            q_args = [
                np.array(self.axis_ratio.cdf_values),
                np.array(self.axis_ratio.x_array),
                np.array(self.axis_ratio.conditioned_y_array),
            ]
        except:
            # SIS case
            q_args = [
                self.axis_ratio.info["q_min"],
                self.axis_ratio.info["q_max"],
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
        """
        Function to sample the axis rotation angle of the elliptical lens galaxy from a uniform distribution.

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        get_attribute : `bool`
            if True, returns a function that can be called with size as input

        Returns
        -------
        phi : `numpy.ndarray`
            axis rotation angle of the elliptical lens galaxy

        Examples
        --------

        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_param_samplers=dict(axis_rotation_angle="axis_rotation_angle_uniform"))
        >>> print(od.axis_rotation_angle_uniform(size=10))
        """

        identifier_dict = {'name': "axis_rotation_angle_uniform"}
        param_dict = self.available_lens_prior_list_and_its_params["axis_rotation_angle"]["axis_rotation_angle_uniform"]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        low = param_dict["phi_min"]
        high = param_dict["phi_max"]
        phi_rvs = njit(lambda size: np.random.uniform(
            low=low,
            high=high,
            size=size,
        ))
        if param_dict["phi_max"] == param_dict["phi_min"]:
            phi_pdf = lambda phi: 1. if phi == param_dict["phi_min"] else 0.
        else:
            phi_pdf = lambda phi: 1/(param_dict["phi_max"]-param_dict["phi_min"])

        phi_object =FunctionConditioning(
            function=None,
            x_array=None,
            param_dict_given=param_dict,
            create_rvs=phi_rvs,
            create_pdf=phi_pdf,
            callback='rvs',
        )

        return phi_object if get_attribute else phi_object.rvs(size)
    
    def axis_ratio_uniform(self, size, get_attribute=False, **kwargs):
        """
        Function to sample the axis ratio of the elliptical lens galaxy from a uniform distribution.

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        get_attribute : `bool`
            if True, returns a function that can be called with size as input

        Returns
        -------
        q : `numpy.ndarray`
            axis ratio of the elliptical lens galaxy

        Examples
        --------

        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_param_samplers=dict(axis_ratio="axis_ratio_uniform"))
        >>> print(od.axis_ratio_uniform(size=10))
        """

        identifier_dict = {'name': "axis_ratio_uniform"}
        param_dict = self.available_lens_prior_list_and_its_params["axis_ratio"]["axis_ratio_uniform"]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        low = param_dict["q_min"]
        high = param_dict["q_max"]
        q_rvs = njit(lambda size: np.random.uniform(
            low=low,
            high=high,
            size=size,
        ))
        if param_dict["q_max"] == param_dict["q_min"]:
            q_pdf = lambda q: 1. if q == param_dict["q_min"] else 0.
        else:
            q_pdf = lambda q: 1/(param_dict["q_max"]-param_dict["q_min"])

        q_object =FunctionConditioning(
            function=None,
            x_array=None,
            param_dict_given=param_dict,
            create_rvs=q_rvs,
            create_pdf=q_pdf,
            callback='rvs',
        )

        return q_object if get_attribute else q_object.rvs(size)
        
    def external_shear_normal(self, size, get_attribute=False, **kwargs):
        """
        Function to sample the external shear parameters from a normal distribution.

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        get_attribute : `bool`
            if True, returns a function that can be called with size as input

        Returns
        -------
        gamma_1 : `numpy.ndarray`
            shear component in the x-direction
        gamma_2 : `numpy.ndarray`
            shear component in the y-direction

        Examples
        --------

        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_param_samplers=dict(external_shear="external_shear_normal"))
        >>> print(od.external_shear_normal(size=10))
        """

        identifier_dict = {'name': "external_shear_normal"}
        param = self.available_lens_prior_list_and_its_params["external_shear"]["external_shear_normal"]
        if param:
            param.update(kwargs)
        else:
            param = kwargs
        identifier_dict.update(param)

        mean = param["mean"]
        std = param["std"]
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
            param_dict_given=identifier_dict,
            create_rvs=shear_rvs,
            create_pdf=shear_pdf,
            callback='rvs',
        )

        return shear_object if get_attribute else shear_object.rvs(size)
    
    def external_shear_numerical_hemanta(self, size, get_attribute=False, **kwargs):
        """
        Function to sample the external shear parameters from a normal distribution.

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        get_attribute : `bool`
            if True, returns a function that can be called with size as input

        Returns
        -------
        gamma_1 : `numpy.ndarray`
            shear component in the x-direction
        gamma_2 : `numpy.ndarray`
            shear component in the y-direction

        Examples
        --------

        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_param_samplers=dict(external_shear="external_shear_normal"))
        >>> print(od.external_shear_normal(size=10))
        """

        identifier_dict = {
            'name': 'external_shear_numerical_hemanta', 
            'external_shear_normal': {'mean': 0., 'std': 0.05},
        }

        gamma1, gamma2 = load_txt_from_module('ler', 'lens_galaxy_population.lens_param_data', 'external_shear_sl.txt')

        shear_object =FunctionConditioning(
            x_array=gamma1,
            y_array=gamma2,
            gaussian_kde=True,
            param_dict_given=identifier_dict,
            create_rvs=True,
            create_pdf=True,
            callback='rvs',
        )

        return shear_object if get_attribute else shear_object.rvs(size)
        
    def density_profile_slope_normal(self, size, get_attribute=False, **kwargs):
        """
        Function to sample the lens galaxy density profile slope with normal distribution.

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        get_attribute : `bool`
            if True, returns a function that can be used to sample velocity dispersion
        **kwargs : `dict`
            additional parameters to be passed to the function, 
            e.g. `mean` and `std` for the normal distribution

        Returns
        -------
        slope : `float`
            density profile slope of the lens galaxy

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_param_samplers=dict(density_profile_slope="density_profile_slope_normal"))
        >>> print(od.density_profile_slope_normal(size=10))
        """

        identifier_dict = {'name': "density_profile_slope_normal"}

        param = self.available_lens_prior_list_and_its_params["density_profile_slope"]["density_profile_slope_normal"]
        if param:
            param.update(kwargs)
        else:
            param = kwargs
        identifier_dict.update(param)

        mean = param["mean"]
        std = param["std"]
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
            param_dict_given=identifier_dict,
            create_rvs=slope_rvs,
            create_pdf=slope_pdf,
            callback='rvs',
        )

        if get_attribute:
            return slope_object
        else:
            return slope_rvs(size)
        
    def density_profile_slope_numerical_hemanta(self, size, get_attribute=False, **kwargs):
        """
        Function to sample the lens galaxy density profile slope with normal distribution.

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        get_attribute : `bool`
            if True, returns a function that can be used to sample velocity dispersion
        **kwargs : `dict`
            additional parameters to be passed to the function, 
            e.g. `mean` and `std` for the normal distribution

        Returns
        -------
        slope : `float`
            density profile slope of the lens galaxy

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_param_samplers=dict(density_profile_slope="density_profile_slope_normal"))
        >>> print(od.density_profile_slope_normal(size=10))
        """

        identifier_dict = {
            'name': 'density_profile_slope_numerical_hemanta', 
            'density_profile_slope_normal': {'mean': 1.99, 'std': 0.149},
        }

        gamma = load_txt_from_module('ler', 'lens_galaxy_population.lens_param_data', 'density_profile_slope_sl.txt')

        slope_object =FunctionConditioning(
            x_array=gamma,
            gaussian_kde=True,
            param_dict_given=identifier_dict,
            create_rvs=True,
            create_pdf=True,
            callback='rvs',
        )

        return slope_object if get_attribute else slope_object.rvs(size)
    
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
        identifier_dict = {'name': "lens_redshift_SDSS_catalogue_sis"}
        identifier_dict['z_min'] = self.z_min
        identifier_dict['z_max'] = self.z_max
        identifier_dict['cosmology'] = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        identifier_dict['velocity_dispersion'] = dict(sigma_min = 10., sigma_max = 420., alpha = 2.32, beta = 2.67, phistar = 8.0e-3*self.cosmo.h**3, sigmastar = 161.0)
        identifier_dict['resolution'] = self.create_new_interpolator["lens_redshift"]["resolution"]
        param_dict = self.available_lens_prior_list_and_its_params["lens_redshift"]["lens_redshift_SDSS_catalogue_sis"]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        zl_resolution = identifier_dict['resolution']
        x_array = np.linspace(0., 1., zl_resolution)
        pdf_ = lambda x: 30* x**2 * (1-x)**2  # Haris et al. 2018 (A7)

        zl_object = FunctionConditioning(
            function=pdf_,
            x_array=x_array,
            param_dict_given=identifier_dict,
            directory=self.directory,
            sub_directory="lens_redshift",
            name=identifier_dict['name'],
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

        return zl_object if get_attribute else zl_object.rvs(size, zs)
        
         
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

        identifier_dict = {'name': "velocity_dispersion_gengamma"}
        identifier_dict['z_min'] = self.z_min
        identifier_dict['z_max'] = self.z_max
        identifier_dict['cosmology'] = self.cosmo
        identifier_dict['resolution'] = self.create_new_interpolator["velocity_dispersion"]["resolution"]
        param_dict = self.available_lens_prior_list_and_its_params["velocity_dispersion"]["velocity_dispersion_gengamma"]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        # setting up inputs for the interpolator
        sigma_array = np.linspace(
            identifier_dict['sigma_min'],
            identifier_dict['sigma_max'],
            identifier_dict["resolution"]
        )
        pdf_func_ = lambda sigma_: gengamma.pdf(
            sigma_/identifier_dict['sigmastar'],
            a= identifier_dict['alpha']/identifier_dict['beta'],
            c= identifier_dict['beta'],
        )   # gengamma pdf

        sigma_object = FunctionConditioning(
            function=pdf_func_,
            x_array=sigma_array,
            param_dict_given=identifier_dict,
            directory=self.directory,
            sub_directory="velocity_dispersion",
            name=identifier_dict['name'],
            create_new=self.create_new_interpolator["velocity_dispersion"]['create_new'],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='rvs',
        )
        
        return sigma_object if get_attribute else sigma_object.rvs(size)
    
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

        identifier_dict = {'name': "velocity_dispersion_bernardi"}
        identifier_dict['z_min'] = self.z_min
        identifier_dict['z_max'] = self.z_max
        identifier_dict['cosmology'] = self.cosmo
        identifier_dict['name'] = "velocity_dispersion_bernardi"
        identifier_dict['resolution'] = self.create_new_interpolator["velocity_dispersion"]["resolution"]
        param_dict = self.available_lens_prior_list_and_its_params["velocity_dispersion"]["velocity_dispersion_bernardi"]
        if param_dict: 
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        # setting up inputs for the interpolator
        sigma_array = np.linspace(
            identifier_dict['sigma_min'],
            identifier_dict['sigma_max'],
            identifier_dict["resolution"]
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
            param_dict_given=identifier_dict,
            directory=self.directory,
            sub_directory="velocity_dispersion",
            name=identifier_dict['name'],
            create_new=self.create_new_interpolator["velocity_dispersion"]['create_new'],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='rvs',
        )
        
        return sigma_object if get_attribute else sigma_object.rvs(size)

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

        identifier_dict = {'name': "velocity_dispersion_ewoud"}
        identifier_dict['z_min'] = self.z_min
        identifier_dict['z_max'] = self.z_max
        identifier_dict['cosmology'] = self.cosmo
        identifier_dict['name'] = "velocity_dispersion_ewoud"
        identifier_dict['resolution'] = self.create_new_interpolator["velocity_dispersion"]["resolution"]
        param_dict = self.available_lens_prior_list_and_its_params["velocity_dispersion"]["velocity_dispersion_ewoud"].copy()
        if param_dict: 
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        # setting up inputs for the interpolator
        sigma_array = np.linspace(
            identifier_dict['sigma_min'],
            identifier_dict['sigma_max'],
            identifier_dict["resolution"]
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
            param_dict_given=identifier_dict,
            directory=self.directory,
            sub_directory="velocity_dispersion",
            name=identifier_dict['name'],
            create_new=self.create_new_interpolator["velocity_dispersion"]['create_new'],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='rvs',
        )
        
        return sigma_object if get_attribute else sigma_object.rvs(size, zl)

    ##################
    # Lens functions #
    ##################
    def optical_depth_numerical(self, zs, get_attribute=False, **kwargs):

        identifier_dict = {'name': "optical_depth_numerical"}
        identifier_dict['resolution'] = self.create_new_interpolator["optical_depth"]["resolution"]
        identifier_dict['z_min'] = self.z_min
        identifier_dict['z_max'] = self.z_max
        identifier_dict['cosmology'] = self.cosmo
        identifier_dict['velocity_dispersion'] = self.lens_param_samplers_params['velocity_dispersion']
        if identifier_dict['velocity_dispersion']: # if velocity_dispersion is not None
            identifier_dict['velocity_dispersion']['name'] = str(self.lens_param_samplers['velocity_dispersion'])
        identifier_dict['axis_ratio'] = self.lens_param_samplers_params['axis_ratio']
        if identifier_dict['axis_ratio']: # if axis_ratio is not None
            identifier_dict['axis_ratio']['name'] = str(self.lens_param_samplers['axis_ratio'])
        identifier_dict['axis_rotation_angle'] = self.lens_param_samplers_params['axis_rotation_angle']
        if identifier_dict['axis_rotation_angle']: # if axis_rotation_angle is not None
            identifier_dict['axis_rotation_angle']['name'] = str(self.lens_param_samplers['axis_rotation_angle'])
        identifier_dict['density_profile_slope'] = self.lens_param_samplers_params['density_profile_slope']
        if identifier_dict['density_profile_slope']: # if density_profile_slope is not None
            identifier_dict['density_profile_slope']['name'] = str(self.lens_param_samplers['density_profile_slope'])
        identifier_dict['external_shear'] = self.lens_param_samplers_params['external_shear']
        if identifier_dict['external_shear']: # if external_shear is not None
            identifier_dict['external_shear']['name'] = str(self.lens_param_samplers['external_shear'])
        
        param_dict = self.available_lens_functions_and_its_params["optical_depth"]["optical_depth_numerical"]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        z_min = self.z_min if self.z_min>0. else 0.001
        z_max = self.z_max
        resolution = identifier_dict['resolution']
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
            param_dict_given=identifier_dict,
            directory=self.directory,
            sub_directory="optical_depth",
            name=identifier_dict['name'],
            create_new=self.create_new_interpolator["optical_depth"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='function',
        )

        return tau_object if get_attribute else tau_object.function(zs)
        
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
            Einstein radii of the lens galaxies in radians. Multiply by 

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

        # setting up input parameters for interpolation
        identifier_dict = {'name': 'optical_depth_numerical',
            'resolution': 48,
            'z_min': 0.0,
            'z_max': 10.0,
            'cosmology': LambdaCDM(H0=70, Om0=0.3, Ode0=0.7),
            'velocity_dispersion': {'sigma_min': 100.0,
            'sigma_max': 400.0,
            'alpha': 0.94,
            'beta': 1.85,
            'phistar': 0.02099,
            'sigmastar': 113.78,
            'name': 'velocity_dispersion_ewoud'},
            'axis_ratio': {'q_min': 0.2, 'q_max': 1.0, 'name': 'axis_ratio_rayleigh'},
            'axis_rotation_angle': {'phi_min': 0.0,
            'phi_max': 6.283185307179586,
            'name': 'axis_rotation_angle_uniform'},
            'density_profile_slope': {'mean': 1.99,
            'std': 0.149,
            'name': 'density_profile_slope_normal'},
            'external_shear': {'mean': 0.0, 'std': 0.05, 'name': 'external_shear_normal'}
        }
        
        print("Getting pre-computed optical depth for EPL+Shear lenses...")

        zs_arr = np.geomspace(
            0.001,
            10,
            identifier_dict['resolution']
        )
        tau_arr = load_txt_from_module('ler', 'lens_galaxy_population.lens_param_data', 'optical_depth_epl_shear_vd_ewoud.txt')

        tau_object = FunctionConditioning(
            function=tau_arr,
            x_array=zs_arr,
            conditioned_y_array=None,
            param_dict_given=identifier_dict,
            directory=self.directory,
            sub_directory="optical_depth",
            name=identifier_dict['name'],
            create_new=self.create_new_interpolator["optical_depth"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='function',
        )

        return tau_object if get_attribute else tau_object.function(zs)

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

        # # z to luminosity_distance (luminosity_distance) conversion
        # #Dc = self.z_to_Dc(zs) * 1e-3  # 1e-3 converts Mpc to Gpc
        # splineDc = self.splineDc

        

        # if get_attribute:
        #     return tau
        # else:
        #     return tau(zs)
        identifier_dict = {'name': "optical_depth_sis_haris"}
        identifier_dict['z_min'] = self.z_min if self.z_min>0. else 0.001
        identifier_dict['z_max'] = self.z_max
        identifier_dict['cosmology'] = self.cosmo
        identifier_dict['resolution'] = self.create_new_interpolator["optical_depth"]["resolution"]
        param_dict = self.available_lens_functions_and_its_params["optical_depth"]["optical_depth_sis_haris"]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        zs_arr = np.geomspace(identifier_dict['z_min'], identifier_dict['z_max'], identifier_dict['resolution'])

        def tau(zs):
            Dc = self.comoving_distance.function(zs)
            return (Dc * 1e-3 / 62.2) ** 3
        
        tau_object = FunctionConditioning(
            function=tau,
            x_array=zs_arr,
            conditioned_y_array=None,
            param_dict_given=identifier_dict,
            directory=self.directory,
            sub_directory="optical_depth",
            name=identifier_dict['name'],
            create_new=self.create_new_interpolator["optical_depth"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='function',
        )

        return tau_object if get_attribute else tau_object.function(zs)

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
            print("Creating interpolated cross section function with sklearn NearestNeighbors, and this will be use for cross section calculation")

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
        resolution = self.create_new_interpolator["comoving_distance"]["resolution"]
        create_new = self.create_new_interpolator["comoving_distance"]["create_new"]
        zs = np.geomspace(z_min, z_max, resolution)
        Dc = self.cosmo.comoving_distance(zs).value  # co-moving distance in Mpc
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
            args = self.lens_functions_params["optical_depth"]
            if args is None:
                self._optical_depth = getattr(self, input_function)(zs=None, get_attribute=True)
            else:
                self._optical_depth = getattr(self, input_function)(zs=None, get_attribute=True, **args)
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
    
        if prior in self.available_lens_prior_list_and_its_params["velocity_dispersion"]:
            print(f"using ler available velocity dispersion function : {prior}")
            args = self.lens_param_samplers_params["velocity_dispersion"]
            if prior=="velocity_dispersion_choi":
                prior = "velocity_dispersion_bernardi"
            if args is None:
                self._velocity_dispersion = getattr(self, prior)(size=None, zl=None, get_attribute=True)
            else:
                self._velocity_dispersion = getattr(self, prior)(size=None, zl=None, get_attribute=True, **args)
        elif callable(prior):
            print("using user provided custom velocity_dispersion function")
            self._velocity_dispersion = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
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
            args = self.lens_param_samplers_params["axis_ratio"]
            if args is None:
                self._axis_ratio = getattr(self, prior)(size=None,sigma=None, get_attribute=True)
            else:
                self._axis_ratio = getattr(self, prior)(size=None, sigma=None, get_attribute=True, **args)
        elif callable(prior):
            print("using user provided custom axis_ratio function")
            self._axis_ratio = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
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
            args = self.lens_param_samplers_params["lens_redshift"]
            if args is None:
                self._lens_redshift = getattr(self, prior)(size=None, zs=None, get_attribute=True)
            else:
                self._lens_redshift = getattr(self, prior)(size=None, zs=None, get_attribute=True, **args)
        elif callable(prior):
            print("using user provided custom lens_redshift function")
            self._lens_redshift = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
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
            args = self.lens_param_samplers_params["axis_rotation_angle"]
            if args is None:
                self._axis_rotation_angle = getattr(self, prior)(size=None, get_attribute=True)
            else:
                self._axis_rotation_angle = getattr(self, prior)(size=None, get_attribute=True, **args)
        elif callable(prior):
            print("using user provided custom axis_rotation_angle function")
            self._axis_rotation_angle = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
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
            args = self.lens_param_samplers_params["external_shear"]
            if args is None:
                self._external_shear = getattr(self, prior)(size=None, get_attribute=True)
            else:
                self._external_shear = getattr(self, prior)(size=None, get_attribute=True, **args)
        elif callable(prior):
            print("using user provided custom external_shear function")
            self._external_shear = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
        elif isinstance(prior, object):
            print("using user provided custom external_shear class/object")
            self._external_shear = prior
        else:
            raise ValueError("prior not in available_lens_prior_list_and_its_params['external_shear']")
        
    @property
    def external_shear_sl(self):
        return self._external_shear_sl
    
    @external_shear_sl.setter
    def external_shear_sl(self, prior):
        if prior in self.available_lens_prior_list_and_its_params["external_shear_sl"]:
            print(f"using ler available external_shear_sl function : {prior}")
            args = self.lens_param_samplers_params["external_shear_sl"]
            if args is None:
                self._external_shear_sl = getattr(self, prior)(size=None, get_attribute=True)
            else:
                self._external_shear_sl = getattr(self, prior)(size=None, get_attribute=True, **args)
        elif callable(prior):
            print("using user provided custom external_shear_sl function")
            self._external_shear_sl = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
        elif isinstance(prior, object):
            print("using user provided custom external_shear_sl class/object")
            self._external_shear_sl = prior
        else:
            raise ValueError("prior not in available_lens_prior_list_and_its_params['external_shear_sl']")
        
    @property
    def density_profile_slope(self):
        return self._density_profile_slope
    
    @density_profile_slope.setter
    def density_profile_slope(self, prior):
        if prior in self.available_lens_prior_list_and_its_params["density_profile_slope"]:
            print(f"using ler available density_profile_slope function : {prior}")
            args = self.lens_param_samplers_params["density_profile_slope"]
            if args is None:
                self._density_profile_slope = getattr(self, prior)(size=None, get_attribute=True)
            else:
                self._density_profile_slope = getattr(self, prior)(size=None, get_attribute=True, **args)
        elif callable(prior):
            print("using user provided custom density_profile_slope function")
            self._density_profile_slope = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
        elif isinstance(prior, object):
            print("using user provided custom density_profile_slope class/object")
            self._density_profile_slope = prior
        else:
            raise ValueError("prior not in available_lens_prior_list_and_its_params['density_profile_slope']")
        
    @property
    def density_profile_slope_sl(self):
        return self._density_profile_slope_sl
    
    @density_profile_slope_sl.setter
    def density_profile_slope_sl(self, prior):
        if prior in self.available_lens_prior_list_and_its_params["density_profile_slope_sl"]:
            print(f"using ler available density_profile_slope_sl function : {prior}")
            args = self.lens_param_samplers_params["density_profile_slope_sl"]
            if args is None:
                self._density_profile_slope_sl = getattr(self, prior)(size=None, get_attribute=True)
            else:
                self._density_profile_slope_sl = getattr(self, prior)(size=None, get_attribute=True, **args)
        elif callable(prior):
            print("using user provided custom density_profile_slope_sl function")
            self._density_profile_slope_sl = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
        elif isinstance(prior, object):
            print("using user provided custom density_profile_slope_sl class/object")
            self._density_profile_slope_sl = prior
        else:
            raise ValueError("prior not in available_lens_prior_list_and_its_params['density_profile_slope_sl']")
        

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
                    sigma_min=100., sigma_max=400., alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78
                ),
                velocity_dispersion_choi = dict(
                    sigma_min=100., sigma_max=400., alpha = 2.32, beta = 2.67, phistar = 8.0e-3*self.cosmo.h**3, sigmastar = 161.0
                ),
                velocity_dispersion_bernardi = dict(
                    sigma_min=100., sigma_max=400., alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78
                ),
                velocity_dispersion_ewoud = dict(
                    sigma_min=100., sigma_max=400., alpha=0.94, beta=1.85, phistar=2.099e-2*(self.cosmo.h/0.7)**3, sigmastar=113.78
                ),
            ),
            axis_ratio = dict(
                axis_ratio_rayleigh = dict(q_min=0.2, q_max=1.0),
                axis_ratio_padilla_strauss = dict(q_min=0.2, q_max=1.0),
                axis_ratio_uniform = dict(q_min=0.2, q_max=1.0),
            ),
            axis_rotation_angle = dict(
                axis_rotation_angle_uniform = dict(
                    phi_min=0.0, phi_max=2 * np.pi
                ),
            ),
            external_shear = dict(
                external_shear_normal = dict(mean=0., std=0.05),
            ),
            external_shear_sl = dict(
                external_shear_normal = dict(mean=0., std=0.05, name="external_shear_normal_sl"),
                external_shear_sl_numerical_hemanta = dict(external_shear_normal=dict(mean=0., std=0.05), name="external_shear_sl_numerical_hemanta"),
            ),
            density_profile_slope = dict(
                density_profile_slope_normal=dict(mean=1.99, std=0.149),
            ),
            density_profile_slope_sl = dict(
                density_profile_slope_normal=dict(mean=2.091, std=0.133, name="density_profile_slope_normal_sl"),
                density_profile_slope_sl_numerical_hemanta=dict(density_profile_slope_normal=dict(mean=1.99, std=0.149), name="density_profile_slope_sl_numerical_hemanta"),
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
                rjs_with_cross_section_mp=None,
            ),
            optical_depth=dict(
                optical_depth_sis_haris=None,
                optical_depth_sie_hemanta=None,
                optical_depth_epl_shear_hemanta=None,
                optical_depth_numerical=None,
            ),
            param_sampler_type=dict(
                sample_all_routine_sie_sl=None,
                sample_all_routine_sis_sl=None,
                sample_all_routine_epl_shear_sl=None,
            ),
            cross_section=dict(
                cross_section_sie_feixu=None,
                cross_section_sis=None,
                cross_section_caustic_area=None,
                interpolated_cross_section_function=None,
            ),
        )

        return self._available_lens_functions_and_its_params









