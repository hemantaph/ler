from numba import njit
from multiprocessing import Pool
import numpy as np
from scipy.integrate import quad
from scipy.stats import gengamma
from scipy.interpolate import CubicSpline
from astropy.cosmology import LambdaCDM
from tqdm import tqdm

from ..utils import (
    cubic_spline_interpolator,
    inverse_transform_sampler,
    cubic_spline_interpolator2d_array,
    save_json,
    load_json,
    interpolator_json_path,
    FunctionConditioning,
    inverse_transform_sampler2d,
    pdf_cubic_spline_interpolator2d_array,
    normal_pdf,
    normal_pdf_2d,
    comoving_distance,
    angular_diameter_distance,
    angular_diameter_distance_z1z2,
    differential_comoving_volume,
    is_njitted,
    redshift_optimal_spacing,
)

from .jit_functions import (
    phi_cut_SIE,
    phi,
    phi_loc_bernardi,
    phi_q2_ellipticity_hemanta,
    axis_ratio_rayleigh_pdf,
)

from .mp import cross_section_mp


class OpticalDepth:
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
        Directory where the interpolators are saved (default is './interpolator_json').
        If True, creates a new interpolator (default is False).
    verbose : bool, optional
        If True, prints additional information during initialization (default is False).

    Raises
    ------
    ValueError
        If `lens_type` is not in ['sie_galaxy', 'epl_shear_galaxy', 'sis_galaxy'].
    """

    def __init__(
        self,
        npool=4,
        z_min=0.0,
        z_max=10.0,
        cosmology=None,
        lens_type="epl_shear_galaxy",
        lens_functions=None,
        lens_functions_params=None,
        lens_param_samplers=None,
        lens_param_samplers_params=None,
        directory="./interpolator_json",
        create_new_interpolator=False,
        verbose=False,
    ):
    
        print("\nInitializing OpticalDepth class\n")
        if lens_type not in ["sie_galaxy", "epl_shear_galaxy", "sis_galaxy"]:
            raise ValueError(
                "lens_type not in ['sie_galaxy', 'epl_shear_galaxy', 'sis_galaxy']"
            )
        self.lens_type = lens_type

        self.npool = npool
        self.z_min = z_min
        self.z_max = z_max
        self.cosmo = cosmology if cosmology else LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        self.directory = directory

        # decision dict fot creating interpolator
        self.initialize_decision_dictionary(create_new_interpolator, lens_type)
        # default lens functions and samplers
        (
            self.lens_param_samplers,
            self.lens_param_samplers_params,
            self.lens_functions,
            self.lens_functions_params,
        ) = self.default_lens_samplers_and_functions(lens_type)

        # change default lens functions and samplers if user input is given
        # if lens functions and samplers are given but not parameters, then use default parameters
        # it updates self.lens_param_samplers, self.lens_param_samplers_params, self.lens_functions, self.lens_functions_params
        self.lens_functions_and_sampler_categorization(
            lens_param_samplers,
            lens_param_samplers_params,
            lens_functions,
            lens_functions_params,
        )

        # Initialize cosmological functions
        self.comoving_distance = comoving_distance(z_min=self.z_min, z_max=self.z_max, cosmo=self.cosmo, directory=self.directory, create_new=self.create_new_interpolator["comoving_distance"]["create_new"], resolution=self.create_new_interpolator["comoving_distance"]["resolution"], get_attribute=True)
        self.angular_diameter_distance = angular_diameter_distance(z_min=self.z_min, z_max=self.z_max, cosmo=self.cosmo, directory=self.directory, create_new=self.create_new_interpolator["angular_diameter_distance"]["create_new"], resolution=self.create_new_interpolator["angular_diameter_distance"]["resolution"], get_attribute=True)
        self.angular_diameter_distance_z1z2 = angular_diameter_distance_z1z2(z_min=self.z_min, z_max=self.z_max, cosmo=self.cosmo, directory=self.directory, create_new=self.create_new_interpolator["angular_diameter_distance_z1z2"]["create_new"], resolution=self.create_new_interpolator["angular_diameter_distance_z1z2"]["resolution"], get_attribute=True)
        self.differential_comoving_volume = differential_comoving_volume(z_min=self.z_min, z_max=self.z_max, cosmo=self.cosmo, directory=self.directory, create_new=self.create_new_interpolator["differential_comoving_volume"]["create_new"], resolution=self.create_new_interpolator["differential_comoving_volume"]["resolution"], get_attribute=True)

        # lens sampler initialization
        self.velocity_dispersion = self.lens_param_samplers["velocity_dispersion"]
        self.axis_ratio = self.lens_param_samplers["axis_ratio"]
        self.axis_rotation_angle = self.lens_param_samplers["axis_rotation_angle"]
        self.density_profile_slope = self.lens_param_samplers["density_profile_slope"]
        self.external_shear = self.lens_param_samplers["external_shear"]

        # cross section initialization 
        self.cross_section = self.lens_functions["cross_section"]
        
        # lens redshift initialization
        self.lens_redshift = self.lens_param_samplers["lens_redshift"]
        self.lens_redshift_intrinsic = self.lens_redshift_intrinsic(
            size=None, get_attribute=True
        )
        # lens parameter initialization 
        # density profile slope initialization
        self.density_profile_slope_sl = self.lens_param_samplers[
            "density_profile_slope_sl"
        ]
        # external shear initialization
        self.external_shear_sl = self.lens_param_samplers["external_shear_sl"]
        # lens function initialization
        self.optical_depth = self.lens_functions["optical_depth"]


    #-----------------
    # Initialization 
    #-----------------
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
        lens_param_samplers_ : `dict`
            dictionary of priors
        lens_param_samplers_params_ : `dict`
            dictionary of priors parameters
        lens_sampler_names_ : `dict`
            dictionary of sampler names
        lens_functions_ : `dict`
            dictionary of lens functions
        """

        if lens_type == "epl_shear_galaxy":
            lens_param_samplers = dict(
                source_redshift_sl="strongly_lensed_source_redshifts",
                lens_redshift="lens_redshift_strongly_lensed_numerical",
                velocity_dispersion="velocity_dispersion_ewoud",
                axis_ratio="axis_ratio_rayleigh",
                axis_rotation_angle="axis_rotation_angle_uniform",
                external_shear="external_shear_normal",
                density_profile_slope="density_profile_slope_normal",
                external_shear_sl="external_shear_normal",
                density_profile_slope_sl="density_profile_slope_normal",
            )
            lens_param_samplers_params = dict(
                source_redshift_sl=None,
                lens_redshift=dict(integration_size=20000),
                velocity_dispersion=dict(
                    sigma_min=100.0,
                    sigma_max=400.0,
                    alpha=0.94,
                    beta=1.85,
                    phistar=2.099e-2 * (self.cosmo.h / 0.7) ** 3,
                    sigmastar=113.78,
                ),
                axis_ratio=dict(q_min=0.2, q_max=1.0),
                axis_rotation_angle=dict(phi_min=0.0, phi_max=2 * np.pi),
                external_shear=dict(mean=0.0, std=0.05),
                density_profile_slope=dict(mean=1.99, std=0.149),
                external_shear_sl=dict(mean=0.0, std=0.05),
                density_profile_slope_sl=dict(mean=2.078, std=0.16),
            )
            lens_functions = dict(
                param_sampler_type="sample_all_routine_epl_shear_sl",
                cross_section_based_sampler="importance_sampling_with_cross_section",
                # cross_section_based_sampler="rejection_sampling_with_cross_section",
                optical_depth="optical_depth_numerical",
                cross_section="cross_section_epl_shear_interpolation",
            )
            lens_functions_params = dict(
                param_sampler_type=None,
                cross_section_based_sampler=dict(n_prop=200),
                # cross_section_based_sampler=dict(saftey_factor=1.2),
                optical_depth=None,
                cross_section=None,
            )
        elif lens_type == "sie_galaxy":
            lens_param_samplers = dict(
                source_redshift_sl="strongly_lensed_source_redshifts",
                lens_redshift="lens_redshift_strongly_lensed_numerical",
                velocity_dispersion="velocity_dispersion_ewoud",
                axis_ratio="axis_ratio_rayleigh",
                axis_rotation_angle="axis_rotation_angle_uniform",
                external_shear="external_shear_normal",
                density_profile_slope="density_profile_slope_normal",
                external_shear_sl="external_shear_normal",
                density_profile_slope_sl="density_profile_slope_normal",
            )
            lens_param_samplers_params = dict(
                source_redshift_sl=None,
                lens_redshift=dict(integration_size=20000),
                velocity_dispersion=dict(
                    sigma_min=100.0,
                    sigma_max=400.0,
                    alpha=0.94,
                    beta=1.85,
                    phistar=2.099e-2 * (self.cosmo.h / 0.7) ** 3,
                    sigmastar=113.78,
                ),
                axis_ratio=dict(q_min=0.2, q_max=1.0),
                axis_rotation_angle=dict(phi_min=0.0, phi_max=2 * np.pi),
                external_shear=dict(mean=0.0, std=0.0),
                density_profile_slope=dict(mean=2.0, std=0.0),
                external_shear_sl=dict(mean=0.0, std=0.0),
                density_profile_slope_sl=dict(mean=2.0, std=0.0),
            )
            lens_functions = dict(
                param_sampler_type="sample_all_routine_epl_shear_sl",
                cross_section_based_sampler="importance_sampling_with_cross_section",
                optical_depth="optical_depth_numerical",
                cross_section="cross_section_sie_feixu",
            )
            lens_functions_params = dict(
                param_sampler_type=None,
                cross_section_based_sampler=dict(n_prop=200),
                optical_depth=dict(interpolated_cross_section=True),
                cross_section=None,
            )
        elif lens_type == "sis_galaxy":
            lens_param_samplers = dict(
                source_redshift_sl="strongly_lensed_source_redshifts",
                lens_redshift="lens_redshift_sis_haris",
                velocity_dispersion="velocity_dispersion_choi",
                axis_ratio="axis_ratio_uniform",
                axis_rotation_angle="axis_rotation_angle_uniform",
                external_shear="external_shear_normal",
                density_profile_slope="density_profile_slope_normal",
                external_shear_sl="external_shear_normal",
                density_profile_slope_sl="density_profile_slope_normal",
            )
            lens_param_samplers_params = dict(
                source_redshift_sl=None,
                lens_redshift=dict(integration_size=20000),
                velocity_dispersion=dict(
                    sigma_min=100.0,
                    sigma_max=400.0,
                    alpha=2.32,
                    beta=2.67,
                    phistar=8.0e-3 * self.cosmo.h**3,
                    sigmastar=161.0,
                ),
                axis_ratio=dict(q_min=1.0, q_max=1.0),
                axis_rotation_angle=dict(phi_min=0.0, phi_max=0.0),
                external_shear=dict(mean=0.0, std=0.0),
                density_profile_slope=dict(mean=2.0, std=0.0),
                external_shear_sl=dict(mean=0.0, std=0.0),
                density_profile_slope_sl=dict(mean=2.0, std=0.0),
            )
            lens_functions = dict(
                param_sampler_type="sample_all_routine_epl_shear_sl",
                cross_section_based_sampler="importance_sampling_with_cross_section",
                optical_depth="optical_depth_sis_haris",
                cross_section="cross_section_sis",
            )
            lens_functions_params = dict(
                param_sampler_type=None,
                cross_section_based_sampler=dict(n_prop=200),
                optical_depth=dict(interpolated_cross_section=True),
                cross_section=None,
            )
        else:
            raise ValueError("lens_type not recognized")

        return (
            lens_param_samplers,
            lens_param_samplers_params,
            lens_functions,
            lens_functions_params,
        )

    def initialize_decision_dictionary(self, create_new_interpolator, lens_type):
        """
        Function to initialize decision dictionary for creating interpolator

        Parameters
        ----------
        create_new_interpolator : `dict` or `bool`
            dictionary to create new interpolator for velocity dispersion and optical depth.
        """

        # initialize the interpolator's parameters
        self.create_new_interpolator = dict(
            velocity_dispersion=dict(create_new=False, resolution=500, zl_resolution=48),
            axis_ratio=dict(create_new=False, resolution=500, sigma_resolution=48),
            lens_redshift=dict(create_new=False, resolution=48, zl_resolution=48),
            lens_redshift_intrinsic=dict(create_new=False, resolution=500),
            optical_depth=dict(create_new=False, resolution=48),
            comoving_distance=dict(create_new=False, resolution=500),
            angular_diameter_distance=dict(create_new=False, resolution=500),
            angular_diameter_distance_z1z2=dict(create_new=False, resolution=500),
            differential_comoving_volume=dict(create_new=False, resolution=500),
            density_profile_slope=dict(create_new=False, resolution=100),
            lens_parameters_kde_sl=dict(create_new=False, resolution=5000),
        )
        if lens_type == "sis_galaxy":
            self.create_new_interpolator.update(
                cross_section=dict(create_new=False, resolution=[5, 5, 5, 5, 5]),
            )
        elif lens_type == "sie_galaxy":
            self.create_new_interpolator.update(
                cross_section=dict(create_new=False, resolution=[25, 25, 5, 5, 5]),
            )
        elif lens_type == "epl_shear_galaxy":
            self.create_new_interpolator.update(
                cross_section=dict(create_new=False, resolution=[25, 25, 45, 15, 15]),
            )

        if isinstance(create_new_interpolator, dict):
            self.create_new_interpolator.update(create_new_interpolator)
        # if create_new_interpolator is True, create new interpolator for all
        elif create_new_interpolator is True:
            for key in self.create_new_interpolator.keys():
                self.create_new_interpolator[key]["create_new"] = True

    def lens_functions_and_sampler_categorization(
        self,
        lens_param_samplers,
        lens_param_samplers_params,
        lens_functions,
        lens_functions_params,
    ):
        """
        Function to initialize velocity dispersion sampler with it's settings. The reason I am seperating this from lens_param_samplers_categorization is only a specific parameters needs special attention.

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

        sampler_prior_names = [
            "lens_redshift",
            "velocity_dispersion",
            "axis_ratio",
            "axis_rotation_angle",
            "external_shear",
            "density_profile_slope",
            "external_shear_sl",
            "density_profile_slope_sl",
        ]

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
                    dict_ = self.available_lens_prior_list_and_its_params[
                        name
                    ]  # e.g. {'axis_ratio_padilla_strauss': {'q_min': 0.2, 'q_max': 1.0}, ....}
                    if sampler_name in dict_:  # e.g. 'axis_ratio_padilla_strauss'
                        param_dict = dict_[
                            sampler_name
                        ]  # e.g. {'q_min': 0.2, 'q_max': 1.0}
                        if (lens_param_samplers_params is None) or (
                            lens_param_samplers_params[name] is None
                        ):  # not a dictionary
                            self.lens_param_samplers_params[name] = param_dict
                        else:  # if there is user input lens_param_samplers_params
                            param_dict.update(lens_param_samplers_params[name])
                            self.lens_param_samplers_params[name].update(
                                param_dict
                            )  # user inputs will override the default values
                    else:
                        raise ValueError(
                            f"{name} sampler {sampler_name} not available.\n Available {name} samplers and its parameters are: {dict_[name]}"
                        )
                elif not callable(lens_param_samplers[name]):
                    raise ValueError(
                        f"Given {name} sampler should be either a string name of available sampler or a function"
                    )

        lens_function_names = ["optical_depth", "cross_section"]

        # if there is user input function, update the sampler priors
        for name in lens_function_names:
            if (lens_functions is not None) and (name in lens_functions):
                function_name = lens_functions[name]
                if isinstance(function_name, str):
                    # available lens functions for name e.g. 'optical_depth'
                    dict_ = self.available_lens_functions_and_its_params[name]
                    if function_name in dict_:
                        param_dict = dict_[function_name]
                        if (lens_functions_params is None) or (
                            lens_functions_params[name] is None
                        ):  # not a dictionary
                            self.lens_functions_params[name] = param_dict
                        else:  # if there is user input lens_functions_params
                            param_dict.update(lens_functions_params[name])
                            self.lens_functions_params[name].update(
                                param_dict
                            )  # user inputs will override the default values
                    else:
                        raise ValueError(
                            f"{name} function {function_name} not available.\n Available {name} functions and its parameters are: {dict_[name]}"
                        )
                elif not callable(lens_functions[name]):
                    raise ValueError(
                        f"Given {name} function should be either a string name of available function or a function"
                    )

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
        >>> print(self.axis_ratio(sigma=200.))
        """

        identifier_dict = {"name": "axis_ratio_rayleigh"}

        identifier_dict["velocity_dispersion"] = {}
        if (
            self.velocity_dispersion.info is not None
        ):  # if velocity_dispersion is not None
            for key, value in self.velocity_dispersion.info.items():
                identifier_dict["velocity_dispersion"][key] = value

        identifier_dict["resolution"] = self.create_new_interpolator["axis_ratio"][
            "resolution"
        ]
        identifier_dict["sigma_resolution"] = self.create_new_interpolator["axis_ratio"][
            "sigma_resolution"
        ]
        param_dict = self.available_lens_prior_list_and_its_params["axis_ratio"][
            "axis_ratio_rayleigh"
        ]
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
            identifier_dict["velocity_dispersion"]["sigma_min"],
            identifier_dict["velocity_dispersion"]["sigma_max"],
            identifier_dict["sigma_resolution"],
        )

        q_pdf = lambda q, sigma: axis_ratio_rayleigh_pdf(  # noqa: E731
            q=q,
            sigma=sigma,
            q_min=identifier_dict["q_min"],
            q_max=identifier_dict["q_max"],
        )

        q_object = FunctionConditioning(
            function=q_pdf,
            x_array=q_array,
            conditioned_y_array=sigma_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="axis_ratio",
            name="axis_ratio_rayleigh",
            create_new=self.create_new_interpolator["axis_ratio"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
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
        >>> print(self.axis_ratio(size=10))
        """

        identifier_dict = {"name": "axis_ratio_padilla_strauss"}
        identifier_dict["resolution"] = self.create_new_interpolator["axis_ratio"][
            "resolution"
        ]
        param_dict = self.available_lens_prior_list_and_its_params["axis_ratio"][
            "axis_ratio_padilla_strauss"
        ]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        # Using Padilla and Strauss 2008 distribution for axis ratio
        q_array = np.array(
            [
                0.04903276402927845,
                0.09210526315789469,
                0.13596491228070173,
                0.20789473684210524,
                0.2899703729522482,
                0.3230132450331126,
                0.35350877192982455,
                0.37946148483792264,
                0.4219298245614036,
                0.4689525967235971,
                0.5075026141512723,
                0.5226472638550018,
                0.5640350877192983,
                0.6096491228070177,
                0.6500000000000001,
                0.6864848379226213,
                0.7377192982456142,
                0.7787295224817011,
                0.8007581038689441,
                0.822786685256187,
                0.8668438480306729,
                0.8973684210526317,
                0.9254385964912283,
            ]
        )
        pdf = np.array(
            [
                0.04185262687135349,
                0.06114520695141845,
                0.096997499638376,
                0.1932510900336828,
                0.39547914337673706,
                0.49569751276216234,
                0.6154609137685201,
                0.7182049959882812,
                0.920153741243567,
                1.1573982157399754,
                1.3353263628106684,
                1.413149656448315,
                1.5790713532948977,
                1.7280185150744938,
                1.8132994441344819,
                1.8365803753840484,
                1.8178662203211204,
                1.748929843583365,
                1.688182592496342,
                1.6274353414093188,
                1.4948487090314488,
                1.402785526832393,
                1.321844068356993,
            ]
        )

        q_object = FunctionConditioning(
            function=pdf,
            x_array=q_array,
            conditioned_y_array=None,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="axis_ratio",
            name="axis_ratio_padilla_strauss",
            create_new=self.create_new_interpolator["axis_ratio"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
        )

        return q_object if get_attribute else q_object(size)

    def lens_redshift_strongly_lensed_numerical(
        self, size=1000, zs=None, get_attribute=False, **kwargs
    ):
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
        >>> od = OpticalDepth(lens_param_samplers=dict(lens_redshift="lens_redshift_strongly_lensed_numerical"))
        >>> print(self.lens_redshift(size=10, zs=1.0))
        """

        identifier_dict = {}
        identifier_dict["name"] = (
            "lens_redshift_strongly_lensed_numerical_" + self.lens_type
        )
        identifier_dict["resolution"] = self.create_new_interpolator["lens_redshift"][
            "resolution"
        ]
        identifier_dict["zl_resolution"] = self.create_new_interpolator["lens_redshift"][
            "zl_resolution"
        ]
        identifier_dict["integration_size"] = self.lens_param_samplers_params[
            "lens_redshift"
        ]["integration_size"]

        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo

        name_list = [
            "z_min",
            "z_max",
            "cosmology",
            "velocity_dispersion",
            "axis_ratio",
            "axis_rotation_angle",
            "density_profile_slope",
            "external_shear",
            # "resolution",
        ]

        identifier_dict["cross_section"] = {}
        identifier_dict["cross_section"]["name"] = self.lens_functions["cross_section"]
        if self.lens_functions_params["cross_section"] is not None:
            for key, value in self.lens_functions_params["cross_section"].items():
                if key not in name_list:
                    identifier_dict["cross_section"][key] = value

        identifier_dict["velocity_dispersion"] = {}
        if (
            self.velocity_dispersion.info is not None
        ):  # if velocity_dispersion is not None
            for key, value in self.velocity_dispersion.info.items():
                if key not in name_list:
                    identifier_dict["velocity_dispersion"][key] = value

        identifier_dict["axis_ratio"] = {}
        if self.axis_ratio.info is not None:  # if axis_ratio is not None
            for key, value in self.axis_ratio.info.items():
                if key not in name_list:
                    identifier_dict["axis_ratio"][key] = value

        identifier_dict["axis_rotation_angle"] = {}
        if (
            self.axis_rotation_angle.info is not None
        ):  # if axis_rotation_angle is not None
            for key, value in self.axis_rotation_angle.info.items():
                if key not in name_list:
                    identifier_dict["axis_rotation_angle"][key] = value

        identifier_dict["density_profile_slope"] = {}
        if (
            self.density_profile_slope.info is not None
        ):  # if density_profile_slope is not None
            for key, value in self.density_profile_slope.info.items():
                if key not in name_list:
                    identifier_dict["density_profile_slope"][key] = value

        identifier_dict["external_shear"] = {}
        if self.external_shear.info is not None:  # if external_shear is not None
            for key, value in self.external_shear.info.items():
                if key not in name_list:
                    identifier_dict["external_shear"][key] = value

        param_dict = self.available_lens_prior_list_and_its_params["lens_redshift"][
            "lens_redshift_strongly_lensed_numerical"
        ]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        print("Numerically solving the lens redshift distribution...")
        zs_resolution = identifier_dict["resolution"]
        zs_min = self.z_min + 0.001 if self.z_min == 0.0 else self.z_min
        zs_max = self.z_max
        zs_array = redshift_optimal_spacing(zs_min, zs_max, zs_resolution)

        _, it_exist = interpolator_json_path(
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="lens_redshift",
            interpolator_name=identifier_dict["name"],
        )

        # zl dependent velocity dispersion distribution
        zl_scaled = []
        zl_min = 0.0001
        for i, zs_ in enumerate(zs_array):
            zl_resolution = identifier_dict["zl_resolution"]
            buffer_ = np.linspace(zl_min, zs_ - zl_min, zl_resolution)
            zl_scaled.append(buffer_ / zs_)
        zl_scaled = np.array(zl_scaled)

        create_new = self.create_new_interpolator["lens_redshift"]["create_new"]
        if not it_exist or create_new:
            number_density = self._helper_number_density_calculation(zl_scaled, zs_array)
            # number density is zero for zl=0 and infinite for zl=zs
            # Adding zero at the first element of each row
            zl_scaled = np.hstack((np.zeros((zl_scaled.shape[0], 1)), zl_scaled))
            number_density = np.hstack(
                (np.zeros((number_density.shape[0], 1)), number_density)
            )
            # Adding one at the last element of each row of zl_scaled
            zl_scaled = np.hstack((zl_scaled, np.ones((zl_scaled.shape[0], 1))))
            # Adding zero at the last element of each row of density
            number_density = np.hstack(
                (number_density, np.zeros((number_density.shape[0], 1)))
            )
        else:
            # if it_exist, it will use the saved interpolator
            number_density = None

        zl_object = FunctionConditioning(
            function=number_density,
            x_array=zl_scaled,
            conditioned_y_array=zs_array,
            non_zero_function=True,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="lens_redshift",
            name=identifier_dict["name"],
            create_new=create_new,
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
        )

        # un-scaled lens redshift
        cdf_values = zl_object.cdf_values
        x_array = zl_object.x_array
        y_array = zl_object.conditioned_y_array
        function_spline = zl_object.function_spline
        pdf_norm_const = zl_object.pdf_norm_const

        zl_object.function = njit(
            lambda x, y: cubic_spline_interpolator2d_array(
                x / y, y, function_spline, x_array, y_array
            )
        )

        zl_object.pdf = njit(
            lambda x, y: pdf_cubic_spline_interpolator2d_array(
                x / y, y, pdf_norm_const, function_spline, x_array, y_array
            )
            / y
        )

        zl_object.rvs = njit(
            lambda size, y: inverse_transform_sampler2d(
                size, y, cdf_values, x_array, y_array
            )
            * y
        )

        return zl_object if get_attribute else zl_object(size, zs)

    def _helper_number_density_calculation(self, zl_scaled, zs):
        """
        Compute the lens redshift distribution using multithreaded njitted function. This function is used when the velocity dispersion distribution depends on lens redshift.

        Parameters
        ----------
        zl_scaled : array_like`
            1D array of lens redshifts, scaled by the source redshift.
        zs : array_like
            1D array of source redshifts.
        zl_distribution_name : str
            Name of the lens redshift distribution to compute.

        Returns
        -------
        density_array : array_like
            2D array of the lens redshift distribution.
        """

        # Check if the function can be JIT compiled
        print("Checking if the function can be JIT compiled")
        _, use_njit_sampler = self._lens_redshift_multithreaded_njit(np.array([zl_scaled[0]]), np.array([zs[0]]))

        # if use_njit_sampler:
        print("Using multithreaded njit")
        density_array, _ = self._lens_redshift_multithreaded_njit(zl_scaled, zs)
        # else:
        # print("Using multiprocessing")
        # density_array = self._lens_redshift_multiprocessing(zl_scaled, zs)

        return density_array

    def _lens_redshift_multithreaded_njit(self, zl_scaled, zs):

        # Get random variable samplers from the initialized objects
        sigma_rvs_ = self.velocity_dispersion.rvs
        q_rvs_ = self.axis_ratio.rvs
        phi_rvs = self.axis_rotation_angle.rvs
        gamma_rvs = self.density_profile_slope.rvs
        shear_rvs = self.external_shear.rvs

        number_density_ = self.velocity_dispersion.function

        cross_section_ = self.cross_section
        # check the number of inputs for cross_section function
        # and JIT compile the wrapper function if necessary
        if cross_section_.__code__.co_argcount == 4:
            cross_section = njit(lambda zs, zl, sigma, q, phi, gamma, gamma1, gamma2: cross_section_(zs, zl, sigma, q))
        elif cross_section_.__code__.co_argcount == 3:
            cross_section = njit(lambda zs, zl, sigma, q, phi, gamma, gamma1, gamma2: cross_section_(zs, zl, sigma))
        else:
            cross_section = cross_section_
        
        # JIT compile the random variable samplers and PDFs if they are not already JIT compiled
        # This is required for compatibility with the JIT compiled cross_section_based_sampler
        if self.velocity_dispersion.conditioned_y_array is None:
            if is_njitted(sigma_rvs_):
                sigma_rvs = njit(lambda size, zl: sigma_rvs_(size))
                number_density = njit(lambda sigma, zl: number_density_(sigma))
            else:
                sigma_rvs = lambda size, zl: sigma_rvs_(size)  # noqa: E731
                number_density = lambda sigma, zl: number_density_(sigma)  # noqa: E731
        else:
            sigma_rvs = sigma_rvs_
            number_density = number_density_
            
        if self.axis_ratio.conditioned_y_array is None:
            if is_njitted(q_rvs_):
                q_rvs = njit(lambda size, sigma: q_rvs_(size))
            else:
                q_rvs = lambda size, sigma: q_rvs_(size)  # noqa: E731
        else:
            q_rvs = q_rvs_

        # List of functions to check if they are JIT compiled
        dict_ = {
            'sigma_rvs': sigma_rvs,
            'q_rvs': q_rvs,
            'phi_rvs': phi_rvs,
            'gamma_rvs': gamma_rvs,
            'shear_rvs': shear_rvs,
            'number_density': number_density,
            'cross_section': cross_section
        }

        # Check if all functions are JIT compiled
        use_njit_sampler = True
        for key, value in dict_.items():
            if not is_njitted(value):
                print(f"Warning: {key} is not njitted.")
                use_njit_sampler = False
            
        sigma_min = self.lens_param_samplers_params['velocity_dispersion']['sigma_min']
        sigma_max = self.lens_param_samplers_params['velocity_dispersion']['sigma_max']

        # integration_size
        integration_size = (
            self.lens_param_samplers_params["lens_redshift"]["integration_size"]
            if "integration_size" in self.lens_param_samplers_params["lens_redshift"]
            else 20000
        )

        # dVcdz_function
        dVcdz_function = self.differential_comoving_volume.function

        from .mp import lens_redshift_strongly_lensed_njit

        density_array = lens_redshift_strongly_lensed_njit(
            zs, # 1D
            zl_scaled, # 2D
            sigma_min,
            sigma_max,
            q_rvs,
            phi_rvs,
            gamma_rvs,
            shear_rvs,
            number_density,
            cross_section,
            dVcdz_function,
            integration_size,
        )

        return density_array, use_njit_sampler

    def _lens_redshift_multiprocessing(self, zl_scaled, zs):
        """
            Compute the lens redshift distribution using multiprocessing.
            
            Parameters
            ----------
            zl_scaled2d : array_like
                2D array of lens redshifts, scaled by the source redshift.
            zs : array_like
                1D array of source redshifts.
            
            Returns
            -------
            density_array : array_like
                2D array of the lens redshift distribution.
        """
   
        # 0. zs
        # 1. zl_scaled

        # 2. sigma_args; velocity_dispersion = self.velocity_dispersion.rvs and density_function = self.velocity_dispersion.function
        sigma_args = self.lens_param_samplers_params["velocity_dispersion"]
        # importance sampling in the backend requires sampling from a uniform distribution within the range
        sigma_args = [
            sigma_args["sigma_min"],
            sigma_args["sigma_max"],
        ]
        # The following is used for prior reweighting in the backend
        sigma_args += [
            # if zl dependent vd: 2d np.ndarray, shape=(n_sigma, n_zl)
            # if zl independent vd: 1d np.ndarray, shape=(n_sigma)
            np.array(self.velocity_dispersion.x_array),
            # if zl dependent vd: 2d np.ndarray, shape=(n_sigma, n_zl)
            # if zl independent vd: None
            np.array(self.velocity_dispersion.conditioned_y_array)
            if self.velocity_dispersion.conditioned_y_array is not None
            else None,
        ]
        if self.velocity_dispersion.function_spline is not None:
            # if zl dependent vd: 3d np.ndarray, shape=(n_sigma, n_zl, n_spline_coeff)
            # if zl independent vd: 2d np.ndarray, shape=(n_sigma, n_spline_coeff)
            sigma_args += [np.array(self.velocity_dispersion.function_spline)]

        # 3. q_args; axis_ratio = self.axis_ratio.rvs
        if self.lens_param_samplers["axis_ratio"] == "axis_ratio_uniform":
            q_args = [
                self.lens_param_samplers["axis_ratio"],
                self.lens_param_samplers_params["axis_ratio"]["q_min"],
                self.lens_param_samplers_params["axis_ratio"]["q_max"],
            ]
        else:
            q_args = [
                self.lens_param_samplers["axis_ratio"],
                np.array(self.axis_ratio.cdf_values),
                np.array(self.axis_ratio.x_array),
                np.array(self.axis_ratio.conditioned_y_array)
                if self.axis_ratio.conditioned_y_array is not None
                else None,
            ]

        # 4. Da_args; angular_diameter_distance = self.angular_diameter_distance.function
        Da_args = [
            np.array(self.angular_diameter_distance.function_spline),
            np.array(self.angular_diameter_distance.x_array),
        ]

        # 5. dVcdz_args; differential_comoving_volume = self.differential_comoving_volume.function
        dVcdz_args = [
            np.array(self.differential_comoving_volume.function_spline),
            np.array(self.differential_comoving_volume.x_array),
        ]

        # 7. cross section
        # cs_fn = self.cross_section_function
        cs_args = [self.lens_functions["cross_section"]]
        if self.lens_functions["cross_section"] == "cross_section_epl_shear_interpolation":
            cs_args += [self.cs_data_points_path]

        # 8. phi; axis_rotation_angle = self.axis_rotation_angle.rvs
        if (
            self.lens_param_samplers["axis_rotation_angle"]
            == "axis_rotation_angle_uniform"
        ):
            phi_args = [
                self.lens_param_samplers["axis_rotation_angle"],
                self.lens_param_samplers_params["axis_rotation_angle"]["phi_min"],
                self.lens_param_samplers_params["axis_rotation_angle"]["phi_max"],
            ]
        else:
            phi_args = [
                self.lens_param_samplers["axis_rotation_angle"],
                np.array(self.axis_rotation_angle.cdf_values),
                np.array(self.axis_rotation_angle.x_array),
            ]

        # 9. shear; external_shear = self.external_shear.rvs
        if self.lens_param_samplers["external_shear"] == "external_shear_normal":
            shear_args = [
                self.lens_param_samplers["external_shear"],
                self.lens_param_samplers_params["external_shear"]["mean"],
                self.lens_param_samplers_params["external_shear"]["std"],
            ]
        else:
            shear_args = [
                self.lens_param_samplers["external_shear"],
                np.array(self.external_shear.cdf_values),
                np.array(self.external_shear.x_array),
                np.array(self.external_shear.conditioned_y_array)
                if self.external_shear.conditioned_y_array is not None
                else None,
            ]

        # 10. slope; density_profile_slope = self.density_profile_slope.rvs
        if (
            self.lens_param_samplers["density_profile_slope"]
            == "density_profile_slope_normal"
        ):
            slope_args = [
                self.lens_param_samplers["density_profile_slope"],
                self.lens_param_samplers_params["density_profile_slope"]["mean"],
                self.lens_param_samplers_params["density_profile_slope"]["std"],
            ]

        # 11. integration_size
        integration_size = (
            self.lens_param_samplers_params["lens_redshift"]["integration_size"]
            if "integration_size" in self.lens_param_samplers_params["lens_redshift"]
            else 25000
        )

        # 6. idx; index to order the results
        idx = np.arange(len(zs))

        input_params = np.array(
            [
                (
                    zs[i],
                    zl_scaled[i],
                    sigma_args,
                    q_args,
                    Da_args,
                    dVcdz_args,
                    idx[i],
                    cs_args,
                    phi_args,
                    shear_args,
                    slope_args,
                    integration_size,
                )
                for i in range(len(zs))
            ],
            dtype=object,
        )

        print("Computing lens redshift distribution with multiprocessing...")
        from .mp import lens_redshift_strongly_lensed_mp
        density_array = np.zeros_like(zl_scaled)
        # with Pool(processes=self.npool) as pool:
        #     for result in tqdm(
        #         pool.imap_unordered(lens_redshift_strongly_lensed_mp, input_params),
        #         total=len(zs),
        #         ncols=100,
        #         disable=False,
        #     ):
        #         # print(result)
        #         (
        #             iter_i,
        #             density_,
        #         ) = result

        #         density_array[iter_i] = density_
        # use for-loop instead
        for i in tqdm(range(len(zs))):
            result = lens_redshift_strongly_lensed_mp(input_params[i])
            (
                iter_i,
                density_,
            ) = result

            density_array[iter_i] = density_

        return density_array

            
    # def _lens_redshift_multiprocessing(self, zl_scaled, zs):
    #     """
    #     Compute the lens redshift distribution using multiprocessing. This function is used when the velocity dispersion distribution depends on lens redshift.

    #     Parameters
    #     ----------
    #     zl_scaled : array_like
    #         1D array of lens redshifts, scaled by the source redshift.
    #     zs : array_like
    #         1D array of source redshifts.
    #     zl_distribution_name : str
    #         Name of the lens redshift distribution to compute.

    #     Returns
    #     -------
    #     density_array : array_like
    #         2D array of the lens redshift distribution.
    #     """

    #     # 0. zs
    #     # 1. zl_scaled

    #     # 2. sigma_args; velocity_dispersion = self.velocity_dispersion.rvs and density_function = self.velocity_dispersion.function
    #     sigma_args = self.lens_param_samplers_params["velocity_dispersion"]
    #     # importance sampling in the backend requires sampling from a uniform distribution within the range
    #     sigma_args = [
    #         sigma_args["sigma_min"],
    #         sigma_args["sigma_max"],
    #     ]
    #     # The following is used for prior reweighting in the backend
    #     sigma_args += [
    #         # if zl dependent vd: 2d np.ndarray, shape=(n_sigma, n_zl)
    #         # if zl independent vd: 1d np.ndarray, shape=(n_sigma)
    #         np.array(self.velocity_dispersion.x_array),
    #         # if zl dependent vd: 2d np.ndarray, shape=(n_sigma, n_zl)
    #         # if zl independent vd: None
    #         np.array(self.velocity_dispersion.conditioned_y_array)
    #         if self.velocity_dispersion.conditioned_y_array is not None
    #         else None,
    #     ]
    #     if self.velocity_dispersion.function_spline is not None:
    #         # if zl dependent vd: 3d np.ndarray, shape=(n_sigma, n_zl, n_spline_coeff)
    #         # if zl independent vd: 2d np.ndarray, shape=(n_sigma, n_spline_coeff)
    #         sigma_args += [np.array(self.velocity_dispersion.function_spline)]

    #     # 3. q_args; axis_ratio = self.axis_ratio.rvs
    #     if self.lens_param_samplers["axis_ratio"] == "axis_ratio_uniform":
    #         q_args = [
    #             self.lens_param_samplers["axis_ratio"],
    #             self.lens_param_samplers_params["axis_ratio"]["q_min"],
    #             self.lens_param_samplers_params["axis_ratio"]["q_max"],
    #         ]
    #     else:
    #         q_args = [
    #             self.lens_param_samplers["axis_ratio"],
    #             np.array(self.axis_ratio.cdf_values),
    #             np.array(self.axis_ratio.x_array),
    #             np.array(self.axis_ratio.conditioned_y_array)
    #             if self.axis_ratio.conditioned_y_array is not None
    #             else None,
    #         ]

    #     # 4. Da_args; angular_diameter_distance = self.angular_diameter_distance.function
    #     Da_args = [
    #         np.array(self.angular_diameter_distance.function_spline),
    #         np.array(self.angular_diameter_distance.x_array),
    #     ]

    #     # 5. dVcdz_args; differential_comoving_volume = self.differential_comoving_volume.function
    #     dVcdz_args = [
    #         np.array(self.differential_comoving_volume.function_spline),
    #         np.array(self.differential_comoving_volume.x_array),
    #     ]

    #     # 7. cross section
    #     # cs_fn = self.cross_section_function
    #     cs_args = [self.lens_functions["cross_section"]]
    #     if self.lens_functions["cross_section"] == "cross_section_epl_shear_interpolation":
    #         cs_args += [self.cs_data_points_path]

    #     # 8. phi; axis_rotation_angle = self.axis_rotation_angle.rvs
    #     if (
    #         self.lens_param_samplers["axis_rotation_angle"]
    #         == "axis_rotation_angle_uniform"
    #     ):
    #         phi_args = [
    #             self.lens_param_samplers["axis_rotation_angle"],
    #             self.lens_param_samplers_params["axis_rotation_angle"]["phi_min"],
    #             self.lens_param_samplers_params["axis_rotation_angle"]["phi_max"],
    #         ]
    #     else:
    #         phi_args = [
    #             self.lens_param_samplers["axis_rotation_angle"],
    #             np.array(self.axis_rotation_angle.cdf_values),
    #             np.array(self.axis_rotation_angle.x_array),
    #         ]

    #     # 9. shear; external_shear = self.external_shear.rvs
    #     if self.lens_param_samplers["external_shear"] == "external_shear_normal":
    #         shear_args = [
    #             self.lens_param_samplers["external_shear"],
    #             self.lens_param_samplers_params["external_shear"]["mean"],
    #             self.lens_param_samplers_params["external_shear"]["std"],
    #         ]
    #     else:
    #         shear_args = [
    #             self.lens_param_samplers["external_shear"],
    #             np.array(self.external_shear.cdf_values),
    #             np.array(self.external_shear.x_array),
    #             np.array(self.external_shear.conditioned_y_array)
    #             if self.external_shear.conditioned_y_array is not None
    #             else None,
    #         ]

    #     # 10. slope; density_profile_slope = self.density_profile_slope.rvs
    #     if (
    #         self.lens_param_samplers["density_profile_slope"]
    #         == "density_profile_slope_normal"
    #     ):
    #         slope_args = [
    #             self.lens_param_samplers["density_profile_slope"],
    #             self.lens_param_samplers_params["density_profile_slope"]["mean"],
    #             self.lens_param_samplers_params["density_profile_slope"]["std"],
    #         ]

    #     # 11. integration_size
    #     integration_size = (
    #         self.lens_param_samplers_params["lens_redshift"]["integration_size"]
    #         if "integration_size" in self.lens_param_samplers_params["lens_redshift"]
    #         else 20000
    #     )

    #     from .mp import lens_redshift_strongly_lensed_njit

    #     density_array = lens_redshift_strongly_lensed_njit(
    #         zs,
    #         zl_scaled,
    #         sigma_args,
    #         q_args,
    #         Da_args,
    #         dVcdz_args,
    #         cs_args,
    #         phi_args,
    #         shear_args,
    #         slope_args,
    #         integration_size,
    #     )

    #     return density_array


        # if isinstance(zs, (float, int)):
        #     idx = 0  # dummy index for single zs

        #     input_params = np.array(
        #         [
        #             zs,
        #             zl_scaled,
        #             sigma_args,
        #             q_args,
        #             Da_args,
        #             dVcdz_args,
        #             idx,
        #             cs_args,
        #             phi_args,
        #             shear_args,
        #             slope_args,
        #             integration_size,
        #         ],
        #         dtype=object,
        #     )

        #     mp_fn = self._helper_sl_disribution_mp
        #     iter_i, density_ = mp_fn(input_params)
        #     density_array = (
        #         density_  # density_ is already the full array for all zl_scaled values
        #     )

        # elif isinstance(zs, np.ndarray) or isinstance(zs, list):
        #     # 6. idx; index to order the results
        #     idx = np.arange(len(zs))

        #     input_params = np.array(
        #         [
        #             (
        #                 zs[i],
        #                 zl_scaled[i],
        #                 sigma_args,
        #                 q_args,
        #                 Da_args,
        #                 dVcdz_args,
        #                 idx[i],
        #                 cs_args,
        #                 phi_args,
        #                 shear_args,
        #                 slope_args,
        #                 integration_size,
        #             )
        #             for i in range(len(zs))
        #         ],
        #         dtype=object,
        #     )

        #     # print("Computing lens redshift distribution with multiprocessing...")
        #     # # with tqdm
        #     # mp_fn = self._helper_sl_disribution_mp
        #     # density_array = np.zeros_like(zl_scaled)
        #     # with Pool(processes=self.npool) as pool:
        #     #     for result in tqdm(
        #     #         pool.imap_unordered(mp_fn, input_params),
        #     #         total=len(zs),
        #     #         ncols=100,
        #     #         disable=False,
        #     #     ):
        #     #         # print(result)
        #     #         (
        #     #             iter_i,
        #     #             density_,
        #     #         ) = result

        #     #         density_array[iter_i] = density_

        #     # return density_array
        #     # using for loop instead
        #     density_array = np.zeros_like(zl_scaled)
        #     for i in range(len(zs)):
        #         density_array[i] = self._helper_sl_disribution_mp(
        #             input_params[i]
        #         )

        #     return density_array    

    def lens_redshift_intrinsic(
        self, size=1000, zs=None, get_attribute=False, **kwargs
    ):
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

        identifier_dict = {"name": "lens_redshift_intrinsic"}
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo

        name_list = ["z_min", "z_max", "cosmology", "velocity_dispersion", "resolution"]

        identifier_dict["velocity_dispersion"] = {}
        if (
            self.velocity_dispersion.info is not None
        ):  # if velocity_dispersion is not None
            for key, value in self.velocity_dispersion.info.items():
                if key not in name_list:
                    identifier_dict["velocity_dispersion"][key] = value

        identifier_dict["resolution"] = self.create_new_interpolator[
            "lens_redshift_intrinsic"
        ]["resolution"]

        _, it_exist = interpolator_json_path(
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="lens_redshift_intrinsic",
            interpolator_name=identifier_dict["name"],
        )

        create_new = self.create_new_interpolator["lens_redshift_intrinsic"][
            "create_new"
        ]

        if not it_exist or create_new:
            z_min = 0.0001 if self.z_min == 0 else self.z_min
            z_max = self.z_max
            resolution = identifier_dict["resolution"]
            zl_array = redshift_optimal_spacing(z_min, z_max, resolution)

            # create number density
            if self.velocity_dispersion.conditioned_y_array is None:
                integrand = (  # noqa: E731
                    lambda sigma, z: self.velocity_dispersion.function(
                        np.array([sigma])
                    )[0]
                    * self.differential_comoving_volume.function(np.array([z]))[0]
                )
            else:
                integrand = (  # noqa: E731
                    lambda sigma, z: self.velocity_dispersion.function(
                        np.array([sigma]), np.array([z])
                    )[0]
                    * self.differential_comoving_volume.function(np.array([z]))[0]
                )
            # fix zl and integrate over sigma
            integral = [
                quad(
                    integrand,
                    identifier_dict["velocity_dispersion"]["sigma_min"],
                    identifier_dict["velocity_dispersion"]["sigma_max"],
                    args=(zl),
                )[0]
                for zl in zl_array
            ]
            # create spline fit for lens redshift
            number_density_spline = CubicSpline(zl_array, integral)

            # zl dependent velocity dispersion distribution
            zl_scaled = []
            number_density = []
            zl_min = 0.0001
            for i, zl_ in enumerate(zl_array):
                zl_size = 48  # number of points for each zl
                buffer_ = np.linspace(zl_min, zl_ - zl_min, zl_size)
                zl_scaled.append(buffer_ / zl_)
                number_density.append(number_density_spline(buffer_))
            zl_scaled = np.array(zl_scaled)
            number_density = np.array(number_density)

        else:
            zl_array = None
            number_density = None
            zl_scaled = None

        zl_object = FunctionConditioning(
            function=number_density,
            x_array=zl_scaled,
            conditioned_y_array=zl_array,
            non_zero_function=True,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="lens_redshift_intrinsic",
            name=identifier_dict["name"],
            create_new=create_new,
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
        )

        # un-scaled lens redshift
        cdf_values = zl_object.cdf_values
        x_array = zl_object.x_array
        y_array = zl_object.conditioned_y_array
        function_spline = zl_object.function_spline
        pdf_norm_const = zl_object.pdf_norm_const

        zl_object.function = njit(
            lambda x, y: cubic_spline_interpolator2d_array(
                x / y, y, function_spline, x_array, y_array
            )
        )

        zl_object.pdf = njit(
            lambda x, y: pdf_cubic_spline_interpolator2d_array(
                x / y, y, pdf_norm_const, function_spline, x_array, y_array
            )
            / y
        )

        zl_object.rvs = njit(
            lambda size, y: inverse_transform_sampler2d(
                size, y, cdf_values, x_array, y_array
            )
            * y
        )

        return zl_object if get_attribute else zl_object(size, zs)

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
        >>> print(self.axis_rotation_angle_uniform(size=10))
        """

        identifier_dict = {"name": "axis_rotation_angle_uniform"}
        param_dict = self.available_lens_prior_list_and_its_params[
            "axis_rotation_angle"
        ]["axis_rotation_angle_uniform"]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        low = param_dict["phi_min"]
        high = param_dict["phi_max"]
        phi_rvs = njit(
            lambda size: np.random.uniform(
                low=low,
                high=high,
                size=size,
            )
        )
        if param_dict["phi_max"] == param_dict["phi_min"]:
            phi_pdf = lambda phi: 1.0 if phi == param_dict["phi_min"] else 0.0
        else:
            phi_pdf = lambda phi: 1 / (param_dict["phi_max"] - param_dict["phi_min"])

        phi_object = FunctionConditioning(
            function=None,
            x_array=None,
            identifier_dict=param_dict,
            create_rvs=phi_rvs,
            create_pdf=phi_pdf,
            callback="rvs",
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
        >>> print(self.axis_ratio_uniform(size=10))
        """

        identifier_dict = {"name": "axis_ratio_uniform"}
        param_dict = self.available_lens_prior_list_and_its_params["axis_ratio"][
            "axis_ratio_uniform"
        ]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        low = param_dict["q_min"]
        high = param_dict["q_max"]
        q_rvs = njit(
            lambda size: np.random.uniform(
                low=low,
                high=high,
                size=size,
            )
        )
        if param_dict["q_max"] == param_dict["q_min"]:
            q_pdf = lambda q: 1.0 if q == param_dict["q_min"] else 0.0
        else:
            q_pdf = lambda q: 1 / (param_dict["q_max"] - param_dict["q_min"])

        q_object = FunctionConditioning(
            function=None,
            x_array=None,
            identifier_dict=param_dict,
            create_rvs=q_rvs,
            create_pdf=q_pdf,
            callback="rvs",
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
        >>> print(self.external_shear_normal(size=10))
        """

        identifier_dict = {"name": "external_shear_normal"}
        param = self.available_lens_prior_list_and_its_params["external_shear"][
            "external_shear_normal"
        ]
        if param:
            param.update(kwargs)
        else:
            param = kwargs
        identifier_dict.update(param)

        mean = param["mean"]
        std = param["std"]
        shear_rvs = njit(
            lambda size: np.random.normal(
                loc=mean,
                scale=std,
                size=(2, size),
            )
        )
        shear_pdf = njit(
            lambda shear1, shear2: normal_pdf_2d(
                x=shear1,
                y=shear2,
                mean_x=mean,
                mean_y=mean,
                std_x=std,
                std_y=std,
            )
        )

        shear_object = FunctionConditioning(
            function=None,
            x_array=None,
            identifier_dict=identifier_dict,
            create_rvs=shear_rvs,
            create_pdf=shear_pdf,
            callback="rvs",
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
        >>> print(self.density_profile_slope_normal(size=10))
        """

        identifier_dict = {"name": "density_profile_slope_normal"}

        param = self.available_lens_prior_list_and_its_params["density_profile_slope"][
            "density_profile_slope_normal"
        ]
        if param:
            param.update(kwargs)
        else:
            param = kwargs
        identifier_dict.update(param)

        mean = param["mean"]
        std = param["std"]
        slope_rvs = njit(
            lambda size: np.random.normal(
                loc=mean,
                scale=std,
                size=size,
            )
        )
        slope_pdf = njit(
            lambda slope: normal_pdf(
                x=slope,
                mean=mean,
                std=std,
            )
        )

        slope_object = FunctionConditioning(
            function=None,
            x_array=None,
            identifier_dict=identifier_dict,
            create_rvs=slope_rvs,
            create_pdf=slope_pdf,
            callback="rvs",
        )

        if get_attribute:
            return slope_object
        else:
            return slope_rvs(size)

    def lens_redshift_sis_haris(self, size, zs, get_attribute=False, **kwargs):
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
        >>> lens.lens_redshift_sis_haris(zs=1.0)
        """

        # # Old way
        # splineDc = self.comoving_distance.function_spline  # spline coefficients for the comoving distance and redshifts
        # splineDcInv = self.comoving_distance.function_inverse_spline # spline coefficients for the redshifts and comoving distance
        # u = np.linspace(0, 1, 500)
        # cdf = (10 * u**3 - 15 * u**4 + 6 * u**5)  # See the integral of Eq. A7 of https://arxiv.org/pdf/1807.07062.pdf (cdf)
        # zs = np.array([zs]).reshape(-1)

        # New way
        identifier_dict = {"name": "lens_redshift_sis_haris"}
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        identifier_dict["velocity_dispersion"] = dict(
            sigma_min=100.0,
            sigma_max=400.0,
            alpha=2.32,
            beta=2.67,
            phistar=8.0e-3 * self.cosmo.h**3,
            sigmastar=161.0,
        )
        identifier_dict["resolution"] = self.create_new_interpolator["lens_redshift"][
            "resolution"
        ]
        param_dict = self.available_lens_prior_list_and_its_params["lens_redshift"][
            "lens_redshift_sis_haris"
        ]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        zl_resolution = identifier_dict["resolution"]
        x_array = np.linspace(0.0, 1.0, zl_resolution)
        pdf_ = lambda x: 30 * x**2 * (1 - x) ** 2  # Haris et al. 2018 (A7)

        zl_object = FunctionConditioning(
            function=pdf_,
            x_array=x_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="lens_redshift",
            name=identifier_dict["name"],
            create_new=self.create_new_interpolator["lens_redshift"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
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
            return (
                cubic_spline_interpolator(r, function_spline_zl, x_array_zl)
                / pdf_norm_const
            )

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
        >>> print(self.velocity_dispersion(size=10))
        """

        identifier_dict = {"name": "velocity_dispersion_gengamma"}
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["resolution"] = self.create_new_interpolator[
            "velocity_dispersion"
        ]["resolution"]
        param_dict = self.available_lens_prior_list_and_its_params[
            "velocity_dispersion"
        ]["velocity_dispersion_gengamma"]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        # setting up inputs for the interpolator
        sigma_array = np.linspace(
            identifier_dict["sigma_min"],
            identifier_dict["sigma_max"],
            identifier_dict["resolution"],
        )
        pdf_func_ = lambda sigma_: gengamma.pdf(
            sigma_ / identifier_dict["sigmastar"],
            a=identifier_dict["alpha"] / identifier_dict["beta"],
            c=identifier_dict["beta"],
        )  # gengamma pdf

        sigma_object = FunctionConditioning(
            function=pdf_func_,
            x_array=sigma_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="velocity_dispersion",
            name=identifier_dict["name"],
            create_new=self.create_new_interpolator["velocity_dispersion"][
                "create_new"
            ],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
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
        >>> print(self.velocity_dispersion(size=10))
        """

        identifier_dict = {"name": "velocity_dispersion_bernardi"}
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["name"] = "velocity_dispersion_bernardi"
        identifier_dict["resolution"] = self.create_new_interpolator[
            "velocity_dispersion"
        ]["resolution"]
        param_dict = self.available_lens_prior_list_and_its_params[
            "velocity_dispersion"
        ]["velocity_dispersion_bernardi"]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        # setting up inputs for the interpolator
        sigma_array = np.linspace(
            identifier_dict["sigma_min"],
            identifier_dict["sigma_max"],
            identifier_dict["resolution"],
        )
        number_density_function = lambda sigma: phi_loc_bernardi(
            sigma,
            alpha=self.lens_param_samplers_params["velocity_dispersion"]["alpha"],
            beta=self.lens_param_samplers_params["velocity_dispersion"]["beta"],
            phistar=self.lens_param_samplers_params["velocity_dispersion"]["phistar"],
            sigmastar=self.lens_param_samplers_params["velocity_dispersion"][
                "sigmastar"
            ],
        )

        sigma_object = FunctionConditioning(
            function=number_density_function,
            x_array=sigma_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="velocity_dispersion",
            name=identifier_dict["name"],
            create_new=self.create_new_interpolator["velocity_dispersion"][
                "create_new"
            ],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
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
        >>> print(self.velocity_dispersion(size=10, zl=0.5))
        """

        identifier_dict = {"name": "velocity_dispersion_ewoud"}
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["name"] = "velocity_dispersion_ewoud"
        identifier_dict["resolution"] = self.create_new_interpolator[
            "velocity_dispersion"
        ]["resolution"]
        identifier_dict["zl_resolution"] = self.create_new_interpolator[
            "velocity_dispersion"
        ]["zl_resolution"]
        param_dict = self.available_lens_prior_list_and_its_params[
            "velocity_dispersion"
        ]["velocity_dispersion_ewoud"].copy()
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        # setting up inputs for the interpolator
        sigma_array = np.linspace(
            identifier_dict["sigma_min"],
            identifier_dict["sigma_max"],
            identifier_dict["resolution"],
        )
        number_density_function = lambda sigma, zl: phi(  # noqa: E731
            sigma,
            zl,
            alpha=self.lens_param_samplers_params["velocity_dispersion"]["alpha"],
            beta=self.lens_param_samplers_params["velocity_dispersion"]["beta"],
            phistar=self.lens_param_samplers_params["velocity_dispersion"]["phistar"],
            sigmastar=self.lens_param_samplers_params["velocity_dispersion"][
                "sigmastar"
            ],
        )

        z_min = self.z_min + 0.001 if self.z_min == 0.0 else self.z_min
        z_max = self.z_max
        z_resolution = identifier_dict["zl_resolution"]
        zl_array = redshift_optimal_spacing(z_min, z_max, z_resolution)

        sigma_object = FunctionConditioning(
            function=number_density_function,
            x_array=sigma_array,
            conditioned_y_array=zl_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="velocity_dispersion",
            name=identifier_dict["name"],
            create_new=self.create_new_interpolator["velocity_dispersion"][
                "create_new"
            ],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
        )

        return sigma_object if get_attribute else sigma_object.rvs(size, zl)

    ##################
    # Lens functions #
    ##################
    def optical_depth_numerical(self, zs, get_attribute=False, **kwargs):
        identifier_dict = {"name": "optical_depth_numerical"}
        identifier_dict["resolution"] = self.create_new_interpolator["optical_depth"][
            "resolution"
        ]
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["velocity_dispersion"] = self.lens_param_samplers_params[
            "velocity_dispersion"
        ]
        if identifier_dict["velocity_dispersion"]:  # if velocity_dispersion is not None
            identifier_dict["velocity_dispersion"]["name"] = str(
                self.lens_param_samplers["velocity_dispersion"]
            )
        identifier_dict["axis_ratio"] = self.lens_param_samplers_params["axis_ratio"]
        if identifier_dict["axis_ratio"]:  # if axis_ratio is not None
            identifier_dict["axis_ratio"]["name"] = str(
                self.lens_param_samplers["axis_ratio"]
            )
        identifier_dict["axis_rotation_angle"] = self.lens_param_samplers_params[
            "axis_rotation_angle"
        ]
        if identifier_dict["axis_rotation_angle"]:  # if axis_rotation_angle is not None
            identifier_dict["axis_rotation_angle"]["name"] = str(
                self.lens_param_samplers["axis_rotation_angle"]
            )
        identifier_dict["density_profile_slope"] = self.lens_param_samplers_params[
            "density_profile_slope"
        ]
        if identifier_dict[
            "density_profile_slope"
        ]:  # if density_profile_slope is not None
            identifier_dict["density_profile_slope"]["name"] = str(
                self.lens_param_samplers["density_profile_slope"]
            )
        identifier_dict["external_shear"] = self.lens_param_samplers_params[
            "external_shear"
        ]
        if identifier_dict["external_shear"]:  # if external_shear is not None
            identifier_dict["external_shear"]["name"] = str(
                self.lens_param_samplers["external_shear"]
            )

        param_dict = self.available_lens_functions_and_its_params["optical_depth"][
            "optical_depth_numerical"
        ]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        z_min = self.z_min if self.z_min>0. else 0.0001
        z_max = self.z_max
        resolution = identifier_dict['resolution']
        zs_array = np.geomspace(
            z_min,
            z_max,
            resolution
        )
        # z_min = self.z_min + 0.001 if self.z_min == 0.0 else self.z_min
        # z_max = self.z_max
        # z_resolution = identifier_dict["resolution"]
        # zs_array = redshift_optimal_spacing(z_min, z_max, z_resolution)

        def tau(zs):
            # self.lens_redshift.function gives cross-section
            integrand = lambda zl_, zs_: self.lens_redshift.function(np.array([zl_]), np.array([zs_]))[0] # noqa: E731
            integral = [quad(integrand, 0.0, z, args=(z))[0] for z in zs]
            return integral

        tau_object = FunctionConditioning(
            function=tau,
            x_array=zs_array,
            conditioned_y_array=None,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="optical_depth",
            name=identifier_dict["name"],
            create_new=self.create_new_interpolator["optical_depth"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="function",
        )

        return tau_object if get_attribute else tau_object.function(zs)

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
        >>> print(self.optical_depth_sis_haris(zs=1.0))
        """

        identifier_dict = {"name": "optical_depth_sis_haris"}
        identifier_dict["z_min"] = self.z_min if self.z_min > 0.0 else 0.0001
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["resolution"] = self.create_new_interpolator["optical_depth"][
            "resolution"
        ]
        param_dict = self.available_lens_functions_and_its_params["optical_depth"][
            "optical_depth_sis_haris"
        ]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        z_min = identifier_dict["z_min"]
        z_max = identifier_dict["z_max"]
        z_resolution = identifier_dict["resolution"]
        zs_arr = redshift_optimal_spacing(z_min, z_max, z_resolution)

        def tau(zs):
            Dc = self.comoving_distance.function(zs)
            return 4.192e-15 * Dc**3

        tau_object = FunctionConditioning(
            function=tau,
            x_array=zs_arr,
            conditioned_y_array=None,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="optical_depth",
            name=identifier_dict["name"],
            create_new=self.create_new_interpolator["optical_depth"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="function",
        )

        return tau_object if get_attribute else tau_object.function(zs)

    ###########################
    # cross section functions #
    ###########################
    def cross_section_sis(self, zs=None, zl=None, sigma=None, get_attribute=False, **kwargs):
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
        >>> print(self.cross_section_sis(sigma=200., zl=0.5, zs=1.0))
        """

        # theta_E = self.compute_einstein_radii(sigma, zl, zs)
        # Compute the angular diameter distances
        _Da = self.angular_diameter_distance.function
        _Da_z1z2 = self.angular_diameter_distance_z1z2.function

        # Compute the Einstein radii
        if get_attribute:
            @njit
            def _cross_section_sis(zs, zl, sigma):
                Ds = _Da(zs)
                Dls = _Da_z1z2(zl, zs)
                theta_E = 4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / (Ds)
                return np.pi * theta_E**2
            return _cross_section_sis
        else:
            Ds = _Da(zs)
            Dls = _Da_z1z2(zl, zs)
            theta_E = 4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / (Ds)
            return np.pi * theta_E**2   

    def cross_section_sie_feixu(self, zs=None, zl=None, sigma=None, q=None, get_attribute=False, **kwargs):
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
        >>> print(self.cross_section_sie_feixu(sigma=200., zl=0.5, zs=1.0, q=1.0))
        """

        if get_attribute:
            _cross_section_sis = self.cross_section_sis(sigma=None, zl=None, zs=None, get_attribute=True)
            @njit
            def _cross_section_sie_feixu(zs, zl, sigma, q):
                return phi_cut_SIE(q) * _cross_section_sis(zs, zl, sigma)
            return _cross_section_sie_feixu
        else:
            _cross_section_sis = self.cross_section_sis(zs=zs, zl=zl, sigma=sigma, get_attribute=False)
            return phi_cut_SIE(q) * _cross_section_sis

    def cross_section_epl_shear_numerical(
        self,
        theta_E,
        gamma,
        gamma1,
        gamma2,
        q=None,
        phi=None,
        e1=None,
        e2=None,
        verbose=False,
        **kwargs,
    ):
        """
        Function to compute the strong lensing cross-section numerically for EPL + external shear lenses.

        Parameters
        ----------
        theta_E : `numpy.ndarray`
            Einstein radii of the lens galaxies in radians
        gamma : `numpy.ndarray`
            external shear magnitudes of the lens galaxies
        gamma1 : `numpy.ndarray`
            external shear component 1 of the lens galaxies
        gamma2 : `numpy.ndarray`
            external shear component 2 of the lens galaxies
        q : `numpy.ndarray`
            axis ratios of the lens galaxies
        phi : `numpy.ndarray`
            axis rotation angles of the lens galaxies in radians
        e1 : `numpy.ndarray`
            ellipticity component 1 of the lens galaxies
        e2 : `numpy.ndarray`
            ellipticity component 2 of the lens galaxies

        Returns
        -------
        cross_section : `numpy.ndarray`
            strong lensing cross-section of the lens galaxies in square radians

        """

        # Transform the axis ratio and the axis rotation angle, to ellipticities e1, e2, using lenstronomy
        if e1 is None or e2 is None:
            if q is None or phi is None:
                raise ValueError("Either (e1, e2) or (q, phi) must be provided.")
            e1, e2 = phi_q2_ellipticity_hemanta(phi, q)

        size = len(theta_E)
        iter_ = np.arange(size)
        points = np.column_stack((theta_E, e1, e2, gamma, gamma1, gamma2, iter_))

        area_list = np.ones(size)
        with Pool(processes=self.npool) as pool:
            if verbose:
                for result in tqdm(
                    pool.imap_unordered(cross_section_mp, points),
                    total=size,
                    ncols=100,
                    disable=False,
                ):
                    (
                        iter_i,
                        area,
                    ) = result
                    area_list[int(iter_i)] = area
            else:
                for result in pool.map(cross_section_mp, points):
                    (
                        iter_i,
                        area,
                    ) = result
                    area_list[int(iter_i)] = area

        return np.array(area_list)

    def create_parameter_grid(self, size_list=[25, 25, 45, 15, 15]):
        size = 1000000
        zl = np.random.uniform(0.1, 10.0, size)
        if self.velocity_dispersion.conditioned_y_array is not None:
            sigma = self.velocity_dispersion.rvs(size, zl)
        else:
            sigma = self.velocity_dispersion.rvs(size)

        if self.axis_ratio.conditioned_y_array is not None:
            q = self.axis_ratio.rvs(size, sigma)
        else:
            q = self.axis_ratio.rvs(size)
        q_min, q_max = q.min(), q.max()
        del_q = 2 * (q_max - q_min) / size_list[0]
        if (q_min - del_q) > 0:
            q_min = q_min - del_q
        if (q_max + del_q) < 1:
            q_max = q_max + del_q
        if del_q == 0:
            q_min = max(0.2, q_min - 0.01)
            q_max = min(1.0, q_max + 0.01)

        phi = self.axis_rotation_angle.rvs(size)
        phi_min, phi_max = phi.min(), phi.max()
        del_phi = 2 * (phi_max - phi_min) / size_list[1]
        if del_phi == 0:
            phi_min = max(0.0, phi_min - 0.01)
            phi_max = min(2.0 * np.pi, phi_max + 0.01)

        e1, e2 = phi_q2_ellipticity_hemanta(phi, q)

        gamma = self.density_profile_slope.rvs(size)
        gamma_min, gamma_max = gamma.min(), gamma.max()
        del_gamma = 2 * (gamma_max - gamma_min) / size_list[2]
        gamma_min = gamma_min - del_gamma
        gamma_max = gamma_max + del_gamma
        if del_gamma == 0:
            gamma_min = gamma_min - 0.01
            gamma_max = gamma_max + 0.01

        gamma1, gamma2 = self.external_shear.rvs(size)
        gamma1_min, gamma1_max = gamma1.min(), gamma1.max()
        gamma2_min, gamma2_max = gamma2.min(), gamma2.max()
        del_gamma1 = 2 * (gamma1_max - gamma1_min) / size_list[3]
        del_gamma2 = 2 * (gamma2_max - gamma2_min) / size_list[4]
        gamma1_min = gamma1_min - del_gamma1
        gamma1_max = gamma1_max + del_gamma1
        gamma2_min = gamma2_min - del_gamma2
        gamma2_max = gamma2_max + del_gamma2
        if del_gamma1 == 0:
            gamma1_min = gamma1_min - 0.0001
            gamma1_max = gamma1_max + 0.0001
        if del_gamma2 == 0:
            gamma2_min = gamma2_min - 0.0001
            gamma2_max = gamma2_max + 0.0001

        ###########################
        # sampling the parameters #
        ###########################
        q = np.linspace(q_min, q_max, size_list[0])
        phi = np.linspace(phi_min, phi_max, size_list[1])
        e1, e2 = phi_q2_ellipticity_hemanta(phi, q)
        e1 = np.sort(e1)
        e2 = np.sort(e2)

        gamma = np.linspace(gamma_min, gamma_max, size_list[2])

        gamma1 = np.linspace(gamma1_min, gamma1_max, size_list[3])
        gamma2 = np.linspace(gamma2_min, gamma2_max, size_list[4])

        return q, phi, e1, e2, gamma, gamma1, gamma2

    def cross_section_epl_shear_interpolation_init(self, file_path, size_list):
        print(f"Cross section interpolation data points will be created at {file_path}")

        #----------------------------------------------------
        # cross section unit to cross section ratio fitting 
        #----------------------------------------------------
        # if self.lens_type == "epl_shear_galaxy":
        size = 10000
        zs = np.random.uniform(self.z_min, self.z_max, size)
        zl = np.zeros(size)
        for i in range(size):
            zl[i] = np.random.uniform(0.001, zs[i] - 0.001)

        sigma = np.random.uniform(
            self.lens_param_samplers_params["velocity_dispersion"]["sigma_min"],
            self.lens_param_samplers_params["velocity_dispersion"]["sigma_max"],
            size,
        )
        theta_E = self.compute_einstein_radii(sigma, zl, zs)
        theta_E_unit = np.ones_like(theta_E)

        # Sample individual parameter values (not grid)
        try:
            q = self.axis_ratio.rvs(size, sigma)
        except TypeError:
            q = self.axis_ratio.rvs(size)
        phi = self.axis_rotation_angle.rvs(size)
        e1, e2 = phi_q2_ellipticity_hemanta(phi, q)
        gamma = self.density_profile_slope.rvs(size)
        gamma1, gamma2 = self.external_shear.rvs(size)

        cs = self.cross_section_epl_shear_numerical(
            theta_E=theta_E, e1=e1, e2=e2, gamma=gamma, gamma1=gamma1, gamma2=gamma2
        )
        cs_unit = self.cross_section_epl_shear_numerical(
            theta_E=theta_E_unit,
            e1=e1,
            e2=e2,
            gamma=gamma,
            gamma1=gamma1,
            gamma2=gamma2,
        )

        x = np.pi * theta_E**2

        idx_cs = (cs_unit > 0.0) & (cs > 0.0) & (x > 0.0)
        ratio = cs[idx_cs] / cs_unit[idx_cs]
        x = x[idx_cs]

        # linear fitting
        from scipy.stats import linregress

        cs_csunit_slope, cs_csunit_intercept, _, _, _ = linregress(x, ratio)

        # cs data points for interpolation
        size = np.prod(size_list)
        _, _, e1, e2, gamma, gamma1, gamma2 = self.create_parameter_grid(
            size_list=size_list
        )

        e1_arr, e2_arr, gamma_arr, gamma1_arr, gamma2_arr = np.meshgrid(
            e1, e2, gamma, gamma1, gamma2, indexing="ij"
        )
        e1_arr = e1_arr.flatten()
        e2_arr = e2_arr.flatten()
        gamma_arr = gamma_arr.flatten()
        gamma1_arr = gamma1_arr.flatten()
        gamma2_arr = gamma2_arr.flatten()

        cs = self.cross_section_epl_shear_numerical(
            theta_E=np.ones(size),
            gamma=gamma_arr,
            gamma1=gamma1_arr,
            gamma2=gamma2_arr,
            e1=e1_arr,
            e2=e2_arr,
            verbose=True,
        )

        cs = cs.reshape(size_list)

        # save cs data
        save_json(
            file_path,
            [
            e1,
            e2,
            gamma,
            gamma1,
            gamma2,
            cs,
            cs_csunit_slope,
            cs_csunit_intercept,
            ],
        )

        return e1, e2, gamma, gamma1, gamma2, cs, cs_csunit_slope, cs_csunit_intercept

    def cross_section_epl_shear_interpolation(
        self,
        zs,
        zl,
        sigma,
        q,
        phi,
        gamma,
        gamma1,
        gamma2,
        get_attribute=False,
        size_list=[25, 25, 45, 15, 15],
        **kwargs,
    ):
        """
        Function to compute the cross-section correction factor
        """

        from .cross_section_interpolator import make_cross_section_reinit

        param_dict_given = dict(
            z_min=self.z_min,
            z_max=self.z_max,
            cosmology=self.cosmo,
            resolution=self.create_new_interpolator["cross_section"]["resolution"],
            size_list=size_list,
            velocity_dispersion=self.lens_param_samplers_params["velocity_dispersion"],
            axis_ratio=self.lens_param_samplers_params["axis_ratio"],
            axis_rotation_angle=self.lens_param_samplers_params["axis_rotation_angle"],
            density_profile_slope=self.lens_param_samplers_params[
                "density_profile_slope"
            ],
            external_shear=self.lens_param_samplers_params["external_shear"],
        )

        file_path, it_exist = interpolator_json_path(
            identifier_dict=param_dict_given,
            directory=self.directory,
            sub_directory="cross_section_function",
            interpolator_name="cross_section_function",
        )

        if (self.create_new_interpolator["cross_section"]["create_new"] is True) or (
            it_exist is False
        ):
            (
                e1_grid,
                e2_grid,
                gamma_grid,
                gamma1_grid,
                gamma2_grid,
                cs_spline_coeff_grid,
                _, # cs_csunit_slope
                _, # cs_csunit_intercept
            ) = self.cross_section_epl_shear_interpolation_init(file_path, size_list)
        else:

            print(f"Cross section interpolation data points loaded from {file_path}")
            (
                e1_grid,
                e2_grid,
                gamma_grid,
                gamma1_grid,
                gamma2_grid,
                cs_spline_coeff_grid,
                _, # cs_csunit_slope
                _, # cs_csunit_intercept
            ) = load_json(file_path)

        self.cs_data_points_path = file_path

        # Create the interpolator instance
        cs_caculator = make_cross_section_reinit(
            e1_grid=np.array(e1_grid),
            e2_grid=np.array(e2_grid),
            gamma_grid=np.array(gamma_grid),
            gamma1_grid=np.array(gamma1_grid),
            gamma2_grid=np.array(gamma2_grid),
            cs_spline_coeff_grid=np.array(cs_spline_coeff_grid),
            Da_instance=self.angular_diameter_distance.function,
            csunit_to_cs_slope=0.31830988618379075,
            csunit_to_cs_intercept=-3.2311742677852644e-27,
        )
        # # define conversion function: cs_unit to cs
        # self.cs_unit_to_cs = lambda area_sis: cs_csunit_slope * area_sis + cs_csunit_intercept

        return (
            cs_caculator
            if get_attribute
            else cs_caculator(
                zs=zs,
                zl=zl,
                sigma=sigma,
                q=q,
                phi=phi,
                gamma=gamma,
                gamma1=gamma1,
                gamma2=gamma2,
            )
        )

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
        >>> print(self.optical_depth(np.array([0.1,0.2,0.3])))
        """

        return self._optical_depth

    @optical_depth.setter
    def optical_depth(self, input_function):
        if (
            input_function
            in self.available_lens_functions_and_its_params["optical_depth"]
        ):
            print(f"using ler available optical depth function : {input_function}")
            args = self.lens_functions_params["optical_depth"]
            if args is None:
                self._optical_depth = getattr(self, input_function)(
                    zs=None, get_attribute=True
                )
            else:
                self._optical_depth = getattr(self, input_function)(
                    zs=None, get_attribute=True, **args
                )
        elif callable(input_function):
            print(f"using user provided custom optical depth function")
            self._optical_depth = FunctionConditioning(
                function=None, x_array=None, create_function=input_function
            )
        elif isinstance(input_function, object):
            print(f"using user provided custom optical depth class/object")
            self._optical_depth = input_function
        else:
            raise ValueError(
                "input_function not in available_lens_functions_and_its_params['optical_depth']"
            )

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
        >>> print(self.velocity_dispersion(size=10))
        """

        return self._velocity_dispersion

    @velocity_dispersion.setter
    def velocity_dispersion(self, prior):
        if (
            prior
            in self.available_lens_prior_list_and_its_params["velocity_dispersion"]
        ):
            print(f"using ler available velocity dispersion function : {prior}")
            args = self.lens_param_samplers_params["velocity_dispersion"]
            if prior == "velocity_dispersion_choi":
                prior = "velocity_dispersion_bernardi"
            if args is None:
                self._velocity_dispersion = getattr(self, prior)(
                    size=None, zl=None, get_attribute=True
                )
            else:
                self._velocity_dispersion = getattr(self, prior)(
                    size=None, zl=None, get_attribute=True, **args
                )

        elif isinstance(prior, object):
            print("using user provided custom velocity_dispersion class/object")
            self._velocity_dispersion = prior

        elif callable(prior):
            raise ValueError(
                "custom fuction for velocity dispersion should be provided as a class/object. See ler.utils.FunctionConditioning for more details."
            )

        else:
            raise ValueError(
                f"velocity_dispersion prior not in available_lens_prior_list_and_its_params['velocity_dispersion']"
            )

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
        >>> print(self.axis_ratio(sigma=200.))
        """

        return self._axis_ratio

    @axis_ratio.setter
    def axis_ratio(self, prior):
        if prior in self.available_lens_prior_list_and_its_params["axis_ratio"]:
            print(f"using ler available axis_ratio function : {prior}")
            args = self.lens_param_samplers_params["axis_ratio"]
            if args is None:
                self._axis_ratio = getattr(self, prior)(
                    size=None, sigma=None, get_attribute=True
                )
            else:
                self._axis_ratio = getattr(self, prior)(
                    size=None, sigma=None, get_attribute=True, **args
                )
        elif callable(prior):
            print("using user provided custom axis_ratio function")
            self._axis_ratio = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior
            )
        elif isinstance(prior, object):
            print("using user provided custom axis_ratio class/object")
            self._axis_ratio = prior
        else:
            raise ValueError(
                "prior not in available_lens_prior_list_and_its_params['axis_ratio']"
            )

    @property
    def lens_redshift(self):
        return self._lens_redshift

    @lens_redshift.setter
    def lens_redshift(self, prior):
        if prior in self.available_lens_prior_list_and_its_params["lens_redshift"]:
            print(f"using ler available lens_redshift function : {prior}")
            args = self.lens_param_samplers_params["lens_redshift"]
            if args is None:
                self._lens_redshift = getattr(self, prior)(
                    size=None, zs=None, get_attribute=True
                )
            else:
                self._lens_redshift = getattr(self, prior)(
                    size=None, zs=None, get_attribute=True, **args
                )
        elif callable(prior):
            print("using user provided custom lens_redshift function")
            self._lens_redshift = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior
            )
        elif isinstance(prior, object):
            print("using user provided custom lens_redshift class/object")
            self._lens_redshift = prior
        else:
            raise ValueError(
                "prior not in available_lens_prior_list_and_its_params['lens_redshift']"
            )

    @property
    def axis_rotation_angle(self):
        return self._axis_rotation_angle

    @axis_rotation_angle.setter
    def axis_rotation_angle(self, prior):
        if (
            prior
            in self.available_lens_prior_list_and_its_params["axis_rotation_angle"]
        ):
            print(f"using ler available axis_rotation_angle function : {prior}")
            args = self.lens_param_samplers_params["axis_rotation_angle"]
            if args is None:
                self._axis_rotation_angle = getattr(self, prior)(
                    size=None, get_attribute=True
                )
            else:
                self._axis_rotation_angle = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif callable(prior):
            print("using user provided custom axis_rotation_angle function")
            self._axis_rotation_angle = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior
            )
        elif isinstance(prior, object):
            print("using user provided custom axis_rotation_angle class/object")
            self._axis_rotation_angle = prior
        else:
            raise ValueError(
                "prior not in available_lens_prior_list_and_its_params['axis_rotation_angle']"
            )

    @property
    def external_shear(self):
        return self._external_shear

    @external_shear.setter
    def external_shear(self, prior):
        if prior in self.available_lens_prior_list_and_its_params["external_shear"]:
            print(f"using ler available external_shear function : {prior}")
            args = self.lens_param_samplers_params["external_shear"]
            if args is None:
                self._external_shear = getattr(self, prior)(
                    size=None, get_attribute=True
                )
            else:
                self._external_shear = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif callable(prior):
            print("using user provided custom external_shear function")
            self._external_shear = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior
            )
        elif isinstance(prior, object):
            print("using user provided custom external_shear class/object")
            self._external_shear = prior
        else:
            raise ValueError(
                "prior not in available_lens_prior_list_and_its_params['external_shear']"
            )

    @property
    def external_shear_sl(self):
        return self._external_shear_sl

    @external_shear_sl.setter
    def external_shear_sl(self, prior):
        if prior in self.available_lens_prior_list_and_its_params["external_shear_sl"]:
            print(f"using ler available external_shear_sl function : {prior}")
            args = self.lens_param_samplers_params["external_shear_sl"]
            if args is None:
                self._external_shear_sl = getattr(self, prior)(
                    size=None, get_attribute=True
                )
            else:
                self._external_shear_sl = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif callable(prior):
            print("using user provided custom external_shear_sl function")
            self._external_shear_sl = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior
            )
        elif isinstance(prior, object):
            print("using user provided custom external_shear_sl class/object")
            self._external_shear_sl = prior
        else:
            raise ValueError(
                "prior not in available_lens_prior_list_and_its_params['external_shear_sl']"
            )

    @property
    def density_profile_slope(self):
        return self._density_profile_slope

    @density_profile_slope.setter
    def density_profile_slope(self, prior):
        if (
            prior
            in self.available_lens_prior_list_and_its_params["density_profile_slope"]
        ):
            print(f"using ler available density_profile_slope function : {prior}")
            args = self.lens_param_samplers_params["density_profile_slope"]
            if args is None:
                self._density_profile_slope = getattr(self, prior)(
                    size=None, get_attribute=True
                )
            else:
                self._density_profile_slope = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif callable(prior):
            print("using user provided custom density_profile_slope function")
            self._density_profile_slope = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior
            )
        elif isinstance(prior, object):
            print("using user provided custom density_profile_slope class/object")
            self._density_profile_slope = prior
        else:
            raise ValueError(
                "prior not in available_lens_prior_list_and_its_params['density_profile_slope']"
            )

    @property
    def density_profile_slope_sl(self):
        return self._density_profile_slope_sl

    @density_profile_slope_sl.setter
    def density_profile_slope_sl(self, prior):
        if (
            prior
            in self.available_lens_prior_list_and_its_params["density_profile_slope_sl"]
        ):
            print(f"using ler available density_profile_slope_sl function : {prior}")
            args = self.lens_param_samplers_params["density_profile_slope_sl"]
            if args is None:
                self._density_profile_slope_sl = getattr(self, prior)(
                    size=None, get_attribute=True
                )
            else:
                self._density_profile_slope_sl = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif callable(prior):
            print("using user provided custom density_profile_slope_sl function")
            self._density_profile_slope_sl = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior
            )
        elif isinstance(prior, object):
            print("using user provided custom density_profile_slope_sl class/object")
            self._density_profile_slope_sl = prior
        else:
            raise ValueError(
                "prior not in available_lens_prior_list_and_its_params['density_profile_slope_sl']"
            )

    @property
    def cross_section(self):
        '''
        Lensing cross section for individual lensing events

        Parameters
        ----------
        

        Returns
        -------
        cross_section : float
            Lensing cross section for individual lensing events
        '''
        return self._cross_section

    @cross_section.setter
    def cross_section(self, cross_section):
        if cross_section == "cross_section_epl_shear_interpolation":
            # this will initialize the cross section interpolator
            size_list = self.create_new_interpolator["cross_section"]["resolution"]
            self._cross_section = self.cross_section_epl_shear_interpolation(
                zs=None,
                zl=None,
                sigma=None,
                q=None,
                phi=None,
                gamma=None,
                gamma1=None,
                gamma2=None,
                get_attribute=True,
                size_list=size_list,
            )
        elif cross_section in ["cross_section_sie_feixu", "cross_section_sis"]:
            self._cross_section = getattr(self, cross_section)(get_attribute=True)
        elif callable(cross_section):
            self._cross_section = cross_section
        elif isinstance(cross_section, str) and hasattr(self, cross_section):
            self._cross_section = getattr(self, cross_section)
        else:
            raise ValueError(
                "cross_section should be a string, callable, or an instance of a class"
            )

    ####################################
    # Available samplers and functions #
    ####################################
    @property
    def available_lens_prior_list_and_its_params(self):
        """
        Dictionary with list all the available priors and it's corresponding parameters. This is an immutable instance attribute.
        """

        self._available_lens_prior_list_and_its_params = dict(
            source_redshift_sl=dict(strongly_lensed_source_redshifts=None),
            lens_redshift=dict(
                lens_redshift_sis_haris=None,
                lens_redshift_strongly_lensed_numerical=dict(integration_size=25000),
                lens_redshift_strongly_lensed_hemanta=None,
            ),
            velocity_dispersion=dict(
                velocity_dispersion_gengamma=dict(
                    sigma_min=100.0,
                    sigma_max=400.0,
                    alpha=0.94,
                    beta=1.85,
                    phistar=2.099e-2 * (self.cosmo.h / 0.7) ** 3,
                    sigmastar=113.78,
                ),
                velocity_dispersion_choi=dict(
                    sigma_min=100.0,
                    sigma_max=400.0,
                    alpha=2.32,
                    beta=2.67,
                    phistar=8.0e-3 * self.cosmo.h**3,
                    sigmastar=161.0,
                ),
                velocity_dispersion_bernardi=dict(
                    sigma_min=100.0,
                    sigma_max=400.0,
                    alpha=0.94,
                    beta=1.85,
                    phistar=2.099e-2 * (self.cosmo.h / 0.7) ** 3,
                    sigmastar=113.78,
                ),
                velocity_dispersion_ewoud=dict(
                    sigma_min=100.0,
                    sigma_max=400.0,
                    alpha=0.94,
                    beta=1.85,
                    phistar=2.099e-2 * (self.cosmo.h / 0.7) ** 3,
                    sigmastar=113.78,
                ),
            ),
            axis_ratio=dict(
                axis_ratio_rayleigh=dict(q_min=0.2, q_max=1.0),
                axis_ratio_padilla_strauss=dict(q_min=0.2, q_max=1.0),
                axis_ratio_uniform=dict(q_min=0.2, q_max=1.0),
            ),
            axis_rotation_angle=dict(
                axis_rotation_angle_uniform=dict(phi_min=0.0, phi_max=2 * np.pi),
            ),
            external_shear=dict(
                external_shear_normal=dict(mean=0.0, std=0.05),
            ),
            external_shear_sl=dict(
                external_shear_normal=dict(
                    mean=0.0, std=0.05
                ),
                external_shear_sl_numerical_hemanta=dict(
                    external_shear_normal=dict(mean=0.0, std=0.05)
                ),
            ),
            density_profile_slope=dict(
                density_profile_slope_normal=dict(mean=1.99, std=0.149),
            ),
            density_profile_slope_sl=dict(
                density_profile_slope_normal=dict(
                    mean=2.091, std=0.133
                ),
                density_profile_slope_sl_numerical_hemanta=dict(
                    density_profile_slope_normal=dict(mean=1.99, std=0.149)
                ),
            ),
            source_parameters=dict(sample_gw_parameters=None),
        )

        return self._available_lens_prior_list_and_its_params

    @property
    def available_lens_functions_and_its_params(self):
        """
        Dictionary with list all the available lens functions. This is an immutable instance attribute.
        """

        self._available_lens_functions_and_its_params = dict(
            cross_section_based_sampler=dict(
                # rejection_sampling_with_cross_section_sie_feixu=None,
                # rejection_sampling_with_cross_section_sis=None,
                rejection_sampling_with_cross_section=dict(saftey_factor=1.2),
                importance_sampling_with_cross_section=dict(n_prop=200),
            ),
            optical_depth=dict(
                optical_depth_sis_haris=None,
                optical_depth_epl_shear_hemanta=None,
                optical_depth_numerical=None,
            ),
            param_sampler_type=dict(
                sample_all_routine_epl_shear_sl=None,
            ),
            cross_section=dict(
                cross_section_sie_feixu=None,
                cross_section_sis=None,
                cross_section_epl_shear_numerical=None,
                cross_section_epl_shear_interpolation=None,
            ),
        )

        return self._available_lens_functions_and_its_params
