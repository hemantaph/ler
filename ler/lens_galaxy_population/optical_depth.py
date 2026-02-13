"""
Module for optical depth and lens parameter distribution calculations.

This module provides the ``OpticalDepth`` class for computing strong lensing
optical depth, cross-section, sampling velocity dispersion, axis ratio, and other lens galaxy
population parameters. It supports multiple lens models including SIS, SIE,
and EPL + external shear.

Key features: \n
- Optical depth computation for strong gravitational lensing \n
- Velocity dispersion sampling with multiple models \n
- Lens redshift distribution sampling \n
- Cross-section calculations for various lens models \n

Copyright (C) 2024 Hemantakumar Phurailatpam. Distributed under MIT License.
"""

from numba import njit
from multiprocessing import Pool
import numpy as np
from scipy.integrate import quad
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
    redshift_optimal_spacing,
)

from .lens_functions import (
    phi_cut_SIE,
    phi_q2_ellipticity,
    cross_section,
)

from .mp import cross_section_mp


class OpticalDepth:
    """
    Class for computing optical depth and lens galaxy population parameters.

    This class calculates strong lensing optical depth, velocity dispersion,
    axis ratio, and other parameters for a lens galaxy population. It supports
    SIS, SIE, and EPL + external shear lens models with customizable samplers
    and interpolators for efficient computation.

    Key Features: \n
    - Multiple lens model support (SIS, SIE, EPL + shear) \n
    - Configurable velocity dispersion distributions \n
    - Cached interpolators for fast optical depth computation \n
    - Flexible parameter sampling with user-defined priors \n

    Parameters
    ----------
    npool : ``int``
        Number of processors for multiprocessing. \n
        default: 4
    z_min : ``float``
        Minimum redshift of the lens galaxy population. \n
        default: 0.0
    z_max : ``float``
        Maximum redshift of the lens galaxy population. \n
        default: 10.0
    cosmology : ``astropy.cosmology`` or ``None``
        Cosmology object for distance calculations. \n
        default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
    lens_type : ``str``
        Type of lens galaxy model. \n
        Options: \n
        - 'epl_shear_galaxy': Elliptical power-law with external shear \n
        - 'sie_galaxy': Singular isothermal ellipsoid \n
        - 'sis_galaxy': Singular isothermal sphere \n
        default: 'epl_shear_galaxy'
    lens_functions : ``dict`` or ``None``
        Dictionary with lens-related functions. \n
        default: None (uses defaults for lens_type)
    lens_functions_params : ``dict`` or ``None``
        Dictionary with parameters for lens-related functions. \n
        default: None
    lens_param_samplers : ``dict`` or ``None``
        Dictionary of sampler functions for lens parameters. \n
        default: None (uses defaults for lens_type)
    lens_param_samplers_params : ``dict`` or ``None``
        Dictionary with parameters for the samplers. \n
        default: None
    directory : ``str``
        Directory where interpolators are saved. \n
        default: './interpolator_json'
    create_new_interpolator : ``bool`` or ``dict``
        Whether to create new interpolators. \n
        default: False
    verbose : ``bool``
        If True, prints additional information. \n
        default: False

    Examples
    --------
    Basic usage:

    >>> from ler.lens_galaxy_population import OpticalDepth
    >>> od = OpticalDepth()
    >>> tau = od.optical_depth(zs=np.array([1.0, 2.0]))


    Instance Methods
    ----------
    OpticalDepth has the following instance methods: \n
    +-----------------------------------------------------+----------------------------------------------------------+
    | Method                                              | Description                                              |
    +=====================================================+==========================================================+
    | :meth:`~axis_ratio_rayleigh`                        | Sample axis ratio from Rayleigh distribution             |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~axis_ratio_padilla_strauss`                 | Sample axis ratio from Padilla & Strauss 2008            |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~axis_ratio_uniform`                         | Sample axis ratio from uniform distribution              |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~axis_rotation_angle_uniform`                | Sample axis rotation angle from uniform distribution     |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~lens_redshift_strongly_lensed_numerical`    | Sample lens redshift for strong lensing                  |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~lens_redshift_strongly_lensed_sis_haris`                    | Sample SIS lens redshift (Haris et al. 2018)             |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~velocity_dispersion_gengamma`               | Sample velocity dispersion from gengamma distribution    |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~velocity_dispersion_bernardi`               | Sample velocity dispersion (Bernardi et al. 2010)        |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~velocity_dispersion_ewoud`                  | Sample redshift-dependent velocity dispersion            |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~external_shear_normal`                      | Sample external shear from normal distribution           |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~density_profile_slope_normal`               | Sample density profile slope from normal distribution    |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~optical_depth_sis_analytic`                    | Compute SIS optical depth (Haris et al. 2018)            |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~cross_section_sis`                          | Compute SIS cross-section                                |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~cross_section_sie_feixu`                    | Compute SIE cross-section (Fei Xu et al. 2021)           |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~cross_section_epl_shear_numerical`          | Compute EPL+shear cross-section numerically              |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~cross_section_epl_shear_interpolation`      | Compute EPL+shear cross-section via interpolation        |
    +-----------------------------------------------------+----------------------------------------------------------+

    Instance Attributes
    ----------
    OpticalDepth has the following instance attributes: \n
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | Attribute                                      | Type                         | Unit  | Description                                              |
    +================================================+==============================+=======+==========================================================+
    | :attr:`~npool`                                 | ``int``                      |       | Number of processors for multiprocessing                 |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | :attr:`~z_min`                                 | ``float``                    |       | Minimum redshift                                         |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | :attr:`~z_max`                                 | ``float``                    |       | Maximum redshift                                         |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | :attr:`~cosmo`                                 | ``astropy.cosmology``        |       | Cosmology object                                         |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | :attr:`~lens_type`                             | ``str``                      |       | Type of lens galaxy model                                |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | :attr:`~directory`                             | ``str``                      |       | Directory for interpolator storage                       |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | :attr:`~optical_depth`                         | ``FunctionConditioning``     |       | Optical depth calculator                                 |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | :attr:`~velocity_dispersion`                   | ``FunctionConditioning``     | km/s  | Velocity dispersion sampler                              |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | :attr:`~axis_ratio`                            | ``FunctionConditioning``     |       | Axis ratio sampler                                       |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | :attr:`~axis_rotation_angle`                   | ``FunctionConditioning``     | rad   | Axis rotation angle sampler                              |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | :attr:`~lens_redshift`                         | ``FunctionConditioning``     |       | Lens redshift sampler                                    |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | :attr:`~external_shear`                        | ``FunctionConditioning``     |       | External shear sampler                                   |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | :attr:`~density_profile_slope`                 | ``FunctionConditioning``     |       | Density profile slope sampler                            |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | :attr:`~cross_section`                         | ``callable``                 | rad²  | Cross-section calculator                                 |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | :attr:`~available_lens_samplers`               | ``dict``                     |       | Available lens parameter samplers                        |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | :attr:`~available_lens_functions`              | ``dict``                     |       | Available lens functions                                 |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
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

        # Initialize decision dictionary for creating interpolators
        self._initialize_decision_dictionary(create_new_interpolator, lens_type)
        # Set default lens functions and samplers
        (
            self.lens_param_samplers,
            self.lens_param_samplers_params,
            self.lens_functions,
            self.lens_functions_params,
        ) = self._default_lens_samplers_and_functions(lens_type)

        # Update lens functions and samplers with user input
        self._lens_functions_and_sampler_categorization(
            lens_param_samplers,
            lens_param_samplers_params,
            lens_functions,
            lens_functions_params,
        )

        # Initialize cosmological functions
        self.comoving_distance = comoving_distance(
            z_min=self.z_min,
            z_max=self.z_max,
            cosmo=self.cosmo,
            directory=self.directory,
            create_new=self.create_new_interpolator["comoving_distance"]["create_new"],
            resolution=self.create_new_interpolator["comoving_distance"]["resolution"],
            get_attribute=True,
        )
        self.angular_diameter_distance = angular_diameter_distance(
            z_min=self.z_min,
            z_max=self.z_max,
            cosmo=self.cosmo,
            directory=self.directory,
            create_new=self.create_new_interpolator["angular_diameter_distance"][
                "create_new"
            ],
            resolution=self.create_new_interpolator["angular_diameter_distance"][
                "resolution"
            ],
            get_attribute=True,
        )
        self.angular_diameter_distance_z1z2 = angular_diameter_distance_z1z2(
            z_min=self.z_min,
            z_max=self.z_max,
            cosmo=self.cosmo,
            directory=self.directory,
            create_new=self.create_new_interpolator["angular_diameter_distance_z1z2"][
                "create_new"
            ],
            resolution=self.create_new_interpolator["angular_diameter_distance_z1z2"][
                "resolution"
            ],
            get_attribute=True,
        )
        self.differential_comoving_volume = differential_comoving_volume(
            z_min=self.z_min,
            z_max=self.z_max,
            cosmo=self.cosmo,
            directory=self.directory,
            create_new=self.create_new_interpolator["differential_comoving_volume"][
                "create_new"
            ],
            resolution=self.create_new_interpolator["differential_comoving_volume"][
                "resolution"
            ],
            get_attribute=True,
        )

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

    # --------------------
    # Initialization helpers
    # --------------------
    def _default_lens_samplers_and_functions(self, lens_type):
        """
        Helper to get default lens samplers and functions based on lens type.

        Parameters
        ----------
        lens_type : ``str``
            Type of lens galaxy model.

        Returns
        -------
        lens_param_samplers : ``dict``
            Dictionary of sampler names.
        lens_param_samplers_params : ``dict``
            Dictionary of sampler parameters.
        lens_functions : ``dict``
            Dictionary of lens function names.
        lens_functions_params : ``dict``
            Dictionary of lens function parameters.
        """

        if lens_type == "epl_shear_galaxy":
            lens_param_samplers = dict(
                source_redshift_sl="strongly_lensed_source_redshift",
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
                lens_redshift=dict(integration_size=25000, use_multiprocessing=False),
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
                # cross_section_based_sampler=dict(safety_factor=1.2),
                optical_depth=None,
                cross_section=None,
            )
        elif lens_type == "sie_galaxy":
            lens_param_samplers = dict(
                source_redshift_sl="strongly_lensed_source_redshift",
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
                lens_redshift=dict(integration_size=25000, use_multiprocessing=False),
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
                optical_depth=None,
                cross_section=None,
            )
        elif lens_type == "sis_galaxy":
            lens_param_samplers = dict(
                source_redshift_sl="strongly_lensed_source_redshift",
                lens_redshift="lens_redshift_strongly_lensed_sis_haris",
                velocity_dispersion="velocity_dispersion_bernardi",
                axis_ratio="axis_ratio_uniform",
                axis_rotation_angle="axis_rotation_angle_uniform",
                external_shear="external_shear_normal",
                density_profile_slope="density_profile_slope_normal",
                external_shear_sl="external_shear_normal",
                density_profile_slope_sl="density_profile_slope_normal",
            )
            lens_param_samplers_params = dict(
                source_redshift_sl=None,
                lens_redshift=dict(integration_size=25000, use_multiprocessing=False),
                velocity_dispersion=dict(
                    sigma_min=100.0,
                    sigma_max=400.0,
                    alpha=0.94,
                    beta=1.85,
                    phistar=2.099e-2 * self.cosmo.h**3,
                    sigmastar=113.78,
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
                optical_depth="optical_depth_sis_analytic",
                cross_section="cross_section_sis",
            )
            lens_functions_params = dict(
                param_sampler_type=None,
                cross_section_based_sampler=dict(n_prop=200),
                optical_depth=None,
                cross_section=None,
            )
        else:
            raise ValueError(
                "lens_type should be 'epl_shear_galaxy' or 'sie_galaxy' or 'sis_galaxy'"
            )

        return (
            lens_param_samplers,
            lens_param_samplers_params,
            lens_functions,
            lens_functions_params,
        )

    def _initialize_decision_dictionary(self, create_new_interpolator, lens_type):
        """
        Helper to initialize decision dictionary for creating interpolators.

        Parameters
        ----------
        create_new_interpolator : ``dict`` or ``bool``
            Configuration for creating new interpolators.
        lens_type : ``str``
            Type of lens galaxy model.
        """

        # initialize the interpolator's parameters
        self.create_new_interpolator = dict(
            velocity_dispersion=dict(
                create_new=False, resolution=500, zl_resolution=48
            ),
            axis_ratio=dict(create_new=False, resolution=500, sigma_resolution=48),
            lens_redshift=dict(create_new=False, resolution=48, zl_resolution=48),
            lens_redshift_intrinsic=dict(create_new=False, resolution=500),
            optical_depth=dict(create_new=False, resolution=48),
            comoving_distance=dict(create_new=False, resolution=500),
            angular_diameter_distance=dict(create_new=False, resolution=500),
            angular_diameter_distance_z1z2=dict(create_new=False, resolution=500),
            differential_comoving_volume=dict(create_new=False, resolution=500),
            # density_profile_slope=dict(create_new=False, resolution=100),
            lens_parameters_kde_sl=dict(create_new=False, resolution=5000),
            source_redshift_sl=dict(create_new=False, resolution=500),
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
        elif create_new_interpolator:
            for key in self.create_new_interpolator.keys():
                self.create_new_interpolator[key]["create_new"] = True

    def _lens_functions_and_sampler_categorization(
        self,
        lens_param_samplers,
        lens_param_samplers_params,
        lens_functions,
        lens_functions_params,
    ):
        """
        Helper to categorize and update lens functions and samplers.

        Parameters
        ----------
        lens_param_samplers : ``dict`` or ``None``
            User-provided sampler names or functions.
        lens_param_samplers_params : ``dict`` or ``None``
            User-provided sampler parameters.
        lens_functions : ``dict`` or ``None``
            User-provided lens function names or functions.
        lens_functions_params : ``dict`` or ``None``
            User-provided lens function parameters.
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
                    dict_ = self.available_lens_samplers[
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
                    dict_ = self.available_lens_functions[name]
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

    # ---------------------------------
    # Lens redshift sampler functions
    # ---------------------------------
    def lens_redshift_strongly_lensed_numerical(
        self, size=1000, zs=None, get_attribute=False, **kwargs
    ):
        """
        Sample lens redshifts conditioned on strong lensing (numerical method).

        This method computes the lens redshift distribution by numerically
        integrating over the velocity dispersion distribution (galaxy density distribution wrt), cross-section and differential comoving volume.

        Parameters
        ----------
        size : ``int``
            Number of samples to generate. \\n
            default: 1000
        zs : ``numpy.ndarray``
            Source redshifts.
        get_attribute : ``bool``
            If True, returns the sampler object instead of samples. \\n
            default: False
        **kwargs : ``dict``
            Additional parameters.

        Returns
        -------
        zl : ``numpy.ndarray`` or ``FunctionConditioning``
            Lens redshift samples or sampler object.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> zl = od.lens_redshift(size=100, zs=np.ones(100)*2.0)
        """

        identifier_dict = {}
        identifier_dict["name"] = (
            "lens_redshift_strongly_lensed_numerical_" + self.lens_type
        )
        identifier_dict["resolution"] = self.create_new_interpolator["lens_redshift"][
            "resolution"
        ]
        identifier_dict["zl_resolution"] = self.create_new_interpolator[
            "lens_redshift"
        ]["zl_resolution"]
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

        param_dict = self.available_lens_samplers["lens_redshift"][
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
            number_density = self._helper_number_density_calculation(
                zl_scaled, zs_array
            )
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
        Helper to compute lens redshift distribution using njitted functions.

        Parameters
        ----------
        zl_scaled : ``numpy.ndarray``
            2D array of lens redshifts scaled by source redshift.
        zs : ``numpy.ndarray``
            1D array of source redshifts.

        Returns
        -------
        density_array : ``numpy.ndarray``
            2D array of the lens redshift distribution.
        """

        # Check if the function can be JIT compiled
        # Get random variable samplers from initialized objects
        sigma_rvs_ = self.velocity_dispersion.rvs
        q_rvs_ = self.axis_ratio.rvs
        phi_rvs = self.axis_rotation_angle.rvs
        gamma_rvs = self.density_profile_slope.rvs
        shear_rvs = self.external_shear.rvs
        sigma_pdf_ = self.velocity_dispersion.pdf
        number_density_ = self.velocity_dispersion.function
        cross_section_ = self.cross_section

        from .sampler_functions import _njit_checks

        use_njit_sampler, sampler_dict = _njit_checks(
            sigma_rvs_,
            q_rvs_,
            phi_rvs,
            gamma_rvs,
            shear_rvs,
            sigma_pdf_,
            number_density_,
            cross_section_,
        )

        use_multiprocessing = self.lens_param_samplers_params["lens_redshift"][
            "use_multiprocessing"
        ]

        if use_multiprocessing or not use_njit_sampler:
            print("Using multiprocessing")
            density_array = self._lens_redshift_multiprocessing(
                zl_scaled, zs, sampler_dict
            )
        else:
            print("Using multithreaded njit")
            density_array = self._lens_redshift_multithreaded_njit(
                zl_scaled, zs, sampler_dict
            )

        return density_array

    def _lens_redshift_multithreaded_njit(self, zl_scaled, zs, sampler_dict):

        q_rvs = sampler_dict["q_rvs"]
        phi_rvs = sampler_dict["phi_rvs"]
        gamma_rvs = sampler_dict["gamma_rvs"]
        shear_rvs = sampler_dict["shear_rvs"]
        number_density = sampler_dict["number_density"]
        cross_section_function = sampler_dict["cross_section_function"]

        sigma_min = self.lens_param_samplers_params["velocity_dispersion"]["sigma_min"]
        sigma_max = self.lens_param_samplers_params["velocity_dispersion"]["sigma_max"]

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
            zs,  # 1D
            zl_scaled,  # 2D
            sigma_min,
            sigma_max,
            q_rvs,
            phi_rvs,
            gamma_rvs,
            shear_rvs,
            number_density,
            cross_section_function,
            dVcdz_function,
            integration_size,
        )

        return density_array

    def _lens_redshift_multiprocessing(self, zl_scaled, zs, sampler_dict):
        """
        Compute the lens redshift distribution using multiprocessing.

        Parameters
        ----------
        zl_scaled : ``numpy.ndarray``
            2D array of lens redshifts scaled by source redshift.
        zs : ``numpy.ndarray``
            1D array of source redshifts.

        Returns
        -------
        density_array : ``numpy.ndarray``
            2D array of the lens redshift distribution.
        """

        # 0. zs
        # 1. zl_scaled

        number_density = sampler_dict["number_density"]
        q_rvs = sampler_dict["q_rvs"]
        phi_rvs = sampler_dict["phi_rvs"]
        gamma_rvs = sampler_dict["gamma_rvs"]
        shear_rvs = sampler_dict["shear_rvs"]
        cross_section_function = sampler_dict["cross_section_function"]

        sigma_args = [
            self.velocity_dispersion.info["sigma_min"],
            self.velocity_dispersion.info["sigma_max"],
            number_density,
        ]

        dVcdz_args = self.differential_comoving_volume.function

        integration_size = (
            self.lens_param_samplers_params["lens_redshift"]["integration_size"]
            if "integration_size" in self.lens_param_samplers_params["lens_redshift"]
            else 20000
        )

        # save input params to be passed to multiprocessing
        from ler.utils import save_pickle

        input_params = np.array(
            [
                sigma_args,
                q_rvs,
                dVcdz_args,
                cross_section_function,
                phi_rvs,
                shear_rvs,
                gamma_rvs,
                integration_size,
            ],
            dtype=object,
        )
        save_pickle("input_params_mp.pkl", input_params)

        # setup input params for multiprocessing (zs, zl_scaled, worker_index)
        input_params = np.array(
            [
                (zs[i], zl_scaled[i], i)  # worker index for result ordering
                for i in range(len(zs))
            ],
            dtype=object,
        )

        print("Computing lens redshift distribution with multiprocessing...")
        from .mp import lens_redshift_strongly_lensed_mp

        density_array = np.zeros_like(zl_scaled)

        with Pool(processes=self.npool) as pool:
            for result in tqdm(
                pool.imap_unordered(lens_redshift_strongly_lensed_mp, input_params),
                total=len(zs),
                ncols=100,
                disable=False,
            ):
                (
                    iter_i,
                    density_,
                ) = result

                density_array[iter_i] = density_

        # # simple for loop check
        # for i in tqdm(range(len(zs)), ncols=100):
        #     result = lens_redshift_strongly_lensed_mp(input_params[i])
        #     iter_i, density_ = result
        #     density_array[iter_i] = density_

        # delete the pickle file
        import os

        os.remove("input_params_mp.pkl")

        return density_array

    def lens_redshift_intrinsic(
        self, size=1000, zs=None, get_attribute=False, **kwargs
    ):
        """
        Sample intrinsic lens redshifts (numerical method).

         This method computes the lens redshift distribution by numerically
         integrating over the velocity dispersion distribution (galaxy density distribution wrt) and differential comoving volume.

         Parameters
         ----------
         size : ``int``
             Number of samples to generate. \n
             default: 1000
         zs : ``numpy.ndarray`` or ``None``
             Source redshifts (not used, kept for API compatibility).
         get_attribute : ``bool``
             If True, returns the sampler object instead of samples. \n
             default: False

         Returns
         -------
         zl : ``numpy.ndarray`` or ``FunctionConditioning``
             Lens redshift samples or sampler object.
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
                integrand = (
                    lambda sigma, z: self.velocity_dispersion.function(
                        np.array([sigma])
                    )[0]
                    * self.differential_comoving_volume.function(np.array([z]))[0]
                )
            else:
                integrand = (
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

    def lens_redshift_strongly_lensed_sis_haris(
        self, size, zs, get_attribute=False, **kwargs
    ):
        """
        Sample SIS lens redshifts using Haris et al. (2018) distribution.

        Parameters
        ----------
        size : ``int``
            Number of samples to generate.
        zs : ``numpy.ndarray``
            Source redshifts.
        get_attribute : ``bool``
            If True, returns the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Additional parameters.

        Returns
        -------
        zl : ``numpy.ndarray`` or ``FunctionConditioning``
            Lens redshift samples or sampler object.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_type='sis_galaxy')
        >>> zl = od.lens_redshift(size=100, zs=np.ones(100)*2.0)
        """

        # # Old way
        # splineDc = self.comoving_distance.function_spline  # spline coefficients for the comoving distance and redshifts
        # splineDcInv = self.comoving_distance.function_inverse_spline # spline coefficients for the redshifts and comoving distance
        # u = np.linspace(0, 1, 500)
        # cdf = (10 * u**3 - 15 * u**4 + 6 * u**5)  # See the integral of Eq. A7 of https://arxiv.org/pdf/1807.07062.pdf (cdf)
        # zs = np.array([zs]).reshape(-1)

        # New way
        identifier_dict = {"name": "lens_redshift_strongly_lensed_sis_haris"}
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        # Harris et al. (2018) choice of velocity dispersion is from Choi et al. (2007)
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
        param_dict = self.available_lens_samplers["lens_redshift"][
            "lens_redshift_strongly_lensed_sis_haris"
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
        # pdf_norm_const = zl_object.pdf_norm_const
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
                # / pdf_norm_const # pdf_norm_const is 1
            )

        # @njit
        # def zl_function(zl, zs):
        #     r = zl / zs
        #     return cubic_spline_interpolator(r, function_spline_zl, x_array_zl)

        # zl_object.function = zl_function
        zl_object.function = None
        zl_object.pdf = zl_pdf
        zl_object.rvs = zl_rvs

        return zl_object if get_attribute else zl_object.rvs(size, zs)

    # ---------------------------------------
    # Velocity dispersion sampler functions
    # ---------------------------------------
    def velocity_dispersion_gengamma(self, size, get_attribute=False, **kwargs):
        """
        Sample velocity dispersion from generalized gamma distribution.

        Parameters
        ----------
        size : ``int``
            Number of samples to generate.
        get_attribute : ``bool``
            If True, returns the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Additional parameters. \n
            +-------------------+-------------------+---------------------------------------+
            | Parameter         | Unit              | Description                           |
            +===================+===================+=======================================+
            | sigma_min         | km/s              | minimum velocity dispersion           |
            +-------------------+-------------------+---------------------------------------+
            | sigma_max         | km/s              | maximum velocity dispersion           |
            +-------------------+-------------------+---------------------------------------+
            | alpha             | dimensionless     | Power-law index governing the slope   |
            |                   |                   | of the distribution at low velocities |
            +-------------------+-------------------+---------------------------------------+
            | beta              | dimensionless     | Exponential parameter determining the |
            |                   |                   | sharpness of the high-velocity cutoff |
            +-------------------+-------------------+---------------------------------------+
            | phistar           | h^3 Mpc^-3        | Normalization constant representing   |
            |                   |                   | the comoving number density of galaxy |
            +-------------------+-------------------+---------------------------------------+
            | sigmastar         | km/s              | Characteristic velocity scale marking |
            |                   |                   | the "knee" transition in the VDF      |
            +-------------------+-------------------+---------------------------------------+

        Returns
        -------
        sigma : ``numpy.ndarray`` or ``FunctionConditioning``
            Velocity dispersion samples (km/s) or sampler object.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_param_samplers=dict(
        ...     velocity_dispersion="velocity_dispersion_gengamma"))
        >>> sigma = od.velocity_dispersion(size=100)
        """

        identifier_dict = {"name": "velocity_dispersion_gengamma"}
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["resolution"] = self.create_new_interpolator[
            "velocity_dispersion"
        ]["resolution"]
        param_dict = self.available_lens_samplers["velocity_dispersion"][
            "velocity_dispersion_gengamma"
        ]
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

        from .sampler_functions import velocity_dispersion_gengamma_density_function

        density_func_ = lambda sigma_: velocity_dispersion_gengamma_density_function(
            sigma=sigma_,
            sigma_min=identifier_dict["sigma_min"],
            sigma_max=identifier_dict["sigma_max"],
            alpha=identifier_dict["alpha"],
            beta=identifier_dict["beta"],
            phistar=identifier_dict["phistar"],
            sigmastar=identifier_dict["sigmastar"],
        )

        sigma_object = FunctionConditioning(
            function=density_func_,
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
        Sample velocity dispersion from Bernardi et al. (2010) distribution.

        Uses inverse transform sampling on the velocity dispersion function.

        Parameters
        ----------
        size : ``int``
            Number of samples to generate.
        get_attribute : ``bool``
            If True, returns the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Additional parameters. \n
            +-------------------+-------------------+---------------------------------------+
            | Parameter         | Unit              | Description                           |
            +===================+===================+=======================================+
            | sigma_min         | km/s              | minimum velocity dispersion           |
            +-------------------+-------------------+---------------------------------------+
            | sigma_max         | km/s              | maximum velocity dispersion           |
            +-------------------+-------------------+---------------------------------------+
            | alpha             | dimensionless     | Power-law index governing the slope   |
            |                   |                   | of the distribution at low velocities |
            +-------------------+-------------------+---------------------------------------+
            | beta              | dimensionless     | Exponential parameter determining the |
            |                   |                   | sharpness of the high-velocity cutoff |
            +-------------------+-------------------+---------------------------------------+
            | phistar           | h^3 Mpc^-3        | Normalization constant representing   |
            |                   |                   | the comoving number density of galaxy |
            +-------------------+-------------------+---------------------------------------+
            | sigmastar         | km/s              | Characteristic velocity scale marking |
            |                   |                   | the "knee" transition in the VDF      |
            +-------------------+-------------------+---------------------------------------+

        Returns
        -------
        sigma : ``numpy.ndarray`` or ``FunctionConditioning``
            Velocity dispersion samples (km/s) or sampler object.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_param_samplers=dict(
        ...     velocity_dispersion="velocity_dispersion_bernardi"))
        >>> sigma = od.velocity_dispersion(size=100)
        """

        identifier_dict = {"name": "velocity_dispersion_bernardi"}
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["resolution"] = self.create_new_interpolator[
            "velocity_dispersion"
        ]["resolution"]
        param_dict = self.available_lens_samplers["velocity_dispersion"][
            "velocity_dispersion_bernardi"
        ]
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

        from .sampler_functions import velocity_dispersion_bernardi_denisty_function

        number_density_function = (
            lambda sigma: velocity_dispersion_bernardi_denisty_function(
                sigma,
                alpha=self.lens_param_samplers_params["velocity_dispersion"]["alpha"],
                beta=self.lens_param_samplers_params["velocity_dispersion"]["beta"],
                phistar=self.lens_param_samplers_params["velocity_dispersion"][
                    "phistar"
                ],
                sigmastar=self.lens_param_samplers_params["velocity_dispersion"][
                    "sigmastar"
                ],
            )
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
        Sample redshift-dependent velocity dispersion from Wempe et al. (2022).

        Uses inverse transform sampling with redshift-dependent velocity
        dispersion function.

        Parameters
        ----------
        size : ``int``
            Number of samples to generate.
        zl : ``numpy.ndarray``
            Lens redshifts.
        get_attribute : ``bool``
            If True, returns the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Additional parameters. \n
            +-------------------+-------------------+---------------------------------------+
            | Parameter         | Unit              | Description                           |
            +===================+===================+=======================================+
            | sigma_min         | km/s              | minimum velocity dispersion           |
            +-------------------+-------------------+---------------------------------------+
            | sigma_max         | km/s              | maximum velocity dispersion           |
            +-------------------+-------------------+---------------------------------------+
            | alpha             | dimensionless     | Power-law index governing the slope   |
            |                   |                   | of the distribution at low velocities |
            +-------------------+-------------------+---------------------------------------+
            | beta              | dimensionless     | Exponential parameter determining the |
            |                   |                   | sharpness of the high-velocity cutoff |
            +-------------------+-------------------+---------------------------------------+
            | phistar           | h^3 Mpc^-3        | Normalization constant representing   |
            |                   |                   | the comoving number density of galaxy |
            +-------------------+-------------------+---------------------------------------+
            | sigmastar         | km/s              | Characteristic velocity scale marking |
            |                   |                   | the "knee" transition in the VDF      |
            +-------------------+-------------------+---------------------------------------+


        Returns
        -------
        sigma : ``numpy.ndarray`` or ``FunctionConditioning``
            Velocity dispersion samples (km/s) or sampler object.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> sigma = od.velocity_dispersion(size=100, zl=np.ones(100)*0.5)
        """

        identifier_dict = {"name": "velocity_dispersion_ewoud"}
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["resolution"] = self.create_new_interpolator[
            "velocity_dispersion"
        ]["resolution"]
        identifier_dict["zl_resolution"] = self.create_new_interpolator[
            "velocity_dispersion"
        ]["zl_resolution"]
        param_dict = self.available_lens_samplers["velocity_dispersion"][
            "velocity_dispersion_ewoud"
        ].copy()
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
        from .sampler_functions import velocity_dispersion_ewoud_denisty_function

        number_density_function = (
            lambda sigma, zl: velocity_dispersion_ewoud_denisty_function(
                sigma,
                zl,
                alpha=self.lens_param_samplers_params["velocity_dispersion"]["alpha"],
                beta=self.lens_param_samplers_params["velocity_dispersion"]["beta"],
                phistar=self.lens_param_samplers_params["velocity_dispersion"][
                    "phistar"
                ],
                sigmastar=self.lens_param_samplers_params["velocity_dispersion"][
                    "sigmastar"
                ],
            )
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

    # ------------------------------
    # Axis ratio sampler functions
    # ------------------------------
    def axis_ratio_rayleigh(self, size, sigma, get_attribute=False, **kwargs):
        """
        Sample axis ratio from Rayleigh distribution conditioned on velocity dispersion.

        Parameters
        ----------
        size : ``int``
            Number of samples to generate.
        sigma : ``numpy.ndarray``
            Velocity dispersion of the lens galaxy (km/s).
        get_attribute : ``bool``
            If True, returns the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Additional parameters (q_min, q_max).

        Returns
        -------
        q : ``numpy.ndarray`` or ``FunctionConditioning``
            Axis ratio samples or sampler object.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_param_samplers=dict(axis_ratio="axis_ratio_rayleigh"))
        >>> q = od.axis_ratio(size=100, sigma=np.ones(100)*200.)
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
        identifier_dict["sigma_resolution"] = self.create_new_interpolator[
            "axis_ratio"
        ]["sigma_resolution"]
        param_dict = self.available_lens_samplers["axis_ratio"]["axis_ratio_rayleigh"]
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

        from .sampler_functions import axis_ratio_rayleigh_pdf

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
        Sample axis ratio from Padilla & Strauss (2008) distribution.

        Parameters
        ----------
        size : ``int``
            Number of samples to generate. \n
            default: 1000
        get_attribute : ``bool``
            If True, returns the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Additional parameters (q_min, q_max).

        Returns
        -------
        q : ``numpy.ndarray`` or ``FunctionConditioning``
            Axis ratio samples or sampler object.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_param_samplers=dict(axis_ratio="axis_ratio_padilla_strauss"))
        >>> q = od.axis_ratio(size=100)
        """

        identifier_dict = {"name": "axis_ratio_padilla_strauss"}
        identifier_dict["resolution"] = self.create_new_interpolator["axis_ratio"][
            "resolution"
        ]
        param_dict = self.available_lens_samplers["axis_ratio"][
            "axis_ratio_padilla_strauss"
        ]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        # Using Padilla and Strauss 2008 distribution for axis ratio
        from .sampler_functions import axis_ratio_padilla_strauss_pdf

        q_array = np.linspace(
            identifier_dict["q_min"],
            identifier_dict["q_max"],
            identifier_dict["resolution"],
        )
        pdf = axis_ratio_padilla_strauss_pdf(q_array)

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

    def axis_rotation_angle_uniform(self, size, get_attribute=False, **kwargs):
        """
        Sample axis rotation angle from uniform distribution.

        Parameters
        ----------
        size : ``int``
            Number of samples to generate.
        get_attribute : ``bool``
            If True, returns the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Additional parameters (phi_min, phi_max).

        Returns
        -------
        phi : ``numpy.ndarray`` or ``FunctionConditioning``
            Axis rotation angle samples (rad) or sampler object.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> phi = od.axis_rotation_angle(size=100)
        """

        identifier_dict = {"name": "axis_rotation_angle_uniform"}
        param_dict = self.available_lens_samplers["axis_rotation_angle"][
            "axis_rotation_angle_uniform"
        ]
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
        Sample axis ratio from uniform distribution.

        Parameters
        ----------
        size : ``int``
            Number of samples to generate.
        get_attribute : ``bool``
            If True, returns the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Additional parameters (q_min, q_max).

        Returns
        -------
        q : ``numpy.ndarray`` or ``FunctionConditioning``
            Axis ratio samples or sampler object.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_param_samplers=dict(axis_ratio="axis_ratio_uniform"))
        >>> q = od.axis_ratio(size=100)
        """

        identifier_dict = {"name": "axis_ratio_uniform"}
        param_dict = self.available_lens_samplers["axis_ratio"]["axis_ratio_uniform"]
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
        Sample external shear parameters from 2D normal distribution.

        Parameters
        ----------
        size : ``int``
            Number of samples to generate.
        get_attribute : ``bool``
            If True, returns the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Additional parameters (mean, std).

        Returns
        -------
        shear : ``numpy.ndarray`` or ``FunctionConditioning``
            Array of shape (2, size) with gamma1, gamma2 or sampler object.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> gamma1, gamma2 = od.external_shear(size=100)
        """

        identifier_dict = {"name": "external_shear_normal"}
        param_dict = self.available_lens_samplers["external_shear"][
            "external_shear_normal"
        ]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        mean = param_dict["mean"]
        std = param_dict["std"]
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
        Sample density profile slope from normal distribution.

        Parameters
        ----------
        size : ``int``
            Number of samples to generate.
        get_attribute : ``bool``
            If True, returns the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Additional parameters (mean, std).

        Returns
        -------
        gamma : ``numpy.ndarray`` or ``FunctionConditioning``
            Density profile slope samples or sampler object.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> gamma = od.density_profile_slope(size=100)
        """

        identifier_dict = {"name": "density_profile_slope_normal"}

        param_dict = self.available_lens_samplers["density_profile_slope"][
            "density_profile_slope_normal"
        ]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        mean = param_dict["mean"]
        std = param_dict["std"]
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

        return slope_object if get_attribute else slope_object.rvs(size)

    # --------------------
    # Lens functions
    # --------------------
    def optical_depth_numerical(self, zs, get_attribute=False, **kwargs):
        """
        Helper to compute optical depth numerically by integrating lens redshift.

        Parameters
        ----------
        zs : ``numpy.ndarray``
            Source redshifts.
        get_attribute : ``bool``
            If True, returns the function object instead of values. \\n
            default: False
        **kwargs : ``dict``
            Additional parameters.

        Returns
        -------
        tau : ``numpy.ndarray`` or ``FunctionConditioning``
            Optical depth values or function object.
        """
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

        param_dict = self.available_lens_functions["optical_depth"][
            "optical_depth_numerical"
        ]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        z_min = self.z_min if self.z_min > 0.0 else 0.0001
        z_max = self.z_max
        resolution = identifier_dict["resolution"]
        # zs_array = np.geomspace(z_min, z_max, resolution)
        # z_min = self.z_min + 0.001 if self.z_min == 0.0 else self.z_min
        # z_max = self.z_max
        # z_resolution = identifier_dict["resolution"]
        zs_array = redshift_optimal_spacing(z_min, z_max, resolution)

        def tau(zs):
            # self.lens_redshift.function gives cross-section
            integrand = lambda zl_, zs_: self.lens_redshift.function(
                np.array([zl_]), np.array([zs_])
            )[0]
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

    def optical_depth_sis_analytic(self, zs, get_attribute=False, **kwargs):
        """
        Function to compute the strong lensing optical depth (SIS).
        LambdaCDM(H0=70, Om0=0.3, Ode0=0.7) was used to derive the following equation. This is the analytic version of optical depth from z=0 to z=zs.

        Parameters
        ----------
        zs : ``numpy.ndarray``
            Source redshifts.
        get_attribute : ``bool``
            If True, returns the function object instead of values. \\n
            default: False
        **kwargs : ``dict``
            Additional parameters.

        Returns
        -------
        tau : ``numpy.ndarray`` or ``FunctionConditioning``
            Optical depth values or function object.
        """

        identifier_dict = {"name": "optical_depth_sis_analytic"}
        identifier_dict["z_min"] = self.z_min if self.z_min > 0.0 else 0.0001
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["resolution"] = self.create_new_interpolator["optical_depth"][
            "resolution"
        ]
        param_dict = self.available_lens_functions["optical_depth"][
            "optical_depth_sis_analytic"
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

    # ------------------------
    # cross section functions
    # ------------------------
    def cross_section_sis(
        self, zs=None, zl=None, sigma=None, get_attribute=False, **kwargs
    ):
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

    def cross_section_sie_feixu(
        self, zs=None, zl=None, sigma=None, q=None, get_attribute=False, **kwargs
    ):
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
            _cross_section_sis = self.cross_section_sis(
                sigma=None, zl=None, zs=None, get_attribute=True
            )

            @njit
            def _cross_section_sie_feixu(zs, zl, sigma, q):
                return phi_cut_SIE(q) * _cross_section_sis(zs, zl, sigma)

            return _cross_section_sie_feixu
        else:
            _cross_section_sis = self.cross_section_sis(
                zs=zs, zl=zl, sigma=sigma, get_attribute=False
            )
            return phi_cut_SIE(q) * _cross_section_sis

    def cross_section_epl_shear_numerical(
        self,
        zs,
        zl,
        sigma,
        q,
        phi,
        gamma,
        gamma1,
        gamma2,
    ):
        """
        Function to compute the strong lensing cross-section numerically for EPL + external shear lenses.

        Parameters
        ----------
        zs : `numpy.ndarray`
            redshift of the source galaxies
        zl : `numpy.ndarray`
            redshift of the lens galaxies
        sigma : `numpy.ndarray`
            velocity dispersion of the lens galaxies
        q : `numpy.ndarray`
            axis ratios of the lens galaxies
        phi : `numpy.ndarray`
            axis rotation angles of the lens galaxies in radians
        gamma : `numpy.ndarray`
            external shear magnitudes of the lens galaxies
        gamma1 : `numpy.ndarray`
            external shear component 1 of the lens galaxies
        gamma2 : `numpy.ndarray`
            external shear component 2 of the lens galaxies

        Returns
        -------
        cross_section : `numpy.ndarray`
            strong lensing cross-section of the lens galaxies in square radians
        """

        size = len(zs)
        e1, e2 = phi_q2_ellipticity(phi, q)
        theta_E = self.compute_einstein_radii(sigma, zl, zs)

        area_list = np.ones(size)
        for i in range(size):
            area_list[i] = cross_section(
                theta_E[i], e1[i], e2[i], gamma[i], gamma1[i], gamma2[i]
            )

        return area_list

    def cross_section_epl_shear_numerical_mp(
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
            e1, e2 = phi_q2_ellipticity(phi, q)

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
        """
        Create a parameter grid for lens galaxies.

        Parameters
        ----------
        size_list : list
            List of sizes for each parameter grid.

        Returns
        -------
        zl : numpy.ndarray
            Lens redshifts.
        sigma : numpy.ndarray
            Velocity dispersions.
        q : numpy.ndarray
            Axis ratios.
        """

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

        e1, e2 = phi_q2_ellipticity(phi, q)

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
        e1, e2 = phi_q2_ellipticity(phi, q)
        e1 = np.sort(e1)
        e2 = np.sort(e2)

        gamma = np.linspace(gamma_min, gamma_max, size_list[2])

        gamma1 = np.linspace(gamma1_min, gamma1_max, size_list[3])
        gamma2 = np.linspace(gamma2_min, gamma2_max, size_list[4])

        return q, phi, e1, e2, gamma, gamma1, gamma2

    def cross_section_epl_shear_interpolation_init(self, file_path, size_list):
        print(f"Cross section interpolation data points will be created at {file_path}")

        # ----------------------------------
        # cross section unit to cross section ratio fitting
        # ----------------------------------
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

        # Sample individual parameter values (not grid)
        try:
            q = self.axis_ratio.rvs(size, sigma)
        except TypeError:
            q = self.axis_ratio.rvs(size)
        phi = self.axis_rotation_angle.rvs(size)
        e1, e2 = phi_q2_ellipticity(phi, q)
        gamma = self.density_profile_slope.rvs(size)
        gamma1, gamma2 = self.external_shear.rvs(size)

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

        cs = self.cross_section_epl_shear_numerical_mp(
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
            ],
        )

        return e1, e2, gamma, gamma1, gamma2, cs

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
            velocity_dispersion=dict(
                sigma_min=self.lens_param_samplers_params["velocity_dispersion"][
                    "sigma_min"
                ],
                sigma_max=self.lens_param_samplers_params["velocity_dispersion"][
                    "sigma_max"
                ],
            ),
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
        Strong lensing optical depth calculator.

        Returns
        -------
        optical_depth : ``FunctionConditioning``
            Function object with `.function(zs)` method that returns \\n
            optical depth for given source redshifts.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> tau = od.optical_depth.function(np.array([1.0, 2.0]))
        """

        return self._optical_depth

    @optical_depth.setter
    def optical_depth(self, input_function):
        # Handle available lens functions (need to map to private method)
        available_functions_mapping = {
            "optical_depth_numerical": "optical_depth_numerical",
        }

        if input_function in self.available_lens_functions["optical_depth"]:
            print(f"using ler available optical depth function : {input_function}")
            # Map to private method if it exists
            method_name = available_functions_mapping.get(
                input_function, input_function
            )
            args = self.lens_functions_params["optical_depth"]
            if args is None:
                self._optical_depth = getattr(self, method_name)(
                    zs=None, get_attribute=True
                )
            else:
                self._optical_depth = getattr(self, method_name)(
                    zs=None, get_attribute=True, **args
                )
        elif isinstance(input_function, FunctionConditioning):
            print(
                "using user provided custom optical depth class/object of type ler.utils.FunctionConditioning"
            )
            self._optical_depth = input_function
        elif callable(input_function):
            print("using user provided custom optical depth function")
            self._optical_depth = FunctionConditioning(
                function=None, x_array=None, create_function=input_function
            )
        else:
            raise ValueError(
                "optical depth function should be string in available_lens_functions['optical_depth'] or callable or class object of 'ler.utils.FunctionConditioning'"
            )

    @property
    def velocity_dispersion(self):
        """
        Velocity dispersion sampler object.

        Returns a ``FunctionConditioning`` object with methods: \n
        - ``rvs(size, zl)``: Sample velocity dispersion values \n
        - ``pdf(sigma, zl)``: Get probability density \n
        - ``function(sigma, zl)``: Get number density function \n

        Returns
        -------
        velocity_dispersion : ``FunctionConditioning``
            Sampler object for velocity dispersion (km/s).

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> sigma = od.velocity_dispersion(size=100, zl=np.ones(100)*0.5)
        """

        return self._velocity_dispersion

    @velocity_dispersion.setter
    def velocity_dispersion(self, prior):
        if prior in self.available_lens_samplers["velocity_dispersion"]:
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
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user provided custom velocity_dispersion class/object of type ler.utils.FunctionConditioning"
            )
            self._velocity_dispersion = prior
        elif callable(prior):
            print("using user provided custom velocity_dispersion function")
            identifier_dict = {"name": "velocity_dispersion_custom"}
            identifier_dict["sigma_min"] = self.lens_param_samplers_params[
                "velocity_dispersion"
            ]["sigma_min"]
            identifier_dict["sigma_max"] = self.lens_param_samplers_params[
                "velocity_dispersion"
            ]["sigma_max"]
            identifier_dict["resolution"] = self.create_new_interpolator[
                "velocity_dispersion"
            ]["resolution"]

            # setting up inputs for the interpolator
            sigma_array = np.linspace(
                identifier_dict["sigma_min"],
                identifier_dict["sigma_max"],
                identifier_dict["resolution"],
            )

            if prior.__code__.co_argcount == 2:
                identifier_dict["z_min"] = self.z_min
                identifier_dict["z_max"] = self.z_max
                identifier_dict["zl_resolution"] = self.create_new_interpolator[
                    "velocity_dispersion"
                ]["zl_resolution"]
                z_min = self.z_min + 0.001 if self.z_min == 0.0 else self.z_min
                z_max = self.z_max
                z_resolution = identifier_dict["zl_resolution"]
                zl_array = redshift_optimal_spacing(z_min, z_max, z_resolution)

                number_density_function = lambda sigma, zl: prior(
                    sigma,
                    zl,
                )
            else:
                number_density_function = lambda sigma, zl: prior(
                    sigma,
                )
                zl_array = None

            self._velocity_dispersion = FunctionConditioning(
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
        else:
            raise ValueError(
                "velocity_dispersion should be sampler name from available_lens_samplers['velocity_dispersion'] or class object of 'ler.utils.FunctionConditioning' or callable function"
            )

    @property
    def axis_ratio(self):
        """
        Axis ratio sampler object.

        Returns a ``FunctionConditioning`` object with methods: \n
        - ``rvs(size, sigma)``: Sample axis ratio values \n
        - ``pdf(q, sigma)``: Get probability density \n
        - ``function(q, sigma)``: Get distribution function \n

        Returns
        -------
        axis_ratio : ``FunctionConditioning``
            Sampler object for axis ratio.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> q = od.axis_ratio(size=100, sigma=np.ones(100)*200.)
        """

        return self._axis_ratio

    @axis_ratio.setter
    def axis_ratio(self, prior):
        if prior in self.available_lens_samplers["axis_ratio"]:
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
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user provided custom axis_ratio class/object of type ler.utils.FunctionConditioning"
            )
            self._axis_ratio = prior
        elif callable(prior):
            print("using user provided custom axis_ratio sampler function")
            self._axis_ratio = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior
            )
        else:
            raise ValueError(
                "axis_ratio should be string in available_lens_samplers['axis_ratio'] or class object of 'ler.utils.FunctionConditioning' or callable function with input argument 'size'"
            )

    @property
    def lens_redshift(self):
        """
        Lens redshift sampler object.

        Returns a ``FunctionConditioning`` object with methods: \n
        - ``rvs(size, zs)``: Sample lens redshifts given source redshifts \n
        - ``pdf(zl, zs)``: Get probability density \n
        - ``function(zl, zs)``: Get effective lensing cross-section \n

        Returns
        -------
        lens_redshift : ``FunctionConditioning``
            Sampler object for lens redshift.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> zl = od.lens_redshift(size=100, zs=np.ones(100)*2.0)
        """

        return self._lens_redshift

    @lens_redshift.setter
    def lens_redshift(self, prior):
        if prior in self.available_lens_samplers["lens_redshift"]:
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
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user provided custom lens_redshift class/object of type ler.utils.FunctionConditioning"
            )
            self._lens_redshift = prior
        else:
            raise ValueError(
                "lens_redshift should be string in available_lens_samplers['lens_redshift'] or class object of 'ler.utils.FunctionConditioning'"
            )

    @property
    def axis_rotation_angle(self):
        """
        Axis rotation angle sampler object.

        Returns a ``FunctionConditioning`` object with methods: \n
        - ``rvs(size)``: Sample axis rotation angles \n
        - ``pdf(phi)``: Get probability density \n

        Returns
        -------
        axis_rotation_angle : ``FunctionConditioning``
            Sampler object for axis rotation angle (rad).

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> phi = od.axis_rotation_angle(size=100)
        """

        return self._axis_rotation_angle

    @axis_rotation_angle.setter
    def axis_rotation_angle(self, prior):
        if prior in self.available_lens_samplers["axis_rotation_angle"]:
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
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user provided custom axis_rotation_angle class/object of type ler.utils.FunctionConditioning"
            )
            self._axis_rotation_angle = prior
        elif callable(prior):
            print("using user provided custom axis_rotation_angle sampler function")
            self._axis_rotation_angle = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior
            )
        else:
            raise ValueError(
                "axis_rotation_angle should be string in available_lens_samplers['axis_rotation_angle'] or class object of 'ler.utils.FunctionConditioning' or callable function with input argument 'size'"
            )

    @property
    def external_shear(self):
        """
        External shear sampler object.

        Returns a ``FunctionConditioning`` object with methods: \n
        - ``rvs(size)``: Sample shear components (gamma1, gamma2) \n
        - ``pdf(gamma1, gamma2)``: Get probability density \n

        Returns
        -------
        external_shear : ``FunctionConditioning``
            Sampler object for external shear.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> gamma1, gamma2 = od.external_shear(size=100)
        """

        return self._external_shear

    @external_shear.setter
    def external_shear(self, prior):
        if prior in self.available_lens_samplers["external_shear"]:
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
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user provided custom external_shear class/object of type ler.utils.FunctionConditioning"
            )
            self._external_shear = prior
        elif callable(prior):
            print("using user provided custom external_shear sampler function")
            self._external_shear = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior
            )
        else:
            raise ValueError(
                "external_shear should be string in available_lens_samplers['external_shear'] or class object of 'ler.utils.FunctionConditioning' or callable function with input argument 'size'"
            )

    @property
    def external_shear_sl(self):
        """
        External shear sampler object (strong lensing conditioned).

        Returns a ``FunctionConditioning`` object with methods: \n
        - ``rvs(size)``: Sample shear components (gamma1, gamma2) \n
        - ``pdf(gamma1, gamma2)``: Get probability density \n

        Returns
        -------
        external_shear_sl : ``FunctionConditioning``
            Sampler object for external shear (strong lensing).

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> gamma1, gamma2 = od.external_shear_sl(size=100)
        """

        return self._external_shear_sl

    @external_shear_sl.setter
    def external_shear_sl(self, prior):
        if prior in self.available_lens_samplers["external_shear_sl"]:
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
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user provided custom external_shear_sl class/object of type ler.utils.FunctionConditioning"
            )
            self._external_shear_sl = prior
        elif callable(prior):
            print("using user provided custom external_shear_sl sampler function")
            self._external_shear_sl = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior
            )
        else:
            raise ValueError(
                "external_shear_sl should be string in available_lens_samplers['external_shear_sl'] or class object of 'ler.utils.FunctionConditioning' or callable function with input argument 'size'"
            )

    @property
    def density_profile_slope(self):
        """
        Density profile slope sampler object.

        Returns a ``FunctionConditioning`` object with methods: \n
        - ``rvs(size)``: Sample density profile slope values \n
        - ``pdf(gamma)``: Get probability density \n

        Returns
        -------
        density_profile_slope : ``FunctionConditioning``
            Sampler object for density profile slope.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> gamma = od.density_profile_slope(size=100)
        """

        return self._density_profile_slope

    @density_profile_slope.setter
    def density_profile_slope(self, prior):
        if prior in self.available_lens_samplers["density_profile_slope"]:
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
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user provided custom density_profile_slope class/object of type ler.utils.FunctionConditioning"
            )
            self._density_profile_slope = prior
        elif callable(prior):
            print("using user provided custom density_profile_slope sampler function")
            self._density_profile_slope = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior
            )
        else:
            raise ValueError(
                "density_profile_slope should be string in available_lens_samplers['density_profile_slope'] or class object of 'ler.utils.FunctionConditioning' or callable function with input argument 'size'"
            )

    @property
    def density_profile_slope_sl(self):
        """
        Density profile slope sampler object (strong lensing conditioned).

        Returns a ``FunctionConditioning`` object with methods: \n
        - ``rvs(size)``: Sample density profile slope values \n
        - ``pdf(gamma)``: Get probability density \n

        Returns
        -------
        density_profile_slope_sl : ``FunctionConditioning``
            Sampler object for density profile slope (strong lensing).

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> gamma = od.density_profile_slope_sl(size=100)
        """

        return self._density_profile_slope_sl

    @density_profile_slope_sl.setter
    def density_profile_slope_sl(self, prior):
        if prior in self.available_lens_samplers["density_profile_slope_sl"]:
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
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user provided custom density_profile_slope_sl class/object of type ler.utils.FunctionConditioning"
            )
            self._density_profile_slope_sl = prior
        elif callable(prior):
            print(
                "using user provided custom density_profile_slope_sl sampler function"
            )
            self._density_profile_slope_sl = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior
            )
        else:
            raise ValueError(
                "density_profile_slope_sl should be string in available_lens_samplers['density_profile_slope_sl'] or class object of 'ler.utils.FunctionConditioning' or callable function with input argument 'size'"
            )

    @property
    def cross_section(self):
        """
        Lensing cross-section calculator.

        Returns a callable that computes lensing cross-section for individual \n
        lensing events. Input parameters depend on lens type: \n
        - EPL+shear: zs, zl, sigma, q, phi, gamma, gamma1, gamma2 \n
        - SIE: zs, zl, sigma, q \n
        - SIS: zs, zl, sigma \n

        Returns
        -------
        cross_section : ``callable``
            Cross-section function (rad² units).

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> cs = od.cross_section(zs=zs, zl=zl, sigma=sigma, ...)
        """
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

    # --------------------
    # Available samplers and functions
    # --------------------
    @property
    def available_lens_samplers(self):
        """
        Dictionary of available lens parameter samplers and their default parameters.

        Returns
        -------
        available_lens_samplers : ``dict``
            Dictionary with sampler names and default parameters.
        """

        self._available_lens_samplers = dict(
            source_redshift_sl=dict(
                strongly_lensed_source_redshift=dict(
                    tau_approximation=True,
                ),
                # strongly_lensed_source_redshift_rjs=dict(
                #     tau_approximation=True,
                # ),
            ),
            lens_redshift=dict(
                lens_redshift_strongly_lensed_sis_haris=None,
                lens_redshift_strongly_lensed_numerical=dict(
                    integration_size=25000, use_multiprocessing=False
                ),
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
                external_shear_normal=dict(mean=0.0, std=0.05),
                external_shear_sl_numerical_hemanta=dict(
                    external_shear_normal=dict(mean=0.0, std=0.05)
                ),
            ),
            density_profile_slope=dict(
                density_profile_slope_normal=dict(mean=1.99, std=0.149),
            ),
            density_profile_slope_sl=dict(
                density_profile_slope_normal=dict(mean=2.091, std=0.133),
                density_profile_slope_sl_numerical_hemanta=dict(
                    density_profile_slope_normal=dict(mean=1.99, std=0.149)
                ),
            ),
            source_parameters=dict(sample_gw_parameters=None),
        )

        return self._available_lens_samplers

    @property
    def available_lens_functions(self):
        """
        Dictionary of available lens functions and their default parameters.

        Returns
        -------
        available_lens_functions : ``dict``
            Dictionary with function names and default parameters.
        """

        self._available_lens_functions = dict(
            cross_section_based_sampler=dict(
                # rejection_sampling_with_cross_section_sie_feixu=None,
                # rejection_sampling_with_cross_section_sis=None,
                rejection_sampling_with_cross_section=dict(safety_factor=1.2),
                importance_sampling_with_cross_section=dict(n_prop=200),
            ),
            optical_depth=dict(
                optical_depth_sis_analytic=None,
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

        return self._available_lens_functions

    # -------------------------------------
    # Simple instance attribute properties
    # -------------------------------------
    @property
    def npool(self):
        """
        Number of processors for multiprocessing.

        Returns
        -------
        npool : ``int``
            Number of parallel processors.
        """
        return self._npool

    @npool.setter
    def npool(self, value):
        self._npool = value

    @property
    def z_min(self):
        """
        Minimum redshift of the lens galaxy population.

        Returns
        -------
        z_min : ``float``
            Minimum redshift.
        """
        return self._z_min

    @z_min.setter
    def z_min(self, value):
        self._z_min = value

    @property
    def z_max(self):
        """
        Maximum redshift of the lens galaxy population.

        Returns
        -------
        z_max : ``float``
            Maximum redshift.
        """
        return self._z_max

    @z_max.setter
    def z_max(self, value):
        self._z_max = value

    @property
    def cosmo(self):
        """
        Cosmology object for distance calculations.

        Returns
        -------
        cosmo : ``astropy.cosmology``
            Cosmology object.
        """
        return self._cosmo

    @cosmo.setter
    def cosmo(self, value):
        self._cosmo = value

    @property
    def lens_type(self):
        """
        Type of lens galaxy model.

        Returns
        -------
        lens_type : ``str``
            Lens type ('epl_shear_galaxy', 'sie_galaxy', or 'sis_galaxy').
        """
        return self._lens_type

    @lens_type.setter
    def lens_type(self, value):
        self._lens_type = value

    @property
    def directory(self):
        """
        Directory for interpolator storage.

        Returns
        -------
        directory : ``str``
            Path to interpolator JSON files.
        """
        return self._directory

    @directory.setter
    def directory(self, value):
        self._directory = value
