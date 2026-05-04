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

from numba import njit, get_num_threads, set_num_threads
from multiprocessing import Pool
import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from astropy.cosmology import LambdaCDM
from tqdm import tqdm

from ..utils import (
    cubic_spline_interpolator,
    inverse_transform_sampler,
    cubic_spline_interpolator2d,
    save_json,
    load_json,
    interpolator_json_path,
    FunctionConditioning,
    inverse_transform_sampler2d,
    pdf_cubic_spline_interpolator2d,
    comoving_distance,
    angular_diameter_distance,
    angular_diameter_distance_z1z2,
    differential_comoving_volume,
    generate_mixed_grid,
)

from .lens_functions import (
    phi_cut_SIE,
    cross_section,
)

from ..image_properties.cross_section_njit import phi_q2_ellipticity

from .mp import cross_section_mp


@njit()
def _seed_numba_rng(seed):
    """Seed Numba's RNG so njit samplers are reproducible across calls."""
    np.random.seed(seed)


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
    lens_priors : ``dict`` or ``None``
        Dictionary of sampler functions for lens parameters. \n
        default: None (uses defaults for lens_type)
    lens_priors_params : ``dict`` or ``None``
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
    | :meth:`~rayleigh`                        | Sample axis ratio from Rayleigh distribution             |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~axis_ratio_padilla_strauss`                 | Sample axis ratio from Padilla & Strauss 2008            |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~axis_ratio_uniform`                         | Sample axis ratio from uniform distribution              |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~uniform`                | Sample axis rotation angle from uniform distribution     |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~lens_redshift_strongly_lensed_numerical`    | Sample lens redshift for strong lensing                  |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~lens_redshift_strongly_lensed_sis_analytical`    | Sample SIS lens redshift (Haris et al. 2018)             |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~gengamma`               | Sample velocity dispersion from gengamma distribution    |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~velocity_dispersion_bernardi`               | Sample velocity dispersion (Bernardi et al. 2010)        |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~velocity_dispersion_ewoud`                  | Sample redshift-dependent velocity dispersion            |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~external_shear_normal`                      | Sample external shear from normal distribution           |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~normal`                          | Sample from normal distribution                          |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~optical_depth_sis_analytic`                 | Compute SIS optical depth (Haris et al. 2018)            |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~cross_section_sis`                          | Compute SIS cross-section                                |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~cross_section_sie_feixu`                    | Compute SIE cross-section (Fei Xu et al. 2021)           |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~cross_section_epl_shear_numerical`          | Compute EPL+shear cross-section numerically              |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~cross_section_epl_shear_interpolation`      | Compute EPL+shear cross-section via interpolation        |
    +-----------------------------------------------------+----------------------------------------------------------+
    | :meth:`~cross_section_epl_shear_njit`               | Compute EPL+shear cross-section using Numba njit         |
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
    | :attr:`~external_shear1`                       | ``FunctionConditioning``     |       | External shear 1 (gamma1) sampler                        |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | :attr:`~external_shear2`                       | ``FunctionConditioning``     |       | External shear 2 (gamma2) sampler                        |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | :attr:`~density_profile_slope`                 | ``FunctionConditioning``     |       | Density profile slope sampler                            |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | :attr:`~cross_section`                         | ``callable``                 | rad²  | Cross-section calculator                                 |
    +------------------------------------------------+------------------------------+-------+----------------------------------------------------------+
    | :attr:`~available_lens_priors`               | ``dict``                     |       | Available lens parameter samplers                        |
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
        lens_priors=None,
        lens_priors_params=None,
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
        self.cosmo = (
            cosmology
            if cosmology
            else LambdaCDM(
                H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0
            )
        )
        self.directory = directory

        # Initialize decision dictionary for creating interpolators
        self._initialize_decision_dictionary(create_new_interpolator, lens_type)

        # Update lens functions and samplers with user input
        self._lens_functions_and_sampler_categorization(
            lens_type,
            lens_priors,
            lens_priors_params,
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
        self.velocity_dispersion = self.lens_priors["velocity_dispersion"]
        self.axis_ratio = self.lens_priors["axis_ratio"]
        self.axis_rotation_angle = self.lens_priors["axis_rotation_angle"]
        self.density_profile_slope = self.lens_priors["density_profile_slope"]
        self.external_shear1 = self.lens_priors["external_shear1"]
        self.external_shear2 = self.lens_priors["external_shear2"]

        # cross section initialization
        self.cross_section = self.lens_functions["cross_section"]

        # # lens redshift initialization
        self.lens_redshift_sl = self.lens_priors["lens_redshift_sl"]
        self.lens_redshift = self.lens_priors["lens_redshift"]
        

        # lens function initialization
        self.optical_depth = self.lens_functions["optical_depth"]

    # --------------------
    # Initialization helpers
    # --------------------
    def _default_lens_priors_and_functions(self, lens_type):
        """
        Helper to get default lens samplers and functions based on lens type.

        Parameters
        ----------
        lens_type : ``str``
            Type of lens galaxy model.

        Returns
        -------
        lens_priors : ``dict``
            Dictionary of sampler names.
        lens_priors_params : ``dict``
            Dictionary of sampler parameters.
        lens_functions : ``dict``
            Dictionary of lens function names.
        lens_functions_params : ``dict``
            Dictionary of lens function parameters.
        """

        if lens_type == "epl_shear_galaxy":
            lens_priors = dict(
                zs_sl="strongly_lensed_source_redshift",
                lens_redshift_sl="lens_redshift_strongly_lensed_numerical",
                lens_redshift="lens_redshift_intrinsic_numerical",
                velocity_dispersion="velocity_dispersion_ewoud",
                axis_ratio="rayleigh",
                axis_rotation_angle="uniform",
                external_shear1="normal",
                external_shear2="normal",
                density_profile_slope="normal",
            )
            lens_priors_params = dict(
                zs_sl=None,
                lens_redshift_sl=dict(
                    param_name = "lens_redshift_sl",
                    sampler_type = "lens_redshift_strongly_lensed_numerical",
                    lens_type = lens_type,
                    integration_size=50000,
                    use_multiprocessing=False,
                    cross_section_epl_shear_interpolation=False,
                ),
                lens_redshift=None,
                velocity_dispersion=dict(
                    param_name='velocity_dispersion', sampler_type='velocity_dispersion_ewoud',
                    sigma_min=100.0,
                    sigma_max=400.0,
                    alpha=0.94,
                    beta=1.85,
                    phistar=2.099e-2 * (self.cosmo.h/0.7)**3,
                    sigmastar=113.78,
                ),
                axis_ratio=dict(
                    param_name='axis_ratio', sampler_type='rayleigh',
                    q_min=0.2, q_max=1.0
                ),
                axis_rotation_angle=dict(
                    param_name='axis_rotation_angle', sampler_type='uniform',
                    x_min=0.0, x_max=2 * np.pi
                ),
                external_shear1=dict(
                    param_name='external_shear1', sampler_type='normal',
                    mu=0.0, sigma=0.05
                ),
                external_shear2=dict(
                    param_name='external_shear2', sampler_type='normal',
                    mu=0.0, sigma=0.05
                ),
                density_profile_slope=dict(
                    param_name='density_profile_slope', sampler_type='normal',
                    mu=1.99, sigma=0.149
                ),
            )
            lens_functions = dict(
                param_sampler_type="epl_shear_sl_parameters_rvs",
                cross_section_based_sampler="importance_sampler_partial",
                optical_depth="optical_depth_numerical",
                cross_section="cross_section_epl_shear_njit",
            )
            lens_functions_params = dict(
                param_sampler_type=None,
                cross_section_based_sampler=dict(
                    n_prop=50,
                    threshold_factor=1e-4,
                    sigma_min=100.0,
                    sigma_max=400.0,
                    q_min=0.2,
                    q_max=1.0,
                    phi_min=0.0,
                    phi_max=2 * np.pi,
                    gamma_min=1.4,
                    gamma_max=2.7,
                    shear_min=-0.22,
                    shear_max=0.20,
                ),
                optical_depth=dict(param_name="optical_depth", function_type="optical_depth_numerical"),
                cross_section=dict(
                    num_th=500, maginf=-100.0
                ),
            )
        elif lens_type == "sie_galaxy":
            lens_priors = dict(
                zs_sl="strongly_lensed_source_redshift",
                lens_redshift_sl="lens_redshift_strongly_lensed_numerical",
                lens_redshift="lens_redshift_intrinsic_numerical",
                velocity_dispersion="velocity_dispersion_ewoud",
                axis_ratio="rayleigh",
                axis_rotation_angle="uniform",
                external_shear1="constant_values_n_size",
                external_shear2="constant_values_n_size",
                density_profile_slope="constant_values_n_size",
            )
            lens_priors_params = dict(
                zs_sl=None,
                lens_redshift_sl=dict(
                    param_name = "lens_redshift_sl",
                    sampler_type = "lens_redshift_strongly_lensed_numerical",
                    lens_type = lens_type,
                    integration_size=50000,
                    use_multiprocessing=False,
                    cross_section_epl_shear_interpolation=False,
                ),
                lens_redshift=None,
                velocity_dispersion=dict(
                    param_name='velocity_dispersion', sampler_type='velocity_dispersion_ewoud',
                    sigma_min=100.0,
                    sigma_max=400.0,
                    alpha=0.94,
                    beta=1.85,
                    phistar=2.099e-2 * (self.cosmo.h/0.7)**3,
                    sigmastar=113.78,
                ),
                axis_ratio=dict(
                    param_name='axis_ratio', sampler_type='rayleigh',
                    q_min=0.2, q_max=1.0
                ),
                axis_rotation_angle=dict(
                    param_name='axis_rotation_angle', sampler_type='uniform',
                    x_min=0.0, x_max=2 * np.pi
                ),
                external_shear1=dict(param_name='external_shear1', sampler_type= 'constant_values_n_size', value=0.0),
                external_shear2=dict(param_name='external_shear2', sampler_type= 'constant_values_n_size', value=0.0),
                density_profile_slope=dict(param_name='density_profile_slope', sampler_type='constant_values_n_size', value=2.0),
            )
            lens_functions = dict(
                param_sampler_type="epl_shear_sl_parameters_rvs",
                cross_section_based_sampler="importance_sampler_partial",
                optical_depth="optical_depth_numerical",
                cross_section="cross_section_sie_feixu",
            )
            lens_functions_params = dict(
                param_sampler_type=None,
                cross_section_based_sampler=dict(
                    n_prop=50,
                    threshold_factor=0.0,
                    sigma_min=100.0,
                    sigma_max=400.0,
                    q_min=0.2,
                    q_max=1.0,
                    phi_min=0.0,
                    phi_max=2 * np.pi,
                    gamma_min=2.0,
                    gamma_max=2.0,
                    shear_min=0.0,
                    shear_max=0.0,
                ),
                optical_depth=dict(param_name="optical_depth", function_type="optical_depth_numerical"),
                cross_section=None,
            )
        elif lens_type == "sis_galaxy":
            lens_priors = dict(
                zs_sl="strongly_lensed_source_redshift",
                lens_redshift_sl="lens_redshift_strongly_lensed_numerical",
                lens_redshift="lens_redshift_intrinsic_numerical",
                velocity_dispersion="velocity_dispersion_ewoud",
                axis_ratio="constant_values_n_size",
                axis_rotation_angle="constant_values_n_size",
                external_shear1="constant_values_n_size",
                external_shear2="constant_values_n_size",
                density_profile_slope="constant_values_n_size",
            )
            lens_priors_params = dict(
                zs_sl=None,
                lens_redshift_sl=dict(
                    param_name = "lens_redshift_sl",
                    sampler_type = "lens_redshift_strongly_lensed_numerical",
                    lens_type = lens_type,
                    integration_size=50000,
                    use_multiprocessing=False,
                    cross_section_epl_shear_interpolation=False,
                ),
                lens_redshift=None,
                velocity_dispersion=dict(
                    param_name='velocity_dispersion', sampler_type='velocity_dispersion_ewoud',
                    sigma_min=100.0,
                    sigma_max=400.0,
                    alpha=0.94,
                    beta=1.85,
                    phistar=2.099e-2 * (self.cosmo.h/0.7)**3,
                    sigmastar=113.78,
                ),
                axis_ratio=dict(
                    param_name='axis_ratio', sampler_type='constant_values_n_size',
                    value=1.0
                ),
                axis_rotation_angle=dict(
                    param_name='axis_rotation_angle', sampler_type='constant_values_n_size',
                    value=0.0
                ),
                external_shear1=dict(
                    param_name='external_shear1', sampler_type='constant_values_n_size',
                    value=0.0
                ),
                external_shear2=dict(
                    param_name='external_shear2', sampler_type='constant_values_n_size',
                    value=0.0
                ),
                density_profile_slope=dict(
                    param_name='density_profile_slope', sampler_type='constant_values_n_size',
                    value=2.0
                ),
            )
            lens_functions = dict(
                param_sampler_type="epl_shear_sl_parameters_rvs",
                cross_section_based_sampler="importance_sampler_partial",
                optical_depth="optical_depth_numerical",
                cross_section="cross_section_sis",
            )
            lens_functions_params = dict(
                param_sampler_type=None,
                cross_section_based_sampler=dict(
                    n_prop=50,
                    threshold_factor=0.0,
                    sigma_min=100.0,
                    sigma_max=400.0,
                    q_min=1.0,
                    q_max=1.0,
                    phi_min=0.0,
                    phi_max=0.0,
                    gamma_min=2.0,
                    gamma_max=2.0,
                    shear_min=0.0,
                    shear_max=0.0,
                ),
                optical_depth=dict(param_name="optical_depth", function_type="optical_depth_numerical"),
                cross_section=None, # if interpolation, you should provide dict(num_th=500, maginf=-1000.0)
            )
        else:
            raise ValueError(
                "lens_type should be 'epl_shear_galaxy' or 'sie_galaxy' or 'sis_galaxy'"
            )

        return (
            lens_priors,
            lens_priors_params,
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

        # preserve user input and initialize defaults
        user_create_new_interpolator = create_new_interpolator
        create_new_interpolator = dict(
            velocity_dispersion=dict(
                create_new=False, resolution=200, zl_resolution=48, cdf_size=400
            ),
            axis_ratio=dict(
                create_new=False, resolution=200, sigma_resolution=48, cdf_size=400
            ),
            lens_redshift_sl=dict(
                create_new=False, resolution=16, zl_resolution=48, cdf_size=400
            ),
            lens_redshift=dict(
                create_new=False, resolution=16, zl_resolution=48, cdf_size=400
            ),
            optical_depth=dict(create_new=False, resolution=16, cdf_size=500),
            comoving_distance=dict(create_new=False, resolution=500),
            angular_diameter_distance=dict(create_new=False, resolution=500),
            angular_diameter_distance_z1z2=dict(create_new=False, resolution=500),
            differential_comoving_volume=dict(create_new=False, resolution=500),
            # density_profile_slope=dict(create_new=False, resolution=100),
            zs_sl=dict(create_new=False, resolution=200, cdf_size=500),
        )
        # the following resolution elements are for the cross section
        # interpolator, [25, 25, 45, 15, 15] corresponds to the resolution for
        # e1, e2, gamma, gamma1, gamma2
        spacing_config = {
            "e1": {
                "mode": "two_sided_mixed_grid", 
                "power_law_part": "lower",
                "spacing_trend": "increasing",
                "power": 2.5,
                "value_transition_fraction": 0.7,
                "num_transition_fraction": 0.9,
                "auto_match_slope": True
            },
            "e2": {
                "mode": "two_sided_mixed_grid", 
                "power_law_part": "lower",
                "spacing_trend": "increasing",
                "power": 2.5,
                "value_transition_fraction": 0.7,
                "num_transition_fraction": 0.9,
                "auto_match_slope": True
            },
            "gamma": {
                "mode": "two_sided_mixed_grid", 
                "power_law_part": "lower",
                "spacing_trend": "increasing",
                "power": 2.5,
                "value_transition_fraction": 0.7,
                "num_transition_fraction": 0.9,
                "auto_match_slope": True
            },
            "gamma1": {
                "mode": "two_sided_mixed_grid", 
                "power_law_part": "lower",
                "spacing_trend": "increasing",
                "power": 2.5,
                "value_transition_fraction": 0.7,
                "num_transition_fraction": 0.9,
                "auto_match_slope": True
            },
            "gamma2": {
                "mode": "two_sided_mixed_grid", 
                "power_law_part": "lower",
                "spacing_trend": "increasing",
                "power": 2.5,
                "value_transition_fraction": 0.7,
                "num_transition_fraction": 0.9,
                "auto_match_slope": True
            },
        }

        if lens_type == "sis_galaxy":
            create_new_interpolator.update(
                cross_section=dict(create_new=False, resolution=[5, 5, 5, 5, 5], spacing_config=spacing_config),
            )
        elif lens_type == "sie_galaxy":
            create_new_interpolator.update(
                cross_section=dict(create_new=False, resolution=[25, 25, 5, 5, 5], spacing_config=spacing_config ),
            )
        elif lens_type == "epl_shear_galaxy":
            create_new_interpolator.update(
                cross_section=dict(create_new=False, resolution=[25, 25, 45, 15, 15], spacing_config=spacing_config),
            )

        if isinstance(user_create_new_interpolator, dict):
            for key, value in user_create_new_interpolator.items():
                if (
                    key in create_new_interpolator
                    and isinstance(create_new_interpolator[key], dict)
                    and isinstance(value, dict)
                ):
                    create_new_interpolator[key].update(value)
                else:
                    create_new_interpolator[key] = value
        # if create_new_interpolator is True, create new interpolator for all
        elif user_create_new_interpolator:
            for key in create_new_interpolator.keys():
                create_new_interpolator[key]["create_new"] = True

        # update the create_new_interpolator
        try:
            if self.create_new_interpolator is None:
                self.create_new_interpolator = create_new_interpolator
            else:
                self.create_new_interpolator.update(create_new_interpolator)
        except AttributeError:
            self.create_new_interpolator = create_new_interpolator

    def _lens_functions_and_sampler_categorization(
        self,
        lens_type,
        lens_priors,
        lens_priors_params,
        lens_functions,
        lens_functions_params,
    ):
        """
        Helper to categorize, set and update lens functions and samplers.

        There are three levels.
        1. Set default lens functions and samplers based on lens type. \n
        2. If there is user input, update the lens functions and samplers with user input. \n
        3. Check the validity of user input and update the lens functions and samplers. User input can be different from the level (1). It can either be one of the internal functions or a user-defined function. If it is one of the internal functions, the parameters will be updated (missing parameters will be taken care of) with user input if given, otherwise it will use the default parameters. If it is a user-defined function, it will be directly used without any parameter update. \n

        Parameters
        ----------
        lens_priors : ``dict`` or ``None``
            User-provided sampler names or sampler instances. For custom samplers,
            pass an instance that provides both ``.rvs(...)`` and ``.pdf(...)``
            (e.g. a ``ler.utils.FunctionConditioning`` object). Plain callable
            sampler functions are not supported.
        lens_priors_params : ``dict`` or ``None``
            User-provided sampler parameters.
            Nested dicts merge into existing defaults key-by-key (outer keys unchanged), so overriding e.g.\ ``cross_section_epl_shear_interpolation`` does not drop other fields for that sampler.
        lens_functions : ``dict`` or ``None``
            User-provided lens function names or functions.
        lens_functions_params : ``dict`` or ``None``
            User-provided lens function parameters.
        """

        # 1. Set default lens functions and samplers based on lens type.

        (
            self.lens_priors,
            self.lens_priors_params,
            self.lens_functions,
            self.lens_functions_params,
        ) = self._default_lens_priors_and_functions(lens_type)

        # 2. If there is user input, update the lens functions and samplers with user input.

        # update the priors if input is given
        if lens_priors:
            self.lens_priors.update(lens_priors)
        if lens_priors_params:
            # Nested merge per sampler category: shallow dict.update replaces whole
            # inner dicts and would silently drop unrelated keys unless the caller
            # repeats every default field (easy to miss ``cross_section_epl_shear_*``).
            for key, value in lens_priors_params.items():
                if isinstance(value, dict) and isinstance(
                    self.lens_priors_params.get(key), dict
                ):
                    self.lens_priors_params[key].update(value)
                else:
                    self.lens_priors_params[key] = value
        if lens_functions:
            self.lens_functions.update(lens_functions)
        if lens_functions_params:
            self.lens_functions_params.update(lens_functions_params)

        # 3. Check the validity of user input and update the lens functions and samplers.

        # EPL + Shear default maginf=-100.0
        if self.lens_functions['cross_section'] == 'cross_section_epl_shear_interpolation':
            if self.lens_type == 'sis_galaxy':
                self.lens_functions_params["cross_section"].update(dict(maginf=-10000.0)) # maginf for SIS is set wrt to analytic cross section
            elif self.lens_type == 'sie_galaxy':
                self.lens_functions_params['cross_section'].update(dict(maginf=-10000.0)) # maginf for SIE is set wrt to analytic cross section

        sampler_prior_names = [
            "zs_sl",
            "lens_redshift_sl",
            "lens_redshift",
            "velocity_dispersion",
            "axis_ratio",
            "axis_rotation_angle",
            "external_shear1",
            "external_shear2",
            "density_profile_slope",
        ]

        def _has_rvs_and_pdf(obj):
            return (
                obj is not None
                and hasattr(obj, "rvs")
                and callable(getattr(obj, "rvs"))
                and hasattr(obj, "pdf")
                and callable(getattr(obj, "pdf"))
            )

        # Reconcile sampler parameter dicts: ``available_lens_priors`` defaults merged with ``self``.
        # ``lens_priors`` selects optional overrides; callers may pass only nested ``lens_priors_params``.
        for name in sampler_prior_names:  # e.g. name='axis_ratio'
            if lens_priors is not None and name in lens_priors:
                sampler_candidate = lens_priors[name]
            else:
                sampler_candidate = self.lens_priors.get(name)

            # Custom sampler object: no template lookup from ``available_lens_priors``.
            if isinstance(sampler_candidate, FunctionConditioning) or _has_rvs_and_pdf(
                sampler_candidate
            ):
                continue

            if sampler_candidate is None:
                continue

            if not isinstance(sampler_candidate, str):
                raise ValueError(
                    f"Given sampler for `{name}` must be (a) a string name from "
                    "``available_lens_priors``, or (b) a ``FunctionConditioning`` "
                    "(or comparable) providing ``.rvs`` and ``.pdf``."
                )

            sampler_type = sampler_candidate
            dict_ = self.available_lens_priors[
                name
            ]  # e.g. {'axis_ratio_padilla_strauss': {'q_min': 0.2, 'q_max': 1.0}, ....}
            if sampler_type in dict_:  # e.g. 'axis_ratio_padilla_strauss'
                template = dict_[sampler_type]
                param_dict = template.copy() if isinstance(template, dict) else template
                stored = self.lens_priors_params.get(name)

                # Merge template defaults with values already accumulated on ``self``
                # (nested ``lens_priors_params`` update in step 2).
                if isinstance(param_dict, dict):
                    if isinstance(stored, dict):
                        param_dict.update(stored)
                    self.lens_priors_params[name] = param_dict
                elif stored is not None:
                    self.lens_priors_params[name] = stored
                else:
                    self.lens_priors_params[name] = param_dict
            else:
                raise ValueError(
                    f"{name} sampler {sampler_type} not available.\n Available {name} samplers and its parameters are: {dict_[name]}"
                )

        lens_function_names = ["cross_section_based_sampler", "optical_depth", "cross_section"]

        # if there is user input function, update the sampler priors
        for name in lens_function_names:
            if (lens_functions is not None) and (name in lens_functions):
                function_name = lens_functions[name]
                if isinstance(function_name, str):
                    # available lens functions for name e.g. 'optical_depth'
                    dict_ = self.available_lens_functions[name]
                    if function_name in dict_:
                        param_dict = dict_[function_name]
                        param_dict = param_dict.copy() if isinstance(param_dict, dict) else param_dict
                        if (lens_functions_params is None) or (
                            lens_functions_params.get(name) is None
                        ):  # not a dictionary
                            self.lens_functions_params[name] = param_dict
                        else:  # if there is user input lens_functions_params
                            param_dict.update(lens_functions_params[name])
                            self.lens_functions_params[name] = (
                                param_dict
                            )  # user inputs override default values
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
        self,
        size=1000,
        zs=None,
        get_attribute=False,
        **kwargs,
    ):
        """
        Sample lens redshifts conditioned on strong lensing (numerical method).

        This method computes the lens redshift distribution by numerically
        integrating over the velocity-dispersion number density, lensing
        cross-section, and differential comoving volume:

        .. math::

            \\frac{d\\tau}{dz_l}(z_s) \\propto
            \\frac{dV_c}{dz_l}\\int d\\sigma\\,
            n(\\sigma, z_l)\\,\\sigma_{\\rm SL}(\\sigma, z_l, z_s).

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
            Additional parameters forwarded from ``lens_priors_params`` (typically merged onto a fresh
            copy of the sampler defaults listed in ``available_lens_priors`` so shared templates are
            never mutated).
        cross_section_epl_shear_interpolation : ``bool``
            If True, uses :meth:`~cross_section_epl_shear_interpolation` as the cross-section
            function (with a lens-type-dependent interpolation grid) during the numerical
            integration for lens redshift distribution.
            default: False

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

        identifier_dict = dict(
            param_name = "lens_redshift_sl",
            sampler_type = "lens_redshift_strongly_lensed_numerical",
            lens_type = self.lens_type,
        )
        identifier_dict["resolution"] = self.create_new_interpolator[identifier_dict["param_name"]][
            "resolution"
        ]
        identifier_dict["zl_resolution"] = self.create_new_interpolator[
            identifier_dict["param_name"]
        ]["zl_resolution"]
        identifier_dict["cdf_size"] = self.create_new_interpolator[
            identifier_dict["param_name"]
        ]["cdf_size"]
        identifier_dict["integration_size"] = self.lens_priors_params[
            identifier_dict["param_name"]
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
            "external_shear1",
            "external_shear2",
        ]

        # cross_section cannot be FunctionConditioning object
        identifier_dict["cross_section"] = {}
        identifier_dict["cross_section"]["name"] = self.lens_functions["cross_section"]
        if self.lens_functions_params["cross_section"] is not None:
            for key, value in self.lens_functions_params["cross_section"].items():
                identifier_dict["cross_section"][key] = value
        identifier_dict["cross_section"]["resolution"] = self.create_new_interpolator[
            "cross_section"
        ]["resolution"]

        # All other parameters can be FunctionConditioning objects
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

        identifier_dict["external_shear1"] = {}
        identifier_dict["external_shear2"] = {}
        if self.external_shear1.info is not None:  # if external_shear1 is not None
            for key, value in self.external_shear1.info.items():
                if key not in name_list:
                    identifier_dict["external_shear1"][key] = value
        if self.external_shear2.info is not None:  # if external_shear2 is not None
            for key, value in self.external_shear2.info.items():
                if key not in name_list:
                    identifier_dict["external_shear2"][key] = value

        # Do not mutate shared entries inside ``available_lens_priors``: copy then merge kwargs.
        template = self.available_lens_priors[identifier_dict["param_name"]][
            identifier_dict["sampler_type"]
        ]
        if isinstance(template, dict):
            if template:
                merged_sampler_params = dict(template)
                merged_sampler_params.update(kwargs)
            else:
                merged_sampler_params = dict(kwargs)
        elif template is None:
            merged_sampler_params = dict(kwargs)
        else:
            merged_sampler_params = dict(kwargs)

        identifier_dict.update(merged_sampler_params)

        print("Numerically solving the lens redshift distribution...")
        zs_resolution = identifier_dict["resolution"]
        zs_min = self.z_min + 0.001 if self.z_min == 0.0 else self.z_min
        zs_max = self.z_max
        zs_array = generate_mixed_grid(zs_min, zs_max, zs_resolution)

        _, it_exist = interpolator_json_path(
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            interpolator_name=identifier_dict["sampler_type"],
        )

        # zl dependent velocity dispersion distribution
        zl_scaled = []
        zl_min = 0.0001
        for i, zs_ in enumerate(zs_array):
            zl_resolution = identifier_dict["zl_resolution"]
            buffer_ = np.linspace(zl_min, zs_ - zl_min, zl_resolution)
            zl_scaled.append(buffer_ / zs_)
        zl_scaled = np.array(zl_scaled)

        create_new = self.create_new_interpolator[identifier_dict["param_name"]]["create_new"]
        if not it_exist or create_new:
            number_density = self._helper_number_density_calculation(
                zl_scaled,
                zs_array,
                cross_section_epl_shear_interpolation=identifier_dict.get(
                    "cross_section_epl_shear_interpolation", False
                ),
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
            non_negative_function=True,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["sampler_type"],
            create_new=create_new,
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
            cdf_size=identifier_dict["cdf_size"],
        )

        # un-scaled lens redshift
        cdf_values = zl_object.cdf_values
        x_array = zl_object.x_array
        x_array_cdf = zl_object.x_array_cdf
        y_array = zl_object.conditioned_y_array
        function_spline = zl_object.function_spline
        pdf_norm_const = zl_object.pdf_norm_const

        @njit
        def zl_function(x, y):
            result = cubic_spline_interpolator2d(
                x / y, y, function_spline, x_array, y_array
            )
            # The differential cross-section is non-negative. Clamp tiny
            # negative values introduced by spline interpolation/roundoff.
            idx = result < 0.0
            result[idx] = 0.0
            return result

        zl_object.function = zl_function

        @njit
        def zl_pdf(x, y):
            result = (
                pdf_cubic_spline_interpolator2d(
                    x / y, y, pdf_norm_const, function_spline, x_array, y_array
                )
                / y
            )
            # PDF is non-negative; clamp tiny interpolation artifacts.
            idx = result < 0.0
            result[idx] = 0.0
            return result

        zl_object.pdf = zl_pdf

        @njit
        def zl_rvs(size, y):
            return (
                inverse_transform_sampler2d(size, y, cdf_values, x_array_cdf, y_array) * y
            )

        zl_object.rvs = zl_rvs

        return zl_object if get_attribute else zl_object(size, zs)

    def _helper_number_density_calculation(
        self, zl_scaled, zs, cross_section_epl_shear_interpolation=False
    ):
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
        shear1_rvs = self.external_shear1.rvs
        shear2_rvs = self.external_shear2.rvs
        # sigma_pdf_ = self.velocity_dispersion.pdf
        number_density_ = self.velocity_dispersion.function
        if cross_section_epl_shear_interpolation:
            # Choose interpolation grid resolution based on lens model type
            if self.lens_type == "sis_galaxy":
                size_list = [5, 5, 5, 5, 5]
            elif self.lens_type == "sie_galaxy":
                size_list = [25, 25, 5, 5, 5]
            else:  # "epl_shear_galaxy" (default)
                size_list = [25, 25, 45, 15, 15]

            cross_section_ = self.cross_section_epl_shear_interpolation(
                get_attribute=True,
                size_list=size_list,
            )
        else:
            cross_section_ = self.cross_section
        dVcdz_function = self.differential_comoving_volume.function

        from .sampler_functions import _njit_checks

        use_njit_sampler, sampler_dict = _njit_checks(
            sigma_rvs=sigma_rvs_,
            q_rvs=q_rvs_,
            phi_rvs=phi_rvs,
            gamma_rvs=gamma_rvs,
            shear1_rvs=shear1_rvs,
            shear2_rvs=shear2_rvs,
            # sigma_pdf=sigma_pdf_,
            number_density_function=number_density_,
            cross_section_function=cross_section_,
            dVcdz_function=dVcdz_function,
            create_njit_sampler=True,
        )

        use_multiprocessing = self.lens_priors_params["lens_redshift_sl"][
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
        shear1_rvs = sampler_dict["shear1_rvs"]
        shear2_rvs = sampler_dict["shear2_rvs"]
        number_density_function = sampler_dict["number_density_function"]
        cross_section_function = sampler_dict["cross_section_function"]

        sigma_min = self.lens_priors_params["velocity_dispersion"]["sigma_min"]
        sigma_max = self.lens_priors_params["velocity_dispersion"]["sigma_max"]

        # integration_size
        integration_size = (
            self.lens_priors_params["lens_redshift_sl"]["integration_size"]
            if "integration_size" in self.lens_priors_params["lens_redshift_sl"]
            else 20000
        )

        # dVcdz_function
        dVcdz_function = self.differential_comoving_volume.function

        from .mp import lens_redshift_strongly_lensed_njit
        from numba import set_num_threads

        # Set Numba threads to npool before entering prange
        set_num_threads(self.npool)

        density_array = lens_redshift_strongly_lensed_njit(
            zs,  # 1D
            zl_scaled,  # 2D
            sigma_min,
            sigma_max,
            q_rvs,
            phi_rvs,
            gamma_rvs,
            shear1_rvs,
            shear2_rvs,
            number_density_function,
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

        number_density = sampler_dict["number_density_function"]
        q_rvs = sampler_dict["q_rvs"]
        phi_rvs = sampler_dict["phi_rvs"]
        gamma_rvs = sampler_dict["gamma_rvs"]
        shear1_rvs = sampler_dict["shear1_rvs"]
        shear2_rvs = sampler_dict["shear2_rvs"]
        cross_section_function = sampler_dict["cross_section_function"]

        sigma_min = self.velocity_dispersion.info["sigma_min"]
        sigma_max = self.velocity_dispersion.info["sigma_max"]

        dVcdz_function = self.differential_comoving_volume.function

        integration_size = (
            self.lens_priors_params["lens_redshift_sl"]["integration_size"]
            if "integration_size" in self.lens_priors_params["lens_redshift_sl"]
            else 20000
        )

        # setup input params for multiprocessing (zs, zl_scaled, worker_index)
        input_params = np.array(
            [
                (zs[i], zl_scaled[i], i)  # worker index for result ordering
                for i in range(len(zs))
            ],
            dtype=object,
        )

        print("Computing lens redshift distribution with multiprocessing...")
        from .mp import lens_redshift_strongly_lensed_mp, _init_lens_redshift_worker

        density_array = np.zeros_like(zl_scaled)

        with Pool(
            processes=self.npool,
            initializer=_init_lens_redshift_worker,
            initargs=(
                sigma_min,
                sigma_max,
                number_density,
                q_rvs,
                phi_rvs,
                gamma_rvs,
                shear1_rvs,
                shear2_rvs,
                dVcdz_function,
                cross_section_function,
                integration_size,
            ),
        ) as pool:
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

        return density_array

    def lens_redshift_intrinsic_numerical(
        self, size=1000, zs=None, get_attribute=False, **kwargs
    ):
        """
        Sample intrinsic lens redshifts (numerical method).

         This method computes the intrinsic lens redshift distribution by
         numerically integrating the velocity-dispersion number density and
         differential comoving volume.

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

        identifier_dict = dict(
            param_name="lens_redshift",
            sampler_type="lens_redshift_intrinsic_numerical",
        )
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["resolution"] = self.create_new_interpolator[
            identifier_dict["param_name"]
        ]["resolution"]
        identifier_dict["cdf_size"] = self.create_new_interpolator[
            identifier_dict["param_name"]
        ]["cdf_size"]
        create_new = self.create_new_interpolator[identifier_dict["param_name"]][
            "create_new"
        ]

        name_list = ["z_min", "z_max", "cosmology", "velocity_dispersion"]

        identifier_dict["velocity_dispersion"] = {}
        if (
            self.velocity_dispersion.info is not None
        ):  # if velocity_dispersion is not None
            for key, value in self.velocity_dispersion.info.items():
                if key not in name_list:
                    identifier_dict["velocity_dispersion"][key] = value

        _, it_exist = interpolator_json_path(
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict['param_name'],
            interpolator_name=identifier_dict["sampler_type"],
        )

        if not it_exist or create_new:
            # calculate for zs=zs_max
            z_min = 0.0001 if self.z_min == 0 else self.z_min
            z_max = self.z_max
            resolution = identifier_dict["resolution"]
            zl_array = generate_mixed_grid(z_min, z_max, resolution)

            # create number density
            if self.velocity_dispersion.conditioned_y_array is None:

                def integrand(sigma, z):
                    return (
                        self.velocity_dispersion.function(np.array([sigma]))[0]
                        * self.differential_comoving_volume.function(np.array([z]))[0]
                    )

            else:

                def integrand(sigma, z):
                    return (
                        self.velocity_dispersion.function(
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
            # create spline fit for lens redshift, with zs=zs_max
            number_density_spline = CubicSpline(zl_array, integral)

            # zl dependent velocity dispersion distribution
            zl_scaled = []
            number_density = []
            zl_min = 0.0001
            zl_size = identifier_dict["resolution"] = self.create_new_interpolator[identifier_dict["param_name"]]["zl_resolution"]

            for i, zl_ in enumerate(zl_array):
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
            non_negative_function=True,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict['param_name'],
            name=identifier_dict["sampler_type"],
            create_new=create_new,
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
            cdf_size=identifier_dict["cdf_size"],
        )

        # un-scaled lens redshift
        cdf_values = zl_object.cdf_values
        x_array = zl_object.x_array
        x_array_cdf = zl_object.x_array_cdf
        y_array = zl_object.conditioned_y_array
        function_spline = zl_object.function_spline
        pdf_norm_const = zl_object.pdf_norm_const

        @njit
        def zl_function(x, y):
            return cubic_spline_interpolator2d(
                x / y, y, function_spline, x_array, y_array
            )

        zl_object.function = zl_function

        @njit
        def zl_pdf(x, y):
            return (
                pdf_cubic_spline_interpolator2d(
                    x / y, y, pdf_norm_const, function_spline, x_array, y_array
                )
                / y
            )

        zl_object.pdf = zl_pdf

        @njit
        def zl_rvs(size, y):
            return (
                inverse_transform_sampler2d(size, y, cdf_values, x_array_cdf, y_array) * y
            )

        zl_object.rvs = zl_rvs

        return zl_object if get_attribute else zl_object(size, zs)

    def lens_redshift_strongly_lensed_sis_analytical(
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
        # identifier_dict = {"name": "lens_redshift_strongly_lensed_sis_analytical"}
        identifier_dict = dict(
            param_name="lens_redshift_sl",
            sampler_type="lens_redshift_strongly_lensed_sis_analytical",
            lens_type=self.lens_type,
        )
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
        identifier_dict["resolution"] = self.create_new_interpolator[identifier_dict["param_name"]][
            "resolution"
        ]
        identifier_dict["cdf_size"] = self.create_new_interpolator[
            identifier_dict["param_name"]
        ]["cdf_size"]
        param_dict = self.available_lens_priors[identifier_dict["param_name"]][
            "lens_redshift_strongly_lensed_sis_analytical"
        ]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        zl_resolution = identifier_dict["resolution"]
        x_array = np.linspace(0.0, 1.0, zl_resolution)

        def pdf_(x):
            return 30 * x**2 * (1 - x) ** 2  # Haris et al. 2018 (A7)

        zl_object = FunctionConditioning(
            function=pdf_,
            x_array=x_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["sampler_type"],
            create_new=self.create_new_interpolator[identifier_dict["param_name"]]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
            cdf_size=identifier_dict["cdf_size"],
        )

        cdf_values_zl = zl_object.cdf_values
        x_array_cdf_zl = zl_object.x_array_cdf
        function_spline_zl = zl_object.function_spline
        # pdf_norm_const = zl_object.pdf_norm_const
        x_array_zl = zl_object.x_array

        spline_Dc = self.comoving_distance.function_spline
        x_array_Dc = self.comoving_distance.x_array
        inverse_spline_Dc = self.comoving_distance.function_inverse_spline
        z_array_Dc = self.comoving_distance.z_array

        @njit
        def zl_rvs(size, zs):
            r = inverse_transform_sampler(size, cdf_values_zl, x_array_cdf_zl)
            zs_Dc = cubic_spline_interpolator(zs, spline_Dc, x_array_Dc)
            zl_Dc = zs_Dc * r
            return cubic_spline_interpolator(zl_Dc, inverse_spline_Dc, z_array_Dc)

        @njit
        def zl_pdf(zl, zs):
            zs_Dc = cubic_spline_interpolator(zs, spline_Dc, x_array_Dc)
            zl_Dc = cubic_spline_interpolator(zl, spline_Dc, x_array_Dc)
            r = zl_Dc / zs_Dc

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
    def gengamma(self, size, get_attribute=False, **kwargs):
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
            | param_name        |                   | Name of the parameter                 |
            +-------------------+-------------------+---------------------------------------+
            | sampler_type      |                   | Type of the sample                    |
            +-------------------+-------------------+---------------------------------------+
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
        >>> od = OpticalDepth(lens_priors=dict(
        ...     velocity_dispersion="gengamma"))
        >>> sigma = od.velocity_dispersion(size=100)
        """

        identifier_dict = {
            "param_name": "velocity_dispersion",
            "sampler_type": "gengamma",
        }

        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["resolution"] = self.create_new_interpolator[
            "velocity_dispersion"
        ]["resolution"]
        identifier_dict["cdf_size"] = self.create_new_interpolator[
            "velocity_dispersion"
        ]["cdf_size"]
        param_dict = self.available_lens_priors["velocity_dispersion"][
            "gengamma"
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

        from .sampler_functions import gengamma_function

        def density_func_(sigma_):
            return gengamma_function(
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
            non_negative_function=True,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["sampler_type"],
            create_new=self.create_new_interpolator["velocity_dispersion"][
                "create_new"
            ],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
            cdf_size=identifier_dict["cdf_size"],
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
            | param_name        |                   | Name of the parameter                 |
            +-------------------+-------------------+---------------------------------------+
            | sampler_type      |                   | Type of the sample                    |
            +-------------------+-------------------+---------------------------------------+
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
        >>> od = OpticalDepth(lens_priors=dict(
        ...     velocity_dispersion="velocity_dispersion_bernardi"))
        >>> sigma = od.velocity_dispersion(size=100)
        """

        identifier_dict = {
            "param_name": "velocity_dispersion",
            "sampler_type": "velocity_dispersion_bernardi",
        }
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["resolution"] = self.create_new_interpolator[
            "velocity_dispersion"
        ]["resolution"]
        identifier_dict["cdf_size"] = self.create_new_interpolator[
            "velocity_dispersion"
        ]["cdf_size"]
        param_dict = self.available_lens_priors["velocity_dispersion"][
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

        def number_density_function(sigma):
            return velocity_dispersion_bernardi_denisty_function(
                sigma,
                alpha=self.lens_priors_params["velocity_dispersion"]["alpha"],
                beta=self.lens_priors_params["velocity_dispersion"]["beta"],
                phistar=self.lens_priors_params["velocity_dispersion"][
                    "phistar"
                ],
                sigmastar=self.lens_priors_params["velocity_dispersion"][
                    "sigmastar"
                ],
            )

        sigma_object = FunctionConditioning(
            function=number_density_function,
            x_array=sigma_array,
            non_negative_function=True,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["sampler_type"],
            create_new=self.create_new_interpolator["velocity_dispersion"][
                "create_new"
            ],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
            cdf_size=identifier_dict["cdf_size"],
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
            | param_name        |                   | Name of the parameter                 |
            +-------------------+-------------------+---------------------------------------+
            | sampler_type      |                   | Type of the sample                    |
            +-------------------+-------------------+---------------------------------------+
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

        identifier_dict = {"param_name": "velocity_dispersion", "sampler_type": "velocity_dispersion_ewoud"}
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["resolution"] = self.create_new_interpolator[
            "velocity_dispersion"
        ]["resolution"]
        identifier_dict["zl_resolution"] = self.create_new_interpolator[
            "velocity_dispersion"
        ]["zl_resolution"]
        identifier_dict["cdf_size"] = self.create_new_interpolator[
            "velocity_dispersion"
        ]["cdf_size"]
        param_dict = self.available_lens_priors["velocity_dispersion"][
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

        def number_density_function(sigma, zl):
            return velocity_dispersion_ewoud_denisty_function(
                sigma,
                zl,
                alpha=self.lens_priors_params["velocity_dispersion"]["alpha"],
                beta=self.lens_priors_params["velocity_dispersion"]["beta"],
                phistar=self.lens_priors_params["velocity_dispersion"][
                    "phistar"
                ],
                sigmastar=self.lens_priors_params["velocity_dispersion"][
                    "sigmastar"
                ],
            )

        z_min = self.z_min + 0.001 if self.z_min == 0.0 else self.z_min
        z_max = self.z_max
        z_resolution = identifier_dict["zl_resolution"]
        zl_array = generate_mixed_grid(z_min, z_max, z_resolution)

        sigma_object = FunctionConditioning(
            function=number_density_function,
            x_array=sigma_array,
            non_negative_function=True,
            conditioned_y_array=zl_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["sampler_type"],
            create_new=self.create_new_interpolator["velocity_dispersion"][
                "create_new"
            ],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
            cdf_size=identifier_dict["cdf_size"],
        )

        return sigma_object if get_attribute else sigma_object.rvs(size, zl)

    # ------------------------------
    # Axis ratio sampler functions
    # ------------------------------
    def rayleigh(self, size, sigma, get_attribute=False, **kwargs):
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
            Additional parameters (param_name, sampler_type, q_min, q_max).

        Returns
        -------
        q : ``numpy.ndarray`` or ``FunctionConditioning``
            Axis ratio samples or sampler object.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_priors=dict(axis_ratio="rayleigh"))
        >>> q = od.axis_ratio(size=100, sigma=np.ones(100)*200.)
        """

        identifier_dict = {"param_name": "axis_ratio", "sampler_type": "rayleigh"}

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
        identifier_dict["cdf_size"] = self.create_new_interpolator["axis_ratio"][
            "cdf_size"
        ]
        param_dict = self.available_lens_priors["axis_ratio"]["rayleigh"]
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

        from .sampler_functions import rayleigh_pdf

        def q_pdf(q, sigma):
            return rayleigh_pdf(
                q=q,
                sigma=sigma,
                q_min=identifier_dict["q_min"],
                q_max=identifier_dict["q_max"],
            )

        q_object = FunctionConditioning(
            function=q_pdf,
            x_array=q_array,
            non_negative_function=True,
            conditioned_y_array=sigma_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["sampler_type"],
            create_new=self.create_new_interpolator["axis_ratio"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
            cdf_size=identifier_dict["cdf_size"],
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
            Additional parameters (param_name, sampler_type, q_min, q_max).

        Returns
        -------
        q : ``numpy.ndarray`` or ``FunctionConditioning``
            Axis ratio samples or sampler object.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth(lens_priors=dict(axis_ratio="axis_ratio_padilla_strauss"))
        >>> q = od.axis_ratio(size=100)
        """

        identifier_dict = {"param_name": "axis_ratio", "sampler_type": "axis_ratio_padilla_strauss"}
        identifier_dict["resolution"] = self.create_new_interpolator["axis_ratio"][
            "resolution"
        ]
        param_dict = self.available_lens_priors["axis_ratio"][
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
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["sampler_type"],
            create_new=self.create_new_interpolator["axis_ratio"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
            cdf_size=identifier_dict["cdf_size"],
        )

        return q_object if get_attribute else q_object(size)

    def uniform(self, size, get_attribute=False, **kwargs):
        """
        Sample values from uniform distribution.

        Parameters
        ----------
        size : ``int``
            Number of samples to draw.
        get_attribute : ``bool``
            If True, return the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Model parameters: \n
            - x_min: Minimum value, default: 0.0 \n
            - x_max: Maximum value, default: 1.0

        Returns
        -------
        values : ``numpy.ndarray``
            Array of uniformly distributed values in range [x_min, x_max].
        """

        identifier_dict = dict(
            param_name = "unknown_parameter",
            sampler_type = "uniform",
        )
        param_dict = dict(x_min=0.0, x_max=1.0)
        param_dict.update(kwargs)
        identifier_dict.update(param_dict)

        x_min = identifier_dict["x_min"]
        x_max = identifier_dict["x_max"]

        @njit
        def pdf_(x):
            # 1.0 / (x_max - x_min) * np.ones(len(x))
            return np.where((x >= x_min) & (x <= x_max), 1.0 / (x_max - x_min), 0.0)

        @njit
        def rvs_(size):
            return np.random.uniform(x_min, x_max, size=size)

        object_ = FunctionConditioning(
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["sampler_type"],
            create_function_inverse=False,
            create_function=False,
            create_pdf=pdf_,
            create_rvs=rvs_,
            callback="rvs",
        )

        if get_attribute:
            return object_
        else:
            return object_.rvs(size)

    def constant_values_n_size(self, size=100, get_attribute=False, **kwargs):
        """
        Return array of constant values.

        Parameters
        ----------
        size : ``int``
            Number of values to return. \n
            default: 100
        get_attribute : ``bool``
            If True, return the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Model parameters: \n
            - value: Constant value to return, default: 0.0

        Returns
        -------
        values : ``numpy.ndarray``
            Array of constant values.
        """

        identifier_dict = dict(
            param_name = "unknown_parameter",
            sampler_type = "constant_values_n_size",
        )
        param_dict = dict(value=0.0)
        param_dict.update(kwargs)
        identifier_dict.update(param_dict)

        value = identifier_dict["value"]

        # pdf_, zero everywhere except at value
        @njit
        def pdf_(x):
            return np.where(x == value, 1.0, 0.0)

        # rvs_, return value
        @njit
        def rvs_(size):
            return np.full(size, value)

        object_ = FunctionConditioning(
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["sampler_type"],
            create_function_inverse=False,
            create_function=False,
            create_pdf=pdf_,
            create_rvs=rvs_,
            callback="rvs",
        )

        if get_attribute:
            return object_
        else:
            return object_.rvs(size)

    def normal(self, size, get_attribute=False, **kwargs):
        """
        Sample with a normal distribution.

        Parameters
        ----------
        size : ``int``
            Number of samples to draw.
        get_attribute : ``bool``
            If True, return the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Model parameters: \n
            - mu: Mean of the normal distribution, default: 1.99 \n
            - sigma: Standard deviation of the normal distribution, default: 0.149 \n

        Returns
        -------
        x : ``numpy.ndarray``
            Array of values drawn from the normal distribution.
        """

        # Get parameters
        identifier_dict = dict(
            param_name = "unknown_parameter",
            sampler_type = "normal",
        )
        param_dict = dict(mu=1.99, sigma=0.149)

        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        mu = identifier_dict["mu"]
        sigma = identifier_dict["sigma"]

        @njit
        def normal_pdf_(x):
            return (1.0 / (sigma * np.sqrt(2.0 * np.pi))) * np.exp(-0.5 * ((x - mu) / sigma) ** 2)

        @njit
        def normal_rvs_(size):
            return np.random.normal(mu, sigma, size)

        object_ = FunctionConditioning(
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["sampler_type"],
            create_function=False,
            create_pdf=normal_pdf_,
            create_rvs=normal_rvs_,
            callback='rvs',
        )

        return object_ if get_attribute else object_.rvs(size)

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
            If True, returns the function object instead of values. \n
            default: False
        **kwargs : ``dict``
            Additional parameters.

        Returns
        -------
        tau : ``numpy.ndarray`` or ``FunctionConditioning``
            Optical depth values or function object.

        Notes
        -----
        The optical depth is a cumulative integral of a differential cross-section
        and is therefore non-negative. When optical depth values are cached via
        spline interpolation, tiny negative values can appear near boundaries due
        to numerical interpolation artifacts. This routine clamps such values to 0
        by constructing the cached interpolator with ``non_negative_function=True``.
        """
        identifier_dict = dict(
            param_name="optical_depth",
            function_type="optical_depth_numerical",
        )

        identifier_dict["resolution"] = self.create_new_interpolator[identifier_dict["param_name"]][
            "resolution"
        ]
        identifier_dict["cdf_size"] = self.create_new_interpolator[
            identifier_dict["param_name"]
        ]["cdf_size"]
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo

        name_list = [
            "z_min",
            "z_max",
            "cosmology",
        ]

        identifier_dict["lens_redshift_sl"] = {}
        if self.lens_redshift_sl.info is not None:
            for key, value in self.lens_redshift_sl.info.items():
                if key not in name_list:
                    identifier_dict["lens_redshift_sl"][key] = value

        param_dict = self.available_lens_functions[identifier_dict["param_name"]][
            identifier_dict["function_type"]
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
        zs_array = generate_mixed_grid(z_min, z_max, resolution)

        def tau(zs):
            # self.lens_redshift_sl.function gives differential cross-section
            def integrand(zl_, zs_):
                return self.lens_redshift_sl.function(np.array([zl_]), np.array([zs_]))[0]

            integral = [quad(integrand, 0.0, z, args=(z))[0] for z in zs]
            return integral

        tau_object = FunctionConditioning(
            function=tau,
            x_array=zs_array,
            conditioned_y_array=None,
            non_negative_function=True,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["function_type"],
            create_new=self.create_new_interpolator[identifier_dict["param_name"]]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="function",
            cdf_size=identifier_dict["cdf_size"],
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
        Ds = self.angular_diameter_distance.function(zs)
        Dls = self.angular_diameter_distance_z1z2.function(zl, zs)
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

        Notes
        -----
        The analytic SIS optical depth is non-negative. When cached via spline
        interpolation, tiny negative values can occur at the edges due to
        interpolation/roundoff. The cached interpolator is therefore created with
        ``non_negative_function=True`` to clamp such artifacts to 0.
        """

        identifier_dict = dict(
            param_name="optical_depth",
            function_type="optical_depth_sis_analytic",
        )
        identifier_dict["z_min"] = self.z_min if self.z_min > 0.0 else 0.0001
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["resolution"] = self.create_new_interpolator[identifier_dict["param_name"]][
            "resolution"
        ]
        identifier_dict["cdf_size"] = self.create_new_interpolator[
            identifier_dict["param_name"]
        ]["cdf_size"]
        param_dict = self.available_lens_functions[identifier_dict["param_name"]][
            identifier_dict["function_type"]
        ]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        z_min = identifier_dict["z_min"]
        z_max = identifier_dict["z_max"]
        z_resolution = identifier_dict["resolution"]
        zs_arr = generate_mixed_grid(z_min, z_max, z_resolution)

        def tau(zs):
            Dc = self.comoving_distance.function(zs)
            return 4.192e-15 * Dc**3

        tau_object = FunctionConditioning(
            function=tau,
            x_array=zs_arr,
            conditioned_y_array=None,
            non_negative_function=True,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["function_type"],
            create_new=self.create_new_interpolator[identifier_dict["param_name"]]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="function",
            cdf_size=identifier_dict["cdf_size"],
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
        >>> print(od.cross_section_sis(sigma=200., zl=0.5, zs=1.0))
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
        >>> print(od.cross_section_sie_feixu(sigma=200., zl=0.5, zs=1.0, q=1.0))
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

    def _two_sided_mixed_grid(
        self,
        x_min,
        x_max,
        resolution,
        power_law_part='lower',
        spacing_trend='increasing',
        power=2.3,
        value_transition_fraction=0.3,
        num_transition_fraction=0.6,
        auto_match_slope=True,
    ):
        # n1 chosen so that n1 + (n1-1) = 2*n1-1 >= resolution.
        # The -1 is because x2 drops the shared midpoint (0.5) via [1:].
        n1 = (resolution + 2) // 2

        # generate grid from 0 to 0.5
        x1 = generate_mixed_grid(
                x_min=0., 
                x_max=0.5, 
                resolution=n1,
                power_law_part=power_law_part,
                spacing_trend=spacing_trend,
                power=power,
                value_transition_fraction=value_transition_fraction,
                num_transition_fraction=num_transition_fraction,
                auto_match_slope=auto_match_slope
            )
        # Mirror x1 onto [0.5, 1.0]; [1:] removes the shared endpoint at 0.5.
        x2 = np.sort(1 - x1)[1:]
        u_grid = np.concatenate([x1, x2])[:resolution]
        x = x_min + u_grid * (x_max - x_min)
        
        return x

    def _parameter_spacing(self, x_min, x_max, resolution, config=None):
        """
        Create 1D grid axis based on spacing config.

        Supported modes:
        - ``linear`` (default)
        - ``geom_mixed`` (power-law + linear)
        - ``gaussian_mixed`` (Gaussian concentration + linear remainder)
        """
        if config is None:
            return np.linspace(x_min, x_max, resolution)

        mode = str(config.get("mode", "linear")).lower()

        if mode == "powerlaw_mixed":
            return generate_mixed_grid(
                x_min=x_min, 
                x_max=x_max, 
                resolution=resolution,
                power_law_part=config.get("power_law_part", "lower"),
                spacing_trend=config.get("spacing_trend", "increasing"),
                power=config.get("power", 2.5),
                value_transition_fraction=config.get("value_transition_fraction", 0.6),
                num_transition_fraction=config.get("num_transition_fraction", 0.3),
                auto_match_slope=config.get("auto_match_slope", True)
            )
        if mode == "two_sided_mixed_grid":
            return self._two_sided_mixed_grid(
                x_min=x_min, 
                x_max=x_max, 
                resolution=resolution,
                power_law_part=config.get("power_law_part", "lower"),
                spacing_trend=config.get("spacing_trend", "increasing"),
                power=config.get("power", 2.5),
                value_transition_fraction=config.get("value_transition_fraction", 0.6),
                num_transition_fraction=config.get("num_transition_fraction", 0.3),
                auto_match_slope=config.get("auto_match_slope", True)
            )

        return np.linspace(x_min, x_max, resolution)

    def _compute_parameter_bounds(self, size_list=[25, 25, 45, 15, 15]):
        """
        Sample lens parameters with a fixed internal seed and return min/max
        bounds used for grid construction.

        The global NumPy random state is saved before seeding and restored
        afterwards, so no external draws are affected.

        Parameters
        ----------
        size_list : list
            Grid axis sizes ``[n_e1, n_e2, n_gamma, n_gamma1, n_gamma2]``.
            Used to set the padding around the sampled extremes.

        Returns
        -------
        bounds : dict
            Keys: ``q_min``, ``q_max``, ``phi_min``, ``phi_max``,
            ``e1_min``, ``e1_max``, ``e2_min``, ``e2_max``,
            ``gamma_min``, ``gamma_max``,
            ``gamma1_min``, ``gamma1_max``,
            ``gamma2_min``, ``gamma2_max``.
        """

        size = 1000000

        _rng_state = np.random.get_state()
        _numba_threads = get_num_threads()
        try:
            np.random.seed(42)
            _seed_numba_rng(42)
            # Keep njit+prange samplers deterministic while computing bounds.
            set_num_threads(1)

            q_array = getattr(self.axis_ratio, "x_array", None)
            if q_array is not None:
                q_min = np.min(q_array)
                q_max = np.max(q_array)
            else:
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

            phi_info = getattr(self.axis_rotation_angle, "info", None)
            if isinstance(phi_info, dict):
                phi_min = phi_info.get("x_min", 0.0)
                phi_max = phi_info.get("x_max", 2.0 * np.pi)
            else:
                phi_min = 0.0
                phi_max = 2.0 * np.pi

            if phi_min == phi_max:
                phi = self.axis_rotation_angle.rvs(size)
                phi_min, phi_max = phi.min(), phi.max()

            del_phi = 2 * (phi_max - phi_min) / size_list[1]
            if del_phi == 0:
                phi_min = max(0.0, phi_min - 0.01)
                phi_max = min(2.0 * np.pi, phi_max + 0.01)

            # Build deterministic e1/e2 bounds from a fixed q-phi grid.
            n_q = max(256, int(size_list[0]) * 16)
            n_phi = max(256, int(size_list[1]) * 16)
            q_grid = np.linspace(q_min, q_max, n_q)
            phi_grid = np.linspace(phi_min, phi_max, n_phi)
            q_eval = np.repeat(q_grid, n_phi)
            phi_eval = np.tile(phi_grid, n_q)

            e1, e2 = phi_q2_ellipticity(phi_eval, q_eval)
            e1_min, e1_max = e1.min(), e1.max()
            e2_min, e2_max = e2.min(), e2.max()
            del_e1 = 2 * (e1_max - e1_min) / size_list[0]
            del_e2 = 2 * (e2_max - e2_min) / size_list[1]
            e1_min = max(-0.9999, e1_min - del_e1)
            e1_max = min(0.9999, e1_max + del_e1)
            e2_min = max(-0.9999, e2_min - del_e2)
            e2_max = min(0.9999, e2_max + del_e2)
            if del_e1 == 0:
                e1_min = max(-0.9999, e1_min - 0.01)
                e1_max = min(0.9999, e1_max + 0.01)
            if del_e2 == 0:
                e2_min = max(-0.9999, e2_min - 0.01)
                e2_max = min(0.9999, e2_max + 0.01)

            gamma = self.density_profile_slope.rvs(size)
            gamma_min, gamma_max = gamma.min(), gamma.max()
            del_gamma = 2 * (gamma_max - gamma_min) / size_list[2]
            gamma_min = gamma_min - del_gamma
            gamma_max = gamma_max + del_gamma
            if del_gamma == 0:
                gamma_min = gamma_min - 0.01
                gamma_max = gamma_max + 0.01

            gamma1 = self.external_shear1.rvs(size)
            gamma2 = self.external_shear2.rvs(size)

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
        finally:
            np.random.set_state(_rng_state)
            set_num_threads(_numba_threads)

        bounds = dict(
            q_min=q_min, q_max=q_max,
            phi_min=phi_min, phi_max=phi_max,
            e1_min=e1_min, e1_max=e1_max,
            e2_min=e2_min, e2_max=e2_max,
            gamma_min=gamma_min, gamma_max=gamma_max,
            gamma1_min=gamma1_min, gamma1_max=gamma1_max,
            gamma2_min=gamma2_min, gamma2_max=gamma2_max,
        )

        return {key: float(np.round(value, 4)) for key, value in bounds.items()}

    def create_parameter_grid(self, size_list=[25, 25, 45, 15, 15], spacing_config=None, bounds=None):
        """
        Create a parameter grid for lens galaxies.

        Parameters
        ----------
        size_list : list
            List of sizes for each parameter grid.
        spacing_config : dict or None
            Optional per-parameter spacing configuration. Keys can include
            ``q``, ``phi``, ``e1``, ``e2``, ``gamma``, ``gamma1``, ``gamma2`` and values are
            dictionaries, e.g. ``{"mode": "two_sided_mixed_grid", "power_law_part": "lower", "spacing_trend": "increasing", "power": 2.5, "value_transition_fraction": 0.6, "num_transition_fraction": 0.3, "auto_match_slope": True}``.
        bounds : dict or None
            Pre-computed parameter bounds from ``_compute_parameter_bounds``. If
            None, bounds are computed internally.

        Returns
        -------
        zl : numpy.ndarray
            Lens redshifts.
        sigma : numpy.ndarray
            Velocity dispersions.
        q : numpy.ndarray
            Axis ratios.
        """

        if bounds is None:
            bounds = self._compute_parameter_bounds(size_list)
        e1_min = bounds["e1_min"]
        e1_max = bounds["e1_max"]
        e2_min = bounds["e2_min"]
        e2_max = bounds["e2_max"]
        gamma_min = bounds["gamma_min"]
        gamma_max = bounds["gamma_max"]
        gamma1_min = bounds["gamma1_min"]
        gamma1_max = bounds["gamma1_max"]
        gamma2_min = bounds["gamma2_min"]
        gamma2_max = bounds["gamma2_max"]

        if spacing_config is None:
            spacing_config = {}

        ###########################
        # sampling the parameters #
        ###########################
        # q = self._parameter_spacing(
        #     q_min, q_max, size_list[0], spacing_config.get("q")
        # )
        # phi = self._parameter_spacing(
        #     phi_min, phi_max, size_list[1], spacing_config.get("phi")
        # )
        # e1, e2 = phi_q2_ellipticity(phi, q)

        # Build interpolation axes directly in e1/e2 space.
        # This avoids distortions introduced by elementwise (q, phi) -> (e1, e2)
        # pairing when nonuniform spacing is used.
        e1 = self._parameter_spacing(
            e1_min,
            e1_max,
            size_list[0],
            spacing_config.get("e1", spacing_config.get("e1")),
        )
        e2 = self._parameter_spacing(
            e2_min,
            e2_max,
            size_list[1],
            spacing_config.get("e2", spacing_config.get("e2")),
        )

        e1 = np.sort(e1)
        e2 = np.sort(e2)
        # print(f"e1 : {e1}")
        # print(f"e2 : {e2}")

        gamma = self._parameter_spacing(
            gamma_min, gamma_max, size_list[2], spacing_config.get("gamma")
        )

        gamma1 = self._parameter_spacing(
            gamma1_min, gamma1_max, size_list[3], spacing_config.get("gamma1")
        )
        gamma2 = self._parameter_spacing(
            gamma2_min, gamma2_max, size_list[4], spacing_config.get("gamma2")
        )

        # return q, phi, e1, e2, gamma, gamma1, gamma2
        return e1, e2, gamma, gamma1, gamma2

    def cross_section_epl_shear_interpolation_init(
        self, file_path, size_list, spacing_config=None, batch_size=50000, bounds=None
    ):
        print(f"Cross section interpolation data points will be created at {file_path}")

        # ----------------------------------
        # cross section unit to cross section ratio fitting
        # ----------------------------------

        # # --------------------------
        # # warm up if cross_section_epl_shear_numerical_mp is used
        # # --------------------------
        # size = 2 # small size for warm up
        # zs = np.random.uniform(self.z_min, self.z_max, size)
        # zl = np.zeros(size)
        # for i in range(size):
        #     zl[i] = np.random.uniform(0.001, zs[i] - 0.001)

        # sigma = np.random.uniform(
        #     self.lens_priors_params["velocity_dispersion"]["sigma_min"],
        #     self.lens_priors_params["velocity_dispersion"]["sigma_max"],
        #     size,
        # )

        # # Sample individual parameter values (not grid)
        # try:
        #     q = self.axis_ratio.rvs(size, sigma)
        # except TypeError:
        #     q = self.axis_ratio.rvs(size)
        # phi = self.axis_rotation_angle.rvs(size)
        # e1, e2 = phi_q2_ellipticity(phi, q)
        # gamma = self.density_profile_slope.rvs(size)
        # gamma1 = self.external_shear1.rvs(size)
        # gamma2 = self.external_shear2.rvs(size)
        # size = np.prod(size_list)
        # # --------------------------

        # cs data points for interpolation
        e1, e2, gamma, gamma1, gamma2 = self.create_parameter_grid(
            size_list=size_list, spacing_config=spacing_config, bounds=bounds
        )

        e1_arr, e2_arr, gamma_arr, gamma1_arr, gamma2_arr = np.meshgrid(
            e1, e2, gamma, gamma1, gamma2, indexing="ij"
        )
        e1_arr = e1_arr.flatten()
        e2_arr = e2_arr.flatten()
        gamma_arr = gamma_arr.flatten()
        gamma1_arr = gamma1_arr.flatten()
        gamma2_arr = gamma2_arr.flatten()

        # # ---------
        # cs_unit = self.cross_section_epl_shear_numerical_mp(
        #     theta_E=np.ones(size),
        #     gamma=gamma_arr,
        #     gamma1=gamma1_arr,
        #     gamma2=gamma2_arr,
        #     e1=e1_arr,
        #     e2=e2_arr,
        #     verbose=True,
        # )
        # # ---------

        # --------------------------
        # cross_section_epl_shear_unit njit
        # --------------------------
        from ..image_properties.cross_section_njit import cross_section_epl_shear_unit

        # Set Numba threads to npool before entering prange
        from numba import set_num_threads

        set_num_threads(self.npool)

        args = self.lens_functions_params["cross_section"]

        if args is not None and 'num_th' in args:
            num_th = args['num_th']
        else:
            num_th = 500
        if args is not None and 'maginf' in args:
            maginf = args['maginf']
        else:
            maginf = 100.0
        # print(f"\n num_th: {num_th}, maginf: {maginf} \n")

        # warm up the Numba function
        cross_section_epl_shear_unit(
            e1=e1_arr[:2],
            e2=e2_arr[:2],
            gamma=gamma_arr[:2],
            gamma1=gamma1_arr[:2],
            gamma2=gamma2_arr[:2],
            num_th=num_th,
            maginf=maginf,
        )

        """Memory-efficient wrapper - reuses the same small work arrays."""
        n = len(e1_arr)
        cs_unit = np.empty(n)

        print(
            f"Computing n={int(n)} cross-section data points for interpolation with ler.image_properties.cross_section_njit.cross_section_epl_shear_unit..."
        )
        
        for start in range(0, n, batch_size):
            end = min(start + batch_size, n)
            cs_unit[start:end] = cross_section_epl_shear_unit(
                e1=e1_arr[start:end], e2=e2_arr[start:end], gamma=gamma_arr[start:end],
                gamma1=gamma1_arr[start:end], gamma2=gamma2_arr[start:end],
                num_th=num_th, maginf=maginf
            )

        # --------------------------

        # # find slope and intercept for cs_unit to cs conversion
        # print("Finding the conversion factor from cross-section unit to actual cross-section...")
        # # choose random gamma, gamma1, gamma2 for testing cs_unit to cs conversion
        # size = min(10000, len(gamma_arr))
        # idx = np.random.choice(len(gamma_arr), size, replace=False)
        # theta_E = np.random.uniform(1.0e-10, 1.0e-8, size)
        # if self.lens_type in ["epl_shear_galaxy", "sie_galaxy"]:
        #     cs = cross_section_epl_shear_unit(
        #         e1=e1_arr[idx],
        #         e2=e2_arr[idx],
        #         gamma=gamma_arr[idx],
        #         gamma1=gamma1_arr[idx],
        #         gamma2=gamma2_arr[idx],
        #         theta_E=theta_E,
        #     )
        # elif self.lens_type == "sis_galaxy":
        #     cs = np.pi * theta_E * theta_E * np.ones(size)
        # else:
        #     raise ValueError(f"Unsupported lens type {self.lens_type} for cross-section interpolation.")
        # # cs = self.cross_section_epl_shear_numerical_mp(
        # #     theta_E=theta_E,
        # #     gamma=gamma_arr[idx],
        # #     gamma1=gamma1_arr[idx],
        # #     gamma2=gamma2_arr[idx],
        # #     e1=e1_arr[idx],
        # #     e2=e2_arr[idx],
        # #     verbose=False,
        # # )
        # y_ = cs / cs_unit[idx]
        # x_ = np.pi * theta_E * theta_E  # SIS cross-section as reference
        # valid = np.isfinite(y_) & np.isfinite(x_) & (y_ > 0.0) & (cs_unit[idx] > 0.0)

        # if np.count_nonzero(valid) >= 2:
        #     cs_csunit_slope, cs_csunit_intercept = np.polyfit(x_[valid], y_[valid], 1)
        # else:
        #     cs_csunit_slope = 0.31830988618379075
        #     cs_csunit_intercept = -3.2311742677852644e-27
        cs_csunit_slope = 0.31830988618379075
        cs_csunit_intercept = -3.2311742677852644e-27

        # # calculate relative error and adjust the offset
        # size = min(100, len(gamma_arr))
        # idx = np.random.choice(len(gamma_arr), size, replace=False)
        # theta_E = np.random.uniform(1.0e-12, 1.0e-5, size)

        # cs_pred_unit = cross_section_epl_shear_unit(
        #     e1=e1_arr[idx],
        #     e2=e2_arr[idx],
        #     gamma=gamma_arr[idx],
        #     gamma1=gamma1_arr[idx],
        #     gamma2=gamma2_arr[idx],
        # )
        # cs_sis = np.pi * theta_E * theta_E
        # cs_pred = cs_pred_unit * (cs_csunit_intercept + cs_csunit_slope * cs_sis)

        # if self.lens_type in ["epl_shear_galaxy", "sie_galaxy"]:
        #     cs_true = self.cross_section_epl_shear_numerical_mp(
        #         theta_E=theta_E,
        #         gamma=gamma_arr[idx],
        #         gamma1=gamma1_arr[idx],
        #         gamma2=gamma2_arr[idx],
        #         e1=e1_arr[idx],
        #         e2=e2_arr[idx],
        #         verbose=False,
        #     )
        # elif self.lens_type == "sis_galaxy":
        #     cs_true = np.pi * theta_E * theta_E * np.ones(size)

        # idx_nonzero = cs_true > 0.0
        # relative_error = np.median((cs_pred[idx_nonzero] - cs_true[idx_nonzero]) / cs_true[idx_nonzero])
        # # correct the intercept to minimize the median relative error
        # cs_csunit_slope = cs_csunit_slope / (1 + relative_error)
        # cs_csunit_intercept = cs_csunit_intercept / (1 + relative_error)

        # ----------
        cs_unit_grid = cs_unit.reshape(size_list)

        # save cs data
        save_json(
            file_path,
            [
                e1,
                e2,
                gamma,
                gamma1,
                gamma2,
                cs_unit_grid,
                cs_csunit_slope,
                cs_csunit_intercept,
            ],
        )

        return e1, e2, gamma, gamma1, gamma2, cs_unit_grid, cs_csunit_slope, cs_csunit_intercept

    def cross_section_epl_shear_interpolation(
        self,
        zs=None,
        zl=None,
        sigma=None,
        q=None,
        phi=None,
        gamma=None,
        gamma1=None,
        gamma2=None,
        get_attribute=False,
        size_list=[25, 25, 45, 15, 15],
        spacing_config=None,
        **kwargs,
    ):
        """
        Function to compute the cross-section correction factor
        """

        from .cross_section_interpolator import make_cross_section_area_reinit

        # spacing_config = kwargs.get("grid_spacing_config", None)
        if spacing_config is None:
            # spacing_config = {'e1': {'mode': 'mixed', 'geom_fraction': 0.6, 'transition_fraction': 0.3}, 'e2': {'mode': 'mixed', 'geom_fraction': 0.6, 'transition_fraction': 0.3}, 'gamma': {'mode': 'linear'}, 'gamma1': {'mode': 'linear'}, 'gamma2': {'mode': 'linear'}}
            spacing_config = self.create_new_interpolator["cross_section"]["spacing_config"]

        bounds = self._compute_parameter_bounds(size_list)

        param_dict_given = dict(
            z_min=self.z_min,
            z_max=self.z_max,
            resolution=self.create_new_interpolator["cross_section"]["resolution"],
            size_list=size_list,
            velocity_dispersion=dict(
                sigma_min=self.lens_priors_params["velocity_dispersion"]["sigma_min"],
                sigma_max=self.lens_priors_params["velocity_dispersion"]["sigma_max"],
            ),
            axis_ratio=dict(
                q_min=bounds["q_min"],
                q_max=bounds["q_max"],
            ),
            axis_rotation_angle=dict(
                phi_min=bounds["phi_min"],
                phi_max=bounds["phi_max"],
            ),
            density_profile_slope=dict(
                gamma_min=bounds["gamma_min"],
                gamma_max=bounds["gamma_max"],
            ),
            external_shear1=dict(
                gamma1_min=bounds["gamma1_min"],
                gamma1_max=bounds["gamma1_max"],
            ),
            external_shear2=dict(
                gamma2_min=bounds["gamma2_min"],
                gamma2_max=bounds["gamma2_max"],
            ),
            grid_spacing_config=spacing_config,
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
                cs_unit_grid,
                cs_csunit_slope,
                cs_csunit_intercept,
            ) = self.cross_section_epl_shear_interpolation_init(
                file_path, size_list, spacing_config=spacing_config, bounds=bounds
            )
        else:

            print(f"Cross section interpolation data points loaded from {file_path}")
            (
                e1_grid,
                e2_grid,
                gamma_grid,
                gamma1_grid,
                gamma2_grid,
                cs_unit_grid,
                cs_csunit_slope,
                cs_csunit_intercept,
            ) = load_json(file_path)

        self.cs_data_points_path = file_path

        # Create the interpolator instance
        cs_caculator = make_cross_section_area_reinit(
            e1_grid=np.array(e1_grid),
            e2_grid=np.array(e2_grid),
            gamma_grid=np.array(gamma_grid),
            gamma1_grid=np.array(gamma1_grid),
            gamma2_grid=np.array(gamma2_grid),
            cs_unit_grid=np.array(cs_unit_grid),
            Da_instance=self.angular_diameter_distance.function,
            csunit_to_cs_slope=cs_csunit_slope,
            csunit_to_cs_intercept=cs_csunit_intercept,
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

    # ----------------------------------------------------
    # EPL + external shear cross section using Numba njit
    # ----------------------------------------------------
    def cross_section_epl_shear_njit(
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
        num_th=500,
        maginf=-100.0,
        **kwargs,
    ):
        """
        Compute the lensing cross-section for EPL+shear lens model using Numba njit.

        Parameters
        ----------
        zs : ``float`` or ``numpy.ndarray``
            Source redshift(s).
        zl : ``float`` or ``numpy.ndarray``
            Lens redshift(s).
        sigma : ``float`` or ``numpy.ndarray``
            Velocity dispersion(s) in km/s.
        q : ``float`` or ``numpy.ndarray``
            Axis ratio(s) (b/a).
        phi : ``float`` or ``numpy.ndarray``
            Axis rotation angle(s) in radians.
        gamma : ``float`` or ``numpy.ndarray``
            Density profile slope(s).
        gamma1 : ``float`` or ``numpy.ndarray``
            External shear component 1.
        gamma2 : ``float`` or ``numpy.ndarray``
            External shear component 2.
        get_attribute : ``bool``, optional
            If True, return the interpolator instance instead of the cross-section.
            Default is False.
        num_th : ``int``, optional
            Number of theta values to use for the spline interpolation.
            Default is 500.
        maginf : ``float``, optional
            Magnitude limit for the cross-section calculation.
            Default is -100.0.
        **kwargs
            Additional keyword arguments to pass to the interpolator.

        Returns
        -------
        cs_caculator : ``ler.image_properties.cross_section_njit.CrossSectionNjit``
            The interpolator instance (if get_attribute=True).
        cross_section : ``float`` or ``numpy.ndarray``
            Lensing cross-section in arcsec^2.
        """

        from ..image_properties.cross_section_njit import make_cross_section_area_reinit

        # Create the interpolator instance
        Da_instance = self.angular_diameter_distance.function
        cs_caculator = make_cross_section_area_reinit(
            Da_instance=Da_instance,
            num_th=num_th,
            maginf=maginf,
        )

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

    # -------------------
    # properties
    # -------------------
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
        if prior in self.available_lens_priors["velocity_dispersion"]:
            print(f"using ler available velocity dispersion function : {prior}")
            args = self.lens_priors_params["velocity_dispersion"]
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
        elif (
            hasattr(prior, "rvs")
            and callable(getattr(prior, "rvs"))
            and hasattr(prior, "pdf")
            and callable(getattr(prior, "pdf"))
            and hasattr(prior, "function")
            and callable(getattr(prior, "function"))
        ):
            self._velocity_dispersion = prior
        else:
            raise ValueError(
                "velocity_dispersion must be either a sampler name from available_lens_priors['velocity_dispersion'] "
                "or an instance providing both `.rvs(...)` and `.pdf(...)` (e.g. `ler.utils.FunctionConditioning`)."
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
        if prior in self.available_lens_priors["axis_ratio"]:
            print(f"using ler available axis_ratio function : {prior}")
            args = self.lens_priors_params["axis_ratio"]
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
        elif (
            hasattr(prior, "rvs")
            and callable(getattr(prior, "rvs"))
            and hasattr(prior, "pdf")
            and callable(getattr(prior, "pdf"))
        ):
            self._axis_ratio = prior
        else:
            raise ValueError(
                "axis_ratio must be either a string in available_lens_priors['axis_ratio'] "
                "or an instance providing both `.rvs(...)` and `.pdf(...)` (e.g. `ler.utils.FunctionConditioning`)."
            )

    @property
    def zs_sl(self):
        """
        Strongly-lensed source redshift sampler object.

        Returns
        -------
        zs_sl : ``FunctionConditioning``
            Sampler object for source redshift conditioned on strong lensing.

        Notes
        -----
        The default strongly-lensed source redshift sampler,
        ``strongly_lensed_source_redshift``, is implemented in
        :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution`.
        This property exists in ``OpticalDepth`` so that derived classes can
        configure/override the prior using the common lens-prior mechanism.
        """
        return self._zs_sl

    @zs_sl.setter
    def zs_sl(self, prior):
        if prior in self.available_lens_priors["zs_sl"]:
            print(f"using ler available zs_sl function : {prior}")
            args = self.lens_priors_params["zs_sl"]
            if args is None:
                self._zs_sl = getattr(self, prior)(
                    size=None, get_attribute=True
                )
            else:
                self._zs_sl = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user provided custom zs_sl class/object of type ler.utils.FunctionConditioning"
            )
            self._zs_sl = prior
        else:
            raise ValueError(
                "zs_sl should be a string in available_lens_priors['zs_sl'] "
                "or an instance of `ler.utils.FunctionConditioning`."
            )

    @property
    def lens_redshift_sl(self):
        """
        Strongly lensed lens redshift sampler object.

        Returns a ``FunctionConditioning`` object with methods: \n
        - ``rvs(size, zs)``: Sample lens redshifts given source redshifts \n
        - ``pdf(zl, zs)``: Get probability density \n
        - ``function(zl, zs)``: Get effective lensing cross-section \n

        Returns
        -------
        lens_redshift_sl : ``FunctionConditioning``
            Sampler object for strongly lensed lens redshift.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> zl = od.lens_redshift_sl(size=100, zs=np.ones(100)*2.0)
        """

        return self._lens_redshift_sl

    @lens_redshift_sl.setter
    def lens_redshift_sl(self, prior):
        if prior in self.available_lens_priors["lens_redshift_sl"]:
            print(f"using ler available lens_redshift_sl function : {prior}")
            args = self.lens_priors_params["lens_redshift_sl"]
            if args is None:
                self._lens_redshift_sl = getattr(self, prior)(
                    size=None, zs=None, get_attribute=True
                )
            else:
                self._lens_redshift_sl = getattr(self, prior)(
                    size=None, zs=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user provided custom lens_redshift_sl class/object of type ler.utils.FunctionConditioning"
            )
            self._lens_redshift_sl = prior
        else:
            raise ValueError(
                "lens_redshift_sl should be string in available_lens_priors['lens_redshift_sl'] or class object of 'ler.utils.FunctionConditioning'"
            )

    @property
    def lens_redshift(self):
        """
        Lens redshift (intrinsic) sampler object.

        Returns a ``FunctionConditioning`` object with methods: \n
        - ``rvs(size)``: Sample lens redshifts \n
        - ``pdf(zl)``: Get probability density \n

        Returns
        -------
        lens_redshift : ``FunctionConditioning``
            Sampler object for lens redshift.

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> zl = od.lens_redshift(size=100)
        """

        return self._lens_redshift

    @lens_redshift.setter
    def lens_redshift(self, prior):
        if prior in self.available_lens_priors["lens_redshift"]:
            print(f"using ler available lens_redshift function : {prior}")
            args = self.lens_priors_params["lens_redshift"]
            if args is None:
                self._lens_redshift = getattr(self, prior)(
                    size=None, get_attribute=True
                )
            else:
                self._lens_redshift = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user provided custom lens_redshift class/object of type ler.utils.FunctionConditioning"
            )
            self._lens_redshift = prior
        elif (
            hasattr(prior, "rvs")
            and callable(getattr(prior, "rvs"))
            and hasattr(prior, "pdf")
            and callable(getattr(prior, "pdf"))
        ):
            self._lens_redshift = prior
        else:
            raise ValueError(
                "lens_redshift must be either a string in available_lens_priors['lens_redshift'] "
                "or an instance providing both `.rvs(...)` and `.pdf(...)` (e.g. `ler.utils.FunctionConditioning`)."
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
        if prior in self.available_lens_priors["axis_rotation_angle"]:
            print(f"using ler available axis_rotation_angle function : {prior}")
            args = self.lens_priors_params["axis_rotation_angle"]
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
        elif (
            hasattr(prior, "rvs")
            and callable(getattr(prior, "rvs"))
            and hasattr(prior, "pdf")
            and callable(getattr(prior, "pdf"))
        ):
            self._axis_rotation_angle = prior
        else:
            raise ValueError(
                "axis_rotation_angle must be either a string in available_lens_priors['axis_rotation_angle'] "
                "or an instance providing both `.rvs(...)` and `.pdf(...)` (e.g. `ler.utils.FunctionConditioning`)."
            )

    @property
    def external_shear1(self):
        """
        External shear sampler object (component 1).

        Returns a ``FunctionConditioning`` object with methods: \n
        - ``rvs(size)``: Sample shear component 1 (gamma1) \n
        - ``pdf(gamma1)``: Get probability density \n

        Returns
        -------
        external_shear1 : ``FunctionConditioning``
            Sampler object for external shear (component 1).

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> gamma1 = od.external_shear1(size=100)
        """

        return self._external_shear1

    @external_shear1.setter
    def external_shear1(self, prior):
        if prior in self.available_lens_priors["external_shear1"]:
            print(f"using ler available external_shear1 function : {prior}")
            args = self.lens_priors_params["external_shear1"]
            if args is None:
                self._external_shear1 = getattr(self, prior)(
                    size=None, get_attribute=True
                )
            else:
                self._external_shear1 = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user provided custom external_shear1 class/object of type ler.utils.FunctionConditioning"
            )
            self._external_shear1 = prior
        elif (
            hasattr(prior, "rvs")
            and callable(getattr(prior, "rvs"))
            and hasattr(prior, "pdf")
            and callable(getattr(prior, "pdf"))
        ):
            self._external_shear1 = prior
        else:
            raise ValueError(
                "external_shear1 must be either a string in available_lens_priors['external_shear1'] "
                "or an instance providing both `.rvs(...)` and `.pdf(...)` (e.g. `ler.utils.FunctionConditioning`)."
            )

    @property
    def external_shear2(self):
        """
        External shear sampler object (component 2).

        Returns a ``FunctionConditioning`` object with methods: \n
        - ``rvs(size)``: Sample shear component 2 (gamma2) \n
        - ``pdf(gamma2)``: Get probability density \n

        Returns
        -------
        external_shear2 : ``FunctionConditioning``
            Sampler object for external shear (component 2).

        Examples
        --------
        >>> from ler.lens_galaxy_population import OpticalDepth
        >>> od = OpticalDepth()
        >>> gamma2 = od.external_shear2(size=100)
        """

        return self._external_shear2

    @external_shear2.setter
    def external_shear2(self, prior):
        if prior in self.available_lens_priors["external_shear2"]:
            print(f"using ler available external_shear2 function : {prior}")
            args = self.lens_priors_params["external_shear2"]
            if args is None:
                self._external_shear2 = getattr(self, prior)(
                    size=None, get_attribute=True
                )
            else:
                self._external_shear2 = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user provided custom external_shear2 class/object of type ler.utils.FunctionConditioning"
            )
            self._external_shear2 = prior
        elif (
            hasattr(prior, "rvs")
            and callable(getattr(prior, "rvs"))
            and hasattr(prior, "pdf")
            and callable(getattr(prior, "pdf"))
        ):
            self._external_shear2 = prior
        else:
            raise ValueError(
                "external_shear2 must be either a string in available_lens_priors['external_shear2'] "
                "or an instance providing both `.rvs(...)` and `.pdf(...)` (e.g. `ler.utils.FunctionConditioning`)."
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
        if prior in self.available_lens_priors["density_profile_slope"]:
            print(f"using ler available density_profile_slope function : {prior}")
            args = self.lens_priors_params["density_profile_slope"]
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
        elif (
            hasattr(prior, "rvs")
            and callable(getattr(prior, "rvs"))
            and hasattr(prior, "pdf")
            and callable(getattr(prior, "pdf"))
        ):
            self._density_profile_slope = prior
        else:
            raise ValueError(
                "density_profile_slope must be either a string in available_lens_priors['density_profile_slope'] "
                "or an instance providing both `.rvs(...)` and `.pdf(...)` (e.g. `ler.utils.FunctionConditioning`)."
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
            print(f"using ler available cross_section function : {cross_section}")
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
        elif cross_section == "cross_section_epl_shear_njit":
            print(f"using ler available cross_section function : {cross_section}")
            args = self.lens_functions_params["cross_section"]
            if args is None:
                self._cross_section = self.cross_section_epl_shear_njit(
                    zs=None,
                    zl=None,
                    sigma=None,
                    q=None,
                    phi=None,
                    gamma=None,
                    gamma1=None,
                    gamma2=None,
                    get_attribute=True,
                )
            else:
                self._cross_section = self.cross_section_epl_shear_njit(
                    zs=None,
                    zl=None,
                    sigma=None,
                    q=None,
                    phi=None,
                    gamma=None,
                    gamma1=None,
                    gamma2=None,
                    get_attribute=True,
                    **args,
                )
        elif cross_section in ["cross_section_sie_feixu", "cross_section_sis"]:
            print(f"using ler available cross_section function : {cross_section}")
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
    def available_lens_priors(self):
        """
        Dictionary of available lens parameter samplers and their default parameters.

        Returns
        -------
        available_lens_priors : ``dict``
            Dictionary with sampler names and default parameters.
        """

        self._available_lens_priors = dict(
            zs_sl=dict(
                strongly_lensed_source_redshift=dict(
                    tau_approximation=True,
                ),
            ),
            lens_redshift_sl=dict(
                lens_redshift_strongly_lensed_sis_analytical=None,
                lens_redshift_strongly_lensed_numerical=dict(
                    param_name = "lens_redshift_sl",
                    sampler_type = "lens_redshift_strongly_lensed_numerical",
                    lens_type = "epl_shear_galaxy",
                    integration_size=50000,
                    use_multiprocessing=False,
                    cross_section_epl_shear_interpolation=False,
                ),
            ),
            lens_redshift=dict(
                lens_redshift_intrinsic_numerical=None,
            ),
            velocity_dispersion=dict(
                gengamma=dict(
                    param_name="velocity_dispersion",
                    sampler_type="gengamma",
                    sigma_min=100.0,
                    sigma_max=400.0,
                    alpha=0.94,
                    beta=1.85,
                    phistar=2.099e-2 * (self.cosmo.h/0.7)**3,
                    sigmastar=113.78,
                ),
                velocity_dispersion_choi=dict(
                    param_name="velocity_dispersion",
                    sampler_type="velocity_dispersion_choi",
                    sigma_min=100.0,
                    sigma_max=400.0,
                    alpha=2.32,
                    beta=2.67,
                    phistar=8.0e-3 * self.cosmo.h**3,
                    sigmastar=161.0,
                ),
                velocity_dispersion_bernardi=dict(
                    param_name="velocity_dispersion",
                    sampler_type="velocity_dispersion_bernardi",
                    sigma_min=100.0,
                    sigma_max=400.0,
                    alpha=0.94,
                    beta=1.85,
                    phistar=2.099e-2 * (self.cosmo.h/0.7)**3,
                    sigmastar=113.78,
                ),
                velocity_dispersion_ewoud=dict(
                    param_name="velocity_dispersion",
                    sampler_type="velocity_dispersion_ewoud",
                    sigma_min=100.0,
                    sigma_max=400.0,
                    alpha=0.94,
                    beta=1.85,
                    phistar=2.099e-2 * (self.cosmo.h/0.7)**3,
                    sigmastar=113.78,
                ),
            ),
            axis_ratio=dict(
                rayleigh=dict(param_name="axis_ratio", sampler_type="rayleigh", q_min=0.2, q_max=1.0),
                axis_ratio_padilla_strauss=dict(param_name="axis_ratio", sampler_type="axis_ratio_padilla_strauss", q_min=0.2, q_max=1.0),
                uniform=dict(param_name="axis_ratio", sampler_type="uniform", x_min=0.2, x_max=1.0),
                constant_values_n_size=dict(param_name="axis_ratio", sampler_type="constant_values_n_size", value=1.0),
            ),
            axis_rotation_angle=dict(
                uniform=dict(param_name="axis_rotation_angle", sampler_type="uniform", x_min=0.0, x_max=2 * np.pi),
                constant_values_n_size=dict(param_name="axis_rotation_angle", sampler_type="constant_values_n_size", value=0.0)
            ),
            density_profile_slope=dict[str, dict[str, str | float]](
                normal=dict(param_name="density_profile_slope", sampler_type="normal", mu=1.99, sigma=0.149),
                constant_values_n_size=dict(param_name="density_profile_slope", sampler_type="constant_values_n_size", value=2.0),
            ),
            external_shear1=dict(
                normal=dict(param_name="external_shear1", sampler_type="normal", mu=0.0, sigma=0.05),
                constant_values_n_size=dict(param_name="external_shear1", sampler_type="constant_values_n_size", value=0.0)
            ),
            external_shear2=dict(
                normal=dict(param_name="external_shear2", sampler_type="normal", mu=0.0, sigma=0.05),
                constant_values_n_size=dict(param_name="external_shear2", sampler_type="constant_values_n_size", value=0.0)
            ),
            # source_parameters=dict(gw_parameters_rvs=None),
        )

        return self._available_lens_priors

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
                rejection_sampler_full=dict(
                    n_prop=10000,
                    threshold_factor=1e-4,
                    zs_min=0.001,
                    zs_max=10.0,
                    zl_min=0.0001,
                    zl_max=None,
                    sigma_min=100.0,
                    sigma_max=400.0,
                    q_min=0.2,
                    q_max=1.0,
                    phi_min=0.0,
                    phi_max=2 * np.pi,
                    gamma_min=1.4,
                    gamma_max=2.7,
                    shear_min=-0.2,
                    shear_max=0.2,
                ),
                importance_sampler_full=dict(
                    n_prop=500,
                    threshold_factor=1e-4,
                    zs_min=0.001,
                    zs_max=10.0,
                    zl_min=0.0001,
                    zl_max=None,
                    sigma_min=100.0,
                    sigma_max=400.0,
                    q_min=0.2,
                    q_max=1.0,
                    phi_min=0.0,
                    phi_max=2 * np.pi,
                    gamma_min=1.4,
                    gamma_max=2.7,
                    shear_min=-0.2,
                    shear_max=0.2,
                ),
                rejection_sampler_partial=dict(
                    n_prop=10000,
                    threshold_factor=1e-4,
                    sigma_min=100.0,
                    sigma_max=400.0,
                    q_min=0.2,
                    q_max=1.0,
                    phi_min=0.0,
                    phi_max=2 * np.pi,
                    gamma_min=1.4,
                    gamma_max=2.7,
                    shear_min=-0.2,
                    shear_max=0.2,
                ),
                importance_sampler_partial=dict(
                    n_prop=400,
                    threshold_factor=1e-4,
                    sigma_min=100.0,
                    sigma_max=400.0,
                    q_min=0.2,
                    q_max=1.0,
                    phi_min=0.0,
                    phi_max=2 * np.pi,
                    gamma_min=1.4,
                    gamma_max=2.7,
                    shear_min=-0.2,
                    shear_max=0.2,
                ),
            ),
            optical_depth=dict(
                optical_depth_sis_analytic=dict(param_name="optical_depth", function_type="optical_depth_sis_analytic"),
                optical_depth_numerical=dict(param_name="optical_depth", function_type="optical_depth_numerical"),
            ),
            param_sampler_type=dict(
                epl_shear_sl_parameters_rvs=None,
            ),
            cross_section=dict(
                cross_section_sie_feixu=None,
                cross_section_sis=None,
                cross_section_epl_shear_numerical=None,
                cross_section_epl_shear_interpolation=dict(
                    num_th=500, maginf=-100.0
                ),
                cross_section_epl_shear_njit=dict(
                    num_th=500, maginf=-100.0
                ),
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
        from numba import set_num_threads
        set_num_threads(value)

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
