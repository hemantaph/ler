# -*- coding: utf-8 -*-
"""
Module for sampling compact binary coalescence (CBC) source parameters.

This module provides the ``CBCSourceParameterDistribution`` class for generating
complete sets of intrinsic and extrinsic gravitational wave parameters for
compact binary sources (BBH, BNS, NSBH). It includes mass distributions,
spin parameters, sky positions, and other parameters needed for GW analysis.

Inheritance hierarchy: \n
- :class:`~ler.gw_source_population.CBCSourceRedshiftDistribution` \n

Usage:
    Basic workflow example:

    >>> from ler.gw_source_population import CBCSourceParameterDistribution
    >>> cbc = CBCSourceParameterDistribution(event_type='BBH')
    >>> params = cbc.gw_parameters_rvs(size=1000)
    >>> print(list(params.keys()))
"""

import warnings

warnings.filterwarnings("ignore")
import numpy as np
from numba import njit

# for redshift to luminosity distance conversion
from astropy.cosmology import LambdaCDM

# from gwcosmo import priors as p
# scipy.integrate.quad removed - was unused

# for multiprocessing
# Import helper routines
from ..utils import FunctionConditioning

# import redshift distribution sampler
from .cbc_source_redshift_distribution import CBCSourceRedshiftDistribution

chunk_size = 10000  # chunk size for rejection sampling


class CBCSourceParameterDistribution(CBCSourceRedshiftDistribution):
    """
    Class for sampling compact binary coalescence source parameters.

    This class generates complete sets of intrinsic and extrinsic gravitational
    wave parameters for compact binary sources including masses, spins, sky
    positions, and orbital parameters. It supports BBH, BNS, NSBH, and primordial
    black hole populations with configurable prior distributions.

    Key Features: \n
    - Multiple mass distribution models (PowerLaw+Gaussian, lognormal, bimodal) \n
    - Configurable spin priors (zero, aligned, precessing) \n
    - Isotropic sky position and orientation sampling \n
    - Built-in support for population III and primordial black holes \n

    Parameters
    ----------
    npool : ``int``
        Number of processors to use for multiprocessing and numba threads. \n
        default: 4
    z_min : ``float``
        Minimum redshift of the source population. \n
        default: 0.0
    z_max : ``float``
        Maximum redshift of the source population. \n
        default: 10.0
    event_type : ``str``
        Type of compact binary event to generate. \n
        Options: \n
        - 'BBH': Binary black hole (Population I/II) \n
        - 'BNS': Binary neutron star \n
        - 'NSBH': Neutron star-black hole \n
        - 'BBH_popIII': Population III binary black hole \n
        - 'BBH_primordial': Primordial binary black hole \n
        default: 'BBH'
    gw_priors : ``dict`` or ``None``
        Dictionary of prior sampler functions for each parameter. \n
        If None, uses default priors based on event_type. \n
        default: None
    gw_priors_params : ``dict`` or ``None``
        Dictionary of parameters for each prior sampler function. \n
        If None, uses default parameters based on event_type. \n
        default: None
    cosmology : ``astropy.cosmology`` or ``None``
        Cosmology to use for distance calculations. \n
        default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)
    spin_zero : ``bool``
        If True, spin parameters are set to zero (no spin sampling). \n
        default: False
    spin_precession : ``bool``
        If True (and spin_zero=False), sample precessing spin parameters. \n
        If False (and spin_zero=False), sample aligned/anti-aligned spins. \n
        default: False
    directory : ``str``
        Directory to store interpolator JSON files. \n
        default: './interpolator_json'
    create_new_interpolator : ``dict`` or ``bool``
        Configuration for creating new interpolators. \n
        If bool, applies to all interpolators. \n
        default: False

    Examples
    --------
    >>> from ler.gw_source_population import CBCSourceParameterDistribution
    >>> cbc = CBCSourceParameterDistribution(event_type='BBH')
    >>> params = cbc.gw_parameters_rvs(size=1000)
    >>> print(list(params.keys()))

    Instance Methods
    ----------
    CBCSourceParameterDistribution has the following methods: \n
    +-----------------------------------------------------+------------------------------------------------+
    | Method                                              | Description                                    |
    +=====================================================+================================================+
    | :meth:`~gw_parameters_rvs`                       | Sample all GW parameters for compact binaries  |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~powerlaw_plus_peak`.       | Sample BBH masses with PowerLaw+PEAK model     |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~binary_masses_BBH_popIII_lognormal`         | Sample pop III BBH masses from lognormal       |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~binary_masses_BBH_primordial_lognormal`     | Sample primordial BBH masses from lognormal    |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~binary_masses_NSBH_broken_powerlaw`         | Sample NSBH masses from broken powerlaw        |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~binary_masses_uniform`                      | Sample masses from uniform distribution        |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~binary_masses_BNS_bimodal`                  | Sample BNS masses from bimodal distribution    |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~constant_values_n_size`                     | Return array of constant values                |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~uniform`                            | Sample from uniform distribution               |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~sampler_cosine`                             | Sample from cosine distribution                |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~sampler_sine`                               | Sample from sine distribution                  |
    +-----------------------------------------------------+------------------------------------------------+

    Instance Attributes
    ----------
    CBCSourceParameterDistribution has the following attributes: \n
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | Attribute                                      | Type                   | Unit  | Description                                    |
    +================================================+========================+=======+================================================+
    | :attr:`~z_min`                                 | ``float``              |       | Minimum redshift of source population          |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~z_max`                                 | ``float``              |       | Maximum redshift of source population          |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~cosmo`                                 | ``astropy.cosmology``  |       | Cosmology for distance calculations            |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~spin_zero`                             | ``bool``               |       | Whether to ignore spin parameters              |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~spin_precession`                       | ``bool``               |       | Whether to use precessing spins                |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~directory`                             | ``str``                |       | Directory for interpolator files               |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~gw_param_samplers`                     | ``dict``               |       | Dictionary of parameter sampler functions      |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~gw_param_samplers_params`              | ``dict``               |       | Dictionary of sampler function parameters      |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~available_gw_prior`                    | ``dict``               |       | Available prior distributions                  |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~mass_1_source`                   | ``callable``           |       | Sampler for source frame masses                |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~zs`                                    | ``callable``           |       | Sampler for source redshift                    |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~geocent_time`                          | ``callable``           |       | Sampler for geocentric time                    |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~ra`                                    | ``callable``           |       | Sampler for right ascension                    |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~dec`                                   | ``callable``           |       | Sampler for declination                        |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~phase`                                 | ``callable``           |       | Sampler for coalescence phase                  |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~psi`                                   | ``callable``           |       | Sampler for polarization angle                 |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~theta_jn`                              | ``callable``           |       | Sampler for inclination angle                  |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~a_1`                                   | ``callable``           |       | Sampler for spin1 magnitude                    |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~a_2`                                   | ``callable``           |       | Sampler for spin2 magnitude                    |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~tilt_1`                                | ``callable``           |       | Sampler for tilt1 angle                        |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~tilt_2`                                | ``callable``           |       | Sampler for tilt2 angle                        |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~phi_12`                                | ``callable``           |       | Sampler for phi_12 angle                       |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    | :attr:`~phi_jl`                                | ``callable``           |       | Sampler for phi_jl angle                       |
    +------------------------------------------------+------------------------+-------+------------------------------------------------+
    """

    # Attributes
    z_min = None
    """``float`` \n
    Minimum redshift of the source population
    """

    z_max = None
    """``float`` \n
    Maximum redshift of the source population
    """

    event_type = None
    """``str`` \n
    Type of event to generate. \n
    e.g. 'BBH', 'BNS', 'NSBH'
    """

    gw_priors = None
    """``dict`` \n
    Dictionary of prior sampler functions.
    """

    gw_priors_params = None
    """``dict`` \n
    Dictionary of prior sampler functions' input parameters.
    """

    cosmo = LambdaCDM(
        H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0
    )
    """``astropy.cosmology`` \n
    Cosmology to use.
    """

    spin_zero = None
    """``bool`` \n
    If True, spin prior is set to zero.
    """

    def __init__(
        self,
        npool=4,
        z_min=0.0,
        z_max=10.0,
        event_type="BBH",
        gw_functions=None,
        gw_functions_params=None,
        gw_priors=None,
        gw_priors_params=None,
        cosmology=None,
        spin_zero=False,
        spin_precession=False,
        directory="./interpolator_json",
        create_new_interpolator=False,
        verbose=False,
    ):
        # set attributes
        self.z_min = z_min
        self.z_max = z_max
        self.cosmo = (
            cosmology
            if cosmology
            else LambdaCDM(
                H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0
            )
        )
        # note that self.cosmo is initialized in the super class
        self.spin_zero = spin_zero
        self.spin_precession = spin_precession
        self.directory = directory

        # setting up the interpolator creation parameters
        create_new_interpolator = self._setup_decision_dictionary_gw_params(
            create_new_interpolator
        )

        # dealing with prior functions and categorization
        (
            self.gw_param_samplers,
            self.gw_param_samplers_params,
            self.gw_functions,
            self.gw_functions_params,
        ) = self._gw_priors_categorization(
            event_type, gw_priors, gw_priors_params, gw_functions, gw_functions_params, 
        )

        # initialize the CBCSourceRedshiftDistribution parent class
        # for redshift distribution
        # instance attribute source_redshift is initialized here
        super().__init__(
            npool=npool,
            z_min=z_min,
            z_max=z_max,
            event_type=event_type,
            merger_rate_density=self.gw_functions["merger_rate_density"],
            merger_rate_density_param=self.gw_functions_params[
                "merger_rate_density"
            ],
            cosmology=cosmology,
            directory=directory,
            create_new_interpolator=create_new_interpolator,
        )

        # Function to sample lens parameters
        self.gw_parameters_rvs_routine = getattr(
            self, self.gw_functions["param_sampler_type"]
        )

        print("\nInitializing CBCSourceParameterDistribution class...\n")
        # initializing samplers
        # it goes through the setter functions and assign the sampler functions
        # self.merger_rate_density_based_source_redshift is already initialized in the super class
        self.zs = self.gw_param_samplers["zs"]
        self.mass_1_source = self.gw_param_samplers["mass_1_source"]
        if self.gw_param_samplers["mass_2_source"] is not None:
            self.mass_2_source = self.gw_param_samplers["mass_2_source"]
            self.mass_ratio = None
        else:
            self.mass_ratio = self.gw_param_samplers["mass_ratio"]
            self.mass_2_source = None
            
        self.geocent_time = self.gw_param_samplers["geocent_time"]
        self.ra = self.gw_param_samplers["ra"]
        self.dec = self.gw_param_samplers["dec"]
        self.phase = self.gw_param_samplers["phase"]
        self.psi = self.gw_param_samplers["psi"]
        self.theta_jn = self.gw_param_samplers["theta_jn"]
        # initialize the spin prior attribute
        self.a_1 = self.gw_param_samplers["a_1"]
        self.a_2 = self.gw_param_samplers["a_2"]
        self.tilt_1 = self.gw_param_samplers["tilt_1"]
        self.tilt_2 = self.gw_param_samplers["tilt_2"]
        self.phi_12 = self.gw_param_samplers["phi_12"]
        self.phi_jl = self.gw_param_samplers["phi_jl"]


        # Function to sample gw parameters
        print("Initializing GW parameter samplers...\n")
        self.sample_gw_parameters_routine = self._initialization_gw_parameters_sampler()

    def _setup_decision_dictionary_gw_params(self, create_new_interpolator):
        """
        Helper function to set up decision dictionary for interpolator creation.

        Parameters
        ----------
        create_new_interpolator : ``dict`` or ``bool``
            If dict, contains boolean values and resolution for each interpolator. \n
            If bool, applies to all interpolators.

        Returns
        -------
        create_new_interpolator_ : ``dict``
            Dictionary with keys for each parameter and values containing \n
            create_new (bool) and resolution (int) settings.
        """
        create_new_interpolator_ = dict(
            mass_1_source=dict(create_new=False, resolution=200, cdf_size=500),
            mass_ratio=dict(create_new=False, resolution=100, m1_resolution=200, cdf_size=200),
            mass_2_source=dict(create_new=False, resolution=200, cdf_size=500),
            geocent_time=dict(create_new=False, resolution=200, cdf_size=500),
            ra=dict(create_new=False, resolution=200, cdf_size=500),
            dec=dict(create_new=False, resolution=200, cdf_size=500),
            phase=dict(create_new=False, resolution=200, cdf_size=500),
            psi=dict(create_new=False, resolution=200, cdf_size=500),
            theta_jn=dict(create_new=False, resolution=200, cdf_size=500),
            a_1=dict(create_new=False, resolution=200, cdf_size=500),
            a_2=dict(create_new=False, resolution=200, cdf_size=500),
            tilt_1=dict(create_new=False, resolution=200, cdf_size=500),
            tilt_2=dict(create_new=False, resolution=200, tilt_1_resolution=100, cdf_size=200),
            phi_12=dict(create_new=False, resolution=200, cdf_size=500),
            phi_jl=dict(create_new=False, resolution=200, cdf_size=500),
            merger_rate_density=dict(create_new=False, resolution=200, cdf_size=500),
            redshift_distribution=dict(create_new=False, resolution=200, cdf_size=500),
            luminosity_distance=dict(create_new=False, resolution=500),
            differential_comoving_volume=dict(create_new=False, resolution=500),
        )
        if isinstance(create_new_interpolator, dict):
            for key, value in create_new_interpolator.items():
                if (
                    key in create_new_interpolator_
                    and isinstance(create_new_interpolator_[key], dict)
                    and isinstance(value, dict)
                ):
                    create_new_interpolator_[key].update(value)
                else:
                    create_new_interpolator_[key] = value
        elif create_new_interpolator is True:
            for key in create_new_interpolator_:
                create_new_interpolator_[key]["create_new"] = True

        return create_new_interpolator_

    def _gw_priors_categorization(
        self, event_type, gw_priors, gw_priors_params, gw_functions, gw_functions_params
    ):
        """
        Helper function to categorize event priors/functions and parameters.

        There are three levels.
        1. Set default gw functions and samplers based on event type. \n
        2. If there is user input, update the gw functions and samplers with user input. \n
        3. Check the validity of user input and update the gw functions and samplers. User input can be different from the level (1). It can either be one of the internal functions or a user-defined function. If it is one of the internal functions, the parameters will be updated (missing parameters will be taken care of) with user input if given, otherwise it will use the default parameters. If it is a user-defined function, it will be directly used without any parameter update. \n

        Parameters
        ----------
        event_type : ``str``
            Type of compact binary event. \n
            Options: 'BBH', 'BNS', 'NSBH', 'BBH_popIII', 'BBH_primordial'
        gw_priors : ``dict`` or ``None``
            User-provided prior sampler functions for each parameter.
        gw_priors_params : ``dict`` or ``None``
            User-provided parameters for each prior sampler.

        Returns
        -------
        gw_priors_ : ``dict``
            Dictionary of prior sampler functions for each parameter.
        gw_priors_params_ : ``dict``
            Dictionary of sampler parameters for each parameter.
        """

        # 1. Set default gw functions and samplers based on event type.

        # for BBH
        if event_type == "BBH":
            merger_rate_density_function = (
                "merger_rate_density_madau_dickinson_belczynski_ng"
            )
            merger_rate_density_function_params = dict(
                param_name = "merger_rate_density",
                function_type = "merger_rate_density_madau_dickinson_belczynski_ng",
                R0=19 * 1e-9, alpha_F=2.57, beta_F=5.83, c_F=3.36
            )
            param_sampler_type = "gw_parameters_rvs"
            param_sampler_type_params = None
            mass_1_source_priors = "broken_powerlaw_plus_2peaks"
            mass_1_source_priors_params = dict(
                param_name = "mass_1_source",
                sampler_type = "broken_powerlaw_plus_2peaks",
                lam_0=0.361, 
                lam_1=0.586, 
                mpp_1=9.764, 
                sigpp_1=0.649, 
                mpp_2=32.763, 
                sigpp_2=3.918, 
                mlow_1=5.059, 
                delta_m_1=4.321, 
                break_mass=35.622, 
                alpha_1=1.728, 
                alpha_2=4.512, 
                mmax=300.0
            )
            mass_ratio_prior = "powerlaw_with_smoothing"
            mass_ratio_prior_params = dict(
                param_name = "mass_ratio",
                sampler_type = "powerlaw_with_smoothing",
                q_min=0.01, q_max=1.0,beta=1.171, mmin=3.551, delta_m=4.910
            )
            mass_2_source_priors = None
            mass_2_source_priors_params = None
            a_max = 1.0

        elif event_type == "BNS":
            merger_rate_density_function = "merger_rate_density_madau_dickinson2014"
            merger_rate_density_function_params = dict(
                param_name = "merger_rate_density",
                function_type = "merger_rate_density_madau_dickinson2014",
                R0=89 * 1e-9, a=0.015, b=2.7, c=2.9, d=5.6
            )
            param_sampler_type = "gw_parameters_rvs"
            param_sampler_type_params = None
            mass_1_source_priors = "truncated_normal"
            mass_1_source_priors_params = dict(
                param_name = "mass_1_source",
                sampler_type = "truncated_normal",
                mu=1.4, sigma=0.68, x_min=1.0, x_max=2.5
            )
            mass_ratio_prior = None
            mass_ratio_prior_params = None
            mass_2_source_priors = "truncated_normal"
            mass_2_source_priors_params = dict(
                param_name = "mass_2_source",
                sampler_type = "truncated_normal",
                mu=1.4, sigma=0.68, x_min=1.0, x_max=2.5
            )
            a_max = 0.4

        elif event_type == "NSBH":
            merger_rate_density_function = "merger_rate_density_madau_dickinson2014"
            merger_rate_density_function_params = dict(
                param_name = "merger_rate_density",
                function_type = "merger_rate_density_madau_dickinson2014",
                R0=23 * 1e-9, a=0.015, b=2.7, c=2.9, d=5.6
            )
            param_sampler_type = "gw_parameters_rvs"
            param_sampler_type_params = None
            mass_1_source_priors = "uniform"
            mass_1_source_priors_params = dict(
                param_name = "mass_1_source",
                sampler_type = "uniform",
                x_min=3.0, x_max=60.0
            )
            mass_ratio_prior = None
            mass_ratio_prior_params = None
            mass_2_source_priors = "truncated_normal"
            mass_2_source_priors_params = dict(
                param_name = "mass_2_source",
                sampler_type = "truncated_normal",
                mu=1.4, sigma=0.68, x_min=1.0, x_max=2.5
            )
            a_max = 1.0

        # elif event_type == "BBH_popIII":
        #     merger_rate_density_function = "merger_rate_density_bbh_popIII_ken2022"
        #     merger_rate_density_function_params = dict(
        #         n0=19.2 * 1e-9, aIII=0.66, bIII=0.3, zIII=11.6
        #     )
        #     param_sampler_type = "sample_gw_parameters"
        #     param_sampler_type_params = None
        #     mass_1_source_priors = "binary_masses_BBH_popIII_lognormal"
        #     mass_1_source_priors_params = dict(
        #         m_min=5.0, m_max=150.0, Mc=30.0, sigma=0.3, chunk_size=chunk_size
        #     )
        #     a_max = 1.0

        # elif event_type == "BBH_primordial":
        #     merger_rate_density_function = "merger_rate_density_bbh_primordial_ken2022"
        #     merger_rate_density_function_params = dict(
        #         n0=0.044 * 1e-9, t0=13.786885302009708
        #     )
        #     param_sampler_type = "sample_gw_parameters"
        #     param_sampler_type_params = None
        #     mass_1_source_priors = "binary_masses_BBH_primordial_lognormal"
        #     mass_1_source_priors_params = dict(
        #         m_min=1.0, m_max=100.0, Mc=20.0, sigma=0.3, chunk_size=chunk_size
        #     )
        #     a_max = 1.0

        else:
            raise ValueError(f"event_type {event_type} is not recognized")

        # setting the priors and its parameters
        gw_priors_ = dict(
            zs="merger_rate_density_based_source_redshift",
            mass_1_source=mass_1_source_priors,
            mass_ratio=mass_ratio_prior,
            mass_2_source=mass_2_source_priors,
            geocent_time="uniform",
            ra="uniform",
            dec="sampler_cosine",
            phase="uniform",
            psi="uniform",
            theta_jn="sampler_sine",
        )

        gw_priors_params_ = dict(
            zs=None,
            mass_1_source=mass_1_source_priors_params,
            mass_ratio=mass_ratio_prior_params,
            mass_2_source=mass_2_source_priors_params,
            geocent_time=dict(x_min=1238166018, x_max=1269702018),
            ra=dict(x_min=0.0, x_max=2.0 * np.pi),
            dec=None,  # dict(x_min=-np.pi/2, x_max=np.pi/2),
            phase=dict(x_min=0.0, x_max=2.0 * np.pi),
            psi=dict(x_min=0.0, x_max=np.pi),
            theta_jn=None,  # dict(x_min=0., x_max=np.pi),
        )

        gw_functions_ = dict(
            merger_rate_density=merger_rate_density_function,
            param_sampler_type=param_sampler_type,
        )

        gw_functions_params_ = dict(
            merger_rate_density=merger_rate_density_function_params,
            param_sampler_type=param_sampler_type_params,
        )

        # spin
        if not self.spin_zero:

            if self.spin_precession:
                gw_priors_["a_1"] = "truncated_normal"
                gw_priors_["a_2"] = "truncated_normal"
                gw_priors_["tilt_1"] = "gaussian_plus_isotropic"
                gw_priors_["tilt_2"] = "gaussian_plus_isotropic_joint"
                gw_priors_["phi_12"] = "uniform"
                gw_priors_["phi_jl"] = "uniform"

                gw_priors_params_["a_1"] = dict(x_min=0.0, x_max=a_max, mu=0.085, sigma=0.330)
                gw_priors_params_["a_2"] = dict(x_min=0.0, x_max=a_max, mu=0.085, sigma=0.330)
                gw_priors_params_["tilt_1"] = dict(tilt_1_min=0.0, tilt_1_max=np.pi, mu_t=0.426, sigma_t=1.222, zeta=0.652)
                gw_priors_params_["tilt_2"] = dict(tilt_2_min=0.0, tilt_2_max=np.pi, mu_t=0.426, sigma_t=1.222, zeta=0.652)
                gw_priors_params_["phi_12"] = dict(x_min=0, x_max=2 * np.pi)
                gw_priors_params_["phi_jl"] = dict(x_min=0, x_max=2 * np.pi)
            else:
                gw_priors_["a_1"] = "uniform"
                gw_priors_["a_2"] = "uniform"
                gw_priors_["tilt_1"] = None
                gw_priors_["tilt_2"] = None
                gw_priors_["phi_12"] = None
                gw_priors_["phi_jl"] = None

                gw_priors_params_["a_1"] = dict(x_min=-a_max, x_max=a_max)
                gw_priors_params_["a_2"] = dict(x_min=-a_max, x_max=a_max)
                gw_priors_params_["tilt_1"] = None
                gw_priors_params_["tilt_2"] = None
                gw_priors_params_["phi_12"] = None
                gw_priors_params_["phi_jl"] = None
        else:   
            gw_priors_["a_1"] = None
            gw_priors_["a_2"] = None
            gw_priors_["tilt_1"] = None
            gw_priors_["tilt_2"] = None
            gw_priors_["phi_12"] = None
            gw_priors_["phi_jl"] = None

            gw_priors_params_["a_1"] = None
            gw_priors_params_["a_2"] = None
            gw_priors_params_["tilt_1"] = None
            gw_priors_params_["tilt_2"] = None
            gw_priors_params_["phi_12"] = None
            gw_priors_params_["phi_jl"] = None


        # 2. If there is user input, update the gw functions and samplers with user input.

        # update the priors if input is given
        if gw_priors:
            # If a built-in prior name is overridden without explicit parameter
            # overrides, discard the old event-type defaults so the new prior can
            # start from its own default parameter set.
            for key, value in gw_priors.items():
                if (
                    key in gw_priors_params_
                    and isinstance(value, str)
                    and (not gw_priors_params or key not in gw_priors_params)
                    and value != gw_priors_.get(key)
                ):
                    gw_priors_params_[key] = None
            gw_priors_.update(gw_priors)
        if gw_priors_params:
            gw_priors_params_.update(gw_priors_params)
        if gw_functions:
            gw_functions_.update(gw_functions)
        if gw_functions_params:
            gw_functions_params_.update(gw_functions_params)

        # # taking care of gw_priors_params from the available_gw_prior
        # for key, value in gw_priors_.items():
        #     if isinstance(value, str):
        #         dict_ = self.available_gw_prior[
        #             key
        #         ]  # e.g. all mass_1_source_priors function names and its parameters
        #         if value in dict_:
        #             param_dict = dict_[value]
        #             if gw_priors_params_[key] is None:
        #                 gw_priors_params_[key] = param_dict
        #             else:
        #                 param_dict.update(gw_priors_params_[key])
        #                 gw_priors_params_[key] = param_dict
        #         else:
        #             raise ValueError(
        #                 f"gw_priors_params_['{key}'] is not in available_gw_prior"
        #             )

        # taking care of gw_priors_params from the available_gw_prior
        for key, value in gw_priors_.items():
            if isinstance(value, str):
                dict_ = self.available_gw_prior[
                    key
                ]  # e.g. all source_frame_masses_prior function names and its parameters
                if value in dict_:
                    param_dict = dict_[value]
                    if gw_priors_params_[key] is None:
                        gw_priors_params_[key] = param_dict
                    else:
                        param_dict.update(gw_priors_params_[key])
                        gw_priors_params_[key] = param_dict
                else:
                    raise ValueError(
                        f"gw_priors_params_['{key}'] is not in available_gw_prior"
                    )
            elif not callable(value):
                if key in ["mass_ratio", "mass_2_source", 'a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl']:
                    if value is not None:
                        raise ValueError(
                            f"gw_priors_params_['{key}'] should be either a string name of available sampler or a function. If you don't want to use this parameter, set it to None."
                        )

                else:
                    raise ValueError(
                        f"gw_priors_params_['{key}'] should be either a string name of available sampler or a function"
                    )

        # 3. Check the validity of user input and update the gw functions and samplers. User input can be different from the level (1).

        # taking care of gw_functions_params from the available_gw_functions
        for key, value in gw_functions_.items():
            if value not in ["sample_gw_parameters"]:
                # if string check if it's in the available_gw_functions and get the parameters
                if isinstance(value, str):
                    dict_ = self.available_gw_functions[
                        key
                    ]  # e.g. all merger_rate_density function names and its parameters
                    # if the function is in the available_gw_functions, get the parameters and update the gw_functions_params_ dictionary
                    if value in dict_:
                        # get the default parameters
                        param_dict = dict_[value]
                        if gw_functions_params_[key] is None:
                            gw_functions_params_[key] = param_dict
                        else:
                            param_dict.update(gw_functions_params_[key])
                            gw_functions_params_[key] = param_dict
                    # if the function is not in the available_gw_prior, raise an error
                    else:
                        raise ValueError(
                            f"gw_functions_params_['{key}'] of '{value}' is not in available_gw_functions"
                        )
                # if it's not a string, it should be a function, check if it's callable
                elif not callable(value):
                    raise ValueError(
                        f"gw_functions_params_['{key}'] should be either a string name of available function or a function"
                    )

        return (gw_priors_, gw_priors_params_, gw_functions_, gw_functions_params_)

    def sample_gw_parameters(self, size=1000, param=None):
        """
        Sample all gravitational wave parameters for compact binaries.

        This method samples GW parameters using the specified routine and
        returns a dictionary with redshift, masses, luminosity distance, and
        orientation parameters.

        Parameters
        ----------
        size : ``int``
            Number of samples to draw.
            default: 1000
        param : ``dict`` or ``None``
            Dictionary of fixed parameter values. \n
            Parameters in this dict will not be sampled. \n
            default: None

        Returns
        -------
        gw_parameters : ``dict``
            Dictionary of sampled GW parameters. The included parameters and their units are as follows (for default settings):\n
            +--------------------+--------------+--------------------------------------+
            | Parameter          | Units        | Description                          |
            +====================+==============+======================================+
            | zs                 |              | redshift of the source               |
            +--------------------+--------------+--------------------------------------+
            | geocent_time       | s            | GPS time of coalescence              |
            +--------------------+--------------+--------------------------------------+
            | ra                 | rad          | right ascension                      |
            +--------------------+--------------+--------------------------------------+
            | dec                | rad          | declination                          |
            +--------------------+--------------+--------------------------------------+
            | phase              | rad          | phase of GW at reference frequency   |
            +--------------------+--------------+--------------------------------------+
            | psi                | rad          | polarization angle                   |
            +--------------------+--------------+--------------------------------------+
            | theta_jn           | rad          | inclination angle                    |
            +--------------------+--------------+--------------------------------------+
            | a_1                |              | spin_1 of the compact binary         |
            +--------------------+--------------+--------------------------------------+
            | a_2                |              | spin of the secondary compact binary |
            +--------------------+--------------+--------------------------------------+
            | luminosity_distance| Mpc          | luminosity distance                  |
            +--------------------+--------------+--------------------------------------+
            | mass_1_source      | Msun         | mass_1 of the compact binary         |
            |                    |              | (source frame)                       |
            +--------------------+--------------+--------------------------------------+
            | mass_2_source      | Msun         | mass_2 of the compact binary         |
            |                    |              | (source frame)                       |
            +--------------------+--------------+--------------------------------------+
            | mass_1             | Msun         | mass_1 of the compact binary         |
            |                    |              | (detector frame)                     |
            +--------------------+--------------+--------------------------------------+
            | mass_2             | Msun         | mass_2 of the compact binary         |
            |                    |              | (detector frame)                     |
            +--------------------+--------------+--------------------------------------+
        """

        if param is None:
            param = {}

        # Sample gw parameters using the specified routine
        result_dict = self.gw_parameters_rvs_routine(size, **param)

        return result_dict

    def gw_parameters_rvs(self, size=1000, **kwargs):
        """
        Sample all gravitational wave parameters for compact binaries.

        This method samples GW parameters using the specified routine and returns a dictionary with redshift, masses, luminosity distance, and
        orientation parameters. 

        Parameters
        ----------
        size : ``int``
            Number of samples to draw.
            default: 1000
        **kwargs : ``dict``
            Dictionary of fixed parameter values. \n
            Parameters in this dict will not be sampled. \n
            e.g. mass_1_source, mass_2_source, mass_ratio, zs, geocent_time, ra, dec, phase, psi, theta_jn, a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl    

       Returns
        -------
        gw_parameters : ``dict``
            Dictionary of sampled GW parameters. The included parameters and their units are as follows (for default settings):\n
            +--------------------+--------------+--------------------------------------+
            | Parameter          | Units        | Description                          |
            +====================+==============+======================================+
            | zs                 |              | redshift of the source               |
            +--------------------+--------------+--------------------------------------+
            | geocent_time       | s            | GPS time of coalescence              |
            +--------------------+--------------+--------------------------------------+
            | ra                 | rad          | right ascension                      |
            +--------------------+--------------+--------------------------------------+
            | dec                | rad          | declination                          |
            +--------------------+--------------+--------------------------------------+
            | phase              | rad          | phase of GW at reference frequency   |
            +--------------------+--------------+--------------------------------------+
            | psi                | rad          | polarization angle                   |
            +--------------------+--------------+--------------------------------------+
            | theta_jn           | rad          | inclination angle                    |
            +--------------------+--------------+--------------------------------------+
            | a_1                |              | spin_1 of the compact binary         |
            +--------------------+--------------+--------------------------------------+
            | a_2                |              | spin of the secondary compact binary |
            +--------------------+--------------+--------------------------------------+
            | luminosity_distance| Mpc          | luminosity distance                  |
            +--------------------+--------------+--------------------------------------+
            | mass_1_source      | Msun         | mass_1 of the compact binary         |
            |                    |              | (source frame)                       |
            +--------------------+--------------+--------------------------------------+
            | mass_2_source      | Msun         | mass_2 of the compact binary         |
            |                    |              | (source frame)                       |
            +--------------------+--------------+--------------------------------------+
            | mass_1             | Msun         | mass_1 of the compact binary         |
            |                    |              | (detector frame)                     |
            +--------------------+--------------+--------------------------------------+
            | mass_2             | Msun         | mass_2 of the compact binary         |
            |                    |              | (detector frame)                     |
            +--------------------+--------------+--------------------------------------+
        """

        result_dict = dict()

        # ----------------------------------------
        if 'zs' not in kwargs:
            result_dict['zs'] = self.zs.rvs(size=size)
        else:
            result_dict['zs'] = kwargs['zs']

        if 'mass_1_source' not in kwargs:
            if self.mass_1_source.rvs.__code__.co_argcount == 1:
                result_dict['mass_1_source'] = self.mass_1_source.rvs(size)  
            else:
                result_dict['mass_1_source'] = self.mass_1_source.rvs(size, result_dict['zs'])
        else:
            result_dict['mass_1_source'] = kwargs['mass_1_source']
        # ----------------------------------------
        if 'mass_2_source' not in kwargs:
            if self.mass_2_source is not None:
                if self.mass_2_source.rvs.__code__.co_argcount == 1:
                    result_dict['mass_2_source'] = self.mass_2_source.rvs(size)  
                else:
                    result_dict['mass_2_source'] = self.mass_2_source.rvs(size, result_dict['mass_1_source'])
            elif self.mass_ratio is not None:
                if self.mass_ratio.rvs.__code__.co_argcount == 1:
                    result_dict['mass_ratio'] = self.mass_ratio.rvs(size)  
                else:
                    result_dict['mass_ratio'] = self.mass_ratio.rvs(size, result_dict['mass_1_source'])
                result_dict['mass_2_source'] = result_dict['mass_1_source'] * result_dict['mass_ratio']
            else:
                raise ValueError("Either mass_2_source or mass_ratio should be defined.")
        else:
            result_dict['mass_2_source'] = kwargs['mass_2_source']
        # ----------------------------------------
        if 'geocent_time' not in kwargs:
            result_dict['geocent_time'] = self.geocent_time.rvs(size=size)
        else:
            result_dict['geocent_time'] = kwargs['geocent_time']
        # ----------------------------------------
        if 'ra' not in kwargs:
            result_dict['ra'] = self.ra.rvs(size=size)
        else:            
            result_dict['ra'] = kwargs['ra']
        # ----------------------------------------
        if 'dec' not in kwargs:
            result_dict['dec'] = self.dec.rvs(size=size)
        else:
            result_dict['dec'] = kwargs['dec']
        # ----------------------------------------
        if 'phase' not in kwargs:
            result_dict['phase'] = self.phase.rvs(size=size)
        else:
            result_dict['phase'] = kwargs['phase']
        # ----------------------------------------
        if 'psi' not in kwargs:
            result_dict['psi'] = self.psi.rvs(size=size)
        else:
            result_dict['psi'] = kwargs['psi']
        # ----------------------------------------
        if 'theta_jn' not in kwargs:
            result_dict['theta_jn'] = self.theta_jn.rvs(size=size)
        else:
            result_dict['theta_jn'] = kwargs['theta_jn']
        # ----------------------------------------
        if not self.spin_zero:
            # ----------------------------------------
            if 'a_1' not in kwargs:
                result_dict['a_1'] = self.a_1.rvs(size=size)
            else:
                result_dict['a_1'] = kwargs['a_1']
            # ----------------------------------------
            if 'a_2' not in kwargs:
                result_dict['a_2'] = self.a_2.rvs(size=size)
            else:
                result_dict['a_2'] = kwargs['a_2']
            # ----------------------------------------
            if self.spin_precession:
                # ----------------------------------------
                if 'tilt_1' not in kwargs:
                    result_dict['tilt_1'] = self.tilt_1.rvs(size=size)
                else:
                    result_dict['tilt_1'] = kwargs['tilt_1']
                # ----------------------------------------
                if 'tilt_2' not in kwargs:
                    if self.tilt_2.rvs.__code__.co_argcount == 1:
                        result_dict['tilt_2'] = self.tilt_2.rvs(size=size)
                    else:
                        result_dict['tilt_2'] = self.tilt_2.rvs(size, result_dict['tilt_1'])
                else:
                    result_dict['tilt_2'] = kwargs['tilt_2']
                # ----------------------------------------
                if 'phi_12' not in kwargs:
                    result_dict['phi_12'] = self.phi_12.rvs(size=size)
                else:
                    result_dict['phi_12'] = kwargs['phi_12']
                # ----------------------------------------
                if 'phi_jl' not in kwargs:
                    result_dict['phi_jl'] = self.phi_jl.rvs(size=size)
                else:
                    result_dict['phi_jl'] = kwargs['phi_jl']
                # ----------------------------------------

        # Compute derived quantities needed by pdet_finder
        zs = result_dict['zs']
        result_dict['luminosity_distance'] = self.luminosity_distance(zs)
        result_dict['mass_1'] = result_dict['mass_1_source'] * (1 + zs)
        result_dict['mass_2'] = result_dict['mass_2_source'] * (1 + zs)

        return result_dict

    def gw_parameters_rvs_njit(self, size=1000, **kwargs):
        """
        Sample all gravitational wave parameters for compact binaries.

        This method samples GW parameters using the specified routine and
        returns a dictionary with redshift, masses, luminosity distance, and
        orientation parameters.

        Parameters
        ----------
        size : ``int``
            Number of samples to draw.
            default: 1000

        Returns
        -------
        gw_parameters : ``dict``
            Dictionary of sampled GW parameters. The included parameters and their units are as follows (for default settings):\n
            +--------------------+--------------+--------------------------------------+
            | Parameter          | Units        | Description                          |
            +====================+==============+======================================+
            | zs                 |              | redshift of the source               |
            +--------------------+--------------+--------------------------------------+
            | geocent_time       | s            | GPS time of coalescence              |
            +--------------------+--------------+--------------------------------------+
            | ra                 | rad          | right ascension                      |
            +--------------------+--------------+--------------------------------------+
            | dec                | rad          | declination                          |
            +--------------------+--------------+--------------------------------------+
            | phase              | rad          | phase of GW at reference frequency   |
            +--------------------+--------------+--------------------------------------+
            | psi                | rad          | polarization angle                   |
            +--------------------+--------------+--------------------------------------+
            | theta_jn           | rad          | inclination angle                    |
            +--------------------+--------------+--------------------------------------+
            | a_1                |              | spin_1 of the compact binary         |
            +--------------------+--------------+--------------------------------------+
            | a_2                |              | spin of the secondary compact binary |
            +--------------------+--------------+--------------------------------------+
            | luminosity_distance| Mpc          | luminosity distance                  |
            +--------------------+--------------+--------------------------------------+
            | mass_1_source      | Msun         | mass_1 of the compact binary         |
            |                    |              | (source frame)                       |
            +--------------------+--------------+--------------------------------------+
            | mass_2_source      | Msun         | mass_2 of the compact binary         |
            |                    |              | (source frame)                       |
            +--------------------+--------------+--------------------------------------+
            | mass_1             | Msun         | mass_1 of the compact binary         |
            |                    |              | (detector frame)                     |
            +--------------------+--------------+--------------------------------------+
            | mass_2             | Msun         | mass_2 of the compact binary         |
            |                    |              | (detector frame)                     |
            +--------------------+--------------+--------------------------------------+
        """

        # Sample gw parameters using the specified routine
        gw_parameters = self.sample_gw_parameters_routine(size=size)

        zs = gw_parameters[0]
        m1 = gw_parameters[1]
        m2 = gw_parameters[2]

        result_dict = dict(
            zs=zs,
            mass_1_source=m1,
            mass_2_source=m2,
            mass_1=m1 * (1 + zs),
            mass_2=m2 * (1 + zs),
            luminosity_distance=self.luminosity_distance(zs),
            geocent_time=gw_parameters[3],
            ra=gw_parameters[4],
            dec=gw_parameters[5],
            phase=gw_parameters[6],
            psi=gw_parameters[7],
            theta_jn=gw_parameters[8],
        )

        if not self.spin_zero:
            result_dict["a_1"] = gw_parameters[9]
            result_dict["a_2"] = gw_parameters[10]
            if self.spin_precession:
                result_dict["tilt_1"] = gw_parameters[11]
                result_dict["tilt_2"] = gw_parameters[12]
                result_dict["phi_12"] = gw_parameters[13]
                result_dict["phi_jl"] = gw_parameters[14]

        return result_dict

    def _initialization_gw_parameters_sampler(self):
        """
        Helper function to initialize the gw parameter samplers.
        """

        from .prior_functions import _njit_checks, create_gw_parameters_sampler

        zs_rvs = self.zs.rvs
        m1_rvs = self.mass_1_source.rvs
        if self.mass_2_source is not None:
            q_rvs = None
            m2_rvs = self.mass_2_source.rvs
        else:
            q_rvs = self.mass_ratio.rvs
            m2_rvs = None
        tc_rvs = self.geocent_time.rvs
        ra_rvs = self.ra.rvs
        dec_rvs = self.dec.rvs
        phase_rvs = self.phase.rvs
        psi_rvs = self.psi.rvs
        theta_jn_rvs = self.theta_jn.rvs
        if not self.spin_zero:
            a_1_rvs = self.a_1.rvs
            a_2_rvs = self.a_2.rvs
            if self.spin_precession:
                tilt_1_rvs = self.tilt_1.rvs
                tilt_2_rvs = self.tilt_2.rvs
                phi_12_rvs = self.phi_12.rvs
                phi_jl_rvs = self.phi_jl.rvs
            else:
                tilt_1_rvs = None
                tilt_2_rvs = None
                phi_12_rvs = None
                phi_jl_rvs = None
        else:
                a_1_rvs = None
                a_2_rvs = None
                tilt_1_rvs = None
                tilt_2_rvs = None
                phi_12_rvs = None
                phi_jl_rvs = None

        use_njit_sampler, dict_ = _njit_checks(
            zs_rvs,
            m1_rvs,
            q_rvs,
            m2_rvs,
            tc_rvs,
            ra_rvs,
            dec_rvs,
            phase_rvs,
            psi_rvs,
            theta_jn_rvs,
            a_1_rvs,
            a_2_rvs,
            tilt_1_rvs,
            tilt_2_rvs,
            phi_12_rvs,
            phi_jl_rvs,
            spin_zero=self.spin_zero,
            spin_precession=self.spin_precession,
        )

        gw_sampling_routine = create_gw_parameters_sampler(
            use_njit_sampler=use_njit_sampler,
            zs_rvs=dict_["zs_rvs"],
            m1_rvs=dict_["m1_rvs"],
            q_rvs=dict_["q_rvs"],
            m2_rvs=dict_["m2_rvs"],
            tc_rvs=dict_["tc_rvs"],
            ra_rvs=dict_["ra_rvs"],
            dec_rvs=dict_["dec_rvs"],
            phase_rvs=dict_["phase_rvs"],
            psi_rvs=dict_["psi_rvs"],
            theta_jn_rvs=dict_["theta_jn_rvs"],
            a_1_rvs=dict_["a_1_rvs"],
            a_2_rvs=dict_["a_2_rvs"],
            tilt_1_rvs=dict_["tilt_1_rvs"],
            tilt_2_rvs=dict_["tilt_2_rvs"],
            phi_12_rvs=dict_["phi_12_rvs"],
            phi_jl_rvs=dict_["phi_jl_rvs"],
            spin_zero=self.spin_zero,
            spin_precession=self.spin_precession,
        )

        return gw_sampling_routine

    # -------------------------------------
    # BBH Primary mass distribution models
    # -------------------------------------

    def broken_powerlaw_plus_2peaks(
        self,
        size,
        get_attribute=False,
        **kwargs,
    ):
        """
        Sample primary source mass with broken power-law + two peaks model for BBH.

        Implements the mass distribution model combining a broken power-law
        with two Gaussian peaks for modeling specific mass populations.

        Parameters
        ----------
        size : ``int``
            Number of samples to draw.
        get_attribute : ``bool``
            If True, return the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Model parameters: \n
            - lam_0: Fraction of broken power-law component, default: 0.361 \n
            - lam_1: Fraction of first Gaussian peak, default: 0.586 \n
            - mpp_1: Mean of first Gaussian peak (Msun), default: 9.764 \n
            - sigpp_1: Std dev of first Gaussian peak (Msun), default: 0.649 \n
            - mpp_2: Mean of second Gaussian peak (Msun), default: 32.763 \n
            - sigpp_2: Std dev of second Gaussian peak (Msun), default: 3.918 \n
            - mlow_1: Minimum primary mass (Msun), default: 5.059 \n
            - delta_m_1: Low-mass tapering range (Msun), default: 4.321 \n
            - break_mass: Power-law break mass (Msun), default: 35.622 \n
            - alpha_1: Power-law index below break, default: 1.728 \n
            - alpha_2: Power-law index above break, default: 4.512 \n
            - mmax: Maximum primary mass (Msun), default: 300.0 \n
            - beta: Mass ratio power-law index, default: 1.171 \n
            - mlow_2: Minimum secondary mass (Msun), default: 3.551 \n
            - delta_m_2: Secondary mass tapering range (Msun), default: 4.910

        Returns
        -------
        mass_1_source : ``numpy.ndarray``
            Array of primary masses in source frame (Msun).
        mass_2_source : ``numpy.ndarray``
            Array of secondary masses in source frame (Msun).

        Examples
        --------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution()
        >>> m1_src, m2_src = cbc.broken_powerlaw_plus_2peaks(size=1000)
        """

        from .broken_powerlaw_plus_2peaks import broken_powerlaw_plus_2peaks_function, broken_powerlaw_plus_2peaks_pdf

        # Get parameters
        identifier_dict = dict(
            param_name = "mass_1_source",
            sampler_type = "broken_powerlaw_plus_2peaks",
        )
        identifier_dict["resolution"] = self.create_new_interpolator[
            identifier_dict["param_name"]
        ]["resolution"]
        identifier_dict["cdf_size"] = self.create_new_interpolator[
            identifier_dict["param_name"]
        ]["cdf_size"]
        param_dict = self.available_gw_prior[identifier_dict["param_name"]][
            identifier_dict["sampler_type"]
        ].copy()
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        def m1_fn(m1):
            return broken_powerlaw_plus_2peaks_function(
                m1,
                lam_0=identifier_dict["lam_0"],
                lam_1=identifier_dict["lam_1"],
                mpp_1=identifier_dict["mpp_1"],
                sigpp_1=identifier_dict["sigpp_1"],
                mpp_2=identifier_dict["mpp_2"],
                sigpp_2=identifier_dict["sigpp_2"],
                mlow_1=identifier_dict["mlow_1"],
                delta_m_1=identifier_dict["delta_m_1"],
                break_mass=identifier_dict["break_mass"],
                alpha_1=identifier_dict["alpha_1"],
                alpha_2=identifier_dict["alpha_2"],
                mmax=identifier_dict["mmax"],
                normalization_size=identifier_dict["normalization_size"],
            )
        
        m1_arr = np.geomspace(
            identifier_dict["mlow_1"], identifier_dict["mmax"], identifier_dict["resolution"]
        )

        m1_object = FunctionConditioning(
            function=m1_fn,
            x_array=m1_arr,
            non_negative_function=True,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["sampler_type"],
            create_new=self.create_new_interpolator[identifier_dict["param_name"]][
                "create_new"
            ],
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
            cdf_size=identifier_dict["cdf_size"],
        )

        # return m1_object if get_attribute else m1_object.rvs(size)
        if get_attribute:
            m1_object.original_pdf_function = broken_powerlaw_plus_2peaks_pdf
            return m1_object
        else:
            return m1_object.rvs(size)

    def powerlaw_plus_peak(
        self,
        size,
        get_attribute=False,
        **kwargs,
    ):
        """
        Sample source masses with PowerLaw+PEAK model for Population I/II BBH.

        GWTC-3 parameters: mminbh=4.98, mmaxbh=112.5, alpha=3.78, mu_g=32.27, sigma_g=3.88, lambda_peak=0.03, delta_m=4.8, beta=0.81

        Implements the mass distribution model from LIGO-Virgo population analyses
        combining a power-law with a Gaussian peak component. Uses separate
        m1 and q interpolators with joint composite sampling.

        Parameters
        ----------
        size : ``int``
            Number of samples to draw.
        get_attribute : ``bool``
            If True, return the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Model parameters: \n
            - mminbh: Minimum BH mass (Msun), default: 4.98 \n
            - mmaxbh: Maximum BH mass (Msun), default: 112.5 \n
            - alpha: Power-law spectral index, default: 3.78 \n
            - mu_g: Gaussian peak mean (Msun), default: 32.27 \n
            - sigma_g: Gaussian peak width (Msun), default: 3.88 \n
            - lambda_peak: Fraction in Gaussian component, default: 0.03 \n
            - delta_m: Low-mass tapering range (Msun), default: 4.8 \n
            - beta: Mass ratio power-law index, default: 0.81

        Returns
        -------
        mass_1_source : ``numpy.ndarray``
            Array of primary masses in source frame (Msun).
        mass_2_source : ``numpy.ndarray``
            Array of secondary masses in source frame (Msun).

        Examples
        --------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution()
        >>> m1_src, m2_src = cbc.powerlaw_plus_peak(size=1000)
        """

        from .powerlaw_plus_peak import powerlaw_plus_peak_function

        # Get parameters
        identifier_dict = dict(
            param_name="mass_1_source",
            sampler_type="powerlaw_plus_peak",
        )
        identifier_dict["resolution"] = self.create_new_interpolator[
            identifier_dict["param_name"]
        ]["resolution"]
        identifier_dict["cdf_size"] = self.create_new_interpolator[
            identifier_dict["param_name"]
        ]["cdf_size"]
        param_dict = self.available_gw_prior[identifier_dict["param_name"]][
            identifier_dict["sampler_type"]
        ].copy()
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        def m1_fn(m1):
            return powerlaw_plus_peak_function(
                m1,
                mminbh=identifier_dict["mminbh"],
                mmaxbh=identifier_dict["mmaxbh"],
                alpha=identifier_dict["alpha"],
                mu_g=identifier_dict["mu_g"],
                sigma_g=identifier_dict["sigma_g"],
                lambda_peak=identifier_dict["lambda_peak"],
                delta_m=identifier_dict["delta_m"],
                normalization_size=identifier_dict["normalization_size"],
            )
        
        m1_arr = np.geomspace(
            identifier_dict["mminbh"], identifier_dict["mmaxbh"], identifier_dict["resolution"]
        )

        m1_object = FunctionConditioning(
            function=m1_fn,
            x_array=m1_arr,
            non_negative_function=True,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["sampler_type"],
            create_new=self.create_new_interpolator[identifier_dict["param_name"]][
                "create_new"
            ],
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
            cdf_size=identifier_dict["cdf_size"],
        )

        return m1_object if get_attribute else m1_object.rvs(size)


    # ---------------------------------
    # Mass ratio distribution models
    # ---------------------------------

    def powerlaw_with_smoothing(
        self,
        size,
        get_attribute=False,
        **kwargs,
    ):
        """
        Sample mass ratio with power-law distribution with smoothing near q=1.

        Implements a power-law distribution for the mass ratio with an additional
        smoothing function to account for the observed excess of equal-mass binaries.

        Parameters
        ----------
        size : ``int``
            Number of samples to draw.
        get_attribute : ``bool``
            If True, return the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Model parameters: \n
            - q_min: Minimum mass ratio, default: 0.01 \n
            - q_max: Maximum mass ratio, default: 1.0 \n
            - beta: Power-law index, default: 1.171 \n
            - mlow_2: Minimum secondary mass (Msun), default: 3.551 \n
            - delta_m_2: Tapering range for secondary mass (Msun), default: 4.910

        Returns
        -------
        q : ``numpy.ndarray``
            Array of mass ratios (m2/m1) for the binary systems.

        Examples
        --------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution()
        >>> m1 = cbc.broken_powerlaw_plus_2peaks(1000)
        >>> q_samples = cbc.powerlaw_with_smoothing(1000, m1)
        """

        from .prior_functions import powerlaw_with_smoothing
        # powerlaw_with_smoothing(q, m2, mmin, beta, delta_m)

        # Get parameters
        identifier_dict = dict(
            param_name = "mass_ratio",
            sampler_type = "powerlaw_with_smoothing",
        )
        identifier_dict["resolution"] = self.create_new_interpolator[identifier_dict["param_name"]][
            "resolution"
        ]
        identifier_dict["cdf_size"] = self.create_new_interpolator[
            identifier_dict["param_name"]
        ]["cdf_size"]
        identifier_dict["m1_resolution"] = self.create_new_interpolator[identifier_dict["param_name"]][
            "m1_resolution"
        ]
        param_dict = self.available_gw_prior[identifier_dict["param_name"]][
            identifier_dict["sampler_type"]
        ].copy()

        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        if 'mlow_2' in identifier_dict:
            mmin = identifier_dict["mlow_2"]
        elif 'mminbh' in identifier_dict:
            mmin = identifier_dict["mminbh"]
        elif 'mmin' in identifier_dict:
            mmin = identifier_dict["mmin"]
        elif 'mlow_1' in identifier_dict:
            mmin = identifier_dict["mlow_1"]
        else:
            mmin = 3.551

        def powerlaw_with_smoothing_(q, m):

            return powerlaw_with_smoothing(
                q,
                q*m,
                mmin=mmin,
                beta=identifier_dict["beta"],
                delta_m=identifier_dict["delta_m"],
            )

        # create pdf 2D array
        size1 = identifier_dict["m1_resolution"]
        size2 = identifier_dict["resolution"]

        # m1 1D array
        identifier_dict_m1 = self.gw_param_samplers_params["mass_1_source"]
        if 'mlow_1' in identifier_dict_m1:
            mlow_1 = identifier_dict_m1["mlow_1"]
        elif 'mminbh' in identifier_dict_m1:
            mlow_1 = identifier_dict_m1["mminbh"]
        elif 'mmin' in identifier_dict_m1:
            mlow_1 = identifier_dict_m1["mmin"]
        else:
            mlow_1 = 5.059

        if 'mmax' in identifier_dict_m1:
            mmax = identifier_dict_m1["mmax"]
        elif 'mmaxbh' in identifier_dict_m1:
            mmax = identifier_dict_m1["mmaxbh"]
        else:
            mmax = 300.0
        m1_arr = np.geomspace(mlow_1, mmax, size1)

        if 'mlow_2' in identifier_dict:
            mlow_2 = identifier_dict["mlow_2"]
        else:            
            mlow_2 = 3.551 

        # q 2D array
        q_arr = np.linspace(
            mlow_2/m1_arr, 1.0, size2).T

        q_object = FunctionConditioning(
            function=powerlaw_with_smoothing_,  # can also be an array of function values
            x_array=q_arr,
            conditioned_y_array=m1_arr,  # if this is not none, 2D interpolation will be used
            non_negative_function=True,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["sampler_type"],
            create_new=self.create_new_interpolator[identifier_dict["param_name"]][
                "create_new"
            ],
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='rvs',
            cdf_size=identifier_dict["cdf_size"],
        )

        return q_object if get_attribute else q_object.rvs(size, m1_arr)
    
    def gaussian_plus_isotropic(self, size, get_attribute=False, **kwargs):
        """
        Sample from tilt_1 distribution for spin tilt angle of primary mass.

        Samples values in range [-1, 1] and converts to tilt angles in radians. Combines a truncated Gaussian component with an isotropic distribution.  

        Parameters
        ----------
        size : ``int``
            Number of samples to draw.
        get_attribute : ``bool``
            If True, return the sampler object instead of samples. \n
            default: False

        Returns
        -------
        values : ``numpy.ndarray``
            Array of tilt_1 angles.
        """

        from .prior_functions import gaussian_plus_isotropic_pdf

        # Get parameters
        identifier_dict = dict(
            param_name = "tilt_1",
            sampler_type = "gaussian_plus_isotropic",
        )
        identifier_dict["resolution"] = self.create_new_interpolator[identifier_dict["param_name"]][
            "resolution"
        ]
        identifier_dict["cdf_size"] = self.create_new_interpolator[
            identifier_dict["param_name"]
        ]["cdf_size"]
        param_dict = self.available_gw_prior[identifier_dict["param_name"]][
            identifier_dict["sampler_type"]
        ].copy()

        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        mu_t = identifier_dict["mu_t"]
        sigma_t = identifier_dict["sigma_t"]
        zeta = identifier_dict["zeta"]
        tilt_1_min = identifier_dict["tilt_1_min"]
        tilt_1_max = identifier_dict["tilt_1_max"]

        def tilt_1_pdf_(x):
            return gaussian_plus_isotropic_pdf(
                np.cos(x), mu_t=mu_t, sigma_t=sigma_t, zeta=zeta
            ) * np.sin(x) # np.sin(x) is the Jacobian for converting from cos(tilt) to tilt

        tilt_1_arr = np.linspace(tilt_1_min, tilt_1_max, identifier_dict["resolution"])

        tilt_1_object = FunctionConditioning(
            function=tilt_1_pdf_,  # can also be an array of function values
            x_array=tilt_1_arr,
            non_negative_function=True,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["sampler_type"],
            create_new=self.create_new_interpolator[identifier_dict["param_name"]][
                "create_new"
            ],
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='rvs',
            cdf_size=identifier_dict["cdf_size"],
        )

        return tilt_1_object if get_attribute else tilt_1_object.rvs(size)

    def gaussian_plus_isotropic_joint(self, size, get_attribute=False, **kwargs):
        """
        Sample from tilt_2 distribution for spin tilt angle of secondary mass, given the tilt_1 distribution for primary mass.

        Samples values in range [-1, 1] and converts to tilt angles in radians. Combines a truncated Gaussian component with an isotropic distribution.  

        Parameters
        ----------
        size : ``int``
            Number of samples to draw.
        get_attribute : ``bool``
            If True, return the sampler object instead of samples. \n
            default: False

        Returns
        -------
        values : ``numpy.ndarray``
            Array of tilt_2 angles.
        """

        from .prior_functions import gaussian_plus_isotropic_joint_pdf

        # Get parameters
        identifier_dict = dict(
            param_name = "tilt_2",
            sampler_type = "gaussian_plus_isotropic_joint",
        )
        identifier_dict["resolution"] = self.create_new_interpolator[identifier_dict["param_name"]][
            "resolution"
        ]
        identifier_dict["cdf_size"] = self.create_new_interpolator[
            identifier_dict["param_name"]
        ]["cdf_size"]
        param_dict = self.available_gw_prior[identifier_dict["param_name"]][
            identifier_dict["sampler_type"]
        ].copy()

        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        mu_t = identifier_dict["mu_t"]
        sigma_t = identifier_dict["sigma_t"]
        zeta = identifier_dict["zeta"]
        tilt_2_min = identifier_dict["tilt_2_min"]
        tilt_2_max = identifier_dict["tilt_2_max"]

        def tilt_2_pdf(tilt_2, tilt_1):
            return gaussian_plus_isotropic_joint_pdf(
                np.cos(tilt_2), np.cos(tilt_1), mu_t=mu_t, sigma_t=sigma_t, zeta=zeta
            ) * np.sin(tilt_2) # np.sin(tilt_2) is the Jacobian for converting from cos(tilt) to tilt

        tilt_1_params = self.gw_param_samplers_params.get("tilt_1")
        if tilt_1_params is None:
            tilt_1_params = self.available_gw_prior["tilt_1"][
                "gaussian_plus_isotropic"
            ]
        tilt_1_min = tilt_1_params.get("tilt_1_min", 0.0)
        tilt_1_max = tilt_1_params.get("tilt_1_max", np.pi)

        tilt_2_arr = np.linspace(tilt_2_min, tilt_2_max, identifier_dict["resolution"])
        tilt_1_arr = np.linspace(
            tilt_1_min,
            tilt_1_max,
            self.create_new_interpolator["tilt_2"]["tilt_1_resolution"],
        )

        tilt_2_object = FunctionConditioning(
            function=tilt_2_pdf,  # can also be an array of function values
            x_array=tilt_2_arr,
            conditioned_y_array=tilt_1_arr,
            non_negative_function=True,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["sampler_type"],
            create_new=self.create_new_interpolator[identifier_dict["param_name"]][
                "create_new"
            ],
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='rvs',
            cdf_size=identifier_dict["cdf_size"],
        )

        return tilt_2_object if get_attribute else tilt_2_object.rvs(size)

    # -------------------------------------
    # BNS Primary mass distribution models
    # -------------------------------------

    def bimodal(
        self,
        size,
        get_attribute=False,
        **kwargs,
    ):
        """
        Sample BNS primary mass from bimodal Gaussian distribution.

        Based on Will M. Farr et al. 2020 Eq. 6 for neutron star mass
        distribution combining two Gaussian peaks.

        .. math::

            p(m) = w\\,\\mathcal{N}(m \\mid \\mu_L, \\sigma_L)
            + (1-w)\\,\\mathcal{N}(m \\mid \\mu_R, \\sigma_R),

        truncated to ``[mmin, mmax]``.

        Parameters
        ----------
        size : ``int``
            Number of samples to draw.
        get_attribute : ``bool``
            If True, return the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Model parameters: \n
            - w: Weight of left peak, default: 0.643 \n
            - muL: Mean of left peak (Msun), default: 1.352 \n
            - sigmaL: Width of left peak (Msun), default: 0.08 \n
            - muR: Mean of right peak (Msun), default: 1.88 \n
            - sigmaR: Width of right peak (Msun), default: 0.3 \n
            - mmin: Minimum mass (Msun), default: 1.0 \n
            - mmax: Maximum mass (Msun), default: 2.3

        Returns
        -------
        mass_1_source : ``numpy.ndarray``
            Array of primary masses in source frame (Msun).

        Examples
        --------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution(event_type='BNS')
        >>> m1_src = cbc.bimodal(size=1000)
        """

        from .prior_functions import bimodal_pdf

        identifier_dict = dict(
            param_name="mass_1_source",
            sampler_type="bimodal",
         )
        param_dict = self.available_gw_prior[identifier_dict["param_name"]][
            identifier_dict["sampler_type"]
        ].copy()
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        identifier_dict["resolution"] = self.create_new_interpolator[
            identifier_dict["param_name"]
        ]["resolution"]
        identifier_dict["cdf_size"] = self.create_new_interpolator[
            identifier_dict["param_name"]
        ]["cdf_size"]

        m_arr = np.linspace(
            identifier_dict["mmin"], identifier_dict["mmax"], identifier_dict["resolution"]
        )

        # mass function for BNS
        def pdf_(m):
            return bimodal_pdf(
                m=m,
                w=identifier_dict["w"],
                muL=identifier_dict["muL"],
                sigmaL=identifier_dict["sigmaL"],
                muR=identifier_dict["muR"],
                sigmaR=identifier_dict["sigmaR"],
                mmin=identifier_dict["mmin"],
                mmax=identifier_dict["mmax"],
            )

        mass_object = FunctionConditioning(
            function=pdf_,
            x_array=m_arr,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["sampler_type"],
            create_new=self.create_new_interpolator[identifier_dict["param_name"]][
                "create_new"
            ],
            create_function_inverse=False,
            create_function=False,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
            cdf_size=identifier_dict["cdf_size"],
        )

        return mass_object if get_attribute else mass_object.rvs(size)

    # -------------------------------------
    # Pop III BBH mass distribution models
    # -------------------------------------
    # def lognormal(
    #     self,
    #     size,
    #     get_attribute=False,
    #     **kwargs,
    # ):
    #     """
    #     Sample source masses for Population III BBH from lognormal distribution.

    #     Based on Eqn. 1 and 4 of Ng et al. 2022 for Population III black holes.

    #     Parameters
    #     ----------
    #     size : ``int``
    #         Number of samples to draw.
    #     get_attribute : ``bool``
    #         If True, return the sampler object instead of samples. \n
    #         default: False
    #     **kwargs : ``dict``
    #         Model parameters: \n
    #         - m_min: Minimum BH mass (Msun), default: 5.0 \n
    #         - m_max: Maximum BH mass (Msun), default: 150.0 \n
    #         - Mc: Central mass scale (Msun), default: 30.0 \n
    #         - sigma: Distribution width, default: 0.3

    #     Returns
    #     -------
    #     mass_1_source : ``numpy.ndarray``
    #         Array of primary masses in source frame (Msun).

    #     """

    #     from .prior_functions import binary_masses_BBH_popIII_lognormal_rvs

    #     identifier_dict = {"name": "binary_masses_BBH_popIII_lognormal"}
    #     param_dict = self.available_gw_prior["mass_1_source"][
    #         "binary_masses_BBH_popIII_lognormal"
    #     ].copy()
    #     if param_dict:
    #         param_dict.update(kwargs)
    #     else:
    #         param_dict = kwargs
    #     identifier_dict.update(param_dict)

    #     def rvs_(size):
    #         return binary_masses_BBH_popIII_lognormal_rvs(
    #             size,
    #             m_min=identifier_dict["m_min"],
    #             m_max=identifier_dict["m_max"],
    #             Mc=identifier_dict["Mc"],
    #             sigma=identifier_dict["sigma"],
    #             chunk_size=chunk_size,
    #         )

    #     mass_object = FunctionConditioning(
    #         function=None,
    #         x_array=None,
    #         identifier_dict=identifier_dict,
    #         directory=self.directory,
    #         sub_directory="mass_1_source",
    #         name=identifier_dict["name"],
    #         create_new=self.create_new_interpolator["mass_1_source"][
    #             "create_new"
    #         ],
    #         create_function_inverse=False,
    #         create_function=False,
    #         create_pdf=False,
    #         create_rvs=rvs_,
    #         callback="rvs",
    #     )

    #     return mass_object if get_attribute else mass_object.rvs(size)

    # def binary_masses_BBH_primordial_lognormal(
    #     self,
    #     size,
    #     get_attribute=False,
    #     **kwargs,
    # ):
    #     """
    #     Sample source masses for primordial BBH from lognormal distribution.

    #     Based on Eqn. 1 and 4 of Ng et al. 2022 for primordial black holes.

    #     Parameters
    #     ----------
    #     size : ``int``
    #         Number of samples to draw.
    #     get_attribute : ``bool``
    #         If True, return the sampler object instead of samples. \n
    #         default: False
    #     **kwargs : ``dict``
    #         Model parameters: \n
    #         - m_min: Minimum BH mass (Msun), default: 1.0 \n
    #         - m_max: Maximum BH mass (Msun), default: 100.0 \n
    #         - Mc: Central mass scale (Msun), default: 20.0 \n
    #         - sigma: Distribution width, default: 0.3

    #     Returns
    #     -------
    #     mass_1_source : ``numpy.ndarray``
    #         Array of primary masses in source frame (Msun).
    #     mass_2_source : ``numpy.ndarray``
    #         Array of secondary masses in source frame (Msun).
    #     """

    #     from .prior_functions import binary_masses_BBH_primordial_lognormal_rvs

    #     identifier_dict = {"name": "binary_masses_BBH_primordial_lognormal"}
    #     param_dict = self.available_gw_prior["mass_1_source"][
    #         "binary_masses_BBH_primordial_lognormal"
    #     ].copy()
    #     if param_dict:
    #         param_dict.update(kwargs)
    #     else:
    #         param_dict = kwargs
    #     identifier_dict.update(param_dict)

    #     def rvs_(size):
    #         return binary_masses_BBH_primordial_lognormal_rvs(
    #             size,
    #             m_min=identifier_dict["m_min"],
    #             m_max=identifier_dict["m_max"],
    #             Mc=identifier_dict["Mc"],
    #             sigma=identifier_dict["sigma"],
    #             chunk_size=chunk_size,
    #         )

    #     mass_object = FunctionConditioning(
    #         function=None,
    #         x_array=None,
    #         identifier_dict=identifier_dict,
    #         directory=self.directory,
    #         sub_directory="mass_1_source",
    #         name=identifier_dict["name"],
    #         create_new=self.create_new_interpolator["mass_1_source"][
    #             "create_new"
    #         ],
    #         create_function_inverse=False,
    #         create_function=False,
    #         create_pdf=False,
    #         create_rvs=rvs_,
    #         callback="rvs",
    #     )

    #     return mass_object if get_attribute else mass_object.rvs(size)

    # -------------------------------------
    # NSBH mass distribution models
    # -------------------------------------

    def broken_powerlaw(self, size, get_attribute=False, **kwargs):
        """
        Sample source masses for NSBH from broken power-law distribution.

        Uses gwcosmo-style broken power-law for black hole mass and power-law
        for neutron star mass.

        Parameters
        ----------
        size : ``int``
            Number of samples to draw.
        get_attribute : ``bool``
            If True, return the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Model parameters: \n
            - mminbh: Minimum BH mass (Msun), default: 26 \n
            - mmaxbh: Maximum BH mass (Msun), default: 125 \n
            - alpha_1: Primary power-law index, default: 6.75 \n
            - alpha_2: Secondary power-law index, default: 6.75 \n
            - b: Break point, default: 0.5 \n
            - delta_m: Tapering range (Msun), default: 5 \n

        Returns
        -------
        mass_1_source : ``numpy.ndarray``
            Array of BH masses in source frame (Msun).
        """

        from .prior_functions import broken_powerlaw_pdf

        identifier_dict = dict(
            param_name="mass_1_source",
            sampler_type="broken_powerlaw",
        )
        identifier_dict["resolution"] = self.create_new_interpolator[
            identifier_dict["param_name"]
        ]["resolution"]
        identifier_dict["cdf_size"] = self.create_new_interpolator[
            identifier_dict["param_name"]
        ]["cdf_size"]
        param_dict = self.available_gw_prior[identifier_dict["param_name"]][
            identifier_dict["sampler_type"]
        ].copy()
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        # mass function for NSBH
        def pdf_(size):
            return broken_powerlaw_pdf(
                size,
                mminbh=identifier_dict["mminbh"],
                mmaxbh=identifier_dict["mmaxbh"],
                alpha_1=identifier_dict["alpha_1"],
                alpha_2=identifier_dict["alpha_2"],
                b=identifier_dict["b"],
                delta_m=identifier_dict["delta_m"],
            )

        m1_arr = np.geomspace(
            identifier_dict["mminbh"],
            identifier_dict["mmaxbh"],
            identifier_dict["resolution"],
        )
        

        m1_object = FunctionConditioning(
            function=pdf_,
            x_array=m1_arr,
            non_negative_function=True,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["sampler_type"],
            create_new=self.create_new_interpolator[identifier_dict["param_name"]][
                "create_new"
            ],
            create_function=False,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
            cdf_size=identifier_dict["cdf_size"],
        )

        return m1_object if get_attribute else m1_object.rvs(size)
    
    # ---------------------------------------
    # General purpose samplers and functions
    # ---------------------------------------

    def truncated_normal(self, size, get_attribute=False, **kwargs):
        """
        Sample with a truncated normal distribution.

        Parameters
        ----------
        size : ``int``
            Number of samples to draw.
        get_attribute : ``bool``
            If True, return the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Model parameters: \n
            - x_min: Minimum value, default: 1.0 \n
            - x_max: Maximum value, default: 2.5 \n
            - mu: Mean of the truncated normal distribution, default: 1.4 \n
            - sigma: Standard deviation of the truncated normal distribution, default: 0.68 \n
            
        Returns
        -------
        x : ``numpy.ndarray``
            Array of values drawn from the truncated normal distribution.
        """

        from .prior_functions import truncated_normal_pdf, truncated_normal_rvs

        # Get parameters
        identifier_dict = dict(
            param_name = "unknown_parameter",
            sampler_type = "truncated_normal",
        )
        param_dict = dict(mu=1.4, sigma=0.68, x_min=1.0, x_max=2.5)

        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        identifier_dict["resolution"] = self.create_new_interpolator[identifier_dict["param_name"]][
            "resolution"
        ]

        x_min = identifier_dict["x_min"]
        x_max = identifier_dict["x_max"]
        mu = identifier_dict["mu"]
        sigma = identifier_dict["sigma"]

        @njit
        def truncated_normal_pdf_(x):
            return truncated_normal_pdf(
                x, mu=mu, sigma=sigma, x_min=x_min, x_max=x_max
            )

        @njit
        def truncated_normal_rvs_(size):
            return truncated_normal_rvs(
                size, mu=mu, sigma=sigma, x_min=x_min, x_max=x_max
            )

        object_ = FunctionConditioning(
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory=identifier_dict["param_name"],
            name=identifier_dict["sampler_type"],
            create_new=self.create_new_interpolator[identifier_dict["param_name"]][
                "create_new"
            ],
            create_function=False,
            create_pdf=truncated_normal_pdf_,
            create_rvs=truncated_normal_rvs_,
            callback='rvs',
        )

        return object_ if get_attribute else object_.rvs(size)

    def powerlaw(self, size, get_attribute=False, **kwargs):
        """
        Sample values from a normalized power-law distribution.
        p(x) ∝ x^{-alpha}, x in [x_min, x_max]

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
            - x_max: Maximum value, default: 1.0 \n
            - alpha: Power-law index, default: 2.35

        Returns
        -------
        values : ``numpy.ndarray``
            Array of values drawn from the normalized power-law distribution.
        """

        identifier_dict = dict(
            param_name = "unknown_parameter",
            sampler_type = "powerlaw",
        )
        param_dict = dict(x_min=0.0, x_max=1.0, alpha=2.35)
        param_dict.update(kwargs)
        identifier_dict.update(param_dict)
        
        x_min = identifier_dict["x_min"]
        x_max = identifier_dict["x_max"]
        alpha = identifier_dict["alpha"]

        from .prior_functions import powerlaw_pdf, powerlaw_rvs

        @njit
        def pdf_(x):
            return powerlaw_pdf(x, x_min=x_min, x_max=x_max, alpha=alpha)

        @njit
        def rvs_(size):
            return powerlaw_rvs(size, x_min=x_min, x_max=x_max, alpha=alpha)   

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

    def sampler_cosine(self, size, get_attribute=False, **kwargs):
        """
        Sample from cosine distribution for declination.

        Samples values in range [-pi/2, pi/2] following a cosine distribution,
        appropriate for isotropic sky position declination.

        Parameters
        ----------
        size : ``int``
            Number of samples to draw.
        get_attribute : ``bool``
            If True, return the sampler object instead of samples. \n
            default: False

        Returns
        -------
        values : ``numpy.ndarray``
            Array of values in range [-pi/2, pi/2] (rad).
        """

        identifier_dict = dict(
            param_name = "unknown_parameter",
            sampler_type = "sampler_cosine",
        )
        identifier_dict.update(kwargs)

        x_min = identifier_dict.get("x_min", -np.pi/2)
        x_max = identifier_dict.get("x_max", np.pi/2)

        @njit
        def pdf_(x):
            # return 0.5 * np.cos(x)
            return np.where((x >= x_min) & (x <= x_max), 0.5 * np.cos(x), 0.0)

        @njit
        def rvs_(size):
            return np.arcsin((np.random.uniform(0, 1, size=size) * 2 - 1))

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
            return object_(size)

    def sampler_sine(self, size, get_attribute=False, **kwargs):
        """
        Sample from sine distribution for inclination angles.

        Samples values in range [0, pi] following a sine distribution,
        appropriate for isotropic orientation angles.

        Parameters
        ----------
        size : ``int``
            Number of samples to draw.
        get_attribute : ``bool``
            If True, return the sampler object instead of samples. \n
            default: False

        Returns
        -------
        values : ``numpy.ndarray``
            Array of values in range [0, pi] (rad).
        """

        identifier_dict = dict(
            param_name = "unknown_parameter",
            sampler_type = "sampler_sine",
        )
        identifier_dict.update(kwargs)
        x_min = identifier_dict.get("x_min", 0.0)
        x_max = identifier_dict.get("x_max", np.pi)

        @njit
        def pdf_(x):
            return np.where((x >= x_min) & (x <= x_max), 0.5 * np.sin(x), 0.0)

        @njit
        def rvs_(size):
            return np.arccos((np.random.uniform(0, 1, size=size) - 0.5) * 2)

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
            return object_(size)

    @property
    def mass_1_source(self):
        """
        Class object (of FunctionConditioning) for source frame primary mass, with rvs/sampler as callback. Can also be a user defined callable sampler. \n

        The class object contains the following attribute methods: \n
        - `function(m1)`: returns the primary mass distribution function in terms of the merger rate at fixed redshift. \n
        - `pdf(m1)`: returns the probability density function of the primary mass distribution. \n
        - `rvs(size)`: returns random samples from the primary mass distribution.

        Returns
        -------
        mass_1_source : ``numpy.ndarray``
            Array of mass_1_source values in solar masses.

        Examples
        --------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc_source_param_dist = CBCSourceParameterDistribution()
        >>> cbc_source_param_dist.mass_1_source(size=10)
        """

        return self._mass_1_source

    @mass_1_source.setter
    def mass_1_source(self, prior):
        """
        Prior be a string corresponding to a built-in prior in `available_gw_prior`, or a user defined callable function or a FunctionConditioning object.
        """

        if prior in self.available_gw_prior["mass_1_source"]:
            print(f"using ler available mass_1_source function : {prior}")
            args = self.gw_param_samplers_params["mass_1_source"]
            if args is None:
                self._mass_1_source = getattr(self, prior)(
                    size=None, get_attribute=True
                )
            else:
                # following should return a sampler function with only one argument (size)
                self._mass_1_source = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user defined custom mass_1_source class/object of type ler.utils.FunctionConditioning"
            )
            self._mass_1_source = prior
        elif callable(prior):
            print("using user defined custom mass_1_source function")
            self._mass_1_source = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior, callback='rvs'
            )
        else:
            raise ValueError(
                "mass_1_source prior not available in available_gw_prior. Must be a string or a callable function."
            )

    @property
    def mass_ratio(self):
        """
        Class object (of FunctionConditioning) for mass ratio, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        -  `rvs(size)` or `rvs(size, m1)`: returns random samples from the mass ratio distribution \n
        - `pdf(q)` or `pdf(q, m1)`: returns the probability density function of the mass ratio distribution \n

        Returns
        -------
        mass_ratio : ``numpy.ndarray``
            Array of mass ratio values.
        """

        return self._mass_ratio

    @mass_ratio.setter
    def mass_ratio(self, prior):
        if prior in self.available_gw_prior["mass_ratio"]:
            print(f"using ler available mass_ratio function : {prior}")
            args = self.gw_param_samplers_params["mass_ratio"]
            if args is None:
                self._mass_ratio = getattr(self, prior)(size=None, get_attribute=True)
            else:
                # following should return a sampler function with only one argument (size)
                self._mass_ratio = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user defined custom mass_ratio class/object of type ler.utils.FunctionConditioning"
            )
            self._mass_ratio = prior
        elif callable(prior):
            print("using user defined custom mass_ratio function")
            self._mass_ratio = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior, callback='rvs'
            )
        elif prior is None:
            print("No mass_ratio prior provided. mass_ratio = mass_2_source / mass_1_source")
            self._mass_ratio = None
        else:
            raise ValueError(
                "mass_ratio prior not available in available_gw_prior. Must be a string or a callable function."
            )

    @property
    def mass_2_source(self):
        """
        Class object (of FunctionConditioning) for source frame secondary mass, with rvs/sampler as callback. Can also be a user defined callable sampler. \n

        The class object contains the following attribute methods: \n
        - `function(m1)`: returns the secondary mass distribution function in terms of the merger rate at fixed redshift. \n
        - `pdf(m1)`: returns the probability density function of the secondary mass distribution. \n
        - `rvs(size)`: returns random samples from the secondary mass distribution.

        Returns
        -------
        mass_2_source : ``numpy.ndarray``
            Array of mass_2_source values in solar masses.

        Examples
        --------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc_source_param_dist = CBCSourceParameterDistribution()
        >>> cbc_source_param_dist.mass_2_source(size=10)
        """

        return self._mass_2_source

    @mass_2_source.setter
    def mass_2_source(self, prior):
        """
        Prior be a string corresponding to a built-in prior in `available_gw_prior`, or a user defined callable function or a FunctionConditioning object.
        """

        if prior in self.available_gw_prior["mass_2_source"]:
            print(f"using ler available mass_2_source function : {prior}")
            args = self.gw_param_samplers_params["mass_2_source"]
            if args is None:
                self._mass_2_source = getattr(self, prior)(
                    size=None, get_attribute=True
                )
            else:
                # following should return a sampler function with only one argument (size)
                self._mass_2_source = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user defined custom mass_2_source class/object of type ler.utils.FunctionConditioning"
            )
            self._mass_2_source = prior
        elif callable(prior):
            print("using user defined custom mass_2_source function")
            self._mass_2_source = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior, callback='rvs'
            )
        elif prior is None:
            print("No mass_2_source prior provided. mass_2_source = mass_1_source * mass_ratio")
            self._mass_2_source = None
        else:
            raise ValueError(
                "mass_2_source prior not available in available_gw_prior. Must be a string or a callable function."
            )

    @property
    def zs(self):
        """
        Class object (of FunctionConditioning) for source redshift, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the redshift distribution \n
        - `pdf`: returns the probability density function of the redshift distribution \n
        - `function`: returns the redshift distribution function.

        Returns
        -------
        zs : ``numpy.ndarray``
            Array of redshift values.
        """

        return self._zs

    @zs.setter
    def zs(self, prior):
        if prior in self.available_gw_prior["zs"]:
            print(f"using ler available zs function : {prior}")
            self._zs = getattr(self, prior)
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user defined custom zs class/object of type ler.utils.FunctionConditioning"
            )
            self._zs = prior
        elif callable(prior):
            print("using user defined custom zs function")
            self._zs = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior, callback='rvs'
            )
        else:
            raise ValueError(
                "zs prior not available in available_gw_prior. Must be a string or a callable function."
            )

    @property
    def geocent_time(self):
        """
        Class object (of FunctionConditioning) for geocentric time, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the geocentric time distribution \n
        - `pdf`: returns the probability density function of the geocentric time distribution \n
        - `function`: returns the geocentric time distribution function.

        Returns
        -------
        geocent_time : ``numpy.ndarray``
            Array of geocentric time values.
        """

        return self._geocent_time

    @geocent_time.setter
    def geocent_time(self, prior):
        if prior in self.available_gw_prior["geocent_time"]:
            print(f"using ler available geocent_time function : {prior}")
            args = self.gw_param_samplers_params["geocent_time"]
            if args is None:
                self._geocent_time = getattr(self, prior)(size=None, get_attribute=True)
            else:
                # following should return a sampler function with only one argument (size)
                self._geocent_time = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user defined custom geocent_time class/object of type ler.utils.FunctionConditioning"
            )
            self._geocent_time = prior
        elif callable(prior):
            print("using user defined custom geocent_time function")
            self._geocent_time = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior, callback='rvs'
            )
        else:
            raise ValueError(
                "geocent_time prior not available in available_gw_prior. Must be a string or a callable function."
            )

    @property
    def ra(self):
        """
        Class object (of FunctionConditioning) for right ascension, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the right ascension distribution \n
        - `pdf`: returns the probability density function of the right ascension distribution \n
        - `function`: returns the right ascension distribution function.

        Returns
        -------
        ra : ``numpy.ndarray``
            Array of right ascension values.
        """

        return self._ra

    @ra.setter
    def ra(self, prior):
        if prior in self.available_gw_prior["ra"]:
            print(f"using ler available ra function : {prior}")
            args = self.gw_param_samplers_params["ra"]
            if args is None:
                self._ra = getattr(self, prior)(size=None, get_attribute=True)
            else:
                # following should return a sampler function with only one argument (size)
                self._ra = getattr(self, prior)(size=None, get_attribute=True, **args)
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user defined custom ra class/object of type ler.utils.FunctionConditioning"
            )
            self._ra = prior
        elif callable(prior):
            print("using user defined custom ra function")
            self._ra = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior, callback='rvs'
            )
        else:
            raise ValueError(
                "ra prior not available in available_gw_prior. Must be a string or a callable function."
            )

    @property
    def dec(self):
        """
        Class object (of FunctionConditioning) for declination, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the declination distribution \n
        - `pdf`: returns the probability density function of the declination distribution \n
        - `function`: returns the declination distribution function.

        Returns
        -------
        dec : ``numpy.ndarray``
            Array of declination values.
        """

        return self._dec

    @dec.setter
    def dec(self, prior):
        if prior in self.available_gw_prior["dec"]:
            print(f"using ler available dec function : {prior}")
            args = self.gw_param_samplers_params["dec"]
            if args is None:
                self._dec = getattr(self, prior)(size=None, get_attribute=True)
            else:
                # following should return a sampler function with only one argument (size)
                self._dec = getattr(self, prior)(size=None, get_attribute=True, **args)
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user defined custom dec class/object of type ler.utils.FunctionConditioning"
            )
            self._dec = prior
        elif callable(prior):
            print("using user defined custom dec function")
            self._dec = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior, callback='rvs'
            )
        else:
            raise ValueError(
                "dec prior not available in available_gw_prior. Must be a string or a callable function."
            )

    @property
    def phase(self):
        """
        Class object (of FunctionConditioning) for coalescence phase, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the coalescence phase distribution \n
        - `pdf`: returns the probability density function of the coalescence phase distribution \n
        - `function`: returns the coalescence phase distribution function.

        Returns
        -------
        phase : ``numpy.ndarray``
            Array of coalescence phase values.
        """

        return self._phase

    @phase.setter
    def phase(self, prior):
        if prior in self.available_gw_prior["phase"]:
            print(f"using ler available phase function : {prior}")
            args = self.gw_param_samplers_params["phase"]
            if args is None:
                self._phase = getattr(self, prior)(size=None, get_attribute=True)
            else:
                # following should return a sampler function with only one argument (size)
                self._phase = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user defined custom phase class/object of type ler.utils.FunctionConditioning"
            )
            self._phase = prior
        elif callable(prior):
            print("using user defined custom phase function")
            self._phase = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior, callback='rvs'
            )
        else:
            raise ValueError(
                "phase prior not available in available_gw_prior. Must be a string or a callable function."
            )

    @property
    def psi(self):
        """
        Class object (of FunctionConditioning) for polarization angle, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the polarization angle distribution \n
        - `pdf`: returns the probability density function of the polarization angle distribution \n
        - `function`: returns the polarization angle distribution function.

        Returns
        -------
        psi : ``numpy.ndarray``
            Array of polarization angle values.
        """

        return self._psi

    @psi.setter
    def psi(self, prior):
        if prior in self.available_gw_prior["psi"]:
            print(f"using ler available psi function : {prior}")
            args = self.gw_param_samplers_params["psi"]
            if args is None:
                self._psi = getattr(self, prior)(size=None, get_attribute=True)
            else:
                # following should return a sampler function with only one argument (size)
                self._psi = getattr(self, prior)(size=None, get_attribute=True, **args)
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user defined custom psi class/object of type ler.utils.FunctionConditioning"
            )
            self._psi = prior
        elif callable(prior):
            print("using user defined custom psi function")
            self._psi = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior, callback='rvs'
            )
        else:
            raise ValueError(
                "psi prior not available in available_gw_prior. Must be a string or a callable function."
            )

    @property
    def theta_jn(self):
        """
        Class object (of FunctionConditioning) for inclination angle, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the inclination angle distribution \n
        - `pdf`: returns the probability density function of the inclination angle distribution \n
        - `function`: returns the inclination angle distribution function.

        Returns
        -------
        theta_jn : ``numpy.ndarray``
            Array of inclination angle values, i.e. the angle between the line of sight and the orbital angular momentum (rad).
        """

        return self._theta_jn

    @theta_jn.setter
    def theta_jn(self, prior):
        if prior in self.available_gw_prior["theta_jn"]:
            print(f"using ler available theta_jn function : {prior}")
            args = self.gw_param_samplers_params["theta_jn"]
            if args is None:
                self._theta_jn = getattr(self, prior)(size=None, get_attribute=True)
            else:
                # following should return a sampler function with only one argument (size)
                self._theta_jn = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user defined custom theta_jn class/object of type ler.utils.FunctionConditioning"
            )
            self._theta_jn = prior
        elif callable(prior):
            print("using user defined custom theta_jn function")
            self._theta_jn = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior, callback='rvs'
            )
        else:
            raise ValueError(
                "theta_jn prior not available in available_gw_prior. Must be a string or a callable function."
            )

    @property
    def a_1(self):
        """
        Class object (of FunctionConditioning) for spin1 magnitude, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the spin1 magnitude distribution \n
        - `pdf`: returns the probability density function of the spin1 magnitude distribution \n
        - `function`: returns the spin1 magnitude distribution function.\n
        - `None` : if spin_zero=True.

        Returns
        -------
        a_1 : ``numpy.ndarray``
            Array of spin magnitude values for the primary body.
        """

        return self._a_1

    @a_1.setter
    def a_1(self, prior):
        if prior in self.available_gw_prior["a_1"]:
            print(f"using ler available a_1 function : {prior}")
            args = self.gw_param_samplers_params["a_1"]
            if args is None:
                self._a_1 = getattr(self, prior)(size=None, get_attribute=True)
            else:
                # following should return a sampler function with only one argument (size)
                self._a_1 = getattr(self, prior)(size=None, get_attribute=True, **args)
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user defined custom a_1 class/object of type ler.utils.FunctionConditioning"
            )
            self._a_1 = prior
        elif callable(prior):
            print("using user defined custom a_1 function")
            self._a_1 = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior, callback='rvs'
            )
        elif prior is None:
            print("No a_1 prior provided. a_1 prior will be set to None")
            self._a_1 = prior
        else:
            raise ValueError(
                "a_1 prior not available in available_gw_prior. Must be a string or a callable function."
            )

    @property
    def a_2(self):
        """
        Class object (of FunctionConditioning) for spin2 magnitude, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the spin2 magnitude distribution \n
        - `pdf`: returns the probability density function of the spin2 magnitude distribution \n
        - `function`: returns the spin2 magnitude distribution function.
        - `None` : if spin_zero=True.

        Returns
        -------
        a_2 : ``numpy.ndarray``
            Array of spin magnitude values for the secondary body.
        """

        return self._a_2

    @a_2.setter
    def a_2(self, prior):
        if prior in self.available_gw_prior["a_2"]:
            print(f"using ler available a_2 function : {prior}")
            args = self.gw_param_samplers_params["a_2"]
            if args is None:
                self._a_2 = getattr(self, prior)(size=None, get_attribute=True)
            else:
                # following should return a sampler function with only one argument (size)
                self._a_2 = getattr(self, prior)(size=None, get_attribute=True, **args)
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user defined custom a_2 class/object of type ler.utils.FunctionConditioning"
            )
            self._a_2 = prior
        elif callable(prior):
            print("using user defined custom a_2 function")
            self._a_2 = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior, callback='rvs'
            )
        elif prior is None:
            print("No a_2 prior provided. a_2 prior will be set to None")
            self._a_2 = prior
        else:
            raise ValueError(
                "a_2 prior not available in available_gw_prior. Must be a string or a callable function."
            )

    @property
    def tilt_1(self):
        """
        Class object (of FunctionConditioning) for tilt1 angle, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the tilt1 angle distribution \n
        - `pdf`: returns the probability density function of the tilt1 angle distribution \n
        - `function`: returns the tilt1 angle distribution function.
        - `None` : if spin_zero=True or spin_precession=False.

        Returns
        -------
        tilt_1 : ``numpy.ndarray``
            Array of the spin tilt angle of the primary body, i.e. the angle between the spin vector and the orbital angular momentum for the primary body (rad).
        """

        return self._tilt_1

    @tilt_1.setter
    def tilt_1(self, prior):
        if prior in self.available_gw_prior["tilt_1"]:
            print(f"using ler available tilt_1 function : {prior}")
            args = self.gw_param_samplers_params["tilt_1"]
            if args is None:
                self._tilt_1 = getattr(self, prior)(size=None, get_attribute=True)
            else:
                self._tilt_1 = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user defined custom tilt_1 class/object of type ler.utils.FunctionConditioning"
            )
            self._tilt_1 = prior
        elif callable(prior):
            print("using user defined custom tilt_1 function")
            self._tilt_1 = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior, callback='rvs'
            )
        elif prior is None:
            print("No tilt_1 prior provided. tilt_1 prior will be set to None")
            self._tilt_1 = prior
        else:
            raise ValueError(
                "tilt_1 prior not available in available_gw_prior. Must be a string or a callable function."
            )

    @property
    def tilt_2(self):
        """
        Class object (of FunctionConditioning) for tilt2 angle, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the tilt2 angle distribution \n
        - `pdf`: returns the probability density function of the tilt2 angle distribution \n
        - `function`: returns the tilt2 angle distribution function.
        - `None` : if spin_zero=True or spin_precession=False.

        Returns
        -------
        tilt_2 : ``numpy.ndarray``
            Array of the spin tilt angle of the secondary body, i.e. the angle between the spin vector and the orbital angular momentum for the secondary body (rad).
        """

        return self._tilt_2

    @tilt_2.setter
    def tilt_2(self, prior):
        if prior in self.available_gw_prior["tilt_2"]:
            print(f"using ler available tilt_2 function : {prior}")
            args = self.gw_param_samplers_params["tilt_2"]
            if args is None:
                self._tilt_2 = getattr(self, prior)(size=None, get_attribute=True)
            else:
                self._tilt_2 = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user defined custom tilt_2 class/object of type ler.utils.FunctionConditioning"
            )
            self._tilt_2 = prior
        elif callable(prior):
            print("using user defined custom tilt_2 function")
            self._tilt_2 = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior, callback='rvs'
            )
        elif prior is None:
            print("No tilt_2 prior provided. tilt_2 prior will be set to None")
            self._tilt_2 = prior
        else:
            raise ValueError(
                "tilt_2 prior not available in available_gw_prior. Must be a string or a callable function."
            )

    @property
    def phi_12(self):
        """
        Class object (of FunctionConditioning) for phi_12 angle, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the phi_12 angle distribution \n
        - `pdf`: returns the probability density function of the phi_12 angle distribution \n
        - `function`: returns the phi_12 angle distribution function.
        - `None` : if spin_zero=True or spin_precession=False.

        Returns
        -------
        phi_12 : ``numpy.ndarray``
            Array of the spin tilt angle between the two spins, i.e., angle between the projections of the two spins onto the orbital plane (rad).
        """

        return self._phi_12

    @phi_12.setter
    def phi_12(self, prior):
        if prior in self.available_gw_prior["phi_12"]:
            print(f"using ler available phi_12 function : {prior}")
            args = self.gw_param_samplers_params["phi_12"]
            if args is None:
                self._phi_12 = getattr(self, prior)(size=None, get_attribute=True)
            else:
                self._phi_12 = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user defined custom phi_12 class/object of type ler.utils.FunctionConditioning"
            )
            self._phi_12 = prior
        elif callable(prior):
            print("using user defined custom phi_12 function")
            self._phi_12 = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior, callback='rvs'
            )
        elif prior is None:
            print("No phi_12 prior provided. phi_12 prior will be set to None")
            self._phi_12 = prior
        else:
            raise ValueError(
                "phi_12 prior not available in available_gw_prior. Must be a string or a callable function."
            )

    @property
    def phi_jl(self):
        """
        Class object (of FunctionConditioning) for phi_jl angle, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the phi_jl angle distribution \n
        - `pdf`: returns the probability density function of the phi_jl angle distribution \n
        - `function`: returns the phi_jl angle distribution function.
        - `None` : if spin_zero=True or spin_precession=False.

        Returns
        -------
        phi_jl : ``numpy.ndarray``
            Array of the angle values between the orientation of the total angular momentum around the orbital angular momentum (rad).
        """
        return self._phi_jl

    @phi_jl.setter
    def phi_jl(self, prior):
        if prior in self.available_gw_prior["phi_jl"]:
            print(f"using ler available phi_jl function : {prior}")
            args = self.gw_param_samplers_params["phi_jl"]
            if args is None:
                self._phi_jl = getattr(self, prior)(size=None, get_attribute=True)
            else:
                self._phi_jl = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user defined custom phi_jl class/object of type ler.utils.FunctionConditioning"
            )
            self._phi_jl = prior
        elif callable(prior):
            print("using user defined custom phi_jl function")
            self._phi_jl = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior, callback='rvs'
            )
        elif prior is None:
            print("No phi_jl prior provided. phi_jl prior will be set to None")
            self._phi_jl = prior
        else:
            raise ValueError(
                "phi_jl prior not available in available_gw_prior. Must be a string or a callable function."
            )

    @property
    def available_gw_prior(self):
        """
        Dictionary of all available prior distributions and their parameters. \n
        This is a dynamically generated dictionary containing available samplers
        for each GW parameter type and their default parameter values.

        Returns
        -------
        available_gw_prior : ``dict``
            Nested dictionary organized by parameter type (e.g., 'mass_1_source', \n
            'geocent_time', etc.) with sampler names and default parameters.
        """

        self._available_gw_prior = dict(
            zs=dict(
                merger_rate_density_based_source_redshift=None,
            ),
            mass_1_source=dict(
                # BBH
                broken_powerlaw_plus_2peaks=dict(
                    param_name="mass_1_source",
                    sampler_type="broken_powerlaw_plus_2peaks",
                    lam_0=0.361,
                    lam_1=0.586,
                    mpp_1=9.764,
                    sigpp_1=0.649,
                    mpp_2=32.763,
                    sigpp_2=3.918,
                    mlow_1=5.059,
                    delta_m_1=4.321,
                    break_mass=35.622,
                    alpha_1=1.728,
                    alpha_2=4.512,
                    mmax=300.0,
                    normalization_size=500,
                ),
                powerlaw_plus_peak=dict(
                    param_name="mass_1_source",
                    sampler_type="powerlaw_plus_peak",
                    mminbh=4.98,
                    mmaxbh=112.5,
                    alpha=3.78,
                    mu_g=32.27,
                    sigma_g=3.88,
                    lambda_peak=0.03,
                    delta_m=4.8,
                    normalization_size=500,
                ),
                # BNS
                truncated_normal=dict(
                    param_name="mass_1_source",
                    sampler_type="truncated_normal",
                    mu=1.4, sigma=0.68, x_min=1.0, x_max=2.5
                ),
                uniform=dict(
                    param_name="mass_1_source",
                    sampler_type="uniform",
                    x_min=1.0, x_max=2.5),
                powerlaw=dict(
                    param_name="mass_1_source",
                    sampler_type="powerlaw",
                    x_min=1.0, x_max=2.5, alpha=7.7
                ),
                bimodal=dict(
                    param_name="mass_1_source",
                    sampler_type="bimodal",
                    w=0.643,
                    muL=1.352,
                    sigmaL=0.08,
                    muR=1.88,
                    sigmaR=0.3,
                    mmin=1.0,
                    mmax=2.3,
                ),
                # binary_masses_BBH_popIII_lognormal=dict(
                #     m_min=5.0, m_max=150.0, Mc=30.0, sigma=0.3
                # ),
                # binary_masses_BBH_primordial_lognormal=dict(
                #     m_min=1.0, m_max=100.0, Mc=20.0, sigma=0.3
                # ),
                # NSBH
                broken_powerlaw=dict(
                    mminbh=26,
                    mmaxbh=125,
                    alpha_1=6.75,
                    alpha_2=6.75,
                    b=0.5,
                    delta_m=5,
                    # mminns=1.0,
                    # mmaxns=3.0,
                    # alphans=0.0,
                ),
            ),
            mass_ratio=dict(
                powerlaw_with_smoothing=dict(
                    param_name="mass_ratio",
                    sampler_type="powerlaw_with_smoothing",
                    q_min=0.01, 
                    q_max=1.0,
                    mlow_2=3.551, 
                    mmax=300.0, 
                    beta=1.171, 
                    delta_m=4.910
                ),
            ),
            mass_2_source=dict(
                truncated_normal=dict(
                    param_name="mass_2_source",
                    sampler_type="truncated_normal",
                    mu=1.4, sigma=0.68, x_min=1.0, x_max=2.5
                ),
                uniform=dict(
                    param_name="mass_2_source",
                    sampler_type="uniform",
                    x_min=1.0, x_max=2.5),
                powerlaw=dict(
                    param_name="mass_2_source",
                    sampler_type="powerlaw",
                    x_min=1.0, x_max=2.5, alpha=7.7
                ),
                bimodal=dict(
                    param_name="mass_2_source",
                    sampler_type="bimodal",
                    w=0.643,
                    muL=1.352,
                    sigmaL=0.08,
                    muR=1.88,
                    sigmaR=0.3,
                    mmin=1.0,
                    mmax=2.3,
                ),
                # NSBH
                # powerlaw=dict(
                #     param_name="mass_ratio",
                #     sampler_type="powerlaw",
                #     x_min=1.0, x_max=3.0, alpha=0.0
                # ),
            ),
            a_1=(
                dict(
                    constant_values_n_size=dict(
                        param_name="a_1", sampler_type="constant_values_n_size",
                        value=0.0
                    ),
                    uniform=dict(param_name="a_1", sampler_type="uniform", x_min=0.0, x_max=0.8),
                    truncated_normal=dict(param_name="a_1", sampler_type="truncated_normal", x_min=0.0, x_max=1.0, mu=0.085, sigma=0.330),
                )
            ),
            a_2=(
                dict(
                    constant_values_n_size=dict(
                        param_name="a_2", sampler_type="constant_values_n_size",
                        value=0.0
                    ),
                    uniform=dict(param_name="a_2", sampler_type="uniform", x_min=0.0, x_max=0.8),
                    truncated_normal=dict(param_name="a_2", sampler_type="truncated_normal", x_min=0.0, x_max=1.0, mu=0.085, sigma=0.330),
                )
            ),
            tilt_1=dict(
                constant_values_n_size=dict(param_name="tilt_1", sampler_type="constant_values_n_size", value=0.0),
                sampler_sine=dict(param_name="tilt_1", sampler_type="sampler_sine", tilt_1_min=0.0, tilt_1_max=np.pi),
                gaussian_plus_isotropic=dict(param_name="tilt_1", sampler_type="gaussian_plus_isotropic", tilt_1_min=0.0, tilt_1_max=np.pi, mu_t=0.426, sigma_t=1.222, zeta=0.652),
            ),
            tilt_2=dict(
                constant_values_n_size=dict(param_name="tilt_2", sampler_type="constant_values_n_size", value=0.0),
                sampler_sine=dict(param_name="tilt_2", sampler_type="sampler_sine", tilt_2_min=0.0, tilt_2_max=np.pi),
                gaussian_plus_isotropic_joint=dict(param_name="tilt_2", sampler_type="gaussian_plus_isotropic_joint", tilt_2_min=0.0, tilt_2_max=np.pi, mu_t=0.426, sigma_t=1.222, zeta=0.652),
            ),
            phi_12=dict(
                constant_values_n_size=dict(param_name="phi_12", sampler_type="constant_values_n_size", value=0.0),
                uniform=dict(param_name="phi_12", sampler_type="uniform", x_min=0.0, x_max=2 * np.pi),
            ),
            phi_jl=dict(
                constant_values_n_size=dict(param_name="phi_jl", sampler_type="constant_values_n_size", value=0.0),
                uniform=dict(param_name="phi_jl", sampler_type="uniform", x_min=0.0, x_max=2 * np.pi),
            ),
            geocent_time=dict(
                uniform=dict(param_name="geocent_time", sampler_type="uniform", x_min=1238166018, x_max=1238166018 + 31557600.0),
                constant_values_n_size=dict(param_name="geocent_time", sampler_type="constant_values_n_size", value=1238166018),
            ),
            ra=dict(
                uniform=dict(param_name="ra", sampler_type="uniform", x_min=0.0, x_max=2 * np.pi),
                constant_values_n_size=dict(param_name="ra", sampler_type="constant_values_n_size", value=0.0),
            ),
            dec=dict(
                sampler_cosine=dict(param_name="dec", sampler_type="sampler_cosine", x_min=-np.pi / 2, x_max=np.pi / 2),
                constant_values_n_size=dict(param_name="dec", sampler_type="constant_values_n_size", value=0.0),
                uniform=dict(param_name="dec", sampler_type="uniform", x_min=-np.pi / 2, x_max=np.pi / 2),
            ),
            phase=dict(
                uniform=dict(param_name="phase", sampler_type="uniform", x_min=0.0, x_max=2 * np.pi),
                constant_values_n_size=dict(param_name="phase", sampler_type="constant_values_n_size", value=0.0),
            ),
            psi=dict(
                uniform=dict(param_name="psi", sampler_type="uniform", x_min=0.0, x_max=np.pi),
                constant_values_n_size=dict(param_name="psi", sampler_type="constant_values_n_size", value=0.0),
            ),
            theta_jn=dict(
                sampler_sine=dict(param_name="theta_jn", sampler_type="sampler_sine", x_min=0.0, x_max=np.pi),
                constant_values_n_size=dict(param_name="theta_jn", sampler_type="constant_values_n_size", value=0.0),
                uniform=dict(param_name="theta_jn", sampler_type="uniform", x_min=0.0, x_max=np.pi),
            ),
        )

        return self._available_gw_prior

    @property
    def available_gw_functions(self):
        """
        Dictionary of available gravitational wave related functions and their default parameters.

        Returns
        -------
        available_gw_functions : ``dict``
            Dictionary with function names and default parameters.
        """

        self._available_gw_functions = dict(
            merger_rate_density=self.available_merger_rate_density_model,
            param_sampler_type=dict(
                gw_parameters_rvs=None,
                gw_parameters_rvs_njit=None,
            ),
        )

        return self._available_gw_functions
