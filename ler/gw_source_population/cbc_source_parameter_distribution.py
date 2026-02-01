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
    >>> params = cbc.sample_gw_parameters(size=1000)
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

chunk_size = 10000 # chunk size for rejection sampling


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
    source_priors : ``dict`` or ``None``
        Dictionary of prior sampler functions for each parameter. \n
        If None, uses default priors based on event_type. \n
        default: None
    source_priors_params : ``dict`` or ``None``
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
    >>> params = cbc.sample_gw_parameters(size=1000)
    >>> print(list(params.keys()))

    Instance Methods
    ----------
    CBCSourceParameterDistribution has the following methods: \n
    +-----------------------------------------------------+------------------------------------------------+
    | Method                                              | Description                                    |
    +=====================================================+================================================+
    | :meth:`~sample_gw_parameters`                       | Sample all GW parameters for compact binaries  |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~binary_masses_BBH_powerlaw_gaussian`| Sample BBH masses with PowerLaw+PEAK model     |
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
    | :meth:`~sampler_uniform`                            | Sample from uniform distribution               |
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
    | :attr:`~source_frame_masses`                   | ``callable``           |       | Sampler for source frame masses                |
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

    source_priors = None
    """``dict`` \n
    Dictionary of prior sampler functions.
    """

    source_priors_params = None
    """``dict`` \n
    Dictionary of prior sampler functions' input parameters.
    """

    cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)
    """``astropy.cosmology`` \n
    Cosmology to use.
    """

    spin_zero = None
    """``bool`` \n
    If True, spin prior is set to zero.
    """

    def __init__(
        self,
        z_min=0.0,
        z_max=10.0,
        event_type="BBH",
        source_priors=None,
        source_priors_params=None,
        cosmology=None,
        spin_zero=False,
        spin_precession=False,
        directory="./interpolator_json",
        create_new_interpolator=False,
    ):
        # set attributes
        self.z_min = z_min
        self.z_max = z_max
        self.cosmo = cosmology if cosmology else LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0)
        # note that self.cosmo is initialized in the super class
        self.spin_zero = spin_zero
        self.spin_precession = spin_precession
        self.directory = directory

        # setting up the interpolator creation parameters
        create_new_interpolator = self._setup_decision_dictionary_gw_params(create_new_interpolator)

        # dealing with prior functions and categorization
        (
            self.gw_param_samplers,
            self.gw_param_samplers_params,
        ) = self._source_priors_categorization(
            event_type, source_priors, source_priors_params
        )

        # initialize the SourceGalaxyPopulationModel mother class
        # for redshift distribution
        # instance attribute source_redshift is initialized here
        super().__init__(
            z_min=z_min,
            z_max=z_max,
            event_type=event_type,
            merger_rate_density=self.gw_param_samplers["merger_rate_density"],
            merger_rate_density_param=self.gw_param_samplers_params[
                "merger_rate_density"
            ],
            cosmology=cosmology,
            directory=directory,
            create_new_interpolator=create_new_interpolator,
        )

        print("\nInitializing CBCSourceParameterDistribution class...\n")
        # initializing samplers
        # it goes through the setter functions and assign the sampler functions
        # self.source_redshift is already initialized in the super class
        self.zs = self.gw_param_samplers["zs"]
        self.source_frame_masses = self.gw_param_samplers["source_frame_masses"]
        self.geocent_time = self.gw_param_samplers["geocent_time"]
        self.ra = self.gw_param_samplers["ra"]
        self.dec = self.gw_param_samplers["dec"]
        self.phase = self.gw_param_samplers["phase"]
        self.psi = self.gw_param_samplers["psi"]
        self.theta_jn = self.gw_param_samplers["theta_jn"]
        # initialize the spin prior attribute
        # remove spin prior if spin_zero is True
        if not spin_zero:
            self.a_1 = self.gw_param_samplers["a_1"]
            self.a_2 = self.gw_param_samplers["a_2"]
            if spin_precession:
                self.tilt_1 = self.gw_param_samplers["tilt_1"]
                self.tilt_2 = self.gw_param_samplers["tilt_2"]
                self.phi_12 = self.gw_param_samplers["phi_12"]
                self.phi_jl = self.gw_param_samplers["phi_jl"]


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
            source_frame_masses=dict(create_new=False, resolution=500),
            geocent_time=dict(create_new=False, resolution=500),
            ra=dict(create_new=False, resolution=500),
            dec=dict(create_new=False, resolution=500),
            phase=dict(create_new=False, resolution=500),
            psi=dict(create_new=False, resolution=500),
            theta_jn=dict(create_new=False, resolution=500),
            a_1=dict(create_new=False, resolution=500),
            a_2=dict(create_new=False, resolution=500),
            tilt_1=dict(create_new=False, resolution=500),
            tilt_2=dict(create_new=False, resolution=500),
            phi_12=dict(create_new=False, resolution=500),
            phi_jl=dict(create_new=False, resolution=500),
            merger_rate_density=dict(create_new=False, resolution=500),
            redshift_distribution=dict(create_new=False, resolution=500), 
            luminosity_distance=dict(create_new=False, resolution=500), 
            differential_comoving_volume=dict(create_new=False, resolution=500),
        )
        if isinstance(create_new_interpolator, dict):
            create_new_interpolator_.update(create_new_interpolator)
        elif create_new_interpolator is True:
            for key in create_new_interpolator_:
                create_new_interpolator_[key]["create_new"] = True

        return create_new_interpolator_

    def _source_priors_categorization(
        self, event_type, source_priors, source_prior_params
    ):
        """
        Helper function to categorize event priors and parameters.

        Sets up default prior distributions based on event type and merges
        with any user-provided priors.

        Parameters
        ----------
        event_type : ``str``
            Type of compact binary event. \n
            Options: 'BBH', 'BNS', 'NSBH', 'BBH_popIII', 'BBH_primordial'
        source_priors : ``dict`` or ``None``
            User-provided prior sampler functions for each parameter.
        source_prior_params : ``dict`` or ``None``
            User-provided parameters for each prior sampler.

        Returns
        -------
        source_priors_ : ``dict``
            Dictionary of prior sampler functions for each parameter.
        source_prior_params_ : ``dict``
            Dictionary of sampler parameters for each parameter.
        """

        # for BBH
        if event_type == "BBH":
            merger_rate_density_prior = "merger_rate_density_madau_dickinson_belczynski_ng"
            merger_rate_density_prior_params = dict(
                R0=19 * 1e-9, alpha_F=2.57, beta_F=5.83, c_F=3.36  
            )
            source_frame_masses_prior = "binary_masses_BBH_powerlaw_gaussian"
            source_frame_masses_prior_params = dict(
                mminbh=4.98,
                mmaxbh=112.5,
                alpha=3.78,
                mu_g=32.27,
                sigma_g=3.88,
                lambda_peak=0.03,
                delta_m=4.8,
                beta=0.81,
            )
            a_max = 0.8

        elif event_type == "BNS":
            merger_rate_density_prior = "merger_rate_density_madau_dickinson2014"
            merger_rate_density_prior_params = dict(
                R0=89 * 1e-9, a=0.015, b=2.7, c=2.9, d=5.6
            )
            source_frame_masses_prior = "binary_masses_BNS_bimodal"
            source_frame_masses_prior_params = dict(
                w=0.643,
                muL=1.352,
                sigmaL=0.08,
                muR=1.88,
                sigmaR=0.3,
                mmin=1.0,
                mmax=2.3,
            )
            a_max = 0.05

        elif event_type == "NSBH":
            merger_rate_density_prior = "merger_rate_density_madau_dickinson2014"
            merger_rate_density_prior_params = dict(
                R0=23 * 1e-9, a=0.015, b=2.7, c=2.9, d=5.6
            )
            source_frame_masses_prior = "binary_masses_NSBH_broken_powerlaw"
            source_frame_masses_prior_params = dict(
                mminbh=26,
                mmaxbh=125,
                alpha_1=6.75,
                alpha_2=6.75,
                b=0.5,
                delta_m=5,
                mminns=1.0,
                mmaxns=3.0,
                alphans=0.0,
            )
            a_max = 0.8

        elif event_type == "BBH_popIII":
            merger_rate_density_prior = "merger_rate_density_bbh_popIII_ken2022"
            merger_rate_density_prior_params = dict(
                n0=19.2 * 1e-9, aIII=0.66, bIII=0.3, zIII=11.6
            )
            source_frame_masses_prior = "binary_masses_BBH_popIII_lognormal"
            source_frame_masses_prior_params = dict(
                m_min=5.0, m_max=150.0, Mc=30.0, sigma=0.3, chunk_size=chunk_size
            )
            a_max = 0.8

        elif event_type == "BBH_primordial":
            merger_rate_density_prior = "merger_rate_density_bbh_primordial_ken2022"
            merger_rate_density_prior_params = dict(
                n0=0.044 * 1e-9, t0=13.786885302009708
            )
            source_frame_masses_prior = "binary_masses_BBH_primordial_lognormal"
            source_frame_masses_prior_params = dict(
                m_min=1.0, m_max=100.0, Mc=20.0, sigma=0.3, chunk_size=chunk_size
            )
            a_max = 0.8

        else:
            raise ValueError("event_type is not recognized")

        # setting the priors and its parameters
        source_priors_ = dict(
            merger_rate_density=merger_rate_density_prior,
            zs='source_redshift',
            source_frame_masses=source_frame_masses_prior,
            geocent_time="sampler_uniform",
            ra="sampler_uniform",
            dec="sampler_cosine",
            phase="sampler_uniform",
            psi="sampler_uniform",
            theta_jn="sampler_sine",
        )
        source_prior_params_ = dict(
            merger_rate_density=merger_rate_density_prior_params,
            zs=None,
            source_frame_masses=source_frame_masses_prior_params,
            geocent_time=dict(xmin=1238166018, xmax=1269702018),
            ra=dict(xmin=0., xmax=2.*np.pi),
            dec=None,  # dict(xmin=-np.pi/2, xmax=np.pi/2),
            phase=dict(xmin=0., xmax=2.*np.pi),
            psi=dict(xmin=0., xmax=np.pi),
            theta_jn=None,  # dict(xmin=0., xmax=np.pi),
        )

        # spin
        if not self.spin_zero:
            source_priors_["a_1"] = "sampler_uniform"
            source_priors_["a_2"] = "sampler_uniform"

            if self.spin_precession:
                # Precessing spins: magnitude from 0 to a_max
                source_prior_params_["a_1"] = dict(xmin=0.0, xmax=a_max)
                source_prior_params_["a_2"] = dict(xmin=0.0, xmax=a_max)
                source_priors_["tilt_1"] = "sampler_sine"
                source_prior_params_["tilt_1"] = None
                source_priors_["tilt_2"] = "sampler_sine"
                source_prior_params_["tilt_2"] = None
                source_priors_["phi_12"] = "sampler_uniform"
                source_prior_params_["phi_12"] = dict(xmin=0, xmax=2 * np.pi)
                source_priors_["phi_jl"] = "sampler_uniform"
                source_prior_params_["phi_jl"] = dict(xmin=0, xmax=2 * np.pi)
            else:
                # Aligned/anti-aligned spins: magnitude from -a_max to a_max
                source_prior_params_["a_1"] = dict(xmin=-a_max, xmax=a_max)
                source_prior_params_["a_2"] = dict(xmin=-a_max, xmax=a_max)

        # update the priors if input is given
        if source_priors:
            source_priors_.update(source_priors)
        if source_prior_params:
            source_prior_params_.update(source_prior_params)

        # taking care of source_prior_params from the available_gw_prior
        for key, value in source_priors_.items():
            if isinstance(value, str):
                dict_ = self.available_gw_prior[key]  # e.g. all source_frame_masses_prior function names and its parameters
                if value in dict_:
                    param_dict = dict_[value]
                    if source_prior_params_[key] is None:
                        source_prior_params_[key] = param_dict
                    else:
                        param_dict.update(source_prior_params_[key])
                        source_prior_params_[key] = param_dict
                else:
                    raise ValueError(
                        f"source_prior_params_['{key}'] is not in available_gw_prior"
                    )
            elif not callable(value):
                raise ValueError(
                    f"source_prior_params_['{key}'] should be either a string name of available sampler or a function"
                )

        return (source_priors_, source_prior_params_)
    
    def sample_gw_parameters(self, size=1000, param=None):
        """
        Sample all gravitational wave parameters for compact binaries.

        Generates a complete set of intrinsic and extrinsic parameters including
        masses, redshift, luminosity distance, sky position, orientation, and 
        optionally spin parameters.

        Parameters
        ----------
        size : ``int``
            Number of samples to draw. \n
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
            | a_2                |              | spin_2 of the compact binary         |
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

        Examples
        --------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution()
        >>> params = cbc.sample_gw_parameters(size=1000)
        >>> print(list(params.keys()))
        """

        # check for input parameters
        # allow some of the parameters to be fixed 
        if param is None:
            param = {}  # empty
            param_keys = param.keys()  # empty
        else:
            param_keys = param.keys()

        # sample parameters
        param_names = list(self.gw_param_samplers.keys())
        del param_names[0]  # remove merger_rate_density

        gw_parameters = {}  # initialize dictionary to store parameters
        for sampler_name in param_names:
            # print(name)
            if sampler_name not in param_keys:
                # Sample the parameter using the specified sampler function
                # try:
                gw_parameters[sampler_name] = getattr(self, str(sampler_name))(size)
                # except:
                #     raise Exception(f"Sampler {sampler_name} is not defined.")
            else:
                # Use the provided value from kwargs
                gw_parameters[sampler_name] = param[sampler_name]

        # calculate luminosity distance
        zs = gw_parameters["zs"]
        gw_parameters["luminosity_distance"] = self.luminosity_distance(zs)  # Mpc
        
        # mass1 and mass2
        m1, m2 = gw_parameters["source_frame_masses"]  # Msun
        # exchange m1 and m2 if m1 < m2 in the fastest way
        idx = m1 < m2
        m1[idx], m2[idx] = m2[idx], m1[idx]
        gw_parameters["mass_1_source"], gw_parameters["mass_2_source"] = m1, m2
        gw_parameters["mass_1"], gw_parameters["mass_2"] = m1 * (1 + zs), m2 * (
            1 + zs
        )  # Msun
        del gw_parameters["source_frame_masses"]

        return gw_parameters

    def binary_masses_BBH_powerlaw_gaussian(
        self,
        size,
        get_attribute=False,
        **kwargs,
    ):
        """
        Sample source masses with PowerLaw+PEAK model for Population I/II BBH.

        Implements the mass distribution model from LIGO-Virgo population analyses
        combining a power-law with a Gaussian peak component.

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
        >>> m1_src, m2_src = cbc.binary_masses_BBH_powerlaw_gaussian(size=1000)
        """

        from .prior_functions import binary_masses_BBH_powerlaw_gaussian_rvs

        identifier_dict = {'name': "binary_masses_BBH_powerlaw_gaussian"}
        param_dict = self.available_gw_prior["source_frame_masses"]["binary_masses_BBH_powerlaw_gaussian"].copy()
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)


        # mass function
        rvs_ = lambda size: binary_masses_BBH_powerlaw_gaussian_rvs(  
            size=size,
            mminbh=identifier_dict["mminbh"],
            mmaxbh=identifier_dict["mmaxbh"],
            alpha=identifier_dict["alpha"],
            mu_g=identifier_dict["mu_g"],
            sigma_g=identifier_dict['sigma_g'],
            lambda_peak=identifier_dict['lambda_peak'],
            delta_m=identifier_dict['delta_m'],
            beta=identifier_dict['beta'],
        )

        mass_object = FunctionConditioning(
            function=None,
            x_array=None,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="source_frame_masses",
            name=identifier_dict['name'],
            create_new=self.create_new_interpolator["source_frame_masses"]["create_new"],
            create_function_inverse=False,
            create_function=False,
            create_pdf=False,
            create_rvs=rvs_,
            callback='rvs',
        )

        return mass_object if get_attribute else mass_object.rvs(size)

    def binary_masses_BBH_popIII_lognormal(
        self,
        size,
        get_attribute=False,
        **kwargs,
    ):
        """
        Sample source masses for Population III BBH from lognormal distribution.

        Based on Eqn. 1 and 4 of Ng et al. 2022 for Population III black holes.

        Parameters
        ----------
        size : ``int``
            Number of samples to draw.
        get_attribute : ``bool``
            If True, return the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Model parameters: \n
            - m_min: Minimum BH mass (Msun), default: 5.0 \n
            - m_max: Maximum BH mass (Msun), default: 150.0 \n
            - Mc: Central mass scale (Msun), default: 30.0 \n
            - sigma: Distribution width, default: 0.3

        Returns
        -------
        mass_1_source : ``numpy.ndarray``
            Array of primary masses in source frame (Msun).
        mass_2_source : ``numpy.ndarray``
            Array of secondary masses in source frame (Msun).

        Examples
        --------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution(event_type='BBH_popIII')
        >>> m1_src, m2_src = cbc.binary_masses_BBH_popIII_lognormal(size=1000)
        """

        from .prior_functions import binary_masses_BBH_popIII_lognormal_rvs

        identifier_dict = {'name': "binary_masses_BBH_popIII_lognormal"}
        param_dict = self.available_gw_prior["source_frame_masses"]["binary_masses_BBH_popIII_lognormal"].copy()
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        rvs_ = lambda size: binary_masses_BBH_popIII_lognormal_rvs(  
            size,
            m_min=identifier_dict["m_min"],
            m_max=identifier_dict["m_max"],
            Mc=identifier_dict["Mc"],
            sigma=identifier_dict["sigma"],
            chunk_size=chunk_size,
        )

        mass_object = FunctionConditioning(
            function=None,
            x_array=None,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="source_frame_masses",
            name=identifier_dict['name'],
            create_new=self.create_new_interpolator["source_frame_masses"]["create_new"],
            create_function_inverse=False,
            create_function=False,
            create_pdf=False,
            create_rvs=rvs_,
            callback='rvs',
        )

        return mass_object if get_attribute else mass_object.rvs(size)

    def binary_masses_BBH_primordial_lognormal(
        self,
        size,
        get_attribute=False,
        **kwargs,
    ):
        """
        Sample source masses for primordial BBH from lognormal distribution.

        Based on Eqn. 1 and 4 of Ng et al. 2022 for primordial black holes.

        Parameters
        ----------
        size : ``int``
            Number of samples to draw.
        get_attribute : ``bool``
            If True, return the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Model parameters: \n
            - m_min: Minimum BH mass (Msun), default: 1.0 \n
            - m_max: Maximum BH mass (Msun), default: 100.0 \n
            - Mc: Central mass scale (Msun), default: 20.0 \n
            - sigma: Distribution width, default: 0.3

        Returns
        -------
        mass_1_source : ``numpy.ndarray``
            Array of primary masses in source frame (Msun).
        mass_2_source : ``numpy.ndarray``
            Array of secondary masses in source frame (Msun).
        """

        from .prior_functions import binary_masses_BBH_primordial_lognormal_rvs

        identifier_dict = {'name': "binary_masses_BBH_primordial_lognormal"}    
        param_dict = self.available_gw_prior["source_frame_masses"]["binary_masses_BBH_primordial_lognormal"].copy()
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        rvs_ = lambda size: binary_masses_BBH_primordial_lognormal_rvs(  
            size,
            m_min=identifier_dict["m_min"],
            m_max=identifier_dict["m_max"],
            Mc=identifier_dict["Mc"],
            sigma=identifier_dict["sigma"],
            chunk_size=chunk_size,
        )

        mass_object = FunctionConditioning(
            function=None,
            x_array=None,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="source_frame_masses",
            name=identifier_dict['name'],
            create_new=self.create_new_interpolator["source_frame_masses"]["create_new"],
            create_function_inverse=False,
            create_function=False,
            create_pdf=False,
            create_rvs=rvs_,
            callback='rvs',
        )

        return mass_object if get_attribute else mass_object.rvs(size)
        
    def binary_masses_NSBH_broken_powerlaw(self, size, get_attribute=False, **kwargs):
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
            - mminns: Minimum NS mass (Msun), default: 1.0 \n
            - mmaxns: Maximum NS mass (Msun), default: 3.0 \n
            - alphans: NS mass power-law index, default: 0.0

        Returns
        -------
        mass_1_source : ``numpy.ndarray``
            Array of BH masses in source frame (Msun).
        mass_2_source : ``numpy.ndarray``
            Array of NS masses in source frame (Msun).

        Examples
        --------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution(event_type='NSBH')
        >>> m1_src, m2_src = cbc.binary_masses_NSBH_broken_powerlaw(size=1000)
        """

        from .prior_functions import binary_masses_NSBH_broken_powerlaw_rvs

        identifier_dict = {'name': "binary_masses_NSBH_broken_powerlaw"}
        param_dict = self.available_gw_prior["source_frame_masses"]["binary_masses_NSBH_broken_powerlaw"].copy()
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        # mass function for NSBH
        rvs_ = lambda size: binary_masses_NSBH_broken_powerlaw_rvs(
            size,
            mminbh=identifier_dict["mminbh"],
            mmaxbh=identifier_dict["mmaxbh"],
            alpha_1=identifier_dict["alpha_1"],
            alpha_2=identifier_dict["alpha_2"],
            b=identifier_dict["b"],
            delta_m=identifier_dict["delta_m"],
            mminns=identifier_dict["mminns"],
            mmaxns=identifier_dict["mmaxns"],
            alphans=identifier_dict["alphans"],
        )

        mass_object = FunctionConditioning(
            function=None,
            x_array=None,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="source_frame_masses",
            name=identifier_dict['name'],
            create_new=self.create_new_interpolator["source_frame_masses"]["create_new"],
            create_function_inverse=False,
            create_function=False,
            create_pdf=False,
            create_rvs=rvs_,
            callback='rvs',
        )

        return mass_object if get_attribute else mass_object.rvs(size)

    def binary_masses_uniform(
        self,
        size,
        get_attribute=False,
        **kwargs,
    ):
        """
        Sample source masses from uniform distribution.

        Parameters
        ----------
        size : ``int``
            Number of samples to draw.
        get_attribute : ``bool``
            If True, return the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Model parameters: \n
            - m_min: Minimum mass (Msun), default: 1.0 \n
            - m_max: Maximum mass (Msun), default: 3.0

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
        >>> m1_src, m2_src = cbc.binary_masses_uniform(size=1000)
        """

        identifier_dict = {'name': "binary_masses_uniform"}
        param_dict = self.available_gw_prior["source_frame_masses"]["binary_masses_uniform"].copy()
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        m_min = identifier_dict["m_min"]
        m_max = identifier_dict["m_max"]

        @njit
        def rvs_(size):
            # generate random masses
            mass_1_source = np.random.uniform(m_min, m_max, size)
            mass_2_source = np.random.uniform(m_min, m_max, size)
            # swap if mass_2_source > mass_1_source
            idx = np.where(mass_2_source > mass_1_source)
            mass_1_source[idx], mass_2_source[idx] = mass_2_source[idx], mass_1_source[idx] 
            return mass_1_source, mass_2_source

        mass_object = FunctionConditioning(
            function=None,
            x_array=None,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="source_frame_masses",
            name=identifier_dict['name'],
            create_new=self.create_new_interpolator["source_frame_masses"]["create_new"],
            create_function_inverse=False,
            create_function=False,
            create_pdf=False,
            create_rvs=rvs_,
            callback='rvs',
        )

        return mass_object if get_attribute else mass_object.rvs(size)

    def binary_masses_BNS_bimodal(
        self,
        size,
        get_attribute=False,
        **kwargs,
    ):
        """
        Sample BNS masses from bimodal Gaussian distribution.

        Based on Will M. Farr et al. 2020 Eqn. 6 for neutron star mass
        distribution combining two Gaussian peaks.

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
        mass_2_source : ``numpy.ndarray``
            Array of secondary masses in source frame (Msun).

        Examples
        --------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution(event_type='BNS')
        >>> m1_src, m2_src = cbc.binary_masses_BNS_bimodal(size=1000)
        """

        from .prior_functions import binary_masses_BNS_bimodal_rvs

        identifier_dict = {'name': "binary_masses_BNS_bimodal"}
        identifier_dict['resolution'] = self.create_new_interpolator["source_frame_masses"]["resolution"]
        param_dict = self.available_gw_prior["source_frame_masses"]["binary_masses_BNS_bimodal"].copy()
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        # mass function for BNS 
        rvs_ = lambda size: binary_masses_BNS_bimodal_rvs(  
            size,
            w=identifier_dict["w"],
            muL=identifier_dict["muL"],
            sigmaL=identifier_dict["sigmaL"],
            muR=identifier_dict["muR"],
            sigmaR=identifier_dict["sigmaR"],
            mmin=identifier_dict["mmin"],
            mmax=identifier_dict["mmax"],
            resolution=identifier_dict["resolution"],
        )
        
        mass_object = FunctionConditioning(
            function=None,
            x_array=None,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="source_frame_masses",
            name=identifier_dict['name'],
            create_new=self.create_new_interpolator["source_frame_masses"]["create_new"],
            create_function_inverse=False,
            create_function=False,
            create_pdf=False,
            create_rvs=rvs_,
            callback='rvs',
        )

        return mass_object if get_attribute else mass_object.rvs(size)

    def constant_values_n_size(
        self, size=100, get_attribute=False, **kwargs
    ):
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

        # if param:
        #     value = param["value"]
        identifier_dict = {'name': "constant_values_n_size"}
        param_dict = dict(value=0.0)
        param_dict.update(kwargs)
        identifier_dict.update(param_dict)

        value = identifier_dict["value"]
        # pdf_, zero everywhere except at value
        pdf_ = njit(lambda x: np.where(x == value, 1.0, 0.0))
        # rvs_, return value
        rvs_ = njit(lambda size: np.ones(size) * value)

        object_ = FunctionConditioning(
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="custom_functions",
            name=identifier_dict['name'],
            create_function_inverse=False,
            create_function=False,
            create_pdf=pdf_,
            create_rvs=rvs_,
            callback='rvs',
        )

        if get_attribute:
            return object_
        else:
            return object_(size)

    def sampler_uniform(
        self, size, get_attribute=False, **kwargs
    ):
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
            - xmin: Minimum value, default: 0.0 \n
            - xmax: Maximum value, default: 1.0

        Returns
        -------
        values : ``numpy.ndarray``
            Array of uniformly distributed values in range [xmin, xmax].
        """
        
        identifier_dict = {'name': "sampler_uniform"}
        param_dict = dict(xmin=0.0, xmax=1.0)
        param_dict.update(kwargs)
        identifier_dict.update(param_dict)

        xmin = identifier_dict['xmin']
        xmax = identifier_dict['xmax']

        pdf_ = njit(lambda x: 1.0 / (xmax - xmin) * np.ones(len(x)))
        rvs_ = njit(lambda size: np.random.uniform(xmin, xmax, size=size))

        object_ = FunctionConditioning(
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="custom_functions",
            name=identifier_dict['name'],
            create_function_inverse=False,
            create_function=False,
            create_pdf=pdf_,
            create_rvs=rvs_,
            callback='rvs',
        )

        if get_attribute:
            return object_
        else:
            return object_(size)

    def sampler_cosine(
        self, size, get_attribute=False, **kwargs
    ):
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

        identifier_dict = {}
        identifier_dict['name'] = "sampler_cosine"
        identifier_dict.update(kwargs)

        pdf_ = njit(lambda x: 0.5 * np.cos(x))
        rvs_ = njit(lambda size: np.arcsin((np.random.uniform(0, 1, size=size) * 2 - 1)))

        object_ = FunctionConditioning(
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="custom_functions",
            name=identifier_dict['name'],
            create_function_inverse=False,
            create_function=False,
            create_pdf=pdf_,
            create_rvs=rvs_,
            callback='rvs',
        )

        if get_attribute:
            return object_
        else:
            return object_(size)

    def sampler_sine(
        self, size, get_attribute=False, **kwargs
    ):
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

        identifier_dict = {}
        identifier_dict['name'] = "sampler_sine"
        identifier_dict.update(kwargs)

        pdf_ = njit(lambda x: 0.5 * np.sin(x))
        rvs_ = njit(lambda size: np.arccos((np.random.uniform(0, 1, size=size) - 0.5) * 2))

        object_ = FunctionConditioning(
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="custom_functions",
            name=identifier_dict['name'],
            create_function_inverse=False,
            create_function=False,
            create_pdf=pdf_,
            create_rvs=rvs_,
            callback='rvs',
        )   

        if get_attribute:
            return object_
        else:
            return object_(size)

    @property
    def source_frame_masses(self):
        """
        Class object (of FunctionConditioning) for source frame masses, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the density profile slope distribution

        Returns
        -------
        mass_1_source : ``numpy.ndarray``
            Array of mass_1_source values in solar masses.
        mass_2_source : ``numpy.ndarray``
            Array of mass_2_source values in solar masses.

        Examples
        --------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc_source_param_dist = CBCSourceParameterDistribution()
        >>> cbc_source_param_dist.source_frame_masses(size=10)
        """

        return self._source_frame_masses

    @source_frame_masses.setter
    def source_frame_masses(self, prior):
        if prior in self.available_gw_prior["source_frame_masses"]:
            print(f"using ler available source_frame_masses function : {prior}")
            args = self.gw_param_samplers_params["source_frame_masses"]
            if args is None:
                self._source_frame_masses = getattr(self, prior)(
                size=None, get_attribute=True
            )
            else:
                # follwing should return a sampler function with only one argument (size)
                self._source_frame_masses = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print("using user defined custom source_frame_masses class/object of type ler.utils.FunctionConditioning")
            self._source_frame_masses = prior
        elif callable(prior):
            print("using user defined custom source_frame_masses function")
            self._source_frame_masses = prior
        else:
            raise ValueError(
                "source_frame_masses prior not available in available_gw_prior. Must be a string or a callable function."
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
            print("using user defined custom zs class/object of type ler.utils.FunctionConditioning")
            self._zs = prior
        elif callable(prior):
            print("using user defined custom zs function")
            self._zs = prior
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
                self._geocent_time = getattr(self, prior)(
                size=None, get_attribute=True
            )
            else:
                # follwing should return a sampler function with only one argument (size)
                self._geocent_time = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print("using user defined custom geocent_time class/object of type ler.utils.FunctionConditioning")
            self._geocent_time = prior
        elif callable(prior):
            print("using user defined custom geocent_time function")
            self._geocent_time = prior
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
                # follwing should return a sampler function with only one argument (size)
                self._ra = getattr(self, prior)(size=None, get_attribute=True, **args)
        elif isinstance(prior, FunctionConditioning):
            print("using user defined custom ra class/object of type ler.utils.FunctionConditioning")
            self._ra = prior
        elif callable(prior):
            print("using user defined custom ra function")
            self._ra = prior
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
                # follwing should return a sampler function with only one argument (size)
                self._dec = getattr(self, prior)(size=None, get_attribute=True, **args)
        elif isinstance(prior, FunctionConditioning):
            print("using user defined custom dec class/object of type ler.utils.FunctionConditioning")
            self._dec = prior
        elif callable(prior):
            print("using user defined custom dec function")
            self._dec = prior
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
                # follwing should return a sampler function with only one argument (size)
                self._phase = getattr(self, prior)(size=None, get_attribute=True, **args)
        elif isinstance(prior, FunctionConditioning):
            print("using user defined custom phase class/object of type ler.utils.FunctionConditioning")
            self._phase = prior
        elif callable(prior):
            print("using user defined custom phase function")
            self._phase = prior
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
        geocent_time : ``numpy.ndarray``
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
                # follwing should return a sampler function with only one argument (size)
                self._psi = getattr(self, prior)(size=None, get_attribute=True, **args)
        elif isinstance(prior, FunctionConditioning):
            print("using user defined custom psi class/object of type ler.utils.FunctionConditioning")
            self._psi = prior
        elif callable(prior):
            print("using user defined custom psi function")
            self._psi = prior
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
                self._theta_jn = getattr(self, prior)(
                size=None, get_attribute=True
            )
            else:
                # follwing should return a sampler function with only one argument (size)
                self._theta_jn = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print("using user defined custom theta_jn class/object of type ler.utils.FunctionConditioning")
            self._theta_jn = prior
        elif callable(prior):
            print("using user defined custom theta_jn function")
            self._theta_jn = prior
        else:
            raise ValueError("theta_jn prior not available in available_gw_prior. Must be a string or a callable function.")

    @property
    def a_1(self):
        """
        Class object (of FunctionConditioning) for spin1 magnitude, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the spin1 magnitude distribution \n
        - `pdf`: returns the probability density function of the spin1 magnitude distribution \n
        - `function`: returns the spin1 magnitude distribution function.

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
                self._a_1 = getattr(self, prior)(
                size=None, get_attribute=True
            )
            else:
                # follwing should return a sampler function with only one argument (size)
                self._a_1 = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print("using user defined custom a_1 class/object of type ler.utils.FunctionConditioning")
            self._a_1 = prior
        elif callable(prior):
            print("using user defined custom a_1 function")
            self._a_1 = prior
        else:
            raise ValueError("a_1 prior not available in available_gw_prior. Must be a string or a callable function.")

    @property
    def a_2(self):
        """
        Class object (of FunctionConditioning) for spin2 magnitude, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the spin2 magnitude distribution \n
        - `pdf`: returns the probability density function of the spin2 magnitude distribution \n
        - `function`: returns the spin2 magnitude distribution function.

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
                self._a_2 = getattr(self, prior)(
                size=None, get_attribute=True
            )
            else:
                # follwing should return a sampler function with only one argument (size)
                self._a_2 = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print("using user defined custom a_2 class/object of type ler.utils.FunctionConditioning")
            self._a_2 = prior
        elif callable(prior):
            print("using user defined custom a_2 function")
            self._a_2 = prior
        else:
            raise ValueError("a_2 prior not available in available_gw_prior. Must be a string or a callable function.")

    @property
    def tilt_1(self):
        """
        Class object (of FunctionConditioning) for tilt1 angle, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the tilt1 angle distribution \n
        - `pdf`: returns the probability density function of the tilt1 angle distribution \n
        - `function`: returns the tilt1 angle distribution function.

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
                self._tilt_1 = getattr(self, prior)(
                size=None, get_attribute=True
            )
            else:
                self._tilt_1 = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print("using user defined custom tilt_1 class/object of type ler.utils.FunctionConditioning")
            self._tilt_1 = prior
        elif callable(prior):
            print("using user defined custom tilt_1 function")
            self._tilt_1 = prior
        else:
            raise ValueError("tilt_1 prior not available in available_gw_prior. Must be a string or a callable function.")

    @property
    def tilt_2(self):
        """
        Class object (of FunctionConditioning) for tilt2 angle, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the tilt2 angle distribution \n
        - `pdf`: returns the probability density function of the tilt2 angle distribution \n
        - `function`: returns the tilt2 angle distribution function.

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
                self._tilt_2 = getattr(self, prior)(
                size=None, get_attribute=True
            )
            else:
                self._tilt_2 = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print("using user defined custom tilt_2 class/object of type ler.utils.FunctionConditioning")
            self._tilt_2 = prior
        elif callable(prior):
            print("using user defined custom tilt_2 function")
            self._tilt_2 = prior
        else:
            raise ValueError("tilt_2 prior not available in available_gw_prior. Must be a string or a callable function.")
        
    @property
    def phi_12(self):
        """
        Class object (of FunctionConditioning) for phi_12 angle, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the phi_12 angle distribution \n
        - `pdf`: returns the probability density function of the phi_12 angle distribution \n
        - `function`: returns the phi_12 angle distribution function.

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
                self._phi_12 = getattr(self, prior)(
                    size=None, get_attribute=True
                )
            else:
                self._phi_12 = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print("using user defined custom phi_12 class/object of type ler.utils.FunctionConditioning")
            self._phi_12 = prior
        elif callable(prior):
            print("using user defined custom phi_12 function")
            self._phi_12 = prior
        else:
            raise ValueError("phi_12 prior not available in available_gw_prior. Must be a string or a callable function.")

    @property
    def phi_jl(self):
        """
        Class object (of FunctionConditioning) for phi_jl angle, with rvs/sampler as callback. Can also be a user defined callable sampler. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the phi_jl angle distribution \n
        - `pdf`: returns the probability density function of the phi_jl angle distribution \n
        - `function`: returns the phi_jl angle distribution function.

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
                self._phi_jl = getattr(self, prior)(
                size=None, get_attribute=True
            )
            else:
                self._phi_jl = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print("using user defined custom phi_jl class/object of type ler.utils.FunctionConditioning")
            self._phi_jl = prior
        elif callable(prior):
            print("using user defined custom phi_jl function")
            self._phi_jl = prior
        else:
            raise ValueError("phi_jl prior not available in available_gw_prior. Must be a string or a callable function.")

    @property
    def available_gw_prior(self):
        """
        Dictionary of all available prior distributions and their parameters. \n
        This is a dynamically generated dictionary containing available samplers
        for each GW parameter type and their default parameter values.

        Returns
        -------
        available_gw_prior : ``dict``
            Nested dictionary organized by parameter type (e.g., 'source_frame_masses', \n
            'geocent_time', etc.) with sampler names and default parameters.
        """

        self._available_gw_prior = dict(
            merger_rate_density=self.available_merger_rate_density_model,
            zs=dict(
                source_redshift=None,
            ),
            source_frame_masses=dict(
                binary_masses_BBH_powerlaw_gaussian=dict(
                    mminbh=4.98,
                    mmaxbh=112.5,
                    alpha=3.78,
                    mu_g=32.27,
                    sigma_g=3.88,
                    lambda_peak=0.03,
                    delta_m=4.8,
                    beta=0.81,
                ),
                binary_masses_BBH_popIII_lognormal=dict(m_min=5.0, m_max=150.0, Mc=30.0, sigma=0.3),
                binary_masses_BBH_primordial_lognormal=dict(
                    m_min=1.0, m_max=100.0, Mc=20.0, sigma=0.3
                ),
                binary_masses_NSBH_broken_powerlaw=dict(
                    mminbh=26,
                    mmaxbh=125,
                    alpha_1=6.75,
                    alpha_2=6.75,
                    b=0.5,
                    delta_m=5,
                    mminns=1.0,
                    mmaxns=3.0,
                    alphans=0.0,
                ),
                binary_masses_uniform=dict(m_min=1.0, m_max=3.0),
                binary_masses_BNS_bimodal=dict(
                    w=0.643,
                    muL=1.352,
                    sigmaL=0.08,
                    muR=1.88,
                    sigmaR=0.3,
                    mmin=1.0,
                    mmax=2.3,
                ),
            ),
            a_1=dict(
                constant_values_n_size=dict(value=0.0),
                sampler_uniform=dict(xmin=-0.8, xmax=0.8),
            ) if not self.spin_precession else dict(
                constant_values_n_size=dict(value=0.0),
                sampler_uniform=dict(xmin=0.0, xmax=0.8),
            ),
            a_2=dict(
                constant_values_n_size=dict(value=0.0),
                sampler_uniform=dict(xmin=-0.8, xmax=0.8),
            ) if not self.spin_precession else dict(
                constant_values_n_size=dict(value=0.0),
                sampler_uniform=dict(xmin=0.0, xmax=0.8),
            ),
            tilt_1=dict(
                constant_values_n_size=dict(value=0.0),
                sampler_sine=None,
            ),
            tilt_2=dict(
                constant_values_n_size=dict(value=0.0),
                sampler_sine=None,
            ),
            phi_12=dict(
                constant_values_n_size=dict(value=0.0),
                sampler_uniform=dict(xmin=0.0, xmax=2 * np.pi),
            ),
            phi_jl=dict(
                constant_values_n_size=dict(value=0.0),
                sampler_uniform=dict(xmin=0.0, xmax=2 * np.pi),
            ),
            geocent_time=dict(
                sampler_uniform=dict(
                    xmin=1238166018, xmax=1238166018 + 31557600.0
                ),
                constant_values_n_size=dict(value=1238166018),
            ),
            ra=dict(
                sampler_uniform=dict(xmin=0.0, xmax=2 * np.pi),
                constant_values_n_size=dict(value=0.0),
            ),
            dec=dict(
                sampler_cosine=None,
                constant_values_n_size=dict(value=0.0),
                sampler_uniform=dict(xmin=-np.pi / 2, xmax=np.pi / 2),
            ),
            phase=dict(
                sampler_uniform=dict(xmin=0.0, xmax=2 * np.pi),
                constant_values_n_size=dict(value=0.0),
            ),
            psi=dict(
                sampler_uniform=dict(xmin=0.0, xmax=np.pi),
                constant_values_n_size=dict(value=0.0),
            ),
            theta_jn=dict(
                sampler_sine=None,
                constant_values_n_size=dict(value=0.0),
                sampler_uniform=dict(xmin=0.0, xmax=np.pi),
            ),
        )

        return self._available_gw_prior