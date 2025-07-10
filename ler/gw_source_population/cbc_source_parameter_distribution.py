# -*- coding: utf-8 -*-
"""
This module contains the classes to generate source compact binary's intrinsic and extrinsic gravitational waves parameters.
"""

import warnings

warnings.filterwarnings("ignore")
import numpy as np
from numba import njit

# for redshift to luminosity distance conversion
from astropy.cosmology import LambdaCDM

# from gwcosmo import priors as p
from scipy.integrate import quad

# for multiprocessing
# Import helper routines
from ..utils import FunctionConditioning

# import redshift distribution sampler
from .cbc_source_redshift_distribution import CBCSourceRedshiftDistribution

from .jit_functions import lognormal_distribution_2D, bns_bimodal_pdf, inverse_transform_sampler_m1m2, sample_powerlaw_gaussian_source_bbh_masses, sample_broken_powerlaw_nsbh_masses

chunk_size = 10000


class CBCSourceParameterDistribution(CBCSourceRedshiftDistribution):
    """Class to generate a population of compact binaries. It helps sample all the intrinsic and extrinsic parameters of compact binaries. This daughter class inherits from :class:`~ler.ler.CBCSourceRedshiftDistribution` class.

    Parameters
    ----------
    z_min : `float`
        Minimum redshift of the source population
        default: 0.001
    z_max : `float`
        Maximum redshift of the source population
        default: 10.
    event_type : `str`
        Type of event to generate.
        e.g. 'BBH', 'BNS', 'NSBH'
    source_priors, source_priors_params : `dict`, `dict`
        Dictionary of prior sampler functions and its input parameters.
        Check for available priors and corresponding input parameters by running,
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution()
        >>> cbc.available_gw_prior_list_and_its_params()
        # To check the current chosen priors and its parameters, run,
        >>> print("default priors=",cbc.gw_param_samplers)
        >>> print("default priors's parameters=",cbc.gw_param_samplers_params)
    cosmology : `astropy.cosmology`
        Cosmology to use
        default: None/astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)
    spin_zero : `bool`
        If True, spin parameters are completely ignore in the sampling.
        default: True
    spin_precession : `bool`
        If spin_zero=False and spin_precession=True, spin parameters are sampled for precessing binaries.
        if spin_zero=False and spin_precession=False, spin parameters are sampled for aligned/anti-aligned spin binaries.
        default: False
    directory : `str`
        Directory to store the interpolator pickle files
        default: './interpolator_pickle'
    create_new_interpolator : `dict`
        Dictionary of boolean values and resolution to create new interpolator.
        default: dict(redshift_distribution=dict(create_new=False, resolution=500), z_to_luminosity_distance=dict(create_new=False, resolution=500), differential_comoving_volume=dict(create_new=False, resolution=500))

    Examples
    ----------
    >>> from ler.gw_source_population import CBCSourceParameterDistribution
    >>> cbc = CBCSourceParameterDistribution()
    >>> params = cbc.gw_parameters(size=1000)
    >>> print("sampled parameters=",list(params.keys()))

    Instance Attributes
    ----------
    CBCSourceParameterDistribution has the following instance attributes:\n
    +-------------------------------------+----------------------------------+
    | Atrributes                          | Type                             |
    +=====================================+==================================+
    |:attr:`~z_min`                       | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`~z_max`                       | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`~event_type`                  | `str`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~source_priors`               | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~source_priors_params`        | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~cosmo`                       | `astropy.cosmology`              |
    +-------------------------------------+----------------------------------+
    |:attr:`~spin_zero`                   | `bool`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~spin_precession`             | `bool`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~directory`                   | `str`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~create_new_interpolator`     | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~available_gw_prior_list_and_its_params`                            |
    +-------------------------------------+----------------------------------+
    |                                     | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~gw_param_samplers`           | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~gw_param_samplers_params`    | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~sampler_names`               | `dict`                           |
    +-------------------------------------+----------------------------------+

    Instance Methods
    ----------
    CBCSourceParameterDistribution has the following instance methods:\n
    +-------------------------------------+----------------------------------+
    | Methods                             | Type                             |
    +=====================================+==================================+
    |:meth:`~source_priors_categorization`                                   |
    +-------------------------------------+----------------------------------+
    |                                     | Function to categorize the event |
    |                                     | priors and its parameters        |
    +-------------------------------------+----------------------------------+
    |:meth:`~lookup_table_luminosity_distance`                               |
    |                                     | Function to create a lookup      |
    |                                     | table for converting redshift    |
    |                                     | to luminosity distance           |
    +-------------------------------------+----------------------------------+
    |:meth:`~gw_parameters`        | Function to sample all the       |
    |                                     | intrinsic and extrinsic          |
    |                                     | parameters of compact binaries   |
    +-------------------------------------+----------------------------------+
    |:meth:`~source_frame_masses`  | Function to sample source mass1  |
    |                                     | and mass2                        |
    +-------------------------------------+----------------------------------+
    |:meth:`~geocent_time`         | Function to sample geocent time  |
    +-------------------------------------+----------------------------------+
    |:meth:`~zs`                   | Function to sample source        |
    |                                     | redshift                         |
    +-------------------------------------+----------------------------------+
    |:meth:`~ra`                   | Function to sample right         |
    |                                     | ascension (sky position)         |
    +-------------------------------------+----------------------------------+
    |:meth:`~dec`                  | Function to sample declination   |
    |                                     | (sky position)                   |
    +-------------------------------------+----------------------------------+
    |:meth:`~phase`                | Function to sample coalescence   |
    |                                     | phase                            |
    +-------------------------------------+----------------------------------+
    |:meth:`~psi`                  | Function to sample polarization  |
    |                                     | angle                            |
    +-------------------------------------+----------------------------------+
    |:meth:`~theta_jn`             | Function to sample inclination   |
    |                                     | angle                            |
    +-------------------------------------+----------------------------------+
    |:meth:`~a_1`                   | Function to sample spin1         |
    |                                     | magnitude                        |
    +-------------------------------------+----------------------------------+
    |:meth:`~a_2`                   | Function to sample spin2         |
    |                                     | magnitude                        |
    +-------------------------------------+----------------------------------+
    |:meth:`~tilt_1`               | Function to sample tilt1 angle   |
    +-------------------------------------+----------------------------------+
    |:meth:`~tilt_2`               | Function to sample tilt2 angle   |
    +-------------------------------------+----------------------------------+
    |:meth:`~phi_12`               | Function to sample phi12 angle   |
    +-------------------------------------+----------------------------------+
    |:meth:`~phi_jl`               | Function to sample phi_jl angle  |
    +-------------------------------------+----------------------------------+
    |:meth:`~binary_masses_BBH_popI_II_powerlaw_gaussian`                    |
    +-------------------------------------+----------------------------------+
    |                                     | Function to sample source mass1  |
    |                                     | and mass2 with PowerLaw+PEAK     |
    |                                     | model                            |
    +-------------------------------------+----------------------------------+
    |:meth:`~binary_masses_BBH_popIII_lognormal`                             |
    +-------------------------------------+----------------------------------+
    |                                     | Function to sample source mass1  |
    |                                     | and mass2 with popIII orgin from |
    |                                     | lognormal distribution. Refer to |
    |                                     | Ng et al. 2022. Eqn. 1 and 4     |
    +-------------------------------------+----------------------------------+
    |:meth:`~binary_masses_BBH_primordial_lognormal`                         |
    +-------------------------------------+----------------------------------+
    |                                     | Function to sample source mass1  |
    |                                     | and mass2 with primordial orgin  |
    |                                     | from lognormal distribution.     |
    |                                     | Refer to Ng et al. 2022. Eqn. 1  |
    |                                     | and 4                            |
    +-------------------------------------+----------------------------------+
    |:meth:`~binary_masses_BNS_bimodal`   | Function to sample source mass1  |
    |                                     | and mass2 from bimodal           |
    |                                     | distribution. Refer to           |
    |                                     | Will M. Farr et al. 2020 Eqn. 6  |
    +-------------------------------------+----------------------------------+
    |:meth:`~constant_values_n_size`      | Function to return array of      |
    |                                     | constant values of size n        |
    +-------------------------------------+----------------------------------+
    |:meth:`~sampler_uniform`             | Function to sample from uniform  |
    |                                     | distribution                     |
    +-------------------------------------+----------------------------------+
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

    cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
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
        spin_zero=True,
        spin_precession=False,
        directory="./interpolator_pickle",
        create_new_interpolator=False,
    ):
        # set attributes
        self.z_min = z_min
        self.z_max = z_max
        self.cosmo = cosmology if cosmology else LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        # note that self.cosmo is initialized in the super class
        self.spin_zero = spin_zero
        self.spin_precession = spin_precession
        self.directory = directory

        # setting up the interpolator creation parameters
        create_new_interpolator = self.setup_decision_dictionary_gw_params(create_new_interpolator)

        # dealing with prior functions and categorization
        (
            self.gw_param_samplers,
            self.gw_param_samplers_params,
        ) = self.source_priors_categorization(
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

        print("\nInitializing CBCSourceParameterDistribution...\n")
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


    def setup_decision_dictionary_gw_params(self, create_new_interpolator):
        """
        Method to set up a decision dictionary for interpolator creation.

        Parameters
        ----------
        create_new_interpolator : `dict`, `bool`
            If `dict`, dictionary of boolean values and resolution to create new interpolator.
            If `bool`, boolean value to create new interpolator for all quantities.

        Returns
        -------
        create_new_interpolator_ : `dict`
            Dictionary of boolean values and resolution to create new interpolator.
            e.g. dict(redshift_distribution=dict(create_new=False, resolution=1000), luminosity_distance=dict(create_new=False, resolution=1000), differential_comoving_volume=dict(create_new=False, resolution=1000))
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

    def source_priors_categorization(
        self, event_type, source_priors, source_prior_params
    ):
        """
        Function to categorize the event priors and its parameters.

        Parameters
        ----------
        event_type : `str`
            Type of event to generate.
            e.g. 'BBH', 'BNS', 'BBH_popIII', 'BBH_primordial', 'NSBH'
        source_priors : `dict`
            Dictionary of prior sampler functions for each parameter
        source_prior_params : `dict`
            Dictionary of sampler parameters for each GW parameter

        Returns
        ----------
        source_priors_ : `dict`
            Dictionary of prior sampler functions for each parameter
        source_prior_params_ : `dict`
            Dictionary of sampler parameters for each parameter
        sampler_names_ : `dict`
            Dictionary of sampler names with description

        Examples
        ----------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution()
        >>> source_priors, source_prior_params, sampler_names = cbc.source_priors_categorization(event_type='BBH', source_priors=None, source_prior_params=None)
        >>> print(source_priors.keys())
        >>> print(source_prior_params.keys())
        >>> print(sampler_names.keys())
        """

        # for BBH
        if event_type == "BBH":
            merger_rate_density_prior = "merger_rate_density_bbh_popI_II_oguri2018"
            merger_rate_density_prior_params = dict(
                R0=23.9 * 1e-9, b2=1.6, b3=2.1, b4=30  # 
            )
            source_frame_masses_prior = "binary_masses_BBH_popI_II_powerlaw_gaussian"
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
            merger_rate_density_prior = "merger_rate_density_bbh_popI_II_oguri2018"
            merger_rate_density_prior_params = dict(
                R0=105.5 * 1e-9, b2=1.6, b3=2.1, b4=30
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
            merger_rate_density_prior = "merger_rate_density_bbh_popI_II_oguri2018"
            merger_rate_density_prior_params = dict(
                R0=27.0 * 1e-9, b2=1.6, b3=2.1, b4=30
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
            geocent_time=dict(min_=1238166018, max_=1269702018),
            ra=dict(min_=0., max_=2.*np.pi),
            dec=None,  # dict(min_=-np.pi/2, max_=np.pi/2),
            phase=dict(min_=0., max_=2.*np.pi),
            psi=dict(min_=0., max_=np.pi),
            theta_jn=None,  # dict(min_=0., max_=np.pi),
        )

        # spin
        if not self.spin_zero:
            source_priors_["a_1"] = "sampler_uniform"
            source_prior_params_["a_1"] = dict(min_=-a_max, max_=a_max)
            source_priors_["a_2"] = "sampler_uniform"
            source_prior_params_["a_2"] = dict(min_=-a_max, max_=a_max)

            if self.spin_precession:
                source_priors_["a_1"] = "sampler_uniform"
                source_prior_params_["a_1"] = dict(min_=0.0, max_=a_max)
                source_priors_["a_2"] = "sampler_uniform"
                source_prior_params_["a_2"] = dict(min_=0.0, max_=a_max)
                source_priors_["tilt_1"] = "sampler_sine"
                source_prior_params_["tilt_1"] = None

                source_priors_["tilt_2"] = "sampler_sine"
                source_prior_params_["tilt_2"] = None

                source_priors_["phi_12"] = "sampler_uniform"
                source_prior_params_["phi_12"] = dict(min_=0, max_=2 * np.pi)
                source_priors_["phi_jl"] = "sampler_uniform"
                source_prior_params_["phi_jl"] = dict(min_=0, max_=2 * np.pi)

        # update the priors if input is given
        if source_priors:
            source_priors_.update(source_priors)
        if source_prior_params:
            source_prior_params_.update(source_prior_params)

        # taking care of source_prior_params from the available_gw_prior_list_and_its_params
        for key, value in source_priors_.items():
            if isinstance(value, str):
                dict_ = self.available_gw_prior_list_and_its_params[key]  # e.g. all source_frame_masses_prior function names and its parameters
                if value in dict_:
                    param_dict = dict_[value]
                    if source_prior_params_[key] is None:
                        source_prior_params_[key] = param_dict
                    else:
                        param_dict.update(source_prior_params_[key])
                        source_prior_params_[key] = param_dict
                else:
                    raise ValueError(
                        f"source_prior_params_['{key}'] is not in available_gw_prior_list_and_its_params"
                    )
            elif not callable(value):
                raise ValueError(
                    f"source_prior_params_['{key}'] should be either a string name of available sampler or a function"
                )

        return (source_priors_, source_prior_params_)
    
    def sample_gw_parameters(self, size=1000, param=None):
        """
        Function to sample BBH/BNS/NSBH intrinsic and extrinsics parameters.

        Parameters
        ----------
        size : `int`
            Number of samples to draw

        Returns
        ----------
        gw_parameters : `dict`
            Dictionary of sampled parameters
            gw_parameters.keys() = ['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'theta_jn', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl']

        Examples
        ----------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution()
        >>> params = cbc.gw_parameters(size=1000)
        >>> print("sampled parameters=",list(params.keys()))
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

    def binary_masses_BBH_popI_II_powerlaw_gaussian(
        self,
        size,
        get_attribute=False,
        **kwargs,
    ):
        """
        Function to sample source mass1 and mass2 with PowerLaw+PEAK model

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        mminbh : `float`
            Minimum mass of the black hole (Msun)
            default: 4.98
        mmaxbh : `float`
            Maximum mass of the black hole (Msun)
            default: 86.22
        alpha : `float`
            Spectral index for the powerlaw of the primary mass distribution
            default: 2.63
        mu_g : `float`
            Mean of the Gaussian component in the primary mass distribution
            default: 33.07
        sigma_g : `float`
            Width of the Gaussian component in the primary mass distribution
            default: 5.69
        lambda_peak : `float`
            Fraction of the model in the Gaussian component
            default: 0.10
        delta_m : `float`
            Range of mass tapering on the lower end of the mass distribution
            default: 4.82
        beta : `float`
            Spectral index for the powerlaw of the mass ratio distribution
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. param = dict(mminbh=4.98, mmaxbh=86.22, alpha=2.63, mu_g=33.07, sigma_g=5.69, lambda_peak=0.10, delta_m=4.82, beta=1.26)

        Returns
        ----------
        mass_1_source : `numpy.ndarray` (1D array of floats)
            Array of mass1 in source frame (Msun)
        mass_2_source : `numpy.ndarray` (1D array of floats)
            Array of mass2 in source frame (Msun)

        Examples
        ----------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution()
        >>> m1_src, m2_src = cbc.binary_masses_BBH_popI_II_powerlaw_gaussian(size=1000)
        """

        identifier_dict = {'name': "binary_masses_BBH_popI_II_powerlaw_gaussian"}
        param_dict = self.available_gw_prior_list_and_its_params["source_frame_masses"]["binary_masses_BBH_popI_II_powerlaw_gaussian"].copy()
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)


        # mass function
        rvs_ = lambda size: sample_powerlaw_gaussian_source_bbh_masses(
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
            param_dict_given=identifier_dict,
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
        # m_min=5.0,
        # m_max=150.0,
        # Mc=30.0,
        # sigma=0.3,
        # chunk_size=10000,
        get_attribute=False,
        **kwargs,
    ):
        """
        Function to sample source mass1 and mass2 with pop III origin. Refer to Eqn. 1 and 4 of Ng et al. 2022

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        m_min : `float`
            Minimum mass of the black hole (popIII) (Msun)
            default: 10.
        m_max : `float`
            Maximum mass of the black hole (popIII) (Msun)
            default: 100.
        Mc    : `float`
            Mass scale; the distribution is centered around Mc
            default: 30.0
        sigma : `float`
            Width of the distribution
            default: 0.3
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. param = dict(m_min=10., m_max=100., Mc=30.0, sigma=0.3)

        Returns
        ----------
        mass_1_source : `numpy.ndarray` (1D array of floats)
            Array of mass1 in source frame (Msun)
        mass_2_source : `numpy.ndarray` (1D array of floats)
            Array of mass2 in source frame (Msun)

        Examples
        ----------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution()
        >>> m1_src, m2_src = cbc.binary_masses_BBH_popIII_lognormal(size=1000)
        """

        identifier_dict = {'name': "binary_masses_BBH_popIII_lognormal"}
        param_dict = self.available_gw_prior_list_and_its_params["source_frame_masses"]["binary_masses_BBH_popIII_lognormal"].copy()
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        rvs_ = lambda size: lognormal_distribution_2D(
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
            param_dict_given=identifier_dict,
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
        # m_min=1.0,
        # m_max=100.0,
        # Mc=20.0,
        # sigma=0.3,
        # chunk_size=10000,
        get_attribute=False,
        **kwargs,
    ):
        """
        Function to sample source mass1 and mass2 with primordial origin. Refer to Eqn. 1 and 4 of Ng et al. 2022

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        m_min : `float`
            Minimum mass of the black hole (primordial) (Msun)
            default: 10.
        m_max : `float`
            Maximum mass of the black hole (primordial) (Msun)
            default: 100.
        Mc, sigma : `float`
            Fitting parameters
            default: Mc=30.0, sigma=0.3
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. param = dict(m_min=10., m_max=100., Mc=30.0, sigma=0.3)

        Returns
        ----------
        mass_1_source : `numpy.ndarray` (1D array of floats)
            Array of mass1 in source frame (Msun)
        mass_2_source : `numpy.ndarray` (1D array of floats)
            Array of mass2 in source frame (Msun)
        """

        identifier_dict = {'name': "binary_masses_BBH_primordial_lognormal"}    
        param_dict = self.available_gw_prior_list_and_its_params["source_frame_masses"]["binary_masses_BBH_primordial_lognormal"].copy()
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        rvs_ = lambda size: lognormal_distribution_2D(
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
            param_dict_given=identifier_dict,
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

    # def binary_masses_BNS_gwcosmo(
    #     self, size, get_attribute=False, **kwargs
    # ):
    #     """
    #     Function to calculate source mass1 and mass2 of BNS from powerlaw distribution (gwcosmo)

    #     Parameters
    #     ----------
    #     size : `int`
    #         Number of samples to draw
    #     mminns : `float`
    #         Minimum mass of the BNS (Msun)
    #         default: 1.0
    #     mmaxns : `float`
    #         Maximum mass of the BNS (Msun)
    #         default: 3.0
    #     alphans : `float`
    #         Power law index
    #         default: 0.0

    #     Returns
    #     ----------
    #     mass_1_source : `numpy.ndarray` (1D array of floats)
    #         Array of mass1 in source frame (Msun)
    #     mass_2_source : `numpy.ndarray` (1D array of floats)
    #         Array of mass2 in source frame (Msun)

    #     Examples
    #     ----------
    #     >>> from ler.gw_source_population import CBCSourceParameterDistribution
    #     >>> cbc = CBCSourceParameterDistribution()
    #     >>> m1_src, m2_src = cbc.binary_masses_BNS_gwcosmo(size=1000)
    #     """

    #     # if param:
    #     #     mminns = param["mminns"]
    #     #     mmaxns = param["mmaxns"]
    #     #     alphans = param["alphans"]

    #     identifier_dict = {'name': "binary_masses_BNS_gwcosmo"}
    #     param_dict = self.available_gw_prior_list_and_its_params["source_frame_masses"]["binary_masses_BNS_gwcosmo"].copy()
    #     if param_dict is not None:
    #         param_dict.update(kwargs)
    #     else:
    #         param_dict = kwargs
    #     identifier_dict.update(param_dict)

    #     # mass function for BNS
    #     model = p.BNS(
    #         mminns=identifier_dict["mminns"],
    #         mmaxns=identifier_dict["mmaxns"],
    #         alphans=identifier_dict["alphans"],
    #     )
    #     rvs_ = lambda size: model.sample(Nsample=size)

    #     mass_object = FunctionConditioning(
    #         function=None,
    #         x_array=None,
    #         param_dict_given=identifier_dict,
    #         directory=self.directory,
    #         sub_directory="source_frame_masses",
    #         name=identifier_dict['name'],
    #         create_new=self.create_new_interpolator["source_frame_masses"]["create_new"],
    #         create_function_inverse=False,
    #         create_function=False,
    #         create_pdf=False,
    #         create_rvs=rvs_,
    #         callback='rvs',
    #     )

    #     return mass_object if get_attribute else mass_object(size)
        
    def binary_masses_NSBH_broken_powerlaw(self, size, get_attribute=False, **kwargs):
        """
        Function to calculate source mass1 and mass2 of NSBH from powerlaw distribution (gwcosmo). Parameters are mminbh=26,mmaxbh=125,alpha_1=6.75,alpha_2=6.75,b=0.5,delta_m=5,mminns=1.0,mmaxns=3.0,alphans=0.0.

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        mminbh : `float`
            Minimum mass of the black hole (Msun)
            default: 26
        mmaxbh : `float`
            Maximum mass of the black hole (Msun)
            default: 125
        alpha_1 : `float`
            Power law index for the primary mass distribution
            default: 6.75
        alpha_2 : `float`
            Power law index for the secondary mass distribution
            default: 6.75
        b : `float`
            Break point of the power law
            default: 0.5
        delta_m : `float`
            Range of mass tapering on
            default: 5
        mminns : `float`
            Minimum mass of the neutron star (Msun)
            default: 1.0
        mmaxns : `float`
            Maximum mass of the neutron star (Msun)
            default: 3.0
        alphans : `float`
            Power law index for the neutron star mass distribution
            default: 0.0
        get_attribute : `bool`
            If True, return a sampler function with size as the only input where parameters are fixed to the given values.
        param : `dict`
            Allows to pass in above parameters as dict.

        Returns
        ----------
        mass_1_source : `numpy.ndarray` (1D array of floats)
            Array of mass1 in source frame (Msun)
        mass_2_source : `numpy.ndarray` (1D array of floats)
            Array of mass2 in source frame (Msun)

        Examples
        ----------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution()
        >>> m1_src, m2_src = cbc.binary_masses_NSBH_broken_powerlaw(size=1000)
        """

        identifier_dict = {'name': "binary_masses_NSBH_broken_powerlaw"}
        param_dict = self.available_gw_prior_list_and_its_params["source_frame_masses"]["binary_masses_NSBH_broken_powerlaw"].copy()
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        # mass function for NSBH
        rvs_ = lambda size: sample_broken_powerlaw_nsbh_masses(
            Nsample=size,
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
            param_dict_given=identifier_dict,
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
        Function to sample source mass1 and mass2 from uniform distribution.

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        m_min : `float`
            Minimum mass of the BNS
            default: 1.0
        m_max : `float`
            Maximum mass of the BNS
            default: 3.0
        get_attribute : `bool`
            If True, return a sampler function with size as the only input where parameters are fixed to the given values.
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. param = dict(m_min=1.0, m_max=3.0)

        Returns
        ----------
        mass_1_source : `numpy.ndarray` (1D array of floats)
            Array of mass1 in source frame (Msun)
        mass_2_source : `numpy.ndarray` (1D array of floats)
            Array of mass2 in source frame (Msun)

        Examples
        ----------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution()
        >>> m1_src, m2_src = cbc.binary_masses_uniform(size=1000)
        """

        identifier_dict = {'name': "binary_masses_uniform"}
        param_dict = self.available_gw_prior_list_and_its_params["source_frame_masses"]["binary_masses_uniform"].copy()
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        m_min = identifier_dict["m_min"]
        m_max = identifier_dict["m_max"]
        rvs_ = njit(lambda size: (np.random.uniform(m_min, m_max, size), np.random.uniform(m_min, m_max, size)))

        mass_object = FunctionConditioning(
            function=None,
            x_array=None,
            param_dict_given=identifier_dict,
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
        # w=0.643,
        # muL=1.352,
        # sigmaL=0.08,
        # muR=1.88,
        # sigmaR=0.3,
        # mmin=1.0,
        # mmax=2.3,
        # resolution=500,
        get_attribute=False,
        **kwargs,
    ):
        """
        Function to sample source mass1 and mass2 from bimodal distribution. Refer to Will M. Farr et al. 2020 Eqn. 6, https://arxiv.org/pdf/2005.00032.pdf .

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        w : `float`
            Weight of the left peak
            default: 0.643
        muL : `float`
            Mean of the left peak
            default: 1.352
        sigmaL : `float`
            Width of the left peak
            default: 0.08
        muR : `float`
            Mean of the right peak
            default: 1.88
        sigmaR : `float`
            Width of the right peak
            default: 0.3
        mmin : `float`
            Minimum mass of the BNS
            default: 1.0
        mmax : `float`
            Maximum mass of the BNS
            default: 2.3
        resolution : `int`
            Number of points to sample
            default: 500
        create_new : `bool`
            If True, create new interpolator
            default: False
        get_attribute : `bool`
            If True, return a sampler function with size as the only input where parameters are fixed to the given values.
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. param = dict(w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1.0, mmax=2.3, resolution=500)

        Returns
        ----------
        mass_1_source : `numpy.ndarray` (1D array of floats)
            Array of mass1 in source frame (Msun)
        mass_2_source : `numpy.ndarray` (1D array of floats)
            Array of mass2 in source frame (Msun)

        Examples
        ----------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution()
        >>> m1_src, m2_src = cbc.binary_masses_BNS_bimodal(size=1000)
        """

        identifier_dict = {'name': "binary_masses_BNS_bimodal"}
        identifier_dict['resolution'] = self.create_new_interpolator["source_frame_masses"]["resolution"]
        param_dict = self.available_gw_prior_list_and_its_params["source_frame_masses"]["binary_masses_BNS_bimodal"].copy()
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        # mass function for BNS
        mass = np.linspace(identifier_dict["mmin"], identifier_dict["mmax"], identifier_dict["resolution"]) 
        model = lambda mass: bns_bimodal_pdf(
            mass,
            w=identifier_dict["w"],
            muL=identifier_dict["muL"],
            sigmaL=identifier_dict["sigmaL"],
            muR=identifier_dict["muR"],
            sigmaR=identifier_dict["sigmaR"],
            mmin=identifier_dict["mmin"],
            mmax=identifier_dict["mmax"],
        )

        mass_object = FunctionConditioning(
            function=model,
            x_array=mass,
            param_dict_given=identifier_dict,
            directory=self.directory,
            sub_directory="source_frame_masses",
            name=identifier_dict['name'],
            create_new=self.create_new_interpolator["source_frame_masses"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='rvs',
        )

        cdf_values = mass_object.cdf_values
        x_array = mass_object.x_array

        mass_object.rvs = njit(lambda size: inverse_transform_sampler_m1m2(
            size, 
            cdf_values, 
            x_array,)
        )

        if get_attribute:
            return mass_object
        else:
            return mass_object(size)

    def constant_values_n_size(
        self, size=100, get_attribute=False, **kwargs
    ):
        """
        Function to sample constant values of size n.

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        get_attribute : `bool`
            If True, return the njitted sampler function with size as the only input where parameters are fixed to the given values.
        kwargs : `keyword arguments`
            Additional parameters to pass to the function


        Returns
        ----------
        values : `numpy.ndarray` (1D array of floats)
            Array of constant values

        Examples
        ----------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution()
        >>> value = cbc.constant_values_n_size(size=1000)
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
        Function to sample values from uniform distribution.

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        get_attribute : `bool`
            If True, return the njitted sampler function with size as the only input where parameters are fixed to the given values.
        param : `dict`
            Allows to pass in above parameters as dict.

        Returns
        ----------
        values : `numpy.ndarray` (1D array of floats)
            Array of uniformly distributed values in the range of [min_, max_]

        Examples
        ----------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution()
        >>> value = cbc.sampler_uniform(size=1000)
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
            param_dict_given=identifier_dict,
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
        Function to sample from sine distribution at the limit of [-np.pi/2, np.pi/2]

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        get_attribute : `bool`
            If True, return the njitted sampler function with size as the only input where parameters are fixed to the given values.
        param : None
            This parameter is not used. It is only here to make the function signature consistent with other samplers.

        Returns
        ----------
        sine : `numpy.ndarray` (1D array of floats)
            Array of values in the range of [-np.pi/2, np.pi/2]
        """

        identifier_dict = {}
        identifier_dict['name'] = "sampler_cosine"
        identifier_dict.update(kwargs)

        pdf_ = njit(lambda x: 0.5 * np.cos(x))
        rvs_ = njit(lambda size: np.arcsin((np.random.uniform(0, 1, size=size) * 2 - 1)))

        object_ = FunctionConditioning(
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
        Function to sample from sine distribution at the limit of [0, np.pi]

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        get_attribute : `bool`
            If True, return the njitted sampler function with size as the only input where parameters are fixed to the given values.
        param : None
            This parameter is not used. It is only here to make the function signature consistent with other samplers.

        Returns
        ----------
        sine : `numpy.ndarray` (1D array of floats)
            Array of values in the range of [0, np.pi]
        """

        identifier_dict = {}
        identifier_dict['name'] = "sampler_sine"
        identifier_dict.update(kwargs)

        pdf_ = njit(lambda x: 0.5 * np.sin(x))
        rvs_ = njit(lambda size: np.arccos((np.random.uniform(0, 1, size=size) - 0.5) * 2))

        object_ = FunctionConditioning(
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
        Function to sample source frame masses (mass1_source, mass2_source) with the initialized prior.

        Parameters
        ----------
        size : `int`
            Number of samples to draw

        Returns
        ----------
        mass_1_source : `numpy.ndarray` (1D array of floats)
            Array of mass1 in source frame
        mass_2_source : `numpy.ndarray` (1D array of floats)
            Array of mass2 in source frame
        """

        return self._source_frame_masses

    @source_frame_masses.setter
    def source_frame_masses(self, prior):
        if prior in self.available_gw_prior_list_and_its_params["source_frame_masses"]:
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
        elif callable(prior):
            print("using user defined custom source_frame_masses function")
            self._source_frame_masses = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
        elif isinstance(prior, object):
            print("using user defined custom source_frame_masses class/object")
            self._source_frame_masses = prior
        else:
            raise ValueError(
                "source_frame_masses prior not available in available_gw_prior_list_and_its_params. Must be a string or a callable function."
            )
        
    @property
    def zs(self):
        """
        Function to sample source redshift with the initialized prior.

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        
        Returns
        ----------
        zs : `numpy.ndarray` (1D array of floats)
            Array of source redshift
        """

        return self._zs
    
    @zs.setter
    def zs(self, prior):
        if prior in self.available_gw_prior_list_and_its_params["zs"]:
            print(f"using ler available zs function : {prior}")
            self._zs = getattr(self, prior)
        elif callable(prior):
            print("using user defined custom zs function")
            self._zs = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
        elif isinstance(prior, object):
            print("using user defined custom zs class/object")
        else:
            raise ValueError(
                "zs prior not available in available_gw_prior_list_and_its_params. Must be a string or a callable function."
            )

    @property
    def geocent_time(self):
        """
        Function to sample geocent time with the initialized prior.

        Parameters
        ----------
        size : `int`
            Number of samples to draw

        Returns
        ----------
        geocent_time : `numpy.ndarray` (1D array of floats)
            Array of geocent_time or time of coalescence
        """

        return self._geocent_time

    @geocent_time.setter
    def geocent_time(self, prior):
        if prior in self.available_gw_prior_list_and_its_params["geocent_time"]:
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
        elif callable(prior):
            print("using user defined custom geocent_time function")
            self._geocent_time = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
        elif isinstance(prior, object):
            print("using user defined custom geocent_time class/object")
            self._geocent_time = prior
        else:
            raise ValueError(
                "geocent_time prior not available in available_gw_prior_list_and_its_params. Must be a string or a callable function."
            )

    @property
    def ra(self):
        """
        Function to sample right ascension of sky position with the initialized prior.

        Parameters
        ----------
        size : `int`
            Number of samples to draw

        Returns
        ----------
        ra : `numpy.ndarray` (1D array of floats)
            Array of right ascension of sky position
        """

        return self._ra

    @ra.setter
    def ra(self, prior):
        if prior in self.available_gw_prior_list_and_its_params["ra"]:
            args = self.gw_param_samplers_params["ra"]
            if args is None:
                self._ra = getattr(self, prior)(size=None, get_attribute=True)
            else:
                # follwing should return a sampler function with only one argument (size)
                self._ra = getattr(self, prior)(size=None, get_attribute=True, **args)
        elif callable(prior):
            print("using user defined custom ra function")
            self._ra = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
        elif isinstance(prior, object):
            print("using user defined custom ra class/object")
            self._ra = prior
        else:
            raise ValueError(
                "ra prior not available in available_gw_prior_list_and_its_params. Must be a string or a callable function."
            )

    @property
    def dec(self):
        """
        Function to sample declination of sky position with the initialized prior.

        Parameters
        ----------
        size : `int`
            Number of samples to draw

        Returns
        ----------
        dec : `numpy.ndarray` (1D array of floats)
            Array of declination of sky position
        """

        return self._dec

    @dec.setter
    def dec(self, prior):
        if prior in self.available_gw_prior_list_and_its_params["dec"]:
            print(f"using ler available dec function : {prior}")
            args = self.gw_param_samplers_params["dec"]
            if args is None:
                self._dec = getattr(self, prior)(size=None, get_attribute=True)
            else:
                # follwing should return a sampler function with only one argument (size)
                self._dec = getattr(self, prior)(size=None, get_attribute=True, **args)
        elif callable(prior):
            print("using user provided custom dec function")
            self._dec = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
        elif isinstance(prior, object):
            print("using user provided custom dec class/object")
            self._dec = prior
        else:
            raise ValueError(
                "dec prior not available in available_gw_prior_list_and_its_params. Must be a string or a callable function."
            )

    @property
    def phase(self):
        """
        Function to sample coalescence phase with the initialized prior.

        Parameters
        ----------
        size : `int`
            Number of samples to draw

        Returns
        ----------
        phase : `numpy.ndarray` (1D array of floats)
            Array of coalescence phase
        """

        return self._phase

    @phase.setter
    def phase(self, prior):
        if prior in self.available_gw_prior_list_and_its_params["phase"]:
            print(f"using ler available phase function : {prior}")
            args = self.gw_param_samplers_params["phase"]
            if args is None:
                self._phase = getattr(self, prior)(size=None, get_attribute=True)
            else:
                # follwing should return a sampler function with only one argument (size)
                self._phase = getattr(self, prior)(size=None, get_attribute=True, **args)
        elif callable(prior):
            print("using user provided custom phase function")
            self._phase = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
        elif isinstance(prior, object):
            print("using user provided custom phase class/object")
            self._phase = prior
        else:
            raise ValueError(
                "phase prior not available in available_gw_prior_list_and_its_params. Must be a string or a callable function."
            )

    @property
    def psi(self):
        """
        Function to sample polarization angle with the initialized prior.

        Parameters
        ----------
        size : `int`
            Number of samples to draw

        Returns
        ----------
        psi : `numpy.ndarray` (1D array of floats)
            Array of polarization angle
        """

        return self._psi

    @psi.setter
    def psi(self, prior):
        if prior in self.available_gw_prior_list_and_its_params["psi"]:
            print(f"using ler available psi function : {prior}")
            args = self.gw_param_samplers_params["psi"]
            if args is None:
                self._psi = getattr(self, prior)(size=None, get_attribute=True)
            else:
                # follwing should return a sampler function with only one argument (size)
                self._psi = getattr(self, prior)(size=None, get_attribute=True, **args)
        elif callable(prior):
            print("using user provided custom psi function")
            self._psi = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
        elif isinstance(prior, object):
            print("using user provided custom psi class/object")
            self._psi = prior
        else:
            raise ValueError(
                "psi prior not available in available_gw_prior_list_and_its_params. Must be a string or a callable function."
            )

    @property
    def theta_jn(self):
        """
        Function to sample theta_jn with the initialized prior.

        Parameters
        ----------
        size : `int`
            Number of samples to draw

        Returns
        ----------
        theta_jn : `numpy.ndarray` (1D array of floats)
            Array of theta_jn
        """
        return self._theta_jn

    @theta_jn.setter
    def theta_jn(self, prior):
        if prior in self.available_gw_prior_list_and_its_params["theta_jn"]:
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
        elif callable(prior):
            print("using user provided custom theta_jn function")
            self._theta_jn = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
        elif isinstance(prior, object):
            print("using user provided custom theta_jn class/object")
            self._theta_jn = prior
        else:
            raise ValueError("theta_jn prior not available in available_gw_prior_list_and_its_params. Must be a string or a callable function.")

    @property
    def a_1(self):
        """
        Function to sample spin magnitude of the compact binaries (body1) with the initialized prior.

        Parameters
        ----------
        size : `int`
            Number of samples to draw

        Returns
        ----------
        a_1 : `numpy.ndarray` (1D array of floats)
            Array of spin magnitude of the compact binaries (body1)
        """
        return self._a_1

    @a_1.setter
    def a_1(self, prior):
        if prior in self.available_gw_prior_list_and_its_params["a_1"]:
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
        elif callable(prior):
            print("using user provided custom a_1 function")
            self._a_1 = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
        elif isinstance(prior, object):
            print("using user provided custom a_1 class/object")
            self._a_1 = prior
        else:
            raise ValueError("a_1 prior not available in available_gw_prior_list_and_its_params. Must be a string or a callable function.")

    @property
    def a_2(self):
        """
        Function to sample spin magnitude of the compact binaries (body2) with the initialized prior.

        Parameters
        ----------
        size : `int`
            Number of samples to draw

        Returns
        ----------
        a_2 : `numpy.ndarray` (1D array of floats)
            Array of spin magnitude of the compact binaries (body2)
        """
        return self._a_2

    @a_2.setter
    def a_2(self, prior):
        if prior in self.available_gw_prior_list_and_its_params["a_2"]:
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
        elif callable(prior):
            print("using user provided custom a_2 function")
            self._a_2 = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
        elif isinstance(prior, object):
            print("using user provided custom a_2 class/object")
            self._a_2 = prior
        else:
            raise ValueError("a_2 prior not available in available_gw_prior_list_and_its_params. Must be a string or a callable function.")

    @property
    def tilt_1(self):
        """
        Function to sample tilt angle of the compact binaries (body1) with the initialized prior.

        Parameters
        ----------
        size : `int`
            Number of samples to draw

        Returns
        ----------
        tilt_1 : `numpy.ndarray` (1D array of floats)
            Array of tilt angle of the compact binaries (body1)
        """
        return self._tilt_1

    @tilt_1.setter
    def tilt_1(self, prior):
        if prior in self.available_gw_prior_list_and_its_params["tilt_1"]:
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
        elif callable(prior):
            print("using user provided custom tilt_1 function")
            self._tilt_1 = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
        elif isinstance(prior, object):
            print("using user provided custom tilt_1 class/object")
            self._tilt_1 = prior
        else:
            raise ValueError("tilt_1 prior not available in available_gw_prior_list_and_its_params. Must be a string or a callable function.")

    @property
    def tilt_2(self):
        """
        Function to sample tilt angle of the compact binaries (body2) with the initialized prior.

        Parameters
        ----------
        size : `int`
            Number of samples to draw

        Returns
        ----------
        tilt_2 : `numpy.ndarray` (1D array of floats)
            Array of tilt angle of the compact binaries (body2)
        """

        return self._tilt_2

    @tilt_2.setter
    def tilt_2(self, prior):
        if prior in self.available_gw_prior_list_and_its_params["tilt_2"]:
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
        elif callable(prior):
            print("using user provided custom tilt_2 function")
            self._tilt_2 = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
        elif isinstance(prior, object):
            print("using user provided custom tilt_2 class/object")
            self._tilt_2 = prior
        else:
            raise ValueError("tilt_2 prior not available in available_gw_prior_list_and_its_params. Must be a string or a callable function.")
        
    @property
    def phi_12(self):
        """
        Function to sample azimuthal angle between the two spins with the initialized prior.

        Parameters
        ----------
        size : `int`
            Number of samples to draw

        Returns
        ----------
        phi_12 : `numpy.ndarray` (1D array of floats)
            Array of azimuthal angle between the two spins
        """

        return self._phi_12

    @phi_12.setter
    def phi_12(self, prior):
        if prior in self.available_gw_prior_list_and_its_params["phi_12"]:
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
        elif callable(prior):
            print("using user provided custom phi_12 function")
            self._phi_12 = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
        elif isinstance(prior, object):
            print("using user provided custom phi_12 class/object")
            self._phi_12 = prior
        else:
            raise ValueError("phi_12 prior not available in available_gw_prior_list_and_its_params. Must be a string or a callable function.")

    @property
    def phi_jl(self):
        """
        Function to sample azimuthal angle between the total angular momentum and the orbital angular momentum with the initialized prior.

        Parameters
        ----------
        size : `int`
            Number of samples to draw

        Returns
        ----------
        phi_jl : `numpy.ndarray` (1D array of floats)
            Array of azimuthal angle between the total angular momentum and the orbital angular momentum
        """
        return self._phi_jl

    @phi_jl.setter
    def phi_jl(self, prior):
        if prior in self.available_gw_prior_list_and_its_params["phi_jl"]:
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
        elif callable(prior):
            print("using user provided custom phi_jl function")
            self._phi_jl = FunctionConditioning(function=None, x_array=None, create_rvs=prior)
        elif isinstance(prior, object):
            print("using user provided custom phi_jl class/object")
            self._phi_jl = prior
        else:
            raise ValueError("phi_jl prior not available in available_gw_prior_list_and_its_params. Must be a string or a callable function.")

    @property
    def available_gw_prior_list_and_its_params(self):
        """
        Dictionary with list all the available priors and it's corresponding parameters. This is an immutable instance attribute.

        Examples
        ----------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution()
        >>> priors = cbc.available_gw_prior_list_and_its_params
        >>> priors.keys()  # type of priors
        dict_keys(['merger_rate_density', 'source_frame_masses', 'spin', 'geocent_time', 'ra', 'phase', 'psi', 'theta_jn'])
        >>> priors['source_frame_masses'].keys()  # type of source_frame_masses priors
        dict_keys(['binary_masses_BBH_popI_II_powerlaw_gaussian', 'binary_masses_BBH_popIII_lognormal', 'binary_masses_BBH_primordial_lognormal', 'binary_masses_BNS_bimodal'])
        >>> priors['source_frame_masses']['binary_masses_BBH_popI_II_powerlaw_gaussian'].keys()  # parameters of binary_masses_BBH_popI_II_powerlaw_gaussian
        dict_keys(['mminbh', 'mmaxbh', 'alpha', 'mu_g', 'sigma_g', 'lambda_peak', 'delta_m', 'beta'])
        """

        self._available_gw_prior_list_and_its_params = dict(
            merger_rate_density=self.merger_rate_density_model_list,
            zs=dict(
                source_redshift=None,
            ),
            source_frame_masses=dict(
                binary_masses_BBH_popI_II_powerlaw_gaussian=dict(
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
                    Mc=30.0, sigma=0.3, beta=1.1
                ),
                binary_masses_NSBH_broken_powerlaw=dict(
                    mminbh=26,
                    mmaxbh=125,
                    alpha_1=6.75,
                    alpha_2=0.0,
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
            ) if self.spin_zero else dict(
                constant_values_n_size=dict(value=0.0),
                sampler_uniform=dict(xmin=0.0, xmax=0.8),
            ),
            a_2=dict(
                constant_values_n_size=dict(value=0.0),
                sampler_uniform=dict(xmin=-0.8, xmax=0.8),
            ) if self.spin_zero else dict(
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

        return self._available_gw_prior_list_and_its_params