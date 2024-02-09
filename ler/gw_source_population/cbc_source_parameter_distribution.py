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

# for generating mass distribution
from gwcosmo import priors as p

# for multiprocessing
# Import helper routines
from ..utils import interpolator_from_pickle, cubic_spline_interpolator

# import redshift distribution sampler
from .cbc_source_redshift_distribution import CBCSourceRedshiftDistribution

from .jit_functions import lognormal_distribution_2D, inverse_transform_sampler_m1m2


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
        >>> from ler.gw_source_population import CompactBinaryPopulation
        >>> cbc = CompactBinaryPopulation()
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
        If spin_zero=True and spin_precession=True, spin parameters are sampled for precessing binaries.
        if spin_zero=True and spin_precession=False, spin parameters are sampled for aligned/anti-aligned spin binaries.
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
    >>> params = cbc.sample_gw_parameters(size=1000)
    >>> print("sampled parameters=",list(params.keys()))

    Instance Attributes
    ----------
    CompactBinaryPopulation has the following instance attributes:\n
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
    CompactBinaryPopulation has the following instance methods:\n
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
    |:meth:`~sample_gw_parameters`        | Function to sample all the       |
    |                                     | intrinsic and extrinsic          |
    |                                     | parameters of compact binaries   |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_source_frame_masses`  | Function to sample source mass1  |
    |                                     | and mass2                        |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_geocent_time`         | Function to sample geocent time  |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_zs`                   | Function to sample source        |
    |                                     | redshift                         |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_ra`                   | Function to sample right         |
    |                                     | ascension (sky position)         |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_dec`                  | Function to sample declination   |
    |                                     | (sky position)                   |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_phase`                | Function to sample coalescence   |
    |                                     | phase                            |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_psi`                  | Function to sample polarization  |
    |                                     | angle                            |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_theta_jn`             | Function to sample inclination   |
    |                                     | angle                            |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_a1`                   | Function to sample spin1         |
    |                                     | magnitude                        |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_a2`                   | Function to sample spin2         |
    |                                     | magnitude                        |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_tilt_1`               | Function to sample tilt1 angle   |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_tilt_2`               | Function to sample tilt2 angle   |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_phi_12`               | Function to sample phi12 angle   |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_phi_jl`               | Function to sample phi_jl angle  |
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
    |:meth:`~binary_masses_BNS_gwcosmo`                                      |
    +-------------------------------------+----------------------------------+
    |                                     | Function to sample source mass1  |
    |                                     | and mass2 from powerlaw          |
    |                                     | distribution.                    |
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
        # initialize the interpolator's parameters
        self.create_new_interpolator = dict(
            redshift_distribution=dict(create_new=False, resolution=500),
            z_to_luminosity_distance=dict(create_new=False, resolution=500),
            differential_comoving_volume=dict(create_new=False, resolution=500),
        )
        if isinstance(create_new_interpolator, dict):
            self.create_new_interpolator.update(create_new_interpolator)
        elif create_new_interpolator is True:
            self.create_new_interpolator = dict(
                redshift_distribution=dict(create_new=True, resolution=500),
                z_to_luminosity_distance=dict(create_new=True, resolution=500),
                differential_comoving_volume=dict(create_new=True, resolution=500),
            )

        # dealing with prior functions and categorization
        (
            self.gw_param_samplers,
            self.gw_param_samplers_params,
            self.sampler_names,
        ) = self.source_priors_categorization(
            event_type, source_priors, source_priors_params
        )

        if self.gw_param_samplers["zs"] == "sample_source_redshift":
            # initialize the SourceGalaxyPopulationModel mother class
            # for redshift distribution
            # instance attribute sample_source_redshift is initialized here
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
                create_new_interpolator=self.create_new_interpolator,
            )
        else:
            # if you already have the redshift distribution function, you can pass it in.
            # super class will not be initialized anymore,
            # get z_to_luminosity_distance
            self.lookup_table_luminosity_distance(z_min, z_max, directory)

        # initializing samplers
        # it goes through the setter functions and assign the sampler functions
        self.sample_source_frame_masses = self.gw_param_samplers["source_frame_masses"]
        self.sample_geocent_time = self.gw_param_samplers["geocent_time"]
        self.sample_zs = self.gw_param_samplers["zs"]
        self.sample_ra = self.gw_param_samplers["ra"]
        self.sample_dec = self.gw_param_samplers["dec"]
        self.sample_phase = self.gw_param_samplers["phase"]
        self.sample_psi = self.gw_param_samplers["psi"]
        self.sample_theta_jn = self.gw_param_samplers["theta_jn"]
        # initialize the spin prior attribute
        # remove spin prior if spin_zero is True
        if not spin_zero:
            self.sample_a_1 = self.gw_param_samplers["a_1"]
            self.sample_a_2 = self.gw_param_samplers["a_2"]
            if spin_precession:
                self.sample_tilt_1 = self.gw_param_samplers["tilt_1"]
                self.sample_tilt_2 = self.gw_param_samplers["tilt_2"]
                self.sample_phi_12 = self.gw_param_samplers["phi_12"]
                self.sample_phi_jl = self.gw_param_samplers["phi_jl"]

    def lookup_table_luminosity_distance(self, z_min, z_max, directory):
        """
        Function to create a lookup table for the differential comoving volume
        and luminosity distance wrt redshift.

        Parameters
        ----------
        z_min : `float`
            Minimum redshift of the source population
        z_max : `float`
            Maximum redshift of the source population

        Attributes
        ----------
        z_to_luminosity_distance : `scipy.interpolate.interpolate`
            Function to convert redshift to luminosity distance
        differential_comoving_volume : `scipy.interpolate.interpolate`
            Function to calculate the differential comoving volume
        """

        # initialing cosmological functions for fast calculation through interpolation
        resolution = self.c_n_i["z_to_luminosity_distance"]["resolution"]
        create_new = self.c_n_i["z_to_luminosity_distance"]["create_new"]
        spline1 = interpolator_from_pickle(
            param_dict_given=dict(
                z_min=z_min, z_max=z_max, cosmology=self.cosmo, resolution=resolution
            ),
            directory=directory,
            sub_directory="z_to_luminosity_distance",
            name="z_to_luminosity_distance",
            x=np.linspace(z_min, z_max, resolution),
            pdf_func=lambda z_: self.cosmo.luminosity_distance(z_).value,
            conditioned_y=None,
            dimension=1,
            category="function",
            create_new=create_new,
        )
        self.z_to_luminosity_distance = njit(
            lambda z_: cubic_spline_interpolator(z_, spline1[0], spline1[1])
        )
        self.z_to_luminosity_distance.__doc__ = """
        Function to convert redshift to luminosity distance.
        
        Parameters
        ----------
        zs : `numpy.ndarray`
            1D array of floats
            Source redshifts

        Returns
        ----------
        luminosity_distance : `numpy.ndarray`
            1D array of floats
            luminosity distance
        """

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
        >>> params = cbc.sample_gw_parameters(size=1000)
        >>> print("sampled parameters=",list(params.keys()))
        """

        # check for input parameters
        if param is None:
            param = {}
            param_keys = param.keys()
        else:
            param_keys = param.keys()

        # sample parameters
        param_names = list(self.gw_param_samplers.keys())
        del param_names[0]  # remove merger_rate_density
        # make sure the order is correct
        sampler_names = list(self.sampler_names.keys())

        gw_parameters = {}  # initialize dictionary to store parameters
        for name, sampler in zip(param_names, sampler_names):
            # print(name)
            if name not in param_keys:
                # Sample the parameter using the specified sampler function
                gw_parameters[name] = getattr(self, sampler)(size)
            else:
                # Use the provided value from kwargs
                gw_parameters[name] = param[name]

        # calculate luminosity distance
        zs = gw_parameters["zs"]
        gw_parameters["luminosity_distance"] = self.z_to_luminosity_distance(zs)  # Mpc
        
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
        mminbh=4.98,
        mmaxbh=112.5,
        alpha=3.78,
        mu_g=32.27,
        sigma_g=3.88,
        lambda_peak=0.03,
        delta_m=4.8,
        beta=0.81,
        get_attribute=False,
        param=None,
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

        if param:
            mminbh = param["mminbh"]
            mmaxbh = param["mmaxbh"]
            alpha = param["alpha"]
            mu_g = param["mu_g"]
            sigma_g = param["sigma_g"]
            lambda_peak = param["lambda_peak"]
            delta_m = param["delta_m"]
            beta = param["beta"]

        # check if the cdf exist

        # mass function
        model = p.BBH_powerlaw_gaussian(
            mminbh=mminbh,
            mmaxbh=mmaxbh,
            alpha=alpha,
            mu_g=mu_g,
            sigma_g=sigma_g,
            lambda_peak=lambda_peak,
            delta_m=delta_m,
            beta=beta,
        )
        if get_attribute:
            sampler_function = lambda size_: model.sample(Nsample=size_)
            return sampler_function
        else:
            # sample mass1 and mass2
            mass_1_source, mass_2_source = model.sample(Nsample=size)

            return (mass_1_source, mass_2_source)

    def binary_masses_BBH_popIII_lognormal(
        self,
        size,
        m_min=5.0,
        m_max=150.0,
        Mc=30.0,
        sigma=0.3,
        chunk_size=10000,
        get_attribute=False,
        param=None,
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

        if param:
            m_min = param["m_min"]
            m_max = param["m_max"]
            Mc = param["Mc"]
            sigma = param["sigma"]

        if get_attribute:
            return njit(
                lambda size_: lognormal_distribution_2D(
                    size_,
                    m_min=m_min,
                    m_max=m_max,
                    Mc=Mc,
                    sigma=sigma,
                    chunk_size=chunk_size,
                )
            )
        else:
            return lognormal_distribution_2D(
                size,
                m_min=m_min,
                m_max=m_max,
                Mc=Mc,
                sigma=sigma,
                chunk_size=chunk_size,
            )

    def binary_masses_BBH_primordial_lognormal(
        self,
        size,
        m_min=1.0,
        m_max=100.0,
        Mc=20.0,
        sigma=0.3,
        chunk_size=10000,
        get_attribute=False,
        param=None,
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

        if param:
            m_min = param["m_min"]
            m_max = param["m_max"]
            Mc = param["Mc"]
            sigma = param["sigma"]

        if get_attribute:
            return njit(
                lambda size_: lognormal_distribution_2D(
                    size_,
                    m_min=m_min,
                    m_max=m_max,
                    Mc=Mc,
                    sigma=sigma,
                    chunk_size=chunk_size,
                )
            )
        else:
            return lognormal_distribution_2D(
                size,
                m_min=m_min,
                m_max=m_max,
                Mc=Mc,
                sigma=sigma,
                chunk_size=chunk_size,
            )

    def binary_masses_BNS_gwcosmo(
        self, size, mminns=1.0, mmaxns=3.0, alphans=0.0, get_attribute=False, param=None
    ):
        """
        Function to calculate source mass1 and mass2 of BNS from powerlaw distribution (gwcosmo)

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        mminns : `float`
            Minimum mass of the BNS (Msun)
            default: 1.0
        mmaxns : `float`
            Maximum mass of the BNS (Msun)
            default: 3.0
        alphans : `float`
            Power law index
            default: 0.0

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
        >>> m1_src, m2_src = cbc.binary_masses_BNS_gwcosmo(size=1000)
        """

        if param:
            mminns = param["mminns"]
            mmaxns = param["mmaxns"]
            alphans = param["alphans"]

        # mass function for BNS
        model = p.BNS(mminns=mminns, mmaxns=mmaxns, alphans=alphans)

        if get_attribute:
            sampler_function = lambda size_: model.sample(Nsample=size_)
            return sampler_function
        else:
            # sampling
            mass_1_source, mass_2_source = model.sample(Nsample=size)

            return (mass_1_source, mass_2_source)

    def binary_masses_BNS_bimodal(
        self,
        size,
        w=0.643,
        muL=1.352,
        sigmaL=0.08,
        muR=1.88,
        sigmaR=0.3,
        mmin=1.0,
        mmax=2.3,
        resolution=500,
        create_new=False,
        get_attribute=False,
        param=None,
    ):
        """
        Function to sample source mass1 and mass2 from bimodal distribution. Refer to Will M. Farr et al. 2020 Eqn. 6

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

        if param:
            w, muL, sigmaL, muR, sigmaR, mmin, mmax = (
                param["w"],
                param["muL"],
                param["sigmaL"],
                param["muR"],
                param["sigmaR"],
                param["mmin"],
                param["mmax"],
            )
        # left and right peak
        pdf_unnormL = lambda m: np.exp(-((m - muL) ** 2) / (2 * sigmaL**2))
        normL = quad(pdf_unnormL, mmin, mmax)[0]  # normalization constant
        pdf_unnormR = lambda m: np.exp(-((m - muR) ** 2) / (2 * sigmaR**2))
        normR = quad(pdf_unnormR, mmin, mmax)[0]  # normalization constant
        # total pdf
        pdf = lambda m: w * pdf_unnormL(m) / normL + (1 - w) * pdf_unnormR(m) / normR

        # find inverse cdf
        inv_cdf = interpolator_from_pickle(
            param_dict_given=dict(
                w=0.643,
                muL=1.352,
                sigmaL=0.08,
                muR=1.88,
                sigmaR=0.3,
                mmin=1.0,
                mmax=2.3,
                resolution=500,
            ),
            directory=self.directory,
            sub_directory="binary_masses_BNS_bimodal",
            name="binary_masses_BNS_bimodal",
            x=np.linspace(mmin, mmax, resolution),
            pdf_func=pdf,
            conditioned_y=None,
            dimension=1,
            category="inv_cdf",
            create_new=create_new,
        )

        # sample from inverse cdf
        if get_attribute:
            return njit(
                lambda size_: inverse_transform_sampler_m1m2(
                    size_, inv_cdf[0], inv_cdf[1]
                )
            )
        else:
            return inverse_transform_sampler_m1m2(size, inv_cdf[0], inv_cdf[1])

    def constant_values_n_size(
        self, size=100, value=0.0, get_attribute=False, param=None
    ):
        """
        Function to sample constant values of size n.

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        value : `float`
            Constant value
            default: 0.0
        get_attribute : `bool`
            If True, return the njitted sampler function with size as the only input where parameters are fixed to the given values.
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. param = dict(value=0.0)

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

        if param:
            value = param["value"]

        if get_attribute:
            return njit(lambda size_: np.ones(size_) * value)
        else:
            return np.ones(size) * value

    def sampler_uniform(
        self, size, min_=0, max_=np.pi, get_attribute=False, param=None
    ):
        """
        Function to sample values from uniform distribution.

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        start_time : `float`
            Start time of the uniform distribution
            default: 1238166018
        end_time : `float`
            End time of the uniform distribution
            default: 1238166018 + 31536000
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

        if param:
            min_ = param["min_"]
            max_ = param["max_"]

        if get_attribute:
            return njit(lambda size_: np.random.uniform(min_, max_, size=size_))
        else:
            return np.random.uniform(min_, max_, size=size)

    def sampler_cosine(
        self, size, get_attribute=False, param=None
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

        if get_attribute:
            return njit(
                lambda size_: np.arcsin((np.random.uniform(0, 1, size=size_) * 2 - 1))
            )
        else:
            return np.arcsin((np.random.uniform(0, 1, size=size) * 2 - 1))

    def sampler_sine(
        self, size, get_attribute=False, param=None
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

        if get_attribute:
            return njit(
                lambda size_: np.arccos((np.random.uniform(0, 1, size=size_) - 0.5) * 2)
            )
        else:
            return np.arccos((np.random.uniform(0, 1, size=size) - 0.5) * 2)

    @property
    def available_gw_prior_list_and_its_params(self):
        """
        Dictionary with list all the available priors and it's corresponding parameters. This is an immutable instance attribute.

        Examples
        ----------
        >>> from ler.gw_source_population import CompactBinaryPopulation
        >>> cbc = CompactBinaryPopulation()
        >>> priors = cbc.available_gw_prior_list_and_its_params
        >>> priors.keys()  # type of priors
        dict_keys(['merger_rate_density', 'source_frame_masses', 'spin', 'geocent_time', 'ra', 'phase', 'psi', 'theta_jn'])
        >>> priors['source_frame_masses'].keys()  # type of source_frame_masses priors
        dict_keys(['binary_masses_BBH_popI_II_powerlaw_gaussian', 'binary_masses_BBH_popIII_lognormal', 'binary_masses_BBH_primordial_lognormal', 'binary_masses_BNS_gwcosmo', 'binary_masses_BNS_bimodal'])
        >>> priors['source_frame_masses']['binary_masses_BBH_popI_II_powerlaw_gaussian'].keys()  # parameters of binary_masses_BBH_popI_II_powerlaw_gaussian
        dict_keys(['mminbh', 'mmaxbh', 'alpha', 'mu_g', 'sigma_g', 'lambda_peak', 'delta_m', 'beta'])
        """

        self._available_gw_prior_list_and_its_params = dict(
            merger_rate_density=self.merger_rate_density_model_list,
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
                binary_masses_BBH_popIII_lognormal=dict(Mc=30.0, sigma=0.3, beta=1.1),
                binary_masses_BBH_primordial_lognormal=dict(
                    Mc=30.0, sigma=0.3, beta=1.1
                ),
                binary_masses_BNS_gwcosmo=dict(mminns=1.0, mmaxns=3.0, alphans=0.0),
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
            zs=dict(
                sample_source_redshift=dict(zs=None),
            ),
            spin=dict(
                constant_values_n_size=dict(value=0.0),
                binary_spin_BBH_bilby=None,
                binary_spin_BNS_bilby=None,
            ),
            geocent_time=dict(
                geocent_time_uniform=dict(
                    start_time=1238166018, end_time=1238166018 + 31536000
                )
            ),
            ra=dict(ra_uniform_bilby=None),
            phase=dict(phase_uniform_bilby=None),
            psi=dict(psi_uniform_bilby=None),
            theta_jn=dict(theta_jn_uniform_bilby=None),
        )

        return self._available_gw_prior_list_and_its_params

    def source_priors_categorization(
        self, event_type, source_priors, event_prior_params
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
        event_prior_params : `dict`
            Dictionary of sampler parameters for each GW parameter

        Returns
        ----------
        source_priors_ : `dict`
            Dictionary of prior sampler functions for each parameter
        event_prior_params_ : `dict`
            Dictionary of sampler parameters for each parameter
        sampler_names_ : `dict`
            Dictionary of sampler names with description

        Examples
        ----------
        >>> from ler.gw_source_population import CBCSourceParameterDistribution
        >>> cbc = CBCSourceParameterDistribution()
        >>> source_priors, event_prior_params, sampler_names = cbc.source_priors_categorization(event_type='BBH', source_priors=None, event_prior_params=None)
        >>> print(source_priors.keys())
        >>> print(event_prior_params.keys())
        >>> print(sampler_names.keys())
        """

        # for BBH
        if event_type == "BBH":
            merger_rate_density_prior = "merger_rate_density_bbh_popI_II_oguri2018"
            merger_rate_density_prior_params = dict(
                R0=23.9 * 1e-9, b2=1.6, b3=2.0, b4=30
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
                R0=105.5 * 1e-9, b2=1.6, b3=2.0, b4=30
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
                R0=27.0 * 1e-9, b2=1.6, b3=2.0, b4=30
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
                m_min=5.0, m_max=150.0, Mc=30.0, sigma=0.3, chunk_size=10000
            )
            a_max = 0.8

        elif event_type == "BBH_primordial":
            merger_rate_density_prior = "merger_rate_density_bbh_primordial_ken2022"
            merger_rate_density_prior_params = dict(
                n0=0.044 * 1e-9, t0=13.786885302009708
            )
            source_frame_masses_prior = "binary_masses_BBH_primordial_lognormal"
            source_frame_masses_prior_params = dict(
                m_min=1.0, m_max=100.0, Mc=20.0, sigma=0.3, chunk_size=10000
            )
            a_max = 0.8

        else:
            raise ValueError("event_type is not recognized")

        # setting the priors and its parameters
        source_priors_ = dict(
            merger_rate_density=merger_rate_density_prior,
            source_frame_masses=source_frame_masses_prior,
            zs="sample_source_redshift",
            geocent_time="sampler_uniform",
            ra="sampler_uniform",
            dec="sampler_cosine",
            phase="sampler_uniform",
            psi="sampler_uniform",
            theta_jn="sampler_sine",
        )
        event_prior_params_ = dict(
            merger_rate_density=merger_rate_density_prior_params,
            source_frame_masses=source_frame_masses_prior_params,
            zs=None,
            geocent_time=dict(min_=1238166018, max_=1269702018),
            ra=dict(min_=0., max_=2.*np.pi),
            dec=None,
            phase=dict(min_=0., max_=2.*np.pi),
            psi=dict(min_=0., max_=np.pi),
            theta_jn=None,
        )

        # dict of sampler names with description
        sampler_names_ = dict(
            sample_source_frame_masses="samples mass1 and mass2 of the compact binaries",
            sample_zs="samples source redshifts",
            sample_geocent_time="samples geocent_time",
            sample_ra="samples right ascension of sky position",
            sample_dec="samples declination of sky position",
            sample_phase="samples coalescence phase",
            sample_psi="samples polarization angle",
            sample_theta_jn="samples inclination angle",
        )

        # spin
        if not self.spin_zero:
            source_priors_["a_1"] = "sampler_uniform"
            event_prior_params_["a_1"] = dict(min_=-a_max, max_=a_max)
            sampler_names_[
                "sample_a_1"
            ] = "samples spin magnitude of the compact binaries (body1)"
            source_priors_["a_2"] = "sampler_uniform"
            event_prior_params_["a_2"] = dict(min_=-a_max, max_=a_max)
            sampler_names_[
                "sample_a_2"
            ] = "samples spin magnitude of the compact binaries (body2)"

            if self.spin_precession:
                source_priors_["a_1"] = "sampler_uniform"
                event_prior_params_["a_1"] = dict(min_=0.0, max_=a_max)
                source_priors_["a_2"] = "sampler_uniform"
                event_prior_params_["a_2"] = dict(min_=0.0, max_=a_max)
                source_priors_["tilt_1"] = "sampler_sine"
                event_prior_params_["tilt_1"] = None
                sampler_names_[
                    "sample_tilt_1"
                ] = "samples tilt angle of the compact binaries (body1)"

                source_priors_["tilt_2"] = "sampler_sine"
                event_prior_params_["tilt_2"] = None
                sampler_names_[
                    "sample_tilt_2"
                ] = "samples tilt angle of the compact binaries (body2)"

                source_priors_["phi_12"] = "sampler_uniform"
                event_prior_params_["phi_12"] = dict(min_=0, max_=2 * np.pi)
                sampler_names_[
                    "sample_phi_12"
                ] = "samples azimuthal angle between the two spins"
                source_priors_["phi_jl"] = "sampler_uniform"
                event_prior_params_["phi_jl"] = dict(min_=0, max_=2 * np.pi)
                sampler_names_[
                    "sample_phi_jl"
                ] = "samples azimuthal angle between the total angular momentum and the orbital angular momentum"

        # update the priors if input is given
        if source_priors:
            source_priors_.update(source_priors)
        if event_prior_params:
            event_prior_params_.update(event_prior_params)

        return (source_priors_, event_prior_params_, sampler_names_)

    @property
    def sample_zs(self):
        """
        Function to sample redshifts with the initialized prior.

        Parameters
        ----------
        size : `int`
            Number of samples to draw

        Returns
        ----------
        zs : `numpy.ndarray` (1D array of floats)
            Array of redshifts
        """

        return self._sample_zs

    @sample_zs.setter
    def sample_zs(self, prior):
        try:
            try:
                self._sample_zs = getattr(self, prior)
            except:
                args = self.gw_param_samplers_params["zs"]
                self._sample_zs = getattr(self, prior)(
                    size=None, get_attribute=True, param=args
                )
        except:
            self._sample_zs = prior

    @property
    def sample_source_frame_masses(self):
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

        return self._sample_source_frame_masses

    @sample_source_frame_masses.setter
    def sample_source_frame_masses(self, prior):
        args = self.gw_param_samplers_params["source_frame_masses"]
    
        # Check if 'prior' is a callable attribute of self
        args = self.gw_param_samplers_params["source_frame_masses"]
        try:
            # If prior is a method of the class
            self._sample_source_frame_masses = getattr(self, prior)(
                size=None, get_attribute=True, param=args,
            )
        except:
            # If prior is a standalone function
            try:
                self._sample_source_frame_masses = lambda size: prior(size, param=args)
            except:
                raise ValueError("given source_frame_masses function should follow the signature of the default ones")

    @property
    def sample_geocent_time(self):
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

        return self._sample_geocent_time

    @sample_geocent_time.setter
    def sample_geocent_time(self, prior):
        try:
            args = self.gw_param_samplers_params["geocent_time"]
            # follwing should return a sampler function with only one argument (size)
            self._sample_geocent_time = getattr(self, prior)(
                size=None, get_attribute=True, param=args,
            )
        except:
            self._sample_geocent_time = prior

    @property
    def sample_ra(self):
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

        return self._sample_ra

    @sample_ra.setter
    def sample_ra(self, prior):
        try:
            args = self.gw_param_samplers_params["ra"]
            self._sample_ra = getattr(self, prior)(
                size=None, get_attribute=True, param=args
            )
        except:
            self._sample_ra = prior

    @property
    def sample_dec(self):
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

        return self._sample_dec

    @sample_dec.setter
    def sample_dec(self, prior):
        try:
            args = self.gw_param_samplers_params["dec"]
            self._sample_dec = getattr(self, prior)(
                size=None, get_attribute=True, param=args
            )
        except:
            self._sample_dec = prior

    @property
    def sample_phase(self):
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

        return self._sample_phase

    @sample_phase.setter
    def sample_phase(self, prior):
        try:
            args = self.gw_param_samplers_params["phase"]
            self._sample_phase = getattr(self, prior)(
                size=None, get_attribute=True, param=args
            )
        except:
            self._sample_phase = prior

    @property
    def sample_psi(self):
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

        return self._sample_psi

    @sample_psi.setter
    def sample_psi(self, prior):
        try:
            args = self.gw_param_samplers_params["psi"]
            self._sample_psi = getattr(self, prior)(
                size=None, get_attribute=True, param=args
            )
        except:
            self._sample_psi = prior

    @property
    def sample_theta_jn(self):
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
        return self._sample_theta_jn

    @sample_theta_jn.setter
    def sample_theta_jn(self, prior):
        try:
            args = self.gw_param_samplers_params["theta_jn"]
            self._sample_theta_jn = getattr(self, prior)(
                size=None, get_attribute=True, param=args
            )
        except:
            self._sample_theta_jn = prior

    @property
    def sample_a_1(self):
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
        return self._sample_a_1

    @sample_a_1.setter
    def sample_a_1(self, prior):
        try:
            args = self.gw_param_samplers_params["a_1"]
            self._sample_a_1 = getattr(self, prior)(
                size=None, get_attribute=True, param=args
            )
        except:
            self._sample_a_1 = prior

    @property
    def sample_a_2(self):
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
        return self._sample_a_2

    @sample_a_2.setter
    def sample_a_2(self, prior):
        try:
            args = self.gw_param_samplers_params["a_2"]
            self._sample_a_2 = getattr(self, prior)(
                size=None, get_attribute=True, param=args
            )
        except:
            self._sample_a_2 = prior

    @property
    def sample_tilt_1(self):
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
        return self._sample_tilt_1

    @sample_tilt_1.setter
    def sample_tilt_1(self, prior):
        try:
            args = self.gw_param_samplers_params["tilt_1"]
            self._sample_tilt_1 = getattr(self, prior)(
                size=None, get_attribute=True, param=args
            )
        except:
            self._sample_tilt_1 = prior

    @property
    def sample_tilt_2(self):
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

        return self._sample_tilt_2

    @sample_tilt_2.setter
    def sample_tilt_2(self, prior):
        try:
            args = self.gw_param_samplers_params["tilt_2"]
            self._sample_tilt_2 = getattr(self, prior)(
                size=None, get_attribute=True, param=args
            )
        except:
            self._sample_tilt_2 = prior

    @property
    def sample_phi_12(self):
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

        return self._sample_phi_12

    @sample_phi_12.setter
    def sample_phi_12(self, prior):
        try:
            args = self.gw_param_samplers_params["phi_12"]
            self._sample_phi_12 = getattr(self, prior)(
                size=None, get_attribute=True, param=args
            )
        except:
            self._sample_phi_12 = prior

    @property
    def sample_phi_jl(self):
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
        return self._sample_phi_jl

    @sample_phi_jl.setter
    def sample_phi_jl(self, prior):
        try:
            args = self.gw_param_samplers_params["phi_jl"]
            self._sample_phi_jl = getattr(self, prior)(
                size=None, get_attribute=True, param=args
            )
        except:
            self._sample_phi_jl = prior
