# -*- coding: utf-8 -*-
"""
Module for compact binary coalescence (CBC) source redshift distribution.

This module provides the :class:`~CBCSourceRedshiftDistribution` class for
generating redshift distributions of compact binary sources (BBH, BNS, NSBH)
based on various merger rate density models. It supports multiple astrophysical
merger rate density prescriptions including PopI/II, PopIII, and primordial
black hole models.

Inheritance hierarchy:

- :class:`~ler.gw_source_population.CBCSourceParameterDistribution` \n
  ↳ inherits from this class

Usage:
    Basic workflow example:

    >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
    >>> cbc = CBCSourceRedshiftDistribution(z_min=0.001, z_max=10)
    >>> zs_samples = cbc.source_redshift(size=1000)

Copyright (C) 2024 Hemantakumar Phurailatpam. Distributed under MIT License.
"""

import warnings

warnings.filterwarnings("ignore")

import numpy as np
from multiprocessing import Pool
from tqdm import tqdm
from scipy.interpolate import CubicSpline

# for redshift to luminosity distance conversion
from astropy.cosmology import LambdaCDM

from ..utils import (
    FunctionConditioning,
    interpolator_json_path,
    luminosity_distance,
    differential_comoving_volume,
)


class CBCSourceRedshiftDistribution(object):
    """
    Class for generating compact binary coalescence source redshift distributions.

    This class generates source redshift distributions for compact binary
    coalescence events (BBH, BNS, NSBH) using various astrophysical merger rate
    density models. It provides interpolated functions for efficient sampling of
    source redshifts weighted by the merger rate density in the detector frame.

    Key Features: \n
    - Multiple merger rate density models (PopI/II, PopIII, Primordial) \n
    - Configurable cosmology for distance calculations \n
    - Cached interpolators for computational efficiency \n
    - Support for user-defined merger rate density functions \n

    Parameters
    ----------
    npool : ``int``
        Number of processors to use for multiprocessing. \n
        default: 4
    z_min : ``float``
        Minimum redshift of the source population. \n
        default: 0.001
    z_max : ``float``
        Maximum redshift of the source population. \n
        default: 10.0
    event_type : ``str``
        Type of compact binary event. \n
        Options: \n
        - 'BBH': Binary black hole \n
        - 'BNS': Binary neutron star \n
        - 'NSBH': Neutron star-black hole \n
        default: 'BBH'
    merger_rate_density : ``str`` or ``callable`` or ``None``
        Merger rate density model to use. \n
        Options: \n
        - 'merger_rate_density_bbh_oguri2018': PopI/II BBH (Oguri 2018) \n
        - 'merger_rate_density_bbh_popIII_ken2022': PopIII BBH (Ng 2022) \n
        - 'merger_rate_density_bbh_primordial_ken2022': Primordial BBH (Ng 2022) \n
        - callable: User-defined function f(z) -> rate density \n
        default: None (uses 'merger_rate_density_bbh_oguri2018')
    merger_rate_density_param : ``dict`` or ``None``
        Parameters for the merger rate density function. \n
        default: None (uses dict(R0=19 * 1e-9, b2=1.6, b3=2.1, b4=30))
    cosmology : ``astropy.cosmology`` or ``None``
        Cosmology for distance calculations. \n
        default: None (uses LambdaCDM(H0=70, Om0=0.3, Ode0=0.7))
    directory : ``str``
        Directory to store interpolator JSON files. \n
        default: './interpolator_json'
    create_new_interpolator : ``dict`` or ``bool``
        Control interpolator creation. \n
        If ``bool``: Apply to all interpolators. \n
        If ``dict``: Per-quantity settings with keys 'create_new' and 'resolution'. \n
        default: False

    Examples
    --------
    Basic usage:

    >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
    >>> cbc = CBCSourceRedshiftDistribution(z_min=0.001, z_max=10)
    >>> zs_samples = cbc.source_redshift(size=1000)
    >>> rate = cbc.merger_rate_density(zs=0.5)

    Instance Methods
    ----------
    CBCSourceRedshiftDistribution has the following methods: \n
    +-----------------------------------------------------+----------------------------------------------------+
    | Method                                              | Description                                        |
    +=====================================================+====================================================+
    | :meth:`~merger_rate_density_detector_frame`         | Compute merger rate density in detector frame      |
    +-----------------------------------------------------+----------------------------------------------------+
    | :meth:`~merger_rate_density_bbh_oguri2018`  | PopI/II merger rate density (Oguri 2018)           |
    +-----------------------------------------------------+----------------------------------------------------+
    | :meth:`~sfr_madau_dickinson2014`                    | Star formation rate (Madau & Dickinson 2014)       |
    +-----------------------------------------------------+----------------------------------------------------+
    | :meth:`~sfr_with_time_delay`                                | SFR with time delay convolution                    |
    +-----------------------------------------------------+----------------------------------------------------+
    | :meth:`~merger_rate_density_bbh_popIII_ken2022`     | PopIII merger rate density (Ng 2022)               |
    +-----------------------------------------------------+----------------------------------------------------+
    | :meth:`~merger_rate_density_bbh_primordial_ken2022` | Primordial BBH merger rate density (Ng 2022)       |
    +-----------------------------------------------------+----------------------------------------------------+

    Instance Attributes
    ----------
    CBCSourceRedshiftDistribution has the following attributes: \n
    +------------------------------------------------+---------------------------+-------+----------------------------------------------+
    | Attribute                                      | Type                      | Unit  | Description                                  |
    +================================================+===========================+=======+==============================================+
    | :attr:`~z_min`                                 | ``float``                 |       | Minimum source redshift                      |
    +------------------------------------------------+---------------------------+-------+----------------------------------------------+
    | :attr:`~z_max`                                 | ``float``                 |       | Maximum source redshift                      |
    +------------------------------------------------+---------------------------+-------+----------------------------------------------+
    | :attr:`~event_type`                            | ``str``                   |       | Type of CBC event (BBH/BNS/NSBH)             |
    +------------------------------------------------+---------------------------+-------+----------------------------------------------+
    | :attr:`~cosmo`                                 | ``astropy.cosmology``     |       | Cosmology for calculations                   |
    +------------------------------------------------+---------------------------+-------+----------------------------------------------+
    | :attr:`~directory`                             | ``str``                   |       | Path for storing interpolators               |
    +------------------------------------------------+---------------------------+-------+----------------------------------------------+
    | :attr:`~merger_rate_density_param`             | ``dict``                  |       | Merger rate density parameters               |
    +------------------------------------------------+---------------------------+-------+----------------------------------------------+
    | :attr:`~normalization_pdf_z`                   | ``float``                 |       | Normalization constant for p(z)              |
    +------------------------------------------------+---------------------------+-------+----------------------------------------------+
    | :attr:`~merger_rate_density`                   | ``callable``              |       | Merger rate density function R(z)            |
    +------------------------------------------------+---------------------------+-------+----------------------------------------------+
    | :attr:`~merger_rate_density_model_list`        | ``dict``                  |       | Available merger rate density models         |
    +------------------------------------------------+---------------------------+-------+----------------------------------------------+
    | :attr:`~source_redshift`                       | ``FunctionConditioning``  |       | Source redshift sampler                      |
    +------------------------------------------------+---------------------------+-------+----------------------------------------------+
    | :attr:`~luminosity_distance`                   | ``FunctionConditioning``  |       | Luminosity distance interpolator             |
    +------------------------------------------------+---------------------------+-------+----------------------------------------------+
    | :attr:`~differential_comoving_volume`          | ``FunctionConditioning``  |       | dVc/dz interpolator                          |
    +------------------------------------------------+---------------------------+-------+----------------------------------------------+
    """

    def __init__(
        self,
        npool=4,
        z_min=0.001,
        z_max=10.0,
        event_type="BBH",
        merger_rate_density=None,
        merger_rate_density_param=None,
        cosmology=None,
        directory="./interpolator_json",
        create_new_interpolator=False,
    ):
        print("\nInitializing CBCSourceRedshiftDistribution class...\n")
        # set attributes
        self.npool = npool
        self.z_min = z_min
        self.z_max = z_max
        self.directory = directory
        self.event_type = event_type
        # if None is passed, use the default cosmology
        self.cosmo = cosmology if cosmology else LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

        # setting up the interpolator creation parameters
        self.create_new_interpolator = self._setup_decision_dictionary(
            create_new_interpolator, merger_rate_density
        )

        # Initialize cosmological functions
        self.luminosity_distance = luminosity_distance(
            z_min=self.z_min,
            z_max=self.z_max,
            cosmo=self.cosmo,
            directory=self.directory,
            create_new=self.create_new_interpolator["luminosity_distance"][
                "create_new"
            ],
            resolution=self.create_new_interpolator["luminosity_distance"][
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

        # function initialization
        merger_rate_density, self.merger_rate_density_param = (
            self._merger_rate_density_priors_categorization(
                event_type, merger_rate_density, merger_rate_density_param
            )
        )
        self.merger_rate_density = (
            merger_rate_density  # this is an initialization, not a variable assignment
        )

        # source redshift distribution initialization
        self.merger_rate_density_detector_frame = (
            self.merger_rate_density_detector_frame(zs=None, get_attribute=True)
        )
        self.source_redshift = FunctionConditioning(
            function=self.merger_rate_density_detector_frame.z_array,
            x_array=self.merger_rate_density_detector_frame.x_array,
            identifier_dict=self.merger_rate_density_detector_frame.info,
            directory=self.directory,
            sub_directory="source_redshift",
            name="source_redshift",
            create_new=self.create_new_interpolator["merger_rate_density"][
                "create_new"
            ],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
        )

        # Normalization of the pdf p(z)
        self.normalization_pdf_z = float(self.source_redshift.pdf_norm_const)

    def _setup_decision_dictionary(self, create_new_interpolator, merger_rate_density):
        """
        Helper to set up a decision dictionary for interpolator creation.

        Parameters
        ----------
        create_new_interpolator : ``dict`` or ``bool``
            If ``dict``: Per-quantity settings with 'create_new' and 'resolution'.
            If ``bool``: Apply to all interpolators.
        merger_rate_density : ``str`` or ``None``
            Merger rate density model name for resolution adjustment.

        Returns
        -------
        create_new_interpolator_ : ``dict``
            Dictionary with interpolator creation settings.
        """
        create_new_interpolator_ = dict(
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

        if isinstance(merger_rate_density, str):
            if merger_rate_density == "sfr_with_time_delay":
                create_new_interpolator_["merger_rate_density"]["resolution"] = 48

        return create_new_interpolator_

    def _merger_rate_density_priors_categorization(
        self, event_type, merger_rate_density, merger_rate_density_param
    ):
        """
        Helper to categorize merger rate density and set default parameters.

        Parameters
        ----------
        event_type : ``str``
            Type of CBC event ('BBH', 'BNS', 'NSBH').
        merger_rate_density : ``str`` or ``callable`` or ``None``
            Merger rate density model name or user function.
        merger_rate_density_param : ``dict`` or ``None``
            User-provided parameters to override defaults.

        Returns
        -------
        merger_rate_density_ : ``str`` or ``callable``
            Validated merger rate density model or function.
        merger_rate_density_param_ : ``dict``
            Parameters for the merger rate density function.
        """

        if event_type == "BBH":
            merger_rate_density_ = "merger_rate_density_madau_dickinson_belczynski_ng"
            merger_rate_density_param_ = dict(
                R0=19 * 1e-9, alpha_F=2.57, beta_F=5.83, c_F=3.36
            )
        elif event_type == "BNS":
            merger_rate_density_ = "merger_rate_density_madau_dickinson2014"
            merger_rate_density_param_ = dict(
                R0=89 * 1e-9, a=0.015, b=2.7, c=2.9, d=5.6
            )
        elif event_type == "NSBH":
            merger_rate_density_ = "merger_rate_density_madau_dickinson2014"
            merger_rate_density_param_ = dict(
                R0=23 * 1e-9, a=0.015, b=2.7, c=2.9, d=5.6
            )

        if merger_rate_density:
            merger_rate_density_ = merger_rate_density

            if isinstance(merger_rate_density, str):
                if merger_rate_density in self.merger_rate_density_model_list:
                    merger_rate_density_param_ = self.merger_rate_density_model_list[
                        merger_rate_density
                    ]
                else:
                    raise ValueError(
                        f"'merger rate density' sampler '{merger_rate_density}' not available.\n Available 'merger rate density' samplers and its parameters are: {self.merger_rate_density_model_list}"
                    )
            elif callable(merger_rate_density):
                print("using user provided custom merger rate density function")
                merger_rate_density_param_ = {}
            else:
                raise ValueError("merger_rate_density must be a string or callable")

        if merger_rate_density_param:
            merger_rate_density_param_.update(merger_rate_density_param)

        return merger_rate_density_, merger_rate_density_param_

    def merger_rate_density_detector_frame(self, zs, get_attribute=False, **kwargs):
        """
        Compute the merger rate density in the detector frame.

        The output is unnormalized and accounts for the (1+z) time dilation
        factor and differential comoving volume.

        Parameters
        ----------
        zs : ```numpy.ndarray``
            Source redshift(s) at which to evaluate the rate density.
        get_attribute : ``bool``
            If True, return the FunctionConditioning object.
            If False, return evaluated values.
            default: False
        **kwargs : ``dict``
            Additional keyword arguments (unused).

        Returns
        -------
        rate_density : ``numpy.ndarray`` or ``FunctionConditioning``
            Merger rate density in detector frame (units: Mpc^-3 yr^-1 sr^-1).

        Examples
        --------
        >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
        >>> cbc = CBCSourceRedshiftDistribution()
        >>> rate_density = cbc.merger_rate_density_detector_frame(zs=0.1)
        """

        identifier_dict = self.merger_rate_density_param.copy()
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["event_type"] = self.event_type
        identifier_dict["name"] = "merger_rate_density_detector_frame"
        identifier_dict["resolution"] = self.create_new_interpolator[
            "merger_rate_density"
        ]["resolution"]

        zs_array = np.linspace(
            identifier_dict["z_min"],
            identifier_dict["z_max"],
            identifier_dict["resolution"],
        )
        Pzs = (
            lambda z: self.merger_rate_density(z)
            / (1 + z)
            * self.differential_comoving_volume(z)
        )

        Pzs_object = FunctionConditioning(
            function=Pzs,
            x_array=zs_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="merger_rate_density",
            name=identifier_dict["name"],
            create_new=self.create_new_interpolator["merger_rate_density"][
                "create_new"
            ],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="function",
        )

        return Pzs_object if get_attribute else Pzs_object(zs)

    def merger_rate_density_bbh_oguri2018(self, zs, get_attribute=False, **kwargs):
        """
        Compute PopI/II BBH merger rate density (Oguri et al. 2018).

        Returns the source-frame merger rate density following the
        Oguri et al. (2018) prescription for PopI/II stellar populations.

        Parameters
        ----------
        zs : ```numpy.ndarray``
            Source redshift(s) at which to evaluate.
        get_attribute : ``bool``
            If True, return the FunctionConditioning object.
            default: False
        **kwargs : ``dict``
            Override default fitting parameters:
            R0=19e-9, b2=1.6, b3=2.1, b4=30.

        Returns
        -------
        rate_density : ```numpy.ndarray`` or ``FunctionConditioning``
            Merger rate density in source frame (units: Mpc^-3 yr^-1).

        Examples
        --------
        >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
        >>> cbc = CBCSourceRedshiftDistribution(merger_rate_density="merger_rate_density_bbh_oguri2018")
        >>> rate = cbc.merger_rate_density(zs=0.5)
        """

        from .prior_functions import merger_rate_density_bbh_oguri2018_function

        identifier_dict = {}
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["event_type"] = self.event_type
        identifier_dict["name"] = "merger_rate_density_bbh_oguri2018"
        identifier_dict["resolution"] = self.create_new_interpolator[
            "merger_rate_density"
        ]["resolution"]
        param_dict = self.merger_rate_density_model_list[
            "merger_rate_density_bbh_oguri2018"
        ].copy()
        param_dict.update(kwargs)
        identifier_dict.update(param_dict)

        zs_array = np.linspace(
            identifier_dict["z_min"],
            identifier_dict["z_max"],
            identifier_dict["resolution"],
        )
        Rzs = lambda zs: merger_rate_density_bbh_oguri2018_function(
            zs=zs,
            R0=identifier_dict["R0"],
            b2=identifier_dict["b2"],
            b3=identifier_dict["b3"],
            b4=identifier_dict["b4"],
        )

        Rzs_object = FunctionConditioning(
            function=Rzs,
            x_array=zs_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="merger_rate_density",
            name=identifier_dict["name"],
            create_new=self.create_new_interpolator["merger_rate_density"][
                "create_new"
            ],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="function",
        )

        return Rzs_object if get_attribute else Rzs_object(zs)

    def sfr_with_time_delay(self, zs, get_attribute=False, **kwargs):
        """
        Compute merger rate density with time delay convolution.

        Convolves the star formation rate with a time delay distribution
        to compute the merger rate density. Uses multiprocessing for
        numerical integration (Borhanian & Sathyaprakash 2024).

        Parameters
        ----------
        zs : ```numpy.ndarray``
            Source redshift(s) at which to evaluate.
        get_attribute : ``bool``
            If True, return the FunctionConditioning object.
            default: False
        **kwargs : ``dict``
            Override default SFR and time delay parameters.

        Returns
        -------
        rate_density : ```numpy.ndarray`` or ``FunctionConditioning``
            Merger rate density (units: Mpc^-3 yr^-1).
        """

        identifier_dict = {}
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["event_type"] = self.event_type
        identifier_dict["name"] = "sfr_with_time_delay"
        identifier_dict["resolution"] = self.create_new_interpolator[
            "merger_rate_density"
        ]["resolution"]
        param_dict = self.merger_rate_density_model_list["sfr_with_time_delay"].copy()
        param_dict.update(kwargs)
        identifier_dict.update(param_dict)

        print("Numerically solving the merger_rate_density with time delay")
        zs_resolution = identifier_dict["resolution"]
        zs_array = (
            np.geomspace(self.z_min + 0.001, self.z_max, zs_resolution)
            if self.z_min == 0
            else np.geomspace(self.z_min, self.z_max, zs_resolution)
        )

        _, it_exist = interpolator_json_path(
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="merger_rate_density",
            interpolator_name=identifier_dict["name"],
        )

        create_new = self.create_new_interpolator["merger_rate_density"]["create_new"]
        if not it_exist or create_new:
            rate_density = self._helper_rate_density_multiprocessing(
                zs_array, identifier_dict
            )
        else:
            rate_density = None

        rate_density_object = FunctionConditioning(
            function=rate_density,
            x_array=zs_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="merger_rate_density",
            name=identifier_dict["name"],
            create_new=create_new,
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="function",
        )

        return rate_density_object if get_attribute else rate_density_object(zs)

    def _helper_rate_density_multiprocessing(self, zs_array, identifier_dict):
        """
        Helper to compute merger rate density with time delay using multiprocessing.

        Parameters
        ----------
        zs_array : ``numpy.ndarray``
            1D array of source redshifts.
        identifier_dict : ``dict``
            Parameters including cosmology (H0, Om0, Ode0),
            time delay bounds (td_min, td_max), and SFR parameters (a, b, c, d).

        Returns
        -------
        rate_density_array : ``numpy.ndarray``
            Computed merger rate density normalized to R0.
        """

        from .sfr_with_time_delay import sfr_with_time_delay_function

        size = len(zs_array)
        input_args = np.array(
            [
                zs_array,  # source redshifts
                np.arange(size),  # index
                identifier_dict["td_min"] * np.ones(size),  # time delay minimum
                identifier_dict["td_max"] * np.ones(size),  # time delay maximum
                identifier_dict["cosmology"].H0.value * np.ones(size),  # H0
                identifier_dict["cosmology"].Om0 * np.ones(size),  # Omega_M
                identifier_dict["cosmology"].Ode0 * np.ones(size),  # Omega_Lambda
                # SFR params
                identifier_dict["a"] * np.ones(size),
                identifier_dict["b"] * np.ones(size),
                identifier_dict["c"] * np.ones(size),
                identifier_dict["d"] * np.ones(size),
            ]
        ).T

        print(
            "Computing merger rate density distribution (with time delay to SFR) using multiprocessing..."
        )
        # with tqdm
        rate_density_array = np.zeros(size)
        with Pool(processes=self.npool) as pool:
            for result in tqdm(
                pool.imap_unordered(sfr_with_time_delay_function, input_args),
                total=size,
                ncols=100,
                disable=False,
            ):
                # print(result)
                (
                    iter_i,
                    density_,
                ) = result

                rate_density_array[iter_i] = density_

        rm_spline = CubicSpline(zs_array, rate_density_array, extrapolate=True)
        # divide
        rate_density_array = rate_density_array / rm_spline(0.0) * identifier_dict["R0"]

        return rate_density_array

    def merger_rate_density_madau_dickinson2014(
        self, zs, get_attribute=False, **kwargs
    ):
        """
        Compute star formation rate following Madau & Dickinson (2014).

        Returns the cosmic star formation rate density as given in
        Equation 15 of Madau & Dickinson (2014).

        Parameters
        ----------
        zs : ```numpy.ndarray``
            Source redshift(s) at which to evaluate.
        get_attribute : ``bool``
            If True, return the FunctionConditioning object.
            default: False
        **kwargs : ``dict``
            Override default fitting parameters: R0=19 * 1e-9, a=0.015, b=2.7, c=2.9, d=5.6.

        Returns
        -------
        rate_density : ```numpy.ndarray`` or ``FunctionConditioning``
            Star formation rate density (units: M_sun yr^-1 Mpc^-3).

        Examples
        --------
        >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
        >>> cbc = CBCSourceRedshiftDistribution(merger_rate_density="merger_rate_density_madau_dickinson2014")
        >>> sfr = cbc.merger_rate_density(zs=2.0)
        """

        from .prior_functions import merger_rate_density_madau_dickinson2014_function

        identifier_dict = {}
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["event_type"] = self.event_type
        identifier_dict["name"] = "merger_rate_density_madau_dickinson2014"
        identifier_dict["resolution"] = self.create_new_interpolator[
            "merger_rate_density"
        ]["resolution"]
        param_dict = self.merger_rate_density_model_list[
            "merger_rate_density_madau_dickinson2014"
        ].copy()
        param_dict.update(kwargs)
        identifier_dict.update(param_dict)

        zs_array = np.linspace(
            identifier_dict["z_min"],
            identifier_dict["z_max"],
            identifier_dict["resolution"],
        )

        Rzs = lambda zs: merger_rate_density_madau_dickinson2014_function(
            zs=zs,
            R0=identifier_dict["R0"],
            a=identifier_dict["a"],
            b=identifier_dict["b"],
            c=identifier_dict["c"],
            d=identifier_dict["d"],
        )

        Rzs_object = FunctionConditioning(
            function=Rzs,
            x_array=zs_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="merger_rate_density",
            name=identifier_dict["name"],
            create_new=self.create_new_interpolator["merger_rate_density"][
                "create_new"
            ],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="function",
        )

        return Rzs_object if get_attribute else Rzs_object(zs)

    def merger_rate_density_madau_dickinson_belczynski_ng(
        self, zs, get_attribute=False, **kwargs
    ):
        """
        Compute BBH merger rate density following Ng et al. (2021).

        This model uses a Madau-Dickinson-like functional form to fit the
        merger rate density of field BHs, accounting for time delays and
        metallicity effects.

        density(zs) ∝ (1 + zs) ** alpha_F / (1 + ((1 + zs) / c_F) ** beta_F)

        Parameters
        ----------
        zs : ```numpy.ndarray``
            Source redshift(s) at which to evaluate.
        get_attribute : ``bool``
            If True, return the FunctionConditioning object.
            default: False
        **kwargs : ``dict``
            Override default fitting parameters: R0, alpha_F, beta_F, c_F.

        Returns
        -------
        rate_density : ```numpy.ndarray`` or ``FunctionConditioning``
            Star formation rate density (units: M_sun yr^-1 Mpc^-3).

        Examples
        --------
        >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
        >>> cbc = CBCSourceRedshiftDistribution(merger_rate_density="merger_rate_density_madau_dickinson_belczynski_ng")
        >>> sfr = cbc.merger_rate_density(zs=2.0)
        """

        from .prior_functions import (
            merger_rate_density_madau_dickinson_belczynski_ng_function,
        )

        identifier_dict = {}
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["event_type"] = self.event_type
        identifier_dict["name"] = "merger_rate_density_madau_dickinson_belczynski_ng"
        identifier_dict["resolution"] = self.create_new_interpolator[
            "merger_rate_density"
        ]["resolution"]
        param_dict = self.merger_rate_density_model_list[
            "merger_rate_density_madau_dickinson_belczynski_ng"
        ].copy()
        param_dict.update(kwargs)
        identifier_dict.update(param_dict)

        zs_array = np.linspace(
            identifier_dict["z_min"],
            identifier_dict["z_max"],
            identifier_dict["resolution"],
        )

        Rzs = merger_rate_density_madau_dickinson_belczynski_ng_function(
            zs=zs_array,
            R0=identifier_dict["R0"],
            alpha_F=identifier_dict["alpha_F"],
            beta_F=identifier_dict["beta_F"],
            c_F=identifier_dict["c_F"],
        )

        Rzs_object = FunctionConditioning(
            function=Rzs,
            x_array=zs_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="merger_rate_density",
            name=identifier_dict["name"],
            create_new=self.create_new_interpolator["merger_rate_density"][
                "create_new"
            ],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="function",
        )

        return Rzs_object if get_attribute else Rzs_object(zs)

    def merger_rate_density_bbh_popIII_ken2022(self, zs, get_attribute=False, **kwargs):
        """
        Compute PopIII BBH merger rate density (Ng et al. 2022).

        Returns the merger rate density for Population III binary black
        holes following the Ng et al. (2022) prescription.

        Parameters
        ----------
        zs : ```numpy.ndarray``
            Source redshift(s) at which to evaluate.
        get_attribute : ``bool``
            If True, return the FunctionConditioning object.
            default: False
        **kwargs : ``dict``
            Override default fitting parameters:
            R0=19.2e-9, aIII=0.66, bIII=0.3, zIII=11.6.

        Returns
        -------
        rate_density : ```numpy.ndarray`` or ``FunctionConditioning``
            Merger rate density (units: Mpc^-3 yr^-1).

        Examples
        --------
        >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
        >>> cbc = CBCSourceRedshiftDistribution(
        ...     z_min=5, z_max=40,
        ...     merger_rate_density="merger_rate_density_bbh_popIII_ken2022"
        ... )
        >>> rate = cbc.merger_rate_density(zs=10)
        """

        from .prior_functions import merger_rate_density_bbh_popIII_ken2022_function

        identifier_dict = {}
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["event_type"] = self.event_type
        identifier_dict["name"] = "merger_rate_density_bbh_popIII_ken2022"
        identifier_dict["resolution"] = self.create_new_interpolator[
            "merger_rate_density"
        ]["resolution"]
        param_dict = self.merger_rate_density_model_list[
            "merger_rate_density_bbh_popIII_ken2022"
        ].copy()
        param_dict.update(kwargs)
        identifier_dict.update(param_dict)

        zs_array = np.linspace(
            identifier_dict["z_min"],
            identifier_dict["z_max"],
            identifier_dict["resolution"],
        )
        Rzs = lambda zs: merger_rate_density_bbh_popIII_ken2022_function(
            zs=zs,
            R0=identifier_dict["R0"],
            aIII=identifier_dict["aIII"],
            bIII=identifier_dict["bIII"],
            zIII=identifier_dict["zIII"],
        )

        Rzs_object = FunctionConditioning(
            function=Rzs,
            x_array=zs_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="merger_rate_density",
            name=identifier_dict["name"],
            create_new=self.create_new_interpolator["merger_rate_density"][
                "create_new"
            ],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="function",
        )

        return Rzs_object if get_attribute else Rzs_object(zs)

    def merger_rate_density_bbh_primordial_ken2022(
        self, zs, get_attribute=False, **kwargs
    ):
        """
        Compute primordial BBH merger rate density (Ng et al. 2022).

        Returns the merger rate density for primordial binary black holes
        following the Ng et al. (2022) prescription.

        Parameters
        ----------
        zs : ```numpy.ndarray``
            Source redshift(s) at which to evaluate.
        get_attribute : ``bool``
            If True, return the FunctionConditioning object.
            default: False
        **kwargs : ``dict``
            Override default fitting parameters:
            R0=0.044e-9, t0=13.786885302009708.

        Returns
        -------
        rate_density : ```numpy.ndarray`` or ``FunctionConditioning``
            Merger rate density (units: Mpc^-3 yr^-1).

        Examples
        --------
        >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
        >>> cbc = CBCSourceRedshiftDistribution(
        ...     z_min=5, z_max=40,
        ...     merger_rate_density="merger_rate_density_bbh_primordial_ken2022"
        ... )
        >>> rate = cbc.merger_rate_density(zs=10)
        """

        from .prior_functions import merger_rate_density_bbh_primordial_ken2022_function

        identifier_dict = {}
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["event_type"] = self.event_type
        identifier_dict["name"] = "merger_rate_density_bbh_primordial_ken2022"
        identifier_dict["resolution"] = self.create_new_interpolator[
            "merger_rate_density"
        ]["resolution"]
        param_dict = self.merger_rate_density_model_list[
            "merger_rate_density_bbh_primordial_ken2022"
        ].copy()
        param_dict.update(kwargs)
        identifier_dict.update(param_dict)

        zs_array = np.linspace(
            identifier_dict["z_min"],
            identifier_dict["z_max"],
            identifier_dict["resolution"],
        )
        Rzs = lambda zs: merger_rate_density_bbh_primordial_ken2022_function(
            zs=zs, R0=identifier_dict["R0"], t0=identifier_dict["t0"]
        )

        Rzs_object = FunctionConditioning(
            function=Rzs,
            x_array=zs_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="merger_rate_density",
            name=identifier_dict["name"],
            create_new=self.create_new_interpolator["merger_rate_density"][
                "create_new"
            ],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="function",
        )

        return Rzs_object if get_attribute else Rzs_object(zs)

    @property
    def merger_rate_density(self):
        """
        Source-frame merger rate density object. \n

        Returns a ``FunctionConditioning`` object with methods: \n
        - ``function(zs)``: Get merger rate density in source frame \n
        - ``rvs(size)``: Sample source redshifts in source frame \n
        - ``pdf(zs)``: Get probability density \n

        Returns
        -------
        merger_rate_density : ``FunctionConditioning``
            Callable that accepts redshift(s) and returns merger rate density
            in source frame (units: Mpc^-3 yr^-1).
        """
        return self._merger_rate_density

    @merger_rate_density.setter
    def merger_rate_density(self, function):

        if function in self.merger_rate_density_model_list:
            print(f"using ler available merger rate density model: {function}")
            args = self.merger_rate_density_param
            if args is None:
                self._merger_rate_density = getattr(self, function)(
                    zs=None, get_attribute=True
                )
            else:
                self._merger_rate_density = getattr(self, function)(
                    zs=None, get_attribute=True, **args
                )

        elif isinstance(function, FunctionConditioning):
            print("using user provided FunctionConditioning object")
            self._merger_rate_density = function

        elif callable(function):
            print("using user provided custom merger rate density function")
            identifier_dict = {}
            identifier_dict["z_min"] = self.z_min
            identifier_dict["z_max"] = self.z_max
            identifier_dict["name"] = "merger_rate_density_custom"
            identifier_dict["resolution"] = self.create_new_interpolator[
                "merger_rate_density"
            ]["resolution"]
            zs_array = np.linspace(
                identifier_dict["z_min"],
                identifier_dict["z_max"],
                identifier_dict["resolution"],
            )

            self._merger_rate_density = FunctionConditioning(
                function=function,
                x_array=zs_array,
                identifier_dict=identifier_dict,
                directory=self.directory,
                sub_directory="merger_rate_density",
                name=identifier_dict["name"],
                create_new=self.create_new_interpolator["merger_rate_density"][
                    "create_new"
                ],
                create_function_inverse=False,
                create_function=True,
                create_pdf=True,
                create_rvs=True,
                callback="function",
            )

        else:
            raise ValueError(
                "merger_rate_density must be a function, FunctionConditioning object, or a string from the available 'merger_rate_density_model_list'"
            )

    @property
    def npool(self):
        """
        Number of processors for multiprocessing. \n

        Returns
        -------
        npool : ``int``
            Number of parallel processes to use. \n
            default: 4
        """
        return self._npool

    @npool.setter
    def npool(self, value):
        self._npool = value

    @property
    def z_min(self):
        """
        Minimum source redshift. \n

        Returns
        -------
        z_min : ``float``
            Lower bound of the redshift range. \n
            default: 0.001
        """
        return self._z_min

    @z_min.setter
    def z_min(self, value):
        self._z_min = value

    @property
    def z_max(self):
        """
        Maximum source redshift. \n

        Returns
        -------
        z_max : ``float``
            Upper bound of the redshift range. \n
            default: 10.0
        """
        return self._z_max

    @z_max.setter
    def z_max(self, value):
        self._z_max = value

    @property
    def directory(self):
        """
        Directory path for storing interpolator JSON files. \n

        Returns
        -------
        directory : ``str``
            Path to the interpolator storage directory. \n
            default: './interpolator_json'
        """
        return self._directory

    @directory.setter
    def directory(self, value):
        self._directory = value

    @property
    def event_type(self):
        """
        Type of compact binary coalescence event. \n

        Returns
        -------
        event_type : ``str``
            CBC event type. \n
            Options: \n
            - 'BBH': Binary black hole \n
            - 'BNS': Binary neutron star \n
            - 'NSBH': Neutron star-black hole \n
            default: 'BBH'
        """
        return self._event_type

    @event_type.setter
    def event_type(self, value):
        self._event_type = value

    @property
    def cosmo(self):
        """
        Astropy cosmology object for distance calculations. \n

        Returns
        -------
        cosmo : ``astropy.cosmology``
            Cosmology used for redshift-distance conversions. \n
            default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        """
        return self._cosmo

    @cosmo.setter
    def cosmo(self, value):
        self._cosmo = value

    @property
    def create_new_interpolator(self):
        """
        Dictionary controlling interpolator creation settings. \n

        Returns
        -------
        create_new_interpolator : ``dict``
            Dictionary with that controls the creation of new interpolators.
            Default: {'merger_rate_density': {'create_new': False, 'resolution': 100}, 'luminosity_distance': {'create_new': False, 'resolution': 100}, 'differential_comoving_volume': {'create_new': False, 'resolution': 100}}
        """
        return self._create_new_interpolator

    @create_new_interpolator.setter
    def create_new_interpolator(self, value):
        self._create_new_interpolator = value

    @property
    def merger_rate_density_param(self):
        """
        Parameters for the merger rate density function. \n

        Returns
        -------
        merger_rate_density_param : ``dict``
            Dictionary of parameters for the selected merger rate density model.
        """
        return self._merger_rate_density_param

    @merger_rate_density_param.setter
    def merger_rate_density_param(self, value):
        self._merger_rate_density_param = value

    @property
    def luminosity_distance(self):
        """
        Class object (of FunctionConditioning) for the luminosity distance, with function as callback, which converts redshift to luminosity distance (in Mpc) for the selected cosmology. \n
        The class object contains the following attribute methods: \n
        - `function`: returns the luminosity distance distribution function. \n
        - `function_inverse`: returns the inverse luminosity distance distribution function, which converts luminosity distance (in Mpc) to redshift.

        Returns
        -------
        luminosity_distance : ``numpy.ndarray``
            Array of luminosity distances (in Mpc).
        """
        return self._luminosity_distance

    @luminosity_distance.setter
    def luminosity_distance(self, value):
        self._luminosity_distance = value

    @property
    def differential_comoving_volume(self):
        """
        Class object (of FunctionConditioning) for the differential comoving volume function, with function as callback, which returns dVc/dz (in Mpc^3 sr^-1) for the selected cosmology. \n
        The class object contains the following attribute methods: \n
        - `function`: returns the differential comoving volume distribution function.

        Returns
        -------
        differential_comoving_volume : ``numpy.ndarray``
            Array of differential comoving volumes (in Mpc^3 sr^-1).
        """
        return self._differential_comoving_volume

    @differential_comoving_volume.setter
    def differential_comoving_volume(self, value):
        self._differential_comoving_volume = value

    @property
    def source_redshift(self):
        """
        Class object (of FunctionConditioning) for the source redshift sampler, with rvs/sampler as callback, which samples source redshifts from p(z) ∝ R(z)/(1+z) dVc/dz , where p(z) is the redshift probability distribution, R(z) is the merger rate density, and dVc/dz is the differential comoving volume. \n
        The class object contains the following attribute methods: \n
        - `rvs`: returns random samples from the source redshift distribution. \n
        - `pdf`: returns the source redshift probability density function. \n
        - `function`: returns the source redshift distribution function.

        Returns
        -------
        source_redshift : ``numpy.ndarray``
            Array of source redshifts (detector frame)
        """
        return self._source_redshift

    @source_redshift.setter
    def source_redshift(self, value):
        self._source_redshift = value

    @property
    def normalization_pdf_z(self):
        """
        Normalization constant for the redshift probability distribution. \n

        Returns
        -------
        normalization_pdf_z : ``float``
            Integral of the unnormalized p(z) over [z_min, z_max].
        """
        return self._normalization_pdf_z

    @normalization_pdf_z.setter
    def normalization_pdf_z(self, value):
        self._normalization_pdf_z = value

    @property
    def merger_rate_density_model_list(self):
        """
        Dictionary of available merger rate density models and default parameters. \n

        Returns
        -------
        merger_rate_density_model_list : ``dict``
            Dictionary with model names as keys and parameter dicts as values. \n
            Available models: \n
            - 'merger_rate_density_bbh_oguri2018' \n
            - 'sfr_with_time_delay' \n
            - 'merger_rate_density_bbh_popIII_ken2022' \n
            - 'merger_rate_density_bbh_primordial_ken2022'
        """
        # Return cached version if available
        if (
            hasattr(self, "_merger_rate_density_model_list")
            and self._merger_rate_density_model_list is not None
        ):
            return self._merger_rate_density_model_list

        self._merger_rate_density_model_list = dict(
            merger_rate_density_bbh_oguri2018=dict(R0=19 * 1e-9, b2=1.6, b3=2.1, b4=30),
            merger_rate_density_madau_dickinson2014=dict(
                R0=19 * 1e-9, a=0.015, b=2.7, c=2.9, d=5.6
            ),
            merger_rate_density_madau_dickinson_belczynski_ng=dict(
                R0=19 * 1e-9, alpha_F=2.57, beta_F=5.83, c_F=3.36
            ),
            sfr_with_time_delay=dict(
                R0=19 * 1e-9, a=0.01, b=2.6, c=3.2, d=6.2, td_min=10e-3, td_max=10.0
            ),
            merger_rate_density_bbh_popIII_ken2022=dict(
                R0=19.2 * 1e-9, aIII=0.66, bIII=0.3, zIII=11.6
            ),
            merger_rate_density_bbh_primordial_ken2022=dict(
                R0=0.044 * 1e-9, t0=13.786885302009708
            ),
        )

        return self._merger_rate_density_model_list
