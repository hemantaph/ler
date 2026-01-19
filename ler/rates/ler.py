# -*- coding: utf-8 -*-
"""
Module for calculating detection rates of gravitational wave events.

This module contains the main ``LeR`` class for calculating the rates of
detectable gravitational wave events, both lensed and unlensed. The class
inherits from :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution`
for source and lens parameters sampling, and utilizes image property calculations.

The inheritance hierarchy is as follows:

- :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution` \n
  - :class:`~ler.lens_galaxy_population.OpticalDepth` \n
  - :class:`~ler.image_properties.ImageProperties` \n
  - :class:`~ler.gw_source_population.CBCSourceParameterDistribution` \n
  - :class:`~ler.gw_source_population.CBCSourceRedshiftDistribution` \n
- Uses the ``gwsnr`` package for pdet calculation.

Usage:
    Basic workflow for rate calculation:

    >>> from ler.rates import LeR
    >>> ler = LeR()
    >>> unlensed_params = ler.unlensed_cbc_statistics()
    >>> ler.unlensed_rate()
    >>> lensed_params = ler.lensed_cbc_statistics()
    >>> ler.lensed_rate()
    >>> ler.rate_ratio()

Copyright (C) 2026 Phurailatpam Hemantakumar. Distributed under MIT License.
"""

import os

# os.environ['OMP_NESTED'] = 'FALSE'
import warnings
import pathlib
import zipfile
from importlib_resources import files as resources_files

warnings.filterwarnings("ignore")
import logging

logging.getLogger("numexpr.utils").setLevel(logging.ERROR)
import contextlib
import numpy as np
from astropy.cosmology import LambdaCDM
from ..lens_galaxy_population import LensGalaxyParameterDistribution
from ..utils import (
    load_json,
    append_json,
    get_param_from_json,
    batch_handler,
    remove_file,
)


class LeR(LensGalaxyParameterDistribution):
    """
    Class to sample lensed and unlensed GW events and calculate their detection rates.

    This class provides functionality for sampling gravitational wave source parameters,
    detection probabilities, and computing detection rates for both lensed and unlensed
    compact binary coalescence events.
    Parameters of simulated events are stored in JSON files (not as class attributes)
    to conserve RAM memory.

    Key Features: \n
    - Sampling of unlensed and lensed CBC event parameters \n
    - Detection probability calculation using ``gwsnr`` package or custom functions \n
    - Rate calculation for detectable events \n
    - Batch processing for memory efficiency \n
    - JSON-based parameter storage for reproducibility \n

    Parameters
    ----------
    npool : ``int``
        Number of cores to use for parallel processing. \n
        default: 4
    z_min : ``float``
        Minimum redshift of the source population. \n
        default: 0.0
    z_max : ``float``
        Maximum redshift of the source population. \n
        default: 10.0
    event_type : ``str``
        Type of event to generate. source_priors and source_priors_params will be set accordingly. \n
        Options: \n
        - 'BBH': Binary Black Hole \n
        - 'BNS': Binary Neutron Star \n
        - 'NSBH': Neutron Star-Black Hole \n
        default: 'BBH'
    lens_type : ``str``
        Type of lens model to use. lens_functions, lens_functions_params, lens_param_samplers and lens_param_samplers_params will be set accordingly. \n
        Options: \n
        - 'epl_shear_galaxy': Exponential Power Law Shear Galaxy \n
        - 'sie_galaxy': Singular Isothermal Ellipsoid Galaxy \n
        - 'sis_galaxy': Singular Isothermal Sphere Galaxy \n
        default: 'epl_shear_galaxy'
    cosmology : ``astropy.cosmology``
        Cosmology to use for the calculation. \n
        default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
    pdet_finder : ``function`` or ``None``
        Custom detection probability finder function. \n
        If None, uses gwsnr's pdet calculator. \n
        The function should follow the signature: \n
        ``def pdet_finder(gw_param_dict): return pdet_net_dict`` \n
        where pdet_net_dict.keys = ['pdet_net']. \n
        default: None
    json_file_names : ``dict``
        Names of the JSON files to store the necessary parameters. \n
        default: dict(
            ler_params="ler_params.json",
            unlensed_param="unlensed_param.json",
            unlensed_param_detectable="unlensed_param_detectable.json",
            lensed_param="lensed_param.json",
            lensed_param_detectable="lensed_param_detectable.json"
        )
    interpolator_directory : ``str``
        Directory to store the interpolators. \n
        default: './interpolator_json'
    create_new_interpolator : ``bool`` or ``dict``
        Whether to create new interpolators. Look at :meth:`~ler.ler_rates.LER.create_new_interpolator` for details. \n
        Options: \n
        - True: Create all interpolators anew \n
        - False: Load existing interpolators if available \n
        - dict: Specify which interpolators to create new \n
        default: False
    ler_directory : ``str``
        Directory to store the output parameters. \n
        default: './ler_data'
    verbose : ``bool``
        If True, print all chosen parameters during initialization. \n
        default: True
    **kwargs : ``dict``
        Additional keyword arguments passed to parent classes: \n
        :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution`, \n
        :class:`~ler.gw_source_population.CBCSourceParameterDistribution`, \n
        :class:`~ler.image_properties.ImageProperties`, and \n
        :class:`~gwsnr.GWSNR` (if snr_finder='gwsnr').

    Examples
    --------
    Basic usage:

    >>> from ler import LeR
    >>> ler = LeR()
    >>> unlensed_params = ler.unlensed_cbc_statistics()
    >>> ler.unlensed_rate()
    >>> lensed_params = ler.lensed_cbc_statistics()
    >>> ler.lensed_rate()
    >>> ler.rate_ratio()


    Instance Methods
    ----------
    LeR class has the following methods: \n
    +-----------------------------------------------------+------------------------------------------------+
    | Method                                              | Description                                    |
    +=====================================================+================================================+
    | :meth:`~unlensed_cbc_statistics`                    | Generate unlensed GW source parameters         |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~unlensed_sampling_routine`                  | Generate unlensed parameters with batching     |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~unlensed_rate`                              | Calculate the unlensed detection rate          |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~lensed_cbc_statistics`                      | Generate lensed GW source parameters           |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~lensed_sampling_routine`                    | Generate lensed parameters with batching       |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~lensed_rate`                                | Calculate the lensed detection rate            |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~rate_function`                              | General helper for rate calculation            |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~rate_ratio`                                 | Calculate lensed/unlensed rate ratio           |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~rate_comparison_with_rate_calculation`      | Calculate and compare lensed/unlensed rates    |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~selecting_n_unlensed_detectable_events`     | Select n unlensed detectable events            |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~selecting_n_lensed_detectable_events`       | Select n lensed detectable events              |
    +-----------------------------------------------------+------------------------------------------------+


    Instance Attributes
    ----------
    LeR class has the following attributes: \n
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | Attribute                                      | Type             | Unit  | Description                                    |
    +================================================+==================+=======+================================================+
    | :meth:`~npool`                                 | ``int``          |       | Number of parallel processing cores            |
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | :meth:`~z_min`                                 | ``float``        |       | Minimum source redshift                        |
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | :meth:`~z_max`                                 | ``float``        |       | Maximum source redshift                        |
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | :meth:`~event_type`                            | ``str``          |       | Type of CBC event (BBH, BNS, NSBH)             |
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | :meth:`~lens_type`                             | ``str``          |       | Type of lens galaxy model                      |
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | :meth:`~cosmo`                                 | ``Cosmology``    |       | Astropy cosmology object                       |
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | :meth:`~json_file_names`                       | ``dict``         |       | JSON file names for parameter storage          |
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | :meth:`~interpolator_directory`                | ``str``          |       | Directory for interpolator files               |
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | :meth:`~ler_directory`                         | ``str``          |       | Directory for output parameter files           |
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | :meth:`~list_of_detectors`                     | ``list``         |       | List of detector names                         |
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | :meth:`~pdet_finder`                           | ``callable``     |       | Detection probability finder function          |
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | :meth:`~ler_args`                              | ``dict``         |       | All LeR initialization arguments               |
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | :meth:`~create_new_interpolator`               | ``dict``         |       | Interpolator creation settings                 |
    +------------------------------------------------+------------------+-------+------------------------------------------------+

    Notes
    -----
    - ``LeR`` class inherits from :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution`. \n
      Refer to that class for additional inherited attributes and methods. \n
    - Parameters are stored in JSON files for memory efficiency and reproducibility. \n
    - For stable rate estimates, use size >= 1e6 samples. \n
    """

    def __init__(
        self,
        npool=int(4),
        z_min=0.0,
        z_max=10.0,
        event_type="BBH",
        lens_type="epl_shear_galaxy",
        cosmology=None,
        pdet_finder=None,  # if not given, gwsnr's pdet calculator will be used
        json_file_names=None,
        interpolator_directory="./interpolator_json",
        create_new_interpolator=False,
        ler_directory="./ler_data",
        verbose=True,
        **kwargs,
    ):

        # getting interpolator data from the package
        # first check if the interpolator directory './interpolator_json' exists
        if not pathlib.Path(interpolator_directory).exists():
            # Get the path to the zip resource using importlib_resources
            zip_resource = resources_files('ler.rates').joinpath('ler_data', 'interpolator_json.zip')
            with zip_resource.open('rb') as zip_file:
                print("Extracting interpolator data from package to the current working directory.")

                # Define destination path (current working directory)
                dest_path = pathlib.Path.cwd()

                # Extract the zip file, skipping __MACOSX metadata
                with zipfile.ZipFile(zip_file, 'r') as zip_ref:
                    for member in zip_ref.namelist():
                        # Skip __MACOSX directory and its contents
                        if member.startswith('__MACOSX'):
                            continue
                        zip_ref.extract(member, dest_path)

        print("\nInitializing LeR class...\n")
        # init ler attributes
        self.npool = npool
        self.z_min = z_min
        self.z_max = z_max
        self.event_type = event_type
        self.lens_type = lens_type
        self.cosmo = cosmology if cosmology else LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

        # init json file names where datas will be stored
        self.json_file_names = dict(
            ler_params="ler_params.json",
            unlensed_param="unlensed_param.json",
            unlensed_param_detectable="unlensed_param_detectable.json",
            lensed_param="lensed_param.json",
            lensed_param_detectable="lensed_param_detectable.json",
        )
        if json_file_names:
            self.json_file_names.update(json_file_names)

        # init interpolator directory
        self.interpolator_directory = interpolator_directory
        # kwargs will be passed as input for the parent class
        kwargs["create_new_interpolator"] = create_new_interpolator
        self.ler_directory = ler_directory
        # create directory if not exists
        if not os.path.exists(ler_directory):
            os.makedirs(ler_directory)

        # parent class initialization
        self._parent_class_initialization(
            params=kwargs, pdet_finder=pdet_finder, verbose=verbose
        )

    def _parent_class_initialization(self, params=None, pdet_finder=None, verbose=True):
        """
        Function to initialize the parent classes.

        Parameters
        ----------
        params : ``dict``
            dictionary of parameters to initialize the parent classes
        """

        def initialization():
            # initialization of parent class
            self._parent_initialization_helper(params=params)
            # initialization self.snr and self.pdet_finder from GWSNR class
            if not pdet_finder:
                self.pdet_finder = self._gwsnr_initialization(params=params)
            else:
                self.pdet_finder = pdet_finder
                self.list_of_detectors = None

            # store all the ler input parameters
            self._store_ler_params(output_jsonfile=self.json_file_names["ler_params"])

        # if not verbose, prevent anything from printing
        if verbose:
            initialization()
            self._print_all_init_args()
        else:
            with contextlib.redirect_stdout(None):
                initialization()

    def _parent_initialization_helper(self, params=None):
        """
        Function to initialize the parent classes.

        Parameters
        ----------
        params : ``dict``
            dictionary of parameters to initialize the parent classes
        """

        # initialization of LensGalaxyParameterDistribution class
        # it also initializes the CBCSourceParameterDistribution and ImageProperties classes
        input_params = dict(
            # LensGalaxyParameterDistribution class params
            #    OpticalDepth class params
            npool=self.npool,
            z_min=self.z_min,
            z_max=self.z_max,
            cosmology=self.cosmo,
            lens_type=self.lens_type,
            lens_functions=None,
            lens_functions_params=None,
            lens_param_samplers=None,
            lens_param_samplers_params=None,
            directory=self.interpolator_directory,
            create_new_interpolator=False,
            # ImageProperties class params
            n_min_images=2,
            n_max_images=4,
            time_window=365 * 24 * 3600 * 20,
            lens_model_list=["EPL_NUMBA", "SHEAR"],
            # CBCSourceParameterDistribution class params
            event_type=self.event_type,
            source_priors=None,
            source_priors_params=None,
            spin_zero=False,
            spin_precession=False,
        )
        # update input_params with params. This will include create_new_interpolator.
        if params:
            for key, value in params.items():
                if key in input_params:
                    input_params[key] = value

        # initialization of parent class
        LensGalaxyParameterDistribution.__init__(
            self,
            # LensGalaxyParameterDistribution class params
            #    OpticalDepth class params
            npool=input_params["npool"],
            z_min=input_params["z_min"],
            z_max=input_params["z_max"],
            cosmology=input_params["cosmology"],
            lens_type=input_params["lens_type"],
            lens_functions=input_params["lens_functions"],
            lens_functions_params=input_params["lens_functions_params"],
            lens_param_samplers=input_params["lens_param_samplers"],
            lens_param_samplers_params=input_params["lens_param_samplers_params"],
            directory=input_params["directory"],
            create_new_interpolator=input_params["create_new_interpolator"],
            # ImageProperties class params
            n_min_images=input_params["n_min_images"],
            n_max_images=input_params["n_max_images"],
            time_window=input_params["time_window"],
            lens_model_list=input_params["lens_model_list"],
            # CBCSourceParameterDistribution class params
            event_type=input_params["event_type"],
            source_priors=input_params["source_priors"],
            source_priors_params=input_params["source_priors_params"],
            spin_zero=input_params["spin_zero"],
            spin_precession=input_params["spin_precession"],
        )

        # some of the None values will have default values after initialization
        input_params["source_priors"] = self.gw_param_samplers.copy()
        input_params["source_priors_params"] = self.gw_param_samplers_params.copy()
        input_params["lens_param_samplers"] = self.lens_param_samplers.copy()
        input_params["lens_param_samplers_params"] = (
            self.lens_param_samplers_params.copy()
        )
        input_params["lens_functions"] = self.lens_functions.copy()
        input_params["lens_functions_params"] = self.lens_functions_params.copy()
        input_params["create_new_interpolator"] = self.create_new_interpolator

        # save input_params to self.ler_args
        self.ler_args = input_params

    def _gwsnr_initialization(self, params=None):
        """
        Function to initialize the GWSNR class from the `gwsnr` package.

        Parameters
        ----------
        params : ``dict``
            dictionary of parameters to initialize the gwsnr class
        """
        from gwsnr import GWSNR

        # initialization of GWSNR class
        if "mminbh" in self.gw_param_samplers_params["source_frame_masses"]:
            min_bh_mass = self.gw_param_samplers_params["source_frame_masses"]["mminbh"]
        else:
            min_bh_mass = 2.0

        if "mmaxbh" in self.gw_param_samplers_params["source_frame_masses"]:
            max_bh_mass = self.gw_param_samplers_params["source_frame_masses"]["mmaxbh"]
        else:
            max_bh_mass = 200.0
        input_params = dict(
            # General settings
            npool=self.npool,
            snr_method="interpolation_aligned_spins",
            snr_type="optimal_snr",
            gwsnr_verbose=True,
            multiprocessing_verbose=True,
            pdet_kwargs=None,
            # Settings for interpolation grid
            mtot_min=min_bh_mass * 2,
            mtot_max=(
                max_bh_mass * 2 * (1 + self.z_max)
                if max_bh_mass * 2 * (1 + self.z_max) < 500.0
                else 500.0
            ),
            ratio_min=0.1,
            ratio_max=1.0,
            spin_max=0.99,
            mtot_resolution=200,
            ratio_resolution=20,
            spin_resolution=10,
            batch_size_interpolation=1000000,
            interpolator_dir=self.directory,
            create_new_interpolator=False,
            # GW signal settings
            sampling_frequency=2048.0,
            waveform_approximant="IMRPhenomD",
            frequency_domain_source_model="lal_binary_black_hole",
            minimum_frequency=20.0,
            reference_frequency=None,
            duration_max=None,
            duration_min=None,
            fixed_duration=None,
            mtot_cut=False,
            # Detector settings
            psds=None,
            ifos=None,
            noise_realization=None,  # not implemented yet
            # ANN settings
            ann_path_dict=None,
            # Hybrid SNR recalculation settings
            snr_recalculation=False,
            snr_recalculation_range=[6, 14],
            snr_recalculation_waveform_approximant="IMRPhenomXPHM",
        )
        # update input_params with params. This will include create_new_interpolator.
        if params:
            for key, value in params.items():
                if key in input_params:
                    input_params[key] = value
        self.ler_args["pdet_args"] = input_params

        # dealing with create_new_interpolator param
        if isinstance(input_params["create_new_interpolator"], bool):
            pass
        elif isinstance(input_params["create_new_interpolator"], dict):
            # check input_params["gwsnr"] exists
            if "gwsnr" in input_params["create_new_interpolator"]:
                if isinstance(input_params["create_new_interpolator"]["gwsnr"], bool):
                    input_params["create_new_interpolator"] = input_params[
                        "create_new_interpolator"
                    ]["gwsnr"]
                else:
                    raise ValueError(
                        "create_new_interpolator['gwsnr'] should be a boolean."
                    )
            else:
                input_params["create_new_interpolator"] = False

        # initialization of GWSNR class
        gwsnr = GWSNR(
            # General settings
            npool=input_params["npool"],
            snr_method=input_params["snr_method"],
            snr_type=input_params["snr_type"],
            gwsnr_verbose=input_params["gwsnr_verbose"],
            multiprocessing_verbose=input_params["multiprocessing_verbose"],
            pdet_kwargs=input_params["pdet_kwargs"],
            # Settings for interpolation grid
            mtot_min=input_params["mtot_min"],
            mtot_max=input_params["mtot_max"],
            ratio_min=input_params["ratio_min"],
            ratio_max=input_params["ratio_max"],
            spin_max=input_params["spin_max"],
            mtot_resolution=input_params["mtot_resolution"],
            ratio_resolution=input_params["ratio_resolution"],
            spin_resolution=input_params["spin_resolution"],
            batch_size_interpolation=input_params["batch_size_interpolation"],
            interpolator_dir=input_params["interpolator_dir"],
            create_new_interpolator=input_params["create_new_interpolator"],
            # GW signal settings
            sampling_frequency=input_params["sampling_frequency"],
            waveform_approximant=input_params["waveform_approximant"],
            frequency_domain_source_model=input_params["frequency_domain_source_model"],
            minimum_frequency=input_params["minimum_frequency"],
            reference_frequency=input_params["reference_frequency"],
            duration_max=input_params["duration_max"],
            duration_min=input_params["duration_min"],
            fixed_duration=input_params["fixed_duration"],
            mtot_cut=input_params["mtot_cut"],
            # Detector settings
            psds=input_params["psds"],
            ifos=input_params["ifos"],
            noise_realization=input_params["noise_realization"],
            # ANN settings
            ann_path_dict=input_params["ann_path_dict"],
            # Hybrid SNR recalculation settings
            snr_recalculation=input_params["snr_recalculation"],
            snr_recalculation_range=input_params["snr_recalculation_range"],
            snr_recalculation_waveform_approximant=input_params[
                "snr_recalculation_waveform_approximant"
            ],
        )

        self.ler_args["list_of_detectors"] = gwsnr.detector_list
        self.list_of_detectors = self.ler_args["list_of_detectors"]

        self.ler_args["pdet_args"]["pdet_kwargs"] = gwsnr.pdet_kwargs
        self.ler_args["pdet_args"]["psds_list"] = gwsnr.psds_list

        return gwsnr.pdet

    def _print_all_init_args(self):
        """
        Function to print all the parameters.
        """

        # print all relevant functions and sampler priors
        print("#-------------------------------------")
        print("# LeR initialization input arguments:")
        print("#-------------------------------------\n")
        print("    # LeR set up input arguments:")
        print(f"    npool = {self.npool},")
        print(f"    z_min = {self.z_min},")
        print(f"    z_max = {self.z_max},")
        print(f"    event_type = '{self.event_type}',")
        print(f"    lens_type = '{self.lens_type}',")
        print(f"    cosmology = {self.cosmo},")
        print(f"    pdet_finder = {self.pdet_finder},")
        print("    json_file_names = dict(")
        for key, value in self.json_file_names.items():
            (
                print(f"        {key} = '{value}',")
                if isinstance(value, str)
                else print(f"        {key} = {value},")
            )
        print("    ),")
        print(f"    interpolator_directory = '{self.interpolator_directory}',")
        print(f"    ler_directory = '{self.ler_directory}',")
        print("    create_new_interpolator = dict(")
        for key, value in self.ler_args["create_new_interpolator"].items():
            (
                print(f"        {key} = '{value}',")
                if isinstance(value, str)
                else print(f"        {key} = {value},")
            )
        print("    ),")

        print(
            "\n    # LeR also takes other CBCSourceParameterDistribution class input arguments as kwargs, as follows:"
        )
        print("    source_priors = dict(")
        for key, value in self.ler_args["source_priors"].items():
            (
                print(f"        {key} = '{value}',")
                if isinstance(value, str)
                else print(f"        {key} = {value},")
            )
        print("    ),")
        print("    source_priors_params = dict(")
        for key, value in self.ler_args["source_priors_params"].items():
            (
                print(f"        {key} = '{value}',")
                if isinstance(value, str)
                else print(f"        {key} = {value},")
            )
        print("    ),")
        print(f"    spin_zero = {self.ler_args['spin_zero']},")
        print(f"    spin_precession = {self.ler_args['spin_precession']},")

        print(
            "\n    # LeR also takes other LensGalaxyParameterDistribution class input arguments as kwargs, as follows:"
        )
        print("    lens_functions = dict(")
        for key, value in self.ler_args["lens_functions"].items():
            (
                print(f"        {key} = '{value}',")
                if isinstance(value, str)
                else print(f"        {key} = {value},")
            )
        print("    ),")
        print("    lens_functions_params = dict(")
        for key, value in self.ler_args["lens_functions_params"].items():
            (
                print(f"        {key} = '{value}',")
                if isinstance(value, str)
                else print(f"        {key} = {value},")
            )
        print("    ),")
        print("    lens_param_samplers = dict(")
        for key, value in self.ler_args["lens_param_samplers"].items():
            (
                print(f"        {key} = '{value}',")
                if isinstance(value, str)
                else print(f"        {key} = {value},")
            )
        print("    ),")
        print("    lens_param_samplers_params = dict(")
        for key, value in self.ler_args["lens_param_samplers_params"].items():
            (
                print(f"        {key} = '{value}',")
                if isinstance(value, str)
                else print(f"        {key} = {value},")
            )
        print("    ),")

        print(
            "\n    # LeR also takes other ImageProperties class input arguments as kwargs, as follows:"
        )
        print(f"    n_min_images = {self.ler_args['n_min_images']},")
        print(f"    n_max_images = {self.ler_args['n_max_images']},")
        print(f"    time_window = {self.ler_args['time_window']},")
        print(f"    lens_model_list = {self.ler_args['lens_model_list']},")

        if "pdet_args" in self.ler_args:
            print(
                "\n    # LeR also takes other gwsnr.GWSNR input arguments as kwargs, as follows:"
            )
            for key, value in self.ler_args["pdet_args"].items():
                if isinstance(value, dict):
                    print(f"    {key} = dict(")
                    for k, v in value.items():
                        (
                            print(f"        {k} = '{v}',")
                            if isinstance(v, str)
                            else print(f"        {k} = {v},")
                        )
                    print("    ),")
                else:
                    (
                        print(f"    {key} = '{value}',")
                        if isinstance(value, str)
                        else print(f"    {key} = {value},")
                    )

    def _store_ler_params(self, output_jsonfile):
        """
        Function to store the all the necessary parameters. This is useful for reproducing the results. All the parameters stored are in string format to make it json compatible.

        Parameters
        ----------
        output_jsonfile : ``str``
            name of the json file to store the parameters
        """

        # ler params
        param_sampler_dict = self.ler_args.copy()
        # convert all dict values to str
        for key, value in param_sampler_dict.items():
            param_sampler_dict[key] = str(value)

        file_name = self.json_file_names["ler_params"]
        append_json(
            self.ler_directory + "/" + file_name, param_sampler_dict, replace=True
        )

    def unlensed_cbc_statistics(
        self,
        size=100000,
        batch_size=50000,
        resume=True,
        save_batch=False,
        output_jsonfile=None,
    ):
        """
        Generate unlensed GW source parameters.

        This function calls the unlensed_sampling_routine function to generate
        the parameters in batches. The generated parameters are stored in a JSON
        file; and if save_batch=True, it keeps updating the file in batches.

        Parameters
        ----------
        size : ``int``
            Number of samples to generate. \n
            default: 100000
        batch_size : ``int``
            Batch size for sampling. \n
            default: 50000
        resume : ``bool``
            If True, the function will resume from the last batch. \n
            default: True
        save_batch : ``bool``
            If True, saves parameters in batches during sampling. \n
            If False, saves all parameters at the end (faster). \n
            default: False
        output_jsonfile : ``str``
            JSON file name for storing the parameters. \n
            default: None (uses self.json_file_names["unlensed_param"])

        Returns
        -------
        unlensed_param : ``dict``
            Dictionary of unlensed GW source parameters. The included parameters and their units are as follows (for default settings):\n
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
            | pdet_L1            |              | pdet of L1                           |
            +--------------------+--------------+--------------------------------------+
            | pdet_H1            |              | pdet of H1                           |
            +--------------------+--------------+--------------------------------------+
            | pdet_V1            |              | pdet of V1                           |
            +--------------------+--------------+--------------------------------------+
            | pdet_net           |              | pdet of the network                  |
            +--------------------+--------------+--------------------------------------+


        Examples
        --------
        >>> from ler import LeR
        >>> ler = LeR()
        >>> unlensed_param = ler.unlensed_cbc_statistics()
        """

        # Note: size must be provided as argument, no default fallback
        output_jsonfile = output_jsonfile or self.json_file_names["unlensed_param"]
        self.json_file_names["unlensed_param"] = output_jsonfile
        output_path = os.path.join(self.ler_directory, output_jsonfile)
        print(f"unlensed params will be stored in {output_path}")

        unlensed_param = batch_handler(
            size=size,
            batch_size=batch_size,
            sampling_routine=self.unlensed_sampling_routine,
            output_jsonfile=output_path,
            save_batch=save_batch,
            resume=resume,
            param_name="unlensed parameters",
        )

        return unlensed_param

    def unlensed_sampling_routine(
        self, size, output_jsonfile, resume=True, save_batch=True
    ):
        """
        Generate unlensed GW source parameters for a single batch.

        This is the core sampling routine called by unlensed_cbc_statistics.
        It samples GW source parameters and calculates detection probabilities.

        Parameters
        ----------
        size : ``int``
            Number of samples to generate.
        output_jsonfile : ``str``
            JSON file name for storing the parameters.
        resume : ``bool``
            If True, appends new samples to existing JSON file. \n
            default: True
        save_batch : ``bool``
            If True, saves parameters in batches during sampling. \n
            default: True

        Returns
        -------
        unlensed_param : ``dict``
            Dictionary of unlensed GW source parameters.
        """

        # get gw params
        print("sampling gw source params...")
        unlensed_param = self.sample_gw_parameters(size=size)

        # Get pdet
        print("calculating pdet...")
        pdet = self.pdet_finder(gw_param_dict=unlensed_param)
        unlensed_param.update(pdet)

        return unlensed_param

    def unlensed_rate(
        self,
        unlensed_param=None,
        pdet_threshold=0.5,
        pdet_type="boolean",
        output_jsonfile=None,
    ):
        """
        Function to calculate the unlensed rate.

        This function calculates the detection rate for unlensed events and stores
        the parameters of the detectable events in a JSON file.

        Parameters
        ----------
        unlensed_param : ``dict`` or ``str``
            Dictionary of GW source parameters or JSON file name. \n
            default: None (uses self.json_file_names["unlensed_param"])
        pdet_threshold : ``float``
            Threshold for detection probability. \n
            default: 0.5
        pdet_type : ``str``
            Detectability condition type. \n
            Options: \n
            - 'boolean': Binary detection based on pdet_threshold \n
            - 'probability_distribution': Uses pdet values directly \n
            default: 'boolean'
        output_jsonfile : ``str``
            JSON file name for storing the parameters of the detectable events. \n
            default: None (uses self.json_file_names["unlensed_param_detectable"])

        Returns
        -------
        total_rate : ``float``
            Total unlensed rate (yr^-1).
        unlensed_param : ``dict``
            Dictionary of unlensed GW source parameters of the detectable events.

        Examples
        --------
        >>> from ler import LeR
        >>> ler = LeR()
        >>> ler.unlensed_cbc_statistics()
        >>> total_rate, unlensed_param_detectable = ler.unlensed_rate()
        """

        unlensed_param = self._load_param(unlensed_param, param_type="unlensed")
        total_events = len(unlensed_param["zs"])

        # find index of detectable events
        pdet = unlensed_param["pdet_net"]
        idx_detectable = pdet > pdet_threshold

        if pdet_type == "boolean":
            detectable_events = np.sum(idx_detectable)
        elif pdet_type == "probability_distribution":
            detectable_events = np.sum(pdet)
        else:
            raise ValueError("pdet_type not recognized")
        # montecarlo integration``
        # The total rate R = norm <Theta(rho-rhoc)>
        total_rate = self.rate_function(
            detectable_events, total_events, param_type="unlensed"
        )

        # store all detectable params in json file
        self._save_detectable_params(
            output_jsonfile,
            unlensed_param,
            idx_detectable,
            key_file_name="unlensed_param_detectable",
            nan_to_num=False,
            verbose=True,
            replace_jsonfile=True,
        )

        # append ler_param and save it
        self._append_ler_param(total_rate, pdet_type=pdet_type, param_type="unlensed")

        return total_rate, unlensed_param

    def _load_param(self, param, param_type="unlensed"):
        """
        Helper function to load or copy unlensed/lensed parameters.

        Parameters
        ----------
        param : ``dict`` or ``str``
            dictionary of unlensed/lensed parameters or json file name.
        param_type : ``str``
            type of parameters.
            default param_type = 'unlensed'. Other options is 'lensed'.

        Returns
        ----------
        param : ``dict``
            dictionary of unlensed/lensed parameters.
        """

        param_type = param_type + "_param"
        if param is None:
            param = self.json_file_names[param_type]
        if isinstance(param, str):
            path_ = self.ler_directory + "/" + param
            print(f"Getting {param_type} from json file {path_}...")
            return get_param_from_json(path_)
        else:
            print(f"Using provided {param_type} dict...")
            return param.copy()

    def rate_function(
        self, detectable_size, total_size, param_type="unlensed", verbose=True
    ):
        """
        Calculate the detection rate for unlensed or lensed events.

        This is a general helper function that computes the rate based on
        Monte Carlo integration using the ratio of detectable to total events.

        Parameters
        ----------
        detectable_size : ``int`` or ``float``
            Number of detectable events (or sum of pdet values).
        total_size : ``int``
            Total number of simulated events.
        param_type : ``str``
            Type of parameters. \n
            Options: \n
            - 'unlensed': Use unlensed normalization \n
            - 'lensed': Use lensed normalization \n
            default: 'unlensed'
        verbose : ``bool``
            If True, print rate information. \n
            default: True

        Returns
        -------
        rate : ``float``
            Event rate (yr^-1).

        Examples
        --------
        >>> from ler import LeR
        >>> ler = LeR()
        >>> rate = ler.rate_function(detectable_size=100, total_size=1000)
        """

        if param_type == "unlensed":
            normalization = self.normalization_pdf_z
        elif param_type == "lensed":
            normalization = self.normalization_pdf_z_lensed
        rate = float(normalization * detectable_size / total_size)
        # print(f"\nnormalization factor for {param_type} rate: {normalization}")
        # print(f"detectable size: {detectable_size}")
        # print(f"total size: {total_size}")
        # print(f"calculated {param_type} rate: {rate}\n")

        if verbose:
            print(f"total {param_type} rate (yr^-1): {rate}")
            print(
                f"number of simulated {param_type} detectable events: {detectable_size}"
            )
            print(f"number of simulated all {param_type} events: {total_size}")

        return rate

    def _save_detectable_params(
        self,
        output_jsonfile,
        param,
        idx_detectable,
        key_file_name="unlensed_param_detectable",
        nan_to_num=False,
        verbose=True,
        replace_jsonfile=True,
    ):
        """
        Helper function to save the detectable parameters in json file.

        Parameters
        ----------
        output_jsonfile : ``str``
            json file name for storing the parameters of the detectable events. This is stored in the self.ler_directory.
        param : ``dict``
            dictionary of GW source parameters.
        idx_detectable : ``numpy.ndarray``
            index of detectable events.
        key_file_name : ``str``
            key name for the json file to be added in self.json_file_names.
        nan_to_num : ``bool``
            if True, it will replace nan with 0.
            default nan_to_num = False.
        verbose : ``bool``
            if True, it will print the path of the json file.
            default verbose = True.
        replace_jsonfile : ``bool``
            if True, it will replace the json file. If False, it will append the json file.
        """

        # store all detectable params in json file
        if nan_to_num:
            for key, value in param.items():
                param[key] = np.nan_to_num(value[idx_detectable])
        else:
            for key, value in param.items():
                param[key] = value[idx_detectable]

        # store all detectable params in json file
        if output_jsonfile is None:
            output_jsonfile = self.json_file_names[key_file_name]
        else:
            self.json_file_names[key_file_name] = output_jsonfile

        output_path = self.ler_directory + "/" + output_jsonfile
        if verbose:
            print(f"storing detectable params in {output_path}")
        append_json(output_path, param, replace=replace_jsonfile)

    def _append_ler_param(self, total_rate, pdet_type="boolean", param_type="unlensed"):
        """
        Helper function to append the final results, total_rate, in the json file.

        Parameters
        ----------
        total_rate : ``float``
            total rate.
        param_type : ``str``
            type of parameters.
            default param_type = 'unlensed'. Other options is 'lensed'.
        """

        data = load_json(self.ler_directory + "/" + self.json_file_names["ler_params"])
        # write the results
        data[f"detectable_{param_type}_rate_per_year"] = total_rate
        data[f"pdet_type_{param_type}"] = pdet_type
        append_json(
            self.ler_directory + "/" + self.json_file_names["ler_params"],
            data,
            replace=True,
        )

    def lensed_cbc_statistics(
        self,
        size=100000,
        batch_size=50000,
        save_batch=False,
        resume=True,
        output_jsonfile=None,
    ):
        """
        Generate lensed GW source parameters.

        This function calls the lensed_sampling_routine function to generate
        the parameters in batches. The generated parameters are stored in a JSON
        file; and if save_batch=True, it keeps updating the file in batches.

        Parameters
        ----------
        size : ``int``
            Number of samples to generate. \n
            default: 100000
        batch_size : ``int``
            Batch size for sampling. \n
            default: 50000
        save_batch : ``bool``
            If True, saves parameters in batches during sampling. \n
            If False, saves all parameters at the end (faster). \n
            default: True
        resume : ``bool``
            If True, the function will resume from the last batch. \n
            default: True
        output_jsonfile : ``str``
            JSON file name for storing the parameters. \n
            default: None (uses self.json_file_names["lensed_param"])

        Returns
        -------
        lensed_param : ``dict``
            Dictionary of lensed GW source parameters. The included parameters and their units are as follows (for default settings):\n
            +------------------------------+-----------+-------------------------------------------------------+
            | Parameter                    | Units     | Description                                           |
            +==============================+===========+=======================================================+
            | zl                           |           | redshift of the lens                                  |
            +------------------------------+-----------+-------------------------------------------------------+
            | zs                           |           | redshift of the source                                |
            +------------------------------+-----------+-------------------------------------------------------+
            | sigma                        | km s^-1   | velocity dispersion                                   |
            +------------------------------+-----------+-------------------------------------------------------+
            | q                            |           | axis ratio                                            |
            +------------------------------+-----------+-------------------------------------------------------+
            | theta_E                      | arcsec    | Einstein radius                                       |
            +------------------------------+-----------+-------------------------------------------------------+
            | phi                          | rad       | axis rotation angle. counter-clockwise from the       |
            |                              |           | positive x-axis (RA-like axis) to the major axis of   |
            |                              |           | the projected mass distribution.                      |
            +------------------------------+-----------+-------------------------------------------------------+
            | gamma                        |           | density profile slope of EPL galaxy                   |
            +------------------------------+-----------+-------------------------------------------------------+
            | gamma1                       |           | external shear component in the x-direction           |
            |                              |           | (RA-like axis)                                        |
            +------------------------------+-----------+-------------------------------------------------------+
            | gamma2                       |           | external shear component in the y-direction           |
            |                              |           | (Dec-like axis)                                       |
            +------------------------------+-----------+-------------------------------------------------------+
            | geocent_time                 | s         | geocent time                                          |
            +------------------------------+-----------+-------------------------------------------------------+
            | ra                           | rad       | right ascension                                       |
            +------------------------------+-----------+-------------------------------------------------------+
            | dec                          | rad       | declination                                           |
            +------------------------------+-----------+-------------------------------------------------------+
            | phase                        | rad       | phase of GW at reference freq                         |
            +------------------------------+-----------+-------------------------------------------------------+
            | psi                          | rad       | polarization angle                                    |
            +------------------------------+-----------+-------------------------------------------------------+
            | theta_jn                     | rad       | inclination angle                                     |
            +------------------------------+-----------+-------------------------------------------------------+
            | a_1                          |           | spin of the primary compact binary                    |
            +------------------------------+-----------+-------------------------------------------------------+
            | a_2                          |           | spin of the secondary compact binary                  |
            +------------------------------+-----------+-------------------------------------------------------+
            | mass_1_source                | Msun      | mass of the primary compact binary (source frame)     |
            +------------------------------+-----------+-------------------------------------------------------+
            | mass_2_source                | Msun      | mass of the secondary compact binary (source frame)   |
            +------------------------------+-----------+-------------------------------------------------------+
            | mass_1                       | Msun      | mass of the primary compact binary (detector frame)   |
            +------------------------------+-----------+-------------------------------------------------------+
            | mass_2                       | Msun      | mass of the secondary compact binary (detector frame) |
            +------------------------------+-----------+-------------------------------------------------------+
            | x0_image_positions           | arcsec    | x-coordinate (RA-like axis) of the images             |
            +------------------------------+-----------+-------------------------------------------------------+
            | x1_image_positions           | arcsec    | y-coordinate (Dec-like axis) of the images            |
            +------------------------------+-----------+-------------------------------------------------------+
            | magnifications               |           | magnifications                                        |
            +------------------------------+-----------+-------------------------------------------------------+
            | time_delays                  |           | time delays                                           |
            +------------------------------+-----------+-------------------------------------------------------+
            | image_type                   |           | image type                                            |
            +------------------------------+-----------+-------------------------------------------------------+
            | n_images                     |           | number of images                                      |
            +------------------------------+-----------+-------------------------------------------------------+
            | x_source                     | arcsec    | x-coordinate (RA-like axis) of the source             |
            +------------------------------+-----------+-------------------------------------------------------+
            | y_source                     | arcsec    | y-coordinate (Dec-like axis) of the source            |
            +------------------------------+-----------+-------------------------------------------------------+
            | effective_luminosity_distance| Mpc       | effective luminosity distance of the images           |
            +------------------------------+-----------+-------------------------------------------------------+
            | effective_geocent_time       | s         | effective GPS time of coalescence of the images       |
            +------------------------------+-----------+-------------------------------------------------------+
            | pdet_L1                      |           | detection probability of L1                           |
            +------------------------------+-----------+-------------------------------------------------------+
            | pdet_H1                      |           | detection probability of H1                           |
            +------------------------------+-----------+-------------------------------------------------------+
            | pdet_V1                      |           | detection probability of V1                           |
            +------------------------------+-----------+-------------------------------------------------------+
            | pdet_net                     |           | detection probability of the network                  |
            +------------------------------+-----------+-------------------------------------------------------+


        Examples
        --------
        >>> from ler import LeR
        >>> ler = LeR()
        >>> lensed_param = ler.lensed_cbc_statistics()
        """

        # Note: size must be provided as argument, no default fallback
        output_jsonfile = output_jsonfile or self.json_file_names["lensed_param"]
        self.json_file_names["lensed_param"] = output_jsonfile
        output_path = os.path.join(self.ler_directory, output_jsonfile)
        print(f"lensed params will be stored in {output_path}")

        lensed_param = batch_handler(
            size=size,
            batch_size=batch_size,
            sampling_routine=self.lensed_sampling_routine,
            output_jsonfile=output_path,
            save_batch=save_batch,
            resume=resume,
            param_name="lensed parameters",
        )

        return lensed_param

    def lensed_sampling_routine(
        self, size, output_jsonfile, save_batch=True, resume=True
    ):
        """
        Generate lensed GW source parameters for a single batch.

        This is the core sampling routine called by lensed_cbc_statistics.
        It samples lens parameters, calculates image properties, and computes
        detection probabilities for the images of lensed events.

        Parameters
        ----------
        size : ``int``
            Number of samples to generate.
        output_jsonfile : ``str``
            JSON file name for storing the parameters.
        save_batch : ``bool``
            If True, saves parameters in batches during sampling. \n
            default: True
        resume : ``bool``
            If True, appends new samples to existing JSON file. \n
            default: True

        Returns
        -------
        lensed_param : ``dict``
            Dictionary of lensed GW source parameters.
        """

        print("sampling lensed params...")
        lensed_param = {}

        # Some of the sample lensed events may not satisfy the strong lensing condition
        # In that case, we will resample those events and replace the values with the corresponding indices
        while True:
            # get lensed params
            lensed_param_ = self.sample_lens_parameters(size=size)
            # now get (strongly lensed) image paramters along with lens parameters
            lensed_param_ = self.image_properties(lensed_param_)

            if len(lensed_param) == 0:  # if empty
                lensed_param = lensed_param_
            else:
                # below will not be used in the first iteration
                # replace the values of the keys
                # idx defines the position that does not satisfy the strong lensing condition
                for key, value in lensed_param_.items():
                    lensed_param[key][idx] = value  # noqa: F821

            # check for invalid samples
            idx = (lensed_param["n_images"] < self.n_min_images) | (
                lensed_param["n_images"] > self.n_max_images
            )

            if np.sum(idx) == 0:
                break
            else:
                print(
                    f"Invalid sample found. Resampling {np.sum(idx)} lensed events..."
                )
                size = np.sum(idx)

        print("calculating pdet...")
        pdet, lensed_param = self.get_lensed_snrs(
            lensed_param=lensed_param,
            list_of_detectors=self.list_of_detectors,
            pdet_calculator=self.pdet_finder,
        )
        lensed_param.update(pdet)

        return lensed_param

    def lensed_rate(
        self,
        lensed_param=None,
        pdet_threshold=[0.5, 0.5],
        num_img=[1, 1],
        output_jsonfile=None,
        nan_to_num=True,
        pdet_type="boolean",
    ):
        """
        Function to calculate the lensed rate.

        This function calculates the detection rate for lensed events and stores
        the parameters of the detectable events in a JSON file.

        Parameters
        ----------
        lensed_param : ``dict`` or ``str``
            Dictionary of lensed GW source parameters or JSON file name. \n
            default: None (uses self.json_file_names["lensed_param"])
        pdet_threshold : ``float`` or ``list``
            Threshold for detection probability. \n
            default: [0.5, 0.5]
        num_img : ``int`` or ``list``
            Number of images corresponding to the pdet_threshold. \n
            Together with pdet_threshold = [0.5, 0.5], it means that two images with pdet > 0.5. \n
            Same condition can also be represented by pdet_threshold = 0.5 and num_img = 2. \n
            default: [1, 1]
        output_jsonfile : ``str``
            JSON file name for storing the parameters of the detectable events. \n
            default: None (uses self.json_file_names["lensed_param_detectable"])
        nan_to_num : ``bool``
            If True, NaN values will be converted to 0. \n
            default: True
        pdet_type : ``str``
            Detectability condition type. \n
            Options: \n
            - 'boolean': Binary detection based on pdet_threshold \n
            - 'probability_distribution': Uses pdet values directly \n
            default: 'boolean'

        Returns
        -------
        total_rate : ``float``
            Total lensed rate (yr^-1).
        lensed_param : ``dict``
            Dictionary of lensed GW source parameters of the detectable events. The included parameters and their units are as follows (for default settings):\n
            +------------------------------+-----------+-------------------------------------------------------+
            | Parameter                    | Units     | Description                                           |
            +==============================+===========+=======================================================+
            | zl                           |           | redshift of the lens                                  |
            +------------------------------+-----------+-------------------------------------------------------+
            | zs                           |           | redshift of the source                                |
            +------------------------------+-----------+-------------------------------------------------------+
            | sigma                        | km s^-1   | velocity dispersion                                   |
            +------------------------------+-----------+-------------------------------------------------------+
            | q                            |           | axis ratio                                            |
            +------------------------------+-----------+-------------------------------------------------------+
            | theta_E                      | arcsec    | Einstein radius                                       |
            +------------------------------+-----------+-------------------------------------------------------+
            | phi                          | rad       | axis rotation angle. counter-clockwise from the       |
            |                              |           | positive x-axis (RA-like axis) to the major axis of   |
            |                              |           | the projected mass distribution.                      |
            +------------------------------+-----------+-------------------------------------------------------+
            | gamma                        |           | density profile slope of EPL galaxy                   |
            +------------------------------+-----------+-------------------------------------------------------+
            | gamma1                       |           | external shear component in the x-direction           |
            |                              |           | (RA-like axis)                                        |
            +------------------------------+-----------+-------------------------------------------------------+
            | gamma2                       |           | external shear component in the y-direction           |
            |                              |           | (Dec-like axis)                                       |
            +------------------------------+-----------+-------------------------------------------------------+
            | geocent_time                 | s         | geocent time                                          |
            +------------------------------+-----------+-------------------------------------------------------+
            | ra                           | rad       | right ascension                                       |
            +------------------------------+-----------+-------------------------------------------------------+
            | dec                          | rad       | declination                                           |
            +------------------------------+-----------+-------------------------------------------------------+
            | phase                        | rad       | phase of GW at reference freq                         |
            +------------------------------+-----------+-------------------------------------------------------+
            | psi                          | rad       | polarization angle                                    |
            +------------------------------+-----------+-------------------------------------------------------+
            | theta_jn                     | rad       | inclination angle                                     |
            +------------------------------+-----------+-------------------------------------------------------+
            | a_1                          |           | spin of the primary compact binary                    |
            +------------------------------+-----------+-------------------------------------------------------+
            | a_2                          |           | spin of the secondary compact binary                  |
            +------------------------------+-----------+-------------------------------------------------------+
            | mass_1_source                | Msun      | mass of the primary compact binary (source frame)     |
            +------------------------------+-----------+-------------------------------------------------------+
            | mass_2_source                | Msun      | mass of the secondary compact binary (source frame)   |
            +------------------------------+-----------+-------------------------------------------------------+
            | mass_1                       | Msun      | mass of the primary compact binary (detector frame)   |
            +------------------------------+-----------+-------------------------------------------------------+
            | mass_2                       | Msun      | mass of the secondary compact binary (detector frame) |
            +------------------------------+-----------+-------------------------------------------------------+
            | x0_image_positions           | arcsec    | x-coordinate (RA-like axis) of the images             |
            +------------------------------+-----------+-------------------------------------------------------+
            | x1_image_positions           | arcsec    | y-coordinate (Dec-like axis) of the images            |
            +------------------------------+-----------+-------------------------------------------------------+
            | magnifications               |           | magnifications                                        |
            +------------------------------+-----------+-------------------------------------------------------+
            | time_delays                  |           | time delays                                           |
            +------------------------------+-----------+-------------------------------------------------------+
            | image_type                   |           | image type                                            |
            +------------------------------+-----------+-------------------------------------------------------+
            | n_images                     |           | number of images                                      |
            +------------------------------+-----------+-------------------------------------------------------+
            | x_source                     | arcsec    | x-coordinate (RA-like axis) of the source             |
            +------------------------------+-----------+-------------------------------------------------------+
            | y_source                     | arcsec    | y-coordinate (Dec-like axis) of the source            |
            +------------------------------+-----------+-------------------------------------------------------+
            | effective_luminosity_distance| Mpc       | effective luminosity distance of the images           |
            +------------------------------+-----------+-------------------------------------------------------+
            | effective_geocent_time       | s         | effective GPS time of coalescence of the images       |
            +------------------------------+-----------+-------------------------------------------------------+
            | pdet_L1                      |           | detection probability of L1                           |
            +------------------------------+-----------+-------------------------------------------------------+
            | pdet_H1                      |           | detection probability of H1                           |
            +------------------------------+-----------+-------------------------------------------------------+
            | pdet_V1                      |           | detection probability of V1                           |
            +------------------------------+-----------+-------------------------------------------------------+
            | pdet_net                     |           | detection probability of the network                  |
            +------------------------------+-----------+-------------------------------------------------------+


        Examples
        --------
        >>> from ler import LeR
        >>> ler = LeR()
        >>> ler.lensed_cbc_statistics()
        >>> total_rate, lensed_param_detectable = ler.lensed_rate()
        """

        # load lensed parameters
        lensed_param = self._load_param(lensed_param, param_type="lensed")

        # re-analyse the provided pdet_threshold and num_img
        pdet_threshold, num_img = self._check_pdet_threshold_lensed(
            pdet_threshold, num_img
        )

        # get size of the lensed_param for a parameter
        total_events = len(lensed_param["zs"])

        # find index of detectable events
        pdet_hit, pdet_prod = self._find_detectable_index_lensed(
            lensed_param, pdet_threshold, num_img, pdet_type
        )

        # montecarlo integration
        pdet_processed = pdet_hit if pdet_type == "boolean" else pdet_prod
        total_rate = self.rate_function(
            np.sum(pdet_processed), total_events, param_type="lensed"
        )

        # store all detectable params in json file
        self._save_detectable_params(
            output_jsonfile,
            lensed_param,
            pdet_hit,
            key_file_name="lensed_param_detectable",
            nan_to_num=nan_to_num,
            verbose=True,
            replace_jsonfile=True,
        )

        # append ler_param and save it
        self._append_ler_param(total_rate, pdet_type=pdet_type, param_type="lensed")

        return total_rate, lensed_param

    def _check_pdet_threshold_lensed(self, pdet_threshold, num_img):
        """
        Helper function to check the pdet_threshold and num_img for lensed events.

        Parameters
        ----------
        pdet_threshold : ``float``
            threshold for detection signal to noise ratio.
            default pdet_threshold = [0.5,0.5].
        num_img : `int`
            number of images corresponding to the pdet_threshold.
            default num_img = [1,1]. Together with pdet_threshold = [0.5,0.5], it means that two images with pdet>0.5. Same condition can also be represented by pdet_threshold = 0.5 and num_img = 2.

        Returns
        ----------
        pdet_threshold : ``float``
            threshold for detection probability.
        num_img : `int`
            number of images corresponding to the snr_threshold.
        """
        # check for images with snr above threshold
        # convert to array
        pdet_threshold_ = np.array([pdet_threshold]).reshape(-1)
        num_img_ = np.array([num_img]).reshape(-1)
        # get descending sorted idx of snr_threshold
        idx = np.argsort(-pdet_threshold_)
        pdet_threshold = pdet_threshold_[idx]
        num_img = num_img_[idx]

        return pdet_threshold, num_img

    def _find_detectable_index_lensed(
        self, lensed_param, pdet_threshold, num_img, pdet_type
    ):
        """
        Helper function to find the index of detectable events based on SNR or p_det.

        Parameters
        ----------
        lensed_param : ``dict``
            dictionary of lensed GW source parameters.
        pdet_threshold : ``float`` or ``list``
            threshold for detection probability.
            default pdet_threshold = 0.5.
        num_img : ``int`` or ``list``
            number of images corresponding to the pdet_threshold.
            default num_img = [1,1].
        pdet_type : ``str``
            detectability condition.
            default pdet_type = 'boolean'.
            other options are 'pdet'.

        Returns
        ----------
        pdet_hit : ``bool``
            boolean array to store the result of the threshold condition.
        """

        if "pdet_net" not in lensed_param:
            raise ValueError("'pdet_net' not in lensed parm dict provided")

        # print(f"given pdet_type == {pdet_type}")
        if pdet_type == "boolean":

            pdet_param = lensed_param["pdet_net"]
            pdet_param = -np.sort(-pdet_param, axis=1)  # sort pdet in descending order
            pdet_hit = np.full(
                len(pdet_param), True
            )  # boolean array to store the result of the threshold condition

            # for each row: choose a threshold and check if the number of images above threshold. Sum over the images. If sum is greater than num_img, then snr_hit = True
            # algorithm:
            # i) consider pdet_threshold=[0.5,0.5] and num_img=[2,1] and first row of pdet_param[0]=[0.6,0.4,0.3,0.2]. Note that the pdet_param is sorted in descending order.
            # ii) for loop runs wrt pdet_threshold. idx_max = idx_max + num_img[i]
            # iii) First iteration: pdet_threshold=0.5 and num_img=2. In pdet_param, column index 0 and 1 (i.e. 0:num_img[0]) are considered. The sum of pdet_param[0, 0:2] > 0.5 is checked. If True, then pdet_hit = True.
            # v) Second iteration: pdet_threshold=0.5 and num_img=1. In pdet_param, column index 2 (i.e. num_img[0]:num_img[1]) is considered. The sum of pdet_param[0, 0:1] > 0.5 is checked. If True, then pdet_hit = True.
            j = 0
            idx_max = 0
            for i, pdet_th in enumerate(pdet_threshold):
                idx_max = idx_max + num_img[i]
                pdet_hit = pdet_hit & (
                    np.sum((pdet_param[:, j:idx_max] > pdet_th), axis=1) >= num_img[i]
                )
                # select according to time delays
                j = idx_max

            pdet_prod = None

        elif pdet_type == "probability_distribution":

            pdet = lensed_param["pdet_net"]
            # sort pdet in descending order
            pdet = -np.sort(-pdet, axis=1)
            # column index beyong np.sum(num_img)-1 are not considered
            pdet = pdet[:, : np.sum(num_img)]

            pdet_prod = np.prod(pdet, axis=1)

            pdet_hit = pdet_prod >= pdet_threshold


        return pdet_hit, pdet_prod

    def rate_comparison_with_rate_calculation(
        self,
        unlensed_param=None,
        pdet_threshold_unlensed=0.5,
        output_jsonfile_unlensed=None,
        lensed_param=None,
        pdet_threshold_lensed=[0.5, 0.5],
        num_img_lensed=[1, 1],
        output_jsonfile_lensed=None,
        nan_to_num=True,
        pdet_type="boolean",
    ):
        """
        Calculate and compare unlensed and lensed detection rates.

        This function calculates both unlensed and lensed rates and computes
        their ratio. It stores the parameters of the detectable events in JSON
        files. Using this function eliminates the need to call unlensed_rate
        and lensed_rate separately.

        Parameters
        ----------
        unlensed_param : ``dict`` or ``str``
            Dictionary of GW source parameters or JSON file name. \n
            default: None (uses self.json_file_names["unlensed_param"])
        pdet_threshold_unlensed : ``float``
            Detection probability threshold for unlensed events. \n
            default: 0.5
        output_jsonfile_unlensed : ``str``
            JSON file name for storing detectable unlensed parameters. \n
            default: None
        lensed_param : ``dict`` or ``str``
            Dictionary of lensed GW source parameters or JSON file name. \n
            default: None (uses self.json_file_names["lensed_param"])
        pdet_threshold_lensed : ``float`` or ``list``
            Detection probability threshold for lensed events. \n
            default: [0.5, 0.5]
        num_img_lensed : ``list``
            Number of images for lensed events. \n
            default: [1, 1]
        output_jsonfile_lensed : ``str``
            JSON file name for storing detectable lensed parameters. \n
            default: None
        nan_to_num : ``bool``
            If True, NaN values will be converted to 0. \n
            default: True
        pdet_type : ``str``
            Detectability condition type. \n
            Options: \n
            - 'boolean': Binary detection based on pdet_threshold \n
            - 'probability_distribution': Uses pdet values directly \n
            default: 'boolean'

        Returns
        -------
        rate_ratio : ``float``
            Ratio of unlensed rate to lensed rate.
        unlensed_param_detectable : ``dict``
            Dictionary of detectable unlensed GW source parameters.
        lensed_param_detectable : ``dict``
            Dictionary of detectable lensed GW source parameters.

        Examples
        --------
        >>> from ler import LeR
        >>> ler = LeR()
        >>> ler.unlensed_cbc_statistics()
        >>> ler.lensed_cbc_statistics()
        >>> rate_ratio, unlensed_param, lensed_param = ler.rate_comparison_with_rate_calculation()
        """

        # get unlensed rate
        unlensed_rate, unlensed_param_detectable = self.unlensed_rate(
            unlensed_param=unlensed_param,
            pdet_threshold=pdet_threshold_unlensed,
            pdet_type=pdet_type,
            output_jsonfile=output_jsonfile_unlensed,
        )
        # get lensed rate
        lensed_rate, lensed_param_detectable = self.lensed_rate(
            lensed_param=lensed_param,
            pdet_threshold=pdet_threshold_lensed,
            num_img=num_img_lensed,
            output_jsonfile=output_jsonfile_lensed,
            nan_to_num=nan_to_num,
            pdet_type=pdet_type,
        )
        # calculate rate ratio
        rate_ratio = self.rate_ratio()

        return (
            unlensed_rate,
            lensed_rate,
            rate_ratio,
            unlensed_param_detectable,
            lensed_param_detectable,
        )

    def rate_ratio(self):
        """
        Calculate and display the unlensed to lensed merger rate ratio.

        This function retrieves the unlensed_rate and lensed_rate from the
        JSON file specified in self.json_file_names["ler_params"] and computes
        their ratio.

        Returns
        -------
        rate_ratio : ``float``
            Ratio of unlensed rate to lensed rate.

        Examples
        --------
        >>> from ler import LeR
        >>> ler = LeR()
        >>> ler.unlensed_cbc_statistics()
        >>> ler.lensed_cbc_statistics()
        >>> ler.unlensed_rate()
        >>> ler.lensed_rate()
        >>> ler.rate_ratio()
        """

        # call json_file_ler_param and add the results
        data = load_json(self.ler_directory + "/" + self.json_file_names["ler_params"])

        try:
            unlensed_rate = data["detectable_unlensed_rate_per_year"]
            lensed_rate = data["detectable_lensed_rate_per_year"]
        except KeyError:
            raise ValueError(
                f"detectable_unlensed_rate_per_year or 'detectable_lensed_rate_per_year' not found in {self.json_file_names['ler_params']} json file. \nRun the functions 'unlensed_rate' and 'lensed_rate' to calculate and save the rates."
            )

        rate_ratio = unlensed_rate / lensed_rate
        # append the results
        data["rate_ratio"] = rate_ratio
        # write the results
        append_json(
            self.ler_directory + "/" + self.json_file_names["ler_params"],
            data,
            replace=True,
        )

        print(f"unlensed_rate: {unlensed_rate}")
        print(f"lensed_rate: {lensed_rate}")
        print(f"ratio: {rate_ratio}")

        return unlensed_rate / lensed_rate

    def selecting_n_unlensed_detectable_events(
        self,
        size=100,
        batch_size=50000,
        stopping_criteria=dict(
            relative_diff_percentage=0.5,
            number_of_last_batches_to_check=4,
        ),
        pdet_threshold=0.5,
        resume=True,
        output_jsonfile="n_unlensed_param_detectable.json",
        meta_data_file="meta_unlensed.json",
        pdet_type="boolean",
        trim_to_size=False,
    ):
        """
        Generate a target number of detectable unlensed events by sampling in batches, with the option to stop once the cumulative rate has stabilized.

        This function samples unlensed parameters and saves only the detectable
        events in a JSON file. It also records metadata including the total
        number of events and the cumulative rate.

        Parameters
        ----------
        size : ``int``
            Target number of detectable samples to collect. \n
            default: 100
        batch_size : ``int``
            Batch size for sampling. \n
            default: 50000
        stopping_criteria : ``dict`` or ``None``
            Criteria for stopping sample collection (but will not stop until n>size). \n
            Keys: \n
            - 'relative_diff_percentage': Maximum relative difference in rate (float) \n
            - 'number_of_last_batches_to_check': Number of batches for comparison (int) \n
            If None, stops when detectable events exceed size. \n
            default: dict(relative_diff_percentage=0.5, number_of_last_batches_to_check=4)
        pdet_threshold : ``float``
            Detection probability threshold. \n
            default: 0.5
        resume : ``bool``
            If True, resumes from last saved batch. \n
            default: True
        output_jsonfile : ``str``
            JSON file name for storing detectable parameters. \n
            default: 'n_unlensed_param_detectable.json'
        meta_data_file : ``str``
            JSON file name for storing metadata. \n
            default: 'meta_unlensed.json'
        pdet_type : ``str``
            Detectability condition type. \n
            Options: \n
            - 'boolean': Binary detection based on pdet_threshold \n
            - 'probability_distribution': Uses pdet values directly \n
            default: 'boolean'
        trim_to_size : ``bool``
            If True, trims final result to exactly size events. \n
            default: False

        Returns
        -------
        param_final : ``dict``
            Dictionary of unlensed GW source parameters of detectable events.

        Examples
        --------
        >>> from ler import LeR
        >>> ler = LeR()
        >>> unlensed_param = ler.selecting_n_unlensed_detectable_events(size=100)
        """

        if stopping_criteria is not None:
            print(
                f"stopping criteria set to when relative difference of total rate for the last {stopping_criteria['number_of_last_batches_to_check']} cumulative batches is less than {stopping_criteria['relative_diff_percentage']}%."
            )
            print(
                "sample collection will stop when the stopping criteria is met and number of detectable events exceeds the specified size."
            )
        else:
            print(
                "stopping criteria not set. sample collection will stop when number of detectable events exceeds the specified size."
            )

        # initial setup
        (
            n,
            events_total,
            output_path,
            meta_data_path,
            buffer_file,
            continue_condition,
        ) = self._initial_setup_for_n_event_selection(
            size, meta_data_file, output_jsonfile, resume, stopping_criteria
        )

        # loop until n==size samples are collected
        while continue_condition:
            # disable print statements
            with contextlib.redirect_stdout(None):
                self.dict_buffer = None  # this is used to store the sampled unlensed_param in batches when running the sampling_routine
                unlensed_param = self.unlensed_sampling_routine(
                    size=batch_size,
                    output_jsonfile=buffer_file,
                    save_batch=False,
                    resume=False,
                )

            total_events_in_this_iteration = len(unlensed_param["zs"])
            # below is use when the snr is calculated with 'ann' method of `gwsnr`

            # find index of detectable events
            pdet = unlensed_param["pdet_net"]
            idx_detectable = pdet > pdet_threshold

            # store all params in json file
            self._save_detectable_params(
                output_jsonfile,
                unlensed_param,
                idx_detectable,
                key_file_name="n_unlensed_detectable_events",
                nan_to_num=False,
                verbose=False,
                replace_jsonfile=False,
            )

            if pdet_type == "boolean":
                n += np.sum(idx_detectable)
            else:
                n += np.sum(pdet)
            events_total += total_events_in_this_iteration
            total_rate = self.rate_function(
                n, events_total, param_type="unlensed", verbose=False
            )

            # bookmark
            buffer_dict = self._append_meta_data(
                meta_data_path, n, events_total, total_rate
            )

            continue_condition = self._continue_condition_check(
                size,
                buffer_dict,
                stopping_criteria,
                initial_continue_condition=continue_condition,
            )

        print(f"stored detectable unlensed params in {output_path}")
        print(f"stored meta data in {meta_data_path}")

        if trim_to_size:
            param_final, total_rate = self._trim_results_to_size(
                size, output_path, meta_data_path
            )
        else:
            param_final = get_param_from_json(output_path)
            meta_data = get_param_from_json(meta_data_path)
            total_rate = meta_data["total_rate"][-1]

        # call self.json_file_names["ler_param"] and for adding the final results
        data = load_json(self.ler_directory + "/" + self.json_file_names["ler_params"])
        # write the results
        try:
            data["detectable_unlensed_rate_per_year"] = total_rate
            data["pdet_type_unlensed"] = pdet_type
        except:
            meta = get_param_from_json(meta_data_path)
            data["detectable_unlensed_rate_per_year"] = meta["total_rate"][-1]
            data["pdet_type_unlensed"] = pdet_type

        append_json(
            self.ler_directory + "/" + self.json_file_names["ler_params"],
            data,
            replace=True,
        )

        return total_rate,param_final

    def selecting_n_lensed_detectable_events(
        self,
        size=100,
        stopping_criteria=dict(
            relative_diff_percentage=2,
            number_of_last_batches_to_check=4,
        ),
        batch_size=50000,
        pdet_threshold=[0.5, 0.5],
        num_img=[1, 1],
        resume=True,
        pdet_type="boolean",
        output_jsonfile="n_lensed_params_detectable.json",
        meta_data_file="meta_lensed.json",
        trim_to_size=False,
        nan_to_num=True,
    ):
        """
        Generate a target number of detectable lensed events by sampling in batches, with the option to stop once the cumulative rate has stabilized.

        This function samples lensed parameters and saves only the detectable
        events in a JSON file. It also records metadata including the total
        number of events and the cumulative rate.

        Parameters
        ----------
        size : ``int``
            Target number of detectable samples to collect. \n
            default: 100
        stopping_criteria : ``dict`` or ``None``
            Criteria for stopping sample collection (but will not stop until n>size). \n
            Keys: \n
            - 'relative_diff_percentage': Maximum relative difference in rate (float) \n
            - 'number_of_last_batches_to_check': Number of batches for comparison (int) \n
            If None, stops when detectable events exceed size. \n
            default: dict(relative_diff_percentage=2, number_of_last_batches_to_check=4)
        batch_size : ``int``
            Batch size for sampling. \n
            default: 50000
        pdet_threshold : ``float`` or ``list``
            Detection probability threshold. \n
            default: [0.5, 0.5]
        num_img : ``list``
            Number of images corresponding to each pdet_threshold. \n
            default: [1, 1]
        resume : ``bool``
            If True, resumes from last saved batch. \n
            default: True
        pdet_type : ``str``
            Detectability condition type. \n
            Options: \n
            - 'boolean': Binary detection based on pdet_threshold \n
            - 'probability_distribution': Uses pdet values directly \n
            default: 'boolean'
        output_jsonfile : ``str``
            JSON file name for storing detectable parameters. \n
            default: 'n_lensed_params_detectable.json'
        meta_data_file : ``str``
            JSON file name for storing metadata. \n
            default: 'meta_lensed.json'
        trim_to_size : ``bool``
            If True, trims final result to exactly size events. \n
            default: False
        nan_to_num : ``bool``
            If True, NaN values will be converted to 0. \n
            default: False

        Returns
        -------
        param_final : ``dict``
            Dictionary of lensed GW source parameters of detectable events.

        Examples
        --------
        >>> from ler import LeR
        >>> ler = LeR()
        >>> lensed_param = ler.selecting_n_lensed_detectable_events(size=100)
        """

        if stopping_criteria is not None:
            print(
                f"stopping criteria set to when relative difference of total rate for the last {stopping_criteria['number_of_last_batches_to_check']} cumulative batches is less than {stopping_criteria['relative_diff_percentage']}%."
            )
            print(
                "sample collection will stop when the stopping criteria is met and number of detectable events exceeds the specified size."
            )
        else:
            print(
                "stopping criteria not set. sample collection will stop when number of detectable events exceeds the specified size."
            )

        # initial setup
        (
            n_cumulative,
            events_total,
            output_path,
            meta_data_path,
            buffer_file,
            continue_condition,
        ) = self._initial_setup_for_n_event_selection(
            size, meta_data_file, output_jsonfile, resume, stopping_criteria
        )

        # re-analyse the provided pdet_threshold and num_img
        pdet_threshold, num_img = self._check_pdet_threshold_lensed(
            pdet_threshold, num_img
        )

        # loop until n==size samples are collected
        while continue_condition:
            # disable print statements
            with contextlib.redirect_stdout(None):
                self.dict_buffer = None  # this is used to store the sampled lensed_param in batches when running the sampling_routine
                lensed_param = self.lensed_sampling_routine(
                    size=batch_size, output_jsonfile=buffer_file, resume=False
                )  # Dimensions are (size, n_max_images)

            total_events_in_this_iteration = len(lensed_param["zs"])

            # find index of detectable events
            pdet_hit, pdet_prod = self._find_detectable_index_lensed(
                lensed_param, pdet_threshold, num_img, pdet_type
            )

            pdet_processed = (
                pdet_hit if pdet_type == "boolean" else pdet_prod
            )

            # store all params in json file
            self._save_detectable_params(
                output_jsonfile,
                lensed_param,
                pdet_hit,
                key_file_name="n_lensed_detectable_events",
                nan_to_num=nan_to_num,
                verbose=False,
                replace_jsonfile=False,
            )

            n_cumulative += np.sum(pdet_processed)
            events_total += total_events_in_this_iteration
            total_rate = self.rate_function(
                n_cumulative, events_total, param_type="lensed", verbose=False
            )

            # save meta data
            buffer_dict = self._append_meta_data(
                meta_data_path, n_cumulative, events_total, total_rate
            )

            continue_condition = self._continue_condition_check(
                size,
                buffer_dict,
                stopping_criteria,
                initial_continue_condition=continue_condition,
            )

        print(f"storing detectable lensed params in {output_path}")
        print(f"storing meta data in {meta_data_path}")

        if trim_to_size:
            param_final, total_rate = self._trim_results_to_size(
                size, output_path, meta_data_path, param_type="lensed"
            )
        else:
            param_final = get_param_from_json(output_path)
            meta_data = get_param_from_json(meta_data_path)
            total_rate = meta_data["total_rate"][-1]

        # call self.json_file_names["ler_param"] and for adding the final results
        data = load_json(self.ler_directory + "/" + self.json_file_names["ler_params"])
        # write the results
        try:
            data["detectable_lensed_rate_per_year"] = total_rate
            data["pdet_type_lensed"] = pdet_type
        except:
            meta = get_param_from_json(meta_data_path)
            data["detectable_lensed_rate_per_year"] = meta["total_rate"][-1]
            data["pdet_type_lensed"] = pdet_type
        buffer_dict = append_json(
            self.ler_directory + "/" + self.json_file_names["ler_params"],
            data,
            replace=True,
        )

        return total_rate, param_final

    def _initial_setup_for_n_event_selection(
        self, size, meta_data_file, output_jsonfile, resume, stopping_criteria
    ):
        """Helper function for selecting_n_unlensed_detectable_events and selecting_n_lensed_detectable_events functions.

        Parameters
        ----------
        size : `int`
            number of samples to be collected.
        meta_data_file : ``str``
            json file name for storing the metadata.
        output_jsonfile : ``str``
            json file name for storing the parameters of the detectable events.
        resume : ``bool``
            resume = False (default) or True.
            if True, the function will resume from the last batch.
        batch_size : ``int``
            batch size for sampling.
            default batch_size = 50000.

        Returns
        ----------
        n : ``int``
            iterator.
        events_total : ``int``
            total number of events.
        output_path : ``str``
            path to the output json file.
        meta_data_path : ``str``
            path to the metadata json file.
        buffer_file : ``str``
            path to the buffer json file.
        """
        continue_condition = True
        meta_data_path = self.ler_directory + "/" + meta_data_file
        output_path = self.ler_directory + "/" + output_jsonfile
        if meta_data_path == output_path:
            raise ValueError("meta_data_file and output_jsonfile cannot be same.")

        if not resume:
            n = 0  # iterator
            events_total = 0
            remove_file(output_path)
            remove_file(meta_data_path)
        else:
            # get sample size as size from json file
            if os.path.exists(meta_data_path):
                param_final = get_param_from_json(output_path)
                n_collected = len(param_final["zs"])
                meta_data = get_param_from_json(meta_data_path)
                n = meta_data["detectable_events"][-1]
                events_total = meta_data["events_total"][-1]

                if not n_collected == n:
                    print(
                        "Number of collected events does not match with the number of events in the meta data file."
                    )

                else:
                    print(f"Resuming from {n} detectable events.")
                    # check if stopping criteria is met
                    continue_condition = self._continue_condition_check(
                        size,
                        meta_data,
                        stopping_criteria,
                        initial_continue_condition=continue_condition,
                    )

                if continue_condition is False:
                    print(
                        "Stopping criteria met. There will be no more samples collected."
                    )
            else:
                n = 0
                events_total = 0

        buffer_file = "params_buffer.json"
        print("collected number of detectable events = ", n)

        return (
            n,
            events_total,
            output_path,
            meta_data_path,
            buffer_file,
            continue_condition,
        )

    def _continue_condition_check(
        self,
        size_to_collect,
        param_dict,
        stopping_criteria,
        initial_continue_condition=True,
    ):

        continue_condition = initial_continue_condition
        already_collected_size = param_dict["detectable_events"][-1]
        # check if stopping criteria is met
        if isinstance(stopping_criteria, dict):
            total_rates = np.array(param_dict["total_rate"])
            limit = stopping_criteria["relative_diff_percentage"]
            num_a = stopping_criteria["number_of_last_batches_to_check"]

            if len(total_rates) > num_a:
                num_a = int(-1 * (num_a))
                percentage_diff = (
                    np.abs((total_rates[num_a:] - total_rates[-1]) / total_rates[-1])
                    * 100
                )
                print(
                    f"percentage difference of total rate for the last {abs(num_a)} cumulative batches = {percentage_diff}"
                )
                if np.any(percentage_diff > limit):
                    continue_condition &= True
                else:
                    print(
                        rf"stopping criteria of rate relative difference of {limit}% for the last {abs(num_a)} cumulative batches reached."
                    )
                    continue_condition &= False

        if isinstance(size_to_collect, int):
            if already_collected_size < size_to_collect:
                continue_condition |= True
            else:
                print(f"Given size={size_to_collect} reached\n")
                continue_condition |= False
                if stopping_criteria is None:
                    continue_condition &= False

        return continue_condition

    def _trim_results_to_size(
        self, size, output_path, meta_data_path, param_type="unlensed"
    ):
        """
        Helper function of 'selecting_n_unlensed_detectable_events' and 'selecting_n_lensed_detectable_events' functions. Trims the data in the output file to the specified size and updates the metadata accordingly.

        Parameters
        ----------
        size : `int`
            number of samples to be selected.
        output_path : ``str``
            path to the output json file.
        meta_data_path : ``str``
            path to the metadata json file.
        param_type : ``str``
            type of parameters.
            default param_type = "unlensed".
            other options are "lensed".

        Returns
        ----------
        param_final : ``dict``
            dictionary of GW source parameters of the detectable events. Refer to :meth:`~unlensed_param` or :meth:`~lensed_param` for details.
        new_total_rate : ``float``
            total rate (yr^-1).
        """

        print(f"\n trmming final result to size={size}")
        param_final = get_param_from_json(output_path)
        # randomly select size number of samples
        len_ = len(list(param_final.values())[0])
        idx = np.random.choice(len_, size, replace=False)
        # trim the final param dictionary, randomly, without repeating
        for key, value in param_final.items():
            param_final[key] = value[idx]

        # change meta data
        meta_data = load_json(meta_data_path)
        old_events_total = meta_data["events_total"][-1]
        old_detectable_events = meta_data["detectable_events"][-1]

        # adjust the meta data
        # following is to keep rate the same
        new_events_total = np.round(size * old_events_total / old_detectable_events)
        new_total_rate = self.rate_function(
            size, new_events_total, param_type=param_type, verbose=False
        )
        meta_data["events_total"][-1] = new_events_total
        meta_data["detectable_events"][-1] = size
        meta_data["total_rate"][-1] = new_total_rate

        print("collected number of detectable events = ", size)
        print("total number of events = ", new_events_total)
        print(f"total {param_type} event rate (yr^-1): {new_total_rate}")

        # save the meta data
        append_json(meta_data_path, meta_data, replace=True)
        # save the final param dictionary
        append_json(output_path, param_final, replace=True)

        return param_final, new_total_rate

    def _append_meta_data(self, meta_data_path, n, events_total, total_rate):
        """
        Helper function for appending meta data json file.

        Parameters
        ----------
        meta_data_path : ``str``
            path to the metadata json file.
        n : `int`
            iterator.
        events_total : `int`
            total number of events.
        total_rate : ``float``
            total rate (yr^-1).
        """

        # save meta data
        meta_data = dict(
            events_total=[events_total],
            detectable_events=[float(n)],
            total_rate=[total_rate],
        )

        if os.path.exists(meta_data_path):
            try:
                dict_ = append_json(meta_data_path, meta_data, replace=False)
            except:
                dict_ = append_json(meta_data_path, meta_data, replace=True)
        else:
            dict_ = append_json(meta_data_path, meta_data, replace=True)

        batch_n = (dict_["detectable_events"][-1]-dict_["detectable_events"][-2]) if len(dict_["detectable_events"]) > 1 else n

        print("collected number of detectable events (batch) = ", batch_n)
        print("collected number of detectable events (cumulative) = ", n)
        print("total number of events = ", events_total)
        print(f"total rate (yr^-1): {total_rate}")

        return dict_

    @property
    def create_new_interpolator(self):
        """
        Configuration dictionary for interpolator creation settings.

        Returns
        -------
        create_new_interpolator : ``dict``
            Dictionary specifying which interpolators to create. \n
            Each key is an interpolator name, and values are dicts with: \n
            - 'create_new': bool - Whether to create new interpolator \n
            - 'resolution': int or list - Grid resolution for interpolation \n
            Special key 'gwsnr' is a bool for GWSNR interpolator creation.
            Default: dict(
                merger_rate_density = {'create_new': False, 'resolution': 500},
                redshift_distribution = {'create_new': False, 'resolution': 500},
                luminosity_distance = {'create_new': False, 'resolution': 500},
                differential_comoving_volume = {'create_new': False, 'resolution': 500},
                source_frame_masses = {'create_new': False, 'resolution': 500},
                geocent_time = {'create_new': False, 'resolution': 500},
                ra = {'create_new': False, 'resolution': 500},
                dec = {'create_new': False, 'resolution': 500},
                phase = {'create_new': False, 'resolution': 500},
                psi = {'create_new': False, 'resolution': 500},
                theta_jn = {'create_new': False, 'resolution': 500},
                a_1 = {'create_new': False, 'resolution': 500},
                a_2 = {'create_new': False, 'resolution': 500},
                tilt_1 = {'create_new': False, 'resolution': 500},
                tilt_2 = {'create_new': False, 'resolution': 500},
                phi_12 = {'create_new': False, 'resolution': 500},
                phi_jl = {'create_new': False, 'resolution': 500},
                velocity_dispersion = {'create_new': False, 'resolution': 500, 'zl_resolution': 48},
                axis_ratio = {'create_new': False, 'resolution': 500, 'sigma_resolution': 48},
                lens_redshift = {'create_new': False, 'resolution': 48, 'zl_resolution': 48},
                lens_redshift_intrinsic = {'create_new': False, 'resolution': 500},
                optical_depth = {'create_new': False, 'resolution': 48},
                comoving_distance = {'create_new': False, 'resolution': 500},
                angular_diameter_distance = {'create_new': False, 'resolution': 500},
                angular_diameter_distance_z1z2 = {'create_new': False, 'resolution': 500},
                density_profile_slope = {'create_new': False, 'resolution': 100},
                lens_parameters_kde_sl = {'create_new': False, 'resolution': 5000},
                cross_section = {'create_new': False, 'resolution': [25, 25, 45, 15, 15]},
                gwsnr = False,
            )
        """
        return self._create_new_interpolator

    @create_new_interpolator.setter
    def create_new_interpolator(self, value):
        self._create_new_interpolator = value

    @property
    def npool(self):
        """
        Number of parallel processing cores.

        Returns
        -------
        npool : ``int``
            Number of logical cores to use for multiprocessing. \n
            default: 4
        """
        return self._npool

    @npool.setter
    def npool(self, value):
        self._npool = value

    @property
    def z_min(self):
        """
        Minimum redshift of the source population.

        Returns
        -------
        z_min : ``float``
            Minimum source redshift for sampling. \n
            default: 0.0
        """
        return self._z_min

    @z_min.setter
    def z_min(self, value):
        self._z_min = value

    @property
    def z_max(self):
        """
        Maximum redshift of the source population.

        Returns
        -------
        z_max : ``float``
            Maximum source redshift for sampling. \n
            default: 10.0
        """
        return self._z_max

    @z_max.setter
    def z_max(self, value):
        self._z_max = value

    @property
    def event_type(self):
        """
        Type of compact binary coalescence event.

        Returns
        -------
        event_type : ``str``
            Type of CBC event. \n
            Options: \n
            - 'BBH': Binary Black Hole \n
            - 'BNS': Binary Neutron Star \n
            - 'NSBH': Neutron Star-Black Hole \n
            default: 'BBH'
        """
        return self._event_type

    @event_type.setter
    def event_type(self, value):
        self._event_type = value

    @property
    def lens_type(self):
        """
        Type of lens galaxy model.

        Returns
        -------
        lens_type : ``str``
            Type of lens model. \n
            Options: \n
            - 'epl_shear_galaxy': Elliptical Power Law with external shear \n
            - 'sie_galaxy': Singular Isothermal Ellipsoid \n
            - 'sis_galaxy': Singular Isothermal Sphere \n
            default: 'epl_shear_galaxy'
        """
        return self._lens_type

    @lens_type.setter
    def lens_type(self, value):
        self._lens_type = value

    @property
    def cosmo(self):
        """
        Astropy cosmology object for distance calculations.

        Returns
        -------
        cosmo : ``astropy.cosmology``
            Cosmology used for luminosity distance and comoving volume calculations. \n
            default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        """
        return self._cosmo

    @cosmo.setter
    def cosmo(self, value):
        self._cosmo = value

    @property
    def json_file_names(self):
        """
        Dictionary of JSON file names for parameter storage.

        Returns
        -------
        json_file_names : ``dict``
            Dictionary with keys: \n
            - 'ler_params': LeR initialization parameters \n
            - 'unlensed_param': Unlensed event parameters \n
            - 'unlensed_param_detectable': Detectable unlensed events \n
            - 'lensed_param': Lensed event parameters \n
            - 'lensed_param_detectable': Detectable lensed events \n
        """
        return self._json_file_names

    @json_file_names.setter
    def json_file_names(self, value):
        self._json_file_names = value

    @property
    def interpolator_directory(self):
        """
        Directory path for interpolator JSON files.

        Returns
        -------
        interpolator_directory : ``str``
            Path to directory containing interpolator data files. \n
            default: './interpolator_json'
        """
        return self._interpolator_directory

    @interpolator_directory.setter
    def interpolator_directory(self, value):
        self._interpolator_directory = value

    @property
    def ler_directory(self):
        """
        Directory path for LeR output files.

        Returns
        -------
        ler_directory : ``str``
            Path to directory for storing output parameter files. \n
            default: './ler_data'
        """
        return self._ler_directory

    @ler_directory.setter
    def ler_directory(self, value):
        self._ler_directory = value

    @property
    def list_of_detectors(self):
        """
        List of gravitational wave detector names.

        Returns
        -------
        list_of_detectors : ``list``
            List of detector identifiers used for pdet calculations. \n
            Typically set from gwsnr.detector_list during initialization.
        """
        return self._list_of_detectors

    @list_of_detectors.setter
    def list_of_detectors(self, value):
        self._list_of_detectors = value

    @property
    def pdet_finder(self):
        """
        Detection probability finder function.

        Returns
        -------
        pdet_finder : ``callable``
            Function that calculates detection probability for GW events. \n
            The function signature should be: \n
            ``pdet_finder(gw_param_dict) -> dict`` with key 'pdet_net'.
        """
        return self._pdet_finder

    @pdet_finder.setter
    def pdet_finder(self, value):
        self._pdet_finder = value

    @property
    def ler_args(self):
        """
        Dictionary of all LeR initialization arguments.

        Returns
        -------
        ler_args : ``dict``
            Dictionary containing all parameters used to initialize LeR and \n
            its parent classes, useful for reproducibility.
        """
        return self._ler_args

    @ler_args.setter
    def ler_args(self, value):
        self._ler_args = value
