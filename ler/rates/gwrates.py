# -*- coding: utf-8 -*-
"""
Module for calculating detection rates of gravitational wave events.

This module contains the main ``GWRATES`` class for calculating the rates of
detectable gravitational wave events. The class inherits from
:class:`~ler.gw_source_population.CBCSourceParameterDistribution` for source
parameters sampling and utilizes the ``gwsnr`` package for detection probability
calculation.

Inheritance hierarchy:

- :class:`~ler.gw_source_population.CBCSourceParameterDistribution` \n
  - :class:`~ler.gw_source_population.CBCSourceRedshiftDistribution` \n
- Uses the ``gwsnr`` package for pdet calculation.

Usage:
    Basic workflow for rate calculation:

    >>> from ler.rates import GWRATES
    >>> gwrates = GWRATES()
    >>> gw_params = gwrates.gw_cbc_statistics()
    >>> gwrates.gw_rate()

Copyright (C) 2026 Phurailatpam Hemantakumar. Distributed under MIT License.
"""

import os
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
from ..gw_source_population import CBCSourceParameterDistribution
from ..utils import (
    load_json,
    append_json,
    get_param_from_json,
    batch_handler,
    remove_file,
)


class GWRATES(CBCSourceParameterDistribution):
    """
    Class to sample GW events and calculate their detection rates.

    This class provides functionality for sampling gravitational wave source
    parameters, detection probabilities, and computing detection rates for
    compact binary coalescence events. Parameters of simulated events are
    stored in JSON files (not as class attributes) to conserve RAM memory.

    Key Features: \n
    - Sampling of GW event parameters \n
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
        Type of event to generate. \n
        Options: \n
        - 'BBH': Binary Black Hole \n
        - 'BNS': Binary Neutron Star \n
        - 'NSBH': Neutron Star-Black Hole \n
        default: 'BBH'
    cosmology : ``astropy.cosmology`` or ``None``
        Cosmology to use for the calculation. \n
        default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
    pdet_finder : ``callable`` or ``None``
        Custom detection probability finder function. \n
        If None, uses gwsnr's pdet calculator. \n
        The function should follow the signature: \n
        ``def pdet_finder(gw_param_dict): return pdet_net_dict`` \n
        where pdet_net_dict.keys = ['pdet_net']. \n
        default: None
    json_file_names : ``dict`` or ``None``
        Names of the JSON files to store the necessary parameters. \n
        default: dict(gwrates_params="gwrates_params.json", gw_param="gw_param.json", gw_param_detectable="gw_param_detectable.json")
    interpolator_directory : ``str``
        Directory to store the interpolators. \n
        default: './interpolator_json'
    create_new_interpolator : ``bool`` or ``dict``
        Whether to create new interpolators. \n
        default: False
    ler_directory : ``str``
        Directory to store the output parameters. \n
        default: './ler_data'
    verbose : ``bool``
        If True, print all chosen parameters during initialization. \n
        default: True
    **kwargs : ``dict``
        Additional keyword arguments passed to parent classes: \n
        :class:`~ler.gw_source_population.CBCSourceParameterDistribution` and \n
        :class:`~gwsnr.GWSNR` (if pdet_finder is not provided).

    Examples
    --------
    Basic usage:

    >>> from ler.rates import GWRATES
    >>> gwrates = GWRATES()
    >>> gw_params = gwrates.gw_cbc_statistics()
    >>> gwrates.gw_rate()

    Instance Methods
    ----------
    GWRATES class has the following methods: \n
    +-----------------------------------------------------+------------------------------------------------+
    | Method                                              | Description                                    |
    +=====================================================+================================================+
    | :meth:`~gw_cbc_statistics`                          | Generate GW source parameters                  |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~gw_rate`                                    | Calculate the GW detection rate                |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~selecting_n_gw_detectable_events`           | Select n GW detectable events                  |
    +-----------------------------------------------------+------------------------------------------------+

    Instance Attributes
    ----------
    GWRATES class has the following attributes: \n
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | Attribute                                      | Type             | Unit  | Description                                    |
    +================================================+==================+=======+================================================+
    | :meth:`~npool`                                 | ``int``          |       | Number of parallel processing cores            |
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | :meth:`~z_min`                                 | ``float``        |       | Minimum source redshift                        |
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | :meth:`~z_max`                                 | ``float``        |       | Maximum source redshift                        |
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | :meth:`~event_type`                            | ``str``          |       | Type of CBC event                              |
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | :meth:`~cosmo`                                 | ``astropy.cosmology`` |  | Cosmology for calculations                     |
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
    | :meth:`~gwrates_args`                          | ``dict``         |       | All GWRATES initialization arguments           |
    +------------------------------------------------+------------------+-------+------------------------------------------------+

    Notes
    -----
    - ``GWRATES`` class inherits from :class:`~ler.gw_source_population.CBCSourceParameterDistribution`. \n
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
            zip_resource = resources_files("ler.rates").joinpath(
                "ler_data", "interpolator_json.zip"
            )
            with zip_resource.open("rb") as zip_file:
                print(
                    "Extracting interpolator data from package to the current working directory."
                )

                # Define destination path (current working directory)
                dest_path = pathlib.Path.cwd()

                # Extract the zip file, skipping __MACOSX metadata
                with zipfile.ZipFile(zip_file, "r") as zip_ref:
                    for member in zip_ref.namelist():
                        # Skip __MACOSX directory and its contents
                        if member.startswith("__MACOSX"):
                            continue
                        zip_ref.extract(member, dest_path)

        print("\nInitializing GWRATES class...\n")
        # init gwrates attributes
        self.npool = npool
        self.z_min = z_min
        self.z_max = z_max
        self.event_type = event_type
        self.cosmo = cosmology if cosmology else LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

        # init json file names where datas will be stored
        self.json_file_names = dict(
            gwrates_params="gwrates_params.json",
            gw_param="gw_param.json",
            gw_param_detectable="gw_param_detectable.json",
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
        Helper function to initialize the parent classes.

        Parameters
        ----------
        params : ``dict`` or ``None``
            Dictionary of parameters to initialize the parent classes.
        pdet_finder : ``callable`` or ``None``
            Custom detection probability finder function.
        verbose : ``bool``
            If True, print initialization parameters. \n
            default: True
        """

        def initialization():
            # initialization of parent class
            self._parent_initialization_helper(params=params)
            # initialization self.pdet_finder from GWSNR class
            if not pdet_finder:
                self.pdet_finder = self._gwsnr_initialization(params=params)
            else:
                self.pdet_finder = pdet_finder

            # store all the gwrates input parameters
            self._store_gwrates_params(
                output_jsonfile=self.json_file_names["gwrates_params"]
            )

        # if not verbose, prevent anything from printing
        if verbose:
            initialization()
            self._print_all_init_args()
        else:
            with contextlib.redirect_stdout(None):
                initialization()

    def _parent_initialization_helper(self, params=None):
        """
        Helper function to initialize CBCSourceParameterDistribution parent class.

        Parameters
        ----------
        params : ``dict`` or ``None``
            Dictionary of parameters to initialize the parent classes.
        """

        # initialization of CBCSourceParameterDistribution class
        # it also initializes the CBCSourceRedshiftDistribution class
        input_params = dict(
            npool=self.npool,
            z_min=self.z_min,
            z_max=self.z_max,
            cosmology=self.cosmo,
            event_type=self.event_type,
            source_priors=None,
            source_priors_params=None,
            spin_zero=False,
            spin_precession=False,
            directory=self.interpolator_directory,
            create_new_interpolator=False,
        )
        # update input_params with params. This will include create_new_interpolator.
        if params:
            for key, value in params.items():
                if key in input_params:
                    input_params[key] = value
        # initialization of parent class
        CBCSourceParameterDistribution.__init__(
            self,
            z_min=input_params["z_min"],
            z_max=input_params["z_max"],
            cosmology=input_params["cosmology"],
            event_type=input_params["event_type"],
            source_priors=input_params["source_priors"],
            source_priors_params=input_params["source_priors_params"],
            spin_zero=input_params["spin_zero"],
            spin_precession=input_params["spin_precession"],
            directory=input_params["directory"],
            create_new_interpolator=input_params["create_new_interpolator"],
        )

        # save input_params to self.gwrates_args
        # some of the None values will have default values after initialization
        input_params["source_priors"] = self.gw_param_samplers.copy()
        input_params["source_priors_params"] = self.gw_param_samplers_params.copy()
        input_params["create_new_interpolator"] = self.create_new_interpolator
        self.gwrates_args = input_params

    def _gwsnr_initialization(self, params=None):
        """
        Helper function to initialize the GWSNR class from the ``gwsnr`` package.

        Parameters
        ----------
        params : ``dict`` or ``None``
            Dictionary of parameters to initialize the GWSNR class.

        Returns
        -------
        pdet : ``callable``
            Detection probability function from GWSNR.
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
        self.gwrates_args["pdet_args"] = input_params

        # dealing with create_new_interpolator param
        if isinstance(input_params["create_new_interpolator"], dict):
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

        self.gwrates_args["list_of_detectors"] = gwsnr.detector_list
        self.list_of_detectors = self.gwrates_args["list_of_detectors"]

        self.gwrates_args["pdet_args"]["pdet_kwargs"] = gwsnr.pdet_kwargs
        self.gwrates_args["pdet_args"]["psds_list"] = gwsnr.psds_list

        return gwsnr.pdet

    def _print_all_init_args(self):
        """
        Helper function to print all initialization parameters.
        """

        # print all relevant functions and sampler priors
        print("#----------------------------------------")
        print("# GWRATES initialization input arguments:")
        print("#----------------------------------------\n")
        print("    # GWRATES set up input arguments:")
        print(f"    npool = {self.npool},")
        print(f"    z_min = {self.z_min},")
        print(f"    z_max = {self.z_max},")
        print(f"    event_type = '{self.event_type}',")
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

        print(
            "\n    # GWRATES also takes other CBCSourceParameterDistribution class input arguments as kwargs, as follows:"
        )
        print("    source_priors = dict(")
        for key, value in self.gwrates_args["source_priors"].items():
            (
                print(f"        {key} = '{value}',")
                if isinstance(value, str)
                else print(f"        {key} = {value},")
            )
        print("    ),")
        print("    source_priors_params = dict(")
        for key, value in self.gwrates_args["source_priors_params"].items():
            (
                print(f"        {key} = '{value}',")
                if isinstance(value, str)
                else print(f"        {key} = {value},")
            )
        print("    ),")
        print(f"    spin_zero = {self.gwrates_args['spin_zero']},")
        print(f"    spin_precession = {self.gwrates_args['spin_precession']},")

        if "pdet_args" in self.gwrates_args:
            print(
                "\n    # GWRATES also takes other gwsnr.GWSNR input arguments as kwargs, as follows:"
            )
            for key, value in self.gwrates_args["pdet_args"].items():
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

    def _store_gwrates_params(self, output_jsonfile=None):
        """
        Helper function to store all initialization parameters to JSON file.

        This is useful for reproducing results. All parameters are stored
        in string format for JSON compatibility.

        Parameters
        ----------
        output_jsonfile : ``str`` or ``None``
            Name of the JSON file to store the parameters. \n
            default: None (uses self.json_file_names["gwrates_params"])
        """

        # gwrates params
        param_sampler_dict = self.gwrates_args.copy()
        # convert all dict values to str
        for key, value in param_sampler_dict.items():
            param_sampler_dict[key] = str(value)

        file_name = (
            output_jsonfile
            if output_jsonfile
            else self.json_file_names["gwrates_params"]
        )
        append_json(
            os.path.join(self.ler_directory, file_name),
            param_sampler_dict,
            replace=True,
        )

    def gw_cbc_statistics(
        self,
        size=100000,
        batch_size=50000,
        resume=True,
        save_batch=False,
        output_jsonfile=None,
    ):
        """
        Generate GW source parameters with detection probabilities.

        This function calls the _gw_sampling_routine function to generate
        the parameters in batches. The generated parameters are stored in
        a JSON file; if save_batch=True, it updates the file after each batch.

        Parameters
        ----------
        size : ``int``
            Number of samples to generate. \n
            default: 100000
        batch_size : ``int``
            Batch size for sampling. \n
            default: 50000
        resume : ``bool``
            If True, resume from the last batch. \n
            default: True
        save_batch : ``bool``
            If True, save parameters after each batch. \n
            If False, save all parameters at the end (faster). \n
            default: False
        output_jsonfile : ``str`` or ``None``
            JSON file name for storing the parameters. \n
            default: 'gw_param.json' (stored in ler_directory)

        Returns
        -------
        gw_param : ``dict``
            dictionary of GW source parameters. The included parameters and their units are as follows (for default settings):\n
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
            | a_1                |              | spin of the primary compact binary         |
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
        >>> from ler.rates import GWRATES
        >>> gwrates = GWRATES()
        >>> param = gwrates.gw_cbc_statistics()
        """

        output_jsonfile = output_jsonfile or self.json_file_names["gw_param"]
        self.json_file_names["gw_param"] = output_jsonfile
        output_path = os.path.join(self.ler_directory, output_jsonfile)
        print(f"Simulated GW params will be stored in {output_path}")

        gw_param = batch_handler(
            size=size,
            batch_size=batch_size,
            sampling_routine=self._gw_sampling_routine,
            output_jsonfile=output_path,
            save_batch=save_batch,
            resume=resume,
            param_name="gw parameters",
        )

        return gw_param

    def _gw_sampling_routine(
        self, size, output_jsonfile, resume=False, save_batch=True
    ):
        """
        Private helper function to generate GW source parameters for a single batch.

        Parameters
        ----------
        size : ``int``
            Number of samples to generate.
        output_jsonfile : ``str``
            JSON file name for storing the parameters.
        resume : ``bool``
            If True, appends new samples to existing file. \n
            default: False
        save_batch : ``bool``
            If True, save parameters in batches. \n
            default: True

        Returns
        -------
        gw_param : ``dict``
            Dictionary of GW source parameters with detection probabilities.
        """
        print("sampling gw source params...")
        gw_param = self.sample_gw_parameters(size=size)

        print("calculating pdet...")
        pdet = self.pdet_finder(gw_param_dict=gw_param.copy())
        gw_param.update(pdet)

        return gw_param

    def gw_rate(
        self,
        gw_param=None,
        pdet_threshold=0.5,
        pdet_type="boolean",
        output_jsonfile=None,
    ):
        """
        Calculate the GW detection rate.

        This function calculates the detection rate and stores the parameters
        of detectable events in a JSON file.

        Parameters
        ----------
        gw_param : ``dict`` or ``str`` or ``None``
            Dictionary of GW source parameters or JSON file name. \n
            default: None (uses self.json_file_names["gw_param"])
        pdet_threshold : ``float``
            Threshold for detection probability. \n
            default: 0.5
        pdet_type : ``str``
            Detectability condition type. \n
            Options: \n
            - 'boolean': Binary detection based on pdet_threshold \n
            - 'probability_distribution': Uses pdet values directly \n
            default: 'boolean'
        output_jsonfile : ``str`` or ``None``
            JSON file name for storing detectable event parameters. \n
            default: 'gw_param_detectable.json'

        Returns
        -------
        total_rate : ``float``
            Total GW detection rate (yr^-1).
        gw_param : ``dict``
            dictionary of GW source parameters of the detectable events. The included parameters and their units are as follows (for default settings):\n
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
            | a_1                |              | spin of the primary compact binary         |
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
        >>> from ler.rates import GWRATES
        >>> gwrates = GWRATES()
        >>> gwrates.gw_cbc_statistics()
        >>> total_rate, gw_param = gwrates.gw_rate()
        """

        gw_param = self._load_param(gw_param, param_type="gw")
        total_events = len(gw_param["zs"])

        # find index of detectable events based on pdet
        pdet_net = gw_param["pdet_net"]
        idx_detectable = pdet_net >= pdet_threshold

        if pdet_type == "boolean":
            detectable_events = np.sum(idx_detectable)
        elif pdet_type == "probability_distribution":
            detectable_events = np.sum(pdet_net)
        else:
            raise ValueError("pdet_type not recognized")
        # Monte Carlo integration: R = norm * <Theta(pdet - pdet_threshold)>
        total_rate = self.rate_function(detectable_events, total_events)

        # store all detectable params in json file
        self._save_detectable_params(
            output_jsonfile,
            gw_param,
            idx_detectable,
            key_file_name="gw_param_detectable",
            nan_to_num=False,
            verbose=True,
            replace_jsonfile=True,
        )

        # append gwrates_param and save it
        self._append_gwrates_param(total_rate)

        return total_rate, gw_param

    def _load_param(self, param, param_type="gw"):
        """
        Helper function to load or copy GW parameters.

        Parameters
        ----------
        param : ``dict`` or ``str`` or ``None``
            Dictionary of GW parameters or JSON file name.
        param_type : ``str``
            Type of parameter for logging. \n
            default: "gw"

        Returns
        -------
        param : ``dict``
            Dictionary of GW parameters.
        """

        if param is None:
            param = self.json_file_names["gw_param"]
        if isinstance(param, str):
            path_ = os.path.join(self.ler_directory, param)
            print(f"Getting {param_type} parameters from json file {path_}...")
            return get_param_from_json(path_)
        else:
            print(f"Using provided {param_type} dict...")
            return param.copy()

    def rate_function(self, detectable_size, total_size, verbose=True):
        """
        Helper function to calculate the detection rate via Monte Carlo integration.

        Parameters
        ----------
        detectable_size : ``int`` or ``float``
            Number of detectable events.
        total_size : ``int``
            Total number of simulated events.
        verbose : ``bool``
            If True, print rate information. \n
            default: True

        Returns
        -------
        rate : ``float``
            Detection rate (yr^-1).
        """
        normalization = self.normalization_pdf_z
        rate = float(normalization * detectable_size / total_size)

        if verbose:
            print(f"total GW event rate (yr^-1): {rate}")
            print(f"number of simulated GW detectable events: {detectable_size}")
            print(f"number of simulated all GW events: {total_size}")

        return rate

    def _save_detectable_params(
        self,
        output_jsonfile,
        param,
        idx_detectable,
        key_file_name="gw_param_detectable",
        nan_to_num=False,
        verbose=True,
        replace_jsonfile=True,
    ):
        """
        Helper function to save detectable event parameters to JSON file.

        Parameters
        ----------
        output_jsonfile : ``str`` or ``None``
            JSON file name for storing detectable event parameters.
        param : ``dict``
            Dictionary of GW source parameters.
        idx_detectable : ``numpy.ndarray``
            Boolean array indicating detectable events.
        key_file_name : ``str``
            Key name for self.json_file_names dictionary. \n
            default: "gw_param_detectable"
        nan_to_num : ``bool``
            If True, replace NaN values with 0. \n
            default: False
        verbose : ``bool``
            If True, print the output file path. \n
            default: True
        replace_jsonfile : ``bool``
            If True, replace existing file. If False, append to file. \n
            default: True
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

        output_path = os.path.join(self.ler_directory, output_jsonfile)
        if verbose:
            print(f"storing detectable params in {output_path}")
        append_json(output_path, param, replace=replace_jsonfile)

    def _append_gwrates_param(self, total_rate):
        """
        Helper function to append the detection rate to the gwrates params JSON file.

        Parameters
        ----------
        total_rate : ``float``
            Total detection rate (yr^-1).
        """

        gwrates_params_path = os.path.join(
            self.ler_directory, self.json_file_names["gwrates_params"]
        )
        data = load_json(gwrates_params_path)
        # write the results
        data["detectable_gw_rate_per_year"] = total_rate
        append_json(gwrates_params_path, data, replace=True)

    def selecting_n_gw_detectable_events(
        self,
        size=100,
        batch_size=50000,
        stopping_criteria=dict(
            relative_diff_percentage=0.5, number_of_last_batches_to_check=4
        ),
        pdet_threshold=0.5,
        resume=True,
        output_jsonfile="gw_params_n_detectable.json",
        meta_data_file="meta_gw.json",
        pdet_type="boolean",
        trim_to_size=False,
    ):
        """
        Generate a target number of detectable GW events by iterative batch sampling.

        This function samples GW parameters in batches and saves only the
        detectable events to a JSON file. It optionally stops when the
        cumulative rate has stabilized based on the stopping criteria.

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
            Threshold for detection probability. \n
            default: 0.5
        resume : ``bool``
            If True, resume from the last batch. \n
            default: True
        output_jsonfile : ``str``
            JSON file name for storing detectable event parameters. \n
            default: 'gw_params_n_detectable.json'
        meta_data_file : ``str``
            JSON file name for storing metadata. \n
            default: 'meta_gw.json'
        pdet_type : ``str``
            Detectability condition type. \n
            Options: \n
            - 'boolean': Binary detection based on pdet_threshold \n
            - 'probability_distribution': Uses pdet values directly \n
            default: 'boolean'
        trim_to_size : ``bool``
            If True, trim final result to exactly the specified size. \n
            default: False

        Returns
        -------
        param_final : ``dict``
            dictionary of GW source parameters of the detectable events. The included parameters and their units are as follows (for default settings):\n
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
            | a_1                |              | spin of the primary compact binary   |
            +--------------------+--------------+--------------------------------------+
            | a_2                |              | spin of the secondary compact binary |
            +--------------------+--------------+--------------------------------------+
            | luminosity_distance| Mpc          | luminosity distance                  |
            +--------------------+--------------+--------------------------------------+
            | mass_1_source      | Msun         | mass of the primary compact binary   |
            |                    |              | (source frame)                       |
            +--------------------+--------------+--------------------------------------+
            | mass_2_source      | Msun         | mass of the secondary compact binary |
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
        >>> from ler.rates import GWRATES
        >>> gwrates = GWRATES()
        >>> gw_param = gwrates.selecting_n_gw_detectable_events(size=100)
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
                self.dict_buffer = None  # this is used to store the sampled GW_param in batches when running the sampling_routine
                gw_param = self._gw_sampling_routine(
                    size=batch_size,
                    output_jsonfile=buffer_file,
                    save_batch=False,
                    resume=False,
                )

            total_events_in_this_iteration = len(gw_param["zs"])
            # below is use when the snr is calculated with 'ann' method of `gwsnr`

            # find index of detectable events
            pdet = gw_param["pdet_net"]
            idx_detectable = pdet > pdet_threshold

            # store all params in json file
            self._save_detectable_params(
                output_jsonfile,
                gw_param,
                idx_detectable,
                key_file_name="n_gw_detectable_events",
                nan_to_num=False,
                verbose=False,
                replace_jsonfile=False,
            )

            if pdet_type == "boolean":
                n += np.sum(idx_detectable)
            else:
                n += np.sum(pdet)
            events_total += total_events_in_this_iteration
            total_rate = self.rate_function(n, events_total, verbose=False)

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

        print(f"stored detectable GW params in {output_path}")
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
        data = load_json(
            self.ler_directory + "/" + self.json_file_names["gwrates_params"]
        )
        # write the results
        try:
            data["detectable_gw_rate_per_year"] = total_rate
            data["pdet_type"] = pdet_type
        except:
            meta = get_param_from_json(meta_data_path)
            data["detectable_gw_rate_per_year"] = meta["total_rate"][-1]
            data["pdet_type"] = pdet_type

        append_json(
            self.ler_directory + "/" + self.json_file_names["gwrates_params"],
            data,
            replace=True,
        )

        return total_rate, param_final

    def _initial_setup_for_n_event_selection(
        self, size, meta_data_file, output_jsonfile, resume, stopping_criteria
    ):
        """
        Helper function for selecting_n_gw_detectable_events initialization.

        Parameters
        ----------
        size : ``int``
            Target number of samples to collect.
        meta_data_file : ``str``
            JSON file name for storing metadata.
        output_jsonfile : ``str``
            JSON file name for storing detectable event parameters.
        resume : ``bool``
            If True, resume from the last batch.
        stopping_criteria : ``dict`` or ``None``
            Stopping criteria for sample collection.

        Returns
        -------
        n : ``int``
            Current count of detectable events.
        events_total : ``int``
            Total number of simulated events.
        output_path : ``str``
            Path to the output JSON file.
        meta_data_path : ``str``
            Path to the metadata JSON file.
        buffer_file : ``str``
            Path to the buffer JSON file.
        continue_condition : ``bool``
            Whether to continue sampling.
        """

        continue_condition = True
        meta_data_path = os.path.join(self.ler_directory, meta_data_file)
        output_path = os.path.join(self.ler_directory, output_jsonfile)
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
        """
        Helper function to check if sampling should continue.

        Parameters
        ----------
        size_to_collect : ``int``
            Target number of detectable events.
        param_dict : ``dict``
            Dictionary containing metadata (detectable_events, total_rate).
        stopping_criteria : ``dict`` or ``None``
            Stopping criteria configuration.
        initial_continue_condition : ``bool``
            Initial value for continue condition. \n
            default: True

        Returns
        -------
        continue_condition : ``bool``
            Whether to continue sampling.
        """

        continue_condition = initial_continue_condition
        already_collected_size = param_dict["detectable_events"][-1]
        stopping_criteria_met = False
        size_reached = False

        # check if stopping criteria is met
        if isinstance(stopping_criteria, dict):
            total_rates = np.array(param_dict["total_rate"])
            limit = stopping_criteria["relative_diff_percentage"]
            num_a = stopping_criteria["number_of_last_batches_to_check"]

            if len(total_rates) > num_a:
                num_a_neg = int(-1 * num_a)
                percentage_diff = (
                    np.abs(
                        (total_rates[num_a_neg:] - total_rates[-1]) / total_rates[-1]
                    )
                    * 100
                )
                print(
                    f"percentage difference of total rate for the last {num_a} cumulative batches = {percentage_diff}"
                )
                if not np.any(percentage_diff > limit):
                    print(
                        rf"stopping criteria of rate relative difference of {limit}% for the last {num_a} cumulative batches reached."
                    )
                    stopping_criteria_met = True

        if isinstance(size_to_collect, int):
            if already_collected_size >= size_to_collect:
                print(f"Given size={size_to_collect} reached\n")
                size_reached = True

        # Determine continue condition
        if stopping_criteria is None:
            # Stop only when size is reached
            continue_condition = not size_reached
        else:
            # Stop when both stopping criteria is met AND size is reached
            continue_condition = not (stopping_criteria_met and size_reached)

        return continue_condition

    def _trim_results_to_size(self, size, output_path, meta_data_path):
        """
        Helper function to trim results to the specified size.

        Parameters
        ----------
        size : ``int``
            Target number of samples.
        output_path : ``str``
            Path to the output JSON file.
        meta_data_path : ``str``
            Path to the metadata JSON file.

        Returns
        -------
        param_final : ``dict``
            Dictionary of trimmed GW source parameters.
        new_total_rate : ``float``
            Updated total rate (yr^-1).
        """

        print(f"\n trimming final result to size={size}")
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
        new_total_rate = self.rate_function(size, new_events_total, verbose=False)
        meta_data["events_total"][-1] = new_events_total
        meta_data["detectable_events"][-1] = size
        meta_data["total_rate"][-1] = new_total_rate

        print("collected number of detectable events = ", size)
        print("total number of events = ", new_events_total)
        print(f"total GW event rate (yr^-1): {new_total_rate}")

        # save the meta data
        append_json(meta_data_path, meta_data, replace=True)
        # save the final param dictionary
        append_json(output_path, param_final, replace=True)

        return param_final, new_total_rate

    def _append_meta_data(self, meta_data_path, n, events_total, total_rate):
        """
        Helper function for appending metadata to JSON file.

        Parameters
        ----------
        meta_data_path : ``str``
            Path to the metadata JSON file.
        n : ``int`` or ``float``
            Number of detectable events collected.
        events_total : ``int``
            Total number of simulated events.
        total_rate : ``float``
            Total detection rate (yr^-1).

        Returns
        -------
        dict_ : ``dict``
            Updated metadata dictionary.
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

        print("collected number of detectable events = ", n)
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
            - 'gwrates_params': GWRATES initialization parameters \n
            - 'gw_param': GW event parameters \n
            - 'gw_param_detectable': Detectable GW events \n
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
    def gwrates_args(self):
        """
        Dictionary of all GWRATES initialization arguments.

        Returns
        -------
        gwrates_args : ``dict``
            Dictionary containing all parameters used to initialize GWRATES and \n
            its parent classes, useful for reproducibility.
        """
        return self._gwrates_args

    @gwrates_args.setter
    def gwrates_args(self, value):
        self._gwrates_args = value
