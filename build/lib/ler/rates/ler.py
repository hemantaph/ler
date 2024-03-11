# -*- coding: utf-8 -*-
"""
This module contains the main class for calculating the rates of detectable gravitational waves events. The class inherits the :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution` class for source parameters and lens parameters sampling. It also finds the image properties. :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution` inherits the :class:`~ler.gw_source_population.CBCSourceParameterDistribution`, :class:`~ler.image_properties.ImageProperties` and uses `gwsnr` package for SNR calculation. 
"""

import os
import warnings
warnings.filterwarnings("ignore")
import contextlib
import numpy as np
from scipy.stats import norm
from astropy.cosmology import LambdaCDM
from ..lens_galaxy_population import LensGalaxyParameterDistribution
from ..utils import load_json, append_json, get_param_from_json, batch_handler


class LeR(LensGalaxyParameterDistribution):
    """Class to calculate both the rates of lensed and unlensed events. Please note that parameters of the simulated events are stored in json file but not as an attribute of the class. This saves RAM memory. 

    Parameters
    ----------
    npool : `int`
        number of cores to use.
        default npool = 4.
    z_min : `float`
        minimum redshift.
        default z_min = 0.
        for popI_II, popIII, primordial, BNS z_min = 0., 5., 5., 0. respectively.
    z_max : `float`
        maximum redshift.
        default z_max = 10.
        for popI_II, popIII, primordial, BNS z_max = 10., 40., 40., 2. respectively.
    size : `int`
        number of samples for sampling.
        default size = 100000.
    batch_size : `int`
        batch size for SNR calculation.
        default batch_size = 25000.
        reduce the batch size if you are getting memory error.
        recommended batch_size = 50000, if size = 1000000.
    snr_finder : `str`
        default snr_finder = 'gwsnr'.
        if 'gwsnr', the SNR will be calculated using the gwsnr package.
        if 'custom', the SNR will be calculated using a custom function.
        The custom function should have input and output as given in GWSNR.snr method.
    json_file_names: `dict`
        names of the json files to strore the necessary parameters.
        default json_file_names = {'ler_param': './LeR_params.json', 'unlensed_param': './unlensed_param.json', 'unlensed_param_detectable': './unlensed_param_detectable.json'}.\n
    kwargs : `keyword arguments`
        Note : kwargs takes input for initializing the :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution`, :meth:`~gwsnr_intialization`.

    Examples
    ----------
    >>> from ler.rates import LeR
    >>> ler = LeR()
    >>> ler.unlensed_cbc_statistics();
    >>> ler.unlensed_rate();
        
    Instance Attributes
    ----------
    LeR class has the following attributes, \n
    +-------------------------------------+----------------------------------+
    | Atrributes                          | Type                             |
    +=====================================+==================================+
    |:attr:`~z_min`                       | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`~z_max`                       | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`~event_type`                  | `str`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~cosmo`                       | `astropy.cosmology`              |
    +-------------------------------------+----------------------------------+
    |:attr:`~size`                        | `int`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~batch_size`                  | `int`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~json_file_names`             | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~directory`                   | `str`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~gw_param_sampler_dict`       | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~snr_calculator_dict`         | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~list_of_detectors`           | `list`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~unlensed_param`              | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~unlensed_param_detectable`   | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~lensed_param`                | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~lensed_param_detectable`     | `dict`                           |
    +-------------------------------------+----------------------------------+

    Instance Methods
    ----------
    LeR class has the following methods, \n
    +-------------------------------------+----------------------------------+
    | Methods                             | Description                      |
    +=====================================+==================================+
    |:meth:`~class_initialization`        | Function to initialize the       |
    |                                     | parent classes                   |
    +-------------------------------------+----------------------------------+
    |:meth:`~gwsnr_intialization`         | Function to initialize the       |
    |                                     | gwsnr class                      |
    +-------------------------------------+----------------------------------+
    |:meth:`~snr`                         | Function to get the snr with the |
    |                                     | given parameters.                |
    +-------------------------------------+----------------------------------+
    |:meth:`~store_ler_params`            | Function to store the all the    |
    |                                     | necessary parameters.            |
    +-------------------------------------+----------------------------------+
    |:meth:`~unlensed_cbc_statistics`     | Function to generate unlensed    |
    |                                     | GW source parameters.            |
    +-------------------------------------+----------------------------------+
    |:meth:`~unlensed_sampling_routine`   | Function to generate unlensed    |
    |                                     | GW source parameters.            |
    +-------------------------------------+----------------------------------+
    |:meth:`~unlensed_rate`               | Function to calculate the        |
    |                                     | unlensed rate.                   |
    +-------------------------------------+----------------------------------+
    |:meth:`~selecting_n_unlensed_detectable_events`                         |
    +-------------------------------------+----------------------------------+
    |                                     | Function to select n unlensed    |
    |                                     | detectable events.               |
    +-------------------------------------+----------------------------------+
    |:meth:`~lensed_cbc_statistics`       | Function to generate lensed      |
    |                                     | GW source parameters.            |
    +-------------------------------------+----------------------------------+
    |:meth:`~lensed_sampling_routine`     | Function to generate lensed      |
    |                                     | GW source parameters.            |
    +-------------------------------------+----------------------------------+
    |:meth:`~lensed_rate`                 | Function to calculate the        |
    |                                     | lensed rate.                     |
    +-------------------------------------+----------------------------------+
    |:meth:`~rate_ratio`                  | Function to calculate the rate   |
    |                                     | ratio.                           |
    +-------------------------------------+----------------------------------+
    |:meth:`~rate_comparision_with_rate_calculation                          |
    +-------------------------------------+----------------------------------+
    |                                     | Function to compare the rates    |
    |                                     | calculated using LeR between     |
    |                                     | unlensed and lensed events.      |
    +-------------------------------------+----------------------------------+
    |:meth:`~param_plot`                  | Function to plot the             |
    |                                     | distribution of various          |
    |                                     | parameters.                      |
    +-------------------------------------+----------------------------------+
    |:meth:`~relative_mu_dt_lensed`       | Function to calculate the        |
    |                                     | relative magnification and       |
    |                                     | relative time-delay of lensed    |
    |                                     | events.                          |
    +-------------------------------------+----------------------------------+
    |:meth:`~relative_mu_dt_unlensed`     | Function to calculate the        |
    |                                     | relative magnification and       |
    |                                     | relative time-delay of unlensed  |
    |                                     | events.                          |
    +-------------------------------------+----------------------------------+
    |:meth:`~ mu_vs_dt_plot`              | Function to plot the             |
    |                                     | relative magnification vs        |
    |                                     | relative time-delay.             |
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

    cosmo = None
    """``astropy.cosmology`` \n
    Cosmology to use for the calculation.
    """

    size = None
    """``int`` \n
    Number of samples for sampling.
    """

    batch_size = None
    """``int`` \n
    Batch size for sampling.
    """

    json_file_names = None
    """``dict`` \n
    Names of the json files to strore the necessary parameters.
    """

    directory = None
    """``str`` \n
    Directory to store the interpolators.
    """

    gw_param_sampler_dict = None
    """``dict`` \n
    Dictionary of parameters to initialize the ``CBCSourceParameterDistribution`` class.
    """

    lens_param_sampler_dict = None
    """``dict`` \n
    Dictionary of parameters to initialize the ``LensGalaxyParameterDistribution`` class.
    """

    snr_calculator_dict = None
    """``dict`` \n
    Dictionary of parameters to initialize the ``GWSNR`` class.
    """

    list_of_detectors = None
    """``list`` \n
    List of detectors.
    """

    def __init__(
        self,
        npool=int(4),
        z_min=0.0,
        z_max=10.0,
        event_type="BBH",
        size=100000,
        batch_size=25000,
        cosmology=None,
        snr_finder="gwsnr",
        json_file_names=None,
        directory="./interpolator_pickle",
        verbose=True,
        **kwargs,
    ):
        
        self.npool = npool
        self.z_min = z_min
        self.z_max = z_max
        self.event_type = event_type
        self.cosmo = cosmology if cosmology else LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        self.size = size
        self.batch_size = batch_size
        self.json_file_names = dict(ler_param="./ler_params.json", unlensed_param="./unlensed_param.json", unlensed_param_detectable="./unlensed_param_detectable.json", lensed_param="./lensed_param.json", lensed_param_detectable="./lensed_param_detectable.json")
        if json_file_names:
            self.json_file_names.update(json_file_names)
        self.directory = directory

        def initialization():
            # initialization of parent class
            self.class_initialization(params=kwargs)
            if snr_finder == "gwsnr":
                # initialization self.snr and self.pdet from GWSNR class
                self.gwsnr_intialization(params=kwargs)
            else:
                self.snr = snr_finder
            # store all the ler input parameters
            self.store_ler_params(output_jsonfile=self.json_file_names["ler_param"])

        # if verbose, prevent anything from printing
        if verbose:
            initialization()
            self.print_all_params()
        else:
            with contextlib.redirect_stdout(None):
                initialization()
            
    def print_all_params(self):
        """
        Function to print all the parameters.
        """

        # print all relevant functions and sampler priors
        print("\n LeR set up params:")
        print("npool = ", self.npool)
        print("z_min = ", self.z_min)
        print("z_max = ", self.z_max)
        print("event_type = ", self.event_type)
        print("size = ", self.size)
        print("batch_size = ", self.batch_size)
        print("cosmology = ", self.cosmo)
        print("snr_finder = ", self.snr)
        print("json_file_names = ", self.json_file_names)
        print("directory = ", self.directory)
        
        print("\n LeR also takes CBCSourceParameterDistribution params as kwargs, as follows:")
        print("source_priors=", self.gw_param_sampler_dict["source_priors"])
        print("source_priors_params=", self.gw_param_sampler_dict["source_priors_params"])
        print("spin_zero=", self.gw_param_sampler_dict["spin_zero"])
        print("spin_precession=", self.gw_param_sampler_dict["spin_precession"])
        print("create_new_interpolator=", self.gw_param_sampler_dict["create_new_interpolator"])

        print("\n LeR also takes LensGalaxyParameterDistribution params as kwargs, as follows:")
        print("lens_type = ", self.gw_param_sampler_dict["lens_type"])
        print("lens_functions = ", self.gw_param_sampler_dict["lens_functions"])
        print("lens_priors = ", self.gw_param_sampler_dict["lens_priors"])
        print("lens_priors_params = ", self.gw_param_sampler_dict["lens_priors_params"])

        print("\n Image properties:")
        print("n_min_images = ", self.n_min_images)
        print("n_max_images = ", self.n_max_images)
        print("geocent_time_min = ", self.geocent_time_min)
        print("geocent_time_max = ", self.geocent_time_max)
        print("lens_model_list = ", self.lens_model_list)

        print("\n LeR also takes gwsnr.GWSNR params as kwargs, as follows:")
        print("mtot_min = ", self.snr_calculator_dict["mtot_min"])
        print("mtot_max = ", self.snr_calculator_dict["mtot_max"])
        print("ratio_min = ", self.snr_calculator_dict["ratio_min"])
        print("ratio_max = ", self.snr_calculator_dict["ratio_max"])
        print("mtot_resolution = ", self.snr_calculator_dict["mtot_resolution"])
        print("ratio_resolution = ", self.snr_calculator_dict["ratio_resolution"])
        print("sampling_frequency = ", self.snr_calculator_dict["sampling_frequency"])
        print("waveform_approximant = ", self.snr_calculator_dict["waveform_approximant"])
        print("minimum_frequency = ", self.snr_calculator_dict["minimum_frequency"])
        print("snr_type = ", self.snr_calculator_dict["snr_type"])
        print("psds = ", self.snr_calculator_dict["psds"])
        print("ifos = ", self.snr_calculator_dict["ifos"])
        print("interpolator_dir = ", self.snr_calculator_dict["interpolator_dir"])
        print("create_new_interpolator = ", self.snr_calculator_dict["create_new_interpolator"])
        print("gwsnr_verbose = ", self.snr_calculator_dict["gwsnr_verbose"])
        print("multiprocessing_verbose = ", self.snr_calculator_dict["multiprocessing_verbose"])
        print("mtot_cut = ", self.snr_calculator_dict["mtot_cut"])

        print("\n For reference, the chosen source parameters are listed below:")
        print("merger_rate_density = ", self.gw_param_samplers["merger_rate_density"])
        print("merger_rate_density_params = ", self.gw_param_samplers_params["merger_rate_density"])
        print("source_frame_masses = ", self.gw_param_samplers["source_frame_masses"])
        print("source_frame_masses_params = ", self.gw_param_samplers_params["source_frame_masses"])
        print("geocent_time = ", self.gw_param_samplers["geocent_time"])
        print("geocent_time_params = ", self.gw_param_samplers_params["geocent_time"])
        print("ra = ", self.gw_param_samplers["ra"])
        print("ra_params = ", self.gw_param_samplers_params["ra"])
        print("dec = ", self.gw_param_samplers["dec"])
        print("dec_params = ", self.gw_param_samplers_params["dec"])
        print("phase = ", self.gw_param_samplers["phase"])
        print("phase_params = ", self.gw_param_samplers_params["phase"])
        print("psi = ", self.gw_param_samplers["psi"])
        print("psi_params = ", self.gw_param_samplers_params["psi"])
        print("theta_jn = ", self.gw_param_samplers["theta_jn"])
        print("theta_jn_params = ", self.gw_param_samplers_params["theta_jn"])
        if self.spin_zero==False:
            print("a_1 = ", self.gw_param_samplers["a_1"])
            print("a_1_params = ", self.gw_param_samplers_params["a_1"])
            print("a_2 = ", self.gw_param_samplers["a_2"])
            print("a_2_params = ", self.gw_param_samplers_params["a_2"])
            if self.spin_precession==True:
                print("tilt_1 = ", self.gw_param_samplers["tilt_1"])
                print("tilt_1_params = ", self.gw_param_samplers_params["tilt_1"])
                print("tilt_2 = ", self.gw_param_samplers["tilt_2"])
                print("tilt_2_params = ", self.gw_param_samplers_params["tilt_2"])
                print("phi_12 = ", self.gw_param_samplers["phi_12"])
                print("phi_12_params = ", self.gw_param_samplers_params["phi_12"])
                print("phi_jl = ", self.gw_param_samplers["phi_jl"])
                print("phi_jl_params = ", self.gw_param_samplers_params["phi_jl"])

        print("\n For reference, the chosen lens related parameters and functions are listed below:")
        print("lens_redshift = ", self.lens_param_samplers["lens_redshift"])
        print("lens_redshift_params = ", self.lens_param_samplers_params["lens_redshift"])
        print("velocity_dispersion = ", self.lens_param_samplers["velocity_dispersion"])
        print("velocity_dispersion_params = ", self.lens_param_samplers_params["velocity_dispersion"])
        print("axis_ratio = ", self.lens_param_samplers["axis_ratio"])
        print("axis_ratio_params = ", self.lens_param_samplers_params["axis_ratio"])
        print("axis_rotation_angle = ", self.lens_param_samplers["axis_rotation_angle"])
        print("axis_rotation_angle_params = ", self.lens_param_samplers_params["axis_rotation_angle"])
        print("shear = ", self.lens_param_samplers["shear"])
        print("shear_params = ", self.lens_param_samplers_params["shear"])
        print("mass_density_spectral_index = ", self.lens_param_samplers["mass_density_spectral_index"])
        print("mass_density_spectral_index_params = ", self.lens_param_samplers_params["mass_density_spectral_index"])
        # lens functions
        print("Lens functions:")
        print("strong_lensing_condition = ", self.lens_functions["strong_lensing_condition"])
        print("optical_depth = ", self.lens_functions["optical_depth"])

    @property
    def snr(self):
        """
        Function to get the snr with the given parameters.

        Parameters
        ----------
        gw_param_dict : `dict`
            dictionary of GW source parameters.
            mass_1 : `numpy.ndarray` or `float`
                mass_1 of the compact binary (detector frame) (Msun).
            mass_2 : `numpy.ndarray` or `float`
                mass_2 of the compact binary (detector frame) (Msun).
            luminosity_distance : `numpy.ndarray` or `float`
                luminosity distance of the source (Mpc).
            theta_jn : `numpy.ndarray` or `float`
                inclination angle of the source (rad).
            psi : `numpy.ndarray` or `float`
                polarization angle of the source (rad).
            phase : `numpy.ndarray` or `float`
                phase of GW at reference frequency  (rad).
            geocent_time : `numpy.ndarray` or `float`
                GPS time of coalescence (s).
            ra : `numpy.ndarray` or `float`
                right ascension of the source (rad).
            dec : `numpy.ndarray` or `float`
                declination of the source (rad).
            a_1 : `numpy.ndarray` or `float`
                dimensionless spin magnitude of the more massive object.
            a_2 : `numpy.ndarray` or `float`
                dimensionless spin magnitude of the less massive object.
            tilt_1 : `numpy.ndarray` or `float`
                tilt angle of the more massive object spin.
            tilt_2 : `numpy.ndarray` or `float`
                tilt angle of the less massive object spin.
            phi_12 : `numpy.ndarray` or `float`
                azimuthal angle between the two spin vectors.
            phi_jl : `numpy.ndarray` or `float`
                azimuthal angle between total angular momentum and the orbital angular momentum.

        Returns
        ----------
        optimal_snr_list : `list`
            e.g. [optimal_snr_net, 'L1', 'H1', 'V1']
            optimal_snr_net : `numpy.ndarray` or `float`
                optimal snr of the network.
            'H1' : `numpy.ndarray` or `float`
                optimal snr of H1.
            'L1' : `numpy.ndarray` or `float`
                optimal snr of L1.
            'V1' : `numpy.ndarray` or `float`
                optimal snr of V1.
        """

        return self._snr
    
    @snr.setter
    def snr(self, snr_finder):
        self._snr = snr_finder

    @property
    def unlensed_param(self):
        """
        Function to get data from the json file self.json_file_names["unlensed_param"].

        Returns
        ----------
        unlensed_param : `dict`
            dictionary of unlensed GW source parameters.
        """

        return get_param_from_json(self.json_file_names["unlensed_param"])
    
    @property
    def unlensed_param_detectable(self):
        """
        Function to get data from the json file self.json_file_names["unlensed_param_detectable"].

        Returns
        ----------
        unlensed_param_detectable : `dict`
            dictionary of unlensed GW source parameters.
        """

        return get_param_from_json(self.json_file_names["unlensed_param_detectable"])
    
    @property
    def lensed_param(self):
        """
        Function to get data from the json file self.json_file_names["lensed_param"].

        Returns
        ----------
        lensed_param : `dict`
            dictionary of lensed GW source parameters.
        """

        return get_param_from_json(self.json_file_names["lensed_param"])
    
    @property
    def lensed_param_detectable(self):
        """
        Function to get data from the json file self.json_file_names["lensed_param_detectable"].

        Returns
        ----------
        lensed_param_detectable : `dict`
            dictionary of lensed GW source parameters.
        """

        return get_param_from_json(self.json_file_names["lensed_param_detectable"])

    def class_initialization(self, params=None):
        """
        Function to initialize the parent classes. List of relevant initialized instances, \n
        1. self.sample_source_redshift
        2. self.sample_unlensed_parameters
        3. self.normalization_pdf_z
        4. self.sample_lens_parameters
        5. self.normalization_pdf_z_lensed
        6. self.image_properties
        7. self.get_lensed_snrs

        Parameters
        ----------
        params : `dict`
            dictionary of parameters to initialize the parent classes
        """

        # initialization of LensGalaxyParameterDistribution class
        # it also initializes the CBCSourceParameterDistribution and ImageProperties classes
        input_params = dict(
            z_min=self.z_min,
            z_max=self.z_max,
            cosmology=self.cosmo,
            event_type=self.event_type,
            lens_type="epl_galaxy",
            lens_functions= None,
            lens_priors=None,
            lens_priors_params=None,
            source_priors=None,
            source_priors_params=None,
            spin_zero=True,
            spin_precession=False,
            directory=self.directory,
            create_new_interpolator=False,
        )
        if params:
            for key, value in params.items():
                if key in input_params:
                    input_params[key] = value
        self.gw_param_sampler_dict = input_params
        # initialization of clasess
        LensGalaxyParameterDistribution.__init__(
            self,
            z_min=input_params["z_min"],
            z_max=input_params["z_max"],
            cosmology=input_params["cosmology"],
            event_type=input_params["event_type"],
            lens_type=input_params["lens_type"],
            lens_functions=input_params["lens_functions"],
            lens_priors=input_params["lens_priors"],
            lens_priors_params=input_params["lens_priors_params"],
            source_priors=input_params["source_priors"],
            source_priors_params=input_params["source_priors_params"],
            spin_zero=input_params["spin_zero"],
            spin_precession=input_params["spin_precession"],
            directory=input_params["directory"],
            create_new_interpolator=input_params["create_new_interpolator"],
        )

        self.gw_param_sampler_dict["source_priors"]=self.gw_param_samplers.copy()
        self.gw_param_sampler_dict["source_priors_params"]=self.gw_param_samplers_params.copy()
        self.gw_param_sampler_dict["lens_priors"]=self.lens_param_samplers.copy()
        self.gw_param_sampler_dict["lens_priors_params"]=self.lens_param_samplers_params.copy()
        self.gw_param_sampler_dict["lens_functions"]=self.lens_functions.copy()

    def gwsnr_intialization(self, params=None):
        """
        Function to initialize the gwsnr class

        Parameters
        ----------
        params : `dict`
            dictionary of parameters to initialize the gwsnr class
        """
        from gwsnr import GWSNR

        # initialization of GWSNR class
        input_params = dict(
            npool=self.npool,
            mtot_min=2.0,
            mtot_max=200,
            ratio_min=0.1,
            ratio_max=1.0,
            mtot_resolution=500,
            ratio_resolution=50,
            sampling_frequency=2048.0,
            waveform_approximant="IMRPhenomD",
            minimum_frequency=20.0,
            snr_type="interpolation",
            psds=None,
            ifos=None,
            interpolator_dir=self.directory,
            create_new_interpolator=False,
            gwsnr_verbose=True,
            multiprocessing_verbose=True,
            mtot_cut=True,
        )
        # if self.event_type == "BNS":
        #     input_params["mtot_max"]= 18.
        if params:
            for key, value in params.items():
                if key in input_params:
                    input_params[key] = value
        self.snr_calculator_dict = input_params
        gwsnr = GWSNR(
                    npool=input_params["npool"],
                    mtot_min=input_params["mtot_min"],
                    mtot_max=input_params["mtot_max"],
                    ratio_min=input_params["ratio_min"],
                    ratio_max=input_params["ratio_max"],
                    mtot_resolution=input_params["mtot_resolution"],
                    ratio_resolution=input_params["ratio_resolution"],
                    sampling_frequency=input_params["sampling_frequency"],
                    waveform_approximant=input_params["waveform_approximant"],
                    minimum_frequency=input_params["minimum_frequency"],
                    snr_type=input_params["snr_type"],
                    psds=input_params["psds"],
                    ifos=input_params["ifos"],
                    interpolator_dir=input_params["interpolator_dir"],
                    # create_new_interpolator=input_params["create_new_interpolator"],
                    gwsnr_verbose=input_params["gwsnr_verbose"],
                    multiprocessing_verbose=input_params["multiprocessing_verbose"],
                    mtot_cut=input_params["mtot_cut"],
                )

        self.snr = gwsnr.snr
        self.list_of_detectors = gwsnr.detector_list
        #self.pdet = gwsnr.pdet

    def store_ler_params(self, output_jsonfile="./ler_params.json"):
        """
        Function to store the all the necessary parameters. This is useful for reproducing the results. All the parameters stored are in string format to make it json compatible.

        Parameters
        ----------
        output_jsonfile : `str`
            name of the json file to store the parameters
        """

        # store param_sampler_dict, lensed_param_sampler_dict and snr_calculator_dict
        parameters_dict = dict(
            npool=str(self.npool),
            z_min=str(self.z_min),
            z_max=str(self.z_max),
            size=str(self.size),
            batch_size=str(self.batch_size),
            cosmology=str(self.cosmo),
            snr_finder=str(self.snr),
            json_file_names=str(self.json_file_names),
            directory=str(self.directory),
        )

        # cbc params
        param_sampler_dict = self.gw_param_sampler_dict.copy()
        # convert all dict values to str
        for key, value in param_sampler_dict.items():
            param_sampler_dict[key] = str(value)
        parameters_dict.update({"gw_param_sampler_dict": param_sampler_dict})

        # snr calculator params
        try:
            snr_calculator_dict = self.snr_calculator_dict.copy()
            for key, value in snr_calculator_dict.items():
                snr_calculator_dict[key] = str(value)
            parameters_dict.update({"snr_calculator_dict": snr_calculator_dict})

            file_name = output_jsonfile
            append_json(file_name, parameters_dict, replace=True)
        except:
            # if snr_calculator is custom function
            pass

    def unlensed_cbc_statistics(
        self, size=None, resume=False, output_jsonfile=None,
    ):
        """
        Function to generate unlensed GW source parameters. This function also stores the parameters in json file.

        Parameters
        ----------
        size : `int`
            number of samples.
            default size = 100000.
        resume : `bool`
            resume = False (default) or True.
            if True, the function will resume from the last batch.
        output_jsonfile : `str`
            json file name for storing the parameters.
            default output_jsonfile = './unlensed_params.json'.

        Returns
        ----------
        unlensed_param : `dict`
            dictionary of unlensed GW source parameters.
            unlensed_param.keys() = ['zs', 'geocent_time', 'ra', 'dec', 'phase', 'psi', 'theta_jn', 'luminosity_distance', 'mass_1_source', 'mass_2_source', 'mass_1', 'mass_2', 'optimal_snr_net', 'L1', 'H1', 'V1']

        Examples
        ----------
        >>> from ler.rates import LeR
        >>> ler = LeR()
        >>> unlensed_param = ler.unlensed_cbc_statistics()
        """

        # gw parameter sampling
        if size is None:
            size = self.size

        # get json file name
        if output_jsonfile is None:
            output_jsonfile = self.json_file_names["unlensed_param"]
        else:
            self.json_file_names["unlensed_param"] = output_jsonfile
        print(f"unlensed params will be store in {output_jsonfile}")

        # sampling in batches
        batch_handler(
            size=size,
            batch_size=self.batch_size,
            sampling_routine=self.unlensed_sampling_routine,
            output_jsonfile=output_jsonfile,
            resume=resume,
        )

        unlensed_param = get_param_from_json(output_jsonfile)
        return unlensed_param
    
    def unlensed_sampling_routine(self, size, output_jsonfile, resume=False,
    ):
        """
        Function to generate unlensed GW source parameters. This function also stores the parameters in json file.

        Parameters
        ----------
        size : `int`
            number of samples.
            default size = 100000.
        resume : `bool`
            resume = False (default) or True.
            if True, the function will resume from the last batch.
        output_jsonfile : `str`
            json file name for storing the parameters.
            default output_jsonfile = './unlensed_params.json'.

        Returns
        ----------
        unlensed_param : `dict`
            dictionary of unlensed GW source parameters.
            unlensed_param.keys() = ['zs', 'geocent_time', 'ra', 'dec', 'phase', 'psi', 'theta_jn', 'luminosity_distance', 'mass_1_source', 'mass_2_source', 'mass_1', 'mass_2', 'optimal_snr_net', 'L1', 'H1', 'V1']
        """

        # get gw params
        print("sampling gw source params...")
        unlensed_param = self.sample_gw_parameters(size=size)
        # Get all of the signal to noise ratios
        print("calculating snrs...")
        snrs = self.snr(gw_param_dict=unlensed_param)
        unlensed_param.update(snrs)

        # store all params in json file
        append_json(file_name=output_jsonfile, dictionary=unlensed_param, replace=not (resume))

        return unlensed_param

    def unlensed_rate(
        self,
        unlensed_param=None,
        snr_threshold=8.0,
        output_jsonfile=None,
        detectability_condition="step_function",
    ):
        """
        Function to calculate the unlensed rate. This function also stores the parameters of the detectable events in json file.

        Parameters
        ----------
        unlensed_param : `dict` or `str`
            dictionary of GW source parameters or json file name.
            default unlensed_param = './unlensed_params.json'.
        snr_threshold : `float`
            threshold for detection signal to noise ratio.
            e.g. snr_threshold = 8.
        output_jsonfile : `str`
            json file name for storing the parameters of the detectable events.
            default output_jsonfile = './unlensed_params_detectable.json'.
        detectability_condition : `str`
            detectability condition. 
            default detectability_condition = 'step_function'.
            other options are 'pdet'.

        Returns
        ----------
        total_rate : `float`
            total unlensed rate (Mpc^-3 yr^-1).
        unlensed_param : `dict`
            dictionary of unlensed GW source parameters of the detectable events.
            unlensed_param.keys() = ['zs', 'geocent_time', 'ra', 'dec', 'phase', 'psi', 'theta_jn', 'luminosity_distance', 'mass_1_source', 'mass_2_source', 'mass_1', 'mass_2', 'optimal_snr_net', 'L1', 'H1', 'V1']

        Examples
        ----------
        >>> from ler.rates import LeR
        >>> ler = LeR()
        >>> ler.unlensed_cbc_statistics();
        >>> total_rate, unlensed_param_detectable = ler.unlensed_rate()
        """
        
        # call self.json_file_names["ler_param"] and for adding the final results
        data = load_json(self.json_file_names["ler_param"])
        
        # get gw params from json file if not provided
        if unlensed_param is None:
            unlensed_param = self.json_file_names["unlensed_param"]
        if type(unlensed_param) == str:
            self.json_file_names["unlensed_param"] = unlensed_param
            print(f"getting unlensed_params from json file {unlensed_param}...")
            unlensed_param = get_param_from_json(unlensed_param)
        else:
            print("using provided unlensed_param dict...")
            unlensed_param = unlensed_param.copy()
            
        if detectability_condition == "step_function":
            param = unlensed_param["optimal_snr_net"]
            threshold = snr_threshold

        elif detectability_condition == "pdet":
            # check if pdet is provided in unlensed_param dict
            if "pdet_net" in unlensed_param.keys():
                param = unlensed_param["pdet_net"]
            else:
                if "optimal_snr_net" in unlensed_param.keys():
                    param = 1 - norm.cdf(snr_threshold - unlensed_param["optimal_snr_net"])
                    unlensed_param["pdet_net"] = param
                else:
                    print("pdet or optimal_snr_net not provided in unlensed_param dict. Exiting...")
                    return None
            threshold = 0.5

        idx_detectable = param > threshold
        print(f"number of simulated unlensed detectable events: {np.sum(idx_detectable)}")
        print(f"number of all simulated unlensed events: {len(idx_detectable)}")
        # montecarlo integration
        # The total rate R = norm <Theta(rho-rhoc)>
        total_rate = self.normalization_pdf_z * np.mean(idx_detectable)
        print(f"total unlensed rate (yr^-1) (with step function): {total_rate}")

        # store all detectable params in json file
        for key, value in unlensed_param.items():
            unlensed_param[key] = value[idx_detectable]

        # store all detectable params in json file
        if output_jsonfile is None:
            output_jsonfile = self.json_file_names["unlensed_param_detectable"]
        else:
            self.json_file_names["unlensed_param_detectable"] = output_jsonfile
        print(f"storing detectable unlensed params in {output_jsonfile}")
        append_json(output_jsonfile, unlensed_param, replace=True)

        # write the results
        data['detectable_unlensed_rate_per_year'] = total_rate
        data["detectability_condition"] = detectability_condition
        append_json(self.json_file_names["ler_param"], data, replace=True)
        
        return total_rate, unlensed_param
    
    def lensed_cbc_statistics(
        self, size=None, resume=False, output_jsonfile=None,
    ):
        """
        Function to generate lensed GW source parameters. This function also stores the parameters in json file.

        Parameters
        ----------
        size : `int`
            number of samples.
            default size = 100000.
        resume : `bool`
            resume = False (default) or True.
            if True, the function will resume from the last batch.
        output_jsonfile : `str`
            json file name for storing the parameters.
            default output_jsonfile = './lensed_params.json'.

        Returns
        ----------
        lensed_param : `dict`
            dictionary of lensed GW source parameters.
            lensed_param.keys() = 

        Examples
        ----------
        >>> from ler.rates import LeR
        >>> ler = LeR()
        >>> lensed_param = ler.lensed_cbc_statistics()
        """

        # gw parameter sampling
        if size is None:
            size = self.size

        # get json file name
        if output_jsonfile is None:
            output_jsonfile = self.json_file_names["lensed_param"]
        else:
            self.json_file_names["lensed_param"] = output_jsonfile
        print(f"lensed params will be store in {output_jsonfile}")

        # sampling in batches
        batch_handler(
            size=size,
            batch_size=self.batch_size,
            sampling_routine=self.lensed_sampling_routine,
            output_jsonfile=output_jsonfile,
            resume=resume,
        )

        lensed_param = get_param_from_json(output_jsonfile)
        return lensed_param
    
    def lensed_sampling_routine(self, size, output_jsonfile, resume=False):
        """
        Function to generate lensed GW source parameters. This function also stores the parameters in json file.

        Parameters
        ----------
        size : `int`
            number of samples.
            default size = 100000.
        resume : `bool`
            resume = False (default) or True.
            if True, the function will resume from the last batch.
        output_jsonfile : `str`
            json file name for storing the parameters.
            default output_jsonfile = './lensed_params.json'.

        Returns
        ----------
        lensed_param : `dict`
            dictionary of lensed GW source parameters.
            lensed_param.keys() = 
        """

        # get lensed params
        print("sampling lensed params...")
        lensed_param = self.sample_lens_parameters(size=size)
        # now get (strongly lensed) image paramters along with lens parameters
        lensed_param = self.image_properties(lensed_param)
        # Get all of the signal to noise ratios
        print("calculating snrs...")
        snrs, lensed_param = self.get_lensed_snrs(
            snr_calculator=self.snr,
            list_of_detectors=self.list_of_detectors,
            lensed_param=lensed_param,
        )
        lensed_param.update(snrs)

        # store all params in json file
        append_json(file_name=output_jsonfile, dictionary=lensed_param, replace=not (resume))

        return lensed_param

    def lensed_rate(
        self,
        lensed_param=None,
        snr_threshold=[8.0,8.0],
        num_img=[1,1],
        output_jsonfile=None,
        nan_to_num=True,
        detectability_condition="step_function",
    ):
        """
        Function to calculate the lensed rate. This function also stores the parameters of the detectable events in json file.

        Parameters
        ----------
        lensed_param : `dict` or `str`
            dictionary of GW source parameters or json file name.
            default lensed_param = './lensed_params.json'.
        snr_threshold : `float`
            threshold for detection signal to noise ratio.
            default snr_threshold = [8.0,8.0].
        num_img : `int`
            number of images.
            default num_img = [1,1].
        output_jsonfile : `str`
            json file name for storing the parameters of the detectable events.
            default output_jsonfile = './lensed_params_detectable.json'.
        nan_to_num : `bool`
            if True, nan values will be converted to 0.
            default nan_to_num = True.
        detectability_condition : `str`
            detectability condition.
            default detectability_condition = 'step_function'.
            other options are 'pdet'.

        Returns
        ----------
        total_rate : `float`
            total lensed rate (Mpc^-3 yr^-1).
        lensed_param : `dict`
            dictionary of lensed GW source parameters of the detectable events.
            lensed_param.keys() =

        Examples
        ----------
        >>> from ler.rates import LeR
        >>> ler = LeR()
        >>> ler.lensed_cbc_statistics();
        >>> total_rate, lensed_param_detectable = ler.lensed_rate()
        """
        
        # call self.json_file_names["ler_param"] and for adding the final results
        data = load_json(self.json_file_names["ler_param"])

        # get gw params from json file if not provided
        if lensed_param is None:
            lensed_param = self.json_file_names["lensed_param"]
        if type(lensed_param) == str:
            self.json_file_names["lensed_param"] = lensed_param
            print(f"getting lensed_params from json file {lensed_param}...")
            lensed_param = get_param_from_json(lensed_param)
        else:
            print("using provided lensed_param dict...")
            lensed_param = lensed_param.copy()

        # check for images with snr above threshold
        # convert to array
        snr_threshold = np.array([snr_threshold]).reshape(-1)  
        num_img = np.array([num_img]).reshape(-1)
        # get descending sorted idx of snr_threshold
        idx = np.argsort(-snr_threshold)
        snr_threshold = snr_threshold[idx]
        num_img = num_img[idx]

        # get size of the lensed_param for a parameter
        size = len(lensed_param["zs"])

        if detectability_condition == "step_function":
            snr_hit = np.full(size, True)  # boolean array to store the result of the threshold condition
            try:
                snr_param = lensed_param["optimal_snr_net"]
                snr_param = -np.sort(-snr_param, axis=1)  # sort snr in descending order
            except:
                print("snr not provided in lensed_param dict. Exiting...")
                return None
            
            # for each row: choose a threshold and check if the number of images above threshold. Sum over the images. If sum is greater than num_img, then snr_hit = True 
            j = 0
            idx_max = 0
            for i in range(len(snr_threshold)):
                idx_max = idx_max + num_img[i]
                snr_hit = snr_hit & (np.sum((snr_param[:,j:idx_max] > snr_threshold[i]), axis=1) >= num_img[i])
                j = idx_max

            # montecarlo integration
            total_rate = self.normalization_pdf_z_lensed * np.mean(snr_hit)
            print("total lensed rate (yr^-1) (with step function): {}".format(total_rate))

        elif detectability_condition == "pdet":
            # check if pdet is provided in unlensed_param dict
            if "pdet_net" in lensed_param.keys():
                pdet = lensed_param["pdet_net"]
                pdet = -np.sort(-pdet, axis=1)  # sort pdet in descending order
                pdet = pdet[:,:np.sum(num_img)]  # use only num_img images
            else:
                if "optimal_snr_net" in lensed_param.keys():
                    # pdet dimension is (size, n_max_images)
                    snr_param = lensed_param["optimal_snr_net"]
                    snr_param = -np.sort(-snr_param, axis=1)  # sort snr in descending order

                    pdet = np.ones(np.shape(snr_param))
                    j = 0
                    idx_max = 0
                    for i in range(len(snr_threshold)):
                        idx_max = idx_max + num_img[i]
                        pdet[:,j:idx_max] = 1 - norm.cdf(snr_threshold[i] - snr_param[:,j:idx_max])
                        j = idx_max
                else:
                    print("pdet or optimal_snr_net not provided in lensed_param dict. Exiting...")
                    return None
                
            snr_hit = np.prod(pdet, axis=1)>0.5
            # montecarlo integration
            total_rate = self.normalization_pdf_z_lensed * np.mean(snr_hit)
            print(f"total lensed rate (yr^-1) (with pdet function): {total_rate}")

        print(f"number of simulated lensed detectable events: {np.sum(snr_hit)}")
        print(f"number of simulated all lensed events: {size}")
        # store all detectable params in json file
        if nan_to_num:
            for key, value in lensed_param.items():
                lensed_param[key] = np.nan_to_num(value[snr_hit])
        else:
            for key, value in lensed_param.items():
                lensed_param[key] = value[snr_hit]
            
        # store all detectable params in json file
        if output_jsonfile is None:
            output_jsonfile = self.json_file_names["lensed_param_detectable"]
        else:
            self.json_file_names["lensed_param_detectable"] = output_jsonfile
        print(f"storing detectable lensed params in {output_jsonfile}")
        append_json(output_jsonfile, lensed_param, replace=True)

        # write the results
        data["detectable_lensed_rate_per_year"] = total_rate
        data["detectability_condition"] = detectability_condition
        append_json(self.json_file_names["ler_param"], data, replace=True)

        return total_rate, lensed_param
     
    def rate_ratio(self):
        """
        Function to calculate and display unlensed and lensed merger rate ratio. 
        It will get the unlensed_rate and lensed_rate from self.json_file_ler_param

        Returns
        ----------
        rate_ratio : `float`
            rate ratio.

        Examples
        ----------
        >>> from ler.rates import LeR
        >>> ler = LeR()
        >>> ler.unlensed_cbc_statistics();
        >>> ler.lensed_cbc_statistics();
        >>> ler.unlensed_rate();
        >>> ler.lensed_rate();
        >>> ler.rate_ratio()
        """

        # call json_file_ler_param and add the results
        data = load_json(self.json_file_names["ler_param"])

        try:
            unlensed_rate = data["detectable_unlensed_rate_per_year"]
            lensed_rate = data["detectable_lensed_rate_per_year"]
        except:
            raise ValueError(f"unlensed_rate_step or lensed_rate_step not found in {self.json_file_names['ler_param']} json file. Exiting...")
        
        rate_ratio = unlensed_rate / lensed_rate
        # append the results
        data['rate_ratio'] = rate_ratio
        # write the results
        append_json(self.json_file_names["ler_param"], data, replace=True)
        
        print(f"unlensed_rate: {unlensed_rate}")
        print(f"lensed_rate: {lensed_rate}")
        print(f"ratio: {rate_ratio}")

        return unlensed_rate / lensed_rate
    
    def rate_comparision_with_rate_calculation(
        self,
        unlensed_param=None,
        snr_threshold_unlensed=8.0,
        output_jsonfile_unlensed=None,
        lensed_param=None,
        snr_threshold_lensed=[8.0,8.0],
        num_img=[1,1],
        output_jsonfile_lensed=None,
        nan_to_num=True,
        detectability_condition="step_function",
    ):
        """
        Function to calculate the unlensed and lensed rate and compare by computing the ratio. This function also stores the parameters of the detectable events in json file.

        Parameters
        ----------
        unlensed_param : `dict` or `str`
            dictionary of GW source parameters or json file name.
            default unlensed_param = './unlensed_params.json'.
        snr_threshold_unlensed : `float`
            threshold for detection signal to noise ratio.
            e.g. snr_threshold = 8.
        output_jsonfile_unlensed : `str`
            json file name for storing the parameters of the detectable events.
            default output_jsonfile = './unlensed_params_detectable.json'.
        lensed_param : `dict` or `str`
            dictionary of GW source parameters or json file name.
            default lensed_param = './lensed_params.json'.
        snr_threshold_lensed : `float`
            threshold for detection signal to noise ratio.
            e.g. snr_threshold = 8.
        output_jsonfile_lensed : `str`
            json file name for storing the parameters of the detectable events.
            default output_jsonfile = './lensed_params_detectable.json'.
        detectability_condition : `str`
            detectability condition. 
            default detectability_condition = 'step_function'.
            other options are 'pdet'.

        Returns
        ----------
        rate_ratio : `float`
            rate ratio.
        unlensed_param : `dict`
            dictionary of unlensed GW source parameters of the detectable events.
        lensed_param : `dict`
            dictionary of lensed GW source parameters of the detectable events.

        Examples
        ----------
        >>> from ler.rates import LeR
        >>> ler = LeR()
        >>> ler.unlensed_cbc_statistics();
        >>> ler.lensed_cbc_statistics();
        >>> rate_ratio, unlensed_param, lensed_param = ler.rate_comparision_with_rate_calculation()
        """

        # call json_file_ler_param and add the results
        data = load_json(self.json_file_names["ler_param"])

        # get unlensed rate
        unlensed_rate, unlensed_param = self.unlensed_rate(
            unlensed_param=unlensed_param,
            snr_threshold=snr_threshold_unlensed,
            output_jsonfile=output_jsonfile_unlensed,
            detectability_condition=detectability_condition,
        )
        # get lensed rate
        lensed_rate, lensed_param = self.lensed_rate(
            lensed_param=lensed_param,
            snr_threshold=snr_threshold_lensed,
            num_img=num_img,
            output_jsonfile=output_jsonfile_lensed,
            nan_to_num=nan_to_num,
            detectability_condition=detectability_condition,
        )
        # calculate rate ratio
        rate_ratio = unlensed_rate / lensed_rate
        # append the results
        data['detectable_unlensed_rate_per_year'] = unlensed_rate
        data['detectable_lensed_rate_per_year'] = lensed_rate
        data['rate_ratio'] = rate_ratio
        data['detectability_condition'] = detectability_condition
        # write the results
        append_json(self.json_file_names["ler_param"], data, replace=True)
        
        print(f"unlensed_rate: {unlensed_rate}")
        print(f"lensed_rate: {lensed_rate}")
        print(f"ratio: {rate_ratio}")

        return rate_ratio, unlensed_param, lensed_param

    def selecting_n_unlensed_detectable_events(
        self,
        size=100,
        batch_size=None,
        snr_threshold=8.0,
        resume=False,
        output_jsonfile="./n_unlensed_param_detectable.json",
        meta_data_file="meta_unlensed.json",
        event_batch_limit=500,
        trim_to_size=True,
    ):
        """
        Function to select n unlensed detectable events.

        Parameters
        ----------
        size : `int`
            number of samples to be selected.
            default size = 100.
        snr_threshold : `float`
            threshold for detection signal to noise ratio.
            e.g. snr_threshold = 8.
        resume : `bool`
            if True, it will resume the sampling from the last batch.
            default resume = False.
        output_jsonfile : `str`
            json file name for storing the parameters.
            default output_jsonfile = './n_unlensed_params_detectable.json'.

        Returns
        ----------
        param_final : `dict`
            dictionary of unlensed GW source parameters of the detectable events.
            param_final.keys() = ['zs', 'geocent_time', 'ra', 'dec', 'phase', 'psi', 'theta_jn', 'luminosity_distance', 'mass_1_source', 'mass_2_source', 'mass_1', 'mass_2', 'optimal_snr_net', 'L1', 'H1', 'V1']

        Examples
        ----------
        >>> from ler.rates import LeR
        >>> ler = LeR()
        >>> unlensed_param_final = ler.selecting_n_unlensed_detectable_events(size=500)
        """

        if batch_size is None:
            batch_size = self.batch_size
        else:
            self.batch_size = batch_size

        if not resume:
            n = 0  # iterator
            events_total = 0
            try:
                os.remove(output_jsonfile)
            except:
                pass
        else:
            # get sample size as size from json file
            try:
                param_final = get_param_from_json(output_jsonfile)
                n = len(param_final["zs"])
                events_total = load_json(meta_data_file)["events_total"]
                del param_final
            except:
                n = 0
                events_total = 0

        buffer_file = "./unlensed_params_buffer.json"
        print("collected number of events = ", n)
        
        while n < size:
            # disable print statements
            with contextlib.redirect_stdout(None):
                unlensed_param = self.unlensed_sampling_routine(
                    size=batch_size, output_jsonfile=buffer_file, resume=False
                )

            # get snr
            snr = unlensed_param["optimal_snr_net"]
            # snr_ = self.snr(gw_param_dict=unlensed_param)["optimal_snr_net"]
            # test = snr == snr_
            # index of detectable events
            idx = snr > snr_threshold

            # store all params in json file
            if np.sum(idx) < event_batch_limit:
                for key, value in unlensed_param.items():
                    unlensed_param[key] = value[idx]
                append_json(output_jsonfile, unlensed_param, replace=False)

                n += np.sum(idx)
                events_total += len(idx)                
                total_rate = self.normalization_pdf_z * n / events_total

                # save meta data
                meta_data = dict(events_total=events_total, total_rate=total_rate)
                try:
                    append_json(meta_data_file, meta_data, replace=True)
                except:
                    append_json(meta_data_file, meta_data, replace=False)
            
            # else:
            #     snr_param = self.snr(gw_param_dict=unlensed_param)['optimal_snr_net']
            #     print(snr_param)

            print("collected number of events = ", n)
            print("total number of events = ", events_total)
            print(f"total unlensed rate (yr^-1): {total_rate}")

        print(f"storing detectable unlensed params in {output_jsonfile}")

        if trim_to_size:
            # trim the final param dictionary
            print(f"\n trmming final result to size={size}")
            param_final = get_param_from_json(output_jsonfile)
            # trim the final param dictionary, randomly, without repeating
            # idx_ = np.random.choice(len(param_final["zs"]), size, replace=False)
            # for key, value in param_final.items():
            #     param_final[key] = param_final[key][idx_]
            for key, value in param_final.items():
                param_final[key] = param_final[key][:size]

            # save the final param dictionary
            
            append_json(output_jsonfile, param_final, replace=True)

        return param_final
    
    def selecting_n_lensed_detectable_events(
        self,
        size=100,
        batch_size=None,
        snr_threshold=8.0,
        num_img=2,
        resume=False,
        detectability_condition="step_function",
        output_jsonfile="./n_lensed_params_detectable.json",
        meta_data_file="./meta_lensed.json",
        event_batch_limit=500,
        trim_to_size=True,
        nan_to_num=False,
    ):
        """
        Function to select n lensed detectable events.

        Parameters
        ----------
        size : `int`
            number of samples to be selected.
            default size = 100.
        snr_threshold : `float`
            threshold for detection signal to noise ratio.
            e.g. snr_threshold = 8.
        num_img : `int`
            number of images.
            default num_img = 2.
        resume : `bool`
            if True, it will resume the sampling from the last batch.
            default resume = False.
        detectability_condition : `str`
            detectability condition.
            default detectability_condition = 'step_function'.
            other options are 'pdet'.
        output_jsonfile : `str`
            json file name for storing the parameters.

        Returns
        ----------
        param_final : `dict`
            dictionary of lensed GW source parameters of the detectable events.
            param_final.keys() =

        Examples
        ----------
        >>> from ler.rates import LeR
        >>> ler = LeR()
        >>> lensed_param_final = ler.selecting_n_lensed_detectable_events(size=500)
        """
        
        if batch_size is None:
            batch_size = self.batch_size
        else:
            self.batch_size = batch_size

        if not resume:
            n = 0  # iterator
            events_total = 0
            try:
                os.remove(output_jsonfile)
            except:
                pass
        else:
            # get sample size as size from json file
            try:
                param_final = get_param_from_json(output_jsonfile)
                n = len(param_final["zs"])
                events_total = load_json(meta_data_file)["events_total"]
                del param_final
            except:
                n = 0
                events_total = 0

        # check for images with snr above threshold
        # convert to array
        snr_threshold = np.array([snr_threshold]).reshape(-1)  
        num_img = np.array([num_img]).reshape(-1)
        # get descending sorted idx of snr_threshold
        idx = np.argsort(-snr_threshold)
        snr_threshold = snr_threshold[idx]
        num_img = num_img[idx]

        buffer_file = "./lensed_params_buffer.json"
        print("collected number of events = ", n)
        while n < size:
            # disable print statements
            with contextlib.redirect_stdout(None):
                lensed_param = self.lensed_sampling_routine(
                    size=self.batch_size, output_jsonfile=buffer_file, resume=False
                )  # Dimensions are (size, n_max_images)

            if detectability_condition == "step_function":
                snr_hit = np.full(len(lensed_param["zs"]), True)  # boolean array to store the result of the threshold condition
                try:
                    snr_param = lensed_param["optimal_snr_net"]
                    snr_param = -np.sort(-snr_param, axis=1)  # sort snr in descending order
                except:
                    raise ValueError("snr not provided in lensed_param dict. Exiting...")
                
                # for each row: choose a threshold and check if the number of images above threshold. Sum over the images. If sum is greater than num_img, then snr_hit = True 
                j = 0
                idx_max = 0
                for i in range(len(snr_threshold)):
                    idx_max = idx_max + num_img[i]
                    snr_hit = snr_hit & (np.sum((snr_param[:,j:idx_max] > snr_threshold[i]), axis=1) >= num_img[i])
                    j = idx_max

            elif detectability_condition == "pdet":
                # check if pdet is provided in unlensed_param dict
                if "pdet_net" in lensed_param.keys():
                    pdet = lensed_param["pdet_net"]
                    pdet = -np.sort(-pdet, axis=1)  # sort pdet in descending order
                    pdet = pdet[:,:np.sum(num_img)]  # use only num_img images
                else:
                    if "optimal_snr_net" in lensed_param.keys():
                        # pdet dimension is (size, n_max_images)
                        snr_param = lensed_param["optimal_snr_net"]
                        snr_param = -np.sort(-snr_param, axis=1)  # sort snr in descending order

                        pdet = np.ones(np.shape(snr_param))
                        j = 0
                        idx_max = 0
                        for i in range(len(snr_threshold)):
                            idx_max = idx_max + num_img[i]
                            pdet[:,j:idx_max] = 1 - norm.cdf(snr_threshold[i] - snr_param[:,j:idx_max])
                            j = idx_max
                    else:
                        print("pdet or optimal_snr_net not provided in lensed_param dict. Exiting...")
                        return None
                    
                snr_hit = np.prod(pdet, axis=1)>0.5
                    
            # store all params in json file
            if np.sum(snr_hit)<event_batch_limit:
                if nan_to_num:
                    for key, value in lensed_param.items():
                        lensed_param[key] = np.nan_to_num(value[snr_hit])
                else:
                    for key, value in lensed_param.items():
                        lensed_param[key] = value[snr_hit]
                try:
                    append_json(output_jsonfile, lensed_param, replace=False)
                except:
                    append_json(output_jsonfile, lensed_param, replace=True)

                n += np.sum(snr_hit)
                events_total += len(snr_hit)

                # save meta data
                total_rate = self.normalization_pdf_z_lensed * n / events_total
                meta_data = dict(events_total=events_total, total_rate=total_rate)
                append_json(meta_data_file, meta_data, replace=True)

            print("collected number of events = ", n)
            print("total number of events = ", events_total)
            print(f"total lensed rate (yr^-1): {total_rate}")

        print(f"storing detectable lensed params in {output_jsonfile}")

        if trim_to_size:
            # trim the final param dictionary
            print(f"\n trmming final result to size={size}")
            param_final = get_param_from_json(output_jsonfile)
            # trim the final param dictionary
            #idx = np.random.choice(len(param_final["zs"]), size, replace=False)
            for key, value in param_final.items():
                param_final[key] = param_final[key][:size]

            # save the final param dictionary
            append_json(output_jsonfile, param_final, replace=True)

        return param_final



