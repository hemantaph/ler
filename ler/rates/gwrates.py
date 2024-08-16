# -*- coding: utf-8 -*-
"""
This module contains the main class for calculating the rates of detectable gravitational waves events. The class inherits the :class:`~ler.gw_source_population.CBCSourceParameterDistribution` class for source parameters sampling and uses `gwsnr` package for SNR calculation. 
"""

import os
import warnings
warnings.filterwarnings("ignore")
import contextlib
import numpy as np
from scipy.stats import norm
from astropy.cosmology import LambdaCDM
from ..gw_source_population import CBCSourceParameterDistribution
from ..utils import load_json, append_json, get_param_from_json, batch_handler


class GWRATES(CBCSourceParameterDistribution):
    """Class to calculate both the rates of lensed and gw events. Please note that parameters of the simulated events are stored in json file but not as an attribute of the class. This saves RAM memory. 

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
        for popI_II, popIII, primordial, BNS z_max = 10., 40., 40., 5. respectively.
    event_type : `str`
        type of event to generate.
        default event_type = 'BBH'. Other options are 'BNS', 'NSBH'.
    size : `int`
        number of samples for sampling.
        default size = 100000. To get stable rates, size should be large (>=1e6).
    batch_size : `int`
        batch size for SNR calculation.
        default batch_size = 50000.
        reduce the batch size if you are getting memory error.
        recommended batch_size = 200000, if size = 1000000.
    cosmology : `astropy.cosmology`
        cosmology to use for the calculation.
        default cosmology = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7).
    snr_finder : `str` or `function`
        default snr_finder = 'gwsnr'.
        if None, the SNR will be calculated using the gwsnr package.
        if custom snr finder function is provided, the SNR will be calculated using a custom function. The custom function should follow the following signature:
        def snr_finder(gw_param_dict):
            ...
            return optimal_snr_dict
        where optimal_snr_dict.keys = ['optimal_snr_net']. Refer to `gwsnr` package's GWSNR.snr attribute for more details.
    pdet_finder : `function`
        default pdet_finder = None.
        The rate calculation uses either the pdet_finder or the snr_finder to calculate the detectable events. The custom pdet finder function should follow the following signature:
        def pdet_finder(gw_param_dict):
            ...
            return pdet_net_dict
        where pdet_net_dict.keys = ['pdet_net']. For example uses, refer to [GRB pdet example](https://ler.readthedocs.io/en/latest/examples/rates/grb%20detection%20rate.html).
    json_file_names: `dict`
        names of the json files to strore the necessary parameters.
        default json_file_names = {'gwrates_params':'gwrates_params.json', 'gw_param': 'gw_param.json', 'gw_param_detectable': 'gw_param_detectable.json'}.
    interpolator_directory : `str`
        directory to store the interpolators.
        default interpolator_directory = './interpolator_pickle'. This is used for storing the various interpolators related to `ler` and `gwsnr` package.
    ler_directory : `str`
        directory to store the parameters.
        default ler_directory = './ler_data'. This is used for storing the parameters of the simulated events.
    verbose : `bool`
        default verbose = True.
        if True, the function will print all chosen parameters.
        Choose False to prevent anything from printing.
    kwargs : `keyword arguments`
        Note : kwargs takes input for initializing the :class:`~ler.gw_source_population.CBCSourceParameterDistribution` and :class:`~ler.gw_source_population.CBCSourceRedshiftDistribution` classes. If snr_finder='gwsnr', then kwargs also takes input for initializing the :class:`~gwsnr.GWSNR` class. Please refer to the respective classes for more details.

    Examples
    ----------
    >>> from ler.rates import GWRATES
    >>> ler = GWRATES()
    >>> ler.gw_cbc_statistics();
    >>> ler.gw_rate();
        
    Instance Attributes
    ----------
    LeR class has the following attributes, \n
    +-------------------------------------+----------------------------------+
    | Atrributes                          | Type                             |
    +=====================================+==================================+
    |:attr:`~npool`                       | `int`                            |
    +-------------------------------------+----------------------------------+
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
    |:attr:`~interpolator_directory`      | `str`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~ler_directory`               | `str`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~gwsnr`                       | `bool`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~gw_param_sampler_dict`       | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~snr_calculator_dict`         | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~gw_param`                    | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~gw_param_detectable`         | `dict`                           |
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
    |:meth:`~store_gwrates_params`        | Function to store the all the    |
    |                                     | necessary parameters.            |
    +-------------------------------------+----------------------------------+
    |:meth:`~gw_cbc_statistics`           | Function to generate gw          |
    |                                     | GW source parameters.            |
    +-------------------------------------+----------------------------------+
    |:meth:`~gw_sampling_routine`         | Function to generate gw          |
    |                                     | GW source parameters.            |
    +-------------------------------------+----------------------------------+
    |:meth:`~gw_rate`                     | Function to calculate the        |
    |                                     | gw rate.                         |
    +-------------------------------------+----------------------------------+
    |:meth:`~selecting_n_gw_detectable_events`                               |
    +-------------------------------------+----------------------------------+
    |                                     | Function to select n gw    |
    |                                     | detectable events.               |
    +-------------------------------------+----------------------------------+
    |:meth:`~gw_param_plot`               | Function to plot the             |
    |                                     | distribution of the GW source    |
    |                                     | parameters.                      |
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
    Names of the json files to store the necessary parameters.
    """

    interpolator_directory = None
    """``str`` \n
    Directory to store the interpolators.
    """

    ler_directory = None
    """``str`` \n
    Directory to store the parameters.
    """

    gwsnr = None
    """``bool`` \n
    If True, the SNR will be calculated using the gwsnr package.
    """

    gw_param_sampler_dict = None
    """``dict`` \n
    Dictionary of parameters to initialize the ``CBCSourceParameterDistribution`` class.
    """

    snr_calculator_dict = None
    """``dict`` \n
    Dictionary of parameters to initialize the ``GWSNR`` class.
    """

    gw_param = None
    """``dict`` \n
    Dictionary of GW source parameters. The included parameters and their units are as follows (for default settings):\n
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
    | L1                 |              | optimal snr of L1                    |
    +--------------------+--------------+--------------------------------------+
    | H1                 |              | optimal snr of H1                    |
    +--------------------+--------------+--------------------------------------+
    | V1                 |              | optimal snr of V1                    |
    +--------------------+--------------+--------------------------------------+
    | optimal_snr_net    |              | optimal snr of the network           |
    +--------------------+--------------+--------------------------------------+
    """

    gw_param_detectable = None
    """``dict`` \n
    Dictionary of detectable GW source parameters. It includes the same parameters as the :attr:`~gw_param` attribute.
    """

    def __init__(
        self,
        npool=int(4),
        z_min=0.0,
        z_max=10.0,
        event_type="BBH",
        size=100000,
        batch_size=50000,
        cosmology=None,
        snr_finder=None,
        pdet_finder=None,
        list_of_detectors=None,
        json_file_names=None,
        interpolator_directory="./interpolator_pickle",
        ler_directory="./ler_data",
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
        self.json_file_names = dict(gwrates_params="gwrates_params.json", gw_param="gw_param.json", gw_param_detectable="gw_param_detectable.json",)
        if json_file_names:
            self.json_file_names.update(json_file_names)
        self.interpolator_directory = interpolator_directory
        self.ler_directory = ler_directory
        # create directory if not exists
        if not os.path.exists(ler_directory):
            os.makedirs(ler_directory)

        def initialization():
            # initialization of parent class
            self.class_initialization(params=kwargs)
            # initialization self.snr and self.pdet from GWSNR class
            if not snr_finder and not pdet_finder:
                self.gwsnr_intialization(params=kwargs)
                self.gwsnr = True
                self.pdet = pdet_finder
            else:
                self.snr = snr_finder
                self.pdet = pdet_finder
                self.gwsnr = False
                self.list_of_detectors = list_of_detectors

            # store all the gwrates input parameters
            self.store_gwrates_params(output_jsonfile=self.json_file_names["gwrates_params"])

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
        print("\n GWRATES set up params:")
        print(f'npool = {self.npool},')
        print(f'z_min = {self.z_min},')
        print(f'z_max = {self.z_max},')
        print(f"event_type = '{self.event_type}',")
        print(f'size = {self.size},')
        print(f'batch_size = {self.batch_size},')
        print(f'cosmology = {self.cosmo},')
        if self.snr:
            print(f'snr_finder = {self.snr},')
        if self.pdet:
            print(f'pdet_finder = {self.pdet},')
        print(f'json_file_names = {self.json_file_names},')
        print(f'interpolator_directory = {self.interpolator_directory},')
        print(f'ler_directory = {self.ler_directory},')

        print("\n GWRATES also takes CBCSourceParameterDistribution params as kwargs, as follows:")
        print(f"source_priors = {self.gw_param_sampler_dict['source_priors']},")
        print(f"source_priors_params = {self.gw_param_sampler_dict['source_priors_params']},")
        print(f"spin_zero = {self.gw_param_sampler_dict['spin_zero']},")
        print(f"spin_precession = {self.gw_param_sampler_dict['spin_precession']},")
        print(f"create_new_interpolator = {self.gw_param_sampler_dict['create_new_interpolator']},")

        if self.gwsnr:
            print("\n LeR also takes gwsnr.GWSNR params as kwargs, as follows:")
            print(f"mtot_min = {self.snr_calculator_dict['mtot_min']},")
            print(f"mtot_max = {self.snr_calculator_dict['mtot_max']},")
            print(f"ratio_min = {self.snr_calculator_dict['ratio_min']},")
            print(f"ratio_max = {self.snr_calculator_dict['ratio_max']},")
            print(f"mtot_resolution = {self.snr_calculator_dict['mtot_resolution']},")
            print(f"ratio_resolution = {self.snr_calculator_dict['ratio_resolution']},")
            print(f"sampling_frequency = {self.snr_calculator_dict['sampling_frequency']},")
            print(f"waveform_approximant = '{self.snr_calculator_dict['waveform_approximant']}',")
            print(f"minimum_frequency = {self.snr_calculator_dict['minimum_frequency']},")
            print(f"snr_type = '{self.snr_calculator_dict['snr_type']}',")
            print(f"psds = {self.snr_calculator_dict['psds']},")
            print(f"ifos = {self.snr_calculator_dict['ifos']},")
            print(f"interpolator_dir = '{self.snr_calculator_dict['interpolator_dir']}',")
            print(f"create_new_interpolator = {self.snr_calculator_dict['create_new_interpolator']},")
            print(f"gwsnr_verbose = {self.snr_calculator_dict['gwsnr_verbose']},")
            print(f"multiprocessing_verbose = {self.snr_calculator_dict['multiprocessing_verbose']},")
            print(f"mtot_cut = {self.snr_calculator_dict['mtot_cut']},")
        # del self.gwsnr

        print("\n For reference, the chosen source parameters are listed below:")
        print(f"merger_rate_density = '{self.gw_param_samplers['merger_rate_density']}'")
        print("merger_rate_density_params = ", self.gw_param_samplers_params["merger_rate_density"])
        print(f"source_frame_masses = '{self.gw_param_samplers['source_frame_masses']}'")
        print("source_frame_masses_params = ", self.gw_param_samplers_params["source_frame_masses"])
        print(f"geocent_time = '{self.gw_param_samplers['geocent_time']}'")
        print("geocent_time_params = ", self.gw_param_samplers_params["geocent_time"])
        print(f"ra = '{self.gw_param_samplers['ra']}'")
        print("ra_params = ", self.gw_param_samplers_params["ra"])
        print(f"dec = '{self.gw_param_samplers['dec']}'")
        print("dec_params = ", self.gw_param_samplers_params["dec"])
        print(f"phase = '{self.gw_param_samplers['phase']}'")
        print("phase_params = ", self.gw_param_samplers_params["phase"])
        print(f"psi = '{self.gw_param_samplers['psi']}'")
        print("psi_params = ", self.gw_param_samplers_params["psi"])
        print(f"theta_jn = '{self.gw_param_samplers['theta_jn']}'")
        print("theta_jn_params = ", self.gw_param_samplers_params["theta_jn"])
        if self.spin_zero is False:
            print(f"a_1 = '{self.gw_param_samplers['a_1']}'")
            print("a_1_params = ", self.gw_param_samplers_params["a_1"])
            print(f"a_2 = '{self.gw_param_samplers['a_2']}'")
            print("a_2_params = ", self.gw_param_samplers_params["a_2"])
            if self.spin_precession is True:
                print(f"tilt_1 = '{self.gw_param_samplers['tilt_1']}'")
                print("tilt_1_params = ", self.gw_param_samplers_params["tilt_1"])
                print(f"tilt_2 = '{self.gw_param_samplers['tilt_2']}'")
                print("tilt_2_params = ", self.gw_param_samplers_params["tilt_2"])
                print(f"phi_12 = '{self.gw_param_samplers['phi_12']}'")
                print("phi_12_params = ", self.gw_param_samplers_params["phi_12"])
                print(f"phi_jl = '{self.gw_param_samplers['phi_jl']}'")
                print("phi_jl_params = ", self.gw_param_samplers_params["phi_jl"])

        
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
    def gw_param(self):
        """
        Function to get data from the json file self.json_file_names["gw_param"].

        Returns
        ----------
        gw_param : `dict`
            dictionary of gw GW source parameters.
        """

        return get_param_from_json(self.ler_directory+"/"+self.json_file_names["gw_param"])
    
    @property
    def gw_param_detectable(self):
        """
        Function to get data from the json file self.json_file_names["gw_param_detectable"].

        Returns
        ----------
        gw_param_detectable : `dict`
            dictionary of gw GW source parameters.
        """

        return get_param_from_json(self.ler_directory+"/"+self.json_file_names["gw_param_detectable"])

    def class_initialization(self, params=None):
        """
        Function to initialize the parent classes. 

        Parameters
        ----------
        params : `dict`
            dictionary of parameters to initialize the parent classes
        """

        # initialization of CompactBinaryPopulation class
        # it also initializes the CBCSourceRedshiftDistribution class
        input_params = dict(
            z_min=self.z_min,
            z_max=self.z_max,
            cosmology=self.cosmo,
            event_type=self.event_type,
            source_priors=None,
            source_priors_params=None,
            spin_zero=True,
            spin_precession=False,
            directory=self.interpolator_directory,
            create_new_interpolator=False,
        )
        if params:
            for key, value in params.items():
                if key in input_params:
                    input_params[key] = value
        self.gw_param_sampler_dict = input_params
        # initialization of clasess
        CBCSourceParameterDistribution.__init__(
            self,
            z_min=input_params["z_min"],
            z_max=input_params["z_max"],
            event_type=input_params["event_type"],
            source_priors=input_params["source_priors"],
            source_priors_params=input_params["source_priors_params"],
            cosmology=input_params["cosmology"],
            spin_zero=input_params["spin_zero"],
            spin_precession=input_params["spin_precession"],
            directory=input_params["directory"],
            create_new_interpolator=input_params["create_new_interpolator"],
        )

        self.gw_param_sampler_dict["source_priors"]=self.gw_param_samplers.copy()
        self.gw_param_sampler_dict["source_priors_params"]=self.gw_param_samplers_params.copy()

    def gwsnr_intialization(self, params=None):
        """
        Function to initialize the GWSNR class from the `gwsnr` package.

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
            interpolator_dir=self.interpolator_directory,
            create_new_interpolator=False,
            gwsnr_verbose=False,
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
        self.snr_bilby = gwsnr.compute_bilby_snr
        self.snr_calculator_dict["mtot_max"] = gwsnr.mtot_max
        self.snr_calculator_dict["psds"] = gwsnr.psds_list
        #self.pdet = gwsnr.pdet

    def store_gwrates_params(self, output_jsonfile="gwrates_params.json"):
        """
        Function to store the all the necessary parameters. This is useful for reproducing the results. All the parameters stored are in string format to make it json compatible.

        Parameters
        ----------
        output_jsonfile : `str`
            name of the json file to store the parameters
        """

        # store gw_param_sampler_dict, lensed_param_sampler_dict and snr_calculator_dict
        parameters_dict = dict(
            npool=str(self.npool),
            z_min=str(self.z_min),
            z_max=str(self.z_max),
            size=str(self.size),
            batch_size=str(self.batch_size),
            cosmology=str(self.cosmo),
            snr_finder=str(self.snr),
            json_file_names=str(self.json_file_names),
            interpolator_directory=str(self.interpolator_directory),
        )

        # cbc params
        gw_param_sampler_dict = self.gw_param_sampler_dict.copy()
        # convert all dict values to str
        for key, value in gw_param_sampler_dict.items():
            gw_param_sampler_dict[key] = str(value)
        parameters_dict.update({"gw_param_sampler_dict": gw_param_sampler_dict})

        # snr calculator params
        try:
            snr_calculator_dict = self.snr_calculator_dict.copy()
            for key, value in snr_calculator_dict.items():
                snr_calculator_dict[key] = str(value)
            parameters_dict.update({"snr_calculator_dict": snr_calculator_dict})

            file_name = output_jsonfile
            append_json(self.ler_directory+"/"+file_name, parameters_dict, replace=True)
        except:
            # if snr_calculator is custom function
            pass

    def gw_cbc_statistics(
        self, size=None, resume=False, save_batch=False, output_jsonfile=None,
    ):
        """
        Function to generate gw GW source parameters. This function calls the gw_sampling_routine function to generate the parameters in batches. The generated parameters are stored in a json file; and if save_batch=True, it keeps updating the file in batches.

        Parameters
        ----------
        size : `int`
            number of samples.
            default size = 100000.
        resume : `bool`
            resume = False (default) or True.
            if True, the function will resume from the last batch.
        save_batch : `bool`
            if True, the function will save the parameters in batches. if False, the function will save all the parameters at the end of sampling. save_batch=False is faster.
        output_jsonfile : `str`
            json file name for storing the parameters.
            default output_jsonfile = 'gw_params.json'. Note that this file will be stored in the self.ler_directory.

        Returns
        ----------
        gw_param : `dict`
            dictionary of gw GW source parameters. Refer to :attr:`~gw_param` for details.
        
        Examples
        ----------
        >>> from ler.rates import GWRATES
        >>> ler = GWRATES()
        >>> param = ler.gw_cbc_statistics()
        """

        size = size or self.size
        output_jsonfile = output_jsonfile or self.json_file_names["gw_param"]
        self.json_file_names["gw_param"] = output_jsonfile
        output_path = os.path.join(self.ler_directory, output_jsonfile)
        print(f"Simulated GW params will be stored in {output_path}")

        # sampling in batches
        if resume and os.path.exists(output_path):
            # get sample from json file
            self.dict_buffer = get_param_from_json(output_path)
        else:
            self.dict_buffer = None

        batch_handler(
            size=size,
            batch_size=self.batch_size,
            sampling_routine=self.gw_sampling_routine,
            output_jsonfile=output_path,
            save_batch=save_batch,
            resume=resume,
        )

        if save_batch:
            gw_param = get_param_from_json(output_path)
        else:
            # this if condition is required if there is nothing to save
            if self.dict_buffer:
                gw_param = self.dict_buffer.copy()
                # store all params in json file
                print(f"saving all gw_params in {output_path}...")
                append_json(output_path, gw_param, replace=True)
            else:
                print("gw_params already sampled.")
                gw_param = get_param_from_json(output_path)
        self.dict_buffer = None  # save memory

        return gw_param
    
    def gw_sampling_routine(self, size, output_jsonfile, resume=False, save_batch=True):
        """
        Function to generate GW source parameters. This function also stores the parameters in json file in the current batch if save_batch=True.

        Parameters
        ----------
        size : `int`
            number of samples.
            default size = 100000.
        output_jsonfile : `str`
            json file name for storing the parameters.
            default output_jsonfile = 'gw_params.json'. Note that this file will be stored in the self.ler_directory.
        resume : `bool`
            resume = False (default) or True. 
            if True, it appends the new samples to the existing json file.
        save_batch : `bool`
            if True, the function will save the parameters in batches. if False, the function will save all the parameters at the end of sampling. save_batch=False is faster.

        Returns
        ----------
        gw_param : `dict`
            dictionary of gw GW source parameters. Refer to :attr:`~gw_param` for details.
        """

        # get gw params
        print("sampling gw source params...")
        gw_param = self.sample_gw_parameters(size=size)
        # Get all of the signal to noise ratios
        if self.snr:
            print("calculating snrs...")
            snrs = self.snr(gw_param_dict=gw_param)
            gw_param.update(snrs)
        elif self.pdet:
            print("calculating pdet...")
            pdet = self.pdet(gw_param_dict=gw_param)
            gw_param.update(pdet)

        # adding batches
        if not save_batch:
            if self.dict_buffer is None:
                self.dict_buffer = gw_param
            else:
                for key, value in gw_param.items():
                    self.dict_buffer[key] = np.concatenate((self.dict_buffer[key], value))
        else:
            # store all params in json file
            self.dict_buffer = append_json(file_name=output_jsonfile, new_dictionary=gw_param,  old_dictionary=self.dict_buffer, replace=not (resume))

        return gw_param

    def gw_rate(
        self,
        gw_param=None,
        snr_threshold=8.0,
        pdet_threshold=0.5,
        output_jsonfile=None,
        detectability_condition="step_function",
        snr_recalculation=False,
        snr_threshold_recalculation=[4, 20],
    ):
        """
        Function to calculate the GW rate. This function also stores the parameters of the detectable events in json file. There are two conditions for detectability: 'step_function' and 'pdet'.

        1. 'step_function': If two images have SNR>8.0, then the event is detectable. This is a step function. This is with the assumption that SNR function is provided and not None. 
        2. 'pdet':
            i) If self.pdet is None and self.snr is not None, then it will calculate the pdet from the snr. There is no hard cut for this pdet and can have value ranging from 0 to 1 near the threshold.
            ii) If self.pdet is not None, then it will use the generated pdet.

        Parameters
        ----------
        gw_param : `dict` or `str`
            dictionary of GW source parameters or json file name.
            default gw_param = self.json_file_names["gw_param"]
        snr_threshold : `float`
            threshold for detection signal to noise ratio.
            e.g. snr_threshold = 8.
        pdet_threshold : `float`
            threshold for detection probability.
            e.g. pdet_threshold = 0.5.
        output_jsonfile : `str`
            json file name for storing the parameters of the detectable events.
            default output_jsonfile = 'gw_params_detectable.json'.
        detectability_condition : `str`
            detectability condition. 
            default detectability_condition = 'step_function'.
            other options are 'pdet'.
        snr_recalculation : `bool`
            if True, the SNR of centain events (snr>snr_threshold_recalculation)will be recalculate with 'inner-product' method. This is useful when the snr is calculated with 'ann' method.
            default snr_recalculation = False.
        snr_threshold_recalculation : `list`
            lower and upper threshold for recalculation of detection signal to noise ratio.
            default snr_threshold_recalculation = [4, 20].

        Returns
        ----------
        total_rate : `float`
            total gw rate (Mpc^-3 yr^-1).
        gw_param : `dict`
            dictionary of gw GW source parameters of the detectable events. Refer to :attr:`~gw_param` for details.

        Examples
        ----------
        >>> from ler.rates import GWRATES
        >>> ler = GWRATES()
        >>> ler.gw_cbc_statistics();    
        >>> total_rate, gw_param = ler.gw_rate()
        """
        
        gw_param = self._load_param(gw_param)
        total_events = len(gw_param["zs"])

        # below is use when the snr is calculated with 'ann' method of `gwsnr`
        if snr_recalculation:
            gw_param = self._recalculate_snr(gw_param, snr_threshold_recalculation)

        # find index of detectable events
        idx_detectable = self._find_detectable_index(gw_param, snr_threshold, pdet_threshold, detectability_condition)

        detectable_events = np.sum(idx_detectable)
        # montecarlo integration
        # The total rate R = norm <Theta(rho-rhoc)>
        total_rate = self.rate_function(detectable_events, total_events)

        # store all detectable params in json file
        self._save_detectable_params(output_jsonfile, gw_param, idx_detectable, key_file_name="gw_param_detectable", nan_to_num=False, verbose=True, replace_jsonfile=True)

        # append ler_param and save it
        self._append_ler_param(total_rate, detectability_condition)
        
        return total_rate, gw_param

    def _load_param(self, param):
        """
        Helper function to load or copy GW parameters.
        
        Parameters
        ----------
        param : `dict` or `str`
            dictionary of GW parameters or json file name.

        Returns
        ----------
        param : `dict`
            dictionary of GW parameters.
        """

        if param is None:
            param = self.json_file_names["gw_param"]
        if isinstance(param, str):
            path_ = self.ler_directory + "/" + param
            print(f"Getting GW parameters from json file {path_}...")
            return get_param_from_json(path_)
        else:
            print("Using provided {param_type} dict...")
            return param.copy()

    def _recalculate_snr(self, gw_param, snr_threshold_recalculation):
        """
        Recalculates SNR for events where the initial SNR is above a given threshold.
        
        Parameters
        ---------- 
        gw_param : `dict`
            dictionary of GW source parameters.
        snr_threshold_recalculation : `list`
            lower and upper threshold for recalculation of detection signal to noise ratio.
            default snr_threshold_recalculation = [4, 20].

        Returns
        ----------
        gw_param : `dict`
            dictionary of GW source parameters.
        """

        snr_param = gw_param["optimal_snr_net"]
        idx_detectable = (snr_param > snr_threshold_recalculation[0]) & (snr_param < snr_threshold_recalculation[1])
        # reduce the size of the dict
        for key, value in gw_param.items():
            gw_param[key] = value[idx_detectable]
        # recalculate more accurate snrs 
        snrs = self.snr_bilby(gw_param_dict=gw_param)
        gw_param.update(snrs)
        return gw_param

    def _find_detectable_index(self, gw_param, snr_threshold, pdet_threshold, detectability_condition):
        """
        Find the index of detectable events based on SNR or p_det.
        
        Parameters
        ----------
        gw_param : `dict`
            dictionary of GW source parameters.
        snr_threshold : `float`
            threshold for detection signal to noise ratio.
        pdet_threshold : `float`
            threshold for detection probability.
        detectability_condition : `str`
            detectability condition.
            default detectability_condition = 'step_function'.
            other options are 'pdet'.

        Returns
        ----------
        idx_detectable : `numpy.ndarray`
            index of detectable events.
        """

        if self.snr:
            if "optimal_snr_net" not in gw_param:
                raise ValueError("'optimal_snr_net' not in gw param dict provided")
            if detectability_condition == "step_function":
                print("given detectability_condition == 'step_function'")
                param = gw_param["optimal_snr_net"]
                threshold = snr_threshold
            elif detectability_condition == "pdet":
                print("given detectability_condition == 'pdet'")
                param = 1 - norm.cdf(snr_threshold - gw_param["optimal_snr_net"])
                gw_param["pdet_net"] = param
                threshold = pdet_threshold
        elif self.pdet:
            if "pdet_net" in gw_param:
                print("given detectability_condition == 'pdet'")
                param = gw_param["pdet_net"]
                threshold = pdet_threshold
            else:
                raise ValueError("'pdet_net' not in gw param dict provided")

        idx_detectable = param > threshold
        return idx_detectable

    def rate_function(self, detectable_size, total_size, verbose=True):
        """
        General helper function to calculate the rate for GW events.

        Parameters
        ----------
        detectable_size : `int`
            number of detectable events.
        total_size : `int`
            total number of events.
        param_type : `str`
            type of parameters.

        Returns
        ----------
        rate : `float`
            rate of the events.

        Examples
        ----------
        >>> from ler.rates import LeR
        >>> ler = LeR()
        >>> rate = ler.rate_function(detectable_size=100, total_size=1000)
        """

        normalization = self.normalization_pdf_z
        rate = normalization * detectable_size / total_size

        if verbose:
            print(f"total GW event rate (yr^-1): {rate}")
            print(f"number of simulated GW detectable events: {detectable_size}")
            print(f"number of simulated all GW events: {total_size}")

        return rate

    def _save_detectable_params(self,
        output_jsonfile,
        param, 
        idx_detectable, 
        key_file_name="gw_param_detectable",
        nan_to_num=False,
        verbose=True,
        replace_jsonfile=True,
    ):
        """
        Helper function to save the detectable parameters in json file.

        Parameters
        ----------
        output_jsonfile : `str`
            json file name for storing the parameters of the detectable events. This is stored in the self.ler_directory.
        param : `dict`
            dictionary of GW source parameters.
        idx_detectable : `numpy.ndarray`
            index of detectable events.
        key_file_name : `str`
            key name for the json file to be added in self.json_file_names.
        nan_to_num : `bool`
            if True, it will replace nan with 0.
            default nan_to_num = False.
        verbose : `bool`
            if True, it will print the path of the json file.
            default verbose = True.
        replace_jsonfile : `bool`
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

        output_path = self.ler_directory+"/"+output_jsonfile
        if verbose:
            print(f"storing detectable params in {output_path}")
        append_json(output_path, param, replace=replace_jsonfile)

    def _append_ler_param(self, total_rate, detectability_condition):
        """
        Helper function to append the final results, total_rate, in the json file.

        Parameters
        ----------
        total_rate : `float`
            total rate.
        detectability_condition : `str`
            detectability condition.
        """

        data = load_json(self.ler_directory+"/"+self.json_file_names["gwrates_params"])
        # write the results
        data["detectable_gw_rate_per_year"] = total_rate
        data["detectability_condition"] = detectability_condition
        append_json(self.ler_directory+"/"+self.json_file_names["gwrates_params"], data, replace=True)

    def selecting_n_gw_detectable_events(
        self,
        size=100,
        batch_size=None,
        snr_threshold=8.0,
        pdet_threshold=0.5,
        resume=False,
        output_jsonfile="gw_params_n_detectable.json",
        meta_data_file="meta_gw.json",
        detectability_condition="step_function",
        trim_to_size=True,
        snr_recalculation=False,
        snr_threshold_recalculation=[4, 12],
    ):
        """
        Function to generate n GW detectable events. This fuction samples the GW parameters and save only the detectable events in json file. It also records metadata in the JSON file, which includes the total number of events and the cumulative rate of events. This functionality is particularly useful for generating a fixed or large number of detectable events until the event rates stabilize.

        Parameters
        ----------
        size : `int`
            number of samples to be selected.
            default size = 100.
        batch_size : `int`
            batch size for sampling.
            default batch_size = 50000.
        snr_threshold : `float`
            threshold for detection signal to noise ratio.
            e.g. snr_threshold = 8.
        pdet_threshold : `float`
            threshold for detection probability.
            default pdet_threshold = 0.5.
        resume : `bool`
            resume = False (default) or True.
            if True, the function will resume from the last batch.
        output_jsonfile : `str`
            json file name for storing the parameters of the detectable events.
            default output_jsonfile = 'n_gw_param_detectable.json'.
        meta_data_file : `str`
            json file name for storing the metadata.
            default meta_data_file = 'meta_gw.json'.
        detectability_condition : `str`
            detectability condition.
            default detectability_condition = 'step_function'.
            other options are 'pdet'.
        trim_to_size : `bool`
            if True, the final result will be trimmed to size.
            default trim_to_size = True.
        snr_recalculation : `bool`
            if True, the SNR of centain events (snr>snr_threshold_recalculation)will be recalculate with 'inner-product' method. This is useful when the snr is calculated with 'ann' method of `gwsnr`.
            default snr_recalculation = False.
        snr_threshold_recalculation : `list`
            lower and upper threshold for recalculation of detection signal to noise ratio.
            default snr_threshold_recalculation = [4, 12].

        Returns
        ----------
        param_final : `dict`
            dictionary of gw GW source parameters of the detectable events. Refer to :attr:`~gw_param` for details.

        Examples
        ----------
        >>> from ler.rates import LeR
        >>> ler = LeR()
        >>> gw_param = ler.selecting_n_gw_detectable_events(size=100)
        """

        # initial setup
        n, events_total, output_path, meta_data_path, buffer_file = self._initial_setup_for_n_event_selection(meta_data_file, output_jsonfile, resume, batch_size)

        # loop until n samples are collected
        while n < size:
            # disable print statements
            with contextlib.redirect_stdout(None):
                self.dict_buffer = None  # this is used to store the sampled gw_param in batches when running the sampling_routine
                gw_param = self.gw_sampling_routine(
                    size=batch_size, output_jsonfile=buffer_file, save_batch=False,resume=False
                )

            total_events_in_this_iteration = len(gw_param["zs"])
            # below is use when the snr is calculated with 'ann' method of `gwsnr`
            if snr_recalculation:
                # select only above centain snr threshold
                gw_param = self._recalculate_snr(gw_param, snr_threshold_recalculation)

            # find index of detectable events
            idx_detectable = self._find_detectable_index(gw_param, snr_threshold, pdet_threshold, detectability_condition)

            # store all params in json file
            self._save_detectable_params(output_jsonfile, gw_param, idx_detectable, key_file_name="n_gw_detectable_events", nan_to_num=False, verbose=False, replace_jsonfile=False)

            n += np.sum(idx_detectable)
            events_total += total_events_in_this_iteration              
            total_rate = self.rate_function(n, events_total, verbose=False)

            # bookmark
            self._append_meta_data(meta_data_path, n, events_total, total_rate)

        print(f"stored detectable gw params in {output_path}")
        print(f"stored meta data in {meta_data_path}")

        if trim_to_size:
            param_final, total_rate = self._trim_results_to_size(size, output_path, meta_data_path)
        else:
            param_final = get_param_from_json(output_path)

        # call self.json_file_names["ler_param"] and for adding the final results
        data = load_json(self.ler_directory+"/"+self.json_file_names["gwrates_params"])
        # write the results
        try:
            data["detectable_gw_rate_per_year"] = total_rate
            data["detectability_condition"] = detectability_condition
        except:
            meta = get_param_from_json(meta_data_path)
            data["detectable_gw_rate_per_year"] = meta["total_rate"][-1]
            data["detectability_condition"] = detectability_condition

        append_json(self.ler_directory+"/"+self.json_file_names["gwrates_params"], data, replace=True)

        return param_final


    def _initial_setup_for_n_event_selection(self, meta_data_file, output_jsonfile, resume, batch_size):
        """Helper function for selecting_n_gw_detectable_events and selecting_n_lensed_detectable_events functions. 

        Parameters
        ----------
        meta_data_file : `str`
            json file name for storing the metadata.
        output_jsonfile : `str`
            json file name for storing the parameters of the detectable events.
        resume : `bool`
            resume = False (default) or True.
            if True, the function will resume from the last batch.
        batch_size : `int`
            batch size for sampling.
            default batch_size = 50000.

        Returns
        ----------
        n : `int`
            iterator.
        events_total : `int`
            total number of events.
        output_path : `str`
            path to the output json file.
        meta_data_path : `str`
            path to the metadata json file.
        buffer_file : `str`
            path to the buffer json file.
        """

        meta_data_path = self.ler_directory+"/"+meta_data_file
        output_path = self.ler_directory+"/"+output_jsonfile
        if meta_data_path==output_path:
            raise ValueError("meta_data_file and output_jsonfile cannot be same.")
            
        if batch_size is None:
            batch_size = self.batch_size
        else:
            self.batch_size = batch_size

        if not resume:
            n = 0  # iterator
            events_total = 0
            if os.path.exists(output_path):
                os.remove(output_path)
            if os.path.exists(meta_data_path):
                os.remove(meta_data_path)
        else:
            # get sample size as size from json file
            if os.path.exists(output_path):
                param_final = get_param_from_json(output_path)
                n = len(param_final["zs"])
                events_total = load_json(meta_data_path)["events_total"][-1]
                del param_final
            else:
                n = 0
                events_total = 0

        buffer_file = "params_buffer.json"
        print("collected number of detectable events = ", n)

        return n, events_total, output_path, meta_data_path, buffer_file

    def _trim_results_to_size(self, size, output_path, meta_data_path):
        """
        Helper function of 'selecting_n_gw_detectable_events' and 'selecting_n_lensed_detectable_events' functions. Trims the data in the output file to the specified size and updates the metadata accordingly.

        Parameters
        ----------
        size : `int`
            number of samples to be selected.
        output_path : `str`
            path to the output json file.
        meta_data_path : `str`
            path to the metadata json file.

        Returns
        ----------
        param_final : `dict`
            dictionary of GW source parameters of the detectable events. Refer to :attr:`~gw_param`
        new_total_rate : `float`
            total rate (Mpc^-3 yr^-1).
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
        new_events_total = np.round(size*old_events_total/old_detectable_events)
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
        Helper function for appending meta data json file.

        Parameters
        ----------
        meta_data_path : `str`
            path to the metadata json file.
        n : `int`
            iterator.
        events_total : `int`
            total number of events.
        total_rate : `float`
            total rate (Mpc^-3 yr^-1).
        """

        # save meta data
        meta_data = dict(events_total=[events_total], detectable_events=[float(n)], total_rate=[total_rate])

        if os.path.exists(meta_data_path):
            try:
                append_json(meta_data_path, meta_data, replace=False)
            except:
                append_json(meta_data_path, meta_data, replace=True)
        else:
            append_json(meta_data_path, meta_data, replace=True)

        print("collected number of detectable events = ", n)
        print("total number of events = ", events_total)
        print(f"total rate (yr^-1): {total_rate}")