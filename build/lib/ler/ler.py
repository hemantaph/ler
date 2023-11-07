# -*- coding: utf-8 -*-
"""
This module contains the main class for calculating the rates of lensed and unlensed events.
"""

import os
import json
import random
import contextlib
import numpy as np
from gwsnr import GWSNR
from scipy.stats import norm, gaussian_kde
from scipy.interpolate import interp1d
from astropy.cosmology import Planck18
from ler.lens_galaxy_population import LensGalaxyPopulation
from ler.source_population import CompactBinaryPopulation
from ler.helperroutines import append_json, get_param_from_json

# Conversions from SI units to CGS units
C = 299792458.0  # m/s
G = 6.67408 * 1e-11  # m^3/kg/s^2


class LeR:
    """Class to calculate both the rates of lensed and unlensed events.

    Parameters
    ----------
    nsamples : `int`
        number of samples for sampling.
        default nsamples = 100000.
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
    batch_size : `int`
        batch size for SNR calculation.
        default batch_size = 25000.
        reduce the batch size if you are getting memory error.
    snr_finder : `str`
        default snr_finder = 'gwsnr'.
        if 'gwsnr', the SNR will be calculated using the gwsnr package.
        if 'custom', the SNR will be calculated using a custom function.
    json_file_ler_param: `str`
        default json_file_ler_param = 'ler_param.json'.
        json file containing the parameters for initializing the :class:`~ler.LeR` class, :class:`~ler.CompactBinaryPopulation` class, :class:`~ler.LensGalaxyPopulation` class, :class:`~gwsnr.GWSNR` class.
    kwargs : `keyword arguments`
        Note : kwargs takes input for initializing the :class:`~ler.CompactBinaryPopulation`, :class:`LensGalaxyPopulation`, :meth:`~gwsnr_intialization`.

    Examples
    ----------
    - class initialization
    - ``ler`` needs `gwsnr <https://github.com/hemantaph/gwsnr/>`_.
    - generation of ``gwsnr`` snr interpolator will take time at the first initialization. The interpolator will be stored in the working dir.
    - ``m_min``, ``m_max`` were used for initializing the ``CompactBinaryPopulation`` class. ``waveform_approximant`` was used for initializing the ``snr_calculator`` (``gwsnr``) class. ``min_lensed_images`` was used for initializing the ``LensGalaxyPopulation`` class.

    >>> from ler import LeR
    >>> ler_ = LeR(nsamples=100000, npool=int(4), z_min=0., z_max=10., batch_size=25000, snr_finder='gwsnr', m_min=4.59, m_max=86.22, waveform_approximant='IMRPhenomD', min_lensed_images=2)
    Given: IMR waveform
    psds not given. Choosing bilby's default psds
    getting stored interpolator...
    In case if you need regeneration of interpolator of the given gwsnr param, please delete this file, ./interpolator_pickle/halfSNR_dict_0.pickle

    Instance Attributes
    ----------
    LeR class has the following attributes, \n
    +-------------------------------------+----------------------------------+
    | Atrributes                          | Type                             |
    +=====================================+==================================+
    |:attr:`~gw_param`                    |`dict`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~gw_param_detectable`         |`dict`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~lensed_param`                |`dict`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~lensed_param_detectable`     |`dict`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~gw_param_sampler_dict`       |`dict`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~lensed_param_sampler_dict`   |`dict`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~snr_calculator_dict`         |`dict`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~z_to_Dc`                     |`scipy.interpolate.interp1d`      |
    +-------------------------------------+----------------------------------+
    |:attr:`~Dc_to_z`                     |`scipy.interpolate.interp1d`      |
    +-------------------------------------+----------------------------------+
    |:attr:`~z_to_luminosity_distance`    |`scipy.interpolate.interp1d`      |
    +-------------------------------------+----------------------------------+
    |:attr:`~differential_comoving_volume`|`scipy.interpolate.interp1d`      |
    +-------------------------------------+----------------------------------+
    |:attr:`~compact_binary_pop`          |`CompactBinaryPopulation class`   |
    +-------------------------------------+----------------------------------+
    |:attr:`~lens_galaxy_pop`             |`LensGalaxyPopulation class`      |
    +-------------------------------------+----------------------------------+
    | :attr:`~snr`                        |``gwsnr`` `package`               |
    +-------------------------------------+----------------------------------+

    Instance Methods
    ----------
    LeR class has the following method(s), \n
    +------------------------------------+-------------------------------------+
    | Method(s)                          | Description                         |
    +====================================+=====================================+
    |:meth:`~gwsnr_intialization`        |Function for initializing the        |
    |                                    |``gwsnr`` package.                   |
    +------------------------------------+-------------------------------------+
    |:meth:`~create_lookup_tables`       |To creating lookup tables for fast   |
    |                                    |calculation for the following        |
    |                                    |conversion operations,               |
    |                                    |redshift to co-moving distance.      |
    |                                    |co-moving distance to redshift.      |
    |                                    |redshift to luminosity distance.     |
    +------------------------------------+-------------------------------------+
    |:meth:`~unlensed_cbc_statistics`    |Function to generate unlensed GW     |
    |                                    |source parameters.                   |
    +------------------------------------+-------------------------------------+
    |:meth:`~unlensed_rate`              |Function to calculate unlensed       |
    |                                    |merger rate.                         |
    +------------------------------------+-------------------------------------+
    |:meth:`~lensed_cbc_statistics`      |Function to generate lensed GW       |
    |                                    |source parameters.                   |
    +------------------------------------+-------------------------------------+
    |:meth:`~lensed_rate`                |Function to calculate lensed         |
    |                                    |merger rate.                         |
    +------------------------------------+-------------------------------------+
    |:meth:`~batch_handler`              |Function to handle the batch size.   |
    +------------------------------------+-------------------------------------+
    |:meth:`~store_ler_params`           |Fuction to store the parameters of   |
    |                                    |the LER model.                       |
    +------------------------------------+-------------------------------------+

    """

    # Attributes
    gw_param_sampler_dict = None
    """``dict`` \n
    dictionary of params for initializing ``CompactBinaryPopulation`` class \n
    this will be used for GW unlensed parameters sampling \n
    gw_param_sampler_dict.keys() = ['nsamples', 'm_min', 'm_max', 'z_min', 'z_max', 'event_type', 'src_model_params']
    """

    lensed_param_sampler_dict = None
    """``dict`` \n
    dictionary of params for initializing ``LensGalaxyPopulation`` class \n
    this will be used for GW lensed parameters sampling \n
    lensed_param_sampler_dict.keys() = ['nsamples', 'min_lensed_images', 'max_lensed_images', 'lensModelList']
    """

    snr_calculator_dict = None
    """``dict`` \n
    dictionary of params for initializing ``snr_calculator`` (``gwsnr``) class \n
    this will be used for SNR calculation \n
    snr_calculator_dict.keys() = ['mtot_min', 'mtot_max', 'nsamples_mtot', 'nsamples_mass_ratio', 'sampling_frequency', 'waveform_approximant', 'minimum_frequency', 'snr_type', 'waveform_inspiral_must_be_above_fmin', 'psds', 'psd_file', 'ifos']
    """

    z_to_Dc = None
    """``scipy.interpolate.interp1d`` \n
    redshift to co-moving distance.
    """

    Dc_to_z = None
    """``scipy.interpolate.interp1d`` \n
    co-moving distance to redshift.
    """

    z_to_luminosity_distance = None
    """``scipy.interpolate.interp1d`` \n
    redshift to luminosity distance.
    """

    differential_comoving_volume = None
    """``scipy.interpolate.interp1d`` \n
    differential comoving volume.
    """

    compact_binary_pop = None
    """``CompactBinaryPopulation class`` \n
    class for sampling GW parameters.
    """

    lens_galaxy_pop = None
    """``LensGalaxyPopulation class`` \n
    class for sampling lensed GW parameters.
    """

    snr = None
    """``gwsnr package`` \n
    class for calculating SNR.
    """

    def __init__(
        self,
        nsamples=100000,
        npool=int(4),
        z_min=0.0001,
        z_max=10.0,
        batch_size=25000,
        snr_finder="gwsnr",
        json_file_ler_param="./LeR_params.json",
        **kwargs,
    ):
        self.z_min = z_min
        self.z_max = z_max
        self.npool = npool
        # batch size for parameters sampling and snr calculation
        # try keeping it below 50000
        # reduce the batch size if you are getting memory error
        self.batch_size = batch_size
        self.gw_param = None
        self.gw_param_detectable = None
        self.lensed_param = None
        self.lensed_param_detectable = None
        self.json_file_ler_param = json_file_ler_param

        # dictionary of params for sampler
        # for unlened case (source params)
        # defualt for 'src_model_params' is set for popI_II PowerLaw+PEAK model
        # for other models, please change the 'src_model_params' accordingly
        self.gw_param_sampler_dict = {
            "nsamples": nsamples,
            "m_min": 4.59,
            "m_max": 86.22,
            "z_min": z_min,
            "z_max": z_max,
            "event_type": "BBH",
            "category": "popI_II",
            "sub_category": "gwcosmo",
            "redshift_event_type": None,
            "redshift_category": None,
            "merger_rate_density_fn": None,
            "merger_rate_density_param": None,
            "src_model_params": None,
            "mass_constant": False,
            "redshift_constant": False,
            "spin_constant": False,
        }
        # for lensed case
        # set 'min_lensed_images' = 2 for double image lensed case
        self.lensed_param_sampler_dict = {
            "nsamples": nsamples,
            "min_lensed_images": 2,
            "max_lensed_images": 4,
            "lensModelList": ["EPL_NUMBA", "SHEAR"],
        }

        # for snr_calculator
        # for 'waveform_approximant' other than IMRPhenomD or TaylorF2, please set 'snr_type' = 'inner_product'
        # you will get accurate results if you set 'nsamples_mtot': 200, 'nsamples_mass_ratio': 500., 'sampling_frequency': 4096.
        self.snr_calculator_dict = {
            "mtot_min": 2.0,
            "mtot_max": 439.6,
            "nsamples_mtot": 100,
            "nsamples_mass_ratio": 50,
            "sampling_frequency": 2048.0,
            "waveform_approximant": "IMRPhenomD",
            "minimum_frequency": 20.0,
            "snr_type": "interpolation",
            "waveform_inspiral_must_be_above_fmin": False,
            "psds": None,
            "psd_file": False,
            "ifos": None,
            "interpolator_dir": "./interpolator_pickle",
        }

        # update dict from kwargs
        keys1 = self.gw_param_sampler_dict.keys()
        keys2 = self.lensed_param_sampler_dict.keys()
        keys3 = self.snr_calculator_dict.keys()
        for key, value in kwargs.items():
            if key in keys1:
                self.gw_param_sampler_dict[key] = value
            if key in keys2:
                self.lensed_param_sampler_dict[key] = value
            if key in keys3:
                self.snr_calculator_dict[key] = value

        # initialization of clasess, for both  CompactBinaryPopulation and LensGalaxyPopulation
        # CompactBinaryPopulation already inherits from Source_Galaxy_Population_Model class form source_population.py
        # LensGalaxyPopulation already inherits from Lens_Galaxy_Population_Model class form lens_galaxy_population.py
        self.class_initialization()

        # initializing function for fast SNR calculation (in case of gwsnr)
        if snr_finder == "gwsnr":
            # default
            self.snr = self.gwsnr_intialization(kwargs)
        else:
            # custom SNR function
            print(f"Input your custom SNR finder function.\n self.snr=custom(gw_param_dict) .\n For details see the documentation.")
            self.snr = None
        # extra note on how to change snr finder function
        # self.snr = custom_snr_finder_function()

        # Create lookup tables
        self.create_lookup_tables(z_min, z_max)

        self.store_ler_params(json_file=json_file_ler_param)

        return None

    # Attributes
    @property
    def gw_param(self):
        """``bool``, ``dict`` \n
        gw_param is a dictionary of unlensed parameters (source parameters) \n
        it will be populated when unlened_cbc_statistics() is called \n
        if unavailable, the unlensed parameters will be sampled when unlensed_rate() is called \n
        gw_param.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time'] \n
        """
        # if file name
        if isinstance(self._gw_param, str):
            f = open(self._gw_param, "r", encoding="utf-8")
            self._gw_param = json.loads(f.read())
        return self._gw_param

    @gw_param.setter
    def gw_param(self, value):
        self._gw_param = value

    @property
    def gw_param_detectable(self):
        """``bool``, ``dict`` \n
        gw_param_detectable is a dictionary of unlensed parameters (source parameters) \n
        it will be populated when unlened_cbc_statistics() is called \n
        if unavailable, the unlensed parameters will be sampled when unlensed_rate() is called \n
        gw_param_detectable.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time'] \n
        """
        # if file name
        if isinstance(self._gw_param_detectable, str):
            f = open(self._gw_param_detectable, "r", encoding="utf-8")
            self._gw_param_detectable = json.loads(f.read())
        return self._gw_param_detectable

    @gw_param_detectable.setter
    def gw_param_detectable(self, value):
        self._gw_param_detectable = value

    @property
    def lensed_param(self):
        """``bool``, ``dict`` \n
        lensed_param is a dictionary of lensed parameters \n
        it will be populated when lensed_cbc_statistics() is called \n
        if unavailable, the lensed parameters will be sampled when lensed_rate() is called \n
        lensed_param.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time', 'lensed_images'] \n
        """
        # if file name
        if isinstance(self._lensed_param, str):
            f = open(self._lensed_param, "r", encoding="utf-8")
            self._lensed_param = json.loads(f.read())
        return self._lensed_param

    @lensed_param.setter
    def lensed_param(self, value):
        self._lensed_param = value

    @property
    def lensed_param_detectable(self):
        """``bool``, ``dict`` \n
        lensed_param_detectable is a dictionary of lensed parameters \n
        it will be populated when lensed_cbc_statistics() is called \n
        if unavailable, the lensed parameters will be sampled when lensed_rate() is called \n
        lensed_param_detectable.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time', 'lensed_images'] \n
        """
        # if file name
        if isinstance(self._lensed_param_detectable, str):
            f = open(self._lensed_param_detectable, "r", encoding="utf-8")
            self._lensed_param_detectable = json.loads(f.read())
        return self._lensed_param_detectable

    @lensed_param_detectable.setter
    def lensed_param_detectable(self, value):
        self._lensed_param_detectable = value

    # CompactBinaryPopulation class and LensGalaxyPopulation class initialization
    def class_initialization(self):
        """
        Function for initializing the ``CompactBinaryPopulation`` and ``LensGalaxyPopulation`` classes.
        """
        self.compact_binary_pop = CompactBinaryPopulation(
            z_min=self.z_min,
            z_max=self.z_max,
            m_min=self.gw_param_sampler_dict["m_min"],
            m_max=self.gw_param_sampler_dict["m_max"],
            event_type=self.gw_param_sampler_dict["event_type"],
            category=self.gw_param_sampler_dict["category"],
            sub_category=self.gw_param_sampler_dict["sub_category"],
            redshift_event_type=self.gw_param_sampler_dict["redshift_event_type"],
            redshift_category=self.gw_param_sampler_dict["redshift_category"],
            merger_rate_density_fn=self.gw_param_sampler_dict[
                "merger_rate_density_fn"],
            merger_rate_density_param=self.gw_param_sampler_dict[
                "merger_rate_density_param"],
            src_model_params=self.gw_param_sampler_dict["src_model_params"],
            mass_constant=self.gw_param_sampler_dict["mass_constant"],
            redshift_constant=self.gw_param_sampler_dict["redshift_constant"],
            spin_constant=self.gw_param_sampler_dict["spin_constant"],
        )
        self.lens_galaxy_pop = LensGalaxyPopulation(self.compact_binary_pop)

        return None

    def store_ler_params(self, json_file="./LeR_params.json"):
        """
        Fuction to store the parameters of the LER model. This is useful for reproducing the results.
        """
        # store gw_param_sampler_dict, lensed_param_sampler_dict and snr_calculator_dict
        parameters_dict = {}

        # cbc params
        gw_param_sampler_dict = self.gw_param_sampler_dict.copy()
        gw_param_sampler_dict["merger_rate_density_fn"] = str(
            gw_param_sampler_dict["merger_rate_density_fn"]
        )
        parameters_dict.update({"gw_param_sampler_dict": gw_param_sampler_dict})

        # lensed params
        parameters_dict.update(
            {"lensed_param_sampler_dict": self.lensed_param_sampler_dict}
        )

        # snr calculator params
        snr_calculator_dict = self.snr_calculator_dict.copy()
        snr_calculator_dict["ifos"] = str(snr_calculator_dict["ifos"])
        parameters_dict.update({"snr_calculator_dict": snr_calculator_dict})

        file_name = json_file
        append_json(file_name, parameters_dict, replace=True)

        return None

    def gwsnr_intialization(self, kwargs_dict):
        """
        Function for initializing the `gwsnr <https://github.com/hemantaph/gwsnr/>`_ package.

        Parameters
        ----------
        kwargs_dict : 'dict'
            keyword arguments for the initialization of the `gwsnr` package.
            kwargs_dict.keys() \n
            ``nsamples_mtot`` : `int`
                nsamples_mtot = 200 (recommended for accurate results)
            ``nsamples_mass_ratio`` : `int`
                nsamples_mass_ratio = 500 (recommended for accurate results)
            ``sampling_frequency`` : `float`
                sampling_frequency = 4096. (recommended for accurate results)
            ``waveform_approximant`` : `str`
                waveform_approximant = "IMRPhenomD" (for BBH) or "TaylorF2" (for BNS)
                if you want to use other approximants, please set ``snr_type`` = 'inner_product'
            ``minimum_frequency`` : `float`
                minimum_frequency = 20. (for O3 and O4 runs) or 10. (for 3G detectors)
            ``snr_type`` : `str`
                snr_type = 'interpolation' (for fast results) or 'inner_product' (for bilby like results)
            ``waveform_inspiral_must_be_above_fmin`` : `bool`
                False if dont want minimum frequency cut-off as higher mass BBH can merger below that frequency.
            ``psds`` : `bool` or `dict` or `str` (txt file)
                e.g. For O4 design sensitivity \n
                    psds = {'L1':'aLIGOaLIGODesignSensitivityT1800044',\n
                    'H1':'aLIGOaLIGODesignSensitivityT1800044',\n
                    'V1':'AdvVirgo'}
            ``psd_file`` : `bool`, `list`
                psd_file = False (if ASD) or True (if PSD file)
                psd_file = [False,True] if psds[0] is a asd and psds[1] is a psd
            ``ifos`` : `list`
                interferometer object name list
                ifos = ['L1', 'H1', 'V1'] (for O4 design sensitivity)

        Returns
        ----------
        snr_ : `the gwsnr object`
            gwsnr object is used to calculate the SNR and pdet (probability of detection)
        """
        gwsnr_param_dict = self.snr_calculator_dict

        # update initialization dict from kwargs
        keys_ = gwsnr_param_dict.keys()
        for key, value in kwargs_dict.items():
            if key in keys_:
                gwsnr_param_dict[key] = value

        snr_ = GWSNR(
            npool=self.npool,
            mtot_min=gwsnr_param_dict["mtot_min"],
            mtot_max=gwsnr_param_dict["mtot_max"],
            nsamples_mtot=gwsnr_param_dict["nsamples_mtot"],
            nsamples_mass_ratio=gwsnr_param_dict["nsamples_mass_ratio"],
            sampling_frequency=gwsnr_param_dict["sampling_frequency"],
            waveform_approximant=gwsnr_param_dict["waveform_approximant"],
            minimum_frequency=gwsnr_param_dict["minimum_frequency"],
            snr_type=gwsnr_param_dict["snr_type"],
            waveform_inspiral_must_be_above_fmin=gwsnr_param_dict[
                "waveform_inspiral_must_be_above_fmin"
            ],
            psds=gwsnr_param_dict["psds"],
            psd_file=gwsnr_param_dict["psd_file"],
            ifos=gwsnr_param_dict["ifos"],
            interpolator_dir=gwsnr_param_dict["interpolator_dir"],
        )

        return snr_

    def create_lookup_tables(self, z_min, z_max):
        """
        To creating lookup tables for fast calculation for the following conversion operations,

        #. redshift to co-moving distance.
        #. co-moving distance to redshift.
        #. redshift to luminosity distance.

        Parameters
        ----------
        z_min : `float`
            minimum redshift.
            for popI_II, popIII, primordial, BNS z_min = 0., 5., 5., 0. respectively.
        z_max : `float`
            maximum redshift.
            for popI_II, popIII, primordial, BNS z_max = 10., 40., 40., 2. respectively.

        Attributes
        ----------
        z_to_Dc : `scipy.interpolate.interp1d`
            redshift to co-moving distance.
        Dc_to_z : `scipy.interpolate.interp1d`
            co-moving distance to redshift.
        z_to_luminosity_distance : `scipy.interpolate.interp1d`
            redshift to luminosity distance.
        differential_comoving_volume : `scipy.interpolate.interp1d`
            differential comoving volume.
        """

        # initialing cosmological functions for fast calculation through interpolation
        z = np.linspace(z_min, z_max, 500)  # red-shifts
        # co-moving distance in Mpc
        Dc = Planck18.comoving_distance(z).value
        # luminosity distance in Mpc
        luminosity_distance = Planck18.luminosity_distance(z).value

        # generating interpolators
        self.z_to_Dc = interp1d(z, Dc, kind="cubic")
        self.Dc_to_z = interp1d(Dc, z, kind="cubic")
        self.z_to_luminosity_distance = interp1d(z, luminosity_distance, kind="cubic")

        # Lookup table for differential comoving distance
        differential_comoving_volume = (
            Planck18.differential_comoving_volume(z).value * 4 * np.pi
        )  # differential comoving volume in Mpc^3
        self.differential_comoving_volume = interp1d(
            z, differential_comoving_volume, kind="cubic"
        )

        return None

    def batch_handler(self, nsamples, sampling_routine, json_file, resume=False):
        """
        Function to handle the batch size.

        Parameters
        ----------
        nsamples : `int`
            number of samples.
        sampling_routine : `function`
            function to sample the parameters.
            e.g. unlensed_sampling_routine() or lensed_sampling_routine()
        json_file : `str`
            name of the json file to store the parameters.
        resume : `bool`
            if True, it will resume the sampling from the last batch.
            default resume = False.
        """
        batch_size = self.batch_size
        # if nsamples is multiple of batch_size
        if nsamples % batch_size == 0:
            num_batches = nsamples // batch_size
        # if nsamples is not multiple of batch_size
        else:
            num_batches = nsamples // batch_size + 1

        print(
            f"chosen batch size = {batch_size}. If you want to change batch size, self.batch_size = new_size"
        )
        print(f"There will be {num_batches} batche(s)")

        # note frac_batches+(num_batches-1)*batch_size = nsamples
        if nsamples > batch_size:
            frac_batches = nsamples - (num_batches - 1) * batch_size
        # if nsamples is less than batch_size
        else:
            frac_batches = nsamples
        track_batches = 0

        if not resume:
            track_batches = track_batches + 1
            print(f"Batch no. {track_batches}")
            # new first batch with the frac_batches
            sampling_routine(nsamples=frac_batches, file_name=json_file)
        else:
            # check where to resume from
            try:
                print(f"resuming from {json_file}")
                with open(json_file, "r", encoding="utf-8") as f:
                    data = json.load(f)
                    track_batches = (len(data["zs"]) - frac_batches) // batch_size + 1
            except:
                track_batches = track_batches + 1
                print(f"Batch no. {track_batches}")
                # new first batch with the frac_batches
                sampling_routine(nsamples=frac_batches, file_name=json_file)

        # ---------------------------------------------------#
        min_, max_ = track_batches, num_batches
        for i in range(min_, max_):
            track_batches = track_batches + 1
            print(f"Batch no. {track_batches}")
            sampling_routine(nsamples=batch_size, file_name=json_file, resume=True)
        # ---------------------------------------------------#

        return None

    def unlensed_sampling_routine(self, nsamples, file_name, resume=False):
        """
        Function to generate unlensed GW source parameters.

        Parameters
        ----------
        nsamples : `int`
            number of samples.
            default nsamples = 100000.
        file_name : `str`
            name of the json file to store the parameters.
        resume : `bool`
            if True, it will resume the sampling from the last batch.
            default resume = False.
        """
        # get gw params
        print("sampling gw source params...")
        gw_param = self.compact_binary_pop.sample_gw_parameters(nsamples=nsamples)
        # Get all of the signal to noise ratios
        print("calculating snrs...")
        snrs = self.snr.snr(GWparam_dict=gw_param)
        gw_param.update(snrs)

        # store all params in json file
        append_json(file_name=file_name, dictionary=gw_param, replace=not (resume))
        gw_param = None  # to free memory

    def unlensed_cbc_statistics(
        self, nsamples=None, resume=False, json_file="./gw_params.json", **kwargs
    ):
        """
        Function to generate unlensed GW source parameters.

        Parameters
        ----------
        nsamples : `int`
            number of samples.
            default nsamples = 100000.
        resume : `bool`
            resume = False (default) or True.
            if True, the function will resume from the last batch.
        json_file : `str`
            json file name for storing the parameters.
            default json_file = './gw_params.json'.
        kwargs : `dict`
            key word arguments for initializing the ``CompactBinaryPopulation`` class. \n
            This initialization is either done at the time of class initialization or at the time of calling this function. \n
            Following parameters can be provided, \n
            ``m_min`` : `float`
                minimum mass of the compact binary (single).
            ``m_max`` : `float`
                maximum mass of the compact binary (single).
            ``event_type`` : `str`
                event_type = 'popI_II' or `popIII` or `primordial`.
            ``src_model_params`` : `dict`
                src_model_params = {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82,\n
                'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08,\n
                'mu_g': 33.07, 'sigma_g': 5.69}}

        Returns
        ----------
        unlensed_gw_params : `dict`
            dictionary of unlensed GW source parameters.
            unlensed_gw_params.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time']

        """

        gw_sampler_dict = self.gw_param_sampler_dict

        # gw parameter sampling
        if nsamples:
            gw_sampler_dict["nsamples"] = nsamples
        else:
            nsamples = gw_sampler_dict["nsamples"]

        try:
            # check if kwargs is empty
            if kwargs:
                for key, value in kwargs.items():
                    if key in gw_sampler_dict:
                        gw_sampler_dict[key] = value
                # re-initializing classes with new params
                self.class_initialization()
        except:
            pass

        # sampling in batches
        self.batch_handler(
            nsamples=nsamples,
            sampling_routine=self.unlensed_sampling_routine,
            json_file=json_file,
            resume=resume,
        )

        gw_param = get_param_from_json(json_file)

        return gw_param

    def unlensed_rate(
        self,
        gw_param="./gw_params.json",
        snr_threshold=8.0,
        jsonfile="./gw_params_detectable.json",
        detectability_condition="step_function",
    ):
        """
        Function to calculate unlensed merger rate.

        .. math::
            R_U = \\mathcal{N}^U\\int dz_s R_o^U(z_s)\\bigg\\{\\Theta[\\rho(z_s,\\theta)-\\rho_{th}] P(\\theta) d\\theta \\bigg\\}

        - where :math:`\\mathcal{N}^U` is the normalization factor of the unlensed merger rate distribution wrt redshift.

        Parameters
        ----------
        gw_param : `dict` or `str` for json file name.
            dictionary of unlensed GW source parameters.
            default gw_param = './gw_params.json'.
        snr_threshold : `float`
            SNR threshold for detection.
            default snr_threshold = 8.
        jsonfile : `str`
            json file name for storing the detectable parameters.
            default jsonfile = './gw_params_detectable.json'.

        Returns
        ----------
        unlensed_rate : (`float`,`float`)
            unlensed merger rate in a year
            unlensed_rate[0] = total unlensed rate with step function
            unlensed_rate[1] = total unlensed rate with pdet function
        gw_param_detectable : `dict`
            dictionary of detectable unlensed GW source parameters.
            gw_param_detectable.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time']

        """
        # get gw params from json file if not provided
        if type(gw_param) == str:
            print(f"getting gw_params from json file {gw_param}...")
            gw_param = get_param_from_json(gw_param)

        # call json_file_ler_param and for adding the final results
        with open(self.json_file_ler_param, 'r', encoding='utf-8') as f:
            data = json.load(f)

        # if detectability_condition == "step_function" or snr is not None:
        if detectability_condition == "step_function":
            try:
                # get snr
                snr = gw_param["opt_snr_net"]
            except:
                # snr not provided
                print("snr not provided in gw_param dict. Exiting...")
                return None
            
            # selecting only detectable
            idx_detectable = snr > snr_threshold
            
            # montecarlo integration
            # The total rate is Eq. A4 of https://arxiv.org/pdf/2106.06303.pdf
            # R = C0 int Theta(rho-rhoc) p(z) p(theta) dtheta dz_s, where C0 = int R(zs)/(1+zs) dVc/dzs dzs is the normalization constant for p(z)
            # Thus R = C0 <Theta(rho-rhoc)>
            c0 = self.compact_binary_pop.normalization_pdf_z
            total_rate = c0 * np.mean(idx_detectable)
            print(f"total unlensed rate (yr^-1) (with step function): {total_rate}")
            # append the results
            data['unlensed_rate_step'] = total_rate

        else:
            # check if pdet is provided
            try:
                pdet = gw_param["pdet_net"]
            except:
                try:
                    snr = gw_param["opt_snr_net"]
                except:
                    print("pdet not provided in gw_param dict. Exiting...")
                    return None

            pdet = 1 - norm.cdf(snr_threshold - snr)
            gw_param["pdet_net"] = pdet

            # selecting only detectable
            idx_detectable = pdet > 0.5

            # with pdet
            # montecarlo integration
            # The total rate is Eq. A4 of https://arxiv.org/pdf/2106.06303.pdf
            # R = C0 int pdet p(z) p(theta) dtheta dz_s, where C0 = int R(zs)/(1+zs) dVc/dzs dzs is the normalization constant for p(z)
            # Thus R = C0 <pdet>
            c0 = self.compact_binary_pop.normalization_pdf_z
            total_rate = c0 * np.mean(pdet)
            print(f"total unlensed rate (yr^-1) (with pdet function): {total_rate}")
            # append the results
            data['unlensed_rate_pdet'] = total_rate

        # store all detectable params in json file
        for key, value in gw_param.items():
            gw_param[key] = value[idx_detectable]

        # store all detectable params in json file
        print(f"storing detectable unlensed params in {jsonfile}")
        append_json(jsonfile, gw_param, replace=True)

        # write the results
        append_json(self.json_file_ler_param, data, replace=True)

        # get the detectable params and return
        gw_param_detectable = get_param_from_json(jsonfile)

        return (total_rate, gw_param_detectable)

    def lensed_sampling_routine(self, nsamples, file_name, resume=False):
        """
        Function to generate lensed GW source parameters, lens galaxy parameters and image paramters.

        Parameters
        ----------
        nsamples : `int`
            number of samples.
        file_name : `str`
            name of the json file to store the parameters.
        resume : `bool`
            if True, it will resume the sampling from the last batch.
            default resume = False.
        """

        # get lensed params
        print("sampling lensed params...")
        lensed_param = self.lens_galaxy_pop.sample_lens_parameters(size=nsamples)
        # now get (strongly lensed) image paramters along with lens parameters
        lensed_param = self.lens_galaxy_pop.get_image_properties(
            n_min_images=self.lensed_param_sampler_dict["min_lensed_images"],
            n_max_images=self.lensed_param_sampler_dict["max_lensed_images"],
            lens_parameters=lensed_param,
            lensModelList=self.lensed_param_sampler_dict["lensModelList"],
            npool=self.npool,
        )
        # Get all of the signal to noise ratios
        print("calculating snrs...")
        snrs = self.lens_galaxy_pop.get_lensed_snrs(
            snr_calculator=self.snr,
            lensed_param=lensed_param,
            n_max_images=self.lensed_param_sampler_dict["max_lensed_images"],
        )
        lensed_param.update(snrs)

        # store all params in json file
        append_json(file_name=file_name, dictionary=lensed_param, replace=not (resume))
        lensed_param = None  # to free up memory

    def lensed_cbc_statistics(
        self, nsamples=None, resume=False, json_file="./lensed_params.json", **kwargs
    ):
        """
        Function to generate lensed GW source parameters, lens galaxy parameters and image paramters.

        Parameters
        ----------
        nsamples : `int`
            number of samples.
            default nsamples = 100000.
        resume : `bool`
            resume = False (default) or True.
            if True, the function will resume from the last batch.
        json_file : `str`
            json file name for storing the parameters.
            default json_file = './lensed_params.json'.
        kwargs : `dict`
            key word arguments for initializing the ``LensGalaxyPopulation`` class. \n
            This initialization is either done at the time of class initialization or at the time of calling this function. \n
            Following parameters can be provided, \n
            ``min_lensed_images`` : `int`
                minimum number of lensed images.
            ``max_lensed_images`` : `int`
                maximum number of lensed images.
            ``lensModelList`` : `list`
                list of lens models.
                e.g. lensModelList = ['EPL_NUMBA', 'SHEAR'].

        Returns
        ----------
        lensed_param : `dict`
            dictionary of lensed GW source parameters, lens galaxy parameters and image paramters.
            lensed_param.keys() = ['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2', 'Dl',
            'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source',
            'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'n_images',
            'x0_image_positions', 'x1_image_positions', 'magnifications', 'time_delays', 'traces',
            'determinants', 'image_type', 'weights', 'opt_snr_net', 'L1', 'H1', 'V1']
        """

        lens_sampler_dict = self.lensed_param_sampler_dict

        # if new paramteres are provided
        if nsamples:
            lens_sampler_dict["nsamples"] = nsamples
        else:
            nsamples = lens_sampler_dict["nsamples"]

        try:
            # check if kwargs is empty
            if kwargs:
                for key, value in kwargs.items():
                    if key in lens_sampler_dict:
                        lens_sampler_dict[key] = value

                # re-initializing classes with new params
                self.class_initialization()
        except:
            pass

        # gw_param will not be kept same as that of unlensed case. So, it is sampled newly
        # sampling in batches
        self.batch_handler(
            nsamples=nsamples,
            sampling_routine=self.lensed_sampling_routine,
            json_file=json_file,
            resume=resume,
        )

        lensed_param = get_param_from_json(json_file)

        return lensed_param

    def lensed_rate(
        self,
        lensed_param="./lensed_params.json",
        snr_threshold=8.0,
        num_img=2,
        jsonfile="./lensed_params_detectable.json",
        none_as_nan=True,
        detectability_condition="step_function",
    ):
        """
        Function to calculate lensed merger rate.

        .. math::
            R_L = \\mathcal{N}^L\\int dz_s R_o^L(z_s)\\bigg\\{\\Theta[\\rho(z_s,\\theta)-\\rho_{th}] P(\\theta) d\\theta \\bigg\\}

        - where :math:`\\mathcal{N}^L` is the normalization factor of the lensed merger rate distribution wrt redshift.

        Parameters
        ----------
        lensed_param : `dict` or `str` for the json file name.
            dictionary of lensed GW source parameters, lens galaxy parameters and image paramters.
            lensed_param.keys() = ['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2', 'Dl',
            'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source',
            'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'n_images',
            'x0_image_positions', 'x1_image_positions', 'magnifications', 'time_delays', 'traces',
            'determinants', 'image_type', 'weights', 'opt_snr_net', 'L1', 'H1', 'V1']
        snr_threshold : `float`
            threshold for detection signal to noise ratio.
            e.g. snr_threshold = 8.
        num_img : `int`
            number of images.
            e.g. num_img = 2.
        jsonfile : `str`
            json file name for storing the parameters.
            default jsonfile = './lensed_params_detectable.json'.
        none_as_nan : `bool`
            if True, replace None with np.nan in the lensed_param dictionary.
            default none_as_nan = True.

        Returns
        ----------
        lensed_rate : `float`
            lensed merger rate in a year.
            lensed_rate[0] = total lensed rate with step function
            lensed_rate[1] = total lensed rate with pdet function
        detectable_lensed_param : `dict`
            dictionary of detectable lensed GW source parameters, lens galaxy parameters and image paramters.
            detectable_lensed_param.keys() = ['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2',
            'Dl', 'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source',
            'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'n_images',
            'x0_image_positions', 'x1_image_positions', 'magnifications', 'time_delays', 'traces',
            'determinants', 'image_type', 'weights', 'opt_snr_net', 'L1', 'H1', 'V1']
        """

        # get lensed params from json file if not provided
        if type(lensed_param) == str:
            print(f"getting lensed_param from json file {lensed_param}...")
            lensed_param = get_param_from_json(lensed_param)

        # get size of the lensed_param for a parameter
        size = len(lensed_param["zs"])

        # call json_file_ler_param and for adding the final results
        with open(self.json_file_ler_param, 'r', encoding='utf-8') as f:
            data = json.load(f)

        # check for images with snr above threshold
        snr_threshold, num_img = np.array([snr_threshold]).reshape(-1), np.array(
            [num_img]).reshape(-1)  # convert to array

        # weights=1 if set minimum number of image is 2
        weights = lensed_param["weights"]
        # rejection sample wrt to weights
        not_rejected = np.random.uniform(0, 1, size) < weights

        # if detectability_condition == "step_function" or snr is not None:
        if detectability_condition == "step_function":
            try:
                # get snr, dimensions are (nsamples, n_max_images)
                snr = lensed_param["opt_snr_net"]
            except:
                # snr not provided
                print("snr not provided in lensed_param dict. Exiting...")
                return None

            snr_hit = np.full(size, True)  # boolean array to store the result of the threshold condition
            # for each row: choose a threshold and check if the number of images above threshold. Sum over the images. If sum is greater than num_img, then snr_hit = True 
            for i in range(len(snr_threshold)):
                snr_hit = snr_hit & (np.sum((snr > snr_threshold[i]), axis=1) >= num_img[i])
            snr_hit = snr_hit & not_rejected

            # montecarlo integration
            # The total rate is Eq. A4 of https://arxiv.org/pdf/2106.06303.pdf
            # R = C0 int Theta(rho-rhoc) p(z) p(theta) dtheta dz_s, where C0 = int R(zs)/(1+zs) dVc/dzs tau(zs) dzs is the normalization constant for p(z)
            # Thus R = C0 <Theta(rho-rhoc)>
            total_rate = self.lens_galaxy_pop.normalization_pdf_z * np.mean(snr_hit)
            print("total lensed rate (yr^-1) (with step function): {}".format(total_rate))
            # append the results
            data['lensed_rate_step'] = total_rate

        else:
            # get pdet, dimensions are (nsamples, n_max_images)
            # check if pdet is provided
            try:
                pdet = lensed_param["pdet_net"]
            except:
                try:
                    # dimensions are (nsamples, n_max_images)
                    pdet = []
                    snr = lensed_param["opt_snr_net"]
                    for i in range(len(snr_threshold)):
                        # store pdet for each threshold condition
                        pdet.append(1 - norm.cdf(snr_threshold[i] - snr))
                    pdet = np.array(pdet)
                except:
                    print("pdet not provided in lensed_param dict. Exiting...")
                    return None
                
            # product of pdet for all images
            # think of 3d layer. x-y correspons to snr, and z wrt snr_threshold
            pdet_combined = np.full(size, 1.0)
            for i in range(len(snr_threshold)):
                # sort pdet and use only num_img images
                pdet_combined = pdet_combined * np.prod(-np.sort(-pdet[i], axis=1)[:, :num_img[i]])
            lensed_param["pdet_net"] = pdet_combined

            # snr_hit
            snr_hit = pdet_combined > 0.5
            snr_hit = snr_hit & not_rejected

            # with pdet
            # montecarlo integration
            # The total rate is Eq. A4 of https://arxiv.org/pdf/2106.06303.pdf
            # R = C0 int Theta(rho-rhoc) p(z) p(theta) dtheta dz_s, where C0 = int R(zs)/(1+zs) dVc/dzs tau(zs) dzs is the normalization constant for p(z)
            # Thus R = C0 <pdet>
            c0 = self.compact_binary_pop.normalization_pdf_z
            total_rate = c0 * np.mean(pdet_combined * weights)
            print(f"total lensed rate (yr^-1) (with pdet function): {total_rate}")
            # append the results
            data['lensed_rate_pdet'] = total_rate

        # store all params in json file
        if none_as_nan == True:
            for key, value in lensed_param.items():
                lensed_param[key] = value[snr_hit]
        else:
            for key, value in lensed_param.items():
                lensed_param[key] = np.nan_to_num(value[snr_hit])

        # store all detectable params in json file
        print(f"storing detectable lensed params in {jsonfile}...")
        append_json(jsonfile, lensed_param, replace=True)

        # write the results
        append_json(self.json_file_ler_param, data, replace=True)
        # get the detectable params and return
        lensed_param_detectable = get_param_from_json(jsonfile)

        return (total_rate, lensed_param_detectable)
    
    def rate_comparision(self, detectability_condition="step_function"):
        """
        Function to calculate unlensed and lensed merger rate and their ratio. 
        It will get the unlensed_rate and lensed_rate from json_file_ler_param="./LeR_params.json"

        Parameters
        ----------
        detectability_condition : `str`
            detectability condition, either "step_function" or "pdet_function"

        Returns
        -------
        unlensed_rate : `float`
            unlensed merger rate
        lensed_rate : `float`
            lensed merger rate
        ratio : `float`
            ratio of lensed_rate and unlensed_rate

        """

        # call json_file_ler_param and add the results
        with open(self.json_file_ler_param, 'r', encoding='utf-8') as f:
            data = json.load(f)

        if detectability_condition == "step_function":
            try:
                unlensed_rate = data['unlensed_rate_step']
                lensed_rate = data['lensed_rate_step']
            except:
                print("unlensed_rate_step or lensed_rate_step not found in json file. Exiting...")
                return None
            rate_ratio = unlensed_rate / lensed_rate
            # append the results
            data['rate_ratio_step'] = rate_ratio

        elif detectability_condition == "pdet":
            try:
                unlensed_rate = data['unlensed_rate_pdet']
                lensed_rate = data['lensed_rate_pdet']
            except:
                print("unlensed_rate_pdet or lensed_rate_pdet not found in json file. Exiting...")
                return None
            rate_ratio = unlensed_rate / lensed_rate
            # append the results
            data['rate_ratio_pdet'] = rate_ratio
        else:
            print("detectability_condition should be either step_function or pdet_function")
            return None
        
        
        print(f"unlensed_rate: {unlensed_rate}")
        print(f"lensed_rate: {lensed_rate}")
        print(f"ratio: {rate_ratio}")

        # write the results
        append_json(self.json_file_ler_param, data, replace=True)

        return (unlensed_rate, lensed_rate, rate_ratio)

    def rate_comparision_with_rate_calculation(
        self,
        snr_threshold_unlensed=8.0,
        unlened_param="./gw_params.json",
        snr_threshold_lensed=8.0,
        num_img=2,
        lensed_param="./lensed_params.json",
        jsonfile_unlensed="./gw_params_detectable.json",
        jsonfile_lensed="./lensed_params_detectable.json",
        detectability_condition="step_function",
    ):
        """
        Function to calculate unlensed and lensed merger rate and their ratio.

        Parameters
        ----------
        snr_threshold_unlensed : `float`
            threshold for detection signal to noise ratio for unlensed case.
            e.g. snr_threshold_unlensed = 8.
        unlened_param : `dict`
            dictionary of unlensed GW source parameters.
            unlened_param.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time']
        snr_threshold_lensed : `float`
            threshold for detection signal to noise ratio for lensed case.
            e.g. snr_threshold_lensed = 8.
        num_img : `int`
            number of images crossing the threshold.
            e.g. num_img = 2.
        lensed_param : `dict`
            dictionary of lensed GW source parameters, lens galaxy parameters and image paramters.
            lensed_param.keys() = ['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2', 'Dl',
            'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source',
            'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'n_images',
            'x0_image_positions', 'x1_image_positions', 'magnifications', 'time_delays', 'traces',
            'determinants', 'image_type', 'weights', 'opt_snr_net', 'L1', 'H1', 'V1']
        jsonfile_unlensed : `str`
            json file name for storing the parameters for unlensed detectable case.
            default jsonfile_unlensed = './gw_params_detectable.json'.
        jsonfile_lensed : `str`
            json file name for storing the parameters for lensed detectable case.
            default jsonfile_lensed = './lensed_params_detectable.json'.

        Returns
        ----------
        unlensed_rate : (`float`,`float`)
            unlensed merger rate in a year
            unlensed_rate[0] = total unlensed rate with step function
            unlensed_rate[1] = total unlensed rate with pdet function
        lensed_rate : (`float`,`float`)
            lensed merger rate in a year
            lensed_rate[0] = total lensed rate with step function
            lensed_rate[1] = total lensed rate with pdet function
        rate_ratio : (`float`,`float`)
            unlensed/lensed rate ratio
            rate_ratio[0] = total unlensed/lensed rate ratio with step function
            rate_ratio[1] = total unlensed/lensed rate ratio with pdet function
        """

        # calculate unlensed rate
        # print(f'getting unlened_param from json file {unlened_param}...')
        unlensed_rate = self.unlensed_rate(
            gw_param=unlened_param,
            snr_threshold=snr_threshold_unlensed,
            jsonfile=jsonfile_unlensed,
            detectability_condition=detectability_condition,
        )[0]

        # calculate lensed rate
        # print(f'getting lensed_param from json file {lensed_param}...')
        lensed_rate = self.lensed_rate(
            lensed_param=lensed_param,
            snr_threshold=snr_threshold_lensed,
            num_img=num_img,
            jsonfile=jsonfile_lensed,
            detectability_condition=detectability_condition,
        )[0]

        # rate ratio
        rate_ratio = unlensed_rate / lensed_rate
        
        print("unlensed/lensed rate ratio = ", rate_ratio)
        # call json_file_ler_param and add the results
        with open(self.json_file_ler_param, 'r', encoding='utf-8') as f:
            data = json.load(f)
        # append the results
        if detectability_condition == "step_function":
            data['rate_ratio_step'] = rate_ratio
        elif detectability_condition == "pdet":
            data['rate_ratio_pdet'] = rate_ratio
        # write the results
        append_json(self.json_file_ler_param, data, replace=True)

        return (unlensed_rate, lensed_rate, rate_ratio)

    # ---------------------------------------------------#
    # functions for selecting n lensed detectable events #
    # ---------------------------------------------------#
    def selecting_n_lensed_detectable_events(
        self,
        nsamples=100,
        snr_threshold=8.0,
        num_img=2,
        resume=False,
        json_file="./lensed_params_detectable.json",
    ):
        """
        Function to select n lensed detectable events.

        Parameters
        ----------
        nsamples : `int`
            number of samples to be selected.
            default size = 100.
        snr_threshold : `float`
            threshold for detection signal to noise ratio.
            e.g. snr_threshold = 8. or [8.,6.]
        num_img : `int`
            number of images crossing the threshold.
            e.g. num_img = 2 or [1,1]
        resume : `bool`
            if True, it will resume the sampling from the last batch.
            default resume = False.
        json_file : `str`
            json file name for storing the parameters.
            default json_file = './lensed_params_detectable.json'.

        Returns
        ----------
        param_final : `dict`
            dictionary of lensed GW source parameters, lens galaxy parameters and image paramters.
            param_final.keys() = ['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2',
            'Dl', 'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source',
            'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'n_images',
            'x0_image_positions', 'x1_image_positions', 'magnifications', 'time_delays', 'image_type',
            'weights', 'opt_snr_net', 'L1', 'H1', 'V1']

        """
        try:
            self.batch_size > 1000
        except:
            self.batch_size = 1000

        if not resume:
            n = 0  # iterator
            try:
                os.remove(json_file)
            except:
                pass
        else:
            # get sample size as nsamples from json file
            param_final = get_param_from_json(json_file)
            n = len(param_final["zs"])
            del param_final

        buffer_file = "./lensed_params_buffer.json"
        print("collected number of events = ", n)
        while n < nsamples:
            # disable print statements
            with contextlib.redirect_stdout(None):
                self.lensed_sampling_routine(
                    nsamples=self.batch_size, file_name=buffer_file, resume=False
                )

                # Dimensions are (nsamples, n_max_images)
                lensed_param = get_param_from_json(buffer_file)

                # get snr
                snr = lensed_param["opt_snr_net"]
                size = len(snr)

                # dealing with snr_threshold and num_img
                snr_threshold, num_img = np.array([snr_threshold]).reshape(
                    -1
                ), np.array([num_img]).reshape(
                    -1
                )  # convert to array
                # sort in descending order of each row
                arg_th = (-snr_threshold).argsort()
                sorted_snr = -np.sort(-snr, axis=1)
                num1 = 0  # tracks the number of images for the current threshold
                num2 = 0  # tracks the column number of the already sorted snr 2D array
                # boolean array to store the result of the threshold condition
                snr_hit = np.full(len(snr), True)
                for i in arg_th:
                    num1 = num_img[i]
                    for j in range(num1):
                        # snr_hit step function case
                        snr_hit = snr_hit & (sorted_snr[:, num2] > snr_threshold[i])
                        num2 += 1

                weights = lensed_param["weights"]
                # rejection sample wrt to weights
                not_rejected = np.random.uniform(0, 1, size) < weights
                snr_hit = snr_hit & not_rejected

                # store all params in json file
                for key, value in lensed_param.items():
                    lensed_param[key] = np.nan_to_num(value[snr_hit])
                append_json(json_file, lensed_param, replace=False)

                n += np.sum(snr_hit)
            print("collected number of events = ", n)

        # trim the final param dictionary
        print(f"trmming final result to size={nsamples}")
        param_final = get_param_from_json(json_file)
        # trim the final param dictionary
        idx = np.random.choice(len(param_final["zs"]), nsamples, replace=False)
        for key, value in param_final.items():
            param_final[key] = param_final[key][idx]

        # save the final param dictionary
        append_json(json_file, param_final, replace=True)

        return param_final

    def relative_mu_dt_lensed(self, lensed_param, snr_threshold=[8.0, 8.0]):
        """
        Function to classify the lensed images wrt to the morse phase difference.

        Parameters
        ----------
        lensed_param : `dict`
            dictionary of lensed GW source parameters, lens galaxy parameters and image paramters.
            lensed_param.keys() = ['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2', 'Dl',
            'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source',
            'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'n_images',
            'x0_image_positions', 'x1_image_positions', 'magnifications', 'time_delays', 'traces',
            'determinants', 'image_type', 'weights', 'opt_snr_net', 'L1', 'H1', 'V1']
        snr_threshold : `float`
            threshold for detection signal to noise ratio.
            e.g. snr_threshold = [8.,8.] or [8.,6.] for subthreshold

        Returns
        ----------
        mu_rel0 : `float.array`
            relative magnification for 0 degree phase difference.
        dt_rel0 : `float.array`
            relative time delay for 0 degree phase difference.
        mu_rel90 : `float.array`
            relative magnification for 90 degree phase difference.
        dt_rel90 : `float.array`
            relative time delay for 90 degree phase difference.
        """

        # get magnifications, time_delays and snr
        mu = np.nan_to_num(lensed_param["magnifications"])
        dt = np.nan_to_num(lensed_param["time_delays"])
        snr = np.nan_to_num(lensed_param["opt_snr_net"])

        # for 0 degree phase difference
        # get the index of the image which cross the threshold
        # get snr_threshold sorted first in descending order
        snr_threshold = -np.sort(-np.array(snr_threshold))
        # for type I
        snr1 = -np.sort(-snr[:, [0, 1]], axis=1)
        # for type II
        snr2 = -np.sort(-snr[:, [2, 3]], axis=1)

        # checking for zero values
        # check for threshold condition
        idx1, idx2 = [], []
        for i in range(len(snr)):
            if (
                any(x != 0.0 for x in snr1[i])
                and snr1[i][0] > snr_threshold[0]
                and snr1[i][1] > snr_threshold[1]
            ):
                idx1.append(i)
            if (
                any(x != 0.0 for x in snr2[i])
                and snr2[i][0] > snr_threshold[0]
                and snr2[i][1] > snr_threshold[1]
            ):
                idx2.append(i)

        # combine magnifications and time_delays
        mu_ = np.concatenate((mu[idx1][:, [0, 1]], mu[idx2][:, [2, 3]]), axis=0)
        dt_ = np.concatenate((dt[idx1][:, [0, 1]], dt[idx2][:, [2, 3]]), axis=0) / (
            60 * 60 * 24
        )  # to days

        # relative magnification
        mu_rel0 = np.abs(mu_[:, 1] / mu_[:, 0])
        # relative time delay
        dt_rel0 = np.abs(dt_[:, 1] - dt_[:, 0])

        # for 90 degree phase difference
        # for type I
        snr1 = -np.sort(-snr[:, [0, 2]], axis=1)
        # for type II
        snr2 = -np.sort(-snr[:, [1, 3]], axis=1)

        # checking for zero values
        # check for threshold condition
        idx1, idx2 = [], []
        for i in range(len(snr)):
            if (
                any(x != 0.0 for x in snr1[i])
                and snr1[i][0] > snr_threshold[0]
                and snr1[i][1] > snr_threshold[1]
            ):
                idx1.append(i)
            if (
                any(x != 0.0 for x in snr2[i])
                and snr2[i][0] > snr_threshold[0]
                and snr2[i][1] > snr_threshold[1]
            ):
                idx2.append(i)

        # combine magnifications and time_delays
        mu_ = np.concatenate((mu[idx1][:, [0, 2]], mu[idx2][:, [1, 3]]), axis=0)
        dt_ = np.concatenate((dt[idx1][:, [0, 2]], dt[idx2][:, [1, 3]]), axis=0) / (
            60 * 60 * 24
        )  # in days

        # relative magnification
        mu_rel90 = np.abs(mu_[:, 1] / mu_[:, 0])
        # relative time delay
        dt_rel90 = np.abs(dt_[:, 1] - dt_[:, 0])

        return (mu_rel0, dt_rel0, mu_rel90, dt_rel90)

    def mu_vs_dt_plot(
        self,
        x_array,
        y_array,
        savefig=False,
        ax=None,
        colors="blue",
        linestyles="-",
        origin="upper",
        alpha=0.6,
        extent=[1e-2, 5e2, 1e-2, 1e2],
        contour_levels=[0.10, 0.40, 0.68, 0.95],
    ):
        """
        Function to generate 2D KDE and plot the relative magnification vs time delay difference for lensed samples.

        Parameters
        ----------
        x_array : `float.array`
            x array.
        y_array : `float.array`
            y array.
        xlabel : `str`
            x label.
        ylabel : `str`
            y label.
        title : `str`
            title.
        savefig : `bool`
            if True, it will save the figure.
            default savefig = False.
        ax : `matplotlib.axes`
            matplotlib axes.
            default ax = None.
        colors : `str`
            color of the plot.
            default colors = 'blue'.
        linestyles : `str`
            linestyle of the plot.
            default linestyles = '-'.
        origin : `str`
            origin of the plot.
            default origin = 'upper'.
        alpha : `float`
            alpha of the plot.
            default alpha = 0.6.
        extent : `list`
            extent of the plot.
            default extent = [1e-2,5e2,1e-2,1e2].
        contour_levels : `list`
            contour levels of the plot.
            default contour_levels = [0.10,0.40,0.68,0.95] which corresponds to 1,2,3,4 sigma.

        Returns
        ----------
        None

        """
        # applying cutt-off
        idx = (
            (x_array > extent[0])
            & (x_array < extent[1])
            & (y_array > extent[2])
            & (y_array < extent[3])
        )
        x_array = x_array[idx]
        y_array = y_array[idx]

        xu = np.log10(x_array)
        yu = np.log10(y_array)

        xmin = np.log10(1e-2)
        xmax = np.log10(5e2)
        ymin = np.log10(1e-2)
        ymax = np.log10(1e2)

        xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
        positions = np.vstack([xx.ravel(), yy.ravel()])
        values = np.vstack([xu, yu])
        kernel = gaussian_kde(values)
        ff = np.reshape(kernel(positions).T, xx.shape)

        zsort = -np.sort(-ff.flatten())

        cumz = np.cumsum(zsort) / np.sum(zsort)
        spl = interp1d(cumz, zsort)

        levels = []
        for i in contour_levels:
            levels.append(spl(i))
        levels = np.array(levels)[::-1]

        ax.contour(
            np.rot90(ff),
            levels,
            colors=colors,
            linestyles=linestyles,
            origin=origin,
            alpha=alpha,
            extent=np.log10(extent),
        )

        # labels
        ax.xlabel(r"$log_{10}\Delta t$ (days)")
        ax.ylabel(r"$\Delta log_{10}\mu$")
        ax.title(r"relative magnification vs relative time delay")

        # save figure
        if savefig:
            ax.savefig("mu_vs_dt.png", dpi=300, bbox_inches="tight")

        return None

    def selecting_n_unlensed_detectable_events(
        self,
        nsamples=100,
        snr_threshold=8.0,
        resume=False,
        json_file="./gw_params_detectable.json",
    ):
        """
        Function to select n unlensed detectable events.

        Parameters
        ----------
        nsamples : `int`
            number of samples to be selected.
            default size = 100.
        snr_threshold : `float`
            threshold for detection signal to noise ratio.
            e.g. snr_threshold = 8.
        resume : `bool`
            if True, it will resume the sampling from the last batch.
            default resume = False.
        json_file : `str`
            json file name for storing the parameters.
            default json_file = './gw_params_detectable.json'.

        Returns
        ----------
        param_final : `dict`
            dictionary of unlensed GW source parameters.
            param_final.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time']

        """

        try:
            self.batch_size > 1000
        except:
            self.batch_size = 1000

        if not resume:
            n = 0  # iterator
            try:
                os.remove(json_file)
            except:
                pass
        else:
            # get sample size as nsamples from json file
            param_final = get_param_from_json(json_file)
            n = len(param_final["zs"])
            del param_final

        buffer_file = "./gw_params_buffer.json"
        print("collected number of events = ", n)
        while n < nsamples:
            # disable print statements
            with contextlib.redirect_stdout(None):
                self.unlensed_sampling_routine(
                    nsamples=self.batch_size, file_name=buffer_file, resume=False
                )

                # get unlensed params
                unlensed_param = get_param_from_json(buffer_file)

                # get snr
                snr = unlensed_param["opt_snr_net"]
                # index of detectable events
                idx = snr > snr_threshold

                # store all params in json file
                for key, value in unlensed_param.items():
                    unlensed_param[key] = value[idx]
                append_json(json_file, unlensed_param, replace=False)

                n += np.sum(idx)
            print("collected number of events = ", n)

        # trim the final param dictionary
        print(f"trmming final result to size={nsamples}")
        param_final = get_param_from_json(json_file)
        # trim the final param dictionary
        idx = np.random.choice(len(param_final["zs"]), nsamples, replace=False)
        for key, value in param_final.items():
            param_final[key] = param_final[key][idx]

        # save the final param dictionary
        append_json(json_file, param_final, replace=True)

        return param_final

    def relative_mu_dt_unlensed(self, param, size=100):
        """
        Function to generate relative magnification vs time delay difference for unlensed samples.

        Parameters
        ----------
        param : `dict`
            dictionary of unlensed GW source parameters.
            unlensed_param.keys() = ['m1', 'm2', 'z', 'snr', 'theta_jn', 'ra', 'dec', 'psi', 'phase', 'geocent_time']

        Returns
        ----------
        dmu : `float.array`
            relative magnification.
        dt : `float.array`
            relative time delay.

        """

        t = param["geocent_time"]
        mu = param["luminosity_distance"]

        len_ = len(t)
        t_ = []
        mu_ = []
        while len(t_) < size:
            idx1, idx2 = random.sample(range(len_), 2)
            t_.append(t[idx2] - t[idx1])
            mu_.append(mu[idx2] / mu[idx1])

        dt = np.abs(np.array(t_)) / (60 * 60 * 24)  # in days
        dmu = np.sqrt(np.abs(np.array(mu_)))

        return (dmu, dt)
