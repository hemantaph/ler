# -*- coding: utf-8 -*-
"""
Functions and objects related to the main ler .
"""
import os
import json
import contextlib
import numpy as np
from gwsnr import GWSNR
from scipy.stats import norm
from scipy.interpolate import interp1d
from scipy.integrate import quad
from astropy.cosmology import Planck18
from ler.lens_galaxy_population import LensGalaxyPopulation
from ler.source_population import CompactBinaryPopulation
from ler.helperroutines import NumpyEncoder
from ler.helperroutines import append_json, dict_list_to_ndarray, trim_dictionary

# Conversions from SI units to CGS units
C = 299792458.  # m/s
G = 6.67408*1e-11  # m^3/kg/s^2 


class LeR():
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
        +---------------------------------------+--------------------------------------------------+
        | Atrributes                            | Type                                             |
        +=======================================+==================================================+
        | :attr:`~gw_param`                     | ``dict``                                         |
        +---------------------------------------+--------------------------------------------------+
        | :attr:`~gw_param_detectable`          | ``dict``                                         |
        +---------------------------------------+--------------------------------------------------+
        | :attr:`~lensed_param`                 | ``dict``                                         |
        +---------------------------------------+--------------------------------------------------+
        | :attr:`~lensed_param_detectable`      | ``dict``                                         |
        +---------------------------------------+--------------------------------------------------+
        | :attr:`~gw_param_sampler_dict`        | ``dict``                                         |
        +---------------------------------------+--------------------------------------------------+
        | :attr:`~lensed_param_sampler_dict`    | ``dict``                                         |
        +---------------------------------------+--------------------------------------------------+
        | :attr:`~snr_calculator_dict`          | ``dict``                                         |
        +---------------------------------------+--------------------------------------------------+
        | :attr:`~z_to_Dc`                      | ``scipy.interpolate.interp1d``                   |
        +---------------------------------------+--------------------------------------------------+
        | :attr:`~Dc_to_z`                      | ``scipy.interpolate.interp1d``                   |
        +---------------------------------------+--------------------------------------------------+
        | :attr:`~z_to_luminosity_distance`     | ``scipy.interpolate.interp1d``                   |
        +---------------------------------------+--------------------------------------------------+
        | :attr:`~differential_comoving_volume` | ``scipy.interpolate.interp1d``                   |
        +---------------------------------------+--------------------------------------------------+
        | :attr:`~compact_binary_pop`           | ``CompactBinaryPopulation class``                |
        +---------------------------------------+--------------------------------------------------+
        | :attr:`~lens_galaxy_pop`              | ``LensGalaxyPopulation class``                   |
        +---------------------------------------+--------------------------------------------------+
        | :attr:`~snr`                          | ``gwsnr package``                                |
        +---------------------------------------+--------------------------------------------------+

    Instance Methods
    ----------
    LeR class has the following method(s), \n
        +---------------------------------------+--------------------------------------------------+
        | Method(s)                             | Description                                      |
        +=======================================+==================================================+
        | :meth:`~gwsnr_intialization`          | Function for initializing the ``gwsnr`` package. |
        +---------------------------------------+--------------------------------------------------+
        | :meth:`~create_lookup_tables`         | To creating lookup tables for fast calculation   |
        |                                       | for the following conversion operations,         |
        |                                       | redshift to co-moving distance.                  |
        |                                       | co-moving distance to redshift.                  |
        |                                       | redshift to luminosity distance.                 |
        +---------------------------------------+--------------------------------------------------+
        | :meth:`~unlensed_cbc_statistics`      | Function to generate unlensed GW source          |
        |                                       | parameters.                                      |
        +---------------------------------------+--------------------------------------------------+
        | :meth:`~unlensed_rate`                | Function to calculate unlensed merger rate.      |
        +---------------------------------------+--------------------------------------------------+
        | :meth:`~lensed_cbc_statistics`        | Function to generate lensed GW source            |
        |                                       | parameters.                                      |
        +---------------------------------------+--------------------------------------------------+
        | :meth:`~lensed_rate`                  | Function to calculate lensed merger rate.        |
        +---------------------------------------+--------------------------------------------------+
    
    """

    gw_param_sampler_dict = None
    """``dict`` \n
    dictionary of params for initializing ``CompactBinaryPopulation`` class \n
    this will be used for GW unlensed parameters sampling \n
    gw_param_sampler_dict.keys() = ['nsamples', 'm_min', 'm_max', 'z_min', 'z_max', 'event_type', 'model_pars']
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

    def __init__(self, nsamples=100000, npool=int(4), z_min=0., z_max=10., batch_size=25000,
                 snr_finder='gwsnr', **kwargs):

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

        # dictionary of params for sampler
        # for unlened case (source params)
        # defualt for 'model_pars' is set for popI_II PowerLaw+PEAK model
        # for other models, please change the 'model_pars' accordingly
        self.gw_param_sampler_dict = {'nsamples': nsamples, 'm_min': 4.59, 'm_max': 86.22, 'z_min': z_min, 'z_max': z_max,
                                      'event_type': 'popI_II', 'model_pars': {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82,
                                                                                      'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08,
                                                                                      'mu_g': 33.07, 'sigma_g': 5.69}}
        # for lensed case
        # set 'min_lensed_images' = 2 for double image lensed case
        self.lensed_param_sampler_dict = {'nsamples': nsamples, 'min_lensed_images': 2, 'max_lensed_images': 4,
                                          'lensModelList': ['EPL_NUMBA', 'SHEAR']}

        # for snr_calculator
        # for 'waveform_approximant' other than IMRPhenomD or TaylorF2, please set 'snr_type' = 'inner_product'
        # you will get accurate results if you set 'nsamples_mtot': 200, 'nsamples_mass_ratio': 500., 'sampling_frequency': 4096.
        self.snr_calculator_dict = {'mtot_min':2., 'mtot_max':439.6, 'nsamples_mtot': 100, 'nsamples_mass_ratio': 50,
                                    'sampling_frequency': 2048.,
                                    'waveform_approximant': "IMRPhenomD", 'minimum_frequency': 20.,
                                    'snr_type': 'interpolation', 'waveform_inspiral_must_be_above_fmin': False,
                                    'psds': False, 'psd_file': False, 'ifos': False}

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
        if snr_finder == 'gwsnr':
            # default
            self.snr = self.gwsnr_intialization(kwargs)
        else:
            # custom SNR function
            self.snr = None
        # extra note on how to change snr finder function
        # self.snr = custom_snr_finder_function()

        # Create lookup tables
        self.create_lookup_tables(z_min, z_max)

        self.store_ler_params()

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
        if self._gw_param=='default':
            f = open ('gw_params.json', "r", encoding='utf-8')
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
        if self._gw_param_detectable=='default':
            f = open ('gw_params_detectable.json', "r", encoding='utf-8')
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
        if self._lensed_param=='default':
            f = open ('lensed_params.json', "r", encoding='utf-8')
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
        if self._lensed_param_detectable=='default':
            f = open ('lensed_params_detectable.json', "r", encoding='utf-8')
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
        self.compact_binary_pop = CompactBinaryPopulation(z_min=self.z_min, z_max=self.z_max,
                                                          m_min=self.gw_param_sampler_dict['m_min'],
                                                          m_max=self.gw_param_sampler_dict['m_max'],
                                                          event_type=self.gw_param_sampler_dict['event_type'],
                                                          model_pars=self.gw_param_sampler_dict['model_pars'])
        self.lens_galaxy_pop = LensGalaxyPopulation(self.compact_binary_pop)

        return None

    def store_ler_params(self):
        """
        Fuction to store the parameters of the LER model. This is useful for reproducing the results.
        """
        # store gw_param_sampler_dict, lensed_param_sampler_dict and snr_calculator_dict
        parameters_dict = {}
        parameters_dict.update(
            {'gw_param_sampler_dict': self.gw_param_sampler_dict})
        parameters_dict.update(
            {'lensed_param_sampler_dict': self.lensed_param_sampler_dict})
        
        snr_calculator_dict = self.snr_calculator_dict.copy()
        del snr_calculator_dict['ifos']
        parameters_dict.update(
            {'snr_calculator_dict': snr_calculator_dict})
        
        file_name = './LeR_params.json'
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

        snr_ = GWSNR(npool=self.npool, 
                     mtot_min=gwsnr_param_dict['mtot_min'],
                     mtot_max=gwsnr_param_dict['mtot_max'],
                     nsamples_mtot=gwsnr_param_dict['nsamples_mtot'],
                     nsamples_mass_ratio=gwsnr_param_dict['nsamples_mass_ratio'],
                     sampling_frequency=gwsnr_param_dict['sampling_frequency'],
                     waveform_approximant=gwsnr_param_dict['waveform_approximant'],
                     minimum_frequency=gwsnr_param_dict['minimum_frequency'],
                     snr_type=gwsnr_param_dict['snr_type'],
                     waveform_inspiral_must_be_above_fmin=gwsnr_param_dict[
                         'waveform_inspiral_must_be_above_fmin'],
                     psds=gwsnr_param_dict['psds'],
                     psd_file=gwsnr_param_dict['psd_file'],
                     ifos=gwsnr_param_dict['ifos'])

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
        self.z_to_Dc = interp1d(z, Dc, kind='cubic')
        self.Dc_to_z = interp1d(Dc, z, kind='cubic')
        self.z_to_luminosity_distance = interp1d(z, luminosity_distance, kind='cubic')

        # Lookup table for differential comoving distance
        differential_comoving_volume = Planck18.differential_comoving_volume(
            z).value*4*np.pi  # differential comoving volume in Mpc^3
        self.differential_comoving_volume = interp1d(
            z, differential_comoving_volume, kind='cubic')
        
        return None
    
    def batch_handler(self,nsamples,sampling_routine,json_file,resume=False):
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
        if nsamples%batch_size == 0:  
            num_batches = nsamples//batch_size
        # if nsamples is not multiple of batch_size
        else:
            num_batches = nsamples//batch_size +1

        print(f"chosen batch size = {batch_size}. If you want to change batch size, self.batch_size = new_size")
        print(f"There will be {num_batches} batche(s)")

        # note frac_batches+(num_batches-1)*batch_size = nsamples
        if nsamples > batch_size:
            frac_batches = nsamples-(num_batches-1)*batch_size
        # if nsamples is less than batch_size
        else:
            frac_batches = nsamples
        track_batches = 0

        if not resume:
            track_batches = track_batches+1
            print(f"Batch no. {track_batches}")
            # new first batch with the frac_batches
            sampling_routine(nsamples=frac_batches,file_name=json_file)
        else:
            # check where to resume from
            try:
                print(f"resuming from {json_file}")
                with open(json_file, "r", encoding='utf-8') as f:
                    data = json.load(f)
                    track_batches = (len(data['zs'])-frac_batches)//batch_size + 1
            except:
                print(f"no json file found with name {json_file}. Set resume=False")
                return None
            
        # ---------------------------------------------------#
        min_, max_ = track_batches, num_batches
        for i in range(min_,max_):
            track_batches = track_batches+1
            print(f"Batch no. {track_batches}")
            sampling_routine(nsamples=batch_size,file_name=json_file,resume=True)     
        # ---------------------------------------------------#

        return None

    def unlensed_sampling_routine(self, nsamples,file_name,resume=False):
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
        print('sampling gw source params...')
        gw_param = self.compact_binary_pop.sample_gw_parameters(
            nsamples=nsamples)
        # Get all of the signal to noise ratios
        print('calculating snrs...')
        snrs = self.snr.snr(GWparam_dict=gw_param)
        gw_param.update(snrs)

        # store all params in json file
        append_json(file_name=file_name, dictionary=gw_param, replace=not(resume))
        gw_param = None # to free memory

    def unlensed_cbc_statistics(self, nsamples=None, resume=False, json_file='./gw_params.json', **kwargs):
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
            ``model_pars`` : `dict`
                model_pars = {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82,\n
                'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08,\n
                'mu_g': 33.07, 'sigma_g': 5.69}}
        
        Returns
        ----------
        unlensed_gw_params : `dict`
            dictionary of unlensed GW source parameters.
            unlensed_gw_params.keys() = ['m1', 'm2', 'z', 'Dc', 'Dl', 'Dl_obs', 'snr', 'pdet', 'event_type']

        """
        
        gw_sampler_dict = self.gw_param_sampler_dict

        # gw parameter sampling
        if nsamples:
            gw_sampler_dict['nsamples'] = nsamples
        else:
            nsamples = gw_sampler_dict['nsamples']

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
        self.batch_handler(nsamples=nsamples,sampling_routine=self.unlensed_sampling_routine,
                           json_file=json_file,resume=resume)
        
        self.gw_param = 'default'
        gw_param = self.gw_param.copy()
        self.gw_param = None  # to free up memory 

        return gw_param
                           
    def unlensed_rate(self, gw_param='default', snr_threshold=8., jsonfile='./gw_params_detectable.json'):
        """
        Function to calculate unlensed merger rate. 

        .. math::
            R_U = \\mathcal{N}^U\\int dz_s R_o^U(z_s)\\bigg\\{\\Theta[\\rho(z_s,\\theta)-\\rho_{th}] P(\\theta) d\\theta \\bigg\\}

        - where :math:`\\mathcal{N}^U` is the normalization factor of the unlensed merger rate distribution wrt redshift.

        Parameters
        ----------
        size : `int`
            number of samples.
        snr_threshold : `float`
            threshold for detection signal to noise ratio.
            e.g. snr_threshold = 8.
        jsonfile : `bool`
            if True, store all gravitational waves source parameters in json file \n
            (for detected sources).

        Returns
        ----------
        unlensed_rate : `float`
            unlensed merger rate in a year.
        """
        # get gw params from json file if not provided
        if gw_param == 'default':
            print('getting gw_params from json file...')
            self.gw_param = 'default'
            gw_param = self.gw_param.copy()
            self.gw_param = None  # to free up memory

        gw_param = dict_list_to_ndarray(gw_param)
        # get snr
        snr = gw_param['opt_snr_net']
        pdet = 1 - norm.cdf(snr_threshold - snr)
        gw_param['pdet_net'] = pdet

        # selecting only detectable
        idx_detectable = snr > snr_threshold
        # store all detectable params in json file
        for key, value in gw_param.items():
            gw_param[key] = value[idx_detectable]

        # store all detectable params in json file
        append_json(jsonfile, gw_param, replace=True)

        # montecarlo integration
        # The total rate is Eq. A4 of https://arxiv.org/pdf/2106.06303.pdf
        # R = C0 int Theta(rho-rhoc) p(z) p(theta) dtheta dz_s, where C0 = int R(zs)/(1+zs) dVc/dzs dzs is the normalization constant for p(z)
        # Thus R = C0 <Theta(rho-rhoc)>
        c0 = self.compact_binary_pop.normalization_pdf_z
        total_rate_step = c0 * np.mean(idx_detectable)
        print(f"total unlensed rate with step function: {total_rate_step}")
        # with pdet
        total_rate_pdet = c0 * np.mean(pdet)
        print(f"total unlensed rate with pdet function: {total_rate_pdet}")

        self.gw_param_detectable = 'default'
        gw_param_detectable = self.gw_param_detectable.copy()
        self.gw_param_detectable = None  # to free up memory

        return([total_rate_step,total_rate_pdet], gw_param_detectable)

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
        print('sampling lensed params...')
        lensed_param = self.lens_galaxy_pop.sample_lens_parameters(
            size=nsamples)
        # now get (strongly lensed) image paramters along with lens parameters
        lensed_param = \
            self.lens_galaxy_pop.get_image_properties(n_min_images=self.lensed_param_sampler_dict['min_lensed_images'], 
                                                      n_max_images=self.lensed_param_sampler_dict['max_lensed_images'], 
                                                      lens_parameters=lensed_param, 
                                                      lensModelList=self.lensed_param_sampler_dict['lensModelList'], 
                                                      npool=self.npool)
        # Get all of the signal to noise ratios
        print('calculating snrs...')
        snrs = self.lens_galaxy_pop.get_lensed_snrs(
            snr_calculator=self.snr, lensed_param=lensed_param, 
            n_max_images=self.lensed_param_sampler_dict['max_lensed_images'])
        lensed_param.update(snrs)

        # store all params in json file
        append_json(file_name=file_name, dictionary=lensed_param, replace=not(resume))
        lensed_param = None  # to free up memory


    def lensed_cbc_statistics(self, nsamples=None, resume=False, json_file='./lensed_params.json', **kwargs):
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
            lens_sampler_dict['nsamples'] = nsamples
        else:
            nsamples = lens_sampler_dict['nsamples']

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
        self.batch_handler(nsamples=nsamples,sampling_routine=self.lensed_sampling_routine,
                            json_file=json_file,resume=resume)
        
        self.lensed_param = 'default'
        lensed_param = self.lensed_param.copy()
        self.lensed_param = None  # to free up memory

        return lensed_param

    def lensed_rate(self, lensed_param='default', snr_threshold=8., num_img=2, 
                    jsonfile='./lensed_params_detectable.json', none_as_nan=True):
        """
        Function to calculate lensed merger rate.

        .. math::
            R_L = \\mathcal{N}^L\\int dz_s R_o^L(z_s)\\bigg\\{\\Theta[\\rho(z_s,\\theta)-\\rho_{th}] P(\\theta) d\\theta \\bigg\\}

        - where :math:`\\mathcal{N}^L` is the normalization factor of the lensed merger rate distribution wrt redshift.

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
            with step function and pdet function.
        detectable_lensed_param : `dict`
            dictionary of detectable lensed GW source parameters, lens galaxy parameters and image paramters.
            detectable_lensed_param.keys() = ['zl', 'zs', 'sigma', 'q', 'e1', 'e2', 'gamma1', 'gamma2', 
            'Dl', 'Ds', 'Dls', 'theta_E', 'gamma', 'mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 
            'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'n_images', 
            'x0_image_positions', 'x1_image_positions', 'magnifications', 'time_delays', 'traces', 
            'determinants', 'image_type', 'weights', 'opt_snr_net', 'L1', 'H1', 'V1']
        """

        # get lensed params from json file if not provided
        if lensed_param == 'default':
            print('getting lensed_params from json file...')
            self.lensed_param = 'default'
            lensed_param = self.lensed_param.copy()
            self.lensed_param = None

        lensed_param = dict_list_to_ndarray(lensed_param)

        # Dimensions are (nsamples, n_max_images)
        snr = np.array(lensed_param['opt_snr_net'])
        size = len(snr)
        snr_threshold, num_img = np.array(
            [snr_threshold]).reshape(-1), np.array([num_img]).reshape(-1)  # convert to array
        # sort in descending order of each row
        arg_th = (-snr_threshold).argsort()
        sorted_snr = -np.sort(-snr, axis=1)
        num1 = 0  # tracks the number of images for the current threshold
        num2 = 0  # tracks the column number of the already sorted snr 2D array
        # boolean array to store the result of the threshold condition
        snr_hit = np.full(len(snr), True)
        # array to store the result of the pdet condition
        pdet_combined = np.full(len(snr), 1.)
        for i in arg_th:
            num1 = num_img[i]
            for j in range(num1):
                # snr_hit step function case
                snr_hit = snr_hit & (sorted_snr[:, num2] > snr_threshold[i])
                # pdet for probability of detection
                pdet = 1 - \
                    norm.cdf(snr_threshold[i] -
                             np.nan_to_num(sorted_snr[:, num2]))
                pdet_combined = pdet_combined*pdet
                num2 += 1
        lensed_param['pdet_net'] = pdet

        weights = lensed_param['weights']
        # rejection sample wrt to weights
        not_rejected = np.random.uniform(0,1,size)<weights
        snr_hit = snr_hit&not_rejected
        
        # store all params in json file
        if none_as_nan == True:
            for key, value in lensed_param.items():
                lensed_param[key] = value[snr_hit]
        else:
            for key, value in lensed_param.items():
                lensed_param[key] = np.nan_to_num(np.array(value[snr_hit]))

        # store all detectable params in json file
        append_json(jsonfile, lensed_param, replace=True)

        # montecarlo integration
        # The total rate is Eq. A4 of https://arxiv.org/pdf/2106.06303.pdf
        # R = C0 int Theta(rho-rhoc) p(z) p(theta) dtheta dz_s, where C0 = int R(zs)/(1+zs) dVc/dzs tau(zs) dzs is the normalization constant for p(z)
        # Thus R = C0 <Theta(rho-rhoc)>
        c0 = self.lens_galaxy_pop.normalization_pdf_z
        total_rate_step = c0 * np.mean(snr_hit)
        print("total lensed rate with step function: {}".format(total_rate_step))

        # with pdet 
        total_rate_pdet = c0 * np.mean(pdet_combined*weights)
        print("total lensed rate with pdet function: {}".format(total_rate_pdet))

        self.lensed_param_detectable = 'default'
        lensed_param_detectable = self.lensed_param_detectable.copy()
        self.lensed_param_detectable = None  # to free up memory

        return([total_rate_step,total_rate_pdet], lensed_param_detectable)
    

    # ---------------------------------------------------#
    # functions for selecting n lensed detectable events #
    # ---------------------------------------------------#

    def selecting_n_lensed_detectable_events(self, size=100, snr_threshold=8., num_img=2, 
                                             none_as_nan=False, lenstype='any',
                                             resume=False, json_file='./lensed_params_detectable.json'):
        """
        Function to select n lensed detectable events.

        Parameters
        ----------
        size : `int`
            number of samples to be selected.
            default size = 100.
        snr_threshold : `float`
            threshold for detection signal to noise ratio.
            e.g. snr_threshold = 8. or [8.,6.]
        num_img : `int`
            number of images crossing the threshold.
            e.g. num_img = 2 or [1,1]
        none_as_nan : `bool`
            if True, replace None with np.nan in the lensed_param dictionary.
            default none_as_nan = True.
        lenstype : `str`
            lens type.
            e.g. lenstype = 'I' or 'II' or 'any'.
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

        if not resume:
            os.remove(json_file) 

        n = 0
        while (n < size):
            # disable print statements
            with contextlib.redirect_stdout(None):
                self.sampling_routine(nsamples=self.batch_size,file_name=json_file,resume=False)
                self.lensed_param = 'default'

                lensed_param = dict_list_to_ndarray(self.lensed_param.copy())
                # Dimensions are (nsamples, n_max_images)
                snr = np.array(lensed_param['opt_snr_net'])
                size = len(snr)

                # lens type selection dictionary
                lens_type_dict = {'I': [0, 1], 'II': [2, 3]}
                if lenstype == 'I':
                    snr = snr[:, lens_type_dict[lenstype]
                              [0]:lens_type_dict[lenstype][1]+1]
                elif lenstype == 'II':
                    snr = snr[:, lens_type_dict[lenstype]
                              [0]:lens_type_dict[lenstype][1]+1]
                elif lenstype == 'any':
                    pass
                else:
                    print('lens type not found. Please choose from',
                          lens_type_dict.keys())
                    return None
                
                # dealing with snr_threshold and num_img
                snr_threshold, num_img = np.array(
                    [snr_threshold]).reshape(-1), np.array([num_img]).reshape(-1)  # convert to array
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

                weights = lensed_param['weights']
                # rejection sample wrt to weights
                not_rejected = np.random.uniform(0,1,size)<weights
                snr_hit = snr_hit&not_rejected
                
                # store all params in json file
                if none_as_nan == True:
                    for key, value in lensed_param.items():
                        lensed_param[key] = value[snr_hit]
                else:
                    for key, value in lensed_param.items():
                        lensed_param[key] = np.nan_to_num(np.array(value[snr_hit]))

                append_json(jsonfile, lensed_param, replace=False)        

                n += np.sum(snr_hit) 
            print('collected number of events = ', n)

        # trim the final param dictionary
        self.lensed_param_detectable = 'default'
        param_final = self.lensed_param_detectable.copy()
        self.lensed_param_detectable = None  # to free up memory
        param_final = trim_dictionary(param_final, size)

        return param_final


    def rate_comparision(self, size=False, snr_threshold=8., num_img=2, jsonfile=True, none_as_nan=True):
        '''
        Function to compare the detectable lensed merger rate with the unlensed merger rate
        Intput Parameters:
            size (int): number of samples
            snr_threshold (float/array): threshold for detection signal to noise ratio
            num_img (int/array): number of images
                                e.g. For Sub-thershold events, snr_threshold=[8.,6.], num_img=[1,1]
                                The event will contain 1 image with snr>8 and 1 image with snr>6
            jsonfile (bool): if True, store all gravitational waves source parameters in json file
            none_as_nan (bool): if True,  no value is kept as np.nan
                                if False, no value is kept as 0.
        Output Parameters:
            unlened_rate (float): unlensed merger rate in yr^-1
            lensed_rate (float): lensed merger rate in yr^-1
            rate_ratio (float): lensed/unlensed merger rate ratio
        '''
        # calculate unlensed rate
        if self.gw_param == False:
            print('unlensed_param not sampled beforehand. Sampling now...')
            self.unlensed_cbc_statistics(nsamples=size, jsonfile=jsonfile)
        unlensed_rate = self.unlensed_rate(
            size=size, snr_threshold=np.max([snr_threshold]), jsonfile=jsonfile)

        # calculate lensed rate
        if self.lensed_param == False:
            print('lensed_param not sampled beforehand. Sampling now...')
            self.lensed_cbc_statistics(nsamples=size, jsonfile=jsonfile)
        lensed_rate = self.lensed_rate(
            size=size, snr_threshold=snr_threshold, num_img=num_img, jsonfile=jsonfile)

        rate_ratio = (unlensed_rate[0]/lensed_rate[0],
                      unlensed_rate[1]/lensed_rate[1])
        print('unlensed/lensed rate ratio = ', rate_ratio)
        return (unlensed_rate, lensed_rate, rate_ratio)
