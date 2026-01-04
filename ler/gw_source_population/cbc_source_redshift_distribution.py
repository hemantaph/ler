# -*- coding: utf-8 -*-
"""
This module contains the classes to generate source galaxy population model
and compact binary population model. The compact binary population model
inherits from the source galaxy population model. The source galaxy population
model is used to generate the source galaxy population model. The compact binary
population model is used to generate the compact binary population model.
"""

import warnings
warnings.filterwarnings("ignore")

import numpy as np
from multiprocessing import Pool
from tqdm import tqdm
from scipy.interpolate import CubicSpline

# for redshift to luminosity distance conversion
from astropy.cosmology import LambdaCDM

from ..utils import  FunctionConditioning, interpolator_json_path, luminosity_distance, differential_comoving_volume
from .jit_functions import merger_rate_density_bbh_popI_II_oguri2018, sfr_madau_dickinson2014, merger_rate_density_bbh_popIII_ken2022, merger_rate_density_bbh_primordial_ken2022
from .sfr_with_time_delay import sfr_with_time_delay


class CBCSourceRedshiftDistribution(object):
    """Class to generate a population of source galaxies.
    This class is inherited by :class:`~ler.ler.CBCSourceParameterDistribution` and :class:`~ler.ler.LensGalaxyParameterDistribution` class.

    Parameters
    ----------
    z_min : `float`
        Minimum redshift of the source population
        default: 0.
    z_max : `float`
        Maximum redshift of the source population
        default: 10.
    event_type : `str`
        Type of event to generate.
        e.g. 'BBH', 'BNS', 'NSBH'
    cosmology : `astropy.cosmology`
        Cosmology to use
        default: None/astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)
    merger_rate_density : `str` or `function`
        Type of merger rate density function to use
        default: 'merger_rate_density_popI_II_oguri2018'
        for others see instance method in :class:`~ler.ler.merger_rate_density_model_list`
    merger_rate_density_param : `dict`
        Dictionary of merger rate density function parameters
        default: None/dict(R0=25 * 1e-9, b2=1.6, b3=2.1, b4=30)
    directory : `str`
        Directory to store the interpolator pickle files
        default: './interpolator_json'
    create_new_interpolator : `dict`
        Dictionary of interpolator creation parameters
        default: None/dict(redshift_distribution=dict(create_new=False, resolution=1000), z_to_luminosity_distance=dict(create_new=False, resolution=1000), differential_comoving_volume=dict(create_new=False, resolution=1000))

    Examples
    ----------
    >>> from ler.gw_source_population import CBCSourceRedshiftDistribution
    >>> cbc = CBCSourceRedshiftDistribution(z_min=0.001, z_max=10, merger_rate_density="merger_rate_density_bbh_popI_II_oguri2018")
    >>> cbc.merger_rate_density(zs=0.0001) # local merger rate density at low redshift

    Instance Attributes
    ----------
    SourceGalaxyPopulationModel has the following instance attributes:\n
    +-------------------------------------+----------------------------------+
    | Atrributes                          | Type                             |
    +=====================================+==================================+
    |:attr:`~z_min`                       | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`~z_max`                       | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`~directory`                   | `str`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~event_type`                  | `str`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~cosmo`                       | `astropy.cosmology`              |
    +-------------------------------------+----------------------------------+
    |:attr:`~merger_rate_density_param`   | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~normalization_pdf_z`         | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`~merger_rate_density_model_list`                                 |
    +-------------------------------------+----------------------------------+
    |                                     | List of available                |
    |                                     | merger rate density functions    |
    |                                     | and its parameters               |
    +-------------------------------------+----------------------------------+
    |:attr:`~merger_rate_density`         | `class object`                   |
    +-------------------------------------+----------------------------------+
    |:attr:`~source_redshift`             | `class object`                   |
    +-------------------------------------+----------------------------------+
    |:attr:`~luminosity_distance`         | `class object`                   |
    +-------------------------------------+----------------------------------+
    |:attr:`~differential_comoving_volume`| `class object`                   |
    +-------------------------------------+----------------------------------+


    Instance Methods
    ----------
    SourceGalaxyPopulationModel has the following instance methods:\n
    +-------------------------------------+----------------------------------+
    | Methods                             | Type                             |
    +=====================================+==================================+
    |:meth:`~merger_rate_density_detector_frame`                             |
    +-------------------------------------+----------------------------------+
    |                                     | Function to compute the merger   |
    |                                     | rate density (detector frame)    |
    +-------------------------------------+----------------------------------+
    |:meth:`~merger_rate_density_bbh_popI_II_oguri2018`                      |
    +-------------------------------------+----------------------------------+
    |                                     | Function to compute the merger   |
    |                                     | rate density (PopI/PopII)        |
    |                                     | from Oguri et al. (2018)         |
    +-------------------------------------+----------------------------------+
    |:meth:`~sfr_madau_dickinson2014`                        |
    +-------------------------------------+----------------------------------+
    |                                     | Function to compute star         |
    |                                     | formation rate as given in       |
    |                                     | Eqn. 15 Madau & Dickinson (2014) |
    +-------------------------------------+----------------------------------+
    |:meth:`~merger_rate_density_bbh_popIII_ken2022`                         |
    +-------------------------------------+----------------------------------+
    |                                     | Function to compute the merger   |
    |                                     | rate density (PopIII)            |
    +-------------------------------------+----------------------------------+
    |:meth:`~merger_rate_density_bbh_primordial_ken2022`                     |
    +-------------------------------------+----------------------------------+
    |                                     | Function to compute the merger   |
    |                                     | rate density (Primordial)        |
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
    Cosmology to use for the redshift distribution. \n
    e.g. Planck18, WMAP9, FlatLambdaCDM(H0=70, Om0=0.3) etc.
    """

    merger_rate_density_param = None
    """``dict`` \n
    Dictionary of merger rate density function input parameters
    """

    create_new_interpolator = None
    """``dict`` \n
    Dictionary of interpolator creation parameters. \n
    e.g. dict(redshift_distribution=dict(create_new=False, resolution=1000), z_to_luminosity_distance=dict(create_new=False, resolution=1000), differential_comoving_volume=dict(create_new=False, resolution=1000))
    """

    normalization_pdf_z = None
    """``float`` \n
    Normalization constant of the pdf p(z)
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
        print("\nInitializing CBCSourceRedshiftDistribution...\n")
        # set attributes
        self.npool = npool
        self.z_min = z_min
        self.z_max = z_max
        self.directory = directory
        self.event_type = event_type
        # if None is passed, use the default cosmology
        self.cosmo = cosmology if cosmology else LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

        # setting up the interpolator creation parameters
        self.create_new_interpolator = self.setup_decision_dictionary(create_new_interpolator, merger_rate_density)

        # Initialize cosmological functions
        self.luminosity_distance = luminosity_distance(z_min=self.z_min, z_max=self.z_max, cosmo=self.cosmo, directory=self.directory, create_new=self.create_new_interpolator["luminosity_distance"]["create_new"], resolution=self.create_new_interpolator["luminosity_distance"]["resolution"], get_attribute=True)
        self.differential_comoving_volume = differential_comoving_volume(z_min=self.z_min, z_max=self.z_max, cosmo=self.cosmo, directory=self.directory, create_new=self.create_new_interpolator["differential_comoving_volume"]["create_new"], resolution=self.create_new_interpolator["differential_comoving_volume"]["resolution"], get_attribute=True)

        # function initialization
        merger_rate_density, self.merger_rate_density_param = self.merger_rate_density_priors_categorization(event_type, merger_rate_density, merger_rate_density_param)
        self.merger_rate_density = merger_rate_density # this is an initialization, not a variable assignment

        # source redshift distribution initialization
        self.merger_rate_density_detector_frame = self.merger_rate_density_detector_frame(zs=None, get_attribute=True)
        self.source_redshift = FunctionConditioning(
            function=self.merger_rate_density_detector_frame.z_array,
            x_array=self.merger_rate_density_detector_frame.x_array,
            identifier_dict=self.merger_rate_density_detector_frame.info,
            directory=self.directory,
            sub_directory="source_redshift",
            name="source_redshift",
            create_new=self.create_new_interpolator["merger_rate_density"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='rvs',
        )

        # Normalization of the pdf p(z)
        self.normalization_pdf_z = self.source_redshift.pdf_norm_const

    def setup_decision_dictionary(self, create_new_interpolator, merger_rate_density):
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
            if merger_rate_density=="sfr_with_td":
                create_new_interpolator_["merger_rate_density"]["resolution"] = 48

        return create_new_interpolator_
    
    def merger_rate_density_priors_categorization(self, event_type, merger_rate_density, merger_rate_density_param):
        """
        Function to categorize the merger rate density and its parameters.

        Parameters
        ----------
        event_type : `str`
            Type of event to generate.
            e.g. 'BBH', 'BNS', 'BBH_popIII', 'BBH_primordial', 'NSBH'
        merger_rate_density : `str` or `callable`
            Merger rate density function name or function itself.
            If `str`, it must be one of the available merger rate density functions.
            If `callable`, it must accept a single argument, the redshift.
        merger_rate_density_param : `dict`
            Dictionary of merger rate density function parameters.
            If `None`, use the default parameters for the chosen merger rate density function.
            If not `None`, must contain the following parameters:
                R0 : `float`
                    Normalization constant of the merger rate density.
                b2 : `float`
                    Power law exponent of the merger rate density.
                b3 : `float`
                    Power law exponent of the merger rate density.
                b4 : `float`
                    Power law exponent of the merger rate density.

        Returns
        -------
        merger_rate_density_ : `str` or `callable`
            Merger rate density function name or function itself.
        merger_rate_density_param_ : `dict`
            Dictionary of merger rate density function parameters.

        Notes
        -----
        If `merger_rate_density` is a string, it must be one of the available merger rate density functions.
        If `merger_rate_density` is a callable, it must accept a single argument, the redshift.
        If `merger_rate_density_param` is `None`, use the default parameters for the chosen merger rate density function.
        If `merger_rate_density_param` is not `None`, it must contain the following parameters: R0, b2, b3, b4.
        """

        # update the merger rate density and its parameters with user provided values
        if isinstance(merger_rate_density, str):
            if merger_rate_density in self.merger_rate_density_model_list:
                merger_rate_density_ = merger_rate_density
                # you can't provide merger_rate_density_param and not merger_rate_density
                merger_rate_density_param_ = self.merger_rate_density_model_list[merger_rate_density]
            else:
                raise ValueError(f"'merger rate density' sampler '{merger_rate_density}' not available.\n Available 'merger rate density' samplers and its parameters are: {self.merger_rate_density_model_list}")
        elif callable(merger_rate_density):
            print("using user provided custom merger rate density function")
            merger_rate_density_ = merger_rate_density
        else:
            merger_rate_density_ = "merger_rate_density_bbh_popI_II_oguri2018"

        # Oguri et al. (2018) merger rate density
        if merger_rate_density_ == "merger_rate_density_bbh_popI_II_oguri2018":
            if event_type == "BBH":
                merger_rate_density_param_ = dict(R0=23.9 * 1e-9, b2=1.6, b3=2.1, b4=30)
            elif event_type == "BNS":
                merger_rate_density_param_ = dict(R0=105.5 * 1e-9, b2=1.6, b3=2.1, b4=30)
            elif event_type == "NSBH":
                merger_rate_density_param_ = dict(R0=45.0 * 1e-9, b2=1.6, b3=2.1, b4=30)
            else:
                raise ValueError("event_type must be 'BBH', 'BNS' or 'NSBH'")
        
        # merger rate with time delay; Borhanian & Sathyaprakash (2024)
        list_ = ['sfr_with_td',
                  #'sfr_madau_fragos2017_with_bbh_dt', 
                  # 'sfr_madau_dickinson2014_with_bbh_dt',
                  # 'sfr_madau_fragos2017_with_bns_dt', 
                  # 'sfr_madau_dickinson2014_with_bns_dt'
                ]
        if merger_rate_density_ in list_:
            if event_type == "BBH":
                merger_rate_density_param_ = dict(R0=23.9 * 1e-9, a=0.01, b=2.6, c=3.2, d=6.2, td_min=10e-3, td_max=10.0)
            elif event_type == "BNS":
                merger_rate_density_param_ = dict(R0=105.5 * 1e-9, a=0.01, b=2.6, c=3.2, d=6.2, td_min=20e-3, td_max=10.0)
            else:
                raise ValueError("event_type must be 'BBH', 'BNS'")
            
        # if callable(merger_rate_density), below code will not matter
        if isinstance(merger_rate_density_param, dict):  # merger_rate_density_param is user provided
            merger_rate_density_param_.update(merger_rate_density_param)

        return merger_rate_density_, merger_rate_density_param_

    def merger_rate_density_detector_frame(self, zs, get_attribute=False, **kwargs):
        """
        Function to compute the merger rate density (detector frame). The output is in detector frame and is unnormalized.

        Parameters
        ----------
        zs : `float` or `numpy.ndarray` (1D array of floats)
            Source redshifts
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. if the merger_rate_density is merger_rate_density_bbh_popI_II_oguri2018
            param = dict(R0=23.9*1e-9, b2=1.6, b3=2.1, b4=30)

        Returns
        ----------
        rate_density : `numpy.ndarray`
            1D array of floats
            merger rate density (detector frame) (Mpc^-3 yr^-1)

        Examples
        ----------
        >>> from ler.gw_source_population import SourceGalaxyPopulationModel
        >>> cbc = SourceGalaxyPopulationModel()
        >>> rate_density = cbc.merger_rate_density_detector_frame(zs=0.1)
        """

        identifier_dict = self.merger_rate_density_param.copy()
        identifier_dict['z_min'] = self.z_min
        identifier_dict['z_max'] = self.z_max
        identifier_dict['cosmology'] = self.cosmo
        identifier_dict['event_type'] = self.event_type
        identifier_dict['name'] = "merger_rate_density_detector_frame"
        identifier_dict['resolution'] = self.create_new_interpolator["merger_rate_density"]["resolution"]

        zs = np.linspace(identifier_dict['z_min'], identifier_dict['z_max'], identifier_dict['resolution'])
        Pzs = lambda z: self.merger_rate_density(z)/(1+z) * self.differential_comoving_volume(z)

        Pzs_object = FunctionConditioning(
            function=Pzs,
            x_array=zs,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="merger_rate_density",
            name=identifier_dict['name'],
            create_new=self.create_new_interpolator["merger_rate_density"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='function',
        )

        if get_attribute:
            return Pzs_object
        else:
            return Pzs_object(zs)

    def merger_rate_density_bbh_popI_II_oguri2018(self, zs, get_attribute=False, **kwargs):
        """
        Function to compute the merger rate density (PopI/PopII). Reference: Oguri et al. (2018). The output is in source frame and is unnormalized.

        Parameters
        ----------
        zs : `float` or `numpy.ndarray` (nD array of floats)
            Source redshifts
        get_attribute : `bool`
            If True, returns the merger rate density function instead of the value
            default: False
        kwargs : `dict`
            Dictionary of merger rate density function fitting parameters. 
            default: R0=23.9*1e-9, b2=1.6, b3=2.1, b4=30
            R0 is the local merger rate density at low redshift in Mpc^-3 yr^-1

        Returns
        ----------
        rate_density : `float` or `numpy.ndarray` (nD array of floats)
            merger rate density in source frame (Mpc^-3 yr^-1)

        Examples
        ----------
        >>> from ler.gw_source_population import SourceGalaxyPopulationModel
        >>> pop = SourceGalaxyPopulationModel(z_min=0.0, z_max=10, event_type = "BBH", merger_rate_density="merger_rate_density_bbh_popI_II_oguri2018")
        >>> rate_density = pop.merger_rate_density(zs=0.1)
        """

        identifier_dict = {}
        identifier_dict['z_min'] = self.z_min
        identifier_dict['z_max'] = self.z_max
        identifier_dict['cosmology'] = self.cosmo
        identifier_dict['event_type'] = self.event_type
        identifier_dict['name'] = "merger_rate_density_bbh_popI_II_oguri2018"
        identifier_dict['resolution'] = self.create_new_interpolator["merger_rate_density"]["resolution"]
        param_dict = self.merger_rate_density_model_list["merger_rate_density_bbh_popI_II_oguri2018"]
        param_dict.update(kwargs)
        identifier_dict.update(param_dict)
        
        zs_array = np.linspace(identifier_dict['z_min'], identifier_dict['z_max'], identifier_dict['resolution'])
        Rzs = lambda zs: merger_rate_density_bbh_popI_II_oguri2018(
            zs=zs, 
            R0=identifier_dict['R0'], 
            b2=identifier_dict['b2'], 
            b3=identifier_dict['b3'], 
            b4=identifier_dict['b4']
        )

        Rzs_object = FunctionConditioning(
            function=Rzs,
            x_array=zs_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="merger_rate_density",
            name=identifier_dict['name'],
            create_new=self.create_new_interpolator["merger_rate_density"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='function',
        )

        return Rzs_object if get_attribute else Rzs_object(zs)
    
    def sfr_with_td(self, zs, get_attribute=False, **kwargs):
        """
        """

        identifier_dict = {}
        identifier_dict['z_min'] = self.z_min
        identifier_dict['z_max'] = self.z_max
        identifier_dict['cosmology'] = self.cosmo
        identifier_dict['event_type'] = self.event_type
        identifier_dict['name'] = "sfr_with_td"
        identifier_dict['resolution'] = self.create_new_interpolator["merger_rate_density"]["resolution"]
        param_dict = self.merger_rate_density_model_list["sfr_with_td"]
        param_dict.update(kwargs)
        identifier_dict.update(param_dict)
        print(identifier_dict)
    
        print("Numerically solving the merger_rate_density with time delay")
        zs_resolution = identifier_dict['resolution']
        zs_array = np.geomspace(self.z_min+0.001, self.z_max, zs_resolution) if self.z_min==0 else np.geomspace(self.z_min, self.z_max, zs_resolution)

        _, it_exist = interpolator_json_path(
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="merger_rate_density",
            interpolator_name=identifier_dict['name'],
        )

        create_new = self.create_new_interpolator["merger_rate_density"]["create_new"]
        if not it_exist or create_new:
            rate_density = self._helper_rate_density_multiprocessing(zs_array, identifier_dict)
        else:
            rate_density=None

        rate_density_object = FunctionConditioning(
            function=rate_density,
            x_array=zs_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="merger_rate_density",
            name=identifier_dict['name'],
            create_new=create_new,
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='function',
        )

        return rate_density_object if get_attribute else rate_density_object(zs)
    
    def _helper_rate_density_multiprocessing(self, zs_array, identifier_dict):
        """
        Computes the merger rate density distribution with a time delay applied to the star formation rate (SFR),
        utilizing multiprocessing for improved performance.

        Parameters
        ----------
        zs_array : `numpy.ndarray`
            1D array of source redshifts for which the merger rate density is calculated.
        identifier_dict : `dict`
            Dictionary containing parameters for the computation, including cosmological parameters
            (H0, Omega_M, Omega_Lambda), time delay bounds (td_min, td_max), and SFR parameters (a, b, c, d).

        Returns
        -------
        rate_density_array : `numpy.ndarray`
            1D array representing the computed merger rate density distribution normalized to the local
            merger rate density at low redshift (R0).
        """

        size = len(zs_array)
        input_args = np.array([
            zs_array, # source redshifts
            np.arange(size), # index
            identifier_dict['td_min']*np.ones(size), # time delay minimum
            identifier_dict['td_max']*np.ones(size), # time delay maximum
            identifier_dict['cosmology'].H0.value*np.ones(size), # H0
            identifier_dict['cosmology'].Om0*np.ones(size), # Omega_M
            identifier_dict['cosmology'].Ode0*np.ones(size), # Omega_Lambda
            # SFR params
            identifier_dict['a']*np.ones(size),
            identifier_dict['b']*np.ones(size),
            identifier_dict['c']*np.ones(size),
            identifier_dict['d']*np.ones(size),
        ]).T

        print("Computing merger rate density distribution (with time delay to SFR) using multiprocessing...")
        # with tqdm
        rate_density_array = np.zeros(size)
        with Pool(processes=self.npool) as pool:            
            for result in tqdm(
                pool.imap_unordered(sfr_with_time_delay, input_args),
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
        rate_density_array = rate_density_array/rm_spline(0.) * identifier_dict['R0']
    
        return rate_density_array

    def sfr_madau_dickinson2014(self, zs, get_attribute=False, **kwargs):
        """
        Formation rate as given in Eqn. 15 Madau & Dickinson (2014). The output is in detector frame and is unnormalized.

        Parameters
        ----------
        zs : `float` or `numpy.ndarray` (nd.array of floats)
            Source redshifts
        get_attribute : `bool`
            If True, returns the merger rate density function instead of the value
            default: False
        kwargs : `dict`
            Dictionary of star formation rate function fitting parameters.
            default: af=2.7, bf=5.6, cf=2.9
            
        Returns
        ----------
        rate_density : `float` or `numpy.ndarray` (nD array of floats)
            merger rate density in detector frame (Mpc^-3 yr^-1)

        Examples
        ----------
        >>> from ler.gw_source_population import SourceGalaxyPopulationModel
        >>> pop = SourceGalaxyPopulationModel(z_min=5., z_max=40., event_type = "BBH", merger_rate_density="sfr_madau_dickinson2014")
        >>> rate_density = pop.merger_rate_density(zs=10)
        """

        identifier_dict = {}
        identifier_dict['z_min'] = self.z_min
        identifier_dict['z_max'] = self.z_max
        identifier_dict['cosmology'] = self.cosmo
        identifier_dict['event_type'] = self.event_type
        identifier_dict['name'] = "sfr_madau_dickinson2014"
        identifier_dict['resolution'] = self.create_new_interpolator["merger_rate_density"]["resolution"]
        param_dict = self.merger_rate_density_model_list["sfr_madau_dickinson2014"]
        param_dict.update(kwargs)
        identifier_dict.update(param_dict)
        
        zs_array = np.linspace(identifier_dict['z_min'], identifier_dict['z_max'], identifier_dict['resolution'])
        Rzs = lambda zs: sfr_madau_dickinson2014(
            zs=zs, 
            a=identifier_dict['a'], 
            b=identifier_dict['b'], 
            c=identifier_dict['c'],
        )

        Rzs_object = FunctionConditioning(
            function=Rzs,
            x_array=zs_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="merger_rate_density",
            name=identifier_dict['name'],
            create_new=self.create_new_interpolator["merger_rate_density"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='function',
        )

        if get_attribute:
            return Rzs_object
        else:
            return Rzs_object(zs)

    def merger_rate_density_bbh_popIII_ken2022(self, zs, get_attribute=False, **kwargs):
        """
        Merger rate density (PopIII). Reference: Ng et al. 2022. The output is in detector frame and is unnormalized.

        Parameters
        ----------
        zs : `float` or `numpy.ndarray` 
            Source redshifts
        get_attribute : `bool`
            If True, returns the merger rate density function instead of the value
            default: False
        kwargs : `dict`
            Dictionary of merger rate density function fitting parameters.
            default: n0=19.2*1e-9, aIII=0.66, bIII=0.3, zIII=11.6
            n0 is the local merger rate density at low redshift in Mpc^-3 yr^-1

        Returns
        ----------
        rate_density : `float` or `numpy.ndarray`
            merger rate density in detector frame (Mpc^-3 yr^-1)

        Examples
        ----------
        >>> from ler.gw_source_population import SourceGalaxyPopulationModel
        >>> pop = SourceGalaxyPopulationModel(z_min=5, z_max=40, event_type = "BBH", merger_rate_density="merger_rate_density_popIII_ken2022")
        >>> rate_density = pop.merger_rate_density(zs=10)
        >>> rate_density  # Mpc^-3 yr^-1
        1.5107979464621443e-08
        """
        
        identifier_dict = {}
        identifier_dict['z_min'] = self.z_min
        identifier_dict['z_max'] = self.z_max
        identifier_dict['cosmology'] = self.cosmo
        identifier_dict['event_type'] = self.event_type
        identifier_dict['name'] = "merger_rate_density_bbh_popIII_ken2022"
        identifier_dict['resolution'] = self.create_new_interpolator["merger_rate_density"]["resolution"]
        param_dict = self.merger_rate_density_model_list["merger_rate_density_bbh_popIII_ken2022"]
        param_dict.update(kwargs)
        identifier_dict.update(param_dict)

        zs_array = np.linspace(identifier_dict['z_min'], identifier_dict['z_max'], identifier_dict['resolution'])
        Rzs = lambda zs: merger_rate_density_bbh_popIII_ken2022(
            zs=zs,
            n0=identifier_dict['n0'],
            aIII=identifier_dict['aIII'],
            bIII=identifier_dict['bIII'],
            zIII=identifier_dict['zIII']
        )

        Rzs_object = FunctionConditioning(
            function=Rzs,
            x_array=zs_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="merger_rate_density",
            name=identifier_dict['name'],
            create_new=self.create_new_interpolator["merger_rate_density"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='function',
        )

        if get_attribute:
            return Rzs_object
        else:
            return Rzs_object(zs)

    
    def merger_rate_density_bbh_primordial_ken2022(self, zs, get_attribute=False, **kwargs):
        """
        Function to compute the merger rate density (Primordial). Reference: Ng et al. 2022. The output is in detector frame and is unnormalized.

        Parameters
        ----------
        zs : `float` or `numpy.ndarray` 
            Source redshifts
        n0 : `float`
            normalization constant
            default: 0.044*1e-9
        t0 : `float`
            Present age of the Universe in Gyr
            default: 13.786885302009708
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. param = dict(t0=13.786885302009708)

        Returns
        ----------
        rate_density : `float`
            merger rate density

        Examples
        ----------
        >>> from ler.gw_source_population import SourceGalaxyPopulationModel
        >>> pop = SourceGalaxyPopulationModel(z_min=5, z_max=40, event_type = "BBH", merger_rate_density="merger_rate_density_bbh_primordial_ken2022")
        >>> rate_density = pop.merger_rate_density(zs=10)
        >>> rate_density  # Mpc^-3 yr^-1
        9.78691173794454e-10
        """

        identifier_dict = {}
        identifier_dict['z_min'] = self.z_min
        identifier_dict['z_max'] = self.z_max
        identifier_dict['cosmology'] = self.cosmo
        identifier_dict['event_type'] = self.event_type
        identifier_dict['name'] = "merger_rate_density_bbh_primordial_ken2022"
        identifier_dict['resolution'] = self.create_new_interpolator["merger_rate_density"]["resolution"]
        param_dict = self.merger_rate_density_model_list["merger_rate_density_bbh_primordial_ken2022"]
        param_dict.update(kwargs)
        identifier_dict.update(param_dict)

        zs_array = np.linspace(identifier_dict['z_min'], identifier_dict['z_max'], identifier_dict['resolution'])
        Rzs = lambda zs: merger_rate_density_bbh_primordial_ken2022(
            zs=zs,
            n0=identifier_dict['n0'],
            t0=identifier_dict['t0']
        )

        Rzs_object = FunctionConditioning(
            function=Rzs,
            x_array=zs_array,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="merger_rate_density",
            name=identifier_dict['name'],
            create_new=self.create_new_interpolator["merger_rate_density"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='function',
        )

        if get_attribute:
            return Rzs_object
        else:
            return Rzs_object(zs)
    

    @property
    def merger_rate_density(self):
        """
        Source frame merger rate density function wrt redshift.

        Parameters
        ----------
        zs : `float`
            1D array of floats
            Source redshifts

        Returns
        ----------
        merger_rate_density : `float`
            merger rate density in detector frame (Mpc^-3 yr^-1)

        Examples
        ----------
        >>> from ler.gw_source_population import SourceGalaxyPopulationModel
        >>> cbc = SourceGalaxyPopulationModel()
        >>> merger_rate_density = cbc.merger_rate_density(zs=0.1)
        """
        return self._merger_rate_density
    
    @merger_rate_density.setter
    def merger_rate_density(self, function):

        if function in self.merger_rate_density_model_list:
            print(f"using ler available merger rate density model: {function}")
            args = self.merger_rate_density_param
            if args is None:
                self._merger_rate_density = getattr(self, function)(zs=None, get_attribute=True)
            else:
                self._merger_rate_density = getattr(self, function)(zs=None, get_attribute=True, **args)

        elif callable(function):
            print("using user provided custom merger rate density function")
            self._merger_rate_density = FunctionConditioning(function=None, x_array=None,create_function=function)

        elif isinstance(function, object):
            print("using user provided custom merger rate density class/object")
            self._merger_rate_density = function

        else:
            raise ValueError("merger_rate_density must be a function or a string from the available merger rate density model list")

    @property
    def merger_rate_density_model_list(self):
        """
        Dictionary of available merger rate density functions and its parameters.
        """

        self._merger_rate_density_model_list = dict(
            merger_rate_density_bbh_popI_II_oguri2018=dict(
                R0=23.9 * 1e-9, b2=1.6, b3=2.1, b4=30
            ),
            sfr_madau_fragos2017=dict(
                a=0.01, b=2.6, c=3.2, d=6.2
            ),
            sfr_madau_dickinson2014=dict(
                a=0.015, b=2.7, c=2.9, d=5.6
            ),
            sfr_with_td=dict(
                R0=23.9 * 1e-9, a=0.01, b=2.6, c=3.2, d=6.2, td_min=10e-3, td_max=10.0
            ),
            merger_rate_density_bbh_popIII_ken2022=dict(
                n0=19.2 * 1e-9, aIII=0.66, bIII=0.3, zIII=11.6
            ),
            merger_rate_density_bbh_primordial_ken2022=dict(
                n0=0.044 * 1e-9, t0=13.786885302009708
            ),
            sfr_madau_fragos2017_with_bbh_dt=dict(
                R0=23.9 * 1e-9,
            ),
            sfr_madau_dickinson2014_with_bbh_dt=dict(
                R0=23.9 * 1e-9,
            ),
            sfr_madau_fragos2017_with_bns_dt=dict(
                R0=105.5 * 1e-9,
            ),
            sfr_madau_dickinson2014_with_bns_dt=dict(
                R0=105.5 * 1e-9,
            ),
        )

        return self._merger_rate_density_model_list