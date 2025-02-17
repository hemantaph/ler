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

# for redshift to luminosity distance conversion
from astropy.cosmology import LambdaCDM

from ..utils import  FunctionConditioning
from .jit_functions import merger_rate_density_bbh_popI_II_oguri2018, star_formation_rate_madau_dickinson2014, merger_rate_density_bbh_popIII_ken2022, merger_rate_density_bbh_primordial_ken2022


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
        default: None/dict(R0=25 * 1e-9, b2=1.6, b3=2.0, b4=30)
    directory : `str`
        Directory to store the interpolator pickle files
        default: './interpolator_pickle'
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
    |:meth:`~create_lookup_table`         | Function to create a lookup      |
    |                                     | table for the differential       |
    |                                     | comoving volume and luminosity   |
    |                                     | distance wrt redshift            |
    +-------------------------------------+----------------------------------+
    |:meth:`~merger_rate_density_bbh_popI_II_oguri2018`                      |
    +-------------------------------------+----------------------------------+
    |                                     | Function to compute the merger   |
    |                                     | rate density (PopI/PopII)        |
    |                                     | from Oguri et al. (2018)         |
    +-------------------------------------+----------------------------------+
    |:meth:`~star_formation_rate_madau_dickinson2014`                        |
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
        z_min=0.001,
        z_max=10.0,
        event_type="BBH",
        merger_rate_density=None,
        merger_rate_density_param=None,
        cosmology=None,
        directory="./interpolator_pickle",
        create_new_interpolator=False,
    ):
        print("\nInitializing CBCSourceRedshiftDistribution...\n")
        # set attributes
        self.z_min = z_min
        self.z_max = z_max
        self.directory = directory
        self.event_type = event_type
        # if None is passed, use the default cosmology
        self.cosmo = cosmology if cosmology else LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

        # setting up the interpolator creation parameters
        self.create_new_interpolator = self.setup_decision_dictionary(create_new_interpolator)
        
        # creating of interpolators for redshift dependent quantities
        self.create_lookup_table()

        # function initialization
        merger_rate_density, self.merger_rate_density_param = self.merger_rate_density_priors_categorization(event_type, merger_rate_density, merger_rate_density_param)
        self.merger_rate_density = merger_rate_density # this is an initialization, not a variable assignment

        # source redshift distribution initialization
        self.source_redshift = self.merger_rate_density_detector_frame(zs=None, get_attribute=True)
        self.normalization_pdf_z = self.source_redshift.pdf_norm_const


    def setup_decision_dictionary(self, create_new_interpolator):
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
        if (isinstance(merger_rate_density, str)) and (merger_rate_density in self.merger_rate_density_model_list):
            merger_rate_density_ = merger_rate_density
            # you can't provide merger_rate_density_param and not merger_rate_density
            merger_rate_density_param_ = self.merger_rate_density_model_list[merger_rate_density]
        elif callable(merger_rate_density):
            print("using user provided custom merger rate density function")
            merger_rate_density_ = merger_rate_density
        else:
            merger_rate_density_ = "merger_rate_density_bbh_popI_II_oguri2018"

        if merger_rate_density_ == "merger_rate_density_bbh_popI_II_oguri2018":
            if event_type == "BBH":
                merger_rate_density_param_ = dict(R0=23.9 * 1e-9, b2=1.6, b3=2.0, b4=30)
            elif event_type == "BNS":
                merger_rate_density_ = "merger_rate_density_bbh_popI_II_oguri2018"
                merger_rate_density_param_ = dict(R0=105.5 * 1e-9, b2=1.6, b3=2.0, b4=30)
            elif event_type == "NSBH":
                merger_rate_density_ = "merger_rate_density_bbh_popI_II_oguri2018"
                merger_rate_density_param_ = dict(R0=45.0 * 1e-9, b2=1.6, b3=2.0, b4=30)
            else:
                raise ValueError("event_type must be 'BBH', 'BNS' or 'NSBH'")
            
        # if callable(merger_rate_density), below code will not matter
        if isinstance(merger_rate_density_param, dict):
            merger_rate_density_param_ = merger_rate_density_param_.update(merger_rate_density_param)

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
            param = dict(R0=23.9*1e-9, b2=1.6, b3=2.0, b4=30)

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
            param_dict_given=identifier_dict,
            directory=self.directory,
            sub_directory="merger_rate_density",
            name=identifier_dict['name'],
            create_new=self.create_new_interpolator["merger_rate_density"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback='rvs',
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
            default: R0=23.9*1e-9, b2=1.6, b3=2.0, b4=30
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

        identifier_dict = self.merger_rate_density_param.copy()
        identifier_dict['z_min'] = self.z_min
        identifier_dict['z_max'] = self.z_max
        identifier_dict['cosmology'] = self.cosmo
        identifier_dict['event_type'] = self.event_type
        identifier_dict['name'] = "merger_rate_density_bbh_popI_II_oguri2018"
        identifier_dict['resolution'] = self.create_new_interpolator["merger_rate_density"]["resolution"]
        identifier_dict.update(kwargs)
        
        zs = np.linspace(identifier_dict['z_min'], identifier_dict['z_max'], identifier_dict['resolution'])
        Rzs = lambda zs: merger_rate_density_bbh_popI_II_oguri2018(
            zs=zs, 
            R0=identifier_dict['R0'], 
            b2=identifier_dict['b2'], 
            b3=identifier_dict['b3'], 
            b4=identifier_dict['b4']
        )

        Rzs_object = FunctionConditioning(
            function=Rzs,
            x_array=zs,
            param_dict_given=identifier_dict,
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
    
    def star_formation_rate_madau_dickinson2014(self, zs, get_attribute=False, **kwargs):
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
        >>> pop = SourceGalaxyPopulationModel(z_min=5., z_max=40., event_type = "BBH", merger_rate_density="star_formation_rate_madau_dickinson2014")
        >>> rate_density = pop.merger_rate_density(zs=10)
        """

        identifier_dict = self.merger_rate_density_param.copy()
        identifier_dict['z_min'] = self.z_min
        identifier_dict['z_max'] = self.z_max
        identifier_dict['cosmology'] = self.cosmo
        identifier_dict['event_type'] = self.event_type
        identifier_dict['name'] = "star_formation_rate_madau_dickinson2014"
        identifier_dict['resolution'] = self.create_new_interpolator["merger_rate_density"]["resolution"]
        identifier_dict.update(kwargs)
        
        zs = np.linspace(identifier_dict['z_min'], identifier_dict['z_max'], identifier_dict['resolution'])
        Rzs = lambda zs: star_formation_rate_madau_dickinson2014(
            zs=zs, 
            af=identifier_dict['af'], 
            bf=identifier_dict['bf'], 
            cf=identifier_dict['cf'],
        )

        Rzs_object = FunctionConditioning(
            function=Rzs,
            x_array=zs,
            param_dict_given=identifier_dict,
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
        
        identifier_dict = self.merger_rate_density_param.copy()
        identifier_dict['z_min'] = self.z_min
        identifier_dict['z_max'] = self.z_max
        identifier_dict['cosmology'] = self.cosmo
        identifier_dict['event_type'] = self.event_type
        identifier_dict['name'] = "merger_rate_density_bbh_popIII_ken2022"
        identifier_dict['resolution'] = self.create_new_interpolator["merger_rate_density"]["resolution"]
        identifier_dict.update(kwargs)

        zs = np.linspace(identifier_dict['z_min'], identifier_dict['z_max'], identifier_dict['resolution'])
        Rzs = lambda zs: merger_rate_density_bbh_popIII_ken2022(
            zs=zs,
            n0=identifier_dict['n0'],
            aIII=identifier_dict['aIII'],
            bIII=identifier_dict['bIII'],
            zIII=identifier_dict['zIII']
        )

        Rzs_object = FunctionConditioning(
            function=Rzs,
            x_array=zs,
            param_dict_given=identifier_dict,
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

        identifier_dict = self.merger_rate_density_param.copy()
        identifier_dict['z_min'] = self.z_min
        identifier_dict['z_max'] = self.z_max
        identifier_dict['cosmology'] = self.cosmo
        identifier_dict['event_type'] = self.event_type
        identifier_dict['name'] = "merger_rate_density_bbh_primordial_ken2022"
        identifier_dict['resolution'] = self.create_new_interpolator["merger_rate_density"]["resolution"]
        identifier_dict.update(kwargs)

        zs = np.linspace(identifier_dict['z_min'], identifier_dict['z_max'], identifier_dict['resolution'])
        Rzs = lambda zs: merger_rate_density_bbh_primordial_ken2022(
            zs=zs,
            n0=identifier_dict['n0'],
            t0=identifier_dict['t0']
        )

        Rzs_object = FunctionConditioning(
            function=Rzs,
            x_array=zs,
            param_dict_given=identifier_dict,
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
    
    def create_lookup_table(self):
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

        z_min = 0.001 if self.z_min == 0. else self.z_min
        z_max = self.z_max

        resolution = self.create_new_interpolator["luminosity_distance"]["resolution"]
        create_new = self.create_new_interpolator["luminosity_distance"]["create_new"]
        zs = np.geomspace(0.001, z_max, 1000)
        Dl = self.cosmo.luminosity_distance(zs).value
        
        self.luminosity_distance = FunctionConditioning(
            function=Dl,
            x_array=zs,
            conditioned_y_array=None,
            param_dict_given=dict(z_min=z_min, z_max=z_max, cosmology=self.cosmo, resolution=resolution, details="luminosity_distance from astropy.cosmology"),
            directory=self.directory,
            sub_directory="luminosity_distance",
            name="luminosity_distance",
            create_new=create_new,
            create_function_inverse=True,
            create_function=True,
            create_pdf=False,
            create_rvs=False,
            callback='function',
        )
        self.luminosity_distance.__doc__ = """
        Redshift to luminosity distance conversion.

        Parameters
        ----------
        zs : `numpy.ndarray` or `float`
            Source redshifts

        Returns
        ----------
        luminosity_distance : `numpy.ndarray`
            luminosity distance in Mpc

        Examples
        ----------
        >>> from ler.gw_source_population import SourceGalaxyPopulationModel
        >>> ler = SourceGalaxyPopulationModel()  # with default LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        >>> luminosity_distance = ler.luminosity_distance(1.)
        >>> luminosity_distance = ler.luminosity_distance.function(np.array([1., 2.]))
        >>> redshift = ler.luminosity_distance.function_inverse(np.array([100., 200.]))
        """

        # get differential co-moving volume interpolator
        resolution = self.create_new_interpolator["differential_comoving_volume"]["resolution"]
        create_new = self.create_new_interpolator["differential_comoving_volume"]["create_new"]
        zs = np.geomspace(z_min, z_max, resolution)
        dVcdz = self.cosmo.differential_comoving_volume(zs).value * 4 * np.pi  # volume of shell in Mpc^3
        self.differential_comoving_volume = FunctionConditioning(
            function=dVcdz,
            x_array=zs,
            conditioned_y_array=None,
            param_dict_given=dict(z_min=z_min, z_max=z_max, cosmology=self.cosmo, resolution=resolution, details="differential_comoving_volume from astropy.cosmology"),
            directory=self.directory,
            sub_directory="differential_comoving_volume",
            name="differential_comoving_volume",
            create_new=create_new,
            create_function_inverse=False,
            create_function=True,
            create_pdf=False,
            create_rvs=False,
            callback='function',
        )
        self.differential_comoving_volume.__doc__ = """
        Redshift to differential comoving volume conversion.

        Parameters
        ----------
        zs : `numpy.ndarray` or `float`
            Source redshifts

        Returns
        ----------
        differential_comoving_volume : `numpy.ndarray`
            differential comoving volume in Mpc^3

        Examples
        ----------
        >>> from ler.len_galaxy_population import OpticalDepth
        >>> ler = OpticalDepth()  # with default LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        >>> differential_comoving_volume = ler.differential_comoving_volume(1.)
        >>> differential_comoving_volume = ler.differential_comoving_volume.function(np.array([1., 2.]))
        """

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
            self._merger_rate_density = getattr(self, function)(zs=None, get_attribute=True)
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
                R0=23.9 * 1e-9, b2=1.6, b3=2.0, b4=30
            ),
            star_formation_rate_madau_dickinson2014=dict(af=2.7, bf=5.6, cf=2.9),
            merger_rate_density_bbh_popIII_ken2022=dict(
                n0=19.2 * 1e-9, aIII=0.66, bIII=0.3, zIII=11.6
            ),
            merger_rate_density_bbh_primordial_ken2022=dict(
                n0=0.044 * 1e-9, t0=13.786885302009708
            ),
        )

        return self._merger_rate_density_model_list
