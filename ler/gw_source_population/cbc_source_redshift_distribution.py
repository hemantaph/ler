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
from numba import njit

# from gwcosmo import priors as p
from scipy.integrate import quad

# for redshift to luminosity distance conversion
from astropy.cosmology import LambdaCDM

from ..utils import  interpolator_from_pickle, cubic_spline_interpolator, inverse_transform_sampler
from .jit_functions import merger_rate_density_bbh_popI_II_oguri2018, star_formation_rate_madau_dickinson2014, merger_rate_density_bbh_popIII_ken2022, merger_rate_density_bbh_primordial_ken2022


class CBCSourceRedshiftDistribution(object):
    """Class to generate a population of source galaxies.
    This class is inherited by :class:`~ler.ler.CompactBinaryPopulation` and :class:`~ler.ler.LensGalaxyParameterDistribution` class.

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
        default: None/dict(redshift_distribution=dict(create_new=False, resolution=500), z_to_luminosity_distance=dict(create_new=False, resolution=500), differential_comoving_volume=dict(create_new=False, resolution=500))

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
    |:attr:`~sample_source_redshift`      | Function to sample source        |
    |                                     | redshifts (source frame)         |
    +-------------------------------------+----------------------------------+

    Instance Methods
    ----------
    SourceGalaxyPopulationModel has the following instance methods:\n
    +-------------------------------------+----------------------------------+
    | Methods                             | Type                             |
    +=====================================+==================================+
    |:attr:`~merger_rate_density`         | `function`                       |
    +-------------------------------------+----------------------------------+
    |:attr:`~z_to_luminosity_distance`    | `function`                       |
    +-------------------------------------+----------------------------------+
    |:attr:`~differential_comoving_volume`| `function`                       |
    +-------------------------------------+----------------------------------+
    |:meth:`~pdf_z`                       | Function to compute the pdf      |
    |                                     | p(z)                             |
    +-------------------------------------+----------------------------------+
    |:meth:`~merger_rate_density_src_frame`                                  |
    +-------------------------------------+----------------------------------+
    |                                     | Function to compute the merger   |
    |                                     | rate density (source frame)      |
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

    c_n_i = None
    """``dict`` \n
    c_n_i stands for 'create new interpolator'. Dictionary of interpolator creation parameters. \n
    e.g. dict(redshift_distribution=dict(create_new=False, resolution=500), z_to_luminosity_distance=dict(create_new=False, resolution=500), differential_comoving_volume=dict(create_new=False, resolution=500))
    """

    normalization_pdf_z = None
    """``float`` \n
    Normalization constant of the pdf p(z)
    """

    z_to_luminosity_distance = None
    """``scipy.interpolate.interpolate`` \n
    Function to convert redshift to luminosity distance
    """

    differential_comoving_volume = None
    """``scipy.interpolate.interpolate`` \n
    Function to calculate the differential comoving volume
    """

    def __init__(
        self,
        z_min=0.001,
        z_max=10.0,
        event_type="BBH",
        merger_rate_density="merger_rate_density_bbh_popI_II_oguri2018",
        merger_rate_density_param=None,
        cosmology=None,
        directory="./interpolator_pickle",
        create_new_interpolator=False,
    ):
        # set attributes
        self.z_min = z_min
        self.z_max = z_max
        self.event_type = event_type
        # if None is passed, use the default cosmology
        self.cosmo = cosmology if cosmology else LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        # setting up the interpolator creation parameters
        self.c_n_i = dict(
            redshift_distribution=dict(create_new=False, resolution=500), z_to_luminosity_distance=dict(create_new=False, resolution=500), differential_comoving_volume=dict(create_new=False, resolution=500))
        if create_new_interpolator:
            self.c_n_i.update(create_new_interpolator)
        # creating of interpolators for redshift dependent quantities
        self.create_lookup_table(z_min, z_max, directory)

        # Define the merger-rate density function/method instances
        # self.merger_rate_density takes function or function name as input 
        self.merger_rate_density = merger_rate_density  # function initialization
        # if None is passed, use the default merger_rate_density_param
        if merger_rate_density_param is None and merger_rate_density== "merger_rate_density_bbh_popI_II_oguri2018":
            if self.event_type == "BBH":
                R0 = 25 * 1e-9
            elif self.event_type == "BNS":
                R0 = 170.0 * 1e-9
            elif self.event_type == "NSBH":
                R0 = 27.0 * 1e-9
            else:
                raise ValueError("event_type must be one of 'BBH', 'BNS', 'NSBH'")
            self.merger_rate_density_param = dict(R0=R0, b2=1.6, b3=2.0, b4=30)

        # To find the normalization constant of the pdf p(z)
        # merger_rate_density_src_frame is njit function and takes array as input but quad integrand function takes only float as input
        merger_rate_density_src_frame = lambda z: self.merger_rate_density_src_frame(np.array([z]))[0]
        # this normalization is important to find the correct pdf p(z)
        self.normalization_pdf_z = quad(
            merger_rate_density_src_frame,
            z_min,
            z_max
        )[0]

        # generate inverse cdf for inverse transform sampling
        # create sampler using the pdf p(z)
        resolution = self.c_n_i["redshift_distribution"]["resolution"]
        create_new = self.c_n_i["redshift_distribution"]["create_new"]
        zs_inv_cdf = interpolator_from_pickle(
                param_dict_given= dict(z_min=z_min, 
                                       z_max=z_max, 
                                       cosmology=self.cosmo,
                                       event_type=event_type,
                                       merger_rate_density=merger_rate_density,
                                        merger_rate_density_param=merger_rate_density_param,
                                        resolution=resolution,
                                        ),
                directory=directory,
                sub_directory=merger_rate_density,
                name=merger_rate_density,
                x = np.linspace(z_min, z_max, resolution),
                pdf_func= lambda z_: self.pdf_z(zs=z_, param=merger_rate_density_param),
                conditioned_y=None,
                dimension=1,
                category="inv_cdf",
                create_new=create_new,
            )

        # redshift sampler using inverse transform sampling
        self.sample_source_redshift = njit( lambda size, param=None: inverse_transform_sampler(size, zs_inv_cdf[0], zs_inv_cdf[1]) )

    def pdf_z(self, zs, param=None):
        """
        Function to compute the pdf p(z). The output is in source frame and is normalized.

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
        pdf : `numpy.ndarray`
            1D array of floats
            pdf p(z)

        Examples
        ----------
        >>> from ler.gw_source_population import SourceGalaxyPopulationModel
        >>> cbc = SourceGalaxyPopulationModel()
        >>> pdf = cbc.pdf_z(zs=0.1)
        """

        return self.merger_rate_density_src_frame(zs=zs,param=param)/self.normalization_pdf_z

    def merger_rate_density_src_frame(self, zs, param=None):
        """
        Function to compute the merger rate density (source frame). The output is in source frame and is unnormalized.

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
            merger rate density (source frame) (Mpc^-3 yr^-1)

        Examples
        ----------
        >>> from ler.gw_source_population import SourceGalaxyPopulationModel
        >>> cbc = SourceGalaxyPopulationModel()
        >>> rate_density = cbc.merger_rate_density_src_frame(zs=0.1)
        """

        zs = np.array([zs]).reshape(-1)
        # Define the merger-rate density function
        rate_density = (
            self.merger_rate_density(zs, param=param)
            / (1 + zs)
            * (self.differential_comoving_volume(zs))
        )

        return rate_density

    def merger_rate_density_bbh_popI_II_oguri2018(self,
        zs, R0=23.9 * 1e-9, b2=1.6, b3=2.0, b4=30, param=None
    ):
        """
        Function to compute the merger rate density (PopI/PopII). Reference: Oguri et al. (2018). The output is in detector frame and is unnormalized.

        Parameters
        ----------
        zs : `float` or `numpy.ndarray` (nD array of floats)
            Source redshifts
        R0 : `float`
            local merger rate density at low redshift
            default: 23.9*1e-9 Mpc^-3 yr^-1
        b2 : `float`
            Fitting paramters
            default: 1.6
        b3 : `float`
            Fitting paramters
            default: 2.0
        b4 : `float`
            Fitting paramters
            default: 30
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. param = dict(R0=23.9*1e-9, b2=1.6, b3=2.0, b4=30)
            default: None

        Returns
        ----------
        rate_density : `float` or `numpy.ndarray` (nD array of floats)
            merger rate density

        Examples
        ----------
        >>> from ler.gw_source_population import SourceGalaxyPopulationModel
        >>> pop = SourceGalaxyPopulationModel(z_min=0.0, z_max=10, event_type = "BBH", merger_rate_density="merger_rate_density_bbh_popI_II_oguri2018")
        >>> rate_density = pop.merger_rate_density(zs=0.1)
        """

        if param:
            R0 = param["R0"]
            b2 = param["b2"]
            b3 = param["b3"]
            b4 = param["b4"]

        return merger_rate_density_bbh_popI_II_oguri2018(zs=zs, R0=R0, b2=b2, b3=b3, b4=b4)
    
    def star_formation_rate_madau_dickinson2014(self,
        zs, af=2.7, bf=5.6, cf=2.9, param=None
    ):
        """
        Function to compute star formation rate as given in Eqn. 15 Madau & Dickinson (2014). The output is in detector frame and is unnormalized.

        Parameters
        ----------
        zs : `float` or `numpy.ndarray` (nD array of floats)
            Source redshifts
        af : `float`
            Fitting paramters
            default: 2.7
        bf : `float`
            Fitting paramters
            default: 5.6
        cf : `float`
            Fitting paramters
            default: 2.9
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. param = dict(af=2.7, bf=5.6, cf=2.9)
            default: None
            
        Returns
        ----------
        rate_density : `float` or `numpy.ndarray` (nD array of floats)
            merger rate density

        Examples
        ----------
        >>> from ler.gw_source_population import SourceGalaxyPopulationModel
        >>> pop = SourceGalaxyPopulationModel(z_min=5., z_max=40., event_type = "BBH", merger_rate_density="star_formation_rate_madau_dickinson2014")
        >>> rate_density = pop.merger_rate_density(zs=10)
        """

        if param:
            af = param["af"]
            bf = param["bf"]
            cf = param["cf"]

        return star_formation_rate_madau_dickinson2014(zs=zs, af=af, bf=bf, cf=cf)

    def merger_rate_density_bbh_popIII_ken2022(self,
        zs, n0=19.2 * 1e-9, aIII=0.66, bIII=0.3, zIII=11.6, param=None
    ):
        """
        Function to compute the unnormalized merger rate density (PopIII). Reference: Ng et al. 2022. The output is in detector frame and is unnormalized.

        Parameters
        ----------
        zs : `float` or `numpy.ndarray` (nD array of floats)
            Source redshifts
        n0 : `float`
            normalization constant
            default: 19.2*1e-9
        aIII : `float`
            Fitting paramters
            default: 0.66
        bIII : `float`
            Fitting paramters
            default: 0.3
        zIII : `float`
            Fitting paramters
            default: 11.6
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. param = dict(aIII=0.66, bIII=0.3, zIII=11.6)
            default: None

        Returns
        ----------
        rate_density : `float` or `numpy.ndarray` (nD array of floats)
            merger rate density

        Examples
        ----------
        >>> from ler.gw_source_population import SourceGalaxyPopulationModel
        >>> pop = SourceGalaxyPopulationModel(z_min=5, z_max=40, event_type = "BBH", merger_rate_density="merger_rate_density_popIII_ken2022")
        >>> rate_density = pop.merger_rate_density(zs=10)
        >>> rate_density  # Mpc^-3 yr^-1
        1.5107979464621443e-08
        """
        if param:
            n0 = param["n0"]
            aIII = param["aIII"]
            bIII = param["bIII"]
            zIII = param["zIII"]

        return merger_rate_density_bbh_popIII_ken2022(zs=zs, n0=n0, aIII=aIII, bIII=bIII, zIII=zIII)
    
    def merger_rate_density_bbh_primordial_ken2022(self,
        zs, n0=0.044 * 1e-9, t0=13.786885302009708, param=None
    ):
        """
        Function to compute the merger rate density (Primordial). Reference: Ng et al. 2022. The output is in detector frame and is unnormalized.

        Parameters
        ----------
        zs : `float`
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

        if param:
            n0 = param["n0"]
            t0 = param["t0"]

        return merger_rate_density_bbh_primordial_ken2022(zs=zs, cosmology=self.cosmo, n0=n0, t0=t0)
    
    def create_lookup_table(self, z_min, z_max, directory):
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
            param_dict_given= dict(z_min=z_min, z_max=z_max, cosmology=self.cosmo, resolution=resolution),
            directory=directory,
            sub_directory="z_to_luminosity_distance",
            name="z_to_luminosity_distance",
            x = np.linspace(z_min, z_max, resolution),
            pdf_func= lambda z_: self.cosmo.luminosity_distance(z_).value, 
            conditioned_y=None, 
            dimension=1,
            category="function",
            create_new=create_new,
        )
        self.z_to_luminosity_distance = njit(lambda z_: cubic_spline_interpolator(z_, spline1[0], spline1[1]))
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

        # Create a lookup table for the differential comoving volume
        # get differential co-moving volume interpolator
        resolution = self.c_n_i["differential_comoving_volume"]["resolution"]
        create_new = self.c_n_i["differential_comoving_volume"]["create_new"]
        spline2 = interpolator_from_pickle(
            param_dict_given= dict(z_min=z_min, z_max=z_max, cosmology=self.cosmo, resolution=resolution), 
            directory=directory,
            sub_directory="differential_comoving_volume", 
            name="differential_comoving_volume",
            x = np.linspace(z_min, z_max, resolution),
            pdf_func= lambda z_: self.cosmo.differential_comoving_volume(z_).value * 4 * np.pi,  # volume of shell
            conditioned_y=None, 
            dimension=1,
            category="function",
            create_new=create_new,
        )
        self.differential_comoving_volume = njit(lambda z_: cubic_spline_interpolator(z_, spline2[0], spline2[1]))
        self.differential_comoving_volume.__doc__ = """
        Function to calculate the differential comoving volume.

        Parameters
        ----------
        zs : `numpy.ndarray`
            1D array of floats
            Source redshifts

        Returns
        ----------
        differential_comoving_volume : `float`
            1D array of floats
            differential comoving volume
        """

    @property
    def sample_source_redshift(self):
        """
        Function to sample source redshifts (source frame) between z_min and z_max from the source galaxy population

        Parameters
        ----------
        size : `int`
            Number of samples to generate

        Returns
        ----------
        zs : 
            Source redshifts
        """
        return self._sample_source_redshift
    
    @sample_source_redshift.setter
    def sample_source_redshift(self, prior):
        # prior has to be a function
        if callable(prior):
            self._sample_source_redshift = prior
    
    @property
    def merger_rate_density(self):
        """
        Function to get the merger rate density function wrt redshift. The output is in detector frame and is unnormalized.

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
    def merger_rate_density(self, merger_rate_density):
        error_msg = ValueError(f"merger_rate_density must be one of {self.merger_rate_density_model_list}")
        # check if it is a string
        if isinstance(merger_rate_density, str):
            try:
                self._merger_rate_density = getattr(self, merger_rate_density)
            except:
                raise error_msg
        # check if it is a function
        elif callable(merger_rate_density):
            self._merger_rate_density = merger_rate_density
        else:
            raise error_msg
        
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
