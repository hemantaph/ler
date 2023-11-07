# -*- coding: utf-8 -*-
"""
This module contains the classes to generate source galaxy population model
and compact binary population model. The compact binary population model
inherits from the source galaxy population model. The source galaxy population
model is used to generate the source galaxy population model. The compact binary
population model is used to generate the compact binary population model.

"""

import numpy as np
import bilby
from scipy.stats import randint

# from gwcosmo import priors as p
from scipy.interpolate import interp1d
from scipy.integrate import quad

# for redshift to luminosity distance conversion
from astropy.cosmology import Planck18

# for generating mass distribution
from gwcosmo import priors as p

# for multiprocessing
# Import helper routines
from ler.helperroutines import rejection_sample, rejection_sample2d


class SourceGalaxyPopulationModel:
    """Class to generate a population of source galaxies.
    This class is inherited by :class:`~ler.ler.CompactBinaryPopulation` class.

    Parameters
    ----------
    z_min : `float`
        Minimum redshift of the source population
        default: 0.
    z_max : `float`
        Maximum redshift of the source population
        default: 10.
    event_type : `str`
        Type of event to generate
        e.g. 'popI_II', 'BNS', 'popIII', 'primordial', 'popI_II_Madau_Dickinson'
        default: 'popI_II'

    Examples
    ----------
    >>> from ler import SourceGalaxyPopulationModel
    >>> pop = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, event_type = "popI_II")
    >>> zs = pop.sample_source_redshifts(size=1000)
    >>> zs
    array([0.0001, 0.0001, 0.0001, ..., 9.9999, 9.9999, 9.9999])

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
    |:attr:`~normalization_pdf_z`         | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`~z_to_luminosity_distance`    | `scipy.interpolate.interpolate`  |
    +-------------------------------------+----------------------------------+
    |:attr:`~differential_comoving_volume`| `scipy.interpolate.interpolate`  |
    +-------------------------------------+----------------------------------+

    Instance Methods
    ----------
    SourceGalaxyPopulationModel has the following instance methods:\n
    +-------------------------------------+----------------------------------+
    | Methods                             | Type                             |
    +=====================================+==================================+
    |:meth:`~create_lookup_table`         | Function to create a lookup      |
    |                                     | table for the differential       |
    |                                     | comoving volume and luminosity   |
    |                                     | distance wrt redshift            |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_source_redshifts`     | Function to sample source        |
    |                                     | redshifts from the source        |
    |                                     | galaxy population model          |
    +-------------------------------------+----------------------------------+
    |:meth:`~merger_rate_density_popI_II` | Function to compute the merger   |
    |                                     | rate density (PopI/PopII)        |
    +-------------------------------------+----------------------------------+
    |:meth:`~merger_rate_density_popI_II_Madau_Dickinson`                    |
    +-------------------------------------+----------------------------------+
    |                                     | Function to compute the          |
    |                                     | merger rate density (PopI/PopII) |
    |                                     | from Madau & Dickinson (2014)    |
    +-------------------------------------+----------------------------------+
    |:meth:`~merger_rate_density_popIII`  | Function to compute the merger   |
    |                                     | rate density (PopIII)            |
    +-------------------------------------+----------------------------------+
    |:meth:`~merger_rate_density_primordial`                                 |
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
            z_min=0.0001,
            z_max=10.0,
            event_type="popI_II",
            category="Oguri",
            merger_rate_density_fn=None,
            merger_rate_density_param=None):

        # set attributes
        self.z_min = z_min
        self.z_max = z_max
        self.create_lookup_table(z_min, z_max)
        self.merger_rate_density_param = merger_rate_density_param

        # To find the normalization constant of the pdf p(z)
        # Define the merger-rate density function
        if merger_rate_density_fn is not None:
            self.merger_rate_density = merger_rate_density_fn
        else:
            try:
                if category == None:
                    category = ""
                else:
                    category = "_" + category
                self.merger_rate_density = getattr(
                    self, "merger_rate_density_" + event_type + category
                )
            except:
                self.merger_rate_density = getattr(self, "merger_rate_density_" + "popI_II_Oguri")

        merger_rate_density_source_frame = lambda z: self.merger_rate_density(z) / (
            1 + z
        )

        # Define the pdf p(z)
        pdf_unnormalized = lambda z: merger_rate_density_source_frame(
            z
        ) * self.differential_comoving_volume(z)
        # Normalize the pdf
        # this normalization factor is common no matter what you choose for z_min and z_max
        self.normalization_pdf_z = quad(pdf_unnormalized, z_min, z_max)[0]
        return None

    def create_lookup_table(self, z_min, z_max):
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
        z = np.linspace(z_min, z_max, 500)  # redshift
        luminosity_distance = Planck18.luminosity_distance(
            z
        ).value  # luminosity distance in Mpc
        self.z_to_luminosity_distance = interp1d(z, luminosity_distance, kind="cubic")

        # Create a lookup table for the differential comoving volume
        dVcdz = Planck18.differential_comoving_volume(z).value * 4 * np.pi
        self.differential_comoving_volume = interp1d(
            z, dVcdz, kind="linear", fill_value="extrapolate"
        )

        return None

    def sample_source_redshifts(self, size=1000, z_min=0.0, z_max=10.0):
        """
        Function to sample source redshifts from the source galaxy population
        model

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        z_min : `float`
            Minimum redshift of the source population
        z_max : `float`
            Maximum redshift of the source population

        Returns
        ----------
        zs : `array`
            Array of sampled redshifts

        Examples
        ----------
        >>> from ler import SourceGalaxyPopulationModel
        >>> pop = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, event_type = "popI_II_Oguri")
        >>> zs = pop.sample_source_redshifts(size=1000)
        >>> zs
        array([0.0001, 0.0001, 0.0001, ..., 9.9999, 9.9999, 9.9999])

        """
        # Define the merger-rate density function
        # print(f"z_max = {z_max}")
        merger_rate_density_source_frame = lambda z: self.merger_rate_density(z) / (
            1 + z
        )
        # Define the pdf p(z)
        pdf_unnormalized = lambda z: merger_rate_density_source_frame(
            z
        ) * self.differential_comoving_volume(z)
        # Normalize the pdf
        normalization = self.normalization_pdf_z
        pdf = lambda z: pdf_unnormalized(z) / normalization
        # Sample the redshifts using rejection sampling
        zs = rejection_sample(pdf, z_min, z_max, size=size)
        return zs

    def merger_rate_density_popI_II_Oguri(self, zs, R0=25.0 * 1e-9, b2=1.6, b3=2.0, b4=30):
        """
        Function to compute the merger rate density (PopI/PopII)

        Parameters
        ----------
        zs : `float`
            Source redshifts
        R0 : `float`
            Normalization constant
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

        Returns
        ----------
        rate_density : `float`
            merger rate density

        Examples
        ----------
        >>> from ler import SourceGalaxyPopulationModel
        >>> pop = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, event_type = "popI_II_Oguri")
        >>> rate_density = pop.merger_rate_density_popI_II_Oguri(zs=0.1)
        >>> rate_density
        2.7848018586883885e-08

        """
        # if self.event_type == "BNS":
        #     R0 = 170.0 * 1e-9
        # if self.event_type == "NSBH":
        #     R0 = 27.0 * 1e-9

        #print("\n merger_rate_density_popI_II_Oguri \n")
        # replace values with self.merger_rate_density_param, if given
        param = dict(R0=R0, b2=b2, b3=b3, b4=b4)
        if self.merger_rate_density_param is not None:
            keys_ = param.keys()
            for key, value in self.merger_rate_density_param.items():
                if key in keys_:
                    param[key] = value
                
        # rate_density = R0 * (b4 + 1) * np.exp(b2 * zs) / (b4 + np.exp(b3 * zs))
        rate_density = param["R0"] * (param["b4"] + 1) * np.exp(
            param["b2"] * zs
        ) / (param["b4"] + np.exp(param["b3"] * zs))


        return rate_density

    def merger_rate_density_popI_II_Madau_Dickinson(self, zs, af=2.7, bf=5.6, cf=1.9):
        """
        Function to compute the unormalized merger rate density (PopI/PopII) from Madau & Dickinson (2014)

        Parameters
        ----------
        zs : `float`
            Source redshifts
        af : `float`
            Fitting paramters
            default: 2.7
        bf : `float`
            Fitting paramters
            default: 5.6
        cf : `float`
            Fitting paramters
            default: 1.9

        Returns
        ----------
        rate_density : `float`
            merger rate density

        Examples
        ----------
        >>> from ler import SourceGalaxyPopulationModel
        >>> pop = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, event_type = "popI_II_Madau_Dickinson")
        >>> rate_density = pop.merger_rate_density_popI_II_Madau_Dickinson(zs=0.1)
        >>> rate_density
        1.2355851838964846

        """
        # print("\n merger_rate_density_popI_II_Madau_Dickinson \n")
        # replace values with self.merger_rate_density_param, if given
        param = dict(af=af, bf=bf, cf=cf)
        if self.merger_rate_density_param is not None:
            keys_ = param.keys()
            for key, value in self.merger_rate_density_param.items():
                if key in keys_:
                    param[key] = value

        rate_density = 0.015 * (1 + zs) ** param["af"] / (1 + ((1 + zs) / param["cf"]) ** param["bf"])

        return rate_density

    def merger_rate_density_popIII_Ken(self, zs, aIII=0.66, bIII=0.3, zIII=11.6):
        """
        Function to compute the unnormalized merger rate density (PopIII)

        Parameters
        ----------
        zs : `float`
            Source redshifts
        aIII : `float`
            Fitting paramters
            default: 0.66
        bIII : `float`
            Fitting paramters
            default: 0.3
        zIII : `float`
            Fitting paramters
            default: 11.6

        Returns
        ----------
        rate_density : `float`
            merger rate density

        Examples
        ----------
        >>> from ler import SourceGalaxyPopulationModel
        >>> pop = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, event_type = "popIII")
        >>> rate_density = pop.merger_rate_density_popIII(zs=0.1)
        >>> rate_density
        0.00010000000000000002
        """
        #print("\n ken \n")
        # replace values with self.merger_rate_density_param, if given
        param = dict(aIII=aIII, bIII=bIII, zIII=zIII)
        if self.merger_rate_density_param is not None:
            keys_ = param.keys()
            for key, value in self.merger_rate_density_param.items():
                if key in keys_:
                    param[key] = value

        rate_density = np.exp(param["aIII"] * (zs - param["zIII"])) / ( param["bIII"] + param["aIII"] * np.exp((param["aIII"] + param["bIII"]) * (zs - param["zIII"]))
        )

        return rate_density

    def merger_rate_density_primordial_Ken(self, zs, t0=13.786885302009708):
        """
        Function to compute the merger rate density (Primordial)

        Parameters
        ----------
        zs : `float`
            Source redshifts
        t0 : `float`
            Present ge of the Universe in Gyr
            default: 13.786885302009708

        Returns
        ----------
        rate_density : `float`
            merger rate density

        Examples
        ----------
        >>> from ler import SourceGalaxyPopulationModel
        >>> pop = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, event_type = "primordial")
        >>> rate_density = pop.merger_rate_density_primordial(zs=0.1)
        >>> rate_density
        0.00010000000000000002
        """
        # replace values with self.merger_rate_density_param, if given
        param = dict(t0=t0)
        if self.merger_rate_density_param is not None:
            keys_ = param.keys()
            for key, value in self.merger_rate_density_param.items():
                if key in keys_:
                    param[key] = value

        rate_density = (Planck18.age(z=zs).value / t0) ** (-34 / 37)

        return rate_density


class CompactBinaryPopulation(SourceGalaxyPopulationModel):
    """Class to generate a population of compact binaries. Inherits from :class:`~ler.ler.SourceGalaxyPopulationModel` class.

    Parameters
    ----------
    z_min : `float`
        Minimum redshift of the source population
    z_max : `float`
        Maximum redshift of the source population
    m_min : `float`
        Minimum mass of the BBHs
    m_max : `float`
        Maximum mass of the BBHs
    event_type : `str`
        Type of event to generate.
        e.g. 'popI_II', 'BNS', 'popIII', 'primordial', 'popI_II_Madau_Dickinson'
    src_model_params : `dict`
        Dictionary of model parameters.
        e.g. for popI_II: {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, 'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, 'mu_g': 33.07, 'sigma_g': 5.69}

    Examples
    ----------
    >>> from ler import CompactBinaryPopulation
    >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popI_II")
    >>> gw_parameters = pop.sample_gw_parameters(nsamples=1000)
    >>> gw_parameters.keys()
    dict_keys(['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl'])

    Instance Attributes
    ----------
    CompactBinaryPopulation has the following instance attributes:\n
    +-------------------------------------+----------------------------------+
    | Atrributes                          | Type                             |
    +=====================================+==================================+
    |:attr:`~z_min`                       | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`~z_max`                       | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`~m_min`                       | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`~m_max`                       | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`~event_type`                  | `str`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~src_model_params`                  | `dict`                           |
    +-------------------------------------+----------------------------------+

    Instance Methods
    ----------
    CompactBinaryPopulation has the following instance methods:\n
    +-------------------------------------+----------------------------------+
    | Methods                             | Type                             |
    +=====================================+==================================+
    |:meth:`~sample_gw_parameters`        | Function for sampling GW         |
    |                                     | parameters from the source       |
    |                                     | galaxy population model          |
    +-------------------------------------+----------------------------------+
    |:meth:`~binary_masses_popI_II`       | Function to calculate source     |
    |                                     | mass1 and mass2 with             |
    |                                     | PowerLaw+PEAK model              |
    +-------------------------------------+----------------------------------+
    |:meth:`~binary_masses_popIII`        | Function to calculate source     |
    |                                     | mass1 and mass2 with pop III     |
    |                                     | origin                           |
    +-------------------------------------+----------------------------------+
    |:meth:`~binary_masses_primordial`    | Function to calculate source     |
    |                                     | mass1 and mass2 for primordial   |
    |                                     | BBHs                             |
    +-------------------------------------+----------------------------------+
    |:meth:`~binary_masses_BNS`           | Function to calculate source     |
    |                                     | mass1 and mass2 of BNS           |
    +-------------------------------------+----------------------------------+
    |:meth:`~mass_ratio`                  | Function to calculate mass ratio |
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

    m_min = None
    """``float`` \n
    Minimum mass of the BBHs
    """

    m_max = None
    """``float`` \n
    Maximum mass of the BBHs
    """

    event_type = None
    """``str`` \n
    Type of event to generate. \n
    e.g. 'popI_II', 'BNS', 'popIII', 'primordial', 'popI_II_Madau_Dickinson'
    """

    src_model_params = None
    """``dict`` \n
    Dictionary of model parameters. \n
    e.g. for popI_II: {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, 'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, 'mu_g': 33.07, 'sigma_g': 5.69} \n
    for popI_II_Madau_Dickinson: {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, 'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, 'mu_g': 33.07, 'sigma_g': 5.69} \n
    for popIII: None \n
    for primordial: {'Mc':30.,'sigma':0.3,'beta':1.1} \n
    for BNS: None

    Examples
    ----------
    >>> from ler import CompactBinaryPopulation
    >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popI_II")
    >>> method_list = [method for method in dir(pop) if method.startswith('__') is False]
    >>> print(method_list)
    ['create_lookup_table', 'differential_comoving_volume', 'merger_rate_density', 'merger_rate_density_popIII', 'merger_rate_density_popI_II', 'merger_rate_density_popI_II_Madau_Dickinson', 'merger_rate_density_primordial', 'normalization_pdf_z', 'sample_source_redshifts', 'z_max', 'z_min', 'z_to_luminosity_distance']

    """

    def __init__(
        self,
        z_min=0.0001,
        z_max=10,
        m_min=None,
        m_max=None,
        event_type="BBH",
        category=None,
        sub_category=None,
        redshift_event_type=None,
        redshift_category=None,
        merger_rate_density_fn=None,
        merger_rate_density_param=None,
        src_model_params=None,
        redshift_constant = False,
        mass_constant=False,
        spin_constant=0.,
    ):
        # initialized SourceGalaxyPopulationModel mother class
        # self.z_min and self.z_max will inherit from SourceGalaxyPopulationModel
        if redshift_event_type is None:
            redshift_event_type = category

        super().__init__(z_min=z_min, z_max=z_max, event_type=redshift_event_type, category=redshift_category, merger_rate_density_fn=merger_rate_density_fn, merger_rate_density_param=merger_rate_density_param)

        # check redshift
        self.redshift_constant = redshift_constant
        if redshift_constant:
            self.sample_source_redshifts = None

        # checking for event-type, category, sub-categoru
        # in BBH_popI_II_gwcosmo, "BBH": event type, "popI_II": category, "gwcosmo": sub-category
        if category==None:
            category = "_popI_II"
        else:
            category = "_"+category
        if sub_category==None:
            sub_category = "_gwcosmo"
        else:
            sub_category = "_"+sub_category

        # check mass
        self.mass_constant = mass_constant
        if mass_constant:
            self.source_binary_masses = None
        else:
            self.m_min = m_min
            self.m_max = m_max
            try:
                self.source_binary_masses = getattr(self, "binary_masses_" + event_type + category + sub_category)
            except:
                pass

        # check spin
        if spin_constant or spin_constant==0.:
            self.spin_constant = spin_constant
        else:
            try:
                self.binary_spin = getattr(self, "binary_spin_" + event_type)
            except:
                self.binary_spin = getattr(self, "binary_spin_BBH")

        # set attribute
        self.src_model_params = src_model_params

        return None

    def sample_gw_parameters(self, nsamples=1000, verbose=False, **kwargs):
        """
        Function to sample BBH parameters from the source galaxy population
        model

        Parameters
        ----------
        nsamples : `int`
            Number of samples to draw
        kwargs : `dict`
            Keyword arguments to pass in parameter values
            e.g. zs = np.array([0.1,0.2,0.3])

        Returns
        ----------
        gw_parameters : `dict`
            Dictionary of sampled parameters
            gw_parameters.keys() = ['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl']

        Examples
        ----------
        >>> from ler import CompactBinaryPopulation
        >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popI_II")
        >>> gw_parameters = pop.sample_gw_parameters(nsamples=1000)
        >>> gw_parameters.keys()
        dict_keys(['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl'])

        """
        if verbose:
            print(
                f"refer to https://arxiv.org/abs/1908.06068 for parameter definitions, priors and sampling methods"
            )

        def warning(param_):
            try:
                len(param_) == nsamples
            except:
                print(f"Error: len({param_}) != nsamples")
                return 0

        # sampling redshifts and luminosity distances
        if self.redshift_constant:
            if type(self.redshift_constant)==float or type(self.redshift_constant)==int:
                zs = self.redshift_constant*np.ones(nsamples)
            else:
                # warning
                print("Error: redshift_constant should be a float or false")
        else:
            try:
                zs = kwargs["zs"]
            except:
                zs = self.sample_source_redshifts(size=nsamples, z_min=self.z_min, z_max=self.z_max)
        luminosity_distance = self.z_to_luminosity_distance(zs)  # Mpc

        # sampling mass1 and mass2
        if self.mass_constant:
            if type(self.mass_constant)==list:
                mass_1_source = self.mass_constant[0]*np.ones(nsamples)
                mass_2_source = self.mass_constant[1]*np.ones(nsamples)
            elif type(self.mass_constant)==float or type(self.mass_constant)==int:
                mass_1_source = self.mass_constant*np.ones(nsamples)
                mass_2_source = self.mass_constant*np.ones(nsamples)
            else:
                # warning
                print("Error: mass_constant should be a float/list or false")
        else:
            try:
                kwargs["mass_1"]
                kwargs["mass_2"]
            except:
                mass_1_source, mass_2_source = self.source_binary_masses(
                    size=nsamples)
            else:
                mass_1 = kwargs["mass_1"]
                warning(mass_1)
                mass_2 = kwargs["mass_2"]
                warning(mass_1)
                mass_1_source, mass_2_source = self.source_binary_masses(
                    size=nsamples)

        # Sample all other parameters
        # use bilby priors
        bilby.core.utils.logger.disabled = True
        prior_default = bilby.gw.prior.BBHPriorDict()
        # draw associated angles
        ra = prior_default["ra"].sample(nsamples)
        dec = prior_default["dec"].sample(nsamples)
        psi = prior_default["psi"].sample(nsamples)
        theta_jn = prior_default["theta_jn"].sample(nsamples)
        phase = prior_default["phase"].sample(nsamples)
        # sampling spin parameters
        if self.spin_constant==0. or self.spin_constant==True:
            a_1 = np.zeros(nsamples)
            a_2 = np.zeros(nsamples)
            tilt_1 = np.zeros(nsamples)
            tilt_2 = np.zeros(nsamples)
            phi_12 = np.zeros(nsamples)
            phi_jl = np.zeros(nsamples)
        elif type(self.spin_constant)==list:
            a_1 = self.spin_constant[0]*np.ones(nsamples)
            a_2 = self.spin_constant[1]*np.ones(nsamples)
            tilt_1 = self.spin_constant[2]*np.ones(nsamples)
            tilt_2 = self.spin_constant[3]*np.ones(nsamples)
            phi_12 = self.spin_constant[4]*np.ones(nsamples)
            phi_jl = self.spin_constant[5]*np.ones(nsamples)
        else:
            a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl = self.binary_spin(nsamples)

        # compute GPS time
        geocent_time = randint.rvs(1238166018, 1238166018 + 31536000, size=nsamples)
        mass_1, mass_2 = mass_1_source * (1 + zs), mass_2_source * (1 + zs)

        gw_parameters = {
            "mass_1": mass_1,
            "mass_2": mass_2,
            "mass_1_source": mass_1_source,
            "mass_2_source": mass_2_source,
            "zs": zs,
            "luminosity_distance": luminosity_distance,
            "iota": theta_jn,
            "psi": psi,
            "phase": phase,
            "geocent_time": geocent_time,
            "ra": ra,
            "dec": dec,
            "a_1": a_1,
            "a_2": a_2,
            "tilt_1": tilt_1,
            "tilt_2": tilt_2,
            "phi_12": phi_12,
            "phi_jl": phi_jl,
        }

        return gw_parameters

    def binary_masses_BBH_popI_II_gwcosmo(self, size, alpha=3.63, beta=1.26, delta_m=4.82, mmin=4.59, mmax=86.22, lambda_peak=0.08, mu_g=33.07, sigma_g=5.69):
        """
        Function to calculate source mass1 and mass2 with PowerLaw+PEAK model

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        src_model_params : `dict`
            Dictionary of model parameters
            e.g. {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, 'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, 'mu_g': 33.07, 'sigma_g': 5.69}

        Returns
        ----------
        mass_1_source : `array`
            Array of mass1 in source frame
        mass_2_source : `array`
            Array of mass2 in source frame

        Examples
        ----------
        >>> from ler import CompactBinaryPopulation
        >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popI_II")
        >>> src_model_params = {'alpha': 3.63, 'beta': 1.26, 'delta_m': 4.82, 'mmin': 4.59, 'mmax': 86.22, 'lambda_peak': 0.08, 'mu_g': 33.07, 'sigma_g': 5.69}
        >>> mass_1_source, mass_2_source = pop.binary_masses_popI_II(size=1000, src_model_params=src_model_params)

        """
        # replace values with self.src_model_params, if given
        param = dict(alpha=alpha, beta=beta, delta_m=delta_m, mmin=mmin, mmax=mmax, lambda_peak=lambda_peak, mu_g=mu_g, sigma_g=sigma_g)
        if self.src_model_params is not None:
            keys_ = param.keys()
            for key, value in self.src_model_params.items():
                if key in keys_:
                    param[key] = value

        if self.m_min is not None:
            param["mmin"] = self.m_min
        if self.m_max is not None:
            param["mmax"] = self.m_max

        model = p.mass_prior("BBH-powerlaw-gaussian", param)
        mass_1_source, mass_2_source = model.sample(Nsample=size)
        while np.any(mass_2_source > mass_1_source):
            mass_1_source, mass_2_source = model.sample(Nsample=size)

        return (mass_1_source, mass_2_source)

    def binary_masses_BBH_popIII_gwcosmo(self, size, Mc=30.0, sigma=0.3, beta=1.1):
        """
        Function to calculate source mass1 and mass2 with pop III origin

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        Mc, sigma, beta : `float`
            Fitting parameters
            default: Mc=30.0, sigma=0.3, beta=1.1

        Returns
        ----------
        mass_1_source : `array`
            Array of mass1 in source frame
        mass_2_source : `array`
            Array of mass2 in source frame

        Examples
        ----------
        >>> from ler import CompactBinaryPopulation
        >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popIII")
        >>> mass_1_source, mass_2_source = pop.binary_masses_popIII(size=1000)

        """
        if self.m_min is not None:
            m_min = self.m_min
        else:
            m_min = 4.59
        if self.m_max is not None:
            m_max = self.m_max
        else:
            m_max = 86.22

        # replace values with self.src_model_params, if given
        param = dict(Mc=Mc, sigma=sigma, beta=beta)
        if self.src_model_params is not None:
            keys_ = param.keys()
            for key, value in self.src_model_params.items():
                if key in keys_:
                    param[key] = value

        
        psi = lambda m: np.exp(-np.log(m / param["Mc"]) ** 2 / (2 * param["sigma"]**2)) / (
            np.sqrt(2 * np.pi) * param["sigma"] * m
        )
    
        pdf = lambda m1,m2: (m1+m2)**(36/37) * (m1*m2)**(32/37) * psi(m1) * psi(m2)
        
        mass_1_source, mass_2_source = rejection_sample2d(pdf=pdf, xmin=m_min, xmax=m_max, ymin=m_min, ymax=m_max, size=size)
        
        idx = mass_2_source>mass_1_source
        mass_1_source[idx], mass_2_source[idx] = mass_2_source[idx], mass_1_source[idx]
        
        return (mass_1_source, mass_2_source)

    def binary_masses_BBH_primordial_lognormal(self, size, Mc=30.0, sigma=0.3, beta=1.1):
        """
        Function to calculate source mass1 and mass2 for primordial BBHs

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        Mc, sigma, beta : `float`
            Fitting parameters
            default: Mc=30.0, sigma=0.3, beta=1.1

        Returns
        ----------
        mass_1_source : `array`
            Array of mass1 in source frame
        mass_2_source : `array`
            Array of mass2 in source frame

        Examples
        ----------
        >>> from ler import CompactBinaryPopulation
        >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "primordial")
        >>> mass_1_source, mass_2_source = pop.binary_masses_primordial(size=1000)

        """
        if self.m_min is not None:
            m_min = self.m_min
        else:
            m_min = 4.59
        if self.m_max is not None:
            m_max = self.m_max
        else:
            m_max = 86.22

        # replace values with self.src_model_params, if given
        param = dict(Mc=Mc, sigma=sigma, beta=beta)
        if self.src_model_params is not None:
            keys_ = param.keys()
            for key, value in self.src_model_params.items():
                if key in keys_:
                    param[key] = value

        
        psi = lambda m: np.exp(-np.log(m / param["Mc"]) ** 2 / (2 * param["sigma"]**2)) / (
            np.sqrt(2 * np.pi) * param["sigma"] * m
        )
    
        pdf = lambda m1,m2: (m1+m2)**(36/37) * (m1*m2)**(32/37) * psi(m1) * psi(m2)
        
        mass_1_source, mass_2_source = rejection_sample2d(pdf=pdf, xmin=m_min, xmax=m_max, ymin=m_min, ymax=m_max, size=size)
        
        idx = mass_2_source>mass_1_source
        mass_1_source[idx], mass_2_source[idx] = mass_2_source[idx], mass_1_source[idx]
        
        return (mass_1_source, mass_2_source)

    def binary_masses_BNS_popI_II_gwcosmo(self, size, mminns=1.0, mmaxns=3.0, alphans=0.0):
        """
        Function to calculate source mass1 and mass2 of BNS (gwcosmo)

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        mminns : `float`
            Minimum mass of the BNS
            default: 1.0
        mmaxns : `float`
            Maximum mass of the BNS
            default: 3.0
        alphans : `float`
            Power law index
            default: 0.0

        Returns
        ----------
        mass_1_source : `array`
            Array of mass1 in source frame
        mass_2_source : `array`
            Array of mass2 in source frame
        """

        # replace values with self.src_model_params, if given
        param = dict(mminns=1.0, mmaxns=3.0, alphans=0.0)
        if self.src_model_params is not None:
            keys_ = param.keys()
            for key, value in self.src_model_params.items():
                if key in keys_:
                    param[key] = value

        if self.m_min is not None:
            param["mminns"] = self.m_min
        if self.m_max is not None:
            param["mmaxns"] = self.m_max

        model = p.mass_prior("BNS", param)
        mass_1_source, mass_2_source = model.sample(Nsample=size)
        while np.any(mass_2_source > mass_1_source):
            mass_1_source, mass_2_source = model.sample(Nsample=size)

        return (mass_1_source, mass_2_source)
    
    def binary_masses_BNS_popI_II_Alsing(self, size, w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1., mmax=2.3):
        """
        Function to calculate source mass1 and mass2 of BNS (Alsing)

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        w, muL, sigmaL, muR, sigmaR : `float`
            Fitting parameters
            default: w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3
        mmin : `float`
            Minimum mass of the BNS
            default: 1.0
        mmax : `float`
            Maximum mass of the BNS
            default: 3.0

        Returns
        ----------
        mass_1_source : `array`
            Array of mass1 in source frame
        mass_2_source : `array`
            Array of mass2 in source frame
        """
        # print("\n BNS alsing \n")
        # replace values with self.src_model_params, if given
        param = dict(w=0.643, muL=1.352, sigmaL=0.08, muR=1.88, sigmaR=0.3, mmin=1., mmax=2.3)
        if self.src_model_params is not None:
            keys_ = param.keys()
            for key, value in self.src_model_params.items():
                if key in keys_:
                    param[key] = value

        if self.m_min is not None:
            param["mmin"] = self.m_min
        if self.m_max is not None:
            param["mmax"] = self.m_max

        pdf_unnormL = lambda m: np.exp(-(m-muL)**2 / (2*sigmaL**2))  
        normL = quad(pdf_unnormL, mmin, mmax)[0]
        pdf_unnormR = lambda m: np.exp(-(m-muR)**2 / (2*sigmaR**2))  
        normR = quad(pdf_unnormR, mmin, mmax)[0]

        pdf = lambda m: w*pdf_unnormL(m) / normL + (1-w)*pdf_unnormR(m) / normR

        mass_1_source = rejection_sample(pdf, mmin, mmax, size=size)
        mass_2_source = rejection_sample(pdf, mmin, mmax, size=size)
        # swap mass_1_source and mass_2_source if mass_2_source > mass_1_source
        idx = mass_2_source>mass_1_source
        mass_1_source[idx], mass_2_source[idx] = mass_2_source[idx], mass_1_source[idx] 

        return(mass_1_source, mass_2_source)
        
    def mass_ratio(self, size, beta=1.1):
        """
        Function to calculate mass ratio with power law distribution

        Parameters
        ----------
        size : `int`
            Number of samples
        beta : `float`
            Power law index

        Returns
        ----------
        q : `array`
            Array of mass ratio

        Examples
        ----------
        >>> from ler import CompactBinaryPopulation
        >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=1.0, m_max=3.0, event_type = "BNS")
        >>> q = pop.mass_ratio(size=1000, beta=1.1)

        """
        pdf = lambda q: q**beta
        q = rejection_sample(pdf, 0, 1, size=size)
        return q
    
    def binary_spin_BBH(self, size):
        """
        Function to calculate spin parameters with PowerLaw+PEAK model

        Parameters
        ----------
        size : `int`
            Number of samples to draw

        Returns
        ----------
        a_1 : `array`
            Array of spin1
        a_2 : `array`
            Array of spin2
        tilt_1 : `array`
            Array of tilt1
        tilt_2 : `array`
            Array of tilt2
        phi_12 : `array`
            Array of phi12
        phi_jl : `array`
            Array of phi_jl

        Examples
        ----------
        >>> from ler import CompactBinaryPopulation
        >>> pop = CompactBinaryPopulation(z_min=0.0001, z_max=10, m_min=4.59, m_max=86.22, event_type = "popI_II")
        >>> a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl = pop.binary_spin_BBH(size=1000)

        """
        bilby.core.utils.logger.disabled = True
        prior_default = bilby.gw.prior.BBHPriorDict()
        a_1 = prior_default["a_1"].sample(size)
        a_2 = prior_default["a_2"].sample(size)
        tilt_1 = prior_default["tilt_1"].sample(size)
        tilt_2 = prior_default["tilt_2"].sample(size)
        phi_12 = prior_default["phi_12"].sample(size)
        phi_jl = prior_default["phi_jl"].sample(size)

        return (a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl)
    
