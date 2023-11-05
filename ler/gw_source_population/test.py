import warnings

warnings.filterwarnings("ignore")
import numpy as np

# import pycbc
import bilby

# from scipy.stats import randint

# from gwcosmo import priors as p
from scipy.interpolate import interp1d
from scipy.integrate import quad

# for redshift to luminosity distance conversion
from astropy.cosmology import Planck18

# for generating mass distribution
from gwcosmo import priors as p

# for multiprocessing
# Import helper routines
from ..utils import rejection_sample, rejection_sample2d, update_dict


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
        Type of event to generate.
        e.g. 'BBH', 'BNS', 'NSBH'
    merger_rate_density : `str`
        Type of merger rate density function to use
        default: None/'merger_rate_density_popI_II_oguri2018'
        for others see instance method in :class:`~ler.ler.SourceGalaxyPopulationModel`
    merger_rate_density_param : `dict`
        Dictionary of merger rate density function parameters
        default: dict(R0=23.9 * 1e-9, b2=1.6, b3=2.0, b4=30)
    cosmology : `astropy.cosmology`
        Cosmology to use
        default: Planck18

    Examples
    ----------
    >>> from ler.gw_source_population import SourceGalaxyPopulationModel
    >>> cbc = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, merger_rate_density="merger_rate_density_bbh_popI_II_oguri2018")
    >>> zs = cbc.sample_redshift_of_source(size=1000)
    >>> zs[:5]
    array([2.9613628 , 1.18360022, 2.47637065, 2.51401502, 4.22868975])

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
    |:attr:`~merger_rate_density`         | `function`                       |
    +-------------------------------------+----------------------------------+
    |:attr:`~merger_rate_density_param`   | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~normalization_pdf_z`         | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`~z_to_luminosity_distance`    | `scipy.interpolate.interpolate`  |
    +-------------------------------------+----------------------------------+
    |:attr:`~differential_comoving_volume`| `scipy.interpolate.interpolate`  |
    +-------------------------------------+----------------------------------+
    |:attr:`~merger_rate_density_model_list`                                |
    +-------------------------------------+----------------------------------+
    |                                     | Function to list available       |
    |                                     | merger rate density functions    |
    |                                     | and its parameters               |
    +-------------------------------------+----------------------------------+
    |:attr:`~sample_redshift_of_source`   | Function to sample source        |
    |                                     | redshifts (source frame)         |
    +-------------------------------------+----------------------------------+

    Instance Methods
    ----------
    SourceGalaxyPopulationModel has the following instance methods:\n
    +-------------------------------------+----------------------------------+
    | Methods                             | Type                             |
    +=====================================+==================================+
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
    |:meth:`~merger_rate_density_primordial_ken2022`                         |
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
    Dictionary of merger rate density function parameters
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
        z_min=0.0,
        z_max=10.0,
        event_type="BBH",
        merger_rate_density="merger_rate_density_bbh_popI_II_oguri2018",
        merger_rate_density_param=dict(R0=23.9 * 1e-9, b2=1.6, b3=2.0, b4=30),
        cosmology=None,
    ):
        # set attributes
        self.z_min = z_min
        self.z_max = z_max
        self.event_type = event_type
        self.cosmo = cosmology if cosmology else Planck18
        self.create_lookup_table(z_min, z_max)

        # Define the merger-rate density function/method instances
        self.merger_rate_density = merger_rate_density  # function initialization
        self.merger_rate_density_param = merger_rate_density_param  # dict of parameters corresponding to the merger_rate_density function 

        # To find the normalization constant of the pdf p(z)
        self.normalization_pdf_z = quad(
            self.merger_rate_density_src_frame,
            z_min,
            z_max,
            args=(merger_rate_density_param,),
        )[0]

        # Inverse transform sampling
        # create sampler using the pdf p(z)
        # Dont be fooled by the '=' sign. self.pdf_z generates probability density function but self.sample_redshift_of_source is a sampler.
        self.sample_redshift_of_source = self.pdf_z

        return None
    
    @property
    def merger_rate_density(self):
        """
        Function to get the merger rate density function wrt redshift.
        """
        return self._merger_rate_density
    
    @merger_rate_density.setter
    def merger_rate_density(self, merger_rate_density):
        error_msg = ValueError(f"merger_rate_density must be one of {self.merger_rate_density_model_list}")
        # check if it is a string
        if isinstance(merger_rate_density, str):
            try:
                self.merger_rate_density = getattr(self, merger_rate_density)
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
            merger_rate_density_primordial_ken2022=dict(
                n0=0.044 * 1e-9, t0=13.786885302009708
            ),
        )

        return self._merger_rate_density_model_list
    
    @property
    def sample_redshift_of_source(self):
        """
        Function to sample source redshifts (source frame) between z_min and z_max from the source galaxy population
        """
        return self._sample_redshift_of_source
    
    @sample_redshift_of_source.setter
    def sample_redshift_of_source(self, pdf_z):
        # Inverse transform sampling
        if self.z_min==0.:
            z = [0.]+np.geomspace(0.001, self.z_max, 99).tolist()
        else:
            z = np.geomspace(self.z_min, self.z_max, 100).tolist()

        cdf_arr = []
        for i in z:
            cdf_arr.append(quad(pdf_z,
                        self.z_min,
                        i,
                    )[0])
            
        inv_cdf = interp1d(cdf_arr, z, kind='cubic')
        self._sample_redshift_of_source = lambda size: inv_cdf(np.random.uniform(0, 1, size=size))

    def pdf_z(self, zs, param=None):
        """
        Function to compute the pdf p(z). The output is in source frame and is normalized.

        Parameters
        ----------
        zs : `float`
            Source redshifts
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. if the merger_rate_density is merger_rate_density_bbh_popI_II_oguri2018
            param = dict(R0=23.9*1e-9, b2=1.6, b3=2.0, b4=30)

        Returns
        ----------
        pdf : `float`
            pdf p(z)
        """

        return self.merger_rate_density_src_frame(zs=zs,param=param)/self.normalization_pdf_z

    def merger_rate_density_src_frame(self, zs, param=None):
        """
        Function to compute the merger rate density (source frame). The output is in source frame and is unnormalized.

        Parameters
        ----------
        zs : `float`
            Source redshifts
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. if the merger_rate_density is merger_rate_density_bbh_popI_II_oguri2018
            param = dict(R0=23.9*1e-9, b2=1.6, b3=2.0, b4=30)

        Returns
        ----------
        rate_density : `float`
            merger rate density
        """

        # Define the merger-rate density function
        rate_density = (
            self.merger_rate_density(zs, param=param)
            / (1 + zs)
            * self.differential_comoving_volume(zs)
        )

        return rate_density

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
        luminosity_distance = self.cosmo.luminosity_distance(
            z
        ).value  # luminosity distance in Mpc
        self.z_to_luminosity_distance = interp1d(z, luminosity_distance, kind="cubic")

        # Create a lookup table for the differential comoving volume
        dVcdz = self.cosmo.differential_comoving_volume(z).value * 4 * np.pi
        self.differential_comoving_volume = interp1d(
            z, dVcdz, kind="linear", fill_value="extrapolate"
        )

        return None

    def merger_rate_density_bbh_popI_II_oguri2018(
        self, zs, R0=23.9 * 1e-9, b2=1.6, b3=2.0, b4=30, param=None
    ):
        """
        Function to compute the merger rate density (PopI/PopII). Reference: Oguri et al. (2018). The output is in detector frame and is unnormalized.

        Parameters
        ----------
        zs : `float`
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
        rate_density : `float`
            merger rate density

        Examples
        ----------
        >>> from ler.gw_source_population import SourceGalaxyPopulationModel
        >>> cbc = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, merger_rate_density="merger_rate_density_bbh_popI_II_oguri2018")
        >>> rate_density = cbc.merger_rate_density(zs=0.0001) # local merger rate density at low redshift
        >>> rate_density  # Mpc^-3 yr^-1
        2.3903670073287287e-08
        """

        if self.event_type == "BNS":
            R0 = 170.0 * 1e-9
        if self.event_type == "NSBH":
            R0 = 27.0 * 1e-9

        if param:
            R0 = param["R0"]
            b2 = param["b2"]
            b3 = param["b3"]
            b4 = param["b4"]

        # rate_density = R0 * (b4 + 1) * np.exp(b2 * zs) / (b4 + np.exp(b3 * zs))
        rate_density = R0 * (b4 + 1) * np.exp(b2 * zs) / (b4 + np.exp(b3 * zs))

        return rate_density

    def star_formation_rate_madau_dickinson2014(
        self, zs, af=2.7, bf=5.6, cf=2.9, param=None
    ):
        """
        Function to compute star formation rate as given in Eqn. 15 Madau & Dickinson (2014).

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
            default: 2.9
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. param = dict(af=2.7, bf=5.6, cf=2.9)
            default: None

        Returns
        ----------
        rate_density : `float`
            merger rate density

        Examples
        ----------
        >>> from ler.gw_source_population import SourceGalaxyPopulationModel
        >>> cbc = SourceGalaxyPopulationModel(z_min=0.0001, z_max=10, merger_rate_density="star_formation_rate_madau_dickinson2014")
        >>> rate_density = cbc.merger_rate_density(zs=0.0001) # local merger rate density at low redshift
        >>> rate_density  # Mpc^-3 yr^-1
        0.014965510855362926
        """
        if param:
            af = param["af"]
            bf = param["bf"]
            cf = param["cf"]

        # rate density
        rate_density = 0.015 * (1 + zs) ** af / (1 + ((1 + zs) / cf) ** bf)

        return rate_density

    def merger_rate_density_popIII_ken2022(
        self, zs, n0=19.2 * 1e-9, aIII=0.66, bIII=0.3, zIII=11.6, param=None
    ):
        """
        Function to compute the unnormalized merger rate density (PopIII). Reference: Ng et al. 2022. The output is in detector frame and is unnormalized.

        Parameters
        ----------
        zs : `float`
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
        rate_density : `float`
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
            aIII = param["aIII"]
            bIII = param["bIII"]
            zIII = param["zIII"]

        # rate density
        rate_density = (
            n0
            * np.exp(aIII * (zs - zIII))
            / (bIII + aIII * np.exp((aIII + bIII) * (zs - zIII)))
        )

        return rate_density

    def merger_rate_density_primordial_ken2022(
        self, zs, n0=0.044 * 1e-9, t0=13.786885302009708, param=None
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
        >>> pop = SourceGalaxyPopulationModel(z_min=5, z_max=40, event_type = "BBH", merger_rate_density="merger_rate_density_primordial_ken2022")
        >>> rate_density = pop.merger_rate_density(zs=10)
        >>> rate_density  # Mpc^-3 yr^-1
        9.78691173794454e-10
        """
        if param:
            t0 = param["t0"]

        # rate density
        rate_density = n0 * (self.cosmo.age(z=zs).value / t0) ** (-34 / 37)

        return rate_density