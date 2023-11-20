# -*- coding: utf-8 -*-
"""
This module contains the LensGalaxyPopulation class, which is used to sample lens galaxy parameters, source parameters conditioned on the source being strongly lensed. \n
The class inherits from the ImageProperties class, which is used calculate image properties (magnification, timedelays, source position, image position, morse phase). \n
Either the class takes in initialized CompactBinaryPopulation class as input or inherits the CompactBinaryPopulation class with default params (if no input) \n
"""

import warnings

warnings.filterwarnings("ignore")
import numpy as np
from scipy.stats import gengamma, rayleigh, norm
from scipy.interpolate import interp1d
from scipy.integrate import quad
from lenstronomy.Util.param_util import phi_q2_ellipticity

# for redshift to luminosity distance conversion
from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
from astropy import constants as const

# the following .py file will be called if they are not given in the class initialization
from ..gw_source_population import CompactBinaryPopulation
from ..image_properties import ImageProperties
from ..utils import add_dictionaries_together, trim_dictionary


class LensGalaxyParameterDistribution(CompactBinaryPopulation, ImageProperties):
    """
    Class to sample lens galaxy parameters

    Parameters
    ----------
    CompactBinaryPopulation_ : CompactBinaryPopulation class
        This is an already initialized class that contains a function (CompactBinaryPopulation.sample_gw_parameters) that actually samples the source parameters. \n
        :class:`~ler.source_population.CompactBinaryPopulation`

    Examples
    --------

    Instance Attributes
    ----------
    LensGalaxyPopulation class has the following instance attributes:\n
    +-------------------------------------+----------------------------------+
    | Atrributes                          | Type                             |
    +=====================================+==================================+
    |:attr:`~cbc_pop`                     | CompactBinaryPopulation class    |
    +-------------------------------------+----------------------------------+
    |:attr:`~z_min`                       | float                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~z_max`                       | float                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~m_min`                       | float                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~m_max`                       | float                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~normalization_pdf_z`         | float                            |
    +-------------------------------------+----------------------------------+

    Instance Methods
    ----------
    LensGalaxyPopulation class has the following instance methods:\n
    +-------------------------------------+----------------------------------+
    | Methods                             | Type                             |
    +=====================================+==================================+
    |:meth:`~create_lookup_table`         | Function to create a lookup      |
    |                                     | table for the differential       |
    |                                     | comoving volume and luminosity   |
    |                                     | distance wrt redshift            |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_lens_parameters`      | Function to sample lens galaxy   |
    |                                     | parameters                       |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_lens_parameters_routine`                                 |
    +-------------------------------------+----------------------------------+
    |                                     | Function to sample lens galaxy   |
    |                                     | parameters                       |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_sl_source_parameters`                       |
    +-------------------------------------+----------------------------------+
    |                                     | Function to sample source        |
    |                                     | parameters conditioned on the    |
    |                                     | source being strongly lensed     |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_lens_redshifts`       | Function to sample lens redshifts|
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_velocity_dispersion_axis_ratio`                          |
    +-------------------------------------+----------------------------------+
    |                                     | Function to sample velocity      |
    |                                     | dispersion and axis ratio of the |
    |                                     | lens galaxy                      |
    +-------------------------------------+----------------------------------+
    |:meth:`~compute_einstein_radii`      | Function to compute the Einstein |
    |                                     | radii of the lens galaxies       |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_axis_ratio_angle_phi` | Function to sample the axis      |
    |                                     | rotation angle of the elliptical |
    |                                     | lens galaxy                      |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_galaxy_shear`         | Function to sample the lens      |
    |                                     | galaxy shear                     |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_gamma`                | Function to sample the lens      |
    |                                     | galaxy spectral index of the     |
    |                                     | density profile                  |
    +-------------------------------------+----------------------------------+
    |:meth:`~rejection_sample_lensing_probability`                           |
    +-------------------------------------+----------------------------------+
    |                                     | Function to conduct rejection    |
    |                                     | sampling wrt einstein radius     |
    +-------------------------------------+----------------------------------+
    |:meth:`~strong_lensing_optical_depth`| Function to compute the strong   |
    |                                     | lensing optical depth            |
    +-------------------------------------+----------------------------------+
    |:meth:`~get_image_properties`        | Function to get the image        |
    |                                     | properties e.g. image positions, |
    |                                     | magnifications, time delays, etc.|
    +-------------------------------------+----------------------------------+
    |:meth:`~get_lensed_snrs`             | Function to get the lensed SNRs  |
    +-------------------------------------+----------------------------------+

    """

    # Attributes
    cbc_pop = None
    """:class:`~CompactBinaryPopulation` class\n
    This is an already initialized class that contains a function (CompactBinaryPopulation.sample_gw_parameters) that actually samples the source parameters. 
    """

    z_min = None
    """`float`\n
    minimum redshift
    """
    z_max = None
    """`float`\n
    maximum redshift
    """

    m_min = None
    """`float`\n
    minimum mass in detector frame
    """

    m_max = None
    """`float`\n
    maximum mass in detector frame
    """

    normalization_pdf_z = None
    """`float`\n
    normalization constant of the pdf p(z)
    """

    def __init__(
        self,
        npool=4,
        z_min=0.001,
        z_max=10.0,
        cosmology=None,
        event_type="BBH",
        CompactBinaryPopulation_=False,
        lens_type="epl_galaxy",
        lens_functions= None,
        sampler_priors=None,
        sampler_priors_params=None,
        **kwargs
    ):
        
        self.npool = npool
        self.z_min = z_min
        self.z_max = z_max
        self.cosmo = cosmology if cosmology else cosmo
        # function initialization
        self.sample_lens_parameters_routine = getattr(
            self, lens_functions["param_sampler_type"]
        )
        self.rejection_sample_sl = getattr(
            self, lens_functions["strong_lensing_condition"]
        )  # SL: Strong Lensing
        # self.strong_lensing_optical_depth = getattr(
        #     self, lens_functions["optical_depth"]
        # )
        # self.rejection_sample_sl = self.rj_on_cross_section
        # dealing with prior functions and categorization
        (
            self.lens_param_samplers,
            self.lens_param_samplers_params,
            self.lens_sampler_names,
            self.lens_functions,
        ) = self.lens_priors_categorization(
            lens_type, sampler_priors, sampler_priors_params, lens_functions,
        )

        # initialization of CompactBinaryPopulation class
        # it also initializes the CBCSourceRedshiftDistribution class
        # list of relevant initialized instances,
        # 1. self.sample_source_redshift
        # 2. self.sample_gw_parameters
        if CompactBinaryPopulation_ == False:
            input_params = dict(
                z_min=self.z_min,
                z_max=self.z_max,
                event_type=event_type,
                event_priors=None,
                event_priors_params=None,
                spin_zero=True,
                cosmology=self.cosmo,
            )
            input_params.update(kwargs)
            # initialization of clasess
            CompactBinaryPopulation.__init__(
                self,
                z_min=self.z_min,
                z_max=self.z_max,
                event_type=input_params["event_type"],
                event_priors=input_params["event_priors"],
                event_priors_params=input_params["event_priors_params"],
                spin_zero=input_params["spin_zero"],
                cosmology=self.cosmo,
            )
        else:
            # if the classes are already initialized, then just use them
            self.__bases__ = (CompactBinaryPopulation_,)

        # initialize the image properties class
        input_params_image = dict(
            npool=npool,
            n_min_images=2,
            n_max_images=4,
            lens_model_list=["EPL_NUMBA", "SHEAR"],
        )
        input_params_image.update(kwargs)
        ImageProperties.__init__(
            self,
            npool=input_params_image["npool"],
            n_min_images=input_params_image["n_min_images"],
            n_max_images=input_params_image["n_max_images"],
            lens_model_list=input_params_image["lens_model_list"],
        )

        # cosmolgy related functions, for fast calculation through interpolation
        self.create_lookup_table_lensing(self.z_max)

        # initializing samplers
        self.sample_source_redshift_sl = self.lens_param_samplers["source_redshift_sl"]
        self.sample_lens_redshift = self.lens_param_samplers["lens_redshift"]
        self.sample_velocity_dispersion = self.lens_param_samplers[
            "velocity_dispersion"
        ]
        self.sample_axis_ratio = self.lens_param_samplers["axis_ratio"]
        self.sample_axis_rotation_angle = self.lens_param_samplers[
            "axis_rotation_angle"
        ]
        self.sample_shear = self.lens_param_samplers["shear"]
        self.sample_mass_density_spectral_index = self.lens_param_samplers[
            "mass_density_spectral_index"
        ]
        self.sample_source_parameters = self.lens_param_samplers["source_parameters"]

        return None

    def lens_priors_categorization(
        self, lens_type, sampler_priors=None, sampler_priors_params=None, lens_functions=None,
    ):
        """
        Function to categorize the lens priors/samplers

        Parameters
        ----------
            lens_type : `str`
                lens type
                e.g. 'epl_galaxy' for elliptical power-law galaxy
            sampler_priors : `dict`
                dictionary of priors
            sampler_priors_params : `dict`
                dictionary of priors parameters
        """

        if lens_type == "epl_galaxy":
            sampler_priors_ = dict(
                source_redshift_sl="strongly_lensed_source_redshifts",
                lens_redshift="lens_redshift_SDSS_catalogue",
                velocity_dispersion="velocity_dispersion_gengamma",
                axis_ratio="axis_ratio_rayleigh",
                axis_rotation_angle="axis_rotation_angle_uniform",
                shear="shear_norm",
                mass_density_spectral_index="mass_density_spectral_index_normal",
                source_parameters="source_parameters",
            )
            sampler_priors_params_ = dict(
                source_redshift_sl=None,
                lens_redshift=None,
                velocity_dispersion=dict(a=2.32 / 2.67, c=2.67),
                axis_ratio=dict(q_min=0.2),
                axis_rotation_angle=dict(phi_min=0.0, phi_max=2 * np.pi),
                shear=dict(scale=0.05),
                mass_density_spectral_index=dict(mean=2.0, std=0.2),
                source_parameters=None,
            )
            lens_functions_ = dict(
                strong_lensing_condition="rjs_with_einstein_radius",
                optical_depth="optical_depth_SIS_haris",
                param_sampler_type="sample_all_routine1",
            ),
        else:
            raise ValueError("lens_type not recognized")

        # update the priors if input is given
        if sampler_priors:
            sampler_priors_.update(sampler_priors)
        if sampler_priors_params:
            sampler_priors_params_.update(sampler_priors_params)
        if lens_functions:
            lens_functions_.update(lens_functions)

        # dict of sampler names with description
        lens_sampler_names_ = dict(
            sample_source_redshift_sl="source parameters conditioned on the source being strongly lensed",
            sample_lens_redshift="lens redshift",
            sample_velocity_dispersion="velocity dispersion of elliptical galaxy",
            sample_axis_ratio="axis ratio of elliptical galaxy",
            sample_axis_rotation_angle="axis rotation angle of elliptical galaxy    ",
            sample_shear="shear of elliptical galaxy",
            sample_mass_density_spectral_index="mass density spectral index of elliptical power-law galaxy",
            sample_source_parameters="source parameters other than redshift",
        )

        return sampler_priors_, sampler_priors_params_, lens_sampler_names_, lens_functions_

    def sample_lens_parameters(
        self,
        size=1000,
        lens_parameters_input=None,
    ):
        """
        Function to call the specific galaxy lens parameters sampler.
        """

        return self.sample_lens_parameters_routine(
            size=size, lens_parameters_input=lens_parameters_input
        )

    def sample_all_routine1(self, size=1000, lens_parameters_input=None):
        """
        Function to sample galaxy lens parameters along with the source parameters.

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        lens_parameters_input : `dict`
            dictionary of lens parameters to sample

        Returns
        -------
        lens_parameters : `dict`
            dictionary of lens parameters and source parameters (lens conditions applied): \n
            zl: lens redshifts \n
            zs: source redshifts, lensed condition applied\n
            sigma: velocity dispersions \n
            q: axis ratios \n
            theta_E: Einstein radii \n
            phi: axis rotation angle \n
            e1: ellipticity component 1 \n
            e2: ellipticity component 2 \n
            gamma1: shear component 1 \n
            gamma2: shear component 2 \n
            gamma: spectral index of the mass density distribution \n
            geocent_time: time of arrival of the unlensed signal\n
            phase: phase of the unlensed signal\n
            psi: polarization angle of the unlensed signal\n
            theta_jn: inclination angle of the unlensed signal\n
            luminosity_distance: luminosity distance of the source\n
            mass_1_source: mass 1 (larger) of the source\n
            mass_2_source: mass 2 (smaller) of the source\n
            ra: right ascension of the source\n
            dec: declination of the source\n
        """

        if lens_parameters_input is None:
            lens_parameters_input = dict()
        samplers_params = self.lens_param_samplers_params.copy()

        # Sample source redshifts from the source population
        # rejection sampled with optical depth
        zs = self.sample_source_redshift_sl(size=size)

        # Sample lens redshifts
        zl = self.sample_lens_redshift(zs=zs)

        # Sample velocity dispersions
        sigma = self.sample_velocity_dispersion(
            len(zs), param=samplers_params["velocity_dispersion"]
        )

        # Sample axis ratios
        q = self.sample_axis_ratio(sigma, param=samplers_params["axis_ratio"])

        # Compute the Einstein radii
        theta_E = self.compute_einstein_radii(sigma, zl, zs)

        # Create a dictionary of the lens parameters
        lens_parameters = {
            "zl": zl,
            "zs": zs,
            "sigma": sigma,
            "q": q,
            "theta_E": theta_E,
        }

        # Rejection sample based on the lensing probability, that is, rejection sample wrt theta_E
        lens_parameters = self.rjs_with_einstein_radius(
            lens_parameters
        )  # proportional to pi theta_E^2

        # Add the lensing parameter dictionaries together
        lens_parameters = add_dictionaries_together(
            lens_parameters, lens_parameters_input
        )

        # check the size of the lens parameters
        if len(lens_parameters["zl"]) < size:
            # Run iteratively until we have the right number of lensing parmaeters
            print("current sampled size", len(lens_parameters["zl"]))
            return self.sample_all_routine1(
                size=size, lens_parameters_input=lens_parameters
            )
        else:
            # Trim dicitionary to right size
            lens_parameters = trim_dictionary(lens_parameters, size)

            # Sample the axis rotation angle
            lens_parameters["phi"] = self.sample_axis_rotation_angle(
                size=size, param=samplers_params["axis_rotation_angle"],
            )

            # Transform the axis ratio and the angle, to ellipticities e1, e2, using lenstronomy
            lens_parameters["e1"], lens_parameters["e2"] = phi_q2_ellipticity(
                lens_parameters["phi"], lens_parameters["q"]
            )

            # Sample shears
            lens_parameters["gamma1"], lens_parameters["gamma2"] = self.sample_shear(
                size=size, param=samplers_params["shear"]
            )

            # Sample the spectral index of the mass density distribution
            lens_parameters["gamma"] = self.sample_mass_density_spectral_index(
                size=size, param=samplers_params["mass_density_spectral_index"]
            )

            # sample gravitional waves source parameter
            param = dict(zs=np.array(zs))
            if samplers_params["source_parameters"]:
                param.update(self.sample_gw_parameters(size=size))
            gw_param = self.sample_source_parameters(size=size, param=param)

            # Add source params strongly lensed to the lens params
            lens_parameters.update(gw_param)
            del (
                lens_parameters["mass_1"],
                lens_parameters["mass_2"],
            )  # remove detector frame masses

            return lens_parameters

    def strongly_lensed_source_redshifts(self, size=1000):
        """
        Function to sample source redshifts and other parameters, conditioned on the source being strongly lensed.

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample

        Returns
        -------
        redshifts : `float`
            source redshifts conditioned on the source being strongly lensed

        """

        z_max = self.z_max

        def zs_function(zs_sl):
            # get zs
            # self.sample_source_redshifts from CBCSourceRedshiftDistribution class
            zs = self.sample_source_redshifts(
                size
            )  # this function is from CompactBinaryPopulation class
            # put strong lensing condition with optical depth
            tau = self.strong_lensing_optical_depth(zs)
            tau_max = self.strong_lensing_optical_depth(z_max)[0] # tau increases with z
            r = np.random.uniform(0, tau_max, size=len(zs)) 
            # Add the strongly lensed source redshifts to the list
            # pick strongly lensed sources
            zs_sl += list(zs[r < tau])  # list concatenation

            # Check if the zs_sl are larger than requested size
            if len(zs_sl) >= size:
                # Trim list to right size
                zs_sl = zs_sl[:size]
                return zs_sl
            else:
                # Run iteratively until we have the right number of lensing parmaeters
                return zs_function(zs_sl)

        zs_sl = []

        return np.array(zs_function(zs_sl))

    def source_parameters(self, size, param=None):
        """
        Function to sample gw source parameters

        Parameters
        ----------
        size : `int`
            Number of samples to draw
        param : `dict`
            Allows to pass in parameters as dict.
            param =

        Returns
        ----------
        source_parameters : `dict`
            Dictionary of source parameters
            source_parameters.keys() = ['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'inclination', 'polarization_angle', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a_2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl']
        """

        # sample gravitional waves source parameter
        return self.sample_gw_parameters(size=size, param=param)

    def lens_redshift_SDSS_catalogue(self, zs):
        """
        Function to sample lens redshifts, conditioned on the lens being strongly lensed

        Parameters
        ----------
        zs : `float`
            source redshifts

        Returns
        -------
        zl : `float`
            lens redshifts
        """

        # lens redshift distribution
        u = njit(np.random.uniform(0, 1, size=len(zs)))
        r = (10 * u**3 - 15 * u**4 + 6 * u**5)  # See the integral of Eq. A7 of https://arxiv.org/pdf/1807.07062.pdf (cdf)
        # comoving distance to the lens galaxy
        # on the condition that lens lie between the source and the observer
        lens_galaxy_Dc = (
            self.z_to_Dc(zs) * r
        )  # corresponding element-wise multiplication between 2 arrays

        # lens redshifts
        return self.Dc_to_z(lens_galaxy_Dc)

    def velocity_dispersion_gengamma(self, size, a=2.32 / 2.67, c=2.67, param=None):
        """
        Function to sample velocity dispersion from gengamma distribution

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        a,c : `float`
            parameters of gengamma distribution
            refer to https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.gengamma.html

        Returns
        -------
        sigma : `array`
            velocity dispersion of the lens galaxy
        """

        if param:
            a = param["a"]
            c = param["c"]

        # sample velocity dispersion from gengamma distribution
        return gengamma.rvs(a, c, size=size)  # km/s

    def axis_ratio_rayleigh(self, sigma, q_min=0.2, param=None):
        """
        Function to sample axis ratio from rayleigh distribution with given velocity dispersion.

        Parameters
        ----------
        sigma : `float: array`
            velocity dispersion of the lens galaxy

        Returns
        -------
        q : `float: array`
            axis ratio of the lens galaxy
        """

        if param:
            q_min = param["q_min"]

        size = len(sigma)
        a = sigma / 161.0
        q = np.ones(size)
        idx = np.arange(size)  # idx tracker
        size_ = size

        while size_ != 0:
            # Draw the axis ratio see Appendix of https://arxiv.org/pdf/1807.07062.pdf
            s = 0.38 - 0.09177 * a[idx]
            b = rayleigh.rvs(scale=s, size=size_)
            q_ = 1.0 - b

            # Weed out axis ratios that have axis ratio below q_min
            idx2 = q_ > q_min
            q[idx[idx2]] = q_[idx2]

            # remaining idx from the original array
            # that still not have axis ratio above q_min
            idx = idx[q <= q_min]
            size_ = len(idx)

        return q

    def axis_rotation_angle_uniform(
        self, size=1000, phi_min=0.0, phi_max=2 * np.pi, param=None
    ):
        """
        Function to sample the axis rotation angle of the elliptical lens galaxy from a uniform distribution

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample

        Returns
        -------
        phi : `float`
            axis rotation angle of the elliptical lens galaxy
        """

        if param:
            phi_min = param["phi_min"]
            phi_max = param["phi_max"]

        # Draw the angles from a uniform distribution
        return np.random.uniform(phi_min, phi_max, size=size)

    def shear_norm(self, size, scale=0.05, param=None):
        """
        Function to sample the elliptical lens galaxy shear from a normal distribution

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample

        Returns
        -------
        gamma_1 : `float`
            shear component in the x-direction
        gamma_2 : `float`
            shear component in the y-direction

        """

        if param:
            scale = param["scale"]

        # Draw an external shear from a normal distribution
        return norm.rvs(size=size, scale=scale), norm.rvs(size=size, scale=scale)

    def mass_density_spectral_index_normal(
        self, size=1000, mean=2.0, std=0.2, param=None
    ):
        """
        Function to sample the lens galaxy spectral index of the mass density profile from a normal distribution

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample

        Returns
        -------
        gamma : `float`
            spectral index of the density profile

        """

        if param:
            mean = param["mean"]
            std = param["std"]

        # Draw the spectral index from a normal distribution
        return np.random.normal(loc=mean, scale=std, size=size)

    def compute_einstein_radii(self, sigma, zl, zs):
        """
        Function to compute the Einstein radii of the lens galaxies

        Parameters
        ----------
        sigma : `float`
            velocity dispersion of the lens galaxy
        zl : `float`
            lens redshifts
        zs : `float`
            source redshifts

        Returns
        -------
        theta_E : `float`
            Einstein radii of the lens galaxies in radian
        """

        # Compute the angular diameter distances
        Ds = self.angular_diameter_distance(zs)
        Dls = self.angular_diameter_distance_z1z2(zl, zs)
        # Compute the Einstein radii
        theta_E = (
            4.0 * np.pi * (sigma / const.c.to("km/s").value) ** 2 * Dls / (Ds)
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        return theta_E

    def rjs_with_einstein_radius(self, param_dict):
        """
        Function to conduct rejection sampling wrt einstein radius

        Parameters
        ----------
        param_dict : `dict`
            dictionary of lens parameters

        Returns
        -------
        lens_params : `dict`
            dictionary of lens parameters after rejection sampling
        """

        theta_E = param_dict["theta_E"]
        size = len(theta_E)
        theta_E_max = np.max(theta_E)  # maximum einstein radius
        u = np.random.uniform(0, theta_E_max**2, size=size)
        mask = u < theta_E**2

        # return the dictionary with the mask applied
        return {key: val[mask] for key, val in param_dict.items()}

    def optical_depth_SIS_haris(self, zs):
        """
        Function to compute the strong lensing optical depth (SIS). \n
        LambdaCDM(H0=70, Om0=0.3, Ode0=0.7) was used to derive the following equation.

        Parameters
        ----------
        zs : `float`
            source redshifts

        Returns
        -------
        tau : `float`
            strong lensing optical depth
        """

        # z to luminosity_distance (luminosity_distance) conversion
        Dc = self.z_to_Dc(zs) * 1e-3  # 1e-3 converts Mpc to Gpc

        return (Dc / 62.2) ** 3  # 62.2 is the critical density in Gpc

    def optical_depth_SIE(self, zs):
        """
        Function to compute the strong lensing optical depth SIE

        Parameters
        ----------
            zs : `float`
                source redshifts

        Returns
        -------
            tau : `float`
                strong lensing optical depth

        """

        # for SIE model
        # dedine a function to calculate the cross section number
        # here we already assume that we have initialized the crossection spline
        def getcrosssect_num(theta_E, q):
            fid_b_I = 10.0  # constant
            idx = q > 0.999
            q[idx] = 0.999
            idx = q < 0.1
            q[idx] = 0.1
            b_I = theta_E * np.sqrt(q)
            return self.cross_sect_spl(q) * (b_I / fid_b_I) ** 2

        # theta_E=bsis=einstein_radius_SIS  # refer to compute_einstein_radii

        return getcrosssect_num(theta_E, q) / (4 * np.pi)

    def create_lookup_table_lensing(self, z_max):

        # create a lookup table for the lens redshift draws
        r = np.linspace(0, 1, num=100)
        # inverse of cdf of the lens redshift distribution
        u = (
            10 * r**3 - 15 * r**4 + 6 * r**5
        )  # See the integral of Eq. A7 of https://arxiv.org/pdf/1807.07062.pdf (cdf)
        self.lens_redshift_sampler_helper_function = interp1d(
            u, r, kind="cubic"
        )  # Computes r(u)

        return None

    @property
    def sample_source_redshift_sl(self):
        """
        Function to sample source redshifts conditioned on the source being strongly lensed
        """
        return self._sample_source_redshift_sl

    @sample_source_redshift_sl.setter
    def sample_source_redshift_sl(self, prior):
        try:
            self._sample_source_redshift_sl = getattr(self, prior)
        except:
            self._sample_source_redshift_sl = prior

    @property
    def sample_lens_redshift(self):
        """
        Function to sample lens redshifts, conditioned on the lens being strongly lensed
        """
        return self._sample_lens_redshift

    @sample_lens_redshift.setter
    def sample_lens_redshift(self, prior):
        try:
            self._sample_lens_redshift = getattr(self, prior)
        except:
            self._sample_lens_redshift = prior

    @property
    def sample_velocity_dispersion(self):
        """
        Function to sample velocity dispersion from gengamma distribution
        """
        return self._sample_velocity_dispersion

    @sample_velocity_dispersion.setter
    def sample_velocity_dispersion(self, prior):
        try:
            self._sample_velocity_dispersion = getattr(self, prior)
        except:
            self._sample_velocity_dispersion = prior

    @property
    def sample_axis_ratio(self):
        """
        Function to sample axis ratio from rayleigh distribution with given velocity dispersion.
        """
        return self._sample_axis_ratio

    @sample_axis_ratio.setter
    def sample_axis_ratio(self, prior):
        try:
            self._sample_axis_ratio = getattr(self, prior)
        except:
            self._sample_axis_ratio = prior

    @property
    def sample_axis_rotation_angle(self):
        """
        Function to sample the axis rotation angle of the elliptical lens galaxy from a uniform distribution
        """
        return self._sample_axis_rotation_angle

    @sample_axis_rotation_angle.setter
    def sample_axis_rotation_angle(self, prior):
        try:
            self._sample_axis_rotation_angle = getattr(self, prior)
        except:
            self._sample_axis_rotation_angle = prior

    @property
    def sample_shear(self):
        """
        Function to sample the elliptical lens galaxy shear from a normal distribution
        """
        return self._sample_shear

    @sample_shear.setter
    def sample_shear(self, prior):
        try:
            self._sample_shear = getattr(self, prior)
        except:
            self._sample_shear = prior

    @property
    def sample_mass_density_spectral_index(self):
        """
        Function to sample the lens galaxy spectral index of the mass density profile from a normal distribution
        """
        return self._sample_mass_density_spectral_index

    @sample_mass_density_spectral_index.setter
    def sample_mass_density_spectral_index(self, prior):
        try:
            self._sample_mass_density_spectral_index = getattr(self, prior)
        except:
            self._sample_mass_density_spectral_index = prior

    @property
    def sample_source_parameters(self):
        """
        Function to sample source parameters conditioned on the source being strongly lensed
        """
        return self._sample_source_parameters

    @sample_source_parameters.setter
    def sample_source_parameters(self, prior):
        try:
            self._sample_source_parameters = getattr(self, prior)
        except:
            self._sample_source_parameters = prior
