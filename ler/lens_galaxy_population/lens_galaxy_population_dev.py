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
from astropy.cosmology import Planck18
from astropy import constants as const

# the following .py file will be called if they are not given in the class initialization
from ..gw_source_population import CompactBinaryPopulation
from ..image_properties import ImageProperties
from ..utils import add_dictionaries_together, trim_dictionary


class LensGalaxyPopulation(CompactBinaryPopulation, ImageProperties):
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
    |:meth:`~sample_strongly_lensed_source_parameters`                       |
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
        CompactBinaryPopulation_=False,
        lens_type="epl_galaxy",
        lens_functions=dict(
            strong_lensing_condition="rjs_with_einstein_radius",
            optical_depth="optical_depth_SIS",
            param_sampler_type="sample_all_routine1",
        ),
        sampler_priors=None,
        sampler_priors_params=None,
        **kwargs
    ):
        self.npool = npool
        self.sample_lens_parameters_routine = getattr(
            self, lens_functions["param_sampler_type"]
        )
        self.rejection_sample_sl = getattr(
            self, lens_functions["strong_lensing_condition"]
        )  # SL: Strong Lensing
        self.strong_lensing_optical_depth = getattr(
            self, lens_functions["optical_depth"]
        )
        # self.rejection_sample_SL = self.rj_on_cross_section

        if CompactBinaryPopulation_ == False:
            input_params = dict(
                z_min=0.0,
                z_max=10.0,
                event_type="BBH",
                event_priors=None,
                event_priors_params=None,
                spin_zero=True,
                cosmology=None,
            )
            input_params.update(kwargs)
            # initialization of clasess
            CompactBinaryPopulation.__init__(
                self,
                z_min=input_params["z_min"],
                z_max=input_params["z_max"],
                event_type=input_params["event_type"],
                event_priors=input_params["event_priors"],
                event_priors_params=input_params["event_priors_params"],
                spin_zero=input_params["spin_zero"],
                cosmology=input_params["cosmology"],
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

        # dict to strore lens parameters
        self.lens_parameters = dict()
        # dealing with prior functions and categorization
        (
            self.lens_param_samplers,
            self.lens_param_samplers_params,
            self.lens_sampler_names,
        ) = self.lens_priors_categorization(
            lens_type, sampler_priors, sampler_priors_params
        )

        # cosmolgy related functions, for fast calculation through interpolation
        self.create_lookup_table_lensing(self.z_max)

        # initializing samplers
        self.sample_source_parameters = self.lens_param_samplers["source_parameters"]
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

        return None

    def lens_priors_categorization(
        self, lens_type, sampler_priors=None, sampler_priors_params=None
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
                source_parameters="sample_strongly_lensed_source_parameters",
                lens_redshift="lens_redshift_SDSS_catalogue",
                velocity_dispersion="velocity_dispersion_gengamma",
                axis_ratio="axis_ratio_rayleigh",
                axis_rotation_angle="axis_rotation_angle_uniform",
                shear="shear_norm",
                mass_density_spectral_index="mass_density_spectral_index_normal",
            )
            sampler_priors_params_ = dict(
                source_parameters=None,
                lens_redshift=dict(sampled_zs=True),
                velocity_dispersion=dict(a=2.32 / 2.67, c=2.67),
                axis_ratio=dict(sampled_sigma=True),
                axis_rotation_angle=dict(phi_min=0.0, phi_max=2 * np.pi),
                shear=dict(scale=0.05),
                mass_density_spectral_index=dict(mean=2.0, std=0.2),
            )
        else:
            raise ValueError("lens_type not recognized")

        # update the priors if input is given
        if sampler_priors:
            sampler_priors_.update(sampler_priors)
        if sampler_priors_params:
            sampler_priors_params_.update(sampler_priors_params)

        # dict of sampler names with description
        lens_sampler_names_ = dict(
            sample_source_parameters="source parameters conditioned on the source being strongly lensed",
            sample_lens_redshift="lens redshift",
            sample_velocity_dispersion="velocity dispersion of elliptical galaxy",
            sample_axis_ratio="axis ratio of elliptical galaxy",
            sample_axis_rotation_angle="axis rotation angle of elliptical galaxy    ",
            sample_shear="shear of elliptical galaxy",
            sample_mass_density_spectral_index="mass density spectral index of elliptical power-law galaxy",
        )

        return sampler_priors_, sampler_priors_params_, lens_sampler_names_

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
            dictionary of lens parameters and source parameters (lens conditions applied)
        """

        if lens_parameters_input is None:
            lens_parameters_input = dict()

        # Sample source redshifts from the source population
        # rejection sampled with optical depth
        source_params_strongly_lensed = self.sample_strongly_lensed_source_parameters(
            size=size
        )
        zs = source_params_strongly_lensed["zs"]

        # Sample lens redshifts
        zl = self.lens_redshift_SDSS_catalogue(size=len(zs), zs=zs)
        # Sample velocity dispersions
        sigma = self.velocity_dispersion_gengamma(len(zs))
        # Sample axis ratios
        q = self.axis_ratio_rayleigh(len(zs), sigma)

        # Compute the Einstein radii
        # print("sampling einstein radius...")
        theta_E = self.compute_einstein_radii(sigma, zl, zs)

        # Sample the axis ratio angle
        # print("sampling axis ratio and angle...")
        axis_ratio_angle_phi = self.axis_ratio_rayleigh(size=size, sigma=sigma)

        # Transform the axis ratio and the angle, to ellipticities e1, e2, using lenstronomy
        e1, e2 = phi_q2_ellipticity(axis_ratio_angle_phi, q)

        # Sample shears
        # print("sampling external shears...")
        gamma1, gamma2 = self.shear_norm(size=size)

        # Sample the spectral index of the mass density distribution
        # print("sampling spectral index...")
        gamma = self.mass_density_spectral_index_normal(size=size)

        # Create a dictionary of the lens parameters
        lens_parameters = {
            "zl": zl,
            "zs": zs,
            "sigma": sigma,
            "q": q,
            "e1": e1,
            "e2": e2,
            "gamma1": gamma1,
            "gamma2": gamma2,
            "theta_E": theta_E,
            "gamma": gamma,
        }

        # Add source params strongly lensed to the lens params
        lens_parameters.update(source_params_strongly_lensed)

        # Rejection sample based on the lensing probability, that is, rejection sample wrt theta_E
        lens_parameters = self.rjs_with_einstein_radius(
            lens_parameters
        )  # proportional to pi theta_E^2

        # Add the lensing parameter dictionaries together
        # Note: This is done after rejection sampling, so that we don't have to worry about the size of the dictionary
        lens_parameters = add_dictionaries_together(
            lens_parameters, lens_parameters_input
        )

        # Check if the lens are larger than requested size
        # len_ = len(lens_parameters['zl'])
        # print(f'len(zs) = {len_}')
        if len(lens_parameters["zl"]) >= size:
            # Trim dicitionary to right size
            lens_parameters = trim_dictionary(lens_parameters, size)
            return lens_parameters
        else:
            # Run iteratively until we have the right number of lensing parmaeters
            # print("current sampled size", len(lens_parameters['zl']))
            return self.sample_all_routine1(
                size=size, lens_parameters_input=lens_parameters
            )

    def sample_all_routine2(self, size=1000, lens_parameters_input=None, **kwargs):
        """
        Under construction!!
        Function to sample galaxy lens parameters along with the source parameters.

        Parameters
        ----------
        size : `int`
            number of samples in each lens parameter
        lens_parameters_input : `dict`
            dictionary of lens parameters which keeps updating as the function runs. This is used for rejection sampling.
        **kwargs : `dict`
            user provided lens parameters. The size of the array should be equal to size input. \n
            e.g. sigma = np.array([100, 200, 300, 400, 500]) \n
            allowed parameters: zs, zl, sigma, q, e1, e2, gamma1, gamma2, gamma, theta_E

        Returns
        -------
        lens_parameters : `dict`
            dictionary of lens parameters and source parameters (lens conditions applied)
        """

        if lens_parameters_input is None:
            lens_parameters_input = dict()  # initialize dictionary to store parameters

        param_names = list(self.lens_param_samplers.keys())
        samplers_params = list(self.lens_param_samplers_params.values())
        # make sure the order is correct
        sampler_names = list(self.lens_sampler_names.keys())
        param_dict = self.lens_parameters  # initialize dictionary to store parameters

        for name, sampler, param in zip(param_names, sampler_names, samplers_params):
            if name not in kwargs:
                # Sample the parameter using the specified sampler function
                # sample stored in param_dict object
                getattr(self, sampler)(size=size, param=param);
            else:
                # Use the provided value from kwargs
                param_dict[name] = kwargs[name]

        # add other lensing parameters
        param_dict["theta_E"] = self.compute_einstein_radii(
            param_dict["sigma"], param_dict["zl"], param_dict["zs"]
        )
        param_dict["e1"], param_dict["e2"] = phi_q2_ellipticity(
            param_dict["axis_rotation_angle"], param_dict["q"]
        )

        # Rejection sample based on the lensing probability
        param_dict = self.rejection_sample_sl(param_dict)

        # Add the lensing parameter dictionaries together
        param_dict = add_dictionaries_together(lens_parameters_input, param_dict)
        # Check if the lens are larger than requested size
        if len(param_dict["zl"]) >= size:
            # Trim dicitionary to right size
            param_dict = trim_dictionary(param_dict, size)
            return param_dict
        else:
            # Run iteratively until we have the right number of lensing parmaeters
            print("current sampled size", len(param_dict['zl']))
            return self.sample_all_routine2(
                size=size, lens_parameters_input=param_dict
            )


    def sample_strongly_lensed_source_parameters(self, size=1000, param=None):
        """
        Function to sample source redshifts and other parameters, conditioned on the source being strongly lensed.

        Parameters
        ----------
            size : `int`
                number of lens parameters to sample

        Returns
        -------
            gw_param_strongly_lensed : `dict`
                dictionary of source parameters. `zs` is sampled considering the merger rate density at source frame, comoving volume and strong lensing optical depth. \n
                e.g. gw_param_strongly_lensed.keys() = ['mass_1', 'mass_2', 'mass_1_source', 'mass_2_source', 'zs', 'luminosity_distance', 'iota', 'psi', 'phase', 'geocent_time', 'ra', 'dec', 'a_1', 'a2', 'tilt1', 'tilt2', 'phi12', 'phi_jl']

        """

        z_max = self.z_max

        def zs_strongly_lensed_function(zs_strongly_lensed):
            # get zs
            zs = self.sample_source_redshifts(
                size
            )  # this function is from CompactBinaryPopulation class
            # put strong lensing condition with optical depth
            tau = self.strong_lensing_optical_depth(zs)
            tau_max = np.max(
                self.strong_lensing_optical_depth(z_max)
            )  # tau increases with z
            r = np.random.uniform(0, tau_max, size=len(zs))
            pick_strongly_lensed = r < tau  # pick strongly lensed sources
            # Add the strongly lensed source redshifts to the list
            zs_strongly_lensed += list(zs[pick_strongly_lensed])  # list concatenation

            # Check if the zs_strongly_lensed are larger than requested size
            if len(zs_strongly_lensed) >= size:
                # Trim list to right size
                zs_strongly_lensed = zs_strongly_lensed[:size]
                return zs_strongly_lensed
            else:
                # Run iteratively until we have the right number of lensing parmaeters
                return zs_strongly_lensed_function(zs_strongly_lensed)

        zs_strongly_lensed = []
        zs_ = np.array(zs_strongly_lensed_function(zs_strongly_lensed))
        # gravitional waves source parameter sampling
        gw_param_strongly_lensed = self.sample_gw_parameters(
            size=size, zs=np.array(zs_)
        )

        self.lens_parameters.update(gw_param_strongly_lensed)

        return gw_param_strongly_lensed

    def lens_redshift_SDSS_catalogue(self, size, zs=None, sampled_zs=True, param=None):
        """
        Function to sample lens redshifts, conditioned on the lens being strongly lensed
        Input parameters:
            zs : source redshifts
        Output parameters:
            zl : lens redshifts
        """

        if param:
            try:
                zs = self.lens_parameters["zs"]
            except:
                raise print("zs is not yet sampled")
        else:
            if zs is None:
                try:
                    zs = self.lens_parameters["zs"]
                    print("zs is not given, using the sampled zs from self.lens_parameters")
                except:
                    raise ValueError("zs is not given")

        # lens redshift distribution
        r = self.lens_redshift_sampler_helper_function(
            np.random.uniform(0, 1, size=len(zs))
        )
        # comoving distance to the lens galaxy
        # on the condition that lens lie between the source and the observer
        lens_galaxy_Dc = (
            self.z_to_Dc(zs) * r
        )  # corresponding element-wise multiplication between 2 arrays
        # lens redshift
        zl = self.Dc_to_z(lens_galaxy_Dc)[:size]  # 2D array

        self.lens_parameters["zl"] = zl
        return zl

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
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. param = dict(a=2.32 / 2.67, c=2.67)

        Returns
        -------
        sigma : `array`
            velocity dispersion of the lens galaxy
        """

        if param:
            a = param["a"]
            c = param["c"]

        sigma = gengamma.rvs(a, c, size=size)
        self.lens_parameters["sigma"] = sigma
        return sigma

    def axis_ratio_rayleigh(self, size, sigma=None, sampled_sigma=True, param=None):
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
            try:
                sigma = self.lens_parameters["sigma"]
            except:
                raise print("sigma is not yet sampled")
        else:
            if sigma is None:
                try:
                    sigma = self.lens_parameters["sigma"]
                    print("sigma is not given, using the sampled sigma from self.lens_parameters")
                except:
                    raise ValueError("sigma is not given")
                
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

            # Weed out sigmas and axis ratios that have axis ratio below 0.2
            idx2 = q_ > 0.2
            q[idx[idx2]] = q_[idx2]

            # remaining idx from the original array
            # that still not have axis ratio above 0.2
            idx = idx[q < 0.2]
            size_ = len(idx)

        self.lens_parameters["q"] = q
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
        # Draw the angles
        phi = np.random.uniform(phi_min, phi_max, size=size)
        self.lens_parameters["axis_rotation_angle"] = phi
        return phi

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

        # Draw an external shear
        gamma_1 = norm.rvs(size=size, scale=scale)
        gamma_2 = norm.rvs(size=size, scale=scale)
        self.lens_parameters["gamma_1"] = gamma_1
        self.lens_parameters["gamma_2"] = gamma_2
        return gamma_1, gamma_2

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
        gamma = np.random.normal(loc=mean, scale=std, size=size)
        self.lens_parameters["gamma"] = gamma
        return gamma

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
        theta_E_max = np.max(theta_E)
        u = np.random.uniform(0, theta_E_max**2, size=size)
        mask = u < theta_E**2

        lens_params = {key: val[mask] for key, val in param_dict.items()}
        return lens_params

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

    def optical_depth_SIS(self, zs):
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
        # For SIS model
        # z to luminosity_distance (luminosity_distance) conversion
        Dc = self.z_to_Dc(zs) * 1e-3  # 1e-3 converts Mpc to Gpc

        return (Dc / 62.2) ** 3

    def create_lookup_table_lensing(self, z_max):
        """
        Functions to create lookup tables
        1. Redshift to co-moving distance.
        2. Co-moving distance to redshift.
        3. Redshift to angular diameter distance.
        4. Lens redshift sampler helper function.

        Parameters
        ----------
        z_max : `float`
            maximum redshift

        Attributes
        ----------
        z_to_Dc : `scipy.interpolate.interpolate`
            redshift to co-moving distance
        Dc_to_z : `scipy.interpolate.interpolate`
            co-moving distance to redshift
        angular_diameter_distance_z1z2 : `lambda function`
            redshift to angular diameter distance (between two redshifts)
        angular_diameter_distance : `lambda function`
            redshift to angular diameter distance
        lens_redshift_sampler_helper_function : `scipy.interpolate.interpolate`
            lens redshift sampler helper function
        """

        # initialing cosmological functions for fast calculation through interpolation
        # z_min is fixed to 0.0, as the lens redshifts are drawn between 0.0 and z_max
        z = np.linspace(0.0, z_max, 500)  # red-shift
        Dc = self.cosmo.comoving_distance(z).value  # co-moving distance in Mpc
        self.z_to_Dc = interp1d(z, Dc, kind="cubic")
        self.Dc_to_z = interp1d(Dc, z, kind="cubic")

        # for angular diameter distance
        quad_ = []  # refers to integration (with quad algorithm) from scipy
        for zl in z:
            quad_.append(
                quad(
                    self.cosmo._inv_efunc_scalar,
                    0.0,
                    zl,
                    args=self.cosmo._inv_efunc_scalar_args,
                )[0]
            )
        quad_ = np.array(quad_)
        quad_int = interp1d(z, np.array(quad_), kind="cubic")

        H0d = self.cosmo._hubble_distance.value
        self.angular_diameter_distance_z1z2 = (
            lambda zl0, zs0: H0d * (quad_int(zs0) - quad_int(zl0)) / (zs0 + 1.0)
        )
        self.angular_diameter_distance = lambda zs0: H0d * quad_int(zs0) / (zs0 + 1.0)

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
