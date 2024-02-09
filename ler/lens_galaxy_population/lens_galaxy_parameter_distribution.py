# -*- coding: utf-8 -*-
"""
This module contains the LensGalaxyPopulation class, which is used to sample lens galaxy parameters, source parameters conditioned on the source being strongly lensed. \n
The class inherits from the ImageProperties class, which is used calculate image properties (magnification, timedelays, source position, image position, morse phase). \n
Either the class takes in initialized CompactBinaryPopulation class as input or inherits the CompactBinaryPopulation class with default params (if no input) \n
"""

import warnings
warnings.filterwarnings("ignore")
import numpy as np
from numba import njit
from scipy.integrate import quad
from lenstronomy.Util.param_util import phi_q2_ellipticity

# for redshift to luminosity distance conversion
from astropy.cosmology import LambdaCDM

# the following .py file will be called if they are not given in the class initialization
from ..gw_source_population import CBCSourceParameterDistribution
from .optical_depth import OpticalDepth
from ..image_properties import ImageProperties
from ..utils import add_dictionaries_together, trim_dictionary
from .jit_functions import phi_cut_SIE, velocity_dispersion_z_dependent, lens_redshift_SDSS_catalogue


class LensGalaxyParameterDistribution(CBCSourceParameterDistribution, ImageProperties, OpticalDepth):
    """
    Class to sample lens galaxy parameters, source parameters conditioned on the source being strongly lensed, and image properties \n

    Parameters
    ----------
    npool : `int`
        number of processors to use
    z_min : `float`
        minimum redshift
    z_max : `float`
        maximum redshift    
    cosmology : `astropy.cosmology`
        Cosmology to use
        default: None/astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)
    event_type : `str`
        Type of event to generate.
        e.g. 'BBH', 'BNS', 'NSBH'
        default: 'BBH'
    lens_type : `str`
        Type of lens galaxy to generate.
        default: 'epl_galaxy'
    lens_functions, lens_priors, lens_priors_params : `dict`, `dict`, `dict`
        dictionary of lens functions, priors, and priors parameters
        Check for default/available lens functions, priors and corresponding input parameters by running,\n
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> print(lens.lens_functions)
        >>> print(lens.lens_priors)
        >>> print(lens.lens_priors_params)
    directory : `str`
        directory to store the interpolators
        default: './interpolator_pickle'
    **kwargs : 
        keyword arguments to pass to the parent classes

    Examples
    --------
    >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
    >>> lens = LensGalaxyParameterDistribution()
    >>> lensed_params = lens.sample_lens_parameters(size=1000)
    >>> lensed_params.keys()

    Instance Attributes
    ----------
    LensGalaxyPopulation class has the following instance attributes:\n
    +-------------------------------------+----------------------------------+
    | Atrributes                          | Type                             |
    +=====================================+==================================+
    |:attr:`~npool`                       | `int`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~z_min`                       | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`~z_max`                       | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`~cosmo`                       | `astropy.cosmology`              |
    +-------------------------------------+----------------------------------+
    |:attr:`~event_type`                  | `str`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~directory`                   | `str`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`~create_new_interpolator`     | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~lens_param_samplers`         | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~lens_param_samplers_params`  | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~lens_sampler_names`          | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~lens_functions`              | `dict`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`~normalization_pdf_z_lensed`  | `float`                          |
    +-------------------------------------+----------------------------------+

    Instance Methods
    ----------
    LensGalaxyPopulation class has the following instance methods:\n
    +-------------------------------------+----------------------------------+
    | Methods                             | Type                             |
    +=====================================+==================================+
    |:meth:`~sample_lens_parameters`      | Function to call the specific    |
    |                                     | galaxy lens parameters sampler   |
    |                                     | routine.                         |
    +-------------------------------------+----------------------------------+
    |:meth:`~sample_all_routine`          | Function to sample galaxy lens   |
    |                                     | parameters along with the source |
    |                                     | parameters.                      |
    +-------------------------------------+----------------------------------+
    |:meth:`~strongly_lensed_source_redshifts`                               |
    +-------------------------------------+----------------------------------+
    |                                     | Function to sample source        |
    |                                     | redshifts conditioned on the     |
    |                                     | source being strongly lensed     |
    +-------------------------------------+----------------------------------+
    |:meth:`~source_parameters`           | Function to sample gw source     |
    |                                     | parameters                       |
    +-------------------------------------+----------------------------------+
    |:meth:`~lens_redshift_SDSS_catalogue`| Function to sample lens          |
    |                                     | redshifts, conditioned on the    |
    |                                     | lens being strongly lensed       |
    +-------------------------------------+----------------------------------+
    |:meth:`~axis_rotation_angle_uniform` | Function to sample the axis      |
    |                                     | rotation angle of the elliptical |
    |                                     | lens galaxy from a uniform       |
    |                                     | distribution                     |
    +-------------------------------------+----------------------------------+
    |:meth:`~shear_norm`                  | Function to sample the           |
    |                                     | elliptical lens galaxy shear     |
    |                                     | from a normal distribution       |
    +-------------------------------------+----------------------------------+
    |:meth:`~mass_density_spectral_index_normal`                             |
    +-------------------------------------+----------------------------------+
    |                                     | Function to sample the lens      |
    |                                     | galaxy spectral index of the     |
    |                                     | mass density profile from a      |
    |                                     | normal distribution              |
    +-------------------------------------+----------------------------------+
    |:meth:`~compute_einstein_radii`      | Function to compute the Einstein |
    |                                     | radii of the lens galaxies       |
    +-------------------------------------+----------------------------------+
    |:meth:`~rjs_with_cross_section_SIE`  | Function to conduct rejection    |
    |                                     | sampling wrt einstein radius     |
    +-------------------------------------+----------------------------------+
    |:meth:`~rjs_with_cross_section_SIE`  | Function to conduct rejection    |
    |                                     | sampling wrt cross_section       |
    +-------------------------------------+----------------------------------+
    |:attr:`~rejection_sample_sl`         | Function to conduct rejection    |
    |                                     | sampling with the given rejection|
    |                                     | sampling function                |
    +-------------------------------------+----------------------------------+
    |:attr:`~sample_source_redshift_sl`   | Function to sample source        |
    |                                     | redshifts conditioned on the     |
    |                                     | source being strongly lensed     |
    +-------------------------------------+----------------------------------+
    |:attr:`~sample_lens_redshift`        | Function to sample lens          |
    |                                     | redshifts, conditioned on the    |
    |                                     | lens being strongly lensed       |
    +-------------------------------------+----------------------------------+
    |:attr:`~sample_axis_rotation_angle`  | Function to sample the axis      |
    |                                     | rotation angle of the elliptical |
    |                                     | lens galaxy from a uniform       |
    |                                     | distribution                     |
    +-------------------------------------+----------------------------------+
    |:attr:`~sample_shear`                | Function to sample the           |
    |                                     | elliptical lens galaxy shear     |
    |                                     | from a normal distribution       |
    +-------------------------------------+----------------------------------+
    |:attr:`~sample_mass_density_spectral_index`                             |
    +-------------------------------------+----------------------------------+
    |                                     | Function to sample the lens      |
    |                                     | galaxy spectral index of the     |
    |                                     | mass density profile from a      |
    |                                     | normal distribution              |
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
        z_min=0.0,
        z_max=10.0,
        cosmology=None,
        event_type="BBH",
        lens_type="epl_galaxy",
        lens_functions= None,
        lens_priors=None,
        lens_priors_params=None,
        directory="./interpolator_pickle",
        create_new_interpolator=False,
        **kwargs
    ):
        
        self.npool = npool
        self.z_min = z_min
        self.z_max = z_max
        self.cosmo = cosmology if cosmology else LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        self.event_type = event_type
        self.directory = directory
        # initialize the interpolator's parameters
        self.create_new_interpolator = dict(
            redshift_distribution=dict(create_new=False, resolution=500),
            z_to_luminosity_distance=dict(create_new=False, resolution=500),
            velocity_dispersion=dict(create_new=False, resolution=500),
            axis_ratio=dict(create_new=False, resolution=500),
            optical_depth=dict(create_new=False, resolution=100),
            z_to_Dc=dict(create_new=False, resolution=500),
            Dc_to_z=dict(create_new=False, resolution=500),
            angular_diameter_distance=dict(create_new=False, resolution=500),
            differential_comoving_volume=dict(create_new=False, resolution=500),
            Dl_to_z=dict(create_new=False, resolution=500),
        )
        if isinstance(create_new_interpolator, dict):
            self.create_new_interpolator.update(create_new_interpolator)
        elif create_new_interpolator is True:
            self.create_new_interpolator = dict(
                redshift_distribution=dict(create_new=True, resolution=500),
                z_to_luminosity_distance=dict(create_new=True, resolution=500),
                velocity_dispersion=dict(create_new=True, resolution=500),
                axis_ratio=dict(create_new=True, resolution=500),
                optical_depth=dict(create_new=True, resolution=100),
                z_to_Dc=dict(create_new=True, resolution=500),
                Dc_to_z=dict(create_new=True, resolution=500),
                angular_diameter_distance=dict(create_new=True, resolution=500),
                differential_comoving_volume=dict(create_new=True, resolution=500),
                Dl_to_z=dict(create_new=True, resolution=500),
            )

        # dealing with prior functions and categorization
        self.lens_param_samplers, self.lens_param_samplers_params, self.lens_sampler_names, self.lens_functions = self.lens_priors_categorization(lens_type, lens_priors,
        lens_priors_params, lens_functions)
        
        # function initialization
        self.sample_lens_parameters_routine = getattr(self, self.lens_functions['param_sampler_type'])
        self.rejection_sample_sl = getattr(self, self.lens_functions['strong_lensing_condition'])

        # initializing parent classes
        self.class_initialization_lens(params=kwargs);

        # initializing samplers
        # self.sample_velocity_dispersion and self.sample_axis_ratio are initialized in OpticalDepth class
        self.sample_source_redshift_sl = self.lens_param_samplers["source_redshift_sl"]
        self.sample_lens_redshift = self.lens_param_samplers["lens_redshift"]
        # self.sample_axis_ratio = self.lens_param_samplers["axis_ratio"]
        self.sample_axis_rotation_angle = self.lens_param_samplers[
            "axis_rotation_angle"
        ]
        self.sample_shear = self.lens_param_samplers["shear"]
        self.sample_mass_density_spectral_index = self.lens_param_samplers[
            "mass_density_spectral_index"
        ]
        self.sample_source_parameters = self.lens_param_samplers["source_parameters"]

        # extra care to the velocity dispersion sampler
        if self.lens_param_samplers["velocity_dispersion"]  == "velocity_dispersion_ewoud":
            vd_inv_cdf = self.vd_inv_cdf
            zl_list = self.zl_list
            self.sample_velocity_dispersion = lambda size, zl: velocity_dispersion_z_dependent(size=size, zl=zl, zl_list=zl_list, vd_inv_cdf=vd_inv_cdf)

        # To find the normalization constant of the pdf p(z)
        # this under the assumption that the event is strongly lensed
        # Define the merger-rate density function
        pdf_unnormalized_ = lambda z: self.merger_rate_density_src_frame(np.array([z])) * self.strong_lensing_optical_depth(np.array([z]))
        pdf_unnormalized = lambda z: pdf_unnormalized_(z)[0]

        self.normalization_pdf_z_lensed = quad(
            pdf_unnormalized,
            self.z_min,
            self.z_max
        )[0]

    def class_initialization_lens(self, params=None):
        """
        Function to initialize the parent classes

        Parameters
        ----------
        params : `dict`
            dictionary of parameters to initialize the parent classes
        """

        # initialization of CompactBinaryPopulation class
        # it also initializes the CBCSourceRedshiftDistribution class
        # list of relevant initialized instances,
        # 1. self.sample_source_redshift
        # 2. self.sample_gw_parameters
        input_params = dict(
            source_priors=None,
            source_priors_params=None,
            spin_zero=True,
            spin_precession=False,
        )
        input_params.update(params)
        # initialization of clasess
        CBCSourceParameterDistribution.__init__(
            self,
            z_min=self.z_min,
            z_max=self.z_max,
            event_type=self.event_type,
            source_priors=input_params["source_priors"],
            source_priors_params=input_params["source_priors_params"],
            spin_zero=input_params["spin_zero"],
            cosmology=self.cosmo,
            spin_precession=input_params["spin_precession"],
            directory=self.directory,
            create_new_interpolator=self.create_new_interpolator,
        )

        # initialize the optical depth class
        # follwing attributes are initialized
        # 1. self.strong_lensing_optical_depth
        # 2. self.sample_velocity_dispersion
        # 3. self.sample_axis_ratio
        OpticalDepth.__init__(
            self,
            npool=self.npool,
            z_min=self.z_min,
            z_max=self.z_max,
            optical_depth_function=self.lens_functions["optical_depth"],
            sampler_priors=dict(
                velocity_dispersion=self.lens_param_samplers["velocity_dispersion"],
                axis_ratio=self.lens_param_samplers["axis_ratio"],
            ),
            sampler_priors_params=dict(
                velocity_dispersion=self.lens_param_samplers_params[
                    "velocity_dispersion"],
                axis_ratio=self.lens_param_samplers_params["axis_ratio"],
            ),
            cosmology=self.cosmo,
            directory=self.directory,
            create_new_interpolator=self.create_new_interpolator,
        )

        # initialize the image properties class
        input_params_image = dict(
            n_min_images=2,
            n_max_images=4,
            geocent_time_min=1126259462.4,
            geocent_time_max=1126259462.4+365*24*3600*100,
            lens_model_list=["EPL_NUMBA", "SHEAR"],
        )
        input_params_image.update(params)
        ImageProperties.__init__(
            self,
            npool=self.npool,
            z_min=self.z_min,
            z_max=self.z_max,
            n_min_images=input_params_image["n_min_images"],
            n_max_images=input_params_image["n_max_images"],
            lens_model_list=input_params_image["lens_model_list"],
            cosmology=self.cosmo,
            geocent_time_min=input_params_image["geocent_time_min"],
            geocent_time_max=input_params_image["geocent_time_max"],
            spin_zero=input_params["spin_zero"],
            spin_precession=input_params["spin_precession"],
            directory=self.directory,
            create_new_interpolator=self.create_new_interpolator,
        )

    def lens_priors_categorization(
        self, lens_type, lens_priors=None, lens_priors_params=None, lens_functions=None,
    ):
        """
        Function to categorize the lens priors/samplers

        Parameters
        ----------
            lens_type : `str`
                lens type
                e.g. 'epl_galaxy' for elliptical power-law galaxy
            lens_priors : `dict`
                dictionary of priors
            lens_priors_params : `dict`
                dictionary of priors parameters
            lens_functions : `dict`
                dictionary of lens functions

        Returns
        -------
        lens_priors_ : `dict`
            dictionary of priors
        lens_priors_params_ : `dict`
            dictionary of priors parameters
        lens_sampler_names_ : `dict`
            dictionary of sampler names
        lens_functions_ : `dict`
            dictionary of lens functions
        """

        if lens_type == "epl_galaxy":
            lens_priors_ = dict(
                source_redshift_sl="strongly_lensed_source_redshifts",
                lens_redshift="lens_redshift_SDSS_catalogue",
                velocity_dispersion="velocity_dispersion_ewoud",
                axis_ratio="axis_ratio_rayleigh",
                axis_rotation_angle="axis_rotation_angle_uniform",
                shear="shear_norm",
                mass_density_spectral_index="mass_density_spectral_index_normal",
                source_parameters="sample_gw_parameters",
            )
            lens_priors_params_ = dict(
                source_redshift_sl=None,
                lens_redshift=None,
                velocity_dispersion=None,
                axis_ratio=dict(q_min=0.2, q_max=1.),
                axis_rotation_angle=dict(phi_min=0.0, phi_max=2 * np.pi),
                shear=dict(scale=0.05),
                mass_density_spectral_index=dict(mean=2.0, std=0.2),
                source_parameters=None,
            )
            lens_functions_ = dict(
                strong_lensing_condition="rjs_with_cross_section_SIE",
                optical_depth="optical_depth_SIE_hemanta",
                param_sampler_type="sample_all_routine",
            )
        else:
            raise ValueError("lens_type not recognized")

        # update the priors if input is given
        if lens_priors:
            lens_priors_.update(lens_priors)
        if lens_priors_params:
            lens_priors_params_.update(lens_priors_params)
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

        return(lens_priors_, lens_priors_params_, lens_sampler_names_, lens_functions_)

    def sample_lens_parameters(
        self,
        size=1000,
        lens_parameters_input=None,
    ):
        """
        Function to call the specific galaxy lens parameters sampler routine.
        """

        return self.sample_lens_parameters_routine(
            size=size, lens_parameters_input=lens_parameters_input
        )

    def sample_all_routine(self, size=1000, lens_parameters_input=None):
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

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> lens.sample_all_routine(size=1000)
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
        try:
            sigma = self.sample_velocity_dispersion(len(zs))
        except:
            sigma = self.sample_velocity_dispersion(len(zs), zl)

        # Sample axis ratios
        try:
            q = self.sample_axis_ratio(sigma)
        except:
            q = self.sample_axis_ratio(len(sigma))

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
        lens_parameters = self.rejection_sample_sl(
            lens_parameters
        )  # proportional to pi theta_E^2

        # Add the lensing parameter dictionaries together
        lens_parameters = add_dictionaries_together(
            lens_parameters, lens_parameters_input
        )

        # check the size of the lens parameters
        if len(lens_parameters["zl"]) < size:
            # Run iteratively until we have the right number of lensing parmaeters
            # print("current sampled size", len(lens_parameters["zl"]))
            return self.sample_all_routine(
                size=size, lens_parameters_input=lens_parameters
            )
        else:
            # Trim dicitionary to right size
            lens_parameters = trim_dictionary(lens_parameters, size)

            # Sample the axis rotation angle
            lens_parameters["phi"] = self.sample_axis_rotation_angle(size=size)

            # Transform the axis ratio and the angle, to ellipticities e1, e2, using lenstronomy
            lens_parameters["e1"], lens_parameters["e2"] = phi_q2_ellipticity(
                lens_parameters["phi"], lens_parameters["q"]
            )

            # Sample shears
            lens_parameters["gamma1"], lens_parameters["gamma2"] = self.sample_shear(
                size=size)

            # Sample the spectral index of the mass density distribution
            lens_parameters["gamma"] = self.sample_mass_density_spectral_index(
                size=size)

            # sample gravitional waves source parameter
            param = dict(zs=np.array(zs))
            if samplers_params["source_parameters"]:
                param.update(self.sample_gw_parameters(size=size))
            gw_param = self.sample_source_parameters(size=size, param=param)

            # Add source params strongly lensed to the lens params
            lens_parameters.update(gw_param)

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

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> lens.strongly_lensed_source_redshifts(size=1000)
        """

        z_max = self.z_max

        def zs_function(zs_sl):
            # get zs
            # self.sample_source_redshifts from CBCSourceRedshiftDistribution class
            zs = self.sample_zs(size)  # this function is from CompactBinaryPopulation class
            # put strong lensing condition with optical depth
            tau = self.strong_lensing_optical_depth(zs)
            tau_max = self.strong_lensing_optical_depth(np.array([z_max]))[0] # tau increases with z
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

    def source_parameters(self, size, get_attribute=False, param=None):
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

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> lens.source_parameters(size=1000)
        """

        if get_attribute:
            return lambda size: self.sample_gw_parameters(size=size, param=param)
        else:
            # sample gravitional waves source parameter
            return self.sample_gw_parameters(size=size, param=param)
    
    def lens_redshift_SDSS_catalogue(self, zs, get_attribute=False, param=None):
        """
        Function to sample lens redshifts, conditioned on the lens being strongly lensed

        Parameters
        ----------
        zs : `float`
            source redshifts
        get_attribute : `bool`
            If True, returns a function that can be called with zs as input

        Returns
        -------
        zl : `float`
            lens redshifts

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> lens.lens_redshift_SDSS_catalogue(zs=1.0)
        """

        splineDc = self.splineDc  # spline coefficients for the comoving distance and redshifts
        splineDcInv = self.splineDcInv  # spline coefficients for the redshifts and comoving distance
        u = np.linspace(0, 1, 500)
        cdf = (10 * u**3 - 15 * u**4 + 6 * u**5)  # See the integral of Eq. A7 of https://arxiv.org/pdf/1807.07062.pdf (cdf)
        zs = np.array([zs]).reshape(-1)

        # lens redshifts
        #return self.Dc_to_z(lens_galaxy_Dc)
        if get_attribute:
            return njit(lambda zs: lens_redshift_SDSS_catalogue(zs, splineDc, splineDcInv, u, cdf))
        else:
            return lens_redshift_SDSS_catalogue(zs, splineDc, splineDcInv, u, cdf)

    def axis_rotation_angle_uniform(
        self, size=1000, phi_min=0.0, phi_max=2 * np.pi, get_attribute=False,param=None
    ):
        """
        Function to sample the axis rotation angle of the elliptical lens galaxy from a uniform distribution.

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        phi_min : `float`
            minimum axis rotation angle of the elliptical lens galaxy
        phi_max : `float`
            maximum axis rotation angle of the elliptical lens galaxy
        get_attribute : `bool`
            If True, returns a function that can be called with size as input
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. param = dict(phi_min=0.0, phi_max=2 * np.pi)

        Returns
        -------
        phi : `float`
            axis rotation angle of the elliptical lens galaxy

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> lens.axis_rotation_angle_uniform(size=1000)
        """

        if param:
            phi_min = param["phi_min"]
            phi_max = param["phi_max"]

        if get_attribute:
            return njit(lambda size: np.random.uniform(phi_min, phi_max, size=size))
        else:
            # Draw the angles from a uniform distribution
            return np.random.uniform(phi_min, phi_max, size=size)

    def shear_norm(self, size, scale=0.05, get_attribute=False, param=None):
        """
        Function to sample the elliptical lens galaxy shear from a normal distribution

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        scale : `float`
            standard deviation of the normal distribution
        get_attribute : `bool`
            If True, returns a function that can be called with size as input
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. param = dict(scale=0.05)

        Returns
        -------
        gamma_1 : `float`
            shear component in the x-direction
        gamma_2 : `float`
            shear component in the y-direction

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> gamma_1, gamma_2 = lens.shear_norm(size=1000)
        """

        if param:
            scale = param["scale"]

        if get_attribute:
            return njit(lambda size: np.random.normal(loc=0, scale=scale,size=(2,size)))
        else:
            # Draw an external shear from a normal distribution
            return np.random.normal(loc=0., scale=scale,size=(2,size))

    def mass_density_spectral_index_normal(
        self, size=1000, mean=2.0, std=0.2, get_attribute=False, param=None
    ):
        """
        Function to sample the lens galaxy spectral index of the mass density profile from a normal distribution

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample
        mean : `float`
            mean of the normal distribution
        std : `float`
            standard deviation of the normal distribution
        param : `dict`
            Allows to pass in above parameters as dict.
            e.g. param = dict(mean=2.0, std=0.2)

        Returns
        -------
        gamma : `float`
            spectral index of the density profile

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> lens.mass_density_spectral_index_normal(size=1000)
        """

        if param:
            mean = param["mean"]
            std = param["std"]

        if get_attribute:
            return njit(lambda size: np.random.normal(loc=mean, scale=std, size=size))
        else:
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

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> sigma = 200.0
        >>> zl = 0.5
        >>> zs = 1.0
        >>> lens.compute_einstein_radii(sigma, zl, zs)
        """

        # Compute the angular diameter distances
        Ds = self.angular_diameter_distance(zs)
        Dls = self.angular_diameter_distance_z1z2(zl, zs)
        # Compute the Einstein radii
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / (Ds)
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        return theta_E

    def rjs_with_cross_section_SIS(self, param_dict):
        """
        Function to conduct rejection sampling wrt einstein radius

        Parameters
        ----------
        param_dict : `dict`
            dictionary of lens parameters and source parameters

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

    def rjs_with_cross_section_SIE(self, param_dict):
        """
        Function to conduct rejection sampling wrt cross_section

        Parameters
        ----------
        param_dict : `dict`
            dictionary of lens parameters and source parameters

        Returns
        -------
        lens_params : `dict`
            dictionary of lens parameters after rejection sampling
        """

        theta_E = param_dict["theta_E"]
        q = param_dict["q"]
        phi_cut = phi_cut_SIE(q)
        size = len(theta_E)
        cross_section = theta_E**2 * phi_cut
        max_ = np.max(cross_section)  # maximum einstein radius
        u = np.random.uniform(0, max_, size=size)
        mask = u < cross_section

        # return the dictionary with the mask applied
        return {key: val[mask] for key, val in param_dict.items()}

    @property
    def sample_source_redshift_sl(self):
        """
        Function to sample source redshifts conditioned on the source being strongly lensed

        Parameters
        ----------
        size : `int`
            number samples to draw

        Returns
        -------
        zs : `numpy.ndarray` (1D array of floats)
            source redshifts conditioned on the source being strongly lensed

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> lens.sample_source_redshift_sl(size=1000)
        """

        return self._sample_source_redshift_sl

    @sample_source_redshift_sl.setter
    def sample_source_redshift_sl(self, prior):
        try:
            self._sample_source_redshift_sl = getattr(self, prior)
        except:
            self._sample_source_redshift_sl = prior

    @property
    def sample_source_parameters(self):
        """
        Function to sample source parameters conditioned on the source being strongly lensed

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample

        Returns
        -------
        source_parameters : `dict`
            dictionary of source parameters conditioned on the source being strongly lensed
        """
        return self._sample_source_parameters
    
    @sample_source_parameters.setter
    def sample_source_parameters(self, prior):
        try:
            args = self.lens_param_samplers_params["source_parameters"]
            self._sample_source_parameters = getattr(self, prior)(size=None, get_attribute=True, param=args)
        except:
            self._sample_source_parameters = prior

    @property
    def sample_lens_redshift(self):
        """
        Function to sample lens redshifts, conditioned on the lens being strongly lensed

        Parameters
        ----------
        zs : `numpy.ndarray` (1D array of floats)
            source redshifts

        Returns
        -------
        zl : `numpy.ndarray` (1D array of floats)
            lens redshifts corresponding to the source redshifts

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> zs = lens.sample_source_redshift_sl(size=1000)
        >>> lens.sample_lens_redshift(zs=zs)
        """

        return self._sample_lens_redshift

    @sample_lens_redshift.setter
    def sample_lens_redshift(self, prior):
        try:
            args = self.lens_param_samplers_params["lens_redshift"]
            self._sample_lens_redshift = getattr(self, prior)(zs=None, get_attribute=True, param=args)
        except:
            self._sample_lens_redshift = prior

    @property
    def sample_axis_rotation_angle(self):
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

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> lens.sample_axis_rotation_angle(size=1000)
        """
        return self._sample_axis_rotation_angle

    @sample_axis_rotation_angle.setter
    def sample_axis_rotation_angle(self, prior):
        try:
            args = self.lens_param_samplers_params["axis_rotation_angle"]
            self._sample_axis_rotation_angle = getattr(self, prior)(size=None, get_attribute=True, param=args)
        except:
            self._sample_axis_rotation_angle = prior

    @property
    def sample_shear(self):
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

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> gamma_1, gamma_2 = lens.shear_norm(size=1000)
        """

        return self._sample_shear

    @sample_shear.setter
    def sample_shear(self, prior):
        try:
            args = self.lens_param_samplers_params["shear"]
            self._sample_shear = getattr(self, prior)(size=None, get_attribute=True, param=args)
        except:
            self._sample_shear = prior

    @property
    def sample_mass_density_spectral_index(self):
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

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> lens.mass_density_spectral_index_normal(size=1000)
        """

        return self._sample_mass_density_spectral_index

    @sample_mass_density_spectral_index.setter
    def sample_mass_density_spectral_index(self, prior):
        try:
            args = self.lens_param_samplers_params["mass_density_spectral_index"]
            self._sample_mass_density_spectral_index = getattr(self, prior)(size=None, get_attribute=True, param=args)
        except:
            self._sample_mass_density_spectral_index = prior

    @property
    def sample_source_parameters(self):
        """
        Function to sample source parameters conditioned on the source being strongly lensed

        Parameters
        ----------
        size : `int`
            number of lens parameters to sample

        Returns
        -------
        source_parameters : `dict`
            dictionary of source parameters conditioned on the source being strongly lensed
        """
        return self._sample_source_parameters

    @sample_source_parameters.setter
    def sample_source_parameters(self, prior):
        try:
            self._sample_source_parameters = getattr(self, prior)
        except:
            self._sample_source_parameters = prior


    @property
    def available_lens_prior_list_and_its_params(self):
        """
        Dictionary with list all the available priors and it's corresponding parameters. This is an immutable instance attribute.
        """
        
        self._available_lens_prior_list_and_its_params = dict(
            source_redshift_sl=dict(strongly_lensed_source_redshifts=None),
            lens_redshift=dict(lens_redshift_SDSS_catalogue=None),
            velocity_dispersion=self.available_velocity_dispersion_list_and_its_params,
            axis_ratio=self.available_axis_ratio_list_and_its_params,
            axis_rotation_angle=dict(axis_rotation_angle_uniform=dict(phi_min=0.0, phi_max=2 * np.pi)),
            shear=dict(shear_norm=dict(scale=0.05)),
            mass_density_spectral_index=dict(mass_density_spectral_index_normal=dict(mean=2.0, std=0.2)),
            source_parameters=dict(sample_gw_parameters=None),
        )

        return self._available_lens_prior_list_and_its_params
    
    @property
    def available_lens_functions(self):
        """
        Dictionary with list all the available lens functions. This is an immutable instance attribute.
        """

        self._available_lens_functions = dict(
            strong_lensing_condition=["rjs_with_cross_section_SIE", "rjs_with_cross_section_SIE"],
            optical_depth=["SIS", "optical_depth_SIS_haris","optical_depth_SIS_hemanta", "SIE", "optical_depth_SIE_hemanta"],
            param_sampler_type=["sample_all_routine"],
        )

        return self._available_lens_functions


