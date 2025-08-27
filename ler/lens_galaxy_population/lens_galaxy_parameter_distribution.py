# -*- coding: utf-8 -*-
"""
This module contains the LensGalaxyPopulation class, which is used to sample lens galaxy parameters, source parameters conditioned on the source being strongly lensed. \n
The class inherits from the ImageProperties class, which is used calculate image properties (magnification, timedelays, source position, image position, morse phase). \n
Either the class takes in initialized CBCSourceParameterDistribution class as input or inherits the CBCSourceParameterDistribution class with default params (if no input) \n
"""

import warnings
warnings.filterwarnings("ignore")
import numpy as np
from scipy.integrate import quad
# from lenstronomy.Util.param_util import phi_q2_ellipticity
from multiprocessing import Pool
from tqdm import tqdm

# the following .py file will be called if they are not given in the class initialization
from ..gw_source_population import CBCSourceParameterDistribution
from .optical_depth import OpticalDepth
from ..image_properties import ImageProperties
from ..utils import add_dictionaries_together, trim_dictionary, FunctionConditioning, interpolator_pickle_path
from .jit_functions import phi_cut_SIE, phi_q2_ellipticity_hemanta
from .mp import cross_section_mp, rjs_sie_mp

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
            default: 'epl_shear_galaxy'
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
        |:meth:`~sample_all_routine_sie`      | Function to sample galaxy lens   |
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
        |:meth:`~density_profile_slope_normal`                             |
        +-------------------------------------+----------------------------------+
        |                                     | Function to sample the lens      |
        |                                     | galaxy spectral index of the     |
        |                                     | mass density profile from a      |
        |                                     | normal distribution              |
        +-------------------------------------+----------------------------------+
        |:meth:`~compute_einstein_radii`      | Function to compute the Einstein |
        |                                     | radii of the lens galaxies       |
        +-------------------------------------+----------------------------------+
        |:meth:`~rjs_with_cross_section_sis`  | Function to conduct rejection    |
        |                                     | sampling wrt einstein radius     |
        +-------------------------------------+----------------------------------+
        |:meth:`~rjs_with_cross_section_sie`  | Function to conduct rejection    |
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
        |:attr:`~sample_density_profile_slope`                             |
        +-------------------------------------+----------------------------------+
        |                                     | Function to sample the lens      |
        |                                     | galaxy spectral index of the     |
        |                                     | mass density profile from a      |
        |                                     | normal distribution              |
        +-------------------------------------+----------------------------------+
    """

    # Attributes
    cbc_pop = None
    """:class:`~CBCSourceParameterDistribution` class\n
    This is an already initialized class that contains a function (CBCSourceParameterDistribution.sample_gw_parameters) that actually samples the source parameters. 
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
        lens_type="epl_shear_galaxy",
        lens_functions= None,
        lens_functions_params=None,
        lens_param_samplers=None,
        lens_param_samplers_params=None,
        directory="./interpolator_pickle",
        create_new_interpolator=False,
        buffer_size=1000,
        **kwargs  # for initialization of CBCSourceParameterDistribution class and ImageProperties class
    ):
        print("\nInitializing LensGalaxyParameterDistribution class...\n")          
        self.event_type = event_type  # needed for the source population
        self.buffer_size = buffer_size  # buffer size for sampling lens parameters        

        # initializing parent classes
        self.class_initialization_lens(
            npool,
            z_min,
            z_max,
            cosmology,
            lens_type,
            lens_functions,
            lens_functions_params,
            lens_param_samplers, 
            lens_param_samplers_params, 
            directory,
            create_new_interpolator,
            params=kwargs,  # related parameters for CBCSourceParameterDistribution, ImageProperties classes
        );

        # function to sample source redshifts conditioned on the source being strongly lensed 
        self.sample_source_redshift_sl = getattr(self,self.lens_param_samplers["source_redshift_sl"])

        # interpolated cross section function is not very accurate for rejection sampling
        if self.lens_functions["cross_section"] == "interpolated_cross_section_function":
            self.cross_section = self.cross_section_caustic_area
        # function to sample lens parameters
        self.sample_lens_parameters_routine = getattr(self, self.lens_functions['param_sampler_type'])  
        # function to rejection sample wrt lens cross section
        self.rejection_sample_sl = getattr(self, self.lens_functions['strong_lensing_condition']) 

        # To find the normalization constant of the pdf p(z)
        # this under the assumption that the event is strongly lensed
        # Define the merger-rate density function
        pdf_unnormalized_ = lambda z: self.merger_rate_density_detector_frame(np.array([z])) * self.optical_depth.function(np.array([z]))
        pdf_unnormalized = lambda z: pdf_unnormalized_(z)[0]

        self.normalization_pdf_z_lensed = quad(
            pdf_unnormalized,
            self.z_min,
            self.z_max
        )[0]

    def class_initialization_lens(self, npool, z_min, z_max, cosmology, lens_type, lens_functions, lens_functions_params, lens_param_samplers,  lens_param_samplers_params,  directory, create_new_interpolator, params):
        """
            Initialize the LensGalaxyParameterDistribution class.

            Parameters
            ----------
            npool : `int`
                number of processors to use for sampling
            z_min : `float`
                minimum redshift of the lens galaxy
            z_max : `float`
                maximum redshift of the lens galaxy
            cosmology : `astropy.cosmology`
                cosmology object
            lens_type : `str`
                type of the lens galaxy
            lens_functions : `dict`
                dictionary with the lens related functions
            lens_functions_params : `dict`
                dictionary with the parameters for the lens related functions
            lens_param_samplers : `dict`
                dictionary with the priors for the sampler
            lens_param_samplers_params : `dict`
                dictionary with the parameters for the priors of the sampler
            directory : `str`
                directory where the interpolators are saved
            create_new_interpolator : `bool`
                if True, creates a new interpolator
            params : `dict`
                additional parameters for the CBCSourceParameterDistribution and ImageProperties classes

        """

        # initialize the optical depth class
        # this also initializes the lens related parameter samplers and functions:
        # self.lens_param_samplers, self.lens_param_samplers_params, self.lens_functions, self.lens_functions_params
        OpticalDepth.__init__(
            self,
            npool = npool,
            z_min = z_min,
            z_max = z_max,
            cosmology=cosmology,
            lens_type = lens_type,
            lens_functions = lens_functions,
            lens_functions_params = lens_functions_params,
            lens_param_samplers=lens_param_samplers,
            lens_param_samplers_params=lens_param_samplers_params,
            directory=directory,
            create_new_interpolator=create_new_interpolator,
        )

        # initialization of CBCSourceParameterDistribution class
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

        # initialize the image properties class
        input_params_image = dict(
            n_min_images=2,
            n_max_images=4,
            time_window=365*24*3600*20,
            # geocent_time_min=1126259462.4,
            # geocent_time_max=1126259462.4+365*24*3600*2,
            lens_model_list=["EPL_NUMBA", "SHEAR"],
        )
        input_params_image.update(params)

        # print("input_params_image", input_params_image)
        ImageProperties.__init__(
            self,
            npool=self.npool,
            z_min=self.z_min,
            z_max=self.z_max,
            n_min_images=input_params_image["n_min_images"],
            n_max_images=input_params_image["n_max_images"],
            lens_model_list=input_params_image["lens_model_list"],
            cosmology=self.cosmo,
            time_window=input_params_image["time_window"],
            # geocent_time_min=input_params_image["geocent_time_min"],
            # geocent_time_max=input_params_image["geocent_time_max"],
            spin_zero=input_params["spin_zero"],
            spin_precession=input_params["spin_precession"],
            directory=self.directory,
            create_new_interpolator=self.create_new_interpolator,
        )

    def sample_lens_parameters(self, size=1000,):
        """
            Function to sample galaxy lens parameters along with the source parameters, conditioned on the source being strongly lensed.

            Parameters
            ----------
            size : `int`
                number of lens parameters to sample

            Returns
            -------
            lens_parameters : `dict`
                dictionary of sampled lens parameters and source parameters. \n
                keys: \n
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
                gamma: density profile slope distribution \n
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
            >>> od = LensGalaxyParameterDistribution(lens_param_samplers=dict(velocity_dispersion="velocity_dispersion_ewoud"))
            >>> print(od.sample_lens_parameters(size=10))
        """

        print(f"sampling lens parameters with {self.lens_functions['param_sampler_type']}...")

        # get sample lens parameters, strongly lensed
        lens_parameters = self.sample_lens_parameters_routine(
            size=size
        )

        # sample gravitional waves source parameter
        param = dict(zs=lens_parameters["zs"])
        # zs won't be sampling again
        gw_param = self.sample_gw_parameters(size=size, param=param)
        # Add source params strongly lensed to the lens params
        lens_parameters.update(gw_param)

        return lens_parameters
    
    def sample_all_routine_sis_sl(self, size=1000):
        """
        Function to sample galaxy lens parameters. SIS cross section is used for rejection sampling.
        """

        buffer_size = self.buffer_size

        # Sample source redshifts from the source population
        # rejection sampled with optical depth
        zs = self.sample_source_redshift_sl(size)

        # # Sample lens redshifts
        zl = self.lens_redshift.rvs(size, zs)

        sigma = np.zeros(size)
        theta_E = np.zeros(size)
        
        sigma_max = self.velocity_dispersion.info['sigma_max']

        for i in tqdm(range(size), ncols=100, disable=False):
            zs_ = zs[i]*np.ones(buffer_size)
            zl_ = zl[i]*np.ones(buffer_size)

            # cross_section_max calculation
            theta_E_max = self.compute_einstein_radii(np.array([sigma_max]), np.array([zl[i]]), np.array([zs[i]]))[0]
            cross_section_max = np.pi*theta_E_max**2

            while True:
                # Create a dictionary of the lens parameters; sigma, theta_E, q, phi, e1, e2
                lens_parameters_ = self.sampling_routine_sie_nsl(zl_, zs_, size=buffer_size)

                # Rejection sample based on the lensing probability, that is, rejection sample wrt theta_E
                lens_parameters_, mask, cross_section_max_ = self.rjs_with_cross_section_sis(
                    lens_parameters_, cross_section_max
                )  # proportional to pi theta_E^2

                if cross_section_max_>cross_section_max:
                    cross_section_max = cross_section_max_

                if np.sum(mask) > 0:
                    break
            
            sigma[i] = lens_parameters_["sigma"][0]
            theta_E[i] = lens_parameters_["theta_E"][0]
            
        # sample additional lens parameters
        # P(q|SL), P(gamma|SL), P(gamma1, gamma2|SL)
        q = self.axis_ratio.rvs(size)
        gamma = self.density_profile_slope_sl.rvs(size)
        gamma1, gamma2 = self.external_shear_sl.rvs(size)

        phi = self.axis_rotation_angle.rvs(size)
        e1, e2 = phi_q2_ellipticity_hemanta(phi, q)

        # Create a dictionary of the lens parameters
        lens_parameters = {
            "zl": zl,
            "zs": zs,
            "sigma": sigma,
            "theta_E": theta_E,
            "q": q,
            "phi": phi,
            "e1": e1,
            "e2": e2,
            "gamma": gamma,
            "gamma1": gamma1,
            "gamma2": gamma2,
        }

        return lens_parameters

    def sample_all_routine_sie_sl(self, size=1000):
        """
        Function to sample galaxy lens parameters. SIE cross section is used for rejection sampling.

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
            gamma: density profile slope distribution \n

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> lens.sample_all_routine_sie(size=1000)
        """

        buffer_size = self.buffer_size

        # Sample source redshifts from the source population
        # rejection sampled with optical depth
        zs = self.sample_source_redshift_sl(size)

        # # Sample lens redshifts
        zl = self.lens_redshift.rvs(size, zs)

        sigma = np.zeros(size)
        theta_E = np.zeros(size)
        q = np.zeros(size)

        # # set-up multiprocessing args
        # q_args_cdf_values = self.axis_ratio.cdf_values
        # q_args_x_array = self.axis_ratio.x_array
        # q_args_conditioned_y_array = self.axis_ratio.conditioned_y_array

        # sigma_args_cdf_values = self.velocity_dispersion.cdf_values
        # sigma_args_x_array = self.velocity_dispersion.x_array
        # sigma_args_conditioned_y_array = self.velocity_dispersion.conditioned_y_array

        # da_args_function_spline = self.angular_diameter_distance.function_spline
        # da_args_x_array = self.angular_diameter_distance.x_array

        # idx = np.arange(size)
        # input_params = np.array([(idx[i], zs[i], zl[i], sigma_args_cdf_values, sigma_args_x_array, sigma_args_conditioned_y_array, q_args_cdf_values, q_args_x_array, q_args_conditioned_y_array, da_args_function_spline, da_args_x_array, buffer_size) for i in range(size)], dtype=object)

        # with Pool(processes=self.npool) as pool:
        #     for result in tqdm(
        #         pool.imap_unordered(rjs_sie_mp, input_params),
        #         total=size,
        #         ncols=100,
        #         disable=False,
        #     ):
        #         (
        #             idx_,
        #             sigma_,
        #             theta_E_,
        #             q_,
        #         ) = result

        #         sigma[idx_] = sigma_
        #         theta_E[idx_] = theta_E_
        #         q[idx_] = q_
        # phi = self.axis_rotation_angle.rvs(size)
        # e1, e2 = phi_q2_ellipticity_hemanta(phi, q)
        
        sigma_max = self.velocity_dispersion.info['sigma_max']

        for i in tqdm(range(size), ncols=100, disable=False):
            zs_ = zs[i]*np.ones(buffer_size)
            zl_ = zl[i]*np.ones(buffer_size)

            # cross_section_max calculation
            theta_E_max = self.compute_einstein_radii(np.array([sigma_max]), np.array([zl[i]]), np.array([zs[i]]))[0]
            sie_factor=1.
            cross_section_max = sie_factor*np.pi*theta_E_max**2

            while True:
                # Create a dictionary of the lens parameters; sigma, theta_E, q, phi, e1, e2
                lens_parameters_ = self.sampling_routine_sie_nsl(zl_, zs_, size=buffer_size)

                # Rejection sample based on the lensing probability, that is, rejection sample wrt theta_E
                lens_parameters_, mask, cross_section_max_ = self.rjs_with_cross_section_sie_feixu(
                    lens_parameters_, cross_section_max
                )  # proportional to pi theta_E^2

                if cross_section_max_>cross_section_max:
                    cross_section_max = cross_section_max_

                if np.sum(mask) > 0:
                    break
            
            sigma[i] = lens_parameters_["sigma"][0]
            theta_E[i] = lens_parameters_["theta_E"][0]
            q[i] = lens_parameters_["q"][0]
            
        # sample additional lens parameters
        # P(gamma|SL), P(gamma1, gamma2|SL)
        gamma = self.density_profile_slope_sl.rvs(size)
        gamma1, gamma2 = self.external_shear_sl.rvs(size)

        phi = self.axis_rotation_angle.rvs(size)
        e1, e2 = phi_q2_ellipticity_hemanta(phi, q)

        # Create a dictionary of the lens parameters
        lens_parameters = {
            "zl": zl,
            "zs": zs,
            "sigma": sigma,
            "theta_E": theta_E,
            "q": q,
            "phi": phi,
            "e1": e1,
            "e2": e2,
            "gamma": gamma,
            "gamma1": gamma1,
            "gamma2": gamma2,
        }

        return lens_parameters
        
    def sample_all_routine_epl_shear_sl(self, size=1000):
        """
            Function to sample galaxy lens parameters along. EPL shear cross section is used for rejection sampling.

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
                gamma: density profile slope distribution \n

            Examples
            --------
            >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
            >>> lens = LensGalaxyParameterDistribution()
            >>> lens.sample_all_routine_sie(size=1000)
        """

        buffer_size = self.buffer_size

        # Sample source redshifts from the source population
        # rejection sampled with optical depth
        zs = self.sample_source_redshift_sl(size)

        # # Sample lens redshifts
        zl = self.lens_redshift.rvs(size, zs)

        sigma = np.zeros(size)
        theta_E = np.zeros(size)
        q = np.zeros(size)
        phi = np.zeros(size)
        gamma1 = np.zeros(size)
        gamma2 = np.zeros(size)
        gamma = np.zeros(size)

        sigma_max = self.velocity_dispersion.info['sigma_max']

        # this will take some time
        for i in tqdm(range(size), ncols=100, disable=False):
            zs_ = zs[i]*np.ones(buffer_size)
            zl_ = zl[i]*np.ones(buffer_size)

            # cross_section_max calculation            
            theta_E_max = self.compute_einstein_radii(np.array([sigma_max]), np.array([zl[i]]), np.array([zs[i]]))[0]
            epl_factor=5.
            cross_section_max = epl_factor*np.pi*theta_E_max**2

            while True:
                # Create a dictionary of the lens parameters, gamma, gamma1, gamma2, sigma, theta_E, q, phi, e1, e2
                lens_parameters_ = self.sampling_routine_epl_shear_nsl(zl_, zs_, size=buffer_size)

                # Rejection sample based on the lensing probability, that is, rejection sample wrt theta_E
                lens_parameters_, mask, cross_section_max_ = self.rejection_sample_sl(
                    lens_parameters_, cross_section_max
                )  # proportional to pi theta_E^2
                
                if cross_section_max_>cross_section_max:
                    cross_section_max = cross_section_max_

                # check the size of the lens parameters
                if np.sum(mask) > 0:
                    break

            sigma[i] = lens_parameters_["sigma"][0]
            theta_E[i] = lens_parameters_["theta_E"][0]
            q[i] = lens_parameters_["q"][0]
            phi[i] = lens_parameters_["phi"][0]
            gamma[i] = lens_parameters_["gamma"][0]
            gamma1[i] = lens_parameters_["gamma1"][0]
            gamma2[i] = lens_parameters_["gamma2"][0]
        e1, e2 = phi_q2_ellipticity_hemanta(phi, q)

        # Create a dictionary of the lens parameters
        lens_parameters = {
            "zl": zl,
            "zs": zs,
            "sigma": sigma,
            "theta_E": theta_E,
            "q": q,
            "phi": phi,
            "e1": e1,
            "e2": e2,
            "gamma": gamma,
            "gamma1": gamma1,
            "gamma2": gamma2,
        }

        return lens_parameters
        
    def sampling_routine_sis_nsl(self, zl, zs, size=1000):
        """
        Function to sample SIS lens related parameters.

        Parameters
        ----------
        zl : `float`
            lens redshifts
        zs : `float`
            source redshifts
        size : `int`
            number of lens parameters to sample

        Returns
        -------
        lens_parameters : `dict`
            dictionary of sampled lens parameters.
            keys: sigma, theta_E
        """
        try:
            sigma = self.velocity_dispersion.rvs(size, zl)
        except:
            sigma = self.velocity_dispersion.rvs(size)

        # Compute the Einstein radii
        theta_E = self.compute_einstein_radii(sigma, zl, zs)

        # Create a dictionary of the lens parameters
        lens_parameters = {
            "zs": zs,
            "zl": zl,
            "sigma": sigma,
            "theta_E": theta_E,
        }

        return lens_parameters
    
    def sampling_routine_sie_nsl(self, zl, zs, size=1000):
        """
        Function to sample SIE lens related parameters.

        Parameters
        ----------
        zl : `float`
            lens redshifts
        zs : `float`
            source redshifts
        size : `int`
            number of lens parameters to sample

        Returns
        -------
        lens_parameters : `dict`
            dictionary of sampled lens parameters.
            keys: sigma, q, phi
        """

        lens_parameters = self.sampling_routine_sis_nsl(zl, zs, size=size)

        # Sample axis ratios
        try:
            q = self.axis_ratio.rvs(size, lens_parameters["sigma"])
        except:
            q = self.axis_ratio.rvs(size)

        # Sample the axis rotation angle
        phi = self.axis_rotation_angle.rvs(size)

        # Transform the axis ratio and the angle, to ellipticities e1, e2, using lenstronomy
        # e1, e2 = phi_q2_ellipticity_hemanta(phi, q)

        # Add the lensing parameter dictionaries together
        lens_parameters["q"] = q
        lens_parameters["phi"] = phi
        # lens_parameters["e1"] = e1
        # lens_parameters["e2"] = e2

        return lens_parameters
    
    def sampling_routine_epl_shear_nsl(self, zl, zs, size=1000):
        """
        Function to sample EPL and shear related parameters.

        Parameters
        ----------
        zl : `float`
            lens redshifts
        zs : `float`
            source redshifts
        size : `int`
            number of lens parameters to sample

        Returns
        -------
        lens_parameters : `dict`
            dictionary of sampled lens parameters.
            keys: sigma, q, phi, gamma, gamma1, gamma2
        """

        lens_parameters = self.sampling_routine_sie_nsl(zl, zs, size=size)

        # Sample the density profile slope distribution
        gamma = self.density_profile_slope(size)

        # Sample shears
        gamma1, gamma2 = self.external_shear.rvs(size)

        # Add the lensing parameter dictionaries together
        lens_parameters["gamma"] = gamma
        lens_parameters["gamma1"] = gamma1
        lens_parameters["gamma2"] = gamma2

        return lens_parameters

    def strongly_lensed_source_redshifts(self, size=1000):
        """
        Function to sample source redshifts, conditioned on the source being strongly lensed.

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
            zs = self.zs(size)  # this function is from CBCSourceParameterDistribution class
            # put strong lensing condition with optical depth
            tau = self.optical_depth(zs)
            tau_max = self.optical_depth(np.array([z_max]))[0] # tau increases with z
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
        
    def rjs_with_cross_section_sis(self, param_dict, cross_section_max=0.):
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
        cross_section = np.pi * theta_E**2

        max_ = np.max(cross_section)
        if cross_section_max>max_:
            max_ = cross_section_max

        mask = np.random.uniform(size=size) < (cross_section / max_)

        # return the dictionary with the mask applied
        dict_ = {key: val[mask] for key, val in param_dict.items()}
        dict_["cross_section"] = cross_section[mask]
        return dict_, mask, max_

    def rjs_with_cross_section_sie_feixu(self, param_dict, cross_section_max=0.):
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

        #print("rjs_with_cross_section_sie_feixu")
        theta_E = param_dict["theta_E"]
        q = param_dict["q"]
        size = len(theta_E)
        cross_section = np.pi * theta_E**2 * phi_cut_SIE(q)

        max_ = np.max(cross_section)
        if cross_section_max>max_:
            max_ = cross_section_max

        u = np.random.uniform(0, max_, size=size)
        mask = u < cross_section

        # return the dictionary with the mask applied
        dict_ = {key: val[mask] for key, val in param_dict.items()}
        dict_["cross_section"] = cross_section[mask]
        return dict_, mask, max_
    
    def rjs_with_cross_section(self, param_dict, cross_section_max=0.):
        """
            Function to conduct rejection sampling wrt cross_section of EPL+Shear lens

            Parameters
            ----------
            param_dict : `dict`
                dictionary of lens parameters and source parameters

            Returns
            -------
            lens_params : `dict`
                dictionary of lens parameters after rejection sampling
        """

        theta_E_cut = 2.9243287409459857e-08 # this is numerically found

        # Pre-filter param_dict directly
        size_original = len(param_dict["theta_E"])
        idx_ = param_dict["theta_E"] > theta_E_cut
        param_dict = {key: val[idx_]
                    for key, val in param_dict.items()}

        size = len(param_dict["theta_E"])  # Update size after filtering
        theta_E = param_dict["theta_E"]
        e1, e2 = phi_q2_ellipticity_hemanta(param_dict["phi"], param_dict["q"])
        gamma = param_dict["gamma"]
        gamma1 = param_dict["gamma1"]
        gamma2 = param_dict["gamma2"]

        cross_section = self.cross_section(theta_E, e1, e2, gamma, gamma1, gamma2)

        max_ = np.max(cross_section)
        if cross_section_max>max_:
            max_ = cross_section_max

        u = np.random.uniform(0, max_, size=size)
        mask = u < cross_section

        # return the dictionary with the mask applied
        dict_ = {key: val[mask] for key, val in param_dict.items()}
        dict_["cross_section"] = cross_section[mask]

        mask_complete = np.zeros(size_original, dtype=bool)
        mask_complete[idx_] = mask
        return dict_, mask_complete, max_

    
    def rjs_with_cross_section_mp(self, param_dict, cross_section_max=0.):
        """
            Function to conduct rejection sampling wrt cross_section, multiprocessing

            Parameters
            ----------
            param_dict : `dict`
                dictionary of lens parameters and source parameters

            Returns
            -------
            lens_params : `dict`
                dictionary of lens parameters after rejection sampling
        """

        theta_E_cut = 1e-09 # this is numerically found
        theta_E = param_dict["theta_E"]
        size_original = len(theta_E)

        # Pre-filter param_dict directly
        
        idx_ = theta_E > theta_E_cut
        if np.sum(idx_) > 0:
        
            param_dict = {key: val[idx_] 
                        for key, val in param_dict.items()}
            size = np.sum(idx_)

            e1, e2 = phi_q2_ellipticity_hemanta(param_dict["phi"], param_dict["q"])

            # Combine parameters into a single array for multiprocessing
            params = np.array([
                param_dict["theta_E"],
                e1,
                e2,
                param_dict["gamma"],
                param_dict["gamma1"],
                param_dict["gamma2"],
                np.arange(size, dtype=int)
            ]).T

            cross_section = np.zeros(size)  # Directly create filtered array
            # with Pool(processes=self.npool) as pool:
            #     for result in tqdm(pool.imap_unordered(cross_section_mp, params), total=size,
            #                     ncols=100):  # Remove total from tqdm
            #         idx_, tau_ = result
            #         cross_section[idx_] = tau_
            # without tqdm. Use map
            # Use multiprocessing to calculate cross-section values
            with Pool(processes=self.npool) as pool:
                for idx, tau in pool.imap_unordered(cross_section_mp, params):
                    cross_section[idx] = tau

            # Perform rejection sampling
            max_ = np.max(cross_section)

            if cross_section_max>max_:
                max_ = cross_section_max
            random_values = np.random.uniform(0, max_, size=size)
            mask = random_values < cross_section
        else:
            print("No valid values found")
            mask = np.zeros(size, dtype=bool)

        # Return the dictionary with the mask applied
        dict_ = {key: val[mask] for key, val in param_dict.items()}
        dict_["cross_section"] = cross_section[mask]

        mask_complete = np.zeros(size_original, dtype=bool)
        mask_complete[idx_] = mask
        return dict_, mask_complete, max_

    # @property
    # def sample_source_parameters(self):
    #     """
    #     Function to sample source parameters conditioned on the source being strongly lensed

    #     Parameters
    #     ----------
    #     size : `int`
    #         number of lens parameters to sample

    #     Returns
    #     -------
    #     source_parameters : `dict`
    #         dictionary of source parameters conditioned on the source being strongly lensed
    #     """
    #     return self._sample_source_parameters

    # @sample_source_parameters.setter
    # def sample_source_parameters(self, prior):
    #     try:
    #         self._sample_source_parameters = getattr(self, prior)
    #     except:
    #         self._sample_source_parameters = prior

