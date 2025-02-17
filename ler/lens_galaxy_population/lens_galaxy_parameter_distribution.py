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

# the following .py file will be called if they are not given in the class initialization
from ..gw_source_population import CBCSourceParameterDistribution
from .optical_depth import OpticalDepth
from ..image_properties import ImageProperties
from ..utils import add_dictionaries_together, trim_dictionary, FunctionConditioning, interpolator_pickle_path
from .jit_functions import phi_cut_SIE, phi_q2_ellipticity_hemanta
from .mp import cross_section_mp

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
        buffer_size=int(2e5),
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
        # create kde for strong lensing paramters: gamma, gamma1, gamma2
        # this is done to speed up the sampling lens parameters
        self.create_kde_sl();
        # function to sample lens parameters
        self.sample_lens_parameters_routine = getattr(self, self.lens_functions['param_sampler_type'])  
        # function to rejection sample wrt lens cross section
        self.rejection_sample_sl = getattr(self, self.lens_functions['strong_lensing_condition']) 
        

        # To find the normalization constant of the pdf p(z)
        # this under the assumption that the event is strongly lensed
        # Define the merger-rate density function
        pdf_unnormalized_ = lambda z: self.merger_rate_density_detector_frame(np.array([z]), param=self.merger_rate_density_param) * self.optical_depth.function(np.array([z]))
        pdf_unnormalized = lambda z: pdf_unnormalized_(z)[0]

        self.normalization_pdf_z_lensed = quad(
            pdf_unnormalized,
            self.z_min,
            self.z_max
        )[0]

    def create_kde_sl(self):
        """
            Creates kernel density estimates (KDE) for strong lensing parameters: gamma, gamma1, and gamma2.
            
            This function checks if the necessary interpolator files exist for the parameters 
            associated with strongly lensed gamma and gamma1, gamma2. If they do not exist, 
            or if a new creation is requested, it generates the interpolators by sampling 
            lens parameters under strong lensing conditions. The KDEs are computed to speed 
            up the sampling of these lens parameters.
            
            The function utilizes various class attributes to obtain parameter information 
            and buffer size. It temporarily modifies certain settings and reverts them back 
            after interpolator creation.

            Attributes
            ----------
            density_profile_slope_sl : callable
                A kernel density estimator for the density profile slope under strong lensing.
            external_shear_sl : callable
                A kernel density estimator for the external shear parameters under strong lensing.

            Raises
            ------
            FileNotFoundError
                If the necessary pickle files do not exist and interpolator creation fails.
        """

        self.cross_section = self.cross_section_caustic_area
        self.rejection_sample_sl = self.rjs_with_cross_section_mp

        param_dict_given_ = self.external_shear.info 
        param_dict_given_.update(self.density_profile_slope.info)
        param_dict_given_.update(self.axis_ratio.info)
        param_dict_given_.update(self.density_profile_slope.info)
        param_dict_given_.update(self.velocity_dispersion.info)
        param_dict_given_['buffer_size'] = self.buffer_size

        gamma_dict = param_dict_given_.copy()
        gamma_dict['name'] = 'density_profile_slope_sl'
        gamma12_dict = param_dict_given_.copy()
        gamma12_dict['name'] = 'external_shear_sl'

        # check first whether the directory, subdirectory and pickle exist
        _, it_exist_gamma = interpolator_pickle_path(
            param_dict_given=gamma_dict,
            directory=self.directory,
            sub_directory='density_profile_slope',
            interpolator_name=gamma_dict['name'],
        )
        _, it_exist_gamma12 = interpolator_pickle_path(
            param_dict_given=gamma12_dict,
            directory=self.directory,
            sub_directory='external_shear',
            interpolator_name=gamma12_dict['name'],
        )

        create_new = self.create_new_interpolator['lens_parameters_kde_sl']['create_new']
        if (it_exist_gamma is False) or (it_exist_gamma12 is False) or create_new:
            print("Creating interpolator for strongly lensed gamma and gamma1, gamma2...")
            # create the interpolator
            buffer_size_ = self.buffer_size
            self.buffer_size = int(2e5)
            # sample lens parameters, strongly lensed
            size = self.create_new_interpolator['lens_parameters_kde_sl']['resolution']
            print(f"Sampling lens parameters for interpolator creation with size: {size}")
            lens_params = self.sample_all_routine_epl_shear_sl(size)
            # revert back to the original settings
            self.buffer_size = buffer_size_

            gamma = lens_params['gamma']
            gamma1 = lens_params['gamma1']
            gamma2 = lens_params['gamma2']
        else:
            print("Loading interpolator for strongly lensed gamma and gamma1, gamma2...")
            gamma = None
            gamma1 = None
            gamma2 = None

        # create kde for gamma
        self.density_profile_slope_sl = getattr(self, 'density_profile_slope_sl_sampler')(size=None, get_attribute=True, gamma=gamma)
        self.external_shear_sl = getattr(self, 'external_shear_sl_sampler')(size=None, get_attribute=True, gamma1=gamma1, gamma2=gamma2)


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
        # this also initializes the lens related parameter samplers and functions
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
            geocent_time_min=1126259462.4,
            geocent_time_max=1126259462.4+365*24*3600*10,
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
            geocent_time_min=input_params_image["geocent_time_min"],
            geocent_time_max=input_params_image["geocent_time_max"],
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

        buffer_size = self.buffer_size if self.buffer_size is not None else size
        
        lens_parameters = dict()

        # Sample source redshifts from the source population
        # rejection sampled with optical depth
        zs = self.sample_source_redshift_sl(buffer_size)

        # # Sample lens redshifts
        zl = self.lens_redshift.rvs(buffer_size, zs)

        while True:

            # Create a dictionary of the lens parameters; q, phi, e1, e2
            lens_parameters_ = self.rjs_with_cross_section_sis(zl, zs, size=buffer_size)

            # Rejection sample based on the lensing probability, that is, rejection sample wrt theta_E
            lens_parameters_ = self.rejection_sample_sis(
                lens_parameters_
            )  # proportional to pi theta_E^2

            # Add the lensing parameter dictionaries together
            lens_parameters = add_dictionaries_together(
                lens_parameters, lens_parameters_
            )

            # check the size of the lens parameters
            size_now = len(lens_parameters["sigma"])
            #print(f"current sampled size: {size_now}")
            if size_now > size:
                break

        # trim to the right size
        lens_parameters = trim_dictionary(lens_parameters, size)
        lens_parameters["zs"] = zs[:size]
        lens_parameters["zl"] = zl[:size]
        lens_parameters["q"] = self.axis_ratio.rvs(size)
        lens_parameters["phi"] = self.axis_rotation_angle.rvs(size)
        
        # sample additional lens parameters
        if self.lens_type == "epl_shear_galaxy":
            lens_parameters["gamma"] = self.density_profile_slope_sl.rvs(size)
            lens_parameters["gamma1"], lens_parameters["gamma2"] = self.external_shear_sl.rvs(size)
        else:
            lens_parameters["gamma"] = self.density_profile_slope.rvs(size)
            lens_parameters["gamma1"], lens_parameters["gamma2"] = self.external_shear.rvs(size)

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

        buffer_size = self.buffer_size if self.buffer_size is not None else size
        
        lens_parameters = dict()

        # Sample source redshifts from the source population
        # rejection sampled with optical depth
        zs = self.sample_source_redshift_sl(buffer_size)

        # # Sample lens redshifts
        zl = self.lens_redshift.rvs(buffer_size, zs)

        while True:

            # Create a dictionary of the lens parameters; q, phi, e1, e2
            lens_parameters_ = self.sampling_routine_sie_nsl(zl, zs, size=buffer_size)

            # Rejection sample based on the lensing probability, that is, rejection sample wrt theta_E
            lens_parameters_ = self.rjs_with_cross_section_sie_feixu(
                lens_parameters_
            )  # proportional to pi theta_E^2

            # Add the lensing parameter dictionaries together
            lens_parameters = add_dictionaries_together(
                lens_parameters, lens_parameters_
            )

            # check the size of the lens parameters
            size_now = len(lens_parameters["sigma"])
            #print(f"current sampled size: {size_now}")
            if size_now > size:
                break

        # trim to the right size
        lens_parameters = trim_dictionary(lens_parameters, size)
        lens_parameters["zs"] = zs[:size]
        lens_parameters["zl"] = zl[:size]   
        
        # sample additional lens parameters
        if self.lens_type == "epl_shear_galaxy":
            lens_parameters["gamma"] = self.density_profile_slope_sl.rvs(size)
            lens_parameters["gamma1"], lens_parameters["gamma2"] = self.external_shear_sl.rvs(size)
        else:
            lens_parameters["gamma"] = self.density_profile_slope.rvs(size)
            lens_parameters["gamma1"], lens_parameters["gamma2"] = self.external_shear.rvs(size)

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

        buffer_size = self.buffer_size if self.buffer_size is not None else size
        
        lens_parameters = dict()

        # Sample source redshifts from the source population
        # rejection sampled with optical depth
        zs = self.sample_source_redshift_sl(buffer_size)

        # # Sample lens redshifts
        zl = self.lens_redshift.rvs(buffer_size, zs)

        while True:

            # Create a dictionary of the lens parameters
            lens_parameters_ = self.sampling_routine_epl_shear_nsl(zl, zs, size=buffer_size)

            # Rejection sample based on the lensing probability, that is, rejection sample wrt theta_E
            lens_parameters_ = self.rejection_sample_sl(
                lens_parameters_
            )  # proportional to pi theta_E^2

            # Add the lensing parameter dictionaries together
            lens_parameters = add_dictionaries_together(
                lens_parameters_, lens_parameters
            )

            # check the size of the lens parameters
            size_now = len(lens_parameters["sigma"])
            print(f"current sampled size: {size_now}")
            if size_now > size:
                break

        # Trim dicitionary to right size
        lens_parameters = trim_dictionary(lens_parameters, size)
        lens_parameters["zs"] = zs
        lens_parameters["zl"] = zl

        # trim to the right size
        for key in lens_parameters.keys():
            lens_parameters[key] = lens_parameters[key][:size]

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
            keys: sigma, q, phi, e1, e2
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
        e1, e2 = phi_q2_ellipticity_hemanta(phi, q)

        # Add the lensing parameter dictionaries together
        lens_parameters["q"] = q
        lens_parameters["phi"] = phi
        lens_parameters["e1"] = e1
        lens_parameters["e2"] = e2

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
            keys: sigma, q, phi, e1, e2, gamma, gamma1, gamma2
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
            zs = self.sample_zs(size)  # this function is from CBCSourceParameterDistribution class
            # put strong lensing condition with optical depth
            tau = self.optical_depth(zs)
            tau_max = self.optical_depth(z_max)[0] # tau increases with z
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

        
    def density_profile_slope_sl_sampler(self, size=1000, get_attribute=False, **kwargs):
        """
            Function to sample the lens galaxy density profile slope profile with strong condition applied.

            Parameters
            ----------
            size : `int`
                number of lens parameters to sample

            Returns
            -------
            gamma : `float`
                spectral index of the density profile

        """

        if 'gamma' in kwargs:
            gamma = kwargs['gamma']
        else:
            gamma = None

        param_dict_given_ = self.external_shear.info 
        param_dict_given_.update(self.density_profile_slope.info)
        param_dict_given_.update(self.axis_ratio.info)
        param_dict_given_.update(self.density_profile_slope.info)
        param_dict_given_.update(self.velocity_dispersion.info)
        param_dict_given_['buffer_size'] = self.buffer_size
        param_dict_given_['name'] = "density_profile_slope_sl"

        gamma_object = FunctionConditioning(
            x_array=gamma,
            y_array=None,
            gaussian_kde=True,
            create_rvs=True,
            create_pdf=True,
            callback='rvs',
            create_new=self.create_new_interpolator['lens_parameters_kde_sl']['create_new'],
            param_dict_given=param_dict_given_,
            directory=self.directory,
            sub_directory='density_profile_slope',
            name=param_dict_given_['name'],
        )

        if get_attribute:
            return gamma_object
        else:
            return gamma_object.rvs(size)
        
    def external_shear_sl_sampler(self, size=1000, get_attribute=False, **kwargs):
        """
            Function to sample the lens galaxy shear with strong condition applied.

            Parameters
            ----------
            size : `int`
                number of lens parameters to sample

            Returns
            -------
            gamma : `float`
                spectral index of the density profile

        """

        if ('gamma1' in kwargs) and ('gamma2' in kwargs):
            gamma1 = kwargs['gamma1']
            gamma2 = kwargs['gamma2']
        else:
            gamma1 = None
            gamma2 = None

        param_dict_given_ = self.external_shear.info 
        param_dict_given_.update(self.density_profile_slope.info)
        param_dict_given_.update(self.axis_ratio.info)
        param_dict_given_.update(self.density_profile_slope.info)
        param_dict_given_.update(self.velocity_dispersion.info)
        param_dict_given_['buffer_size'] = self.buffer_size
        param_dict_given_['name'] = "external_shear_sl"

        gamma12_object = FunctionConditioning(
            x_array=gamma1,
            y_array=gamma2,
            gaussian_kde=True,
            create_rvs=True,
            create_pdf=True,
            callback='rvs',
            create_new=self.create_new_interpolator['lens_parameters_kde_sl']['create_new'],
            param_dict_given=param_dict_given_,
            directory=self.directory,
            sub_directory="external_shear",
            name="external_shear_sl",
        )

        if get_attribute:
            return gamma12_object
        else:
            return gamma12_object.rvs(size)
        
    def rjs_with_cross_section_sis(self, param_dict):
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
        cross_section_max = np.max(cross_section)  # maximum einstein radius
        mask = np.random.uniform(size=size) < (cross_section / cross_section_max)

        # return the dictionary with the mask applied
        dict_ = {key: val[mask] for key, val in param_dict.items()}
        dict_["cross_section"] = cross_section[mask]
        return dict_

    def rjs_with_cross_section_sie_feixu(self, param_dict):
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
        max_ = np.max(cross_section)  # maximum einstein radius
        u = np.random.uniform(0, max_, size=size)
        mask = u < cross_section

        # return the dictionary with the mask applied
        dict_ = {key: val[mask] for key, val in param_dict.items()}
        dict_["cross_section"] = cross_section[mask]
        return dict_
    
    def rjs_with_cross_section(self, param_dict):
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

        theta_E_cut = 2.9243287409459857e-08

        # Pre-filter param_dict directly
        idx = param_dict["theta_E"] > theta_E_cut
        param_dict = {key: val[idx] 
                    for key, val in param_dict.items()}

        size = len(param_dict["theta_E"])  # Update size after filtering
        theta_E = param_dict["theta_E"]
        e1 = param_dict["e1"]
        e2 = param_dict["e2"]
        gamma = param_dict["gamma"]
        gamma1 = param_dict["gamma1"]
        gamma2 = param_dict["gamma2"]

        cross_section = self.cross_section_function(theta_E, e1, e2, gamma, gamma1, gamma2)

        max_ = np.max(cross_section)  # maximum einstein radius
        u = np.random.uniform(0, max_, size=size)
        mask = u < cross_section

        # return the dictionary with the mask applied
        dict_ = {key: val[mask] for key, val in param_dict.items()}
        dict_["cross_section"] = cross_section[mask]
        return dict_

    
    def rjs_with_cross_section_mp(self, param_dict):
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


        theta_E_cut = 2.9243287409459857e-08

        # Pre-filter param_dict directly
        idx = param_dict["theta_E"] > theta_E_cut
        param_dict = {key: val[idx] 
                    for key, val in param_dict.items()}

        # Update size after filtering
        size = len(param_dict["theta_E"])

        # Combine parameters into a single array for multiprocessing
        params = np.array([
            param_dict["theta_E"],
            param_dict["e1"],
            param_dict["e2"],
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
        max_cross_section = np.max(cross_section)
        random_values = np.random.uniform(0, max_cross_section, size=size)
        mask = random_values < cross_section

        # Return the dictionary with the mask applied
        dict_ = {key: val[mask] for key, val in param_dict.items()}
        dict_["cross_section"] = cross_section[mask]

        return dict_

    @property
    def density_profile_slope_sl(self):
        """
        Function to sample the lens galaxy density profile slope profile from a normal distribution

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
        >>> lens.density_profile_slope_normal(size=1000)
        """

        return self._sample_density_profile_slope

    @density_profile_slope_sl.setter
    def density_profile_slope_sl(self, prior):
        if isinstance(prior, str):
            args = self.lens_param_samplers_params["density_profile_slope"]
            self._sample_density_profile_slope = getattr(self, prior)(size=None, get_attribute=True, param=args)
        elif callable(prior):
            self._sample_density_profile_slope = prior
        else:
            raise ValueError("Invalid input for sample_density_profile_slope. Must be a string or a callable function.")

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

