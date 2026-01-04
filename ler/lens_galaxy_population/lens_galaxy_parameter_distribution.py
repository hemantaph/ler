"""
This module contains the ``LensGalaxyParameterDistribution`` class.

The class is used to sample lens galaxy parameters and source parameters conditioned on
the source being strongly lensed.
It inherits from the ``ImageProperties`` class to calculate image properties (magnification,
time delays, source position, image position, morse phase).
It also inherits from ``CBCSourceParameterDistribution``.

Copyright (c) 2026 Phurailatpam Hemantakumar
License: MIT
"""

import warnings
import numpy as np
from scipy.integrate import quad
from numba import njit

# the following .py file will be called if they are not given in the class initialization
from ..gw_source_population import CBCSourceParameterDistribution
from .optical_depth import OpticalDepth
from ..image_properties import ImageProperties
from ..utils import is_njitted
warnings.filterwarnings("ignore")

class LensGalaxyParameterDistribution(CBCSourceParameterDistribution, ImageProperties, OpticalDepth):
    """
    Class to sample lens galaxy parameters and source parameters conditioned on the source being strongly lensed.

    This class deals with the distribution of lens galaxy parameters, such as velocity dispersion,
    axis ratio, axis rotation angle, shear, and density profile slope. It also handles the
    sampling of source parameters conditioned on the source being strongly lensed.

    Parameters
    ----------
    npool : int, optional
        Number of processors to use.
        Default is 4.
    z_min : float, optional
        Minimum redshift.
        Default is 0.0.
    z_max : float, optional
        Maximum redshift.
        Default is 10.0.
    cosmology : astropy.cosmology, optional
        Cosmology to use.
        Default is None, which falls back to ``astropy.cosmology.FlatLambdaCDM(H0=70, Om0=0.3)``.
    event_type : str, optional
        Type of event to generate. e.g. 'BBH', 'BNS', 'NSBH'.
        Default is 'BBH'.
    lens_type : str, optional
        Type of lens galaxy to generate.
        Default is 'epl_shear_galaxy'.
    lens_functions : dict, optional
        Dictionary of lens functions.
    lens_functions_params : dict, optional
        Dictionary of parameters for lens functions.
    lens_param_samplers : dict, optional
        Dictionary of lens parameter samplers.
    lens_param_samplers_params : dict, optional
        Dictionary of parameters for lens parameter samplers.
    directory : str, optional
        Directory to store the interpolators.
        Default is './interpolator_json'.
    create_new_interpolator : bool, optional
        If True, creates a new interpolator.
        Default is False.
    buffer_size : int, optional
        Buffer size for sampling lens parameters.
        Default is 1000.
    **kwargs
        Keyword arguments to pass to the parent classes.

    Attributes
    ----------
    npool : int
        Number of processors to use.
    z_min : float
        Minimum redshift.
    z_max : float
        Maximum redshift.
    cosmo : astropy.cosmology
        Cosmology object.
    event_type : str
        Type of event to generate.
    directory : str
        Directory to store the interpolators.
    create_new_interpolator : dict
        Dictionary to check if new interpolator is created.
    lens_param_samplers : dict
        Dictionary of lens parameter samplers.
    lens_param_samplers_params : dict
        Dictionary of lens parameter sampler parameters.
    lens_functions : dict
        Dictionary of lens functions.
    normalization_pdf_z_lensed : float
        Normalization constant of the pdf p(z) for lensed events.

    Examples
    --------
    >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
    >>> lens = LensGalaxyParameterDistribution()
    >>> lensed_params = lens.sample_lens_parameters(size=1000)
    >>> print(lensed_params.keys())
    """

    # Attributes
    cbc_pop = None
    """:class:`~ler.gw_source_population.CBCSourceParameterDistribution`: Inherited class for sampling source parameters."""

    z_min = None
    """float: Minimum redshift."""
    z_max = None
    """float: Maximum redshift."""

    m_min = None
    """float: Minimum mass in detector frame."""

    m_max = None
    """float: Maximum mass in detector frame."""

    normalization_pdf_z = None
    """float: Normalization constant of the pdf p(z)."""

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
        directory="./interpolator_json",
        create_new_interpolator=False,
        buffer_size=1000,
        **kwargs  # for initialization of CBCSourceParameterDistribution class and ImageProperties class
    ):
        """
        Construct an :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution`.

        Parameters
        ----------
        npool : int, optional
            Number of processors to use.
            Default is 4.
        z_min : float, optional
            Minimum redshift.
            Default is 0.0.
        z_max : float, optional
            Maximum redshift.
            Default is 10.0.
        cosmology : astropy.cosmology, optional
            Cosmology to use.
            Default is None.
        event_type : str, optional
            Type of event to generate.
            Default is "BBH".
        lens_type : str, optional
            Type of lens galaxy to generate.
            Default is "epl_shear_galaxy".
        lens_functions : dict, optional
            Dictionary of lens functions.
            Default is None.
        lens_functions_params : dict, optional
            Dictionary of parameters for lens functions.
            Default is None.
        lens_param_samplers : dict, optional
            Dictionary of lens parameter samplers.
            Default is None.
        lens_param_samplers_params : dict, optional
            Dictionary of parameters for lens parameter samplers.
            Default is None.
        directory : str, optional
            Directory to store the interpolators.
            Default is "./interpolator_json".
        create_new_interpolator : bool, optional
            If True, creates a new interpolator.
            Default is False.
        buffer_size : int, optional
            Buffer size for sampling lens parameters.
            Default is 1000.
        **kwargs
            Keyword arguments to pass to the parent classes.
        """

        print("\nInitializing LensGalaxyParameterDistribution class...\n")          
        self.event_type = event_type  # needed for the source population    

        # initializing parent classes
        # This will separate initialization for:
        # 1. OpticalDepth: lens galaxy population parameters and functions
        # 2. CBCSourceParameterDistribution: GW source population
        # 3. ImageProperties: strong lensing image properties
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

        # # interpolated cross section function is not very accurate for rejection sampling
        # if self.lens_functions["cross_section"] == "cross_section_epl_shear_interpolation":
        #     self.cross_section = self.cross_section_epl_shear_numerical
        # function to sample lens parameters
        self.sample_lens_parameters_routine = getattr(self, self.lens_functions['param_sampler_type']) 

        # function to select lens parameters conditioned on the source being strongly lensed
        self.cross_section_based_sampler = self._initialization_cross_section_sampler()

        # To find the normalization constant of the pdf p(z)
        # this under the assumption that the event is strongly lensed
        # Define the merger-rate density function
        def pdf_unnormalized_(z):
            return self.merger_rate_density_detector_frame(np.array([z])) * self.optical_depth.function(np.array([z]))

        def pdf_unnormalized(z):
            return pdf_unnormalized_(z)[0]

        self.normalization_pdf_z_lensed = quad(
            pdf_unnormalized,
            self.z_min,
            self.z_max
        )[0]
        """float: Normalization constant of the pdf p(z) for lensed events."""

    def _initialization_cross_section_sampler(self):
        """
        Initialize the cross-section based sampler.

        Returns
        -------
        cross_section_based_sampler : callable
            A function that samples lens parameters based on the cross-section.
        """
        # Get random variable samplers from the initialized objects
        sigma_rvs_ = self.velocity_dispersion.rvs
        q_rvs_ = self.axis_ratio.rvs
        phi_rvs = self.axis_rotation_angle.rvs
        gamma_rvs = self.density_profile_slope.rvs
        shear_rvs = self.external_shear.rvs

        sigma_pdf_ = self.velocity_dispersion.pdf

        cross_section_ = self.cross_section
        # check the number of inputs for cross_section function
        # and JIT compile the wrapper function if necessary
        if cross_section_.__code__.co_argcount == 4:
            cross_section = njit(lambda zs, zl, sigma, q, phi, gamma, gamma1, gamma2: cross_section_(zs, zl, sigma, q))
        elif cross_section_.__code__.co_argcount == 3:
            cross_section = njit(lambda zs, zl, sigma, q, phi, gamma, gamma1, gamma2: cross_section_(zs, zl, sigma))
        else:
            cross_section = cross_section_
        
        # JIT compile the random variable samplers and PDFs if they are not already JIT compiled
        # This is required for compatibility with the JIT compiled cross_section_based_sampler
        if self.velocity_dispersion.conditioned_y_array is None:
            if is_njitted(sigma_rvs_):
                sigma_rvs = njit(lambda size, zl: sigma_rvs_(size))
                sigma_pdf = njit(lambda sigma, zl: sigma_pdf_(sigma))
            else:
                sigma_rvs = lambda size, zl: sigma_rvs_(size)  # noqa: E731
                sigma_pdf = lambda sigma, zl: sigma_pdf_(sigma)  # noqa: E731
        else:
            sigma_rvs = sigma_rvs_
            sigma_pdf = sigma_pdf_
            
        if self.axis_ratio.conditioned_y_array is None:
            if is_njitted(q_rvs_):
                q_rvs = njit(lambda size, sigma: q_rvs_(size))
            else:
                q_rvs = lambda size, sigma: q_rvs_(size)  # noqa: E731
        else:
            q_rvs = q_rvs_

        # List of functions to check if they are JIT compiled
        list_ = [sigma_rvs, q_rvs, phi_rvs, gamma_rvs, shear_rvs, sigma_pdf, cross_section]

        # Check if all functions are JIT compiled
        use_njit_sampler = True
        for func in list_:
            if not is_njitted(func):
                print(f"Warning: {func.__name__} is not njitted.")
                use_njit_sampler = False

        sigma_min = self.lens_param_samplers_params['velocity_dispersion']['sigma_min']
        sigma_max = self.lens_param_samplers_params['velocity_dispersion']['sigma_max']

        # Choose the sampling method based on the configuration
        if self.lens_functions['cross_section_based_sampler']=='rejection_sampling_with_cross_section':
            # Use rejection sampling
            from .lens_parameter_sampler import create_rejection_sampler
            saftey_factor = self.lens_functions_params['cross_section_based_sampler']['saftey_factor']

            cross_section_based_sampler = create_rejection_sampler(
                sigma_max=sigma_max,
                sigma_rvs=sigma_rvs,
                q_rvs=q_rvs,
                phi_rvs=phi_rvs,
                gamma_rvs=gamma_rvs,
                shear_rvs=shear_rvs,
                cross_section=cross_section,
                saftey_factor=saftey_factor,
                use_njit_sampler=use_njit_sampler,
            )

        elif self.lens_functions['cross_section_based_sampler']=='importance_sampling_with_cross_section':
            # Use importance sampling
            from .lens_parameter_sampler import create_importance_sampler
            n_prop = self.lens_functions_params['cross_section_based_sampler']['n_prop']
            cross_section_based_sampler = create_importance_sampler(
                sigma_min=sigma_min, 
                sigma_max=sigma_max,
                q_rvs=q_rvs,
                phi_rvs=phi_rvs,
                gamma_rvs=gamma_rvs,
                shear_rvs=shear_rvs,
                sigma_pdf=sigma_pdf,
                cross_section=cross_section,
                n_prop=n_prop,
                use_njit_sampler=use_njit_sampler,
            )
        else:
            raise ValueError("Invalid cross_section_based_sampler")

        # # initialize the jitted cross section based sampler
        # size = 2
        # # sample zs source redshift
        # zs = self.sample_source_redshift_sl(size) 
        # # sample zl lens redshiftsd
        # zl = self.lens_redshift.rvs(size, zs)
        # cross_section_based_sampler(zl, zs);

        return cross_section_based_sampler  
        

    def class_initialization_lens(self, npool, z_min, z_max, cosmology, lens_type, lens_functions, lens_functions_params, lens_param_samplers,  lens_param_samplers_params,  directory, create_new_interpolator, params):
        """
        Initialize the LensGalaxyParameterDistribution class.

        Parameters
        ----------
        npool : int
            Number of processors to use for sampling.
        z_min : float
            Minimum redshift of the lens galaxy.
        z_max : float
            Maximum redshift of the lens galaxy.
        cosmology : astropy.cosmology
            Cosmology object.
        lens_type : str
            Type of the lens galaxy.
        lens_functions : dict
            Dictionary with the lens related functions.
        lens_functions_params : dict
            Dictionary with the parameters for the lens related functions.
        lens_param_samplers : dict
            Dictionary with the priors for the sampler.
        lens_param_samplers_params : dict
            Dictionary with the parameters for the priors of the sampler.
        directory : str
            Directory where the interpolators are saved.
        create_new_interpolator : bool
            If True, creates a new interpolator.
        params : dict
            Additional parameters for the ``CBCSourceParameterDistribution`` and ``ImageProperties`` classes.
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
        Sample lens galaxy parameters along with the source parameters, conditioned on the source being strongly lensed.

        Parameters
        ----------
        size : int, optional
            Number of lens parameters to sample.
            Default is 1000.

        Returns
        -------
        lens_parameters : dict
            Dictionary of sampled lens parameters and source parameters.
            Keys include ``zl``, ``zs``, ``sigma``, ``q``, ``theta_E``, ``phi``, ``e1``, ``e2``,
            ``gamma1``, ``gamma2``, ``gamma``, ``geocent_time``, ``phase``, ``psi``, ``theta_jn``,
            ``luminosity_distance``, ``mass_1_source``, ``mass_2_source``, ``ra``, ``dec``.

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
        
    def sample_all_routine_epl_shear_sl(self, size=1000):
        """
        Sample galaxy lens parameters. EPL shear cross section is used for rejection sampling.

        Parameters
        ----------
        size : int, optional
            Number of lens parameters to sample.
            Default is 1000.

        Returns
        -------
        lens_parameters : dict
            Dictionary of lens parameters and source parameters (lens conditions applied).
            Keys include ``zl``, ``zs``, ``sigma``, ``q``, ``theta_E``, ``phi``, ``e1``, ``e2``,
            ``gamma1``, ``gamma2``, ``gamma``.

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> lens.sample_all_routine_epl_shear_sl(size=1000)
        """

        # sample zs source redshift
        zs = self.sample_source_redshift_sl(size) 

        # sample zl lens redshifts
        # this is conditioned on the source redshift, as the lens must be at a lower redshift than the source
        zl = self.lens_redshift.rvs(size, zs)

        # sample other lens parameters conditioned on the source and lens redshifts
        # this uses the initialized cross-section based sampler (rejection or importance sampling)
        sigma, q, phi, gamma, gamma1, gamma2 = self.cross_section_based_sampler(zs, zl)

        # Create a dictionary of the le ns parameters
        lens_parameters = {
            "zl": zl,
            "zs": zs,
            "sigma": sigma,
            "theta_E": self.compute_einstein_radii(sigma, zl, zs),
            "q": q,
            "phi": phi,
            "gamma": gamma,
            "gamma1": gamma1,
            "gamma2": gamma2,
        }

        return lens_parameters

    def strongly_lensed_source_redshifts(self, size=1000):
        """
        Sample source redshifts, conditioned on the source being strongly lensed.

        Parameters
        ----------
        size : int, optional
            Number of lens parameters to sample.
            Default is 1000.

        Returns
        -------
        redshifts : numpy.ndarray
            Source redshifts conditioned on the source being strongly lensed.

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> lens.strongly_lensed_source_redshifts(size=1000)
        """

        z_max = self.z_max

        def zs_function(zs_sl):
            # Sample source redshifts from the source population
            # self.zs is from CBCSourceParameterDistribution class
            zs = self.zs(size)
            
            # Apply strong lensing condition using optical depth
            # Calculate optical depth for each source redshift
            tau = self.optical_depth(zs)
            # Calculate max optical depth at z_max (tau increases with z)
            tau_max = self.optical_depth(np.array([z_max]))[0]
            
            # Rejection sampling step
            # Generate random numbers for rejection criterion
            r = np.random.uniform(0, tau_max, size=len(zs)) 
            
            # Pick strongly lensed sources where random number < optical depth
            zs_sl += list(zs[r < tau])  # list concatenation

            # Check if we have enough samples
            if len(zs_sl) >= size:
                # Trim list to the requested size
                zs_sl = zs_sl[:size]
                return zs_sl
            else:
                # Recursive call to get more samples
                return zs_function(zs_sl)

        zs_sl = []

        return np.array(zs_function(zs_sl))

    def sample_all_routine_epl_shear_intrinsic(self, size=1000):
        """
        Sample galaxy lens parameters. EPL shear cross section is used for rejection sampling.

        Parameters
        ----------
        size : int, optional
            Number of lens parameters to sample.
            Default is 1000.

        Returns
        -------
        lens_parameters : dict
            Dictionary of lens parameters and source parameters (lens conditions applied).
            Keys include ``zl``, ``zs``, ``sigma``, ``q``, ``theta_E``, ``phi``, ``e1``, ``e2``,
            ``gamma1``, ``gamma2``, ``gamma``.

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> lens.sample_all_routine_epl_shear_intrinsic(size=1000)
        """
        
        # sample zs source redshift
        zs = self.zs(size)

        # sample zl lens redshifts
        zl = self.lens_redshift_intrinsic.rvs(size, zs)

        # sample velocity dispersion
        if self.velocity_dispersion is not None:
            sigma = self.velocity_dispersion.rvs(size, zl)
        else:
            sigma = self.velocity_dispersion.rvs(size)

        # sample q
        if self.axis_ratio is not None:
            q = self.axis_ratio.rvs(size, sigma)
        else:
            q = self.axis_ratio.rvs(size)

        # sample phi
        phi = self.axis_rotation_angle.rvs(size)

        # sample gamma
        gamma = self.density_profile_slope.rvs(size)

        # sample gamma1 and gamma2
        gamma1, gamma2 = self.external_shear.rvs(size)

        # Create a dictionary of the lens parameters
        lens_parameters = {
            "zl": zl,
            "zs": zs,
            "sigma": sigma,
            "q": q,
            "phi": phi,
            "gamma": gamma,
            "gamma1": gamma1,
            "gamma2": gamma2,
        }

        return lens_parameters


