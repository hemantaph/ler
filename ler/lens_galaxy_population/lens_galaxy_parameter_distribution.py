# -*- coding: utf-8 -*-
"""
Module for lens galaxy parameter distributions.

This module contains the ``LensGalaxyParameterDistribution`` class which samples
lens galaxy parameters and source parameters conditioned on the source being
strongly lensed. It handles the distribution of lens parameters such as velocity
dispersion, axis ratio, rotation angle, shear, and density profile slope.

Inheritance hierarchy:

- :class:`~ler.gw_source_population.CBCSourceParameterDistribution` \n
- :class:`~ler.image_properties.ImageProperties` \n
- :class:`~ler.lens_galaxy_population.OpticalDepth` \n

Usage:
    Basic workflow example:

    >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
    >>> lens = LensGalaxyParameterDistribution()
    >>> lensed_params = lens.sample_lens_parameters(size=1000)
    >>> print(lensed_params.keys())

Copyright (C) 2026 Phurailatpam Hemantakumar. Distributed under MIT License.
"""

import warnings
import numpy as np
from scipy.integrate import quad
from numba import njit

from .optical_depth import OpticalDepth
from .sampler_functions import _njit_checks

from ..gw_source_population import CBCSourceParameterDistribution
from ..image_properties import ImageProperties
warnings.filterwarnings("ignore")


class LensGalaxyParameterDistribution(CBCSourceParameterDistribution, ImageProperties, OpticalDepth):
    """
    Sample lens galaxy parameters conditioned on strong lensing.

    This class handles the distribution of lens galaxy parameters such as velocity
    dispersion, axis ratio, axis rotation angle, shear, and density profile slope.
    It samples source parameters conditioned on the source being strongly lensed
    using cross-section based rejection or importance sampling.

    Key Features: \n
    - Samples lens parameters using EPL+shear galaxy model \n
    - Supports rejection and importance sampling based on cross-section \n
    - Computes optical depth weighted source redshift distributions \n
    - Integrates with GW source population and image property calculations \n

    Parameters
    ----------
    npool : ``int``
        Number of processors to use for parallel sampling. \n
        default: 4
    z_min : ``float``
        Minimum redshift for source and lens populations. \n
        default: 0.0
    z_max : ``float``
        Maximum redshift for source and lens populations. \n
        default: 10.0
    cosmology : ``astropy.cosmology`` or ``None``
        Cosmology object for distance calculations. \n
        default: None (uses FlatLambdaCDM with H0=70, Om0=0.3)
    event_type : ``str``
        Type of compact binary coalescence event. \n
        Options: \n
        - 'BBH': Binary black hole \n
        - 'BNS': Binary neutron star \n
        - 'NSBH': Neutron star-black hole \n
        default: 'BBH'
    lens_type : ``str``
        Type of lens galaxy model to use. \n
        default: 'epl_shear_galaxy'
    lens_functions : ``dict`` or ``None``
        Dictionary specifying lens-related functions. \n
        default: None (uses defaults from OpticalDepth)
    lens_functions_params : ``dict`` or ``None``
        Parameters for lens functions. \n
        default: None
    lens_param_samplers : ``dict`` or ``None``
        Dictionary specifying lens parameter sampling functions. \n
        default: None (uses defaults from OpticalDepth)
    lens_param_samplers_params : ``dict`` or ``None``
        Parameters for lens parameter samplers. \n
        default: None
    directory : ``str``
        Directory for storing interpolator files. \n
        default: './interpolator_json'
    create_new_interpolator : ``bool``
        If True, recreates interpolators even if files exist. \n
        default: False
    buffer_size : ``int``
        Buffer size for batch sampling of lens parameters. \n
        default: 1000
    **kwargs : ``dict``
        Additional keyword arguments passed to parent classes: \n
        :class:`~ler.gw_source_population.CBCSourceParameterDistribution`, \n
        :class:`~ler.image_properties.ImageProperties`, \n
        :class:`~ler.lens_galaxy_population.OpticalDepth`.

    Examples
    --------
    Basic usage:

    >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
    >>> lens = LensGalaxyParameterDistribution()
    >>> lensed_params = lens.sample_lens_parameters(size=1000)
    >>> print(lensed_params.keys())


    Instance Methods
    ----------
    LensGalaxyParameterDistribution has the following methods: \n
    +-----------------------------------------------------+------------------------------------------------+
    | Method                                              | Description                                    |
    +=====================================================+================================================+
    | :meth:`~sample_lens_parameters`                     | Sample lens and source parameters              |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~sample_all_routine_epl_shear_intrinsic`    | Sample EPL+shear lens parameters from intrinsic |
    |                                                     | distributions                                  |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~sample_all_routine_epl_shear_sl`            | Sample EPL+shear lens parameters with strong   |
    |                                                     | lensing condition                              |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~strongly_lensed_source_redshifts`           | Sample source redshifts with lensing condition |
    +-----------------------------------------------------+------------------------------------------------+

    Instance Attributes
    ----------
    LensGalaxyParameterDistribution has the following attributes: \n
    +------------------------------------------------+----------------------+-------+------------------------------------------------+
    | Attribute                                      | Type                 | Unit  | Description                                    |
    +================================================+======================+=======+================================================+
    | :attr:`~npool`                                 | ``int``              |       | Number of processors for parallel computation  |
    +------------------------------------------------+----------------------+-------+------------------------------------------------+
    | :attr:`~z_min`                                 | ``float``            |       | Minimum redshift                               |
    +------------------------------------------------+----------------------+-------+------------------------------------------------+
    | :attr:`~z_max`                                 | ``float``            |       | Maximum redshift                               |
    +------------------------------------------------+----------------------+-------+------------------------------------------------+
    | :attr:`~cosmo`                                 | ``astropy.cosmology``|       | Cosmology object for calculations              |
    +------------------------------------------------+----------------------+-------+------------------------------------------------+
    | :attr:`~event_type`                            | ``str``              |       | Type of CBC event (BBH, BNS, NSBH)             |
    +------------------------------------------------+----------------------+-------+------------------------------------------------+
    | :attr:`~directory`                             | ``str``              |       | Path to interpolator storage directory         |
    +------------------------------------------------+----------------------+-------+------------------------------------------------+
    | :attr:`~lens_param_samplers`                   | ``dict``             |       | Dictionary of lens parameter sampler names     |
    +------------------------------------------------+----------------------+-------+------------------------------------------------+
    | :attr:`~lens_param_samplers_params`            | ``dict``             |       | Parameters for lens parameter samplers         |
    +------------------------------------------------+----------------------+-------+------------------------------------------------+
    | :attr:`~lens_functions`                        | ``dict``             |       | Dictionary of lens function names              |
    +------------------------------------------------+----------------------+-------+------------------------------------------------+
    | :attr:`~normalization_pdf_z_lensed`            | ``float``            |       | Normalization constant for lensed source z pdf |
    +------------------------------------------------+----------------------+-------+------------------------------------------------+
    """

    def __init__(
        self,
        npool=4,
        z_min=0.0,
        z_max=10.0,
        cosmology=None,
        event_type="BBH",
        lens_type="epl_shear_galaxy",
        lens_functions=None,
        lens_functions_params=None,
        lens_param_samplers=None,
        lens_param_samplers_params=None,
        directory="./interpolator_json",
        create_new_interpolator=False,
        buffer_size=1000,
        **kwargs
    ):
        print("\nInitializing LensGalaxyParameterDistribution class...\n")
        self.event_type = event_type

        # Initialize parent classes
        self._class_initialization_lens(
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
            params=kwargs,
        )

        # Function to sample source redshifts conditioned on strong lensing
        self.sample_source_redshift_sl = getattr(self, self.lens_param_samplers["source_redshift_sl"])

        # Function to sample lens parameters
        self.sample_lens_parameters_routine = getattr(self, self.lens_functions["param_sampler_type"])

        # Function to select lens parameters conditioned on strong lensing
        self.cross_section_based_sampler = self._initialization_cross_section_sampler()

        # Compute normalization constant for lensed event pdf
        def pdf_unnormalized_(z):
            return self.merger_rate_density_detector_frame(np.array([z])) * self.optical_depth.function(np.array([z]))

        def pdf_unnormalized(z):
            return pdf_unnormalized_(z)[0]

        self.normalization_pdf_z_lensed = float(quad(
            pdf_unnormalized,
            self.z_min,
            self.z_max
        )[0])

    def _initialization_cross_section_sampler(self):
        """
        Initialize the cross-section based lens parameter sampler.

        Returns
        -------
        cross_section_based_sampler : ``callable``
            Function that samples lens parameters weighted by cross-section.
        """
        # Get random variable samplers from initialized objects
        sigma_rvs = self.velocity_dispersion.rvs
        q_rvs = self.axis_ratio.rvs
        phi_rvs = self.axis_rotation_angle.rvs
        gamma_rvs = self.density_profile_slope.rvs
        shear_rvs = self.external_shear.rvs
        sigma_pdf = self.velocity_dispersion.pdf
        number_density = self.velocity_dispersion.function
        cross_section_function = self.cross_section


        if cross_section_function.__code__.co_argcount <= 4: # sie or sis cross-section
            gamma_rvs = self.density_profile_slope_sl.rvs
            shear_rvs = self.external_shear_sl.rvs
            # not yet implemented
            # q_rvs = self.axis_ratio_sl.rvs
            # phi_rvs = self.axis_rotation_angle_sl.rvs
            
        use_njit_sampler, dict_ = _njit_checks(sigma_rvs, q_rvs, phi_rvs, gamma_rvs, shear_rvs, sigma_pdf, number_density, cross_section_function)

        sigma_rvs = dict_['sigma_rvs']
        sigma_pdf = dict_['sigma_pdf']
        q_rvs = dict_['q_rvs']
        phi_rvs = dict_['phi_rvs']
        gamma_rvs = dict_['gamma_rvs']
        shear_rvs = dict_['shear_rvs']
        cross_section_function = dict_['cross_section_function']

        sigma_max = self.lens_param_samplers_params["velocity_dispersion"]["sigma_max"]

        # Choose sampling method based on configuration
        if self.lens_functions["cross_section_based_sampler"] == "rejection_sampling_with_cross_section":
            from .sampler_functions import create_rejection_sampler
            safety_factor = self.lens_functions_params["cross_section_based_sampler"]["safety_factor"]

            cross_section_based_sampler = create_rejection_sampler(
                sigma_max=sigma_max,
                sigma_rvs=sigma_rvs,
                q_rvs=q_rvs,
                phi_rvs=phi_rvs,
                gamma_rvs=gamma_rvs,
                shear_rvs=shear_rvs,
                cross_section=cross_section_function,
                safety_factor=safety_factor,
                use_njit_sampler=use_njit_sampler,
            )

        elif self.lens_functions["cross_section_based_sampler"] == "importance_sampling_with_cross_section":
            from .sampler_functions import create_importance_sampler
            sigma_min = self.lens_param_samplers_params["velocity_dispersion"]["sigma_min"]
            n_prop = self.lens_functions_params["cross_section_based_sampler"]["n_prop"]
            cross_section_based_sampler = create_importance_sampler(
                sigma_min=sigma_min,
                sigma_max=sigma_max,
                q_rvs=q_rvs,
                phi_rvs=phi_rvs,
                gamma_rvs=gamma_rvs,
                shear_rvs=shear_rvs,
                sigma_pdf=sigma_pdf,
                cross_section=cross_section_function,
                n_prop=n_prop,
                use_njit_sampler=use_njit_sampler,
            )
        else:
            raise ValueError("Invalid cross_section_based_sampler")

        return cross_section_based_sampler

    def _class_initialization_lens(
        self,
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
        params,
    ):
        """
        Initialize parent classes for lens galaxy parameter distribution.

        Parameters
        ----------
        npool : ``int``
            Number of processors for parallel sampling.
        z_min : ``float``
            Minimum redshift.
        z_max : ``float``
            Maximum redshift.
        cosmology : ``astropy.cosmology`` or ``None``
            Cosmology object.
        lens_type : ``str``
            Type of lens galaxy model.
        lens_functions : ``dict`` or ``None``
            Dictionary of lens function names.
        lens_functions_params : ``dict`` or ``None``
            Parameters for lens functions.
        lens_param_samplers : ``dict`` or ``None``
            Dictionary of lens parameter sampler names.
        lens_param_samplers_params : ``dict`` or ``None``
            Parameters for lens parameter samplers.
        directory : ``str``
            Directory for interpolator storage.
        create_new_interpolator : ``bool``
            Whether to recreate interpolators.
        params : ``dict``
            Additional parameters for parent classes.
        """
        # Initialize OpticalDepth class
        OpticalDepth.__init__(
            self,
            npool=npool,
            z_min=z_min,
            z_max=z_max,
            cosmology=cosmology,
            lens_type=lens_type,
            lens_functions=lens_functions,
            lens_functions_params=lens_functions_params,
            lens_param_samplers=lens_param_samplers,
            lens_param_samplers_params=lens_param_samplers_params,
            directory=directory,
            create_new_interpolator=create_new_interpolator,
        )

        # Initialize CBCSourceParameterDistribution class
        input_params = dict(
            source_priors=None,
            source_priors_params=None,
            spin_zero=True,
            spin_precession=False,
        )
        input_params.update(params)

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

        # Initialize ImageProperties class
        input_params_image = dict(
            n_min_images=2,
            n_max_images=4,
            time_window=365 * 24 * 3600 * 2,
            lens_model_list=["EPL_NUMBA", "SHEAR"],
            effective_params_in_output=True,
        )
        input_params_image.update(params)

        ImageProperties.__init__(
            self,
            npool=self.npool,
            n_min_images=input_params_image["n_min_images"],
            n_max_images=input_params_image["n_max_images"],
            lens_model_list=input_params_image["lens_model_list"],
            cosmology=self.cosmo,
            time_window=input_params_image["time_window"],
            spin_zero=self.spin_zero,
            spin_precession=self.spin_precession,
            effective_params_in_output=input_params_image["effective_params_in_output"],
        )

    def sample_lens_parameters(self, size=1000):
        """
        Sample lens galaxy and source parameters conditioned on strong lensing.

        This method samples both lens galaxy parameters (velocity dispersion, axis
        ratio, shear, etc.) and gravitational wave source parameters, with the
        source redshift distribution weighted by strong lensing optical depth.

        Parameters
        ----------
        size : ``int``
            Number of lens-source parameter sets to sample. \n
            default: 1000

        Returns
        -------
        lens_parameters : ``dict``
            Dictionary containing sampled lens and source parameters. \n
            The included parameters and their units are as follows (for default settings):\n
            +------------------------------+-----------+-------------------------------------------------------+
            | Parameter                    | Units     | Description                                           |
            +==============================+===========+=======================================================+
            | zl                           |           | redshift of the lens                                  |
            +------------------------------+-----------+-------------------------------------------------------+
            | zs                           |           | redshift of the source                                |
            +------------------------------+-----------+-------------------------------------------------------+
            | sigma                        | km s^-1   | velocity dispersion                                   |
            +------------------------------+-----------+-------------------------------------------------------+
            | q                            |           | axis ratio                                            |
            +------------------------------+-----------+-------------------------------------------------------+
            | theta_E                      | radian    | Einstein radius                                       |
            +------------------------------+-----------+-------------------------------------------------------+
            | phi                          | rad       | axis rotation angle. counter-clockwise from the       |
            |                              |           | positive x-axis (RA-like axis) to the major axis of   |
            |                              |           | the projected mass distribution.                      |
            +------------------------------+-----------+-------------------------------------------------------+
            | gamma                        |           | density profile slope of EPL galaxy                   |
            +------------------------------+-----------+-------------------------------------------------------+
            | gamma1                       |           | external shear component in the x-direction           |
            |                              |           | (RA-like axis)                                        |
            +------------------------------+-----------+-------------------------------------------------------+
            | gamma2                       |           | external shear component in the y-direction           |
            |                              |           | (Dec-like axis)                                       |
            +------------------------------+-----------+-------------------------------------------------------+
            | geocent_time                 | s         | geocent time                                          |
            +------------------------------+-----------+-------------------------------------------------------+
            | ra                           | rad       | right ascension                                       |
            +------------------------------+-----------+-------------------------------------------------------+
            | dec                          | rad       | declination                                           |
            +------------------------------+-----------+-------------------------------------------------------+
            | phase                        | rad       | phase of GW at reference freq                         |
            +------------------------------+-----------+-------------------------------------------------------+
            | psi                          | rad       | polarization angle                                    |
            +------------------------------+-----------+-------------------------------------------------------+
            | theta_jn                     | rad       | inclination angle                                     |
            +------------------------------+-----------+-------------------------------------------------------+
            | a_1                          |           | spin of the primary compact binary                    |
            +------------------------------+-----------+-------------------------------------------------------+
            | a_2                          |           | spin of the secondary compact binary                  |
            +------------------------------+-----------+-------------------------------------------------------+
            | luminosity_distance          | Mpc       | luminosity distance of the source                      |
            +------------------------------+-----------+-------------------------------------------------------+
            | mass_1_source                | Msun      | mass of the primary compact binary (source frame)     |
            +------------------------------+-----------+-------------------------------------------------------+
            | mass_2_source                | Msun      | mass of the secondary compact binary (source frame)   |
            +------------------------------+-----------+-------------------------------------------------------+
            | mass_1                       | Msun      | mass of the primary compact binary (detector frame)   |
            +------------------------------+-----------+-------------------------------------------------------+
            | mass_2                       | Msun      | mass of the secondary compact binary (detector frame) |
            +------------------------------+-----------+-------------------------------------------------------+
            
        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> params = lens.sample_lens_parameters(size=1000)
        >>> print(params.keys())
        """
        print(f"sampling lens parameters with {self.lens_functions['param_sampler_type']}...")

        # Sample lens parameters with strong lensing condition
        lens_parameters = self.sample_lens_parameters_routine(size=size)

        # Sample GW source parameters (zs already sampled)
        param = dict(zs=lens_parameters["zs"])
        gw_param = self.sample_gw_parameters(size=size, param=param)

        # Combine lens and source parameters
        lens_parameters.update(gw_param)

        return lens_parameters

    def sample_all_routine_epl_shear_sl(self, size=1000):
        """
        Sample EPL+shear galaxy lens parameters with strong lensing condition.

        Parameters
        ----------
        size : ``int``
            Number of lens parameters to sample. \n
            default: 1000

        Returns
        -------
        lens_parameters : ``dict``
            Dictionary of sampled lens parameters. \n
            The included parameters and their units are as follows (for default settings):\n
            +------------------------------+-----------+-------------------------------------------------------+
            | Parameter                    | Units     | Description                                           |
            +==============================+===========+=======================================================+
            | zl                           |           | redshift of the lens                                  |
            +------------------------------+-----------+-------------------------------------------------------+
            | zs                           |           | redshift of the source                                |
            +------------------------------+-----------+-------------------------------------------------------+
            | sigma                        | km s^-1   | velocity dispersion                                   |
            +------------------------------+-----------+-------------------------------------------------------+
            | q                            |           | axis ratio                                            |
            +------------------------------+-----------+-------------------------------------------------------+
            | theta_E                      | radian    | Einstein radius                                       |
            +------------------------------+-----------+-------------------------------------------------------+
            | phi                          | rad       | axis rotation angle. counter-clockwise from the       |
            |                              |           | positive x-axis (RA-like axis) to the major axis of   |
            |                              |           | the projected mass distribution.                      |
            +------------------------------+-----------+-------------------------------------------------------+
            | gamma                        |           | density profile slope of EPL galaxy                   |
            +------------------------------+-----------+-------------------------------------------------------+
            | gamma1                       |           | external shear component in the x-direction           |
            |                              |           | (RA-like axis)                                        |
            +------------------------------+-----------+-------------------------------------------------------+
            | gamma2                       |           | external shear component in the y-direction           |
            |                              |           | (Dec-like axis)                                       |
            +------------------------------+-----------+-------------------------------------------------------+
        """
        # Sample source redshift with strong lensing weighting
        zs = self.sample_source_redshift_sl(size)

        # Sample lens redshift conditioned on source redshift
        zl = self.lens_redshift.rvs(size, zs)

        # Sample other lens parameters using cross-section based sampler
        sigma, q, phi, gamma, gamma1, gamma2 = self.cross_section_based_sampler(zs, zl)

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
        Sample source redshifts conditioned on strong lensing.

        Uses rejection sampling to generate source redshifts from the CBC source
        population weighted by the optical depth, which increases with redshift.

        Parameters
        ----------
        size : ``int``
            Number of redshifts to sample. \n
            default: 1000

        Returns
        -------
        redshifts : ``numpy.ndarray``
            Array of source redshifts conditioned on strong lensing.

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> zs = lens.strongly_lensed_source_redshifts(size=1000)
        >>> print(f"Mean source redshift: {zs.mean():.2f}")
        """
        z_max = self.z_max

        def zs_function(zs_sl):
            # Sample from source population
            zs = self.zs(size)

            # Apply optical depth weighting
            tau = self.optical_depth(zs)
            tau_max = self.optical_depth(np.array([z_max]))[0]

            # Rejection sampling
            r = np.random.uniform(0, tau_max, size=len(zs))
            zs_sl += list(zs[r < tau])

            if len(zs_sl) >= size:
                return zs_sl[:size]
            else:
                return zs_function(zs_sl)

        zs_sl = []
        return np.array(zs_function(zs_sl))

    def sample_all_routine_epl_shear_intrinsic(self, size=1000):
        """
        Sample EPL+shear galaxy lens parameters from intrinsic distributions.

        Samples lens parameters from their intrinsic distributions without
        applying strong lensing cross-section weighting.

        Parameters
        ----------
        size : ``int``
            Number of lens parameters to sample. \n
            default: 1000

        Returns
        -------
        lens_parameters : ``dict``
            Dictionary of sampled lens parameters. \n
            The included parameters and their units are as follows (for default settings):\n
            +------------------------------+-----------+-------------------------------------------------------+
            | Parameter                    | Units     | Description                                           |
            +==============================+===========+=======================================================+
            | zl                           |           | redshift of the lens                                  |
            +------------------------------+-----------+-------------------------------------------------------+
            | zs                           |           | redshift of the source                                |
            +------------------------------+-----------+-------------------------------------------------------+
            | sigma                        | km s^-1   | velocity dispersion                                   |
            +------------------------------+-----------+-------------------------------------------------------+
            | q                            |           | axis ratio                                            |
            +------------------------------+-----------+-------------------------------------------------------+
            | theta_E                      | radian    | Einstein radius                                       |
            +------------------------------+-----------+-------------------------------------------------------+
            | phi                          | rad       | axis rotation angle. counter-clockwise from the       |
            |                              |           | positive x-axis (RA-like axis) to the major axis of   |
            |                              |           | the projected mass distribution.                      |
            +------------------------------+-----------+-------------------------------------------------------+
            | gamma                        |           | density profile slope of EPL galaxy                   |
            +------------------------------+-----------+-------------------------------------------------------+
            | gamma1                       |           | external shear component in the x-direction           |
            |                              |           | (RA-like axis)                                        |
            +------------------------------+-----------+-------------------------------------------------------+
            | gamma2                       |           | external shear component in the y-direction           |
            |                              |           | (Dec-like axis)                                       |
            +------------------------------+-----------+-------------------------------------------------------+
        """
        # Sample source redshift from intrinsic distribution
        zs = self.zs(size)

        # Sample lens redshift
        zl = self.lens_redshift_intrinsic.rvs(size, self.z_max * np.ones(size))

        # Sample velocity dispersion
        if self.velocity_dispersion.conditioned_y_array is not None:
            sigma = self.velocity_dispersion.rvs(size, zl)
        else:
            sigma = self.velocity_dispersion.rvs(size)

        # Sample axis ratio
        if self.axis_ratio.conditioned_y_array is not None:
            q = self.axis_ratio.rvs(size, sigma)
        else:
            q = self.axis_ratio.rvs(size)

        # Sample remaining parameters
        phi = self.axis_rotation_angle.rvs(size)
        gamma = self.density_profile_slope.rvs(size)
        gamma1, gamma2 = self.external_shear.rvs(size)

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

    # -------------------------------------------------------------------------
    # Properties
    # -------------------------------------------------------------------------
    @property
    def lens_param_samplers(self):
        """
        Dictionary of lens parameter sampler function names.

        Returns
        -------
        lens_param_samplers : ``dict``
            Dictionary mapping parameter names to sampler function names. \n
            Keys include: 'source_redshift_sl', 'lens_redshift', \n
            'velocity_dispersion', 'axis_ratio', 'axis_rotation_angle', \n
            'external_shear', 'density_profile_slope'.
        """
        return self._lens_param_samplers

    @lens_param_samplers.setter
    def lens_param_samplers(self, value):
        self._lens_param_samplers = value

    @property
    def lens_param_samplers_params(self):
        """
        Dictionary of parameters for lens parameter samplers.

        Returns
        -------
        lens_param_samplers_params : ``dict``
            Dictionary with sampler parameters. \n
            Each key corresponds to a sampler in lens_param_samplers.
        """
        return self._lens_param_samplers_params

    @lens_param_samplers_params.setter
    def lens_param_samplers_params(self, value):
        self._lens_param_samplers_params = value

    @property
    def lens_functions(self):
        """
        Dictionary of lens-related function names.

        Returns
        -------
        lens_functions : ``dict``
            Dictionary mapping function types to function names. \n
            Keys include: 'param_sampler_type', 'cross_section_based_sampler', \n
            'optical_depth', 'cross_section'.
        """
        return self._lens_functions

    @lens_functions.setter
    def lens_functions(self, value):
        self._lens_functions = value

    @property
    def normalization_pdf_z_lensed(self):
        """
        Normalization constant for the lensed source redshift pdf.

        This constant is used to normalize the probability distribution \n
        of source redshifts conditioned on strong lensing. It is computed \n
        by integrating the merger rate density times optical depth.

        Returns
        -------
        normalization_pdf_z_lensed : ``float``
            Normalization constant for lensed redshift distribution.
        """
        return self._normalization_pdf_z_lensed

    @normalization_pdf_z_lensed.setter
    def normalization_pdf_z_lensed(self, value):
        self._normalization_pdf_z_lensed = value

