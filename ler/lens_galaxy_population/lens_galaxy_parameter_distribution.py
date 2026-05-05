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
from astropy.cosmology import LambdaCDM

from .optical_depth import OpticalDepth
from ..gw_source_population import CBCSourceParameterDistribution
from ..image_properties import ImageProperties
from ..utils import (
    FunctionConditioning,
    generate_mixed_grid,
)

warnings.filterwarnings("ignore")


class LensGalaxyParameterDistribution(
    CBCSourceParameterDistribution, ImageProperties, OpticalDepth
):
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
    lens_priors : ``dict`` or ``None``
        Dictionary specifying lens parameter sampling functions. \n
        default: None (uses defaults from OpticalDepth)
    lens_priors_params : ``dict`` or ``None``
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
    | :meth:`~sample_all_routine_epl_shear_intrinsic`     | Sample EPL+shear lens parameters from          |
    |                                                     | intrinsic distributions                        |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~epl_shear_sl_parameters_rvs`            | Sample EPL+shear lens parameters with strong   |
    |                                                     | lensing condition                              |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~strongly_lensed_source_redshift`            | Sample source redshifts with lensing condition |
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
    | :attr:`~lens_priors`                   | ``dict``             |       | Dictionary of lens parameter sampler names     |
    +------------------------------------------------+----------------------+-------+------------------------------------------------+
    | :attr:`~lens_priors_params`            | ``dict``             |       | Parameters for lens parameter samplers         |
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
        lens_priors=None,
        lens_priors_params=None,
        directory="./interpolator_json",
        create_new_interpolator=False,
        buffer_size=1000,
        **kwargs,
    ):
        print("\nInitializing LensGalaxyParameterDistribution class...\n")
        cosmology = (
            cosmology
            if cosmology
            else LambdaCDM(
                H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0
            )
        )

        # Initialize parent classes
        self._class_initialization_lens(
            npool,
            z_min,
            z_max,
            cosmology,
            event_type,
            lens_type,
            lens_functions,
            lens_functions_params,
            lens_priors,
            lens_priors_params,
            directory,
            create_new_interpolator,
            params=kwargs,
        )

        # Function to sample source redshifts conditioned on strong lensing
        self.zs_sl = self.lens_priors["zs_sl"]

        # Function to sample lens parameters
        self.sample_lens_parameters_routine = getattr(
            self, self.lens_functions["param_sampler_type"]
        )

        # Function to select lens parameters conditioned on strong lensing
        self.cross_section_based_sampler = self._initialization_cross_section_sampler()

        # Compute normalization constant for lensed event pdf
        def pdf_unnormalized_(z):
            return self.merger_rate_density_detector_frame(
                np.array([z])
            ) * self.optical_depth.function(np.array([z]))

        def pdf_unnormalized(z):
            return pdf_unnormalized_(z)[0]

        self.normalization_pdf_z_lensed = float(
            quad(pdf_unnormalized, self.z_min, self.z_max)[0]
        )

    def _initialization_cross_section_sampler(self):
        """
        Initialize the cross-section based lens parameter sampler.

        ``_njit_checks`` always validates the intrinsic ``zs`` / ``lens_redshift``
        PDFs (as used in *full* modes). *Partial* modes draw ``(zs,zl)`` from
        ``zs_sl`` / ``lens_redshift_sl`` separately; conditional weights over
        ``(sigma, lambda)`` omit redshift and ``dV/dz`` factors (constant given
        fixed ``(zs,zl)``).

        Returns
        -------
        cross_section_based_sampler : ``callable``
            Function that samples lens parameters weighted by cross-section.
        """
        sampler_type = self.lens_functions["cross_section_based_sampler"]
        sampler_params = self.lens_functions_params["cross_section_based_sampler"]

        range_kwargs = dict(
            z_min=self.z_min,
            z_max=self.z_max,
            sigma_min=sampler_params.get(
                "sigma_min",
                self.lens_priors_params["velocity_dispersion"].get("sigma_min", 100.0),
            ),
            sigma_max=sampler_params.get(
                "sigma_max",
                self.lens_priors_params["velocity_dispersion"].get("sigma_max", 400.0),
            ),
            q_min=sampler_params.get(
                "q_min", self.lens_priors_params["axis_ratio"].get("q_min", 0.2)
            ),
            q_max=sampler_params.get(
                "q_max", self.lens_priors_params["axis_ratio"].get("q_max", 1.0)
            ),
            phi_min=sampler_params.get("phi_min", 0.0),
            phi_max=sampler_params.get("phi_max", 2 * np.pi),
            gamma_min=sampler_params.get("gamma_min", 1.5),
            gamma_max=sampler_params.get("gamma_max", 2.5),
            shear_min=sampler_params.get("shear_min", -0.2),
            shear_max=sampler_params.get("shear_max", 0.2),
        )

        # PDFs used in ``_njit_checks`` and full cross-section weights:
        # * Full samplers: intrinsic ``P(zs)``, ``dV/dz_l``, and Schechter-style
        #   ``n(sigma,zl)`` (no extra ``P(zl|zs)``: ``z_l`` is implicit in the
        #   uniform-in-interval proposal vs comoving-volume × number density).
        # * Partial samplers fix ``(zs,zl)`` from ``zs_sl`` and ``lens_redshift_sl``;
        #   weights over ``(sigma,lambda)`` use only ``n(sigma,zl)`` and shape PDFs
        #   (no ``P(zs|SL)``, ``P(zl|zs,SL)``, or ``dV/dz`` — constant given fixed redshifts).
        zs_pdf = self.zs.pdf
        zl_pdf = self.lens_redshift.pdf

        number_density = self.velocity_dispersion.function
        q_pdf=self.axis_ratio.pdf
        phi_pdf=self.axis_rotation_angle.pdf
        gamma_pdf=self.density_profile_slope.pdf
        shear1_pdf=self.external_shear1.pdf
        shear2_pdf=self.external_shear2.pdf
        cross_section=self.cross_section
        dVdz=self.differential_comoving_volume.function

        from .sampler_functions import _njit_checks

        use_njit_sampler, dict_ = _njit_checks(
            zs_pdf=zs_pdf,  
            zl_pdf=zl_pdf,
            number_density_function=number_density,
            q_pdf=q_pdf,
            phi_pdf=phi_pdf,
            gamma_pdf=gamma_pdf,
            shear1_pdf=shear1_pdf,
            shear2_pdf=shear2_pdf,
            cross_section_function=cross_section,
            dVcdz_function=dVdz,
            create_njit_sampler=True,
        )

        if use_njit_sampler:
            try:
                from numba import threading_layer

                if threading_layer() == "workqueue":
                    # Numba's workqueue backend is not safe for this parallel
                    # sampler path and can abort the interpreter. Keep the
                    # scalar njit callables, but use the interpreted sampler.
                    use_njit_sampler = False
            except ValueError:
                # Threading layer not initialized yet; leave the normal njit
                # path available and let Numba initialize the selected backend.
                pass
        
        from .sampler_functions import create_sampler

        _full = sampler_type in (
            "rejection_sampler_full",
            "importance_sampler_full",
        )

        return create_sampler(
            zs_pdf=dict_["zs_pdf"] if _full else None,
            number_density=dict_["number_density_function"],
            q_pdf=dict_["q_pdf"],
            phi_pdf=dict_["phi_pdf"],
            gamma_pdf=dict_["gamma_pdf"],
            shear1_pdf=dict_["shear1_pdf"],
            shear2_pdf=dict_["shear2_pdf"],
            cross_section=dict_["cross_section_function"],
            dVdz=dict_["dVcdz_function"],
            use_njit_sampler=use_njit_sampler,
            sampler_type=sampler_type,
            threshold_factor=sampler_params.get("threshold_factor", 1e-4),
            n_prop=sampler_params.get("n_prop", 50),
            lens_type=self.lens_type,
            npool=self.npool,
            **range_kwargs
        )


    def _class_initialization_lens(
        self,
        npool,
        z_min,
        z_max,
        cosmology,
        event_type,
        lens_type,
        lens_functions,
        lens_functions_params,
        lens_priors,
        lens_priors_params,
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
        lens_priors : ``dict`` or ``None``
            Dictionary of lens parameter sampler names.
        lens_priors_params : ``dict`` or ``None``
            Parameters for lens parameter samplers.
        directory : ``str``
            Directory for interpolator storage.
        create_new_interpolator : ``bool``
            Whether to recreate interpolators.
        params : ``dict``
            Additional parameters for parent classes.
        """

        # # print input params
        # print("z_min = ", z_min)
        # print("z_max = ", z_max)
        # print("cosmology = ", cosmology)
        # print("event_type = ", event_type)
        # print("lens_type = ", lens_type)
        # print("lens_functions = ", lens_functions)
        # print("lens_functions_params = ", lens_functions_params)
        # print("lens_priors = ", lens_priors)
        # print("lens_priors_params = ", lens_priors_params)
        # print("directory = ", directory)
        # print("create_new_interpolator = ", create_new_interpolator)

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
            lens_priors=lens_priors,
            lens_priors_params=lens_priors_params,
            directory=directory,
            create_new_interpolator=create_new_interpolator,
        )

        # Initialize CBCSourceParameterDistribution class
        input_params = dict(
            gw_priors=None,
            gw_priors_params=None,
            gw_functions=None,
            gw_functions_params=None,
            spin_zero=False,
            spin_precession=False,
        )
        input_params.update(params)

        CBCSourceParameterDistribution.__init__(
            self,
            npool=npool,
            z_min=z_min,
            z_max=z_max,
            event_type=event_type,
            gw_priors=input_params["gw_priors"],
            gw_priors_params=input_params["gw_priors_params"],
            gw_functions=input_params["gw_functions"],
            gw_functions_params=input_params["gw_functions_params"],
            spin_zero=input_params["spin_zero"],
            cosmology=cosmology,
            spin_precession=input_params["spin_precession"],
            directory=directory,
            create_new_interpolator=create_new_interpolator,
        )

        # Initialize ImageProperties class
        input_params_image = dict(
            n_min_images=2,
            n_max_images=4,
            time_window=365 * 24 * 3600 * 1,
            lens_model_list=["EPL_NUMBA", "SHEAR"],
            image_properties_function="image_properties_epl_shear_njit",
            image_properties_function_params=None,
            include_effective_parameters=False,
            multiprocessing_verbose=True,
            include_redundant_parameters=False,
        )
        input_params_image.update(params)

        ImageProperties.__init__(
            self,
            npool=npool,
            n_min_images=input_params_image["n_min_images"],
            n_max_images=input_params_image["n_max_images"],
            lens_model_list=input_params_image["lens_model_list"],
            image_properties_function=input_params_image["image_properties_function"],
            image_properties_function_params=input_params_image[
                "image_properties_function_params"
            ],
            cosmology=cosmology,
            directory=directory,
            z_min=z_min,
            z_max=z_max,
            create_new_interpolator=create_new_interpolator,
            time_window=input_params_image["time_window"],
            spin_zero=input_params["spin_zero"],
            spin_precession=input_params["spin_precession"],
            include_effective_parameters=input_params_image[
                "include_effective_parameters"
            ],
            multiprocessing_verbose=input_params_image["multiprocessing_verbose"],
            include_redundant_parameters=input_params_image[
                "include_redundant_parameters"
            ],
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
                | luminosity_distance          | Mpc       | luminosity distance of the source                     |
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

        print(
            f"sampling lens parameters with {self.lens_functions['param_sampler_type']}..."
        )

        # Sample lens parameters with strong lensing condition
        lens_parameters = self.sample_lens_parameters_routine(size=size)

        # Sample GW source parameters (zs already sampled)
        param = dict(zs=lens_parameters["zs"])
        gw_param = self.sample_gw_parameters(size=size, param=param)

        # Combine lens and source parameters
        lens_parameters.update(gw_param)

        return lens_parameters

    def epl_shear_sl_parameters_rvs(self, size=1000):
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

        sampler_type = self.lens_functions["cross_section_based_sampler"]

        # Full modes resample zs, zl, sigma, and lens nuisance parameters together.
        if sampler_type in [
            "rejection_sampler_full",
            "importance_sampler_full",
        ]:
            zs, zl, sigma, q, phi, gamma, gamma1, gamma2 = (
                self.cross_section_based_sampler(size)
            )

        # Partial modes keep the redshifts fixed and sample all lens parameters
        # from the configured uniform proposal ranges.
        elif sampler_type in [
            "rejection_sampler_partial",
            "importance_sampler_partial",
        ]:
            zs = self.zs_sl.rvs(size)
            zl = self.lens_redshift_sl.rvs(size, zs)
            sigma, q, phi, gamma, gamma1, gamma2 = self.cross_section_based_sampler(
                size, zs, zl
            )

        else:
            raise ValueError(
                f"Invalid cross_section_based_sampler: {sampler_type}. "
                "Available options are 'rejection_sampler_full', "
                "'importance_sampler_full', 'rejection_sampler_partial', "
                "and 'importance_sampler_partial'."
            )

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
        zl = self.lens_redshift.rvs(size, self.z_max * np.ones(size))

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
        gamma1 = self.external_shear1.rvs(size)
        gamma2 = self.external_shear2.rvs(size)

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

    # ---------------------------------
    # Source redshift sampler functions
    # ---------------------------------
    def strongly_lensed_source_redshift(self, size, get_attribute=False, **kwargs):
        """
        Sample source redshifts conditioned on strong lensing.

        Notes
        -----
        This sampler needs a merger-rate-weighted intrinsic source redshift model via
        ``self.merger_rate_density_based_source_redshift``. That attribute is provided by
        :class:`~ler.gw_source_population.CBCSourceParameterDistribution` and therefore is
        available in :class:`~ler.lens_galaxy_population.LensGalaxyParameterDistribution`
        (which inherits both).

        Parameters
        ----------
        size : ``int``
            Number of samples to generate. \n
            default: 1000
        get_attribute : ``bool``
            If True, returns the sampler object instead of samples. \n
            default: False
        **kwargs : ``dict``
            Additional parameters.

        Returns
        -------
        redshifts : ``numpy.ndarray`` or ``ler.utils.FunctionConditioning``
            Array of source redshifts conditioned on strong lensing, or the sampler object.

        Examples
        --------
        >>> from ler.lens_galaxy_population import LensGalaxyParameterDistribution
        >>> lens = LensGalaxyParameterDistribution()
        >>> zs = lens.strongly_lensed_source_redshift(size=1000)
        >>> print(f"strongly lensed source redshift: {zs.mean():.2f}")
        """

        if not hasattr(self, "merger_rate_density_based_source_redshift"):
            raise AttributeError(
                "Missing attribute `merger_rate_density_based_source_redshift`. "
                "This sampler must be called from a class that also inherits "
                "`ler.gw_source_population.CBCSourceParameterDistribution` "
                "(e.g. `LensGalaxyParameterDistribution`)."
            )

        identifier_dict = {"name": "strongly_lensed_source_redshift"}
        identifier_dict["resolution"] = self.create_new_interpolator["zs_sl"][
            "resolution"
        ]
        identifier_dict["cdf_size"] = self.create_new_interpolator["zs_sl"][
            "cdf_size"
        ]
        identifier_dict["z_min"] = self.z_min
        identifier_dict["z_max"] = self.z_max
        identifier_dict["cosmology"] = self.cosmo
        identifier_dict["resolution"] = self.create_new_interpolator["zs_sl"][
            "resolution"
        ]
        identifier_dict["optical_depth"] = self.optical_depth.info  # can be None
        identifier_dict["merger_rate_density_based_source_redshift"] = (
            self.merger_rate_density_based_source_redshift.info
        )

        param_dict = self.available_lens_priors["zs_sl"][
            "strongly_lensed_source_redshift"
        ]
        if param_dict:
            param_dict.update(kwargs)
        else:
            param_dict = kwargs
        identifier_dict.update(param_dict)

        z_min = self.z_min if self.z_min > 0.0 else 0.0001
        z_max = self.z_max
        resolution = identifier_dict["resolution"]
        zs_array = generate_mixed_grid(z_min, z_max, resolution)

        if param_dict["tau_approximation"]:

            def n_zs_sl(zs):
                # number density of strongly lensed sources
                P_sl_zs = self.optical_depth(zs)
                return P_sl_zs * self.merger_rate_density_based_source_redshift.function(
                    zs
                )

        else:

            def n_zs_sl(zs):
                # number density of strongly lensed sources
                tau = self.optical_depth(zs)
                P_sl_zs = tau * np.exp(-tau)
                return P_sl_zs * self.merger_rate_density_based_source_redshift.function(
                    zs
                )

        zs_sl_object = FunctionConditioning(
            function=n_zs_sl,
            x_array=zs_array,
            conditioned_y_array=None,
            identifier_dict=identifier_dict,
            directory=self.directory,
            sub_directory="zs_sl",
            name=identifier_dict["name"],
            create_new=self.create_new_interpolator["zs_sl"]["create_new"],
            create_function_inverse=False,
            create_function=True,
            create_pdf=True,
            create_rvs=True,
            callback="rvs",
            cdf_size=identifier_dict["cdf_size"],
        )

        return zs_sl_object if get_attribute else zs_sl_object.rvs(size)

    # -------------------------------------------------------------------------
    # Properties
    # -------------------------------------------------------------------------
    @property
    def zs_sl(self):
        """
        Function to sample source redshifts conditioned on strong lensing.

        Returns
        -------
        zs_sl : ``FunctionConditioning``
            Sampler object for strongly-lensed source redshift.
        """
        return self._zs_sl

    @zs_sl.setter
    def zs_sl(self, prior):
        if prior in self.available_lens_priors["zs_sl"]:
            print(f"using ler available zs_sl function : {prior}")
            args = self.lens_priors_params["zs_sl"]
            if args is None:
                self._zs_sl = getattr(self, prior)(
                    size=None, get_attribute=True
                )
            else:
                self._zs_sl = getattr(self, prior)(
                    size=None, get_attribute=True, **args
                )
        elif isinstance(prior, FunctionConditioning):
            print(
                "using user provided custom zs_sl class/object of type ler.utils.FunctionConditioning"
            )
            self._zs_sl = prior
        elif callable(prior):
            print("using user provided custom zs_sl sampler function")
            self._zs_sl = FunctionConditioning(
                function=None, x_array=None, create_rvs=prior
            )
        else:
            raise ValueError(
                "zs_sl should be string in available_lens_priors['zs_sl'] "
                "or class object of 'ler.utils.FunctionConditioning' "
                "or callable function with input argument 'size'"
            )

    @property
    def cross_section_based_sampler(self):
        """
        Cross-section based lens parameter sampler function. This is an initialized function of sampler_functions.create_importance_sampler or sampler_functions.create_rejection_sampler based on the configuration.

        This function samples lens parameters weighted by the strong lensing
        cross-section.

        Parameters
        ----------
        size : ``int``
            Number of strongly-lensed parameter sets to sample.

        Returns
        -------
        zs : ``numpy.ndarray``
            Array of source redshifts.
        zl : ``numpy.ndarray``
            Array of lens redshifts.
        sigma : ``numpy.ndarray``
            Array of velocity dispersions.
        q : ``numpy.ndarray``
            Array of axis ratios.
        phi : ``numpy.ndarray``
            Array of axis rotation angles.
        gamma : ``numpy.ndarray``
            Array of density profile slopes.
        gamma1 : ``numpy.ndarray``
            Array of external shear components in the x-direction.
        gamma2 : ``numpy.ndarray``
            Array of external shear components in the y-direction.

        """
        return self._cross_section_based_sampler

    @cross_section_based_sampler.setter
    def cross_section_based_sampler(self, value):
        self._cross_section_based_sampler = value

    @property
    def lens_priors(self):
        """
        Dictionary of lens parameter sampler function names.

        Returns
        -------
        lens_priors : ``dict``
            Dictionary mapping parameter names to sampler function names. \n
            Keys include: 'zs_sl', 'lens_redshift', \n
            'velocity_dispersion', 'axis_ratio', 'axis_rotation_angle', \n
            'external_shear1', 'external_shear2', 'density_profile_slope'.
        """
        return self._lens_priors

    @lens_priors.setter
    def lens_priors(self, value):
        self._lens_priors = value

    @property
    def lens_priors_params(self):
        """
        Dictionary of parameters for lens parameter samplers.

        Returns
        -------
        lens_priors_params : ``dict``
            Dictionary with sampler parameters. \n
            Each key corresponds to a sampler in lens_priors.
        """
        return self._lens_priors_params

    @lens_priors_params.setter
    def lens_priors_params(self, value):
        self._lens_priors_params = value

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
