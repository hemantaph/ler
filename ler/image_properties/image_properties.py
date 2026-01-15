# -*- coding: utf-8 -*-
"""
Module for computing image properties of strongly lensed gravitational wave events.

This module contains the ImageProperties class, which computes image positions,
magnifications, time delays, and other lensing-related quantities for gravitational
wave sources that are strongly lensed by intervening galaxies. The class uses
multiprocessing to efficiently solve lens equations for large samples.

Usage:
    Basic workflow example:

    >>> from ler.image_properties import ImageProperties
    >>> ip = ImageProperties()
    >>> lens_parameters = dict(zs=np.array([2.0]), zl=np.array([0.5]), ...)
    >>> result = ip.image_properties(lens_parameters)

Copyright (C) 2026 Phurailatpam Hemanta Kumar. Distributed under MIT License.
"""

import warnings
warnings.filterwarnings("ignore")
import logging
logging.getLogger('numexpr.utils').setLevel(logging.ERROR)
# for multiprocessing
from multiprocessing import Pool
from tqdm import tqdm

import numpy as np
from numba import njit

from astropy.cosmology import LambdaCDM
cosmo = LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)

# the following .py file will be called if they are not given in the class initialization
from .multiprocessing_routine import solve_lens_equation
from ..lens_galaxy_population.lens_functions import phi_q2_ellipticity


class ImageProperties():
    """
    Class to compute image properties of strongly lensed gravitational wave events.

    This class solves the lens equation to find image positions, magnifications,
    time delays, and image types (morse phase) for strongly lensed sources. It uses
    multiprocessing for efficient computation of large samples.

    Key Features: \n
    - Solves lens equations using multiprocessing for efficiency \n
    - Computes image positions, magnifications, and time delays \n
    - Classifies image types using morse phase \n
    - Calculates detection probabilities for lensed images \n

    Parameters
    ----------
    npool : ``int``
        Number of processes for multiprocessing. \n
        default: 4
    n_min_images : ``int``
        Minimum number of images required for a valid lensing event. \n
        default: 2
    n_max_images : ``int``
        Maximum number of images to consider per event. \n
        default: 4
    time_window : ``float``
        Time window for lensed events (units: seconds). \n
        default: 365*24*3600*20 (20 years)
    lens_model_list : ``list``
        List of lens models to use. \n
        default: ['EPL_NUMBA', 'SHEAR']
    cosmology : ``astropy.cosmology`` or ``None``
        Cosmology for distance calculations. \n
        If None, uses default LambdaCDM. \n
        default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
    spin_zero : ``bool``
        If True, spin parameters are set to zero (no spin sampling). \n
        default: False
    spin_precession : ``bool``
        If True (and spin_zero=False), sample precessing spin parameters. \n
        If False (and spin_zero=False), sample aligned/anti-aligned spins. \n
        default: False

    Examples
    --------
    Basic usage:

    >>> from ler.image_properties import ImageProperties
    >>> ip = ImageProperties()
    >>> lens_parameters = dict(
    ...     zs=np.array([2.0]),
    ...     zl=np.array([0.5]),
    ...     gamma1=np.array([0.0]),
    ...     gamma2=np.array([0.0]),
    ...     phi=np.array([0.0]),
    ...     q=np.array([0.8]),
    ...     gamma=np.array([2.0]),
    ...     theta_E=np.array([1.0])
    ... )
    >>> result = ip.image_properties(lens_parameters)
    >>> print(result.keys())
    

    Instance Methods
    ----------
    ImageProperties has the following methods: \n
    +-----------------------------------------------------+------------------------------------------------+
    | Method                                              | Description                                    |
    +=====================================================+================================================+
    | :meth:`~image_properties`                           | Compute image properties for lensed events     |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~get_lensed_snrs`                            | Compute detection probability for lensed images|
    +-----------------------------------------------------+------------------------------------------------+

    Instance Attributes
    ----------
    ImageProperties has the following attributes: \n
    +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
    | Attribute                                           | Type                      | Unit     | Description                                    |
    +=====================================================+===========================+==========+================================================+
    | :attr:`~npool`                                      | ``int``                   |          | Number of multiprocessing workers              |
    +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
    | :attr:`~n_min_images`                               | ``int``                   |          | Minimum number of images required              |
    +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
    | :attr:`~n_max_images`                               | ``int``                   |          | Maximum number of images per event             |
    +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
    | :attr:`~time_window`                                | ``float``                 | s        | Time window for lensed events                  |
    +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
    | :attr:`~lens_model_list`                            | ``list``                  |          | List of lens models                            |
    +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
    | :attr:`~cosmo`                                      | ``astropy.cosmology``     |          | Cosmology for calculations                     |
    +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
    | :attr:`~spin_zero`                                  | ``bool``                  |          | Flag for zero spin assumption                  |
    +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
    | :attr:`~spin_precession`                            | ``bool``                  |          | Flag for spin precession                       |
    +-----------------------------------------------------+---------------------------+----------+------------------------------------------------+
    """

    def __init__(self, 
                 npool=4,
                 n_min_images=2, 
                 n_max_images=4,
                 lens_model_list=['EPL_NUMBA', 'SHEAR'],
                 cosmology=None,
                 time_window=365*24*3600*20,
                 spin_zero=True,
                 spin_precession=False,
        ):

        self.npool = npool
        self.n_min_images = n_min_images
        self.n_max_images = n_max_images
        self.lens_model_list = lens_model_list  # list of lens models
        self.spin_zero = spin_zero
        self.spin_precession = spin_precession
        self.time_window = time_window
        self.cosmo = cosmology if cosmology else cosmo

    def image_properties(self, lens_parameters):
        """
        Compute image properties for strongly lensed events.

        Solves the lens equation using multiprocessing to find image positions,
        magnifications, time delays, and image types for each lensing event.

        Parameters
        ----------
        lens_parameters : ``dict``
            Dictionary containing lens and source parameters with keys: \n
            - 'zs': source redshift (array) \n
            - 'zl': lens redshift (array) \n
            - 'gamma1': external shear component 1 (array) \n
            - 'gamma2': external shear component 2 (array) \n
            - 'phi': position angle of lens ellipticity (array) \n
            - 'q': axis ratio of lens (array) \n
            - 'gamma': power-law slope of mass density (array) \n
            - 'theta_E': Einstein radius in radians (array) \n

        Returns
        -------
        lens_parameters : ``dict``
            Updated dictionary with additional image properties: \n
            - 'x0_image_positions': x-coordinates of images (shape: size x n_max_images) \n
            - 'x1_image_positions': y-coordinates of images (shape: size x n_max_images) \n
            - 'magnifications': magnification factors (shape: size x n_max_images) \n
            - 'time_delays': time delays relative to first image (shape: size x n_max_images, units: s) \n
            - 'image_type': morse phase classification (1=minimum, 2=saddle, 3=maximum) \n
            - 'n_images': number of images per event (array) \n
            - 'x_source': source x-position (array) \n
            - 'y_source': source y-position (array) \n
        """

        zs = lens_parameters["zs"]
        size = len(zs)
        zl = lens_parameters["zl"]
        # external shear params to the 'PEMD' galaxy lens
        gamma1, gamma2 = lens_parameters["gamma1"], lens_parameters["gamma2"]
        # ellipticity of the galaxy lens
        phi = lens_parameters["phi"]
        q = lens_parameters["q"]
        e1, e2 = phi_q2_ellipticity(phi, q)

        gamma = lens_parameters["gamma"]
        einstein_radius = lens_parameters["theta_E"]
        # Create the lens model list (note: can be a different lens model for different samples)
        lensModelList = np.array(self.lens_model_list) * np.ones(
            (size, len(self.lens_model_list)), dtype=object
        )
        min_img_arr = self.n_min_images * np.ones((size), dtype=int)

        # get image properties (with Multiprocessing)
        iterations = np.arange(size)
        input_arguments = np.array(
            [
                min_img_arr,
                e1,
                e2,
                gamma,
                gamma1,
                gamma2,
                zl,
                zs,
                einstein_radius,
                iterations,
            ],
            dtype=object,
        ).T
        input_arguments = np.concatenate((input_arguments, lensModelList), axis=1)
        # Initialize the image positions and lens argument list here.
        x0_image_positions = np.ones((size, self.n_max_images)) * np.nan
        x1_image_positions = np.ones((size, self.n_max_images)) * np.nan
        magnifications = np.ones((size, self.n_max_images)) * np.nan
        time_delays = np.ones((size, self.n_max_images)) * np.nan
        determinants = np.ones((size, self.n_max_images)) * np.nan
        traces = np.ones((size, self.n_max_images)) * np.nan
        n_images = np.ones(size, dtype=int)
        x_source, y_source = np.ones(size) * np.nan, np.ones(size) * np.nan

        # Solve the lens equation
        print("solving lens equations...")
        if self.n_min_images == 2:
            solve_lens_equation_ = solve_lens_equation
        elif self.n_min_images > 2:
            print("n_min_images > 2 is not supported yet")
            raise NotImplementedError
        else:
            raise ValueError("n_min_images should be greater than 1")
        with Pool(processes=self.npool) as pool:
            # call the same function with different data in parallel
            # imap->retain order in the list, while map->doesn't
            for result in tqdm(
                pool.imap_unordered(solve_lens_equation_, input_arguments),
                total=len(input_arguments),
                ncols=100,
                disable=False,
            ):
                (
                    x_source_i,
                    y_source_i,
                    x0_image_position_i,
                    x1_image_position_i,
                    magnifications_i,
                    time_delays_i,
                    n_image_i,
                    determinant_i,
                    trace_i,
                    iter_i,
                ) = result

                n_image_i = min(n_image_i, self.n_max_images)
                n_images[iter_i] = n_image_i
                x0_image_position = np.ones(self.n_max_images) * np.nan
                x1_image_position = np.ones(self.n_max_images) * np.nan
                x0_image_position[:n_image_i] = x0_image_position_i[:n_image_i]
                x1_image_position[:n_image_i] = x1_image_position_i[:n_image_i]
                x0_image_positions[
                    iter_i
                ] = x0_image_position  # shape = (size, n_max_images)
                x1_image_positions[
                    iter_i
                ] = x1_image_position  # shape = (size, n_max_images)
                magnification = np.ones(self.n_max_images) * np.nan
                time_delay = np.ones(self.n_max_images) * np.nan
                determinant = np.ones(self.n_max_images) * np.nan
                trace = np.ones(self.n_max_images) * np.nan
                magnification[:n_image_i] = magnifications_i[:n_image_i]
                time_delay[:n_image_i] = time_delays_i[:n_image_i]
                determinant[:n_image_i] = determinant_i[:n_image_i]
                trace[:n_image_i] = trace_i[:n_image_i]
                # Add the magnifications, time delays, determinants, and traces to their respective arrays
                magnifications[iter_i] = magnification
                time_delays[iter_i] = time_delay
                determinants[iter_i] = determinant
                traces[iter_i] = trace
                x_source[iter_i] = x_source_i
                y_source[iter_i] = y_source_i

        # time-delays: convert to positive values
        # time-delays will be relative to the first arrived signal of an lensed event
        time_delays = time_delays - np.array([np.sort(time_delays, axis=1)[:, 0]]).T # this is alright if time delays are already sorted
        
        # image type classification (morse phase)
        image_type = np.zeros((size, self.n_max_images))
        for i in range(size):
            for j in range(self.n_max_images):
                if determinants[i, j] < 0:
                    image_type[i, j] = 2
                elif traces[i, j] > 0:
                    image_type[i, j] = 1
                elif traces[i, j] < 0:
                    image_type[i, j] = 3
                else:
                    image_type[i, j] = np.nan

        # Return a dictionary with all of the lens information but also the BBH parameters from gw_param
        image_parameters = {
            "x0_image_positions": x0_image_positions,
            "x1_image_positions": x1_image_positions,
            "magnifications": magnifications,
            "time_delays": time_delays,
            "image_type": image_type,
        }

        # sorting wrt time delays
        idx_sort = np.argsort(time_delays, axis=1) # sort each row
        # idx_sort has the shape (size, n_max_images)
        for key, value in image_parameters.items():
            # sort each row
            image_parameters[key] = np.array([value[i,idx_sort[i]] for i in range(size)])
        lens_parameters.update(image_parameters)
        lens_parameters["n_images"] = n_images
        lens_parameters["x_source"] = x_source
        lens_parameters["y_source"] = y_source

        return lens_parameters

    def get_lensed_snrs(self, lensed_param, pdet_calculator, list_of_detectors=None):
        """
        Compute detection probability for each lensed image.

        Calculates the effective luminosity distance, geocent time, and phase
        for each image accounting for magnification and morse phase, then
        computes detection probabilities using the provided calculator.

        Parameters
        ----------
        lensed_param : ``dict``
            Dictionary containing lensed source and image parameters with keys: \n
            - 'mass_1', 'mass_2': detector-frame masses (array) \n
            - 'luminosity_distance' or 'effective_luminosity_distance': distance (array) \n
            - 'geocent_time' or 'effective_geocent_time': GPS time (array) \n
            - 'phase' or 'effective_phase': coalescence phase (array) \n
            - 'theta_jn', 'psi', 'ra', 'dec': orientation and position (arrays) \n
            - 'magnifications': image magnifications (shape: size x n_max_images) \n
            - 'time_delays': image time delays (shape: size x n_max_images) \n
            - 'image_type': morse phase type (shape: size x n_max_images) \n
        pdet_calculator : ``callable``
            Function that computes detection probability given GW parameters.
        list_of_detectors : ``list`` or ``None``
            List of detector names (e.g., ['H1', 'L1', 'V1']) for per-detector results. \n
            default: None

        Returns
        -------
        result_dict : ``dict``
            Dictionary containing: \n
            - 'pdet_net': network detection probability (shape: size x n_max_images) \n
            - Individual detector probabilities if list_of_detectors provided \n
        lensed_param : ``dict``
            Updated dictionary with effective parameters: \n
            - 'effective_luminosity_distance': magnification-corrected distance \n
            - 'effective_geocent_time': time-delay-corrected GPS time \n
            - 'effective_phase': morse-phase-corrected coalescence phase \n
        """
        # needed to calculate effective luminosity distance and effective time delay
        magnifications = lensed_param["magnifications"]
        time_delays = lensed_param["time_delays"]
        image_type = lensed_param["image_type"].copy()  # copy to avoid modifying original
        size = len(magnifications)

        # image type to morse phase
        image_type[image_type==1.] = 0.
        image_type[image_type==2.] = np.pi/2
        image_type[image_type==3.] = np.pi
        # Get the binary parameters
        mass_1, mass_2, theta_jn, psi, ra, dec, phase, a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl = (
            lensed_param["mass_1"],
            lensed_param["mass_2"],
            lensed_param["theta_jn"],
            lensed_param["psi"],
            lensed_param["ra"],
            lensed_param["dec"],
            lensed_param["phase"],
            np.zeros(size),
            np.zeros(size),
            np.zeros(size),
            np.zeros(size),
            np.zeros(size),
            np.zeros(size),
        )
        if not self.spin_zero:
            a_1, a_2 = (lensed_param["a_1"], lensed_param["a_2"])
            if self.spin_precession:
                tilt_1, tilt_2, phi_12, phi_jl = (
                    lensed_param["tilt_1"],
                    lensed_param["tilt_2"],
                    lensed_param["phi_12"],
                    lensed_param["phi_jl"],
                )
        
        # setting up pdet dictionary
        result_dict = dict()
        result_dict["pdet_net"] = (
            np.ones((size, self.n_max_images)) * np.nan
        )
        
        # if detector list are provided for pdet calculation
        if list_of_detectors:
            for detector in list_of_detectors:
                result_dict[detector] = (
                    np.ones((size, self.n_max_images)) * np.nan
                )
        
        # for updating the lensed_param
        if "luminosity_distance" in lensed_param:
            luminosity_distance_ = lensed_param["luminosity_distance"]
            lensed_param["effective_luminosity_distance"] = np.ones((size, self.n_max_images)) * np.nan
            dl_eff_present = False
        elif "effective_luminosity_distance" in lensed_param:
            dl_eff = lensed_param["effective_luminosity_distance"]
            dl_eff_present = True
        else:
            raise ValueError("luminosity_distance or effective_luminosity_distance not given")
        
        if "geocent_time" in lensed_param:
            geocent_time = lensed_param["geocent_time"]
            lensed_param["effective_geocent_time"] = np.ones((size, self.n_max_images)) * np.nan
            time_eff_present = False
        elif "effective_geocent_time" in lensed_param:
            time_eff = lensed_param["effective_geocent_time"]
            time_eff_present = True
        else:
            raise ValueError("geocent_time or effective_geocent_time not given")
        
        if "phase" in lensed_param:
            lensed_param["effective_phase"] = np.ones((size, self.n_max_images)) * np.nan
            phase_eff_present = False
        elif "effective_phase" in lensed_param:
            phase_eff = lensed_param["effective_phase"]
            phase_eff_present = True
        else:
            raise ValueError("phase or effective_phase not given")


        # Get the optimal signal to noise ratios for each image
        # iterate over the image type (column)
        geocent_time_min = np.min(geocent_time)
        geocent_time_max = geocent_time_min + self.time_window

        for i in range(self.n_max_images):
            
            # get the effective time for each image type
            if not time_eff_present:
                effective_geocent_time = geocent_time + time_delays[:, i]
            else:
                effective_geocent_time = time_eff[:, i]
                
            # choose only the events that are within the time range and also not nan
            idx = (effective_geocent_time <= geocent_time_max) & (effective_geocent_time >= geocent_time_min)

            # get the effective luminosity distance for each image type
            if not dl_eff_present:
                effective_luminosity_distance = luminosity_distance_ / np.sqrt(
                    np.abs(magnifications[:, i])
                )
            else:
                effective_luminosity_distance = dl_eff[:, i]
            # get the effective phase for each image type
            if not phase_eff_present:
                effective_phase = phase - image_type[:, i]  # morse phase correction
            else:
                effective_phase = phase_eff[:, i]
            
            # check for nan values
            idx = idx & ~np.isnan(effective_luminosity_distance) & ~np.isnan(effective_geocent_time) & ~np.isnan(effective_phase)

            # Each image has their own effective luminosity distance and effective geocent time
            if sum(idx) != 0:
                # Returns a dictionary
                pdet = pdet_calculator(
                    gw_param_dict= dict(
                        mass_1=mass_1[idx],
                        mass_2=mass_2[idx],
                        luminosity_distance=effective_luminosity_distance[idx],
                        theta_jn=theta_jn[idx],
                        psi=psi[idx],
                        phase= effective_phase[idx],
                        geocent_time=effective_geocent_time[idx],
                        ra=ra[idx],
                        dec=dec[idx],
                        a_1=a_1[idx],
                        a_2=a_2[idx],
                        tilt_1=tilt_1[idx],
                        tilt_2=tilt_2[idx],
                        phi_12=phi_12[idx],
                        phi_jl=phi_jl[idx],
                    ),
                )
                result_dict["pdet_net"][idx, i] = pdet["pdet_net"]

                if list_of_detectors:
                    for detector in list_of_detectors:
                        if detector in pdet:
                            result_dict[detector][idx, i] = pdet[detector]

            # Update lensed_param with effective values only if they weren't already present
            if not dl_eff_present:
                lensed_param["effective_luminosity_distance"][:, i] = effective_luminosity_distance
            if not time_eff_present:
                lensed_param["effective_geocent_time"][:, i] = effective_geocent_time
            if not phase_eff_present:
                lensed_param["effective_phase"][:, i] = effective_phase

        return result_dict, lensed_param

    # -------------
    # Properties
    # -------------

    @property
    def npool(self):
        """
        Number of multiprocessing workers.

        Returns
        -------
        npool : ``int``
            Number of processes for multiprocessing. \n
            default: 4
        """
        return self._npool

    @npool.setter
    def npool(self, value):
        self._npool = value

    @property
    def n_min_images(self):
        """
        Minimum number of images required for a valid lensing event.

        Returns
        -------
        n_min_images : ``int``
            Minimum number of images required. \n
            default: 2
        """
        return self._n_min_images

    @n_min_images.setter
    def n_min_images(self, value):
        self._n_min_images = value

    @property
    def n_max_images(self):
        """
        Maximum number of images per event.

        Returns
        -------
        n_max_images : ``int``
            Maximum number of images to consider per event. \n
            default: 4
        """
        return self._n_max_images

    @n_max_images.setter
    def n_max_images(self, value):
        self._n_max_images = value

    @property
    def time_window(self):
        """
        Time window for lensed events.

        Returns
        -------
        time_window : ``float``
            Time window for lensed events (units: s). \n
            default: 365*24*3600*20 (20 years)
        """
        return self._time_window

    @time_window.setter
    def time_window(self, value):
        self._time_window = value

    @property
    def lens_model_list(self):
        """
        List of lens models to use.

        Returns
        -------
        lens_model_list : ``list``
            List of lens model names. \n
            default: ['EPL_NUMBA', 'SHEAR']
        """
        return self._lens_model_list

    @lens_model_list.setter
    def lens_model_list(self, value):
        self._lens_model_list = value

    @property
    def cosmo(self):
        """
        Astropy cosmology object for calculations.

        Returns
        -------
        cosmo : ``astropy.cosmology``
            Cosmology used for distance calculations. \n
            default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        """
        return self._cosmo

    @cosmo.setter
    def cosmo(self, value):
        self._cosmo = value

    @property
    def spin_zero(self):
        """
        Flag for zero spin assumption.

        Returns
        -------
        spin_zero : ``bool``
            Whether to assume zero spin for compact objects. \n
            If True, spin parameters are set to zero (no spin sampling). \n
            If False, spin parameters are sampled. \n
            default: False
        """
        return self._spin_zero

    @spin_zero.setter
    def spin_zero(self, value):
        self._spin_zero = value

    @property
    def spin_precession(self):
        """
        Flag for spin precession.

        Returns
        -------
        spin_precession : ``bool``
            Whether to include spin precession effects. \n
            If True (and spin_zero=False), sample precessing spin parameters. \n
            If False (and spin_zero=False), sample aligned/anti-aligned spins. \n
            default: False
        """
        return self._spin_precession

    @spin_precession.setter
    def spin_precession(self, value):
        self._spin_precession = value
