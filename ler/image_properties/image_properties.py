# -*- coding: utf-8 -*-
"""
This module contains the LensGalaxyPopulation class, which is used to sample lens galaxy parameters, source parameters conditioned on the source being strongly lensed, image properties, and lensed SNRs. \n
The class inherits from the CompactBinaryPopulation class, which is used to sample source parameters. \n
"""
import warnings
warnings.filterwarnings("ignore")
import numpy as np

# the following .py file will be called if they are not given in the class initialization
from .multiprocessing_routine import solve_lens_equation1, solve_lens_equation2

# for multiprocessing
from multiprocessing import Pool
from tqdm import tqdm


class ImageProperties():
    """
    Class to find the image properties of a lensed event. Image properties include image positions, magnifications, time delays, etc.

    Parameters
    ----------
        npool : `int`
            number of processes to use
            default: 4
        n_min_images : `int`
            minimum number of images to consider
            default: 2
        n_max_images : `int`
            maximum number of images to consider
            default: 4
        lens_model_list : `list`
            list of lens models
            default: ['EPL_NUMBA', 'SHEAR']

    """

    def __init__(self, npool=4, n_min_images=2, n_max_images=4, lens_model_list=['EPL_NUMBA', 'SHEAR']):

        self.npool = 4
        self.n_min_images = n_min_images
        self.n_max_images = n_max_images
        self.lens_model_list = lens_model_list  # list of lens models

    def image_properties(
            self,
            lens_parameters
    ):
        """
        Function to get the image properties e.g. image positions, magnifications, time delays, etc.

        Parameters
        ----------
            lens_parameters : `dict`
                dictionary of lens parameters
                e.g. lens_parameters.keys() = ['zs', 'zl', 'gamma1', 'gamma2', 'e1', 'e2', 'gamma', 'theta_E']

        Returns
        -------
            lens_parameters : `dict`
                dictionary of lens parameters and image properties
                e.g. lens_parameters contains the following keys:\n
                lens related=>['zs': source redshift, 'zl': lens redshift, 'gamma1': shear component in the x-direction, 'gamma2': shear component in the y-direction, 'e1': ellipticity component in the x-direction, 'e2': ellipticity component in the y-direction, 'gamma': spectral index of the mass density distribution, 'theta_E': einstein radius in radian]\n
                source related=>['mass_1': mass in detector frame (mass1>mass2), 'mass_2': mass in detector frame, 'mass_1_source':mass in source frame, 'mass_2_source':mass source frame, 'luminosity_distance': luminosity distance, 'iota': inclination angle, 'psi': polarization angle, 'phase': coalesence phase, 'geocent_time': coalensence GPS time at geocenter, 'ra': right ascension, 'dec': declination, 'a_1': spin magnitude of the more massive black hole, 'a2': spin magnitude of the less massive black hole, 'tilt_1': tilt angle of the more massive black hole, 'tilt_2': tilt angle of the less massive black hole, 'phi_12': azimuthal angle between the two spins, 'phi_jl': azimuthal angle between the total angular momentum and the orbital angular momentum]\n
                image related=>['x_source': source position in the x-direction, 'y_source': source position in the y-direction, 'x0_image_position': image position in the x-direction, 'x1_image_position': image position in the y-direction, 'magnifications': magnifications, 'time_delays': time delays, 'n_images': number of images formed, 'determinant': determinants, 'trace': traces, 'iteration': to keep track of the iteration number, 'weights': weights for the caustic considered]

        """

        npool = self.npool
        n_min_images = self.n_min_images
        n_max_images = self.n_max_images
        
        zs = lens_parameters["zs"]
        size = len(zs)
        zl = lens_parameters["zl"]
        # external shear params to the 'PEMD' galaxy lens
        gamma1, gamma2 = lens_parameters["gamma1"], lens_parameters["gamma2"]
        # ellipticity of the galaxy lens
        e1, e2 = lens_parameters["e1"], lens_parameters["e2"]
        gamma = lens_parameters["gamma"]
        einstein_radius = lens_parameters["theta_E"]
        # Create the lens model list (note: can be a different lens model for different samples)
        lensModelList = np.array(lensModelList) * np.ones(
            (size, len(lensModelList)), dtype=object
        )
        min_img_arr = n_min_images * np.ones((size), dtype=int)

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
        x0_image_positions = np.ones((size, n_max_images)) * np.nan
        x1_image_positions = np.ones((size, n_max_images)) * np.nan
        magnifications = np.ones((size, n_max_images)) * np.nan
        time_delays = np.ones((size, n_max_images)) * np.nan
        determinants = np.ones((size, n_max_images)) * np.nan
        traces = np.ones((size, n_max_images)) * np.nan
        n_images = np.ones(size, dtype=int)
        x_source, y_source = np.ones(size) * np.nan, np.ones(size) * np.nan
        weights = np.ones(size) * np.nan

        # Solve the lens equation
        print("solving lens equations...")
        if n_min_images == 2:
            solve_lens_equation = solve_lens_equation1
        elif n_min_images > 2:
            solve_lens_equation = solve_lens_equation2
        else:
            raise ValueError("n_min_images should be greater than 1")
        with Pool(processes=npool) as pool:
            # call the same function with different data in parallel
            # imap->retain order in the list, while map->doesn't
            for result in tqdm(
                pool.imap(solve_lens_equation, input_arguments),
                total=len(input_arguments),
                ncols=100,
                disable=False,
            ):
                # print(result)
                """
                for i in tqdm(range(size)):
                    result = self.solve_lens_equation(input_arguments[i])
                """
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
                    weights_i,
                ) = result

                n_image_i = min(n_image_i, n_max_images)
                n_images[iter_i] = n_image_i
                x0_image_position = np.ones(n_max_images) * np.nan
                x1_image_position = np.ones(n_max_images) * np.nan
                x0_image_position[:n_image_i] = x0_image_position_i[:n_image_i]
                x1_image_position[:n_image_i] = x1_image_position_i[:n_image_i]
                x0_image_positions[
                    iter_i
                ] = x0_image_position  # shape = (size, n_max_images)
                x1_image_positions[
                    iter_i
                ] = x1_image_position  # shape = (size, n_max_images)
                magnification = np.ones(n_max_images) * np.nan
                time_delay = np.ones(n_max_images) * np.nan
                determinant = np.ones(n_max_images) * np.nan
                trace = np.ones(n_max_images) * np.nan
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
                weights[iter_i] = weights_i

        # time-delays: convert to positive values
        # time-delays will be relative to the first arrived signal of an lensed event
        time_delays = time_delays - np.array([np.sort(time_delays, axis=1)[:, 0]]).T

        # select only strongly lensed events are selected
        assert np.all(n_images >= 2), "There are events with no images!"

        # image type classification (morse phase)
        number_of_lensed_events = size
        image_type = np.zeros((number_of_lensed_events, n_max_images))
        image_type[traces < 0] = 3
        image_type[traces > 0] = 1
        image_type[determinants < 0] = 2

        # Return a dictionary with all of the lens information but also the BBH parameters from gw_param
        image_parameters = {
            "n_images": n_images,
            "x0_image_positions": x0_image_positions,
            "x1_image_positions": x1_image_positions,
            "magnifications": magnifications,
            "time_delays": time_delays,
            "image_type": image_type,
            "weights": weights,
        }
        lens_parameters.update(image_parameters)

        return lens_parameters

    def get_lensed_snrs(self, snr_calculator, lensed_param, n_max_images=4):
        """
        Function to calculate the signal to noise ratio for each image in each event.

        Parameters
        ----------
            snr_calculator : `class`
                snr_calculator class
                this is an already initialized class that contains a function (snr_calculator.snr) that actually calculates snr with the given gw_params.\n
                Luminosity distance and time delay are modified to be effective luminosity distance and effective time delay, respectively, for each image using the magnifications and time delays.\n
            lensed_param : `dict`
                dictionary containing the both already lensed source paramters and image parameters.
                e.g. lensed_param.keys() = ['mass_1', 'mass_2', 'zs', 'luminosity_distance', 'iota', 'psi', 'phi', 'ra', 'dec', 'geocent_time', 'phase', 'a_1', 'a2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl', 'magnifications', 'time_delays']
            n_max_images : `int`
                maximum number of images to consider
                default: 4

        Returns
        -------
            snrs : `dict`
                signal to noise ratio for each image in each event.
                (dictionary containing 'H1', 'L1', ..., and 'opt_snr_net', which is the network snr, for each image as an array with dimensions (number_of_lensed_events,n_max_images) )

        """
        # needed to calculate effective luminosity distance and effective time delay
        magnifications = lensed_param["magnifications"]
        time_delays = lensed_param["time_delays"]

        # Get the binary parameters
        number_of_lensed_events = len(magnifications)
        mass_1, mass_2, luminosity_distance, iota, psi, ra, dec, geocent_time, phase, a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl  = (
            lensed_param["mass_1"],
            lensed_param["mass_2"],
            lensed_param["luminosity_distance"],
            lensed_param["iota"],
            lensed_param["psi"],
            lensed_param["ra"],
            lensed_param["dec"],
            lensed_param["geocent_time"],
            lensed_param["phase"],
            lensed_param["a_1"],
            lensed_param["a_2"],
            lensed_param["tilt_1"],
            lensed_param["tilt_2"],
            lensed_param["phi_12"],
            lensed_param["phi_jl"],
        )

        # setting up snr dictionary
        detectors = snr_calculator.list_of_detectors
        optimal_snrs = dict()
        optimal_snrs["opt_snr_net"] = (
            np.ones((number_of_lensed_events, n_max_images)) * np.nan
        )
        for detector in detectors:
            optimal_snrs[detector] = (
                np.ones((number_of_lensed_events, n_max_images)) * np.nan
            )

        # LALSimulation cannot handle NaN
        if snr_calculator.snr_type == "inner_product":
            print("There will be {} progress bar iteration".format(n_max_images))

        for i in range(n_max_images):
            # Get the optimal signal to noise ratios for each image
            buffer = magnifications[:, i]
            idx = ~np.isnan(buffer)  # index of not-nan
            effective_luminosity_distance = luminosity_distance[idx] / np.sqrt(
                np.abs(buffer[idx])
            )
            effective_geocent_time = geocent_time[idx] + time_delays[idx, i]
            # if GPS time is negative, shift it
            # by a year until it is positive
            effective_geocent_time[
                effective_geocent_time < 0
            ] += 31556952  # number of seconds in a year

            # Each image has their own effective luminosity distance and effective geocent time
            if len(effective_luminosity_distance) != 0:
                # Returns a dictionary
                optimal_snr = snr_calculator.snr(
                    mass_1[idx],
                    mass_2[idx],
                    effective_luminosity_distance,
                    iota[idx],
                    psi[idx],
                    phase[idx],
                    effective_geocent_time,
                    ra[idx],
                    dec[idx],
                    a_1[idx],
                    a_2[idx],
                    tilt_1[idx],
                    tilt_2[idx],
                    phi_12[idx],
                    phi_jl[idx],
                    jsonFile=False,
                )

                optimal_snrs["opt_snr_net"][idx, i] = optimal_snr["opt_snr_net"]
                for detector in detectors:
                    optimal_snrs[detector][idx, i] = optimal_snr[detector]

        return optimal_snrs