# -*- coding: utf-8 -*-
"""
This module contains the LensGalaxyPopulation class, which is used to sample lens galaxy parameters, source parameters conditioned on the source being strongly lensed, image properties, and lensed SNRs. \n
The class inherits from the CBCSourceParameterDistribution class, which is used to sample source parameters. \n
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



from ..utils import  interpolator_from_pickle, cubic_spline_interpolator

class ImageProperties():
    """
    Class to find the image properties of a lensed event. Image properties include image positions, magnifications, time delays, etc.

    Parameters
    ----------
    npool : `int`
        number of processes to use
        default: 4
    z_min : `float`
        minimum redshift to consider
        default: 0.0
    z_max : `float`
        maximum redshift to consider
        default: 10.0
    n_min_images : `int`
        minimum number of images to consider
        default: 2
    n_max_images : `int`
        maximum number of images to consider
        default: 4
    geocent_time_min : `float`
        minimum geocent time to consider
        default: 1126259462.4 , which is the GPS time of the first GW detection
    geocent_time_max : `float`
        maximum geocent time to consider
        default: 1126259462.4+365*24*3600*100 , which is the GPS time of the first GW detection + 100 years. Some time delays can be very large.
    lens_model_list : `list`
        list of lens models
        default: ['EPL_NUMBA', 'SHEAR']
    cosmology : `astropy.cosmology`
        cosmology
        default: None/astropy.cosmology.LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
    spin_zero : `bool`
        whether to assume spin zero or not
        default: True
    spin_precession : `bool`
        whether to assume spin precession or not
        default: False
    directory : `str`
        directory to save the interpolator pickle files
        default: "./interpolator_pickle"
    create_new_interpolator : `dict`
        dictionary to create new interpolator pickle files
        default: dict(luminosity_distance_to_z=dict(create_new=False, resolution=1000))

    Examples
    --------
    >>> from ler.image_properties import ImageProperties
    >>> image_properties = ImageProperties()
    >>> lens_parameters = dict(zs=2.0, zl=0.5, gamma1=0.0, gamma2=0.0, e1=0.0, e2=0.0, gamma=2.0, theta_E=1.0)
    >>> lens_parameters = image_properties.image_properties(lens_parameters)
    >>> print(lens_parameters.keys())

    Instance Attributes
    ----------
    ImageProperties has the following instance attributes:\n
    +-------------------------------------+----------------------------------+
    | Atrributes                          | Type                             |
    +=====================================+==================================+
    |:attr:`npool`                        | `int`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`z_min`                        | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`z_max`                        | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`n_min_images`                 | `int`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`n_max_images`                 | `int`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`geocent_time_min`             | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`geocent_time_max`             | `float`                          |
    +-------------------------------------+----------------------------------+
    |:attr:`lens_model_list`              | `list`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`cosmo`                        | `astropy.cosmology`              |
    +-------------------------------------+----------------------------------+
    |:attr:`spin_zero`                    | `bool`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`spin_precession`              | `bool`                           |
    +-------------------------------------+----------------------------------+
    |:attr:`directory`                    | `str`                            |
    +-------------------------------------+----------------------------------+
    |:attr:`create_new_interpolator`      | `dict`                           |
    +-------------------------------------+----------------------------------+

    """

    def __init__(self, 
                 npool=4,
                 z_min=0.0,
                 z_max=10,
                 n_min_images=2, 
                 n_max_images=4,
                 time_window=365*24*3600*20,
                 lens_model_list=['EPL_NUMBA', 'SHEAR'],
                 cosmology=None,
                 spin_zero=True,
                 spin_precession=False,
                 directory="./interpolator_pickle",
                 create_new_interpolator=False,
        ):

        self.npool = npool
        self.n_min_images = n_min_images
        self.n_max_images = n_max_images
        self.lens_model_list = lens_model_list  # list of lens models
        self.spin_zero = spin_zero
        self.spin_precession = spin_precession
        self.time_window = time_window
        # self.geocent_time_min = geocent_time_min
        # self.geocent_time_max = geocent_time_max
        self.cosmo = cosmology if cosmology else cosmo
        
        # initialize the interpolator's parameters
        self.create_new_interpolator = dict(
            luminosity_distance_to_z=dict(create_new=False, resolution=1000),
        )
        if isinstance(create_new_interpolator, dict):
            self.create_new_interpolator.update(create_new_interpolator)
        elif create_new_interpolator is True:
            self.create_new_interpolator = dict(
                luminosity_distance_to_z=dict(create_new=True, resolution=1000)
            )

        resolution = self.create_new_interpolator["luminosity_distance_to_z"]["resolution"]
        create_new = self.create_new_interpolator["luminosity_distance_to_z"]["create_new"]
        spline1 = interpolator_from_pickle(
            param_dict_given= dict(z_min=z_min, z_max=z_max, cosmology=self.cosmo, resolution=resolution),
            directory=directory,
            sub_directory="luminosity_distance_to_z",
            name="luminosity_distance_to_z",
            x = np.linspace(z_min, z_max, resolution),
            pdf_func= lambda z_: cosmo.luminosity_distance(z_).value, 
            conditioned_y=None, 
            dimension=1,
            category="function_inverse",
            create_new=create_new,
        )
        self.luminosity_distance_to_z = njit(lambda z_: cubic_spline_interpolator(z_, spline1[0], spline1[1]))

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
                source related=>['mass_1': mass in detector frame (mass1>mass2), 'mass_2': mass in detector frame, 'mass_1_source':mass in source frame, 'mass_2_source':mass source frame, 'luminosity_distance': luminosity distance, 'theta_jn': inclination angle, 'psi': polarization angle, 'phase': coalesence phase, 'geocent_time': coalensence GPS time at geocenter, 'ra': right ascension, 'dec': declination, 'a_1': spin magnitude of the more massive black hole, 'a2': spin magnitude of the less massive black hole, 'tilt_1': tilt angle of the more massive black hole, 'tilt_2': tilt angle of the less massive black hole, 'phi_12': azimuthal angle between the two spins, 'phi_jl': azimuthal angle between the total angular momentum and the orbital angular momentum]\n
                image related=>['x_source': source position in the x-direction, 'y_source': source position in the y-direction, 'x0_image_position': image position in the x-direction, 'x1_image_position': image position in the y-direction, 'magnifications': magnifications, 'time_delays': time delays, 'n_images': number of images formed, 'determinant': determinants, 'trace': traces, 'iteration': to keep track of the iteration number
        """

        npool = self.npool
        n_min_images = self.n_min_images
        n_max_images = self.n_max_images
        lensModelList = self.lens_model_list
        
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

        # Solve the lens equation
        print("solving lens equations...")
        if n_min_images == 2:
            solve_lens_equation_ = solve_lens_equation
        elif n_min_images > 2:
            print("n_min_images > 2 is not supported yet")
            raise NotImplementedError
        else:
            raise ValueError("n_min_images should be greater than 1")
        with Pool(processes=npool) as pool:
            # call the same function with different data in parallel
            # imap->retain order in the list, while map->doesn't
            for result in tqdm(
                pool.imap_unordered(solve_lens_equation_, input_arguments),
                total=len(input_arguments),
                ncols=100,
                disable=False,
            ):
                # print(result)
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

        # time-delays: convert to positive values
        # time-delays will be relative to the first arrived signal of an lensed event
        time_delays = time_delays - np.array([np.sort(time_delays, axis=1)[:, 0]]).T # this is alright if time delays are already sorted
        
        # select only strongly lensed events are selected
        # assert np.all(n_images >= 2), "There are events with no images!"

        # image type classification (morse phase)
        number_of_lensed_events = size
        image_type = np.zeros((number_of_lensed_events, n_max_images))
        # image_type[traces < 0] = 3
        # image_type[traces > 0] = 1
        # image_type[determinants < 0] = 2
        for i in range(number_of_lensed_events):
            for j in range(n_max_images):
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

        ############ Important ############
        # sorting wrt time delays
        idx_sort = np.argsort(time_delays, axis=1) # sort each row
        # idx_sort has the shape (number_of_lensed_events, n_max_images)
        for key, value in image_parameters.items():
            # sort each row
            image_parameters[key] = np.array([value[i,idx_sort[i]] for i in range(number_of_lensed_events)])
        ####################################
        lens_parameters.update(image_parameters)
        lens_parameters["n_images"] = n_images
        lens_parameters["x_source"] = x_source
        lens_parameters["y_source"] = y_source

        return lens_parameters

    def get_lensed_snrs(self, lensed_param, list_of_detectors=None, snr_calculator=None, pdet_calculator=None):
        """
        Function to calculate the signal to noise ratio for each image in each event.

        Parameters
        ----------
            snr_calculator : `function`
                snr function, as describe in the :class:`~ler.rates.GWRATES` class.
            list_of_detectors : `list`
                list of detectors
                e.g. ['H1', 'L1', 'V1']
            lensed_param : `dict`
                dictionary containing the both already lensed source paramters and image parameters.
                e.g. lensed_param.keys() = ['mass_1', 'mass_2', 'zs', 'luminosity_distance', 'theta_jn', 'psi', 'phi', 'ra', 'dec', 'geocent_time', 'phase', 'a_1', 'a2', 'tilt_1', 'tilt_2', 'phi_12', 'phi_jl', 'magnifications', 'time_delays']
            n_max_images : `int`
                maximum number of images to consider
                default: 4

        Returns
        -------
            snrs : `dict`
                signal to noise ratio for each image in each event.
                (dictionary containing 'H1', 'L1', ..., and 'snr_net', which is the network snr, for each image as an array with dimensions (number_of_lensed_events,n_max_images) )

        """
        # needed to calculate effective luminosity distance and effective time delay
        n_max_images = self.n_max_images
        magnifications = lensed_param["magnifications"]
        time_delays = lensed_param["time_delays"]
        imgage_type = lensed_param["image_type"]
        size = len(magnifications)

        # image type to morse phase
        imgage_type[imgage_type==1.] = 0.
        imgage_type[imgage_type==2.] = np.pi/2

        # Get the binary parameters
        number_of_lensed_events = len(magnifications)
        mass_1, mass_2, theta_jn, psi, ra, dec, phase, a_1, a_2, tilt_1, tilt_2, phi_12, phi_jl  = (
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
        
        # setting up snr dictionary
        result_dict = dict()
        if snr_calculator:
            result_dict["snr_net"] = (
                np.ones((number_of_lensed_events, n_max_images)) * np.nan
            )
        # setting up pdet dictionary
        elif pdet_calculator:
            result_dict["pdet_net"] = (
                np.ones((number_of_lensed_events, n_max_images)) * np.nan
            )
        else:
            raise ValueError("snr_calculator or pdet_calculator not given")
        
        # if detector list are provided for snr or pdet calculation
        if list_of_detectors:
            for detector in list_of_detectors:
                result_dict[detector] = (
                    np.ones((number_of_lensed_events, n_max_images)) * np.nan
                )
        
        # for updating the lensed_param
        if "luminosity_distance" in lensed_param:
            luminosity_distance = lensed_param["luminosity_distance"]
            lensed_param["effective_luminosity_distance"] = np.ones((number_of_lensed_events, n_max_images)) * np.nan
            dl_eff_present = False
        elif "effective_luminosity_distance" in lensed_param:
            dl_eff = lensed_param["effective_luminosity_distance"]
            dl_eff_present = True
        else:
            raise ValueError("luminosity_distance or effective_luminosity_distance not given")
        
        if "geocent_time" in lensed_param:
            geocent_time = lensed_param["geocent_time"]
            lensed_param["effective_geocent_time"] = np.ones((number_of_lensed_events, n_max_images)) * np.nan
            time_eff_present = False
        elif "effective_geocent_time" in lensed_param:
            time_eff = lensed_param["effective_geocent_time"]
            time_eff_present = True
        else:
            raise ValueError("geocent_time or effective_geocent_time not given")
        
        if "phase" in lensed_param:
            phase = lensed_param["phase"]
            lensed_param["effective_phase"] = np.ones((number_of_lensed_events, n_max_images)) * np.nan
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

        for i in range(n_max_images):
            
            # get the effective time for each image type
            if not time_eff_present:
                effective_geocent_time = geocent_time + time_delays[:, i]
            else:
                effective_geocent_time = time_eff[:, i]
                
            # choose only the events that are within the time range and also not nan
            idx = (effective_geocent_time <= geocent_time_max) & (effective_geocent_time >= geocent_time_min)

            # get the effective luminosity distance for each image type
            if not dl_eff_present:
                effective_luminosity_distance = luminosity_distance / np.sqrt(
                    np.abs(magnifications[:, i])
                )
            else:
                effective_luminosity_distance = dl_eff[:, i]
            # get the effective phase for each image type
            if not phase_eff_present:
                effective_phase = phase - imgage_type[:, i]  # morse phase correction
            else:
                effective_phase = phase_eff[:, i]
            
            # check for nan values
            idx = idx & ~np.isnan(effective_luminosity_distance) & ~np.isnan(effective_geocent_time) & ~np.isnan(effective_phase)

            # Each image has their own effective luminosity distance and effective geocent time
            if sum(idx) != 0:
                # Returns a dictionary
                if snr_calculator:
                    optimal_snr = snr_calculator(
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
                        output_jsonfile=False,
                    )
                    result_dict["snr_net"][idx, i] = optimal_snr["snr_net"]

                    if list_of_detectors:
                        for detector in list_of_detectors:
                            if detector in optimal_snr:
                                result_dict[detector][idx, i] = optimal_snr[detector]

                elif pdet_calculator:
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
                

                lensed_param["effective_luminosity_distance"][:, i] = effective_luminosity_distance
                lensed_param["effective_geocent_time"][:, i] = effective_geocent_time
                lensed_param["effective_phase"][:, i] = effective_phase

        if dl_eff_present:
            del lensed_param["effective_luminosity_distance"]
        if time_eff_present:
            del lensed_param["effective_geocent_time"]
        if phase_eff_present:
            del lensed_param["effective_phase"]

        return result_dict, lensed_param
