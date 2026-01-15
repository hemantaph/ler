# -*- coding: utf-8 -*-
"""
Module for solving lens equations using multiprocessing.

This sub-module contains functions to solve the lens equation for a given set
of lens parameters. The lens equation is solved using the analytical solver in
lenstronomy. These functions are used in the multiprocessing routine within
the ImageProperties class.

Usage:
    Basic workflow example:

    >>> from ler.image_properties.multiprocessing_routine import solve_lens_equation
    >>> import numpy as np
    >>> from multiprocessing import Pool
    >>> lens_params = np.array([2, 0.02, -0.01, 1.9, 0.1, 0.09, 0.25, 0.94, 1e-6, 0, 'EPL_NUMBA', 'SHEAR'], dtype=object)
    >>> result = solve_lens_equation(lens_params)

Copyright (C) 2026 Phurailatpam Hemanta Kumar. Distributed under MIT License.
"""

import numpy as np
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
from lenstronomy.LensModel.Solver.epl_shear_solver import caustics_epl_shear

# For sampling from caustic
from shapely.geometry import Polygon
import pointpats

# Constants
MAX_RETRIES = 100
MIN_MAGNIFICATION = 0.01


def _create_nan_result(iteration):
    """
    Create a NaN-filled result tuple for invalid lens configurations.

    Parameters
    ----------
    iteration : ``int``
        The iteration number to include in the result.

    Returns
    -------
    result : ``tuple``
        A tuple with NaN values for all output parameters.
    """
    nan_array = np.array([np.nan, np.nan])
    return (
        np.nan,
        np.nan,
        nan_array,
        nan_array,
        nan_array,
        nan_array,
        0,
        nan_array,
        nan_array,
        iteration,
    )


def solve_lens_equation(lens_parameters):
    """
    Solve the lens equation to find image properties.

    Uses the analytical solver from lenstronomy to find image positions,
    magnifications, time delays, and hessian properties for strongly
    lensed sources. Source positions are sampled from within the caustic
    region to ensure multiple imaging.

    Parameters
    ----------
    lens_parameters : ``numpy.ndarray``
        Array of lens configuration parameters with the following structure: \n
        - [0]: n_min_images - minimum number of images required \n
        - [1]: e1 - ellipticity component 1 \n
        - [2]: e2 - ellipticity component 2 \n
        - [3]: gamma - power-law slope of mass density \n
        - [4]: gamma1 - external shear component 1 \n
        - [5]: gamma2 - external shear component 2 \n
        - [6]: zl - lens redshift \n
        - [7]: zs - source redshift \n
        - [8]: einstein_radius - Einstein radius (units: arcsec) \n
        - [9]: iteration - iteration index for tracking \n
        - [10:]: lens_model_list - lens model names (e.g., 'EPL_NUMBA', 'SHEAR') \n

    Returns
    -------
    x_source : ``float``
        Source x-position (units: arcsec).
    y_source : ``float``
        Source y-position (units: arcsec).
    x0_image_position : ``numpy.ndarray``
        Image x-positions (units: arcsec).
    x1_image_position : ``numpy.ndarray``
        Image y-positions (units: arcsec).
    magnifications : ``numpy.ndarray``
        Magnification factors for each image.
    time_delays : ``numpy.ndarray``
        Time delays for each image (units: seconds).
    nImages : ``int``
        Number of images formed.
    determinant : ``numpy.ndarray``
        Determinant of the lensing Jacobian for each image.
    trace : ``numpy.ndarray``
        Trace of the lensing Jacobian for each image.
    iteration : ``int``
        Iteration index passed through for tracking.

    Examples
    --------
    >>> from ler.image_properties.multiprocessing_routine import solve_lens_equation
    >>> import numpy as np
    >>> from multiprocessing import Pool
    >>> lens_parameters1 = np.array([2, 0.024, -0.016, 1.89, 0.10, 0.09, 0.25, 0.94, 2.5e-06, 0, 'EPL_NUMBA', 'SHEAR'], dtype=object)
    >>> lens_parameters2 = np.array([2, -0.040, -0.014, 2.00, 0.08, -0.01, 1.09, 2.55, 1.0e-06, 1, 'EPL_NUMBA', 'SHEAR'], dtype=object)
    >>> input_arguments = np.vstack((lens_parameters1, lens_parameters2))
    >>> with Pool(2) as p:
    ...     result = p.map(solve_lens_equation, input_arguments)
    >>> print(f"Number of images: {result[0][6]}")
    """
    n_min_images = int(lens_parameters[0])
    zl = lens_parameters[6]
    zs = lens_parameters[7]
    einstein_radius = lens_parameters[8]
    iteration = lens_parameters[9]

    # Initialize lens model for image position, magnification, and time-delay calculations
    lensModel = LensModel(
        lens_model_list=lens_parameters[10:].tolist(), z_lens=zl, z_source=zs
    )
    lens_eq_solver = LensEquationSolver(lensModel)

    # Build lens model parameters
    factor = 1.0
    kwargs_lens = [
        {
            "theta_E": factor,
            "e1": lens_parameters[1],
            "e2": lens_parameters[2],
            "gamma": lens_parameters[3],
            "center_x": 0.0,
            "center_y": 0.0,
        },
        {
            "gamma1": lens_parameters[4],
            "gamma2": lens_parameters[5],
            "ra_0": 0,
            "dec_0": 0,
        },
    ]

    # Get the caustic curve for double-imaging region
    caustic_double_points = caustics_epl_shear(
        kwargs_lens, return_which="double", maginf=-100
    )
    caustic_is_valid = not np.isnan(caustic_double_points).any()

    if not caustic_is_valid:
        return _create_nan_result(iteration)

    # Define region where 2+ images form
    caustic_double = Polygon(caustic_double_points.T)

    # Sample source positions until strong lensing condition is satisfied
    strongly_lensed = False
    retry_count = 0
    while not strongly_lensed:
        x_source, y_source = pointpats.random.poisson(caustic_double, size=1)
        try:
            (
                x0_image_position,
                x1_image_position,
            ) = lens_eq_solver.image_position_from_source(
                sourcePos_x=x_source,
                sourcePos_y=y_source,
                kwargs_lens=kwargs_lens,
                solver="analytical",
                magnification_limit=MIN_MAGNIFICATION,
                arrival_time_sort=True,
            )
        except ValueError:
            return _create_nan_result(iteration)

        nImages = len(x0_image_position)
        if nImages >= n_min_images:
            strongly_lensed = True

        if retry_count > MAX_RETRIES:
            return _create_nan_result(iteration)
        retry_count += 1

    # Compute magnifications and time delays
    theta_E_nImages = einstein_radius * np.ones(nImages)
    radian_to_arcseconds = 180.0 / np.pi * 3600.0
    days_to_seconds = 24.0 * 3600.0

    magnifications = lensModel.magnification(
        x0_image_position, x1_image_position, kwargs_lens
    )
    time_delays = (
        lensModel.arrival_time(x0_image_position, x1_image_position, kwargs_lens)
        * (theta_E_nImages * radian_to_arcseconds) ** 2
        * days_to_seconds
    )

    # Compute hessian properties for image type classification
    hessian = lensModel.hessian(x0_image_position, x1_image_position, kwargs_lens)
    determinant = np.array(
        (1 - hessian[0]) * (1 - hessian[3]) - hessian[1] * hessian[2]
    )
    trace = np.array(2 - hessian[0] - hessian[3])

    # Scale positions to physical units
    x_source = x_source * einstein_radius
    y_source = y_source * einstein_radius
    x0_image_position = x0_image_position * einstein_radius
    x1_image_position = x1_image_position * einstein_radius

    return (
        x_source,
        y_source,
        x0_image_position,
        x1_image_position,
        magnifications,
        time_delays,
        nImages,
        determinant,
        trace,
        iteration,
    )
