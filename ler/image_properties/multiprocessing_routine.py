# -*- coding: utf-8 -*-
"""
This sub-module contains the functions to solve the lens equation for a given set of lens parameters. The lens equation is solved using the analytical solver in lenstronomy. The functions in this sub-module are used in the multiprocessing routine to solve the lens equation for a given set of lens parameters.
"""

import numpy as np
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
from lenstronomy.LensModel.Solver.epl_shear_solver import caustics_epl_shear

# For sampling from caustic
from shapely.geometry import Polygon
import pointpats


def solve_lens_equation(lens_parameters):
    """
    Function to solve the lens equation (min_image = 2)

    Parameters
    ----------
    lens_parameters : `list`
        a list of parameters
        lens_parameters[0] = min_images : minimum number of images
        lens_parameters[1] = e1 : ellipticity
        lens_parameters[2] = e2 : ellipticity
        lens_parameters[3] = gamma : power-law index
        lens_parameters[4] = gamma1 : shear
        lens_parameters[5] = gamma2 : shear
        lens_parameters[6] = zl : redshift of the lens
        lens_parameters[7] = zs : redshift of the source
        lens_parameters[8] = einstein_radius : Einstein radius
        lens_parameters[9] = iteration : iteration number
        lens_parameters[10:] = lens_model_list : numpy array of lens models

    Returns
    -------
    x_source : `float`
        x position of the source in the source plane
    y_source : `float`
        y position of the source in the source plane
    x0_image_position : `float`
        x position of the images in the source plane
    x1_image_position : `float`
        y position of the images in the source plane
    magnifications : `float`
        magnification of the images
    time_delays : `float`
        time-delay of the images
    nImages : `int`
        number of images
    determinant : `float`
        determinant of the hessian matrix
    trace : `float`
        trace of the hessian matrix
    iteration : `int`
        iteration number

    Examples
    --------
    >>> from ler.multiprocessing_routine import solve_lens_equation1
    >>> import numpy as np
    >>> from multiprocessing import Pool
    >>> # lens parameters input contains 12 parameters [e1, e2, gamma, gamma1, gamma2, zl, zs, einstein_radius, iteration, lens_model_list]
    >>> lens_parameters1 = np.array([2, 0.024069457093642648, -0.016002190961948142, 1.8945414936459974, 0.10117465203892329, 0.09600089396968613, 0.2503743800068136, 0.9418211055453296, 2.5055790287104725e-06, 0, 'EPL_NUMBA', 'SHEAR'], dtype=object)
    >>> lens_parameters2 = np.array([2, -0.04030088581646998, -0.01419438113690042, 2.0068239327017, 0.08482718989370612, -0.015393332086560785, 1.0952303138971118, 2.5534097159384417, 1.0125570159563301e-06, 1, 'EPL_NUMBA', 'SHEAR'], dtype=object)
    >>> input_arguments = np.vstack((lens_parameters1, lens_parameters2))
    >>> # solve the lens equation for each set of lens parameters
    >>> with Pool(2) as p:
    ...     result = p.map(solve_lens_equation1, input_arguments)
    >>> # result is a list of tuples
    >>> # each tuple contains the output parameters of the function
    >>> # each output parameter contains x_source, y_source, x0_image_position, x1_image_position, magnifications, time_delays, nImages, determinant, trace, iteration
    >>> print(f"magnification of images with lens parameters 'lens_parameters1' is {result[0][6]}")
    magnification of images with lens parameters 'lens_parameters1' is [ 2.18973765 -1.27542831]

    """
    n_min_images = int(lens_parameters[0])
    zl = lens_parameters[6]
    zs = lens_parameters[7]
    einstein_radius = lens_parameters[8]
    iteration = lens_parameters[9]
    # lensModel parameters are the same for the three functions used for image param calculation
    # 1. x-y position of images in the source plane, 2. magnifications, 3. time-delays (relative)
    lensModel = LensModel(
        lens_model_list=lens_parameters[10:].tolist(), z_lens=zl, z_source=zs
    )

    lens_eq_solver = LensEquationSolver(lensModel)

    factor = 1.0
    # ---------------------------------------------------#
    #     x-y position of images in the source plane
    # ---------------------------------------------------#
    # Get the caustic curve cut by the lens
    # First check if there is any nan in the caustic points
    while True:
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
        caustic_double_points = caustics_epl_shear(
            kwargs_lens, return_which="double", maginf=-100
        )
        caustic = np.logical_not(np.isnan(caustic_double_points).any())
        # If there is a nan, caustic=False, draw a new gamma
        if caustic:
            break
        else:
            lens_parameters[3] = np.random.normal(loc=2.0, scale=0.2, size=1)[0]
    caustic_double = Polygon(caustic_double_points.T)

    # check for strong lensed condition
    strongly_lensed = False
    while strongly_lensed == False:
        # Draw random points within the caustic
        # sometimes x_source, y_source positions are at same location and the solver fails
        # so we use a try-except block to catch the error and draw a new point
        try:
            x_source, y_source = pointpats.random.poisson(caustic_double, size=1)
            # Solve the lens equation
            (
                x0_image_position,
                x1_image_position,
            ) = lens_eq_solver.image_position_from_source(
                sourcePos_x=x_source,
                sourcePos_y=y_source,
                kwargs_lens=kwargs_lens,
                solver="analytical",
                magnification_limit=1.0 / 1000.0,
            )
            nImages = len(x0_image_position)  # shows how many images
            if nImages >= n_min_images:
                strongly_lensed = True
        except:
            pass

    # ---------------------------------------------------#
    #          magnification and time-delay
    # ---------------------------------------------------#
    # theta_E is in arcsec
    theta_E_nImages = einstein_radius * np.ones(nImages)
    radian_to_arcseconds = 180.0 / np.pi * 3600.0
    days_to_seconds = 24.0 * 3600.0
    # can have multiple magnification
    magnifications = lensModel.magnification(
        x0_image_position, x1_image_position, kwargs_lens
    )
    time_delays = (
        lensModel.arrival_time(x0_image_position, x1_image_position, kwargs_lens)
        * (theta_E_nImages * radian_to_arcseconds) ** 2
        * days_to_seconds
    )

    # ---------------------------------------------------#
    #     Params needed for image-type classification
    # ---------------------------------------------------#
    # it is faster to use numpy array operation to do image classification
    # return: f_xx, f_xy, f_yx, f_yy components
    hessian = lensModel.hessian(x0_image_position, x1_image_position, kwargs_lens)
    determinant = np.array(
        (1 - hessian[0]) * (1 - hessian[3]) - hessian[1] * hessian[2]
    )
    trace = np.array(2 - hessian[0] - hessian[3])

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
        lens_parameters[3],
        iteration,
    )
