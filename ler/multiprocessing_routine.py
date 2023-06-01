import numpy as np
# For caustic manipulation
from lenstronomy.LensModel.lens_model import LensModel
from lenstronomy.LensModel.Solver.lens_equation_solver import LensEquationSolver
from lenstronomy.LensModel.Solver.epl_shear_solver import caustics_epl_shear
from shapely.geometry import Polygon
import pointpats
from ler.helperroutines import cart2pol

from time import sleep

def solve_lens_equation1(lens_parameters):
    '''
    Function to solve the lens equation (min_image = 2)
    Input parameters:
        lens_parameters : a list of parameters
                        lens_parameters[0] = e1 : ellipticity 
                        lens_parameters[1] = e2 : ellipticity
                        lens_parameters[2] = gamma : power-law index
                        lens_parameters[3] = gamma1 : shear
                        lens_parameters[4] = gamma2 : shear
                        lens_parameters[5] = zl : redshift of the lens
                        lens_parameters[6] = zs : redshift of the source
                        lens_parameters[7] = einstein_radius : Einstein radius
                        lens_parameters[8] = iteration : iteration number
                        lens_parameters[9:] = lens_model_list : list of lens models
    Output parameters:
        x_source : x position of the source in the source plane
        y_source : y position of the source in the source plane
        eta : polar coordinate of the source in the source plane
        phi : polar coordinate of the source in the source plane
        x0_image_position : x position of the images in the source plane
        x1_image_position : y position of the images in the source plane
        magnifications : magnification of the images
        time_delays : time-delay of the images
        nImages : number of images
        determinant : determinant of the hessian matrix
        trace : trace of the hessian matrix
        iteration : iteration number
        weights : weights for the caustic
    '''
    n_min_images = int(lens_parameters[0])
    zl = lens_parameters[6]
    zs = lens_parameters[7]
    einstein_radius = lens_parameters[8]
    iteration = lens_parameters[9]
    # lensModel parameters are the same for the three functions used for image param calculation
    # 1. x-y position of images in the source plane, 2. magnifications, 3. time-delays (relative)
    lensModel = LensModel(lens_model_list = lens_parameters[10:].tolist(), 
                          z_lens = zl,
                          z_source = zs )

    lens_eq_solver = LensEquationSolver(lensModel)

    factor = 1.0
    #---------------------------------------------------#
    #     x-y position of images in the source plane
    #---------------------------------------------------#
    # Get the caustic curve cut by the lens
    # First check if there is any nan in the caustic points
    while True:
        kwargs_lens = [{'theta_E': factor, 'e1': lens_parameters[1], 'e2': lens_parameters[2], 'gamma': lens_parameters[3], \
                      'center_x': 0.0, 'center_y': 0.0}, {'gamma1': lens_parameters[4], 'gamma2': lens_parameters[5], 'ra_0': 0, 'dec_0':0}]
        caustic_double_points = caustics_epl_shear(kwargs_lens, return_which='double', maginf=-100)
        caustic = np.logical_not(np.isnan(caustic_double_points).any()) 
        # If there is a nan, caustic=False, draw a new gamma
        if caustic:
            break
        else:
            lens_parameters[3] = np.random.normal(loc = 2., scale = 0.2, size = 1)[0]
    caustic_double = Polygon(caustic_double_points.T)

    # check for strong lensed condition
    strongly_lensed = False
    while strongly_lensed == False:
        # Draw random points within the caustic
        x_source, y_source = pointpats.random.poisson(caustic_double, size=1)
        # Transform to polar coordinates
        eta, phi = cart2pol(x_source, y_source)
        # Solve the lens equation
        x0_image_position,x1_image_position = \
        lens_eq_solver.image_position_from_source( sourcePos_x = x_source, sourcePos_y = y_source, kwargs_lens = kwargs_lens, \
                                                  solver='analytical', magnification_limit = 1./100.)
        nImages = len(x0_image_position) # shows how many images
        if nImages >= n_min_images:
            strongly_lensed = True

    #---------------------------------------------------#
    #          magnification and time-delay
    #---------------------------------------------------#
    # theta_E is in arcsec
    theta_E_nImages = einstein_radius*np.ones(nImages)
    radian_to_arcseconds = 180.0/np.pi*3600.0
    days_to_seconds = 24.0*3600.0
    # can have multiple magnification
    magnifications = lensModel.magnification(x0_image_position, x1_image_position, kwargs_lens)
    time_delays = lensModel.arrival_time(x0_image_position, x1_image_position, kwargs_lens)*(theta_E_nImages*radian_to_arcseconds)**2*days_to_seconds

    #---------------------------------------------------#
    #     Params needed for image-type classification
    #---------------------------------------------------#
    # it is faster to use numpy array operation to do image classification
    # return: f_xx, f_xy, f_yx, f_yy components
    hessian = lensModel.hessian(x0_image_position, x1_image_position, kwargs_lens)
    determinant = np.array( (1 - hessian[0])*(1 - hessian[3]) - hessian[1]*hessian[2] )
    trace = np.array(2 - hessian[0] - hessian[3])
    weights = 1.0

    return(x_source, y_source, eta, phi, x0_image_position, x1_image_position, magnifications,time_delays, nImages , \
           determinant, trace, iteration, weights)

def solve_lens_equation2(lens_parameters):
    '''
    Function to solve the lens equation (min_image > 2)
    Input parameters:
        lens_parameters : a list of parameters
                        lens_parameters[0] = e1 : ellipticity 
                        lens_parameters[1] = e2 : ellipticity
                        lens_parameters[2] = gamma : power-law index
                        lens_parameters[3] = gamma1 : shear
                        lens_parameters[4] = gamma2 : shear
                        lens_parameters[5] = zl : redshift of the lens
                        lens_parameters[6] = zs : redshift of the source
                        lens_parameters[7] = einstein_radius : Einstein radius
                        lens_parameters[8] = iteration : iteration number
                        lens_parameters[9:] = lens_model_list : list of lens models
    Output parameters:
        x_source : x position of the source in the source plane
        y_source : y position of the source in the source plane
        eta : polar coordinate of the source in the source plane
        phi : polar coordinate of the source in the source plane
        x0_image_position : x position of the images in the source plane
        x1_image_position : y position of the images in the source plane
        magnifications : magnification of the images
        time_delays : time-delay of the images
        nImages : number of images
        determinant : determinant of the hessian matrix
        trace : trace of the hessian matrix
        iteration : iteration number
        weights : weights for the caustic
    '''
    n_min_images = int(lens_parameters[0])
    zl = lens_parameters[6]
    zs = lens_parameters[7]
    einstein_radius = lens_parameters[8]
    iteration = lens_parameters[9]
    # lensModel parameters are the same for the three functions used for image param calculation
    # 1. x-y position of images in the source plane, 2. magnifications, 3. time-delays (relative)
    lensModel = LensModel(lens_model_list = lens_parameters[10:].tolist(), 
                          z_lens = zl,
                          z_source = zs )

    lens_eq_solver = LensEquationSolver(lensModel)

    factor = 1.0
    #---------------------------------------------------#
    #     x-y position of images in the source plane
    #---------------------------------------------------#
    # Get the caustic curve cut by the lens
    # First check if there is any nan in the caustic points
    while True:
        kwargs_lens = [{'theta_E': factor, 'e1': lens_parameters[1], 'e2': lens_parameters[2], 'gamma': lens_parameters[3], \
                      'center_x': 0.0, 'center_y': 0.0}, {'gamma1': lens_parameters[4], 'gamma2': lens_parameters[5], 'ra_0': 0, 'dec_0':0}]
        caustic_double_points = caustics_epl_shear(kwargs_lens, return_which='double', maginf=-100)
        caustic_quad_points = caustics_epl_shear(kwargs_lens, return_which='caustic', maginf=-100)
        caustic = not(np.isnan(caustic_double_points).any()) 
        caustic &= not(np.isnan(caustic_quad_points).any())
        # If there is a nan, caustic=False, draw a new gamma
        if caustic:
            break
        else:
            #print('\n ERROR1 \n')
            lens_parameters[3] = np.random.normal(loc = 2., scale = 0.2, size = 1)[0]

    caustic_double = Polygon(caustic_double_points.T)
    caustic_quad = Polygon(caustic_quad_points.T)
    # check for strong lensed condition
    strongly_lensed = False
    while strongly_lensed == False:
        # Draw random points within the caustic
        try:
            x_source, y_source = pointpats.random.poisson(caustic_quad, size=1)
            # Solve the lens equation
            x0_image_position,x1_image_position = \
            lens_eq_solver.image_position_from_source( sourcePos_x = x_source, sourcePos_y = y_source, kwargs_lens = kwargs_lens, \
                                                      solver='analytical', magnification_limit = 1./100.)
            nImages = len(x0_image_position) # shows how many images
            if nImages >= n_min_images:
                strongly_lensed = True
        except:
            #print(lens_parameters[1], lens_parameters[2], lens_parameters[3], lens_parameters[4], lens_parameters[5], lens_parameters[5], lens_parameters[6])
            # print(x_source,y_source)
            x_source, y_source = pointpats.random.poisson(caustic_quad, size=1)
            # Solve the lens equation
            
            try:
                x0_image_position,x1_image_position = \
                lens_eq_solver.image_position_from_source( sourcePos_x = x_source, sourcePos_y = y_source, kwargs_lens = kwargs_lens, \
                                                          solver='analytical', magnification_limit = 1./100.)
                nImages = len(x0_image_position) # shows how many images
                if nImages >= n_min_images:
                    strongly_lensed = True
            except:
                print('\n ERROR in sampling from caustic \n')
                pass

    
    # Transform to polar coordinates
    eta, phi = cart2pol(x_source, y_source)
    #---------------------------------------------------#
    #          magnification and time-delay
    #---------------------------------------------------#
    # theta_E is in arcsec
    theta_E_nImages = einstein_radius*np.ones(nImages)
    radian_to_arcseconds = 180.0/np.pi*3600.0
    days_to_seconds = 24.0*3600.0
    # can have multiple magnification
    magnifications = lensModel.magnification(x0_image_position, x1_image_position, kwargs_lens)
    time_delays = lensModel.arrival_time(x0_image_position, x1_image_position, kwargs_lens)*(theta_E_nImages*radian_to_arcseconds)**2*days_to_seconds

    #---------------------------------------------------#
    #     Params needed for image-type classification
    #---------------------------------------------------#
    # it is faster to use numpy array operation to do image classification
    # return: f_xx, f_xy, f_yx, f_yy components
    hessian = lensModel.hessian(x0_image_position, x1_image_position, kwargs_lens)
    determinant = np.array( (1 - hessian[0])*(1 - hessian[3]) - hessian[1]*hessian[2] )
    trace = np.array(2 - hessian[0] - hessian[3])
    weights = caustic_quad.area/caustic_double.area

    return(x_source, y_source, eta, phi, x0_image_position, x1_image_position, magnifications,time_delays, nImages , \
           determinant, trace, iteration, weights)
