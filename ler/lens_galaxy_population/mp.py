import numpy as np
from scipy.integrate import quad
from numba import njit, prange
from lenstronomy.LensModel.Solver.epl_shear_solver import caustics_epl_shear
from shapely.geometry import Polygon

from .jit_functions import phi_cut_SIE, phi_q2_ellipticity_hemanta
from ..utils import inverse_transform_sampler, cubic_spline_interpolator, cubic_spline_interpolator2d_array, inverse_transform_sampler2d, load_json

from .cross_section_interpolator import make_cross_section_reinit

@njit(parallel=True)
def lens_redshift_strongly_lensed_njit(
        zs_array, # 1D
        zl_scaled, # 2D
        sigma_min,
        sigma_max,
        q_rvs,
        phi_rvs,
        gamma_rvs,
        shear_rvs,
        number_density,
        cross_section,
        dVcdz_function,
        integration_size,
    ):

    size_zs = zs_array.shape[0]
    size_zl = zl_scaled.shape[1]
    result_array = np.zeros((size_zs, size_zl))

    for i in prange(size_zs):

        # redshift of the lens between the source and the observer
        zl_array = zl_scaled[i]*zs_array[i] # 1D array, unscaled
        
        for j in range(size_zl):

            # monte carlo integration 
            zs_ = zs_array[i]*np.ones(integration_size)
            zl = zl_array[j]
            zl_ = zl*np.ones(integration_size)

            # velocity dispersion distribution
            sigma = np.random.uniform(sigma_min, sigma_max, integration_size)
            q = q_rvs(integration_size, sigma)
            phi = phi_rvs(integration_size)
            gamma = gamma_rvs(integration_size)
            gamma1, gamma2 = shear_rvs(integration_size)
            
            area_array = cross_section(
                zs_, zl_, sigma, q, phi, gamma, gamma1, gamma2
            )


            # # this is to avoid nan, inf, -inf
            # # idx = (area_array>=2.4906932608826035e-28) & (area_array<=1.972789211435565e-08)
            # idx = np.logical_not(np.isinf(area_array))
            # idx &= (area_array>0)
            # # idx = np.logical_not(np.isnan(area_array))
            idx = np.logical_not(np.isinf(area_array))
            idx &= (area_array>0)
            if idx.sum() == 0:
                result_array[i, j] = 0.
                continue

            phi_sigma = number_density(sigma, zl_)
                
            # diffrential co-moving volume 
            dVcdz = dVcdz_function(np.array([zl]))[0]

            # print('dVcdz: ', dVcdz)

            # cross-section
            result = (sigma_max-sigma_min)*np.average( area_array[idx] * phi_sigma[idx] * dVcdz )/(4*np.pi)
            result_array[i, j] = result

            # if (zs_array[i]==0.001) and (zl_array[j]==0.0001):
            #     print('zs: ', zs_array[i])
            #     print('zl: ', zl_array[j])
            #     print(f'sigma: {sigma}')
            #     print(f'q: {q}')
            #     print('area_array: ', area_array)
            #     print('phi_sigma: ', phi_sigma)
            #     print('dVcdz: ', dVcdz)
            #     print('result: ', result)
        
    return(result_array)

def lens_redshift_strongly_lensed_mp(params):

    size = params[11]  # integration size
    zs = params[0] # float
    # angular diameter distance
    Da_zs = cubic_spline_interpolator(np.array([zs]), params[4][0], params[4][1])[0] # float

    # redshift of the lens between the source and the observer
    zl_array = params[1]*zs # 1D array, unscaled

    # initialize interpolator
    if params[7][0] == 'cross_section_epl_shear_interpolation':
        cs_grid, cs_csunit_slope, cs_csunit_intercept, e1_grid, e2_grid, gamma_grid, gamma1_grid, gamma2_grid = load_json(params[7][1])

        # Create the interpolator instance
        cs_caculator = make_cross_section_reinit(
            e1_grid=e1_grid, 
            e2_grid=e2_grid, 
            gamma_grid=gamma_grid, 
            gamma1_grid=gamma1_grid, 
            gamma2_grid=gamma2_grid, 
            cs_grid=cs_grid,
            slope=0.31830988618379075,
            intercept=-3.2311742677852644e-27
        )

    result_array = []
    for i, zl in enumerate(zl_array):
        # monte carlo integration 
        # (sigma_max-sigma_min)*np.average( area_array[idx] * phi_sigma[idx] * dVcdz )/(4*np.pi)
        # normalization factor (sigma_max-sigma_min)/(4*np.pi)
        # cross-section part: area_array[idx]
        # number density part: phi_sigma[idx]
        # zs is fixed for the entire function call
        # zl is fixed for each iteration of the for loop

        # velocity dispersion distribution
        sigma_min, sigma_max = params[2][0], params[2][1] 
        sigma = np.random.uniform(sigma_min, sigma_max, size)

        ## Einstein radius calculation
        # angular diameter distance
        Da_zl = cubic_spline_interpolator(np.array([zl]), params[4][0], params[4][1])[0]
        Dls = (Da_zs*(1+zs) - Da_zl*(1+zl))/(1+zs)

        # einstein radius
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / Da_zs
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        #######
        # SIS #
        #######
        if params[7][0] == 'cross_section_sis':
            # # print('\nUsing SIS cross-section')
            area_array = np.pi*theta_E**2

            # # print('area_array: ', area_array)

        else: 
            # axis ratio distribution
            if params[3][0] == 'axis_ratio_uniform':
                q = np.random.uniform(params[3][1], params[3][2], size)
            else:
                if params[3][3] is not None:
                    q = inverse_transform_sampler2d(size, sigma, params[3][1], params[3][2], params[3][3])
                else:
                    q = inverse_transform_sampler(size, params[3][1], params[3][2])
            
            #######
            # SIE #
            #######
            if params[7][0] == 'cross_section_sie_feixu':
                area_array = phi_cut_SIE(q) * np.pi*theta_E**2

            else:
                # axis rotation angle
                if params[8][0] == 'axis_rotation_angle_uniform':
                    phi = np.random.uniform(params[8][1], params[8][2], size)
                else:
                    phi = inverse_transform_sampler(size, params[8][1], params[8][2])

                # Transform the axis ratio and the axis rotation angle, to ellipticities e1, e2, using lenstronomy
                e1, e2 = phi_q2_ellipticity_hemanta(phi, q)

                # Sample shears
                # gamma1, gamma2 = np.random.normal(loc=params[9][0], scale=params[9][1], size=(2,size))
                if params[9][0] == 'external_shear_normal':
                    gamma1, gamma2 = np.random.normal(loc=params[9][1], scale=params[9][2], size=(2,size))
                else:
                    raise ValueError('params[9][0], for gamma1 and gamma2, is not recognized. Only external_shear_normal is implemented.')


                # Sample the density profile slop of the mass density distribution
                # gamma = np.random.normal(loc=params[10][0], scale=params[10][1], size=size)
                if params[10][0] == 'density_profile_slope_normal':
                    gamma = np.random.normal(loc=params[10][1], scale=params[10][2], size=size)
                else:
                    gamma = inverse_transform_sampler(size, params[10][1], params[10][2])

                #############
                # EPL+Shear #
                #############
                # cross-section
                # option 1: use nearest-neighbor interpolation (fast)
                # option 2: direct cross-section calculation with lenstronomy
                if params[7][0] == 'cross_section_epl_shear_interpolation':
                    area_array = cs_caculator(theta_E, q, phi, gamma, gamma1, gamma2)
                elif params[7][0] == 'cross_section_epl_shear_numerical':
                    area_array = []
                    for i in range(size):
                        area_array.append(cross_section(theta_E[i], e1[i], e2[i], gamma[i], gamma1[i], gamma2[i]))
                    area_array = np.array(area_array)
                else:
                    raise(ValueError('params[7][0] not recognized. Should be cross_section_epl_shear_interpolation, cross_section_epl_shear_numerical, cross_section_sis, or cross_section_sie_feixu.'))

        # # this is to avoid nan, inf, -inf
        # # idx = (area_array>=2.4906932608826035e-28) & (area_array<=1.972789211435565e-08)
        # idx = np.logical_not(np.isinf(area_array))
        # idx &= (area_array>0)
        # # idx = np.logical_not(np.isnan(area_array))
        idx = np.logical_not(np.isinf(area_array))
        idx &= (area_array>0)
        if idx.sum() == 0:
            result_array.append(0.)
            continue

        # velocity dispersion distribution
        zl_ = zl*np.ones(size)
        if params[2][3] is None:
            phi_sigma = cubic_spline_interpolator(
                sigma,
                params[2][4],
                params[2][2]
            )   
        else:
            phi_sigma = cubic_spline_interpolator2d_array(
                sigma, 
                zl_,
                params[2][4],
                params[2][2],
                params[2][3]
            )
            
        # diffrential co-moving volume 
        dVcdz = cubic_spline_interpolator(np.array([zl]), params[5][0], params[5][1])[0]

        # cross-section
        result = (sigma_max-sigma_min)*np.average( area_array[idx] * phi_sigma[idx] * dVcdz )/(4*np.pi)
        result_array.append(result)

        if (zs==0.001) and (zl==0.0001):
            print('zs: ', zs)
            print('zl: ', zl)
            print(f'sigma: {sigma}')
            print(f'q: {q}')
            print('area_array: ', area_array)
            print('phi_sigma: ', phi_sigma)
            print('dVcdz: ', dVcdz)
            print('result: ', result)
            print('theta_E: ', theta_E)
    
    return(params[6], np.array(result_array))

# def lens_redshift_strongly_lensed_mp(params):
#     """
#     input_params = np.array(
#         [
#             zs,
#             zl_scaled,
#             sigma_args,
#             q_args,
#             Da_args,
#             dVcdz_args,
#             idx,
#             cs_args,
#             phi_args,
#             shear_args,
#             slope_args,
#             integration_size,
#         ],
#         dtype=object,
#     )
#     """

#     size = params[11]  # integration size
#     zs = params[0] # float
#     # angular diameter distance
#     Da_function = njit(lambda z: cubic_spline_interpolator(z, params[4][0], params[4][1]))
#     Da_zs = Da_function(np.array([zs]))[0]

#     # redshift of the lens between the source and the observer
#     zl_array = params[1]*zs # 1D array, unscaled

#     # initialize interpolator
#     if params[7][0] == 'cross_section_epl_shear_interpolation':
#         (
#             e1_grid,
#             e2_grid,
#             gamma_grid,
#             gamma1_grid,
#             gamma2_grid,
#             cs_spline_coeff_grid,
#             _, # cs_csunit_slope
#             _, # cs_csunit_intercept
#         ) = load_json(params[7][1])

#         # Create the interpolator instance
#         cs_caculator = make_cross_section_reinit(
#             e1_grid=e1_grid, 
#             e2_grid=e2_grid, 
#             gamma_grid=gamma_grid, 
#             gamma1_grid=gamma1_grid, 
#             gamma2_grid=gamma2_grid, 
#             cs_spline_coeff_grid=cs_spline_coeff_grid,
#             slope=0.31830988618379075,
#             intercept=-3.2311742677852644e-27
#         )

#     result_array = []
#     for i, zl in enumerate(zl_array):
#         # monte carlo integration 
#         # (sigma_max-sigma_min)*np.average( area_array[idx] * phi_sigma[idx] * dVcdz )/(4*np.pi)
#         # normalization factor (sigma_max-sigma_min)/(4*np.pi)
#         # cross-section part: area_array[idx]
#         # number density part: phi_sigma[idx]
#         # zs is fixed for the entire function call
#         # zl is fixed for each iteration of the for loop

#         # velocity dispersion distribution
#         sigma_min, sigma_max = params[2][0], params[2][1] 
#         sigma = np.random.uniform(sigma_min, sigma_max, size)

#         # print('sigma: ', sigma)

#         #######
#         # SIS #
#         #######
#         if params[7][0] == 'cross_section_sis':
#             # # print('\nUsing SIS cross-section')
#             area_array = np.pi*theta_E**2

#             # # print('area_array: ', area_array)

#         else: 
#             # axis ratio distribution
#             if params[3][0] == 'axis_ratio_uniform':
#                 q = np.random.uniform(params[3][1], params[3][2], size)
#             else:
#                 if params[3][3] is not None:
#                     q = inverse_transform_sampler2d(size, sigma, params[3][1], params[3][2], params[3][3])
#                 else:
#                     q = inverse_transform_sampler(size, params[3][1], params[3][2])
            
#             #######
#             # SIE #
#             #######
#             if params[7][0] == 'cross_section_sie_feixu':
#                 area_array = phi_cut_SIE(q) * np.pi*theta_E**2

#             else:
#                 # axis rotation angle
#                 if params[8][0] == 'axis_rotation_angle_uniform':
#                     phi = np.random.uniform(params[8][1], params[8][2], size)
#                 else:
#                     phi = inverse_transform_sampler(size, params[8][1], params[8][2])

#                 # Transform the axis ratio and the axis rotation angle, to ellipticities e1, e2, using lenstronomy
#                 e1, e2 = phi_q2_ellipticity_hemanta(phi, q)

#                 # Sample shears
#                 # gamma1, gamma2 = np.random.normal(loc=params[9][0], scale=params[9][1], size=(2,size))
#                 if params[9][0] == 'external_shear_normal':
#                     gamma1, gamma2 = np.random.normal(loc=params[9][1], scale=params[9][2], size=(2,size))
#                 else:
#                     raise ValueError('params[9][0], for gamma1 and gamma2, is not recognized. Only external_shear_normal is implemented.')


#                 # Sample the density profile slop of the mass density distribution
#                 # gamma = np.random.normal(loc=params[10][0], scale=params[10][1], size=size)
#                 if params[10][0] == 'density_profile_slope_normal':
#                     gamma = np.random.normal(loc=params[10][1], scale=params[10][2], size=size)
#                 else:
#                     gamma = inverse_transform_sampler(size, params[10][1], params[10][2])

#                 #############
#                 # EPL+Shear #
#                 #############
#                 # cross-section
#                 if params[7][0] == 'cross_section_epl_shear_interpolation':
#                     area_array = cs_caculator(theta_E, q, phi, gamma, gamma1, gamma2)






#                 elif params[7][0] == 'cross_section_epl_shear_numerical':

#                     area_array = []
#                     for i in range(size):
#                         area_array.append(cross_section(theta_E[i], e1[i], e2[i], gamma[i], gamma1[i], gamma2[i]))
#                     area_array = np.array(area_array)
#                 else:
#                     raise(ValueError('params[7][0] not recognized. Should be cross_section_epl_shear_interpolation, cross_section_epl_shear_numerical, cross_section_sis, or cross_section_sie_feixu.'))

#         # # this is to avoid nan, inf, -inf
#         # # idx = (area_array>=2.4906932608826035e-28) & (area_array<=1.972789211435565e-08)
#         # idx = np.logical_not(np.isinf(area_array))
#         # idx &= (area_array>0)
#         # # idx = np.logical_not(np.isnan(area_array))
#         idx = np.logical_not(np.isinf(area_array))
#         idx &= (area_array>0)
#         if idx.sum() == 0:
#             result_array.append(0.)
#             continue

#         # velocity dispersion distribution
#         zl_ = zl*np.ones(size)
#         if params[2][3] is None:
#             phi_sigma = cubic_spline_interpolator(
#                 sigma,
#                 params[2][4],
#                 params[2][2]
#             )   
#         else:
#             phi_sigma = cubic_spline_interpolator2d_array(
#                 sigma, 
#                 zl_,
#                 params[2][4],
#                 params[2][2],
#                 params[2][3]
#             )

#         # print('phi_sigma: ', phi_sigma[idx])
            
#         # diffrential co-moving volume 
#         dVcdz = cubic_spline_interpolator(np.array([zl]), params[5][0], params[5][1])[0]

#         # print('dVcdz: ', dVcdz)

#         # cross-section
#         result = (sigma_max-sigma_min)*np.average( area_array[idx] * phi_sigma[idx] * dVcdz )/(4*np.pi)
#         result_array.append(result)

#         # print('result: ', result)
    
#     return(params[6], np.array(result_array))

def cross_section_unit_mp(params):

    kwargs_lens = [
        {
            "theta_E": 1.0,
            "e1": params[0],
            "e2": params[1],
            "gamma": params[2],
            "center_x": 0.0,
            "center_y": 0.0,
        },
        {
            "gamma1": params[3],
            "gamma2": params[4],
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
        area_ = (Polygon(caustic_double_points.T).area)
    else:
        area_ = 0.0
        
    return params[5], area_

def cross_section_mp(params):

    theta_E = params[0]
    e1 = params[1]
    e2 = params[2]
    gamma = params[3]
    gamma1 = params[4]
    gamma2 = params[5]
    idx = params[6]

    kwargs_lens = [
        {
            "theta_E": theta_E,
            "e1": e1,
            "e2": e2,
            "gamma": gamma,
            "center_x": 0.0,
            "center_y": 0.0,
        },
        {
            "gamma1": gamma1,
            "gamma2": gamma2,
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
        area = Polygon(caustic_double_points.T).area
    else:
        area = 0.0

    return int(idx) , area

# def cross_section_interpolated_mp(params):

#     # theta_E = params[0]
#     # e1 = np.array(params[1])
#     # e2 = np.array(params[2])
#     # gamma = np.array(params[3])
#     # gamma1 = np.array(params[4])
#     # gamma2 = np.array(params[5])
#     # idx = params[6]
#     # nbrs = params[7]
#     # values = params[8]
#     # cross_section_spline = params[9]
#     # sis_area_array = params[10]
#     theta_E, e1, e2, gamma, gamma1, gamma2, idx, nbrs, values, cross_section_spline, sis_area_array = params
    
#     points_ = np.array([e1, e2, gamma, gamma1, gamma2]).T
#     _, indices = nbrs.kneighbors(points_)
#     new_values_nn1 = np.mean(values[indices], axis=1)

#     theta_E_correction = cubic_spline_interpolator(np.pi*theta_E**2, cross_section_spline, sis_area_array)
#     area = new_values_nn1*theta_E_correction

#     return idx, area

def cross_section(theta_E, e1, e2, gamma, gamma1, gamma2):

    kwargs_lens = [
        {
            "theta_E": theta_E,
            "e1": e1,
            "e2": e2,
            "gamma": gamma,
            "center_x": 0.0,
            "center_y": 0.0,
        },
        {
            "gamma1": gamma1,
            "gamma2": gamma2,
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
        area = Polygon(caustic_double_points.T).area
    else:
        area = 0.0

    return area