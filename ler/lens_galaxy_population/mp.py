import numpy as np
from scipy.integrate import quad
from numba import njit
from lenstronomy.LensModel.Solver.epl_shear_solver import caustics_epl_shear
from shapely.geometry import Polygon

from .jit_functions import phi_cut_SIE, phi_q2_ellipticity_hemanta
from ..utils import inverse_transform_sampler, cubic_spline_interpolator, cubic_spline_interpolator2d_array, inverse_transform_sampler2d

# for testing
from .jit_functions import phi_loc_bernardi, phi
# phi(s, z, alpha, beta, phistar, sigmastar)


def lens_redshift_sis1_mp(params):

    size = 20000  # integration size
    zs = params[0] # float

    # redshift of the lens between the source and the observer
    zl_array = params[1]*zs # 1D array, unscaled

    result_array = []
    for i, zl in enumerate(zl_array):
        # velocity dispersion distribution
        sigma_min, sigma_max = params[2][0], params[2][1] # min_max[0]
        sigma = np.random.uniform(sigma_min, sigma_max, size)
    
        ## Einstein radius calculation
        # angular diameter distance
        Da_zs = cubic_spline_interpolator(np.array([zs]), params[4][0], params[4][1])[0] # float
        Da_zl = cubic_spline_interpolator(np.array([zl]), params[4][0], params[4][1])[0]

        # einstein radius 
        Dls = (Da_zs*(1+zs) - Da_zl*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / Da_zs
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        # cross-section
        area_array = np.pi*theta_E**2
        # non-nan index
        idx = np.logical_not(np.isnan(area_array))

        # velocity dispersion distribution
        phi_sigma = phi_loc_bernardi(
            sigma=sigma,
            alpha=params[2][2], 
            beta=params[2][3],
            phistar=params[2][4],
            sigmastar=params[2][5],
        )

        # diffrential co-moving volume 
        dVcdz = cubic_spline_interpolator(np.array([zl]), params[5][0], params[5][1])[0]

        # cross-section
        result_array.append((sigma_max-sigma_min)*np.average( area_array[idx] * phi_sigma[idx] * dVcdz )/(4*np.pi))
    
    return(params[6], np.array(result_array))

def lens_redshift_sis2_mp(params):

    size = 20000  # integration size
    zs = params[0] # float

    # redshift of the lens between the source and the observer
    zl_array = params[1]*zs # 1D array, unscaled

    result_array = []
    for i, zl in enumerate(zl_array):
        # velocity dispersion distribution
        sigma_min, sigma_max = params[2][0], params[2][1] # min_max[0]
        sigma = np.random.uniform(sigma_min, sigma_max, size)
    
        ## Einstein radius calculation
        # angular diameter distance
        Da_zs = cubic_spline_interpolator(np.array([zs]), params[4][0], params[4][1])[0] # float
        Da_zl = cubic_spline_interpolator(np.array([zl]), params[4][0], params[4][1])[0]

        # einstein radius 
        Dls = (Da_zs*(1+zs) - Da_zl*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / Da_zs
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        # cross-section
        area_array = np.pi*theta_E**2
        # non-nan index
        idx = np.logical_not(np.isnan(area_array))

        # velocity dispersion distribution
        zl_ = zl*np.ones(size)
        phi_sigma = phi(
            s=sigma,
            z=zl_, 
            alpha=params[2][2], 
            beta=params[2][3],
            phistar=params[2][4],
            sigmastar=params[2][5],
        )

        # diffrential co-moving volume 
        dVcdz = cubic_spline_interpolator(np.array([zl]), params[5][0], params[5][1])[0]

        # cross-section
        result_array.append((sigma_max-sigma_min)*np.average( area_array[idx] * phi_sigma[idx] * dVcdz )/(4*np.pi))
    
    return(params[6], np.array(result_array))

def lens_redshift_sie1_mp(params):

    size = 20000  # integration size
    zs = params[0] # float

    # redshift of the lens between the source and the observer
    zl_array = params[1]*zs # 1D array, unscaled

    result_array = []
    for i, zl in enumerate(zl_array):
        # velocity dispersion distribution
        sigma_min, sigma_max = params[2][0], params[2][1] # min_max[0]
        sigma = np.random.uniform(sigma_min, sigma_max, size)

        # axis ratio distribution
        if len(params[3]) == 3:
            print
            q = inverse_transform_sampler2d(size, sigma, params[3][0], params[3][1], params[3][2])
        else:
            if isinstance(params[3][0], float):
                q = np.random.uniform(params[3][0], params[3][1], size)
            else:
                q = inverse_transform_sampler(size, params[3][0], params[3][1])
    
        ## Einstein radius calculation
        # angular diameter distance
        Da_zs = cubic_spline_interpolator(np.array([zs]), params[4][0], params[4][1])[0] # float
        Da_zl = cubic_spline_interpolator(np.array([zl]), params[4][0], params[4][1])[0]

        # einstein radius 
        Dls = (Da_zs*(1+zs) - Da_zl*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / Da_zs
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        # cross-section
        area_array = phi_cut_SIE(q) * np.pi*theta_E**2
        # area_array_sis = np.pi*theta_E**2
        # logic = area_array_sis < area_array
        # # area_array should be less than area_array_sis
        # if logic.any():
        #     raise(ValueError('area_array should be less than area_array_sis'))

        # non-nan index
        idx = np.logical_not(np.isnan(area_array))

        # velocity dispersion distribution
        phi_sigma = phi_loc_bernardi(
            sigma=sigma,
            alpha=params[2][2], 
            beta=params[2][3],
            phistar=params[2][4],
            sigmastar=params[2][5],
        )

        # diffrential co-moving volume 
        dVcdz = cubic_spline_interpolator(np.array([zl]), params[5][0], params[5][1])[0]

        # cross-section
        result_array.append((sigma_max-sigma_min)*np.average( area_array[idx] * phi_sigma[idx] * dVcdz )/(4*np.pi))
    
    return(params[6], np.array(result_array))


def lens_redshift_sie2_mp(params):

    size = 20000  # integration size
    zs = params[0] # float

    # redshift of the lens between the source and the observer
    zl_array = params[1]*zs # 1D array, unscaled

    result_array = []
    for i, zl in enumerate(zl_array):
        # velocity dispersion distribution
        sigma_min, sigma_max = params[2][0], params[2][1] # min_max[0]
        sigma = np.random.uniform(sigma_min, sigma_max, size)

        # axis ratio distribution
        if len(params[3]) == 3:
            q = inverse_transform_sampler2d(size, sigma, params[3][0], params[3][1], params[3][2])
        else:
            if isinstance(params[3][0], float):
                q = np.random.uniform(params[3][0], params[3][1], size)
            else:
                q = inverse_transform_sampler(size, params[3][0], params[3][1])
    
        ## Einstein radius calculation
        # angular diameter distance
        Da_zs = cubic_spline_interpolator(np.array([zs]), params[4][0], params[4][1])[0] # float
        Da_zl = cubic_spline_interpolator(np.array([zl]), params[4][0], params[4][1])[0]

        # einstein radius 
        Dls = (Da_zs*(1+zs) - Da_zl*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / Da_zs
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        # cross-section
        area_array = phi_cut_SIE(q) * np.pi*theta_E**2
        # non-nan index
        idx = np.logical_not(np.isnan(area_array))

        # velocity dispersion distribution
        zl_ = zl*np.ones(size)
        phi_sigma = phi(
            s=sigma,
            z=zl_, 
            alpha=params[2][2], 
            beta=params[2][3],
            phistar=params[2][4],
            sigmastar=params[2][5],
        )

        # diffrential co-moving volume 
        dVcdz = cubic_spline_interpolator(np.array([zl]), params[5][0], params[5][1])[0]

        # cross-section
        result_array.append((sigma_max-sigma_min)*np.average( area_array[idx] * phi_sigma[idx] * dVcdz )/(4*np.pi))
    
    return(params[6], np.array(result_array))


def lens_redshift_sie3_mp(params):

    size = 20000  # integration size
    zs = params[0] # float

    # redshift of the lens between the source and the observer
    zl_array = params[1]*zs # 1D array, unscaled

    result_array = []
    for i, zl in enumerate(zl_array):
        # velocity dispersion distribution
        sigma_min, sigma_max = params[2][0], params[2][1] # min_max[0]
        sigma = np.random.uniform(sigma_min, sigma_max, size)

        # axis ratio distribution
        if len(params[3]) == 3:
            q = inverse_transform_sampler2d(size, sigma, params[3][0], params[3][1], params[3][2])
        else:
            if isinstance(params[3][0], float):
                q = np.random.uniform(params[3][0], params[3][1], size)
            else:
                q = inverse_transform_sampler(size, params[3][0], params[3][1])

        # axis rotation angle
        psi = np.random.uniform(params[8][0], params[8][1], size)

        # Transform the axis ratio and the angle, to ellipticities e1, e2, using lenstronomy
        e1, e2 = phi_q2_ellipticity_hemanta(psi, q)

        # Sample shears
        gamma1, gamma2 = np.zeros((2,size))

        # Sample the density profile slop of the mass density distribution
        gamma = np.ones(size)*2.0
    
        ## Einstein radius calculation
        # angular diameter distance
        Da_zs = cubic_spline_interpolator(np.array([zs]), params[4][0], params[4][1])[0] # float
        Da_zl = cubic_spline_interpolator(np.array([zl]), params[4][0], params[4][1])[0]

        # einstein radius 
        Dls = (Da_zs*(1+zs) - Da_zl*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / Da_zs
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        # cross-section
        nbrs = params[7][0]
        values = params[7][1]
        points_ = np.array([e1, e2, gamma, gamma1, gamma2]).T
        _, indices = nbrs.kneighbors(points_)
        new_values_nn1 = np.mean(values[indices], axis=1)
        theta_E_correction = cubic_spline_interpolator(np.pi*theta_E**2, params[7][2], params[7][3])
        area_array = new_values_nn1*theta_E_correction
        idx = (area_array>=2.4906932608826035e-28) & (area_array<=1.972789211435565e-08)

        # velocity dispersion distribution
        phi_sigma = phi_loc_bernardi(
            sigma=sigma,
            alpha=params[2][2], 
            beta=params[2][3],
            phistar=params[2][4],
            sigmastar=params[2][5],
        )

        # diffrential co-moving volume 
        dVcdz = cubic_spline_interpolator(np.array([zl]), params[5][0], params[5][1])[0]

        # cross-section
        result_array.append((sigma_max-sigma_min)*np.average( area_array[idx] * phi_sigma[idx] * dVcdz )/(4*np.pi))
    
    return(params[6], np.array(result_array))


def lens_redshift_sie4_mp(params):

    size = 20000  # integration size
    zs = params[0] # float

    # redshift of the lens between the source and the observer
    zl_array = params[1]*zs # 1D array, unscaled

    result_array = []
    for i, zl in enumerate(zl_array):
        # velocity dispersion distribution
        sigma_min, sigma_max = params[2][0], params[2][1] # min_max[0]
        sigma = np.random.uniform(sigma_min, sigma_max, size)

        # axis ratio distribution
        if len(params[3]) == 3:
            q = inverse_transform_sampler2d(size, sigma, params[3][0], params[3][1], params[3][2])
        else:
            if isinstance(params[3][0], float):
                q = np.random.uniform(params[3][0], params[3][1], size)
            else:
                q = inverse_transform_sampler(size, params[3][0], params[3][1])

        # axis rotation angle
        psi = np.random.uniform(params[8][0], params[8][1], size)

        # Transform the axis ratio and the angle, to ellipticities e1, e2, using lenstronomy
        e1, e2 = phi_q2_ellipticity_hemanta(psi, q)

        # Sample shears
        gamma1, gamma2 = np.zeros((2,size))

        # Sample the density profile slop of the mass density distribution
        gamma = np.ones(size)*2.0
    
        ## Einstein radius calculation
        # angular diameter distance
        Da_zs = cubic_spline_interpolator(np.array([zs]), params[4][0], params[4][1])[0] # float
        Da_zl = cubic_spline_interpolator(np.array([zl]), params[4][0], params[4][1])[0]

        # einstein radius 
        Dls = (Da_zs*(1+zs) - Da_zl*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / Da_zs
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        # cross-section
        nbrs = params[7][0]
        values = params[7][1]
        points_ = np.array([e1, e2, gamma, gamma1, gamma2]).T
        _, indices = nbrs.kneighbors(points_)
        new_values_nn1 = np.mean(values[indices], axis=1)
        theta_E_correction = cubic_spline_interpolator(np.pi*theta_E**2, params[7][2], params[7][3])
        area_array = new_values_nn1*theta_E_correction
        idx = (area_array>=2.4906932608826035e-28) & (area_array<=1.972789211435565e-08)

        # velocity dispersion distribution
        zl_ = zl*np.ones(size)
        phi_sigma = phi(
            s=sigma,
            z=zl_, 
            alpha=params[2][2], 
            beta=params[2][3],
            phistar=params[2][4],
            sigmastar=params[2][5],
        )

        # diffrential co-moving volume 
        dVcdz = cubic_spline_interpolator(np.array([zl]), params[5][0], params[5][1])[0]

        # cross-section
        result_array.append((sigma_max-sigma_min)*np.average( area_array[idx] * phi_sigma[idx] * dVcdz )/(4*np.pi))
    
    return(params[6], np.array(result_array))


def lens_redshift_epl_shear1_mp(params):

    size = 20000  # integration size
    zs = params[0] # float

    # redshift of the lens between the source and the observer
    zl_array = params[1]*zs # 1D array, unscaled

    result_array = []
    for i, zl in enumerate(zl_array):
        # velocity dispersion distribution
        sigma_min, sigma_max = params[2][0], params[2][1] # min_max[0]
        sigma = np.random.uniform(sigma_min, sigma_max, size)

        # axis ratio distribution
        if len(params[3]) == 3:
            q = inverse_transform_sampler2d(size, sigma, params[3][0], params[3][1], params[3][2])
        else:
            if isinstance(params[3][0], float):
                q = np.random.uniform(params[3][0], params[3][1], size)
            else:
                q = inverse_transform_sampler(size, params[3][0], params[3][1])

        # axis rotation angle
        psi = np.random.uniform(params[8][0], params[8][1], size)

        # Transform the axis ratio and the angle, to ellipticities e1, e2, using lenstronomy
        e1, e2 = phi_q2_ellipticity_hemanta(psi, q)

        # Sample shears
        gamma1, gamma2 = np.random.normal(loc=params[9][0], scale=params[9][1], size=(2,size))

        # Sample the density profile slop of the mass density distribution
        gamma = np.random.normal(loc=params[10][0], scale=params[10][1], size=size)
    
        ## Einstein radius calculation
        # angular diameter distance
        Da_zs = cubic_spline_interpolator(np.array([zs]), params[4][0], params[4][1])[0] # float
        Da_zl = cubic_spline_interpolator(np.array([zl]), params[4][0], params[4][1])[0]

        # einstein radius 
        Dls = (Da_zs*(1+zs) - Da_zl*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / Da_zs
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        # cross-section
        nbrs = params[7][0]
        if nbrs is not None:
            values = params[7][1]
            points_ = np.array([e1, e2, gamma, gamma1, gamma2]).T
            _, indices = nbrs.kneighbors(points_)
            new_values_nn1 = np.mean(values[indices], axis=1)
            theta_E_correction = cubic_spline_interpolator(np.pi*theta_E**2, params[7][2], params[7][3])
            area_array = new_values_nn1*theta_E_correction
        else:
            area_array = []
            for i in range(size):
                area_array.append(cross_section(theta_E[i], e1[i], e2[i], gamma[i], gamma1[i], gamma2[i]))
            area_array = np.array(area_array)
        idx = (area_array>=2.4906932608826035e-28) & (area_array<=1.972789211435565e-08)

        # velocity dispersion distribution
        phi_sigma = phi_loc_bernardi(
            sigma=sigma,
            alpha=params[2][2], 
            beta=params[2][3],
            phistar=params[2][4],
            sigmastar=params[2][5],
        )

        # diffrential co-moving volume 
        dVcdz = cubic_spline_interpolator(np.array([zl]), params[5][0], params[5][1])[0]

        # cross-section
        result_array.append((sigma_max-sigma_min)*np.average( area_array[idx] * phi_sigma[idx] * dVcdz )/(4*np.pi))
    
    return(params[6], np.array(result_array))


def lens_redshift_epl_shear2_mp(params):


    size = 20000  # integration size
    zs = params[0] # float

    # redshift of the lens between the source and the observer
    zl_array = params[1]*zs # 1D array, unscaled

    result_array = []
    for i, zl in enumerate(zl_array):
        # velocity dispersion distribution
        sigma_min, sigma_max = params[2][0], params[2][1] # min_max[0]
        sigma = np.random.uniform(sigma_min, sigma_max, size)

        # axis ratio distribution
        if len(params[3]) == 3:
            q = inverse_transform_sampler2d(size, sigma, params[3][0], params[3][1], params[3][2])
        else:
            if isinstance(params[3][0], float):
                q = np.random.uniform(params[3][0], params[3][1], size)
            else:
                q = inverse_transform_sampler(size, params[3][0], params[3][1])

        # axis rotation angle
        psi = np.random.uniform(params[8][0], params[8][1], size)

        # Transform the axis ratio and the axis rotation angle, to ellipticities e1, e2, using lenstronomy
        e1, e2 = phi_q2_ellipticity_hemanta(psi, q)

        # Sample shears
        gamma1, gamma2 = np.random.normal(loc=params[9][0], scale=params[9][1], size=(2,size))

        # Sample the density profile slop of the mass density distribution
        gamma = np.random.normal(loc=params[10][0], scale=params[10][1], size=size)
    
        ## Einstein radius calculation
        # angular diameter distance
        Da_zs = cubic_spline_interpolator(np.array([zs]), params[4][0], params[4][1])[0] # float
        Da_zl = cubic_spline_interpolator(np.array([zl]), params[4][0], params[4][1])[0]

        # einstein radius 
        Dls = (Da_zs*(1+zs) - Da_zl*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / Da_zs
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        # cross-section
        nbrs = params[7][0]
        if nbrs is not None:
            values = params[7][1]
            points_ = np.array([e1, e2, gamma, gamma1, gamma2]).T
            _, indices = nbrs.kneighbors(points_)
            new_values_nn1 = np.mean(values[indices], axis=1)
            theta_E_correction = cubic_spline_interpolator(np.pi*theta_E**2, params[7][2], params[7][3])
            area_array = new_values_nn1*theta_E_correction
        else:
            area_array = []
            for i in range(size):
                area_array.append(cross_section(theta_E[i], e1[i], e2[i], gamma[i], gamma1[i], gamma2[i]))
            area_array = np.array(area_array)

        # this is to avoid nan, inf, -inf
        idx = (area_array>=2.4906932608826035e-28) & (area_array<=1.972789211435565e-08)

        # velocity dispersion distribution
        zl_ = zl*np.ones(size)
        phi_sigma = phi(
            s=sigma,
            z=zl_, 
            alpha=params[2][2], 
            beta=params[2][3],
            phistar=params[2][4],
            sigmastar=params[2][5],
        )

        # diffrential co-moving volume 
        dVcdz = cubic_spline_interpolator(np.array([zl]), params[5][0], params[5][1])[0]

        # cross-section
        result_array.append((sigma_max-sigma_min)*np.average( area_array[idx] * phi_sigma[idx] * dVcdz )/(4*np.pi))
    
    return(params[6], np.array(result_array))

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

def cross_section_interpolated_mp(params):

    # theta_E = params[0]
    # e1 = np.array(params[1])
    # e2 = np.array(params[2])
    # gamma = np.array(params[3])
    # gamma1 = np.array(params[4])
    # gamma2 = np.array(params[5])
    # idx = params[6]
    # nbrs = params[7]
    # values = params[8]
    # cross_section_spline = params[9]
    # sis_area_array = params[10]
    theta_E, e1, e2, gamma, gamma1, gamma2, idx, nbrs, values, cross_section_spline, sis_area_array = params
    
    points_ = np.array([e1, e2, gamma, gamma1, gamma2]).T
    _, indices = nbrs.kneighbors(points_)
    new_values_nn1 = np.mean(values[indices], axis=1)

    theta_E_correction = cubic_spline_interpolator(np.pi*theta_E**2, cross_section_spline, sis_area_array)
    area = new_values_nn1*theta_E_correction

    return idx, area

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

def rjs_sie_mp(params):

    size = 1000  # integration size
    zs = params[0] # float

    # redshift of the lens between the source and the observer
    zl_array = params[1]

    while True:
        # velocity dispersion distribution
        # njit might give problem
        if len(params[2]) == 2:
            sigma = inverse_transform_sampler2d(size, zl_array, params[2][0], params[2][1], params[2][2])
        else:
            sigma = inverse_transform_sampler(size, params[2][0], params[2][1])

        # axis ratio distribution
        if len(params[3]) == 3:
            q = inverse_transform_sampler2d(size, sigma, params[3][0], params[3][1], params[3][2])
        else:
            q = inverse_transform_sampler(size, params[3][0], params[3][1])
    
        ## Einstein radius calculation
        # angular diameter distance
        Da_zs = cubic_spline_interpolator(np.array([zs]), params[4][0], params[4][1])[0] # float
        Da_zl = cubic_spline_interpolator(np.array([zl]), params[4][0], params[4][1])[0]

        # einstein radius 
        Dls = (Da_zs*(1+zs) - Da_zl*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / Da_zs
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        # cross-section
        cross_section = phi_cut_SIE(q) * np.pi*theta_E**2

        # rejection sampling
        max_ = np.max(cross_section)  # maximum einstein radius
        u = np.random.uniform(0, max_, size=size)
        mask = u < cross_section
        if np.sum(mask) > 0:
            break

    return(params[6], sigma[mask][0], q[mask][0], theta_E[mask][0])

def rjs_sie_mp(args):

    idx, zs, zl, sigma_args_cdf_values, sigma_args_x_array, sigma_args_conditioned_y_array, q_args_cdf_values, q_args_x_array, q_args_conditioned_y_array, da_args_function_spline, da_args_x_array, buffer_size = args

    zl_ = zl*np.ones(buffer_size) # need for sigma computation
    max_value = 0.
    while True:
        # velocity dispersion distribution
        if sigma_args_conditioned_y_array is not None:
            sigma = inverse_transform_sampler2d(buffer_size, zl_, sigma_args_cdf_values, sigma_args_x_array, sigma_args_conditioned_y_array)
        else:
            sigma = inverse_transform_sampler(buffer_size, zl_, sigma_args_cdf_values, sigma_args_x_array)

        # axis ratio distribution
        if q_args_conditioned_y_array is not None:
            q = inverse_transform_sampler2d(buffer_size, sigma, q_args_cdf_values, q_args_x_array, q_args_conditioned_y_array)
        else:
            q = inverse_transform_sampler(buffer_size, sigma, q_args_cdf_values, q_args_x_array)
    
        ## Einstein radius calculation
        # angular diameter distance
        Da_zs = cubic_spline_interpolator(np.array([zs]), da_args_function_spline, da_args_x_array)[0] # float
        Da_zl = cubic_spline_interpolator(np.array([zl]), da_args_function_spline, da_args_x_array)[0]

        # einstein radius 
        Dls = (Da_zs*(1+zs) - Da_zl*(1+zl))/(1+zs)
        theta_E = (
            4.0 * np.pi * (sigma / 299792.458) ** 2 * Dls / Da_zs
        )  # Note: km/s for sigma; Dls, Ds are in Mpc

        # cross-section
        cross_section = phi_cut_SIE(q) * np.pi*theta_E**2

        # rejection sampling
        max_ = np.max(cross_section)  # maximum einstein radius
        if max_ > max_value:
            max_value = max_
        u = np.random.uniform(0, max_value, size=buffer_size)
        mask = u < cross_section
        if np.sum(mask) > 0:
            break

    sigma_final = sigma[mask][0]
    theta_E_final = theta_E[mask][0]
    q_final = q[mask][0]
    # print(idx)

    return(idx, sigma_final, theta_E_final, q_final)