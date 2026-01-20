# -*- coding: utf-8 -*-
"""
Multiprocessing helper functions for lens galaxy population calculations.

This module provides parallelized functions for computing optical depth 
and cross-sections for strong gravitational lensing. These functions are 
designed to be used with Python's multiprocessing module for efficient 
Monte Carlo integration over lens parameters.

Key Features: \n
- JIT-compiled parallel function for lens redshift sampling \n
- Multiprocessing-compatible wrapper functions for optical depth \n
- Cross-section calculation using EPL+Shear lens models \n

Copyright (C) 2026 Phurailatpam Hemantakumar. Distributed under MIT License.
"""

import numpy as np
from numba import njit, prange

from .lens_functions import cross_section

from ..utils import load_pickle

# Constants for cross-section unit scaling (1/pi and numerical offset)
CS_UNIT_SLOPE = 0.31830988618379075
CS_UNIT_INTERCEPT = -3.2311742677852644e-27


@njit(parallel=True)
def lens_redshift_strongly_lensed_njit(
        zs_array,
        zl_scaled,
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
    """
    JIT-compiled parallel computation of lens redshift optical depth.

    Computes the differential optical depth for strong lensing as a function 
    of source and lens redshifts using Monte Carlo integration with parallel 
    execution over source redshifts.

    Parameters
    ----------
    zs_array : ``numpy.ndarray``
        1D array of source redshifts.
    zl_scaled : ``numpy.ndarray``
        2D array of scaled lens redshifts (zl/zs).
    sigma_min : ``float``
        Minimum velocity dispersion (km/s).
    sigma_max : ``float``
        Maximum velocity dispersion (km/s).
    q_rvs : ``callable``
        Function to sample axis ratios given size and sigma.
    phi_rvs : ``callable``
        Function to sample axis rotation angles.
    gamma_rvs : ``callable``
        Function to sample density profile slopes.
    shear_rvs : ``callable``
        Function to sample external shear components (gamma1, gamma2).
    number_density : ``callable``
        Function to compute velocity dispersion number density.
    cross_section : ``callable``
        Function to compute lensing cross-section.
    dVcdz_function : ``callable``
        Function to compute differential comoving volume.
    integration_size : ``int``
        Number of Monte Carlo samples per (zs, zl) pair.

    Returns
    -------
    result_array : ``numpy.ndarray``
        2D array of optical depth values with shape (len(zs_array), len(zl_scaled[0])).
    """
    size_zs = zs_array.shape[0]
    size_zl = zl_scaled.shape[1]
    result_array = np.zeros((size_zs, size_zl))

    for i in prange(size_zs):
        # Unscale lens redshifts for this source redshift
        zl_array = zl_scaled[i] * zs_array[i]
        
        for j in range(size_zl):
            # Monte Carlo integration setup
            zs_ = zs_array[i] * np.ones(integration_size)
            zl = zl_array[j]
            zl_ = zl * np.ones(integration_size)

            # Sample lens parameters from their distributions
            sigma = np.random.uniform(sigma_min, sigma_max, integration_size)
            q = q_rvs(integration_size, sigma)
            phi = phi_rvs(integration_size)
            gamma = gamma_rvs(integration_size)
            gamma1, gamma2 = shear_rvs(integration_size)
            
            # Compute cross-sections
            area_array = cross_section(
                zs_, zl_, sigma, q, phi, gamma, gamma1, gamma2
            )

            # Filter out invalid values (inf, non-positive)
            idx = np.logical_not(np.isinf(area_array))
            idx &= (area_array > 0)
            if idx.sum() == 0:
                result_array[i, j] = 0.
                continue

            # Compute number density weights
            phi_sigma = number_density(sigma, zl_)
                
            # Differential comoving volume
            dVcdz = dVcdz_function(np.array([zl]))[0]

            # Compute optical depth contribution
            result = (sigma_max - sigma_min) * np.average(area_array[idx] * phi_sigma[idx] * dVcdz) / (4 * np.pi)
            result_array[i, j] = result
        
    return result_array


def lens_redshift_strongly_lensed_mp(params):
    """
    Multiprocessing worker for lens redshift optical depth calculation.

    Computes the differential optical depth for a single source redshift 
    across multiple lens redshifts. Designed to be called via multiprocessing 
    Pool.map() for parallel computation.

    Parameters
    ----------
    params : ``tuple``
        Packed parameters tuple containing: \n
        - params[0]: Source redshift (float) \n
        - params[1]: Scaled lens redshift array (1D array) \n
        - params[2]: Worker index (int) \n

    Returns
    -------
    idx : ``int``
        Worker index for result ordering.
    result_array : ``numpy.ndarray``
        1D array of optical depth values for each lens redshift.
    """

    # Load shared data from pickle file (set once per process)
    [
        sigma_args,
        q_rvs,
        dVcdz_function,
        cs_function,
        phi_rvs,
        shear_rvs,
        gamma_rvs,
        integration_size,
    ] = load_pickle('input_params_mp.pkl')
    sigma_min, sigma_max, sigma_function = sigma_args[0], sigma_args[1], sigma_args[2]

    zs = params[0]
    worker_idx = params[2]

    # Unscale lens redshifts
    zl_array = params[1] * zs

    result_array = np.zeros(len(zl_array))
    for i, zl in enumerate(zl_array):
        # Sample velocity dispersion (uniform for importance sampling)
        sigma = np.random.uniform(sigma_min, sigma_max, integration_size)

        # Sample axis ratio
        q = q_rvs(integration_size, sigma)

        # Sample axis rotation angle
        phi = phi_rvs(integration_size)

        # Sample density slope
        gamma = gamma_rvs(integration_size)

        # Sample external shear components (gamma1, gamma2)
        gamma1, gamma2 = shear_rvs(integration_size)

        # Compute cross-section
        zs_arr = zs * np.ones(integration_size)
        zl_arr = zl * np.ones(integration_size)
        area_array = cs_function(zs_arr, zl_arr, sigma, q, phi, gamma, gamma1, gamma2)

        # Filter out invalid values (inf, non-positive)
        valid_idx = np.logical_not(np.isinf(area_array))
        valid_idx &= (area_array > 0)
        if valid_idx.sum() == 0:
            result_array[i] = 0.
            continue

        # Compute number density weights (velocity dispersion distribution)
        phi_sigma = sigma_function(sigma, zl_arr)
            
        # Differential comoving volume 
        dVcdz = dVcdz_function(np.array([zl]))[0]

        # Compute optical depth contribution (importance sampling correction)
        result = (sigma_max - sigma_min) * np.average(area_array[valid_idx] * phi_sigma[valid_idx] * dVcdz) / (4 * np.pi)
        result_array[i] = result
    
    return worker_idx, np.array(result_array)


def cross_section_unit_mp(params):
    """
    Multiprocessing worker for unit Einstein radius cross-section.

    Computes the lensing cross-section for a lens with unit Einstein radius 
    (theta_E = 1). Used for building cross-section interpolation grids.

    Parameters
    ----------
    params : ``tuple``
        Packed parameters (e1, e2, gamma, gamma1, gamma2, idx).

    Returns
    -------
    idx : ``int``
        Worker index for result ordering.
    area : ``float``
        Cross-section area in square arcseconds.
    """
    e1, e2, gamma, gamma1, gamma2, idx = params
    area = cross_section(1.0, e1, e2, gamma, gamma1, gamma2)
    return idx, area


def cross_section_mp(params):
    """
    Multiprocessing worker for cross-section calculation.

    Computes the lensing cross-section for given lens parameters. 
    Designed to be called via multiprocessing Pool.map().

    Parameters
    ----------
    params : ``tuple``
        Packed parameters (theta_E, e1, e2, gamma, gamma1, gamma2, idx).

    Returns
    -------
    idx : ``int``
        Worker index for result ordering.
    area : ``float``
        Cross-section area in square arcseconds.
    """
    theta_E, e1, e2, gamma, gamma1, gamma2, idx = params
    area = cross_section(theta_E, e1, e2, gamma, gamma1, gamma2)
    return int(idx), area