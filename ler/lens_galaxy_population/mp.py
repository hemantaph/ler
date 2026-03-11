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
from numba import njit, prange, set_num_threads

from .lens_functions import cross_section

# Global variable to hold shared data in worker processes.
# This avoids pickling large data for each work item — only once per worker.
_worker_shared_data = {}

# Constants for cross-section unit scaling (1/pi and numerical offset)
CS_UNIT_SLOPE = 0.31830988618379075
CS_UNIT_INTERCEPT = -3.2311742677852644e-27


@njit(parallel=True, cache=True)
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
    JIT-compiled parallel computation of differential optical dept (lens redshift).

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
            area_array = cross_section(zs_, zl_, sigma, q, phi, gamma, gamma1, gamma2)

            # Filter out invalid values (inf, non-positive)
            idx = np.logical_not(np.isinf(area_array))
            idx &= area_array > 0
            if idx.sum() == 0:
                result_array[i, j] = 0.0
                continue

            # Compute number density weights
            phi_sigma = number_density(sigma, zl_)

            # Differential comoving volume
            dVcdz = dVcdz_function(np.array([zl]))[0]

            # Compute optical depth contribution
            result = (
                (sigma_max - sigma_min)
                * np.average(area_array[idx] * phi_sigma[idx] * dVcdz)
                / (4 * np.pi)
            )
            result_array[i, j] = result

    return result_array


def _init_lens_redshift_worker(
    sigma_min,
    sigma_max,
    sigma_function,
    q_rvs,
    phi_rvs,
    gamma_rvs,
    shear_rvs,
    dVcdz_function,
    cs_function,
    integration_size,
):
    """
    Initialize worker process with shared data for lens redshift.

    Called once per worker process when the pool is created.
    Shared data is stored in the module-level ``_worker_shared_data``
    global so that ``lens_redshift_strongly_lensed_mp`` can access it
    without any per-task pickling overhead.

    Parameters
    ----------
    sigma_min : ``float``
        Minimum velocity dispersion (km/s).
    sigma_max : ``float``
        Maximum velocity dispersion (km/s).
    sigma_function : ``callable``
        Velocity dispersion number density function.
    q_rvs : ``callable``
        Axis ratio sampler.
    phi_rvs : ``callable``
        Axis rotation angle sampler.
    gamma_rvs : ``callable``
        Density slope sampler.
    shear_rvs : ``callable``
        External shear sampler.
    dVcdz_function : ``callable``
        Differential comoving volume function.
    cs_function : ``callable``
        Cross-section function.
    integration_size : ``int``
        Number of Monte Carlo samples per (zs, zl) pair.
    """
    global _worker_shared_data
    _worker_shared_data["sigma_min"] = sigma_min
    _worker_shared_data["sigma_max"] = sigma_max
    _worker_shared_data["sigma_function"] = sigma_function
    _worker_shared_data["q_rvs"] = q_rvs
    _worker_shared_data["phi_rvs"] = phi_rvs
    _worker_shared_data["gamma_rvs"] = gamma_rvs
    _worker_shared_data["shear_rvs"] = shear_rvs
    _worker_shared_data["dVcdz_function"] = dVcdz_function
    _worker_shared_data["cs_function"] = cs_function
    _worker_shared_data["integration_size"] = integration_size


def lens_redshift_strongly_lensed_mp(params):
    """
    Multiprocessing worker for computation of differential optical dept (lens redshift).

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

    # Re-seed RNG from OS entropy so forked workers don't share the parent's state
    np.random.seed()

    # Limit to 1 Numba thread per worker — parallelism comes from the process pool
    set_num_threads(1)

    # Retrieve shared data set by _init_lens_redshift_worker
    sigma_min = _worker_shared_data["sigma_min"]
    sigma_max = _worker_shared_data["sigma_max"]
    sigma_function = _worker_shared_data["sigma_function"]
    q_rvs = _worker_shared_data["q_rvs"]
    phi_rvs = _worker_shared_data["phi_rvs"]
    gamma_rvs = _worker_shared_data["gamma_rvs"]
    shear_rvs = _worker_shared_data["shear_rvs"]
    dVcdz_function = _worker_shared_data["dVcdz_function"]
    cs_function = _worker_shared_data["cs_function"]
    integration_size = _worker_shared_data["integration_size"]

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
        valid_idx &= area_array > 0
        if valid_idx.sum() == 0:
            result_array[i] = 0.0
            continue

        # Compute number density weights (velocity dispersion distribution)
        phi_sigma = sigma_function(sigma, zl_arr)

        # Differential comoving volume
        dVcdz = dVcdz_function(np.array([zl]))[0]

        # Compute optical depth contribution (importance sampling correction)
        result = (
            (sigma_max - sigma_min)
            * np.average(area_array[valid_idx] * phi_sigma[valid_idx] * dVcdz)
            / (4 * np.pi)
        )
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
        Cross-section area (dimensionless, for theta_E = 1).
    """
    set_num_threads(1)

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
        Packed parameters (theta_E, e1, e2, gamma, gamma1, gamma2, idx). theta_E is in radians.

    Returns
    -------
    idx : ``int``
        Worker index for result ordering.
    area : ``float``
        Cross-section area (radian^2).
    """
    set_num_threads(1)

    theta_E, e1, e2, gamma, gamma1, gamma2, idx = params
    area = cross_section(theta_E, e1, e2, gamma, gamma1, gamma2)
    return int(idx), area
