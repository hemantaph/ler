# -*- coding: utf-8 -*-
"""
Module for computing star formation rates with time delays.

This module provides functions for computing the star formation rate at a given
redshift, accounting for time delays between formation and observation. The time
delay distribution follows a 1/t power-law form, and the formation redshift is
computed using the cosmological age-redshift relation.

Key Features: \n
- Time-delayed star formation rate computation \n
- Integration with Madau & Fragos (2017) SFR model \n
- Monte Carlo integration for time delay averaging \n
- Cosmological calculations using Astropy \n

Copyright (C) 2026 Hemanta Ph. Distributed under MIT License.
"""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import fsolve
from astropy.cosmology import LambdaCDM

from .prior_functions  import sfr_madau_fragos2017


def sfr_with_time_delay_function(input_args):
    """
    Compute star formation rate at observed redshift with time delay.

    The star formation rate is time-delayed relative to the observed redshift,
    with a time delay uniformly distributed between td_min and td_max. The
    formation redshift is computed using the cosmological age-redshift relation.

    Parameters
    ----------
    input_args : ``list``
        List containing the following elements in order: \n
        - z (``float``): Observed redshift \n
        - idx (``int``): Index identifier for the computation \n
        - td_min (``float``): Minimum time delay (Gyr) \n
        - td_max (``float``): Maximum time delay (Gyr) \n
        - H0 (``float``): Hubble constant (km/s/Mpc) \n
        - Omega_M (``float``): Matter density parameter \n
        - Omega_Lambda (``float``): Dark energy density parameter \n
        - a (``float``): Madau-Fragos SFR normalization parameter \n
        - b (``float``): Madau-Fragos low-z power-law slope \n
        - c (``float``): Madau-Fragos turnover parameter \n
        - d (``float``): Madau-Fragos high-z power-law slope \n

    Returns
    -------
    idx : ``int``
        Index identifier (same as input).
    result : ``float``
        Time-averaged star formation rate at observed redshift z.

    Examples
    --------
    >>> from ler.gw_source_population.sfr_with_time_delay import sfr_with_time_delay
    >>> args = [0.5, 0, 0.02, 13.0, 70.0, 0.3, 0.7, 0.01, 2.6, 3.2, 6.2]
    >>> idx, sfr = sfr_with_time_delay(args)
    """
    z = input_args[0]
    idx = input_args[1]
    td_min = input_args[2]
    td_max = input_args[3]
    H0 = input_args[4]
    Omega_M = input_args[5]
    Omega_Lambda = input_args[6]
    a = input_args[7]
    b = input_args[8]
    c = input_args[9]
    d = input_args[10]

    def _E(z_prime):
        """Compute dimensionless Hubble parameter E(z)."""
        return np.sqrt(Omega_M * (1 + z_prime)**3 + Omega_Lambda)

    def _integrand(z_prime):
        """Compute integrand for lookback time calculation."""
        return 1 / (H0 * (1 + z_prime) * _E(z_prime)) * 977.813

    def _time_delay(zform, z):
        """Compute time delay between formation and observation redshifts."""
        integral, _ = quad(_integrand, z, zform)
        return integral

    def _equation_to_solve(zform, z, td):
        """Residual equation for finding formation redshift."""
        return td - _time_delay(zform, z)

    def _find_zform(z, td):
        """Find formation redshift for given observation redshift and time delay."""
        zform_solution = fsolve(_equation_to_solve, z, args=(z, td))
        return zform_solution

    def _integrand_rates(z, size=100000, zform_max=1000.):
        """
        Compute time-averaged SFR using Monte Carlo integration.

        Parameters
        ----------
        z : ``float``
            Observed redshift.
        size : ``int``
            Number of Monte Carlo samples. \n
            default: 100000
        zform_max : ``float``
            Maximum formation redshift cutoff. \n
            default: 1000.0

        Returns
        -------
        integral : ``float``
            Time-averaged star formation rate.
        """
        td = np.random.uniform(td_min, td_max, size)
        td_max_allowed = _time_delay(zform_max, z)
        idx_valid = np.where(td < td_max_allowed)[0]
        P_td = np.zeros_like(td)
        P_td[idx_valid] = 1 / (np.log(td_max / td_min) * td[idx_valid])

        zform = np.zeros_like(td)
        for idx_ in idx_valid:
            zform[idx_] = _find_zform(z, td[idx_])

        psi = sfr_madau_fragos2017(zform, a, b, c, d)

        integral = 1 / (td_max - td_min) * np.sum(P_td * psi)
        return integral

    result = _integrand_rates(z)
    return int(idx), result