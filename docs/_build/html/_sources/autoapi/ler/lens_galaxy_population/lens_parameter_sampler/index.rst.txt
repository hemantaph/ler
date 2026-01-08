:py:mod:`ler.lens_galaxy_population.lens_parameter_sampler`
===========================================================

.. py:module:: ler.lens_galaxy_population.lens_parameter_sampler

.. autoapi-nested-parse::

   Module for lens parameter sampling using rejection and importance sampling methods.

   This module provides functions to sample lens galaxy parameters (velocity dispersion,
   axis ratio, orientation, power-law index, and external shear) weighted by the
   gravitational lensing cross section. Two sampling strategies are implemented:

   - Rejection sampling: Uses accept/reject method with cross section as the target density

   - Importance sampling: Uses uniform proposal with cross section-based importance weights


   Both samplers can be optionally JIT-compiled using Numba for improved performance.

   Copyright (C) 2024 LeR Team. Distributed under MIT License.

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.lens_parameter_sampler.create_rejection_sampler
   ler.lens_galaxy_population.lens_parameter_sampler.create_importance_sampler



.. py:function:: create_rejection_sampler(sigma_max, sigma_rvs, q_rvs, phi_rvs, gamma_rvs, shear_rvs, cross_section, safety_factor=1.2, use_njit_sampler=True)

   
   Create a rejection sampler for lens parameters weighted by cross section.

   Returns a callable that samples lens parameters using rejection sampling,
   optionally JIT-compiled for improved performance.

   :Parameters:

       **sigma_max** : ``float``
           Maximum velocity dispersion (km/s) for computing upper bound.

       **sigma_rvs** : ``callable``
           Function to sample velocity dispersion: sigma_rvs(n, zl) -> array.

       **q_rvs** : ``callable``
           Function to sample axis ratio: q_rvs(n, sigma) -> array.

       **phi_rvs** : ``callable``
           Function to sample orientation angle: phi_rvs(n) -> array.

       **gamma_rvs** : ``callable``
           Function to sample power-law index: gamma_rvs(n) -> array.

       **shear_rvs** : ``callable``
           Function to sample external shear: shear_rvs(n) -> (gamma1, gamma2).

       **cross_section** : ``callable``
           Function to compute cross section.

       **safety_factor** : ``float``
           Multiplicative safety factor for the upper bound.

           default: 1.2

       **use_njit_sampler** : ``bool``
           If True, uses Numba JIT compilation for faster execution.

           default: True

   :Returns:

       **rejection_sampler_wrapper** : ``callable``
           Function with signature (zs, zl) -> (sigma, q, phi, gamma, gamma1, gamma2).













   ..
       !! processed by numpydoc !!

.. py:function:: create_importance_sampler(sigma_min, sigma_max, q_rvs, phi_rvs, gamma_rvs, shear_rvs, sigma_pdf, cross_section, n_prop, use_njit_sampler=True)

   
   Create an importance sampler for lens parameters weighted by cross section.

   Returns a callable that samples lens parameters using importance sampling
   with uniform proposal, optionally JIT-compiled for improved performance.

   :Parameters:

       **sigma_min** : ``float``
           Minimum velocity dispersion (km/s) for uniform proposal.

       **sigma_max** : ``float``
           Maximum velocity dispersion (km/s) for uniform proposal.

       **q_rvs** : ``callable``
           Function to sample axis ratio: q_rvs(n, sigma) -> array.

       **phi_rvs** : ``callable``
           Function to sample orientation angle: phi_rvs(n) -> array.

       **gamma_rvs** : ``callable``
           Function to sample power-law index: gamma_rvs(n) -> array.

       **shear_rvs** : ``callable``
           Function to sample external shear: shear_rvs(n) -> (gamma1, gamma2).

       **sigma_pdf** : ``callable``
           PDF of velocity dispersion: sigma_pdf(sigma, zl) -> array.

       **cross_section** : ``callable``
           Function to compute cross section.

       **n_prop** : ``int``
           Number of proposal samples per lens.

       **use_njit_sampler** : ``bool``
           If True, uses Numba JIT compilation for faster execution.

           default: True

   :Returns:

       **importance_sampler_wrapper** : ``callable``
           Function with signature (zs, zl) -> (sigma, q, phi, gamma, gamma1, gamma2).













   ..
       !! processed by numpydoc !!

