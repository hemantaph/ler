:py:mod:`ler.gw_source_population.sfr_with_time_delay`
======================================================

.. py:module:: ler.gw_source_population.sfr_with_time_delay

.. autoapi-nested-parse::

   Module for computing star formation rates with time delays.

   This module provides functions for computing the star formation rate at a given
   redshift, accounting for time delays between formation and observation. The time
   delay distribution follows a 1/t power-law form, and the formation redshift is
   computed using the cosmological age-redshift relation.

   Key Features:

   - Time-delayed star formation rate computation

   - Integration with Madau & Fragos (2017) SFR model

   - Monte Carlo integration for time delay averaging

   - Cosmological calculations using Astropy


   Copyright (C) 2026 Hemanta Ph. Distributed under MIT License.

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.gw_source_population.sfr_with_time_delay.sfr_with_time_delay_function



.. py:function:: sfr_with_time_delay_function(input_args)

   
   Compute star formation rate at observed redshift with time delay.

   The star formation rate is time-delayed relative to the observed redshift,
   with a time delay uniformly distributed between td_min and td_max. The
   formation redshift is computed using the cosmological age-redshift relation.

   :Parameters:

       **input_args** : ``list``
           List containing the following elements in order:

           - z (``float``): Observed redshift

           - idx (``int``): Index identifier for the computation

           - td_min (``float``): Minimum time delay (Gyr)

           - td_max (``float``): Maximum time delay (Gyr)

           - H0 (``float``): Hubble constant (km/s/Mpc)

           - Omega_M (``float``): Matter density parameter

           - Omega_Lambda (``float``): Dark energy density parameter

           - a (``float``): Madau-Fragos SFR normalization parameter

           - b (``float``): Madau-Fragos low-z power-law slope

           - c (``float``): Madau-Fragos turnover parameter

           - d (``float``): Madau-Fragos high-z power-law slope

   :Returns:

       **idx** : ``int``
           Index identifier (same as input).

       **result** : ``float``
           Time-averaged star formation rate at observed redshift z.










   .. rubric:: Examples

   >>> from ler.gw_source_population.sfr_with_time_delay import sfr_with_time_delay
   >>> args = [0.5, 0, 0.02, 13.0, 70.0, 0.3, 0.7, 0.01, 2.6, 3.2, 6.2]
   >>> idx, sfr = sfr_with_time_delay(args)



   ..
       !! processed by numpydoc !!

