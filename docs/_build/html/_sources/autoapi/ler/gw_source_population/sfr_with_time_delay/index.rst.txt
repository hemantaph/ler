:py:mod:`ler.gw_source_population.sfr_with_time_delay`
======================================================

.. py:module:: ler.gw_source_population.sfr_with_time_delay


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.gw_source_population.sfr_with_time_delay.sfr_with_time_delay



Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.gw_source_population.sfr_with_time_delay.cosmo


.. py:data:: cosmo

   

.. py:function:: sfr_with_time_delay(input_args)

   
   Compute the star formation rate at redshift z, given parameters a, b, c, and d,
   and cosmological parameters H0, Omega_M, and Omega_Lambda.
   The star formation rate is time-delayed relative to the observed redshift,
   with a time delay uniformly distributed between td_min and td_max.
   The time delay is computed using the cosmology provided by astropy.


   :Parameters:

       **input_args** : list
           z : float
               observed redshift
           idx : int
               index of the galaxy
           td_min : float
               minimum time delay in Gyr
           td_max : float
               maximum time delay in Gyr
           H0 : float
               Hubble constant in km/s/Mpc
           Omega_M : float
               matter density parameter
           Omega_Lambda : float
               dark energy density parameter
           a : float
               parameter of the Madau-Fragos star formation rate
           b : float
               parameter of the Madau-Fragos star formation rate
           c : float
               parameter of the Madau-Fragos star formation rate
           d : float
               parameter of the Madau-Fragos star formation rate

   :Returns:

       **idx** : int
           index of the galaxy

       **result** : float
           star formation rate at observed redshift z













   ..
       !! processed by numpydoc !!

