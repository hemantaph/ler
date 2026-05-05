:py:mod:`ler.utils.cosmological_coversions`
===========================================

.. py:module:: ler.utils.cosmological_coversions


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.utils.cosmological_coversions.generate_mixed_grid
   ler.utils.cosmological_coversions.luminosity_distance
   ler.utils.cosmological_coversions.differential_comoving_volume
   ler.utils.cosmological_coversions.comoving_distance
   ler.utils.cosmological_coversions.angular_diameter_distance
   ler.utils.cosmological_coversions.angular_diameter_distance_z1z2



.. py:function:: generate_mixed_grid(x_min, x_max, resolution, power_law_part='lower', geomspace_part=False, spacing_trend='increasing', power=2.3, value_transition_fraction=0.6, num_transition_fraction=0.8, auto_match_slope=True)

   
   Generalized mixed spacing grid generator. Safely handles negative ranges.


   :Parameters:

       **x_min** : float
           Minimum value of the grid.

       **x_max** : float
           Maximum value of the grid.

       **resolution** : int
           Total number of grid points.

       **power_law_part** : str, optional
           Which part of the grid should follow the power-law spacing. Options: 'lower' or 'upper'. Default is 'lower'.

       **geomspace_part** : bool or str, optional
           If `False`, keep the existing linear + power-law behavior. If `'lower'` or `'upper'`,
           replace that segment with geometric spacing while keeping the other segment linear.
           Geometric spacing is only used when the selected segment endpoints are strictly positive;
           otherwise the function falls back to the standard mixed-grid construction. Default is `False`.

       **spacing_trend** : str, optional
           Whether the power-law spacing should be increasing or decreasing. Options: 'increasing' or 'decreasing'. Default is 'increasing'.

       **power** : float, optional
           The power-law exponent. Higher values lead to more extreme spacing. Default is 2.3.

       **value_transition_fraction** : float, optional
           The fraction of the total value range at which to transition from linear to power-law spacing. Must be between 0 and 1. Default is 0.6.

       **num_transition_fraction** : float, optional
           The fraction of the total number of points at which to transition from linear to power-law spacing. Must be between 0 and 1. Default is 0.8.

       **auto_match_slope** : bool, optional
           Whether to automatically adjust the power-law exponent to match the slope of the linear spacing at the transition point. Default is True.
           This is ignored for the geometric-spacing segment when `geomspace_part` is used.

   :Returns:

       numpy.ndarray
           The generated grid points.










   .. rubric:: Examples

   from ler.utils.cosmological_conversions import generate_mixed_grid

   resolution=20
   # linear+power-law with power-law in the upper segment and decreasing step sizes
   x = generate_mixed_grid(
       x_min=0.0, x_max=10.0, resolution=resolution,
       power_law_part='upper',
       spacing_trend='decreasing',  # Forces largest steps near z_trans
       power=2.5,
       value_transition_fraction=0.6,
       num_transition_fraction=0.3,
       auto_match_slope=True       # We accept the kink to control the exact power
   )
   # powerlaw+linear with power-law in the lower segment and increasing step sizes
   x = generate_mixed_grid(
       x_min=0.0, x_max=10.0, resolution=resolution,
       power_law_part='lower',
       spacing_trend='increasing',  # Forces largest steps near z_trans
       power=2.5,
       value_transition_fraction=0.3,
       num_transition_fraction=0.6,
       auto_match_slope=True       # We accept the kink to control the exact power
   )
   # linear+geomspace with geometric spacing in the upper segment
   x = generate_mixed_grid(
       x_min=0.1, x_max=10.0, resolution=resolution,
       geomspace_part='lower',
       value_transition_fraction=0.3,
       num_transition_fraction=0.6,
   )



   ..
       !! processed by numpydoc !!

.. py:function:: luminosity_distance(z=None, z_min=0.001, z_max=10.0, cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0), directory='./interpolator_json', create_new=False, resolution=500, get_attribute=True)

   
   Function to create a lookup table for the luminosity distance wrt redshift.


   :Parameters:

       **z** : `numpy.ndarray` or `float`
           Source redshifts

       **z_min** : `float`
           Minimum redshift of the source population

       **z_max** : `float`
           Maximum redshift of the source population












   :Attributes:

       **z_to_luminosity_distance** : `ler.utils.FunctionConditioning`
           Object of FunctionConditioning class containing the luminosity distance wrt redshift


   ..
       !! processed by numpydoc !!

.. py:function:: differential_comoving_volume(z=None, z_min=0.001, z_max=10.0, cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0), directory='./interpolator_json', create_new=False, resolution=500, get_attribute=True)


.. py:function:: comoving_distance(z=None, z_min=0.001, z_max=10.0, cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0), directory='./interpolator_json', create_new=False, resolution=500, get_attribute=True)


.. py:function:: angular_diameter_distance(z=None, z_min=0.001, z_max=10.0, cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0), directory='./interpolator_json', create_new=False, resolution=500, get_attribute=True)


.. py:function:: angular_diameter_distance_z1z2(z1=None, z2=None, z_min=0.001, z_max=10.0, cosmo=LambdaCDM(H0=70, Om0=0.3, Ode0=0.7, Tcmb0=0.0, Neff=3.04, m_nu=None, Ob0=0.0), directory='./interpolator_json', create_new=False, resolution=500, get_attribute=True)


