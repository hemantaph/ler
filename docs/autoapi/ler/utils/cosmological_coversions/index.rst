:py:mod:`ler.utils.cosmological_coversions`
===========================================

.. py:module:: ler.utils.cosmological_coversions


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.utils.cosmological_coversions.redshift_optimal_spacing
   ler.utils.cosmological_coversions.luminosity_distance
   ler.utils.cosmological_coversions.differential_comoving_volume
   ler.utils.cosmological_coversions.comoving_distance
   ler.utils.cosmological_coversions.angular_diameter_distance
   ler.utils.cosmological_coversions.angular_diameter_distance_z1z2



.. py:function:: redshift_optimal_spacing(z_min, z_max, resolution)


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


