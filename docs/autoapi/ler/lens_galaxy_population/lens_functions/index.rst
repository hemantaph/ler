:py:mod:`ler.lens_galaxy_population.lens_functions`
===================================================

.. py:module:: ler.lens_galaxy_population.lens_functions

.. autoapi-nested-parse::

   Lens functions for gravitational lensing calculations.

   This module provides utility functions for computing lens galaxy velocity
   dispersion functions, ellipticity conversions, and strong lensing cross-sections.
   These functions support the lens galaxy population modeling in the ``ler`` package.

   The velocity dispersion functions follow the models from Oguri et al. (2018b)
   and Bernardi et al. (2010). Cross-section calculations use lenstronomy for
   EPL+Shear lens models.

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.lens_functions.phi_cut_SIE
   ler.lens_galaxy_population.lens_functions.phi_q2_ellipticity
   ler.lens_galaxy_population.lens_functions.cross_section



.. py:function:: phi_cut_SIE(q)

   
   Calculate cross-section scaling factor for SIE lens galaxy from SIS.

   Computes the ratio of the SIE (Singular Isothermal Ellipsoid) cross-section
   to the SIS (Singular Isothermal Sphere) cross-section for a given axis ratio.

   :Parameters:

       **q** : ``numpy.ndarray``
           Axis ratio of the lens galaxy (0 < q <= 1).

   :Returns:

       **result** : ``numpy.ndarray``
           Scaling factor (normalized to pi).

           For q -> 1 (spherical): returns 1.0

           For q -> 0 (highly elliptical): returns ~0.













   ..
       !! processed by numpydoc !!

.. py:function:: phi_q2_ellipticity(phi, q)

   
   Convert position angle and axis ratio to ellipticity components.


   :Parameters:

       **phi** : ``numpy.ndarray``
           Position angle of the major axis (radians).

       **q** : ``numpy.ndarray``
           Axis ratio (0 < q <= 1).

   :Returns:

       **e1** : ``numpy.ndarray``
           First ellipticity component.

       **e2** : ``numpy.ndarray``
           Second ellipticity component.













   ..
       !! processed by numpydoc !!

.. py:function:: cross_section(theta_E, e1, e2, gamma, gamma1, gamma2)

   
   Compute the strong lensing cross-section for an EPL+Shear lens.

   Uses lenstronomy to compute the caustic structure and returns the
   area enclosed by the double-image (outer) caustic.

   :Parameters:

       **theta_E** : ``float``
           Einstein radius (arcseconds).

       **e1** : ``float``
           First ellipticity component.

       **e2** : ``float``
           Second ellipticity component.

       **gamma** : ``float``
           Power-law density profile slope.

       **gamma1** : ``float``
           First external shear component.

       **gamma2** : ``float``
           Second external shear component.

   :Returns:

       **area** : ``float``
           Cross-section area (arcseconds^2).

           Returns 0.0 if caustic computation fails.













   ..
       !! processed by numpydoc !!

