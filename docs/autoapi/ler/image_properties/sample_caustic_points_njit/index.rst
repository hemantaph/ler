:py:mod:`ler.image_properties.sample_caustic_points_njit`
=========================================================

.. py:module:: ler.image_properties.sample_caustic_points_njit

.. autoapi-nested-parse::

   Module for sampling source positions from EPL + shear caustic regions.

   Provides numba-accelerated routines for exact fan-triangulation sampling
   inside arbitrary star-convex polygons, with a convenience wrapper to
   sample from the double caustic of an EPL + external shear lens model.

   Usage:
       Basic workflow example:

       >>> from ler.image_properties.sample_caustic_points_njit import sample_source_from_double_caustic
       >>> xs, ys = sample_source_from_double_caustic(
       ...     theta_E=1.0, q=0.8, phi=0.0, gamma=2.0, gamma1=0.03, gamma2=-0.01
       ... )

   Copyright (C) 2026 Phurailatpam Hemantakumar. Distributed under MIT License.

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.image_properties.sample_caustic_points_njit.sample_source_from_double_caustic



.. py:function:: sample_source_from_double_caustic(theta_E, q, phi, gamma, gamma1, gamma2, num_th=500, maginf=-100.0)

   
   Sample a single source position from the double caustic region.

   Precomputes the caustic boundary once (the expensive part) and
   draws one uniform sample via exact fan-triangulation from the origin.

   :Parameters:

       **theta_E** : ``float``
           Einstein radius.

       **q** : ``float``
           Axis ratio.

       **phi** : ``float``
           Position angle of the lens major axis in radians.

       **gamma** : ``float``
           EPL slope exponent.

       **gamma1** : ``float``
           First shear component.

       **gamma2** : ``float``
           Second shear component.

       **num_th** : ``int``
           Number of angular samples for caustic boundary.

           default: 500

       **maginf** : ``float``
           Magnification cut threshold.

           default: -100.0

   :Returns:

       **xs** : ``float``
           Sampled source x-coordinate (``NaN`` if caustic is invalid).

       **ys** : ``float``
           Sampled source y-coordinate (``NaN`` if caustic is invalid).










   .. rubric:: Examples

   >>> xs, ys = sample_source_from_double_caustic(
   ...     theta_E=1.0, q=0.8, phi=0.0, gamma=2.0, gamma1=0.03, gamma2=-0.01
   ... )



   ..
       !! processed by numpydoc !!

