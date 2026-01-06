:py:mod:`ler.lens_galaxy_population.lens_parameter_sampler`
===========================================================

.. py:module:: ler.lens_galaxy_population.lens_parameter_sampler


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.lens_parameter_sampler.rejection_sampler
   ler.lens_galaxy_population.lens_parameter_sampler.create_rejection_sampler
   ler.lens_galaxy_population.lens_parameter_sampler.sigma_proposal_uniform
   ler.lens_galaxy_population.lens_parameter_sampler.weighted_choice_1d
   ler.lens_galaxy_population.lens_parameter_sampler.importance_sampler
   ler.lens_galaxy_population.lens_parameter_sampler.create_importance_sampler



Attributes
~~~~~~~~~~

.. autoapisummary::

   ler.lens_galaxy_population.lens_parameter_sampler.C_LIGHT


.. py:data:: C_LIGHT
   :value: '299792.458'

   

.. py:function:: rejection_sampler(zs, zl, sigma_max, sigma_rvs, q_rvs, phi_rvs, gamma_rvs, shear_rvs, cross_section, saftey_factor=1.2)


.. py:function:: create_rejection_sampler(sigma_max, sigma_rvs, q_rvs, phi_rvs, gamma_rvs, shear_rvs, cross_section, saftey_factor=1.2, use_njit_sampler=True)


.. py:function:: sigma_proposal_uniform(n, sigma_min, sigma_max)


.. py:function:: weighted_choice_1d(weights)

   
   Draw an index with probability proportional to 'weights' (assumed >=0).
   Numba-safe replacement for np.random.choice(n, p=weights).
















   ..
       !! processed by numpydoc !!

.. py:function:: importance_sampler(zs, zl, sigma_min, sigma_max, q_rvs, phi_rvs, gamma_rvs, shear_rvs, sigma_pdf, cross_section, n_prop)


.. py:function:: create_importance_sampler(sigma_min, sigma_max, q_rvs, phi_rvs, gamma_rvs, shear_rvs, sigma_pdf, cross_section, n_prop, use_njit_sampler=True)


