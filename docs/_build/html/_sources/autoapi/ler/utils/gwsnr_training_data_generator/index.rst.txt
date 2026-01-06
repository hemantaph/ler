:py:mod:`ler.utils.gwsnr_training_data_generator`
=================================================

.. py:module:: ler.utils.gwsnr_training_data_generator


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   ler.utils.gwsnr_training_data_generator.TrainingDataGenerator




.. py:class:: TrainingDataGenerator(npool=4, z_min=0.0, z_max=5.0, verbose=True, **kwargs)


   .. py:attribute:: npool
      :value: '4'

      

   .. py:attribute:: z_min
      :value: '0.0'

      

   .. py:attribute:: z_max
      :value: '5.0'

      

   .. py:attribute:: verbose
      :value: 'True'

      

   .. py:attribute:: ler_init_args

      

   .. py:method:: gw_parameters_generator(size, batch_size=100000, snr_recalculation=True, trim_to_size=False, verbose=True, replace=False, data_distribution_range=[0, 2, 4, 6, 8, 10, 12, 14, 16, 100], output_jsonfile='gw_parameters.json')


   .. py:method:: helper_data_distribution(gw_param, data_distribution_range)


   .. py:method:: combine_dicts(file_name_list=None, path_list=None, detector='L1', parameter_list=['mass_1', 'mass_2', 'luminosity_distance', 'theta_jn', 'psi', 'geocent_time', 'ra', 'dec', 'a_1', 'a_2', 'tilt_1', 'tilt_2'], output_jsonfile='combined_data.json')


   .. py:method:: delete_json_file(path_list)



