:py:mod:`ler.helperroutines`
============================

.. py:module:: ler.helperroutines

.. autoapi-nested-parse::

   Helper functions for various tasks used in LER, for example combining dictionaries together. Really this is a place for routines which don't seem to fit into anywhere else.

   ..
       !! processed by numpydoc !!


Module Contents
---------------


Functions
~~~~~~~~~

.. autoapisummary::

   ler.helperroutines.trim_dictionary
   ler.helperroutines.trim_dictionary_by_indices
   ler.helperroutines.add_dictionaries_together
   ler.helperroutines.rejection_sample
   ler.helperroutines.combine_lens_parameter_dictionaries
   ler.helperroutines.save_dictionary_to_numpy_txt_file



.. py:function:: trim_dictionary(dictionary, size)

   
   Filters an event dictionary to only contain the size.
















   ..
       !! processed by numpydoc !!

.. py:function:: trim_dictionary_by_indices(dictionary, indices)

   
   Filters an event dictionary to only contain the indices.
















   ..
       !! processed by numpydoc !!

.. py:function:: add_dictionaries_together(dictionary1, dictionary2)

   
   Adds two dictionaries with the same keys together.
















   ..
       !! processed by numpydoc !!

.. py:function:: rejection_sample(pdf, xmin, xmax, size=100)

   
   Helper function for rejection sampling from a pdf with maximum and minimum arguments.
   Input parameters:
       pdf: the pdf to sample from
       xmin: the minimum argument of the pdf
       xmax: the maximum argument of the pdf
       size: the number of samples to draw
   Output:
       samples: the samples drawn from the pdf
















   ..
       !! processed by numpydoc !!

.. py:function:: combine_lens_parameter_dictionaries(lensed_parameters, lensed_parameters_draw, idx, n_images)

   
   Adds lensed_parameters_draw to lensed_parameters dictionary for selected events idx and for n_images.

   Input parameters:
   lensed_parameters (dict): Dictionary of lensed parameters
   lensed_parameters_draw (dict): Dictionary of lensed parameters to be added to lensed_parameters
   idx (int): Index of the events to be added to lensed_parameters
   n_images (int): Number of images

   Output parameters:
   lensed_parameters (dict): Dictionary of lensed parameters















   ..
       !! processed by numpydoc !!

.. py:function:: save_dictionary_to_numpy_txt_file(detectable_lensed_event_parameters, fname='detectable_lensed_event_parameters.txt')

   
   Saves a dictionary to a numpy txt file.

   :param detectable_lensed_event_parameters: Dictionary to be saved
   :type detectable_lensed_event_parameters: dict
   :param fname: Name of the file to be saved
   :type fname: str

   Example:
   from ler import helperroutines as hr
   # Save the detectable lensed event parameters dictionary as numpy txt file
   hr.save_dictionary_to_numpy_txt_file(detectable_lensed_event_parameters, fname= 'detectable_lensed_event_parameters.txt' )
   # Load the detectable lensed event parameters dictionary
   data = np.genfromtxt('detectable_lensed_4_image_event_parameters.txt', names=True)
   names = data.dtype.names















   ..
       !! processed by numpydoc !!

