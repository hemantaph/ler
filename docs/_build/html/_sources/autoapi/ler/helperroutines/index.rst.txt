:py:mod:`ler.helperroutines`
============================

.. py:module:: ler.helperroutines

.. autoapi-nested-parse::

   This module contains helper routines for other modules in the ler package.

   ..
       !! processed by numpydoc !!


Module Contents
---------------

Classes
~~~~~~~

.. autoapisummary::

   ler.helperroutines.NumpyEncoder



Functions
~~~~~~~~~

.. autoapisummary::

   ler.helperroutines.append_json
   ler.helperroutines.get_param_from_json
   ler.helperroutines.rejection_sample
   ler.helperroutines.add_dictionaries_together
   ler.helperroutines.trim_dictionary



.. py:class:: NumpyEncoder(*, skipkeys=False, ensure_ascii=True, check_circular=True, allow_nan=True, sort_keys=False, indent=None, separators=None, default=None)

   Bases: :py:obj:`json.JSONEncoder`

   .. autoapi-inheritance-diagram:: ler.helperroutines.NumpyEncoder
      :parts: 1

   
   class for storing a numpy.ndarray or any nested-list composition as JSON file


   :Parameters:

       **json.JSONEncoder** : `class`
           class for encoding JSON file

   :Returns:

       **json.JSONEncoder.default** : `function`
           function for encoding JSON file













   ..
       !! processed by numpydoc !!
   .. py:method:: default(obj)

      
      Implement this method in a subclass such that it returns
      a serializable object for ``o``, or calls the base implementation
      (to raise a ``TypeError``).

      For example, to support arbitrary iterators, you could
      implement default like this::

          def default(self, o):
              try:
                  iterable = iter(o)
              except TypeError:
                  pass
              else:
                  return list(iterable)
              # Let the base class default method raise the TypeError
              return JSONEncoder.default(self, o)















      ..
          !! processed by numpydoc !!


.. py:function:: append_json(file_name, dictionary, replace=False)

   
   Append and update a json file with a dictionary.


   :Parameters:

       **file_name** : `str`
           json file name for storing the parameters.

       **dictionary** : `dict`
           dictionary to be appended to the json file.

       **replace** : `bool`, optional
           If True, replace the json file with the dictionary. Default is False.














   ..
       !! processed by numpydoc !!

.. py:function:: get_param_from_json(json_file)

   
   Function to get the parameters from json file.


   :Parameters:

       **json_file** : `str`
           json file name for storing the parameters.

   :Returns:

       **param** : `dict`
           ..













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

.. py:function:: add_dictionaries_together(dictionary1, dictionary2)

   
   Adds two dictionaries with the same keys together.
















   ..
       !! processed by numpydoc !!

.. py:function:: trim_dictionary(dictionary, size)

   
   Filters an event dictionary to only contain the size.
















   ..
       !! processed by numpydoc !!

