Welcome to ler's documentation!
===============================

.. image:: ../lerlogo.png
   :align: center
   :width: 30%
   :alt: ler logo

==================
ler
==================

``ler`` (/ˈɛlɚ/): Lensed (or un-lensed) gravitational waves Event (compact-binaries) Rate calculator and simulator.

``ler`` is a statistical based python package whose core function is to calculate detectable rates of gravitational waves events. Description available in the `Summary <https://ler.readthedocs.io/en/latest/Summary.html>`_ section.

| The code is available at `github.ler <https://github.com/hemantaph/ler>`_.
| For reaching out to the developer, please raise issue in `github.ler.issue <https://github.com/hemantaph/ler/issues>`_.

| ``ler`` main developer: `Hemanta Ph. <https://newww.phy.cuhk.edu.hk/postgraduate/phurailatpam-hemantakumar>`_
| ``ler`` developer and analyst: Harsh Narola, Anupreeta More, Ng Chung Yin (Leo), Justin Janquart, Chris Van Den Broeck, Otto Akseli HANNUKSELA

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Installation2
   Summary
   code_overview
   GW_events
   GW_equations
   Lensed_events
   Lensed_equations
   
.. toctree::
   :maxdepth: 1
   :caption: API:

   autoapi/ler/rates/index.rst
   autoapi/ler/gw_source_population/index.rst
   autoapi/ler/lens_galaxy_population/index.rst
   autoapi/ler/image_properties/index.rst
   autoapi/ler/utils/index.rst
   
.. toctree::
   :maxdepth: 2
   :caption: Examples:
   
   examples/rates/LeR complete examples
   examples/rates/GWRATES complete exmaples
   examples/rates/grb detection rate
   examples/rates/ler bns example
   examples/rates/rates_with_3G_detectors
   examples/source population/merger rate density evolution with redshift
   examples/source population/sample compact binary parameters
   examples/source population/sample compact binary parameters with redshift
   examples/source population/statistical study of changing mass model params
   examples/source population/statistical study of cosmological friction
   examples/lens_parameters/sample lens parameters
   examples/optical_depth/validation_SIE
   examples/optical_depth/validation_SIS
   examples/image_properties/dt_vs_dmu


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`