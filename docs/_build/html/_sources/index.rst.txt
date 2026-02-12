Welcome to :red:`ler`'s documentation!
=============================================

.. raw:: html

   <div style="text-align: center; margin-top: 0px; margin-bottom: 10px; padding-top: 0px; padding-bottom: 10px;">
      <img src="_static/lerlogo.png" width="220rem" alt="ler logo">
   </div>

:red:`ler` : :red_first:`LVK` (LIGO-Virgo-KAGRA) :red_first:`Event` (compact-binary mergers) :red_first:`Rate` calculator and simulator
----------------------------------------------------------------------------------------

``ler`` (/ˈɛlɚ/) is a statistical-based Python package whose core function is to calculate detectable rates of gravitational wave events (both lensed and unlensed). It provides a comprehensive framework for simulating compact binary coalescence events and their detection by ground-based gravitational wave detectors.

``ler`` is closely integrated with the ``gwsnr`` package (`see gwsnr documentation <https://gwsnr.hemantaph.com>`_), which provides efficient gravitational-wave Signal-to-Noise Ratio calculation. This allows researchers to perform fast $P_{\rm det}$ calculations for population studies and hierarchical Bayesian inference with selection effects.

For a detailed technical overview, see the :doc:`Summary` section, and browse the topics under :blue:`CONTENTS` in the sidebar (or top-left menu) of your browser.


Quick Start
-----------

Install the package using pip:

.. code-block:: bash

   pip install ler

Then, calculate gravitational wave event rates:

.. code-block:: python

   from ler import LeR

   # Initialize the rate calculator
   ler = LeR()

   # Generate unlensed and lensed event parameters
   unlensed_param = ler.unlensed_cbc_statistics(size=1000)
   lensed_param = ler.lensed_cbc_statistics(size=1000)

   # Calculate detectable rates
   unlensed_rate, unlensed_param_detected = ler.unlensed_rate()
   lensed_rate, lensed_param_detected = ler.lensed_rate()

   # compare rates
   ratio = ler.rate_ratio()

.. note::

   ``ler`` supports Python 3.10+ (but 3.11 recommended) and utilizes multi-core CPUs and multi-threaded CPUs when available. Refer to the :doc:`Installation` section for detailed setup instructions.

About the Project
-----------

* **Source Code:** `github.com/hemantaph/ler <https://github.com/hemantaph/ler>`_
* **Issue Tracker:** `Report an issue <https://github.com/hemantaph/ler/issues>`_
* **Main Developer:** `Hemanta Ph. <https://www.hemantaph.com>`_
* **Contributors:** Anupreeta More, Harsh Narola, Ng Chung Yin (Leo), Justin Janquart, Chris Van Den Broeck, Otto Akseli Hannuksela, Neha Singh, David Keitel
* **Citation:** If you use ``ler`` in your research, please cite the `ler paper <https://arxiv.org/abs/2306.03827>`_.



.. _glossary:

Glossary
========

.. glossary::

   Gravitational waves

      Ripples in the fabric of spacetime, first predicted by Albert Einstein in his theory of General Relativity in 1916. They are created by some of the most violent and energetic events in the universe, such as the collision of black holes, the merging of neutron stars, or supernova explosions. These waves travel outward from their source at the speed of light, carrying information about their origins and the nature of gravity itself.

      For a century, gravitational waves remained a theoretical prediction. It wasn't until 2015 that the LIGO and Virgo collaborations made the first direct detection, an achievement that earned the 2017 Nobel Prize in Physics and opened an entirely new way of observing the cosmos.

      .. raw:: html

         <div style="text-align:center;">
         <img src="_static/gw.gif" width="480px" alt="Animation of GW propagation">
          <div style="text-align:left; max-width:480px; margin: 0 auto;">
            <p style="font-size: 0.9em; font-family: Arial, sans-serif; line-height: 1.5em;">
               Animation showing the propagation of gravitational waves from inspiraling binary black holes. As the waves travel, they stretch and squeeze spacetime in their path. <em>Source: <a href="https://community.wolfram.com/groups/-/m/t/790989">Jeffrey Bryant, Wolfram | Alpha, LLC.</a>.</em>
            </p>
         </div>
         </div>

   Lensing of gravitational waves

      A process similar to the gravitational lensing of light, where gravitational waves emitted from distant astrophysical events are bent and split (strong-lensing case) into multiple images by the gravity of intervening massive objects, such as galaxies and galaxy clusters. This can magnify and change the arrival time of the gravitational waves.

      .. raw:: html

         <div style="text-align: center;">
         <iframe src="_static/gwlensing.html"
            width="90%"
            height="500"
            frameborder="0"></iframe>
         <div style="text-align:left; max-width:600px; margin: 0 auto;">
            <p style="font-size: 0.9em; font-family: Arial, sans-serif; line-height: 1.5em;">
               Interactive animation showing the lensing of gravitational waves by a massive object.
            </p>
         </div>
         </div>

      .. raw:: html

         <div style="text-align:center;">
         <img src="_static/Light_Bending.gif" width="480px" alt="Animation of light bending">
          <div style="text-align:left; max-width:480px; margin: 0 auto;">
            <p style="font-size: 0.9em; font-family: Arial, sans-serif; line-height: 1.5em;">
               Animation showing the bending of light by a massive object. <em>Source: <a href="https://science.nasa.gov/universe/how-gravity-warps-light/">NASA, ESA, and Goddard Space Flight Center/K. Jackson.</a>.</em>
            </p>
         </div>
         </div>

   Rate calculation and statistics

      In the context of astrophysics and gravitational wave science, rate calculation involves predicting the frequency of events such as compact-binary mergers based on observed data and theoretical models. Statistical methods are used to analyze and interpret the data, estimate parameters, and test hypotheses.

      .. raw:: html

         <div style="text-align:center;">
         <img src="_static/BBH_Merger_rate_density_and_PDF_of_redshift.png" width="480px" alt="BBH Merger rate density">
          <div style="text-align:left; max-width:480px; margin: 0 auto;">
            <p style="font-size: 0.9em; font-family: Arial, sans-serif; line-height: 1.5em;">
               Figure showing the merger rate density of binary black hole (BBH) mergers as a function of redshift. <em>Source: <a href="https://ler.hemantaph.comGW_events.html">generated by the LeR package.</a>.</em>
            </p>
         </div>
         </div>

      The ``ler`` package efficiently computes detectable event rates for large-scale astrophysical simulations, supporting both lensed and unlensed gravitational wave events.


.. toctree::
   :maxdepth: 2
   :caption: Contents:

   Installation
   Summary
   code_overview
   analytical_formulation_unlensed
   analytical_formulation_lensed
   
.. toctree::
   :maxdepth: 2
   :caption: API:

   autoapi/ler/rates/index
   autoapi/ler/gw_source_population/index
   autoapi/ler/lens_galaxy_population/index
   autoapi/ler/image_properties/index
   autoapi/ler/utils/index
   
.. toctree::
   :maxdepth: 2
   :caption: Examples:
   
   examples/LeR_short_examples
   examples/GWRATES_complete_examples
   examples/LeR_custom_functions
   examples/LeR_advanced_sampling


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`