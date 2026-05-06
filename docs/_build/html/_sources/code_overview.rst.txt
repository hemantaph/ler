====================
Code overview
====================

The diagram below shows the main internal workflow of ``ler``: how the
public classes connect to the source-population samplers, lens-population
samplers, image solver, detectability calculation, rate estimators, and
JSON outputs.

.. image:: _static/ler_flowchart.svg
    :align: center
    :width: 100%
    :alt: Flowchart of the internal workflow of ler

Main entry points
-----------------

``GWRATES`` handles unlensed compact-binary population sampling and detection
rates. ``LeR`` extends that workflow to strongly lensed events by adding lens
galaxy sampling, optical-depth weighting, image-property calculations, and
image-level detectability.

The same source-population machinery is used in both paths. The lensed path
adds lens parameters and transforms each lensed image into effective GW
parameters before applying the same detection-probability and rate machinery.
