============
Installation
============

.. note::

    For arm64 architecture processor (e.g. apple silicon) you might need to install some of the dependencies with conda (:ref:`dependencies`). For package development and contribution refer here (:ref:`development`).

.. tabs::
        
     .. code-tab:: console pip

        pip install ler


This will also install the dependencies needed by the lastest ``ler`` version.  

.. _development:
ler for development
======================

To install ``ler`` for development purposes use `github.ler <https://github.com/hemantaph/ler/>`_. Use conda environment to avoid dependency error. 

    
.. tabs::

     .. code-tab:: Console with new conda env

        git clone https://github.com/hemantaph/ler.git
        cd ler
        conda env create -f ler.yml
        conda activate ler
        pip install -e .
        
     .. code-tab:: with existing conda env
     
        git clone https://github.com/hemantaph/ler.git
        cd ler
        conda env update --file ler.yml
        pip install -e .
    
.. _dependencies:
Installation of numba with conda
=======================

.. code-block:: console

    conda install -c conda-forge numba