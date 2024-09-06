============
Installation
============

.. note::

    For arm64 architecture processor (e.g. apple silicon) you might need to install some specific dependencies with conda (:ref:`dependencies`). For package development and contribution refer here (:ref:`development`). There can be error related to multiprocessing on either arm64 or x86_64 architecture. This can be resolved by using `main()` function as shown in this `example <https://github.com/hemantaph/ler/blob/main/examples/rates/LeR_short_example.py>`_. 

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