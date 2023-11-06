============
Installation
============

.. note::

    For arm64 architecture processor (e.g. apple silicon) you might need to install some of the dependencies with conda (:ref:`dependencies`). For package development and contribution refer here (:ref:`development`).

.. tabs::
        
     .. code-tab:: console pip

        pip install ler gwcosmo@git+https://git.ligo.org/lscsoft/gwcosmo.git@v1.0.0


This will also install the dependencies needed by the lastest ``ler`` version.  

.. _development:
ler for development
======================

To install ``ler`` for development purposes use `github.ler <https://github.com/hemantaph/ler/>`_. Use conda environment to avoid dependency error. 

    
.. tabs::

     .. code-tab:: console with new conda env

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
Installation of numba and healpy with conda
=======================

.. code-block:: console

    conda install -c conda-forge numba
    conda install -c conda-forge healpy
    
.. _gwcosmo:
Installation from the source for gwcosmo
=======================
    
.. tabs::

     .. code-tab:: console pip

        pip install git+https://git.ligo.org/lscsoft/gwcosmo.git@v1.0.0
        
     .. code-tab:: console git
     
        git clone https://git.ligo.org/lscsoft/gwcosmo.git
        cd gwcosmo
        pip install -r requirements.txt
        pip install .