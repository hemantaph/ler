============
Installation
============

.. note::

    Install numba, healpy seperately with conda if you are with ARM processor, e.g. apple silicon. 

.. tabs::

     .. code-tab:: console pip

        pip install ler


This will also install the dependencies needed by the lastest ``ler`` version.  


ler for development
======================

To install ``ler`` for development purposes use `github.ler <https://github.com/hemantaph/ler/>`_.

.. code-block:: console

    git clone https://github.com/hemantaph/ler.git
    cd ler
    pip install . -e
    

Installation of numba and healpy with conda
=======================

.. code-block:: console

    conda install -c conda-forge numba
    conda install -c conda-forge healpy