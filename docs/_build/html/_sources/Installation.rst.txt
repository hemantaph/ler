============
Installation
============

.. note::

    ``ler`` supports Python 3.10+ (but 3.11 recommended) and utilizes multi-core CPUs and multi-threaded CPUs when available.
    
    For package development and contribution refer here (:ref:`development`).

.. tabs::

   .. code-tab:: bash pip (standard)

        pip install ler

   .. code-tab:: bash pip (with gwsnr JAX backend)

      pip install ler jax jaxlib
      pip install -U "jax[cuda12]" # optional, for Nvidia GPU support

   .. code-tab:: bash pip (with gwsnr MLX backend; Apple Silicon)

      pip install ler mlx

   .. code-tab:: bash pip (with tensorflow based ANN)

      pip install ler scikit-learn tensorflow
      pip install --upgrade ml-dtypes # optional, for compatibility


This will also install the dependencies needed by the lastest ``ler`` version.  

.. _development:
ler for development
======================

To install ``ler`` for development purposes use `github.ler <https://github.com/hemantaph/ler/>`_. Use conda environment to avoid dependency error. 

    
.. tabs::

     .. code-tab:: bash with new conda env

        git clone https://github.com/hemantaph/ler.git
        cd ler
        conda env create -f ler.yml
        conda activate ler
        pip install -e .
        
     .. code-tab:: bash with existing conda env
     
        git clone https://github.com/hemantaph/ler.git
        cd ler
        conda env update --file ler.yml
        pip install -e .
    
.. _dependencies:
Installation of numba with conda
================================

.. code-block:: bash

    conda install -c conda-forge numba

.. note::

    For arm64 architecture processor (e.g. Apple Silicon) you might need to install some specific dependencies with conda. There can be errors related to multiprocessing on either arm64 or x86_64 architecture. This can be resolved by using the `main()` function as shown in this `example <https://github.com/hemantaph/ler/blob/main/examples/rates/LeR_short_example.py>`_.