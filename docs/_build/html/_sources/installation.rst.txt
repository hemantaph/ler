============
Installation
============

.. note::

    It is recommended to install ``gwcosmo`` and ``gwsnr`` before installing ``ler``. See the `gwsnr documentation <https://github.com/hemantaph/gwsnr/>`_ . You can install ``gwsnr`` with pip but ``gwcosmo`` needs custom :ref:`my-reference-label`.

.. tabs::

     .. code-tab:: console pip

        pip install ler gwsnr


This will install the dependencies needed but the version of gwsnr will not necessarily be the correct version for your system.


ler for development
======================

To install ``ler`` for development purposes use <https://github.com/hemantaph/ler/>`_.

.. code-block:: console

    git clone https://github.com/hemantaph/ler.git
    cd ler
    pip install . -e
    

.. _my-reference-label:
Installation from the source for gwcosmo
=======================

.. code-block:: console

    git clone https://git.ligo.org/lscsoft/gwcosmo.git
    cd gwcosmo
    pip install -r requirements.txt
    pip install .