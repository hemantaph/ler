============
Installation
============

.. note::

    It is recommended to install ``gwsnr`` and ``gwcosmo`` before installing ``ler``. See the `gwsnr documentation <https://github.com/hemantaph/gwsnr/>`_ . Other dependencies are listed `here <https://github.com/hemantaph/ler/blob/main/setup.py>`_ in the setup.py, which will be automatically installed with ``ler``. 

.. tabs::

     .. code-tab:: console pip

        pip install ler gwsnr


This will install the dependencies needed but the version of gwsnr will not necessarily be the correct version for your system.


ler for development
======================

To install ``ler`` for development purposes use <https://github.com/hemantaph/ler/>`_.

Installation of gwcosmo
=======================
See the ```gwcosmo`` documentation <https://git.ligo.org/lscsoft/gwcosmo/>`_ . It is use for generating component masses of PopI/PopII binary blackholes.

.. code-block:: console

   $ git clone https://git.ligo.org/lscsoft/gwcosmo.git
   $ cd gwcosmo/
   $ pip install -r requirements.txt
   $ pip install .
   
Installing gwsnr from the source
=======================

.. code-block:: console

   $ git clone git@github.com:hemantaph/gwsnr.git
   $ cd gwsnr/
   $ pip install -r requirements.txt
   $ pip install .
   
   
ler for development
======================

.. code-block:: console

   $ git clone git@github.com:hemantaph/gwsnr.git
   $ cd gwsnr/
   $ pip install -r requirements.txt
   $ pip install .
   
   
ler for development
======================
To install ``ler`` for development purposes use <https://github.com/hemantaph/ler/>`_.