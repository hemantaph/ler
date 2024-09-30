"""
LeR
"""

# mypackage/cli.py
import argparse
# import subprocess, os, sys, signal, warnings

## import pycbc
import os
import multiprocessing as mp

def set_multiprocessing_start_method():
    if os.name == 'posix':  # posix indicates the program is run on Unix/Linux/Mac
      print("Setting multiprocessing start method to 'fork'")
      try:
         mp.set_start_method('fork', force=True)
      except RuntimeError:
         # The start method can only be set once and must be set before any process starts
         pass
    else:
      print("Setting multiprocessing start method to 'spawn'")
      # For Windows and other operating systems, use 'spawn'
      try:
         mp.set_start_method('spawn', force=True)
      except RuntimeError:
         pass

set_multiprocessing_start_method()

# try:
#    mp.set_start_method('fork', force=True)
# except RuntimeError:
#    pass

# if sys.platform == 'darwin':
#     HAVE_OMP = False

#     # MacosX after python3.7 switched to 'spawn', however, this does not
#     # preserve common state information which we have relied on when using
#     # multiprocessing based pools.
#     import multiprocessing
#     if multiprocessing.get_start_method(allow_none=True) is None:
#         if hasattr(multiprocessing, 'set_start_method'):
#             multiprocessing.set_start_method('fork')
#     elif multiprocessing.get_start_method() != 'fork':
#         warnings.warn("PyCBC requires the use of the 'fork' start method for multiprocessing, it is currently set to {}".format(multiprocessing.get_start_method()))
# else:
#     HAVE_OMP = True

__author__ = 'hemanta_ph <hemantaphurailatpam@gmail.com>'

# The version as used in the setup.py
__version__ = "0.4.1"

# add __file__
import os
__file__ = os.path.abspath(__file__)

from . import rates, gw_source_population, lens_galaxy_population, image_properties, utils

from .rates import ler, gwrates
from .gw_source_population import cbc_source_parameter_distribution, cbc_source_redshift_distribution
from .lens_galaxy_population import lens_galaxy_parameter_distribution, optical_depth
from .image_properties import image_properties
from .utils import utils, plots


