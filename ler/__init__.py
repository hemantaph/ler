"""
LeR
"""

import sys
import multiprocessing as mp

try:
   mp.set_start_method('spawn', force=True)
except RuntimeError:
   pass
#from . import rates, gw_source_population, lens_galaxy_population, image_properties, utils

__author__ = 'hemanta_ph <hemantaphurailatpam@gmail.com>'

__version__ = "0.2.7"


