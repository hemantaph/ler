"""
LeR
"""

#import pycbc
import multiprocessing as mp

try:
   mp.set_start_method('spawn', force=True)
except RuntimeError:
   pass
#from . import rates, gw_source_population, lens_galaxy_population, image_properties, utils

__author__ = 'hemanta_ph <hemantaphurailatpam@gmail.com>'

__version__ = "0.3.3"

from . import rates, gw_source_population, lens_galaxy_population, image_properties, utils

from .rates import ler, gwrates
from .gw_source_population import cbc_source_parameter_distribution, cbc_source_redshift_distribution
from .lens_galaxy_population import lens_galaxy_parameter_distribution, optical_depth
from .image_properties import image_properties
from .utils import utils, plots


