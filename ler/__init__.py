"""
ler: LVK (LIGO-Virgo-KAGRA collaboration) Event (compact-binary mergers) Rate calculator and simulator
"""
# In your package __init__.py
import os
import multiprocessing as mp
import warnings
# disable OpenMP warnings
os.environ['KMP_WARNINGS'] = 'off'

def set_multiprocessing_start_method():
    """
    Set the multiprocessing start method based on OS and environment variables.
    By default, sets 'fork' on POSIX systems unless overridden by LER_USE_SPAWN=True.
    """
    # Only set if not already set
    method = None

    # POSIX = Linux, Mac
    if os.name == 'posix':
        # User override: LER_USE_SPAWN=True
        use_spawn = os.environ.get("LER_USE_SPAWN", "False").lower() == "true"
        if use_spawn:
            method = "spawn"
            print("ler: Setting multiprocessing start method to 'spawn' per environment variable.")
        else:
            method = "fork"
            # print(
            #     "ler: Setting multiprocessing start method to 'fork'.\n"
            #     "If you need to use the 'spawn' method (in case error or warning due to other library dependencies),\n"
            #     "set the environment variable LER_USE_SPAWN=True *before* running your script."
            #     "\n"
            #     "Command line (single line):\n"
            #     "    LER_USE_SPAWN=True python yourscript.py\n"
            #     "In a Python script (before importing ler):\n"
            #     "    import os\n"
            #     "    os.environ['LER_USE_SPAWN'] = 'True'\n"
            #     "    import ler\n"
            # )
        try:
            mp.set_start_method(method, force=True)
        except RuntimeError:
            # Already set (possibly by another library, or called too late)
            warnings.warn(
                f"ler: Could not set multiprocessing start method to '{method}'. "
                "This is usually because the start method was already set by another library or your script, "
                "or because it is being set too late.\n"
                "If you need to control the multiprocessing start method, set it at the very beginning of your script.\n"
                "To use the 'spawn' method instead, set the environment variable LER_USE_SPAWN=True before running your script.\n"
                "Command line (single line):\n"
                "    LER_USE_SPAWN=True python yourscript.py\n"
                "In a Python script (before importing ler):\n"
                "    import os\n"
                "    os.environ['LER_USE_SPAWN'] = 'True'\n"
                "    import ler\n"
            )
    else:
        # For Windows, default is already 'spawn', nothing to do.
        pass

# Call the function on package import
set_multiprocessing_start_method()

import warnings
import logging
from . import rates
from .rates import LeR
from .rates import GWRATES
from ._version import __version__

# Package metadata
__all__ = ['LeR', 'rates']
__author__ = 'Hemantakumar Phurailatpam <hemantaphurailatpam@gmail.com>'
__license__ = "MIT"
__email__ = "hemantaphurailatpam@gmail.com"
__maintainer__ = "Hemantakumar Phurailatpam"
__status__ = "Development"
__url__ = "https://github.com/hemantaph/ler"
__description__ = "ler is a statistical based python package whose core function is to calculate detectable rates of gravitational waves events (both lensed and un-lensed events)."
__version_info__ = tuple(map(int, __version__.split('.')))

warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")
logging.getLogger(__name__).addHandler(logging.NullHandler())


