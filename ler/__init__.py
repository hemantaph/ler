"""
ler: LVK (LIGO-Virgo-KAGRA collaboration) Event (compact-binary mergers) Rate calculator and simulator
"""

# In your package __init__.py
import os
import multiprocessing as mp
import warnings

# disable OpenMP warnings
os.environ["KMP_WARNINGS"] = "off"

# Set default thread counts before numba/MKL/OpenBLAS are imported.
# NUMBA_NUM_THREADS: controls Numba prange parallelism (default: all cores).
#   Dynamically adjusted to npool via numba.set_num_threads() at call sites,
#   and set to 1 inside multiprocessing workers to avoid CPU oversubscription.
# OMP/MKL/OPENBLAS: default to 1 so BLAS libraries don't compete with
#   Numba's own thread pool or multiprocessing workers.
os.environ.setdefault("NUMBA_NUM_THREADS", str(os.cpu_count()))
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")


def set_multiprocessing_start_method():
    """
    Set the multiprocessing start method based on OS and environment variables.
    Defaults:
    - macOS: 'spawn' (safer with threaded native libraries)
    - Linux/other POSIX: 'fork'

    Overrides:
    - LER_USE_SPAWN=True forces 'spawn'
    - LER_USE_FORK=True forces 'fork'
    """
    method = None

    if os.name == "posix":
        use_spawn = os.environ.get("LER_USE_SPAWN", "False").lower() == "true"
        use_fork = os.environ.get("LER_USE_FORK", "False").lower() == "true"

        if use_spawn and use_fork:
            warnings.warn(
                "ler: Both LER_USE_SPAWN=True and LER_USE_FORK=True are set. "
                "Using 'spawn'."
            )
            method = "spawn"
            print(
                "ler: Setting multiprocessing start method to 'spawn' per environment variable."
            )
        elif use_spawn:
            method = "spawn"
            print(
                "ler: Setting multiprocessing start method to 'spawn' per environment variable."
            )
        elif use_fork:
            method = "fork"
            print(
                "ler: Setting multiprocessing start method to 'fork' per environment variable."
            )
        else:
            method = "spawn" if os.uname().sysname == "Darwin" else "fork"

        current = mp.get_start_method(allow_none=True)
        if current is None:
            try:
                mp.set_start_method(method)
            except RuntimeError:
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
        elif current != method:
            warnings.warn(
                f"ler: Multiprocessing start method is already '{current}'. "
                f"Requested '{method}' was not applied."
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
__all__ = ["LeR", "rates"]
__author__ = "Hemantakumar Phurailatpam <hemantaphurailatpam@gmail.com>"
__license__ = "MIT"
__email__ = "hemantaphurailatpam@gmail.com"
__maintainer__ = "Hemantakumar Phurailatpam"
__status__ = "Development"
__url__ = "https://github.com/hemantaph/ler"
__description__ = "ler is a statistical based python package whose core function is to calculate detectable rates of gravitational waves events (both lensed and un-lensed events)."
__version_info__ = tuple(map(int, __version__.split(".")))

warnings.filterwarnings("ignore", "Wswiglal-redir-stdio")
logging.getLogger(__name__).addHandler(logging.NullHandler())


"""
Where multiprocessing is used:
1. ler.lens_galaxy_population.optical_depth.cross_section_epl_shear_numerical_mp
  - uses lenstronomy.LensModel.Solver.epl_shear_solver.caustics_epl_shear for generate double caustic boundary
2. ler.lens_galaxy_population.optical_depth._lens_redshift_multiprocessing
  - uses njitted sampling functions (see ler.lens_galaxy_population.mp.lens_redshift_strongly_lensed_mp)

Where njitted prange is used:
1. ler.image_properties.cross_section_njit.make_cross_section_reinit
  - lenstronomy style njitted cross section calculation
2. ler.image_properties.epl_shear_njit.create_epl_shear_solver
  - lenstronomy style njitted lens equation solver
3. ler.image_properties.sample_caustic_points_njit._points_in_poly_precomp
  - 
4. ler.lens_galaxy_population.cross_section_interpolator._map_coordinates_5d_cubic_nearest
  - cross section interpolator 
5. ler.lens_galaxy_population.mp.lens_redshift_strongly_lensed_njit
  - JIT-compiled parallel computation of differential optical dept (lens redshift)
6. ler.lens_galaxy_population.sampler_functions.importance_sampler
  - importance sampling for lens galaxy parameters

There can be nested prange in the code, but I believe that Numba would silently serialize nested prange.
"""
