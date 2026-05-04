"""
Gravitational-wave lensed-event rate code for LVK compact binary coalescences.

``import ler`` only configures lightweight threading defaults; the main classes
``LeR`` and ``GWRATES`` are loaded on first access (see ``__getattr__``).
"""

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
    Set ``multiprocessing`` start method once per process when explicitly called.

    Default choices: ``spawn`` on macOS, ``fork`` on other POSIX systems. Windows
    is left unchanged (``spawn``).

    Environment overrides (POSIX): ``LER_USE_SPAWN=True`` or ``LER_USE_FORK=True``
    (if both are set, ``spawn`` is used and a warning is issued).

    Returns
    -------
    None
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


import logging
from ._version import __version__

# Package metadata
__all__ = [
    "LeR",
    "GWRATES",
    "rates",
    "set_multiprocessing_start_method",
    "__version__",
]
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

def __getattr__(name):
    """
    Resolve ``LeR``, ``GWRATES``, or the ``rates`` subpackage on first use.

    Parameters
    ----------
    name : str
        Attribute name on the ``ler`` package.

    Returns
    -------
    type or module
        ``LeR`` or ``GWRATES`` from ``ler.rates``, or the ``ler.rates`` package.

    Raises
    ------
    AttributeError
        If ``name`` is not one of the names exported via this hook.
    """
    if name == "rates":
        from . import rates as _rates
        globals()[name] = _rates
        return _rates
    if name == "LeR":
        from .rates import LeR as _LeR
        globals()[name] = _LeR
        return _LeR
    if name == "GWRATES":
        from .rates import GWRATES as _GWRATES
        globals()[name] = _GWRATES
        return _GWRATES
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


# Multiprocessing entry points (see package internals):
# - ler.lens_galaxy_population.optical_depth.cross_section_epl_shear_numerical_mp
#   (lenstronomy EPL+shear caustics for double-image region).
# - ler.lens_galaxy_population.optical_depth._lens_redshift_multiprocessing
#   (ler.lens_galaxy_population.mp.lens_redshift_strongly_lensed_mp).
# Numba prange is used in, e.g., ler.image_properties.cross_section_njit,
# epl_shear_njit, sample_caustic_points_njit; cross_section_interpolator;
# lens_galaxy_population.mp / sampler_functions. Nested prange may serialize
# inside Numba depending on build and runtime.
