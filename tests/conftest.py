"""
Pytest configuration and shared fixtures for the ler test suite.

This module sets up session-scoped fixtures that are shared across all tests.
Most importantly, it ensures that the bundled interpolator data is extracted
to a session cache directory once per test session, avoiding repeated costly
interpolator construction. The cache directory is intentionally persistent so
that `interpolator_json` is not deleted after the test run.

Wall-clock / JIT comparison tests, heavy cross-section samplers, and the LeR
speed benchmark subprocess checks are marked ``@pytest.mark.slow`` and are
deselected from the default collection.
"""

import os
import re
from pathlib import Path
import zipfile
import pytest
from importlib_resources import files as resources_files


_TEST_CACHE_DIR = Path(__file__).resolve().parent / ".cache"
_NUMBA_CACHE_DIR = _TEST_CACHE_DIR / "numba"
_MPL_CACHE_DIR = _TEST_CACHE_DIR / "matplotlib"
_NUMBA_CACHE_DIR.mkdir(parents=True, exist_ok=True)
_MPL_CACHE_DIR.mkdir(parents=True, exist_ok=True)

# Some optional backends imported during collection enable Numba caching at
# import time. Provide cache locations before test modules import those backends.
os.environ.setdefault("NUMBA_CACHE_DIR", str(_NUMBA_CACHE_DIR))
os.environ.setdefault("MPLCONFIGDIR", str(_MPL_CACHE_DIR))


def pytest_addoption(parser):
    """Register ``--run-slow`` so default runs can omit ``slow`` tests."""
    parser.addoption(
        "--run-slow",
        action="store_true",
        default=False,
        help=(
            "Include @pytest.mark.slow tests in the default suite (wall-clock / heavy benchmarks). "
            "Not needed when selecting them explicitly (e.g. ``pytest -m slow``)."
        ),
    )


def _mark_expression_requests_slow(markexpr):
    """
    Return True if ``-m markexpr`` is intended to select ``slow`` tests.

    Default (empty expression) returns False so we deselect ``slow`` tests.
    Expressions containing ``not slow`` return False.
    """
    m = (markexpr or "").strip()
    if not m:
        return False
    if re.search(r"\bnot\s+slow\b", m):
        return False
    return bool(re.search(r"\bslow\b", m))


@pytest.fixture(scope="session")
def interpolator_directory(tmp_path_factory):
    """
    Session-scoped fixture that extracts the bundled interpolator data to a
    cache directory. The directory is reused for the entire test session and
    kept on disk after the tests finish.

    Returns
    -------
    interpolator_directory : str
        Absolute path to the extracted interpolator_json directory.
    """
    # Use a persistent cache directory by default so interpolator_json is not
    # deleted when the pytest temporary directory is cleaned up.
    #
    # Users can override the location if needed (e.g. on shared runners).
    cache_root = os.environ.get("LER_TEST_INTERPOLATOR_CACHE_DIR", "").strip()
    if cache_root:
        cache_dir = Path(cache_root).expanduser().resolve()
    else:
        cache_dir = Path(__file__).resolve().parent / ".cache" / "interpolators"
    cache_dir.mkdir(parents=True, exist_ok=True)

    interpolator_json_dir = cache_dir / "interpolator_json"
    if interpolator_json_dir.exists():
        return str(interpolator_json_dir)

    # NOTE: use the top-level package as the anchor to avoid importing
    # heavy optional dependencies that may be pulled in by `ler.rates`.
    zip_resource = resources_files("ler").joinpath(
        "rates", "ler_data", "interpolator_json.zip"
    )
    with zip_resource.open("rb") as zip_file:
        with zipfile.ZipFile(zip_file, "r") as zip_ref:
            for member in zip_ref.namelist():
                if not member.startswith("__MACOSX"):
                    zip_ref.extract(member, cache_dir)

    return str(interpolator_json_dir)


@pytest.fixture(scope="session")
def ler_directory(tmp_path_factory):
    """
    Session-scoped output directory shared by unit and integration tests.

    Returns
    -------
    ler_directory : str
        Absolute path to a temporary ler_data directory.
    """
    return str(tmp_path_factory.mktemp("ler_data", numbered=False))


@pytest.fixture(scope="session")
def interpolator_dir(interpolator_directory):
    """Backward-compatible alias for older tests."""
    return interpolator_directory


def pytest_collection_modifyitems(session, config, items):
    """Deterministic module order; default runs omit ``slow`` tests."""
    module_order = {
        "test_cbc_source_redshift_distribution.py": 1,
        "test_cbc_source_parameter_distribution.py": 2,
        "test_gwrates.py": 3,
        "test_optical_depth.py": 4,
        "test_lens_galaxy_parameter_distribution.py": 5,
        "test_image_properties.py": 6,
        "test_ler.py": 7,
    }

    def _sort_key(item):
        file_name = item.fspath.basename
        return (module_order.get(file_name, 100), str(item.fspath), item.name)

    items.sort(key=_sort_key)

    if config.getoption("--run-slow"):
        return
    markexpr = getattr(config.option, "markexpr", None) or ""
    if _mark_expression_requests_slow(markexpr):
        return

    deselected = [item for item in items if "slow" in item.keywords]
    items[:] = [item for item in items if "slow" not in item.keywords]
    if deselected:
        config.hook.pytest_deselected(items=deselected)
