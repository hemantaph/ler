"""
Integration benchmark: JIT ``LeR`` (parallel ``npool``, intent ``~6``, clamped per host) vs no-JIT subprocess.

The no-JIT branch mirrors the slow backend configuration used by the speed
comparison example:

- Subprocess ``NUMBA_DISABLE_JIT=1``, ``npool=1``, ``rejection_sampler_partial``,
  ``cross_section_epl_shear_njit``, ``image_properties_epl_shear_lenstronomy``.
- Timing (each branch): three passes; **mean** wall time across those three for
  **unlensed** and **lensed** separately, with::

    unlensed_cbc_statistics(size=100, batch_size=100, ...)
    lensed_cbc_statistics(size=10, batch_size=10, ...)

Fast branch: default ``LeR`` with ``npool=6``, same ``size`` / ``batch_size``.
Mock ``pdet_finder`` only. Reports averages; requires **lensed** no-JIT mean wall time
>= JIT (the heavy lenstronomy / no-JIT path). **Unlensed** times are printed but not
ordered: with ``npool=6``, JIT pays multiprocessing overhead while ``NUMBA_DISABLE_JIT``
can make the unlensed draw+mock-``pdet`` path faster at small sizes, so
``no_jit >= jit`` is not reliable there.

The parent pytest process **never** imports ``LeR`` / Numba for the timed work.
Both the JIT (``npool=6``) and no-JIT paths run as **fresh subprocesses**. If
``NUMBA_THREADING_LAYER`` is unset, a short subprocess probe picks the first
loadable backend among **tbb** -> **omp**. With ``NPOOL_JIT > 1``, **workqueue**
is **not** auto-selected: initializing it avoids the earlier ``omp`` ImportError-style
failure, but the benchmark then often hits Numba's "Concurrent access has been detected"
abort alongside ``ler`` multiprocessing. The benchmark requires a loadable
**tbb** or **omp** threading backend, or an explicit ``NUMBA_THREADING_LAYER``.

Each subprocess writes two **stdout** lines with tokens ``LER_SPEED_BENCHMARK_UNLENSED_MEAN_S``
and ``LER_SPEED_BENCHMARK_LENSED_MEAN_S`` (``flush=True``). ``ler`` progress is redirected
mostly to stderr, so timings are parsed from subprocess stdout rather than stderr.

Each subprocess also gets per-run ``NUMBA_CACHE_DIR`` and ``MPLCONFIGDIR`` values under
``ler_directory/ler_speed_benchmark``. That keeps import-time cache discovery and
Matplotlib cache writes away from user-level locations during the benchmark.

Marked ``slow`` and deselected from the default pytest collection.
"""

import os
import subprocess
import sys

import numpy as np
import pytest

from tests_utils import clamp_npool_for_numba


@pytest.mark.slow
class TestLeRSpeedBenchmark:
    """Mean wall times (3 passes) for unlensed vs lensed, JIT vs no-JIT subprocesses."""

    UNLENSED_SIZE = 100
    UNLENSED_BATCH = 100
    LENSED_SIZE = 10
    LENSED_BATCH = 10
    N_AVG = 3
    NPOOL_JIT = clamp_npool_for_numba(6)

    @classmethod
    def _numba_threading_layer_for_subprocess_env(cls):
        """
        Pick a Numba threading layer before ``LeR`` calls ``set_num_threads``.

        ``LeR`` sets ``npool`` via ``numba.set_num_threads``, which requires a
        loadable threading layer.

        When ``NPOOL_JIT > 1``, **workqueue** is not considered: probes pass, but the
        real benchmark subprocess often exits with concurrent-access abort.

        Returns
        -------
        str or None
            Layer name to set on the subprocess environment, or ``None`` if
            ``NUMBA_THREADING_LAYER`` is already defined (respect user / CI).
        """
        if os.environ.get("NUMBA_THREADING_LAYER"):
            return None
        # Only tbb / omp when JIT uses multiprocessing pools (see module docstring).
        candidates = ("tbb", "omp") if cls.NPOOL_JIT > 1 else ("tbb", "omp", "workqueue")
        for cand in candidates:
            probe = (
                "import os\n"
                f"os.environ['NUMBA_THREADING_LAYER'] = {cand!r}\n"
                "from numba import set_num_threads\n"
                "set_num_threads(2)\n"
            )
            proc = subprocess.run(
                [sys.executable, "-c", probe],
                capture_output=True,
                text=True,
                timeout=120,
            )
            if proc.returncode == 0:
                return cand
        pytest.skip(
            (
                "No loadable Numba threading layer "
                + ("(tried tbb, omp, workqueue). " if cls.NPOOL_JIT <= 1 else "(tried tbb, omp). ")
                + (
                    "`workqueue` is skipped for JIT npool>1 — it tends to fatal-abort with "
                    "ler multiprocessing. " if cls.NPOOL_JIT > 1 else ""
                )
                + "Provide a Numba-backed threading layer or set NUMBA_THREADING_LAYER "
                "(see Numba threading-layer docs)."
            )
        )

    def test_ler_unlensed_lensed_wall_time_jit_vs_numba_disable_jit_subprocess(
        self,
        interpolator_directory,
        ler_directory,
    ):
        """Two subprocesses (JIT ``npool=6`` vs no-JIT backend kwargs); compare means."""
        base_out = os.path.join(ler_directory, "ler_speed_benchmark")
        jit_dir = os.path.join(base_out, "jit_npool_subprocess")
        sub_dir = os.path.join(base_out, "no_jit_subprocess")
        numba_cache_dir = os.path.join(base_out, "numba_cache")
        mpl_cache_dir = os.path.join(base_out, "matplotlib_cache")
        os.makedirs(jit_dir, exist_ok=True)
        os.makedirs(sub_dir, exist_ok=True)
        os.makedirs(numba_cache_dir, exist_ok=True)
        os.makedirs(mpl_cache_dir, exist_ok=True)

        jit_script = (
            """import time
import numpy as np

def mock_pdet_finder(gw_param_dict):
    n = len(list(gw_param_dict.values())[0])
    return dict(pdet_net=np.ones(n))

from ler.rates import LeR

ler = LeR(
    npool=""" + str(self.NPOOL_JIT) + """,
    pdet_finder=mock_pdet_finder,
    interpolator_directory=""" + repr(interpolator_directory) + """,
    ler_directory=""" + repr(jit_dir) + """,
    event_type="BBH",
    lens_type="epl_shear_galaxy",
    create_new_interpolator=False,
    verbose=False,
    spin_zero=False,
    spin_precession=False,
    multiprocessing_verbose=False,
)

times_u = []
times_l = []
for _ in range(""" + str(self.N_AVG) + """):
    t0 = time.perf_counter()
    _ = ler.unlensed_cbc_statistics(
        size=""" + str(self.UNLENSED_SIZE) + """,
        batch_size=""" + str(self.UNLENSED_BATCH) + """,
        resume=False,
        output_jsonfile=False,
    )
    times_u.append(time.perf_counter() - t0)

    t0 = time.perf_counter()
    _ = ler.lensed_cbc_statistics(
        size=""" + str(self.LENSED_SIZE) + """,
        batch_size=""" + str(self.LENSED_BATCH) + """,
        resume=False,
        save_batch=False,
        output_jsonfile=False,
    )
    times_l.append(time.perf_counter() - t0)

# Timings on stdout; ``ler`` chatty output is typically on stderr.
print("LER_SPEED_BENCHMARK_UNLENSED_MEAN_S", float(np.mean(np.asarray(times_u))), flush=True)
print("LER_SPEED_BENCHMARK_LENSED_MEAN_S", float(np.mean(np.asarray(times_l))), flush=True)
"""
        )

        nojit_script = (
            """import time
import numpy as np

def mock_pdet_finder(gw_param_dict):
    n = len(list(gw_param_dict.values())[0])
    return dict(pdet_net=np.ones(n))

from ler.rates import LeR

ler = LeR(
    npool=1,
    pdet_finder=mock_pdet_finder,
    interpolator_directory=""" + repr(interpolator_directory) + """,
    ler_directory=""" + repr(sub_dir) + """,
    event_type="BBH",
    lens_type="epl_shear_galaxy",
    create_new_interpolator=False,
    verbose=False,
    spin_zero=False,
    spin_precession=False,
    multiprocessing_verbose=False,
    lens_functions=dict(
        cross_section_based_sampler="rejection_sampler_partial",
        cross_section="cross_section_epl_shear_njit",
    ),
    image_properties_function="image_properties_epl_shear_lenstronomy",
)

times_u = []
times_l = []
for _ in range(""" + str(self.N_AVG) + """):
    t0 = time.perf_counter()
    _ = ler.unlensed_cbc_statistics(
        size=""" + str(self.UNLENSED_SIZE) + """,
        batch_size=""" + str(self.UNLENSED_BATCH) + """,
        resume=False,
        output_jsonfile=False,
    )
    times_u.append(time.perf_counter() - t0)

    t0 = time.perf_counter()
    _ = ler.lensed_cbc_statistics(
        size=""" + str(self.LENSED_SIZE) + """,
        batch_size=""" + str(self.LENSED_BATCH) + """,
        resume=False,
        save_batch=False,
        output_jsonfile=False,
    )
    times_l.append(time.perf_counter() - t0)

# Timings on stdout; ``ler`` chatty output is typically on stderr.
print("LER_SPEED_BENCHMARK_UNLENSED_MEAN_S", float(np.mean(np.asarray(times_u))), flush=True)
print("LER_SPEED_BENCHMARK_LENSED_MEAN_S", float(np.mean(np.asarray(times_l))), flush=True)
"""
        )

        layer = self._numba_threading_layer_for_subprocess_env()
        env_jit = dict(os.environ)
        env_jit.pop("NUMBA_DISABLE_JIT", None)
        env_jit.setdefault("NUMBA_CACHE_DIR", numba_cache_dir)
        env_jit.setdefault("MPLCONFIGDIR", mpl_cache_dir)
        if layer is not None:
            env_jit["NUMBA_THREADING_LAYER"] = layer

        env_nojit = dict(os.environ)
        env_nojit["NUMBA_DISABLE_JIT"] = "1"
        env_nojit.setdefault("NUMBA_CACHE_DIR", numba_cache_dir)
        env_nojit.setdefault("MPLCONFIGDIR", mpl_cache_dir)
        if layer is not None:
            env_nojit["NUMBA_THREADING_LAYER"] = layer

        proc_jit = subprocess.run(
            [sys.executable, "-c", jit_script],
            capture_output=True,
            text=True,
            env=env_jit,
            timeout=7200,
        )
        assert proc_jit.returncode == 0, (
            f"JIT subprocess failed (returncode={proc_jit.returncode}):\n"
            f"STDOUT:\n{proc_jit.stdout}\nSTDERR:\n{proc_jit.stderr}"
        )

        proc_nojit = subprocess.run(
            [sys.executable, "-c", nojit_script],
            capture_output=True,
            text=True,
            env=env_nojit,
            timeout=7200,
        )
        assert proc_nojit.returncode == 0, (
            f"no-JIT subprocess failed (returncode={proc_nojit.returncode}):\n"
            f"STDOUT:\n{proc_nojit.stdout}\nSTDERR:\n{proc_nojit.stderr}"
        )

        def _parse_benchmark_stdout(proc, label):
            unl, lens = None, None
            for ln in proc.stdout.splitlines():
                parts = ln.split(maxsplit=1)
                if len(parts) == 2 and parts[0] == "LER_SPEED_BENCHMARK_UNLENSED_MEAN_S":
                    unl = float(parts[1])
                elif len(parts) == 2 and parts[0] == "LER_SPEED_BENCHMARK_LENSED_MEAN_S":
                    lens = float(parts[1])
            assert unl is not None and lens is not None, (
                f"{label}: missing LER_SPEED_BENCHMARK_* lines on stdout:\n{proc.stdout!r}"
            )
            return unl, lens

        avg_unl_jit, avg_len_jit = _parse_benchmark_stdout(proc_jit, "JIT subprocess")
        avg_unl_nojit, avg_len_nojit = _parse_benchmark_stdout(proc_nojit, "no-JIT subprocess")

        for lbl, av in (
            ("jit unlensed", avg_unl_jit),
            ("jit lensed", avg_len_jit),
            ("nojit unlensed", avg_unl_nojit),
            ("nojit lensed", avg_len_nojit),
        ):
            assert np.isfinite(av) and av > 0.0, f"invalid avg {lbl}: {av}"

        print(
            f"\nLeR benchmark (mean of {self.N_AVG} passes, subprocess per branch):\n"
            f"  unlensed  size={self.UNLENSED_SIZE} batch={self.UNLENSED_BATCH}: "
            f"jit(npool={self.NPOOL_JIT})={avg_unl_jit:.4f}s  "
            f"nojit(npool=1,NUMBA_DISABLE_JIT)={avg_unl_nojit:.4f}s\n"
            f"  lensed    size={self.LENSED_SIZE} batch={self.LENSED_BATCH}: "
            f"jit(npool={self.NPOOL_JIT})={avg_len_jit:.4f}s  "
            f"nojit(npool=1,NUMBA_DISABLE_JIT)={avg_len_nojit:.4f}s\n",
        )

        assert avg_len_nojit >= avg_len_jit, (
            f"expected no-JIT lensed mean >= JIT (jit={avg_len_jit:.6f}s, "
            f"nojit={avg_len_nojit:.6f}s)"
        )
