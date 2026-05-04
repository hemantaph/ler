"""
Unit test utilities for the ler test suite.

This module provides shared helper classes and functions used across
multiple unit test files. It mirrors the pattern used in gwsnr for shared test helpers.
"""

import inspect
import time

import numpy as np

from ler.utils import FunctionConditioning


# ---------------------------------------------------------------------------
# Expected parameter keys (shared across unit + integration tests)
# ---------------------------------------------------------------------------
# Minimal GW columns before detection (unlensed or intrinsic CBC draw).
EXPECTED_GW_PARAM_KEYS = [
    "zs",
    "luminosity_distance",
    "mass_1_source",
    "mass_2_source",
    "mass_1",
    "mass_2",
    "geocent_time",
    "ra",
    "dec",
    "phase",
    "psi",
    "theta_jn",
]

EXPECTED_SPIN_PRECESSING_KEYS = [
    "a_1",
    "a_2",
    "tilt_1",
    "tilt_2",
    "phi_12",
    "phi_jl",
]

# Spin columns for default BBH with ``spin_precession=False``: only the magnitude
# pair at the start of ``EXPECTED_SPIN_PRECESSING_KEYS``.
EXPECTED_ALIGNED_MODE_SPIN_KEYS = list(EXPECTED_SPIN_PRECESSING_KEYS[:2])

# BBH with ``spin_zero=True``: these keys must be absent from ``gw_parameters_rvs``.
EXPECTED_SPIN_ZERO_FORBIDDEN_BBH_KEYS = EXPECTED_ALIGNED_MODE_SPIN_KEYS

# Intrinsic GW columns for ``unlensed_cbc_statistics`` / ``gw_cbc_statistics`` (before
# ``pdet_net``): choose the list that matches ``spin_zero`` / ``spin_precession``.
EXPECTED_UNLENSED_GW_KEYS_NO_SPIN = list(EXPECTED_GW_PARAM_KEYS)
EXPECTED_UNLENSED_GW_KEYS_ALIGNED_SPIN = (
    EXPECTED_GW_PARAM_KEYS + EXPECTED_ALIGNED_MODE_SPIN_KEYS
)
EXPECTED_UNLENSED_GW_KEYS_PRECESSING_SPIN = (
    EXPECTED_GW_PARAM_KEYS + EXPECTED_SPIN_PRECESSING_KEYS
)

# 1-D per-event lens columns after ``lensed_cbc_statistics`` when redundant
# mass/distance fields have been stripped (default ``include_redundant_parameters=False``).
EXPECTED_LENSED_PARAM_KEYS = [
    "zl",
    "zs",
    "sigma",
    "q",
    "phi",
    "gamma",
    "gamma1",
    "gamma2",
    "mass_ratio",
    "geocent_time",
    "ra",
    "dec",
    "phase",
    "psi",
    "theta_jn",
    "mass_1",
    "mass_2",
]

# Same as ``EXPECTED_LENSED_PARAM_KEYS`` plus aligned-spin magnitudes ``a_1``, ``a_2``
# (default BBH with ``spin_zero=False``, ``spin_precession=False``).
EXPECTED_LENSED_PARAM_KEYS_ALIGNED_SPIN = (
    EXPECTED_LENSED_PARAM_KEYS + EXPECTED_ALIGNED_MODE_SPIN_KEYS
)

# Per-image columns: shape ``(n_events, n_max_images)``; pads may be NaN (use
# ``_assert_param_dict_valid(..., nan_to_num=True)`` in tests).
EXPECTED_IMAGE_KEYS = [
    "x0_image_positions",
    "x1_image_positions",
    "magnifications",
    "time_delays",
    "image_type",
]

# Source position in the lens plane (not per-image columns); 1-D length ``n_events``.
EXPECTED_SOURCE_POS_KEYS = ["x_source", "y_source"]

# Present only when ``include_redundant_parameters=True`` (retained through the pipeline).
EXPECTED_REDUNDANT_KEYS = [
    "theta_E",
    "n_images",
    "mass_1_source",
    "mass_2_source",
    "luminosity_distance",
]

# Per-image effective GW parameters when ``include_effective_parameters=True`` (2-D).
EXPECTED_EFFECTIVE_KEYS = [
    "effective_luminosity_distance",
    "effective_geocent_time",
    "effective_phase",
    "effective_ra",
    "effective_dec",
]

EXPECTED_PDET_NET_KEYS = ["pdet_net"]

# Unlensed CBC output = intrinsic GW columns + scalar ``pdet_net`` per event.
EXPECTED_UNLENSED_PARAM_KEYS = EXPECTED_GW_PARAM_KEYS + EXPECTED_PDET_NET_KEYS


def median_call_time(fn, repeats=3):
    """Median wall time in seconds for ``fn()`` over ``repeats`` calls.

    Used by performance-style tests to reduce sensitivity to CPU noise.
    """
    times = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        fn()
        times.append(time.perf_counter() - t0)
    return float(np.median(np.asarray(times)))


class CommonTestUtils:
    """Shared utilities and validators for ler unit tests."""

    def _assert_array_valid(
        self,
        arr,
        name="array",
        size=None,
        positive=False,
        finite=True,
        lo=None,
        hi=None,
        nan_to_num=False,
    ):
        """Assert that a numpy array satisfies common validity checks.

        Parameters
        ----------
        arr : array-like
            The array to validate.
        name : str
            Label used in assertion messages.
        size : int or tuple of int or None
            Expected 1-D length (``int``) or full ``shape`` tuple for ND arrays.
        positive : bool
            If True, assert all values are >= 0.
        finite : bool
            If True, assert all values are finite (no NaN or Inf).
        lo : float or None
            Lower bound (inclusive) for array values.
        hi : float or None
            Upper bound (inclusive) for array values.
        nan_to_num : bool
            If True, apply :func:`numpy.nan_to_num` before other checks.
        """

        arr = np.asarray(arr)
        if nan_to_num:
            arr = np.nan_to_num(arr)

        assert arr.ndim >= 1, f"{name}: expected at least 1-D array"
        if size is not None:
            if isinstance(size, int):
                assert arr.shape == (size,), f"{name}: expected length {size}, got {arr.shape}"
            else:
                assert arr.shape == size, f"{name}: expected shape {size}, got {arr.shape}"
        if finite:
            assert np.all(np.isfinite(arr)), f"{name}: contains NaN or Inf"
        if positive:
            assert np.all(arr >= 0), f"{name}: contains negative values"
        if lo is not None:
            assert np.all(arr >= lo), f"{name}: values below lower bound {lo}"
        if hi is not None:
            assert np.all(arr <= hi), f"{name}: values above upper bound {hi}"

    def _assert_param_dict_valid(
        self,
        param_dict,
        expected_keys,
        size=None,
        *,
        finite=True,
        nan_to_num=False,
        positive=False,
        lo=None,
        hi=None,
    ):
        """Validate a GW / lensed parameter dictionary output.

        Parameters
        ----------
        param_dict : dict
            Dictionary returned by a parameter sampler.
        expected_keys : list of str
            Keys that must be present in the dictionary.
        size : int or tuple of int or None
            Expected shape: ``int`` → 1-D length ``(size,)``; ``tuple`` → exact
            ``ndarray.shape``. If ``None``, shape is not checked.
        finite : bool, optional
            If True (default), all values must be finite (after ``nan_to_num`` if set).
        nan_to_num : bool, optional
            If True, run :func:`numpy.nan_to_num` on each array before finiteness checks
            (typical for padded per-image columns and ``pdet_net``).
        positive : bool, optional
            Passed to :meth:`_assert_array_valid` (e.g. ``pdet_net`` >= 0).
        lo, hi : float or None, optional
            Inclusive bounds passed to :meth:`_assert_array_valid` (e.g. ``pdet_net`` in ``[0, 1]``).
        """
        assert isinstance(param_dict, dict), "output must be a dict"
        for key in expected_keys:
            assert key in param_dict, f"missing key: {key}"
            self._assert_array_valid(
                param_dict[key],
                name=key,
                size=size,
                finite=finite,
                nan_to_num=nan_to_num,
                positive=positive,
                lo=lo,
                hi=hi,
            )

    def _assert_unlensed_cbc_outputs(self, param_dict, size, *, expected_gw_keys=None):
        """Intrinsic GW columns and 1-D ``pdet_net`` in ``[0, 1]``.

        Pass one of :data:`EXPECTED_UNLENSED_GW_KEYS_NO_SPIN`,
        :data:`EXPECTED_UNLENSED_GW_KEYS_ALIGNED_SPIN`, or
        :data:`EXPECTED_UNLENSED_GW_KEYS_PRECESSING_SPIN`.  Default is no-spin
        (``EXPECTED_GW_PARAM_KEYS`` only).
        """
        if expected_gw_keys is None:
            gw_keys = EXPECTED_UNLENSED_GW_KEYS_NO_SPIN
        else:
            gw_keys = list(expected_gw_keys)
        self._assert_param_dict_valid(param_dict, gw_keys, size)
        self._assert_param_dict_valid(
            param_dict,
            EXPECTED_PDET_NET_KEYS,
            size,
            positive=True,
            lo=0.0,
            hi=1.0,
        )

    def _assert_lensed_pdet_net_2d(self, param_dict, shape2):
        """2-D ``pdet_net`` with NaN padding; probability values in ``[0, 1]`` after ``nan_to_num``."""
        self._assert_param_dict_valid(
            param_dict,
            EXPECTED_PDET_NET_KEYS,
            shape2,
            nan_to_num=True,
            positive=True,
            lo=0.0,
            hi=1.0,
        )

    def _assert_valid_object(
        self,
        object,
        x_array=None,
        y_array=None,
        name="object",
        size=None,
        positive=False,
        finite=True,
        lo=None,
        hi=None,
        pdf=False,
        rvs=False,
        function=False,
    ):
        """Validate a ``FunctionConditioning`` object interface.

        This helper is used across unit tests to ensure that a sampler/interpolator
        object exposes the expected callables and that they return outputs with
        basic invariants.

        Parameters
        ----------
        object : ``FunctionConditioning``
            Object under test.
        x_array, y_array : array-like or None
            Evaluation points. If ``y_array`` is provided, the callable is assumed
            to be conditioned (2D case), e.g. ``pdf(x, y)``.
        name : str
            Label used in assertion messages.
        size : int or None
            Expected output size. If not given, inferred from ``x_array``.
        positive : bool
            If True, require outputs from selected callables to be >= 0.
            Note: this must be False for distributions with support on negative
            values (e.g. normal-distributed shear components).
        finite : bool
            If True, require outputs to be finite.
        lo, hi : float or None
            Optional inclusive bounds check for ``rvs`` outputs.
        pdf, rvs, function : bool
            Which callables to validate.
        """

        # if size is not None, size==len(x_array) or size==len(y_array)
        if size is not None:
            if x_array is not None:
                assert len(x_array) == size, f"{name}: expected length {size}, got {len(x_array)}"
            if y_array is not None:
                assert len(y_array) == size, f"{name}: expected length {size}, got {len(y_array)}"
        else:
            if x_array is not None:
                size = len(x_array)
                self._assert_array_valid(x_array, name=f"{name}.x_array", size=size)
            else:
                raise ValueError(f"{name}: size is None and x_array is None")

            if y_array is not None:
                self._assert_array_valid(y_array, name=f"{name}.y_array", size=size)

        # check if object is an instance of FunctionConditioning
        assert isinstance(object, FunctionConditioning), f"{name}: expected instance of FunctionConditioning"

        if y_array is None:
            if pdf:
                self._assert_array_valid(object.pdf(x_array), name=f"{name}.pdf", size=size, positive=positive, finite=finite)
            if rvs:
                self._assert_array_valid(
                    object.rvs(size),
                    name=f"{name}.rvs",
                    size=size,
                    positive=positive,
                    finite=finite,
                    lo=lo,
                    hi=hi,
                )
            if function:
                self._assert_array_valid(
                    object.function(x_array),
                    name=f"{name}.function",
                    size=size,
                    positive=positive,
                    finite=finite,
                )

        else:
            if pdf:
                self._assert_array_valid(object.pdf(x_array, y_array), name=f"{name}.pdf", size=size, positive=positive, finite=finite)
            if rvs:
                self._assert_array_valid(
                    object.rvs(size, y_array),
                    name=f"{name}.rvs",
                    size=size,
                    positive=positive,
                    finite=finite,
                    lo=lo,
                    hi=hi,
                )
            if function:
                self._assert_array_valid(
                    object.function(x_array, y_array),
                    name=f"{name}.function",
                    size=size,
                    positive=positive,
                    finite=finite,
                )

    def _arg_count(self, function):
        """Return the number of positional parameters accepted by `function`.

        Supports regular Python callables and numba-jitted callables (via `py_func`).
        """
        if function is None:
            return 0

        fn = function.py_func if hasattr(function, "py_func") else function
        try:
            sig = inspect.signature(fn)
        except (TypeError, ValueError):
            if hasattr(fn, "__code__"):
                return int(fn.__code__.co_argcount)
            return 0

        count = 0
        for p in sig.parameters.values():
            if p.kind in (inspect.Parameter.POSITIONAL_ONLY, inspect.Parameter.POSITIONAL_OR_KEYWORD):
                count += 1
        return count

    def _assert_param_dicts_equal(self, d1, d2, label=""):
        """Compare two parameter dicts whose values are numpy arrays."""
        prefix = f"{label}: " if label else ""
        assert d1.keys() == d2.keys(), f"{prefix}key mismatch: {d1.keys()} vs {d2.keys()}"
        for k in d1:
            np.testing.assert_array_equal(d1[k], d2[k], err_msg=f"{prefix}key '{k}' differs")
