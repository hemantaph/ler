# -*- coding: utf-8 -*-
"""
Module for function interpolation and conditioning using cubic splines and Gaussian KDE.

This module provides the FunctionConditioning class, which implements flexible function \n
interpolation for both 1D and 2D cases with support for inverse functions, probability \n
density functions (PDFs), and random variable sampling. It supports both cubic spline \n
interpolation and Gaussian Kernel Density Estimation (KDE).

Key Features: \n
- 1D and 2D cubic spline interpolation via JSON caching \n
- Automatic inverse function generation and normalization \n
- PDF computation with proper normalization constants \n
- Inverse transform sampling for random variable generation \n
- Gaussian KDE as an alternative to spline interpolation \n
- Numba JIT compilation for performance optimization \n

Copyright (C) 2026 Hemanta Ph. Distributed under MIT License.
"""

from .utils import (
    interpolator_json_path,
    cubic_spline_interpolator,
    pdf_cubic_spline_interpolator,
    cubic_spline_interpolator2d,
    inverse_transform_sampler,
    pdf_cubic_spline_interpolator2d,
    save_json,
    load_json,
    inverse_transform_sampler2d,
    cumulative_trapezoid,
    inverse_transform_sampler2d_spline,
    cumulative_spline,
)
import numpy as np
from scipy.interpolate import CubicSpline
from scipy.stats import gaussian_kde
from numba import njit


class FunctionConditioning:
    """
    Cubic spline interpolation and conditioning for functions with optional PDF and sampling.

    This class creates efficiently interpolated representations of scalar or multi-dimensional \n
    functions, supporting inverse functions, probability density functions, and random \n
    variable sampling. Interpolators are cached as JSON files to avoid recomputation. \n
    Both 1D and 2D conditioning scenarios are supported.

    Parameters
    ----------
    function : ``callable`` or ``numpy.ndarray`` or ``None``
        Function to interpolate or array of function values at x_array points. \n
        If callable, will be evaluated at x_array (and conditioned_y_array if provided). \n
        default: None
    x_array : ``numpy.ndarray``
        X-axis values where the function is defined. \n
        For 2D case, can be 1D (replicated for each y) or 2D (shape: len(conditioned_y_array), len(x)). \n
        default: None
    conditioned_y_array : ``numpy.ndarray`` or ``None``
        Y-axis values for 2D conditioning. If None, 1D interpolation is used. \n
        default: None
    y_array : ``numpy.ndarray`` or ``None``
        Y values for 2D Gaussian KDE (used only when gaussian_kde=True). \n
        default: None
    non_negative_function : ``bool``
        If True, clamps negative function values to 0. Useful for PDF-like functions. \n
        default: False
    gaussian_kde : ``bool``
        If True, uses Gaussian KDE instead of cubic spline interpolation. \n
        default: False
    gaussian_kde_kwargs : ``dict``
        Keyword arguments passed to scipy.stats.gaussian_kde constructor. \n
        default: None
    identifier_dict : ``dict``
        Dictionary to uniquely identify the interpolator in the cache. \n
        Used in directory structure for saving/loading. \n
        default: None
    directory : ``str``
        Root directory for storing interpolator JSON files. \n
        default: './interpolator_json'
    sub_directory : ``str``
        Subdirectory within directory for organizing interpolators. \n
        default: 'default'
    name : ``str``
        Name/identifier of this specific interpolator. \n
        default: 'default'
    create_new : ``bool``
        If True, forces creation of new interpolator, ignoring cached versions. \n
        default: False
    create_function : ``bool`` or ``callable``
        If True, compiles interpolated function via Numba JIT. \n
        If callable, uses provided function directly. \n
        default: False
    create_function_inverse : ``bool`` or ``callable``
        If True, compiles inverse interpolated function via Numba JIT. \n
        If callable, uses provided function directly. \n
        default: False
    create_pdf : ``bool`` or ``callable``
        If True, compiles PDF function (normalized) via Numba JIT. \n
        If callable, uses provided function directly. \n
        default: False
    create_rvs : ``bool`` or ``callable``
        If True, compiles random variable sampler via Numba JIT. \n
        If callable, uses provided function directly. \n
        default: False
    multiprocessing_function : ``bool``
        If True, assumes function accepts vectorized input for 2D evaluation. \n
        default: False
    callback : ``str`` or ``None``
        Name of the method to call when instance is called (e.g., 'function', 'pdf', 'rvs'). \n
        default: None
    use_spline_ppf : ``bool``
        If True, uses spline-based inverse-CDF sampling when creating random-variable samplers. \n
        default: False
    cdf_size : ``int``
        Number of points to use when computing CDF for inverse transform sampling. \n

    Examples
    --------
    1D cubic spline interpolation:

    >>> import numpy as np
    >>> from ler.utils import FunctionConditioning
    >>> x = np.linspace(0, 10, 100)
    >>> def f(x): return np.sin(x)
    >>> fc = FunctionConditioning(
    ...     function=f,
    ...     x_array=x,
    ...     identifier_dict={'param': 'sin_function'},
    ...     create_function=True,
    ...     callback='function'
    ... )
    >>> result = fc.function(5.0)

    2D conditioning with PDF:

    >>> y = np.linspace(0, 5, 50)
    >>> def f2d(x, y): return np.exp(-x**2/y)
    >>> fc2d = FunctionConditioning(
    ...     function=f2d,
    ...     x_array=x,
    ...     conditioned_y_array=y,
    ...     identifier_dict={'param': 'exp_function_2d'},
    ...     create_pdf=True,
    ...     create_rvs=True
    ... )

    Instance Methods
    ----------
    FunctionConditioning has the following public methods: \n
    +-------------------------------------------------+------------------------------------------------+
    | Method                                          | Description                                    |
    +=================================================+================================================+
    | :meth:`~__call__`                               | Call the callback method with arguments        |
    +-------------------------------------------------+------------------------------------------------+

    Instance Attributes
    ----------
    FunctionConditioning has the following attributes: \n
    +-----------------------------------------------+------------------+-------+------------------------------------------------+
    | Attribute                                     | Type             | Unit  | Description                                    |
    +===============================================+==================+=======+================================================+
    | :attr:`~info`                                 | ``dict``         |       | Identifier dictionary for the interpolator     |
    +-----------------------------------------------+------------------+-------+------------------------------------------------+
    | :attr:`~callback`                             | ``str``          |       | Name of method to invoke when calling instance |
    +-----------------------------------------------+------------------+-------+------------------------------------------------+
    | :attr:`~function`                             | ``callable``     |       | Interpolated function (if created)             |
    +-----------------------------------------------+------------------+-------+------------------------------------------------+
    | :attr:`~function_inverse`                     | ``callable``     |       | Inverse interpolated function (if created)     |
    +-----------------------------------------------+------------------+-------+------------------------------------------------+
    | :attr:`~pdf`                                  | ``callable``     |       | Normalized PDF function (if created)           |
    +-----------------------------------------------+------------------+-------+------------------------------------------------+
    | :attr:`~rvs`                                  | ``callable``     |       | Random variable sampler (if created)           |
    +-----------------------------------------------+------------------+-------+------------------------------------------------+
    | :attr:`~x_array`                              | ``numpy.ndarray``|       | X-axis points where function is evaluated      |
    +-----------------------------------------------+------------------+-------+------------------------------------------------+
    | :attr:`~z_array`                              | ``numpy.ndarray``|       | Function values at x_array points              |
    +-----------------------------------------------+------------------+-------+------------------------------------------------+
    | :attr:`~conditioned_y_array`                  | ``numpy.ndarray``|       | Y-axis conditioning points (2D case)           |
    +-----------------------------------------------+------------------+-------+------------------------------------------------+
    | :attr:`~function_spline`                      | ``numpy.ndarray``|       | Cubic spline coefficients for function         |
    +-----------------------------------------------+------------------+-------+------------------------------------------------+
    | :attr:`~function_inverse_spline`              | ``numpy.ndarray``|       | Cubic spline coefficients for inverse          |
    +-----------------------------------------------+------------------+-------+------------------------------------------------+
    | :attr:`~pdf_norm_const`                       | ``float`` or ``numpy.ndarray`` |   | PDF normalization constant(s)       |
    +-----------------------------------------------+------------------+-------+------------------------------------------------+
    | :attr:`~cdf_values`                           | ``numpy.ndarray``|       | CDF values for inverse transform sampling      |
    +-----------------------------------------------+------------------+-------+------------------------------------------------+
    | :attr:`~y_array`                              | ``numpy.ndarray``|       | Y values (for KDE case)                        |
    +-----------------------------------------------+------------------+-------+------------------------------------------------+
    | :attr:`~kde_object`                           | ``scipy.stats.gaussian_kde`` |   | KDE object (if gaussian_kde=True)  |
    +-----------------------------------------------+------------------+-------+------------------------------------------------+

    Notes
    -----
    - Interpolators are cached as JSON for efficiency; use create_new=True to regenerate. \n
    - For PDF-like functions, set non_negative_function=True to avoid negative values. \n
    - Numba JIT compilation requires ~0.5-1 second on first call for startup. \n
    - 2D interpolation requires arrays sorted by conditioned_y_array values internally. \n
    """

    def __init__(
        self,
        function=None,  # can also be an array of function values
        x_array=None,
        conditioned_y_array=None,  # if this is not none, 2D interpolation will be used
        y_array=None,
        non_negative_function=False,
        gaussian_kde=False,
        gaussian_kde_kwargs=None,
        identifier_dict=None,
        directory="./interpolator_json",
        sub_directory="default",
        name="default",
        create_new=False,
        create_function=False,
        create_function_inverse=False,
        create_pdf=False,
        create_rvs=False,
        multiprocessing_function=False,
        callback=None,
        use_spline_ppf=False,
        cdf_size=500,
    ):
        """
        Initialize FunctionConditioning with function interpolation setup.

        Loads or creates an interpolator based on provided function and parameters. \n
        Handles caching via JSON files and supports both cubic spline and KDE approaches.

        Parameters
        ----------
        See class docstring for detailed parameter descriptions.

        Raises
        ------
        ValueError
            If function is 2D array but conditioned_y_array is None.
        """
        self.cdf_size = cdf_size

        if gaussian_kde_kwargs is None:
            gaussian_kde_kwargs = {}
        if identifier_dict is None:
            identifier_dict = {}

        self.use_spline_ppf = use_spline_ppf

        # Determine whether to create a new interpolator or use provided callables
        create = self._create_decision_function(
            create_function, create_function_inverse, create_pdf, create_rvs
        )
        self.info = identifier_dict.copy()
        self.callback = callback

        if create:
            # create_interpolator input list
            input_list = [
                function,
                x_array,
                conditioned_y_array,
                create_function_inverse,
                create_pdf,
                create_rvs,
                multiprocessing_function,
            ]
            input_list_kde = [x_array, y_array, gaussian_kde_kwargs]

            # check first whether the directory, subdirectory and pickle exist
            path_inv_cdf, it_exist = interpolator_json_path(
                identifier_dict=identifier_dict,
                directory=directory,
                sub_directory=sub_directory,
                interpolator_name=name,
            )

            # if the interpolator exists, load it
            if create_new:
                it_exist = False

            if it_exist:
                print(f"{name} interpolator will be loaded from {path_inv_cdf}")
                # load the interpolator, where the interpolator is a dictionary
                # interpolator = {
                #     'x_array': x_array,
                #     'z_array': z_array,
                #     'conditioned_y_array': conditioned_y_array,
                #     'function_spline': function_spline,
                #     'function_inverse_spline': function_inverse_spline,
                #     'pdf_norm_const': pdf_norm_const,
                #     'cdf_values': cdf_values,
                # }
                interpolator = load_json(path_inv_cdf)
            else:
                print(f"{name} interpolator will be generated at {path_inv_cdf}")
                if not gaussian_kde:
                    interpolator = self._create_interpolator(*input_list)
                    # save the interpolator
                    save_json(path_inv_cdf, interpolator)
                else:
                    # gaussian KDE
                    interpolator = self._create_gaussian_kde(*input_list_kde)

                    # ------------
                    # check scipy gaussian kde can be saved as json
                    # ------------
                    # save_json(path_inv_cdf, interpolator)

            if not gaussian_kde:
                # Convert loaded JSON lists to NumPy arrays for Numba compatibility
                x_array = np.array(interpolator["x_array"])
                z_array = np.array(interpolator["z_array"])
                conditioned_y_array = (
                    np.array(interpolator["conditioned_y_array"])
                    if interpolator["conditioned_y_array"] is not None
                    else None
                )
                y_array = None
                function_spline = np.array(interpolator["function_spline"])
                function_inverse_spline = (
                    np.array(interpolator["function_inverse_spline"])
                    if interpolator["function_inverse_spline"] is not None
                    else None
                )
                pdf_norm_const = (
                    np.array(interpolator["pdf_norm_const"])
                    if interpolator["pdf_norm_const"] is not None
                    else None
                )
                cdf_values = (
                    np.array(interpolator["cdf_values"])
                    if interpolator["cdf_values"] is not None
                    else None
                )
                # x_array_cdf = (
                #     np.array(interpolator["x_array_cdf"])
                #     if interpolator["x_array_cdf"] is not None
                #     else None
                # )
                if "x_array_cdf" in interpolator and interpolator["x_array_cdf"] is not None:
                    x_array_cdf = np.array(interpolator["x_array_cdf"])
                else:
                    x_array_cdf = None


                if conditioned_y_array is None:
                    # function is 1D
                    # njit(lambda x: cubic_spline_interpolator(x, function_spline, x_array))
                    def function_any(x):
                        return cubic_spline_interpolator(x, function_spline, x_array)

                    def function_non_zero(x):
                        result = cubic_spline_interpolator(x, function_spline, x_array)
                        idx = result < 0.0
                        result[idx] = 0.0
                        return result

                    function_final = (
                        function_non_zero if non_negative_function else function_any
                    )

                    self.function = (
                        njit(function_final) if create_function else None
                    )

                    # inverse function is 1D
                    def function_any(x):
                        return cubic_spline_interpolator(
                            x, function_inverse_spline, z_array
                        )

                    def function_non_zero(x):
                        result = cubic_spline_interpolator(
                            x, function_inverse_spline, z_array
                        )
                        idx = result < 0.0
                        result[idx] = 0.0
                        return result

                    function_final = (
                        function_non_zero if non_negative_function else function_any
                    )

                    self.function_inverse = (
                        njit(function_final)
                        if create_function_inverse
                        else None
                    )

                    # pdf is 1D
                    @njit
                    def pdf_wrapper(x):
                        return pdf_cubic_spline_interpolator(
                            x, pdf_norm_const, function_spline, x_array
                        )

                    self.pdf = pdf_wrapper if create_pdf else None

                    if self.use_spline_ppf:
                        @njit
                        def rvs_wrapper(size):
                            return inverse_transform_sampler(
                                size=size, cdf=cdf_values, x=x_array_cdf
                            )
                    else:
                        # sampler is 1D
                        @njit
                        def rvs_wrapper(size):
                            return inverse_transform_sampler(
                                size=size, cdf=cdf_values, x=x_array_cdf
                            )

                    self.rvs = rvs_wrapper if create_rvs else None

                else:
                    # function is 2D
                    def function_any(x, y):
                        return cubic_spline_interpolator2d(
                            xnew_array=x, ynew_array=y, coefficients=function_spline, x=x_array, y=conditioned_y_array
                        )

                    def function_non_zero(x, y):
                        result = cubic_spline_interpolator2d(
                            xnew_array=x, ynew_array=y, coefficients=function_spline, x=x_array, y=conditioned_y_array
                        )
                        idx = result < 0.0
                        result[idx] = 0.0
                        return result

                    function_final = (
                        function_non_zero if non_negative_function else function_any
                    )
                    self.function = (
                        njit(function_final) if create_function else None
                    )

                    # inverse function is 2D
                    def function_any(x, y):
                        return cubic_spline_interpolator2d(
                            xnew_array=x, ynew_array=y, coefficients=function_inverse_spline, x=z_array, y=conditioned_y_array
                        )

                    def function_non_zero(x, y):
                        result = cubic_spline_interpolator2d(
                            xnew_array=x, ynew_array=y, coefficients=function_inverse_spline, x=z_array, y=conditioned_y_array
                        )
                        idx = result < 0.0
                        result[idx] = 0.0
                        return result

                    function_final = (
                        function_non_zero if non_negative_function else function_any
                    )
                    self.function_inverse = (
                        njit(function_final)
                        if create_function_inverse
                        else None
                    )

                    @njit
                    def pdf_wrapper(x, y):
                        return pdf_cubic_spline_interpolator2d(
                            xnew_array=x,
                            ynew_array=y,
                            norm_array=pdf_norm_const,
                            coefficients=function_spline,
                            x=x_array,
                            y=conditioned_y_array,
                        )

                    self.pdf = pdf_wrapper if create_pdf else None

                    if self.use_spline_ppf:
                        @njit
                        def rvs_wrapper(size, y):
                            return inverse_transform_sampler2d_spline(
                                size=size,
                                conditioned_y=y,
                                cdf2d=cdf_values,
                                x2d=x_array_cdf,
                                y1d=conditioned_y_array,
                            )
                    else:
                        @njit
                        def rvs_wrapper(size, y):
                            return inverse_transform_sampler2d(
                                size=size,
                                conditioned_y=y,
                                cdf2d=cdf_values,
                                x2d=x_array_cdf,
                                y1d=conditioned_y_array,
                            )

                    self.rvs = rvs_wrapper if create_rvs else None

                self.x_array = x_array
                self.z_array = z_array
                self.conditioned_y_array = conditioned_y_array
                self.function_spline = function_spline
                self.function_inverse_spline = function_inverse_spline
                self.pdf_norm_const = pdf_norm_const
                self.cdf_values = cdf_values
                self.x_array_cdf = x_array_cdf
            else:
                x_array = interpolator["x_array"]
                y_array = interpolator["y_array"]
                kde_object = interpolator["kde_object"]

                def kde_pdf(x):
                    return kde_object.pdf(x)

                self.pdf = kde_pdf if create_pdf else None
                if y_array is None:

                    def kde_rvs(size):
                        return kde_object.resample(size)[0]

                    self.rvs = kde_rvs if create_rvs else None
                else:

                    def kde_rvs(size):
                        return kde_object.resample(size)

                    self.rvs = kde_rvs if create_rvs else None

                self.x_array = x_array
                self.y_array = y_array
                self.kde_object = kde_object
        else:
            self.conditioned_y_array = None
            self.x_array = None
            self.y_array = None
            self.kde_object = None

    def __call__(self, *args):
        """
        Call the instance with arguments forwarded to the callback method.

        Converts float arguments to reshaped arrays for compatibility with \n
        interpolation functions and forwards to the method specified by callback.

        Parameters
        ----------
        *args : ``float`` or ``numpy.ndarray``
            Arguments to pass to callback method. Float values are converted to \n
            reshaped (-1,) arrays for interpolation compatibility.

        Returns
        -------
        result : ``float`` or ``numpy.ndarray``
            Output from the callback method.
        """
        args = [
            np.array(arg).reshape(-1) if isinstance(arg, float) else arg for arg in args
        ]
        return getattr(self, self.callback)(*args)

    def _create_decision_function(
        self, create_function, create_function_inverse, create_pdf, create_rvs
    ):
        """
        Assign user-provided functions or determine if interpolation is needed.

        Checks if any of the create_* parameters are callable. If so, assigns them \n
        directly as attributes and returns False to skip interpolator creation. \n
        Otherwise returns True to proceed with interpolator generation.

        Parameters
        ----------
        create_function : ``bool`` or ``callable``
            Whether to create function or callable to use directly.
        create_function_inverse : ``bool`` or ``callable``
            Whether to create inverse function or callable to use directly.
        create_pdf : ``bool`` or ``callable``
            Whether to create PDF or callable to use directly.
        create_rvs : ``bool`` or ``callable``
            Whether to create sampler or callable to use directly.

        Returns
        -------
        decision_bool : ``bool``
            True if interpolator creation is needed, False if user provided callables.
        """
        decision_bool = True
        if not isinstance(create_function, bool) and callable(create_function):
            self.function = create_function
            decision_bool = False
        if not isinstance(create_function_inverse, bool) and callable(
            create_function_inverse
        ):
            self.function_inverse = create_function_inverse
            decision_bool = False
        if not isinstance(create_pdf, bool) and callable(create_pdf):
            self.pdf = create_pdf
            decision_bool = False
        if not isinstance(create_rvs, bool) and callable(create_rvs):
            self.rvs = create_rvs
            decision_bool = False

        return decision_bool

    def _create_gaussian_kde(self, x_array, y_array, gaussian_kde_kwargs):
        """
        Create a Gaussian Kernel Density Estimator interpolator.

        Constructs a KDE object for 1D or 2D data using scipy.stats.gaussian_kde.

        Parameters
        ----------
        x_array : ``numpy.ndarray``
            X-axis data points (1D array).
        y_array : ``numpy.ndarray`` or ``None``
            Y-axis data points (1D array). If None, performs 1D KDE.
        gaussian_kde_kwargs : ``dict``
            Keyword arguments for gaussian_kde constructor (e.g., bw_method).

        Returns
        -------
        interpolator : ``dict``
            Dictionary containing 'x_array', 'y_array', and 'kde_object' keys.
        """
        # 1d KDE
        if y_array is None:
            kde = gaussian_kde(x_array, **gaussian_kde_kwargs)
        else:
            data = np.vstack([x_array, y_array])
            kde = gaussian_kde(data, **gaussian_kde_kwargs)

        return {
            "x_array": x_array,
            "y_array": y_array,
            "kde_object": kde,
        }

    def _create_interpolator(
        self,
        function,
        x_array,
        conditioned_y_array,
        create_function_inverse,
        create_pdf,
        create_rvs,
        multiprocessing_function,
    ):
        """
        Create cubic spline interpolator data structure.

        Evaluates function on provided grid, computes spline coefficients, and \n
        optionally generates inverse spline, PDF normalization, and CDF values for \n
        inverse transform sampling.

        Parameters
        ----------
        function : ``callable`` or ``numpy.ndarray``
            Function to interpolate or array of function values.
        x_array : ``numpy.ndarray``
            1D or 2D array of x-axis points.
        conditioned_y_array : ``numpy.ndarray`` or ``None``
            1D array of conditioning y values (1D interpolation if None).
        create_function_inverse : ``bool``
            Whether to generate inverse function spline.
        create_pdf : ``bool``
            Whether to compute PDF normalization and CDF for sampling.
        create_rvs : ``bool``
            Whether to generate CDF values for inverse transform sampling.
        multiprocessing_function : ``bool``
            Whether function accepts vectorized 2D evaluation.

        Returns
        -------
        interpolator : ``dict``
            Dictionary with keys: 'x_array', 'z_array', 'conditioned_y_array', \n
            'function_spline', 'function_inverse_spline', 'pdf_norm_const', 'cdf_values'.
        """
        # function can be numpy array or callable
        # x_array, z_array are 2D arrays if conditioned_y_array is not None
        x_array, z_array, conditioned_y_array = self._create_z_array(
            x_array,
            function,
            conditioned_y_array,
            create_pdf,
            create_rvs,
            multiprocessing_function,
        )
        del function

        function_spline = self._function_spline_generator(
            x_array, z_array, conditioned_y_array
        )

        if create_function_inverse:
            if conditioned_y_array is None:
                idx_sort = np.argsort(z_array)
            else:
                idx_sort = np.argsort(z_array, axis=1)
            x_array_ = x_array[idx_sort]
            z_array_ = z_array[idx_sort]

            # check z_array is strictly increasing
            # if (not np.all(np.diff(z_array) > 0)) or (not np.all(np.diff(z_array) < 0)):
            #     raise ValueError("z_array must be strictly increasing")

            function_inverse_spline = self._function_spline_generator(
                z_array_, x_array_, conditioned_y_array
            )
        else:
            function_inverse_spline = None

        if create_pdf or create_rvs:
            # cannot have -ve pdf
            pdf_norm_const = self._pdf_norm_const_generator(
                x_array, function_spline, conditioned_y_array
            )

            if create_rvs:
                cdf_values, x_array_cdf = self._cdf_values_generator(
                    x_array, z_array, conditioned_y_array, function_spline
                )
            else:
                cdf_values = None
                x_array_cdf = None
        else:
            pdf_norm_const = None
            cdf_values = None
            x_array_cdf = None

        return {
            "x_array": x_array,
            "z_array": z_array,
            "conditioned_y_array": conditioned_y_array,
            "function_spline": function_spline,
            "function_inverse_spline": function_inverse_spline,
            "pdf_norm_const": pdf_norm_const,
            "cdf_values": cdf_values,
            "x_array_cdf": x_array_cdf,
        }

    def _create_z_array(
        self,
        x_array,
        function,
        conditioned_y_array,
        create_pdf,
        create_rvs,
        multiprocessing_function,
    ):
        """
        Evaluate function on grid and remove NaN values.

        Handles both callable and array inputs, processes 1D and 2D cases, and \n
        clamps negative values to zero if creating PDF or sampler.

        Parameters
        ----------
        x_array : ``numpy.ndarray``
            1D or 2D array of x-axis evaluation points.
        function : ``callable`` or ``numpy.ndarray``
            Function to evaluate or pre-computed values.
        conditioned_y_array : ``numpy.ndarray`` or ``None``
            Conditioning parameter; None for 1D case.
        create_pdf : ``bool``
            Set floor at 0 if True (PDF non-negative constraint).
        create_rvs : ``bool``
            Set floor at 0 if True (sampler non-negative constraint).
        multiprocessing_function : ``bool``
            Whether function supports vectorized 2D evaluation.

        Returns
        -------
        x_array : ``numpy.ndarray``
            Cleaned and sorted x-axis points.
        z_array : ``numpy.ndarray``
            Function values at x_array, cleaned and normalized.
        conditioned_y_array : ``numpy.ndarray`` or ``None``
            Sorted conditioning parameter if provided.
        """
        if callable(function):
            # 1D
            if conditioned_y_array is None:
                z_array = function(x_array)
                # remove nan values
                idx = np.argwhere(np.isnan(z_array))
                x_array = np.delete(x_array, idx)
                z_array = np.delete(z_array, idx)
            # 2D
            else:
                # check if x_array is 2D, if not, make it 2D of shape (len(conditioned_y_array), len(x_array))
                if x_array.ndim == 1:
                    x_array = np.array([x_array] * len(conditioned_y_array))

                idx = np.argsort(conditioned_y_array)
                conditioned_y_array = conditioned_y_array[idx]
                # x_array is 2D here, each row corresponds to a different y value
                x_array = x_array[idx]
                # sort each row of x_array
                x_array = np.sort(x_array, axis=1)

                if multiprocessing_function:
                    z_array = function(x_array, conditioned_y_array)
                else:
                    z_list = []
                    for i, y in enumerate(conditioned_y_array):
                        try:
                            z_list.append(
                                function(x_array[i], y * np.ones_like(x_array[i]))
                            )
                        except Exception:
                            # print(x_array[i], y)
                            z_list.append(function(x_array[i], y))
                    z_array = np.array(z_list)
        else:
            if conditioned_y_array is None:
                z_array = function
                # remove nan values
                idx = np.argwhere(np.isnan(z_array))
                x_array = np.delete(x_array, idx)
                z_array = np.delete(z_array, idx)
            else:
                if x_array.ndim == 1:
                    x_array = np.array([x_array] * len(conditioned_y_array))
                if function.ndim == 1:
                    raise ValueError(
                        "function must be 2D array if conditioned_y_array is not None"
                    )
                # row sort
                idx = np.argsort(conditioned_y_array)
                conditioned_y_array = conditioned_y_array[idx]
                # x_array is 2D here, each row corresponds to a different y value
                x_array = x_array[idx]
                z_array = function[idx]

                z_list = []
                x_list = []
                for i in range(len(conditioned_y_array)):
                    # column sort
                    idx = np.argsort(x_array[i])
                    x_list.append(x_array[i][idx])
                    z_list.append(z_array[i][idx])
                x_array = np.array(x_list)
                z_array = np.array(z_list)

        # cannot have -ve pdf
        if create_pdf or create_rvs:
            z_array[z_array < 0.0] = 0.0

        return x_array, z_array, conditioned_y_array

    def _upsample_x_array(self, x_array, size):
        """
        Increase the number of x points without replacing the original spacing profile.

        Parameters
        ----------
        x_array : ``numpy.ndarray``
            Monotonic 1D x-axis points.
        size : ``int``
            Target number of points.

        Returns
        -------
        x_array_new : ``numpy.ndarray``
            Upsampled x-axis points.
        """
        if len(x_array) >= size:
            return x_array

        index_array = np.arange(len(x_array), dtype=float)
        index_array_new = np.linspace(0.0, len(x_array) - 1, size)
        return np.interp(index_array_new, index_array, x_array)

    def _cdf_values_generator(self, x_array, z_array, conditioned_y_array, function_spline):
        """
        Compute normalized cumulative distribution function from PDF values.

        Integrates z_array (PDF) along x_array to produce normalized CDF for \n
        use in inverse transform sampling. Handles both 1D and 2D cases.

        Parameters
        ----------
        x_array : ``numpy.ndarray``
            1D or 2D array of x-axis points.
        z_array : ``numpy.ndarray``
            1D or 2D array of PDF values at x_array points.
        conditioned_y_array : ``numpy.ndarray`` or ``None``
            Conditioning parameter; None for 1D case.

        Returns
        -------
        cdf_values : ``numpy.ndarray``
            Normalized cumulative integration of z_array values.
        """
        # 1D case
        cdf_size = self.cdf_size
        if conditioned_y_array is None:
            z_array[z_array<0.] = 0. # already done
            if self.use_spline_ppf:
                cdf_values, x_array_cdf, _ = cumulative_spline(y=z_array, x=x_array, initial=0.0)
            else:
                # create more sampled x_array for better integration accuracy
                size = np.max([len(x_array), cdf_size])
                # use the spline to evaluate z_array at the new x_array
                x_array_new = self._upsample_x_array(x_array, size)
                z_array_new = cubic_spline_interpolator(x_array_new, function_spline, x_array)
                # cdf_values are already normalized internally
                cdf_values, x_array_cdf, _ = cumulative_trapezoid(y=z_array_new, x=x_array_new, initial=0.0)
        # 2D case
        else:
            cdf_values = []
            x_array_cdf = []
            for i, y in enumerate(conditioned_y_array):
                z_array_ = z_array[i]
                if self.use_spline_ppf:
                    z_array_[z_array_ < 0.0] = 0.0
                    cdfs_, x_array_cdf_, _ = cumulative_spline(y=z_array_, x=x_array[i], initial=0.0)
                else:
                    # create more sampled x_array for better integration accuracy
                    size = np.max([len(x_array[i]), cdf_size])
                    x_array_new = self._upsample_x_array(x_array[i], size)
                    z_array_new = cubic_spline_interpolator(x_array_new, function_spline[i], x_array[i])
                    # set negative values to 0
                    z_array_new[z_array_new < 0.0] = 0.0

                    cdfs_, x_array_cdf_, _ = cumulative_trapezoid(y=z_array_new, x=x_array_new, initial=0.0)

                cdf_values.append(cdfs_)
                x_array_cdf.append(x_array_cdf_)

        return np.array(cdf_values), np.array(x_array_cdf)

    def _pdf_norm_const_generator(self, x_array, function_spline, conditioned_y_array):
        """
        Compute PDF normalization constant via numerical integration.

        Integrates the spline-interpolated function over the domain to determine \n
        the normalization constant needed to make the PDF integrate to 1. \n
        Handles both 1D and 2D cases.

        Parameters
        ----------
        x_array : ``numpy.ndarray``
            1D or 2D array of x-axis points.
        function_spline : ``numpy.ndarray``
            Cubic spline coefficients from CubicSpline.
        conditioned_y_array : ``numpy.ndarray`` or ``None``
            Conditioning parameter; None for 1D case.

        Returns
        -------
        norm : ``float`` or ``numpy.ndarray``
            Normalization constant(s) for the PDF.
        """
        cdf_size = self.cdf_size
        # 1D case
        if conditioned_y_array is None:
            # Use a dense grid (like in CDF generation) so trapezoidal integration
            # is accurate even when the original x_array resolution is low.
            size = np.max([len(x_array), cdf_size])
            x_array_new = self._upsample_x_array(x_array, size)
            z_array_new = cubic_spline_interpolator(x_array_new, function_spline, x_array)
            z_array_new[z_array_new < 0.0] = 0.0
            # `np.trapezoid` is preferred, but older NumPy may only provide `np.trapz`.
            integrate = getattr(np, "trapezoid", np.trapz)
            return float(integrate(z_array_new, x_array_new))
        # 2D case
        else:
            norm = []
            for i, y in enumerate(conditioned_y_array):
                size = np.max([len(x_array[i]), cdf_size])
                x_array_new = self._upsample_x_array(x_array[i], size)
                z_array_new = cubic_spline_interpolator(x_array_new, function_spline[i], x_array[i])
                z_array_new[z_array_new < 0.0] = 0.0
                integrate = getattr(np, "trapezoid", np.trapz)
                norm.append(float(integrate(z_array_new, x_array_new)))

            return np.array(norm)

    def _function_spline_generator(self, x_array, z_array, conditioned_y_array):
        """
        Generate cubic spline coefficients from data points.

        Creates scipy CubicSpline object coefficients for smooth interpolation. \n
        Returns spline coefficients in format compatible with Numba JIT compilation. \n
        Handles both 1D and 2D cases independently.

        Parameters
        ----------
        x_array : ``numpy.ndarray``
            1D or 2D array of x-axis evaluation points.
        z_array : ``numpy.ndarray``
            1D or 2D array of function values.
        conditioned_y_array : ``numpy.ndarray`` or ``None``
            Conditioning parameter; None for 1D case.

        Returns
        -------
        function_spline : ``numpy.ndarray``
            Cubic spline polynomial coefficients (either 2D for 1D case or \n
            3D array of 2D arrays for 2D case).
        """
        # 1D case
        if conditioned_y_array is None:
            function_spline = CubicSpline(x_array, z_array, extrapolate=True).c
        # 2D case
        else:
            function_spline = []
            for i, y in enumerate(conditioned_y_array):
                function_spline.append(CubicSpline(x_array[i], z_array[i], extrapolate=True).c)

        return np.array(function_spline)
