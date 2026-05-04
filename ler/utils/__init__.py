# -*- coding: utf-8 -*-
"""
Shared numerics, plotting, I/O, splines, and cosmology helpers.

Public symbols are resolved lazily so importing one helper does not also import
plotting, ``gwsnr``, or the full interpolation/cosmology stack.
"""

from importlib import import_module

_SUBMODULES = {
    "utils": ".utils",
    "plots": ".plots",
    "gwsnr_training_data_generator": ".gwsnr_training_data_generator",
    "function_interpolation": ".function_interpolation",
    "cosmological_conversions": ".cosmological_conversions",
}

_SYMBOL_TO_MODULE = {
    "is_njitted": ".utils",
    "remove_file": ".utils",
    "load_pickle": ".utils",
    "save_pickle": ".utils",
    "load_hdf5": ".utils",
    "save_hdf5": ".utils",
    "NumpyEncoder": ".utils",
    "load_json": ".utils",
    "save_json": ".utils",
    "append_json": ".utils",
    "get_param_from_json": ".utils",
    "add_dictionaries_together": ".utils",
    "trim_dictionary": ".utils",
    "load_txt_from_module": ".utils",
    "rejection_sample": ".utils",
    "rejection_sample2d": ".utils",
    "create_func": ".utils",
    "create_func_inv": ".utils",
    "create_pdf": ".utils",
    "create_inv_cdf_array": ".utils",
    "create_conditioned_pdf": ".utils",
    "create_conditioned_inv_cdf_array": ".utils",
    "interpolator_from_json": ".utils",
    "interpolator_json_path": ".utils",
    "monte_carlo_integration": ".utils",
    "cubic_spline_interpolator": ".utils",
    "cubic_spline_interpolator_scalar": ".utils",
    "pdf_cubic_spline_interpolator": ".utils",
    "cubic_hermite_interpolation": ".utils",
    "cubic_spline_interpolator2d": ".utils",
    "pdf_cubic_spline_interpolator2d": ".utils",
    "inverse_transform_sampler2d": ".utils",
    "inverse_transform_sampler2d_spline": ".utils",
    "inverse_transform_sampler": ".utils",
    "inverse_transform_sampler_spline": ".utils",
    "trim_cdf": ".utils",
    "normal_pdf": ".utils",
    "normal_pdf_2d": ".utils",
    "cumulative_spline": ".utils",
    "cumulative_trapezoid": ".utils",
    "get_pchip_spline_coeffs": ".utils",
    "batch_handler": ".utils",
    "create_batch_params": ".utils",
    "KStest": ".utils",
    "param_plot": ".plots",
    "relative_mu_dt_unlensed": ".plots",
    "relative_mu_dt_lensed": ".plots",
    "mu_vs_dt_plot": ".plots",
    "TrainingDataGenerator": ".gwsnr_training_data_generator",
    "FunctionConditioning": ".function_interpolation",
    "generate_mixed_grid": ".cosmological_conversions",
    "luminosity_distance": ".cosmological_conversions",
    "differential_comoving_volume": ".cosmological_conversions",
    "comoving_distance": ".cosmological_conversions",
    "angular_diameter_distance": ".cosmological_conversions",
    "angular_diameter_distance_z1z2": ".cosmological_conversions",
}

__all__ = sorted([*_SUBMODULES, *_SYMBOL_TO_MODULE])


def __getattr__(name):
    if name in _SUBMODULES:
        module = import_module(_SUBMODULES[name], __name__)
        globals()[name] = module
        return module

    module_name = _SYMBOL_TO_MODULE.get(name)
    if module_name is None:
        raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

    value = getattr(import_module(module_name, __name__), name)
    globals()[name] = value
    return value


def __dir__():
    return sorted([*globals(), *__all__])
