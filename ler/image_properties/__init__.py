# -*- coding: utf-8 -*-
"""
Lensing image geometry: magnifications, time delays, EPL+shear numerics, caustics.

Public symbols are resolved lazily so importing this package does not load every
Numba lensing kernel up front.
"""

from importlib import import_module

_SUBMODULES = {
    "image_properties": ".image_properties",
    "multiprocessing_routine_epl_shear": ".multiprocessing_routine_epl_shear",
    "epl_shear_njit": ".epl_shear_njit",
    "cross_section_njit": ".cross_section_njit",
    "sample_caustic_points_njit": ".sample_caustic_points_njit",
}

_SYMBOL_TO_MODULE = {
    "ImageProperties": ".image_properties",
    "solve_lens_equation": ".multiprocessing_routine_epl_shear",
    "pol_to_ell": ".epl_shear_njit",
    "omega_scalar": ".epl_shear_njit",
    "lensing_diagnostics_scalar": ".epl_shear_njit",
    "fermat_potential_scalar": ".epl_shear_njit",
    "image_position_analytical_njit": ".epl_shear_njit",
    "create_epl_shear_solver": ".epl_shear_njit",
    "phi_q2_ellipticity": ".cross_section_njit",
    "ellipticity2phi_q": ".cross_section_njit",
    "omega": ".cross_section_njit",
    "cdot": ".cross_section_njit",
    "pol_to_cart": ".cross_section_njit",
    "cart_to_pol": ".cross_section_njit",
    "caustic_points_epl_shear": ".cross_section_njit",
    "caustic_area_epl_shear": ".cross_section_njit",
    "half_symmetric_polygon_area": ".cross_section_njit",
    "polygon_area": ".cross_section_njit",
    "make_cross_section_area_reinit": ".cross_section_njit",
    "cross_section_epl_shear_unit": ".cross_section_njit",
    "caustic_points_epl_shear_half": ".sample_caustic_points_njit",
    "sample_source_from_double_caustic": ".sample_caustic_points_njit",
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
