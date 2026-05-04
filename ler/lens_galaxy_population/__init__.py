# -*- coding: utf-8 -*-
"""
Lens galaxy population: parameter distributions, optical depth, sampling, cross-section grids.

Public symbols are resolved lazily so importing this package does not load every
lensing sampler and numerical backend up front.
"""

from importlib import import_module

_SUBMODULES = {
    "lens_galaxy_parameter_distribution": ".lens_galaxy_parameter_distribution",
    "optical_depth": ".optical_depth",
    "mp": ".mp",
    "sampler_functions": ".sampler_functions",
    "lens_functions": ".lens_functions",
    "cross_section_interpolator": ".cross_section_interpolator",
}

_SYMBOL_TO_MODULE = {
    "LensGalaxyParameterDistribution": ".lens_galaxy_parameter_distribution",
    "OpticalDepth": ".optical_depth",
    "lens_redshift_strongly_lensed_njit": ".mp",
    "lens_redshift_strongly_lensed_mp": ".mp",
    "cross_section_unit_mp": ".mp",
    "cross_section_mp": ".mp",
    "lens_redshift_strongly_lensed_sis_analytical_pdf": ".sampler_functions",
    "lens_redshift_strongly_lensed_sis_analytical_rvs": ".sampler_functions",
    "velocity_dispersion_ewoud_denisty_function": ".sampler_functions",
    "velocity_dispersion_bernardi_denisty_function": ".sampler_functions",
    "gengamma_function": ".sampler_functions",
    "gengamma_pdf": ".sampler_functions",
    "gengamma_rvs": ".sampler_functions",
    "rayleigh_rvs": ".sampler_functions",
    "rayleigh_pdf": ".sampler_functions",
    "axis_ratio_padilla_strauss_rvs": ".sampler_functions",
    "axis_ratio_padilla_strauss_pdf": ".sampler_functions",
    "rejection_sampler_full": ".sampler_functions",
    "rejection_sampler_partial": ".sampler_functions",
    "importance_sampler_partial": ".sampler_functions",
    "importance_sampler_full": ".sampler_functions",
    "create_sampler": ".sampler_functions",
    "phi_cut_SIE": ".lens_functions",
    "einstein_radius": ".lens_functions",
    "cross_section": ".lens_functions",
    "make_cross_section_area_reinit": ".cross_section_interpolator",
}

__all__ = sorted([*_SUBMODULES, *_SYMBOL_TO_MODULE, "available_sampler_list"])


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


def available_sampler_list():
    """
    Return string keys for lens-parameter samplers and related helpers.

    Each entry is the name of a callable (PDF, inverse CDF draw, or helper) in
    this subpackage; pass the string where API expects a ``sampler`` or
    ``_function`` key (see ``LensGalaxyParameterDistribution`` / ``OpticalDepth``).

    Returns
    -------
    list of str
        Registered sampler-related names.

    Examples
    --------
    >>> from ler.lens_galaxy_population import available_sampler_list
    >>> names = available_sampler_list()
    >>> "importance_sampler" in names
    True
    """
    return [
        "lens_redshift_strongly_lensed_sis_analytical_pdf",
        "lens_redshift_strongly_lensed_sis_analytical_rvs",
        "velocity_dispersion_ewoud_denisty_function",
        "velocity_dispersion_bernardi_denisty_function",
        "gengamma_function",
        "gengamma_pdf",
        "gengamma_rvs",
        "rayleigh_rvs",
        "rayleigh_pdf",
        "axis_ratio_padilla_strauss_rvs",
        "axis_ratio_padilla_strauss_pdf",
        "bounded_normal_sample",
        "rejection_sampler",
        "importance_sampler",
        "importance_sampler_mp",
    ]
