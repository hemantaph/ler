# -*- coding: utf-8 -*-
"""
Compact-binary source population: masses, redshift, merger-rate priors, SFR models.

Public symbols are resolved lazily so importing this package does not load every
population model and dependency up front.
"""

from importlib import import_module

_SUBMODULES = {
    "cbc_source_redshift_distribution": ".cbc_source_redshift_distribution",
    "cbc_source_parameter_distribution": ".cbc_source_parameter_distribution",
    "prior_functions": ".prior_functions",
    "sfr_with_time_delay": ".sfr_with_time_delay",
    "broken_powerlaw_plus_2peaks": ".broken_powerlaw_plus_2peaks",
    "powerlaw_plus_peak": ".powerlaw_plus_peak",
}

_SYMBOL_TO_MODULE = {
    "CBCSourceRedshiftDistribution": ".cbc_source_redshift_distribution",
    "CBCSourceParameterDistribution": ".cbc_source_parameter_distribution",
    "merger_rate_density_bbh_oguri2018_function": ".prior_functions",
    "merger_rate_density_bbh_popIII_ken2022_function": ".prior_functions",
    "merger_rate_density_madau_dickinson2014_function": ".prior_functions",
    "merger_rate_density_madau_dickinson_belczynski_ng_function": ".prior_functions",
    "merger_rate_density_bbh_primordial_ken2022_function": ".prior_functions",
    "sfr_madau_fragos2017_with_bbh_td": ".prior_functions",
    "sfr_madau_dickinson2014_with_bbh_td": ".prior_functions",
    "sfr_madau_fragos2017_with_bns_td": ".prior_functions",
    "sfr_madau_dickinson2014_with_bns_td": ".prior_functions",
    "sfr_madau_fragos2017": ".prior_functions",
    "sfr_madau_dickinson2014": ".prior_functions",
    "ng2022_lognormal_joint_pdf": ".prior_functions",
    "binary_masses_BBH_popIII_lognormal_rvs": ".prior_functions",
    "binary_masses_BBH_primordial_lognormal_rvs": ".prior_functions",
    "bimodal_pdf": ".prior_functions",
    "bimodal_cdf": ".prior_functions",
    "binary_masses_BNS_bimodal_rvs": ".prior_functions",
    "broken_powerlaw_pdf": ".prior_functions",
    "gaussian_plus_isotropic_pdf": ".prior_functions",
    "gaussian_plus_isotropic_joint_pdf": ".prior_functions",
    "powerlaw_pdf": ".prior_functions",
    "powerlaw_rvs": ".prior_functions",
    "truncated_normal_pdf": ".prior_functions",
    "truncated_normal_rvs": ".prior_functions",
    "powerlaw_with_smoothing": ".prior_functions",
    "create_gw_parameters_sampler": ".prior_functions",
    "sfr_with_time_delay_function": ".sfr_with_time_delay",
    "broken_powerlaw_plus_2peaks_pdf": ".broken_powerlaw_plus_2peaks",
    "broken_powerlaw_plus_2peaks_function": ".broken_powerlaw_plus_2peaks",
    "broken_powerlaw_plus_2peaks_cdf": ".broken_powerlaw_plus_2peaks",
    "broken_powerlaw_plus_2peaks_rvs": ".broken_powerlaw_plus_2peaks",
    "mass_ratio_powerlaw_with_smoothing_scalar": ".broken_powerlaw_plus_2peaks",
    "mass_ratio_powerlaw_with_smoothing_pdf": ".broken_powerlaw_plus_2peaks",
    "mass_ratio_powerlaw_with_smoothing_cdf": ".broken_powerlaw_plus_2peaks",
    "mass_ratio_powerlaw_with_smoothing_rvs": ".broken_powerlaw_plus_2peaks",
    "powerlaw_plus_peak_unnormalized": ".powerlaw_plus_peak",
    "powerlaw_plus_peak_pdf": ".powerlaw_plus_peak",
    "powerlaw_plus_peak_cdf": ".powerlaw_plus_peak",
    "powerlaw_plus_peak_function": ".powerlaw_plus_peak",
    "powerlaw_plus_peak_rvs": ".powerlaw_plus_peak",
}

__all__ = sorted([*_SUBMODULES, *_SYMBOL_TO_MODULE, "available_prior_list"])


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


def available_prior_list():
    """
    Return string keys for merger-rate models, mass priors, and related callables.

    Each entry names a function or closure defined in this subpackage; use the
    string as a ``source_priors`` (or related) key in ``CBCSourceParameterDistribution``
    / rate drivers. Mathematical forms live in the underlying module docstrings,
    not in this registry.

    Returns
    -------
    list of str
        Registered prior-related names.

    Examples
    --------
    >>> from ler.gw_source_population import available_prior_list
    >>> names = available_prior_list()
    >>> "powerlaw_rvs" in names
    True
    """
    return [
        'merger_rate_density_bbh_oguri2018_function',
        'merger_rate_density_bbh_popIII_ken2022_function',
        'merger_rate_density_madau_dickinson2014_function',
        'merger_rate_density_madau_dickinson_belczynski_ng_function',
        'merger_rate_density_bbh_primordial_ken2022_function',
        'sfr_madau_fragos2017_with_bbh_td',
        'sfr_madau_dickinson2014_with_bbh_td',
        'sfr_madau_fragos2017_with_bns_td',
        'sfr_madau_dickinson2014_with_bns_td',
        'sfr_madau_fragos2017',
        'sfr_madau_dickinson2014',
        'ng2022_lognormal_joint_pdf',
        'binary_masses_BBH_popIII_lognormal_rvs',
        'binary_masses_BBH_primordial_lognormal_rvs',
        'bimodal_pdf',
        'binary_masses_BNS_bimodal_rvs',
        'broken_powerlaw_pdf',
        'gaussian_plus_isotropic_pdf',
        'gaussian_plus_isotropic_joint_pdf',
        'powerlaw_pdf',
        'powerlaw_rvs',
        'truncated_normal_pdf',
        'truncated_normal_rvs',
        'powerlaw_with_smoothing',
        'powerlaw_plus_peak_pdf',
        'powerlaw_plus_peak_function',
        'powerlaw_plus_peak_rvs',
        'broken_powerlaw_plus_2peaks_pdf',
        'broken_powerlaw_plus_2peaks_function',
        'broken_powerlaw_plus_2peaks_rvs',
        'mass_ratio_powerlaw_with_smoothing_pdf',
        'mass_ratio_powerlaw_with_smoothing_rvs',
    ]
