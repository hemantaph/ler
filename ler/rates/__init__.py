# -*- coding: utf-8 -*-
"""
Rate drivers: ``LeR`` (lensed + unlensed CBC) and ``GWRATES`` (unlensed CBC).

Public symbols are resolved lazily so importing ``ler.rates`` does not initialize
the full lensed/rate stack until ``LeR`` or ``GWRATES`` is requested.
"""

from importlib import import_module

_SUBMODULES = {
    "ler": ".ler",
    "gwrates": ".gwrates",
}

_SYMBOL_TO_MODULE = {
    "LeR": ".ler",
    "GWRATES": ".gwrates",
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
