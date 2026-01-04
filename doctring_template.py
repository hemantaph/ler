"""
<one-line summary of the module.>

<longer module summary. Explain the purpose, scope, and key concepts. If there are
multiple public entry points, list them in prose (not as a table).

This module is intended to be rendered by Sphinx + AutoAPI. It uses reST-friendly
docstrings with NumPy-style sections and Sphinx roles for cross-references.

Copyright (c) 2026 Phurailatpam Hemantakumar
License: MIT
"""

from __future__ import annotations

from typing import Any, Dict, Optional

# Public imports should be explicit to help AutoAPI and readers.
# Avoid "import *" and avoid importing heavy deps at import time if possible.
from parent_class import ParentClass


class ExampleClass(ParentClass):
    r"""
    <one-line summary of the class.>

    <longer class summary. Mention what it represents, how it is typically used,
    and any important invariants.>

    Parameters
    ----------
    param1 : type
        Parameter description.
    param2 : dict
        Parameter description.

        The dictionary may contain the following keys (optional):

        +-----------+--------+-----------------------------+
        | Key       | Type   | Meaning                     |
        +===========+========+=============================+
        | ``key1``  | type   | description                 |
        +-----------+--------+-----------------------------+
        | ``key2``  | type   | description                 |
        +-----------+--------+-----------------------------+

    Attributes
    ----------
    param1 : type
        Attribute description.
    param2 : dict
        Attribute description.

    See Also
    --------
    :class:`parent_class.ParentClass`
        Brief explanation of the relationship/inheritance.

    Notes
    -----
    - Keep cross-references in Sphinx roles so AutoAPI can link them:
      :class:`~package.module.Class`, :meth:`~package.module.Class.method`,
      :attr:`~package.module.Class.attribute`, :func:`~package.module.func`.
    - Prefer listing *public* attributes/methods in prose or in the dedicated
      sections above. Avoid “Instance Methods/Attributes” tables because AutoAPI
      already generates member listings and tables can break easily in reST.

    Examples
    --------
    Basic usage:

    >>> from ler.docstring import ExampleClass
    >>> obj = ExampleClass(param1=1, param2={"key1": 2})
    >>> out = obj.method1()
    """
    # --- Class-level attributes (constants / defaults) ---
    # For “real” instance attributes, document them in class docstring under
    # Attributes and set them in __init__ (preferred).
    DEFAULT_PARAM2: Dict[str, Any] = {}
    """dict : Default value used when ``param2`` is not provided."""

    def __init__(self, param1: Any, param2: Optional[Dict[str, Any]] = None) -> None:
        """
        Construct an :class:`~ler.docstring.ExampleClass`.

        Parameters
        ----------
        param1 : type
            Parameter description.
        param2 : dict, optional
            Parameter description. If ``None``, a default is used.

        Notes
        -----
        Keep ``__init__`` concise. Longer conceptual docs belong in the class
        docstring, not here.
        """
        self.param1 = param1
        self.param2 = self.DEFAULT_PARAM2 if param2 is None else param2

    def _method1_helper(self) -> None:
        """
        Internal helper for :meth:`~ler.docstring.ExampleClass.method1`.

        Notes
        -----
        Private methods typically do not need extensive docstrings unless the
        internal API is nontrivial. Keep it short and focused.
        """
        return None

    def method1(self, param1: Any = None, param2: Any = None) -> Any:
        r"""
        <one-line summary of what the method does.>

        <longer explanation. Mention algorithmic intent, assumptions, and edge
        cases. If applicable, note units and conventions.>

        Parameters
        ----------
        param1 : type, optional
            Parameter description. If ``None``, a default is inferred from
            :attr:`~ler.docstring.ExampleClass.param1`.
        param2 : type, optional
            Parameter description. If ``None``, a default is inferred from
            :attr:`~ler.docstring.ExampleClass.param2`.

        Returns
        -------
        type
            Return description.

        Raises
        ------
        ValueError
            Explanation of when/why it is raised.
        RuntimeError
            Explanation of when/why it is raised.

        See Also
        --------
        :meth:`~ler.docstring.ExampleClass.method2`
            Related method.

        Examples
        --------
        >>> from ler.docstring import ExampleClass
        >>> obj = ExampleClass(param1=1, param2={"key1": 2})
        >>> obj.method1()
        """
        self._method1_helper()
        return None

    def method2(self, param1: Any = None, param2: Any = None) -> Any:
        r"""
        <one-line summary of what the method does.>

        Parameters
        ----------
        param1 : type, optional
            Parameter description.
        param2 : type, optional
            Parameter description.

        Returns
        -------
        type
            Return description.

        Notes
        -----
        Add implementation notes only if they matter for correct use.

        Examples
        --------
        >>> from ler.docstring import ExampleClass
        >>> obj = ExampleClass(param1=1, param2={"key1": 2})
        >>> obj.method2()
        """
        return None
