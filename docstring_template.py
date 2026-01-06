# -*- coding: utf-8 -*-
"""
Docstring Template based on LeR style conventions.

This template follows NumPy/SciPy docstring conventions with reStructuredText 
(RST) formatting for Sphinx documentation. It demonstrates the standard format 
used throughout the ``ler`` package.

Style Guidelines:
    - Use double backticks (``type``) for type annotations in docstrings
    - Use ``\n`` for explicit line breaks in RST output
    - Use 4-column format for Instance Attributes tables (Attribute, Type, Unit, Description)
    - Use 2-column format for Instance Methods tables (Method, Description)
    - Private methods (prefixed with ``_``) should not appear in Instance Methods tables
    - Properties should have docstrings with Returns section only (not Parameters)

Copyright (C) 2026 Author Name. Distributed under MIT License.

AI prompt to use this template: 
1. Check for any inconsistencies and redundancies in the code (don't touch the docstrings). Please point them out without making any changes.
2. Rewrite docstrings (including inline docstrings) of *.py using docstring_template.py. Convert some of the Instance Methods to private methods with underscore prefix and remove them from the class docstring. Use '\n' for newline and don't use '\\n'. Dont add unnecessary docstrings like the following:
    # =============================================================================
    # MODULE DOCSTRING TEMPLATE
    # =============================================================================
3. Create properties of the Instance Attributes if it is not there. Properties should be at the end.
"""

# =============================================================================
# MODULE DOCSTRING TEMPLATE
# =============================================================================
"""
Module for [brief description].

[Extended description explaining the module's purpose, main classes/functions,
and how it fits into the larger package structure.]

Inheritance hierarchy (if applicable):

- :class:`~package.ParentClass` \n
  - :class:`~package.GrandparentClass` \n

Usage:
    Basic workflow example:

    >>> from package import ClassName
    >>> obj = ClassName()
    >>> result = obj.some_method()

Copyright (C) 2026 Author Name. Distributed under MIT License.
"""


# =============================================================================
# CLASS DOCSTRING TEMPLATE
# =============================================================================
class TemplateClass:
    """
    One-line summary of the class.

    Extended description of the class explaining its purpose, functionality,
    and how it relates to other classes in the package. This should be 2-4
    sentences that provide context.

    Key Features: \n
    - Feature one description \n
    - Feature two description \n
    - Feature three description \n

    Parameters
    ----------
    param1 : ``int``
        Description of param1. \n
        default: 4
    param2 : ``float``
        Description of param2. \n
        default: 0.0
    param3 : ``str``
        Description of param3 with options. \n
        Options: \n
        - 'option_a': Description of option A \n
        - 'option_b': Description of option B \n
        - 'option_c': Description of option C \n
        default: 'option_a'
    param4 : ``dict`` or ``None``
        Description of dict parameter. \n
        default: dict(key1="value1", key2="value2")
    param5 : ``callable`` or ``None``
        Description of callable parameter. \n
        If None, uses default behavior. \n
        The function should follow the signature: \n
        ``def func(arg1, arg2): return result`` \n
        default: None
    **kwargs : ``dict``
        Additional keyword arguments passed to parent classes: \n
        :class:`~package.ParentClass1`, \n
        :class:`~package.ParentClass2`.

    Examples
    --------
    Basic usage:

    >>> from package import TemplateClass
    >>> obj = TemplateClass()
    >>> result = obj.some_method()


    Instance Methods
    ----------
    TemplateClass has the following methods: \n
    +-----------------------------------------------------+------------------------------------------------+
    | Method                                              | Description                                    |
    +=====================================================+================================================+
    | :meth:`~method_one`                                 | Brief description of method one                |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~method_two`                                 | Brief description of method two                |
    +-----------------------------------------------------+------------------------------------------------+
    | :meth:`~method_three`                               | Brief description of method three              |
    +-----------------------------------------------------+------------------------------------------------+

    Instance Attributes
    ----------
    TemplateClass has the following attributes: \n
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | Attribute                                      | Type             | Unit  | Description                                    |
    +================================================+==================+=======+================================================+
    | :attr:`~attr1`                                 | ``int``          |       | Description of attr1                           |
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | :attr:`~attr2`                                 | ``float``        | Hz    | Description of attr2 with units                |
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | :attr:`~attr3`                                 | ``str``          |       | Description of attr3                           |
    +------------------------------------------------+------------------+-------+------------------------------------------------+
    | :attr:`~attr4`                                 | ``dict``         |       | Description of attr4                           |
    +------------------------------------------------+------------------+-------+------------------------------------------------+

    Notes
    -----
    - Note 1 about the class. \n
    - Note 2 providing additional context. \n
    - References to related classes or documentation. \n
    """

    def __init__(self, param1=4, param2=0.0, param3="option_a", param4=None, **kwargs):
        self.param1 = param1
        self.param2 = param2
        self.param3 = param3
        self.param4 = param4 or {}


# =============================================================================
# METHOD DOCSTRING TEMPLATES
# =============================================================================

    def public_method(self, arg1, arg2=None, arg3=0.5):
        """
        One-line summary of the method.

        Extended description explaining what the method does, any important
        behavior, and how parameters affect the output.

        Parameters
        ----------
        arg1 : ``dict`` or ``str``
            Description of arg1. \n
            default: None (uses self.some_attribute)
        arg2 : ``float``
            Description of arg2. \n
            default: 0.5
        arg3 : ``str``
            Description with options. \n
            Options: \n
            - 'option1': Description of option1 \n
            - 'option2': Description of option2 \n
            default: 'option1'

        Returns
        -------
        result1 : ``float``
            Description of first return value (with units if applicable).
        result2 : ``dict``
            Description of second return value.

        Examples
        --------
        >>> from package import TemplateClass
        >>> obj = TemplateClass()
        >>> result1, result2 = obj.public_method(arg1=data)
        """
        pass

    def _private_helper_method(self, param):
        """
        Helper function to [brief description].
        
        Parameters
        ----------
        param : ``dict`` or ``str``
            Description of the parameter.
        
        Returns
        -------
        result : ``dict``
            Description of return value.
        """
        pass


# =============================================================================
# PROPERTY DOCSTRING TEMPLATES
# =============================================================================

    @property
    def simple_property(self):
        """
        One-line description of the property.

        Returns
        -------
        simple_property : ``int``
            Description of what this property represents. \n
            default: 4
        """
        return self._simple_property

    @simple_property.setter
    def simple_property(self, value):
        self._simple_property = value

    @property
    def property_with_units(self):
        """
        Property representing a physical quantity.

        Returns
        -------
        property_with_units : ``float``
            Description (units: Hz, Mpc, etc.). \n
            default: 1.0
        """
        return self._property_with_units

    @property_with_units.setter
    def property_with_units(self, value):
        self._property_with_units = value

    @property
    def enum_property(self):
        """
        Property with enumerated string options.

        Returns
        -------
        enum_property : ``str``
            Description of the property. \n
            Options: \n
            - 'option_a': Description of option A \n
            - 'option_b': Description of option B \n
            - 'option_c': Description of option C \n
            default: 'option_a'
        """
        return self._enum_property

    @enum_property.setter
    def enum_property(self, value):
        self._enum_property = value

    @property
    def dict_property(self):
        """
        Configuration dictionary property.

        Returns
        -------
        dict_property : ``dict``
            Dictionary with keys: \n
            - 'key1': Description of key1 \n
            - 'key2': Description of key2 \n
            - 'key3': Description of key3 \n
        """
        return self._dict_property

    @dict_property.setter
    def dict_property(self, value):
        self._dict_property = value

    @property
    def list_property(self):
        """
        List of objects property.

        Returns
        -------
        list_property : ``list``
            List of items (e.g., detector names, parameter keys).
        """
        return self._list_property

    @list_property.setter
    def list_property(self, value):
        self._list_property = value

    @property
    def callable_property(self):
        """
        Callable function property.

        Returns
        -------
        callable_property : ``callable``
            Function that performs some operation. \n
            The function signature should be: \n
            ``func(arg1, arg2) -> result``
        """
        return self._callable_property

    @callable_property.setter
    def callable_property(self, value):
        self._callable_property = value

    @property
    def cosmology_property(self):
        """
        Astropy cosmology object for calculations.

        Returns
        -------
        cosmology_property : ``astropy.cosmology``
            Cosmology used for distance and volume calculations. \n
            default: LambdaCDM(H0=70, Om0=0.3, Ode0=0.7)
        """
        return self._cosmology_property

    @cosmology_property.setter
    def cosmology_property(self, value):
        self._cosmology_property = value


# =============================================================================
# FORMATTING REFERENCE
# =============================================================================
"""
QUICK REFERENCE:

1. Types: Use double backticks
   - ``int``, ``float``, ``str``, ``bool``, ``dict``, ``list``
   - ``callable``, ``numpy.ndarray``, ``astropy.cosmology``
   - Combined types: ``dict`` or ``str``, ``float`` or ``list``

2. Line breaks: Use \n at end of lines for RST rendering
   description line 1. \n
   description line 2. \n

3. Options format:
   Options: \n
   - 'option1': Description \n
   - 'option2': Description \n
   default: 'option1'

4. Defaults format:
   default: value
   default: None (uses self.attribute_name)
   default: dict(key1="val1", key2="val2")

5. Cross-references:
   :class:`~package.ClassName`
   :meth:`~method_name`
   :attr:`~attribute_name`

6. Section headers (7 dashes for Returns, more for sections):
   Parameters
   ----------
   
   Returns
   -------
   
   Examples
   --------

7. RST Tables (for Instance Methods and Attributes):
   +-----------------------------------------------------+------------------------------------------------+
   | Column1                                             | Column2                                        |
   +=====================================================+================================================+
   | :meth:`~method_name`                                | Description                                    |
   +-----------------------------------------------------+------------------------------------------------+
"""
