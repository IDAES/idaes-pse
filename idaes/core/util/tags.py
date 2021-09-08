#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""This module containts utility classes that allow users to tag model quntities
and put them in groups, for easy display formatting and model input.
"""

import pyomo.environ as pyo
from pyomo.core.base.indexed_component_slice import IndexedComponent_slice

__Author__ = "John Eslick"


class ModelTag:
    """The purpose of this class is to facilitate a simpler method of accessing,
    displaying, and reporting model quantities. The structure of IDAES models is
    a complex hiarachy. The class allows quanties to be accessed more directly
    and provides more control over how they are reported."""

    __slots__ = [
        "_format",
        "_expression",
        "_doc",
        "_display_units",
        "_cache_validation_value",
        "_cache_display_value",
        "_name",
        "_root",
        "_index",
        "_group",
        "_str_units",
        "_str_index",
        "_set_in_display_units",
    ]

    def __init__(self, expr, format_string="{}", doc="", display_units=None):
        """initialize a model tag instance.
        Args:
            expr: A Pyomo Var, Expression, Param, Reference, or unnamed
                expression to tag. This can be a scalar or indexed.
            format_string: A formating string to print an elememnt of the tagged
                expression (e.g. '{:.3f}' prints a float with 3 decimal places).
            doc: A description of the tagged qunatity.
            display_units: Pyomo units to display the qunatity in.
        """
        super().__init__()
        self._format = format_string  # format string for printing expression
        self._expression = expr  # tag expression (can be unnamed)
        self._doc = doc  # documentation for a tag
        self._display_units = display_units  # unit to display value in
        self._cache_validation_value = {}  # value when converted value stored
        self._cache_display_value = {}  # value to display after unit conversions
        self._name = None  # tag name (just used to claify error messages)
        self._root = None  # use this to cache scalar tags in indexed parent
        self._index = None  # index to get cached converted value from parent
        self._group = None  # Group object if this is a member of a group
        self._str_units = True  # include units to stringify the tag
        self._str_index = 0  # use this index to stringify the tag
        self._set_in_display_units = False  # if a value is rovided with no units
        # assume display units for fix and set

    def __getitem__(self, k):
        """Returns a new model tag for a scalar element of a tagged indexed
        quantity"""
        try:
            tag = ModelTag(
                expr=self.expression[k],
                format_string=self._format,
                display_units=self._display_units,
                doc=self._doc,
            )
        except KeyError as key_err:
            if self._name is None:
                raise KeyError(f"{k} is not a valid index for tag") from key_err
            raise KeyError(
                f"{k} is not a valid index for tag {self._name}"
            ) from key_err
        tag._str_units = self._str_units
        tag._str_index = self._str_index
        tag._set_in_display_units = self._set_in_display_units
        if (
            self._root is None or self.is_slice
        ):  # cache the unit conversion in root object
            tag._root = self
        else:
            tag._root = self._root
        tag._index = k
        tag._group = self._group
        return tag

    def __str__(self):
        """Returns the default string representation of a tagged quantity. If
        the tagged expression is indexed this uses the default index.  This can
        be handy for things like the time index."""
        if self._group is None:
            include_units = self._str_units
            expr_index = self._str_index  # this is ignored if not indexed
        else:
            include_units = self._group._str_units
            expr_index = self._group._str_index
        return self.display(units=include_units, index=expr_index)

    def __call__(self, *args, **kwargs):
        """Calling an instance of a tagged quantitiy gets the display string see
        display()"""
        return self.display(*args, **kwargs)

    def display(self, units=True, format_string=None, index=None):
        """Get a string representation of the tagged quantity

        Args:
            units (bool): Include units of measure in the string
            format_string (str): Formatting string, if supplied overrides the default
            index: If the tagged quantity is indexed, an index for the element
                to display is required, or the default index is used.

        Returns:
            str
        """
        if self.is_indexed and index is None:
            if self._group is not None:
                index = self._group.str_default_index
            else:
                index = self._str_index
        val = self.get_display_value(index=index)
        if format_string is None:
            return self.get_format(units=units, index=index).format(val)
        if units:
            return self._join_units(format_string=format_string, index=index).format(
                val
            )
        return format_string.format(val)

    def _join_units(self, index=None, format_string=None):
        """Private method to join the format string with the units of measure
        string.
        """
        if format_string is None:
            format_string = self._format
        units = self.get_unit_str(index=index)
        if units == "None" or units == "" or units is None:
            return format_string
        if units == "%":
            return "".join([format_string, units])
        return " ".join([format_string, units])

    @property
    def expression(self):
        """The tagged expression"""
        return self._expression

    @property
    def doc(self):
        """Tag documentation string"""
        return self._doc

    def get_format(self, units=True, index=None):
        """Get the formatting string.

        Args:
            units: if True include units of measure
            index: index for indexed expressions

        Returns
            str:
        """
        if units:
            return self._join_units(index=index)
        return self._format

    def get_display_value(self, index=None, convert=True):
        """Get the value of the expression to display.  Do unit conversion if
        needed.  This caches the unit conversion, to save time if this is called
        repeatededly.  The unconverted value is used to ensure the cached
        converted value is still valid.

        Args:
            index: index of value to display if expression is indexed
            convert: if False don't do unit conversion

        Returns:
            numeric expression value
        """
        if self.is_indexed:
            try:
                expr = self.expression[index]
            except KeyError as key_err:
                if self._name is None:
                    raise KeyError(f"{index} is not a valid index for tag") from key_err
                raise KeyError(
                    f"{index} is not a valid index for tag {self._name}"
                ) from key_err
        else:
            expr = self.expression
        if self._display_units is None:
            # no display units, so display in original units
            return pyo.value(expr)
        if not convert or isinstance(self._display_units, str):
            # display units can be a string, if you have a model without units
            # attached but still want to display units, you can specify them as
            # a string, just for printing
            return pyo.value(expr)
        val = pyo.value(expr)
        if self._root is None:
            cache_validate = self._cache_validation_value
            cache_value = self._cache_display_value
        else:
            cache_validate = self._root._cache_validation_value
            cache_value = self._root._cache_display_value
            if not self.is_indexed:
                index = self._index
        if val == cache_validate.get(index, None):
            return cache_value[index]
        cache_validate[index] = val
        cache_value[index] = pyo.value(pyo.units.convert(expr, self._display_units))
        return cache_value[index]

    def get_unit_str(self, index=None):
        """String representation of the tagged quantity's units of measure"""
        if self._display_units is None:
            if self.is_indexed:
                return str(pyo.units.get_units(self.expression[index]))
            return str(pyo.units.get_units(self.expression))
        return str(self._display_units)

    @property
    def is_var(self):
        """Whether the tagged expression is a Pyomo Var. Tagged variables
        can be fixed or set, while expressions cannot.

        Args:
            None

        Returns:
            True if tagged expression is a variable
        """
        try:
            return issubclass(self._expression.ctype, pyo.Var)
        except AttributeError:
            return False

    @property
    def is_slice(self):
        """Whether the tagged expression is a Pyomo slice.

        Args:
            None

        Returns:
            True if tagged expression is a Pyomo component slice
        """
        try:
            return isinstance(self._expression, IndexedComponent_slice)
        except AttributeError:
            return False

    @property
    def is_indexed(self):
        """Returns whether the tagged expression is a indexed.

        Args:
            None

        Returns:
            True if tagged expression is a variable
        """
        try:
            return self.expression.is_indexed()
        except AttributeError:
            return False

    @property
    def indexes(self):
        """The index set of the tagged quantity"""
        if self.is_indexed:
            return list(self.expression.keys())
        return None

    @property
    def indices(self):
        """The index set of the tagged qunatity"""
        if self.is_indexed:
            return list(self.expression.keys())
        return None

    @property
    def str_include_units(self):
        """Set whether to include units by default in the tag's string
        representation"""
        if self._group is None:
            return self._str_units
        return self._group._str_units

    @str_include_units.setter
    def str_include_units(self, val):
        if self._group is not None:
            raise RuntimeError("str_include_units is superseded by the group property.")
        self._str_units = val

    @property
    def str_default_index(self):
        """Default index to use in the tag's string representation, this
        is required for indexed quntities if you want to automatically convert
        to string. An example use it for a time indexed tag, to display a
        specific time point."""
        if self._group is None:
            return self._str_index
        return self._group._str_index

    @str_default_index.setter
    def str_default_index(self, index):
        """Default index to use in the tag's string representation, this
        is required for indexed quntities if you want to automatically convert
        to string. An example use it for a time indexed tag, to display a
        specific time point."""
        if self._group is not None:
            raise RuntimeError("str_default_index is superseded by the group property.")
        self._str_index = index

    @property
    def set_in_display_units(self):
        """Default index to use in the tag's string representation, this
        is required for indexed quntities if you want to automatically convert
        to string. An example use it for a time indexed tag, to display a
        specific time point."""
        if self._group is None:
            return self._set_in_display_units
        return self._group.set_in_display_units

    @set_in_display_units.setter
    def set_in_display_units(self, b):
        """Default index to use in the tag's string representation, this
        is required for indexed quntities if you want to automatically convert
        to string. An example use it for a time indexed tag, to display a
        specific time point."""
        if self._group is not None:
            raise RuntimeError(
                "set_in_display_units is superseded by the group property."
            )
        self._set_in_display_units = b

    def set(self, val, in_display_units=None):
        """Set the value of a tagged variable.

        Args:
            v: value
            in_display_units: if true assume the value is in the display units

        Returns:
            None
        """
        if in_display_units is None:
            in_display_units = self.set_in_display_units

        if in_display_units and pyo.units.get_units(val) is None:
            if self._display_units is not None:
                val *= self._display_units

        try:
            self.expression.value = val
        except AttributeError as attr_err:
            if self._name is None:
                raise AttributeError(
                    f"Tagged expression {self._name}, has no attribute 'value'."
                ) from attr_err
            raise AttributeError(
                "Tagged expression has no attribute 'value'."
            ) from attr_err

    def fix(self, v=None, in_display_units=None):
        """Fix the value of a tagged variable.

        Args:
            val: value, if None fix without setting a value
            in_display_units: if true assume the value is in the display units

        Returns:
            None
        """
        if in_display_units is None:
            in_display_units = self.set_in_display_units

        if in_display_units and pyo.units.get_units(v) is None:
            if not (
                isinstance(self._display_units, str) or self._display_units is None
            ):
                v *= self._display_units

        def _fix(comp, v):
            try:
                if v is None:
                    comp.fix()
                else:
                    comp.fix(v)
            except AttributeError as attr_err:
                raise AttributeError(
                    f"Tagged expression has no fix method ({self._name})"
                ) from attr_err

        if self.is_slice:
            for comp in self.expression:
                _fix(comp, v)
        else:
            _fix(self.expression, v)

    def unfix(self):
        """Unfix the value of a tagged variable.

        Args:
            v: value, if None fix without setting a value

        Returns:
            None
        """

        def _unfix(c):
            try:
                c.unfix()
            except AttributeError as attr_err:
                raise AttributeError(
                    f"Tagged expression {self._name}, has no unfix method."
                ) from attr_err

        if self.is_slice:
            for c in self.expression:
                _unfix(c)
        else:
            _unfix(self.expression)

    @property
    def var(self):
        """Get the tagged variable if the tag is not a variable, raise TypeError"""
        if not self.is_var:
            raise TypeError(
                "Can only return a variable if the expression is a variable."
            )
        return self.expression


class ModelTagGroup(dict):
    """This a basically a dictionary of ModelTag objects. This is used to group
    sets of tags, and contains methods to operate on sets of tags. It inherits dict
    so dictionary methods can be used."""

    __slots__ = [
        "_str_index",
        "_str_units",
        "_set_in_display_units",
    ]

    def __init__(self):
        super().__init__()
        self._str_index = 0
        self._str_units = True
        self._set_in_display_units = False

    def __setitem__(self, key, val):
        if isinstance(val, ModelTag):
            super().__setitem__(key, val)
            self[key]._name = key
            self[key]._group = self
        else:
            raise TypeError("Only ModelTag objects can be part of a ModelTagGroup.")

    def add(self, name, expr, **kwargs):
        """Add a new model tag to the group"""
        if isinstance(expr, ModelTag):
            self[name] = expr
        else:
            self[name] = ModelTag(expr=expr, **kwargs)

    def expr_dict(self):
        """Get a dictionary of expressions with tag keys."""
        d = {}
        for name, v in self.items():
            d[name] = v.expression
        return d

    def format_dict(self, units=True):
        """Get a dictionary of format strings with tag keys."""
        d = {}
        for n, v in self.items():
            d[n] = v.format(units=units)
        return d

    @property
    def str_include_units(self):
        """When converting a tag in this group directly to a string, include units
        or not.
        """
        return self._str_units

    @str_include_units.setter
    def str_include_units(self, value):
        """When converting a tag in this group directly to a string, include units
        or not.
        """
        self._str_units = value

    @property
    def set_in_display_units(self):
        """When this is True, and set() or fix() are called on a tag in the group,
        with a value that doesn't include units, assume the display units.
        """
        return self._set_in_display_units

    @set_in_display_units.setter
    def set_in_display_units(self, value):
        """When this is True, and set() or fix() are called on a tag in the group,
        with a value that doesn't include units, assume the display units.
        """
        self._set_in_display_units = value

    @property
    def str_default_index(self):
        """When converting a tag in this group directly to a string, use index as
        the default index for indexed varaibles. This can be useful when the tag
        group where indexed expressions have a consitent index set.  For example,
        where all the indexed expressions are indexed by time, and you want to
        display results at a specific time point.
        """
        return self._str_index

    @str_default_index.setter
    def str_default_index(self, value):
        """When converting a tag in this group directly to a string, use index as
        the default index for indexed varaibles. This can be useful when the tag
        group where indexed expressions have a consitent index set.  For example,
        where all the indexed expressions are indexed by time, and you want to
        display results at a specific time point.
        """
        self._str_index = value

    def names_with_units_dict(self, format_string="{name} [{units}]"):
        """Returns a dictionary with tag names as keys and tags with units for
        values. This may be useful for tabulating results of multiple model
        results among other things.

        Args:
            format_string: The format of the tag with units string by default the formating
                is '{name} [{units}]'

        """
        d = {}
        for name, v in self.items():
            d[name] = format_string.format(name, v.unit_str)
        return d
