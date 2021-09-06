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

__Author__ = "John Eslick"

import pyomo.environ as pyo
from pyomo.core.base.indexed_component_slice import IndexedComponent_slice


class ModelTag(object):
    """The purpose of this class is to facilitate a simpler method of accessing
    displaying and reporting model quantities. The structure of IDAES models is
    a complex hiarachy. The class allows quanties to be accessed more directly
    and provides more control over how they are reported."""

    def __init__(self, expr, format="{}", doc="", name=None, display_units=None):
        super().__init__()
        self._format = format  # format string for printing expression
        self._expression = expr  # tag expression (can be unnamed)
        self._doc = doc  # documentation for a tag
        self._display_units = display_units  # unit to display value in
        self._cache_validation_value = {}  # value when converted value stored
        self._cache_display_value = {}  # value to display after unit conversions
        self._name = name  # tag name (just used to claify error messages)
        self._str_units = True  # include units to stringify the tag
        self._str_index = 0  # use this index to stringify the tag
        self._root = None  # use this to cache scalar tags in indexed parent
        self._index = None  # index to get cached converted value from parent

    def __getitem__(self, k):
        """Returns a new model tag for a scalar element of a tagged indexed
        quantity"""
        try:
            tag = ModelTag(
                expr=self.expression[k],
                format=self._format,
                display_units=self._display_units,
                doc=self._doc,
            )
        except KeyError:
            if self._name is None:
                raise KeyError(f"{index} is not a valid index for tag")
            else:
                raise KeyError(f"{index} is not a valid index for tag {self._name}")
        tag._str_units = self._str_units
        tag._str_index = self._str_index
        if (
            self._root is None or self.is_slice
        ):  # cache the unit conversion in root object
            tag._root = self
        else:
            tag._root = self._root
        tag._index = k
        return tag

    def __str__(self):
        """Returns the default string representation of a tagged quantitiy"""
        return self.display(units=self._str_units, index=self._str_index)

    def __call__(self, *args, **kwargs):
        """Calling an instance of a tagged quantitiy gets the display string see
        display()"""
        return self.display(*args, **kwargs)

    def display(self, units=True, format=None, index=None):
        """Get a string representation of the tagged quantity

        Args:
            units (bool): Include units of measure in the string
            format (str): Formatting string, if supplied overrides the default
            index: If the tagged quantity is indexed, an index for the element
                to display is required, or the default index is used.

        Returns:
            str
        """
        if self.is_indexed and index is None:
            index = self._str_index
        v = self.get_display_value(index=index)
        if format is None:
            return self.get_format(units=units, index=index).format(v)
        else:
            if units:
                return self._join_units(format=format, index=index).format(v)
            else:
                return format.format(v)

    def _join_units(self, index=None, format=None):
        """Private method to join the format string with the units of measure
        string.
        """
        if format is None:
            format = self._format
        u = self.get_unit_str(index=index)
        if u == "None" or u == "" or u is None:
            return format
        elif u == "%":
            return "".join([format, u])
        else:
            return " ".join([format, u])

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
        else:
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
                e = self.expression[index]
            except KeyError:
                if self._name is None:
                    raise KeyError(f"{index} is not a valid index for tag")
                else:
                    raise KeyError(f"{index} is not a valid index for tag {self._name}")
        else:
            e = self.expression
        if self._display_units is None:
            # no display units, so display in original units
            return pyo.value(e)
        elif not convert or isinstance(self._display_units, str):
            # display units can be a string, if you have a model without units
            # attached but still want to display units, you can specify them as
            # a string, just for printing
            return pyo.value(e)
        else:
            v = pyo.value(e)
            if self._root is None:
                cv = self._cache_validation_value
                cd = self._cache_display_value
            else:
                cv = self._root._cache_validation_value
                cd = self._root._cache_display_value
                if not self.is_indexed:
                    index = self._index
            if v == cv.get(index, None):
                return cd[index]
            else:
                cv[index] = v
                cd[index] = pyo.value(pyo.units.convert(e, self._display_units))
                return cd[index]

    def get_unit_str(self, index=None):
        """String representation of the tagged quantity's units of measure"""
        if self._display_units is None:
            if self.is_indexed:
                return str(pyo.units.get_units(self.expression[index]))
            else:
                return str(pyo.units.get_units(self.expression))
        else:
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
        else:
            return None

    @property
    def indices(self):
        """The index set of the tagged qunatity"""
        if self.is_indexed:
            return list(self.expression.keys())
        else:
            return None

    def str_include_units(self, b=True):
        """Set weather to include units by default in the tag's string
        representation"""
        self._str_units = b

    def str_default_index(self, index=0):
        """Default index to use in the tag's string representation, this
        is required for indexed quntities if you want to automatically convert
        to string. An example use it for a time indexed tag, to display a
        specific time point."""
        self._str_index = index

    def set(self, v, in_display_units=False):
        """Set the value of a tagged variable.

        Args:
            v: value
            in_display_units: if true assume the value is in the display units

        Returns:
            None
        """
        if in_display_units:
            v *= self._display_units

        def _setv(c, v):
            try:
                c.value = v
            except AttributeError:
                raise AttributeError(
                    f"Tagged expression {self._name}, has no attribute 'value'."
                )

        if self.is_slice:
            for c in self.expression:
                _setv(c, v)
        else:
            _setv(self.expression, v)

    def fix(self, v=None, in_display_units=False):
        """Fix the value of a tagged variable.

        Args:
            v: value, if None fix without setting a value
            in_display_units: if true assume the value is in the display units

        Returns:
            None
        """
        if in_display_units:
            v *= self._display_units

        def _fix(c, v):
            try:
                if v is None:
                    c.fix()
                else:
                    c.fix(v)
            except AttributeError:
                raise AttributeError(
                    f"Tagged expression {self._name}, has no fix method."
                )

        if self.is_slice:
            for c in self.expression:
                _fix(c, v)
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
            except AttributeError:
                raise AttributeError(
                    f"Tagged expression {self._name}, has no unfix method."
                )

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

    def __init__(self):
        super().__init__()

    def add(self, name, expr, format="{}", display_units=None, doc=""):
        if isinstance(expr, ModelTag):
            self[name] = expr
            self[name]._name = name
        else:
            self[name] = ModelTag(
                expr=expr,
                format=format,
                display_units=display_units,
                name=name,
                doc=doc,
            )

    def expr_dict(self):
        d = {}
        for n, v in self.items():
            d[name] = v.expression
        return d

    def format_dict(self, units=True):
        d = {}
        for n, v in self.items():
            d[n] = v.format(units=units)
        return d

    def str_include_units(self, b=True):
        """When convertig a tag in this group directly to a string, include units
        or not.

        Args:
            b (bool): if True include units in the string
        """
        for v in self.values():
            v.str_include_units(b)

    def str_default_index(self, index=0):
        """When convertig a tag in this group directly to a string, use index as
        the default index for indexed varaibles. This can be useful when the tag
        group where indexed expressions have a consitent index set.  For example,
        where all the indexed expressions are indexed by time, and you want to
        display results at a specific time point.

        Args:
            index: default index for indexed expressions when stringifying a tag

        Returns:
            None
        """
        for v in self.values():
            v.str_default_index(index)

    def names_with_units_dict(self, format="{name} [{units}]"):
        """Returns a dictionary with tag names as keys and tags with units for
        values. This may be useful for tabulating results of multiple model
        results among other things.

        Args:
            format: The format of the tag with units string by default the formating
                is '{name} [{units}]'

        """
        d = {}
        for n, v in self.items():
            d[name] = format.format(name=n, units=v.unit_str)
        return d
