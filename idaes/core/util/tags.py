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
    displaying and reporting model quantities. The structure of IDAES models can be
    a cmplex hiarachy. The class allows quanties to be accessed more directly and
    provides control over how they are reported."""

    def __init__(self, expr, format="{}", doc="", units=None, name=None):
        super().__init__()
        self._format = format  # format string for printing expression
        self._expression = expr  # tag expression (can be unnamed)
        self._units = units  # user supplied or otherwise obtained from pyomo
        self._doc = doc  # documentation for a tag
        self._name = name  # tag name (just used to claify error messages)
        self._str_units = True  # include units to stringify the tag
        self._str_index = 0  # use this index to stringify the tag

    def __getitem__(self, k):
        """Returns a new model tag for a scalar element of a tagged indexed
        quantity"""
        try:
            tag = ModelTag(
                expr=self.expression[k],
                format=self._format,
                units=self._units,
                doc=self._doc,
            )
        except KeyError:
            if self._name is None:
                raise KeyError(f"{index} is not a valid index for tag")
            else:
                raise KeyError(f"{index} is not a valid index for tag {self._name}")
        tag._str_units = self._str_units
        # here you may expect for copy the default str index, but the new tag
        # object shouldn't be indexed anymore, so it doesn't matter.
        return tag

    def __str__(self):
        """Returns the default string representation of a tagged quantitiy"""
        return self.display(units=self._str_units, index=self._str_index)

    def __call__(self, *args, **kwargs):
        """Calling an instance of a tagged quantitiy get the display string see
        display()"""
        return self.display(*args, **kwargs)

    def display(self, units=True, format=None, index=None):
        """Get a string representation of the tagged quantity

        Args:
            units (bool): Include units of measure in the string
            format (str): Formatting string, if supplied overrides the default
            index: If the tagged quantity is indexed, and index for the element
                to display is required.

        Returns:
            str
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
        if format is None:
            return self.get_format(units=units).format(pyo.value(e))
        else:
            if units:
                return self._join_units(format=format).format(pyo.value(e))
            else:
                return format.format(pyo.value(e))

    def _join_units(self, format=None):
        """Private method to join the format string with the units of measure
        string.
        """
        if format is None:
            format = self._format
        u = self.unit_str
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

    def get_format(self, units=True):
        """Get the formatting string.

        Args:
            units: if True include units of measure

        Returns
            str:
        """
        if units:
            return self._join_units()
        else:
            return self._format

    @property
    def unit_str(self):
        """String representation of the tagged quantity's units of measure"""
        if self._units is not None:
            return self._units
        else:
            if self.is_indexed: # issue with references assume units are same
                first = self.expression.index_set().first()
                self._units = str(pyo.units.get_units(self.expression[first]))
            else:
                self._units = str(pyo.units.get_units(self.expression))
            #raise TypeError(self._units)
            return self._units

    @property
    def is_var(self):
        """Returns whether the tagged expression is a Pyomo Var. Tagged variables
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
        """Returns whether the tagged expression is a Pyomo slice.

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
        to string."""
        self._str_index = index

    def set(self, v, index=None):
        """Set the value of a tagged variable.

        Args:
            v: value
            index: if variable is indexed specifies the index to set

        Returns:
            None
        """

        def _setv(c, v):
            try:
                c.value = v
            except AttributeError:
                raise AttributeError(
                    f"Tagged expression {self._name}, has attribute 'value'."
                )

        if self.is_slice:
            for c in self.expression:
                _setv(c, v)
        else:
            _setv(self.expression, v)

    def fix(self, v=None):
        """Fix the value of a tagged variable.

        Args:
            v: value, if None fix without setting a value

        Returns:
            None
        """

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

    def add(self, name, expr, format="{}", units=None, doc=""):
        if isinstance(expr, ModelTag):
            self[name] = expr
            self[name]._name = name
        else:
            self[name] = ModelTag(
                expr=expr, format=format, units=units, name=name, doc=doc
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

    def str_index(self, index=0):
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


if __name__ == "__main__":
    """This contains some examples"""
    m = pyo.ConcreteModel()
    m.x = pyo.Var([1, 2], initialize={1: 3, 2: 4}, units=pyo.units.kg)
    m.y = pyo.Var(initialize=5, units=pyo.units.s)
    m.e = pyo.Expression(expr=m.x[1] / m.y)
    m.z = pyo.Var(initialize=5)
    m.w = pyo.Var(initialize=5, units=pyo.units.kg)

    g = ModelTagGroup()

    g.add("e", ModelTag(expr=m.e, format="{:.3f}"))
    g.add("z", expr=m.z, format="{:.1f}")
    g.add("x", expr=m.x, format="{}")
    g.add("y", expr=m.y, format="{}")
    g.add("f", expr=m.x[1] / m.w * 100, units="%")

    for t, v in g.items():
        print(f"{t} is var: {v.is_var}")
        print(f"{t} is indexed: {v.is_indexed}")
        print(f"{t} indexes: {v.indexes}")

    # Various ways to display a tag

    # indexed
    try:
        print(f"{g['x']}")
    except KeyError:
        print("I didn't provide a valid index")
    g["x"].str_default_index(1)
    print(f"{g['x']}")
    print(f"{g['x'][1]}")
    print(f"{g['x'](index=1)}")
    print(f"{g['x'][1](units=False)}")
    print(f"{g['x'][1].display(units=False, index=2)}")

    # not indexed
    print(f"{g['e']}")
    print(f"{g['e'].display(units=False)}")
    print(f"{g['e'].display(units=False, format='{:.1e}')}")

    # no units
    print(f"{g['y']}")
    print(f"{g['y'].display(format='{:.1e}')}")

    # user supplied units
    print(f"{g['f']}")

    # use tags to set value

    g["x"][1].fix(1)
    g["x"][2].set(4)
    g["x"].expression.display()

    # another way
    g["x"][1].var.fix(2)
    g["x"][2].var.value = 6
    g["x"].var.display()

    g["x"][:].fix(1)
    g["x"].var.display()

    g["x"][:].set(2)
    g["x"].var.display()

    g["x"][:].unfix()
    g["x"].var.display()
