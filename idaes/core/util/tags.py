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
"""This module containts utility classes that allow users to tag model quantities
and group them, for easy display, formatting, and input.
"""
import xml.dom.minidom
import collections

import pyomo.environ as pyo
from pyomo.core.base.indexed_component_slice import IndexedComponent_slice

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

__Author__ = "John Eslick"


class ModelTag:
    """The purpose of this class is to facilitate a simpler method of accessing,
    displaying, and reporting model quantities. The structure of IDAES models is
    a complex hierarchy. This class allows quantities to be accessed more directly
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
        "_set_in_display_units",
    ]

    def __init__(self, expr, format_string="{}", doc="", display_units=None):
        """initialize a model tag instance.
        Args:
            expr: A Pyomo Var, Expression, Param, Reference, or unnamed
                expression to tag. This can be a scalar or indexed.
            format_string: A formatting string used to print an element of the
                tagged expression (e.g. '{:.3f}').
            doc: A description of the tagged quantity.
            display_units: Pyomo units to display the quantity in. If a string
                is provided it will be used to display as the unit, but will not
                be used to convert units. If None, use native units of the
                quantity.
        """
        super().__init__()
        self._format = format_string  # format string for printing expression
        if isinstance(expr, IndexedComponent_slice):
            self._expression = pyo.Reference(expr)  # tag expression (can be unnamed)
        else:
            self._expression = expr
        self._doc = doc  # documentation for a tag
        self._display_units = display_units  # unit to display value in
        self._cache_validation_value = {}  # value when converted value stored
        self._cache_display_value = {}  # value to display after unit conversions
        self._name = None  # tag name (just used to clarify error messages)
        self._root = None  # use this to cache scalar tags in indexed parent
        self._index = None  # index to get cached converted value from parent
        self._group = None  # Group object if this is a member of a group
        self._str_units = True  # include units to stringify the tag
        # if _set_in_display_units is True and no units are provided for for
        # set, fix, setub, and setlb, the units will be assumed to be the
        # display units.  If it is false and no units are provided, the units are
        # assumed to be the native units of the quantity
        self._set_in_display_units = False

    def __getitem__(self, k):
        """Returns a new ModelTag for a scalar element of a tagged indexed
        quantity or a ModelTag with a slice as the expression."""
        if k is None and not self.is_indexed:
            return self
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
        if self._root is None:  # cache the unit conversion in root object
            tag._root = self
        else:
            tag._root = self._root
        tag._index = k
        return tag

    def __len__(self):
        """Number of elements in a tag"""
        return len(self.indexes)

    def keys(self):
        """Iterator of keys in the tag.  If scalar, is a single None."""
        indx_set = self.indexes
        for i in indx_set:
            yield i

    def values(self):
        """Iterator for scalar elements in a tag (The elements are scalar tags)"""
        indx_set = self.indexes
        for i in indx_set:
            yield self[i]

    def items(self):
        """Iterator for key scalar elements pairs in a tag"""
        indx_set = self.indexes
        for i in indx_set:
            yield i, self[i]

    def __str__(self):
        """Returns the default string representation of a tagged quantity. If
        the tagged expression is indexed this uses the default index.  This can
        be handy for things like the time index."""
        return self.display(units=self.str_include_units)

    def __call__(self, *args, **kwargs):
        """Calling an instance of a tagged quantity gets the display string see
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
        val = self.get_display_value(index=index)
        if format_string is None:
            format_string = self.get_format(units=False, index=index)
        if callable(format_string):
            format_string = format_string(val)
        if units:
            format_string = self._join_units(index=index, format_string=format_string)
        try:
            return format_string.format(val)
        except ValueError:
            # Probably trying to put a string through the number format
            # for various reasons, I'll allow it.
            return str(val)
        except TypeError:
            # Probably trying to put None through the numeric format.  This
            # can happen for example when variables don't have values.  I'll
            # allow 'None' to be printed.  It's not uncommon to happen, and I
            # don't want to raise an exception.
            return str(val)

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
    def value(self):
        """The numeric value of the tagged expression in the display units"""
        if not self.is_indexed:
            return self.get_display_value()
        else:
            val_dict = dict.fromkeys(self.indexes)
            for k, v in self.items():
                val_dict[k] = v.get_display_value()
            return val_dict

    @property
    def native_value(self):
        """The numeric value of the tagged expression in the native units"""
        if not self.is_indexed:
            return self.get_display_value(convert=False)
        else:
            val_dict = dict.fromkeys(self.indexes)
            for k, v in self.items():
                val_dict[k] = v.get_display_value(convert=False)
            return val_dict

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
        repeatedly.  The unconverted value is used to ensure the cached
        converted value is still valid.

        Args:
            index: index of value to display if expression is indexed
            convert: if False don't do unit conversion

        Returns:
            numeric expression value
        """
        if self._root is not None:
            if not self.is_indexed:
                index = self._index
            return self._root.get_display_value(index=index, convert=convert)

        try:
            expr = self.expression[index]
        except KeyError as key_err:
            if self._name is None:
                raise KeyError(f"{index} not a valid key for tag") from key_err
            raise KeyError(f"{index} not a valid key for tag {self._name}") from key_err
        except TypeError:
            expr = self.expression

        if (
            self._display_units is None
            or isinstance(self._display_units, str)
            or not convert
        ):
            # no display units, so display in original units, no convert opt
            try:
                return pyo.value(expr, exception=False)
            except ZeroDivisionError:
                return "ZeroDivisionError"
            except ValueError:  # no value
                return None

        try:
            val = pyo.value(expr, exception=False)
        except ZeroDivisionError:
            return "ZeroDivisionError"
        except ValueError:  # no value
            return None

        cache_validate = self._cache_validation_value
        cache_value = self._cache_display_value
        if index in cache_validate and val == cache_validate[index]:
            return cache_value[index]
        cache_validate[index] = val
        if self._display_units is None:
            cache_value[index] = pyo.value(expr, exception=False)
        else:
            try:
                cache_value[index] = pyo.value(
                    pyo.units.convert(expr, self._display_units), exception=False
                )
            except ValueError:
                return None
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
        """
        try:
            return issubclass(self._expression.ctype, pyo.Var)
        except AttributeError:
            return False

    @property
    def fixed(self):
        """Get the tagged variable if the tag is not a variable, raise TypeError"""
        if not self.is_var:
            return False
        return self.expression.fixed

    @property
    def is_indexed(self):
        """Returns whether the tagged expression is an indexed."""
        try:
            return self.expression.is_indexed()
        except AttributeError:
            return False

    @property
    def indexes(self):
        """The index set of the tagged quantity"""
        if self.is_indexed:
            return list(self.expression.keys())
        return (None,)

    @property
    def group(self):
        """The ModelTagGroup object that this belongs to, if any."""
        if self._root is not None:
            return self._root.group
        return self._group

    @group.setter
    def group(self, val):
        """The ModelTagGroup object that this belongs to, if any."""
        if self._root is not None:
            raise RuntimeError("group is superseded by the root property.")
        self._group = val

    @property
    def str_include_units(self):
        """Set whether to include units by default in the tag's string
        representation"""
        if self.group is not None:
            return self.group.str_include_units
        if self._root is not None:
            return self._root.str_include_units
        return self._str_units

    @str_include_units.setter
    def str_include_units(self, val):
        if self._group is not None:
            raise RuntimeError("str_include_units is superseded by the group property.")
        if self._root is not None:
            raise RuntimeError("str_include_units is superseded by the root property.")
        self._str_units = val

    @property
    def set_in_display_units(self):
        """Default index to use in the tag's string representation, this
        is required for indexed quantities if you want to automatically convert
        to string. An example use it for a time indexed tag, to display a
        specific time point."""
        if self.group is not None:
            return self.group.set_in_display_units
        if self._root is not None:
            return self._root.set_in_display_units
        return self._set_in_display_units

    @set_in_display_units.setter
    def set_in_display_units(self, val):
        """Default index to use in the tag's string representation, this
        is required for indexed quantities if you want to automatically convert
        to string. An example use it for a time indexed tag, to display a
        specific time point."""
        if self.group is not None:
            raise RuntimeError(
                "set_in_display_units is superseded by the group property."
            )
        if self._root is not None:
            raise RuntimeError(
                "set_in_display_units is superseded by the root property."
            )
        self._set_in_display_units = val

    def set(self, val, in_display_units=None):
        """Set the value of a tagged variable.

        Args:
            val: value
            in_display_units: if true assume the value is in the display units

        Returns:
            None
        """
        # If tag is a common immutable then allow it to be set
        if self.expression is None or isinstance(
            self.expression, (str, float, int, tuple)
        ):
            self._expression = val
            return

        if in_display_units is None:
            in_display_units = self.set_in_display_units

        if in_display_units and pyo.units.get_units(val) is None:
            if self._display_units is not None:
                val *= self._display_units

        try:
            try:
                self.expression.set_value(val)
            except ValueError:  # it's an indexed expression or slice
                for index in self.expression:
                    self.expression[index].set_value(val)
        except AttributeError as attr_err:
            if self._name:
                raise AttributeError(
                    f"Tagged expression {self._name}, has no set_value()."
                ) from attr_err
            raise AttributeError("Tagged expression has no set_value().") from attr_err

    def setlb(self, val, in_display_units=None):
        """Set the lower bound of a tagged variable.

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
            try:
                self.expression.setlb(val)
            except ValueError:  # it's an indexed expression or slice
                for index in self.expression:
                    self.expression[index].lb(val)
        except AttributeError as attr_err:
            if self._name:
                raise AttributeError(
                    f"Tagged expression {self._name}, has no setlb()."
                ) from attr_err
            raise AttributeError("Tagged expression has no setlb().") from attr_err

    def setub(self, val, in_display_units=None):
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
            try:
                self.expression.setub(val)
            except ValueError:  # it's an indexed expression or slice
                for index in self.expression:
                    self.expression[index].setub(val)
        except AttributeError as attr_err:
            if self._name:
                raise AttributeError(
                    f"Tagged expression {self._name}, has no setub()."
                ) from attr_err
            raise AttributeError("Tagged expression has no setub().") from attr_err

    def fix(self, val=None, in_display_units=None):
        """Fix the value of a tagged variable.

        Args:
            val: value, if None fix without setting a value
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
            if val is None:
                self.expression.fix()
            else:
                self.expression.fix(val)
        except AttributeError as attr_err:
            if self._name:
                raise AttributeError(f"Tag {self._name} has no fix.") from attr_err
            raise AttributeError("Tag has no fix.") from attr_err

    def unfix(self):
        """Unfix the value of a tagged variable.

        Args:
            v: value, if None fix without setting a value

        Returns:
            None
        """
        try:
            self.expression.unfix()
        except AttributeError as attr_err:
            if self._name:
                raise AttributeError(f"Tag, {self._name}, has no unfix.") from attr_err
            raise AttributeError("Tag has no unfix.") from attr_err

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
        "_str_units",
        "_set_in_display_units",
    ]

    def __init__(self):
        super().__init__()
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

    def _table_tagkey_index_lists(self, tags=None):
        """PRIVATE utility function to take tags (list of tag keys or 2-element
        list of a tag key and index) and turn it into a list of tag keys and list
        of indexes. Indexed tags with no index will be expanded to include all
        indexes.
        """

        if tags is None:
            tags = list(self.keys())

        tag_list = []
        indexes = []
        for i, tag in enumerate(tags):
            if not isinstance(tag, collections.abc.Hashable):
                if len(tag) == 2:
                    tag_list.append(tag[0])
                    indexes.append(tag[1])
                else:
                    raise ValueError("Key-index pairs should be a list of length 2")
            else:
                if not self[tag].is_indexed:
                    tag_list.append(tag)
                    indexes.append(None)
                else:
                    for k in self[tag].keys():
                        tag_list.append(tag)
                        indexes.append(k)
        return tag_list, indexes

    def table_heading(self, tags=None, units=True):
        """Create a table heading with a given set of tags, for tabulating model
        results.

        Args:
            tags: List (not tuple, since a tuple can be a key) of tag keys or
                2-element list of tags key and an index. If a key is for an
                indexed value and no index is given with the key, it will be
                flattened to create a column for each index.
            units: If True, include the units of measure in the column heading
                names.

        Returns:
            list of column headings
        """
        tag_list, indexes = self._table_tagkey_index_lists(tags)

        row = [None] * len(tag_list)
        for i, tag in enumerate(tag_list):
            if indexes[i] is None:
                name = tag
            else:
                name = f"{tag}[{indexes[i]}]"
            if units:
                row[i] = f"{name} ({self[tag].get_unit_str(index=indexes[i])})"
            else:
                row[i] = name
        return row

    def table_row(self, tags=None, units=True, numeric=False):
        """Create a table row with a given set of tags, for tabulating model
        results.  The row contains values of tagged quantities either as strings
        of numeric values.

        Args:
            tags: List (not tuple, since a tuple can be a key) of tag keys or
                2-element list of tags key and an index. If a key is for an
                indexed value and no index is given with the key, it will be
                flattened to create a column for each index.
            units: If true include the units of measure in the value, if the
                values are not numeric.
            numeric: If true provide numeric values rather than a formatted
                string.

        Returns:
            list of values for a table row either numeric or as formatted strings
        """
        tag_list, indexes = self._table_tagkey_index_lists(tags)

        row = [None] * len(tag_list)
        for i, tag in enumerate(tag_list):
            if numeric:
                row[i] = self[tag][indexes[i]].value
            else:
                row[i] = self[tag][indexes[i]].display(units=units)
        return row

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


def svg_tag(
    svg=None,
    tag_group=None,
    outfile=None,
    idx=None,
    tag_map=None,
    show_tags=False,
    byte_encoding="utf-8",
):
    """
    Replace text in a SVG with tag values for the model. This works by looking
    for text elements in the SVG with IDs that match the tags or are in tag_map.

    Args:
        svg: a file pointer or a string containing svg contents
        tag_group: a ModelTagGroup with tags to display in the SVG
        outfile: a file name to save the results, if None don't save
        idx: if None not indexed, otherwise an index in the indexing set of the
            reference
        tag_map: dictionary with svg id keys and tag values, to map svg ids to
            tags, used in cases where tags contain characters that cannot be used
            in the svg's xml
        show_tags: Put tag labels of the diagram instead of numbers
        byte_encoding: If svg is given as a byte-array, use this encoding to
            convert it to a string.

    Returns:
        SVG String
    """
    if svg is None:
        raise RuntimeError("svg string or file-like object required.")

    # get SVG content string
    if isinstance(svg, str):  # already a string
        pass
    elif isinstance(svg, bytes):  # bytes to string
        svg = svg.decode(byte_encoding)  # file-like to string
    elif hasattr(svg, "read"):
        svg = svg.read()
    else:  # Can't handle whatever this is.
        raise TypeError("SVG must either be a string or a file-like object")

    # Make tag map here because the tags may not make valid XML IDs if no
    # tag_map provided we'll go ahead and handle XML @ (maybe more in future)
    if tag_map is None:
        tag_map = dict()
        for tag in tag_group:
            new_tag = tag.replace("@", "_")
            new_tag = new_tag.replace(" ", "_")
            tag_map[new_tag] = tag

    # Ture SVG string into XML document
    doc = xml.dom.minidom.parseString(svg)
    # Get the text elements of the SVG
    texts = doc.getElementsByTagName("text")

    # Add some text
    for t in texts:
        id = t.attributes["id"].value
        if id in tag_map:
            # if it's multiline change last line
            try:
                tspan = t.getElementsByTagName("tspan")[-1]
            except IndexError:
                _log.warning(f"Text object but no tspan for tag {tag_map[id]}.")
                _log.warning(f"Skipping output for {tag_map[id]}.")
                continue
            try:
                tspan = tspan.childNodes[0]
            except IndexError:
                # No child node means there is a line with no text, so add some.
                tspan.appendChild(doc.createTextNode(""))
                tspan = tspan.childNodes[0]

            if show_tags:
                val = tag_map[id]
            else:
                if tag_group[tag_map[id]].is_indexed:
                    val = tag_group[tag_map[id]][idx]
                else:
                    val = tag_group[tag_map[id]]

            tspan.nodeValue = str(val)

    new_svg = doc.toxml()
    # If outfile is provided save to a file
    if outfile is not None:
        with open(outfile, "w") as f:
            f.write(new_svg)
    # Return the SVG as a string.  This lets you take several passes at adding
    # output without saving and loading files.
    return new_svg
