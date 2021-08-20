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

"""
This module contains miscellaneous utility functions for use in IDAES models.
"""
import xml.dom.minidom

import pyomo.environ as pyo
from pyomo.core.base.expression import _GeneralExpressionData
from pyomo.core.base.component import ModelComponentFactory
from pyomo.core.base.indexed_component import (
    UnindexedComponent_set, )
from pyomo.core.base.disable_methods import disable_methods
from pyomo.common.config import ConfigBlock

import idaes.logger as idaeslog
import idaes.core.solvers

_log = idaeslog.getLogger(__name__)


# Author: Andrew Lee
def get_solver(solver=None, options=None):
    """
    General method for getting a solver object which defaults to the standard
    IDAES solver (defined in the IDAES configuration).

    Args:
        solver: string name for desired solver. Default=None, use default solver
        options: dict of solver options to use, overwrites any settings
                 provided by IDAES configuration. Default = None, use default
                 solver options.

    Returns:
        A Pyomo solver object
    """
    if solver is None:
        solver = "default"
    solver_obj = idaes.core.solvers.SolverWrapper(solver, register=False)()

    if options is not None:
        solver_obj.options.update(options)

    return solver_obj


# Author: Andrew Lee
def add_object_reference(self, local_name, remote_object):
    """
    Method to create a reference in the local model to a remote Pyomo object.
    This method should only be used where Pyomo Reference objects are not
    suitable (such as for referencing scalar Pyomo objects where the None
    index is undesirable).

    Args:
        local_name : name to use for local reference (str)
        remote_object : object to make a reference to

    Returns:
        None
    """
    try:
        object.__setattr__(self, local_name, remote_object)
    except AttributeError:
        raise AttributeError(
            "{} failed to construct reference to {} - remote "
            "object does not exist.".format(self.name, remote_object)
        )


# Author: Jaffer Ghouse
def extract_data(data_dict):
    """
    General method that returns a rule to extract data from a python
    dictionary. This method allows the param block to have a database for
    a parameter but extract a subset of this data to initialize a Pyomo
    param object.
    """

    def _rule_initialize(m, *args):
        if len(args) > 1:
            return data_dict[args]
        else:
            return data_dict[args[0]]

    return _rule_initialize


# Author: John Eslick
def TagReference(s, description=""):
    """
    Create a Pyomo reference with an added description string attribute to
    describe the reference. The intended use for these references is to create
    a time-indexed reference to variables in a model corresponding to plant
    measurment tags.

    Args:
        s: Pyomo time slice of a variable or expression
        description (str): A description the measurment

    Returns:
        A Pyomo Reference object with an added doc attribute
    """
    r = pyo.Reference(s)
    r.description = description
    return r


# Author John Eslick
def svg_tag(
    tags,
    svg,
    outfile=None,
    idx=None,
    tag_map=None,
    show_tags=False,
    byte_encoding="utf-8",
    tag_format=None,
    tag_format_default="{:.4e}"
):
    """
    Replace text in a SVG with tag values for the model. This works by looking
    for text elements in the SVG with IDs that match the tags or are in tag_map.

    Args:
        tags: A dictionary where the key is the tag and the value is a Pyomo
            Reference.  The reference could be indexed. In typical IDAES
            applications the references would be indexed by time.
        svg: a file pointer or a string continaing svg contents
        outfile: a file name to save the results, if None don't save
        idx: if None not indexed, otherwise an index in the indexing set of the
            reference
        tag_map: dictionary with svg id keys and tag values, to map svg ids to
            tags
        show_tags: Put tag labels of the diagram instead of numbers
        byte_encoding: If svg is given as a byte-array, use this encoding to
            convert it to a string.
        tag_format: A dictionary of formatting strings.  If the formatting
            string is a callable, it should be a function that takes the value
            to display and returns a formatting string.
        tag_format_default: The default formatting if not explicitly by
            tag_format. If the formatting string is a callable, it should be a
            function that takes the value to display and returns a formatting
            string.

    Returns:
        SVG String
    """
    if tag_format is None:
        tag_format = {}

    if isinstance(svg, str):  # assume this is svg content string
        pass
    elif isinstance(svg, bytes):
        svg = svg.decode(byte_encoding)
    elif hasattr(svg, "read"):  # file-like object to svg
        svg = svg.read()
    else:
        raise TypeError("SVG must either be a string or a file-like object")
    # Make tag map here because the tags may not make valid XML IDs if no
    # tag_map provided we'll go ahead and handle XML @ (maybe more in future)
    if tag_map is None:
        tag_map = dict()
        for tag in tags:
            new_tag = tag.replace("@", "_")
            tag_map[new_tag] = tag
    # Search for text in the svg that has an id in tags
    doc = xml.dom.minidom.parseString(svg)
    texts = doc.getElementsByTagName("text")
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
            try:
                if show_tags:
                    val = tag_map[id]
                elif idx is None:
                    val = pyo.value(tags[tag_map[id]], exception=False)
                else:
                    val = pyo.value(tags[tag_map[id]][idx], exception=False)
            except ZeroDivisionError:
                val = "Divide_by_0"
            tf = tag_format.get(tag_map[id], tag_format_default)
            try:
                if callable(tf): # conditional formatting
                    tspan.nodeValue = tf(val).format(val)
                else:
                    tspan.nodeValue = tf.format(val)
            except ValueError:
                # whatever it is, it doesn't match the format.  Usually this
                # happens when a string is given, but it is using a default
                # number format
                tspan.nodeValue = val

    new_svg = doc.toxml()
    # If outfile is provided save to a file
    if outfile is not None:
        with open(outfile, "w") as f:
            f.write(new_svg)
    # Return the SVG as a string.  This lets you take several passes at adding
    # output without saving and loading files.
    return new_svg


# Author: John Eslick
def copy_port_values(destination=None, source=None, arc=None,
        direction="forward"):
    """
    Moved to idaes.core.util.initialization.propagate_state.
    Leaving redirection function here for deprecation warning.
    """
    _log.warning("DEPRECATED: copy_port_values has been deprecated. "
            "The same functionality can be found in "
            "idaes.core.util.initialization.propagate_state.")
    from idaes.core.util.initialization import propagate_state
    propagate_state(destination=destination, source=source, arc=arc,
            direction=direction)


def set_param_from_config(b, param, config=None, index=None):
    """
    Utility method to set parameter value from a config block. This allows for
    converting units if required. This method directly sets the value of the
    parameter.

    Args:
        b - block on which parameter and config block are defined
        param - name of parameter as str. Used to find param and config arg
        units - units of param object (used if conversion required)
        config - (optional) config block to get parameter data from. If
                unset, assumes b.config.
        index - (optional) used for pure component properties where a single
                property may have multiple parameters associated with it.

    Returns:
        None
    """
    if config is None:
        try:
            config = b.config
        except AttributeError:
            raise AttributeError(
                "{} - set_param_from_config method was not provided with a "
                "config argument, but no default Config block exists. Please "
                "specify the Config block to use via the config argument."
                .format(b.name))

    # Check that config is an instance of a Config Block
    if not isinstance(config, ConfigBlock):
        raise TypeError(
            "{} - set_param_from_config - config argument provided is not an "
            "instance of a Config Block.".format(b.name))

    if index is None:
        try:
            param_obj = getattr(b, param)
        except AttributeError:
            raise AttributeError(
                "{} - set_param_from_config method was provided with param "
                "argument {}, but no attribute of that name exists."
                .format(b.name, param))

        try:
            p_data = config.parameter_data[param]
        except (KeyError, AttributeError):
            raise KeyError(
                "{} - set_param_from_config method was provided with param "
                "argument {}, but the config block does not contain a "
                "value for this parameter.".format(b.name, param))
    else:
        try:
            param_obj = getattr(b, param+"_"+index)
        except AttributeError:
            raise AttributeError(
                "{} - set_param_from_config method was provided with param and"
                " index arguments {} {}, but no attribute with that "
                "combination ({}_{}) exists."
                .format(b.name, param, index, param, index))

        try:
            p_data = config.parameter_data[param][index]
        except (KeyError, AttributeError):
            raise KeyError(
                "{} - set_param_from_config method was provided with param and"
                " index arguments {} {}, but the config block does not contain"
                " a value for this parameter and index."
                .format(b.name, param, index))

    units = param_obj.get_units()

    if isinstance(p_data, tuple):
        # 11 Dec 2020 - There is currently a bug in Pyomo where trying to
        # convert the units of a unitless quantity results in a TypeError.
        # To avoid this, we check here for cases where both the parameter and
        # user provided value are unitless and bypass unit conversion.
        if ((units is None or units is pyo.units.dimensionless) and
                (p_data[1] is None or p_data[1] is pyo.units.dimensionless)):
            param_obj.value = p_data[0]
        else:
            param_obj.value = pyo.units.convert_value(
                p_data[0], from_units=p_data[1], to_units=units)
    else:
        _log.debug("{} no units provided for parameter {} - assuming default "
                   "units".format(b.name, param))
        param_obj.value = p_data


# -----------------------------------------------------------------------------
# Creating a Component derived from Pyomo's Expression to use in cases
# where an Expression could be mistaken for a Var.
# Author: Andrew Lee
class _GeneralVarLikeExpressionData(_GeneralExpressionData):
    """
    An object derived from _GeneralExpressionData which implements methods for
    common APIs on Vars.

    Constructor Arguments:
        expr: The Pyomo expression stored in this expression.

        component: The Expression object that owns this data.

    Public Class Attributes:
        expr: The expression owned by this data.

    Private class attributes:
        _component: The expression component.
    """

    # Define methods for common APIs on Vars in case user mistakes
    # an Expression for a Var
    @property
    def value(self):
        # Overload value so it behaves like a Var
        return pyo.value(self.expr)

    @value.setter
    def value(self, expr):
        # Overload value seter to prevent users changing the expression body
        raise TypeError(
            "%s is an Expression and does not have a value which can be set."
            % (self.name))

    def setlb(self, val=None):
        raise TypeError(
            "%s is an Expression and can not have bounds. "
            "Use an inequality Constraint instead."
            % (self.name))

    def setub(self, val=None):
        raise TypeError(
            "%s is an Expression and can not have bounds. "
            "Use an inequality Constraint instead."
            % (self.name))

    def fix(self, val=None):
        raise TypeError(
            "%s is an Expression and can not be fixed. "
            "Use an equality Constraint instead."
            % (self.name))

    def unfix(self):
        raise TypeError(
            "%s is an Expression and can not be unfixed."
            % (self.name))

@ModelComponentFactory.register(
    "Named expressions that can be used in places of variables.")
class VarLikeExpression(pyo.Expression):
    """
    A shared var-like expression container, which may be defined over a index.

    Constructor Arguments:
        initialize: A Pyomo expression or dictionary of expressions used
        to initialize this object.

        expr: A synonym for initialize.

        rule: A rule function used to initialize this object.
    """

    _ComponentDataClass = _GeneralVarLikeExpressionData
    NoConstraint    = (1000,)
    Skip            = (1000,)

    def __new__(cls, *args, **kwds):
        if cls is not VarLikeExpression:
            return super(VarLikeExpression, cls).__new__(cls)
        if not args or (args[0] is UnindexedComponent_set and len(args) == 1):
            return super(VarLikeExpression, cls).__new__(
                AbstractSimpleVarLikeExpression)
        else:
            return super(VarLikeExpression, cls).__new__(
                IndexedVarLikeExpression)


class SimpleVarLikeExpression(_GeneralVarLikeExpressionData,
                              VarLikeExpression):

    def __init__(self, *args, **kwds):
        _GeneralVarLikeExpressionData.__init__(self, expr=None, component=self)
        VarLikeExpression.__init__(self, *args, **kwds)

    #
    # From Pyomo: Leaving this method for backward compatibility reasons.
    # (probably should be removed)
    # Note: Doesn't seem to work without it
    #
    def add(self, index, expr):
        """Add an expression with a given index."""
        if index is not None:
            raise KeyError(
                "SimpleExpression object '%s' does not accept "
                "index values other than None. Invalid value: %s"
                % (self.name, index))
        if (type(expr) is tuple) and \
           (expr == pyo.Expression.Skip):
            raise ValueError(
                "Expression.Skip can not be assigned "
                "to an Expression that is not indexed: %s"
                % (self.name))
        self.set_value(expr)
        return self


@disable_methods({'set_value', 'is_constant', 'is_fixed', 'expr'})
class AbstractSimpleVarLikeExpression(SimpleVarLikeExpression):
    pass


class IndexedVarLikeExpression(VarLikeExpression):

    #
    # From Pyomo: Leaving this method for backward compatibility reasons
    # Note: It allows adding members outside of self._index.
    #       This has always been the case. Not sure there is
    #       any reason to maintain a reference to a separate
    #       index set if we allow this.
    #
    def add(self, index, expr):
        """Add an expression with a given index."""
        if (type(expr) is tuple) and (expr == pyo.Expression.Skip):
            return None
        cdata = _GeneralVarLikeExpressionData(expr, component=self)
        self._data[index] = cdata
        return cdata

    # Define methods for common APIs on Vars in case user mistakes
    # an Expression for a Var
    def setlb(self, val=None):
        raise TypeError(
            "%s is an Expression and can not have bounds. "
            "Use inequality Constraints instead."
            % (self.name))

    def setub(self, val=None):
        raise TypeError(
            "%s is an Expression and can not have bounds. "
            "Use inequality Constraints instead."
            % (self.name))

    def fix(self, val=None):
        raise TypeError(
            "%s is an Expression and can not be fixed. "
            "Use equality Constraints instead."
            % (self.name))

    def unfix(self):
        raise TypeError(
            "%s is an Expression and can not be unfixed."
            % (self.name))
