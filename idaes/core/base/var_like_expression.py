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
Creating a Component derived from Pyomo's Expression to use in cases
where an Expression could be mistaken for a Var.
"""

import pyomo.environ as pyo
from pyomo.core.base.expression import _GeneralExpressionData
from pyomo.core.base.component import ModelComponentFactory
from pyomo.core.base.indexed_component import (
    UnindexedComponent_set,
)
from pyomo.core.base.disable_methods import disable_methods


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
    def set_value(self, value, force=False):
        """
        Overload set_value method to provide meaningful error if user attempts
        to set the value of the Expression. In order to support changing the
        expression (and setting it originally), if self._expr is None or
        force=True, the value of the expression will be updated, otherwise a
        TypeError will be raised.

        Args:
            value: value to set for _expr
            force: force updating of _expr if True (default = False)

        Returns:
            None

        Raises:
            TypeError if _expr is not None and force=False
        """
        if self._expr is None or force:
            super().set_value(value)
        else:
            raise TypeError(
                f"{self.name} is an Expression and does not have a value "
                f"which can be set."
            )

    @property
    def value(self):
        raise TypeError(
            f"{self.name} is an Expression and does not have a value "
            f"attribute. Use the 'value()' method instead."
        )

    @value.setter
    def value(self, expr):
        raise TypeError(
            "%s is an Expression and does not have a value which can be set."
            % (self.name)
        )

    def setlb(self, val=None):
        raise TypeError(
            "%s is an Expression and can not have bounds. "
            "Use an inequality Constraint instead." % (self.name)
        )

    def setub(self, val=None):
        raise TypeError(
            "%s is an Expression and can not have bounds. "
            "Use an inequality Constraint instead." % (self.name)
        )

    def fix(self, val=None):
        raise TypeError(
            "%s is an Expression and can not be fixed. "
            "Use an equality Constraint instead." % (self.name)
        )

    def unfix(self):
        raise TypeError("%s is an Expression and can not be unfixed." % (self.name))


@ModelComponentFactory.register(
    "Named expressions that can be used in places of variables."
)
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
    NoConstraint = (1000,)
    Skip = (1000,)

    def __new__(cls, *args, **kwds):
        if cls is not VarLikeExpression:
            return super(VarLikeExpression, cls).__new__(cls)
        if not args or (args[0] is UnindexedComponent_set and len(args) == 1):
            return super(VarLikeExpression, cls).__new__(
                AbstractSimpleVarLikeExpression
            )
        else:
            return super(VarLikeExpression, cls).__new__(IndexedVarLikeExpression)


class SimpleVarLikeExpression(_GeneralVarLikeExpressionData, VarLikeExpression):
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
                "index values other than None. Invalid value: %s" % (self.name, index)
            )
        if (type(expr) is tuple) and (expr == pyo.Expression.Skip):
            raise ValueError(
                "Expression.Skip can not be assigned "
                "to an Expression that is not indexed: %s" % (self.name)
            )
        self.set_value(expr)
        return self


@disable_methods({"set_value", "is_constant", "is_fixed", "expr"})
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
            "Use inequality Constraints instead." % (self.name)
        )

    def setub(self, val=None):
        raise TypeError(
            "%s is an Expression and can not have bounds. "
            "Use inequality Constraints instead." % (self.name)
        )

    def fix(self, val=None):
        raise TypeError(
            "%s is an Expression and can not be fixed. "
            "Use equality Constraints instead." % (self.name)
        )

    def unfix(self):
        raise TypeError("%s is an Expression and can not be unfixed." % (self.name))
