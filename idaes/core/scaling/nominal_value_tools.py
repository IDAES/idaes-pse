#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
The NominalValueExtractionVisitor and supporting functions.

Authors: Andrew Lee, Douglas Allan
"""

import math

from pyomo.environ import (
    Binary,
    Boolean,
    NegativeIntegers,
    NegativeReals,
    NonNegativeIntegers,
    NonNegativeReals,
    NonPositiveIntegers,
    NonPositiveReals,
    PositiveIntegers,
    PositiveReals,
    value,
)
from pyomo.core.base.var import VarData
from pyomo.core.base.expression import ExpressionData
from pyomo.core.base.param import ParamData
from pyomo.core import expr as EXPR
from pyomo.core.base.units_container import _PyomoUnit

from pyomo.common.numeric_types import native_types

from idaes.core.scaling.util import get_scaling_factor
from idaes.core.util.exceptions import BurntToast
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

TAB = " " * 4


def get_nominal_value(component):
    """
    Get the signed nominal value for a VarData or ParamData component.

    For Params, the current value of the component will be returned.

    For Vars, the nominal value is determined using the assigned scaling factor
    and the sign determined based on the bounds and domain of the variable (defaulting to
    positive). If no scaling factor is set, then the current value will be used if set,
    otherwise the absolute nominal value will be equal to 1.

    Args:
        component: component to determine nominal value for

    Returns:
        signed float with nominal value

    Raises:
        TypeError if component is not instance of VarData or ParamData
    """
    # Determine if Var or Param
    if isinstance(component, VarData):
        # Get scaling factor for Var
        sf = get_scaling_factor(component)
        if sf is None:
            # No scaling factor - see if Var has a value
            if component.value is not None:
                # If it has a value, use that as the nominal value
                # As we are using the actual value, do not need to determine sign
                return value(component)
            else:
                # Otherwise assign a nominal value of 1
                sf = 1

        # Try to determine expected sign of node
        ub = component.ub
        lb = component.lb
        domain = component.domain

        # To avoid NoneType errors, assign dummy values in place of None
        if ub is None:
            # No upper bound, take a positive value
            ub = 1000
        if lb is None:
            # No lower bound, take a negative value
            lb = -1000

        if lb >= 0 or domain in [
            NonNegativeReals,
            PositiveReals,
            PositiveIntegers,
            NonNegativeIntegers,
            Boolean,
            Binary,
        ]:
            # Strictly positive
            sign = 1
        elif ub <= 0 or domain in [
            NegativeReals,
            NonPositiveReals,
            NegativeIntegers,
            NonPositiveIntegers,
        ]:
            # Strictly negative
            sign = -1
        else:
            # Unbounded, see if there is a current value
            # Assume positive until proven otherwise
            sign = 1
            if component.value is not None:
                val = value(component)
                if val < 0:
                    # Assigned negative value, assume value will remain negative
                    sign = -1

        return sign / sf

    elif isinstance(component, ParamData):
        # Nominal value of a parameter is always its value
        return value(component)
    else:
        # Not a Var or Param - invalid component type
        raise TypeError(
            f"get_nominal_value - {component.name} is not an instance of a Var or Param."
        )


class NominalValueExtractionVisitor(EXPR.StreamBasedExpressionVisitor):
    """
    Expression walker for collecting scaling factors in an expression and determining the
    expected value of the expression using the scaling factors as nominal inputs.

    By default, the get_nominal_value method is used to determine the nominal value for
    all Vars and Params in the expression, however this can be changed by setting the
    nominal_value_callback argument.

    Returns a list of expected values for each additive term in the expression.

    In order to properly assess the expected value of terms within functions, the sign
    of each term is maintained throughout thus returned values may be negative. Functions
    using this walker should handle these appropriately.
    """

    def __init__(self, nominal_value_callback=get_nominal_value):
        """
        Visitor class used to determine nominal values of all terms in an expression based on
        scaling factors assigned to the associated variables. Do not use this class directly.

        Args:
            nominal_value_callback - method to use to get nominal value of root nodes.

        Notes
        -----
        This class inherits from the :class:`StreamBasedExpressionVisitor` to implement
        a walker that returns the nominal value corresponding to all additive terms in an
        expression.
        There are class attributes (dicts) that map the expression node type to the
        particular method that should be called to return the nominal value of the node based
        on the nominal value of its child arguments. This map is used in exitNode.
        """
        super().__init__()

        self._nominal_value_callback = nominal_value_callback

    def _get_magnitude_base_type(self, node):
        try:
            return [self._nominal_value_callback(node)]
        except TypeError:
            # Not a Var or Param - something went wrong
            raise BurntToast(
                "NominalValueExtractionVisitor found root node that was not a Var or Param. "
                "This should never happen - please contact the developers with this bug."
            )

    def _get_nominal_value_for_sum_subexpression(self, child_nominal_values):
        return sum(i for i in child_nominal_values)

    def _get_nominal_value_for_sum(self, node, child_nominal_values):
        # For sums, collect all child values into a list
        mag = []
        for i in child_nominal_values:
            for j in i:
                mag.append(j)
        return mag

    def _get_nominal_value_for_product(self, node, child_nominal_values):
        mag = []
        for i in child_nominal_values[0]:
            for j in child_nominal_values[1]:
                mag.append(i * j)
        return mag

    def _get_nominal_value_for_division(self, node, child_nominal_values):
        numerator = self._get_nominal_value_for_sum_subexpression(
            child_nominal_values[0]
        )
        denominator = self._get_nominal_value_for_sum_subexpression(
            child_nominal_values[1]
        )
        if denominator == 0:
            # Assign a nominal value of 1 so that we can continue
            denominator = 1
            # Log a warning for the user
            _log.warning(
                "Nominal value of 0 found in denominator of division expression. "
                "Assigning a value of 1. You should check you scaling factors and models to "
                "ensure there are no values of 0 that can appear in these functions."
            )
        return [numerator / denominator]

    def _get_nominal_value_for_power(self, node, child_nominal_values):
        # Use the absolute value of the base term to avoid possible complex numbers
        base = abs(
            self._get_nominal_value_for_sum_subexpression(child_nominal_values[0])
        )
        exponent = self._get_nominal_value_for_sum_subexpression(
            child_nominal_values[1]
        )

        return [base**exponent]

    def _get_nominal_value_single_child(self, node, child_nominal_values):
        return child_nominal_values[0]

    def _get_nominal_value_abs(self, node, child_nominal_values):
        return [abs(i) for i in child_nominal_values[0]]

    def _get_nominal_value_negation(self, node, child_nominal_values):
        return [-i for i in child_nominal_values[0]]

    def _get_nominal_value_for_unary_function(self, node, child_nominal_values):
        func_name = node.getname()
        func_nominal = self._get_nominal_value_for_sum_subexpression(
            child_nominal_values[0]
        )
        func = getattr(math, func_name)
        try:
            return [func(func_nominal)]
        except ValueError:
            raise ValueError(
                f"Evaluation error occurred when getting nominal value in {func_name} "
                f"expression with input {func_nominal}. You should check you scaling factors "
                f"and model to address any numerical issues or scale this constraint manually."
            )

    def _get_nominal_value_expr_if(self, node, child_nominal_values):
        return child_nominal_values[1] + child_nominal_values[2]

    def _get_nominal_value_external_function(self, node, child_nominal_values):
        # First, need to get expected magnitudes of input terms, which may be sub-expressions
        input_mag = []
        for i in child_nominal_values:
            if isinstance(i[0], str):
                # Sometimes external functions might have string arguments
                # Check here, and return the string if true
                input_mag.append(i[0])
            else:
                input_mag.append(self._get_nominal_value_for_sum_subexpression(i))

        # Next, create a copy of the external function with expected magnitudes as inputs
        newfunc = node.create_node_with_local_data(input_mag)

        # Evaluate new function and return the absolute value
        return [value(newfunc)]

    node_type_method_map = {
        EXPR.EqualityExpression: _get_nominal_value_for_sum,
        EXPR.InequalityExpression: _get_nominal_value_for_sum,
        EXPR.RangedExpression: _get_nominal_value_for_sum,
        EXPR.SumExpression: _get_nominal_value_for_sum,
        EXPR.NPV_SumExpression: _get_nominal_value_for_sum,
        EXPR.ProductExpression: _get_nominal_value_for_product,
        EXPR.MonomialTermExpression: _get_nominal_value_for_product,
        EXPR.NPV_ProductExpression: _get_nominal_value_for_product,
        EXPR.DivisionExpression: _get_nominal_value_for_division,
        EXPR.NPV_DivisionExpression: _get_nominal_value_for_division,
        EXPR.PowExpression: _get_nominal_value_for_power,
        EXPR.NPV_PowExpression: _get_nominal_value_for_power,
        EXPR.NegationExpression: _get_nominal_value_negation,
        EXPR.NPV_NegationExpression: _get_nominal_value_negation,
        EXPR.AbsExpression: _get_nominal_value_abs,
        EXPR.NPV_AbsExpression: _get_nominal_value_abs,
        EXPR.UnaryFunctionExpression: _get_nominal_value_for_unary_function,
        EXPR.NPV_UnaryFunctionExpression: _get_nominal_value_for_unary_function,
        EXPR.Expr_ifExpression: _get_nominal_value_expr_if,
        EXPR.ExternalFunctionExpression: _get_nominal_value_external_function,
        EXPR.NPV_ExternalFunctionExpression: _get_nominal_value_external_function,
        EXPR.LinearExpression: _get_nominal_value_for_sum,
    }

    def beforeChild(self, node, child, child_idx):
        """
        Callback for :class:`pyomo.core.current.StreamBasedExpressionVisitor`. This method
        is called before entering a child node. If we encounter a named Expression with
        a scaling hint, then we use that scaling hint instead of descending further into
        the expression tree.
        """
        if isinstance(child, ExpressionData):
            sf = get_scaling_factor(child, warning=False)
            if sf is not None:
                # Crude way to determine sign of expression. Maybe fbbt could be used here?
                try:
                    val = value(child)
                except ValueError:
                    # Some variable isn't defined, etc.
                    val = 1
                if val < 0:
                    return (False, [-1 / sf])
                else:
                    return (False, [1 / sf])
        return (True, None)

    def exitNode(self, node, data):
        """Callback for :class:`pyomo.core.current.StreamBasedExpressionVisitor`. This
        method is called when moving back up the tree in a depth first search."""

        # first check if the node is a leaf
        nodetype = type(node)

        if nodetype in native_types:
            return [node]

        node_func = self.node_type_method_map.get(nodetype, None)
        if node_func is not None:
            return node_func(self, node, data)

        elif not node.is_expression_type():
            # this is a leaf, but not a native type
            if nodetype is _PyomoUnit:
                return [1]
            else:
                return self._get_magnitude_base_type(node)
                # might want to add other common types here

        # not a leaf - check if it is a named expression
        if (
            hasattr(node, "is_named_expression_type")
            and node.is_named_expression_type()
        ):
            return self._get_nominal_value_single_child(node, data)

        raise TypeError(
            f"An unhandled expression node type: {str(nodetype)} was encountered while "
            f"retrieving the nominal value of expression {str(node)}"
        )
