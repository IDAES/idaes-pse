# -*- coding: utf-8 -*-
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
This module contains a a visitor for analyzing constraint and identifying
numerically problematic terms.
"""

__author__ = "Alexander Dowling, Douglas Allan, Andrew Lee, Robby Parker, Ben Knueven"

from math import log10
from itertools import combinations, chain

from pyomo.environ import (
    ComponentMap,
    value,
)
from pyomo.core.base.var import VarData
from pyomo.core import expr as EXPR
from pyomo.common.numeric_types import native_types
from pyomo.core.base.units_container import _PyomoUnit

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


class ConstraintTermAnalysisVisitor(EXPR.StreamBasedExpressionVisitor):
    """
    Expression walker for checking Constraints for problematic terms.

    This walker will walk the expression and look for summation terms
    with mismatched magnitudes or potential cancellations.

    Args:
        term_mismatch_tolerance: tolerance to use when determining mismatched
            terms
        term_cancellation_tolerance: tolerance to use when identifying
            possible cancellation of terms
        term_zero_tolerance: tolerance for considering terms equal to zero
        max_canceling_terms: maximum number of terms to consider when looking
            for canceling combinations (None = consider all possible combinations)
        max_cancellations_per_node: maximum number of cancellations to collect
            for a single node. Collection will terminate when this many cancellations
            have been identified (None = collect all cancellations)

    Returns:
        list of values for top-level summation terms
        list of terms with mismatched magnitudes
        list of terms with potential cancellations
        bool indicating whether expression is a constant
    """

    def __init__(
        self,
        term_mismatch_tolerance: float = 1e6,
        term_cancellation_tolerance: float = 1e-4,
        term_zero_tolerance: float = 1e-10,
        max_canceling_terms: int = 4,
        max_cancellations_per_node: int = 5,
    ):
        super().__init__()

        # Tolerance attributes
        self._log_mm_tol = log10(term_mismatch_tolerance)
        self._sum_tol = term_cancellation_tolerance
        self._zero_tolerance = term_zero_tolerance
        self._max_canceling_terms = max_canceling_terms
        self._max_cancellations_per_node = max_cancellations_per_node

        # Placeholders for collecting results
        self.canceling_terms = ComponentMap()
        self.mismatched_terms = ComponentMap()

        # Flag for if cancellation collection hit limit
        self._cancellation_tripped = False

    def _get_value_for_sum_subexpression(self, child_data):
        # child_data is a tuple, with the 0-th element being the node values
        if isinstance(child_data[0][0], str):
            # Values may be a list containing a string in some cases (e.g. external functions)
            # Return the string in this case
            return child_data[0][0]
        return sum(i for i in child_data[0])

    def _generate_combinations(self, inputs, equality=False):
        # We want to test all combinations of terms for cancellation

        # The number of combinations we check depends on whether this is an (in)equality
        # expression or a sum node deeper in the expression tree.
        # We expect (in)equalities to generally sum to 0 (0 == expr) thus we want to
        # check if any subset of the sum terms sum to zero (i.e. are any terms unnecessary).
        # For other sum nodes, we need to check for any combination of terms.

        # Maximum number of terms to include in combinations
        max_comb = len(inputs)
        if equality:
            # Subtract 1 if (in)equality node
            max_comb += -1
        # We also have a limit on the maximum number of terms to consider
        if self._max_canceling_terms is not None:
            max_comb = min(max_comb, self._max_canceling_terms)

        # Single terms cannot cancel, thus we want all combinations of length 2 to max terms
        # Note the need for +1 due to way range works
        combo_range = range(2, max_comb + 1)

        # Collect combinations of terms in an iterator
        # In order to identify the terms in each cancellation, we will pair each value
        # with its index in the input set as a tuple using enumerate
        for i in chain.from_iterable(
            combinations(enumerate(inputs), r) for r in combo_range
        ):
            # Yield each combination in the set
            yield i

    def _check_sum_cancellations(self, values_list, equality=False):
        # First, strip any terms with value 0 as they do not contribute to cancellation
        # We do this to keep the number of possible combinations as small as possible
        stripped = [i for i in values_list if abs(i) >= self._zero_tolerance]

        cancellations = []

        if len(stripped) == 0:
            # If the stripped list is empty, there are no non-zero terms
            # We can stop here and return False as there are no possible cancellations
            return cancellations

        # For scaling of tolerance, we want to compare to the largest absolute value of
        # the input values
        max_value = abs(max(stripped, key=abs))

        for i in self._generate_combinations(stripped, equality):
            # Generate combinations will return a list of combinations of input terms
            # each element of the list will be a 2-tuple representing a term in the
            # input list with the first value being the position in the input set and
            # the second being the value

            # Check if the sum of values in the combination are below tolerance
            if abs(sum(j[1] for j in i)) <= self._sum_tol * max_value:
                # If so, record combination as canceling
                cancellations.append(i)

            # Terminate loop if we have reached the max cancellations to collect
            if len(cancellations) >= self._max_cancellations_per_node:
                self._cancellation_tripped = True
                break

        return cancellations

    def _perform_checks(self, node, child_data):
        # Perform checks for problematic expressions
        # First, need to check to see if any child data is a list
        # This indicates a sum expression
        const = True

        for d in child_data:
            # We will check for canceling terms here, rather than the sum itself, to handle special cases
            # We want to look for cases where a sum term results in a value much smaller
            # than the terms of the sum
            # Each element of child_data is a tuple where the 0-th element is the node values
            if isinstance(d[0][0], str):
                # Values may be a list containing a string in some cases (e.g. external functions)
                # Skip if this is the case
                pass
            else:
                for c in self._check_sum_cancellations(d[0]):
                    if node in self.canceling_terms.keys():
                        self.canceling_terms[node].append(c)
                    else:
                        self.canceling_terms[node] = [c]

            # Expression is not constant if any child is not constant
            # Element 1 is a bool indicating if the child node is constant
            if not d[1]:
                const = False

        # Return any problematic terms found
        return const

    def _check_base_type(self, node):
        if isinstance(node, VarData):
            const = node.fixed
        else:
            const = True
        return [value(node)], const, False

    def _check_equality_expression(self, node, child_data):
        # (In)equality expressions are a special case of sum expressions
        # child_data has two elements; 0 is the LHS of the (in)equality and 1 is the RHS
        # Each of these then contains three elements; 0 is a list of values for the sum components,
        # 1 is a bool indicating if the child node term is constant, and 2 is a bool indicating if
        # the child ode is a sum expression.

        # First, to check for cancellation we need to negate one side of the expression
        # mdata will contain the new set of child_data with negated values
        mdata = []
        # child_data[0][0] has the values of the LHS of the (in)equality, and we will negate these
        vals = []
        for j in child_data[0][0]:
            vals.append(-j)
        # Append the negated values along with whether the node is constant to mdata
        mdata.append((vals, child_data[0][1]))
        # child_data[1] is the RHS, so we take this as it appears and append to mdata
        mdata.append(child_data[1])

        # Next, call the method to check the sum expression
        vals, const, _ = self._check_sum_expression(node, mdata)

        # Next, we need to check for canceling terms.
        # child_data[x][2] indicates if a node is a sum expression or not
        if not child_data[0][2] and not child_data[1][2]:
            # If both sides are not sum expressions, we have the form a == b
            # Simple lLinking constraints do not need to be checked
            pass
        # Next, we can ignore any term that has already been flagged as mismatched
        elif node in self.mismatched_terms.keys():
            pass
        # We can also ignore any case where one side of the (in)equality is constant
        # I.e. if either child_node[x][1] is True
        elif any(d[1] for d in mdata):
            pass
        else:
            # Check for cancellation
            # First, collect terms from both sides
            # Note: outer loop comes first in list comprehension
            t = [i for d in mdata for i in d[0]]

            # Then check for cancellations
            for c in self._check_sum_cancellations(t, equality=True):
                if node in self.canceling_terms.keys():
                    self.canceling_terms[node].append(c)
                else:
                    self.canceling_terms[node] = [c]

        return vals, const, False

    def _check_product_expr(self, node, child_data):
        # We expect child_data to be a tuple of len 2 (a*b)
        # If both or neither a and b are sums (xor), handle like other expressions
        if not child_data[0][2] ^ child_data[1][2]:
            return self._check_general_expr(node, child_data)
        else:
            # Here we have the case of a sum and a multiplier
            # For this case, we will make the result look like a sum with the
            # multiplier applied
            # This is important for scaling of the sum terms, as cancellation should be
            # Checked on the scaled value

            # First, check if both terms are constants - if so we can just get the value of
            # this node and move on.
            if child_data[0][1] and child_data[1][1]:
                return self._check_general_expr(node, child_data)

            # The length of the values (child_data[i][0]) of the multiplier will be 1
            # We can just iterate over all terms in both value lists and multiply
            # each pair
            vals = [i * j for i in child_data[0][0] for j in child_data[1][0]]

            # Return the list of values, not constant, is a sum expression (apparent)
            return vals, False, True

    def _check_div_expr(self, node, child_data):
        # We expect child_data to be a tuple of len 2 (a/b)
        # If the numerator is not a sum, handle like general expression
        # child_data[0] is numerator, child_data[1] is denominator
        if not child_data[0][2]:
            return self._check_general_expr(node, child_data)
        else:
            # If the numerator is a sum, we will treat this as if the division
            # were applied to each term in the numerator separately
            # This is important for scaling of the sum terms, as cancellation should be
            # Checked on the scaled value

            # First, check if the numerator is constant - if so we can just get the value of
            # this node and move on.
            if child_data[0][1]:
                return self._check_general_expr(node, child_data)

            # Next, we need to get the value of the denominator
            denom = self._get_value_for_sum_subexpression(child_data[1])

            try:
                vals = [i / denom for i in child_data[0][0]]
            except ZeroDivisionError:
                raise ZeroDivisionError(
                    f"Error in ConstraintTermAnalysisVisitor: found division with denominator of 0 "
                    f"({str(node)})."
                )

            # Return the list of values, not constant, is a sum expression (apparent)
            return vals, False, True

    def _check_negation_expr(self, node, child_data):
        # Here we need to defer checking of cancellations too due to different ways
        # these can appear in an expression.
        # We will simply negate all values on the child node values (child_data[0][0])
        # and pass on the rest.
        vals = [-i for i in child_data[0][0]]

        # Return the list of values, not constant, is a sum expression (apparent)
        return vals, child_data[0][1], child_data[0][2]

    def _check_general_expr(self, node, child_data):
        const = self._perform_checks(node, child_data)

        try:
            # pylint: disable-next=protected-access
            val = node._apply_operation(
                list(map(self._get_value_for_sum_subexpression, child_data))
            )
        except ValueError:
            raise ValueError(
                f"Error in ConstraintTermAnalysisVisitor: error evaluating {str(node)}."
            )
        except ZeroDivisionError:
            raise ZeroDivisionError(
                f"Error in ConstraintTermAnalysisVisitor: found division with denominator of 0 "
                f"({str(node)})."
            )
        except Exception as err:
            # Catch and re-raise any other error that occurs
            _log.exception(f"Unexpected {err=}, {type(err)=}")
            raise

        return [val], const, False

    def _check_other_expression(self, node, child_data):
        const = self._perform_checks(node, child_data)

        # First, need to get value of input terms, which may be sub-expressions
        input_mag = []
        for i in child_data:
            input_mag.append(self._get_value_for_sum_subexpression(i))

        # Next, create a copy of the function with expected magnitudes as inputs
        newfunc = node.create_node_with_local_data(input_mag)

        # Evaluate new function and return the value along with check results
        return [value(newfunc)], const, False

    def _check_ranged_expression(self, node, child_data):
        # child_data should have 3 elements, LHS, middle term, and RHS
        lhs_vals, lhs_const, _ = self._check_equality_expression(node, child_data[:2])
        rhs_vals, rhs_const, _ = self._check_equality_expression(node, child_data[1:])

        # Constant is a bit vague in terms ranged expressions.
        # We will call the term constant only if all parts are constant
        const = lhs_const and rhs_const

        # For values, we need to avoid double counting the middle term
        # Also for sign convention, we will negate the outer terms
        vals = lhs_vals + [-rhs_vals[1]]

        return vals, const, False

    def _check_sum_expression(self, node, child_data):
        # Sum expressions need special handling
        # For sums, collect all child values into a list
        vals = []
        # We will check for cancellation in this node at the next level
        # Pyomo is generally good at simplifying compound sums
        const = True
        # Collect data from child nodes
        for d in child_data:
            vals.append(self._get_value_for_sum_subexpression(d))

            # Expression is not constant if any child is not constant
            # Element 1 is a bool indicating if node is constant
            if not d[1]:
                const = False

        # Check for mismatched terms
        if len(vals) > 1:
            absvals = [abs(v) for v in vals if abs(v) >= self._zero_tolerance]
            if len(absvals) > 0:
                vl = max(absvals)
                vs = min(absvals)
                if vl != vs:
                    if vs == 0:
                        # This branch is reachable only if the user
                        # set _zero_tolerance to 0 or a negative number
                        diff = log10(vl)
                    else:
                        diff = log10(vl / vs)

                    if diff >= self._log_mm_tol:
                        self.mismatched_terms[node] = (vl, vs)

        return vals, const, True

    node_type_method_map = {
        EXPR.EqualityExpression: _check_equality_expression,
        EXPR.InequalityExpression: _check_equality_expression,
        EXPR.RangedExpression: _check_ranged_expression,
        EXPR.SumExpression: _check_sum_expression,
        EXPR.NPV_SumExpression: _check_sum_expression,
        EXPR.ProductExpression: _check_product_expr,
        EXPR.MonomialTermExpression: _check_product_expr,
        EXPR.NPV_ProductExpression: _check_product_expr,
        EXPR.DivisionExpression: _check_div_expr,
        EXPR.NPV_DivisionExpression: _check_div_expr,
        EXPR.PowExpression: _check_general_expr,
        EXPR.NPV_PowExpression: _check_general_expr,
        EXPR.NegationExpression: _check_negation_expr,
        EXPR.NPV_NegationExpression: _check_negation_expr,
        EXPR.AbsExpression: _check_general_expr,
        EXPR.NPV_AbsExpression: _check_general_expr,
        EXPR.UnaryFunctionExpression: _check_general_expr,
        EXPR.NPV_UnaryFunctionExpression: _check_general_expr,
        EXPR.Expr_ifExpression: _check_other_expression,
        EXPR.ExternalFunctionExpression: _check_other_expression,
        EXPR.NPV_ExternalFunctionExpression: _check_other_expression,
        EXPR.LinearExpression: _check_sum_expression,
    }

    def exitNode(self, node, data):
        """
        Method to call when exiting node to check for potential issues.
        """
        # Return [node values], constant (bool), sum expression (bool)
        # first check if the node is a leaf
        nodetype = type(node)

        if nodetype in native_types:
            return [node], True, False

        node_func = self.node_type_method_map.get(nodetype, None)
        if node_func is not None:
            return node_func(self, node, data)

        if not node.is_expression_type():
            # this is a leaf, but not a native type
            if nodetype is _PyomoUnit:
                # Unit have no value, so return 1 as a placeholder
                return [1], True, False

            # Var or Param
            return self._check_base_type(node)
            # might want to add other common types here

        # not a leaf - check if it is a named expression
        if (
            hasattr(node, "is_named_expression_type")
            and node.is_named_expression_type()
        ):
            return self._check_other_expression(node, data)

        raise TypeError(
            f"An unhandled expression node type: {str(nodetype)} was encountered while "
            f"analyzing constraint terms {str(node)}"
        )

    def walk_expression(self, expr):
        """
        Main method to call to walk an expression and return analysis.

        Args:
            expr - expression to be analyzed

        Returns:
            list of values of top-level additive terms
            ComponentSet containing any mismatched terms
            ComponentSet containing any canceling terms
            Bool indicating whether expression is a constant
        """
        # Create new holders for collected terms
        self.canceling_terms = ComponentMap()
        self.mismatched_terms = ComponentMap()

        # Call parent walk_expression method
        vals, const, _ = super().walk_expression(expr)

        # Return results
        return (
            vals,
            self.mismatched_terms,
            self.canceling_terms,
            const,
            self._cancellation_tripped,
        )
