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
This module contains an expression walker which looks for potential numerical evaluation errors.
"""

__author__ = "Alexander Dowling, Douglas Allan, Andrew Lee, Robby Parker, Ben Knueven"

from math import inf, isfinite
from typing import List

from pyomo.environ import (
    Binary,
    Integers,
)
from pyomo.core.expr.numeric_expr import (
    DivisionExpression,
    NPV_DivisionExpression,
    PowExpression,
    NPV_PowExpression,
    UnaryFunctionExpression,
    NPV_UnaryFunctionExpression,
    NumericExpression,
)
from pyomo.core.base.var import VarData
from pyomo.repn.standard_repn import (  # pylint: disable=no-name-in-module
    generate_standard_repn,
)
from pyomo.common.collections import ComponentSet
from pyomo.common.config import (
    ConfigDict,
)
from pyomo.core.expr.visitor import StreamBasedExpressionVisitor
from pyomo.contrib.fbbt.fbbt import compute_bounds_on_expr


def _get_bounds_with_inf(node: NumericExpression):
    lb, ub = compute_bounds_on_expr(node)
    if lb is None:
        lb = -inf
    if ub is None:
        ub = inf
    return lb, ub


def _check_eval_error_division(
    node: NumericExpression, warn_list: List[str], config: ConfigDict
):
    lb, ub = _get_bounds_with_inf(node.args[1])
    if (config.warn_for_evaluation_error_at_bounds and (lb <= 0 <= ub)) or (
        lb < 0 < ub
    ):
        msg = f"Potential division by 0 in {node}; Denominator bounds are ({lb}, {ub})"
        warn_list.append(msg)


def _check_eval_error_pow(
    node: NumericExpression, warn_list: List[str], config: ConfigDict
):
    arg1, arg2 = node.args
    lb1, ub1 = _get_bounds_with_inf(arg1)
    lb2, ub2 = _get_bounds_with_inf(arg2)

    integer_domains = ComponentSet([Binary, Integers])

    integer_exponent = False
    # if the exponent is an integer, there should not be any evaluation errors
    if isinstance(arg2, VarData) and arg2.domain in integer_domains:
        # The exponent is an integer variable
        # check if the base can be zero
        integer_exponent = True
    if lb2 == ub2 and lb2 == round(lb2):
        # The exponent is fixed to an integer
        integer_exponent = True
    repn = generate_standard_repn(arg2, quadratic=True)
    if (
        repn.nonlinear_expr is None
        and repn.constant == round(repn.constant)
        and all(i.domain in integer_domains for i in repn.linear_vars)
        and all(i[0].domain in integer_domains for i in repn.quadratic_vars)
        and all(i[1].domain in integer_domains for i in repn.quadratic_vars)
        and all(i == round(i) for i in repn.linear_coefs)
        and all(i == round(i) for i in repn.quadratic_coefs)
    ):
        # The exponent is a linear or quadratic expression containing
        # only integer variables with integer coefficients
        integer_exponent = True

    if integer_exponent and (
        (lb1 > 0 or ub1 < 0)
        or (not config.warn_for_evaluation_error_at_bounds and (lb1 >= 0 or ub1 <= 0))
    ):
        # life is good; the exponent is an integer and the base is nonzero
        return None
    elif integer_exponent and lb2 >= 0:
        # life is good; the exponent is a nonnegative integer
        return None

    # if the base is positive, there should not be any evaluation errors
    if lb1 > 0 or (not config.warn_for_evaluation_error_at_bounds and lb1 >= 0):
        return None
    if lb1 >= 0 and lb2 >= 0:
        return None

    msg = f"Potential evaluation error in {node}; "
    msg += f"base bounds are ({lb1}, {ub1}); "
    msg += f"exponent bounds are ({lb2}, {ub2})"
    warn_list.append(msg)


def _check_eval_error_log(
    node: NumericExpression, warn_list: List[str], config: ConfigDict
):
    lb, ub = _get_bounds_with_inf(node.args[0])
    if (config.warn_for_evaluation_error_at_bounds and lb <= 0) or lb < 0:
        msg = f"Potential log of a non-positive number in {node}; Argument bounds are ({lb}, {ub})"
        warn_list.append(msg)


def _check_eval_error_tan(
    node: NumericExpression, warn_list: List[str], config: ConfigDict
):
    lb, ub = _get_bounds_with_inf(node)
    if not (isfinite(lb) and isfinite(ub)):
        msg = f"{node} may evaluate to -inf or inf; Argument bounds are {_get_bounds_with_inf(node.args[0])}"
        warn_list.append(msg)


def _check_eval_error_asin(
    node: NumericExpression, warn_list: List[str], config: ConfigDict
):
    lb, ub = _get_bounds_with_inf(node.args[0])
    if lb < -1 or ub > 1:
        msg = f"Potential evaluation of asin outside [-1, 1] in {node}; Argument bounds are ({lb}, {ub})"
        warn_list.append(msg)


def _check_eval_error_acos(
    node: NumericExpression, warn_list: List[str], config: ConfigDict
):
    lb, ub = _get_bounds_with_inf(node.args[0])
    if lb < -1 or ub > 1:
        msg = f"Potential evaluation of acos outside [-1, 1] in {node}; Argument bounds are ({lb}, {ub})"
        warn_list.append(msg)


def _check_eval_error_sqrt(
    node: NumericExpression, warn_list: List[str], config: ConfigDict
):
    lb, ub = _get_bounds_with_inf(node.args[0])
    if lb < 0:
        msg = f"Potential square root of a negative number in {node}; Argument bounds are ({lb}, {ub})"
        warn_list.append(msg)


_unary_eval_err_handler = dict()
_unary_eval_err_handler["log"] = _check_eval_error_log
_unary_eval_err_handler["log10"] = _check_eval_error_log
_unary_eval_err_handler["tan"] = _check_eval_error_tan
_unary_eval_err_handler["asin"] = _check_eval_error_asin
_unary_eval_err_handler["acos"] = _check_eval_error_acos
_unary_eval_err_handler["sqrt"] = _check_eval_error_sqrt


def _check_eval_error_unary(
    node: NumericExpression, warn_list: List[str], config: ConfigDict
):
    if node.getname() in _unary_eval_err_handler:
        _unary_eval_err_handler[node.getname()](node, warn_list, config)


_eval_err_handler = dict()
_eval_err_handler[DivisionExpression] = _check_eval_error_division
_eval_err_handler[NPV_DivisionExpression] = _check_eval_error_division
_eval_err_handler[PowExpression] = _check_eval_error_pow
_eval_err_handler[NPV_PowExpression] = _check_eval_error_pow
_eval_err_handler[UnaryFunctionExpression] = _check_eval_error_unary
_eval_err_handler[NPV_UnaryFunctionExpression] = _check_eval_error_unary


class EvalErrorWalker(StreamBasedExpressionVisitor):
    """
    Expression walker that looks for potential numerical evaluation errors.
    """

    def __init__(self, config: ConfigDict):
        super().__init__()
        self._warn_list = list()
        self._config = config

    def exitNode(self, node, data):
        """
        callback to be called as the visitor moves from the leaf
        nodes back to the root node.

        Args:
            node: a pyomo expression node
            data: not used in this walker
        """
        if type(node) in _eval_err_handler:
            _eval_err_handler[type(node)](node, self._warn_list, self._config)
        return self._warn_list
