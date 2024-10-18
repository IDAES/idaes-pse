#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
This module contains functions for formulating the inner problem of the 
flexibility test problem
"""
from typing import MutableMapping, Tuple, Optional, Mapping, Union
from pyomo.contrib.fbbt.fbbt import compute_bounds_on_expr
import pyomo.environ as pe
from pyomo.core.base.block import _BlockData
from pyomo.core.base.var import _GeneralVarData, ScalarVar
from pyomo.core.expr.numeric_expr import ExpressionBase
from pyomo.contrib.solver.util import get_objective
from .var_utils import (
    BoundsManager,
    _apply_var_bounds,
    get_used_unfixed_variables,
    _remove_var_bounds,
)
from .indices import _ConIndex, _VarIndex


def _get_g_bounds(
    m: _BlockData, original_vars, valid_var_bounds: Mapping
) -> MutableMapping[Union[_ConIndex, _VarIndex], Tuple[float, float]]:

    key_list = list()
    g_list = list()

    for c in list(
        m.component_data_objects(pe.Constraint, descend_into=True, active=True)
    ):
        if c.lower is not None:
            key = _ConIndex(c, "lb")
            g = c.lower - c.body
            key_list.append(key)
            g_list.append(g)
        if c.upper is not None:
            key = _ConIndex(c, "ub")
            g = c.body - c.upper
            key_list.append(key)
            g_list.append(g)

    for v in original_vars:
        if v.is_integer():
            raise ValueError("Original problem must be continuous")
        if v.lb is not None:
            key = _VarIndex(v, "lb")
            g = v.lb - v
            key_list.append(key)
            g_list.append(g)
        if v.ub is not None:
            key = _VarIndex(v, "ub")
            g = v - v.ub
            key_list.append(key)
            g_list.append(g)

    bounds_manager = BoundsManager(m)
    bounds_manager.save_bounds()

    _remove_var_bounds(m)
    _apply_var_bounds(valid_var_bounds)

    g_bounds = dict()
    for key, g in zip(key_list, g_list):
        g_bounds[key] = compute_bounds_on_expr(g)

    bounds_manager.pop_bounds()
    return g_bounds


def _add_total_violation_disjunctions(
    m: _BlockData, g: ExpressionBase, key, g_bounds: MutableMapping
):
    m.zero_violation_cons[key] = (
        None,
        (
            m.constraint_violation[key]
            - (1 - m.zero_violation[key]) * m.violation_disjunction_BigM[key]
        ),
        0,
    )
    m.nonzero_violation_cons[key] = (
        None,
        (
            m.constraint_violation[key]
            - g
            - (1 - m.nonzero_violation[key]) * m.violation_disjunction_BigM[key]
        ),
        0,
    )
    m.violation_disjunction_cons[key] = (
        m.zero_violation[key] + m.nonzero_violation[key],
        1,
    )
    lb, ub = g_bounds[key]
    m.violation_disjunction_BigM[key].value = max(abs(lb), abs(ub))


def _process_constraint(
    m: _BlockData,
    g: ExpressionBase,
    key: Union[_ConIndex, _VarIndex],
    unique_constraint_violations: bool,
    total_violation: bool,
    total_violation_disjunctions: bool,
    g_bounds: MutableMapping,
):
    m.ineq_violation_set.add(key)
    if total_violation:
        m.constraint_violation[key].setlb(0)
        if total_violation_disjunctions:
            _add_total_violation_disjunctions(m, g, key, g_bounds)
        else:
            m.ineq_violation_cons[key] = (
                None,
                g - m.constraint_violation[key],
                0,
            )
    elif unique_constraint_violations:
        m.ineq_violation_cons[key] = (
            g - m.constraint_violation[key],
            0,
        )
    else:
        m.ineq_violation_cons[key] = (
            None,
            g - m.max_constraint_violation,
            0,
        )


def _build_inner_problem(
    m: _BlockData,
    enforce_equalities: bool,
    unique_constraint_violations: bool,
    valid_var_bounds: Optional[MutableMapping[_GeneralVarData, Tuple[float, float]]],
    total_violation: bool = False,
    total_violation_disjunctions: bool = False,
):
    """
    If enfoce equalities is True and unique_constraint_violations is False, then this
    function converts

        min f(x)
        s.t.
            c(x) = 0
            g(x) <= 0

    to

        min u
        s.t.
            c(x) = 0
            g(x) <= u

    If enfoce equalities is False and unique_constraint_violations is False, then this
    function converts

        min f(x)
        s.t.
            c(x) = 0
            g(x) <= 0

    to

        min u
        s.t.
            c(x) <= u
            -c(x) <= u
            g(x) <= u

    If enfoce equalities is True and unique_constraint_violations is True, then this
    function converts

        min f(x)
        s.t.
            c(x) = 0
            g(x) <= 0

    to

        min u
        s.t.
            c(x) = 0
            g_i(x) == u_i
            u = sum(u_i * y_i)
            sum(y_i) = 1
    Of course, the nonlinear constraint u = sum(u_i * y_i) gets reformulated.

    If total_violation is True, then unique_constraint_violations is ignored and
    this function converts

        min f(x)
        s.t.
            c(x) = 0
            g(x) <= 0

    to (if enforce_equalities is True)

        min sum_{i} u_i
        s.t.
            c(x) = 0
            g_i(x) <= u_i
            u_i >= 0

    or to (if enforce_equalities is False)

        min sum_{i} u_i + sum_{j} u_j + sum_{k} u_k
        s.t.
            c_j(x) <= u_j
            -c_k(x) <= u_k
            g_i(x) <= u_i
            u_i >= 0
            u_j >= 0
            u_k >= 0

    or to (if total_violation_disjunctions is True)

        max sum_{i} u_i
        s.t.
            c(x) = 0
            u_i >= 0
            [u_i <= 0] v [u_i <= g_i(x)]

    This function will also modify valid_var_bounds to include any new variables
    """
    if total_violation_disjunctions:
        assert total_violation

    obj = get_objective(m)
    if obj is not None:
        obj.deactivate()

    for v in m.unc_param_vars.values():
        v.fix()
    original_vars = list(get_used_unfixed_variables(m))
    for v in m.unc_param_vars.values():
        v.unfix()

    if valid_var_bounds is None:
        g_bounds = dict()
    else:
        g_bounds = _get_g_bounds(m, original_vars, valid_var_bounds)

    m.ineq_violation_set = pe.Set()
    if not total_violation_disjunctions:
        m.ineq_violation_cons = pe.Constraint(m.ineq_violation_set)

    if not total_violation:
        m.max_constraint_violation = ScalarVar()

    if unique_constraint_violations or total_violation:
        m.constraint_violation = pe.Var(m.ineq_violation_set)

    if total_violation_disjunctions:
        m.zero_violation = pe.Var(m.ineq_violation_set, domain=pe.Binary)
        m.nonzero_violation = pe.Var(m.ineq_violation_set, domain=pe.Binary)
        m.zero_violation_cons = pe.Constraint(m.ineq_violation_set)
        m.nonzero_violation_cons = pe.Constraint(m.ineq_violation_set)
        m.violation_disjunction_cons = pe.Constraint(m.ineq_violation_set)
        m.violation_disjunction_BigM = pe.Param(m.ineq_violation_set, mutable=True)

    for c in list(
        m.component_data_objects(pe.Constraint, descend_into=True, active=True)
    ):
        if c.equality and enforce_equalities:
            continue
        if c.lower is not None:
            key = _ConIndex(c, "lb")
            g = c.lower - c.body
            _process_constraint(
                m,
                g,
                key,
                unique_constraint_violations,
                total_violation,
                total_violation_disjunctions,
                g_bounds,
            )
        if c.upper is not None:
            key = _ConIndex(c, "ub")
            g = c.body - c.upper
            _process_constraint(
                m,
                g,
                key,
                unique_constraint_violations,
                total_violation,
                total_violation_disjunctions,
                g_bounds,
            )

    for v in original_vars:
        if v.is_integer():
            raise ValueError("Original problem must be continuous")
        if v.lb is not None:
            key = _VarIndex(v, "lb")
            g = v.lb - v
            _process_constraint(
                m,
                g,
                key,
                unique_constraint_violations,
                total_violation,
                total_violation_disjunctions,
                g_bounds,
            )
        if v.ub is not None:
            key = _VarIndex(v, "ub")
            g = v - v.ub
            _process_constraint(
                m,
                g,
                key,
                unique_constraint_violations,
                total_violation,
                total_violation_disjunctions,
                g_bounds,
            )

    if total_violation:
        m.total_constraint_violation_obj = pe.Objective(
            expr=sum(m.constraint_violation.values())
        )
    else:
        m.min_constraint_violation_obj = pe.Objective(expr=m.max_constraint_violation)

    for key in m.ineq_violation_set:
        if isinstance(key, _ConIndex):
            key.con.deactivate()
        else:
            key.var.setlb(None)
            key.var.setub(None)
            key.var.domain = pe.Reals

    if (
        total_violation or unique_constraint_violations
    ) and valid_var_bounds is not None:
        for key in m.ineq_violation_set:
            lb, ub = g_bounds[key]
            v = m.constraint_violation[key]
            valid_var_bounds[v] = (min(lb, 0), max(ub, 0))

    if unique_constraint_violations and not total_violation:
        # max_constraint_violation = sum(constraint_violation[i] * y[i])
        # sum(y[i]) == 1
        # reformulate as
        # max_constraint_violation = sum(u_hat[i])
        # u_hat[i] = constraint_violation[i] * y[i]
        # and use mccormick for the last constraint

        m.max_violation_selector = pe.Var(
            m.ineq_violation_set, domain=pe.Binary
        )  # y[i]
        m.one_max_violation = pe.Constraint(
            expr=sum(m.max_violation_selector.values()) == 1
        )
        m.u_hat = pe.Var(m.ineq_violation_set)
        m.max_violation_sum = pe.Constraint(
            expr=m.max_constraint_violation == sum(m.u_hat.values())
        )
        m.u_hat_cons = pe.ConstraintList()
        for key in m.ineq_violation_set:
            violation_var = m.constraint_violation[key]
            viol_lb, viol_ub = g_bounds[key]
            valid_var_bounds[m.u_hat[key]] = (min(viol_lb, 0), max(viol_ub, 0))
            valid_var_bounds[m.max_violation_selector[key]] = (0, 1)
            y_i = m.max_violation_selector[key]
            m.u_hat_cons.add(m.u_hat[key] <= viol_ub * y_i)
            m.u_hat_cons.add(m.u_hat[key] >= viol_lb * y_i)
            m.u_hat_cons.add(m.u_hat[key] <= violation_var + viol_lb * y_i - viol_lb)
            m.u_hat_cons.add(m.u_hat[key] >= viol_ub * y_i + violation_var - viol_ub)

    if valid_var_bounds is not None and not total_violation:
        valid_var_bounds[m.max_constraint_violation] = (
            min(i[0] for i in g_bounds.values()),
            max(i[1] for i in g_bounds.values()),
        )
