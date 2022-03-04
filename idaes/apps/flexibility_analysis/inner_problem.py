import pyomo.environ as pe
from pyomo.core.base.block import _BlockData
from pyomo.core.base.var import _GeneralVarData, ScalarVar
from pyomo.common.dependencies import attempt_import
coramin, coramin_available = attempt_import('coramin', 'coramin is required for flexibility analysis')
from .var_utils import get_all_unfixed_variables, BoundsManager, _apply_var_bounds
from .indices import _ConIndex, _VarIndex
from typing import MutableMapping, Tuple, Optional, Mapping, Sequence
import math
from pyomo.contrib.fbbt.fbbt import compute_bounds_on_expr
from pyomo.core.expr.calculus.diff_with_pyomo import reverse_sd


def _get_bounds_on_max_constraint_violation(
    m: _BlockData, valid_var_bounds: Mapping[_GeneralVarData, Sequence[float]]
):
    bounds_manager = BoundsManager(m)
    bounds_manager.save_bounds()

    _apply_var_bounds(valid_var_bounds)

    min_constraint_violation = math.inf
    max_constraint_violation = -math.inf
    m.max_constraint_violation.fix(0)
    for c in m.ineq_violation_cons.values():
        exprs = list()
        assert c.lower is None
        assert c.upper is not None
        e = c.body - c.upper
        ders = reverse_sd(e)
        assert ders[m.max_constraint_violation] == -1
        _lb, _ub = compute_bounds_on_expr(e)
        if _lb < min_constraint_violation:
            min_constraint_violation = _lb
        if _ub > max_constraint_violation:
            max_constraint_violation = _ub

    bounds_manager.pop_bounds()
    m.max_constraint_violation.unfix()

    return min_constraint_violation, max_constraint_violation


def _get_constraint_violation_bounds(
    m: _BlockData, valid_var_bounds: Mapping[_GeneralVarData, Sequence[float]]
) -> MutableMapping[_GeneralVarData, Tuple[float, float]]:
    bounds_manager = BoundsManager(m)
    bounds_manager.save_bounds()

    _apply_var_bounds(valid_var_bounds)

    constraint_violation_bounds = pe.ComponentMap()
    for key in m.ineq_violation_set:
        v = m.constraint_violation[key]
        v.fix(0)

        c = m.ineq_violation_cons[key]
        assert c.equality
        e = c.body - c.upper
        ders = reverse_sd(e)
        assert ders[v] == -1
        _lb, _ub = compute_bounds_on_expr(e)
        constraint_violation_bounds[v] = (_lb, _ub)
        v.unfix()

    bounds_manager.pop_bounds()

    return constraint_violation_bounds


def _build_inner_problem(
    m: _BlockData,
    enforce_equalities: bool,
    unique_constraint_violations: bool,
    valid_var_bounds: Optional[MutableMapping[_GeneralVarData, Tuple[float, float]]],
):
    """
    If enfoce equalities is True and unique_constraint_violations is False, then this function converts

        min f(x)
        s.t.
            c(x) = 0
            g(x) <= 0

    to

        min u
        s.t.
            c(x) = 0
            g(x) <= u

    If enfoce equalities is False and unique_constraint_violations is False, then this function converts

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

    If enfoce equalities is True and unique_constraint_violations is True, then this function converts

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

    This function will also modify valid_var_bounds to include any new variables
    """
    obj = coramin.utils.get_objective(m)
    if obj is not None:
        obj.deactivate()

    for v in m.unc_param_vars.values():
        v.fix()
    original_vars = list(get_all_unfixed_variables(m))
    for v in m.unc_param_vars.values():
        v.unfix()

    m.max_constraint_violation = ScalarVar()
    m.min_constraint_violation_obj = pe.Objective(expr=m.max_constraint_violation)

    m.ineq_violation_set = pe.Set()
    m.ineq_violation_cons = pe.Constraint(m.ineq_violation_set)

    if unique_constraint_violations:
        m.constraint_violation = pe.Var(m.ineq_violation_set)

    for c in list(
        m.component_data_objects(pe.Constraint, descend_into=True, active=True)
    ):
        if c.equality and enforce_equalities:
            continue
        if c.lower is not None:
            key = _ConIndex(c, "lb")
            m.ineq_violation_set.add(key)
            if unique_constraint_violations:
                m.ineq_violation_cons[key] = (
                    c.lower - c.body - m.constraint_violation[key],
                    0,
                )
            else:
                m.ineq_violation_cons[key] = (
                    None,
                    c.lower - c.body - m.max_constraint_violation,
                    0,
                )
        if c.upper is not None:
            key = _ConIndex(c, "ub")
            m.ineq_violation_set.add(key)
            if unique_constraint_violations:
                m.ineq_violation_cons[key] = (
                    c.body - c.upper - m.constraint_violation[key],
                    0,
                )
            else:
                m.ineq_violation_cons[key] = (
                    None,
                    c.body - c.upper - m.max_constraint_violation,
                    0,
                )

    for v in original_vars:
        if v.is_integer():
            raise ValueError("Original problem must be continuous")
        if v.lb is not None:
            key = _VarIndex(v, "lb")
            m.ineq_violation_set.add(key)
            if unique_constraint_violations:
                m.ineq_violation_cons[key] = (v.lb - v - m.constraint_violation[key], 0)
            else:
                m.ineq_violation_cons[key] = (
                    None,
                    v.lb - v - m.max_constraint_violation,
                    0,
                )
        if v.ub is not None:
            key = _VarIndex(v, "ub")
            m.ineq_violation_set.add(key)
            if unique_constraint_violations:
                m.ineq_violation_cons[key] = (v - v.ub - m.constraint_violation[key], 0)
            else:
                m.ineq_violation_cons[key] = (
                    None,
                    v - v.ub - m.max_constraint_violation,
                    0,
                )

    for key in m.ineq_violation_set:
        if isinstance(key, _ConIndex):
            key.con.deactivate()
        else:
            key.var.setlb(None)
            key.var.setub(None)
            key.var.domain = pe.Reals

    if unique_constraint_violations:
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
        constraint_violation_bounds = _get_constraint_violation_bounds(
            m, valid_var_bounds
        )
        m.u_hat_cons = pe.ConstraintList()
        for key in m.ineq_violation_set:
            violation_var = m.constraint_violation[key]
            viol_lb, viol_ub = constraint_violation_bounds[violation_var]
            y_i = m.max_violation_selector[key]
            m.u_hat_cons.add(m.u_hat[key] <= viol_ub * y_i)
            m.u_hat_cons.add(m.u_hat[key] >= viol_lb * y_i)
            m.u_hat_cons.add(m.u_hat[key] <= violation_var + viol_lb * y_i - viol_lb)
            m.u_hat_cons.add(m.u_hat[key] >= viol_ub * y_i + violation_var - viol_ub)
        valid_var_bounds.update(constraint_violation_bounds)
        valid_var_bounds[m.max_constraint_violation] = (
            min(i[0] for i in constraint_violation_bounds.values()),
            max(i[1] for i in constraint_violation_bounds.values()),
        )
        for key in m.ineq_violation_set:
            valid_var_bounds[m.max_violation_selector[key]] = (0, 1)
            valid_var_bounds[m.u_hat[key]] = (
                min(0.0, constraint_violation_bounds[m.constraint_violation[key]][0]),
                max(0.0, constraint_violation_bounds[m.constraint_violation[key]][1]),
            )
    else:
        if valid_var_bounds is not None:
            valid_var_bounds[
                m.max_constraint_violation
            ] = _get_bounds_on_max_constraint_violation(
                m=m, valid_var_bounds=valid_var_bounds
            )
