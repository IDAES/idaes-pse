import pyomo.environ as pe
from pyomo.core.expr.calculus.diff_with_pyomo import reverse_sd
from pyomo.common.dependencies import attempt_import

from pyomo.core.base.block import _BlockData
from pyomo.contrib.fbbt.fbbt import fbbt
from .var_utils import (
    get_used_unfixed_variables,
    _apply_var_bounds,
)
from typing import Sequence, Mapping
from pyomo.core.base.var import _GeneralVarData
from pyomo.core.expr.sympy_tools import sympyify_expression, sympy2pyomo_expression
from pyomo.contrib.solver.util import get_objective
from .indices import _VarIndex, _ConIndex


def _simplify(expr):
    cm, sympy_expr = sympyify_expression(expr)
    sympy_expr = sympy_expr.simplify()
    pyomo_expr = sympy2pyomo_expression(sympy_expr, cm)
    return pyomo_expr


def _add_grad_lag_constraints(m: _BlockData) -> _BlockData:
    primal_vars = get_used_unfixed_variables(m)

    m.duals_eq_set = pe.Set()
    m.duals_eq = pe.Var(m.duals_eq_set)

    m.duals_ineq_set = pe.Set()
    m.duals_ineq = pe.Var(m.duals_ineq_set, bounds=(0, None))

    obj = get_objective(m)
    assert obj.sense == pe.minimize
    if obj is None:
        lagrangian = 0
    else:
        lagrangian = obj.expr

    for c in m.component_data_objects(pe.Constraint, active=True, descend_into=True):
        if c.equality:
            key = _ConIndex(c, "eq")
            m.duals_eq_set.add(key)
            lagrangian += m.duals_eq[key] * (c.body - c.upper)
        else:
            if c.upper is not None:
                key = _ConIndex(c, "ub")
                m.duals_ineq_set.add(key)
                lagrangian += m.duals_ineq[key] * (c.body - c.upper)
            if c.lower is not None:
                key = _ConIndex(c, "lb")
                m.duals_ineq_set.add(key)
                lagrangian += m.duals_ineq[key] * (c.lower - c.body)

    for v in primal_vars:
        assert v.is_continuous()
        if v.ub is not None:
            key = _VarIndex(v, "ub")
            m.duals_ineq_set.add(key)
            lagrangian += m.duals_ineq[key] * (v - v.ub)
        if v.lb is not None:
            key = _VarIndex(v, "lb")
            m.duals_ineq_set.add(key)
            lagrangian += m.duals_ineq[key] * (v.lb - v)

    grad_lag = reverse_sd(lagrangian)

    m.grad_lag_set = pe.Set()
    m.grad_lag = pe.Constraint(m.grad_lag_set)
    for v in primal_vars:
        if v in grad_lag and (type(grad_lag[v]) != float or grad_lag[v] != 0):
            key = _VarIndex(v, None)
            m.grad_lag_set.add(key)
            m.grad_lag[key] = _simplify(grad_lag[v]) == 0

    return m


def _introduce_inequality_slacks(m) -> _BlockData:
    m.slacks = pe.Var(m.duals_ineq_set, bounds=(0, None))
    m.ineq_cons_with_slacks = pe.Constraint(m.duals_ineq_set)

    for key in m.duals_ineq_set:
        s = m.slacks[key]
        bnd = key.bound

        if isinstance(key, _ConIndex):
            e = key.con.body
            lb = key.con.lower
            ub = key.con.upper
        else:
            assert isinstance(key, _VarIndex)
            e = key.var
            lb = e.lb
            ub = e.ub

        if bnd == "ub":
            m.ineq_cons_with_slacks[key] = s + e - ub == 0
        else:
            assert bnd == "lb"
            m.ineq_cons_with_slacks[key] = s + lb - e == 0

    for key in m.duals_ineq_set:
        if isinstance(key, _ConIndex):
            key.con.deactivate()
        else:
            key.var.setlb(None)
            key.var.setub(None)
            key.var.domain = pe.Reals

    return m


def _do_fbbt(m, uncertain_params):
    p_bounds = pe.ComponentMap()
    for p in uncertain_params:
        p.unfix()
        p_bounds[p] = (p.lb, p.ub)
    fbbt(m)
    for p in uncertain_params:
        if p.lb > p_bounds[p][0] + 1e-6 or p.ub < p_bounds[p][1] - 1e-6:
            raise RuntimeError(
                "The bounds provided in valid_var_bounds were proven to "
                "be invalid for some values of the uncertain parameters."
            )
        p.fix()


def add_kkt_with_milp_complementarity_conditions(
    m: _BlockData,
    uncertain_params: Sequence[_GeneralVarData],
    valid_var_bounds: Mapping[_GeneralVarData, Sequence[float]],
    default_M=None,
) -> _BlockData:
    for v in uncertain_params:
        v.fix()

    _add_grad_lag_constraints(m)
    obj = get_objective(m)
    obj.deactivate()

    _apply_var_bounds(valid_var_bounds)
    _do_fbbt(m, uncertain_params)

    _introduce_inequality_slacks(m)

    _do_fbbt(m, uncertain_params)

    m.active_indicator = pe.Var(m.duals_ineq_set, domain=pe.Binary)
    m.dual_ineq_0_if_not_active = pe.Constraint(m.duals_ineq_set)
    m.slack_0_if_active = pe.Constraint(m.duals_ineq_set)
    m.dual_M = pe.Param(m.duals_ineq_set, mutable=True)
    m.slack_M = pe.Param(m.duals_ineq_set, mutable=True)
    for key in m.duals_ineq_set:
        if m.duals_ineq[key].ub is None:
            if default_M is None:
                raise RuntimeError(
                    f"could not compute upper bound on multiplier for inequality {key}."
                )
            else:
                dual_M = default_M
        else:
            dual_M = m.duals_ineq[key].ub
        if m.slacks[key].ub is None:
            if default_M is None:
                raise RuntimeError(
                    f"could not compute upper bound on slack for inequality {key}"
                )
            else:
                slack_M = default_M
        else:
            slack_M = m.slacks[key].ub
        m.dual_M[key].value = dual_M
        m.slack_M[key].value = slack_M
        m.dual_ineq_0_if_not_active[key] = (
            m.duals_ineq[key] <= m.active_indicator[key] * m.dual_M[key]
        )
        m.slack_0_if_active[key] = (
            m.slacks[key] <= (1 - m.active_indicator[key]) * m.slack_M[key]
        )

    for v in uncertain_params:
        v.unfix()

    return m
