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
This module contains a tool for determining whether a model is ill-conditioned.
"""

__author__ = "Alexander Dowling, Douglas Allan, Andrew Lee, Robby Parker, Ben Knueven"


from pyomo.environ import (
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Expression,
    Objective,
    RangeSet,
    SolverFactory,
    value,
    Var,
)

from idaes.core.scaling.util import (
    get_jacobian,
)
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def compute_ill_conditioning_certificate(
    model,
    target_feasibility_tol: float = 1e-06,
    ratio_cutoff: float = 1e-04,
    direction: str = "row",
):
    """
    Finds constraints (rows) or variables (columns) in the model Jacobian that
    may be contributing to ill conditioning.

    This method is based on work published in:

    Klotz, E., Identification, Assessment, and Correction of Ill-Conditioning and
    Numerical Instability in Linear and Integer Programs, Informs 2014, pgs. 54-108
    https://pubsonline.informs.org/doi/epdf/10.1287/educ.2014.0130

    Args:
        model: model to be analysed
        target_feasibility_tol: target tolerance for solving ill conditioning problem
        ratio_cutoff: cut-off for reporting ill conditioning
        direction: 'row' (default, constraints) or 'column' (variables)

    Returns:
        list of strings reporting ill-conditioned variables/constraints and their
        associated y values
    """
    _log.warning(
        "Ill conditioning checks are a beta capability. Please be aware that "
        "the name, location, and API for this may change in future releases."
    )
    # Thanks to B. Knueven for this implementation

    if direction not in ["row", "column"]:
        raise ValueError(
            f"Unrecognised value for direction ({direction}). "
            "Must be 'row' or 'column'."
        )

    jac, nlp = get_jacobian(model)

    inverse_target_kappa = 1e-16 / target_feasibility_tol

    # Set up the components we will analyze, either row or column
    if direction == "row":
        components = nlp.get_pyomo_constraints()
        components_set = RangeSet(0, len(components) - 1)
        results_set = RangeSet(0, nlp.n_primals() - 1)
        jac = jac.transpose().tocsr()
    else:  # direction == "column"
        components = nlp.get_pyomo_variables()
        components_set = RangeSet(0, len(components) - 1)
        results_set = RangeSet(0, nlp.n_constraints() - 1)
        jac = jac.tocsr()

    # Build test problem
    inf_prob = ConcreteModel()

    inf_prob.y_pos = Var(components_set, bounds=(0, None))
    inf_prob.y_neg = Var(components_set, bounds=(0, None))
    inf_prob.y = Expression(components_set, rule=lambda m, i: m.y_pos[i] - m.y_neg[i])

    inf_prob.res_pos = Var(results_set, bounds=(0, None))
    inf_prob.res_neg = Var(results_set, bounds=(0, None))
    inf_prob.res = Expression(
        results_set, rule=lambda m, i: m.res_pos[i] - m.res_neg[i]
    )

    def b_rule(b, i):
        lhs = 0.0

        row = jac.getrow(i)
        for j, val in zip(row.indices, row.data):
            lhs += val * b.y[j]

        return lhs == b.res[i]

    inf_prob.by = Constraint(results_set, rule=b_rule)

    # Normalization of y
    inf_prob.normalize = Constraint(
        expr=1 == sum(inf_prob.y_pos.values()) - sum(inf_prob.y_neg.values())
    )

    inf_prob.y_norm = Var()
    inf_prob.y_norm_constr = Constraint(
        expr=inf_prob.y_norm
        == sum(inf_prob.y_pos.values()) + sum(inf_prob.y_neg.values())
    )

    inf_prob.res_norm = Var()
    inf_prob.res_norm_constr = Constraint(
        expr=inf_prob.res_norm
        == sum(inf_prob.res_pos.values()) + sum(inf_prob.res_neg.values())
    )

    # Objective -- minimize residual
    inf_prob.min_res = Objective(expr=inf_prob.res_norm)

    solver = SolverFactory("cbc")  # TODO: Consider making this an option

    # tighten tolerances  # TODO: If solver is an option, need to allow user options
    solver.options["primalT"] = target_feasibility_tol * 1e-1
    solver.options["dualT"] = target_feasibility_tol * 1e-1

    results = solver.solve(inf_prob, tee=False)
    if not check_optimal_termination(results):
        # TODO: maybe we should tighten tolerances first?
        raise RuntimeError("Ill conditioning diagnostic problem infeasible")

    result_norm = inf_prob.res_norm.value
    if result_norm < 0.0:
        # TODO: try again with tighter tolerances?
        raise RuntimeError(
            "Ill conditioning diagnostic problem has numerically troublesome solution"
        )
    if result_norm >= inverse_target_kappa:
        return []

    # find an equivalent solution which minimizes y_norm
    inf_prob.min_res.deactivate()
    inf_prob.res_norm.fix()

    inf_prob.min_y = Objective(expr=inf_prob.y_norm)

    # if this problem is numerically infeasible, we can still report something to the user
    results = solver.solve(inf_prob, tee=False, load_solutions=False)
    if check_optimal_termination(results):
        inf_prob.solutions.load_from(results)

    ill_cond = []
    slist = sorted(
        inf_prob.y, key=lambda dict_key: abs(value(inf_prob.y[dict_key])), reverse=True
    )
    cutoff = None
    for i in slist:
        if cutoff is None:
            cutoff = abs(value(inf_prob.y[i])) * ratio_cutoff
        val = value(inf_prob.y[i])
        if abs(val) < cutoff:
            break
        ill_cond.append((components[i], val))

    return ill_cond
