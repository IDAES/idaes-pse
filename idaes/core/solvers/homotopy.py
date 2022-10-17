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
IDAES Homotopy meta-solver routine.
"""

__author__ = "Andrew Lee"

import logging

from pyomo.environ import Block, SolverFactory, TerminationCondition
from pyomo.core.base.var import _VarData
from pyomo.contrib.parmest.utils.ipopt_solver_wrapper import ipopt_solve_with_stats

from idaes.core.util.model_serializer import to_json, from_json
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def homotopy(
    model,
    variables,
    targets,
    max_solver_iterations=50,
    max_solver_time=10,
    step_init=0.1,
    step_cut=0.5,
    iter_target=4,
    step_accel=0.5,
    max_step=1,
    min_step=0.05,
    max_eval=200,
):
    """
    Homotopy meta-solver routine using Ipopt as the non-linear solver. This
    routine takes a model along with a list of fixed variables in that model
    and a list of target values for those variables. The routine then tries to
    iteratively move the values of the fixed variables to their target values
    using an adaptive step size.

    Args:
        model : model to be solved
        variables : list of Pyomo Var objects to be varied using homotopy.
                    Variables must be fixed.
        targets : list of target values for each variable
        max_solver_iterations : maximum number of solver iterations per
                    homotopy step (default=50)
        max_solver_time : maximum cpu time for the solver per homotopy step
                    (default=10)
        step_init : initial homotopy step size (default=0.1)
        step_cut : factor by which to reduce step size on failed step
                    (default=0.5)
        step_accel : acceleration factor for adjusting step size on successful
                     step (default=0.5)
        iter_target : target number of solver iterations per homotopy step
                    (default=4)
        max_step : maximum homotopy step size (default=1)
        min_step : minimum homotopy step size (default=0.05)
        max_eval : maximum number of homotopy evaluations (both successful and
                   unsuccessful) (default=200)

    Returns:
        Termination Condition : A Pyomo TerminationCondition Enum indicating
            how the meta-solver terminated (see documentation)
        Solver Progress : a fraction indication how far the solver progressed
            from the initial values to the target values
        Number of Iterations : number of homotopy evaluations before solver
            terminated
    """
    eps = 1e-3  # Tolerance for homotopy step convergence to 1

    # Get model logger
    _log = logging.getLogger(__name__)

    # Validate model is an instance of Block
    if not isinstance(model, Block):
        raise TypeError(
            "Model provided was not a valid Pyomo model object "
            "(instance of Block). Please provide a valid model."
        )
    if degrees_of_freedom(model) != 0:
        raise ConfigurationError(
            "Degrees of freedom in model are not equal to zero. Homotopy "
            "should not be used on probelms which are not well-defined."
        )

    # Validate variables and targets
    if len(variables) != len(targets):
        raise ConfigurationError("Number of variables and targets do not match.")
    for i in range(len(variables)):
        v = variables[i]
        t = targets[i]

        if not isinstance(v, _VarData):
            raise TypeError(
                "Variable provided ({}) was not a valid Pyomo Var "
                "component.".format(v)
            )

        # Check that v is part of model
        parent = v.parent_block()
        while parent != model:
            if parent is None:
                raise ConfigurationError("Variable {} is not part of model".format(v))
            parent = parent.parent_block()

        # Check that v is fixed
        if not v.fixed:
            raise ConfigurationError(
                "Homotopy metasolver provided with unfixed variable {}."
                "All variables must be fixed.".format(v.name)
            )

        # Check bounds on v (they don't really matter, but check for sanity)
        if v.ub is not None:
            if v.value > v.ub:
                raise ConfigurationError(
                    "Current value for variable {} is greater than the "
                    "upper bound for that variable. Please correct this "
                    "before continuing.".format(v.name)
                )
            if t > v.ub:
                raise ConfigurationError(
                    "Target value for variable {} is greater than the "
                    "upper bound for that variable. Please correct this "
                    "before continuing.".format(v.name)
                )
        if v.lb is not None:
            if v.value < v.lb:
                raise ConfigurationError(
                    "Current value for variable {} is less than the "
                    "lower bound for that variable. Please correct this "
                    "before continuing.".format(v.name)
                )
            if t < v.lb:
                raise ConfigurationError(
                    "Target value for variable {} is less than the "
                    "lower bound for that variable. Please correct this "
                    "before continuing.".format(v.name)
                )

    # TODO : Should we be more restrictive on these values to avoid users
    # TODO : picking numbers that are less likely to solve (but still valid)?
    # Validate homotopy parameter selections
    if not 0.05 <= step_init <= 0.8:
        raise ConfigurationError(
            "Invalid value for step_init ({}). Must lie "
            "between 0.05 and 0.8.".format(step_init)
        )
    if not 0.1 <= step_cut <= 0.9:
        raise ConfigurationError(
            "Invalid value for step_cut ({}). Must lie "
            "between 0.1 and 0.9.".format(step_cut)
        )
    if step_accel < 0:
        raise ConfigurationError(
            "Invalid value for step_accel ({}). Must be "
            "greater than or equal to 0.".format(step_accel)
        )
    if iter_target < 1:
        raise ConfigurationError(
            "Invalid value for iter_target ({}). Must be "
            "greater than or equal to 1.".format(iter_target)
        )
    if not isinstance(iter_target, int):
        raise ConfigurationError(
            "Invalid value for iter_target ({}). Must be "
            "an an integer.".format(iter_target)
        )
    if not 0.05 <= max_step <= 1:
        raise ConfigurationError(
            "Invalid value for max_step ({}). Must lie "
            "between 0.05 and 1.".format(max_step)
        )
    if not 0.01 <= min_step <= 0.1:
        raise ConfigurationError(
            "Invalid value for min_step ({}). Must lie "
            "between 0.01 and 0.1.".format(min_step)
        )
    if not min_step <= max_step:
        raise ConfigurationError(
            "Invalid argumnets: step_min must be less " "or equal to step_max."
        )
    if not min_step <= step_init <= max_step:
        raise ConfigurationError(
            "Invalid arguments: step_init must lie " "between min_step and max_step."
        )
    if max_eval < 1:
        raise ConfigurationError(
            "Invalid value for max_eval ({}). Must be "
            "greater than or equal to 1.".format(step_accel)
        )
    if not isinstance(max_eval, int):
        raise ConfigurationError(
            "Invalid value for max_eval ({}). Must be "
            "an an integer.".format(iter_target)
        )

    # Create solver object
    solver_obj = SolverFactory("ipopt")

    # Perform initial solve of model to confirm feasible initial solution
    results, solved, sol_iter, sol_time, sol_reg = ipopt_solve_with_stats(
        model, solver_obj, max_solver_iterations, max_solver_time
    )

    if not solved:
        _log.exception("Homotopy Failed - initial solution infeasible.")
        return TerminationCondition.infeasible, 0, 0
    elif sol_reg != "-":
        _log.warning("Homotopy - initial solution converged with regularization.")
        return TerminationCondition.other, 0, 0
    else:
        _log.info("Homotopy - initial point converged")

    # Set up homotopy variables
    # Get initial values and deltas for all variables
    v_init = []
    for i in range(len(variables)):
        v_init.append(variables[i].value)

    n_0 = 0.0  # Homotopy progress variable
    s = step_init  # Set step size to step_init
    iter_count = 0  # Counter for homotopy iterations

    # Save model state to dict
    # TODO : for very large models, it may be necessary to dump this to a file
    current_state = to_json(model, return_dict=True)

    while n_0 < 1.0:
        iter_count += 1  # Increase iter_count regardless of success or failure

        # Calculate next n value given current step size
        if n_0 + s >= 1.0 - eps:
            n_1 = 1.0
        else:
            n_1 = n_0 + s

        _log.info(
            "Homotopy Iteration {}. Next Step: {} (Current: {})".format(
                iter_count, n_1, n_0
            )
        )

        # Update values for all variables using n_1
        for i in range(len(variables)):
            variables[i].fix(targets[i] * n_1 + v_init[i] * (1 - n_1))

        # Solve model at new state
        results, solved, sol_iter, sol_time, sol_reg = ipopt_solve_with_stats(
            model, solver_obj, max_solver_iterations, max_solver_time
        )

        # Check solver output for convergence
        if solved:
            # Step succeeded - accept current state
            current_state = to_json(model, return_dict=True)

            # Update n_0 to accept current step
            n_0 = n_1

            # Check solver iterations and calculate next step size
            s_proposed = s * (1 + step_accel * (iter_target / sol_iter - 1))

            if s_proposed > max_step:
                s = max_step
            elif s_proposed < min_step:
                s = min_step
            else:
                s = s_proposed
        else:
            # Step failed - reload old state
            from_json(model, current_state)

            # Try to cut back step size
            if s > min_step:
                # Step size can be cut
                s = max(min_step, s * step_cut)
            else:
                # Step is already at minimum size, terminate homotopy
                _log.exception(
                    "Homotopy failed - could not converge at minimum step "
                    "size. Current progress is {}".format(n_0)
                )
                return TerminationCondition.minStepLength, n_0, iter_count

        if iter_count >= max_eval:  # Use greater than or equal to to be safe
            _log.exception(
                "Homotopy failed - maximum homotopy iterations "
                "exceeded. Current progress is {}".format(n_0)
            )
            return TerminationCondition.maxEvaluations, n_0, iter_count

    if sol_reg == "-":
        _log.info(
            "Homotopy successful - converged at target values in {} "
            "iterations.".format(iter_count)
        )
        return TerminationCondition.optimal, n_0, iter_count
    else:
        _log.exception(
            "Homotopy failed - converged at target values with "
            "regularization in {} iterations.".format(iter_count)
        )
        return TerminationCondition.other, n_0, iter_count
