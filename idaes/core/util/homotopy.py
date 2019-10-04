##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
IDAES Homotopy meta-solver routine.
"""

__author__ = "Andrew Lee"

import logging

from pyomo.environ import (Block,
                           SolverFactory,
                           SolverStatus,
                           TerminationCondition)

from idaes.core.util.model_serializer import to_json, from_json
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import ConfigurationError


def homotopy(model, variables, targets, solver='ipopt', solver_options={},
             step_init=0.1, step_cut = 0.5, iter_target = 10, step_accel = 1,
             max_step = 1, min_step = 0.05, max_iter = 200):
    """
    Homotopy meta-solver routine.

    Args:
        model : model to be solved
        variables : list of Pyomo Var objects to be varied using homotopy. "
                    "Variables must be fixed.
        targets : list of target values for each variable
        solver : solver to use to converge model at each step.
        solver_options : dict of option arguments to be passed to solver object
        step_init : initial homotopy step size
        step_cut : factor to reduce step size by on failed step
        step_accel : acceleration factor for adjusting step size on successful
                     step
        iter_target : target number of solver iterations per homotopy step
        max_step : maximum homotopy step size
        min_step : minimum homotopy step size
        max_iter : maximum number of homotopy iterations (but successful and
                   unsuccessful)

    Returns:
        A ComponentSet including all Block components in block (including block
        itself)
    """
    # Add a setting to display debugging output - don't expose to the user
    DEBUG = True
    eps = 1e-6  # Tolerance for homotopy step convergence to 1

    # Get model logger
    _log = logging.getLogger(__name__)

    # Validate model is an instance of Block
    if not isinstance(model, Block):
        raise TypeError("Model provided was not a valid Pyomo model object "
                        "(instance of Block). Please provide a valid model.")
    if degrees_of_freedom(model) != 0:
        raise ConfigurationError(
                "Degrees of freedom in model are not equal to zero. Homotopy "
                "should not be used on probelms which are not well-defined.")

    # Validate variables and targets
    if len(variables) != len(targets):
        raise ConfigurationError(
                "Number of variables and targets do not match.")
    for i in range(len(variables)):
        v = variables[i]
        t = targets[i]
        if not v.fixed:
            raise ConfigurationError(
                    "Homotopy metasolver provided with unfixed variable {}."
                    "All variables must be fixed.".format(v.name))
        if v.ub is not None:
            if v.value > v.ub:
                raise ConfigurationError(
                        "Current value for variable {} is greater than the "
                        "upper bound for that variable. Please correct this "
                        "before continuing.".format(v.name))
            if t > v.ub:
                raise ConfigurationError(
                        "Target value for variable {} is greater than the "
                        "upper bound for that variable. Please correct this "
                        "before continuing.".format(v.name))
        if v.lb is not None:
            if v.value < v.lb:
                raise ConfigurationError(
                        "Current value for variable {} is less than the "
                        "lower bound for that variable. Please correct this "
                        "before continuing.".format(v.name))
            if t < v.lb:
                raise ConfigurationError(
                        "Target value for variable {} is less than the "
                        "lower bound for that variable. Please correct this "
                        "before continuing.".format(v.name))

    # TODO : Should we be more restrictive on these values to avoid users
    # TODO : picking numbers that are less likely to solve (but still valid)?
    # Validate homotopy parameter selections
    if not 0 < step_init <= 1:
        raise ConfigurationError("Invalid value for step_init ({}). Must lie "
                                 "between 0 and 1.".format(step_init))
    if not 0.1 <= step_cut <= 0.9:
        raise ConfigurationError("Invalid value for step_cut ({}). Must lie "
                                 "between 0.1 and 0.9.".format(step_cut))
    if step_accel < 1:
        raise ConfigurationError("Invalid value for step_accel ({}). Must be "
                                 "greater than or equal to 1."
                                 .format(step_accel))
    if iter_target < 1:
        raise ConfigurationError("Invalid value for iter_target ({}). Must be "
                                 "greater than or equal to 1."
                                 .format(iter_target))
    if not 0.05 < max_step <= 1:
        raise ConfigurationError("Invalid value for max_step ({}). Must lie "
                                 "between 0.05 and 1."
                                 .format(max_step))
    if not 0.01 < min_step <= 0.1:
        raise ConfigurationError("Invalid value for min_step ({}). Must lie "
                                 "between 0.01 and 0.1."
                                 .format(min_step))
    if not min_step <= max_step:
        raise ConfigurationError("Invalid argumnets: step_min must be less "
                                 "or equal to step_max.")
    if not min_step <= step_init <= max_step:
        raise ConfigurationError("Invalid arguments: step_init must lie "
                                 "between min_step and max_step.")

    # Create solver object
    solver_obj = SolverFactory(solver)
    solver_obj.option = solver_options

    # Perform initial solve of model to confirm feasible initial solution
    results = solver_obj.solve(model, tee=DEBUG)
    
    if not (results.solver.termination_condition == 
            TerminationCondition.optimal and
            results.solver.status == SolverStatus.ok):
        _log.exception("Homotopy Failed - initial solution infeasible.")
    
    print(results)

    # TODO : Check that number of iterations can be accessed from results object

    # Set up homotopy variables
    # Get initial values and deltas for all variables
    v_init = []
    delta = []
    for i in range(len(variables)):
        v_init.append(variables[i].value)
        delta.append(targets[i] - variables[i].value)

    if DEBUG:
        print(v_init, delta)

    n_0 = 0.0  # Homotopy progress variable
    s = step_init  # Set step size to step_init
    iter_count = 0  # Counter for homotopy iterations

    # Save model state to dict
    # TODO : for very large models, it may be necessary to dump this to a file
    current_state = to_json(model, return_dict=True)

    while n_0 < 1.0:
        iter_count += 1  # Increase iter_count regardless of success or failure

        # Calculate next n value given current step size
        if n_0 + s >= 1.0-eps:
            n_1 = 1.0
        else:
            n_1 = n_0 + s

        _log.info("Homotopy Iteration {}. Next Step {} (current step {})"
                  .format(iter_count, n_1, n_0))

        # Update values for all variables using n_1
        for i in range(len(variables)):
            variables[i].fix(v_init[i] + delta[i]*n_1)
            if DEBUG:
                variables[i].display()

        # Solve model at new state
        results = solver_obj.solve(model, tee=DEBUG)

        # Check solver output for convergence
        if (results.solver.termination_condition == 
            TerminationCondition.optimal and
            results.solver.status == SolverStatus.ok):
            # Step succeeded - accept current state
            current_state = to_json(model, return_dict=True)

            # Update n_0 to accept current step
            n_0 = n_1

            # Check solver iterations and calculate next step size
            solver_iterations = 8  # TODO : for now picking a value of iterations
            
            s_proposed = s*(1 + step_accel*(iter_target/solver_iterations-1))
            
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
            if s != min_step:
                # Step size can be cut
                s = max(min_step, s*step_cut)
            else:
                # Step is already at minimum size, terminate homotopy
                _log.exception(
                        "Homotopy failed - could not converge at minimum step "
                        "size. Current progress is {}".format(n_0))
        
        if iter_count >= max_iter:  # Use greater than or equal to to be safe
                _log.exception(
                        "Homotopy failed - maximum homotopy iterations "
                        "exceeded. Current progress is {}".format(n_0))

    _log.info("Homotopy successful - converged at target values in {} "
              "iterations.".format(iter_count))
