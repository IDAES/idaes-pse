# -*- coding: utf-8 -*-
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
A module of helper functions for working with flattened DAE models.
"""

from pyomo.environ import (
    check_optimal_termination,
    Constraint,
    Var,
    SolverFactory,
)
from pyomo.common.collections import ComponentSet, ComponentMap
from pyomo.dae.flatten import flatten_dae_components
from pyomo.dae.set_utils import is_in_block_indexed_by
from pyomo.core.expr.visitor import identify_variables
from pyomo.core.base.constraint import _ConstraintData
from pyomo.core.base.block import _BlockData

from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.dyn_utils import (
    get_activity_dict,
    deactivate_model_at,
    get_implicit_index_of_set,
    get_fixed_dict,
    deactivate_constraints_unindexed_by,
)
from idaes.apps.caprese.common.config import NoiseBoundOption
import idaes.logger as idaeslog

import random

__author__ = "Robert Parker"


class CachedVarsContext(object):
    """This is a simple class to cache the values of variables,
    for instance while a solve is performed, so that they may be loaded
    later.

    For example, to set values of disturbances different than those
    already loaded at t0 when solving for the steady state setpoint:

    >>> ctrl = nmpc.controller
    >>> disturbances = [ctrl.FIXED_BLOCK[0].var, ctrl.FIXED_BLOCK[1].var]
    >>> with CachedVarsContext(disturbances, [t0]):
    >>>     for d in disturbances:
    >>>         d[t0].set_value(1.5)
    >>>     ctrl.solve_setpoint(solver)

    """

    def __init__(self, varlist, tlist):
        if type(tlist) is not list:
            tlist = [tlist]
        self.n_t = len(tlist)
        self.vars = varlist
        self.tlist = tlist
        self.cache = [[None for j in range(self.n_t)] for i in range(len(self.vars))]

    def __enter__(self):
        for i in range(len(self.vars)):
            for j, t in enumerate(self.tlist):
                self.cache[i][j] = self.vars[i][t].value
        return self

    def __exit__(self, a, b, c):
        for i in range(len(self.vars)):
            for j, t in enumerate(self.tlist):
                self.vars[i][t].set_value(self.cache[i][j])


def initialize_by_element_in_range(
    model,
    time,
    t_start,
    t_end,
    time_linking_vars=[],
    dae_vars=[],
    max_linking_range=0,
    **kwargs
):
    """Function for solving a square model, time element-by-time element,
    between specified start and end times.

    Args:
        model : Flowsheet model to solve
        t_start : Beginning of timespan over which to solve
        t_end : End of timespan over which to solve

    Kwargs:
        solver : Solver option used to solve portions of the square model
        outlvl : idaes.logger output level
    """
    solver = kwargs.pop("solver", SolverFactory("ipopt"))
    outlvl = kwargs.pop("outlvl", idaeslog.NOTSET)
    init_log = idaeslog.getInitLogger("nmpc", outlvl)
    solver_log = idaeslog.getSolveLogger("nmpc", outlvl)
    solve_initial_conditions = kwargs.pop("solve_initial_conditions", False)

    # TODO: Move to docstring
    # Variables that will be fixed for time points outside the finite element
    # when constraints for a finite element are activated.
    # For a "normal" process, these should just be differential variables
    # (and maybe derivative variables). For a process with a (PID) controller,
    # these should also include variables used by the controller.
    # If these variables are not specified,

    # Timespan over which these variables will be fixed, counting backwards
    # from the first time point in the finite element (which will always be
    # fixed)
    # Should I specify max_linking_range as an integer number of finite
    # elements, an integer number of time points, or a float in actual time
    # units? Go with latter for now.

    assert t_start in time.get_finite_elements()
    assert t_end in time.get_finite_elements()
    # assert degrees_of_freedom(model) == 0
    # No need to check dof here as we will check right before each solve

    # dae_vars = kwargs.pop('dae_vars', [])
    if not dae_vars:
        scalar_vars, dae_vars = flatten_dae_components(model, time, ctype=Var)
        for var in scalar_vars:
            var.fix()
        deactivate_constraints_unindexed_by(model, time)

    ncp = time.get_discretization_info()["ncp"]

    fe_in_range = [
        i
        for i, fe in enumerate(time.get_finite_elements())
        if fe >= t_start and fe <= t_end
    ]
    t_in_range = [t for t in time if t >= t_start and t <= t_end]

    fe_in_range.pop(0)
    n_fe_in_range = len(fe_in_range)

    was_originally_active = get_activity_dict(model)
    was_originally_fixed = get_fixed_dict(model)

    # Deactivate model
    if not solve_initial_conditions:
        time_list = [t for t in time]
        deactivated = deactivate_model_at(model, time, time_list, outlvl=idaeslog.ERROR)
    else:
        time_list = [t for t in time if t != time.first()]
        deactivated = deactivate_model_at(model, time, time_list, outlvl=idaeslog.ERROR)

        assert degrees_of_freedom(model) == 0
        with idaeslog.solver_log(solver_log, level=idaeslog.DEBUG) as slc:
            results = solver.solve(model, tee=slc.tee)
        if check_optimal_termination(results):
            pass
        else:
            raise ValueError("Failed to solve for consistent initial conditions.")

        deactivated[time.first()] = deactivate_model_at(
            model, time, time.first(), outlvl=idaeslog.ERROR
        )[time.first()]

    # "Integration" loop
    for i in fe_in_range:
        t_prev = time.at((i - 1) * ncp + 1)

        fe = [time.at(k) for k in range((i - 1) * ncp + 2, i * ncp + 2)]

        con_list = []
        for t in fe:
            # These will be fixed vars in constraints at t
            # Probably not necessary to record at what t
            # they occur
            for comp in deactivated[t]:
                if was_originally_active[id(comp)]:
                    comp.activate()
                    if not time_linking_vars:
                        if isinstance(comp, _ConstraintData):
                            con_list.append(comp)
                        elif isinstance(comp, _BlockData):
                            # Active here should be independent of whether block
                            # was active
                            con_list.extend(
                                list(
                                    comp.component_data_objects(Constraint, active=True)
                                )
                            )

        if not time_linking_vars:
            fixed_vars = []
            for con in con_list:
                for var in identify_variables(con.expr, include_fixed=False):
                    # use var_locator/ComponentMap to get index somehow
                    t_idx = get_implicit_index_of_set(var, time)
                    if t_idx is None:
                        assert not is_in_block_indexed_by(var, time)
                        continue
                    if t_idx <= t_prev:
                        fixed_vars.append(var)
                        var.fix()
        else:
            fixed_vars = []
            time_range = [
                t for t in time if t_prev - t <= max_linking_range and t <= t_prev
            ]
            time_range = [t_prev]
            for _slice in time_linking_vars:
                for t in time_range:
                    # if not _slice[t].fixed:
                    _slice[t].fix()
                    fixed_vars.append(_slice[t])

        # Here I assume that the only variables that can appear in
        # constraints at a different (later) time index are derivatives
        # and differential variables (they do so in the discretization
        # equations) and that they only participate at t_prev.
        #
        # This is not the case for, say, PID controllers, in which case
        # I should pass in a list of "complicating variables," then fix
        # them at all time points outside the finite element.
        #
        # Alternative solution is to identify_variables in each constraint
        # that is activated and fix those belonging to a previous finite
        # element. (Should not encounter variables belonging to a future
        # finite element.)
        # ^ This option is easier, less efficient
        #
        # In either case need to record whether variable was previously fixed
        # so I know if I should unfix it or not.

        for t in fe:
            for _slice in dae_vars:
                if not _slice[t].fixed:
                    # Fixed DAE variables are time-dependent disturbances,
                    # whose values should not be altered by this function.
                    _slice[t].set_value(_slice[t_prev].value)

        assert degrees_of_freedom(model) == 0

        with idaeslog.solver_log(solver_log, level=idaeslog.DEBUG) as slc:
            results = solver.solve(model, tee=slc.tee)
        if check_optimal_termination(results):
            pass
        else:
            raise ValueError("Failed to solve for finite element %s" % i)

        for t in fe:
            for comp in deactivated[t]:
                comp.deactivate()

        for var in fixed_vars:
            if not was_originally_fixed[id(var)]:
                var.unfix()

    for t in time:
        for comp in deactivated[t]:
            if was_originally_active[id(comp)]:
                comp.activate()


def get_violated_bounds(val, bounds):
    """This function tests a value against a lower and an upper bound,
    returning which if either is violated, as well as a direction
    that the value needs to move for feasibility.

    Arguments:
        val: Value to be tested
        bounds: Tuple containing the lower, upper bounds

    """
    lower = bounds[0]
    upper = bounds[1]
    if upper is not None:
        if val > upper:
            return (upper, -1)
    if lower is not None:
        if val < lower:
            return (lower, 1)
    return (None, 0)


class MaxDiscardError(Exception):
    pass


def apply_noise(val_list, noise_params, noise_function):
    """
    Applies noise to each value in a list of values and returns the result.
    Noise is generated by a user-provided function that maps a value and
    parameters to a random value.
    """
    result = []
    for val, params in zip(val_list, noise_params):
        if type(params) is not tuple:
            # better be a scalar
            params = (params,)
        result.append(noise_function(val, *params))
    return result


def apply_bounded_noise_discard(
    val, params, noise_function, bounds, max_number_discards
):
    i = 0
    while i <= max_number_discards:
        newval = noise_function(val, *params)

        violated_bound, direction = get_violated_bounds(newval, bounds)
        if violated_bound is None:
            return newval
        else:
            # Discard.
            i += 1

    # This could be caught by the caller to raise a more useful
    # error that includes the variable whose noise violates a
    # bound.
    raise MaxDiscardError("Max number of discards exceeded when applying noise.")


def apply_bounded_noise_push(val, params, noise_function, bounds, bound_push):
    newval = noise_function(val, *params)
    violated_bound, direction = get_violated_bounds(newval, bounds)
    if not violated_bound:
        return newval
    return violated_bound + bound_push * direction


def apply_bounded_noise_fail(val, params, noise_function, bounds):
    newval = noise_function(val, *params)
    violated_bound, direction = get_violated_bounds(newval, bounds)
    if violated_bound:
        raise RuntimeError("Applying noise caused a bound to be violated")
    return newval


def apply_noise_with_bounds(
    val_list,
    noise_params,
    noise_function,
    bound_list,
    bound_option=NoiseBoundOption.DISCARD,
    max_number_discards=5,
    bound_push=0.0,
):
    result = []
    for val, params, bounds in zip(val_list, noise_params, bound_list):
        if type(params) is not tuple:
            # better be a scalar
            # better check: if type(params) not in {sequence_types}...
            params = (params,)

        if bound_option == NoiseBoundOption.DISCARD:
            newval = apply_bounded_noise_discard(
                val, params, noise_function, bounds, max_number_discards
            )
        elif bound_option == NoiseBoundOption.PUSH:
            newval = apply_bounded_noise_push(
                val, params, noise_function, bounds, bound_push
            )
        elif bound_option == NoiseBoundOption.FAIL:
            newval = apply_bounded_noise_fail(val, params, noise_function, bounds)
        else:
            raise RuntimeError("Bound violation option not recognized")

        result.append(newval)
    return result
