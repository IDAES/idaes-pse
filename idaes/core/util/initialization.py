# -*- coding: UTF-8 -*-
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
This module contains utility functions for initialization of IDAES models.
"""

from pyomo.environ import (
    Block,
    check_optimal_termination,
    Constraint,
    value,
)
from pyomo.network import Arc
from pyomo.dae import ContinuousSet
from pyomo.core.expr.visitor import identify_variables

from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.dyn_utils import (
    get_activity_dict,
    deactivate_model_at,
    deactivate_constraints_unindexed_by,
    fix_vars_unindexed_by,
    get_derivatives_at,
    copy_values_at_time,
    get_implicit_index_of_set,
)
import idaes.logger as idaeslog
from idaes.core.solvers import get_solver

__author__ = "Andrew Lee, John Siirola, Robert Parker"


def fix_state_vars(blk, state_args=None):
    """
    Method for fixing state variables within StateBlocks. Method takes an
    optional argument of values to use when fixing variables.

    Args:
        blk : An IDAES StateBlock object in which to fix the state variables
        state_args : a dict containing values to use when fixing state
                variables. Keys must match with names used in the
                define_state_vars method, and indices of any variables must
                agree.

    Returns:
        A dict keyed by block index, state variable name (as defined by
        define_state_variables) and variable index indicating the fixed status
        of each variable before the fix_state_vars method was applied.
    """
    # For sanity, handle cases where state_args is None
    if state_args is None:
        state_args = {}

    flags = {}
    for k in blk.keys():
        for n, v in blk[k].define_state_vars().items():
            for i in v:
                flags[k, n, i] = v[i].is_fixed()

                # If not fixed, fix at either guess provided or current value
                if not v[i].is_fixed():
                    if n in state_args:
                        # Try to get initial guess from state_args
                        try:
                            if i is None:
                                val = state_args[n]
                            else:
                                val = state_args[n][i]
                        except KeyError:
                            raise ConfigurationError(
                                "Indexes in state_args did not agree with "
                                "those of state variable {}. Please ensure "
                                "that indexes for initial guesses are correct.".format(
                                    n
                                )
                            )
                        v[i].fix(val)
                    else:
                        # No guess, try to use current value
                        if v[i].value is not None:
                            v[i].fix()
                        else:
                            # No initial value - raise Exception before this
                            # gets to a solver.
                            raise ConfigurationError(
                                "State variable {} does not have a value "
                                "assigned. This usually occurs when a Var "
                                "is not assigned an initial value when it is "
                                "created. Please ensure all variables have "
                                "valid values before fixing them.".format(v.name)
                            )

    return flags


def revert_state_vars(blk, flags):
    """
    Method to revert the fixed state of the state variables within an IDAES
    StateBlock based on a set of flags of the previous state.

    Args:
        blk : an IDAES StateBlock
        flags : a dict of bools indicating previous state with keys in the form
                (StateBlock index, state variable name (as defined by
                define_state_vars), var indices).

    Returns:
        None
    """
    for k in blk.keys():
        for n, v in blk[k].define_state_vars().items():
            for i in v:
                try:
                    if not flags[k, n, i]:
                        v[i].unfix()
                except KeyError:
                    raise ConfigurationError(
                        "Indices of flags proved do not match with indices of"
                        "the StateBlock. Please make sure you are using the "
                        "correct StateBlock."
                    )


def propagate_state(
    destination=None, source=None, arc=None, direction="forward", overwrite_fixed=False
):
    """
    This method propagates values between Ports along Arcs. Values can be
    propagated in either direction using the direction argument.

    Args:
        destination (Port): Port to copy values to or None if specifying arc
        source (Port): Port to copy values from or None if specifying arc
        arc (Arc): If arc is provided, use arc to define source and destination
        direction (str): Direction in which to propagate values.
                Default = 'forward' Valid values: 'forward', 'backward'.
        overwrite_fixed (bool): If True overwrite fixed values, otherwise do not
            overwrite fixed values.

    Returns:
        None
    """
    _log = idaeslog.getLogger(__name__)
    # Allow an arc to be passed as a positional arg
    if destination is not None and source is None and destination.ctype is Arc:
        arc = destination
        if source in ("forward", "backward"):
            direction = source
        destination = None
        source = None
    # Check that only arc or source and destination are passed
    if arc is None and (destination is None or source is None):
        raise RuntimeError(
            "In propagate_state(), must provide source and destination or arc"
        )
    if arc is not None:
        if destination is not None or source is not None:
            raise RuntimeError(
                "In propagate_state(), provide only arc or " "source and destination"
            )
        try:
            source = arc.src
            destination = arc.dest
        except AttributeError:
            _log.error(
                "In propagate_state(), unexpected type of arc "
                "argument. Value must be a Pyomo Arc"
            )
            raise

    if direction == "forward":
        pass
    elif direction == "backward":
        source, destination = destination, source
    else:
        raise ValueError(
            "Unexpected value for direction argument: ({}). "
            "Value must be either 'forward' or 'backward'.".format(direction)
        )

    try:
        destination_vars = destination.vars
    except AttributeError:
        _log.error(
            "In propagate_state(), unexpected type of destination "
            "port argument. Value must be a Pyomo Port"
        )
        raise
    try:
        source_vars = source.vars
    except AttributeError:
        _log.error(
            "In propagate_state(), unexpected type of destination "
            "source argument. Value must be a Pyomo Port"
        )
        raise

    # Copy values
    for k, v in destination_vars.items():
        try:
            for i in v:
                if v[i].is_variable_type() and ((not v[i].fixed) or overwrite_fixed):
                    v[i].value = value(source_vars[k][i])
                elif not v[i].is_variable_type():
                    raise TypeError(
                        f"propagate_state() is can only change the value of "
                        f"variables and cannot set a {v[i].ctype}.  This "
                        f"likely indicates either a malformed port or a misuse "
                        f"of propagate_state."
                    )
        except KeyError as e:
            raise KeyError(
                "In propagate_state, variables have incompatible index sets"
            ) from e


# HACK, courtesy of J. Siirola
def solve_indexed_blocks(solver, blocks, **kwds):
    """
    This method allows for solving of Indexed Block components as if they were
    a single Block. A temporary Block object is created which is populated with
    the contents of the objects in the blocks argument and then solved.

    Args:
        solver : a Pyomo solver object to use when solving the Indexed Block
        blocks : an object which inherits from Block, or a list of Blocks
        kwds : a dict of argumnets to be passed to the solver

    Returns:
        A Pyomo solver results object
    """
    # Check blocks argument, and convert to a list of Blocks
    if isinstance(blocks, Block):
        blocks = [blocks]

    try:
        # Create a temporary Block
        tmp = Block(concrete=True)

        nBlocks = len(blocks)

        # Iterate over indexed objects
        for i, b in enumerate(blocks):
            # Check that object is a Block
            if not isinstance(b, Block):
                raise TypeError(
                    "Trying to apply solve_indexed_blocks to "
                    "object containing non-Block objects"
                )
            # Append components of BlockData to temporary Block
            try:
                tmp._decl["block_%s" % i] = i
                tmp._decl_order.append((b, i + 1 if i < nBlocks - 1 else None))
            except Exception:
                raise Exception(
                    "solve_indexed_blocks method failed adding "
                    "components to temporary block."
                )

        # Set ctypes on temporary Block
        tmp._ctypes[Block] = [0, nBlocks - 1, nBlocks]

        # Solve temporary Block
        results = solver.solve(tmp, **kwds)

    finally:
        # Clean up temporary Block contents so they are not removed when Block
        # is garbage collected.
        tmp._decl = {}
        tmp._decl_order = []
        tmp._ctypes = {}

    # Return results
    return results


def initialize_by_time_element(fs, time, **kwargs):
    """
    Function to initialize Flowsheet fs element-by-element along
    ContinuousSet time. Assumes sufficient initialization/correct degrees
    of freedom such that the first finite element can be solved immediately
    and each subsequent finite element can be solved by fixing differential
    and derivative variables at the initial time point of that finite element.

    Args:
        fs : Flowsheet to initialize
        time : Set whose elements will be solved for individually
        solver : Pyomo solver object initialized with user's desired options
        outlvl : IDAES logger outlvl
        ignore_dof : Bool. If True, checks for square problems will be skipped.

    Returns:
        None
    """
    if not fs.is_flowsheet():
        raise TypeError("First arg must be a FlowsheetBlock")
    if not isinstance(time, ContinuousSet):
        raise TypeError("Second arg must be a ContinuousSet")

    if time.get_discretization_info() == {}:
        raise ValueError("ContinuousSet must be discretized")

    scheme = time.get_discretization_info()["scheme"]
    fep_list = time.get_finite_elements()
    nfe = time.get_discretization_info()["nfe"]

    if scheme == "LAGRANGE-RADAU":
        ncp = time.get_discretization_info()["ncp"]
    elif scheme == "LAGRANGE-LEGENDRE":
        msg = "Initialization does not support collocation with Legendre roots"
        raise NotImplementedError(msg)
    elif scheme == "BACKWARD Difference":
        ncp = 1
    elif scheme == "FORWARD Difference":
        ncp = 1
        msg = "Forward initialization (explicit Euler) has not yet been implemented"
        raise NotImplementedError(msg)
    elif scheme == "CENTRAL Difference":
        msg = "Initialization does not support central finite difference"
        raise NotImplementedError(msg)
    else:
        msg = "Unrecognized discretization scheme. "
        "Has the model been discretized along the provided ContinuousSet?"
        raise ValueError(msg)
    # Disallow Central/Legendre discretizations.
    # Neither of these seem to be square by default for multi-finite element
    # initial value problems.

    # Create logger objects
    outlvl = kwargs.pop("outlvl", idaeslog.NOTSET)
    init_log = idaeslog.getInitLogger(__name__, level=outlvl)
    solver_log = idaeslog.getSolveLogger(__name__, level=outlvl)

    ignore_dof = kwargs.pop("ignore_dof", False)
    solver = kwargs.pop("solver", get_solver())
    fix_diff_only = kwargs.pop("fix_diff_only", True)
    # This option makes the assumption that the only variables that
    # link constraints to previous points in time (which must be fixed)
    # are the derivatives and differential variables. Not true if a controller
    # is being present, but should be a good assumption otherwise, and is
    # significantly faster than searching each constraint for time-linking
    # variables.

    if not ignore_dof:
        if degrees_of_freedom(fs) != 0:
            msg = (
                "Original model has nonzero degrees of freedom. This was "
                "unexpected. Use keyword arg igore_dof=True to skip this "
                "check."
            )
            init_log.error(msg)
            raise ValueError("Nonzero degrees of freedom.")

    # Get dict telling which constraints/blocks are already inactive:
    # dict: id(compdata) -> bool (is active?)
    was_originally_active = get_activity_dict(fs)

    # Deactivate flowsheet except at t0, solve to ensure consistency
    # of initial conditions.
    non_initial_time = [t for t in time]
    non_initial_time.remove(time.first())
    deactivated = deactivate_model_at(fs, time, non_initial_time, outlvl=idaeslog.ERROR)

    if not ignore_dof:
        if degrees_of_freedom(fs) != 0:
            msg = (
                "Model has nonzero degrees of freedom at initial conditions."
                " This was unexpected. Use keyword arg igore_dof=True to skip"
                " this check."
            )
            init_log.error(msg)
            raise ValueError("Nonzero degrees of freedom.")

    init_log.info(
        "Model is inactive except at t=0. Solving for consistent initial conditions."
    )
    with idaeslog.solver_log(solver_log, level=idaeslog.DEBUG) as slc:
        results = solver.solve(fs, tee=slc.tee)
    if check_optimal_termination(results):
        init_log.info("Successfully solved for consistent initial conditions")
    else:
        init_log.error("Failed to solve for consistent initial conditions")
        raise ValueError("Solver failed in initialization")

    deactivated[time.first()] = deactivate_model_at(
        fs, time, time.first(), outlvl=idaeslog.ERROR
    )[time.first()]

    # Here, deactivate non-time-indexed components. Do this after solve
    # for initial conditions in case these were used to specify initial
    # conditions
    con_unindexed_by_time = deactivate_constraints_unindexed_by(fs, time)
    var_unindexed_by_time = fix_vars_unindexed_by(fs, time)

    # Now model is completely inactive

    # For each timestep, we need to
    # 1. Activate model at points we're solving for
    # 2. Fix initial conditions (differential variables at previous timestep)
    #    of finite element
    # 3. Solve the (now) square system
    # 4. Revert the model to its prior state

    # This will make use of the following dictionaries mapping
    # time points -> time derivatives and time-differential variables
    derivs_at_time = get_derivatives_at(fs, time, [t for t in time])
    dvars_at_time = {
        t: [d.parent_component().get_state_var()[d.index()] for d in derivs_at_time[t]]
        for t in time
    }

    # Perform a solve for 1 -> nfe; i is the index of the finite element
    init_log.info(
        "Flowsheet has been deactivated. Beginning element-wise initialization"
    )
    for i in range(1, nfe + 1):
        t_prev = time.at((i - 1) * ncp + 1)
        # Non-initial time points in the finite element:
        fe = [time.at(k) for k in range((i - 1) * ncp + 2, i * ncp + 2)]

        init_log.info(f"Entering step {i}/{nfe} of initialization")

        # Activate components of model that were active in the presumably
        # square original system
        for t in fe:
            for comp in deactivated[t]:
                if was_originally_active[id(comp)]:
                    comp.activate()

        # Get lists of derivative and differential variables
        # at initial time point of finite element
        init_deriv_list = derivs_at_time[t_prev]
        init_dvar_list = dvars_at_time[t_prev]

        # Variables that were originally fixed
        fixed_vars = []
        if fix_diff_only:
            for drv in init_deriv_list:
                # Cannot fix variables with value None.
                # Any variable with value None was not solved for
                # (either stale or not included in previous solve)
                # and we don't want to fix it.
                if not drv.fixed:
                    fixed_vars.append(drv)
                if not drv.value is None:
                    drv.fix()
            for dv in init_dvar_list:
                if not dv.fixed:
                    fixed_vars.append(dv)
                if not dv.value is None:
                    dv.fix()
        else:
            for con in fs.component_data_objects(Constraint, active=True):
                for var in identify_variables(con.expr, include_fixed=False):
                    t_idx = get_implicit_index_of_set(var, time)
                    if t_idx is None:
                        continue
                    if t_idx <= t_prev:
                        fixed_vars.append(var)
                        var.fix()

        # Initialize finite element from its initial conditions
        for t in fe:
            copy_values_at_time(
                fs, fs, t, t_prev, copy_fixed=False, outlvl=idaeslog.ERROR
            )

        # Log that we are solving finite element {i}
        init_log.info(f"Solving finite element {i}")

        if not ignore_dof:
            if degrees_of_freedom(fs) != 0:
                msg = (
                    f"Model has nonzero degrees of freedom at finite element"
                    " {i}. This was unexpected. "
                    "Use keyword arg igore_dof=True to skip this check."
                )
                init_log.error(msg)
                raise ValueError("Nonzero degrees of freedom")

        with idaeslog.solver_log(solver_log, level=idaeslog.DEBUG) as slc:
            results = solver.solve(fs, tee=slc.tee)
        if check_optimal_termination(results):
            init_log.info(f"Successfully solved finite element {i}")
        else:
            init_log.error(f"Failed to solve finite element {i}")
            raise ValueError("Failure in initialization solve")

        # Deactivate components that may have been activated
        for t in fe:
            for comp in deactivated[t]:
                comp.deactivate()

        # Unfix variables that have been fixed
        for var in fixed_vars:
            var.unfix()

        # Log that initialization step {i} has been finished
        init_log.info(f"Initialization step {i} complete")

    # Reactivate components of the model that were originally active
    for t in time:
        for comp in deactivated[t]:
            if was_originally_active[id(comp)]:
                comp.activate()

    for con in con_unindexed_by_time:
        con.activate()
    for var in var_unindexed_by_time:
        var.unfix()

    # Logger message that initialization is finished
    init_log.info("Initialization completed. Model has been reactivated")
