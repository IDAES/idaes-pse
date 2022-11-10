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
""" Block-like object meant for controller models.
"""

import idaes.logger as idaeslog
from idaes.apps.caprese.util import initialize_by_element_in_range
from idaes.apps.caprese.common.config import (
    ControlPenaltyType,
)
from idaes.apps.caprese.common.config import VariableCategory as VC
from idaes.apps.caprese.categorize import (
    categorize_dae_variables,
    CATEGORY_TYPE_MAP,
)
from idaes.apps.caprese.nmpc_var import (
    NmpcVar,
    DiffVar,
    AlgVar,
    InputVar,
    DerivVar,
    FixedVar,
    MeasuredVar,
)
from idaes.apps.caprese.dynamic_block import (
    _DynamicBlockData,
    IndexedDynamicBlock,
    DynamicBlock,
)
from idaes.core.util.model_statistics import degrees_of_freedom

from pyomo.environ import (
    Objective,
    check_optimal_termination,
    Constraint,
    Block,
)
from pyomo.core.base.block import _BlockData
from pyomo.common.collections import ComponentMap
from pyomo.core.base.range import remainder
from pyomo.dae.set_utils import deactivate_model_at
from pyomo.dae.flatten import flatten_dae_components
from pyomo.core.base.indexed_component import UnindexedComponent_set


def pwc_rule(ctrl, i, t):
    time = ctrl.time
    sp_set = set(ctrl.sample_points)
    if t in sp_set:
        # No need to check for time.first() as it is a sample point
        return Constraint.Skip
    t_next = time.next(t)
    return ctrl.vectors.input[i, t_next] == ctrl.vectors.input[i, t]


class _ControllerBlockData(_DynamicBlockData):
    """This class adds methods useful for working with dynamic
    models to be used by a controller. These include methods for
    calculating a setpoint and adding objective functions.
    """

    def solve_setpoint(self, solver, **kwargs):
        """This method performs a "real time optimization-type"
        solve on the model, using only the variables and constraints
        at the first point in time. The purpose is to calculate a
        (possibly steady state) setpoint.

        Parameters:
            solver: A Pyomo solver object that will be used to solve
                    the model with only variables and constraints at t0
                    active.
            require_steady: Whether derivatives should be fixed to zero
                            for the solve. Default is `True`.

        """
        require_steady = kwargs.pop("require_steady", True)
        config = self.CONFIG(kwargs)
        model = self.mod
        time = self.time
        t0 = time.first()

        was_originally_active = ComponentMap(
            [
                (comp, comp.active)
                for comp in model.component_data_objects((Constraint, Block))
            ]
        )
        non_initial_time = list(time)[1:]
        deactivated = deactivate_model_at(
            model,
            time,
            non_initial_time,
            allow_skip=True,
            suppress_warnings=True,
        )
        was_fixed = ComponentMap()

        # Cache "important" values to re-load after solve
        init_input = (
            list(self.vectors.input[:, t0].value) if VC.INPUT in self.categories else []
        )
        init_meas = (
            list(self.vectors.measurement[:, t0].value)
            if VC.MEASUREMENT in self.categories
            else []
        )

        # Fix/unfix variables as appropriate
        # Order matters here. If a derivative is used as an IC, we still want
        # it to be fixed if steady state is required.
        if VC.MEASUREMENT in self.categories:
            self.vectors.measurement[:, t0].unfix()

        if VC.INPUT in self.categories:
            input_vars = self.vectors.input
            was_fixed = ComponentMap((var, var.fixed) for var in input_vars[:, t0])
            input_vars[:, t0].unfix()
        if VC.DERIVATIVE in self.categories:
            if require_steady == True:
                self.vectors.derivative[:, t0].fix(0.0)

        self.setpoint_objective.activate()

        # Solve single-time point optimization problem
        #
        # If these sets do not necessarily exist, then I don't think
        # I can make any assertion about the number of degrees of freedom
        # dof = degrees_of_freedom(model)
        # if require_steady:
        #    assert dof == len(self.INPUT_SET)
        # else:
        #    assert dof == (len(self.INPUT_SET) +
        #            len(self.DIFFERENTIAL_SET))
        results = solver.solve(self, tee=config.tee)
        if check_optimal_termination(results):
            pass
        else:
            msg = "Failed to solve for full state setpoint values"
            raise RuntimeError(msg)

        self.setpoint_objective.deactivate()

        # Revert changes. Again, order matters
        if VC.DERIVATIVE in self.categories:
            if require_steady == True:
                self.vectors.derivative[:, t0].unfix()
        if VC.MEASUREMENT in self.categories:
            self.vectors.measurement[:, t0].fix()

        # Reactivate components that were deactivated
        for t, complist in deactivated.items():
            for comp in complist:
                if was_originally_active[comp]:
                    comp.activate()

        # Fix inputs that were originally fixed
        if VC.INPUT in self.categories:
            for var in self.vectors.input[:, t0]:
                if was_fixed[var]:
                    var.fix()

        setpoint_ctype = (DiffVar, AlgVar, InputVar, FixedVar, DerivVar, MeasuredVar)
        for var in self.component_objects(setpoint_ctype):
            var.setpoint = var[t0].value

        # Restore cached values
        if VC.INPUT in self.categories:
            self.vectors.input.values = init_input
        if VC.MEASUREMENT in self.categories:
            self.vectors.measurement.values = init_meas

    def add_setpoint_objective(
        self,
        setpoint,
        weights,
    ):
        """This method uses user-provided setpoints and weights
        to directly construct a least-squares objective for the
        specified variables at the first time point. The purpose
        is to then use the first time point to solve for a "proper
        setpoint" that provides a value for every variable.

        Parameters:
            setpoint: List of vardata, value tuples describing the
                      setpoint values of these specified variables.
            weights: List of vardata, value tuples describing the
                     weight for each variable's term in the objective.

        """
        vardata_map = self.vardata_map
        for vardata, weight in weights:
            nmpc_var = vardata_map[vardata]
            nmpc_var.weight = weight

        weight_vector = []
        for vardata, sp in setpoint:
            nmpc_var = vardata_map[vardata]
            if nmpc_var.weight is None:
                self.logger.warning("Weight not supplied for %s" % var.name)
                nmpc_var.weight = 1.0
            weight_vector.append(nmpc_var.weight)

        obj_expr = sum(
            weight_vector[i] * (var - sp) ** 2 for i, (var, sp) in enumerate(setpoint)
        )
        self.setpoint_objective = Objective(expr=obj_expr)

    def add_tracking_objective(
        self,
        weights,
        control_penalty_type=ControlPenaltyType.ERROR,
        state_ctypes=DiffVar,
        # TODO: Option for user to provide a setpoint here.
        #       (Should ignore setpoint attrs)
        state_weight=1.0,  # These are liable to get confused with other weights
        input_weight=1.0,
        objective_weight=1.0,
    ):
        """This method constructs a tracking objective based on existing
        `setpoint` attributes of this block's `NmpcVar`s and adds it to the
        block. Only specific ctypes listed in `state_ctypes` are directly
        penalized in the objective. These should be the differential variables
        or the measurements, for traditional NMPC.

        Parameters:
            weights: A list of vardata-value tuples describing the weight to
                     be given to the objective term containing each variable.
            control_penalty_type: An entry of the `ControlPenaltyType` enum,
                                  used to determine if control variable error
                                  or action will be penalized.
            state_ctypes: A ctype for the state variables that should be
                          included in the objective function. This should be
                          one of the types from `nmpc_var.py`.
            state_weight: A scalar weight for the state terms in the objective.
            input_weight: A scalar weight for the input terms in the objective.
            objective_weight: A scalar weight for the entire objective.
        """
        samples = self.sample_points
        # Since t0 ~is~ a sample point, will have to iterate
        # over samples[1:]
        n_sample_points = len(samples)

        # First set the weight of each variable specified
        vardata_map = self.vardata_map
        for vardata, weight in weights:
            var = vardata_map[vardata]
            var.weight = weight

        if not (
            control_penalty_type == ControlPenaltyType.ERROR
            or control_penalty_type == ControlPenaltyType.ACTION
            or control_penalty_type == ControlPenaltyType.NONE
        ):
            raise ValueError(
                "control_penalty_type argument must be 'ACTION', 'ERROR', " "or 'NONE'."
            )

        states = list(self.component_objects(state_ctypes))
        inputs = list(self.component_objects(InputVar))

        state_term = sum(
            state.weight * (state[t] - state.setpoint) ** 2
            for state in states
            if state.weight is not None and state.setpoint is not None
            for t in samples[1:]
        )

        if control_penalty_type == ControlPenaltyType.ERROR:
            input_term = sum(
                var.weight * (var[t] - var.setpoint) ** 2
                for var in inputs
                if var.weight is not None and var.setpoint is not None
                for t in samples[1:]
            )
            obj_expr = objective_weight * (
                state_weight * state_term + input_weight * input_term
            )
        elif control_penalty_type == ControlPenaltyType.ACTION:
            input_term = sum(
                var.weight * (var[samples[k]] - var[samples[k - 1]]) ** 2
                for var in inputs
                if var.weight is not None and var.setpoint is not None
                for k in range(1, n_sample_points)
            )
            obj_expr = objective_weight * (
                state_weight * state_term + input_weight * input_term
            )
        elif control_penalty_type == ControlPenaltyType.NONE:
            obj_expr = objective_weight * state_weight * state_term

        self.tracking_objective = Objective(expr=obj_expr)

    def constrain_control_inputs_piecewise_constant(self):
        """Adds an indexed constraint `pwc_constraint` that
        sets each control variable equal to itself at the previous
        time point. This constraint is skipped at sample points.
        """
        time = self.time
        if VC.INPUT in self.categories:
            input_set = self.INPUT_SET
            self.pwc_constraint = Constraint(input_set, time, rule=pwc_rule)
        else:
            raise RuntimeError(
                "Trying to add constraints on inputs but no input "
                "variables have been specified."
            )


class ControllerBlock(DynamicBlock):
    """This is a user-facing class to be instantiated when one
    wants to work with a block capable of the methods of
    `_ControllerBlockData`.
    """

    _ComponentDataClass = _ControllerBlockData

    def __new__(cls, *args, **kwds):
        # Decide what class to allocate
        if cls != ControllerBlock:
            target_cls = cls
        elif not args or (args[0] is UnindexedComponent_set and len(args) == 1):
            target_cls = SimpleControllerBlock
        else:
            target_cls = IndexedControllerBlock
        return super(ControllerBlock, cls).__new__(target_cls)


class SimpleControllerBlock(_ControllerBlockData, ControllerBlock):
    def __init__(self, *args, **kwds):
        _ControllerBlockData.__init__(self, component=self)
        ControllerBlock.__init__(self, *args, **kwds)

    # Pick up the display() from Block and not BlockData
    display = ControllerBlock.display


class IndexedControllerBlock(ControllerBlock):
    def __init__(self, *args, **kwargs):
        ControllerBlock.__init__(self, *args, **kwargs)
