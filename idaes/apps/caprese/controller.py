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
from idaes.apps.caprese.dynamic_var import (
        DynamicVar,
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
        TerminationCondition,
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
    """ This class adds methods useful for working with dynamic
    models to be used by a controller. These include methods for
    calculating a setpoint and adding objective functions.
    """

    def _construct(self):
        super(_ControllerBlockData, self)._construct()
        self.has_estimator = False

    def solve_setpoint(self, solver, **kwargs):
        ic_type = kwargs.pop("ic_type", "measurement_var")
        self.solve_single_time_optimization(
            solver, ic_type=ic_type, load_setpoints=True, **kwargs
        )

    def add_tracking_objective(self,
            weights,
            control_penalty_type=ControlPenaltyType.ERROR,
            state_ctypes=DiffVar,
            # TODO: Option for user to provide a setpoint here.
            #       (Should ignore setpoint attrs)
            state_weight=1.0, # These are liable to get confused with other weights
            input_weight=1.0,
            objective_weight=1.0,
            ):
        """ This method constructs a tracking objective based on existing
        `setpoint` attributes of this block's `DynamicVar`s and adds it to the
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

        if not (control_penalty_type == ControlPenaltyType.ERROR or
                control_penalty_type == ControlPenaltyType.ACTION or
                control_penalty_type == ControlPenaltyType.NONE):
            raise ValueError(
                "control_penalty_type argument must be 'ACTION', 'ERROR', "
                "or 'NONE'."
                )

        states = list(self.component_objects(state_ctypes))
        inputs = list(self.component_objects(InputVar))

        state_term = sum(
                state.weight*(state[t] - state.setpoint)**2
                for state in states 
                if state.weight is not None and state.setpoint is not None
                for t in samples[1:]
                )

        if control_penalty_type == ControlPenaltyType.ERROR:
            input_term = sum(
                    var.weight*(var[t] - var.setpoint)**2
                    for var in inputs
                    if var.weight is not None and var.setpoint is not None
                    for t in samples[1:]
                    )
            obj_expr = objective_weight*(
                    state_weight*state_term + input_weight*input_term)
        elif control_penalty_type == ControlPenaltyType.ACTION:
            input_term = sum(
                    var.weight*(var[samples[k]] - var[samples[k-1]])**2
                    for var in inputs
                    if var.weight is not None and var.setpoint is not None
                    for k in range(1, n_sample_points)
                    )
            obj_expr = objective_weight*(
                    state_weight*state_term + input_weight*input_term)
        elif control_penalty_type == ControlPenaltyType.NONE:
            obj_expr = objective_weight*state_weight*state_term

        self.tracking_objective = Objective(expr=obj_expr)

    def constrain_control_inputs_piecewise_constant(self):
        """ Adds an indexed constraint `pwc_constraint` that
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

    def load_estimates(self, estimates):
        tp = self.time.first()
        for var, val in zip(self.DIFFERENTIAL_BLOCK[:].var, estimates):
            var[tp].fix(val)

    def load_initial_conditions(self, ics, type_of_ics = None):
        '''
        Load initial conditions for the NMPC controller. 
        If MHE exists, initial conditions should be the estimates from MHE.
        If MHE doesn't exist, initial conditions are the measurements directly
        from the plant.

        Parameters
        ----------
        ics : list of initial conditions.
        '''

        if type_of_ics is None:
            print("Loading type of initial conditions is not declared. \n"
                  "Determine the type by checking whether MHE exists or not.")
            if self.has_estimator:
                type_of_ics = "estimate"
            else:
                type_of_ics = "measurement"

        if type_of_ics == "estimate":
            self.load_estimates(ics)
        elif type_of_ics == "measurement":
            self.load_measurements(
                ics, timepoint = self.time.first()
            )
        else:
            raise RuntimeError("Declared type of initial conditions does not support."
                               "Please use either 'measurement' or 'estimate'.")

    def load_measurements(self, measured, timepoint=None):
        super(_ControllerBlockData, self).load_measurements(
            measured, "measurement", timepoint
        )


class ControllerBlock(DynamicBlock):
    """ This is a user-facing class to be instantiated when one
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
