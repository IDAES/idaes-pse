import time as time_module
from collections import OrderedDict

from six import iteritems

import idaes.logger as idaeslog
from idaes.apps.caprese.base_class import DynamicBase
from idaes.apps.caprese.util import (
        initialize_by_element_in_range,
        NMPCVarLocator,
        cuid_from_timeslice,
        )
from idaes.apps.caprese.common.config import (
        VariableCategory,
        ControlPenaltyType,
        ControlInitOption,
        )
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
        )
from idaes.apps.caprese.model import (
        _DynamicBlockData,
        IndexedDynamicBlock,
        DynamicBlock,
        )
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.dyn_utils import (
        find_comp_in_block_at_time,
        )
from pyomo.environ import (
        value,
        Objective,
        TerminationCondition,
        Constraint,
        Block,
        Reference,
        TransformationFactory,
        Var,
        Set,
        )
from pyomo.core.base.block import _BlockData, declare_custom_block
from pyomo.common.collections import ComponentSet, ComponentMap
from pyomo.common.modeling import unique_component_name
from pyomo.common.timing import ConstructionTimer
from pyomo.core.base.range import remainder
from pyomo.dae.set_utils import deactivate_model_at
from pyomo.dae.flatten import flatten_dae_components
from pyomo.core.base.indexed_component import UnindexedComponent_set


class _ControllerBlockData(_DynamicBlockData):

    def solve_setpoint(self, solver, require_steady=True):
        model = self.mod
        time = self.time
        t0 = time.first()

        was_originally_active = ComponentMap([(comp, comp.active) for comp in 
                model.component_data_objects((Constraint, Block))])
        non_initial_time = list(time)[1:]
        deactivated = deactivate_model_at(
                model,
                time,
                non_initial_time,
                allow_skip=True,
                suppress_warnings=True,
                )
        was_fixed = ComponentMap()

        # Fix/unfix variables as appropriate
        # Order matters here. If a derivative is used as an IC, we still want
        # it to be fixed if steady state is required.
        self.vectors.measurement[:,t0].unfix()

        input_vars = self.vectors.input
        was_fixed = ComponentMap(
                (var, var.fixed) for var in input_vars[:,t0]
                )
        input_vars[:,t0].unfix()
        if require_steady == True:
            self.vectors.derivative[:,t0].fix(0.)

        self.setpoint_objective.activate()

        # Solve single-time point optimization problem
        dof = degrees_of_freedom(model)
        if require_steady:
            assert dof == len(self.INPUT_SET)
        else:
            assert dof == (len(self.INPUT_SET) +
                    len(self.DIFFERENTIAL_SET))
        results = solver.solve(self, tee=True)
        if results.solver.termination_condition == TerminationCondition.optimal:
            pass
        else:
            msg = 'Failed to solve for full state setpoint values'
            raise RuntimeError(msg)

        self.setpoint_objective.deactivate()

        # Revert changes. Again, order matters
        if require_steady == True:
            self.vectors.derivative[:,t0].unfix()
        self.vectors.measurement[:,t0].fix()

        # Reactivate components that were deactivated
        for t, complist in deactivated.items():
            for comp in complist:
                if was_originally_active[comp]:
                    comp.activate()

        # Fix inputs that were originally fixed
        for var in self.vectors.input[:,t0]:
            if was_fixed[var]:
                var.fix()

        setpoint_ctype = (DiffVar, AlgVar, InputVar, FixedVar, DerivVar)
        for var in self.component_objects(setpoint_ctype):
            var.setpoint = var[t0].value

    def add_setpoint_objective(self, 
            setpoint,
            weights,
            ):
        """
        """
        vardata_map = self.vardata_map
        for vardata, weight in weights:
            nmpc_var = vardata_map[vardata]
            nmpc_var.weight = weight

        weight_vector = []
        for vardata, sp in setpoint:
            nmpc_var = vardata_map[vardata]
            if nmpc_var.weight is None:
                # TODO: config with outlvl, use logger here.
                print('WARNING: weight not supplied for %s' % var.name)
                nmpc_var.weight = 1.0
            weight_vector.append(nmpc_var.weight)

        obj_expr = sum(
            weight_vector[i]*(var - sp)**2 for
            i, (var, sp) in enumerate(setpoint))
        self.setpoint_objective = Objective(expr=obj_expr)

    def add_tracking_objective(self,
            weights,
            control_penalty_type=ControlPenaltyType.ERROR,
            state_categories=[
                VariableCategory.DIFFERENTIAL,
                ],
            # TODO: Option for user to provide a setpoint here.
            #       (Should ignore setpoint attrs)
            state_weight=1.0, # These are liable to get confused with other weights
            input_weight=1.0,
            objective_weight=1.0,
            ):
        """
        """
        namespace = self.namespace
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
                "or 'None'."
                )

        category_dict = self.category_dict

        states = [var for c in state_categories for var in category_dict[c]]
        inputs = self.input_vars

        state_term = sum(
                state.weight*(state[t] - state.setpoint)**2
                for state in states 
                if state.weight is not None and state.setpoint is not None
                for t in samples[1:]
                )

        if control_penalty_type == ControlPenaltyType.ERROR:
            input_term = sum(
                    var.weight*(var[t] - var.setpoint)**2
                    for var in inputs[:]
                    if var.weight is not None and var.setpoint is not None
                    for t in samples[1:]
                    )
            obj_expr = objective_weight*(
                    state_weight*state_term + input_weight*input_term)
        elif control_penalty_type == ControlPenaltyType.ACTION:
            input_term = sum(
                    var.weight*(var[samples[k]] - var[samples[k-1]])**2
                    for var in inputs[:]
                    if var.weight is not None and var.setpoint is not None
                    for k in range(1, n_sample_points)
                    )
            obj_expr = objective_weight*(
                    state_weight*state_term + input_weight*input_term)
        elif control_penalty_type == ControlPenaltyType.NONE:
            obj_expr = objective_weight*state_weight*state_term

        namespace.tracking_objective = Objective(expr=obj_expr)

    def constrain_control_inputs_piecewise_constant(self):
        namespace = self.namespace
        time = self.time
        input_set = namespace.INPUT_SET

        def pwc_rule(ns, i, t):
            time = self.time
            sp_set = set(self.sample_points)
            if t in sp_set:
                # No need to check for time.first() as it is a sample point
                return Constraint.Skip
            t_next = time.next(t)
            inputs = self.input_vars
            var = inputs[i]
            return var[t_next] == var[t]

        namespace.pwc_constraint = Constraint(input_set, time, rule=pwc_rule)

    def initialize_last_sample(self,
            categories=(
                VariableCategory.DIFFERENTIAL,
                VariableCategory.DERIVATIVE,
                VariableCategory.ALGEBRAIC,
                VariableCategory.INPUT,
                VariableCategory.FIXED,
                ),
            # TODO: include option for how to initialize
            # In particular, how should a known disturbance be initialized?
            tolerance=1e-8,
            ):
        time = self.time
        sample_time = self.sample_time
        t = time.last() - sample_time
        idx = time.find_nearest_index(t, tolerance)
        n = len(time)
        category_dict = self.category_dict
        if type(categories) is VariableCategory:
            categories = (categories,)
        for categ in categories:
            # Would be nice to use component_objects here
            for v in category_dict[categ]:
                for i in range(idx+1, n+1):
                    # Do not initialize the first point in this sample
                    # I do not consider it to be "in the sample"
                    t = time[i]
                    # All variables should have a setpoint attribute,
                    # regardless of category. May not be meaningful
                    # for fixed variables.
                    v[t].set_value(v.setpoint)


class ControllerBlock(DynamicBlock):
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
