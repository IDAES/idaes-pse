import time as time_module
from collections import OrderedDict

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
from idaes.apps.caprese.categorize import categorize_dae_variables
from idaes.apps.caprese.nmpc_var import NmpcVar
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
from pyomo.core.base.range import remainder
from pyomo.dae.set_utils import deactivate_model_at
from pyomo.dae.flatten import flatten_dae_components

@declare_custom_block('DynamicBlock')
class _DynamicBlockData(_BlockData):
    # TODO: Any advantage to inheriting from ConcreteModel or Block?
    #       Would not have to access the "model" attribute
    #       As a Block, could "contain" the newly constructed components?
    # TODO: This class should probably give the option to clone
    # the user's model.

    namespace_name = '_CAPRESE_NAMESPACE'

    def __init__(self, model=None, time=None, inputs=None, **kwargs):
        # TODO: Figure out how to properly get the arguments I want into
        # this custom block class.
        super(_DynamicBlockData, self).__init__(**kwargs)
        self.model = model
        self.time = time

        self.add_namespace()
        self.add_time_to_namespace()

        scalar_vars, dae_vars = flatten_dae_components(
                model,
                time,
                Var,
                new_ctype=NmpcVar,
                )
        self.dae_vars = dae_vars
        category_dict = categorize_dae_variables(dae_vars, time, inputs)
        self.category_dict = category_dict
        self.measured_vars = category_dict.pop(VariableCategory.MEASUREMENT)
        # The categories in category_dict now form a partition of the
        # time-indexed variables.

        # Maps each vardata (of a time-indexed var) to the NmpcVar
        # that contains it.
        self.vardata_map = ComponentMap((var[t], var) 
                for varlist in category_dict.values()
                for var in varlist
                for t in time
                )
        # NOTE: looking up var[t] instead of iterating over values() 
        # appears to be ~ 5x faster

        self.add_category_blocks_to_namespace()
        self.add_category_references()

    _var_name = 'var'

    def add_category_blocks_to_namespace(self):
        category_dict = self.category_dict
        namespace = self.namespace
        var_name = self._var_name
        for categ, varlist in category_dict.items():
            categ_name = str(categ).split('.')[1]
            block_name = categ_name + '_BLOCK'
            set_name = categ_name + '_SET'
            set_range = range(len(varlist))

            category_set = Set(initialize=set_range)
            namespace.add_component(set_name, category_set)

            category_block = Block(category_set)
            namespace.add_component(block_name, category_block)
            # Don't want these blocks sent to any solver.
            category_block.deactivate()

            for i, var in enumerate(varlist):
                # Add references to new blocks
                category_block[i].add_component(var_name, var)

    def add_category_references(self):
        namespace = self.namespace
        # TODO: Don't hardcode name of `.var` here
        self.deriv_vars = Reference(namespace.DERIVATIVE_BLOCK[:].var)
        self.diff_vars = Reference(namespace.DIFFERENTIAL_BLOCK[:].var)
        self.alg_vars = Reference(namespace.ALGEBRAIC_BLOCK[:].var)
        self.input_vars = Reference(namespace.INPUT_BLOCK[:].var)
        self.fixed_vars = Reference(namespace.FIXED_BLOCK[:].var)

    def add_namespace(self):
        model = self.model
        namespace_name = self.namespace_name
        namespace_name = unique_component_name(model, namespace_name)
        self.namespace_name = namespace_name

        model.add_component(namespace_name, Block())
        self.namespace = getattr(model, namespace_name)

    def add_time_to_namespace(self):
        # Do this because I can't add a reference to a set
        super(_BlockData, self.namespace).__setattr__('time', self.time)

    def find_components(self, model, comps, time):
        t0_src = time.first()
        t0_tgt = self.time.first()
        tgt_comps = []
        for comp in comps:
            src_vardata = comp[t0_src]
            tgt_vardata = find_comp_in_block_at_time(
                    self.model,
                    model,
                    src_vardata,
                    self.time,
                    t0_tgt,
                    )
            tgt_comps.append(self.vardata_map[tgt_vardata])
        return tgt_comps

    def set_sample_time(self, sample_time, tolerance=1e-8):
        self.validate_sample_time(sample_time, tolerance)
        self.sample_time = sample_time

    def validate_sample_time(self, sample_time, tolerance=1e-8):
        """Makes sure sample points, or integer multiple of sample time-offsets
        from time.first(), lie on finite element boundaries, and that the 
        horizon of each model is an integer multiple of sample time. Assembles 
        a list of sample points and a dictionary mapping sample points to the 
        number of finite elements in the preceding sampling period, and adds 
        them as attributes to _NMPC_NAMESPACE.

        Args:
            sample_time: Sample time to check

        """
        namespace = self.namespace
        time = namespace.time
        horizon_length = time.last() - time.first()

        # TODO: This should probably be a DAE utility
        min_spacing = horizon_length
        for t in time:
            if t == time.first():
                continue
            prev = time.prev(t)
            if t - prev < min_spacing:
                min_spacing = t - prev
        # Sanity check:
        assert min_spacing > 0
        # Required so only one point can satisfy equality to tolerance
        if tolerance >= min_spacing/2:
            raise ValueError(
                'ContinuousSet tolerance is larger than half the minimum '
                'spacing. An element of this set will not necessarily be '
                'unique within this tolerance.')

        off_by = abs(remainder(horizon_length, sample_time))
        if off_by > tolerance:
            raise ValueError(
                'Sampling time must be an integer divider of '
                'horizon length within tolerance %f' % tolerance)
        n_samples = round(horizon_length/sample_time)
        self.samples_per_horizon = n_samples

        finite_elements = time.get_finite_elements()
        sample_points = [time.first()]
        sample_no = 1
        fe_per = 0
        fe_per_sample_dict = {}
        for t in finite_elements:
            if t == time.first():
                continue
            fe_per += 1
            time_since = t - time.first()
            sp = sample_no*sample_time
            diff = abs(sp-time_since)
            if diff < tolerance:
                sample_points.append(t)
                sample_no += 1
                fe_per_sample_dict[sample_no] = fe_per
                fe_per = 0
            if time_since > sp:
                raise ValueError(
                        'Could not find a time point for the %ith '
                        'sample point' % sample_no)
        assert len(sample_points) == n_samples + 1
        self.fe_per_sample = fe_per_sample_dict
        self.sample_points = sample_points

    def initialize(self, option):
        namespace = self.namespace
        time = self.time

        if option == ControlInitOption.FROM_PREVIOUS:
            pass
        elif option == ControlInitOption.BY_TIME_ELEMENT:
            pass
        elif option == ControlInitOption.FROM_INITIAL_CONDITIONS:
            pass
        elif option == ControlInitOption.SETPOINT:
            pass

    def initialize_from_setpoint(self,
            categories=[
                VariableCategory.DIFFERENTIAL,
                VariableCategory.ALGEBRAIC,
                VariableCategory.DERIVATIVE,
                VariableCategory.INPUT,
                ],
            ):
        namespace = self.namespace
        time = self.time
        t0 = time.first()
        category_dict = namespace.category_dict
        for categ, varlist in categories.items():
            # Would like:
            # vector[:,t1:].set_value(vector[:].setpoint)
            for var in varlist:
                for t in time:
                    if t == t0:
                        continue
                    var[t].set_value(var.setpoint)

    def initialize_from_initial_conditions(self,
            categories=[
                VariableCategory.DERIVATIVE,
                VariableCategory.DIFFERENTIAL,
                VariableCategory.ALGEBRAIC,
                ],
            ):
        """ 
        Set values of differential, algebraic, and derivative variables to
        their values at the initial conditions.
        An implicit assumption here is that the initial conditions are
        consistent.

        Args:
            model : Flowsheet model whose variables are initialized
            categories : List of VariableCategory enum items to 
                         initialize. Default contains DERIVATIVE, DIFFERENTIAL,
                         and ALGEBRAIC.

        """
        namespace = self.namespace
        time = self.time
        t0 = time.first()
        cat_dict = self.category_dict
        for categ in categories:
            for v in cat_dict[categ]:
                v[:].set_value(v[t0].value)
        # for var in component_objects(categories):
        #     var[:].set_value(var[t0].value)
        # This should work without a custom var class.

    def initialize_by_solving_elements(self, solver, fix_inputs=False):
        namespace = self.namespace
        time = self.time
        model = self.model
        strip_bounds = TransformationFactory('contrib.strip_var_bounds')
        strip_bounds.apply_to(model, reversible=True)

        if fix_inputs:
            input_vars = self.input_vars
            for var in input_vars[:]:
                sp = var.setpoint
                # Would like:
                # var[t1:].fix(sp)
                for t in time:
                    if t != time.first():
                        var[t].fix(sp)
                    else:
                        var[t].fix()
        
        if hasattr(namespace, 'tracking_objective'):
            namespace.tracking_objective.deactivate()
        if hasattr(namespace, 'pwc_constraint'):
            namespace.pwc_constraint.deactivate()

        initialize_by_element_in_range(
                model,
                time,
                time.first(),
                time.last(),
                dae_vars=self.dae_vars,
                time_linking_vars=list(self.diff_vars[:]),
                outlvl=idaeslog.DEBUG,
                solver=solver,
                )

        if hasattr(namespace, 'tracking_objective'):
            namespace.tracking_objective.activate()
        if hasattr(namespace, 'pwc_constraint'):
            namespace.pwc_constraint.activate()

        if fix_inputs:
            for var in self.input_vars[:]:
                for t in time:
                    if t == time.first():
                        continue
                    var[t].unfix()

        strip_bounds.revert(model)
        # TODO: Can check for violated bounds if I'm really worried about it.
        # Should have a general method to deal with violated bounds that I
        # can use for noise as well.

    def initialize_sample_by_element(self, solver, ts, fix_inputs=False):
        namespace = self.namespace
        time = self.time
        model = self.model
        strip_bounds = TransformationFactory('contrib.strip_var_bounds')
        strip_bounds.apply_to(model, reversible=True)

        t0 = ts - self.sample_time
        idx_0 = time.find_nearest_index(t0)
        t0 = time[idx_0]
        idx_s = time.find_nearest_index(ts)

        if fix_inputs:
            input_vars = self.input_vars
            for var in input_vars[:]:
                sp = var.setpoint
                # Would like:
                # var[t1:].fix(sp)
                for i in range(idx_0+1, idx_s+1):
                    t = time[i]
                    var[t].fix(sp)
        
        if hasattr(namespace, 'tracking_objective'):
            namespace.tracking_objective.deactivate()
        if hasattr(namespace, 'pwc_constraint'):
            namespace.pwc_constraint.deactivate()

        initialize_by_element_in_range(
                model,
                time,
                t0,
                ts,
                dae_vars=self.dae_vars,
                time_linking_vars=list(self.diff_vars[:]),
                outlvl=idaeslog.DEBUG,
                solver=solver,
                )

        if hasattr(namespace, 'tracking_objective'):
            namespace.tracking_objective.activate()
        if hasattr(namespace, 'pwc_constraint'):
            namespace.pwc_constraint.activate()

        if fix_inputs:
            for var in self.input_vars[:]:
                for i in range(idx_0+1, idx_s+1):
                    t = time[i]
                    var[t].unfix()

        strip_bounds.revert(model)

    def generate_inputs_at_time(self, t):
        for val in self.input_vars[:][t].value:
            yield val

    def generate_measurements_at_time(self, t):
        for var in self.measured_vars:
            yield var[t].value

    def inject_inputs(self, inputs):
        # To simulate computational delay, this function would 
        # need an argument for the start time of inputs.
        for var, val in zip(self.input_vars[:], inputs):
            # Would like:
            # self.input_vars[:,:].fix(inputs)
            # This is an example of setting a matrix from a vector.
            # Could even aspire towards:
            # self.input_vars[:,t0:t1].fix(inputs[t1])
            var[:].fix(val)

    def load_measurements(self, measured):
        t0 = self.time.first()
        # Want: self.measured_vars[:,t0].fix(measured)
        for var, val in zip(self.measured_vars, measured):
            var[t0].fix(val)

    def cycle_by(self,
            t_cycle,
            categories=(
                VariableCategory.DIFFERENTIAL,
                VariableCategory.DERIVATIVE,
                VariableCategory.ALGEBRAIC,
                VariableCategory.INPUT,
                VariableCategory.FIXED,
                ),
            tolerance=1e-8,
            ):
        time = self.time
        category_dict = self.category_dict
        for t in time:
            ts = t + t_cycle
            idx = time.find_nearest_index(ts, tolerance)
            if idx is None:
                # t + sample_time is outside the model's "horizon"
                continue
            ts = time[idx]
            for categ in categories:
                # Would like:
                # for var in component_objects(categories):
                for var in category_dict[categ]:
                    var[t].set_value(var[ts].value)

    def cycle(self,
            categories=(
                VariableCategory.DIFFERENTIAL,
                VariableCategory.DERIVATIVE,
                VariableCategory.ALGEBRAIC,
                VariableCategory.INPUT,
                VariableCategory.FIXED,
                ),
            tolerance=1e-8,
            ):
        #horizon = self.time.last() - self.time.last()
        sample_time = self.sample_time
        self.cycle_by(
                sample_time,
                categories=categories,
                tolerance=tolerance,
                )

    def generate_time_in_sample(self,
            ts,
            t0=None,
            tolerance=1e-8,
            ):
        time = self.time
        idx_s = time.find_nearest_index(ts, tolerance=tolerance)
        ts = time[idx_s]
        if t0 is None:
            t0 = ts - self.sample_time
        idx_0 = time.find_nearest_index(t0, tolerance=tolerance)
        for i in range(idx_0+1, idx_s+1):
            # Don't want to include first point in sample
            yield time[i]

    def get_data(self,
            ts,
            variables=(
                VariableCategory.DIFFERENTIAL,
                VariableCategory.INPUT,
                ),
            tolerance=1e-8,
            ):
        time = self.time
        sample_time = self.sample_time
        category_dict = self.category_dict
        vardata_map = self.vardata_map

        data = OrderedDict()
        queue = list(variables)
        for var in queue:
            if var in VariableCategory:
                category = var
                varlist = category_dict[category]
                queue.extend(var[ts] for var in varlist)
                continue
            _slice = vardata_map[var]
            cuid = cuid_from_timeslice(_slice, time)
            data[cuid] = [_slice[t].value for t in
                    self.generate_time_in_sample(ts, tolerance=tolerance)]

        return data


@declare_custom_block('ControllerBlock')
class _ControllerBlockData(_DynamicBlockData):

    def solve_setpoint(self, solver, require_steady=True):
        controller = self.model
        time = self.time
        t0 = time.first()
        namespace = self.namespace
        category_dict = self.category_dict

        was_originally_active = ComponentMap([(comp, comp.active) for comp in 
                controller.component_data_objects((Constraint, Block))])
        non_initial_time = list(time)[1:]
        deactivated = deactivate_model_at(
                controller,
                time,
                non_initial_time,
                allow_skip=True,
                suppress_warnings=True,
                )
        was_fixed = ComponentMap()

        # Fix/unfix variables as appropriate
        # Order matters here. If a derivative is used as an IC, we still want
        # it to be fixed if steady state is required.
        for var in self.measured_vars:
            # Would like:
            # self.measured_vars[:, t0].unfix()
            # Not possible yet because measured_vars are not a reference
            var[t0].unfix()

        input_vars = self.input_vars
        was_fixed = ComponentMap(
                (var[t0], var[t0].fixed) for var in input_vars[:]
                )
        input_vars[:][t0].unfix()
        if require_steady == True:
            self.deriv_vars[:][t0].fix(0.)

        namespace.setpoint_objective.activate()

        # Solve single-time point optimization problem
        dof = degrees_of_freedom(controller)
        if require_steady:
            assert dof == len(namespace.INPUT_SET)
        else:
            assert dof == (len(namespace.INPUT_SET) +
                    len(namespace.DIFFERENTIAL_SET))
        results = solver.solve(controller, tee=True)
        if results.solver.termination_condition == TerminationCondition.optimal:
            pass
        else:
            msg = 'Failed to solve for full state setpoint values'
            raise RuntimeError(msg)

        namespace.setpoint_objective.deactivate()

        # Revert changes. Again, order matters
        if require_steady == True:
            self.deriv_vars[:][t0].unfix()
        for var in self.measured_vars:
            # Again, would like this to be one line
            var[t0].fix()

        # Reactivate components that were deactivated
        for t, complist in deactivated.items():
            for comp in complist:
                if was_originally_active[comp]:
                    comp.activate()

        # Fix inputs that were originally fixed
        for var in input_vars[:]:
        # for var in input_vars[:, t0]:... would be nicer
            if was_fixed[var[t0]]:
                var[t0].fix()

        setpoint_categories = [
                self.diff_vars,
                self.alg_vars,
                self.input_vars,
                self.fixed_vars,
                self.deriv_vars,
                ]
        for vector in setpoint_categories:
            # Would like:
            # vector.set_setpoint(vector[:,0].value)
            # Needs a ctype that's aware of this vectorization
            for var in vector[:]:
                var.setpoint = var[t0].value

        # for vector in component_objects(NmpcVector):
        #     vector.set_setpoint(vector[:,0].value)
        # or
        #     vector.setpoint = vector[:,0].value
        # and this has to actually change the 

    def add_setpoint_objective(self, 
            setpoint,
            weights,
            ):
        """
        """
        namespace = self.namespace
        vardata_map = self.vardata_map
        for vardata, weight in weights:
            nmpc_var = vardata_map[vardata]
            nmpc_var.weight = weight

        weight_vector = []
        for vardata, sp in setpoint:
            nmpc_var = vardata_map[vardata]
            if nmpc_var.weight is None:
                print('WARNING: weight not supplied for %s' % var.name)
                nmpc_var.weight = 1.0
            weight_vector.append(nmpc_var.weight)

        obj_expr = sum(
            weight_vector[i]*(var - sp)**2 for
            i, (var, sp) in enumerate(setpoint))
        namespace.setpoint_objective = Objective(expr=obj_expr)
        self.setpoint_objective = namespace.setpoint_objective

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
