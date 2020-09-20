from idaes.apps.caprese.base_class import DynamicBase
from idaes.apps.caprese.util import initialize_by_element_in_range
from idaes.apps.caprese.common.config import (
        VariableCategory,
        ControlPenaltyType,
        ControlInitOption,
        )
from idaes.core.util.model_statistics import degrees_of_freedom
from pyomo.environ import (
        value,
        Objective,
        TerminationCondition,
        Constraint,
        Block,
        Reference,
        TransformationFactory,
        )
from pyomo.kernel import ComponentMap
from pyomo.core.base.range import remainder
from pyomo.dae.set_utils import deactivate_model_at

ATTRIBUTES = (
        'dae_vars',
        'diff_vars',
        'deriv_vars',
        'input_vars',
        'alg_vars',
        'fixed_vars',
        'scalar_vars',
        'ic_vars',
        'variables_categorized',
        'ncp',
        'category_dict',
        'var_locator',
        )

class DynamicModelHelper(object):
    # TODO: Any advantage to inheriting from ConcreteModel or Block?
    #       Would not have to access the "model" attribute
    #       As a Block, could "contain" the newly constructed components?
    # TODO: This class should probably give the option to clone
    # the user's model.

    def __init__(self, model, time, inputs_at_t0):
        self.model = model
        self.time = time

        # TODO: These methods should be moved to this class
        DynamicBase.add_namespace_to(model, time)
        namespace = getattr(model, DynamicBase.namespace_name)
        self.namespace = namespace
        DynamicBase.categorize_variables(model, inputs_at_t0)
        category_dict = {
                VariableCategory.DIFFERENTIAL: namespace.diff_vars,
                VariableCategory.DERIVATIVE: namespace.deriv_vars,
                VariableCategory.ALGEBRAIC: namespace.alg_vars,
                VariableCategory.SCALAR: namespace.scalar_vars,
                VariableCategory.INPUT: namespace.input_vars,
                VariableCategory.FIXED: namespace.fixed_vars,
                }
        namespace.category_dict = category_dict
        DynamicBase.build_variable_locator(
                model, 
                category_dict, 
                ic_vars=namespace.ic_vars)
        for attr in ATTRIBUTES:
            ns_obj = getattr(self.namespace, attr)
            setattr(self, attr, ns_obj)

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
        time = namespace.get_time()
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
        namespace.samples_per_horizon = n_samples

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
        namespace.fe_per_sample = fe_per_sample_dict
        namespace.sample_points = sample_points

    def solve_setpoint(self, solver, require_steady=True):
        controller = self.model
        time = self.time
        t0 = time.first()
        namespace = self.namespace
        category_dict = namespace.category_dict

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
        for var in namespace.ic_vars:
            var[t0].unfix()
        for var in category_dict[VariableCategory.INPUT]:
            was_fixed[var[t0]] = var[t0].fixed
            var[t0].unfix()
        if require_steady == True:
            for var in category_dict[VariableCategory.DERIVATIVE]:
                var[t0].fix(0.0)

        namespace.setpoint_objective.activate()

        # Solve single-time point optimization problem
        dof = degrees_of_freedom(controller)
        if require_steady:
            assert dof == namespace.n_input_vars
        else:
            assert dof == namespace.n_input_vars + namespace.n_diff_vars
        results = solver.solve(controller, tee=True)
        if results.solver.termination_condition == TerminationCondition.optimal:
            pass
        else:
            msg = 'Failed to solve for full state setpoint values'
            raise RuntimeError(msg)

        namespace.setpoint_objective.deactivate()

        # Revert changes. Again, order matters
        if require_steady == True:
            for var in category_dict[VariableCategory.DERIVATIVE]:
                var[t0].unfix()
        for var in namespace.ic_vars:
            var[t0].fix()

        # Reactivate components that were deactivated
        for t, complist in deactivated.items():
            for comp in complist:
                if was_originally_active[comp]:
                    comp.activate()

        # Fix inputs that were originally fixed
        for var in category_dict[VariableCategory.INPUT]:
            if was_fixed[var[t0]]:
                var[t0].fix()

        setpoint_categories = [
                VariableCategory.DIFFERENTIAL,
                VariableCategory.ALGEBRAIC,
                VariableCategory.INPUT,
                VariableCategory.FIXED,
                VariableCategory.DERIVATIVE,
                ]
        for categ in setpoint_categories:
            group = category_dict[categ]
            for i, var in enumerate(group):
                group.set_setpoint(i, var[t0].value)
                # use value attribute instead of function
                # so this doesn't fail if value == None

    def add_setpoint_objective(self, 
            setpoint,
            weights,
            ):
        """
        """
        namespace = self.namespace
        category_dict = namespace.category_dict
        for var, weight in weights:
            locator = namespace.var_locator[var]
            categ = locator.category
            location = locator.location
            category_dict[categ].weights[location] = weight

        weight_vector = []
        for var, sp in setpoint:
            locator = namespace.var_locator[var]
            categ = locator.category
            location = locator.location
            group = category_dict[categ]
            _slice = group[location]
            if group.weights[location] is None:
                print('WARNING: weight not supplied for %s' % var.name)
                group.weights[location] = 1.0
            weight_vector.append(group.weights[location])

        obj_expr = sum(
            weight_vector[i]*(var - sp)**2 for
            i, (var, sp) in enumerate(setpoint))
        namespace.setpoint_objective = Objective(expr=obj_expr)
        self.setpoint_objective = namespace.setpoint_objective

    def add_tracking_objective(self, weights,
            control_penalty_type=ControlPenaltyType.ERROR,
            state_categories=[
                VariableCategory.DIFFERENTIAL,
                ],
            state_weight=1.0, # These are liable to get confused with other weights
            input_weight=1.0,
            objective_weight=1.0,
            ):
        """
        """
        namespace = self.namespace
        sample_points = namespace.sample_points
        # Since t0 ~is~ a sample point, will have to iterate
        # over sample_points[1:]
        n_sample_points = len(sample_points)

        # First set the weight of each variable specified
        for var, weight in weights:
            info = namespace.var_locator[var]
            group = info.group
            loc = info.location
            group.weights[loc] = weight

        if not (control_penalty_type == ControlPenaltyType.ERROR or
                control_penalty_type == ControlPenaltyType.ACTION or
                control_penalty_type == ControlPenaltyType.NONE):
            raise ValueError(
                "control_penalty_type argument must be 'ACTION' or 'ERROR'")

        category_dict = namespace.category_dict
        

        states = [var for c in state_categories 
                for var in category_dict[c].varlist]
        state_weights = [w for c in state_categories 
                for w in category_dict[c].weights]
        state_setpoint = [sp for c in state_categories 
                for sp in category_dict[c].setpoint]
        n_states = len(states)

        inputs = namespace.input_vars.varlist
        input_weights = namespace.input_vars.weights
        input_setpoint = namespace.input_vars.setpoint
        n_inputs = len(inputs)

        state_term = sum(
                state_weights[i]*(states[i][t] - state_setpoint[i])**2
                for i in range(n_states) if (state_weights[i] is not None 
                    and state_setpoint[i] is not None)
                for t in sample_points[1:]
                )

        if control_penalty_type == ControlPenaltyType.ERROR:
            input_term = sum(
                    input_weights[i]*(inputs[i][t] - input_setpoint[i])**2
                    for i in range(n_inputs) if (input_weights[i] is not None
                        and input_setpoint[i] is not None)
                    for t in sample_points[1:]
                    )
            obj_expr = objective_weight*(
                    state_weight*state_term + input_weight*input_term)
        elif control_penalty_type == ControlPenaltyType.ACTION:
            input_term = sum(
                    input_weights[i]*(inputs[i][sample_points[k]] -
                        inputs[i][sample_points[k-1]])**2
                    for i in range(n_inputs) if input_weights[i] is not None
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
        input_indices = list(range(len(namespace.input_vars)))
        def pwc_rule(ns, t, i):
            time = ns.get_time()
            sp_set = set(ns.sample_points)
            if t in sp_set:
                # No need to check for time.first() as it is a sample point
                return Constraint.Skip
            t_next = time.next(t)
            inputs = ns.input_vars
            _slice = inputs[i]
            return _slice[t_next] == _slice[t]

        namespace.pwc_constraint = Constraint(
                time,
                input_indices,
                rule=pwc_rule,
                )
        namespace.pwc_constraint_list = [
                Reference(namespace.pwc_constraint[:, i])
                for i in input_indices
                ]

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
        for categ in categories:
            group = category_dict[categ]
            for _slice, sp in zip(group, group.setpoint):
                for t in time:
                    if t == t0:
                        continue
                    _slice[t].set_value(sp)

    def initialize_from_initial_conditions(self,
            categories=[
                VariableCategory.DERIVATIVE,
                VariableCategory.DIFFERENTIAL,
                VariableCategory.ALGEBRAIC
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
        time = namespace.get_time()
        cat_dict = namespace.category_dict
        for categ in categories:
            varlist = cat_dict[categ].varlist
            for v in varlist:
                v[:].set_value(v[0].value)

    def initialize_by_solving_elements(self, solver):
        namespace = self.namespace
        time = self.time
        model = self.model
        strip_bounds = TransformationFactory('contrib.strip_var_bounds')
        strip_bounds.apply_to(model, reversible=True)

        input_vars = namespace.input_vars
        for var, sp in zip(input_vars, input_vars.setpoint):
            for t in time:
                if t != time.first():
                    var[t].fix(sp)
                else:
                    var[t].fix()
        
        namespace.tracking_objective.deactivate()
        namespace.pwc_constraint.deactivate()

        initialize_by_element_in_range(
                model,
                time,
                time.first(),
                time.last(),
                dae_vars=namespace.dae_vars,
                time_linking_variables=namespace.diff_vars,
                )

        namespace.tracking_objective.activate()
        namespace.pwc_constraint.activate()

        for var in namespace.input_vars:
            for t in time:
                if t == time.first():
                    continue
                var[t].unfix()

        strip_bounds.revert(model)
        # TODO: Can check for violated bounds if I'm really worried about it.
