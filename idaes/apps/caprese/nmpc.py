# -*- coding: utf-8 -*-
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
Class for performing NMPC simulations of IDAES flowsheets
"""

from pyomo.environ import (
        Block, 
        Constraint, 
        Var, 
        TerminationCondition,
        SolverFactory, 
        Objective, 
        NonNegativeReals, 
        Reals, 
        TransformationFactory, 
        Reference, 
        value,
        )
from pyomo.core.base.range import remainder
from pyomo.common.collections import ComponentMap
from pyomo.dae.initialization import (
        solve_consistent_initial_conditions,
        get_inconsistent_initial_conditions,
        )
from pyutilib.misc.config import ConfigDict, ConfigValue

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import (
        degrees_of_freedom, 
        activated_equalities_generator,
        )
from idaes.core.util.dyn_utils import (
        deactivate_model_at,
        path_from_block, 
        find_comp_in_block, 
        find_comp_in_block_at_time,
        )
from idaes.apps.caprese.common.config import (
        ControlInitOption,
        ElementInitializationInputOption,
        TimeResolutionOption,
        ControlPenaltyType,
        VariableCategory)
from idaes.apps.caprese.util import (
        initialize_by_element_in_range,
        find_slices_in_model, 
        NMPCVarGroup, 
        NMPCVarLocator, 
        copy_values_at_time, 
        add_noise_at_time,
        validate_list_of_vardata, 
        validate_list_of_vardata_value_tuples, 
        validate_solver,
        find_point_in_continuousset,
        get_violated_bounds_at_time)
from idaes.apps.caprese.base_class import DynamicBase
import idaes.logger as idaeslog

__author__ = "Robert Parker and David Thierry"


class NMPCSim(DynamicBase):
    """
    Main class for NMPC simulations of Pyomo models.
    """
    # pyomo.common.config.add_docstring_list
    CONFIG = DynamicBase.CONFIG

    # TODO: How to document config values?
    CONFIG.declare(
            'control_init_option',
            ConfigValue(
                default=ControlInitOption.FROM_INITIAL_CONDITIONS,
                domain=ControlInitOption.from_enum_or_string,
                doc='Option for how to initialize the controller model'
                )
            )
    CONFIG.declare(
            'element_initialization_input_option',
            ConfigValue(
                default=ElementInitializationInputOption.SET_POINT,
                domain=ElementInitializationInputOption.from_enum_or_string,
                doc=('Option for how to fix inputs when initializing '
                    'by time element')
                )
            )
    CONFIG.declare(
            'time_resolution_option',
            ConfigValue(
                default=TimeResolutionOption.SAMPLE_POINTS,
                domain=TimeResolutionOption.from_enum_or_string,
                doc=('Option for specifying a time resolution in the '
                    'objective function')
                )
            )
    CONFIG.declare(
            'calculate_error',
            ConfigValue(
                default=True,
                domain=bool,
                doc=('Flag for whether or not to calculate set-point-error '
                    'when simulating plant')
                )
            )
    CONFIG.declare(
            'state_objective_weight_matrix_diagonal',
            ConfigValue(
                default=True,
                domain=bool,
                doc='Flag for whether state objective weights are diagonal'
                )
            )
    CONFIG.declare(
            'control_objective_weight_matrix_diagonal',
            ConfigValue(
                default=True,
                domain=bool,
                doc='Flag for whether control objective weights are diagonal'
                )
            )
    CONFIG.declare(
            'control_penalty_type',
            ConfigValue(
                default=ControlPenaltyType.ERROR,
                domain=ControlPenaltyType.from_enum_or_string,
                doc=('Type of control penalty that will be normed in '
                    'objective functions')
                )
            )
    # TODO: Should I combine these into one config argument, then just override
    # for each's function if they need to change?
    CONFIG.declare(
            'add_plant_noise',
            ConfigValue(
                default=True,
                domain=bool,
                doc='Flag for whether to add noise to state loaded from plant'
                )
            )
    CONFIG.declare(
            'add_input_noise',
            ConfigValue(
                default=True,
                domain=bool,
                doc=('Flag for whether to add noise to inputs injected '
                    'into plant')
                )
            )
    CONFIG.declare(
            'noise_weights',
            ConfigValue(
                default=[],
                domain=list,
                doc=('List of weights to override weights for variance '
                    'in noise function')
                )
            )
    # ^ TODO: Really this should be a list of vardata, value tuples
    CONFIG.declare(
            'max_noise_weight',
            ConfigValue(
                default=1e6,
                domain=float,
                doc='Maximum value by which noise variance can be weighted'
                )
            )
    CONFIG.declare(
            'noise_arguments',
            ConfigValue(
                default={},
                domain=dict,
                doc='Extra arguments for noise function')
            )
    CONFIG.declare(
            'noise_sigma_0',
            ConfigValue(
                default=0.05,
                domain=float,
                doc=('Nominal value of variance that will be scaled by weights '
                    'for each state')
                )
            )
    CONFIG.declare(
            'setpoint',
            ConfigValue(
                default=[],
                domain=validate_list_of_vardata_value_tuples,
                doc=('User-specified list of VarDatas and their corresponding '
                    'setpoints')
                )
            )
    CONFIG.declare('objective_weight_tolerance',
            ConfigValue(
                default=1e-6,
                domain=float,
                doc=('Minimum delta between nominal and set-point that will '
                    'be used to calculate objective function weights')
                )
            )
    CONFIG.declare('objective_weight_override',
            ConfigValue(
                default=[],
                domain=validate_list_of_vardata_value_tuples,
                doc=('User-specified objective weight values for given '
                    'variables that take precedence over calculated values')
                )
            )
    CONFIG.declare('objective_state_categories',
            ConfigValue(
                default=[VariableCategory.DIFFERENTIAL],
                domain=list,
                doc=('Variable categories that will be penalized as '
                    'states in the objective function'),
                )
            )
    CONFIG.declare('sample_time',
            ConfigValue(
                default=1,
                domain=float,
                doc='Time period over which inputs will be held'
                )
            )
    CONFIG.declare('inputs_at_t0',
            ConfigValue(
                default=[],
                domain=validate_list_of_vardata,
                doc=('List of VarData objects corresponding to the inputs '
                    'at time.first() in the plant model')
                )
            )
    CONFIG.declare('user_objective_name',
            ConfigValue(
                default='user_objective',
                domain=str,
                doc=('Name for the objective function created from the '
                    'set-point provided by the user')
                )
            )
    CONFIG.declare('full_state_objective_name',
            ConfigValue(
                default='tracking_objective',
                domain=str,
                doc=('Name for full-state objective function calculated '
                    'from that provided by the user')
                )
            )


    namespace_name = '_NMPC_NAMESPACE'
    
    @classmethod
    def get_namespace_name(cls):
        return cls.namespace_name
    

    def __init__(self, plant_model=None, plant_time_set=None, 
        controller_model=None, controller_time_set=None, inputs_at_t0=None,
        sample_time=None, **kwargs):
        """Constructor method. Accepts plant and controller models needed for 
        NMPC simulation, as well as time sets (Pyomo Sets) in each model
        Inputs at the first time point in the plant model are also required.
        Models provided are added to the NMPCSim instance as attributes.
        This constructor solves for consistent initial conditions 
        in the plant and controller and performs categorization into lists of
        differential, derivative, algebraic, input, fixed, and scalar variables,
        which are added as attributes to a _NMPC_NAMESPACE Block on each model.

        Args:
            plant_model : Plant Pyomo model, NMPC of which will be 
                          simulated. Currently this must contain the entire 
                          timespan it is desired to simulate.
            plant_time_set : Set to treat as time in the plant model
            controller_model : Model to be used to calculate control inputs
                               for the plant. Control inputs in controller
                               must exist in the plant, and initial condition
                               variables in the plant must exist in the 
                               controller.
            controller_time_set : Set to treat as time in the controller model
            inputs_at_t0 : List of VarData objects containing the variables
                             to be treated as control inputs, at time.first().
            solver : Solver to be used for verification of consistent initial 
                     conditions, will also be used as the default solver if
                     another is not provided for initializing or solving the 
                     optimal control problem.
            outlvl : IDAES logger output level. Default is idaes.logger.INFO.
                     To see solver output, use idaes.logger.DEBUG.
            sample_time : Length of time each control input will be held for.
                          This must be an integer multiple of the (finite
                          element) discretization spacing in both the plant
                          and controller models. Default is to use the 
                          controller model's discretization spacing.

        """
        self.config = self.CONFIG(kwargs)
        super(NMPCSim, self).__init__(plant_model, plant_time_set,
                controller_model, controller_time_set, inputs_at_t0,
                **kwargs)

        # Should I provide solver and outlvl as explicit args here?
        self.config.sample_time = sample_time
        self.config.inputs_at_t0 = inputs_at_t0
        # Maybe include a kwarg for require_steady - if False, set-point is not
        # forced to be a steady state

        # TODO: validate_time_set function

        init_log = idaeslog.getInitLogger('nmpc', level=self.config.outlvl)

        # Only need to manipulate bounds of controller model. Assume the 
        # bounds in the plant model should remain in place for simulation.
        # (Should probably raise a warning if bounds are present...)
        for categ, vargroup in self.controller._NMPC_NAMESPACE.category_dict.items():
            self.set_bounds_from_initial(vargroup)
        # ^ This may be removed in favor of strip_bounds transformation

        # Validate inputs in the plant model and initial conditions
        # in the control model.
        # TODO: allow user to specify this if names don't match
        self.plant._NMPC_NAMESPACE.controller_ic_vars = find_slices_in_model(
                self.plant, self.plant_time,
                self.controller, self.controller_time,
                self.plant._NMPC_NAMESPACE.var_locator,
                self.controller._NMPC_NAMESPACE.ic_vars)
        self.controller._NMPC_NAMESPACE.plant_input_vars = find_slices_in_model(
                self.controller, self.controller_time,
                self.plant, self.plant_time,
                self.controller._NMPC_NAMESPACE.var_locator,
                self.plant._NMPC_NAMESPACE.input_vars.varlist)

        self.validate_fixedness(self.plant, self.controller)

        self.sample_time = self.config.sample_time
        self.validate_sample_time(self.sample_time, 
                self.controller, self.plant)

        # Flag for whether controller has been initialized
        # by a previous solve
        self.controller_solved = False

        # Maps sample times in plant model to the normalized state error
        # This error will be defined by:
        # <(x_pred-x_meas), Q(x_pred-x_meas)>
        # where Q is the positive semi-definite matrix defining the norm
        # used in the objective function.
        #
        # Currently only diagonal matrices Q are supported, and values of None
        # are interpreted as zeros
        self.state_error = {}
        # Should I set state_error[0] = 0? Probably not, in case there is for
        # instance some measurement noise. 
        # Remember: Need to calculate weight matrices before populating this. 

        self.current_plant_time = 0

    @classmethod
    def add_namespace_to(cls, model, time):
        """Adds the _NMPC_NAMESPACE block a model with a given time set.
        All necessary model-specific attributes, including constraints
        and objectives, will be added to this block.

        Args:
            model : Model to which to add the namespace
            time : Set to treat as time in the given model

        """
        name = '_NMPC_NAMESPACE'
        # Not _CAPRESE_NAMESPACE as I might want to add a similar 
        # namespace for MHE
        if hasattr(model, name):
            raise ValueError('%s already exists on model. Please fix this.'
                             % name)
        model.add_component(name, Block())
        super(NMPCSim, cls).add_namespace_to(model, time)
        # TODO: Should this method call _populate_namespace?
        # Otherwise this namespace doesn't have access to the get_time function.
        # Don't want to require that a model is categorized just to get_time,
        # which is what I'm doing right now, unless I call get_time on the base
        # namespace.

    def validate_sample_time(self, sample_time, *models, **kwargs):
        """Makes sure sample points, or integer multiple of sample time-offsets
        from time.first(), lie on finite element boundaries, and that the 
        horizon of each model is an integer multiple of sample time. Assembles 
        a list of sample points and a dictionary mapping sample points to the 
        number of finite elements in the preceding sampling period, and adds 
        them as attributes to _NMPC_NAMESPACE.

        Args:
            sample_time: Sample time to check
            models: List of flowsheet models to check

        """
        config = self.config(kwargs)
        tolerance = config.continuous_set_tolerance
        for model in models:
            time = model._NMPC_NAMESPACE.get_time()
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
            model._NMPC_NAMESPACE.samples_per_horizon = n_samples

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
            model._NMPC_NAMESPACE.fe_per_sample = fe_per_sample_dict
            model._NMPC_NAMESPACE.sample_points = sample_points


    def validate_slices(self, tgt_model, src_model, src_time, src_slices):
        """
        Given list of time-only slices in a source model, attempts to find
        each of them in the target model and returns a list of the found 
        slices in the same order.
        Expects to find a var_locator ComponentMap attribute in the 
        _NMPC_NAMESPACE of the target model.

        Args:
            tgt_model : Model to search for time-slices
            src_model : Model containing the slices to search for
            src_slices : List of time-only slices of variables in the source
                         model

        Returns:
            List of time-only slices to same-named variables in the target 
            model
        """
        t0 = src_time.first()
        tgt_slices = []
        namespace = getattr(tgt_model, self.namespace_name)
        locator = namespace.var_locator
        for _slice in src_slices:
            init_vardata = _slice[t0]
            # FIXME
            # This assumes that t0 is a valid time point for the target
            # model, even it is taken from the source.
            # Should use find_comp_in_block_at_time, which essentially
            # does the work of constructing a CUID with wildcard, but
            # here is tied to the target model's time set.
            # A better validation method might be:
            #     src_comp -> cuid w/ wildcard -> target_comp...
            # this would still be the same amount of work/code in this
            # method. Would be nice to go straight from the source
            # Reference to the CUID, and from the CUID to the/a ref-
            # to-slice. cuid.to_reference() would be nice.
            tgt_vardata = find_comp_in_block(tgt_model, 
                                             src_model, 
                                             init_vardata)
            tgt_container = locator[tgt_vardata].group.varlist
            location = locator[tgt_vardata].location
            tgt_slices.append(tgt_container[location])
        return tgt_slices


    def validate_fixedness(self, *models):
        """
        Makes sure that assumptions regarding fixedness for different points
        in time are valid. Differential, algebraic, and derivative variables
        may be fixed only at t0, only if they are initial conditions.
        Fixed variables must be fixed at all points in time, except possibly
        initial conditions. 

        Expects to find "alg," "diff," "deriv," and "fixed" vars on each
        model's _NMPC_NAMESPACE, as well as a var_locator ComponentMap.

        Args:
            models: Models for which to validate fixedness

        """
        for model in models:
            time = model._NMPC_NAMESPACE.get_time()
            t0 = time.first()
            locator = model._NMPC_NAMESPACE.var_locator

            # Appropriate for this function to have categories specified
            for _slice in (model._NMPC_NAMESPACE.alg_vars.varlist + 
                           model._NMPC_NAMESPACE.diff_vars.varlist + 
                           model._NMPC_NAMESPACE.deriv_vars.varlist):
                var0 = _slice[t0]
                if locator[var0].is_ic:
                    assert var0.fixed
                    for t in time:
                        if t == t0:
                            continue
                        assert not _slice[t].fixed
                else:
                    for t in time:
                        assert not _slice[t].fixed

            for var in model._NMPC_NAMESPACE.fixed_vars.varlist:
                for t in time:
                    # Fixed vars, e.g. those used in boundary conditions,
                    # may "overlap" with initial conditions. It is up to the user
                    # to make sure model has appropriate number of degrees of
                    # freedom
                    if t == t0:
                        continue
                    assert var[t].fixed
                    

    def transfer_current_plant_state_to_controller(self, t_plant, **kwargs):
        """Transfers values of the initial condition variables at a specified
        time in the plant model to the initial time point of the controller
        model, adding noise if desired.

        Args:
            t_plant: Time point in plant model whose values will be transferred

        """
        # Would like to pass "noise_args" in as a bundle here. This can
        # probably be done with config blocks somehow.
        # TODO: allow specification of noise args
        config = self.config(kwargs)

        time = self.controller_time
        t0 = time.first()

        copy_values_at_time(self.controller._NMPC_NAMESPACE.ic_vars,
                            self.plant._NMPC_NAMESPACE.controller_ic_vars,
                            t0,
                            t_plant)

        # Apply noise to new initial conditions
        add_noise = config.add_plant_noise

        noise_weights = config.noise_weights
        noise_sig_0 = config.noise_sigma_0
        noise_args = config.noise_arguments
        max_noise_weight = config.max_noise_weight

        locator = self.controller._NMPC_NAMESPACE.var_locator
        if add_noise:
            if not noise_weights:
                noise_weights = []
                for var in self.controller._NMPC_NAMESPACE.ic_vars:
                    info = locator[var[t0]]
                    loc = info.location
                    obj_weight = info.group.weights[loc]

                    if obj_weight is not None and obj_weight != 0:
                        noise_weights.append(min(1/obj_weight, 
                                                 max_noise_weight))
                    else:
                        noise_weights.append(None)

            add_noise_at_time(self.controller._NMPC_NAMESPACE.ic_vars,
                              t0,
                              weights=noise_weights,
                              sigma_0=noise_sig_0,
                              **noise_args)


    def inject_control_inputs_into_plant(self, t_plant, **kwargs):
        """Injects input variables from the first sampling time in the 
        controller model to the sampling period in the plant model that
        starts at the specified time, adding noise if desired.

        Args:
            t_plant : First time point in plant model where inputs will be
                      applied.
            
        """
        # config args for control_input_noise
        config = self.config(kwargs)
        tolerance = config.continuous_set_tolerance
        sample_time = self.config.sample_time

        # Send inputs to plant that were calculated for the end
        # of the first sample
        t_controller = find_point_in_continuousset(
                self.controller_time.first() + sample_time, 
                self.controller_time, tolerance)
        assert t_controller in self.controller_time

        time = self.plant_time
        plant_sample_end = find_point_in_continuousset(
                t_plant + sample_time, 
                time, tolerance)
        assert plant_sample_end in time
        plant_sample = [t for t in time if t > t_plant and t<= plant_sample_end]
        assert plant_sample_end in plant_sample
        # len(plant_sample) should be ncp*nfe_per_sample, assuming the expected
        # sample_time is passed in

        add_noise = config.add_input_noise
        noise_weights = config.noise_weights
        noise_sig_0 = config.noise_sigma_0
        noise_args = config.noise_arguments
        max_noise_weight = config.max_noise_weight

        # Need to get proper weights for plant's input vars
        locator = self.controller._NMPC_NAMESPACE.var_locator
        if add_noise:
            if not noise_weights:
                noise_weights = []
                for var in self.controller._NMPC_NAMESPACE.plant_input_vars:
                    info = locator[var[t_controller]]
                    loc = info.location
                    obj_weight = info.group.weights[loc]
                    if obj_weight is not None and obj_weight != 0:
                        noise_weights.append(min(1/obj_weight, max_noise_weight))
                    else:
                        # By default, if state is not penalized in objective,
                        # noise will not be applied to it here.
                        # This may be incorrect, but user will have to override,
                        # by providing their own weights, as I don't see a good
                        # way of calculating a weight
                        noise_weights.append(None)

            add_noise_at_time(self.controller._NMPC_NAMESPACE.plant_input_vars,
                              t_controller,
                              weights=noise_weights,
                              sigma_0=noise_sig_0,
                              **noise_args)
            #add_noise_at_time(self.plant.input_vars,
            #                  t_plant+sample_time,
            #                  weights=noise_weights,
            #                  sigma_0=noise_sig_0,
            #                  **noise_args)
            # Slight bug in logic here: noise is applied to plant variables,
            # but only controller variables have bounds.
            # Alternatives: add bounds to plant variables (undesirable)  
            #               apply noise to controller variables (maybe okay...)
            #                ^ can always record nominal values, then revert
            #                  noise after it's copied into plant...
            # Right now I apply noise to controller model, and don't revert

        copy_values_at_time(self.plant._NMPC_NAMESPACE.input_vars.varlist,
                            self.controller._NMPC_NAMESPACE.plant_input_vars,
                            plant_sample,
                            t_controller)

    def has_consistent_initial_conditions(self, model, **kwargs):
        """
        Finds constraints at time.first() that are violated by more than
        tolerance. Returns True if any are found.
        """
        # This will raise an error if any constraints at t0 cannot be
        # evaluated, i.e. contain a variable of value None.
        namespace = getattr(model, self.get_namespace_name())
        time = namespace.get_time()
        config = self.config(kwargs)
        tolerance = config.tolerance
        inconsistent = get_inconsistent_initial_conditions(
                model, 
                time, 
                tol=tolerance,
                suppress_warnings=True)
        return not inconsistent

    def solve_consistent_initial_conditions(self, model, **kwargs):
        """
        Uses pyomo.dae.initialization solve_consistent_initial_conditions
        function to solve for consistent initial conditions. Inputs are
        fixed at time.first() in attempt to eliminate degrees of freedom.
        """
        namespace = getattr(model, self.get_namespace_name())
        time = namespace.get_time()
        strip_bounds = kwargs.pop('strip_bounds', True)
        config = self.config(kwargs)
        outlvl = config.outlvl
        solver = config.solver
        solver_log = idaeslog.getSolveLogger('nmpc', level=outlvl)
        t0 = time.first()

        previously_fixed = ComponentMap()
        for var in namespace.input_vars:
            var0 = var[t0]
            previously_fixed[var0] = var0.fixed
            var0.fix()

        if strip_bounds:
            strip_var_bounds = TransformationFactory(
                                           'contrib.strip_var_bounds')
            strip_var_bounds.apply_to(model, reversible=True)

        with idaeslog.solver_log(solver_log, level=idaeslog.DEBUG) as slc:
            result = solve_consistent_initial_conditions(model, time, solver,
                    tee=slc.tee)

        if strip_bounds:
            strip_var_bounds.revert(model)

        for var, was_fixed in previously_fixed.items():
            if not was_fixed:
                var.unfix()

        return result


    def calculate_full_state_setpoint(self, 
            setpoint, 
            require_steady=True, 
            allow_inconsistent=True,
            **kwargs):
        """Given a user-defined setpoint, i.e. a list of VarData, value tuples,
        calculates a full-state setpoint to be used in the objective function
        of the dynamic optimization problem. This is done by solving a single-
        time point optimization problem with the user's setpoint in the 
        objective function.

        The solve is performed in the first time point blocks/constraints of the
        controller model. The procedure is:

            i. Check for inconsistent initial conditions. Warn user if found.
            ii. Populate controller setpoint attributes with user-defined 
                values.
            iii. Populate reference attributes with (now consistent) initial
                 conditions.
            iv. Calculate weights for variables specified.
            v. Add objective function based on these weights and setpoint
                 values.
            vi. Solve for setpoint.
            vii. Deactivate just-added objective function.

        Args:
            setpoint : List of VarData, value tuples to be used in the objective
                       function of the single-time point optimization problem
            require_steady : Bool telling whether or not to fix derivatives to
                             zero when performing optimization

        """
        config = self.config(kwargs)
        solver = config.solver
        outlvl = config.outlvl
        tolerance = config.tolerance
        init_log = idaeslog.getInitLogger('nmpc', outlvl)
        user_objective_name = config.user_objective_name

        controller = self.controller
        time = self.controller_time
        t0 = time.first()
        category_dict = controller._NMPC_NAMESPACE.category_dict
        locator = controller._NMPC_NAMESPACE.var_locator

        # User should have already solved for consistent initial conditions if
        # they want them.
        inconsistent = get_inconsistent_initial_conditions(
                controller, 
                time,
                tol=tolerance,
                suppress_warnings=True)
        if inconsistent:
            msg = ('Initial conditions are inconistent. Weights in the '
            'setpoint optimization problem may not be reasonable. Use the '
            'solve_consistent_initial_conditions before calling '
            'calculate_full_state_setpoint to remedy')
            if allow_inconsistent:
                init_log.warning(msg)
            else:
                raise RuntimeError(msg)

        # Categories of variables whose set point values will be added to controller
        # TODO: maybe this should be an argument
        categories = [VariableCategory.DIFFERENTIAL,
                      VariableCategory.ALGEBRAIC,
                      VariableCategory.DERIVATIVE,
                      VariableCategory.INPUT,
                      VariableCategory.SCALAR]
    
        # Clear any previous existing setpoint values in these variables
        for categ in categories:
            group = category_dict[categ]
            for i in range(group.n_vars):
                group.set_setpoint(i, None)

        # Populate appropriate setpoint values from argument
        for vardata, val in setpoint:
            info = locator[vardata]
            categ = info.category
            loc = info.location
            group = category_dict[categ]
            group.set_setpoint(loc, val)

        # Calculate objective weights for all variables.
        for categ, vargroup in category_dict.items():
            if categ == VariableCategory.SCALAR:
                for i, var in enumerate(vargroup):
                    vargroup.set_reference(i, var.value)
            else:
                for i, var in enumerate(vargroup):
                    vargroup.set_reference(i, var[t0].value)

        override = config.objective_weight_override
        tolerance = config.objective_weight_tolerance

        self.construct_objective_weights(
                controller,
                objective_weight_override=override,
                objective_weight_tolerance=tolerance,
                categories=[VariableCategory.DIFFERENTIAL,
                            VariableCategory.ALGEBRAIC,
                            VariableCategory.DERIVATIVE,
                            VariableCategory.INPUT])
                
        # Save user setpoint and weights as attributes of namespace
        # in case they are required later
        user_setpoint = []
        user_setpoint_vars = []
        user_sp_weights = []
        for var, val in setpoint:
            user_setpoint.append(val)
            user_setpoint_vars.append(var)
            loc = locator[var].location
            user_sp_weights.append(locator[var].group.weights[loc])
        controller._NMPC_NAMESPACE.user_setpoint_weights = user_sp_weights
        controller._NMPC_NAMESPACE.user_setpoint = user_setpoint
        controller._NMPC_NAMESPACE.user_setpoint_vars = user_setpoint_vars

        # Add an objective function that only involves variables at t0
        self.add_objective_function(controller,
                control_penalty_type=ControlPenaltyType.ERROR,
                name=user_objective_name,
                objective_state_categories=[
                    VariableCategory.DIFFERENTIAL,
                    VariableCategory.ALGEBRAIC,
                    ],
                time_resolution_option=TimeResolutionOption.INITIAL_POINT)
        temp_objective = getattr(controller._NMPC_NAMESPACE, user_objective_name)

        self.solve_setpoint(
                categories=categories,
                require_steady=require_steady,
                **kwargs)

        # Deactivate objective that was just created
        temp_objective.deactivate()

        # Transfer setpoint values and reset initial values
        for categ in categories:
            vargroup = category_dict[categ]
            if categ == VariableCategory.SCALAR:
                for i, var in enumerate(vargroup):
                    vargroup.set_setpoint(i, var.value)
                    var.set_value(vargroup.reference[i])
            else:
                for i, var in enumerate(vargroup):
                    vargroup.set_setpoint(i, var[t0].value)
                    var[t0].set_value(vargroup.reference[i])


    def solve_setpoint(self, 
            categories = [VariableCategory.DIFFERENTIAL,
                          VariableCategory.ALGEBRAIC,
                          VariableCategory.DERIVATIVE,
                          VariableCategory.INPUT,
                          VariableCategory.SCALAR],
            require_steady=True,    
            **kwargs):
        config = self.config(kwargs)
        solver = config.solver
        outlvl = config.outlvl
        init_log = idaeslog.getInitLogger('nmpc', outlvl)
        solver_log = idaeslog.getSolveLogger('nmpc', outlvl)
        controller = self.controller
        time = self.controller_time
        t0 = time.first()
        namespace = getattr(controller, self.namespace_name)
        category_dict = namespace.category_dict

        was_originally_active = ComponentMap([(comp, comp.active) for comp in 
                controller.component_data_objects((Constraint, Block))])
        non_initial_time = list(time)[1:]
        deactivated = deactivate_model_at(controller, time, non_initial_time, outlvl)
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

        # Solve single-time point optimization problem
        dof = degrees_of_freedom(controller)
        if require_steady:
            assert dof == namespace.n_input_vars
        else:
            assert dof == namespace.n_input_vars + namespace.n_diff_vars
        init_log.info('Solving for full-state setpoint values')
        with idaeslog.solver_log(solver_log, level=idaeslog.DEBUG) as slc:
            results = solver.solve(controller, tee=slc.tee)
        if results.solver.termination_condition == TerminationCondition.optimal:
            init_log.info(
                    'Successfully solved for full state setpoint values')
        else:
            msg = 'Failed to solve for full state setpoint values'
            init_log.error(msg)
            raise RuntimeError(msg)

        # Revert changes. Again, order matters
        if require_steady == True:
            for var in category_dict[VariableCategory.DERIVATIVE]:
                var[t0].unfix()
        for var in controller._NMPC_NAMESPACE.ic_vars:
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


    def add_setpoint_to_controller(self, objective_name='tracking_objective',
            **kwargs):
        """User-facing function for the addition of a setpoint to the 
        controller. Assumes the controller model's setpoint attributes have
        been populated with desired values. This function first calculates
        weights, then adds an objective function based on those weights 
        and existing setpoint values.

        Args:
            objective_name : Name to use for the objective function added
        """
        # TODO: allow user to specify a steady state to use without having
        # called create_steady_state_setpoint
        config = self.config(kwargs)
        weight_override = config.objective_weight_override
        weight_tolerance = config.objective_weight_tolerance
        objective_state_categories = config.objective_state_categories
        time_resolution_option = config.time_resolution_option
        outlvl = config.outlvl
        #objective_name = config.full_state_objective_name

        self.construct_objective_weights(self.controller,
                objective_weight_override=weight_override,
                objective_weight_tolerance=weight_tolerance,
                categories=[
                    VariableCategory.DIFFERENTIAL,
                    VariableCategory.ALGEBRAIC,
                    VariableCategory.DERIVATIVE,
                    VariableCategory.INPUT,
                    ])

        # TODO: set point changes.
        self.add_objective_function(self.controller,
#                control_penalty_type=ControlPenaltyType.ACTION,
# NOTE: Leaving this commented here in case this breaks something
                control_penalty_type=config.control_penalty_type,    
                objective_state_categories=objective_state_categories,
                time_resolution_option=time_resolution_option,
                name=objective_name)


    def set_reference_values_from_initial(self, vargroup, t0=None):
        """Sets the values in the reference list of an NMPCVarGroup from the
        values of the group's variables at t0

        Args:
            vargroup : NMPCVarGroup instance whose reference values to set
            t0 : Point in time at which variable values will be used to set
                 reference values

        """
        if vargroup.is_scalar:
            raise ValueError(
                'No way to get initial conditions for a scalar component')
        else:
            if t0 is None:
                t0 = vargroup.t0
        for i in range(vargroup.n_vars):
            vargroup.reference[i] = vargroup.varlist[i][t0].value


    def construct_objective_weights(self, model,
            categories=[VariableCategory.DIFFERENTIAL,
                        VariableCategory.ALGEBRAIC,
                        VariableCategory.DERIVATIVE,
                        VariableCategory.INPUT], 
            **kwargs):
        """Constructs the objective weight values for the specified variable
        categories of a specified model. Weights are calculated for each 
        variable in each group by taking the difference between the initial
        value and the setpoint value, making sure it is above a tolerance,
        and taking its reciprocal. Weights can be overridden by a list
        of VarData, value tuples passed in as the "objective_weight_override"
        config argument.

        Args:
            model : Model whose variables will be accessed to calculate weights,
                    and whose weight attributes will be set.
            categories : List of VariableCategory enum items for which to 
                         calculate weights. Default is DIFFERENTIAL, ALGEBRAIC,
                         DERIVATIVE, and INPUT

        """
        config = self.config(kwargs)
        override = config.objective_weight_override
        tol = config.objective_weight_tolerance

        # Variables to override must be VarData objects in the model
        # for whose objective function we are calculating weights
        category_dict = model._NMPC_NAMESPACE.category_dict
        
        weights_to_override = {}
        for ow_tpl in override:
            locator = model._NMPC_NAMESPACE.var_locator[ow_tpl[0]]
            weights_to_override[(locator.category, locator.location)] = \
                    ow_tpl[1]

        # Given a vardata here, need to know its location so I know which
        # weight to override

        # Attempt to construct weight for each type of setpoint

        for categ in categories:
            vargroup = category_dict[categ]
            reference = vargroup.reference
            setpoint = vargroup.setpoint
            weights = vargroup.weights
            # construct the diagonal matrix (list).
            for loc, sp_value in enumerate(setpoint):
    
                # This assumes the vardata in sp is the same one
                # provided by the user. But these could differ by time
                # index...
                # Need to check by location, category here
                if (categ, loc) in weights_to_override:
                    weights[loc] = weights_to_override[categ, loc]
                    continue
    
                # If value is None, but variable was provided as override,
                # weight can still be non-None. This is okay.
                if sp_value is None or reference[loc] is None:
                    weights[loc] = None
                    continue
    
                diff = abs(reference[loc] - sp_value)
                if diff > tol:
                    weight = 1./diff
                else:
                    weight = 1./tol
                weights[loc] = weight
    

    def add_objective_function(self, model, name='objective', state_weight=1,
            control_weight=1, **kwargs):
        """Adds an objective function based on already calculated weights
        and setpoint values to the _NMPC_NAMESPACE of a model.

        Args:
            model : Model to which to add objective function
            name : Name of objective function to add
            state_weight : Additional weight factor to apply to each state
                           term in the objective function. Intended for a user
                           that wants to weigh states and controls differently
            control_weight : Addtional weight factor to apply to each control
                             term in the objective function. Intended for a user
                             that wants to weigh states and controls differently

        """
        config = self.config(kwargs)
        outlvl = config.outlvl
        init_log = idaeslog.getInitLogger('nmpc', level=outlvl)
        time_resolution = config.time_resolution_option
        state_categories = config.objective_state_categories

        # Q and R are p.s.d. matrices that weigh the state and
        # control norms in the objective function
        Q_diagonal = config.state_objective_weight_matrix_diagonal
        R_diagonal = config.control_objective_weight_matrix_diagonal

        # User may want to penalize control action, i.e. ||u_i - u_{i-1}||,
        # or control error (from set point), i.e. ||u_i - u*||
        # Valid values are ACTION or ERROR
        control_penalty_type = config.control_penalty_type
        if not (control_penalty_type == ControlPenaltyType.ERROR or
                control_penalty_type == ControlPenaltyType.ACTION or
                control_penalty_type == ControlPenaltyType.NONE):
            raise ValueError(
                "control_penalty_type argument must be 'ACTION' or 'ERROR'")

        if not Q_diagonal or not R_diagonal:
            raise NotImplementedError('Q and R must be diagonal for now.')
        
        category_dict = model._NMPC_NAMESPACE.category_dict 
        states = []
        Q_entries = []
        sp_states = []
        for categ in state_categories:
            if (categ == VariableCategory.INPUT and 
                control_penalty_type != ControlPenaltyType.NONE):
                raise ValueError(
        '''INPUT variable cannot be penalized as both states and controls.
        Either set control_penalty_type to ControlPenaltyType.NONE or
        omit VariableCategory.INPUT from objective_state_categories.'''
        )
            vargroup = category_dict[categ]
            states += vargroup.varlist
            Q_entries += vargroup.weights
            sp_states += vargroup.setpoint

        input_group = category_dict[VariableCategory.INPUT]
        controls = input_group.varlist
        R_entries = input_group.weights
        sp_controls = input_group.setpoint

        mod_time = model._NMPC_NAMESPACE.get_time()
        t0 = mod_time.first()
        # NOTE: t0 is now omitted from objective function, unless
        # INITIAL_POINT option is used
        if time_resolution == TimeResolutionOption.COLLOCATION_POINTS:
            time = [t for t in mod_time if t != mod_time.first()]
        if time_resolution == TimeResolutionOption.FINITE_ELEMENTS:
            time = [t for t in mod_time.get_finite_elements() 
                    if t != mod_time.first()]
        if time_resolution == TimeResolutionOption.SAMPLE_POINTS:
            sample_time = self.sample_time
            time = [t for t in model._NMPC_NAMESPACE.sample_points
                    if t != mod_time.first()]
        if time_resolution == TimeResolutionOption.INITIAL_POINT:
            time = [t0]

        state_term = sum(Q_entries[i]*(states[i][t] - sp_states[i])**2
                for i in range(len(states)) if (Q_entries[i] is not None
                                            and sp_states[i] is not None)
                         for t in time)
        # TODO: With what time resolution should states/controls be penalized?
        #       I think they should be penalized every sample point

        if control_penalty_type == ControlPenaltyType.ERROR:
            control_term = sum(R_entries[i]*(controls[i][t] - sp_controls[i])**2
                    for i in range(len(controls)) if (R_entries[i] is not None
                                                and sp_controls[i] is not None)
                                                for t in time)
        elif control_penalty_type == ControlPenaltyType.ACTION:
            # Override time list to be the list of sample points,
            # as these are the only points control action can be 
            # nonzero
            action_time = model._NMPC_NAMESPACE.sample_points[1:]
            time_len = len(action_time)
            if time_len == 1:
                init_log.warning(
                        'Warning: Control action penalty specfied '
                        'for a model with a single time point.'
                        'Control term in objective function will be empty.')
            control_term = sum(
                    R_entries[i]*(controls[i][action_time[k]] - 
                    controls[i][action_time[k-1]])**2
                    for i in range(len(controls)) 
                    if (R_entries[i] is not None
                    and sp_controls[i] is not None)
                    for k in range(1, time_len)
                    )
        elif control_penalty_type == ControlPenaltyType.NONE:
            control_term = 0
            # Note: This term is only non-zero at the boundary between sampling
            # times. Could use this info to make the expression more compact

        obj_expr = state_term + control_term

        # TODO: Deactivate existing objectives
        obj = Objective(expr=obj_expr)
        model._NMPC_NAMESPACE.add_component(name, obj)


    def set_bounds_from_initial(self, vargroup):
        """
        Builds lists of lower bound, upper bound tuples as attributes of the 
        input model, based on the current bounds (and domains) of
        differential, algebraic, and input variables.

        Args:
            model : Model whose variables will be checked for bounds.

        """

        varlist = vargroup.varlist
        if not vargroup.is_scalar:
            t0 = vargroup.index_set.first()
        for i, var in enumerate(varlist):
            if not vargroup.is_scalar:
                # Just assume these (t0) are the bounds/domain I want
                lb = var[t0].lb
                ub = var[t0].ub
                domain = var[t0].domain
            else:
                lb = var.lb
                ub = var.ub
                domain = var.domain
            if (domain == NonNegativeReals and lb is None):
                lb = 0
            elif (domain == NonNegativeReals and lb < 0):
                lb = 0
            vargroup.set_lb(i, lb)
            vargroup.set_ub(i, ub)


    def constrain_control_inputs_piecewise_constant(self,
            **kwargs):
        """Function to add piecewise constant (PWC) constraints to controller
        model. Requires model's _NMPC_NAMESPACE to know about input vars
        and to have as an attribute a sample points list.

        """
        config = self.config(kwargs)
        sample_time = config.sample_time
        outlvl = config.outlvl
        init_log = idaeslog.getInitLogger('nmpc', outlvl)
        init_log.info('Adding piecewise-constant constraints')

        model = self.controller

        # If sample_time is overwritten here, assume that the 
        # provided sample_time should be used going forward
        # (in input injection, plant simulation, and controller initialization)
        if sample_time != self.config.sample_time:
            self.validate_sample_time(sample_time, 
                    self.controller, self.plant)
            self.config.sample_time = sample_time

        time = model._NMPC_NAMESPACE.get_time()

        # This rule will not be picklable as it is not declared
        # at module namespace
        # Can access sample_time as attribute of namespace block,
        # then rule can be located outside of class
        input_indices = [i for i in range(model._NMPC_NAMESPACE.input_vars.n_vars)]
        def pwc_rule(ns, t, i):
            # Unless t is at the boundary of a sample, require
            # input[t] == input[t_next]
            time = ns.get_time()
            if t in ns.sample_points or t == time.first():
                return Constraint.Skip
            t_next = time.next(t)
            inputs = ns.input_vars.varlist
            _slice = inputs[i]
            return _slice[t_next] == _slice[t]

        name = 'pwc_constraint'
        pwc_constraint = Constraint(time, input_indices, 
                rule=pwc_rule)
        model._NMPC_NAMESPACE.add_component(name, pwc_constraint)

        pwc_constraint_list = [Reference(pwc_constraint[:, i])
                           for i in input_indices]
        model._NMPC_NAMESPACE.pwc_constraint_list = pwc_constraint_list


    def initialize_control_problem(self, **kwargs):
        """Function to initialize the controller model before solving the
        optimal control problem. Possible strategies are to use the initial
        conditions, to perform a simulation, or to use the results of the 
        previous solve. Initialization from a previous (optimization)
        solve can only be done if an optimization solve has been performed
        since the last initialization. The strategy may be passed in as
        the control_init_option keyword (config) argument, otherwise the
        default will be used.

        """
        config = self.config(kwargs)
        strategy = config.control_init_option
        solver = config.solver

        input_type = config.element_initialization_input_option

        time = self.controller_time

        if strategy == ControlInitOption.FROM_PREVIOUS:
            self.initialize_from_previous_sample(self.controller, **kwargs)

        elif strategy == ControlInitOption.BY_TIME_ELEMENT:
            self.initialize_by_solving_elements(self.controller, self.controller_time,
                    input_type=input_type, **kwargs)

        elif strategy == ControlInitOption.FROM_INITIAL_CONDITIONS:
            self.initialize_from_initial_conditions(self.controller, **kwargs)
        
        # Add check that initialization did not violate bounds/equalities?

        self.controller_solved = False


    def initialize_by_solving_elements(self, model, time,
            input_type=ElementInitializationInputOption.SET_POINT,
            objective_name='tracking_objective',
            **kwargs):
        """Initializes the controller model by solving (a square simulation
        for) each time element.

        Args:
            model : Model to initialize
            time : Set to treat as time
            input_type : ElementInitializationInputOption enum item 
                         telling how to fix the inputs for the simulation

        """
        config = self.config(kwargs)
        tol = config.tolerance
        outlvl = config.outlvl
        objective = getattr(model._NMPC_NAMESPACE, 
                            objective_name)
        namespace = model._NMPC_NAMESPACE

        # Strip bounds before simulation as square solves will be performed
        strip_controller_bounds = TransformationFactory(
                                      'contrib.strip_var_bounds')
        strip_controller_bounds.apply_to(model, reversible=True)

        input_vars = model._NMPC_NAMESPACE.input_vars
        if input_type == ElementInitializationInputOption.SET_POINT:
            for i, _slice in enumerate(input_vars.varlist):
                for t in time:
                    if t != time.first():
                        _slice[t].fix(input_vars.setpoint[i])
                    else:
                        _slice[t].fix()
        elif input_type == ElementInitializationInputOption.INITIAL:
            for i, _slice in enumerate(input_vars.varlist):
                t0 = time.first()
                for t in time:
                    _slice[t].fix(_slice[t0].value)
        else:
            raise ValueError('Unrecognized input option')
        # The above should ensure that all inputs are fixed and the 
        # model has no dof upon simulation

        # Deactivate objective function
        # Here I assume the name of the objective function.
        # TODO: ObjectiveType Enum and objective_dict
        objective.deactivate()
        model._NMPC_NAMESPACE.pwc_constraint.deactivate()

        initialize_by_element_in_range(self.controller, self.controller_time, 
                    time.first(), 
                    time.last(),
                    outlvl=outlvl,
                    dae_vars=self.controller._NMPC_NAMESPACE.dae_vars,
                    time_linking_variables=self.controller._NMPC_NAMESPACE.diff_vars)

        objective.activate()
        model._NMPC_NAMESPACE.pwc_constraint.activate()

        for _slice in self.controller._NMPC_NAMESPACE.input_vars:
            for t in time:
                if t != time.first():
                    # Don't want to unfix inputs at time.first()
                    _slice[t].unfix()

        strip_controller_bounds.revert(self.controller)

        timelist = list(time)
        for cat, group in namespace.category_dict.items():
            if (cat == VariableCategory.FIXED or cat == VariableCategory.INPUT
                    or cat == VariableCategory.SCALAR):
                continue
            violated = get_violated_bounds_at_time(group, timelist, tol)
            if violated:
                raise ValueError(
                    'Bounds violated after solving elements: %s'
                    % str(violated))


    def initialize_from_previous_sample(self, model,
            categories=[VariableCategory.DIFFERENTIAL,
                        VariableCategory.ALGEBRAIC,
                        VariableCategory.DERIVATIVE,
                        VariableCategory.INPUT],
            **kwargs):
        """Re-initializes values of variables in model to the values one 
        sampling time in the future. Values for the last sampling time are 
        currently set to values in the steady state model, assumed to be the 
        set point.

        Args:
            model : Flowsheet model to initialize
            categories : List of VariableCategory enum items to initialize.
                         Default contains DIFFERENTIAL, ALGEBRAIC, DERIVATIVE,
                         and INPUT

        """
        # Should only do this if controller is initialized
        # from a prior solve.
        if not self.controller_solved:
            raise RuntimeError(
                    'Cannot initialize from previous if the control '
                    'problem has not previously been solved.'
                    )

        config = self.config(kwargs)
        sample_time = config.sample_time
        tolerance = config.continuous_set_tolerance

        # TODO
        # Should initialize dual variables here too.

        time = model._NMPC_NAMESPACE.get_time()
        category_dict = model._NMPC_NAMESPACE.category_dict
        # TODO: have some attribute for steady time 
        # Or better yet, don't use a steady_model at all

        for categ in categories:
            varlist = category_dict[categ].varlist
            for i, _slice in enumerate(varlist):
                for t in time:
                    # If not in last sample:
                    if (time.last() - t) >= sample_time:
                        t_next = find_point_in_continuousset(
                                t + sample_time,
                                time, tolerance=tolerance)
                        _slice[t].set_value(_slice[t_next].value)
                    else:
                        _slice[t].set_value(category_dict[categ].setpoint[i])


    def initialize_from_initial_conditions(self, model, 
            categories=[VariableCategory.DERIVATIVE,
                        VariableCategory.DIFFERENTIAL,
                        VariableCategory.ALGEBRAIC],
            **kwargs):
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
        config = self.config(kwargs)
        time = model._NMPC_NAMESPACE.get_time()
        cat_dict = model._NMPC_NAMESPACE.category_dict
        for categ in categories:
            varlist = cat_dict[categ].varlist
            for v in varlist:
                v[:].set_value(v[0].value)

    
    def solve_control_problem(self, **kwargs):
        """Function for solving optimal control problem, which calculates
        control inputs for the plant.

        """
        config = self.config(kwargs)
        solver = config.solver
        outlvl = config.outlvl
        init_log = idaeslog.getInitLogger('nmpc', level=outlvl)
        s_log = idaeslog.getSolveLogger('nmpc', level=outlvl)

        time = self.controller_time
        for _slice in self.controller._NMPC_NAMESPACE.input_vars:
            for t in time:
                if t == time.first():
                    _slice[t].fix()
                else:
                    _slice[t].unfix()

        assert (degrees_of_freedom(self.controller) == 
                self.controller._NMPC_NAMESPACE.n_input_vars*
                (self.controller._NMPC_NAMESPACE.samples_per_horizon))

        with idaeslog.solver_log(s_log, idaeslog.DEBUG) as slc:
            results = solver.solve(self.controller, tee=slc.tee)
        if results.solver.termination_condition == TerminationCondition.optimal:
            init_log.info('Successfully solved optimal control problem')
            self.controller_solved = True
        else:
            init_log.error('Failed to solve optimal control problem')
            raise ValueError


    def simulate_plant(self, t_start, **kwargs):
        """Function for simulating plant model for one sampling period after
        inputs have been assigned from solve of controller model.

        Args:
            t_start : Beginning of timespan over which to simulate

        """
        config = self.config(kwargs)

        sample_time = self.config.sample_time
        # ^ Use self.config here, as I don't want user to override sample_time
        #   at this point. How to throw an error if they do? - use immutable param
        # TODO
        calculate_error = config.calculate_error
        outlvl = config.outlvl
        init_log = idaeslog.getInitLogger('nmpc', level=outlvl)
        tol = config.continuous_set_tolerance

        t_end = t_start + sample_time 
        assert t_start in self.plant_time
        t_end = find_point_in_continuousset(t_end, self.plant_time, tol)
        assert t_end in self.plant_time

        initialize_by_element_in_range(self.plant, self.plant_time, t_start, t_end, 
                dae_vars=self.plant._NMPC_NAMESPACE.dae_vars, 
                time_linking_vars=self.plant._NMPC_NAMESPACE.diff_vars,
                outlvl=outlvl)

        self.current_plant_time = t_end

        msg = ('Successfully simulated plant over the sampling period '
                'through ' + str(self.current_plant_time))
        init_log.info(msg)

        tc1 = self.controller_time.first() + sample_time

        if self.controller_solved and calculate_error:
            self.state_error[t_end] = self.calculate_error_between_states(
                    self.controller, self.plant, tc1, t_end)


    def calculate_error_between_states(self, mod1, mod2, t1, t2, 
            Q_matrix=[],
            categories=[VariableCategory.DIFFERENTIAL],
            **kwargs):
        """
        Calculates the normalized (by the weighting matrix already calculated)
        error between the differential variables in different models and at
        different points in time.

        Args:
            mod1 : First flowsheet model
            mod2 : Second flowsheet model (may be same as the first)
            t1 : Time point of interest in first model
            t2 : Time point of interest in second model
            Q_matrix : List of weights by which to weigh the error for
                       each state. Default is to use the same weights calculated
                       for the controller objective function.

        """
        config = self.config(kwargs)

        Q_diagonal = config.state_objective_weight_matrix_diagonal
        if not Q_diagonal:
            raise ValueError('Only diagonal weighting matrices are supported')
        # Grab the weighting matrix from the controller model regardless of what
        # mod1 and mod2 are. This can be overwritten if desired.

        # TODO: allow option to override weights
        # As the default, weights are taken from model 1

        # Used to specify variables other than differential to use for
        # error calculation
        
        varlist_1 = []
        varlist_2 = []

        weight_matrix_provided = bool(Q_matrix)
        for categ in categories:
            varlist_1 += mod1._NMPC_NAMESPACE.category_dict[categ].varlist
            varlist_2 += mod2._NMPC_NAMESPACE.category_dict[categ].varlist
            if not weight_matrix_provided:
                Q_matrix += self.controller._NMPC_NAMESPACE.category_dict[categ].weights
        assert len(varlist_1) == len(varlist_2)
        n = len(varlist_1)

        assert t1 in mod1._NMPC_NAMESPACE.get_time()
        assert t2 in mod2._NMPC_NAMESPACE.get_time()

        error = sum(Q_matrix[i]*(varlist_1[i][t1].value - 
                                 varlist_2[i][t2].value)**2
                    for i in range(n) if Q_matrix[i] is not None)

        return error

