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
from idaes.apps.caprese.model import (
        DynamicModelHelper,
        ControllerHelper,
        )
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
#    CONFIG = DynamicBase.CONFIG
    CONFIG = ConfigDict()

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
                default=ElementInitializationInputOption.SETPOINT,
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

        self.plant = DynamicModelHelper(
                plant_model,
                plant_time_set,
                inputs_at_t0,
                )

        init_controller_inputs = [
                find_comp_in_block(controller_model, plant_model, comp)
                for comp in inputs_at_t0
                ]
        self.controller = ControllerHelper(
                controller_model,
                controller_time_set,
                init_controller_inputs,
                )

        controller_measurements = self.controller.measured_vars
        plant_measurements = self.plant.find_components(
                self.controller.model,
                controller_measurements,
                self.controller.time,
                )

        plant_inputs = self.plant.category_dict[VariableCategory.INPUT]
        controller_inputs = self.controller.find_components(
                self.plant.model,
                plant_inputs,
                self.plant.time,
                )

        self.controller_measurement_map = ComponentMap(
                zip(plant_measurements, controller_measurements))

        self.plant_input_map = ComponentMap(
                zip(controller_inputs, plant_inputs))

        self.controller.validate_sample_time(sample_time)
        self.plant.validate_sample_time(sample_time)
        self.sample_time = sample_time

        # NOTE: This is probably unnecessary:
#        self.validate_fixedness(self.plant, self.controller)

        self.current_plant_time = 0

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

    def inject_inputs_into_plant(self):
        time = self.controller.time
        ts = time.first() + self.sample_time
        for p_var, c_var in self.plant_input_map.items():
            # TODO: Could apply noise here, or could do it elsewhere
            p_var[:].fix(c_var[ts].value)

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

