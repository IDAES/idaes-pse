# -*- coding: UTF-8 -*-
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

from pyomo.environ import (Block, Constraint, Var, TerminationCondition,
        SolverFactory, Objective, NonNegativeReals, Reals, 
        TransformationFactory, Reference)
from pyomo.core.base.var import _GeneralVarData
from pyomo.kernel import ComponentSet, ComponentMap
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.dae.flatten import flatten_dae_variables
from pyomo.opt.solver import SystemCallSolver
from pyutilib.misc.config import ConfigDict, ConfigValue

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.dyn_utils import (get_activity_dict, deactivate_model_at,
        path_from_block, find_comp_in_block, find_comp_in_block_at_time)
from idaes.core.util.initialization import initialize_by_time_element
from idaes.dynamic.caprese.util import (initialize_by_element_in_range,
        find_slices_in_model, NMPCVarLocator, copy_values_at_time, 
        add_noise_at_time, ElementInitializationInputOption, 
        TimeResolutionOption, ControlInitOption, ControlPenaltyType,
        VariableCategory, validate_list_of_vardata, 
        validate_list_of_vardata_value_tuples, validate_solver,
        NMPCVarGroup)
import idaes.logger as idaeslog

from collections import OrderedDict
import time as timemodule
import enum
import pdb

__author__ = "Robert Parker and David Thierry"


# TODO tags:
# - KWARGS
# - CONFIG
# - NAMESPACE
# - RENAME
# - REMOVECONSTANT
# - GENERALIZE


class NMPCSim(object):
    """
    Main class for NMPC simulations of IDAES flowsheets. 
    """
    CONFIG = ConfigDict()

    # What is the right way to configure logger outlvl?
    # Should I use a default (Init/Solver) IDAES logger or make my own?
    CONFIG.declare(
            'outlvl',
            ConfigValue(
                default=idaeslog.INFO, 
                doc='Output verbosity level for IDAES logger'
                )
            )
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
            'continuous_set_tolerance',
            ConfigValue(
                default=1e-8,
                domain=float,
                doc=('Tolerance used for determining whether a float is a '
                    'member of a ContinuousSet')
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
    # These args, I just declare once each... perhaps because I don't expect
    # them to be overwritten?
    CONFIG.declare(
            'noise_weights',
            ConfigValue(
                default=[],
                domain=list,
                doc=('List of weights to override weights for variance '
                    'in noise function')
                )
            )
    # ^ Really this should be a list of vardata, value tuples
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
    # Still unsure of the purpose of Config.
    # - Is it a way to pass arguments with default values and validation?
    #   If so, should every argument have a ConfigValue?
    # - Is it a way to set default (or constructor-level) arguments, then
    #   safely and selectively override them if desired during the algorithm?
    #   This seems like the most practical use to me.

    # How should I handle solver config? What is a good domain?
    # Since each solver could have a different interface for options,
    # I'll refrain from adding a ConfigValue for solver options.
    # Plus having a "config"-inside-a-config seems kind of wonky.
    CONFIG.declare('solver',
            ConfigValue(
                default=SolverFactory('ipopt'),
                domain=validate_solver,
                doc='Pyomo solver object to be used to solve generated NLPs'
                )
            )
    CONFIG.declare('tolerance',
            ConfigValue(
                default=1e-8,
                domain=float,
                doc='Tolerance for checking constraint violation'
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
    # If this is intended to be declared only in each function (so it can
    # be declared with variables from the right model), does it make sense
    # to declare it here? I don't really want the user passing in values at the
    # constructor level, but I do want the default value...
    CONFIG.declare('objective_weight_override',
            ConfigValue(
                default=[],
                domain=validate_list_of_vardata_value_tuples,
                doc=('User-specified objective weight values for given '
                    'variables that take precedence over calculated values')
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
    

    def __init__(self, plant_model=None, plant_time_set=None, 
        controller_model=None, controller_time_set=None, inputs_at_t0=None,
        sample_time=None, outlvl=idaeslog.NOTSET, solver=None):
        """Constructor method. Accepts plant and controller models needed for 
        NMPC simulation, as well as inputs at the first time point in the 
        plant model. Models provided are added to the self as attributes.
        This constructor solves for consistent initial conditions 
        in the plant and controller and performs categorization into lists of
        differential, derivative, algebraic, input, fixed, and scalar variables,
        which are added as attributes to the provided models.

        Args:
            plant_model : Plant flowsheet model, NMPC of which will be 
                          simulated. Currently this must contain the entire 
                          timespan it is desired to simulate.
            controller_model : Model to be used to calculate control inputs
                               for the plant. Control inputs in controller
                               must exist in the plant, and initial condition
                               variables in the plant must exist in the 
                               controller.
            initial_inputs : List of VarData objects containing the variables
                             to be treated as control inputs, at t = 0.

        Kwargs:
            solver : Solver to be used for verification of consistent initial 
                     conditions, will also be used as the default solver if
                     another is not provided for initializing or solving the 
                     optimal control problem.
            outlvl : IDAES logger output level. Default is idaes.logger.NOTSET.
                     To see solver traces, use idaes.logger.DEBUG.
            sample_time : Length of time each control input will be held for.
                          This must be an integer multiple of the (finite
                          element) discretization spacing in both the plant
                          and controller models. Default is to use the 
                          controller model's discretization spacing.
        """
        self.config = self.CONFIG
        self.config.sample_time = sample_time
        self.config.inputs_at_t0 = inputs_at_t0
        self.config.outlvl = outlvl
        self.config.solver = solver
        # Maybe include a kwarg for require_steady - if False, set-point is not
        # forced to be a steady state

        # TODO: handle arguments (solver, outlvl, ...) via config blocks
        #       separate solvers for each solve, in init
        # KWARGS, CONFIG, NAMESPACE - remove kwargs and set these up here

#        solver = kwargs.pop('solver', SolverFactory('ipopt'))
#        self.default_solver = solver

        # TODO: validate_time_set function

        # Logger properties:
        # outlvl added as an attribute so it can be referenced later
#        self.outlvl = kwargs.pop('outlvl', idaeslog.NOTSET)
#        self.outlvl = outlvl

        init_log = idaeslog.getInitLogger('nmpc', level=self.config.outlvl)

        # Set up attributes
        self.p_mod = plant_model
        self.p_mod_time = plant_time_set
        self.c_mod = controller_model
        self.c_mod_time = controller_time_set

        self.add_namespace_to(self.p_mod, self.p_mod_time)
        self.add_namespace_to(self.c_mod, self.c_mod_time)

        # Validate models
        self.validate_models(self.p_mod, self.c_mod)

        # Solve for consistent initial conditions

        # Categorize variables in plant model
        init_log.info('Categorizing variables in plant model') 
        self.categorize_variables(self.p_mod, self.p_mod_time, inputs_at_t0)

        self.solve_initial_conditions(self.p_mod, self.p_mod_time)
        # TODO: move into own function
        #       possibly a DAE utility
        #       - check for consistency of initial conditions
        #       - if not consistent, tell user to go run dae.solve_initial_conditions

        # TODO: - add attributes as private "_initial_inputs" - make block private
        #       - check that I don't override anything
        #       - store in a block: util.modeling.unique_component_name
        #       - cleanup if desired
        self.p_mod.initial_inputs = inputs_at_t0
        # ^ Add these as an attribute so they can be used later to validate
        # steady state model
        # ^ What did I mean by this? In what way does steady state model need
        # to resemble plant model?
        # Not sure this is still necessary

        self.p_mod._NMPC_NAMESPACE.category_dict = {
                VariableCategory.DIFFERENTIAL:
                        self.p_mod._NMPC_NAMESPACE.diff_vars,
                VariableCategory.DERIVATIVE:
                        self.p_mod._NMPC_NAMESPACE.deriv_vars,
                VariableCategory.ALGEBRAIC:
                        self.p_mod._NMPC_NAMESPACE.alg_vars,
                VariableCategory.INPUT: 
                        self.p_mod._NMPC_NAMESPACE.input_vars,
                VariableCategory.FIXED: 
                        self.p_mod._NMPC_NAMESPACE.fixed_vars,
                VariableCategory.SCALAR: 
                        self.p_mod._NMPC_NAMESPACE.scalar_vars,
                }
        self.build_variable_locator(self.p_mod, self.p_mod_time, 
                self.p_mod._NMPC_NAMESPACE.category_dict,
                ic_vars=self.p_mod._NMPC_NAMESPACE.ic_vars)
        # Now adding a locator to the plant model so I can find plant model
        # variables corresponding to the controller's initial conditions

        # Categorize variables in controller model
        init_controller_inputs = self.validate_initial_inputs(self.c_mod,
                self.p_mod, self.p_mod_time, self.p_mod.initial_inputs)
        init_log.info('Categorizing variables in the controller model')
        self.categorize_variables(self.c_mod, self.c_mod_time, init_controller_inputs)
        # TODO: Would like to use my variable category enum here
        # As simple as replacing 'differential' by 'DIFFERENTIAL', etc.?
        # Could construct a dict: enum -> var list, but that might be more work than necessary?
        # Then constructing ComponentMap is much easier...
        self.c_mod._NMPC_NAMESPACE.category_dict = {
                VariableCategory.DIFFERENTIAL:
                        self.c_mod._NMPC_NAMESPACE.diff_vars,
                VariableCategory.DERIVATIVE:
                        self.c_mod._NMPC_NAMESPACE.deriv_vars,
                VariableCategory.ALGEBRAIC:
                        self.c_mod._NMPC_NAMESPACE.alg_vars,
                VariableCategory.INPUT:
                        self.c_mod._NMPC_NAMESPACE.input_vars,
                VariableCategory.FIXED:
                        self.c_mod._NMPC_NAMESPACE.fixed_vars,
                VariableCategory.SCALAR:
                        self.c_mod._NMPC_NAMESPACE.scalar_vars,
                }
        self.build_variable_locator(self.c_mod, self.c_mod_time,
                self.c_mod._NMPC_NAMESPACE.category_dict,
                ic_vars=self.c_mod._NMPC_NAMESPACE.ic_vars)

        # Only need to manipulate bounds of controller model. Assume the 
        # bounds in the plant model should remain in place for simulation.
        # (Should probably raise a warning if bounds are present...)
        for categ, vargroup in self.c_mod._NMPC_NAMESPACE.category_dict.items():
            self.set_bounds_from_initial(vargroup)
        # ^ This may be removed in favor of strip_bounds transformation

        # Validate inputs in the plant model and initial conditions
        # in the control model.
        # TODO: allow user to specify this if names don't match
        self.p_mod._NMPC_NAMESPACE.controller_ic_vars = find_slices_in_model(
                self.p_mod,
                self.c_mod,
                self.p_mod._NMPC_NAMESPACE.var_locator,
                self.c_mod._NMPC_NAMESPACE.ic_vars)
        self.c_mod._NMPC_NAMESPACE.plant_input_vars = find_slices_in_model(
                self.c_mod,
                self.p_mod,
                self.c_mod._NMPC_NAMESPACE.var_locator,
                self.p_mod._NMPC_NAMESPACE.input_vars.varlist)

        self.validate_fixedness(self.p_mod, self.p_mod_time)
        self.validate_fixedness(self.c_mod, self.c_mod_time)

# TODO: remove. Place in solve_initial_conditions method if it exists.
#       If desired ('strict mode') check for consistency.
#####################
        copy_values_at_time(self.c_mod._NMPC_NAMESPACE.ic_vars,
                            self.p_mod._NMPC_NAMESPACE.controller_ic_vars,
                            self.c_mod_time.first(),
                            self.p_mod_time.first())
        # May eventually want to "load_initial_conditions" here so I can apply 
        # noise, but for now this is more direct.

        # Should strip bounds before this IC solve, since the controller
        # model should have bounds
        self.strip_controller_bounds = TransformationFactory(
                                       'contrib.strip_var_bounds')
        self.strip_controller_bounds.apply_to(self.c_mod, reversible=True)

        # Controller model has already been categorized... No need 
        # to provide init_controller_inputs
        self.solve_initial_conditions(self.c_mod, self.c_mod_time)
        # ^ move into initialize_control_problem
        self.strip_controller_bounds.revert(self.c_mod)
#####################

#        self.sample_time = kwargs.pop('sample_time', 
#                self.c_mod.time.get_finite_elements()[1] -
#                    self.c_mod.time.get_finite_elements()[0])
        self.sample_time = self.config.sample_time
        self.validate_sample_time(self.sample_time, 
                self.c_mod, self.p_mod)

        # TODO: validate_discretization_scheme option?
        scheme = self.c_mod.time.get_discretization_info()['scheme']
        if scheme == 'LAGRANGE-RADAU':
            self.c_mod._NMPC_NAMESPACE.ncp = \
                    self.c_mod.time.get_discretization_info()['ncp']
        elif scheme == 'BACKWARD Difference':
            self.c_mod._NMPC_NAMESPACE.ncp = 1
        else:
            raise NotImplementedError

        scheme = self.p_mod.time.get_discretization_info()['scheme']
        if scheme == 'LAGRANGE-RADAU':
            self.p_mod._NMPC_NAMESPACE.ncp = \
                    self.p_mod.time.get_discretization_info()['ncp']
        elif scheme == 'BACKWARD Difference':
            self.p_mod._NMPC_NAMESPACE.ncp = 1
        else:
            raise NotImplementedError

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


    def add_namespace_to(self, model, time):
        name = '_NMPC_NAMESPACE'
        # Not _CAPRESE_NAMESPACE as I might want to add a similar 
        # namespace for MHE
        if hasattr(model, name):
            raise ValueError('%s already exists on model. Please fix this.'
                             % name)
        model.add_component(name, Block())
        namespace = getattr(model, name)
        assert time.model() == model.model()

        def get_time():
            return time
        # Will not be bound. Not sure if that matters
        namespace.get_time = get_time


    def validate_sample_time(self, sample_time, *models):
        """Makes sure sample time is an integer multiple of discretization
        spacing in each model, and that horizon of each model is an integer
        multiple of sample time.

        Args:
            sample_time: Sample time to check
            models: List of flowsheet models to check
        """
        for model in models:
            time = model._NMPC_NAMESPACE.get_time()
            fe_spacing = (time.get_finite_elements()[1] -
                          time.get_finite_elements()[0])
            if fe_spacing > sample_time:
                raise ValueError(
                    'Sampling time must be at least as long '
                    'as a finite element in time')
            # TODO: how to handle roundoff here?
            if (sample_time / fe_spacing) % 1 != 0:
                raise ValueError(
                    'Sampling time must be an integer multiple of '
                    'finite element spacing for every model')
            else:
                model._NMPC_NAMESPACE.fe_per_sample = \
                        int(sample_time/fe_spacing)

            horizon_length = time.last() - time.first()
            if (horizon_length / sample_time) % 1 != 0:
                raise ValueError(
                    'Sampling time must be an integer divider of '
                    'horizon length')
            else:
                model._NMPC_NAMESPACE.samples_per_horizon = \
                        int(horizon_length/sample_time)


    def validate_slices(self, tgt_model, src_model, src_slices):
        """
        Given list of time-only slices in a source model, attempts to find
        each of them in the target model and returns a list of the found 
        slices in the same order.
        Expects to find a var_locator dictionary attribute in the target
        model.

        Args:
            tgt_model : Model to search for time-slices
            src_model : Model containing the slices to search for
            src_slices : List of time-only slices of variables in the source
                         model

        Returns:
            List of time-only slices to same-named variables in the target 
            model
        """
        # What I need to do is actually a little more complicated...
        # Because I can't just find the slice that I want in the target model.
        # That slice doesn't exist in the model, it's essentially just a nice
        # literary device for iterating over some related variables.
        # So for every slice provided, extract the first VarData, find it in
        # the target model, then use the NMPCVarLocator to find the corresponding
        # slice in the target model...
        t0 = src_model.time.first()
        tgt_slices = []
        locator = tgt_model._NMPC_NAMESPACE.var_locator
        for _slice in src_slices:
            init_vardata = _slice[t0]
            tgt_vardata = find_comp_in_block(tgt_model, 
                                             src_model, 
                                             init_vardata)
            tgt_container = locator[tgt_vardata].group.varlist
            location = locator[tgt_vardata].location
            tgt_slices.append(tgt_container[location])
        return tgt_slices


    def validate_fixedness(self, model, time):
        """
        Makes sure that assumptions regarding fixedness for different points
        in time are valid. Differential, algebraic, and derivative variables
        may be fixed only at t0, only if they are initial conditions.
        Fixed variables must be fixed at all points in time, except possibly
        initial conditions.
        """
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
                    

    # TODO: option to skip this step by user specification of input pairs
    def validate_initial_inputs(self, tgt_model, src_model, src_time, 
            src_inputs=None, **kwargs):
        """Uses initial inputs in the source model to find variables of the
        same name in a target model.
        
        Args:
           tgt_model : Flowsheet model to search for input variables
           src_model : Flowsheet model containing inputs to search for
           src_inputs : List of input variables at the initial time point
                        to find in target model. If not provided, the
                        initial_inputs attribute will be used.
        """
        config = self.config(kwargs)
        outlvl = config.outlvl
        
        log = idaeslog.getInitLogger('nmpc', level=outlvl)

        # src_time only necessary to find inputs if not provided?

        if not src_inputs:
            # If source inputs are not specified, assume they
            # already exist in src_model
            try:
                t0 = src_time.first()
                src_inputs = [v[t0] for v in src_model.input_vars]
            except AttributeError:
                msg = ('Error validating inputs. Either provide src_inputs '
                      'or categorize_inputs in the source model first')
                idaeslog.error(msg)
                raise

        tgt_inputs = []
        for inp in src_inputs:
            local_parent = tgt_model
            for r in path_from_block(inp, src_model, include_comp=True):
                try:
                    local_parent = getattr(local_parent, r[0])[r[1]]
                except AttributeError:
                    msg = (f'Error validating input {inp.name}.'
                           'Could not find component {r[0]} in block '
                           '{local_parent.name}.')
                    log.error(msg)
                    raise
                except KeyError:
                    msg = (f'Error validating {inp.name}.'
                           'Could not find key {r[1]} in component '
                           '{getattr(local_parent, r[0]).name}.')
                    log.error(msg)
                    raise
            tgt_inputs.append(local_parent)
        return tgt_inputs


    def validate_models(self, m1, m2):
        """
        Makes sure the two models are FlowsheetBlocks and are
        distinct.

        Args:
            m1 : First flowsheet model
            m2 : Second flowsheet model
        """
        if not (isinstance(m1, Block) and
                isinstance(m2, Block)):
            raise ValueError(
                    'Provided models must be Blocks')
        if m1.model() is m2.model():
            raise ValueError(
                    'Provided models must not live in the same top-level'
                    'ConcreteModel')
        return True


    def transfer_current_plant_state_to_controller(self, t_plant, 
            config=ConfigDict(), **kwargs):
        # Would like to pass "noise_args" in as a bundle here. This can
        # probably be done with config blocks somehow.

        # TESTME: this function is not tested.
        # copy_values and add_noise are tested, however
        config = self.config(config)

        # appropriate to hard-code c_mod_time here
        time = self.c_mod_time
        t0 = time.first()

        copy_values_at_time(self.c_mod._NMPC_NAMESPACE.ic_vars,
                            self.p_mod._NMPC_NAMESPACE.controller_ic_vars,
                            t0,
                            t_plant)

        # Apply noise to new initial conditions
        add_noise = config.add_plant_noise

        # TODO: config for noise args
        noise_weights = config.noise_weights
        noise_sig_0 = config.noise_sigma_0
        noise_args = config.noise_arguments
        max_noise_weight = config.max_noise_weight

        locator = self.c_mod._NMPC_NAMESPACE.var_locator
        if add_noise:
            if not noise_weights:
                noise_weights = []
                for var in self.c_mod._NMPC_NAMESPACE.ic_vars:
                    info = locator[var[t0]]
#                    category = info.category
                    loc = info.location
                    # TODO: Should be a dict mapping: category -> correct weights
                    #       var_info[category].weights (set point, bounds, etc.)
                    #       This would be moot if I could map directly from 
                    #       variable to weights
                    obj_weight = info.group.weights[loc]
#                    if category == VariableCategory.DIFFERENTIAL:
#                        obj_weight = self.c_mod._NMPC_NAMESPACE.diff_weights[loc]
#                    elif category == VariableCategory.DERIVATIVE:
#                        # Is it okay to weigh derivative noise by 
#                        # the state variable's weight?
#                        obj_weight = self.c_mod._NMPC_NAMESPACE.diff_weights[loc]
#                    else:
#                        # FIXME: Won't have obj weight for algebraic equations
#                        # NAMESPACE: empty list should be added automatically in
#                        # namespace block
#                        # Really should calculate these weights, even if they
#                        # aren't used
#                        obj_weight = None

                    if obj_weight is not None and obj_weight != 0:
                        noise_weights.append(min(1/obj_weight, 
                                                 max_noise_weight))
                    else:
                        noise_weights.append(None)

            add_noise_at_time(self.c_mod._NMPC_NAMESPACE.ic_vars,
                              t0,
                              weights=noise_weights,
                              sigma_0=noise_sig_0,
                              **noise_args)


    def inject_control_inputs_into_plant(self, t_plant, config=ConfigDict(), 
            **kwargs):
        """Injects input variables from the first sampling time in the 
        controller model to the sampling period in the plant model that
        starts at the specified time.

        Args:
            t_plant : First time point in plant model where inputs will be
                      applied.
            
        Kwargs:
            sample_time : Sample time in the plant over which the inputs
                          will be applied. Default is to use the sample
                          time assigned by the constructor or overwritten
                          by the creation of PWC constraints.
            add_noise : Bool telling whether or not to apply noise
            noise_weights : List of weights for each state's standard deviation
            noise_sig_0 : Standard deviation for a state with unit weight
            noise_args : Additional kwargs to pass apply_noise_at_time

        """
        # KWARGS, CONFIG - handle all of these through config. Pass in noise args
        # as a bundle
        # config args for control_input_noise
        config = self.config(config)

        sample_time = self.config.sample_time

        # Send inputs to plant that were calculated for the end
        # of the first sample
        t_controller = self.c_mod_time.first() + sample_time
        assert t_controller in self.c_mod_time

        time = self.p_mod_time
        plant_sample = [t for t in time if t > t_plant and t<= t_plant+sample_time]
        assert t_plant+sample_time in plant_sample
        # len(plant_sample) should be ncp*nfe_per_sample, assuming the expected
        # sample_time is passed in
        # Should have some logic to correct if assert fails due to a rounding
        # error.

        add_noise = config.add_input_noise
        noise_weights = config.noise_weights
        noise_sig_0 = config.noise_sigma_0
        noise_args = config.noise_arguments
        max_noise_weight = config.max_noise_weight

        # Need to get proper weights for plant's input vars
        locator = self.c_mod._NMPC_NAMESPACE.var_locator
        if add_noise:
            if not noise_weights:
                noise_weights = []
                for var in self.c_mod._NMPC_NAMESPACE.plant_input_vars:
                    info = locator[var[t_controller]]
#                    category = locator.category
                    loc = info.location
#                    # TODO: allow other variable categories here?
#                    # GENERALIZE
#                    if category == VariableCategory.INPUT:
#                        obj_weight = self.c_mod._NMPC_NAMESPACE.input_weights[loc]
#                    elif category == VariableCategory.DIFFERENTIAL:
#                        obj_weight = self.c_mod._NMPC_NAMESPACE.diff_weights[loc]
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

            add_noise_at_time(self.c_mod._NMPC_NAMESPACE.plant_input_vars,
                              t_controller,
                              weights=noise_weights,
                              sigma_0=noise_sig_0,
                              **noise_args)
            #add_noise_at_time(self.p_mod.input_vars,
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

        copy_values_at_time(self.p_mod._NMPC_NAMESPACE.input_vars.varlist,
                            self.c_mod._NMPC_NAMESPACE.plant_input_vars,
                            plant_sample,
                            t_controller)


    def solve_initial_conditions(self, model, time, **kwargs):
        """Function to solve for consistent initial conditions in
        the provided flowsheet model.

        Args:
            model : Flowsheet model whose initial conditions are solved

        Kwargs:
            solver : Solver object to use
            outlvl : idaes.logger output level
        """
        # Record which Constraints/Variables are initially inactive
        # deactivate model except at t=0
        # fix initial inputs - raise error if no value
        # raise error if not square
        # re-activate model
        #
        # Later include option to skip solve for consistent initial conditions
        #
        # Will only work as written for "True" initial conditions since 
        # it doesn't try to deactivate discretization equations or fix
        # derivative/differential variables.

        config = self.config(kwargs)

        # TODO: Yes. Should activity_dict be a ComponentMap?
        was_originally_active = get_activity_dict(model)
        solver = config.solver
        outlvl = config.outlvl
        init_log = idaeslog.getInitLogger('nmpc', level=outlvl)
        solver_log = idaeslog.getSolveLogger('nmpc', level=outlvl)

        toplevel = model.model()
        t0 = time.first()

        non_initial_time = [t for t in time]
        non_initial_time.remove(t0)
        deactivated = deactivate_model_at(model, time, non_initial_time,
                outlvl=idaeslog.ERROR)

        # Proper to specify input_vars here
        for _slice in model._NMPC_NAMESPACE.input_vars.varlist:
            vardata = _slice[t0]
            if vardata.model() is not toplevel:
                raise ValueError(
                        f"Trying to fix an input that does not belong to model"
                        " {toplevel.name}. Are 'initial_inputs' arguments"
                        " contained in the proper model?")
            if vardata.value == None:
                raise ValueError(
                        "Trying to solve for consistent initial conditions with "
                        "an input of value None. Inputs must have valid values "
                        "to solve for consistent initial conditions")
            else:
                vardata.fix()

        if degrees_of_freedom(model) != 0:
            raise ValueError(
                    f"Non-zero degrees of freedom in model {toplevel.name}"
                    " when solving for consistent initial conditions. Have the "
                    "right number of initial conditions been fixed?")

        with idaeslog.solver_log(solver_log, level=idaeslog.DEBUG) as slc:
            results = solver.solve(model, tee=slc.tee)

        if results.solver.termination_condition == TerminationCondition.optimal:
            init_log.info(
                    'Successfully solved for consistent initial conditions')
        else:
            init_log.error('Failed to solve for consistent initial conditions')
            raise ValueError(
            f'Falied to solve {toplevel.name} for consistent initial conditions')

        if was_originally_active is not None:
            for t in non_initial_time:
                for comp in deactivated[t]:
                    if was_originally_active[id(comp)]:
                        comp.activate()


    def categorize_variables(self, model, time, initial_inputs):
        """Function to create lists of time-only-slices of the different 
        types of variables in a model, given knowledge of which are inputs. 
        These lists are added as attributes to the model.

        Args:
            model : Model whose variables will be flattened and categorized
            initial_inputs : List of VarData objects that are input variables
                             at the initial time point

        """

        # Would probably be much faster to do this categorization
        # during variable flattening
        #time = model.time

        # Works for steady state models as time will be an ordered
        # (although not continuous) set:
        t0 = time.first()
        try:
            t1 = time.get_finite_elements()[1]
        except AttributeError:
            t1 = t0

        # TODO: subblock
        deriv_vars = []
        diff_vars = []
        input_vars = []
        alg_vars = []
        fixed_vars = []

        ic_vars = []

        # Create list of time-only-slices of time indexed variables
        # (And list of VarData objects for scalar variables)
        scalar_vars, dae_vars = flatten_dae_variables(model, time)

        dae_map = ComponentMap([(v[t0], v) for v in dae_vars])
        t0_vardata = list(dae_map.keys())
        #model.dae_vars = list(dae_map.values())
        model._NMPC_NAMESPACE.dae_vars = list(dae_map.values())
        #model.scalar_vars = list(ComponentMap([(v, v) for v in scalar_vars]).values())
        model._NMPC_NAMESPACE.scalar_vars = \
            NMPCVarGroup(
                list(ComponentMap([(v, v) for v in scalar_vars]).values()),
                index_set=None, is_scalar=True)
        model._NMPC_NAMESPACE.n_scalar_vars = \
                model._NMPC_NAMESPACE.scalar_vars.n_vars
        input_set = ComponentSet(initial_inputs)
        updated_input_set = ComponentSet(initial_inputs)
        diff_set = ComponentSet()

        # Iterate over initial vardata, popping from dae map when an input,
        # derivative, or differential var is found.
        for var0 in t0_vardata:
            if var0 in updated_input_set:
                input_set.remove(var0)
                time_slice = dae_map.pop(var0)
                input_vars.append(time_slice)
             
            parent = var0.parent_component()
            if not isinstance(parent, DerivativeVar):
                continue
            if not time in ComponentSet(parent.get_continuousset_list()):
                continue
            index0 = var0.index()
            var1 = dae_map[var0][t1]
            index1 = var1.index()
            state = parent.get_state_var()

            if state[index1].fixed:
                # Assume state var is fixed everywhere, so derivative
                # 'isn't really' a derivative.
                # Should be safe to remove state from dae_map here
                state_slice = dae_map.pop(state[index0])
                fixed_vars.append(state_slice)
                continue
            if state[index0] in input_set:
                # If differential variable is an input, then this DerivativeVar
                # is 'not really a derivative'
                continue

            deriv_slice = dae_map.pop(var0)
            if var1.fixed:
                # Assume derivative has been fixed everywhere.
                # Add to list of fixed variables, and don't remove its state variable.
                fixed_vars.append(deriv_slice)
            elif var0.fixed:
                # In this case the derivative has been used as an initial condition. 
                # Still want to include it in the list of derivatives.
                ic_vars.append(deriv_slice)
                state_slice = dae_map.pop(state[index0])
                if state[index0].fixed:
                    ic_vars.append(state_slice)
                deriv_vars.append(deriv_slice)
                diff_vars.append(state_slice)
            else:
                # Neither is fixed. This should be the most common case.
                state_slice = dae_map.pop(state[index0])
                if state[index0].fixed:
                    ic_vars.append(state_slice)
                deriv_vars.append(deriv_slice)
                diff_vars.append(state_slice)

        if not updated_input_set:
            raise RuntimeError('Not all inputs could be found')
        assert len(deriv_vars) == len(diff_vars)

        for var0, time_slice in dae_map.items():
            var1 = time_slice[t1]
            # If the variable is still in the list of time-indexed vars,
            # it must either be fixed (not a var) or be an algebraic var
            if var1.fixed:
                fixed_vars.append(time_slice)
            else:
                if var0.fixed:
                    ic_vars.append(time_slice)
                alg_vars.append(time_slice)

        model._NMPC_NAMESPACE.deriv_vars = NMPCVarGroup(deriv_vars, time)
        model._NMPC_NAMESPACE.diff_vars = NMPCVarGroup(diff_vars, time)
        model._NMPC_NAMESPACE.n_diff_vars = len(diff_vars)
        model._NMPC_NAMESPACE.n_deriv_vars = len(deriv_vars)
        assert (model._NMPC_NAMESPACE.n_diff_vars == 
                model._NMPC_NAMESPACE.n_deriv_vars)
                
        # ic_vars will not be stored as a NMPCVarGroup - don't want to store
        # all the info twice
        model._NMPC_NAMESPACE.ic_vars = ic_vars
        model._NMPC_NAMESPACE.n_ic_vars = len(ic_vars)
        #assert model.n_dv == len(ic_vars)
        # Would like this to be true, but accurately detecting differential
        # variables that are not implicitly fixed (by fixing some input)
        # is difficult

        model._NMPC_NAMESPACE.input_vars = NMPCVarGroup(input_vars, time)
        model._NMPC_NAMESPACE.n_input_vars = len(input_vars)

        model._NMPC_NAMESPACE.alg_vars = NMPCVarGroup(alg_vars, time)
        model._NMPC_NAMESPACE.n_alg_vars = len(alg_vars)

        model._NMPC_NAMESPACE.fixed_vars = NMPCVarGroup(fixed_vars, time)
        model._NMPC_NAMESPACE.n_fixed_vars = len(fixed_vars)


    def build_variable_locator(self, model, time, category_dict,
            ic_vars=[]):
        """Constructs a dictionary mapping the id of each VarData object
        to a NMPCVarLocator object. This dictionary is added as an attribute to
        the model.

        Args:
            model : Flowsheet model containing the variables provided

        Kwargs:
            algebraic : List of algebraic variable time-slices
            differential : List of differential variable time-slices
            derivative : List of derivative variable time-slices
            input : List of input variable time-slices
            fixed : List of fixed variable time-slices
            scalar : List of non-time-indexed variables
            ic : List of differential, algebraic, or derivative variables
                 that will be fixed as initial conditions
        """
        ic_list = ic_vars

        locator = ComponentMap()
        for categ, vargroup in category_dict.items():
            varlist = vargroup.varlist
            if categ == VariableCategory.SCALAR:
                for i, var in enumerate(varlist):
                    locator[var] = NMPCVarLocator(categ, vargroup, i)
            else:
                for i, var in enumerate(varlist):
                    for t in time:
                        locator[var[t]] = NMPCVarLocator(categ,
                                vargroup, i)

        # Since these variables already have NMPCVarLocator objects,
        # just set the desired attribute.
        for i, _slice in enumerate(ic_list):
            for t in model.time:
                locator[_slice[t]].is_ic = True

        model._NMPC_NAMESPACE.var_locator = locator


    def solve_steady_state_setpoint(self, set_point, steady_model, **kwargs):
        config = self.config(kwargs)
        outlvl = config.outlvl
        weight_override = config.objective_weight_override
        weight_tolerance = config.objective_weight_tolerance

        if not steady_model:
            raise ValueError(
               "'steady_model' required to validate set point")
        self.s_mod = steady_model
        self.add_namespace_to(self.s_mod, steady_model.time)
        self.validate_models(self.s_mod, self.p_mod)
        self.validate_steady_setpoint(set_point, self.s_mod,
                                      outlvl=outlvl,
                                      objective_weight_override=weight_override,
                                      objective_weight_tolerance=weight_tolerance)
        # ^ result should be that controller now has set point attributes


    def add_setpoint_to_controller(self, **kwargs):
        """User-facing function for the addition of a set point to the 
        controller.
        """
        # TODO: allow user to specify a steady state to use without having
        # called create_steady_state_setpoint
        config = self.config(kwargs)
        weight_override = self.config.objective_weight_override
        weight_tolerance = self.config.objective_weight_tolerance
        outlvl = config.outlvl

        self.construct_objective_weights(self.c_mod,
                objective_weight_override=weight_override,
                objective_weight_tolerance=weight_tolerance)

        self.add_objective_function(self.c_mod,
                control_penalty_type=ControlPenaltyType.ACTION,
                name='tracking_objective')


    def validate_steady_setpoint(self, set_point, steady_model, **kwargs):

        config = self.config(kwargs)
        solver = config.solver
        outlvl = config.outlvl
        
        init_log = idaeslog.getInitLogger('nmpc', level=outlvl)
        solver_log = idaeslog.getSolveLogger('nmpc', level=outlvl)
        weight_override = config.objective_weight_override
        weight_tolerance = config.objective_weight_tolerance

        # The following loop will create steady-state variable lists in
        # proper order, initialize steady state model to initial conditions
        # of controller model, and make sure that the fixed status of 
        # variables is the same across these models.

        # Assume that variable with time-index t0 exists in steady model
        t0 = self.c_mod_time.first()
        # Compared fixed-ness as t1 so we're not thrown off by fixed
        # initial conditions
        t1 = self.c_mod_time.get_finite_elements()[1]

        steady_cat_dict = {}
        for categ, vargroup in self.c_mod._NMPC_NAMESPACE.category_dict.items():
            if categ == VariableCategory.DERIVATIVE:
                continue
            varlist = []
            for _slice in vargroup.varlist:
                # TODO: need some way to know which category a variable is
                # (ComponentMap?)
                if not vargroup.is_scalar:
                    vardata_t0 = _slice[t0]
                    vardata_t1 = _slice[t1]
                else:
                    vardata_t0 = _slice
                    vardata_t1 = _slice

                local_parent = steady_model
                # TODO: replace with helper function
                # Is there a UID function I could use to get this?
                for r in path_from_block(vardata_t0, self.c_mod,
                                         include_comp=True):
                    try:
                        local_parent = getattr(local_parent, r[0])[r[1]]
                    except AttributeError:
                        init_log.error(
                             f'Error initializing steady state model: '
                             'Could not find {r[0]} in {local_parent.name}. '
                             'Was the steady state model constructed with '
                             'has_holdups=True?')
                        raise
                    except KeyError:
                        init_log.error(
                            f'KeyError while initializing steady state model. '
                            'Was steady state model constructed with same '
                            'spacial discretization and component/phase lists '
                            '(if applicable)? Does time start at 0 in '
                            'controller model?')
                        raise

                var_steady = local_parent
                var_steady.set_value(vardata_t0.value)
                
                if not vargroup.is_scalar:
                    varlist.append({t0: var_steady})
                    # Append dict to list here so that steady state
                    # VarDatas for "DAE vars" are accessed with the
                    # same syntax as in the dynamic case.
                else:
                    varlist.append(var_steady)

                # Copy fixed status from controller model at t1
                if not var_steady.fixed == vardata_t1.fixed:
                    var_steady.fixed = vardata_t1.fixed

            steady_cat_dict[categ] = NMPCVarGroup(varlist,
                    steady_model.time, is_scalar=vargroup.is_scalar)
        steady_model._NMPC_NAMESPACE.category_dict = steady_cat_dict
        steady_model._NMPC_NAMESPACE.diff_vars = \
                steady_cat_dict[VariableCategory.DIFFERENTIAL].varlist
        steady_model._NMPC_NAMESPACE.alg_vars = \
                steady_cat_dict[VariableCategory.ALGEBRAIC].varlist
        steady_model._NMPC_NAMESPACE.input_vars = \
                steady_cat_dict[VariableCategory.INPUT].varlist
        steady_model._NMPC_NAMESPACE.fixed_vars = \
                steady_cat_dict[VariableCategory.FIXED].varlist
        steady_model._NMPC_NAMESPACE.scalar_vars = \
                steady_cat_dict[VariableCategory.SCALAR].varlist

        # TODO: these don't exist yet
        assert (len(steady_model._NMPC_NAMESPACE.diff_vars) == 
                self.c_mod._NMPC_NAMESPACE.n_diff_vars)
        assert (len(steady_model._NMPC_NAMESPACE.alg_vars) == 
                self.c_mod._NMPC_NAMESPACE.n_alg_vars)
        assert (len(steady_model._NMPC_NAMESPACE.input_vars) == 
                self.c_mod._NMPC_NAMESPACE.n_input_vars)
        assert (len(steady_model._NMPC_NAMESPACE.fixed_vars) == 
                self.c_mod._NMPC_NAMESPACE.n_fixed_vars)
        assert (len(steady_model._NMPC_NAMESPACE.scalar_vars) == 
                self.c_mod._NMPC_NAMESPACE.n_scalar_vars)

        # Now that steady state model has been validated and initialized,
        # construct complete lists for the *_sp attributes
# TODO: these are already None
#        diff_var_sp = [None for var in steady_model._NMPC_NAMESPACE.diff_vars]
#        alg_var_sp = [None for var in steady_model._NMPC_NAMESPACE.alg_vars]
#        input_var_sp = [None for var in steady_model._NMPC_NAMESPACE.input_vars]
        # ^ list entries are not actually vars, they are dictionaries...
        # maybe this is too confusing

        # This is where I map user values for set points into lists that
        # I can use to build the objective function (and weight matrices).
        for vardata, value in set_point:
            # set_point variables should be members of controller model
            info = self.c_mod._NMPC_NAMESPACE.var_locator[vardata]
            category = info.category
            location = info.location
            group = info.group

            # Only allow set point variables from following categories:
            if (category != VariableCategory.DIFFERENTIAL and
                category != VariableCategory.ALGEBRAIC and
                category != VariableCategory.INPUT):
                raise ValueError('Set point variables (data objects) must be'
                                 'differential, algebraic, or inputs.')
            # Maybe this restriction should stand, but I should probably map
            # category to the corresponding set-point list
            group.set_setpoint(location, value)

        # Add attributes (with proper names!) to steady model
        # so these values can be accessed later
# These attributes no longer exist and the info has been properly recorded
#        steady_model._NMPC_NAMESPACE.diff_setpoints = diff_var_sp
#        steady_model._NMPC_NAMESPACE.alg_setpoints = alg_var_sp
#        steady_model._NMPC_NAMESPACE.input_setpoints = input_var_sp

# Has already been created.
#        category_dict = {
#        VariableCategory.DIFFERENTIAL: steady_model._NMPC_NAMESPACE.diff_vars,
#        VariableCategory.ALGEBRAIC: steady_model._NMPC_NAMESPACE.alg_vars,
#        VariableCategory.INPUT: steady_model._NMPC_NAMESPACE.input_vars,
#        VariableCategory.FIXED: steady_model._NMPC_NAMESPACE.fixed_vars,
#        VariableCategory.SCALAR: steady_model._NMPC_NAMESPACE.scalar_vars,
#        }
        self.build_variable_locator(steady_model, 
                steady_model._NMPC_NAMESPACE.get_time(), steady_cat_dict)

        # Set values of reference variables
        for categ, group in steady_cat_dict.items():
            if categ != VariableCategory.SCALAR:
                self.set_reference_values_from_initial(group)
        for categ, group in self.c_mod._NMPC_NAMESPACE.category_dict.items():
            if categ != VariableCategory.SCALAR:
                self.set_reference_values_from_initial(group)

        self.construct_objective_weights(steady_model,
                categories=[VariableCategory.ALGEBRAIC,
                           VariableCategory.DIFFERENTIAL,
                           VariableCategory.INPUT],
                objective_weight_override=weight_override,
                objective_weight_tolerance=weight_tolerance)

        # Add objective to steady_model
        # Add bounds to steady_model (taken from control model)
        #
        # Make sure inputs are unfixed - validate degrees of freedom
        #
        # solve steady model for set point
        # add sp attributes to control model

        self.add_objective_function(steady_model,
                            state_categories=[VariableCategory.DIFFERENTIAL,
                                              VariableCategory.ALGEBRAIC],
                            control_penalty_type=ControlPenaltyType.ERROR,
                            name='user_setpoint_objective')

        # Transfer bounds to steady state model
        for categ, vargroup in steady_cat_dict.items():
            controller_group = self.c_mod._NMPC_NAMESPACE.category_dict[categ]
            self.transfer_bounds(vargroup, controller_group)

        # Unfix inputs for solve
        for var in steady_model._NMPC_NAMESPACE.input_vars:
            for t in steady_model.time:
                var[t].unfix()

        # Verify proper number of degrees of freedom
        # TODO: abstract this into a verify_dof function
        assert (degrees_of_freedom(steady_model) ==
                len(steady_model._NMPC_NAMESPACE.input_vars)*len(steady_model.time))

        # Solve steady state model 
        init_log.info('Solving for steady state set-point.')
        with idaeslog.solver_log(solver_log, level=idaeslog.DEBUG) as slc:
            results = solver.solve(steady_model, tee=slc.tee)
        if results.solver.termination_condition == TerminationCondition.optimal:
            init_log.info(
                    'Successfully solved for steady state set-point')
        else:
            msg = 'Failed to solve for steady state setpoint'
            init_log.error(msg)
            raise RuntimeError(msg)

        for categ, steady_group in steady_cat_dict.items():
            controller_group = self.c_mod._NMPC_NAMESPACE.category_dict[categ]
            for i, var in enumerate(steady_group):
                controller_group.setpoint[i] = var[t0].value
        # TODO: Change if economic NMPC (solve for "consistent derivatives")
        for i in range(self.c_mod._NMPC_NAMESPACE.deriv_vars.n_vars):
            self.c_mod._NMPC_NAMESPACE.deriv_vars.setpoint[i] = 0


    def set_reference_values_from_initial(self, vargroup, t0=None):
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
        """
        Do I even need the model? 
        ^ it makes things slightly easier to provide only the model, I think...
        Need model to get the var_locator

        Do I want to allow user to provide a setpoint here?
        ^ No, but want to allow them to override weights if they want

        Here each setpoint is a list of (VarData, value) tuples.

        Overwrite is a lists of (VarData, value) tuples, where 
        this value will directly override the weight of the
        corresponding VarData. (For all time..., what if user
        provides multiple weights for VarDatas that only differ
        by time? Raise warning and use last weight provided)
        """
        config = self.config(kwargs)
#        override = kwargs.pop('override', [])
        override = config.objective_weight_override

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

        tol = config.objective_weight_tolerance

        # Need t0 to get the 'reference' value to compare with set point value
        # for weight calculation
#        t0 = vargroup.index_set.first()
# list of initial values should have already been populated

        # Attempt to construct weight for each type of setpoint
        # that we could allow

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
                if sp_value is None:
                    weights[loc] = None
                    continue
    
                # This line works for steady state variables thanks to
                # duck typing
                diff = abs(reference[loc] - sp_value)
                if diff > tol:
                    weight = 1./diff
                else:
                    weight = 1./tol
                weights[loc] = weight
    

    def add_objective_function(self, model, name='objective', state_weight=1,
            control_weight=1, 
            state_categories=[VariableCategory.DIFFERENTIAL], 
            **kwargs):
        """
        Assumes that model has already been populated with set point 
        and weights.
        Need to include state?
        Can't access ss vars in same manner as vars in dynamic model - 
        entries in varlists are not slices, they are already VarDatas...
        Solution would be to either add dynamic/ss flag, or to modify ss varlists
        to look like those dynamic. (Can't actually categorize because no derivs,
        but could flatten into slices, then assemble into lists based on names. 
        This seems like a lot of extra work though.)
        """
        config = self.config(kwargs)
        name = kwargs.pop('name', 'objective')
        outlvl = config.outlvl
        init_log = idaeslog.getInitLogger('nmpc', level=outlvl)
        time_resolution = config.time_resolution_option

        # Q and R are p.s.d. matrices that weigh the state and
        # control norms in the objective function
        Q_diagonal = config.state_objective_weight_matrix_diagonal
        R_diagonal = config.control_objective_weight_matrix_diagonal

        # User may want to penalize control action, i.e. ||u_i - u_{i-1}||,
        # or control error (from set point), i.e. ||u_i - u*||
        # Valid values are ACTION or ERROR
        control_penalty_type = config.control_penalty_type
        if not (control_penalty_type == ControlPenaltyType.ERROR or
                control_penalty_type == ControlPenaltyType.ACTION):
            raise ValueError(
                "control_penalty_type argument must be 'ACTION' or 'ERROR'")

        if not Q_diagonal or not R_diagonal:
            raise NotImplementedError('Q and R must be diagonal for now.')
        
        category_dict = model._NMPC_NAMESPACE.category_dict 
        states = []
        Q_entries = []
        sp_states = []
        for categ in state_categories:
            vargroup = category_dict[categ]
            states += vargroup.varlist
            Q_entries += vargroup.weights
            sp_states += vargroup.setpoint

        input_group = category_dict[VariableCategory.INPUT]
        controls = input_group.varlist
        R_entries = input_group.weights
        sp_controls = input_group.setpoint

        # TODO: Make this a list of time points, replace with list of
        # time-finite elements or sampling points depending on user preference.
        time = model.time

        state_term = sum(Q_entries[i]*(states[i][t] - sp_states[i])**2
                for i in range(len(states)) if (Q_entries[i] is not None
                                            and sp_states[i] is not None)
                         for t in time)
        # TODO: With what time resolution should states/controls be penalized?
        #       I think they should be penalized every sample point

        # Should be some check that Q_entries == None <=> sp_states == None
        # Not necessarily true, Q can have an entry due to an override

        if control_penalty_type == ControlPenaltyType.ERROR:
            control_term = sum(R_entries[i]*(controls[i][t] - sp_controls[i])**2
                    for i in range(len(controls)) if (R_entries[i] is not None
                                                and sp_controls[i] is not None)
                                                for t in time)
        elif control_penalty_type == ControlPenaltyType.ACTION:
            if len(time) == 1:
                init_log.warning(
                        'Warning: Control action penalty specfied '
                        'for a model with a single time point.'
                        'Control term in objective function will be empty.')
            control_term = sum(R_entries[i]*
                                  (controls[i][time[k]] - controls[i][time[k-1]])**2
                for i in range(len(controls)) if (R_entries[i] is not None
                                            and sp_controls[i] is not None)
                                            for k in range(2, len(time)+1))
        # TODO: allow a control penalty type of None to turn off this term.
            # Note: This term is only non-zero at the boundary between sampling
            # times. Could use this info to make the expression more compact
            # (TODO)
            # Populate time as either the list of finite element of the list
            # of sample points. Change time to be a list, not ContinuousSet

        obj_expr = state_term + control_term

        # TODO: namespace block
        obj = Objective(expr=obj_expr)
        model.add_component(name, obj)


    def set_bounds_from_initial(self, vargroup):
        """
        Builds lists of lower bound, upper bound tuples as attributes of the 
        input model, based on the current bounds (and domains) of
        differential, algebraic, and input variables.

        Args:
            model : Model whose variables will be checked for bounds.

        Returns:
            None
        """

#        for categ, vargroup in model._NMPC_NAMESPACE.category_dict.items():
        varlist = vargroup.varlist
        if not vargroup.is_scalar:
            t0 = vargroup.index_set.first()
        for i, var in enumerate(varlist):
            if not vargroup.is_scalar:
                # Just assume these are the bounds/domain I want
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


    def transfer_bounds(self, tgt_group, src_group):
        """
        Transfers bounds from source model's bound lists
        to target model's differential, algebraic, and input
        variables, and sets domain to Reals.

        Args:
            tgt_model : Model whose variables bounds will be transferred to
            src_model : Model whose bound lists will be used to set bounds.

        Returns:
            None
        """
        n_vars = tgt_group.n_vars
        for i in range(n_vars):
            tgt_group.set_lb(i, src_group.lb[i])
            tgt_group.set_ub(i, src_group.ub[i])
            tgt_group.set_domain(i, Reals)


    def constrain_control_inputs_piecewise_constant(self,
            **kwargs):
        """Function to add piecewise constant (PWC) constraints to
        model. Inputs and sample time are already known, so no arguments are
        necessary.

        Kwargs:
            model : Model to which PWC constraints are added. Default is   
                    controller model
            sample_time : Duration for which inputs will be forced constant
            outlvl : idaes.logger output level
        """
        config = self.config(kwargs)
        sample_time = config.sample_time
        outlvl = config.outlvl
        init_log = idaeslog.getInitLogger('nmpc', outlvl)
        init_log.info('Adding piecewise-constant constraints')

        model = self.c_mod

        # If sample_time is overwritten here, assume that the 
        # provided sample_time should be used going forward
        # (in input injection, plant simulation, and controller initialization)
        if sample_time != self.config.sample_time:
            self.validate_sample_time(sample_time, 
                    self.c_mod, self.p_mod)
            self.config.sample_time = sample_time

        time = model.time
        t0 = model.time.get_finite_elements()[0]
        t1 = model.time.get_finite_elements()[1]
        nfe_spacing = t1-t0
        nfe_per_sample = sample_time / nfe_spacing 
        if nfe_per_sample - int(nfe_per_sample) != 0:
            raise ValueError
        nfe_per_sample = int(nfe_per_sample)

        # This rule will not be picklable as it is not declared
        # at module namespace
        # Can access sample_time as attribute of namespace block,
        # then rule can be located outside of class
        input_indices = [i for i in range(model._NMPC_NAMESPACE.input_vars.n_vars)]
        def pwc_rule(m, t, i):
            # Unless t is at the boundary of a sample, require
            # input[t] == input[t_next]
            # TODO: Do this to tolerance
            # pyomo.core.base.range::Remainder
            # (or just math.remainder -- nope, added in 3.7)
            # if abs(remainder(t-time.first(), sample_time)) < 1e-8
            # set collocation_tolerance = 1/2(min spacing) during
            # construction, reference it here.
            if (t - time.first()) % sample_time == 0:
                return Constraint.Skip
            t_next = time.next(t)
            inputs = m._NMPC_NAMESPACE.input_vars.varlist
            _slice = inputs[i]
            return _slice[t_next] == _slice[t]

        name = 'pwc_input'
        pwc_constraint = Constraint(model.time, input_indices, 
                rule=pwc_rule)
        model.add_component(name, pwc_constraint)

        pwc_constraint_list = [Reference(pwc_constraint[:, i])
                           for i in input_indices]
        model._NMPC_NAMESPACE.pwc_constraint_list = pwc_constraint_list


    def initialize_control_problem(self, **kwargs):
        """Function to initialize the controller model before solving the
        optimal control problem. Possible strategies are to use the initial
        conditions, to perform a simulation, or to use the results of the 
        previous solve. Initialization from a previous (optimization)
        solve can only be done if an optimization solve has been performed
        since the last initialization.

        Kwargs:
            strategy : String describing the initialization strategy. Possible
                       values are 'from_previous', 'from_simulation', and
                       'initial_conditions'. Default is 'from_previous'.
            solver : Solver object to be used for initialization from simulation
            solver_options : Dictionary of options to pass to the solver
        """
        config = self.config(kwargs)
        strategy = config.control_init_option
        solver = config.solver

        input_type = config.element_initialization_input_option

        time = self.c_mod_time

        if strategy == ControlInitOption.FROM_PREVIOUS:
            self.initialize_from_previous_sample(self.c_mod)

        elif strategy == ControlInitOption.BY_TIME_ELEMENT:
            self.initialize_by_solving_elements(self.c_mod, self.c_mod_time,
                    input_type=element_initialization_input_option)

        elif strategy == ControlInitOption.FROM_INITIAL_CONDITIONS:
            self.initialize_from_initial_conditions(self.c_mod)
        
        # Add check that initialization did not violate bounds/equalities?

        self.controller_solved = False


    def initialize_by_solving_elements(self, model, time,
            input_type=ElementInitializationInputOption.SET_POINT):
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
        elif input_type == ElementInitializationInputOption.INITIAL_CONDITIONS:
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
        model._NMPC_NAMESPACE.tracking_objective.deactivate()
        model._NMPC_NAMESPACE.pwc_constraint.deactivate()

        initialize_by_element_in_range(self.c_mod, self.c_mod_time, 
                    time.first(), time.last(),
                    dae_vars=self.c_mod._NMPC_NAMESPACE.dae_vars,
                    time_linking_variables=self.c_mod._NMPC_NAMESPACE.diff_vars)

        # Reactivate objective, pwc constraints, bounds
        # TODO: safer objective name
        self.c_mod._NMPC_NAMESPACE.tracking_objective.activate()
        model.pwc_constraint.activate()

        for _slice in self.c_mod._NMPC_NAMESPACE.input_vars:
            for t in time:
                _slice[t].unfix()

        strip_controller_bounds.revert(self.c_mod)


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

        Kwargs: 
            sample_time : Length of time by which to shift variable values.
                          Default uses the sample time provided to the 
                          constructor or overwritten by the PWC constraints.
            attr_list : List of attribute names containing variables whose
                        values should be re-initialized. Default is 
                        'diff_vars', 'alg_vars', 'deriv_vars', and
                        'input_vars'.
        """
        # Should only do this if controller is initialized
        # from a prior solve.
        if not self.controller_solved:
            raise ValueError

        config = self.config(kwargs)
        sample_time = config.sample_time

        # TODO
        # Should initialize dual variables here too.

        time = model._NMPC_NAMESPACE.get_time()
        category_dict = model._NMPC_NAMESPACE.category_dict
        # TODO: have some attribute for steady time 
        # Or better yet, don't use a steady_model at all
        # Here I use steady model, assuming it holds the desired set point values
        # Really, I should access set point lists...

        for categ in categories:
            varlist = category_dict[categ].varlist
            for i, _slice in enumerate(varlist):
                for t in time:
                    # If not in last sample:
                    if (time.last() - t) >= sample_time:
                        t_next = t + sample_time

                        # Performing addition on CtsSet indices can result in
                        # rounding errors. Round to 8th decimal place here:
                        # TODO: config.continuous_set_tolerance
                        t_next = int(round(t_next*1e8))/1e8

                        assert t_next in time
                        _slice[t].set_value(_slice[t_next].value)

                    # Otherwise initialize to final steady state
                    # (Should provide some other options here,
                    # like copy inputs from set point then simulate)
                    else:
#                        if not attrname == 'deriv_vars':
                        _slice[t].set_value(category_dict[categ].setpoint[i])
#                        else:
#                            _slice[t].set_value(0)


    def initialize_from_initial_conditions(self, model, 
            categories=[VariableCategory.DERIVATIVE,
                        VariableCategory.DIFFERENTIAL,
                        VariableCategory.ALGEBRAIC]):
        """ 
        Set values of differential, algebraic, and derivative variables to
        their values at the initial conditions.
        An implicit assumption here is that the initial conditions are
        consistent.

        Args:
            model : Flowsheet model whose variables are initialized

        Kwargs:
            attr_list : List of names of attributes that contain variables
                        whose values should be initialized
        """
        # could use copy_values_at_time here
#        attr_list = kwargs.pop('attr_list', ['deriv_vars', 'diff_vars', 'alg_vars'])
        # TODO: user provides their own list of variables to initialize, or
        # I grab all these attributes from the model and add them up here.
        time = model._NMPC_NAMESPACE.get_time()
        cat_dict = model._NMPC_NAMESPACE.category_dict
        for categ in categories:
            varlist = cat_dict[categ].varlist
            for v in varlist:
                v[:].set_value(v[0].value)

    
    def solve_control_problem(self, **kwargs):
        """Function for solving optimal control problem, which calculates
        control inputs for the plant.

        Kwargs:
            solver : Solver object to be used, already loaded with user's
                     desired options. Default is that provided to the 
                     constructor.
            outlvl : idaes.logger output level. Default is that provided
                     to the constructor.
        """
        config = self.config(kwargs)
        solver = config.solver
        outlvl = config.outlvl
        init_log = idaeslog.getInitLogger('nmpc', level=outlvl)
        s_log = idaeslog.getSolveLogger('nmpc', level=outlvl)

        time = self.c_mod_time
        for _slice in self.c_mod._NMPC_NAMESPACE.input_vars:
            for t in time:
                _slice[t].unfix()
        # ^ Maybe should fix inputs at time.first()?
        # Also, this should be redundant as inputs have been unfixed
        # after initialization

        assert (degrees_of_freedom(self.c_mod) == 
                self.c_mod._NMPC_NAMESPACE.n_input_vars*
                (self.c_mod._NMPC_NAMESPACE.samples_per_horizon + 1))

        with idaeslog.solver_log(s_log, idaeslog.DEBUG) as slc:
            results = solver.solve(self.c_mod, tee=slc.tee)
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

        Kwargs:
            sample_time : Length of timespan to simulate. Default is the sample
                          time provided to the constructor or overwritten by
                          PWC constraints.
            outlvl : idaes.logger output level
        """
        config = self.config(kwargs)

        sample_time = self.config.sample_time
        # ^ Use self.config here, as I don't want user to override sample_time
        #   at this point. How to throw an error if they do? - use immutable param
        calculate_error = config.calculate_error
        outlvl = config.outlvl
        init_log = idaeslog.getInitLogger('nmpc', level=outlvl)

        t_end = t_start + sample_time 
        assert t_start in self.p_mod_time
        # TODO: add tolerance here, as t_end could change due to roundoff
        #       Then change t_end s.t. it is a point in time
        # A helper function may be useful - is_in_to_tolerance
        # Need to adjust t_end 
        assert t_end in self.p_mod_time

        initialize_by_element_in_range(self.p_mod, self.p_mod_time, t_start, t_end, 
                dae_vars=self.p_mod._NMPC_NAMESPACE.dae_vars, 
                time_linking_vars=self.p_mod._NMPC_NAMESPACE.diff_vars,
                outlvl=outlvl)
        msg = ('Successfully simulated plant over the sampling period '
                'beginning at ' + str(t_start))
        init_log.info(msg)

        tc1 = self.c_mod.time.first() + sample_time

        if self.controller_solved and calculate_error:
            self.state_error[t_end] = self.calculate_error_between_states(
                    self.c_mod, self.p_mod, tc1, t_end)


    def calculate_error_between_states(self, mod1, mod2, t1, t2, 
            state_attrname='diff_vars', 
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

        Kwargs:
            Q_diagonal : Flag for whether weighting matrix is diagonal. Default
                         True. False is not supported for now.
            Q_matrix : Weighting "matrix." For now just a list of values to 
                       weight the error between each state. Default is to use
                       the same weights calculated for controller objective
                       function.
            state_attrname : Name of the list attribute containing the
                             variables whose errors will be calculated.
            state_attrname_1 : Name of variable list in first model if
                               different
            state_attrname_2 : Name of variable list in second model if
                               different
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
        Q_matrix = []
        for categ in categories:
            varlist_1 += mod1._NMPC_NAMESPACE.category_dict[categ].varlist
            Q_matrix += mod1._NMPC_NAMESPACE.category_dict[categ].weights
            varlist_2 += mod2._NMPC_NAMESPACE.category_dict[categ].varlist
        assert len(varlist_1) == len(varlist_2)
        n = len(varlist_1)

        assert t1 in mod1._NMPC_NAMESPACE.get_time()
        assert t2 in mod2._NMPC_NAMESPACE.get_time()

        error = sum(Q_matrix[i]*(varlist_1[i][t1].value - 
                                 varlist_2[i][t2].value)**2
                    for i in range(n))

        return error

