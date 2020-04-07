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
from pyomo.kernel import ComponentSet, ComponentMap
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.dae.flatten import flatten_dae_variables

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.dyn_utils import (get_activity_dict, deactivate_model_at,
        path_from_block, find_comp_in_block, find_comp_in_block_at_time)
from idaes.core.util.initialization import initialize_by_time_element
from idaes.dynamic.caprese.util import (initialize_by_element_in_range,
        find_slices_in_model, VarLocator, copy_values_at_time, 
        add_noise_at_time)
import idaes.logger as idaeslog

from collections import OrderedDict
import time as timemodule
import pdb

__author__ = "Robert Parker and David Thierry"


# TODO tags:
# - KWARGS
# - CONFIG
# - NAMESPACE
# - RENAME
# - REMOVECONSTANT


class NMPCSim(object):
    """
    Main class for NMPC simulations of IDAES flowsheets. 
    """

    def __init__(self, plant_model, controller_model, initial_inputs, **kwargs):
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
        # Maybe include a kwarg for require_steady - if False, set-point is not
        # forced to be a steady state

        # TODO: handle arguments (solver, outlvl, ...) via config blocks
        #       separate solvers for each solve, in init
        # KWARGS, CONFIG, NAMESPACE - remove kwargs and set these up here

        solver = kwargs.pop('solver', SolverFactory('ipopt'))
        self.default_solver = solver

        # Logger properties:
        # outlvl added as an attribute so it can be referenced later
        self.outlvl = kwargs.pop('outlvl', idaeslog.NOTSET)
        init_log = idaeslog.getInitLogger('nmpc', level=self.outlvl)

        # Set up attributes
        self.p_mod = plant_model
        self.c_mod = controller_model

        self.add_namespace_to(self.p_mod)
        self.add_namespace_to(self.c_mod)

        # Validate models
        self.validate_models(self.p_mod, self.c_mod)

        # Solve for consistent initial conditions
#        self.p_mod._originally_active = get_activity_dict(self.p_mod)

        # Categorize variables in plant model
        init_log.info('Categorizing variables in plant model') 
        self.categorize_variables(self.p_mod, initial_inputs)

        self.solve_initial_conditions(self.p_mod,
                solver=solver)
        # TODO: move into own function
        #       possibly a DAE utility
        #       - check for consistency of initial conditions
        #       - if not consistent, tell user to go run dae.solve_initial_conditions

        # TODO: - add attributes as private "_initial_inputs" - make block private
        #       - check that I don't overwrite anything
        #       - store in a block: util.modeling.unique_component_name
        #       - cleanup if desired
        self.p_mod.initial_inputs = initial_inputs
        # ^ Add these as an attribute so they can be used later to validate
        # steady state model
        # ^ What did I mean by this? In what way does steady state model need
        # to resemble plant model?
        # Not sure this is still necessary

        self.build_variable_locator(self.p_mod, 
                differential=self.p_mod.diff_vars,
                derivative=self.p_mod.deriv_vars,
                algebraic=self.p_mod.alg_vars,
                inputs=self.p_mod.input_vars,
                fixed=self.p_mod.fixed_vars,
                ic=self.p_mod.ic_vars)
        # Now adding a locator to the plant model so I can find plant model
        # variables corresponding to the controller's initial conditions

        # Categorize variables in controller model
        init_controller_inputs = self.validate_initial_inputs(controller_model,
                plant_model, initial_inputs)
        init_log.info('Categorizing variables in the controller model')
        self.categorize_variables(self.c_mod, init_controller_inputs)
        self.build_variable_locator(self.c_mod, 
                differential=self.c_mod.diff_vars,
                derivative=self.c_mod.deriv_vars,
                algebraic=self.c_mod.alg_vars,
                inputs=self.c_mod.input_vars,
                fixed=self.c_mod.fixed_vars,
                ic=self.c_mod.ic_vars)
        # Only expecting user arguments (set point) in form of controller
        # variables, so only build locator for controller model
        # for now.
        # ^ Not true, see above.
        #
        # Not convinced that having a variable locator like this is the
        # best thing to do, but will go with it for now.

        # Only need to manipulate bounds of controller model. Assume the 
        # bounds in the plant model should remain in place for simulation.
        # (Should probably raise a warning if bounds are present...)
        self.build_bound_lists(self.c_mod)
        # ^ This may be removed in favor of strip_bounds transformation

        # Validate inputs in the plant model and initial conditions
        # in the control model.
        # TODO: allow user to specify this if names don't match
        self.p_mod.controller_ic_vars = find_slices_in_model(
                self.p_mod,
                self.c_mod,
                self.p_mod.var_locator,
                self.c_mod.ic_vars)
        self.c_mod.plant_input_vars = find_slices_in_model(
                self.c_mod,
                self.p_mod,
                self.c_mod.var_locator,
                self.p_mod.input_vars)

        self.validate_fixedness(self.p_mod)
        self.validate_fixedness(self.c_mod)

# TODO: remove. Place in solve_initial_conditions method if it exists.
#       If desired ('strict mode') check for consistency.
#####################
        copy_values_at_time(self.c_mod.ic_vars,
                            self.p_mod.controller_ic_vars,
                            self.c_mod.time.first(),
                            self.p_mod.time.first())
        # May eventually want to "load_initial_conditions" here so I can apply 
        # noise, but for now this is more direct.

        # Should strip bounds before this IC solve, since the controller
        # model should have bounds
        self.strip_controller_bounds = TransformationFactory(
                                       'contrib.strip_var_bounds')
        self.strip_controller_bounds.apply_to(self.c_mod, reversible=True)

        # Controller model has already been categorized... No need 
        # to provide init_controller_inputs
        self.solve_initial_conditions(self.c_mod, solver=solver)
        # ^ move into initialize_control_problem
        self.strip_controller_bounds.revert(self.c_mod)
#####################

        self.sample_time = kwargs.pop('sample_time', 
                self.c_mod.time.get_finite_elements()[1] -
                    self.c_mod.time.get_finite_elements()[0])
        self.validate_sample_time(self.sample_time, 
                self.c_mod, self.p_mod)

        scheme = self.c_mod.time.get_discretization_info()['scheme']
        if scheme == 'LAGRANGE-RADAU':
            self.c_mod._ncp = self.c_mod.time.get_discretization_info()['ncp']
        elif scheme == 'BACKWARD Difference':
            self.c_mod._ncp = 1
        else:
            raise NotImplementedError

        scheme = self.p_mod.time.get_discretization_info()['scheme']
        if scheme == 'LAGRANGE-RADAU':
            self.p_mod._ncp = self.p_mod.time.get_discretization_info()['ncp']
        elif scheme == 'BACKWARD Difference':
            self.p_mod._ncp = 1
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


    def add_namespace_to(self, model):
        name = '_NMPC_NAMESPACE'
        # Not _CAPRESE_NAMESPACE as I might want to add a similar 
        # namespace for MHE
        if hasattr(model, name):
            raise ValueError('%s already exists on model. Please fix this.'
                             % name)
        model.add_component(name, Block())
        namespace = getattr(model, name)

        namespace.input_vars = []
        namespace.n_input_vars = 0
        namespace.input_weights = []
        namespace.input_setpoints = []
        namespace.input_bounds = []

        namespace.alg_vars = []
        namespace.n_alg_vars = 0
        namespace.alg_weights = []
        namespace.alg_setpoints = []
        namespace.alg_bounds = []

        namespace.diff_vars = []
        namespace.n_diff_vars = 0
        namespace.diff_weights = []
        namespace.diff_setpoints = []
        namespace.diff_bounds = []

        namespace.deriv_vars = []
        namespace.n_deriv_vars = 0
        namespace.deriv_weights = []
        namespace.deriv_setpoints = []
        namespace.deriv_bounds = []

        namespace.fixed_vars = []
        namespace.n_fixed_vars = 0
        namespace.fixed_weights = []
        namespace.fixed_setpoints = []
        namespace.fixed_bounds = []

        namespace.scalar_vars = []
        namespace.n_scalar_vars = 0
        namespace.scalar_weights = []
        namespace.scalar_setpoints = []
        namespace.scalar_bounds = []

        namespace.ic_vars = []
        namespace.n_ic_vars = 0

        namespace.dae_vars = []
        namespace.n_dae_vars = 0

        namespace.var_locator = ComponentMap()


    def validate_sample_time(self, sample_time, *models):
        """Makes sure sample time is an integer multiple of discretization
        spacing in each model, and that horizon of each model is an integer
        multiple of sample time.

        Args:
            sample_time: Sample time to check
            models: List of flowsheet models to check
        """
        for model in models:
            time = model.time
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
                model._fe_per_sample = int(sample_time/fe_spacing)

            horizon_length = time.last() - time.first()
            if (horizon_length / sample_time) % 1 != 0:
                raise ValueError(
                    'Sampling time must be an integer divider of '
                    'horizon length')
            else:
                model._samples_per_horizon = int(horizon_length/sample_time)


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
        # the target model, then use the VarLocator to find the corresponding
        # slice in the target model...
        t0 = src_model.time.first()
        tgt_slices = []
        for _slice in src_slices:
            init_vardata = _slice[t0]
            tgt_vardata = find_comp_in_block(tgt_model, 
                                             src_model, 
                                             init_vardata)
            tgt_slice = tgt_model.var_locator[id(tgt_vardata)].container
            location = tgt_model.var_locator[id(tgt_vardata)].location
            tgt_slices.append(tgt_slice[location])
        return tgt_slices


    def validate_fixedness(self, model):
        """
        Makes sure that assumptions regarding fixedness for different points
        in time are valid. Differential, algebraic, and derivative variables
        may be fixed only at t0, only if they are initial conditions.
        Fixed variables must be fixed at all points in time, except possibly
        initial conditions.
        """
        time = model.time
        t0 = time.first()
        locator = model.var_locator

        for _slice in model.alg_vars + model.diff_vars + model.deriv_vars:
            var0 = _slice[t0]
            if locator[id(var0)].is_ic:
                assert var0.fixed
                for t in time:
                    if t == t0:
                        continue
                    assert not _slice[t].fixed
            else:
                for t in time:
                    assert not _slice[t].fixed

        for var in model.fixed_vars:
            for t in time:
                # Fixed vars, e.g. those used in boundary conditions,
                # may "overlap" with initial conditions. It is up to the user
                # to make sure model has appropriate number of degrees of
                # freedom
                if t == t0:
                    continue
                assert var[t].fixed
                    

    def validate_initial_inputs(self, tgt_model, src_model, src_inputs=None,
            outlvl=idaeslog.NOTSET):
        """Uses initial inputs in the source model to find variables of the
        same name in a target model.
        
        Args:
           tgt_model : Flowsheet model to search for input variables
           src_model : Flowsheet model containing inputs to search for
           src_inputs : List of input variables at the initial time point
                        to find in target model. If not provided, the
                        initial_inputs attribute will be used.
        """
        # TODO: set level through config
        log = idaeslog.getInitLogger('nmpc', level=outlvl)

        if not src_inputs:
            # If source inputs are not specified, assume they
            # already exist in src_model
            try:
                t0 = src_model.time.first()
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
        if not (isinstance(m1, FlowsheetBlock) and
                isinstance(m2, FlowsheetBlock)):
            raise ValueError(
                    'Provided models must be FlowsheetBlocks')
        if m1.model() is m2.model():
            raise ValueError(
                    'Provided models must not live in the same top-level'
                    'ConcreteModel')
        return True


    def transfer_current_plant_state_to_controller(self, t_plant, **kwargs):
        # Would like to pass "noise_args" in as a bundle here. This can
        # probably be done with config blocks somehow.

        # TESTME: this function is not tested.
        # copy_values and add_noise are tested, however

        time = self.c_mod.time
        t0 = time.first()

        copy_values_at_time(self.c_mod.ic_vars,
                            self.p_mod.controller_ic_vars,
                            t0,
                            t_plant)

        # Apply noise to new initial conditions
        add_noise = kwargs.pop('add_noise', False)
        noise_weights = kwargs.pop('noise_weights', [])
        noise_sig_0 = kwargs.pop('noise_sig_0', 0.05)
        noise_args = kwargs.pop('noise_args', {})

        if add_noise:
            if not noise_weights:
                noise_weights = []
                for var in self.c_mod.ic_vars:
                    locator = self.c_mod.var_locator[id(var[t0])]
                    category = locator.category
                    loc = locator.location
                    # REMOVECONSTANT
                    # There should be an enum for variable category
                    if category == 'differential':
                        obj_weight = self.c_mod.diff_weights[loc]
                    elif category == 'derivative':
                        # Is it okay to weigh derivative noise by 
                        # the state variable's weight?
                        obj_weight = self.c_mod.diff_weights[loc]
                    else:
                        # FIXME: Won't have obj weight for algebraic equations
                        # NAMESPACE: empty list should be added automatically in
                        # namespace block
                        obj_weight = None

                    if obj_weight is not None and obj_weight != 0:
                        # TODO: make an option for max noise weight
                        noise_weights.append(min(1/obj_weight, 1e6))
                    else:
                        noise_weights.append(None)

            add_noise_at_time(self.c_mod.ic_vars,
                              t0,
                              weights=noise_weights,
                              sigma_0=noise_sig_0,
                              **noise_args)


    def inject_control_inputs_into_plant(self, t_plant, **kwargs):
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
            src_attrname : Name (string) of the list attribute containing
                           the variables in the controller model corresponding
                           to the plant model's inputs. (These are not
                           necessarily the control model's inputs!)
            add_noise : Bool telling whether or not to apply noise
            noise_weights : List of weights for each state's standard deviation
            noise_sig_0 : Standard deviation for a state with unit weight
            noise_args : Additional kwargs to pass apply_noise_at_time

        """
        # KWARGS, CONFIG - handle all of these through config. Pass in noise args
        # as a bundle
        # config args for control_input_noise

        sample_time = kwargs.pop('sample_time', self.sample_time)
        src_attrname = kwargs.pop('src_attrname', 'plant_input_vars')

        # Send inputs to plant that were calculated for the end
        # of the first sample
        t_controller = self.c_mod.time.first() + sample_time
        assert t_controller in self.c_mod.time

        time = self.p_mod.time
        plant_sample = [t for t in time if t > t_plant and t<= t_plant+sample_time]
        assert t_plant+sample_time in plant_sample
        # len(plant_sample) should be ncp*nfe_per_sample, assuming the expected
        # sample_time is passed in

        add_noise = kwargs.pop('add_noise', False)
        noise_weights = kwargs.pop('noise_weights', [])
        noise_sig_0 = kwargs.pop('noise_sig_0', 0.05)
        noise_args = kwargs.pop('noise_args', {})

        # Need to get proper weights for plant's input vars
        if add_noise:
            if not noise_weights:
                noise_weights = []
                for var in self.c_mod.plant_input_vars:
                    locator = self.c_mod.var_locator[id(var[t_controller])]
                    category = locator.category
                    loc = locator.location
                    if category == 'input':
                        obj_weight = self.c_mod.input_weights[loc]
                    elif category == 'differential':
                        obj_weight = self.c_mod.diff_weights[loc]
                    if obj_weight is not None and obj_weight != 0:
                        noise_weights.append(min(1/obj_weight, 1e6))
                    else:
                        # By default, if state is not penalized in objective,
                        # noise will not be applied to it here.
                        # This may be incorrect, but user will have to override,
                        # by providing their own weights, as I don't see a good
                        # way of calculating a weight
                        noise_weights.append(None)

            add_noise_at_time(self.c_mod.plant_input_vars,
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

        copy_values_at_time(self.p_mod.input_vars,
                            self.c_mod.plant_input_vars,
                            plant_sample,
                            t_controller)


    def solve_initial_conditions(self, model, **kwargs):
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

        was_originally_active = get_activity_dict(model)
        solver = kwargs.pop('solver', self.default_solver)
        outlvl = kwargs.pop('outlvl', self.outlvl)
        init_log = idaeslog.getInitLogger('nmpc', level=outlvl)
        solver_log = idaeslog.getSolveLogger('nmpc', level=outlvl)

        toplevel = model.model()
        t0 = model.time.first()

        non_initial_time = [t for t in model.time]
        non_initial_time.remove(t0)
        deactivated = deactivate_model_at(model, model.time, non_initial_time, 
                outlvl=idaeslog.ERROR)

        for _slice in model.input_vars:
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


    def categorize_variables(self, model, initial_inputs):
        """Function to create lists of time-only-slices of the different 
        types of variables in a model, given knowledge of which are inputs. 
        These lists are added as attributes to the model.

        Args:
            model : Model whose variables will be flattened and categorized
            initial_inputs : List of VarData objects that are input variables
                             at the initial time point

        Kwargs:
            is_steady : Flag to set True if the model is a steady state model,
                        and thus its time set is a singleton.
        """

        # Would probably be much faster to do this categorization
        # during variable flattening
        time = model.time

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
        model.dae_vars = list(dae_map.values())
        model.scalar_vars = list(ComponentMap([(v, v) for v in scalar_vars]).values())
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

        # TODO: attribute of subblock
        model.deriv_vars = deriv_vars
        model.diff_vars = diff_vars
        model.n_dv = len(diff_vars)
        assert model.n_dv == len(deriv_vars)

        model.ic_vars = ic_vars
        #assert model.n_dv == len(ic_vars)
        # Would like this to be true, but accurately detecting differential
        # variables that are not implicitly fixed (by fixing some input)
        # is difficult

        model.input_vars = input_vars
        model.n_iv = len(input_vars)

        model.alg_vars = alg_vars
        model.n_av = len(alg_vars)

        model.fixed_vars = fixed_vars 
        model.n_nv = len(fixed_vars)


    def build_variable_locator(self, model, **kwargs):
        """Constructs a dictionary mapping the id of each VarData object
        to a VarLocator object. This dictionary is added as an attribute to
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
        alg_list = kwargs.pop('algebraic', [])
        diff_list = kwargs.pop('differential', [])
        deriv_list = kwargs.pop('derivative', [])
        input_list = kwargs.pop('inputs', []) # note the s. input is reserved
        fixed_list = kwargs.pop('fixed', [])
        scalar_list = kwargs.pop('scalar', [])
        ic_list = kwargs.pop('ic', [])

        locator = {}

        for i, _slice in enumerate(alg_list):
            for t in model.time:
                locator[id(_slice[t])] = VarLocator('algebraic', alg_list, i)

        for i, _slice in enumerate(diff_list):
            for t in model.time:
                locator[id(_slice[t])] = VarLocator('differential', diff_list, i)

        for i, _slice in enumerate(deriv_list):
            for t in model.time:
                locator[id(_slice[t])] = VarLocator('derivative', deriv_list, i)

        for i, _slice in enumerate(input_list):
            for t in model.time:
                locator[id(_slice[t])] = VarLocator('input', input_list, i)

        for i, _slice in enumerate(fixed_list):
            for t in model.time:
                locator[id(_slice[t])] = VarLocator('fixed', fixed_list, i)

        for i, vardata in enumerate(scalar_list):
            for t in model.time:
                locator[id(vardata)] = VarLocator('scalar', scalar_list, i)

        # Since these variables have already VarLocator objects,
        # just set the desired attribute.
        for i, _slice in enumerate(ic_list):
            for t in model.time:
                locator[id(_slice[t])].is_ic = True

        model.var_locator = locator


    def add_setpoint(self, set_point, **kwargs):
        """User-facing function for the addition of a set point to the 
        controller. Will validate th
        """
        # TODO: maybe split into two functions
        skip_validation = kwargs.pop('skip_validation', False)
        outlvl = kwargs.pop('outlvl', idaeslog.NOTSET)
        weight_overwrite = kwargs.pop('steady_weight_overwrite', [])
        weight_tolerance = kwargs.pop('steady_weight_tolerance', 1e-6)
        dynamic_weight_overwrite = kwargs.pop('dynamic_weight_overwrite', [])
        dynamic_weight_tol = kwargs.pop('dynamic_weight_tol', 1e-6)

        if skip_validation:
            raise NotImplementedError(
                    'Maybe one day...')
            # In the implementation to 'skip validation', should still
            # get an error if a set_point var is not present in controller
            # model.
            #
            # Result should be that setpoint attributes have been added to 
            # controller model
        else:
            steady_model = kwargs.pop('steady_model', None)
            if not steady_model:
                raise ValueError(
                   "'steady_model' required to validate set point")
            self.s_mod = steady_model
            self.validate_models(self.s_mod, self.p_mod)
            self.validate_steady_setpoint(set_point, self.s_mod,
                                          outlvl=outlvl,
                                          weight_overwrite=weight_overwrite,
                                          weight_tolerance=weight_tolerance)
            # ^ result should be that controller now has set point attributes

        self.construct_objective_weight_matrices(self.c_mod,
                weight_overwrite=dynamic_weight_overwrite,
                tol=dynamic_weight_tol)

        self.add_objective_function(self.c_mod,
                control_penalty_type='action',
                name='tracking_objective')


    def validate_steady_setpoint(self, set_point, steady_model, **kwargs):
        solver = kwargs.pop('solver', self.default_solver)
        outlvl = kwargs.pop('outlvl', self.outlvl)
        init_log = idaeslog.getInitLogger('nmpc', level=outlvl)
        solver_log = idaeslog.getSolveLogger('nmpc', level=outlvl)
        # Above can be setup through common config
        weight_overwrite = kwargs.pop('weight_overwrite', [])
        weight_tolerance = kwargs.pop('weight_tolerance', 1e-6)

        # The following loop will create steady-state variable lists in
        # proper order, initialize steady state model to initial conditions
        # of controller model, and make sure that the fixed status of 
        # variables is the same across these models.

        # Assume that variable with time-index t0 exists in steady model
        t0 = self.c_mod.time.first()
        # Compared fixed-ness as t1 so we're not thrown off by fixed
        # initial conditions
        t1 = self.c_mod.time.get_finite_elements()[1]
        for attr_name in ['diff_vars', 'alg_vars', 'input_vars', 
                          'fixed_vars', 'scalar_vars']:
            # TODO: subblock for naming
            setattr(steady_model, attr_name, [])
            for _slice in getattr(self.c_mod, attr_name):
                if attr_name != 'scalar_vars':
                    vardata_t0 = _slice[t0]
                    vardata_t1 = _slice[t1]
                else:
                    vardata_t0 = _slice
                    vardata_t1 = _slice

                local_parent = steady_model
                # TODO: replace with helper function
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
                
                varlist = getattr(steady_model, attr_name)
                if attr_name != 'scalar_vars':
                    varlist.append({t0: var_steady})
                    # Append dict to list here so that steady state
                    # VarDatas for "DAE vars" are accessed with the
                    # same syntax as in the dynamic case.
                else:
                    varlist.append(var_steady)

                # Copy fixed status from controller model at t1
                if not var_steady.fixed == vardata_t1.fixed:
                    var_steady.fixed = vardata_t1.fixed

        assert len(steady_model.diff_vars) == len(self.c_mod.diff_vars)
        assert len(steady_model.alg_vars) == len(self.c_mod.alg_vars)
        assert len(steady_model.input_vars) == len(self.c_mod.input_vars)
        assert len(steady_model.fixed_vars) == len(self.c_mod.fixed_vars)
        assert len(steady_model.scalar_vars) == len(self.c_mod.scalar_vars)

        # Now that steady state model has been validated and initialized,
        # construct complete lists for the *_sp attributes
        diff_var_sp = [None for var in steady_model.diff_vars]
        alg_var_sp = [None for var in steady_model.alg_vars]
        input_var_sp = [None for var in steady_model.input_vars]
        # ^ list entries are not actually vars, they are dictionaries...
        # maybe this is too confusing

        # This is where I map user values for set points into lists that
        # I can use to build the objective function (and weight matrices).
        for vardata, value in set_point:
            # set_point variables should be members of controller model
            info = self.c_mod.var_locator[id(vardata)]
            category = info.category
            location = info.location

            # Only allow set point variables from following categories:
            if category == 'differential':
                diff_var_sp[location] = value
            elif category == 'algebraic':
                alg_var_sp[location] = value
            elif category == 'input':
                input_var_sp[location] = value
            else:
                raise ValueError('Set point variables (data objects) must be'
                                 'differential, algebraic, or inputs.')

        # Add attributes (with proper names!) to steady model
        # so these values can be accessed later
        # TODO: subblock for naming
        steady_model.diff_sp = diff_var_sp
        steady_model.alg_sp = alg_var_sp
        steady_model.input_sp = input_var_sp

        # TODO: ComponentMap
        self.build_variable_locator(steady_model,
                differential=steady_model.diff_vars,
                algebraic=steady_model.alg_vars,
                inputs=steady_model.input_vars,
                fixed=steady_model.fixed_vars,
                scalar=steady_model.scalar_vars)
        self.construct_objective_weight_matrices(steady_model,
                overwrite=weight_overwrite,
                tol=weight_tolerance)

        assert hasattr(steady_model, 'diff_weights')
        assert hasattr(steady_model, 'alg_weights')
        assert hasattr(steady_model, 'input_weights')

        assert len(steady_model.diff_weights) == len(steady_model.diff_vars)
        assert len(steady_model.alg_weights) == len(steady_model.alg_vars)
        assert len(steady_model.input_weights) == len(steady_model.input_vars)

        # Add objective to steady_model
        # Add bounds to steady_model (taken from control model)
        #
        # Make sure inputs are unfixed - validate degrees of freedom
        #
        # solve steady model for set point
        # add sp attributes to control model

        self.add_objective_function(steady_model,
                                    control_penalty_type='error',
                                    name='user_setpoint_objective')
        self.transfer_bounds(steady_model, self.c_mod)

        for var in steady_model.input_vars:
            for t in steady_model.time:
                var[t].unfix()

        assert (degrees_of_freedom(steady_model) ==
                len(steady_model.input_vars)*len(steady_model.time))

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

        c_diff_sp = []
        for var in steady_model.diff_vars:
            c_diff_sp.append(var[t0].value)
        self.c_mod.diff_sp = c_diff_sp

        c_input_sp = []
        for var in steady_model.input_vars:
            c_input_sp.append(var[t0].value)
        self.c_mod.input_sp = c_input_sp

        # TODO: add alg var sp values


    def construct_objective_weight_matrices(self, model, **kwargs):
        """
        Do I even need the model? 
        ^ it makes things slightly easier to provide only the model, I think...

        Do I want to allow user to provide a setpoint here?
        ^ No, but want to allow them to overwrite weights if they want

        Here each setpoint is a list of (VarData, value) tuples.

        Overwrite is a lists of (VarData, value) tuples, where 
        this value will directly overwrite the weight of the
        corresponding VarData. (For all time..., what if user
        provides multiple weights for VarDatas that only differ
        by time? Raise warning and use last weight provided)
        """
        overwrite = kwargs.pop('overwrite', [])

        # Variables to overwrite must be VarData objects in the model
        # for whose objective function we are calculating weights
        
        weights_to_overwrite = {}
        for ow_tpl in overwrite:
            locator = model.var_locator[id(ow_tpl[0])]
            weights_to_overwrite[(locator.category, locator.location)] = \
                    ow_tpl[1]
        # dictionary mapping id of VarDatas to overwrite to specified weight
        # Is the proper way to do this to use suffixes?
        # But which vardata (time index) is the right one to have the suffix?

        # Given a vardata here, need to know its location so I know which
        # weight to overwrite

        tol = kwargs.pop('tol', 1e-6)

        # Need t0 to get the 'reference' value to compare with set point value
        # for weight calculation
        t0 = model.time.first()

        # Attempt to construct weight for each type of setpoint
        # that we could allow

        # TODO: subblock for naming
        for name in ['diff', 'alg', 'input']:
            sp_attrname = name + '_sp'
            if not hasattr(model, sp_attrname):
                continue
            sp = getattr(model, sp_attrname)

            varlist_attrname = name + '_vars'
            varlist = getattr(model, varlist_attrname)

            # construct the (diagonal) matrix (list).
            matrix = []
            for loc, value in enumerate(sp):

                # This assumes the vardata in sp is the same one
                # provided by the user. But these could differ by time
                # index...
                # Need to check by location, category here
                if name == 'diff':
                    if ('differential', loc) in weights_to_overwrite:
                        weight = weights_to_overwrite['differential', loc]
                        matrix.append(weight)
                        continue
                elif name == 'alg':
                    if ('algebraic', loc) in weights_to_overwrite:
                        weight = weights_to_overwrite['algebraic', loc]
                        matrix.append(weight)
                        continue
                elif name == 'input':
                    if ('input', loc) in weights_to_overwrite:
                        weight = weights_to_overwrite['input', loc]
                        matrix.append(weight)
                        continue

                # If value is None, but variable was provided as overwrite,
                # weight can still be non-None. This is okay.
                if value is None:
                    weight = None
                    matrix.append(weight)
                    continue

                diff = abs(varlist[loc][t0].value - value)
                if diff > tol:
                    weight = 1./diff
                else:
                    weight = 1./tol
                matrix.append(weight)

            # TODO: subblock for naming
            weight_attrname = name + '_weights'
            setattr(model, weight_attrname, matrix)


    def add_objective_function(self, model, **kwargs):
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
        # TODO: kwargs ???
        name = kwargs.pop('name', 'objective')
        outlvl = kwargs.pop('outlvl', idaeslog.NOTSET)
        init_log = idaeslog.getInitLogger('nmpc', level=outlvl)

        # Q and R are p.s.d. matrices that weigh the state and
        # control norms in the objective function
        Q_diagonal = kwargs.pop('Q_diagonal', True)
        R_diagonal = kwargs.pop('R_diagonal', True)

        state_weight = kwargs.pop('state_weight', 1)
        control_weight = kwargs.pop('control_weight', 1)
        # The names of these arguments could get confusing...

        # User may want to penalize control action, i.e. ||u_i - u_{i-1}||,
        # or control error (from set point), i.e. ||u_i - u*||
        # Valid values are 'action' or 'error'
        control_penalty_type = kwargs.pop('control_penalty_type', 'action')
        if not (control_penalty_type == 'action' or
                control_penalty_type == 'error'):
            raise ValueError(
                "control_penalty_type argument must be 'action' or 'error'")

        # TODO: rename Q/R more descriptive, or comment well.
        if not Q_diagonal or not R_diagonal:
            raise NotImplementedError('Q and R must be diagonal for now.')

        # This implementation assumes diff_weights/diff_sp are just values
        # TODO: square this with implementation above
        if hasattr(model, 'alg_sp') and hasattr(model, 'alg_weights'):
            # Lists of time-indexed 'variables' - this is fine
            states = model.diff_vars + model.alg_vars
            Q_entries = model.diff_weights + model.alg_weights
            sp_states = model.diff_sp + model.alg_sp
        else:
            states = model.diff_vars
            Q_entries = model.diff_weights
            sp_states = model.diff_sp
        # TODO: Setup all data structures (namespace attributes) in constructor
        # and document there to avoid lots of if statements like this.

        controls = model.input_vars
        R_entries = model.input_weights
        sp_controls = model.input_sp

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
        # Not necessarily true, Q can have an entry due to an overwrite

        if control_penalty_type == 'error':
            control_term = sum(R_entries[i]*(controls[i][t] - sp_controls[i])**2
                    for i in range(len(controls)) if (R_entries[i] is not None
                                                and sp_controls[i] is not None)
                               for t in time)
        elif control_penalty_type == 'action':
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
            # Note: This term is only non-zero at the boundary between sampling
            # times. Could use this info to make the expression more compact
            # (TODO)
            # Populate time as either the list of finite element of the list
            # of sample points. Change time to be a list, not ContinuousSet

        obj_expr = state_term + control_term

        # TODO: namespace block
        obj = Objective(expr=obj_expr)
        model.add_component(name, obj)


    def build_bound_lists(self, model):
        """
        Builds lists of lower bound, upper bound tuples as attributes of the 
        input model, based on the current bounds (and domains) of
        differential, algebraic, and input variables.

        Args:
            model : Model whose variables will be checked for bounds.

        Returns:
            None
        """
        time = model.time
        t0 = time.first()

        model.diff_var_bounds = []
        for var in model.diff_vars:
            # Assume desired bounds are given at t0
            lb = var[t0].lb
            ub = var[t0].ub
            if (var[t0].domain == NonNegativeReals and lb is None):
                lb = 0
            elif (var[t0].domain == NonNegativeReals and lb < 0):
                lb = 0
            for t in time:
                var[t].setlb(lb)
                var[t].setub(ub)
                var[t].domain = Reals
            model.diff_var_bounds.append((lb, ub))

        model.alg_var_bounds = []
        for var in model.alg_vars:
            # Assume desired bounds are given at t0
            lb = var[t0].lb
            ub = var[t0].ub
            if (var[t0].domain == NonNegativeReals and lb is None):
                lb = 0
            elif (var[t0].domain == NonNegativeReals and lb < 0):
                lb = 0
            for t in time:
                var[t].setlb(lb)
                var[t].setub(ub)
                var[t].domain = Reals
            model.alg_var_bounds.append((lb, ub))

        model.input_var_bounds = []
        for var in model.input_vars:
            # Assume desired bounds are given at t0
            lb = var[t0].lb
            ub = var[t0].ub
            if (var[t0].domain == NonNegativeReals and lb is None):
                lb = 0
            elif (var[t0].domain == NonNegativeReals and lb < 0):
                lb = 0
            for t in time:
                var[t].setlb(lb)
                var[t].setub(ub)
                var[t].domain = Reals
            model.input_var_bounds.append((lb, ub))


    def transfer_bounds(self, tgt_model, src_model):
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
        time = tgt_model.time

        diff_var_bounds = src_model.diff_var_bounds
        for i, var in enumerate(tgt_model.diff_vars):
            for t in time:
                var[t].setlb(diff_var_bounds[i][0])
                var[t].setub(diff_var_bounds[i][1])
                var[t].domain = Reals

        alg_var_bounds = src_model.alg_var_bounds
        for i, var in enumerate(tgt_model.alg_vars):
            for t in time:
                var[t].setlb(alg_var_bounds[i][0])
                var[t].setub(alg_var_bounds[i][1])
                var[t].domain = Reals

        input_var_bounds = src_model.input_var_bounds
        for i, var in enumerate(tgt_model.input_vars):
            for t in time:
                var[t].setlb(input_var_bounds[i][0])
                var[t].setub(input_var_bounds[i][1])
                var[t].domain = Reals


# RENAME
#    def add_pwc_constraints(self, **kwargs):
    def constrain_control_inputs_piecewise_constant(self, **kwargs):
        """Function to add piecewise constant (PWC) constraints to
        model. Inputs and sample time are already known, so no arguments are
        necessary.

        Kwargs:
            model : Model to which PWC constraints are added. Default is   
                    controller model
            sample_time : Duration for which inputs will be forced constant
            outlvl : idaes.logger output level
        """
        # KWARGS
        # CONFIG: all of these options can/should be handled through config
        model = kwargs.pop('model', self.c_mod)
        sample_time = kwargs.pop('sample_time', self.sample_time)
        outlvl = kwargs.pop('outlvl', self.outlvl)
        init_log = idaeslog.getInitLogger('nmpc', outlvl)
        init_log.info('Adding piecewise-constant constraints')

        # If sample_time is overwritten here, assume that the 
        # provided sample_time should be used going forward
        # (in input injection, plant simulation, and controller initialization)
        if sample_time != self.sample_time:
            self.validate_sample_time(sample_time, self.c_mod, self.p_mod)
            self.sample_time = sample_time

        time = model.time
        t0 = model.time.get_finite_elements()[0]
        t1 = model.time.get_finite_elements()[1]
        nfe_spacing = t1-t0
        nfe_per_sample = sample_time / nfe_spacing 
        if nfe_per_sample - int(nfe_per_sample) != 0:
            raise ValueError
        nfe_per_sample = int(nfe_per_sample)

        pwc_constraints = []

        # This rule will not be picklable as it is not declared
        # at module namespace
        # Can access sample_time as attribute of namespace block,
        # then rule can be located outside of class
        input_indices = [i for i in range(len(model.input_vars))]
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
            # NAMESPACE: input_vars. Also, sample_time
            # probably make sample time an attribute of NMPCSim, but still
            # provide as a config. CONFIG
            inputs = m.input_vars
            _slice = inputs[i]
            return _slice[t_next] == _slice[t]
        # NAMESPACE: pwc_input
        name = '_pwc_input'
        pwc_constraint = Constraint(model.time, input_indices, 
                rule=pwc_rule)
        # Maybe a check here to delete existing pwc constraints
        model.add_component(name, pwc_constraint)

        pwc_constraints = [Reference(pwc_constraint[:, i])
                           for i in input_indices]
        # TODO: namespace block
        model._pwc_constraints = pwc_constraints


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
        strategy = kwargs.pop('strategy', 'from_previous')
        solver = kwargs.pop('solver', self.default_solver)
        solver_options = kwargs.pop('solver_options', None)
        if solver_options:
            solver.options = solver_options

        input_type = kwargs.pop('input_type', 'set_point')

        time = self.c_mod.time

        # TODO: remove embedded strings. 
        # Handle these options in config block? Here, and in constructor - 'ephemeral'
        # ^ domain maps string to enum, or just takes enum
        # Here, options should be enums
        # Read documentation in pyutilib.misc.config
        if strategy == 'from_previous':
            self.initialize_from_previous_sample(self.c_mod)

        elif strategy == 'from_simulation':
            # TODO: 'by_time_element'
            self.initialize_by_solving_elements(self.c_mod)

        elif strategy == 'initial_conditions':
            self.initialize_from_initial_conditions(self.c_mod)
        
        # Add check that initialization did not violate bounds/equalities?

        self.controller_solved = False


    def initialize_by_solving_elements(self, model, input_type='set_point'):
        time = model.time 
        # Strip bounds before simulation as square solves will be performed
        strip_controller_bounds = TransformationFactory(
                                      'contrib.strip_var_bounds')
        strip_controller_bounds.apply_to(model, reversible=True)

        if input_type == 'set_point':
            for i, _slice in enumerate(model.input_vars):
                for t in time:
                    if t != time.first():
                        _slice[t].fix(model.input_sp[i])
                    else:
                        _slice[t].fix()

        elif input_type == 'initial':
            for i, _slice in enumerate(model.input_vars):
                t0 = time.first()
                for t in time:
                    _slice[t].fix(model.input_vars[i][t0].value)

        else:
            raise ValueError('Unrecognized input type')
        # The above should ensure that all inputs are fixed and the 
        # model has no dof upon simulation

        # Deactivate objective function
        # Here I assume the name of the objective function.
        model.tracking_objective.deactivate()
        # Deactivate pwc constraints (as inputs are fixed)
        #for con in model._pwc_constraints:
        #    con.deactivate()
        model._pwc_input.deactivate()

        initialize_by_element_in_range(self.c_mod, time.first(), time.last(),
                            dae_vars=self.c_mod.dae_vars,
                            time_linking_variables=self.c_mod.diff_vars)
        # TODO: add function that is wrapper around initialize_by_element_in_range
        #       then map option to function

        # Reactivate objective, pwc constraints, bounds
        self.c_mod.tracking_objective.activate()
        #for con in self.c_mod._pwc_constraints:
        #    con.activate()
        model._pwc_input.activate()

        for _slice in self.c_mod.input_vars:
            for t in time:
                _slice[t].unfix()

        strip_controller_bounds.revert(self.c_mod)


    def initialize_from_previous_sample(self, model, **kwargs):
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

        sample_time = kwargs.pop('sample_time', self.sample_time)

        attr_list = kwargs.pop('attr_list', 
                ['diff_vars', 'alg_vars', 'deriv_vars', 'input_vars'])
        # Replace this argument with a single list of flattened variables.
        # Use as default a non_fixed_variable_list, which was defined in constructor.
        # ^ or pull off namespace block if arg not provided, knowing what these are
        # called

        # ^ why not just iterate over dae_vars here?
        # If fixed vars need to be changed (time-varying known disturbance),
        # they should be changed manually. (This is my current opinion)

        # TODO
        # Should initialize dual variables here too.

        time = model.time
        steady_t = self.s_mod.time.first()

        for attrname in attr_list:
            varlist = getattr(model, attrname)
            if not attrname == 'deriv_vars':
                steady_varlist = getattr(self.s_mod, attrname)
            for i, _slice in enumerate(varlist):
                for t in time:
                    # If not in last sample:
                    if (time.last() - t) >= sample_time:
                        t_next = t + sample_time

                        # Performing addition on CtsSet indices can result in
                        # rounding errors. Round to 8th decimal place here:
                        t_next = int(round(t_next*1e8))/1e8

                        assert t_next in time
                        _slice[t].set_value(_slice[t_next].value)

                    # Otherwise initialize to final steady state
                    # (Should provide some other options here,
                    # like copy inputs from set point then simulate)
                    else:
                        if not attrname == 'deriv_vars':
                            _slice[t].set_value(
                                    steady_varlist[i][steady_t].value)
                        else:
                            _slice[t].set_value(0)


    def initialize_from_initial_conditions(self, model, **kwargs):
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
        attr_list = kwargs.pop('attr_list', ['deriv_vars', 'diff_vars', 'alg_vars'])
        # TODO: user provides their own list of variables to initialize, or
        # I grab all these attributes from the model and add them up here.
        time = model.time
        for attrname in attr_list:
            varlist = getattr(model, attrname)
            for v in varlist:
                for t in time:
                    v[t].set_value(v[0].value)

    
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
        solver = kwargs.pop('solver', self.default_solver)
        outlvl = kwargs.pop('outlvl', self.outlvl) 
        init_log = idaeslog.getInitLogger('nmpc', level=outlvl)
        s_log = idaeslog.getSolveLogger('nmpc', level=outlvl)

        time = self.c_mod.time
        for _slice in self.c_mod.input_vars:
            for t in time:
                _slice[t].unfix()
        # ^ Maybe should fix inputs at time.first()?
        # Also, this should be redundant as inputs have been unfixed
        # after initialization

        assert (degrees_of_freedom(self.c_mod) == 
                self.c_mod.n_iv*(self.c_mod._samples_per_horizon + 1))

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
        sample_time = kwargs.pop('sample_time', self.sample_time)
        calculate_error = kwargs.pop('calculate_errors', True)
        outlvl = kwargs.pop('outlvl', self.outlvl)
        init_log = idaeslog.getInitLogger('nmpc', level=outlvl)

        t_end = t_start + sample_time 
        assert t_start in self.p_mod.time
        # TODO: add tolerance here, as t_end could change due to roundoff
        #       Then change t_end s.t. it is a point in time
        # A helper function may be useful - is_in_to_tolerance
        assert t_end in self.p_mod.time

        initialize_by_element_in_range(self.p_mod, t_start, t_end, 
                dae_vars=self.p_mod.dae_vars, 
                time_linking_vars=self.p_mod.diff_vars,
                outlvl=outlvl)
        msg = ('Successfully simulated plant over the sampling period '
                'beginning at ' + str(t_start))
        init_log.info(msg)

        tc1 = self.c_mod.time.first() + sample_time

        if self.controller_solved and calculate_error:
            self.state_error[t_end] = self.calculate_error_between_states(
                    self.c_mod, self.p_mod, tc1, t_end)


    def calculate_error_between_states(self, mod1, mod2, t1, t2, **kwargs):
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

        Q_diagonal = kwargs.pop('Q_diagonal', True)
        if not Q_diagonal:
            raise ValueError('Only diagonal weighting matrices are supported')
        Q_matrix = kwargs.pop('Q_matrix', self.c_mod.diff_weights)
        # Grab the weighting matrix from the controller model regardless of what
        # mod1 and mod2 are. This can be overwritten if desired.

        # Used to specify variables other than differential to use for
        # error calculation
        state_attrname = kwargs.pop('state_attrname', 'diff_vars')
        state_attrname_1 = kwargs.pop('state_attrname_1', state_attrname)
        state_attrname_2 = kwargs.pop('state_attrname_2', state_attrname)
        # ^ To be used if 

        if (not hasattr(mod1, state_attrname_1) or 
                not hasattr(mod2, state_attrname_1)):
            raise ValueError(
                    'Cannot calculate error in model with no '
                    + state_attrname + ' attribute')

        varlist_1 = getattr(mod1, state_attrname_1)
        varlist_2 = getattr(mod2, state_attrname_2)
        assert len(varlist_1) == len(varlist_2)
        n = len(varlist_1)

        assert t1 in mod1.time
        assert t2 in mod2.time

        error = sum(Q_matrix[i]*(varlist_1[i][t1].value - 
                                 varlist_2[i][t2].value)**2
                    for i in range(n))

        return error

