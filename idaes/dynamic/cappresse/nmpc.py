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
        TransformationFactory)
from pyomo.kernel import ComponentSet
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.dae.flatten import flatten_dae_variables

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.dyn_utils import (get_activity_dict, deactivate_model_at,
        path_from_block)
from idaes.core.util.initialization import initialize_by_time_element
import idaes.logger as idaeslog

from collections import OrderedDict
import time as timemodule
import pdb

__author__ = "Robert Parker and David Thierry"


# TODO: move this to dyn_utils module
def find_comp_in_block(tgt_block, src_block, src_comp, **kwargs):
    """
    Returns:
        Component with the same name in the target block
    """
    outlvl = kwargs.pop('outlvl', idaeslog.NOTSET)
    init_log = idaeslog.getInitLogger(__name__, outlvl)

    allow_miss = kwargs.pop('allow_miss', False)

    local_parent = tgt_block
    for r in path_from_block(src_comp, src_block, include_comp=True):
        # Better name for include_comp might be include_leaf
        try:
            local_parent = getattr(local_parent, r[0])[r[1]]
        except AttributeError:
            init_log.warning('Warning.')
            if not allow_miss:
                raise
        except KeyError:
            init_log.warning('Warning.')
            if not allow_miss:
                raise
    tgt_comp = local_parent
    return tgt_comp


class VarLocator(object):
    """
    Class for storing information used to locate a VarData object.
    Used because I want to allow the user to supply set-point in terms
    of any variables they want. I then need to find these variables in the
    proper container.
    """

    def __init__(self, category, container, location, is_ic=False):
        # Should this class store the time index of the variable?
        # probably not (a. might not exist, b. should already be known)
        if type(category) is not str:
            raise TypeError(
            'category argument must be a string')
        self.category = category

        if type(container) is not list:
            raise TypeError(
            'varlist argument must be a list')
        self.container = container

        if type(location) is not int:
            raise TypeError(
            'location argument must be an integer index')
        if location >= len(container):
            raise ValueError(
            'location must be a valid index for the container') 
        self.location = location

        if type(is_ic) is not bool:
            raise ValueError()
        self.is_ic = is_ic


class NMPCSim(object):
    """
    Main class for NMPC simulations of IDAES flowsheets. 
    """

    def __init__(self, plant_model, controller_model, initial_inputs, **kwargs):
        # Maybe include a kwarg for require_steady - if False, set-point is not
        # forced to be a steady state

        solver = kwargs.pop('solver', SolverFactory('ipopt'))
        self.default_solver = solver

        # Logger properties:
        # outlvl added as an attribute so it can be referenced later
        self.outlvl = kwargs.pop('outlvl', idaeslog.NOTSET)
        init_log = idaeslog.getInitLogger('nmpc', level=self.outlvl)

        # Set up attributes
        self.p_mod = plant_model
        self.c_mod = controller_model

        # Validate models
        self.validate_models(self.p_mod, self.c_mod)

        # Solve for consistent initial conditions
#        self.p_mod._originally_active = get_activity_dict(self.p_mod)

        # Categorize variables in plant model
        init_log.info('Categorizing variables in plant model') 
        self.categorize_variables(self.p_mod, initial_inputs)

        self.solve_initial_conditions(self.p_mod,
                solver=solver)

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
        t1 = timemodule.time() 
        self.categorize_variables(self.c_mod, init_controller_inputs)
        t2 = timemodule.time()
        print(f'categorizing controller variables took {t2-t1} seconds')
        self.build_variable_locator(self.c_mod, 
                differential=self.c_mod.diff_vars,
                derivative=self.c_mod.deriv_vars,
                algebraic=self.c_mod.alg_vars,
                inputs=self.c_mod.input_vars,
                fixed=self.c_mod.fixed_vars,
                ic=self.c_mod.ic_vars)
        t3 = timemodule.time()
        print(f'building locator took {t3-t2} seconds')
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

        # Validate inputs in the plant model and initial conditions
        # in the control model.
        self.p_mod.controller_ic_vars = self.validate_slices(
                self.p_mod,
                self.c_mod,
                self.c_mod.ic_vars)
        self.c_mod.plant_input_vars = self.validate_slices(
                self.c_mod,
                self.p_mod,
                self.p_mod.input_vars)

        self.validate_fixedness(self.p_mod)
        self.validate_fixedness(self.c_mod)

        self.load_initial_conditions_into(self.c_mod, self.p_mod,
                src_attrname='controller_ic_vars')
#        self.c_mod._originally_active = get_activity_dict(self.c_mod)

        # Should strip bounds before this IC solve, since the controller
        # model should have bounds
        self.strip_controller_bounds = TransformationFactory(
                                       'contrib.strip_var_bounds')
        self.strip_controller_bounds.apply_to(self.c_mod, reversible=True)

        # Controller model has already been categorized... No need 
        # to provide init_controller_inputs
        self.solve_initial_conditions(self.c_mod, solver=solver)
        #strip_bounds.revert(self.c_mod)

        self.sample_time = kwargs.pop('sample_time', 
                self.c_mod.time.get_finite_elements()[1] -
                    self.c_mod.time.get_finite_elements()[0])
        self.validate_sample_time(self.sample_time, 
                self.c_mod, self.p_mod)

        scheme = self.c_mod.time.get_discretization_info()['scheme']
        if scheme == 'LAGRANGE-RADAU':
            self.c_mod._ncp = time.get_discretization_info()['ncp']
        elif scheme == 'BACKWARD Difference':
            self.c_mod._ncp = 1
        else:
            raise NotImplementedError

        # Flag for whether controller has been initialized
        # by a previous solve
        self.controller_solved = False


    def validate_sample_time(self, sample_time, *models):
        for model in models:
            time = model.time
            fe_spacing = (time.get_finite_elements()[1] -
                          time.get_finite_elements()[0])
            if fe_spacing > sample_time:
                raise ValueError(
                    'Sampling time must be at least as long '
                    'as a finite element in time')
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
        Fixed variables must be fixed at all points in time.
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
        # Should probably merge with validate_vars, but the functionality is 
        # slightly different. Here...
        # Actually seems to do exactly what I want...
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


    def load_initial_conditions_into(self, tgt_model, src_model,
            src_attrname='ic_vars', **kwargs):
        t_src = kwargs.pop('t_src', src_model.time.first())
        t0 = tgt_model.time.first()

        # src_attrname provided so we can use controller_ic_vars
        src_ic_vars = getattr(src_model, src_attrname)
        for i, var_tgt in enumerate(tgt_model.ic_vars):
            var_tgt[t0].fix(src_ic_vars[i][t_src].value)


    def inject_inputs_into(self, tgt_model, src_model, **kwargs):
        """
        Sets the values of the input variables in the target model 
        over a sampling range to the values of those in the source 
        model at the specified point in time.

        Args:
            tgt_model : The target model
            src_model : The source model
            src_attrname : The attribute in the source model that contains
                           variables corresponding to the plant's input 
                           variables. Default is input_vars, but may want
                           to use plant_input_vars so as not to assume that
                           input variables are the same between models.
            t_src : The point in time from which input values in the source
                    model will be extracted.
            t_tgt : The start of the range into which input values will be 
                    copied.
            sample_time : The length of the range into which input values will
                          be copied.
        """
        src_attrname = kwargs.pop('src_attrname', 'input_vars')
        t1 = src_model.time.get_finite_elements()[1]
        t0 = src_model.time.get_finite_elements()[0]
        t_src = kwargs.pop('t_src', t1)
        t_tgt = kwargs.pop('t_tgt', tgt_model.time.first())
        sample_time = kwargs.pop('sample_time', self.sample_time)

        src_inputs = getattr(src_model, src_attrname)
        tgt_inputs = tgt_model.input_vars
        for t in tgt_model.time:
            if t <= t_tgt or t > t_tgt + sample_time:
                continue
            for i, inp in enumerate(tgt_inputs):
                inp[t].set_value(src_inputs[i][t_src].value)
                # Set values instead of fixing here - inputs should
                # be fixed already, or at least fixed before solve.


    def inject_inputs_into_plant(self, t_plant, **kwargs):

        sample_time = kwargs.pop('sample_time', self.sample_time)
        src_attrname = kwargs.pop('src_attrname', 'plant_input_vars')

        # Send inputs to plant that were calculated for the end
        # of the first sample
        t_controller = self.c_mod.time.first() + sample_time
        assert t_controller in self.c_mod.time

        self.inject_inputs_into(self.p_mod,
                                self.c_mod,
                                t_tgt=t_plant,
                                t_src=t_controller,
                                sample_time=sample_time,
                                src_attrname=src_attrname)


    def solve_initial_conditions(self, model, **kwargs):
        # Record which Constraints/Variables are initially inactive
        # deactivate model except at t=0
        # fix initial inputs - raise error if no value
        # raise error if not square
        # re-activate model
        #
        # Later include option to skip solve for consistent initial conditions
        #
        # Will only work as written for "True" initial conditions since 
        # it doesn't try to deactivate discretization equations.

        was_originally_active = get_activity_dict(model)
        solver = kwargs.pop('solver', SolverFactory('ipopt'))
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


    def categorize_variables(self, model, initial_inputs, **kwargs):
        # Would probably be much faster to do this categorization
        # during variable flattening
        is_steady = kwargs.pop('is_steady', False)
        time = model.time

        # Works for steady state models as time will be an ordered
        # (although not continuous) set:
        t0 = time.first()
        if is_steady:
            t1 = t0
        else:
            t1 = time.get_finite_elements()[1]

        deriv_vars = []
        diff_vars = []
        input_vars = []
        alg_vars = []
        fixed_vars = []

        ic_vars = []

        scalar_vars, dae_vars = flatten_dae_variables(model, time)

        # Remove duplicates from these var lists
        dae_no_dups = list(OrderedDict([(id(v[t0]), v)
                           for v in dae_vars]).values())
        scalar_no_dups = list(OrderedDict([(id(v), v)
                              for v in scalar_vars]).values())

        model.scalar_vars = scalar_no_dups

        # Remove duplicates from variable lists
        model.dae_vars = dae_no_dups.copy()

        # Find input variables
        temp_inps = initial_inputs.copy()
        for i, var in reversed(list(enumerate(dae_no_dups))):
            for j, ini in enumerate(temp_inps):
                if var[t0] is ini:
                    assert ini is temp_inps.pop(j)
                    assert var is dae_no_dups.pop(i)
                    input_vars.append(var)
                    break

        inputs_remaining = [v.name for v in temp_inps]
        if inputs_remaining:
            raise ValueError(
                f'Could not find variables for initial inputs '
                '{inputs_remaining}')

        # Find derivative variables
        init_dv_list = []
        for i, var in reversed(list(enumerate(dae_no_dups))):
            parent = var[t0].parent_component()
            index0 = var[t0].index()
            index1 = var[t1].index()
            if not isinstance(parent, DerivativeVar):
                continue
            if time not in ComponentSet(parent.get_continuousset_list()):
                continue
            assert var is dae_no_dups.pop(i)

            # Do The Correct Thing depending on which variables aren't fixed
            # Should have some robust tests for this behavior
            if parent.get_state_var()[index1].fixed:
                # Assume state var has been fixed everywhere.
                # Then derivative is 'not really a derivative'
                # in the sense that it doesn't specify a differential
                # variable through some discretization equations.
                # In this case do nothing.
                continue

            if var[t1].fixed:
                # Assume derivative has been fixed everywhere.
                # It will be added to the list of fixed variables,
                # and don't want to search for its state variable later.
                fixed_vars.append(var)
            elif var[t0].fixed:
                # In this case the derivative has been used as an
                # initial condition. Still want to include it in the
                # list of derivatives, and want to look for its state var.
                ic_vars.append(var)
                deriv_vars.append(var)
                init_dv_list.append(parent.get_state_var()[index0])
            else:
                # Neither is fixed. This should be the most common case.
                deriv_vars.append(var)
                init_dv_list.append(parent.get_state_var()[index0])

        # Find differential variables
        for i, diffvar in reversed(list(enumerate(dae_no_dups))):
            for j, dv in enumerate(init_dv_list):
                # Don't need to reverse here as we intend to break the loop
                # after the first pop
                if diffvar[t0] is dv:
                    if diffvar[t0].fixed and not diffvar[t1].fixed:
                        ic_vars.append(diffvar)
                    assert dv is init_dv_list.pop(j)
                    assert diffvar is dae_no_dups.pop(i)
                    diff_vars.append(diffvar)
                    break

        dvs_remaining = [v.name for v in init_dv_list]
        if dvs_remaining:
            raise ValueError(
                f'Could not find variables for initial states '
                '{dvs_remaining}')

        # Find algebraic variables
        for i, var in reversed(list(enumerate(dae_no_dups))):
            # If the variable is still in the list of time-indexed vars,
            # it must either be fixed (not a var) or be an algebraic var
            # - - - 
            # Check at t1 instead of t0 as algebraic vars might be fixed
            # at t0 as initial conditions instead of some differential vars
            if var[t1].fixed:
                fixed_vars.append(dae_no_dups.pop(i))
            else:
                if var[t0].fixed:
                    ic_vars.append(var)
                alg_vars.append(dae_no_dups.pop(i))


        model.deriv_vars = deriv_vars
        model.diff_vars = diff_vars
        model.n_dv = len(diff_vars)
        assert model.n_dv == len(deriv_vars)

        model.ic_vars = ic_vars
        assert model.n_dv == len(ic_vars)

        model.input_vars = input_vars
        model.n_iv = len(input_vars)

        model.alg_vars = alg_vars
        model.n_av = len(alg_vars)

        model.fixed_vars = fixed_vars 
        model.n_nv = len(fixed_vars)


    def build_variable_locator(self, model, **kwargs):
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
        skip_validation = kwargs.pop('skip_validation', False)
        outlvl = kwargs.pop('outlvl', idaeslog.NOTSET)
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
            self.validate_steady_setpoint(set_point, self.s_mod)
            # ^ result should be that controller now has set point attributes
            # return result <- what did I mean by this?

        self.construct_objective_weight_matrices(self.c_mod,
                weight_overwrite=dynamic_weight_overwrite,
                tol=dynamic_weight_tol)

        self.add_objective_function(self.c_mod,
                control_penalty_type='action',
                name='tracking_objective')
#
#        ### TESTING PURPOSES ONLY ###
#        time = self.c_mod.time
#        for inp in self.c_mod.input_vars:
#            for t in time:
#                if t != time.first():
#                    inp[t].unfix()
#
#        assert (degrees_of_freedom(self.c_mod) == 
#                time.get_discretization_info()['nfe']*self.c_mod.n_iv)
#
#        results = self.default_solver.solve(self.c_mod, tee=True)
#
#        pdb.set_trace()


    def validate_steady_setpoint(self, set_point, steady_model, **kwargs):
        solver = kwargs.pop('solver', self.default_solver)
        outlvl = kwargs.pop('outlvl', self.outlvl)
        init_log = idaeslog.getInitLogger('nmpc', level=outlvl)
        solver_log = idaeslog.getSolveLogger('nmpc', level=outlvl)
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
            setattr(steady_model, attr_name, [])
            for _slice in getattr(self.c_mod, attr_name):
                if attr_name != 'scalar_vars':
                    vardata_t0 = _slice[t0]
                    vardata_t1 = _slice[t1]
                else:
                    vardata_t0 = _slice
                    vardata_t1 = _slice

                local_parent = steady_model
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
        steady_model.diff_sp = diff_var_sp
        steady_model.alg_sp = alg_var_sp
        steady_model.input_sp = input_var_sp

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
            init_log.error('Failed to solve for steady state set-point')
            raise ValueError

        c_diff_sp = []
        for var in steady_model.diff_vars:
            c_diff_sp.append(var[t0].value)
        self.c_mod.diff_sp = c_diff_sp

        c_input_sp = []
        for var in steady_model.input_vars:
            c_input_sp.append(var[t0].value)
        self.c_mod.input_sp = c_input_sp


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

        matrices = []
        # Attempt to construct weight for each type of setpoint
        # that we could allow
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
                    weight = 1/diff
                else:
                    weight = 1/tol
                matrix.append(weight)

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

        controls = model.input_vars
        R_entries = model.input_weights
        sp_controls = model.input_sp

        time = model.time

        state_term = sum(sum(Q_entries[i]*(states[i][t] - sp_states[i])**2
                for i in range(len(states)) if (Q_entries[i] is not None
                                            and sp_states[i] is not None))
                         for t in time)
        # Should be some check that Q_entries == None <=> sp_states == None
        # Not necessarily true, Q can have an entry due to an overwrite

        if control_penalty_type == 'error':
            control_term = sum(sum(R_entries[i]*(controls[i][t] - sp_controls[i])**2
                    for i in range(len(controls)) if (R_entries[i] is not None
                                                and sp_controls[i] is not None))
                               for t in time)
        elif control_penalty_type == 'action':
            if len(time) == 1:
                init_log.warning(
                        'Warning: Control action penalty specfied '
                        'for a model with a single time point.'
                        'Control term in objective function will be empty.')
            control_term = sum(sum(R_entries[i]*
                                  (controls[i][time[k]] - controls[i][time[k-1]])**2
                for i in range(len(controls)) if (R_entries[i] is not None
                                            and sp_controls[i] is not None))
                               for k in range(2, len(time)+1))

        obj_expr = state_term + control_term

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


    def add_pwc_constraints(self, **kwargs):
        model = kwargs.pop('model', self.c_mod)
        sample_time = kwargs.pop('sample_time', self.sample_time)
        outlvl = kwargs.pop('outlvl', self.outlvl)
        init_log = idaeslog.getInitLogger('nmpc', outlvl)
        init_log.info('Adding piecewise-constant constraints')

        time = model.time
        t0 = model.time.get_finite_elements()[0]
        t1 = model.time.get_finite_elements()[1]
        nfe_spacing = t1-t0
        nfe_per_sample = sample_time / nfe_spacing 
        if nfe_per_sample - int(nfe_per_sample) != 0:
            raise ValueError
        nfe_per_sample = int(nfe_per_sample)

        pwc_constraints = []

        for i, _slice in enumerate(model.input_vars):
            def pwc_rule(m, t):
                # Unless t is at the boundary of a sample, require
                # input[t] == input[t_next]
                if (t - time.first()) % sample_time == 0:
                    return Constraint.Skip

                t_next = time.next(t)
                return _slice[t_next] == _slice[t]

            name = '_pwc_input_' + str(i)
            pwc_constraint = Constraint(model.time, rule=pwc_rule)
            model.add_component(name, pwc_constraint)
            pwc_constraints.append(getattr(model, name))

        model._pwc_constraints = pwc_constraints


    def initialize_control_problem(self, **kwargs):
        strategy = kwargs.pop('strategy', 'from_previous')
        solver = kwargs.pop('solver', self.default_solver)
        solver_options = kwargs.pop('solver_options', None)
        if solver_options:
            solver.options = solver_options

        time = self.c_mod.time

        if strategy == 'from_previous':
            self.initialize_from_previous()

        elif strategy == 'simulate':
            input_type = kwargs.pop('inputs', 'set_point')
            if input_type == 'set_point':
                for i, _slice in enumerate(self.c_mod.input_vars):
                    for t in time:
                        if t != time.first():
                            _slice[t].fix(self.c_mod.input_sp[i])

            elif input_type == 'initial':
                for i, _slice in enumerate(self.c_mod.input_vars):
                    t0 = time.first()
                    for t in time:
                        _slice[t].fix(self.c_mod.input_vars[t0].value)

            # Deactivate objective function
            # strip bounds
            self.c_mod.tracking_objective.deactivate()
            for con in self.c_mod._pwc_constraints:
                con.deactivate()

            initialize_by_time_element(self.c_mod, self.c_mod.time,
                    solver=self.default_solver, outlvl=idaeslog.DEBUG)

            # Reactivate objective, pwc constraints, bounds

            self.c_mod.tracking_objective.activate()
            for con in self.c_mod._pwc_constraints:
                con.activate()

            for _slice in self.c_mod.input_vars:
                for t in time:
                    _slice[t].unfix()

            self.strip_controller_bounds.revert(self.c_mod)

        elif strategy == 'initial_conditions':
            self.initialize_from_initial_conditions(self)


    def initialize_from_previous_sample(self, model, sample_time, **kwargs):
        # Should only do this if controller is initialized
        # from a prior solve.
        if not self.controller_solved:
            raise ValueError

        attr_list = kwargs.pop('attr_list', 
                ['diff_vars', 'alg_vars', 'deriv_vars', 'input_vars'])
        # ^ why not just iterate over dae_vars here?
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
                        assert t_next in time
                        _slice[t].set_value(_slice[t_next].value)
                    # Othersie initialize to final steady state
                    # (Should provide some other options here,
                    # like copy inputs from set point then simulate)
                    else:
                        if not attrname == 'deriv_vars':
                            _slice[t].set_value(steady_varlist[i][steady_t])
                        else:
                            _slice[t].set_value(0)


    def initialize_from_simulation(self):
        raise NotImplementedError


    def initialize_from_initial_conditions(self):
        raise NotImplementedError

    
    def solve_control_problem(self, **kwargs):
        solver = kwargs.pop('solver', self.default_solver)
        outlvl = kwargs.pop('outlvl', self.outlvl) 
        s_log = idaeslog.getSolveLogger('nmpc', level=outlvl)

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


    def simulate_over_range(self, model, t_start, t_end, **kwargs):
        # should use knowledge of variables...

        # deactivate model except at t_start
        # deactivate disc. equations at t_start
        # solve for consistent 'initial' conditions at t_start

        # for each finite element in range... which are these?
        # get indices of finite elements that lie in range
        # t_start, t_end should correspond to finite elements

        # end by reactivating parts of model that were originally active
        outlvl = kwargs.pop('outlvl', idaeslog.NOTSET)
        init_log = idaeslog.getInitLogger('nmpc', outlvl)
        solver_log = idaeslog.getSolveLogger('nmpc', outlvl)

        time = model.time
        assert t_start in time.get_finite_elements()
        assert t_end in time.get_finite_elements()
        assert degrees_of_freedom(model) == 0
        ncp = model._ncp

        fe_in_range = [i for i, fe in enumerate(time.get_finite_elements())
                                if fe >= t_start and fe <= t_end]
        fe_in_range.pop(0)
        n_fe_in_range = len(fe_in_range)

        was_originally_active = get_activity_dict(model)

        non_initial_time = [t for t in time]
        non_initial_time.remove(t_start)
        deactivated = deactivate_model_at(model, time, non_initial_time,
                outlvl=idaeslog.ERROR)

        for i in fe_in_range:
            t_prev = time[(i-1)*ncp+1]

            fe = [time[k] for k in range((i-1)*ncp+1, i*ncp+2)]

            for t in fe:
                for comp in deactivated[t]:
                    if was_originally_active[id(comp)]:
                        comp.activate()

            for drv in model.deriv_vars:
                drv[t_prev].fix()
            for dv in model.diff_vars:
                dv[t_prev].fix()

            for t in fe:
                for _slice in model.dae_vars:
                    if not _slice[t].fixed:
                        _slice[t].set_value(_slice[t_prev].value)

            assert degrees_of_freedom(model) == 0

            with idaeslog.solver_log(solver_log, level=idaeslog.DEBUG) as slv:
                results = solver.solve(fs, tee=slc.tee)
            if results.solver.termination_condition == TerminationCondition.optimal:
                pass
            else:
                raise ValueError

            for t in fe:
                for comp in deactivated[t]:
                    comp.deactivate()

            for drv in model.deriv_vars:
                drv[t_prev].unfix()
            for dv in model.diff_vars:
                dv[t_prev].unfix()

        for t in time:
            for comp in deactivated[t]:
                if was_originally_active[id(comp)]:
                    comp.activate()
