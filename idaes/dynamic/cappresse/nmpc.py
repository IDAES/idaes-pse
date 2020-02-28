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
        SolverFactory)
from pyomo.kernel import ComponentSet
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.dae.flatten import flatten_dae_variables

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.dyn_utils import (get_activity_dict, deactivate_model_at,
        path_from_block)
import idaes.logger as idaeslog

from collections import OrderedDict
import pdb

__author__ = "Robert Parker and David Thierry"


class NMPCSim(object):
    """
    Main class for NMPC simulations of IDAES flowsheets. 
    """

    def __init__(self, plant_model, controller_model, initial_inputs, **kwargs):
        # Maybe include a kwarg for require_steady - if False, set-point is not
        # forced to be a steady state

        solver = kwargs.pop('solver', SolverFactory('ipopt'))

        # Logger properties:
        outlvl = kwargs.pop('outlvl', idaeslog.NOTSET)

        # Set up attributes
        self.p_mod = plant_model
        self.c_mod = controller_model

        # Validate models
        self.validate_models(self.p_mod, self.c_mod)

        # Solve for consistent initial conditions
        self.p_mod._originally_active = get_activity_dict(self.p_mod)
        self.solve_initial_conditions(self.p_mod, initial_inputs,
                activity_dict=self.p_mod._originally_active,
                solver=solver)

        # Categorize variables in plant model
        self.categorize_variables(self.p_mod, initial_inputs)
        self.p_mod.initial_inputs = initial_inputs
        # ^ Add these as an attribute so they can be used later to validate
        # steady state model

        # Categorize variables in controller model
        init_controller_inputs = self.validate_inputs(controller_model,
                plant_model, initial_inputs)
        self.categorize_variables(self.c_mod, init_controller_inputs)


    def validate_inputs(self, tgt_model, src_model, src_inputs=None,
            outlvl=idaeslog.NOTSET):
        log = idaeslog.getInitLogger(__name__, level=outlvl)

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
        if not (isinstance(m1, FlowsheetBlock) and 
                isinstance(m2, FlowsheetBlock)):
            raise ValueError(
                    'Provided models must be FlowsheetBlocks')
        if m1.model() is m2.model():
            raise ValueError(
                    'Provided models must not live in the same top-level'
                    'ConcreteModel')
        return True


    def solve_initial_conditions(self, model, initial_inputs, **kwargs):
        # Record which Constraints/Variables are initially inactive
        # deactivate model except at t=0
        # fix initial inputs - raise error if no value
        # raise error if not square
        # re-activate model
        #
        # Later include option to skip solve for consistent initial conditions

        was_originally_active = kwargs.pop('activity_dict', None)
        solver = kwargs.pop('solver', SolverFactory('ipopt'))
        outlvl = kwargs.pop('outlvl', idaeslog.NOTSET)
        solver_log = idaeslog.getSolveLogger(__name__, level=outlvl)
        init_log = idaeslog.getInitLogger(__name__, level=outlvl)

        toplevel = model.model()

        non_initial_time = [t for t in model.time]
        non_initial_time.remove(model.time.first())
        deactivated = deactivate_model_at(model, model.time, non_initial_time, 
                outlvl=idaeslog.ERROR)

        for vardata in initial_inputs:
            if vardata.model() is not toplevel:
                raise ValueError(
                        f"Trying to fix an input that does not belong to model"
                        " {toplevel.name}. Are 'initial_inputs' arguments"
                        " contained in the plant model?")
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
                    "right number of initial conditions in the plant model been"
                    " fixed?")

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
        not_vars = []

        scalar_vars, dae_vars = flatten_dae_variables(model, time)
        # Remove duplicates from these var lists
        dae_no_dups = list(OrderedDict([(id(v[t0]), v) for v in dae_vars]).values())
        scalar_no_dups = list(OrderedDict([(id(v), v) for v in scalar_vars]).values())

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
            index = var[t0].index()
            if not isinstance(parent, DerivativeVar):
                continue
            if time not in ComponentSet(parent.get_continuousset_list()):
                continue
            assert var is dae_no_dups.pop(i)
            deriv_vars.append(var)
            init_dv_list.append(parent.get_state_var()[index])

        # Find differential variables
        for i, diffvar in reversed(list(enumerate(dae_no_dups))):
            for j, dv in enumerate(init_dv_list):
                # Don't need to reverse here as we intend to break the loop
                # after the first pop
                if diffvar[t0] is dv:
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
                not_vars.append(dae_no_dups.pop(i))
            else:
                alg_vars.append(dae_no_dups.pop(i))

        model.deriv_vars = deriv_vars
        model.diff_vars = diff_vars
        model.n_dv = len(diff_vars)
        assert model.n_dv == len(deriv_vars)

        model.input_vars = input_vars
        model.n_iv = len(input_vars)

        model.alg_vars = alg_vars
        model.n_av = len(alg_vars)

        model.not_vars = not_vars
        model.n_nv = len(not_vars)


    def add_setpoint(self, set_point, **kwargs):
        skip_validation = kwargs.pop('skip_validation', False)
        if skip_validation:
            raise NotImplementedError(
                    'Maybe one day...')
        else:
            steady_model = kwargs.pop('steady_model', None)
            if not steady_model:
                raise ValueError(
                   "'steady_model' required to validate set point")
            self.steady_model = steady_model
            self.validate_models(self.steady_model, self.p_mod)
            # validate set point
            # return result


    def validate_tracking_setpoint(self, set_point, steady_model):
        steady_inputs = self.validate_inputs(steady_model, self.p_mod)
        self.categorize_variables(steady_model, steady_inputs)
        # put set_point into the form of lists of proper dimension
        # will require map vardata -> (category, location)
        pass


    def construct_objective_function(self, **kwargs):
        # Read set point from arguments, if provided
        # State and control set poitns 
        set_point_state = kwargs.pop('set_point_state',
                kwargs.pop('set_point', None))
        set_point_control = kwargs.pop('set_point_control', None)

        # If no set point is provided,
        if not set_point_state:
            pass
        else:
            assert len(set_point_state) == self.c_mod.n_dv

        if not set_point_control:
            pass
        else:
            assert len(set_point_control) == self.c_mod.n_iv

        # Q and R are p.s.d. matrices that weigh the state and
        # control norms in the objective function
        Q_diagonal = kwargs.pop('Q_diagonal', True)
        R_diagonal = kwargs.pop('R_diagonal', True)

        Q_matrix = kwargs.pop('Q_matrix', None) 
        R_matrix = kwargs.pop('R_matrix', None)

        state_weight = kwargs.pop('state_weight', 1)
        control_weight = kwargs.pop('control_weight', 1)

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

        if not Q_matrix:
            # call function to assemble Q matrix from initial conditions 
            # and set point
            pass
        else:
            assert len(Q_matrix) == self.c_mod.n_dv
            # Q_matrix is an iterable of proper dimension.
            # Still require that we can use it to construct a Pyomo expression

        # Add objective to the model
