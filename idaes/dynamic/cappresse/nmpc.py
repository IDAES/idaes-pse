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
        SolverFactory, Objective)
from pyomo.kernel import ComponentSet
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.dae.flatten import flatten_dae_variables

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.dyn_utils import (get_activity_dict, deactivate_model_at,
        path_from_block)
import idaes.logger as idaeslog

from collections import OrderedDict
import time as timemodule
import pdb

__author__ = "Robert Parker and David Thierry"


class VarLocator(object):
    """
    Class for storing information used to locate a VarData object.
    Used because I want to allow the user to supply set-point in terms
    of any variables they want. I then need to find these variables in the
    proper container.
    """

    def __init__(self, category, container, location):
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
        t1 = timemodule.time() 
        self.categorize_variables(self.c_mod, init_controller_inputs)
        t2 = timemodule.time()
        print(f'categorizing controller variables took {t2-t1} seconds')
        self.build_variable_locator(self.c_mod, 
                differential=self.c_mod.diff_vars,
                derivative=self.c_mod.deriv_vars,
                algebraic=self.c_mod.alg_vars,
                inputs=self.c_mod.input_vars,
                fixed=self.c_mod.fixed_vars)
        t3 = timemodule.time()
        print(f'building locator took {t3-t2} seconds')
        # Only expecting user arguments (set point) in form of control
        # variables, so only build locator for control variable model
        # for now.
        #
        # Not convinced that having a variable locator like this is the
        # best thing to do, but will go with it for now.


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
        fixed_vars = []

        scalar_vars, dae_vars = flatten_dae_variables(model, time)
        # Remove duplicates from these var lists
        dae_no_dups = list(OrderedDict([(id(v[t0]), v) for v in dae_vars]).values())
        scalar_no_dups = list(OrderedDict([(id(v), v) for v in scalar_vars]).values())

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
                fixed_vars.append(dae_no_dups.pop(i))
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

        model.fixed_vars = fixed_vars 
        model.n_nv = len(fixed_vars)


    def build_variable_locator(self, model, **kwargs):
        alg_list = kwargs.pop('algebraic', [])
        diff_list = kwargs.pop('differential', [])
        deriv_list = kwargs.pop('derivative', [])
        input_list = kwargs.pop('inputs', []) # note the s. input is reserved
        fixed_list = kwargs.pop('fixed', [])
        scalar_list = kwargs.pop('scalar', [])

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

        for i, _slice in enumerate(scalar_list):
            for t in model.time:
                locator[id(_slice[t])] = VarLocator('scalar', scalar_list, i)

        model.var_locator = locator


    def add_setpoint(self, set_point, **kwargs):
        skip_validation = kwargs.pop('skip_validation', False)
        if skip_validation:
            raise NotImplementedError(
                    'Maybe one day...')
            # In the implementation to 'skip validation', should still
            # get an error if a set_point var is not present in controller
            # model.
        else:
            steady_model = kwargs.pop('steady_model', None)
            if not steady_model:
                raise ValueError(
                   "'steady_model' required to validate set point")
            self.s_mod = steady_model
            self.validate_models(self.s_mod, self.p_mod)
            self.validate_steady_setpoint(set_point, self.s_mod)
            # return result


#    def _location_in_model(self, tgt_model, src_model, vardata, category=None):
#        """
#        Given a VarData object in a source model, attempts to find an object
#        of the same name in a target model. Then returns the object's index
#        in its target model container.
#
#        This uses target model's var_locator, so this must have been defined.
#        It is intended that the VarDatas in the source model have already been
#        categorized, so only the index is necessary to locate it in the target
#        model.
#
#        This can then be used to sort source model containers (e.g. alg_vars)
#        in the same order as the target model.
#
#        Has the dual purpose of validating the source model against the target
#        and sorting the source model's variable lists.
#        """
#        if vardata.model() is not src_model:
#            raise ValueError(
#            f'{vardata.name} is not a member of source model {src_model.name}')
#
#        local_parent = tgt_model
#        for r in path_from_block(vardata, src_model, include_comp=True):
#            local_parent = getattr(local_parent, r[0])[r[1]]
#        var_tgt = local_parent
#        info = tgt_model.var_locator[id(var_tgt)]
#
#        if category and info.category != category:
#            raise ValueError(
#                f'{var_tgt.name} is not in the expected category {category} '
#                'in model {tgt_model}')
#
#        return info.location


    def validate_steady_setpoint(self, set_point, steady_model, **kwargs):
        solver = kwargs.pop('solver', SolverFactory('ipopt'))
        outlvl = kwargs.pop('outlvl', idaeslog.NOTSET)
        init_log = idaeslog.getInitLogger(__name__, level=outlvl)
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
        # construct complete (vardata, value) lists for the set_point
        # argument in construct_objective_weight_matrices.
        
        diff_var_sp = [(var[t0], None) for var in steady_model.diff_vars]
        alg_var_sp = [(var[t0], None) for var in steady_model.alg_vars]
        input_var_sp = [(var[t0], None) for var in steady_model.input_vars]
        # ^ list entries are not actually vars, they are dictionaries...
        # maybe this is too confusing
        # If I know that these are ordered properly, why do I need to
        # provide the VarData object? I don't think I do...

        # This is where I map user values for set points into lists that
        # I can use to build the objective function (and weight matrices).
        for vardata, value in set_point:
            info = self.c_mod.var_locator[id(vardata)]
            category = info.category
            location = info.location
            if category == 'differential':
                diff_var_sp[location] = (steady_model.diff_vars[location],
                                         value)
            elif category == 'algebraic':
                alg_var_sp[location] = (steady_model.alg_vars[location],
                                         value)
            elif category == 'input':
                input_var_sp[location] = (steady_model.input_vars[location],
                                         value)

        # If I have the model and the model has attributes for set points
        # and variables (from which I can get initial values),
        # why do I need to add anything else to this function.
        # Still need an option to overwrite for the ambitious user...
        # If I'm not providing other arguments though, the format for
        # weight_overwrite is a bit more flexible
        diff_weights, alg_weights, input_weights = \
                self.construct_objective_weight_matrices(
                        diff_var_sp, alg_var_sp, input_var_sp,
                        overwrite=weight_overwrite,
                        tol=weight_tolerance)
        # Need to provide these arguments in order to calculate weights
        # for a variable number of varlists...
        # Alternative is to have a flag for include_algebraic
        # Or can just hasattr(alg_var_sp)...
        # (And also to just add weights as attributes to model,
        # not to return them...)

        # Add attributes (with proper names!) to steady model
        # so these values can be accessed later
        steady_model.diff_sp = diff_var_sp
        steady_model.alg_sp = alg_var_sp
        steady_model.input_sp = input_var_sp

        steady_model.diff_weights = diff_weights
        steady_model.alg_weights = alg_weights
        steady_model.input_weights = input_weights

        pdb.set_trace()


    def construct_objective_weight_matrices(self, *setpoints, **kwargs):
        """
        Do I even need the model? 
        Do I want to allow user to provide a setpoint here?

        Here each setpoint is a list of (VarData, value) tuples

        Overwrite is a lists of (VarData, value) tuples, where 
        this value will directly overwrite the weight of the
        corresponding VarData.
        """
        overwrite = kwargs.pop('overwrite', [])
        vars_to_overwrite = {id(tpl[0]): tpl[1] for tpl in overwrite}
        # dictionary mapping id of VarDatas to overwrite to specified weight

        tol = kwargs.pop('tol', 1e-6)

        matrices = []
        for sp in setpoints:
            # construct the (diagonal) matrix (list).
            matrix = [] 
            for vardata, value in sp:

                if id(vardata) in vars_to_overwrite:
                    weight = vars_to_overwrite[id(vardata)]
                    matrix.append(weight)
                    continue

                if value is None:
                    weight = None
                    matrix.append(weight)
                    continue

                diff = abs(vardata.value - value)
                if diff > tol:
                    weight = 1/diff
                else:
                    weight = 1/tol
                matrix.append(weight)

            matrices.append(matrix)
            
        return (matrix for matrix in matrices)


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
        name = kwargs.pop(name, 'objective')
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

        state_term = sum(sum(Q_entries[i]*(states[i][t] - sp_states[i])**2
                for i in range(len(states)) if (Q_entries[i] is not None 
                                            and sp_states[i] is not None))
                         for t in model.time)
        # Should be some check that Q_entries == None <=> sp_states == None

        control_term = sum(sum(R_entries[i]*(controls[i][t] = sp_controls[i])**2
                for i in range(len(controls)) if (R_entries[i] is not None
                                            and sp_controls[i] is not None))
                           for t in model.time)

        obj_expr = state_term + control_term

        obj = Objective(expr=obj_expr)
        model.add_component(name, obj)


