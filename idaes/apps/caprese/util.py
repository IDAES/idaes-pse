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
A module of helper functions for working with flattened DAE models.
"""

from pyomo.environ import (Block, Constraint, Var, TerminationCondition,
        SolverFactory, Objective, NonNegativeReals, Reals, 
        TransformationFactory)
from pyomo.common.collections import ComponentSet, ComponentMap
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.dae.flatten import flatten_dae_components
from pyomo.dae.set_utils import is_in_block_indexed_by
from pyomo.core.expr.visitor import identify_variables
from pyomo.core.base.constraint import _ConstraintData
from pyomo.core.base.block import _BlockData
from pyomo.core.base.var import _GeneralVarData
from pyomo.opt.solver import SystemCallSolver

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.dyn_utils import (get_activity_dict, 
                                       deactivate_model_at,
                                       path_from_block, 
                                       find_comp_in_block_at_time, 
                                       get_implicit_index_of_set,
                                       get_fixed_dict, 
                                       deactivate_constraints_unindexed_by, 
                                       find_comp_in_block)
from idaes.core.util.initialization import initialize_by_time_element
from idaes.apps.caprese.common.config import VariableCategory
import idaes.logger as idaeslog

from collections import OrderedDict
import random
import time as timemodule
import enum
import pdb

__author__ = "Robert Parker and David Thierry"


# TODO: clean up this file - add license, remove solver_available
# See if ipopt is available and set up solver
solver_available = SolverFactory('ipopt').available()
if solver_available:
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8,
                      'halt_on_ampl_error': 'yes'}
else:
    solver = None


class CachedVarsContext(object):
    def __init__(self, varlist, nvars, tlist):
        if type(tlist) is not list:
            tlist = [tlist]
        self.n_t = len(tlist)
        self.vars = varlist
        self.nvars = nvars
        self.tlist = tlist
        self.cache = [[None for j in range(self.n_t)] 
                for i in range(self.nvars)]

    def __enter__(self):
        for i in range(self.nvars):
            for j, t in enumerate(self.tlist):
                self.cache[i][j] = self.vars[i][t].value
        return self

    def __exit__(self, a, b, c):
        for i in range(self.nvars):
            for j, t in enumerate(self.tlist):
                self.vars[i][t].set_value(self.cache[i][j])


# NMPC Var could inherit from "DAE Var" - just add setpoint
# Inherit from Var? Advantage is that it could behave like Var
# fix(), setub() etc.
# Don't want to duplicate the VarData though - i.e. have them show
# up twice in component_data_objects
class NMPCVar(object):
    def __init__(self, _slice, category):
        self.timeslice = _slice
        # ^ _data dict: t -> vardata
        self.setpoint = None
        self.bounds = (None, None)
        self.category = category
        self.is_initial_condition = False


class NMPCVarGroup(object):
    # TODO: implement __iter__ that iterates over varlist
    def __init__(self, varlist, index_set, is_scalar=False):
        if type(varlist) is not list:
            raise TypeError(
                    'varlist argument must be a list')

        self.n_vars = len(varlist)
        self.varlist = varlist
        self.is_scalar = is_scalar
        if is_scalar:
            if index_set is not None:
                print('Warning. index_set provided for a list of scalars. '
                      'Setting to None.')
            self.index_set = None
            self.t0 = None
        else:
            self.index_set = self.validate_index_set(index_set)
            self.t0 = self.index_set.first()

        self.lb = self.n_vars*[None]
        self.ub = self.n_vars*[None]
        self.setpoint = self.n_vars*[None]
        self.reference = self.n_vars*[None]
        self.weights = self.n_vars*[None]
        # Can I get time set here? From an element (slice) of varlist
        # Complicated by the fact that varlist could be scalar_vars

    def __iter__(self):
        for var in self.varlist:
            yield var

    def __len__(self):
        return self.n_vars

    def __getitem__(self, i):
        return self.varlist[i]

    def validate_index_set(self, index_set):
        for var in self.varlist:
            # Hack so this doesn't fail for dicts that act as wrappers around
            # steady state model variables.
            # Should be able remove when I stop using a steady state model.
            if type(var) is dict:
                continue
            # For lack of a better error message, just assert:
            assert var.index_set() is index_set
        return index_set

    def validate_index(self, i):
        if i >= self.n_vars:
            raise ValueError(
                'Var index %i out of range. n_vars = %i' %(i, self.n_vars))
        if self.n_vars == 0:
            raise ValueError(
                'Index %i is invalid for an empty NMPCVarGroup' % i)

    def set_setpoint(self, i, val):
        self.validate_index(i)
        self.setpoint[i] = val

    def set_ub(self, i, val):
        self.validate_index(i)
        self.ub[i] = val
        if self.is_scalar:
            self.varlist[i].setub(val)
        else:
            for t in self.index_set:
                self.varlist[i][t].setub(val)
            #self.varlist[i][:].setub(val)
            # Reverted from using slices here as they don't work with
            # the steady model, which I (stupidly) store as dicts.
            # TODO: Change back one I get rid of steady model.

    def set_lb(self, i, val):
        self.validate_index(i)
        self.lb[i] = val
        if self.is_scalar:
            self.varlist[i].setlb(val)
        else:
            for t in self.index_set:
                self.varlist[i][t].setlb(val)
            #self.varlist[i][:].setlb(val)

    def set_domain(self, i, domain):
        self.validate_index(i)
        if self.is_scalar:
            self.varlist[i].domain = domain
        else:
            for t in self.index_set:
                self.varlist[i][t].domain = domain

    def set_value(self, i, val):
        self.validate_index(i)
        if self.is_scalar:
            self.varlist[i].set_value(val)
        else:
            for t in self.index_set:
                self.varlist[i][t].set_value(val)
            #self.varlist[i][:].set_value(val)

    def set_reference(self, i, val):
        self.validate_index(i)
        self.reference[i] = val

    def set_weight(self, i, val):
        self.validate_index(i)
        self.weights[i] = val


# This function is used as the domain for the user-provided
# list of inputs at time.first().
def validate_list_of_vardata(varlist):
    if not isinstance(varlist, list):
        raise TypeError('Not a list of VarData')
    for var in varlist:
        if not isinstance(var, _GeneralVarData):
            raise TypeError('Not a list of VarData')
    return varlist


def validate_list_of_vardata_value_tuples(varvaluelist):
    if not isinstance(varvaluelist, list):
        raise TypeError('Not a list')
    for item in varvaluelist:
        if not isinstance(item, tuple):
            raise TypeError('Item in list is not a tuple')
        if not len(item) == 2:
            raise ValueError('Tuple in list does not have correct length')
        if not isinstance(item[0], _GeneralVarData):
            raise TypeError('First entry is not a VarData')
        item = (item[0], float(item[1]))
    return varvaluelist


def validate_solver(solver):
    if type(solver) is str:
        solver = SolverFactory(solver)
    if not hasattr(solver, 'solve'):
        raise TypeError(
            'Solver does not implement solve method')
    return solver


class NMPCVarLocator(object):
    """
    Class for storing information used to locate a VarData object.
    Used because I want to allow the user to supply set-point in terms
    of any variables they want. I then need to find these variables in the
    proper container.
    """

    def __init__(self, category, group, location, is_ic=False, 
            is_measurement=False):
        """Constructor method. Assigns attributes based on arguments.

        Args:
            category : String describing a type of variable. Should be
                       something like 'differential' or 'algebraic'
            container : The list object that contains a VarData's
                        time slice.
            location : The index within the container where the variable
                       lives.
            is_ic : True if the variable is used as an initial condition.
        """
        # Should this class store the time index of the variable?
        # probably not (a. might not exist, b. should already be known)
        if category not in VariableCategory:
            raise TypeError(
            'category argument must be a valid VariableCategory enum item')
        self.category = category

        #if type(container) is not list:
        #    raise TypeError(
        #    'varlist argument must be a list')
        self.group = group

        if type(location) is not int:
            raise TypeError(
            'location argument must be an integer index')
        if location >= group.n_vars:
            raise ValueError(
            'location must be a valid index for the container') 
        self.location = location

        if type(is_ic) is not bool:
            raise TypeError(
                    'is_ic must be a bool')
        self.is_ic = is_ic

        if type(is_measurement) is not bool:
            raise TypeError(
                    'is_measurement must be a bool')
        self.is_measurement = is_measurement


def get_violated_bounds_at_time(group, timepoints, tolerance=1e-8):
    if type(timepoints) is not list:
        timepoints = [timepoints]
    violated = []
    for i, var in enumerate(group):
        ub = group.ub[i]
        lb = group.lb[i]
        if ub is not None:
            for t in timepoints:
                if var[t].value - ub > tolerance:
                    violated.append(var[t])
                    continue
        if lb is not None:
            for t in timepoints:
                if lb - var[t].value > tolerance:
                    violated.append(var[t])
                    continue
    return violated


def find_point_in_continuousset(point, cset, tolerance=1e-8):
    for t in cset:
        diff = abs(point-t)
        if diff < tolerance:
            return t
        if t > point:
            break
    return None


def copy_weights(tgt_group, src_group):
    assert tgt_group.n_vars == src_group.n_vars
    for i in range(tgt_group.n_vars):
        tgt_group.set_weight(i, src_group.weights[i])


def copy_values_at_time(varlist_tgt, varlist_src, t_tgt, t_src):
    """Copies values from time-indexed variables in one list, at one point
    in time to another list, at another point in time

    Args:
        varlist_tgt : List containing variables whose values will be set
        varlist_src : List containing variables whose values will be copied
        t_tgt : Point in time, or list of points in time, at which 
                variable values will be copied over
        t_src : Point in time from which variable values will be copied

    """
    # Downside to passing varlists as arguments directly is that I can't
    # validate that time points are valid for each model's time set
    # without just trying to access the VarDatas
    # ^ This could be circumvented by passing vargroups, then accessing
    # the groups' index_set attributes
    assert len(varlist_tgt) == len(varlist_src)

    if not isinstance(t_tgt, list):
        t_tgt = [t_tgt]

    # TODO? This could be made more compact with vargroups...
    for i, tgt_slice in enumerate(varlist_tgt):
        src_slice = varlist_src[i]

        try:
            src_value = src_slice[t_src].value
        except KeyError:
            raise KeyError(
                f'{t_src} does not seem to be a valid time index '
                'for the source variables')

        for t in t_tgt:
            try:
                var_tgt = tgt_slice[t]
            except KeyError:
                raise KeyError(
                    f'{t} does not seem to be a valid time index '
                    'for the target variables')

            var_tgt.set_value(src_value)


def find_slices_in_model(tgt_model, tgt_time, src_model, src_time, 
        tgt_locator, src_slices):
    """
    Given list of time-only slices in a source model and dictionary mapping
    VarData ids to VarLocator objects, attempts to find each slice
    in the target model and returns a list of the found slices in the same 
    order. 

    Args:
        tgt_model : Model to search for time-slices
        src_model : Model containing the slices to search for
        src_slices : List of time-only slices of variables in the source
                     model

    Returns:
        List of time-only slices to same-named variables in the target 
        model
    """
    # src_slice -> src_vardata -> name/route -> tgt_vardata -> 
    # tgt_slice (via locator) -> append to list 
    # (need both models to find_comp_in_block, and to get some point
    # in time at which to access slices)
    t0_src = src_time.first()
    t0_tgt = tgt_time.first()
    tgt_slices = []
    for _slice in src_slices:
        init_var = _slice[t0_src]
        tgt_vardata = find_comp_in_block_at_time(tgt_model,
                                                 src_model,
                                                 init_var,
                                                 tgt_time,
                                                 t0_tgt)
        # This logic is bad because tgt_var might not be explicitly
        # time-indexed, or it may be indexed by things other than time
        #tgt_vardata = tgt_var[t0_tgt]

        try:
            # locator is now a ComponentMap and takes in componentdatas
            tgt_container = tgt_locator[tgt_vardata].group.varlist
        except KeyError:
            raise KeyError(
                'Locator does not seem to know about ' + 
                tgt_vardata.name)

        location = tgt_locator[tgt_vardata].location
        tgt_slices.append(tgt_container[location])
    return tgt_slices


def initialize_by_element_in_range(model, time, t_start, t_end, 
        time_linking_vars=[],
        dae_vars=[],
        max_linking_range=0,
        **kwargs):
    """Function for solving a square model, time element-by-time element,
    between specified start and end times.

    Args:
        model : Flowsheet model to solve
        t_start : Beginning of timespan over which to solve
        t_end : End of timespan over which to solve

    Kwargs:
        solver : Solver option used to solve portions of the square model
        outlvl : idaes.logger output level
    """
    # TODO: How to handle config arguments here? Should this function
    # be moved to be a method of NMPC? Have a module-level config?
    # CONFIG, KWARGS: handle these kwargs through config

    solver = kwargs.pop('solver', SolverFactory('ipopt'))
    outlvl = kwargs.pop('outlvl', idaeslog.NOTSET)
    init_log = idaeslog.getInitLogger('nmpc', outlvl)
    solver_log = idaeslog.getSolveLogger('nmpc', outlvl)
    solve_initial_conditions = kwargs.pop('solve_initial_conditions', False)

    #TODO: Move to docstring
    # Variables that will be fixed for time points outside the finite element
    # when constraints for a finite element are activated.
    # For a "normal" process, these should just be differential variables
    # (and maybe derivative variables). For a process with a (PID) controller,
    # these should also include variables used by the controller.
    # If these variables are not specified, 

    # Timespan over which these variables will be fixed, counting backwards
    # from the first time point in the finite element (which will always be
    # fixed)
    # Should I specify max_linking_range as an integer number of finite
    # elements, an integer number of time points, or a float in actual time
    # units? Go with latter for now.

    # TODO: Should I fix scalar vars? Intuition is that they should already
    # be fixed.

    #time = model.time
    assert t_start in time.get_finite_elements()
    assert t_end in time.get_finite_elements()
    assert degrees_of_freedom(model) == 0

    #dae_vars = kwargs.pop('dae_vars', [])
    if not dae_vars:
        scalar_vars, dae_vars = flatten_dae_components(model, time, Var)
        for var in scalar_vars:
            var.fix()
        deactivate_constraints_unindexed_by(model, time)

    ncp = time.get_discretization_info()['ncp']

    fe_in_range = [i for i, fe in enumerate(time.get_finite_elements())
                            if fe >= t_start and fe <= t_end]
    t_in_range = [t for t in time if t >= t_start and t <= t_end]

    fe_in_range.pop(0)
    n_fe_in_range = len(fe_in_range)

    was_originally_active = get_activity_dict(model)
    was_originally_fixed = get_fixed_dict(model)

    # Deactivate model
    if not solve_initial_conditions:
        time_list = [t for t in time]
        deactivated = deactivate_model_at(model, time, time_list,
                outlvl=idaeslog.ERROR)
    else:
        time_list = [t for t in time if t != time.first()]
        deactivated = deactivate_model_at(model, time, time_list,
                outlvl=idaeslog.ERROR)

        assert degrees_of_freedom(model) == 0
        with idaeslog.solver_log(solver_log, level=idaeslog.DEBUG) as slc:
            results = solver.solve(model, tee=slc.tee)
        if results.solver.termination_condition == TerminationCondition.optimal:
            pass
        else:
            raise ValueError

        deactivated[time.first()] = deactivate_model_at(model, time, 
                time.first(),
                outlvl=idaeslog.ERROR)[time.first()]

    # "Integration" loop
    for i in fe_in_range:
        t_prev = time[(i-1)*ncp+1]

        fe = [time[k] for k in range((i-1)*ncp+2, i*ncp+2)]

        con_list = []
        for t in fe:
            # These will be fixed vars in constraints at t
            # Probably not necessary to record at what t
            # they occur
            for comp in deactivated[t]:
                if was_originally_active[id(comp)]:
                   comp.activate()
                   if not time_linking_vars:
                       if isinstance(comp, _ConstraintData):
                           con_list.append(comp)
                       elif isinstance(comp, _BlockData):
                           # Active here should be independent of whether block
                           # was active
                           con_list.extend(
                               list(comp.component_data_objects(Constraint,
                                                                 active=True)))

        if not time_linking_vars:
            fixed_vars = []
            for con in con_list:
                for var in identify_variables(con.expr,
                                              include_fixed=False):
                    # use var_locator/ComponentMap to get index somehow
                    t_idx = get_implicit_index_of_set(var, time)
                    if t_idx is None:
                        assert not is_in_block_indexed_by(var, time)
                        continue
                    if t_idx <= t_prev:
                        fixed_vars.append(var)
                        var.fix()
        else:
            fixed_vars = []
            time_range = [t for t in time 
                          if t_prev - t <= max_linking_range
                          and t <= t_prev]
            time_range = [t_prev]
            for _slice in time_linking_vars:
                for t in time_range:
                    #if not _slice[t].fixed:
                    _slice[t].fix()
                    fixed_vars.append(_slice[t])

        # Here I assume that the only variables that can appear in 
        # constraints at a different (later) time index are derivatives
        # and differential variables (they do so in the discretization
        # equations) and that they only participate at t_prev.
        #
        # This is not the case for, say, PID controllers, in which case
        # I should pass in a list of "complicating variables," then fix
        # them at all time points outside the finite element.
        #
        # Alternative solution is to identify_variables in each constraint
        # that is activated and fix those belonging to a previous finite
        # element. (Should not encounter variables belonging to a future
        # finite element.)
        # ^ This option is easier, less efficient
        #
        # In either case need to record whether variable was previously fixed
        # so I know if I should unfix it or not.

        for t in fe:
            for _slice in dae_vars:
                if not _slice[t].fixed:
                    # Fixed DAE variables are time-dependent disturbances,
                    # whose values should not be altered by this function.
                    _slice[t].set_value(_slice[t_prev].value)

        assert degrees_of_freedom(model) == 0

        with idaeslog.solver_log(solver_log, level=idaeslog.DEBUG) as slc:
            results = solver.solve(model, tee=slc.tee)
        if results.solver.termination_condition == TerminationCondition.optimal:
            pass
        else:
            raise ValueError

        for t in fe:
            for comp in deactivated[t]:
                comp.deactivate()

        for var in fixed_vars:
            if not was_originally_fixed[id(var)]:
                var.unfix()

    for t in time:
        for comp in deactivated[t]:
            if was_originally_active[id(comp)]:
                comp.activate()

    assert degrees_of_freedom(model) == 0


def add_noise_at_time(varlist, t_list, **kwargs):
    """Function to add random noise to the values of variables at a particular
    point in time.

    Args:
        varlist : List of (only) time-indexed variables (or References) to
                  which to add noise
        t_list : Point in time, or list of points in time, at which to add noise

    Kwargs:
        random_function : Function that will be called to get the random noise 
                          added to each variable
        args_function : Function that maps index (location) of a variable in
                      varlist to a list of arguments that can be passed to the
                      random function
        weights : List of weights for random distribution arguments such as 
                  standard deviation (Gaussian), radius (uniform), or lambda
                  (Laplacian)
        sigma_0 : Value of standard deviation for a variable with unit weight.
                  Default is 0.05.
        random_arg_dict : Dictionary containing other values users may want
                          to use in their argument function
        bound_strategy : String describing strategy for case in which a bound
                         is violated. Options are 'discard' (default) or 'push'.
        discard_limit : Number of discarded random values after which an
                        exception will be raised. Default is 5
        bound_push : Distance from bound if a push strategy is used for bound
                     violate. Default is 0

    Returns:
        A dictionary mapping each t in t_list to the 
    """
    n = len(varlist)
    rand_fcn = kwargs.pop('random_function', random.gauss)
    weights = kwargs.pop('weights', [1 for i in range(n)])
    random_arg_dict = kwargs.pop('random_arg_dict', {})
    assert len(weights) == n
    sig_0 = kwargs.pop('sigma_0', 0.05)
    sig = [w*sig_0 if w is not None else None for w in weights]

    args_fcn = kwargs.pop('args_function',
                          lambda i, val, **kwargs: [val, sig[i]] 
                                 if sig[i] is not None else None)

    bound_strategy = kwargs.pop('bound_strategy', 'discard')
    discard_limit = kwargs.pop('discard_limit', 5)
    bound_push = kwargs.pop('bound_push', 0)
    assert bound_push >= 0
    assert discard_limit >= 0

    if type(t_list) is not list:
        t_list = [t_list]

    nom_values = {t: [var[t].value for var in varlist] for t in t_list}
    for t in t_list:
        if any([val is None for val in nom_values[t]]):
            raise ValueError(
                    'Cannot apply noise to an uninitialized variable')

    def violated_bounds(var, t, val):
        if var[t].ub is not None:
            if val > var[t].ub:
                return ('upper', var[t].ub)
        if var[t].lb is not None:
            if val < var[t].lb:
                return ('lower', var[t].lb)
        return None

    for i, var in enumerate(varlist):
        for t in t_list:
            rand_args = args_fcn(i, var[t].value, **random_arg_dict)
            if not rand_args:
                # If a certain value is to be skipped,
                # args_fcn should return None
                continue
            newval = rand_fcn(*rand_args)

            violated = violated_bounds(var, t, newval)
            if not violated:
                var[t].set_value(newval)
                continue
            if bound_strategy == 'discard':
                for count in range(0, discard_limit):
                    newval = rand_fcn(*rand_args) 
                    if not violated_bounds(var, t, newval):
                        break
                if violated_bounds(var, t, newval):
                    raise ValueError(
                        'Discard limit exceeded when trying to apply noise to '
                        + var[t].name + ' with arguments ' + str(rand_args) +
                        '. Please adjust bounds or tighten distribution.')
            elif bound_strategy == 'push':
                if violated[0] == 'upper':
                    newval = violated[1] - bound_push
                elif violated[0] == 'lower':
                    newval = violated[1] + bound_push
                if violated_bounds(var, t, newval):
                    raise ValueError(
                            'Value after noise violates bounds even after '
                            'push. Please use a smaller bound push.')
            var[t].set_value(newval) 

    return nom_values
