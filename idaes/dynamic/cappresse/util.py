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
A module of helper functions for working with flattened DAE models.
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
        path_from_block, find_comp_in_block_at_time)
from idaes.core.util.initialization import initialize_by_time_element
from idaes.dynamic.cappresse.nmpc import find_comp_in_block
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

    def __init__(self, category, container, location, is_ic=False):
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
    assert len(varlist_tgt) == len(varlist_src)

    if not isinstance(t_tgt, list):
        t_tgt = [t_tgt]

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


def find_slices_in_model(tgt_model, src_model, tgt_locator, src_slices):
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
    t0_src = src_model.time.first()
    t0_tgt = tgt_model.time.first()
    tgt_slices = []
    for _slice in src_slices:
        init_var = _slice[t0_src]
#        tgt_var = find_comp_in_block(tgt_model, 
#                                     src_model, 
#                                     init_var)
        tgt_vardata = find_comp_in_block_at_time(tgt_model,
                                                 src_model,
                                                 init_var,
                                                 tgt_model.time,
                                                 t0_tgt)
        # This logic is bad because tgt_var might not be explicitly
        # time-indexed, or it may be indexed by things other than time
        #tgt_vardata = tgt_var[t0_tgt]

        try:
            tgt_container = tgt_locator[id(tgt_vardata)].container
        except KeyError:
            raise KeyError(
                'Locator does not seem to know about ' + 
                tgt_vardata.name)

        location = tgt_locator[id(tgt_vardata)].location
        tgt_slices.append(tgt_container[location])
    return tgt_slices


def simulate_over_range(self, model, t_start, t_end, **kwargs):
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
    # deactivate model except at t_start
    # deactivate disc. equations at t_start
    # solve for consistent 'initial' conditions at t_start

    # for each finite element in range... which are these?
    # get indices of finite elements that lie in range
    # t_start, t_end should correspond to finite elements

    # end by reactivating parts of model that were originally active

    # What variables does this function need knowledge of?

    solver = kwargs.pop('solver', self.default_solver)
    outlvl = kwargs.pop('outlvl', self.outlvl)
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

    # Deactivate model
    time_list = [t for t in time]
    deactivated = deactivate_model_at(model, time, time_list,
            outlvl=idaeslog.ERROR)

    for i in fe_in_range:
        t_prev = time[(i-1)*ncp+1]

        fe = [time[k] for k in range((i-1)*ncp+2, i*ncp+2)]

        for t in fe:
            for comp in deactivated[t]:
                if was_originally_active[id(comp)]:
                    comp.activate()

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

        was_fixed = {}
        for drv in model.deriv_vars:
            was_fixed[id(drv[t_prev])] = drv[t_prev].fixed
            drv[t_prev].fix()
        for dv in model.diff_vars:
            was_fixed[id(dv[t_prev])] = dv[t_prev].fixed
            dv[t_prev].fix()

        for t in fe:
            for _slice in model.dae_vars:
                if not _slice[t].fixed:
                # Fixed DAE variables are time-dependent disturbances,
                # whose values should not be altered by this function.
                #
                # ^ Exception is (true) initial conditions, but t should never
                # be time.first() because first point is always excluded from
                # fe. This assumption may need to change when/if forward
                # integration schemes are considered.
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

        for drv in model.deriv_vars:
            if not was_fixed[id(drv[t_prev])]:
                drv[t_prev].unfix()
        for dv in model.diff_vars:
            if not was_fixed[id(dv[t_prev])]:
                dv[t_prev].unfix()

    for t in time:
        for comp in deactivated[t]:
            if was_originally_active[id(comp)]:
                comp.activate()
