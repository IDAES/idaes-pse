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
'''
Some functions that are useful for initializing and otherwise working with
dynamic models

Author: Robert Parker
'''
from pyomo.environ import (Var, Block, Constraint, Set, Reference, value,
        Reals)
from pyomo.dae import DerivativeVar
from pyomo.dae.misc import get_index_information
from pyomo.core.kernel.component_set import ComponentSet
from pyomo.core.expr.visitor import identify_components
from pyomo.core.expr.visitor import identify_variables
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import (activated_equalities_generator,
        degrees_of_freedom)
from idaes.core.util.dyn_utils import (get_index_set_except,
        is_explicitly_indexed_by)
import pdb

def write_violated_equalities(m, tol=1e-5, filename=None):
    # m is any block...
    if not isinstance(filename, str):
        print('- - -\n Constraints violated:\n- - -')
        for con in activated_equalities_generator(m):
            upper_infeas = value(con.body) - value(con.upper)
            lower_infeas = value(con.lower) - value(con.body)
            infeas = max(upper_infeas, lower_infeas)
            if infeas > tol:
                print(con.name, infeas)
    else:
        with open(filename, 'w') as f:
            f.write('- - -\n Constraints violated:\n- - -\n')
            for con in activated_equalities_generator(m):
                upper_infeas = value(con.body) - value(con.upper)
                lower_infeas = value(con.lower) - value(con.body)
                infeas = max(upper_infeas, lower_infeas)
                if infeas > tol:
                    f.write(con.name + ': ' +  str(infeas) + '\n')

def remove_component_with_suffix(clist, suffix):
    lsuffix = len(suffix)
    to_remove = []
    for i, comp in enumerate(clist):
        if comp.name[-lsuffix:] == suffix:
           to_remove.append(i)
    for i in sorted(to_remove, reverse=True):
        del clist[i]

def remove_bounds_from(m):
    for var in ComponentSet(m.component_objects(Var)):
        if not isinstance(var, DerivativeVar):
            var.setub(None)
            var.setlb(None)
            var.domain = Reals

def get_differential_equations(m, time, t0=None):
    # Rename time cset
    '''
    Find/identify non-disc equations containing a derivative wrt time,
    call calculate_variable_from_constraint on corresponding DerivativeVar,
    equation at t = t0

    Assumes semi-explicit form of DAE, i.e. only one time derivative will
    appear in each differential equation, time derivative will appear linearly

    Does not distinguish between different orders of time derivatives or mixed
    time derivatives -- they are all assumed explicit and updated.
    '''
    diff_eqn_list = []
    deriv_list = []
    # Aside: can component_data_objects have repeats?
    for eqn in ComponentSet(m.component_data_objects(Constraint)):
        epc = eqn.parent_component()
        if eqn.name[-8:] == '_disc_eq':
            continue
        expr = eqn.expr
        for dv in identify_variables(expr):
            parent = dv.parent_component()
            if not isinstance(parent, DerivativeVar):
                continue
            if not time in ComponentSet(parent.get_continuousset_list()):
                continue
            if t0 is None:
                diff_eqn_list.append(eqn)
                deriv_list.append(dv)
            elif get_index_of_set(dv, time) == t0:
                diff_eqn_list.append(eqn)
                deriv_list.append(dv)

    return diff_eqn_list, deriv_list

def get_derivatives_at(b, time, t):
    '''
    Returns a list of VarDatas in block b that have derivatives
    with respect to time and time index t. No distinction made for
    multiple derivatives or mixed partials
    '''
    dvlist = []
    for var in ComponentSet(b.component_objects(Var)):
        # for DerivativeVars wrt time:
        if not isinstance(var, DerivativeVar):
            continue
        if time not in ComponentSet(var.get_continuousset_list()):
            continue
        n = var.index_set().dimen
        if n == 1:
            dvlist.append(var[t])
        elif n >= 2:
            info = get_index_set_except(var, time)        
            non_time_set = info['set_except']
            index_getter = info['index_getter']
            for non_time_index in non_time_set:
                index = index_getter(non_time_index, t)
                dvlist.append(var[index])
    return dvlist

def fix_initial_conditions(m, time, t0=0):
    '''
    Fix initial conditions of model to their current values.
    arguments are m, a block to operate on, and time,
    the set whose initial ( = t0) conditions will be fixed
    
    Want a way to skip at boundaries...
    does fixing differential variables specify the inputs? - no
    but if DVs at inlet boundaries (which are not actually differential
    variables) are specified, then the inputs must not be specified
    
    Can figure out what to not fix given a space domain set and a
    flow direction
    '''
    for var in ComponentSet(m.component_objects(Var)):
        # for DerivativeVars wrt time:
        if not isinstance(var, DerivativeVar):
            continue
        if time not in ComponentSet(var.get_continuousset_list()):
            continue
        state = var.get_state_var()

        n = state.index_set().dimen
        if n == 1:
            state[t0].fix()
        elif n >= 2:
            index_info = get_index_set_except(state, time)
            non_time_index_set = index_info['set_except']
            index_getter = index_info['index_getter']
            for non_time_index in non_time_index_set:
                index = index_getter(non_time_index, t0)
                state[index].fix()

