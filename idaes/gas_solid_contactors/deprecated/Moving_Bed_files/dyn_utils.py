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
from pyomo.environ import (Var, Block, value, Reals)
from pyomo.dae import DerivativeVar
from pyomo.dae.misc import get_index_information
from pyomo.core.kernel.component_set import ComponentSet
from idaes.core.util.model_statistics import activated_equalities_generator
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

def remove_bounds_from(m):
    for var in ComponentSet(m.component_objects(Var)):
        if not isinstance(var, DerivativeVar):
            var.setub(None)
            var.setlb(None)
            var.domain = Reals

def get_ik_from_index(s, val, tol=1e-8):
    # Gets fin. element and col. point of float index val
    # in discretized ContinuousSet s. 
    # Checks equality to within tol.
    # Used to facilitate using get_index_information when
    # index value, rather than (i,k) is known.
    info = s.get_discretization_info()
    nfe = info['nfe']
    ncp = info['ncp']
    for i in range(nfe):
        for k in range(ncp):
            location = i*ncp + k + 1
            if abs(s[location]-val) < 1e-8:
                return (i,k)
    if abs(s.last()-val) < 1e-8:
        return (nfe, 0)

    # return empty tuple if val could not be found within tolerance
    return ()

def is_indexed_by(comp, s):
    # Returns True if component comp is indexed by set s
    if not comp.is_indexed():
        return False
    n = comp.index_set().dimen
    if n == 1:
        if comp.index_set() is s:
            return True
        else:
            return False
    if n >= 2:
        if s in comp.index_set().set_tuple:
            return True
        else:
            return False

def is_implicitly_indexed_by(comp, s, stop_at=None):
    # Returns True if component comp or any of its parent blocks
    # are indexed by set s. Works by recursively checking parent
    # blocks.
    #
    # If block stop_at is provided, function will return False
    # if stop_at is reached, regardless of whether stop_at is 
    # indexed by s. Meant to be an "upper bound" for blocks to 
    # check, like a flowsheet.
    if is_indexed_by(comp, s):
        return True
    parent = comp.parent_block()
    # parent block could be block, blockdata, or none
    # want to check the block for a time index

    # assume here that top-level block cannot be indexed
    while not (parent is None):
        parent = parent.parent_component()
        if parent is stop_at:
            return False
        if is_indexed_by(parent, s):
            return True
        else:
            parent = parent.parent_block()
    return False

def get_index_set_except(comp, s):
    '''
    Returns a dictionary:
      'set_except'   -> Pyomo Set or SetProduct indexing comp, with s omitted.
      'index_getter' -> Function to return an index for comp given an index
                        from set_except and a value from set s.
                        Won't check if value is in s, so can be used to get
                        and index for a component that has a different s-set
    User should have already confirmed that s is an indexing set of comp.
    If combined with an interpolate function, could be used to project comp
    onto a hyperplane of constant s.
    '''
    info = {}
    # If indexed by one set:
    if comp.dim() == 1:
        info['set_except'] = [None] 
        info['index_getter'] = lambda incomplete_index, newval: newval
        return info

    # Otherwise find location of s within comp's index set
    count = 0
    location = 0
    other_ind_sets = []
    for ind_set in comp.index_set().set_tuple:
        if ind_set is s:
            location = count 
        else:
            other_ind_sets.append(ind_set)
            count += 1

    # other_ind_sets should not be empty
    # Should have length var.dim()-1

    if len(other_ind_sets) == 1:
        set_except = other_ind_sets[0]
    elif len(other_ind_sets) >= 2:
        set_except = other_ind_sets[0].cross(*other_ind_sets[1:])

    index_getter = (lambda incomplete_index, newval: 
            complete_index(location, incomplete_index, newval))

    info['set_except'] = set_except
    info['index_getter'] = index_getter
    return info
            
def complete_index(loc, index, newval):
    if not isinstance(index, tuple):
        index = (index,)
    return index[0:loc] + (newval,) + index[loc:]

def path_to(comp, include_comp=False):
    # Returns list of (local_name, index) tuples required
    # to locate comp in its top-level model
    #
    # Use this instead of find_component because find_component
    # differentiates between, for example, 0 and 0.0
    parent_data = comp.parent_block()
    route = []
    while parent_data.parent_block() != None:
        parent_obj = parent_data.parent_component()
        route = [(parent_obj.local_name, parent_data.index())] + route
        parent_data = parent_data.parent_block()
    if include_comp:
        route.append((comp.parent_component().local_name, comp.index()))
    return route

def path_from_block(comp, blk, include_comp=False):
    # Returns a list of tuples with (local_name, index) pairs required
    # to locate comp from blk
    #
    # blk should be a BlockData object
    parent_data = comp.parent_block()
    route = []
    while parent_data != blk:
        parent_obj = parent_data.parent_component()
        route = [(parent_obj.local_name, parent_data.index())] + route
        if parent_data.parent_block() == None:
            break
        parent_data = parent_data.parent_block()
    if include_comp:
        route.append((comp.parent_component().local_name, comp.index()))
    return route

def fix_initial_conditions(m, time):
    # Fix initial conditions of model to their current values.
    # arguments are m, a block to operate on, and time,
    # the set whose initial ( = 0) conditions will be fixed
    #
    # Could be extended to fix values at any point given
    # location of finite element and collocation point
    for var in ComponentSet(m.component_objects(Var)):
        # for DerivativeVars wrt time:
        if not isinstance(var, DerivativeVar):
            continue
        if time not in var.get_continuousset_list():
            continue
        state = var.get_state_var()

        n = state.index_set().dimen
        if n == 1:
            state[0].fix()
        elif n >= 2:
            index_info = get_index_information(state, time)
            non_time_index_set = index_info['non_ds']
            index_getter = index_info['index function']
            for non_time_index in non_time_index_set:
                index = index_getter(non_time_index, 0, 0)
                state[index].fix()

def copy_non_time_indexed_values(fs_tgt, fs_src, copy_fixed=True):
    # Copies values of non-time-indexed variables from src to tgt flowsheets.
    # These variables are both i) not explicitly indexed by time
    #                     and ii) not contained in a block indexed by time
    # May or may not want to copy over values that have been fixed
    # in the target flowsheet.
    time_tgt = fs_tgt.time
    for var_tgt in ComponentSet(fs_tgt.component_objects(Var,
                                            descend_into=False)):
        if is_indexed_by(var_tgt, time_tgt):
            continue
        var_src = fs_src.find_component(var_tgt.local_name)
        # ^ this find_component is fine because var_tgt is a Var not VarData
        # and its local_name is used
        for index in var_tgt:
            if not copy_fixed and var_tgt[index].fixed:
                continue
            var_tgt[index].set_value(var_src.value)

    for blk_tgt in ComponentSet(fs_tgt.component_objects(Block)):
        if is_implicitly_indexed_by(blk_tgt, time_tgt):
            continue
        # block is not even implicitly indexed by time
        for var_tgt in ComponentSet(blk_tgt.component_objects(Var,
                                                 descend_into=False)):
            if is_indexed_by(var_tgt, time_tgt):
                continue

            # can't used find_component(local_name) here because I need
            # the name of a variable with respect to a specific block
            local_parent = fs_src
            for r in path_from_block(var_tgt, fs_tgt):
                local_parent = getattr(local_parent, r[0])[r[1]]
            var_src = getattr(local_parent, var_tgt.local_name)

            for index in var_tgt:
                if not copy_fixed and var_tgt[index].fixed:
                    continue
                var_tgt[index].set_value(var_src[index].value)

def copy_values_at_time(fs_tgt, fs_src, t_target, t_source, copy_fixed=True):
    '''
    For all variables in fs_tgt (implicitly or explicitly) indexed by time,
    sets the value at t_target to that of the same variable in fs_src 
    at t_source

    Currently relies on fact that input blocks are flowsheets.
    Could extend to apply to non-flowsheet blocks, but would need to
    pass in the time set as well - could then apply to any set.


    '''
    time_target = fs_tgt.time
    for var_target in ComponentSet(fs_tgt.component_objects(Var)):
        # iterate over Vars (rather than VarDatas) to avoid repeating work

        # in this implementation, we don't care about variables
        # that aren't indexed by time
        #
        # note that indiscriminantly copying vars that aren't indexed by time
        # would lead to the incorrect copying of variables in blocks
        # that are indexed by time
        if not is_indexed_by(var_target, time_target):
            continue
        n = var_target.index_set().dimen

        local_parent = fs_src
        for r in path_from_block(var_target, fs_tgt):
            local_parent = getattr(local_parent, r[0])[r[1]]
        try:
            var_source = getattr(local_parent, var_target.local_name)
        except AttributeError:
            # if the variable does not exist in the source model, continue
            print()
            print('Warning while copying values:\nVariable: ' + var_target.name + 
                  '\ndoes not exist in source block: ' + fs_src.name)
            continue

        if n == 1:
            if not copy_fixed and var_target[t_target].fixed:
                continue
            var_target[t_target].set_value(var_source[t_source].value)
        elif n >= 2:
            index_info = get_index_set_except(var_target, time_target)
            non_time_index_set = index_info['set_except']
            index_getter = index_info['index_getter']
            for non_time_index in non_time_index_set:
                source_index = index_getter(non_time_index, t_source)
                target_index = index_getter(non_time_index, t_target)
                if not copy_fixed and var_target[target_index].fixed:
                    continue
                var_target[target_index].set_value(var_source[source_index].value)

    for blk_target in ComponentSet(fs_tgt.component_objects(Block)):
        if not is_indexed_by(blk_target, time_target):
            continue
        n = blk_target.index_set().dimen

        local_parent = fs_src 
        for r in path_from_block(blk_target, fs_tgt):
            local_parent = getattr(local_parent, r[0])[r[1]]
        blk_source = getattr(local_parent, blk_target.local_name)

        if n == 1:
            target_index = t_target
            source_index = t_source
            for var_target in ComponentSet(
                    blk_target[target_index].component_data_objects(Var)):
                if not copy_fixed and var_target.fixed:
                    continue
                local_parent = blk_source[source_index]
                for r in path_from_block(var_target, blk_target[target_index]):
                    local_parent = getattr(local_parent, r[0])[r[1]]
                var_source = getattr(local_parent, 
                        var_target.parent_component().local_name)[
                                var_target.index()]
                var_target.set_value(var_source.value)
        elif n >= 2:
            index_info = get_index_set_except(blk_target, time_target)
            non_time_index_set = index_info['set_except']
            index_getter = index_info['index_getter']
            for non_time_index in non_time_index_set:
                source_index = index_getter(non_time_index, t_source)
                target_index = index_getter(non_time_index, t_target)
                for var_target in ComponentSet(
                        blk_target[target_index].component_data_objects(Var)):
                    if not copy_fixed and var_target.fixed:
                        continue
                    local_parent = blk_source[source_index]
                    for r in path_from_block(var_target, blk_target[target_index]):
                        local_parent = getattr(local_parent, r[0])[r[1]]
                    var_source = getattr(local_parent, 
                            var_target.parent_component().local_name)[
                                    var_target.index()]
                    var_target.set_value(var_source.value)

def integrate_flowsheet(fs):
    time = fs.time
    tfe_list = time.get_finite_elements()
