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
This module contains utility functions for dynamic IDAES models.
"""

from pyomo.environ import Block, Constraint, Var
from pyomo.dae import ContinuousSet, DerivativeVar
from pyomo.dae.set_utils import (is_explicitly_indexed_by,
        is_implicitly_indexed_by, get_index_set_except)

from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import ComponentSet
from idaes.core.util.exceptions import ConfigurationError
import idaes.logger as idaeslog
import pdb

__author__ = "Robert Parker"

# Set up logger
_log = idaeslog.getLogger(__name__)


def get_activity_dict(b):
    """
    args:
        b: a block to be searched for active components
    returns:
        a dictionary mapping id of constraints and blocks
        to a bool indicating if they are fixed
    """
    return {id(con): con.active 
                     for con in b.component_data_objects((Constraint, Block))}

def deactivate_model_at(b, cset, pt, outlvl=idaeslog.NOTSET):
    """
    Finds any block or constraint in block b, indexed explicitly or implicitly
    by cset, and deactivates each instance at pt
    Args:
        b - Block to search
        cset - ContinuousSet of interest
        pt - value, in ContinuousSet, to deactivate at
    Returns:
        deactivated - a list of component data that have been deactivated
    """
    if not pt in cset:
        msg = str(pt) + ' is not in ContinuousSet ' + cset.name
        raise ValueError(msg)
    deactivated = []
    
    for block in b.component_objects(Block):
        if is_explicitly_indexed_by(block, cset):
            info = get_index_set_except(block, cset)
            non_cset_set = info['set_except']
            index_getter = info['index_getter']
            for non_cset_index in non_cset_set:
                index = index_getter(non_cset_index, pt)
                try:
                    block[index].deactivate()
                    deactivated.append(block[index])
                except KeyError:
                    # except KeyError to allow Block.Skip
                    # TODO: use logger to give a warning here
                    msg = (block.name + ' has no index ' + str(index))
                    init_log = idaeslog.getInitLogger(b.name, outlvl)
                    init_log.warning(msg)
                    continue                

    for con in b.component_objects(Constraint):
        if (is_explicitly_indexed_by(con, cset) and
                not is_implicitly_indexed_by(con, cset)):
            info = get_index_set_except(con, cset)
            non_cset_set = info['set_except']
            index_getter = info['index_getter']
            for non_cset_index in non_cset_set:
                index = index_getter(non_cset_index, pt)
                try:
                    con[index].deactivate()
                    deactivated.append(con[index])
                except KeyError:
                    # except KeyError to allow Constraint.Skip
                    # TODO: use logger to give a warning here
                    msg = (con.name + ' has no index ' + str(index))
                    init_log = idaeslog.getInitLogger(b.name, outlvl)
                    init_log.warning(msg)
                    continue
                 
    return deactivated

# TODO: get_time_component_dict
# TODO: get_time_derivative_dict
#       will get components in one pass so I don't have to look through 
#       component_objects at each step of the integrator


def get_derivatives_at(b, time, t):
    """
    Finds derivatives with respect to time at point t.
    No distinction made for multiple derivatives or mixed partials.
    Args:
        b - Block to search for derivatives
        time - ContinuousSet to look for derivatives with respect to 
        t - point at which to return derivatives
    Returns:
        dvlist - list of derivatives found
    """
    dvlist = [] 
    for var in ComponentSet(b.component_objects(Var)):
        if not isinstance(var, DerivativeVar):
            continue
        if time not in set(var.get_continuousset_list()):
            continue

        info = get_index_set_except(var, time)
        non_time_set = info['set_except']
        index_getter = info['index_getter']
        for non_time_index in non_time_set:
            index = index_getter(non_time_index, t)
            dvlist.append(var[index])
    return dvlist

# TODO: should be able to replace this function everywhere
#       with getname and find_component
#       ^ not true. Cannot call find_component from a BlockData object
#                   Or on a name containing a decimal index
#       component looks like a similar substitue for BlockDatas, but
#       cannot seem to call on names including indices at all
def path_from_block(comp, blk, include_comp=False):
    """
    Returns a list of tuples with (local_name, index) pairs required
    to locate comp from blk
    Args:
        comp - Component(Data) to locate
        blk - Block(Data) to locate comp from
        include_comp - bool of whether or not to include the
                       local_name, index of the component itself
    Returns:
        route - A list of string, index tuples
    """
    parent_data = comp.parent_block()
    route = []
    # Walk up the hierarchy tree until blk is reached
    while parent_data != blk:
        parent_obj = parent_data.parent_component()
        # Record the local name and index required to locate the current block
        # from the parent_component of its parent_block
        # Pre-pend to the existing (name, index) list
        route = [(parent_obj.local_name, parent_data.index())] + route
        # If top-levelmodel has been reached, break
        # (This should not happen)
        if parent_data.parent_block() == None:
            break
        parent_data = parent_data.parent_block()

    # Append comp's name and index to the list if desired:
    if include_comp:
        route.append((comp.parent_component().local_name, comp.index()))
    return route

def copy_values_at_time(fs_tgt, fs_src, t_target, t_source, 
        copy_fixed=True, outlvl=idaeslog.NOTSET):
    """
    For all variables in fs_tgt (implicitly or explicitly) indexed by time,
    sets the value at t_target to that of the same variable in fs_src
    at t_source.

    Currently relies on fact that input blocks are flowsheets.
    Could extend to apply to non-flowsheet blocks, but would need to
    pass in the time set as well - could then apply to any set.

    Args:
        fs_tgt - target flowsheet
        fs_src - source flowsheet
        t_target - target time point
        t_source - source time point
        copy_fixed - bool of whether or not to copy over 
                     fixed variables in target model
    """
    time_target = fs_tgt.time
    for var_target in ComponentSet(fs_tgt.component_objects(Var)):
        if not is_explicitly_indexed_by(var_target, time_target):
            continue
        n = var_target.index_set().dimen

        local_parent = fs_src

        varname = var_target.getname(fully_qualified=True, relative_to=fs_tgt)
        # Calling find_component here makes the assumption that fs_src 
        # is not indexed
        var_source = fs_src.find_component(varname)
        if var_source is None:
            # Log a warning
            msg = varname + ' does not exist in ' + fs_src.name
            init_log = idaeslog.getInitLogger('copy_values logger', outlvl)
            init_log.warning(msg)
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
                var_target[target_index].set_value(
                        var_source[source_index].value)

    for blk_target in ComponentSet(fs_tgt.component_objects(Block)):
        if not is_explicitly_indexed_by(blk_target, time_target):
            continue
        n = blk_target.index_set().dimen

        blkname = blk_target.getname(fully_qualified=True, relative_to=fs_tgt)
        blk_source = fs_src.find_component(blkname)
        if blk_source is None:
            # log warning
            msg = blkname + ' does not exist in ' + fs_src.name
            init_log = idaeslog.getInitLogger('copy_values logger', outlvl)
            init_log.warning(msg)
            continue

        if n == 1:
            target_index = t_target
            source_index = t_source
            for var_target in ComponentSet(
                    blk_target[target_index].component_data_objects(Var)):
                if not copy_fixed and var_target.fixed:
                    continue
                # Here, find_component will not work from BlockData object
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
                    for r in path_from_block(var_target,
                                             blk_target[target_index]):
                        local_parent = getattr(local_parent, r[0])[r[1]]
                    var_source = getattr(local_parent,
                            var_target.parent_component().local_name)[
                                    var_target.index()]
                    var_target.set_value(var_source.value)

