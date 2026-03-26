#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
This module contains tools for setting and checking variable bounds based on model metadata.
"""

__author__ = "Alexander Dowling, Douglas Allan, Andrew Lee, Robby Parker, Ben Knueven"

from pyomo.environ import (
    Param,
    value,
    Var,
)
from pyomo.core.base.block import BlockData

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def get_valid_range_of_component(component):
    """
    Return the valid range for a component as specified in the model metadata.

    Args:
        component: Pyomo component to get valid range for

    Returns:
        valid range for component if found. This will either be a 2-tuple (low, high) or None.

    Raises:
        AttributeError if metadata object not found

    """
    # Get metadata for component
    parent = component.parent_block()

    try:
        if hasattr(parent, "params"):
            meta = parent.params.get_metadata().properties
        else:
            meta = parent.get_metadata().properties
    except AttributeError:
        raise AttributeError(f"Could not find metadata for component {component.name}")

    # Get valid range from metadata
    try:
        n, i = meta.get_name_and_index(component.parent_component().local_name)
        cmeta = getattr(meta, n)[i]
        valid_range = cmeta.valid_range
    except ValueError:
        # Assume no metadata for this property
        _log.debug(f"No metadata entry for component {component.name}; returning None")
        valid_range = None

    return valid_range


def set_bounds_from_valid_range(component, descend_into=True):
    """
    Set bounds on Pyomo components based on valid range recorded in model metadata.
    WARNING - this function will overwrite any bounds already set on the component/model.

    This function will iterate over component data objects in Blocks and indexed components.

    Args:
        component: Pyomo component to set bounds on. This can be a Block, Var or Param.
        descend_into: (optional) Whether to descend into components on child Blocks (default=True)

    Returns:
         None

    """
    if component.is_indexed():
        for k in component:
            set_bounds_from_valid_range(component[k])
    elif isinstance(component, BlockData):
        for i in component.component_data_objects(
            ctype=[Var, Param], descend_into=descend_into
        ):
            set_bounds_from_valid_range(i)
    elif not hasattr(component, "bounds"):
        raise TypeError(
            f"Component {component.name} does not have bounds. Only Vars and Params have bounds."
        )
    else:
        valid_range = get_valid_range_of_component(component)

        if valid_range is None:
            valid_range = (None, None)

        component.setlb(valid_range[0])
        component.setub(valid_range[1])


def list_components_with_values_outside_valid_range(component, descend_into=True):
    """
    Return a list of component objects with values outside the valid range specified in the model
    metadata.

    This function will iterate over component data objects in Blocks and indexed components.

    Args:
        component: Pyomo component to search for component outside of range on.
            This can be a Block, Var or Param.
        descend_into: (optional) Whether to descend into components on child Blocks (default=True)

    Returns:
         list of component objects found with values outside the valid range.
    """
    comp_list = []

    if component.is_indexed():
        for k in component:
            comp_list.extend(
                list_components_with_values_outside_valid_range(component[k])
            )
    elif isinstance(component, BlockData):
        for i in component.component_data_objects(
            ctype=[Var, Param], descend_into=descend_into
        ):
            comp_list.extend(list_components_with_values_outside_valid_range(i))
    else:
        valid_range = get_valid_range_of_component(component)

        if valid_range is not None:
            cval = value(component)
            if cval is not None and (cval < valid_range[0] or cval > valid_range[1]):
                comp_list.append(component)

    return comp_list
