#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Utility functions for scaling.

Author: Andrew Lee
"""

from copy import deepcopy
import sys

import json

from pyomo.environ import Block, Constraint, Suffix, value, Var
from pyomo.core.base.block import BlockData

import idaes.logger as idaeslog


_log = idaeslog.getLogger(__name__)

TAB = " " * 4


def get_scaling_suffix(component):
    """
    Get scaling suffix for component.

    If component is not a Block, gets scaling suffix from parent block.
    Creates a new suffix if one is not found.

    Args:
        component: component to get suffix for

    Returns:
        Pyomo scaling Suffix

    Raises:
        TypeError is component is an IndexedBlock
    """
    if isinstance(component, BlockData):
        blk = component
    elif isinstance(component, Block):
        raise TypeError(
            "IndexedBlocks cannot have scaling factors attached to them. "
            "Please assign scaling factors to the elements of the IndexedBlock."
        )
    else:
        blk = component.parent_block()

    try:
        sfx = blk.scaling_factor
    except AttributeError:
        # No existing suffix, create one
        _log.debug(f"Created new scaling suffix for {blk.name}")
        sfx = blk.scaling_factor = Suffix(direction=Suffix.EXPORT)

    return sfx


def scaling_factors_to_dict(blk_or_suffix, descend_into: bool = True):
    """
    Write scaling suffixes to a serializable json dict.

    Component objects are replaced by their local names so they can be
    serialized.

    Args:
        blk_or_suffix - Pyomo Block or Suffix object to covert to dict
        descend_into - for Blocks, whether to descend into any child blocks

    Returns
        dict containing scaling factors indexed by component names

    Raises:
        TypeError if blk_or_suffix is not an instance of Block or Suffix

    """
    # First, determine what type of component we have
    if isinstance(blk_or_suffix, Suffix):
        # Suffix
        sdict = _suffix_to_dict(blk_or_suffix)
        blk = blk_or_suffix.parent_block()
    elif isinstance(blk_or_suffix, BlockData):
        # Scalar block or element of indexed block
        sdict = _collect_block_suffixes(blk_or_suffix, descend_into=descend_into)
        blk = blk_or_suffix
    elif isinstance(blk_or_suffix, Block):
        # Indexed block
        blk = blk_or_suffix
        sdict = {}
        sdict["block_datas"] = {}
        for bd in blk_or_suffix.values():
            sdict["block_datas"][bd.name] = _collect_block_suffixes(
                bd, descend_into=descend_into
            )
    else:
        # Not a block or suffix
        raise TypeError(
            f"{blk_or_suffix.name} is not an instance of a Block of Suffix."
        )

    # Attach block name for future verification
    sdict["block_name"] = blk.name

    return sdict


def scaling_factors_from_dict(
    blk_or_suffix, json_dict: dict, overwrite: bool = False, verify_names: bool = True
):
    """
    Set scaling factors based on values in a serializable json dict.

    This method expects components to be referenced by their local names.

    Args:
        blk_or_suffix - Pyomo Block or Suffix object to set scaling factors on
        json_dict - dict of scaling factors to load into model
        overwrite - (bool) whether to overwrite existing scaling factors or not
        verify_names - (bool) whether to verify that all names in dict exist on block

    Returns
        None

    Raises:
        TypeError if blk_or_suffix is not an instance of Block or Suffix

    """
    # First, copy json_dict so we do not mutate original
    sdict = deepcopy(json_dict)
    # Pop block name for verification
    block_name = sdict.pop("block_name")

    # Next, determine what type of component we have
    if isinstance(blk_or_suffix, Suffix):
        # Suffix
        if verify_names and block_name != blk_or_suffix.parent_block().name:
            raise ValueError(
                f"Name of parent block ({blk_or_suffix.parent_block().name}) does "
                f"not match that recorded in json_dict ({block_name})"
            )
        _suffix_from_dict(
            blk_or_suffix, sdict, overwrite=overwrite, verify_names=verify_names
        )
    elif isinstance(blk_or_suffix, BlockData):
        # Scalar block or element of indexed block
        if verify_names and block_name != blk_or_suffix.name:
            raise ValueError(
                f"Block name ({blk_or_suffix.name}) does "
                f"not match that recorded in json_dict ({block_name})"
            )
        _set_block_suffixes_from_dict(
            blk_or_suffix, sdict, overwrite=overwrite, verify_names=verify_names
        )
    elif isinstance(blk_or_suffix, Block):
        # Indexed block
        if verify_names and block_name != blk_or_suffix.name:
            raise ValueError(
                f"Block name ({blk_or_suffix.name}) does "
                f"not match that recorded in json_dict ({block_name})"
            )
        for bd_name, bd_dict in sdict["block_datas"].items():
            bd = blk_or_suffix.parent_block().find_component(bd_name)
            _set_block_suffixes_from_dict(
                bd, bd_dict, overwrite=overwrite, verify_names=verify_names
            )
    else:
        # Not a block or suffix
        raise TypeError(
            f"{blk_or_suffix.name} is not an instance of a Block of Suffix."
        )


def scaling_factors_to_json_file(blk_or_suffix, filename: str):
    """
    Serialize scaling factors to file in json format.

    Args:
        blk_of_suffix - Block or Suffix to save scaling factors for
        filename - name of file to write to as string

    Returns:
        None

    Raises:
        TypeError if blk_or_suffix is not an instance of Block or Suffix
    """
    with open(filename, "w") as fd:
        json.dump(scaling_factors_to_dict(blk_or_suffix), fd, indent=3)


def scaling_factors_from_json_file(
    blk_or_suffix, filename: str, overwrite: bool = False, verify_names: bool = True
):
    """
    Load scaling factors from json file.

    Args:
        blk_of_suffix - Block or Suffix to load scaling factors for
        filename - name of file to load as string
        overwrite - (bool) whether to overwrite existing scaling factors or not
        verify_names - (bool) whether to verify that all names in dict exist on block

    Returns:
        None

    Raises:
        TypeError if blk_or_suffix is not an instance of Block or Suffix
    """
    with open(filename, "r") as f:
        scaling_factors_from_dict(
            blk_or_suffix, json.load(f), overwrite=overwrite, verify_names=verify_names
        )
    f.close()


def _collect_block_suffixes(block_data, descend_into=True):
    suffix = get_scaling_suffix(block_data)
    sdict = _suffix_to_dict(suffix)

    if descend_into:
        sdict["subblock_suffixes"] = {}
        for sb in block_data.component_data_objects(Block, descend_into=False):
            sdict["subblock_suffixes"][sb.local_name] = _collect_block_suffixes(
                sb, descend_into
            )

    return sdict


def _set_block_suffixes_from_dict(
    block_data, json_dict, verify_names=True, overwrite=False
):
    # First, copy dict so we can take it apart
    sdict = deepcopy(json_dict)

    # Pop any subblock suffixes
    sb_dict = sdict.pop("subblock_suffixes", None)

    # Set local suffix values
    suffix = get_scaling_suffix(block_data)
    _suffix_from_dict(suffix, sdict, verify_names=verify_names, overwrite=overwrite)

    if sb_dict is not None:
        # Get each subblock and apply function recursively
        for sb, sb_dict_value in sb_dict.items():
            subblock = block_data.find_component(sb)

            if subblock is not None:
                _set_block_suffixes_from_dict(
                    subblock,
                    sb_dict_value,
                    verify_names=verify_names,
                    overwrite=overwrite,
                )
            elif verify_names:
                raise AttributeError(
                    f"Block {block_data.name} does not have a subblock named {sb}."
                )


def _suffix_to_dict(suffix):
    sdict = {}

    for k, v in suffix.items():
        # Record components by their local name so we can use
        # find_Component to retrieve them later
        sdict[k.local_name] = v

    return sdict


def _suffix_from_dict(suffix, json_dict, verify_names=True, overwrite=False):
    parent_block = suffix.parent_block()

    for k, v in json_dict.items():
        # Safety catch in case this gets left in by other functions
        if k == "parent_name":
            continue

        comp = parent_block.find_component(k)
        if comp is not None:
            if overwrite or comp not in suffix:
                suffix[comp] = v
        elif verify_names:
            raise ValueError(
                f"Could not find component {k} on block {parent_block.name}."
            )


def get_scaling_factor(component):
    """
    Get scaling factor for component.

    Args:
        component: component to get scaling factor for

    Returns:
        float - scaling factor

    Raises:
        TypeError if component is a Block
    """
    if isinstance(component, (Block, BlockData)):
        raise TypeError("Blocks cannot have scaling factors.")

    try:
        sfx = get_scaling_suffix(component)
        return sfx[component]
    except (AttributeError, KeyError):
        # No scaling factor found, return None
        return None


def set_scaling_factor(component, scaling_factor: float, overwrite: bool = False):
    """
    Set scaling factor for component.

    Scaling factors must be positive, non-zero floats.

    Args:
        component - component to set scaling factor for
        scaling_factor - scaling factor to assign
        overwrite - (bool) whether to overwrite existing scaling factor

    Returns:
        None

    Raises:
        ValueError is scaling_factor is 0 or negative
    """
    # Cast to float to catch any garbage inputs
    scaling_factor = float(scaling_factor)

    # Check for negative or 0 scaling factors
    if scaling_factor < 0:
        raise ValueError(
            f"scaling factor for {component.name} is negative ({scaling_factor}). "
            "Scaling factors must be strictly positive."
        )
    elif scaling_factor == 0:
        raise ValueError(
            f"scaling factor for {component.name} is zero. "
            "Scaling factors must be strictly positive."
        )

    # Get suffix and assign scaling factor
    sfx = get_scaling_suffix(component)

    if not overwrite and component in sfx:
        _log.debug(
            f"Existing scaling factor for {component.name} found and overwrite=False. "
            "Scaling factor unchanged."
        )
    else:
        sfx[component] = scaling_factor


def del_scaling_factor(component, delete_empty_suffix: bool = False):
    """
    Delete scaling factor for component.

    Args:
        component - component to delete scaling factor for
        delete_empty_suffix - (bool) whether to delete scaling Suffix if it is
          empty after deletion.
    """
    # Get suffix
    parent = component.parent_block()
    sfx = get_scaling_suffix(parent)

    # Delete entry for component if it exists
    # Pyomo handles case where value does not exist in suffix with a no-op
    sfx.clear_value(component)

    if delete_empty_suffix:
        # Check if Suffix is empty (i.e. length 0)
        if len(sfx) == 0:
            # If so, delete suffix from parent block of component
            _log.debug(f"Deleting empty scaling suffix from {parent.name}")
            parent.del_component(sfx)


def report_scaling_factors(
    blk: Block, ctype=None, descend_into: bool = False, stream=None
):
    """
    Write the scaling factors for all components in a Block to a stream.

    Args:
        blk - Block to get scaling factors from.
        ctype - None, Var or Constraint. Type of component to show scaling factors for
          (if None, shows both Vars and Constraints).
        descend_into - whether to show scaling factors for components in sub-blocks.
        stream - StringIO object to write results to. If not provided, writes to stdout.

    Raises:
        TypeError if blk is not a Pyomo Block.
        ValueError is ctype is not None, Var or Constraint.
    """
    if stream is None:
        stream = sys.stdout

    if ctype not in [None, Var, Constraint]:
        raise ValueError(
            f"report_scaling_factors only supports None, Var or Constraint for argument ctype: "
            f"received {ctype}."
        )

    if not isinstance(blk, (Block, BlockData)):
        raise TypeError(
            "report_scaling_factors: blk must be an instance of a Pyomo Block."
        )

    stream.write(f"Scaling Factors for {blk.name}\n")

    # We will report Vars and Constraint is separate sections for clarity - iterate separately
    if ctype != Constraint:
        # Collect Vars
        vdict = {}
        for blkdata in blk.values():
            for vardata in blkdata.component_data_objects(
                Var, descend_into=descend_into
            ):
                val = vardata.value
                sf = get_scaling_factor(vardata)

                if sf is not None:
                    sfstr = "{:.3E}".format(sf)
                else:
                    sfstr = "None     "

                if val is not None:
                    vstr = "{:.3E}".format(val)
                    if sf is not None:
                        sval = "{:.3E}".format(value(vardata * sf))
                    else:
                        sval = vstr
                else:
                    vstr = "None     "
                    sval = "None"

                vdict[vardata.name] = (sfstr, vstr, sval)

        # Write Var section - skip if no Vars
        if len(vdict) > 0:
            # Get longest var name
            header = "Variable"
            maxname = len(max(vdict.keys(), key=len))
            if maxname < len(header):
                maxname = len(header)

            stream.write(
                f"\n{header}{' '*(maxname-len(header))}{TAB}Scaling Factor{TAB}Value{' '*4}{TAB}Scaled Value\n"
            )

            for n, i in vdict.items():
                # Pad name to length
                stream.write(
                    f"{n + ' '*(maxname-len(n))}{TAB}{i[0]}{' '*5}{TAB}{i[1]}{TAB}{i[2]}\n"
                )

    if ctype != Var:
        # Collect Constraints
        for blkdata in blk.values():
            cdict = {}
            for condata in blkdata.component_data_objects(
                Constraint, descend_into=descend_into
            ):
                sf = get_scaling_factor(condata)

                if sf is not None:
                    sfstr = "{:.3E}".format(sf)
                else:
                    sfstr = "None"

                cdict[condata.name] = sfstr

        # Write Constraint section - skip if no Constraints
        if len(cdict) > 0:
            # Get longest con name
            header = "Constraint"
            maxname = len(max(cdict.keys(), key=len))
            if maxname < len(header):
                maxname = len(header)

            stream.write(
                f"\n{header}{' ' * (maxname - len(header))}{TAB}Scaling Factor\n"
            )

            for n, i in cdict.items():
                # Pad name to length
                stream.write(f"{n + ' ' * (maxname - len(n))}{TAB}{i}\n")
