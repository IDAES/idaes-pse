# -*- coding: utf-8 -*-
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
This module contains utility functions for reporting structural statistics of
IDAES models.
"""

__author__ = "Andrew Lee"

import sys

from pyomo.environ import Block, Constraint, Expression, Objective, Var, value
from pyomo.dae import DerivativeVar
from pyomo.core.expr import identify_variables
from pyomo.common.collections import ComponentMap, ComponentSet
from pyomo.common.deprecation import deprecation_warning
from pyomo.contrib.pynumero.interfaces.external_grey_box import ExternalGreyBoxBlock
from pyomo.contrib.pynumero.interfaces.external_grey_box_constraint import (
    ExternalGreyBoxConstraint,
)

import idaes.logger as idaeslog
from idaes.core.scaling import get_scaling_factor

_log = idaeslog.getLogger(__name__)


# -------------------------------------------------------------------------
# Generator to handle cases where the input is an indexed Block
# Indexed blocks do not have component_data_objects, so we need to iterate
# over the indexed block first.
def _iter_indexed_block_data_objects(block, ctype, active, descend_into):
    if block.is_indexed():
        for bd in block.values():
            for c in bd.component_data_objects(
                ctype=ctype, active=active, descend_into=descend_into
            ):
                yield c
    else:
        for c in block.component_data_objects(
            ctype=ctype, active=active, descend_into=descend_into
        ):
            yield c


# -------------------------------------------------------------------------
# Block methods
def total_blocks_set(block):
    """
    Method to return a ComponentSet of all Block components in a model.

    Args:
        block : model to be studied

    Returns:
        A ComponentSet including all Block components in block (including block
        itself)
    """
    total_blocks_set = ComponentSet(
        _iter_indexed_block_data_objects(
            block, ctype=Block, active=None, descend_into=True
        )
    )
    total_blocks_set.add(block)
    return total_blocks_set


def number_total_blocks(block):
    """
    Method to return the number of Block components in a model.

    Args:
        block : model to be studied

    Returns:
        Number of Block components in block (including block itself)
    """
    # +1 to include main model
    return (
        sum(
            1
            for _ in _iter_indexed_block_data_objects(
                block, ctype=Block, active=None, descend_into=True
            )
        )
        + 1
    )


def activated_blocks_set(block):
    """
    Method to return a ComponentSet of all activated Block components in a
    model.

    Args:
        block : model to be studied

    Returns:
        A ComponentSet including all activated Block components in block
        (including block itself)
    """
    block_set = ComponentSet()
    if block.active:
        block_set.add(block)
        for b in _iter_indexed_block_data_objects(
            block, ctype=Block, active=True, descend_into=True
        ):
            block_set.add(b)
    return block_set


def greybox_block_set(block):
    """
    Function to return ComponentSet of all Greybox Blocks components in a
    model.

    Args:
        block : model to be studied

    Returns:
        A ComponentSet including all GreyBox Block components in block
        (including block itself)
    """
    block_set = ComponentSet()
    for grey_box in activated_block_component_generator(
        block, ctype=ExternalGreyBoxBlock
    ):
        block_set.add(grey_box)

    return block_set


def activated_greybox_block_set(block):
    """
    Function to return ComponentSet of activated Greybox Blocks components in a
    model.

    Args:
        block : model to be studied

    Returns:
        A ComponentSet including all GreyBox Block components in block
        (including block itself)
    """
    block_set = ComponentSet()
    for grey_box in greybox_block_set(block):
        if grey_box.active:
            block_set.add(grey_box)

    return block_set


def deactivated_greybox_block_set(block):
    """
    Function to return ComponentSet of deactivated Greybox Blocks components in a
    model.

    Args:
        block : model to be studied

    Returns:
        A ComponentSet including all GreyBox Block components in block
        (including block itself)
    """
    return greybox_block_set(block) - activated_greybox_block_set(block)


def number_deactivated_greybox_block(block):
    """
    Function to return a Number of deactivated Greybox Blocks components in a
    model.

    Args:
        block : model to be studied

    Returns:
        number of deactivated greybox blocks
    """
    return len(deactivated_greybox_block_set(block))


def number_greybox_blocks(block):
    """
    Function to return a number of Greybox Blocks components in a
    model.

    Args:
        block : model to be studied

    Returns:
        number of activated greybox blocks
    """
    return len(greybox_block_set(block))


def number_activated_greybox_blocks(block):
    """
    Function to return a Number of activated Greybox Blocks components in a
    model.

    Args:
        block : model to be studied

    Returns:
        number of activated greybox blocks
    """
    return len(activated_greybox_block_set(block))


def number_activated_blocks(block):
    """
    Method to return the number of activated Block components in a model.

    Args:
        block : model to be studied

    Returns:
        Number of activated Block components in block (including block itself)
    """
    b = 0
    if block.active:
        b = 1
        b += sum(
            1
            for _ in _iter_indexed_block_data_objects(
                block, ctype=Block, active=True, descend_into=True
            )
        )
    return b


def deactivated_blocks_set(block):
    """
    Method to return a ComponentSet of all deactivated Block components in a
    model.

    Args:
        block : model to be studied

    Returns:
        A ComponentSet including all deactivated Block components in block
        (including block itself)
    """
    # component_data_objects active=False does not seem to work as expected
    # Use difference of total and active block sets
    return total_blocks_set(block) - activated_blocks_set(block)


def number_deactivated_blocks(block):
    """
    Method to return the number of deactivated Block components in a model.

    Args:
        block : model to be studied

    Returns:
        Number of deactivated Block components in block (including block
        itself)
    """
    # component_data_objects active=False does not seem to work as expected
    # Use difference of total and active block sets
    return number_total_blocks(block) - number_activated_blocks(block)


# -------------------------------------------------------------------------
# Basic Constraint methods
def total_constraints_set(block, include_greybox=True):
    """
    Method to return a ComponentSet of all Constraint components in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A ComponentSet including all Constraint components in block
    """
    if include_greybox:
        ctypes = (Constraint, ExternalGreyBoxConstraint)
    else:
        ctypes = (Constraint,)
    return ComponentSet(
        activated_block_component_generator(
            block, ctype=ctypes, include_greybox=include_greybox
        )
    )


def number_total_constraints(block, include_greybox=True):
    """
    Method to return the total number of Constraint components in a model.
    This will include the number of constraints provided by Greybox models using
    the number_activated_greybox_equalities function.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        Number of Constraint components in block
    """
    return sum(1 for _ in total_constraints_set(block, include_greybox=include_greybox))


def activated_constraints_generator(block, include_greybox=True):
    """
    Generator which returns all activated Constraint components in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A generator which returns all activated Constraint components block
    """
    if include_greybox:
        ctypes = (Constraint, ExternalGreyBoxConstraint)
    else:
        ctypes = (Constraint,)
    for c in activated_block_component_generator(
        block, ctype=ctypes, include_greybox=include_greybox
    ):
        if c.active:
            yield c


def activated_constraints_set(block, include_greybox=True):
    """
    Method to return a ComponentSet of all activated Constraint components in a
    model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A ComponentSet including all activated Constraint components in block
    """
    return ComponentSet(
        activated_constraints_generator(block, include_greybox=include_greybox)
    )


def number_activated_constraints(block, include_greybox=True):
    """
    Method to return the number of activated Constraint components in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        Number of activated Constraint components in block
    """
    return sum(
        1
        for _ in activated_constraints_generator(block, include_greybox=include_greybox)
    )


def deactivated_constraints_generator(block, include_greybox=True):
    """
    Generator which returns all deactivated Constraint components in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A generator which returns all deactivated Constraint components block
    """
    # GreyBox constraints are a special case, as their activity depends on the ExternalGreyBoxBlock
    # Using activated_block_component_generator will filter these out, so we need to handle them separately
    for c in activated_block_component_generator(block, ctype=Constraint):
        if not c.active:
            yield c

    if include_greybox:
        # we could potentially save time by looking only for deactivated greybox blocks here,
        # but for robustness we will be thorough
        for b in _iter_indexed_block_data_objects(
            block, ctype=ExternalGreyBoxBlock, active=None, descend_into=True
        ):
            for c in b.component_data_objects(
                ctype=ExternalGreyBoxConstraint, active=None, descend_into=False
            ):
                if not c.active:
                    yield c


def deactivated_constraints_set(block, include_greybox=True):
    """
    Method to return a ComponentSet of all deactivated Constraint components in
    a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A ComponentSet including all deactivated Constraint components in block
    """
    return ComponentSet(
        deactivated_constraints_generator(block, include_greybox=include_greybox)
    )


def number_deactivated_constraints(block, include_greybox=True):
    """
    Method to return the number of deactivated Constraint components in a
    model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        Number of deactivated Constraint components in block
    """
    return sum(
        1
        for _ in deactivated_constraints_generator(
            block, include_greybox=include_greybox
        )
    )


# -------------------------------------------------------------------------
# Equality Constraints
def total_equalities_generator(block, include_greybox=True):
    """
    Generator which returns all equality Constraint components in a model.

    Args:
        block : model to be studied

    Returns:
        A generator which returns all equality Constraint components block
    """
    for c in activated_block_component_generator(block, ctype=Constraint):
        if c.upper is not None and c.lower is not None and c.upper == c.lower:
            yield c

    # GreyBox constraints are a special case, as they are always equalities
    if include_greybox:
        for b in _iter_indexed_block_data_objects(
            block, ctype=ExternalGreyBoxBlock, active=None, descend_into=True
        ):
            for c in b.component_data_objects(
                ctype=ExternalGreyBoxConstraint, active=None, descend_into=False
            ):
                yield c


def total_equalities_set(block, include_greybox=True):
    """
    Method to return a ComponentSet of all equality Constraint components in a
    model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A ComponentSet including all equality Constraint components in block
    """
    return ComponentSet(
        total_equalities_generator(block, include_greybox=include_greybox)
    )


def number_total_equalities(block, include_greybox=True):
    """
    Method to return the total number of equality Constraint components in a
    model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        Number of equality Constraint components in block
    """
    return sum(
        1 for _ in total_equalities_generator(block, include_greybox=include_greybox)
    )


def activated_equalities_generator(block, include_greybox=True):
    """
    Generator which returns all activated equality Constraint components in a
    model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A generator which returns all activated equality Constraint components
        block
    """
    for c in _iter_indexed_block_data_objects(
        block, Constraint, active=True, descend_into=True
    ):
        if (
            c.upper is not None
            and c.lower is not None
            and value(c.upper) == value(c.lower)
        ):
            yield c

    # GreyBox constraints are a special case, as they are always equalities
    if include_greybox:
        for b in _iter_indexed_block_data_objects(
            block, ctype=ExternalGreyBoxBlock, active=True, descend_into=True
        ):
            for c in b.component_data_objects(
                ctype=ExternalGreyBoxConstraint, active=True, descend_into=False
            ):
                yield c


def activated_equalities_set(block, include_greybox=True):
    """
    Method to return a ComponentSet of all activated equality Constraint
    components in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A ComponentSet including all activated equality Constraint components
        in block
    """
    return ComponentSet(
        activated_equalities_generator(block, include_greybox=include_greybox)
    )


def number_activated_equalities(block, include_greybox=True):
    """
    Method to return the number of activated equality Constraint components in
    a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        Number of activated equality Constraint components in block
    """
    return sum(
        1
        for _ in activated_equalities_generator(block, include_greybox=include_greybox)
    )


def number_activated_greybox_equalities(block) -> int:
    """
    Function to compute total number of equality constraints for all GreyBox objects in this block.

    Args:
        block : pyomo concrete model or pyomo block

    Returns:
        Number of equality constraints in all GreyBox objects on the provided block
    """
    equalities = 0
    for grey_box in activated_greybox_block_set(block):
        for _ in grey_box.component_data_objects(
            ctype=ExternalGreyBoxConstraint, active=None, descend_into=False
        ):
            equalities += 1
    return equalities


def number_deactivated_greybox_equalities(block) -> int:
    """
    Function to compute total number of equality constraints for all GreyBox objects in this block.

    Args:
        block : pyomo concrete model or pyomo block

    Returns:
        Number of equality constraints in all GreyBox objects on the provided block
    """
    equalities = 0
    for grey_box in deactivated_greybox_block_set(block):
        for _ in grey_box.component_data_objects(
            ctype=ExternalGreyBoxConstraint, active=None, descend_into=False
        ):
            equalities += 1
    return equalities


def deactivated_equalities_generator(block, include_greybox=True):
    """
    Generator which returns all deactivated equality Constraint components in a
    model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A generator which returns all deactivated equality Constraint
        components block
    """
    for c in total_equalities_generator(block, include_greybox=include_greybox):
        if not c.active:
            yield c


def deactivated_equalities_set(block, include_greybox=True):
    """
    Method to return a ComponentSet of all deactivated equality Constraint
    components in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A ComponentSet including all deactivated equality Constraint components
        in block
    """
    return ComponentSet(
        deactivated_equalities_generator(block, include_greybox=include_greybox)
    )


def number_deactivated_equalities(block, include_greybox=True):
    """
    Method to return the number of deactivated equality Constraint components
    in a model. This will include the number of deactivated equality constraints in Greybox models.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        Number of deactivated equality Constraint components in block
    """
    return sum(
        1
        for _ in deactivated_equalities_generator(
            block, include_greybox=include_greybox
        )
    )


# -------------------------------------------------------------------------
# Inequality Constraints
def total_inequalities_generator(block):
    """
    Generator which returns all inequality Constraint components in a
    model.

    Args:
        block : model to be studied

    Returns:
        A generator which returns all inequality Constraint components block
    """
    for c in activated_block_component_generator(block, ctype=Constraint):
        if c.upper is None or c.lower is None:
            yield c


def total_inequalities_set(block):
    """
    Method to return a ComponentSet of all inequality Constraint components in
    a model.

    Args:
        block : model to be studied

    Returns:
        A ComponentSet including all inequality Constraint components in block
    """
    return ComponentSet(total_inequalities_generator(block))


def number_total_inequalities(block):
    """
    Method to return the total number of inequality Constraint components in a
    model.

    Args:
        block : model to be studied

    Returns:
        Number of inequality Constraint components in block
    """
    return sum(1 for _ in total_inequalities_generator(block))


def activated_inequalities_generator(block):
    """
    Generator which returns all activated inequality Constraint components in a
    model.

    Args:
        block : model to be studied

    Returns:
        A generator which returns all activated inequality Constraint
        components block
    """
    for c in _iter_indexed_block_data_objects(
        block, Constraint, active=True, descend_into=True
    ):
        if c.upper is None or c.lower is None:
            yield c


def activated_inequalities_set(block):
    """
    Method to return a ComponentSet of all activated inequality Constraint
    components in a model.

    Args:
        block : model to be studied

    Returns:
        A ComponentSet including all activated inequality Constraint components
        in block
    """
    return ComponentSet(activated_inequalities_generator(block))


def number_activated_inequalities(block):
    """
    Method to return the number of activated inequality Constraint components
    in a model.

    Args:
        block : model to be studied

    Returns:
        Number of activated inequality Constraint components in block
    """
    return sum(1 for _ in activated_inequalities_generator(block))


def deactivated_inequalities_generator(block):
    """
    Generator which returns all deactivated inequality Constraint components in
    a model.

    Args:
        block : model to be studied

    Returns:
        A generator which returns all deactivated equality Constraint
        components block
    """
    for c in total_inequalities_generator(block):
        if not c.active:
            yield c


def deactivated_inequalities_set(block):
    """
    Method to return a ComponentSet of all deactivated inequality Constraint
    components in a model.

    Args:
        block : model to be studied

    Returns:
        A ComponentSet including all deactivated inequality Constraint
        components in block
    """
    return ComponentSet(deactivated_inequalities_generator(block))


def number_deactivated_inequalities(block):
    """
    Method to return the number of deactivated inequality Constraint components
    in a model.

    Args:
        block : model to be studied

    Returns:
        Number of deactivated inequality Constraint components in block
    """
    return sum(1 for _ in deactivated_inequalities_generator(block))


# -------------------------------------------------------------------------
# Basic Variable Methods
# Always use ComponentSets for Vars to avoid duplication of References
# i.e. number methods should always use the ComponentSet, not a generator
def variables_generator(block, include_greybox=True, active=True):
    """
    Generator which returns all Var components in a model.

    Args:
        block : model to be studied
        include_greybox : whether to include greybox variables
        active : whether to include only active sub-blocks (default = True)

    Returns:
        A generator which returns all Var components in block
    """
    for var in _iter_indexed_block_data_objects(
        block, ctype=Var, active=active, descend_into=True
    ):
        yield var

    if include_greybox:
        for egb in _iter_indexed_block_data_objects(
            block, ctype=ExternalGreyBoxBlock, active=None, descend_into=True
        ):
            for v in egb.component_data_objects(
                ctype=Var, active=None, descend_into=False
            ):
                yield v


def variables_set(block, include_greybox=True, active=True):
    """
    Method to return a ComponentSet of all Var components in a model.

    Args:
        block : model to be studied
        include_greybox : whether to include greybox variables
        active : whether to include only active sub-blocks (default = True)

    Returns:
        A ComponentSet including all Var components in block
    """
    var_set = ComponentSet()
    for var in variables_generator(
        block, include_greybox=include_greybox, active=active
    ):
        var_set.add(var)
    return var_set


def number_variables(block, include_greybox=True, active=True):
    """
    Method to return the number of Var components in a model.

    Args:
        block : model to be studied
        include_greybox : whether to include greybox variables
        active : whether to include only active sub-blocks (default = True)

    Returns:
        Number of Var components in block
    """
    return len(variables_set(block, include_greybox=include_greybox, active=active))


def fixed_variables_generator(block, include_greybox=True):
    """
    Generator which returns all fixed Var components in a model.

    Args:
        block : model to be studied
        include_greybox : whether to include greybox variables

    Returns:
        A generator which returns all fixed Var components block
    """
    for v in variables_generator(block, include_greybox=include_greybox):
        if v.fixed:
            yield v


def fixed_variables_set(block, include_greybox=True):
    """
    Method to return a ComponentSet of all fixed Var components in a model.

    Args:
        block : model to be studied
        include_greybox : whether to include greybox variables

    Returns:
        A ComponentSet including all fixed Var components in block
    """
    return ComponentSet(
        fixed_variables_generator(block, include_greybox=include_greybox)
    )


def number_fixed_variables(block, include_greybox=True):
    """
    Method to return the number of fixed Var components in a model.

    Args:
        block : model to be studied
        include_greybox : whether to include greybox variables

    Returns:
        Number of fixed Var components in block
    """
    return len(fixed_variables_set(block, include_greybox=include_greybox))


def unfixed_variables_generator(block, include_greybox=True):
    """
    Generator which returns all unfixed Var components in a model.

    Args:
        block : model to be studied
        include_greybox : whether to include greybox variables

    Returns:
        A generator which returns all unfixed Var components block
    """
    for v in variables_generator(block, include_greybox=include_greybox):
        if not v.fixed:
            yield v


def unfixed_variables_set(block, include_greybox=True):
    """
    Method to return a ComponentSet of all unfixed Var components in a model.

    Args:
        block : model to be studied
        include_greybox : whether to include greybox variables

    Returns:
        A ComponentSet including all unfixed Var components in block
    """
    return ComponentSet(
        unfixed_variables_generator(block, include_greybox=include_greybox)
    )


def number_unfixed_variables(block, include_greybox=True):
    """
    Method to return the number of unfixed Var components in a model.

    Args:
        block : model to be studied
        include_greybox : whether to include greybox variables

    Returns:
        Number of unfixed Var components in block
    """
    return len(unfixed_variables_set(block, include_greybox=include_greybox))


def variables_near_bounds_generator(
    block,
    tol=None,
    relative=None,
    skip_lb=False,
    skip_ub=False,
    abs_tol=1e-4,
    rel_tol=1e-4,
    include_greybox=True,
    apply_scaling=True,
):
    """
    Generator which returns all Var components in a model which have a value
    within tol (default: relative) of a bound.

    Args:
        block : model to be studied
        abs_tol : absolute tolerance for inclusion in generator (default = 1e-4)
        rel_tol : relative tolerance for inclusion in generator (default = 1e-4)
        skip_lb: Boolean to skip lower bound (default = False)
        skip_ub: Boolean to skip upper bound (default = False)
        include_greybox : whether to include greybox variables (default = True)
        apply_scaling: whether to apply scaling factors to the variable when determining if it is near a bound (default = True)

    Returns:
        A generator which returns all Var components block that are close to a
        bound
    """
    # Check for deprecated arguments
    if relative is not None:
        msg = (
            "variables_near_bounds_generator has deprecated the relative argument. "
            "Please set abs_tol and rel_tol arguments instead."
        )
        deprecation_warning(msg=msg, logger=_log, version="2.2.0", remove_in="2.11.0")
    if tol is not None:
        msg = (
            "variables_near_bounds_generator has deprecated the tol argument. "
            "Please set abs_tol and rel_tol arguments instead."
        )
        deprecation_warning(msg=msg, logger=_log, version="2.2.0", remove_in="2.11.0")
        # Set tolerances using the provided value
        abs_tol = tol
        rel_tol = tol

    for v in variables_generator(block, include_greybox=include_greybox):
        # To avoid errors, check that v has a value
        if v.value is None:
            continue
        sf = get_scaling_factor(v, default=1, warning=False) if apply_scaling else 1
        # First, magnitude of variable
        if v.ub is not None and v.lb is not None:
            # Both upper and lower bounds, apply tol to (upper - lower)
            mag = value(v.ub - v.lb)
        elif v.ub is not None:
            # Only upper bound, apply tol to bound value
            mag = abs(value(v.ub))
        elif v.lb is not None:
            # Only lower bound, apply tol to bound value
            mag = abs(value(v.lb))
        else:
            mag = 0

        # Calculate largest tolerance from absolute and relative
        tol = max(abs_tol / sf, mag * rel_tol)

        if v.ub is not None and not skip_ub and value(v.ub - v.value) <= tol:
            yield v
        elif v.lb is not None and not skip_lb and value(v.value - v.lb) <= tol:
            yield v


def variables_near_bounds_set(
    block,
    tol=None,
    relative=None,
    skip_lb=False,
    skip_ub=False,
    abs_tol=1e-4,
    rel_tol=1e-4,
    include_greybox=True,
    apply_scaling=True,
):
    """
    Method to return a ComponentSet of all Var components in a model which have
    a value within tolerance of a bound.

    Args:
        block : model to be studied
        abs_tol : absolute tolerance for inclusion in generator (default = 1e-4)
        rel_tol : relative tolerance for inclusion in generator (default = 1e-4)
        skip_lb: Boolean to skip lower bound (default = False)
        skip_ub: Boolean to skip upper bound (default = False)
        include_greybox : whether to include greybox variables (default = True)
        apply_scaling: whether to apply scaling factors to the variable when determining
        if it is near a bound (default = True)

    Returns:
        A ComponentSet including all Var components block that are close to a
        bound
    """
    return ComponentSet(
        variables_near_bounds_generator(
            block,
            tol,
            relative,
            skip_lb,
            skip_ub,
            abs_tol,
            rel_tol,
            include_greybox=include_greybox,
            apply_scaling=apply_scaling,
        )
    )


def number_variables_near_bounds(
    block, tol=None, abs_tol=1e-4, rel_tol=1e-4, include_greybox=True
):
    """
    Method to return the number of all Var components in a model which have
    a value within tol (relative) of a bound.

    Args:
        block : model to be studied
        abs_tol : absolute tolerance for inclusion in generator (default = 1e-4)
        rel_tol : relative tolerance for inclusion in generator (default = 1e-4)
        include_greybox : whether to include greybox variables (default = True)

    Returns:
        Number of components block that are close to a bound
    """
    return len(
        variables_near_bounds_set(
            block,
            tol=tol,
            abs_tol=abs_tol,
            rel_tol=rel_tol,
            include_greybox=include_greybox,
        )
    )


def variables_violating_bounds_generator(
    block: Block,
    abs_tol: float = 1e-4,
    rel_tol: float = 1e-4,
    include_greybox: bool = True,
    apply_scaling: bool = True,
):
    """
    Generator which returns all Var components in a model which have a value
    that violates their bounds.

    Args:
        block : model to be studied
        abs_tol : absolute tolerance for violation of bounds
        rel_tol : relative tolerance for violation of bounds
        include_greybox : whether to include greybox variables (default = True)
        apply_scaling: whether to apply scaling factors to the variable when determining if it is violating bounds (default = True)

    Returns:
        A generator which returns all Var components block that are violating
        their bounds
    """
    for v in variables_generator(block, include_greybox=include_greybox):
        if v.value is not None:
            sf = get_scaling_factor(v, default=1, warning=False) if apply_scaling else 1
            tolerance = max(abs_tol, abs(value(v)) * rel_tol)
            if v.lb is not None and sf * value(v) <= sf * v.lb - tolerance:
                yield v
            elif v.ub is not None and sf * value(v) >= sf * v.ub + tolerance:
                yield v


def variables_violating_bounds_set(
    block: Block,
    abs_tol: float = 1e-4,
    rel_tol: float = 1e-4,
    include_greybox: bool = True,
    apply_scaling: bool = True,
):
    """
    Method to return a ComponentSet of all Var components in a model which have
    a value that violates their bounds.

    Args:
        block : model to be studied
        abs_tol : absolute tolerance for violation of bounds
        rel_tol : relative tolerance for violation of bounds
        include_greybox : whether to include greybox variables (default = True)
        apply_scaling: whether to apply scaling factors to the variable when determining if it is violating bounds (default = True)
    Returns:
        A ComponentSet including all Var components block that are violating
        their bounds
    """
    return ComponentSet(
        variables_violating_bounds_generator(
            block,
            abs_tol=abs_tol,
            rel_tol=rel_tol,
            include_greybox=include_greybox,
            apply_scaling=apply_scaling,
        )
    )


def number_variables_violating_bounds(
    block: Block,
    abs_tol: float = 1e-4,
    rel_tol: float = 1e-4,
    include_greybox: bool = True,
    apply_scaling: bool = True,
):
    """
    Method to return the number of all Var components in a model which have
    a value that violates their bounds.

    Args:
        block : model to be studied
        abs_tol : absolute tolerance for violation of bounds
        rel_tol : relative tolerance for violation of bounds
        include_greybox : whether to include greybox variables (default = True)
        apply_scaling: whether to apply scaling factors to the variable when determining if it is violating bounds (default = True)
    Returns:
        Number of components block that are violating their bounds
    """
    return len(
        variables_violating_bounds_set(
            block,
            abs_tol=abs_tol,
            rel_tol=rel_tol,
            include_greybox=include_greybox,
            apply_scaling=apply_scaling,
        )
    )


# -------------------------------------------------------------------------
# Variables in Constraints
def variables_in_activated_constraints_set(block, include_greybox=True):
    """
    Method to return a ComponentSet of all Var components which appear within a
    Constraint in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A ComponentSet including all Var components which appear within
        activated Constraints in block
    """
    var_set = ComponentSet()
    for c in _iter_indexed_block_data_objects(
        block, ctype=Constraint, active=True, descend_into=True
    ):
        for v in identify_variables(c.body):
            var_set.add(v)

    if include_greybox:
        for egb in _iter_indexed_block_data_objects(
            block, ctype=ExternalGreyBoxBlock, active=True, descend_into=True
        ):
            # For simplicity and efficiency, we will assume all variables in greybox blocks are included
            # For this to be False, the grey box itself would have to have a dangling variable
            for v in egb.component_data_objects(
                ctype=Var, active=None, descend_into=False
            ):
                var_set.add(v)
    return var_set


def number_variables_in_activated_constraints(block, include_greybox=True):
    """
    Method to return the number of Var components that appear within active
    Constraints in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        Number of Var components which appear within active Constraints in
        block
    """
    return len(
        variables_in_activated_constraints_set(block, include_greybox=include_greybox)
    )


def variables_not_in_activated_constraints_set(block, include_greybox=True):
    """
    Method to return a ComponentSet of all Var components which do not appear within a
    Constraint in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A ComponentSet including all Var components which do not appear within
        activated Constraints in block
    """
    all_vars = variables_set(block, include_greybox=include_greybox)
    active_vars = variables_in_activated_constraints_set(
        block, include_greybox=include_greybox
    )

    var_set = ComponentSet()
    for v in all_vars:
        if v not in active_vars:
            var_set.add(v)
    return var_set


def number_variables_not_in_activated_constraints(block, include_greybox=True):
    """
    Method to return the number of Var components that do not appear within active
    Constraints in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        Number of Var components which do not appear within active Constraints in
        block
    """
    return len(
        variables_not_in_activated_constraints_set(
            block, include_greybox=include_greybox
        )
    )


def variables_in_activated_equalities_set(block, include_greybox=True):
    """
    Method to return a ComponentSet of all Var components which appear within
    an equality Constraint in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A ComponentSet including all Var components which appear within
        activated equality Constraints in block
    """
    var_set = ComponentSet()
    # identify_variables does not work on ExternalGreyBoxConstraints, so we need to handle these separately
    for c in activated_equalities_generator(block, include_greybox=False):
        for v in identify_variables(c.body):
            var_set.add(v)

    if include_greybox:
        for egb in _iter_indexed_block_data_objects(
            block, ctype=ExternalGreyBoxBlock, active=True, descend_into=True
        ):
            # For simplicity and efficiency, we will assume all variables in greybox blocks are included
            # For this to be False, the grey box itself would have to have a dangling variable
            # All grey box constraints are equalities, so we do not need to check for inequality constraints here
            for v in egb.component_data_objects(
                ctype=Var, active=None, descend_into=False
            ):
                var_set.add(v)
    return var_set


def number_variables_in_activated_equalities(block, include_greybox=True):
    """
    Method to return the number of Var components which appear within activated
    equality Constraints in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

        block : model to be studied

    Returns:
        Number of Var components which appear within activated equality
        Constraints in block
    """
    return len(
        variables_in_activated_equalities_set(block, include_greybox=include_greybox)
    )


def variables_in_activated_inequalities_set(block):
    """
    Method to return a ComponentSet of all Var components which appear within
    an inequality Constraint in a model.

    Args:
        block : model to be studied

    Returns:
        A ComponentSet including all Var components which appear within
        activated inequality Constraints in block
    """
    # GreyBox constraints are all equalities, so we do not need to check for them here
    var_set = ComponentSet()
    for c in activated_inequalities_generator(block):
        for v in identify_variables(c.body):
            var_set.add(v)
    return var_set


def number_variables_in_activated_inequalities(block):
    """
    Method to return the number of Var components which appear within activated
    inequality Constraints in a model.

    Args:
        block : model to be studied

    Returns:
        Number of Var components which appear within activated inequality
        Constraints in block
    """
    return len(variables_in_activated_inequalities_set(block))


def variables_only_in_inequalities(block):
    """
    Method to return a ComponentSet of all Var components which appear only
    within inequality Constraints in a model.

    Args:
        block : model to be studied

    Returns:
        A ComponentSet including all Var components which appear only within
        inequality Constraints in block
    """
    return variables_in_activated_inequalities_set(
        block
    ) - variables_in_activated_equalities_set(block)


def number_variables_only_in_inequalities(block):
    """
    Method to return the number of Var components which appear only within
    activated inequality Constraints in a model.

    Args:
        block : model to be studied

    Returns:
        Number of Var components which appear only within activated inequality
        Constraints in block
    """
    return len(variables_only_in_inequalities(block))


# -------------------------------------------------------------------------
# Fixed Variables in Constraints
def fixed_variables_in_activated_equalities_set(block, include_greybox=True):
    """
    Method to return a ComponentSet of all fixed Var components which appear
    within an equality Constraint in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A ComponentSet including all fixed Var components which appear within
        activated equality Constraints in block
    """
    var_set = ComponentSet()
    for v in variables_in_activated_equalities_set(
        block, include_greybox=include_greybox
    ):
        if v.fixed:
            var_set.add(v)
    return var_set


def number_fixed_variables_in_activated_equalities(block, include_greybox=True):
    """
    Method to return the number of fixed Var components which appear within
    activated equality Constraints in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        Number of fixed Var components which appear within activated equality
        Constraints in block
    """
    return len(
        fixed_variables_in_activated_equalities_set(
            block, include_greybox=include_greybox
        )
    )


def unfixed_variables_in_activated_equalities_set(block, include_greybox=True):
    """
    Method to return a ComponentSet of all unfixed Var components which appear
    within an activated equality Constraint in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A ComponentSet of all unfixed Var components which appear within
        activated equality Constraints in block
    """
    var_set = ComponentSet()
    for v in variables_in_activated_equalities_set(
        block, include_greybox=include_greybox
    ):
        if not v.fixed:
            var_set.add(v)
    return var_set


def unfixed_greybox_variables(block):
    """
    Function to return a ComponentSet of all unfixed Var in GreyBoxModels

    Args:
        block : model to be studied

    Returns:
        A ComponentSet of all unfixed Var components which appear in Greybox models
    """
    # GreyBox variable should never be fixed - this function is useful to confirm that
    var_set = ComponentSet()
    for var in greybox_variables(block):
        if not var.fixed:
            var_set.add(var)
    return var_set


def greybox_variables(block):
    """
    Function to return a ComponentSet of all Var in GreyBoxModels

    Args:
        block : model to be studied

    Returns:
        A ComponentSet of all Var components which appear within
        activated Greybox model blocks
    """
    var_set = ComponentSet()
    for grey_box in activated_greybox_block_set(block):
        for v in grey_box.component_data_objects(
            ctype=Var, active=None, descend_into=False
        ):
            var_set.add(v)
    return var_set


def number_of_unfixed_greybox_variables(block):
    """
    Function to return a number of unfixed variables in grey box
    Args:
        block : model to be studied

    Returns:
        number of unfixed greybox variables
    """

    return len(unfixed_greybox_variables(block))


def number_of_greybox_variables(block):
    """
    Function to return a number of variables in grey box
    Args:
        block : model to be studied

    Returns:
        number of greybox variables
    """

    return len(greybox_variables(block))


def number_unfixed_variables_in_activated_equalities(block, include_greybox=True):
    """
    Method to return the number of unfixed Var components which appear within
    activated equality Constraints in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        Number of unfixed Var components which appear within activated equality
        Constraints in block
    """
    return len(
        unfixed_variables_in_activated_equalities_set(
            block, include_greybox=include_greybox
        )
    )


def fixed_variables_only_in_inequalities(block):
    """
    Method to return a ComponentSet of all fixed Var components which appear
    only within activated inequality Constraints in a model.

    Args:
        block : model to be studied

    Returns:
        A ComponentSet including all fixed Var components which appear only
        within activated inequality Constraints in block
    """
    var_set = ComponentSet()
    for v in variables_only_in_inequalities(block):
        if v.fixed:
            var_set.add(v)
    return var_set


def number_fixed_variables_only_in_inequalities(block):
    """
    Method to return the number of fixed Var components which only appear
    within activated inequality Constraints in a model.

    Args:
        block : model to be studied

    Returns:
        Number of fixed Var components which only appear within activated
        inequality Constraints in block
    """
    return len(fixed_variables_only_in_inequalities(block))


def external_variables_set(block: Block, include_greybox: bool = True) -> ComponentSet:
    """
    Method to return a ComponentSet of all Var components which appear within a
    Constraint in a model, but do not appear in the variables_set of the model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A ComponentSet including all Var components which appear within
        activated Constraints in block, but do not appear in the variables_set of
        the model
    """
    ext_vars = ComponentSet()
    local_var_set = variables_set(block, include_greybox=include_greybox, active=None)
    conc_var_set = variables_in_activated_constraints_set(
        block, include_greybox=include_greybox
    )
    for v in conc_var_set:
        if v not in local_var_set:
            ext_vars.add(v)
    return ext_vars


def number_external_variables(block: Block, include_greybox: bool = True) -> int:
    """
    Method to return the number of Var components which appear within a
    Constraint in a model, but do not appear in the variables_set of the model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        Number of Var components which appear within activated Constraints in block,
        but do not appear in the variables_set of the model
    """
    return len(external_variables_set(block, include_greybox=include_greybox))


def variables_fixed_to_zero_set(
    block: Block, include_greybox: bool = True
) -> ComponentSet:
    """
    Method to return a ComponentSet of all Var components which are fixed to zero in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A ComponentSet including all Var components which are fixed to zero in block
    """
    var_set = ComponentSet()
    for v in variables_set(block, include_greybox=include_greybox, active=None):
        if v.fixed and value(v) == 0:
            var_set.add(v)
    return var_set


def number_variables_fixed_to_zero(block: Block, include_greybox: bool = True) -> int:
    """
    Method to return the number of Var components which are fixed to zero in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        Number of Var components which are fixed to zero in block
    """
    return len(variables_fixed_to_zero_set(block, include_greybox=include_greybox))


def variables_near_zero_set(
    block: Block,
    tol: float = 1e-4,
    include_greybox: bool = True,
    scale_variables: bool = True,
) -> ComponentSet:
    """
    Method to return a ComponentSet of all Var components which have a value within tol of zero in a model.

    Args:
        block : model to be studied
        tol: Tolerance for inclusion in generator (default = 1e-4)
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)
        scale_variables: Boolean to scale variables by their scaling factor when determining if they are near zero (default = True)

    Returns:
        A ComponentSet including all Var components which have a value within tol of zero in block
    """
    var_set = ComponentSet()
    for v in variables_generator(block, include_greybox=include_greybox, active=None):
        sf = get_scaling_factor(v, default=1, warning=False) if scale_variables else 1
        if v.value is not None and abs(value(v) * sf) <= tol:
            var_set.add(v)
    return var_set


def number_variables_near_zero(
    block: Block,
    tol: float = 1e-4,
    include_greybox: bool = True,
    scale_variables: bool = True,
) -> int:
    """
    Method to return the number of Var components which have a value within tol of zero in a model.

    Args:
        block : model to be studied
        tol: Tolerance for inclusion in generator (default = 1e-4)
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)
        scale_variables: Boolean to scale variables by their scaling factor when determining if they are near zero (default = True)

    Returns:
        Number of Var components which have a value within tol of zero in block
    """
    return len(
        variables_near_zero_set(
            block,
            tol=tol,
            include_greybox=include_greybox,
            scale_variables=scale_variables,
        )
    )


def variables_with_none_value_set(block, include_greybox=True):
    """
    Method to return a ComponentSet of all Var components which have a value of None in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A ComponentSet including all Var components which have a value of None in block
    """
    var_set = ComponentSet()
    for v in variables_generator(block, include_greybox=include_greybox, active=None):
        if v.value is None:
            var_set.add(v)
    return var_set


def number_variables_with_none_value(block, include_greybox=True):
    """
    Method to return the number of Var components which have a value of None in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        Number of Var components which have a value of None in block
    """
    return len(variables_with_none_value_set(block, include_greybox=include_greybox))


def variables_with_extreme_values_set(
    block: Block,
    large: float,
    small: float,
    zero: float,
    include_greybox: bool = True,
    apply_scaling: bool = True,
):
    """Method to return a ComponentSet of all Var components which have a value with magnitude larger than large or smaller than small (but larger than zero) in a model.

    Args:
        block : model to be studied
        large: Threshold for large values
        small: Threshold for small values
        zero: Threshold for zero values (small values must be larger than this value to be included)
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)
        apply_scaling: Boolean to apply scaling factors to the variable when determining if it has an extreme value (default = True)
    Returns:
        A ComponentSet including all Var components which have a value with magnitude larger than large or smaller than small (but larger than zero) in block
    """
    extreme_vars = ComponentSet()
    for v in variables_generator(block, include_greybox=include_greybox, active=None):
        sf = get_scaling_factor(v, default=1, warning=False) if apply_scaling else 1
        if v.value is not None:
            mag = sf * abs(value(v))
            if mag > abs(large):
                extreme_vars.add(v)
            elif mag < abs(small) and mag > abs(zero):
                extreme_vars.add(v)

    return extreme_vars


def number_variables_with_extreme_values(
    block: Block,
    large: float,
    small: float,
    zero: float,
    include_greybox: bool = True,
    apply_scaling: bool = True,
):
    """Method to return the number of Var components which have a value with magnitude larger than large or smaller than small (but larger than zero) in a model.

    Args:
        block : model to be studied
        large: Threshold for large values
        small: Threshold for small values
        zero: Threshold for zero values (small values must be larger than this value to be included)
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)
        apply_scaling: Boolean to apply scaling factors to the variable when determining if it has an extreme value (default = True)
    Returns:
        Number of Var components which have a value with magnitude larger than large or smaller than small (but larger than zero) in block
    """
    return len(
        variables_with_extreme_values_set(
            block,
            large,
            small,
            zero,
            include_greybox=include_greybox,
            apply_scaling=apply_scaling,
        )
    )


# -------------------------------------------------------------------------
# Unused and un-Transformed Variables
def unused_variables_set(block, include_greybox=True):
    """
    Method to return a ComponentSet of all Var components which do not appear
    within any activated Constraint in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A ComponentSet including all Var components which do not appear within
        any Constraints in block
    """
    return variables_set(
        block, include_greybox=include_greybox
    ) - variables_in_activated_constraints_set(block, include_greybox=include_greybox)


def number_unused_variables(block, include_greybox=True):
    """
    Method to return the number of Var components which do not appear within
    any activated Constraint in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        Number of Var components which do not appear within any activated
        Constraints in block
    """
    return len(unused_variables_set(block, include_greybox=include_greybox))


def fixed_unused_variables_set(block, include_greybox=True):
    """
    Method to return a ComponentSet of all fixed Var components which do not
    appear within any activated Constraint in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A ComponentSet including all fixed Var components which do not appear
        within any Constraints in block
    """
    var_set = ComponentSet()
    for v in unused_variables_set(block, include_greybox=include_greybox):
        if v.fixed:
            var_set.add(v)
    return var_set


def number_fixed_unused_variables(block, include_greybox=True):
    """
    Method to return the number of fixed Var components which do not appear
    within any activated Constraint in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        Number of fixed Var components which do not appear within any activated
        Constraints in block
    """
    return len(fixed_unused_variables_set(block, include_greybox=include_greybox))


def derivative_variables_set(block):
    """
    Method to return a ComponentSet of all DerivativeVar components which
    appear in a model. Users should note that DerivativeVars are converted to
    ordinary Vars when a DAE transformation is applied. Thus, this method is
    useful for detecting any DerivativeVars which were do transformed.

    Args:
        block : model to be studied

    Returns:
        A ComponentSet including all DerivativeVar components which appear in
        block
    """
    return ComponentSet(
        _iter_indexed_block_data_objects(
            block, ctype=DerivativeVar, active=True, descend_into=True
        )
    )


def number_derivative_variables(block):
    """
    Method to return the number of DerivativeVar components which
    appear in a model. Users should note that DerivativeVars are converted to
    ordinary Vars when a DAE transformation is applied. Thus, this method is
    useful for detecting any DerivativeVars which were do transformed.

    Args:
        block : model to be studied

    Returns:
        Number of DerivativeVar components which appear in block
    """
    return len(derivative_variables_set(block))


# -------------------------------------------------------------------------
# Objective methods
def total_objectives_generator(block):
    """
    Generator which returns all Objective components in a model.

    Args:
        block : model to be studied

    Returns:
        A generator which returns all Objective components block
    """
    for o in activated_block_component_generator(block, ctype=Objective):
        yield o


def total_objectives_set(block):
    """
    Method to return a ComponentSet of all Objective components which appear
    in a model.

    Args:
        block : model to be studied

    Returns:
        A ComponentSet including all Objective components which appear in block
    """
    return ComponentSet(total_objectives_generator(block))


def number_total_objectives(block):
    """
    Method to return the number of Objective components which appear in a model

    Args:
        block : model to be studied

    Returns:
        Number of Objective components which appear in block
    """
    return sum(1 for _ in total_objectives_generator(block))


def activated_objectives_generator(block):
    """
    Generator which returns all activated Objective components in a model.

    Args:
        block : model to be studied

    Returns:
        A generator which returns all activated Objective components block
    """
    for o in activated_block_component_generator(block, ctype=Objective):
        if o.active:
            yield o


def activated_objectives_set(block):
    """
    Method to return a ComponentSet of all activated Objective components which
    appear in a model.

    Args:
        block : model to be studied

    Returns:
        A ComponentSet including all activated Objective components which
        appear in block
    """
    return ComponentSet(activated_objectives_generator(block))


def number_activated_objectives(block):
    """
    Method to return the number of activated Objective components which appear
    in a model.

    Args:
        block : model to be studied

    Returns:
        Number of activated Objective components which appear in block
    """
    return sum(1 for _ in activated_objectives_generator(block))


def deactivated_objectives_generator(block):
    """
    Generator which returns all deactivated Objective components in a model.

    Args:
        block : model to be studied

    Returns:
        A generator which returns all deactivated Objective components block
    """
    for o in activated_block_component_generator(block, ctype=Objective):
        if not o.active:
            yield o


def deactivated_objectives_set(block):
    """
    Method to return a ComponentSet of all deactivated Objective components
    which appear in a model.

    Args:
        block : model to be studied

    Returns:
        A ComponentSet including all deactivated Objective components which
        appear in block
    """
    return ComponentSet(deactivated_objectives_generator(block))


def number_deactivated_objectives(block):
    """
    Method to return the number of deactivated Objective components which
    appear in a model.

    Args:
        block : model to be studied

    Returns:
        Number of deactivated Objective components which appear in block
    """
    return sum(1 for _ in deactivated_objectives_generator(block))


# -------------------------------------------------------------------------
# Expression methods
# Always use ComponentsSets here to avoid duplication of References
def expressions_set(block):
    """
    Method to return a ComponentSet of all Expression components which appear
    in a model.

    Args:
        block : model to be studied

    Returns:
        A ComponentSet including all Expression components which  appear in
        block
    """
    return ComponentSet(
        _iter_indexed_block_data_objects(
            block, ctype=Expression, active=True, descend_into=True
        )
    )


def number_expressions(block):
    """
    Method to return the number of Expression components which appear in a
    model.

    Args:
        block : model to be studied

    Returns:
        Number of Expression components which  appear in block
    """
    return len(expressions_set(block))


# -------------------------------------------------------------------------
# Other model statistics
def degrees_of_freedom(block, include_greybox=True):
    """
    Method to return the degrees of freedom of a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include GreyBox models (default = True)

    Returns:
        Number of degrees of freedom in block.
    """
    return number_unfixed_variables_in_activated_equalities(
        block, include_greybox=include_greybox
    ) - number_activated_equalities(block, include_greybox=include_greybox)


def large_residuals_set(
    block,
    tol=1e-5,
    return_residual_values=False,
    include_greybox=True,
    apply_scaling=True,
):
    """
    Method to return a ComponentSet of all Constraint components with a
    residual greater than a given threshold which appear in a model.

    Args:
        block : model to be studied
        tol : residual threshold for inclusion in ComponentSet
        return_residual_values: boolean, if true return dictionary with
            residual values
        include_greybox: boolean, if true include GreyBox constraints (default = True)
        apply_scaling: boolean, if true apply scaling factors to the residuals when determining
        if they are large (default = True)

    Returns:
        large_residual_set: A ComponentSet including all Constraint components
        with a residual greater than tol which appear in block (if
        return_residual_values is false) residual_values: dictionary with
        constraint as key and residual (float) as value (if
        return_residual_values is true)
    """
    large_residuals_set = ComponentSet()
    residual_values = ComponentMap()

    for c in activated_constraints_generator(block, include_greybox=include_greybox):
        # Grey Box constraints cannot be scaled in the normal fashion
        sf = 1
        if apply_scaling and not isinstance(c, ExternalGreyBoxConstraint):
            sf = get_scaling_factor(c, default=1, warning=False)

        try:
            val = value(c.body)
        except ValueError:
            val = None
        if val is not None:
            if c.lb is None:
                r = 0
            else:
                r = max(c.lb - val, 0)

            if c.ub is not None:
                r = max(r, val - c.ub)

            if r * sf > tol:
                large_residuals_set.add(c)
                if return_residual_values:
                    residual_values[c] = r * sf
        else:
            large_residuals_set.add(c)

            if return_residual_values:
                residual_values[c] = None

    if return_residual_values:
        return residual_values
    else:
        return large_residuals_set


def number_large_residuals(block, tol=1e-5, include_greybox=True):
    """
    Method to return the number Constraint components with a residual greater
    than a given threshold which appear in a model.

    Args:
        block : model to be studied
        tol : residual threshold for inclusion in ComponentSet
        include_greybox: boolean, if true include GreyBox constraints (default = True)

    Returns:
        Number of Constraint components with a residual greater than tol which
        appear in block
    """
    return len(large_residuals_set(block, tol=tol, include_greybox=include_greybox))


def active_variables_in_deactivated_blocks_set(block, include_greybox=True):
    """
    Method to return a ComponentSet of any Var components which appear within
    an active Constraint but belong to a deactivated Block in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A ComponentSet including any Var components which belong to a
        deactivated Block but appear in an activate Constraint in block
    """
    var_set = ComponentSet()

    block_set = activated_blocks_set(block)
    if include_greybox:
        for c in activated_greybox_block_set(block):
            block_set.add(c)

    for v in variables_in_activated_constraints_set(
        block, include_greybox=include_greybox
    ):
        if v.parent_block() not in block_set:
            var_set.add(v)
    return var_set


def number_active_variables_in_deactivated_blocks(block, include_greybox=True):
    """
    Method to return the number of Var components which appear within an active
    Constraint but belong to a deactivated Block in a model.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        Number of Var components which belong to a deactivated Block but appear
        in an activate Constraint in block
    """
    return len(
        active_variables_in_deactivated_blocks_set(
            block, include_greybox=include_greybox
        )
    )


def variables_with_none_value_in_activated_equalities_set(block, include_greybox=True):
    """
    Method to return a ComponentSet of all Var components which
    have a value of None in the set of activated constraints.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        A ComponentSet including all Var components which
        have a value of None in the set of activated constraints.
    """
    var_set = ComponentSet()
    for v in variables_in_activated_equalities_set(
        block, include_greybox=include_greybox
    ):
        if v.value is None:
            var_set.add(v)
    return var_set


def number_variables_with_none_value_in_activated_equalities(
    block, include_greybox=True
):
    """
    Method to return the number of Var components which
    have a value of None in the set of activated constraints.

    Args:
        block : model to be studied
        include_greybox: Boolean to include implicit constraints from GreyBox
        models (default = True)

    Returns:
        Number of Var components which
        have a value of None in the set of activated constraints.
    """

    return len(
        variables_with_none_value_in_activated_equalities_set(
            block, include_greybox=include_greybox
        )
    )


# -------------------------------------------------------------------------
# Reporting methods
# TODO: This appears to be duplicated in the diagnostics module - we should consider whether to merge these or keep them separate
def report_statistics(block, ostream=None):
    """
    Method to print a report of the model statistics for a Pyomo Block

    Args:
        block : the Block object to report statistics from
        ostream : output stream for printing (defaults to sys.stdout)

    Returns:
        Printed output of the model statistics
    """
    if ostream is None:
        ostream = sys.stdout

    tab = " " * 4
    header = "=" * 72

    if block.name == "unknown":
        name_str = ""
    else:
        name_str = f"-  {block.name}"

    ostream.write("\n")
    ostream.write(header + "\n")
    ostream.write(f"Model Statistics  {name_str} \n")
    ostream.write("\n")
    ostream.write(f"Degrees of Freedom: " f"{degrees_of_freedom(block)} \n")
    ostream.write("\n")
    ostream.write(f"Total No. Variables: " f"{number_variables(block)} \n")
    ostream.write(
        f"{tab}No. Fixed Variables: " f"{number_fixed_variables(block)}" f"\n"
    )
    ostream.write(
        f"{tab}No. Unused Variables: "
        f"{number_unused_variables(block)} (Fixed):"
        f"{number_fixed_unused_variables(block)})"
        f"\n"
    )
    nv_alias = number_variables_only_in_inequalities
    nfv_alias = number_fixed_variables_only_in_inequalities
    ostream.write(
        f"{tab}No. Variables only in Inequalities:"
        f" {nv_alias(block)}"
        f" (Fixed: {nfv_alias(block)}) \n"
    )
    ostream.write("\n")
    ostream.write(f"Total No. Constraints: " f"{number_total_constraints(block)} \n")
    ostream.write(
        f"{tab}No. Equality Constraints: "
        f"{number_total_equalities(block)}"
        f" (Deactivated: "
        f"{number_deactivated_equalities(block)})"
        f"\n"
    )
    ostream.write(
        f"{tab}No. Inequality Constraints: "
        f"{number_total_inequalities(block)}"
        f" (Deactivated: "
        f"{number_deactivated_inequalities(block)})"
        f"\n"
    )
    ostream.write("\n")
    ostream.write(
        f"No. Objectives: "
        f"{number_total_objectives(block)}"
        f" (Deactivated: "
        f"{number_deactivated_objectives(block)})"
        f"\n"
    )
    ostream.write("\n")
    ostream.write(
        f"No. Blocks: {number_total_blocks(block)}"
        f" (Deactivated: "
        f"{number_deactivated_blocks(block)}) \n"
    )
    ostream.write(f"No. Expressions: " f"{number_expressions(block)} \n")
    if number_activated_greybox_blocks(block) != 0:
        ostream.write(
            f"No. Activated GreyBox Blocks: {number_activated_greybox_blocks(block)} \n"
        )
        ostream.write(f"No. GreyBox Variables: {number_of_greybox_variables(block)} \n")
        ostream.write(
            f"No. Fixed GreyBox Variables: {number_of_greybox_variables(block)-number_of_unfixed_greybox_variables(block)} \n"
        )
        ostream.write(
            f"No. GreyBox Equalities: {number_activated_greybox_equalities(block)} \n"
        )
    ostream.write(header + "\n")
    ostream.write("\n")


# -------------------------------------------------------------------------
# Common sub-methods
def activated_block_component_generator(block, ctype, include_greybox=False):
    """
    Generator which returns all the components of a given ctype which exist in
    activated Blocks within a model.

    Args:
        block : model to be studied
        ctype : type of Pyomo component to be returned by generator.
        include_greybox : whether to include greybox components

    Returns:
        A generator which returns all components of ctype which appear in
        activated Blocks in block
    """
    # Yield local components first
    for c in _iter_indexed_block_data_objects(
        block, ctype=ctype, active=None, descend_into=False
    ):
        yield c

    # Then yield components in active sub-blocks
    for b in _iter_indexed_block_data_objects(
        block, ctype=Block, active=True, descend_into=True
    ):
        for c in b.component_data_objects(ctype=ctype, active=None, descend_into=False):
            yield c

    # If include_greybox is True, also yield components in active greybox sub-blocks
    if include_greybox:
        for b in _iter_indexed_block_data_objects(
            block, ctype=ExternalGreyBoxBlock, active=True, descend_into=True
        ):
            for c in b.component_data_objects(
                ctype=ctype, active=None, descend_into=False
            ):
                yield c
