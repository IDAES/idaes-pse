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
This module contains utility functions for reporting structural statistics of
IDAES models.
"""

__author__ = "Andrew Lee"

from pyomo.environ import Block, Constraint, Expression, Objective, Var
from pyomo.dae import DerivativeVar
from pyomo.core.expr.current import identify_variables
from pyomo.core.kernel.component_set import ComponentSet

from idaes.core.util.exceptions import ConfigurationError


# TODO : Need to deal with variables in inactive constraints

def enumerate_component(block, ctype,
                        deactivated_blocks=True, descend_into=True):

    set_components = ComponentSet()

    if deactivated_blocks is False and block.active is False:
        return set_components

    # Enumerate local components
    for o in block.component_data_objects(ctype=ctype, descend_into=False):
        set_components.add(o)

    if descend_into:
        for b in block.component_data_objects(ctype=Block, descend_into=False):
            sc = enumerate_component(b,
                                     ctype,
                                     deactivated_blocks=deactivated_blocks,
                                     descend_into=True)

            set_components.update(sc)

    return set_components


def enumerate_blocks(block, deactivated_blocks=False, descend_into=True):
    """
    Method to enumerate the Block components appearing within a Pyomo Block.

    Args:
        block - the Block to be enumerated
        deactivated_blocks - indicates whether deactivated Blocks (and
                sub-Blocks) should be included (deafult = False)
        descend_into - indicates whether only local Blocks should be enumerated
                (False), or if child Blocks should be included (True, default)

    Returns:
        set_blocks - a ComponentSet containing all the Blocks that
            appear within block
    """
    set_blocks = ComponentSet()

    if deactivated_blocks is False and block.active is False:
        return set_blocks

    # Enumerate local components
    for o in block.component_data_objects(ctype=Block, descend_into=False):
        if deactivated_blocks is True or o.active is True:
            set_blocks.add(o)

    if descend_into:
        for b in block.component_data_objects(ctype=Block, descend_into=False):
            sb = enumerate_blocks(b,
                                  deactivated_blocks=deactivated_blocks,
                                  descend_into=True)

            set_blocks.update(sb)

    return set_blocks


def enumerate_variables(block, deactivated_blocks=False, descend_into=True):
    """
    Method to enumerate the Variables appearing within a model.

    Args:
        block - the Block object to be enumerated
        deactivated_blocks - indicates whether Vars in deactivated Blocks (and
                sub-Blocks) should be included (deafult = False)
        descend_into - indicates whether only local Vars within Blocks should
                be enumerated (False), or if Vars in child Blocks should be
                included (True, default)

    Returns:
        set_vars - a ComponentSet containing all Vars that appear within block
    """
    return enumerate_component(block,
                               ctype=Var,
                               deactivated_blocks=deactivated_blocks,
                               descend_into=descend_into)


def enumerate_derivative_variables(block,
                                   deactivated_blocks=False,
                                   descend_into=True):
    """
    Method to enumerate the DerivativeVars appearing within a model. Users
    should note that applying a DAE transformation converts DerivativeVars into
    ordinary Vars. Thus, this method is useful for identifying any
    DerivativeVars that have not been trasnformed.

    Args:
        block - the Block object to be enumerated
        deactivated_blocks - indicates whether DerivativeVars in deactivated
                Blocks (and sub-Blocks) should be included (deafult = False)
        descend_into - indicates whether only local DerivativeVars within
                Blocks should be enumerated (False), or if DerivativeVars in
                child Blocks should be included (True, default)

    Returns:
        set_vars - a ComponentSet containing all DerivativeVars that appear
                within block
    """
    return enumerate_component(block,
                               ctype=DerivativeVar,
                               deactivated_blocks=deactivated_blocks,
                               descend_into=descend_into)


def enumerate_expressions(block, deactivated_blocks=False, descend_into=True):
    """
    Method to enumerate the Expression components appearing within a model.

    Args:
        block - the Block object to be enumerated
        deactivated_blocks - indicates whether Expressions in deactivated
                Blocks (and sub-Blocks) should be included (deafult = False)
        descend_into - indicates whether only local Expressions within
                Blocks should be enumerated (False), or if Expressions in
                child Blocks should be included (True, default)

    Returns:
        set_expressions - a ComponentSet containing all Expressions that appear
            within block
    """
    return enumerate_component(block,
                               ctype=Expression,
                               deactivated_blocks=deactivated_blocks,
                               descend_into=descend_into)


def enumerate_objectives(block, deactivated_blocks=False, descend_into=True):
    """
    Method to enumerate the Objective components appearing within a model.

    Args:
        block - the Block object to be enumerated
        deactivated_blocks - indicates whether Objective in deactivated
                Blocks (and sub-Blocks) should be included (deafult = False)
        descend_into - indicates whether only local Objectives within
                Blocks should be enumerated (False), or if Objectives in
                child Blocks should be included (True, default)

    Returns:
        set_expressions - a ComponentSet containing all Objectives that appear
            within block
    """
    return enumerate_component(block,
                               ctype=Objective,
                               deactivated_blocks=deactivated_blocks,
                               descend_into=descend_into)


def enumerate_constraints(block, deactivated_blocks=False, descend_into=True):
    """
    Method to enumerate the Constraint components appearing within a model.

    Args:
        block - the Block object to be enumerated
        deactivated_blocks - indicates whether Constraints in deactivated
                Blocks (and sub-Blocks) should be included (deafult = False)
        descend_into - indicates whether only local Constraints within
                Blocks should be enumerated (False), or if Constraints in
                child Blocks should be included (True, default)

    Returns:
        set_expressions - a ComponentSet containing all Constraints that appear
            within block
    """
    return enumerate_component(block,
                               ctype=Constraint,
                               deactivated_blocks=deactivated_blocks,
                               descend_into=descend_into)


def enumerate_equality_constraints(block_or_set,
                                   deactivated_blocks=False,
                                   descend_into=True):
    """
    Method to enumerate the equality Constraint components appearing within a
    model.

    Args:
        block_or_set - the Block object to be enumerated or an existing
                ComponentSet of Constraints to be enumerated
        deactivated_blocks - used when block_or_set is a Block object and
                indicates whether Constraints in deactivated Blocks (and
                sub-Blocks) should be included (deafult = False)
        descend_into - used when block_or_set is a Block object and indicates
                whether only local Constraints within Blocks should be
                enumerated (False), or if Constraints in child Blocks should be
                included (True, default)

    Returns:
        set_equalites - a ComponentSet containing all equality Constraints that
            appear within block_or_set
    """
    if not isinstance(block_or_set, ComponentSet):
        block_or_set = enumerate_constraints(
                            block_or_set,
                            deactivated_blocks=deactivated_blocks,
                            descend_into=descend_into)

    set_equalities = ComponentSet()

    for c in block_or_set:
        if c.upper is not None and c.lower is not None and c.upper == c.lower:
            # Constraint is an equality constraint
            set_equalities.add(c)

    return set_equalities


def enumerate_inequality_constraints(block_or_set,
                                     deactivated_blocks=False,
                                     descend_into=True):
    """
    Method to enumerate the inequality Constraint components appearing within a
    model.

    Args:
        block_or_set - the Block object to be enumerated or an existing
                ComponentSet of Constraints to be enumerated
        deactivated_blocks - used when block_or_set is a Block object and
                indicates whether Constraints in deactivated Blocks (and
                sub-Blocks) should be included (deafult = False)
        descend_into - used when block_or_set is a Block object and indicates
                whether only local Constraints within Blocks should be
                enumerated (False), or if Constraints in child Blocks should be
                included (True, default)

    Returns:
        set_inequalites - a ComponentSet containing all inequality Constraints
            that appear within block_or_set
    """
    if not isinstance(block_or_set, ComponentSet):
        block_or_set = enumerate_constraints(
                            block_or_set,
                            deactivated_blocks=deactivated_blocks,
                            descend_into=descend_into)

    set_inequalities = ComponentSet()

    for c in block_or_set:
        if c.upper is None or c.lower is None:
            # Constraint is an inequality constraint
            set_inequalities.add(c)

    return set_inequalities


def enumerate_activated_components(component_set):
    """
    Method to enumerate the activated components appearing within a
    set of components.

    Args:
        block_or_set - a ComponentSet of components to be enumerated

    Returns:
        set_activated - a ComponentSet containing all activated components
            that appear within component_set
    """
    set_activated = ComponentSet()

    for c in component_set:
        if c.active:
            set_activated.add(c)

    return set_activated


def enumerate_variables_in_constraints(block_or_set,
                                       deactivated_blocks=False,
                                       descend_into=True):
    """
    Method to enumerate the Variables which occur within Constraints within a
    model.

    Args:
        block_or_set - the Block object containing the Constraints to be
                enumerated or an existing ComponentSet of Constraints
        deactivated_blocks - used when block_or_set is a Block object and
                indicates whether Constraints in deactivated Blocks (and
                sub-Blocks) should be included (deafult = False)
        descend_into - used when block_or_set is a Block object and indicates
                whether only local Constraints within Blocks should be
                enumerated (False), or if Constraints in child Blocks should be
                included (True, default)

    Returns:
        set_vars_in_constraints - a ComponentSet containing all Var components
            which appear in the Constraints in block_or_set
    """
    if not isinstance(block_or_set, ComponentSet):
        block_or_set = enumerate_constraints(
                            block_or_set,
                            deactivated_blocks=deactivated_blocks,
                            descend_into=descend_into)

    set_vars_in_constraints = ComponentSet()

    for c in block_or_set:
        for v in identify_variables(c.body):
            set_vars_in_constraints.add(v)

    return set_vars_in_constraints


def enumerate_fixed_variables(block_or_set,
                              deactivated_blocks=False,
                              descend_into=True):
    """
    Method to enumerate the fixed Variables which within a model.

    Args:
        block_or_set - the Block object containing the Vars to be
                enumerated or an existing ComponentSet of Vars
        deactivated_blocks - used when block_or_set is a Block object and
                indicates whether Vars in deactivated Blocks (and
                sub-Blocks) should be included (deafult = False)
        descend_into - used when block_or_set is a Block object and indicates
                whether only local Vars within Blocks should be
                enumerated (False), or if Vars in child Blocks should be
                included (True, default)

    Returns:
        set_fixed_vars - a ComponentSet containing all fixed Var components
            which appear in block_or_set
    """
    if not isinstance(block_or_set, ComponentSet):
        block_or_set = enumerate_variables(
                            block_or_set,
                            deactivated_blocks=deactivated_blocks,
                            descend_into=descend_into)

    set_fixed_vars = ComponentSet()

    for v in block_or_set:
        if v.fixed:
            set_fixed_vars.add(v)

    return set_fixed_vars


def enumerate_active_varaibles_in_deactived_blocks(block):
    """
    This method enumerates any Vars which belong to a deactivated Block (or
    a child of a deactiveted Block) but appear in an active Constraint within
    a model. This can cause problems with the Pyomo solver writers, and is
    often an indication of a mistake when setting up the model.

    Args:
        block - the Block object to be enumerated

    Returns:
        act_vars_in_deact_blocks - a ComponentSet containing all Vars from
            deactived Blocks which appear in active Constraints within block.
    """
    # Check that block is activated
    if block.active is False:
        raise ConfigurationError(
                f"Cannot enumerate active variables in deactivated "
                "constraints. The Block provided ({block.name}) is deactivated"
                )

    # Get all Constraints in active Blocks
    all_constraints = enumerate_constraints(block,
                                            deactivated_blocks=False,
                                            descend_into=True)

    # Get activated Constraints from all_constraints
    act_constraints = enumerate_activated_components(all_constraints)

    # Get all Vars that appear in activated Constraints
    vars_in_act_constraints = enumerate_variables_in_constraints(
            act_constraints)

    # Get all active Blocks
    act_blocks = enumerate_blocks(block,
                                  deactivated_blocks=False,
                                  descend_into=True)

    act_vars_in_deact_blocks = ComponentSet()

    for v in vars_in_act_constraints:
        if (v.parent_block() not in act_blocks and
                v.parent_block() is not block):
            act_vars_in_deact_blocks.add(v)

    return act_vars_in_deact_blocks


def calculate_degrees_of_freedom(block):
    equalities = enumerate_equality_constraints(block,
                                                deactivated_blocks=False,
                                                descend_into=True)
    act_equalities = enumerate_activated_components(equalities)

    vars_in_act_equals = enumerate_variables_in_constraints(act_equalities)

    fixed_vars_in_act_equals = enumerate_fixed_variables(vars_in_act_equals)

    dof = len(vars_in_act_equals-fixed_vars_in_act_equals)-len(act_equalities)

    return dof


def report_model_statistics(block,
                            deactivated_blocks=False,
                            descend_into=True):
    """
    Method to print a report of the model statistics for a Pyomo Block

    Args:
        block - the Block object to report statistics from
        deactivated_blocks - indicates whether to include deactivated Blocks
                (and sub-Blocks) in report (deafult = False)
        descend_into - indicates whether only components within block should be
                included (False), or if components in child Blocks should be
                included (True, default)

    Returns:
        Printed output of the model statistics
    """
    total_vars = enumerate_variables(block,
                                     deactivated_blocks=deactivated_blocks,
                                     descend_into=descend_into)
    fixed_vars = enumerate_fixed_variables(total_vars)

    total_cons = enumerate_constraints(block,
                                       deactivated_blocks=deactivated_blocks,
                                       descend_into=descend_into)
    eq_cons = enumerate_equality_constraints(total_cons)
    act_eq_cons = enumerate_activated_components(eq_cons)
    ineq_cons = enumerate_inequality_constraints(total_cons)
    act_ineq_cons = enumerate_activated_components(ineq_cons)

    total_objs = enumerate_objectives(block,
                                      deactivated_blocks=deactivated_blocks,
                                      descend_into=descend_into)
    act_objs = enumerate_activated_components(total_objs)

    total_exprs = enumerate_expressions(block,
                                        deactivated_blocks=deactivated_blocks,
                                        descend_into=descend_into)

    total_blocks = enumerate_blocks(block,
                                    deactivated_blocks=deactivated_blocks,
                                    descend_into=descend_into)
    act_blocks = enumerate_activated_components(total_blocks)

    vars_in_act_equals = enumerate_variables_in_constraints(act_eq_cons)
    fixed_vars_in_act_equals = enumerate_fixed_variables(vars_in_act_equals)

    vars_in_inequals = enumerate_variables_in_constraints(act_ineq_cons)
    vars_only_in_inequals = ComponentSet()
    for v in vars_in_inequals:
        if v not in vars_in_act_equals:
            vars_only_in_inequals.add(v)
    fixed_vars_only_in_inequals = \
        enumerate_fixed_variables(vars_only_in_inequals)

    unused_vars = total_vars-vars_in_act_equals
    fixed_unused_vars = enumerate_fixed_variables(unused_vars)

    dof = len(vars_in_act_equals-fixed_vars_in_act_equals)-len(act_eq_cons)

    pad = " "*4
    header = '='*64

    if block.name == "unknown":
        name_str = ""
    else:
        name_str = f"-  ({block.name})"

    print()
    print(header)
    print(f"Model Statistics  {name_str}")
    print()
    print(f"Degrees of Freedom: {dof}")
    print()
    print(f"Total No. Variables: {len(total_vars)}")
    print(f"{pad}No. Fixed Variables: {len(fixed_vars)}")
    print(f"{pad}No. Unused Variables: {len(unused_vars)}"
          f" (Fixed: {len(fixed_unused_vars)})")
    print(f"{pad}No. Variables only in Inequalities:"
          f" {len(vars_only_in_inequals)}"
          f" (Fixed: {len(fixed_vars_only_in_inequals)})")
    print()
    print(f"Total No. Constraints: {len(total_cons)}")
    print(f"{pad}No. Equality Constraints: {len(eq_cons)}"
          f" (Deactivated: {len(eq_cons-act_eq_cons)})")
    print(f"{pad}No. Inequality Constraints: {len(ineq_cons)}"
          f" (Deactivated: {len(ineq_cons-act_ineq_cons)})")
    print()
    print(f"No. Objectives: {len(total_objs)}"
          f" (Deactivated: {len(total_objs-act_objs)})")
    print()
    print(f"No. Blocks: {len(total_blocks)}"
          f" (Deactivated: {len(total_blocks-act_blocks)})")
    print(f"No. Expressions: {len(total_exprs)}")
    print(header)
    print()
