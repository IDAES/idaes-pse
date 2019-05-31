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

from pyomo.environ import Block, Constraint, Expression, Objective, Var, value
from pyomo.dae import DerivativeVar
from pyomo.core.expr.current import identify_variables
from pyomo.core.kernel.component_set import ComponentSet

from idaes.core.util.exceptions import ConfigurationError


def _create_component_set(block, ctype, active_blocks=True, descend_into=True):
    set_components = ComponentSet()

    if active_blocks is True and block.active is False:
        return set_components

    set_components.update(block.component_data_objects(ctype=ctype,
                                                       active=None,
                                                       descend_into=False))

    if descend_into:
        set_blocks = block_set(block, active=active_blocks, descend_into=True)

        for b in set_blocks:
            set_components.update(b.component_data_objects(ctype=ctype,
                                                           active=None,
                                                           descend_into=False))

    return set_components


def block_set(block, active=True, descend_into=True):
    """
    Method to return a set of the Block components appearing within a Pyomo
    Block.

    Args:
        block - the Block to be enumerated
        active - indicates whether components should be included based on their
                `active` status. Can have the follwoing values:
                * None - inlcude all components
                * True - include only components with active = True (default)
                * False - include only components with active = False
        descend_into - indicates whether only local Blocks should be enumerated
                (False), or if child Blocks should be included (True, default)

    Returns:
        set_blocks - a ComponentSet containing all the Blocks that
            appear within block
    """
    set_blocks = ComponentSet()

    set_blocks.update(block.component_data_objects(ctype=Block,
                                                   active=active,
                                                   descend_into=descend_into))

    return set_blocks


def variable_set(block, active=True, descend_into=True):
    """
    Method to return a set of the Variables appearing within a model.

    Args:
        block - the Block object to be enumerated
        active - indicates whether Vars in deactivated Blocks (and
                sub-Blocks) should be included. Valid option are:
                * None - include Vars in all Blocks
                * True - only include Vars in Blocks with active=True (default)
                * False - only include Vars in Blocks with active=False
        descend_into - indicates whether only local Vars within Blocks should
                be enumerated (False), or if Vars in child Blocks should be
                included (True, default)

    Returns:
        set_vars - a ComponentSet containing all Vars that appear within block
    """
    set_vars = ComponentSet()
    set_vars.update(block.component_data_objects(ctype=Var,
                                                 active=active,
                                                 descend_into=descend_into))
    return set_vars


def derivative_variables_set(block, active_blocks=True, descend_into=True):
    """
    Method to return a set of the DerivativeVars appearing within a model.
    Users should note that applying a DAE transformation converts
    DerivativeVars into ordinary Vars. Thus, this method is useful for
    identifying any DerivativeVars that have not been transformed.

    Args:
        block - the Block object to be enumerated
        active_blocks - indicates whether DerivativeVars in deactivated Blocks
                (and sub-Blocks) should be included. Valid option are:
                * None - include Vars in all Blocks
                * True - only include Vars in Blocks with active=True (default)
                * False - only include Vars in Blocks with active=False
        descend_into - indicates whether only local DerivativeVars within
                Blocks should be enumerated (False), or if DerivativeVars in
                child Blocks should be included (True, default)

    Returns:
        set_vars - a ComponentSet containing all DerivativeVars that appear
                within block
    """
    return _create_component_set(block,
                                 ctype=DerivativeVar,
                                 active_blocks=active_blocks,
                                 descend_into=descend_into)


def expression_set(block, active_blocks=True, descend_into=True):
    """
    Method to return a set of the Expression components appearing within a
    model.

    Args:
        block - the Block object to be enumerated
        active_blocks - indicates whether Expressions in deactivated Blocks
                (and sub-Blocks) should be included. Valid option are:
                * None - include Vars in all Blocks
                * True - only include Vars in Blocks with active=True (default)
                * False - only include Vars in Blocks with active=False
        descend_into - indicates whether only local Expressions within
                Blocks should be enumerated (False), or if Expressions in
                child Blocks should be included (True, default)

    Returns:
        set_expressions - a ComponentSet containing all Expressions that appear
            within block
    """
    return _create_component_set(block,
                                 ctype=Expression,
                                 active_blocks=active_blocks,
                                 descend_into=descend_into)


def objective_set(block, active_blocks=True, descend_into=True):
    """
    Method to return a set of the Objective components appearing within a
    model.

    Args:
        block - the Block object to be enumerated
        active_blocks - indicates whether Objectives in deactivated Blocks
                (and sub-Blocks) should be included. Valid option are:
                * None - include Vars in all Blocks
                * True - only include Vars in Blocks with active=True (default)
                * False - only include Vars in Blocks with active=False
        descend_into - indicates whether only local Objectives within
                Blocks should be enumerated (False), or if Objectives in
                child Blocks should be included (True, default)

    Returns:
        set_expressions - a ComponentSet containing all Objectives that appear
            within block
    """
    return _create_component_set(block,
                                 ctype=Objective,
                                 active_blocks=active_blocks,
                                 descend_into=descend_into)


def constraint_set(block, active_blocks=True, descend_into=True):
    """
    Method to return a set of the Constraint components appearing within a
    model.

    Args:
        block - the Block object to be enumerated
        active_blocks - indicates whether Constraints in deactivated Blocks
                (and sub-Blocks) should be included. Valid option are:
                * None - include Vars in all Blocks
                * True - only include Vars in Blocks with active=True (default)
                * False - only include Vars in Blocks with active=False
        descend_into - indicates whether only local Constraints within
                Blocks should be enumerated (False), or if Constraints in
                child Blocks should be included (True, default)

    Returns:
        set_expressions - a ComponentSet containing all Constraints that appear
            within block
    """
    return _create_component_set(block,
                                 ctype=Constraint,
                                 active_blocks=active_blocks,
                                 descend_into=descend_into)


def equality_constraint_set(block_or_set,
                            active_blocks=True,
                            descend_into=True):
    """
    Method to return a set of the equality Constraint components appearing
    within a model.

    Args:
        block_or_set - the Block object to be enumerated or an existing
                ComponentSet of Constraints to be enumerated
        active_blocks - used when block_or_set is a Block object and indicates
                whether Constraints in deactivated Blocks (and sub-Blocks)
                should be included. Valid option are:
                * None - include Vars in all Blocks
                * True - only include Vars in Blocks with active=True (default)
                * False - only include Vars in Blocks with active=False
        descend_into - used when block_or_set is a Block object and indicates
                whether only local Constraints within Blocks should be
                enumerated (False), or if Constraints in child Blocks should be
                included (True, default)

    Returns:
        set_equalites - a ComponentSet containing all equality Constraints that
            appear within block_or_set
    """
    if not isinstance(block_or_set, ComponentSet):
        block_or_set = constraint_set(
                            block_or_set,
                            active_blocks=active_blocks,
                            descend_into=descend_into)

    set_equalities = ComponentSet()

    for c in block_or_set:
        if c.upper is not None and c.lower is not None and c.upper == c.lower:
            # Constraint is an equality constraint
            set_equalities.add(c)

    return set_equalities


def inequality_constraints_set(block_or_set,
                               active_blocks=True,
                               descend_into=True):
    """
    Method to return a set of the inequality Constraint components appearing
    within a model.

    Args:
        block_or_set - the Block object to be enumerated or an existing
                ComponentSet of Constraints to be enumerated
        active_blocks - used when block_or_set is a Block object and indicates
                whether Constraints in deactivated Blocks (and sub-Blocks)
                should be included. Valid option are:
                * None - include Vars in all Blocks
                * True - only include Vars in Blocks with active=True (default)
                * False - only include Vars in Blocks with active=False
        descend_into - used when block_or_set is a Block object and indicates
                whether only local Constraints within Blocks should be
                enumerated (False), or if Constraints in child Blocks should be
                included (True, default)

    Returns:
        set_inequalites - a ComponentSet containing all inequality Constraints
            that appear within block_or_set
    """
    if not isinstance(block_or_set, ComponentSet):
        block_or_set = constraint_set(
                            block_or_set,
                            active_blocks=active_blocks,
                            descend_into=descend_into)

    set_inequalities = ComponentSet()

    for c in block_or_set:
        if c.upper is None or c.lower is None:
            # Constraint is an inequality constraint
            set_inequalities.add(c)

    return set_inequalities


def large_residual_set(block_or_set,
                       tol=1e-5,
                       active_blocks=True,
                       descend_into=True):
    """
    Method to return a set of the Constraint components with large residuals
    appearing within a model.

    Args:
        block_or_set - the Block object to be enumerated or an existing
                ComponentSet of Constraints to be enumerated
        tol - show Constraints with residuals greated than tol
        active_blocks - used when block_or_set is a Block object and indicates
                whether Constraints in deactivated Blocks (and sub-Blocks)
                should be included. Valid option are:
                * None - include Vars in all Blocks
                * True - only include Vars in Blocks with active=True (default)
                * False - only include Vars in Blocks with active=False
        descend_into - used when block_or_set is a Block object and indicates
                whether only local Constraints within Blocks should be
                enumerated (False), or if Constraints in child Blocks should be
                included (True, default)

    Returns:
        set_equalites - a ComponentSet containing all equality Constraints that
            appear within block_or_set
    """
    if not isinstance(block_or_set, ComponentSet):
        block_or_set = constraint_set(
                            block_or_set,
                            active_blocks=active_blocks,
                            descend_into=descend_into)

    large_residuals = ComponentSet()

    for c in block_or_set:
        if c.active and value(c.lower - c.body()) > tol:
            large_residuals.add(c)
        elif c.active and value(c.body() - c.upper) > tol:
            large_residuals.add(c)

    return large_residuals


def activated_component_set(component_set):
    """
    Method to return a set of the activated components appearing within a
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


def variables_in_constraints_set(block_or_set,
                                 active_blocks=True,
                                 descend_into=True):
    """
    Method to return a set of the Vars which occur within Constraints within a
    model.

    Args:
        block_or_set - the Block object containing the Constraints to be
                enumerated or an existing ComponentSet of Constraints
        active_blocks - used when block_or_set is a Block object and indicates
                whether Constraints in deactivated Blocks (and sub-Blocks)
                should be included. Valid option are:
                * None - include Vars in all Blocks
                * True - only include Vars in Blocks with active=True (default)
                * False - only include Vars in Blocks with active=False
        descend_into - used when block_or_set is a Block object and indicates
                whether only local Constraints within Blocks should be
                enumerated (False), or if Constraints in child Blocks should be
                included (True, default)

    Returns:
        set_vars_in_constraints - a ComponentSet containing all Var components
            which appear in the Constraints in block_or_set
    """
    if not isinstance(block_or_set, ComponentSet):
        block_or_set = constraint_set(
                            block_or_set,
                            active_blocks=active_blocks,
                            descend_into=descend_into)

    set_vars_in_constraints = ComponentSet()

    for c in block_or_set:
        for v in identify_variables(c.body):
            set_vars_in_constraints.add(v)

    return set_vars_in_constraints


def fixed_variable_set(block_or_set,
                       active=True,
                       descend_into=True):
    """
    Method to enumerate the fixed Variables which within a model.

    Args:
        block_or_set - the Block object containing the Vars to be
                enumerated or an existing ComponentSet of Vars
        active - indicates whether Vars in deactivated Blocks (and
                sub-Blocks) should be included. Valid option are:
                * None - include Vars in all Blocks
                * True - only include Vars in Blocks with active=True (default)
                * False - only include Vars in Blocks with active=False
        descend_into - used when block_or_set is a Block object and indicates
                whether only local Vars within Blocks should be
                enumerated (False), or if Vars in child Blocks should be
                included (True, default)

    Returns:
        set_fixed_vars - a ComponentSet containing all fixed Var components
            which appear in block_or_set
    """
    if not isinstance(block_or_set, ComponentSet):
        block_or_set = variable_set(
                            block_or_set,
                            active=active,
                            descend_into=descend_into)

    set_fixed_vars = ComponentSet()

    for v in block_or_set:
        if v.fixed:
            set_fixed_vars.add(v)

    return set_fixed_vars


def unfixed_variable_set(block_or_set,
                         active=True,
                         descend_into=True):
    """
    Method to enumerate the unfixed Variables which within a model.

    Args:
        block_or_set - the Block object containing the Vars to be
                enumerated or an existing ComponentSet of Vars
        active - indicates whether Vars in deactivated Blocks (and
                sub-Blocks) should be included. Valid option are:
                * None - include Vars in all Blocks
                * True - only include Vars in Blocks with active=True (default)
                * False - only include Vars in Blocks with active=False
        descend_into - used when block_or_set is a Block object and indicates
                whether only local Vars within Blocks should be
                enumerated (False), or if Vars in child Blocks should be
                included (True, default)

    Returns:
        set_fixed_vars - a ComponentSet containing all unfixed Var components
            which appear in block_or_set
    """
    if not isinstance(block_or_set, ComponentSet):
        block_or_set = variable_set(
                            block_or_set,
                            active=active,
                            descend_into=descend_into)

    set_fixed_vars = ComponentSet()

    for v in block_or_set:
        if not v.fixed:
            set_fixed_vars.add(v)

    return set_fixed_vars


def active_variables_in_deactived_blocks_set(block):
    """
    This method returns a set of any Vars which belong to a deactivated Block
    (or a child of a deactivated Block) but appear in an active Constraint
    within a model. This can cause problems with the Pyomo solver writers, and
    is often an indication of a mistake when setting up the model.

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
    all_constraints = constraint_set(block,
                                     active_blocks=True,
                                     descend_into=True)

    # Get activated Constraints from all_constraints
    act_constraints = activated_component_set(all_constraints)

    # Get all Vars that appear in activated Constraints
    vars_in_act_constraints = variables_in_constraints_set(
            act_constraints)

    # Get all active Blocks
    act_blocks = block_set(block, active=True, descend_into=True)

    act_vars_in_deact_blocks = ComponentSet()

    for v in vars_in_act_constraints:
        if (v.parent_block() not in act_blocks and
                v.parent_block() is not block):
            act_vars_in_deact_blocks.add(v)

    return act_vars_in_deact_blocks


def calculate_degrees_of_freedom(block):
    equalities = equality_constraint_set(block,
                                         active_blocks=True,
                                         descend_into=True)
    act_equalities = activated_component_set(equalities)

    vars_in_act_equals = variables_in_constraints_set(act_equalities)

    fixed_vars_in_act_equals = fixed_variable_set(vars_in_act_equals)

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
    if deactivated_blocks:
        active = None
    else:
        active = True

    total_vars = variable_set(block, active=active, descend_into=descend_into)
    fixed_vars = fixed_variable_set(total_vars)

    total_cons = constraint_set(block,
                                active_blocks=active,
                                descend_into=descend_into)
    eq_cons = equality_constraint_set(total_cons)
    act_eq_cons = activated_component_set(eq_cons)
    ineq_cons = inequality_constraints_set(total_cons)
    act_ineq_cons = activated_component_set(ineq_cons)

    total_objs = objective_set(block,
                               active_blocks=active,
                               descend_into=descend_into)
    act_objs = activated_component_set(total_objs)

    total_exprs = expression_set(block,
                                 active_blocks=active,
                                 descend_into=descend_into)

    total_blocks = block_set(block, active=active, descend_into=descend_into)
    act_blocks = activated_component_set(total_blocks)

    vars_in_act_equals = variables_in_constraints_set(act_eq_cons)
    fixed_vars_in_act_equals = fixed_variable_set(vars_in_act_equals)

    vars_in_inequals = variables_in_constraints_set(act_ineq_cons)
    vars_only_in_inequals = ComponentSet()
    for v in vars_in_inequals:
        if v not in vars_in_act_equals:
            vars_only_in_inequals.add(v)
    fixed_vars_only_in_inequals = \
        fixed_variable_set(vars_only_in_inequals)

    unused_vars = total_vars-vars_in_act_equals
    fixed_unused_vars = fixed_variable_set(unused_vars)

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
