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
This module contains utility functions for reporting model diagnostics.
"""

__author__ = "Alexander Dowling, Douglas Allan, Andrew Lee, Robby Parker, Ben Knueven"

from pyomo.common.collections import ComponentSet

from idaes.core.util.model_statistics import (
    activated_blocks_set,
    deactivated_blocks_set,
    activated_equalities_set,
    deactivated_equalities_set,
    activated_inequalities_set,
    deactivated_inequalities_set,
    activated_objectives_set,
    deactivated_objectives_set,
    variables_in_activated_constraints_set,
    number_activated_greybox_equalities,
    number_deactivated_greybox_equalities,
    activated_greybox_block_set,
    deactivated_greybox_block_set,
    greybox_block_set,
    unfixed_greybox_variables,
    greybox_variables,
)
from idaes.core.util.diagnostics_tools.utils import (
    var_in_block,
)

MAX_STR_LENGTH = 84
TAB = " " * 4


def collect_model_statistics(model):
    """
    Collects statistics about the model that are relevant for diagnostics tools.
    """
    vars_in_constraints = variables_in_activated_constraints_set(model)
    fixed_vars_in_constraints = ComponentSet()
    free_vars_in_constraints = ComponentSet()
    free_vars_lb = ComponentSet()
    free_vars_ub = ComponentSet()
    free_vars_lbub = ComponentSet()
    ext_fixed_vars_in_constraints = ComponentSet()
    ext_free_vars_in_constraints = ComponentSet()
    for v in vars_in_constraints:
        if v.fixed:
            fixed_vars_in_constraints.add(v)
            if not var_in_block(v, model):
                ext_fixed_vars_in_constraints.add(v)
        else:
            free_vars_in_constraints.add(v)
            if not var_in_block(v, model):
                ext_free_vars_in_constraints.add(v)
            if v.lb is not None:
                if v.ub is not None:
                    free_vars_lbub.add(v)
                else:
                    free_vars_lb.add(v)
            elif v.ub is not None:
                free_vars_ub.add(v)

    # Generate report
    # TODO: Binary and boolean vars
    stats = []
    stats.append(
        f"{TAB}Activated Blocks: {len(activated_blocks_set(model))} "
        f"(Deactivated: {len(deactivated_blocks_set(model))})"
    )
    stats.append(
        f"{TAB}Free Variables in Activated Constraints: "
        f"{len(free_vars_in_constraints)} "
        f"(External: {len(ext_free_vars_in_constraints)})"
    )
    stats.append(f"{TAB * 2}Free Variables with only lower bounds: {len(free_vars_lb)}")
    stats.append(f"{TAB * 2}Free Variables with only upper bounds: {len(free_vars_ub)}")
    stats.append(
        f"{TAB * 2}Free Variables with upper and lower bounds: "
        f"{len(free_vars_lbub)}"
    )
    stats.append(
        f"{TAB}Fixed Variables in Activated Constraints: "
        f"{len(fixed_vars_in_constraints)} "
        f"(External: {len(ext_fixed_vars_in_constraints)})"
    )
    stats.append(
        f"{TAB}Activated Equality Constraints: {len(activated_equalities_set(model))+number_activated_greybox_equalities(model)} "
        f"(Deactivated: {len(deactivated_equalities_set(model))+number_deactivated_greybox_equalities(model)})"
    )
    stats.append(
        f"{TAB}Activated Inequality Constraints: {len(activated_inequalities_set(model))} "
        f"(Deactivated: {len(deactivated_inequalities_set(model))})"
    )
    stats.append(
        f"{TAB}Activated Objectives: {len(activated_objectives_set(model))} "
        f"(Deactivated: {len(deactivated_objectives_set(model))})"
    )

    # Only show graybox info if they are present
    if len(greybox_block_set(model)) != 0:
        stats.append(f"{TAB}GreyBox Statistics")
        stats.append(
            f"{TAB* 2}Activated GreyBox models: {len(activated_greybox_block_set(model))} "
            f"(Deactivated: {len(deactivated_greybox_block_set(model))})"
        )
        stats.append(
            f"{TAB* 2}Activated GreyBox Equalities: {number_activated_greybox_equalities(model)} "
            f"(Deactivated: {number_deactivated_greybox_equalities(model)})"
        )
        stats.append(
            f"{TAB* 2}Free Variables in Activated GreyBox Equalities: {len(unfixed_greybox_variables(model))} (Fixed: {len(greybox_variables(model)-unfixed_greybox_variables(model))})"
        )

    return stats


def write_report_section(
    stream,
    lines_list,
    title=None,
    line_if_empty=None,
    end_line=None,
    header="-",
    footer=None,
):
    """
    Writes output in standard format for report and display methods.

    Args:
        stream: stream to write to
        lines_list: list containing lines to be written in body of report
        title: title to be put at top of report
        line_if_empty: line to be written if lines_list is empty
        end_line: line to be written at end of report
        header: character to use to write header separation line
        footer: character to use to write footer separation line

    Returns:
        None

    """
    stream.write(f"{header * MAX_STR_LENGTH}\n")
    if title is not None:
        stream.write(f"{title}\n\n")
    if len(lines_list) > 0:
        for i in lines_list:
            stream.write(f"{TAB}{i}\n")
    elif line_if_empty is not None:
        stream.write(f"{TAB}{line_if_empty}\n")
    stream.write("\n")
    if end_line is not None:
        stream.write(f"{end_line}\n")
    if footer is not None:
        stream.write(f"{footer * MAX_STR_LENGTH}\n")
