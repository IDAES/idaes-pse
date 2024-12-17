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
This module contains functions that build unit commitment-type
constraints: startup/shutdown, uptime/downtime constraints,
capacity limit constraints, and ramping constraints.
"""

from typing import Union, Any
from pyomo.environ import Block, Constraint, Var, RangeSet


def startup_shutdown_constraints(
    blk: Block,
    op_blocks: dict,
    install_unit: Union[int, Var],
    up_time: int,
    down_time: int,
    set_time: RangeSet,
):
    """
    Appends startup and shutdown constraints for a given unit/process
    """

    @blk.Constraint(set_time)
    def binary_relationship_con(_, t):
        if t == 1:
            return Constraint.Skip

        return (
            op_blocks[t].op_mode - op_blocks[t - 1].op_mode
            == op_blocks[t].startup - op_blocks[t].shutdown
        )

    @blk.Constraint(set_time)
    def minimum_up_time_con(_, t):
        if t < up_time:
            return Constraint.Skip

        return (
            sum(op_blocks[i].startup for i in range(t - up_time + 1, t + 1))
            <= op_blocks[t].op_mode
        )

    @blk.Constraint(set_time)
    def minimum_down_time_con(_, t):
        if t < down_time:
            return Constraint.Skip

        return (
            sum(op_blocks[i].shutdown for i in range(t - down_time + 1, t + 1))
            <= install_unit - op_blocks[t].op_mode
        )


def capacity_limits(
    blk: Block,
    op_blocks: dict,
    commodity: str,
    limits: tuple,
    set_time: RangeSet,
):
    """
    Appends capacity limit constraints
    """
    commodity = {t: getattr(blk, commodity) for t, blk in op_blocks.items()}

    @blk.Constraint(set_time)
    def capacity_low_limit_con(_, t):
        return limits[0] * op_blocks[t].op_mode <= commodity[t]

    @blk.Constraint(set_time)
    def capacity_high_limit_con(_, t):
        return commodity[t] <= limits[1] * op_blocks[t].op_mode


def ramping_limits(
    blk: Block,
    op_blocks: dict,
    ramping_var: str,
    startup_rate: Any,
    shutdown_rate: Any,
    rampup_rate: Any,
    rampdown_rate: Any,
    set_time: RangeSet,
):
    """
    Appends ramping constraints
    """
    ramping_var = {t: getattr(blk, ramping_var) for t, blk in op_blocks.items()}

    @blk.Constraint(set_time)
    def ramp_up_con(_, t):
        if t == 1:
            return Constraint.Skip

        return (
            ramping_var[t] - ramping_var[t - 1]
            <= startup_rate * op_blocks[t].startup
            + rampup_rate * op_blocks[t - 1].op_mode
        )

    @blk.Constraint(set_time)
    def ramp_down_con(_, t):
        if t == 1:
            return Constraint.Skip

        return (
            ramping_var[t - 1] - ramping_var[t]
            <= shutdown_rate * op_blocks[t].shutdown
            + rampdown_rate * op_blocks[t].op_mode
        )
