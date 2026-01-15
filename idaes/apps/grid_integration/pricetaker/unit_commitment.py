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
This module contains functions that build unit commitment-type
constraints: startup/shutdown, uptime/downtime constraints,
capacity limit constraints, and ramping constraints.

The unit commitment model is taken from: 
Knueven, Bernard, James Ostrowski, and Jean-Paul Watson. 
"On mixed-integer programming formulations for the unit commitment problem." 
INFORMS Journal on Computing 32, no. 4 (2020): 857-876.
"""

from typing import Union, Optional
from pyomo.common.config import ConfigDict, ConfigValue
from pyomo.environ import Block, Constraint, Var, Param, RangeSet, Binary
from idaes.core.util.config import ConfigurationError, is_in_range


class UnitCommitmentData:
    """Dataclass to store startup, shutdown, and ramp rates"""

    def __init__(self, blk_name: str, commodity_name: str, **kwargs):
        self.blk_name = blk_name
        self.commodity_name = commodity_name
        self.config = self._get_config()
        self.update(**kwargs)

    def _get_config(self):
        config = ConfigDict()
        config.declare(
            "startup_rate",
            ConfigValue(
                domain=is_in_range(0, 1),
                doc="Startup rate as a fraction of design capacity",
            ),
        )
        config.declare(
            "shutdown_rate",
            ConfigValue(
                domain=is_in_range(0, 1),
                doc="Shutdown rate as a fraction of design capacity",
            ),
        )
        config.declare(
            "rampup_rate",
            ConfigValue(
                domain=is_in_range(0, 1),
                doc="Rampup rate as a fraction of design capacity",
            ),
        )
        config.declare(
            "rampdown_rate",
            ConfigValue(
                domain=is_in_range(0, 1),
                doc="Rampdown rate as a fraction of design capacity",
            ),
        )
        config.declare(
            "op_range_lb",
            ConfigValue(
                domain=is_in_range(0, 1),
                doc="Minimum stable operation range as a fraction of design capacity",
            ),
        )
        config.declare(
            "capacity",
            ConfigValue(
                doc=(
                    "Parameter/variable denoting the maximum "
                    "capacity of the commodity"
                ),
            ),
        )

        return config

    def update(self, **kwargs):
        """Updates the attribute values"""
        self.config.set_value(kwargs)
        self.assert_startup_rate_validity()
        self.assert_shutdown_rate_validity()

    def assert_startup_rate_validity(self):
        """Raises an error if startup rate value is not valid."""
        if None in (self.config.op_range_lb, self.config.startup_rate):
            # One of the values is not provided, so skip the check
            return

        if self.config.startup_rate < self.config.op_range_lb:
            raise ConfigurationError(
                f"For commidity {self.commodity_name} in operational "
                f"block {self.blk_name}, \n\tthe startup rate is less than "
                f"the minimum stable operation value."
            )

    def assert_shutdown_rate_validity(self):
        """Raises an error if shutdown rate value is not valid"""
        if None in (self.config.op_range_lb, self.config.shutdown_rate):
            # One of the values is not provided, so skip the check
            return

        if self.config.shutdown_rate < self.config.op_range_lb:
            raise ConfigurationError(
                f"For commidity {self.commodity_name} in operational "
                f"block {self.blk_name}, \n\tthe shutdown rate is less than "
                f"the minimum stable operation value."
            )

    def assert_ramping_args_present(self):
        """
        Raises an error if any of the arguments needed for the ramping
        constraints is missing.
        """
        cf = self.config
        if (
            None
            in (
                cf.startup_rate,
                cf.shutdown_rate,
                cf.rampup_rate,
                cf.rampdown_rate,
            )
            or cf.capacity is None
        ):
            raise ConfigurationError(
                "Necessary arguments needed for the ramping constraints are missing."
            )


def startup_shutdown_constraints(
    blk: Block,
    op_blocks: dict,
    install_unit: Union[int, Var],
    minimum_up_time: int,
    minimum_down_time: int,
    set_time: RangeSet,
    startup_transition_time: Optional[dict] = None,
):
    """
    Appends startup and shutdown constraints for a given unit/process.
    Supports multiples types of startup.

    Args:
        startup_transition_time (dict or None): A dictionary with keys as startup types and values are the time of startup transition.

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
        if t < minimum_up_time:
            return Constraint.Skip

        return (
            sum(op_blocks[i].startup for i in range(t - minimum_up_time + 1, t + 1))
            <= op_blocks[t].op_mode
        )

    @blk.Constraint(set_time)
    def minimum_down_time_con(_, t):
        if t < minimum_down_time:
            return Constraint.Skip

        return (
            sum(op_blocks[i].shutdown for i in range(t - minimum_down_time + 1, t + 1))
            <= install_unit - op_blocks[t].op_mode
        )

    if startup_transition_time is None or len(startup_transition_time) == 0:
        # if there is only one startup type, return
        # startup_transition_time can be None.
        # This is a double insurance check, that empty dict or None will skip the following.
        return

    # multiple startup types
    if startup_transition_time:
        # there will be at least two types of startup
        startup_names = list(startup_transition_time.keys())

        # assume the first should be max(min_down_time, startup_transition_time["hot"])
        startup_transition_time[startup_names[0]] = max(
            minimum_down_time, startup_transition_time[startup_names[0]]
        )

        # this is necessary, because we have updated the startup_transition_time.
        blk.startup_duration = Param(startup_names, initialize=startup_transition_time)

    @blk.Constraint(set_time)
    def tot_startup_type_rule(_, t):
        """
        Eq 55 in Knueven et.al.
        """

        return (
            sum(op_blocks[t].startup_type_vars[k] for k in startup_names)
            == op_blocks[t].startup
        )

    # add the startup type constraints for each type of startup
    for idx, key in enumerate(startup_names):
        if idx == 0:
            prev_key = key
            continue

        def startup_type_rule(_, t, key=key, prev_key=prev_key):
            """
            Eq 54 in Knueven et.al.
            """
            if t < blk.startup_duration[key]:
                return Constraint.Skip
            return op_blocks[t].startup_type_vars[prev_key] <= sum(
                op_blocks[t - i].shutdown
                for i in range(
                    blk.startup_duration[prev_key], blk.startup_duration[key]
                )
            )

        setattr(
            blk,
            f"Startup_Type_Constraint_{key}",
            Constraint(set_time, rule=startup_type_rule),
        )
        prev_key = key


def capacity_limits(
    blk: Block,
    op_blocks: dict,
    uc_data: UnitCommitmentData,
    set_time: RangeSet,
):
    """
    Appends capacity limit constraints
    """
    commodity = uc_data.commodity_name
    limits = (
        uc_data.config.op_range_lb * uc_data.config.capacity,
        uc_data.config.capacity,
    )
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
    uc_data: UnitCommitmentData,
    set_time: RangeSet,
):
    """
    Appends ramping constraints
    """
    _ramping_var = uc_data.commodity_name
    startup_rate = uc_data.config.startup_rate * uc_data.config.capacity
    shutdown_rate = uc_data.config.shutdown_rate * uc_data.config.capacity
    rampup_rate = uc_data.config.rampup_rate * uc_data.config.capacity
    rampdown_rate = uc_data.config.rampdown_rate * uc_data.config.capacity
    ramping_var = {t: getattr(blk, _ramping_var) for t, blk in op_blocks.items()}

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
