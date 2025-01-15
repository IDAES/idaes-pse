#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
0D Control Volume class with support for isothermal energy balance.
"""

__author__ = "Andrew Lee"

# Import IDAES cores
from idaes.core.base.control_volume0d import ControlVolume0DBlockData
from idaes.core import declare_process_block_class
from idaes.core.util.exceptions import ConfigurationError

import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


@declare_process_block_class(
    "ExtendedControlVolume0DBlock",
    doc="""
    ExtendedControlVolume0DBlock is an extension of the ControlVolume0D
    block with support for isothermal conditions in place of a formal
    energy balance.""",
)
class ExtendedControlVolume0DBlockData(ControlVolume0DBlockData):
    """
    Extended 0-Dimensional (Non-Discretized) ControlVolume Class

    This class extends the existing ControlVolume0DBlockData class
    with support for isothermal energy balances.
    """

    def add_isothermal_constraint(
        self,
        has_heat_of_reaction=False,
        has_heat_transfer=False,
        has_work_transfer=False,
        has_enthalpy_transfer=False,
        custom_term=None,
    ):
        """
        This method constructs an isothermal constraint for the control volume.

        Arguments are supported for compatibility with other forms but must be False
        or None otherwise an Exception is raised.

        Args:
            has_heat_of_reaction: whether terms for heat of reaction should
                be included in enthalpy balance
            has_heat_transfer: whether terms for heat transfer should be
                included in enthalpy balances
            has_work_transfer: whether terms for work transfer should be
                included in enthalpy balances
            has_enthalpy_transfer: whether terms for enthalpy transfer due to
                mass transfer should be included in enthalpy balance. This
                should generally be the same as the has_mass_transfer
                argument in the material balance methods
            custom_term: a Python method which returns Pyomo expressions representing
                custom terms to be included in enthalpy balances.
                Method should accept time and phase list as arguments.

        Returns:
            Constraint object representing enthalpy balances
        """
        if has_heat_transfer:
            raise ConfigurationError(
                f"{self.name}: isothermal energy balance option requires that has_heat_transfer is False. "
                "If you are trying to solve for heat duty to achieve isothermal operation, please use "
                "a full energy balance as add a constraint to equate inlet and outlet temperatures."
            )
        if has_work_transfer:
            raise ConfigurationError(
                f"{self.name}: isothermal energy balance option requires that has_work_transfer is False. "
                "If you are trying to solve for work under isothermal operation, please use "
                "a full energy balance as add a constraint to equate inlet and outlet temperatures."
            )
        if has_enthalpy_transfer:
            raise ConfigurationError(
                f"{self.name}: isothermal energy balance option does not support enthalpy transfer."
            )
        if has_heat_of_reaction:
            raise ConfigurationError(
                f"{self.name}: isothermal energy balance option requires that has_heat_of_reaction is False. "
                "If you are trying to solve for heat duty to achieve isothermal operation, please use "
                "a full energy balance as add a constraint to equate inlet and outlet temperatures."
            )
        if custom_term is not None:
            raise ConfigurationError(
                f"{self.name}: isothermal energy balance option does not support custom terms."
            )

        # Add isothermal constraint
        @self.Constraint(self.flowsheet().time, doc="Energy balances")
        def enthalpy_balances(b, t):
            return b.properties_in[t].temperature == b.properties_out[t].temperature
