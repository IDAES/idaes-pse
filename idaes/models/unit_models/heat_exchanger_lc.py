#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Extension of the IDAES heat exchanger model for transient simulations
"""

__author__ = "Rusty Gentile, John Eslick, Andrew Lee"

from pyomo.environ import Var, Param
from pyomo.dae import DerivativeVar
from pyomo.common.config import ConfigValue
from idaes.core import declare_process_block_class, useDefault
from idaes.core.util.config import DefaultBool
from idaes.core.util.exceptions import ConfigurationError, IdaesError
from .heat_exchanger import HeatExchangerData


@declare_process_block_class(
    "HeatExchangerLumpedCapacitance", doc="0D heat exchanger for transient simulations"
)
class HeatExchangerLumpedCapacitanceData(HeatExchangerData):

    CONFIG = HeatExchangerData.CONFIG(implicit=True)

    CONFIG.declare(
        "dynamic_heat_balance",
        ConfigValue(
            default=useDefault,
            domain=DefaultBool,
            doc="""Indicates whether heat holdup in the wall material should 
be included in the overall energy balance,
**default** - useDefault.
**Valid values:** {
**useDefault** - get flag from parent (default = False),
**True** - include wall material heat holdup,
**False** - do not include wall material heat holdup.}""",
        ),
    )

    def _add_wall_variables(self):

        # Use the hot side as a reference
        s1_metadata = self.config.hot_side.property_package.get_metadata()

        # Unit system
        temp_units = s1_metadata.get_derived_units("temperature")
        energy_units = s1_metadata.get_derived_units("energy")
        u_units = s1_metadata.get_derived_units("heat_transfer_coefficient")
        area_units = s1_metadata.get_derived_units("area")
        ua_units = u_units * area_units
        r_units = 1 / ua_units

        self.temperature_wall = Var(
            self.flowsheet().config.time,
            initialize=300.0,
            doc="Average wall temperature",
            units=temp_units,
        )
        self.heat_capacity_wall = Param(
            initialize=1.0,
            mutable=True,
            doc="Total heat capacity of heat exchanger material",
            units=energy_units / temp_units,
        )
        self.ua_cold_side = Var(
            self.flowsheet().config.time,
            initialize=1000.0,
            doc="Overall heat transfer coefficient from the cold side",
            units=ua_units,
        )
        self.ua_hot_side = Var(
            self.flowsheet().config.time,
            initialize=1000.0,
            doc="Overall heat transfer coefficient from the hot side",
            units=ua_units,
        )
        self.thermal_resistance_wall = Param(
            initialize=0.0,
            mutable=True,
            doc="Total thermal resistance of heat exchanger material",
            units=r_units,
        )
        self.thermal_fouling_cold_side = Param(
            initialize=0.0,
            mutable=True,
            doc="Total thermal resistance due to fouling on the cold side",
            units=r_units,
        )
        self.thermal_fouling_hot_side = Param(
            initialize=0.0,
            mutable=True,
            doc="Total thermal resistance due to fouling on the hot side",
            units=r_units,
        )
        self.ua_hot_side_to_wall = Var(
            self.flowsheet().config.time,
            initialize=1000.0,
            doc="Overall heat transfer coefficient between the hot side and "
            "the center of the wall",
            units=ua_units,
        )

    def _add_wall_variable_constraints(self):
        @self.Constraint(
            self.flowsheet().config.time,
            doc="Overall heat transfer coefficient for wall temperature equation",
        )
        def ua_hot_side_to_wall_eq(b, t):
            return b.ua_hot_side_to_wall[t] == 1 / (
                1 / b.ua_hot_side[t]
                + b.thermal_fouling_hot_side
                + 0.5 * b.thermal_resistance_wall
            )

        @self.Constraint(
            self.flowsheet().config.time,
            doc="Overall heat transfer coefficient equation",
        )
        def ua_total_eq(b, t):
            return b.overall_heat_transfer_coefficient[t] * b.area == 1 / (
                1 / b.ua_hot_side[t]
                + b.thermal_fouling_hot_side
                + b.thermal_resistance_wall
                + b.thermal_fouling_cold_side
                + 1 / b.ua_cold_side[t]
            )

        @self.Constraint(self.flowsheet().config.time, doc="Wall temperature equation")
        def wall_temperature_eq(b, t):
            return (
                b.temperature_wall[t]
                == 0.5
                * (
                    b.hot_side.properties_in[t].temperature
                    + b.hot_side.properties_out[t].temperature
                )
                + b.hot_side.heat[t] / b.ua_hot_side_to_wall[t]
            )

    def activate_dynamic_heat_eq(self):
        """
        Activates the heat holdup term in the overall energy balance for
        transient simulations. Should only be used with dynamic flowsheets
        and if ``dynamic_heat_balance`` is True.

        Args:
            None
        Returns:
            None
        """

        if not self.config.dynamic_heat_balance:
            raise IdaesError(
                "{} heat holdup term cannot be activated "
                "when `dynamic_heat_balance=False`".format(self.name)
            )

        self.unit_heat_balance.deactivate()
        self.dynamic_heat_balance.activate()

    def deactivate_dynamic_heat_eq(self):
        """
        Removes the heat holdup term from the overall energy balance,
        restoring the behavior of the base 0D HeatExchanger model. This may
        be useful for more complex models which need to be solved in several
        stages.

        Args:
            None
        Returns:
            None
        """
        self.dynamic_heat_balance.deactivate()
        self.unit_heat_balance.activate()

    def build(self):
        """
        Building model

        Args:
            None
        Returns:
            None
        """

        super().build()
        self._add_wall_variables()
        self._add_wall_variable_constraints()

        if self.config.dynamic_heat_balance:

            s1_metadata = self.config.hot_side.property_package.get_metadata()
            temp_units = s1_metadata.get_derived_units("temperature")
            time_units = s1_metadata.get_derived_units("time")

            if not self.flowsheet().config.dynamic:
                raise ConfigurationError(
                    "{} dynamic heat balance cannot be "
                    "used with a steady-state "
                    "flowsheet".format(self.name)
                )

            self.dT_wall_dt = DerivativeVar(
                self.temperature_wall,
                wrt=self.flowsheet().config.time,
                doc="Derivative of wall temperature with respect to time",
                units=temp_units / time_units,
            )

            @self.Constraint(
                self.flowsheet().time,
                doc="Unit heat balance equation with the added capacitance term",
            )
            def dynamic_heat_balance(b, t):
                return (
                    b.hot_side.heat[t] + b.cold_side.heat[t]
                    == -b.dT_wall_dt[t] * b.heat_capacity_wall
                )

            self.activate_dynamic_heat_eq()
            self.dT_wall_dt[self.flowsheet().config.time.first()].fix(0)
            self.dT_wall_dt[:].value = 0

    def initialize(self, *args, **kwargs):

        if self.config.dynamic_heat_balance:
            # If the time derivative terms are defined, deactivate them first
            # and reactivate after calling the base initialization method
            self.deactivate_dynamic_heat_eq()
            super().initialize(*args, **kwargs)
            self.activate_dynamic_heat_eq()
        else:
            # Otherwise, just initialize as normal
            super().initialize(*args, **kwargs)
