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
from pyomo.common.config import ConfigValue
from pyomo.core.base.constraint import Constraint
from pyomo.core.base.var import Var
from pyomo.core.base.units_container import units as pyunits

from idaes.core import (
    declare_process_block_class,
    MaterialBalanceType,
    UnitModelBlockData,
)
from idaes.core.util.config import is_physical_parameter_block

"""
A simple unit model to calculate power required to achieve isothermal
gas compression. 

Data sources:
    [1] Stochastic Optimal Control Model for Natural Gas Network
        Operations. V. Zavala, 2014, Comp. Chem. Eng.
"""


@declare_process_block_class("IsothermalCompressor")
class IsothermalCompressorData(UnitModelBlockData):

    CONFIG = UnitModelBlockData.CONFIG()

    CONFIG.declare(
        "property_package",
        ConfigValue(default=None, domain=is_physical_parameter_block),
    )

    def build(self):
        super(IsothermalCompressorData, self).build()

        time = self.flowsheet().time
        property_package = self.config.property_package

        if len(property_package.phase_list) != 1:
            raise ValueError(
                "%s can only be constructed with a "
                "single-phase property package.\n"
                "Got phases %s."
                % (self.__class__, [p for p in property_package.phase_list])
            )
        self.phase = next(iter(property_package.phase_list))
        if self.phase != "Vap":
            raise ValueError(
                '%s can only be constructed with a single phase, "Vap".'
                "Got phase %s." % (self.__class__, self.phase)
            )

        inlet_config = {"defined_state": True}
        self.inlet_state = property_package.build_state_block(time, **inlet_config)
        self.outlet_state = property_package.build_state_block(time)

        # A little annoying that add_port assumes the state block is indexed
        # by exactly time. If we supported unindexed state blocks here,
        # this entire unit model could be agnostic of time.
        # (This is not strictly true because add_state_material_balances
        # writes time-indexed balance equations, but we could get around this
        # by writing our own balance equations, which is easy enough.)
        self.add_port(
            name="inlet_port",
            block=self.inlet_state,
            doc="The inlet to the compressor",
        )
        self.add_port(
            name="outlet_port",
            block=self.outlet_state,
            doc="The outlet from the compressor",
        )

        self.add_state_material_balances(
            # TODO: Does componentPhase vs. componentTotal matter here?
            MaterialBalanceType.componentPhase,
            self.inlet_state,
            self.outlet_state,
        )

        self.add_state_isothermal_equation(self.inlet_state, self.outlet_state)

        # Add boost pressure and pressure equation
        self.boost_pressure = Var(
            time,
            initialize=0.0,
            units=pyunits.bar,
            # TODO: Should this variable not be initialized to its bound?
            bounds=(0.0, None),
            doc="Increase in pressure provided by this compressor",
        )
        self.add_pressure_change_equation(self.inlet_state, self.outlet_state)

        # Add compression coefficient
        self.beta = Var(
            time,
            initialize=0.28,  # Approximately (1.4 - 1.0)/1.4
            units=pyunits.dimensionless,
            bounds=(0.0, None),
            doc="Compression coefficient, beta.",
        )
        self.add_beta_equation(self.inlet_state)

        # Add power requirement
        self.power = Var(
            time,
            initialize=0.0,
            units=pyunits.kW,
            bounds=(0.0, None),
            doc="Work required to achieve the desired compression",
        )
        self.add_power_equation(self.inlet_state)

    def add_state_isothermal_equation(self, state1, state2):
        time = self.flowsheet().time
        # NOTE: We assume our states are only indexed by time and that
        # temperature is unindexed.
        def isothermal_rule(b, t):
            return state1[t].temperature == state2[t].temperature

        self.state_isothermal_eqn = Constraint(time, rule=isothermal_rule)

    def add_pressure_change_equation(self, inlet, outlet):
        time = self.flowsheet().time

        def pressure_rule(b, t):
            # Here we are using some quantities from the state blocks
            # and some from the unit model, so we must convert units.
            return pyunits.convert(inlet[t].pressure, pyunits.bar) + pyunits.convert(
                self.boost_pressure[t], pyunits.bar
            ) == pyunits.convert(outlet[t].pressure, pyunits.bar)

        self.pressure_change_eqn = Constraint(time, rule=pressure_rule)

    def add_beta_equation(self, inlet_state):
        time = self.flowsheet().time
        p = next(iter(self.config.property_package.phase_list))

        def beta_rule(b, t):
            # We know that we are constructing this property package with a
            # single phase, but we use heat_capacity_ratio_phase nonetheless
            # for compatibility with "generic" property packages.
            gamma = inlet_state[t].heat_capacity_ratio_phase[p]
            return b.beta[t] == (gamma - 1.0) / gamma

        self.beta_eqn = Constraint(time, rule=beta_rule)

    def add_power_equation(self, inlet_state):
        """
        This is Equation 2.10 in [1]
        """
        time = self.flowsheet().time

        def power_rule(b, t):
            cp_mass = inlet_state[t].cp_mol / inlet_state[t].mw
            power_expr = (
                inlet_state[t].temperature
                * cp_mass
                * inlet_state[t].flow_mass
                * (
                    (
                        (
                            pyunits.convert(inlet_state[t].pressure, pyunits.bar)
                            + pyunits.convert(self.boost_pressure[t], pyunits.bar)
                        )
                        / pyunits.convert(inlet_state[t].pressure, pyunits.bar)
                    )
                    ** b.beta[t]
                    - 1.0
                )
            )
            return b.power[t] == pyunits.convert(power_expr, pyunits.kW)

        self.power_eqn = Constraint(time, rule=power_rule)
