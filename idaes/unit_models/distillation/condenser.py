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
Condenser model for distillation
"""

__author__ = "Jaffer Ghouse"

import logging
from enum import Enum

# Import Pyomo libraries
from pyomo.common.config import ConfigBlock, ConfigValue, In
from pyomo.network import Port
from pyomo.environ import Reference, Expression, Var, Constraint, \
    TerminationCondition

# Import IDAES cores
from idaes.core import (ControlVolume0DBlock,
                        declare_process_block_class,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        MaterialBalanceType,
                        UnitModelBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import add_object_reference

_log = logging.getLogger(__name__)


class CondenserType(Enum):
    totalCondenser = 0
    partialCondenser = 1


@declare_process_block_class("Condenser")
class CondenserData(UnitModelBlockData):
    """
    Condenser unit for distillation model.
    Unit model to condense (total/partial) the vapor from the top tray of
    the distillation column.
    """
    CONFIG = UnitModelBlockData.CONFIG()
    CONFIG.declare("condenser_type", ConfigValue(
        default=CondenserType.totalCondenser,
        domain=In(CondenserType),
        description="Type of condenser flag",
        doc="""Indicates what type of condenser should be constructed,
**default** - CondenserType.totalCondenser.
**Valid values:** {
**CondenserType.totalCondenser** - Incoming vapor from top tray is condensed
to all liquid,
**CondenserType.partialCondenser** - Incoming vapor from top tray is
partially condensed to a vapor and liquid stream.}"""))
    CONFIG.declare("material_balance_type", ConfigValue(
        default=MaterialBalanceType.componentPhase,
        domain=In(MaterialBalanceType),
        description="Material balance construction flag",
        doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.componentPhase.
**Valid values:** {
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}"""))
    CONFIG.declare("energy_balance_type", ConfigValue(
        default=EnergyBalanceType.enthalpyTotal,
        domain=In(EnergyBalanceType),
        description="Energy balance construction flag",
        doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.enthalpyTotal.
**Valid values:** {
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}"""))
    CONFIG.declare("momentum_balance_type", ConfigValue(
        default=MomentumBalanceType.pressureTotal,
        domain=In(MomentumBalanceType),
        description="Momentum balance construction flag",
        doc="""Indicates what type of momentum balance should be constructed,
**default** - MomentumBalanceType.pressureTotal.
**Valid values:** {
**MomentumBalanceType.none** - exclude momentum balances,
**MomentumBalanceType.pressureTotal** - single pressure balance for material,
**MomentumBalanceType.pressurePhase** - pressure balances for each phase,
**MomentumBalanceType.momentumTotal** - single momentum balance for material,
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}"""))
    CONFIG.declare("has_pressure_change", ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}"""))
    CONFIG.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}"""))
    CONFIG.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))

    def build(self):
        """Build the model.

        Args:
            None
        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(CondenserData, self).build()

        # Add Control Volume for the condenser
        self.control_volume = ControlVolume0DBlock(default={
            "dynamic": self.config.dynamic,
            "has_holdup": self.config.has_holdup,
            "property_package": self.config.property_package,
            "property_package_args": self.config.property_package_args})

        self.control_volume.add_state_blocks(
            has_phase_equilibrium=True)

        self.control_volume.add_material_balances(
            balance_type=self.config.material_balance_type,
            has_phase_equilibrium=True)

        self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=True)

        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change)

        self._make_ports()

        if self.config.condenser_type == CondenserType.totalCondenser:

            self._make_splits_total_condenser()

            # Set condition for total condenser (T_cond = T_bubble) (Option 1)

            def rule_total_cond(self, t):
                return self.control_volume.properties_out[t].temperature == \
                    self.control_volume.properties_out[t].temperature_bubble
            self.eq_total_cond_spec = Constraint(self.flowsheet().time,
                                                 rule=rule_total_cond)
        else:
            self._make_splits_partial_condenser()

        # Add object reference to variables of the control volume
        # Reference to the heat duty
        add_object_reference(self, "heat_duty", self.control_volume.heat)
        # Reference to the pressure drop
        add_object_reference(self, "deltaP", self.control_volume.deltaP)

    def _make_ports(self):

        # Add Ports for the condenser
        # Inlet port (the vapor from the top tray)
        self.add_inlet_port()

        # Outlet ports that always exist irrespective of condenser type
        self.reflux = Port(noruleinit=True, doc="Reflux stream that is"
                           " returned to the top tray.")
        self.distillate = Port(noruleinit=True, doc="Reflux stream that is"
                               " returned to the top tray.")

        if self.config.condenser_type == CondenserType.partialCondenser:
            self.vapor_outlet = Port(noruleinit=True,
                                     doc="Vapor outlet port from a "
                                     "partial condenser")
        # Add codnenser specific variables
        self.reflux_ratio = Var(initialize=1, doc="reflux ratio for "
                                "the condenser")

    def _make_splits_total_condenser(self):
        # Get dict of Port members and names
        member_list = self.control_volume.\
            properties_out[0].define_port_members()

        if self.config.condenser_type == CondenserType.totalCondenser:
            # Create references and populate the reflux, distillate ports
            for k in member_list:
                # Create references and populate the intensive variables
                if "flow" not in k:
                    if not member_list[k].is_indexed():
                        var = self.control_volume.properties_out[:].\
                            component(member_list[k].local_name)
                    else:
                        var = self.control_volume.properties_out[:].\
                            component(member_list[k].local_name)[...]

                    # add the reference and variable name to the reflux port
                    self.reflux.add(Reference(var), k)

                    # add the reference and variable name to the distillate port
                    self.distillate.add(Reference(var), k)
                else:
                    # Create references and populate the extensive variables
                    # This is for vars that are not indexed
                    if not member_list[k].is_indexed():
                        # Expression for reflux flow and relation to the
                        # reflux_ratio variable
                        def rule_reflux_flow(self, t):
                            return self.control_volume.properties_out[t].\
                                component(member_list[k].local_name) * \
                                (self.reflux_ratio / (1 + self.reflux_ratio))
                        self.e_reflux_flow = Expression(self.flowsheet().time,
                                                        rule=rule_reflux_flow)
                        self.reflux.add(self.e_reflux_flow, k)

                        # Expression for distillate flow and relation to the
                        # reflux_ratio variable
                        def rule_distillate_flow(self, t):
                            return self.control_volume.properties_out[t].\
                                component(member_list[k].local_name) / \
                                (1 + self.reflux_ratio)
                        self.e_distillate_flow = Expression(
                            self.flowsheet().time, rule=rule_distillate_flow)
                        self.distillate.add(self.e_distillate_flow, k)
                    else:
                        # Create references and populate the extensive variables
                        # This is for vars that are indexed
                        index_set = member_list[k].index_set()

                        def rule_reflux_flow(self, t, *args):
                            return self.control_volume.properties_out[t].\
                                component(member_list[k].local_name)[args] * \
                                (self.reflux_ratio / (1 + self.reflux_ratio))
                        self.e_reflux_flow = Expression(
                            (self.flowsheet().time, index_set),
                            rule=rule_reflux_flow)
                        self.reflux.add(self.e_reflux_flow, k)

                        # Create references and populate the extensive variables
                        # This is for vars that are indexed
                        index_set = member_list[k].index_set()

                        def rule_distillate_flow(self, t, *args):
                            return self.control_volume.properties_out[t].\
                                component(member_list[k].local_name)[args] / \
                                (1 + self.reflux_ratio)
                        self.e_distillate_flow = Expression(
                            (self.flowsheet().time, index_set),
                            rule=rule_distillate_flow)
                        self.distillate.add(self.e_distillate_flow, k)

    def _make_splits_partial_condenser(self):
        # Get dict of Port members and names
        member_list = self.control_volume.\
            properties_out[0].define_port_members()
        if self.config.condenser_type == CondenserType.partialCondenser:
            # Create references and populate the reflux, distillate ports
            for k in member_list:
                # Create references and populate the intensive variables
                if "flow" not in k and "frac" not in k and "enth" not in k:
                    if not member_list[k].is_indexed():
                        var = self.control_volume.properties_out[:].\
                            component(member_list[k].local_name)
                    else:
                        var = self.control_volume.properties_out[:].\
                            component(member_list[k].local_name)[...]

                    # add the reference and variable name to the reflux port
                    self.reflux.add(Reference(var), k)

                    # add the reference and variable name to the distillate port
                    self.distillate.add(Reference(var), k)

                    # add the reference and variable name to the
                    # vapor outlet port
                    self.vapor_outlet.add(Reference(var), k)
                elif "mass_frac" in k or "mole_frac" in k:
                    # Mole/mass frac is indexed
                    index_set = member_list[k].index_set()

                    # if state var is not mole/mass frac by phase
                    if "phase" not in k:
                        local_name = str(member_list[k].local_name) + \
                            "_phase"

                        def rule_liq_frac(self, t, i):
                            return self.control_volume.properties_out[t].\
                                component(local_name)["Liq", i]
                        self.e_liq_frac = Expression(
                            self.flowsheet().time, index_set,
                            rule=rule_liq_frac)

                        def rule_vap_frac(self, t, i):
                            return self.control_volume.properties_out[t].\
                                component(local_name)["Vap", i]
                        self.e_vap_frac = Expression(
                            self.flowsheet().time, index_set,
                            rule=rule_vap_frac)

                        # add the reference and variable name to the reflux port
                        self.reflux.add(self.e_liq_frac, k)

                        # add the reference and variable name to the
                        # distillate port
                        self.distillate.add(self.e_liq_frac, k)

                        # add the reference and variable name to the
                        # distillate port
                        self.vapor_outlet.add(self.e_vap_frac, k)
                    else:

                        var = self.control_volume.properties_out[:].\
                            component(member_list[k].local_name)[...]

                        # add the reference and variable name to the reflux port
                        self.reflux.add(Reference(var), k)

                        # add the reference and variable name to the distillate port
                        self.distillate.add(Reference(var), k)
                elif "flow" in k:
                    if "comp" not in k or "phase" not in k:

                        # if state var is not mole/mass frac by phase
                        local_name = str(member_list[k].local_name) + \
                            "_phase"

                        def rule_vap_flow(self, t):
                            return self.control_volume.properties_out[t].\
                                component(local_name)["Vap"]
                        self.e_vap_flow = Expression(
                            self.flowsheet().time,
                            rule=rule_vap_flow)

                        def rule_reflux_flow(self, t):
                            return self.control_volume.properties_out[t].\
                                component(local_name)["Liq"] * \
                                (self.reflux_ratio / (1 + self.reflux_ratio))
                        self.e_reflux_flow = Expression(
                            self.flowsheet().time,
                            rule=rule_reflux_flow)

                        def rule_distillate_flow(self, t):
                            return self.control_volume.properties_out[t].\
                                component(local_name)["Liq"] / \
                                (1 + self.reflux_ratio)
                        self.e_distillate_flow = Expression(
                            self.flowsheet().time,
                            rule=rule_distillate_flow)

                        # add the reference and variable name to the reflux port
                        self.reflux.add(self.e_reflux_flow, k)

                        # add the reference and variable name to the
                        # distillate port
                        self.distillate.add(self.e_distillate_flow, k)

                        # add the reference and variable name to the
                        # distillate port
                        self.vapor_outlet.add(self.e_vap_flow, k)

    def initialize(self, solver=None, outlvl=None):

        # TODO: Fix the inlets to the condenser to the vapor flow from
        # the top tray or take it as an argument to this method.
        if self.config.condenser_type == CondenserType.totalCondenser:
            self.eq_total_cond_spec.deactivate()

        # Initialize the inlet and outlet state blocks
        self.control_volume.initialize(outlvl=outlvl)

        # Activate the total condenser spec
        if self.config.condenser_type == CondenserType.totalCondenser:
            self.eq_total_cond_spec.activate()

        if solver is not None:
            if outlvl > 2:
                tee = True
            else:
                tee = False

            solver_output = solver.solve(self, tee=tee)
            if solver_output.solver.termination_condition == \
                    TerminationCondition.optimal:
                _log.info('{} Condenser Initialisation Complete.'
                          .format(self.name))
