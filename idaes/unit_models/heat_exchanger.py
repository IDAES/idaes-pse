##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Heat Exchanger Models.
"""
from __future__ import division

# Import Pyomo libraries
from pyomo.environ import Reals, Var
from pyomo.common.config import ConfigBlock, ConfigValue, In

# Import IDAES cores
from idaes.core import (ControlVolume0D,
                        declare_process_block_class,
                        EnergyBalanceType,
                        MomentumBalanceType,
                        MaterialBalanceType,
                        UnitBlockData,
                        useDefault)
from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.misc import add_object_reference

__author__ = "John Eslick"

def _make_heater_control_volume(o, name, config):
    """
    This is seperated from the main model so it can be reused to different types
    of heat exchangers and heater/cooler units.
    """
    control_volume = ControlVolume0D(default={
            "dynamic": config.dynamic,
            "property_package": config.property_package,
            "property_package_args": config.property_package_args})
    # we have to attach this control volume to the model for the rest of
    # the steps to work
    setattr(o, name, control_volume)
    # Add inlet and outlet state blocks to control volume
    control_volume.add_state_blocks()
    # Add material balance
    control_volume.add_material_balances(
        balance_type=config.material_balance_type,
        has_phase_equilibrium=config.calculate_phase_equilibrium)
    # add energy balance
    control_volume.add_energy_balances(
        balance_type=config.energy_balance_type,
        has_heat_transfer=config.has_heat_transfer)
    # add momentum balance
    control_volume.add_momentum_balances(
        balance_type=config.momentum_balance_type,
        has_pressure_change=config.has_pressure_change)
    return control_volume

def _make_heater_config_block(config):
    """
    Declare options of a HeaterData unit.
    """
    config.declare("dynamic", ConfigValue(
        domain=In([True, False]),
        default=False,
        description="Dynamic model flag",
        doc="Indicates whether the model is dynamic."))
    config.declare("has_holdup", ConfigValue(
        default=useDefault,
        domain=In([useDefault, True, False]),
        description="Holdup construction flag",
        doc="""Indicates whether holdup terms should be constructed or not.
Must be True if dynamic = True,
**default** - False.
**Valid values:** {
**True** - construct holdup terms,
**False** - do not construct holdup terms}"""))
    config.declare("material_balance_type", ConfigValue(
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
    config.declare("energy_balance_type", ConfigValue(
        default=EnergyBalanceType.enthalpyTotal,
        domain=In(EnergyBalanceType),
        description="Energy balance construction flag",
        doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.enthalpyTotal.
**Valid values:** {
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single ethalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - ethalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}"""))
    config.declare("momentum_balance_type", ConfigValue(
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
    config.declare("has_heat_transfer", ConfigValue(
        default=True,
        domain=In([True, False]),
        description="Heat transfer term construction flag",
        doc="""Indicates whether terms for heat transfer should be constructed,
**default** - False.
**Valid values:** {
**True** - include heat transfer terms,
**False** - exclude heat transfer terms.}"""))
    config.declare("calculate_phase_equilibrium", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Calculate phase equilibrium in mixed stream",
        doc="""Argument indicating whether phase equilibrium should be
calculated for the resulting mixed stream,
**default** - False.
**Valid values:** {
**True** - calculate phase equilibrium in mixed stream,
**False** - do not calculate equilibrium in mixed stream.}"""))
    config.declare("has_pressure_change", ConfigValue(
        default=False,
        domain=In([True, False]),
        description="Pressure change term construction flag",
        doc="""Indicates whether terms for pressure change should be
constructed,
**default** - False.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}"""))
    config.declare("property_package", ConfigValue(
        default=useDefault,
        domain=is_physical_parameter_block,
        description="Property package to use for control volume",
        doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}"""))
    config.declare("property_package_args", ConfigBlock(
        implicit=True,
        description="Arguments to use for constructing property packages",
        doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}"""))


@declare_process_block_class("Heater", doc="Simple 0D heater/cooler model.")
class HeaterData(UnitBlockData):
    """
    Simple 0D heat exchange unit.

    Unit model to add or remove heat from a material.
    """
    CONFIG = ConfigBlock()
    _make_heater_config_block(CONFIG)

    def build(self):
        """
        Begin building model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(HeaterData, self).build()
        # Add Control Volume
        _make_heater_control_volume(self, "control_volume", self.config)
        # Add Ports
        self.add_inlet_port()
        self.add_outlet_port()
        # Add a convienient reference to heat duty.
        add_object_reference(self, "heat_duty", self.control_volume.heat)
