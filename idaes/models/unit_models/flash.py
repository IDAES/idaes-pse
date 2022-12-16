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
Standard IDAES flash model.
"""
# Import Python libraries
import logging
from pandas import DataFrame

# Import Pyomo libraries
from pyomo.environ import Constraint, Reference
from pyomo.common.config import ConfigBlock, ConfigValue, In, Bool
from pyomo.network import Port

# Import IDAES cores
from idaes.core import (
    ControlVolume0DBlock,
    declare_process_block_class,
    MaterialBalanceType,
    EnergyBalanceType,
    MomentumBalanceType,
    UnitModelBlockData,
    useDefault,
)
from idaes.models.unit_models.separator import (
    Separator,
    SplittingType,
    EnergySplittingType,
)

from idaes.core.util.config import is_physical_parameter_block
from idaes.core.util.units_of_measurement import report_quantity


__author__ = "Andrew Lee, Jaffer Ghouse"


# Set up logger
logger = logging.getLogger("idaes.unit_model")


@declare_process_block_class("Flash")
class FlashData(UnitModelBlockData):
    """
    Standard Flash Unit Model Class
    """

    CONFIG = ConfigBlock()
    CONFIG.declare(
        "dynamic",
        ConfigValue(
            domain=In([False]),
            default=False,
            description="Dynamic model flag - must be False",
            doc="""Indicates whether this model will be dynamic or not,
**default** = False. Flash units do not support dynamic behavior.""",
        ),
    )
    CONFIG.declare(
        "has_holdup",
        ConfigValue(
            default=False,
            domain=In([False]),
            description="Holdup construction flag - must be False",
            doc="""Indicates whether holdup terms should be constructed or not.
**default** - False. Flash units do not have defined volume, thus
this must be False.""",
        ),
    )
    CONFIG.declare(
        "material_balance_type",
        ConfigValue(
            default=MaterialBalanceType.useDefault,
            domain=In(MaterialBalanceType),
            description="Material balance construction flag",
            doc="""Indicates what type of mass balance should be constructed,
**default** - MaterialBalanceType.useDefault.
**Valid values:** {
**MaterialBalanceType.useDefault - refer to property package for default
balance type
**MaterialBalanceType.none** - exclude material balances,
**MaterialBalanceType.componentPhase** - use phase component balances,
**MaterialBalanceType.componentTotal** - use total component balances,
**MaterialBalanceType.elementTotal** - use total element balances,
**MaterialBalanceType.total** - use total material balance.}""",
        ),
    )
    CONFIG.declare(
        "energy_balance_type",
        ConfigValue(
            default=EnergyBalanceType.useDefault,
            domain=In(EnergyBalanceType),
            description="Energy balance construction flag",
            doc="""Indicates what type of energy balance should be constructed,
**default** - EnergyBalanceType.useDefault.
**Valid values:** {
**EnergyBalanceType.useDefault - refer to property package for default
balance type
**EnergyBalanceType.none** - exclude energy balances,
**EnergyBalanceType.enthalpyTotal** - single enthalpy balance for material,
**EnergyBalanceType.enthalpyPhase** - enthalpy balances for each phase,
**EnergyBalanceType.energyTotal** - single energy balance for material,
**EnergyBalanceType.energyPhase** - energy balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "momentum_balance_type",
        ConfigValue(
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
**MomentumBalanceType.momentumPhase** - momentum balances for each phase.}""",
        ),
    )
    CONFIG.declare(
        "energy_split_basis",
        ConfigValue(
            default=EnergySplittingType.equal_temperature,
            domain=EnergySplittingType,
            description="Type of constraint to write for energy splitting",
            doc="""Argument indicating basis to use for splitting energy this is
not used for when ideal_separation == True.
**default** - EnergySplittingType.equal_temperature.
**Valid values:** {
**EnergySplittingType.equal_temperature** - outlet temperatures equal inlet
**EnergySplittingType.equal_molar_enthalpy** - oulet molar enthalpies equal
inlet,
**EnergySplittingType.enthalpy_split** - apply split fractions to enthalpy
flows.}""",
        ),
    )
    CONFIG.declare(
        "ideal_separation",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Ideal splitting flag",
            doc="""Argument indicating whether ideal splitting should be used.
Ideal splitting assumes perfect separation of material, and attempts to
avoid duplication of StateBlocks by directly partitioning outlet flows to
ports,
**default** - True.
**Valid values:** {
**True** - use ideal splitting methods. Cannot be combined with
has_phase_equilibrium = True,
**False** - use explicit splitting equations with split fractions.}""",
        ),
    )
    CONFIG.declare(
        "has_heat_transfer",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Heat transfer term construction flag",
            doc="""Indicates whether terms for heat transfer should be constructed,
**default** - False.
**Valid values:** {
**True** - include heat transfer terms,
**False** - exclude heat transfer terms.}""",
        ),
    )
    CONFIG.declare(
        "has_pressure_change",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Pressure change term construction flag",
            doc="""Indicates whether terms for pressure change should be
constructed,
**default** - True.
**Valid values:** {
**True** - include pressure change terms,
**False** - exclude pressure change terms.}""",
        ),
    )
    CONFIG.declare(
        "property_package",
        ConfigValue(
            default=useDefault,
            domain=is_physical_parameter_block,
            description="Property package to use for control volume",
            doc="""Property parameter object used to define property calculations,
**default** - useDefault.
**Valid values:** {
**useDefault** - use default package from parent model or flowsheet,
**PropertyParameterObject** - a PropertyParameterBlock object.}""",
        ),
    )
    CONFIG.declare(
        "property_package_args",
        ConfigBlock(
            implicit=True,
            description="Arguments to use for constructing property packages",
            doc="""A ConfigBlock with arguments to be passed to a property block(s)
and used when constructing these,
**default** - None.
**Valid values:** {
see property package for documentation.}""",
        ),
    )

    def build(self):
        """
        Begin building model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        # Call UnitModel.build to setup dynamics
        super(FlashData, self).build()

        # Build Control Volume
        self.control_volume = ControlVolume0DBlock(
            dynamic=self.config.dynamic,
            has_holdup=self.config.has_holdup,
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
        )

        self.control_volume.add_state_blocks(has_phase_equilibrium=True)

        self.control_volume.add_material_balances(
            balance_type=self.config.material_balance_type, has_phase_equilibrium=True
        )

        self.control_volume.add_energy_balances(
            balance_type=self.config.energy_balance_type,
            has_heat_transfer=self.config.has_heat_transfer,
        )

        self.control_volume.add_momentum_balances(
            balance_type=self.config.momentum_balance_type,
            has_pressure_change=self.config.has_pressure_change,
        )

        # Add Ports
        self.add_inlet_port()

        split_map = {}
        for p in self.control_volume.properties_in.phase_list:
            p_obj = self.config.property_package.get_phase(p)
            if p_obj.is_vapor_phase():
                # Vapor leaves through Vap outlet
                split_map[p] = "Vap"
            else:
                # All other phases leave through Liq outlet
                split_map[p] = "Liq"

        self.split = Separator(
            property_package=self.config.property_package,
            property_package_args=self.config.property_package_args,
            outlet_list=["Vap", "Liq"],
            split_basis=SplittingType.phaseFlow,
            ideal_separation=self.config.ideal_separation,
            ideal_split_map=split_map,
            mixed_state_block=self.control_volume.properties_out,
            has_phase_equilibrium=not self.config.ideal_separation,
            energy_split_basis=self.config.energy_split_basis,
        )
        if not self.config.ideal_separation:

            def split_frac_rule(b, t, o):
                return b.split.split_fraction[t, o, o] == 1

            self.split_fraction_eq = Constraint(
                self.flowsheet().time, self.split.outlet_idx, rule=split_frac_rule
            )

        self.vap_outlet = Port(extends=self.split.Vap)
        self.liq_outlet = Port(extends=self.split.Liq)

        # Add references
        if (
            self.config.has_heat_transfer is True
            and self.config.energy_balance_type != EnergyBalanceType.none
        ):
            self.heat_duty = Reference(self.control_volume.heat[:])
        if (
            self.config.has_pressure_change is True
            and self.config.momentum_balance_type != MomentumBalanceType.none
        ):
            self.deltaP = Reference(self.control_volume.deltaP[:])

    def _get_performance_contents(self, time_point=0):
        var_dict = {}
        if hasattr(self, "heat_duty"):
            var_dict["Heat Duty"] = self.heat_duty[time_point]
        if hasattr(self, "deltaP"):
            var_dict["Pressure Change"] = self.deltaP[time_point]

        return {"vars": var_dict}

    def _get_stream_table_contents(self, time_point=0):
        stream_attributes = {}
        stream_attributes["Units"] = {}

        sblocks = {
            "Inlet": self.control_volume.properties_in,
        }
        if not self.config.ideal_separation:
            # If not using ideal separation, we can get outlet state directly
            # from the state blocks
            sblocks["Vapor Outlet"] = self.split.Vap_state
            sblocks["Liquid Outlet"] = self.split.Liq_state

        for n, v in sblocks.items():
            dvars = v[time_point].define_display_vars()

            stream_attributes[n] = {}

            for k in dvars:
                for i in dvars[k].keys():
                    stream_key = k if i is None else f"{k} {i}"

                    quant = report_quantity(dvars[k][i])

                    stream_attributes[n][stream_key] = quant.m
                    stream_attributes["Units"][stream_key] = quant.u

        if self.config.ideal_separation:
            # If using ideal separation, get values from Ports and hope they map
            # to names in Inlet
            # TODO: Add a better way to map these if necessary
            for n, v in {
                "Vapor Outlet": "vap_outlet",
                "Liquid Outlet": "liq_outlet",
            }.items():
                port_obj = getattr(self, v)

                stream_attributes[n] = {}

                for k in port_obj.vars:
                    for i in port_obj.vars[k].keys():
                        if isinstance(i, float):
                            quant = report_quantity(port_obj.vars[k][time_point])
                            stream_attributes[n][k] = quant.m
                            stream_attributes["Units"][k] = quant.u
                        else:
                            if len(i) == 2:
                                kname = str(i[1])
                            else:
                                kname = str(i[1:])
                            quant = report_quantity(port_obj.vars[k][time_point, i[1:]])
                            stream_attributes[n][k + " " + kname] = quant.m
                            stream_attributes["Units"][k + " " + kname] = quant.u

        return DataFrame.from_dict(stream_attributes, orient="columns")
