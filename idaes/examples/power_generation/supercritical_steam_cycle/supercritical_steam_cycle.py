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

# Import Python libraries
import os
import logging
import sys
import csv
from collections import OrderedDict
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import argparse

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.opt import TerminationCondition
from pyomo.common.fileutils import this_file_dir
from pyomo.network import Arc, Port
from pyomo.core.kernel.component_set import ComponentSet

# Import IDAES
from idaes.core import FlowsheetBlock, MaterialBalanceType
from idaes.core.util import model_serializer as ms
from idaes.core.util import copy_port_values as _set_port
from idaes.core.util.misc import svg_tag
from idaes.core.util.tables import create_stream_table_dataframe
from idaes.unit_models.power_generation import TurbineMultistage, FWH0D
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.unit_models import (
    Mixer,
    HeatExchanger,
    PressureChanger,
    MomentumMixingType,
    Heater,
)
from idaes.unit_models.heat_exchanger import delta_temperature_underwood_callback
from idaes.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.property_models import iapws95
from idaes.core.util import model_serializer as ms


def create_model():
    """Create the flowsheet and add unit models. Fixing model inputs is done
    in a separate function, to try to keep this fairly clean and easy to follow.

    Args:
        None

    Returns:
        (ConcreteModel) Steam cycle model
    """
    ############################################################################
    #  Flowsheet and Properties                                                #
    ############################################################################
    m = pyo.ConcreteModel(name="Steam Cycle Model")
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.prop_water = iapws95.Iapws95ParameterBlock(
        default={"phase_presentation": iapws95.PhaseType.LG}
    )
    m.fs.prop_water_tpx = iapws95.Iapws95ParameterBlock(
        default={
            "phase_presentation": iapws95.PhaseType.LG,
            "state_vars": iapws95.StateVars.TPX,
        }
    )
    ############################################################################
    #  Turbine with fill-in reheat constraints                                 #
    ############################################################################
    m.fs.turb = TurbineMultistage(
        default={
            "property_package": m.fs.prop_water,
            "num_parallel_inlet_stages": 4,  # partial arc admission
            "num_hp": 7,
            "num_ip": 10,
            "num_lp": 11,
            "hp_split_locations": [4, 7],  # hp steam extraction locations
            "ip_split_locations": [5, 10],  # ip steam extraction locations
            "lp_split_locations": [4, 8, 10, 11],  # lp steam extraction locations
            "hp_disconnect": [7],  # disconnect hp from ip to insert reheater
            "ip_split_num_outlets": {10: 3},
        }
    )  # default is 2, so only set this
    # This is only the steam cycle, so don't have the flue gas side of the
    # reheater. To fill in I'll just add a few constraints for the flow,
    # pressure drop, and outlet temperature. Can insert a unit model later
    @m.fs.turb.Constraint(m.fs.time)
    def constraint_reheat_flow(b, t):
        return b.ip_stages[1].inlet.flow_mol[t] == b.hp_split[7].outlet_1.flow_mol[t]

    m.fs.turb.reheat_delta_p = pyo.Var(m.fs.time)
    m.fs.turb.reheat_delta_p.fix(0)

    @m.fs.turb.Constraint(m.fs.time)
    def constraint_reheat_press(b, t):
        return (
            b.ip_stages[1].inlet.pressure[t]
            == b.hp_split[7].outlet_1.pressure[t] + b.reheat_delta_p[t]
        )

    m.fs.turb.reheat_out_T = pyo.Var(m.fs.time, initialize=866)
    m.fs.turb.reheat_out_T.fix()

    @m.fs.turb.Constraint(m.fs.time)
    def constraint_reheat_temp(b, t):
        return (
            b.ip_stages[1].control_volume.properties_in[t].temperature
            == b.reheat_out_T[t]
        )

    ############################################################################
    #  Add Condenser/hotwell/condensate pump                                   #
    ############################################################################
    # condenser hx
    m.fs.condenser_mix = Mixer(
        default={
            "momentum_mixing_type": MomentumMixingType.none,
            "material_balance_type": MaterialBalanceType.componentTotal,
            "inlet_list": ["main", "bfpt"],
            "property_package": m.fs.prop_water,
        }
    )
    # Add an extra port to the mixer outlet that will let me hook it to the
    # condenser heatexchanger which uses T,P,x state variables because it is
    # easier solve the condenser pressure calculation with those vars.
    m.fs.condenser_mix._outlet_temperature_ref = pyo.Reference(
        m.fs.condenser_mix.mixed_state[:].temperature
    )
    m.fs.condenser_mix._outlet_vapor_fraction_ref = pyo.Reference(
        m.fs.condenser_mix.mixed_state[:].vapor_frac
    )
    m.fs.condenser_mix.outlet_tpx = Port(
        initialize={
            "flow_mol": m.fs.condenser_mix._outlet_flow_mol_ref,
            "temperature": m.fs.condenser_mix._outlet_temperature_ref,
            "pressure": m.fs.condenser_mix._outlet_pressure_ref,
            "vapor_frac": m.fs.condenser_mix._outlet_vapor_fraction_ref,
        }
    )

    @m.fs.condenser_mix.Constraint(m.fs.time)
    def mixer_pressure_constraint(b, t):
        return b.main_state[t].pressure == b.mixed_state[t].pressure

    m.fs.condenser = HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_underwood_callback,
            "shell": {
                "property_package": m.fs.prop_water_tpx,
                "material_balance_type": MaterialBalanceType.componentTotal,
            },
            "tube": {
                "property_package": m.fs.prop_water,
                "material_balance_type": MaterialBalanceType.componentTotal,
            },
        }
    )
    m.fs.condenser.delta_temperature_out.fix(5)
    # I know everything condenses so want to calculate pressure
    m.fs.condenser.shell.properties_out[:].eq_complementarity.deactivate()

    m.fs.condenser.pressure_over_sat = pyo.Var(m.fs.time,
        doc="Pressure added to Psat in the condeser. This is to account for"
            "some subcooling. (Pa)")
    m.fs.condenser.pressure_over_sat.fix(2000)
    @m.fs.condenser.Constraint(m.fs.time)
    def eq_pressure(b, t):
        return (
            b.shell.properties_out[t].pressure
            == b.shell.properties_out[t].pressure_sat + b.pressure_over_sat[t]
        )

    # m.fs.condenser.heat_transfer_equation.deactivate()
    m.fs.condenser.shell.properties_out[:].vapor_frac.fix(0)
    # Extra port on condenser to hook back up to pressure-enthalpy properties
    m.fs.condenser._outlet_1_enth_mol_ref = pyo.Reference(
        m.fs.condenser.shell.properties_out[:].enth_mol
    )
    m.fs.condenser.outlet_1_ph = Port(
        initialize={
            "flow_mol": m.fs.condenser._outlet_1_flow_mol_ref,
            "pressure": m.fs.condenser._outlet_1_pressure_ref,
            "enth_mol": m.fs.condenser._outlet_1_enth_mol_ref,
        }
    )
    # hotwell
    m.fs.hotwell = Mixer(
        default={
            "momentum_mixing_type": MomentumMixingType.none,
            "material_balance_type": MaterialBalanceType.componentTotal,
            "inlet_list": ["condensate", "makeup"],
            "property_package": m.fs.prop_water,
        }
    )

    @m.fs.hotwell.Constraint(m.fs.time)
    def mixer_pressure_constraint(b, t):
        return b.condensate_state[t].pressure == b.mixed_state[t].pressure

    # condensate pump
    m.fs.cond_pump = PressureChanger(
        default={
            "property_package": m.fs.prop_water,
            "material_balance_type": MaterialBalanceType.componentTotal,
            "thermodynamic_assumption": ThermodynamicAssumption.pump,
        }
    )
    ############################################################################
    #  Add low pressure feedwater heaters                                      #
    ############################################################################
    # Need to set the material balance type on all the feedwater heaters, so
    # define fwh_config to condense the code down a bit.
    fwh_config = {
        "delta_temperature_callback": delta_temperature_underwood_callback,
        "shell": {"material_balance_type": MaterialBalanceType.componentTotal},
        "tube": {"material_balance_type": MaterialBalanceType.componentTotal},
    }
    # fwh1
    m.fs.fwh1 = FWH0D(
        default={
            "has_desuperheat": False,
            "has_drain_cooling": False,
            "has_drain_mixer": True,
            "property_package": m.fs.prop_water,
            "condense": fwh_config,
        }
    )
    # pump for fwh1 condensate, to pump it ahead and mix with feedwater
    m.fs.fwh1_pump = PressureChanger(
        default={
            "property_package": m.fs.prop_water,
            "material_balance_type": MaterialBalanceType.componentTotal,
            "thermodynamic_assumption": ThermodynamicAssumption.pump,
        }
    )
    # Mix the FWH1 drain back into the feedwater
    m.fs.fwh1_return = Mixer(
        default={
            "momentum_mixing_type": MomentumMixingType.none,
            "material_balance_type": MaterialBalanceType.componentTotal,
            "inlet_list": ["feedwater", "fwh1_drain"],
            "property_package": m.fs.prop_water,
        }
    )

    @m.fs.fwh1_return.Constraint(m.fs.time)
    def mixer_pressure_constraint(b, t):
        return b.feedwater_state[t].pressure == b.mixed_state[t].pressure

    m.fs.fwh2 = FWH0D(
        default={
            "has_desuperheat": True,
            "has_drain_cooling": True,
            "has_drain_mixer": True,
            "property_package": m.fs.prop_water,
            "desuperheat": fwh_config,
            "cooling": fwh_config,
            "condense": fwh_config,
        }
    )
    m.fs.fwh3 = FWH0D(
        default={
            "has_desuperheat": True,
            "has_drain_cooling": True,
            "has_drain_mixer": True,
            "property_package": m.fs.prop_water,
            "desuperheat": fwh_config,
            "cooling": fwh_config,
            "condense": fwh_config,
        }
    )
    m.fs.fwh4 = FWH0D(
        default={
            "has_desuperheat": True,
            "has_drain_cooling": True,
            "has_drain_mixer": False,
            "property_package": m.fs.prop_water,
            "desuperheat": fwh_config,
            "cooling": fwh_config,
            "condense": fwh_config,
        }
    )
    ############################################################################
    #  Add deaerator and boiler feed pump (BFP)                                #
    ############################################################################
    m.fs.fwh5_da = Mixer(
        default={
            "momentum_mixing_type": MomentumMixingType.none,
            "material_balance_type": MaterialBalanceType.componentTotal,
            "inlet_list": ["steam", "drain", "feedwater"],
            "property_package": m.fs.prop_water,
        }
    )

    @m.fs.fwh5_da.Constraint(m.fs.time)
    def mixer_pressure_constraint(b, t):
        # Not sure about deaerator pressure, so assume same as feedwater inlet
        return b.feedwater_state[t].pressure == b.mixed_state[t].pressure

    m.fs.bfp = PressureChanger(
        default={
            "property_package": m.fs.prop_water,
            "material_balance_type": MaterialBalanceType.componentTotal,
            "thermodynamic_assumption": ThermodynamicAssumption.pump,
        }
    )
    m.fs.bfpt = PressureChanger(
        default={
            "property_package": m.fs.prop_water,
            "compressor": False,
            "material_balance_type": MaterialBalanceType.componentTotal,
            "thermodynamic_assumption": ThermodynamicAssumption.isentropic,
        }
    )

    @m.fs.Constraint(m.fs.time)
    def constraint_out_pressure(b, t):
        return (
            b.bfpt.control_volume.properties_out[t].pressure
            == b.condenser_mix.mixed_state[t].pressure
        )

    @m.fs.Constraint(m.fs.time)
    def constraint_bfp_power(b, t):
        return 0 == b.bfp.control_volume.work[t] + b.bfpt.control_volume.work[t]

    ############################################################################
    #  Add high pressure feedwater heaters                                     #
    ############################################################################
    m.fs.fwh6 = FWH0D(
        default={
            "has_desuperheat": True,
            "has_drain_cooling": True,
            "has_drain_mixer": True,
            "property_package": m.fs.prop_water,
            "desuperheat": fwh_config,
            "cooling": fwh_config,
            "condense": fwh_config,
        }
    )
    m.fs.fwh7 = FWH0D(
        default={
            "has_desuperheat": True,
            "has_drain_cooling": True,
            "has_drain_mixer": True,
            "property_package": m.fs.prop_water,
            "desuperheat": fwh_config,
            "cooling": fwh_config,
            "condense": fwh_config,
        }
    )
    m.fs.fwh8 = FWH0D(
        default={
            "has_desuperheat": True,
            "has_drain_cooling": True,
            "has_drain_mixer": False,
            "property_package": m.fs.prop_water,
            "desuperheat": fwh_config,
            "cooling": fwh_config,
            "condense": fwh_config,
        }
    )
    ############################################################################
    #  Additional Constraints/Expressions                                      #
    ############################################################################

    # Add a constraint to set the amount of pressure lost from the boiler feed
    # pump to the turbine inlet.  This isn't great, but it's an okay
    # approximation given that this particular model doesn't contain a boiler
    # model.
    m.fs.boiler_pressure_drop_fraction = pyo.Var(
        m.fs.time,
        initialize=0.01,
        doc="Fraction of pressure lost from boiler feed pump and turbine inlet",
    )

    @m.fs.Constraint(m.fs.time)
    def boiler_pressure_drop(b, t):
        return (
            m.fs.bfp.control_volume.properties_out[t].pressure
            * (1 - b.boiler_pressure_drop_fraction[t])
            == m.fs.turb.inlet_split.mixed_state[t].pressure
        )

    # Again since the boiler is missing, set the flow of steam into the turbine
    # equal to the flow of feedwater out of the last feedwater heater.
    @m.fs.Constraint(m.fs.time)
    def close_flow(b, t):
        return (
            m.fs.bfp.control_volume.properties_out[t].flow_mol
            == m.fs.turb.inlet_split.mixed_state[t].flow_mol
        )

    # Add some expressions to calculate efficiency.
    @m.fs.Expression(m.fs.time)
    def boiler_heat(b, t):
        return (
            b.turb.inlet_split.mixed_state[t].enth_mol
            * b.turb.inlet_split.mixed_state[t].flow_mol
            - b.fwh8.desuperheat.tube.properties_out[t].enth_mol
            * b.fwh8.desuperheat.tube.properties_out[t].flow_mol
            + b.turb.ip_stages[1].control_volume.properties_in[t].enth_mol
            * b.turb.ip_stages[1].control_volume.properties_in[t].flow_mol
            - b.turb.hp_split[7].outlet_1.enth_mol[t]
            * b.turb.hp_split[7].outlet_1.flow_mol[t]
        )

    @m.fs.Expression(m.fs.time)
    def steam_cycle_eff(b, t):
        return -100 * b.turb.power[t] / b.boiler_heat[t]

    ############################################################################
    #  Create the stream Arcs and return the model                             #
    ############################################################################
    _create_arcs(m)
    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)
    return m


def _create_arcs(m):
    ############################################################################
    #  Connect turbine and condenser units                                     #
    ############################################################################
    m.fs.turb_to_condenser = Arc(
        source=m.fs.turb.outlet_stage.outlet, destination=m.fs.condenser_mix.main
    )
    m.fs.condenser_mix_to_condenser = Arc(
        source=m.fs.condenser_mix.outlet_tpx, destination=m.fs.condenser.inlet_1
    )
    m.fs.condenser_to_hotwell = Arc(
        source=m.fs.condenser.outlet_1_ph, destination=m.fs.hotwell.condensate
    )
    m.fs.hotwell_to_cond_pump = Arc(
        source=m.fs.hotwell.outlet, destination=m.fs.cond_pump.inlet
    )
    ############################################################################
    #  Low pressure FWHs                                                       #
    ############################################################################
    # fwh1
    m.fs.turb_to_fwh1 = Arc(
        source=m.fs.turb.lp_split[11].outlet_2, destination=m.fs.fwh1.drain_mix.steam
    )
    m.fs.cond_pump_to_fwh1 = Arc(
        source=m.fs.cond_pump.outlet, destination=m.fs.fwh1.condense.inlet_2
    )
    m.fs.fwh1_drain_to_pump = Arc(
        source=m.fs.fwh1.condense.outlet_1, destination=m.fs.fwh1_pump.inlet
    )
    m.fs.fwh1_pump_to_return = Arc(
        source=m.fs.fwh1_pump.outlet, destination=m.fs.fwh1_return.fwh1_drain
    )
    m.fs.fwh1_to_retrun = Arc(
        source=m.fs.fwh1.condense.outlet_2, destination=m.fs.fwh1_return.feedwater
    )
    # fwh2
    m.fs.fwh1_return_to_fwh2 = Arc(
        source=m.fs.fwh1_return.outlet, destination=m.fs.fwh2.cooling.inlet_2
    )
    m.fs.fwh2_drain_to_fwh1 = Arc(
        source=m.fs.fwh2.cooling.outlet_1, destination=m.fs.fwh1.drain_mix.drain
    )
    m.fs.turb_to_fwh2 = Arc(
        source=m.fs.turb.lp_split[10].outlet_2,
        destination=m.fs.fwh2.desuperheat.inlet_1,
    )
    # fwh3
    m.fs.fwh2_to_fwh3 = Arc(
        source=m.fs.fwh2.desuperheat.outlet_2, destination=m.fs.fwh3.cooling.inlet_2
    )
    m.fs.fwh3_drain_to_fwh2 = Arc(
        source=m.fs.fwh3.cooling.outlet_1, destination=m.fs.fwh2.drain_mix.drain
    )
    m.fs.turb_to_fwh3 = Arc(
        source=m.fs.turb.lp_split[8].outlet_2, destination=m.fs.fwh3.desuperheat.inlet_1
    )
    # fwh4
    m.fs.fwh3_to_fwh4 = Arc(
        source=m.fs.fwh3.desuperheat.outlet_2, destination=m.fs.fwh4.cooling.inlet_2
    )
    m.fs.fwh4_drain_to_fwh3 = Arc(
        source=m.fs.fwh4.cooling.outlet_1, destination=m.fs.fwh3.drain_mix.drain
    )
    m.fs.turb_to_fwh4 = Arc(
        source=m.fs.turb.lp_split[4].outlet_2, destination=m.fs.fwh4.desuperheat.inlet_1
    )
    ############################################################################
    #  FWH5 (Deaerator) and boiler feed pump (BFP)                             #
    ############################################################################
    m.fs.fwh4_to_fwh5 = Arc(
        source=m.fs.fwh4.desuperheat.outlet_2, destination=m.fs.fwh5_da.feedwater
    )
    m.fs.turb_to_fwh5 = Arc(
        source=m.fs.turb.ip_split[10].outlet_2, destination=m.fs.fwh5_da.steam
    )
    m.fs.fwh5_to_bfp = Arc(source=m.fs.fwh5_da.outlet, destination=m.fs.bfp.inlet)
    m.fs.turb_to_bfpt = Arc(
        source=m.fs.turb.ip_split[10].outlet_3, destination=m.fs.bfpt.inlet
    )
    m.fs.bfpt_to_condenser = Arc(
        source=m.fs.bfpt.outlet, destination=m.fs.condenser_mix.bfpt
    )
    ############################################################################
    #  High pressure feedwater heaters                                         #
    ############################################################################
    # fwh6
    m.fs.bfp_to_fwh6 = Arc(
        source=m.fs.bfp.outlet, destination=m.fs.fwh6.cooling.inlet_2
    )
    m.fs.fwh6_drain_to_fwh5 = Arc(
        source=m.fs.fwh6.cooling.outlet_1, destination=m.fs.fwh5_da.drain
    )
    m.fs.turb_to_fwh6 = Arc(
        source=m.fs.turb.ip_split[5].outlet_2, destination=m.fs.fwh6.desuperheat.inlet_1
    )
    # fwh7
    m.fs.fwh6_to_fwh7 = Arc(
        source=m.fs.fwh6.desuperheat.outlet_2, destination=m.fs.fwh7.cooling.inlet_2
    )
    m.fs.fwh7_drain_to_fwh6 = Arc(
        source=m.fs.fwh7.cooling.outlet_1, destination=m.fs.fwh6.drain_mix.drain
    )
    m.fs.turb_to_fwh7 = Arc(
        source=m.fs.turb.hp_split[7].outlet_2, destination=m.fs.fwh7.desuperheat.inlet_1
    )
    # fwh8
    m.fs.fwh7_to_fwh8 = Arc(
        source=m.fs.fwh7.desuperheat.outlet_2, destination=m.fs.fwh8.cooling.inlet_2
    )
    m.fs.fwh8_drain_to_fwh7 = Arc(
        source=m.fs.fwh8.cooling.outlet_1, destination=m.fs.fwh7.drain_mix.drain
    )
    m.fs.turb_to_fwh8 = Arc(
        source=m.fs.turb.hp_split[4].outlet_2, destination=m.fs.fwh8.desuperheat.inlet_1
    )


def _stream_dict(m):
    """Adds _streams to m, which contains a dictionary of streams for display

    Args:
        m (ConcreteModel): A Pyomo model from create_model()

    Returns:
        None
    """

    m._streams = OrderedDict(
        [
            ("STEAM_MAIN", m.fs.turb.inlet_split.mixed_state),
            ("THRTL1", m.fs.turb.inlet_stage[1].control_volume.properties_in),
            ("THRTL2", m.fs.turb.inlet_stage[2].control_volume.properties_in),
            ("THRTL3", m.fs.turb.inlet_stage[3].control_volume.properties_in),
            ("THRTL4", m.fs.turb.inlet_stage[4].control_volume.properties_in),
            ("RHT_COLD", m.fs.turb.hp_split[7].outlet_1_state),
            ("RHT_HOT", m.fs.turb.ip_stages[1].control_volume.properties_in),
            ("STEAM_LP", m.fs.turb.ip_split[10].outlet_1_state),
            ("EXTR_HP4", m.fs.turb_to_fwh8),
            ("EXTR_HP7", m.fs.turb_to_fwh7),
            ("EXTR_IP5", m.fs.turb_to_fwh6),
            ("EXTR_IP10", m.fs.turb_to_fwh5),
            ("EXTR_BFPT_A", m.fs.turb.ip_split[10].outlet_3_state),
            ("EXTR_LP4", m.fs.turb_to_fwh4),
            ("EXTR_LP8", m.fs.turb_to_fwh3),
            ("EXTR_LP10", m.fs.turb_to_fwh2),
            ("EXTR_LP11", m.fs.turb_to_fwh1),
            ("EXHST_MAIN", m.fs.turb_to_condenser),
            ("EXHST_BFPT", m.fs.bfpt_to_condenser),
            ("CW01", m.fs.condenser.tube.properties_in),
            ("CW02", m.fs.condenser.tube.properties_out),
            ("MAKEUP_01", m.fs.hotwell.makeup_state),
            ("COND_01", m.fs.condenser_to_hotwell),
            ("COND_02", m.fs.hotwell_to_cond_pump),
            ("COND_03", m.fs.cond_pump_to_fwh1),
            ("FW01A", m.fs.fwh1_to_retrun),
            ("FW01B", m.fs.fwh1_return_to_fwh2),
            ("FW02", m.fs.fwh2_to_fwh3),
            ("FW03", m.fs.fwh3_to_fwh4),
            ("FW04", m.fs.fwh4_to_fwh5),
            ("FW05A", m.fs.fwh5_to_bfp),
            ("FW05B", m.fs.bfp_to_fwh6),
            ("FW06", m.fs.fwh6_to_fwh7),
            ("FW07", m.fs.fwh7_to_fwh8),
            ("FW08", m.fs.fwh8.desuperheat.tube.properties_out),
            ("FWH1_DRN1", m.fs.fwh1_drain_to_pump),
            ("FWH1_DRN2", m.fs.fwh1_pump_to_return),
            ("FWH2_DRN", m.fs.fwh2_drain_to_fwh1),
            ("FWH3_DRN", m.fs.fwh3_drain_to_fwh2),
            ("FWH4_DRN", m.fs.fwh4_drain_to_fwh3),
            ("FWH6_DRN", m.fs.fwh6_drain_to_fwh5),
            ("FWH7_DRN", m.fs.fwh7_drain_to_fwh6),
            ("FWH8_DRN", m.fs.fwh8_drain_to_fwh7),
        ]
    )


def set_model_input(m):
    """ Fix some variables and set values. Generally get the model ready to run
    in simulation mode (0 degrees of freedom).

    Args:
        m (ConcreteModel): A Pyomo model from create_model()

    Returns:
        None
    """
    ############################################################################
    #  Turbine input                                                           #
    ############################################################################
    main_steam_pressure = 2.423e7
    m.fs.turb.turbine_inlet_cf_fix(2.7e-4)
    m.fs.turb.turbine_outlet_cf_fix(0.15)
    # Set the turbine steam inlet conditions and flow guess for init
    m.fs.turb.inlet_split.inlet.enth_mol.fix(62710)
    m.fs.turb.inlet_split.inlet.pressure.fix(main_steam_pressure)
    m.fs.boiler_pressure_drop_fraction.fix(0.1)
    m.fs.turb.inlet_split.inlet.flow_mol[:].value = 26000
    m.fs.turb.inlet_split.inlet.flow_mol.unfix()  # Pressure driven
    m.fs.turb.inlet_mix.use_equal_pressure_constraint()  # Pressure driven mix
    # Set throttle valve
    m.fs.turb.throttle_cv_fix(0.0009)
    m.fs.turb.throttle_valve[1].valve_opening.fix(0.9)
    m.fs.turb.throttle_valve[2].valve_opening.fix(0.9)
    m.fs.turb.throttle_valve[3].valve_opening.fix(0.9)
    m.fs.turb.throttle_valve[4].valve_opening.fix(0.9)
    # Set the efficiency and pressure ratios of stages other than inlet and outlet
    for i, s in m.fs.turb.hp_stages.items():
        s.ratioP.fix(0.80)
        s.efficiency_isentropic.fix(0.9)
    for i, s in m.fs.turb.ip_stages.items():
        s.ratioP.fix(0.78)
        s.efficiency_isentropic.fix(0.9)
    for i, s in m.fs.turb.lp_stages.items():
        s.ratioP[:].fix(0.76)
        s.efficiency_isentropic[:].fix(0.9)
    ############################################################################
    #  Condenser section inputs                                                #
    ############################################################################
    m.fs.condenser.inlet_2.flow_mol.fix(2000000)
    m.fs.condenser.inlet_2.enth_mol.fix(1600)
    m.fs.condenser.inlet_2.pressure.fix(500000)
    m.fs.condenser.delta_temperature_out.fix(3)
    m.fs.condenser.area.fix(2.2e4)
    m.fs.condenser.area.unfix()  # to make solving easier specifying approach
    m.fs.condenser.overall_heat_transfer_coefficient.fix(2000)
    m.fs.hotwell.makeup.flow_mol[:].value = 1
    m.fs.hotwell.makeup.enth_mol.fix(2500)
    m.fs.hotwell.makeup.pressure.fix(101325)
    m.fs.cond_pump.efficiency_pump.fix(0.80)
    m.fs.cond_pump.deltaP.fix(1e6)
    ############################################################################
    #  Low pressure FWH section inputs                                         #
    ############################################################################
    # fwh1
    m.fs.fwh1.condense.area.fix(400)
    m.fs.fwh1.condense.overall_heat_transfer_coefficient.fix(2000)
    # fwh1 pump
    m.fs.fwh1_pump.efficiency_pump.fix(0.80)
    m.fs.fwh1_pump.deltaP.fix(1.2e6)
    # fwh2
    m.fs.fwh2.condense.area.fix(150)
    m.fs.fwh2.condense.overall_heat_transfer_coefficient.fix(2000)
    m.fs.fwh2.desuperheat.area.fix(75)
    m.fs.fwh2.desuperheat.overall_heat_transfer_coefficient.fix(600)
    m.fs.fwh2.cooling.area.fix(75)
    m.fs.fwh2.cooling.overall_heat_transfer_coefficient.fix(300)
    # fwh3
    m.fs.fwh3.condense.area.fix(100)
    m.fs.fwh3.condense.overall_heat_transfer_coefficient.fix(2000)
    m.fs.fwh3.desuperheat.area.fix(50)
    m.fs.fwh3.desuperheat.overall_heat_transfer_coefficient.fix(600)
    m.fs.fwh3.cooling.area.fix(50)
    m.fs.fwh3.cooling.overall_heat_transfer_coefficient.fix(300)
    # fwh4
    m.fs.fwh4.condense.area.fix(100)
    m.fs.fwh4.condense.overall_heat_transfer_coefficient.fix(2000)
    m.fs.fwh4.desuperheat.area.fix(50)
    m.fs.fwh4.desuperheat.overall_heat_transfer_coefficient.fix(600)
    m.fs.fwh4.cooling.area.fix(50)
    m.fs.fwh4.cooling.overall_heat_transfer_coefficient.fix(300)
    ############################################################################
    #  Deaerator and boiler feed pump (BFP) Input                              #
    ############################################################################
    m.fs.bfp.efficiency_pump.fix(0.80)
    m.fs.bfp.outlet.pressure[:].value = main_steam_pressure * 1.1
    m.fs.bfpt.efficiency_isentropic.fix(0.80)
    # m.fs.bfpt.ratioP.fix(0.5)
    ############################################################################
    #  High pressure feedwater heater                                          #
    ############################################################################
    # fwh6
    m.fs.fwh6.condense.area.fix(300)
    m.fs.fwh6.condense.overall_heat_transfer_coefficient.fix(2000)
    m.fs.fwh6.desuperheat.area.fix(150)
    m.fs.fwh6.desuperheat.overall_heat_transfer_coefficient.fix(600)
    m.fs.fwh6.cooling.area.fix(150)
    m.fs.fwh6.cooling.overall_heat_transfer_coefficient.fix(300)
    # fwh7
    m.fs.fwh7.condense.area.fix(200)
    m.fs.fwh7.condense.overall_heat_transfer_coefficient.fix(2000)
    m.fs.fwh7.desuperheat.area.fix(100)
    m.fs.fwh7.desuperheat.overall_heat_transfer_coefficient.fix(600)
    m.fs.fwh7.cooling.area.fix(100)
    m.fs.fwh7.cooling.overall_heat_transfer_coefficient.fix(300)
    # fwh8
    m.fs.fwh8.condense.area.fix(200)
    m.fs.fwh8.condense.overall_heat_transfer_coefficient.fix(2000)
    m.fs.fwh8.desuperheat.area.fix(100)
    m.fs.fwh8.desuperheat.overall_heat_transfer_coefficient.fix(600)
    m.fs.fwh8.cooling.area.fix(100)
    m.fs.fwh8.cooling.overall_heat_transfer_coefficient.fix(300)


def initialize(m, fileoutput=None, fileinput=None):
    """ Initialize a mode from create_model(), set model inputs before
    initializing.

    Args:
        m (ConcreteModel): A Pyomo model from create_model()
        fileoutput (str|None): File path to save initialized model state. If
            None, don't save anything
        fileinput (str|None): File to load initialized model state from. If a
            file is supplied skip initialization routine. If None, initialize.

    Returns:
        solver: A Pyomo solver object, that can be used to solve the model.

    """
    if len(m.fs.time) > 1:
        m.fs.condenser.heat_transfer_equation.deactivate()
    solver = pyo.SolverFactory("ipopt")
    solver.options = {
        "tol": 1e-6,
        #'linear_solver': "mumps",
        "max_iter": 300,
    }
    if fileinput is not None:
        ms.from_json(m, fname=fileinput)
        return solver
    ############################################################################
    #  Initialize turbine                                                      #
    ############################################################################
    # Extraction rates are calculated from the feedwater heater models, so to
    # initialize the tubine fix some initial guesses. They get unfixed after
    # solving the turbine
    m.fs.turb.outlet_stage.control_volume.properties_out[:].pressure.fix(9000)
    m.fs.turb.lp_split[11].split_fraction[:, "outlet_2"].fix(0.04403)
    m.fs.turb.lp_split[10].split_fraction[:, "outlet_2"].fix(0.04025)
    m.fs.turb.lp_split[8].split_fraction[:, "outlet_2"].fix(0.04362)
    m.fs.turb.lp_split[4].split_fraction[:, "outlet_2"].fix(0.08025)
    m.fs.turb.ip_split[10].split_fraction[:, "outlet_2"].fix(0.045)
    m.fs.turb.ip_split[10].split_fraction[:, "outlet_3"].fix(0.04)
    m.fs.turb.ip_split[5].split_fraction[:, "outlet_2"].fix(0.05557)
    m.fs.turb.hp_split[7].split_fraction[:, "outlet_2"].fix(0.09741)
    m.fs.turb.hp_split[4].split_fraction[:, "outlet_2"].fix(0.0740)
    # Put in a rough initial guess for the ip section inlet, since it is
    # disconnected from the hp section for the reheater.
    ip1_pin = 5.35e6
    ip1_hin = iapws95.htpx(T=866, P=ip1_pin)
    ip1_fin = pyo.value(m.fs.turb.inlet_split.inlet.flow_mol[0])
    m.fs.turb.ip_stages[1].inlet.enth_mol[:].value = ip1_hin
    m.fs.turb.ip_stages[1].inlet.flow_mol[:].value = ip1_fin
    m.fs.turb.ip_stages[1].inlet.pressure[:].value = ip1_pin
    # initialize turbine
    m.fs.turb.initialize(outlvl=1)
    # solve with pressure driven flow to get a bit better solution before
    # initializing rest of units
    assert degrees_of_freedom(m.fs.turb) == 0
    res = solver.solve(m.fs.turb, tee=True)
    m.fs.turb.outlet_stage.control_volume.properties_out[:].pressure.unfix()
    # Extraction rates are calculated once the feedwater heater models are
    # add so unfix the splits.
    m.fs.turb.lp_split[11].split_fraction[:, "outlet_2"].unfix()
    m.fs.turb.lp_split[10].split_fraction[:, "outlet_2"].unfix()
    m.fs.turb.lp_split[8].split_fraction[:, "outlet_2"].unfix()
    m.fs.turb.lp_split[4].split_fraction[:, "outlet_2"].unfix()
    m.fs.turb.ip_split[10].split_fraction[:, "outlet_3"].unfix()
    m.fs.turb.ip_split[5].split_fraction[:, "outlet_2"].unfix()
    m.fs.turb.hp_split[7].split_fraction[:, "outlet_2"].unfix()
    m.fs.turb.hp_split[4].split_fraction[:, "outlet_2"].unfix()
    # initialize the boiler feed pump turbine.
    _set_port(m.fs.bfpt.inlet, m.fs.turb.ip_split[10].outlet_3)
    m.fs.bfpt.ratioP.fix(0.4)
    m.fs.bfpt.initialize(outlvl=1)
    # m.fs.bfpt.ratioP.unfix()
    m.fs.constraint_out_pressure.deactivate()
    ############################################################################
    #  Condenser section                                                       #
    ############################################################################
    # initialize condenser mixer
    _set_port(m.fs.condenser_mix.main, m.fs.turb.outlet_stage.outlet)
    _set_port(m.fs.condenser_mix.bfpt, m.fs.bfpt.outlet)
    m.fs.condenser_mix.initialize(outlvl=1)
    # initialize condenser hx
    _set_port(m.fs.condenser.inlet_1, m.fs.condenser_mix.outlet_tpx)
    _set_port(m.fs.condenser.outlet_2, m.fs.condenser.inlet_2)
    # This still has the outlet vapor fraction fixed at 0, so if that isn't
    # true for the initial pressure guess this could fail to initialize
    m.fs.condenser.inlet_1.fix()
    m.fs.condenser.inlet_1.pressure.unfix()
    solver.solve(m.fs.condenser, tee=True)
    m.fs.condenser.inlet_1.unfix()
    # m.fs.condenser.eq_pressure.deactivate()
    # m.fs.condenser.initialize(outlvl=1)
    # m.fs.condenser.eq_pressure.activate()
    # initialize hotwell
    _set_port(m.fs.hotwell.condensate, m.fs.condenser.outlet_1_ph)
    m.fs.hotwell.initialize(outlvl=1)
    m.fs.hotwell.condensate.unfix()
    # initialize condensate pump
    _set_port(m.fs.cond_pump.inlet, m.fs.hotwell.outlet)
    m.fs.cond_pump.initialize(outlvl=1)
    ############################################################################
    #  Low pressure FWH section                                                #
    ############################################################################
    # fwh1
    m.fs.fwh1.drain_mix.drain.flow_mol[:] = 1000
    m.fs.fwh1.drain_mix.drain.pressure[:] = 1e5
    m.fs.fwh1.drain_mix.drain.enth_mol[:] = 6117
    _set_port(m.fs.fwh1.condense.inlet_2, m.fs.cond_pump.outlet)
    _set_port(m.fs.fwh1.drain_mix.steam, m.fs.turb.lp_split[11].outlet_2)
    m.fs.fwh1.initialize(outlvl=1)
    # initialize fwh1 drain pump
    _set_port(m.fs.fwh1_pump.inlet, m.fs.fwh1.condense.outlet_1)
    m.fs.fwh1_pump.initialize(outlvl=1)
    # initialize mixer to add fwh1 drain to feedwater
    _set_port(m.fs.fwh1_return.feedwater, m.fs.fwh1.condense.outlet_2)
    _set_port(m.fs.fwh1_return.fwh1_drain, m.fs.fwh1.condense.outlet_1)
    m.fs.fwh1_return.initialize(outlvl=1)
    m.fs.fwh1_return.feedwater.unfix()
    m.fs.fwh1_return.fwh1_drain.unfix()
    # fwh2
    m.fs.fwh2.drain_mix.drain.flow_mol[:] = 100
    m.fs.fwh2.drain_mix.drain.pressure[:] = 1.5e5
    m.fs.fwh2.drain_mix.drain.enth_mol[:] = 7000
    _set_port(m.fs.fwh2.cooling.inlet_2, m.fs.fwh1_return.outlet)
    _set_port(m.fs.fwh2.desuperheat.inlet_1, m.fs.turb.lp_split[10].outlet_2)
    m.fs.fwh2.initialize(outlvl=1)
    # fwh3
    m.fs.fwh3.drain_mix.drain.flow_mol[:] = 100
    m.fs.fwh3.drain_mix.drain.pressure[:] = 2.5e5
    m.fs.fwh3.drain_mix.drain.enth_mol[:] = 8000
    _set_port(m.fs.fwh3.cooling.inlet_2, m.fs.fwh2.desuperheat.outlet_2)
    _set_port(m.fs.fwh3.desuperheat.inlet_1, m.fs.turb.lp_split[8].outlet_2)
    m.fs.fwh3.initialize(outlvl=1)
    # fwh4
    _set_port(m.fs.fwh4.cooling.inlet_2, m.fs.fwh3.desuperheat.outlet_2)
    _set_port(m.fs.fwh4.desuperheat.inlet_1, m.fs.turb.lp_split[4].outlet_2)
    m.fs.fwh4.initialize(outlvl=1)
    ############################################################################
    #  boiler feed pump and deaerator                                          #
    ############################################################################
    _set_port(m.fs.fwh5_da.feedwater, m.fs.fwh4.desuperheat.outlet_2)
    _set_port(m.fs.fwh5_da.steam, m.fs.turb.ip_split[10].outlet_2)
    m.fs.fwh5_da.drain.flow_mol[:] = 2000
    m.fs.fwh5_da.drain.pressure[:] = 3e6
    m.fs.fwh5_da.drain.enth_mol[:] = 9000
    m.fs.fwh5_da.initialize(outlvl=1)
    _set_port(m.fs.bfp.inlet, m.fs.fwh5_da.outlet)
    m.fs.bfp.control_volume.properties_out[:].pressure.fix()
    m.fs.bfp.initialize(outlvl=1)
    m.fs.bfp.control_volume.properties_out[:].pressure.unfix()
    ############################################################################
    #  High-pressure feedwater heaters                                         #
    ############################################################################
    # fwh6
    m.fs.fwh6.drain_mix.drain.flow_mol[:] = 1000
    m.fs.fwh6.drain_mix.drain.pressure[:] = 1e7
    m.fs.fwh6.drain_mix.drain.enth_mol[:] = 9500
    _set_port(m.fs.fwh6.cooling.inlet_2, m.fs.bfp.outlet)
    _set_port(m.fs.fwh6.desuperheat.inlet_1, m.fs.turb.ip_split[5].outlet_2)
    m.fs.fwh6.initialize(outlvl=1)
    # fwh7
    m.fs.fwh7.drain_mix.drain.flow_mol[:] = 2000
    m.fs.fwh7.drain_mix.drain.pressure[:] = 1e7
    m.fs.fwh7.drain_mix.drain.enth_mol[:] = 9500
    _set_port(m.fs.fwh7.cooling.inlet_2, m.fs.fwh6.desuperheat.outlet_2)
    _set_port(m.fs.fwh7.desuperheat.inlet_1, m.fs.turb.hp_split[7].outlet_2)
    m.fs.fwh7.initialize(outlvl=1)
    # fwh8
    _set_port(m.fs.fwh8.cooling.inlet_2, m.fs.fwh7.desuperheat.outlet_2)
    _set_port(m.fs.fwh8.desuperheat.inlet_1, m.fs.turb.hp_split[4].outlet_2)
    m.fs.fwh8.initialize(outlvl=1)

    ############################################################################
    #  Save and return solver                                                  #
    ############################################################################
    if fileoutput is not None:
        raise NotImplementedError("Saving initialization not finished")
    return solver


def pfd_result(outfile, m, df):
    tags = {}
    for i in df.index:
        tags[i + "_F"] = df.loc[i, "Molar Flow (mol/s)"]
        tags[i + "_T"] = df.loc[i, "T (K)"]
        tags[i + "_P"] = df.loc[i, "P (Pa)"]
        tags[i + "_X"] = df.loc[i, "Vapor Fraction"]
    tags["gross_power"] = -pyo.value(m.fs.turb.power[0])
    tags["gross_power_mw"] = -pyo.value(m.fs.turb.power[0]) * 1e-6
    tags["steam_mass_flow"] = df.loc["STEAM_MAIN", "Mass Flow (kg/s)"]
    tags["sc_eff"] = pyo.value(m.fs.steam_cycle_eff[0])
    tags["boiler_heat"] = pyo.value(m.fs.boiler_heat[0]) * 1e-6
    tags["steam_pressure"] = df.loc["STEAM_MAIN", "P (Pa)"] / 1000.0
    tags["cond_pressure"] = df.loc["EXHST_MAIN", "P (Pa)"] / 1000.0
    tags["bfp_power"] = pyo.value(m.fs.bfp.work_mechanical[0])
    tags["bfp_eff"] = pyo.value(m.fs.bfp.efficiency_pump[0]) * 100
    tags["bfpt_power"] = pyo.value(m.fs.bfpt.work_mechanical[0])
    tags["bfpt_eff"] = pyo.value(m.fs.bfpt.efficiency_isentropic[0]) * 100
    original_svg_file = os.path.join(this_file_dir(), "supercritical_steam_cycle.svg")
    with open(original_svg_file, "r") as f:
        s = svg_tag(tags, f, outfile=outfile)


def main(initialize_from_file=None, store_initialization=None):
    """ Create and initalize a model and solver

    Args:
        None

    Returns:
        A tuple of a model and solver
    """
    m = create_model()
    _stream_dict(m)
    set_model_input(m)
    if initialize_from_file is None:
        solver = initialize(m)
    else:
        solver = initialize(m, fileinput=initialize_from_file)
    solver.solve(m, tee=True)
    if store_initialization is not None:
        ms.to_json(m, fname=store_initialization)
    return m, solver


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--initialize_from_file", help="File from which to load"
                       " initialized values. If specified, the initialization"
                       " proceedure will be skipped.", default=None)
    parser.add_argument("--store_initialization", help="If specified, store"
                        " the initialized model values, to reload later.",
                        default=None)
    args = parser.parse_args()
    m, solver = main(initialize_from_file=args.initialize_from_file,
                     store_initialization=args.store_initialization)
    df = create_stream_table_dataframe(streams=m._streams, orient="index")
    pfd_result("result.svg", m, df)
