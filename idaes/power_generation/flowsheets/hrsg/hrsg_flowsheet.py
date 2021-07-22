##############################################################################
# EMRE/DAC CRADA, NETL Collaboration agreement
# ToDo: complete section
##############################################################################
##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
NGCC HRSG Subsystem for a 690MWe plant
Three Pressure level heat exchanger network:

Water/Steam Route:
* LP Steam: FWH - LP ECON, ECON out + Nat Gas Pre heat - LP EVAP - LP SH
* IP Steam: LP EVAP to IP pump -> IP ECON1 -> IP Splitter1 -> IP ECON 2 -> IP EVAP -> IP SH1 -> IP SH1 + Cold reheat to IP Mixer -> IP SH2/RH1 -> IP SH3/RH2
                                                           -> Nat Gas Preheat
* HP Steam: LP EVAP to HP pump -> HP ECON1 -> ... -> HP ECON5 -> HP EVAP -> HP SH1 -> HP SH2 -> HP SH3 -> HP SH4


Flue Gas Route:
Nat Gas Turbine -> HP SH4 -> IP SH3 -> HP SH3 -> HP SH2 -> IP SH2 -> HP SH1 ->*
* HP EVAP -> HP ECON5 -> IP SH1 -> HP ECON4 -> HP ECON3 -> LP SH -> IP EVAP ->*
* IP ECON2 -> HP ECON 2 -> IP ECON1 -> HP ECON1 -> LP EVAP -> LP ECON

ECON - Economizer
SH - Super Heater
EVAP - Evaporator
FWH - Feed Water Heater
"""

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.network import Arc

# Import IDAES core
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util import copy_port_values as _set_port
from idaes.core import FlowsheetBlock
import idaes.logger as idaeslog

# Import IDAES standard unit model
from idaes.generic_models.unit_models import Mixer
from idaes.generic_models.properties import iapws95
from idaes.power_generation.properties import FlueGasParameterBlock

# Import IDAES power generation unit models
from idaes.power_generation.unit_models.helm import (
    HelmMixer,
    MomentumMixingType,
    HelmSplitter,
    HelmIsentropicCompressor as WaterPump
    )
from idaes.generic_models.unit_models.heat_exchanger import (
    HeatExchanger,
    HeatExchangerFlowPattern,
    delta_temperature_underwood_callback)
import idaes.core.util.scaling as iscale
from idaes.power_generation.unit_models.boiler_heat_exchanger import (
    BoilerHeatExchanger, TubeArrangement, DeltaTMethod)
from idaes.power_generation.unit_models.helm.phase_separator import \
    HelmPhaseSeparator

from idaes.core.util.misc import svg_tag
import os
from pyomo.common.fileutils import this_file_dir
import idaes.core.util.tables as ta
from idaes.core.util import model_serializer as ms
from idaes.generic_models.unit_models import Separator as Splitter

__author__ = "M. Zamarripa"


def add_unit_models(m):
    """
    Function to add unit operation models to the flowsheet
    """
    # Models will be added to the boiler sub-flowsheet of the main flowsheet
    fs = m.fs

    # Water and ideal gas properties defined on main flowsheet
    prop_water = m.fs.prop_water

    # -------------------------------------------------------------------------
    # Low Pressure System: 1 economizer, 1 evaporator, 1 superheater
    fs.LP_ECON = BoilerHeatExchanger(default={
        "side_1_property_package": fs.prop_water,
        "side_2_property_package": fs.prop_gas,
        "has_pressure_change": False,
        "has_holdup": False,
        "delta_T_method": DeltaTMethod.counterCurrent,
        "tube_arrangement": TubeArrangement.inLine,
        "side_1_water_phase": "Liq",
        "has_radiation": False,
        "underwood":True,})

    # Mixer for LP_ECON outlet and Preheater streams
    fs.Mixer1 = HelmMixer(
        default={"dynamic": False,
                 "property_package": prop_water,
                 "momentum_mixing_type": MomentumMixingType.minimize,
                 "inlet_list": ["LP_ECON", "Preheater"]})

    # LP drum - flash separator
    m.fs.LP_DRUM = HelmPhaseSeparator(default={'property_package':
                                               m.fs.prop_water})

    fs.LP_EVAP = HeatExchanger(default={
        "shell": {"property_package": m.fs.prop_gas},
        "tube": {"property_package": m.fs.prop_water},
        "delta_temperature_callback":delta_temperature_underwood_callback,
        "flow_pattern": HeatExchangerFlowPattern.countercurrent})

    m.fs.LP_FGsplit = Splitter(default={"property_package": m.fs.prop_gas,
                               "ideal_separation": False,
                               "outlet_list": ["toLP_SH", "toMixer"]})

    m.fs.LP_Mixer2 = Mixer(default={"property_package": m.fs.prop_gas,
                           "inlet_list": ["fromLP_SH", "bypass"]})

    # LP SH = Low pressure superheater
    fs.LP_SH = BoilerHeatExchanger(default={
        "side_1_property_package": fs.prop_water,
        "side_2_property_package": fs.prop_gas,
        "has_pressure_change": False,
        "has_holdup": False,
        "delta_T_method": DeltaTMethod.counterCurrent,
        "tube_arrangement": TubeArrangement.inLine,
        "side_1_water_phase": "Vap",
        "has_radiation": False,
        "underwood":True,})

    # LP DRUM LIQUID Outlet - Liquid Stream to IP and HP pumps
    m.fs.Splitter1 = HelmSplitter(
        default={
            "property_package": m.fs.prop_water,
            "outlet_list": ["toIP", "toHP"]
        }
    )

    # Intermediate Pressure Pump
    m.fs.IP_pump = WaterPump(
        default={
            "property_package": m.fs.prop_water,
        }
    )

    # High Pressure Pump
    m.fs.HP_pump = WaterPump(
        default={
            "property_package": m.fs.prop_water,
        }
    )
    # -------------------------------------------------------------------------
    # IP Section: IP economizer, IP Drum, IP Evaporator, IP Superheater 1 and 2
    # IP ECON
    # Low Pressure System: 1 economizer, 1 evaporator, 1 superheater
    fs.IP_ECON1 = BoilerHeatExchanger(default={
        "side_1_property_package": fs.prop_water,
        "side_2_property_package": fs.prop_gas,
        "has_pressure_change": True,
        "has_holdup": False,
        "delta_T_method": DeltaTMethod.counterCurrent,
        "tube_arrangement": TubeArrangement.inLine,
        "side_1_water_phase": "Liq",
        "has_radiation": False,
        "underwood":True,})

    # IP ECON Splitter 1, IP ECON outlet to IP ECON 2 and NG preheater (NGPH)
    m.fs.IP_Splitter1 = HelmSplitter(
        default={
            "property_package": m.fs.prop_water,
            "outlet_list": ["toIP_ECON2", "toNGPH"]
        }
    )
    # IP ECON2 (from splitter)
    fs.IP_ECON2 = BoilerHeatExchanger(default={
        "side_1_property_package": fs.prop_water,
        "side_2_property_package": fs.prop_gas,
        "has_pressure_change": True,
        "has_holdup": False,
        "delta_T_method": DeltaTMethod.counterCurrent,
        "tube_arrangement": TubeArrangement.inLine,
        "side_1_water_phase": "Liq",
        "has_radiation": False,
        "underwood":True,})

    fs.IP_EVAP = HeatExchanger(default={
        "shell": {"property_package": m.fs.prop_gas},
        "tube": {"property_package": m.fs.prop_water,
                 "has_pressure_change": True},
        "delta_temperature_callback":delta_temperature_underwood_callback,
        "flow_pattern": HeatExchangerFlowPattern.countercurrent})

    # IP SH = Low pressure superheater
    fs.IP_SH1 = BoilerHeatExchanger(default={
        "side_1_property_package": fs.prop_water,
        "side_2_property_package": fs.prop_gas,
        "has_pressure_change": True,
        "has_holdup": False,
        "delta_T_method": DeltaTMethod.counterCurrent,
        "tube_arrangement": TubeArrangement.inLine,
        "side_1_water_phase": "Vap",
        "has_radiation": False,
        "underwood":True,})

    # Mixer for LP_ECON outlet and Preheater streams
    fs.IP_Mixer1 = HelmMixer(
        default={"dynamic": False,
                 "property_package": prop_water,
                 "momentum_mixing_type": MomentumMixingType.minimize,
                 "inlet_list": ["IP_SH1", "Cold_reheat"]})

    # IP Splitter 2, HP Exit: Steam air ejector, attemp1, attemp 2
    m.fs.IP_Splitter2 = HelmSplitter(
        default={
            "property_package": m.fs.prop_water,
            "outlet_list": ["Cold_reheat", "toEjector", "toReclaimer",
                            "toDryer"]
        }
    )

    # IP SH2 = Intermediate pressure superheater2 and reheater 1
    fs.IP_SH2 = BoilerHeatExchanger(default={
        "side_1_property_package": fs.prop_water,
        "side_2_property_package": fs.prop_gas,
        "has_pressure_change": True,
        "has_holdup": False,
        "delta_T_method": DeltaTMethod.counterCurrent,
        "tube_arrangement": TubeArrangement.inLine,
        "side_1_water_phase": "Vap",
        "has_radiation": False,
        "underwood":True,})

    # IP SH3 = Intermediate pressure superheater3 and reheater 2
    fs.IP_SH3 = BoilerHeatExchanger(default={
        "side_1_property_package": fs.prop_water,
        "side_2_property_package": fs.prop_gas,
        "has_pressure_change": True,
        "has_holdup": False,
        "delta_T_method": DeltaTMethod.counterCurrent,
        "tube_arrangement": TubeArrangement.inLine,
        "side_1_water_phase": "Vap",
        "has_radiation": False,
        "underwood":True,})
    # -------------------------------------------------------------------------
    # High Pressure System ====================================================
    # inlet from # High Pressure Pump (m.fs.HP_pump.outlet)
    # HP_ECON1 = ECONOMIZER 1
    fs.HP_ECON1 = BoilerHeatExchanger(default={
        "side_1_property_package": fs.prop_water,
        "side_2_property_package": fs.prop_gas,
        "has_pressure_change": True,
        "has_holdup": False,
        "delta_T_method": DeltaTMethod.counterCurrent,
        "tube_arrangement": TubeArrangement.inLine,
        "side_1_water_phase": "Liq",
        "has_radiation": False,
        "underwood":True,})

    # HP_ECON2 = ECONOMIZER 2
    fs.HP_ECON2 = BoilerHeatExchanger(default={
        "side_1_property_package": fs.prop_water,
        "side_2_property_package": fs.prop_gas,
        "has_pressure_change": True,
        "has_holdup": False,
        "delta_T_method": DeltaTMethod.counterCurrent,
        "tube_arrangement": TubeArrangement.inLine,
        "side_1_water_phase": "Liq",
        "has_radiation": False,
        "underwood":True,})

    # HP_ECON3 = ECONOMIZER 3
    fs.HP_ECON3 = BoilerHeatExchanger(default={
        "side_1_property_package": fs.prop_water,
        "side_2_property_package": fs.prop_gas,
        "has_pressure_change": True,
        "has_holdup": False,
        "delta_T_method": DeltaTMethod.counterCurrent,
        "tube_arrangement": TubeArrangement.inLine,
        "side_1_water_phase": "Liq",
        "has_radiation": False,
        "underwood":True,})

    # HP_ECON4 = ECONOMIZER 4
    fs.HP_ECON4 = BoilerHeatExchanger(default={
        "side_1_property_package": fs.prop_water,
        "side_2_property_package": fs.prop_gas,
        "has_pressure_change": True,
        "has_holdup": False,
        "delta_T_method": DeltaTMethod.counterCurrent,
        "tube_arrangement": TubeArrangement.inLine,
        "side_1_water_phase": "Liq",
        "has_radiation": False,
        "underwood":True,})

    # HP_ECON5 = ECONOMIZER 5
    fs.HP_ECON5 = BoilerHeatExchanger(default={
        "side_1_property_package": fs.prop_water,
        "side_2_property_package": fs.prop_gas,
        "has_pressure_change": True,
        "has_holdup": False,
        "delta_T_method": DeltaTMethod.counterCurrent,
        "tube_arrangement": TubeArrangement.inLine,
        "side_1_water_phase": "Liq",
        "has_radiation": False,
        "underwood":True,})

    fs.HP_EVAP = HeatExchanger(default={
        "shell": {"property_package": m.fs.prop_gas},
        "tube": {"property_package": m.fs.prop_water},
        "delta_temperature_callback":delta_temperature_underwood_callback,
        "flow_pattern": HeatExchangerFlowPattern.countercurrent})

    # HP_SH1 = superheater 1
    fs.HP_SH1 = BoilerHeatExchanger(default={
        "side_1_property_package": fs.prop_water,
        "side_2_property_package": fs.prop_gas,
        "has_pressure_change": True,
        "has_holdup": False,
        "delta_T_method": DeltaTMethod.counterCurrent,
        "tube_arrangement": TubeArrangement.inLine,
        "side_1_water_phase": "Vap",
        "has_radiation": False,
        "underwood":True,})

    # HP_SH2 = superheater 2
    fs.HP_SH2 = BoilerHeatExchanger(default={
        "side_1_property_package": fs.prop_water,
        "side_2_property_package": fs.prop_gas,
        "has_pressure_change": True,
        "has_holdup": False,
        "delta_T_method": DeltaTMethod.counterCurrent,
        "tube_arrangement": TubeArrangement.inLine,
        "side_1_water_phase": "Vap",
        "has_radiation": False,
        "underwood":True,})

    # HP_SH3 = superheater 3
    fs.HP_SH3 = BoilerHeatExchanger(default={
        "side_1_property_package": fs.prop_water,
        "side_2_property_package": fs.prop_gas,
        "has_pressure_change": True,
        "has_holdup": False,
        "delta_T_method": DeltaTMethod.counterCurrent,
        "tube_arrangement": TubeArrangement.inLine,
        "side_1_water_phase": "Vap",
        "has_radiation": False,
        "underwood":True,})

    # HP_SH4 = superheater 4
    fs.HP_SH4 = BoilerHeatExchanger(default={
        "side_1_property_package": fs.prop_water,
        "side_2_property_package": fs.prop_gas,
        "has_pressure_change": True,
        "has_holdup": False,
        "delta_T_method": DeltaTMethod.counterCurrent,
        "tube_arrangement": TubeArrangement.inLine,
        "side_1_water_phase": "Vap",
        "has_radiation": False,
        "underwood":True,})

    # Low Pressure Evaporator Performance must be fixed to ~18 % of water inlet
    m.fs.LP_EVAP.vapor_frac_control = pyo.Var(initialize=0.18,
                                              doc="parameter to determine"
                                              "vapor flowrate")
    m.fs.LP_EVAP.vapor_frac_control.fix(0.18)
    @m.fs.LP_EVAP.Constraint(m.fs.time)
    def vap_fraceq(b, t):
        return b.tube.properties_out[0].vapor_frac == \
            m.fs.LP_EVAP.vapor_frac_control
    m.fs.LP_EVAP.heat_transfer_equation.deactivate()
    # unfixing U to calculate UA lumped parameter
    #m.fs.LP_EVAP.overall_heat_transfer_coefficient.unfix()
    #m.fs.LP_EVAP.overall_heat_transfer_coefficient.setub(250)

    # IP Inlet must be vaporized
    @m.fs.IP_EVAP.Constraint(m.fs.config.time)
    def sat_vap_eqn(b, t):
        return b.tube.properties_out[t].enth_mol == \
            b.tube.properties_out[t].enth_mol_sat_phase["Vap"] + 30
    # deactivate heat transfer to close DOF
    m.fs.IP_EVAP.heat_transfer_equation.deactivate()

    # HP Inlet must be vaporized
    @m.fs.HP_EVAP.Constraint(m.fs.config.time)
    def sat_vap_eqn2(b, t):
        return b.tube.properties_out[t].enth_mol == \
            b.tube.properties_out[t].enth_mol_sat_phase["Vap"] + 30
    return m
    # deactivate heat transfer to close DOF
    m.fs.HP_EVAP.heat_transfer_equation.deactivate()

    # Pressure Drop
    # Since tube geometry is an estimate of a real plant, off design conditions
    # can be calculated for overall heat transfer. However, pressure drop
    # cannot be estimated for this design. Therefore, pressure drop is fixed
    # across the system to match NETL Baseline Report Case B31B rev 4.
    # LP pressure drop has been ignored
    # Intermediate Pressure System - pressure drop
    m.fs.IP_ECON1.deltaP_tube_eqn.deactivate()
    m.fs.IP_ECON1.deltaP_tube.fix(-1.79E+05)
    m.fs.IP_ECON2.deltaP_tube_eqn.deactivate()
    m.fs.IP_ECON2.deltaP_tube.fix(-1.63E+05)
    m.fs.IP_EVAP.tube.deltaP.fix(-1.62E+05)
    m.fs.IP_SH1.deltaP_tube_eqn.deactivate()
    m.fs.IP_SH1.deltaP_tube.fix(-1.55E+05)
    m.fs.IP_SH2.deltaP_tube_eqn.deactivate()
    m.fs.IP_SH2.deltaP_tube.fix(-7.45E+04)
    m.fs.IP_SH3.deltaP_tube_eqn.deactivate()
    m.fs.IP_SH3.deltaP_tube.fix(-7.31E+04)
    # High Pressure System - pressure drop
    m.fs.HP_ECON1.deltaP_tube_eqn.deactivate()
    m.fs.HP_ECON1.deltaP_tube.fix(-9.79E+05)
    m.fs.HP_ECON2.deltaP_tube_eqn.deactivate()
    m.fs.HP_ECON2.deltaP_tube.fix(-9.38E+05)
    m.fs.HP_ECON3.deltaP_tube_eqn.deactivate()
    m.fs.HP_ECON3.deltaP_tube.fix(-8.96E+05)
    m.fs.HP_ECON4.deltaP_tube_eqn.deactivate()
    m.fs.HP_ECON4.deltaP_tube.fix(-8.62E+05)
    m.fs.HP_ECON5.deltaP_tube_eqn.deactivate()
    m.fs.HP_ECON5.deltaP_tube.fix(-8.27E+05)
    m.fs.HP_SH1.deltaP_tube_eqn.deactivate()
    m.fs.HP_SH1.deltaP_tube.fix(-8.27E+05)
    m.fs.HP_SH2.deltaP_tube_eqn.deactivate()
    m.fs.HP_SH2.deltaP_tube.fix(-7.93E+05)
    m.fs.HP_SH3.deltaP_tube_eqn.deactivate()
    m.fs.HP_SH3.deltaP_tube.fix(-5.52E+05)
    m.fs.HP_SH4.deltaP_tube_eqn.deactivate()
    m.fs.HP_SH4.deltaP_tube.fix(-5.45E+05)


def set_inputs(m):
    """
    Function to set the basic inputs for the unit models
    including geometry and fixed design and operating variables
    """
    # ============== LP Section ==============================================
    # ------------- LP ECONOMIZER --------------------------------------------
    # economizer design variables and parameters
    ITM = 0.0254  # inch to meter conversion
    m.fs.LP_ECON.tube_di.fix(0.038)  # (2-2*0.188)*ITM
    m.fs.LP_ECON.tube_thickness.fix(0.003)  # 0.188*ITM
    m.fs.LP_ECON.pitch_x.fix(0.09)  # 3.5*ITM
    m.fs.LP_ECON.pitch_y.fix(0.09)  # 5.03*ITM
    m.fs.LP_ECON.tube_length.fix(94.08/2)  # 53.41*105*ITM
    m.fs.LP_ECON.tube_nrow.fix(4)  # 31
    m.fs.LP_ECON.tube_ncol.fix(78)
    m.fs.LP_ECON.nrow_inlet.fix(2)
    m.fs.LP_ECON.delta_elevation.fix(1)
    # parameters
    # heat transfer resistance due to tube side fouling (water scales)
    m.fs.LP_ECON.tube_r_fouling = 0.000176
    # heat transfer resistance due to tube shell fouling (ash deposition)
    m.fs.LP_ECON.shell_r_fouling = 0.00088
    if m.fs.LP_ECON.config.has_radiation is True:
        m.fs.LP_ECON.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.LP_ECON.fcorrection_htc.fix(1.5)
    # correction factor for pressure drop calc tube side
    m.fs.LP_ECON.fcorrection_dp_tube.fix(1.0)
    # correction factor for pressure drop calc shell side
    m.fs.LP_ECON.fcorrection_dp_shell.fix(1.0)
    # prop_gas = m.prop_gas
    # LP Mixer1 no design variables -------------------------------------------
    # LP EVAP no design variables ---------------------------------------------
    # LP_DRUM LIQUID Outlet Splitter ------------------------------------------
    m.fs.Splitter1.split_fraction[0, "toIP"].fix(0.20626)
    # m.fs.Splitter1.split_fraction[0, "toHP"] = 0.704

    # IP Pump design ----------------------------------------------------------
    m.fs.IP_pump.efficiency_isentropic.fix(0.80)
    # HP Pump design ----------------------------------------------------------
    m.fs.HP_pump.efficiency_isentropic.fix(0.80)

    # LP EVAPORATOR -----------------------------------------------------------
    m.fs.LP_EVAP.area.fix(42032.0*0.5)

    # LP SH (Superheater) -----------------------------------------------------
    # Evaporator design variables and parameters
    ITM = 0.0254  # inch to meter conversion
    m.fs.LP_SH.tube_di.fix(0.051)
    m.fs.LP_SH.tube_thickness.fix(0.003)
    m.fs.LP_SH.pitch_x.fix(0.11)
    m.fs.LP_SH.pitch_y.fix(0.11)
    m.fs.LP_SH.tube_length.fix(7*9)
    m.fs.LP_SH.tube_nrow.fix(30)
    m.fs.LP_SH.tube_ncol.fix(10*3)
    m.fs.LP_SH.nrow_inlet.fix(10)
    m.fs.LP_SH.delta_elevation.fix(1)
    # parameters
    # heat transfer resistance due to tube side fouling (water scales)
    m.fs.LP_SH.tube_r_fouling = 0.000176
    # heat transfer resistance due to tube shell fouling (ash deposition)
    m.fs.LP_SH.shell_r_fouling = 0.00088
    if m.fs.LP_SH.config.has_radiation is True:
        m.fs.LP_SH.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.LP_SH.fcorrection_htc.fix(1.0)
    # correction factor for pressure drop calc tube side
    m.fs.LP_SH.fcorrection_dp_tube.fix(1.0)
    # correction factor for pressure drop calc shell side
    m.fs.LP_SH.fcorrection_dp_shell.fix(1.0)

    # -------------------------------------------------------------------------
    # IP Section --------------------------------------------------------------

    # IP_ECON1 (from IP Pump)--------------------------------------------------
    # Evaporator design variables and parameters
    ITM = 0.0254  # inch to meter conversion
    m.fs.IP_ECON1.tube_di.fix(0.038)
    m.fs.IP_ECON1.tube_thickness.fix(0.003)
    m.fs.IP_ECON1.pitch_x.fix(0.09)
    m.fs.IP_ECON1.pitch_y.fix(0.09)
    m.fs.IP_ECON1.tube_length.fix(7*14)
    m.fs.IP_ECON1.tube_nrow.fix(6)
    m.fs.IP_ECON1.tube_ncol.fix(12)
    m.fs.IP_ECON1.nrow_inlet.fix(3)
    m.fs.IP_ECON1.delta_elevation.fix(1.0)
    # parameters
    # heat transfer resistance due to tube side fouling (water scales)
    m.fs.IP_ECON1.tube_r_fouling = 0.000176
    # heat transfer resistance due to tube shell fouling (ash deposition)
    m.fs.IP_ECON1.shell_r_fouling = 0.00088
    if m.fs.IP_ECON1.config.has_radiation is True:
        m.fs.IP_ECON1.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.IP_ECON1.fcorrection_htc.fix(1.9)
    # correction factor for pressure drop calc tube side
    m.fs.IP_ECON1.fcorrection_dp_tube.fix(1.0)
    # correction factor for pressure drop calc shell side
    m.fs.IP_ECON1.fcorrection_dp_shell.fix(1.0)

    # IP_Splitter1 (from IP ECON) ---------------------------------------------
    m.fs.IP_Splitter1.split_fraction[0, "toNGPH"].fix(0.4592)

    # IP_ECON2 (from IP Splitter 1)--------------------------------------------
    # Evaporator design variables and parameters
    ITM = 0.0254  # inch to meter conversion
    m.fs.IP_ECON2.tube_di.fix(0.051)
    m.fs.IP_ECON2.tube_thickness.fix(0.003)
    m.fs.IP_ECON2.pitch_x.fix(0.09)
    m.fs.IP_ECON2.pitch_y.fix(0.09)
    m.fs.IP_ECON2.tube_length.fix(7)
    m.fs.IP_ECON2.tube_nrow.fix(15)
    m.fs.IP_ECON2.tube_ncol.fix(71)
    m.fs.IP_ECON2.nrow_inlet.fix(2)
    m.fs.IP_ECON2.delta_elevation.fix(12)
    # parameters
    # heat transfer resistance due to tube side fouling (water scales)
    m.fs.IP_ECON2.tube_r_fouling = 0.000176
    # heat transfer resistance due to tube shell fouling (ash deposition)
    m.fs.IP_ECON2.shell_r_fouling = 0.00088
    if m.fs.IP_ECON2.config.has_radiation is True:
        m.fs.IP_ECON2.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.IP_ECON2.fcorrection_htc.fix(1.0)
    # correction factor for pressure drop calc tube side
    m.fs.IP_ECON2.fcorrection_dp_tube.fix(1.0)
    # correction factor for pressure drop calc shell side
    m.fs.IP_ECON2.fcorrection_dp_shell.fix(1.0)

    # IP EVAPORATOR ----------------------------------------------------------
    m.fs.IP_EVAP.area.fix(8368.6)
    m.fs.IP_EVAP.overall_heat_transfer_coefficient.fix(150)

    # IP SH1 (Superheater) ----------------------------------------------------
    # Evaporator design variables and parameters
    ITM = 0.0254  # inch to meter conversion
    m.fs.IP_SH1.tube_di.fix(0.051)
    m.fs.IP_SH1.tube_thickness.fix(0.003)
    m.fs.IP_SH1.pitch_x.fix(0.11)
    m.fs.IP_SH1.pitch_y.fix(0.11)
    m.fs.IP_SH1.tube_length.fix(7*4)
    m.fs.IP_SH1.tube_nrow.fix(10)
    m.fs.IP_SH1.tube_ncol.fix(8)
    m.fs.IP_SH1.nrow_inlet.fix(5)
    m.fs.IP_SH1.delta_elevation.fix(1)
    # parameters
    # heat transfer resistance due to tube side fouling (water scales)
    m.fs.IP_SH1.tube_r_fouling = 0.000176
    # heat transfer resistance due to tube shell fouling (ash deposition)
    m.fs.IP_SH1.shell_r_fouling = 0.00088
    if m.fs.IP_SH1.config.has_radiation is True:
        m.fs.IP_SH1.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.IP_SH1.fcorrection_htc.fix(0.85)
    # correction factor for pressure drop calc tube side
    m.fs.IP_SH1.fcorrection_dp_tube.fix(1.0)
    # correction factor for pressure drop calc shell side
    m.fs.IP_SH1.fcorrection_dp_shell.fix(1.0)

    # IP SH2 (Superheater and Reheater 1) -------------------------------------
    # Evaporator design variables and parameters
    ITM = 0.0254  # inch to meter conversion
    m.fs.IP_SH2.tube_di.fix(0.0635)
    m.fs.IP_SH2.tube_thickness.fix(0.003)
    m.fs.IP_SH2.pitch_x.fix(0.11)
    m.fs.IP_SH2.pitch_y.fix(0.11)
    m.fs.IP_SH2.tube_length.fix(7*4)
    m.fs.IP_SH2.tube_nrow.fix(16)
    m.fs.IP_SH2.tube_ncol.fix(37)
    m.fs.IP_SH2.nrow_inlet.fix(4)
    m.fs.IP_SH2.delta_elevation.fix(1)
    # parameters
    # heat transfer resistance due to tube side fouling (water scales)
    m.fs.IP_SH2.tube_r_fouling = 0.000176
    # heat transfer resistance due to tube shell fouling (ash deposition)
    m.fs.IP_SH2.shell_r_fouling = 0.00088
    if m.fs.IP_SH2.config.has_radiation is True:
        m.fs.IP_SH2.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.IP_SH2.fcorrection_htc.fix(0.85)
    # correction factor for pressure drop calc tube side
    m.fs.IP_SH2.fcorrection_dp_tube.fix(1.0)
    # correction factor for pressure drop calc shell side
    m.fs.IP_SH2.fcorrection_dp_shell.fix(1.0)

    # IP SH3 (Superheater 3 and Reheater 2) -----------------------------------
    # Evaporator design variables and parameters
    ITM = 0.0254  # inch to meter conversion
    m.fs.IP_SH3.tube_di.fix(0.0635)
    m.fs.IP_SH3.tube_thickness.fix(0.003)
    m.fs.IP_SH3.pitch_x.fix(0.11)
    m.fs.IP_SH3.pitch_y.fix(0.11)
    m.fs.IP_SH3.tube_length.fix(7*10)
    m.fs.IP_SH3.tube_nrow.fix(16)
    m.fs.IP_SH3.tube_ncol.fix(37)
    m.fs.IP_SH3.nrow_inlet.fix(4)
    m.fs.IP_SH3.delta_elevation.fix(1)
    # parameters
    # heat transfer resistance due to tube side fouling (water scales)
    m.fs.IP_SH3.tube_r_fouling = 0.000176
    # heat transfer resistance due to tube shell fouling (ash deposition)
    m.fs.IP_SH3.shell_r_fouling = 0.00088
    if m.fs.IP_SH3.config.has_radiation is True:
        m.fs.IP_SH3.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.IP_SH3.fcorrection_htc.fix(0.95)
    # correction factor for pressure drop calc tube side
    m.fs.IP_SH3.fcorrection_dp_tube.fix(1.0)
    # correction factor for pressure drop calc shell side
    m.fs.IP_SH3.fcorrection_dp_shell.fix(1.0)

    # HP_ECON1 (from HP Pump)--------------------------------------------------
    # Evaporator design variables and parameters
    ITM = 0.0254  # inch to meter conversion
    m.fs.HP_ECON1.tube_di.fix(0.051)
    m.fs.HP_ECON1.tube_thickness.fix(0.003)
    m.fs.HP_ECON1.pitch_x.fix(0.1)
    m.fs.HP_ECON1.pitch_y.fix(0.1)
    m.fs.HP_ECON1.tube_length.fix(7*10)
    m.fs.HP_ECON1.tube_nrow.fix(8)
    m.fs.HP_ECON1.tube_ncol.fix(61)
    m.fs.HP_ECON1.nrow_inlet.fix(4)
    m.fs.HP_ECON1.delta_elevation.fix(1.0)
    # parameters
    # heat transfer resistance due to tube side fouling (water scales)
    m.fs.HP_ECON1.tube_r_fouling = 0.000176
    # heat transfer resistance due to tube shell fouling (ash deposition)
    m.fs.HP_ECON1.shell_r_fouling = 0.00088
    if m.fs.HP_ECON1.config.has_radiation is True:
        m.fs.HP_ECON1.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.HP_ECON1.fcorrection_htc.fix(1.0)
    # correction factor for pressure drop calc tube side
    m.fs.HP_ECON1.fcorrection_dp_tube.fix(0.05)
    # correction factor for pressure drop calc shell side
    m.fs.HP_ECON1.fcorrection_dp_shell.fix(0.05)

    # HP_ECON2 ----------------------------------------------------------------
    # Evaporator design variables and parameters
    ITM = 0.0254  # inch to meter conversion
    m.fs.HP_ECON2.tube_di.fix(0.051)
    m.fs.HP_ECON2.tube_thickness.fix(0.003)
    m.fs.HP_ECON2.pitch_x.fix(0.1)
    m.fs.HP_ECON2.pitch_y.fix(0.1)
    m.fs.HP_ECON2.tube_length.fix(7*10)
    m.fs.HP_ECON2.tube_nrow.fix(8)
    m.fs.HP_ECON2.tube_ncol.fix(61)
    m.fs.HP_ECON2.nrow_inlet.fix(4)
    m.fs.HP_ECON2.delta_elevation.fix(1.0)
    # parameters
    # heat transfer resistance due to tube side fouling (water scales)
    m.fs.HP_ECON2.tube_r_fouling = 0.000176
    # heat transfer resistance due to tube shell fouling (ash deposition)
    m.fs.HP_ECON2.shell_r_fouling = 0.00088
    if m.fs.HP_ECON2.config.has_radiation is True:
        m.fs.HP_ECON2.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.HP_ECON2.fcorrection_htc.fix(1.0)
    # correction factor for pressure drop calc tube side
    m.fs.HP_ECON2.fcorrection_dp_tube.fix(0.05)
    # correction factor for pressure drop calc shell side
    m.fs.HP_ECON2.fcorrection_dp_shell.fix(0.05)

    # HP_ECON3 ----------------------------------------------------------------
    # Evaporator design variables and parameters
    ITM = 0.0254  # inch to meter conversion
    m.fs.HP_ECON3.tube_di.fix(0.051)
    m.fs.HP_ECON3.tube_thickness.fix(0.003)
    m.fs.HP_ECON3.pitch_x.fix(0.1)
    m.fs.HP_ECON3.pitch_y.fix(0.1)
    m.fs.HP_ECON3.tube_length.fix(7*12)
    m.fs.HP_ECON3.tube_nrow.fix(8)
    m.fs.HP_ECON3.tube_ncol.fix(61)
    m.fs.HP_ECON3.nrow_inlet.fix(4)
    m.fs.HP_ECON3.delta_elevation.fix(1.0)
    # parameters
    # heat transfer resistance due to tube side fouling (water scales)
    m.fs.HP_ECON3.tube_r_fouling = 0.000176
    # heat transfer resistance due to tube shell fouling (ash deposition)
    m.fs.HP_ECON3.shell_r_fouling = 0.00088
    if m.fs.HP_ECON3.config.has_radiation is True:
        m.fs.HP_ECON3.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.HP_ECON3.fcorrection_htc.fix(1.2)
    # correction factor for pressure drop calc tube side
    m.fs.HP_ECON3.fcorrection_dp_tube.fix(0.05)
    # correction factor for pressure drop calc shell side
    m.fs.HP_ECON3.fcorrection_dp_shell.fix(0.05)

    # HP_ECON4 ----------------------------------------------------------------
    # Evaporator design variables and parameters
    ITM = 0.0254  # inch to meter conversion
    m.fs.HP_ECON4.tube_di.fix(0.051)
    m.fs.HP_ECON4.tube_thickness.fix(0.003)
    m.fs.HP_ECON4.pitch_x.fix(0.1)
    m.fs.HP_ECON4.pitch_y.fix(0.1)
    m.fs.HP_ECON4.tube_length.fix(7*10)
    m.fs.HP_ECON4.tube_nrow.fix(8)
    m.fs.HP_ECON4.tube_ncol.fix(61)
    m.fs.HP_ECON4.nrow_inlet.fix(4)
    m.fs.HP_ECON4.delta_elevation.fix(1.0)
    # parameters
    # heat transfer resistance due to tube side fouling (water scales)
    m.fs.HP_ECON4.tube_r_fouling = 0.000176
    # heat transfer resistance due to tube shell fouling (ash deposition)
    m.fs.HP_ECON4.shell_r_fouling = 0.00088
    if m.fs.HP_ECON4.config.has_radiation is True:
        m.fs.HP_ECON4.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.HP_ECON4.fcorrection_htc.fix(1.0)
    # correction factor for pressure drop calc tube side
    m.fs.HP_ECON4.fcorrection_dp_tube.fix(0.05)
    # correction factor for pressure drop calc shell side
    m.fs.HP_ECON4.fcorrection_dp_shell.fix(0.05)

    # HP_ECON5 ----------------------------------------------------------------
    # Evaporator design variables and parameters
    ITM = 0.0254  # inch to meter conversion
    m.fs.HP_ECON5.tube_di.fix(0.051)
    m.fs.HP_ECON5.tube_thickness.fix(0.003)
    m.fs.HP_ECON5.pitch_x.fix(0.1)
    m.fs.HP_ECON5.pitch_y.fix(0.1)
    m.fs.HP_ECON5.tube_length.fix(7*10)
    m.fs.HP_ECON5.tube_nrow.fix(8)
    m.fs.HP_ECON5.tube_ncol.fix(61)
    m.fs.HP_ECON5.nrow_inlet.fix(4)
    m.fs.HP_ECON5.delta_elevation.fix(1.0)
    # parameters
    # heat transfer resistance due to tube side fouling (water scales)
    m.fs.HP_ECON5.tube_r_fouling = 0.000176
    # heat transfer resistance due to tube shell fouling (ash deposition)
    m.fs.HP_ECON5.shell_r_fouling = 0.00088
    if m.fs.HP_ECON5.config.has_radiation is True:
        m.fs.HP_ECON5.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.HP_ECON5.fcorrection_htc.fix(1.0)
    # correction factor for pressure drop calc tube side
    m.fs.HP_ECON5.fcorrection_dp_tube.fix(0.05)
    # correction factor for pressure drop calc shell side
    m.fs.HP_ECON5.fcorrection_dp_shell.fix(0.05)

    # HP SH1 (Superheater 1) -----------------------------------
    # Evaporator design variables and parameters
    ITM = 0.0254  # inch to meter conversion
    m.fs.HP_SH1.tube_di.fix(0.0635)
    m.fs.HP_SH1.tube_thickness.fix(0.003)
    m.fs.HP_SH1.pitch_x.fix(0.11)
    m.fs.HP_SH1.pitch_y.fix(0.11)
    m.fs.HP_SH1.tube_length.fix(7*10)
    m.fs.HP_SH1.tube_nrow.fix(16)
    m.fs.HP_SH1.tube_ncol.fix(37)
    m.fs.HP_SH1.nrow_inlet.fix(4)
    m.fs.HP_SH1.delta_elevation.fix(1)
    # parameters
    # heat transfer resistance due to tube side fouling (water scales)
    m.fs.HP_SH1.tube_r_fouling = 0.000176
    # heat transfer resistance due to tube shell fouling (ash deposition)
    m.fs.HP_SH1.shell_r_fouling = 0.00088
    if m.fs.HP_SH1.config.has_radiation is True:
        m.fs.HP_SH1.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.HP_SH1.fcorrection_htc.fix(1.2)
    # correction factor for pressure drop calc tube side
    m.fs.HP_SH1.fcorrection_dp_tube.fix(1.0)
    # correction factor for pressure drop calc shell side
    m.fs.HP_SH1.fcorrection_dp_shell.fix(1.0)

    # HP SH2 (Superheater 2) -----------------------------------
    # Evaporator design variables and parameters
    ITM = 0.0254  # inch to meter conversion
    m.fs.HP_SH2.tube_di.fix(0.0635)
    m.fs.HP_SH2.tube_thickness.fix(0.003)
    m.fs.HP_SH2.pitch_x.fix(0.11)
    m.fs.HP_SH2.pitch_y.fix(0.11)
    m.fs.HP_SH2.tube_length.fix(7*10)
    m.fs.HP_SH2.tube_nrow.fix(16)
    m.fs.HP_SH2.tube_ncol.fix(37)
    m.fs.HP_SH2.nrow_inlet.fix(4)
    m.fs.HP_SH2.delta_elevation.fix(1)
    # parameters
    # heat transfer resistance due to tube side fouling (water scales)
    m.fs.HP_SH2.tube_r_fouling = 0.000176
    # heat transfer resistance due to tube shell fouling (ash deposition)
    m.fs.HP_SH2.shell_r_fouling = 0.00088
    if m.fs.HP_SH2.config.has_radiation is True:
        m.fs.HP_SH2.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.HP_SH2.fcorrection_htc.fix(0.95)
    # correction factor for pressure drop calc tube side
    m.fs.HP_SH2.fcorrection_dp_tube.fix(1.0)
    # correction factor for pressure drop calc shell side
    m.fs.HP_SH2.fcorrection_dp_shell.fix(1.0)

    # HP SH3 (Superheater 3) -----------------------------------
    # Evaporator design variables and parameters
    ITM = 0.0254  # inch to meter conversion
    m.fs.HP_SH3.tube_di.fix(0.0635)
    m.fs.HP_SH3.tube_thickness.fix(0.003)
    m.fs.HP_SH3.pitch_x.fix(0.11)
    m.fs.HP_SH3.pitch_y.fix(0.11)
    m.fs.HP_SH3.tube_length.fix(7*10)
    m.fs.HP_SH3.tube_nrow.fix(16)
    m.fs.HP_SH3.tube_ncol.fix(37)
    m.fs.HP_SH3.nrow_inlet.fix(4)
    m.fs.HP_SH3.delta_elevation.fix(1)
    # parameters
    # heat transfer resistance due to tube side fouling (water scales)
    m.fs.HP_SH3.tube_r_fouling = 0.000176
    # heat transfer resistance due to tube shell fouling (ash deposition)
    m.fs.HP_SH3.shell_r_fouling = 0.00088
    if m.fs.HP_SH3.config.has_radiation is True:
        m.fs.HP_SH3.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.HP_SH3.fcorrection_htc.fix(0.95)
    # correction factor for pressure drop calc tube side
    m.fs.HP_SH3.fcorrection_dp_tube.fix(1.0)
    # correction factor for pressure drop calc shell side
    m.fs.HP_SH3.fcorrection_dp_shell.fix(1.0)

    # HP SH4 (Superheater 4) -----------------------------------
    # Evaporator design variables and parameters
    m.fs.HP_SH4.tube_di.fix(0.0635)
    m.fs.HP_SH4.tube_thickness.fix(0.003)
    m.fs.HP_SH4.pitch_x.fix(0.11)
    m.fs.HP_SH4.pitch_y.fix(0.11)
    m.fs.HP_SH4.tube_length.fix(7*10)
    m.fs.HP_SH4.tube_nrow.fix(16)
    m.fs.HP_SH4.tube_ncol.fix(37)
    m.fs.HP_SH4.nrow_inlet.fix(4)
    m.fs.HP_SH4.delta_elevation.fix(1)
    # parameters
    # heat transfer resistance due to tube side fouling (water scales)
    m.fs.HP_SH4.tube_r_fouling = 0.000176
    # heat transfer resistance due to tube shell fouling (ash deposition)
    m.fs.HP_SH4.shell_r_fouling = 0.00088
    if m.fs.HP_SH4.config.has_radiation is True:
        m.fs.HP_SH4.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.HP_SH4.fcorrection_htc.fix(0.95)
    # correction factor for pressure drop calc tube side
    m.fs.HP_SH4.fcorrection_dp_tube.fix(1.0)
    # correction factor for pressure drop calc shell side
    m.fs.HP_SH4.fcorrection_dp_shell.fix(1.0)
    return m


def init_function(m):
    """Initialize unit models"""
    iscale.calculate_scaling_factors(m)
    fs = m.fs
    # prop_gas = m.prop_gas
    outlvl = idaeslog.INFO_LOW
    _log = idaeslog.getLogger(fs.name, outlvl, tag="unit")
    solve_log = idaeslog.getSolveLogger(fs.name, outlvl, tag="unit")
    solver = pyo.SolverFactory("ipopt")
    # solver.options = {
    #         "tol": 1e-8,
    #         "linear_solver": "ma27",
    #         "max_iter": 30,
    #         # "halt_on_ampl_error": "yes",
    # }
    options = {
            "tol": 1e-6,
            "linear_solver": "ma27",
            "max_iter": 30,
            # "halt_on_ampl_error": "yes",
    }

    # ============== LP Section ==============================================
    # ------------- LP ECONOMIZER --------------------------------------------
    # from FeedWater Heater3 = 1,399,485 lbm/hr; 634795 kg/hr;
    # Tube side (side 1)
    m.fs.LP_ECON.side_1_inlet.flow_mol[0].fix(9790.55)  # mol/s
    LP_ECON_h = iapws95.htpx(T=357.03*pyo.units.K, P=599844*pyo.units.Pa)
    m.fs.LP_ECON.side_1_inlet.enth_mol[0].fix(LP_ECON_h)
    m.fs.LP_ECON.side_1_inlet.pressure[0].fix(599844)  # Pa

    # FLUE GAS Inlet to LP_ECON (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPa abs)
    FGrate = 38446.11  # mol/s  from Baseline report table 5-22
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.LP_ECON.side_2_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    m.fs.LP_ECON.side_2_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    m.fs.LP_ECON.side_2_inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    m.fs.LP_ECON.side_2_inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    m.fs.LP_ECON.side_2_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    m.fs.LP_ECON.side_2_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    m.fs.LP_ECON.side_2_inlet.temperature[0].fix(433)  # 878.15)  # K
    m.fs.LP_ECON.side_2_inlet.pressure[0].fix(103421)  # Pa
    fs.LP_ECON.initialize(outlvl=outlvl)
    m.fs.LP_ECON.side_1_inlet.flow_mol[0].fix()
    m.fs.LP_ECON.side_1_inlet.enth_mol[0].fix()
    m.fs.LP_ECON.side_1_inlet.pressure[0].fix()
    solver.solve(m.fs.LP_ECON, tee=False)

    # ------------- LP Mixer 1 -----------------------------------------------
    _set_port(m.fs.Mixer1.LP_ECON, m.fs.LP_ECON.side_1_outlet)
    # Stream from Nat Gas Pre-heater:
    # F=58831.83 kg/hr / 18.01 kg/kmol / 3600 s/hr * 1000 mol/kmol,
    # T=333.15 K, 3.509e6Pa
    m.fs.Mixer1.Preheater.flow_mol[0].fix(907.89)  # mol/s
    Mixer1_h = iapws95.htpx(T=333.15*pyo.units.K, P=3.509e6*pyo.units.Pa)
    m.fs.Mixer1.Preheater.enth_mol[0].fix(Mixer1_h)  # j/mol
    m.fs.Mixer1.Preheater.pressure[0].fix(3.509e6)  # Pa
    m.fs.Mixer1.initialize()

    # ------------- LP EVAPORATOR --------------------------------------------
    # Tube side (side 1)
    # m.fs.LP_EVAP.tube_inlet.flow_mol[0].fix(10698)  # mol/s
    # LP_EVAP_h = iapws95.htpx(T=395.03*pyo.units.K, P=599844*pyo.units.Pa)
    # m.fs.LP_EVAP.tube_inlet.enth_mol[0].fix(48959.4365)
    # m.fs.LP_EVAP.tube_inlet.pressure[0].fix(599844.00)  # Pa
    _set_port(m.fs.LP_EVAP.tube_inlet, m.fs.Mixer1.outlet)

    # FLUE GAS Inlet to LP_EVAP (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPa abs)
    FGrate = 38446.11  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.LP_EVAP.shell_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    m.fs.LP_EVAP.shell_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    m.fs.LP_EVAP.shell_inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    m.fs.LP_EVAP.shell_inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    m.fs.LP_EVAP.shell_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    m.fs.LP_EVAP.shell_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    m.fs.LP_EVAP.shell_inlet.temperature[0].fix(543.31)  # K
    m.fs.LP_EVAP.shell_inlet.pressure[0].fix(103421)  # Pa
    m.fs.LP_EVAP.overall_heat_transfer_coefficient.fix(212)
    m.fs.LP_EVAP.vap_fraceq.deactivate()
    m.fs.LP_EVAP.tube_inlet.flow_mol[0].fix()
    m.fs.LP_EVAP.tube_inlet.enth_mol[0].fix()
    m.fs.LP_EVAP.tube_inlet.pressure[0].fix()
    m.fs.LP_EVAP.initialize()
    m.fs.LP_EVAP.vap_fraceq.activate()
    m.fs.LP_EVAP.overall_heat_transfer_coefficient.unfix()
    solver.solve(m.fs.LP_EVAP, tee=False)

    # -------------- LP Evaporator Initialization ----------------------------
    _set_port(m.fs.LP_DRUM.inlet, m.fs.LP_EVAP.tube_outlet)
    m.fs.LP_DRUM.initialize()

    m.fs.LP_DRUM.inlet.flow_mol[0].fix()
    m.fs.LP_DRUM.inlet.enth_mol[0].fix()
    m.fs.LP_DRUM.inlet.pressure[0].fix()
    solver.solve(m.fs.LP_DRUM, tee=False)
    # m.fs.LP_DRUM.display()
    m.fs.LP_DRUM.inlet.flow_mol[0].unfix()
    m.fs.LP_DRUM.inlet.enth_mol[0].unfix()
    m.fs.LP_DRUM.inlet.pressure[0].unfix()

    # -------------- LP_DRUM LIQUID Outlet Splitter --------------------------
    _set_port(m.fs.Splitter1.inlet, m.fs.LP_DRUM.liq_outlet)
    m.fs.Splitter1.initialize()

    # -------------- IP Pump init ------------------------------------------
    _set_port(m.fs.IP_pump.inlet, m.fs.Splitter1.toIP)
    m.fs.IP_pump.outlet.pressure[0].fix(4.385e6)
    m.fs.IP_pump.initialize()

    # -------------- HP Pump init ------------------------------------------
    _set_port(m.fs.HP_pump.inlet, m.fs.Splitter1.toHP)
    m.fs.HP_pump.outlet.pressure[0].fix(2.4359e7)
    m.fs.HP_pump.initialize()

    # flue gas bypass LP SH
    # FLUE GAS Inlet to splitter (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPa abs)
    FGrate = 38446.11  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.LP_FGsplit.inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    m.fs.LP_FGsplit.inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    m.fs.LP_FGsplit.inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    m.fs.LP_FGsplit.inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    m.fs.LP_FGsplit.inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    m.fs.LP_FGsplit.inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    m.fs.LP_FGsplit.inlet.temperature[0].fix(640.15)  # K
    m.fs.LP_FGsplit.inlet.pressure[0].fix(103421)  # Pa
    m.fs.LP_FGsplit.split_fraction[0, 'toLP_SH'].fix(0.5)
    m.fs.LP_FGsplit.initialize()

    # ------------- LP Superheater -------------------------------------------
    # Tube side (side 1)
    # m.fs.LP_SH.side_1_inlet.flow_mol[0].fix(2*4704.64)  # mol/s
    # LP_SH_h = iapws95.htpx(T=357.03*pyo.units.K, P=599844*pyo.units.Pa)
    # m.fs.LP_SH.side_1_inlet.enth_mol[0].fix(48959.4365)
    # m.fs.LP_SH.side_1_inlet.pressure[0].fix(256399.36)  # Pa
    _set_port(m.fs.LP_SH.side_1_inlet, m.fs.LP_DRUM.vap_outlet)
    _set_port(m.fs.LP_SH.side_2_inlet, m.fs.LP_FGsplit.toLP_SH)

    # # FLUE GAS Inlet to LP_EVAP (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPaabs)
    # FGrate = 38446.11*0.5  # mol/s
    # # Use FG molar composition to set component flow rates (baseline report)
    # m.fs.LP_SH.side_2_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    # m.fs.LP_SH.side_2_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    # m.fs.LP_SH.side_2_inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    # m.fs.LP_SH.side_2_inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    # m.fs.LP_SH.side_2_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    # m.fs.LP_SH.side_2_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    # m.fs.LP_SH.side_2_inlet.temperature[0].fix(640.15)  # 648.35 # K
    # m.fs.LP_SH.side_2_inlet.pressure[0].fix(103421)  # Pa

    fs.LP_SH.initialize(outlvl=0, optarg=options)

    m.fs.LP_SH.side_1_inlet.flow_mol[0].fix()
    m.fs.LP_SH.side_1_inlet.enth_mol[0].fix()
    m.fs.LP_SH.side_1_inlet.pressure[0].fix()

    _log.info("Low pressure system initialization - Completed")

    _set_port(m.fs.LP_Mixer2.fromLP_SH, m.fs.LP_SH.side_2_inlet)
    _set_port(m.fs.LP_Mixer2.bypass, m.fs.LP_FGsplit.toMixer)
    m.fs.LP_Mixer2.initialize()

    # IP Section --------------------------------------------------------------
    # IP ECONOMIZER 1 --------------------------------------------
    # Tube side (side 1)
    # m.fs.IP_ECON1.side_1_inlet.flow_mol[0].fix(2*4704.64)  # mol/s
    # IP_ECON1_h = iapws95.htpx(T=357.03*pyo.units.K, P=599844*pyo.units.Pa)
    # m.fs.IP_ECON1.side_1_inlet.enth_mol[0].fix(48959.4365)
    # m.fs.IP_ECON1.side_1_inlet.pressure[0].fix(256399.36)  # Pa
    _set_port(m.fs.IP_ECON1.side_1_inlet, m.fs.IP_pump.outlet)

    # FLUE GAS Inlet to LP_EVAP (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPa abs)
    FGrate = 38446.11  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.IP_ECON1.side_2_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    m.fs.IP_ECON1.side_2_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    m.fs.IP_ECON1.side_2_inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    m.fs.IP_ECON1.side_2_inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    m.fs.IP_ECON1.side_2_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    m.fs.IP_ECON1.side_2_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    m.fs.IP_ECON1.side_2_inlet.temperature[0].fix(579.46)  # K
    m.fs.IP_ECON1.side_2_inlet.pressure[0].fix(103421)  # Pa
    fs.IP_ECON1.initialize()
    m.fs.IP_ECON1.deltaP_tube_eqn.deactivate()
    m.fs.IP_ECON1.deltaP_tube.fix(-1.79E+05)
    m.fs.IP_ECON1.side_1_inlet.flow_mol[0].fix()
    m.fs.IP_ECON1.side_1_inlet.enth_mol[0].fix()
    m.fs.IP_ECON1.side_1_inlet.pressure[0].fix()
    solver.solve(m.fs.IP_ECON1, tee=False)
    m.fs.IP_ECON1.side_1_inlet.flow_mol[0].unfix()
    m.fs.IP_ECON1.side_1_inlet.enth_mol[0].unfix()
    m.fs.IP_ECON1.side_1_inlet.pressure[0].unfix()

    # IP Splitter 1: IP ECON to IP ECON2 and NGPH (nat gas pre heater) -------
    _set_port(m.fs.IP_Splitter1.inlet, m.fs.IP_ECON1.side_1_outlet)
    m.fs.IP_Splitter1.initialize()
    # m.fs.IP_Splitter1.toIP_ECON2.display()

    # IP ECONOMIZER 2 --------------------------------------------
    # Tube side (side 1)
    # m.fs.IP_ECON1.side_1_inlet.flow_mol[0].fix(4704.64)  # mol/s
    # IP_ECON1_h = iapws95.htpx(T=357.03*pyo.units.K, P=599844*pyo.units.Pa)
    # m.fs.IP_ECON1.side_1_inlet.enth_mol[0].fix(48959.4365)
    # m.fs.IP_ECON1.side_1_inlet.pressure[0].fix(256399.36)  # Pa
    _set_port(m.fs.IP_ECON2.side_1_inlet, m.fs.IP_Splitter1.toIP_ECON2)

    # FLUE GAS Inlet to LP_EVAP (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPa abs)
    FGrate = 38446.11  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.IP_ECON2.side_2_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    m.fs.IP_ECON2.side_2_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    m.fs.IP_ECON2.side_2_inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    m.fs.IP_ECON2.side_2_inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    m.fs.IP_ECON2.side_2_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    m.fs.IP_ECON2.side_2_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    m.fs.IP_ECON2.side_2_inlet.temperature[0].fix(607)  # K
    m.fs.IP_ECON2.side_2_inlet.pressure[0].fix(103421)  # Pa
    m.fs.IP_ECON2.initialize()
    m.fs.IP_ECON2.deltaP_tube_eqn.deactivate()
    m.fs.IP_ECON2.deltaP_tube.fix(-1.63E+05)
    m.fs.IP_ECON2.side_1_inlet.flow_mol[0].fix()
    m.fs.IP_ECON2.side_1_inlet.enth_mol[0].fix()
    m.fs.IP_ECON2.side_1_inlet.pressure[0].fix()
    solver.solve(m.fs.IP_ECON2, tee=False)

    _set_port(m.fs.IP_EVAP.tube_inlet, m.fs.IP_ECON2.side_1_outlet)
    m.fs.IP_EVAP.shell_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    m.fs.IP_EVAP.shell_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    m.fs.IP_EVAP.shell_inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    m.fs.IP_EVAP.shell_inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    m.fs.IP_EVAP.shell_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    m.fs.IP_EVAP.shell_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    m.fs.IP_EVAP.shell_inlet.temperature[0].fix(618.0)  # K
    m.fs.IP_EVAP.shell_inlet.pressure[0].fix(103421)  # Pa
    m.fs.IP_EVAP.sat_vap_eqn.deactivate()
    m.fs.IP_EVAP.heat_transfer_equation.activate()
    m.fs.IP_EVAP.tube.deltaP.fix(-1.62E+05)
    m.fs.IP_EVAP.initialize()
    m.fs.IP_EVAP.tube_inlet.flow_mol[0].fix()
    m.fs.IP_EVAP.tube_inlet.pressure[0].fix()
    m.fs.IP_EVAP.tube_inlet.enth_mol[0].fix()
    m.fs.IP_EVAP.overall_heat_transfer_coefficient.fix(150)
    # m.fs.IP_EVAP.tube.deltaP.fix(0)
    m.fs.IP_EVAP.sat_vap_eqn.activate()
    m.fs.IP_EVAP.heat_transfer_equation.deactivate()
    solver.solve(m.fs.IP_EVAP, tee=False)

    # IP Super Heater 1
    _set_port(m.fs.IP_SH1.side_1_inlet, m.fs.IP_EVAP.tube_outlet)
    # FLUE GAS Inlet to LP_EVAP (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPa abs)
    FGrate = 38446.11  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.IP_SH1.side_2_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    m.fs.IP_SH1.side_2_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    m.fs.IP_SH1.side_2_inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    m.fs.IP_SH1.side_2_inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    m.fs.IP_SH1.side_2_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    m.fs.IP_SH1.side_2_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    m.fs.IP_SH1.side_2_inlet.temperature[0].fix(668.15)  # 623)  # K
    m.fs.IP_SH1.side_2_inlet.pressure[0].fix(103421)  # Pa
    fs.IP_SH1.initialize()
    m.fs.IP_SH1.deltaP_tube_eqn.deactivate()
    m.fs.IP_SH1.deltaP_tube.fix(-1.55E+05)
    m.fs.IP_SH1.side_1_inlet.flow_mol[0].fix()
    m.fs.IP_SH1.side_1_inlet.enth_mol[0].fix()
    m.fs.IP_SH1.side_1_inlet.pressure[0].fix()
    solver.solve(m.fs.IP_SH1, tee=False)

    # Tube side HP Outlet Steam (Cold Reheat)
    # cold reheat 486241.9488 kg/hr *(1/18gr/mol)*(1000gr/kg)*(1hr/3600)= mol/s
    m.fs.IP_Splitter2.inlet.flow_mol[0].fix(7503.7337)  # mol/s
    IP_splitter2_h = iapws95.htpx(T=628.03*pyo.units.K, P=3.737e6*pyo.units.Pa)
    m.fs.IP_Splitter2.inlet.enth_mol[0].fix(IP_splitter2_h)
    m.fs.IP_Splitter2.inlet.pressure[0].fix(3.737e6)  # Pa
    m.fs.IP_Splitter2.split_fraction[0, "Cold_reheat"].fix(0.9941)
    m.fs.IP_Splitter2.split_fraction[0, "toEjector"].fix(0.00074)
    m.fs.IP_Splitter2.split_fraction[0, "toDryer"].fix(0.000274)
    # m.fs.IP_Splitter2.split_fraction[0, "toReclaimer"].fix(0.004801)
    ["Cold_reheat", "toEjector", "toReclaimer", "toDryer"]
    m.fs.IP_Splitter2.initialize()
    solver.solve(m.fs.IP_Splitter2, tee=False)

    # IP Mixer 1 (cold reheat + IP_SH1 outlet) ["IP_SH1", "Cold_reheat"]
    _set_port(m.fs.IP_Mixer1.IP_SH1, m.fs.IP_SH1.side_1_outlet)
    _set_port(m.fs.IP_Mixer1.Cold_reheat, m.fs.IP_Splitter2.Cold_reheat)
    m.fs.IP_Mixer1.IP_SH1.flow_mol[0].fix()
    m.fs.IP_Mixer1.IP_SH1.enth_mol[0].fix()
    m.fs.IP_Mixer1.IP_SH1.pressure[0].fix()
    m.fs.IP_Mixer1.Cold_reheat.flow_mol[0].fix()
    m.fs.IP_Mixer1.Cold_reheat.enth_mol[0].fix()
    m.fs.IP_Mixer1.Cold_reheat.pressure[0].fix()
    m.fs.IP_Mixer1.initialize()


    # IP Super Heater 2 and Reheater 1
    _set_port(m.fs.IP_SH2.side_1_inlet, m.fs.IP_Mixer1.outlet)
    # FLUE GAS Inlet to LP_EVAP (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPa abs)
    FGrate = 38446.11  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.IP_SH2.side_2_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    m.fs.IP_SH2.side_2_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    m.fs.IP_SH2.side_2_inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    m.fs.IP_SH2.side_2_inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    m.fs.IP_SH2.side_2_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    m.fs.IP_SH2.side_2_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    m.fs.IP_SH2.side_2_inlet.temperature[0].fix(802.35)  # 878.15)  # K
    m.fs.IP_SH2.side_2_inlet.pressure[0].fix(103421)  # Pa
    m.fs.IP_SH2.initialize()
    m.fs.IP_SH2.deltaP_tube_eqn.deactivate()
    m.fs.IP_SH2.deltaP_tube.fix(-7.45E+04)
    m.fs.IP_SH2.side_1_inlet.flow_mol[0].fix()
    m.fs.IP_SH2.side_1_inlet.enth_mol[0].fix()
    m.fs.IP_SH2.side_1_inlet.pressure[0].fix()
    solver.solve(m.fs.IP_SH2, tee=False)

    # IP Super Heater 3 and Reheater 2
    _set_port(m.fs.IP_SH3.side_1_inlet, m.fs.IP_SH2.side_1_outlet)
    # FLUE GAS Inlet to LP_EVAP (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPa abs)
    FGrate = 38446.11  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.IP_SH3.side_2_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    m.fs.IP_SH3.side_2_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    m.fs.IP_SH3.side_2_inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    m.fs.IP_SH3.side_2_inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    m.fs.IP_SH3.side_2_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    m.fs.IP_SH3.side_2_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    m.fs.IP_SH3.side_2_inlet.temperature[0].fix(868.35)  # 878.15)  # K
    m.fs.IP_SH3.side_2_inlet.pressure[0].fix(103421)  # Pa
    m.fs.IP_SH3.initialize()
    m.fs.IP_SH3.deltaP_tube_eqn.deactivate()
    m.fs.IP_SH3.deltaP_tube.fix(-7.31E+04)
    m.fs.IP_SH3.side_1_inlet.flow_mol[0].fix()
    m.fs.IP_SH3.side_1_inlet.enth_mol[0].fix()
    m.fs.IP_SH3.side_1_inlet.pressure[0].fix()
    solver.solve(m.fs.IP_SH3, tee=False)
    _log.info("Intermediate pressure system initialization - Completed")

    # -------------------------------------------------------------------------
    # High Pressure System ----------------------------------------------------
    # HP ECONOMIZER 1
    _set_port(m.fs.HP_ECON1.side_1_inlet, m.fs.HP_pump.outlet)
    # FLUE GAS Inlet to LP_EVAP (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPa abs)
    FGrate = 38446.11  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.HP_ECON1.side_2_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    m.fs.HP_ECON1.side_2_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    m.fs.HP_ECON1.side_2_inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    m.fs.HP_ECON1.side_2_inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    m.fs.HP_ECON1.side_2_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    m.fs.HP_ECON1.side_2_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    m.fs.HP_ECON1.side_2_inlet.temperature[0].fix(566.35)  # K
    m.fs.HP_ECON1.side_2_inlet.pressure[0].fix(103421)  # Pa
    m.fs.HP_ECON1.initialize()
    m.fs.HP_ECON1.deltaP_tube_eqn.deactivate()
    m.fs.HP_ECON1.deltaP_tube.fix(-9.79E+05)
    m.fs.HP_ECON1.side_1_inlet.flow_mol[0].fix()
    m.fs.HP_ECON1.side_1_inlet.enth_mol[0].fix()
    m.fs.HP_ECON1.side_1_inlet.pressure[0].fix()
    solver.solve(m.fs.HP_ECON1, tee=False)

    # HP ECONOMIZER 2
    _set_port(m.fs.HP_ECON2.side_1_inlet, m.fs.HP_ECON1.side_1_outlet)
    # FLUE GAS Inlet to LP_EVAP (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPa abs)
    FGrate = 38446.11  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.HP_ECON2.side_2_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    m.fs.HP_ECON2.side_2_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    m.fs.HP_ECON2.side_2_inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    m.fs.HP_ECON2.side_2_inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    m.fs.HP_ECON2.side_2_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    m.fs.HP_ECON2.side_2_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    m.fs.HP_ECON2.side_2_inlet.temperature[0].fix(596.35)   # K
    m.fs.HP_ECON2.side_2_inlet.pressure[0].fix(103421)  # Pa
    m.fs.HP_ECON2.initialize()
    m.fs.HP_ECON2.deltaP_tube_eqn.deactivate()
    m.fs.HP_ECON2.deltaP_tube.fix(-9.38E+05)
    m.fs.HP_ECON2.side_1_inlet.flow_mol[0].fix()
    m.fs.HP_ECON2.side_1_inlet.enth_mol[0].fix()
    m.fs.HP_ECON2.side_1_inlet.pressure[0].fix()
    solver.solve(m.fs.HP_ECON2, tee=False)

    # HP ECONOMIZER 3
    _set_port(m.fs.HP_ECON3.side_1_inlet, m.fs.HP_ECON2.side_1_outlet)
    # FLUE GAS Inlet to LP_EVAP (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPa abs)
    FGrate = 38446.11  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.HP_ECON3.side_2_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    m.fs.HP_ECON3.side_2_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    m.fs.HP_ECON3.side_2_inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    m.fs.HP_ECON3.side_2_inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    m.fs.HP_ECON3.side_2_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    m.fs.HP_ECON3.side_2_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    m.fs.HP_ECON3.side_2_inlet.temperature[0].fix(667.35)  # 878.15)  # K
    m.fs.HP_ECON3.side_2_inlet.pressure[0].fix(103421)  # Pa
    m.fs.HP_ECON3.initialize()
    m.fs.HP_ECON3.deltaP_tube_eqn.deactivate()
    m.fs.HP_ECON3.deltaP_tube.fix(-8.96E+05)
    m.fs.HP_ECON3.side_1_inlet.flow_mol[0].fix()
    m.fs.HP_ECON3.side_1_inlet.enth_mol[0].fix()
    m.fs.HP_ECON3.side_1_inlet.pressure[0].fix()
    solver.solve(m.fs.HP_ECON3, tee=False)

    # HP ECONOMIZER 4
    _set_port(m.fs.HP_ECON4.side_1_inlet, m.fs.HP_ECON3.side_1_outlet)
    # FLUE GAS Inlet to LP_EVAP (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPa abs)
    FGrate = 38446.11  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.HP_ECON4.side_2_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    m.fs.HP_ECON4.side_2_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    m.fs.HP_ECON4.side_2_inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    m.fs.HP_ECON4.side_2_inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    m.fs.HP_ECON4.side_2_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    m.fs.HP_ECON4.side_2_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    m.fs.HP_ECON4.side_2_inlet.temperature[0].fix(624.35)  # K
    m.fs.HP_ECON4.side_2_inlet.pressure[0].fix(103421)  # Pa
    m.fs.HP_ECON4.initialize()
    m.fs.HP_ECON4.deltaP_tube_eqn.deactivate()
    m.fs.HP_ECON4.deltaP_tube.fix(-8.62E+05)
    m.fs.HP_ECON4.side_1_inlet.flow_mol[0].fix()
    m.fs.HP_ECON4.side_1_inlet.enth_mol[0].fix()
    m.fs.HP_ECON4.side_1_inlet.pressure[0].fix()
    solver.solve(m.fs.HP_ECON4, tee=False)

    # HP ECONOMIZER 5
    _set_port(m.fs.HP_ECON5.side_1_inlet, m.fs.HP_ECON4.side_1_outlet)
    # FLUE GAS Inlet to LP_EVAP (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPa abs)
    FGrate = 38446.11  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.HP_ECON5.side_2_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    m.fs.HP_ECON5.side_2_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    m.fs.HP_ECON5.side_2_inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    m.fs.HP_ECON5.side_2_inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    m.fs.HP_ECON5.side_2_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    m.fs.HP_ECON5.side_2_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    m.fs.HP_ECON5.side_2_inlet.temperature[0].fix(631.15)  # K
    m.fs.HP_ECON5.side_2_inlet.pressure[0].fix(103421)  # Pa
    m.fs.HP_ECON5.initialize()
    m.fs.HP_ECON5.deltaP_tube_eqn.deactivate()
    m.fs.HP_ECON5.deltaP_tube.fix(-8.27E+05)
    m.fs.HP_ECON5.side_1_inlet.flow_mol[0].fix()
    m.fs.HP_ECON5.side_1_inlet.enth_mol[0].fix()
    m.fs.HP_ECON5.side_1_inlet.pressure[0].fix()
    solver.solve(m.fs.HP_ECON5, tee=False)

    # HP EVAPORATOR ----------------------------------------------------------
    m.fs.HP_EVAP.area.fix(8368.6)
    m.fs.HP_EVAP.overall_heat_transfer_coefficient.fix(150)

    _set_port(m.fs.HP_EVAP.tube_inlet, m.fs.HP_ECON5.side_1_outlet)

    # FLUE GAS Inlet to LP_EVAP (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPa abs)
    FGrate = 38446.11  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.HP_EVAP.shell_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    m.fs.HP_EVAP.shell_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    m.fs.HP_EVAP.shell_inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    m.fs.HP_EVAP.shell_inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    m.fs.HP_EVAP.shell_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    m.fs.HP_EVAP.shell_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    m.fs.HP_EVAP.shell_inlet.temperature[0].fix(729.15)  # K
    m.fs.HP_EVAP.shell_inlet.pressure[0].fix(103421)  # Pa
    m.fs.HP_EVAP.sat_vap_eqn2.deactivate()
    m.fs.HP_EVAP.heat_transfer_equation.activate()
    m.fs.HP_EVAP.initialize()
    m.fs.HP_EVAP.sat_vap_eqn2.activate()
    m.fs.HP_EVAP.heat_transfer_equation.deactivate()
    m.fs.HP_EVAP.tube_inlet.flow_mol[0].fix()
    m.fs.HP_EVAP.tube_inlet.enth_mol[0].fix()
    m.fs.HP_EVAP.tube_inlet.pressure[0].fix()
    solver.solve(m.fs.HP_EVAP, tee=False)

    # HP superheater 1
    _set_port(m.fs.HP_SH1.side_1_inlet, m.fs.HP_EVAP.tube_outlet)
    # FLUE GAS Inlet to LP_EVAP (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPa abs)
    FGrate = 38446.11  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.HP_SH1.side_2_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    m.fs.HP_SH1.side_2_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    m.fs.HP_SH1.side_2_inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    m.fs.HP_SH1.side_2_inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    m.fs.HP_SH1.side_2_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    m.fs.HP_SH1.side_2_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    m.fs.HP_SH1.side_2_inlet.temperature[0].fix(760.35)  # K
    m.fs.HP_SH1.side_2_inlet.pressure[0].fix(103421)  # Pa
    fs.HP_SH1.initialize()
    m.fs.HP_SH1.deltaP_tube_eqn.deactivate()
    m.fs.HP_SH1.deltaP_tube.fix(-8.27E+05)
    m.fs.HP_SH1.side_1_inlet.flow_mol[0].fix()
    m.fs.HP_SH1.side_1_inlet.enth_mol[0].fix()
    m.fs.HP_SH1.side_1_inlet.pressure[0].fix()
    solver.solve(m.fs.HP_SH1, tee=False)

    # HP superheater 2
    _set_port(m.fs.HP_SH2.side_1_inlet, m.fs.HP_SH1.side_1_outlet)
    # FLUE GAS Inlet to LP_EVAP (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPa abs)
    FGrate = 38446.11  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.HP_SH2.side_2_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    m.fs.HP_SH2.side_2_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    m.fs.HP_SH2.side_2_inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    m.fs.HP_SH2.side_2_inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    m.fs.HP_SH2.side_2_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    m.fs.HP_SH2.side_2_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    m.fs.HP_SH2.side_2_inlet.temperature[0].fix(806.35)  # 878.15)  # K
    m.fs.HP_SH2.side_2_inlet.pressure[0].fix(103421)  # Pa
    fs.HP_SH2.initialize()
    m.fs.HP_SH2.deltaP_tube_eqn.deactivate()
    m.fs.HP_SH2.deltaP_tube.fix(-7.93E+05)
    m.fs.HP_SH2.side_1_inlet.flow_mol[0].fix()
    m.fs.HP_SH2.side_1_inlet.enth_mol[0].fix()
    m.fs.HP_SH2.side_1_inlet.pressure[0].fix()
    solver.solve(m.fs.HP_SH2, tee=False)

    # HP superheater 3
    _set_port(m.fs.HP_SH3.side_1_inlet, m.fs.HP_SH2.side_1_outlet)
    # FLUE GAS Inlet to LP_EVAP (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPa abs)
    FGrate = 38446.11  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.HP_SH3.side_2_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    m.fs.HP_SH3.side_2_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    m.fs.HP_SH3.side_2_inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    m.fs.HP_SH3.side_2_inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    m.fs.HP_SH3.side_2_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    m.fs.HP_SH3.side_2_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    m.fs.HP_SH3.side_2_inlet.temperature[0].fix(835.35)  # 878.15)  # K
    m.fs.HP_SH3.side_2_inlet.pressure[0].fix(103421)  # Pa
    m.fs.HP_SH3.initialize()
    m.fs.HP_SH3.deltaP_tube_eqn.deactivate()
    m.fs.HP_SH3.deltaP_tube.fix(-5.52E+05)
    m.fs.HP_SH3.side_1_inlet.flow_mol[0].fix()
    m.fs.HP_SH3.side_1_inlet.enth_mol[0].fix()
    m.fs.HP_SH3.side_1_inlet.pressure[0].fix()
    solver.solve(m.fs.HP_SH3, tee=False)

    # HP superheater 4
    _set_port(m.fs.HP_SH4.side_1_inlet, m.fs.HP_SH3.side_1_outlet)
    # FLUE GAS Inlet to LP_EVAP (F = 138406 kgmol/hr, T = 15 C, P=0.1 MPa abs)
    FGrate = 38446.11  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.HP_SH4.side_2_inlet.flow_mol_comp[0, "H2O"].fix(FGrate*0.0875)
    m.fs.HP_SH4.side_2_inlet.flow_mol_comp[0, "CO2"].fix(FGrate*0.0408)
    m.fs.HP_SH4.side_2_inlet.flow_mol_comp[0, "N2"].fix(FGrate*0.75)
    m.fs.HP_SH4.side_2_inlet.flow_mol_comp[0, "O2"].fix(FGrate*0.12)
    m.fs.HP_SH4.side_2_inlet.flow_mol_comp[0, "NO"].fix(FGrate*0.001)
    m.fs.HP_SH4.side_2_inlet.flow_mol_comp[0, "SO2"].fix(FGrate*0.0007)
    m.fs.HP_SH4.side_2_inlet.temperature[0].fix(878.15)  # 878.15)  # K
    m.fs.HP_SH4.side_2_inlet.pressure[0].fix(103421)  # Pa
    fs.HP_SH4.initialize()
    m.fs.HP_SH1.deltaP_tube_eqn.deactivate()
    m.fs.HP_SH1.deltaP_tube.fix(-5.45E+05)
    m.fs.HP_SH4.side_1_inlet.flow_mol[0].fix()
    m.fs.HP_SH4.side_1_inlet.enth_mol[0].fix()
    m.fs.HP_SH4.side_1_inlet.pressure[0].fix()
    solver.solve(m.fs.HP_SH4, tee=False)
    _log.info("High pressure system initialization - Completed")

    return m


def set_arcs(m):
    # -------------------------------------------------------------------------
    # Arcs ===================================================================
    # LP_ECON to LP DRUM - steam - to LP EVAP, LP EVAP to LP SH
    # LP EVAP - liq outlet - to Pumps (HP and IP)
    # Mixer to LP_ECON
    m.fs.lp02 = Arc(
        source=m.fs.LP_ECON.side_1_outlet, destination=m.fs.Mixer1.LP_ECON
    )

    # LP_ECON to LP_Mixer then LP_Mixer to LP_EVAP
    m.fs.lp04 = Arc(
        source=m.fs.Mixer1.outlet, destination=m.fs.LP_EVAP.tube_inlet
    )
    m.fs.LP_EVAP.tube_inlet.flow_mol[0].unfix()
    m.fs.LP_EVAP.tube_inlet.enth_mol[0].unfix()
    m.fs.LP_EVAP.tube_inlet.pressure[0].unfix()

    # LP_DRUM STEAM OUTLET: inlet to evaporator
    # LP_DRUM to LP_EVAP
    m.fs.lp05 = Arc(
        source=m.fs.LP_EVAP.tube_outlet, destination=m.fs.LP_DRUM.inlet
    )

    # LP_EVAP to LP_SH (superheater)
    m.fs.lp10 = Arc(
        source=m.fs.LP_DRUM.vap_outlet, destination=m.fs.LP_SH.side_1_inlet
    )
    m.fs.LP_SH.side_1_inlet.flow_mol[0].unfix()
    m.fs.LP_SH.side_1_inlet.enth_mol[0].unfix()
    m.fs.LP_SH.side_1_inlet.pressure[0].unfix()
    # LP_DRUM LIQUID OUTLET: inlet to splitters
    # Arc to connect liquid outlet LP evap to splitter
    m.fs.lp06 = Arc(
        source=m.fs.LP_DRUM.liq_outlet, destination=m.fs.Splitter1.inlet
    )

    # arc to connect splitter to IP Pump
    m.fs.lp08 = Arc(
        source=m.fs.Splitter1.toIP, destination=m.fs.IP_pump.inlet
    )

    # arc to connect splitter to HP pump
    m.fs.lp09 = Arc(
        source=m.fs.Splitter1.toHP, destination=m.fs.HP_pump.inlet
    )

    # IP Section --------------------------------------------------------------
    # arc to connect IP Pump to IP_ECON1
    m.fs.ip01 = Arc(
        source=m.fs.IP_pump.outlet, destination=m.fs.IP_ECON1.side_1_inlet
    )
    # IP_ECON1 to IP_splitter 1
    m.fs.ip02 = Arc(
        source=m.fs.IP_ECON1.side_1_outlet, destination=m.fs.IP_Splitter1.inlet
    )

    # IP_splitter_1 to IP_ECON2
    m.fs.ip03 = Arc(
        source=m.fs.IP_Splitter1.toIP_ECON2,
        destination=m.fs.IP_ECON2.side_1_inlet)
    m.fs.IP_ECON2.side_1_inlet.flow_mol[0].unfix()
    m.fs.IP_ECON2.side_1_inlet.enth_mol[0].unfix()
    m.fs.IP_ECON2.side_1_inlet.pressure[0].unfix()

    # IP_ECON2 to IP_EVAP
    m.fs.ip05 = Arc(
        source=m.fs.IP_ECON2.side_1_outlet,
        destination=m.fs.IP_EVAP.tube_inlet)
    m.fs.IP_EVAP.tube_inlet.flow_mol[0].unfix()
    m.fs.IP_EVAP.tube_inlet.pressure[0].unfix()
    m.fs.IP_EVAP.tube_inlet.enth_mol[0].unfix()

    # IP_EVAP to IP_SH1
    m.fs.ip06 = Arc(
        source=m.fs.IP_EVAP.tube_outlet,
        destination=m.fs.IP_SH1.side_1_inlet)
    m.fs.IP_SH1.side_1_inlet.flow_mol[0].unfix()
    m.fs.IP_SH1.side_1_inlet.enth_mol[0].unfix()
    m.fs.IP_SH1.side_1_inlet.pressure[0].unfix()

    # IP_SH1 to IP_Mixer1
    m.fs.ip07 = Arc(
        source=m.fs.IP_SH1.side_1_outlet,
        destination=m.fs.IP_Mixer1.IP_SH1)

    # cold reheat to IP_Mixer1
    m.fs.ip15 = Arc(
        source=m.fs.IP_Splitter2.Cold_reheat,
        destination=m.fs.IP_Mixer1.Cold_reheat)
    m.fs.IP_Mixer1.IP_SH1.flow_mol[0].unfix()
    m.fs.IP_Mixer1.IP_SH1.enth_mol[0].unfix()
    m.fs.IP_Mixer1.IP_SH1.pressure[0].unfix()
    m.fs.IP_Mixer1.Cold_reheat.flow_mol[0].unfix()
    m.fs.IP_Mixer1.Cold_reheat.enth_mol[0].unfix()
    m.fs.IP_Mixer1.Cold_reheat.pressure[0].unfix()

    # IP_Mixer1 to IP_SH2
    m.fs.ip08 = Arc(
        source=m.fs.IP_Mixer1.outlet,
        destination=m.fs.IP_SH2.side_1_inlet)
    m.fs.IP_SH2.side_1_inlet.flow_mol[0].unfix()
    m.fs.IP_SH2.side_1_inlet.enth_mol[0].unfix()
    m.fs.IP_SH2.side_1_inlet.pressure[0].unfix()

    # IP_SH2 to IP_SH3
    m.fs.ip09 = Arc(
        source=m.fs.IP_SH2.side_1_outlet,
        destination=m.fs.IP_SH3.side_1_inlet)
    m.fs.IP_SH3.side_1_inlet.flow_mol[0].unfix()
    m.fs.IP_SH3.side_1_inlet.enth_mol[0].unfix()
    m.fs.IP_SH3.side_1_inlet.pressure[0].unfix()

    # HP_pump to HP_ECON1
    m.fs.hp01 = Arc(
        source=m.fs.HP_pump.outlet,
        destination=m.fs.HP_ECON1.side_1_inlet)
    m.fs.HP_ECON1.side_1_inlet.flow_mol[0].unfix()
    m.fs.HP_ECON1.side_1_inlet.enth_mol[0].unfix()
    m.fs.HP_ECON1.side_1_inlet.pressure[0].unfix()

    # HP_ECON1 to HP_ECON2
    m.fs.hp02 = Arc(
        source=m.fs.HP_ECON1.side_1_outlet,
        destination=m.fs.HP_ECON2.side_1_inlet)
    m.fs.HP_ECON2.side_1_inlet.flow_mol[0].unfix()
    m.fs.HP_ECON2.side_1_inlet.enth_mol[0].unfix()
    m.fs.HP_ECON2.side_1_inlet.pressure[0].unfix()

    # HP_ECON2 to HP_ECON3
    m.fs.hp03 = Arc(
        source=m.fs.HP_ECON2.side_1_outlet,
        destination=m.fs.HP_ECON3.side_1_inlet)
    m.fs.HP_ECON3.side_1_inlet.flow_mol[0].unfix()
    m.fs.HP_ECON3.side_1_inlet.enth_mol[0].unfix()
    m.fs.HP_ECON3.side_1_inlet.pressure[0].unfix()

    # HP_ECON3 to HP_ECON4
    m.fs.hp04 = Arc(
        source=m.fs.HP_ECON3.side_1_outlet,
        destination=m.fs.HP_ECON4.side_1_inlet)
    m.fs.HP_ECON4.side_1_inlet.flow_mol[0].unfix()
    m.fs.HP_ECON4.side_1_inlet.enth_mol[0].unfix()
    m.fs.HP_ECON4.side_1_inlet.pressure[0].unfix()

    # HP_ECON4 to HP_ECON5
    m.fs.hp05 = Arc(
        source=m.fs.HP_ECON4.side_1_outlet,
        destination=m.fs.HP_ECON5.side_1_inlet)
    m.fs.HP_ECON5.side_1_inlet.flow_mol[0].unfix()
    m.fs.HP_ECON5.side_1_inlet.enth_mol[0].unfix()
    m.fs.HP_ECON5.side_1_inlet.pressure[0].unfix()

    # HP_ECON5 to HP_EVAP
    m.fs.hp06 = Arc(
        source=m.fs.HP_ECON5.side_1_outlet, destination=m.fs.HP_EVAP.tube_inlet
    )
    m.fs.HP_EVAP.tube_inlet.flow_mol[0].unfix()
    m.fs.HP_EVAP.tube_inlet.enth_mol[0].unfix()
    m.fs.HP_EVAP.tube_inlet.pressure[0].unfix()

    # HP_EVAP to HP_SH1
    m.fs.hp07 = Arc(
        source=m.fs.HP_EVAP.tube_outlet, destination=m.fs.HP_SH1.side_1_inlet
    )
    m.fs.HP_SH1.side_1_inlet.flow_mol[0].unfix()
    m.fs.HP_SH1.side_1_inlet.enth_mol[0].unfix()
    m.fs.HP_SH1.side_1_inlet.pressure[0].unfix()

    # HP_SH1 to HP_SH2
    m.fs.hp08 = Arc(
        source=m.fs.HP_SH1.side_1_outlet, destination=m.fs.HP_SH2.side_1_inlet
    )
    m.fs.HP_SH2.side_1_inlet.flow_mol[0].unfix()
    m.fs.HP_SH2.side_1_inlet.enth_mol[0].unfix()
    m.fs.HP_SH2.side_1_inlet.pressure[0].unfix()

    # HP_SH2 to HP_SH3
    m.fs.hp09 = Arc(
        source=m.fs.HP_SH2.side_1_outlet, destination=m.fs.HP_SH3.side_1_inlet
    )
    m.fs.HP_SH3.side_1_inlet.flow_mol[0].unfix()
    m.fs.HP_SH3.side_1_inlet.enth_mol[0].unfix()
    m.fs.HP_SH3.side_1_inlet.pressure[0].unfix()

    # HP_SH3 to HP_SH4
    m.fs.hp10 = Arc(
        source=m.fs.HP_SH3.side_1_outlet, destination=m.fs.HP_SH4.side_1_inlet
    )
    m.fs.HP_SH4.side_1_inlet.flow_mol[0].unfix()
    m.fs.HP_SH4.side_1_inlet.enth_mol[0].unfix()
    m.fs.HP_SH4.side_1_inlet.pressure[0].unfix()

    # Flue Gas Route ----------------------------------------------------------
    # HP_SH4 to IP_SH3 or IP_Reheater 2
    m.fs.g09 = Arc(
        source=m.fs.HP_SH4.side_2_outlet, destination=m.fs.IP_SH3.side_2_inlet)
    m.fs.IP_SH3.side_2_inlet.flow_mol_comp[:, :].unfix()
    m.fs.IP_SH3.side_2_inlet.pressure[0].unfix()
    m.fs.IP_SH3.side_2_inlet.temperature[0].unfix()

    # IP_SH3 to HP_SH3
    m.fs.g10 = Arc(
        source=m.fs.IP_SH3.side_2_outlet,
        destination=m.fs.HP_SH3.side_2_inlet)
    m.fs.HP_SH3.side_2_inlet.flow_mol_comp[:, :].unfix()
    m.fs.HP_SH3.side_2_inlet.pressure[0].unfix()
    m.fs.HP_SH3.side_2_inlet.temperature[0].unfix()

    # HP_SH3 to HP_SH2
    m.fs.g11 = Arc(
        source=m.fs.HP_SH3.side_2_outlet,
        destination=m.fs.HP_SH2.side_2_inlet)
    m.fs.HP_SH2.side_2_inlet.flow_mol_comp[:, :].unfix()
    m.fs.HP_SH2.side_2_inlet.pressure[0].unfix()
    m.fs.HP_SH2.side_2_inlet.temperature[0].unfix()

    # HP_SH2 to IP_SH2
    m.fs.g12 = Arc(
        source=m.fs.HP_SH2.side_2_outlet,
        destination=m.fs.IP_SH2.side_2_inlet)
    m.fs.IP_SH2.side_2_inlet.flow_mol_comp[:, :].unfix()
    m.fs.IP_SH2.side_2_inlet.pressure[0].unfix()
    m.fs.IP_SH2.side_2_inlet.temperature[0].unfix()

    # IP_SH2 to HP_SH1
    m.fs.g13 = Arc(
        source=m.fs.IP_SH2.side_2_outlet,
        destination=m.fs.HP_SH1.side_2_inlet)
    m.fs.HP_SH1.side_2_inlet.flow_mol_comp[:, :].unfix()
    m.fs.HP_SH1.side_2_inlet.pressure[0].unfix()
    m.fs.HP_SH1.side_2_inlet.temperature[0].unfix()

    # HP_SH1 to HP_EVAP
    m.fs.g14 = Arc(
        source=m.fs.HP_SH1.side_2_outlet,
        destination=m.fs.HP_EVAP.shell_inlet)

    m.fs.HP_EVAP.shell_inlet.flow_mol_comp[:, :].unfix()
    m.fs.HP_EVAP.shell_inlet.pressure[0].unfix()
    m.fs.HP_EVAP.shell_inlet.temperature[0].unfix()
    # m.fs.HP_EVAP.overall_heat_transfer_coefficient.fix(150)

    # HP_EVAP to HP_ECON5
    m.fs.g15 = Arc(
        source=m.fs.HP_EVAP.shell_outlet,
        destination=m.fs.HP_ECON5.side_2_inlet)
    m.fs.HP_ECON5.side_2_inlet.flow_mol_comp[:, :].unfix()
    m.fs.HP_ECON5.side_2_inlet.pressure[0].unfix()
    m.fs.HP_ECON5.side_2_inlet.temperature[0].unfix()

    # HP_ECON5 to IP_SH1
    m.fs.g16 = Arc(
        source=m.fs.HP_ECON5.side_2_outlet,
        destination=m.fs.IP_SH1.side_2_inlet)
    m.fs.IP_SH1.side_2_inlet.flow_mol_comp[:, :].unfix()
    m.fs.IP_SH1.side_2_inlet.pressure[0].unfix()
    m.fs.IP_SH1.side_2_inlet.temperature[0].unfix()

    # IP_SH1 to HP_ECON4
    m.fs.g17 = Arc(
        source=m.fs.IP_SH1.side_2_outlet,
        destination=m.fs.HP_ECON4.side_2_inlet)
    m.fs.HP_ECON4.side_2_inlet.flow_mol_comp[:, :].unfix()
    m.fs.HP_ECON4.side_2_inlet.pressure[0].unfix()
    m.fs.HP_ECON4.side_2_inlet.temperature[0].unfix()

    # IP_SH1 to HP_ECON3
    m.fs.g18 = Arc(
        source=m.fs.HP_ECON4.side_2_outlet,
        destination=m.fs.HP_ECON3.side_2_inlet)
    m.fs.HP_ECON3.side_2_inlet.flow_mol_comp[:, :].unfix()
    m.fs.HP_ECON3.side_2_inlet.pressure[0].unfix()
    m.fs.HP_ECON3.side_2_inlet.temperature[0].unfix()

    # HP_ECON3 to LP_SH ** fix splitter before LP_SH, Mix after LP_SH
    m.fs.g19 = Arc(
        source=m.fs.HP_ECON3.side_2_outlet,
        destination=m.fs.LP_FGsplit.inlet)
    m.fs.LP_FGsplit.inlet.flow_mol_comp[:, :].unfix()
    m.fs.LP_FGsplit.inlet.pressure[0].unfix()
    m.fs.LP_FGsplit.inlet.temperature[0].unfix()

    m.fs.g20 = Arc(
        source=m.fs.LP_FGsplit.toLP_SH,
        destination=m.fs.LP_SH.side_2_inlet)
    m.fs.LP_SH.side_2_inlet.flow_mol_comp[:, :].unfix()
    m.fs.LP_SH.side_2_inlet.pressure[0].unfix()
    m.fs.LP_SH.side_2_inlet.temperature[0].unfix()

    m.fs.g21 = Arc(
        source=m.fs.LP_SH.side_2_outlet,
        destination=m.fs.LP_Mixer2.fromLP_SH)

    # splitter to mixer
    m.fs.g22 = Arc(
        source=m.fs.LP_FGsplit.toMixer,
        destination=m.fs.LP_Mixer2.bypass)
    m.fs.LP_FGsplit.inlet.flow_mol_comp[:, :].unfix()
    m.fs.LP_FGsplit.inlet.pressure[0].unfix()
    m.fs.LP_FGsplit.inlet.temperature[0].unfix()

    # LP_SH to IP_EVAP splitter and mixer to bypass
    m.fs.g23 = Arc(
        source=m.fs.LP_Mixer2.outlet,
        destination=m.fs.IP_EVAP.shell_inlet)
    m.fs.IP_EVAP.shell_inlet.flow_mol_comp[:, :].unfix()
    m.fs.IP_EVAP.shell_inlet.pressure[0].unfix()
    m.fs.IP_EVAP.shell_inlet.temperature[0].unfix()

    # IP_EVAP to IP_ECON2
    m.fs.g24 = Arc(
        source=m.fs.IP_EVAP.shell_outlet,
        destination=m.fs.IP_ECON2.side_2_inlet)
    m.fs.IP_ECON2.side_2_inlet.flow_mol_comp[:, :].unfix()
    m.fs.IP_ECON2.side_2_inlet.pressure[0].unfix()
    m.fs.IP_ECON2.side_2_inlet.temperature[0].unfix()

    # IP_ECON2 to HP_ECON2
    m.fs.g25 = Arc(
        source=m.fs.IP_ECON2.side_2_outlet,
        destination=m.fs.HP_ECON2.side_2_inlet)
    m.fs.HP_ECON2.side_2_inlet.flow_mol_comp[:, :].unfix()
    m.fs.HP_ECON2.side_2_inlet.pressure[0].unfix()
    m.fs.HP_ECON2.side_2_inlet.temperature[0].unfix()

    # HP_ECON2 to IP_ECON1
    m.fs.g26 = Arc(
        source=m.fs.HP_ECON2.side_2_outlet,
        destination=m.fs.IP_ECON1.side_2_inlet)
    m.fs.IP_ECON1.side_2_inlet.flow_mol_comp[:, :].unfix()
    m.fs.IP_ECON1.side_2_inlet.pressure[0].unfix()
    m.fs.IP_ECON1.side_2_inlet.temperature[0].unfix()

    # IP_ECON1 to HP_ECON1
    m.fs.g27 = Arc(
        source=m.fs.IP_ECON1.side_2_outlet,
        destination=m.fs.HP_ECON1.side_2_inlet)
    m.fs.HP_ECON1.side_2_inlet.flow_mol_comp[:, :].unfix()
    m.fs.HP_ECON1.side_2_inlet.pressure[0].unfix()
    m.fs.HP_ECON1.side_2_inlet.temperature[0].unfix()

    # # HP_ECON1 to LP_EVAP
    m.fs.g28 = Arc(
        source=m.fs.HP_ECON1.side_2_outlet,
        destination=m.fs.LP_EVAP.shell_inlet)
    m.fs.LP_EVAP.shell_inlet.flow_mol_comp[:, :].unfix()
    m.fs.LP_EVAP.shell_inlet.pressure[0].unfix()
    m.fs.LP_EVAP.shell_inlet.temperature[0].unfix()


    # LP_EVAP to LP_ECON1
    m.fs.g29 = Arc(
        source=m.fs.LP_EVAP.shell_outlet,
        destination=m.fs.LP_ECON.side_2_inlet)
    m.fs.LP_ECON.side_2_inlet.flow_mol_comp[:, :].unfix()
    m.fs.LP_ECON.side_2_inlet.pressure[0].unfix()
    m.fs.LP_ECON.side_2_inlet.temperature[0].unfix()

    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs)

    return m


def set_scaling_factors(m):
    """ Set scaling factors for variables and expressions. These are used for
    variable scaling and used by the framework to scale constraints.

    Args:
        m: plant model to set scaling factors for.

    Returns:
        None
    """
    hx_list = [m.fs.LP_ECON,
               m.fs.LP_SH,
               m.fs.IP_ECON1,
               m.fs.IP_ECON2,
               m.fs.IP_SH1,
               m.fs.IP_SH2,
               m.fs.IP_SH3,
               m.fs.HP_ECON1,
               m.fs.HP_ECON2,
               m.fs.HP_ECON3,
               m.fs.HP_ECON4,
               m.fs.HP_ECON5,
               m.fs.HP_SH1,
               m.fs.HP_SH2,
               m.fs.HP_SH3,
               m.fs.HP_SH4]
    for unit in hx_list:
        iscale.set_scaling_factor(unit.heat_transfer_correlation, 1e-8)
        iscale.set_scaling_factor(unit.side_1.heat, 1e-6)
        iscale.set_scaling_factor(unit.side_2.heat, 1e-6)
        iscale.set_scaling_factor(unit.side_1.volume, 1)
        iscale.set_scaling_factor(unit.side_2.volume, 1)


    fs = m.fs

    iscale.set_scaling_factor(fs.LP_EVAP.shell.heat, 1e-8)
    iscale.set_scaling_factor(fs.LP_EVAP.tube.heat, 1e-8)
    iscale.set_scaling_factor(fs.LP_EVAP.tube.heat, 1e-8)
    iscale.set_scaling_factor(fs.LP_EVAP.area, 1e-4)
    iscale.set_scaling_factor(fs.LP_EVAP.overall_heat_transfer_coefficient, 1e-2)

    iscale.set_scaling_factor(fs.IP_EVAP.shell.heat, 1e-8)
    iscale.set_scaling_factor(fs.IP_EVAP.tube.heat, 1e-8)
    iscale.set_scaling_factor(fs.IP_EVAP.tube.heat, 1e-8)
    iscale.set_scaling_factor(fs.IP_EVAP.area, 1e-4)
    iscale.set_scaling_factor(fs.IP_EVAP.overall_heat_transfer_coefficient, 1e-2)

    iscale.set_scaling_factor(fs.IP_pump.control_volume.work, 1e-6)
    iscale.set_scaling_factor(fs.HP_pump.control_volume.work, 1e-7)

    iscale.set_scaling_factor(fs.HP_EVAP.shell.heat, 1e-8)
    iscale.set_scaling_factor(fs.HP_EVAP.tube.heat, 1e-8)
    iscale.set_scaling_factor(fs.HP_EVAP.tube.heat, 1e-8)
    iscale.set_scaling_factor(fs.HP_EVAP.area, 1e-4)
    iscale.set_scaling_factor(fs.HP_EVAP.overall_heat_transfer_coefficient, 1e-2)

    # Calculate calculated scaling factors
    iscale.calculate_scaling_factors(m)


def tag_model(m):
    """
    Create to dictionaries. One is called tags with tags for keys and expressions
    for values. The second is called tag_format with tags for keys and a
    formatting string.  The formating string controls the number format and can
    be used to add additional text, to report units for example. The two
    dictionaries are returned, and also attached to the model.

    Args:
        m: (ConcreteModel) model to tag

    Returns:
        (dict) tags, (dict) tag_format
    """
    tags = {} # dict of with tag keys and expressions for their values
    tag_format = {} # format string for the tags
    def new_tag(name, expr, format):
        # funcion to keep it more compact
        tags[name] = expr
        tag_format[name] = format
    # Create a dict with Arc name keys and state block values
    stream_states = ta.stream_states_dict(
        ta.arcs_to_stream_dict(
            m.fs,
            descend_into=False,
            additional={
                "lp01": m.fs.LP_ECON.side_1_inlet,
                "lp03": m.fs.Mixer1.Preheater,
                "ip04": m.fs.IP_Splitter1.toNGPH,
                "lp11": m.fs.LP_SH.side_1_outlet,
                "hp11": m.fs.HP_SH4.side_1_outlet,
                "ip10": m.fs.IP_SH3.side_1_outlet,
                "ip11": m.fs.IP_Splitter2.inlet,
                "ip15": m.fs.IP_Splitter2.Cold_reheat,
                "ip13": m.fs.IP_Splitter2.toReclaimer,
                "ip12": m.fs.IP_Splitter2.toEjector,
                "ip14": m.fs.IP_Splitter2.toDryer,
                "g30": m.fs.LP_ECON.side_2_outlet,
                "g19": m.fs.HP_ECON3.side_2_outlet,}))
    for i, s in stream_states.items(): # create the tags for steam quantities
        new_tag(f"{i}_Fvol", expr=s.flow_vol, format="{:.1f} m^3/s")
        new_tag(f"{i}_Fmol", expr=s.flow_mol/1000, format="{:.3f} kmol/s")
        new_tag(f"{i}_F", expr=s.flow_mass, format="{:.3f} kg/s")
        new_tag(f"{i}_P", expr=s.pressure/1000, format="{:.2f} kPa")
        new_tag(f"{i}_T", expr=s.temperature, format="{:.2f} K")
        try:
            new_tag(f"{i}_H", expr=s.enth_mol/1000, format="{:.2f} kJ/mol")
        except:
            pass
        try:
            new_tag(f"{i}_X", expr=s.vapor_frac*100, format="{:.2f} %")
        except:
            pass
        try:
            for c in s.mole_frac_comp:
                new_tag(f"{i}_y{c}", expr=s.mole_frac_comp[c]*100, format="{:.3f}%")
        except:
            pass
    if hasattr(m, "tags"):
        m.tags.update(tags)
        m.tag_format.update(tag_format)
    else:
        m.tags = tags
        m.tag_format = tag_format
    return tags, tag_format


def pfd_result(outfile, m):
    # Add tags and data parameters
    template = os.path.join(this_file_dir(), "hrsg_template.svg")
    with open(template, "r") as f:
        svg_tag(
            svg=f,
            tags=m.tags,
            outfile=outfile,
            tag_format=m.tag_format)


def print_results(m):
    print('\nResults')
    print("\n\n Flue Gas Route ---------------")
    print('\n HP SH4')
    print('Flue Gas shell inlet temperature')
    m.fs.HP_SH4.side_2.properties_in[0].temperature.display()
    print('flue gas shell outlet temperature')
    m.fs.HP_SH4.side_2.properties_out[0].temperature.display()
    print('tube steam temperature')
    m.fs.HP_SH4.side_1.properties_in[0].temperature.display()
    m.fs.HP_SH4.side_1.properties_out[0].temperature.display()

    print('\n IP SH3')
    print('tube steam temperature')
    m.fs.IP_SH3.side_1.properties_in[0].temperature.display()
    m.fs.IP_SH3.side_1.properties_out[0].temperature.display()
    print('IP SH3 shell inlet/outlet temperature')
    m.fs.IP_SH3.side_2.properties_in[0].temperature.display()
    m.fs.IP_SH3.side_2.properties_out[0].temperature.display()

    print('\n HP_SH3')
    print('tube steam temperature')
    m.fs.HP_SH3.side_1.properties_in[0].temperature.display()
    m.fs.HP_SH3.side_1.properties_out[0].temperature.display()
    print('HP SH3 shell inlet/outlet temperature')
    m.fs.HP_SH3.side_2.properties_in[0].temperature.display()
    m.fs.HP_SH3.side_2.properties_out[0].temperature.display()

    print('\n HP SH2')
    print('tube steam temperature')
    m.fs.HP_SH2.side_1.properties_in[0].temperature.display()
    m.fs.HP_SH2.side_1.properties_out[0].temperature.display()
    print('HP SH2 shell inlet/outlet temperature')
    m.fs.HP_SH2.side_2.properties_in[0].temperature.display()
    m.fs.HP_SH2.side_2.properties_out[0].temperature.display()

    print('\n IP SH2')
    print('tube steam temperature')
    m.fs.IP_SH2.side_1.properties_in[0].temperature.display()
    m.fs.IP_SH2.side_1.properties_out[0].temperature.display()
    print('IP SH2 shell inlet/outlet temperature')
    m.fs.IP_SH2.side_2.properties_in[0].temperature.display()
    m.fs.IP_SH2.side_2.properties_out[0].temperature.display()

    print('\n HP SH1')
    print('tube steam temperature')
    m.fs.HP_SH1.side_1.properties_in[0].temperature.display()
    m.fs.HP_SH1.side_1.properties_out[0].temperature.display()
    print('HP SH1 shell inlet/outlet temperature')
    m.fs.HP_SH1.side_2.properties_in[0].temperature.display()
    m.fs.HP_SH1.side_2.properties_out[0].temperature.display()

    print('\n HP EVAP')
    print('tube steam temperature')
    m.fs.HP_EVAP.tube.properties_in[0].temperature.display()
    m.fs.HP_EVAP.tube.properties_out[0].temperature.display()
    print('HP EVAP shell inlet/outlet temperature')
    m.fs.HP_EVAP.shell.properties_in[0].temperature.display()
    m.fs.HP_EVAP.shell.properties_out[0].temperature.display()

    print('\n HP ECON5')
    print('tube steam temperature')
    m.fs.HP_ECON5.side_1.properties_in[0].temperature.display()
    m.fs.HP_ECON5.side_1.properties_out[0].temperature.display()
    print('HP ECON5 shell outlet temperature')
    m.fs.HP_ECON5.side_2.properties_in[0].temperature.display()
    m.fs.HP_ECON5.side_2.properties_out[0].temperature.display()

    print('\n IP SH1')
    print('tube steam temperature')
    m.fs.IP_SH1.side_1.properties_in[0].temperature.display()
    m.fs.IP_SH1.side_1.properties_out[0].temperature.display()
    print('IP SH2 shell inlet/outlet temperature')
    m.fs.IP_SH1.side_2.properties_in[0].temperature.display()
    m.fs.IP_SH1.side_2.properties_out[0].temperature.display()

    print('\n HP ECON4')
    print('tube steam temperature')
    m.fs.HP_ECON4.side_1.properties_in[0].temperature.display()
    m.fs.HP_ECON4.side_1.properties_out[0].temperature.display()
    print('HP ECON4 shell outlet temperature')
    m.fs.HP_ECON4.side_2.properties_in[0].temperature.display()
    m.fs.HP_ECON4.side_2.properties_out[0].temperature.display()

    print('\n HP ECON3')
    print('tube steam temperature')
    m.fs.HP_ECON3.side_1.properties_in[0].temperature.display()
    m.fs.HP_ECON3.side_1.properties_out[0].temperature.display()
    print('HP ECON3 shell outlet temperature')
    m.fs.HP_ECON3.side_2.properties_in[0].temperature.display()
    m.fs.HP_ECON3.side_2.properties_out[0].temperature.display()

    print('\n LP SH')
    print('tube steam temperature')
    m.fs.LP_SH.side_1.properties_in[0].temperature.display()
    m.fs.LP_SH.side_1.properties_out[0].temperature.display()
    print('IP SH2 shell inlet/outlet temperature')
    m.fs.LP_SH.side_2.properties_in[0].temperature.display()
    m.fs.LP_SH.side_2.properties_out[0].temperature.display()

    print('\n IP EVAP')
    print('tube steam temperature')
    m.fs.IP_EVAP.tube.properties_in[0].temperature.display()
    m.fs.IP_EVAP.tube.properties_out[0].temperature.display()
    print('HP EVAP shell inlet/outlet temperature')
    m.fs.IP_EVAP.shell.properties_in[0].temperature.display()
    m.fs.IP_EVAP.shell.properties_out[0].temperature.display()

    print('\n IP ECON2')
    print('tube steam temperature')
    m.fs.IP_ECON2.side_1.properties_in[0].temperature.display()
    m.fs.IP_ECON2.side_1.properties_out[0].temperature.display()
    print('IP ECON2 shell outlet temperature')
    m.fs.IP_ECON2.side_2.properties_in[0].temperature.display()
    m.fs.IP_ECON2.side_2.properties_out[0].temperature.display()

    print('\n HP ECON2')
    print('tube steam temperature')
    m.fs.HP_ECON2.side_1.properties_in[0].temperature.display()
    m.fs.HP_ECON2.side_1.properties_out[0].temperature.display()
    print('HP ECON2 shell outlet temperature')
    m.fs.HP_ECON2.side_2.properties_in[0].temperature.display()
    m.fs.HP_ECON2.side_2.properties_out[0].temperature.display()

    print('\n IP ECON1')
    print('tube steam temperature')
    m.fs.IP_ECON1.side_1.properties_in[0].temperature.display()
    m.fs.IP_ECON1.side_1.properties_out[0].temperature.display()
    print('IP ECON1 shell outlet temperature')
    m.fs.IP_ECON1.side_2.properties_in[0].temperature.display()
    m.fs.IP_ECON1.side_2.properties_out[0].temperature.display()

    print('\n HP ECON1')
    print('tube steam temperature')
    m.fs.HP_ECON1.side_1.properties_in[0].temperature.display()
    m.fs.HP_ECON1.side_1.properties_out[0].temperature.display()
    print('shell outlet temperature')
    m.fs.HP_ECON1.side_2.properties_in[0].temperature.display()
    m.fs.HP_ECON1.side_2.properties_out[0].temperature.display()

    print('\n LP EVAP')
    print('tube steam temperature')
    m.fs.LP_EVAP.tube.properties_in[0].temperature.display()
    m.fs.LP_EVAP.tube.properties_out[0].temperature.display()
    print('shell outlet temperature')
    m.fs.LP_EVAP.shell.properties_in[0].temperature.display()
    m.fs.LP_EVAP.shell.properties_out[0].temperature.display()

    print('\n LP ECON')
    print('tube steam temperature')
    m.fs.LP_ECON.side_1.properties_in[0].temperature.display()
    m.fs.LP_ECON.side_1.properties_out[0].temperature.display()
    print('shell outlet temperature')
    m.fs.LP_ECON.side_2.properties_in[0].temperature.display()
    m.fs.LP_ECON.side_2.properties_out[0].temperature.display()


def get_model(m=None, init=True):
    if m is None:
        m = pyo.ConcreteModel()
    if not hasattr(m, "fs"):
        m.fs = FlowsheetBlock(default={"dynamic": False})
    # add property packages to flowsheet library
    m.fs.prop_water = iapws95.Iapws95ParameterBlock()
    m.fs.prop_gas = FlueGasParameterBlock()

    # adding unit models
    m = add_unit_models(m)
    # fixing unit design variables
    m = set_inputs(m)
    set_scaling_factors(m)
    init_fname = 'hrsg_init.json.gz'
    if os.path.exists(init_fname):
        print('loading initial conditions')
        ms.from_json(m, fname=init_fname, wts=ms.StoreSpec(suffix=False))
        m = set_arcs(m)

    else:
        if init:
            # unit by unit initialization
            m = init_function(m)

        solver = pyo.SolverFactory("ipopt")
        print('degrees of freedom before arcs = ' + str(degrees_of_freedom(m)))
        # set arcs and unfix variables after initialization
        m = set_arcs(m)

        print('degrees of freedom = ' + str(degrees_of_freedom(m)))
        solver.options = {
            "tol": 1e-8,
            "max_iter": 20,
            # "halt_on_ampl_error": "yes",
        }
        # Solving full space model
        print('solving fullspace model')
        solver.solve(m, tee=True)
        print('saving init file')
        ms.to_json(m, fname=init_fname)
    tag_model(m)
    return m


def flow_sensitivity(m):
    # sensitivity analysis to test model robustness
    # reducing water flowrate to HRSG by 20%
    # a complete study should reduce flue gas flowrate similarly.
    steam_flow = []
    change = [1, 0.99, 0.95, 0.90, 0.87, 0.85, 0.80]
    for i in range(0, len(change)):
        m.fs.LP_ECON.side_1_inlet.flow_mol[0].fix(9790.55*change[i])  # mol/s
        FGrate = 38446.11  # mol/s
        # Use FG molar composition to set component flow rate (baseline report)
        m.fs.HP_SH4.side_2_inlet.flow_mol_comp[0, "H2O"].\
            fix(FGrate * change[i] * 0.0875)
        m.fs.HP_SH4.side_2_inlet.flow_mol_comp[0, "CO2"].\
            fix(FGrate * change[i] * 0.0408)
        m.fs.HP_SH4.side_2_inlet.flow_mol_comp[0, "N2"].\
            fix(FGrate * change[i] * 0.75)
        m.fs.HP_SH4.side_2_inlet.flow_mol_comp[0, "O2"].\
            fix(FGrate * change[i] * 0.12)
        m.fs.HP_SH4.side_2_inlet.flow_mol_comp[0, "NO"].\
            fix(FGrate * change[i] * 0.001)
        m.fs.HP_SH4.side_2_inlet.flow_mol_comp[0, "SO2"].\
            fix(FGrate * change[i] * 0.0007)
        solver.solve(m, tee=True)
        steam_flow.append(pyo.value(m.fs.LP_DRUM.vap_state[0].flow_mol))
        print(steam_flow)


if __name__ == "__main__":
    m = get_model()
    solver = pyo.SolverFactory("ipopt")
    solver.options = {
        "tol": 1e-8,
        "max_iter": 20,
    }
    set_scaling_factors(m)
    print('degrees of freedom =' + str(degrees_of_freedom(m)))
    res = solver.solve(m, tee=True, symbolic_solver_labels=True)

    print('\nFullspace Problem - Solver Termination Condition')
    print(res.solver.termination_condition)
    print('printing PFD results')
    print('HP steam conditions')
    m.fs.HP_SH4.side_1.properties_out[0].temperature.display()
    m.fs.HP_SH4.side_1_outlet.display()
    print('IP steam conditions')
    m.fs.IP_SH3.side_1.properties_out[0].temperature.display()
    m.fs.IP_SH3.side_1_outlet.display()
