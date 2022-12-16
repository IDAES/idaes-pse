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
# =============================================================================
"""
SCPC Boiler subflowsheet

Inputs:
    BFW - boiler feed water (from Feed water heaters)
    Coal from pulverizers

Main Assumptions:
    Coal flowrate as a function of load, coal HHV is fixed and heat dutty
    split from fire side to water wall and platen superheater is fixed.

    Boiler heat exchanger network:
        Water Flow:
            BFW -> ECONOMIZER -> Water Wall -> Primary SH -> Platen SH -> Finishing Superheate -> HP Turbine -> Reheater -> IP Turbine
        Flue Gas Flow:
            Fire Ball -> Platen SH -> Finishing SH -> Reheater  -> o -> Economizer -> Air Preheater
                                                   -> Primary SH --^

        * HP Turbine, IP Turbine, Air Preheater ==> not included in this release

    Models used:
        - Mixers: Attemperator, Flue gas mix
        - Heater: Platen SH, Fire/Water side (simplified model)
        - BoilerHeatExchanger: Economizer, Primary SH, Finishing SH, Reheater
            + Shell and tube heat exchanger
                - tube side: Steam (side 1 holdup)
                - shell side: flue gas (side 2 holdup)

    Property packages used:
        - IAPWS: Water/steam side
        - IDEAL GAS: Flue Gas side

Created: 1/10/2020 by Boiler subsystem team (M Zamarripa)

"""
__author__ = "Miguel Zamarripa"

# Import Pyomo libraries
from pyomo.environ import (
    ConcreteModel,
    value,
    TransformationFactory,
    units as pyunits,
)
from pyomo.network import Arc
from idaes.core.util.tags import svg_tag

# Import IDAES core
from idaes.core import FlowsheetBlock

# Import Unit Model Modules
from idaes.models.properties import iapws95

# Import Property Modules
from idaes.models_extra.power_generation.properties import FlueGasParameterBlock

# Import Unit Model Modules
from idaes.models.unit_models import Heater, Mixer
from idaes.models_extra.power_generation.unit_models.boiler_heat_exchanger import (
    BoilerHeatExchanger,
    TubeArrangement,
    HeatExchangerFlowPattern,
)
from idaes.models.unit_models.separator import (
    Separator,
    SplittingType,
    EnergySplittingType,
)
from pyomo.common.fileutils import this_file_dir
from collections import OrderedDict
import os
from idaes.core.util.model_statistics import degrees_of_freedom
import logging
from idaes.core.solvers import get_solver


def main():
    """
    Make the flowsheet object, fix some variables, and solve the problem
    """
    # Create a Concrete Model as the top level object
    m = ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.prop_water = iapws95.Iapws95ParameterBlock()

    build_boiler(m.fs)

    TransformationFactory("network.expand_arcs").apply_to(m)

    # Create a solver
    solver = get_solver()
    return (m, solver)


def build_boiler(fs):

    # Add property packages to flowsheet library
    fs.prop_fluegas = FlueGasParameterBlock()

    # Create unit models
    # Boiler Economizer
    fs.ECON = BoilerHeatExchanger(
        cold_side={"property_package": fs.prop_water, "has_pressure_change": True},
        hot_side={"property_package": fs.prop_fluegas, "has_pressure_change": True},
        has_holdup=False,
        flow_pattern=HeatExchangerFlowPattern.countercurrent,
        tube_arrangement=TubeArrangement.inLine,
        cold_side_water_phase="Liq",
        has_radiation=False,
    )
    # Primary Superheater
    fs.PrSH = BoilerHeatExchanger(
        cold_side={"property_package": fs.prop_water, "has_pressure_change": True},
        hot_side={"property_package": fs.prop_fluegas, "has_pressure_change": True},
        has_holdup=False,
        flow_pattern=HeatExchangerFlowPattern.countercurrent,
        tube_arrangement=TubeArrangement.inLine,
        cold_side_water_phase="Vap",
        has_radiation=True,
    )

    # Finishing Superheater
    fs.FSH = BoilerHeatExchanger(
        cold_side={"property_package": fs.prop_water, "has_pressure_change": True},
        hot_side={"property_package": fs.prop_fluegas, "has_pressure_change": True},
        has_holdup=False,
        flow_pattern=HeatExchangerFlowPattern.countercurrent,
        tube_arrangement=TubeArrangement.inLine,
        cold_side_water_phase="Vap",
        has_radiation=True,
    )

    # Reheater
    fs.RH = BoilerHeatExchanger(
        cold_side={"property_package": fs.prop_water, "has_pressure_change": True},
        hot_side={"property_package": fs.prop_fluegas, "has_pressure_change": True},
        has_holdup=False,
        flow_pattern=HeatExchangerFlowPattern.countercurrent,
        tube_arrangement=TubeArrangement.inLine,
        cold_side_water_phase="Vap",
        has_radiation=True,
    )
    # Platen Superheater
    fs.PlSH = Heater(property_package=fs.prop_water)

    # Boiler Water Wall
    fs.Water_wall = Heater(property_package=fs.prop_water)

    # Boiler Splitter (splits FSH flue gas outlet to Reheater and PrSH)
    fs.Spl1 = Separator(
        property_package=fs.prop_fluegas,
        split_basis=SplittingType.totalFlow,
        energy_split_basis=EnergySplittingType.equal_temperature,
    )
    # Flue gas mixer (mixing FG from Reheater and Primary SH, inlet to ECON)
    fs.mix1 = Mixer(
        property_package=fs.prop_fluegas,
        inlet_list=["Reheat_out", "PrSH_out"],
        dynamic=False,
    )

    # Mixer for Attemperator #1 (between PrSH and PlSH)
    fs.ATMP1 = Mixer(
        property_package=fs.prop_water,
        inlet_list=["Steam", "SprayWater"],
        dynamic=False,
    )

    # Build connections (streams)

    # Steam Route (side 1 = tube side = steam/water side)
    # Boiler feed water to Economizer (to be imported in full plant model)
    #    fs.bfw2econ = Arc(source=fs.FWH8.outlet,
    #                           destination=fs.ECON.cold_side_inlet)
    fs.econ2ww = Arc(source=fs.ECON.cold_side_outlet, destination=fs.Water_wall.inlet)
    fs.ww2prsh = Arc(source=fs.Water_wall.outlet, destination=fs.PrSH.cold_side_inlet)
    fs.prsh2plsh = Arc(source=fs.PrSH.cold_side_outlet, destination=fs.PlSH.inlet)
    fs.plsh2fsh = Arc(source=fs.PlSH.outlet, destination=fs.FSH.cold_side_inlet)
    fs.FSHtoATMP1 = Arc(source=fs.FSH.cold_side_outlet, destination=fs.ATMP1.Steam)
    #    fs.fsh2hpturbine=Arc(source=fs.ATMP1.outlet,
    #                           destination=fs.HPTinlet)
    # (to be imported in full plant model)

    # Flue gas route ---------------------------------------------------------
    # water wall connected with boiler block (to fix the heat duty)
    # platen SH connected with boiler block (to fix the heat duty)
    # Finishing superheater connected with a flowsheet level constraint
    fs.fg_fsh2_separator = Arc(source=fs.FSH.hot_side_outlet, destination=fs.Spl1.inlet)
    fs.fg_fsh2rh = Arc(source=fs.Spl1.outlet_1, destination=fs.RH.hot_side_inlet)
    fs.fg_fsh2PrSH = Arc(source=fs.Spl1.outlet_2, destination=fs.PrSH.hot_side_inlet)
    fs.fg_rhtomix = Arc(source=fs.RH.hot_side_outlet, destination=fs.mix1.Reheat_out)
    fs.fg_prsh2mix = Arc(source=fs.PrSH.hot_side_outlet, destination=fs.mix1.PrSH_out)
    fs.fg_mix2econ = Arc(source=fs.mix1.outlet, destination=fs.ECON.hot_side_inlet)


# Set inputs ==========================
def initialize(m):
    # ------------- ECONOMIZER -----------------------------------------------
    # BFW Boiler Feed Water inlet temeperature = 555 F = 563.706 K
    m.fs.ECON.cold_side_inlet.flow_mol[0].fix(24194.177)  # mol/s
    m.fs.ECON.cold_side_inlet.enth_mol[0].fix(14990.97)
    m.fs.ECON.cold_side_inlet.pressure[0].fix(26922222.222)  # Pa

    # FLUE GAS Inlet from Primary Superheater
    FGrate = 21290.6999  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.ECON.hot_side_inlet.flow_mol_comp[0, "H2O"].fix(FGrate * 8.69 / 100)
    m.fs.ECON.hot_side_inlet.flow_mol_comp[0, "CO2"].fix(FGrate * 14.49 / 100)
    m.fs.ECON.hot_side_inlet.flow_mol_comp[0, "N2"].fix(FGrate * 74.34 / 100)
    m.fs.ECON.hot_side_inlet.flow_mol_comp[0, "O2"].fix(FGrate * 2.47 / 100)
    m.fs.ECON.hot_side_inlet.flow_mol_comp[0, "NO"].fix(FGrate * 0.0006)
    m.fs.ECON.hot_side_inlet.flow_mol_comp[0, "SO2"].fix(FGrate * 0.002)
    m.fs.ECON.hot_side_inlet.temperature[0].fix(682.335)  # K
    m.fs.ECON.hot_side_inlet.pressure[0].fix(100145)  # Pa

    # economizer design variables and parameters
    ITM = 0.0254  # inch to meter conversion
    # Based on NETL Baseline Report Rev3
    m.fs.ECON.tube_di.fix((2 - 2 * 0.188) * ITM)  # calc inner diameter
    #                                   (2 = outer diameter, thickness = 0.188)
    m.fs.ECON.tube_thickness.fix(0.188 * ITM)  # tube thickness
    m.fs.ECON.pitch_x.fix(3.5 * ITM)
    # pitch_y = (54.5) gas path transverse width /columns
    m.fs.ECON.pitch_y.fix(5.03 * ITM)
    m.fs.ECON.tube_length.fix(53.41 * 12 * ITM)  # use tube length (53.41 ft)
    m.fs.ECON.tube_nrow.fix(36 * 2.5)  # use to match baseline performance
    m.fs.ECON.tube_ncol.fix(130)  # 130 from NETL
    m.fs.ECON.nrow_inlet.fix(2)
    m.fs.ECON.delta_elevation.fix(50)
    # parameters
    # heat transfer resistance due to tube side fouling (water scales)
    m.fs.ECON.tube_r_fouling = 0.000176
    # heat transfer resistance due to tube shell fouling (ash deposition)
    m.fs.ECON.shell_r_fouling = 0.00088
    if m.fs.ECON.config.has_radiation is True:
        m.fs.ECON.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.ECON.fcorrection_htc.fix(1.5)
    # correction factor for pressure drop calc tube side
    m.fs.ECON.fcorrection_dp_tube.fix(1.0)
    # correction factor for pressure drop calc shell side
    m.fs.ECON.fcorrection_dp_shell.fix(1.0)

    # --------- Primary Superheater ------------
    # Steam from water wall
    hprsh = value(iapws95.htpx(773.15 * pyunits.K, 24865516.722 * pyunits.Pa))
    print(hprsh)
    m.fs.PrSH.cold_side_inlet.flow_mol[0].fix(24190.26)  # mol/s
    m.fs.PrSH.cold_side_inlet.enth_mol[0].fix(hprsh)  # J/mol
    m.fs.PrSH.cold_side_inlet.pressure[0].fix(2.5449e7)  # Pascals

    # FLUE GAS Inlet from Primary Superheater
    FGrate = 21290.6999 * 0.18  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.PrSH.hot_side_inlet.flow_mol_comp[0, "H2O"].fix(FGrate * 8.69 / 100)
    m.fs.PrSH.hot_side_inlet.flow_mol_comp[0, "CO2"].fix(FGrate * 14.49 / 100)
    m.fs.PrSH.hot_side_inlet.flow_mol_comp[0, "N2"].fix(FGrate * 74.34 / 100)
    m.fs.PrSH.hot_side_inlet.flow_mol_comp[0, "O2"].fix(FGrate * 2.47 / 100)
    m.fs.PrSH.hot_side_inlet.flow_mol_comp[0, "NO"].fix(FGrate * 0.0006)
    m.fs.PrSH.hot_side_inlet.flow_mol_comp[0, "SO2"].fix(FGrate * 0.002)
    m.fs.PrSH.hot_side_inlet.temperature[0].fix(1180.335)
    m.fs.PrSH.hot_side_inlet.pressure[0].fix(100145)

    # Primary Superheater
    ITM = 0.0254  # inch to meter conversion
    m.fs.PrSH.tube_di.fix((2.5 - 2 * 0.165) * ITM)
    m.fs.PrSH.tube_thickness.fix(0.165 * ITM)
    m.fs.PrSH.pitch_x.fix(3 * ITM)
    # gas path transverse width 54.78 ft / number of columns
    m.fs.PrSH.pitch_y.fix(54.78 / 108 * 12 * ITM)
    m.fs.PrSH.tube_length.fix(53.13 * 12 * ITM)
    m.fs.PrSH.tube_nrow.fix(20 * 2)
    m.fs.PrSH.tube_ncol.fix(108)
    m.fs.PrSH.nrow_inlet.fix(4)
    m.fs.PrSH.delta_elevation.fix(50)
    m.fs.PrSH.tube_r_fouling = 0.000176  # (0.001 h-ft^2-F/BTU)
    m.fs.PrSH.shell_r_fouling = 0.003131  # (0.03131 - 0.1779 h-ft^2-F/BTU)
    if m.fs.PrSH.config.has_radiation is True:
        m.fs.PrSH.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.PrSH.fcorrection_htc.fix(1.5)
    # correction factor for pressure drop calc tube side
    m.fs.PrSH.fcorrection_dp_tube.fix(1.0)
    # correction factor for pressure drop calc shell side
    m.fs.PrSH.fcorrection_dp_shell.fix(1.0)

    #  ----- Finishing Superheater ----------------------
    # Steam from Platen Supeheater
    hfsh = value(iapws95.htpx(823.15 * pyunits.K, 24790249.01 * pyunits.Pa))
    m.fs.FSH.cold_side_inlet.flow_mol[0].fix(24194.177)  # mol/s
    m.fs.FSH.cold_side_inlet.enth_mol[0].fix(hfsh)  # J/mol
    m.fs.FSH.cold_side_inlet.pressure[0].fix(24790249.01)  # Pascals

    # FLUE GAS Inlet from Primary Superheater
    FGrate = 21290.6999  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.FSH.hot_side_inlet.flow_mol_comp[0, "H2O"].fix(FGrate * 8.69 / 100)
    m.fs.FSH.hot_side_inlet.flow_mol_comp[0, "CO2"].fix(FGrate * 14.49 / 100)
    m.fs.FSH.hot_side_inlet.flow_mol_comp[0, "N2"].fix(FGrate * 74.34 / 100)
    m.fs.FSH.hot_side_inlet.flow_mol_comp[0, "O2"].fix(FGrate * 2.47 / 100)
    m.fs.FSH.hot_side_inlet.flow_mol_comp[0, "NO"].fix(FGrate * 0.0006)
    m.fs.FSH.hot_side_inlet.flow_mol_comp[0, "SO2"].fix(FGrate * 0.002)
    m.fs.FSH.hot_side_inlet.temperature[0].fix(1300.335)
    m.fs.FSH.hot_side_inlet.pressure[0].fix(100145)

    # Finishing Superheater
    ITM = 0.0254  # inch to meter conversion
    m.fs.FSH.tube_di.fix((2.5 - 2 * 0.165) * ITM)
    m.fs.FSH.tube_thickness.fix(0.165 * ITM)
    m.fs.FSH.pitch_x.fix(3 * ITM)
    # gas path transverse width 54.78 ft / number of columns
    m.fs.FSH.pitch_y.fix(54.78 / 108 * 12 * ITM)
    m.fs.FSH.tube_length.fix(53.13 * 12 * ITM)
    m.fs.FSH.tube_nrow.fix(8 * 2)
    m.fs.FSH.tube_ncol.fix(85)
    m.fs.FSH.nrow_inlet.fix(2)
    m.fs.FSH.delta_elevation.fix(50)
    m.fs.FSH.tube_r_fouling = 0.000176  # (0.001 h-ft^2-F/BTU)
    m.fs.FSH.shell_r_fouling = 0.003131  # (0.03131 - 0.1779 h-ft^2-F/BTU)
    if m.fs.FSH.config.has_radiation is True:
        m.fs.FSH.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.FSH.fcorrection_htc.fix(1.0)
    # correction factor for pressure drop calc tube side
    m.fs.FSH.fcorrection_dp_tube.fix(1.0)
    # correction factor for pressure drop calc shell side
    m.fs.FSH.fcorrection_dp_shell.fix(1.0)

    #  ----- Reheater Superheater ----------------------
    #   Steam from HP Turbine outlet
    m.fs.RH.cold_side_inlet.flow_mol[0].fix(21235.27)  # mol/s
    m.fs.RH.cold_side_inlet.enth_mol[0].fix(53942.7569)  # J/mol
    m.fs.RH.cold_side_inlet.pressure[0].fix(3677172.33638)  # Pascals

    # FLUE GAS Inlet from Finishing Superheater
    FGrate = 21290.6999 * 0.85  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.RH.hot_side_inlet.flow_mol_comp[0, "H2O"].fix(FGrate * 8.69 / 100)
    m.fs.RH.hot_side_inlet.flow_mol_comp[0, "CO2"].fix(FGrate * 14.49 / 100)
    m.fs.RH.hot_side_inlet.flow_mol_comp[0, "N2"].fix(FGrate * 74.34 / 100)
    m.fs.RH.hot_side_inlet.flow_mol_comp[0, "O2"].fix(FGrate * 2.47 / 100)
    m.fs.RH.hot_side_inlet.flow_mol_comp[0, "NO"].fix(FGrate * 0.0006)
    m.fs.RH.hot_side_inlet.flow_mol_comp[0, "SO2"].fix(FGrate * 0.002)
    m.fs.RH.hot_side_inlet.temperature[0].fix(1180.335)
    m.fs.RH.hot_side_inlet.pressure[0].fix(100145)

    # Reheater Superheater
    ITM = 0.0254  # inch to meter conversion
    m.fs.RH.tube_di.fix((2.5 - 2 * 0.165) * ITM)
    m.fs.RH.tube_thickness.fix(0.11 * ITM)
    m.fs.RH.pitch_x.fix(3 * ITM)
    # gas path transverse width 54.08 ft / number of columns
    m.fs.RH.pitch_y.fix(54.08 / 108 * 12 * ITM)
    m.fs.RH.tube_length.fix(53.82 * 12 * ITM)
    m.fs.RH.tube_nrow.fix(18 * 2)
    m.fs.RH.tube_ncol.fix(82 + 70)
    m.fs.RH.nrow_inlet.fix(2)
    m.fs.RH.delta_elevation.fix(30)
    m.fs.RH.tube_r_fouling = 0.000176  # (0.001 h-ft^2-F/BTU)
    m.fs.RH.shell_r_fouling = 0.00088  # (0.03131 - 0.1779 h-ft^2-F/BTU)
    if m.fs.RH.config.has_radiation is True:
        m.fs.RH.emissivity_wall.fix(0.7)  # wall emissivity
    # correction factor for overall heat transfer coefficient
    m.fs.RH.fcorrection_htc.fix(1.8)
    # correction factor for pressure drop calc tube side
    m.fs.RH.fcorrection_dp_tube.fix(1.0)
    # correction factor for pressure drop calc shell side
    m.fs.RH.fcorrection_dp_shell.fix(1.0)

    #    Platen Superheater ------------------------------------------
    hpl = value(iapws95.htpx(798.15 * pyunits.K, 24790249.01 * pyunits.Pa))
    m.fs.PlSH.inlet[:].flow_mol.fix(24194.177)
    m.fs.PlSH.inlet[:].enth_mol.fix(hpl)
    m.fs.PlSH.inlet[:].pressure.fix(24790249.01)
    m.fs.PlSH.heat_duty[:].fix(5.5e7)

    #    Water wall Superheater ------------------------------------------
    hww = value(iapws95.htpx(588.15 * pyunits.K, 2.5449e7 * pyunits.Pa))
    m.fs.Water_wall.inlet[:].flow_mol.fix(24194.177)
    m.fs.Water_wall.inlet[:].enth_mol.fix(hww)
    m.fs.Water_wall.inlet[:].pressure.fix(24865516.722)
    m.fs.Water_wall.heat_duty[:].fix(7.51e8)  # 8.76e8

    #   splitter flue gas from Finishing SH to Reheater and Primary SH
    m.fs.Spl1.split_fraction[0, "outlet_1"].fix(0.75)  # 0.85)
    # FLUE GAS Inlet from Primary Superheater
    FGrate = 21290.6999  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.Spl1.inlet.flow_mol_comp[0, "H2O"].fix(FGrate * 8.69 / 100)
    m.fs.Spl1.inlet.flow_mol_comp[0, "CO2"].fix(FGrate * 14.49 / 100)
    m.fs.Spl1.inlet.flow_mol_comp[0, "N2"].fix(FGrate * 74.34 / 100)
    m.fs.Spl1.inlet.flow_mol_comp[0, "O2"].fix(FGrate * 2.47 / 100)
    m.fs.Spl1.inlet.flow_mol_comp[0, "NO"].fix(FGrate * 0.0006)
    m.fs.Spl1.inlet.flow_mol_comp[0, "SO2"].fix(FGrate * 0.002)
    m.fs.Spl1.inlet.temperature[0].fix(1180.335)
    m.fs.Spl1.inlet.pressure[0].fix(100145)

    # mixer (econ inlet) fluegas outlet from reheater and primary superheater
    #   Mixer inlets  ['Reheat_out', 'PrSH_out']
    m.fs.mix1.Reheat_out.flow_mol_comp[0, "H2O"].fix(FGrate * 8.69 / 100)
    m.fs.mix1.Reheat_out.flow_mol_comp[0, "CO2"].fix(FGrate * 14.49 / 100)
    m.fs.mix1.Reheat_out.flow_mol_comp[0, "N2"].fix(FGrate * 74.34 / 100)
    m.fs.mix1.Reheat_out.flow_mol_comp[0, "O2"].fix(FGrate * 2.47 / 100)
    m.fs.mix1.Reheat_out.flow_mol_comp[0, "NO"].fix(FGrate * 0.0006)
    m.fs.mix1.Reheat_out.flow_mol_comp[0, "SO2"].fix(FGrate * 0.002)
    m.fs.mix1.Reheat_out.pressure.fix(100145)
    m.fs.mix1.Reheat_out.temperature.fix(731.5)
    # Fixed SprayWater as small flow and arbitrary conditions
    # (will be connected from FW pump splitter)
    m.fs.mix1.PrSH_out.flow_mol_comp[0, "H2O"].fix(FGrate * 8.69 / 100)
    m.fs.mix1.PrSH_out.flow_mol_comp[0, "CO2"].fix(FGrate * 14.49 / 100)
    m.fs.mix1.PrSH_out.flow_mol_comp[0, "N2"].fix(FGrate * 74.34 / 100)
    m.fs.mix1.PrSH_out.flow_mol_comp[0, "O2"].fix(FGrate * 2.47 / 100)
    m.fs.mix1.PrSH_out.flow_mol_comp[0, "NO"].fix(FGrate * 0.0006)
    m.fs.mix1.PrSH_out.flow_mol_comp[0, "SO2"].fix(FGrate * 0.002)
    m.fs.mix1.PrSH_out.pressure.fix(100145)
    m.fs.mix1.PrSH_out.temperature.fix(731.15)

    # Attemperator inputs
    hatt = value(iapws95.htpx(875.15 * pyunits.K, 24865516.722 * pyunits.Pa))
    m.fs.ATMP1.Steam.flow_mol.fix(24194.177)
    m.fs.ATMP1.Steam.enth_mol.fix(hatt)
    m.fs.ATMP1.Steam.pressure.fix(24865516.722)
    # Fixed SprayWater from FW pump splitter (splitter is needded)
    hatt2 = value(iapws95.htpx(563.15 * pyunits.K, 2.5449e7 * pyunits.Pa))
    m.fs.ATMP1.SprayWater_state[:].flow_mol.fix(0.001)
    m.fs.ATMP1.SprayWater_state[:].pressure.fix(1.22e8)
    m.fs.ATMP1.SprayWater_state[:].enth_mol.fix(hatt2)

    # ---------  Initialization ------------------------------------------

    # Initialize Units
    m.fs.ECON.initialize(outlvl=logging.INFO)
    m.fs.PrSH.initialize(outlvl=logging.INFO)
    m.fs.FSH.initialize(outlvl=logging.INFO)
    m.fs.RH.initialize(outlvl=logging.INFO)
    m.fs.PlSH.initialize(outlvl=logging.INFO)
    m.fs.Water_wall.initialize(outlvl=logging.INFO)
    m.fs.Spl1.initialize(outlvl=logging.INFO)
    m.fs.mix1.initialize(outlvl=logging.INFO)
    m.fs.ATMP1.initialize(outlvl=logging.INFO)
    print("initialization done")


def unfix_inlets(m):

    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.ECON.hot_side_inlet.flow_mol_comp[0, "H2O"].unfix()
    m.fs.ECON.hot_side_inlet.flow_mol_comp[0, "CO2"].unfix()
    m.fs.ECON.hot_side_inlet.flow_mol_comp[0, "N2"].unfix()
    m.fs.ECON.hot_side_inlet.flow_mol_comp[0, "O2"].unfix()
    m.fs.ECON.hot_side_inlet.flow_mol_comp[0, "NO"].unfix()
    m.fs.ECON.hot_side_inlet.flow_mol_comp[0, "SO2"].unfix()
    m.fs.ECON.hot_side_inlet.temperature[0].unfix()
    m.fs.ECON.hot_side_inlet.pressure[0].unfix()

    # PrSH Primary superheater inlets (steam and flue gas) --------------------
    m.fs.PrSH.cold_side_inlet.flow_mol.unfix()
    m.fs.PrSH.cold_side_inlet.enth_mol[0].unfix()
    m.fs.PrSH.cold_side_inlet.pressure[0].unfix()

    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.PrSH.hot_side_inlet.flow_mol_comp[0, "H2O"].unfix()
    m.fs.PrSH.hot_side_inlet.flow_mol_comp[0, "CO2"].unfix()
    m.fs.PrSH.hot_side_inlet.flow_mol_comp[0, "N2"].unfix()
    m.fs.PrSH.hot_side_inlet.flow_mol_comp[0, "O2"].unfix()
    m.fs.PrSH.hot_side_inlet.flow_mol_comp[0, "NO"].unfix()
    m.fs.PrSH.hot_side_inlet.flow_mol_comp[0, "SO2"].unfix()
    m.fs.PrSH.hot_side_inlet.temperature[0].unfix()
    m.fs.PrSH.hot_side_inlet.pressure[0].unfix()

    # WaterWall  water from economizer  ---------------------------------------
    m.fs.Water_wall.inlet[:].flow_mol.unfix()
    m.fs.Water_wall.inlet[:].enth_mol.unfix()
    m.fs.Water_wall.inlet[:].pressure.unfix()

    # Use FG molar composition to set component flow rates (baseline report)
    m.fs.RH.hot_side_inlet.flow_mol_comp[0, "H2O"].unfix()
    m.fs.RH.hot_side_inlet.flow_mol_comp[0, "CO2"].unfix()
    m.fs.RH.hot_side_inlet.flow_mol_comp[0, "N2"].unfix()
    m.fs.RH.hot_side_inlet.flow_mol_comp[0, "O2"].unfix()
    m.fs.RH.hot_side_inlet.flow_mol_comp[0, "NO"].unfix()
    m.fs.RH.hot_side_inlet.flow_mol_comp[0, "SO2"].unfix()
    m.fs.RH.hot_side_inlet.temperature[0].unfix()
    m.fs.RH.hot_side_inlet.pressure[0].unfix()

    # Finishing Superheater (steam from Platen SH)-----------------------------
    m.fs.FSH.cold_side_inlet.flow_mol.unfix()
    m.fs.FSH.cold_side_inlet.enth_mol[0].unfix()
    m.fs.FSH.cold_side_inlet.pressure[0].unfix()

    # Platen SH steam inlet conditions from Water Wall-------------------------
    m.fs.PlSH.inlet[:].flow_mol.unfix()
    m.fs.PlSH.inlet[:].enth_mol.unfix()
    m.fs.PlSH.inlet[:].pressure.unfix()

    # unfix Splitter flue gas inlet to RH and PrSH
    m.fs.Spl1.inlet.flow_mol_comp[0, "H2O"].unfix()
    m.fs.Spl1.inlet.flow_mol_comp[0, "CO2"].unfix()
    m.fs.Spl1.inlet.flow_mol_comp[0, "N2"].unfix()
    m.fs.Spl1.inlet.flow_mol_comp[0, "O2"].unfix()
    m.fs.Spl1.inlet.flow_mol_comp[0, "NO"].unfix()
    m.fs.Spl1.inlet.flow_mol_comp[0, "SO2"].unfix()
    m.fs.Spl1.inlet.temperature[0].unfix()
    m.fs.Spl1.inlet.pressure[0].unfix()

    # unfix mixer inlets (now connected)
    m.fs.mix1.Reheat_out.flow_mol_comp[0, "H2O"].unfix()
    m.fs.mix1.Reheat_out.flow_mol_comp[0, "CO2"].unfix()
    m.fs.mix1.Reheat_out.flow_mol_comp[0, "N2"].unfix()
    m.fs.mix1.Reheat_out.flow_mol_comp[0, "O2"].unfix()
    m.fs.mix1.Reheat_out.flow_mol_comp[0, "NO"].unfix()
    m.fs.mix1.Reheat_out.flow_mol_comp[0, "SO2"].unfix()
    m.fs.mix1.Reheat_out.pressure.unfix()
    m.fs.mix1.Reheat_out.temperature.unfix()
    # PrSH output
    m.fs.mix1.PrSH_out.flow_mol_comp[0, "H2O"].unfix()
    m.fs.mix1.PrSH_out.flow_mol_comp[0, "CO2"].unfix()
    m.fs.mix1.PrSH_out.flow_mol_comp[0, "N2"].unfix()
    m.fs.mix1.PrSH_out.flow_mol_comp[0, "O2"].unfix()
    m.fs.mix1.PrSH_out.flow_mol_comp[0, "NO"].unfix()
    m.fs.mix1.PrSH_out.flow_mol_comp[0, "SO2"].unfix()
    m.fs.mix1.PrSH_out.pressure.unfix()
    m.fs.mix1.PrSH_out.temperature.unfix()

    # unfix attemperator inlet from (Finishing SH)
    m.fs.ATMP1.Steam.flow_mol.unfix()
    m.fs.ATMP1.Steam.enth_mol.unfix()
    m.fs.ATMP1.Steam.pressure.unfix()


def pfd_result(outfile, m, df):
    tags = {}
    for i in df.index:
        tags[i + "_F"] = df.loc[i, "Molar Flow"]
        tags[i + "_T"] = df.loc[i, "T"]
        tags[i + "_P"] = df.loc[i, "P"]
        tags[i + "_X"] = df.loc[i, "Vapor Fraction"]

    tags["FG_2_RH_Fm"] = value(m.fs.RH.side_2.properties_in[0].flow_mass)
    tags["FG_2_RH_T"] = value(m.fs.RH.side_2.properties_in[0].temperature)
    tags["FG_2_RH_P"] = value(m.fs.RH.side_2.properties_in[0].pressure)

    tags["FG_RH_2_Mix_Fm"] = value(m.fs.RH.side_2.properties_out[0].flow_mass)
    tags["FG_RH_2_Mix_T"] = value(m.fs.RH.side_2.properties_out[0].temperature)
    tags["FG_RH_2_Mix_P"] = value(m.fs.RH.side_2.properties_out[0].pressure)

    tags["FG_2_FSH_Fm"] = value(m.fs.FSH.side_2.properties_in[0].flow_mass)
    tags["FG_2_FSH_T"] = value(m.fs.FSH.side_2.properties_in[0].temperature)
    tags["FG_2_FSH_P"] = value(m.fs.FSH.side_2.properties_in[0].pressure)

    tags["FG_2_PrSH_Fm"] = value(m.fs.PrSH.side_2.properties_in[0].flow_mass)
    tags["FG_2_PrSH_T"] = value(m.fs.PrSH.side_2.properties_in[0].temperature)
    tags["FG_2_PrSH_P"] = value(m.fs.PrSH.side_2.properties_in[0].pressure)

    tags["FG_PrSH_2_Mix_Fm"] = value(m.fs.PrSH.side_2.properties_out[0].flow_mass)
    tags["FG_PrSH_2_Mix_T"] = value(m.fs.PrSH.side_2.properties_out[0].temperature)
    tags["FG_PrSH_2_Mix_P"] = value(m.fs.PrSH.side_2.properties_out[0].pressure)

    tags["FG_2_ECON_Fm"] = value(m.fs.ECON.side_2.properties_in[0].flow_mass)
    tags["FG_2_ECON_T"] = value(m.fs.ECON.side_2.properties_in[0].temperature)
    tags["FG_2_ECON_P"] = value(m.fs.ECON.side_2.properties_in[0].pressure)

    tags["FG_2_AIRPH_Fm"] = value(m.fs.ECON.side_2.properties_out[0].flow_mass)
    tags["FG_2_AIRPH_T"] = value(m.fs.ECON.side_2.properties_out[0].temperature)
    tags["FG_2_AIRPH_P"] = value(m.fs.ECON.side_2.properties_out[0].pressure)

    tags["FG_2_STACK_Fm"] = value(m.fs.ECON.side_2.properties_out[0].flow_mass)
    tags["FG_2_STACK_T"] = value(m.fs.ECON.side_2.properties_out[0].temperature)
    tags["FG_2_STACK_P"] = value(m.fs.ECON.side_2.properties_out[0].pressure)

    original_svg_file = os.path.join(this_file_dir(), "Boiler_scpc_PFD.svg")
    with open(original_svg_file, "r") as f:
        s = svg_tag(tags, f, outfile=outfile)


def _stream_dict(m):
    """Adds _streams to m, which contains a dictionary of streams for display

    Args:
        m (ConcreteModel): A Pyomo model from create_model()

    Returns:
        None
    """

    m._streams = OrderedDict(
        [
            ("MS", m.fs.ATMP1.mixed_state),
            ("ATMP_In", m.fs.FSH.side_1.properties_out),
            ("FSH_In", m.fs.FSH.side_1.properties_in),
            ("PrSH_IN", m.fs.PrSH.side_1.properties_in),
            ("RHT_COLD", m.fs.RH.side_1.properties_in),
            ("RHT_HOT", m.fs.RH.side_1.properties_out),
            ("PlatenSH_IN", m.fs.PlSH.control_volume.properties_in),
            ("BFW", m.fs.ECON.side_1.properties_in),
            ("ECON_OUT", m.fs.ECON.side_1.properties_out),
        ]
    )


def print_results(m):
    print()
    print("Results")
    print()

    print("viscosity gas side = ", m.fs.PrSH.side_2.properties_in[0].visc_d.value)
    print(
        "conductivity gas side = ", m.fs.PrSH.side_2.properties_in[0].therm_cond.value
    )
    print("velocity_tube = ", m.fs.PrSH.v_tube[0].value)
    print("velocity_shell = ", m.fs.PrSH.v_shell[0].value)
    print("Re_tube = ", m.fs.PrSH.N_Re_tube[0].value)
    print("Re_shell = ", m.fs.PrSH.N_Re_shell[0].value)
    print("hconv_tube = ", m.fs.PrSH.hconv_tube[0].value)
    print("hconv_shell_rad = ", m.fs.PrSH.hconv_shell_rad[0].value)
    print("hconv_shell_conv = ", m.fs.PrSH.hconv_shell_conv[0].value)
    print("hconv_shell_total = ", m.fs.PrSH.hconv_shell_total[0].value)
    print("driving force = ", m.fs.PrSH.temperature_driving_force[0].value)
    print("dT_inlet = ", m.fs.PrSH.deltaT_1[0].value)
    print("dT_outlet = ", m.fs.PrSH.deltaT_2[0].value)
    print("deltaP tube = ", m.fs.PrSH.deltaP_tube[0].value)
    print("deltaP shell = ", m.fs.PrSH.deltaP_shell[0].value)
    print("mbl = ", value(m.fs.PrSH.mbl))

    if m.fs.PrSH.config.has_radiation is True:
        print("gas emissivity = ", m.fs.PrSH.gas_emissivity[0].value)
        print("gas emissivity div2 = ", m.fs.PrSH.gas_emissivity_div2[0].value)
        print("gas emissivity mul2 = ", m.fs.PrSH.gas_emissivity_mul2[0].value)
        print("gas gray fraction = ", m.fs.PrSH.gas_gray_fraction[0].value)
    print(
        "liquid density in = ",
        value(m.fs.PrSH.side_1.properties_in[0].dens_mass_phase["Liq"]),
    )

    print(
        "liquid density out = ",
        value(m.fs.PrSH.side_1.properties_out[0].dens_mass_phase["Liq"]),
    )
    print("heat transfer area = ", value(m.fs.PrSH.area_heat_transfer))
    print(
        "overal heat transfer = ", value(m.fs.PrSH.overall_heat_transfer_coefficient[0])
    )

    print("\n\n ------------- Economizer   ---------")
    print("liquid temp in = ", value(m.fs.ECON.side_1.properties_in[0].temperature))
    print("liquid temp out = ", value(m.fs.ECON.side_1.properties_out[0].temperature))
    print("gas temp in = ", value(m.fs.ECON.side_2.properties_in[0].temperature))
    print("gas temp out = ", value(m.fs.ECON.side_2.properties_out[0].temperature))

    print("\n\n ------------- water wall  ---------")
    print(
        "liquid temp in = ",
        value(m.fs.Water_wall.control_volume.properties_in[0].temperature),
    )
    print(
        "steam temp out = ",
        value(m.fs.Water_wall.control_volume.properties_out[0].temperature),
    )

    print("\n\n ------------- Primary Superheater  ---------")
    print("steam temp in = ", value(m.fs.PrSH.side_1.properties_in[0].temperature))
    print("steam temp out = ", value(m.fs.PrSH.side_1.properties_out[0].temperature))
    print("gas temp in = ", value(m.fs.PrSH.side_2.properties_in[0].temperature))
    print("gas temp out = ", value(m.fs.PrSH.side_2.properties_out[0].temperature))

    print("\n\n ------------- Platen SH  ---------")
    print(
        "steam temp in = ", value(m.fs.PlSH.control_volume.properties_in[0].temperature)
    )
    print(
        "steam temp out = ",
        value(m.fs.PlSH.control_volume.properties_out[0].temperature),
    )

    print("\n\n ------------- Finishing Superheater  ---------")
    print("steam temp in = ", value(m.fs.FSH.side_1.properties_in[0].temperature))
    print(
        "steam temp out (to attmp) = ",
        value(m.fs.FSH.side_1.properties_out[0].temperature),
    )
    print("gas temp in = ", value(m.fs.FSH.side_2.properties_in[0].temperature))
    print("gas temp out = ", value(m.fs.FSH.side_2.properties_out[0].temperature))

    print("\n\n ------------- Attemperator  ---------")
    print(
        "steam temp in = ",
        value(m.fs.ATMP1.Steam.enth_mol[0]),
        value(m.fs.FSH.side_1.properties_out[0].temperature),
    )
    print(
        "steam temp out (to HP turbine) = ",
        value(m.fs.ATMP1.outlet.enth_mol[0]),
        value(m.fs.ATMP1.mixed_state[0].temperature),
    )
    print()

    print("\n\n ------------- Reheater  ---------")
    print("liquid temp in = ", value(m.fs.RH.side_1.properties_in[0].temperature))
    print("liquid temp out = ", value(m.fs.RH.side_1.properties_out[0].temperature))
    print("gas temp in = ", value(m.fs.RH.side_2.properties_in[0].temperature))
    print("gas temp out = ", value(m.fs.RH.side_2.properties_out[0].temperature))


if __name__ == "__main__":
    m, solver = main()

    solver.options = {
        "tol": 1e-6,
        "linear_solver": "ma27",
        "max_iter": 100,
        "halt_on_ampl_error": "yes",
    }
    # initialize each unit at the time
    initialize(m)
    print("flowsheet degrees of freedom = " + str(degrees_of_freedom(m)))
    # unfix inlets to build arcs at the flowsheet level
    unfix_inlets(m)
    # users can change the following fixed values to obtain certain performance
    # m.fs.ATMP1.outlet.enth_mol[0].fix(62710.01)
    # m.fs.Water_wall.heat_duty.unfix()
    # m.fs.RH.side_1_outlet.enth_mol.fix(66043.35)
    # m.fs.PlSH.heat_duty.unfix()
    #    m.fs.ATMP1.SprayWater.flow_mol[0].unfix()
    print("flowsheet degrees of freedom = " + str(degrees_of_freedom(m)))
    results = solver.solve(m, tee=True, symbolic_solver_labels=True)
    print_results(m)

    #    print results in the PFD file (activate if needed)
    #    _stream_dict(m)
    #    df = create_stream_table_dataframe(streams=m._streams, orient="index")
    #    pfd_result("Boiler_Results.svg", m, df)
