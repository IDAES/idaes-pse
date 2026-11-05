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
            BFW -> ECONOMIZER -> Water Wall -> Primary SH -> Platen SH -> Finishing
            Superheate -> HP Turbine -> Reheater -> IP Turbine
        Flue Gas Flow:
            Fire Ball -> Platen SH -> Finishing SH -> Reheater  -> o -> Economizer ->
            Air Preheater -> Primary SH --^

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

    Numerical scaling approach: Scaling is done by initializing the coupled flowsheet
    sequentially with a reconciliation pass, applying each unit’s default scaler, and
    using AutoScaler to fill remaining gaps based on initialized variable magnitudes
    and Jacobian row norms before the final solve.


Created: 1/10/2020 by Boiler subsystem team (M Zamarripa)
Modified: 10/01/2026 (F Hanna)
"""

# pylint: disable=missing-function-docstring

__author__ = "Miguel Zamarripa"

# import dependencies
from collections import OrderedDict
import os
import logging

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.environ import (
    ConcreteModel,
    value,
    TransformationFactory,
    units as pyunits,
)
from pyomo.network import Arc
from pyomo.common.fileutils import this_file_dir
from pyomo.opt import check_optimal_termination

import idaes.logger as idaeslog

# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.util.tags import svg_tag
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state as _set_port
from idaes.core.scaling import AutoScaler, set_scaling_factor
from idaes.core.solvers import get_solver
from idaes.core.scaling.util import (
    get_scaling_factor,
    list_unscaled_constraints,
    list_unscaled_variables,
)

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

# setup logger
_log = idaeslog.getModelLogger(__name__, logging.INFO)


def main():
    """
    Make the flowsheet object, fix some variables, and solve the problem
    """
    # Create a Concrete Model as the top level object
    _log.info("Creating Concrete Model for Boiler Subflowsheet")
    m = ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(dynamic=False)

    m.fs.prop_water = iapws95.Iapws95ParameterBlock()

    # build the boiler subflowsheet
    _log.info("Building Boiler Subflowsheet")
    build_boiler(m.fs)

    # Create a solver
    _log.info("Creating Solver")
    solver = get_solver()
    return (m, solver)


def boiler_hx(fs, phase, radiation):
    """
    Helper function to create a boiler heat exchanger with the given phase and radiation flag.

    Parameters
    ----------
    fs : FlowsheetBlock
        Flowsheet containing the water and flue-gas property packages.
    phase : str
        The phase of the cold side water ('Liq' or 'Vap').
    radiation : bool
        Whether to include radiation in the heat exchanger.
    """
    return BoilerHeatExchanger(
        cold_side={"property_package": fs.prop_water, "has_pressure_change": True},
        hot_side={"property_package": fs.prop_fluegas, "has_pressure_change": True},
        has_holdup=False,
        flow_pattern=HeatExchangerFlowPattern.countercurrent,
        tube_arrangement=TubeArrangement.inLine,
        cold_side_water_phase=phase,
        has_radiation=radiation,
    )


def build_boiler(fs):
    """
    This function builds the boiler subflowsheet within the given flowsheet block.
    This includes creating the necessary unit models and connecting them with arcs to
    represent the flow of water/steam and flue gas through the boiler system.
    The boiler subflowsheet includes two main flow paths:

    1. The water/steam flow path, which consists of the following components:

       * Economizer
       * Water wall
       * Primary superheater
       * Platen superheater
       * Finishing superheater
       * Reheater

    2. The flue gas flow path, which consists of the following components:

       * Finishing superheater
       * Splitter
       * Reheater
       * Mixer

    Parameters
    ----------
    fs : FlowsheetBlock
        The flowsheet block to which the boiler subflowsheet will be added.

    Returns
    -------
    None
    """
    # Add property packages to flowsheet library
    fs.prop_fluegas = FlueGasParameterBlock()

    # Create unit models
    _log.info("Creating Boiler Unit Models")
    # Boiler Economizer
    fs.ECON = boiler_hx(fs, "Liq", False)
    # Boiler Water Wall
    fs.Water_wall = Heater(property_package=fs.prop_water)
    # Primary Superheater
    fs.PrSH = boiler_hx(fs, "Vap", True)
    # Platen Superheater
    fs.PlSH = Heater(property_package=fs.prop_water)
    # Finishing Superheater
    fs.FSH = boiler_hx(fs, "Vap", True)
    # Reheater
    fs.RH = boiler_hx(fs, "Vap", True)
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
    #    fs.bfw2econ = Arc (source=fs.FWH8.outlet,
    #                       destination=fs.ECON.cold_side_inlet)
    _log.info("Creating Boiler Arcs")
    fs.econ2ww = Arc(source=fs.ECON.cold_side_outlet, destination=fs.Water_wall.inlet)
    fs.ww2prsh = Arc(source=fs.Water_wall.outlet, destination=fs.PrSH.cold_side_inlet)
    fs.prsh2plsh = Arc(source=fs.PrSH.cold_side_outlet, destination=fs.PlSH.inlet)
    fs.plsh2fsh = Arc(source=fs.PlSH.outlet, destination=fs.FSH.cold_side_inlet)
    fs.FSHtoATMP1 = Arc(source=fs.FSH.cold_side_outlet, destination=fs.ATMP1.Steam)
    # The attemperator outlet is out of scope for the boiler subflowsheet
    # it will be connected to the HP turbine inlet in the full plant model
    #    fs.fsh2hpturbine = Arc(source=fs.ATMP1.outlet,
    #                           destination=fs.HPTinlet)

    # Flue gas route
    # water wall connected with boiler block (to fix the heat duty)
    # platen SH connected with boiler block (to fix the heat duty)
    # Finishing superheater connected with a flowsheet level constraint
    fs.fg_fsh2_separator = Arc(source=fs.FSH.hot_side_outlet, destination=fs.Spl1.inlet)
    fs.fg_fsh2rh = Arc(source=fs.Spl1.outlet_1, destination=fs.RH.hot_side_inlet)
    fs.fg_fsh2PrSH = Arc(source=fs.Spl1.outlet_2, destination=fs.PrSH.hot_side_inlet)
    fs.fg_rhtomix = Arc(source=fs.RH.hot_side_outlet, destination=fs.mix1.Reheat_out)
    fs.fg_prsh2mix = Arc(source=fs.PrSH.hot_side_outlet, destination=fs.mix1.PrSH_out)
    fs.fg_mix2econ = Arc(source=fs.mix1.outlet, destination=fs.ECON.hot_side_inlet)

    TransformationFactory("network.expand_arcs").apply_to(fs)


def scale_unit(unit, label):
    """
    Apply the unit model's default scaler, then use AutoScaler to fill any
    remaining variable and constraint scaling gaps based on initialized
    variable magnitudes and Jacobian row norms.

    Parameters
    ----------
    unit : ProcessBlockData
        Initialized IDAES unit model to scale.
    label : str
        Descriptive unit name used in the scaling summary output.

    Returns
    -------
    tuple
        A two-element tuple containing:

        default_scaler
            The unit model's default scaler instance, or ``None`` if the unit
            does not provide one.
        gap_scaler : AutoScaler
            The AutoScaler instance used to fill missing scaling factors.
    """
    # The Helmholtz/IAPWS default scaler requires representative molar-flow
    # factors. Supply these from initialized nominal values before invoking it.
    required_flow_factors = 0
    for variable in unit.component_data_objects(pyo.Var, descend_into=True):
        if variable.parent_component().local_name == "flow_mol":
            nominal = abs(pyo.value(variable))
            if nominal > 1e-12 and get_scaling_factor(variable) is None:
                set_scaling_factor(variable, 1 / nominal, overwrite=False)
                required_flow_factors += 1

    default_scaler = unit.default_scaler(overwrite=False)
    if default_scaler is not None:
        default_scaler.scale_model(unit)

    variables_after_default = list_unscaled_variables(
        unit, descend_into=True, include_fixed=False
    )
    constraints_after_default = list_unscaled_constraints(unit, descend_into=True)

    # Fill only gaps left by the unit's default scaler. Existing factors are
    # protected because overwrite=False.
    gap_scaler = AutoScaler(overwrite=False)
    if variables_after_default:
        gap_scaler.scale_variables_by_magnitude(unit, descend_into=True)
    if constraints_after_default:
        gap_scaler.scale_constraints_by_jacobian_norm(unit, norm=2, descend_into=True)

    remaining_variables = list_unscaled_variables(
        unit, descend_into=True, include_fixed=False
    )
    remaining_constraints = list_unscaled_constraints(unit, descend_into=True)
    scaler_name = (
        type(default_scaler).__name__ if default_scaler is not None else "None"
    )
    print(
        f"{label}: default={scaler_name}, required flow factors={required_flow_factors}, "
        f"gaps before AutoScaler=({len(variables_after_default)} vars, "
        f"{len(constraints_after_default)} cons), remaining="
        f"({len(remaining_variables)} vars, {len(remaining_constraints)} cons)"
    )

    return default_scaler, gap_scaler


def initialize(m):
    """
    This function initializes the boiler subflowsheet model by setting initial values for

    the economizer, water wall, primary superheater, platen superheater, finishing superheater,

    attemperator, splitter, reheater, and mixer units.
    It also propagates the state of the streams between the units and scales the units using the
    scale_unit helper function.


    Parameters
    ----------
    m : ConcreteModel
        The Pyomo model containing the flowsheet to be initialized.

    Returns
    -------
    None
    """
    # ------------- ECONOMIZER -----------------------------------------------
    # BFW Boiler Feed Water inlet temeperature = 555 F = 563.706 K
    _log.info("Setting Economizer Inlet Conditions")
    m.fs.ECON.cold_side_inlet.flow_mol[0].fix(24194.177)  # mol/s
    m.fs.ECON.cold_side_inlet.enth_mol[0].fix(14990.97)
    m.fs.ECON.cold_side_inlet.pressure[0].fix(26922222.222)  # Pa

    # FLUE GAS Inlet from Primary Superheater
    FGrate = 21290.6999  # mol/s
    comp = {
        "H2O": 0.0869,
        "CO2": 0.1449,
        "N2": 0.7434,
        "O2": 0.0247,
        "NO": 0.0006,
        "SO2": 0.002,
    }
    for j, x in comp.items():
        m.fs.ECON.hot_side_inlet.flow_mol_comp[0, j].fix(FGrate * x)
    m.fs.ECON.hot_side_inlet.temperature[0].fix(682.335)  # K
    m.fs.ECON.hot_side_inlet.pressure[0].fix(100145)  # Pa

    # Economizer design variables and parameters
    ITM = 0.0254  # inch to meter conversion
    # Based on NETL Baseline Report Rev3
    # calc inner diameter (2 = outer diameter, thickness = 0.188)
    m.fs.ECON.tube_di.fix((2 - 2 * 0.188) * ITM)
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

    # initialize the economizer unit
    _log.info("Initializing Economizer Unit")
    m.fs.ECON.initialize(outlvl=logging.INFO)
    _log.info("Scaling Economizer Unit")
    scale_unit(m.fs.ECON, "ECON")

    # ------- Water wall Superheater ----------------------------------------
    # propagate the economizer outlet state to the water-wall inlet
    _log.info("Setting Water Wall Inlet Conditions")
    _set_port(arc=m.fs.econ2ww)
    m.fs.Water_wall.heat_duty[:].fix(7.51e8)  # 8.76e8

    # Initialize the water wall unit
    _log.info("Initializing Water Wall Unit")
    m.fs.Water_wall.initialize(outlvl=logging.INFO)
    _log.info("Scaling Water Wall Unit")
    scale_unit(m.fs.Water_wall, "Water_wall")

    # --------- Primary Superheater -----------------------------------------
    # Steam from water wall
    # propagate the water-wall outlet state to the PrSH cold-side inlet
    _log.info("Setting Primary Superheater Inlet Conditions")
    _set_port(arc=m.fs.ww2prsh)

    # FLUE GAS Inlet from Primary Superheater
    FGrate = 21290.6999 * 0.18  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    for j, x in comp.items():
        m.fs.PrSH.hot_side_inlet.flow_mol_comp[0, j].fix(FGrate * x)
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

    # Initialize the PrSH unit
    _log.info("Initializing Primary Superheater Unit")
    m.fs.PrSH.initialize(outlvl=logging.INFO)
    # scale the PrSH unit
    _log.info("Scaling Primary Superheater Unit")
    scale_unit(m.fs.PrSH, "PrSH")

    # --------- Platen Superheater ------------------------------------------
    # propagate the PrSH outlet state to the PlSH inlet
    _log.info("Setting Platen Superheater Inlet Conditions")
    _set_port(arc=m.fs.prsh2plsh)
    # Fix the PlSH heat duty
    m.fs.PlSH.heat_duty[:].fix(5.5e7)
    # initialize the PlSH unit
    _log.info("Initializing Platen Superheater Unit")
    m.fs.PlSH.initialize(outlvl=logging.INFO)
    # scale the PlSH unit
    _log.info("Scaling Platen Superheater Unit")
    scale_unit(m.fs.PlSH, "PlSH")

    #  -------- Finishing Superheater ----------------------------------------
    # propagate the PlSH outlet state to the FSH cold-side inlet
    _log.info("Setting Finishing Superheater Inlet Conditions")
    _set_port(arc=m.fs.plsh2fsh)

    # FLUE GAS Inlet from Primary Superheater
    FGrate = 21290.6999  # mol/s
    # Use FG molar composition to set component flow rates (baseline report)
    for j, x in comp.items():
        m.fs.FSH.hot_side_inlet.flow_mol_comp[0, j].fix(FGrate * x)
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

    # Initialize the FSH unit
    _log.info("Initializing Finishing Superheater Unit")
    m.fs.FSH.initialize(outlvl=logging.INFO)
    # scale the FSH unit
    _log.info("Scaling Finishing Superheater Unit")
    scale_unit(m.fs.FSH, "FSH")

    # --------- Attemperator inputs ------------------------------------------
    # propagate the FSH steam outlet state to the ATMP1 steam inlet
    _log.info("Setting Attemperator Inlet Conditions")
    _set_port(arc=m.fs.FSHtoATMP1)
    # Fixed SprayWater from FW pump splitter (splitter is needded)
    hatt2 = value(iapws95.htpx(563.15 * pyunits.K, 2.5449e7 * pyunits.Pa))
    m.fs.ATMP1.SprayWater_state[:].flow_mol.fix(0.001)
    m.fs.ATMP1.SprayWater_state[:].pressure.fix(1.22e8)
    m.fs.ATMP1.SprayWater_state[:].enth_mol.fix(hatt2)
    # Initialize the ATMP1 unit
    _log.info("Initializing Attemperator Unit")
    m.fs.ATMP1.initialize(outlvl=logging.INFO)
    # scale the ATMP1 unit
    _log.info("Scaling Attemperator Unit")
    scale_unit(m.fs.ATMP1, "ATMP1")

    # --------- Splitter ----------------------------------------------------
    # splitter flue gas from Finishing SH to Reheater and Primary SH
    _log.info("Setting Splitter Fraction to 0.75 for Reheater and 0.25 for Primary SH")
    m.fs.Spl1.split_fraction[0, "outlet_1"].fix(0.75)  # 0.85)
    # FLUE GAS Inlet from Primary Superheater
    # Propagate the FSH outlet state to the splitter inlet.
    _log.info("Setting Splitter Inlet Conditions")
    _set_port(arc=m.fs.fg_fsh2_separator)
    for variable in m.fs.Spl1.inlet.vars.values():
        variable.fix()
    # Initialize the splitter unit
    _log.info("Initializing Splitter Unit")
    m.fs.Spl1.initialize(outlvl=logging.INFO)
    # scale the splitter unit
    _log.info("Scaling Splitter Unit")
    m.fs.Spl1.scaler = scale_unit(m.fs.Spl1, "Spl1")

    # ----------Propagate Splitter outlet -----------------------------------
    # Propagate the splitter results before initializing either gas branch.
    # Propagate the splitter outlet state 1 to the RH hot-side inlet
    _set_port(arc=m.fs.fg_fsh2rh)
    # Propagate the splitter outlet state 2 to the PrSH hot-side inlet
    _set_port(arc=m.fs.fg_fsh2PrSH, overwrite_fixed=True)

    # ----------- Reheater Superheater --------------------------------------
    #   Steam from HP Turbine outlet
    _log.info("Setting Reheater Inlet Conditions")
    m.fs.RH.cold_side_inlet.flow_mol[0].fix(21235.27)  # mol/s
    m.fs.RH.cold_side_inlet.enth_mol[0].fix(53942.7569)  # J/mol
    m.fs.RH.cold_side_inlet.pressure[0].fix(3677172.33638)  # Pascals

    # Hold the propagated splitter outlet while the RH is initialized alone.
    for variable in m.fs.RH.hot_side_inlet.vars.values():
        variable.fix()

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

    # Initialize the RH unit
    _log.info("Initializing Reheater Unit")
    m.fs.RH.initialize(outlvl=logging.INFO)
    # scale the RH unit
    _log.info("Scaling Reheater Unit")
    scale_unit(m.fs.RH, "RH")

    # Reinitialize PrSH with the propagated splitter branch. Its steam-side inlet
    # remains the initialized Water_wall outlet from the earlier steam-path pass.
    u = m.fs.PrSH
    _log.info("Reinitializing Primary Superheater Unit with Splitter Outlet")
    u.initialize(outlvl=logging.INFO)

    # Refresh the downstream steam path because the updated PrSH gas inlet can
    # change its cold-side outlet and, consequently, the FSH outlet states.
    _set_port(arc=m.fs.prsh2plsh)
    m.fs.PlSH.initialize(outlvl=logging.INFO)
    _set_port(arc=m.fs.plsh2fsh)
    m.fs.FSH.initialize(outlvl=logging.INFO)

    # Refresh the splitter and both gas branches once after the steam-path update.
    _set_port(arc=m.fs.fg_fsh2_separator, overwrite_fixed=True)
    m.fs.Spl1.initialize(outlvl=logging.INFO)
    _set_port(arc=m.fs.fg_fsh2rh, overwrite_fixed=True)
    m.fs.RH.initialize(outlvl=logging.INFO)
    _set_port(arc=m.fs.fg_fsh2PrSH, overwrite_fixed=True)
    m.fs.PrSH.initialize(outlvl=logging.INFO)

    # Initialize the gas mixer from the refreshed RH and PrSH outlets.
    _log.info(
        "Setting Mixer Inlet Conditions from Reheater and Primary Superheater Outlets"
    )
    for inlet, outlet in (
        (m.fs.mix1.Reheat_out, m.fs.RH.hot_side_outlet),
        (m.fs.mix1.PrSH_out, m.fs.PrSH.hot_side_outlet),
    ):
        for j in comp:
            inlet.flow_mol_comp[0, j].fix(value(outlet.flow_mol_comp[0, j]))
        inlet.temperature[0].fix(value(outlet.temperature[0]))
        inlet.pressure[0].fix(value(outlet.pressure[0]))
    # Initialize the mixer unit
    _log.info("Initializing Mixer Unit")
    m.fs.mix1.initialize(outlvl=logging.INFO)
    # scale the mixer unit
    _log.info("Scaling Mixer Unit")
    scale_unit(m.fs.mix1, "mix1")

    # Release states determined by connected flue-gas arcs. FSH hot-side inlet
    # remains the external boiler-gas boundary; RH cold-side inlet remains the
    # external HP-turbine-exhaust boundary.
    connected_inlet_ports = (
        m.fs.Spl1.inlet,
        m.fs.RH.hot_side_inlet,
        m.fs.PrSH.hot_side_inlet,
        m.fs.mix1.Reheat_out,
        m.fs.mix1.PrSH_out,
        m.fs.ECON.hot_side_inlet,
    )
    for port in connected_inlet_ports:
        for variable in port.vars.values():
            variable.unfix()

    # Populate consistent starting values along the gas path.
    _set_port(arc=m.fs.fg_fsh2_separator)
    _set_port(arc=m.fs.fg_fsh2rh)
    _set_port(arc=m.fs.fg_fsh2PrSH)
    _set_port(arc=m.fs.fg_rhtomix)
    _set_port(arc=m.fs.fg_prsh2mix)
    _set_port(arc=m.fs.fg_mix2econ)

    # ------------------------------------------------------------------------
    print("initialization done")


def scale_solve(m):
    """
    This function scales the entire flowsheet and solves it.

    Function logic:
    1. Check the degrees of freedom of the flowsheet. If it is not zero, raise an error.
    2. Create an AutoScaler instance to scale the entire flowsheet without overwriting existing factors.
    3. Scale variables by magnitude and constraints by Jacobian norm, without overwriting existing factors.
    4. Solve the flowsheet using the specified solver and print the results.

    Parameters
    ----------
    m : ConcreteModel
        The Pyomo model containing the flowsheet to be scaled and solved.

    Returns
    -------
    results : SolverResults
        The results of the solver after solving the flowsheet.
    """
    _log.info("Scaling and Solving Full Flowsheet")
    # make sure the flowsheet has 0 degrees of freedom before solving
    if degrees_of_freedom(m.fs) != 0:
        raise RuntimeError("Connected flowsheet is not square")
    # scale the entire flowsheet, but do not overwrite existing factors
    fs_scaler = AutoScaler(overwrite=False)
    # scale variables by magnitude, but do not overwrite existing factors
    fs_scaler.scale_variables_by_magnitude(m.fs, descend_into=True)
    # scale constraints by Jacobian norm, but do not overwrite existing factors
    fs_scaler.scale_constraints_by_jacobian_norm(m.fs, norm=2, descend_into=True)
    _log.info(
        "Unscaled active variables: {}".format(
            len(list_unscaled_variables(m.fs, descend_into=True, include_fixed=False))
        )
    )
    _log.info(
        "Unscaled constraints: {}".format(
            len(list_unscaled_constraints(m.fs, descend_into=True))
        )
    )

    # solve the connected flowsheet with the updated flue-gas path and refreshed steam path.
    solver = get_solver()
    results = solver.solve(m.fs, tee=True)
    _log.info(
        "Solver Status: {}, Termination Condition: {}".format(
            results.solver.status, results.solver.termination_condition
        )
    )
    if not check_optimal_termination(results):
        raise RuntimeError("Connected solve failed")

    return results


def pfd_result(outfile, m, df):
    tags = {}
    for i in df.index:
        tags[i + "_F"] = df.loc[i, "Molar Flow"]
        tags[i + "_T"] = df.loc[i, "T"]
        tags[i + "_P"] = df.loc[i, "P"]
        tags[i + "_X"] = df.loc[i, "Vapor Fraction"]

    tags["FG_2_RH_Fm"] = value(m.fs.RH.hot_side.properties_in[0].flow_mass)
    tags["FG_2_RH_T"] = value(m.fs.RH.hot_side.properties_in[0].temperature)
    tags["FG_2_RH_P"] = value(m.fs.RH.hot_side.properties_in[0].pressure)

    tags["FG_RH_2_Mix_Fm"] = value(m.fs.RH.hot_side.properties_out[0].flow_mass)
    tags["FG_RH_2_Mix_T"] = value(m.fs.RH.hot_side.properties_out[0].temperature)
    tags["FG_RH_2_Mix_P"] = value(m.fs.RH.hot_side.properties_out[0].pressure)

    tags["FG_2_FSH_Fm"] = value(m.fs.FSH.hot_side.properties_in[0].flow_mass)
    tags["FG_2_FSH_T"] = value(m.fs.FSH.hot_side.properties_in[0].temperature)
    tags["FG_2_FSH_P"] = value(m.fs.FSH.hot_side.properties_in[0].pressure)

    tags["FG_2_PrSH_Fm"] = value(m.fs.PrSH.hot_side.properties_in[0].flow_mass)
    tags["FG_2_PrSH_T"] = value(m.fs.PrSH.hot_side.properties_in[0].temperature)
    tags["FG_2_PrSH_P"] = value(m.fs.PrSH.hot_side.properties_in[0].pressure)

    tags["FG_PrSH_2_Mix_Fm"] = value(m.fs.PrSH.hot_side.properties_out[0].flow_mass)
    tags["FG_PrSH_2_Mix_T"] = value(m.fs.PrSH.hot_side.properties_out[0].temperature)
    tags["FG_PrSH_2_Mix_P"] = value(m.fs.PrSH.hot_side.properties_out[0].pressure)

    tags["FG_2_ECON_Fm"] = value(m.fs.ECON.hot_side.properties_in[0].flow_mass)
    tags["FG_2_ECON_T"] = value(m.fs.ECON.hot_side.properties_in[0].temperature)
    tags["FG_2_ECON_P"] = value(m.fs.ECON.hot_side.properties_in[0].pressure)

    tags["FG_2_AIRPH_Fm"] = value(m.fs.ECON.hot_side.properties_out[0].flow_mass)
    tags["FG_2_AIRPH_T"] = value(m.fs.ECON.hot_side.properties_out[0].temperature)
    tags["FG_2_AIRPH_P"] = value(m.fs.ECON.hot_side.properties_out[0].pressure)

    tags["FG_2_STACK_Fm"] = value(m.fs.ECON.hot_side.properties_out[0].flow_mass)
    tags["FG_2_STACK_T"] = value(m.fs.ECON.hot_side.properties_out[0].temperature)
    tags["FG_2_STACK_P"] = value(m.fs.ECON.hot_side.properties_out[0].pressure)

    original_svg_file = os.path.join(this_file_dir(), "Boiler_scpc_PFD.svg")
    with open(original_svg_file, "r") as f:
        svg_tag(tags, f, outfile=outfile)


def _stream_dict(m):
    """Adds _streams to m, which contains a dictionary of streams for display

    Parameters
    ----------
        m (ConcreteModel): A Pyomo model from create_model()

    Returns
    -------
        None
    """
    # We control m
    # pylint: disable-next=protected-access
    m._streams = OrderedDict(
        [
            ("MS", m.fs.ATMP1.mixed_state),
            ("ATMP_In", m.fs.FSH.cold_side.properties_out),
            ("FSH_In", m.fs.FSH.cold_side.properties_in),
            ("PrSH_IN", m.fs.PrSH.cold_side.properties_in),
            ("RHT_COLD", m.fs.RH.cold_side.properties_in),
            ("RHT_HOT", m.fs.RH.cold_side.properties_out),
            ("PlatenSH_IN", m.fs.PlSH.control_volume.properties_in),
            ("BFW", m.fs.ECON.cold_side.properties_in),
            ("ECON_OUT", m.fs.ECON.cold_side.properties_out),
        ]
    )


def print_results(m):
    """
    This function prints the results of the boiler subflowsheet model.

    Parameters
    ----------
    m : ConcreteModel
        A Pyomo model from create_model()

    Returns
    -------
    None
    """
    _log.info("Printing Results")
    _log.info("Results")
    _log.info("")

    _log.info(
        "viscosity gas side = %s", m.fs.PrSH.hot_side.properties_in[0].visc_d.value
    )
    _log.info(
        "conductivity gas side = %s",
        m.fs.PrSH.hot_side.properties_in[0].therm_cond.value,
    )
    _log.info("velocity_tube = %s", m.fs.PrSH.v_tube[0].value)
    _log.info("velocity_shell = %s", m.fs.PrSH.v_shell[0].value)
    _log.info("Re_tube = %s", m.fs.PrSH.N_Re_tube[0].value)
    _log.info("Re_shell = %s", m.fs.PrSH.N_Re_shell[0].value)
    _log.info("hconv_tube = %s", m.fs.PrSH.hconv_tube[0].value)
    _log.info("hconv_shell_rad = %s", m.fs.PrSH.hconv_shell_rad[0].value)
    _log.info("hconv_shell_conv = %s", m.fs.PrSH.hconv_shell_conv[0].value)
    _log.info("hconv_shell_total = %s", m.fs.PrSH.hconv_shell_total[0].value)
    _log.info("driving force = %s", value(m.fs.PrSH.delta_temperature[0]))
    _log.info("dT_inlet = %s", m.fs.PrSH.deltaT_1[0].value)
    _log.info("dT_outlet = %s", m.fs.PrSH.deltaT_2[0].value)
    _log.info("deltaP tube = %s", m.fs.PrSH.deltaP_tube[0].value)
    _log.info("deltaP shell = %s", m.fs.PrSH.deltaP_shell[0].value)
    _log.info("mbl = %s", value(m.fs.PrSH.mbl))

    if m.fs.PrSH.config.has_radiation is True:
        _log.info("gas emissivity = %s", m.fs.PrSH.gas_emissivity[0].value)
        _log.info("gas emissivity div2 = %s", m.fs.PrSH.gas_emissivity_div2[0].value)
        _log.info("gas emissivity mul2 = %s", m.fs.PrSH.gas_emissivity_mul2[0].value)
        _log.info("gas gray fraction = %s", m.fs.PrSH.gas_gray_fraction[0].value)
    _log.info(
        "liquid density in = %s",
        value(m.fs.PrSH.cold_side.properties_in[0].dens_mass_phase["Liq"]),
    )

    _log.info(
        "liquid density out = %s",
        value(m.fs.PrSH.cold_side.properties_out[0].dens_mass_phase["Liq"]),
    )
    _log.info("heat transfer area = %s", value(m.fs.PrSH.area))
    _log.info(
        "overall heat transfer = %s",
        value(m.fs.PrSH.overall_heat_transfer_coefficient[0]),
    )

    _log.info("\n\n ------------- Economizer   ---------")
    _log.info(
        "liquid temp in = %s", value(m.fs.ECON.cold_side.properties_in[0].temperature)
    )
    _log.info(
        "liquid temp out = %s", value(m.fs.ECON.cold_side.properties_out[0].temperature)
    )
    _log.info(
        "gas temp in = %s", value(m.fs.ECON.hot_side.properties_in[0].temperature)
    )
    _log.info(
        "gas temp out = %s", value(m.fs.ECON.hot_side.properties_out[0].temperature)
    )

    _log.info("\n\n ------------- water wall  ---------")
    _log.info(
        "liquid temp in = %s",
        value(m.fs.Water_wall.control_volume.properties_in[0].temperature),
    )
    _log.info(
        "steam temp out = %s",
        value(m.fs.Water_wall.control_volume.properties_out[0].temperature),
    )

    _log.info("\n\n ------------- Primary Superheater  ---------")
    _log.info(
        "steam temp in = %s", value(m.fs.PrSH.cold_side.properties_in[0].temperature)
    )
    _log.info(
        "steam temp out = %s", value(m.fs.PrSH.cold_side.properties_out[0].temperature)
    )
    _log.info(
        "gas temp in = %s", value(m.fs.PrSH.hot_side.properties_in[0].temperature)
    )
    _log.info(
        "gas temp out = %s", value(m.fs.PrSH.hot_side.properties_out[0].temperature)
    )

    _log.info("\n\n ------------- Platen SH  ---------")
    _log.info(
        "steam temp in = %s",
        value(m.fs.PlSH.control_volume.properties_in[0].temperature),
    )
    _log.info(
        "steam temp out = %s",
        value(m.fs.PlSH.control_volume.properties_out[0].temperature),
    )

    _log.info("\n\n ------------- Finishing Superheater  ---------")
    _log.info(
        "steam temp in = %s", value(m.fs.FSH.cold_side.properties_in[0].temperature)
    )
    _log.info(
        "steam temp out (to attmp) = %s",
        value(m.fs.FSH.cold_side.properties_out[0].temperature),
    )
    _log.info("gas temp in = %s", value(m.fs.FSH.hot_side.properties_in[0].temperature))
    _log.info(
        "gas temp out = %s", value(m.fs.FSH.hot_side.properties_out[0].temperature)
    )

    _log.info("\n\n ------------- Attemperator  ---------")
    _log.info(
        "steam temp in = %s, temperature = %s",
        value(m.fs.ATMP1.Steam.enth_mol[0]),
        value(m.fs.FSH.cold_side.properties_out[0].temperature),
    )
    _log.info(
        "steam temp out (to HP turbine) = %s, temperature = %s",
        value(m.fs.ATMP1.outlet.enth_mol[0]),
        value(m.fs.ATMP1.mixed_state[0].temperature),
    )
    _log.info("")

    _log.info("\n\n ------------- Reheater  ---------")
    _log.info(
        "liquid temp in = %s", value(m.fs.RH.cold_side.properties_in[0].temperature)
    )
    _log.info(
        "liquid temp out = %s", value(m.fs.RH.cold_side.properties_out[0].temperature)
    )
    _log.info("gas temp in = %s", value(m.fs.RH.hot_side.properties_in[0].temperature))
    _log.info(
        "gas temp out = %s", value(m.fs.RH.hot_side.properties_out[0].temperature)
    )


def unfix_inlets(m):
    """
    This method unfixes the inlet conditions for the economizer, primary superheater, water
    wall, reheater, finishing superheater, platen superheater, and splitter units in the boiler
    subflowsheet model.

    This function is intended to be used in a complete flowsheet where the conditions for the reheater,
    economizer and attemperator are determined by the rest steam cycle flowsheet unit models.

    The remaining parameters are unfixed to ensure that the full model can be solved with the
    appropriate degrees of freedom.

    Parameters
    ----------
    m : ConcreteModel
        A Pyomo model from create_model()

    Returns
    -------
    None
    """
    # Use FG molar composition to set component flow rates (baseline report)
    _log.info("Unfixing hot-side inlet conditions for economizer")
    m.fs.ECON.hot_side_inlet.flow_mol_comp[0, "H2O"].unfix()
    m.fs.ECON.hot_side_inlet.flow_mol_comp[0, "CO2"].unfix()
    m.fs.ECON.hot_side_inlet.flow_mol_comp[0, "N2"].unfix()
    m.fs.ECON.hot_side_inlet.flow_mol_comp[0, "O2"].unfix()
    m.fs.ECON.hot_side_inlet.flow_mol_comp[0, "NO"].unfix()
    m.fs.ECON.hot_side_inlet.flow_mol_comp[0, "SO2"].unfix()
    m.fs.ECON.hot_side_inlet.temperature[0].unfix()
    m.fs.ECON.hot_side_inlet.pressure[0].unfix()

    # PrSH Primary superheater inlets (steam and flue gas) --------------------
    _log.info("Unfixing cold-side inlet conditions for primary superheater")
    m.fs.PrSH.cold_side_inlet.flow_mol.unfix()
    m.fs.PrSH.cold_side_inlet.enth_mol[0].unfix()
    m.fs.PrSH.cold_side_inlet.pressure[0].unfix()

    # Use FG molar composition to set component flow rates (baseline report)
    _log.info("Unfixing hot-side inlet conditions for primary superheater")
    m.fs.PrSH.hot_side_inlet.flow_mol_comp[0, "H2O"].unfix()
    m.fs.PrSH.hot_side_inlet.flow_mol_comp[0, "CO2"].unfix()
    m.fs.PrSH.hot_side_inlet.flow_mol_comp[0, "N2"].unfix()
    m.fs.PrSH.hot_side_inlet.flow_mol_comp[0, "O2"].unfix()
    m.fs.PrSH.hot_side_inlet.flow_mol_comp[0, "NO"].unfix()
    m.fs.PrSH.hot_side_inlet.flow_mol_comp[0, "SO2"].unfix()
    m.fs.PrSH.hot_side_inlet.temperature[0].unfix()
    m.fs.PrSH.hot_side_inlet.pressure[0].unfix()

    # WaterWall  water from economizer  ---------------------------------------
    _log.info("Unfixing inlet conditions for water wall")
    m.fs.Water_wall.inlet[:].flow_mol.unfix()
    m.fs.Water_wall.inlet[:].enth_mol.unfix()
    m.fs.Water_wall.inlet[:].pressure.unfix()

    # Use FG molar composition to set component flow rates (baseline report)
    _log.info("Unfixing hot-side inlet conditions for reheater")
    m.fs.RH.hot_side_inlet.flow_mol_comp[0, "H2O"].unfix()
    m.fs.RH.hot_side_inlet.flow_mol_comp[0, "CO2"].unfix()
    m.fs.RH.hot_side_inlet.flow_mol_comp[0, "N2"].unfix()
    m.fs.RH.hot_side_inlet.flow_mol_comp[0, "O2"].unfix()
    m.fs.RH.hot_side_inlet.flow_mol_comp[0, "NO"].unfix()
    m.fs.RH.hot_side_inlet.flow_mol_comp[0, "SO2"].unfix()
    m.fs.RH.hot_side_inlet.temperature[0].unfix()
    m.fs.RH.hot_side_inlet.pressure[0].unfix()

    # Finishing Superheater (steam from Platen SH)-----------------------------
    _log.info("Unfixing cold-side inlet conditions for finishing superheater")
    m.fs.FSH.cold_side_inlet.flow_mol.unfix()
    m.fs.FSH.cold_side_inlet.enth_mol[0].unfix()
    m.fs.FSH.cold_side_inlet.pressure[0].unfix()

    # Platen SH steam inlet conditions from Water Wall-------------------------
    _log.info("Unfixing inlet conditions for platen superheater")
    m.fs.PlSH.inlet[:].flow_mol.unfix()
    m.fs.PlSH.inlet[:].enth_mol.unfix()
    m.fs.PlSH.inlet[:].pressure.unfix()

    # unfix Splitter flue gas inlet to RH and PrSH
    _log.info("Unfixing inlet conditions for splitter")
    m.fs.Spl1.inlet.flow_mol_comp[0, "H2O"].unfix()
    m.fs.Spl1.inlet.flow_mol_comp[0, "CO2"].unfix()
    m.fs.Spl1.inlet.flow_mol_comp[0, "N2"].unfix()
    m.fs.Spl1.inlet.flow_mol_comp[0, "O2"].unfix()
    m.fs.Spl1.inlet.flow_mol_comp[0, "NO"].unfix()
    m.fs.Spl1.inlet.flow_mol_comp[0, "SO2"].unfix()
    m.fs.Spl1.inlet.temperature[0].unfix()
    m.fs.Spl1.inlet.pressure[0].unfix()

    # unfix mixer inlets (now connected)
    _log.info("Unfixing outlet conditions for mixer")
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
    _log.info("Unfixing steam inlet conditions for attemperator")
    m.fs.ATMP1.Steam.flow_mol.unfix()
    m.fs.ATMP1.Steam.enth_mol.unfix()
    m.fs.ATMP1.Steam.pressure.unfix()


if __name__ == "__main__":
    # generate the model and solver
    m, solver = main()
    # set solver options
    solver.options = {
        "tol": 1e-6,
        "linear_solver": "ma27",
        "max_iter": 100,
        "halt_on_ampl_error": "yes",
    }
    # initialize the model
    initialize(m)
    # Scale and solve the model
    results = scale_solve(m)
    # print the results
    print_results(m)
