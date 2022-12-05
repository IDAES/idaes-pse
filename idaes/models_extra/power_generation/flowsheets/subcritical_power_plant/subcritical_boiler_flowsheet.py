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
Dynamic sub-flowsheet for a subcritical 300MWe boiler system
"""

# Import Pyomo libraries
import pyomo.environ as pyo
from pyomo.network import Arc

# Import IDAES core
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.initialization import propagate_state as _set_port
from idaes.core.solvers import get_solver
from idaes.core import FlowsheetBlock
import idaes.logger as idaeslog

# Import IDAES standard unit model
from idaes.models.unit_models import Mixer
from idaes.models.properties import iapws95
from idaes.models_extra.power_generation.properties import FlueGasParameterBlock

# Import IDAES power generation unit models
from idaes.models_extra.power_generation.unit_models.helm import (
    HelmMixer,
    MomentumMixingType,
    HelmSplitter,
)

from idaes.models_extra.power_generation.unit_models import (
    Drum1D,
    HeatExchangerWith3Streams,
    WaterPipe,
    HeatExchangerCrossFlow2D_Header,
    Downcomer,
    WaterwallSection,
    BoilerFireside,
    SteamHeater,
)

# Import boiler fire-side surrogate model as dictionary
from idaes.models_extra.power_generation.flowsheets.subcritical_power_plant.generic_surrogate_dict import (
    data_dic,
)
import idaes.core.util.scaling as iscale
from idaes.core.util.dyn_utils import copy_values_at_time, copy_non_time_indexed_values

import matplotlib.pyplot as plt

__author__ = "Boiler Subsystem Team (J. Ma, M. Zamarripa)"


def add_unit_models(m):
    """
    Function to add unit operation models to the flowsheet
    """
    # Models will be added to the boiler sub-flowsheet of the main flowsheet
    fs = m.fs_main.fs_blr

    # Water and ideal gas properties defined on main flowsheet
    prop_water = m.fs_main.prop_water
    prop_gas = m.fs_main.prop_gas

    # Define 12 vertical waterwall zones on the fire side of the boiler
    fs.ww_zones = pyo.RangeSet(12)

    # Unit model for boiler fire side based on surrogate
    fs.aBoiler = BoilerFireside(
        dynamic=False,
        property_package=prop_gas,
        calculate_PA_SA_flows=True,
        number_of_zones=12,
        has_platen_superheater=True,
        has_roof_superheater=True,
        surrogate_dictionary=data_dic,
    )

    # Unit model for boiler drum
    fs.aDrum = Drum1D(
        property_package=prop_water,
        has_holdup=True,
        has_heat_transfer=True,
        has_pressure_change=True,
        finite_elements=4,
        drum_inner_diameter=1.8,
        drum_thickness=0.13,
    )

    # Unit model for splitter from drum to downcomers and blowdown
    fs.blowdown_split = HelmSplitter(
        dynamic=False,
        property_package=prop_water,
        outlet_list=["FW_Downcomer", "FW_Blowdown"],
    )

    # Unit model for downcomer
    fs.aDowncomer = Downcomer(
        dynamic=False,
        property_package=prop_water,
        has_holdup=True,
        has_heat_transfer=True,
    )

    # Unit models for 12 waterwall sections
    fs.Waterwalls = WaterwallSection(
        fs.ww_zones,
        has_holdup=True,
        property_package=prop_water,
        has_heat_transfer=True,
        has_pressure_change=True,
    )

    # Unit model for roof superheater
    fs.aRoof = SteamHeater(
        dynamic=False,
        property_package=prop_water,
        has_holdup=True,
        has_heat_transfer=True,
        has_pressure_change=True,
        single_side_only=True,
    )

    # Unit model for platen superheater
    fs.aPlaten = SteamHeater(
        dynamic=False,
        property_package=prop_water,
        has_holdup=True,
        has_heat_transfer=True,
        has_pressure_change=True,
        single_side_only=False,
    )

    # Unit model for 1st reheater
    fs.aRH1 = HeatExchangerCrossFlow2D_Header(
        tube_side={"property_package": prop_water, "has_pressure_change": True},
        shell_side={"property_package": prop_gas, "has_pressure_change": True},
        finite_elements=4,
        flow_type="counter_current",
        tube_arrangement="in-line",
        tube_side_water_phase="Vap",
        has_radiation=True,
        radial_elements=5,
        tube_inner_diameter=2.25 * 0.0254,
        tube_thickness=0.15 * 0.0254,
        has_header=False,
    )

    # Unit model for 2nd reheater
    fs.aRH2 = HeatExchangerCrossFlow2D_Header(
        tube_side={"property_package": prop_water, "has_pressure_change": True},
        shell_side={"property_package": prop_gas, "has_pressure_change": True},
        finite_elements=2,
        flow_type="counter_current",
        tube_arrangement="in-line",
        tube_side_water_phase="Vap",
        has_radiation=True,
        radial_elements=5,
        tube_inner_diameter=2.25 * 0.0254,
        tube_thickness=0.15 * 0.0254,
        has_header=False,
    )

    # Unit model for primary superheater with header
    fs.aPSH = HeatExchangerCrossFlow2D_Header(
        tube_side={"property_package": prop_water, "has_pressure_change": True},
        shell_side={"property_package": prop_gas, "has_pressure_change": True},
        finite_elements=6,
        flow_type="counter_current",
        tube_arrangement="in-line",
        tube_side_water_phase="Vap",
        has_radiation=True,
        radial_elements=5,
        tube_inner_diameter=1.5 * 0.0254,
        tube_thickness=0.16 * 0.0254,
        header_radial_elements=5,
        header_inner_diameter=12 * 0.0254,
        header_wall_thickness=1.35 * 0.0254,
    )

    # Unit model for economizer
    fs.aECON = HeatExchangerCrossFlow2D_Header(
        tube_side={"property_package": prop_water, "has_pressure_change": True},
        shell_side={"property_package": prop_gas, "has_pressure_change": True},
        finite_elements=5,
        flow_type="counter_current",
        tube_arrangement="in-line",
        tube_side_water_phase="Liq",
        has_radiation=False,
        radial_elements=5,
        tube_inner_diameter=1.5 * 0.0254,
        tube_thickness=0.15 * 0.0254,
        has_header=False,
    )

    # Unit model for water pipe from economizer outlet to drum
    fs.aPipe = WaterPipe(
        dynamic=False,
        property_package=prop_water,
        has_holdup=True,
        has_heat_transfer=False,
        has_pressure_change=True,
        water_phase="Liq",
        contraction_expansion_at_end="None",
    )

    # Unit model for a mixer to mix hot primary air with tempering air
    fs.Mixer_PA = Mixer(
        dynamic=False,
        property_package=prop_gas,
        momentum_mixing_type=MomentumMixingType.equality,
        inlet_list=["PA_inlet", "TA_inlet"],
    )

    # Unit model for attemperator for main steam before platen SH
    fs.Attemp = HelmMixer(
        dynamic=False,
        property_package=prop_water,
        momentum_mixing_type=MomentumMixingType.equality,
        inlet_list=["Steam_inlet", "Water_inlet"],
    )

    # Unit model for air preheater as three-stream heat exchanger
    # with heat loss to ambient
    # side_1: flue gas
    # side_2: priamry air
    # side_3: secondry air
    fs.aAPH = HeatExchangerWith3Streams(
        dynamic=False,
        side_1_property_package=prop_gas,
        side_2_property_package=prop_gas,
        side_3_property_package=prop_gas,
        has_heat_transfer=True,
        has_pressure_change=True,
        has_holdup=False,
        flow_type_side_2="counter-current",
        flow_type_side_3="counter-current",
    )
    return m


def set_arcs_and_constraints(m):
    """
    Function to connect streams on the sub-flowsheet using pyomo arc
    and specify constraints on the sub-flowsheet
    """
    # Obtain boiler sub-flowsheet and ideal gas property
    fs = m.fs_main.fs_blr
    prop_gas = m.fs_main.prop_gas

    # Declare arcs
    fs.B001 = Arc(source=fs.aECON.tube_outlet, destination=fs.aPipe.inlet)
    fs.B001b = Arc(source=fs.aPipe.outlet, destination=fs.aDrum.feedwater_inlet)
    fs.B011 = Arc(source=fs.aDrum.liquid_outlet, destination=fs.blowdown_split.inlet)
    fs.B011b = Arc(
        source=fs.blowdown_split.FW_Downcomer, destination=fs.aDowncomer.inlet
    )
    fs.B007 = Arc(source=fs.aDowncomer.outlet, destination=fs.Waterwalls[1].inlet)

    def ww_arc_rule(b, i):
        return (b.Waterwalls[i].outlet, b.Waterwalls[i + 1].inlet)

    fs.ww_arcs = Arc(range(1, 12), rule=ww_arc_rule)

    fs.B008 = Arc(
        source=fs.Waterwalls[12].outlet, destination=fs.aDrum.water_steam_inlet
    )
    fs.B002 = Arc(source=fs.aDrum.steam_outlet, destination=fs.aRoof.inlet)
    fs.B003 = Arc(source=fs.aRoof.outlet, destination=fs.aPSH.tube_inlet)
    fs.B004 = Arc(source=fs.aPSH.tube_outlet, destination=fs.Attemp.Steam_inlet)
    fs.B005 = Arc(source=fs.Attemp.outlet, destination=fs.aPlaten.inlet)
    fs.B012 = Arc(source=fs.aRH1.tube_outlet, destination=fs.aRH2.tube_inlet)
    fs.PA04 = Arc(source=fs.aAPH.side_2_outlet, destination=fs.Mixer_PA.PA_inlet)
    fs.PA05 = Arc(source=fs.Mixer_PA.outlet, destination=fs.aBoiler.primary_air_inlet)
    fs.SA03 = Arc(
        source=fs.aAPH.side_3_outlet, destination=fs.aBoiler.secondary_air_inlet
    )
    fs.FG01 = Arc(source=fs.aBoiler.flue_gas_outlet, destination=fs.aRH2.shell_inlet)
    fs.FG02 = Arc(source=fs.aRH2.shell_outlet, destination=fs.aRH1.shell_inlet)
    fs.FG03 = Arc(source=fs.aRH1.shell_outlet, destination=fs.aPSH.shell_inlet)
    fs.FG04 = Arc(source=fs.aPSH.shell_outlet, destination=fs.aECON.shell_inlet)
    fs.FG05 = Arc(source=fs.aECON.shell_outlet, destination=fs.aAPH.side_1_inlet)

    # Expand arcs. This must be called after discretization call
    pyo.TransformationFactory("network.expand_arcs").apply_to(fs)

    # Follwing are flowsheet level constraints
    #
    # Constraint to set boiler zone heat duty equal to waterwall section duty
    @fs.Constraint(fs.time, fs.ww_zones, doc="boiler zone heat duty")
    def zone_heat_loss_eqn(b, t, izone):
        return (
            1e-6 * b.aBoiler.waterwall_heat[t, izone]
            == 1e-6 * b.Waterwalls[izone].heat_fireside[t]
        )

    # Constraint to set boiler platen heat duty equal to platen model duty
    @fs.Constraint(fs.time, doc="platen SH heat duty")
    def platen_heat_loss_eqn(b, t):
        return 1e-6 * b.aBoiler.platen_heat[t] == 1e-6 * b.aPlaten.heat_fireside[t]

    # Constraint to set boiler roof heat duty equal to roof model duty
    @fs.Constraint(fs.time, doc="roof SH heat duty")
    def roof_heat_loss_eqn(b, t):
        return 1e-6 * b.aBoiler.roof_heat[t] == 1e-6 * b.aRoof.heat_fireside[t]

    # Constraint to set boiler slag layer wall temperature equal to
    # waterwall section model slag wall temperature
    @fs.Constraint(fs.time, fs.ww_zones, doc="zone wall temperature")
    def zone_wall_temp_eqn(b, t, izone):
        return (
            b.aBoiler.wall_temperature_waterwall[t, izone]
            == b.Waterwalls[izone].temp_slag_boundary[t]
        )

    # Constraint to set boiler platen slag wall temperature equal to
    # platen superheater slag wall temperature
    @fs.Constraint(fs.time, doc="platen wall temperature")
    def platen_wall_temp_eqn(b, t):
        return b.aBoiler.wall_temperature_platen[t] == b.aPlaten.temp_slag_boundary[t]

    # Constraint to set boiler roof slag wall temperature equal to
    # roof superheater slag wall temperature
    @fs.Constraint(fs.time, doc="roof wall temperature")
    def roof_wall_temp_eqn(b, t):
        return b.aBoiler.wall_temperature_roof[t] == b.aRoof.temp_slag_boundary[t]

    # Constraint to set RH steam flow as 90% of main steam flow
    # This constraint should be deactivated if this boiler sub-flowsheet
    # is combined with the steam cycle sub-flowsheet to form a full plant model
    @fs.Constraint(fs.time, doc="RH steam flow")
    def flow_mol_steam_rh_eqn(b, t):
        return b.aRH1.tube_inlet.flow_mol[t] == 0.9 * b.aPlaten.outlet.flow_mol[t]

    # Constraint to set pressure drop of PA equal to that of SA
    @fs.Constraint(fs.time, doc="Pressure drop of PA and SA of APH")
    def pressure_drop_of_APH_eqn(b, t):
        return b.aAPH.deltaP_side_2[t] == b.aAPH.deltaP_side_3[t]

    # Constraint to set the temperature of PA before APH equal to
    # the temperature of TA
    @fs.Constraint(fs.time, doc="Same inlet temperature for PA and SA")
    def pa_ta_temperature_identical_eqn(b, t):
        return b.aAPH.side_2_inlet.temperature[t] == b.Mixer_PA.TA_inlet.temperature[t]

    # Constraint to set blowdown water flow rate equal to 2% of feed water flow
    @fs.Constraint(fs.time, doc="Blowdown water flow fraction")
    def blowdown_flow_fraction_eqn(b, t):
        return (
            b.blowdown_split.FW_Blowdown.flow_mol[t]
            == 0.02 * b.aECON.tube_inlet.flow_mol[t]
        )

    # Surrogate model of UA (overall heat transfer coefficient times area)
    # for heat transfer from flue gas to primary air in air preheater
    # as a function of raw coal flow rate, which is related to load
    @fs.Constraint(fs.time, doc="UA for APH of PA")
    def ua_side_2_eqn(b, t):
        return 1e-4 * b.aAPH.ua_side_2[t] == 1e-4 * (
            -150.0 * b.aBoiler.flowrate_coal_raw[t] ** 2
            + 9400.0 * b.aBoiler.flowrate_coal_raw[t]
            + 65000
        )

    # Surrogate model of UA (overall heat transfer coefficient times area)
    # for heat transfer from flue gas to secondary air in air preheater
    # as a function of raw coal flow rate, which is related to load
    @fs.Constraint(fs.time, doc="UA for APH of SA")
    def ua_side_3_eqn(b, t):
        return 1e-5 * b.aAPH.ua_side_3[t] == 1e-5 * (
            30000 * b.aBoiler.flowrate_coal_raw[t] + 100000
        )

    # The follwing three constraints are related to the actual operating
    # curves set by the controllers or by the boiler manufacturers
    #
    # Constraints to set tempering air (unheated) to total PA
    # (heated + unheated) flow ratio as a function of raw coal flow rate
    # This is usually controlled to get a desired mill outlet temperature
    # Currently we don't have a mill model and the controller on the flowsheet
    # Usually more temperatuing air is needed at lower load
    @fs.Constraint(
        fs.time,
        prop_gas.component_list,
        doc="Fraction of tempering air as total PA flow",
    )
    def fraction_of_ta_in_total_pa_eqn(b, t, j):
        return b.Mixer_PA.PA_inlet.flow_mol_comp[t, j] == (
            (
                b.Mixer_PA.PA_inlet.flow_mol_comp[t, j]
                + b.Mixer_PA.TA_inlet.flow_mol_comp[t, j]
            )
            * (
                0.0003 * b.aBoiler.flowrate_coal_raw[t] ** 2
                - 0.0175 * b.aBoiler.flowrate_coal_raw[t]
                + 0.5
            )
        )

    # Constraint to set total PA flow to coal flow ratio as a function of
    # raw coal flow rate
    # This is related to mill curve and burner out of service scheduling
    @fs.Constraint(fs.time, doc="PA to coal ratio")
    def pa_to_coal_ratio_eqn(b, t):
        return b.aBoiler.ratio_PA2coal[t] == (
            0.0018 * b.aBoiler.flowrate_coal_raw[t] ** 2
            - 0.11 * b.aBoiler.flowrate_coal_raw[t]
            + 3.45
        )

    # Constraint to set dry O2 mole percent in flue gas
    # as a function of coal flow rate
    # This is related to scheduling by the controller to adjust to obtain
    # a desired main steam temperature, unburned carbon and NOx
    @fs.Constraint(fs.time, doc="Steady state dry O2 in flue gas")
    def dry_o2_in_flue_gas_eqn(b, t):
        return b.aBoiler.fluegas_o2_pct_dry[t] == (
            -0.00075 * b.aBoiler.flowrate_coal_raw[t] ** 3
            + 0.067 * b.aBoiler.flowrate_coal_raw[t] ** 2
            - 2.0 * b.aBoiler.flowrate_coal_raw[t]
            + 22.95
        )

    # Boiler efficiency based on enthalpy increase of main and RH steams
    @fs.Expression(fs.time, doc="boiler efficiency based on steam")
    def boiler_efficiency_steam(b, t):
        return (
            (b.aPlaten.outlet.flow_mol[t] - b.Attemp.Water_inlet.flow_mol[t])
            * (b.aPlaten.outlet.enth_mol[t] - b.aECON.tube_inlet.enth_mol[t])
            + b.aRH2.tube_outlet.flow_mol[t]
            * (b.aRH2.tube_outlet.enth_mol[t] - b.aRH1.tube_inlet.enth_mol[t])
            + b.Attemp.Water_inlet.flow_mol[t]
            * (b.aPlaten.outlet.enth_mol[t] - b.Attemp.Water_inlet.enth_mol[t])
        ) / (
            b.aBoiler.flowrate_coal_raw[t]
            * (1 - b.aBoiler.mf_H2O_coal_raw[t])
            * b.aBoiler.hhv_coal_dry
        )

    # Boiler efficiency based on heat absorbed
    @fs.Expression(fs.time, doc="boiler efficiency based on heat")
    def boiler_efficiency_heat(b, t):
        return (
            b.aBoiler.heat_total[t]
            + b.aRH2.total_heat[t]
            + b.aRH1.total_heat[t]
            + b.aPSH.total_heat[t]
            + b.aECON.total_heat[t]
        ) / (
            b.aBoiler.flowrate_coal_raw[t]
            * (1 - b.aBoiler.mf_H2O_coal_raw[t])
            * b.aBoiler.hhv_coal_dry
        )

    return m


def set_inputs(m):
    """
    Function to set the basic inputs for the unit models
    including geometry and fixed design and operating variables
    """
    fs = m.fs_main.fs_blr
    # Specify air composition (mole fractions)
    # based on 298.15 K and 0.5 relative humidity
    # It is defined in fire-side boiler model but used for all air inlets
    fs.aBoiler.mole_frac_air["O2"] = 0.206201
    fs.aBoiler.mole_frac_air["N2"] = 0.777811
    fs.aBoiler.mole_frac_air["CO2"] = 0.0003346
    fs.aBoiler.mole_frac_air["H2O"] = 0.0156532
    fs.aBoiler.mole_frac_air["SO2"] = 0.0000001
    fs.aBoiler.mole_frac_air["NO"] = 0.0000001

    # Coal fuel analysis data, always fix in the current model
    fs.aBoiler.mf_C_coal_dry.fix(0.717259)
    fs.aBoiler.mf_H_coal_dry.fix(0.0506301)
    fs.aBoiler.mf_O_coal_dry.fix(0.0789829)
    fs.aBoiler.mf_N_coal_dry.fix(0.0140639)
    fs.aBoiler.mf_S_coal_dry.fix(0.0282403)
    fs.aBoiler.mf_Ash_coal_dry.fix(0.110824)
    fs.aBoiler.hhv_coal_dry.fix(3.05052e007)

    # Assume mill outlet temperature is controlled and fixed
    fs.aBoiler.temperature_coal[:].fix(338.7)

    # Assume 60% of coal moisture is vaporized in mills
    fs.aBoiler.frac_moisture_vaporized[:].fix(0.6)

    # drum inputs
    fs.aDrum.drum_length.fix(16.0)
    fs.aDrum.level[:].fix(0.9)
    fs.aDrum.number_downcomer.fix(8)
    fs.aDrum.downcomer_diameter.fix(0.375)
    fs.aDrum.temperature_ambient[:].fix(300)
    fs.aDrum.insulation_thickness.fix(0.15)

    # blowdown split fraction initially set to a small value
    # it will eventually unfixed due to a flowsheet constraint
    fs.blowdown_split.split_fraction[:, "FW_Blowdown"].fix(0.001)

    # downcomer inputs
    fs.aDowncomer.diameter.fix(0.375)
    fs.aDowncomer.height.fix(45)
    fs.aDowncomer.number_downcomers.fix(8)
    fs.aDowncomer.heat_duty[:].fix(0.0)

    # inputs for 12 waterwall sections
    for i in fs.ww_zones:
        fs.Waterwalls[i].tube_diameter.fix(0.055)
        fs.Waterwalls[i].tube_thickness.fix(0.0055)
        fs.Waterwalls[i].fin_thickness.fix(0.005)
        fs.Waterwalls[i].slag_thickness[:].fix(0.001)
        fs.Waterwalls[i].fin_length.fix(0.013)
        fs.Waterwalls[i].number_tubes.fix(660)
        fs.Waterwalls[i].fcorrection_dp.fix(1.2)

    # water wall section height
    fs.Waterwalls[1].height.fix(7.0)
    fs.Waterwalls[2].height.fix(4.8)
    fs.Waterwalls[3].height.fix(2.6)
    fs.Waterwalls[4].height.fix(2.6)
    fs.Waterwalls[5].height.fix(2.65)
    fs.Waterwalls[6].height.fix(2.7)
    fs.Waterwalls[7].height.fix(2.7)
    fs.Waterwalls[8].height.fix(2.7)
    fs.Waterwalls[9].height.fix(2.75)
    fs.Waterwalls[10].height.fix(3.5)
    fs.Waterwalls[11].height.fix(5.5)
    fs.Waterwalls[12].height.fix(5.5)

    # water wall section projected area
    fs.Waterwalls[1].projected_area.fix(359.738)
    fs.Waterwalls[2].projected_area.fix(249.6)
    fs.Waterwalls[3].projected_area.fix(135.2)
    fs.Waterwalls[4].projected_area.fix(135.2)
    fs.Waterwalls[5].projected_area.fix(137.8)
    fs.Waterwalls[6].projected_area.fix(140.4)
    fs.Waterwalls[7].projected_area.fix(140.4)
    fs.Waterwalls[8].projected_area.fix(140.4)
    fs.Waterwalls[9].projected_area.fix(143)
    fs.Waterwalls[10].projected_area.fix(183.81)
    fs.Waterwalls[11].projected_area.fix(179.3)
    fs.Waterwalls[12].projected_area.fix(179.3)

    # roof
    fs.aRoof.diameter_in.fix(0.055)
    fs.aRoof.tube_thickness.fix(0.0055)
    fs.aRoof.fin_thickness.fix(0.005)
    fs.aRoof.slag_thickness[:].fix(0.001)
    fs.aRoof.fin_length.fix(0.013)
    fs.aRoof.tube_length.fix(8.8)
    fs.aRoof.number_tubes.fix(190)
    fs.aRoof.therm_cond_slag.fix(1.3)

    # platen superheater
    fs.aPlaten.diameter_in.fix(0.042)
    fs.aPlaten.tube_thickness.fix(0.0065)
    fs.aPlaten.fin_thickness.fix(0.005)
    fs.aPlaten.slag_thickness[:].fix(0.001)
    fs.aPlaten.fin_length.fix(0.01)
    fs.aPlaten.tube_length.fix(44)
    fs.aPlaten.number_tubes.fix(14 * 18)
    fs.aPlaten.therm_cond_slag.fix(1.3)

    # RH1
    fs.aRH1.pitch_x.fix(4.55 * 0.0254)
    fs.aRH1.pitch_y.fix(7.0 * 0.0254)
    fs.aRH1.tube_length_seg.fix(350 * 0.0254)
    fs.aRH1.tube_nseg.fix(4)
    fs.aRH1.tube_ncol.fix(83)
    fs.aRH1.tube_inlet_nrow.fix(3)
    fs.aRH1.delta_elevation.fix(0.0)
    fs.aRH1.therm_cond_wall = 43.0
    fs.aRH1.emissivity_wall.fix(0.6)
    fs.aRH1.dens_wall = 7800
    fs.aRH1.cp_wall = 470
    fs.aRH1.Young_modulus = 1.90e5
    fs.aRH1.Possion_ratio = 0.29
    fs.aRH1.coefficient_therm_expansion = 1.3e-5
    fs.aRH1.tube_r_fouling = 0.00017
    fs.aRH1.shell_r_fouling = 0.00088
    fs.aRH1.fcorrection_htc_tube.fix(1.0)
    fs.aRH1.fcorrection_htc_shell.fix(1.0)
    fs.aRH1.fcorrection_dp_tube.fix(5.0)
    fs.aRH1.fcorrection_dp_shell.fix(3.5)

    # RH2
    fs.aRH2.pitch_x.fix(4.55 * 0.0254)
    fs.aRH2.pitch_y.fix(14.0 * 0.0254)
    fs.aRH2.tube_length_seg.fix(420 * 0.0254)
    fs.aRH2.tube_nseg.fix(2)
    fs.aRH2.tube_ncol.fix(41)
    fs.aRH2.tube_inlet_nrow.fix(6)
    fs.aRH2.delta_elevation.fix(0.0)
    fs.aRH2.therm_cond_wall = 43.0
    fs.aRH2.emissivity_wall.fix(0.6)
    fs.aRH2.dens_wall = 7800
    fs.aRH2.cp_wall = 470
    fs.aRH2.Young_modulus = 1.90e5
    fs.aRH2.Possion_ratio = 0.29
    fs.aRH2.coefficient_therm_expansion = 1.3e-5
    fs.aRH2.tube_r_fouling = 0.00017
    fs.aRH2.shell_r_fouling = 0.00088
    fs.aRH2.fcorrection_htc_tube.fix(1.0)
    fs.aRH2.fcorrection_htc_shell.fix(1.0)
    fs.aRH2.fcorrection_dp_tube.fix(5.0)
    fs.aRH2.fcorrection_dp_shell.fix(3.5)

    # PSH
    fs.aPSH.pitch_x.fix(3.8 * 0.0254)
    fs.aPSH.pitch_y.fix(6.5 * 0.0254)
    fs.aPSH.tube_length_seg.fix(350 * 0.0254)
    fs.aPSH.tube_nseg.fix(12)
    fs.aPSH.tube_ncol.fix(90)
    fs.aPSH.tube_inlet_nrow.fix(4)
    fs.aPSH.delta_elevation.fix(5.0)
    fs.aPSH.therm_cond_wall = 49.0  # Carbon steel SA 209 T1
    fs.aPSH.dens_wall = 7800
    fs.aPSH.cp_wall = 470
    fs.aPSH.Young_modulus = 1.90e5
    fs.aPSH.Possion_ratio = 0.29
    fs.aPSH.coefficient_therm_expansion = 1.3e-5
    fs.aPSH.tube_r_fouling = 0.00017
    fs.aPSH.shell_r_fouling = 0.00088
    fs.aPSH.emissivity_wall.fix(0.7)
    fs.aPSH.fcorrection_htc_tube.fix(1.0)
    fs.aPSH.fcorrection_htc_shell.fix(1.0)
    fs.aPSH.fcorrection_dp_tube.fix(5.0)
    fs.aPSH.fcorrection_dp_shell.fix(1.25)
    fs.aPSH.temperature_ambient.fix(350)
    fs.aPSH.head_insulation_thickness.fix(0.025)

    # economizer
    fs.aECON.pitch_x.fix(3.8 * 0.0254)
    fs.aECON.pitch_y.fix(4.25 * 0.0254)
    fs.aECON.tube_length_seg.fix(350 * 0.0254)
    fs.aECON.tube_nseg.fix(18)
    fs.aECON.tube_ncol.fix(138)
    fs.aECON.tube_inlet_nrow.fix(2)
    fs.aECON.delta_elevation.fix(12.0)
    fs.aECON.therm_cond_wall = 43.0
    fs.aECON.dens_wall = 7800
    fs.aECON.cp_wall = 470
    fs.aECON.Young_modulus = 1.90e5
    fs.aECON.Possion_ratio = 0.29
    fs.aECON.coefficient_therm_expansion = 1.3e-5
    fs.aECON.tube_r_fouling = 0.00017
    fs.aECON.shell_r_fouling = 0.00088
    fs.aECON.fcorrection_htc_tube.fix(1.0)
    fs.aECON.fcorrection_htc_shell.fix(1.0)
    fs.aECON.fcorrection_dp_tube.fix(10.0)
    fs.aECON.fcorrection_dp_shell.fix(4.5)

    # APH
    fs.aAPH.ua_side_2[:].fix(170000)
    fs.aAPH.ua_side_3[:].fix(677000)
    fs.aAPH.frac_heatloss.fix(0.15)
    fs.aAPH.deltaP_side_1[:].fix(-1000)
    fs.aAPH.deltaP_side_2[:].fix(-1000)
    fs.aAPH.deltaP_side_3[:].fix(-1000)

    # 138 economizer rising tubes
    fs.aPipe.diameter.fix(0.04)
    fs.aPipe.length.fix(35)
    fs.aPipe.number_of_pipes.fix(138)
    fs.aPipe.elevation_change.fix(20)
    fs.aPipe.fcorrection_dp.fix(1)

    return m


def initialize(m):
    """Initialize unit models"""
    fs = m.fs_main.fs_blr
    prop_gas = m.fs_main.prop_gas
    outlvl = idaeslog.INFO_LOW
    _log = idaeslog.getLogger(fs.name, outlvl, tag="unit")
    solve_log = idaeslog.getSolveLogger(fs.name, outlvl, tag="unit")
    solver = get_solver(
        options={"linear_solver": "ma57", "OF_ma57_automatic_scaling": "yes"}
    )

    # set initial condition to steady-state condition for dynamic flowsheet
    if m.dynamic is True:
        fs.aDrum.set_initial_condition()
        fs.aDowncomer.set_initial_condition()
        for i in fs.ww_zones:
            fs.Waterwalls[i].set_initial_condition()
        fs.aRoof.set_initial_condition()
        fs.aPSH.set_initial_condition()
        fs.aPlaten.set_initial_condition()
        fs.aRH1.set_initial_condition()
        fs.aRH2.set_initial_condition()
        fs.aECON.set_initial_condition()
        fs.aPipe.set_initial_condition()

    # fix operating conditions
    #
    # boiler
    fs.aBoiler.mf_H2O_coal_raw[:].fix(0.1112)
    fs.aBoiler.flowrate_coal_raw[:].fix(29)
    fs.aBoiler.SR[:].fix(1.19)
    fs.aBoiler.SR_lf[:].fix(0.97)
    fs.aBoiler.ratio_PA2coal[:].fix(2.0)
    fs.aBoiler.deltaP.fix(1000)
    fs.aBoiler.wall_temperature_waterwall[:, 1].fix(641)
    fs.aBoiler.wall_temperature_waterwall[:, 2].fix(664)
    fs.aBoiler.wall_temperature_waterwall[:, 3].fix(722)
    fs.aBoiler.wall_temperature_waterwall[:, 4].fix(735)
    fs.aBoiler.wall_temperature_waterwall[:, 5].fix(744)
    fs.aBoiler.wall_temperature_waterwall[:, 6].fix(747)
    fs.aBoiler.wall_temperature_waterwall[:, 7].fix(746)
    fs.aBoiler.wall_temperature_waterwall[:, 8].fix(729)
    fs.aBoiler.wall_temperature_waterwall[:, 9].fix(716)
    fs.aBoiler.wall_temperature_waterwall[:, 10].fix(698)
    fs.aBoiler.wall_temperature_waterwall[:, 11].fix(681)
    fs.aBoiler.wall_temperature_waterwall[:, 12].fix(665)
    fs.aBoiler.wall_temperature_platen[:].fix(790)
    fs.aBoiler.wall_temperature_roof[:].fix(643)
    fs.aBoiler.fcorrection_heat_ww.fix(1.0)
    fs.aBoiler.fcorrection_heat_platen.fix(1.0)
    # boiler inlet air pressure slightly higher than atmospheric pressure
    fs.aBoiler.primary_air_inlet.pressure[:].fix(101825)
    fs.aBoiler.secondary_air_inlet.pressure[:].fix(101825)
    # initial guess for SA temperature, will be calculated later by APH model
    fs.aBoiler.secondary_air_inlet.temperature[:].fix(634)

    # drum level set as drum inner radius initially
    fs.aDrum.level[:].fix(0.9)

    # 12 waterwalls, initial guess, will be calculated by boiler model
    fs.Waterwalls[1].heat_fireside[:].fix(4.42e07)
    fs.Waterwalls[2].heat_fireside[:].fix(5.04e07)
    fs.Waterwalls[3].heat_fireside[:].fix(3.43e07)
    fs.Waterwalls[4].heat_fireside[:].fix(3.62e07)
    fs.Waterwalls[5].heat_fireside[:].fix(3.58e07)
    fs.Waterwalls[6].heat_fireside[:].fix(2.88e07)
    fs.Waterwalls[7].heat_fireside[:].fix(2.54e07)
    fs.Waterwalls[8].heat_fireside[:].fix(2.27e07)
    fs.Waterwalls[9].heat_fireside[:].fix(2.04e07)
    fs.Waterwalls[10].heat_fireside[:].fix(2.35e07)
    fs.Waterwalls[11].heat_fireside[:].fix(1.24e07)
    fs.Waterwalls[12].heat_fireside[:].fix(9.06e06)

    # roof superheater, initial guess, unfixed later
    fs.aRoof.heat_fireside[:].fix(9.92e06)

    # platen superheater, initial guess, unfixed later
    fs.aPlaten.heat_fireside[:].fix(6.69e07)

    # economizer water inlet, will be calculated if connected with steam cycle
    fs.aECON.tube_inlet.flow_mol[:].fix(12000)
    fs.aECON.tube_inlet.pressure[:].fix(1.4e7)
    fs.aECON.tube_inlet.enth_mol[:].fix(18335.7)

    # fixed RH1 inlet, will be calculated if connected with steam cycle
    fs.aRH1.tube_inlet.flow_mol[:].fix(12000)
    fs.aRH1.tube_inlet.enth_mol[:].fix(55879)
    fs.aRH1.tube_inlet.pressure[:].fix(3100000)

    fs.model_check()

    # Initialize unit operation models on the flowsheet
    # since dynamic model is initialized by copy steady-state model results,
    # it will not call those initialize functions
    #
    # tear flue gas stream between PSH and ECON
    # Use FG molar composition to set component flow rates
    fs.aECON.shell_inlet.flow_mol_comp[:, "H2O"].value = 748
    fs.aECON.shell_inlet.flow_mol_comp[:, "CO2"].value = 1054
    fs.aECON.shell_inlet.flow_mol_comp[:, "N2"].value = 5377
    fs.aECON.shell_inlet.flow_mol_comp[:, "O2"].value = 194
    fs.aECON.shell_inlet.flow_mol_comp[:, "SO2"].value = 9
    fs.aECON.shell_inlet.flow_mol_comp[:, "NO"].value = 2.6
    fs.aECON.shell_inlet.temperature[:].value = 861.3
    fs.aECON.shell_inlet.pressure[:].value = 100000
    if m.dynamic is False:
        fs.aECON.initialize(outlvl=outlvl)
        _log.info("Completed economizer initialization")
    if m.dynamic is False:
        _set_port(fs.aPipe.inlet, fs.aECON.tube_outlet)
        fs.aPipe.initialize(outlvl=outlvl)
        _log.info("Completed pipe initialization")

    # PA to APH
    flow_mol_pa = 1410
    for i in prop_gas.component_list:
        fs.aAPH.side_2_inlet.flow_mol_comp[:, i].fix(
            pyo.value(flow_mol_pa * fs.aBoiler.mole_frac_air[i])
        )
    fs.aAPH.side_2_inlet.temperature[:].fix(324.8)
    fs.aAPH.side_2_inlet.pressure[:].fix(102325)

    # SA to APH
    flow_mol_sa = 4716
    for i in prop_gas.component_list:
        fs.aAPH.side_3_inlet.flow_mol_comp[:, i].fix(
            pyo.value(flow_mol_sa * fs.aBoiler.mole_frac_air[i])
        )
    fs.aAPH.side_3_inlet.temperature[:].fix(366.2)
    fs.aAPH.side_3_inlet.pressure[:].fix(102325)

    flow_mol_ta = 721
    for i in prop_gas.component_list:
        fs.Mixer_PA.TA_inlet.flow_mol_comp[:, i].fix(
            pyo.value(flow_mol_ta * fs.aBoiler.mole_frac_air[i])
        )
    fs.Mixer_PA.TA_inlet.temperature[:].value = 324.8
    fs.Mixer_PA.TA_inlet.pressure[:].value = 101325

    if m.dynamic is False:
        _set_port(fs.aAPH.side_1_inlet, fs.aECON.shell_outlet)
        fs.aAPH.initialize(outlvl=outlvl)
        _log.info("Completed air preheater initialization")

    if m.dynamic is False:
        _set_port(fs.Mixer_PA.PA_inlet, fs.aAPH.side_2_outlet)
        fs.Mixer_PA.initialize(outlvl=outlvl)
        _log.info("Completed Mixer_PA initialization")

    if m.dynamic is False:
        _set_port(fs.aBoiler.primary_air_inlet, fs.Mixer_PA.outlet)
        _set_port(fs.aBoiler.secondary_air_inlet, fs.aAPH.side_3_outlet)
        fs.aBoiler.initialize(outlvl=outlvl)
        _log.info("Completed boiler initialization")

    # tear stream from last boiler water wall section
    fs.aDrum.water_steam_inlet.flow_mol[:].value = 133000
    fs.aDrum.water_steam_inlet.pressure[:].value = 12000000
    fs.aDrum.water_steam_inlet.enth_mol[:].value = 27832

    fs.aDrum.water_steam_inlet.flow_mol[:].fix()
    fs.aDrum.water_steam_inlet.pressure[:].fix()
    fs.aDrum.water_steam_inlet.enth_mol[:].fix()

    if m.dynamic is False:
        _set_port(fs.aDrum.feedwater_inlet, fs.aPipe.outlet)
        fs.aDrum.initialize(outlvl=outlvl)
        _log.info("Completed drum initialization")

    if m.dynamic is False:
        _set_port(fs.blowdown_split.inlet, fs.aDrum.liquid_outlet)
        fs.blowdown_split.initialize()
        _log.info("Completed blowdown_split initialization")

    if m.dynamic is False:
        _set_port(fs.aDowncomer.inlet, fs.blowdown_split.FW_Downcomer)
        fs.aDowncomer.initialize(outlvl=outlvl)
        _log.info("Completed downcomer initialization")

    if m.dynamic is False:
        _set_port(fs.Waterwalls[1].inlet, fs.aDowncomer.outlet)
        fs.Waterwalls[1].initialize(outlvl=outlvl)
        _log.info("Completed zone 1 initialization")

    for i in fs.ww_zones:
        if i > 1:
            if m.dynamic is False:
                _set_port(fs.Waterwalls[i].inlet, fs.Waterwalls[i - 1].outlet)
                fs.Waterwalls[i].initialize(outlvl=outlvl)
                _log.info("Completed zone {} initialization".format(i))

    # tear steam between RH1 and RH2
    fs.aRH2.tube_inlet.flow_mol[:].value = 12000
    fs.aRH2.tube_inlet.enth_mol[:].value = 61628
    fs.aRH2.tube_inlet.pressure[:].value = 3000000

    if m.dynamic is False:
        _set_port(fs.aRH2.shell_inlet, fs.aBoiler.flue_gas_outlet)
        # use a lower temperature to avoid convegence issue
        fs.aRH2.shell_inlet.temperature[:].value = 1350
        fs.aRH2.initialize(outlvl=outlvl)
        _log.info("Completed RH2 initialization")

    if m.dynamic is False:
        _set_port(fs.aRH1.shell_inlet, fs.aRH2.shell_outlet)
        # use a lower temperature to avoid convegence issue
        fs.aRH1.shell_inlet.temperature[:].value = 1100
        fs.aRH1.initialize(outlvl=outlvl)
        _log.info("Completed RH1 initialization")

    if m.dynamic is False:
        _set_port(fs.aRoof.inlet, fs.aDrum.steam_outlet)
        fs.aRoof.initialize(outlvl=outlvl)
        _log.info("Completed foof initialization")

    if m.dynamic is False:
        _set_port(fs.aPSH.tube_inlet, fs.aRoof.outlet)
        _set_port(fs.aPSH.shell_inlet, fs.aRH1.shell_outlet)
        # use a lower temperature to avoid convegence issue
        fs.aPSH.shell_inlet.temperature[:].value = 1000
        fs.aPSH.initialize(outlvl=outlvl)
        _log.info("Completed PSH initialization")

    # fixed inlet conditions, to be linked with steam cycle
    fs.Attemp.Water_inlet.flow_mol[:].fix(20)
    fs.Attemp.Water_inlet.pressure[:].value = 1.35e7
    fs.Attemp.Water_inlet.enth_mol[:].fix(12767)

    if m.dynamic is False:
        _set_port(fs.Attemp.Steam_inlet, fs.aPSH.tube_outlet)
        fs.Attemp.initialize(outlvl=outlvl)
        _log.info("Completed attemperator initialization")

    if m.dynamic is False:
        _set_port(fs.aPlaten.inlet, fs.Attemp.outlet)
        fs.aPlaten.initialize(outlvl=outlvl)
        _log.info("Completed platen superheater initialization")

    # Unfix variables after initialization
    # blowdown split fraction
    fs.blowdown_split.split_fraction[:, "FW_Blowdown"].unfix()
    fs.aDrum.water_steam_inlet.unfix()
    # Waterwalls[] heat duty
    for i in fs.ww_zones:
        fs.Waterwalls[i].heat_fireside[:].unfix()
    # heat duty to aPlaten and aRoof
    fs.aPlaten.heat_fireside[:].unfix()
    fs.aRoof.heat_fireside[:].unfix()
    # all wall temperatures
    for i in fs.aBoiler.zones:
        fs.aBoiler.wall_temperature_waterwall[:, i].unfix()
    fs.aBoiler.wall_temperature_platen[:].unfix()
    fs.aBoiler.wall_temperature_roof[:].unfix()
    # air and gas component flows
    for i in prop_gas.component_list:
        # SA at FD fan outlet and aAPH inlet
        fs.aAPH.side_2_inlet.flow_mol_comp[:, i].unfix()
        # PA at aAPH inlet
        fs.aAPH.side_3_inlet.flow_mol_comp[:, i].unfix()
        # Tempering air at Mixer_PA inlet
        fs.Mixer_PA.TA_inlet.flow_mol_comp[:, i].unfix()
    # SA pressure and temperature need to be unfixed
    fs.aBoiler.secondary_air_inlet.pressure[:].unfix()
    fs.aBoiler.secondary_air_inlet.temperature[:].unfix()
    # PA inlet pressure needs to be unfixed
    fs.aBoiler.primary_air_inlet.pressure[:].unfix()
    # fs.aBoiler.primary_air_inlet.temperature[:].unfix()
    # inlet flow based on constraint
    fs.aRH1.tube_inlet.flow_mol[:].unfix()

    # fix feed water pressure and enthalpy but allow flow rate to change
    fs.aECON.tube_inlet.flow_mol[:].unfix()
    fs.aBoiler.SR[:].unfix()

    fs.aBoiler.ratio_PA2coal[:].unfix()
    # due to a constraint to set SA and PA deltP the same
    fs.aAPH.deltaP_side_2[:].unfix()
    fs.aAPH.ua_side_2.unfix()
    fs.aAPH.ua_side_3.unfix()

    df = degrees_of_freedom(fs)
    _log.info("boiler flowsheet degree of freedom = {}".format(df))
    assert df == 0

    if m.dynamic is False:
        _log.info("Solving boiler steady-state problem...")
        with idaeslog.solver_log(solve_log, idaeslog.DEBUG) as slc:
            res = solver.solve(fs, tee=slc.tee)
        _log.info(
            "Solving boiler steady-state problem: {}".format(idaeslog.condition(res))
        )
        _log.info("feed water flow={}".format(fs.aECON.tube_inlet.flow_mol[0].value))
        _log.info("main steam flow={}".format(fs.aPlaten.outlet.flow_mol[0].value))
        _log.info("spray water flow={}".format(fs.Attemp.Water_inlet.flow_mol[0].value))
        _log.info(
            "flowrate_coal_raw={}".format(pyo.value(fs.aBoiler.flowrate_coal_raw[0]))
        )
        _log.info(
            "fraction of moisture in raw coal={}".format(
                pyo.value(fs.aBoiler.mf_H2O_coal_raw[0])
            )
        )
        _log.info(
            "flowrate_fluegas={}".format(pyo.value(fs.aBoiler.flue_gas[0].flow_mass))
        )
        _log.info(
            "tempering air flow={}".format(
                pyo.value(fs.Mixer_PA.TA_inlet_state[0].flow_mass)
            )
        )
        _log.info(
            "total PA flow={}".format(pyo.value(fs.aBoiler.primary_air[0].flow_mass))
        )
        _log.info("SA flow={}".format(pyo.value(fs.aBoiler.secondary_air[0].flow_mass)))
        _log.info("TCA flow={}".format(pyo.value(fs.aBoiler.flow_mass_TCA[0])))
        _log.info("ratio_PA2coal={}".format(pyo.value(fs.aBoiler.ratio_PA2coal[0])))
        _log.info("SR={}".format(fs.aBoiler.SR[0].value))
        _log.info("SR_lf={}".format(fs.aBoiler.SR_lf[0].value))
        _log.info("total ww heat={}".format(pyo.value(fs.aBoiler.heat_total_ww[0])))
        _log.info("total heat={}".format(pyo.value(fs.aBoiler.heat_total_ww[0])))
        _log.info(
            "FEGT={}".format(pyo.value(fs.aBoiler.flue_gas_outlet.temperature[0]))
        )
    return m


def set_scaling_factors(m):
    """Set scaling factors for variables and expressions. These are used for
    variable scaling and used by the framework to scale constraints.

    Args:
        m: plant model to set scaling factors for.

    Returns:
        None
    """

    fs = m.fs_main.fs_blr

    for i, ww in fs.Waterwalls.items():
        iscale.set_scaling_factor(ww.heat_fireside, 1e-7)
        if i < 4:
            iscale.set_scaling_factor(ww.heat_flux_conv, 1e-4)
        else:
            iscale.set_scaling_factor(ww.heat_flux_conv, 1e-5)
        iscale.set_scaling_factor(ww.tube_diameter, 100)
        iscale.set_scaling_factor(ww.control_volume.material_holdup, 1e-4)
        iscale.set_scaling_factor(ww.control_volume.energy_holdup, 1e-8)
        iscale.set_scaling_factor(ww.energy_holdup_slag, 1e-3)
        iscale.set_scaling_factor(ww.energy_holdup_metal, 1e-6)
        iscale.set_scaling_factor(ww.N_Re, 1e-6)
        iscale.set_scaling_factor(ww.pitch, 1e3)
        for j, c in ww.hconv_lo_eqn.items():
            iscale.constraint_scaling_transform(c, 1e-2)

    iscale.set_scaling_factor(fs.aRoof.heat_fireside, 1e-6)
    iscale.set_scaling_factor(fs.aRoof.heat_flux_conv, 1e-4)
    iscale.set_scaling_factor(fs.aRoof.hconv, 1e-3)
    iscale.set_scaling_factor(fs.aRoof.deltaP, 1e-3)
    iscale.set_scaling_factor(fs.aRoof.diameter_in, 100)
    iscale.set_scaling_factor(fs.aRoof.N_Re, 1e-6)

    iscale.set_scaling_factor(fs.aPlaten.heat_fireside, 1e-7)
    iscale.set_scaling_factor(fs.aPlaten.heat_flux_conv, 1e-4)
    iscale.set_scaling_factor(fs.aPlaten.hconv, 1e-3)
    iscale.set_scaling_factor(fs.aPlaten.deltaP, 1e-3)
    iscale.set_scaling_factor(fs.aPlaten.diameter_in, 100)
    iscale.set_scaling_factor(fs.aPlaten.N_Re, 1e-6)
    iscale.set_scaling_factor(fs.aDrum.control_volume.energy_holdup, 1e-10)
    iscale.set_scaling_factor(fs.aDrum.control_volume.material_holdup, 1e-5)
    if m.dynamic:
        for t, c in fs.aDrum.control_volume.energy_accumulation_disc_eq.items():
            iscale.constraint_scaling_transform(c, 1e-4)

    iscale.set_scaling_factor(fs.aDowncomer.control_volume.energy_holdup, 1e-10)
    iscale.set_scaling_factor(fs.aDowncomer.control_volume.material_holdup, 1e-5)
    iscale.set_scaling_factor(fs.aRoof.control_volume.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.aRoof.control_volume.material_holdup, 1e-6)
    iscale.set_scaling_factor(fs.aPlaten.control_volume.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.aPlaten.control_volume.material_holdup, 1e-6)
    iscale.set_scaling_factor(fs.aPipe.control_volume.energy_holdup, 1e-9)
    iscale.set_scaling_factor(fs.aPipe.control_volume.material_holdup, 1e-5)
    iscale.set_scaling_factor(fs.aBoiler.waterwall_heat, 1e-7)
    iscale.set_scaling_factor(fs.aBoiler.platen_heat, 1e-7)
    iscale.set_scaling_factor(fs.aBoiler.roof_heat, 1e-7)
    iscale.set_scaling_factor(fs.aECON.shell._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(fs.aECON.tube._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(fs.aECON.shell.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(fs.aECON.tube.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(fs.aECON.shell.heat, 1e-7)
    iscale.set_scaling_factor(fs.aECON.tube.heat, 1e-7)
    for t, c in fs.aECON.shell.enthalpy_flow_dx_disc_eq.items():
        iscale.constraint_scaling_transform(c, 1e-7)
    for t, c in fs.aECON.tube.enthalpy_flow_dx_disc_eq.items():
        iscale.constraint_scaling_transform(c, 1e-7)

    iscale.set_scaling_factor(fs.aPSH.shell._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(fs.aPSH.tube._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(fs.aPSH.shell.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(fs.aPSH.tube.enthalpy_flow_dx, 1e-7)
    for t, c in fs.aPSH.shell.enthalpy_flow_dx_disc_eq.items():
        iscale.constraint_scaling_transform(c, 1e-7)
    for t, c in fs.aPSH.tube.enthalpy_flow_dx_disc_eq.items():
        iscale.constraint_scaling_transform(c, 1e-7)

    iscale.set_scaling_factor(fs.aRH1.shell._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(fs.aRH1.tube._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(fs.aRH1.shell.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(fs.aRH1.tube.enthalpy_flow_dx, 1e-7)
    for t, c in fs.aRH1.shell.enthalpy_flow_dx_disc_eq.items():
        iscale.constraint_scaling_transform(c, 1e-7)
    for t, c in fs.aRH1.tube.enthalpy_flow_dx_disc_eq.items():
        iscale.constraint_scaling_transform(c, 1e-7)

    iscale.set_scaling_factor(fs.aRH2.shell._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(fs.aRH2.tube._enthalpy_flow, 1e-8)
    iscale.set_scaling_factor(fs.aRH2.shell.enthalpy_flow_dx, 1e-7)
    iscale.set_scaling_factor(fs.aRH2.tube.enthalpy_flow_dx, 1e-7)
    for t, c in fs.aRH2.shell.enthalpy_flow_dx_disc_eq.items():
        iscale.constraint_scaling_transform(c, 1e-7)
    for t, c in fs.aRH2.tube.enthalpy_flow_dx_disc_eq.items():
        iscale.constraint_scaling_transform(c, 1e-7)

    # Calculate calculated scaling factors
    iscale.calculate_scaling_factors(m)


def main_steady_state():
    m_ss = get_model(dynamic=False)
    aBoiler = m_ss.fs_main.fs_blr.aBoiler
    outlvl = idaeslog.DEBUG
    _log = idaeslog.getLogger(aBoiler.name, outlvl, tag="flowsheet")
    _log.debug("PA flow={}".format(pyo.value(aBoiler.primary_air[0].flow_mass)))
    _log.debug("SA flow={}".format(pyo.value(aBoiler.secondary_air[0].flow_mass)))
    _log.debug("O2 in flue gas{}".format(pyo.value(aBoiler.fluegas_o2_pct_dry[0])))
    _log.debug("FEGT={}".format(pyo.value(aBoiler.flue_gas[0].temperature)))
    _log.debug("heat_total={}".format(pyo.value(aBoiler.heat_total[0])))
    _log.debug("total heat_ww={}".format(pyo.value(aBoiler.heat_total_ww[0])))
    _log.debug("UBC={}".format(pyo.value(aBoiler.ubc_in_flyash[0])))
    _log.debug("NOx={}".format(pyo.value(aBoiler.frac_mol_NOx_fluegas[0])))
    _log.debug(
        "main steam flow={}".format(
            pyo.value(m_ss.fs_main.fs_blr.aPlaten.outlet.flow_mol[0])
        )
    )
    _log.debug(
        "main steam pressure={}".format(
            pyo.value(m_ss.fs_main.fs_blr.aPlaten.outlet.pressure[0])
        )
    )
    _log.debug(
        "main steam temperature={}".format(
            pyo.value(
                m_ss.fs_main.fs_blr.aPlaten.control_volume.properties_out[0].temperature
            )
        )
    )
    _log.debug(
        "RH steam flow={}".format(
            pyo.value(m_ss.fs_main.fs_blr.aRH2.tube_outlet.flow_mol[0])
        )
    )
    _log.debug(
        "RH steam pressure={}".format(
            pyo.value(m_ss.fs_main.fs_blr.aRH2.tube_outlet.pressure[0])
        )
    )
    _log.debug(
        "RH steam temperature={}".format(
            pyo.value(m_ss.fs_main.fs_blr.aRH2.tube.properties[0, 0].temperature)
        )
    )
    _log.debug(
        "flue gas econ inlet temp={}".format(
            pyo.value(m_ss.fs_main.fs_blr.aECON.shell.properties[0, 0].temperature)
        )
    )
    _log.debug(
        "flue gas econ out temp={}".format(
            pyo.value(m_ss.fs_main.fs_blr.aECON.shell.properties[0, 1].temperature)
        )
    )
    _log.debug(
        "flue gas aph out temp={}".format(
            pyo.value(m_ss.fs_main.fs_blr.aAPH.side_1.properties_out[0].temperature)
        )
    )
    _log.debug(
        "econ out water temp={}".format(
            pyo.value(m_ss.fs_main.fs_blr.aECON.tube.properties[0, 0].temperature)
        )
    )
    _log.debug(
        "econ out water sat temp={}".format(
            pyo.value(m_ss.fs_main.fs_blr.aECON.tube.properties[0, 0].temperature_sat)
        )
    )
    return m_ss


def main_dynamic():
    m_ss = get_model(dynamic=False)
    # aBoiler = m_ss.fs_main.fs_blr.aBoiler
    m_dyn = get_model(dynamic=True)
    copy_non_time_indexed_values(
        m_dyn.fs_main, m_ss.fs_main, copy_fixed=True, outlvl=idaeslog.ERROR
    )
    for t in m_dyn.fs_main.time:
        copy_values_at_time(
            m_dyn.fs_main, m_ss.fs_main, t, 0.0, copy_fixed=True, outlvl=idaeslog.ERROR
        )

    solver = get_solver()

    dof = degrees_of_freedom(m_dyn.fs_main)
    # solving dynamic model at steady-state
    print("solving dynamic model at steady-state...")
    solver.solve(m_dyn.fs_main, tee=True)

    # run dynamic model
    run_dynamic(m_dyn)
    return m_dyn


def get_model(dynamic=True, init=True):
    m = pyo.ConcreteModel()
    m.dynamic = dynamic
    m.init_dyn = False
    if m.dynamic:
        m.fs_main = FlowsheetBlock(
            dynamic=True, time_set=[0, 60], time_units=pyo.units.s
        )
    else:
        m.fs_main = FlowsheetBlock(dynamic=False)
    # Add property packages to flowsheet library
    m.fs_main.prop_water = iapws95.Iapws95ParameterBlock()
    m.fs_main.prop_gas = FlueGasParameterBlock()
    m.fs_main.fs_blr = FlowsheetBlock(time_units=pyo.units.s)
    m = add_unit_models(m)
    if m.dynamic:
        m.discretizer = pyo.TransformationFactory("dae.finite_difference")
        m.discretizer.apply_to(m, nfe=2, wrt=m.fs_main.time, scheme="BACKWARD")
    m = set_arcs_and_constraints(m)
    m = set_inputs(m)
    set_scaling_factors(m)
    if init:
        m = initialize(m)

    return m


def run_dynamic(m):
    print("//////////////////start dynamic solver///////////////")
    fs = m.fs_main.fs_blr

    # Create a disturbance
    for t in m.fs_main.time:
        if t >= 30:
            fs.aBoiler.flowrate_coal_raw[t].fix(
                fs.aBoiler.flowrate_coal_raw[0].value * 1.025
            )
        else:
            fs.aBoiler.flowrate_coal_raw[t].fix(fs.aBoiler.flowrate_coal_raw[0].value)
    df = degrees_of_freedom(m)
    assert df == 0
    solver = get_solver()

    solver.solve(m, tee=True)


def print_dynamic_results(m):
    fs = m.fs_main.fs_blr

    # ploting results
    time = []
    coal_flow = []
    steam_flow = []
    steam_temp = []
    steam_pres = []
    fegt = []
    drum_level = []

    for t in m.fs_main.time:
        time.append(t)
        coal_flow.append(fs.aBoiler.flowrate_coal_raw[t].value)
        steam_flow.append(fs.aPlaten.outlet.flow_mol[t].value / 1000)
        steam_temp.append(
            pyo.value(fs.aPlaten.control_volume.properties_out[t].temperature)
        )
        steam_pres.append(fs.aPlaten.outlet.pressure[t].value / 1e6)
        fegt.append(fs.aBoiler.flue_gas_outlet.temperature[t].value)
        drum_level.append(fs.aDrum.level[t].value)

    plt.figure(1)
    plt.plot(time, coal_flow)
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Coal Flow Rate [kg/s]")
    plt.show(block=False)

    plt.figure(2)
    plt.plot(time, steam_flow)
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Steam Flow Rate [kmol/s]")
    plt.show(block=False)

    plt.figure(3)
    plt.plot(time, steam_temp)
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Steam Temperature [K]")
    plt.show(block=False)

    plt.figure(4)
    plt.plot(time, steam_pres)
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Steam Pressure [MPa]")
    plt.show(block=False)

    plt.figure(5)
    plt.plot(time, fegt)
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("FEGT [K]")
    plt.show(block=True)


if __name__ == "__main__":
    # This method builds and runs a steady state subcritical coal fired boiler
    # flowsheet. The subcritical flowsheet consists mainly of:
    # water/steam route: Economizer, water pipe, drum, blowdown splitter,
    # downcomer, water wall tubes, drum, roof backpas superheater, primary sh,
    # platen superheater, the main steam is connected with the HP turbine.
    # HP turbine exit enters the Rehater 1, and leaves to IP Turbine.
    # fluegas route: Boiler fireside, reheater, primary superheater, economizer
    # air preheater, and stack.
    # fixed inputs: Water inlet to economizer (from feed water heaters)
    # Steam to Reheaters (from HP turbine exit)
    m = main_steady_state()

    # This method builds a steady state boiler, then builds a dynamic flowsheet
    # copies the steady state solution to initialize the dynamic model and
    # solves the dynamic simulation of a subcritical boiler with a step
    # function for the coal flowrate in the boiler
    # m = main_dynamic()
    # print_dynamic_results(m)
