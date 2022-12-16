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

import os
import time as wall_clock
import matplotlib.pyplot as plt

import pyomo.environ as pyo
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.models.properties import iapws95
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.core.util.scaling as iscale
from idaes.models_extra.power_generation.properties import FlueGasParameterBlock
from idaes.models.control.controller import (
    PIDController,
    ControllerType,
)
import idaes.models_extra.power_generation.flowsheets.subcritical_power_plant.subcritical_boiler_flowsheet as blr
import idaes.models_extra.power_generation.flowsheets.subcritical_power_plant.steam_cycle_flowsheet as stc
from idaes.core.util.dyn_utils import copy_values_at_time, copy_non_time_indexed_values
import idaes.logger as idaeslog
import os
import idaes.core.util.tables as tables
from idaes.core.util.tags import svg_tag, ModelTagGroup
from idaes.core.solvers import get_solver

_log = idaeslog.getLogger(__name__)

__author__ = "J. Ma, M. Wang, M. Zamarripa, J. Eslick"


def set_scaling_factors(m):
    """Set scaling factors for variables and expressions. These are used for
    variable scaling and used by the framework to scale constraints.

    Args:
        m: plant model to set scaling factors for.

    Returns:
        None
    """
    # Set scaling factors for boiler system
    fs = m.fs_main.fs_blr

    for i, ww in fs.Waterwalls.items():
        iscale.set_scaling_factor(ww.heat_fireside, 1e-7)
        if i < 3:
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

    iscale.set_scaling_factor(fs.aRoof.heat_fireside, 1e-7)
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
            iscale.constraint_scaling_transform(c, 1e-5)
        for t, c in fs.aDrum.heat_loss_eqn.items():
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

    # Set scale factors for steam cycle
    fs = m.fs_main.fs_stc

    iscale.set_scaling_factor(fs.condenser.hot_side.heat, 1e-9)
    iscale.set_scaling_factor(fs.condenser.cold_side.heat, 1e-9)

    iscale.set_scaling_factor(fs.aux_condenser.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(fs.aux_condenser.cold_side.heat, 1e-7)

    iscale.set_scaling_factor(fs.hotwell_tank.control_volume.energy_holdup, 1e-10)
    iscale.set_scaling_factor(fs.hotwell_tank.control_volume.material_holdup, 1e-6)
    if m.dynamic:
        for t, c in fs.hotwell_tank.control_volume.energy_accumulation_disc_eq.items():
            iscale.constraint_scaling_transform(c, 1e-6)

    iscale.set_scaling_factor(fs.fwh1.condense.hot_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh1.condense.hot_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh1.condense.cold_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh1.condense.cold_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh1.condense.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(fs.fwh1.condense.cold_side.heat, 1e-7)

    iscale.set_scaling_factor(fs.fwh2.condense.hot_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh2.condense.hot_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh2.condense.cold_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh2.condense.cold_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh2.condense.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(fs.fwh2.condense.cold_side.heat, 1e-7)

    iscale.set_scaling_factor(fs.fwh3.condense.hot_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh3.condense.hot_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh3.condense.cold_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh3.condense.cold_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh3.condense.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(fs.fwh3.condense.cold_side.heat, 1e-7)

    iscale.set_scaling_factor(fs.da_tank.control_volume.energy_holdup, 1e-11)
    iscale.set_scaling_factor(fs.da_tank.control_volume.material_holdup, 1e-6)
    if m.dynamic:
        for t, c in fs.da_tank.control_volume.energy_accumulation_disc_eq.items():
            iscale.constraint_scaling_transform(c, 1e-7)

    iscale.set_scaling_factor(fs.fwh5.condense.hot_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh5.condense.hot_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh5.condense.cold_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh5.condense.cold_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh5.condense.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(fs.fwh5.condense.cold_side.heat, 1e-7)

    iscale.set_scaling_factor(fs.fwh6.condense.hot_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh6.condense.hot_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh6.condense.cold_side.material_holdup, 1e-4)
    iscale.set_scaling_factor(fs.fwh6.condense.cold_side.energy_holdup, 1e-8)
    iscale.set_scaling_factor(fs.fwh6.condense.hot_side.heat, 1e-7)
    iscale.set_scaling_factor(fs.fwh6.condense.cold_side.heat, 1e-7)

    # scaling factor for control valves
    for t in m.fs_main.time:
        iscale.set_scaling_factor(
            fs.spray_valve.control_volume.properties_in[t].flow_mol, 0.05
        )
        iscale.set_scaling_factor(
            fs.makeup_valve.control_volume.properties_in[t].flow_mol, 0.01
        )
        iscale.set_scaling_factor(
            fs.fwh2_valve.control_volume.properties_in[t].flow_mol, 0.001
        )
        iscale.set_scaling_factor(
            fs.fwh3_valve.control_volume.properties_in[t].flow_mol, 0.001
        )
        iscale.set_scaling_factor(
            fs.fwh5_valve.control_volume.properties_in[t].flow_mol, 0.001
        )
        iscale.set_scaling_factor(
            fs.fwh6_valve.control_volume.properties_in[t].flow_mol, 0.001
        )
        iscale.set_scaling_factor(
            fs.bfp_turb_valve.control_volume.properties_in[t].flow_mol, 0.01
        )

    # Calculate calculated scaling factors
    iscale.calculate_scaling_factors(m)


def add_overall_performance_expressions(m):
    @m.fs_main.Expression(m.fs_main.time, doc="Heat input to the steam cycle (W)")
    def boiler_heat(b, t):
        return (
            (
                b.fs_stc.turb.inlet_split.mixed_state[t].flow_mol
                - b.fs_stc.spray_valve.outlet.flow_mol[t]
            )
            * (
                b.fs_stc.turb.inlet_split.mixed_state[t].enth_mol
                - b.fs_stc.fwh6.desuperheat.cold_side.properties_out[t].enth_mol
            )
            + b.fs_stc.spray_valve.outlet.flow_mol[t]
            * (
                b.fs_stc.turb.inlet_split.mixed_state[t].enth_mol
                - b.fs_stc.spray_valve.outlet.enth_mol[t]
            )
            + b.fs_stc.turb.ip_stages[1].control_volume.properties_in[t].flow_mol
            * (
                b.fs_stc.turb.ip_stages[1].control_volume.properties_in[t].enth_mol
                - b.fs_stc.turb.hp_split[14].outlet_1.enth_mol[t]
            )
        )

    # Calculate steam cycle efficiency
    @m.fs_main.Expression(m.fs_main.time)
    def steam_cycle_eff(b, t):
        return -b.fs_stc.turb.power[t] / b.boiler_heat[t]

    # Calculate gross heat rate of the plant.
    @m.fs_main.Expression(m.fs_main.time, doc="Heat rate based on gross power (BTU/MW)")
    def gross_heat_rate(b, t):
        return (
            -b.fs_blr.aBoiler.flowrate_coal_raw[t]
            * (1 - b.fs_blr.aBoiler.mf_H2O_coal_raw[t])
            * b.fs_blr.aBoiler.hhv_coal_dry
            * 3600
            / b.fs_stc.turb.power[t]
            * 1000
            / 1054
        )

    # Calculate the overall efficiency of the plant.
    @m.fs_main.Expression(
        m.fs_main.time, doc="Overall efficency based on gross power (%)"
    )
    def plant_gross_efficiency(b, t):
        return -b.fs_stc.turb.power[t] / (
            b.fs_blr.aBoiler.flowrate_coal_raw[t]
            * (1 - b.fs_blr.aBoiler.mf_H2O_coal_raw[t])
            * b.fs_blr.aBoiler.hhv_coal_dry
        )

    # Calculate total auxillary power based on a surrogate expression
    @m.fs_main.Expression(m.fs_main.time)
    def aux_power(b, t):
        steam_f = m.fs_main.fs_stc.turb.inlet_split.mixed_state[t].flow_mass * 7.937
        air_f = m.fs_main.fs_blr.aBoiler.primary_air[t].flow_mass * 7.937
        cw_f = (
            m.fs_main.fs_stc.condenser.tube.properties_in[t].flow_mass
            / 997
            * 15850
            / 1000
        )
        steam_t = (
            m.fs_main.fs_blr.aPlaten.control_volume.properties_out[t].temperature
            - 273.15
        ) * 9 / 5 + 32
        steam_p = (
            m.fs_main.fs_blr.aPlaten.control_volume.properties_out[t].pressure - 79410
        ) / 6895
        air_p = m.fs_main.fs_blr.aAPH.side_2_inlet.pressure[t] / 1000
        return (
            -7.743210e-2 * cw_f
            + 2.424004e-7 * steam_f**2
            + 1.718702e-6 * steam_f * air_f
            + 5.301540e-11 * steam_f**3
            + 2.806236e-10 * steam_f**2 * steam_t
            + 7.943112e-10 * steam_f**2 * air_f
            + 5.106800e-10 * steam_f * steam_t * air_f
            + 6.356823e-10 * steam_f * air_f**2
            + -2.094381e-7 * steam_p * cw_f**2
            + 5.849115e-10 * steam_t**3
            + 1.049080e-11 * steam_t**2 * air_f
            + 1.913389e-9 * steam_t * air_p**2
            + 15.19445
        ) * 1e6


def main_steady_state():
    m = get_model(dynamic=False)
    return m


def input_profile(t, x0):
    """
    calculate the user input x as a function of time
    x0 is the initial value
    ramp down from 100% to 50% load at 5%/min rate
    hold for half an hour
    ramp up from 50% to 100% load at 5%/min rate
    hold for 20 minutes
    """
    if t < 60:
        x = x0
    elif t < 660:
        x = x0 * (1 - (t - 60) / 600 * 0.5)
    elif t < 2460:
        x = x0 * 0.5
    elif t < 3060:
        x = x0 * (0.5 + (t - 2460) / 600 * 0.5)
    else:
        # Hold for 1200 sec to 4260 sec
        x = x0
    return x


def main_dynamic():
    """Main function to set up and run dynamic flowsheet simulation"""
    start_time = wall_clock.time()
    # Declare dictionary for the result data
    plot_data = {}
    plot_data["time"] = []
    plot_data["coal_flow"] = []
    plot_data["PA_to_APH_flow"] = []
    plot_data["PA_temp_flow"] = []
    plot_data["SA_flow"] = []
    plot_data["SR"] = []
    plot_data["dry_O2_pct_in_flue_gas"] = []
    plot_data["bfpt_opening"] = []
    plot_data["gross_power"] = []
    plot_data["ww_heat"] = []
    plot_data["fegt"] = []
    plot_data["drum_level"] = []
    plot_data["feed_water_flow_sp"] = []
    plot_data["drum_master_ctrl_op"] = []
    plot_data["feed_water_flow"] = []
    plot_data["spray_flow"] = []
    plot_data["main_steam_flow"] = []
    plot_data["rh_steam_flow"] = []
    plot_data["bfpt_flow"] = []
    plot_data["main_steam_temp"] = []
    plot_data["rh_steam_temp"] = []
    plot_data["fw_pres"] = []
    plot_data["drum_pres"] = []
    plot_data["main_steam_pres"] = []
    plot_data["rh_steam_pres"] = []
    plot_data["hw_tank_level"] = []
    plot_data["da_tank_level"] = []
    plot_data["fwh2_level"] = []
    plot_data["fwh3_level"] = []
    plot_data["fwh5_level"] = []
    plot_data["fwh6_level"] = []
    plot_data["makeup_valve_opening"] = []
    plot_data["cond_valve_opening"] = []
    plot_data["fwh2_valve_opening"] = []
    plot_data["fwh3_valve_opening"] = []
    plot_data["fwh5_valve_opening"] = []
    plot_data["fwh6_valve_opening"] = []
    plot_data["spray_valve_opening"] = []
    plot_data["tube_temp_rh2"] = []
    plot_data["temp_fg_econ_exit"] = []
    plot_data["temp_fg_aph_exit"] = []
    plot_data["throttle_opening"] = []
    plot_data["load_demand"] = []
    plot_data["sliding_pressure"] = []
    plot_data["steam_pressure_sp"] = []
    plot_data["deaerator_pressure"] = []
    plot_data["temp_econ_in"] = []
    plot_data["temp_econ_out"] = []
    plot_data["temp_econ_out_sat"] = []
    plot_data["boiler_efficiency_heat"] = []
    plot_data["boiler_efficiency_steam"] = []
    plot_data["steam_cycle_efficiency"] = []
    plot_data["plant_gross_efficiency"] = []
    plot_data["gross_heat_rate"] = []

    # Boiler health related data
    # For drum
    plot_data["drum_inner_wall_temperature"] = []
    plot_data["drum_outer_wall_temperature"] = []
    plot_data["drum_inner_theta_sigma_m"] = []
    plot_data["drum_inner_theta_sigma_t"] = []
    plot_data["drum_inner_vM_sigma"] = []
    plot_data["drum_outer_theta_sigma_m"] = []
    plot_data["drum_outer_theta_sigma_t"] = []
    plot_data["drum_outer_vM_sigma"] = []
    plot_data["drum_delta_sigma_r_theta_outer"] = []
    plot_data["drum_delta_sigma_theta_z_outer"] = []
    plot_data["drum_delta_sigma_z_r_outer"] = []
    # Stress at opening of the drum
    plot_data["drum_circumferential_P1"] = []
    plot_data["drum_circumferential_P2"] = []
    plot_data["drum_von_Mises_P1"] = []
    plot_data["drum_von_Mises_P2"] = []

    # For primary superheater
    plot_data["PSH_flue_gas_flow_disc1"] = []
    plot_data["PSH_steam_flow_disc9"] = []
    plot_data["PSH_steam_temperature_disc1"] = []
    plot_data["PSH_steam_temperature_disc2"] = []
    plot_data["PSH_steam_temperature_disc4"] = []
    plot_data["PSH_steam_temperature_disc8"] = []
    plot_data["PSH_steam_temperature_disc9"] = []
    plot_data["PSH_temp_wall_tube_disc1"] = []
    plot_data["PSH_temp_wall_tube_disc2"] = []
    plot_data["PSH_temp_wall_tube_disc4"] = []
    plot_data["PSH_temp_wall_tube_disc8"] = []
    plot_data["PSH_temp_wall_tube_disc9"] = []
    plot_data["PSH_temp_wall_shell_disc1"] = []
    plot_data["PSH_temp_wall_shell_disc2"] = []
    plot_data["PSH_temp_wall_shell_disc4"] = []
    plot_data["PSH_temp_wall_shell_disc8"] = []
    plot_data["PSH_temp_wall_shell_disc9"] = []
    plot_data["PSH_temp_wall_shell_fouling_disc1"] = []
    plot_data["PSH_temp_wall_shell_fouling_disc2"] = []
    plot_data["PSH_temp_wall_shell_fouling_disc4"] = []
    plot_data["PSH_temp_wall_shell_fouling_disc8"] = []
    plot_data["PSH_temp_wall_shell_fouling_disc9"] = []
    plot_data["PSH_flue_gas_temperature_disc1"] = []
    plot_data["PSH_flue_gas_temperature_disc2"] = []
    plot_data["PSH_flue_gas_temperature_disc4"] = []
    plot_data["PSH_flue_gas_temperature_disc8"] = []
    plot_data["PSH_flue_gas_temperature_disc9"] = []
    plot_data["PSH_von_Mises_inside_disc1"] = []
    plot_data["PSH_von_Mises_inside_disc2"] = []
    plot_data["PSH_von_Mises_inside_disc4"] = []
    plot_data["PSH_von_Mises_inside_disc8"] = []
    plot_data["PSH_von_Mises_inside_disc9"] = []
    plot_data["PSH_delta_s12_inside_disc1"] = []
    plot_data["PSH_delta_s23_inside_disc1"] = []
    plot_data["PSH_delta_s31_inside_disc1"] = []
    plot_data["PSH_mech_circumferential_inside_disc1"] = []
    plot_data["PSH_ther_circumferential_inside_disc1"] = []
    plot_data["PSH_circumferential_inside_disc1"] = []
    plot_data["PSH_mech_circumferential_outside_disc1"] = []
    plot_data["PSH_ther_circumferential_outside_disc1"] = []
    plot_data["PSH_circumferential_outside_disc1"] = []

    # For primary superheater header
    plot_data["Header_inner_temperature"] = []
    plot_data["Header_outside_temperature"] = []
    plot_data["Header_circumferential_P1"] = []
    plot_data["Header_circumferential_P2"] = []
    plot_data["Header_circumferential_body"] = []
    plot_data["Header_von_Mises_P1"] = []
    plot_data["Header_von_Mises_P2"] = []
    plot_data["Header_rupture_time_P1"] = []
    plot_data["Header_rupture_time_P2"] = []

    # Set up and solve a steady-state model first
    m_ss = get_model(dynamic=False)
    # Set up two dynamic models with different time step size
    # Both dynamic models have 2 discretized time steps
    num_step = [2, 2]
    # The time step sizes for the two models are 30 s and 60 s, respectively
    step_size = [30, 60]
    # Set up first dynamic model
    m_dyn0 = get_model(
        dynamic=True, time_set=[0, num_step[0] * step_size[0]], nstep=num_step[0]
    )
    # Set up second dynamic model, if not needed, set it to None
    # m_dyn1 = get_model(dynamic=True,
    #                   time_set=[0, num_step[1]*step_size[1]],
    #                   nstep=num_step[1])
    m_dyn1 = None

    # Copy non-indexed variable values from steady-state to dynamic model
    # Initial condition of dynamic model is steady-state
    if m_dyn0:
        copy_non_time_indexed_values(
            m_dyn0.fs_main, m_ss.fs_main, copy_fixed=True, outlvl=idaeslog.ERROR
        )
    if m_dyn1:
        copy_non_time_indexed_values(
            m_dyn1.fs_main, m_ss.fs_main, copy_fixed=True, outlvl=idaeslog.ERROR
        )

    # model type list with a list of time periods in sequence
    itype_list = []
    # 71 periods of 60 s each (0 to 4260 s)
    for i in range(71):
        itype_list.append(0)
    """
    # Followed by 20 periods of 120 s each (1860 to 4260 s)
    for i in range(20):
        itype_list.append(1)
    # Followed by 30 periods of 60 s each (4260 to 6060 s)
    for i in range(30):
        itype_list.append(0)
    """
    # Total number of period for rolling time windwo simulations
    nperiod = len(itype_list)
    tstart = []
    model_list = []
    t = 0
    for i in range(nperiod):
        tstart.append(t)
        t += step_size[itype_list[i]] * num_step[itype_list[i]]
        if itype_list[i] == 0:
            model_list.append(m_dyn0)
        else:
            model_list.append(m_dyn1)

    # Prepare dynamic model for the first time period
    m_dyn = model_list[0]

    # Copy steady-state results to dynamic model for time-indexed vars
    for t in m_dyn.fs_main.time:
        copy_values_at_time(
            m_dyn.fs_main, m_ss.fs_main, t, 0.0, copy_fixed=True, outlvl=idaeslog.ERROR
        )

    # Reset steady-state bias of controllers to current steady-state values
    # This makes both error and integral error to zero for PID controllers
    # without bounds for manipulated variables
    # For bounded controllers, error is zero and integral error is non-zero
    #
    # Controllers on the steam cycle sub-flowsheet
    t0 = m_dyn.fs_main.time.first()
    m_dyn.fs_main.fs_stc.fwh2_ctrl.mv_ref.value = (
        m_dyn.fs_main.fs_stc.fwh2_valve.valve_opening[t0].value
    )
    m_dyn.fs_main.fs_stc.fwh3_ctrl.mv_ref.value = (
        m_dyn.fs_main.fs_stc.fwh3_valve.valve_opening[t0].value
    )
    m_dyn.fs_main.fs_stc.fwh5_ctrl.mv_ref.value = (
        m_dyn.fs_main.fs_stc.fwh5_valve.valve_opening[t0].value
    )
    m_dyn.fs_main.fs_stc.fwh6_ctrl.mv_ref.value = (
        m_dyn.fs_main.fs_stc.fwh6_valve.valve_opening[t0].value
    )
    m_dyn.fs_main.fs_stc.da_ctrl.mv_ref.value = (
        m_dyn.fs_main.fs_stc.cond_valve.valve_opening[t0].value
    )
    m_dyn.fs_main.fs_stc.makeup_ctrl.mv_ref.value = (
        m_dyn.fs_main.fs_stc.makeup_valve.valve_opening[t0].value
    )
    m_dyn.fs_main.fs_stc.spray_ctrl.mv_ref.value = (
        m_dyn.fs_main.fs_stc.spray_valve.valve_opening[t0].value
    )
    # For two bounded controllers, set the initial integral_of_error
    m_dyn.fs_main.fs_stc.makeup_ctrl.integral_of_error[:].value = pyo.value(
        m_dyn.fs_main.fs_stc.makeup_ctrl.integral_of_error_ref[t0]
    )
    m_dyn.fs_main.fs_stc.spray_ctrl.integral_of_error[:].value = pyo.value(
        m_dyn.fs_main.fs_stc.spray_ctrl.integral_of_error_ref[t0]
    )

    # Controllers on the main flowsheet
    # drum_master_ctrl.mv_ref is always 0 (already fixed to 0)
    m_dyn.fs_main.drum_slave_ctrl.mv_ref.value = (
        m_dyn.fs_main.fs_stc.bfp_turb_valve.valve_opening[t0].value
    )
    m_dyn.fs_main.turbine_master_ctrl.mv_ref.value = (
        m_ss.fs_main.fs_stc.turb.throttle_valve[1].valve_opening[t0].value
    )
    m_dyn.fs_main.turbine_master_ctrl.setpoint[
        :
    ].value = m_ss.fs_main.fs_stc.power_output[t0].value
    m_dyn.fs_main.boiler_master_ctrl.mv_ref.value = (
        m_ss.fs_main.fs_blr.aBoiler.flowrate_coal_raw[t0].value
    )
    m_dyn.fs_main.main_steam_pressure[:].value = (
        m_ss.fs_main.fs_blr.aPlaten.outlet.pressure[t0].value / 1e6
    )
    m_dyn.fs_main.boiler_master_ctrl.setpoint[:].value = (
        m_ss.fs_main.fs_blr.aPlaten.outlet.pressure[t0].value / 1e6
    )

    _log.info(
        "boiler_master_setpoint={}".format(
            m_dyn.fs_main.boiler_master_ctrl.setpoint[0].value
        )
    )
    _log.info(
        "sliding pressure={}".format(pyo.value(m_dyn.fs_main.sliding_pressure[0]))
    )

    # Set solver
    solver = get_solver()

    # Check degree of freedom for the dynamic model
    dof = degrees_of_freedom(m_dyn)
    _log.info("dof of full plant model is {}".format(dof))
    assert dof == 0

    # Solve dynamic model at steady-state
    _log.info("Solving dynamic model at steady-state...")
    solver.solve(m_dyn, tee=True)

    # Start dynamic simulation for 1st period
    _log.info("Solving for time period 0 from {} s".format(tstart[0]))
    # Use steady-state gross power output as basis for load changes
    ss_value = m_ss.fs_main.fs_stc.power_output[0].value
    run_dynamic(m_dyn, ss_value, tstart[0], plot_data, solver)

    # Loop for remaining time periods
    tlast = m_dyn.fs_main.time.last()
    m_prev = m_dyn
    for i in range(1, nperiod):
        m_dyn = model_list[i]
        for t in m_dyn.fs_main.time:
            if itype_list[i] != itype_list[i - 1] or t != tlast:
                # Copy results from previous time period to current period as
                # initial condition and guess
                copy_values_at_time(
                    m_dyn.fs_main,
                    m_prev.fs_main,
                    t,
                    tlast,
                    copy_fixed=True,
                    outlvl=idaeslog.ERROR,
                )
            _log.info(
                "Spray control windup={}".format(
                    m_dyn.fs_main.fs_stc.spray_ctrl.integral_of_error[t].value
                )
            )
            mv_unbounded = pyo.value(m_dyn.fs_main.fs_stc.spray_ctrl.mv_unbounded[t])
            if mv_unbounded < pyo.value(m_dyn.fs_main.fs_stc.spray_ctrl.mv_lb):
                # Reset spray control windup
                if (
                    pyo.value(
                        m_dyn.fs_main.fs_stc.spray_ctrl.integral_of_error[t].value
                    )
                    > 3000
                ):
                    m_dyn.fs_main.fs_stc.spray_ctrl.integral_of_error[t].value = 3000
                    _log.info(
                        "Reached lower bound of manipulated variable and"
                        " maximum windup of 3000. Reset windup to 3000."
                    )
            if mv_unbounded > pyo.value(m_dyn.fs_main.fs_stc.spray_ctrl.mv_ub):
                # Reset spray control windup
                if (
                    pyo.value(
                        m_dyn.fs_main.fs_stc.spray_ctrl.integral_of_error[t].value
                    )
                    < -10000
                ):
                    m_dyn.fs_main.fs_stc.spray_ctrl.integral_of_error[t].value = -10000
                    _log.info(
                        "Reached upper bound of manipulated variable and"
                        " minimum windup of -1e4. Reset windup to -1e4."
                    )
        _log.info("Solving for period number {} from {} s".format(i, tstart[i]))
        # Solve dynamic case for current period of time
        run_dynamic(m_dyn, ss_value, tstart[i], plot_data, solver)
        tlast = m_dyn.fs_main.time.last()
        _log.info(
            "spray flow={}".format(
                m_dyn.fs_main.fs_stc.spray_valve.outlet.flow_mol[tlast].value
            )
        )
        m_prev = m_dyn
    end_time = wall_clock.time()
    time_used = end_time - start_time
    _log.info("simulation time={}".format(time_used))
    write_data_to_txt_file(plot_data)
    plot_results(plot_data)
    return m_dyn


def get_model(dynamic=True, time_set=None, nstep=None, init=True):
    m = pyo.ConcreteModel()
    m.dynamic = dynamic
    if time_set is None:
        time_set = [0, 20, 200]
    if nstep is None:
        nstep = 2
    if m.dynamic:
        m.fs_main = FlowsheetBlock(
            dynamic=True, time_set=time_set, time_units=pyo.units.s
        )
    else:
        m.fs_main = FlowsheetBlock(dynamic=False)

    # Add property packages to flowsheet library
    m.fs_main.prop_water = iapws95.Iapws95ParameterBlock()
    m.fs_main.prop_gas = FlueGasParameterBlock()
    # Declare two sub-flowsheets, one for boiler and the other for steam cycle
    m.fs_main.fs_blr = FlowsheetBlock(time_units=pyo.units.s)
    m.fs_main.fs_stc = FlowsheetBlock(time_units=pyo.units.s)
    # Parameter for sliding pressure slope versus gross power output
    m.fs_main.slope_pslide = pyo.Param(initialize=0.02, doc="slope of sliding pressure")
    m = blr.add_unit_models(m)
    m = stc.add_unit_models(m)
    # Declare a variable for main steam pressure in MPa
    # Used as boiler master process variable
    m.fs_main.main_steam_pressure = pyo.Var(
        m.fs_main.time,
        initialize=13,
        doc="main steam pressure in MPa for boiler master controller",
    )

    # Constraint for main steam pressure in MPa
    @m.fs_main.Constraint(m.fs_main.time, doc="main steam pressure in MPa")
    def main_steam_pressure_eqn(b, t):
        return b.main_steam_pressure[t] == 1e-6 * b.fs_blr.aPlaten.outlet.pressure[t]

    # Add controllers and variables on main flowsheet if dynamic model
    if m.dynamic:
        # master level control output
        # desired feed water flow rate adjustment due to level deviation
        m.fs_main.flow_level_ctrl_output = pyo.Var(
            m.fs_main.time,
            initialize=0,
            doc="feed water flow demand from drum level master controller",
        )
        # PID controllers
        # master of cascading level controller
        m.fs_main.drum_master_ctrl = PIDController(
            process_var=m.fs_main.fs_blr.aDrum.level,
            manipulated_var=m.fs_main.flow_level_ctrl_output,
            type=ControllerType.PI,
            calculate_initial_integral=False,
        )
        # slave of cascading level controller
        m.fs_main.drum_slave_ctrl = PIDController(
            process_var=m.fs_main.fs_stc.bfp.outlet.flow_mol,
            manipulated_var=m.fs_main.fs_stc.bfp_turb_valve.valve_opening,
            type=ControllerType.PI,
            calculate_initial_integral=False,
        )
        # turbine master PID controller to control power output in MW
        # by manipulating throttling valve
        m.fs_main.turbine_master_ctrl = PIDController(
            process_var=m.fs_main.fs_stc.power_output,
            manipulated_var=m.fs_main.fs_stc.turb.throttle_valve[1].valve_opening,
            type=ControllerType.PI,
            calculate_initial_integral=False,
        )
        # boiler master PID controller to control main steam pressure in MPa
        # by manipulating coal feed rate
        m.fs_main.boiler_master_ctrl = PIDController(
            process_var=m.fs_main.main_steam_pressure,
            manipulated_var=m.fs_main.fs_blr.aBoiler.flowrate_coal_raw,
            type=ControllerType.PI,
            calculate_initial_integral=False,
        )

        # Call Pyomo DAE discretizer
        m.discretizer = pyo.TransformationFactory("dae.finite_difference")
        m.discretizer.apply_to(m, nfe=nstep, wrt=m.fs_main.time, scheme="BACKWARD")

        # desired sliding pressure in MPa as a function of power demand in MW
        @m.fs_main.Expression(
            m.fs_main.time,
            doc="Sliding pressure as a function"
            " of power output scheduled for controller",
        )
        def sliding_pressure(b, t):
            return 14.4 + b.slope_pslide * (b.turbine_master_ctrl.setpoint[t] - 320)

        # Constraint for setpoint of the slave controller of
        # three-element drum level controller
        @m.fs_main.Constraint(
            m.fs_main.time, doc="Set point of drum level slave control"
        )
        def drum_level_control_setpoint_eqn(b, t):
            return b.drum_slave_ctrl.setpoint[t] == (
                b.flow_level_ctrl_output[t]
                + b.fs_blr.aPlaten.outlet.flow_mol[t]
                + b.fs_blr.blowdown_split.FW_Blowdown.flow_mol[t]
            )

        # Constraint for setpoint of boiler master
        @m.fs_main.Constraint(m.fs_main.time, doc="Set point of boiler master")
        def boiler_master_setpoint_eqn(b, t):
            return b.boiler_master_ctrl.setpoint[t] == (
                0.02 * (b.turbine_master_ctrl.setpoint[t] - b.fs_stc.power_output[t])
                + b.sliding_pressure[t]
            )

        # Constraint for setpoint of boiler master
        @m.fs_main.Constraint(m.fs_main.time, doc="dry O2 in flue gas in dynamic mode")
        def dry_o2_in_flue_gas_dyn_eqn(b, t):
            return b.fs_blr.aBoiler.fluegas_o2_pct_dry[t] == (
                0.05
                * (b.fs_stc.spray_ctrl.setpoint[t] - b.fs_stc.temperature_main_steam[t])
                - 0.00075 * b.fs_blr.aBoiler.flowrate_coal_raw[t] ** 3
                + 0.067 * b.fs_blr.aBoiler.flowrate_coal_raw[t] ** 2
                - 2.0 * b.fs_blr.aBoiler.flowrate_coal_raw[t]
                + 22.95
            )

        # controller inputs
        # for master level controller
        m.fs_main.drum_master_ctrl.gain_p.fix(40000)
        m.fs_main.drum_master_ctrl.gain_i.fix(100)
        m.fs_main.drum_master_ctrl.setpoint.fix(0.9)
        m.fs_main.drum_master_ctrl.mv_ref.fix(0)
        # for slave level controller
        # note the setpoint is defined by the constraint
        m.fs_main.drum_slave_ctrl.gain_p.fix(2e-2)
        m.fs_main.drum_slave_ctrl.gain_i.fix(2e-4)
        m.fs_main.drum_slave_ctrl.mv_ref.fix(0.5)
        # for turbine master controller, note the setpoint is the power demand
        m.fs_main.turbine_master_ctrl.gain_p.fix(5e-4)
        m.fs_main.turbine_master_ctrl.gain_i.fix(5e-4)
        m.fs_main.turbine_master_ctrl.mv_ref.fix(0.6)
        # for boiler master controller
        # note the setpoint is specified by the constraint
        m.fs_main.boiler_master_ctrl.gain_p.fix(10)
        m.fs_main.boiler_master_ctrl.gain_i.fix(0.25)
        m.fs_main.boiler_master_ctrl.mv_ref.fix(29.0)

        t0 = m.fs_main.time.first()
        m.fs_main.drum_master_ctrl.integral_of_error[t0].fix(0)
        m.fs_main.drum_slave_ctrl.integral_of_error[t0].fix(0)
        m.fs_main.turbine_master_ctrl.integral_of_error[t0].fix(0)
        m.fs_main.boiler_master_ctrl.integral_of_error[t0].fix(0)
        m.fs_main.flow_level_ctrl_output.value = 0
    else:

        @m.fs_main.Constraint(m.fs_main.time, doc="sliding pressure in MPa")
        def sliding_pressure_eqn(b, t):
            return b.main_steam_pressure[t] == 14.4 + b.slope_pslide * (
                b.fs_stc.power_output[t] - 320
            )

    # Set arc connections and inputs for boiler system
    blr.set_arcs_and_constraints(m)
    blr.set_inputs(m)
    # Set arc connections and input for steam cycle system
    stc.set_arcs_and_constraints(m)
    stc.set_inputs(m)
    # Set and calculate scaling factors for the full plant model
    set_scaling_factors(m)
    # Add overall performation expressions
    add_overall_performance_expressions(m)
    if init:
        # Initialize boiler and steam cycle sub-flowsheets
        blr.initialize(m)
        stc.initialize(m)

    # Set arc connections between two sub-flowsheets,
    # deactivate some constraints of two sub-flowsheets
    # and unfix connected streams
    m.fs_main.S001 = Arc(
        source=m.fs_main.fs_blr.aPlaten.outlet,
        destination=m.fs_main.fs_stc.turb.inlet_split.inlet,
    )
    m.fs_main.S005 = Arc(
        source=m.fs_main.fs_stc.turb.hp_split[14].outlet_1,
        destination=m.fs_main.fs_blr.aRH1.tube_inlet,
    )
    m.fs_main.S009 = Arc(
        source=m.fs_main.fs_blr.aRH2.tube_outlet,
        destination=m.fs_main.fs_stc.turb.ip_stages[1].inlet,
    )
    m.fs_main.S042 = Arc(
        source=m.fs_main.fs_stc.fwh6.desuperheat.cold_side_outlet,
        destination=m.fs_main.fs_blr.aECON.tube_inlet,
    )
    m.fs_main.B006 = Arc(
        source=m.fs_main.fs_stc.spray_valve.outlet,
        destination=m.fs_main.fs_blr.Attemp.Water_inlet,
    )
    # Apply Pyomo arc connections
    pyo.TransformationFactory("network.expand_arcs").apply_to(m.fs_main)

    # Deactivate constraints on two sub-flowsheets
    # Boler reheat stream flow based on steam cycle steam extraction
    m.fs_main.fs_blr.flow_mol_steam_rh_eqn.deactivate()
    m.fs_main.fs_stc.turb.constraint_reheat_flow.deactivate()
    # Feed water flow coupled with blowdown and makeup
    m.fs_main.fs_stc.fw_flow_constraint.deactivate()

    # Unfix all connected streams
    m.fs_main.fs_stc.turb.inlet_split.inlet.unfix()
    m.fs_main.fs_stc.turb.hp_split[14].outlet_1.unfix()
    m.fs_main.fs_blr.aRH1.tube_inlet.unfix()
    m.fs_main.fs_stc.turb.ip_stages[1].inlet.unfix()
    m.fs_main.fs_blr.aECON.tube_inlet.unfix()
    m.fs_main.fs_blr.Attemp.Water_inlet.unfix()
    m.fs_main.fs_stc.spray_valve.outlet.unfix()
    m.fs_main.fs_blr.aBoiler.flowrate_coal_raw.unfix()
    # Note that pressure driven model will calculate the steam flow
    # given the feed water pump outlet pressure and the
    # required coal flow rate is determined by steam flow rate

    if m.dynamic is False:
        m.fs_main.fs_stc.turb.throttle_valve[1].valve_opening.unfix()
        _log.info("Solve connected models...")
        _log.info("Degrees of freedom = {}".format(degrees_of_freedom(m)))

        # Solver for solving combined full plant model
        solver = get_solver(options={"max_iter": 50})

        if init:
            res = solver.solve(m, tee=True)
            _log.info("Solved: {}".format(idaeslog.condition(res)))
        # Main performance data
        _log.info(
            "Power output of main turbine={}".format(
                pyo.value(m.fs_main.fs_stc.power_output[0])
            )
        )
        _log.info(
            "Power output of bfp turbine={}".format(
                pyo.value(
                    -m.fs_main.fs_stc.bfp_turb.control_volume.work[0]
                    - m.fs_main.fs_stc.bfp_turb_os.control_volume.work[0]
                )
                * 1e-6
            )
        )
        _log.info(
            "coal flow rate={}".format(
                pyo.value(m.fs_main.fs_blr.aBoiler.flowrate_coal_raw[0])
            )
        )
        _log.info(
            "feed water flow_mol={}".format(
                pyo.value(
                    m.fs_main.fs_stc.fwh6.desuperheat.cold_side_outlet.flow_mol[0]
                )
            )
        )
        _log.info(
            "main steam flow_mol={}".format(
                pyo.value(m.fs_main.fs_stc.turb.inlet_split.inlet.flow_mol[0])
            )
        )
        _log.info(
            "makeup flow={}".format(
                m.fs_main.fs_stc.makeup_valve.outlet.flow_mol[0].value
            )
        )
        _log.info(
            "reheat steam flow_mol={}".format(
                pyo.value(m.fs_main.fs_stc.turb.ip_stages[1].inlet.flow_mol[0])
            )
        )
        _log.info(
            "spray flow={}".format(
                m.fs_main.fs_stc.spray_valve.outlet.flow_mol[0].value
            )
        )
        _log.info(
            "bfpt steam flow={}".format(
                pyo.value(m.fs_main.fs_stc.bfp_turb.inlet.flow_mol[0])
            )
        )
        _log.info(
            "condenser condensate flow_mol={}".format(
                pyo.value(m.fs_main.fs_stc.turb.outlet_stage.outlet.flow_mol[0])
            )
        )
        _log.info(
            "SA temperature={}".format(
                pyo.value(m.fs_main.fs_blr.aBoiler.secondary_air_inlet.temperature[0])
            )
        )
        _log.info(
            "FEGT={}".format(
                pyo.value(m.fs_main.fs_blr.aBoiler.flue_gas_outlet.temperature[0])
            )
        )
        _log.info(
            "main steam temp={}".format(
                pyo.value(
                    m.fs_main.fs_blr.aPlaten.control_volume.properties_out[
                        0
                    ].temperature
                )
            )
        )
        _log.info(
            "reheat steam temp={}".format(
                pyo.value(m.fs_main.fs_blr.aRH2.tube.properties[0, 0].temperature)
            )
        )
        _log.info(
            "ECON inlet temperature={}".format(
                pyo.value(m.fs_main.fs_blr.aECON.tube.properties[0, 1].temperature)
            )
        )
        _log.info(
            "ECON inlet temperature sat={}".format(
                pyo.value(m.fs_main.fs_blr.aECON.tube.properties[0, 1].temperature_sat)
            )
        )
        _log.info(
            "ECON outlet temperature={}".format(
                pyo.value(m.fs_main.fs_blr.aECON.tube.properties[0, 0].temperature)
            )
        )
        _log.info(
            "ECON outlet temperature sat={}".format(
                pyo.value(m.fs_main.fs_blr.aECON.tube.properties[0, 0].temperature_sat)
            )
        )
        _log.info(
            "condenser temp={}".format(
                pyo.value(
                    m.fs_main.fs_stc.turb.outlet_stage.control_volume.properties_out[
                        0
                    ].temperature
                )
            )
        )
        _log.info(
            "FW pressure={}".format(pyo.value(m.fs_main.fs_stc.bfp.outlet.pressure[0]))
        )
        _log.info(
            "ECON inlet pressure={}".format(
                pyo.value(m.fs_main.fs_blr.aECON.tube.properties[0, 1].pressure)
            )
        )
        _log.info(
            "drum pressure={}".format(
                pyo.value(m.fs_main.fs_blr.aRoof.inlet.pressure[0])
            )
        )
        _log.info(
            "main steam pressure={}".format(
                pyo.value(m.fs_main.fs_blr.aPlaten.outlet.pressure[0])
            )
        )
        _log.info(
            "reheat steam pressure={}".format(
                pyo.value(m.fs_main.fs_stc.turb.ip_stages[1].inlet.pressure[0])
            )
        )
        _log.info(
            "main condenser pressure={}".format(
                pyo.value(m.fs_main.fs_stc.turb.outlet_stage.outlet.pressure[0])
            )
        )
        _log.info(
            "Attemp water pressure={}".format(
                pyo.value(m.fs_main.fs_blr.Attemp.Water_inlet.pressure[0])
            )
        )
        _log.info("Cv fwh2={}".format(m.fs_main.fs_stc.fwh2_valve.Cv.value))
        _log.info("Cv fwh3={}".format(m.fs_main.fs_stc.fwh3_valve.Cv.value))
        _log.info("Cv fwh5={}".format(m.fs_main.fs_stc.fwh5_valve.Cv.value))
        _log.info("Cv fwh6={}".format(m.fs_main.fs_stc.fwh6_valve.Cv.value))
        _log.info("Cv cond_valve={}".format(m.fs_main.fs_stc.cond_valve.Cv.value))
        _log.info("Cv makeup_valve={}".format(m.fs_main.fs_stc.makeup_valve.Cv.value))
        _log.info("Cv spray_valve={}".format(m.fs_main.fs_stc.spray_valve.Cv.value))
        _log.info(
            "Cv bfp_turb_valve={}".format(m.fs_main.fs_stc.bfp_turb_valve.Cv.value)
        )
        _log.info(
            "Cv throttle valve={}".format(
                m.fs_main.fs_stc.turb.throttle_valve[1].Cv.value
            )
        )
        _log.info(
            "valve opening fwh2={}".format(
                m.fs_main.fs_stc.fwh2_valve.valve_opening[0].value
            )
        )
        _log.info(
            "valve opening fwh3={}".format(
                m.fs_main.fs_stc.fwh3_valve.valve_opening[0].value
            )
        )
        _log.info(
            "valve opening fwh5={}".format(
                m.fs_main.fs_stc.fwh5_valve.valve_opening[0].value
            )
        )
        _log.info(
            "valve opening fwh6={}".format(
                m.fs_main.fs_stc.fwh6_valve.valve_opening[0].value
            )
        )
        _log.info(
            "valve opening cond_valve={}".format(
                m.fs_main.fs_stc.cond_valve.valve_opening[0].value
            )
        )
        _log.info(
            "valve opening makeup_valve={}".format(
                m.fs_main.fs_stc.makeup_valve.valve_opening[0].value
            )
        )
        _log.info(
            "valve opening spray_valve={}".format(
                m.fs_main.fs_stc.spray_valve.valve_opening[0].value
            )
        )
        _log.info(
            "valve opening bfp_turb_valve={}".format(
                m.fs_main.fs_stc.bfp_turb_valve.valve_opening[0].value
            )
        )
        _log.info(
            "valve opening throttle={}".format(
                m.fs_main.fs_stc.turb.throttle_valve[1].valve_opening[0].value
            )
        )
        _log.info("fwh2 level={}".format(m.fs_main.fs_stc.fwh2.condense.level[0].value))
        _log.info("fwh3 level={}".format(m.fs_main.fs_stc.fwh3.condense.level[0].value))
        _log.info("fwh5 level={}".format(m.fs_main.fs_stc.fwh5.condense.level[0].value))
        _log.info("fwh6 level={}".format(m.fs_main.fs_stc.fwh6.condense.level[0].value))
        _log.info(
            "hotwell tank level={}".format(
                m.fs_main.fs_stc.hotwell_tank.tank_level[0].value
            )
        )
        _log.info(
            "da tank level={}".format(m.fs_main.fs_stc.da_tank.tank_level[0].value)
        )
        _log.info("gross heat rate={}".format(pyo.value(m.fs_main.gross_heat_rate[0])))

        # set load to 320 MW
        _log.info("Set power output is 320 MW...")
        m.fs_main.fs_stc.bfp.outlet.pressure.unfix()
        m.fs_main.fs_blr.aBoiler.flowrate_coal_raw.unfix()
        m.fs_main.fs_stc.power_output.fix(320)
        m.fs_main.fs_stc.spray_valve.valve_opening.unfix()
        m.fs_main.fs_stc.temperature_main_steam.fix(810)
        if init:
            res = solver.solve(m, tee=True)

        # Main performance data
        _log.info(
            "Power output of main turbine={}".format(
                pyo.value(m.fs_main.fs_stc.power_output[0])
            )
        )
        _log.info(
            "Power output of bfp turbine={}".format(
                pyo.value(
                    -m.fs_main.fs_stc.bfp_turb.control_volume.work[0]
                    - m.fs_main.fs_stc.bfp_turb_os.control_volume.work[0]
                )
                * 1e-6
            )
        )
        _log.info(
            "coal flow rate={}".format(
                pyo.value(m.fs_main.fs_blr.aBoiler.flowrate_coal_raw[0])
            )
        )
        _log.info(
            "feed water flow_mol={}".format(
                pyo.value(
                    m.fs_main.fs_stc.fwh6.desuperheat.cold_side_outlet.flow_mol[0]
                )
            )
        )
        _log.info(
            "main steam flow_mol={}".format(
                pyo.value(m.fs_main.fs_stc.turb.inlet_split.inlet.flow_mol[0])
            )
        )
        _log.info(
            "makeup flow={}".format(
                m.fs_main.fs_stc.makeup_valve.outlet.flow_mol[0].value
            )
        )
        _log.info(
            "reheat steam flow_mol={}".format(
                pyo.value(m.fs_main.fs_stc.turb.ip_stages[1].inlet.flow_mol[0])
            )
        )
        _log.info(
            "spray flow={}".format(
                m.fs_main.fs_stc.spray_valve.outlet.flow_mol[0].value
            )
        )
        _log.info(
            "bfpt steam flow={}".format(
                pyo.value(m.fs_main.fs_stc.bfp_turb.inlet.flow_mol[0])
            )
        )
        _log.info(
            "condenser condensate flow_mol={}".format(
                pyo.value(m.fs_main.fs_stc.turb.outlet_stage.outlet.flow_mol[0])
            )
        )
        _log.info(
            "SA temperature={}".format(
                pyo.value(m.fs_main.fs_blr.aBoiler.secondary_air_inlet.temperature[0])
            )
        )
        _log.info(
            "FEGT={}".format(
                pyo.value(m.fs_main.fs_blr.aBoiler.flue_gas_outlet.temperature[0])
            )
        )
        _log.info(
            "main steam temp={}".format(
                pyo.value(
                    m.fs_main.fs_blr.aPlaten.control_volume.properties_out[
                        0
                    ].temperature
                )
            )
        )
        _log.info(
            "reheat steam temp={}".format(
                pyo.value(m.fs_main.fs_blr.aRH2.tube.properties[0, 0].temperature)
            )
        )
        _log.info(
            "ECON inlet temperature={}".format(
                pyo.value(m.fs_main.fs_blr.aECON.tube.properties[0, 1].temperature)
            )
        )
        _log.info(
            "ECON inlet temperature sat={}".format(
                pyo.value(m.fs_main.fs_blr.aECON.tube.properties[0, 1].temperature_sat)
            )
        )
        _log.info(
            "ECON outlet temperature={}".format(
                pyo.value(m.fs_main.fs_blr.aECON.tube.properties[0, 0].temperature)
            )
        )
        _log.info(
            "ECON outlet temperature sat={}".format(
                pyo.value(m.fs_main.fs_blr.aECON.tube.properties[0, 0].temperature_sat)
            )
        )
        _log.info(
            "condenser temp={}".format(
                pyo.value(
                    m.fs_main.fs_stc.turb.outlet_stage.control_volume.properties_out[
                        0
                    ].temperature
                )
            )
        )
        _log.info(
            "FW pressure={}".format(pyo.value(m.fs_main.fs_stc.bfp.outlet.pressure[0]))
        )
        _log.info(
            "ECON inlet pressure={}".format(
                pyo.value(m.fs_main.fs_blr.aECON.tube.properties[0, 1].pressure)
            )
        )
        _log.info(
            "drum pressure={}".format(
                pyo.value(m.fs_main.fs_blr.aRoof.inlet.pressure[0])
            )
        )
        _log.info(
            "main steam pressure={}".format(
                pyo.value(m.fs_main.fs_blr.aPlaten.outlet.pressure[0])
            )
        )
        _log.info(
            "reheat steam pressure={}".format(
                pyo.value(m.fs_main.fs_stc.turb.ip_stages[1].inlet.pressure[0])
            )
        )
        _log.info(
            "main condenser pressure={}".format(
                pyo.value(m.fs_main.fs_stc.turb.outlet_stage.outlet.pressure[0])
            )
        )
        _log.info(
            "Attemp water pressure={}".format(
                pyo.value(m.fs_main.fs_blr.Attemp.Water_inlet.pressure[0])
            )
        )
        _log.info("Cv fwh2={}".format(m.fs_main.fs_stc.fwh2_valve.Cv.value))
        _log.info("Cv fwh3={}".format(m.fs_main.fs_stc.fwh3_valve.Cv.value))
        _log.info("Cv fwh5={}".format(m.fs_main.fs_stc.fwh5_valve.Cv.value))
        _log.info("Cv fwh6={}".format(m.fs_main.fs_stc.fwh6_valve.Cv.value))
        _log.info("Cv cond_valve={}".format(m.fs_main.fs_stc.cond_valve.Cv.value))
        _log.info("Cv makeup_valve={}".format(m.fs_main.fs_stc.makeup_valve.Cv.value))
        _log.info("Cv spray_valve={}".format(m.fs_main.fs_stc.spray_valve.Cv.value))
        _log.info(
            "Cv bfp_turb_valve={}".format(m.fs_main.fs_stc.bfp_turb_valve.Cv.value)
        )
        _log.info(
            "Cv throttle valve={}".format(
                m.fs_main.fs_stc.turb.throttle_valve[1].Cv.value
            )
        )
        _log.info(
            "valve opening fwh2={}".format(
                m.fs_main.fs_stc.fwh2_valve.valve_opening[0].value
            )
        )
        _log.info(
            "valve opening fwh3={}".format(
                m.fs_main.fs_stc.fwh3_valve.valve_opening[0].value
            )
        )
        _log.info(
            "valve opening fwh5={}".format(
                m.fs_main.fs_stc.fwh5_valve.valve_opening[0].value
            )
        )
        _log.info(
            "valve opening fwh6={}".format(
                m.fs_main.fs_stc.fwh6_valve.valve_opening[0].value
            )
        )
        _log.info(
            "valve opening cond_valve={}".format(
                m.fs_main.fs_stc.cond_valve.valve_opening[0].value
            )
        )
        _log.info(
            "valve opening makeup_valve={}".format(
                m.fs_main.fs_stc.makeup_valve.valve_opening[0].value
            )
        )
        _log.info(
            "valve opening spray_valve={}".format(
                m.fs_main.fs_stc.spray_valve.valve_opening[0].value
            )
        )
        _log.info(
            "valve opening bfp_turb_valve={}".format(
                m.fs_main.fs_stc.bfp_turb_valve.valve_opening[0].value
            )
        )
        _log.info(
            "valve opening throttle={}".format(
                m.fs_main.fs_stc.turb.throttle_valve[1].valve_opening[0].value
            )
        )
        _log.info("fwh2 level={}".format(m.fs_main.fs_stc.fwh2.condense.level[0].value))
        _log.info("fwh3 level={}".format(m.fs_main.fs_stc.fwh3.condense.level[0].value))
        _log.info("fwh5 level={}".format(m.fs_main.fs_stc.fwh5.condense.level[0].value))
        _log.info("fwh6 level={}".format(m.fs_main.fs_stc.fwh6.condense.level[0].value))
        _log.info(
            "hotwell tank level={}".format(
                m.fs_main.fs_stc.hotwell_tank.tank_level[0].value
            )
        )
        _log.info(
            "da tank level={}".format(m.fs_main.fs_stc.da_tank.tank_level[0].value)
        )
        _log.info("gross heat rate={}".format(pyo.value(m.fs_main.gross_heat_rate[0])))
    else:
        m.fs_main.fs_blr.dry_o2_in_flue_gas_eqn.deactivate()
        t0 = m.fs_main.time.first()
        m.fs_main.fs_stc.fwh2.condense.level[t0].fix()
        m.fs_main.fs_stc.fwh3.condense.level[t0].fix()
        m.fs_main.fs_stc.fwh5.condense.level[t0].fix()
        m.fs_main.fs_stc.fwh6.condense.level[t0].fix()
        m.fs_main.fs_stc.hotwell_tank.tank_level[t0].fix()
        m.fs_main.fs_stc.da_tank.tank_level[t0].fix()
        m.fs_main.fs_stc.temperature_main_steam[t0].unfix()
        m.fs_main.fs_stc.spray_valve.valve_opening[t0].fix()
        m.fs_main.fs_blr.aDrum.level.unfix()
        m.fs_main.fs_blr.aDrum.level[t0].fix()
        m.fs_main.flow_level_ctrl_output.unfix()
        m.fs_main.flow_level_ctrl_output[t0].fix()
        m.fs_main.fs_stc.bfp.outlet.pressure.unfix()
        m.fs_main.fs_blr.aBoiler.flowrate_coal_raw.unfix()
        m.fs_main.fs_blr.aBoiler.flowrate_coal_raw[t0].fix()
        m.fs_main.fs_stc.turb.throttle_valve[1].valve_opening.unfix()
        m.fs_main.fs_stc.turb.throttle_valve[1].valve_opening[t0].fix()
        m.fs_main.turbine_master_ctrl.setpoint.fix()

    return m


def run_dynamic(m, x0, t0, pd, solver):
    """Function to solve dynamic model for a given period of time"""
    # Set input profile of load demand
    i = 0
    for t in m.fs_main.time:
        power_demand = input_profile(t0 + t, x0)
        m.fs_main.turbine_master_ctrl.setpoint[t].value = power_demand
        _log.info("Point {}, time {}, power demand={}".format(i, t + t0, power_demand))
        i += 1
    # Solve the model
    results = solver.solve(m, tee=True)
    # Print results
    print(results)

    # append results
    for t in m.fs_main.time:
        if t == 0 and len(pd["time"]) > 0:
            continue
        pd["time"].append(t + t0)
        pd["coal_flow"].append(m.fs_main.fs_blr.aBoiler.flowrate_coal_raw[t].value)
        pd["PA_to_APH_flow"].append(
            pyo.value(m.fs_main.fs_blr.aAPH.side_2.properties_in[t].flow_mol)
        )
        pd["PA_temp_flow"].append(
            pyo.value(m.fs_main.fs_blr.Mixer_PA.TA_inlet_state[t].flow_mol)
        )
        pd["SA_flow"].append(
            pyo.value(m.fs_main.fs_blr.aAPH.side_3.properties_in[t].flow_mol)
        )
        pd["SR"].append(m.fs_main.fs_blr.aBoiler.SR[t].value)
        pd["dry_O2_pct_in_flue_gas"].append(
            m.fs_main.fs_blr.aBoiler.fluegas_o2_pct_dry[t].value
        )
        pd["bfpt_opening"].append(
            m.fs_main.fs_stc.bfp_turb_valve.valve_opening[t].value
        )
        pd["gross_power"].append(pyo.value(m.fs_main.fs_stc.power_output[t]))
        pd["ww_heat"].append(pyo.value(m.fs_main.fs_blr.aBoiler.heat_total_ww[t]) / 1e6)
        pd["fegt"].append(m.fs_main.fs_blr.aBoiler.flue_gas_outlet.temperature[t].value)
        pd["drum_level"].append(m.fs_main.fs_blr.aDrum.level[t].value)
        pd["feed_water_flow_sp"].append(m.fs_main.drum_slave_ctrl.setpoint[t].value)
        pd["drum_master_ctrl_op"].append(m.fs_main.flow_level_ctrl_output[t].value)
        pd["feed_water_flow"].append(
            m.fs_main.fs_blr.aECON.tube_inlet.flow_mol[t].value
        )
        pd["main_steam_flow"].append(m.fs_main.fs_blr.aPlaten.outlet.flow_mol[t].value)
        pd["rh_steam_flow"].append(m.fs_main.fs_blr.aRH2.tube_outlet.flow_mol[t].value)
        pd["bfpt_flow"].append(m.fs_main.fs_stc.bfp_turb_valve.outlet.flow_mol[t].value)
        pd["spray_flow"].append(m.fs_main.fs_stc.spray_valve.outlet.flow_mol[t].value)
        pd["main_steam_temp"].append(m.fs_main.fs_stc.temperature_main_steam[t].value)
        pd["rh_steam_temp"].append(
            pyo.value(m.fs_main.fs_blr.aRH2.tube.properties[t, 0].temperature)
        )
        pd["fw_pres"].append(m.fs_main.fs_stc.bfp.outlet.pressure[t].value / 1e6)
        pd["drum_pres"].append(
            m.fs_main.fs_blr.aDrum.water_steam_inlet.pressure[t].value / 1e6
        )
        pd["main_steam_pres"].append(m.fs_main.main_steam_pressure[t].value)
        pd["rh_steam_pres"].append(
            m.fs_main.fs_blr.aRH2.tube_outlet.pressure[t].value / 1e6
        )
        pd["hw_tank_level"].append(m.fs_main.fs_stc.hotwell_tank.tank_level[t].value)
        pd["da_tank_level"].append(m.fs_main.fs_stc.da_tank.tank_level[t].value)
        pd["fwh2_level"].append(m.fs_main.fs_stc.fwh2.condense.level[t].value)
        pd["fwh3_level"].append(m.fs_main.fs_stc.fwh3.condense.level[t].value)
        pd["fwh5_level"].append(m.fs_main.fs_stc.fwh5.condense.level[t].value)
        pd["fwh6_level"].append(m.fs_main.fs_stc.fwh6.condense.level[t].value)
        pd["makeup_valve_opening"].append(
            m.fs_main.fs_stc.makeup_valve.valve_opening[t].value
        )
        pd["cond_valve_opening"].append(
            m.fs_main.fs_stc.cond_valve.valve_opening[t].value
        )
        pd["fwh2_valve_opening"].append(
            m.fs_main.fs_stc.fwh2_valve.valve_opening[t].value
        )
        pd["fwh3_valve_opening"].append(
            m.fs_main.fs_stc.fwh3_valve.valve_opening[t].value
        )
        pd["fwh5_valve_opening"].append(
            m.fs_main.fs_stc.fwh5_valve.valve_opening[t].value
        )
        pd["fwh6_valve_opening"].append(
            m.fs_main.fs_stc.fwh6_valve.valve_opening[t].value
        )
        pd["spray_valve_opening"].append(
            m.fs_main.fs_stc.spray_valve.valve_opening[t].value
        )
        pd["tube_temp_rh2"].append(
            m.fs_main.fs_blr.aRH2.shell_wall_temperature[t, 0].value
        )
        pd["temp_fg_econ_exit"].append(
            m.fs_main.fs_blr.aECON.shell_outlet.temperature[t].value
        )
        pd["temp_fg_aph_exit"].append(
            m.fs_main.fs_blr.aAPH.side_1_outlet.temperature[t].value
        )
        pd["throttle_opening"].append(
            m.fs_main.fs_stc.turb.throttle_valve[1].valve_opening[t].value
        )
        pd["load_demand"].append(m.fs_main.turbine_master_ctrl.setpoint[t].value)
        pd["sliding_pressure"].append(pyo.value(m.fs_main.sliding_pressure[t]))
        pd["steam_pressure_sp"].append(m.fs_main.boiler_master_ctrl.setpoint[t].value)
        pd["deaerator_pressure"].append(
            m.fs_main.fs_stc.da_tank.inlet.pressure[t].value / 1e6
        )
        pd["temp_econ_in"].append(
            pyo.value(m.fs_main.fs_blr.aECON.tube.properties[t, 1].temperature)
        )
        pd["temp_econ_out"].append(
            pyo.value(m.fs_main.fs_blr.aECON.tube.properties[t, 0].temperature)
        )
        pd["temp_econ_out_sat"].append(
            pyo.value(m.fs_main.fs_blr.aECON.tube.properties[t, 0].temperature_sat)
        )
        pd["boiler_efficiency_heat"].append(
            pyo.value(m.fs_main.fs_blr.boiler_efficiency_heat[t])
        )
        pd["boiler_efficiency_steam"].append(
            pyo.value(m.fs_main.fs_blr.boiler_efficiency_steam[t])
        )
        pd["steam_cycle_efficiency"].append(pyo.value(m.fs_main.steam_cycle_eff[t]))
        pd["plant_gross_efficiency"].append(
            pyo.value(m.fs_main.plant_gross_efficiency[t])
        )
        pd["gross_heat_rate"].append(pyo.value(m.fs_main.gross_heat_rate[t]))

        # Boiler health related data
        # For drum
        pd["drum_inner_wall_temperature"].append(
            pyo.value(
                m.fs_main.fs_blr.aDrum.drum_wall_temperature[
                    t, m.fs_main.fs_blr.aDrum.radial_domain.first()
                ]
            )
        )
        pd["drum_outer_wall_temperature"].append(
            pyo.value(
                m.fs_main.fs_blr.aDrum.drum_wall_temperature[
                    t, m.fs_main.fs_blr.aDrum.radial_domain.last()
                ]
            )
        )
        pd["drum_inner_theta_sigma_m"].append(
            pyo.value(
                m.fs_main.fs_blr.aDrum.mech_sigma_theta[
                    t, m.fs_main.fs_blr.aDrum.radial_domain.first()
                ]
            )
        )
        pd["drum_inner_theta_sigma_t"].append(
            pyo.value(
                m.fs_main.fs_blr.aDrum.therm_sigma_theta[
                    t, m.fs_main.fs_blr.aDrum.radial_domain.first()
                ]
            )
        )
        pd["drum_inner_vM_sigma"].append(
            pyo.value(
                m.fs_main.fs_blr.aDrum.sigma_von_Mises[
                    t, m.fs_main.fs_blr.aDrum.radial_domain.first()
                ]
            )
        )
        pd["drum_outer_theta_sigma_m"].append(
            pyo.value(
                m.fs_main.fs_blr.aDrum.mech_sigma_theta[
                    t, m.fs_main.fs_blr.aDrum.radial_domain.last()
                ]
            )
        )
        pd["drum_outer_theta_sigma_t"].append(
            pyo.value(
                m.fs_main.fs_blr.aDrum.therm_sigma_theta[
                    t, m.fs_main.fs_blr.aDrum.radial_domain.last()
                ]
            )
        )
        pd["drum_outer_vM_sigma"].append(
            pyo.value(
                m.fs_main.fs_blr.aDrum.sigma_von_Mises[
                    t, m.fs_main.fs_blr.aDrum.radial_domain.last()
                ]
            )
        )
        pd["drum_delta_sigma_r_theta_outer"].append(
            pyo.value(
                m.fs_main.fs_blr.aDrum.delta_sigma_r_theta[
                    t, m.fs_main.fs_blr.aDrum.radial_domain.last()
                ]
            )
        )
        pd["drum_delta_sigma_theta_z_outer"].append(
            pyo.value(
                m.fs_main.fs_blr.aDrum.delta_sigma_theta_z[
                    t, m.fs_main.fs_blr.aDrum.radial_domain.last()
                ]
            )
        )
        pd["drum_delta_sigma_z_r_outer"].append(
            pyo.value(
                m.fs_main.fs_blr.aDrum.delta_sigma_z_r[
                    t, m.fs_main.fs_blr.aDrum.radial_domain.last()
                ]
            )
        )
        pd["drum_circumferential_P1"].append(
            pyo.value(m.fs_main.fs_blr.aDrum.sigma_theta_P1[t])
        )
        pd["drum_circumferential_P2"].append(
            pyo.value(m.fs_main.fs_blr.aDrum.sigma_theta_P2[t])
        )
        pd["drum_von_Mises_P1"].append(
            pyo.value(m.fs_main.fs_blr.aDrum.sigma_eff_P1[t])
        )
        pd["drum_von_Mises_P2"].append(
            pyo.value(m.fs_main.fs_blr.aDrum.sigma_eff_P2[t])
        )

        # For primary superheater
        pd["PSH_flue_gas_flow_disc1"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.shell.properties[t, 0].flow_mol)
        )
        pd["PSH_steam_flow_disc9"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.tube.properties[t, 1].flow_mol)
        )
        pd["PSH_steam_temperature_disc1"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.tube.properties[t, 0].temperature)
        )
        pd["PSH_steam_temperature_disc2"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.tube.properties[t, 0.166667].temperature)
        )
        pd["PSH_steam_temperature_disc4"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.tube.properties[t, 0.5].temperature)
        )
        pd["PSH_steam_temperature_disc8"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.tube.properties[t, 0.833333].temperature)
        )
        pd["PSH_steam_temperature_disc9"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.tube.properties[t, 1].temperature)
        )
        pd["PSH_temp_wall_tube_disc1"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.tube_wall_temperature[
                    t, 0, m.fs_main.fs_blr.aPSH.r.first()
                ]
            )
        )
        pd["PSH_temp_wall_tube_disc2"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.tube_wall_temperature[
                    t, 0.166667, m.fs_main.fs_blr.aPSH.r.first()
                ]
            )
        )
        pd["PSH_temp_wall_tube_disc4"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.tube_wall_temperature[
                    t, 0.5, m.fs_main.fs_blr.aPSH.r.first()
                ]
            )
        )
        pd["PSH_temp_wall_tube_disc8"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.tube_wall_temperature[
                    t, 0.833333, m.fs_main.fs_blr.aPSH.r.first()
                ]
            )
        )
        pd["PSH_temp_wall_tube_disc9"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.tube_wall_temperature[
                    t, 1, m.fs_main.fs_blr.aPSH.r.first()
                ]
            )
        )
        pd["PSH_temp_wall_shell_disc1"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.tube_wall_temperature[
                    t, 0, m.fs_main.fs_blr.aPSH.r.last()
                ]
            )
        )
        pd["PSH_temp_wall_shell_disc2"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.tube_wall_temperature[
                    t, 0.166667, m.fs_main.fs_blr.aPSH.r.last()
                ]
            )
        )
        pd["PSH_temp_wall_shell_disc4"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.tube_wall_temperature[
                    t, 0.5, m.fs_main.fs_blr.aPSH.r.last()
                ]
            )
        )
        pd["PSH_temp_wall_shell_disc8"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.tube_wall_temperature[
                    t, 0.833333, m.fs_main.fs_blr.aPSH.r.last()
                ]
            )
        )
        pd["PSH_temp_wall_shell_disc9"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.tube_wall_temperature[
                    t, 1, m.fs_main.fs_blr.aPSH.r.last()
                ]
            )
        )
        pd["PSH_temp_wall_shell_fouling_disc1"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.shell_wall_temperature[t, 0])
        )
        pd["PSH_temp_wall_shell_fouling_disc2"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.shell_wall_temperature[t, 0.166667])
        )
        pd["PSH_temp_wall_shell_fouling_disc4"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.shell_wall_temperature[t, 0.5])
        )
        pd["PSH_temp_wall_shell_fouling_disc8"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.shell_wall_temperature[t, 0.833333])
        )
        pd["PSH_temp_wall_shell_fouling_disc9"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.shell_wall_temperature[t, 1])
        )
        pd["PSH_flue_gas_temperature_disc1"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.shell.properties[t, 0].temperature)
        )
        pd["PSH_flue_gas_temperature_disc2"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.shell.properties[t, 0.166667].temperature)
        )
        pd["PSH_flue_gas_temperature_disc4"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.shell.properties[t, 0.5].temperature)
        )
        pd["PSH_flue_gas_temperature_disc8"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.shell.properties[t, 0.833333].temperature)
        )
        pd["PSH_flue_gas_temperature_disc9"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.shell.properties[t, 1].temperature)
        )
        pd["PSH_von_Mises_inside_disc1"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.sigma_von_Mises[
                    t, 0, m.fs_main.fs_blr.aPSH.r.first()
                ]
            )
        )
        pd["PSH_von_Mises_inside_disc2"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.sigma_von_Mises[
                    t, 0.166667, m.fs_main.fs_blr.aPSH.r.first()
                ]
            )
        )
        pd["PSH_von_Mises_inside_disc4"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.sigma_von_Mises[
                    t, 0.5, m.fs_main.fs_blr.aPSH.r.first()
                ]
            )
        )
        pd["PSH_von_Mises_inside_disc8"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.sigma_von_Mises[
                    t, 0.833333, m.fs_main.fs_blr.aPSH.r.first()
                ]
            )
        )
        pd["PSH_von_Mises_inside_disc9"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.sigma_von_Mises[
                    t, 1, m.fs_main.fs_blr.aPSH.r.first()
                ]
            )
        )
        pd["PSH_delta_s12_inside_disc1"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.delta_sigma_r_theta[
                    t, 0, m.fs_main.fs_blr.aPSH.r.first()
                ]
            )
        )
        pd["PSH_delta_s23_inside_disc1"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.delta_sigma_theta_z[
                    t, 0, m.fs_main.fs_blr.aPSH.r.first()
                ]
            )
        )
        pd["PSH_delta_s31_inside_disc1"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.delta_sigma_z_r[
                    t, 0, m.fs_main.fs_blr.aPSH.r.first()
                ]
            )
        )
        pd["PSH_mech_circumferential_inside_disc1"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.mech_sigma_theta[
                    t, 0, m.fs_main.fs_blr.aPSH.r.first()
                ]
            )
        )
        pd["PSH_ther_circumferential_inside_disc1"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.therm_sigma_theta[
                    t, 0, m.fs_main.fs_blr.aPSH.r.first()
                ]
            )
        )
        pd["PSH_circumferential_inside_disc1"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.sigma_theta[t, 0, m.fs_main.fs_blr.aPSH.r.first()]
            )
        )
        pd["PSH_mech_circumferential_outside_disc1"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.mech_sigma_theta[
                    t, 0, m.fs_main.fs_blr.aPSH.r.last()
                ]
            )
        )
        pd["PSH_ther_circumferential_outside_disc1"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.therm_sigma_theta[
                    t, 0, m.fs_main.fs_blr.aPSH.r.last()
                ]
            )
        )
        pd["PSH_circumferential_outside_disc1"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.sigma_theta[t, 0, m.fs_main.fs_blr.aPSH.r.last()]
            )
        )

        # For primary superheater header
        pd["Header_inner_temperature"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.header_wall_temperature[
                    t, m.fs_main.fs_blr.aPSH.head_r.first()
                ]
            )
        )
        pd["Header_outside_temperature"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.header_wall_temperature[
                    t, m.fs_main.fs_blr.aPSH.head_r.last()
                ]
            )
        )
        pd["Header_circumferential_P1"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.sigma_theta_P1[t])
        )
        pd["Header_circumferential_P2"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.sigma_theta_P2[t])
        )
        pd["Header_circumferential_body"].append(
            pyo.value(
                m.fs_main.fs_blr.aPSH.sigma_theta_header[
                    t, m.fs_main.fs_blr.aPSH.head_r.first()
                ]
            )
        )
        pd["Header_von_Mises_P1"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.sigma_eff_P1[t])
        )
        pd["Header_von_Mises_P2"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.sigma_eff_P2[t])
        )
        pd["Header_rupture_time_P1"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.rupture_time_crotch_corner[t])
        )
        pd["Header_rupture_time_P2"].append(
            pyo.value(m.fs_main.fs_blr.aPSH.rupture_time_P2[t])
        )


def plot_results(pd):
    # ploting responses
    plt.figure(1)
    plt.plot(pd["time"], pd["coal_flow"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Coal Flow Rate [kg/s]")
    plt.show(block=False)

    plt.figure(2)
    plt.plot(pd["time"], pd["bfpt_opening"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("BFPT Valve Opening")
    plt.show(block=False)

    plt.figure(3)
    plt.plot(pd["time"], pd["gross_power"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Gross Power Output [MW]")
    plt.show(block=False)

    plt.figure(4)
    plt.plot(pd["time"], pd["ww_heat"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Waterwall Heat [MW]")
    plt.show(block=False)

    plt.figure(5)
    plt.plot(pd["time"], pd["fegt"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("FEGT [K]")
    plt.show(block=False)

    plt.figure(6)
    plt.plot(pd["time"], pd["drum_level"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Drum Level [m]")
    plt.show(block=False)

    plt.figure(7)
    plt.plot(pd["time"], pd["feed_water_flow"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Feed Water Flow [kmol/s]")
    plt.show(block=False)

    plt.figure(8)
    plt.plot(pd["time"], pd["main_steam_flow"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Main Steam Flow [kmol/s]")
    plt.show(block=False)

    plt.figure(9)
    plt.plot(pd["time"], pd["rh_steam_flow"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("RH Steam Flow [kmol/s]")
    plt.show(block=False)

    plt.figure(10)
    plt.plot(pd["time"], pd["bfpt_flow"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("BFPT Flow [mol/s]")
    plt.show(block=False)

    plt.figure(11)
    plt.plot(pd["time"], pd["spray_flow"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Water Spray Flow [mol/s]")
    plt.show(block=False)

    plt.figure(12)
    plt.plot(pd["time"], pd["main_steam_temp"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Main Steam Temperature [K]")
    plt.show(block=False)

    plt.figure(13)
    plt.plot(pd["time"], pd["rh_steam_temp"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("RH Steam Temperature [K]")
    plt.show(block=False)

    plt.figure(14)
    plt.plot(pd["time"], pd["fw_pres"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Feed Water Pressure [MPa]")
    plt.show(block=False)

    plt.figure(15)
    plt.plot(pd["time"], pd["drum_pres"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Drum Pressure [MPa]")
    plt.show(block=False)

    plt.figure(16)
    plt.plot(pd["time"], pd["main_steam_pres"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Main Steam Pressure [MPa]")
    plt.show(block=False)

    plt.figure(17)
    plt.plot(pd["time"], pd["rh_steam_pres"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("RH Steam Pressure [MPa]")
    plt.show(block=False)

    plt.figure(18)
    plt.plot(pd["time"], pd["hw_tank_level"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Hotwell Tank Level [m]")
    plt.show(block=False)

    plt.figure(19)
    plt.plot(pd["time"], pd["da_tank_level"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("DA Tank Level [m]")
    plt.show(block=False)

    plt.figure(20)
    plt.plot(pd["time"], pd["fwh2_level"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("FWH2 Level [m]")
    plt.show(block=False)

    plt.figure(21)
    plt.plot(pd["time"], pd["fwh3_level"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("FWH3 Level [m]")
    plt.show(block=False)

    plt.figure(22)
    plt.plot(pd["time"], pd["fwh5_level"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("FWH5 Level [m]")
    plt.show(block=False)

    plt.figure(23)
    plt.plot(pd["time"], pd["fwh6_level"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("FWH6 Level [m]")
    plt.show(block=False)

    plt.figure(24)
    plt.plot(pd["time"], pd["makeup_valve_opening"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Makeup Valve Opening")
    plt.show(block=False)

    plt.figure(25)
    plt.plot(pd["time"], pd["cond_valve_opening"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Condensate Valve Opening")
    plt.show(block=False)

    plt.figure(26)
    plt.plot(pd["time"], pd["fwh2_valve_opening"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("FWH2 Valve Opening")
    plt.show(block=False)

    plt.figure(27)
    plt.plot(pd["time"], pd["fwh3_valve_opening"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("FWH3 Valve Opening")
    plt.show(block=False)

    plt.figure(28)
    plt.plot(pd["time"], pd["fwh5_valve_opening"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("FWH5 Valve Opening")
    plt.show(block=False)

    plt.figure(29)
    plt.plot(pd["time"], pd["fwh6_valve_opening"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("FWH6 Valve Opening")
    plt.show(block=False)

    plt.figure(30)
    plt.plot(pd["time"], pd["spray_valve_opening"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Spray Valve Opening")
    plt.show(block=False)

    plt.figure(31)
    plt.plot(pd["time"], pd["tube_temp_rh2"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("RH Tube Temperature [K]")
    plt.show(block=False)

    plt.figure(32)
    plt.plot(pd["time"], pd["temp_fg_econ_exit"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Flue Gas T at Econ Exit [K]")
    plt.show(block=False)

    plt.figure(33)
    plt.plot(pd["time"], pd["temp_fg_aph_exit"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Flue Gas T at APH Exit [K]")
    plt.show(block=False)

    plt.figure(34)
    plt.plot(pd["time"], pd["throttle_opening"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Throttle Valve Opening")
    plt.show(block=False)

    plt.figure(35)
    plt.plot(pd["time"], pd["load_demand"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Load Demand [MW]")
    plt.show(block=False)

    plt.figure(36)
    plt.plot(pd["time"], pd["sliding_pressure"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Desired Sliding Pressure [MPa]")
    plt.show(block=False)

    plt.figure(37)
    plt.plot(pd["time"], pd["gross_heat_rate"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("Gross Heat Rate [BTU/kW-hr]")
    plt.show(block=False)

    plt.figure(38)
    plt.plot(pd["time"], pd["deaerator_pressure"])
    plt.grid()
    plt.xlabel("Time [s]")
    plt.ylabel("DA pressure [MPa]")
    plt.show(block=True)


def write_data_to_txt_file(plot_data):
    # write data to file
    ntime = len(plot_data["time"])
    ncount = len(plot_data)
    icount = 0
    with open("case_5pct_result.txt", "w") as fout:
        # write headings
        for k in plot_data:
            icount += 1
            fout.write(k)
            if icount < ncount:
                fout.write("\t")
            else:
                fout.write("\n")
        # write values
        for i in range(ntime):
            icount = 0
            for k in plot_data:
                icount += 1
                fout.write(str(plot_data[k][i]))
                if icount < ncount:
                    fout.write("\t")
                else:
                    fout.write("\n")


def print_pfd_results(m):
    streams = tables.arcs_to_stream_dict(
        m.fs_main,
        descend_into=False,
        additional={
            "S002": m.fs_main.fs_stc.turb.inlet_stage[1].inlet,
            "S012": m.fs_main.fs_stc.turb.lp_stages[1].inlet,
            "S050b": m.fs_main.fs_stc.condenser_hotwell.makeup,
            "S050": m.fs_main.fs_stc.makeup_valve.inlet,
            "S018": m.fs_main.fs_stc.condenser.tube.properties_in,
            "S020": m.fs_main.fs_stc.condenser.tube.properties_out,
            "S048": m.fs_main.fs_stc.aux_condenser.tube.properties_in,
            "S049": m.fs_main.fs_stc.aux_condenser.tube.properties_out,
            "FG06": m.fs_main.fs_blr.aAPH.side_1_outlet,
            "FG01": m.fs_main.fs_blr.aRH2.shell_inlet,
            "FG04": m.fs_main.fs_blr.aECON.shell_inlet,
            "B001": m.fs_main.fs_blr.aECON.tube_outlet,
        },
    )
    tables.arcs_to_stream_dict(m.fs_main.fs_stc, s=streams, descend_into=False)
    sd = tables.stream_states_dict(streams, 0)
    sdf = tables.generate_table(
        blocks=sd,
        attributes=[
            "flow_mass",
            "flow_mol",
            "temperature",
            "pressure",
            "enth_mol",
            "vapor_frac",
            ("flow_mol_comp", "O2"),
            ("flow_mol_comp", "N2"),
            ("flow_mol_comp", "NO"),
            ("flow_mol_comp", "CO2"),
            ("flow_mol_comp", "H2O"),
            ("flow_mol_comp", "SO2"),
        ],
        exception=False,
    )
    sdf.sort_index(inplace=True)
    sdf.to_csv("streams.csv")
    tags = {}
    tag_formats = {}
    for i, s in sd.items():
        tags[i + "_Fmass"] = s.flow_mass
        tag_formats[i + "_Fmass"] = (
            lambda x: "{:.1f} kg/s" if x >= 10 else "{:.2f} kg/s"
        )
        tags[i + "_F"] = s.flow_mol
        tag_formats[i + "_F"] = "{:,.0f} mol/s"
        tags[i + "_T"] = s.temperature
        tag_formats[i + "_T"] = "{:,.0f} K"
        tags[i + "_P_kPa"] = s.pressure / 1000
        tag_formats[i + "_P_kPa"] = (
            lambda x: "{:,.0f} kPa" if x >= 100 else "{:.2f} kPa"
        )
        tags[i + "_P"] = s.pressure / 1000
        tag_formats[i + "_P"] = "{:,.0f} Pa"
        try:
            tags[i + "_hmass"] = s.enth_mass / 1000.0
            tag_formats[i + "_hmass"] = "{:,.0f} kJ/kg"
            tags[i + "_h"] = s.enth_mol
            tag_formats[i + "_h"] = "{:,.0f} J/mol"
        except AttributeError:
            pass
        try:
            tags[i + "_x"] = s.vapor_frac
            tag_formats[i + "_x"] = "{:.3f}"
        except AttributeError:
            pass
        try:
            tags[i + "_yN2"] = s.mole_frac_comp["N2"]
            tags[i + "_yO2"] = s.mole_frac_comp["O2"]
            tags[i + "_yNO"] = s.mole_frac_comp["NO"]
            tags[i + "_yCO2"] = s.mole_frac_comp["CO2"]
            tags[i + "_yH2O"] = s.mole_frac_comp["H2O"]
            tags[i + "_ySO2"] = s.mole_frac_comp["SO2"]
        except AttributeError:
            pass

    tag_group = ModelTagGroup()
    for t, v in tags.items():
        try:
            formatter = tag_formats[t]
        except KeyError:
            formatter = "{:.3f}"
        tag_group.add(t, v, format_string=formatter)

    dirpath = os.path.dirname(__file__)
    svgpath = os.path.join(dirpath, "plant_pfd.svg")
    with open(svgpath, "r") as f:
        svg_tag(
            svg=f,
            tag_group=tag_group,
            outfile="plant_pfd_result.svg",
        )
    os.remove("streams.csv")


if __name__ == "__main__":
    """
    Main function to to run simulation
    To run steady-state model, call main_steady()
    to run dynamic model, call main_dyn()
    """
    # This method builds and runs a subcritical coal-fired power plant
    # dynamic simulation. The simulation consists of 5%/min ramping down from
    # full load to 50% load, holding for 30 minutes and then ramping up
    # to 100% load and holding for 20 minutes.
    # uncomment the code (line 1821) to run this simulation,
    # note that this simulation takes around ~60 minutes to complete
    # m_dyn = main_dynamic()

    # This method builds and runs a steady state subcritical coal-fired power
    # plant, the simulation consists of a typical base load case.
    m_ss = main_steady_state()
    print_pfd_results(m_ss)
