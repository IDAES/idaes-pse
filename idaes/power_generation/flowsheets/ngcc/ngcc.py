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

import os
import csv
from collections import OrderedDict
import numpy as np

import pyomo.environ as pyo
from pyomo.environ import units as pyunits
from pyomo.network import Arc

import idaes
from idaes.core.solvers import use_idaes_solver_configuration_defaults
from idaes.power_generation.flowsheets.gas_turbine import gas_turbine
from idaes.power_generation.flowsheets.hrsg import hrsg_flowsheet as hrsg_module
from idaes.power_generation.flowsheets.ngcc import steam_turbine as sturb_module
from idaes.power_generation.flowsheets.ngcc.ngcc_costing import (
    get_ngcc_costing,
    build_ngcc_OM_costs
)
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.generic_models.unit_models as gum
import idaes.core.util.misc as imisc
import idaes.core.util.model_serializer as ms
import idaes.core.util.scaling as iscale
from idaes.generic_models.unit_models.heat_exchanger import (
    delta_temperature_underwood_callback
)
import idaes.core.util.tables as tables
from pyomo.environ import units as pyunits

def get_model():
    expand_arcs = pyo.TransformationFactory("network.expand_arcs")
    use_idaes_solver_configuration_defaults()
    idaes.cfg.ipopt["options"]["nlp_scaling_method"] = "user-scaling"
    idaes.cfg.ipopt["options"]["bound_push"] = 1e-10
    idaes.cfg.ipopt["options"]['ma27_pivtol'] = 0.05
    idaes.cfg.ipopt["options"]['ma27_pivtolmax'] = 0.9
    idaes.cfg.ipopt["options"]["max_iter"] = 50
    comps = { # components present
        "CH4", "C2H6", "C3H8", "C4H10", "O2", "H2O", "CO2", "N2", "Ar"}
    rxns = { # reactions and key components for conversion
        "ch4_cmb":"CH4",
        "c2h6_cmb":"C2H6",
        "c3h8_cmb":"C3H8",
        "c4h10_cmb":"C4H10"}
    phases = ["Vap"]
    air_comp = {
        "CH4":0.0,
        "C2H6":0.0,
        "C3H8":0.0,
        "C4H10":0.0,
        "O2":0.2074,
        "H2O":0.0099,
        "CO2":0.0003,
        "N2":0.7732,
        "Ar":0.0092}
    ng_comp = {
        "CH4":0.931,
        "C2H6":0.0320,
        "C3H8":0.007,
        "C4H10":0.004,
        "O2":0.0,
        "H2O":0.0,
        "CO2":0.01,
        "N2":0.0160,
        "Ar":0.0}

    m = pyo.ConcreteModel("NGCC")
    if os.path.isfile("init_ngcc.json.gz"):
        m, solver = gas_turbine.main(
            m=m,
            initialize=False,
            comps=comps,
            rxns=rxns,
            phases=phases,
            air_comp=air_comp,
            ng_comp=ng_comp)
    else:
        m, solver = gas_turbine.main(
            m=m,
            comps=comps,
            rxns=rxns,
            phases=phases,
            air_comp=air_comp,
            ng_comp=ng_comp)
        gas_turbine.run_full_load(m, solver)

    if not os.path.isfile("init_ngcc.json.gz"):
        print("Adjust GT")
        m.fs.gt_power[0].fix(-477e6)
        m.fs.exhaust_1.temperature.fix(898)
        m.fs.feed_fuel1.temperature.fix(460)
        m.fs.gts3.control_volume.properties_out[0].pressure.fix(103421)
        solver.solve(m, tee=True)

    m.fs.net_power_mw = pyo.Var(m.fs.config.time, initialize=600)

    hrsg_module.get_model(m)
    hrsg_module.set_scaling_factors(m)

    if not os.path.isfile("init_ngcc.json.gz"):
        print("Solve Disconnected GT and HRSG (to be sure its all good)")
        solver.solve(m, tee=True)

    m.fs.fg_translate = gum.Translator(default={
        "inlet_property_package": m.fs.gas_prop_params,
        "outlet_property_package": m.fs.prop_gas})
    m.fs.g08t = Arc(
        source=m.fs.exhaust_1.inlet, destination=m.fs.fg_translate.inlet)
    m.fs.g08t_tofg = Arc(
        source=m.fs.fg_translate.outlet, destination=m.fs.HP_SH4.side_2_inlet)
    expand_arcs.apply_to(m)

    m.fs.HP_SH4.side_2_inlet.unfix()
    imisc.copy_port_values(
        source=m.fs.exhaust_1.inlet,
        destination=m.fs.fg_translate.inlet)
    @m.fs.fg_translate.Constraint(m.fs.time, m.fs.prop_gas.component_list)
    def mol_frac_eqn(b, t, i):
        return b.inlet.flow_mol[t]*b.inlet.mole_frac_comp[t, i] == \
            b.outlet.flow_mol_comp[t, i]
    for (t, i) in m.fs.fg_translate.mol_frac_eqn:
        iscale.constraint_scaling_transform(m.fs.fg_translate.mol_frac_eqn[t, i], 1e-2)
    @m.fs.fg_translate.Constraint(m.fs.time)
    def temperature_eqn(b, t):
        return b.inlet.temperature[t] == b.outlet.temperature[t]
    for t in m.fs.fg_translate.temperature_eqn:
        iscale.constraint_scaling_transform(m.fs.fg_translate.temperature_eqn[t], 1e-2)
    @m.fs.fg_translate.Constraint(m.fs.time)
    def pressure_eqn(b, t):
        return b.inlet.pressure[t] == b.outlet.pressure[t]
    for t in m.fs.fg_translate.pressure_eqn:
        iscale.constraint_scaling_transform(m.fs.fg_translate.pressure_eqn[t], 1e-5)

    iscale.calculate_scaling_factors(m)
    if not os.path.isfile("init_ngcc.json.gz"):
        m.fs.fg_translate.initialize()
        solver.solve(m, tee=True)

    if os.path.isfile("init_ngcc.json.gz"):
        ms.from_json(m, fname="init_ngcc.json.gz", wts=ms.StoreSpec(suffix=False))
    else:
        ms.to_json(m, fname="init_ngcc.json.gz")

    #Add steam turbine
    if not os.path.isfile("init_ngcc_st.json.gz"):
        solver.solve(m, tee=True)

    if not os.path.isfile("init_ngcc_st.json.gz"):
        sturb_module.get_model(m)
    else:
        sturb_module.get_model(m, init=False)

    strip_bounds = pyo.TransformationFactory("contrib.strip_var_bounds")
    strip_bounds.apply_to(m, reversible=False)

    print("Connect LP Steam")
    m.fs.t05 = Arc(
        source=m.fs.LP_SH.side_1_outlet,
        destination=m.fs.steam_turbine_lp_mix.hrsg)
    imisc.copy_port_values(arc=m.fs.t05)
    expand_arcs.apply_to(m)
    m.fs.steam_turbine_lp_mix.hrsg.unfix()
    if not os.path.isfile("init_ngcc_st.json.gz"):
        solver.solve(m, tee=True)

    print("connect the reheater")
    m.fs.t02_dummy.deactivate()
    m.fs.t03_dummy.deactivate()
    m.fs.dummy_reheat.deactivate()
    m.fs.t02 = Arc(
        source=m.fs.steam_turbine.hp_stages[7].outlet,
        destination=m.fs.IP_Splitter2.inlet)
    m.fs.t03 = Arc(
        source=m.fs.IP_SH3.side_1_outlet,
        destination=m.fs.steam_turbine.ip_stages[1].inlet)
    expand_arcs.apply_to(m)
    m.fs.IP_Splitter2.inlet.unfix()
    imisc.copy_port_values(arc=m.fs.t02)
    imisc.copy_port_values(arc=m.fs.t03)
    if not os.path.isfile("init_ngcc_st.json.gz"):
        solver.solve(m, tee=True)

    print("Tie throttle valves together and set delta P")
    dp = pyo.value(-5e5)
    m.fs.steam_turbine.throttle_valve[1].deltaP[0].fix(dp)
    m.fs.steam_turbine.throttle_valve[1].pressure_flow_equation.deactivate()
    m.fs.steam_turbine.throttle_valve[2].pressure_flow_equation.deactivate()
    m.fs.steam_turbine.throttle_valve[3].pressure_flow_equation.deactivate()
    m.fs.steam_turbine.throttle_valve[4].pressure_flow_equation.deactivate()
    @m.fs.steam_turbine.Constraint([2,3,4], m.fs.time)
    def same_throttle_pos_eqn(b, i, t):
        return b.throttle_valve[i].deltaP[t] == b.throttle_valve[1].deltaP[t]
    if not os.path.isfile("init_ngcc_st.json.gz"):
        solver.solve(m, tee=True)
    for i in m.fs.steam_turbine.same_throttle_pos_eqn:
        iscale.constraint_scaling_transform(
            m.fs.steam_turbine.same_throttle_pos_eqn[i], 1e-5)


    print("Connect t11")
    m.fs.t11 = Arc(
        source=m.fs.return_mix.outlet,
        destination=m.fs.LP_ECON.side_1_inlet)
    m.fs.LP_ECON.side_1_inlet.unfix()
    expand_arcs.apply_to(m)
    m.fs.cond_pump.control_volume.properties_out[0].pressure.fix(
        pyo.value(m.fs.LP_ECON.side_1_inlet.pressure[0]))
    m.fs.cond_pump.deltaP.unfix()
    imisc.copy_port_values(arc=m.fs.t11)
    if not os.path.isfile("init_ngcc_st.json.gz"):
        solver.solve(m, tee=True)

    print("connect HP steam")
    m.fs.t01 = Arc(
        source=m.fs.HP_SH4.side_1_outlet,
        destination=m.fs.steam_turbine.inlet_split.inlet)
    expand_arcs.apply_to(m)
    m.fs.hp_steam_temperature.unfix()
    m.fs.hotwell.makeup.flow_mol[:].unfix()
    m.fs.steam_turbine.inlet_split.inlet.unfix()
    imisc.copy_port_values(arc=m.fs.t01)
    if not os.path.isfile("init_ngcc_st.json.gz"):
        solver.solve(m, tee=True)

    print("Adjust HRSG Splits")
    m.fs.LP_EVAP.vapor_frac_control.fix(0.16)
    #m.fs.Splitter1.split_fraction[0, "toIP"].fix(0.16)
    #m.fs.LP_FGsplit.split_fraction[0, 'toLP_SH'].fix(0.4)
    #m.fs.steam_turbine.throttle_valve[1].deltaP[0].fix(-5e5)
    #m.fs.steam_turbine.outlet_stage.flow_coeff.fix(0.12)
    #m.fs.IP_Splitter1.split_fraction[0, "toNGPH"]
    if not os.path.isfile("init_ngcc_st.json.gz"):
        solver.solve(m, tee=True)
    m.fs.LP_EVAP.vapor_frac_control.fix(0.15)
    if not os.path.isfile("init_ngcc_st.json.gz"):
        solver.solve(m, tee=True)
    m.fs.LP_EVAP.vapor_frac_control.fix(0.14)
    if not os.path.isfile("init_ngcc_st.json.gz"):
        solver.solve(m, tee=True)
    m.fs.LP_EVAP.vapor_frac_control.fix(0.13)
    if not os.path.isfile("init_ngcc_st.json.gz"):
        solver.solve(m, tee=True)
    m.fs.LP_EVAP.vapor_frac_control.fix(0.12)
    if not os.path.isfile("init_ngcc_st.json.gz"):
        solver.solve(m, tee=True)

    if os.path.isfile("init_ngcc_st.json.gz"):
        ms.from_json(m, fname="init_ngcc_st.json.gz", wts=ms.StoreSpec(suffix=False))
    else:
        ms.to_json(m, fname="init_ngcc_st.json.gz")

    # Add the NG preheater
    m.fs.ng_preheater = gum.HeatExchanger(default={
        "delta_temperature_callback": delta_temperature_underwood_callback,
        "shell": {"property_package": m.fs.prop_water},
        "tube": {"property_package": m.fs.gas_prop_params}})

    m.fs.ng_preheater.area.fix(5000)
    m.fs.ng_preheater.overall_heat_transfer_coefficient.fix(100)

    #m.fs.ng_preheater.shell_inlet.fix()
    m.fs.ng_preheater.tube_inlet.fix()
    m.fs.Mixer1.Preheater.unfix()
    m.fs.feed_fuel1.flow_mol.unfix()
    m.fs.feed_fuel1.pressure.unfix()
    m.fs.feed_fuel1.temperature.unfix()
    m.fs.feed_fuel1.mole_frac_comp.unfix()
    m.fs.feed_fuel1.pressure.unfix()
    m.fs.feed_fuel1.mole_frac_comp.unfix()

    imisc.copy_port_values(
        source=m.fs.feed_fuel1.outlet,
        destination=m.fs.ng_preheater.tube_inlet)
    imisc.copy_port_values(
        source=m.fs.IP_Splitter1.toNGPH,
        destination=m.fs.ng_preheater.shell_inlet)


    m.fs.fuel01dupe = Arc(
        source=m.fs.ng_preheater.tube_outlet,
        destination=m.fs.inject1.gas)
    m.fs.st01 = Arc(
        source=m.fs.IP_Splitter1.toNGPH,
        destination=m.fs.ng_preheater.shell_inlet)
    m.fs.st02 = Arc(
        source=m.fs.ng_preheater.shell_outlet,
        destination=m.fs.Mixer1.Preheater)

    # add the preheater stream to tags
    st = { # streams that are half in HRSG, and may not have arcs
        "st01":m.fs.ng_preheater.shell.properties_in[0],
        "st02":m.fs.ng_preheater.shell.properties_out[0],
        "fuel02":m.fs.ng_preheater.tube.properties_in[0]}
    def new_tag(name, expr, format):
        # funcion to keep it more compact
        m.tags[name] = expr
        m.tag_format[name] = format
    for i, s in st.items(): # create the tags for steam quantities
        new_tag(f"{i}_Fvol", expr=s.flow_vol, format="{:.1f} m^3/s")
        new_tag(f"{i}_Fmol", expr=s.flow_mol/1000, format="{:.3f} kmol/s")
        new_tag(f"{i}_F", expr=s.flow_mass, format="{:.3f} kg/s")
        new_tag(f"{i}_P", expr=s.pressure/1000, format="{:.3f} kPa")
        new_tag(f"{i}_T", expr=s.temperature, format="{:.2f} K")
        new_tag(f"{i}_H", expr=s.enth_mol/1000, format="{:.2f} kJ/mol")
        if hasattr(s, "mole_frac_comp"):
            for c in s.mole_frac_comp:
                new_tag(
                    f"{i}_y{c}", expr=s.mole_frac_comp[c]*100, format="{:.3f}%")

    expand_arcs.apply_to(m)

    # Set scaling factors
    iscale.set_scaling_factor(m.fs.ng_preheater.shell.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.ng_preheater.tube.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.ng_preheater.area, 1e-2)
    iscale.set_scaling_factor(
        m.fs.ng_preheater.overall_heat_transfer_coefficient, 1e-2)
    iscale.set_scaling_factor(m.fs.hot_reheat_temperature, 1e-2)
    iscale.set_scaling_factor(m.fs.net_power_mw, 1e-2)
    #iscale.set_scaling_factor(m.fs.natural_gas, 1e-4)
    iscale.set_scaling_factor(m.fs.HP_pump.control_volume.deltaP, 1e-8)
    iscale.set_scaling_factor(m.fs.LP_Mixer2.minimum_pressure, 1e-5)
    iscale.set_scaling_factor(m.fs.main_condenser.shell.heat, 1e-8)
    iscale.set_scaling_factor(m.fs.main_condenser.tube.heat, 1e-8)
    for t in m.fs.main_condenser.unit_heat_balance:
        iscale.constraint_scaling_transform(m.fs.main_condenser.unit_heat_balance[t], 1e-8, overwrite=True)
    for t in m.fs.main_condenser.heat_transfer_equation:
        iscale.constraint_scaling_transform(m.fs.main_condenser.heat_transfer_equation[t], 1e-8, overwrite=True)

    iscale.set_scaling_factor(m.fs.main_condenser.shell.properties_in[0.0].pressure, 1e-3)
    iscale.set_scaling_factor(m.fs.main_condenser.shell.properties_out[0.0].pressure, 1e-3)
    iscale.set_scaling_factor(m.fs.cond_pump.control_volume.properties_in[0.0].pressure, 1e-3)
    iscale.set_scaling_factor(m.fs.hotwell.condensate_state[0.0].pressure, 1e-3)
    iscale.set_scaling_factor(m.fs.steam_turbine.outlet_stage.control_volume.properties_out[0.0].pressure, 1e-3)

    iscale.calculate_scaling_factors(m)
    m#.fs.IP_Splitter1.split_fraction[0, "toNGPH"].unfix()
    m.fs.ng_preheater.tube_inlet.flow_mol.unfix()
    #m.fs.feed_fuel1.temperature.fix(460)
    m.fs.ng_preheater.tube_inlet.temperature.fix(311)
    if not os.path.isfile("init_preheat.json.gz"):
        m.fs.ng_preheater.initialize()
        imisc.copy_port_values(
            source=m.fs.ng_preheater.shell_outlet,
            destination=m.fs.Mixer1.Preheater)
        solver.solve(m, tee=True)

    m.fs.cond_pump.control_volume.properties_out[0].pressure.fix(655000)
    m.fs.cmp1.efficiency_isentropic.fix(0.85)
    for c in [m.fs.HP_SH1.fcorrection_htc, m.fs.HP_SH2.fcorrection_htc,
        m.fs.HP_SH3.fcorrection_htc, m.fs.HP_SH4.fcorrection_htc]:
        c.value = pyo.value(c * 0.75)
    for c in [m.fs.IP_SH1.fcorrection_htc, m.fs.IP_SH2.fcorrection_htc,
        m.fs.IP_SH3.fcorrection_htc]:
        c.value = pyo.value(c * 1.45)
    for c in [m.fs.LP_SH.fcorrection_htc]:
        c.value = pyo.value(c * 1.1)
    for i, s in m.fs.steam_turbine.lp_stages.items():
        s.ratioP[:] = 0.754
    if not os.path.isfile("init_preheat.json.gz"):
        solver.solve(m, tee=True)

    if os.path.isfile("init_preheat.json.gz"):
        ms.from_json(m, fname="init_preheat.json.gz", wts=ms.StoreSpec(suffix=False))
    else:
        ms.to_json(m, fname="init_preheat.json.gz")

    # add a heater for the reboiler.
    print("Adding reboiler heater")
    m.fs.reboiler = gum.Heater(default={"property_package": m.fs.prop_water})
    @m.fs.reboiler.Constraint(m.fs.config.time)
    def reboiler_condense_eqn(b, t):
        return b.control_volume.properties_out[t].enth_mol == \
            b.control_volume.properties_out[t].enth_mol_sat_phase["Liq"] - 100
    m.fs.return_mix.reboiler.unfix()
    imisc.copy_port_values(
        source=m.fs.return_mix.reboiler,
        destination=m.fs.reboiler.outlet)
    imisc.copy_port_values(
        source=m.fs.steam_turbine_lp_split.reboiler,
        destination=m.fs.reboiler.inlet)
    m.fs.t17 = Arc(
        source=m.fs.reboiler.outlet,
        destination=m.fs.return_mix.reboiler)
    m.fs.t13 = Arc(
        source=m.fs.steam_turbine_lp_split.reboiler,
        destination=m.fs.reboiler.inlet)
    expand_arcs.apply_to(m)
    iscale.set_scaling_factor(m.fs.reboiler.control_volume.heat, 1e-8)
    iscale.calculate_scaling_factors(m)
    if not os.path.isfile("init_reboiler.json.gz"):
        solver.solve(m, tee=True)

    if os.path.isfile("init_reboiler.json.gz"):
        ms.from_json(m, fname="init_reboiler.json.gz", wts=ms.StoreSpec(suffix=False))
    else:
        ms.to_json(m, fname="init_reboiler.json.gz")


    # Aux power expressions
    @m.fs.Expression(m.fs.config.time)
    def gross_power(b, t):
        return b.gt_power[t] + b.steam_turbine.power[t]

    @m.fs.Expression(m.fs.config.time)
    def aux_cooling(b, t):
        return 1e3*(4580 + 2370)

    # Aux power expressions
    @m.fs.Expression(m.fs.config.time)
    def aux_combustion(b, t):
        return 1e3*1020

    @m.fs.Expression(m.fs.config.time)
    def aux_capture(b, t): #scale to flue gas flow
        return 1e3*10600* \
            (b.gts2.control_volume.properties_out[t].flow_mass/1090.759)

    @m.fs.Expression(m.fs.config.time)
    def aux_compression(b, t): #scale to flue gas flow
        return 1e3*17090* \
            (b.gts2.control_volume.properties_out[0].flow_mol_comp["CO2"]/1547.75)

    @m.fs.Expression(m.fs.config.time)
    def aux_transformer(b, t): # scale to gross power
        return 1e3*2200 * \
            (b.gross_power[t]/687.0e6)

    @m.fs.Expression(m.fs.config.time)
    def aux_misc(b, t):
        return 1e3*1000.

    @m.fs.Expression(m.fs.config.time)
    def net_power(b, t):
        return (
            b.gt_power[t] +
            b.steam_turbine.power[t] +
            b.cond_pump.work[t] +
            b.HP_pump.work[t] +
            b.IP_pump.work[t] +
            b.aux_cooling[t] +
            b.aux_combustion[t] +
            b.aux_capture[t] +
            b.aux_compression[t] +
            b.aux_transformer[t] +
            b.aux_misc[t])

    m.fs.fuel_lhv = pyo.Var(initialize=47.2e6,
                            units=pyunits.J/pyunits.kg) # J/kg
    m.fs.fuel_lhv.fix(47.2e6)

    @m.fs.Expression(m.fs.config.time)
    def fuel_thermal_in_mbtu(b, t):
        return pyunits.convert(m.fs.fuel_lhv*m.fs.inject1.gas_state[t].flow_mass, pyunits.MBtu/pyunits.hr)

    @m.fs.Expression(m.fs.config.time)
    def lhv_efficiency(b, t):
        return -b.net_power[t]/b.inject1.gas_state[t].flow_mass/b.fuel_lhv

    @m.fs.Expression(m.fs.config.time)
    def reboiler_duty_expr(b, t): #scale to flue gas flow
        return -176.35920313360697e6 * \
            (b.gts2.control_volume.properties_out[t].flow_mass/1090.759)

    print("Control steam to maintain temp")
    #m.fs.steam_turbine.throttle_valve[1].deltaP[0].unfix()
    m.fs.HP_pump.outlet.pressure[0].unfix()
    m.fs.hp_steam_temperature.fix(855)
    if not os.path.isfile("init_control_steam.json.gz"):
        solver.solve(m, tee=True)

    if os.path.isfile("init_control_steam.json.gz"):
        ms.from_json(m, fname="init_control_steam.json.gz", wts=ms.StoreSpec(suffix=False))
    else:
        ms.to_json(m, fname="init_control_steam.json.gz")

    print("Reboiler Heat Constraint")
    @m.fs.Constraint(m.fs.config.time)
    def reboiler_duty_eqn(b, t):
        return b.reboiler_duty_expr[t] == m.fs.reboiler.control_volume.heat[t]
    for i in m.fs.reboiler_duty_eqn:
        iscale.constraint_scaling_transform(m.fs.reboiler_duty_eqn[i], 1e-8)
    m.fs.steam_turbine_lp_split.reboiler.flow_mol.unfix()
    if not os.path.isfile("init_reboiler_heat.json.gz"):
        solver.solve(m, tee=True)
    if os.path.isfile("init_reboiler_heat.json.gz"):
        ms.from_json(m, fname="init_reboiler_heat.json.gz", wts=ms.StoreSpec(suffix=False))
    else:
        ms.to_json(m, fname="init_reboiler_heat.json.gz")

    print("Add costing and net power")
    m.fs.LP_FGsplit.split_fraction[:, "toLP_SH"].fix(0.55)
    @m.fs.Constraint(m.fs.config.time)
    def net_power_constraint(b, t):
        return b.net_power_mw[t]/100.0 == -b.net_power[t]/1e6/100.0
    # adding Total Plant Cost and Fixed and Variable O&M costs
    get_ngcc_costing(m, evaluate_cost=True)
    build_ngcc_OM_costs(m)
    @m.fs.Constraint(m.fs.time)
    def eq1(c, t):
        return m.fs.natural_gas[t] == pyunits.convert(
            m.fs.inject1.gas_state[t].flow_mass*m.fs.fuel_lhv,
            pyunits.MBtu/pyunits.day)
    for t in m.fs.eq1:
        iscale.constraint_scaling_transform(m.fs.eq1[t], 1e-4)
    m.fs.natural_gas.unfix()
    @m.fs.Constraint(m.fs.time)
    def eq2(c, t):
        return m.fs.net_power_cost[t] == m.fs.net_power_mw[t]
    for t in m.fs.eq2:
        iscale.constraint_scaling_transform(m.fs.eq2[t], 1e-2)
    m.fs.net_power_cost.unfix()

    if not os.path.isfile("init_add_costing.json.gz"):
        solver.solve(m, tee=True)
    if os.path.isfile("init_add_costing.json.gz"):
        ms.from_json(m, fname="init_add_costing.json.gz", wts=ms.StoreSpec(suffix=False))
    else:
        ms.to_json(m, fname="init_add_costing.json.gz")


    iscale.set_scaling_factor(m.fs.costing.total_TPC, 1e-2)
    iscale.set_scaling_factor(m.fs.natural_gas, 1e-4)
    iscale.set_scaling_factor(m.fs.net_power_cost, 1e-2)
    iscale.set_scaling_factor(m.fs.HP_ECON1.side_1.deltaP[0.0], 1e-7)
    iscale.set_scaling_factor(m.fs.HP_ECON2.side_1.deltaP[0.0], 1e-7)
    iscale.set_scaling_factor(m.fs.HP_ECON3.side_1.deltaP[0.0], 1e-7)
    iscale.set_scaling_factor(m.fs.HP_ECON4.side_1.deltaP[0.0], 1e-7)
    iscale.set_scaling_factor(m.fs.HP_ECON5.side_1.deltaP[0.0], 1e-7)
    iscale.set_scaling_factor(m.fs.HP_ECON1.side_2.deltaP[0.0], 1)
    iscale.set_scaling_factor(m.fs.HP_ECON2.side_2.deltaP[0.0], 1)
    iscale.set_scaling_factor(m.fs.HP_ECON3.side_2.deltaP[0.0], 1)
    iscale.set_scaling_factor(m.fs.HP_ECON4.side_2.deltaP[0.0], 1)
    iscale.set_scaling_factor(m.fs.HP_ECON5.side_2.deltaP[0.0], 1)
    # Set some constraint scaling factors that I want to override calc


    for i, c in m.fs.cond_pump.eq_work.items():
        iscale.constraint_scaling_transform(c, 1e-4, overwrite=True)
    for i, c in m.fs.IP_EVAP.sat_vap_eqn.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for i, c in m.fs.HP_EVAP.sat_vap_eqn2.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for i, c in m.fs.HP_ECON1.side_1.pressure_balance.items():
        iscale.constraint_scaling_transform(c, 1e-6, overwrite=True)
    for i, c in m.fs.HP_ECON2.side_1.pressure_balance.items():
        iscale.constraint_scaling_transform(c, 1e-6, overwrite=True)
    for i, c in m.fs.HP_ECON3.side_1.pressure_balance.items():
        iscale.constraint_scaling_transform(c, 1e-6, overwrite=True)
    for i, c in m.fs.HP_ECON4.side_1.pressure_balance.items():
        iscale.constraint_scaling_transform(c, 1e-6, overwrite=True)
    for i, c in m.fs.HP_ECON5.side_1.pressure_balance.items():
        iscale.constraint_scaling_transform(c, 1e-6, overwrite=True)
    for i, c in m.fs.HP_pump.control_volume.pressure_balance.items():
        iscale.constraint_scaling_transform(c, 1e-7, overwrite=True)
    for i, c in m.fs.reheat_steam_temperature_eqn.items():
        iscale.constraint_scaling_transform(c, 1e-2)
    for i, c in m.fs.LP_Mixer2.minimum_pressure_constraint.items():
        iscale.constraint_scaling_transform(c, 1e-5)
    for i, c in m.fs.LP_Mixer2.mixture_pressure.items():
        iscale.constraint_scaling_transform(c, 1e-5)
    for i, c in m.fs.LP_Mixer2.enthalpy_mixing_equations.items():
        iscale.constraint_scaling_transform(c, 1e-8)
    iscale.constraint_scaling_transform(m.fs.costing.annual_labor_cost_rule, 1e-6)
    for t in m.fs.reboiler.reboiler_condense_eqn:
        iscale.constraint_scaling_transform(m.fs.reboiler.reboiler_condense_eqn[t], 1e-3)
    for t, c in m.fs.costing.total_TPC_eq.items():
        iscale.constraint_scaling_transform(c, 1e-2)

    return m, solver

def tabulated_output_dict(m, add=[]):
    d = OrderedDict([
        ("Net Power (MW)", -m.fs.net_power[0]/1e6),
        ("Fuel Flow (kg/s)", m.fs.inject1.gas_state[0].flow_mass),
        ("Gas Turbine Power (MW)", -m.fs.gt_power[0]/1e6),
        ("Steam Turbine Power (MW)", -m.fs.steam_turbine.power[0]/1e6),
        ("LHV Efficiency (%)", m.fs.lhv_efficiency[0]*100),
        ("Capture Reboiler Duty (MW)", m.fs.reboiler_duty_expr[0]/1e6),
        ("Total Plant Cost ($/hr)", m.fs.costing.total_TPC/365/24*1e6),
        ("Fixed O&M Cost ($/hr)", m.fs.costing.total_fixed_OM_cost/365/24*1e6),
        ("Variable O&M Cost ($/hr)", m.fs.costing.total_variable_OM_cost[0]*(-m.fs.net_power[0]/1e6)),
    ])
    for a in add:
        d[a.getname()] = a
    return d

def result_summary(m, add=[]):
    d = tabulated_output_dict(m, add=add)
    print("")
    for k, v, in d.items():
        print(f"{k}: {pyo.value(v)}")
    print("")

def write_csv_header(filename, d):
    with open(filename, "w", newline='') as f:
        w = csv.writer(f)
        w.writerow(list(d.keys()))


def write_csv_row(filename, d):
    with open(filename, "a", newline='') as f:
        w = csv.writer(f)
        w.writerow([pyo.value(x) for x in d.values()])


def write_pfds(m, aname=""):
    gas_turbine.write_pfd_results(f"pfd_results_gt{aname}.svg", m.tags, m.tag_format, infilename="gas_turbine.svg")
    sturb_module.write_pfd_results(f"pfd_results_st{aname}.svg", m.tags, m.tag_format)
    hrsg_module.pfd_result(f"pfd_results_hrsg{aname}.svg", m)


def run_series(m):
    m.fs.gt_power.unfix()
    netpow = np.linspace(650, 150, 50).tolist()
    d = tabulated_output_dict(m)
    solver = pyo.SolverFactory("ipopt")
    write_csv_header("series.csv", d)
    for p in netpow:
        p = int(p)
        print(f"Start Net Power {p} Run")
        m.fs.net_power_mw.fix(p)
        sf = f"state_ngcc_state_{p}.json.gz"
        if not os.path.isfile(sf):
            solver.solve(m, tee=True)
            ms.to_json(m, fname=sf)
        else:
            ms.from_json(m, fname=sf, wts=ms.StoreSpec(suffix=False))
        write_csv_row("series.csv", d)
        result_summary(m)
        write_pfds(m, aname=f"_{p}")


def run_opt_series(m):
    m.fs.gt_power.unfix()
    #m.obj = pyo.Objective(expr=m.fs.costing.total_variable_OM_cost[0])
    m.obj = pyo.Objective(expr=-m.fs.lhv_efficiency[0]*100)


    dvar_list = []
    def make_dvar(v, l, u):
        v.unfix()
        v.setlb(l)
        v.setub(u)
        dvar_list.append(v)

    #make_dvar(m.fs.exhaust_1.temperature[0], 890, 910)
    m.fs.upper_reheat_temp = pyo.Constraint(expr=
        m.fs.IP_SH3.side_1.properties_out[0].temperature <= 860)
    #m.fs.lower_reheat_temp = pyo.Constraint(expr=
    #    m.fs.IP_SH3.side_1.properties_out[0].temperature >= 830)
    #make_dvar(m.fs.hp_steam_temperature[0], 845, 860)
    make_dvar(m.fs.LP_EVAP.vapor_frac_control, 0.10, 0.15)
    make_dvar(m.fs.IP_Splitter1.split_fraction[0, "toNGPH"], 0.35, 0.55)
    make_dvar(m.fs.LP_FGsplit.split_fraction[0, "toLP_SH"], 0.3, 0.7)
    make_dvar(m.fs.Splitter1.split_fraction[0, "toIP"], 0.15, 0.25)
    solver = pyo.SolverFactory("ipopt")
    solver.options["max_iter"] = 100
    #solver.options["halt_on_ampl_error"] = "yes"
    netpow = np.linspace(650, 150, 50).tolist()
    d = tabulated_output_dict(m, add=dvar_list)
    write_csv_header("opt_series.csv", d)
    linear_eliminate = pyo.TransformationFactory("simple_equality_eliminator")
    for p in netpow:
        p = int(p)
        print(f"Start Net Power {p} Run")
        m.fs.net_power_mw.fix(p)
        sf = f"state_ngcc_ostate_{p}.json.gz"
        if not os.path.isfile(sf):
            linear_eliminate.apply_to(m, max_iter=15, reversible=True)
            solver.solve(m, tee=True, symbolic_solver_labels=True)
            linear_eliminate.revert()
            ms.to_json(m, fname=sf)
        else:
            ms.from_json(m, fname=sf, wts=ms.StoreSpec(suffix=False))
        result_summary(m, add=dvar_list)
        write_csv_row("opt_series.csv", d)
        write_pfds(m, aname=f"_o{p}")



def check_scaling(m):
    jac, nlp = iscale.get_jacobian(m, scaled=True)
    print("Extreme Jacobian entries:")
    for i in iscale.extreme_jacobian_entries(jac=jac, nlp=nlp, small=1e-6, large=100):
        print(f"    {i[0]:.2e}, [{i[1]}, {i[2]}]")
    print("Unscaled constraints:")
    for c in iscale.unscaled_constraints_generator(m):
        print(f"    {c}")
    print("Scaled constraints by factor:")
    for c, s in iscale.constraints_with_scale_factor_generator(m):
        print(f"    {c}, {s}")
    print("Badly scaled variables:")
    for v, sv in iscale.badly_scaled_var_generator(m, large=1e2, small=1e-2, zero=1e-12):
        print(f"    {v} -- {sv} -- {iscale.get_scaling_factor(v)}")
    print(f"Jacobian Condition Number: {iscale.jacobian_cond(jac=jac):.2e}")


if __name__ == "__main__":
    m, solver = get_model()
    solver.solve(m, tee=True)
    write_pfds(m)
    result_summary(m)
    check_scaling(m)
    #run_opt_series(m)
