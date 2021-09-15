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

__author__ = "John Eslick"

import os
import csv

import numpy as np
import matplotlib.pyplot as plt

import pyomo.environ as pyo
from pyomo.network import Arc, Port
from pyomo.common.fileutils import this_file_dir
import pyomo.common.errors

from idaes.core import FlowsheetBlock
from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock,
)
from idaes.generic_models.properties.core.generic.generic_reaction import (
    GenericReactionParameterBlock,
)
from idaes.power_generation.unit_models.helm import (
    HelmMixer,
    MomentumMixingType,
    HelmSplitter,
    HelmIsentropicCompressor
)
import idaes.generic_models.unit_models as gum  # generic unit models
import idaes.power_generation.unit_models as pum  # power unit models
import idaes.core.util as iutil
import idaes.core.util.tables as tables
import idaes.core.util.scaling as iscale
import idaes.core.util.initialization as iinit

import idaes.core.plugins
from idaes.power_generation.properties.natural_gas_PR import get_prop, get_rxn, EosType
from idaes.core.solvers import use_idaes_solver_configuration_defaults
from idaes.power_generation.properties import FlueGasParameterBlock
from idaes.generic_models.properties import iapws95

from idaes.generic_models.unit_models.heat_exchanger import (
    delta_temperature_underwood_callback,
    delta_temperature_lmtd_callback,
    delta_temperature_lmtd2_callback,
    delta_temperature_lmtd3_callback,
)

from idaes.generic_models.properties.helmholtz.helmholtz import (
    HelmholtzThermoExpressions as ThermoExpr,
)
import idaes.logger as idaeslog


def _set_port(port, F, T, P, comp, fix=True):
    if fix:
        port.flow_mol.fix(F)
        port.temperature.fix(T)
        port.pressure.fix(P)
        for k, v in comp.items():
            port.mole_frac_comp[:, k].fix(v)
    else:
        port.flow_mol[:].value = F
        port.temperature[:].value = T
        port.pressure[:] = P
        for k, v in comp.items():
            port.mole_frac_comp[:, k].value = v


def add_flowsheet(m=None, name="SOEC Module"):
    if m is None:
        m = pyo.ConcreteModel(name)
    if not hasattr(m, "fs"):
        m.fs = FlowsheetBlock(default={"dynamic": False})
    return m


def add_properties(m):
    comps = {  # components present
        "CH4",
        "C2H6",
        "C3H8",
        "C4H10",
        "O2",
        "H2O",
        "CO2",
        "N2",
        "Ar",
    }
    m.fs.fg_prop = GenericParameterBlock(
        default=get_prop(components=comps, phases=["Vap"], eos=EosType.IDEAL)
    )
    m.fs.fg_prop.set_default_scaling("mole_frac_comp", 10)
    m.fs.fg_prop.set_default_scaling("mole_frac_phase_comp", 10)
    rxns = {  # reactions and key components for conversion
        "ch4_cmb": "CH4",
        "c2h6_cmb": "C2H6",
        "c3h8_cmb": "C3H8",
        "c4h10_cmb": "C4H10",
    }
    m.rxns = rxns
    m.fs.fg_combust = GenericReactionParameterBlock(default=get_rxn(m.fs.fg_prop, rxns))
    m.fs.water_prop = iapws95.Iapws95ParameterBlock()
    m.fs.h2_prop = GenericParameterBlock(
        default=get_prop(components={"H2", "H2O"}, phases=["Vap"], eos=EosType.IDEAL)
    )
    m.fs.o2_prop = GenericParameterBlock(
        default=get_prop(components={"O2", "H2O"}, phases=["Vap"], eos=EosType.IDEAL)
    )


def add_preheater(m):
    m.fs.air_preheater = gum.HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_underwood_callback,
            "shell": {"property_package": m.fs.fg_prop},
            "tube": {"property_package": m.fs.fg_prop},
        }
    )
    m.fs.ng_preheater = gum.HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_underwood_callback,
            "shell": {"property_package": m.fs.fg_prop},
            "tube": {"property_package": m.fs.fg_prop},
        }
    )
    m.fs.preheat_split = gum.Separator(
        default={
            "property_package": m.fs.fg_prop,
            "outlet_list": ["air", "ng"],
        }
    )
    m.fs.fg04 = Arc(
        source=m.fs.preheat_split.ng, destination=m.fs.ng_preheater.shell_inlet
    )
    m.fs.fg05 = Arc(
        source=m.fs.preheat_split.air, destination=m.fs.air_preheater.shell_inlet
    )


def add_combustor(m):
    m.fs.cmb_mix = gum.Mixer(default={
        "property_package": m.fs.fg_prop,
        "inlet_list":["ng", "air"],
        "momentum_mixing_type":gum.MomentumMixingType.none})
    m.fs.cmb = gum.StoichiometricReactor(default={
        "property_package": m.fs.fg_prop,
        "reaction_package": m.fs.fg_combust,
        "has_pressure_change": False})
    @m.fs.cmb_mix.Constraint(m.fs.time)
    def pressure_eqn(b, t):
        return b.mixed_state[t].pressure == b.air_state[t].pressure
    @m.fs.cmb.Constraint(m.fs.time, m.rxns.keys())
    def reaction_extent(b, t, r):
        k = m.rxns[r]
        prp = b.control_volume.properties_in[t]
        stc = -m.fs.fg_combust.rate_reaction_stoichiometry[r, "Vap", k]
        extent = b.rate_reaction_extent[t, r]
        return extent == prp.flow_mol*prp.mole_frac_comp[k]/stc
    m.fs.ba03 = Arc(
        source=m.fs.air_preheater.tube_outlet, destination=m.fs.cmb_mix.air
    )
    m.fs.bng03 = Arc(
        source=m.fs.ng_preheater.tube_outlet, destination=m.fs.cmb_mix.ng
    )
    m.fs.bng04 = Arc(
        source=m.fs.cmb_mix.outlet, destination=m.fs.cmb.inlet
    )


def add_aux_boiler_steam(m):
    m.fs.bhx2 = gum.HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_underwood_callback,
            "shell": {"property_package": m.fs.fg_prop},
            "tube": {"property_package": m.fs.water_prop},
        }
    )
    m.fs.main_steam_split = HelmSplitter(
        default={
            "property_package": m.fs.water_prop,
            "outlet_list": ["h_side", "o_side"]
        }
    )
    m.fs.aux_boiler_feed_pump = HelmIsentropicCompressor(
        default={"property_package": m.fs.water_prop}
    )
    m.fs.bhx1 = gum.HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_underwood_callback,
            "shell": {"property_package": m.fs.fg_prop},
            "tube": {"property_package": m.fs.water_prop},
        }
    )
    m.fs.recover_split = HelmSplitter(
        default={
            "property_package": m.fs.water_prop,
            "outlet_list": ["h_side", "o_side"]
        }
    )
    m.fs.fg01 = Arc(source=m.fs.cmb.outlet, destination=m.fs.bhx2.shell_inlet)
    m.fs.s02 = Arc(
        source=m.fs.aux_boiler_feed_pump.outlet,
        destination=m.fs.bhx1.tube_inlet
    )
    m.fs.fg02 = Arc(source=m.fs.bhx2.shell_outlet, destination=m.fs.bhx1.shell_inlet)
    m.fs.s03 = Arc(
        source=m.fs.bhx1.tube_outlet,
        destination=m.fs.recover_split.inlet
    )
    m.fs.s09 = Arc(
        source=m.fs.bhx2.tube_outlet,
        destination=m.fs.main_steam_split.inlet
    )


def add_recovery_hx(m):
    m.fs.hxh2 = gum.HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_underwood_callback,
            "shell": {"property_package": m.fs.h2_prop},
            "tube": {"property_package": m.fs.water_prop},
        }
    )
    m.fs.hxo2 = gum.HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_underwood_callback,
            "shell": {"property_package": m.fs.o2_prop},
            "tube": {"property_package": m.fs.water_prop},
        }
    )
    m.fs.smix1 = HelmMixer(
        default={
            "momentum_mixing_type": MomentumMixingType.none,
            "inlet_list": ["hxh2", "hxo2"],
            "property_package": m.fs.water_prop,
        }
    )
    m.fs.h02 = Arc(source=m.fs.spltf1.out, destination=m.fs.hxh2.shell_inlet)
    m.fs.s05 = Arc(source=m.fs.recover_split.h_side, destination=m.fs.hxh2.tube_inlet)
    m.fs.o02 = Arc(source=m.fs.splta1.out, destination=m.fs.hxo2.shell_inlet)
    m.fs.s04 = Arc(source=m.fs.recover_split.o_side, destination=m.fs.hxo2.tube_inlet)
    m.fs.s06 = Arc(source=m.fs.hxo2.tube_outlet, destination=m.fs.smix1.hxo2)
    m.fs.s07 = Arc(source=m.fs.hxh2.tube_outlet, destination=m.fs.smix1.hxh2)

    @m.fs.smix1.Constraint(m.fs.time)
    def pressure_eqn(b, t):
        return b.mixed_state[t].pressure == b.hxh2_state[t].pressure


def add_soec_unit(m):
    m.fs.soec = pum.IsothermalSofc(
        default={
            "nz": 20,
            "nxfe": 10,
            "nxae": 10,
            "soec": True,
            "air_side_comp_list": ["H2O", "O2"],
            "fuel_side_comp_list": ["H2O", "H2"],
            "air_side_stoich": {"H2O": 0, "O2": -0.25},
        }
    )
    m.fs.spltf1 = gum.Separator(
        default={
            "property_package": m.fs.h2_prop,
            "outlet_list": ["out", "recycle"],
        }
    )
    m.fs.splta1 = gum.Separator(
        default={
            "property_package": m.fs.o2_prop,
            "outlet_list": ["out", "recycle"],
        }
    )
    m.fs.h01 = Arc(source=m.fs.soec.outlet_fc_mult, destination=m.fs.spltf1.inlet)
    m.fs.o01 = Arc(source=m.fs.soec.outlet_ac_mult, destination=m.fs.splta1.inlet)

    m.fs.soec.E_cell.fix(1.28)  # unfix after initialize
    m.fs.soec.el.thickness.fix(9e-6)
    m.fs.soec.fe.thickness.fix(1e-3)
    m.fs.soec.ae.thickness.fix(20e-6)
    m.fs.soec.length.fix(0.05)
    m.fs.soec.width.fix(0.05)
    m.fs.soec.k_ae.fix(26.1e7)
    m.fs.soec.eact_ae.fix(120000)
    m.fs.soec.alpha_ae.fix(0.4)
    m.fs.soec.k_fe.fix(1.35e10)
    m.fs.soec.eact_fe.fix(110000)
    m.fs.soec.alpha_fe.fix(0.5)
    m.fs.soec.fe.k_res.fix(2.98e-5)
    m.fs.soec.fe.E_res.fix(-1392)
    m.fs.soec.ae.k_res.fix(8.114e-5)
    m.fs.soec.ae.E_res.fix(600)
    m.fs.soec.el.k_res.fix(2.94e-5)
    m.fs.soec.el.E_res.fix(10350)
    m.fs.soec.fc.thickness.fix(0.002)
    m.fs.soec.ac.thickness.fix(0.002)
    m.fs.soec.fe.porosity.fix(0.48)
    m.fs.soec.fe.tortuosity.fix(5.4)
    m.fs.soec.ae.porosity.fix(0.48)
    m.fs.soec.ae.tortuosity.fix(5.4)
    temperature = 1073.15
    m.fs.soec.el.temperature.fix(temperature)
    m.fs.soec.fc.temperature.fix(temperature)
    m.fs.soec.ac.temperature.fix(temperature)
    m.fs.soec.fe.temperature.fix(temperature)
    m.fs.soec.ae.temperature.fix(temperature)


def add_soec_inlet_mix(m):
    m.fs.mxf1 = gum.Mixer(
        default={
            "property_package": m.fs.h2_prop,
            "inlet_list": ["water", "recycle"],
            "momentum_mixing_type": gum.MomentumMixingType.none,
        }
    )
    m.fs.mxa1 = gum.Mixer(
        default={
            "property_package": m.fs.o2_prop,
            "inlet_list": ["water", "recycle"],
            "momentum_mixing_type": gum.MomentumMixingType.none,
        }
    )

    @m.fs.mxf1.Constraint(m.fs.time)
    def fmxpress_eqn(b, t):
        return b.mixed_state[t].pressure == b.water_state[t].pressure

    @m.fs.mxa1.Constraint(m.fs.time)
    def amxpress_eqn(b, t):
        return b.mixed_state[t].pressure == b.water_state[t].pressure

    # Add ports to connet pure steam to steam + h2 or steam + o2
    m.fs.main_steam_split._temperature_h_side_ref = pyo.Reference(
        m.fs.main_steam_split.h_side_state[:].temperature
    )
    m.fs.main_steam_split._temperature_o_side_ref = pyo.Reference(
        m.fs.main_steam_split.o_side_state[:].temperature
    )
    @m.fs.main_steam_split.Expression(m.fs.time, m.fs.soec.fc.config.comp_list)
    def h_side_mole_frac_expr(b, t, i):
        if i == "H2O":
            return 1
        else:
            return 0
    @m.fs.main_steam_split.Expression(m.fs.time, m.fs.soec.ac.config.comp_list)
    def o_side_mole_frac_expr(b, t, i):
        if i == "H2O":
            return 1
        else:
            return 0

    m.fs.main_steam_split.h_side_adapt = Port(
        rule=lambda b: {
            "flow_mol": m.fs.main_steam_split._flow_mol_h_side_ref,
            "pressure": m.fs.main_steam_split._pressure_h_side_ref,
            "temperature": m.fs.main_steam_split._temperature_h_side_ref,
            "mole_frac_comp": m.fs.main_steam_split.h_side_mole_frac_expr,
        }
    )
    m.fs.main_steam_split.o_side_adapt = Port(
        rule=lambda b: {
            "flow_mol": m.fs.main_steam_split._flow_mol_o_side_ref,
            "pressure": m.fs.main_steam_split._pressure_o_side_ref,
            "temperature": m.fs.main_steam_split._temperature_o_side_ref,
            "mole_frac_comp": m.fs.main_steam_split.o_side_mole_frac_expr,
        }
    )

    m.fs.s10 = Arc(
        source=m.fs.main_steam_split.h_side_adapt,
        destination=m.fs.mxf1.water,
    )
    m.fs.s11 = Arc(
        source=m.fs.main_steam_split.o_side_adapt,
        destination=m.fs.mxa1.water,
    )
    m.fs.hr01 = Arc(
        source=m.fs.spltf1.recycle,
        destination=m.fs.mxf1.recycle,
    )
    m.fs.or01 = Arc(
        source=m.fs.splta1.recycle,
        destination=m.fs.mxa1.recycle,
    )
    m.fs.s12 = Arc(source=m.fs.mxf1.outlet, destination=m.fs.soec.inlet_fc_mult)
    m.fs.s13 = Arc(source=m.fs.mxa1.outlet, destination=m.fs.soec.inlet_ac_mult)


def add_more_hx_connections(m):
    m.fs.fg03 = Arc(
        source=m.fs.bhx1.shell_outlet,
        destination=m.fs.preheat_split.inlet
    )
    m.fs.s08 = Arc(
        source=m.fs.smix1.outlet,
        destination=m.fs.bhx2.tube_inlet
    )


def add_constraints(m):
    m.fs.soec_heat_duty = pyo.Var(m.fs.time, units=pyo.units.W)
    @m.fs.Constraint(m.fs.time)
    def heat_duty_soec_zero_eqn(b, t):
        return b.soec.heat_duty[t] == b.soec_heat_duty[t]

    m.fs.soec_cmb_temperature = pyo.Var(m.fs.time, initialize=2000, units=pyo.units.K)
    @m.fs.Constraint(m.fs.time)
    def soec_cmb_temperature_eqn(b, t):
        return m.fs.cmb.outlet.temperature[t] == m.fs.soec_cmb_temperature[t]

    m.fs.soec_steam_temperature = pyo.Var(m.fs.time, initialize=1073.15, units=pyo.units.K)
    @m.fs.Constraint(m.fs.time)
    def soec_steam_temperature_eqn(b, t):
        return m.fs.bhx2.tube.properties_out[t].temperature == m.fs.soec_steam_temperature[t]

def set_guess(m):
    fg_comp_guess = {
        "CH4": 0.0,
        "C2H6": 0.0,
        "C3H8": 0.0,
        "C4H10": 0.0,
        "O2": 0.049,
        "H2O": 0.2,
        "CO2": 0.2,
        "N2": 0.55,
        "Ar": 0.001,
    }
    _set_port(
        m.fs.preheat_split.inlet, F=650, T=550, P=1.04e5, comp=fg_comp_guess, fix=True
    )


def set_inputs(m):
    m.fs.soec_cmb_temperature.fix(2000)
    m.fs.soec_steam_temperature.fix(1073.15)
    m.fs.air_preheater.area.fix(3000)
    m.fs.air_preheater.overall_heat_transfer_coefficient.fix(100)
    m.fs.ng_preheater.area.fix(300)
    m.fs.ng_preheater.overall_heat_transfer_coefficient.fix(100)
    m.fs.bhx2.area.fix(3000)
    m.fs.bhx2.overall_heat_transfer_coefficient.fix(100)
    m.fs.bhx1.area.fix(500)
    m.fs.bhx1.overall_heat_transfer_coefficient.fix(100)
    m.fs.hxh2.area.fix(4000)
    m.fs.hxh2.overall_heat_transfer_coefficient.fix(100)
    m.fs.hxo2.area.fix(4000)
    m.fs.hxo2.overall_heat_transfer_coefficient.fix(100)


    m.fs.preheat_split.split_fraction[:, "air"].fix(0.9)

    air_comp = {
        "CH4": 0.0,
        "C2H6": 0.0,
        "C3H8": 0.0,
        "C4H10": 0.0,
        "O2": 0.2074,
        "H2O": 0.0099,
        "CO2": 0.0003,
        "N2": 0.7732,
        "Ar": 0.0092,
    }
    ng_comp = {
        "CH4": 0.931,
        "C2H6": 0.0320,
        "C3H8": 0.007,
        "C4H10": 0.004,
        "O2": 0.0,
        "H2O": 0.0,
        "CO2": 0.01,
        "N2": 0.0160,
        "Ar": 0.0,
    }
    _set_port(
        m.fs.air_preheater.tube_inlet, F=3500, T=330, P=1.04e5, comp=air_comp, fix=True
    )
    _set_port(
        m.fs.ng_preheater.tube_inlet, F=280, T=330, P=1.04e5, comp=ng_comp, fix=True
    )
    m.fs.bhx2.tube_inlet.flow_mol.fix(5000)
    m.fs.bhx2.tube_inlet.enth_mol.fix(iapws95.htpx(T=950*pyo.units.K, P=20.6e5*pyo.units.Pa))
    m.fs.bhx2.tube_inlet.pressure.fix(20.6e5)

    m.fs.main_steam_split.split_fraction[:, "h_side"].fix(0.5)
    m.fs.recover_split.split_fraction[:, "h_side"].fix(0.5)

    m.fs.aux_boiler_feed_pump.inlet.flow_mol.fix(5000)
    m.fs.aux_boiler_feed_pump.inlet.enth_mol.fix(
        iapws95.htpx(T=310*pyo.units.K, P=101325*pyo.units.Pa)
    )
    m.fs.aux_boiler_feed_pump.inlet.pressure.fix(101325)
    m.fs.aux_boiler_feed_pump.outlet.pressure.fix(1.1e5)#20.6e5)
    m.fs.aux_boiler_feed_pump.efficiency_isentropic.fix(0.85)

    m.fs.spltf1.split_fraction[:, "out"].fix(0.98)
    m.fs.splta1.split_fraction[:, "out"].fix(0.85)

    m.fs.soec.n_cells.fix(300e6)
    m.fs.soec_heat_duty.fix(0) # going for the theroneutral point here

    m.fs.soec.fc.flow_mol[:, 0].fix(8e-6)
    m.fs.soec.fc.pressure[:, 0].fix(1e5)
    m.fs.soec.fc.mole_frac_comp[:, 0, "H2O"].fix(0.90)
    m.fs.soec.fc.mole_frac_comp[:, 0, "H2"].fix(0.10)

    m.fs.soec.ac.flow_mol[:, 0].fix(1e-5)
    m.fs.soec.ac.pressure[:, 0].fix(1e5)
    m.fs.soec.ac.mole_frac_comp[:, 0, "O2"].fix(0.1)
    m.fs.soec.ac.mole_frac_comp[:, 0, "H2O"].fix(0.9)


def do_initialize(m, solver):
    m.fs.preheat_split.initialize(outlvl=idaeslog.DEBUG)
    iinit.propagate_state(m.fs.fg04)
    iinit.propagate_state(m.fs.fg05)
    m.fs.air_preheater.initialize(outlvl=idaeslog.DEBUG)
    m.fs.ng_preheater.initialize(outlvl=idaeslog.DEBUG)
    iinit.propagate_state(m.fs.ba03)
    iinit.propagate_state(m.fs.bng03)
    m.fs.cmb_mix.initialize(outlvl=idaeslog.DEBUG)
    iinit.propagate_state(m.fs.bng04)
    m.fs.cmb.initialize(outlvl=idaeslog.DEBUG)
    iinit.propagate_state(m.fs.fg01)
    m.fs.bhx2.initialize(outlvl=idaeslog.DEBUG)

    m.fs.aux_boiler_feed_pump.initialize()
    iinit.propagate_state(m.fs.fg02)
    iinit.propagate_state(m.fs.s02)
    m.fs.bhx1.initialize()

    m.fs.soec.initialize()
    iinit.propagate_state(m.fs.h01)
    iinit.propagate_state(m.fs.o01)
    m.fs.spltf1.initialize()
    m.fs.splta1.initialize()

    iinit.propagate_state(m.fs.s03)
    iinit.propagate_state(m.fs.s09)
    m.fs.recover_split.initialize()
    m.fs.main_steam_split.initialize()

    iinit.propagate_state(m.fs.h02)
    iinit.propagate_state(m.fs.s05)
    iinit.propagate_state(m.fs.o02)
    iinit.propagate_state(m.fs.s04)
    m.fs.hxh2.initialize()
    m.fs.hxo2.initialize()

    iinit.propagate_state(m.fs.s06)
    iinit.propagate_state(m.fs.s07)
    m.fs.smix1.initialize()

    iinit.propagate_state(m.fs.s10)
    iinit.propagate_state(m.fs.s11)
    iinit.propagate_state(m.fs.hr01)
    iinit.propagate_state(m.fs.or01)
    m.fs.mxf1.initialize()
    m.fs.mxa1.initialize()


    iinit.propagate_state(m.fs.fg03)
    iinit.propagate_state(m.fs.s08)

    m.fs.bhx2.tube_inlet.unfix()
    m.fs.preheat_split.inlet.unfix()

    m.fs.soec.E_cell.unfix()
    m.fs.air_preheater.tube_inlet.flow_mol.unfix()
    m.fs.ng_preheater.tube_inlet.flow_mol.unfix()

    m.fs.s12_expanded.deactivate()
    m.fs.s13_expanded.deactivate()
    solver.solve(m, tee=True)

    m.fs.s12_expanded.activate()
    m.fs.s13_expanded.activate()

    m.fs.soec.fc.pressure[:, 0].unfix()
    m.fs.soec.fc.temperature[:, 0].unfix()
    m.fs.soec.fc.mole_frac_comp[:, 0, "H2O"].unfix()

    m.fs.soec.ac.pressure[:, 0].unfix()
    m.fs.soec.ac.temperature[:, 0].unfix()
    m.fs.soec.ac.mole_frac_comp[:, 0, "H2O"].unfix()
    m.fs.main_steam_split.split_fraction[:, "h_side"].unfix()
    m.fs.aux_boiler_feed_pump.inlet.flow_mol.unfix()

    m.fs.spltf1.split_fraction[:, "out"].unfix()
    m.fs.splta1.split_fraction[:, "out"].unfix()

    solver.solve(m, tee=True)



def get_solver():
    use_idaes_solver_configuration_defaults()
    idaes.cfg.ipopt["options"]["nlp_scaling_method"] = "user-scaling"
    idaes.cfg.ipopt["options"]["tol"] = 1e-7
    # due to a lot of component mole fractions being on their lower bound of 0
    # bound push result in much longer solve times, so set it low.
    idaes.cfg.ipopt["options"]["bound_push"] = 1e-9
    idaes.cfg.ipopt["options"]["linear_solver"] = "ma27"
    idaes.cfg.ipopt["options"]["max_iter"] = 400
    #idaes.cfg.ipopt["options"]["ma27_pivtol"] = 1e-1
    #idaes.cfg.ipopt["options"]["ma57_pivtol"] = 1e-1
    return pyo.SolverFactory("ipopt")


def tag_inputs_opt_vars(m):
    tags = iutil.ModelTagGroup()
    tags["single_cell_h2_side_inlet_flow"] = iutil.ModelTag(
        expr=m.fs.soec.fc.flow_mol[0, 0],
        format_string="{:.3f}",
        display_units=pyo.units.micromol/pyo.units.s,
        doc="Single cell H2 side inlet flow (feed)",
    )
    tags["single_cell_sweep_flow"] = iutil.ModelTag(
        expr=m.fs.soec.ac.flow_mol[0, 0],
        format_string="{:.3f}",
        display_units=pyo.units.micromol/pyo.units.s,
        doc="Single cell O2 side inlet flow (sweep)",
    )
    tags["feed_h2_frac"] = iutil.ModelTag(
        expr=m.fs.soec.fc.mole_frac_comp[0, 0, "H2"],
        format_string="{:.3f}",
        display_units=None,
        doc="H2 side inlet H2 mole frac (from recycle)",
    )
    tags["sweep_o2_frac"] = iutil.ModelTag(
        expr=m.fs.soec.ac.mole_frac_comp[0, 0, "O2"],
        format_string="{:.3f}",
        display_units=None,
        doc="O2 side inlet O2 mole frac (from recycle)",
    )
    tags["single_cell_heat_required"] = iutil.ModelTag(
        expr=m.fs.soec_heat_duty[0],
        format_string="{:.3f}",
        display_units=pyo.units.W,
        doc="Heat duty of a singel SOEC cell",
    )
    tags["combustor_temperature"] = iutil.ModelTag(
        expr=m.fs.soec_cmb_temperature[0],
        format_string="{:.3f}",
        display_units=pyo.units.K,
        doc="Combustor temperature",
    )
    tags["steam_temperature"] = iutil.ModelTag(
        expr=m.fs.soec_steam_temperature[0],
        format_string="{:.3f}",
        display_units=pyo.units.K,
        doc="Steam temperature (should be same as SOEC)",
    )
    tags["preheat_fg_split_to_air"] = iutil.ModelTag(
        expr=m.fs.preheat_split.split_fraction[0, "air"],
        format_string="{:.3f}",
        display_units=None,
        doc="Split fraction of flue gas to air preheater, rest goes to NG heater",
    )
    tags["recover_split_to_hxh2"] = iutil.ModelTag(
        expr=m.fs.recover_split.split_fraction[0, "h_side"],
        format_string="{:.3f}",
        display_units=None,
        doc="Split fraction of steam to H2 recovery HX, rest goes to O2 HX",
    )
    tags["n_cells"] = iutil.ModelTag(
        expr=m.fs.soec.n_cells,
        format_string="{:,.0f}",
        display_units=None,
        doc="Number of SOEC cells",
    )
    tags["cell_pressure"] = iutil.ModelTag(
        expr=m.fs.aux_boiler_feed_pump.outlet.pressure[0],
        format_string="{:.3f}",
        display_units=pyo.units.kPa,
        doc="Steam and SOEC pressure",
    )
    tags["air_preheater_area"] = iutil.ModelTag(
        expr=m.fs.air_preheater.area,
        format_string="{:.3f}",
        display_units=pyo.units.m**2,
        doc="Air preheater area",
    )
    tags["ng_preheater_area"] = iutil.ModelTag(
        expr=m.fs.ng_preheater.area,
        format_string="{:.3f}",
        display_units=pyo.units.m**2,
        doc="NG preheater area",
    )
    tags["bhx1_area"] = iutil.ModelTag(
        expr=m.fs.bhx1.area,
        format_string="{:.3f}",
        display_units=pyo.units.m**2,
        doc="bhx1 area",
    )
    tags["bhx2_area"] = iutil.ModelTag(
        expr=m.fs.bhx2.area,
        format_string="{:.3f}",
        display_units=pyo.units.m**2,
        doc="bhx2 area",
    )
    tags["hxh2_area"] = iutil.ModelTag(
        expr=m.fs.hxh2.area,
        format_string="{:.3f}",
        display_units=pyo.units.m**2,
        doc="hxh2 area",
    )
    tags["hxo2_area"] = iutil.ModelTag(
        expr=m.fs.hxo2.area,
        format_string="{:.3f}",
        display_units=pyo.units.m**2,
        doc="hxo2 area",
    )
    m.tag_input = tags
    display_input_tags(m)


def tag_for_pfd_and_tables(m):
    tag_group = iutil.ModelTagGroup()
    stream_states = tables.stream_states_dict(
        tables.arcs_to_stream_dict(
            m.fs,
            descend_into=False,
            additional={
                "fg06":m.fs.ng_preheater.shell_outlet,
                "fg07":m.fs.air_preheater.shell_outlet,
                "s01":m.fs.aux_boiler_feed_pump.inlet,
                "h03":m.fs.hxh2.shell_outlet,
                "o03":m.fs.hxo2.shell_outlet,
            },
        )
    )
    for i, s in stream_states.items():  # create the tags for steam quantities
        tag_group[f"{i}_Fmol"] = iutil.ModelTag(
            expr=s.flow_mol,
            format_string="{:.3f}",
            display_units=pyo.units.kmol/pyo.units.s
        )
        tag_group[f"{i}_Fmass"] = iutil.ModelTag(
            expr=s.flow_mass,
            format_string="{:.3f}",
            display_units=pyo.units.kg/pyo.units.s
        )
        tag_group[f"{i}_P"] = iutil.ModelTag(
            expr=s.pressure,
            format_string="{:.1f}",
            display_units=pyo.units.kPa
        )
        tag_group[f"{i}_T"] = iutil.ModelTag(
            expr=s.temperature,
            format_string="{:.2f}",
            display_units=pyo.units.K
        )
        try:
            tag_group[f"{i}_vf"] = iutil.ModelTag(
                expr=s.phase_frac["Vap"],
                format_string="{:.3f}",
                display_units=None
        )
        except (KeyError, AttributeError):
            pass
        try:
            for c in s.mole_frac_comp:
                tag_group[f"{i}_y{c}"] = iutil.ModelTag(
                    expr=s.mole_frac_comp[c]*100,
                    format_string="{:.3f}",
                    display_units="%"
                )
        except (KeyError, AttributeError):
            pass
        try:
            tag_group[f"{i}_y"] = iutil.ModelTag(
                expr=s.mole_frac_comp,
                format_string="{:.3f}",
                display_units=None
            )
        except (KeyError, AttributeError):
            pass

    tag_group["soec_power"] = iutil.ModelTag(
        expr=m.fs.soec.total_power[0],
        format_string="{:.2f}",
        display_units=pyo.units.MW
    )
    tag_group["soec_n_cells"] = iutil.ModelTag(
        expr=m.fs.soec.n_cells,
        format_string="{:,.0f}",
        display_units=None
    )
    tag_group["E_cell"] = iutil.ModelTag(
        expr=m.fs.soec.E_cell[0],
        format_string="{:.4f}",
        display_units=pyo.units.V
    )

    m.tag_pfd = tag_group


def display_input_tags(m):
    # this is special for this model.  The input tags are not indexed
    print("")
    for key, tag in m.tag_input.items():
        print(key)
        print(f"    {tag.doc}")
        print(f"    display units: {tag._display_units}")
        print(f"    native units: {pyo.units.get_units(tag.expression)}")
        print(f"    value {tag}, fixed: {tag.expression.fixed}")
        print("")



def check_scaling(m):
    jac, nlp = iscale.get_jacobian(m, scaled=True)
    print("Extreme Jacobian entries:")
    for i in iscale.extreme_jacobian_entries(jac=jac, nlp=nlp, large=100):
        print(f"    {i[0]:.2e}, [{i[1]}, {i[2]}]")
    print("Unscaled constraints:")
    for c in iscale.unscaled_constraints_generator(m):
        print(f"    {c}")
    print("Scaled constraints by factor:")
    for c, s in iscale.constraints_with_scale_factor_generator(m):
        print(f"    {c}, {s}")
    print("Badly scaled variables:")
    for v, sv in iscale.badly_scaled_var_generator(
        m, large=1e2, small=1e-2, zero=1e-12
    ):
        print(f"    {v} -- {sv} -- {iscale.get_scaling_factor(v)}")
    print(f"Jacobian Condition Number: {iscale.jacobian_cond(jac=jac):.2e}")


def write_pfd_results(m, filename, infilename=None):
    """
    Write simulation results in a template PFD in svg format and save as
    filename.

    Args:
        filename: (str) file namd for output
        tags: (dict) tag keys and expression values
        tag_format: (dict) tag keys and format string values
        infilename: input file name, if you want to use an alternative diagram

    Returns:
        None
    """
    if infilename is None:
        infilename = os.path.join(this_file_dir(), "soec.svg")
    with open(infilename, "r") as f:
        iutil.svg_tag(svg=f, tag_group=m.tag_pfd, outfile=filename)


def get_model(m=None, name="SOEC Module"):
    m = add_flowsheet(m)
    add_properties(m)
    add_preheater(m)
    add_combustor(m)
    add_aux_boiler_steam(m)
    add_soec_unit(m)
    add_recovery_hx(m)
    add_more_hx_connections(m)
    add_soec_inlet_mix(m)
    add_constraints(m)
    expand_arcs = pyo.TransformationFactory("network.expand_arcs")
    expand_arcs.apply_to(m)

    tag_inputs_opt_vars(m)
    tag_for_pfd_and_tables(m)
    set_guess(m)
    set_inputs(m)
    solver = get_solver()
    do_initialize(m, solver)
    display_input_tags(m)
    return m



if __name__ == "__main__":
    m = get_model()
    write_pfd_results(m, "soec_init.svg")
    check_scaling(m)
