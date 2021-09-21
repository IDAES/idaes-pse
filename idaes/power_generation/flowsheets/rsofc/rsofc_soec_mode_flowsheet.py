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
from idaes.core.util.model_statistics import degrees_of_freedom, variables_near_bounds_generator


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
from idaes.generic_models.unit_models.pressure_changer import \
    ThermodynamicAssumption
from idaes.generic_models.unit_models.separator import SplittingType

from idaes.generic_models.properties.helmholtz.helmholtz import (
    HelmholtzThermoExpressions as ThermoExpr,
)
import idaes.logger as idaeslog


# Check algebraic jacobian
from pyomo.util.subsystems import (
        TemporarySubsystemManager,
        )
from pyomo.contrib.incidence_analysis import (
        IncidenceGraphInterface,
        solve_strongly_connected_components,
        )
from pyomo.common.collections import ComponentSet, ComponentMap
from pyomo.core.expr.visitor import identify_variables

# Imports necessary to check DAE index
from pyomo.dae.flatten import flatten_dae_components

# sys.path.append("C:/cokoli/IDAES_2/workspace/idaes-2/research")
# Currently depedenent on Robby's idaes-pse fork, capr-categorize branch
from idaes.apps.caprese.categorize import (
        categorize_dae_variables_and_constraints,
        VariableCategory,
        ConstraintCategory,
        )
VC = VariableCategory
CC = ConstraintCategory
from pyomo.contrib.incidence_analysis import (
        IncidenceGraphInterface,
        maximum_matching,
        )
from pyomo.contrib.incidence_analysis.interface import (
        get_structural_incidence_matrix,
        )
from pyomo.contrib.incidence_analysis.dulmage_mendelsohn import (
        dulmage_mendelsohn,
        )
from scipy.sparse.linalg import splu


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
    fuel_comp = {  # components present
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
    air_comp = {  # components present
        "O2",
        "H2O",
        "CO2",
        "N2",
        "Ar",
    }
    m.fs.fg_prop = GenericParameterBlock(
        default=get_prop(components=fuel_comp, phases=["Vap"], eos=EosType.IDEAL)
    )
    m.fs.fg_prop.set_default_scaling("mole_frac_comp", 1e2)
    m.fs.fg_prop.set_default_scaling("mole_frac_phase_comp", 1e2)
    # m.fs.fg_prop.set_default_scaling("flow_mol", 1)
    # m.fs.fg_prop.set_default_scaling("flow_mol_phase", 1)
    # m.fs.fg_prop.set_default_scaling("temperature", 1e-3)
    # m.fs.fg_prop.set_default_scaling("pressure", 1e-2)
    # m.fs.fg_prop.set_default_scaling("enth_mol_phase", 1e-4)
    # m.fs.fg_prop.set_default_scaling("entr_mol_phase", 1e-2)    
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
    m.fs.air_prop = GenericParameterBlock(
        default=get_prop(components=air_comp, phases=["Vap"], eos=EosType.IDEAL)
    )
    m.fs.air_prop.set_default_scaling("mole_frac_comp", 1e2)
    m.fs.air_prop.set_default_scaling("mole_frac_phase_comp", 1e2)
    # m.fs.air_prop.set_default_scaling("flow_mol", 1)
    # m.fs.air_prop.set_default_scaling("flow_mol_phase", 1)
    # m.fs.air_prop.set_default_scaling("temperature", 1e-3)
    # m.fs.air_prop.set_default_scaling("pressure", 1e-2)
    # m.fs.air_prop.set_default_scaling("enth_mol_phase", 1e-4)
    # m.fs.air_prop.set_default_scaling("entr_mol_phase", 1e-2)

def add_asu(m):   
    # build ASU, oxycombustor
    m.fs.air_compressor_s1 = gum.PressureChanger(
        default={"compressor": True,
                 "property_package": m.fs.air_prop,
                 "thermodynamic_assumption":
                     ThermodynamicAssumption.isentropic})

    m.fs.intercooler_s1 = gum.Heater(
        default={"property_package": m.fs.air_prop,
                 "has_pressure_change": True})

    m.fs.air_compressor_s2 = gum.PressureChanger(
        default={"compressor": True,
                 "property_package": m.fs.air_prop,
                 "thermodynamic_assumption":
                     ThermodynamicAssumption.isentropic})

    m.fs.intercooler_s2 = gum.Heater(
        default={"property_package": m.fs.air_prop,
                 "has_pressure_change": True})

    m.fs.ASU = gum.Separator(
        default={"outlet_list": ["N2_outlet", "O2_outlet"],
                 "split_basis": SplittingType.componentFlow,
                 "property_package": m.fs.air_prop})

    # m.fs.ASU_O2_outlet = gum.Heater(
    #     default={"has_pressure_change": True,
    #              "property_package": m.fs.fg_prop})

    # arcs for ASU, oxycombustor and CPU
    m.fs.STAGE_1_OUT = Arc(
        source=m.fs.air_compressor_s1.outlet,
        destination=m.fs.intercooler_s1.inlet)

    m.fs.IC_1_OUT = Arc(
        source=m.fs.intercooler_s1.outlet,
        destination=m.fs.air_compressor_s2.inlet)

    m.fs.STAGE_2_OUT = Arc(
        source=m.fs.air_compressor_s2.outlet,
        destination=m.fs.intercooler_s2.inlet)

    m.fs.TO_ASU = Arc(
        source=m.fs.intercooler_s2.outlet,
        destination=m.fs.ASU.inlet)


def add_preheater(m):
    m.fs.air_preheater = gum.HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_underwood_callback,
            "shell": {"property_package": m.fs.air_prop},
            "tube": {"property_package": m.fs.air_prop},
        }
    )
    m.fs.ng_preheater = gum.HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_underwood_callback,
            "shell": {"property_package": m.fs.air_prop},
            "tube": {"property_package": m.fs.fg_prop},
        }
    )
    m.fs.preheat_split = gum.Separator(
        default={
            "property_package": m.fs.air_prop,
            "outlet_list": ["air", "ng"],
        }
    )
    m.fs.fg04 = Arc(
        source=m.fs.preheat_split.ng,
        destination=m.fs.ng_preheater.shell_inlet
    )
    m.fs.fg05 = Arc(
        source=m.fs.preheat_split.air,
        destination=m.fs.air_preheater.shell_inlet
    )
    m.fs.ba02 = Arc(
        source=m.fs.ASU.O2_outlet,
        destination=m.fs.air_preheater.tube_inlet
        )

def add_combustor(m):
    m.fs.pre_oxycombustor_translator = gum.Translator(
        default={"outlet_state_defined": True,
                  "inlet_property_package": m.fs.air_prop,
                  "outlet_property_package": m.fs.fg_prop})

    # Additional constraints to specify the translator block
    @m.fs.pre_oxycombustor_translator.Constraint(m.fs.time)
    def pre_oxycombustor_translator_F(b, t):
        return b.inlet.flow_mol[t] == b.outlet.flow_mol[t]

    @m.fs.pre_oxycombustor_translator.Constraint(m.fs.time)
    def pre_oxycombustor_translator_T(b, t):
        return b.inlet.temperature[t] == b.outlet.temperature[t]

    @m.fs.pre_oxycombustor_translator.Constraint(m.fs.time)
    def pre_oxycombustor_translator_P(b, t):
        return b.inlet.pressure[t] == b.outlet.pressure[t]

    @m.fs.pre_oxycombustor_translator.Constraint(m.fs.time,
                                             m.fs.air_prop.component_list)
    def pre_oxycombustor_translator_x(b, t, j):
        return b.inlet.mole_frac_comp[t, j] == b.outlet.mole_frac_comp[t, j]

    for j in m.fs.fg_prop.component_list:
        if j not in m.fs.air_prop.component_list:
            m.fs.pre_oxycombustor_translator.outlet.mole_frac_comp[0, j].fix(0)

    m.fs.cmb_mix = gum.Mixer(default={
        "property_package": m.fs.fg_prop,
        "inlet_list":["ng", "air"],
        "momentum_mixing_type":gum.MomentumMixingType.none})
    m.fs.cmb = gum.StoichiometricReactor(default={
        "property_package": m.fs.fg_prop,
        "reaction_package": m.fs.fg_combust,
        "has_pressure_change": False,
        # "has_heat_of_reaction": True
        })
    @m.fs.cmb_mix.Constraint(m.fs.time)
    def pressure_eqn(b, t):
        return b.mixed_state[t].pressure == b.air_state[t].pressure
    @m.fs.cmb.Constraint(m.fs.time, m.rxns.keys())
    def reaction_extent(b, t, r):
        k = m.rxns[r]
        prp = b.control_volume.properties_in[t]
        stc = -m.fs.fg_combust.rate_reaction_stoichiometry[r, "Vap", k]
        extent = b.rate_reaction_extent[t, r]
        return extent*stc == prp.flow_mol*prp.mole_frac_comp[k]

    m.fs.post_oxycombustor_translator = gum.Translator(
        default={"outlet_state_defined": True,
                  "inlet_property_package": m.fs.fg_prop,
                  "outlet_property_package": m.fs.air_prop})

    # Additional constraints to specify the translator block
    @m.fs.post_oxycombustor_translator.Constraint(m.fs.time)
    def post_oxycombustor_translator_F(b, t):
        return b.inlet.flow_mol[t] == b.outlet.flow_mol[t]

    @m.fs.post_oxycombustor_translator.Constraint(m.fs.time)
    def post_oxycombustor_translator_T(b, t):
        return b.inlet.temperature[t] == b.outlet.temperature[t]

    @m.fs.post_oxycombustor_translator.Constraint(m.fs.time)
    def post_oxycombustor_translator_P(b, t):
        return b.inlet.pressure[t] == b.outlet.pressure[t]

    @m.fs.post_oxycombustor_translator.Constraint(m.fs.time,
                                             m.fs.air_prop.component_list)
    def post_oxycombustor_translator_x(b, t, j):
        return b.inlet.mole_frac_comp[t, j] == b.outlet.mole_frac_comp[t, j]

    # for j in m.fs.fg_prop.component_list:
    #     if j not in m.fs.air_prop.component_list:
    #         m.fs.post_oxycombustor_translator.outlet.mole_frac_comp[0, j].fix(0)
    m.fs.cmb_mix_in = Arc(  # TODO - rename arc
        source=m.fs.air_preheater.tube_outlet,
        destination=m.fs.pre_oxycombustor_translator.inlet
    )
    m.fs.ba03 = Arc(
        source=m.fs.pre_oxycombustor_translator.outlet,
        destination=m.fs.cmb_mix.air
    )
    m.fs.bng03 = Arc(
        source=m.fs.ng_preheater.tube_outlet, destination=m.fs.cmb_mix.ng
    )
    m.fs.bng04 = Arc(
        source=m.fs.cmb_mix.outlet, destination=m.fs.cmb.inlet
    )
    m.fs.cmb_out = Arc(  # TODO - rename arc
        source=m.fs.cmb.outlet,
        destination=m.fs.post_oxycombustor_translator.inlet
    )

def add_aux_boiler_steam(m):
    m.fs.bhx2 = gum.HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_underwood_callback,
            "shell": {"property_package": m.fs.air_prop},
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
            "shell": {"property_package": m.fs.air_prop},
            "tube": {"property_package": m.fs.water_prop},
        }
    )
    m.fs.recover_split = HelmSplitter(
        default={
            "property_package": m.fs.water_prop,
            "outlet_list": ["h_side", "o_side"]
        }
    )
    m.fs.fg01 = Arc(source=m.fs.post_oxycombustor_translator.outlet,
                    destination=m.fs.bhx2.shell_inlet)
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
    # m.fs.hxo2 = gum.HeatExchanger(
    #     default={
    #         "delta_temperature_callback": delta_temperature_underwood_callback,
    #         "shell": {"property_package": m.fs.air_prop},
    #         "tube": {"property_package": m.fs.water_prop},
    #     }
    # )
    # m.fs.smix1 = HelmMixer(
    #     default={
    #         "momentum_mixing_type": MomentumMixingType.none,
    #         "inlet_list": ["hxh2", "hxo2"],
    #         "property_package": m.fs.water_prop,
    #     }
    # )
    m.fs.h02 = Arc(source=m.fs.spltf1.out, destination=m.fs.hxh2.shell_inlet)
    m.fs.s05 = Arc(source=m.fs.recover_split.h_side, destination=m.fs.hxh2.tube_inlet)
    # m.fs.o02 = Arc(source=m.fs.splta1.out, destination=m.fs.hxo2.shell_inlet)
    # m.fs.s04 = Arc(source=m.fs.recover_split.o_side, destination=m.fs.hxo2.tube_inlet)
    # m.fs.s06 = Arc(source=m.fs.hxo2.tube_outlet, destination=m.fs.smix1.hxo2)
    # m.fs.s07 = Arc(source=m.fs.hxh2.tube_outlet, destination=m.fs.smix1.hxh2)

    # @m.fs.smix1.Constraint(m.fs.time)
    # def pressure_eqn(b, t):
    #     return b.mixed_state[t].pressure == b.hxh2_state[t].pressure


def add_soec_unit(m):
    m.fs.soec = pum.IsothermalSofc(
        default={
            "nz": 20,
            "nxfe": 10,
            "nxae": 10,
            "soec": True,
            "air_side_comp_list": ["O2", "H2O",
                                   "CO2", "N2", "Ar"],
            "fuel_side_comp_list": ["H2O", "H2"],
            "air_side_stoich": {"O2": -0.25, "H2O": 0,
                                "CO2": 0, "N2": 0, "Ar": 0},
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
            "property_package": m.fs.air_prop,
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
            "property_package": m.fs.air_prop,
            "inlet_list": ["air", "recycle"],
            "momentum_mixing_type": gum.MomentumMixingType.none,
        }
    )

    @m.fs.mxf1.Constraint(m.fs.time)
    def fmxpress_eqn(b, t):
        return b.mixed_state[t].pressure == b.water_state[t].pressure

    # TODO - removed as air stream now goes to mxa1
    # @m.fs.mxa1.Constraint(m.fs.time)
    # def amxpress_eqn(b, t):
    #     return b.mixed_state[t].pressure == b.water_state[t].pressure

    # Add ports to connet pure steam to steam + h2 or steam + o2
    m.fs.main_steam_split._temperature_h_side_ref = pyo.Reference(
        m.fs.main_steam_split.h_side_state[:].temperature
    )
    # TODO - removed as air stream now goes to mxa1
    # m.fs.main_steam_split._temperature_o_side_ref = pyo.Reference(
    #     m.fs.main_steam_split.o_side_state[:].temperature
    # )

    @m.fs.main_steam_split.Expression(m.fs.time, m.fs.soec.fc.config.comp_list)
    def h_side_mole_frac_expr(b, t, i):
        if i == "H2O":
            return 1
        else:
            return 0

    # TODO - removed as air stream now goes to mxa1
    # @m.fs.main_steam_split.Expression(m.fs.time, m.fs.soec.ac.config.comp_list)
    # def o_side_mole_frac_expr(b, t, i):
    #     if i == "H2O":
    #         return 1
    #     else:
    #         return 0

    m.fs.main_steam_split.h_side_adapt = Port(
        rule=lambda b: {
            "flow_mol": m.fs.main_steam_split._flow_mol_h_side_ref,
            "pressure": m.fs.main_steam_split._pressure_h_side_ref,
            "temperature": m.fs.main_steam_split._temperature_h_side_ref,
            "mole_frac_comp": m.fs.main_steam_split.h_side_mole_frac_expr,
        }
    )

    # TODO - removed as air stream now goes to mxa1
    # m.fs.main_steam_split.o_side_adapt = Port(
    #     rule=lambda b: {
    #         "flow_mol": m.fs.main_steam_split._flow_mol_o_side_ref,
    #         "pressure": m.fs.main_steam_split._pressure_o_side_ref,
    #         "temperature": m.fs.main_steam_split._temperature_o_side_ref,
    #         "mole_frac_comp": m.fs.main_steam_split.o_side_mole_frac_expr,
    #     }
    # )

    m.fs.s10 = Arc(
        source=m.fs.main_steam_split.h_side_adapt,
        destination=m.fs.mxf1.water,
    )
    # TODO - removed as air stream now goes to mxa1
    # m.fs.s11 = Arc(
    #     source=m.fs.main_steam_split.o_side_adapt,
    #     destination=m.fs.mxa1.water,
    # )
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
        source=m.fs.hxh2.tube_outlet,
        destination=m.fs.bhx2.tube_inlet
    )


def add_constraints(m):
    m.fs.soec_heat_duty = pyo.Var(m.fs.time, units=pyo.units.W)
    @m.fs.Constraint(m.fs.time)
    def heat_duty_soec_zero_eqn(b, t):
        return b.soec.heat_duty[t] == b.soec_heat_duty[t]

    # TODO - why is this eqn relevant? cmb outlet temperature is not 2,000
    m.fs.soec_cmb_temperature = pyo.Var(m.fs.time, initialize=2000, units=pyo.units.K)
    @m.fs.Constraint(m.fs.time)
    def soec_cmb_temperature_eqn(b, t):
        return m.fs.cmb.outlet.temperature[t] == m.fs.soec_cmb_temperature[t]

    m.fs.soec_steam_temperature = pyo.Var(m.fs.time, initialize=1073.15, units=pyo.units.K)
    @m.fs.Constraint(m.fs.time)
    def soec_steam_temperature_eqn(b, t):
        return m.fs.bhx2.tube.properties_out[t].temperature == m.fs.soec_steam_temperature[t]

def set_guess(m):
    fg_comp_guess = { # air_prop is used as assumption is all fuel is combusted
        "O2": 0.049,
        "H2O": 0.2,
        "CO2": 0.2,
        "N2": 0.55,
        "Ar": 0.001,
    }
    _set_port(
        m.fs.preheat_split.inlet, F=1650, T=550, P=1.04e5,
        comp=fg_comp_guess, fix=True
    )


def set_inputs(m):
    m.fs.soec_cmb_temperature.fix(2000)
    m.fs.soec_steam_temperature.fix(1073.15)
    # TODO - reduced as choice of the air_preheater area affects convergence of air_preheater
    m.fs.air_preheater.area.fix(2000)
    m.fs.air_preheater.overall_heat_transfer_coefficient.fix(100)
    m.fs.ng_preheater.area.fix(300)
    m.fs.ng_preheater.overall_heat_transfer_coefficient.fix(100)
    # TODO - reduced as choice of the bxh2 area affects convergence of bxh2
    m.fs.bhx2.area.fix(4000)
    m.fs.bhx2.overall_heat_transfer_coefficient.fix(100)
    m.fs.bhx1.area.fix(500)
    m.fs.bhx1.overall_heat_transfer_coefficient.fix(100)
    m.fs.hxh2.area.fix(4000)
    m.fs.hxh2.overall_heat_transfer_coefficient.fix(100)
    # m.fs.hxo2.area.fix(4000)
    # m.fs.hxo2.overall_heat_transfer_coefficient.fix(100)


    m.fs.preheat_split.split_fraction[:, "air"].fix(0.9)

    air_comp = {
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
        "O2": 1e-5,
        "H2O": 1e-5,
        "CO2": 0.01,
        "N2": 0.0160,
        "Ar": 1e-5,
    }
    _set_port( # TODO - relate NG feed to air feed requirement
        m.fs.air_compressor_s1.inlet, F=3500, T=330, P=1.04e5, comp=air_comp, fix=True
    )
    _set_port(
        m.fs.ng_preheater.tube_inlet, F=280, T=330, P=1.04e5, comp=ng_comp, fix=True
    )
    _set_port(  # TODO - increased flow to 10K to avoid log eval error in solve step 2
        m.fs.mxa1.air, F=10000, T=1073.15, P=1.11e5,
        comp={"O2": 0.2074, "H2O": 0.0099, "CO2": 0.0003,
              "N2": 0.7732, "Ar": 0.0092,}, fix=True
    )

    m.fs.bhx2.tube_inlet.flow_mol.fix(5000)
    # TODO - should this pressure be set to 20.6e5 Pa. Why not the expected 1.1e5
    m.fs.bhx2.tube_inlet.enth_mol.fix(iapws95.htpx(T=950*pyo.units.K, P=20.6e5*pyo.units.Pa))
    m.fs.bhx2.tube_inlet.pressure.fix(20.6e5)

    # TODO - Set to 0.9999 as air stream now goes to mxa1
    m.fs.main_steam_split.split_fraction[:, "h_side"].fix(0.9999)
    m.fs.recover_split.split_fraction[:, "h_side"].fix(0.9999)

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
    # TODO - fixed N2, Ar and CO2 for air feed
    m.fs.soec.ac.mole_frac_comp[:, 0, "N2"].fix(0)
    m.fs.soec.ac.mole_frac_comp[:, 0, "Ar"].fix(0)
    m.fs.soec.ac.mole_frac_comp[:, 0, "CO2"].fix(0)

    # air compressors and intercoolers
    m.fs.air_compressor_s1.outlet.pressure.fix(111422)  # Pa (34 psia)
    m.fs.air_compressor_s1.efficiency_isentropic.fix(0.84)
    m.fs.intercooler_s1.outlet.temperature.fix(310.93)  # K (100 F)
    m.fs.intercooler_s1.deltaP.fix(-3447)  # Pa (-0.5 psi)
    m.fs.air_compressor_s2.outlet.pressure.fix(130686)  # Pa (79 psia)
    m.fs.air_compressor_s2.efficiency_isentropic.fix(0.84)
    m.fs.intercooler_s2.outlet.temperature.fix(310.93)  # K (100 F)
    m.fs.intercooler_s2.deltaP.fix(-3447)  # Pa (-0.5 psi)

    # air seperation unit
    m.fs.ASU.split_fraction[0, "O2_outlet", "CO2"].fix(0)
    m.fs.ASU.split_fraction[0, "O2_outlet", "H2O"].fix(0)
    m.fs.ASU.split_fraction[0, "O2_outlet", "N2"].fix(0.0005)
    m.fs.ASU.split_fraction[0, "O2_outlet", "O2"].fix(0.9691)
    m.fs.ASU.split_fraction[0, "O2_outlet", "Ar"].fix(0.0673)


def do_initialize(m, solver):

    m.fs.preheat_split.initialize(outlvl=idaeslog.DEBUG)
    iinit.propagate_state(m.fs.fg04)
    iinit.propagate_state(m.fs.fg05)

    m.fs.air_compressor_s1.initialize(outlvl=idaeslog.DEBUG)
    iinit.propagate_state(m.fs.STAGE_1_OUT)
    m.fs.intercooler_s1.initialize(outlvl=idaeslog.DEBUG)
    iinit.propagate_state(m.fs.IC_1_OUT)
    m.fs.air_compressor_s2.initialize(outlvl=idaeslog.DEBUG)
    iinit.propagate_state(m.fs.STAGE_2_OUT)
    m.fs.intercooler_s2.initialize(outlvl=idaeslog.DEBUG)
    iinit.propagate_state(m.fs.TO_ASU)
    m.fs.ASU.initialize(outlvl=idaeslog.DEBUG)
    iinit.propagate_state(m.fs.ba02)
    m.fs.air_preheater.initialize(outlvl=idaeslog.DEBUG)
    iinit.propagate_state(m.fs.ba03)
    m.fs.pre_oxycombustor_translator.initialize()
    iinit.propagate_state(m.fs.cmb_mix_in)
    


    m.fs.ng_preheater.initialize(outlvl=idaeslog.DEBUG)
    iinit.propagate_state(m.fs.bng03)
    m.fs.cmb_mix.initialize(outlvl=idaeslog.DEBUG)
    iinit.propagate_state(m.fs.bng04)
    m.fs.cmb.initialize(outlvl=idaeslog.DEBUG)
    iinit.propagate_state(m.fs.cmb_out)
    # residual_checker(m.fs.cmb)
    # m.fs.cmb.inlet.display()
    # m.fs.cmb.outlet.display()
    # print('display cmb')
    # m.fs.cmb.display()
    # print('dof = ', degrees_of_freedom(m.fs.cmb))   # check_scaling(m.fs.cmb)
    # import sys
    # sys.exit('stop')
    m.fs.pre_oxycombustor_translator.initialize()
    iinit.propagate_state(m.fs.fg01)
    m.fs.bhx2.initialize(outlvl=idaeslog.DEBUG)

    m.fs.aux_boiler_feed_pump.initialize()
    iinit.propagate_state(m.fs.fg02)
    iinit.propagate_state(m.fs.s02)
    m.fs.bhx1.initialize()
    # m.fs.bhx2.tube_outlet.display()
    # import sys
    # sys.exit('stop')
    m.fs.soec.initialize()
    iinit.propagate_state(m.fs.h01)
    iinit.propagate_state(m.fs.o01)
    m.fs.spltf1.initialize()
    m.fs.splta1.initialize()
    # import sys
    # sys.exit('stop')
    iinit.propagate_state(m.fs.s03)
    iinit.propagate_state(m.fs.s09)
    m.fs.recover_split.initialize()
    m.fs.main_steam_split.initialize()

    iinit.propagate_state(m.fs.h02)
    iinit.propagate_state(m.fs.s05)
    # iinit.propagate_state(m.fs.o02)
    # iinit.propagate_state(m.fs.s04)
    m.fs.hxh2.initialize()
    # m.fs.hxo2.initialize()

    # iinit.propagate_state(m.fs.s06)
    # iinit.propagate_state(m.fs.s07)
    # m.fs.smix1.initialize()

    iinit.propagate_state(m.fs.s10)
    # TODO - removed as air stream now goes to mxa1
    # iinit.propagate_state(m.fs.s11)
    iinit.propagate_state(m.fs.hr01)
    iinit.propagate_state(m.fs.or01)
    m.fs.mxf1.initialize()
    m.fs.mxa1.initialize()


    iinit.propagate_state(m.fs.fg03)
    iinit.propagate_state(m.fs.s08)

    m.fs.bhx2.tube_inlet.unfix()
    m.fs.preheat_split.inlet.unfix()
    
    # TODO - air_preheater inlet now replaced with air_compressor_s1 inlet
    # m.fs.air_preheater.tube_inlet.flow_mol.unfix()
    m.fs.air_compressor_s1.inlet.flow_mol.unfix()

    m.fs.ng_preheater.tube_inlet.flow_mol.unfix()
    m.fs.soec.E_cell.unfix()
    m.fs.s12_expanded.deactivate()
    m.fs.s13_expanded.deactivate()

    iscale.calculate_scaling_factors(m)    
    # TODO - this solve isn't robust
    # strip_bounds=True # strip bounds
    # print('dof = ', degrees_of_freedom(m))
    solver.solve(m, tee=True,
                  options={"max_iter":400}
                  )
    
    # check_scaling(m)
    # residual_checker(m)
    # import sys
    # sys.exit('stop')

    m.fs.s12_expanded.activate()
    m.fs.s13_expanded.activate()

    m.fs.soec.fc.pressure[:, 0].unfix()
    m.fs.soec.fc.temperature[:, 0].unfix()
    m.fs.soec.fc.mole_frac_comp[:, 0, "H2O"].unfix()

    m.fs.soec.ac.pressure[:, 0].unfix()
    m.fs.soec.ac.temperature[:, 0].unfix()
    m.fs.soec.ac.mole_frac_comp[:, 0, "H2O"].unfix()
    # TODO - unfixed N2, Ar and CO2 for air feed
    m.fs.soec.ac.mole_frac_comp[:, 0, "N2"].unfix()
    m.fs.soec.ac.mole_frac_comp[:, 0, "Ar"].unfix()
    m.fs.soec.ac.mole_frac_comp[:, 0, "CO2"].unfix()
    m.fs.aux_boiler_feed_pump.inlet.flow_mol.unfix()
    # TODO - Removed as air stream now goes to mxa1
    # m.fs.main_steam_split.split_fraction[:, "h_side"].unfix()
    

    m.fs.spltf1.split_fraction[:, "out"].unfix()
    m.fs.splta1.split_fraction[:, "out"].unfix()

    solver.solve(m, tee=True, symbolic_solver_labels=True)


def additional_scaling(m):
    # for t, c in m.fs.smix1.pressure_eqn.items():
    #     iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in m.fs.air_preheater.heat_transfer_equation.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in m.fs.ng_preheater.heat_transfer_equation.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in m.fs.bhx2.heat_transfer_equation.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in m.fs.bhx1.heat_transfer_equation.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in m.fs.hxh2.heat_transfer_equation.items():
        iscale.constraint_scaling_transform(c, 1e-6, overwrite=True)
    # for t, c in m.fs.hxo2.heat_transfer_equation.items():
    #     iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in m.fs.mxf1.enthalpy_mixing_equations.items():
        iscale.constraint_scaling_transform(c, 1e-6, overwrite=True)
    for t, c in m.fs.mxa1.enthalpy_mixing_equations.items():
        iscale.constraint_scaling_transform(c, 1e-6, overwrite=True)
    for t, c in m.fs.air_preheater.unit_heat_balance.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in m.fs.ng_preheater.unit_heat_balance.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in m.fs.bhx2.unit_heat_balance.items():
        iscale.constraint_scaling_transform(c, 1e-6, overwrite=True)
    for t, c in m.fs.bhx1.unit_heat_balance.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in m.fs.hxh2.unit_heat_balance.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    # for t, c in m.fs.hxo2.unit_heat_balance.items():
    #     iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in m.fs.mxf1.fmxpress_eqn.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    # for t, c in m.fs.mxa1.amxpress_eqn.items():
    #     iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in m.fs.cmb_mix.pressure_eqn.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)

    for t, c in m.fs.cmb_mix.enthalpy_mixing_equations.items():
        iscale.constraint_scaling_transform(c, 1e-6, overwrite=True)
    for t, c in m.fs.cmb.control_volume.enthalpy_balances.items():
        iscale.constraint_scaling_transform(c, 1e-6, overwrite=True)
    for t, c in m.fs.cmb.control_volume.properties_in[t].total_flow_balance.items():
        iscale.constraint_scaling_transform(c, 1e-1, overwrite=True)
    for t, c in m.fs.cmb.control_volume.properties_in[0].component_flow_balances.items():
        iscale.constraint_scaling_transform(c, 1e3, overwrite=True)
    for (t, j), c in m.fs.cmb.control_volume.material_balances.items():
        iscale.constraint_scaling_transform(c, 1e-1, overwrite=True)
    for (t, p, j), c in m.fs.cmb.control_volume.rate_reaction_stoichiometry_constraint.items():
        iscale.constraint_scaling_transform(c, 1e-1, overwrite=True) 

    for t, c in m.fs.pre_oxycombustor_translator.pre_oxycombustor_translator_F.items():
        iscale.constraint_scaling_transform(c, 1e-1, overwrite=True)
    for t, c in m.fs.pre_oxycombustor_translator.pre_oxycombustor_translator_P.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)

    for t, c in m.fs.post_oxycombustor_translator.post_oxycombustor_translator_F.items():
        iscale.constraint_scaling_transform(c, 1e-1, overwrite=True)
    for t, c in m.fs.post_oxycombustor_translator.post_oxycombustor_translator_P.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    # for t, c in m.fs.air_compressor_s1.isentropic_energy_balance.items():
    #     iscale.constraint_scaling_transform(c, 1e-1, overwrite=True)
    # for t, c in m.fs.air_compressor_s1.control_volume.enthalpy_balances.items():
    #     iscale.constraint_scaling_transform(c, 1e-1, overwrite=True)
    # for t, c in m.fs.air_compressor_s2.isentropic_energy_balance.items():
    #     iscale.constraint_scaling_transform(c, 1e-1, overwrite=True)
    # for t, c in m.fs.air_compressor_s2.control_volume.enthalpy_balances.items():
    #     iscale.constraint_scaling_transform(c, 1e-1, overwrite=True)
    # for t, c in m.fs.intercooler_s1.control_volume.enthalpy_balances.items():
    #     iscale.constraint_scaling_transform(c, 1e-1, overwrite=True)
    # for t, c in m.fs.intercooler_s2.control_volume.enthalpy_balances.items():
    #     iscale.constraint_scaling_transform(c, 1e-1, overwrite=True)


    # for t, v in m.fs.aux_boiler_feed_pump.properties_out[t].enth_mol.items():
    #     iscale.set_scaling_factor(m.fs.aux_boiler_feed_pump.properties_out[t].enth_mol, 1e3) 

    for t, c in m.fs.aux_boiler_feed_pump.eq_work.items():
        iscale.constraint_scaling_transform(c, 1e-10, overwrite=True)
    # for (t, r), c in m.fs.cmb.control_volume.rate_reaction_extent.items():
    #     iscale.constraint_scaling_transform(c, 1e3, overwrite=True)
        # iscale.set_scaling_factor(m.fs.cmb.rate_reaction_extent, 1e3)
        
    # Scale down equations to minimize residual values
    # fs.bhx2.tube.material_balances[0.0,H2O]
    # fs.soec.fc.flow_mol_eqn[0.0,0]
    # fs.soec.ac.mass_balance_eqn[0.0,1,O2]
    # # fs.spltf1.material_splitting_eqn[0.0,out,H2O]
    # # fs.hxh2.heat_transfer_equation[0.0] # already scaled
    # fs.hxh2.shell.material_balances[0.0,H2O]
    # fs.mxa1.material_mixing_equations[0.0,CO2]
    # # fs.mxa1.air_state[0.0].total_flow_balance
    # # fs.s10_expanded.temperature_equality[0.0]
    # # fs.s13_expanded.flow_mol_equality[0.0]

    for (t, j), c in m.fs.bhx2.tube.material_balances.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for (t, n), c in m.fs.soec.fc.flow_mol_eqn.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for (t, n, j), c in m.fs.soec.ac.mass_balance_eqn.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for (t, j), c in m.fs.hxh2.shell.material_balances.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for (t, j), c in m.fs.mxa1.material_mixing_equations.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for (t, n, j), c in m.fs.soec.ac.mass_balance_eqn.items():
        iscale.constraint_scaling_transform(c, 1e-6, overwrite=True)







def get_solver():
    use_idaes_solver_configuration_defaults()
    idaes.cfg.ipopt["options"]["nlp_scaling_method"] = "user-scaling"
    idaes.cfg.ipopt["options"]["tol"] = 1e-4
    # due to a lot of component mole fractions being on their lower bound of 0
    # bound push result in much longer solve times, so set it low.
    idaes.cfg.ipopt["options"]["halt_on_ampl_error"] = 'yes'
    idaes.cfg.ipopt["options"]["bound_push"] = 1e-9
    idaes.cfg.ipopt["options"]["linear_solver"] = "ma57"
    idaes.cfg.ipopt["options"]["max_iter"] = 400
    #idaes.cfg.ipopt["options"]["ma27_pivtol"] = 1e-1
    # idaes.cfg.ipopt["options"]["ma57_pivtol"] = 1e-1
    # idaes.cfg.ipopt["options"]["ma57_pivtolmax"] = 1e-1
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
        expr=m.fs.cmb.outlet.temperature[0],
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
    # tags["hxo2_area"] = iutil.ModelTag(
    #     expr=m.fs.hxo2.area,
    #     format_string="{:.3f}",
    #     display_units=pyo.units.m**2,
    #     doc="hxo2 area",
    # )
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
                # "o03":m.fs.hxo2.shell_outlet,
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
    import idaes.core.util.scaling as iscale
    jac, nlp = iscale.get_jacobian(m, scaled=True)
    # print("Extreme Jacobian entries:")
    sourceFile = open('extreme_jacobian.txt', 'w')
    for i in iscale.extreme_jacobian_entries(jac=jac, nlp=nlp, small=1e-6, large=1e3):
        print(f"    {i[0]:.2e}, [{i[1]}, {i[2]}]", file=sourceFile)
    sourceFile.close()

    # print("Unscaled constraints:")
    sourceFile2 = open('unscaled_constraints.txt', 'w')
    for c in iscale.unscaled_constraints_generator(m):
        print(f"    {c}", file=sourceFile2)
    sourceFile2.close()

    sourceFile3 = open('constraints_with_scale_factor.txt', 'w')
    # print("Scaled constraints by factor:")
    for c, s in iscale.constraints_with_scale_factor_generator(m):
        print(f"    {c}, {s}", file=sourceFile3)
    sourceFile3.close()

    # print("Badly scaled variables:")
    sourceFile4 = open('badly_scaled_var.txt', 'w')
    for v, sv in iscale.badly_scaled_var_generator(m, large=1e3, small=1e-4, zero=1e-12):
        print(f"    {v} -- {sv} -- {iscale.get_scaling_factor(v)}", file=sourceFile4)
    sourceFile.close()
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
    add_asu(m)
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
    # print('dof = ', degrees_of_freedom(m))
    # import sys
    # sys.exit('stop')
    additional_scaling(m)
    iscale.calculate_scaling_factors(m)
    # iscale.scale_arc_constraints(m)
    # print(iscale.get_scaling_factor(m.fs.aux_boiler_feed_pump.eq_work[0]))
    # print(pyo.value(m.fs.aux_boiler_feed_pump.eq_work[0]))

    # import sys
    # sys.exit('stop')
    # iscale.constraint_autoscale_large_jac(m)
    # m.strip_bounds()  # TODO - what is this method currently called?
    do_initialize(m, solver)
    display_input_tags(m)
    return m

def residual_checker(m):
    # Print model constraints - residual checker
    with open('model_residuals.txt','w') as f:
        for v in m.component_data_objects(pyo.Constraint, descend_into=True, active=True):
            residual_value = abs(v.body() - v.lower())
            if residual_value >= 1e-6:
                a = " FAIL"
            else:
                a = " PASS"
            residual_print = [str(v.name), ": value = ", str(residual_value), a, " \n"]
            f.writelines(residual_print)
def variables_above_bounds_generator(
        block, tol=1e-6, relative=True, skip_lb=False, skip_ub=False):
    """
    Generator which returns all Var components in a model which have a value
    within tol (default: relative) of a bound.

    Args:
        block : model to be studied
        tol : (relative) tolerance for inclusion in generator (default = 1e-4)
        relative : Boolean, use relative tolerance (default = True)
        skip_lb: Boolean to skip lower bound (default = False)
        skip_ub: Boolean to skip upper bound (default = False)

    Returns:
        A generator which returns all Var components block that are close to a
        bound
    """
    for v in block.component_data_objects(
            ctype=pyo.Var, active=True, descend_into=True):
        # To avoid errors, check that v has a value
        if v.value is None:
            continue

        if relative:
            # First, determine absolute tolerance to apply to bounds
            if v.ub is not None and v.lb is not None:
                # Both upper and lower bounds, apply tol to (upper - lower)
                atol = pyo.value((v.ub - v.lb)*tol)
            elif v.ub is not None:
                # Only upper bound, apply tol to bound value
                atol = abs(pyo.value(v.ub*tol))
            elif v.lb is not None:
                # Only lower bound, apply tol to bound value
                atol = abs(pyo.value(v.lb*tol))
            else:
                continue
        else:
            atol = tol

        if v.ub is not None and not skip_lb and pyo.value(v.ub - v.value) >= atol:
            yield v
        elif (v.lb is not None and not skip_ub and
              pyo.value(v.value - v.lb) >= atol):
            yield v


def check_algebraic_jacobian():
    """ Main function to set up and run dynamic flowsheet simulation """
    # start_time = time.time()

    # # Number of time steps, we need to solve for at least 20 time steps
    # # Please try to increase num_step to test the convergence issue
    # # We may also want to revise the input_profile() function to test the convergence
    # # Eventually we want to provide the start and end points for the profile and
    # # let the optimizer to find the best profile 
    # num_step = 8

    # # The time step size in second
    # step_size = 30

    # Set up first dynamic model
    m_ss = get_model()

    # Check degree of freedom for the dynamic model
    dof = degrees_of_freedom(m_ss)
    assert dof == 0

    time = m_ss.fs.time
    t1 = 0

    scalar_vars, dae_vars = flatten_dae_components(
            m_ss,
            time,
            pyo.Var
            )
    scalar_cons, dae_cons = flatten_dae_components(
            m_ss,
            time,
            pyo.Constraint,
            )

    var_partition, con_partition = categorize_dae_variables_and_constraints(
            m_ss,
            dae_vars,
            dae_cons,
            time,
            index=0 # index set equal to zero for steady-state systems
            )

    alg_vars = [var[t1] for var in var_partition[VC.ALGEBRAIC]]
    alg_cons = [con[t1] for con in con_partition[CC.ALGEBRAIC]]

    igraph = IncidenceGraphInterface()
    matching = igraph.maximum_matching(alg_vars, alg_cons)

    var_dmp, con_dmp = igraph.dulmage_mendelsohn(alg_vars, alg_cons)

    matched_var_set = ComponentSet(matching.values())

    print("Unmatched constraints:")
    for con in alg_cons:
        if con not in matching:
            print("  %s" % con.name)

    print("Unmatched variables:")
    for var in alg_vars:
        if var not in matched_var_set:
            print("  %s" % var.name)

    # At this point check the sizes of the underconstrained and
    # overconstrained systems; decide which to debug.
    #import pdb; pdb.set_trace()

    print("\n- - -")
    print("Dulmage-Mendelsohn partition:")

    print("Unmatched constraints:")
    for con in con_dmp[0]:
        print(con.name)
    for con, var in zip(con_dmp[1], var_dmp[2]):
        print(con.name, var.name)
if __name__ == "__main__":
    # m = check_algebraic_jacobian()
    m = get_model()
    # write_pfd_results(m, "soec_init.svg")
    print('dof = ', degrees_of_freedom(m))
    # check_scaling(m)
    residual_checker(m)
    generator=variables_above_bounds_generator(m)
    generator2=variables_near_bounds_generator(m, tol=1e-6)
           

    sourceFile = open('model_bounds.txt', 'w')          
    for v in generator:
        print("var_name: %a, var_lb: %a, var_value: %a, var_ub: %a" %(
                v.name, pyo.value(v.lb), v.value, pyo.value(v.ub)),
                file=sourceFile)  
    sourceFile.close()              

    sourceFile1 = open('vars_near_bounds.txt', 'w')          
    for v in generator2:
        print("var_name: %a, var_lb: %a, var_value: %a, var_ub: %a" %(
                v.name, pyo.value(v.lb), v.value, pyo.value(v.ub)),
                file=sourceFile1)  
    sourceFile1.close() 
    
    
    # # Print all variables and constraints in the model - update for BFB
    # with open('model_bounds.txt','w') as f:
    #     for v in m.component_objects(pyo.Var, descend_into=True):
    #         print("var_name:" + v.name,
    #               "var_lb", pyo.value(v.lb),
    #               "var_value", v.value,
    #               "var_ub", pyo.value(v.lb))
            # v.pprint(ostream=f)

        # return ( headers,
        #          self._data.items(),
        #          ( "Lower","Value","Upper","Fixed","Stale","Domain"),
        #          lambda k, v: [ value(v.lb),
        #                         v.value,
        #                         value(v.ub),
        #                         v.fixed,
        #                         v.stale,
        #                         v.domain
        #                         ]
        #          )
    
# def unscaled_variables_generator(blk, descend_into=True, include_fixed=False):
#     """Generator for unscaled variables

#     Args:
#         block

#     Yields:
#         variables with no scale factor
#     """
#     for v in blk.component_data_objects(pyo.Var, descend_into=descend_into):
#         if v.fixed and not include_fixed:
#             continue
#         if get_scaling_factor(v) is None:
#             yield v