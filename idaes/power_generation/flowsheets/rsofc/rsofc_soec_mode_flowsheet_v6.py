###############################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
###############################################################################

__author__ = "Chinedu Okoli", "John Eslick", "Alex Noring"

# Import python modules
import os
import csv
import numpy as np
# import matplotlib.pyplot as plt

# Import pyomo modules
import pyomo.environ as pyo
from pyomo.util.calc_var_value import calculate_variable_from_constraint
from pyomo.network import Arc, Port
from pyomo.common.fileutils import this_file_dir
# import pyomo.common.errors

# Import modules from core
from idaes.core import FlowsheetBlock
import idaes.core.util as iutil
import idaes.core.util.tables as tables
import idaes.core.util.scaling as iscale
import idaes.core.util.initialization as iinit
from idaes.core.util import copy_port_values
import idaes.core.plugins

# Import unit models
from idaes.power_generation.unit_models.helm import (
    # HelmMixer,
    # MomentumMixingType,
    HelmSplitter,
    HelmIsentropicCompressor
)
import idaes.generic_models.unit_models as gum  # generic unit models
import idaes.power_generation.unit_models as pum  # power unit models
from idaes.power_generation.flowsheets.sofc.surrogates.cpu import CPU
from idaes.generic_models.unit_models.heat_exchanger import (
    delta_temperature_underwood_callback
)
from idaes.generic_models.unit_models.pressure_changer import \
    ThermodynamicAssumption
from idaes.generic_models.unit_models.separator import SplittingType

# Import properties
from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock,
)
from idaes.generic_models.properties.core.generic.generic_reaction import (
    GenericReactionParameterBlock,
)
from idaes.power_generation.properties.natural_gas_PR import (
    get_prop, get_rxn, EosType)
from idaes.power_generation.flowsheets.rsofc.properties. \
    CO2_H2O_Ideal_VLE_scaled import configuration as CO2_H2O_VLE_config
from idaes.core.solvers import use_idaes_solver_configuration_defaults
from idaes.generic_models.properties import iapws95
# from idaes.generic_models.properties.helmholtz.helmholtz import (
#     HelmholtzThermoExpressions as ThermoExpr,
# )

# Import costing
from idaes.power_generation.flowsheets.rsofc import rsofc_costing as rsofc_cost

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
    # if m is None:
    #     m = pyo.ConcreteModel(name)
    if not hasattr(m, "soec_fs"):
        m.soec_fs = FlowsheetBlock(default={"dynamic": False})
    return m


def add_properties(fs):
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
    fs.fg_prop = GenericParameterBlock(
        default=get_prop(components=fuel_comp, phases=["Vap"],
                         eos=EosType.IDEAL)
    )
    fs.fg_prop.set_default_scaling("mole_frac_comp", 1e2)
    fs.fg_prop.set_default_scaling("mole_frac_phase_comp", 1e2)
    rxns = {  # reactions and key components for conversion
        "ch4_cmb": "CH4",
        "c2h6_cmb": "C2H6",
        "c3h8_cmb": "C3H8",
        "c4h10_cmb": "C4H10",
    }
    fs.rxns = rxns
    fs.fg_combust = GenericReactionParameterBlock(
        default=get_rxn(fs.fg_prop, rxns))
    fs.water_prop = iapws95.Iapws95ParameterBlock()
    fs.h2_prop = GenericParameterBlock(
        default=get_prop(components={"H2", "H2O"},
                         phases=["Vap"], eos=EosType.IDEAL)
    )
    fs.o2_prop = GenericParameterBlock(
        default=get_prop(components={"O2", "H2O"},
                         phases=["Vap"], eos=EosType.IDEAL)
    )
    fs.air_prop = GenericParameterBlock(
        default=get_prop(components=air_comp,
                         phases=["Vap"], eos=EosType.IDEAL)
    )
    fs.air_prop.set_default_scaling("mole_frac_comp", 1e2)
    fs.air_prop.set_default_scaling("mole_frac_phase_comp", 1e2)
    fs.h2_compress_prop = GenericParameterBlock(
        default=get_prop(components={"H2"}, phases=["Vap"], eos=EosType.PR)
    )
    fs.h2_compress_prop.set_default_scaling("mole_frac_comp", 1e2)
    fs.h2_compress_prop.set_default_scaling("mole_frac_phase_comp", 1e2)

    fs.CO2_H2O_VLE = GenericParameterBlock(default=CO2_H2O_VLE_config)
    fs.CO2_H2O_VLE.set_default_scaling("mole_frac_comp", 1e2)
    fs.CO2_H2O_VLE.set_default_scaling("mole_frac_phase_comp", 1e2)


def add_asu(fs):
    ###########################################################################
    #  Build unit operations
    ###########################################################################
    fs.air_compressor_s1 = gum.PressureChanger(
        default={"compressor": True,
                 "property_package": fs.air_prop,
                 "thermodynamic_assumption":
                     ThermodynamicAssumption.isentropic})

    fs.intercooler_s1 = gum.Heater(
        default={"property_package": fs.air_prop,
                 "has_pressure_change": True})

    fs.air_compressor_s2 = gum.PressureChanger(
        default={"compressor": True,
                 "property_package": fs.air_prop,
                 "thermodynamic_assumption":
                     ThermodynamicAssumption.isentropic})

    fs.intercooler_s2 = gum.Heater(
        default={"property_package": fs.air_prop,
                 "has_pressure_change": True})

    fs.ASU = gum.Separator(
        default={"outlet_list": ["N2_outlet", "O2_outlet"],
                 "split_basis": SplittingType.componentFlow,
                 "property_package": fs.air_prop})

    ###########################################################################
    #  Specify performance variables of unit operations
    ###########################################################################
    # air compressors and intercoolers
    fs.air_compressor_s1.outlet.pressure.fix(111422)  # Pa (34 psia)
    fs.air_compressor_s1.efficiency_isentropic.fix(0.84)
    fs.intercooler_s1.outlet.temperature.fix(310.93)  # K (100 F)
    fs.intercooler_s1.deltaP.fix(-3447)  # Pa (-0.5 psi)
    fs.air_compressor_s2.outlet.pressure.fix(130686)  # Pa (79 psia)
    fs.air_compressor_s2.efficiency_isentropic.fix(0.84)
    fs.intercooler_s2.outlet.temperature.fix(310.93)  # K (100 F)
    fs.intercooler_s2.deltaP.fix(-3447)  # Pa (-0.5 psi)

    # air seperation unit
    fs.ASU.split_fraction[0, "O2_outlet", "CO2"].fix(0)
    fs.ASU.split_fraction[0, "O2_outlet", "H2O"].fix(0)
    fs.ASU.split_fraction[0, "O2_outlet", "N2"].fix(0.0005)
    fs.ASU.split_fraction[0, "O2_outlet", "O2"].fix(0.9691)
    fs.ASU.split_fraction[0, "O2_outlet", "Ar"].fix(0.0673)

    ###########################################################################
    #  Add stream connections
    ###########################################################################
    # Arcs for ASU, oxycombustor and CPU
    fs.STAGE_1_OUT = Arc(
        source=fs.air_compressor_s1.outlet,
        destination=fs.intercooler_s1.inlet)

    fs.IC_1_OUT = Arc(
        source=fs.intercooler_s1.outlet,
        destination=fs.air_compressor_s2.inlet)

    fs.STAGE_2_OUT = Arc(
        source=fs.air_compressor_s2.outlet,
        destination=fs.intercooler_s2.inlet)

    fs.ba01 = Arc(
        source=fs.intercooler_s2.outlet,
        destination=fs.ASU.inlet)


def add_preheater(fs):
    ###########################################################################
    #  Build unit operations
    ###########################################################################
    fs.oxygen_preheater = gum.HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_underwood_callback,
            "shell": {"property_package": fs.air_prop},
            "tube": {"property_package": fs.air_prop},
        }
    )
    fs.ng_preheater = gum.HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_underwood_callback,
            "shell": {"property_package": fs.air_prop},
            "tube": {"property_package": fs.fg_prop},
        }
    )
    fs.preheat_split = gum.Separator(
        default={
            "property_package": fs.air_prop,
            "outlet_list": ["air", "ng"],
        }
    )

    ###########################################################################
    #  Specify performance variables of unit operations
    ###########################################################################
    fs.oxygen_preheater.overall_heat_transfer_coefficient.fix(100)
    fs.oxygen_preheater.delta_temperature_in.fix(30)  # fix DT for pinch side

    fs.ng_preheater.overall_heat_transfer_coefficient.fix(100)
    fs.ng_preheater.delta_temperature_in.fix(30)  # fix DT for pinch side

    fs.preheat_split.split_fraction[:, "air"].fix(0.9)

    ###########################################################################
    #  Add stream connections
    ###########################################################################
    fs.o04 = Arc(
        source=fs.preheat_split.ng,
        destination=fs.ng_preheater.shell_inlet
    )
    fs.o05 = Arc(
        source=fs.preheat_split.air,
        destination=fs.oxygen_preheater.shell_inlet
    )
    fs.ba02 = Arc(
        source=fs.ASU.O2_outlet,
        destination=fs.oxygen_preheater.tube_inlet
        )


def add_combustor(fs):
    ###########################################################################
    #  Build unit operations
    ###########################################################################
    fs.pre_oxycombustor_translator = gum.Translator(
        default={"outlet_state_defined": True,
                 "inlet_property_package": fs.air_prop,
                 "outlet_property_package": fs.fg_prop}
        )
    fs.cmb_mix = gum.Mixer(
        default={
            "property_package": fs.fg_prop,
            "inlet_list": ["ng", "air"],
            "momentum_mixing_type": gum.MomentumMixingType.none}
        )
    fs.cmb = gum.StoichiometricReactor(
        default={
            "property_package": fs.fg_prop,
            "reaction_package": fs.fg_combust,
            "has_pressure_change": False,
            "has_heat_of_reaction": True,
            "has_heat_transfer": True  # Use to add Q to water/steam in tubes
            }
        )
    fs.post_oxycombustor_translator = gum.Translator(
        default={"outlet_state_defined": True,
                 "inlet_property_package": fs.fg_prop,
                 "outlet_property_package": fs.air_prop}
        )

    ###########################################################################
    #  Specify performance variables of unit operations
    ###########################################################################
    # Additional constraints to specify the pre combustor translator block
    @fs.pre_oxycombustor_translator.Constraint(fs.time)
    def pre_oxycombustor_translator_F(b, t):
        return b.inlet.flow_mol[t] == b.outlet.flow_mol[t]

    @fs.pre_oxycombustor_translator.Constraint(fs.time)
    def pre_oxycombustor_translator_T(b, t):
        return b.inlet.temperature[t] == b.outlet.temperature[t]

    @fs.pre_oxycombustor_translator.Constraint(fs.time)
    def pre_oxycombustor_translator_P(b, t):
        return b.inlet.pressure[t] == b.outlet.pressure[t]

    @fs.pre_oxycombustor_translator.Constraint(fs.time,
                                                 fs.air_prop.component_list)
    def pre_oxycombustor_translator_x(b, t, j):
        return b.inlet.mole_frac_comp[t, j] == b.outlet.mole_frac_comp[t, j]

    for j in fs.fg_prop.component_list:
        if j not in fs.air_prop.component_list:
            fs.pre_oxycombustor_translator.outlet.mole_frac_comp[0, j].fix(0)

    # Additional constraints to specify the mixer
    @fs.cmb_mix.Constraint(fs.time)
    def pressure_eqn(b, t):
        return b.mixed_state[t].pressure == b.ng_state[t].pressure

    # Additional constraints to specify the combustor reactions
    @fs.cmb.Constraint(fs.time, fs.rxns.keys())
    def reaction_extent(b, t, r):
        k = fs.rxns[r]
        prp = b.control_volume.properties_in[t]
        stc = -fs.fg_combust.rate_reaction_stoichiometry[r, "Vap", k]
        extent = b.rate_reaction_extent[t, r]
        return extent*stc == prp.flow_mol*prp.mole_frac_comp[k]

    # Additional constraints to specify the post combustor  translator block
    @fs.post_oxycombustor_translator.Constraint(fs.time)
    def post_oxycombustor_translator_F(b, t):
        return b.inlet.flow_mol[t] == b.outlet.flow_mol[t]

    @fs.post_oxycombustor_translator.Constraint(fs.time)
    def post_oxycombustor_translator_T(b, t):
        return b.inlet.temperature[t] == b.outlet.temperature[t]

    @fs.post_oxycombustor_translator.Constraint(fs.time)
    def post_oxycombustor_translator_P(b, t):
        return b.inlet.pressure[t] == b.outlet.pressure[t]

    @fs.post_oxycombustor_translator.Constraint(fs.time,
                                                  fs.air_prop.component_list)
    def post_oxycombustor_translator_x(b, t, j):
        return b.inlet.mole_frac_comp[t, j] == b.outlet.mole_frac_comp[t, j]

    ###########################################################################
    #  Add stream connections
    ###########################################################################
    fs.cmb_mix_in = Arc(  # TODO - rename arc
        source=fs.oxygen_preheater.tube_outlet,
        destination=fs.pre_oxycombustor_translator.inlet
    )
    fs.ba03 = Arc(
        source=fs.pre_oxycombustor_translator.outlet,
        destination=fs.cmb_mix.air
    )
    fs.bng01 = Arc(
        source=fs.ng_preheater.tube_outlet, destination=fs.cmb_mix.ng
    )
    fs.bng02 = Arc(
        source=fs.cmb_mix.outlet, destination=fs.cmb.inlet
    )
    fs.cmb_out = Arc(  # TODO - rename arc
        source=fs.cmb.outlet,
        destination=fs.post_oxycombustor_translator.inlet
    )
    fs.fg01 = Arc(source=fs.post_oxycombustor_translator.outlet,
                    destination=fs.air_preheater_2.shell_inlet)


def add_aux_boiler_steam(fs):
    ###########################################################################
    #  Build unit operations
    ###########################################################################
    fs.bhx2 = gum.Heater(  # in-bed heat exchanger to the oxycombustor
        default={"has_pressure_change": False,
                 "property_package": fs.water_prop})

    fs.main_steam_split = HelmSplitter(
        default={
            "property_package": fs.water_prop,
            "outlet_list": ["h_side", "o_side"]
        }
    )
    fs.aux_boiler_feed_pump = HelmIsentropicCompressor(
        default={"property_package": fs.water_prop}
    )
    fs.bhx1 = gum.HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_underwood_callback,
            "shell": {"property_package": fs.h2_prop},
            "tube": {"property_package": fs.water_prop},
        }
    )

    ###########################################################################
    #  Specify performance variables of unit operations
    ###########################################################################
    # fs.bhx2.shell_outlet.temperature[0].fix(700)
    # fs.bhx2.overall_heat_transfer_coefficient.fix(100)
    # fs.bhx1.area.fix(400)
    # htpx(T=None, P=None, x=None)
    # bhx2 - inbed combustor heat exchanger spec
    # TODO - Rewrite this spec in a more elegant way - Q comes from cmb
    h_bhx2 = pyo.value(iapws95.htpx(1073.15*pyo.units.K,1.1e5*pyo.units.Pa)) # enthalpy outlet
    fs.bhx2.outlet.enth_mol.fix(h_bhx2)  # K (100 F) # unfix after initalize and spec Q from cmb

    fs.bhx1.overall_heat_transfer_coefficient.fix(100)
    fs.bhx1.delta_temperature_out.fix(100)  # fix DT for pinch side

    fs.aux_boiler_feed_pump.outlet.pressure.fix(1.1e5)
    fs.aux_boiler_feed_pump.efficiency_isentropic.fix(0.85)

    fs.main_steam_split.split_fraction[:, "h_side"].fix(0.9999)

    ###########################################################################
    #  Add stream connections
    ###########################################################################
    fs.s01 = Arc(
        source=fs.aux_boiler_feed_pump.outlet,
        destination=fs.bhx1.tube_inlet
    )
    fs.s02 = Arc(
        source=fs.bhx1.tube_outlet,
        destination=fs.bhx2.inlet
    )
    fs.s03 = Arc(
        source=fs.bhx2.outlet,
        destination=fs.main_steam_split.inlet
    )


def add_soec_air_side_units(fs):
    ###########################################################################
    #  Build unit operations
    ###########################################################################
    fs.air_blower = gum.PressureChanger(
        default={'compressor': True,
                 'property_package': fs.air_prop,
                 'thermodynamic_assumption':
                     ThermodynamicAssumption.isentropic}
    )
    fs.air_preheater_1 = gum.HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_underwood_callback,
            "shell": {"property_package": fs.air_prop},
            "tube": {"property_package": fs.air_prop},
        }
    )
    fs.air_preheater_2 = gum.HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_underwood_callback,
            "shell": {"property_package": fs.air_prop},
            "tube": {"property_package": fs.air_prop},
        }
    )
    ###########################################################################
    #  Specify performance variables of unit operations
    ###########################################################################
    fs.air_blower.outlet.pressure.fix(1.1e5)  # Pa (15.3 psi)
    fs.air_blower.efficiency_isentropic.fix(0.8)

    fs.air_preheater_1.overall_heat_transfer_coefficient.fix(100)
    fs.air_preheater_1.delta_temperature_in.fix(50)  # fix DT for pinch side

    fs.air_preheater_2.tube_outlet.temperature[0].fix(1073.15)  # Known value
    fs.air_preheater_2.overall_heat_transfer_coefficient.fix(100)
    # fs.air_preheater_2.area.fix(1000)

    ###########################################################################
    #  Add stream connections
    ###########################################################################
    fs.a01 = Arc(source=fs.air_blower.outlet,
                   destination=fs.air_preheater_1.tube_inlet)
    fs.a02 = Arc(source=fs.air_preheater_1.tube_outlet,
                   destination=fs.air_preheater_2.tube_inlet)


def add_soec_unit(fs):
    ###########################################################################
    #  Build unit operations
    ###########################################################################
    fs.soec = pum.IsothermalSofc(
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
    fs.spltf1 = gum.Separator(
        default={
            "property_package": fs.h2_prop,
            "outlet_list": ["out", "recycle"],
        }
    )
    fs.splta1 = gum.Separator(
        default={
            "property_package": fs.air_prop,
            "outlet_list": ["out", "recycle"],
        }
    )

    ###########################################################################
    #  Specify performance variables of unit operations
    ###########################################################################
    # soec performance variables
    fs.soec.E_cell.fix(1.28)  # TODO - unfix after initialize
    fs.soec.el.thickness.fix(9e-6)
    fs.soec.fe.thickness.fix(1e-3)
    fs.soec.ae.thickness.fix(20e-6)
    fs.soec.length.fix(0.05)
    fs.soec.width.fix(0.05)
    fs.soec.k_ae.fix(26.1e7)
    fs.soec.eact_ae.fix(120000)
    fs.soec.alpha_ae.fix(0.4)
    fs.soec.k_fe.fix(1.35e10)
    fs.soec.eact_fe.fix(110000)
    fs.soec.alpha_fe.fix(0.5)
    fs.soec.fe.k_res.fix(2.98e-5)
    fs.soec.fe.E_res.fix(-1392)
    fs.soec.ae.k_res.fix(8.114e-5)
    fs.soec.ae.E_res.fix(600)
    fs.soec.el.k_res.fix(2.94e-5)
    fs.soec.el.E_res.fix(10350)
    fs.soec.fc.thickness.fix(0.002)
    fs.soec.ac.thickness.fix(0.002)
    fs.soec.fe.porosity.fix(0.48)
    fs.soec.fe.tortuosity.fix(5.4)
    fs.soec.ae.porosity.fix(0.48)
    fs.soec.ae.tortuosity.fix(5.4)
    temperature = 1073.15  # TODO - should temperatures be fixed here?
    fs.soec.el.temperature.fix(temperature)
    fs.soec.fc.temperature.fix(temperature)
    fs.soec.ac.temperature.fix(temperature)
    fs.soec.fe.temperature.fix(temperature)
    fs.soec.ae.temperature.fix(temperature)

    # performance variables for other units
    fs.spltf1.split_fraction[:, "out"].fix(0.98)
    fs.splta1.split_fraction[:, "out"].fix(0.9999)

    ###########################################################################
    #  Add stream connections
    ###########################################################################
    fs.h01 = Arc(source=fs.soec.outlet_fc_mult,
                   destination=fs.spltf1.inlet)
    fs.o01 = Arc(source=fs.soec.outlet_ac_mult,
                   destination=fs.splta1.inlet)


def add_soec_inlet_mix(fs):
    ###########################################################################
    #  Build unit operations
    ###########################################################################
    fs.mxf1 = gum.Mixer(
        default={
            "property_package": fs.h2_prop,
            "inlet_list": ["water", "recycle"],
            "momentum_mixing_type": gum.MomentumMixingType.none,
        }
    )
    fs.mxa1 = gum.Mixer(
        default={
            "property_package": fs.air_prop,
            "inlet_list": ["air", "recycle"],
            "momentum_mixing_type": gum.MomentumMixingType.none,
        }
    )

    ###########################################################################
    #  Specify performance variables of unit operations
    ###########################################################################
    @fs.mxf1.Constraint(fs.time)
    def fmxpress_eqn(b, t):
        return b.mixed_state[t].pressure == b.water_state[t].pressure

    @fs.mxa1.Constraint(fs.time)
    def amxpress_eqn(b, t):
        return b.mixed_state[t].pressure == b.air_state[t].pressure

    ###########################################################################
    #  Add stream connections
    ###########################################################################
    # Add ports to connect pure steam to steam + h2
    fs.main_steam_split._temperature_h_side_ref = pyo.Reference(
        fs.main_steam_split.h_side_state[:].temperature
    )

    @fs.main_steam_split.Expression(fs.time, fs.soec.fc.config.comp_list)
    def h_side_mole_frac_expr(b, t, i):
        if i == "H2O":
            return 1
        else:
            return 0

    fs.main_steam_split.h_side_adapt = Port(
        rule=lambda b: {
            "flow_mol": fs.main_steam_split._flow_mol_h_side_ref,
            "pressure": fs.main_steam_split._pressure_h_side_ref,
            "temperature": fs.main_steam_split._temperature_h_side_ref,
            "mole_frac_comp": fs.main_steam_split.h_side_mole_frac_expr,
        }
    )

    fs.s10 = Arc(
        source=fs.main_steam_split.h_side_adapt,
        destination=fs.mxf1.water,
    )

    fs.hr01 = Arc(
        source=fs.spltf1.recycle,
        destination=fs.mxf1.recycle,
    )

    fs.s12 = Arc(source=fs.mxf1.outlet,
                   destination=fs.soec.inlet_fc_mult)
    fs.s13 = Arc(source=fs.mxa1.outlet,
                   destination=fs.soec.inlet_ac_mult)

    fs.a03 = Arc(
        source=fs.air_preheater_2.tube_outlet,
        destination=fs.mxa1.air
    )


def add_h2_compressor(fs):
    ###########################################################################
    #  Build unit operations
    ###########################################################################
    fs.hcmp_ic01 = gum.Heater(
                default={"property_package": fs.h2_compress_prop})
    fs.hcmp01 = gum.Compressor(
                default={"property_package": fs.h2_compress_prop,
                         "thermodynamic_assumption":
                             ThermodynamicAssumption.isentropic})
    fs.hcmp_ic02 = gum.Heater(
                default={"property_package": fs.h2_compress_prop})
    fs.hcmp02 = gum.Compressor(
                default={"property_package": fs.h2_compress_prop,
                         "thermodynamic_assumption":
                             ThermodynamicAssumption.isentropic})
    fs.hcmp_ic03 = gum.Heater(
                default={"property_package": fs.h2_compress_prop})
    fs.hcmp03 = gum.Compressor(
        default={"property_package": fs.h2_compress_prop,
                 "thermodynamic_assumption":
                     ThermodynamicAssumption.isentropic})
    fs.hcmp_ic04 = gum.Heater(
                default={"property_package": fs.h2_compress_prop})
    fs.hcmp04 = gum.Compressor(
                default={"property_package": fs.h2_compress_prop,
                         "thermodynamic_assumption":
                             ThermodynamicAssumption.isentropic})

    ###########################################################################
    #  Specify performance variables of unit operations
    ###########################################################################
    @fs.bhx1.Expression(fs.time, {"H2"})
    def waterless_mole_frac_expr(b, t, i):
        return 1

    @fs.bhx1.Expression(fs.time)
    def waterless_flow_expr(b, t):
        return (
            fs.bhx1._flow_mol_shell_outlet_ref[t]
            * fs.bhx1._mole_frac_comp_shell_outlet_ref[t, "H2"]
        )

    fs.bhx1.shell_outlet_drop_water = Port(
        rule=lambda b: {
            "flow_mol": fs.bhx1.waterless_flow_expr,
            "pressure": fs.bhx1._pressure_shell_outlet_ref,
            "temperature": fs.bhx1._temperature_shell_outlet_ref,
            "mole_frac_comp": fs.bhx1.waterless_mole_frac_expr,
        }
    )

    fs.hcmp_ic01.outlet.temperature.fix(320)
    fs.hcmp01.outlet.pressure.fix(40e5)
    fs.hcmp01.efficiency_isentropic.fix(0.9)

    fs.hcmp_ic02.outlet.temperature.fix(320)
    fs.hcmp02.outlet.pressure.fix(80e5)
    fs.hcmp02.efficiency_isentropic.fix(0.9)

    fs.hcmp_ic03.outlet.temperature.fix(320)
    fs.hcmp03.outlet.pressure.fix(160e5)
    fs.hcmp03.efficiency_isentropic.fix(0.9)

    fs.hcmp_ic04.outlet.temperature.fix(320)
    fs.hcmp04.outlet.pressure.fix(320e5)
    fs.hcmp04.efficiency_isentropic.fix(0.9)

    ###########################################################################
    #  Add stream connections
    ###########################################################################
    fs.h03 = Arc(
        source=fs.bhx1.shell_outlet_drop_water,
        destination=fs.hcmp_ic01.inlet
    )
    fs.h04 = Arc(source=fs.hcmp_ic01.outlet,
                   destination=fs.hcmp01.inlet)
    fs.h05 = Arc(source=fs.hcmp01.outlet,
                   destination=fs.hcmp_ic02.inlet)
    fs.h06 = Arc(source=fs.hcmp_ic02.outlet,
                   destination=fs.hcmp02.inlet)
    fs.h07 = Arc(source=fs.hcmp02.outlet,
                   destination=fs.hcmp_ic03.inlet)
    fs.h08 = Arc(source=fs.hcmp_ic03.outlet,
                   destination=fs.hcmp03.inlet)
    fs.h09 = Arc(source=fs.hcmp03.outlet,
                   destination=fs.hcmp_ic04.inlet)
    fs.h10 = Arc(source=fs.hcmp_ic04.outlet,
                   destination=fs.hcmp04.inlet)


def add_hrsg_and_cpu(fs):
    ###########################################################################
    #  Build unit operations
    ###########################################################################
    fs.soec_air_HRSG_1 = gum.Heater(
        default={"has_pressure_change": False,
                 "property_package": fs.air_prop})
    fs.soec_air_HRSG_2 = gum.Heater(
        default={"has_pressure_change": False,
                 "property_package": fs.air_prop})
    fs.fluegas_HRSG = gum.Heater(
        default={"has_pressure_change": False,
                 "property_package": fs.air_prop})

    fs.CPU_translator = gum.Translator(
        default={"outlet_state_defined": True,
                 "inlet_property_package": fs.air_prop,
                 "outlet_property_package": fs.CO2_H2O_VLE})
    fs.condenser = gum.Heater(
        default={"has_pressure_change": False,
                 "property_package": fs.CO2_H2O_VLE})

    fs.flash = gum.Flash(
        default={"has_heat_transfer": True,
                 "has_pressure_change": True,
                 "property_package": fs.CO2_H2O_VLE})
    fs.CPU = CPU()

    ###########################################################################
    #  Specify performance variables of unit operations
    ###########################################################################
    fs.soec_air_HRSG_1.outlet.temperature.fix(405)  # K

    fs.soec_air_HRSG_2.outlet.temperature.fix(405)  # K

    fs.fluegas_HRSG.outlet.temperature.fix(405)  # K

    fs.condenser.outlet.temperature.fix(310.9)  # K (100 F)

    fs.flash.control_volume.properties_out[0].temperature.fix(310.9)  # K
    fs.flash.control_volume.properties_out[0].pressure.fix(1013529)  # Pa

    # CPU translator constraints
    @fs.CPU_translator.Constraint(fs.time)
    def CPU_translator_F(b, t):
        return b.inlet.flow_mol[t] == b.outlet.flow_mol[t]

    @fs.CPU_translator.Constraint(fs.time)
    def CPU_translator_T(b, t):
        return b.inlet.temperature[t] == b.outlet.temperature[t]

    @fs.CPU_translator.Constraint(fs.time)
    def CPU_translator_P(b, t):
        return b.inlet.pressure[t] == b.outlet.pressure[t]

    @fs.CPU_translator.Constraint(fs.time, fs.CO2_H2O_VLE.component_list)
    def CPU_translator_x(b, t, j):
        return b.inlet.mole_frac_comp[t, j] == b.outlet.mole_frac_comp[t, j]
    for j in fs.air_prop.component_list:
        if j not in fs.CO2_H2O_VLE.component_list:
            fs.pre_oxycombustor_translator.outlet.mole_frac_comp[0, j].fix(0)

    ###########################################################################
    #  Add stream connections
    ###########################################################################
    fs.a06 = Arc(
        source=fs.oxygen_preheater.shell_outlet,
        destination=fs.soec_air_HRSG_1.inlet
    )
    fs.a07 = Arc(
        source=fs.ng_preheater.shell_outlet,
        destination=fs.soec_air_HRSG_2.inlet
    )
    fs.fg02 = Arc(source=fs.air_preheater_2.shell_outlet,
                    destination=fs.fluegas_HRSG.inlet)
    fs.fg03 = Arc(source=fs.fluegas_HRSG.outlet,
                    destination=fs.CPU_translator.inlet)
    fs.fg04 = Arc(source=fs.CPU_translator.outlet,
                    destination=fs.condenser.inlet)
    fs.fg05 = Arc(source=fs.condenser.outlet,
                    destination=fs.flash.inlet)

    # TODO - copy_port method used instead of propagate_state because of issue
    # TODO - noticed in tags with propagating flash.vap_outlet to CPU.inlet
    # fs.fg06 = Arc(source=fs.flash.vap_outlet,
    #                 destination=fs.CPU.inlet)
    @fs.Constraint(fs.time)
    def CPU_inlet_F(fs, t):
        return (1e-3*fs.CPU.inlet.flow_mol[t] ==
                1e-3*fs.flash.vap_outlet.flow_mol[t])

    @fs.Constraint(fs.time)
    def CPU_inlet_T(fs, t):
        return (1e-3*fs.CPU.inlet.temperature[t] ==
                1e-3*fs.flash.vap_outlet.temperature[t])

    @fs.Constraint(fs.time)
    def CPU_inlet_P(fs, t):
        return (1e-3*fs.CPU.inlet.pressure[t] ==
                1e-3*fs.flash.vap_outlet.pressure[t])

    @fs.Constraint(fs.time, fs.CO2_H2O_VLE.component_list)
    def CPU_inlet_x(fs, t, j):
        return (1e2*fs.CPU.inlet.mole_frac_comp[t, j] ==
                1e2*fs.flash.vap_outlet.mole_frac_comp[t, j])


def add_more_hx_connections(fs):
    fs.h02 = Arc(source=fs.spltf1.out,
                   destination=fs.bhx1.shell_inlet)
    fs.o02 = Arc(source=fs.splta1.out,
                   destination=fs.air_preheater_1.shell_inlet)
    fs.o03 = Arc(
        source=fs.air_preheater_1.shell_outlet,
        destination=fs.preheat_split.inlet
        )


def add_design_constraints(fs):
    fs.soec_heat_duty = pyo.Var(fs.time, units=pyo.units.W)

    @fs.Constraint(fs.time)
    def heat_duty_soec_zero_eqn(b, t):
        return b.soec.heat_duty[t] == b.soec_heat_duty[t]

    fs.cmb_temperature = pyo.Var(fs.time, initialize=2000,
                                   units=pyo.units.K)

    @fs.Constraint(fs.time)
    def cmb_temperature_eqn(b, t):
        return fs.cmb.outlet.temperature[t] == fs.cmb_temperature[t]

    fs.soec_steam_temperature = pyo.Var(fs.time, initialize=1073.15,
                                          units=pyo.units.K)

    @fs.Constraint(fs.time)
    def soec_steam_temperature_eqn(b, t):
        return (fs.bhx2.control_volume.properties_out[t].temperature ==
                fs.soec_steam_temperature[t])

    @fs.Constraint(fs.time)
    def air_to_combustor(b, t):
        XS = 1.09  # excess oxygen

        F_CH4 = (fs.ng_preheater.tube_inlet.flow_mol[t] *
                 fs.ng_preheater.tube_inlet.mole_frac_comp[(t, "CH4")])
        F_C2H6 = (fs.ng_preheater.tube_inlet.flow_mol[t] *
                  fs.ng_preheater.tube_inlet.mole_frac_comp[(t, "C2H6")])
        F_C3H8 = (fs.ng_preheater.tube_inlet.flow_mol[t] *
                  fs.ng_preheater.tube_inlet.mole_frac_comp[(t, "C3H8")])
        F_C4H10 = (fs.ng_preheater.tube_inlet.flow_mol[t] *
                   fs.ng_preheater.tube_inlet.mole_frac_comp[(t, "C4H10")])
        O2_required = 2*F_CH4 + 3.5*F_C2H6 + 5*F_C3H8 + 6.5*F_C4H10

        O2_fed = (fs.air_compressor_s1.inlet.flow_mol[t] *
                  fs.air_compressor_s1.inlet.mole_frac_comp[(t, "O2")])
        return 1e-1*O2_fed == 1e-1*XS*O2_required

    # TODO - add constraint for cmb Q to in-bed heat exchanger (produces steam)
    @fs.Constraint(fs.time)
    def combustor_heat(b, t):
        return fs.cmb.heat_duty[t] == -1 * (fs.bhx2.heat_duty[t])


def add_result_constraints(fs):
    fs.soec_power = pyo.Var(fs.time,
                              units=pyo.units.MW)

    @fs.Constraint(fs.time)
    def soec_power_constraint(b, t):
        return b.soec_power[t] == pyo.units.convert(
            b.soec.total_power[t], pyo.units.MW)

    @fs.Expression(fs.time)
    def hydrogen_product_rate_expr(b, t):
        return (
            b.bhx1.shell_outlet.flow_mol[t]
            * b.bhx1.shell_outlet.mole_frac_comp[t, "H2"]
        )

    fs.hydrogen_product_rate = pyo.Var(fs.time,
                                         units=pyo.units.mol / pyo.units.s)

    @fs.Constraint(fs.time)
    def hydrogen_product_rate_constraint(b, t):
        return b.hydrogen_product_rate[t] == b.hydrogen_product_rate_expr[t]

    @fs.Expression(fs.time)
    def h2_product_rate_mass(b, t):  # kg/s
        return (fs.hydrogen_product_rate[t] *
                0.002 * pyo.units.kg / pyo.units.mol)

    @fs.Expression(fs.time)
    def soec_power_per_h2(b, t):
        return (
            b.soec.total_power[t]
            / b.hydrogen_product_rate_expr[t]
            / (0.002 * pyo.units.kg / pyo.units.mol)
        )

    @fs.Expression(fs.time)
    def h2_compressor_power_expr(b, t):
        return (
            fs.hcmp01.control_volume.work[t]
            + fs.hcmp02.control_volume.work[t]
            + fs.hcmp03.control_volume.work[t]
            + fs.hcmp04.control_volume.work[t]
        )

    fs.h2_compressor_power = pyo.Var(fs.time,
                                       units=pyo.units.MW)

    @fs.Constraint(fs.time)
    def h2_compressor_power_constraint(b, t):
        return b.h2_compressor_power[t] == pyo.units.convert(
            b.h2_compressor_power_expr[t], pyo.units.MW)

    # total heat supplied to HRSG
    fs.HRSG_heat_duty = pyo.Var(fs.time,
                                  initialize=300, units=pyo.units.MW)

    @fs.Constraint(fs.time)
    def HRSG_heat_duty_constraint(fs, t):
        return (
            -1*fs.HRSG_heat_duty[t] ==
            pyo.units.convert(fs.soec_air_HRSG_1.heat_duty[t], pyo.units.MW) +
            pyo.units.convert(fs.soec_air_HRSG_2.heat_duty[t], pyo.units.MW) +
            pyo.units.convert(fs.fluegas_HRSG.heat_duty[t], pyo.units.MW)
            )

    # steam required for ASU
    fs.ASU_HP_steam_heat = pyo.Var(fs.time,
                                     initialize=4, units=pyo.units.MW)

    @fs.Constraint(fs.time)
    def ASU_HP_steam_constraint(fs, t):
        return (fs.ASU_HP_steam_heat[t] ==
                0.005381*pyo.units.MJ/pyo.units.mol *
                fs.ASU.O2_outlet.flow_mol[t])

    # HRSG heat applied toward steam cycle
    fs.steam_cycle_heat = pyo.Var(fs.time,
                                    initialize=300, units=pyo.units.MW)

    @fs.Constraint(fs.time)
    def steam_cycle_heat_constraint(fs, t):
        return (fs.steam_cycle_heat[t] == fs.HRSG_heat_duty[t]
                - fs.ASU_HP_steam_heat[t])

    # power generated by steam cycle
    fs.steam_cycle_efficiency = pyo.Param(initialize=0.381, mutable=True)
    fs.steam_cycle_power = pyo.Var(fs.time,
                                     initialize=100, units=pyo.units.MW)
    fs.steam_cycle_loss = pyo.Var(fs.time,
                                    initialize=200, units=pyo.units.MW)

    @fs.Constraint(fs.time)
    def steam_cycle_power_constraint(fs, t):
        return (fs.steam_cycle_power[t] ==
                fs.steam_cycle_heat[t]*fs.steam_cycle_efficiency)

    @fs.Constraint(fs.time)
    def steam_cycle_loss_constraint(fs, t):
        return (fs.steam_cycle_loss[t] ==
                fs.steam_cycle_heat[t]*(1 - fs.steam_cycle_efficiency))

    # # gross plant power
    # fs.gross_power = pyo.Var(fs.time, initialize=670, units=pyo.units.MW)

    # @fs.Constraint(fs.time)
    # def gross_power_constraint(fs, t):
    #     return (fs.gross_power[t] == fs.stack_power_AC[t] +
    #             fs.steam_cycle_power[t])

    # boiler feedwater pump load - steam cycle not modeled, scaled from BB
    fs.feedwater_pump_work = pyo.Var(fs.time,
                                       initialize=2, units=pyo.units.MW)

    @fs.Constraint(fs.time)
    def feedwater_pump_work_constraint(fs, t):
        ref_steam_cycle_power = 262.8*pyo.units.MW
        ref_pump_work = 4.83*pyo.units.MW
        return (fs.feedwater_pump_work[t]*ref_steam_cycle_power ==
                ref_pump_work*fs.steam_cycle_power[t])

    # condensate pump load - steam cycle not modeled, scaled from BB
    fs.condensate_pump_work = pyo.Var(fs.time,
                                        initialize=0.1, units=pyo.units.MW)

    @fs.Constraint(fs.time)
    def condensate_pump_work_constraint(fs, t):
        ref_steam_cycle_power = 262.8*pyo.units.MW
        ref_pump_work = 0.15*pyo.units.MW
        return (fs.condensate_pump_work[t]*ref_steam_cycle_power ==
                ref_pump_work*fs.steam_cycle_power[t])

    # steam turbine auxiliary load - steam cycle not modeled, scaled from BB
    fs.steam_turbine_auxiliary = pyo.Var(fs.time,
                                           initialize=0.1, units=pyo.units.MW)

    @fs.Constraint(fs.time)
    def steam_turbine_auxiliary_constraint(fs, t):
        ref_steam_cycle_power = 262.8*pyo.units.MW
        ref_turbine_aux = 0.2*pyo.units.MW
        return (fs.steam_turbine_auxiliary[t]*ref_steam_cycle_power ==
                ref_turbine_aux*fs.steam_cycle_power[t])

    # miscellaneous BOP - scaled from BB based on NG mass flow
    fs.misc_BOP_load = pyo.Var(fs.time, initialize=0.5, units=pyo.units.MW)

    @fs.Constraint(fs.time)
    def misc_BOP_load_constraint(fs, t):
        NG_flow = pyo.units.convert(
            fs.ng_preheater.tube.properties_in[t].flow_mass,
            pyo.units.lb/pyo.units.hr)
        ref_NG_flow = 148095*pyo.units.lb/pyo.units.hr
        ref_load = 0.396*pyo.units.MW
        return fs.misc_BOP_load[t] == ref_load*NG_flow/ref_NG_flow

    # calculate cooling water flowrate
    fs.cooling_water_duty = pyo.Var(fs.time,
                                      initialize=450, units=pyo.units.MW)
    fs.cooling_water_flowrate = pyo.Var(fs.time, initialize=3000,
                                          units=pyo.units.lb/pyo.units.s)

    @fs.Constraint(fs.time)
    def cooling_water_duty_constraint(fs, t):
        return (fs.cooling_water_duty[t] ==
                fs.steam_cycle_loss[t] +
                fs.CPU.heat_duty[t]/1e6*pyo.units.MW +
                7.327*pyo.units.MW +  # 25 MMbtu/hr of additional heat
                pyo.units.convert(
                    -1*fs.intercooler_s1.heat_duty[t] +
                    -1*fs.intercooler_s2.heat_duty[t] +
                    -1*fs.hcmp_ic01.heat_duty[t] +
                    -1*fs.hcmp_ic02.heat_duty[t] +
                    -1*fs.hcmp_ic03.heat_duty[t] +
                    -1*fs.hcmp_ic04.heat_duty[t],
                    pyo.units.MW) +
                pyo.units.convert(
                    -1*fs.condenser.heat_duty[t],
                    pyo.units.MW)
                )

    @fs.Constraint(fs.time)
    def cooling_water_flowrate_constraint(fs, t):
        heat_capacity = 75.38*pyo.units.J/pyo.units.mol/pyo.units.K
        delta_T = 36*pyo.units.K
        molar_mass = 18*pyo.units.g/pyo.units.mol
        return (fs.cooling_water_flowrate[t] ==
                pyo.units.convert(
                    fs.cooling_water_duty[t]*molar_mass/heat_capacity/delta_T,
                    pyo.units.lb/pyo.units.s))

    # circulating pump work - not modeled, scaled from NGFC pathways study
    fs.circulating_pump_work = pyo.Var(fs.time,
                                         initialize=1, units=pyo.units.MW)

    @fs.Constraint(fs.time)
    def circulating_pump_work_constraint(fs, t):
        ref_flowrate = 18602.6*pyo.units.lb/pyo.units.s
        ref_load = 2.778*pyo.units.MW
        return (fs.circulating_pump_work[t]*ref_flowrate ==
                ref_load*fs.cooling_water_flowrate[t])

    # cooling tower fan load - not modeled, scaled from NGFC pathways study
    fs.cooling_tower_load = pyo.Var(fs.time,
                                      initialize=1, units=pyo.units.MW)

    @fs.Constraint(fs.time)
    def cooling_tower_load_constraint(fs, t):
        ref_flowrate = 18602.6*pyo.units.lb/pyo.units.s
        ref_load = 1.4372*pyo.units.MW
        return (fs.cooling_tower_load[t]*ref_flowrate ==
                ref_load*fs.cooling_water_flowrate[t])

    # auxiliary load of the plant
    fs.auxiliary_load = pyo.Var(fs.time, initialize=60, units=pyo.units.MW)

    @fs.Constraint(fs.time)
    def auxiliary_load_constraint(fs, t):
        return (fs.auxiliary_load[t] == fs.CPU.work[t]/1e6 +
                fs.feedwater_pump_work[t] + fs.condensate_pump_work[t] +
                fs.steam_turbine_auxiliary[t] + fs.misc_BOP_load[t] +
                fs.circulating_pump_work[t] + fs.cooling_tower_load[t] +
                pyo.units.convert(
                   (fs.air_blower.work_mechanical[t]/0.95 +
                    fs.aux_boiler_feed_pump.work_mechanical[t]/0.95 +
                    fs.air_compressor_s1.work_mechanical[t]/0.96 +
                    fs.air_compressor_s2.work_mechanical[t]/0.96 +
                    fs.hcmp01.work_mechanical[t]/0.96 +
                    fs.hcmp02.work_mechanical[t]/0.96 +
                    fs.hcmp03.work_mechanical[t]/0.96 +
                    fs.hcmp04.work_mechanical[t]/0.96
                    ), pyo.units.MW))

    # TODO - should the soec have transformer losses i.e. no power sent to grid
    # transformer losses
    fs.transformer_losses = pyo.Var(fs.time,
                                      initialize=2, units=pyo.units.MW)

    @fs.Constraint(fs.time)
    def transformer_losses_constraint(fs, t):
        HV_aux = pyo.units.convert(
            fs.air_compressor_s1.work_mechanical[t]/0.96 +
            fs.air_compressor_s2.work_mechanical[t]/0.96 +
            fs.air_compressor_s2.work_mechanical[t]/0.96 +
            fs.hcmp01.work_mechanical[t]/0.96 +
            fs.hcmp02.work_mechanical[t]/0.96 +
            fs.hcmp03.work_mechanical[t]/0.96 +
            fs.hcmp04.work_mechanical[t]/0.96,
            pyo.units.MW) + fs.CPU.work[t]/1e6
        MV_aux = fs.auxiliary_load[t] - HV_aux
        LV_aux = MV_aux*0.15

        HV_gen_loss = ((fs.steam_cycle_power[t] - fs.auxiliary_load[t])
                       + HV_aux)*.003
        HV_loss = HV_aux*.003
        MV_loss = MV_aux*.005
        LV_loss = LV_aux*.005

        return (fs.transformer_losses[t] ==
                HV_gen_loss + HV_loss + MV_loss + LV_loss)

    # net plant power
    fs.net_power = pyo.Var(fs.time, initialize=660, units=pyo.units.MW)

    @fs.Constraint(fs.time)
    def net_power_constraint(fs, t):
        return (fs.net_power[t] == fs.steam_cycle_power[t] -
                fs.soec_power[t] -
                fs.auxiliary_load[t] -
                fs.transformer_losses[t])

    # HHV efficiency
    fs.HHV_efficiency = pyo.Var(fs.time, initialize=0.6)

    # TODO - HHV efficiency should include H2 produced
    @fs.Constraint(fs.time)
    def efficiency_rule(fs, t):
        NG_HHV = 908839.23*pyo.units.J/pyo.units.mol
        H2_HHV = 141.7e6*pyo.units.J/pyo.units.kg  # Ref: https://www.engineeringtoolbox.com/fuels-higher-calorific-values-d_169.html
        return (fs.HHV_efficiency[t] == (fs.net_power[t] +
                pyo.units.convert(
                    (H2_HHV * fs.h2_product_rate_mass[t]),
                    pyo.units.MW))
                /
                pyo.units.convert(
                    (NG_HHV * fs.ng_preheater.tube.properties_in[t].flow_mol),
                    pyo.units.MW))

    # CO2 captured in kg/hr
    fs.CO2_captured = pyo.Var(fs.time, initialize=300,
                                units=pyo.units.kg/pyo.units.hr)

    @fs.Constraint(fs.time)
    def CO2_captured_constraint(fs, t):
        mass_flow = (fs.CPU.pureco2.flow_mol[t] *
                     fs.CPU.pureco2.mole_frac_comp[t, 'CO2'] *
                     0.04401*pyo.units.kg/pyo.units.mol)
        return fs.CO2_captured[t] == pyo.units.convert(
                            mass_flow, pyo.units.kg/pyo.units.hr)
    # CO2 emissions in kg/hr
    fs.CO2_emissions = pyo.Var(fs.time, initialize=300,
                                 units=pyo.units.kg/pyo.units.hr)

    @fs.Constraint(fs.time)
    def CO2_emission_constraint(fs, t):
        mass_flow = (fs.CPU.vent.flow_mol[t] *
                     fs.CPU.vent.mole_frac_comp[t, 'CO2'] *
                     0.04401*pyo.units.kg/pyo.units.mol)
        return fs.CO2_emissions[t] == pyo.units.convert(
                            -1*mass_flow, pyo.units.kg/pyo.units.hr)
    # TODO - include results constraints for water usage/makeup water i.e.
    # TODO - cooling water loss, steam cycle loss, process water for H2 gen


def initialize_results(fs):

    variables = [
        fs.hydrogen_product_rate,
        fs.h2_compressor_power,
        fs.soec_power,
        fs.HRSG_heat_duty,
        fs.ASU_HP_steam_heat,
        fs.steam_cycle_heat,
        fs.steam_cycle_power,
        fs.steam_cycle_loss,
        fs.feedwater_pump_work,
        fs.condensate_pump_work,
        fs.steam_turbine_auxiliary,
        fs.misc_BOP_load,
        fs.cooling_water_duty,
        fs.cooling_water_flowrate,
        fs.circulating_pump_work,
        fs.cooling_tower_load,
        fs.auxiliary_load,
        fs.transformer_losses,
        fs.net_power,
        fs.HHV_efficiency,
        fs.CO2_captured,
        fs.CO2_emissions]

    constraints = [
        fs.hydrogen_product_rate_constraint,
        fs.h2_compressor_power_constraint,
        fs.soec_power_constraint,
        fs.HRSG_heat_duty_constraint,
        fs.ASU_HP_steam_constraint,
        fs.steam_cycle_heat_constraint,
        fs.steam_cycle_power_constraint,
        fs.steam_cycle_loss_constraint,
        fs.feedwater_pump_work_constraint,
        fs.condensate_pump_work_constraint,
        fs.steam_turbine_auxiliary_constraint,
        fs.misc_BOP_load_constraint,
        fs.cooling_water_duty_constraint,
        fs.cooling_water_flowrate_constraint,
        fs.circulating_pump_work_constraint,
        fs.cooling_tower_load_constraint,
        fs.auxiliary_load_constraint,
        fs.transformer_losses_constraint,
        fs.net_power_constraint,
        fs.efficiency_rule,
        fs.CO2_captured_constraint,
        fs.CO2_emission_constraint]

    for v, c in zip(variables, constraints):
        for t in fs.time:
            calculate_variable_from_constraint(v[t], c[t])


def set_guess(fs):
    # Set guess for tear streams/ports (preheat_split.inlet, soec inlet)
    comp_guess = {  # air_prop is used as assumption as all fuel is combusted
        "O2": 0.2074,
        "H2O": 0.0099,
        "CO2": 0.0003,
        "N2": 0.7732,
        "Ar": 0.0092
    }

    _set_port(
        fs.preheat_split.inlet, F=6000, T=700, P=1.04e5,
        comp=comp_guess, fix=True
    )

    # Set guess for temp, pressure and mole frac conditions to initalize soec
    fs.soec.fc.temperature[:, 0].fix(1073.15)
    fs.soec.fc.pressure[:, 0].fix(1e5)
    fs.soec.fc.mole_frac_comp[:, 0, "H2O"].fix(0.90)
    fs.soec.fc.mole_frac_comp[:, 0, "H2"].fix(0.10)

    fs.soec.ac.temperature[:, 0].fix(1073.15)
    fs.soec.ac.pressure[:, 0].fix(1e5)
    fs.soec.ac.mole_frac_comp[:, 0, "O2"].fix(0.2074)
    fs.soec.ac.mole_frac_comp[:, 0, "H2O"].fix(0.0099)
    fs.soec.ac.mole_frac_comp[:, 0, "N2"].fix(0.7732)
    fs.soec.ac.mole_frac_comp[:, 0, "Ar"].fix(0.0092)
    fs.soec.ac.mole_frac_comp[:, 0, "CO2"].fix(0.0003)


def set_inputs(fs):
    # Set combustor outlet temperature and soec steam temperature
    fs.cmb_temperature.fix(1300)
    fs.soec_steam_temperature.fix(1073.15)

    # Set feed conditions of input fuel and utilities
    # TODO - note that flowrate vars for mxa1.air, air_compressor_s1.inlet,
    # aux_boiler_feed_pump.inlet, and fs.ng_preheater.tube_inlet are dof or
    # design specs that will be unfixed at final solve
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
    _set_port(  # TODO - relate NG feed to air feed requirement
        fs.air_compressor_s1.inlet, F=5000, T=330, P=1.01325e5,
        comp=air_comp, fix=True
    )
    _set_port(
        fs.ng_preheater.tube_inlet, F=380, T=330, P=1.04e5,
        comp=ng_comp, fix=True
    )
    _set_port(  # TODO - replaced mxa1 with air blower inlet
        fs.air_blower.inlet, F=3000, T=330, P=1.01325e5,
        comp=air_comp, fix=True
    )
    _set_port(  # TODO - remove alongside mxa1 as recycle is no longer spec'd
        fs.mxa1.recycle, F=1e-3, T=1073.15, P=1.04e5,  # Flow set to tiny val
        comp={"O2": 0.2074, "H2O": 0.0099, "CO2": 0.0003,
              "N2": 0.7732, "Ar": 0.0092}, fix=True
    )

    fs.aux_boiler_feed_pump.inlet.flow_mol.fix(2000)
    fs.aux_boiler_feed_pump.inlet.enth_mol.fix(
        iapws95.htpx(T=310*pyo.units.K, P=101325*pyo.units.Pa)
    )
    fs.aux_boiler_feed_pump.inlet.pressure.fix(101325)

    # Set known fixed (assumed) conditions for soec
    # TODO - increased n_cells to meet 600MW usage (equiv prod from sofc stack)
    fs.soec.n_cells.fix(400e6)
    fs.soec_heat_duty.fix(0)  # going for the thermoneutral point here

    fs.soec.fc.mole_frac_comp[:, 0, "H2"].fix(0.10)  # why not unfixed later
    fs.soec.fc.flow_mol[:, 0].fix(1e-5)
    fs.soec.ac.flow_mol[:, 0].fix(1e-5)


def initialize_plant(fs, solver):
    fs.preheat_split.initialize()
    iinit.propagate_state(fs.o04)
    iinit.propagate_state(fs.o05)

    fs.air_compressor_s1.initialize()
    iinit.propagate_state(fs.STAGE_1_OUT)
    fs.intercooler_s1.initialize()
    iinit.propagate_state(fs.IC_1_OUT)
    fs.air_compressor_s2.initialize()
    iinit.propagate_state(fs.STAGE_2_OUT)
    fs.intercooler_s2.initialize()
    iinit.propagate_state(fs.ba01)
    fs.ASU.initialize()
    iinit.propagate_state(fs.ba02)
    fs.oxygen_preheater.initialize()
    iinit.propagate_state(fs.cmb_mix_in)
    fs.pre_oxycombustor_translator.initialize()
    iinit.propagate_state(fs.ba03)

    fs.ng_preheater.initialize()
    iinit.propagate_state(fs.bng01)
    fs.cmb_mix.initialize()
    iinit.propagate_state(fs.bng02)
    fs.cmb.initialize()
    iinit.propagate_state(fs.cmb_out)

    fs.post_oxycombustor_translator.initialize()
    iinit.propagate_state(fs.fg01)
    fs.soec.initialize()
    iinit.propagate_state(fs.h01)
    iinit.propagate_state(fs.o01)
    fs.spltf1.initialize()
    iinit.propagate_state(fs.h02)
    iinit.propagate_state(fs.hr01)
    fs.splta1.initialize()
    iinit.propagate_state(fs.o02)

    fs.aux_boiler_feed_pump.initialize()
    iinit.propagate_state(fs.s01)
    fs.bhx1.initialize()
    iinit.propagate_state(fs.s02)
    fs.bhx2.initialize()
    iinit.propagate_state(fs.s03)

    fs.main_steam_split.initialize()
    iinit.propagate_state(fs.s10)

    fs.mxf1.initialize()
    iinit.propagate_state(fs.s12)

    fs.air_blower.initialize()
    iinit.propagate_state(fs.a01)
    fs.air_preheater_1.initialize()
    iinit.propagate_state(fs.a02)
    fs.air_preheater_2.initialize()
    iinit.propagate_state(fs.a03)
    iinit.propagate_state(fs.o03)
    fs.mxa1.initialize()
    iinit.propagate_state(fs.s13)

    # Unfix tear/guess streams
    fs.preheat_split.inlet.unfix()

    # Unfix some dof
    # TODO - maybe tie the unfixing of these vars to the add_constraint method
    fs.bhx2.outlet.enth_mol.unfix()  # related to soec_steam_temperature eqn
    fs.air_compressor_s1.inlet.flow_mol.unfix()  #  related to fixing soec.ac.flow_mol
    fs.ng_preheater.tube_inlet.flow_mol.unfix()  # related to fixing soec.fc.flow_mol
    fs.soec.E_cell.unfix()  # related to thermoneutral point i.e. soec_heat_duty.fix(0)

    # Solve soec model alone before flowsheet solve to check that it converges
    # solver.solve(fs.soec, tee=True, symbolic_solver_labels=True)
    # Unfix because cell is solved as thermoneutral (soec_heat_duty.fix(0))
    fs.s12_expanded.deactivate()
    fs.s13_expanded.deactivate()

    # # TODO - this seems to improve the convergence in some cases
    # iscale.constraint_autoscale_large_jac(m)

    # solver.solve(m, tee=True,
    #              # options={"max_iter":100},
    #              symbolic_solver_labels=True
    #              )

    # Connect soec to BOP - activate relavant arcs, unfix relevant variables
    fs.s12_expanded.activate()
    fs.s13_expanded.activate()

    fs.soec.fc.pressure[:, 0].unfix()
    fs.soec.fc.temperature[:, 0].unfix()
    fs.soec.fc.mole_frac_comp[:, 0, "H2O"].unfix()
    # TODO - fs.soec.fc.mole_frac_comp[:, 0, "H2"].unfix() - why not unfixed
    fs.soec.ac.pressure[:, 0].unfix()
    fs.soec.ac.temperature[:, 0].unfix()
    fs.soec.ac.mole_frac_comp[:, 0, "H2O"].unfix()
    fs.soec.ac.mole_frac_comp[:, 0, "O2"].unfix()
    fs.soec.ac.mole_frac_comp[:, 0, "N2"].unfix()
    fs.soec.ac.mole_frac_comp[:, 0, "Ar"].unfix()
    fs.soec.ac.mole_frac_comp[:, 0, "CO2"].unfix()

    fs.spltf1.split_fraction[:, "out"].unfix()

    # TODO - these unfixes are related to the fixed soec flow_mol for fc and ac
    # Unfix soec.fc.flow_mol if aux_boiler_feed_pump.inlet.flow_mol is fixed
    fs.aux_boiler_feed_pump.inlet.flow_mol.unfix()
    # Unfix soec.ac.flow_mol if mxa1.air.flow_mol is fixed
    fs.air_blower.inlet.flow_mol.unfix()

    solver.solve(fs, tee=True,
                 options={
                       "max_iter": 200,
                       "tol": 1e-7,
                       "bound_push": 1e-12,
                       "linear_solver": "ma57"
                           },
                 symbolic_solver_labels=True)


def initialize_bop(fs, solver):
    copy_port_values(source=fs.bhx1.shell_outlet,
                     destination=fs.hcmp_ic01.inlet)
    fs.hcmp_ic01.initialize()
    iinit.propagate_state(fs.h04)
    fs.hcmp01.initialize()
    iinit.propagate_state(fs.h05)
    fs.hcmp_ic02.initialize()
    iinit.propagate_state(fs.h06)
    fs.hcmp02.initialize()
    iinit.propagate_state(fs.h07)
    fs.hcmp_ic03.initialize()
    iinit.propagate_state(fs.h08)
    fs.hcmp03.initialize()
    iinit.propagate_state(fs.h09)
    fs.hcmp_ic04.initialize()
    iinit.propagate_state(fs.h10)
    fs.hcmp04.initialize()

    iinit.propagate_state(fs.a06)
    fs.soec_air_HRSG_1.initialize()
    iinit.propagate_state(fs.a07)
    fs.soec_air_HRSG_2.initialize()
    iinit.propagate_state(fs.fg02)
    fs.fluegas_HRSG.initialize()
    iinit.propagate_state(fs.fg03)
    fs.CPU_translator.initialize()

    iinit.propagate_state(fs.fg04)
    fs.condenser.initialize()
    iinit.propagate_state(fs.fg05)
    fs.flash.initialize()

    # TODO - copy_port method used instead of propagate_state because of issue
    # TODO - noticed in tags with propagating flash.vap_outlet to CPU.inlet
    # iinit.propagate_state(fs.fg06)
    copy_port_values(source=fs.flash.vap_outlet,
                     destination=fs.CPU.inlet)
    calculate_variable_from_constraint(fs.CPU.inlet.flow_mol[0],
                                       fs.CPU_inlet_F[0])
    calculate_variable_from_constraint(fs.CPU.inlet.pressure[0],
                                       fs.CPU_inlet_P[0])
    fs.CPU.initialize()

    # TODO - this seems to improve the convergence in some cases
    iscale.constraint_autoscale_large_jac(fs)

    solver.solve(fs, tee=True,
                 symbolic_solver_labels=True)


def additional_scaling(fs):
    for t, c in fs.oxygen_preheater.heat_transfer_equation.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in fs.ng_preheater.heat_transfer_equation.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in fs.bhx1.heat_transfer_equation.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in fs.air_preheater_1.heat_transfer_equation.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in fs.air_preheater_2.heat_transfer_equation.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)

    for t, c in fs.mxf1.enthalpy_mixing_equations.items():
        iscale.constraint_scaling_transform(c, 1e-6, overwrite=True)
    for t, c in fs.mxa1.enthalpy_mixing_equations.items():
        iscale.constraint_scaling_transform(c, 1e-6, overwrite=True)
    for t, c in fs.oxygen_preheater.unit_heat_balance.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in fs.ng_preheater.unit_heat_balance.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in fs.bhx1.unit_heat_balance.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in fs.air_preheater_1.unit_heat_balance.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in fs.air_preheater_2.unit_heat_balance.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)

    for t, c in fs.mxf1.fmxpress_eqn.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in fs.mxa1.amxpress_eqn.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in fs.cmb_mix.pressure_eqn.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)

    for t, c in fs.cmb_mix.enthalpy_mixing_equations.items():
        iscale.constraint_scaling_transform(c, 1e-6, overwrite=True)
    for t, c in fs.cmb.control_volume.enthalpy_balances.items():
        iscale.constraint_scaling_transform(c, 1e-6, overwrite=True)
    for t, c in fs.cmb.control_volume.properties_in[t].total_flow_balance.items():
        iscale.constraint_scaling_transform(c, 1e-1, overwrite=True)
    for t, c in fs.cmb.control_volume.properties_in[0].component_flow_balances.items():
        iscale.constraint_scaling_transform(c, 1e3, overwrite=True)
    for (t, j), c in fs.cmb.control_volume.material_balances.items():
        iscale.constraint_scaling_transform(c, 1e-1, overwrite=True)
    for (t, p, j), c in fs.cmb.control_volume.rate_reaction_stoichiometry_constraint.items():
        iscale.constraint_scaling_transform(c, 1e-1, overwrite=True)

    for t, c in fs.pre_oxycombustor_translator.pre_oxycombustor_translator_F.items():
        iscale.constraint_scaling_transform(c, 1e-1, overwrite=True)
    for t, c in fs.pre_oxycombustor_translator.pre_oxycombustor_translator_P.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)

    for t, c in fs.post_oxycombustor_translator.post_oxycombustor_translator_F.items():
        iscale.constraint_scaling_transform(c, 1e-1, overwrite=True)
    for t, c in fs.post_oxycombustor_translator.post_oxycombustor_translator_P.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)

    # This scaling term is needed to improve model convergence
    for t, c in fs.aux_boiler_feed_pump.eq_work.items():
        iscale.constraint_scaling_transform(c, 1e-10, overwrite=True)

    for (t, n), c in fs.soec.fc.flow_mol_eqn.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for (t, n, j), c in fs.soec.ac.mass_balance_eqn.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for (t, j), c in fs.mxa1.material_mixing_equations.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for (t, n, j), c in fs.soec.ac.mass_balance_eqn.items():
        iscale.constraint_scaling_transform(c, 1e-6, overwrite=True)

    for (t, r), v in fs.cmb.control_volume.rate_reaction_extent.items():
        iscale.set_scaling_factor(
            fs.cmb.control_volume.rate_reaction_extent, 1e-1)
    for (t, p, j), v in fs.cmb.control_volume.rate_reaction_generation.items():
        iscale.set_scaling_factor(
            fs.cmb.control_volume.rate_reaction_generation, 1e-1)

    iscale.calculate_scaling_factors(fs)


def additional_scaling_2(fs):
    iscale.set_scaling_factor(fs.soec_air_HRSG_1.control_volume.heat, 1e-5)
    iscale.set_scaling_factor(fs.soec_air_HRSG_2.control_volume.heat, 1e-5)
    iscale.set_scaling_factor(fs.fluegas_HRSG.control_volume.heat, 1e-5)
    iscale.set_scaling_factor(fs.condenser.control_volume.heat, 1e-4)
    iscale.set_scaling_factor(fs.flash.control_volume.heat, 1e-4)

    for t, c in fs.hcmp03.isentropic_pressure.items():
        iscale.constraint_scaling_transform(c, 1e-6, overwrite=True)

    for t, c in fs.condenser.control_volume.enthalpy_balances.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for t, c in fs.CPU_translator.CPU_translator_F.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for t, c in fs.CPU_translator.CPU_translator_P.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)

    for (p, q, j), c in fs.CPU_translator.properties_out[t].eq_mole_frac_tdew.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for (p, q, j), c in fs.CPU_translator.properties_out[0.0].equilibrium_constraint.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)

    for (p, q, j), c in fs.condenser.control_volume.properties_in[0.0].eq_mole_frac_tdew.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for (p, q, j), c in fs.condenser.control_volume.properties_out[0.0].eq_mole_frac_tdew.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for (p, q, j), c in fs.condenser.control_volume.properties_in[0.0].equilibrium_constraint.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for (p, q, j), c in fs.condenser.control_volume.properties_out[0.0].equilibrium_constraint.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)

    for (p, q, j), c in fs.flash.control_volume.properties_in[0.0].eq_mole_frac_tdew.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for (p, q, j), c in fs.flash.control_volume.properties_out[0.0].eq_mole_frac_tdew.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for (p, q, j), c in fs.flash.control_volume.properties_in[0.0].equilibrium_constraint.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for (p, q, j), c in fs.flash.control_volume.properties_out[0.0].equilibrium_constraint.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)

    for t, c in fs.fluegas_HRSG.control_volume.enthalpy_balances.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)

    for t, c in fs.CPU.pureco2_pressure_eq.items():
        iscale.constraint_scaling_transform(c, 1e-5, overwrite=True)
    for t, c in fs.CPU.water_pressure_eq.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)
    for t, c in fs.CPU.vent_pressure_eq.items():
        iscale.constraint_scaling_transform(c, 1e-3, overwrite=True)

    iscale.calculate_scaling_factors(fs)


def get_solver():
    use_idaes_solver_configuration_defaults()
    idaes.cfg.ipopt["options"]["nlp_scaling_method"] = "user-scaling"
    idaes.cfg.ipopt["options"]["tol"] = 1e-7
    # due to a lot of component mole fractions being on their lower bound of 0
    # bound push result in much longer solve times, so set it low.
    # idaes.cfg.ipopt["options"]["halt_on_ampl_error"] = 'yes'
    idaes.cfg.ipopt["options"]["bound_push"] = 1e-12
    idaes.cfg.ipopt["options"]["linear_solver"] = "ma57"
    idaes.cfg.ipopt["options"]["max_iter"] = 200
    # idaes.cfg.ipopt["options"]["ma27_pivtol"] = 1e-1
    # idaes.cfg.ipopt["options"]["ma57_pivtol"] = 1e-1
    # idaes.cfg.ipopt["options"]["ma57_pivtolmax"] = 1e-1
    return pyo.SolverFactory("ipopt")


def tag_inputs_opt_vars(fs):
    tags = iutil.ModelTagGroup()
    tags["single_cell_h2_side_inlet_flow"] = iutil.ModelTag(
        expr=fs.soec.fc.flow_mol[0, 0],
        format_string="{:.3f}",
        display_units=pyo.units.micromol/pyo.units.s,
        doc="Single cell H2 side inlet flow (feed)",
    )
    tags["single_cell_sweep_flow"] = iutil.ModelTag(
        expr=fs.soec.ac.flow_mol[0, 0],
        format_string="{:.3f}",
        display_units=pyo.units.micromol/pyo.units.s,
        doc="Single cell O2 side inlet flow (sweep)",
    )
    tags["feed_h2_frac"] = iutil.ModelTag(
        expr=fs.soec.fc.mole_frac_comp[0, 0, "H2"],
        format_string="{:.3f}",
        display_units=None,
        doc="H2 side inlet H2 mole frac (from recycle)",
    )
    tags["feed_water_flow"] = iutil.ModelTag(
        expr=fs.aux_boiler_feed_pump.inlet.flow_mol[0],
        format_string="{:.3f}",
        display_units=None,
        doc="Feed water flow (provides steam to soec fuel side)",
    )
    tags["sweep_air_flow"] = iutil.ModelTag(
        expr=fs.air_compressor_s1.inlet.flow_mol[0],
        format_string="{:.3f}",
        display_units=None,
        doc="Air flow to soec air side",
    )
    tags["single_cell_heat_required"] = iutil.ModelTag(
        expr=fs.soec_heat_duty[0],
        format_string="{:.3f}",
        display_units=pyo.units.W,
        doc="Heat duty of a singel SOEC cell",
    )
    tags["combustor_temperature"] = iutil.ModelTag(
        expr=fs.cmb.outlet.temperature[0],
        format_string="{:.3f}",
        display_units=pyo.units.K,
        doc="Combustor temperature",
    )
    tags["preheat_fg_split_to_air"] = iutil.ModelTag(
        expr=fs.preheat_split.split_fraction[0, "air"],
        format_string="{:.3f}",
        display_units=None,
        doc="Split frac. of soec air to air preheater, rest goes to NG heater",
    )
    tags["n_cells"] = iutil.ModelTag(
        expr=fs.soec.n_cells,
        format_string="{:,.0f}",
        display_units=None,
        doc="Number of SOEC cells",
    )
    tags["cell_pressure"] = iutil.ModelTag(
        expr=fs.aux_boiler_feed_pump.outlet.pressure[0],
        format_string="{:.3f}",
        display_units=pyo.units.kPa,
        doc="Steam and SOEC pressure",
    )
    # TODO - HEX areas aren't optimization variables as DT is fixed.
    tags["oxygen_preheater_area"] = iutil.ModelTag(
        expr=fs.oxygen_preheater.area,
        format_string="{:.3f}",
        display_units=pyo.units.m**2,
        doc="Air preheater area",
    )
    tags["ng_preheater_area"] = iutil.ModelTag(
        expr=fs.ng_preheater.area,
        format_string="{:.3f}",
        display_units=pyo.units.m**2,
        doc="NG preheater area",
    )
    tags["bhx1_area"] = iutil.ModelTag(
        expr=fs.bhx1.area,
        format_string="{:.3f}",
        display_units=pyo.units.m**2,
        doc="bhx1 area",
    )
    tags["air_preheater_1_area"] = iutil.ModelTag(
        expr=fs.air_preheater_1.area,
        format_string="{:.3f}",
        display_units=pyo.units.m**2,
        doc="air_preheater_1 area",
    )
    tags["air_preheater_2_area"] = iutil.ModelTag(
        expr=fs.air_preheater_2.area,
        format_string="{:.3f}",
        display_units=pyo.units.m**2,
        doc="air_preheater_2 area",
    )
    tags["hydrogen_product_rate"] = iutil.ModelTag(
        expr=fs.hydrogen_product_rate[0],
        format_string="{:.3f}",
        display_units=pyo.units.kmol / pyo.units.s,
        doc="Rate of hydrogen production",
    )
    fs.tag_input = tags
    display_input_tags(fs)


def tag_for_pfd_and_tables(fs):
    tag_group = iutil.ModelTagGroup()
    stream_states = tables.stream_states_dict(
        tables.arcs_to_stream_dict(
            fs,
            descend_into=False,
            additional={
                # TODO - include all inlet streams
                "a08": fs.soec_air_HRSG_1.outlet,
                "a09": fs.soec_air_HRSG_2.outlet,
                # "fg06": fs.flash.vap_outlet,
                "boiler feed water": fs.aux_boiler_feed_pump.inlet,
                "H2 product": fs.hcmp04.outlet,
                # "CO2 captured": fs.CPU.pureco2,
                # "CO2 vented": fs.CPU.vent,
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
        expr=fs.soec.total_power[0],
        format_string="{:.2f}",
        display_units=pyo.units.MW
    )
    tag_group["soec_n_cells"] = iutil.ModelTag(
        expr=fs.soec.n_cells,
        format_string="{:,.0f}",
        display_units=None
    )
    tag_group["E_cell"] = iutil.ModelTag(
        expr=fs.soec.E_cell[0],
        format_string="{:.4f}",
        display_units=pyo.units.V
    )

    tag_group["h2_product_rate_mass"] = iutil.ModelTag(
        expr=fs.h2_product_rate_mass[0],
        format_string="{:.3f}",
        display_units=pyo.units.kg / pyo.units.s,
    )
    # tag_group["co2_product_rate"] = iutil.ModelTag(
    #     expr=fs.co2_product_rate[0],
    #     format_string="{:.3f}",
    #     display_units=pyo.units.kmol / pyo.units.s,
    # )
    # tag_group["co2_product_rate_mass"] = iutil.ModelTag(
    #     expr=fs.co2_product_rate_mass[0],
    #     format_string="{:.3f}",
    #     display_units=pyo.units.kg / pyo.units.s,
    # )
    # tag_group["fuel_rate"] = iutil.ModelTag(
    #     expr=fs.ng_preheater.shell.properties_in[0].flow_mol,
    #     format_string="{:.3f}",
    #     display_units=pyo.units.kmol / pyo.units.s,
    # )
    # tag_group["fuel_rate_mass"] = iutil.ModelTag(
    #     expr=fs.ng_preheater.shell.properties_in[0].flow_mass,
    #     format_string="{:.3f}",
    #     display_units=pyo.units.kg / pyo.units.s,
    # )
    # tag_group["electric_power"] = iutil.ModelTag(
    #     expr=fs.electric_power[0],
    #     format_string="{:.3f}",
    #     display_units=pyo.units.MW,
    # )
    # tag_group["electric_power_per_mass_h2"] = iutil.ModelTag(
    #     expr=fs.electric_power_per_mass_h2[0],
    #     format_string="{:.3f}",
    #     display_units=pyo.units.MWh/pyo.units.kg,
    # )

    tag_group["status"] = iutil.ModelTag(expr=None, format_string="{}")
    fs.tag_pfd = tag_group


def display_input_tags(fs):
    # this is special for this model.  The input tags are not indexed
    print("")
    for key, tag in fs.tag_input.items():
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
    for i in iscale.extreme_jacobian_entries(
            jac=jac, nlp=nlp, small=1e-6, large=1e3):
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
    for v, sv in iscale.badly_scaled_var_generator(
            m, large=1e3, small=1e-4, zero=1e-12):
        print(f"    {v} -- {sv} -- {iscale.get_scaling_factor(v)}",
              file=sourceFile4)
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
        infilename = os.path.join(this_file_dir(), "rsofc_soec_mode.svg")
    with open(infilename, "r") as f:
        iutil.svg_tag(svg=f, tag_group=m.soec_fs.tag_pfd, outfile=filename)


def get_model(m=None, name="SOEC Module"):
    add_flowsheet(m)
    add_properties(m.soec_fs)

    add_asu(m.soec_fs)
    add_preheater(m.soec_fs)
    add_soec_air_side_units(m.soec_fs)
    add_combustor(m.soec_fs)
    add_aux_boiler_steam(m.soec_fs)
    add_soec_unit(m.soec_fs)
    add_more_hx_connections(m.soec_fs)
    add_soec_inlet_mix(m.soec_fs)

    add_design_constraints(m.soec_fs)
    expand_arcs = pyo.TransformationFactory("network.expand_arcs")
    expand_arcs.apply_to(m.soec_fs)
    set_guess(m.soec_fs)
    set_inputs(m.soec_fs)
    solver = get_solver()
    additional_scaling(m.soec_fs)
    iscale.scale_arc_constraints(m.soec_fs)

    initialize_plant(m.soec_fs, solver)
    add_h2_compressor(m.soec_fs)
    add_hrsg_and_cpu(m.soec_fs)
    expand_arcs = pyo.TransformationFactory("network.expand_arcs")
    expand_arcs.apply_to(m.soec_fs)
    additional_scaling_2(m.soec_fs)
    initialize_bop(m.soec_fs, solver)

    add_result_constraints(m.soec_fs)
    initialize_results(m.soec_fs)
    tag_inputs_opt_vars(m.soec_fs)
    tag_for_pfd_and_tables(m.soec_fs)

    display_input_tags(m.soec_fs)

    return m, solver

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
def variables_generator(
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
        if v.value is not None:
            yield v

if __name__ == "__main__":
    m = pyo.ConcreteModel()
    m, solver = get_model(m)
    # check_scaling(m)
    # write_pfd_results(m, "soec_init.svg")
  

    # Did not fix area of air preheater 2 as that it uses the outlet T as the dof instead
    # m.soec_fs.air_preheater_2.tube_outlet.temperature[0].unfix()
    m.soec_fs.air_preheater_2.overall_heat_transfer_coefficient.unfix()
    m.soec_fs.air_preheater_2.area.fix()

    # m.soec_fs.bhx1.delta_temperature_out.unfix()
    m.soec_fs.bhx1.overall_heat_transfer_coefficient.unfix()
    m.soec_fs.bhx1.area.fix()

    # m.soec_fs.air_preheater_1.delta_temperature_in.unfix()
    m.soec_fs.air_preheater_1.overall_heat_transfer_coefficient.unfix()
    m.soec_fs.air_preheater_1.area.fix()

    # m.soec_fs.oxygen_preheater.delta_temperature_in.unfix()
    m.soec_fs.oxygen_preheater.overall_heat_transfer_coefficient.unfix() 
    m.soec_fs.oxygen_preheater.area.fix()

    m.soec_fs.ng_preheater.overall_heat_transfer_coefficient.unfix()   
    m.soec_fs.ng_preheater.area.fix()

    print(f"Hydrogen product rate {m.soec_fs.tag_input['hydrogen_product_rate']}.")

    m.soec_fs.tag_input["hydrogen_product_rate"].fix()
    m.soec_fs.tag_input["single_cell_h2_side_inlet_flow"].unfix()
    solver.solve(m, tee=True)

    m.soec_fs.visualize("rSOEC Flowsheet")  # visualize flowsheet

    # rsofc_cost.get_rsofc_capital_costing(m)
    # rsofc_cost.lock_capital_cost(m)
    # rsofc_cost.get_rsofc_fixed_OM_costing(m, 650)
    # rsofc_cost.get_rsofc_soec_variable_OM_costing(m.soec_fs)

    # iscale.calculate_scaling_factors(m.fs.costing)
    # iscale.calculate_scaling_factors(m.soec_fs.costing)

    # solver.solve(m, tee=True)

    # m.fs.costing.display()
    # m.soec_fs.H2_costing.display()

    import sys
    sys.exit('stop')
    # strip_bounds = pyo.TransformationFactory("contrib.strip_var_bounds")
    # strip_bounds.apply_to(m, reversible=False)

    m.soec_fs.sweep_constraint = pyo.Constraint(
        expr=1e2* m.soec_fs.tag_input["single_cell_sweep_flow"].expression
        ==
        1e2* m.soec_fs.tag_input["single_cell_h2_side_inlet_flow"].expression
    )
    m.soec_fs.tag_input["single_cell_sweep_flow"].unfix()
    m.soec_fs.tag_input["sweep_air_flow"].unfix()
    m.soec_fs.tag_input["sweep_air_flow"].setlb(1000)
    m.soec_fs.tag_input["sweep_air_flow"].setub(20000)

    # m.soec_fs.tag_input["feed_water_flow"].unfix()  # redundant as already unfixed
    # m.soec_fs.tag_input["feed_water_flow"].setlb(1000)
    # m.soec_fs.tag_input["feed_water_flow"].setub(10000)

    # # TODO - is the sum of mole frac needed for this constraint?
    # m.soec_fs.tag_input["feed_h2_frac"].unfix()
    # m.soec_fs.tag_input["feed_h2_frac"].setlb(0.001)
    # m.soec_fs.tag_input["feed_h2_frac"].setub(0.98)

    # m.soec_fs.tag_input["preheat_fg_split_to_air"].unfix()
    # m.soec_fs.tag_input["preheat_fg_split_to_air"].setlb(0.85)
    # m.soec_fs.tag_input["preheat_fg_split_to_air"].setub(0.95)

    # m.soec_fs.tag_input["recover_split_to_hxh2"].unfix()
    # m.soec_fs.tag_input["recover_split_to_hxh2"].setlb(0.20)
    # m.soec_fs.tag_input["recover_split_to_hxh2"].setub(0.80)

    m.soec_fs.spltf1.split_fraction[:, "out"].setlb(0.001)
    m.soec_fs.spltf1.split_fraction[:, "out"].setub(0.999)

    # m.soec_fs.splta1.split_fraction[:, "out"].setlb(0.001)
    # m.soec_fs.splta1.split_fraction[:, "out"].setub(0.95)

    # m.soec_fs.aux_boiler_feed_pump.outlet.pressure.fix(20e5)
    # m.soec_fs.aux_boiler_feed_pump.outlet.pressure.setlb(1.1e5)
    # m.soec_fs.aux_boiler_feed_pump.outlet.pressure.setub(40e5)

    # m.soec_fs.obj = pyo.Objective(
    #     expr=m.soec_fs.ng_preheater.tube_inlet.flow_mol[0] / 10)
    # m.soec_fs.obj = pyo.Objective(expr=-m.soec_fs.hxh2.shell_outlet.mole_frac_comp[0, "H2"]*10)
    # m.soec_fs.obj = pyo.Objective(expr=m.soec_fs.H2_costing.total_variable_OM_cost[0])

    # m.soec_fs.tag_pfd["total_variable_OM_cost"] = iutil.ModelTag(
    #     expr=m.soec_fs.H2_costing.total_variable_OM_cost[0],
    #     format_string="{:.3f}",
    #     display_units=pyo.units.USD / pyo.units.kg,
    #     doc="Variable Hydrogen Production Cost",
    # )

    cols_input = ("hydrogen_product_rate",)
    cols_pfd = (
        "status",
        # "total_variable_OM_cost",
        "h2_product_rate_mass",
        # "co2_product_rate",
        # "co2_product_rate_mass",
        # "soec_power",
        # "h2_compressor_power",
        # "feed_pump_power",
        # "fuel_rate",
        # "fuel_rate_mass",
        # "electric_power",
        # "electric_power_per_mass_h2",
    )

    head_1 = m.soec_fs.tag_input.table_heading(tags=cols_input, units=True)
    head_2 = m.soec_fs.tag_pfd.table_heading(tags=cols_pfd, units=True)
    with open("opt_res.csv", "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(head_1 + head_2)
    for h in np.linspace(2.5, 0.1, 25):  # Generate 25 points btw 2.5 to 0.1
        m.soec_fs.tag_input["hydrogen_product_rate"].fix(
            float(h) * pyo.units.kmol / pyo.units.s
        )
        print(f"Hydrogen product rate {m.soec_fs.tag_input['hydrogen_product_rate']}.")
        # idaes.cfg.ipopt["options"]["max_iter"] = 5
        res = solver.solve(m.soec_fs, tee=True, symbolic_solver_labels=True,
                           options={"max_iter":100})

        # generator=variables_generator(m.soec_fs.soec)
        # # generator2=variables_near_bounds_generator(m, tol=1e-6)
               
    
        # sourceFile = open('model_variables.txt', 'w')          
        # for v in generator:
        #     print("var_name: %a, var_lb: %a, var_value: %a, var_ub: %a" %(
        #             v.name, pyo.value(v.lb), v.value, pyo.value(v.ub)),
        #             file=sourceFile)  
        # sourceFile.close()



        # # m.soec_fs.soec.fe.mole_frac_comp.display()
        # # m.soec_fs.soec.ae.mole_frac_comp.display()
        # residual_checker(m.soec_fs.soec)
        # import sys
        # sys.exit('stop')
        
        stat = idaeslog.condition(res)
        m.soec_fs.tag_pfd["status"].set(stat)
        # soec_cost.display_soec_costing(m)
        row_1 = m.soec_fs.tag_input.table_row(tags=cols_input, numeric=True)
        row_2 = m.soec_fs.tag_pfd.table_row(tags=cols_pfd, numeric=True)
        with open("opt_res.csv", "a", newline="") as f:
            w = csv.writer(f)
            w.writerow(row_1 + row_2)
        # write_pfd_results(
        #     m, f"soec_{m.tag_input['hydrogen_product_rate'].display(units=False)}.svg"
        # )