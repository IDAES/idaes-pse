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
    m.fs.fg01 = Arc(source=m.fs.cmb.outlet, destination=m.fs.bhx2.shell_inlet)
    m.fs.s02 = Arc(
        source=m.fs.aux_boiler_feed_pump.outlet,
        destination=m.fs.bhx1.tube_inlet
    )
    m.fs.fg02 = Arc(source=m.fs.bhx2.shell_outlet, destination=m.fs.bhx1.shell_inlet)

def add_soec_unit(m):
    m.fs.soec = pum.IsothermalSofc(
        default={
            "nz": 20,
            "nxfe": 6,
            "nxae": 6,
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
    m.fs.soec.el.thickness.fix(8e-6)
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
    m.fs.air_preheater.area.fix(1000)
    m.fs.air_preheater.overall_heat_transfer_coefficient.fix(100)
    m.fs.ng_preheater.area.fix(1000)
    m.fs.ng_preheater.overall_heat_transfer_coefficient.fix(100)
    m.fs.bhx2.area.fix(200)
    m.fs.bhx2.overall_heat_transfer_coefficient.fix(100)
    m.fs.bhx1.area.fix(200)
    m.fs.bhx1.overall_heat_transfer_coefficient.fix(100)


    m.fs.preheat_split.split_fraction[:, "air"].fix(0.5)

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
        m.fs.air_preheater.tube_inlet, F=600, T=330, P=1.04e5, comp=air_comp, fix=True
    )
    _set_port(
        m.fs.ng_preheater.tube_inlet, F=50, T=330, P=1.04e5, comp=ng_comp, fix=True
    )
    m.fs.bhx2.tube_inlet.flow_mol.fix(500)
    m.fs.bhx2.tube_inlet.enth_mol.fix(iapws95.htpx(T=950*pyo.units.K, P=20.6e5*pyo.units.Pa))
    m.fs.bhx2.tube_inlet.pressure.fix(20.6e5)

    m.fs.aux_boiler_feed_pump.inlet.flow_mol.fix(500)
    m.fs.aux_boiler_feed_pump.inlet.enth_mol.fix(
        iapws95.htpx(T=310*pyo.units.K, P=101325*pyo.units.Pa)
    )
    m.fs.aux_boiler_feed_pump.inlet.pressure.fix(101325)
    m.fs.aux_boiler_feed_pump.outlet.pressure.fix(20.6e5)
    m.fs.aux_boiler_feed_pump.efficiency_isentropic.fix(0.85)

    m.fs.spltf1.split_fraction[:, "out"].fix(0.95)
    m.fs.splta1.split_fraction[:, "out"].fix(0.85)
    m.fs.soec.n_cells.fix(200e6)

    m.fs.soec.fc.flow_mol[:, 0].fix(9e-5/10)
    m.fs.soec.fc.pressure[:, 0].fix(1e5)
    m.fs.soec.fc.mole_frac_comp[:, 0, "H2O"].fix(0.95)
    m.fs.soec.fc.mole_frac_comp[:, 0, "H2"].fix(0.05)

    m.fs.soec.ac.flow_mol[:, 0].fix(1e-4/10)
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


def tag_model(m):
    tag_group = iutil.ModelTagGroup()
    stream_states = tables.stream_states_dict(
        tables.arcs_to_stream_dict(
            m.fs,
            descend_into=False,
            additional={
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
                format_string="{:.1f}",
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

    m.tag_group = tag_group


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
        iutil.svg_tag(svg=f, tag_group=m.tag_group, outfile=filename)


def get_model(m=None, name="SOEC Module"):
    m = add_flowsheet(m)
    add_properties(m)
    add_preheater(m)
    add_combustor(m)
    add_aux_boiler_steam(m)
    add_soec_unit(m)
    expand_arcs = pyo.TransformationFactory("network.expand_arcs")
    expand_arcs.apply_to(m)

    tag_model(m)
    set_guess(m)
    set_inputs(m)
    solver = get_solver()
    do_initialize(m, solver)
    return m



if __name__ == "__main__":
    m = get_model()
    write_pfd_results(m, "soec_init.svg")
