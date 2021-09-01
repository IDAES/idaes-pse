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


def add_properties(m):
    m.fs.o2_side_prop = GenericParameterBlock(
        default=get_prop(components={"O2", "H2O"}, phases=["Vap"], eos=EosType.IDEAL)
    )
    m.fs.h2_side_prop = GenericParameterBlock(
        default=get_prop(components={"H2", "H2O"}, phases=["Vap"], eos=EosType.IDEAL)
    )
    m.fs.o2l_side_prop = GenericParameterBlock(
        default=get_prop(components={"O2", "H2O"}, phases=["Vap", "Liq"], eos=EosType.IDEAL)
    )
    m.fs.h2l_side_prop = GenericParameterBlock(
        default=get_prop(components={"H2", "H2O"}, phases=["Vap", "Liq"], eos=EosType.IDEAL)
    )
    m.fs.fg_side_prop = GenericParameterBlock(
        default=get_prop(components=["N2", "O2", "CO2", "H2O"], phases=["Vap"], eos=EosType.IDEAL)
    )
    m.fs.water_prop = iapws95.Iapws95ParameterBlock()
    m.fs.h2_side_prop.set_default_scaling("mole_frac_comp", 10)
    m.fs.h2_side_prop.set_default_scaling("mole_frac_phase_comp", 10)
    m.fs.h2_side_prop.set_default_scaling("enth_mol_phase", 1e-3)

    m.fs.o2_side_prop.set_default_scaling("mole_frac_comp", 10)
    m.fs.o2_side_prop.set_default_scaling("mole_frac_phase_comp", 10)
    m.fs.o2_side_prop.set_default_scaling("enth_mol_phase", 1e-3)

    m.fs.h2l_side_prop.set_default_scaling("mole_frac_comp", 10)
    m.fs.h2l_side_prop.set_default_scaling("mole_frac_phase_comp", 10)
    m.fs.h2l_side_prop.set_default_scaling("enth_mol_phase", 1e-3)

    m.fs.o2l_side_prop.set_default_scaling("mole_frac_comp", 10)
    m.fs.o2l_side_prop.set_default_scaling("mole_frac_phase_comp", 10)
    m.fs.o2l_side_prop.set_default_scaling("enth_mol_phase", 1e-3)

    m.fs.fg_side_prop.set_default_scaling("mole_frac_comp", 10)
    m.fs.fg_side_prop.set_default_scaling("mole_frac_phase_comp", 10)
    m.fs.fg_side_prop.set_default_scaling("enth_mol_phase", 1e-3)


def add_soec(m):
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
            "property_package": m.fs.h2_side_prop,
            "outlet_list": ["out", "recycle"],
        }
    )
    m.fs.splta1 = gum.Separator(
        default={
            "property_package": m.fs.o2_side_prop,
            "outlet_list": ["out", "recycle"],
        }
    )

    m.fs.soec.E_cell.fix(1.28)  # unfix after initialize
    m.fs.soec.el.thickness.fix(8e-6)
    m.fs.soec.fe.thickness.fix(1e-3)
    m.fs.soec.ae.thickness.fix(20e-6)
    m.fs.soec.length.fix(0.05)
    m.fs.soec.width.fix(0.05)
    m.fs.soec.k_ae.fix(26.1e7)
    m.fs.soec.eact_ae.fix(120000)
    m.fs.soec.k_fe.fix(1.35e10)
    m.fs.soec.eact_fe.fix(110000)
    m.fs.soec.fe.k_res.fix(2.98e-5)
    m.fs.soec.fe.E_res.fix(-1392)
    m.fs.soec.ae.k_res.fix(8.114e-5)
    m.fs.soec.ae.E_res.fix(600)
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


def add_recycle_inlet_mixers(m):
    m.fs.mxf1 = gum.Mixer(
        default={
            "property_package": m.fs.h2_side_prop,
            "inlet_list": ["water", "recycle"],
            "momentum_mixing_type": gum.MomentumMixingType.none,
        }
    )
    m.fs.mxa1 = gum.Mixer(
        default={
            "property_package": m.fs.o2_side_prop,
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


def add_heatexchangers(m):
    m.fs.hxf1 = gum.HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_lmtd2_callback,
            "shell": {"property_package": m.fs.o2l_side_prop},
            "tube": {"property_package": m.fs.water_prop},
        }
    )
    m.fs.hxa1 = gum.HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_lmtd2_callback,
            "shell": {"property_package": m.fs.h2l_side_prop},
            "tube": {"property_package": m.fs.water_prop},
        }
    )
    m.fs.hxf2 = gum.HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_lmtd2_callback,
            "shell": {"property_package": m.fs.fg_side_prop},
            "tube": {"property_package": m.fs.water_prop},
        }
    )
    m.fs.hxa2 = gum.HeatExchanger(
        default={
            "delta_temperature_callback": delta_temperature_lmtd2_callback,
            "shell": {"property_package": m.fs.fg_side_prop},
            "tube": {"property_package": m.fs.water_prop},
        }
    )
    # Add new ports the the heat exchanger outlet for connect the IAPWS property
    # pacakge for steam in tubes to muticomponent in the mixer
    m.fs.hxa2.tube_outlet_adapt = Port()
    m.fs.hxa2._temperature_tube_outlet_ref = pyo.Reference(
        m.fs.hxa2.tube.properties_out[:].temperature
    )

    @m.fs.hxa2.Expression(m.fs.time, m.fs.soec.ac.config.comp_list)
    def atube_mole_frac_expr(b, t, i):
        if i == "H2O":
            return 1
        else:
            return 0

    m.fs.hxa2.tube_outlet_adapt.add(m.fs.hxa2._flow_mol_tube_outlet_ref, "flow_mol")
    m.fs.hxa2.tube_outlet_adapt.add(m.fs.hxa2._pressure_tube_outlet_ref, "pressure")
    m.fs.hxa2.tube_outlet_adapt.add(
        m.fs.hxa2._temperature_tube_outlet_ref, "temperature"
    )
    m.fs.hxa2.tube_outlet_adapt.add(m.fs.hxa2.atube_mole_frac_expr, "mole_frac_comp")
    m.fs.hxf2.tube_outlet_adapt = Port()
    m.fs.hxf2._temperature_tube_outlet_ref = pyo.Reference(
        m.fs.hxf2.tube.properties_out[:].temperature
    )

    @m.fs.hxf2.Expression(m.fs.time, m.fs.soec.fc.config.comp_list)
    def ftube_mole_frac_expr(b, t, i):
        if i == "H2O":
            return 1
        else:
            return 0

    m.fs.hxf2.tube_outlet_adapt.add(m.fs.hxf2._flow_mol_tube_outlet_ref, "flow_mol")
    m.fs.hxf2.tube_outlet_adapt.add(m.fs.hxf2._pressure_tube_outlet_ref, "pressure")
    m.fs.hxf2.tube_outlet_adapt.add(
        m.fs.hxf2._temperature_tube_outlet_ref, "temperature"
    )
    m.fs.hxf2.tube_outlet_adapt.add(m.fs.hxf2.ftube_mole_frac_expr, "mole_frac_comp")


def add_constraints(m):
    @m.fs.Expression(m.fs.time)
    def total_soec_power_expr(b, t):
        return b.soec.power[t] * b.soec.n_cells

    m.fs.total_soec_power = pyo.Var(m.fs.time, initialize=-10e6)

    @m.fs.Constraint(m.fs.time)
    def total_soec_power_eqn(b, t):
        return b.total_soec_power[t] == b.total_soec_power_expr[t]

    te = ThermoExpr(blk=m.fs, parameters=m.fs.water_prop)

    @m.fs.Constraint(m.fs.time)
    def ftemp_in_eqn(b, t):
        # T = b.soec.inlet_fc.temperature[t]
        p = b.hxf2.tube.properties_out[t].pressure
        hhx = b.hxf2.tube.properties_out[t].enth_mol
        return hhx == te.h(T=1073.15, p=p, x=1)

    @m.fs.Constraint(m.fs.time)
    def atemp_in_eqn(b, t):
        # T = b.soec.inlet_ac.temperature[t]
        p = b.hxa2.tube.properties_out[t].pressure
        hhx = b.hxa2.tube.properties_out[t].enth_mol
        return hhx == te.h(T=1073.15, p=p, x=1)

    @m.fs.Constraint(m.fs.time)
    def heat_duty_zero_eqn(b, t):
        return b.soec.heat_duty[t] == 0

    @m.fs.Expression(m.fs.time)
    def h2_flow_out_expr(b, t):
        return (
            m.fs.hxa1.shell_outlet.flow_mol[t]
            * m.fs.hxa1.shell_outlet.mole_frac_comp[t, "H2"]
        )

    m.fs.h2_out_flow = pyo.Var(m.fs.time, initialize=1.0)

    @m.fs.Constraint(m.fs.time)
    def h2_flow_out_eqn(b, t):
        return m.fs.h2_out_flow[t] == m.fs.h2_flow_out_expr[t]

    @m.fs.Expression(m.fs.time)
    def power_per_h2_MW(b, t):
        return (
            m.fs.total_soec_power_expr[t]
            / 1e3
            / 1000
            / 0.002
            / m.fs.hxa1.shell_outlet.flow_mol[t]
            / m.fs.hxa1.shell_outlet.mole_frac_comp[t, "H2"]
        )

    @m.fs.Expression(m.fs.time)
    def heat_per_h2_MW(b, t):
        return (
            (m.fs.hxf2.heat_duty[0] + m.fs.hxa2.heat_duty[t])
            / 1e3
            / 1000
            / 0.002
            / m.fs.hxa1.shell_outlet.flow_mol[t]
            / m.fs.hxa1.shell_outlet.mole_frac_comp[t, "H2"]
        )


def add_arcs(m):
    m.fs.f02 = Arc(source=m.fs.hxf1.tube_outlet, destination=m.fs.hxf2.tube_inlet)
    m.fs.a02 = Arc(source=m.fs.hxa1.tube_outlet, destination=m.fs.hxa2.tube_inlet)
    m.fs.f03 = Arc(source=m.fs.hxf2.tube_outlet_adapt, destination=m.fs.mxf1.water)
    m.fs.a03 = Arc(source=m.fs.hxa2.tube_outlet_adapt, destination=m.fs.mxa1.water)
    m.fs.f04 = Arc(source=m.fs.mxf1.outlet, destination=m.fs.soec.inlet_fc_mult)
    m.fs.a04 = Arc(source=m.fs.mxa1.outlet, destination=m.fs.soec.inlet_ac_mult)
    m.fs.f05 = Arc(source=m.fs.soec.outlet_fc_mult, destination=m.fs.spltf1.inlet)
    m.fs.a05 = Arc(source=m.fs.soec.outlet_ac_mult, destination=m.fs.splta1.inlet)
    m.fs.f06 = Arc(source=m.fs.spltf1.out, destination=m.fs.hxa1.shell_inlet)
    m.fs.a06 = Arc(source=m.fs.splta1.out, destination=m.fs.hxf1.shell_inlet)
    m.fs.f08 = Arc(source=m.fs.spltf1.recycle, destination=m.fs.mxf1.recycle)
    m.fs.a08 = Arc(source=m.fs.splta1.recycle, destination=m.fs.mxa1.recycle)

    expand_arcs = pyo.TransformationFactory("network.expand_arcs")
    expand_arcs.apply_to(m)


def set_inputs(m):
    m.fs.soec.n_cells.fix(200e6)

    m.fs.soec.fc.pressure[:, 0].fix(1e5)
    m.fs.soec.fc.flow_mol[:, 0].fix(9e-5/10)
    m.fs.soec.fc.mole_frac_comp[:, 0, "H2O"].fix(0.95)
    m.fs.soec.fc.mole_frac_comp[:, 0, "H2"].fix(0.05)

    m.fs.soec.ac.pressure[:, 0].fix(1e5)
    m.fs.soec.ac.flow_mol[:, 0].fix(1e-4/10)
    m.fs.soec.ac.mole_frac_comp[:, 0, "O2"].fix(0.1)
    m.fs.soec.ac.mole_frac_comp[:, 0, "H2O"].fix(0.9)

    m.fs.spltf1.split_fraction[:, "out"].fix(0.95)
    m.fs.splta1.split_fraction[:, "out"].fix(0.85)

    m.fs.hxf1.tube_inlet.enth_mol.fix(3000)
    m.fs.hxf1.tube_inlet.pressure.fix(1e5)
    m.fs.hxf1.tube_inlet.flow_mol.fix(1700)

    m.fs.hxa1.tube_inlet.enth_mol.fix(3000)
    m.fs.hxa1.tube_inlet.pressure.fix(1e5)
    m.fs.hxa1.tube_inlet.flow_mol.fix(1600)

    m.fs.hxf2.shell_inlet.mole_frac_comp[:, "CO2"].fix(0.04)
    m.fs.hxf2.shell_inlet.mole_frac_comp[:, "H2O"].fix(0.09)
    m.fs.hxf2.shell_inlet.mole_frac_comp[:, "O2"].fix(0.11)
    m.fs.hxf2.shell_inlet.mole_frac_comp[:, "N2"].fix(0.76)
    m.fs.hxf2.shell_inlet.pressure.fix(1.03e5)
    m.fs.hxf2.shell_inlet.temperature.fix(1200)
    m.fs.hxf2.shell_inlet.flow_mol[0].value = 600

    m.fs.hxa2.shell_inlet.mole_frac_comp[:, "CO2"].fix(0.04)
    m.fs.hxa2.shell_inlet.mole_frac_comp[:, "H2O"].fix(0.09)
    m.fs.hxa2.shell_inlet.mole_frac_comp[:, "O2"].fix(0.11)
    m.fs.hxa2.shell_inlet.mole_frac_comp[:, "N2"].fix(0.76)
    m.fs.hxa2.shell_inlet.pressure.fix(1.03e5)
    m.fs.hxa2.shell_inlet.temperature.fix(1200)
    m.fs.hxa2.shell_inlet.flow_mol[0].value = 400

    # Fix basic heat exchanger dimensions
    m.fs.hxf1.overall_heat_transfer_coefficient.fix(100)
    m.fs.hxa1.overall_heat_transfer_coefficient.fix(100)
    m.fs.hxf2.overall_heat_transfer_coefficient.fix(100)
    m.fs.hxa2.overall_heat_transfer_coefficient.fix(100)
    m.fs.hxf1.area.fix(5000)
    m.fs.hxa1.area.fix(2000)
    m.fs.hxf2.area.fix(5000)
    m.fs.hxa2.area.fix(5000)


def do_scaling(m):
    # flow through a single cell is very small
    iscale.set_scaling_factor(m.fs.total_soec_power, 1e-7)
    iscale.set_scaling_factor(m.fs.hxf1.overall_heat_transfer_coefficient, 1e-2)
    iscale.set_scaling_factor(m.fs.hxa1.overall_heat_transfer_coefficient, 1e-2)
    iscale.set_scaling_factor(m.fs.hxf2.overall_heat_transfer_coefficient, 1e-2)
    iscale.set_scaling_factor(m.fs.hxa2.overall_heat_transfer_coefficient, 1e-2)
    iscale.set_scaling_factor(m.fs.hxf1.delta_temperature_in, 1e-2)
    iscale.set_scaling_factor(m.fs.hxa1.delta_temperature_in, 1e-2)
    iscale.set_scaling_factor(m.fs.hxf2.delta_temperature_in, 1e-2)
    iscale.set_scaling_factor(m.fs.hxa2.delta_temperature_in, 1e-2)
    iscale.set_scaling_factor(m.fs.hxf1.area, 1e-4)
    iscale.set_scaling_factor(m.fs.hxa1.area, 1e-4)
    iscale.set_scaling_factor(m.fs.hxf2.area, 1e-4)
    iscale.set_scaling_factor(m.fs.hxa2.area, 1e-4)
    iscale.set_scaling_factor(m.fs.hxa1.tube.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.hxa1.shell.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.hxa2.tube.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.hxa2.shell.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.hxf1.tube.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.hxf1.shell.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.hxf2.tube.heat, 1e-6)
    iscale.set_scaling_factor(m.fs.hxf2.shell.heat, 1e-6)

    iscale.constraint_scaling_transform(m.fs.ftemp_in_eqn[0], 1e-3)
    iscale.constraint_scaling_transform(m.fs.atemp_in_eqn[0], 1e-3)
    iscale.constraint_scaling_transform(m.fs.mxf1.fmxpress_eqn[0], 1e-5)
    iscale.constraint_scaling_transform(m.fs.mxa1.amxpress_eqn[0], 1e-5)




    for i, c in m.fs.total_soec_power_eqn.items():
        s = iscale.get_scaling_factor(m.fs.total_soec_power[i])
        iscale.constraint_scaling_transform(c, s)
    iscale.set_scaling_factor(m.fs.h2_out_flow, 1e-3)
    for i, c in m.fs.h2_flow_out_eqn.items():
        s = iscale.get_scaling_factor(m.fs.h2_out_flow[i])
        iscale.constraint_scaling_transform(c, s)



    iscale.calculate_scaling_factors(m)


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


def do_initialization(m, solver):

    m.fs.soec.initialize()
    iinit.propagate_state(arc=m.fs.f05)
    iinit.propagate_state(arc=m.fs.a05)
    m.fs.spltf1.initialize()
    m.fs.splta1.initialize()
    iinit.propagate_state(arc=m.fs.f06)
    iinit.propagate_state(arc=m.fs.a06)
    m.fs.hxf1.initialize(outlvl=idaeslog.DEBUG)
    m.fs.hxa1.initialize()
    iinit.propagate_state(arc=m.fs.f02)
    iinit.propagate_state(arc=m.fs.a02)
    m.fs.hxf2.initialize()
    m.fs.hxa2.initialize()
    iinit.propagate_state(arc=m.fs.f03)
    iinit.propagate_state(arc=m.fs.a03)
    iinit.propagate_state(arc=m.fs.f08)
    iinit.propagate_state(arc=m.fs.a08)
    m.fs.mxf1.initialize()
    m.fs.mxa1.initialize()
    m.fs.soec.E_cell.unfix()
    m.fs.f04_expanded.deactivate()
    m.fs.a04_expanded.deactivate()
    solver.solve(m, tee=True)

    m.fs.f04_expanded.activate()
    m.fs.f04_expanded.flow_mol_equality.deactivate()
    m.fs.soec.fc.pressure[:, 0].unfix()
    m.fs.soec.fc.temperature[:, 0].unfix()
    m.fs.soec.fc.mole_frac_comp[:, 0, :].unfix()
    solver.solve(m, tee=True)

    m.fs.f04_expanded.flow_mol_equality.activate()
    m.fs.hxf1.tube_inlet.flow_mol.unfix()
    solver.solve(m, tee=True)

    m.fs.a04_expanded.activate()
    m.fs.a04_expanded.flow_mol_equality.deactivate()
    m.fs.soec.ac.pressure[:, 0].unfix()
    m.fs.soec.ac.temperature[:, 0].unfix()
    m.fs.soec.ac.mole_frac_comp[:, 0, :].unfix()
    solver.solve(m, tee=True)

    m.fs.a04_expanded.flow_mol_equality.activate()
    m.fs.hxa1.tube_inlet.flow_mol.unfix()
    solver.solve(m, tee=True)


def tag_model(m):
    """
    Create two dictionaries, one is called tags with tags for keys and
    expressions for values. The second is called tag_format with tags for keys
    and a formatting string. The formating string controls the number format and
    can be used to add additional text, to report units for example. The two
    dictionaries are returned, and also attached to the model.

    Args:
        m: (ConcreteModel) model to tag

    Returns:
        (dict) tags, (dict) tag_format
    """
    tags = {}  # dict of with tag keys and expressions for their values
    tag_format = {}  # format string for the tags

    def new_tag(name, expr, format):
        # funcion to keep it more compact
        tags[name] = expr
        tag_format[name] = format

    # Create a dict with Arc name keys and state block values
    stream_states = tables.stream_states_dict(
        tables.arcs_to_stream_dict(
            m.fs,
            descend_into=False,
            additional={
                "fg01": m.fs.hxf2.shell_inlet,
                "fg02": m.fs.hxf2.shell_outlet,
                "fg03": m.fs.hxa2.shell_inlet,
                "fg04": m.fs.hxa2.shell_outlet,
                "f01": m.fs.hxf1.tube_inlet,
                "a01": m.fs.hxa1.tube_inlet,
                "a07": m.fs.hxf1.shell_outlet,
                "f07": m.fs.hxa1.shell_outlet,
            },
        )
    )
    for i, s in stream_states.items():  # create the tags for steam quantities
        new_tag(f"{i}_Fmol", expr=s.flow_mol / 1000, format="{:.3f} kmol/s")
        new_tag(f"{i}_F", expr=s.flow_mass, format="{:.3f} kg/s")
        new_tag(f"{i}_P", expr=s.pressure / 1000, format="{:.3f} kPa")
        new_tag(f"{i}_T", expr=s.temperature, format="{:.2f} K")
        try:
            new_tag(f"{i}_x", expr=s.phase_frac["Vap"] * 100, format="{:.1f}%")
        except (KeyError, AttributeError):
            pass
        try:
            for c in s.mole_frac_comp:
                new_tag(f"{i}_y{c}", expr=s.mole_frac_comp[c] * 100, format="{:.3f}%")
        except (KeyError, AttributeError):
            pass

    new_tag("soec_power", expr=m.fs.total_soec_power_expr[0] / 1e6, format="{:.2f} MW")
    new_tag("soec_n_cells", expr=m.fs.soec.n_cells, format="{:.2e}")
    new_tag("E_cell", expr=m.fs.soec.E_cell[0], format="{:.4f} V")
    new_tag(
        "soec_heat_duty",
        expr=m.fs.soec.heat_duty[0] * m.fs.soec.n_cells / 1e6,
        format="{:.4f} MW",
    )
    new_tag("hxf2_heat_duty", expr=m.fs.hxf2.heat_duty[0] / 1e6, format="{:.2f} MW")
    new_tag("hxa2_heat_duty", expr=m.fs.hxa2.heat_duty[0] / 1e6, format="{:.2f} MW")
    new_tag("power_per_h2", expr=m.fs.power_per_h2_MW[0], format="{:.2f} MJ/kg")
    new_tag("heat_per_h2", expr=m.fs.heat_per_h2_MW[0], format="{:.2f} MJ/kg")
    new_tag(
        "spltf1_recycle_frac",
        expr=m.fs.spltf1.split_fraction[0, "recycle"],
        format="{:.3f}",
    )
    new_tag(
        "splta1_recycle_frac",
        expr=m.fs.splta1.split_fraction[0, "recycle"],
        format="{:.3f}",
    )
    new_tag(
        "h2_flow_out",
        expr=m.fs.h2_flow_out_expr[0] / 1000,
        format="{:.3f} kmol/s",
    )
    new_tag(
        "h2_product_h2_pct",
        expr=m.fs.hxa1.shell_outlet.mole_frac_comp[0, "H2"]*100,
        format="{:.1f}%",
    )
    new_tag(
        "h2_product_h2o_pct",
        expr=m.fs.hxa1.shell_outlet.mole_frac_comp[0, "H2O"]*100,
        format="{:.1f}%",
    )

    m.tags = tags
    m.tag_format = tag_format
    return tags, tag_format


def write_csv_header(m, columns, filename="res.csv"):
    # Write CSV header
    with open(filename, "w", newline="") as f:
        cw = csv.writer(f)
        cw.writerow(columns)


def write_csv_row(m, columns, filename="res.csv"):
    # Write CSV header
    with open(filename, "a", newline="") as f:
        cw = csv.writer(f)
        cw.writerow([pyo.value(m.tags[k]) for k in columns])


def write_pfd_results(filename, infilename=None):
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
        iutil.svg_tag(svg=f, tags=m.tags, outfile=filename, tag_format=m.tag_format)


def plot_soec(m):
    z = [
        (m.fs.soec.zset[iz] + m.fs.soec.zset[iz - 1]) / 2.0 for iz in m.fs.soec.izset_cv
    ]
    v = {t: None for t in m.fs.time}
    for t in v:
        v[t] = [
            pyo.value(m.fs.soec.fc.mole_frac_comp[t, iz, "H2"])
            for iz in m.fs.soec.izset_cv
        ]

    for t in v:
        plt.plot(z, v[t], label=f"t = {t}")

    plt.title(
        f"Fuel Channel H$_2$ (E$_{{cell}}$ = {pyo.value(m.fs.soec.E_cell[0]):.4f}V)"
    )
    plt.xlabel("z/L")
    plt.ylabel("$x_{H_2}$")
    plt.legend()
    plt.show()


def get_model(m=None):
    if m is None:
        m = pyo.ConcreteModel("SOEC Module")
    if not hasattr(m, "fs"):
        m.fs = FlowsheetBlock(default={"dynamic": False})
    add_properties(m)
    add_soec(m)
    add_heatexchangers(m)
    add_recycle_inlet_mixers(m)
    add_arcs(m)
    add_constraints(m)

    set_inputs(m)
    do_scaling(m)
    solver = get_solver()
    do_initialization(m, solver)
    tag_model(m)

    return m, solver


if __name__ == "__main__":
    m, solver = get_model()
    # plot_soec(m)
    write_pfd_results("soec_init.svg")

    #xs = np.linspace(1, 1.75, 31)
    xs = np.linspace(1, 0.10, 36)
    ohf = pyo.value(m.fs.h2_out_flow[0])
    m.obj = pyo.Objective(expr=m.fs.heat_per_h2_MW[0] + m.fs.power_per_h2_MW[0])
    m.fs.spltf1.split_fraction[0, "out"].unfix()
    m.fs.spltf1.split_fraction[0, "out"].setlb(0.80)
    m.fs.spltf1.split_fraction[0, "out"].setub(0.98)
    m.fs.spltf1.split_fraction[0, "out"].fix(0.98)
    m.fs.splta1.split_fraction[0, "out"].unfix()
    m.fs.splta1.split_fraction[0, "out"].setlb(0.70)
    m.fs.splta1.split_fraction[0, "out"].setub(0.98)
    m.fs.soec.ac.flow_mol[:, 0].unfix()
    m.fs.sweep_flow_ineq = pyo.Constraint(
        expr=m.fs.soec.ac.flow_mol[0, 0] == 0.5 * m.fs.soec.fc.flow_mol[0, 0]
    )
    m.fs.soec.inlet_fc.flow_mol.unfix()

    cols = [
        "number",
        "stat",
        "h2_flow_out",
        "soec_power",
        "spltf1_recycle_frac",
        "splta1_recycle_frac",
        "power_per_h2",
        "heat_per_h2",
        "h2_product_h2_pct",
        "h2_product_h2o_pct",
    ]
    write_csv_header(m, cols)
    #iscale.constraint_autoscale_large_jac(m)
    #solver.options["halt_on_ampl_error"] = "yes"
    for x in xs:
        m.fs.h2_out_flow.fix(ohf * x)
        m.tags["number"] = x
        try:
            res = solver.solve(m, tee=True)
            write_pfd_results(f"soec_{x:.4f}.svg")
            stat = idaeslog.condition(res)
            m.tags["stat"] = stat
            write_csv_row(m, cols)
        except ValueError:
            m.tags["stat"] = "Solver Error"
            write_csv_row(m, ["number", "stat"])
        except pyomo.common.errors.ApplicationError:
            m.tags["stat"] = "Solver Halted on Eval Error"
            write_csv_row(m, ["number", "stat"])
    check_scaling = True
    if check_scaling:
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
