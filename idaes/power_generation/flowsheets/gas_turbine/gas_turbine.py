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

import pyomo.environ as pyo
from pyomo.network import Arc
from pyomo.common.fileutils import this_file_dir

from idaes.core import FlowsheetBlock
from idaes.generic_models.properties.core.generic.generic_property import (
    GenericParameterBlock)
from idaes.generic_models.properties.core.generic.generic_reaction import (
    GenericReactionParameterBlock)
import idaes.generic_models.unit_models as um # um = unit models
import idaes.core.util as iutil
import idaes.core.util.tables as tables
import idaes.core.util.scaling as iscale
import idaes.core.plugins
from idaes.power_generation.properties.natural_gas_PR import get_prop, get_rxn
from idaes.core.solvers import use_idaes_solver_configuration_defaults


def tag_model(m):
    """
    Create to dictionaries. One is called tags with tags for keys and expressions
    for values. The second is called tag_format with tags for keys and a
    formatting string.  The formating string controls the number format and can
    be used to add additional text, to report units for example. The two
    dictionaries are returned, and also attached to the model.

    Args:
        m: (ConcreteModel) model to tag

    Returns:
        (dict) tags, (dict) tag_format
    """
    tags = {} # dict of with tag keys and expressions for their values
    tag_format = {} # format string for the tags
    def new_tag(name, expr, format):
        # funcion to keep it more compact
        tags[name] = expr
        tag_format[name] = format
    # Create a dict with Arc name keys and state block values
    stream_states = tables.stream_states_dict(
        tables.arcs_to_stream_dict(m.fs, descend_into=False))
    for i, s in stream_states.items(): # create the tags for steam quantities
        new_tag(f"{i}_Fvol", expr=s.flow_vol, format="{:.1f} m^3/s")
        new_tag(f"{i}_Fmol", expr=s.flow_mol/1000, format="{:.3f} kmol/s")
        new_tag(f"{i}_F", expr=s.flow_mass, format="{:.3f} kg/s")
        new_tag(f"{i}_P", expr=s.pressure/1000, format="{:.3f} kPa")
        new_tag(f"{i}_T", expr=s.temperature, format="{:.2f} K")
        for c in s.mole_frac_comp:
            new_tag(f"{i}_y{c}", expr=s.mole_frac_comp[c]*100, format="{:.3f}%")
    # Tags for non-stream things
    new_tag(
        "cmp1_power",
        expr=m.fs.cmp1.control_volume.work[0]/1e6,
        format="{:.2f} MW")
    new_tag(
        "gts1_power",
        expr=m.fs.gts1.control_volume.work[0]/1e6,
        format="{:.2f} MW")
    new_tag(
        "gts2_power",
        expr=m.fs.gts2.control_volume.work[0]/1e6,
        format="{:.2f} MW")
    new_tag(
        "gts3_power",
        expr=m.fs.gts3.control_volume.work[0]/1e6,
        format="{:.2f} MW")
    new_tag(
        "cmp1_head_isen",
        expr=m.fs.cmp1.performance_curve.head_isentropic[0]/1000,
        format="{:.2f} kJ/kg")
    new_tag(
        "gts1_head_isen",
        expr=m.fs.gts1.performance_curve.head_isentropic[0]/1000,
        format="{:.2f} kJ/kg")
    new_tag(
        "gts2_head_isen",
        expr=m.fs.gts2.performance_curve.head_isentropic[0]/1000,
        format="{:.2f} kJ/kg")
    new_tag(
        "gts3_head_isen",
        expr=m.fs.gts3.performance_curve.head_isentropic[0]/1000,
        format="{:.2f} kJ/kg")
    new_tag(
        "cmp1_eff_isen",
        expr=m.fs.cmp1.efficiency_isentropic[0]*100,
        format="{:.2f}%")
    new_tag(
        "gts1_eff_isen",
        expr=m.fs.gts1.efficiency_isentropic[0]*100,
        format="{:.2f}%")
    new_tag(
        "gts2_eff_isen",
        expr=m.fs.gts2.efficiency_isentropic[0]*100,
        format="{:.2f}%")
    new_tag(
        "gts3_eff_isen",
        expr=m.fs.gts3.efficiency_isentropic[0]*100,
        format="{:.2f}%")
    new_tag(
        "valve01_opening",
        expr=m.fs.valve01.valve_opening[0]*100,
        format="{:.2f}%")
    new_tag(
        "valve02_opening",
        expr=m.fs.valve02.valve_opening[0]*100,
        format="{:.2f}%")
    new_tag(
        "valve03_opening",
        expr=m.fs.valve03.valve_opening[0]*100,
        format="{:.2f}%")
    new_tag(
        "gt_total_power",
        expr=m.fs.gt_power[0]/1e6,
        format="{:.2f} MW")
    m.tags = tags
    m.tag_format = tag_format
    return tags, tag_format


def write_pfd_results(filename, tags, tag_format, infilename=None):
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
        infilename = os.path.join(this_file_dir(), "gas_turbine.svg")
    with open(infilename, "r") as f:
        iutil.svg_tag(
            svg=f,
            tags=tags,
            outfile=filename,
            tag_format=tag_format)


def performance_curves(m, flow_scale=0.896):
    """
    Add head and efficenty curves to the turbine stages in a model.

    Args:
        m: model with expcted structure

    Returns:
        None
    """
    # The flow scale variable less you scale the turbine up or down for differnt
    # nominal flows and differnt full load power outputs.
    m.fs.performace_flow_scale = pyo.Var(initialize=flow_scale)
    m.fs.performace_flow_scale.fix()
    fscale = m.fs.performace_flow_scale

    # PYLINT-TODO the names "eff_isen_eqn" and "head_isen_eqn" are reused below to define other constraint functions,
    # causing pylint to report function-redefined errors
    # this likely does not actually cause issues at runtime,
    # but it could be worth to check anyway if the pylint errors can be addressed
    # e.g. by giving unique names to each of the affected functions
    # pylint: disable=function-redefined

    # Efficiency curves for three stages
    @m.fs.gts1.performance_curve.Constraint(m.fs.time)
    def eff_isen_eqn(b, t):
        f = fscale*b.parent_block().control_volume.properties_in[t].flow_vol
        return b.parent_block().efficiency_isentropic[t] == 1.02*(
            1.4469E-14*f**5 -
            6.3333E-11*f**4 +
            6.6179E-08*f**3 -
            3.1728E-05*f**2 +
            7.7846E-03*f +
            1.0724E-01)
    @m.fs.gts2.performance_curve.Constraint(m.fs.time)
    def eff_isen_eqn(b, t):
        f = fscale*b.parent_block().control_volume.properties_in[t].flow_vol
        return b.parent_block().efficiency_isentropic[t] == 1.02*(
            2.6599E-16*f**5 -
            2.5894E-12*f**4 +
            6.0174E-09*f**3 -
            6.4156E-06*f**2 +
            3.5005E-03*f +
            1.0724E-01)
    @m.fs.gts3.performance_curve.Constraint(m.fs.time)
    def eff_isen_eqn(b, t):
        f = fscale*b.parent_block().control_volume.properties_in[t].flow_vol
        return b.parent_block().efficiency_isentropic[t] == 1.02*(
            5.8407E-18*f**5 -
            1.2203E-13*f**4 +
            6.0863E-10*f**3 -
            1.3927E-06*f**2 +
            1.6310E-03*f +
            1.0724E-01)
    # Head curves for three stages
    @m.fs.gts1.performance_curve.Constraint(m.fs.time)
    def head_isen_eqn(b, t):
        f = pyo.log(
            fscale*b.parent_block().control_volume.properties_in[t].flow_vol)
        return b.head_isentropic[t] == -(
            -2085.1*f**3 + 38433*f**2 - 150764*f + 422313)
    @m.fs.gts2.performance_curve.Constraint(m.fs.time)
    def head_isen_eqn(b, t):
        f = pyo.log(
            fscale*b.parent_block().control_volume.properties_in[t].flow_vol)
        return b.head_isentropic[t] == -(
            -1676.3*f**3 + 34916*f**2 - 173801*f + 456957)
    @m.fs.gts3.performance_curve.Constraint(m.fs.time)
    def head_isen_eqn(b, t):
        f = pyo.log(
            fscale*b.parent_block().control_volume.properties_in[t].flow_vol)
        return b.head_isentropic[t] == -(
            -1373.6*f**3 + 31759*f**2 - 188528*f + 500520)


def main(
    comps,
    rxns,
    phases,
    air_comp,
    ng_comp,
    m=None,
    initialize=True,
    flow_scale=0.896):
    """Generate and initialize the gas turbine flowsheet. This method returns
    a model and solver. The flow_scale argument can be used to scale the turnine
    up or down with the same relative performace.

    Args:
        comps: (set of str) a set of componets
        rxns: (dict) reactions with reaction name key and key component for
            conversion calculations
        phases: (list) phases potentially present
        air_comp: (dict) mole-fraction air composition
        ng_comp: (dict) mole-fraction natural gas composition
        m: (ConcreteModel) model to add flowsheet too.  If None create a new one.
        initialize: (bool) If true initialize the flowsheet
        flow_scale: (float) flow scale multiplier to scale the performace curves
            up or down. Original nominal volumetric flow/new nominal volumetric
            flow (default=0.896).

    Returns:
        (ConcreteModel) model, (Solver) solver
    """
    expand_arcs = pyo.TransformationFactory("network.expand_arcs")
    #
    # Model, flowsheet and properties
    #
    if m is None:
        m = pyo.ConcreteModel("Gas Turbine Model")
    if not hasattr(m, "fs"):
        m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.gas_prop_params = GenericParameterBlock(default=get_prop(comps, phases))
    m.fs.gas_prop_params.set_default_scaling("mole_frac_comp", 10)
    m.fs.gas_prop_params.set_default_scaling("mole_frac_phase_comp", 10)
    _mf_scale = {
        "Ar":100,
        "O2":100,
        "H2S":1000,
        "SO2":1000,
        "H2":1000,
        "CO":1000,
        "C2H4":1000,
        "C2H6":1000,
        "C3H8":1000,
        "C4H10":1000,
        "CO2":1000}
    for c, s in _mf_scale.items():
        m.fs.gas_prop_params.set_default_scaling(
            "mole_frac_comp", s, index=c)
        m.fs.gas_prop_params.set_default_scaling(
            "mole_frac_phase_comp", s, index=("Vap", c))
    m.fs.gas_prop_params.set_default_scaling(
        "enth_mol_phase", 1e-3, index="Vap")

    m.fs.gas_combustion = GenericReactionParameterBlock(
        default=get_rxn(m.fs.gas_prop_params, rxns))
    # Variable for mole-fraction of O2 in the flue gas.  To initialize the model
    # will use a fixed value here.
    m.fs.cmbout_o2_mol_frac = pyo.Var(m.fs.time, initialize=0.1157)
    m.fs.cmbout_o2_mol_frac.fix()
    #
    # Unit models
    #
    # inlet/outlet blocks
    m.fs.feed_air1 = um.Feed(default={"property_package": m.fs.gas_prop_params})
    m.fs.feed_fuel1 = um.Feed(default={"property_package": m.fs.gas_prop_params})
    m.fs.exhaust_1 = um.Product(default={"property_package": m.fs.gas_prop_params})
    # Valve roughly approximating variable stator vanes to control air flow.
    m.fs.vsv = um.Valve(default={
        "valve_function_callback":um.ValveFunctionType.linear,
        "property_package": m.fs.gas_prop_params})
    # Comprssor for air
    m.fs.cmp1 = um.Compressor(default={
        "property_package": m.fs.gas_prop_params,
        "support_isentropic_performance_curves":True})
    # Blade cooling air splitter
    m.fs.splt1 = um.Separator(default={
        "property_package": m.fs.gas_prop_params,
        "outlet_list":["air04", "air05", "air07", "air09"]})
    # Three valves for blade cooling air.
    m.fs.valve01 = um.Valve(default={
        "valve_function_callback":um.ValveFunctionType.linear,
        "property_package": m.fs.gas_prop_params})
    m.fs.valve02 = um.Valve(default={
        "valve_function_callback":um.ValveFunctionType.linear,
        "property_package": m.fs.gas_prop_params})
    m.fs.valve03 = um.Valve(default={
        "valve_function_callback":um.ValveFunctionType.linear,
        "property_package": m.fs.gas_prop_params})
    # Mixer for fuel injection. Since this is a steady state model the pressure
    # drop on the fuel control valve is lumped into this mixer
    m.fs.inject1 = um.Mixer(default={
        "property_package": m.fs.gas_prop_params,
        "inlet_list":["gas", "air"],
        "momentum_mixing_type":um.MomentumMixingType.none})
    # Three blade cooling air mixers. The blade cooling isn't explicitly modeled
    # so these mainly just maintain the proper mass and enegy balance.
    m.fs.mx1 = um.Mixer(default={
        "property_package": m.fs.gas_prop_params,
        "inlet_list":["gas", "air"],
        "momentum_mixing_type":um.MomentumMixingType.equality})
    m.fs.mx2 = um.Mixer(default={
        "property_package": m.fs.gas_prop_params,
        "inlet_list":["gas", "air"],
        "momentum_mixing_type":um.MomentumMixingType.equality})
    m.fs.mx3 = um.Mixer(default={
        "property_package": m.fs.gas_prop_params,
        "inlet_list":["gas", "air"],
        "momentum_mixing_type":um.MomentumMixingType.equality})
    # Combustor, for now assuming complete reactions and no NOx
    m.fs.cmb1 = um.StoichiometricReactor(default={
        "property_package": m.fs.gas_prop_params,
        "reaction_package": m.fs.gas_combustion,
        "has_pressure_change": True})
    # Three gas turbine stages, with performance curves for all three added by
    # the performance_curves() fuctions.
    m.fs.gts1 = um.Turbine(default={
        "property_package": m.fs.gas_prop_params,
        "support_isentropic_performance_curves":True})
    m.fs.gts2 = um.Turbine(default={
        "property_package": m.fs.gas_prop_params,
        "support_isentropic_performance_curves":True})
    m.fs.gts3 = um.Turbine(default={
        "property_package": m.fs.gas_prop_params,
        "support_isentropic_performance_curves":True})
    performance_curves(m)
    #
    # Additional constraints/expressions
    #
    # to make this simpler, just deal with VSV valve-approximation pressure drop
    # and forget the opening
    m.fs.vsv.pressure_flow_equation.deactivate()
    # O2 fraction in combustor constraint (to set var)
    @m.fs.Constraint(m.fs.time)
    def o2_flow(b, t):
        return m.fs.cmb1.control_volume.properties_out[t].mole_frac_comp["O2"] \
            == m.fs.cmbout_o2_mol_frac[t]
    # Fuel injector pressure is the air pressure. This lumps in the gas valve
    @m.fs.inject1.Constraint(m.fs.time)
    def mxpress_eqn(b, t):
        return b.mixed_state[t].pressure == b.air_state[t].pressure
    # The pressure drop in the cumbustor is just a fixed 5%.
    @m.fs.cmb1.Constraint(m.fs.time)
    def pressure_drop_eqn(b, t):
        return (0.95 ==
            b.control_volume.properties_out[t].pressure /
            b.control_volume.properties_in[t].pressure)
    # Constraints for complete combustion, use key compoents and 100% conversion
    @m.fs.cmb1.Constraint(m.fs.time, rxns.keys())
    def reaction_extent(b, t, r):
        k = rxns[r]
        prp = b.control_volume.properties_in[t]
        stc = -m.fs.gas_combustion.rate_reaction_stoichiometry[r, "Vap", k]
        extent = b.rate_reaction_extent[t, r]
        return extent == prp.flow_mol*prp.mole_frac_comp[k]/stc
    # Calculate the total gas turnine gross power output.
    @m.fs.Expression(m.fs.time)
    def gt_power_expr(b, t):
        return (
            m.fs.cmp1.control_volume.work[t] +
            m.fs.gts1.control_volume.work[t] +
            m.fs.gts2.control_volume.work[t] +
            m.fs.gts3.control_volume.work[t])
    # Add a varable and constraint for gross power.  This allows fixing power
    # for simulations where a specific power output is desired.
    m.fs.gt_power = pyo.Var(m.fs.time)
    @m.fs.Constraint(m.fs.time)
    def gt_power_eqn(b, t):
        return b.gt_power[t] == b.gt_power_expr[t]
    #
    # Arcs
    #
    # Arc names are imoprtant here becuase they match the PFD and are used for
    # tagging stream quantities
    m.fs.fuel01 = Arc(source=m.fs.feed_fuel1.outlet, destination=m.fs.inject1.gas)
    m.fs.air01 = Arc(source=m.fs.feed_air1.outlet, destination=m.fs.vsv.inlet)
    m.fs.air02 = Arc(source=m.fs.vsv.outlet, destination=m.fs.cmp1.inlet)
    m.fs.air03 = Arc(source=m.fs.cmp1.outlet, destination=m.fs.splt1.inlet)
    m.fs.air04 = Arc(source=m.fs.splt1.air04, destination=m.fs.inject1.air)
    m.fs.air05 = Arc(source=m.fs.splt1.air05, destination=m.fs.valve01.inlet)
    m.fs.air06 = Arc(source=m.fs.valve01.outlet, destination=m.fs.mx1.air)
    m.fs.air07 = Arc(source=m.fs.splt1.air07, destination=m.fs.valve02.inlet)
    m.fs.air08 = Arc(source=m.fs.valve02.outlet, destination=m.fs.mx2.air)
    m.fs.air09 = Arc(source=m.fs.splt1.air09, destination=m.fs.valve03.inlet)
    m.fs.air10 = Arc(source=m.fs.valve03.outlet, destination=m.fs.mx3.air)
    m.fs.g01 = Arc(source=m.fs.inject1.outlet, destination=m.fs.cmb1.inlet)
    m.fs.g02 = Arc(source=m.fs.cmb1.outlet, destination=m.fs.gts1.inlet)
    m.fs.g03 = Arc(source=m.fs.gts1.outlet, destination=m.fs.mx1.gas)
    m.fs.g04 = Arc(source=m.fs.mx1.outlet, destination=m.fs.gts2.inlet)
    m.fs.g05 = Arc(source=m.fs.gts2.outlet, destination=m.fs.mx2.gas)
    m.fs.g06 = Arc(source=m.fs.mx2.outlet, destination=m.fs.gts3.inlet)
    m.fs.g07 = Arc(source=m.fs.gts3.outlet, destination=m.fs.mx3.gas)
    m.fs.g08 = Arc(source=m.fs.mx3.outlet, destination=m.fs.exhaust_1.inlet)
    # Create arc constraints
    expand_arcs.apply_to(m)
    #
    # Tag model
    #
    tags, tag_format = tag_model(m)
    #
    # Set some scaling
    #
    for i in ["air05", "air07", "air09"]:
        iscale.set_scaling_factor(m.fs.splt1.split_fraction[0.0, i], 100)
    iscale.set_scaling_factor(m.fs.valve01.control_volume.work, 1e-8)
    iscale.set_scaling_factor(m.fs.valve02.control_volume.work, 1e-8)
    iscale.set_scaling_factor(m.fs.valve03.control_volume.work, 1e-8)
    iscale.set_scaling_factor(m.fs.vsv.control_volume.work, 1e-8)
    iscale.set_scaling_factor(m.fs.cmp1.control_volume.work, 1e-8)
    iscale.set_scaling_factor(m.fs.gts1.control_volume.work, 1e-8)
    iscale.set_scaling_factor(m.fs.gts2.control_volume.work, 1e-8)
    iscale.set_scaling_factor(m.fs.gts3.control_volume.work, 1e-8)
    iscale.set_scaling_factor(
        m.fs.gts1.control_volume.properties_in[0].flow_mol, 1e-5)
    iscale.set_scaling_factor(
        m.fs.gts2.control_volume.properties_in[0].flow_mol, 1e-5)
    iscale.set_scaling_factor(
        m.fs.gts3.control_volume.properties_in[0].flow_mol, 1e-5)
    iscale.set_scaling_factor(
        m.fs.gts1.control_volume.properties_out[0].flow_mol, 1e-5)
    iscale.set_scaling_factor(
        m.fs.gts2.control_volume.properties_out[0].flow_mol, 1e-5)
    iscale.set_scaling_factor(
        m.fs.gts3.control_volume.properties_out[0].flow_mol, 1e-5)
    for i, v in m.fs.cmb1.control_volume.rate_reaction_extent.items():
        if i[1] == "ch4_cmb":
            iscale.set_scaling_factor(v, 1e-2)
        else:
            iscale.set_scaling_factor(v, 1)
    for i, c in m.fs.cmb1.reaction_extent.items():
            iscale.constraint_scaling_transform(c, 1e-2)
    for i, v in m.fs.cmb1.control_volume.rate_reaction_generation.items():
        if i[2] in ["O2", "CO2", "CH4", "H2O"]:
            iscale.set_scaling_factor(v, 0.001)
        else:
            iscale.set_scaling_factor(v, 0.1)
    for v in m.fs.cmp1.deltaP.values():
        iscale.set_scaling_factor(v, 1e-6)
    for v in m.fs.gts1.deltaP.values():
        iscale.set_scaling_factor(v, 1e-6)
    for v in m.fs.gts2.deltaP.values():
        iscale.set_scaling_factor(v, 1e-6)
    for v in m.fs.gts3.deltaP.values():
        iscale.set_scaling_factor(v, 1e-6)
    for v in m.fs.valve01.deltaP.values():
        iscale.set_scaling_factor(v, 1e-6)
    for v in m.fs.valve02.deltaP.values():
        iscale.set_scaling_factor(v, 1e-6)
    for v in m.fs.valve03.deltaP.values():
        iscale.set_scaling_factor(v, 1e-6)
    for c in m.fs.o2_flow.values():
        iscale.constraint_scaling_transform(c, 10)
    for c in m.fs.gt_power_eqn.values():
        iscale.constraint_scaling_transform(c, 1e-8)
    for c in m.fs.mx1.enthalpy_mixing_equations.values():
        iscale.constraint_scaling_transform(c, 1e-8)
    for c in m.fs.mx2.enthalpy_mixing_equations.values():
        iscale.constraint_scaling_transform(c, 1e-8)
    for c in m.fs.mx3.enthalpy_mixing_equations.values():
        iscale.constraint_scaling_transform(c, 1e-8)
    for c in m.fs.inject1.enthalpy_mixing_equations.values():
        iscale.constraint_scaling_transform(c, 1e-8)
    for c in m.fs.gts1.performance_curve.head_isen_eqn.values():
        iscale.constraint_scaling_transform(c, 1e-5)
    for c in m.fs.gts1.performance_curve.eff_isen_eqn.values():
        iscale.constraint_scaling_transform(c, 2)
    for c in m.fs.gts2.performance_curve.head_isen_eqn.values():
        iscale.constraint_scaling_transform(c, 1e-5)
    for c in m.fs.gts2.performance_curve.eff_isen_eqn.values():
        iscale.constraint_scaling_transform(c, 2)
    for c in m.fs.gts3.performance_curve.head_isen_eqn.values():
        iscale.constraint_scaling_transform(c, 1e-5)
    for c in m.fs.gts3.performance_curve.eff_isen_eqn.values():
        iscale.constraint_scaling_transform(c, 2)
    for c in m.fs.inject1.mxpress_eqn.values():
        iscale.constraint_scaling_transform(c, 1e-5)
    #
    # Calculate scaling factors/scaling transform of constraints
    #
    iscale.calculate_scaling_factors(m)
    #
    # Set basic model inputs for initialization
    #
    # compressor
    m.fs.vsv.deltaP.fix(-100)
    m.fs.cmp1.efficiency_isentropic.fix(0.9)
    m.fs.cmp1.ratioP.fix(19.075)
    # blabe cooling air valves use expected flow to calculate valve flow
    # coefficients, so here set expected flow and after init the opening will
    # be the fixed var.
    m.fs.valve01.control_volume.properties_in[0].flow_mol.fix(2315)
    m.fs.valve02.control_volume.properties_in[0].flow_mol.fix(509)
    m.fs.valve03.control_volume.properties_in[0].flow_mol.fix(354)
    m.fs.valve01.valve_opening.unfix()
    m.fs.valve02.valve_opening.unfix()
    m.fs.valve03.valve_opening.unfix()
    # Feed streams
    # Air at compressor inlet
    m.fs.feed_air1.flow_mol[:] = 34875.9
    m.fs.feed_air1.temperature.fix(288.15)
    m.fs.feed_air1.pressure.fix(101047)
    for i, v in air_comp.items():
        m.fs.feed_air1.mole_frac_comp[:,i].fix(v)
    # Gas at fuel injection
    m.fs.feed_fuel1.flow_mol.fix(1348.8)
    m.fs.feed_fuel1.temperature.fix(310.928)
    m.fs.feed_fuel1.pressure.fix(3.10264e6)
    for i, v in ng_comp.items():
        m.fs.feed_fuel1.mole_frac_comp[:,i].fix(v)
    #
    # Initialization
    #
    solver = pyo.SolverFactory('ipopt')
    if initialize:
        # feeds
        m.fs.feed_air1.initialize()
        m.fs.feed_fuel1.initialize()
        # compressor
        iutil.copy_port_values(m.fs.air01)
        m.fs.vsv.initialize()
        iutil.copy_port_values(m.fs.air02)
        m.fs.cmp1.initialize()
        # splitter
        iutil.copy_port_values(m.fs.air03)
        m.fs.splt1.split_fraction[0, "air05"].fix(0.0916985*0.73)
        m.fs.splt1.split_fraction[0, "air07"].fix(0.0916985*0.27*0.59)
        m.fs.splt1.split_fraction[0, "air09"].fix(0.0916985*0.27*0.41)
        m.fs.splt1.initialize()
        m.fs.splt1.split_fraction[0, "air05"].unfix()
        m.fs.splt1.split_fraction[0, "air07"].unfix()
        m.fs.splt1.split_fraction[0, "air09"].unfix()
        # inject
        iutil.copy_port_values(m.fs.air04)
        iutil.copy_port_values(m.fs.fuel01)
        m.fs.inject1.mixed_state[0].pressure = pyo.value(m.fs.inject1.air.pressure[0])
        m.fs.inject1.initialize()
        # combustor
        iutil.copy_port_values(m.fs.g01)
        m.fs.cmb1.initialize()
        # gas turbine stage 1
        iutil.copy_port_values(m.fs.g02)
        m.fs.gts1.ratioP[0] = 0.7
        m.fs.gts1.initialize()
        # blade cooling air valve01, and calculate a flow coefficent
        iutil.copy_port_values(m.fs.air05)
        m.fs.valve01.Cv = 2
        m.fs.valve01.Cv.unfix()
        m.fs.valve01.valve_opening.fix(0.85)
        m.fs.valve01.control_volume.properties_out[0].pressure.fix(
            pyo.value(m.fs.gts1.control_volume.properties_out[0].pressure))
        m.fs.valve01.initialize()
        m.fs.valve01.control_volume.properties_out[0].pressure.unfix()
        m.fs.valve01.Cv.fix()
        m.fs.valve01.valve_opening.unfix()
        # mixer 1
        iutil.copy_port_values(m.fs.air06)
        iutil.copy_port_values(m.fs.g03)
        m.fs.mx1.mixed_state[0].pressure = pyo.value(m.fs.mx1.gas.pressure[0])
        m.fs.mx1.initialize()
        # gas turbine stage 2
        iutil.copy_port_values(m.fs.g04)
        m.fs.gts2.ratioP[0] = 0.7
        m.fs.gts2.initialize()
        # blade cooling air valve02, and calculate a flow coefficent
        iutil.copy_port_values(m.fs.air07)
        m.fs.valve02.Cv = 2
        m.fs.valve02.Cv.unfix()
        m.fs.valve02.valve_opening.fix(0.85)
        m.fs.valve02.control_volume.properties_out[0].pressure.fix(
            pyo.value(m.fs.gts2.control_volume.properties_out[0].pressure))
        m.fs.valve02.initialize()
        m.fs.valve02.control_volume.properties_out[0].pressure.unfix()
        m.fs.valve02.Cv.fix()
        m.fs.valve02.valve_opening.unfix()
        # mixer 2
        iutil.copy_port_values(m.fs.air08)
        iutil.copy_port_values(m.fs.g05)
        m.fs.mx2.mixed_state[0].pressure = pyo.value(m.fs.mx2.gas.pressure[0])
        m.fs.mx2.initialize()
        # gas turbine stage 3
        iutil.copy_port_values(m.fs.g06)
        m.fs.gts3.ratioP[0] = 0.7
        m.fs.gts3.initialize()
        # blade cooling air valve03, and calculate a flow coefficent
        iutil.copy_port_values(m.fs.air09)
        m.fs.valve03.Cv = 2
        m.fs.valve03.Cv.unfix()
        m.fs.valve03.valve_opening.fix(0.85)
        m.fs.valve03.control_volume.properties_out[0].pressure.fix(
            pyo.value(m.fs.gts3.control_volume.properties_out[0].pressure))
        m.fs.valve03.initialize()
        m.fs.valve03.control_volume.properties_out[0].pressure.unfix()
        m.fs.valve03.Cv.fix()
        m.fs.valve03.valve_opening.unfix()
        # mixer 3
        iutil.copy_port_values(m.fs.air10)
        iutil.copy_port_values(m.fs.g07)
        m.fs.mx3.mixed_state[0].pressure = pyo.value(m.fs.mx3.gas.pressure[0])
        m.fs.mx3.initialize()    # product blocks
        # product blocks
        iutil.copy_port_values(m.fs.g08)
        m.fs.exhaust_1.initialize()
        # Solve
        solver.solve(m, tee=True)
    return m, solver

def run_full_load(m, solver):
    """ Set the model up for off-design pressure driven flow.  Run 480 MW case.

    Args:
        m: (ConcreteModel) model to run
        solver: (Solver)

    Returns:
        None
    """
    m.fs.feed_fuel1.flow_mol.unfix()
    # For initialization, flow was specified and and valve opening calculated
    m.fs.valve01.control_volume.properties_in[0].flow_mol.unfix()
    m.fs.valve02.control_volume.properties_in[0].flow_mol.unfix()
    m.fs.valve03.control_volume.properties_in[0].flow_mol.unfix()
    # deltaP will be whatever is needed to satisfy the power requirment
    m.fs.vsv.deltaP.unfix()
    # The compressor efficiency is a little high since it doesn't include
    # throttling in the valve use to approximate VSV.
    m.fs.cmp1.efficiency_isentropic.fix(0.92)
    m.fs.cmp1.ratioP.fix(17.5) # lowering this ratio, just means less pressure
                               # drop in the VSV valve, decresing throttle loss
    # Exhaust pressure will be a bit over ATM due to HRSG. This will come from
    # HRSG model when coupled to form NGCC model
    m.fs.gts3.control_volume.properties_out[0].pressure.fix(102594)
    # Don't know how much blade cooling air is needed for off desing case, but
    # full load flows were based on WVU model.  For now just leave valves at
    # fixed opening.
    m.fs.valve01.valve_opening.fix()
    m.fs.valve02.valve_opening.fix()
    m.fs.valve03.valve_opening.fix()
    # Feeds
    # Air at compressor inlet
    m.fs.feed_air1.temperature.fix(288.15)
    m.fs.feed_air1.pressure.fix(103421)
    # Gas at fuel injection
    m.fs.feed_fuel1.temperature.fix(310.928)
    m.fs.feed_fuel1.pressure.fix(3.10264e6)
    # It looks like the O2 or air/fuel ratio is actually determined to roughly
    # get a constant exhaust temperature.  In a real gas turbine there are
    # probably a lot of complex factors, but here we fix the exhaust temp.
    m.fs.cmbout_o2_mol_frac.unfix()
    m.fs.exhaust_1.temperature.fix(910)
    # Fix the gross power output
    m.fs.gt_power[0].fix(-480e6) # Watts negative because is power out
    # Solve
    return solver.solve(m, tee=True)

def run_series(m, solver):
    """Run off design at 480 MW then a series of lower loads. Call run_full_load
    first to configure the model properly before calling this function. Saves
    PFD and row in tabulated csv file results for each case.

    Args:
        m: (ConcreteModel) model to run
        solver: (Solver)

    Returns:
        None
    """
    tags = m.tags
    tag_format = m.tag_format
    # Write the results to a CSV with these columns
    columns=[
        "power", "Fuel_flow", "Air_flow", "TI_temp", "TI_pressure",
        "TI_flow", "EX_temp", "EX_pressure", "EX_flow", "FG_O2", "FG_CO2"]

    # Write CSV header
    with open("res.csv", "w", newline='') as f:
        cw = csv.writer(f)
        cw.writerow(columns)

    # Run simulations, write CSV rows, and write flowsheets
    for i, p in enumerate([480 - 48*j for j in range(9)]):
        print(f"Loop {p}") # Print power to track progress
        m.fs.gt_power[0].fix(-p*1e6) # fix gross power
        solver.solve(m, tee=True)
        with open("res.csv", "a", newline='') as f:
            csv.writer(f).writerow([
                p,
                pyo.value(tags["fuel01_F"]), # mass flow
                pyo.value(tags["air01_F"]), # mass flow
                pyo.value(tags["g02_T"]),
                pyo.value(tags["g02_P"]),
                pyo.value(tags["g02_F"]),
                pyo.value(tags["g08_T"]),
                pyo.value(tags["g08_P"]),
                pyo.value(tags["g08_F"]),
                pyo.value(tags["g08_yO2"]),
                pyo.value(tags["g08_yCO2"])])
        write_pfd_results(f"gas_turbine_results_{p}.svg", tags, tag_format)


if __name__ == "__main__":
    #
    # Solver config
    #
    # Use idaes config, because these will apply to all ipopt solvers created
    use_idaes_solver_configuration_defaults()
    idaes.cfg.ipopt["options"]["nlp_scaling_method"] = "user-scaling"
    # due to a lot of component mole fractions being on their lower bound of 0
    # bound push result in much longer solve times, so set it low.
    idaes.cfg.ipopt["options"]["bound_push"] = 1e-10


    comps = { # components present
        "CH4", "C2H6", "C2H4", "CO", "H2S", "H2", "O2", "H2O", "CO2", "N2",
        "Ar", "SO2"}
    rxns = { # reactions and key compoents for conversion
        "ch4_cmb":"CH4",
        "c2h6_cmb":"C2H6",
        "c2h4_cmb":"C2H4",
        "co_cmb":"CO",
        "h2s_cmb":"H2S",
        "h2_cmb":"H2"}
    phases = ["Vap"]
    air_comp = {
        "CH4":0.0,
        "C2H6":0.0,
        "C2H4":0.0,
        "CO":0.0,
        "H2S":0.0,
        "H2":0.0,
        "O2":0.2074,
        "H2O":0.0099,
        "CO2":0.0003,
        "N2":0.7732,
        "Ar":0.0092,
        "SO2":0.0}
    ng_comp = {
        "CH4":0.87,
        "C2H6":0.0846,
        "C2H4":0.0003,
        "CO":0.0009,
        "H2S":0.0004,
        "H2":0.0036,
        "O2":0.0007,
        "H2O":0.0,
        "CO2":0.0034,
        "N2":0.0361,
        "Ar":0.0,
        "SO2":0.0}

    m, solver = main(
        comps=comps,
        rxns=rxns,
        phases=phases,
        air_comp=air_comp,
        ng_comp=ng_comp)
    run_full_load(m, solver)
    #iscale.constraint_autoscale_large_jac(m)
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
    for v, sv in iscale.badly_scaled_var_generator(m, large=1e2, small=1e-2, zero=1e-12):
        print(f"    {v} -- {sv} -- {iscale.get_scaling_factor(v)}")
    print(f"Jacobian Condition Number: {iscale.jacobian_cond(jac=jac):.2e}")
    write_pfd_results("gas_turbine_results.svg", m.tags, m.tag_format)
    #run_series(m, solver)
