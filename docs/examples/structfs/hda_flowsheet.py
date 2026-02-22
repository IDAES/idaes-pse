###############################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
###############################################################################


#
# # HDA Flowsheet Simulation and Optimization
#
# Author: Jaffer Ghouse
# Maintainer: Brandon Paul
# Updated: 2023-06-01
#
# ## Learning outcomes
#
#
# - Construct a steady-state flowsheet using the IDAES unit model library
# - Connecting unit models in a  flowsheet using Arcs
# - Using the SequentialDecomposition tool to initialize a flowsheet with recycle
# - Formulate and solve an optimization problem
#     - Defining an objective function
#     - Setting variable bounds
#     - Adding additional constraints
#
#
# ## Problem Statement
#
# Hydrodealkylation is a chemical reaction that often involves reacting
# an aromatic hydrocarbon in the presence of hydrogen gas to form a
# simpler aromatic hydrocarbon devoid of functional groups. In this
# example, toluene will be reacted with hydrogen gas at high temperatures
#  to form benzene via the following reaction:
#
# **C<sub>6</sub>H<sub>5</sub>CH<sub>3</sub> + H<sub>2</sub> → C<sub>6</sub>H<sub>6</sub> + CH<sub>4</sub>**
#
#
# This reaction is often accompanied by an equilibrium side reaction
# which forms diphenyl, which we will neglect for this example.
#
# This example is based on the 1967 AIChE Student Contest problem as
# present by Douglas, J.M., Chemical  Design of Chemical Processes, 1988,
# McGraw-Hill.
#
# The flowsheet that we will be using for this module is shown below with the stream conditions. We will be processing toluene and hydrogen to produce at least 370 TPY of benzene. As shown in the flowsheet, there are two flash tanks, F101 to separate out the non-condensibles and F102 to further separate the benzene-toluene mixture to improve the benzene purity.  Note that typically a distillation column is required to obtain high purity benzene but that is beyond the scope of this workshop. The non-condensibles separated out in F101 will be partially recycled back to M101 and the rest will be either purged or combusted for power generation.We will assume ideal gas for this flowsheet. The properties required for this module are available in the same directory:
#
# - hda_ideal_VLE.py
# - hda_reaction.py
#
# The state variables chosen for the property package are **flows of component by phase, temperature and pressure**. The components considered are: **toluene, hydrogen, benzene and methane**. Therefore, every stream has 8 flow variables, 1 temperature and 1 pressure variable.
#
# ![](HDA_flowsheet.png)
#
#

# ## Importing required pyomo and idaes components
#
#
# To construct a flowsheet, we will need several components from the pyomo and idaes package. Let us first import the following components from Pyomo:
# - Constraint (to write constraints)
# - Var (to declare variables)
# - ConcreteModel (to create the concrete model object)
# - Expression (to evaluate values as a function of variables defined in the model)
# - Objective (to define an objective function for optimization)
# - SolverFactory (to solve the problem)
# - TransformationFactory (to apply certain transformations)
# - Arc (to connect two unit models)
# - SequentialDecomposition (to initialize the flowsheet in a sequential mode)
#
# For further details on these components, please refer to the pyomo documentation: https://pyomo.readthedocs.io/en/stable/
#


from pyomo.environ import (
    Constraint,
    Var,
    ConcreteModel,
    Expression,
    Objective,
    SolverFactory,
    TerminationCondition,
    TransformationFactory,
    value,
)
from pyomo.network import Arc, SequentialDecomposition

from idaes.core import FlowsheetBlock

from idaes.models.unit_models import (
    PressureChanger,
    Mixer,
    Separator as Splitter,
    Heater,
    StoichiometricReactor,
)

from idaes.models.unit_models import Flash
from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from idaes.core.util.model_statistics import degrees_of_freedom

import idaes.logger as idaeslog
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import InitializationError

from idaes.core.util.structfs.fsrunner import FlowsheetRunner, Context

import hda_ideal_VLE as thermo_props
import hda_reaction as reaction_props


FS = FlowsheetRunner()


@FS.step("build")
def build_model(ctx: Context):
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    add_property_packages(m)
    add_units(m)
    connect_units(m)
    add_expr(m)
    ctx.model = m


# We now need to add the property packages to the flowsheet. Unlike Module 1, where we only had a thermo property package, for this flowsheet we will also need to add a reaction property package.


@FS.substep("build", "add_props")
def add_property_packages(m):
    m.fs.thermo_params = thermo_props.HDAParameterBlock()
    m.fs.reaction_params = reaction_props.HDAReactionParameterBlock(
        property_package=m.fs.thermo_params
    )


@FS.substep("build", "add_units")
def add_units(m):
    """Add the unit models we have imported to the flowsheet.

    Here, we are adding the Mixer (assigned a name M101) and a Heater (assigned a name H101).
    Note that, all unit models need to be given a property package argument.
    """
    m.fs.M101 = Mixer(
        property_package=m.fs.thermo_params,
        inlet_list=["toluene_feed", "hydrogen_feed", "vapor_recycle"],
    )

    m.fs.H101 = Heater(
        property_package=m.fs.thermo_params,
        has_pressure_change=False,
        has_phase_equilibrium=True,
    )

    # Todo: Add reactor with the specifications above
    m.fs.R101 = StoichiometricReactor(
        property_package=m.fs.thermo_params,
        reaction_package=m.fs.reaction_params,
        has_heat_of_reaction=True,
        has_heat_transfer=True,
        has_pressure_change=False,
    )

    m.fs.F101 = Flash(
        property_package=m.fs.thermo_params,
        has_heat_transfer=True,
        has_pressure_change=True,
    )

    m.fs.S101 = Splitter(
        property_package=m.fs.thermo_params,
        ideal_separation=False,
        outlet_list=["purge", "recycle"],
    )

    m.fs.C101 = PressureChanger(
        property_package=m.fs.thermo_params,
        compressor=True,
        thermodynamic_assumption=ThermodynamicAssumption.isothermal,
    )

    m.fs.F102 = Flash(
        property_package=m.fs.thermo_params,
        has_heat_transfer=True,
        has_pressure_change=True,
    )


@FS.substep("build", "create_arcs")
def connect_units(m):
    """Connect Unit Models using Arcs"""
    m.fs.s03 = Arc(source=m.fs.M101.outlet, destination=m.fs.H101.inlet)
    m.fs.s04 = Arc(source=m.fs.H101.outlet, destination=m.fs.R101.inlet)
    m.fs.s05 = Arc(source=m.fs.R101.outlet, destination=m.fs.F101.inlet)
    m.fs.s06 = Arc(source=m.fs.F101.vap_outlet, destination=m.fs.S101.inlet)
    m.fs.s08 = Arc(source=m.fs.S101.recycle, destination=m.fs.C101.inlet)
    m.fs.s09 = Arc(source=m.fs.C101.outlet, destination=m.fs.M101.vapor_recycle)
    m.fs.s10 = Arc(source=m.fs.F101.liq_outlet, destination=m.fs.F102.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)


@FS.substep("build", "add_expressions")
def add_expr(m):
    """Add expressions to compute purity and operating costs"""

    m.fs.purity = Expression(
        expr=m.fs.F102.vap_outlet.flow_mol_phase_comp[0, "Vap", "benzene"]
        / (
            m.fs.F102.vap_outlet.flow_mol_phase_comp[0, "Vap", "benzene"]
            + m.fs.F102.vap_outlet.flow_mol_phase_comp[0, "Vap", "toluene"]
        )
    )
    m.fs.cooling_cost = Expression(
        expr=0.212e-7 * (-m.fs.F101.heat_duty[0]) + 0.212e-7 * (-m.fs.R101.heat_duty[0])
    )
    m.fs.heating_cost = Expression(
        expr=2.2e-7 * m.fs.H101.heat_duty[0] + 1.9e-7 * m.fs.F102.heat_duty[0]
    )
    m.fs.operating_cost = Expression(
        expr=(3600 * 24 * 365 * (m.fs.heating_cost + m.fs.cooling_cost))
    )


@FS.step("set_operating_conditions")
def set_op_cond(ctx):
    m = ctx.model
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(0.30)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "hydrogen"].fix(1e-5)
    m.fs.M101.toluene_feed.flow_mol_phase_comp[0, "Liq", "methane"].fix(1e-5)
    m.fs.M101.toluene_feed.temperature.fix(303.2)
    m.fs.M101.toluene_feed.pressure.fix(350000)

    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "benzene"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "toluene"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "hydrogen"].fix(0.30)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Vap", "methane"].fix(0.02)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "benzene"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "toluene"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "hydrogen"].fix(1e-5)
    m.fs.M101.hydrogen_feed.flow_mol_phase_comp[0, "Liq", "methane"].fix(1e-5)
    m.fs.M101.hydrogen_feed.temperature.fix(303.2)
    m.fs.M101.hydrogen_feed.pressure.fix(350000)

    m.fs.H101.outlet.temperature.fix(600)

    m.fs.R101.conversion = Var(initialize=0.75, bounds=(0, 1))

    m.fs.R101.conv_constraint = Constraint(
        expr=m.fs.R101.conversion
        * m.fs.R101.inlet.flow_mol_phase_comp[0, "Vap", "toluene"]
        == (
            m.fs.R101.inlet.flow_mol_phase_comp[0, "Vap", "toluene"]
            - m.fs.R101.outlet.flow_mol_phase_comp[0, "Vap", "toluene"]
        )
    )

    m.fs.R101.conversion.fix(0.75)
    m.fs.R101.heat_duty.fix(0)

    # Flash 1
    m.fs.F101.vap_outlet.temperature.fix(325.0)
    m.fs.F101.deltaP.fix(0)

    # Flash 2
    m.fs.F102.vap_outlet.temperature.fix(375)
    m.fs.F102.deltaP.fix(-200000)

    # Purge split fraction
    m.fs.S101.split_fraction[0, "purge"].fix(0.2)
    m.fs.C101.outlet.pressure.fix(350000)


@FS.step("initialize")
def initialization(ctx):
    m = ctx.model
    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 3

    # Using the SD tool
    G = seq.create_graph(m)
    heuristic_tear_set = seq.tear_set_arcs(G, method="heuristic")
    order = seq.calculation_order(G)

    tear_guesses = {
        "flow_mol_phase_comp": {
            (0, "Vap", "benzene"): 1e-5,
            (0, "Vap", "toluene"): 1e-5,
            (0, "Vap", "hydrogen"): 0.30,
            (0, "Vap", "methane"): 0.02,
            (0, "Liq", "benzene"): 1e-5,
            (0, "Liq", "toluene"): 0.30,
            (0, "Liq", "hydrogen"): 1e-5,
            (0, "Liq", "methane"): 1e-5,
        },
        "temperature": {0: 303},
        "pressure": {0: 350000},
    }

    # Pass the tear_guess to the SD tool
    seq.set_guesses_for(m.fs.H101.inlet, tear_guesses)

    def init_function(unit):
        try:
            initializer = unit.default_initializer()
            initializer.initialize(unit, output_level=idaeslog.INFO)
        except InitializationError:
            solver = get_solver()
            solver.solve(unit)

    seq.run(m, init_function)


@FS.step("set_solver")
def set_solver(ctx):
    ctx.solver = SolverFactory("ipopt")


@FS.step("solve_initial")
def solve(ctx):
    """Perform the initial model solve."""
    ctx["status"] = results = ctx.solver.solve(ctx.model, tee=ctx["tee"])
    assert results.solver.termination_condition == TerminationCondition.optimal


@FS.step("solve_optimization")
def solve_opt(ctx):
    m = ctx.model
    m.fs.objective = Objective(expr=m.fs.operating_cost)

    m.fs.H101.outlet.temperature.unfix()
    m.fs.R101.heat_duty.unfix()
    m.fs.F101.vap_outlet.temperature.unfix()
    m.fs.F102.vap_outlet.temperature.unfix()

    m.fs.F102.deltaP.unfix()

    assert degrees_of_freedom(m) == 5

    m.fs.H101.outlet.temperature[0].setlb(500)
    m.fs.H101.outlet.temperature[0].setub(600)

    m.fs.R101.outlet.temperature[0].setlb(600)
    m.fs.R101.outlet.temperature[0].setub(800)

    m.fs.F101.vap_outlet.temperature[0].setlb(298.0)
    m.fs.F101.vap_outlet.temperature[0].setub(450.0)
    m.fs.F102.vap_outlet.temperature[0].setlb(298.0)
    m.fs.F102.vap_outlet.temperature[0].setub(450.0)
    m.fs.F102.vap_outlet.pressure[0].setlb(105000)
    m.fs.F102.vap_outlet.pressure[0].setub(110000)

    m.fs.overhead_loss = Constraint(
        expr=m.fs.F101.vap_outlet.flow_mol_phase_comp[0, "Vap", "benzene"]
        <= 0.20 * m.fs.R101.outlet.flow_mol_phase_comp[0, "Vap", "benzene"]
    )

    m.fs.product_flow = Constraint(
        expr=m.fs.F102.vap_outlet.flow_mol_phase_comp[0, "Vap", "benzene"] >= 0.15
    )

    m.fs.product_purity = Constraint(expr=m.fs.purity >= 0.80)

    results = ctx.solver.solve(ctx.model, tee=ctx["tee"])
    ctx["results"] = results
