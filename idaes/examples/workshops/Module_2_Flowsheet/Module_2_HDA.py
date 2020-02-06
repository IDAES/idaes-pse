##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Demonstration of HDA flowsheet with optimization.

"""

# Import Pyomo libraries
from pyomo.environ import (Constraint,
                           ConcreteModel,
                           Expression,
                           Objective,
                           minimize,
                           SolverFactory,
                           TransformationFactory,
                           value,
                           Var)
from pyomo.network import Arc, SequentialDecomposition

# Import IDAES core
from idaes.core import FlowsheetBlock

# Import Unit Model Modules
from . import hda_ideal_VLE as thermo_props
from . import hda_reaction as reaction_props

# Import Unit Model Modules
from idaes.unit_models import (PressureChanger,
                               StoichiometricReactor,
                               Flash,
                               Heater,
                               Mixer,
                               Separator as Splitter)

from idaes.unit_models.pressure_changer import ThermodynamicAssumption

from idaes.core.util.model_statistics import degrees_of_freedom


def main():
    # Create a Concrete Model as the top level object
    m = ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # Add property packages to flowsheet library
    m.fs.thermo_params = thermo_props.HDAParameterBlock()
    m.fs.reaction_params = reaction_props.HDAReactionParameterBlock(
        default={"property_package": m.fs.thermo_params})

    # Create unit models
    m.fs.M101 = Mixer(default={"property_package": m.fs.thermo_params,
                               "inlet_list": ["toluene_feed", "hydrogen_feed",
                                              "vapor_recycle"]})

    m.fs.H101 = Heater(default={"property_package": m.fs.thermo_params,
                                "has_pressure_change": False,
                                "has_phase_equilibrium": True})

    m.fs.R101 = StoichiometricReactor(
            default={"property_package": m.fs.thermo_params,
                     "reaction_package": m.fs.reaction_params,
                     "has_heat_of_reaction": True,
                     "has_heat_transfer": True,
                     "has_pressure_change": False})

    m.fs.F101 = Flash(default={"property_package": m.fs.thermo_params,
                               "has_heat_transfer": True,
                               "has_pressure_change": True})

    m.fs.S101 = Splitter(default={"property_package": m.fs.thermo_params,
                                  "ideal_separation": False,
                                  "outlet_list": ["purge", "recycle"]})

    # This is needed to avoid pressure degeneracy in recylce loop
    m.fs.C101 = PressureChanger(default={
            "property_package": m.fs.thermo_params,
            "compressor": True,
            "thermodynamic_assumption": ThermodynamicAssumption.isothermal})

    m.fs.F102 = Flash(default={"property_package": m.fs.thermo_params,
                               "has_heat_transfer": True,
                               "has_pressure_change": True})

    m.fs.H101.control_volume.scaling_factor_energy = 1e-3
    m.fs.R101.control_volume.scaling_factor_energy = 1e-3
    m.fs.F101.control_volume.scaling_factor_energy = 1e-3
    m.fs.C101.control_volume.scaling_factor_energy = 1e-3
    m.fs.F102.control_volume.scaling_factor_energy = 1e-3

    # Connect units
    m.fs.s03 = Arc(source=m.fs.M101.outlet, destination=m.fs.H101.inlet)
    m.fs.s04 = Arc(source=m.fs.H101.outlet, destination=m.fs.R101.inlet)
    m.fs.s05 = Arc(source=m.fs.R101.outlet, destination=m.fs.F101.inlet)
    m.fs.s06 = Arc(source=m.fs.F101.vap_outlet, destination=m.fs.S101.inlet)
    m.fs.s08 = Arc(source=m.fs.S101.recycle, destination=m.fs.C101.inlet)
    m.fs.s09 = Arc(source=m.fs.C101.outlet,
                   destination=m.fs.M101.vapor_recycle)
    m.fs.s10 = Arc(source=m.fs.F101.liq_outlet, destination=m.fs.F102.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    # Set operating conditions
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
        expr=m.fs.R101.conversion * m.fs.R101.inlet.
        flow_mol_phase_comp[0, "Vap", "toluene"] ==
        (m.fs.R101.inlet.flow_mol_phase_comp[0, "Vap", "toluene"] -
         m.fs.R101.outlet.flow_mol_phase_comp[0, "Vap", "toluene"]))

    m.fs.R101.conversion.fix(0.75)
    m.fs.R101.heat_duty.fix(0)

    m.fs.F101.vap_outlet.temperature.fix(325.0)
    m.fs.F101.deltaP.fix(0)

    m.fs.S101.split_fraction[0, "purge"].fix(0.2)

    m.fs.C101.outlet.pressure.fix(350000)

    m.fs.F102.vap_outlet.temperature.fix(375)
    m.fs.F102.deltaP.fix(-200000)

    # Define expressions
    # Product purity
    m.fs.purity = Expression(
        expr=m.fs.F102.vap_outlet.flow_mol_phase_comp[0, "Vap", "benzene"] /
        (m.fs.F102.vap_outlet.flow_mol_phase_comp[0, "Vap", "benzene"]
         + m.fs.F102.vap_outlet.flow_mol_phase_comp[0, "Vap", "toluene"]))

    # Operating cost ($/yr)
    m.fs.cooling_cost = Expression(expr=0.212e-7 * -m.fs.F101.heat_duty[0] +
                                   0.212e-7 * -m.fs.R101.heat_duty[0])
    m.fs.heating_cost = Expression(expr=2.2e-7 * m.fs.H101.heat_duty[0] +
                                   1.9e-7 * m.fs.F102.heat_duty[0])
    m.fs.operating_cost = Expression(expr=(3600 * 24 * 365 *
                                           (m.fs.heating_cost +
                                            m.fs.cooling_cost)))
    print(degrees_of_freedom(m))

    # Initialize Units
    # Define method for initialising each block
    def function(unit):
        unit.initialize(outlvl=1)

    # Create instance of sequential decomposition tool
    seq = SequentialDecomposition()
    seq.options.select_tear_method = "heuristic"
    seq.options.tear_method = "Wegstein"
    seq.options.iterLim = 5

    # Determine tear stream and calculation order
    G = seq.create_graph(m)
    heu_result = seq.tear_set_arcs(G, method="heuristic")
    order = seq.calculation_order(G)

    # Display tear stream and calculation order
    for o in heu_result:
        print(o.name)
    for o in order:
        for oo in o:
            print(oo.name)

    # Set guesses for tear stream
    tear_guesses = {
        "flow_mol_phase_comp": {
                (0, "Vap", "benzene"): 1e-5,
                (0, "Vap", "toluene"): 1e-5,
                (0, "Vap", "hydrogen"): 0.30,
                (0, "Vap", "methane"): 0.02,
                (0, "Liq", "benzene"): 1e-5,
                (0, "Liq", "toluene"): 0.30,
                (0, "Liq", "hydrogen"): 1e-5,
                (0, "Liq", "methane"): 1e-5},
        "temperature": {0: 303},
        "pressure": {0: 350000}}
    seq.set_guesses_for(m.fs.H101.inlet, tear_guesses)

    # Run sequential initialization
    seq.run(m, function)

    # # Create a solver
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6}
    solver.options = {'tol': 1e-6, 'max_iter': 5000}
    results = solver.solve(m, tee=True)

    # Print results
    print("M101 Outlet")
    m.fs.M101.outlet.display()

    print("H101 Outlet")
    m.fs.H101.outlet.display()

    print("R101 Outlet")
    m.fs.R101.outlet.display()

    print("F101")
    m.fs.F101.liq_outlet.display()
    m.fs.F101.vap_outlet.display()

    print("F102")
    m.fs.F102.liq_outlet.display()
    m.fs.F102.vap_outlet.display()

    print("Purge")
    m.fs.S101.purge.display()

    print("Purity:", value(m.fs.purity))

    # Optimize process
    m.fs.objective = Objective(sense=minimize, expr=m.fs.operating_cost)

    # Decision variables
    m.fs.H101.outlet.temperature.unfix()
    m.fs.R101.heat_duty.unfix()
    m.fs.F101.vap_outlet.temperature.unfix()
    m.fs.F102.vap_outlet.temperature.unfix()
    m.fs.F102.deltaP.unfix()

    # Variable bounds
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

    # Additional Constraints
    m.fs.overhead_loss = Constraint(
        expr=m.fs.F101.vap_outlet.flow_mol_phase_comp[0, "Vap", "benzene"] <=
        0.20 * m.fs.R101.outlet.flow_mol_phase_comp[0, "Vap", "benzene"])

    m.fs.product_flow = Constraint(
        expr=m.fs.F102.vap_outlet.flow_mol_phase_comp[0, "Vap", "benzene"] >=
        0.15)

    m.fs.product_purity = Constraint(expr=m.fs.purity >= 0.80)

    # Create a solver
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6}
    solver.options = {'tol': 1e-6, 'max_iter': 5000}
    results = solver.solve(m, tee=True)

    # Print optimization results
    print()
    print("Optimal Solution")
    m.fs.operating_cost.display()
    m.fs.H101.heat_duty.display()
    m.fs.R101.heat_duty.display()
    m.fs.F101.heat_duty.display()
    m.fs.F102.heat_duty.display()

    # Print results
    print("M101 Outlet")
    m.fs.M101.outlet.display()

    print("H101 Outlet")
    m.fs.H101.outlet.display()

    print("R101 Outlet")
    m.fs.R101.outlet.display()

    print("F101")
    m.fs.F101.liq_outlet.display()
    m.fs.F101.vap_outlet.display()

    print("F102")
    m.fs.F102.liq_outlet.display()
    m.fs.F102.vap_outlet.display()

    print("Purge")
    m.fs.S101.purge.display()

    print("Recycle")
    m.fs.S101.recycle.display()

    print("Purity:", value(m.fs.purity))

    # For testing purposes
    return(m, results)


if __name__ == "__main__":
    main()
