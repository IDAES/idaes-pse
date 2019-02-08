##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Demonstration and test flowsheet for a dynamic flowsheet.

"""
from __future__ import division
from __future__ import print_function

# Import Pyomo libraries
from pyomo.environ import ConcreteModel, SolverFactory, TransformationFactory
from pyomo.network import Arc

# Import IDAES core
from idaes.core import FlowsheetBlock

# Import Unit Model Modules
import idaes.property_models.saponification_thermo as thermo_props
import idaes.property_models.saponification_reactions as reaction_props

# Import Unit Model Modules
from idaes.unit_models import CSTR, Mixer


def main():
    """
    Make the flowsheet object, fix some variables, and solve the problem
    """
    # Create a Concrete Model as the top level object
    m = ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(default={"dynamic": True,
                                   "time_set": [0, 1, 1000]})

    # Add property packages to flowsheet library
    m.fs.thermo_params = thermo_props.PhysicalParameterBlock()
    m.fs.reaction_params = reaction_props.ReactionParameterBlock(default={
                            "property_package": m.fs.thermo_params})

    # Create unit models
    m.fs.mix = Mixer(default={"dynamic": False,
                              "property_package": m.fs.thermo_params})
    m.fs.tank1 = CSTR(default={"property_package": m.fs.thermo_params,
                               "reaction_package": m.fs.reaction_params,
                               "has_holdup": True,
                               "has_equilibrium_reactions": False,
                               "has_heat_transfer": True,
                               "has_pressure_change": False})
    m.fs.tank2 = CSTR(default={"property_package": m.fs.thermo_params,
                               "reaction_package": m.fs.reaction_params,
                               "has_holdup": True,
                               "has_equilibrium_reactions": False,
                               "has_heat_transfer": True,
                               "has_pressure_change": False})

    # Make Streams to connect units
    m.fs.stream1 = Arc(source=m.fs.mix.outlet,
                       destination=m.fs.tank1.inlet)

    m.fs.stream2 = Arc(source=m.fs.tank1.outlet,
                       destination=m.fs.tank2.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    # Discretize time domain
    m.discretizer = TransformationFactory('dae.finite_difference')
    m.discretizer.apply_to(m,
                           nfe=4,
                           wrt=m.fs.time,
                           scheme="BACKWARD")

    # Set inlet and operating conditions, and some initial conditions.
    m.fs.mix.inlet_1[0].flow_vol.fix(0.5)
    m.fs.mix.inlet_1[0].conc_mol_comp["H2O"].fix(55388.0)
    m.fs.mix.inlet_1[0].conc_mol_comp["NaOH"].fix(100.0)
    m.fs.mix.inlet_1[0].conc_mol_comp["EthylAcetate"].fix(0.0)
    m.fs.mix.inlet_1[0].conc_mol_comp["SodiumAcetate"].fix(0.0)
    m.fs.mix.inlet_1[0].conc_mol_comp["Ethanol"].fix(0.0)
    m.fs.mix.inlet_1[0].temperature.fix(303.15)
    m.fs.mix.inlet_1[0].pressure.fix(101325.0)

    m.fs.mix.inlet_2[0].flow_vol.fix(0.5)
    m.fs.mix.inlet_2[0].conc_mol_comp["H2O"].fix(55388.0)
    m.fs.mix.inlet_2[0].conc_mol_comp["NaOH"].fix(0.0)
    m.fs.mix.inlet_2[0].conc_mol_comp["EthylAcetate"].fix(100.0)
    m.fs.mix.inlet_2[0].conc_mol_comp["SodiumAcetate"].fix(0.0)
    m.fs.mix.inlet_2[0].conc_mol_comp["Ethanol"].fix(0.0)
    m.fs.mix.inlet_2[0].temperature.fix(303.15)
    m.fs.mix.inlet_2[0].pressure.fix(101325.0)

    m.fs.tank1.inlet[0].temperature.fix(303.15)
    m.fs.tank1.inlet[0].pressure.fix(101325.0)

    m.fs.tank1.volume.fix(1.0)
    m.fs.tank1.heat_duty.fix(0.0)

    m.fs.tank2.volume.fix(1.0)
    m.fs.tank2.heat_duty.fix(0.0)

    # Initialize Units
    m.fs.mix.initialize()
#
#    m.fs.tank1.initialize(state_args={
#            "flow_vol": 1.0,
#            "conc_mol_comp": {"H2O": 55388.0,
#                              "NaOH": 100.0,
#                              "EthylAcetate": 100.0,
#                              "SodiumAcetate": 0.0,
#                              "Ethanol": 0.0},
#            "temperature": 303.15,
#            "pressure": 101325.0})
#
#    m.fs.tank2.initialize(state_args={
#            "flow_vol": 1.0,
#            "conc_mol_comp": {"H2O": 55388.0,
#                              "NaOH": 100.0,
#                              "EthylAcetate": 100.0,
#                              "SodiumAcetate": 0.0,
#                              "Ethanol": 0.0},
#            "temperature": 303.15,
#            "pressure": 101325.0})

    # Create a solver
    solver = SolverFactory('ipopt')
    results = solver.solve(m.fs.mix, tee=True)

    # Print results
#    print(results)
#    print()
#    print("Results")
#    print()
#    print("Tank 1 Outlet")
#    m.fs.tank1.outlet.display()
#    print()
#    print("Tank 2 Outlet")
#    m.fs.tank2.outlet.display()
    m.fs.mix.outlet.display()

    # For testing purposes
    return(m, results)


if __name__ == "__main__":
    main()
