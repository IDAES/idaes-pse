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
Demonstration and test flowsheet for a dynamic flowsheet.

"""

# Import Pyomo libraries
from pyomo.environ import ConcreteModel, SolverFactory, TransformationFactory
from pyomo.network import Arc

# Import IDAES core
from idaes.core import FlowsheetBlock
from idaes.core.util.model_statistics import report_statistics

# Import Unit Model Modules
import idaes.property_models.examples.saponification_thermo as thermo_props
import idaes.property_models.examples.saponification_reactions as reaction_props

# Import Unit Model Modules
from idaes.unit_models import CSTR


def main():
    """
    Make the flowsheet object, fix some variables, and solve the problem
    """
    # Create a Concrete Model as the top level object
    m = ConcreteModel()

    # Add a flowsheet object to the model
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # Add property packages to flowsheet library
    m.fs.thermo_params = thermo_props.SaponificationParameterBlock()
    m.fs.reaction_params = reaction_props.SaponificationReactionParameterBlock(
        default={"property_package": m.fs.thermo_params})

    # Create unit models
    m.fs.Tank1 = CSTR(default={"property_package": m.fs.thermo_params,
                               "reaction_package": m.fs.reaction_params,
                               "has_equilibrium_reactions": False,
                               "has_heat_of_reaction": True,
                               "has_heat_transfer": True,
                               "has_pressure_change": False})
    m.fs.Tank2 = CSTR(default={"property_package": m.fs.thermo_params,
                               "reaction_package": m.fs.reaction_params,
                               "has_equilibrium_reactions": False,
                               "has_heat_of_reaction": True,
                               "has_heat_transfer": True,
                               "has_pressure_change": False})

    # Make Streams to connect units
    m.fs.stream = Arc(source=m.fs.Tank1.outlet,
                      destination=m.fs.Tank2.inlet)

    TransformationFactory("network.expand_arcs").apply_to(m)

    # Set inlet and operating conditions, and some initial conditions.
    m.fs.Tank1.inlet.flow_vol[0].fix(1.0)
    m.fs.Tank1.inlet.conc_mol_comp[0, "H2O"].fix(55388.0)
    m.fs.Tank1.inlet.conc_mol_comp[0, "NaOH"].fix(100.0)
    m.fs.Tank1.inlet.conc_mol_comp[0, "EthylAcetate"].fix(100.0)
    m.fs.Tank1.inlet.conc_mol_comp[0, "SodiumAcetate"].fix(0.0)
    m.fs.Tank1.inlet.conc_mol_comp[0, "Ethanol"].fix(0.0)

    m.fs.Tank1.inlet.temperature.fix(303.15)
    m.fs.Tank1.inlet.pressure.fix(101325.0)

    m.fs.Tank1.volume.fix(1.0)
    m.fs.Tank1.heat_duty.fix(0.0)

    m.fs.Tank2.volume.fix(1.0)
    m.fs.Tank2.heat_duty.fix(0.0)

    # Initialize Units
    m.fs.Tank1.initialize()
    m.fs.Tank2.initialize(state_args={
            "flow_vol": 1.0,
            "conc_mol_comp": {"H2O": 55388.0,
                              "NaOH": 100.0,
                              "EthylAcetate": 100.0,
                              "SodiumAcetate": 0.0,
                              "Ethanol": 0.0},
            "temperature": 303.15,
            "pressure": 101325.0})

    # Create a solver
    solver = SolverFactory('ipopt')
    results = solver.solve(m, tee=False)

    # Print results
    print(results)
    print()
    print("Results")
    print()
    print("Tank 1 Outlet")
    m.fs.Tank1.outlet.display()
    print()
    print("Tank 2 Outlet")
    m.fs.Tank2.outlet.display()

    # For testing purposes
    return(m, results)


if __name__ == "__main__":
    main()
