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
"""
Flowsheet example of the 0D FixedBed model for an iron-oxide reduction
with methane case study

Created: 05/14/2020

Author: Chinedu Okoli
"""

import time

from pyomo.environ import ConcreteModel, TransformationFactory, value, units as pyunits

from idaes.core import FlowsheetBlock, EnergyBalanceType
from idaes.core.util.initialization import initialize_by_time_element
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver

from idaes.models_extra.gas_solid_contactors.unit_models.fixed_bed_0D import FixedBed0D
from idaes.models_extra.gas_solid_contactors.properties.methane_iron_OC_reduction.gas_phase_thermo import (
    GasPhaseParameterBlock,
)
from idaes.models_extra.gas_solid_contactors.properties.methane_iron_OC_reduction.solid_phase_thermo import (
    SolidPhaseParameterBlock,
)
from idaes.models_extra.gas_solid_contactors.properties.methane_iron_OC_reduction.hetero_reactions import (
    HeteroReactionParameterBlock,
)


# -----------------------------------------------------------------------------
def main(m):
    m.fs = FlowsheetBlock(dynamic=True, time_set=[0, 3600], time_units=pyunits.s)

    m.fs.gas_props = GasPhaseParameterBlock()
    m.fs.solid_props = SolidPhaseParameterBlock()
    m.fs.solid_rxns = HeteroReactionParameterBlock(
        solid_property_package=m.fs.solid_props, gas_property_package=m.fs.gas_props
    )

    m.fs.TGA = FixedBed0D(
        energy_balance_type=EnergyBalanceType.none,
        gas_property_package=m.fs.gas_props,
        solid_property_package=m.fs.solid_props,
        reaction_package=m.fs.solid_rxns,
    )

    # Discretize time domain
    m.discretizer = TransformationFactory("dae.finite_difference")
    m.discretizer.apply_to(m, nfe=100, wrt=m.fs.time, scheme="BACKWARD")

    # Set reactor design conditions
    m.fs.TGA.bed_diameter.fix(1)  # diameter of the TGA reactor [m]
    m.fs.TGA.bed_height.fix(1)  # height of solids in the TGA reactor [m]

    # Set initial conditions of the solid phase
    m.fs.TGA.solids[0].particle_porosity.fix(0.20)
    m.fs.TGA.solids[0].mass_frac_comp["Fe2O3"].fix(0.45)
    m.fs.TGA.solids[0].mass_frac_comp["Fe3O4"].fix(0)
    m.fs.TGA.solids[0].mass_frac_comp["Al2O3"].fix(0.55)
    m.fs.TGA.solids[0].temperature.fix(1273.15)

    # Set conditions of the gas phase (this is all fixed as gas side assumption
    # is excess gas flowrate which means all state variables remain unchanged)
    for t in m.fs.time:
        m.fs.TGA.gas[t].temperature.fix(1273.15)
        m.fs.TGA.gas[t].pressure.fix(1.01325e5)  # 1atm
        m.fs.TGA.gas[t].mole_frac_comp["CO2"].fix(0.4)
        m.fs.TGA.gas[t].mole_frac_comp["H2O"].fix(0.5)
        m.fs.TGA.gas[t].mole_frac_comp["CH4"].fix(0.1)

    # Solver options
    optarg = {"tol": 1e-6}

    t_start = time.time()  # Run start time

    print()
    print("Apply scaling transformation")
    # Scale the model by applying scaling transformation
    # This reduces ill conditioning of the model
    iscale.calculate_scaling_factors(m)

    print()
    print("Initialize the model")
    m.fs.TGA.initialize(optarg=optarg)

    t_initialize = time.time()  # Initialization time

    solver = get_solver("ipopt", optarg)  # create solver

    initialize_by_time_element(m.fs, m.fs.time, solver=solver)
    solver.solve(m, tee=True)

    t_simulation = time.time()  # Simulation time

    print("\n")
    print("----------------------------------------------------------")
    print("Total initialization time: ", value(t_initialize - t_start), " s")
    print("----------------------------------------------------------")

    print("\n")
    print("----------------------------------------------------------")
    print("Total simulation time: ", value(t_simulation - t_start), " s")
    print("----------------------------------------------------------")

    return m


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    m = ConcreteModel()
    m = main(m)
