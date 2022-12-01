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
Flowsheet example of the MB model for a methane combustion with
iron-oxide case study


Author: Chinedu Okoli
"""

import time

from pyomo.environ import ConcreteModel, value

from idaes.core import FlowsheetBlock
from idaes.core.util import scaling as iscale
from idaes.core.solvers import get_solver

# Import IDAES logger
import idaes.logger as idaeslog

# Import MBR unit model
from idaes.models_extra.gas_solid_contactors.unit_models.moving_bed import MBR

# Import property packages
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
def main():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # Set up thermo props and reaction props
    m.fs.gas_properties = GasPhaseParameterBlock()
    m.fs.solid_properties = SolidPhaseParameterBlock()

    m.fs.hetero_reactions = HeteroReactionParameterBlock(
        solid_property_package=m.fs.solid_properties,
        gas_property_package=m.fs.gas_properties,
    )

    m.fs.MB = MBR(
        transformation_method="dae.collocation",
        gas_phase_config={"property_package": m.fs.gas_properties},
        solid_phase_config={
            "property_package": m.fs.solid_properties,
            "reaction_package": m.fs.hetero_reactions,
        },
    )

    # Fix bed geometry variables
    m.fs.MB.bed_diameter.fix(6.5)  # m
    m.fs.MB.bed_height.fix(5)  # m

    # Fix inlet port variables for gas and solid
    m.fs.MB.gas_inlet.flow_mol[0].fix(128.20513)  # mol/s
    m.fs.MB.gas_inlet.temperature[0].fix(298.15)  # K
    m.fs.MB.gas_inlet.pressure[0].fix(2.00e5)  # Pa = 1E5 bar
    m.fs.MB.gas_inlet.mole_frac_comp[0, "CO2"].fix(0.02499)
    m.fs.MB.gas_inlet.mole_frac_comp[0, "H2O"].fix(0.00001)
    m.fs.MB.gas_inlet.mole_frac_comp[0, "CH4"].fix(0.975)

    m.fs.MB.solid_inlet.flow_mass[0].fix(591.4)  # kg/s
    # Particle porosity:
    # The porosity of the OC particle at the inlet is calculated from the
    # known bulk density of the fresh OC particle (3251.75 kg/m3), and the
    # skeletal density of the fresh OC particle (calculated from the known
    # composition of the fresh particle, and the skeletal density of its
    # components [see the solids property package])
    m.fs.MB.solid_inlet.particle_porosity[0].fix(0.27)
    m.fs.MB.solid_inlet.temperature[0].fix(1183.15)  # K
    m.fs.MB.solid_inlet.mass_frac_comp[0, "Fe2O3"].fix(0.45)
    m.fs.MB.solid_inlet.mass_frac_comp[0, "Fe3O4"].fix(1e-9)
    m.fs.MB.solid_inlet.mass_frac_comp[0, "Al2O3"].fix(0.55)

    # Initialize fuel reactor
    t_start = time.time()  # Run start time

    # State arguments for initializing property state blocks
    # Gas phase temperature is initialized at solid
    # temperature because thermal mass of solid >> thermal mass of gas
    # Particularly useful for initialization if reaction takes place
    blk = m.fs.MB
    gas_phase_state_args = {
        "flow_mol": blk.gas_inlet.flow_mol[0].value,
        "temperature": blk.solid_inlet.temperature[0].value,
        "pressure": blk.gas_inlet.pressure[0].value,
        "mole_frac": {
            "CH4": blk.gas_inlet.mole_frac_comp[0, "CH4"].value,
            "CO2": blk.gas_inlet.mole_frac_comp[0, "CO2"].value,
            "H2O": blk.gas_inlet.mole_frac_comp[0, "H2O"].value,
        },
    }
    solid_phase_state_args = {
        "flow_mass": blk.solid_inlet.flow_mass[0].value,
        "particle_porosity": blk.solid_inlet.particle_porosity[0].value,
        "temperature": blk.solid_inlet.temperature[0].value,
        "mass_frac": {
            "Fe2O3": blk.solid_inlet.mass_frac_comp[0, "Fe2O3"].value,
            "Fe3O4": blk.solid_inlet.mass_frac_comp[0, "Fe3O4"].value,
            "Al2O3": blk.solid_inlet.mass_frac_comp[0, "Al2O3"].value,
        },
    }

    print()
    print("Apply scaling transformation")
    # Scale the model by applying scaling transformation
    # This reduces ill conditioning of the model
    iscale.calculate_scaling_factors(m)

    print()
    print("Initialize the model")
    m.fs.MB.initialize(
        outlvl=idaeslog.INFO,
        optarg={"tol": 1e-5},
        gas_phase_state_args=gas_phase_state_args,
        solid_phase_state_args=solid_phase_state_args,
    )

    t_initialize = time.time()  # Initialization time

    print()
    print("Solve the model")
    # Create a solver
    solver = get_solver()
    solver.solve(m.fs.MB, tee=True, options={"tol": 1e-5})

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
    m = main()
    stream_table = m.fs.MB._get_stream_table_contents()
    print(stream_table)
