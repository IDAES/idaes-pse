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
Flowsheet example of the 2-region bubbling fluidized bed (BFB) model
for a methane combustion with iron-oxide case study

Created: 05/04/2020

Author: Chinedu Okoli
"""
# Import python libraries
import time

# Import Pyomo libraries
from pyomo.environ import ConcreteModel, value

# Import IDAES core modules
from idaes.core import FlowsheetBlock
from idaes.core.util import get_solver

# Import IDAES logger
import idaes.logger as idaeslog

# Import BFB unit model
from idaes.gas_solid_contactors.unit_models.bubbling_fluidized_bed \
    import BubblingFluidizedBed

# Import property packages
from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    gas_phase_thermo import GasPhaseParameterBlock
from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    solid_phase_thermo import SolidPhaseParameterBlock
from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    hetero_reactions import HeteroReactionParameterBlock


# -----------------------------------------------------------------------------
def main():

    # ---------------------------------------------------------------------
    # Build model

    # Create a concrete model
    m = ConcreteModel()

    # Create a steady-state flowsheet
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # Set up thermo-physical and reaction properties
    m.fs.gas_properties = GasPhaseParameterBlock()
    m.fs.solid_properties = SolidPhaseParameterBlock()

    m.fs.hetero_reactions = HeteroReactionParameterBlock(
            default={"solid_property_package": m.fs.solid_properties,
                     "gas_property_package": m.fs.gas_properties})

    # Build the BFB in the flowsheet
    m.fs.BFB = BubblingFluidizedBed(
            default={
                    "flow_type": "co_current",
                    "finite_elements": 5,
                    "transformation_method": "dae.collocation",
                    "gas_phase_config":
                    {"property_package": m.fs.gas_properties},
                    "solid_phase_config":
                    {"property_package": m.fs.solid_properties,
                     "reaction_package": m.fs.hetero_reactions
                     }})

    # ---------------------------------------------------------------------
    # Set design and operating variables of the BFB model

    # Fix design variables
    m.fs.BFB.number_orifice.fix(2500)  # [-]
    m.fs.BFB.bed_diameter.fix(6.5)  # m
    m.fs.BFB.bed_height.fix(5)  # m

    # Fix inlet port variables for gas and solid
    m.fs.BFB.gas_inlet.flow_mol[0].fix(272.81)  # mol/s
    m.fs.BFB.gas_inlet.temperature[0].fix(373)  # K
    m.fs.BFB.gas_inlet.pressure[0].fix(1.86)  # bar
    m.fs.BFB.gas_inlet.mole_frac_comp[0, "CO2"].fix(0.4772)
    m.fs.BFB.gas_inlet.mole_frac_comp[0, "H2O"].fix(0.0646)
    m.fs.BFB.gas_inlet.mole_frac_comp[0, "CH4"].fix(0.4582)

    m.fs.BFB.solid_inlet.flow_mass[0].fix(1230)  # kg/s
    # Particle porosity:
    # The porosity of the OC particle at the inlet is calculated from the
    # known bulk density of the fresh OC particle (3251.75 kg/m3), and the
    # skeletal density of the fresh OC particle (calculated from the known
    # composition of the fresh particle, and the skeletal density of its
    # components [see the solids property package])
    m.fs.BFB.solid_inlet.particle_porosity[0].fix(0.27)
    m.fs.BFB.solid_inlet.temperature[0].fix(1186)  # K
    m.fs.BFB.solid_inlet.mass_frac_comp[0, "Fe2O3"].fix(0.45)
    m.fs.BFB.solid_inlet.mass_frac_comp[0, "Fe3O4"].fix(1e-9)
    m.fs.BFB.solid_inlet.mass_frac_comp[0, "Al2O3"].fix(0.55)

    # ---------------------------------------------------------------------
    # Initialize reactor

    t_start = time.time()  # Run start time

    # State arguments for initializing property state blocks
    # Bubble and gas_emulsion temperatures are initialized at solid
    # temperature because thermal mass of solid >> thermal mass of gas
    blk = m.fs.BFB
    gas_phase_state_args = {
            'flow_mol': blk.gas_inlet.flow_mol[0].value,
            'temperature': blk.solid_inlet.temperature[0].value,
            'pressure': blk.gas_inlet.pressure[0].value,
            'mole_frac': {
                'CH4': blk.gas_inlet.mole_frac_comp[0, 'CH4'].value,
                'CO2': blk.gas_inlet.mole_frac_comp[0, 'CO2'].value,
                'H2O': blk.gas_inlet.mole_frac_comp[0, 'H2O'].value}}
    solid_phase_state_args = {
            'flow_mass': blk.solid_inlet.flow_mass[0].value,
            'particle_porosity': blk.solid_inlet.particle_porosity[0].value,
            'temperature': blk.solid_inlet.temperature[0].value,
            'mass_frac': {
                    'Fe2O3': blk.solid_inlet.mass_frac_comp[0, 'Fe2O3'].value,
                    'Fe3O4': blk.solid_inlet.mass_frac_comp[0, 'Fe3O4'].value,
                    'Al2O3': blk.solid_inlet.mass_frac_comp[0, 'Al2O3'].value}}

    m.fs.BFB.initialize(outlvl=idaeslog.INFO,
                        gas_phase_state_args=gas_phase_state_args,
                        solid_phase_state_args=solid_phase_state_args)

    t_initialize = time.time()  # Initialization time

    # ---------------------------------------------------------------------
    # Final solve

    # Create solver
    solver = get_solver()

    solver.solve(m.fs.BFB, tee=True)

    t_simulation = time.time()  # Simulation time

    print("\n")
    print("----------------------------------------------------------")
    print('Total initialization time: ', value(t_initialize - t_start), " s")
    print("----------------------------------------------------------")

    print("\n")
    print("----------------------------------------------------------")
    print('Total simulation time: ', value(t_simulation - t_start), " s")
    print("----------------------------------------------------------")

    return m


# -----------------------------------------------------------------------------
if __name__ == "__main__":
    m = main()

    # Plot some results
    print()
    stream_table = m.fs.BFB._get_stream_table_contents()
    print(stream_table)
