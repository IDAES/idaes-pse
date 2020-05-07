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
from pyomo.environ import ConcreteModel, SolverFactory, value
from pyomo.util.calc_var_value import calculate_variable_from_constraint

# Import IDAES core modules
from idaes.core import FlowsheetBlock, EnergyBalanceType

# Import BFB unit model
from idaes.gas_solid_contactors.unit_models.bubbling_fluidized_bed \
    import BubblingFluidizedBed

# Import property packages
from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    gas_phase_thermo import Gas_Phase_Thermo_ParameterBlock
from idaes.gas_solid_contactors.properties.methane_iron_OC_reduction. \
    solid_phase_thermo import Solid_Phase_Thermo_ParameterBlock
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
    m.fs.gas_properties = Gas_Phase_Thermo_ParameterBlock()
    m.fs.solid_properties = Solid_Phase_Thermo_ParameterBlock()

    m.fs.hetero_reactions = HeteroReactionParameterBlock(
            default={"solid_property_package": m.fs.solid_properties,
                     "gas_property_package": m.fs.gas_properties})

    # Build the BFB in the flowsheet
    m.fs.BFB = BubblingFluidizedBed(
            default={
                    "flow_type": "co_current",
                    "finite_elements": 5,
                    "transformation_method": "dae.collocation",
                    "bubble_config":
                    {"property_package": m.fs.gas_properties},
                    "gas_emulsion_config":
                    {"property_package": m.fs.gas_properties,
                     "has_pressure_change": True},
                    "solid_emulsion_config":
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
    m.fs.BFB.gas_inlet.mole_frac[0, "CO2"].fix(0.4772)
    m.fs.BFB.gas_inlet.mole_frac[0, "H2O"].fix(0.0646)
    m.fs.BFB.gas_inlet.mole_frac[0, "CH4"].fix(0.4582)

    m.fs.BFB.solid_inlet.flow_mass[0].fix(1230)  # kg/s
    m.fs.BFB.solid_inlet.temperature[0].fix(1186)  # K
    m.fs.BFB.solid_inlet.mass_frac[0, "Fe2O3"].fix(0.45)
    m.fs.BFB.solid_inlet.mass_frac[0, "Fe3O4"].fix(1e-9)
    m.fs.BFB.solid_inlet.mass_frac[0, "Al2O3"].fix(0.55)

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
                'CH4': blk.gas_inlet.mole_frac[0, 'CH4'].value,
                'CO2': blk.gas_inlet.mole_frac[0, 'CO2'].value,
                'H2O': blk.gas_inlet.mole_frac[0, 'H2O'].value}}
    solid_phase_state_args = {
            'flow_mass': blk.solid_inlet.flow_mass[0].value,
            'temperature': blk.solid_inlet.temperature[0].value,
            'mass_frac': {
                    'Fe2O3': blk.solid_inlet.mass_frac[0, 'Fe2O3'].value,
                    'Fe3O4': blk.solid_inlet.mass_frac[0, 'Fe3O4'].value,
                    'Al2O3': blk.solid_inlet.mass_frac[0, 'Al2O3'].value}}

    m.fs.BFB.initialize(outlvl=0,
                        gas_phase_state_args=gas_phase_state_args,
                        solid_phase_state_args=solid_phase_state_args)

    t_initialize = time.time()  # Initialization time

    # ---------------------------------------------------------------------
    # Final solve

    # Create a solver
    solver = SolverFactory('ipopt')
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

# %% # ------------------------------------------------------------------------


# Function to print some variables of interest
def print_summary(self):
    if m.fs.BFB.config.solid_emulsion_config.reaction_package is not None:
        print()
        print("Reaction stoichiometry ratio check (expected value is 12) =>")

        mole_gas_reacted = (
                m.fs.BFB.gas_inlet.flow_mol[0].value *
                m.fs.BFB.gas_inlet.mole_frac[0, 'CH4'].value -
                m.fs.BFB.gas_outlet.flow_mol[0].value *
                m.fs.BFB.gas_outlet.mole_frac[0, 'CH4'].value)

        mole_solid_reacted = (
                (m.fs.BFB.solid_inlet.flow_mass[0].value *
                 m.fs.BFB.solid_inlet.mass_frac[0, 'Fe2O3'].value /
                 m.fs.BFB.solid_inlet_block[0]._params.mw['Fe2O3']) -
                (m.fs.BFB.solid_outlet.flow_mass[0].value *
                 m.fs.BFB.solid_outlet.mass_frac[0, 'Fe2O3'].value /
                 m.fs.BFB.solid_outlet_block[0]._params.mw['Fe2O3']))

        print('stoichiometric ratio:', mole_solid_reacted/mole_gas_reacted)

    print()
    print('Mass balance closure check =>')

    calculate_variable_from_constraint(
                m.fs.BFB.gas_inlet_block[0].mw_gas,
                m.fs.BFB.gas_inlet_block[0].mw_gas_eqn)
    calculate_variable_from_constraint(
                m.fs.BFB.gas_outlet_block[0].mw_gas,
                m.fs.BFB.gas_outlet_block[0].mw_gas_eqn)

    mbal_gas = ((m.fs.BFB.gas_inlet.flow_mol[0].value *
                m.fs.BFB.gas_inlet_block[0].mw_gas.value) -
                (m.fs.BFB.gas_outlet_block[0].flow_mol.value *
                m.fs.BFB.gas_outlet_block[0].mw_gas.value))

    mbal_solid = (
            m.fs.BFB.solid_inlet.flow_mass[0].value -
            m.fs.BFB.solid_outlet.flow_mass[0].value)

    mbal_tol = mbal_gas + mbal_solid

    print('Mass Balance Tolerance:', mbal_tol)
    print('Mass balance gas_FR:', mbal_gas)
    print('Mass balance solids_FR:', mbal_solid)

    if m.fs.BFB.config.energy_balance_type != EnergyBalanceType.none:
        print()
        print('Energy balance closure check =>')
        ebal_gas = (
            (m.fs.BFB.gas_inlet.flow_mol[0].value *
             m.fs.BFB.gas_inlet_block[0].enth_mol.value) -
            (m.fs.BFB.gas_outlet.flow_mol[0].value *
             m.fs.BFB.gas_outlet_block[0].enth_mol.value))

        ebal_solid = (
            (m.fs.BFB.solid_inlet.flow_mass[0].value *
             m.fs.BFB.solid_inlet_block[0].enth_mass.value) -
            (m.fs.BFB.solid_outlet.flow_mass[0].value *
             m.fs.BFB.solid_outlet_block[0].enth_mass.value))

        if m.fs.BFB.config.solid_emulsion_config.reaction_package is not None:
            e_reaction = (mole_gas_reacted *
                          m.fs.BFB.solid_emulsion_region.reactions[0, 0].
                          _params.dh_rxn["R1"])
        else:
            e_reaction = 0

        ebal_tol = ebal_gas + ebal_solid - e_reaction

        print('Energy Balance Tolerance_FR:', ebal_tol)
        print('Energy balance gas_FR:', ebal_gas)
        print('Energy balance solids_FR:', ebal_solid)
        print('Total heat due to reaction:', e_reaction)

    if m.fs.BFB.config.solid_emulsion_config.reaction_package is not None:
        print()
        # Oxygen carrier and fuel conversion
        Conv_gas = 1 - (
            ((m.fs.BFB.gas_outlet.flow_mol[0].value *
              m.fs.BFB.gas_outlet_block[0].mole_frac['CH4'].value)) /
            (m.fs.BFB.gas_inlet.flow_mol[0].value *
             m.fs.BFB.gas_inlet_block[0].mole_frac['CH4'].value))

        if m.fs.BFB.config.flow_type == "counter_current":
            Conv_OC = (m.fs.BFB.solid_emulsion_region.reactions[0, 0].
                       OC_conv.value)
        else:
            Conv_OC = (m.fs.BFB.solid_emulsion_region.reactions[0, 1].
                       OC_conv.value)

        print('Fuel and OC conversion =>')
        print("CH4 conversion:", value(Conv_gas)*100, "%")
        print("Fe2O3 conversion:", value(Conv_OC)*100, "%")


# %%        # -----------------------------------------------------------------
if __name__ == "__main__":
    m = main()

    # Plot some results
    print_summary(m)
    m.fs.BFB._results_plot_FR()
    stream_table = m.fs.BFB._get_stream_table_contents()
    print(stream_table)
