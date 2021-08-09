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
Flowsheet example of the MB model for a syngas combustion with
iron-oxide case study

Author: Chinedu Okoli, Anca Ostace
"""
import time
import sys
import os

from pyomo.environ import ConcreteModel, SolverFactory, value

from idaes.core import FlowsheetBlock

# Access parent directory (chemical_looping) of the current directory
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

from unit_models.moving_bed import MovingBed
from property_packages.syngas_iron_OC_reduction.gas_phase_thermo \
    import Gas_Phase_Thermo_ParameterBlock
from property_packages.syngas_iron_OC_reduction.solid_phase_thermo \
    import Solid_Phase_Thermo_ParameterBlock
from property_packages.syngas_iron_OC_reduction.hetero_reactions \
    import HeteroReactionParameterBlock


# -----------------------------------------------------------------------------
def main():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})

    # Set up thermo props and reaction props
    m.fs.gas_properties = Gas_Phase_Thermo_ParameterBlock()
    m.fs.solid_properties = Solid_Phase_Thermo_ParameterBlock()
    m.fs.hetero_reactions = HeteroReactionParameterBlock(
        default={"solid_property_package": m.fs.solid_properties,
                 "gas_property_package": m.fs.gas_properties})

    # Create a custom grid, fe_set
    nfe = 10

    fe_set = [0]
#    fe_a = 1/4.0
#    fe_b = 0.2
#    for i in range(1, nfe+1):
#        if i < nfe*fe_a:
#            fe_set.append(i*fe_b/(nfe*fe_a))
#        elif i == nfe:
#            fe_set.append(1)
#        else:
#            fe_set.append(fe_b + (i-nfe*fe_a)*(1-fe_b)/(nfe*(1-fe_a)))

    m.fs.MB = MovingBed(default={
        "finite_elements": nfe,
        "length_domain_set": fe_set,
        "transformation_method": "dae.collocation",
        "transformation_scheme": "LAGRANGE-RADAU",
        "gas_phase_config":
            {"property_package": m.fs.gas_properties,
             "has_pressure_change": True,
             "pressure_drop_type": "simple_correlation"},
        "solid_phase_config": {"property_package": m.fs.solid_properties,
                               "reaction_package": m.fs.hetero_reactions}})

    # Fix bed geometry variables
    m.fs.MB.bed_diameter.fix(6.5)  # m
    m.fs.MB.bed_height.fix(5)  # m

    # Fix inlet port variables for gas and solid
    m.fs.MB.gas_inlet.flow_mol[0].fix(130)  # 128.20513)  # mol/s
    m.fs.MB.gas_inlet.temperature[0].fix(298.15)  # K
    m.fs.MB.gas_inlet.pressure[0].fix(2.00)  # bar

    m.fs.MB.gas_inlet.mole_frac[0, "CO2"].fix(0.1622)  # 0.02499)
    m.fs.MB.gas_inlet.mole_frac[0, "H2O"].fix(0.0186)  # 0.00001)
    m.fs.MB.gas_inlet.mole_frac[0, "CH4"].fix(0.0217)  # 0.325)
    m.fs.MB.gas_inlet.mole_frac[0, "CO"].fix(0.4631)  # 0.325)
    m.fs.MB.gas_inlet.mole_frac[0, "H2"].fix(0.3344)  # 0.325)

    m.fs.MB.solid_inlet.flow_mass[0].fix(122.4)  # 591.4) # kg/s
    m.fs.MB.solid_inlet.temperature[0].fix(1183.15)  # K
    m.fs.MB.solid_inlet.mass_frac[0, "Fe2O3"].fix(0.45)
    m.fs.MB.solid_inlet.mass_frac[0, "Fe3O4"].fix(1e-9)
    m.fs.MB.solid_inlet.mass_frac[0, "Al2O3"].fix(0.55)

    # Initialize fuel reactor
    t_start = time.time()  # Run start time

    # State arguments for initializing property state blocks
    # Gas phase temperature is initialized at solid
    # temperature because thermal mass of solid >> thermal mass of gas
    blk = m.fs.MB
    gas_phase_state_args = {
        'flow_mol': blk.gas_inlet.flow_mol[0].value,
        'temperature': blk.solid_inlet.temperature[0].value,
        'pressure': blk.gas_inlet.pressure[0].value,
        'mole_frac': {
            'CH4': blk.gas_inlet.mole_frac[0, 'CH4'].value,
            'CO2': blk.gas_inlet.mole_frac[0, 'CO2'].value,
            'H2O': blk.gas_inlet.mole_frac[0, 'H2O'].value,
            'CO': blk.gas_inlet.mole_frac[0, 'CO'].value,
            'H2': blk.gas_inlet.mole_frac[0, 'H2'].value}}
    solid_phase_state_args = {
        'flow_mass': blk.solid_inlet.flow_mass[0].value,
        'temperature': blk.solid_inlet.temperature[0].value,
        'mass_frac': {
            'Fe2O3': blk.solid_inlet.mass_frac[0, 'Fe2O3'].value,
            'Fe3O4': blk.solid_inlet.mass_frac[0, 'Fe3O4'].value,
            'Al2O3': blk.solid_inlet.mass_frac[0, 'Al2O3'].value}}

    # Reaction rates scaled down by a factor of 3 to enable simulation run.
    for x in m.fs.MB.length_domain:
        m.fs.MB.solid_phase.reactions[0, x]._params._scale_factor_rxn = 1e-3

    m.fs.MB.initialize(outlvl=5,
                       gas_phase_state_args=gas_phase_state_args,
                       solid_phase_state_args=solid_phase_state_args)

    t_initialize = time.time()  # Initialization time

    # Create a solver
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-8,
                      'linear_solver': 'ma57',
                      'bound_push': 1e-6,
                      'max_cpu_time': 3600,
                      'max_iter': 3000,
                      'mu_strategy': 'monotone',  # 'adaptive',#
                      'mu_init': 1e-8,
                      'print_level': 5,
                      'print_user_options': 'yes',
                      'file_print_level': 5,
                      'output_file': 'output.log',
                      'halt_on_ampl_error': 'no',
                      'warm_start_init_point': 'yes'}

    solver.solve(m.fs.MB, tee=True)

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


# %% # ------------------------------------------------------------------
    # Print and display some results
def print_summary(self):
    m.fs.MB.solid_phase.properties[0, 0].flow_mass.display()
    m.fs.MB.solid_phase.properties[0, 1].flow_mass.display()
    print()
    ebal_gas = ((m.fs.MB.gas_phase.properties[0, 0].flow_mol.value *
                 m.fs.MB.gas_phase.properties[0, 0].enth_mol.value) -
                (m.fs.MB.gas_phase.properties[0, 1].flow_mol.value *
                 m.fs.MB.gas_phase.properties[0, 1].enth_mol.value))

    ebal_solid = ((m.fs.MB.solid_phase.properties[0, 1].flow_mass.value *
                   m.fs.MB.solid_phase.properties[0, 1].enth_mass.value) -
                  (m.fs.MB.solid_phase.properties[0, 0].flow_mass.value *
                   m.fs.MB.solid_phase.properties[0, 0].enth_mass.value))
    ebal_tol = ebal_gas + ebal_solid

    print('Energy balance gas_FR:', ebal_gas)
    print('Energy balance solids_FR:', ebal_solid)
    print('Energy Balance Tolerance_FR:', ebal_tol)
    print()

    # Oxygen carrier and fuel conversion
    Conv_gas_CH4 = 1 - (
        (m.fs.MB.gas_phase.properties[0, 1].flow_mol.value *
         m.fs.MB.gas_phase.properties[0, 1].mole_frac['CH4'].value) /
        (m.fs.MB.gas_phase.properties[0, 0].flow_mol.value *
         m.fs.MB.gas_phase.properties[0, 0].mole_frac['CH4'].value))
    Conv_gas_CO = 1 - (
        (m.fs.MB.gas_phase.properties[0, 1].flow_mol.value *
         m.fs.MB.gas_phase.properties[0, 1].mole_frac['CO'].value) /
        (m.fs.MB.gas_phase.properties[0, 0].flow_mol.value *
         m.fs.MB.gas_phase.properties[0, 0].mole_frac['CO'].value))
    Conv_gas_H2 = 1 - (
        (m.fs.MB.gas_phase.properties[0, 1].flow_mol.value *
         m.fs.MB.gas_phase.properties[0, 1].mole_frac['H2'].value) /
        (m.fs.MB.gas_phase.properties[0, 0].flow_mol.value *
         m.fs.MB.gas_phase.properties[0, 0].mole_frac['H2'].value))
    Conv_OC = 1 - (
        (m.fs.MB.solid_phase.properties[0, 0].flow_mass.value *
         m.fs.MB.solid_phase.properties[0, 0].mass_frac['Fe2O3'].value) /
        (m.fs.MB.solid_phase.properties[0, 1].flow_mass.value *
         m.fs.MB.solid_phase.properties[0, 1].mass_frac['Fe2O3'].value))
    print("\nCH4 conversion:", value(Conv_gas_CH4)*100, " %")
    print("CO conversion:", value(Conv_gas_CO)*100, " %")
    print("H2 conversion:", value(Conv_gas_H2)*100, " %")
    print("\nFe2O3 conversion:", value(Conv_OC)*100, " %")


# %%       # ------------------------------------------------------------------
if __name__ == "__main__":
    m = main()
    print_summary(m)
