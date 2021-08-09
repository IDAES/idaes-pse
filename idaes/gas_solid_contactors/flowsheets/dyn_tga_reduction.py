#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 15:40:33 2019

@author: alee
"""

import sys
import os

from pyomo.environ import ConcreteModel, SolverFactory, TransformationFactory

from idaes.core import (FlowsheetBlock, MaterialBalanceType, EnergyBalanceType,
        MomentumBalanceType)
#from idaes.core.util.initialization import initialize_by_time_element
from idaes.core.util.model_statistics import degrees_of_freedom

# Access parent directory (chemical_looping) of the current directory
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

from unit_models.fixed_bed_0D import FixedBed0D
import property_packages.methane_iron_OC_reduction.gas_phase_thermo as gas_props
import property_packages.methane_iron_OC_reduction.solid_phase_thermo as solid_props
import property_packages.methane_iron_OC_reduction.hetero_reactions as solid_rxns


m = ConcreteModel()

m.fs = FlowsheetBlock(default={"dynamic": True, "time_set": [0, 3600]})

m.fs.gas_props = gas_props.Gas_Phase_Thermo_ParameterBlock()
m.fs.solid_props = solid_props.Solid_Phase_Thermo_ParameterBlock()
m.fs.solid_rxns = solid_rxns.HeteroReactionParameterBlock(
        default={"solid_property_package": m.fs.solid_props,
                 "gas_property_package": m.fs.gas_props})
#        default={'property_package': m.fs.solid_props})

m.fs.TGA = FixedBed0D(default={"gas_property_package": m.fs.gas_props,
                               "solid_property_package": m.fs.solid_props,
                               "reaction_package": m.fs.solid_rxns,
                               
                               'material_balance_type': MaterialBalanceType.componentTotal,
                               'energy_balance_type': EnergyBalanceType.enthalpyTotal,
                               'momentum_balance_type': MomentumBalanceType.none})

# Discretize time domain
m.discretizer = TransformationFactory('dae.finite_difference')
m.discretizer.apply_to(m,
                       nfe=100,
                       wrt=m.fs.time,
                       scheme="BACKWARD")

# No idea what the proper conditions should be 
#m.fs.TGA.inlet.flow_mol.fix(1)
#m.fs.TGA.inlet.temperature.fix(303.15)
#m.fs.TGA.inlet.pressure.fix(2.0)

# Volume of reactor and initial mass of solids of a TGA are known
m.fs.TGA.volume_reactor.fix(0.2)
m.fs.TGA.mass_solids[0].fix(300)

## Fix value of bed_voidage
#m.fs.TGA.bed_voidage.fix(0.4) # Bed voidage for a fixed bed is usually known

# Fix solid material holdup to arbitrary values for now
#m.fs.TGA.solids_material_holdup[0, 'Fe2O3'].fix(4.5)
#m.fs.TGA.solids_material_holdup[0, 'Fe3O4'].fix(0)
#m.fs.TGA.solids_material_holdup[0, 'Al2O3'].fix(5.5)
#m.fs.TGA.solids[0].temperature.fix(303.15)

m.fs.TGA.solids[0].mass_frac['Fe2O3'].fix(0.45)
m.fs.TGA.solids[0].mass_frac['Fe3O4'].fix(0)
m.fs.TGA.solids[0].mass_frac['Al2O3'].fix(0.55)
m.fs.TGA.solids[0].temperature.fix(303.15)

# solids.flow_mass is included in the property package - need to remove
for t in m.fs.time:
#    m.fs.TGA.solids[t].flow_mass.fix(0)
#   ^ Does not seem to have an effect on the degrees of freedom.
#   Must not be present in any active constraints?
    m.fs.TGA.inlet.flow_mol[t].fix(1)
    m.fs.TGA.inlet.temperature[t].fix(303.15)
    m.fs.TGA.inlet.pressure[t].fix(2.0)

    m.fs.TGA.inlet.mole_frac[t, 'CO2'].fix(0.02499)
    m.fs.TGA.inlet.mole_frac[t, 'H2O'].fix(0.00001)
    m.fs.TGA.inlet.mole_frac[t, 'CH4'].fix(0.975)
    m.fs.TGA.gas_phase.heat[t].fix(0) 

# Set initial conditions - accumulation = 0 at time = 0
m.fs.fix_initial_conditions(state="steady-state")

#m.fs.TGA.initialize()

# Getting degree of freedom error in initialization. dof=101
print('dof:', degrees_of_freedom(m.fs))

solver = SolverFactory('ipopt')
#initialize_by_time_element(m.fs, m.fs.time, solver=solver)
solver.solve(m, tee=True)

#m.fs.TGA.outlet.display()
#m.fs.TGA.solids_material_holdup.display()
#m.fs.TGA.reactions[3600].display()
#m.fs.TGA.solids[0].display()
m.fs.TGA.solids[3600].display()
#m.fs.TGA.mass_solids.display()
