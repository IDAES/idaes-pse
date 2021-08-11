#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 17 15:40:33 2019

@author: alee
"""

import sys
import os

from pyomo.environ import ConcreteModel, SolverFactory, TransformationFactory

from idaes.core import FlowsheetBlock
from idaes.core.util.initialization import initialize_by_time_element
from idaes.core.util import get_solver

# Access parent directory (chemical_looping) of the current directory
sys.path.append(os.path.join(os.path.dirname(__file__), os.pardir))

from idaes.gas_solid_contactors.unit_models.fixed_bed_0D import FixedBed0D
import property_packages.TGA_example_properties.TGA_gas_example as gas_props
import property_packages.TGA_example_properties.TGA_solid_example as solid_props
import property_packages.TGA_example_properties.TGA_reactions_example as solid_rxns


m = ConcreteModel()

m.fs = FlowsheetBlock(default={"dynamic": True, "time_set": [0, 3600]})

m.fs.gas_props = gas_props.AirParameterBlock()
m.fs.solid_props = solid_props.ZeoliteParameterBlock()
m.fs.solid_rxns = solid_rxns.AdsorptionParameterBlock(
        default={"solid_property_package": m.fs.solid_props,
                 "gas_property_package": m.fs.gas_props})

m.fs.TGA = FixedBed0D(default={"gas_property_package": m.fs.gas_props,
                               "solid_property_package": m.fs.solid_props,
                               "reaction_package": m.fs.solid_rxns})

# Discretize time domain
m.discretizer = TransformationFactory('dae.finite_difference')
m.discretizer.apply_to(m,
                       nfe=100,
                       wrt=m.fs.time,
                       scheme="BACKWARD")

m.fs.TGA.inlet.flow_mol.fix(1)
m.fs.TGA.inlet.temperature.fix(303.15)
m.fs.TGA.inlet.pressure.fix(1.2e5)

m.fs.TGA.volume_reactor.fix(0.2)

m.fs.TGA.solids_material_holdup[0, "Zeo"].fix(0.1)
m.fs.TGA.solids[0].loading["N2"].fix(0.0)
m.fs.TGA.solids[0].temperature.fix(303.15)

for t in m.fs.time:
    m.fs.TGA.inlet.mole_frac_comp[t, "N2"].fix(0.79)
    m.fs.TGA.inlet.mole_frac_comp[t, "O2"].fix(0.21)

# Set initial conditions - accumulation = 0 at time = 0
m.fs.fix_initial_conditions(state="steady-state")

#m.fs.TGA.initialize()

opt = get_solver(solver, optarg) # create solver
initialize_by_time_element(m.fs, m.fs.time, solver=solver)
solver.solve(m, tee=True)

#m.fs.TGA.outlet.display()
#m.fs.TGA.solids_material_holdup.display()
#m.fs.TGA.reactions[3600].display()
#m.fs.TGA.solids[0].display()
m.fs.TGA.solids[3600].display()
#m.fs.TGA.mass_solids.display()
import pdb; pdb.set_trace()
