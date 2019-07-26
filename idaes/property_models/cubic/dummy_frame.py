#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 14:36:11 2019

@author: alee
"""

from idaes.core import FlowsheetBlock
from idaes.property_models.cubic import ASU_PR_VLE, cubic_prop_pack_VLE

from pyomo.environ import ConcreteModel, SolverFactory


m = ConcreteModel()

m.fs = FlowsheetBlock(default={'dynamic': False})

m.fs.props = ASU_PR_VLE.ASUParameterBlock(default={'valid_phase': 'Vap'})

#m.fs.props.pprint()

m.fs.state = m.fs.props.state_block_class(default={'parameters': m.fs.props})

m.fs.state.flow_mol.fix(100)
m.fs.state.mole_frac["N2"].fix(0.6)
m.fs.state.mole_frac["O2"].fix(0.35)
#m.fs.state.mole_frac["Ar"].fix(0.05)
m.fs.state.temperature.fix(98.078)
m.fs.state.pressure.fix(500000)

m.fs.state.mole_frac["Ar"] = 0.05

m.fs.state.mole_frac_phase["Vap", "N2"] = 0.6
m.fs.state.mole_frac_phase["Vap", "O2"] = 0.35
m.fs.state.mole_frac_phase["Vap", "Ar"] = 0.05

print(cubic_prop_pack_VLE.cubic_roots_available())

m.fs.state.enth_mol_phase.display()

# Create a solver
solver = SolverFactory('ipopt')
results = solver.solve(m, tee=True)

# Print results
print(results)

#m.fs.state.pprint()
m.fs.state.compress_fact_vap.display()
m.fs.state.fug_coeff_phase.display()
m.fs.state.fug_phase.display()
m.fs.state.dens_mol_phase.display()
m.fs.state.mw_phase.display()
m.fs.state.mw.display()
m.fs.state.enth_mol_phase.display()
