#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 14:36:11 2019

@author: alee
"""

from idaes.core import FlowsheetBlock
from idaes.property_models.cubic import BT_PR_VLE
from idaes.property_models.ideal import BTX_ideal_VLE

from pyomo.environ import ConcreteModel, SolverFactory, value


m = ConcreteModel()

m.fs = FlowsheetBlock(default={'dynamic': False})

m.fs.props = BT_PR_VLE.BTParameterBlock(
        default={'valid_phase': ('Vap', 'Liq')})
m.fs.props_ideal = BTX_ideal_VLE.BTXParameterBlock(
        default={'valid_phase': ('Vap', 'Liq')})

m.fs.state = m.fs.props.state_block_class(default={'parameters': m.fs.props})
m.fs.ideal = m.fs.props_ideal.state_block_class(
        default={'parameters': m.fs.props_ideal})

# -----------------------------------------------------------------------------
m.fs.state.flow_mol.fix(100)
m.fs.state.mole_frac["benzene"].fix(0.5)
m.fs.state.temperature.fix(368)
m.fs.state.pressure.fix(101325)

m.fs.ideal.flow_mol.fix(100)
m.fs.ideal.mole_frac["benzene"].fix(0.5)
m.fs.ideal.temperature.fix(368)
m.fs.ideal.pressure.fix(101325)

# -----------------------------------------------------------------------------
m.fs.state.enth_mol
m.fs.state.entr_mol

m.fs.ideal.enth_mol_phase_comp

m.fs.state.mole_frac["toluene"] = 0.5

m.fs.state.initialize()
m.fs.ideal.initialize()

# -----------------------------------------------------------------------------
# Create a solver
solver = SolverFactory('ipopt')
results = solver.solve(m, tee=True)

# Print results
print(results)

print()
print("Cubic")
m.fs.state.flow_mol_phase.display()
m.fs.state.mole_frac_phase.display()
m.fs.state.enth_mol.display()
m.fs.state.enth_mol_phase.display()
m.fs.state.temperature.display()

m.fs.state.temperature_dew.display()
m.fs.state.temperature_bubble.display()
m.fs.state._teq.display()
m.fs.state._fug_phase_eq.display()

print()
for j in m.fs.props.component_list:
    print(j, value(m.fs.state._enth_mol_comp_ig(j)))

print()
print("ideal")
for j in m.fs.props.component_list:
    print(j, value(m.fs.ideal.enth_mol_phase_comp["Vap", j] - m.fs.ideal.dh_vap[j]))
