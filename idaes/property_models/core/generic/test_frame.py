#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 14:36:11 2019

@author: alee
"""

from idaes.core import FlowsheetBlock
from idaes.property_models.core.generic.test_params import TestParameterBlock

from pyomo.environ import ConcreteModel, SolverFactory


def main():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={'dynamic': False})

    m.fs.props = TestParameterBlock()

    m.fs.state = m.fs.props.state_block_class(default={'parameters': m.fs.props,
                                                       'defined_state': True})

    m.fs.state.flow_mol.fix(100)
    m.fs.state.mole_frac["benzene"].fix(0.5)
    m.fs.state.mole_frac["toluene"].fix(0.5)
    m.fs.state.temperature.fix(368)
    m.fs.state.pressure.fix(101325)

    m.fs.state.temperature_bubble.display()
    m.fs.state.temperature_dew.display()
    m.fs.state.pressure_bubble.display()
    m.fs.state.pressure_dew.display()
    m.fs.state.enth_mol.display()
    m.fs.state.entr_mol.display()
    m.fs.state.gibbs_mol.display()

#    m.fs.state.pprint()
#    print(m.fs.state.define_state_vars())
#    print(m.fs.state.get_material_flow_terms("Vap", "benzene"))

    solver = SolverFactory('ipopt')
    results = solver.solve(m, tee=True)

#    m.fs.state.temperature_bubble.display()
    m.fs.state.display()
    m.fs.state.mw_phase.display()
    m.fs.state.mw.display()
    m.fs.state.dens_mol_phase.display()
    m.fs.state.dens_mass_phase.display()
    m.fs.state.dens_mol.display()
    m.fs.state.dens_mass.display()
    m.fs.state.enth_mol_phase_comp.display()
    m.fs.state.enth_mol.display()
    m.fs.state.entr_mol_phase_comp.display()
    m.fs.state.entr_mol.display()
    m.fs.state.gibbs_mol_phase_comp.display()
    m.fs.state.gibbs_mol.display()
    m.fs.state.fug_coeff.display()
    m.fs.state.fug.display()


if __name__ == "__main__":
    main()
