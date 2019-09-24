#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 14:36:11 2019

@author: alee
"""

from idaes.core import FlowsheetBlock
from idaes.property_models.cubic import BT_PR_VLE_A, cubic_prop_pack_VLE

from pyomo.environ import ConcreteModel, SolverFactory, value


def main():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={'dynamic': False})

    m.fs.props = BT_PR_VLE_A.BTParameterBlock(
            default={'valid_phase': ('Vap', 'Liq')})

    m.fs.state = m.fs.props.state_block_class(default={'parameters': m.fs.props,
                                                       'defined_state': True})

    F = 100
    x_benzene = 0.5
    T = 400
    P = 1e5

    # -----------------------------------------------------------------------------
    m.fs.state.flow_mol.fix(F)
    m.fs.state.mole_frac["benzene"].fix(x_benzene)
    m.fs.state.mole_frac["toluene"].fix(1-x_benzene)
    m.fs.state.temperature.fix(T)
    m.fs.state.pressure.fix(P)

    # -----------------------------------------------------------------------------
    m.fs.state.enth_mol
    m.fs.state.entr_mol

    m.fs.state.initialize(outlvl=5)

    # -----------------------------------------------------------------------------
    # Create a solver
    solver = SolverFactory('ipopt')
    results = solver.solve(m, tee=True)

    # Print results
    print(results)

    print()
    print("Cubic")
#    m.fs.state.temperature_dew.display()
#    m.fs.state.temperature_bubble.display()
    m.fs.state._teq.display()
#    m.fs.state.dens_mol_phase.display()
#    m.fs.state.dens_mass_phase.display()

    m.fs.state.compress_fact.display()
    m.fs.state.enth_mol_phase.display()
    m.fs.state.entr_mol_phase.display()
    m.fs.state.mole_frac_phase.display()
    m.fs.state.fug_phase.display()
    print("Vapor Fraction:", value(m.fs.state.flow_mol_phase["Vap"]/m.fs.state.flow_mol))

    print()
    print(value(cubic_prop_pack_VLE.CubicStateBlockData._enth_mol_comp_ig(m.fs.state, 'benzene')),
          value(cubic_prop_pack_VLE.CubicStateBlockData._enth_mol_comp_ig(m.fs.state, 'toluene')))
    m.fs.state.dadT.display()
    print(value(cubic_prop_pack_VLE.CubicStateBlockData._enth_mol_cubic(m.fs.state, 'Vap')))


if __name__ == "__main__":
    main()
