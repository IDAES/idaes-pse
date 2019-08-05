#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 14:36:11 2019

@author: alee
"""

from idaes.core import FlowsheetBlock
from idaes.property_models.cubic import ASU_PR_VLE, cubic_prop_pack_VLE

from pyomo.environ import ConcreteModel, SolverFactory


def main():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={'dynamic': False})

    m.fs.props = ASU_PR_VLE.ASUParameterBlock(
            default={'valid_phase': ('Vap', 'Liq')})
    #m.fs.props = ASU_PR_VLE.ASUParameterBlock(default={'valid_phase': 'Liq'})

    #m.fs.props.pprint()

    m.fs.state = m.fs.props.state_block_class(default={'parameters': m.fs.props})

    m.fs.state.flow_mol.fix(100)
    m.fs.state.mole_frac["N2"].fix(0.65)
    m.fs.state.mole_frac["O2"].fix(0.3)
    #m.fs.state.mole_frac["Ar"].fix(0.05)
    m.fs.state.temperature.fix(200)
    m.fs.state.pressure.fix(101000)

    m.fs.state.mole_frac["Ar"].value = 0.05

    print(cubic_prop_pack_VLE.cubic_roots_available())

    m.fs.state.initialize()

    # Create a solver
    solver = SolverFactory('ipopt')
    results = solver.solve(m, tee=True)

    # Print results
    print(results)

    m.fs.state.flow_mol_phase.display()
    m.fs.state.mole_frac_phase.display()

    m.fs.state.temperature_dew.display()
    m.fs.state.temperature_bubble.display()
    m.fs.state._teq.display()


if __name__ == "__main__":
    main()
