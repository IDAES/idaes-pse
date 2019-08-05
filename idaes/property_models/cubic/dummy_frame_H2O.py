#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 14:36:11 2019

@author: alee
"""

from idaes.core import FlowsheetBlock
from idaes.property_models.cubic import H2O_PR_VLE

from pyomo.environ import ConcreteModel, SolverFactory


def main():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={'dynamic': False})

    m.fs.props = H2O_PR_VLE.H2OParameterBlock(
            default={'valid_phase': ('Vap', 'Liq')})

    m.fs.state = m.fs.props.state_block_class(default={'parameters': m.fs.props})

    m.fs.state.flow_mol.fix(100)
    m.fs.state.temperature.fix(298.15)
    m.fs.state.pressure.fix(101325)
    m.fs.state.enth_mol

    m.fs.state.mole_frac["H2O"].value = 1

    m.fs.state.initialize()

    # Create a solver
    solver = SolverFactory('ipopt')
    results = solver.solve(m, tee=True)

    # -----------------------------------------------------------------------------
    m.fs.state.temperature.unfix()
    m.fs.state.enth_mol.fix()

    results = solver.solve(m, tee=True)
    #
    #m.fs.state.enth_mol.fix(-45000)
    #
    ## Create a solver
    #results = solver.solve(m, tee=True)

    # Print results
    print(results)

    m.fs.state.flow_mol_phase.display()
    m.fs.state.mole_frac_phase.display()
    m.fs.state.enth_mol.display()
    m.fs.state.enth_mol_phase.display()
    m.fs.state.temperature.display()

    m.fs.state.temperature_dew.display()
    m.fs.state.temperature_bubble.display()
    m.fs.state._teq.display()
    m.fs.state._fug_phase_eq.display()


    #m.fs.state.enth_mol.fix(-40000)
    #
    ## Create a solver
    #results = solver.solve(m, tee=True)
    #
    ## Print results
    #print(results)
    #
    #m.fs.state.flow_mol_phase.display()
    #m.fs.state.mole_frac_phase.display()
    #m.fs.state.enth_mol.display()
    #m.fs.state.enth_mol_phase.display()
    #m.fs.state.temperature.display()
    #
    #m.fs.state.temperature_dew.display()
    #m.fs.state.temperature_bubble.display()
    #m.fs.state._teq.display()
    #m.fs.state._fug_phase_eq.display()


if __name__ == "__main__":
    main()
