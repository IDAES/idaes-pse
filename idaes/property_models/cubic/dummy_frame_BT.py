#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 14:36:11 2019

@author: alee
"""

from idaes.core import FlowsheetBlock
from idaes.property_models.cubic import BT_PR

from pyomo.environ import ConcreteModel, SolverFactory, value


def main():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={'dynamic': False})

    m.fs.props = BT_PR.BTParameterBlock(
            default={'valid_phase': ('Vap', 'Liq')})

    m.fs.state = m.fs.props.state_block_class(default={'parameters': m.fs.props,
                                                       'defined_state': True})

    F = 100
    x_benzene = 0.5
    T = 350
    P = 5e5

    # -----------------------------------------------------------------------------
    m.fs.state.flow_mol.fix(F)
    m.fs.state.mole_frac_comp["benzene"].fix(x_benzene)
    m.fs.state.mole_frac_comp["toluene"].fix(1-x_benzene)
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
    m.fs.state.temperature_dew.display()
    m.fs.state.temperature_bubble.display()
    m.fs.state._teq.display()
#    m.fs.state.dens_mol_phase.display()
#    m.fs.state.dens_mass_phase.display()

    m.fs.state.compress_fact_phase.display()
    m.fs.state.enth_mol_phase.display()
    m.fs.state.entr_mol_phase.display()
    m.fs.state.mole_frac_phase_comp.display()
    m.fs.state.fug_phase_comp.display()
    m.fs.state.fug_coeff_phase_comp.display()
    print("Vapor Fraction:", value(m.fs.state.flow_mol_phase["Vap"]/m.fs.state.flow_mol))
    print()
#    m.fs.state.delta.display()
#    m.fs.state.a.display()
#    m.fs.state.am.display()
#    m.fs.state.b.display()
#    m.fs.state.bm.display()
#    m.fs.state.A.display()
#    m.fs.state.B.display()
#    m.fs.state.compress_fact_phase.display()


if __name__ == "__main__":
    main()
