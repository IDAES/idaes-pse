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
import pytest

from idaes.core import FlowsheetBlock
from idaes.property_models.cubic import BT_PR

from pyomo.environ import (ConcreteModel,
                           Objective,
                           SolverFactory,
                           TerminationCondition,
                           value)


def test_T_sweep():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={'dynamic': False})

    m.fs.props = BT_PR.BTParameterBlock(
            default={'valid_phase': ('Vap', 'Liq')})

    m.fs.state = m.fs.props.state_block_class(
            default={'parameters': m.fs.props,
                     'defined_state': True})

    m.fs.obj = Objective(expr=(m.fs.state.temperature - 510)**2)

    for logP in range(8, 13, 1):
        m.fs.obj.deactivate()

        m.fs.state.flow_mol.fix(100)
        m.fs.state.mole_frac_comp["benzene"].fix(0.5)
        m.fs.state.mole_frac_comp["toluene"].fix(0.5)
        m.fs.state.temperature.fix(300)
        m.fs.state.pressure.fix(10**(0.5*logP))

        m.fs.state.initialize(outlvl=0)

        m.fs.state.temperature.unfix()
        m.fs.obj.activate()

        solver = SolverFactory('ipopt')
        results = solver.solve(m, tee=True)

        assert results.solver.termination_condition == \
            TerminationCondition.optimal
        assert m.fs.state.flow_mol_phase["Liq"].value <= 1e-5


def test_P_sweep():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={'dynamic': False})

    m.fs.props = BT_PR.BTParameterBlock(
            default={'valid_phase': ('Vap', 'Liq')})

    m.fs.state = m.fs.props.state_block_class(
            default={'parameters': m.fs.props,
                     'defined_state': True})

    for T in range(370, 500, 25):
        m.fs.state.flow_mol.fix(100)
        m.fs.state.mole_frac_comp["benzene"].fix(0.5)
        m.fs.state.mole_frac_comp["toluene"].fix(0.5)
        m.fs.state.temperature.fix(T)
        m.fs.state.pressure.fix(1e5)

        m.fs.state.initialize(outlvl=0)

        solver = SolverFactory('ipopt')
        results = solver.solve(m)

        assert results.solver.termination_condition == \
            TerminationCondition.optimal

        while m.fs.state.pressure.value <= 1e6:
            m.fs.state.pressure.value = m.fs.state.pressure.value + 1e5
            solver = SolverFactory('ipopt')
            results = solver.solve(m)
            assert results.solver.termination_condition == \
                TerminationCondition.optimal
            print(T, m.fs.state.pressure.value)


def test_T350_P1_x5():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={'dynamic': False})

    m.fs.props = BT_PR.BTParameterBlock(
            default={'valid_phase': ('Vap', 'Liq')})

    m.fs.state = m.fs.props.state_block_class(
            default={'parameters': m.fs.props,
                     'defined_state': True})

    m.fs.state.flow_mol.fix(100)
    m.fs.state.mole_frac_comp["benzene"].fix(0.5)
    m.fs.state.mole_frac_comp["toluene"].fix(0.5)
    m.fs.state.temperature.fix(350)
    m.fs.state.pressure.fix(1e5)

    # Trigger build of enthalpy and entropy
    m.fs.state.enth_mol_phase
    m.fs.state.entr_mol_phase

    m.fs.state.initialize(outlvl=0)

    solver = SolverFactory('ipopt')
    solver.solve(m)

    assert pytest.approx(value(m.fs.state._teq), abs=1e-1) == 365
    assert pytest.approx(
            value(m.fs.state.compress_fact_phase["Liq"]), 1e-5) == 0.0035346
    assert pytest.approx(
            value(m.fs.state.compress_fact_phase["Vap"]), 1e-5) == 0.966749
    assert pytest.approx(
            value(m.fs.state.fug_coeff_phase_comp["Liq", "benzene"]),
            1e-5) == 0.894676
    assert pytest.approx(
            value(m.fs.state.fug_coeff_phase_comp["Liq", "toluene"]),
            1e-5) == 0.347566
    assert pytest.approx(
            value(m.fs.state.fug_coeff_phase_comp["Vap", "benzene"]),
            1e-5) == 0.971072
    assert pytest.approx(
            value(m.fs.state.fug_coeff_phase_comp["Vap", "toluene"]),
            1e-5) == 0.959791

    assert pytest.approx(
            value(m.fs.state.mole_frac_phase_comp["Liq", "benzene"]),
            1e-5) == 0.5
    assert pytest.approx(
            value(m.fs.state.mole_frac_phase_comp["Liq", "toluene"]),
            1e-5) == 0.5
    assert pytest.approx(
            value(m.fs.state.mole_frac_phase_comp["Vap", "benzene"]),
            1e-5) == 0.70584
    assert pytest.approx(
            value(m.fs.state.mole_frac_phase_comp["Vap", "toluene"]),
            1e-5) == 0.29416

    assert pytest.approx(
            value(m.fs.state.enth_mol_phase["Liq"]), 1e-5) == 38942.8
    assert pytest.approx(
            value(m.fs.state.enth_mol_phase["Vap"]), 1e-5) == 78048.7
    assert pytest.approx(
            value(m.fs.state.entr_mol_phase["Liq"]), 1e-5) == -367.558
    assert pytest.approx(
            value(m.fs.state.entr_mol_phase["Vap"]), 1e-5) == -269.0553


def test_T350_P5_x5():
    m = ConcreteModel()

    m.fs = FlowsheetBlock(default={'dynamic': False})

    m.fs.props = BT_PR.BTParameterBlock(
            default={'valid_phase': ('Vap', 'Liq')})

    m.fs.state = m.fs.props.state_block_class(
            default={'parameters': m.fs.props,
                     'defined_state': True})

    m.fs.state.flow_mol.fix(100)
    m.fs.state.mole_frac_comp["benzene"].fix(0.5)
    m.fs.state.mole_frac_comp["toluene"].fix(0.5)
    m.fs.state.temperature.fix(350)
    m.fs.state.pressure.fix(5e5)

    # Trigger build of enthalpy and entropy
    m.fs.state.enth_mol_phase
    m.fs.state.entr_mol_phase

    m.fs.state.initialize(outlvl=0)

    solver = SolverFactory('ipopt')
    solver.solve(m)

    assert pytest.approx(value(m.fs.state._teq), 1e-5) == 431.47
    assert pytest.approx(
            value(m.fs.state.compress_fact_phase["Liq"]), 1e-5) == 0.01766
    assert pytest.approx(
            value(m.fs.state.compress_fact_phase["Vap"]), 1e-5) == 0.80245
    assert pytest.approx(
            value(m.fs.state.fug_coeff_phase_comp["Liq", "benzene"]),
            1e-5) == 0.181229
    assert pytest.approx(
            value(m.fs.state.fug_coeff_phase_comp["Liq", "toluene"]),
            1e-5) == 0.070601
    assert pytest.approx(
            value(m.fs.state.fug_coeff_phase_comp["Vap", "benzene"]),
            1e-5) == 0.856523
    assert pytest.approx(
            value(m.fs.state.fug_coeff_phase_comp["Vap", "toluene"]),
            1e-5) == 0.799237

    assert pytest.approx(
            value(m.fs.state.mole_frac_phase_comp["Liq", "benzene"]),
            1e-5) == 0.5
    assert pytest.approx(
            value(m.fs.state.mole_frac_phase_comp["Liq", "toluene"]),
            1e-5) == 0.5
    assert pytest.approx(
            value(m.fs.state.mole_frac_phase_comp["Vap", "benzene"]),
            1e-5) == 0.65415
    assert pytest.approx(
            value(m.fs.state.mole_frac_phase_comp["Vap", "toluene"]),
            1e-5) == 0.34585

    assert pytest.approx(
            value(m.fs.state.enth_mol_phase["Liq"]), 1e-5) == 38966.9
    assert pytest.approx(
            value(m.fs.state.enth_mol_phase["Vap"]), 1e-5) == 75150.7
    assert pytest.approx(
            value(m.fs.state.entr_mol_phase["Liq"]), 1e-5) == -367.6064
    assert pytest.approx(
            value(m.fs.state.entr_mol_phase["Vap"]), 1e-5) == -287.3318
