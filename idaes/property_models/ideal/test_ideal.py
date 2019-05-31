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
"""
Tests for ideal state block; tests for construction and solves
Author: Jaffer Ghouse
"""
import pytest
from pyomo.environ import ConcreteModel, SolverFactory, TerminationCondition, \
    SolverStatus, value

from idaes.core import FlowsheetBlock
from idaes.property_models.ideal.BTX_ideal_VLE import BTXParameterBlock
from idaes.core.util.model_statistics import calculate_degrees_of_freedom

# See if ipopt is available and set up solver
if SolverFactory('ipopt').available():
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8}
else:
    solver = None

# -----------------------------------------------------------------------------
m = ConcreteModel()
# Create a flowsheet object to test inlet state blocks
m.fs = FlowsheetBlock(default={"dynamic": False})

# vapor-liquid
m.fs.properties_vl = BTXParameterBlock(default={"valid_phase":
                                       ('Liq', 'Vap')})
m.fs.state_block_vl = m.fs.properties_vl.state_block_class(
    default={"parameters": m.fs.properties_vl,
             "defined_state": True})

# liquid only
m.fs.properties_l = BTXParameterBlock(default={"valid_phase": 'Liq'})
m.fs.state_block_l = m.fs.properties_l.state_block_class(
    default={"parameters": m.fs.properties_l,
             "has_phase_equilibrium": False,
             "defined_state": True})

# vapor only
m.fs.properties_v = BTXParameterBlock(default={"valid_phase": 'Vap'})
m.fs.state_block_v = m.fs.properties_v.state_block_class(
    default={"parameters": m.fs.properties_v,
             "has_phase_equilibrium": False,
             "defined_state": True})


def test_build_inlet_state_blocks():
    assert len(m.fs.properties_vl.config) == 2

    # vapor-liquid
    assert m.fs.properties_vl.config.valid_phase == ('Vap', 'Liq') or \
        m.fs.properties_vl.config.valid_phase == ('Liq', 'Vap')
    assert len(m.fs.properties_vl.phase_list) == 2
    assert m.fs.properties_vl.phase_list == ["Liq", "Vap"]
    assert hasattr(m.fs.state_block_vl, "equilibrium_constraint")
    assert not hasattr(m.fs.state_block_vl, "eq_mol_frac_out")

    # liquid only
    assert m.fs.properties_l.config.valid_phase == "Liq"
    assert len(m.fs.properties_l.phase_list) == 1
    assert m.fs.properties_l.phase_list == ["Liq"]
    assert not hasattr(m.fs.state_block_l, "equilibrium_constraint")
    assert not hasattr(m.fs.state_block_vl, "eq_h_vap")
    assert not hasattr(m.fs.state_block_l, "eq_mol_frac_out")

    # vapor only
    assert m.fs.properties_v.config.valid_phase == "Vap"
    assert len(m.fs.properties_v.phase_list) == 1
    assert m.fs.properties_v.phase_list == ["Vap"]
    assert not hasattr(m.fs.state_block_v, "equilibrium_constraint")
    assert not hasattr(m.fs.state_block_vl, "eq_h_liq")
    assert not hasattr(m.fs.state_block_v, "eq_mol_frac_out")


def test_setInputs_inlet_state_blocks():

    # vapor-liquid
    m.fs.state_block_vl.flow_mol.fix(1)
    m.fs.state_block_vl.temperature.fix(368)
    m.fs.state_block_vl.pressure.fix(101325)
    m.fs.state_block_vl.mole_frac["benzene"].fix(0.5)
    m.fs.state_block_vl.mole_frac["toluene"].fix(0.5)

    assert calculate_degrees_of_freedom(m.fs.state_block_vl) == 0

    # liquid only
    m.fs.state_block_l.flow_mol.fix(1)
    m.fs.state_block_l.temperature.fix(362)
    m.fs.state_block_l.pressure.fix(101325)
    m.fs.state_block_l.mole_frac["benzene"].fix(0.5)
    m.fs.state_block_l.mole_frac["toluene"].fix(0.5)

    assert calculate_degrees_of_freedom(m.fs.state_block_l) == 0

    # vapor only
    m.fs.state_block_v.flow_mol.fix(1)
    m.fs.state_block_v.temperature.fix(375)
    m.fs.state_block_v.pressure.fix(101325)
    m.fs.state_block_v.mole_frac["benzene"].fix(0.5)
    m.fs.state_block_v.mole_frac["toluene"].fix(0.5)

    assert calculate_degrees_of_freedom(m.fs.state_block_v) == 0


@pytest.mark.skipif(solver is None, reason="Solver not available")
def test_solve():
    # vapor-liquid
    m.fs.state_block_vl.initialize()
    results = solver.solve(m.fs.state_block_vl, tee=False)

    # Check for optimal solution
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    # Check for VLE results
    assert value(m.fs.state_block_vl.mole_frac_phase['Liq', 'benzene']) == \
        pytest.approx(0.4215, abs=1e-3)
    assert value(m.fs.state_block_vl.mole_frac_phase['Vap', 'benzene']) == \
        pytest.approx(0.643, abs=1e-3)

    # liquid only
    m.fs.state_block_l.initialize()
    results = solver.solve(m.fs.state_block_l, tee=False)

    # Check for optimal solution
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    # Check for results
    assert value(m.fs.state_block_l.mole_frac_phase['Liq', 'benzene']) == \
        pytest.approx(0.5, abs=1e-3)
    assert value(m.fs.state_block_l.mole_frac_phase['Liq', 'toluene']) == \
        pytest.approx(0.5, abs=1e-3)

    # vapor only
    m.fs.state_block_v.initialize()
    results = solver.solve(m.fs.state_block_v, tee=False)

    # Check for optimal solution
    assert results.solver.termination_condition == TerminationCondition.optimal
    assert results.solver.status == SolverStatus.ok

    # Check for results
    assert value(m.fs.state_block_v.mole_frac_phase['Vap', 'benzene']) == \
        pytest.approx(0.5, abs=1e-3)
    assert value(m.fs.state_block_v.mole_frac_phase['Vap', 'toluene']) == \
        pytest.approx(0.5, abs=1e-3)


def test_bubbleT_inlet_state_blocks():
    assert m.fs.state_block_vl.calculate_bubble_point_temperature() == \
        pytest.approx(365.64, abs=1e-2)


def test_dewT_inlet_state_blocks():
    assert m.fs.state_block_vl.calculate_dew_point_temperature() == \
        pytest.approx(372.31, abs=1e-2)


def test_bubbleP_inlet_state_blocks():
    assert m.fs.state_block_vl.calculate_bubble_point_pressure() == \
        pytest.approx(108600, abs=1e3)


def test_dewP_inlet_state_blocks():
    assert m.fs.state_block_vl.calculate_dew_point_pressure() == \
        pytest.approx(89000, abs=1e3)


# Create a flowsheet object to test outlet state blocks
m.fs1 = FlowsheetBlock(default={"dynamic": False})
# vapor-liquid
m.fs1.properties_vl = BTXParameterBlock(default={"valid_phase":
                                                 ('Liq', 'Vap')})
m.fs1.state_block_vl = m.fs1.properties_vl.state_block_class(
    default={"parameters": m.fs1.properties_vl,
             "defined_state": False})

# liquid only
m.fs1.properties_l = BTXParameterBlock(default={"valid_phase": 'Liq'})
m.fs1.state_block_l = m.fs1.properties_l.state_block_class(
    default={"parameters": m.fs1.properties_l,
             "has_phase_equilibrium": False,
             "defined_state": False})

# vapor only
m.fs1.properties_v = BTXParameterBlock(default={"valid_phase": 'Vap'})
m.fs1.state_block_v = m.fs1.properties_v.state_block_class(
    default={"parameters": m.fs1.properties_v,
             "has_phase_equilibrium": False,
             "defined_state": False})


def test_build_outlet_state_blocks():
    assert len(m.fs1.properties_vl.config) == 2

    # vapor-liquid
    assert m.fs1.properties_vl.config.valid_phase == ('Vap', 'Liq') or \
        m.fs1.properties_vl.config.valid_phase == ('Liq', 'Vap')
    assert len(m.fs1.properties_vl.phase_list) == 2
    assert m.fs1.properties_vl.phase_list == ["Liq", "Vap"]
    assert hasattr(m.fs1.state_block_vl, "equilibrium_constraint")
    assert hasattr(m.fs1.state_block_vl, "sum_mole_frac_out")

    # liquid only
    assert m.fs1.properties_l.config.valid_phase == "Liq"
    assert len(m.fs1.properties_l.phase_list) == 1
    assert m.fs1.properties_l.phase_list == ["Liq"]
    assert not hasattr(m.fs1.state_block_l, "equilibrium_constraint")
    assert not hasattr(m.fs1.state_block_vl, "eq_h_vap")
    assert hasattr(m.fs1.state_block_l, "sum_mole_frac_out")

    # vapor only
    assert m.fs1.properties_v.config.valid_phase == "Vap"
    assert len(m.fs1.properties_v.phase_list) == 1
    assert m.fs1.properties_v.phase_list == ["Vap"]
    assert not hasattr(m.fs1.state_block_v, "equilibrium_constraint")
    assert not hasattr(m.fs1.state_block_vl, "eq_h_liq")
    assert hasattr(m.fs1.state_block_v, "sum_mole_frac_out")


def test_setInputs_outlet_state_blocks():

    # vapor-liquid
    m.fs1.state_block_vl.flow_mol.fix(1)
    m.fs1.state_block_vl.temperature.fix(368)
    m.fs1.state_block_vl.pressure.fix(101325)
    m.fs1.state_block_vl.mole_frac["benzene"].fix(0.5)

    assert calculate_degrees_of_freedom(m.fs.state_block_vl) == 0

    # liquid only
    m.fs1.state_block_l.flow_mol.fix(1)
    m.fs1.state_block_l.temperature.fix(362)
    m.fs1.state_block_l.pressure.fix(101325)
    m.fs1.state_block_l.mole_frac["benzene"].fix(0.5)

    assert calculate_degrees_of_freedom(m.fs.state_block_l) == 0

    # vapor only
    m.fs1.state_block_v.flow_mol.fix(1)
    m.fs1.state_block_v.temperature.fix(375)
    m.fs1.state_block_v.pressure.fix(101325)
    m.fs1.state_block_v.mole_frac["benzene"].fix(0.5)

    assert calculate_degrees_of_freedom(m.fs.state_block_v) == 0
