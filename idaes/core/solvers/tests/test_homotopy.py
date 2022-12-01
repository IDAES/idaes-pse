#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
IDAES Homotopy meta-solver tests.
"""

__author__ = "Andrew Lee"

import pytest

from pyomo.environ import ConcreteModel, Constraint, Param, TerminationCondition, Var

from idaes.core import FlowsheetBlock
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.solvers import get_solver
from idaes.core.solvers.homotopy import homotopy

# Set module level pyest marker
pytestmark = pytest.mark.solver

# -----------------------------------------------------------------------------
# Get default solver for testing
solver = get_solver()


@pytest.fixture()
def model():
    m = ConcreteModel()
    m.x = Var(initialize=1)
    m.y = Var(initialize=1)

    m.c = Constraint(expr=m.y == m.x**2)

    m.x.fix(10)

    return m


# -----------------------------------------------------------------------------
# Test argument validation
@pytest.mark.unit
def test_invalid_model(model):
    o = object()

    with pytest.raises(TypeError):
        homotopy(o, [model.x], [20])


@pytest.mark.unit
def test_not_var(model):
    model.p = Param()

    with pytest.raises(TypeError):
        homotopy(model, [model.p], [20])


@pytest.mark.unit
def test_var_not_in_model(model):
    m2 = ConcreteModel()

    with pytest.raises(ConfigurationError):
        homotopy(m2, [model.x], [20])


@pytest.mark.unit
def test_var_not_fixed(model):
    model.x.unfix()

    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20])


@pytest.mark.unit
def test_current_value_ub(model):
    model.x.setub(5)

    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20])


@pytest.mark.unit
def test_target_value_ub(model):
    model.x.setub(12)

    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20])


@pytest.mark.unit
def test_current_value_lb(model):
    model.x.setlb(12)

    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20])


@pytest.mark.unit
def test_target_value_lb(model):
    model.x.fix(50)
    model.x.setlb(30)

    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20])


@pytest.mark.unit
def test_step_init(model):
    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20], step_init=0.01)

    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20], step_init=0.9)


@pytest.mark.unit
def test_step_cut(model):
    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20], step_cut=0.09)

    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20], step_cut=0.91)


@pytest.mark.unit
def test_step_accel(model):
    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20], step_accel=-1)


@pytest.mark.unit
def test_iter_target(model):
    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20], iter_target=0)

    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20], iter_target=1.7)


@pytest.mark.unit
def test_max_step(model):
    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20], max_step=0.04)

    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20], max_step=1.1)


@pytest.mark.unit
def test_min_step(model):
    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20], min_step=0.009)

    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20], min_step=0.11)


@pytest.mark.unit
def test_min_max_step(model):
    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20], min_step=0.1, max_step=0.05)


@pytest.mark.unit
def test_init_step_min_max(model):
    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20], step_init=0.2, max_step=0.1)

    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20], step_init=0.05, min_step=0.1)


@pytest.mark.unit
def test_max_eval(model):
    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20], max_eval=0)

    with pytest.raises(ConfigurationError):
        homotopy(model, [model.x], [20], max_eval=1.7)


# -----------------------------------------------------------------------------
# Test termination conditions
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.unit
def test_basic(model):
    tc, prog, ni = homotopy(model, [model.x], [20])

    assert model.y.value == 400

    assert tc == TerminationCondition.optimal
    assert prog == 1
    assert ni == 4


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.unit
def test_basic_overshoot(model):
    # Use a big step such that overshoot will occur if not caught
    tc, prog, ni = homotopy(model, [model.x], [20], step_init=0.6)

    assert model.y.value == 400

    assert tc == TerminationCondition.optimal
    assert prog == 1
    assert ni == 2


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.unit
def test_basic_constraint_violation(model):
    # Add a constraint to limit y
    model.c2 = Constraint(expr=model.y <= 300)

    # Try to walk to a value of x that gives an infeasible value of y
    tc, prog, ni = homotopy(model, [model.x], [20])

    assert pytest.approx(model.y.value, 1e-4) == 293.27

    assert tc == TerminationCondition.minStepLength
    assert prog == 0.7125
    assert ni == 12


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.unit
def test_basic_max_iter(model):
    tc, prog, ni = homotopy(model, [model.x], [20], max_eval=2)

    assert pytest.approx(model.y.value, 1e-4) == 182.25

    assert tc == TerminationCondition.maxEvaluations
    assert prog == 0.35
    assert ni == 2


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.unit
def test_basic_infeasible_init(model):
    model.c2 = Constraint(expr=model.y <= 50)

    tc, prog, ni = homotopy(model, [model.x], [20])

    assert tc == TerminationCondition.infeasible
    assert prog == 0
    assert ni == 0


# TODO : need tests for convergence with regularisation
# -----------------------------------------------------------------------------
# Test that parameters have correct effect
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.unit
def test_basic_step_accel(model):
    # With zero acceleration, should take 10 steps
    tc, prog, ni = homotopy(model, [model.x], [20], step_accel=0)

    assert model.y.value == 400

    assert tc == TerminationCondition.optimal
    assert prog == 1
    assert ni == 10


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.unit
def test_basic_step_init(model):
    # With zero acceleration and initial step of 0.05, should take 20 steps
    tc, prog, ni = homotopy(model, [model.x], [20], step_init=0.05, step_accel=0)

    assert model.y.value == 400

    assert tc == TerminationCondition.optimal
    assert prog == 1
    assert ni == 20


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.unit
def test_basic_step_cut(model):
    # Add a constraint to limit y
    model.c2 = Constraint(expr=model.y <= 196)

    # Should take 6 steps
    # 4 steps to reach 14, and 2 to cut back to min_step
    tc, prog, ni = homotopy(
        model,
        [model.x],
        [20],
        step_init=0.1,
        min_step=0.025,
        step_cut=0.25,
        step_accel=0,
    )

    assert model.y.value == 196

    assert tc == TerminationCondition.minStepLength
    assert prog == 0.4
    assert ni == 6


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.unit
def test_basic_step_cut_2(model):
    # Add a constraint to limit y
    model.c2 = Constraint(expr=model.y <= 196)

    # Should take 7 steps
    # 4 steps to reach 14, and 3 to cut back to min_step
    tc, prog, ni = homotopy(
        model,
        [model.x],
        [20],
        step_init=0.1,
        min_step=0.01,
        step_cut=0.25,
        step_accel=0,
    )

    assert model.y.value == 196

    assert tc == TerminationCondition.minStepLength
    assert prog == 0.4
    assert ni == 7


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.unit
def test_basic_iter_target(model):
    # Should take 5 steps
    tc, prog, ni = homotopy(model, [model.x], [20], iter_target=2)

    assert model.y.value == 400

    assert tc == TerminationCondition.optimal
    assert prog == 1
    assert ni == 5


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.unit
def test_basic_max_step(model):
    # With max_step = step_init = 0.1, should take 10 steps
    tc, prog, ni = homotopy(model, [model.x], [20], max_step=0.1)

    assert model.y.value == 400

    assert tc == TerminationCondition.optimal
    assert prog == 1
    assert ni == 10


# -----------------------------------------------------------------------------
# Test a more complex problem
@pytest.fixture()
def model2():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)

    # vapor-liquid (ideal) - FTPz
    m.fs.properties_ideal_vl_FTPz = BTXParameterBlock(
        valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal", state_vars="FTPz"
    )
    m.fs.state_block = m.fs.properties_ideal_vl_FTPz.build_state_block(
        defined_state=True
    )

    m.fs.state_block.flow_mol.fix(1)
    m.fs.state_block.temperature.fix(360)
    m.fs.state_block.pressure.fix(101325)
    m.fs.state_block.mole_frac_comp["benzene"].fix(0.5)
    m.fs.state_block.mole_frac_comp["toluene"].fix(0.5)

    assert degrees_of_freedom(m.fs.state_block) == 0

    m.fs.state_block.initialize()

    return m


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.unit
def test_ideal_prop(model2):
    tc, prog, ni = homotopy(model2, [model2.fs.state_block.temperature], [390])

    assert tc == TerminationCondition.optimal
    assert prog == 1
    assert ni == 5

    # Check for VLE results
    assert model2.fs.state_block.mole_frac_phase_comp[
        "Liq", "benzene"
    ].value == pytest.approx(0.291, abs=1e-3)
    assert model2.fs.state_block.mole_frac_phase_comp[
        "Liq", "toluene"
    ].value == pytest.approx(0.709, abs=1e-3)
    assert model2.fs.state_block.mole_frac_phase_comp[
        "Vap", "benzene"
    ].value == pytest.approx(0.5, abs=1e-5)
    assert model2.fs.state_block.mole_frac_phase_comp[
        "Vap", "toluene"
    ].value == pytest.approx(0.5, abs=1e-5)


# Test max_iter here, as a more complicated model is needed
@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.unit
def test_ideal_prop_max_iter(model2):
    tc, prog, ni = homotopy(
        model2,
        [model2.fs.state_block.temperature],
        [390],
        max_solver_iterations=3,
        min_step=0.01,
        step_accel=0,
    )

    assert tc == TerminationCondition.optimal
    assert prog == 1
    assert ni == 19

    # Check for VLE results
    assert model2.fs.state_block.mole_frac_phase_comp[
        "Liq", "benzene"
    ].value == pytest.approx(0.291, abs=1e-3)
    assert model2.fs.state_block.mole_frac_phase_comp[
        "Liq", "toluene"
    ].value == pytest.approx(0.709, abs=1e-3)
    assert model2.fs.state_block.mole_frac_phase_comp[
        "Vap", "benzene"
    ].value == pytest.approx(0.5, abs=1e-5)
    assert model2.fs.state_block.mole_frac_phase_comp[
        "Vap", "toluene"
    ].value == pytest.approx(0.5, abs=1e-5)
