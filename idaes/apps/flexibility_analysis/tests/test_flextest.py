#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
from idaes.apps.flexibility_analysis import _check_dependencies
import pyomo.environ as pe
from idaes.apps.flexibility_analysis.flextest import (
    build_active_constraint_flextest,
    ActiveConstraintConfig,
)
import unittest
from idaes.apps.flexibility_analysis.indices import _VarIndex
from pyomo.contrib.fbbt import interval
import pytest
from idaes.core.util.testing import _enable_scip_solver_for_testing
from contextlib import contextmanager


@contextmanager
def scip_solver():
    solver = pe.SolverFactory("scip")
    undo_changes = None

    if not solver.available():
        undo_changes = _enable_scip_solver_for_testing()
    if not solver.available():
        pytest.skip(reason="SCIP solver not available")
    yield solver
    if undo_changes is not None:
        undo_changes()


def create_poly_model():
    m = pe.ConcreteModel()
    m.z = pe.Var(initialize=8, bounds=(-10, 15))
    m.theta = pe.Param(mutable=True)

    offset = 1.5

    m.obj = pe.Objective(expr=m.z**2)
    m.c1 = pe.Constraint(
        expr=0.01 * (m.z - offset) ** 4
        - 0.05 * (m.z - offset) ** 3
        - (m.z - offset) ** 2
        - (m.z - offset)
        - 10
        + m.theta
        <= 0
    )
    m.c2 = pe.Constraint(
        expr=(-0.02 * m.theta - 14) * m.z + (1.66 * m.theta - 100) <= 0
    )
    m.c3 = pe.Constraint(
        expr=30 * m.z - 50 - 4 * m.theta + pe.exp(-0.2 * m.theta + 1) <= 0
    )

    nominal_values = pe.ComponentMap()
    nominal_values[m.theta] = 22.5
    param_bounds = pe.ComponentMap()
    param_bounds[m.theta] = (-20, 65)

    return m, nominal_values, param_bounds


def create_hx_network_model():
    """
    The model used for this test comes from

    Grossmann, I. E., & Floudas, C. A. (1987). Active constraint strategy for
    flexibility analysis in chemical processes. Computers & Chemical Engineering,
    11(6), 675-693.
    """
    m = pe.ConcreteModel()

    m.uncertain_temps_set = pe.Set(initialize=[1, 3, 5, 8])
    m.uncertain_temps = pe.Param(
        m.uncertain_temps_set, mutable=True, initialize={1: 620, 3: 388, 5: 583, 8: 313}
    )
    nominal_values = pe.ComponentMap()
    for p in m.uncertain_temps.values():
        nominal_values[p] = p.value

    param_bounds = pe.ComponentMap()
    for p in m.uncertain_temps.values():
        param_bounds[p] = (p.value - 10.0, p.value + 10.0)

    m.variable_temps_set = pe.Set(initialize=[2, 4, 6, 7])
    m.variable_temps = pe.Var(m.variable_temps_set, bounds=(0, 1000))
    m.qc = pe.Var()

    m.balances = pe.Constraint([1, 2, 3, 4])
    m.balances[1] = 1.5 * (m.uncertain_temps[1] - m.variable_temps[2]) == 2 * (
        m.variable_temps[4] - m.uncertain_temps[3]
    )
    m.balances[2] = m.uncertain_temps[5] - m.variable_temps[6] == 2 * (
        563 - m.variable_temps[4]
    )
    m.balances[3] = m.variable_temps[6] - m.variable_temps[7] == 3 * (
        393 - m.uncertain_temps[8]
    )
    m.balances[4] = m.qc == 1.5 * (m.variable_temps[2] - 350)

    m.temp_approaches = pe.Constraint([1, 2, 3, 4])
    m.temp_approaches[1] = m.variable_temps[2] >= m.uncertain_temps[3]
    m.temp_approaches[2] = m.variable_temps[6] >= m.variable_temps[4]
    m.temp_approaches[3] = m.variable_temps[7] >= m.uncertain_temps[8]
    m.temp_approaches[4] = 393 <= m.variable_temps[6]

    m.performance = pe.Constraint(expr=m.variable_temps[7] <= 323)

    return m, nominal_values, param_bounds


@pytest.mark.component
class TestFlexTest(unittest.TestCase):
    def test_poly(self):
        m, nominal_values, param_bounds = create_poly_model()
        var_bounds = pe.ComponentMap()
        var_bounds[m.z] = (-20, 20)
        build_active_constraint_flextest(
            m,
            uncertain_params=list(nominal_values.keys()),
            param_nominal_values=nominal_values,
            param_bounds=param_bounds,
            valid_var_bounds=var_bounds,
        )
        with scip_solver() as opt:
            res = opt.solve(m)
            pe.assert_optimal_termination(res)
            self.assertAlmostEqual(m.max_constraint_violation.value, 48.4649, 4)
            self.assertAlmostEqual(m.z.value, -2.6513, 4)
            ndx = _VarIndex(m.theta, None)
            self.assertAlmostEqual(m.unc_param_vars[ndx].value, 65, 5)

    def test_hx_network(self):
        expected = [4, 8.8]
        enforce_eq = [False, True]
        for exp, ee in zip(expected, enforce_eq):
            m, nominal_values, param_bounds = create_hx_network_model()
            var_bounds = pe.ComponentMap()
            for v in m.variable_temps.values():
                var_bounds[v] = (100, 1000)
            var_bounds[m.qc] = interval.mul(
                1.5, 1.5, *interval.sub(100, 1000, 350, 350)
            )
            config = ActiveConstraintConfig()
            config.enforce_equalities = ee
            build_active_constraint_flextest(
                m,
                uncertain_params=list(nominal_values.keys()),
                param_nominal_values=nominal_values,
                param_bounds=param_bounds,
                valid_var_bounds=var_bounds,
                config=config,
            )
            with scip_solver() as opt:
                res = opt.solve(m, tee=False)
                pe.assert_optimal_termination(res)
                self.assertAlmostEqual(m.max_constraint_violation.value, exp, 4)
