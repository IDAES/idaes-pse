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
from contextlib import contextmanager
import unittest
import pytest
import pyomo.environ as pe
from idaes.apps.flexibility_analysis import (
    _check_dependencies,
    _check_relu_dependencies,
)
import idaes.apps.flexibility_analysis as flex
from idaes.apps.flexibility_analysis.examples import (
    linear_hx_network,
    idaes_hx_network,
    nonlin_hx_network,
)
from idaes.core.util.testing import _enable_scip_solver_for_testing


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


@pytest.mark.component
class TestExamples(unittest.TestCase):
    def test_nonlin_hx_network(self):
        with scip_solver():
            res = nonlin_hx_network.main(
                method=flex.FlexTestMethod.active_constraint,
            )
            self.assertAlmostEqual(res, 0.1474609375, 3)

    def test_linear_hx_network(self):
        with scip_solver():
            res = linear_hx_network.main(
                flex_index=False,
                method=flex.FlexTestMethod.active_constraint,
                plot_history=False,
            )
            self.assertEqual(
                res.termination, flex.FlexTestTermination.found_infeasible_point
            )
            self.assertAlmostEqual(res.max_constraint_violation, 8.8, 5)

            res = linear_hx_network.main(
                flex_index=False,
                method=flex.FlexTestMethod.vertex_enumeration,
                plot_history=False,
            )
            self.assertEqual(
                res.termination, flex.FlexTestTermination.found_infeasible_point
            )
            self.assertAlmostEqual(res.max_constraint_violation, 8.8, 5)

            res = linear_hx_network.main(
                flex_index=False,
                method=flex.FlexTestMethod.linear_decision_rule,
                plot_history=False,
            )
            self.assertEqual(
                res.termination, flex.FlexTestTermination.found_infeasible_point
            )

            res = linear_hx_network.main(
                flex_index=True,
                method=flex.FlexTestMethod.active_constraint,
                plot_history=False,
            )
            self.assertAlmostEqual(res, 0.5, 5)

            res = linear_hx_network.main(
                flex_index=True,
                method=flex.FlexTestMethod.vertex_enumeration,
                plot_history=False,
            )
            self.assertAlmostEqual(res, 0.5, 5)

            res = linear_hx_network.main(
                flex_index=True,
                method=flex.FlexTestMethod.linear_decision_rule,
                plot_history=False,
            )
            self.assertLessEqual(res, 0.5)

            res = linear_hx_network.main(
                flex_index=True,
                method=flex.FlexTestMethod.relu_decision_rule,
                plot_history=False,
                relu_epochs=100,
            )
            self.assertGreaterEqual(res, 0.1)
            self.assertLessEqual(res, 0.5)

    def test_idaes_hx_network(self):
        with scip_solver():
            res = idaes_hx_network.main(flex.FlexTestMethod.sampling)
            self.assertEqual(
                res.termination, flex.FlexTestTermination.found_infeasible_point
            )
            self.assertAlmostEqual(
                res.max_constraint_violation, 0.375170890924453
            )  # regression
