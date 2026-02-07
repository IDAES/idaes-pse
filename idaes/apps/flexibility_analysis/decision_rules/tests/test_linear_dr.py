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
import unittest
import pyomo.environ as pe
from idaes.apps.flexibility_analysis.decision_rules.linear_dr import (
    construct_linear_decision_rule,
    LinearDRConfig,
)
import numpy as np
import pytest


def y1_func(x1, x2):
    return 3 * x1 - 2 * x2 + 5


def y2_func(x1, x2):
    return -x1 + 0.5 * x2


class TestLinearDecisionRule(unittest.TestCase):
    @pytest.mark.component
    def test_construct_linear_dr(self):
        x1_samples = [float(i) for i in np.linspace(-5, 5, 100)]
        x2_samples = [float(i) for i in np.linspace(-5, 5, 100)]
        x1_samples.extend(float(i) for i in np.linspace(-5, 5, 100))
        x2_samples.extend(float(i) for i in np.linspace(5, -5, 100))

        x1_samples = np.array(x1_samples)
        x2_samples = np.array(x2_samples)
        y1_samples = y1_func(x1_samples, x2_samples)
        y2_samples = y2_func(x1_samples, x2_samples)

        m = pe.ConcreteModel()
        m.x1 = pe.Var(initialize=1.7)
        m.x2 = pe.Var(initialize=-3.1)
        m.y1 = pe.Var(initialize=0.2)
        m.y2 = pe.Var(initialize=2.5)

        input_vals = pe.ComponentMap()
        input_vals[m.x1] = [float(i) for i in x1_samples]
        input_vals[m.x2] = [float(i) for i in x2_samples]

        output_vals = pe.ComponentMap()
        output_vals[m.y1] = [float(i) for i in y1_samples]
        output_vals[m.y2] = [float(i) for i in y2_samples]

        opt = pe.SolverFactory("ipopt")
        config = LinearDRConfig()
        config.solver = opt
        m.dr = construct_linear_decision_rule(
            input_vals=input_vals, output_vals=output_vals, config=config
        )

        self.assertEqual(pe.value(m.dr.decision_rule[0].lower), 0)
        self.assertEqual(pe.value(m.dr.decision_rule[0].upper), 0)
        self.assertAlmostEqual(
            pe.value(m.dr.decision_rule[0].body),
            y1_func(m.x1.value, m.x2.value) - m.y1.value,
        )

        self.assertEqual(pe.value(m.dr.decision_rule[1].lower), 0)
        self.assertEqual(pe.value(m.dr.decision_rule[1].upper), 0)
        self.assertAlmostEqual(
            pe.value(m.dr.decision_rule[1].body),
            y2_func(m.x1.value, m.x2.value) - m.y2.value,
        )
