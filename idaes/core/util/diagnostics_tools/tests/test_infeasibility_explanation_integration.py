#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
This module contains tests for the DiagnosticsToolbox interface to Pyomo's infeasibility explanation tool.
"""

from io import StringIO

import pytest

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Var,
)

from idaes.core.util.diagnostics_tools.diagnostics_toolbox import (
    DiagnosticsToolbox,
)

__author__ = "Alex Dowling, Douglas Allan, Andrew Lee"


class TestComputeInfeasibilityExplanation:

    @pytest.fixture(scope="class")
    def model(self):
        # create an infeasible model for demonstration
        m = ConcreteModel()

        m.name = "test_infeas"
        m.x = Var([1, 2], bounds=(0, 1))
        m.y = Var(bounds=(0, 1))

        m.c = Constraint(expr=m.x[1] * m.x[2] == -1)
        m.d = Constraint(expr=m.x[1] + m.y >= 1)

        return m

    @pytest.mark.component
    @pytest.mark.solver
    def test_output(self, model):
        dt = DiagnosticsToolbox(model)

        stream = StringIO()

        dt.compute_infeasibility_explanation(stream=stream)

        expected = """Computed Minimal Intractable System (MIS)!
Constraints / bounds in MIS:
	lb of var x[2]
	lb of var x[1]
	constraint: c
Constraints / bounds in guards for stability:
"""
        assert expected in stream.getvalue()
