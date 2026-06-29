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
Tests for Block Triangularization initialization
"""

import pytest
import re
import types

from pyomo.environ import ConcreteModel, Constraint, Var

from idaes.core.initialization.block_triangularization import (
    BlockTriangularizationInitializer,
)
from idaes.core.initialization.initializer_base import InitializationStatus
from idaes.core.util.exceptions import InitializationError

__author__ = "Andrew Lee"


class TestBTSubMethods:
    @pytest.mark.unit
    def test_config(self):
        initializer = BlockTriangularizationInitializer()

        assert hasattr(initializer, "config")

        assert "constraint_tolerance" in initializer.config

        assert "block_solver" in initializer.config
        assert "block_solver_options" in initializer.config
        assert "block_solver_call_options" in initializer.config
        assert "calculate_variable_options" in initializer.config

    # TODO: Tests for prechecks and initialization_routine stand alone

    @pytest.fixture
    def model(self):
        m = ConcreteModel()

        m.v1 = Var()
        m.v2 = Var()
        m.v3 = Var()
        m.v4 = Var()

        m.c1 = Constraint(expr=m.v1 == m.v2)
        m.c2 = Constraint(expr=2 * m.v2 == m.v3 + m.v4)
        m.c3 = Constraint(expr=m.v3 - m.v4 == 0)

        # Add a dummy method for fixing initialization states
        def fix_initialization_states(blk):
            blk.v1.fix(4)

        m.fix_initialization_states = types.MethodType(fix_initialization_states, m)

        return m

    @pytest.mark.component
    def test_workflow(self, model):
        initializer = BlockTriangularizationInitializer()

        status = initializer.initialize(model)

        assert model.v1.value == 4
        assert model.v2.value == 4
        assert model.v3.value == 4
        assert model.v4.value == 4

        assert not model.v1.fixed

        assert status == InitializationStatus.Ok

    @pytest.mark.component
    def test_workflow_no_presolve(self, model):
        # If the linear presolve is used, it will solve
        # the entire system without passing it to IPOPT.
        # Here we turn it off so IPOPT is called.
        initializer = BlockTriangularizationInitializer(
            block_solver_writer_config={"linear_presolve": False},
        )

        status = initializer.initialize(model)

        assert model.v1.value == 4
        assert model.v2.value == 4
        assert model.v3.value == 4
        assert model.v4.value == 4

        assert not model.v1.fixed

        assert status == InitializationStatus.Ok

    # Smoke test to make sure skip_final_solve option doesn't raise an error
    @pytest.mark.component
    def test_skip_final_solve(self, model):
        initializer = BlockTriangularizationInitializer(skip_final_solve=True)

        status = initializer.initialize(model)

        assert model.v1.value == 4
        assert model.v2.value == 4
        assert model.v3.value == 4
        assert model.v4.value == 4

        assert not model.v1.fixed

        assert status == InitializationStatus.Ok

    @pytest.mark.component
    def test_final_solve_fail(self, model):
        model.v5 = Var(bounds=(5, None))
        model.c4 = Constraint(expr=(model.v4 == model.v5))

        # In order to guarantee that block triangularization fails
        # during the final solve (and not in solve_strongly_connected_components),
        # we need to exploit some detailed properties in how it works.
        # 1) The 1x1 solver calculate_variable_from_constraints ignores variable bounds
        # 2) The linear presolve will recognize this system as being infeasible,
        #    so we have to skip it
        # 3) The final solver will run into the bound for v5 and return infeasible
        initializer = BlockTriangularizationInitializer(
            block_solver_options={"max_iter": 1},
            block_solver_writer_config={"linear_presolve": False},
        )

        with pytest.raises(
            InitializationError,
            match=re.escape(
                "Could not solve unknown after block triangularization finished."
            ),
        ):
            initializer.initialize(model)
