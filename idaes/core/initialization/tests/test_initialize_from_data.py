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
"""
Tests for initializing from data
"""
import os
import pytest
import types

from pyomo.environ import ConcreteModel, Constraint, Var

from idaes.core.initialization.initialize_from_data import (
    FromDataInitializer,
)
from idaes.core.initialization.initializer_base import InitializationStatus


__author__ = "Andrew Lee"


path = os.path.dirname(__file__)
fname = os.path.join(path, "init_example.json")


class TestInitFromDataSubMethods:
    @pytest.mark.unit
    def test_config(self):
        initializer = FromDataInitializer()

        assert hasattr(initializer, "config")

        assert "constraint_tolerance" in initializer.config

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
    def test_workflow_from_json(self, model):
        initializer = FromDataInitializer()

        status = initializer.initialize(model, json_file=fname)

        assert model.v1.value == 4
        assert model.v2.value == 4
        assert model.v3.value == 4
        assert model.v4.value == 4

        assert not model.v1.fixed

        assert status == InitializationStatus.Ok
