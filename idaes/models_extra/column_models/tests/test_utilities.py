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
Tests for column model utitliy functions

Author: Alejandro Garcia-Diego, andrew Lee
"""
import pytest
from types import MethodType

from pyomo.environ import ConcreteModel, Block, Var, Reference, Set
from pyomo.network import Port

from idaes.models_extra.column_models.util import make_phase_split

# Mark module as all unit tests
pytestmark = pytest.mark.unit


@pytest.fixture
def m():
    # Mock up something that looks like a tray model
    m = ConcreteModel()
    m.time = Set(initialize=[1, 2, 3])

    def flowsheet(m):
        return m

    m.flowsheet = MethodType(flowsheet, m)

    m.properties_out = Block(m.time)

    m.port = Port()

    return m


def test_make_phase_split_general(m):
    # Test for variables without special indices
    # Mock define_port_members method
    def define_port_members(m):
        return {"temperature": m.temperature, "indexed": m.indexed}

    # Add variables and define_port_members method to dummy state block
    for i in m.time:
        m.properties_out[i].temperature = Var()
        m.properties_out[i].indexed = Var(["a", "b"])

        m.properties_out[i].define_port_members = MethodType(
            define_port_members, m.properties_out[i]
        )

    # Call make_phase_split function
    make_phase_split(m, m.port)

    # Check for expected results
    for i in m.time:
        assert m.port.temperature[i] is m.properties_out[i].temperature
        for j in ["a", "b"]:
            assert m.port.indexed[i, j] is m.properties_out[i].indexed[j]
