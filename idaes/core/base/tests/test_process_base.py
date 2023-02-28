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
Tests for process_base.

Author: Andrew Lee
"""
import pytest
from io import StringIO
import types
import pandas

from pyomo.environ import Block, ConcreteModel, Set, Var, Param, Expression, units

from idaes.core.base.process_base import ProcessBaseBlock
from idaes.core import FlowsheetBlockData, declare_process_block_class
from idaes.core.initialization import BlockTriangularizationInitializer


@declare_process_block_class("Flowsheet")
class _Flowsheet(FlowsheetBlockData):
    def build(self):
        super(FlowsheetBlockData, self).build()


@pytest.mark.unit
def test_base_attributes():
    # Test flowsheet method
    m = ConcreteModel()
    m.a = Flowsheet()

    assert m.a.initialization_order == [m.a]
    assert m.a.default_initializer is BlockTriangularizationInitializer


@pytest.mark.unit
def test_flowsheet():
    # Test flowsheet method
    m = ConcreteModel()
    m.a = Flowsheet()

    assert m.a.flowsheet() is None

    m.b = Block()
    m.b.c = Flowsheet()
    assert m.b.c.flowsheet() is None

    m.a.d = Flowsheet()
    assert m.a.d.flowsheet() is m.a

    m.a.e = Flowsheet([1, 2])
    assert m.a.e[1].flowsheet() is m.a

    m.a.e[1].f = Flowsheet()
    assert m.a.e[1].f.flowsheet() is m.a.e[1]

    m.a.g = Block()
    m.a.g.h = Flowsheet()
    assert m.a.g.h.flowsheet() is m.a

    m.a.i = Block([1, 2])
    m.a.i[1].j = Flowsheet()
    assert m.a.i[1].j.flowsheet() is m.a


@pytest.mark.unit
def test_get_performance_contents():
    m = ConcreteModel()
    m.b = ProcessBaseBlock()
    assert m.b._get_performance_contents(time_point=0) is None


@pytest.mark.unit
def test_get_stream_table_contents():
    m = ConcreteModel()
    m.b = ProcessBaseBlock()
    assert m.b._get_stream_table_contents(time_point=0) is None


@pytest.mark.unit
def test_report_dof():
    m = ConcreteModel()
    m.b = ProcessBaseBlock()

    stream = StringIO()

    m.b.report(ostream=stream, dof=True)

    expected = """
====================================================================================
Unit : b                                                                   Time: 0.0
====================================================================================
    Local Degrees of Freedom: 0
    Total Variables: 0    Activated Constraints: 0    Activated Blocks: 1
====================================================================================
"""

    assert stream.getvalue().strip() == expected.strip()


@pytest.mark.unit
def test_report_perf_dict():
    m = ConcreteModel()
    m.b = ProcessBaseBlock()

    m.b.s = Set(initialize=[1, 2])

    m.b.sv = Var(initialize=7, units=units.m)
    m.b.iv = Var(m.b.s, initialize=42, units=units.s)

    m.b.p = Param(initialize=4, units=units.mole)

    m.b.e = Expression(rule=72)

    def dummy_performance(self, time_point):
        return {
            "vars": {"Scalar Var": m.b.sv, "Indexed Var 1": m.b.iv[1]},
            "params": {"Param": m.b.p},
            "exprs": {"Expression": m.b.e},
        }

    m.b._get_performance_contents = types.MethodType(dummy_performance, m.b)

    stream = StringIO()

    m.b.report(ostream=stream)

    expected = """
====================================================================================
Unit : b                                                                   Time: 0.0
------------------------------------------------------------------------------------
    Unit Performance

    Variables: 

    Key           : Value  : Units  : Fixed : Bounds
    Indexed Var 1 : 42.000 : second : False : (None, None)
       Scalar Var : 7.0000 :  meter : False : (None, None)

    Expressions: 

    Key        : Value  : Units
    Expression : 72.000 : dimensionless

    Parameters: 

    Key   : Value : Units : Mutable
    Param :     4 :  mole :    True

====================================================================================
"""

    assert stream.getvalue().strip() == expected.strip()


@pytest.mark.unit
def test_report_stream_table():
    m = ConcreteModel()
    m.b = ProcessBaseBlock()

    data = {"x1": [1, 2, 3, 4], "x2": [5, 6, 7, 8], "z1": [10, 20, 30, 40]}
    df = pandas.DataFrame(data)

    def dummy_stream_table(self, time_point):
        return df

    m.b._get_stream_table_contents = types.MethodType(dummy_stream_table, m.b)

    stream = StringIO()

    m.b.report(ostream=stream)

    expected = """
====================================================================================
Unit : b                                                                   Time: 0.0
------------------------------------------------------------------------------------
    Stream Table
       x1  x2  z1
    0   1   5  10
    1   2   6  20
    2   3   7  30
    3   4   8  40
====================================================================================
"""

    assert stream.getvalue().strip() == expected.strip()


@pytest.mark.unit
def test_default_scaling():
    m = ConcreteModel()
    m.b = ProcessBaseBlock()

    assert isinstance(m.b._default_scaling_factors, dict)
    assert m.b.default_scaling_factors is m.b._default_scaling_factors
    assert m.b.default_scaling_factor is m.b._default_scaling_factors


@pytest.mark.unit
def test_set_default_scaling():
    m = ConcreteModel()
    m.b = ProcessBaseBlock()

    m.b.set_default_scaling("a", 11)
    m.b.set_default_scaling("b", 22, index="c")

    assert ("a", None) in m.b._default_scaling_factors
    assert m.b._default_scaling_factors[("a", None)] == 11
    assert ("b", "c") in m.b._default_scaling_factors
    assert m.b._default_scaling_factors[("b", "c")] == 22


@pytest.mark.unit
def test_unset_default_scaling():
    m = ConcreteModel()
    m.b = ProcessBaseBlock()

    m.b.set_default_scaling("a", 11)
    m.b.set_default_scaling("b", 22, index="c")

    m.b.unset_default_scaling("a")
    # Try something that does not exist
    m.b.unset_default_scaling("b", index="d")

    assert ("a", None) not in m.b._default_scaling_factors
    assert ("b", "c") in m.b._default_scaling_factors
    assert m.b._default_scaling_factors[("b", "c")] == 22


@pytest.mark.unit
def test_get_default_scaling():
    m = ConcreteModel()
    m.b = ProcessBaseBlock()

    m.b.set_default_scaling("a", 11)
    m.b.set_default_scaling("b", 22, index="c")

    assert m.b.get_default_scaling("a") == 11
    # Lookup scaling factor for an index that doesn't exist, but has overall scaling
    assert m.b.get_default_scaling("a", index="c") == 11
    assert m.b.get_default_scaling("b", "c") == 22
