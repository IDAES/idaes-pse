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
This module contains model diagnostic utility functions for use in IDAES (Pyomo) models.
"""
from io import StringIO
import math
import numpy as np
import pytest

from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    Expression,
    log,
    tan,
    asin,
    acos,
    sqrt,
    Objective,
    Set,
    SolverFactory,
    Suffix,
    TransformationFactory,
    units,
    value,
    Var,
    Param,
    Integers,
)
from pyomo.common.collections import ComponentSet
from pyomo.contrib.pynumero.asl import AmplInterface
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP

import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.solvers import get_solver
from idaes.core import FlowsheetBlock
from idaes.core.util.testing import PhysicalParameterTestBlock
from unittest import TestCase

# TODO: Add pyomo.dae test case
"""
from pyomo.environ import TransformationFactory
from pyomo.dae import ContinuousSet, DerivativeVar
"""

# Need to update
from idaes.core.util.model_diagnostics import (
    DiagnosticsToolbox,
    SVDToolbox,
    DegeneracyHunter,
    DegeneracyHunter2,
    svd_dense,
    svd_sparse,
    get_valid_range_of_component,
    set_bounds_from_valid_range,
    list_components_with_values_outside_valid_range,
    ipopt_solve_halt_on_error,
    _var_in_block,
    _vars_fixed_to_zero,
    _vars_near_zero,
    _vars_violating_bounds,
    _vars_with_none_value,
    _vars_with_extreme_values,
    _write_report_section,
    _collect_model_statistics,
)
from idaes.core.util.testing import _enable_scip_solver_for_testing


__author__ = "Alex Dowling, Douglas Allan, Andrew Lee"


@pytest.fixture(scope="module")
def scip_solver():
    solver = SolverFactory("scip")
    undo_changes = None

    if not solver.available():
        undo_changes = _enable_scip_solver_for_testing()
    if not solver.available():
        pytest.skip(reason="SCIP solver not available")
    yield solver
    if undo_changes is not None:
        undo_changes()


@pytest.fixture
def model():
    m = ConcreteModel()
    m.b = Block()

    m.v1 = Var()
    m.v2 = Var()
    m.v3 = Var()
    m.v4 = Var()

    m.v1.fix(0)
    m.v2.fix(3)
    m.v3.set_value(0)

    return m


@pytest.mark.unit
def test_var_in_block(model):
    assert _var_in_block(model.v1, model)
    assert not _var_in_block(model.v1, model.b)


@pytest.mark.unit
def test_vars_fixed_to_zero(model):
    zero_vars = _vars_fixed_to_zero(model)
    assert isinstance(zero_vars, ComponentSet)
    assert len(zero_vars) == 1
    for i in zero_vars:
        assert i is model.v1


@pytest.mark.unit
def test_vars_near_zero(model):
    model.v3.set_value(1e-5)

    near_zero_vars = _vars_near_zero(model, variable_zero_value_tolerance=1e-5)
    assert isinstance(near_zero_vars, ComponentSet)
    assert len(near_zero_vars) == 2
    for i in near_zero_vars:
        assert i.local_name in ["v1", "v3"]

    near_zero_vars = _vars_near_zero(model, variable_zero_value_tolerance=1e-6)
    assert isinstance(near_zero_vars, ComponentSet)
    assert len(near_zero_vars) == 1
    for i in near_zero_vars:
        assert i is model.v1


@pytest.mark.unit
def test_vars_with_none_value(model):
    none_value = _vars_with_none_value(model)

    assert isinstance(none_value, ComponentSet)
    assert len(none_value) == 1
    for i in none_value:
        assert i is model.v4


@pytest.mark.unit
def test_vars_with_bounds_issues(model):
    model.v1.setlb(2)
    model.v1.setub(6)
    model.v2.setlb(0)
    model.v2.setub(10)
    model.v4.set_value(10)
    model.v4.setlb(0)
    model.v4.setub(1)

    bounds_issue = _vars_violating_bounds(model, tolerance=0)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 2
    for i in bounds_issue:
        assert i.local_name in ["v1", "v4"]

    m = ConcreteModel()
    m.v = Var(initialize=-1e-8, bounds=(0, 1))

    bounds_issue = _vars_violating_bounds(m, tolerance=1e-6)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0


@pytest.mark.unit
def test_vars_with_extreme_values():
    m = ConcreteModel()
    m.v1 = Var(initialize=1e-12)  # below zero
    m.v2 = Var(initialize=1e-8)  # small
    m.v3 = Var(initialize=1e-4)
    m.v4 = Var(initialize=1e0)
    m.v5 = Var(initialize=1e4)
    m.v6 = Var(initialize=1e8)
    m.v6 = Var(initialize=1e12)  # large

    xvars = _vars_with_extreme_values(m, large=1e9, small=1e-7, zero=1e-10)

    assert len(xvars) == 2
    for i in xvars:
        assert i.name in ["v2", "v6"]


@pytest.mark.unit
def test_write_report_section_all():
    stream = StringIO()

    _write_report_section(
        stream=stream,
        lines_list=["a", "b", "c"],
        title="foo",
        line_if_empty="bar",
        end_line="baz",
        header="-",
        footer="=",
    )

    expected = """------------------------------------------------------------------------------------
foo

    a
    b
    c

baz
====================================================================================
"""
    assert stream.getvalue() == expected


@pytest.mark.unit
def test_write_report_section_no_lines():
    stream = StringIO()

    _write_report_section(
        stream=stream,
        lines_list=[],
        title="foo",
        line_if_empty="bar",
        end_line="baz",
        header="-",
        footer="=",
    )

    expected = """------------------------------------------------------------------------------------
foo

    bar

baz
====================================================================================
"""
    assert stream.getvalue() == expected


@pytest.mark.unit
def test_write_report_section_lines_only():
    stream = StringIO()

    _write_report_section(
        stream=stream,
        lines_list=["a", "b", "c"],
    )

    expected = """------------------------------------------------------------------------------------
    a
    b
    c

"""
    assert stream.getvalue() == expected


@pytest.mark.unit
class TestStatsWriter:
    def test_blocks(self):
        m = ConcreteModel()
        m.b1 = Block()
        m.b1.b2 = Block()
        m.b3 = Block()
        m.b3.b4 = Block()
        m.b5 = Block()
        m.b5.b6 = Block()

        m.b3.b4.deactivate()
        m.b5.deactivate()

        stats = _collect_model_statistics(m)

        assert len(stats) == 9
        assert stats[0] == "    Activated Blocks: 4 (Deactivated: 3)"
        assert (
            stats[1] == "    Free Variables in Activated Constraints: 0 (External: 0)"
        )
        assert stats[2] == "        Free Variables with only lower bounds: 0"
        assert stats[3] == "        Free Variables with only upper bounds: 0"
        assert stats[4] == "        Free Variables with upper and lower bounds: 0"
        assert (
            stats[5] == "    Fixed Variables in Activated Constraints: 0 (External: 0)"
        )
        assert stats[6] == "    Activated Equality Constraints: 0 (Deactivated: 0)"
        assert stats[7] == "    Activated Inequality Constraints: 0 (Deactivated: 0)"
        assert stats[8] == "    Activated Objectives: 0 (Deactivated: 0)"

    def test_constraints(self):
        m = ConcreteModel()
        m.b1 = Block()
        m.b2 = Block()

        m.b1.v1 = Var()
        m.b1.v2 = Var()
        m.b1.v3 = Var()
        m.b1.v4 = Var()

        m.b1.c1 = Constraint(expr=m.b1.v1 == m.b1.v2)
        m.b1.c2 = Constraint(expr=m.b1.v1 >= m.b1.v2)
        m.b1.c3 = Constraint(expr=m.b1.v1 == m.b1.v2)
        m.b1.c4 = Constraint(expr=m.b1.v1 >= m.b1.v2)

        m.b2.c1 = Constraint(expr=m.b1.v1 == m.b1.v2)
        m.b2.c2 = Constraint(expr=m.b1.v1 >= m.b1.v2)

        m.b2.deactivate()
        m.b1.c3.deactivate()
        m.b1.c4.deactivate()

        stats = _collect_model_statistics(m)

        assert len(stats) == 9
        assert stats[0] == "    Activated Blocks: 2 (Deactivated: 1)"
        assert (
            stats[1] == "    Free Variables in Activated Constraints: 2 (External: 0)"
        )
        assert stats[2] == "        Free Variables with only lower bounds: 0"
        assert stats[3] == "        Free Variables with only upper bounds: 0"
        assert stats[4] == "        Free Variables with upper and lower bounds: 0"
        assert (
            stats[5] == "    Fixed Variables in Activated Constraints: 0 (External: 0)"
        )
        assert stats[6] == "    Activated Equality Constraints: 1 (Deactivated: 1)"
        assert stats[7] == "    Activated Inequality Constraints: 1 (Deactivated: 1)"
        assert stats[8] == "    Activated Objectives: 0 (Deactivated: 0)"

    def test_objectives(self):
        m = ConcreteModel()
        m.b1 = Block()
        m.b2 = Block()

        m.b1.o1 = Objective(expr=1)
        m.b1.o2 = Objective(expr=1)

        m.b2.o1 = Objective(expr=1)

        m.b2.deactivate()
        m.b1.o2.deactivate()

        stats = _collect_model_statistics(m)

        assert len(stats) == 9
        assert stats[0] == "    Activated Blocks: 2 (Deactivated: 1)"
        assert (
            stats[1] == "    Free Variables in Activated Constraints: 0 (External: 0)"
        )
        assert stats[2] == "        Free Variables with only lower bounds: 0"
        assert stats[3] == "        Free Variables with only upper bounds: 0"
        assert stats[4] == "        Free Variables with upper and lower bounds: 0"
        assert (
            stats[5] == "    Fixed Variables in Activated Constraints: 0 (External: 0)"
        )
        assert stats[6] == "    Activated Equality Constraints: 0 (Deactivated: 0)"
        assert stats[7] == "    Activated Inequality Constraints: 0 (Deactivated: 0)"
        assert stats[8] == "    Activated Objectives: 1 (Deactivated: 1)"

    def test_free_variables(self):
        m = ConcreteModel()
        m.b1 = Block()

        m.v1 = Var(bounds=(0, 1))

        m.b1.v2 = Var(bounds=(0, None))
        m.b1.v3 = Var(bounds=(None, 0))
        m.b1.v4 = Var(bounds=(0, 1))
        m.b1.v5 = Var()
        m.b1.v6 = Var(bounds=(0, 1))

        m.b1.c1 = Constraint(expr=0 == m.v1 + m.b1.v2 + m.b1.v3 + m.b1.v4 + m.b1.v5)

        stats = _collect_model_statistics(m.b1)

        assert len(stats) == 9
        assert stats[0] == "    Activated Blocks: 1 (Deactivated: 0)"
        assert (
            stats[1] == "    Free Variables in Activated Constraints: 5 (External: 1)"
        )
        assert stats[2] == "        Free Variables with only lower bounds: 1"
        assert stats[3] == "        Free Variables with only upper bounds: 1"
        assert stats[4] == "        Free Variables with upper and lower bounds: 2"
        assert (
            stats[5] == "    Fixed Variables in Activated Constraints: 0 (External: 0)"
        )
        assert stats[6] == "    Activated Equality Constraints: 1 (Deactivated: 0)"
        assert stats[7] == "    Activated Inequality Constraints: 0 (Deactivated: 0)"
        assert stats[8] == "    Activated Objectives: 0 (Deactivated: 0)"

    def test_fixed_variables(self):
        m = ConcreteModel()
        m.b1 = Block()

        m.v1 = Var(bounds=(0, 1))

        m.b1.v2 = Var(bounds=(0, None))
        m.b1.v3 = Var(bounds=(None, 0))
        m.b1.v4 = Var(bounds=(0, 1))
        m.b1.v5 = Var()
        m.b1.v6 = Var(bounds=(0, 1))

        m.b1.c1 = Constraint(expr=0 == m.v1 + m.b1.v2 + m.b1.v3 + m.b1.v4 + m.b1.v5)

        m.v1.fix(0.5)
        m.b1.v2.fix(-0.5)
        m.b1.v5.fix(-0.5)

        stats = _collect_model_statistics(m.b1)

        assert len(stats) == 9
        assert stats[0] == "    Activated Blocks: 1 (Deactivated: 0)"
        assert (
            stats[1] == "    Free Variables in Activated Constraints: 2 (External: 0)"
        )
        assert stats[2] == "        Free Variables with only lower bounds: 0"
        assert stats[3] == "        Free Variables with only upper bounds: 1"
        assert stats[4] == "        Free Variables with upper and lower bounds: 1"
        assert (
            stats[5] == "    Fixed Variables in Activated Constraints: 3 (External: 1)"
        )
        assert stats[6] == "    Activated Equality Constraints: 1 (Deactivated: 0)"
        assert stats[7] == "    Activated Inequality Constraints: 0 (Deactivated: 0)"
        assert stats[8] == "    Activated Objectives: 0 (Deactivated: 0)"


@pytest.mark.solver
class TestDiagnosticsToolbox:
    @pytest.mark.unit
    def test_invalid_model_type(self):
        with pytest.raises(
            TypeError,
            match="model argument must be an instance of a Pyomo BlockData object "
            "\(either a scalar Block or an element of an indexed Block\).",
        ):
            DiagnosticsToolbox(model="foo")

        # Check for indexed Blocks
        m = ConcreteModel()
        m.s = Set(initialize=[1, 2, 3])
        m.b = Block(m.s)

        with pytest.raises(
            TypeError,
            match="model argument must be an instance of a Pyomo BlockData object "
            "\(either a scalar Block or an element of an indexed Block\).",
        ):
            DiagnosticsToolbox(model=m.b)

    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.b = Block()

        m.v1 = Var(units=units.m)  # External variable
        m.b.v2 = Var(units=units.m)
        m.b.v3 = Var(bounds=(0, 5))
        m.b.v4 = Var()
        m.b.v5 = Var(bounds=(0, 5))
        m.b.v6 = Var()
        m.b.v7 = Var(
            units=units.m, bounds=(0, 1)
        )  # Poorly scaled variable with lower bound
        m.b.v8 = Var()  # unused variable

        m.b.c1 = Constraint(expr=m.v1 + m.b.v2 == 10)  # Unit consistency issue
        m.b.c2 = Constraint(expr=m.b.v3 == m.b.v4 + m.b.v5)
        m.b.c3 = Constraint(expr=2 * m.b.v3 == 3 * m.b.v4 + 4 * m.b.v5 + m.b.v6)
        m.b.c4 = Constraint(expr=m.b.v7 == 2e-8 * m.v1)  # Poorly scaled constraint

        m.b.v2.fix(5)
        m.b.v5.fix(2)
        m.b.v6.fix(0)

        solver = get_solver()
        solver.solve(m)

        return m

    @pytest.mark.component
    def test_display_external_variables(self, model):
        dt = DiagnosticsToolbox(model=model.b)

        stream = StringIO()
        dt.display_external_variables(stream)

        expected = """====================================================================================
The following external variable(s) appear in constraints within the model:

    v1

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_display_unused_variables(self, model):
        dt = DiagnosticsToolbox(model=model.b)

        stream = StringIO()
        dt.display_unused_variables(stream)

        expected = """====================================================================================
The following variable(s) do not appear in any activated constraints within the model:

    b.v8

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_display_variables_fixed_to_zero(self, model):
        dt = DiagnosticsToolbox(model=model.b)

        stream = StringIO()
        dt.display_variables_fixed_to_zero(stream)

        expected = """====================================================================================
The following variable(s) are fixed to zero:

    b.v6

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_display_variables_at_or_outside_bounds(self, model):
        dt = DiagnosticsToolbox(model=model.b)

        stream = StringIO()
        dt.display_variables_at_or_outside_bounds(stream)

        expected = """====================================================================================
The following variable(s) have values at or outside their bounds (tol=0.0E+00):

    b.v3 (free): value=0.0 bounds=(0, 5)

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_display_variables_with_none_value(self, model):
        dt = DiagnosticsToolbox(model=model.b)

        stream = StringIO()
        dt.display_variables_with_none_value(stream)

        expected = """====================================================================================
The following variable(s) have a value of None:

    b.v8

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_display_variables_with_value_near_zero(self, model):
        dt = DiagnosticsToolbox(model=model.b)

        stream = StringIO()
        dt.display_variables_with_value_near_zero(stream)

        expected = """====================================================================================
The following variable(s) have a value close to zero (tol=1.0E-08):

    b.v3: value=0.0
    b.v6: value=0

====================================================================================
"""
        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_display_variables_with_extreme_values(self, model):
        dt = DiagnosticsToolbox(model=model.b)

        stream = StringIO()
        dt.display_variables_with_extreme_values(stream)

        expected = """====================================================================================
The following variable(s) have extreme values (<1.0E-04 or > 1.0E+04):

    b.v7: 1.0000939326524314e-07

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_display_variables_near_bounds(self, model):
        dt = DiagnosticsToolbox(model=model.b)

        stream = StringIO()
        dt.display_variables_near_bounds(stream)

        expected = """====================================================================================
The following variable(s) have values close to their bounds (abs=1.0E-04, rel=1.0E-04):

    b.v3: value=0.0 bounds=(0, 5)
    b.v7: value=1.0000939326524314e-07 bounds=(0, 1)

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_display_components_with_inconsistent_units(self, model):
        dt = DiagnosticsToolbox(model=model.b)

        stream = StringIO()
        dt.display_components_with_inconsistent_units(stream)

        expected = """====================================================================================
The following component(s) have unit consistency issues:

    b.c1

For more details on unit inconsistencies, import the assert_units_consistent method
from pyomo.util.check_units
====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_display_constraints_with_large_residuals(self, model):
        dt = DiagnosticsToolbox(model=model.b)

        stream = StringIO()
        dt.display_constraints_with_large_residuals(stream)

        expected = """====================================================================================
The following constraint(s) have large residuals (>1.0E-05):

    b.c2: 6.66667E-01

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_get_dulmage_mendelsohn_partition(self, model):
        # Clone model so we can add some singularities
        m = model.clone()

        # Create structural singularities
        m.b.v2.unfix()
        m.b.v4.fix(2)

        # Add a second set of structural singularities
        m.b.b2 = Block()
        m.b.b2.v1 = Var()
        m.b.b2.v2 = Var()
        m.b.b2.v3 = Var()
        m.b.b2.v4 = Var()

        m.b.b2.c1 = Constraint(expr=m.b.b2.v1 == m.b.b2.v2)
        m.b.b2.c2 = Constraint(expr=2 * m.b.b2.v1 == 3 * m.b.b2.v2)
        m.b.b2.c3 = Constraint(expr=m.b.b2.v3 == m.b.b2.v4)

        m.b.b2.v2.fix(42)

        dt = DiagnosticsToolbox(model=m.b)

        (
            uc_vblocks,
            uc_cblocks,
            oc_vblocks,
            oc_cblocks,
        ) = dt.get_dulmage_mendelsohn_partition()

        assert len(uc_vblocks) == 2
        assert len(uc_vblocks[0]) == 3
        for i in uc_vblocks[0]:
            assert i.name in ["v1", "b.v2", "b.v7"]
        assert len(uc_vblocks[1]) == 2
        for i in uc_vblocks[1]:
            assert i.name in ["b.b2.v3", "b.b2.v4"]

        assert len(uc_cblocks) == 2
        assert len(uc_cblocks[0]) == 2
        for i in uc_cblocks[0]:
            assert i.name in ["b.c1", "b.c4"]
        assert len(uc_cblocks[1]) == 1
        for i in uc_cblocks[1]:
            assert i.name in ["b.b2.c3"]

        assert len(oc_vblocks) == 2
        assert len(oc_vblocks[0]) == 1
        for i in oc_vblocks[0]:
            assert i.name in ["b.v3"]
        assert len(oc_vblocks[1]) == 1
        for i in oc_vblocks[1]:
            assert i.name in ["b.b2.v1"]

        assert len(oc_cblocks) == 2
        assert len(oc_cblocks[0]) == 2
        for i in oc_cblocks[0]:
            assert i.name in ["b.c2", "b.c3"]
        assert len(oc_cblocks[1]) == 2
        for i in oc_cblocks[1]:
            assert i.name in ["b.b2.c1", "b.b2.c2"]

    @pytest.mark.component
    def test_display_underconstrained_set(self, model):
        # Clone model so we can add some singularities
        m = model.clone()

        # Create structural singularities
        m.b.v2.unfix()
        m.b.v4.fix(2)

        # Add a second set of structural singularities
        m.b.b2 = Block()
        m.b.b2.v1 = Var()
        m.b.b2.v2 = Var()
        m.b.b2.v3 = Var()
        m.b.b2.v4 = Var()

        m.b.b2.c1 = Constraint(expr=m.b.b2.v1 == m.b.b2.v2)
        m.b.b2.c2 = Constraint(expr=2 * m.b.b2.v1 == 3 * m.b.b2.v2)
        m.b.b2.c3 = Constraint(expr=m.b.b2.v3 == m.b.b2.v4)

        m.b.b2.v2.fix(42)

        dt = DiagnosticsToolbox(model=m.b)

        stream = StringIO()
        dt.display_underconstrained_set(stream)

        expected = """====================================================================================
Dulmage-Mendelsohn Under-Constrained Set

    Independent Block 0:

        Variables:

            b.v2
            v1
            b.v7

        Constraints:

            b.c1
            b.c4

    Independent Block 1:

        Variables:

            b.b2.v4
            b.b2.v3

        Constraints:

            b.b2.c3

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_display_overconstrained_set(self, model):
        # Clone model so we can add some singularities
        m = model.clone()

        # Create structural singularities
        m.b.v2.unfix()
        m.b.v4.fix(2)

        # Add a second set of structural singularities
        m.b.b2 = Block()
        m.b.b2.v1 = Var()
        m.b.b2.v2 = Var()
        m.b.b2.v3 = Var()
        m.b.b2.v4 = Var()

        m.b.b2.c1 = Constraint(expr=m.b.b2.v1 == m.b.b2.v2)
        m.b.b2.c2 = Constraint(expr=2 * m.b.b2.v1 == 3 * m.b.b2.v2)
        m.b.b2.c3 = Constraint(expr=m.b.b2.v3 == m.b.b2.v4)

        m.b.b2.v2.fix(42)

        dt = DiagnosticsToolbox(model=m.b)

        stream = StringIO()
        dt.display_overconstrained_set(stream)

        expected = """====================================================================================
Dulmage-Mendelsohn Over-Constrained Set

    Independent Block 0:

        Variables:

            b.v3

        Constraints:

            b.c2
            b.c3

    Independent Block 1:

        Variables:

            b.b2.v1

        Constraints:

            b.b2.c1
            b.b2.c2

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_display_variables_with_extreme_jacobians(self):
        model = ConcreteModel()
        model.v1 = Var(initialize=1e-8)
        model.v2 = Var()
        model.v3 = Var()

        model.c1 = Constraint(expr=model.v1 == model.v2)
        model.c2 = Constraint(expr=model.v1 == 1e-8 * model.v3)
        model.c3 = Constraint(expr=1e8 * model.v1 + 1e10 * model.v2 == 1e-6 * model.v3)

        dt = DiagnosticsToolbox(model=model)

        stream = StringIO()
        dt.display_variables_with_extreme_jacobians(stream)

        expected = """====================================================================================
The following variable(s) are associated with extreme Jacobian values (<1.0E-04 or>1.0E+04):

    v2: 1.000E+10
    v1: 1.000E+08
    v3: 1.000E-06

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_display_constraints_with_extreme_jacobians(self):
        model = ConcreteModel()
        model.v1 = Var(initialize=1e-8)
        model.v2 = Var()
        model.v3 = Var()

        model.c1 = Constraint(expr=model.v1 == model.v2)
        model.c2 = Constraint(expr=model.v1 == 1e-8 * model.v3)
        model.c3 = Constraint(expr=1e8 * model.v1 + 1e10 * model.v2 == 1e-6 * model.v3)

        dt = DiagnosticsToolbox(model=model)

        stream = StringIO()
        dt.display_constraints_with_extreme_jacobians(stream)

        expected = """====================================================================================
The following constraint(s) are associated with extreme Jacobian values (<1.0E-04 or>1.0E+04):

    c3: 1.000E+10

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_display_extreme_jacobian_entries(self):
        model = ConcreteModel()
        model.v1 = Var(initialize=1e-8)
        model.v2 = Var()
        model.v3 = Var()

        model.c1 = Constraint(expr=model.v1 == model.v2)
        model.c2 = Constraint(expr=model.v1 == 1e-8 * model.v3)
        model.c3 = Constraint(expr=1e8 * model.v1 + 1e10 * model.v2 == 1e-6 * model.v3)

        dt = DiagnosticsToolbox(model=model)

        stream = StringIO()
        dt.display_extreme_jacobian_entries(stream)

        expected = """====================================================================================
The following constraint(s) and variable(s) are associated with extreme Jacobian
values (<1.0E-04 or>1.0E+04):

    c3, v2: 1.000E+10
    c2, v3: 1.000E-08
    c3, v1: 1.000E+08
    c3, v3: 1.000E-06

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_collect_structural_warnings_base_case(self, model):
        dt = DiagnosticsToolbox(model=model.b)

        warnings, next_steps = dt._collect_structural_warnings()

        assert len(warnings) == 1
        assert warnings == ["WARNING: 1 Component with inconsistent units"]
        assert len(next_steps) == 1
        assert next_steps == ["display_components_with_inconsistent_units()"]

    @pytest.mark.component
    def test_collect_structural_warnings_underconstrained(self, model):
        # Clone model so we can add some singularities
        m = model.clone()

        # Create structural singularities
        m.b.v2.unfix()

        dt = DiagnosticsToolbox(model=m.b)

        warnings, next_steps = dt._collect_structural_warnings()

        assert len(warnings) == 3
        assert "WARNING: 1 Component with inconsistent units" in warnings
        assert "WARNING: 1 Degree of Freedom" in warnings
        assert (
            """WARNING: Structural singularity found
        Under-Constrained Set: 3 variables, 2 constraints
        Over-Constrained Set: 0 variables, 0 constraints"""
            in warnings
        )

        assert len(next_steps) == 2
        assert "display_components_with_inconsistent_units()" in next_steps
        assert "display_underconstrained_set()" in next_steps

    @pytest.mark.component
    def test_collect_structural_warnings_overconstrained(self, model):
        # Clone model so we can add some singularities
        m = model.clone()

        # Fix units
        m.b.del_component(m.b.c1)
        m.b.c1 = Constraint(expr=m.v1 + m.b.v2 == 10 * units.m)

        # Create structural singularities
        m.b.v4.fix(2)

        dt = DiagnosticsToolbox(model=m.b)

        warnings, next_steps = dt._collect_structural_warnings()

        assert len(warnings) == 2
        assert "WARNING: -1 Degree of Freedom" in warnings
        assert (
            """WARNING: Structural singularity found
        Under-Constrained Set: 0 variables, 0 constraints
        Over-Constrained Set: 1 variables, 2 constraints"""
            in warnings
        )

        assert len(next_steps) == 1
        assert "display_overconstrained_set()" in next_steps

    @pytest.mark.component
    def test_collect_structural_cautions(self, model):
        dt = DiagnosticsToolbox(model=model.b)

        cautions = dt._collect_structural_cautions()

        assert len(cautions) == 2
        assert "Caution: 1 variable fixed to 0" in cautions
        assert "Caution: 1 unused variable (0 fixed)" in cautions

    @pytest.mark.component
    def test_collect_numerical_warnings(self, model):
        dt = DiagnosticsToolbox(model=model.b)

        warnings, next_steps = dt._collect_numerical_warnings()

        assert len(warnings) == 2
        assert "WARNING: 1 Constraint with large residuals (>1.0E-05)" in warnings
        assert "WARNING: 1 Variable at or outside bounds (tol=0.0E+00)" in warnings

        assert len(next_steps) == 2
        assert "display_constraints_with_large_residuals()" in next_steps
        assert "display_variables_at_or_outside_bounds()" in next_steps

    @pytest.mark.component
    def test_collect_numerical_warnings_corrected(self, model):
        m = model.clone()

        # Fix numerical issues
        m.b.v3.setlb(-5)
        m.b.v5.setub(10)

        solver = get_solver()
        solver.solve(m)

        dt = DiagnosticsToolbox(model=m.b)

        warnings, next_steps = dt._collect_numerical_warnings()

        assert len(warnings) == 0

        assert len(next_steps) == 0

    @pytest.mark.component
    def test_collect_numerical_warnings_jacobian(self):
        model = ConcreteModel()
        model.v1 = Var(initialize=1e-8)
        model.v2 = Var(initialize=0)
        model.v3 = Var(initialize=0)

        model.c1 = Constraint(expr=model.v1 == model.v2)
        model.c2 = Constraint(expr=model.v1 == 1e-8 * model.v3)
        model.c3 = Constraint(expr=1e8 * model.v1 + 1e10 * model.v2 == 1e-6 * model.v3)

        dt = DiagnosticsToolbox(model=model)

        warnings, next_steps = dt._collect_numerical_warnings()

        assert len(warnings) == 3
        assert (
            "WARNING: 2 Variables with extreme Jacobian values (<1.0E-08 or >1.0E+08)"
            in warnings
        )
        assert (
            "WARNING: 1 Constraint with extreme Jacobian values (<1.0E-08 or >1.0E+08)"
            in warnings
        )
        assert "WARNING: 1 Constraint with large residuals (>1.0E-05)" in warnings

        assert len(next_steps) == 3
        assert "display_variables_with_extreme_jacobians()" in next_steps
        assert "display_constraints_with_extreme_jacobians()" in next_steps
        assert "display_constraints_with_large_residuals()" in next_steps

    @pytest.mark.component
    def test_collect_numerical_cautions(self, model):
        dt = DiagnosticsToolbox(model=model.b)

        cautions = dt._collect_numerical_cautions()
        assert len(cautions) == 5
        assert (
            "Caution: 2 Variables with value close to their bounds (abs=1.0E-04, rel=1.0E-04)"
            in cautions
        )
        assert "Caution: 2 Variables with value close to zero (tol=1.0E-08)" in cautions
        assert "Caution: 1 Variable with None value" in cautions
        assert "Caution: 1 extreme Jacobian Entry (<1.0E-04 or >1.0E+04)" in cautions
        assert (
            "Caution: 1 Variable with extreme value (<1.0E-04 or >1.0E+04)" in cautions
        )

    @pytest.mark.component
    def test_collect_numerical_cautions_jacobian(self):
        model = ConcreteModel()
        model.v1 = Var(initialize=1e-8)
        model.v2 = Var(initialize=0)
        model.v3 = Var(initialize=0)

        model.c1 = Constraint(expr=model.v1 == model.v2)
        model.c2 = Constraint(expr=model.v1 == 1e-8 * model.v3)
        model.c3 = Constraint(expr=1e8 * model.v1 + 1e10 * model.v2 == 1e-6 * model.v3)

        dt = DiagnosticsToolbox(model=model)

        cautions = dt._collect_numerical_cautions()

        assert len(cautions) == 4
        assert "Caution: 3 Variables with value close to zero (tol=1.0E-08)" in cautions
        assert (
            "Caution: 3 Variables with extreme Jacobian values (<1.0E-04 or >1.0E+04)"
            in cautions
        )
        assert (
            "Caution: 1 Constraint with extreme Jacobian values (<1.0E-04 or >1.0E+04)"
            in cautions
        )
        assert "Caution: 4 extreme Jacobian Entries (<1.0E-04 or >1.0E+04)" in cautions

    @pytest.mark.component
    def test_assert_no_structural_warnings(self, model):
        m = model.clone()
        dt = DiagnosticsToolbox(model=m.b)

        with pytest.raises(AssertionError, match="Structural issues found \(1\)."):
            dt.assert_no_structural_warnings()

        # Fix units issue
        m.b.del_component(m.b.c1)
        m.b.c1 = Constraint(expr=m.v1 + m.b.v2 == 10 * units.m)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_assert_no_numerical_warnings(self, model):
        m = model.clone()
        dt = DiagnosticsToolbox(model=m.b)

        with pytest.raises(AssertionError, match="Numerical issues found \(2\)."):
            dt.assert_no_numerical_warnings()

        # Fix numerical issues
        m.b.v3.setlb(-5)

        solver = get_solver()
        solver.solve(m)

        dt = DiagnosticsToolbox(model=m.b)
        dt.assert_no_numerical_warnings()

    @pytest.mark.component
    def test_report_structural_issues(self, model):
        dt = DiagnosticsToolbox(model=model.b)

        stream = StringIO()
        dt.report_structural_issues(stream)

        expected = """====================================================================================
Model Statistics

        Activated Blocks: 1 (Deactivated: 0)
        Free Variables in Activated Constraints: 4 (External: 1)
            Free Variables with only lower bounds: 0
            Free Variables with only upper bounds: 0
            Free Variables with upper and lower bounds: 2
        Fixed Variables in Activated Constraints: 3 (External: 0)
        Activated Equality Constraints: 4 (Deactivated: 0)
        Activated Inequality Constraints: 0 (Deactivated: 0)
        Activated Objectives: 0 (Deactivated: 0)

------------------------------------------------------------------------------------
1 WARNINGS

    WARNING: 1 Component with inconsistent units

------------------------------------------------------------------------------------
2 Cautions

    Caution: 1 variable fixed to 0
    Caution: 1 unused variable (0 fixed)

------------------------------------------------------------------------------------
Suggested next steps:

    display_components_with_inconsistent_units()

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_report_numerical_issues(self, model):
        dt = DiagnosticsToolbox(model=model.b)

        stream = StringIO()
        dt.report_numerical_issues(stream)

        expected = """====================================================================================
Model Statistics

    Jacobian Condition Number: 1.700E+01

------------------------------------------------------------------------------------
2 WARNINGS

    WARNING: 1 Constraint with large residuals (>1.0E-05)
    WARNING: 1 Variable at or outside bounds (tol=0.0E+00)

------------------------------------------------------------------------------------
5 Cautions

    Caution: 2 Variables with value close to their bounds (abs=1.0E-04, rel=1.0E-04)
    Caution: 2 Variables with value close to zero (tol=1.0E-08)
    Caution: 1 Variable with extreme value (<1.0E-04 or >1.0E+04)
    Caution: 1 Variable with None value
    Caution: 1 extreme Jacobian Entry (<1.0E-04 or >1.0E+04)

------------------------------------------------------------------------------------
Suggested next steps:

    display_constraints_with_large_residuals()
    display_variables_at_or_outside_bounds()

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_report_numerical_issues_jacobian(self):
        model = ConcreteModel()
        model.v1 = Var(initialize=1e-8)
        model.v2 = Var(initialize=0)
        model.v3 = Var(initialize=0)

        model.c1 = Constraint(expr=model.v1 == model.v2)
        model.c2 = Constraint(expr=model.v1 == 1e-8 * model.v3)
        model.c3 = Constraint(expr=1e8 * model.v1 + 1e10 * model.v2 == 1e-6 * model.v3)

        dt = DiagnosticsToolbox(model=model)

        stream = StringIO()
        dt.report_numerical_issues(stream)

        expected = """====================================================================================
Model Statistics

    Jacobian Condition Number: 1.407E+18

------------------------------------------------------------------------------------
3 WARNINGS

    WARNING: 1 Constraint with large residuals (>1.0E-05)
    WARNING: 2 Variables with extreme Jacobian values (<1.0E-08 or >1.0E+08)
    WARNING: 1 Constraint with extreme Jacobian values (<1.0E-08 or >1.0E+08)

------------------------------------------------------------------------------------
4 Cautions

    Caution: 3 Variables with value close to zero (tol=1.0E-08)
    Caution: 3 Variables with extreme Jacobian values (<1.0E-04 or >1.0E+04)
    Caution: 1 Constraint with extreme Jacobian values (<1.0E-04 or >1.0E+04)
    Caution: 4 extreme Jacobian Entries (<1.0E-04 or >1.0E+04)

------------------------------------------------------------------------------------
Suggested next steps:

    display_constraints_with_large_residuals()
    display_variables_with_extreme_jacobians()
    display_constraints_with_extreme_jacobians()

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.skipif(
        not AmplInterface.available(), reason="pynumero_ASL is not available"
    )
    @pytest.mark.integration
    def test_prepare_svd_toolbox(self, model):
        dt = DiagnosticsToolbox(model=model.b)
        svd = dt.prepare_svd_toolbox()

        assert isinstance(svd, SVDToolbox)

    @pytest.mark.skipif(
        not AmplInterface.available(), reason="pynumero_ASL is not available"
    )
    @pytest.mark.integration
    def test_prepare_degeneracy_hunter(self, model):
        dt = DiagnosticsToolbox(model=model.b)
        dh = dt.prepare_degeneracy_hunter()

        assert isinstance(dh, DegeneracyHunter2)


def dummy_callback(arg1):
    pass


def dummy_callback2(arg1=None, arg2=None):
    pass


@pytest.mark.skipif(
    not AmplInterface.available(), reason="pynumero_ASL is not available"
)
class TestSVDToolbox:
    @pytest.mark.unit
    def test_svd_callback_domain(self, dummy_problem):
        with pytest.raises(
            ValueError,
            match="SVD callback must be a callable which takes at least two arguments.",
        ):
            SVDToolbox(dummy_problem, svd_callback="foo")

        with pytest.raises(
            ValueError,
            match="SVD callback must be a callable which takes at least two arguments.",
        ):
            SVDToolbox(dummy_problem, svd_callback=dummy_callback)

        svd = SVDToolbox(dummy_problem, svd_callback=dummy_callback2)
        assert svd.config.svd_callback is dummy_callback2

    @pytest.mark.unit
    def test_init(self, dummy_problem):
        svd = SVDToolbox(dummy_problem)

        assert svd._model is dummy_problem
        assert svd.u is None
        assert svd.s is None
        assert svd.v is None

        # Get Jacobian and NLP
        jac = {
            (0, 0): 100.0,
            (1, 1): 1.0,
            (2, 2): 10.0,
            (3, 3): 0.1,
            (4, 4): 5.0,
        }
        for i, j in jac.items():
            assert j == svd.jacobian[i]

        assert isinstance(svd.nlp, PyomoNLP)

    @pytest.mark.unit
    def test_init_small_model(self):
        m = ConcreteModel()
        m.v = Var()
        m.c = Constraint(expr=m.v == 10)

        with pytest.raises(
            ValueError,
            match="Model needs at least 2 equality constraints to perform svd_analysis.",
        ):
            svd = SVDToolbox(m)

    @pytest.mark.unit
    def test_run_svd_analysis(self, dummy_problem):
        svd = SVDToolbox(dummy_problem)

        assert svd.config.svd_callback is svd_dense

        svd.run_svd_analysis()

        np.testing.assert_array_almost_equal(
            svd.u,
            np.array(
                [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 0, 1, 0]]
            ),
        )
        np.testing.assert_array_almost_equal(svd.s, np.array([0.1, 1, 5, 10]))
        np.testing.assert_array_almost_equal(
            svd.v,
            np.array(
                [[0, 0, 0, 1, 0], [0, 1, 0, 0, 0], [0, 0, 0, 0, 1], [0, 0, 1, 0, 0]]
            ).T,
        )

    @pytest.mark.unit
    def test_run_svd_analysis_sparse(self, dummy_problem):
        svd = SVDToolbox(dummy_problem, svd_callback=svd_sparse)
        svd.run_svd_analysis()

        # SVD sparse is not consistent with signs - manually iterate and check abs value
        for i in range(5):
            for j in range(4):
                if (i, j) in [(1, 1), (2, 3), (3, 0), (4, 2)]:
                    assert abs(svd.u[i, j]) == pytest.approx(1, abs=1e-6, rel=1e-6)
                else:
                    assert svd.u[i, j] == pytest.approx(0, abs=1e-6)

        np.testing.assert_array_almost_equal(svd.s, np.array([0.1, 1, 5, 10]))

        for i in range(5):
            for j in range(4):
                if (i, j) in [(1, 1), (2, 3), (3, 0), (4, 2)]:
                    assert abs(svd.v[i, j]) == pytest.approx(1, abs=1e-6, rel=1e-6)
                else:
                    assert svd.v[i, j] == pytest.approx(0, abs=1e-6)

    @pytest.mark.unit
    def test_run_svd_analysis_sparse_limit(self, dummy_problem):
        svd = SVDToolbox(
            dummy_problem, svd_callback=svd_sparse, number_of_smallest_singular_values=2
        )
        svd.run_svd_analysis()

        # SVD sparse is not consistent with signs - manually iterate and check abs value
        for i in range(5):
            for j in range(2):
                if (i, j) in [(1, 1), (3, 0)]:
                    assert abs(svd.u[i, j]) == pytest.approx(1, abs=1e-6, rel=1e-6)
                else:
                    assert svd.u[i, j] == pytest.approx(0, abs=1e-6)

        np.testing.assert_array_almost_equal(svd.s, np.array([0.1, 1]))

        for i in range(5):
            for j in range(2):
                if (i, j) in [(1, 1), (3, 0)]:
                    assert abs(svd.v[i, j]) == pytest.approx(1, abs=1e-6, rel=1e-6)
                else:
                    assert svd.v[i, j] == pytest.approx(0, abs=1e-6)

    @pytest.mark.unit
    def test_display_rank_of_equality_constraints(self, dummy_problem):
        svd = SVDToolbox(dummy_problem)

        stream = StringIO()
        svd.display_rank_of_equality_constraints(stream=stream)

        expected = """====================================================================================

Number of Singular Values less than 1.0E-6 is 0

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.unit
    def test_display_rank_of_equality_constraints(self, dummy_problem):
        svd = SVDToolbox(dummy_problem, singular_value_tolerance=1)

        stream = StringIO()
        svd.display_rank_of_equality_constraints(stream=stream)

        expected = """====================================================================================

Number of Singular Values less than 1.0E+00 is 1

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.unit
    def test_display_underdetermined_variables_and_constraints(self, dummy_problem):
        svd = SVDToolbox(dummy_problem)

        stream = StringIO()
        svd.display_underdetermined_variables_and_constraints(stream=stream)

        expected = """====================================================================================
Constraints and Variables associated with smallest singular values

    Smallest Singular Value 1:

        Variables:

            x[3]

        Constraints:

            dummy_eqn[3]

    Smallest Singular Value 2:

        Variables:

            x[1]

        Constraints:

            dummy_eqn[1]

    Smallest Singular Value 3:

        Variables:

            x[4]

        Constraints:

            dummy_eqn[4]

    Smallest Singular Value 4:

        Variables:

            x[2]

        Constraints:

            dummy_eqn[2]

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.unit
    def test_display_underdetermined_variables_and_constraints_specific(
        self, dummy_problem
    ):
        svd = SVDToolbox(dummy_problem)

        stream = StringIO()
        svd.display_underdetermined_variables_and_constraints(
            singular_values=[1], stream=stream
        )

        expected = """====================================================================================
Constraints and Variables associated with smallest singular values

    Smallest Singular Value 1:

        Variables:

            x[3]

        Constraints:

            dummy_eqn[3]

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.unit
    def test_display_underdetermined_variables_and_constraints(self, dummy_problem):
        svd = SVDToolbox(dummy_problem, size_cutoff_in_singular_vector=1)

        stream = StringIO()
        svd.display_underdetermined_variables_and_constraints(stream=stream)

        expected = """====================================================================================
Constraints and Variables associated with smallest singular values

    Smallest Singular Value 1:

        Variables:


        Constraints:


    Smallest Singular Value 2:

        Variables:


        Constraints:


    Smallest Singular Value 3:

        Variables:


        Constraints:


    Smallest Singular Value 4:

        Variables:


        Constraints:


====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.unit
    def test_display_constraints_including_variable(self):
        m = ConcreteModel()
        m.s = Set(initialize=[1, 2, 3, 4])
        m.v = Var(m.s)

        m.c1 = Constraint(expr=m.v[1] + 2 * m.v[2] == 10)
        m.c2 = Constraint(expr=3 * m.v[2] + 4 * m.v[3] == 20)
        m.c3 = Constraint(expr=5 * m.v[3] + 6 * m.v[4] == 30)
        m.c4 = Constraint(expr=7 * m.v[4] + 8 * m.v[1] == 40)

        svd = SVDToolbox(m)

        stream = StringIO()
        svd.display_constraints_including_variable(variable=m.v[1], stream=stream)

        expected = """====================================================================================
The following constraints involve v[1]:

    c1: 1.000e+00
    c4: 8.000e+00

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.unit
    def test_display_constraints_including_variable_invalid(self):
        m = ConcreteModel()
        m.s = Set(initialize=[1, 2, 3, 4])
        m.v = Var(m.s)

        m.c1 = Constraint(expr=m.v[1] + 2 * m.v[2] == 10)
        m.c2 = Constraint(expr=3 * m.v[2] + 4 * m.v[3] == 20)
        m.c3 = Constraint(expr=5 * m.v[3] + 6 * m.v[4] == 30)
        m.c4 = Constraint(expr=7 * m.v[4] + 8 * m.v[1] == 40)

        svd = SVDToolbox(m)

        with pytest.raises(
            TypeError,
            match="variable argument must be an instance of a Pyomo _VarData "
            "object \(got foo\).",
        ):
            svd.display_constraints_including_variable(variable="foo")

    @pytest.mark.unit
    def test_display_constraints_including_variable_not_in_model(self):
        m = ConcreteModel()
        m.s = Set(initialize=[1, 2, 3, 4])
        m.v = Var(m.s)
        m2 = ConcreteModel()
        m2.y = Var()

        m.c1 = Constraint(expr=m.v[1] + 2 * m.v[2] == 10)
        m.c2 = Constraint(expr=3 * m.v[2] + 4 * m.v[3] == 20)
        m.c3 = Constraint(expr=5 * m.v[3] + 6 * m.v[4] == 30)
        m.c4 = Constraint(expr=7 * m.v[4] + 8 * m.v[1] == 40)

        svd = SVDToolbox(m)

        with pytest.raises(AttributeError, match="Could not find y in model."):
            svd.display_constraints_including_variable(variable=m2.y)

    @pytest.mark.unit
    def test_display_variables_in_constraint(self):
        m = ConcreteModel()
        m.s = Set(initialize=[1, 2, 3, 4])
        m.v = Var(m.s)

        m.c1 = Constraint(expr=m.v[1] + 2 * m.v[2] == 10)
        m.c2 = Constraint(expr=3 * m.v[2] + 4 * m.v[3] == 20)
        m.c3 = Constraint(expr=5 * m.v[3] + 6 * m.v[4] == 30)
        m.c4 = Constraint(expr=7 * m.v[4] + 8 * m.v[1] == 40)

        svd = SVDToolbox(m)

        stream = StringIO()
        svd.display_variables_in_constraint(constraint=m.c1, stream=stream)

        expected = """====================================================================================
The following variables are involved in c1:

    v[1]: 1.000e+00
    v[2]: 2.000e+00

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.unit
    def test_display_variables_in_constraint_invalid(self):
        m = ConcreteModel()
        m.s = Set(initialize=[1, 2, 3, 4])
        m.v = Var(m.s)

        m.c1 = Constraint(expr=m.v[1] + 2 * m.v[2] == 10)
        m.c2 = Constraint(expr=3 * m.v[2] + 4 * m.v[3] == 20)
        m.c3 = Constraint(expr=5 * m.v[3] + 6 * m.v[4] == 30)
        m.c4 = Constraint(expr=7 * m.v[4] + 8 * m.v[1] == 40)

        svd = SVDToolbox(m)

        with pytest.raises(
            TypeError,
            match="constraint argument must be an instance of a Pyomo _ConstraintData "
            "object \(got foo\).",
        ):
            svd.display_variables_in_constraint(constraint="foo")

    @pytest.mark.unit
    def test_display_variables_in_constraint_no_in_model(self):
        m = ConcreteModel()
        m.s = Set(initialize=[1, 2, 3, 4])
        m.v = Var(m.s)

        m.c1 = Constraint(expr=m.v[1] + 2 * m.v[2] == 10)
        m.c2 = Constraint(expr=3 * m.v[2] + 4 * m.v[3] == 20)
        m.c3 = Constraint(expr=5 * m.v[3] + 6 * m.v[4] == 30)
        m.c4 = Constraint(expr=7 * m.v[4] + 8 * m.v[1] == 40)

        c6 = Constraint(expr=m.v[1] == m.v[2])

        svd = SVDToolbox(m)

        with pytest.raises(
            AttributeError, match="Could not find AbstractScalarConstraint in model."
        ):
            svd.display_variables_in_constraint(constraint=c6)


@pytest.mark.skipif(
    not AmplInterface.available(), reason="pynumero_ASL is not available"
)
class TestDegeneracyHunter:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()

        m.I = Set(initialize=[i for i in range(1, 4)])

        m.x = Var(m.I, bounds=(0, 5), initialize=1.0)

        m.con1 = Constraint(expr=m.x[1] + m.x[2] >= 1)
        m.con2 = Constraint(expr=m.x[1] + m.x[2] + m.x[3] == 1)
        m.con3 = Constraint(expr=m.x[2] - 2 * m.x[3] <= 1)
        m.con4 = Constraint(expr=m.x[1] + m.x[3] >= 1)

        m.con5 = Constraint(expr=m.x[1] + m.x[2] + m.x[3] == 1)

        m.obj = Objective(expr=sum(m.x[i] for i in m.I))

        return m

    @pytest.mark.unit
    def test_init(self, model):
        dh = DegeneracyHunter2(model)

        assert dh._model is model

        # Get Jacobian and NLP
        jac = {
            (0, 0): 1.0,
            (0, 1): 1.0,
            (0, 2): 1.0,
            (1, 0): 1.0,
            (1, 1): 1.0,
            (1, 2): 1.0,
        }

        for i, j in jac.items():
            assert j == dh.jacobian[i]

        assert isinstance(dh.nlp, PyomoNLP)

        assert dh.degenerate_set == {}
        assert dh.irreducible_degenerate_sets == []

    @pytest.mark.unit
    def test_get_solver(self, model):
        dh = DegeneracyHunter2(model, solver="ipopt", solver_options={"maxiter": 50})

        solver = dh._get_solver()

        assert solver.options == {"maxiter": 50}

    @pytest.mark.unit
    def test_prepare_candidates_milp(self, model):
        dh = DegeneracyHunter2(model)
        dh._prepare_candidates_milp()

        assert isinstance(dh.candidates_milp, ConcreteModel)

    @pytest.mark.unit
    def test_identify_candidates(self, model):
        dh = DegeneracyHunter2(model)
        dh._prepare_candidates_milp()

        dh.candidates_milp.nu[0].set_value(-1e-05)
        dh.candidates_milp.nu[1].set_value(1e-05)

        dh.candidates_milp.y_pos[0].set_value(0)
        dh.candidates_milp.y_pos[1].set_value(1)

        dh.candidates_milp.y_neg[0].set_value(-0)
        dh.candidates_milp.y_neg[1].set_value(-0)

        dh.candidates_milp.abs_nu[0].set_value(1e-05)
        dh.candidates_milp.abs_nu[1].set_value(1e-05)

        dh._identify_candidates()

        assert dh.degenerate_set == {
            model.con2: -1e-05,
            model.con5: 1e-05,
        }

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_candidates_milp(self, model, scip_solver):
        dh = DegeneracyHunter2(model)
        dh._prepare_candidates_milp()
        dh._solve_candidates_milp()

        assert value(dh.candidates_milp.nu[0]) == pytest.approx(1e-05, rel=1e-5)
        assert value(dh.candidates_milp.nu[1]) == pytest.approx(-1e-05, rel=1e-5)

        assert value(dh.candidates_milp.y_pos[0]) == pytest.approx(0, abs=1e-5)
        assert value(dh.candidates_milp.y_pos[1]) == pytest.approx(0, rel=1e-5)

        assert value(dh.candidates_milp.y_neg[0]) == pytest.approx(0, abs=1e-5)
        assert value(dh.candidates_milp.y_neg[1]) == pytest.approx(1, abs=1e-5)

        assert value(dh.candidates_milp.abs_nu[0]) == pytest.approx(1e-05, rel=1e-5)
        assert value(dh.candidates_milp.abs_nu[1]) == pytest.approx(1e-05, rel=1e-5)

        assert dh.degenerate_set == {
            model.con2: value(dh.candidates_milp.nu[0]),
            model.con5: value(dh.candidates_milp.nu[1]),
        }

    @pytest.mark.unit
    def test_prepare_ids_milp(self, model):
        dh = DegeneracyHunter2(model)
        dh._prepare_ids_milp()

        assert isinstance(dh.ids_milp, ConcreteModel)

    @pytest.mark.unit
    def test_solve_ids_milp(self, model):
        dh = DegeneracyHunter2(model)
        dh._prepare_ids_milp()

        dh.ids_milp.nu[0].set_value(1)
        dh.ids_milp.nu[1].set_value(-1)

        dh.ids_milp.y[0].set_value(1)
        dh.ids_milp.y[1].set_value(1)

        ids_ = dh._get_ids()

        assert ids_ == {
            model.con2: 1,
            model.con5: -1,
        }

    @pytest.mark.solver
    @pytest.mark.component
    def test_solve_ids_milp(self, model, scip_solver):
        dh = DegeneracyHunter2(model)
        dh._prepare_ids_milp()
        ids_ = dh._solve_ids_milp(cons=model.con2)

        assert ids_ == {
            model.con2: 1,
            model.con5: -1,
        }

        assert value(dh.ids_milp.nu[0]) == pytest.approx(1, rel=1e-5)
        assert value(dh.ids_milp.nu[1]) == pytest.approx(-1, rel=1e-5)

        assert value(dh.ids_milp.y[0]) == pytest.approx(1, rel=1e-5)
        assert value(dh.ids_milp.y[1]) == pytest.approx(1, rel=1e-5)

    @pytest.mark.solver
    @pytest.mark.component
    def test_find_irreducible_degenerate_sets(self, model, scip_solver):
        dh = DegeneracyHunter2(model)
        dh.find_irreducible_degenerate_sets()

        assert dh.irreducible_degenerate_sets == [
            {model.con2: 1, model.con5: -1},
            {model.con5: 1, model.con2: -1},
        ]

    @pytest.mark.solver
    @pytest.mark.component
    def test_report_irreducible_degenerate_sets(self, model, scip_solver):
        stream = StringIO()

        dh = DegeneracyHunter2(model)
        dh.report_irreducible_degenerate_sets(stream=stream)

        expected = """====================================================================================
Irreducible Degenerate Sets

    Irreducible Degenerate Set 0
        nu    Constraint Name
        1.0   con2
        -1.0  con5

    Irreducible Degenerate Set 1
        nu    Constraint Name
        -1.0  con2
        1.0   con5

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.solver
    @pytest.mark.component
    def test_report_irreducible_degenerate_sets_none(self, model, scip_solver):
        stream = StringIO()

        # Delete degenerate constraint
        model.del_component(model.con5)

        dh = DegeneracyHunter2(model)
        dh.report_irreducible_degenerate_sets(stream=stream)

        expected = """====================================================================================
Irreducible Degenerate Sets

    No candidate equations. The Jacobian is likely full rank.

====================================================================================
"""

        assert stream.getvalue() == expected


@pytest.fixture()
def dummy_problem():
    m = ConcreteModel()

    m.I = Set(initialize=[i for i in range(5)])

    m.x = Var(m.I, initialize=1.0)

    diag = [100, 1, 10, 0.1, 5]
    out = [1, 1, 1, 1, 1]

    @m.Constraint(m.I)
    def dummy_eqn(b, i):
        return out[i] == diag[i] * m.x[i]

    m.obj = Objective(expr=0)
    return m


@pytest.fixture()
def u_exp():
    return np.array(
        [[0, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1], [1, 0, 0, 0], [0, 0, 1, 0]]
    )


@pytest.fixture()
def s_exp():
    return np.array([0.1, 1, 5, 10])


@pytest.fixture()
def v_exp():
    return np.array(
        [[0, 0, 0, 1, 0], [0, 1, 0, 0, 0], [0, 0, 0, 0, 1], [0, 0, 1, 0, 0]]
    ).T


@pytest.mark.unit
def test_deprecate_degeneracy_hunter(caplog):
    m = ConcreteModel()
    m.v = Var()
    m.o = Objective(expr=m.v)
    dh = DegeneracyHunter(m)

    msg = (
        "DEPRECATED: DegeneracyHunter is being deprecated in favor of the new "
        "DiagnosticsToolbox.  (deprecated in 2.2.0, will be removed in (or after) 3.0.0)"
    )
    assert msg.replace(" ", "") in caplog.records[0].message.replace("\n", "").replace(
        " ", ""
    )


@pytest.mark.skipif(
    not AmplInterface.available(), reason="pynumero_ASL is not available"
)
@pytest.mark.unit
def test_dense_svd(dummy_problem, u_exp, s_exp, v_exp):
    m = dummy_problem
    dh = DegeneracyHunter(m)
    dh.svd_analysis(dense=True)
    assert dh.s == pytest.approx(s_exp, 1e-6, abs=1e-10)
    assert u_exp == pytest.approx(np.abs(dh.u), 1e-6, abs=1e-10)
    assert v_exp == pytest.approx(np.abs(dh.v), 1e-6, abs=1e-10)


@pytest.mark.skipif(
    not AmplInterface.available(), reason="pynumero_ASL is not available"
)
@pytest.mark.unit
def test_sparse_svd(dummy_problem, u_exp, s_exp, v_exp):
    m = dummy_problem
    dh = DegeneracyHunter(m)
    dh.svd_analysis()
    assert dh.s == pytest.approx(s_exp, 1e-6, abs=1e-10)
    assert u_exp == pytest.approx(np.abs(dh.u), 1e-6, abs=1e-10)
    assert v_exp == pytest.approx(np.abs(dh.v), 1e-6, abs=1e-10)


@pytest.mark.skipif(
    not AmplInterface.available(), reason="pynumero_ASL is not available"
)
@pytest.mark.unit
def test_scaling(dummy_problem, u_exp, s_exp, v_exp):
    ssf = iscale.set_scaling_factor
    cst = iscale.constraint_scaling_transform
    m = dummy_problem
    cst(m.dummy_eqn[0], 1e-2)
    ssf(m.x[1], 1e3)
    cst(m.dummy_eqn[1], 1e3)
    cst(m.dummy_eqn[2], 1e-1)
    ssf(m.x[3], 1e-3)
    cst(m.dummy_eqn[3], 1e2)
    cst(m.dummy_eqn[4], 0.2)

    dh = DegeneracyHunter(m)
    dh.svd_analysis(dense=False)
    assert dh.s == pytest.approx(np.ones((4,)), 1e-6)
    dh.svd_analysis(dense=True)
    assert dh.s == pytest.approx(np.ones((4,)), 1e-6)


@pytest.mark.skipif(
    not AmplInterface.available(), reason="pynumero_ASL is not available"
)
@pytest.mark.unit
def test_underdetermined_variables_and_constraints(dummy_problem, capsys):
    m = dummy_problem
    dh = DegeneracyHunter(m)

    dh.svd_analysis(n_sv=3)
    captured = capsys.readouterr()
    assert captured.out == "Computing the 3 smallest singular value(s)\n"

    dh.underdetermined_variables_and_constraints()
    captured = capsys.readouterr()
    assert captured.out == (
        "Column:    Variable\n3: x[3]\n\nRow:    " "Constraint\n3: dummy_eqn[3]\n"
    )

    dh.underdetermined_variables_and_constraints(n_calc=3)
    captured = capsys.readouterr()
    assert captured.out == (
        "Column:    Variable\n4: x[4]\n\nRow:    " "Constraint\n4: dummy_eqn[4]\n"
    )
    with pytest.raises(
        ValueError,
        match="User wanted constraints and variables associated "
        "with the 4-th smallest singular value, "
        "but only 3 small singular values have been "
        "calculated. Run svd_analysis again and specify "
        "n_sv>=4.",
    ):
        dh.underdetermined_variables_and_constraints(n_calc=4)

    dh.underdetermined_variables_and_constraints(tol=2)
    captured = capsys.readouterr()
    assert captured.out == ("Column:    Variable\n\nRow:    Constraint\n")


@pytest.mark.skipif(
    not AmplInterface.available(), reason="pynumero_ASL is not available"
)
@pytest.mark.unit
def test_underdetermined_calls_svd_analysis(dummy_problem, capsys):
    m = dummy_problem
    dh = DegeneracyHunter(m)

    dh.underdetermined_variables_and_constraints(n_calc=1)
    captured = capsys.readouterr()
    assert captured.out == (
        "Computing the 4 smallest singular value(s)\n"
        "Column:    Variable\n3: x[3]\n\nRow:    "
        "Constraint\n3: dummy_eqn[3]\n"
    )


@pytest.mark.skipif(
    not AmplInterface.available(), reason="pynumero_ASL is not available"
)
@pytest.mark.unit
def test_sv_value_error(dummy_problem):
    m = dummy_problem
    dh = DegeneracyHunter(m)
    with pytest.raises(
        ValueError,
        match="For a 5 by 5 system, svd_analysis can compute at most 4 "
        "singular values and vectors, but 5 were called for.",
    ):
        dh.svd_analysis(n_sv=5)
    with pytest.raises(ValueError, match="Nonsense value for n_sv=-1 received."):
        dh.svd_analysis(n_sv=-1)


@pytest.mark.skipif(
    not AmplInterface.available(), reason="pynumero_ASL is not available"
)
@pytest.mark.unit
def test_single_eq_error(capsys):
    m = ConcreteModel()
    m.x = Var(initialize=1)
    m.con = Constraint(expr=(2 * m.x == 1))
    m.obj = Objective(expr=0)

    dh = DegeneracyHunter(m)
    with pytest.raises(
        ValueError,
        match="Model needs at least 2 equality constraints to " "perform svd_analysis.",
    ):
        dh.svd_analysis()

    dh.check_rank_equality_constraints()
    captured = capsys.readouterr()
    assert captured.out == (
        "\nChecking rank of Jacobian of equality constraints...\n"
        "Model contains 1 equality constraints and 1 variables.\n"
        "Only singular value: 2.0\n"
    )


# This was from
# @pytest.fixture()
def problem1():
    m = ConcreteModel()

    m.I = Set(initialize=[i for i in range(5)])

    m.x = Var(m.I, bounds=(-10, 10), initialize=1.0)

    m.con1 = Constraint(expr=m.x[0] + m.x[1] - m.x[3] >= 10)
    m.con2 = Constraint(expr=m.x[0] * m.x[3] + m.x[1] >= 0)
    m.con3 = Constraint(expr=m.x[4] * m.x[3] + m.x[0] * m.x[3] - m.x[4] == 0)

    m.obj = Objective(expr=sum(m.x[i] ** 2 for i in m.I))

    return m


def example2(with_degenerate_constraint=True):
    """Create the Pyomo model for Example 2

    Arguments:
        with_degenerate_constraint: Boolean, if True, include the redundant linear constraint

    Returns:
        m2: Pyomo model
    """

    m2 = ConcreteModel()

    m2.I = Set(initialize=[i for i in range(1, 4)])

    m2.x = Var(m2.I, bounds=(0, 5), initialize=1.0)

    m2.con1 = Constraint(expr=m2.x[1] + m2.x[2] >= 1)
    m2.con2 = Constraint(expr=m2.x[1] + m2.x[2] + m2.x[3] == 1)
    m2.con3 = Constraint(expr=m2.x[2] - 2 * m2.x[3] <= 1)
    m2.con4 = Constraint(expr=m2.x[1] + m2.x[3] >= 1)

    if with_degenerate_constraint:
        m2.con5 = Constraint(expr=m2.x[1] + m2.x[2] + m2.x[3] == 1)

    m2.obj = Objective(expr=sum(m2.x[i] for i in m2.I))

    return m2


def extract_constraint_names(cs):
    """Get constraint names from ComponentSet

    Arguments:
        cs: ComponentSet object

    Return:
        constraint_names: list of constraint names (strings)

    """

    constraint_names = []
    for i in cs:
        constraint_names.append(i.name)
    return constraint_names


# Problem 1
@pytest.mark.skipif(
    not AmplInterface.available(), reason="pynumero_ASL is not available"
)
@pytest.mark.skipif(not SolverFactory("ipopt").available(False), reason="no Ipopt")
@pytest.mark.unit
def test_problem1():
    # Create test problem
    m = problem1()

    # Specify Ipopt as the solver
    opt = SolverFactory("ipopt")

    # Specifying an iteration limit of 0 allows us to inspect the initial point
    opt.options["max_iter"] = 0

    # "Solving" the model with an iteration limit of 0 load the initial point and applies
    # any preprocessors (e.g., enforces bounds)
    opt.solve(m, tee=True)

    # Create Degeneracy Hunter object
    dh = DegeneracyHunter(m)

    # Find constraints with residuals > 0.1
    initial_point_constraints = dh.check_residuals(tol=0.1)

    # Check there are 2 constraints with large residuals
    assert len(initial_point_constraints) == 2

    initial_point_constraint_names = extract_constraint_names(initial_point_constraints)

    # Check first constraint
    assert initial_point_constraint_names[0] == "con1"

    # Check second constraint
    assert initial_point_constraint_names[1] == "con3"

    opt.options["max_iter"] = 50

    # Solve
    opt.solve(m, tee=True)

    # Find constraints with residuals > 0.1
    solution_constraints = dh.check_residuals(tol=1e-6)

    # Check at the solution no constraints are violated
    assert len(solution_constraints) == 0

    # Check no constraints are near their bounds
    solution_bounds = dh.check_variable_bounds(tol=0.1)

    # Check at the solution no constraints are violated
    assert len(solution_bounds) == 0


# Problem 2 without degenerate constraint
@pytest.mark.skipif(
    not AmplInterface.available(), reason="pynumero_ASL is not available"
)
@pytest.mark.skipif(not SolverFactory("ipopt").available(False), reason="no Ipopt")
@pytest.mark.unit
def test_problem2_without_degenerate_constraint():
    # Create test problem instance
    m2 = example2(with_degenerate_constraint=False)

    # Specify Ipopt as the solver
    opt = SolverFactory("ipopt")

    # Specifying an iteration limit of 0 allows us to inspect the initial point
    opt.options["max_iter"] = 0

    # "Solving" the model with an iteration limit of 0 load the initial point and applies
    # any preprocessors (e.g., enforces bounds)
    opt.solve(m2, tee=True)

    # Create Degeneracy Hunter object
    dh2 = DegeneracyHunter(m2)

    # Check for violated constraints at the initial point
    initial_point_constraints = dh2.check_residuals(tol=0.1)

    # Check there are 1 constraints with large residuals
    assert len(initial_point_constraints) == 1

    initial_point_constraint_names = extract_constraint_names(initial_point_constraints)

    # Check first constraint
    assert initial_point_constraint_names[0] == "con2"

    # Resolve
    opt.options["max_iter"] = 500
    opt.solve(m2, tee=True)

    # Check solution
    x_sln = []

    for i in m2.I:
        x_sln.append(m2.x[i]())

    assert pytest.approx(x_sln[0], abs=1e-6) == 1.0
    assert pytest.approx(x_sln[1], abs=1e-6) == 0.0
    assert pytest.approx(x_sln[2], abs=1e-6) == 0.0


# Problem 2 with degenerate constraint
@pytest.mark.skipif(
    not AmplInterface.available(), reason="pynumero_ASL is not available"
)
@pytest.mark.skipif(not SolverFactory("ipopt").available(False), reason="no Ipopt")
@pytest.mark.unit
def test_problem2_with_degenerate_constraint():
    # Create test problem instance
    m2 = example2(with_degenerate_constraint=True)

    # Specify Ipopt as the solver
    opt = SolverFactory("ipopt")

    # Specifying an iteration limit of 0 allows us to inspect the initial point
    opt.options["max_iter"] = 0

    # "Solving" the model with an iteration limit of 0 load the initial point and applies
    # any preprocessors (e.g., enforces bounds)
    opt.solve(m2, tee=True)

    # Create Degeneracy Hunter object
    dh2 = DegeneracyHunter(m2)

    # Check for violated constraints at the initial point
    initial_point_constraints = dh2.check_residuals(tol=0.1)

    # Check there are 2 constraints with large residuals
    assert len(initial_point_constraints) == 2

    initial_point_constraint_names = extract_constraint_names(initial_point_constraints)

    # Check first constraint
    assert initial_point_constraint_names[0] == "con2"

    # Check first constraint
    assert initial_point_constraint_names[1] == "con5"

    # Resolve
    opt.options["max_iter"] = 500
    opt.solve(m2, tee=True)

    # Check solution
    x_sln = []

    for i in m2.I:
        x_sln.append(m2.x[i]())

    assert pytest.approx(x_sln[0], abs=1e-6) == 1.0
    assert pytest.approx(x_sln[1], abs=1e-6) == 0.0
    assert pytest.approx(x_sln[2], abs=1e-6) == 0.0

    # Check the rank
    n_rank_deficient = dh2.check_rank_equality_constraints()

    # Test DH with SVD

    assert n_rank_deficient == 1

    # TODO: Add MILP solver to idaes get-extensions and add more tests


@pytest.mark.unit
def test_get_valid_range_of_component():
    m = ConcreteModel()
    m.fs = FlowsheetBlock()

    m.fs.params = PhysicalParameterTestBlock()
    m.fs.state = m.fs.params.build_state_block(m.fs.time)

    # No valid range set yet, should return None
    assert get_valid_range_of_component(m.fs.state[0].flow_vol) is None

    # Set valid_range for flow_vol
    meta = m.fs.params.get_metadata().properties
    meta.flow_vol[None]._set_valid_range((0, 1))

    assert get_valid_range_of_component(m.fs.state[0].flow_vol) == (0, 1)


@pytest.mark.unit
def test_get_valid_range_of_component_no_metadata():
    m = ConcreteModel()
    m.fs = FlowsheetBlock()

    with pytest.raises(
        AttributeError, match="Could not find metadata for component fs"
    ):
        get_valid_range_of_component(m.fs)


@pytest.mark.unit
def test_get_valid_range_of_component_no_metadata_entry(caplog):
    m = ConcreteModel()
    m.fs = FlowsheetBlock()

    m.fs.params = PhysicalParameterTestBlock()
    m.fs.state = m.fs.params.build_state_block(m.fs.time)

    caplog.set_level(idaeslog.DEBUG, logger="idaes.core.util")
    assert get_valid_range_of_component(m.fs.state[0].test_var) is None

    assert (
        "No metadata entry for component fs.state[0.0].test_var; returning None"
        in caplog.text
    )


@pytest.mark.unit
def test_set_bounds_from_valid_range_scalar():
    m = ConcreteModel()
    m.fs = FlowsheetBlock()

    m.fs.params = PhysicalParameterTestBlock()
    m.fs.state = m.fs.params.build_state_block(m.fs.time)

    # Set valid_range for flow_vol
    meta = m.fs.params.get_metadata().properties
    meta.flow_vol[None]._set_valid_range((0, 1))

    assert m.fs.state[0].flow_vol.bounds == (None, None)

    set_bounds_from_valid_range(m.fs.state[0].flow_vol)
    assert m.fs.state[0].flow_vol.bounds == (0, 1)

    meta.flow_vol[None]._set_valid_range(None)
    set_bounds_from_valid_range(m.fs.state[0].flow_vol)
    assert m.fs.state[0].flow_vol.bounds == (None, None)


@pytest.mark.unit
def test_set_bounds_from_valid_range_indexed():
    m = ConcreteModel()
    m.fs = FlowsheetBlock()

    m.fs.params = PhysicalParameterTestBlock()
    m.fs.state = m.fs.params.build_state_block(m.fs.time)

    # Set valid_range for flow_mol_phase_comp
    meta = m.fs.params.get_metadata().properties
    meta.flow_mol["phase_comp"]._set_valid_range((0, 1))

    for k in m.fs.state[0].flow_mol_phase_comp:
        assert m.fs.state[0].flow_mol_phase_comp[k].bounds == (None, None)

    set_bounds_from_valid_range(m.fs.state[0].flow_mol_phase_comp)
    for k in m.fs.state[0].flow_mol_phase_comp:
        assert m.fs.state[0].flow_mol_phase_comp[k].bounds == (0, 1)

    meta.flow_mol["phase_comp"]._set_valid_range(None)
    set_bounds_from_valid_range(m.fs.state[0].flow_mol_phase_comp)
    for k in m.fs.state[0].flow_mol_phase_comp:
        assert m.fs.state[0].flow_mol_phase_comp[k].bounds == (None, None)


@pytest.mark.unit
def test_set_bounds_from_valid_range_block():
    m = ConcreteModel()
    m.fs = FlowsheetBlock()

    m.fs.params = PhysicalParameterTestBlock()
    m.fs.state = m.fs.params.build_state_block(m.fs.time)

    # Set valid_range
    meta = m.fs.params.get_metadata().properties
    meta.flow_vol[None]._set_valid_range((2, 8))
    meta.flow_mol["phase_comp"]._set_valid_range((0, 1))

    set_bounds_from_valid_range(m.fs.state)
    assert m.fs.state[0].flow_vol.bounds == (2, 8)
    for k in m.fs.state[0].flow_mol_phase_comp:
        assert m.fs.state[0].flow_mol_phase_comp[k].bounds == (0, 1)


@pytest.mark.unit
def test_set_bounds_from_valid_range_invalid_type():
    m = ConcreteModel()
    m.fs = FlowsheetBlock()

    m.fs.foo = Expression(expr=1)

    with pytest.raises(
        TypeError,
        match="Component fs.foo does not have bounds. Only Vars and Params have bounds.",
    ):
        set_bounds_from_valid_range(m.fs.foo)


@pytest.mark.unit
def test_list_components_with_values_outside_valid_range_scalar():
    m = ConcreteModel()
    m.fs = FlowsheetBlock()

    m.fs.params = PhysicalParameterTestBlock()
    m.fs.state = m.fs.params.build_state_block(m.fs.time)

    # Check return if no valid range
    m.fs.state[0].flow_vol.set_value(5)
    clist = list_components_with_values_outside_valid_range(m.fs.state[0].flow_vol)
    assert clist == []

    # Set valid_range for flow_vol
    meta = m.fs.params.get_metadata().properties
    meta.flow_vol[None]._set_valid_range((0, 1))

    # Set value outside range (low)
    m.fs.state[0].flow_vol.set_value(-1)
    clist = list_components_with_values_outside_valid_range(m.fs.state[0].flow_vol)
    assert clist == [m.fs.state[0].flow_vol]

    # Set value outside range (high)
    m.fs.state[0].flow_vol.set_value(10)
    clist = list_components_with_values_outside_valid_range(m.fs.state[0].flow_vol)
    assert clist == [m.fs.state[0].flow_vol]

    # Set value at min range
    m.fs.state[0].flow_vol.set_value(0)
    clist = list_components_with_values_outside_valid_range(m.fs.state[0].flow_vol)
    assert clist == []

    # Set value at max range
    m.fs.state[0].flow_vol.set_value(1)
    clist = list_components_with_values_outside_valid_range(m.fs.state[0].flow_vol)
    assert clist == []


@pytest.mark.unit
def test_list_components_with_values_outside_valid_range_indexed():
    m = ConcreteModel()
    m.fs = FlowsheetBlock()

    m.fs.params = PhysicalParameterTestBlock()
    m.fs.state = m.fs.params.build_state_block(m.fs.time)

    # Set valid_range for flow_mol_phase_comp
    meta = m.fs.params.get_metadata().properties
    meta.flow_mol["phase_comp"]._set_valid_range((0, 1))

    # Set values for each index
    m.fs.state[0].flow_mol_phase_comp["p1", "c1"].set_value(-1)  # low
    m.fs.state[0].flow_mol_phase_comp["p1", "c2"].set_value(0)  # min range
    m.fs.state[0].flow_mol_phase_comp["p2", "c1"].set_value(1)  # max range
    m.fs.state[0].flow_mol_phase_comp["p2", "c2"].set_value(10)  # high

    clist = list_components_with_values_outside_valid_range(
        m.fs.state[0].flow_mol_phase_comp
    )
    assert len(clist) == 2
    for i in clist:
        assert i.name in [
            "fs.state[0.0].flow_mol_phase_comp[p1,c1]",
            "fs.state[0.0].flow_mol_phase_comp[p2,c2]",
        ]


@pytest.mark.unit
def test_list_components_with_values_outside_valid_range_block():
    m = ConcreteModel()
    m.fs = FlowsheetBlock()

    m.fs.params = PhysicalParameterTestBlock()
    m.fs.state = m.fs.params.build_state_block(m.fs.time)

    # Set valid_range for some vars
    meta = m.fs.params.get_metadata().properties
    meta.flow_vol[None]._set_valid_range((0, 1))
    meta.flow_mol["phase_comp"]._set_valid_range((0, 1))

    # Set values for same vars
    m.fs.state[0].flow_vol.set_value(100)  # high
    m.fs.state[0].flow_mol_phase_comp["p1", "c1"].set_value(-1)  # low
    m.fs.state[0].flow_mol_phase_comp["p1", "c2"].set_value(0)  # min range
    m.fs.state[0].flow_mol_phase_comp["p2", "c1"].set_value(1)  # max range
    m.fs.state[0].flow_mol_phase_comp["p2", "c2"].set_value(10)  # high

    clist = list_components_with_values_outside_valid_range(m.fs.state)
    assert len(clist) == 3
    for i in clist:
        assert i.name in [
            "fs.state[0.0].flow_vol",
            "fs.state[0.0].flow_mol_phase_comp[p1,c1]",
            "fs.state[0.0].flow_mol_phase_comp[p2,c2]",
        ]


@pytest.mark.component
def test_ipopt_solve_halt_on_error(capsys):
    m = ConcreteModel()

    m.v = Var(initialize=-5, bounds=(None, -1))
    m.e = Expression(expr=log(m.v))
    m.c = Constraint(expr=m.e == 1)

    try:
        results = ipopt_solve_halt_on_error(m)
    except Exception:  # we expect this to fail
        pass

    captured = capsys.readouterr()
    assert "c: can't evaluate log(-5)." in captured.out


class TestEvalErrorDetection(TestCase):
    @pytest.mark.unit
    def test_div(self):
        m = ConcreteModel()
        m.x = Var(bounds=(1, None))
        m.y = Var()
        m.c = Constraint(expr=m.y == 1 / m.x)
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(0)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w, "c: Potential division by 0 in 1/x; Denominator bounds are (0, inf)"
        )

        dtb.config.warn_for_evaluation_error_at_bounds = False
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(-1)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w, "c: Potential division by 0 in 1/x; Denominator bounds are (-1, inf)"
        )

    @pytest.mark.unit
    def test_pow1(self):
        m = ConcreteModel()
        m.x = Var(bounds=(None, None))
        m.y = Var()
        m.p = Param(initialize=2, mutable=True)
        m.c = Constraint(expr=m.y == m.x**m.p)
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.p.value = 2.5
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation error in x**p; base bounds are (-inf, inf); exponent bounds are (2.5, 2.5)",
        )

        m.x.setlb(1)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

    @pytest.mark.unit
    def test_pow2(self):
        m = ConcreteModel()
        m.x = Var(bounds=(1, None))
        m.y = Var()
        m.p = Var(domain=Integers)
        m.c = Constraint(expr=m.y == m.x**m.p)
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(None)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation error in x**p; base bounds are (-inf, inf); exponent bounds are (-inf, inf)",
        )

    @pytest.mark.unit
    def test_pow3(self):
        m = ConcreteModel()
        m.x = Var(bounds=(0, None))
        m.y = Var()
        m.p = Var(bounds=(0, None))
        m.c = Constraint(expr=m.y == m.x**m.p)
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

    @pytest.mark.unit
    def test_pow4(self):
        m = ConcreteModel()
        m.x = Var(bounds=(0, None))
        m.y = Var()
        m.c = Constraint(expr=m.y == m.x ** (-2))
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation error in x**-2; base bounds are (0, inf); exponent bounds are (-2, -2)",
        )

        dtb.config.warn_for_evaluation_error_at_bounds = False
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(-1)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation error in x**-2; base bounds are (-1, inf); exponent bounds are (-2, -2)",
        )

    @pytest.mark.unit
    def test_pow5(self):
        m = ConcreteModel()
        m.x = Var(bounds=(0, None))
        m.y = Var()
        m.c = Constraint(expr=m.y == m.x ** (-2.5))
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation error in x**-2.5; base bounds are (0, inf); exponent bounds are (-2.5, -2.5)",
        )

        dtb.config.warn_for_evaluation_error_at_bounds = False
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(-1)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation error in x**-2.5; base bounds are (-1, inf); exponent bounds are (-2.5, -2.5)",
        )

    @pytest.mark.unit
    def test_log(self):
        m = ConcreteModel()
        m.x = Var(bounds=(1, None))
        m.y = Var()
        m.c = Constraint(expr=m.y == log(m.x))
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(0)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential log of a non-positive number in log(x); Argument bounds are (0, inf)",
        )

        dtb.config.warn_for_evaluation_error_at_bounds = False
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(-1)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential log of a non-positive number in log(x); Argument bounds are (-1, inf)",
        )

    @pytest.mark.unit
    def test_tan(self):
        m = ConcreteModel()
        m.x = Var(bounds=(-math.pi / 4, math.pi / 4))
        m.y = Var()
        m.c = Constraint(expr=m.y == tan(m.x))
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(-math.pi)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: tan(x) may evaluate to -inf or inf; Argument bounds are (-3.141592653589793, 0.7853981633974483)",
        )

    @pytest.mark.unit
    def test_asin(self):
        m = ConcreteModel()
        m.x = Var(bounds=(-0.5, 0.5))
        m.y = Var()
        m.c = Constraint(expr=m.y == asin(m.x))
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(None)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation of asin outside [-1, 1] in asin(x); Argument bounds are (-inf, 0.5)",
        )

        m.x.setlb(-0.5)
        m.x.setub(None)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation of asin outside [-1, 1] in asin(x); Argument bounds are (-0.5, inf)",
        )

    @pytest.mark.unit
    def test_acos(self):
        m = ConcreteModel()
        m.x = Var(bounds=(-0.5, 0.5))
        m.y = Var()
        m.c = Constraint(expr=m.y == acos(m.x))
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(None)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation of acos outside [-1, 1] in acos(x); Argument bounds are (-inf, 0.5)",
        )

        m.x.setlb(-0.5)
        m.x.setub(None)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential evaluation of acos outside [-1, 1] in acos(x); Argument bounds are (-0.5, inf)",
        )

    @pytest.mark.unit
    def test_sqrt(self):
        m = ConcreteModel()
        m.x = Var(bounds=(1, None))
        m.y = Var()
        m.c = Constraint(expr=m.y == sqrt(m.x))
        dtb = DiagnosticsToolbox(m)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 0)

        m.x.setlb(-1)
        warnings = dtb._collect_potential_eval_errors()
        self.assertEqual(len(warnings), 1)
        w = warnings[0]
        self.assertEqual(
            w,
            "c: Potential square root of a negative number in sqrt(x); Argument bounds are (-1, inf)",
        )

    @pytest.mark.unit
    def test_display(self):
        stream = StringIO()
        m = ConcreteModel()
        m.x = Var()
        m.y = Var()
        m.obj = Objective(expr=m.x**2 + m.y**2.5)
        m.c1 = Constraint(expr=m.y >= log(m.x))
        m.c2 = Constraint(expr=m.y >= (m.x - 1) ** 2.5)
        m.c3 = Constraint(expr=m.x - 1 >= 0)
        dtb = DiagnosticsToolbox(m)
        dtb.display_potential_evaluation_errors(stream=stream)
        expected = "====================================================================================\n3 WARNINGS\n\n    c1: Potential log of a non-positive number in log(x); Argument bounds are (-inf, inf)\n    c2: Potential evaluation error in (x - 1)**2.5; base bounds are (-inf, inf); exponent bounds are (2.5, 2.5)\n    obj: Potential evaluation error in y**2.5; base bounds are (-inf, inf); exponent bounds are (2.5, 2.5)\n\n====================================================================================\n"
        got = stream.getvalue()
        exp_list = expected.split("\n")
        got_list = got.split("\n")
        self.assertEqual(len(exp_list), len(got_list))
        for _exp, _got in zip(exp_list, got_list):
            self.assertEqual(_exp, _got)
