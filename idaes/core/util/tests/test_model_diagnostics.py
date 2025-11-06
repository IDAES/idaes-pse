#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
This module contains model diagnostic utility functions for use in IDAES (Pyomo) models.
"""
from copy import deepcopy
from io import StringIO
import math
import os
from copy import deepcopy
from pyomo.contrib.pynumero.interfaces.external_grey_box import (
    ExternalGreyBoxBlock,
    ExternalGreyBoxModel,
)
import re
from unittest import TestCase

import numpy as np
from pandas import DataFrame
import pytest

from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    Expression,
    Expr_if,
    Objective,
    PositiveIntegers,
    Set,
    SolverFactory,
    units,
    value,
    Var,
    assert_optimal_termination,
    Param,
    Integers,
    exp,
    log,
    log10,
    sin,
    cos,
    tan,
    asin,
    acos,
    atan,
    sqrt,
    sinh,
    cosh,
    tanh,
    asinh,
    acosh,
    atanh,
)
from pyomo.common.collections import ComponentSet
from pyomo.contrib.pynumero.asl import AmplInterface
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
from pyomo.common.fileutils import this_file_dir
from pyomo.common.tempfiles import TempfileManager
from pyomo.core import expr as EXPR

import idaes.core.util.scaling as iscale
import idaes.logger as idaeslog
from idaes.core.solvers import get_solver
from idaes.core import FlowsheetBlock
from idaes.core.util.testing import PhysicalParameterTestBlock
from idaes.models.properties.modular_properties.eos.ceos_common import (
    cubic_roots_available,
    CubicThermoExpressions,
    CubicType as CubicEoS,
)
from idaes.core.util.model_diagnostics import (
    DiagnosticsToolbox,
    SVDToolbox,
    DegeneracyHunter,
    IpoptConvergenceAnalysis,
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
    check_parallel_jacobian,
    compute_ill_conditioning_certificate,
    ConstraintTermAnalysisVisitor,
)
from idaes.core.util.parameter_sweep import (
    SequentialSweepRunner,
    ParameterSweepSpecification,
)
from idaes.core.surrogate.pysmo.sampling import (
    UniformSampling,
)
from idaes.core.util.testing import _enable_scip_solver_for_testing
from idaes.models.properties import iapws95
from idaes.core.scaling import set_scaling_factor

__author__ = "Alex Dowling, Douglas Allan, Andrew Lee"


# TODO: Add pyomo.dae test cases
solver_available = SolverFactory("scip").available()
currdir = this_file_dir()


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

    set_scaling_factor(model.v3, 1e5)
    near_zero_vars = _vars_near_zero(model, variable_zero_value_tolerance=1e-5)
    assert isinstance(near_zero_vars, ComponentSet)
    assert len(near_zero_vars) == 1
    for i in near_zero_vars:
        assert i is model.v1

    near_zero_vars = _vars_near_zero(model, variable_zero_value_tolerance=1)
    assert isinstance(near_zero_vars, ComponentSet)
    assert len(near_zero_vars) == 2
    for i in near_zero_vars:
        assert i.local_name in ["v1", "v3"]


class TestVariablesWithNoneValue:
    def err_msg(self, j):
        return (
            f"Found {j} variables with a value of None in the mathematical "
            "program generated by the model. They must be initialized with non-None "
            "values before numerical analysis can proceed. Run "
            + DiagnosticsToolbox.display_variables_with_none_value_in_activated_constraints.__name__
            + " to display a list of these variables."
        )

    @pytest.fixture(scope="function")
    def model(self):
        m = ConcreteModel()
        m.u = Var()
        m.w = Var(initialize=1)
        m.x = Var(initialize=1)
        m.y = Var(range(3), initialize=[1, None, 2])
        m.z = Var()

        m.con_w = Constraint(expr=m.w == 0)
        m.con_x = Constraint(expr=m.x == 1)

        def rule_con_y(b, i):
            return b.y[i] == 3

        m.con_y = Constraint(range(3), rule=rule_con_y)
        m.con_z = Constraint(expr=m.x + m.z == -1)

        m.block1 = Block()
        m.block1.a = Var()
        m.block1.b = Var(initialize=7)
        m.block1.c = Var(initialize=-5)
        m.block1.d = Var()

        m.block1.con_a = Constraint(expr=m.block1.a**2 + m.block1.b == 1)
        m.block1.con_b = Constraint(expr=m.block1.b + m.block1.c == 3)
        m.block1.con_c = Constraint(expr=m.block1.c + m.block1.a == -4)
        return m

    @pytest.mark.unit
    def test_verify_active_variables_initialized(self, model):
        diag_tbx = DiagnosticsToolbox(model)
        with pytest.raises(
            RuntimeError,
            match=self.err_msg(3),
        ):
            diag_tbx._verify_active_variables_initialized()

        model.block1.deactivate()
        with pytest.raises(
            RuntimeError,
            match=self.err_msg(2),
        ):
            diag_tbx._verify_active_variables_initialized()

        model.con_y[1].deactivate()
        model.con_z.deactivate()
        # No uninitialized variables in remaining constraints
        diag_tbx._verify_active_variables_initialized()

    @pytest.mark.unit
    def test_error_handling(self, model):
        diag_tbx = DiagnosticsToolbox(model)
        methods = [
            "display_variables_with_extreme_jacobians",
            "display_constraints_with_extreme_jacobians",
            "display_extreme_jacobian_entries",
            "display_near_parallel_constraints",
            "display_near_parallel_variables",
            "display_constraints_with_mismatched_terms",
            "display_constraints_with_canceling_terms",
            "display_constraints_with_no_free_variables",
            "report_numerical_issues",
        ]
        for mthd in methods:
            mthd_obj = getattr(diag_tbx, mthd)
            with pytest.raises(
                RuntimeError,
                match=self.err_msg(3),
            ):
                mthd_obj()

    @pytest.mark.unit
    def test_display_problematic_constraint_terms(self, model):
        """
        This method needs to be split into a separate function
        because it takes a constraint as an argument.
        """
        diag_tbx = DiagnosticsToolbox(model)
        with pytest.raises(
            RuntimeError,
            match=self.err_msg(3),
        ):
            diag_tbx.display_problematic_constraint_terms(model.con_y[1])
        # The current implementation shows an error message whether or not
        # a problematic constraint is called.
        with pytest.raises(
            RuntimeError,
            match=self.err_msg(3),
        ):
            diag_tbx.display_problematic_constraint_terms(model.con_y[0])

    @pytest.mark.unit
    def test_display_variables_with_none_value_in_activated_constraints(self, model):
        diag_tbx = DiagnosticsToolbox(model)
        stream = StringIO()

        diag_tbx.display_variables_with_none_value_in_activated_constraints(
            stream=stream,
        )

        expected = """====================================================================================
The following variable(s) have a value of None:

    y[1]
    z
    block1.a

====================================================================================
"""
        assert stream.getvalue() == expected


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
    m.v = Var(initialize=-1e-4, bounds=(0, 1))

    # Just outside lower bound
    bounds_issue = _vars_violating_bounds(m, tolerance=1e-6)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 1

    bounds_issue = _vars_violating_bounds(m, tolerance=1e3)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    set_scaling_factor(m.v, 1e-10, overwrite=True)
    bounds_issue = _vars_violating_bounds(m, tolerance=1e3)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    set_scaling_factor(m.v, 1e10, overwrite=True)
    bounds_issue = _vars_violating_bounds(m, tolerance=1e-6)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 1

    # Just inside lower bound
    m.v.set_value(1e-4)

    set_scaling_factor(m.v, 1, overwrite=True)
    bounds_issue = _vars_violating_bounds(m, tolerance=1e-6)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    bounds_issue = _vars_violating_bounds(m, tolerance=1e3)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    set_scaling_factor(m.v, 1e-10, overwrite=True)
    bounds_issue = _vars_violating_bounds(m, tolerance=1e3)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    set_scaling_factor(m.v, 1e10, overwrite=True)
    bounds_issue = _vars_violating_bounds(m, tolerance=1e-6)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    # Just inside upper bound
    m.v.set_value(1 - 1e-4)

    set_scaling_factor(m.v, 1, overwrite=True)
    bounds_issue = _vars_violating_bounds(m, tolerance=1e-6)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    bounds_issue = _vars_violating_bounds(m, tolerance=1e3)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    set_scaling_factor(m.v, 1e-10, overwrite=True)
    bounds_issue = _vars_violating_bounds(m, tolerance=1e3)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    set_scaling_factor(m.v, 1e10, overwrite=True)
    bounds_issue = _vars_violating_bounds(m, tolerance=1e-6)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    # Just outside upper bound
    m.v.set_value(1 + 1e-4)

    set_scaling_factor(m.v, 1, overwrite=True)
    bounds_issue = _vars_violating_bounds(m, tolerance=1e-6)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 1

    bounds_issue = _vars_violating_bounds(m, tolerance=1e3)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    set_scaling_factor(m.v, 1e-10, overwrite=True)
    bounds_issue = _vars_violating_bounds(m, tolerance=1e3)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 0

    set_scaling_factor(m.v, 1e10, overwrite=True)
    bounds_issue = _vars_violating_bounds(m, tolerance=1e-6)
    assert isinstance(bounds_issue, ComponentSet)
    assert len(bounds_issue) == 1


@pytest.mark.unit
def test_vars_with_extreme_values():
    m = ConcreteModel()
    m.v1 = Var(initialize=1e-12)  # below zero
    m.v2 = Var(initialize=1e-8)  # small
    m.v3 = Var(initialize=1e-4)
    m.v4 = Var(initialize=1e0)
    m.v5 = Var(initialize=1e4)
    m.v6 = Var(initialize=1e8)
    m.v7 = Var(initialize=1e12)  # large

    xvars = _vars_with_extreme_values(m, large=1e9, small=1e-7, zero=1e-10)

    assert len(xvars) == 2
    for i in xvars:
        assert i.name in ["v2", "v7"]

    set_scaling_factor(m.v1, 1e3)
    set_scaling_factor(m.v2, 1e6)
    set_scaling_factor(m.v4, 1e-12)
    set_scaling_factor(m.v6, 1e3)
    set_scaling_factor(m.v7, 1e-12)

    xvars = _vars_with_extreme_values(m, large=1e9, small=1e-7, zero=1e-10)

    assert len(xvars) == 2
    for i in xvars:
        assert i.name in ["v1", "v6"]


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

    def test_with_greybox_variables(self):
        """non functional graybox model added to m fixture, to test DOFs

        GreyBoxModel has 3 inputs and 2 outputs calculated an unknown function,
        input a1 and a2 are bound by equality constraint through internal graybox model
        """

        class BasicGrayBox(ExternalGreyBoxModel):
            def input_names(self):
                return ["a1", "a2", "a3"]

            def output_names(self):
                return ["o1", "o2"]

            def equality_constraint_names(self):
                return ["a_sum"]

            def evaluate_equality_constraints(self):
                a1 = self._input_values[0]
                a2 = self._input_values[1]
                return [a1 * 0.5 + a2]

        m = ConcreteModel()

        m.gb = ExternalGreyBoxBlock(external_model=BasicGrayBox())
        m.gb_inactive = ExternalGreyBoxBlock(external_model=BasicGrayBox())
        m.gb_inactive.deactivate()
        m.a1 = Var(initialize=1)
        m.a1.fix()
        m.gb.inputs["a2"].unfix()
        m.gb.inputs["a3"].fix()
        m.a1_eq = Constraint(expr=m.a1 == m.gb.inputs["a1"])
        m.o1 = Var(initialize=1)
        m.o1_eq = Constraint(expr=m.o1 == m.gb.outputs["o1"])
        m.o1.fix()
        stats = _collect_model_statistics(m)
        for k in stats:
            print(k)
        print(stats, len(stats))
        m.display()
        tab = " " * 4
        assert len(stats) == 13
        assert stats[0] == f"{tab}Activated Blocks: 1 (Deactivated: 0)"
        assert (
            stats[1] == f"{tab}Free Variables in Activated Constraints: 4 (External: 0)"
        )
        assert stats[2] == f"{tab*2}Free Variables with only lower bounds: 0"
        assert stats[3] == f"{tab*2}Free Variables with only upper bounds: 0"
        assert stats[4] == f"{tab*2}Free Variables with upper and lower bounds: 0"
        assert (
            stats[5]
            == f"{tab}Fixed Variables in Activated Constraints: 3 (External: 0)"
        )
        assert stats[6] == f"{tab}Activated Equality Constraints: 5 (Deactivated: 3)"
        assert stats[7] == f"{tab}Activated Inequality Constraints: 0 (Deactivated: 0)"
        assert stats[8] == f"{tab}Activated Objectives: 0 (Deactivated: 0)"
        assert stats[9] == f"{tab}GreyBox Statistics"
        assert stats[10] == f"{tab*2}Activated GreyBox models: 1 (Deactivated: 1)"
        assert stats[11] == f"{tab*2}Activated GreyBox Equalities: 3 (Deactivated: 3)"
        assert (
            stats[12]
            == f"{tab*2}Free Variables in Activated GreyBox Equalities: 4 (Fixed: 1)"
        )


@pytest.mark.solver
class TestDiagnosticsToolbox:
    @pytest.mark.unit
    def test_invalid_model_type(self):
        with pytest.raises(
            TypeError,
            match=re.escape(
                "model argument must be an instance of a Pyomo BlockData object "
                "(either a scalar Block or an element of an indexed Block)."
            ),
        ):
            DiagnosticsToolbox(model="foo")

        # Check for indexed Blocks
        m = ConcreteModel()
        m.s = Set(initialize=[1, 2, 3])
        m.b = Block(m.s)

        with pytest.raises(
            TypeError,
            match=re.escape(
                "model argument must be an instance of a Pyomo BlockData object "
                "(either a scalar Block or an element of an indexed Block)."
            ),
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

        # Presolver identifies problem as trivially infeasible (correctly). Turn off presolve.
        solver = get_solver("ipopt_v2", writer_config={"linear_presolve": False})
        solver.solve(m)

        return m

    @pytest.mark.component
    def test_with_grey_box(self):

        class BasicGrayBox(ExternalGreyBoxModel):
            def input_names(self):
                return ["a1", "a2", "a3"]

            def output_names(self):
                return ["o1", "o2"]

            def equality_constraint_names(self):
                return ["a_sum"]

            def evaluate_equality_constraints(self):
                a1 = self._input_values[0]
                a2 = self._input_values[1]
                return [a1 * 0.5 + a2]

        m = ConcreteModel()

        m.gb = ExternalGreyBoxBlock(external_model=BasicGrayBox())
        with pytest.raises(NotImplementedError):
            DiagnosticsToolbox(model=m)

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
        model.v2 = Var(initialize=1)
        model.v3 = Var(initialize=1)

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
        model.v2 = Var(initialize=1)
        model.v3 = Var(initialize=1)

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
        model.v2 = Var(initialize=1)
        model.v3 = Var(initialize=1)

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
    def test_display_near_parallel_constraints(self):
        model = ConcreteModel()
        model.v1 = Var(initialize=1e-8)
        model.v2 = Var(initialize=1)
        model.v3 = Var(initialize=1)

        model.c1 = Constraint(expr=model.v1 == model.v2)
        model.c2 = Constraint(expr=model.v1 == 1e-8 * model.v3)
        model.c3 = Constraint(expr=1e8 * model.v1 + 1e10 * model.v2 == 1e-6 * model.v3)
        model.c4 = Constraint(expr=-model.v1 == -0.99999 * model.v2)

        dt = DiagnosticsToolbox(model=model)

        stream = StringIO()
        dt.display_near_parallel_constraints(stream)

        expected = """====================================================================================
The following pairs of constraints are nearly parallel:

    c1, c4

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_display_near_parallel_variables(self):
        model = ConcreteModel()
        model.v1 = Var(initialize=1e-8)
        model.v2 = Var(initialize=1)
        model.v3 = Var(initialize=1)
        model.v4 = Var(initialize=1)

        model.c1 = Constraint(expr=1e-8 * model.v1 == 1e-8 * model.v2 - 1e-8 * model.v4)
        model.c2 = Constraint(expr=1e-8 * model.v1 + 1e-8 * model.v4 == model.v3)
        model.c3 = Constraint(
            expr=1e3 * (model.v1 + model.v4) + 1e3 * model.v2 == model.v3
        )

        dt = DiagnosticsToolbox(model=model)

        stream = StringIO()
        dt.display_near_parallel_variables(stream)

        expected = """====================================================================================
The following pairs of variables are nearly parallel:

    v1, v2
    v1, v4
    v2, v4

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_collect_constraint_mismatches(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=2)
        m.v2 = Var(initialize=3)

        # Constraint with no free variables
        m.c1 = Constraint(expr=m.v1 == m.v2)
        m.v1.fix()
        m.v2.fix()

        # Constraint with mismatched terms
        m.v3 = Var(initialize=10)
        m.v4 = Var(initialize=10)
        m.v5 = Var(initialize=1e-6)
        m.c2 = Constraint(expr=m.v3 == m.v4 + m.v5)

        # Constraint with cancellation
        m.v6 = Var(initialize=10)
        m.c3 = Constraint(expr=m.v6 == 10 + m.v3 - m.v4)

        dt = DiagnosticsToolbox(model=m)

        mismatch, cancellation, constant = dt._collect_constraint_mismatches()

        assert mismatch == ["c2: 1 mismatched term(s)"]
        assert cancellation == [
            "c2: 1 potential canceling term(s)",
            "c3: 1 potential canceling term(s)",
        ]
        assert constant == ["c1"]

    @pytest.mark.component
    def test_display_problematic_constraint_terms(self):
        m = ConcreteModel()

        # Constraint with mismatched terms
        m.v3 = Var(initialize=10)
        m.v4 = Var(initialize=10)
        m.v5 = Var(initialize=1e-6)
        m.c2 = Constraint(expr=m.v3 == m.v4 + m.v5)

        dt = DiagnosticsToolbox(model=m)

        stream = StringIO()
        dt.display_problematic_constraint_terms(m.c2, stream=stream)

        expected = """====================================================================================
The following terms in c2 are potentially problematic:

    Mismatch in v4 + v5 (Max 10, Min 1e-06)
    Cancellation in v3  ==  v4 + v5. Terms 1 (-10), 2 (10)

====================================================================================
"""

        # Variable scaling doesn't affect anything
        assert stream.getvalue() == expected

        set_scaling_factor(m.v5, 1e-5, overwrite=True)

        stream = StringIO()
        dt.display_problematic_constraint_terms(m.c2, stream=stream)

        assert stream.getvalue() == expected

        set_scaling_factor(m.v5, 1e5, overwrite=True)

        stream = StringIO()
        dt.display_problematic_constraint_terms(m.c2, stream=stream)

        assert stream.getvalue() == expected

        # Constraint scaling affects only the zero tolerance
        set_scaling_factor(m.c2, 1e-5, overwrite=True)

        stream = StringIO()
        dt.display_problematic_constraint_terms(m.c2, stream=stream)

        assert stream.getvalue() == expected

        set_scaling_factor(m.c2, 1e5, overwrite=True)

        expected = """====================================================================================
The following terms in c2 are potentially problematic:


====================================================================================
"""
        stream = StringIO()
        dt.display_problematic_constraint_terms(m.c2, stream=stream)

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_display_problematic_constraint_terms_limited(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=10)
        m.v2 = Var(initialize=-10)
        m.v4 = Var(initialize=10)

        # Constraint with cancellation
        m.v6 = Var(initialize=10)
        m.c1 = Constraint(expr=m.v6 == 10 - m.v4 + m.v1 + m.v2)

        dt = DiagnosticsToolbox(model=m)

        stream = StringIO()
        dt.display_problematic_constraint_terms(m.c1, stream=stream)

        expected = """====================================================================================
The following terms in c1 are potentially problematic:

    Cancellation in v6  ==  10 - v4 + v1 + v2. Terms 1 (-10), 2 (10)
    Cancellation in v6  ==  10 - v4 + v1 + v2. Terms 1 (-10), 4 (10)
    Cancellation in v6  ==  10 - v4 + v1 + v2. Terms 2 (10), 3 (-10)
    Cancellation in v6  ==  10 - v4 + v1 + v2. Terms 2 (10), 5 (-10)
    Cancellation in v6  ==  10 - v4 + v1 + v2. Terms 3 (-10), 4 (10)

Number of canceling terms per node limited to 5.
====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_display_problematic_constraint_terms_indexed_error(self):
        m = ConcreteModel()
        m.s = Set(initialize=[1, 2])
        m.v1 = Var(m.s, initialize=2)
        m.v2 = Var(m.s, initialize=3)

        def _c_rule(b, i):
            return b.v1[i] == b.v2[i]

        m.c1 = Constraint(m.s, rule=_c_rule)

        dt = DiagnosticsToolbox(model=m)

        # Check for indexed constraints
        with pytest.raises(
            TypeError,
            match=re.escape(
                "c1 is an IndexedConstraint. Please provide "
                "an individual element of c1 (ConstraintData) "
                "to be examined for problematic terms."
            ),
        ):
            dt.display_problematic_constraint_terms(m.c1)

        # Check for not a constraints
        with pytest.raises(
            TypeError, match="v1 is not an instance of a Pyomo Constraint."
        ):
            dt.display_problematic_constraint_terms(m.v1)

    @pytest.mark.component
    def test_display_problematic_constraint_terms_named_expression(self):
        m = ConcreteModel()

        # Constraint with mismatched terms
        m.v1 = Var(initialize=10)
        m.v2 = Var(initialize=10)
        m.v3 = Var(initialize=1e-6)

        m.e1 = Expression(expr=(m.v1 - m.v2))

        m.c2 = Constraint(expr=m.v3 == m.e1**2)

        dt = DiagnosticsToolbox(model=m)

        stream = StringIO()
        dt.display_problematic_constraint_terms(m.c2, stream=stream)

        expected = """====================================================================================
The following terms in c2 are potentially problematic:

    Cancellation in Expression e1. Terms 1 (10), 2 (-10)

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_display_constraints_with_mismatched_terms(self):
        m = ConcreteModel()
        # Constraint with mismatched terms
        m.v3 = Var(initialize=10)
        m.v4 = Var(initialize=10)
        m.v5 = Var(initialize=1e-6)
        m.c2 = Constraint(expr=m.v3 == m.v4 + m.v5)

        dt = DiagnosticsToolbox(model=m)

        stream = StringIO()
        dt.display_constraints_with_mismatched_terms(stream=stream)

        expected = """====================================================================================
The following constraints have mismatched terms:

    c2: 1 mismatched term(s)

Call display_problematic_constraint_terms(constraint) for more information.
====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_display_constraints_with_canceling_terms(self):
        m = ConcreteModel()
        # Constraint with mismatched terms
        m.v3 = Var(initialize=10)
        m.v4 = Var(initialize=10)

        # Constraint with cancellation
        m.v6 = Var(initialize=10)
        m.c3 = Constraint(expr=m.v6 == 10 + m.v3 - m.v4)

        dt = DiagnosticsToolbox(model=m)

        stream = StringIO()
        dt.display_constraints_with_canceling_terms(stream=stream)

        expected = """====================================================================================
The following constraints have canceling terms:

    c3: 1 potential canceling term(s)

Call display_problematic_constraint_terms(constraint) for more information.
====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_display_constraints_with_no_free_variables(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=2)
        m.v2 = Var(initialize=3)

        # Constraint with no free variables
        m.c1 = Constraint(expr=m.v1 == m.v2)
        m.v1.fix()
        m.v2.fix()

        dt = DiagnosticsToolbox(model=m)

        stream = StringIO()
        dt.display_constraints_with_no_free_variables(stream=stream)

        expected = """====================================================================================
The following constraints have no free variables:

    c1

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

        assert len(next_steps) == 3
        assert "display_constraints_with_large_residuals()" in next_steps
        assert "compute_infeasibility_explanation()" in next_steps
        assert "display_variables_at_or_outside_bounds()" in next_steps

    @pytest.mark.component
    def test_collect_numerical_warnings_corrected(self, model):
        m = model.clone()

        # Fix numerical issues
        m.b.v3.setlb(-5)
        m.b.v5.setub(10)

        # Presolver identifies problem as trivially infeasible (correctly). Turn off presolve.
        solver = get_solver("ipopt_v2", writer_config={"linear_presolve": False})
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

        assert len(warnings) == 4
        assert (
            "WARNING: 2 Variables with extreme Jacobian values (<1.0E-08 or >1.0E+08)"
            in warnings
        )
        assert (
            "WARNING: 1 Constraint with extreme Jacobian values (<1.0E-08 or >1.0E+08)"
            in warnings
        )
        assert "WARNING: 1 Constraint with large residuals (>1.0E-05)" in warnings

        assert len(next_steps) == 5
        assert "display_variables_with_extreme_jacobians()" in next_steps
        assert "display_constraints_with_extreme_jacobians()" in next_steps
        assert "display_constraints_with_large_residuals()" in next_steps
        assert "compute_infeasibility_explanation()" in next_steps

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

        with pytest.raises(
            AssertionError, match=re.escape("Structural issues found (1).")
        ):
            dt.assert_no_structural_warnings()

        # Fix units issue
        m.b.del_component(m.b.c1)
        m.b.c1 = Constraint(expr=m.v1 + m.b.v2 == 10 * units.m)
        dt.assert_no_structural_warnings()

    @pytest.mark.component
    def test_assert_no_numerical_warnings(self, model):
        m = model.clone()
        dt = DiagnosticsToolbox(model=m.b)

        with pytest.raises(
            AssertionError, match=re.escape("Numerical issues found (2).")
        ):
            dt.assert_no_numerical_warnings()

        # Fix numerical issues
        m.b.v3.setlb(-5)

        # Presolver identifies problem as trivially infeasible (correctly). Turn off presolve.
        solver = get_solver("ipopt_v2", writer_config={"linear_presolve": False})
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
    def test_report_structural_issues_ok(self):
        m = ConcreteModel()

        m.v1 = Var(initialize=1)
        m.v2 = Var(initialize=2)
        m.v3 = Var(initialize=3)

        m.c1 = Constraint(expr=2 * m.v1 == m.v2)
        m.c2 = Constraint(expr=m.v1 + m.v2 == m.v3)
        m.c3 = Constraint(expr=m.v1 == 1)

        dt = DiagnosticsToolbox(model=m)

        stream = StringIO()
        dt.report_structural_issues(stream)

        expected = """====================================================================================
Model Statistics

        Activated Blocks: 1 (Deactivated: 0)
        Free Variables in Activated Constraints: 3 (External: 0)
            Free Variables with only lower bounds: 0
            Free Variables with only upper bounds: 0
            Free Variables with upper and lower bounds: 0
        Fixed Variables in Activated Constraints: 0 (External: 0)
        Activated Equality Constraints: 3 (Deactivated: 0)
        Activated Inequality Constraints: 0 (Deactivated: 0)
        Activated Objectives: 0 (Deactivated: 0)

------------------------------------------------------------------------------------
0 WARNINGS

    No warnings found!

------------------------------------------------------------------------------------
0 Cautions

    No cautions found!

------------------------------------------------------------------------------------
Suggested next steps:

    Try to initialize/solve your model and then call report_numerical_issues()

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_report_numerical_issues_ok(self):
        m = ConcreteModel()

        m.v1 = Var(initialize=1)
        m.v2 = Var(initialize=2)
        m.v3 = Var(initialize=3)

        m.c1 = Constraint(expr=2 * m.v1 == m.v2)
        m.c2 = Constraint(expr=m.v1 + m.v2 == m.v3)
        m.c3 = Constraint(expr=m.v1 == 1)

        dt = DiagnosticsToolbox(model=m)

        stream = StringIO()
        dt.report_numerical_issues(stream)

        expected = """====================================================================================
Model Statistics

    Jacobian Condition Number: 1.237E+01

------------------------------------------------------------------------------------
0 WARNINGS

    No warnings found!

------------------------------------------------------------------------------------
0 Cautions

    No cautions found!

------------------------------------------------------------------------------------
Suggested next steps:

    If you still have issues converging your model consider:

        prepare_degeneracy_hunter()
        prepare_svd_toolbox()

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_report_numerical_issues_exactly_singular(self):
        m = ConcreteModel()
        m.x = Var([1, 2], initialize=1.0)
        m.eq = Constraint(PositiveIntegers)
        m.eq[1] = m.x[1] * m.x[2] == 1.5
        m.eq[2] = m.x[2] * m.x[1] == 1.5
        m.obj = Objective(expr=m.x[1] ** 2 + 2 * m.x[2] ** 2)

        dt = DiagnosticsToolbox(m)
        dt.report_numerical_issues()

        stream = StringIO()
        dt.report_numerical_issues(stream)

        expected = """====================================================================================
Model Statistics

    Jacobian Condition Number: Undefined (Exactly Singular)

------------------------------------------------------------------------------------
3 WARNINGS

    WARNING: 2 Constraints with large residuals (>1.0E-05)
    WARNING: 1 pair of constraints are parallel (to tolerance 1.0E-08)
    WARNING: 1 pair of variables are parallel (to tolerance 1.0E-08)

------------------------------------------------------------------------------------
0 Cautions

    No cautions found!

------------------------------------------------------------------------------------
Suggested next steps:

    display_constraints_with_large_residuals()
    compute_infeasibility_explanation()
    display_near_parallel_constraints()
    display_near_parallel_variables()

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
    compute_infeasibility_explanation()
    display_variables_at_or_outside_bounds()

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_report_numerical_issues_cancellation(self):
        model = ConcreteModel()

        model.v1 = Var(initialize=1)
        model.v2 = Var(initialize=2)
        model.v3 = Var(initialize=1e-8)

        # Non-problematic constraints
        model.c1 = Constraint(expr=2 * model.v1 == model.v2)

        # Problematic constraints
        model.c10 = Constraint(expr=0 == exp(model.v1 - 0.5 * model.v2))
        model.c11 = Constraint(expr=0 == model.v1 - 0.5 * model.v2 + model.v3)

        dt = DiagnosticsToolbox(model=model)

        stream = StringIO()
        dt.report_numerical_issues(stream)

        expected = """====================================================================================
Model Statistics

    Jacobian Condition Number: Undefined (Exactly Singular)

------------------------------------------------------------------------------------
3 WARNINGS

    WARNING: 1 Constraint with large residuals (>1.0E-05)
    WARNING: 1 pair of constraints are parallel (to tolerance 1.0E-08)
    WARNING: 1 pair of variables are parallel (to tolerance 1.0E-08)

------------------------------------------------------------------------------------
3 Cautions

    Caution: 1 Variable with value close to zero (tol=1.0E-08)
    Caution: 1 Constraint with mismatched terms
    Caution: 1 Constraint with potential cancellation of terms

------------------------------------------------------------------------------------
Suggested next steps:

    display_constraints_with_large_residuals()
    compute_infeasibility_explanation()
    display_near_parallel_constraints()
    display_near_parallel_variables()

====================================================================================
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_report_numerical_issues_jacobian(self):
        model = ConcreteModel()
        model.v1 = Var(initialize=1e-8)
        model.v2 = Var(initialize=0)
        model.v3 = Var(initialize=0)

        model.c1 = Constraint(expr=1e-2 * model.v1 == model.v2)
        model.c2 = Constraint(expr=1e-2 * model.v1 == 1e-8 * model.v3)
        model.c3 = Constraint(expr=1e8 * model.v1 + 1e10 * model.v2 == 1e-6 * model.v3)

        dt = DiagnosticsToolbox(model=model)

        stream = StringIO()
        dt.report_numerical_issues(stream)

        expected = """====================================================================================
Model Statistics

    Jacobian Condition Number: 1.118E+18

------------------------------------------------------------------------------------
4 WARNINGS

    WARNING: 1 Constraint with large residuals (>1.0E-05)
    WARNING: 2 Variables with extreme Jacobian values (<1.0E-08 or >1.0E+08)
    WARNING: 1 Constraint with extreme Jacobian values (<1.0E-08 or >1.0E+08)
    WARNING: 3 pairs of variables are parallel (to tolerance 1.0E-08)

------------------------------------------------------------------------------------
4 Cautions

    Caution: 3 Variables with value close to zero (tol=1.0E-08)
    Caution: 3 Variables with extreme Jacobian values (<1.0E-04 or >1.0E+04)
    Caution: 1 Constraint with extreme Jacobian values (<1.0E-04 or >1.0E+04)
    Caution: 4 extreme Jacobian Entries (<1.0E-04 or >1.0E+04)

------------------------------------------------------------------------------------
Suggested next steps:

    display_constraints_with_large_residuals()
    compute_infeasibility_explanation()
    display_variables_with_extreme_jacobians()
    display_constraints_with_extreme_jacobians()
    display_near_parallel_variables()

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

    @pytest.mark.component
    def test_with_grey_box(self):

        class BasicGrayBox(ExternalGreyBoxModel):
            def input_names(self):
                return ["a1", "a2", "a3"]

            def output_names(self):
                return ["o1", "o2"]

            def equality_constraint_names(self):
                return ["a_sum"]

            def evaluate_equality_constraints(self):
                a1 = self._input_values[0]
                a2 = self._input_values[1]
                return [a1 * 0.5 + a2]

        m = ConcreteModel()

        m.gb = ExternalGreyBoxBlock(external_model=BasicGrayBox())
        with pytest.raises(NotImplementedError):
            SVDToolbox(model=m)

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
            match=re.escape(
                "variable argument must be an instance of a Pyomo VarData "
                "object (got foo)."
            ),
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
            match=re.escape(
                "constraint argument must be an instance of a Pyomo ConstraintData "
                "object (got foo)."
            ),
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

    @pytest.mark.component
    def test_with_grey_box(self):

        class BasicGrayBox(ExternalGreyBoxModel):
            def input_names(self):
                return ["a1", "a2", "a3"]

            def output_names(self):
                return ["o1", "o2"]

            def equality_constraint_names(self):
                return ["a_sum"]

            def evaluate_equality_constraints(self):
                a1 = self._input_values[0]
                a2 = self._input_values[1]
                return [a1 * 0.5 + a2]

        m = ConcreteModel()

        m.gb = ExternalGreyBoxBlock(external_model=BasicGrayBox())
        with pytest.raises(NotImplementedError):
            DegeneracyHunter2(model=m)

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

        assert dh.degenerate_set == {
            model.con2: value(dh.candidates_milp.nu[0]),
            model.con5: value(dh.candidates_milp.nu[1]),
        }

        assert abs(value(dh.candidates_milp.nu[0])) == pytest.approx(1e-05, rel=1e-5)
        assert abs(value(dh.candidates_milp.nu[1])) == pytest.approx(1e-05, rel=1e-5)

        # One must be positive and one must be negative, so produce will be negative
        assert value(
            dh.candidates_milp.nu[0] * dh.candidates_milp.nu[1]
        ) == pytest.approx(-1e-10, rel=1e-5)

        assert (
            value(
                dh.candidates_milp.y_pos[0]
                + dh.candidates_milp.y_pos[1]
                + dh.candidates_milp.y_neg[0]
                + dh.candidates_milp.y_neg[1]
            )
            >= 1
        )

        assert value(dh.candidates_milp.abs_nu[0]) == pytest.approx(1e-05, rel=1e-5)
        assert value(dh.candidates_milp.abs_nu[1]) == pytest.approx(1e-05, rel=1e-5)

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

    # TODO does this test function have the exact same name as the one above?
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


ca_dict = {
    "specification": {
        "inputs": {
            "v2": {
                "pyomo_path": "v2",
                "lower": 2,
                "upper": 6,
            },
        },
        "sampling_method": "UniformSampling",
        "sample_size": [2],
        "samples": {
            "index": [0, 1],
            "columns": ["v2"],
            "data": [[2.0], [6.0]],
            "index_names": [None],
            "column_names": [None],
        },
    },
    "results": {
        0: {
            "success": True,
            "results": 2,
        },
        1: {
            "success": True,
            "results": 6,
        },
    },
}

ca_res = {
    "specification": {
        "inputs": {"v2": {"pyomo_path": "v2", "lower": 2, "upper": 6}},
        "sampling_method": "UniformSampling",
        "sample_size": [2],
        "samples": {
            "index": [0, 1],
            "columns": ["v2"],
            "data": [[2.0], [6.0]],
            "index_names": [None],
            "column_names": [None],
        },
    },
    "results": {
        0: {
            "success": False,
            "results": {
                "iters": 7,
                "iters_in_restoration": 4,
                "iters_w_regularization": 0,
                "time": 0.0,
                "numerical_issues": True,
            },
        },
        1: {
            "success": False,
            "results": {
                "iters": 7,
                "iters_in_restoration": 4,
                "iters_w_regularization": 0,
                "time": 0.0,
                "numerical_issues": True,
            },
        },
    },
}


class TestIpoptConvergenceAnalysis:
    @pytest.fixture
    def model(self):
        m = ConcreteModel()

        m.v1 = Var(bounds=(None, 1))
        m.v2 = Var()
        m.c = Constraint(expr=m.v1 == m.v2)

        m.v2.fix(0)

        return m

    @pytest.mark.component
    def test_with_grey_box(self):

        class BasicGrayBox(ExternalGreyBoxModel):
            def input_names(self):
                return ["a1", "a2", "a3"]

            def output_names(self):
                return ["o1", "o2"]

            def equality_constraint_names(self):
                return ["a_sum"]

            def evaluate_equality_constraints(self):
                a1 = self._input_values[0]
                a2 = self._input_values[1]
                return [a1 * 0.5 + a2]

        m = ConcreteModel()

        m.gb = ExternalGreyBoxBlock(external_model=BasicGrayBox())
        with pytest.raises(NotImplementedError):
            IpoptConvergenceAnalysis(model=m)

    @pytest.mark.unit
    def test_init(self, model):
        ca = IpoptConvergenceAnalysis(model)

        assert ca._model is model
        assert isinstance(ca._psweep, SequentialSweepRunner)
        assert isinstance(ca.results, dict)
        assert ca.config.input_specification is None
        assert ca.config.solver_options is None

    @pytest.mark.unit
    def test_init_solver_options(self, model):
        solver_obj = SolverFactory("ipopt")
        options = solver_obj.options["max_iter"] = 100
        ca = IpoptConvergenceAnalysis(
            model, solver_obj=solver_obj, solver_options=options
        )

        assert ca._model is model
        assert isinstance(ca._psweep, SequentialSweepRunner)
        assert isinstance(ca.results, dict)
        assert ca.config.input_specification is None
        assert ca.config.solver_options is not None

    @pytest.mark.unit
    def test_build_model(self, model):
        ca = IpoptConvergenceAnalysis(model)

        clone = ca._build_model()

        assert clone is not model
        assert isinstance(clone.v1, Var)
        assert isinstance(clone.v2, Var)
        assert isinstance(clone.c, Constraint)

    @pytest.mark.unit
    def test_parse_ipopt_output(self, model):
        ca = IpoptConvergenceAnalysis(model)

        fname = os.path.join(currdir, "ipopt_output.txt")
        iters, restoration, regularization, time = ca._parse_ipopt_output(fname)

        assert iters == 43
        assert restoration == 39
        assert regularization == 4
        assert time == 0.016 + 0.035

    @pytest.mark.component
    @pytest.mark.solver
    def test_run_ipopt_with_stats(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1)
        m.c1 = Constraint(expr=m.v1 == 4)

        ca = IpoptConvergenceAnalysis(m)
        solver = SolverFactory("ipopt")

        (
            status,
            iters,
            iters_in_restoration,
            iters_w_regularization,
            time,
        ) = ca._run_ipopt_with_stats(m, solver)

        assert_optimal_termination(status)
        assert iters == 1
        assert iters_in_restoration == 0
        assert iters_w_regularization == 0
        assert isinstance(time, float)

    @pytest.mark.component
    @pytest.mark.solver
    def test_run_model(self, model):
        ca = IpoptConvergenceAnalysis(model)

        model.v2.fix(0.5)

        solver = SolverFactory("ipopt")

        success, run_stats = ca._run_model(model, solver)

        assert success
        assert value(model.v1) == pytest.approx(0.5, rel=1e-8)

        assert len(run_stats) == 4
        assert run_stats[0] == 1
        assert run_stats[1] == 0
        assert run_stats[2] == 0

    @pytest.mark.unit
    def test_build_outputs(self, model):
        ca = IpoptConvergenceAnalysis(model)

        model.v1.set_value(0.5)
        model.v2.fix(0.5)

        results = ca._build_outputs(model, (1, 2, 3, 4))

        assert results == {
            "iters": 1,
            "iters_in_restoration": 2,
            "iters_w_regularization": 3,
            "time": 4,
            "numerical_issues": False,
        }

    @pytest.mark.unit
    def test_build_outputs_with_warnings(self, model):
        ca = IpoptConvergenceAnalysis(model)

        model.v1.set_value(4)
        model.v2.fix(0.5)

        results = ca._build_outputs(model, (1, 2, 3, 4))

        assert results == {
            "iters": 1,
            "iters_in_restoration": 2,
            "iters_w_regularization": 3,
            "time": 4,
            "numerical_issues": True,
        }

    @pytest.mark.unit
    def test_recourse(self, model):
        ca = IpoptConvergenceAnalysis(model)

        assert ca._recourse(model) == {
            "iters": -1,
            "iters_in_restoration": -1,
            "iters_w_regularization": -1,
            "time": -1,
            "numerical_issues": -1,
        }

    @pytest.mark.integration
    @pytest.mark.solver
    def test_run_convergence_analysis(self, model):
        spec = ParameterSweepSpecification()
        spec.add_sampled_input("v2", lower=0, upper=3)
        spec.set_sampling_method(UniformSampling)
        spec.set_sample_size([4])

        ca = IpoptConvergenceAnalysis(model, input_specification=spec)

        ca.run_convergence_analysis()

        assert isinstance(ca.results, dict)
        assert len(ca.results) == 4

        # Ignore time, as it is too noisy to test
        # Sample 0 should solve cleanly
        assert ca.results[0]["success"]
        assert ca.results[0]["results"]["iters"] == 0
        assert ca.results[0]["results"]["iters_in_restoration"] == 0
        assert ca.results[0]["results"]["iters_w_regularization"] == 0
        assert not ca.results[0]["results"]["numerical_issues"]

        # Sample 1 should solve, but have issues due to bound on v1
        assert ca.results[1]["success"]
        assert ca.results[1]["results"]["iters"] == pytest.approx(3, abs=1)
        assert ca.results[1]["results"]["iters_in_restoration"] == 0
        assert ca.results[1]["results"]["iters_w_regularization"] == 0
        assert ca.results[1]["results"]["numerical_issues"]

        # Other iterations should fail due to bound
        assert not ca.results[2]["success"]
        assert ca.results[2]["results"]["iters"] == pytest.approx(7, abs=1)
        assert ca.results[2]["results"]["iters_in_restoration"] == pytest.approx(
            4, abs=1
        )
        assert ca.results[2]["results"]["iters_w_regularization"] == 0
        assert ca.results[2]["results"]["numerical_issues"]

        assert not ca.results[3]["success"]
        assert ca.results[3]["results"]["iters"] == pytest.approx(8, abs=1)
        assert ca.results[3]["results"]["iters_in_restoration"] == pytest.approx(
            5, abs=1
        )
        assert ca.results[3]["results"]["iters_w_regularization"] == 0
        assert ca.results[3]["results"]["numerical_issues"]

    @pytest.fixture(scope="class")
    def ca_with_results(self):
        spec = ParameterSweepSpecification()
        spec.set_sampling_method(UniformSampling)
        spec.add_sampled_input("v2", 2, 6)
        spec.set_sample_size([2])
        spec.generate_samples()

        ca = IpoptConvergenceAnalysis(
            model=ConcreteModel(),
            input_specification=spec,
        )

        ca._psweep._results = {
            0: {"success": True, "results": 2},
            1: {"success": True, "results": 6},
        }

        return ca

    @pytest.mark.unit
    def test_report_convergence_summary(self):
        stream = StringIO()

        ca = IpoptConvergenceAnalysis(
            model=ConcreteModel(),
        )

        ca._psweep._results = {
            0: {
                "success": True,
                "results": {
                    "iters_in_restoration": 1,
                    "iters_w_regularization": 0,
                    "numerical_issues": 10,
                },
            },
            1: {
                "success": True,
                "results": {
                    "iters_in_restoration": 0,
                    "iters_w_regularization": 5,
                    "numerical_issues": 5,
                },
            },
            2: {
                "success": False,
                "results": {
                    "iters_in_restoration": 0,
                    "iters_w_regularization": 0,
                    "numerical_issues": 0,
                },
            },
        }

        ca.report_convergence_summary(stream)

        expected = """Successes: 2, Failures 1 (66.66666666666667%)
Runs with Restoration: 1
Runs with Regularization: 1
Runs with Numerical Issues: 2
"""

        assert stream.getvalue() == expected

    @pytest.mark.component
    def test_to_dict(self, ca_with_results):
        outdict = ca_with_results.to_dict()
        assert outdict == ca_dict

    @pytest.mark.unit
    def test_from_dict(self):
        ca = IpoptConvergenceAnalysis(
            model=ConcreteModel(),
        )

        ca.from_dict(ca_dict)

        input_spec = ca._psweep.get_input_specification()

        assert isinstance(input_spec, ParameterSweepSpecification)
        assert len(input_spec.inputs) == 1

        assert input_spec.sampling_method is UniformSampling
        assert isinstance(input_spec.samples, DataFrame)
        assert input_spec.sample_size == [2]

        assert isinstance(ca.results, dict)
        assert len(ca.results) == 2

        for i in [0, 1]:
            assert ca.results[i]["success"]
            assert ca.results[i]["results"] == 2 + i * 4

    @pytest.mark.component
    def test_to_json_file(self, ca_with_results):
        temp_context = TempfileManager.new_context()
        tmpfile = temp_context.create_tempfile(suffix=".json")

        ca_with_results.to_json_file(tmpfile)

        with open(tmpfile, "r") as f:
            lines = f.read()
        f.close()

        expected = """{
   "specification": {
      "inputs": {
         "v2": {
            "pyomo_path": "v2",
            "lower": 2,
            "upper": 6
         }
      },
      "sampling_method": "UniformSampling",
      "sample_size": [
         2
      ],
      "samples": {
         "index": [
            0,
            1
         ],
         "columns": [
            "v2"
         ],
         "data": [
            [
               2.0
            ],
            [
               6.0
            ]
         ],
         "index_names": [
            null
         ],
         "column_names": [
            null
         ]
      }
   },
   "results": {
      "0": {
         "success": true,
         "results": 2
      },
      "1": {
         "success": true,
         "results": 6
      }
   }
}"""

        assert lines == expected

        # Check for clean up
        temp_context.release(remove=True)
        assert not os.path.exists(tmpfile)

    @pytest.mark.unit
    def test_load_from_json_file(self):
        fname = os.path.join(currdir, "load_psweep.json")

        ca = IpoptConvergenceAnalysis(
            model=ConcreteModel(),
        )
        ca.from_json_file(fname)

        input_spec = ca._psweep.get_input_specification()

        assert isinstance(input_spec, ParameterSweepSpecification)
        assert len(input_spec.inputs) == 1

        assert input_spec.sampling_method is UniformSampling
        assert isinstance(input_spec.samples, DataFrame)
        assert input_spec.sample_size == [2]

        assert isinstance(ca.results, dict)
        assert len(ca.results) == 2

        for i in [0, 1]:
            assert ca.results[i]["success"]
            assert ca.results[i]["results"] == 2 + i * 4

    @pytest.mark.integration
    @pytest.mark.solver
    def test_run_convergence_analysis_from_dict(self, model):
        ca = IpoptConvergenceAnalysis(
            model=model,
        )
        ca.run_convergence_analysis_from_dict(ca_dict)

        input_spec = ca._psweep.get_input_specification()

        assert isinstance(input_spec, ParameterSweepSpecification)
        assert len(input_spec.inputs) == 1

        assert input_spec.sampling_method is UniformSampling
        assert isinstance(input_spec.samples, DataFrame)
        assert input_spec.sample_size == [2]

        assert isinstance(ca.results, dict)
        assert len(ca.results) == 2

        assert not ca.results[0]["success"]
        assert ca.results[0]["results"]["iters"] == pytest.approx(7, abs=1)
        assert ca.results[0]["results"]["iters_in_restoration"] == pytest.approx(
            4, abs=1
        )
        assert ca.results[0]["results"]["iters_w_regularization"] == 0
        assert ca.results[0]["results"]["numerical_issues"]

        assert not ca.results[1]["success"]
        assert ca.results[1]["results"]["iters"] == pytest.approx(7, abs=1)
        assert ca.results[1]["results"]["iters_in_restoration"] == pytest.approx(
            4, abs=1
        )
        assert ca.results[1]["results"]["iters_w_regularization"] == 0
        assert ca.results[1]["results"]["numerical_issues"]

    @pytest.mark.integration
    @pytest.mark.solver
    def test_run_convergence_analysis_from_file(self, model):
        fname = os.path.join(currdir, "load_psweep.json")

        ca = IpoptConvergenceAnalysis(
            model=model,
        )
        ca.run_convergence_analysis_from_file(fname)

        input_spec = ca._psweep.get_input_specification()

        assert isinstance(input_spec, ParameterSweepSpecification)
        assert len(input_spec.inputs) == 1

        assert input_spec.sampling_method is UniformSampling
        assert isinstance(input_spec.samples, DataFrame)
        assert input_spec.sample_size == [2]

        assert isinstance(ca.results, dict)
        assert len(ca.results) == 2

        assert not ca.results[0]["success"]
        assert ca.results[0]["results"]["iters"] == pytest.approx(7, abs=1)
        assert ca.results[0]["results"]["iters_in_restoration"] == pytest.approx(
            4, abs=1
        )
        assert ca.results[0]["results"]["iters_w_regularization"] == 0
        assert ca.results[0]["results"]["numerical_issues"]

        assert not ca.results[1]["success"]
        assert ca.results[1]["results"]["iters"] == pytest.approx(7, abs=1)
        assert ca.results[1]["results"]["iters_in_restoration"] == pytest.approx(
            4, abs=1
        )
        assert ca.results[1]["results"]["iters_w_regularization"] == 0
        assert ca.results[1]["results"]["numerical_issues"]

    @pytest.fixture(scope="class")
    def conv_anal(self):
        ca = IpoptConvergenceAnalysis(
            model=ConcreteModel(),
        )
        ca.from_dict(ca_res)

        return ca

    @pytest.mark.unit
    def test_compare_results_to_dict_ok(self, conv_anal):
        diffs = conv_anal._compare_results_to_dict(ca_res)

        assert diffs["success"] == []
        assert diffs["iters"] == []
        assert diffs["iters_in_restoration"] == []
        assert diffs["iters_w_regularization"] == []
        assert diffs["numerical_issues"] == []

    @pytest.mark.unit
    def test_compare_results_to_dict_success(self, conv_anal):
        ca_copy = deepcopy(ca_res)
        ca_copy["results"][0]["success"] = True

        diffs = conv_anal._compare_results_to_dict(ca_copy)

        assert diffs["success"] == [0]
        assert diffs["iters"] == []
        assert diffs["iters_in_restoration"] == []
        assert diffs["iters_w_regularization"] == []
        assert diffs["numerical_issues"] == []

    @pytest.mark.unit
    def test_compare_results_to_dict_iters(self, conv_anal):
        ca_copy = deepcopy(ca_res)
        ca_copy["results"][0]["results"]["iters"] = 8
        ca_copy["results"][1]["results"]["iters"] = 9

        diffs = conv_anal._compare_results_to_dict(ca_copy)

        assert diffs["success"] == []
        assert diffs["iters"] == [1]
        assert diffs["iters_in_restoration"] == []
        assert diffs["iters_w_regularization"] == []
        assert diffs["numerical_issues"] == []

        diffs = conv_anal._compare_results_to_dict(ca_copy, abs_tol=0, rel_tol=0)

        assert diffs["success"] == []
        assert diffs["iters"] == [0, 1]
        assert diffs["iters_in_restoration"] == []
        assert diffs["iters_w_regularization"] == []
        assert diffs["numerical_issues"] == []

    @pytest.mark.unit
    def test_compare_results_to_dict_restoration(self, conv_anal):
        ca_copy = deepcopy(ca_res)
        ca_copy["results"][0]["results"]["iters_in_restoration"] = 5
        ca_copy["results"][1]["results"]["iters_in_restoration"] = 6

        diffs = conv_anal._compare_results_to_dict(ca_copy)

        assert diffs["success"] == []
        assert diffs["iters"] == []
        assert diffs["iters_in_restoration"] == [1]
        assert diffs["iters_w_regularization"] == []
        assert diffs["numerical_issues"] == []

        diffs = conv_anal._compare_results_to_dict(ca_copy, abs_tol=0, rel_tol=0)

        assert diffs["success"] == []
        assert diffs["iters"] == []
        assert diffs["iters_in_restoration"] == [0, 1]
        assert diffs["iters_w_regularization"] == []
        assert diffs["numerical_issues"] == []

    @pytest.mark.unit
    def test_compare_results_to_dict_regularization(self, conv_anal):
        ca_copy = deepcopy(ca_res)
        ca_copy["results"][0]["results"]["iters_w_regularization"] = 1
        ca_copy["results"][1]["results"]["iters_w_regularization"] = 2

        diffs = conv_anal._compare_results_to_dict(ca_copy)

        assert diffs["success"] == []
        assert diffs["iters"] == []
        assert diffs["iters_in_restoration"] == []
        assert diffs["iters_w_regularization"] == [1]
        assert diffs["numerical_issues"] == []

        diffs = conv_anal._compare_results_to_dict(ca_copy, abs_tol=0, rel_tol=0)

        assert diffs["success"] == []
        assert diffs["iters"] == []
        assert diffs["iters_in_restoration"] == []
        assert diffs["iters_w_regularization"] == [0, 1]
        assert diffs["numerical_issues"] == []

    @pytest.mark.unit
    def test_compare_results_to_dict_numerical_issues(self, conv_anal):
        ca_copy = deepcopy(ca_res)
        ca_copy["results"][1]["results"]["numerical_issues"] = False

        diffs = conv_anal._compare_results_to_dict(ca_copy)

        assert diffs["success"] == []
        assert diffs["iters"] == []
        assert diffs["iters_in_restoration"] == []
        assert diffs["iters_w_regularization"] == []
        assert diffs["numerical_issues"] == [1]

    @pytest.mark.integration
    @pytest.mark.solver
    def test_compare_convergence_to_baseline(self, model):
        fname = os.path.join(currdir, "convergence_baseline.json")

        ca = IpoptConvergenceAnalysis(
            model=model,
        )

        diffs = ca.compare_convergence_to_baseline(fname)

        # Baseline has incorrect values
        assert diffs == {
            "success": [0],
            "iters": [],
            "iters_in_restoration": [],
            "iters_w_regularization": [1],
            "numerical_issues": [],
        }

    @pytest.mark.integration
    @pytest.mark.solver
    def test_assert_baseline_comparison(self, model):
        fname = os.path.join(currdir, "convergence_baseline.json")

        ca = IpoptConvergenceAnalysis(
            model=model,
        )

        # Baseline has incorrect values
        with pytest.raises(
            AssertionError,
            match="Convergence analysis does not match baseline",
        ):
            ca.assert_baseline_comparison(fname)


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
            "c: Potential evaluation error in x**(-2); base bounds are (0, inf); exponent bounds are (-2, -2)",
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
            "c: Potential evaluation error in x**(-2); base bounds are (-1, inf); exponent bounds are (-2, -2)",
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
            "c: Potential evaluation error in x**(-2.5); base bounds are (0, inf); exponent bounds are (-2.5, -2.5)",
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
            "c: Potential evaluation error in x**(-2.5); base bounds are (-1, inf); exponent bounds are (-2.5, -2.5)",
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


class TestCheckParallelJacobian:
    @pytest.mark.unit
    def test_invalid_direction(self):
        m = ConcreteModel()

        with pytest.raises(
            ValueError,
            match=re.escape(
                "Unrecognised value for direction (foo). " "Must be 'row' or 'column'."
            ),
        ):
            check_parallel_jacobian(m, direction="foo")

    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.v1 = Var(initialize=1e-8)
        m.v2 = Var()
        m.v3 = Var()
        m.v4 = Var()

        m.c1 = Constraint(expr=m.v1 == m.v2 - 0.99999 * m.v4)
        m.c2 = Constraint(expr=m.v1 + 1.00001 * m.v4 == 1e-8 * m.v3)
        m.c3 = Constraint(expr=1e8 * (m.v1 + m.v4) + 1e10 * m.v2 == 1e-6 * m.v3)
        m.c4 = Constraint(expr=-m.v1 == -0.99999 * (m.v2 - m.v4))

        return m

    @pytest.mark.unit
    def test_rows(self, model):
        assert check_parallel_jacobian(model, direction="row") == [(model.c1, model.c4)]
        assert check_parallel_jacobian(model, direction="row", tolerance=0) == []

    @pytest.mark.unit
    def test_columns(self, model):
        pcol = check_parallel_jacobian(model, direction="column")

        expected = [
            ("v1", "v2"),
            ("v1", "v4"),
            ("v2", "v4"),
        ]

        for i in pcol:
            assert tuple(sorted([i[0].name, i[1].name])) in expected


class TestCheckIllConditioning:
    @pytest.mark.unit
    def test_invalid_direction(self):
        m = ConcreteModel()

        with pytest.raises(
            ValueError,
            match=re.escape(
                "Unrecognised value for direction (foo). " "Must be 'row' or 'column'."
            ),
        ):
            compute_ill_conditioning_certificate(m, direction="foo")

    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()

        m.v1 = Var(initialize=1e-8)
        m.v2 = Var()
        m.v3 = Var()
        m.v4 = Var()

        m.c1 = Constraint(expr=m.v1 == m.v2 - 0.99999 * m.v4)
        m.c2 = Constraint(expr=m.v1 + 1.00001 * m.v4 == 1e-8 * m.v3)
        m.c3 = Constraint(expr=1e8 * (m.v1 + m.v4) + 1e10 * m.v2 == 1e-6 * m.v3)
        m.c4 = Constraint(expr=-m.v1 == -0.99999 * (m.v2 - m.v4))

        return m

    @pytest.fixture(scope="class")
    def afiro(self):
        # NETLIB AFIRO example
        m = ConcreteModel()

        # Vars
        m.X01 = Var(initialize=1)
        m.X02 = Var(initialize=1)
        m.X03 = Var(initialize=1)
        m.X04 = Var(initialize=1)
        m.X06 = Var(initialize=1)
        m.X07 = Var(initialize=1)
        m.X08 = Var(initialize=1)
        m.X09 = Var(initialize=1)
        m.X10 = Var(initialize=1)
        m.X11 = Var(initialize=1)
        m.X12 = Var(initialize=1)
        m.X13 = Var(initialize=1)
        m.X14 = Var(initialize=1)
        m.X15 = Var(initialize=1)
        m.X16 = Var(initialize=1)
        m.X22 = Var(initialize=1)
        m.X23 = Var(initialize=1)
        m.X24 = Var(initialize=1)
        m.X25 = Var(initialize=1)
        m.X26 = Var(initialize=1)
        m.X28 = Var(initialize=1)
        m.X29 = Var(initialize=1)
        m.X30 = Var(initialize=1)
        m.X31 = Var(initialize=1)
        m.X32 = Var(initialize=1)
        m.X33 = Var(initialize=1)
        m.X34 = Var(initialize=1)
        m.X35 = Var(initialize=1)
        m.X36 = Var(initialize=1)
        m.X37 = Var(initialize=1)
        m.X38 = Var(initialize=1)
        m.X39 = Var(initialize=1)

        # Constraints

        m.R09 = Constraint(expr=-m.X01 + m.X02 + m.X03 == 0)
        m.R10 = Constraint(expr=-1.06 * m.X01 + m.X04 == 0)
        m.X05 = Constraint(expr=m.X01 <= 80)
        m.X21 = Constraint(expr=-m.X02 + 1.4 * m.X14 <= 0)
        m.R12 = Constraint(expr=-m.X06 - m.X07 - m.X08 - m.X09 + m.X14 + m.X15 == 0)
        m.R13 = Constraint(
            expr=-1.06 * m.X06 - 1.06 * m.X07 - 0.96 * m.X08 - 0.86 * m.X09 + m.X16 == 0
        )
        m.X17 = Constraint(expr=m.X06 - m.X10 <= 80)
        m.X18 = Constraint(expr=m.X07 - m.X11 <= 0)
        m.X19 = Constraint(expr=m.X08 - m.X12 <= 0)
        m.X20 = Constraint(expr=m.X09 - m.X13 <= 0)
        m.R19 = Constraint(expr=-1 * m.X22 + 1 * m.X23 + 1 * m.X24 + 1 * m.X25 == 0)
        m.R20 = Constraint(expr=-0.43 * m.X22 + m.X26 == 0)
        m.X27 = Constraint(expr=m.X22 <= 500)
        m.X44 = Constraint(expr=-m.X23 + 1.4 * m.X36 <= 0)
        m.R22 = Constraint(
            expr=-0.43 * m.X28 - 0.43 * m.X29 - 0.39 * m.X30 - 0.37 * m.X31 + m.X38 == 0
        )
        m.R23 = Constraint(
            expr=1 * m.X28
            + 1 * m.X29
            + 1 * m.X30
            + 1 * m.X31
            - 1 * m.X36
            + 1 * m.X37
            + 1 * m.X39
            == 44
        )
        m.X40 = Constraint(expr=m.X28 - m.X32 <= 500)
        m.X41 = Constraint(expr=m.X29 - m.X33 <= 0)
        m.X42 = Constraint(expr=m.X30 - m.X34 <= 0)
        m.X43 = Constraint(expr=m.X31 - m.X35 <= 0)
        m.X45 = Constraint(
            expr=2.364 * m.X10
            + 2.386 * m.X11
            + 2.408 * m.X12
            + 2.429 * m.X13
            - m.X25
            + 2.191 * m.X32
            + 2.219 * m.X33
            + 2.249 * m.X34
            + 2.279 * m.X35
            <= 0
        )
        m.X46 = Constraint(expr=-m.X03 + 0.109 * m.X22 <= 0)
        m.X47 = Constraint(
            expr=-m.X15 + 0.109 * m.X28 + 0.108 * m.X29 + 0.108 * m.X30 + 0.107 * m.X31
            <= 0
        )
        m.X48 = Constraint(expr=0.301 * m.X01 - m.X24 <= 0)
        m.X49 = Constraint(
            expr=0.301 * m.X06 + 0.313 * m.X07 + 0.313 * m.X08 + 0.326 * m.X09 - m.X37
            <= 0
        )
        m.X50 = Constraint(expr=m.X04 + m.X26 <= 310)
        m.X51 = Constraint(expr=m.X16 + m.X38 <= 300)

        # Degenerate constraint
        m.R09b = Constraint(expr=m.X01 - 0.999999999 * m.X02 - m.X03 == 0)

        return m

    @pytest.mark.unit
    def test_beta(self, model, caplog):

        compute_ill_conditioning_certificate(model)

        expected = (
            "Ill conditioning checks are a beta capability. Please be aware that "
            "the name, location, and API for this may change in future releases."
        )

        assert expected in caplog.text

    @pytest.mark.component
    @pytest.mark.solver
    def test_rows(self, model):
        assert compute_ill_conditioning_certificate(model, direction="row") == [
            (model.c4, pytest.approx(0.50000002, rel=1e-5)),
            (model.c4, pytest.approx(0.49999998, rel=1e-5)),
        ]

    @pytest.mark.component
    @pytest.mark.solver
    def test_rows(self, afiro):
        assert compute_ill_conditioning_certificate(afiro, direction="row") == [
            (afiro.R09, pytest.approx(0.5, rel=1e-5)),
            (afiro.R09b, pytest.approx(0.5, rel=1e-5)),
        ]

    @pytest.mark.component
    @pytest.mark.solver
    def test_columns(self, afiro):
        assert compute_ill_conditioning_certificate(afiro, direction="column") == [
            (afiro.X39, pytest.approx(1.1955465, rel=1e-5)),
            (afiro.X23, pytest.approx(1.0668697, rel=1e-5)),
            (afiro.X25, pytest.approx(-1.0668697, rel=1e-5)),
            (afiro.X09, pytest.approx(-0.95897123, rel=1e-5)),
            (afiro.X13, pytest.approx(-0.95897123, rel=1e-5)),
            (afiro.X06, pytest.approx(0.91651956, rel=1e-5)),
            (afiro.X10, pytest.approx(0.91651956, rel=1e-5)),
            (afiro.X36, pytest.approx(0.76204977, rel=1e-5)),
            (afiro.X31, pytest.approx(-0.39674454, rel=1e-5)),
            (afiro.X35, pytest.approx(-0.39674454, rel=1e-5)),
            (afiro.X16, pytest.approx(0.14679548, rel=1e-5)),
            (afiro.X38, pytest.approx(-0.14679548, rel=1e-5)),
            (afiro.X15, pytest.approx(-0.042451666, rel=1e-5)),
            (afiro.X37, pytest.approx(-0.036752232, rel=1e-5)),
        ]


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


class TestConstraintTermAnalysisVisitor:
    @pytest.mark.unit
    def test_sum_combinations(self):
        # Check method to generate sums of all combinations of terms
        # excludes single term sums
        terms = [1, 2, 3, 4, 5]
        visitor = ConstraintTermAnalysisVisitor(max_canceling_terms=None)
        sums = [i for i in visitor._generate_combinations(terms)]

        print(sums)

        expected = [
            ((0, 1), (1, 2)),
            ((0, 1), (2, 3)),
            ((0, 1), (3, 4)),
            ((0, 1), (4, 5)),
            ((1, 2), (2, 3)),
            ((1, 2), (3, 4)),
            ((1, 2), (4, 5)),
            ((2, 3), (3, 4)),
            ((2, 3), (4, 5)),
            ((3, 4), (4, 5)),
            ((0, 1), (1, 2), (2, 3)),
            ((0, 1), (1, 2), (3, 4)),
            ((0, 1), (1, 2), (4, 5)),
            ((0, 1), (2, 3), (3, 4)),
            ((0, 1), (2, 3), (4, 5)),
            ((0, 1), (3, 4), (4, 5)),
            ((1, 2), (2, 3), (3, 4)),
            ((1, 2), (2, 3), (4, 5)),
            ((1, 2), (3, 4), (4, 5)),
            ((2, 3), (3, 4), (4, 5)),
            ((0, 1), (1, 2), (2, 3), (3, 4)),
            ((0, 1), (1, 2), (2, 3), (4, 5)),
            ((0, 1), (1, 2), (3, 4), (4, 5)),
            ((0, 1), (2, 3), (3, 4), (4, 5)),
            ((1, 2), (2, 3), (3, 4), (4, 5)),
            ((0, 1), (1, 2), (2, 3), (3, 4), (4, 5)),
        ]

        assert sums == expected

    @pytest.mark.unit
    def test_sum_combinations_limited_default(self):
        # Check method to generate sums of all combinations of terms
        # excludes single term sums
        terms = [1, 2, 3, 4, 5]
        visitor = ConstraintTermAnalysisVisitor()
        sums = [i for i in visitor._generate_combinations(terms)]

        expected = expected = [
            ((0, 1), (1, 2)),
            ((0, 1), (2, 3)),
            ((0, 1), (3, 4)),
            ((0, 1), (4, 5)),
            ((1, 2), (2, 3)),
            ((1, 2), (3, 4)),
            ((1, 2), (4, 5)),
            ((2, 3), (3, 4)),
            ((2, 3), (4, 5)),
            ((3, 4), (4, 5)),
            ((0, 1), (1, 2), (2, 3)),
            ((0, 1), (1, 2), (3, 4)),
            ((0, 1), (1, 2), (4, 5)),
            ((0, 1), (2, 3), (3, 4)),
            ((0, 1), (2, 3), (4, 5)),
            ((0, 1), (3, 4), (4, 5)),
            ((1, 2), (2, 3), (3, 4)),
            ((1, 2), (2, 3), (4, 5)),
            ((1, 2), (3, 4), (4, 5)),
            ((2, 3), (3, 4), (4, 5)),
            ((0, 1), (1, 2), (2, 3), (3, 4)),
            ((0, 1), (1, 2), (2, 3), (4, 5)),
            ((0, 1), (1, 2), (3, 4), (4, 5)),
            ((0, 1), (2, 3), (3, 4), (4, 5)),
            ((1, 2), (2, 3), (3, 4), (4, 5)),
            # 5 term combination should be excluded
        ]

        assert sums == expected

    @pytest.mark.unit
    def test_check_sum_cancellations(self):
        terms = [1, -2, 3, -4, 5]
        visitor = ConstraintTermAnalysisVisitor()
        cancellations = visitor._check_sum_cancellations(terms)

        # We expect to canceling combinations
        # Results should be a list with each entry being a tuple containing a
        # canceling combination
        # In turn, each element of the tuple should be a 2-tuple with the
        # position of the term in the input list and its value
        expected = [((0, 1), (2, 3), (3, -4)), ((0, 1), (1, -2), (3, -4), (4, 5))]

        assert cancellations == expected

    @pytest.mark.unit
    def test_int(self):
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=7)

        assert vv == [7]
        assert len(mm) == 0
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.mark.unit
    def test_float(self):
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=7.7)

        assert vv == [7.7]
        assert len(mm) == 0
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.mark.unit
    def test_scalar_param(self):
        m = ConcreteModel()
        m.scalar_param = Param(initialize=1, mutable=True)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=m.scalar_param
        )

        assert vv == [1]
        assert len(mm) == 0
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.mark.unit
    def test_indexed_param(self):
        m = ConcreteModel()
        m.indexed_param = Param(["a", "b"], initialize=1, mutable=True)

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=m.indexed_param["a"]
        )
        assert vv == [1]
        assert len(mm) == 0
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.mark.unit
    def test_scalar_var(self):
        m = ConcreteModel()
        m.scalar_var = Var(initialize=1)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=m.scalar_var
        )

        assert vv == [1]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_indexed_var(self):
        m = ConcreteModel()
        m.indexed_var = Var(["a", "b"], initialize=1)

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=m.indexed_var["a"]
        )
        assert vv == [1]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_scalar_var_fixed(self):
        m = ConcreteModel()
        m.scalar_var = Var(initialize=1)
        m.scalar_var.fix()
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=m.scalar_var
        )

        assert vv == [1]
        assert len(mm) == 0
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.mark.unit
    def test_indexed_var_fixed(self):
        m = ConcreteModel()
        m.indexed_var = Var(["a", "b"], initialize=1)
        m.indexed_var.fix()

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=m.indexed_var["a"]
        )
        assert vv == [1]
        assert len(mm) == 0
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.mark.unit
    def test_equality_expr(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1e-7)
        m.v2 = Var(initialize=1e7)

        expr = m.v1 == m.v2
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1e-7, 1e7]
        assert expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_equality_expr_constant(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1e-7)
        m.v2 = Var(initialize=1e7)

        # Fix v1, not constant yet as v2 still free
        m.v1.fix()

        expr = m.v1 == m.v2
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1e-7, 1e7]
        assert expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert not k
        assert not tr

        # Fix v2, now constant
        m.v2.fix()

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1e-7, 1e7]
        assert expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.mark.unit
    def test_sum_expr(self):
        m = ConcreteModel()
        m.v1 = Var(["a", "b", "c"], initialize=1e7)
        m.v1["a"].set_value(1e-7)

        expr = sum(m.v1[i] for i in m.v1)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [1e-7, 1e7, 1e7]
        assert expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_sum_expr_constant(self):
        m = ConcreteModel()
        m.v1 = Var(["a", "b", "c"], initialize=1e7)
        m.v1["a"].set_value(1e-7)
        m.v1.fix()

        expr = sum(m.v1[i] for i in m.v1)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [1e-7, 1e7, 1e7]
        assert expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.fixture(scope="class")
    def model(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=2)
        m.v2 = Var(initialize=3)
        m.v3 = Var(initialize=5.0000001)
        m.v4 = Var(initialize=5)

        return m

    @pytest.mark.unit
    def test_product_expr(self, model):
        m = ConcreteModel()
        expr = model.v1 * model.v2 * model.v3
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [pytest.approx(30.0000006, rel=1e-8)]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_product_sum_expr(self, model):
        expr = (model.v1 + model.v2) * (model.v3 + model.v4)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [pytest.approx((2 + 3) * (5.0000001 + 5), rel=1e-8)]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_product_sum_expr_w_negation(self, model):
        expr = (model.v1 + model.v2) * (model.v3 - model.v4)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [pytest.approx(0.0000005, rel=1e-8)]
        assert len(mm) == 0
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_division_expr(self, model):
        expr = model.v1 / model.v2 / model.v3
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [pytest.approx(2 / 3 / 5.0000001, rel=1e-8)]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_division_sum_expr(self, model):
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=(model.v1 + model.v2) / (model.v3 + model.v4)
        )

        assert vv == [
            pytest.approx(2 / (5.0000001 + 5), rel=1e-8),
            pytest.approx(3 / (5.0000001 + 5), rel=1e-8),
        ]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_division_sum_expr_w_negation(self, model):
        expr = (model.v1 - model.v2) / (model.v3 - model.v4)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [
            pytest.approx(2 / 0.0000001, rel=1e-8),
            pytest.approx(-3 / 0.0000001, rel=1e-8),
        ]
        assert len(mm) == 0
        # Check for cancellation should be deferred
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_division_sum_expr_w_zero_denominator(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=2)
        m.v2 = Var(initialize=3)
        m.v3 = Var(initialize=5)
        m.v4 = Var(initialize=5)

        expr = (m.v1 - m.v2) / (m.v3 - m.v4)
        with pytest.raises(
            ZeroDivisionError,
            match=re.escape(
                "Error in ConstraintTermAnalysisVisitor: found division with denominator of 0 "
                "((v1 - v2)/(v3 - v4))."
            ),
        ):
            ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

    @pytest.mark.unit
    def test_division_sum_expr_w_negation_deferred(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=2)
        m.v2 = Var(initialize=2)
        m.v3 = Var(initialize=5.0000001)
        m.v4 = Var(initialize=5)

        expr = ((m.v1 - m.v2) / (m.v3 - m.v4)) ** 2
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [pytest.approx(0, abs=1e-8)]
        assert len(mm) == 0
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_division_expr_error(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=2)
        m.v2 = Var(initialize=0)

        with pytest.raises(
            ZeroDivisionError,
            match=re.escape(
                "Error in ConstraintTermAnalysisVisitor: found division with "
                "denominator of 0 (v1/v2)."
            ),
        ):
            ConstraintTermAnalysisVisitor().walk_expression(expr=m.v1 / m.v2)

    @pytest.mark.unit
    def test_pow_expr(self, model):
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=model.v1**model.v2
        )

        assert vv == [8]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_pow_sum_expr(self, model):
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=(model.v1 + model.v2) ** (model.v3 + model.v4)
        )

        assert vv == [pytest.approx(5**10.0000001, rel=1e-8)]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_pow_sum_expr_w_negation(self, model):
        expr = (model.v1 + model.v2) ** (model.v3 - model.v4)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [pytest.approx((2 + 3) ** (0.0000001), rel=1e-8)]
        assert len(mm) == 0
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    @pytest.fixture(scope="class")
    def func_model(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1)
        m.v2 = Var(initialize=0.99999)
        m.v3 = Var(initialize=-100)

        m.v2.fix()

        return m

    @pytest.mark.unit
    def test_negation_expr(self, func_model):
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=-func_model.v1
        )

        assert vv == [-1]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_negation_sum_expr(self, func_model):
        expr = -(func_model.v1 - func_model.v2)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        # Checking the cancellation should be deferred
        # Expect to get two values back
        assert vv == [pytest.approx(-1, rel=1e-8), pytest.approx(0.99999, rel=1e-8)]
        assert len(mm) == 0
        # Check is deferred, so no cancellations identified
        assert len(cc) == 0
        assert not k
        assert not tr

    # acosh has bounds that don't fit with other functions - we will assume we caught enough
    func_list = [
        exp,
        log,
        log10,
        sqrt,
        sin,
        cos,
        tan,
        asin,
        acos,
        atan,
        sinh,
        cosh,
        tanh,
        asinh,
        atanh,
    ]
    func_error_list = [log, log10, sqrt, asin, acos, acosh, atanh]

    @pytest.mark.unit
    @pytest.mark.parametrize("func", func_list)
    def test_func_expr(self, func_model, func):
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=func(func_model.v2)
        )

        assert vv == [pytest.approx(value(func(0.99999)), rel=1e-8)]
        assert len(mm) == 0
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.mark.unit
    @pytest.mark.parametrize("func", func_list)
    def test_func_sum_expr(self, func_model, func):
        expr = func(func_model.v1 - func_model.v2)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [pytest.approx(value(func(0.00001)), rel=1e-8)]
        assert len(mm) == 0
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    @pytest.mark.unit
    @pytest.mark.parametrize("func", func_list)
    def test_func_sum_expr_w_negation(self, func_model, func):
        expr = func(func_model.v1 - func_model.v2)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [pytest.approx(value(func(0.00001)), rel=1e-8)]
        assert len(mm) == 0
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    @pytest.mark.unit
    @pytest.mark.parametrize("func", func_error_list)
    def test_func_expr_error(self, func_model, func):
        with pytest.raises(
            ValueError,
            match=re.escape(
                "Error in ConstraintTermAnalysisVisitor: error evaluating "
            ),
        ):
            ConstraintTermAnalysisVisitor().walk_expression(expr=func(func_model.v3))

    @pytest.mark.unit
    def test_expr_if(self, model):
        model.exprif = Expr_if(
            IF=model.v1,
            THEN=model.v1 + model.v2,
            ELSE=model.v3 - model.v4,
        )

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=model.exprif
        )

        assert vv == [pytest.approx(5, rel=1e-8)]
        assert len(mm) == 0
        assert model.exprif in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    @pytest.mark.unit
    @pytest.mark.skipif(
        not AmplInterface.available(), reason="pynumero_ASL is not available"
    )
    @pytest.mark.skipif(not cubic_roots_available, reason="Cubic roots not available")
    def test_ext_func(self):
        # Use the cubic root external function to test
        m = ConcreteModel()
        m.a = Var(initialize=1)
        m.b = Var(initialize=1)

        m.expr_write = CubicThermoExpressions(m)
        Z = m.expr_write.z_liq(eos=CubicEoS.PR, A=m.a, B=m.b)

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=Z)

        assert vv == [pytest.approx(-2.1149075414767577, rel=1e-8)]
        assert len(mm) == 0
        assert Z in cc
        assert len(cc) == 1
        assert not k
        assert not tr

        # Check that model state did not change
        assert value(m.a) == 1
        assert value(m.b) == 1
        assert value(Z) == pytest.approx(-2.1149075414767577, rel=1e-8)

    @pytest.mark.unit
    def test_equality_sum_expr(self):
        m = ConcreteModel()
        m.v1 = Var(["a", "b", "c"], initialize=1e7)
        m.v2 = Var(initialize=1e-7)

        expr = m.v2 == sum(m.v1[i] for i in m.v1)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1e-7, 3e7]
        assert expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_inequality_sum_expr(self):
        m = ConcreteModel()
        m.v1 = Var(["a", "b", "c"], initialize=1e7)
        m.v2 = Var(initialize=1e-7)

        expr = m.v2 <= sum(m.v1[i] for i in m.v1)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1e-7, 3e7]
        assert expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_compound_equality_expr_1(self):
        m = ConcreteModel()
        m.v1 = Var(["a", "b", "c"], initialize=1e7)
        m.v2 = Var(initialize=1e-7)

        expr = 6 * m.v2 == 8 * sum(m.v1[i] for i in m.v1)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-6e-7, 2.4e8]
        assert expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_ranged_expr(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1e7)
        m.v2 = Var(initialize=1e-7)
        m.v3 = Var(initialize=1e7)

        m.expr = EXPR.RangedExpression(args=(m.v1, m.v2, m.v3), strict=True)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=m.expr)

        assert vv == [-1e7, 1e-7, -1e7]
        assert m.expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert not k
        assert not tr

        # Fix v1 and v2 to make first two terms constant
        m.v1.fix()
        m.v2.fix()
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=m.expr)

        # Should not be flagged as constant due to v3
        assert vv == [-1e7, 1e-7, -1e7]
        assert m.expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert not k
        assert not tr

        # Fix v3 to make all terms constant
        m.v3.fix()
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=m.expr)

        # Should now be constant
        assert vv == [-1e7, 1e-7, -1e7]
        assert m.expr in mm
        assert len(mm) == 1
        assert len(cc) == 0
        assert k
        assert not tr

    @pytest.mark.unit
    def test_compound_equality_expr_2(self):
        m = ConcreteModel()
        m.v1 = Var(["a", "b", "c"], initialize=1e7)
        m.v2 = Var(initialize=1e-7)
        m.v3 = Var(initialize=1e3)

        # Set this small so we get two mismatched warnings
        m.v1["a"].set_value(1e-7)

        expr1 = sum(m.v1[i] for i in m.v1)
        expr = 6 * m.v2 == 8 * expr1 + m.v3
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-6e-7, pytest.approx(8 * (2e7 + 1e-7) + 1000, rel=1e-8)]
        assert expr in mm
        assert expr1 in mm
        assert len(mm) == 2
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_canceling_sum_expr(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=2)
        m.v2 = Var(initialize=2)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=m.v1 - m.v2
        )

        assert vv == [2, -2]
        assert len(mm) == 0
        # We do not check cancellation at the sum, so this should be empty
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_expr_w_canceling_sum(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=2)
        m.v2 = Var(initialize=2)
        m.v3 = Var(initialize=3)

        expr = m.v3 * (m.v1 - m.v2)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [6, -6]
        assert len(mm) == 0
        # Check for cancellation should be deferred
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_expr_w_deferred_canceling_sum(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=2)
        m.v2 = Var(initialize=2)
        m.v3 = Var(initialize=3)

        expr = (m.v3 * (m.v1 - m.v2)) ** 2
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [0]
        assert len(mm) == 0
        # We should get a warning about canceling sums here
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

        # Check for tolerance of sum cancellation
        m.v2.set_value(2.00022)

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor(
            term_cancellation_tolerance=1e-4
        ).walk_expression(expr=expr)

        assert vv == [pytest.approx((3 * -0.00022) ** 2, rel=1e-8)]
        assert len(mm) == 0
        # This should pass as the difference is greater than tol
        assert len(cc) == 0
        assert not k
        assert not tr

        # Change tolerance so it should identify cancellation
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor(
            term_cancellation_tolerance=1e-3
        ).walk_expression(expr=expr)

        assert vv == [pytest.approx((3 * -0.00022) ** 2, rel=1e-8)]
        assert len(mm) == 0
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_canceling_equality_expr_safe(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1e7)
        m.v2 = Var(initialize=1e7)

        # This is a standard constraint form, so should have no issues despite cancellation
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=0 == m.v1 - m.v2
        )

        assert vv == [0, pytest.approx(0, abs=1e-12)]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_canceling_equality_expr_zero_term(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1e7)
        m.v2 = Var(initialize=1e7)
        m.v3 = Var(initialize=0)

        expr = m.v3 == m.v1 - m.v2
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [0, pytest.approx(0, abs=1e-12)]
        assert len(mm) == 0
        # No canceling terms as v3=0 is ignored thus we have a=b form
        assert len(cc) == 0
        assert not k
        assert not tr

        # Set v3 above zero tolerance
        m.v3.set_value(1e-4)

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)
        assert vv == [-1e-4, pytest.approx(0, abs=1e-12)]
        assert len(mm) == 0
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_canceling_equality_expr_compound(self):
        m = ConcreteModel()
        m.v1 = Var(["a", "b"], initialize=5e6)
        m.v2 = Var(initialize=1e7)
        m.v3 = Var(initialize=0)

        expr = m.v3 == sum(v for v in m.v1.values()) - m.v2
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [0, pytest.approx(0, abs=1e-12)]
        assert len(mm) == 0
        # No canceling terms as v3=0 is ignored thus we have a=b form
        assert len(cc) == 0
        assert not k
        assert not tr

        # Set v3 above zero tolerance
        m.v3.set_value(1e-4)

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)
        assert vv == [-1e-4, pytest.approx(0, abs=1e-12)]
        assert len(mm) == 0
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    # Double check for a+eps=c form gets flagged in some way
    @pytest.mark.unit
    def test_canceling_equality_expr_canceling_sides(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1)
        m.v2 = Var(initialize=1)
        m.v3 = Var(initialize=1e-8)

        expr1 = m.v1 + m.v3
        expr = m.v2 == expr1
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1, pytest.approx(1 + 1e-8, abs=1e-8)]
        assert expr1 in mm
        assert len(mm) == 1
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    # Check to make sure simple linking constraints are not flagged as canceling
    @pytest.mark.unit
    def test_linking_equality_expr(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1)
        m.v2 = Var(initialize=1)

        expr = m.v1 == m.v2
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1, 1]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_linking_equality_expr_compound(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1)
        m.v2 = Var(initialize=1)
        m.v3 = Var(initialize=1)

        expr = m.v1 == m.v2 * m.v3
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1, 1]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.component
    def test_external_function_w_string_argument(self):
        m = ConcreteModel()
        m.properties = iapws95.Iapws95ParameterBlock()
        m.state = m.properties.build_state_block([0])

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=m.state[0].enth_mol
        )

        assert vv == [pytest.approx(1.1021387e-2, rel=1e-6)]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

        # Test nested external functions
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(
            expr=m.state[0].temperature
        )

        assert vv == [pytest.approx(270.4877, rel=1e-6)]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_zero_tolerance(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1e-12)
        m.v2 = Var(initialize=1)
        m.v3 = Var(initialize=1e-12)

        expr1 = m.v1 - m.v3
        expr2 = expr1 + 1
        expr = m.v2 == expr2
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1, pytest.approx(1, abs=1e-8)]
        # We expect no mismatches, as smallest terms are below zero tolerance
        assert len(mm) == 0
        # We expect no canceling terms, as v1 and v3 are ignored (zero tolerance), leaving
        # 1 == 1
        assert len(cc) == 0
        assert not k
        assert not tr

        # Set v1 and v3 above zero tolerance
        m.v1.set_value(1e-8)
        m.v3.set_value(1e-8)

        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [-1, pytest.approx(1, abs=1e-8)]
        assert expr2 in mm
        assert len(mm) == 1
        assert expr in cc
        assert len(cc) == 1
        assert not k
        assert not tr

    # Test to make sure scaled constraints are not flagged as issues
    @pytest.mark.unit
    def test_scaled_equality_expr_product(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1)
        m.v2 = Var(initialize=1)
        m.v3 = Var(initialize=2)

        expr = 0 == m.v3 * (m.v1 - m.v2)
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [0, 0]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_scaled_equality_expr_division(self):
        m = ConcreteModel()
        m.v1 = Var(initialize=1)
        m.v2 = Var(initialize=1)
        m.v3 = Var(initialize=2)

        expr = 0 == (m.v1 - m.v2) / m.v3
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [0, 0]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_mole_fraction_constraint(self):
        m = ConcreteModel()
        m.mole_frac_a = Var(initialize=0.5)
        m.flow_a = Var(initialize=100)
        m.flow_b = Var(initialize=100)

        expr = m.mole_frac_a * (m.flow_a + m.flow_b) == m.flow_a
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [pytest.approx(-100, rel=1e-5), pytest.approx(100, rel=1e-5)]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr

    @pytest.mark.unit
    def test_mole_fraction_constraint_trace(self):
        m = ConcreteModel()
        m.mole_frac_a = Var(initialize=0.999810)
        m.flow_a = Var(initialize=122.746)
        m.flow_b = Var(initialize=0.0233239)

        expr = m.mole_frac_a * (m.flow_a + m.flow_b) == m.flow_a
        vv, mm, cc, k, tr = ConstraintTermAnalysisVisitor().walk_expression(expr=expr)

        assert vv == [
            pytest.approx(-122.746, rel=1e-5),
            pytest.approx(122.746, rel=1e-5),
        ]
        assert len(mm) == 0
        assert len(cc) == 0
        assert not k
        assert not tr
