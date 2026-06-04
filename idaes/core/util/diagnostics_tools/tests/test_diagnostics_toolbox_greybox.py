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
This module contains tests for the Diagnostics Toolbox using ExternalGreyBox models.
"""

from io import StringIO
import logging
import re

import pytest

import numpy as np
import pytest
from scipy.sparse import coo_matrix

from pyomo.contrib.pynumero.interfaces.external_grey_box import (
    ExternalGreyBoxBlock,
    ExternalGreyBoxModel,
)
from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    SolverFactory,
    units,
    Var,
)

from idaes.core.util.diagnostics_tools.diagnostics_toolbox import (
    DiagnosticsToolbox,
)

logging.getLogger("cyipopt").setLevel(logging.WARNING)


# TODO: Pyomo NLP does not know how to handle Grey Box components
# TODO: Pyomo MIS does not handle Grey Box components


class OutputGrayBox(ExternalGreyBoxModel):
    def input_names(self):
        return ["v2", "v5", "v6"]

    def output_names(self):
        return ["v1", "v3", "v4", "v7"]

    def set_input_values(self, input_values):
        self._input_values = list(input_values)

    def finalize_block_construction(self, pyomo_block):
        pyomo_block.outputs["v3"].setlb(0)  # v3
        pyomo_block.outputs["v3"].setub(5)  # v3
        pyomo_block.inputs["v5"].setlb(0)  # v5
        pyomo_block.inputs["v5"].setub(5)  # v5

    def evaluate_outputs(self):
        v1 = 10 - self._input_values[0]  # c1
        # Linear combination of c2 and c3 to eliminate v3 and solve for v4
        v4 = -2 * self._input_values[1] - 2 * self._input_values[1]
        v3 = v4 + self._input_values[1]  # c2
        v7 = 2e-8 * v1  # c4
        return [v1, v3, v4, v7]

    def evaluate_jacobian_outputs(self):
        row = np.zeros(6)  # output index
        col = np.zeros(6)  # input index
        data = np.zeros(6)  # jacobian values

        # dv1/dv2
        row[0], col[0], data[0] = (0, 0, -1)
        # dv4/dv5
        row[1], col[1], data[1] = (2, 1, -2)
        # dv4/dv6
        row[2], col[2], data[2] = (2, 2, -1)
        # dv3/dv5
        row[3], col[3], data[3] = (1, 1, -1)
        # dv3/dv6
        row[4], col[4], data[4] = (1, 2, -1)
        # dv7/dv2
        row[5], col[5], data[5] = (3, 0, -2e-8)

        return coo_matrix((data, (row, col)), shape=(4, 3))


class ResidualGrayBox(ExternalGreyBoxModel):
    def input_names(self):
        return ["v1", "v2", "v3", "v4", "v5", "v6", "v7"]

    def output_names(self):
        return []

    def equality_constraint_names(self):
        return ["c1", "c2", "c3", "c4"]

    def set_input_values(self, input_values):
        self._input_values = list(input_values)

    def finalize_block_construction(self, pyomo_block):
        pyomo_block.inputs["v3"].setlb(0)  # v3
        pyomo_block.inputs["v3"].setub(5)  # v3
        pyomo_block.inputs["v5"].setlb(0)  # v5
        pyomo_block.inputs["v5"].setub(5)  # v5

    def evaluate_outputs(self):
        raise NotImplementedError("This grey box only provides equality constraints")

    def evaluate_jacobian_outputs(self):
        raise NotImplementedError("This grey box only provides equality constraints")

    def evaluate_equality_constraints(self):
        # c1: m.v1 + m.b.v2 == 10
        # c2: m.b.v3 == m.b.v4 + m.b.v5
        # c3: 2 * m.b.v3 == 3 * m.b.v4 + 4 * m.b.v5 + m.b.v6
        # c4: m.b.v7 == 2e-8 * m.v1
        r = np.zeros(4)
        r[0] = self._input_values[0] + self._input_values[1] - 10
        r[1] = self._input_values[2] - self._input_values[3] - self._input_values[4]
        r[2] = (
            2 * self._input_values[2]
            - 3 * self._input_values[3]
            - 4 * self._input_values[4]
            - self._input_values[5]
        )
        r[3] = self._input_values[6] - 2e-8 * self._input_values[0]
        return r

    def evaluate_jacobian_equality_constraints(self):
        row = np.zeros(11)  # constraint index
        col = np.zeros(11)  # input index
        data = np.zeros(11)  # jacobian values

        # dc1/dv1
        row[0], col[0], data[0] = (0, 0, 1)
        # dc1/dv2
        row[1], col[1], data[1] = (0, 1, 1)
        # dc2/dv3
        row[2], col[2], data[2] = (1, 2, 1)
        # dc2/dv4
        row[3], col[3], data[3] = (1, 3, -1)
        # dc2/dv5
        row[4], col[4], data[4] = (1, 4, -1)
        # dc3/dv3
        row[5], col[5], data[5] = (2, 2, 2)
        # dc3/dv4
        row[6], col[6], data[6] = (2, 3, -3)
        # dc3/dv5
        row[7], col[7], data[7] = (2, 4, -4)
        # dc3/dv6
        row[8], col[8], data[8] = (2, 5, -1)
        # dc4/dv1
        row[9], col[9], data[9] = (3, 0, -2e-8)
        # dc4/dv7
        row[10], col[10], data[10] = (3, 6, 1)

        return coo_matrix((data, (row, col)), shape=(4, 7))


@pytest.fixture()
def model():
    # Note: This model is deliberately infeasible and poorly scaled to
    # test the diagnostics toolbox capabilities.
    m = ConcreteModel()
    m.b = Block()

    m.v1_1 = Var(units=units.m)  # External variable
    m.b.v2_1 = Var(units=units.m)
    m.b.v3_1 = Var(bounds=(0, 5))
    m.b.v4_1 = Var()
    m.b.v5_1 = Var(bounds=(0, 5))
    m.b.v6_1 = Var()
    m.b.v7_1 = Var(
        units=units.m, bounds=(0, 1)
    )  # Poorly scaled variable with lower bound

    m.v1_2 = Var(units=units.m)  # External variable
    m.b.v2_2 = Var(units=units.m)
    m.b.v3_2 = Var(bounds=(0, 5))
    m.b.v4_2 = Var()
    m.b.v5_2 = Var(bounds=(0, 5))
    m.b.v6_2 = Var()
    m.b.v7_2 = Var(
        units=units.m, bounds=(0, 1)
    )  # Poorly scaled variable with lower bound

    m.b.v8 = Var()  # Unused variable

    m.b.v2_1.fix(5)
    m.b.v5_1.fix(2)
    m.b.v6_1.fix(0)
    m.b.v2_2.fix(5)
    m.b.v5_2.fix(2)
    m.b.v6_2.fix(0)

    m.b.gb1 = ExternalGreyBoxBlock(
        external_model=OutputGrayBox(), build_implicit_constraint_objects=True
    )
    m.b.gb1_link_v1 = Constraint(expr=m.b.gb1.outputs["v1"] == m.v1_1)
    m.b.gb1_link_v2 = Constraint(expr=m.b.gb1.inputs["v2"] == m.b.v2_1)
    m.b.gb1_link_v3 = Constraint(expr=m.b.gb1.outputs["v3"] == m.b.v3_1)
    m.b.gb1_link_v4 = Constraint(expr=m.b.gb1.outputs["v4"] == m.b.v4_1)
    m.b.gb1_link_v5 = Constraint(expr=m.b.gb1.inputs["v5"] == m.b.v5_1)
    m.b.gb1_link_v6 = Constraint(expr=m.b.gb1.inputs["v6"] == m.b.v6_1)
    m.b.gb1_link_v7 = Constraint(expr=m.b.gb1.outputs["v7"] == m.b.v7_1)

    m.b.gb2 = ExternalGreyBoxBlock(
        external_model=ResidualGrayBox(),
        build_implicit_constraint_objects=True,
    )
    m.b.gb2_link_v1 = Constraint(expr=m.b.gb2.inputs["v1"] == m.v1_2)
    m.b.gb2_link_v2 = Constraint(expr=m.b.gb2.inputs["v2"] == m.b.v2_2)
    m.b.gb2_link_v3 = Constraint(expr=m.b.gb2.inputs["v3"] == m.b.v3_2)
    m.b.gb2_link_v4 = Constraint(expr=m.b.gb2.inputs["v4"] == m.b.v4_2)
    m.b.gb2_link_v5 = Constraint(expr=m.b.gb2.inputs["v5"] == m.b.v5_2)
    m.b.gb2_link_v6 = Constraint(expr=m.b.gb2.inputs["v6"] == m.b.v6_2)
    m.b.gb2_link_v7 = Constraint(expr=m.b.gb2.inputs["v7"] == m.b.v7_2)

    # Solve model
    solver = SolverFactory("cyipopt")
    res = solver.solve(m)

    return m


@pytest.fixture()
def diagnostics_toolbox(model):
    return DiagnosticsToolbox(model.b, include_grey_box_blocks=True)


@pytest.mark.unit
def test_display_external_variables(diagnostics_toolbox):
    stream = StringIO()

    diagnostics_toolbox.display_external_variables(stream=stream)

    expected = """====================================================================================
The following external variable(s) appear in constraints within the model:

    v1_1
    v1_2

====================================================================================
"""
    assert stream.getvalue() == expected


@pytest.mark.unit
def test_display_unused_variables(diagnostics_toolbox):
    stream = StringIO()

    diagnostics_toolbox.display_unused_variables(stream=stream)

    expected = """====================================================================================
The following variable(s) do not appear in any activated constraints within the model:

    b.v8

====================================================================================
"""

    assert stream.getvalue() == expected


@pytest.mark.unit
def test_display_variables_fixed_to_zero(diagnostics_toolbox):
    stream = StringIO()

    diagnostics_toolbox._model.gb1.outputs["v3"].fix(0)  # Fix v3 to zero for testing
    diagnostics_toolbox._model.gb2.inputs["v3"].fix(0)  # Fix v3 to zero for testing

    diagnostics_toolbox.display_variables_fixed_to_zero(stream=stream)

    expected = """====================================================================================
The following variable(s) are fixed to zero:

    b.v6_1
    b.v6_2
    b.gb1.outputs[v3]
    b.gb2.inputs[v3]

====================================================================================
"""

    assert stream.getvalue() == expected


@pytest.mark.unit
def test_display_variables_at_or_outside_bounds(diagnostics_toolbox):
    stream = StringIO()

    diagnostics_toolbox.display_variables_at_or_outside_bounds(stream=stream)

    expected = """====================================================================================
The following variable(s) have values at or outside their bounds (tol=0.0E+00):

    b.v3_1 (free): value=-9.997671315624657e-09 bounds=(0, 5)
    b.v3_2 (free): value=-9.9943899712328e-09 bounds=(0, 5)
    b.gb1.outputs[v3] (free): value=-9.998508737996996e-09 bounds=(0, 5)
    b.gb2.inputs[v3] (free): value=-9.994692548467552e-09 bounds=(0, 5)

====================================================================================
"""

    assert stream.getvalue() == expected


@pytest.mark.unit
def test_display_variables_with_none_value(diagnostics_toolbox):
    stream = StringIO()

    # Set some variable values to None for testing
    diagnostics_toolbox._model.gb1.outputs["v3"].set_value(None)
    diagnostics_toolbox._model.gb2.inputs["v3"].set_value(None)

    diagnostics_toolbox.display_variables_with_none_value(stream=stream)

    expected = """====================================================================================
The following variable(s) have a value of None:

    b.v8
    b.gb1.outputs[v3]
    b.gb2.inputs[v3]

====================================================================================
"""

    assert stream.getvalue() == expected


@pytest.mark.unit
def test_display_variables_with_none_value_in_activated_constraints(
    diagnostics_toolbox,
):
    stream = StringIO()

    # Set some variable values to None for testing
    diagnostics_toolbox._model.v4_1.set_value(None)
    diagnostics_toolbox._model.gb1.outputs["v3"].set_value(None)
    diagnostics_toolbox._model.gb2.inputs["v3"].set_value(None)

    diagnostics_toolbox.display_variables_with_none_value_in_activated_constraints(
        stream=stream
    )

    expected = """====================================================================================
The following variable(s) have a value of None:

    b.gb1.outputs[v3]
    b.v4_1
    b.gb2.inputs[v3]

====================================================================================
"""

    assert stream.getvalue() == expected


@pytest.mark.unit
def test_display_variables_with_value_near_zero(diagnostics_toolbox):
    stream = StringIO()

    diagnostics_toolbox.display_variables_with_value_near_zero(stream=stream)

    expected = """====================================================================================
The following variable(s) have a value close to zero (tol=1.0E-08):

    b.v3_1: value=-9.997671315624657e-09
    b.v6_1: value=0
    b.v3_2: value=-9.9943899712328e-09
    b.v6_2: value=0
    b.gb1.outputs[v3]: value=-9.998508737996996e-09
    b.gb2.inputs[v3]: value=-9.994692548467552e-09
    b.gb2.inputs[v6]: value=-6.81818317855146e-13

====================================================================================
"""

    assert stream.getvalue() == expected


@pytest.mark.unit
def test_display_variables_with_extreme_values(diagnostics_toolbox):
    stream = StringIO()

    diagnostics_toolbox.display_variables_with_extreme_values(stream=stream)

    expected = """====================================================================================
The following variable(s) have extreme values (<1.0E-04 or > 1.0E+04):

    b.v7_1: 1.0000000003005368e-07
    b.v7_2: 1.0000000003005368e-07
    b.gb1.outputs[v7]: 1.0000000001502684e-07
    b.gb2.inputs[v7]: 1.0000000001502684e-07

====================================================================================
"""

    assert stream.getvalue() == expected


@pytest.mark.unit
def test_display_variables_near_bounds(diagnostics_toolbox):
    stream = StringIO()

    diagnostics_toolbox.display_variables_near_bounds(stream=stream)

    expected = """====================================================================================
The following variable(s) have values close to their bounds (abs=1.0E-04, rel=1.0E-04):

    b.v3_1: value=-9.997671315624657e-09 bounds=(0, 5)
    b.v7_1: value=1.0000000003005368e-07 bounds=(0, 1)
    b.v3_2: value=-9.9943899712328e-09 bounds=(0, 5)
    b.v7_2: value=1.0000000003005368e-07 bounds=(0, 1)
    b.gb1.outputs[v3]: value=-9.998508737996996e-09 bounds=(0, 5)
    b.gb2.inputs[v3]: value=-9.994692548467552e-09 bounds=(0, 5)

====================================================================================
"""

    assert stream.getvalue() == expected


@pytest.mark.unit
def test_display_components_with_inconsistent_units(diagnostics_toolbox):
    stream = StringIO()

    diagnostics_toolbox.display_components_with_inconsistent_units(stream=stream)

    expected = """====================================================================================
The following component(s) have unit consistency issues:

    b.gb1_link_v1
    b.gb1_link_v2
    b.gb1_link_v7
    b.gb2_link_v1
    b.gb2_link_v2
    b.gb2_link_v7

For more details on unit inconsistencies, import the assert_units_consistent method
from pyomo.util.check_units
====================================================================================
"""

    assert stream.getvalue() == expected


@pytest.mark.unit
def test_display_constraints_with_large_residuals(diagnostics_toolbox):
    stream = StringIO()

    diagnostics_toolbox.display_constraints_with_large_residuals(stream=stream)

    expected = """====================================================================================
The following constraint(s) have large residuals (>1.0E-05):

    b.gb1_link_v5: 1.99380E+00
    b.gb1_link_v6: 1.11198E-02
    b.gb1.v3_constraint: 1.86107E-02
    b.gb2.c2: 6.66667E-01

====================================================================================
"""

    assert stream.getvalue() == expected


@pytest.mark.unit
def test_display_underconstrained_set(diagnostics_toolbox):
    stream = StringIO()

    # Create structural singularities
    blk = diagnostics_toolbox._model
    blk.v2_1.unfix()
    blk.v2_2.unfix()

    diagnostics_toolbox.display_underconstrained_set(stream=stream)

    expected = """====================================================================================
Dulmage-Mendelsohn Under-Constrained Set

    Independent Block 0:

        Variables:

            b.gb1.inputs[v2]
            b.gb1.outputs[v1]
            b.v2_1
            b.gb1.outputs[v7]
            v1_1
            b.v7_1

        Constraints:

            b.gb1.v1_constraint
            b.gb1_link_v2
            b.gb1.v7_constraint
            b.gb1_link_v1
            b.gb1_link_v7

    Independent Block 1:

        Variables:

            b.gb2.inputs[v2]
            b.gb2.inputs[v1]
            b.v2_2
            v1_2
            b.gb2.inputs[v7]
            b.v7_2

        Constraints:

            b.gb2.c1
            b.gb2_link_v2
            b.gb2_link_v1
            b.gb2.c4
            b.gb2_link_v7

====================================================================================
"""

    assert stream.getvalue() == expected


@pytest.mark.unit
def test_display_overconstrained_set(diagnostics_toolbox):
    stream = StringIO()

    # Create structural singularities
    blk = diagnostics_toolbox._model
    blk.v4_1.fix()
    blk.v4_2.fix()

    diagnostics_toolbox.display_overconstrained_set(stream=stream)

    expected = """====================================================================================
Dulmage-Mendelsohn Over-Constrained Set

    Independent Block 0:

        Variables:

            b.gb1.outputs[v4]
            b.gb1.inputs[v5]
            b.gb1.inputs[v6]

        Constraints:

            b.gb1_link_v4
            b.gb1_link_v5
            b.gb1_link_v6
            b.gb1.v4_constraint

    Independent Block 1:

        Variables:

            b.gb2.inputs[v4]
            b.gb2.inputs[v5]
            b.gb2.inputs[v6]
            b.gb2.inputs[v3]

        Constraints:

            b.gb2_link_v4
            b.gb2_link_v5
            b.gb2_link_v6
            b.gb2.c2
            b.gb2.c3

====================================================================================
"""

    assert stream.getvalue() == expected


@pytest.mark.component
def test_display_near_parallel_constraints(diagnostics_toolbox):
    m = diagnostics_toolbox.model

    # Add parallel constraints
    m.parallel_con1 = Constraint(expr=m.gb1.outputs["v1"] == 10 - m.gb1.inputs["v2"])
    m.parallel_con2 = Constraint(expr=m.gb2.inputs["v1"] == 10 - m.gb2.inputs["v2"])

    stream = StringIO()

    diagnostics_toolbox.display_near_parallel_constraints(stream=stream)

    expected = """====================================================================================
The following pairs of constraints are nearly parallel:

    b.parallel_con1, b.gb1.v1_constraint
    b.parallel_con2, b.gb2.c1

====================================================================================
"""

    assert stream.getvalue() == expected


@pytest.mark.component
def test_display_near_parallel_variables(diagnostics_toolbox):
    m = diagnostics_toolbox.model

    # Add parallel vars
    m.pv1 = Var(initialize=m.gb1.outputs["v1"].value)  # Parallel to gb1.outputs[v1]
    m.pv2 = Var(initialize=m.gb1.inputs["v2"].value)  # Parallel to gb1.inputs[v2]
    m.parallel_con1 = Constraint(expr=m.pv1 == 10 - m.pv2)

    stream = StringIO()

    diagnostics_toolbox.display_near_parallel_variables(stream=stream)

    expected = """====================================================================================
The following pairs of variables are nearly parallel:

    b.pv1, b.pv2

====================================================================================
"""

    assert stream.getvalue() == expected


# ====================================================================================
# Reporting methods
@pytest.mark.component
def test_collect_structural_warnings(diagnostics_toolbox):
    # Create structural singularities
    blk = diagnostics_toolbox._model
    blk.v4_1.fix()
    blk.v4_2.fix()

    warnings, next_steps = diagnostics_toolbox._collect_structural_warnings()

    assert len(warnings) == 3
    assert "WARNING: -2 Degrees of Freedom" in warnings
    assert "WARNING: 6 Components with inconsistent units" in warnings
    assert """WARNING: Structural singularity found
        Under-Constrained Set: 0 variables, 0 constraints
        Over-Constrained Set: 7 variables, 9 constraints""" in warnings

    assert len(next_steps) == 2
    assert "display_components_with_inconsistent_units()" in next_steps
    assert "display_overconstrained_set()" in next_steps


@pytest.mark.component
def test_collect_structural_cautions(diagnostics_toolbox):
    # Create structural singularities
    blk = diagnostics_toolbox._model
    blk.v4_1.fix()
    blk.v4_2.fix()

    cautions = diagnostics_toolbox._collect_structural_cautions()

    assert len(cautions) == 2
    assert "Caution: 2 variables fixed to 0" in cautions
    assert "Caution: 1 unused variable (0 fixed)" in cautions


@pytest.mark.component
def test_collect_numerical_warnings(diagnostics_toolbox):
    warnings, next_steps = diagnostics_toolbox._collect_numerical_warnings()

    assert len(warnings) == 2
    assert "WARNING: 4 Constraints with large residuals (>1.0E-05)" in warnings
    assert "WARNING: 4 Variables at or outside bounds (tol=0.0E+00)" in warnings

    assert len(next_steps) == 3
    assert "display_constraints_with_large_residuals()" in next_steps
    assert "compute_infeasibility_explanation()" in next_steps
    assert "display_variables_at_or_outside_bounds()" in next_steps


@pytest.mark.component
def test_collect_numerical_cautions(diagnostics_toolbox):
    cautions = diagnostics_toolbox._collect_numerical_cautions()

    for i in cautions:
        print(i)

    assert len(cautions) == 5
    assert (
        "Caution: 6 Variables with value close to their bounds (abs=1.0E-04, rel=1.0E-04)"
        in cautions
    )
    assert "Caution: 7 Variables with value close to zero (tol=1.0E-08)" in cautions
    assert "Caution: 4 Variables with extreme value (<1.0E-04 or >1.0E+04)" in cautions
    assert "Caution: 1 Variable with None value" in cautions
    assert "Caution: 2 extreme Jacobian Entries (<1.0E-04 or >1.0E+04)" in cautions


@pytest.mark.component
def test_assert_no_structural_warnings(diagnostics_toolbox):
    # Create structural singularities
    blk = diagnostics_toolbox._model
    blk.v4_1.fix()
    blk.v4_2.fix()

    with pytest.raises(AssertionError, match=re.escape("Structural issues found (3).")):
        diagnostics_toolbox.assert_no_structural_warnings()


@pytest.mark.component
def test_assert_no_numerical_warnings(diagnostics_toolbox):
    with pytest.raises(AssertionError, match=re.escape("Numerical issues found (2).")):
        diagnostics_toolbox.assert_no_numerical_warnings()


# ====================================================================================
# Reporting methods
@pytest.mark.component
def test_report_structural_issues(diagnostics_toolbox):
    # Create structural singularities
    blk = diagnostics_toolbox._model
    blk.v4_1.fix()
    blk.v4_2.fix()

    stream = StringIO()
    diagnostics_toolbox.report_structural_issues(stream=stream)

    expected = """====================================================================================
Model Statistics

        Activated Blocks: 1 (Deactivated: 0)
        Free Variables in Activated Constraints: 20 (External: 2)
            Free Variables with only lower bounds: 0
            Free Variables with only upper bounds: 0
            Free Variables with upper and lower bounds: 8
        Fixed Variables in Activated Constraints: 8 (External: 0)
        Activated Equality Constraints: 22 (Deactivated: 0)
        Activated Inequality Constraints: 0 (Deactivated: 0)
        Activated Objectives: 0 (Deactivated: 0)
        GreyBox Statistics
            Activated GreyBox models: 2 (Deactivated: 0)
            Activated GreyBox Equalities: 8 (Deactivated: 0)
            Free Variables in Activated GreyBox Equalities: 14 (Fixed: 0)

------------------------------------------------------------------------------------
3 WARNINGS

    WARNING: -2 Degrees of Freedom
    WARNING: 6 Components with inconsistent units
    WARNING: Structural singularity found
        Under-Constrained Set: 0 variables, 0 constraints
        Over-Constrained Set: 7 variables, 9 constraints

------------------------------------------------------------------------------------
2 Cautions

    Caution: 2 variables fixed to 0
    Caution: 1 unused variable (0 fixed)

------------------------------------------------------------------------------------
Suggested next steps:

    display_components_with_inconsistent_units()
    display_overconstrained_set()

====================================================================================
"""

    assert stream.getvalue() == expected


@pytest.mark.component
def test_report_numerical_issues(diagnostics_toolbox):
    stream = StringIO()
    diagnostics_toolbox.report_numerical_issues(stream=stream)

    print(stream.getvalue())
    expected = """====================================================================================
Model Statistics

    Jacobian Condition Number: 7.849E+01

------------------------------------------------------------------------------------
2 WARNINGS

    WARNING: 4 Constraints with large residuals (>1.0E-05)
    WARNING: 4 Variables at or outside bounds (tol=0.0E+00)

------------------------------------------------------------------------------------
5 Cautions

    Caution: 6 Variables with value close to their bounds (abs=1.0E-04, rel=1.0E-04)
    Caution: 7 Variables with value close to zero (tol=1.0E-08)
    Caution: 4 Variables with extreme value (<1.0E-04 or >1.0E+04)
    Caution: 1 Variable with None value
    Caution: 2 extreme Jacobian Entries (<1.0E-04 or >1.0E+04)

------------------------------------------------------------------------------------
Suggested next steps:

    display_constraints_with_large_residuals()
    compute_infeasibility_explanation()
    display_variables_at_or_outside_bounds()

====================================================================================
"""

    assert stream.getvalue() == expected


# ====================================================================================
# Functionality tests
# These tests are for methods that can/do not interact with Grey Box components,
# but we want to ensure work when Grey Box components are present in the model.
@pytest.mark.unit
def test_collect_constraint_mismatches(diagnostics_toolbox):
    blk = diagnostics_toolbox._model

    blk.b3 = Block()
    blk.b3.v1 = Var(initialize=2)
    blk.b3.v2 = Var(initialize=3)

    # Constraint with no free variables
    blk.b3.c1 = Constraint(expr=blk.b3.v1 == blk.b3.v2)
    blk.b3.v1.fix()
    blk.b3.v2.fix()

    # Constraint with mismatched terms
    blk.b3.v3 = Var(initialize=10)
    blk.b3.v4 = Var(initialize=10)
    blk.b3.v5 = Var(initialize=1e-6)
    blk.b3.c2 = Constraint(expr=blk.b3.v3 == blk.b3.v4 + blk.b3.v5)

    # Constraint with cancellation
    blk.b3.v6 = Var(initialize=10)
    blk.b3.c3 = Constraint(expr=blk.b3.v6 == 10 + blk.b3.v3 - blk.b3.v4)

    mismatch, canceling, constant = diagnostics_toolbox._collect_constraint_mismatches()

    print(mismatch)
    print(canceling)
    print(constant)

    assert mismatch == ["b.b3.c2: 1 mismatched term(s)"]
    assert canceling == [
        "b.b3.c2: 1 potential canceling term(s)",
        "b.b3.c3: 1 potential canceling term(s)",
    ]
    assert constant == ["b.b3.c1"]


@pytest.mark.component
def test_collect_potential_eval_errors(diagnostics_toolbox):
    blk = diagnostics_toolbox._model

    blk.b3 = Block()
    blk.b3.v1 = Var(initialize=2)
    blk.b3.v2 = Var(initialize=0)

    # Constraint with potential division by zero
    blk.b3.c1 = Constraint(expr=blk.b3.v1 == 10 / blk.b3.v2)

    errors = diagnostics_toolbox._collect_potential_eval_errors()

    assert len(errors) == 1
    assert "b.b3.c1: Potential division by 0 in 10/b.b3.v2" in errors[0]


@pytest.mark.component
def test_display_potential_evaluation_errors(diagnostics_toolbox):
    blk = diagnostics_toolbox._model

    blk.b3 = Block()
    blk.b3.v1 = Var(initialize=2)
    blk.b3.v2 = Var(initialize=0)

    # Constraint with potential division by zero
    blk.b3.c1 = Constraint(expr=blk.b3.v1 == 10 / blk.b3.v2)

    stream = StringIO()
    diagnostics_toolbox.display_potential_evaluation_errors(stream=stream)

    expected = """====================================================================================
1 WARNINGS

    b.b3.c1: Potential division by 0 in 10/b.b3.v2; Denominator bounds are (-inf, inf)

====================================================================================
"""

    assert stream.getvalue() == expected


@pytest.mark.component
def test_compute_infeasibility_explanation(diagnostics_toolbox):
    with pytest.raises(
        NotImplementedError,
        match="Model contains Greybox models, which are not supported by the "
        "Infeasibility Explanation tools at the moment",
    ):
        diagnostics_toolbox.compute_infeasibility_explanation()


@pytest.mark.component
def test_prepare_svd_toolbox(diagnostics_toolbox):
    with pytest.raises(
        NotImplementedError,
        match="Model contains Greybox models, which are not supported by the "
        "SVD Toolbox at the moment",
    ):
        diagnostics_toolbox.prepare_svd_toolbox()


@pytest.mark.component
def test_prepare_degeneracy_hunter(diagnostics_toolbox):
    with pytest.raises(
        NotImplementedError,
        match="Model contains Greybox models, which are not supported by "
        "Degeneracy Hunter at the moment",
    ):
        diagnostics_toolbox.prepare_degeneracy_hunter()
