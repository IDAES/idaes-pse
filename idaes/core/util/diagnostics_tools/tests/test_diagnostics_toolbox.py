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
This module contains model diagnostic utility functions for use in IDAES (Pyomo) models.
"""

from io import StringIO
import re

import pytest

from pyomo.environ import (
    Block,
    ConcreteModel,
    Constraint,
    Expression,
    Objective,
    PositiveIntegers,
    Set,
    units,
    Var,
    exp,
)
from pyomo.contrib.pynumero.asl import AmplInterface
from pyomo.contrib.pynumero.interfaces.external_grey_box import (
    ExternalGreyBoxBlock,
    ExternalGreyBoxModel,
)

from idaes.core.solvers import get_solver
from idaes.core.util.diagnostics_tools.diagnostics_toolbox import (
    DiagnosticsToolbox,
)
from idaes.core.util.diagnostics_tools.degeneracy_hunter import (
    DegeneracyHunter,
)
from idaes.core.scaling import set_scaling_factor
from idaes.core.util.diagnostics_tools.svd_toolbox import (
    SVDToolbox,
)

__author__ = "Alex Dowling, Douglas Allan, Andrew Lee"


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
The following variable(s) correspond to Jacobian columns with extreme norms(<1.0E-04 or>1.0E+04):

    v2: 1.000E+10
    v1: 1.000E+08
    v3: 1.000E-06

====================================================================================
"""

        assert stream.getvalue() == expected

        # Test to make sure scaled Jacobian is used
        set_scaling_factor(model.v3, 1e-6)

        stream = StringIO()
        dt.display_variables_with_extreme_jacobians(stream)

        expected = """====================================================================================
The following variable(s) correspond to Jacobian columns with extreme norms(<1.0E-04 or>1.0E+04):

    v2: 1.000E+10
    v1: 1.000E+08

====================================================================================
"""

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
The following constraint(s) correspond to Jacobian rows with extreme norms (<1.0E-04 or>1.0E+04):

    c3: 1.000E+10

====================================================================================
"""

        assert stream.getvalue() == expected

        # Test to make sure scaled Jacobian is used
        set_scaling_factor(model.c3, 1e-8)

        stream = StringIO()
        dt.display_constraints_with_extreme_jacobians(stream)

        expected = """====================================================================================
The following constraint(s) correspond to Jacobian rows with extreme norms (<1.0E-04 or>1.0E+04):


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
entries (<1.0E-04 or>1.0E+04):

    c3, v2: 1.000E+10
    c2, v3: 1.000E-08
    c3, v1: 1.000E+08
    c3, v3: 1.000E-06

====================================================================================
"""

        assert stream.getvalue() == expected

        # Test to make sure scaled Jacobian is used
        set_scaling_factor(model.c3, 1e-8)

        stream = StringIO()
        dt.display_extreme_jacobian_entries(stream)

        expected = """====================================================================================
The following constraint(s) and variable(s) are associated with extreme Jacobian
entries (<1.0E-04 or>1.0E+04):

    c3, v3: 1.000E-14
    c2, v3: 1.000E-08

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
        assert """WARNING: Structural singularity found
        Under-Constrained Set: 3 variables, 2 constraints
        Over-Constrained Set: 0 variables, 0 constraints""" in warnings

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
        assert """WARNING: Structural singularity found
        Under-Constrained Set: 0 variables, 0 constraints
        Over-Constrained Set: 1 variables, 2 constraints""" in warnings

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
            "WARNING: 2 Variables with extreme Jacobian column norms (<1.0E-08 or >1.0E+08)"
            in warnings
        )
        assert (
            "WARNING: 1 Constraint with extreme Jacobian row norms (<1.0E-08 or >1.0E+08)"
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
            "Caution: 3 Variables with extreme Jacobian column norms (<1.0E-04 or >1.0E+04)"
            in cautions
        )
        assert (
            "Caution: 1 Constraint with extreme Jacobian row norms (<1.0E-04 or >1.0E+04)"
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
    WARNING: 2 Variables with extreme Jacobian column norms (<1.0E-08 or >1.0E+08)
    WARNING: 1 Constraint with extreme Jacobian row norms (<1.0E-08 or >1.0E+08)
    WARNING: 3 pairs of variables are parallel (to tolerance 1.0E-08)

------------------------------------------------------------------------------------
4 Cautions

    Caution: 3 Variables with value close to zero (tol=1.0E-08)
    Caution: 3 Variables with extreme Jacobian column norms (<1.0E-04 or >1.0E+04)
    Caution: 1 Constraint with extreme Jacobian row norms (<1.0E-04 or >1.0E+04)
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

        # Test that scaled Jacobian is used
        set_scaling_factor(model.c3, 1e-8)
        set_scaling_factor(model.v3, 1e-6)

        stream = StringIO()
        dt.report_numerical_issues(stream)

        expected = """====================================================================================
Model Statistics

    Jacobian Condition Number: 1.225E+04

------------------------------------------------------------------------------------
0 WARNINGS

    No warnings found!

------------------------------------------------------------------------------------
2 Cautions

    Caution: 3 Variables with value close to zero (tol=1.0E-08)
    Caution: 1 extreme Jacobian Entry (<1.0E-04 or >1.0E+04)

------------------------------------------------------------------------------------
Suggested next steps:

    If you still have issues converging your model consider:

        prepare_degeneracy_hunter()
        prepare_svd_toolbox()

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

        assert isinstance(dh, DegeneracyHunter)
