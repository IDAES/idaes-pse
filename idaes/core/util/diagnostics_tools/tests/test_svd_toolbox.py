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
This module contains tests for the SVD Toolbox.
"""

from io import StringIO
from pyomo.contrib.pynumero.interfaces.external_grey_box import (
    ExternalGreyBoxBlock,
    ExternalGreyBoxModel,
)
import re

import numpy as np
import pytest

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Set,
    Var,
)
from pyomo.contrib.pynumero.asl import AmplInterface
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP

from idaes.core.util.diagnostics_tools.svd_toolbox import (
    SVDToolbox,
    svd_dense,
    svd_sparse,
)
from idaes.core.util.diagnostics_tools.tests.utils import (
    dummy_problem,
)

__author__ = "Alex Dowling, Douglas Allan, Andrew Lee"


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

    1st Smallest Singular Value: 1.000e-01

        Variables:

            x[3]

        Constraints:

            dummy_eqn[3]

    2nd Smallest Singular Value: 1.000e+00

        Variables:

            x[1]

        Constraints:

            dummy_eqn[1]

    3rd Smallest Singular Value: 5.000e+00

        Variables:

            x[4]

        Constraints:

            dummy_eqn[4]

    4th Smallest Singular Value: 1.000e+01

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

    1st Smallest Singular Value: 1.000e-01

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

    1st Smallest Singular Value: 1.000e-01

        Variables:


        Constraints:


    2nd Smallest Singular Value: 1.000e+00

        Variables:


        Constraints:


    3rd Smallest Singular Value: 5.000e+00

        Variables:


        Constraints:


    4th Smallest Singular Value: 1.000e+01

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
