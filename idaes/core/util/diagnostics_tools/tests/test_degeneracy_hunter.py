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
This module contains tests for Degeneracy Hunter.
"""

from io import StringIO

import pytest

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Objective,
    Set,
    SolverFactory,
    value,
    Var,
)
from pyomo.contrib.pynumero.asl import AmplInterface
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
from pyomo.contrib.pynumero.interfaces.external_grey_box import (
    ExternalGreyBoxBlock,
    ExternalGreyBoxModel,
)

from idaes.core.util.diagnostics_tools.degeneracy_hunter import (
    DegeneracyHunter,
)
from idaes.core.util.testing import _enable_scip_solver_for_testing

__author__ = "Alex Dowling, Douglas Allan, Andrew Lee"


# TODO: Add pyomo.dae test cases
solver_available = SolverFactory("scip").available()


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
        dh = DegeneracyHunter(model)

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
            DegeneracyHunter(model=m)

    @pytest.mark.unit
    def test_get_solver(self, model):
        dh = DegeneracyHunter(model, solver="ipopt", solver_options={"maxiter": 50})

        solver = dh._get_solver()

        assert solver.options == {"maxiter": 50}

    @pytest.mark.unit
    def test_prepare_candidates_milp(self, model):
        dh = DegeneracyHunter(model)
        dh._prepare_candidates_milp()

        assert isinstance(dh.candidates_milp, ConcreteModel)

    @pytest.mark.unit
    def test_identify_candidates(self, model):
        dh = DegeneracyHunter(model)
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
        dh = DegeneracyHunter(model)
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
        dh = DegeneracyHunter(model)
        dh._prepare_ids_milp()

        assert isinstance(dh.ids_milp, ConcreteModel)

    @pytest.mark.unit
    def test_solve_ids_milp(self, model):
        dh = DegeneracyHunter(model)
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
        dh = DegeneracyHunter(model)
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
        dh = DegeneracyHunter(model)
        dh.find_irreducible_degenerate_sets()

        assert dh.irreducible_degenerate_sets == [
            {model.con2: 1, model.con5: -1},
            {model.con5: 1, model.con2: -1},
        ]

    @pytest.mark.solver
    @pytest.mark.component
    def test_report_irreducible_degenerate_sets(self, model, scip_solver):
        stream = StringIO()

        dh = DegeneracyHunter(model)
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

        dh = DegeneracyHunter(model)
        dh.report_irreducible_degenerate_sets(stream=stream)

        expected = """====================================================================================
Irreducible Degenerate Sets

    No candidate equations. The Jacobian is likely full rank.

====================================================================================
"""

        assert stream.getvalue() == expected
