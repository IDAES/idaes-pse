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
This module contains tests for the ill conditioning tools.
"""

import re
import pytest

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Var,
)

from idaes.core.util.diagnostics_tools.ill_conditioning import (
    compute_ill_conditioning_certificate,
)

__author__ = "Alex Dowling, Douglas Allan, Andrew Lee"


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
