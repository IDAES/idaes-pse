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
This module contains tests for the ipopt_solve_halt_on_error function.
"""

import pytest

from pyomo.environ import (
    ConcreteModel,
    Constraint,
    Expression,
    Var,
    log,
)

from idaes.core.util.diagnostics_tools.ipopt_halt_on_error import (
    ipopt_solve_halt_on_error,
)


@pytest.mark.component
def test_ipopt_solve_halt_on_error(capsys):
    m = ConcreteModel()

    m.v = Var(initialize=-5, bounds=(None, -1))
    m.e = Expression(expr=log(m.v))
    m.c = Constraint(expr=m.e == 1)

    try:
        _ = ipopt_solve_halt_on_error(m)
    except Exception:  # we expect this to fail
        pass

    captured = capsys.readouterr()
    assert "c: can't evaluate log(-5)." in captured.out
