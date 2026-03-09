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
This module contains utility functions for diagnostics tool tests.
"""

import pytest

from pyomo.environ import (
    ConcreteModel,
    Objective,
    Set,
    Var,
)

__author__ = "Alex Dowling, Douglas Allan, Andrew Lee"


@pytest.fixture()
def dummy_problem():
    """
    Create a dummy problem for testing purposes.

    Returns:
        A dummy Pyomo model.
    """
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
