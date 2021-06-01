#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
from pyomo.core.expr.numeric_expr import MonomialTermExpression
from idaes.apps.matopt.opt.pyomo_modeling import getLB, getUB
from pyomo.environ import *
import pytest


@pytest.mark.unit
def test_getLB():
    mod = ConcreteModel()
    mod.x = Var(within=NonNegativeReals, bounds=(0, 2), initialize=1)
    mod.e = MonomialTermExpression((1, mod.x))
    assert isinstance(mod.e, MonomialTermExpression)
    LB = getLB(mod.e)


@pytest.mark.unit
def test_getUB():
    mod = ConcreteModel()
    mod.x = Var(within=NonNegativeReals, bounds=(0, 2), initialize=1)
    mod.e = MonomialTermExpression((1, mod.x))
    assert isinstance(mod.e, MonomialTermExpression)
    UB = getUB(mod.e)
