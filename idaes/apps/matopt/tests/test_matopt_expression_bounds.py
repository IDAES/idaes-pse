###############################################################################
# ** Copyright Notice **
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021 by the
# software owners: The Regents of the University of California, through Lawrence
# Berkeley National Laboratory,  National Technology & Engineering Solutions of
# Sandia, LLC, Carnegie Mellon University, West Virginia University Research
# Corporation, et al.  All rights reserved.
#
# NOTICE.  This Software was developed under funding from the U.S. Department of
# Energy and the U.S. Government consequently retains certain rights. As such, the
# U.S. Government has been granted for itself and others acting on its behalf a
# paid-up, nonexclusive, irrevocable, worldwide license in the Software to
# reproduce, distribute copies to the public, prepare derivative works, and
# perform publicly and display publicly, and to permit other to do so.
###############################################################################
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
