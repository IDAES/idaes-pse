from pyomo.core.expr.numeric_expr import MonomialTermExpression
from apps.matopt.opt.pyomo_modeling import getLB, getUB
from pyomo.environ import *


def test_getLB():
    mod = ConcreteModel()
    mod.x = Var(within=NonNegativeReals, bounds=(0, 2), initialize=1)
    mod.e = MonomialTermExpression((1, mod.x))
    assert isinstance(mod.e, MonomialTermExpression)
    LB = getLB(mod.e)


def test_getUB():
    mod = ConcreteModel()
    mod.x = Var(within=NonNegativeReals, bounds=(0, 2), initialize=1)
    mod.e = MonomialTermExpression((1, mod.x))
    assert isinstance(mod.e, MonomialTermExpression)
    UB = getUB(mod.e)
