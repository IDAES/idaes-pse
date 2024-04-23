from pyomo.core.expr.numeric_expr import NumericExpression
from pyomo.core.expr.sympy_tools import sympyify_expression, sympy2pyomo_expression
from pyomo.core.expr import is_fixed, value


def simplify_expr(expr: NumericExpression):
    om, se = sympyify_expression(expr)
    se = se.simplify()
    new_expr = sympy2pyomo_expression(se, om)
    if is_fixed(new_expr):
        new_expr = value(new_expr)
    return new_expr
