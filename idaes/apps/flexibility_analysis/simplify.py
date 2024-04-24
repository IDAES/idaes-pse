#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
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
