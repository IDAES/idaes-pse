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
#  ___________________________________________________________________________
#
#  Pyomo: Python Optimization Modeling Objects
#  Copyright 2017 National Technology and Engineering Solutions of Sandia, LLC
#  Under the terms of Contract DE-NA0003525 with National Technology and
#  Engineering Solutions of Sandia, LLC, the U.S. Government retains certain
#  rights in this software.
#  This software is distributed under the 3-clause BSD License.
#  ___________________________________________________________________________

__author__ = "Oluwamayowa Amusat, John Siirola"

from pyomo.core.expr import current as EXPR, native_types
from pyomo.core.expr.numvalue import value

_numpy_available = True
try:
    import numpy

    _functionMap = {
        "exp": numpy.exp,
        "log": numpy.log,
        "log10": numpy.log10,
        "sin": numpy.sin,
        "asin": numpy.arcsin,
        "sinh": numpy.sinh,
        "asinh": numpy.arcsinh,
        "cos": numpy.cos,
        "acos": numpy.arccos,
        "cosh": numpy.cosh,
        "acosh": numpy.arccosh,
        "tan": numpy.tan,
        "atan": numpy.arctan,
        "tanh": numpy.tanh,
        "atanh": numpy.arctanh,
        "ceil": numpy.ceil,
        "floor": numpy.floor,
        "sqrt": numpy.sqrt,
    }
except ImportError:
    _numpy_available = False


class NumpyEvaluator(EXPR.StreamBasedExpressionVisitor):
    def __init__(self, object_map):
        super(NumpyEvaluator, self).__init__()
        self.object_map = object_map

    def exitNode(self, node, values):
        if (
            node.__class__ is EXPR.UnaryFunctionExpression
            or node.__class__ is EXPR.NPV_UnaryFunctionExpression
        ):
            return _functionMap[node._name](values[0])
        if (
            node.__class__ is EXPR.AbsExpression
            or node.__class__ is EXPR.NPV_AbsExpression
        ):
            return numpy.abs(values[0])
        return node._apply_operation(values)

    def beforeChild(self, node, child, child_idx):
        #
        # Don't replace native types
        #
        if type(child) in native_types:
            return False, child
        #
        # We will descend into all expressions...
        #
        if child.is_expression_type():
            return True, None
        #
        # Replace pyomo variables with numpy variables
        #
        if child in self.object_map:
            return False, self.object_map[child]
        #
        # Assume everything else is a constant...
        #
        return False, value(child)
