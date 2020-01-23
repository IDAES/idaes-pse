##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
This module contains utilities to provide variable and expression scaling factors
by providing an expression to calculate them via a suffix
"""

import enum
import pyomo.environ as pyo
from pyomo.contrib.pynumero.interfaces import PyomoNLP
from pyomo.core.expr import current as EXPR
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

class ScalingBasis(enum.Enum):
    """Basis value type for scaling expression calculations."""
    Value = 1 # use the variables current value
    VarScale = 2 # use the variable scale
    InverseVarScale = 3 # use 1/(variable scale factor) most common
    Lower = 4 # use the lower bound
    Upper = 5 # use the upper bound
    Mid = 6 # use the bound mid-point


def _replace(expr, replacement):
    """Replace variables in an expression by the basis value used for
    calculating scale factors.

    Args:
        expr: expression to replace variables in
        replacement: a replacement visitor used to walk the expression tree,
            or None to leave unchanged
    Returns:
        expression"""
    if replacement is None:
        return expr
    else:
        return replacement.dfs_postorder_stack(expr)


def _replacement(m, basis):
    """Create a replacement visitor for model m to replace variables by the
    scaling basis.  The replacemnt visitor walks an expression tree and replaces
    variables, by a value to be used in the scaling calculation.

    Args:
        m (Block): model to collect vars from
        basis (list of ScalingBasis): value type to use as basis for scaling calcs

    Return:
        None or ExpressionReplacementVisitor
    """
    if basis[0] == ScalingBasis.Value:
        return None # no need to replace anything if using value
    else:
        rdict = {}
        for v in m.component_data_objects(pyo.Var):
            val = 1.0
            for b in basis:
                try:
                    if b == ScalingBasis.VarScale:
                        val = v.parent_block().scaling_factor[v]
                        break
                    elif b == ScalingBasis.InverseVarScale:
                        val = 1/v.parent_block().scaling_factor[v]
                        break
                    elif b == ScalingBasis.Value:
                        val = pyo.value(v)
                        break
                    elif b == ScalingBasis.Mid:
                        if v.lb is not None and v.ub is not None:
                            val = (v.ub + v.lb)/2.0
                            break
                    elif b == ScalingBasis.Lower:
                        if v.lb is not None:
                            val = v.lb
                            break
                    elif b == ScalingBasis.Upper:
                        if v.ub is not None:
                            val = v.ub
                            break
                    else:
                        _log.warning("Unknown scaling expression basis {}".format(b))
                except AttributeError:
                    pass
                except KeyError:
                    pass
            rdict[id(v)] = val
        return EXPR.ExpressionReplacementVisitor(substitute=rdict)


def _calculate_scale_factors_from_expr(m, replacement, cls):
    # Calculate scaling factors for each constraint
    for c in m.component_data_objects(cls):
        # Check for a scaling expression.  If there is one, use it to calculate
        # a scaling factor.  If use autoscaling.
        if not hasattr(c.parent_block(), "scaling_expression"):
            continue # no scaling expression supplied
        elif c not in c.parent_block().scaling_expression:
            continue # no scaling expression supplied
        if not hasattr(c.parent_block(), "scaling_factor"):
            # if there is no scaling_factor Suffix yet make one
            c.parent_block().scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)

        #Take scaling expression provided by modeler and put in basis values
        expr = _replace(c.parent_block().scaling_expression[c], replacement)
        #Add constraint scaling factor by evaluating modeler provided scale expr
        c.parent_block().scaling_factor[c] = pyo.value(expr)


def apply_scaling(
    m,
    basis=(
        ScalingBasis.InverseVarScale,
        ScalingBasis.Mid,
        ScalingBasis.Value,
    )
):
    """Set scale factors for variables and constraints from expressions, which
    calcualte them based on supplied variable scale factors, values, or bounds.

    Args:
        basis: (ScalingBasis or List-like of ScalingBasis): Value to use
            when evaluating scaling expression either the current variable value.
            A list can be provided to allow falling back on a differnt value if
            one is not available.

    Returns:
        None
    """
    # Map the scaling expression calculation values to the variables, and get
    # a replacement visitor to swap variable values for basis values
    if isinstance(basis, ScalingBasis):
         basis = (basis, )
    replacement = _replacement(m, basis)

    # Fist calculate variable scale factors, where expressions where provided
    _calculate_scale_factors_from_expr(m, replacement=replacement, cls=pyo.Var)
    # Then constraints
    _calculate_scale_factors_from_expr(m, replacement=replacement, cls=pyo.Constraint)
