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
by providing an expression to calculate them via a suffix.

The main purpose of this code is to use the calculate_scaling_factors function to
calculate scaling factors to be used with the Pyomo scaling transformation or
with solvers. A user can provide a scaling_expression suffix to calculate scale
factors from existing variable scaling factors. This allows scaling factors from
a small set of fundamental variables to be propagated to the rest of the model.

The scaling_expression suffix contains Pyomo expressions with model variables.
The expressions can be evaluated with variable scaling factors in place of
variables to calculate additional scaling factors.
"""

import enum
import pyomo.environ as pyo
from pyomo.core.expr import current as EXPR
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

class ScalingBasis(enum.Enum):
    """Basis value type for scaling expression calculations. These are values
    substituted into the scaling expressions in place of the variables and
    Expressions in the scaling expressions."""
    Value = 1 # use the variables current value
    VarScale = 2 # use the variable scale factor
    InverseVarScale = 3 # use 1/(variable scale factor) most common
    Lower = 4 # use the lower bound
    Upper = 5 # use the upper bound
    Mid = 6 # use the bound mid-point


def _replacement(m, basis):
    """PRIVATE FUNCTION
    Create a replacement visitor. The replacement visitor is used on
    user-provided scaling expressions. These expressions are written
    with model variables, but you generally don't want to calculate
    scaling factors based on the curent value of the model variables,
    you want to use their scaling factors, so the replacment visitor
    takes the user-defined scaling expression and replaces the model
    varible by some scaling factor, and returns a new expression. The
    basis argument can be used to specify the basis to use for scaling.

    Args:
        m (Block): model to collect vars from
        basis (list of ScalingBasis): value type to use as basis for scaling calcs

    Return:
        None or ExpressionReplacementVisitor

    """
    # These long ifs up front find values to replace variables in the scaling
    # expressions with.
    if basis[0] == ScalingBasis.Value:
        return None # no need to replace anything if using value
    else:
        rdict = {}
        for v in m.component_data_objects((pyo.Var)):
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
        for v in m.component_data_objects((pyo.Expression)):
            # check for expression scaling factors, while expressions don't
            # get scaled, the factor can be used in the calculation of other
            # scale factors.
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
                    else: # Expressions don't have bounds
                        continue
                except AttributeError:
                    pass
                except KeyError:
                    pass
            rdict[id(v)] = val
        # Use the substitutions dictionary from above to make a replacemnt visitor
        return EXPR.ExpressionReplacementVisitor(substitute=rdict)

def _calculate_scale_factors_from_nominal(m):
    """PRIVATE FUNCTION
    For variables and expressions, if a nominal value is provided calculate the
    scaling factor.

    Args:
        m (Block): a pyomo block to calculate scaling factors for

    Returns:
        None
    """

    for c in m.component_data_objects((pyo.Var, pyo.Expression)):
        # Check for a scaling expression.  If there is one, use it to calculate
        # a scaling factor otherwise use autoscaling.
        if not hasattr(c.parent_block(), "nominal_value"):
            continue # no scaling expression supplied
        elif c not in c.parent_block().nominal_value:
            continue # no scaling expression supplied

        if not hasattr(c.parent_block(), "scaling_factor"):
            # if there is no scaling_factor Suffix yet make one
            c.parent_block().scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)

        # Add scaling factor from nominal value of variables or expressions
        c.parent_block().scaling_factor[c] = 1/c.parent_block().nominal_value[c]


def _calculate_scale_factors_from_expr(m, replacement, cls):
    """PRIVATE FUNCTION
    Take the expressions from the scaling_expression suffix and use them to
    calculate scaling factors for the scaling_factor suffix that is used by Pyomo
    or the solver to do variable and constraint scaling. The resulting scaling
    factors are put into the scaling factor suffix.

    Args:
        m (Block): a pyomo block to calculate scaling factors for
        replacement (ReplacementVisitor): A pyomo replacment visitor to replace
            the variable in a scaling factor expression from the scaling_factor
            suffix and return a new expression for calculating scaling factors
        cls: The class to calculate scaling factors for Var or Constraint

    Returns:
        None
    """

    # Calculate scaling factors for each constraint
    for c in m.component_data_objects(cls):
        # Check for a scaling expression.  If there is one, use it to calculate
        # a scaling factor otherwise use autoscaling.
        if not hasattr(c.parent_block(), "scaling_expression"):
            continue # no scaling expression supplied
        elif c not in c.parent_block().scaling_expression:
            continue # no scaling expression supplied

        if not hasattr(c.parent_block(), "scaling_factor"):
            # if there is no scaling_factor Suffix yet make one
            c.parent_block().scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)

        # Take scaling expression provided by modeler and put in basis values
        if replacement is None:
            expr = c.parent_block().scaling_expression[c]
        else:
            expr = replacement.dfs_postorder_stack(
                c.parent_block().scaling_expression[c]
            )

        # Add constraint scaling factor by evaluating modeler provided scale expr
        c.parent_block().scaling_factor[c] = pyo.value(expr)


def calculate_scaling_factors(
    m,
    basis=(
        ScalingBasis.InverseVarScale,
        ScalingBasis.Mid,
        ScalingBasis.Value,
    )
):
    """Set scale factors for variables and constraints from expressions stored in
    the scaling_expression suffix. The variables and Expressions in the scaling
    expressions are replaced by the scaling basis values before calculating
    the scaling factor. Variable scale factors are calculated first, and variable
    scaling expressions should be based on variables whose scale factors are
    supplied directly. Constraint scaling expressions can be based on any variables.

    Args:
        m (Block): A Pyomo model or block to apply the scaling expressions to.
        basis: (ScalingBasis or List-like of ScalingBasis): Value to use
            when evaluating scaling expressions. A list-like of ScalingBasis can
            be used to provide fall-back values in the event that the first
            choice is not available.  If none of the bases are available, 1 is used.

    Returns:
        None
    """
    # Map the scaling expression calculation values to the variables, and get
    # a replacement visitor to swap variable values for basis values
    if isinstance(basis, ScalingBasis):
         basis = (basis, )
    replacement = _replacement(m, basis)

    # If nominal values are supplied for Vars or named Expressions, use them
    # to set scaling factors
    _calculate_scale_factors_from_nominal(m)

    # First calculate variable scale factors, where expressions were provided
    _calculate_scale_factors_from_expr(m, replacement=replacement, cls=pyo.Var)
    # Then constraints
    _calculate_scale_factors_from_expr(m, replacement=replacement, cls=pyo.Constraint)


def badly_scaled_var_generator(blk, large=1e4, small=1e-3, zero=1e-10):
    """This provides a rough check for variables with poor scaling based on their
    current scale factors and values. For each potentially poorly scaled variable
    it returns the var and its current scaled value.

    Args:
        large: Magnitude that is considered to be too large
        small: Magnitude that is considered to be too small
        zero: Magnitude that is considered to be zero, variables with a value of
            zero are okay, and not reported.

    Yields:
        variable data object, current absolute value of scaled value
    """
    for v in blk.component_data_objects(pyo.Var):
        if v.fixed:
            continue
        try:
            sf = v.parent_block().scaling_factor.get(v, 1)
        except AttributeError: # no scaling factor suffix
            sf = 1
        sv = abs(pyo.value(v)*sf) # scaled value
        if sv > large:
            yield v, sv
        elif sv < zero:
            continue
        elif sv < small:
            yield v, sv


def grad_fd(c, scaled=False, h=1e-6):
    """Finite difference the gradient for a constraint, objective or named
    expression.  This is only for use in examining scaling.  For faster more
    accurate gradients refer to pynumero.

    Args:
        c: constraint to evaluate
        scaled: if True calculate the scaled grad (default=False)
        h: step size for calculating finite differnced derivatives

    Returns:
        (list of gradient values, list for varibles in the constraint) The order
            of the variables coresoponds to the gradient values.
    """
    try:
        ex = c.body
    except AttributeError:
        ex = c.expr

    vars = list(EXPR.identify_variables(ex))
    grad = [None]*len(vars)

    r = {}
    if scaled: # If you want the scaled grad put in variable scaling tansform
        orig = [pyo.value(v) for v in vars]
        for i, v in enumerate(vars):
            try:
                sf = v.parent_block().scaling_factor.get(v, 1)
            except AttributeError:
                sf = 1
            r[id(v)] = v/sf
            v.value = orig[i]*sf
        vis = EXPR.ExpressionReplacementVisitor(
            substitute=r,
            remove_named_expressions=True,
        )
        e = vis.dfs_postorder_stack(ex)
    else:
        e = ex

    for i, v in enumerate(vars):
        ov = pyo.value(v) # original variable value
        f1 = pyo.value(e)
        v.value = ov + h
        f2 = pyo.value(e)
        v.value = ov
        if scaled:
            try:
                sf = c.parent_block().scaling_factor.get(c, 1)
            except AttributeError:
                sf = 1
            grad[i] = sf*(f2 - f1)/h
        else:
            grad[i] = (f2 - f1)/h

    if scaled:
        for i, v in enumerate(vars):
            v.value = orig[i]

    return grad, vars


def scale_constraint(c, v=None):
    """This transforms a constraint with its scaling factor or a given scaling
    factor value.  If it uses the scaling factor suffix value, the scaling factor
    suffix is set to 1 to avoid double scaling the constraint.  This can be used
    when to scale constraints before sending the model to the solver.

    Args:
        c: Pyomo constraint
        v: Scale factor. If None, use value from scaling factor suffix and set
            suffix value to 1.

    Returns:
        None
    """
    if v is None:
        try:
            v = c.parent_block().scaling_factor[c]
            c.parent_block().scaling_factor[c] = 1
        except:
            v = None

    if v is not None:
        c.set_value((c.lower*v, c.body*v, c.upper*v))


def constraint_fd_autoscale(c, min_scale=1e-6, max_grad=100):
    """Autoscale constraints so that if there maximum partial derivative with
    respect to any variable is greater than max_grad at the current variable
    values, the method will attempt to assign a scaling factor to the constraint
    that makes the maximum derivative max_grad.  The min_scale value provides a
    lower limit allowed for constraint scaling factors.  If the caclulated
    scaling factor to make the maxium derivative max_grad is less than min_scale,
    min_scale is used instead.  Derivatives are approximated using finite
    differnce.

    Args:
        c: constraint object
        max_grad: the largest derivative after scaling subject to min_scale
        min_scale: the minimum scale factor allowed

    Returns:
        None
    """
    g, v = grad_fd(c, scaled=True)
    if not hasattr(c.parent_block(), "scaling_factor"):
        c.parent_block().scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
    s0 = c.parent_block().scaling_factor.get(c, 1)
    maxg = max(map(abs,g))
    ming = min(map(abs,g))
    if maxg > max_grad:
        c.parent_block().scaling_factor[c] = s0*max_grad/maxg
        if c.parent_block().scaling_factor[c] < min_scale:
            c.parent_block().scaling_factor[c] = min_scale


def set_scaling_factor(c, v):
    """Set a scaling factor for a model component.  This function creates the
    scaling_factor suffix if needed.

    Args:
        c: component to supply scaling factor for
        v: scaling factor
    Returns:
        None
    """
    try:
        c.parent_block().scaling_factor[c] = v
    except AttributeError:
        c.parent_block().scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
        c.parent_block().scaling_factor[c] = v
