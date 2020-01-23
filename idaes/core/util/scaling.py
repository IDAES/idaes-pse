import enum
import pyomo.environ as pyo
from pyomo.contrib.pynumero.interfaces import PyomoNLP
from pyomo.core.expr import current as EXPR
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)

class ScalingBasis(enum.Enum):
    """Basis value type for scaling expression calculations."""
    VALUE = 1
    VARSCALE = 2
    INVVARSCALE = 3
    LOWER = 4
    UPPER = 5
    MID = 6


def _replace(expr, replacment):
    """Use the replacment visitor to replace variables in an expression by the
    basis used for calculating scale factors.

    Args:
        expr: expression to replace variables in
        replacment: a replacment visitor, or None to leave unchanged
    Returns:
        expression"""
    if replacment is None:
        return expr
    else:
        return replacment.dfs_postorder_stack(expr)


def _replacment(m, basis):
    """Create a replacment visitor for model m to replace variables by the
    scaling basis.
    Args:
        m (Block): model to collect vars from
        basis (list of ScalingBasis): value type to use as basis for scaling calcs

    Return:
        None or ExpressionReplacementVisitor
    """
    if basis[0] == ScalingBasis.VALUE:
        return None # no need to replace anything if using value
    else:
        rdict = {}
        for v in m.component_data_objects(pyo.Var):
            val = 1.0
            for b in basis:
                try:
                    if b == ScalingBasis.VARSCALE:
                        val = v.parent_block().scaling_factor[v]
                        break
                    elif b == ScalingBasis.INVVARSCALE:
                        val = 1/v.parent_block().scaling_factor[v]
                        break
                    elif b == ScalingBasis.VALUE:
                        val = pyo.value(v)
                        break
                    elif b == ScalingBasis.MID:
                        if v.lb is not None and v.ub is not None:
                            val = (v.ub + v.lb)/2.0
                            break
                    elif b == ScalingBasis.LOWER:
                        if v.lb is not None:
                            val = v.lb
                            break
                    elif b == ScalingBasis.UPPER:
                        if v.ub is not None:
                            val = v.ub
                            break
                    else:
                        _log.warning("Unknown scaling expression basis {}".format(int(b)))
                except AttributeError:
                    pass
                except KeyError:
                    pass
            rdict[id(v)] = val
        return EXPR.ExpressionReplacementVisitor(substitute=rdict)


def _calculate_scale_factors_from_expr(m, replacment, cls):
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
        expr = _replace(c.parent_block().scaling_expression[c], replacment)
        #Add constraint scaling factor by evaluating modeler provided scale expr
        c.parent_block().scaling_factor[c] = pyo.value(expr)


def apply_scaling(
    m,
    basis=(
        ScalingBasis.INVVARSCALE,
        ScalingBasis.MID,
        ScalingBasis.VALUE,
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
    # a replacment visitor to swap vairable values for basis values
    if isinstance(basis, ScalingBasis):
         basis = (basis, )
    replacment = _replacment(m, basis)

    # Fist calculate variable scale factors, where expressions where provided
    _calculate_scale_factors_from_expr(m, replacment=replacment, cls=pyo.Var)
    # Then constraints
    _calculate_scale_factors_from_expr(m, replacment=replacment, cls=pyo.Constraint)
