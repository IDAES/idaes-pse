##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
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
This module contains utilities to provide variable and expression scaling
factors by providing an expression to calculate them via a suffix.

The main purpose of this code is to use the calculate_scaling_factors function
to calculate scaling factors to be used with the Pyomo scaling transformation or
with solvers. A user can provide a scaling_expression suffix to calculate scale
factors from existing variable scaling factors. This allows scaling factors from
a small set of fundamental variables to be propagated to the rest of the model.

The scaling_expression suffix contains Pyomo expressions with model variables.
The expressions can be evaluated with variable scaling factors in place of
variables to calculate additional scaling factors.
"""

__author__ = "John Eslick, Tim Bartholomew"

import pyomo.environ as pyo
from pyomo.core.base.constraint import _ConstraintData
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
from pyomo.common.modeling import unique_component_name
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def __none_mult(x, y):
    """PRIVATE FUNCTION, Multiply x by y and if x is None return None"""
    if x is not None:
        x *= y
    return x


def map_scaling_factor(iter, default=1, warning=True, func=min):
    """Map get_scaling_factor to an iterable of Pyomo compoents, and call func
    on the result.  This could be use, for example, to get the minimum or
    maximum scaling factor of a set of compoents.

    Args:
        iter: Iterable yeilding Pyomo componentes
        default: The default value used when a scaling factor is missing. The
            default is default=1.
        warning: Log a warning for missing scaling factors
        func: The function to call on the resulting iterable of scaling factors.
            The default is min().

    Returns:
        The result of func on the set of scaling factors
    """
    return func(
        map(
            lambda x: get_scaling_factor(x, default=default, warning=warning),
            iter
        )
    )


def min_scaling_factor(iter, default=1, warning=True):
    """Map get_scaling_factor to an iterable of Pyomo compoents, and get the
    minimum caling factor.

    Args:
        iter: Iterable yeilding Pyomo componentes
        default: The default value used when a scaling factor is missing.  If
            None, this will raise an exception when scaling factors are missing.
            The default is default=1.
        warning: Log a warning for missing scaling factors

    Returns:
        Minimum scaling factor of the compoents in iter
    """
    return map_scaling_factor(iter, default=default, warning=warning, func=min)


def propagate_indexed_component_scaling_factors(
    blk,
    typ=(pyo.Var, pyo.Constraint, pyo.Expression, pyo.Param),
    overwrite=False,
    descend_into=True):
    """Use the parent compoent scaling factor to set all component data object
    scaling factors.

    Args:
        blk: The block on which to search for components
        typ: Component type(s) (default=(Var, Constraint, Expression, Param))
        overwrite: if a data object already has a scaling factor should it be
            overwrittten (default=False)
        descend_into: descend into child blocks (default=True)
    """
    for c in blk.component_objects(typ, descend_into=descend_into):
        if get_scaling_factor(c) is not None and c.is_indexed():
            for cdat in c.values():
                if overwrite or get_scaling_factor(cdat) is None:
                    set_scaling_factor(cdat, get_scaling_factor(c))


def calculate_scaling_factors(blk):
    """Look for calculate scaling factor methods this uses a depth first
    ordering, so sub-block scale factors are called first. Scale factor
    calculations should only depend on descendent blocks.
    """
    def cs(blk2):
        """ Recursive function for depth first ordering """
        for b in blk2.component_data_objects(pyo.Block, descend_into=False):
            cs(b)
        if hasattr(blk2, "calculate_scaling_factors"):
            blk2.calculate_scaling_factors()
    cs(blk)
    propagate_indexed_component_scaling_factors(blk)


def set_scaling_factor(c, v, data_objects=True):
    """Set a scaling factor for a model component. This function creates the
    scaling_factor suffix if needed.

    Args:
        c: component to supply scaling factor for
        v: scaling factor
    Returns:
        None
    """
    if isinstance(c, (float, int)):
        # property packages can return 0 for material balance terms on compoents
        # doesn't exist.  This handels the case where you get a constant 0 and
        # need it's scale factor to scale the mass balance.
        return 1
    try:
        suf = c.parent_block().scaling_factor
    except AttributeError:
        c.parent_block().scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
        suf = c.parent_block().scaling_factor

    suf[c] = v
    if data_objects and c.is_indexed():
        for cdat in c.values():
            suf[cdat] = v


def get_scaling_factor(c, default=None, warning=False, exception=False):
    """Get a component scale factor.

    Args:
        c: compoent
        default: value to return if no scale factor exists (default=None)
    """
    try:
        sf = c.parent_block().scaling_factor.get(c, default)
    except AttributeError:
        if warning:
            _log.warning(f"Accessing missing scale factor for {c}")
        if exception:
            _log.error(f"Accessing missing scale factor for {c}")
            raise
        sf = default
    return sf


def unset_scaling_factor(c):
    """Delete a component scaling factor.

    Args:
        c: component

    Returns:
        None
    """
    try:
        del c.parent_block().scaling_factor[c]
    except AttributeError:
        pass # no scaling factor suffix, is fine
    except KeyError:
        pass # no scaling factor is fine


def __set_constarint_tranform_applied_scaling_factor(c, v):
    """PRIVATE FUNCTION Set the scaling factor used to transform a constraint.
    This is used to keep track of scaling tranformations that have been applied
    to constraints.

    Args:
        c: component to supply scaling factor for
        v: scaling factor
    Returns:
        None
    """
    try:
        c.parent_block().constaint_transformed_scaling_factor[c] = v
    except AttributeError:
        c.parent_block().constaint_transformed_scaling_factor = pyo.Suffix(
            direction=pyo.Suffix.LOCAL)
        c.parent_block().constaint_transformed_scaling_factor[c] = v


def get_constarint_tranform_applied_scaling_factor(c, default=None):
    """Get a the scale factor that was used to transform a
    constraint.

    Args:
        c: constraint data object
        default: value to return if no scaling factor exisits (default=None)

    Returns:
        The scaling factor that has been used to transform the constrain or the
        default.
    """
    try:
        sf = c.parent_block().constaint_transformed_scaling_factor.get(c, default)
    except AttributeError:
        sf = default # when there is no suffix
    return sf


def __unset_constarint_tranform_applied_scaling_factor(c):
    """PRIVATE FUNCTION: Delete the recored scale factor that has been used
    to transofrm constraint c.  This is used when undoing a constraint
    transformation.
    """
    try:
        del c.parent_block().constaint_transformed_scaling_factor[c]
    except AttributeError:
        pass # no scaling factor suffix, is fine
    except KeyError:
        pass # no scaling factor is fine


def constraint_scaling_transform(c, s):
    """This transforms a constraint by the argument s.  The scaling factor
    applies to original constraint (e.g. if one where to call this twice in a row
    for a constraint with a scaling factor of 2, the original constraint would
    still, only be scaled by a factor of 2.)

    Args:
        c: Pyomo constraint
        s: scale factor applied to the constraint as origianlly written

    Returns:
        None
    """
    if not isinstance(c, _ConstraintData):
        raise TypeError(f"{c} is not a constraint or is an indexed constraint")
    st = get_constarint_tranform_applied_scaling_factor(c, default=1)
    v = s/st
    c.set_value(
        (__none_mult(c.lower, v), __none_mult(c.body, v), __none_mult(c.upper, v)))
    __set_constarint_tranform_applied_scaling_factor(c, s)


def constraint_scaling_transform_undo(c):
    """The undoes the scaling transforms previously applied to a constaint.

    Args:
        c: Pyomo constraint

    Returns:
        None
    """
    if not isinstance(c, _ConstraintData):
        raise TypeError(f"{c} is not a constraint or is an indexed constraint")
    v = get_constarint_tranform_applied_scaling_factor(c)
    if v is None:
        return # hasn't been transformed, so nothing to do.
    c.set_value(
        (__none_mult(c.lower, 1/v), __none_mult(c.body, 1/v), __none_mult(c.upper, 1/v)))
    __unset_constarint_tranform_applied_scaling_factor(c)


def unscaled_variables_generator(blk, descend_into=True, include_fixed=False):
    """Generator for unscaled variables

    Args:
        block

    Yields:
        variables with no scale factor
    """

    for v in blk.component_data_objects(pyo.Var, descend_into=descend_into):
        if v.fixed and not include_fixed:
            continue
        if get_scaling_factor(v) is None:
            yield v


def unscaled_constraints_generator(blk, descend_into=True):
    """Generator for unscaled constraints

    Args:
        block

    Yields:
        constraints with no scale factor
    """
    for c in blk.component_data_objects(
        pyo.Constraint, active=True, descend_into=descend_into):
        if get_scaling_factor(c) is None:
            yield c


def badly_scaled_var_generator(
    blk, large=1e4, small=1e-3, zero=1e-10, descend_into=True, include_fixed=False):
    """This provides a rough check for variables with poor scaling based on
    their current scale factors and values. For each potentially poorly scaled
    variable it returns the var and its current scaled value.

    Args:
        blk: pyomo block
        large: Magnitude that is considered to be too large
        small: Magnitude that is considered to be too small
        zero: Magnitude that is considered to be zero, variables with a value of
            zero are okay, and not reported.

    Yields:
        variable data object, current absolute value of scaled value
    """
    for v in blk.component_data_objects(pyo.Var, descend_into=descend_into):
        if v.fixed and not include_fixed:
            continue
        val = pyo.value(v, exception=False)
        if val is None:
            continue
        sf = get_scaling_factor(v, default=1)
        sv = abs(val * sf)  # scaled value
        if sv > large:
            yield v, sv
        elif sv < zero:
            continue
        elif sv < small:
            yield v, sv


def constraint_autoscale_large_jac(
    m,
    ignore_constraint_scaling=False,
    ignore_variable_scaling=False,
    max_grad=100,
    min_scale=1e-6,
    no_scale = False
):
    """Automatically scale constraints based on the Jacobian.  This function
    immitates Ipopt's default constraint scaling.  This scales constraints down
    to avoid extremely large values in the Jacobian

    Args:
        m: model to scale
        ignore_constraint_scaling: ignore existing constraint scaling
        ignore_variable_scaling: ignore existing variable scaling
        max_grad: maximum value in Jacobian after scaling, subject to minimum
            scaling factor restriction.
        min_scale: minimum scaling factor allowed, keeps constraints from being
            scaled too much.
        no_scale: just calculate the Jacobian and scaled Jacobian, don't scale
            anything
    """
    # Pynumero requires an objective, but I don't, so let's see if we have one
    n_obj = 0
    for c in m.component_data_objects(pyo.Objective, active=True):
        n_obj += 1
    # Add an objective if there isn't one
    if n_obj == 0:
        dummy_objective_name = unique_component_name(m, "objective")
        setattr(m, dummy_objective_name, pyo.Objective(expr=0))
    # Create NLP and calculate the objective
    nlp = PyomoNLP(m)
    jac = nlp.evaluate_jacobian().tocsr()
    # Get lists of varibles and constraints to translate Jacobian indexes
    clist = nlp.get_pyomo_constraints()
    vlist = nlp.get_pyomo_variables()
    # Create a scaled Jacobian to account for variable scaling, for now ignore
    # constraint scaling
    jac_scaled = jac.copy()
    for i in range(len(clist)):
        for j in jac_scaled[i].indices:
            v = vlist[j]
            if ignore_variable_scaling:
                sv = 1
            else:
                sv = get_scaling_factor(v, default=1)
            jac_scaled[i,j] = jac_scaled[i,j]/sv
    # calculate constraint scale factors
    for i in range(len(clist)):
        c = clist[i]
        sc = get_scaling_factor(c, default=1)
        if not no_scale:
            if (ignore_constraint_scaling or get_scaling_factor(c) is None):
                row = jac_scaled[i]
                for d in row.indices:
                    row[0,d] = abs(row[0,d])
                mg = row.max()
                if mg > max_grad:
                    sc = max(min_scale, max_grad/mg)
                set_scaling_factor(c, sc)
        for j in jac_scaled[i].indices:
            # update the scaled jacobian
            jac_scaled[i,j] = jac_scaled[i,j]*sc
    # delete dummy objective
    if n_obj == 0:
        delattr(m, dummy_objective_name)
    return jac, jac_scaled, nlp


# DEPRECATED functions below.

def scale_single_constraint(c):
    """This transforms a constraint with its scaling factor. If there is no
    scaling factor for the constraint, the constraint is not scaled and a
    message is logged. After transforming the constraint the scaling factor,
    scaling expression, and nomical value are all unset to ensure the constraint
    isn't scaled twice.

    Args:
        c: Pyomo constraint

    Returns:
        None
    """
    _log.warning(
        "DEPRECATED: scale_single_constraint() will be removed and has no "
        "direct replacement")
    if not isinstance(c, _ConstraintData):
        raise TypeError(
            "{} is not a constraint and cannot be the input to "
            "scale_single_constraint".format(c.name))

    v = get_scaling_factor(c)
    if v is None:
        _log.warning(
            f"{c.name} constraint has no scaling factor, so it was not scaled.")
        return
    c.set_value(
        (__none_mult(c.lower, v), __none_mult(c.body, v), __none_mult(c.upper, v)))
    unset_scaling_factor(c)


def scale_constraints(blk, descend_into=True):
    """This scales all constraints with their scaling factor suffix for a model
    or block. After scaling the constraints, the scaling factor and expression
    for each constraint is set to 1 to avoid double scaling the constraints.

    Args:
        blk: Pyomo block
        descend_into: indicates whether to descend into the other blocks on blk.
            (default = True)

    Returns:
        None
    """
    _log.warning(
        "DEPRECATED: scale_single_constraint() will be removed and has no "
        "direct replacement")
    for c in blk.component_data_objects(pyo.Constraint, descend_into=False):
        scale_single_constraint(c)
    if descend_into:
        for b in blk.component_data_objects(pyo.Block, descend_into=True):
            for c in b.component_data_objects(pyo.Constraint, descend_into=False):
                scale_single_constraint(c)
