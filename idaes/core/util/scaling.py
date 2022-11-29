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

__author__ = "John Eslick, Tim Bartholomew, Robert Parker"

from math import log10
import scipy.sparse.linalg as spla
import scipy.linalg as la

import pyomo.environ as pyo
from pyomo.core.expr.visitor import identify_variables
from pyomo.network import Arc
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
from pyomo.common.modeling import unique_component_name
from pyomo.core.base.constraint import _ConstraintData
from pyomo.common.collections import ComponentMap, ComponentSet
from pyomo.dae import DerivativeVar
from pyomo.dae.flatten import slice_component_along_sets
from pyomo.util.calc_var_value import calculate_variable_from_constraint
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


def __none_left_mult(x, y):
    """PRIVATE FUNCTION, If x is None return None, else return x * y"""
    if x is not None:
        return x * y
    return None


def __scale_constraint(c, v):
    """PRIVATE FUNCTION, scale Constraint c to value v"""
    if c.equality:
        c.set_value((c.lower * v, c.body * v))
    else:
        c.set_value(
            (__none_left_mult(c.lower, v), c.body * v, __none_left_mult(c.upper, v))
        )


def scale_arc_constraints(blk):
    """Find Arc constraints in a block and its subblocks.  Then scale them based
    on the minimum scaling factor of the variables in the constraint.

    Args:
        blk: Block in which to look for Arc constraints to scale.

    Returns:
        None
    """
    for arc in blk.component_data_objects(Arc, descend_into=True):
        arc_block = arc.expanded_block
        if arc_block is None:  # arc not expanded or port empty?
            _log.warning(
                f"{arc} has no constraints. Has the Arc expansion transform "
                "been applied?"
            )
            continue
        warning = (
            "Automatic scaling for arc constraints is supported for "
            "only the Equality rule. Variable {name} on Port {port} was "
            "created with a different rule, so the corresponding constraint "
            "on {arc_name} will not be scaled."
        )
        port1 = arc.ports[0]
        port2 = arc.ports[1]
        for name in port1.vars.keys():
            if not port1.is_equality(name):
                _log.warning(
                    warning.format(name=name, port=port1.name, arc_name=arc.name)
                )
                continue
            if not port2.is_equality(name):
                _log.warning(
                    warning.format(name=name, port=port2.name, arc_name=arc.name)
                )
                continue
            con = getattr(arc_block, name + "_equality")
            for i, c in con.items():
                if i is None:
                    sf = min_scaling_factor([port1.vars[name], port2.vars[name]])
                else:
                    sf = min_scaling_factor([port1.vars[name][i], port2.vars[name][i]])
                constraint_scaling_transform(c, sf)


def map_scaling_factor(iter, default=1, warning=False, func=min, hint=None):
    """Map get_scaling_factor to an iterable of Pyomo components, and call func
    on the result.  This could be use, for example, to get the minimum or
    maximum scaling factor of a set of components.

    Args:
        iter: Iterable yielding Pyomo components
        default: The default value used when a scaling factor is missing. The
            default is default=1.
        warning: Log a warning for missing scaling factors
        func: The function to call on the resulting iterable of scaling factors.
            The default is min().
        hint: Paired with warning=True, this is a string to indicate where the
            missing scaling factor was being accessed, to easier diagnose issues.

    Returns:
        The result of func on the set of scaling factors
    """
    return func(
        map(
            lambda x: get_scaling_factor(
                x, default=default, warning=warning, hint=hint
            ),
            iter,
        )
    )


def min_scaling_factor(iter, default=1, warning=True, hint=None):
    """Map get_scaling_factor to an iterable of Pyomo components, and get the
    minimum scaling factor.

    Args:
        iter: Iterable yielding Pyomo components
        default: The default value used when a scaling factor is missing.  If
            None, this will raise an exception when scaling factors are missing.
            The default is default=1.
        warning: Log a warning for missing scaling factors
        hint: Paired with warning=True, this is a string to indicate where the
            missing scaling factor was being accessed, to easier diagnose issues.

    Returns:
        Minimum scaling factor of the components in iter
    """
    return map_scaling_factor(iter, default=default, warning=warning, func=min)


def propagate_indexed_component_scaling_factors(
    blk, typ=None, overwrite=False, descend_into=True
):
    """Use the parent component scaling factor to set all component data object
    scaling factors.

    Args:
        blk: The block on which to search for components
        typ: Component type(s) (default=(Var, Constraint, Expression, Param))
        overwrite: if a data object already has a scaling factor should it be
            overwrittten (default=False)
        descend_into: descend into child blocks (default=True)
    """
    if typ is None:
        typ = (pyo.Var, pyo.Constraint, pyo.Expression)

    for c in blk.component_objects(typ, descend_into=descend_into):
        if get_scaling_factor(c) is not None and c.is_indexed():
            for cdat in c.values():
                if overwrite or get_scaling_factor(cdat) is None:
                    set_scaling_factor(cdat, get_scaling_factor(c))


def calculate_scaling_factors(blk):
    """Look for calculate_scaling_factors methods and run them. This uses a
    recursive function to execute the subblock calculate_scaling_factors
    methods first.
    """

    def cs(blk2):
        """Recursive function for to do subblocks first"""
        for b in blk2.component_data_objects(pyo.Block, descend_into=False):
            cs(b)
        if hasattr(blk2, "calculate_scaling_factors"):
            blk2.calculate_scaling_factors()

    # Call recursive function to run calculate_scaling_factors on blocks from
    # the bottom up.
    cs(blk)
    # If a scale factor is set for an indexed component, propagate it to the
    # component data if a scale factor hasn't already been explicitly set
    propagate_indexed_component_scaling_factors(blk)
    # Use the variable scaling factors to scale the arc constraints.
    scale_arc_constraints(blk)


def set_scaling_factor(c, v, data_objects=True, overwrite=True):
    """Set a scaling factor for a model component. This function creates the
    scaling_factor suffix if needed.

    Args:
        c: component to supply scaling factor for
        v: scaling factor
        data_objects: set scaling factors for indexed data objects (default=True)
        overwrite: whether to overwrite an existing scaling factor
    Returns:
        None
    """
    if isinstance(c, (float, int)):
        # property packages can return 0 for material balance terms on components
        # doesn't exist.  This handles the case where you get a constant 0 and
        # need its scale factor to scale the mass balance.
        return 1
    try:
        suf = c.parent_block().scaling_factor

    except AttributeError:
        c.parent_block().scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
        suf = c.parent_block().scaling_factor

    if not overwrite:
        try:
            tmp = suf[c]
            # Able to access suffix value for c, so return without setting scaling factor
            return
        except KeyError:
            # No value exists yet for c, go ahead and set it
            pass
    suf[c] = v
    if data_objects and c.is_indexed():
        for cdat in c.values():
            if not overwrite:
                try:
                    tmp = suf[cdat]
                    continue
                except KeyError:
                    pass
            suf[cdat] = v


def get_scaling_factor(c, default=None, warning=False, exception=False, hint=None):
    """Get a component scale factor.

    Args:
        c: component
        default: value to return if no scale factor exists (default=None)
        warning: whether to log a warning if a scaling factor is not found
                 (default=False)
        exception: whether to raise an Exception if a scaling factor is not
                   found (default=False)
        hint: (str) a string to add to the warning or exception message to help
            locate the source.

    Returns:
        scaling factor (float)
    """
    try:
        sf = c.parent_block().scaling_factor[c]
    except (AttributeError, KeyError):
        if hint is None:
            h = ""
        else:
            h = f", {hint}"
        if warning:
            if hasattr(c, "is_component_type") and c.is_component_type():
                _log.warning(f"Missing scaling factor for {c}{h}")
            else:
                _log.warning(f"Trying to get scaling factor for unnamed expr {h}")
        if exception and default is None:
            if hasattr(c, "is_component_type") and c.is_component_type():
                _log.error(f"Missing scaling factor for {c}{h}")
            else:
                _log.error(f"Trying to get scaling factor for unnamed expr {h}")
            raise
        sf = default
    return sf


def set_and_get_scaling_factor(c, default, warning=False, exception=False, hint=None):
    """Checks whether a scaling factor exists for a component, sets the scaling factor
    if it doesn't exist, then returns the scaling factor on the component (which is
    the default value if it wasn't set originally).

    Args:
        c: component
        default: default value to use for scaling factor of c if there is none
        warning: whether to log a warning if a scaling factor is not found
                 (default=False)
        exception: whether to raise an Exception if a scaling factor is not
                   found (default=False)
        hint: (str) a string to add to the warning or exception message to help
            locate the source.

    Returns:
        scaling factor (float)
    """
    if c.is_indexed():
        raise AttributeError(
            f"Ambiguous which scaling factor to return for indexed component {c.name}."
        )
    sf = get_scaling_factor(c, warning=warning, exception=exception, hint=hint)
    if sf is None:
        sf = default
        set_scaling_factor(c, sf, data_objects=False)
    return sf


def unset_scaling_factor(c, data_objects=True):
    """Delete a component scaling factor.

    Args:
        c: component

    Returns:
        None
    """
    try:
        del c.parent_block().scaling_factor[c]
    except (AttributeError, KeyError):
        pass  # no scaling factor suffix, is fine
    try:
        if data_objects and c.is_indexed():
            for cdat in c.values():
                del cdat.parent_block().scaling_factor[cdat]
    except (AttributeError, KeyError):
        pass  # no scaling factor suffix, is fine


def populate_default_scaling_factors(c):
    """
    Method to set default scaling factors for a number of common quantities
    based of typical values expressed in SI units. Values are converted to
    those used by the property package using Pyomo's unit conversion tools.
    """
    units = c.get_metadata().derived_units

    si_scale = {
        "temperature": (100 * pyo.units.K, "temperature"),
        "pressure": (1e5 * pyo.units.Pa, "pressure"),
        "dens_mol_phase": (100 * pyo.units.mol / pyo.units.m**3, "density_mole"),
        "enth_mol": (1e4 * pyo.units.J / pyo.units.mol, "energy_mole"),
        "entr_mol": (100 * pyo.units.J / pyo.units.mol / pyo.units.K, "entropy_mole"),
        "fug_phase_comp": (1e4 * pyo.units.Pa, "pressure"),
        "fug_coeff_phase_comp": (1 * pyo.units.dimensionless, None),
        "gibbs_mol": (1e4 * pyo.units.J / pyo.units.mol, "energy_mole"),
        "mole_frac_comp": (0.001 * pyo.units.dimensionless, None),
        "mole_frac_phase_comp": (0.001 * pyo.units.dimensionless, None),
        "mw": (1e-3 * pyo.units.kg / pyo.units.mol, "molecular_weight"),
        "mw_comp": (1e-3 * pyo.units.kg / pyo.units.mol, "molecular_weight"),
        "mw_phase": (1e-3 * pyo.units.kg / pyo.units.mol, "molecular_weight"),
    }

    for p, f in si_scale.items():
        # If a default scaling factor exists, do not over write it
        if p not in c.default_scaling_factor.keys():
            if f[1] is not None:
                v = pyo.units.convert(f[0], to_units=units[f[1]])
            else:
                v = f[0]

            sf = 1 / (10 ** round(log10(pyo.value(v))))

            c.set_default_scaling(p, sf)


def __set_constraint_transform_applied_scaling_factor(c, v):
    """PRIVATE FUNCTION Set the scaling factor used to transform a constraint.
    This is used to keep track of scaling transformations that have been applied
    to constraints.

    Args:
        c: component to supply scaling factor for
        v: scaling factor
    Returns:
        None
    """
    try:
        c.parent_block().constraint_transformed_scaling_factor[c] = v
    except AttributeError:
        c.parent_block().constraint_transformed_scaling_factor = pyo.Suffix(
            direction=pyo.Suffix.LOCAL
        )
        c.parent_block().constraint_transformed_scaling_factor[c] = v


def get_constraint_transform_applied_scaling_factor(c, default=None):
    """Get a the scale factor that was used to transform a
    constraint.

    Args:
        c: constraint data object
        default: value to return if no scaling factor exists (default=None)

    Returns:
        The scaling factor that has been used to transform the constraint or the
        default.
    """
    try:
        sf = c.parent_block().constraint_transformed_scaling_factor.get(c, default)
    except AttributeError:
        sf = default  # when there is no suffix
    return sf


def __unset_constraint_transform_applied_scaling_factor(c):
    """PRIVATE FUNCTION: Delete the recorded scale factor that has been used
    to transform constraint c.  This is used when undoing a constraint
    transformation.
    """
    try:
        del c.parent_block().constraint_transformed_scaling_factor[c]
    except AttributeError:
        pass  # no scaling factor suffix, is fine
    except KeyError:
        pass  # no scaling factor is fine


def constraint_scaling_transform(c, s, overwrite=True):
    """This transforms a constraint by the argument s.  The scaling factor
    applies to original constraint (e.g. if one where to call this twice in a row
    for a constraint with a scaling factor of 2, the original constraint would
    still, only be scaled by a factor of 2.)

    Args:
        c: Pyomo constraint
        s: scale factor applied to the constraint as originally written
        overwrite: overwrite existing scaling factors if present (default=True)

    Returns:
        None
    """
    if not isinstance(c, _ConstraintData):
        raise TypeError(f"{c} is not a constraint or is an indexed constraint")
    st = get_constraint_transform_applied_scaling_factor(c, default=None)

    if not overwrite and st is not None:
        # Existing scaling factor and overwrite False, do nothing
        return

    if st is None:
        # If no existing scaling factor, use value of 1
        st = 1

    v = s / st
    __scale_constraint(c, v)
    __set_constraint_transform_applied_scaling_factor(c, s)


def constraint_scaling_transform_undo(c):
    """The undoes the scaling transforms previously applied to a constraint.

    Args:
        c: Pyomo constraint

    Returns:
        None
    """
    if not isinstance(c, _ConstraintData):
        raise TypeError(f"{c} is not a constraint or is an indexed constraint")
    v = get_constraint_transform_applied_scaling_factor(c)
    if v is None:
        return  # hasn't been transformed, so nothing to do.
    __scale_constraint(c, 1 / v)
    __unset_constraint_transform_applied_scaling_factor(c)


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
        pyo.Constraint, active=True, descend_into=descend_into
    ):
        if (
            get_scaling_factor(c) is None
            and get_constraint_transform_applied_scaling_factor(c) is None
        ):
            yield c


def constraints_with_scale_factor_generator(blk, descend_into=True):
    """Generator for constraints scaled by a scaling factor, may or not have
    been transformed.

    Args:
        block

    Yields:
        constraint with a scale factor, scale factor
    """
    for c in blk.component_data_objects(
        pyo.Constraint, active=True, descend_into=descend_into
    ):
        s = get_scaling_factor(c)
        if s is not None:
            yield c, s


def badly_scaled_var_generator(
    blk, large=1e4, small=1e-3, zero=1e-10, descend_into=True, include_fixed=False
):
    """This provides a rough check for variables with poor scaling based on
    their current scale factors and values. For each potentially poorly scaled
    variable it returns the var and its current scaled value.

    Note that while this method is a reasonable heuristic for nonnegative
    variables like (absolute) temperature and pressure, molar flows, etc., it
    can be misleading for variables like enthalpies and fluxes.

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
    no_scale=False,
    equality_constraints_only=False,
):
    """Automatically scale constraints based on the Jacobian.  This function
    imitates Ipopt's default constraint scaling.  This scales constraints down
    to avoid extremely large values in the Jacobian.  This function also returns
    the unscaled and scaled Jacobian matrixes and the Pynumero NLP which can be
    used to identify the constraints and variables corresponding to the rows and
    comlumns.

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
        equality_constraints_only: Include only the equality constraints in the
            Jacobian

    Returns:
        unscaled Jacobian CSR from, scaled Jacobian CSR from, Pynumero NLP
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
    if equality_constraints_only:
        jac = nlp.evaluate_jacobian_eq().tocsr()
    else:
        jac = nlp.evaluate_jacobian().tocsr()
    # Get lists of varibles and constraints to translate Jacobian indexes
    # save them on the NLP for later, since genrating them seems to take a while
    if equality_constraints_only:
        nlp.clist = clist = nlp.get_pyomo_equality_constraints()
    else:
        nlp.clist = clist = nlp.get_pyomo_constraints()
    nlp.vlist = vlist = nlp.get_pyomo_variables()
    # Create a scaled Jacobian to account for variable scaling, for now ignore
    # constraint scaling
    jac_scaled = jac.copy()
    for i, c in enumerate(clist):
        for j in jac_scaled[i].indices:
            v = vlist[j]
            if ignore_variable_scaling:
                sv = 1
            else:
                sv = get_scaling_factor(v, default=1)
            jac_scaled[i, j] = jac_scaled[i, j] / sv
    # calculate constraint scale factors
    for i, c in enumerate(clist):
        sc = get_scaling_factor(c, default=1)
        if not no_scale:
            if ignore_constraint_scaling or get_scaling_factor(c) is None:
                sc = 1
                row = jac_scaled[i]
                for d in row.indices:
                    row[0, d] = abs(row[0, d])
                mg = row.max()
                if mg > max_grad:
                    sc = max(min_scale, max_grad / mg)
                set_scaling_factor(c, sc)
        for j in jac_scaled[i].indices:
            # update the scaled jacobian
            jac_scaled[i, j] = jac_scaled[i, j] * sc
    # delete dummy objective
    if n_obj == 0:
        delattr(m, dummy_objective_name)
    return jac, jac_scaled, nlp


def get_jacobian(m, scaled=True, equality_constraints_only=False):
    """
    Get the Jacobian matrix at the current model values. This function also
    returns the Pynumero NLP which can be used to identify the constraints and
    variables corresponding to the rows and columns.

    Args:
        m: model to get Jacobian from
        scaled: if True return scaled Jacobian, else get unscaled
        equality_constraints_only: Only include equality constraints in the
            Jacobian calculated and scaled

    Returns:
        Jacobian matrix in Scipy CSR format, Pynumero nlp
    """
    jac, jac_scaled, nlp = constraint_autoscale_large_jac(
        m, no_scale=True, equality_constraints_only=equality_constraints_only
    )
    if scaled:
        return jac_scaled, nlp
    else:
        return jac, nlp


def extreme_jacobian_entries(
    m=None, scaled=True, large=1e4, small=1e-4, zero=1e-10, jac=None, nlp=None
):
    """
    Show very large and very small Jacobian entries.

    Args:
        m: model
        scaled: if true use scaled Jacobian
        large: >= to this value is considered large
        small: <= to this and >= zero is considered small

    Returns:
        (list of tuples), Jacobian entry, Constraint, Variable
    """
    if jac is None or nlp is None:
        jac, nlp = get_jacobian(m, scaled)
    el = []
    for i, c in enumerate(nlp.clist):
        for j in jac[i].indices:
            v = nlp.vlist[j]
            e = abs(jac[i, j])
            if (e <= small and e > zero) or e >= large:
                el.append((e, c, v))
    return el


def extreme_jacobian_rows(
    m=None, scaled=True, large=1e4, small=1e-4, jac=None, nlp=None
):
    """
    Show very large and very small Jacobian rows. Typically indicates a badly-
    scaled constraint.

    Args:
        m: model
        scaled: if true use scaled Jacobian
        large: >= to this value is considered large
        small: <= to this is considered small

    Returns:
        (list of tuples), Row norm, Constraint
    """
    # Need both jac for the linear algebra and nlp for constraint names
    if jac is None or nlp is None:
        jac, nlp = get_jacobian(m, scaled)
    el = []
    for i, c in enumerate(nlp.clist):
        norm = 0
        # Calculate L2 norm
        for j in jac[i].indices:
            norm += jac[i, j] ** 2
        norm = norm**0.5
        if norm <= small or norm >= large:
            el.append((norm, c))
    return el


def extreme_jacobian_columns(
    m=None, scaled=True, large=1e4, small=1e-4, jac=None, nlp=None
):
    """
    Show very large and very small Jacobian columns. A more reliable indicator
    of a badly-scaled variable than badly_scaled_var_generator.

    Args:
        m: model
        scaled: if true use scaled Jacobian
        large: >= to this value is considered large
        small: <= to this is considered small

    Returns:
        (list of tuples), Column norm, Variable
    """
    # Need both jac for the linear algebra and nlp for variable names
    if jac is None or nlp is None:
        jac, nlp = get_jacobian(m, scaled)
    jac = jac.tocsc()
    el = []
    for j, v in enumerate(nlp.vlist):
        norm = 0
        # Calculate L2 norm
        for i in jac.getcol(j).indices:
            norm += jac[i, j] ** 2
        norm = norm**0.5
        if norm <= small or norm >= large:
            el.append((norm, v))
    return el


def jacobian_cond(m=None, scaled=True, ord=None, pinv=False, jac=None):
    """
    Get the condition number of the scaled or unscaled Jacobian matrix of a model.

    Args:
        m: calculate the condition number of the Jacobian from this model.
        scaled: if True use scaled Jacobian, else use unscaled
        ord: norm order, None = Frobenius, see scipy.sparse.linalg.norm for more
        pinv: Use pseudoinverse, works for non-square matrices
        jac: (optional) previously calculated Jacobian

    Returns:
        (float) Condition number
    """
    if jac is None:
        jac, nlp = get_jacobian(m, scaled)
    jac = jac.tocsc()
    if jac.shape[0] != jac.shape[1] and not pinv:
        _log.warning("Nonsquare Jacobian using pseudo inverse")
        pinv = True
    if not pinv:
        jac_inv = spla.inv(jac)
        return spla.norm(jac, ord) * spla.norm(jac_inv, ord)
    else:
        jac_inv = la.pinv(jac.toarray())
        return spla.norm(jac, ord) * la.norm(jac_inv, ord)


def scale_time_discretization_equations(blk, time_set, time_scaling_factor):
    """
    Scales time discretization equations generated via a Pyomo discretization
    transformation. Also scales continuity equations for collocation methods
    of discretization that require them.

    Args:
        blk: Block whose time discretization equations are being scaled
        time_set: Time set object. For an IDAES flowsheet object fs, this is fs.time.
        time_scaling_factor: Scaling factor to use for time

    Returns:
        None
    """

    tname = time_set.local_name

    # Copy and pasted from solvers.petsc.find_discretization_equations then modified
    for var in blk.component_objects(pyo.Var):
        if isinstance(var, DerivativeVar):
            cont_set_set = ComponentSet(var.get_continuousset_list())
            if time_set in cont_set_set:
                if len(cont_set_set) > 1:
                    _log.warning(
                        "IDAES presently does not support automatically scaling discretization equations for "
                        f"second order or higher derivatives like {var.name} that are differentiated at least once with "
                        "respect to time. Please scale the corresponding discretization equation yourself."
                    )
                    continue
                state_var = var.get_state_var()
                parent_block = var.parent_block()

                disc_eq = getattr(parent_block, var.local_name + "_disc_eq")
                # Look for continuity equation, which exists only for collocation with certain sets of polynomials
                try:
                    cont_eq = getattr(
                        parent_block, state_var.local_name + "_" + tname + "_cont_eq"
                    )
                except AttributeError:
                    cont_eq = None

                deriv_dict = dict(
                    (key, pyo.Reference(slc))
                    for key, slc in slice_component_along_sets(var, (time_set,))
                )
                state_dict = dict(
                    (key, pyo.Reference(slc))
                    for key, slc in slice_component_along_sets(state_var, (time_set,))
                )
                disc_dict = dict(
                    (key, pyo.Reference(slc))
                    for key, slc in slice_component_along_sets(disc_eq, (time_set,))
                )
                if cont_eq is not None:
                    cont_dict = dict(
                        (key, pyo.Reference(slc))
                        for key, slc in slice_component_along_sets(cont_eq, (time_set,))
                    )
                for key, deriv in deriv_dict.items():
                    state = state_dict[key]
                    disc = disc_dict[key]
                    if cont_eq is not None:
                        cont = cont_dict[key]
                    for t in time_set:
                        s_state = get_scaling_factor(state[t], default=1, warning=True)
                        set_scaling_factor(
                            deriv[t], s_state / time_scaling_factor, overwrite=False
                        )
                        s_deriv = get_scaling_factor(deriv[t])
                        # Check time index to decide what constraints to scale
                        if cont_eq is None:
                            if t == time_set.first() or t == time_set.last():
                                try:
                                    constraint_scaling_transform(
                                        disc[t], s_deriv, overwrite=False
                                    )
                                except KeyError:
                                    # Discretization and continuity equations may or may not exist at the first or last time
                                    # points depending on the method. Backwards skips first, forwards skips last, central skips
                                    # both (which means the user needs to provide additional equations)
                                    pass
                            else:
                                constraint_scaling_transform(
                                    disc[t], s_deriv, overwrite=False
                                )
                        else:
                            # Lagrange-Legendre is a pain, because it has continuity equations on the edges of finite
                            # instead of discretization equations, but no intermediate continuity equations, so we have
                            # to look for both at every timepoint
                            try:
                                constraint_scaling_transform(
                                    disc[t], s_deriv, overwrite=False
                                )
                            except KeyError:
                                if t != time_set.first():
                                    constraint_scaling_transform(
                                        cont[t], s_state, overwrite=False
                                    )


class CacheVars(object):
    """
    A class for saving the values of variables then reloading them,
    usually after they have been used to perform some solve or calculation.
    """

    def __init__(self, vardata_list):
        self.vars = vardata_list
        self.cache = [None for var in self.vars]

    def __enter__(self):
        for i, var in enumerate(self.vars):
            self.cache[i] = var.value
        return self

    def __exit__(self, ex_type, ex_value, ex_traceback):
        for i, var in enumerate(self.vars):
            var.set_value(self.cache[i])


class FlattenedScalingAssignment(object):
    """
    A class to assist in the calculation of scaling factors when a
    variable-constraint assignment can be constructed, especially when
    the variables and constraints are all indexed by some common set(s).
    """

    def __init__(self, scaling_factor, varconlist=None, nominal_index=()):
        """
        Args:
            scaling_factor: A Pyomo scaling_factor Suffix that will hold all
                            the scaling factors calculated
            varconlist: A list of variable, constraint tuples. These variables
                        and constraints should be indexed by the same sets,
                        so they may need to be references-to-slices along some
                        common sets.
            nominal_index: The index of variables and constraints to access
                           when a calculation needs to be performed using
                           data objects.
        """
        if varconlist is None:
            varconlist = []

        self.scaling_factor = scaling_factor
        self.nominal_index = nominal_index
        if nominal_index is None or nominal_index == ():
            self.dim = 0
        else:
            try:
                self.dim = len(nominal_index)
            except TypeError:
                self.dim = 1

        varlist = []
        conlist = []
        for var, con in varconlist:
            varlist.append(var)
            conlist.append(con)
        self.varlist = varlist
        self.conlist = conlist

        data_getter = self.get_representative_data_object
        var_con_data_list = [
            (data_getter(var), data_getter(con)) for var, con in varconlist
        ]
        con_var_data_list = [
            (data_getter(con), data_getter(var)) for var, con in varconlist
        ]
        self.var2con = ComponentMap(var_con_data_list)
        self.con2var = ComponentMap(con_var_data_list)

    def get_representative_data_object(self, obj):
        """
        Gets a data object from an object of the appropriate dimension
        """
        if self.dim == 0:
            # In this way, obj can be a data object and this class can be
            # used even if the assignment is not between "flattened components"
            return obj
        else:
            nominal_index = self.nominal_index
            return obj[nominal_index]

    def calculate_variable_scaling_factor(self, var, include_fixed=False):
        """
        Calculates the scaling factor of a variable based on the
        constraint assigned to it. Loads each variable in that constraint
        with its nominal value (inverse of scaling factor), calculates
        the value of the target variable from the constraint, then sets
        its scaling factor to the inverse of the calculated value.
        """
        vardata = self.get_representative_data_object(var)
        condata = self.var2con[vardata]
        scaling_factor = self.scaling_factor

        in_constraint = list(
            identify_variables(
                condata.expr,
                include_fixed=include_fixed,
            )
        )
        source_vars = [v for v in in_constraint if v is not vardata]
        nominal_source = [1 / scaling_factor[var] for var in source_vars]

        with CacheVars(in_constraint) as cache:
            for v, nom_val in zip(source_vars, nominal_source):
                v.set_value(nom_val)
            # This assumes that target var is initialized to a somewhat
            # reasonable value
            calculate_variable_from_constraint(vardata, condata)
            nominal_target = vardata.value
        if nominal_target == 0:
            target_factor = 1.0
        else:
            target_factor = abs(1 / nominal_target)

        if self.dim == 0:
            scaling_factor[var] = target_factor
        else:
            for v in var.values():
                scaling_factor[v] = target_factor

    def set_constraint_scaling_factor(self, con):
        """
        Sets the scaling factor of a constraint to that of its assigned variable
        """
        condata = self.get_representative_data_object(con)
        vardata = self.con2var[condata]
        scaling_factor = self.scaling_factor

        var_factor = scaling_factor[vardata]
        if self.dim == 0:
            scaling_factor[con] = var_factor
        else:
            for c in con.values():
                scaling_factor[c] = var_factor

    def set_derivative_factor_from_state(self, deriv, nominal_wrt=1.0):
        """
        Sets the scaling factor for a DerivativeVar equal to the factor for
        its state var at every index. This method needs access to the
        get_state_var method, so deriv must be an actual DerivativeVar,
        not a reference-to-slice.
        """
        scaling_factor = self.scaling_factor
        state_var = deriv.get_state_var()
        for index, dv in deriv.items():
            state_data = state_var[index]
            nominal_state = 1 / scaling_factor[state_data]
            nominal_deriv = nominal_state / nominal_wrt
            scaling_factor[dv] = 1 / nominal_deriv
