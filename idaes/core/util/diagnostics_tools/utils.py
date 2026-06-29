#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
This module contains utility functions used for computing model diagnostics.
"""

__author__ = "Alexander Dowling, Douglas Allan, Andrew Lee, Robby Parker, Ben Knueven"


import numpy as np
from scipy.sparse.linalg import norm

from pyomo.environ import (
    value,
    Var,
)
from pyomo.common.collections import ComponentSet

from idaes.core.scaling.util import (
    get_jacobian,
    get_scaling_factor,
)


def check_parallel_jacobian(
    model,
    tolerance: float = 1e-4,
    direction: str = "row",
    jac=None,
    nlp=None,
):
    """
    Check for near-parallel rows or columns in the Jacobian.

    Near-parallel rows or columns indicate a potential degeneracy in the model,
    as this means that the associated constraints or variables are (near)
    duplicates of each other.

    For efficiency, the ``jac`` and ``nlp`` arguments may be provided if they are
    already available. If these are provided, the provided model is not used. If
    either ``jac`` or ``nlp`` is not provided, a Jacobian and ``PyomoNLP`` are
    computed using the model.

    This method is based on work published in:

    Klotz, E., Identification, Assessment, and Correction of Ill-Conditioning and
    Numerical Instability in Linear and Integer Programs, Informs 2014, pgs. 54-108
    https://pubsonline.informs.org/doi/epdf/10.1287/educ.2014.0130

    Args:
        model: model to be analysed
        tolerance: tolerance to use to determine if constraints/variables are parallel
        direction: 'row' (default, constraints) or 'column' (variables)
        jac: model Jacobian as a ``scipy.sparse.coo_matrix``, optional
        nlp: ``PyomoNLP`` of model, optional

    Returns:
        list of 2-tuples containing parallel Pyomo components

    """
    # Thanks to Robby Parker for the sparse matrix implementation and
    # significant performance improvements

    if direction not in ["row", "column"]:
        raise ValueError(
            f"Unrecognised value for direction ({direction}). "
            "Must be 'row' or 'column'."
        )

    if jac is None or nlp is None:
        jac, nlp = get_jacobian(model)

    # Get vectors that we will check, and the Pyomo components
    # they correspond to.
    if direction == "row":
        components = nlp.get_pyomo_constraints()
        csrjac = jac.tocsr()
        # Make everything a column vector (CSC) for consistency
        vectors = [csrjac[i, :].transpose().tocsc() for i in range(len(components))]
    else:  # direction == "column"
        components = nlp.get_pyomo_variables()
        cscjac = jac.tocsc()
        vectors = [cscjac[:, i] for i in range(len(components))]

    # List to store pairs of parallel components
    parallel = []

    vectors_by_nz = {}
    for vecidx, vec in enumerate(vectors):
        maxval = max(np.abs(vec.data))
        # Construct tuple of sorted col/row indices that participate
        # in this vector (with non-negligible coefficient).
        nz = tuple(
            sorted(
                idx
                for idx, val in zip(vec.indices, vec.data)
                if abs(val) > tolerance and abs(val) / maxval > tolerance
            )
        )
        if nz in vectors_by_nz:
            # Store the index as well so we know what component this
            # correrponds to.
            vectors_by_nz[nz].append((vec, vecidx))
        else:
            vectors_by_nz[nz] = [(vec, vecidx)]

    for vecs in vectors_by_nz.values():
        for idx, (u, uidx) in enumerate(vecs):
            # idx is the "local index", uidx is the "global index"
            # Frobenius norm of the matrix is 2-norm of this column vector
            unorm = norm(u, ord="fro")
            for v, vidx in vecs[idx + 1 :]:
                vnorm = norm(v, ord="fro")

                # Explicitly multiply a row vector * column vector
                prod = u.transpose().dot(v)
                absprod = abs(prod[0, 0])
                diff = abs(absprod - unorm * vnorm)
                if diff <= tolerance or diff <= tolerance * max(unorm, vnorm):
                    parallel.append((uidx, vidx))

    parallel = [(components[uidx], components[vidx]) for uidx, vidx in parallel]
    return parallel


def extreme_jacobian_entries(
    jac,
    nlp,
    large=1e4,
    small=1e-4,
    zero=1e-10,
):
    """
    Show very large and very small Jacobian entries.

    Args:
        jac: already-existing Jacobian matrix
        nlp: already-existing Pynumero NLP object from
            get_jacobian (and thus having vlist and clist
            attributes)
        large: >= to this value is considered large
        small: <= to this and >= zero is considered small
        zero: <= to this value is ignored

    Returns:
        (list of tuples), Jacobian entry, Constraint, Variable
    """
    el = []
    for i, c in enumerate(nlp.clist):
        for j in jac[i].indices:
            v = nlp.vlist[j]
            e = abs(jac[i, j])
            if (e <= small and e > zero) or e >= large:
                el.append((e, c, v))
    return el


def extreme_jacobian_rows(
    jac,
    nlp,
    large=1e4,
    small=1e-4,
):
    """
    Show very large and very small Jacobian rows. Typically indicates a badly-
    scaled constraint.

    Args:
        jac: already-existing Jacobian matrix
        nlp: already-existing Pynumero NLP object from
            get_jacobian (and thus having vlist and clist
            attributes)
        large: >= to this value is considered large
        small: <= to this is considered small

    Returns:
        (list of tuples), Row norm, Constraint
    """
    row_norms = norm(jac, ord=2, axis=1)
    # Array with values of 1 for entries with extreme row norms
    # and values of 0 otherwise
    condition_vector = np.logical_or(row_norms >= large, row_norms <= small)
    # Array of indices for which condition_vector is 1
    extreme_indices = np.nonzero(condition_vector)[0]
    return [(row_norms[k], nlp.clist[k]) for k in extreme_indices]


def extreme_jacobian_columns(
    jac,
    nlp,
    large=1e4,
    small=1e-4,
):
    """
    Show very large and very small Jacobian columns. A more reliable indicator
    of a badly-scaled variable than badly_scaled_var_generator.

    Args:
        jac: already-existing Jacobian matrix
        nlp: already-existing Pynumero NLP object from
            get_jacobian (and thus having vlist and clist
            attributes)
        large: >= to this value is considered large
        small: <= to this is considered small

    Returns:
        (list of tuples), Column norm, Variable
    """
    # Convert to csc to make iterating over columns easier
    jac = jac.tocsc()
    column_norms = norm(jac, ord=2, axis=0)
    # Array with values of 1 for entries with extreme row norms
    # and values of 0 otherwise
    condition_vector = np.logical_or(column_norms >= large, column_norms <= small)
    # Array of indices for which condition_vector is 1
    extreme_indices = np.nonzero(condition_vector)[0]
    return [(column_norms[k], nlp.vlist[k]) for k in extreme_indices]


def var_in_block(var, block):
    """
    Check if a variable is within a specific block.

    Args:
        var: The variable to check.
        block: The block to check against.

    Returns:
        True if the variable is within the block, False otherwise.
    """
    parent = var.parent_block()
    while parent is not None:
        if parent is block:
            return True
        parent = parent.parent_block()
    return False


def vars_fixed_to_zero(model):
    """
    Set of variables fixed to 0.

    Args:
        model: The model to check.

    Returns:
        A set of variables fixed to 0.
    """
    zero_vars = ComponentSet()
    for v in model.component_data_objects(Var, descend_into=True):
        if v.fixed and value(v) == 0:
            zero_vars.add(v)
    return zero_vars


def vars_near_zero(model, variable_zero_value_tolerance):
    """
    Set of variables with value near 0, as determined by the provided tolerance.

    Args:

        model: The model to check.
        variable_zero_value_tolerance: The tolerance to use when determining if a variable is near zero.
        This is applied to the scaled value of the variable, so should be chosen with scaling in mind.

    Returns:
        A set of variables with value near 0.
    """
    near_zero_vars = ComponentSet()
    for v in model.component_data_objects(Var, descend_into=True):
        sf = get_scaling_factor(v, default=1, warning=False)
        if v.value is not None and sf * abs(value(v)) <= variable_zero_value_tolerance:
            near_zero_vars.add(v)
    return near_zero_vars


def vars_violating_bounds(model, tolerance):
    """
    Set of variables with values violating their bounds by more than the provided tolerance.

    Args:
        model: The model to check.
        tolerance: The tolerance to use when determining if a variable is violating its bounds.
        This is applied to the scaled value of the variable, so should be chosen with scaling in mind.

    Returns:
        A set of variables with values violating their bounds by more than the provided tolerance.
    """
    violated_bounds = ComponentSet()
    for v in model.component_data_objects(Var, descend_into=True):
        sf = get_scaling_factor(v, default=1, warning=False)
        if v.value is not None:
            if v.lb is not None and sf * v.value <= sf * v.lb - tolerance:
                violated_bounds.add(v)
            elif v.ub is not None and sf * v.value >= sf * v.ub + tolerance:
                violated_bounds.add(v)

    return violated_bounds


def vars_with_none_value(model):
    """
    Set of variables with value None.

    Args:
        model: The model to check.

    Returns:
        A set of variables with value None.
    """
    none_value = ComponentSet()
    for v in model.component_data_objects(Var, descend_into=True):
        if v.value is None:
            none_value.add(v)

    return none_value


def vars_with_extreme_values(model, large, small, zero):
    """
    Set of variables with very large or very small values, as determined by the provided tolerances.

    Tolerances are applied to the scaled value of the variable, so should be chosen with scaling in mind.

    Args:
        model: The model to check.
        large: The threshold above which a variable is considered to have a very large value.
        small: The threshold below which a variable is considered to have a very small value.
        zero: The threshold below which a variable is considered to have a zero value and thus ignored.
    Returns:
        A set of variables with very large or very small values.
    """
    extreme_vars = ComponentSet()
    for v in model.component_data_objects(Var, descend_into=True):
        sf = get_scaling_factor(v, default=1, warning=False)
        if v.value is not None:
            mag = sf * abs(value(v))
            if mag > abs(large):
                extreme_vars.add(v)
            elif mag < abs(small) and mag > abs(zero):
                extreme_vars.add(v)

    return extreme_vars
