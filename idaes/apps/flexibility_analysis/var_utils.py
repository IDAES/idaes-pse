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
"""
Some utility functions for working with pyomo variables.
"""
from typing import Mapping, Sequence, MutableSet
import pyomo.environ as pe
from pyomo.core.base.block import _BlockData
from pyomo.common.collections import ComponentSet
from pyomo.core.expr.visitor import identify_variables
from pyomo.contrib.solver.util import get_objective
from pyomo.core.base.var import _GeneralVarData


def get_all_unfixed_variables(m: _BlockData) -> MutableSet[_GeneralVarData]:
    """
    Returns a set containing all unfixed variables on the model m.

    Parameters
    ----------
    m: _BlockData
        The model for which to get the variables.

    Returns
    -------
    vset: MutableSet[_GeneralVarData]
        The set of all variables on m.
    """
    return ComponentSet(
        v
        for v in m.component_data_objects(pe.Var, descend_into=True, active=True)
        if not v.is_fixed()
    )


def get_used_unfixed_variables(m: _BlockData) -> MutableSet[_GeneralVarData]:
    """
    Returns a set containing all unfixed variables in any active constraint
    or objective on the model m.

    Parameters
    ----------
    m: _BlockData
        The model for which to get the variables.

    Returns
    -------
    vset: MutableSet[_GeneralVarData]
        The set of all variables in active constraints or objectives.
    """
    res = ComponentSet()
    for c in m.component_data_objects(pe.Constraint, active=True, descend_into=True):
        res.update(v for v in identify_variables(c.body, include_fixed=False))
    obj = get_objective(m)
    if obj is not None:
        res.update(identify_variables(obj.expr, include_fixed=False))
    return res


class BoundsManager(object):
    """
    A class for saving and restoring variable bounds.
    """

    def __init__(self, m: _BlockData):
        # TODO: maybe use get_used_unfixed_variables here?
        self._vars = ComponentSet(m.component_data_objects(pe.Var, descend_into=True))
        self._saved_bounds = list()

    def save_bounds(self):
        """
        Save the variable bounds for later use.
        """
        bnds = pe.ComponentMap()
        for v in self._vars:
            bnds[v] = (v.lb, v.ub)
        self._saved_bounds.append(bnds)

    def pop_bounds(self, ndx=-1):
        """
        Restore the variable bounds that were
        previously saved.

        Parameters
        ----------
        ndx: int
            Indicates which set of bounds to restore.
            By default, this will use the most recently
            saved bounds.
        """
        bnds = self._saved_bounds.pop(ndx)
        for v, _bnds in bnds.items():
            lb, ub = _bnds
            v.setlb(lb)
            v.setub(ub)


def _remove_var_bounds(m: _BlockData):
    for v in get_all_unfixed_variables(m):
        v.setlb(None)
        v.setub(None)
        if v.is_integer():
            raise ValueError("Unwilling to remove domain from integer variable")
        v.domain = pe.Reals


def _apply_var_bounds(bounds: Mapping[_GeneralVarData, Sequence[float]]):
    for v, (lb, ub) in bounds.items():
        if v.lb is None or v.lb < lb:
            v.setlb(lb)
        if v.ub is None or v.ub > ub:
            v.setub(ub)
