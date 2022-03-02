import pyomo.environ as pe
from pyomo.core.base.block import _BlockData
from pyomo.common.collections import ComponentSet
from pyomo.core.expr.visitor import identify_variables
from coramin.utils import get_objective
from typing import Mapping, Sequence
from pyomo.core.base.var import _GeneralVarData


def get_all_unfixed_variables(m: _BlockData):
    return ComponentSet(
        v
        for v in m.component_data_objects(pe.Var, descend_into=True, active=True)
        if not v.is_fixed()
    )


def get_used_unfixed_variables(m: _BlockData):
    res = ComponentSet()
    for c in m.component_data_objects(pe.Constraint, active=True, descend_into=True):
        res.update(v for v in identify_variables(c.body, include_fixed=False))
    obj = get_objective(m)
    if obj is not None:
        res.update(identify_variables(obj.expr, include_fixed=False))
    return res


class BoundsManager(object):
    def __init__(self, m: _BlockData):
        self._vars = ComponentSet(m.component_data_objects(pe.Var, descend_into=True))
        self._saved_bounds = list()

    def save_bounds(self):
        bnds = pe.ComponentMap()
        for v in self._vars:
            bnds[v] = (v.lb, v.ub)
        self._saved_bounds.append(bnds)

    def pop_bounds(self, ndx=-1):
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
