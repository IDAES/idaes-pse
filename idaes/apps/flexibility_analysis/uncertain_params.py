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
This module defines a function for replacing uncertain parameters with variables
"""
from typing import Sequence, Union, Mapping
from pyomo.core.base.block import _BlockData
import pyomo.environ as pe
from pyomo.core.expr.visitor import replace_expressions
from pyomo.core.base.var import _GeneralVarData
from pyomo.core.base.param import _ParamData
from .var_utils import get_all_unfixed_variables
from .indices import _VarIndex, _ConIndex


def _replace_uncertain_params(
    m: _BlockData,
    uncertain_params: Sequence[Union[_GeneralVarData, _ParamData]],
    param_nominal_values: Mapping[Union[_GeneralVarData, _ParamData], float],
    param_bounds: Mapping[Union[_GeneralVarData, _ParamData], Sequence[float]],
) -> _BlockData:
    for v in get_all_unfixed_variables(m):
        if not pe.is_constant(v._lb):  # pylint: disable=protected-access
            raise ValueError(
                f"The lower bound on {str(v)} is not constant. All variable bounds must be constant."
            )
        if not pe.is_constant(v._ub):  # pylint: disable=protected-access
            raise ValueError(
                f"The upper bound on {str(v)} is not constant. All variable bounds must be constant."
            )

    m.unc_param_vars_set = pe.Set()
    m.unc_param_vars = pe.Var(m.unc_param_vars_set)
    sub_map = dict()
    for p in uncertain_params:
        key = _VarIndex(p, None)
        m.unc_param_vars_set.add(key)
        sub_map[id(p)] = m.unc_param_vars[key]
        m.unc_param_vars[key].setlb(param_bounds[p][0])
        m.unc_param_vars[key].setub(param_bounds[p][1])
        m.unc_param_vars[key].value = param_nominal_values[p]

    m.unc_cons_set = pe.Set()
    m.unc_cons = pe.Constraint(m.unc_cons_set)
    for c in list(
        m.component_data_objects(pe.Constraint, descend_into=True, active=True)
    ):
        new_body = replace_expressions(c.body, substitution_map=sub_map)
        new_lower = replace_expressions(c.lower, substitution_map=sub_map)
        new_upper = replace_expressions(c.upper, substitution_map=sub_map)
        if c.equality:
            key = _ConIndex(c, None)
            m.unc_cons_set.add(key)
            m.unc_cons[key] = new_body == new_lower
        else:
            if c.lower is not None:
                key = _ConIndex(c, "lb")
                m.unc_cons_set.add(key)
                m.unc_cons[key] = new_lower <= new_body
            if c.upper is not None:
                key = _ConIndex(c, "ub")
                m.unc_cons_set.add(key)
                m.unc_cons[key] = new_body <= new_upper
        c.deactivate()

    return m
