# -*- coding: utf-8 -*-
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
""" Functions for categorizing variable into e.g. differential and algebraic.
"""

from pyomo.environ import (
        Reference,
        )
from pyomo.dae import DerivativeVar
from pyomo.common.collections import ComponentSet, ComponentMap

import idaes.apps.caprese.nmpc_var as nmpc_var
from idaes.apps.caprese.common.config import VariableCategory

CATEGORY_TYPE_MAP = {
        VariableCategory.DIFFERENTIAL: nmpc_var.DiffVar,
        VariableCategory.ALGEBRAIC: nmpc_var.AlgVar,
        VariableCategory.DERIVATIVE: nmpc_var.DerivVar,
        VariableCategory.INPUT: nmpc_var.InputVar,
        VariableCategory.FIXED: nmpc_var.FixedVar,
        VariableCategory.MEASUREMENT: nmpc_var.MeasuredVar,
        }

def categorize_dae_variables(dae_vars, time, inputs, measurements=None):
    t0 = time.first()
    t1 = time.get_finite_elements()[1]
    deriv_vars = []
    diff_vars = []
    input_vars = []
    alg_vars = []
    fixed_vars = []

    # TODO: give user ability to specify measurements and disturbances
    measured_vars = []

    if measurements is not None:
        infer_measurements = False
        user_measurements = ComponentSet(measurements)
        updated_user_measurements = ComponentSet(measurements)
        user_measured_vars = []
    else:
        infer_measurements = True
        updated_user_measurements = ComponentSet()

    dae_map = ComponentMap([(v[t0], v) for v in dae_vars])
    t0_vardata = list(dae_map.keys())

    if inputs is None:
        inputs = []
    input_set = ComponentSet(inputs)
    updated_input_set = ComponentSet(inputs)

    for var0 in t0_vardata:
        if var0 in input_set:
            updated_input_set.remove(var0)
            time_slice = dae_map.pop(var0)
            input_vars.append(time_slice)

        if var0 in updated_user_measurements:
            updated_user_measurements.remove(var0)
            # Don't pop measured vars. They will be popped elsewhere.
            time_slice = dae_map[var0]
            user_measured_vars.append(time_slice)

        parent = var0.parent_component()
        if not isinstance(parent, DerivativeVar):
            continue
        if not time in ComponentSet(parent.get_continuousset_list()):
            continue
        index0 = var0.index()
        var1 = dae_map[var0][t1]
        index1 = var1.index()
        state = parent.get_state_var()

        if state[index1].fixed:
            # Assume state var is fixed everywhere, so derivative
            # 'isn't really' a derivative.
            # Should be safe to remove state from dae_map here
            state_slice = dae_map.pop(state[index0])
            fixed_vars.append(state_slice)
            continue
        if state[index0] in input_set:
            # If differential variable is an input, then this DerivativeVar
            # is 'not really a derivative'
            continue

        deriv_slice = dae_map.pop(var0)

        if var1.fixed:
            # Assume derivative has been fixed everywhere.
            # Add to list of fixed variables, and don't remove its state variable.
            fixed_vars.append(deriv_slice)
        elif var0.fixed:
            # In this case the derivative has been used as an initial condition. 
            # Still want to include it in the list of derivatives.
            measured_vars.append(deriv_slice)
            state_slice = dae_map.pop(state[index0])
            if state[index0].fixed:
                measured_vars.append(state_slice)
            deriv_vars.append(deriv_slice)
            diff_vars.append(state_slice)
        else:
            # Neither is fixed. This should be the most common case.
            state_slice = dae_map.pop(state[index0])
            if state[index0].fixed:
                measured_vars.append(state_slice)
            deriv_vars.append(deriv_slice)
            diff_vars.append(state_slice)

    if updated_input_set:
        raise RuntimeError('Not all inputs could be found')
    assert len(deriv_vars) == len(diff_vars)

    for var0, time_slice in dae_map.items():
        var1 = time_slice[t1]
        # If the variable is still in the list of time-indexed vars,
        # it must either be fixed (not a var) or be an algebraic var
        if var1.fixed:
            fixed_vars.append(time_slice)
        else:
            if var0.fixed:
                measured_vars.append(time_slice)
            alg_vars.append(time_slice)

    category_list_map = {
            VariableCategory.DERIVATIVE: deriv_vars,
            VariableCategory.DIFFERENTIAL: diff_vars,
            VariableCategory.ALGEBRAIC: alg_vars,
            VariableCategory.INPUT: input_vars,
            VariableCategory.FIXED: fixed_vars,
            VariableCategory.MEASUREMENT: measured_vars,
            }
    if measurements is not None:
        # If the user provided their own measurements,
        # override the inferred measurements. Assume the user
        # will modify the state of their variables appropriately.
        category_list_map[VariableCategory.MEASUREMENT] = user_measured_vars
    category_dict = {
            category: [
                Reference(ref.referent, ctype=ctype)
                for ref in category_list_map[category]
                ]
            for category, ctype in CATEGORY_TYPE_MAP.items()
            }
    return category_dict

