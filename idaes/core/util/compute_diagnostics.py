# -*- coding: utf-8 -*-
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
Compute diagnostics values and return as Pydantic data models.

Usage:

    from idaes.core.util.compute_diagnostics import DiagnosticsData, DiagnosticsToolbox
    import json

    model = build_and_run_model()  # replace with real code
    toolbox = DiagnosticsToolbox(model)
    dd = DiagnosticsData(toolbox)
    # collect data
    data = {
        "variable_issues": dd.variables(),
        "structural_issues": dd.structural_issues()
    }
    # pretty-print it
    print(json.dumps(data, indent=2))

"""
__author__ = "Dan Gunter (LBNL)"

# stdlib
from enum import StrEnum
from math import log

# third party
from pydantic import BaseModel, Field
from pyomo.environ import (
    # Binary,
    # Integers,
    # Block,
    # check_optimal_termination,
    # ComponentMap,
    # ConcreteModel,
    # Constraint,
    # Expression,
    # Objective,
    # Param,
    # RangeSet,
    # Set,
    # SolverFactory,
    value,
    # Var,
)
from pyomo.util.check_units import identify_inconsistent_units

# package
from idaes.core.scaling.util import (
    get_jacobian,
)
from idaes.core.util.model_statistics import (
    degrees_of_freedom,
    # activated_blocks_set,
    # deactivated_blocks_set,
    # activated_equalities_set,
    # deactivated_equalities_set,
    # activated_inequalities_set,
    # deactivated_inequalities_set,
    # activated_objectives_set,
    # deactivated_objectives_set,
    variables_in_activated_constraints_set,
    variables_not_in_activated_constraints_set,
    variables_with_none_value_in_activated_equalities_set,
    # number_activated_greybox_equalities,
    # number_deactivated_greybox_equalities,
    # activated_greybox_block_set,
    # deactivated_greybox_block_set,
    # greybox_block_set,
    # unfixed_greybox_variables,
    # greybox_variables,
    # large_residuals_set,
    variables_near_bounds_set,
)
from .model_diagnostics import (
    DiagnosticsToolbox,
    _extreme_jacobian_columns,
    _var_in_block,
    _vars_fixed_to_zero,
    _vars_near_zero,
    _vars_violating_bounds,
    _vars_with_extreme_values,
    #    _vars_with_none_value
)

# -------------------------------------------------------------------------------
# Pydantic models
# -------------------------------------------------------------------------------


def _block_list_names(blocks) -> list[str]:
    """Retrieve names of blocks, including indexed blocks, into a list."""
    b_items = []
    for b in blocks:
        if hasattr(b, "name"):
            b_items.append(b.name)
        else:  # indexed
            for i in b:
                b_items.append(i.name)
    return b_items


class VariableListData(BaseModel):
    """List of variables that satisfy some condition."""

    #: short description of condition these variables satisfy
    tag: str
    #: longer description of condition these variables satisfy
    description: str
    #: names of variables
    variables: list[str] = Field(default_factory=list)
    #: optional values of variables
    values: list[float | None] = Field(default_factory=list)
    value_format: str = ".6E"
    #: optional descriptive details for variables
    details: list[str] = Field(default_factory=list)


class VCSet(BaseModel):
    """Combined variables and constraints.

    Not returned directly; used by `StructuralWarningsData`.
    """

    variables: list[str]
    constraints: list[str]

    @classmethod
    def from_blocks(cls, var_blocks, const_blocks) -> "VCSet":
        v_items = _block_list_names(var_blocks)
        c_items = _block_list_names(const_blocks)
        return VCSet(variables=v_items, constraints=c_items)


class EvalErrorData(BaseModel):
    component_name: str
    message: str


class StructuralWarningsData(BaseModel):
    """Structural warnings.

    All possibilities are listed, the value will be None if it is not
    an issue for this model.
    """

    dof: int | None = None
    inconsistent_units: list[str] | None = None
    underconstrained_set: VCSet | None = None
    overconstrained_set: VCSet | None = None
    evaluation_errors: list[EvalErrorData] | None = None


class StructuralCautionsData(BaseModel):
    """Structural cautions.

    All possibilities are listed, the value will be None if it is not
    an issue for this model.
    """

    zero_vars: list[str] | None = None
    unused_vars_free: list[str] | None = None
    unused_vars_fixed: list[str] | None = None


class StructuralIssuesData(BaseModel):
    """Structural issues: warnings and cautions."""

    warnings: StructuralWarningsData
    cautions: StructuralCautionsData


# -------------------------------------------------------------------------------
# Interface
# -------------------------------------------------------------------------------


class VariableCondition(StrEnum):
    external = "are external variables that appear in constraints"
    unused = "do not appear in any activated constraints"
    fixed_to_zero = "are fixed to zero"
    at_or_outside_bounds = "have values that fall at or outside their bounds"
    with_none_value = "have a value of none"
    value_near_zero = "have a value near zero"
    extreme_values = "have extreme values"
    near_bounds = "have values close to their bounds"
    extreme_jacobians = "corresponding to Jacobian columns with extreme norms"


class DiagnosticsData:
    """Interface to get diagnostics data"""

    VC = VariableCondition  # alias

    def __init__(self, toolbox: DiagnosticsToolbox):
        self._toolbox = toolbox

    def variables(
        self, conditions: list[VariableCondition] | None = None
    ) -> list[VariableListData]:
        """Compute the list of variables meeting some condition and return as Pydantic object.

        Args:
            conditions: Zero or more conditions. If zero, look for all
                  conditions. If one, just return one. If multiple, return one per condition.

        Returns:
            Selected variables and associated metadata, as a list of length 1 or greater
        """
        results = []
        if conditions is None or len(conditions) == 0:
            conditions = [v for v in VariableCondition]
        for cond in conditions:
            results.append(self._get_variables_for_condition(cond))
        return results

    def structural_issues(
        self, evaluation_errors=True, unit_consistency=True
    ) -> StructuralIssuesData:
        """Compute structural warnings and cautions.

        Args:
            evaluation_errors: Include potential evaluation errors
            unit_consistency: Include unit consistency checks

        Returns:
            Pydantic data model representing found issues
        """
        tbx, model = self._toolbox, self._toolbox.model
        uc = [] if unit_consistency else identify_inconsistent_units(model)
        uc_var, uc_con, oc_var, oc_con = tbx.get_dulmage_mendelsohn_partition()
        w, c = StructuralWarningsData(), StructuralCautionsData()

        # Warnings

        dof = degrees_of_freedom(model)
        if dof != 0:
            w.dof = dof

        if len(uc) > 0:
            w.inconsistent_units = uc
        if len(uc_var) + len(uc_con) > 0:
            uc_set = set(uc_var) + set(uc_con)
            w.underconstrained_set = VCSet.from_blocks(uc_var, uc_con)
        if len(oc_var) + len(oc_con) > 0:
            w.overconstrained_set = VCSet.from_blocks(oc_var, oc_con)

        if not evaluation_errors:
            eval_warnings = tbx._collect_potential_eval_errors()
            if len(eval_warnings) > 0:
                w.evaluation_errors = []
                for ew_raw in eval_warnings:
                    ew_comp, ew_msg = ew_raw.split(":")
                    ee = EvalErrorData(
                        component_name=ew_comp.strip(), message=ew_msg.strip()
                    )
                    w.evaluation_errors.append(ee)

        # Cautions

        zero_vars = _vars_fixed_to_zero(model)
        if len(zero_vars) > 0:
            c.zero_vars = _block_list_names(zero_vars)

        unused_vars = variables_not_in_activated_constraints_set(model)
        if len(unused_vars) > 0:
            uv_free, uv_fixed = [], []
            for v in unused_vars:
                if v.fixed:
                    uv_fixed.append(v)
                else:
                    uv_free.append(v)
            if uv_fixed:
                c.unused_vars_fixed = _block_list_names(uv_fixed)
            if uv_free:
                c.unused_vars_free = _block_list_names(uv_free)

        return StructuralIssuesData(warnings=w, cautions=c)

    def _get_variables_for_condition(self, cond: VariableCondition) -> VariableListData:
        tbx = self._toolbox  # alias
        kwargs = {}  # additional kw for VariableListData
        desc = str(cond)  # default description
        details, values = None, None
        if cond == VariableCondition.external:
            cvars = variables_in_activated_constraints_set(tbx._model)
            variables = [v.name for v in cvars if not _var_in_block(v, tbx._model)]
        elif cond == VariableCondition.unused:
            variables = [
                str(v) for v in variables_not_in_activated_constraints_set(tbx._model)
            ]
        elif cond == VariableCondition.fixed_to_zero:
            variables = [str(v) for v in _vars_fixed_to_zero(tbx._model)]
        elif cond == VariableCondition.at_or_outside_bounds:
            variables, details = [], []
            for v in _vars_violating_bounds(
                tbx._model,
                tolerance=tbx.config.variable_bounds_violation_tolerance,
            ):
                variables.append(f"{v.name} ({'fixed' if v.fixed else 'free'})")
                details.append(f"bounds={v.bounds}")
        elif cond == VariableCondition.with_none_value:
            variables = [
                v.name
                for v in variables_with_none_value_in_activated_equalities_set(
                    tbx._model
                )
            ]
        elif cond == VariableCondition.value_near_zero:
            variables, values = [], []
            for v in _vars_near_zero(
                tbx._model, tbx.config.variable_zero_value_tolerance
            ):
                variables.append(v.name)
                values.append(value(v))
        elif cond == VariableCondition.extreme_values:
            desc += f" (<{tbx.config.variable_small_value_tolerance:.1E} or "
            desc += f"> {tbx.config.variable_large_value_tolerance:.1E})"
            variables, values = [], []
            for v in _vars_with_extreme_values(
                model=tbx._model,
                large=tbx.config.variable_large_value_tolerance,
                small=tbx.config.variable_small_value_tolerance,
                zero=tbx.config.variable_zero_value_tolerance,
            ):
                variables.append(v.name)
                values.append(value(v))
        elif cond == VariableCondition.near_bounds:
            desc += f" (abs={tbx.config.variable_bounds_absolute_tolerance:.1E}, "
            desc += f"rel={tbx.config.variable_bounds_relative_tolerance:.1E})"
            variables, values = [], []
            for v in variables_near_bounds_set(
                tbx._model,
                abs_tol=tbx.config.variable_bounds_absolute_tolerance,
                rel_tol=tbx.config.variable_bounds_relative_tolerance,
            ):
                variables.append(v.name)
                values.append(value(v))
        elif cond == VariableCondition.extreme_jacobians:
            tbx._verify_active_variables_initialized()
            desc += f" (<{tbx.config.jacobian_small_value_caution:.1E} or "
            desc += f">{tbx.config.jacobian_large_value_caution:.1E})"
            # compute the extreme jacobians
            jac, nlp = get_jacobian(tbx._model)
            xjc = _extreme_jacobian_columns(
                jac=jac,
                nlp=nlp,
                large=tbx.config.jacobian_large_value_caution,
                small=tbx.config.jacobian_small_value_caution,
            )
            xjc.sort(key=lambda i: abs(log(i[0])), reverse=True)
            # place in output object
            variables, values = [], []
            for v in xjc:
                variables.append(v[1].name)
                values.append(v[0])
            kwargs["value_format"] = ".3E"

        # Fill out optional parts, if missing
        if details is None:
            details = [""] * len(variables)
        if values is None:
            values = [None] * len(variables)

        # return as Pydantic data object
        return VariableListData(
            tag=cond.value,
            description=desc,
            variables=variables,
            details=details,
            values=values,
            **kwargs,
        )
