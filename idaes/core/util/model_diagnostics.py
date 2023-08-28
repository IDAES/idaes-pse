# -*- coding: utf-8 -*-
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
This module contains a collection of tools for diagnosing modeling issues.
"""

__author__ = "Alexander Dowling, Douglas Allan, Andrew Lee"

from operator import itemgetter
from sys import stdout
from math import log

import numpy as np
from scipy.linalg import svd
from scipy.sparse.linalg import svds, norm
from scipy.sparse import issparse, find

from pyomo.environ import (
    Binary,
    Block,
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Objective,
    Param,
    Set,
    SolverFactory,
    value,
    Var,
)
from pyomo.core.base.block import _BlockData
from pyomo.common.collections import ComponentSet
from pyomo.common.config import ConfigDict, ConfigValue, document_kwargs_from_configdict
from pyomo.util.check_units import identify_inconsistent_units
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.core.expr.visitor import identify_variables
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
from pyomo.contrib.pynumero.asl import AmplInterface
from pyomo.common.deprecation import deprecation_warning

from idaes.core.util.model_statistics import (
    activated_blocks_set,
    deactivated_blocks_set,
    activated_equalities_set,
    deactivated_equalities_set,
    activated_inequalities_set,
    deactivated_inequalities_set,
    activated_objectives_set,
    deactivated_objectives_set,
    variables_in_activated_constraints_set,
    variables_not_in_activated_constraints_set,
    degrees_of_freedom,
    large_residuals_set,
    variables_near_bounds_set,
)
from idaes.core.util.scaling import (
    get_jacobian,
    extreme_jacobian_columns,
    extreme_jacobian_rows,
    extreme_jacobian_entries,
    jacobian_cond,
)
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


MAX_STR_LENGTH = 84
TAB = " " * 4

# TODO: Add suggested steps to cautions - how?

CONFIG = ConfigDict()
CONFIG.declare(
    "variable_bounds_absolute_tolerance",
    ConfigValue(
        default=1e-4,
        domain=float,
        description="Absolute tolerance for considering a variable to be close "
        "to its bounds.",
    ),
)
CONFIG.declare(
    "variable_bounds_relative_tolerance",
    ConfigValue(
        default=1e-4,
        domain=float,
        description="Relative tolerance for considering a variable to be close "
        "to its bounds.",
    ),
)
CONFIG.declare(
    "variable_bounds_violation_tolerance",
    ConfigValue(
        default=0,
        domain=float,
        description="Absolute tolerance for considering a variable to violate its bounds.",
        doc="Absolute tolerance for considering a variable to violate its bounds. "
        "Some solvers relax bounds on variables thus allowing a small violation to be "
        "considered acceptable.",
    ),
)
CONFIG.declare(
    "constraint_residual_tolerance",
    ConfigValue(
        default=1e-5,
        domain=float,
        description="Absolute tolerance to use when checking constraint residuals.",
    ),
)
CONFIG.declare(
    "variable_large_value_tolerance",
    ConfigValue(
        default=1e4,
        domain=float,
        description="Absolute tolerance for considering a value to be large.",
    ),
)
CONFIG.declare(
    "variable_small_value_tolerance",
    ConfigValue(
        default=1e-4,
        domain=float,
        description="Absolute tolerance for considering a value to be small.",
    ),
)
CONFIG.declare(
    "variable_zero_value_tolerance",
    ConfigValue(
        default=1e-8,
        domain=float,
        description="Absolute tolerance for considering a value to be near to zero.",
    ),
)
CONFIG.declare(
    "jacobian_large_value_caution",
    ConfigValue(
        default=1e4,
        domain=float,
        description="Tolerance for raising a caution for large Jacobian values.",
    ),
)
CONFIG.declare(
    "jacobian_large_value_warning",
    ConfigValue(
        default=1e8,
        domain=float,
        description="Tolerance for raising a warning for large Jacobian values.",
    ),
)
CONFIG.declare(
    "jacobian_small_value_caution",
    ConfigValue(
        default=1e-4,
        domain=float,
        description="Tolerance for raising a caution for small Jacobian values.",
    ),
)
CONFIG.declare(
    "jacobian_small_value_warning",
    ConfigValue(
        default=1e-8,
        domain=float,
        description="Tolerance for raising a warning for small Jacobian values.",
    ),
)


@document_kwargs_from_configdict(CONFIG)
class DiagnosticsToolbox:
    """
    The IDAES Model DiagnosticsToolbox.

    To get started:

      1. Create an instance of your model (this does not need to be initialized yet).
      2. Fix variables until you have 0 degrees of freedom. Many of these tools presume
         a square model, and a square model should always be the foundation of any more
         advanced model.
      3. Create an instance of the DiagnosticsToolbox and provide the model to debug as
         the model argument.
      4. Call the ``report_structural_issues()`` method.

    Model diagnostics is an iterative process and you will likely need to run these
    tools multiple times to resolve all issues. After making a change to your model,
    you should always start from the beginning again to ensure the change did not
    introduce any new issues; i.e., always start from the report_structural_issues()
    method.

    Note that structural checks do not require the model to be initialized, thus users
    should start with these. Numerical checks require at least a partial solution to the
    model and should only be run once all structural issues have been resolved.

    Report methods will print a summary containing three parts:

    1. Warnings - these are critical issues that should be resolved before continuing.
       For each warning, a method will be suggested in the Next Steps section to get
       additional information.
    2. Cautions - these are things that could be correct but could also be the source of
       solver issues. Not all cautions need to be addressed, but users should investigate
       each one to ensure that the behavior is correct and that they will not be the source
       of difficulties later. Methods exist to provide more information on all cautions,
       but these will not appear in the Next Steps section.
    3. Next Steps - these are recommended methods to call from the DiagnosticsToolbox to
       get further information on warnings. If no warnings are found, this will suggest
       the next report method to call.

    Args:

        model: model to be diagnosed. The DiagnosticsToolbox does not support indexed Blocks.

    """

    def __init__(self, model: _BlockData, **kwargs):
        # TODO: In future may want to generalise this to accept indexed blocks
        # However, for now some of the tools do not support indexed blocks
        if not isinstance(model, _BlockData):
            raise TypeError(
                "model argument must be an instance of a Pyomo BlockData object "
                "(either a scalar Block or an element of an indexed Block)."
            )

        self._model = model
        self.config = CONFIG(kwargs)

    @property
    def model(self):
        """
        Model currently being diagnosed.
        """
        return self._model

    def display_external_variables(self, stream=stdout):
        """
        Prints a list of variables that appear within activated Constraints in the
        model but are not contained within the model themselves.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        ext_vars = []
        for v in variables_in_activated_constraints_set(self._model):
            if not _var_in_block(v, self._model):
                ext_vars.append(v.name)

        _write_report_section(
            stream=stream,
            lines_list=ext_vars,
            title="The following external variable(s) appear in constraints within the model:",
            header="=",
            footer="=",
        )

    def display_unused_variables(self, stream=stdout):
        """
        Prints a list of variables that do not appear in any activated Constraints.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        _write_report_section(
            stream=stream,
            lines_list=variables_not_in_activated_constraints_set(self._model),
            title="The following variable(s) do not appear in any activated constraints within the model:",
            header="=",
            footer="=",
        )

    def display_variables_fixed_to_zero(self, stream=stdout):
        """
        Prints a list of variables that are fixed to an absolute value of 0.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        _write_report_section(
            stream=stream,
            lines_list=_vars_fixed_to_zero(self._model),
            title="The following variable(s) are fixed to zero:",
            header="=",
            footer="=",
        )

    def display_variables_at_or_outside_bounds(self, stream=stdout):
        """
        Prints a list of variables with values that fall at or outside the bounds
        on the variable.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        _write_report_section(
            stream=stream,
            lines_list=[
                f"{v.name} ({'fixed' if v.fixed else 'free'}): value={value(v)} bounds={v.bounds}"
                for v in _vars_violating_bounds(
                    self._model,
                    tolerance=self.config.variable_bounds_violation_tolerance,
                )
            ],
            title="The following variable(s) have values at or outside their bounds:",
            header="=",
            footer="=",
        )

    def display_variables_with_none_value(self, stream=stdout):
        """
        Prints a list of variables with a value of None.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        _write_report_section(
            stream=stream,
            lines_list=_vars_with_none_value(self._model),
            title="The following variable(s) have a value of None:",
            header="=",
            footer="=",
        )

    def display_variables_with_value_near_zero(self, stream=stdout):
        """
        Prints a list of variables with a value close to zero. The tolerance
        for determining what is close to zero can be set in the class configuration
        options.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        _write_report_section(
            stream=stream,
            lines_list=[
                f"{v.name}: value={value(v)}"
                for v in _vars_near_zero(
                    self._model, self.config.variable_zero_value_tolerance
                )
            ],
            title="The following variable(s) have a value close to zero:",
            header="=",
            footer="=",
        )

    def display_variables_with_extreme_values(self, stream=stdout):
        """
        Prints a list of variables with extreme values.

        Tolerances can be set in the class configuration options.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        _write_report_section(
            stream=stream,
            lines_list=[
                f"{i.name}: {value(i)}"
                for i in _vars_with_extreme_values(
                    model=self._model,
                    large=self.config.variable_large_value_tolerance,
                    small=self.config.variable_small_value_tolerance,
                    zero=self.config.variable_zero_value_tolerance,
                )
            ],
            title="The following variable(s) have extreme values:",
            header="=",
            footer="=",
        )

    def display_variables_near_bounds(self, stream=stdout):
        """
        Prints a list of variables with values close to their bounds. Tolerance can
        be set in the class configuration options.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        _write_report_section(
            stream=stream,
            lines_list=[
                f"{v.name}: value={value(v)} bounds={v.bounds}"
                for v in variables_near_bounds_set(
                    self._model,
                    abs_tol=self.config.variable_bounds_absolute_tolerance,
                    rel_tol=self.config.variable_bounds_relative_tolerance,
                )
            ],
            title="The following variable(s) have values close to their bounds:",
            header="=",
            footer="=",
        )

    def display_components_with_inconsistent_units(self, stream=stdout):
        """
        Prints a list of all Constraints, Expressions and Objectives in the
        model with inconsistent units of measurement.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        _write_report_section(
            stream=stream,
            lines_list=identify_inconsistent_units(self._model),
            title="The following component(s) have unit consistency issues:",
            end_line="For more details on unit inconsistencies, import the "
            "assert_units_consistent method\nfrom pyomo.util.check_units",
            header="=",
            footer="=",
        )

    def display_constraints_with_large_residuals(self, stream=stdout):
        """
        Prints a list of Constraints with residuals greater than a specified tolerance.
        Tolerance can be set in the class configuration options.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        _write_report_section(
            stream=stream,
            lines_list=large_residuals_set(
                self._model, tol=self.config.constraint_residual_tolerance
            ),
            title="The following constraint(s) have large residuals:",
            header="=",
            footer="=",
        )

    def get_dulmage_mendelsohn_partition(self):
        """
        Performs a Dulmage-Mendelsohn partitioning on the model and returns
        the over- and under-constrained sub-problems.

        Returns:
            list-of-lists variables in each independent block of the under-constrained set
            list-of-lists constraints in each independent block of the under-constrained set
            list-of-lists variables in each independent block of the over-constrained set
            list-of-lists constraints in each independent block of the over-constrained set

        """
        igraph = IncidenceGraphInterface(self._model, include_inequality=False)
        var_dm_partition, con_dm_partition = igraph.dulmage_mendelsohn()

        # Collect under- and over-constrained sub-system
        uc_var = var_dm_partition.unmatched + var_dm_partition.underconstrained
        uc_con = con_dm_partition.underconstrained
        oc_var = var_dm_partition.overconstrained
        oc_con = con_dm_partition.overconstrained + con_dm_partition.unmatched

        uc_vblocks, uc_cblocks = igraph.get_connected_components(uc_var, uc_con)
        oc_vblocks, oc_cblocks = igraph.get_connected_components(oc_var, oc_con)

        return uc_vblocks, uc_cblocks, oc_vblocks, oc_cblocks

    def display_underconstrained_set(self, stream=stdout):
        """
        Prints the variables and constraints in the under-constrained sub-problem
        from a Dulmage-Mendelsohn partitioning.

        This can be used to identify the under-defined part of a model and thus
        where additional information (fixed variables or constraints) are required.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        uc_vblocks, uc_cblocks, _, _ = self.get_dulmage_mendelsohn_partition()

        stream.write("=" * MAX_STR_LENGTH + "\n")
        stream.write("Dulmage-Mendelsohn Under-Constrained Set\n\n")

        for i, uc_vblock in enumerate(uc_vblocks):
            stream.write(f"{TAB}Independent Block {i}:\n\n")
            stream.write(f"{2*TAB}Variables:\n\n")
            for v in uc_vblock:
                stream.write(f"{3*TAB}{v.name}\n")

            stream.write(f"\n{2*TAB}Constraints:\n\n")
            for c in uc_cblocks[i]:
                stream.write(f"{3*TAB}{c.name}\n")
            stream.write("\n")

        stream.write("=" * MAX_STR_LENGTH + "\n")

    def display_overconstrained_set(self, stream=stdout):
        """
        Prints the variables and constraints in the over-constrained sub-problem
        from a Dulmage-Mendelsohn partitioning.

        This can be used to identify the over-defined part of a model and thus
        where constraints must be removed or variables unfixed.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        _, _, oc_vblocks, oc_cblocks = self.get_dulmage_mendelsohn_partition()

        stream.write("=" * MAX_STR_LENGTH + "\n")
        stream.write("Dulmage-Mendelsohn Over-Constrained Set\n\n")

        for i, oc_vblock in enumerate(oc_vblocks):
            stream.write(f"{TAB}Independent Block {i}:\n\n")
            stream.write(f"{2*TAB}Variables:\n\n")
            for v in oc_vblock:
                stream.write(f"{3*TAB}{v.name}\n")

            stream.write(f"\n{2*TAB}Constraints:\n\n")
            for c in oc_cblocks[i]:
                stream.write(f"{3*TAB}{c.name}\n")
            stream.write("\n")

        stream.write("=" * MAX_STR_LENGTH + "\n")

    def display_variables_with_extreme_jacobians(self, stream=stdout):
        """
        Prints the variables associated with columns in the Jacobian with extreme
        L2 norms. This often indicates poorly scaled variables.

        Tolerances can be set via the DiagnosticsToolbox config.

        Args:
            stream: an I/O object to write the output to (default = stdout)

        Returns:
            None

        """
        xjc = extreme_jacobian_columns(
            m=self._model,
            scaled=False,
            large=self.config.jacobian_large_value_caution,
            small=self.config.jacobian_small_value_caution,
        )
        xjc.sort(key=lambda i: abs(log(i[0])), reverse=True)

        _write_report_section(
            stream=stream,
            lines_list=[f"{i[1].name}: {i[0]:.3E}" for i in xjc],
            title="The following variable(s) are associated with extreme Jacobian values:",
            header="=",
            footer="=",
        )

    def display_constraints_with_extreme_jacobians(self, stream=stdout):
        """
        Prints the constraints associated with rows in the Jacobian with extreme
        L2 norms. This often indicates poorly scaled constraints.

        Tolerances can be set via the DiagnosticsToolbox config.

        Args:
            stream: an I/O object to write the output to (default = stdout)

        Returns:
            None

        """
        xjr = extreme_jacobian_rows(
            m=self._model,
            scaled=False,
            large=self.config.jacobian_large_value_caution,
            small=self.config.jacobian_small_value_caution,
        )
        xjr.sort(key=lambda i: abs(log(i[0])), reverse=True)

        _write_report_section(
            stream=stream,
            lines_list=[f"{i[1].name}: {i[0]:.3E}" for i in xjr],
            title="The following constraint(s) are associated with extreme Jacobian values:",
            header="=",
            footer="=",
        )

    def display_extreme_jacobian_entries(self, stream=stdout):
        """
        Prints variables and constraints associated with entries in the Jacobian with extreme
        values. This can be indicative of poor scaling, especially for isolated terms (e.g.
        variables which appear only in one term of a single constraint).

        Tolerances can be set via the DiagnosticsToolbox config.

        Args:
            stream: an I/O object to write the output to (default = stdout)

        Returns:
            None

        """
        xje = extreme_jacobian_entries(
            m=self._model,
            scaled=False,
            large=self.config.jacobian_large_value_caution,
            small=self.config.jacobian_small_value_caution,
            zero=0,
        )
        xje.sort(key=lambda i: abs(log(i[0])), reverse=True)

        _write_report_section(
            stream=stream,
            lines_list=[f"{i[1].name}, {i[2].name}: {i[0]:.3E}" for i in xje],
            title="The following constraint(s) and variable(s) are associated with extreme Jacobian\nvalues:",
            header="=",
            footer="=",
        )

    # TODO: Block triangularization analysis
    # Number and size of blocks, polynomial degree of 1x1 blocks, simple pivot check of moderate sized sub-blocks?

    def _collect_structural_warnings(self):
        """
        Runs checks for structural warnings and returns two lists.

        Returns:
            warnings - list of warning messages from structural analysis
            next_steps - list of suggested next steps to further investigate warnings

        """
        uc = identify_inconsistent_units(self._model)
        uc_var, uc_con, oc_var, oc_con = self.get_dulmage_mendelsohn_partition()

        # Collect warnings
        warnings = []
        next_steps = []
        dof = degrees_of_freedom(self._model)
        if dof != 0:
            dstring = "Degrees"
            if abs(dof) == 1:
                dstring = "Degree"
            warnings.append(f"WARNING: {dof} {dstring} of Freedom")
        if len(uc) > 0:
            cstring = "Components"
            if len(uc) == 1:
                cstring = "Component"
            warnings.append(f"WARNING: {len(uc)} {cstring} with inconsistent units")
            next_steps.append(
                self.display_components_with_inconsistent_units.__name__ + "()"
            )
        if any(len(x) > 0 for x in [uc_var, uc_con, oc_var, oc_con]):
            warnings.append(
                f"WARNING: Structural singularity found\n"
                f"{TAB*2}Under-Constrained Set: {len(sum(uc_var, []))} "
                f"variables, {len(sum(uc_con, []))} constraints\n"
                f"{TAB*2}Over-Constrained Set: {len(sum(oc_var, []))} "
                f"variables, {len(sum(oc_con, []))} constraints"
            )

        if any(len(x) > 0 for x in [uc_var, uc_con]):
            next_steps.append(self.display_underconstrained_set.__name__ + "()")
        if any(len(x) > 0 for x in [oc_var, oc_con]):
            next_steps.append(self.display_overconstrained_set.__name__ + "()")

        return warnings, next_steps

    def _collect_structural_cautions(self):
        """
        Runs checks for structural cautions and returns a list.

        Returns:
            cautions - list of caution messages from structural analysis

        """
        # Collect cautions
        cautions = []
        zero_vars = _vars_fixed_to_zero(self._model)
        if len(zero_vars) > 0:
            vstring = "variables"
            if len(zero_vars) == 1:
                vstring = "variable"
            cautions.append(f"Caution: {len(zero_vars)} {vstring} fixed to 0")
        unused_vars = variables_not_in_activated_constraints_set(self._model)
        unused_vars_fixed = 0
        for v in unused_vars:
            if v.fixed:
                unused_vars_fixed += 1
        if len(unused_vars) > 0:
            vstring = "variables"
            if len(unused_vars) == 1:
                vstring = "variable"
            cautions.append(
                f"Caution: {len(unused_vars)} "
                f"unused {vstring} ({unused_vars_fixed} fixed)"
            )

        return cautions

    def _collect_numerical_warnings(self, jac=None, nlp=None):
        """
        Runs checks for numerical warnings and returns two lists.

        Returns:
            warnings - list of warning messages from numerical analysis
            next_steps - list of suggested next steps to further investigate warnings

        """
        if jac is None or nlp is None:
            jac, nlp = get_jacobian(self._model, scaled=False)

        warnings = []
        next_steps = []

        # Large residuals
        large_residuals = large_residuals_set(
            self._model, tol=self.config.constraint_residual_tolerance
        )
        if len(large_residuals) > 0:
            cstring = "Constraints"
            if len(large_residuals) == 1:
                cstring = "Constraint"
            warnings.append(
                f"WARNING: {len(large_residuals)} {cstring} with large residuals"
            )
            next_steps.append(
                self.display_constraints_with_large_residuals.__name__ + "()"
            )

        # Variables outside bounds
        violated_bounds = _vars_violating_bounds(
            self._model, tolerance=self.config.variable_bounds_violation_tolerance
        )
        if len(violated_bounds) > 0:
            cstring = "Variables"
            if len(violated_bounds) == 1:
                cstring = "Variable"
            warnings.append(
                f"WARNING: {len(violated_bounds)} {cstring} at or outside bounds"
            )
            next_steps.append(
                self.display_variables_at_or_outside_bounds.__name__ + "()"
            )

        # Extreme Jacobian rows and columns
        jac_col = extreme_jacobian_columns(
            jac=jac,
            nlp=nlp,
            large=self.config.jacobian_large_value_warning,
            small=self.config.jacobian_small_value_warning,
        )
        if len(jac_col) > 0:
            cstring = "Variables"
            if len(jac_col) == 1:
                cstring = "Variable"
            warnings.append(
                f"WARNING: {len(jac_col)} {cstring} with extreme Jacobian values"
            )
            next_steps.append(
                self.display_variables_with_extreme_jacobians.__name__ + "()"
            )

        jac_row = extreme_jacobian_rows(
            jac=jac,
            nlp=nlp,
            large=self.config.jacobian_large_value_warning,
            small=self.config.jacobian_small_value_warning,
        )
        if len(jac_row) > 0:
            cstring = "Constraints"
            if len(jac_row) == 1:
                cstring = "Constraint"
            warnings.append(
                f"WARNING: {len(jac_row)} {cstring} with extreme Jacobian values"
            )
            next_steps.append(
                self.display_constraints_with_extreme_jacobians.__name__ + "()"
            )

        return warnings, next_steps

    def _collect_numerical_cautions(self, jac=None, nlp=None):
        """
        Runs checks for numerical cautions and returns a list.

        Returns:
            cautions - list of caution messages from numerical analysis

        """
        if jac is None or nlp is None:
            jac, nlp = get_jacobian(self._model, scaled=False)

        cautions = []

        # Variables near bounds
        near_bounds = variables_near_bounds_set(
            self._model,
            abs_tol=self.config.variable_bounds_absolute_tolerance,
            rel_tol=self.config.variable_bounds_relative_tolerance,
        )
        if len(near_bounds) > 0:
            cstring = "Variables"
            if len(near_bounds) == 1:
                cstring = "Variable"
            cautions.append(
                f"Caution: {len(near_bounds)} {cstring} with value close to their bounds"
            )

        # Variables near zero
        near_zero = _vars_near_zero(
            self._model, self.config.variable_zero_value_tolerance
        )
        if len(near_zero) > 0:
            cstring = "Variables"
            if len(near_zero) == 1:
                cstring = "Variable"
            cautions.append(
                f"Caution: {len(near_zero)} {cstring} with value close to zero"
            )

        # Variables with extreme values
        xval = _vars_with_extreme_values(
            model=self._model,
            large=self.config.variable_large_value_tolerance,
            small=self.config.variable_small_value_tolerance,
            zero=self.config.variable_zero_value_tolerance,
        )
        if len(xval) > 0:
            cstring = "Variables"
            if len(xval) == 1:
                cstring = "Variable"
            cautions.append(f"Caution: {len(xval)} {cstring} with extreme value")

        # Variables with value None
        none_value = _vars_with_none_value(self._model)
        if len(none_value) > 0:
            cstring = "Variables"
            if len(none_value) == 1:
                cstring = "Variable"
            cautions.append(f"Caution: {len(none_value)} {cstring} with None value")

        # Extreme Jacobian rows and columns
        jac_col = extreme_jacobian_columns(
            jac=jac,
            nlp=nlp,
            large=self.config.jacobian_large_value_caution,
            small=self.config.jacobian_small_value_caution,
        )
        if len(jac_col) > 0:
            cstring = "Variables"
            if len(jac_col) == 1:
                cstring = "Variable"
            cautions.append(
                f"Caution: {len(jac_col)} {cstring} with extreme Jacobian values"
            )

        jac_row = extreme_jacobian_rows(
            jac=jac,
            nlp=nlp,
            large=self.config.jacobian_large_value_caution,
            small=self.config.jacobian_small_value_caution,
        )
        if len(jac_row) > 0:
            cstring = "Constraints"
            if len(jac_row) == 1:
                cstring = "Constraint"
            cautions.append(
                f"Caution: {len(jac_row)} {cstring} with extreme Jacobian values"
            )

        # Extreme Jacobian entries
        extreme_jac = extreme_jacobian_entries(
            jac=jac,
            nlp=nlp,
            large=self.config.jacobian_large_value_caution,
            small=self.config.jacobian_small_value_caution,
            zero=0,
        )
        if len(extreme_jac) > 0:
            cstring = "Entries"
            if len(extreme_jac) == 1:
                cstring = "Entry"
            cautions.append(f"Caution: {len(extreme_jac)} extreme Jacobian {cstring}")

        return cautions

    def assert_no_structural_warnings(self):
        """
        Checks for structural warnings in the model and raises an AssertionError
        if any are found.

        Raises:
            AssertionError if any warnings are identified by structural analysis.

        """
        warnings, _ = self._collect_structural_warnings()
        if len(warnings) > 0:
            raise AssertionError(f"Structural issues found ({len(warnings)}).")

    def assert_no_numerical_warnings(self):
        """
        Checks for numerical warnings in the model and raises an AssertionError
        if any are found.

        Raises:
            AssertionError if any warnings are identified by numerical analysis.

        """
        warnings, _ = self._collect_numerical_warnings()
        if len(warnings) > 0:
            raise AssertionError(f"Numerical issues found ({len(warnings)}).")

    def report_structural_issues(self, stream=stdout):
        """
        Generates a summary report of any structural issues identified in the model provided
        and suggests next steps for debugging the model.

        This should be the first method called when debugging a model and after any change
        is made to the model. These checks can be run before trying to initialize and solve
        the model.

        Args:
            stream: I/O object to write report to (default = stdout)

        Returns:
            None

        """
        # Potential evaluation errors
        # TODO: High Index?
        stats = _collect_model_statistics(self._model)
        warnings, next_steps = self._collect_structural_warnings()
        cautions = self._collect_structural_cautions()

        _write_report_section(
            stream=stream, lines_list=stats, title="Model Statistics", header="="
        )
        _write_report_section(
            stream=stream,
            lines_list=warnings,
            title=f"{len(warnings)} WARNINGS",
            line_if_empty="No warnings found!",
        )
        _write_report_section(
            stream=stream,
            lines_list=cautions,
            title=f"{len(cautions)} Cautions",
            line_if_empty="No cautions found!",
        )
        _write_report_section(
            stream=stream,
            lines_list=next_steps,
            title="Suggested next steps:",
            line_if_empty="Try to initialize/solve your model and then call report_numerical_issues()",
            footer="=",
        )

    def report_numerical_issues(self, stream=stdout):
        """
        Generates a summary report of any numerical issues identified in the model provided
        and suggest next steps for debugging model.

        Numerical checks should only be performed once all structural issues have been resolved,
        and require that at least a partial solution to the model is available.

        Args:
            stream: I/O object to write report to (default = stdout)

        Returns:
            None

        """
        jac, nlp = get_jacobian(self._model, scaled=False)

        warnings, next_steps = self._collect_numerical_warnings(jac=jac, nlp=nlp)
        cautions = self._collect_numerical_cautions(jac=jac, nlp=nlp)

        stats = []
        stats.append(
            f"Jacobian Condition Number: {jacobian_cond(jac=jac, scaled=False):.3E}"
        )
        _write_report_section(
            stream=stream, lines_list=stats, title="Model Statistics", header="="
        )
        _write_report_section(
            stream=stream,
            lines_list=warnings,
            title=f"{len(warnings)} WARNINGS",
            line_if_empty="No warnings found!",
        )
        _write_report_section(
            stream=stream,
            lines_list=cautions,
            title=f"{len(cautions)} Cautions",
            line_if_empty="No cautions found!",
        )
        _write_report_section(
            stream=stream,
            lines_list=next_steps,
            title="Suggested next steps:",
            line_if_empty=f"If you still have issues converging your model consider:\n"
            f"{TAB*2}svd_analysis(TBA)\n{TAB*2}degeneracy_hunter (TBA)",
            footer="=",
        )


class DegeneracyHunter:
    """
    Degeneracy Hunter is a collection of utility functions to assist in mathematical
    modeling in Pyomo.
    """

    def __init__(self, block_or_jac, solver=None):
        """Initialize Degeneracy Hunter Object

        Args:
            block_or_jac: Pyomo model or Jacobian
            solver: Pyomo SolverFactory

        Notes:
            Passing a Jacobian to Degeneracy Hunter is current untested.

        """
        msg = (
            "DegeneracyHunter is being deprecated in favor of the new "
            "DiagnosticsToolbox."
        )
        deprecation_warning(msg=msg, logger=_log, version="2.2.0", remove_in="3.0.0")

        block_like = False
        try:
            block_like = issubclass(block_or_jac.ctype, Block)
        except AttributeError:
            pass

        if block_like:
            # Add Pyomo model to the object
            self.block = block_or_jac

            # setup pynumero interface
            if not AmplInterface.available():
                raise RuntimeError("Pynumero not available.")
            self.nlp = PyomoNLP(self.block)

            # Get the scaled Jacobian of equality constraints
            self.jac_eq = get_jacobian(self.block, equality_constraints_only=True)[0]

            # Create a list of equality constraint names
            self.eq_con_list = self.nlp.get_pyomo_equality_constraints()
            self.var_list = self.nlp.get_pyomo_variables()

            self.candidate_eqns = None

        elif type(block_or_jac) is np.array:  # pylint: disable=unidiomatic-typecheck
            raise NotImplementedError(
                "Degeneracy Hunter currently only supports analyzing a Pyomo model"
            )

            # # TODO: Need to refactor, document, and test support for Jacobian
            # self.jac_eq = block_or_jac
            # self.eq_con_list = None

        else:
            raise TypeError("Check the type for 'block_or_jac'")

        # number of equality constraints, variables
        self.n_eq = self.jac_eq.shape[0]
        self.n_var = self.jac_eq.shape[1]

        # Initialize solver
        if solver is None:
            # TODO: Test performance with open solvers such as cbc
            self.solver = SolverFactory("gurobi")
            self.solver.options = {"NumericFocus": 3}

        else:
            # TODO: Make this a custom exception following IDAES standards
            # assert type(solver) is SolverFactory, "Argument solver should be type SolverFactory"
            self.solver = solver

        # Create spot to store singular values
        self.s = None

        # Set constants for MILPs
        self.max_nu = 1e5
        self.min_nonzero_nu = 1e-5

    def check_residuals(self, tol=1e-5, print_level=1, sort=True):
        """
        Method to return a ComponentSet of all Constraint components with a
        residual greater than a given threshold which appear in a model.

        Args:
            block: model to be studied
            tol: residual threshold for inclusion in ComponentSet
            print_level: controls to extend of output to the screen:

                * 0 - nothing printed
                * 1 - only name of constraint printed
                * 2 - each constraint is pretty printed
                * 3 - pretty print each constraint, then print value for included variable

            sort: sort residuals in descending order for printing

        Returns:
            A ComponentSet including all Constraint components with a residual
            greater than tol which appear in block

        """

        if print_level > 0:
            residual_values = large_residuals_set(self.block, tol, True)
        else:
            return large_residuals_set(self.block, tol, False)

        print(" ")
        if len(residual_values) > 0:
            print("All constraints with residuals larger than", tol, ":")
            if print_level == 1:
                print("Count\tName\t|residual|")

            if sort:
                residual_values = dict(
                    sorted(residual_values.items(), key=itemgetter(1), reverse=True)
                )

            for i, (c, r) in enumerate(residual_values.items()):
                if print_level == 1:
                    # Basic print statement. count, constraint, residual
                    print(i, "\t", c, "\t", r)
                else:
                    # Pretty print constraint
                    print("\ncount =", i, "\t|residual| =", r)
                    c.pprint()

                if print_level == 2:
                    # print values and bounds for each variable in the constraint
                    print("variable\tlower\tvalue\tupper")
                    for v in identify_variables(c.body):
                        self.print_variable_bounds(v)
        else:
            print("No constraints with residuals larger than", tol, "!")

        return residual_values.keys()

    def check_variable_bounds(
        self, tol=1e-5, relative=False, skip_lb=False, skip_ub=False, verbose=True
    ):
        """
        Return a ComponentSet of all variables within a tolerance of their bounds.

        Args:
            block: model to be studied
            tol: residual threshold for inclusion in ComponentSet (default = 1e-5)
            relative : Boolean, use relative tolerance (default = False)
            skip_lb: Boolean to skip lower bound (default = False)
            skip_ub: Boolean to skip upper bound (default = False)
            verbose: Boolean to toggle on printing to screen (default = True)

        Returns:
            A ComponentSet including all Constraint components with a residual
            greater than tol which appear in block

        """
        vnbs = variables_near_bounds_set(self.block, tol, relative, skip_lb, skip_ub)

        if verbose:
            print(" ")
            if relative:
                s = "(relative)"
            else:
                s = "(absolute)"
            if len(vnbs) > 0:
                print("Variables within", tol, s, "of their bounds:")
                print("variable\tlower\tvalue\tupper")
                for v in vnbs:
                    self.print_variable_bounds(v)
            else:
                print("No variables within", tol, s, "of their bounds.")

        return vnbs

    def check_rank_equality_constraints(self, tol=1e-6, dense=False):
        """
        Method to check the rank of the Jacobian of the equality constraints

        Args:
            tol: Tolerance for smallest singular value (default=1E-6)
            dense: If True, use a dense svd to perform singular value analysis,
                which tends to be slower but more reliable than svds

        Returns:
            Number of singular values less than tolerance (-1 means error)

        """

        print("\nChecking rank of Jacobian of equality constraints...")

        print(
            "Model contains",
            self.n_eq,
            "equality constraints and",
            self.n_var,
            "variables.",
        )

        counter = 0
        if self.n_eq > 1:
            if self.s is None:
                self.svd_analysis(dense=dense)

            n = len(self.s)

            print("Smallest singular value(s):")
            for i in range(n):
                print("%.3E" % self.s[i])
                if self.s[i] < tol:
                    counter += 1
        else:
            print(f"Only singular value: {norm(self.jac_eq,'fro')}")

        return counter

    # TODO: Refactor, this should not be a staticmethod
    @staticmethod
    def _prepare_ids_milp(jac_eq, M=1e5):
        """
        Prepare MILP to compute the irreducible degenerate set

        Args:
            jac_eq Jacobian of equality constraints [matrix]
            M: largest value for nu

        Returns:
            m_dh: Pyomo model to calculate irreducible degenerate sets

        """

        n_eq = jac_eq.shape[0]
        n_var = jac_eq.shape[1]

        # Create Pyomo model for irreducible degenerate set
        m_dh = ConcreteModel()

        # Create index for constraints
        m_dh.C = Set(initialize=range(n_eq))

        m_dh.V = Set(initialize=range(n_var))

        # Add variables
        m_dh.nu = Var(m_dh.C, bounds=(-M, M), initialize=1.0)
        m_dh.y = Var(m_dh.C, domain=Binary)

        # Constraint to enforce set is degenerate
        if issparse(jac_eq):
            m_dh.J = jac_eq.tocsc()

            def eq_degenerate(m_dh, v):
                # Find the columns with non-zero entries
                C_ = find(m_dh.J[:, v])[0]
                return sum(m_dh.J[c, v] * m_dh.nu[c] for c in C_) == 0

        else:
            m_dh.J = jac_eq

            def eq_degenerate(m_dh, v):
                return sum(m_dh.J[c, v] * m_dh.nu[c] for c in m_dh.C) == 0

        m_dh.degenerate = Constraint(m_dh.V, rule=eq_degenerate)

        def eq_lower(m_dh, c):
            return -M * m_dh.y[c] <= m_dh.nu[c]

        m_dh.lower = Constraint(m_dh.C, rule=eq_lower)

        def eq_upper(m_dh, c):
            return m_dh.nu[c] <= M * m_dh.y[c]

        m_dh.upper = Constraint(m_dh.C, rule=eq_upper)

        m_dh.obj = Objective(expr=sum(m_dh.y[c] for c in m_dh.C))

        return m_dh

    # TODO: Refactor, this should not be a staticmethod
    @staticmethod
    def _prepare_find_candidates_milp(jac_eq, M=1e5, m_small=1e-5):
        """
        Prepare MILP to find candidate equations for consider for IDS

        Args:
            jac_eq Jacobian of equality constraints [matrix]
            M: maximum value for nu
            m_small: smallest value for nu to be considered non-zero

        Returns:
            m_fc: Pyomo model to find candidates

        """

        n_eq = jac_eq.shape[0]
        n_var = jac_eq.shape[1]

        # Create Pyomo model for irreducible degenerate set
        m_dh = ConcreteModel()

        # Create index for constraints
        m_dh.C = Set(initialize=range(n_eq))

        m_dh.V = Set(initialize=range(n_var))

        # Specify minimum size for nu to be considered non-zero
        m_dh.m_small = m_small

        # Add variables
        m_dh.nu = Var(m_dh.C, bounds=(-M - m_small, M + m_small), initialize=1.0)
        m_dh.y_pos = Var(m_dh.C, domain=Binary)
        m_dh.y_neg = Var(m_dh.C, domain=Binary)
        m_dh.abs_nu = Var(m_dh.C, bounds=(0, M + m_small))

        m_dh.pos_xor_neg = Constraint(m_dh.C)

        # Constraint to enforce set is degenerate
        if issparse(jac_eq):
            m_dh.J = jac_eq.tocsc()

            def eq_degenerate(m_dh, v):
                if np.sum(np.abs(m_dh.J[:, v])) > 1e-6:
                    # Find the columns with non-zero entries
                    C_ = find(m_dh.J[:, v])[0]
                    return sum(m_dh.J[c, v] * m_dh.nu[c] for c in C_) == 0
                else:
                    # This variable does not appear in any constraint
                    return Constraint.Skip

        else:
            m_dh.J = jac_eq

            def eq_degenerate(m_dh, v):
                if np.sum(np.abs(m_dh.J[:, v])) > 1e-6:
                    return sum(m_dh.J[c, v] * m_dh.nu[c] for c in m_dh.C) == 0
                else:
                    # This variable does not appear in any constraint
                    return Constraint.Skip

        m_dh.pprint()

        m_dh.degenerate = Constraint(m_dh.V, rule=eq_degenerate)

        # When y_pos = 1, nu >= m_small
        # When y_pos = 0, nu >= - m_small
        def eq_pos_lower(m_dh, c):
            return m_dh.nu[c] >= -m_small + 2 * m_small * m_dh.y_pos[c]

        m_dh.pos_lower = Constraint(m_dh.C, rule=eq_pos_lower)

        # When y_pos = 1, nu <= M + m_small
        # When y_pos = 0, nu <= m_small
        def eq_pos_upper(m_dh, c):
            return m_dh.nu[c] <= M * m_dh.y_pos[c] + m_small

        m_dh.pos_upper = Constraint(m_dh.C, rule=eq_pos_upper)

        # When y_neg = 1, nu <= -m_small
        # When y_neg = 0, nu <= m_small
        def eq_neg_upper(m_dh, c):
            return m_dh.nu[c] <= m_small - 2 * m_small * m_dh.y_neg[c]

        m_dh.neg_upper = Constraint(m_dh.C, rule=eq_neg_upper)

        # When y_neg = 1, nu >= -M - m_small
        # When y_neg = 0, nu >= - m_small
        def eq_neg_lower(m_dh, c):
            return m_dh.nu[c] >= -M * m_dh.y_neg[c] - m_small

        m_dh.neg_lower = Constraint(m_dh.C, rule=eq_neg_lower)

        # Absolute value
        def eq_abs_lower(m_dh, c):
            return -m_dh.abs_nu[c] <= m_dh.nu[c]

        m_dh.abs_lower = Constraint(m_dh.C, rule=eq_abs_lower)

        def eq_abs_upper(m_dh, c):
            return m_dh.nu[c] <= m_dh.abs_nu[c]

        m_dh.abs_upper = Constraint(m_dh.C, rule=eq_abs_upper)

        # At least one constraint must be in the degenerate set
        m_dh.degenerate_set_nonempty = Constraint(
            expr=sum(m_dh.y_pos[c] + m_dh.y_neg[c] for c in m_dh.C) >= 1
        )

        # Minimize the L1-norm of nu
        m_dh.obj = Objective(expr=sum(m_dh.abs_nu[c] for c in m_dh.C))

        return m_dh

    # TODO: Refactor, this should not be a staticmethod
    @staticmethod
    def _check_candidate_ids(ids_milp, solver, c, eq_con_list=None, tee=False):
        """Solve MILP to check if equation 'c' is a significant component in an irreducible
        degenerate set

        Args:
            ids_milp: Pyomo model to calculate IDS
            solver: Pyomo solver (must support MILP)
            c: index for the constraint to consider [integer]
            eq_con_list: names of equality constraints. If none, use elements of ids_milp (default=None)
            tee: Boolean, print solver output (default = False)

        Returns:
            ids: either None or dictionary containing the IDS

        """

        # Fix weight on candidate equation
        ids_milp.nu[c].fix(1.0)

        # Solve MILP
        results = solver.solve(ids_milp, tee=tee)

        ids_milp.nu[c].unfix()

        if check_optimal_termination(results):
            # We found an irreducible degenerate set

            # Create empty dictionary
            ids_ = {}

            for i in ids_milp.C:
                # Check if constraint is included
                if ids_milp.y[i]() > 0.9:
                    # If it is, save the value of nu
                    if eq_con_list is None:
                        name = i
                    else:
                        name = eq_con_list[i]
                    ids_[name] = ids_milp.nu[i]()
            return ids_
        else:
            return None

    # TODO: Refactor, this should not be a staticmethod
    @staticmethod
    def _find_candidate_eqs(candidates_milp, solver, eq_con_list=None, tee=False):
        """Solve MILP to generate set of candidate equations

        Arguments:
            candidates_milp: Pyomo model to calculate IDS
            solver: Pyomo solver (must support MILP)
            eq_con_list: names of equality constraints. If none, use elements of ids_milp (default=None)
            tee: Boolean, print solver output (default = False)

        Returns:
            candidate_eqns: either None or list of indices
            degenerate_set: either None or dictionary containing the degenerate_set

        """

        results = solver.solve(candidates_milp, tee=tee)

        if check_optimal_termination(results):
            # We found a degenerate set

            # Create empty dictionary
            ds_ = {}

            # Create empty list
            candidate_eqns = []

            for i in candidates_milp.C:
                # Check if constraint is included
                if candidates_milp.abs_nu[i]() > candidates_milp.m_small * 0.99:
                    # If it is, save the value of nu
                    if eq_con_list is None:
                        name = i
                    else:
                        name = eq_con_list[i]
                    ds_[name] = candidates_milp.nu[i]()
                    candidate_eqns.append(i)

            return candidate_eqns, ds_
        else:
            return None, None

    def svd_analysis(self, n_sv=None, dense=False):
        """
        Perform SVD analysis of the constraint Jacobian

        Args:
            n_sv: number of smallest singular values to compute
            dense: If True, use a dense svd to perform singular value analysis,
                which tends to be slower but more reliable than svds

        Returns:
            None

        Actions:
            Stores SVD results in object

        """
        if n_sv is None:
            # Determine the number of singular values to compute
            # The "-1" is needed to avoid an error with svds
            n_sv = min(10, min(self.n_eq, self.n_var) - 1)

        if self.n_eq > 1:
            if n_sv >= min(self.n_eq, self.n_var):
                raise ValueError(
                    f"For a {self.n_eq} by {self.n_var} system, svd_analysis "
                    f"can compute at most {min(self.n_eq, self.n_var) - 1} "
                    f"singular values and vectors, but {n_sv} were called for."
                )
            if n_sv < 1:
                raise ValueError(f"Nonsense value for n_sv={n_sv} received.")
            print("Computing the", n_sv, "smallest singular value(s)")

            # Perform SVD
            # Recall J is a n_eq x n_var matrix
            # Thus U is a n_eq x n_eq matrix
            # And V is a n_var x n_var
            # (U or V may be smaller in economy mode)
            if dense:
                u, s, vT = svd(self.jac_eq.todense(), full_matrices=False)
                # Reorder singular values and vectors so that the singular
                # values are from least to greatest
                u = np.flip(u[:, -n_sv:], axis=1)
                s = np.flip(s[-n_sv:], axis=0)
                vT = np.flip(vT[-n_sv:, :], axis=0)
            else:
                # svds does not guarantee the order in which it generates
                # singular values, but typically generates them least-to-greatest.
                # Maybe the warning is for singular values of nearly equal
                # magnitude or a similar edge case?
                u, s, vT = svds(self.jac_eq, k=n_sv, which="SM")  # , solver='lobpcg')

            # Save results
            self.u = u
            self.s = s
            self.v = vT.transpose()

        else:
            raise ValueError(
                "Model needs at least 2 equality constraints to perform svd_analysis."
            )

    def underdetermined_variables_and_constraints(self, n_calc=1, tol=0.1, dense=False):
        """
        Determines constraints and variables associated with the smallest
        singular values by having large components in the left and right
        singular vectors, respectively, associated with those singular values.

        Args:
            n_calc: The singular value, as ordered from least to greatest
                starting from 1, to calculate associated constraints and variables
            tol: Size below which to ignore constraints and variables in
                the singular vector
            dense: If True, use a dense svd to perform singular value analysis,
                which tends to be slower but more reliable than svds

        Returns:
            None

        """
        if self.s is None:
            self.svd_analysis(
                n_sv=max(n_calc, min(10, min(self.n_eq, self.n_var) - 1)), dense=dense
            )
        n_sv = len(self.s)
        if n_sv < n_calc:
            raise ValueError(
                f"User wanted constraints and variables associated "
                f"with the {n_calc}-th smallest singular value, "
                f"but only {n_sv} small singular values have been "
                f"calculated. Run svd_analysis again and specify "
                f"n_sv>={n_calc}."
            )
        print("Column:    Variable")
        for i in np.where(abs(self.v[:, n_calc - 1]) > tol)[0]:
            print(str(i) + ": " + self.var_list[i].name)
        print("")
        print("Row:    Constraint")
        for i in np.where(abs(self.u[:, n_calc - 1]) > tol)[0]:
            print(str(i) + ": " + self.eq_con_list[i].name)

    def find_candidate_equations(self, verbose=True, tee=False):
        """
        Solve MILP to find a degenerate set and candidate equations

        Args:
            verbose: Print information to the screen (default=True)
            tee: Print solver output to screen (default=True)

        Returns:
            ds: either None or dictionary of candidate equations

        """

        if verbose:
            print("*** Searching for a Single Degenerate Set ***")
            print("Building MILP model...")
        self.candidates_milp = self._prepare_find_candidates_milp(
            self.jac_eq, self.max_nu, self.min_nonzero_nu
        )

        if verbose:
            print("Solving MILP model...")
        ce, ds = self._find_candidate_eqs(
            self.candidates_milp, self.solver, self.eq_con_list, tee
        )

        if ce is not None:
            self.candidate_eqns = ce

        return ds

    def find_irreducible_degenerate_sets(self, verbose=True, tee=False):
        """
        Compute irreducible degenerate sets

        Args:
            verbose: Print information to the screen (default=True)
            tee: Print solver output to screen (default=True)

        Returns:
            irreducible_degenerate_sets: list of irreducible degenerate sets

        """

        # If there are no candidate equations, find them!
        if not self.candidate_eqns:
            self.find_candidate_equations()

        irreducible_degenerate_sets = []

        # Check if it is empty or None
        if self.candidate_eqns:
            if verbose:
                print("*** Searching for Irreducible Degenerate Sets ***")
                print("Building MILP model...")
            self.dh_milp = self._prepare_ids_milp(self.jac_eq, self.max_nu)

            # Loop over candidate equations
            for i, c in enumerate(self.candidate_eqns):
                if verbose:
                    print("Solving MILP", i + 1, "of", len(self.candidate_eqns), "...")

                # Check if equation 'c' is a major element of an IDS
                ids_ = self._check_candidate_ids(
                    self.dh_milp, self.solver, c, self.eq_con_list, tee
                )

                if ids_ is not None:
                    irreducible_degenerate_sets.append(ids_)

            if verbose:
                for i, s in enumerate(irreducible_degenerate_sets):
                    print("\nIrreducible Degenerate Set", i)
                    print("nu\tConstraint Name")
                    for k, v in s.items():
                        print(v, "\t", k)
        else:
            print("No candidate equations. The Jacobian is likely full rank.")

        return irreducible_degenerate_sets

    ### Helper Functions

    # Note: This makes sense as a static method
    @staticmethod
    def print_variable_bounds(v):
        """
        Print variable, bounds, and value

        Args:
            v: variable

        Returns:
            None

        """
        print(v, "\t\t", v.lb, "\t", v.value, "\t", v.ub)


def get_valid_range_of_component(component):
    """
    Return the valid range for a component as specified in the model metadata.

    Args:
        component: Pyomo component to get valid range for

    Returns:
        valid range for component if found. This will either be a 2-tuple (low, high) or None.

    Raises:
        AttributeError if metadata object not found

    """
    # Get metadata for component
    parent = component.parent_block()

    try:
        if hasattr(parent, "params"):
            meta = parent.params.get_metadata().properties
        else:
            meta = parent.get_metadata().properties
    except AttributeError:
        raise AttributeError(f"Could not find metadata for component {component.name}")

    # Get valid range from metadata
    try:
        n, i = meta.get_name_and_index(component.parent_component().local_name)
        cmeta = getattr(meta, n)[i]
        valid_range = cmeta.valid_range
    except ValueError:
        # Assume no metadata for this property
        _log.debug(f"No metadata entry for component {component.name}; returning None")
        valid_range = None

    return valid_range


def set_bounds_from_valid_range(component, descend_into=True):
    """
    Set bounds on Pyomo components based on valid range recorded in model metadata.
    WARNING - this function will overwrite any bounds already set on the component/model.

    This function will iterate over component data objects in Blocks and indexed components.

    Args:
        component: Pyomo component to set bounds on. This can be a Block, Var or Param.
        descend_into: (optional) Whether to descend into components on child Blocks (default=True)

    Returns:
         None

    """
    if component.is_indexed():
        for k in component:
            set_bounds_from_valid_range(component[k])
    elif isinstance(component, _BlockData):
        for i in component.component_data_objects(
            ctype=[Var, Param], descend_into=descend_into
        ):
            set_bounds_from_valid_range(i)
    elif not hasattr(component, "bounds"):
        raise TypeError(
            f"Component {component.name} does not have bounds. Only Vars and Params have bounds."
        )
    else:
        valid_range = get_valid_range_of_component(component)

        if valid_range is None:
            valid_range = (None, None)

        component.setlb(valid_range[0])
        component.setub(valid_range[1])


def list_components_with_values_outside_valid_range(component, descend_into=True):
    """
    Return a list of component objects with values outside the valid range specified in the model
    metadata.

    This function will iterate over component data objects in Blocks and indexed components.

    Args:
        component: Pyomo component to search for component outside of range on.
            This can be a Block, Var or Param.
        descend_into: (optional) Whether to descend into components on child Blocks (default=True)

    Returns:
         list of component objects found with values outside the valid range.
    """
    comp_list = []

    if component.is_indexed():
        for k in component:
            comp_list.extend(
                list_components_with_values_outside_valid_range(component[k])
            )
    elif isinstance(component, _BlockData):
        for i in component.component_data_objects(
            ctype=[Var, Param], descend_into=descend_into
        ):
            comp_list.extend(list_components_with_values_outside_valid_range(i))
    else:
        valid_range = get_valid_range_of_component(component)

        if valid_range is not None:
            cval = value(component)
            if cval is not None and (cval < valid_range[0] or cval > valid_range[1]):
                comp_list.append(component)

    return comp_list


def ipopt_solve_halt_on_error(model, options=None):
    """
    Run IPOPT to solve model with debugging output enabled.

    This function calls IPOPT to solve the model provided with settings
    to halt on AMPL evaluation errors and report these with symbolic names.

    Args:
        model: Pyomo model to be solved.
        options: solver options to be passed to IPOPT

    Returns:
        Pyomo solver results dict

    """
    if options is None:
        options = {}

    solver = SolverFactory("ipopt")
    solver.options = options
    solver.options["halt_on_ampl_error"] = "yes"

    return solver.solve(
        model, tee=True, symbolic_solver_labels=True, export_defined_variables=False
    )


# -------------------------------------------------------------------------------------------
# Private module functions
def _var_in_block(var, block):
    parent = var.parent_block()
    while parent is not None:
        if parent is block:
            return True
        parent = parent.parent_block()
    return False


def _vars_fixed_to_zero(model):
    # Set of variables fixed to 0
    zero_vars = ComponentSet()
    for v in model.component_data_objects(Var, descend_into=True):
        if v.fixed and value(v) == 0:
            zero_vars.add(v)
    return zero_vars


def _vars_near_zero(model, variable_zero_value_tolerance):
    # Set of variables with values close to 0
    near_zero_vars = ComponentSet()
    for v in model.component_data_objects(Var, descend_into=True):
        if v.value is not None and abs(value(v)) <= variable_zero_value_tolerance:
            near_zero_vars.add(v)
    return near_zero_vars


def _vars_violating_bounds(model, tolerance):
    violated_bounds = ComponentSet()
    for v in model.component_data_objects(Var, descend_into=True):
        if v.value is not None:
            if v.lb is not None and v.value <= v.lb - tolerance:
                violated_bounds.add(v)
            elif v.ub is not None and v.value >= v.ub + tolerance:
                violated_bounds.add(v)

    return violated_bounds


def _vars_with_none_value(model):
    none_value = ComponentSet()
    for v in model.component_data_objects(Var, descend_into=True):
        if v.value is None:
            none_value.add(v)

    return none_value


def _vars_with_extreme_values(model, large, small, zero):
    extreme_vars = ComponentSet()
    for v in model.component_data_objects(Var, descend_into=True):
        if v.value is not None:
            mag = abs(value(v))
            if mag > abs(large):
                extreme_vars.add(v)
            elif mag < abs(small) and mag > abs(zero):
                extreme_vars.add(v)

    return extreme_vars


def _write_report_section(
    stream,
    lines_list,
    title=None,
    line_if_empty=None,
    end_line=None,
    header="-",
    footer=None,
):
    """
    Writes output in standard format for report and display methods.

    Args:
        stream: stream to write to
        lines_list: list containing lines to be written in body of report
        title: title to be put at top of report
        line_if_empty: line to be written if lines_list is empty
        end_line: line to be written at end of report
        header: character to use to write header separation line
        footer: character to use to write footer separation line

    Returns:
        None

    """
    stream.write(f"{header * MAX_STR_LENGTH}\n")
    if title is not None:
        stream.write(f"{title}\n\n")
    if len(lines_list) > 0:
        for i in lines_list:
            stream.write(f"{TAB}{i}\n")
    elif line_if_empty is not None:
        stream.write(f"{TAB}{line_if_empty}\n")
    stream.write("\n")
    if end_line is not None:
        stream.write(f"{end_line}\n")
    if footer is not None:
        stream.write(f"{footer * MAX_STR_LENGTH}\n")


def _collect_model_statistics(model):
    vars_in_constraints = variables_in_activated_constraints_set(model)
    fixed_vars_in_constraints = ComponentSet()
    free_vars_in_constraints = ComponentSet()
    free_vars_lb = ComponentSet()
    free_vars_ub = ComponentSet()
    free_vars_lbub = ComponentSet()
    ext_fixed_vars_in_constraints = ComponentSet()
    ext_free_vars_in_constraints = ComponentSet()
    for v in vars_in_constraints:
        if v.fixed:
            fixed_vars_in_constraints.add(v)
            if not _var_in_block(v, model):
                ext_fixed_vars_in_constraints.add(v)
        else:
            free_vars_in_constraints.add(v)
            if not _var_in_block(v, model):
                ext_free_vars_in_constraints.add(v)
            if v.lb is not None:
                if v.ub is not None:
                    free_vars_lbub.add(v)
                else:
                    free_vars_lb.add(v)
            elif v.ub is not None:
                free_vars_ub.add(v)

    # Generate report
    # TODO: Binary and boolean vars
    stats = []
    stats.append(
        f"{TAB}Activated Blocks: {len(activated_blocks_set(model))} "
        f"(Deactivated: {len(deactivated_blocks_set(model))})"
    )
    stats.append(
        f"{TAB}Free Variables in Activated Constraints: "
        f"{len(free_vars_in_constraints)} "
        f"(External: {len(ext_free_vars_in_constraints)})"
    )
    stats.append(f"{TAB * 2}Free Variables with only lower bounds: {len(free_vars_lb)}")
    stats.append(f"{TAB * 2}Free Variables with only upper bounds: {len(free_vars_ub)}")
    stats.append(
        f"{TAB * 2}Free Variables with upper and lower bounds: "
        f"{len(free_vars_lbub)}"
    )
    stats.append(
        f"{TAB}Fixed Variables in Activated Constraints: "
        f"{len(fixed_vars_in_constraints)} "
        f"(External: {len(ext_fixed_vars_in_constraints)})"
    )
    stats.append(
        f"{TAB}Activated Equality Constraints: {len(activated_equalities_set(model))} "
        f"(Deactivated: {len(deactivated_equalities_set(model))})"
    )
    stats.append(
        f"{TAB}Activated Inequality Constraints: {len(activated_inequalities_set(model))} "
        f"(Deactivated: {len(deactivated_inequalities_set(model))})"
    )
    stats.append(
        f"{TAB}Activated Objectives: {len(activated_objectives_set(model))} "
        f"(Deactivated: {len(deactivated_objectives_set(model))})"
    )

    return stats
