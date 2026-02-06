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
This module contains a collection of tools for diagnosing modeling issues.
"""

__author__ = "Alexander Dowling, Douglas Allan, Andrew Lee, Robby Parker, Ben Knueven"

import sys
from math import log
from typing import List
import logging

import numpy as np
from scipy.sparse.linalg import norm

from pyomo.environ import (
    Constraint,
    Objective,
    Param,
    SolverFactory,
    value,
    Var,
)
from pyomo.core.base.block import BlockData
from pyomo.core.base.constraint import ConstraintData
from pyomo.core.base.expression import ExpressionData
from pyomo.common.collections import ComponentSet
from pyomo.common.config import (
    ConfigDict,
    ConfigValue,
    document_kwargs_from_configdict,
    NonNegativeFloat,
    NonNegativeInt,
)
from pyomo.util.check_units import identify_inconsistent_units
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.contrib.iis import mis

from idaes.core.solvers.get_solver import get_solver
from idaes.core.util.misc import compact_expression_to_string
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
    variables_with_none_value_in_activated_equalities_set,
    number_activated_greybox_equalities,
    number_deactivated_greybox_equalities,
    activated_greybox_block_set,
    deactivated_greybox_block_set,
    greybox_block_set,
    unfixed_greybox_variables,
    greybox_variables,
    degrees_of_freedom,
    large_residuals_set,
    variables_near_bounds_set,
)
from idaes.core.scaling.util import (
    get_jacobian,
    get_scaling_factor,
    jacobian_cond,
)
import idaes.logger as idaeslog

from idaes.core.util.diagnostics_tools.evaluation_error import (
    EvalErrorWalker as _EvalErrorWalker,
)
from idaes.core.util.diagnostics_tools.writer_utils import (
    write_report_section as _write_report_section,
    MAX_STR_LENGTH,
    TAB,
)


# Imports for backward compatibility
from idaes.core.util.diagnostics_tools import (
    ConstraintTermAnalysisVisitor,
    compute_ill_conditioning_certificate,
    # TODO: Rename and redirect once old DegeneracyHunter is removed
    DegeneracyHunter as DegeneracyHunter2,
    DHCONFIG,
    # IpoptConvergenceAnalysis,
    SVDToolbox,
    SVDCONFIG,
)
from idaes.core.util.diagnostics_tools.deprecated.degeneracy_hunter_legacy import DegeneracyHunter


_log = idaeslog.getLogger(__name__)


# TODO: Add suggested steps to cautions - how?


CONFIG = ConfigDict()
CONFIG.declare(
    "variable_bounds_absolute_tolerance",
    ConfigValue(
        default=1e-4,
        domain=NonNegativeFloat,
        description="Absolute tolerance for considering a variable to be close "
        "to its bounds.",
    ),
)
# TODO is a relative tolerance necessary if we're including scaling?
CONFIG.declare(
    "variable_bounds_relative_tolerance",
    ConfigValue(
        default=1e-4,
        domain=NonNegativeFloat,
        description="Relative tolerance for considering a variable to be close "
        "to its bounds.",
    ),
)
CONFIG.declare(
    "variable_bounds_violation_tolerance",
    ConfigValue(
        default=0,
        domain=NonNegativeFloat,
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
        domain=NonNegativeFloat,
        description="Absolute tolerance to use when checking constraint residuals.",
    ),
)
CONFIG.declare(
    "constraint_term_mismatch_tolerance",
    ConfigValue(
        default=1e6,
        domain=NonNegativeFloat,
        description="Magnitude difference to use when checking for mismatched additive terms in constraints.",
    ),
)
CONFIG.declare(
    "constraint_term_cancellation_tolerance",
    ConfigValue(
        default=1e-4,
        domain=NonNegativeFloat,
        description="Absolute tolerance to use when checking for canceling additive terms in constraints.",
    ),
)
CONFIG.declare(
    "max_canceling_terms",
    ConfigValue(
        default=5,
        domain=NonNegativeInt,
        description="Maximum number of terms to consider when looking for canceling combinations in expressions.",
    ),
)
CONFIG.declare(
    "constraint_term_zero_tolerance",
    ConfigValue(
        default=1e-10,
        domain=NonNegativeFloat,
        description="Absolute tolerance to use when determining if a constraint term is equal to zero.",
    ),
)
CONFIG.declare(
    "variable_large_value_tolerance",
    ConfigValue(
        default=1e4,
        domain=NonNegativeFloat,
        description="Absolute tolerance for considering a value to be large.",
    ),
)
CONFIG.declare(
    "variable_small_value_tolerance",
    ConfigValue(
        default=1e-4,
        domain=NonNegativeFloat,
        description="Absolute tolerance for considering a value to be small.",
    ),
)
CONFIG.declare(
    "variable_zero_value_tolerance",
    ConfigValue(
        default=1e-8,
        domain=NonNegativeFloat,
        description="Absolute tolerance for considering a value to be near to zero.",
    ),
)
CONFIG.declare(
    "jacobian_large_value_caution",
    ConfigValue(
        default=1e4,
        domain=NonNegativeFloat,
        description="Tolerance for raising a caution for large Jacobian values.",
    ),
)
CONFIG.declare(
    "jacobian_large_value_warning",
    ConfigValue(
        default=1e8,
        domain=NonNegativeFloat,
        description="Tolerance for raising a warning for large Jacobian values.",
    ),
)
CONFIG.declare(
    "jacobian_small_value_caution",
    ConfigValue(
        default=1e-4,
        domain=NonNegativeFloat,
        description="Tolerance for raising a caution for small Jacobian values.",
    ),
)
CONFIG.declare(
    "jacobian_small_value_warning",
    ConfigValue(
        default=1e-8,
        domain=NonNegativeFloat,
        description="Tolerance for raising a warning for small Jacobian values.",
    ),
)
CONFIG.declare(
    "warn_for_evaluation_error_at_bounds",
    ConfigValue(
        default=True,
        domain=bool,
        description="If False, warnings will not be generated for things like log(x) with x >= 0",
    ),
)
CONFIG.declare(
    "parallel_component_tolerance",
    ConfigValue(
        default=1e-8,
        domain=NonNegativeFloat,
        description="Tolerance for identifying near-parallel Jacobian rows/columns",
    ),
)
CONFIG.declare(
    "absolute_feasibility_tolerance",
    ConfigValue(
        default=1e-6,
        domain=NonNegativeFloat,
        description="Feasibility tolerance for identifying infeasible constraints and bounds",
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

    def __init__(self, model: BlockData, **kwargs):
        # TODO: In future may want to generalise this to accept indexed blocks
        # However, for now some of the tools do not support indexed blocks
        if not isinstance(model, BlockData):
            raise TypeError(
                "model argument must be an instance of a Pyomo BlockData object "
                "(either a scalar Block or an element of an indexed Block)."
            )
        if len(greybox_block_set(model)) != 0:
            raise NotImplementedError(
                "Model contains Greybox models, which are not supported by Diagnostics toolbox at the moment"
            )
        self._model = model
        self.config = CONFIG(kwargs)

    @property
    def model(self):
        """
        Model currently being diagnosed.
        """
        return self._model

    def display_external_variables(self, stream=None):
        """
        Prints a list of variables that appear within activated Constraints in the
        model but are not contained within the model themselves.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        if stream is None:
            stream = sys.stdout

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

    def display_unused_variables(self, stream=None):
        """
        Prints a list of variables that do not appear in any activated Constraints.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        if stream is None:
            stream = sys.stdout

        _write_report_section(
            stream=stream,
            lines_list=variables_not_in_activated_constraints_set(self._model),
            title="The following variable(s) do not appear in any activated constraints within the model:",
            header="=",
            footer="=",
        )

    def display_variables_fixed_to_zero(self, stream=None):
        """
        Prints a list of variables that are fixed to an absolute value of 0.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        if stream is None:
            stream = sys.stdout

        _write_report_section(
            stream=stream,
            lines_list=_vars_fixed_to_zero(self._model),
            title="The following variable(s) are fixed to zero:",
            header="=",
            footer="=",
        )

    def display_variables_at_or_outside_bounds(self, stream=None):
        """
        Prints a list of variables with values that fall at or outside the bounds
        on the variable.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        if stream is None:
            stream = sys.stdout

        _write_report_section(
            stream=stream,
            lines_list=[
                f"{v.name} ({'fixed' if v.fixed else 'free'}): value={value(v)} bounds={v.bounds}"
                for v in _vars_violating_bounds(
                    self._model,
                    tolerance=self.config.variable_bounds_violation_tolerance,
                )
            ],
            title="The following variable(s) have values at or outside their bounds "
            f"(tol={self.config.variable_bounds_violation_tolerance:.1E}):",
            header="=",
            footer="=",
        )

    def display_variables_with_none_value(self, stream=None):
        """
        Prints a list of variables with a value of None.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        if stream is None:
            stream = sys.stdout

        _write_report_section(
            stream=stream,
            lines_list=_vars_with_none_value(self._model),
            title="The following variable(s) have a value of None:",
            header="=",
            footer="=",
        )

    def display_variables_with_none_value_in_activated_constraints(self, stream=None):
        """
        Prints a list of variables with values of None that are present in the
        mathematical program generated to solve the model. This list includes only
        variables in active constraints that are reachable through active blocks.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        if stream is None:
            stream = sys.stdout

        _write_report_section(
            stream=stream,
            lines_list=[
                f"{v.name}"
                for v in variables_with_none_value_in_activated_equalities_set(
                    self._model
                )
            ],
            title="The following variable(s) have a value of None:",
            header="=",
            footer="=",
        )

    def _verify_active_variables_initialized(self, stream=None):
        """
        Validate that all variables are initialized (i.e., have values set to
        something other than None) before doing further numerical analysis.
        Stream argument provided for forward compatibility (in case we want
        to print a list or something).
        """
        n_uninit = len(
            variables_with_none_value_in_activated_equalities_set(self._model)
        )
        if n_uninit > 0:
            raise RuntimeError(
                f"Found {n_uninit} variables with a value of None in the mathematical "
                "program generated by the model. They must be initialized with non-None "
                "values before numerical analysis can proceed. Run "
                + self.display_variables_with_none_value_in_activated_constraints.__name__
                + " to display a list of these variables."
            )

    def display_variables_with_value_near_zero(self, stream=None):
        """
        Prints a list of variables with a value close to zero. The tolerance
        for determining what is close to zero can be set in the class configuration
        options.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        if stream is None:
            stream = sys.stdout

        _write_report_section(
            stream=stream,
            lines_list=[
                f"{v.name}: value={value(v)}"
                for v in _vars_near_zero(
                    self._model, self.config.variable_zero_value_tolerance
                )
            ],
            title=f"The following variable(s) have a value close to zero "
            f"(tol={self.config.variable_zero_value_tolerance:.1E}):",
            header="=",
            footer="=",
        )

    def display_variables_with_extreme_values(self, stream=None):
        """
        Prints a list of variables with extreme values.

        Tolerances can be set in the class configuration options.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        if stream is None:
            stream = sys.stdout

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
            title=f"The following variable(s) have extreme values "
            f"(<{self.config.variable_small_value_tolerance:.1E} or "
            f"> {self.config.variable_large_value_tolerance:.1E}):",
            header="=",
            footer="=",
        )

    def display_variables_near_bounds(self, stream=None):
        """
        Prints a list of variables with values close to their bounds. Tolerance can
        be set in the class configuration options.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        if stream is None:
            stream = sys.stdout

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
            title=f"The following variable(s) have values close to their bounds "
            f"(abs={self.config.variable_bounds_absolute_tolerance:.1E}, "
            f"rel={self.config.variable_bounds_relative_tolerance:.1E}):",
            header="=",
            footer="=",
        )

    def display_components_with_inconsistent_units(self, stream=None):
        """
        Prints a list of all Constraints, Expressions and Objectives in the
        model with inconsistent units of measurement.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        if stream is None:
            stream = sys.stdout

        _write_report_section(
            stream=stream,
            lines_list=identify_inconsistent_units(self._model),
            title="The following component(s) have unit consistency issues:",
            end_line="For more details on unit inconsistencies, import the "
            "assert_units_consistent method\nfrom pyomo.util.check_units",
            header="=",
            footer="=",
        )

    def display_constraints_with_large_residuals(self, stream=None):
        """
        Prints a list of Constraints with residuals greater than a specified tolerance.
        Tolerance can be set in the class configuration options.

        Args:
            stream: an I/O object to write the list to (default = stdout)

        Returns:
            None

        """
        if stream is None:
            stream = sys.stdout

        lrdict = large_residuals_set(
            self._model,
            tol=self.config.constraint_residual_tolerance,
            return_residual_values=True,
        )

        lrs = []
        for k, v in lrdict.items():
            lrs.append(f"{k.name}: {v:.5E}")

        _write_report_section(
            stream=stream,
            lines_list=lrs,
            title=f"The following constraint(s) have large residuals "
            f"(>{self.config.constraint_residual_tolerance:.1E}):",
            header="=",
            footer="=",
        )

    def compute_infeasibility_explanation(self, stream=None, solver=None, tee=False):
        """
        This function attempts to determine why a given model is infeasible. It deploys
        two main algorithms:

        1. Relaxes the constraints of the problem, and reports to the user
           some sets of constraints and variable bounds, which when relaxed, creates a
           feasible model.
        2. Uses the information collected from (1) to attempt to compute a Minimal
           Infeasible System (MIS), which is a set of constraints and variable bounds
           which appear to be in conflict with each other. It is minimal in the sense
           that removing any single constraint or variable bound would result in a
           feasible subsystem.

        Args:
            stream: I/O object to write report to (default = stdout)
            solver: A pyomo solver object or a string for SolverFactory
                (default = get_solver())
            tee:  Display intermediate solves conducted (False)

        Returns:
            None

        """
        if solver is None:
            solver = get_solver("ipopt_v2")
        if stream is None:
            stream = sys.stdout

        h = logging.StreamHandler(stream)
        h.setLevel(logging.INFO)

        l = logging.Logger(name=__name__ + ".compute_infeasibility_explanation")
        l.setLevel(logging.INFO)
        l.addHandler(h)

        mis.compute_infeasibility_explanation(
            self._model,
            solver,
            tee=tee,
            tolerance=self.config.absolute_feasibility_tolerance,
            logger=l,
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

    def display_underconstrained_set(self, stream=None):
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
        if stream is None:
            stream = sys.stdout

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

    def display_overconstrained_set(self, stream=None):
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
        if stream is None:
            stream = sys.stdout

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

    def display_variables_with_extreme_jacobians(self, stream=None):
        """
        Prints the variables corresponding to columns in the Jacobian with extreme
        L2 norms. This often indicates poorly scaled variables.

        Tolerances can be set via the DiagnosticsToolbox config.

        Args:
            stream: an I/O object to write the output to (default = stdout)

        Returns:
            None

        """
        self._verify_active_variables_initialized(stream=stream)

        if stream is None:
            stream = sys.stdout

        jac, nlp = get_jacobian(self._model)

        xjc = _extreme_jacobian_columns(
            jac=jac,
            nlp=nlp,
            large=self.config.jacobian_large_value_caution,
            small=self.config.jacobian_small_value_caution,
        )
        xjc.sort(key=lambda i: abs(log(i[0])), reverse=True)

        _write_report_section(
            stream=stream,
            lines_list=[f"{i[1].name}: {i[0]:.3E}" for i in xjc],
            title=f"The following variable(s) correspond to Jacobian columns with extreme norms"
            f"(<{self.config.jacobian_small_value_caution:.1E} or"
            f">{self.config.jacobian_large_value_caution:.1E}):",
            header="=",
            footer="=",
        )

    def display_constraints_with_extreme_jacobians(self, stream=None):
        """
        Prints the constraints corresponding to rows in the Jacobian with extreme
        L2 norms. This often indicates poorly scaled constraints.

        Tolerances can be set via the DiagnosticsToolbox config.

        Args:
            stream: an I/O object to write the output to (default = stdout)

        Returns:
            None

        """
        self._verify_active_variables_initialized(stream=stream)

        if stream is None:
            stream = sys.stdout

        jac, nlp = get_jacobian(self._model)

        xjr = _extreme_jacobian_rows(
            jac=jac,
            nlp=nlp,
            large=self.config.jacobian_large_value_caution,
            small=self.config.jacobian_small_value_caution,
        )
        xjr.sort(key=lambda i: abs(log(i[0])), reverse=True)

        _write_report_section(
            stream=stream,
            lines_list=[f"{i[1].name}: {i[0]:.3E}" for i in xjr],
            title="The following constraint(s) correspond to Jacobian rows with extreme norms "
            f"(<{self.config.jacobian_small_value_caution:.1E} or"
            f">{self.config.jacobian_large_value_caution:.1E}):",
            header="=",
            footer="=",
        )

    def display_extreme_jacobian_entries(self, stream=None):
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
        self._verify_active_variables_initialized(stream=stream)

        if stream is None:
            stream = sys.stdout

        jac, nlp = get_jacobian(self._model, include_scaling_factors=True)
        xje = _extreme_jacobian_entries(
            jac,
            nlp,
            large=self.config.jacobian_large_value_caution,
            small=self.config.jacobian_small_value_caution,
            zero=0,
        )
        xje.sort(key=lambda i: abs(log(i[0])), reverse=True)

        _write_report_section(
            stream=stream,
            lines_list=[f"{i[1].name}, {i[2].name}: {i[0]:.3E}" for i in xje],
            title="The following constraint(s) and variable(s) are associated with extreme "
            f"Jacobian\nentries (<{self.config.jacobian_small_value_caution:.1E} or"
            f">{self.config.jacobian_large_value_caution:.1E}):",
            header="=",
            footer="=",
        )

    def display_near_parallel_constraints(self, stream=None):
        """
        Display near-parallel (duplicate) constraints in model.

        Args:
            stream: I/O object to write report to (default = stdout)

        Returns:
            None

        """
        self._verify_active_variables_initialized(stream=stream)

        if stream is None:
            stream = sys.stdout

        parallel = [
            f"{i[0].name}, {i[1].name}"
            for i in check_parallel_jacobian(
                model=self._model,
                tolerance=self.config.parallel_component_tolerance,
                direction="row",
            )
        ]

        # Write the output
        _write_report_section(
            stream=stream,
            lines_list=parallel,
            title="The following pairs of constraints are nearly parallel:",
            header="=",
            footer="=",
        )

    def display_near_parallel_variables(self, stream=None):
        """
        Display near-parallel (duplicate) variables in model.

        Args:
            stream: I/O object to write report to (default = stdout)

        Returns:
            None

        """
        self._verify_active_variables_initialized(stream=stream)

        if stream is None:
            stream = sys.stdout

        parallel = [
            f"{i[0].name}, {i[1].name}"
            for i in check_parallel_jacobian(
                model=self._model,
                tolerance=self.config.parallel_component_tolerance,
                direction="column",
            )
        ]

        # Write the output
        _write_report_section(
            stream=stream,
            lines_list=parallel,
            title="The following pairs of variables are nearly parallel:",
            header="=",
            footer="=",
        )

    # TODO: Block triangularization analysis
    # Number and size of blocks, polynomial degree of 1x1 blocks, simple pivot check of moderate sized sub-blocks?

    def _collect_constraint_mismatches(self, descend_into=True):
        """
        Call ConstraintTermAnalysisVisitor on all Constraints in model to walk expression
        tree and collect any instances of sum expressions with mismatched terms or potential
        cancellations.

        Args:
            descend_into: whether to descend_into child_blocks

        Returns:
            List of strings summarising constraints with mismatched terms
            List of strings summarising constraints with cancellations
            List of strings with constraint names where constraint contains no free variables
        """
        walker = ConstraintTermAnalysisVisitor(
            term_mismatch_tolerance=self.config.constraint_term_mismatch_tolerance,
            term_cancellation_tolerance=self.config.constraint_term_cancellation_tolerance,
            term_zero_tolerance=self.config.constraint_term_zero_tolerance,
            # for the high level summary, we only need to know if there are any cancellations,
            # but don't need to find all of them
            max_cancellations_per_node=1,
            max_canceling_terms=self.config.max_canceling_terms,
        )

        mismatch = []
        cancellation = []
        constant = []

        for c in self._model.component_data_objects(
            Constraint, descend_into=descend_into
        ):
            _, expr_mismatch, expr_cancellation, expr_constant, _ = (
                walker.walk_expression(c.expr)
            )

            if len(expr_mismatch) > 0:
                mismatch.append(f"{c.name}: {len(expr_mismatch)} mismatched term(s)")

            if len(expr_cancellation) > 0:
                cancellation.append(
                    f"{c.name}: {len(expr_cancellation)} potential canceling term(s)"
                )

            if expr_constant:
                constant.append(c.name)

        return mismatch, cancellation, constant

    def display_constraints_with_mismatched_terms(self, stream=None):
        """
        Display constraints in model which contain additive terms of significantly different magnitude.

        Args:
            stream: I/O object to write report to (default = stdout)

        Returns:
            None

        """
        self._verify_active_variables_initialized(stream=stream)

        if stream is None:
            stream = sys.stdout

        mismatch, _, _ = self._collect_constraint_mismatches()

        # Write the output
        _write_report_section(
            stream=stream,
            lines_list=mismatch,
            title="The following constraints have mismatched terms:",
            end_line="Call display_problematic_constraint_terms(constraint) for more information.",
            header="=",
            footer="=",
        )

    def display_constraints_with_canceling_terms(self, stream=None):
        """
        Display constraints in model which contain additive terms which potentially cancel each other.

        Note that this method looks a the current state of the constraint, and will flag terms as
        cancelling if you have a form A == B + C where C is significantly smaller than A and B. In some
        cases this behavior is intended, as C is a correction term which happens to be very
        small at the current state. However, you should review these constraints to determine whether
        the correction term is important for the situation you are modeling and consider removing the
        term if it will never be significant.

        Args:
            stream: I/O object to write report to (default = stdout)

        Returns:
            None

        """
        self._verify_active_variables_initialized(stream=stream)

        if stream is None:
            stream = sys.stdout

        _, canceling, _ = self._collect_constraint_mismatches()

        # Write the output
        _write_report_section(
            stream=stream,
            lines_list=canceling,
            title="The following constraints have canceling terms:",
            end_line="Call display_problematic_constraint_terms(constraint) for more information.",
            header="=",
            footer="=",
        )

    def display_problematic_constraint_terms(
        self, constraint, max_cancellations: int = 5, stream=None
    ):
        """
        Display a summary of potentially problematic terms in a given constraint.

        Note that this method looks a the current state of the constraint, and will flag terms as
        cancelling if you have a form A == B + C where C is significantly smaller than A and B. In some
        cases this behavior is intended, as C is a correction term which happens to be very
        small at the current state. However, you should review these constraints to determine whether
        the correction term is important for the situation you are modeling and consider removing the
        term if it will never be significant.

        Args:
            constraint: ConstraintData object to be examined
            max_cancellations: maximum number of cancellations per node before termination.
                None = find all cancellations.
            stream: I/O object to write report to (default = stdout)

        Returns:
            None

        """
        self._verify_active_variables_initialized(stream=stream)

        if stream is None:
            stream = sys.stdout

        # Check that constraint is of correct type to give useful error message
        if not isinstance(constraint, ConstraintData):
            # Wrong type, check if it is an indexed constraint
            if isinstance(constraint, Constraint):
                raise TypeError(
                    f"{constraint.name} is an IndexedConstraint. Please provide "
                    f"an individual element of {constraint.name} (ConstraintData) "
                    "to be examined for problematic terms."
                )
            else:
                # Not a constraint
                raise TypeError(
                    f"{constraint.name} is not an instance of a Pyomo Constraint."
                )
        sf = get_scaling_factor(constraint, default=1, warning=False)

        # Don't need to scale constraint_term_mismatch_tolerance and
        # constraint_term_cancellation_tolerance because they are both
        # relative tolerances. term_zero_tolerance is an absolute
        # tolerance, so it must be scaled.
        walker = ConstraintTermAnalysisVisitor(
            term_mismatch_tolerance=self.config.constraint_term_mismatch_tolerance,
            term_cancellation_tolerance=self.config.constraint_term_cancellation_tolerance,
            term_zero_tolerance=sf * self.config.constraint_term_zero_tolerance,
            max_cancellations_per_node=max_cancellations,
            max_canceling_terms=self.config.max_canceling_terms,
        )

        _, expr_mismatch, expr_cancellation, _, tripped = walker.walk_expression(
            constraint.expr
        )

        # Combine mismatches and cancellations into a summary list
        issues = []
        for k, v in expr_mismatch.items():
            tag = " "
            if isinstance(k, ExpressionData):
                # For clarity, if the problem node is a named Expression, note this in output
                tag = " Expression "
            # Want to show full expression node plus largest and smallest magnitudes
            issues.append(
                f"Mismatch in{tag}{compact_expression_to_string(k)} (Max {v[0]}, Min {v[1]})"
            )
        # Collect summary of cancelling terms for user
        # Walker gives us back a list of nodes with cancelling terms
        for k, v in expr_cancellation.items():
            # Each node may have multiple cancellations, these are given as separate tuples
            tag = " "
            if isinstance(k, ExpressionData):
                # For clarity, if the problem node is a named Expression, note this in output
                tag = " Expression "
            for i in v:
                # For each cancellation, iterate over contributing terms and write a summary
                terms = ""
                for j in i:
                    if len(terms) > 0:
                        terms += ", "
                    # +1 to switch from 0-index to 1-index
                    terms += f"{j[0]+1} ({j[1]})"
                issues.append(
                    f"Cancellation in{tag}{compact_expression_to_string(k)}. Terms {terms}"
                )

        if tripped:
            end_line = (
                f"Number of canceling terms per node limited to {max_cancellations}."
            )
        else:
            end_line = None

        # Write the output
        _write_report_section(
            stream=stream,
            lines_list=issues,
            title=f"The following terms in {constraint.name} are potentially problematic:",
            end_line=end_line,
            header="=",
            footer="=",
        )

    def display_constraints_with_no_free_variables(self, stream=None):
        """
        Display constraints in model which contain no free variables.

        Args:
            stream: I/O object to write report to (default = stdout)

        Returns:
            None

        """
        # Although, in principle, this method doesn't require
        # all variables to be initialized, its current
        # implementation does.
        self._verify_active_variables_initialized(stream=stream)

        if stream is None:
            stream = sys.stdout

        _, _, constant = self._collect_constraint_mismatches()

        # Write the output
        _write_report_section(
            stream=stream,
            lines_list=constant,
            title="The following constraints have no free variables:",
            header="=",
            footer="=",
        )

    def _collect_structural_warnings(
        self, ignore_evaluation_errors=False, ignore_unit_consistency=False
    ):
        """
        Runs checks for structural warnings and returns two lists.

        Args:
            ignore_evaluation_errors - ignore potential evaluation error warnings
            ignore_unit_consistency - ignore unit consistency warnings

        Returns:
            warnings - list of warning messages from structural analysis
            next_steps - list of suggested next steps to further investigate warnings

        """
        if not ignore_unit_consistency:
            uc = identify_inconsistent_units(self._model)
        else:
            uc = []
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

        if not ignore_evaluation_errors:
            eval_warnings = self._collect_potential_eval_errors()
            if len(eval_warnings) > 0:
                warnings.append(
                    f"WARNING: Found {len(eval_warnings)} potential evaluation errors."
                )
                next_steps.append(
                    self.display_potential_evaluation_errors.__name__ + "()"
                )

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

    def _collect_numerical_warnings(
        self, jac=None, nlp=None, ignore_parallel_components=False
    ):
        """
        Runs checks for numerical warnings and returns two lists.

        Args:
            ignore_parallel_components - ignore checks for parallel components

        Returns:
            warnings - list of warning messages from numerical analysis
            next_steps - list of suggested next steps to further investigate warnings

        """
        if jac is None or nlp is None:
            jac, nlp = get_jacobian(self._model)

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
                f"WARNING: {len(large_residuals)} {cstring} with large residuals "
                f"(>{self.config.constraint_residual_tolerance:.1E})"
            )
            next_steps.append(
                self.display_constraints_with_large_residuals.__name__ + "()"
            )
            next_steps.append(self.compute_infeasibility_explanation.__name__ + "()")

        # Variables outside bounds
        violated_bounds = _vars_violating_bounds(
            self._model, tolerance=self.config.variable_bounds_violation_tolerance
        )
        if len(violated_bounds) > 0:
            cstring = "Variables"
            if len(violated_bounds) == 1:
                cstring = "Variable"
            warnings.append(
                f"WARNING: {len(violated_bounds)} {cstring} at or outside bounds "
                f"(tol={self.config.variable_bounds_violation_tolerance:.1E})"
            )
            next_steps.append(
                self.display_variables_at_or_outside_bounds.__name__ + "()"
            )

        # Extreme Jacobian rows and columns
        jac_col = _extreme_jacobian_columns(
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
                f"WARNING: {len(jac_col)} {cstring} with extreme Jacobian column norms "
                f"(<{self.config.jacobian_small_value_warning:.1E} or "
                f">{self.config.jacobian_large_value_warning:.1E})"
            )
            next_steps.append(
                self.display_variables_with_extreme_jacobians.__name__ + "()"
            )

        jac_row = _extreme_jacobian_rows(
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
                f"WARNING: {len(jac_row)} {cstring} with extreme Jacobian row norms "
                f"(<{self.config.jacobian_small_value_warning:.1E} or "
                f">{self.config.jacobian_large_value_warning:.1E})"
            )
            next_steps.append(
                self.display_constraints_with_extreme_jacobians.__name__ + "()"
            )

        # Parallel variables and constraints
        if not ignore_parallel_components:
            partol = self.config.parallel_component_tolerance
            par_cons = check_parallel_jacobian(
                self._model, tolerance=partol, direction="row", jac=jac, nlp=nlp
            )
            par_vars = check_parallel_jacobian(
                self._model, tolerance=partol, direction="column", jac=jac, nlp=nlp
            )
            if par_cons:
                p = "pair" if len(par_cons) == 1 else "pairs"
                warnings.append(
                    f"WARNING: {len(par_cons)} {p} of constraints are parallel"
                    f" (to tolerance {partol:.1E})"
                )
                next_steps.append(
                    self.display_near_parallel_constraints.__name__ + "()"
                )
            if par_vars:
                p = "pair" if len(par_vars) == 1 else "pairs"
                warnings.append(
                    f"WARNING: {len(par_vars)} {p} of variables are parallel"
                    f" (to tolerance {partol:.1E})"
                )
                next_steps.append(self.display_near_parallel_variables.__name__ + "()")

        return warnings, next_steps

    def _collect_numerical_cautions(self, jac=None, nlp=None):
        """
        Runs checks for numerical cautions and returns a list.

        Returns:
            cautions - list of caution messages from numerical analysis

        """
        if jac is None or nlp is None:
            jac, nlp = get_jacobian(self._model)

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
                f"Caution: {len(near_bounds)} {cstring} with value close to their bounds "
                f"(abs={self.config.variable_bounds_absolute_tolerance:.1E}, "
                f"rel={self.config.variable_bounds_absolute_tolerance:.1E})"
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
                f"Caution: {len(near_zero)} {cstring} with value close to zero "
                f"(tol={self.config.variable_zero_value_tolerance:.1E})"
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
            cautions.append(
                f"Caution: {len(xval)} {cstring} with extreme value "
                f"(<{self.config.variable_small_value_tolerance:.1E} or "
                f">{self.config.variable_large_value_tolerance:.1E})"
            )

        # Variables with value None
        none_value = _vars_with_none_value(self._model)
        if len(none_value) > 0:
            cstring = "Variables"
            if len(none_value) == 1:
                cstring = "Variable"
            cautions.append(f"Caution: {len(none_value)} {cstring} with None value")

        # Constraints with possible ill-posed terms
        mismatch, cancellation, constant = self._collect_constraint_mismatches()
        if len(mismatch) > 0:
            cstring = "Constraints"
            if len(mismatch) == 1:
                cstring = "Constraint"
            cautions.append(f"Caution: {len(mismatch)} {cstring} with mismatched terms")
        if len(cancellation) > 0:
            cstring = "Constraints"
            if len(cancellation) == 1:
                cstring = "Constraint"
            cautions.append(
                f"Caution: {len(cancellation)} {cstring} with potential cancellation of terms"
            )
        if len(constant) > 0:
            cstring = "Constraints"
            if len(constant) == 1:
                cstring = "Constraint"
            cautions.append(
                f"Caution: {len(constant)} {cstring} with no free variables"
            )

        # Extreme Jacobian rows and columns
        jac_col = _extreme_jacobian_columns(
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
                f"Caution: {len(jac_col)} {cstring} with extreme Jacobian column norms "
                f"(<{self.config.jacobian_small_value_caution:.1E} or "
                f">{self.config.jacobian_large_value_caution:.1E})"
            )

        jac_row = _extreme_jacobian_rows(
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
                f"Caution: {len(jac_row)} {cstring} with extreme Jacobian row norms "
                f"(<{self.config.jacobian_small_value_caution:.1E} or "
                f">{self.config.jacobian_large_value_caution:.1E})"
            )

        # Extreme Jacobian entries
        extreme_jac = _extreme_jacobian_entries(
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
            cautions.append(
                f"Caution: {len(extreme_jac)} extreme Jacobian {cstring} "
                f"(<{self.config.jacobian_small_value_caution:.1E} or "
                f">{self.config.jacobian_large_value_caution:.1E})"
            )

        return cautions

    def assert_no_structural_warnings(
        self,
        ignore_evaluation_errors: bool = False,
        ignore_unit_consistency: bool = False,
    ):
        """
        Checks for structural warnings in the model and raises an AssertionError
        if any are found.

        Args:
            ignore_evaluation_errors - ignore potential evaluation error warnings
            ignore_unit_consistency - ignore unit consistency warnings

        Raises:
            AssertionError if any warnings are identified by structural analysis.

        """
        warnings, _ = self._collect_structural_warnings(
            ignore_evaluation_errors=ignore_evaluation_errors,
            ignore_unit_consistency=ignore_unit_consistency,
        )
        if len(warnings) > 0:
            raise AssertionError(f"Structural issues found ({len(warnings)}).")

    def assert_no_numerical_warnings(self, ignore_parallel_components=False):
        """
        Checks for numerical warnings in the model and raises an AssertionError
        if any are found.

        Args:
            ignore_parallel_components - ignore checks for parallel components

        Raises:
            AssertionError if any warnings are identified by numerical analysis.

        """
        warnings, _ = self._collect_numerical_warnings(
            ignore_parallel_components=ignore_parallel_components
        )
        if len(warnings) > 0:
            raise AssertionError(f"Numerical issues found ({len(warnings)}).")

    def report_structural_issues(self, stream=None):
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
        if stream is None:
            stream = sys.stdout

        # Potential evaluation errors
        # TODO: High Index?
        if len(greybox_block_set(self._model)) != 0:
            raise NotImplementedError(
                "Model contains Greybox models, which are not supported by Diagnostics toolbox at the moment"
            )
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

    def report_numerical_issues(self, stream=None):
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
        self._verify_active_variables_initialized(stream=stream)

        if stream is None:
            stream = sys.stdout
        jac, nlp = get_jacobian(self._model)

        warnings, next_steps = self._collect_numerical_warnings(jac=jac, nlp=nlp)
        cautions = self._collect_numerical_cautions(jac=jac, nlp=nlp)

        stats = []
        try:
            stats.append(
                f"Jacobian Condition Number: {jacobian_cond(jac=jac, scaled=True):.3E}"
            )
        except RuntimeError as err:
            if "Factor is exactly singular" in str(err):
                _log.info(err)
                stats.append("Jacobian Condition Number: Undefined (Exactly Singular)")
            else:
                raise

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
            f"\n{TAB*2}prepare_degeneracy_hunter()\n{TAB*2}prepare_svd_toolbox()",
            footer="=",
        )

    def _collect_potential_eval_errors(self) -> List[str]:
        warnings = list()
        for con in self.model.component_data_objects(
            Constraint, active=True, descend_into=True
        ):
            walker = _EvalErrorWalker(self.config)
            con_warnings = walker.walk_expression(con.body)
            for msg in con_warnings:
                msg = f"{con.name}: " + msg
                warnings.append(msg)
        for obj in self.model.component_data_objects(
            Objective, active=True, descend_into=True
        ):
            walker = _EvalErrorWalker(self.config)
            obj_warnings = walker.walk_expression(obj.expr)
            for msg in obj_warnings:
                msg = f"{obj.name}: " + msg
                warnings.append(msg)

        return warnings

    def display_potential_evaluation_errors(self, stream=None):
        """
        Prints constraints that may be prone to evaluation errors
        (e.g., log of a negative number) based on variable bounds.

        Args:
            stream: an I/O object to write the output to (default = stdout)

        Returns:
            None
        """
        if stream is None:
            stream = sys.stdout

        warnings = self._collect_potential_eval_errors()
        _write_report_section(
            stream=stream,
            lines_list=warnings,
            title=f"{len(warnings)} WARNINGS",
            line_if_empty="No warnings found!",
            header="=",
            footer="=",
        )

    @document_kwargs_from_configdict(SVDCONFIG)
    def prepare_svd_toolbox(self, **kwargs):
        """
        Create an instance of the SVDToolbox and store as self.svd_toolbox.

        After creating an instance of the toolbox, call
        display_underdetermined_variables_and_constraints().

        Returns:

            Instance of SVDToolbox

        """
        self.svd_toolbox = SVDToolbox(self.model, **kwargs)

        return self.svd_toolbox

    @document_kwargs_from_configdict(DHCONFIG)
    def prepare_degeneracy_hunter(self, **kwargs):
        """
        Create an instance of the DegeneracyHunter and store as self.degeneracy_hunter.

        After creating an instance of the toolbox, call
        report_irreducible_degenerate_sets.

        Returns:

            Instance of DegeneracyHunter

        """
        self.degeneracy_hunter = DegeneracyHunter2(self.model, **kwargs)

        return self.degeneracy_hunter


def _extreme_jacobian_entries(
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


def _extreme_jacobian_rows(
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


def _extreme_jacobian_columns(
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
    elif isinstance(component, BlockData):
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
    elif isinstance(component, BlockData):
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
        sf = get_scaling_factor(v, default=1, warning=False)
        if v.value is not None and sf * abs(value(v)) <= variable_zero_value_tolerance:
            near_zero_vars.add(v)
    return near_zero_vars


def _vars_violating_bounds(model, tolerance):
    violated_bounds = ComponentSet()
    for v in model.component_data_objects(Var, descend_into=True):
        sf = get_scaling_factor(v, default=1, warning=False)
        if v.value is not None:
            if v.lb is not None and sf * v.value <= sf * v.lb - tolerance:
                violated_bounds.add(v)
            elif v.ub is not None and sf * v.value >= sf * v.ub + tolerance:
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
        sf = get_scaling_factor(v, default=1, warning=False)
        if v.value is not None:
            mag = sf * abs(value(v))
            if mag > abs(large):
                extreme_vars.add(v)
            elif mag < abs(small) and mag > abs(zero):
                extreme_vars.add(v)

    return extreme_vars


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
        f"{TAB}Activated Equality Constraints: {len(activated_equalities_set(model))+number_activated_greybox_equalities(model)} "
        f"(Deactivated: {len(deactivated_equalities_set(model))+number_deactivated_greybox_equalities(model)})"
    )
    stats.append(
        f"{TAB}Activated Inequality Constraints: {len(activated_inequalities_set(model))} "
        f"(Deactivated: {len(deactivated_inequalities_set(model))})"
    )
    stats.append(
        f"{TAB}Activated Objectives: {len(activated_objectives_set(model))} "
        f"(Deactivated: {len(deactivated_objectives_set(model))})"
    )

    # Only show graybox info if they are present
    if len(greybox_block_set(model)) != 0:
        stats.append(f"{TAB}GreyBox Statistics")
        stats.append(
            f"{TAB* 2}Activated GreyBox models: {len(activated_greybox_block_set(model))} "
            f"(Deactivated: {len(deactivated_greybox_block_set(model))})"
        )
        stats.append(
            f"{TAB* 2}Activated GreyBox Equalities: {number_activated_greybox_equalities(model)} "
            f"(Deactivated: {number_deactivated_greybox_equalities(model)})"
        )
        stats.append(
            f"{TAB* 2}Free Variables in Activated GreyBox Equalities: {len(unfixed_greybox_variables(model))} (Fixed: {len(greybox_variables(model)-unfixed_greybox_variables(model))})"
        )

    return stats
