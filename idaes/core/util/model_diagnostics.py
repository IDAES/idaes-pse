# -*- coding: utf-8 -*-
#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
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

from operator import itemgetter
import sys
from inspect import signature
from math import log, log10, isclose, inf, isfinite
import json
from typing import List
import logging
from itertools import combinations, chain

import numpy as np
from scipy.linalg import svd
from scipy.sparse.linalg import svds, norm
from scipy.sparse import issparse, find

from pyomo.environ import (
    Binary,
    Integers,
    Block,
    check_optimal_termination,
    ComponentMap,
    ConcreteModel,
    Constraint,
    Expression,
    Objective,
    Param,
    RangeSet,
    Set,
    SolverFactory,
    value,
    Var,
)
from pyomo.core.expr.numeric_expr import (
    DivisionExpression,
    NPV_DivisionExpression,
    PowExpression,
    NPV_PowExpression,
    UnaryFunctionExpression,
    NPV_UnaryFunctionExpression,
    NumericExpression,
)
from pyomo.core.base.block import BlockData
from pyomo.core.base.var import VarData
from pyomo.core.base.constraint import ConstraintData
from pyomo.core.base.expression import ExpressionData
from pyomo.repn.standard_repn import (  # pylint: disable=no-name-in-module
    generate_standard_repn,
)
from pyomo.common.collections import ComponentSet
from pyomo.common.config import (
    ConfigDict,
    ConfigValue,
    document_kwargs_from_configdict,
    PositiveInt,
    NonNegativeFloat,
    NonNegativeInt,
)
from pyomo.util.check_units import identify_inconsistent_units
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.core.expr.visitor import identify_variables, StreamBasedExpressionVisitor
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
from pyomo.contrib.pynumero.asl import AmplInterface
from pyomo.contrib.fbbt.fbbt import compute_bounds_on_expr
from pyomo.contrib.iis import mis
from pyomo.common.deprecation import deprecation_warning
from pyomo.common.errors import PyomoException
from pyomo.common.tempfiles import TempfileManager
from pyomo.core import expr as EXPR
from pyomo.common.numeric_types import native_types
from pyomo.core.base.units_container import _PyomoUnit

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
from idaes.core.util.scaling import (
    extreme_jacobian_columns,
    extreme_jacobian_rows,
    extreme_jacobian_entries,
)
from idaes.core.scaling.util import (
    get_jacobian,
    get_scaling_factor,
    jacobian_cond,
)
from idaes.core.util.parameter_sweep import (
    SequentialSweepRunner,
    ParameterSweepBase,
    is_psweepspec,
)
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


MAX_STR_LENGTH = 84
TAB = " " * 4

# Constants for Degeneracy Hunter
YTOL = 0.9
MMULT = 0.99

# TODO: Add suggested steps to cautions - how?


def svd_callback_validator(val):
    """Domain validator for SVD callbacks

    Args:
        val : value to be checked

    Returns:
        TypeError if val is not a valid callback
    """
    if callable(val):
        sig = signature(val)
        if len(sig.parameters) >= 2:
            return val

    _log.error(
        f"SVD callback {val} must be a callable which takes at least two arguments."
    )
    raise ValueError(
        "SVD callback must be a callable which takes at least two arguments."
    )


def svd_dense(jacobian, number_singular_values):
    """
    Callback for performing SVD analysis using scipy.linalg.svd

    Args:
        jacobian: Jacobian to be analysed
        number_singular_values: number of singular values to compute

    Returns:
        u, s and v numpy arrays

    """
    u, s, vT = svd(jacobian.todense(), full_matrices=False)
    # Reorder singular values and vectors so that the singular
    # values are from least to greatest
    u = np.flip(u[:, -number_singular_values:], axis=1)
    s = np.flip(s[-number_singular_values:], axis=0)
    vT = np.flip(vT[-number_singular_values:, :], axis=0)

    return u, s, vT.transpose()


def svd_sparse(jacobian, number_singular_values):
    """
    Callback for performing SVD analysis using scipy.sparse.linalg.svds

    Args:
        jacobian: Jacobian to be analysed
        number_singular_values: number of singular values to compute

    Returns:
        u, s and v numpy arrays

    """
    u, s, vT = svds(jacobian, k=number_singular_values, which="SM")

    return u, s, vT.transpose()


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


SVDCONFIG = ConfigDict()
SVDCONFIG.declare(
    "number_of_smallest_singular_values",
    ConfigValue(
        domain=PositiveInt,
        description="Number of smallest singular values to compute",
    ),
)
SVDCONFIG.declare(
    "svd_callback",
    ConfigValue(
        default=svd_dense,
        domain=svd_callback_validator,
        description="Callback to SVD method of choice (default = svd_dense)",
        doc="Callback to SVD method of choice (default = svd_dense). "
        "Callbacks should take the Jacobian and number of singular values "
        "to compute as options, plus any method specific arguments, and should "
        "return the u, s and v matrices as numpy arrays.",
    ),
)
SVDCONFIG.declare(
    "svd_callback_arguments",
    ConfigValue(
        default=None,
        domain=dict,
        description="Optional arguments to pass to  SVD callback (default = None)",
    ),
)
SVDCONFIG.declare(
    "singular_value_tolerance",
    ConfigValue(
        default=1e-6,
        domain=NonNegativeFloat,
        description="Tolerance for defining a small singular value",
    ),
)
SVDCONFIG.declare(
    "size_cutoff_in_singular_vector",
    ConfigValue(
        default=0.1,
        domain=NonNegativeFloat,
        description="Size below which to ignore constraints and variables in "
        "the singular vector",
    ),
)


DHCONFIG = ConfigDict()
DHCONFIG.declare(
    "solver",
    ConfigValue(
        default="scip",
        domain=str,
        description="MILP solver to use for finding irreducible degenerate sets.",
    ),
)
DHCONFIG.declare(
    "solver_options",
    ConfigValue(
        domain=None,
        description="Options to pass to MILP solver.",
    ),
)
DHCONFIG.declare(
    "M",  # TODO: Need better name
    ConfigValue(
        default=1e5,
        domain=NonNegativeFloat,
        description="Maximum value for nu in MILP models.",
    ),
)
DHCONFIG.declare(
    "m_small",  # TODO: Need better name
    ConfigValue(
        default=1e-5,
        domain=NonNegativeFloat,
        description="Smallest value for nu to be considered non-zero in MILP models.",
    ),
)
DHCONFIG.declare(
    "trivial_constraint_tolerance",
    ConfigValue(
        default=1e-6,
        domain=NonNegativeFloat,
        description="Tolerance for identifying non-zero rows in Jacobian.",
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
        Prints the variables associated with columns in the Jacobian with extreme
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
            title=f"The following variable(s) are associated with extreme Jacobian values "
            f"(<{self.config.jacobian_small_value_caution:.1E} or"
            f">{self.config.jacobian_large_value_caution:.1E}):",
            header="=",
            footer="=",
        )

    def display_constraints_with_extreme_jacobians(self, stream=None):
        """
        Prints the constraints associated with rows in the Jacobian with extreme
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
            title="The following constraint(s) are associated with extreme Jacobian values "
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
            title="The following constraint(s) and variable(s) are associated with extreme "
            f"Jacobian\nvalues (<{self.config.jacobian_small_value_caution:.1E} or"
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
                f"WARNING: {len(jac_col)} {cstring} with extreme Jacobian values "
                f"(<{self.config.jacobian_small_value_warning:.1E} or "
                f">{self.config.jacobian_large_value_warning:.1E})"
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
                f"WARNING: {len(jac_row)} {cstring} with extreme Jacobian values "
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
                f"Caution: {len(jac_col)} {cstring} with extreme Jacobian values "
                f"(<{self.config.jacobian_small_value_caution:.1E} or "
                f">{self.config.jacobian_large_value_caution:.1E})"
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
                f"Caution: {len(jac_row)} {cstring} with extreme Jacobian values "
                f"(<{self.config.jacobian_small_value_caution:.1E} or "
                f">{self.config.jacobian_large_value_caution:.1E})"
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
                f"Jacobian Condition Number: {jacobian_cond(jac=jac, scaled=False):.3E}"
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


@document_kwargs_from_configdict(SVDCONFIG)
class SVDToolbox:
    """
    Toolbox for performing Singular Value Decomposition on the model Jacobian.

    Used help() for more information on available methods.

    Original code by Doug Allan

    Args:

        model: model to be diagnosed. The SVDToolbox does not support indexed Blocks.

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
        self.config = SVDCONFIG(kwargs)

        self.u = None
        self.s = None
        self.v = None

        # Get Jacobian and NLP
        self.jacobian, self.nlp = get_jacobian(
            self._model, equality_constraints_only=True
        )

        # Get list of equality constraint and variable names
        self._eq_con_list = self.nlp.get_pyomo_equality_constraints()
        self._var_list = self.nlp.get_pyomo_variables()

        if self.jacobian.shape[0] < 2:
            raise ValueError(
                "Model needs at least 2 equality constraints to perform svd_analysis."
            )

    def run_svd_analysis(self):
        """
        Perform SVD analysis of the constraint Jacobian

        Args:

            None

        Returns:

            None

        Actions:
            Stores SVD results in object

        """
        n_eq = self.jacobian.shape[0]
        n_var = self.jacobian.shape[1]

        n_sv = self.config.number_of_smallest_singular_values
        if n_sv is None:
            # Determine the number of singular values to compute
            # The "-1" is needed to avoid an error with svds
            n_sv = min(10, min(n_eq, n_var) - 1)
        elif n_sv >= min(n_eq, n_var):
            raise ValueError(
                f"For a {n_eq} by {n_var} system, svd_analysis "
                f"can compute at most {min(n_eq, n_var) - 1} "
                f"singular values and vectors, but {n_sv} were called for."
            )

        # Get optional arguments for SVD callback
        svd_callback_arguments = self.config.svd_callback_arguments
        if svd_callback_arguments is None:
            svd_callback_arguments = {}

        # Perform SVD
        # Recall J is a n_eq x n_var matrix
        # Thus U is a n_eq x n_eq matrix
        # And V is a n_var x n_var
        # (U or V may be smaller in economy mode)
        u, s, v = self.config.svd_callback(
            self.jacobian,
            number_singular_values=n_sv,
            **svd_callback_arguments,
        )

        # Save results
        self.u = u
        self.s = s
        self.v = v

    def display_rank_of_equality_constraints(self, stream=None):
        """
        Method to display the number of singular values that fall below
        tolerance specified in config block.

        Args:
            stream: I/O object to write report to (default = stdout)

        Returns:
            None

        """
        if stream is None:
            stream = sys.stdout

        if self.s is None:
            self.run_svd_analysis()

        counter = 0
        for e in self.s:
            if e < self.config.singular_value_tolerance:
                counter += 1

        stream.write("=" * MAX_STR_LENGTH + "\n\n")
        stream.write(
            f"Number of Singular Values less than "
            f"{self.config.singular_value_tolerance:.1E} is {counter}\n\n"
        )
        stream.write("=" * MAX_STR_LENGTH + "\n")

    def display_underdetermined_variables_and_constraints(
        self, singular_values=None, stream=None
    ):
        """
        Determines constraints and variables associated with the smallest
        singular values by having large components in the left and right
        singular vectors, respectively, associated with those singular values.

        Args:
            singular_values: List of ints representing singular values to display,
                as ordered from least to greatest starting from 1 (default show all)
            stream: I/O object to write report to (default = stdout)

        Returns:
            None

        """
        if stream is None:
            stream = sys.stdout

        if self.s is None:
            self.run_svd_analysis()

        tol = self.config.size_cutoff_in_singular_vector

        if singular_values is None:
            singular_values = range(1, len(self.s) + 1)

        stream.write("=" * MAX_STR_LENGTH + "\n")
        stream.write(
            "Constraints and Variables associated with smallest singular values\n\n"
        )

        for e in singular_values:
            # First, make sure values are feasible
            if e > len(self.s):
                raise ValueError(
                    f"Cannot display the {e}-th smallest singular value. "
                    f"Only {len(self.s)} small singular values have been "
                    "calculated. You can set the number_of_smallest_singular_values "
                    "config argument and call run_svd_analysis again to get more "
                    "singular values."
                )

            stream.write(f"{TAB}Smallest Singular Value {e}:\n\n")
            stream.write(f"{2 * TAB}Variables:\n\n")
            for v in np.where(abs(self.v[:, e - 1]) > tol)[0]:
                stream.write(f"{3 * TAB}{self._var_list[v].name}\n")

            stream.write(f"\n{2 * TAB}Constraints:\n\n")
            for c in np.where(abs(self.u[:, e - 1]) > tol)[0]:
                stream.write(f"{3 * TAB}{self._eq_con_list[c].name}\n")
            stream.write("\n")

        stream.write("=" * MAX_STR_LENGTH + "\n")

    def display_constraints_including_variable(self, variable, stream=None):
        """
        Display all constraints that include the specified variable and the
        associated Jacobian coefficient.

        Args:
            variable: variable object to get associated constraints for
            stream: I/O object to write report to (default = stdout)

        Returns:
            None

        """
        if stream is None:
            stream = sys.stdout

        # Validate variable argument
        if not isinstance(variable, VarData):
            raise TypeError(
                f"variable argument must be an instance of a Pyomo VarData "
                f"object (got {variable})."
            )

        # Get index of variable in Jacobian
        try:
            var_idx = self.nlp.get_primal_indices([variable])[0]
        except (KeyError, PyomoException):
            raise AttributeError(f"Could not find {variable.name} in model.")

        nonzeros = self.jacobian.getcol(var_idx).nonzero()

        # Build a list of all constraints that include var
        cons_w_var = []
        for r in nonzeros[0]:
            cons_w_var.append(
                f"{self._eq_con_list[r].name}: {self.jacobian[(r, var_idx)]:.3e}"
            )

        # Write the output
        _write_report_section(
            stream=stream,
            lines_list=cons_w_var,
            title=f"The following constraints involve {variable.name}:",
            header="=",
            footer="=",
        )

    def display_variables_in_constraint(self, constraint, stream=None):
        """
        Display all variables that appear in the specified constraint and the
        associated Jacobian coefficient.

        Args:
            constraint: constraint object to get associated variables for
            stream: I/O object to write report to (default = stdout)

        Returns:
            None

        """
        if stream is None:
            stream = sys.stdout

        # Validate variable argument
        if not isinstance(constraint, ConstraintData):
            raise TypeError(
                f"constraint argument must be an instance of a Pyomo ConstraintData "
                f"object (got {constraint})."
            )

        # Get index of variable in Jacobian
        try:
            con_idx = self.nlp.get_constraint_indices([constraint])[0]
        except KeyError:
            raise AttributeError(f"Could not find {constraint.name} in model.")

        nonzeros = self.jacobian[con_idx, :].nonzero()

        # Build a list of all vars in constraint
        vars_in_cons = []
        for c in nonzeros[1]:
            vars_in_cons.append(
                f"{self._var_list[c].name}: {self.jacobian[(con_idx, c)]:.3e}"
            )

        # Write the output
        _write_report_section(
            stream=stream,
            lines_list=vars_in_cons,
            title=f"The following variables are involved in {constraint.name}:",
            header="=",
            footer="=",
        )


def _get_bounds_with_inf(node: NumericExpression):
    lb, ub = compute_bounds_on_expr(node)
    if lb is None:
        lb = -inf
    if ub is None:
        ub = inf
    return lb, ub


def _check_eval_error_division(
    node: NumericExpression, warn_list: List[str], config: ConfigDict
):
    lb, ub = _get_bounds_with_inf(node.args[1])
    if (config.warn_for_evaluation_error_at_bounds and (lb <= 0 <= ub)) or (
        lb < 0 < ub
    ):
        msg = f"Potential division by 0 in {node}; Denominator bounds are ({lb}, {ub})"
        warn_list.append(msg)


def _check_eval_error_pow(
    node: NumericExpression, warn_list: List[str], config: ConfigDict
):
    arg1, arg2 = node.args
    lb1, ub1 = _get_bounds_with_inf(arg1)
    lb2, ub2 = _get_bounds_with_inf(arg2)

    integer_domains = ComponentSet([Binary, Integers])

    integer_exponent = False
    # if the exponent is an integer, there should not be any evaluation errors
    if isinstance(arg2, VarData) and arg2.domain in integer_domains:
        # The exponent is an integer variable
        # check if the base can be zero
        integer_exponent = True
    if lb2 == ub2 and lb2 == round(lb2):
        # The exponent is fixed to an integer
        integer_exponent = True
    repn = generate_standard_repn(arg2, quadratic=True)
    if (
        repn.nonlinear_expr is None
        and repn.constant == round(repn.constant)
        and all(i.domain in integer_domains for i in repn.linear_vars)
        and all(i[0].domain in integer_domains for i in repn.quadratic_vars)
        and all(i[1].domain in integer_domains for i in repn.quadratic_vars)
        and all(i == round(i) for i in repn.linear_coefs)
        and all(i == round(i) for i in repn.quadratic_coefs)
    ):
        # The exponent is a linear or quadratic expression containing
        # only integer variables with integer coefficients
        integer_exponent = True

    if integer_exponent and (
        (lb1 > 0 or ub1 < 0)
        or (not config.warn_for_evaluation_error_at_bounds and (lb1 >= 0 or ub1 <= 0))
    ):
        # life is good; the exponent is an integer and the base is nonzero
        return None
    elif integer_exponent and lb2 >= 0:
        # life is good; the exponent is a nonnegative integer
        return None

    # if the base is positive, there should not be any evaluation errors
    if lb1 > 0 or (not config.warn_for_evaluation_error_at_bounds and lb1 >= 0):
        return None
    if lb1 >= 0 and lb2 >= 0:
        return None

    msg = f"Potential evaluation error in {node}; "
    msg += f"base bounds are ({lb1}, {ub1}); "
    msg += f"exponent bounds are ({lb2}, {ub2})"
    warn_list.append(msg)


def _check_eval_error_log(
    node: NumericExpression, warn_list: List[str], config: ConfigDict
):
    lb, ub = _get_bounds_with_inf(node.args[0])
    if (config.warn_for_evaluation_error_at_bounds and lb <= 0) or lb < 0:
        msg = f"Potential log of a non-positive number in {node}; Argument bounds are ({lb}, {ub})"
        warn_list.append(msg)


def _check_eval_error_tan(
    node: NumericExpression, warn_list: List[str], config: ConfigDict
):
    lb, ub = _get_bounds_with_inf(node)
    if not (isfinite(lb) and isfinite(ub)):
        msg = f"{node} may evaluate to -inf or inf; Argument bounds are {_get_bounds_with_inf(node.args[0])}"
        warn_list.append(msg)


def _check_eval_error_asin(
    node: NumericExpression, warn_list: List[str], config: ConfigDict
):
    lb, ub = _get_bounds_with_inf(node.args[0])
    if lb < -1 or ub > 1:
        msg = f"Potential evaluation of asin outside [-1, 1] in {node}; Argument bounds are ({lb}, {ub})"
        warn_list.append(msg)


def _check_eval_error_acos(
    node: NumericExpression, warn_list: List[str], config: ConfigDict
):
    lb, ub = _get_bounds_with_inf(node.args[0])
    if lb < -1 or ub > 1:
        msg = f"Potential evaluation of acos outside [-1, 1] in {node}; Argument bounds are ({lb}, {ub})"
        warn_list.append(msg)


def _check_eval_error_sqrt(
    node: NumericExpression, warn_list: List[str], config: ConfigDict
):
    lb, ub = _get_bounds_with_inf(node.args[0])
    if lb < 0:
        msg = f"Potential square root of a negative number in {node}; Argument bounds are ({lb}, {ub})"
        warn_list.append(msg)


_unary_eval_err_handler = dict()
_unary_eval_err_handler["log"] = _check_eval_error_log
_unary_eval_err_handler["log10"] = _check_eval_error_log
_unary_eval_err_handler["tan"] = _check_eval_error_tan
_unary_eval_err_handler["asin"] = _check_eval_error_asin
_unary_eval_err_handler["acos"] = _check_eval_error_acos
_unary_eval_err_handler["sqrt"] = _check_eval_error_sqrt


def _check_eval_error_unary(
    node: NumericExpression, warn_list: List[str], config: ConfigDict
):
    if node.getname() in _unary_eval_err_handler:
        _unary_eval_err_handler[node.getname()](node, warn_list, config)


_eval_err_handler = dict()
_eval_err_handler[DivisionExpression] = _check_eval_error_division
_eval_err_handler[NPV_DivisionExpression] = _check_eval_error_division
_eval_err_handler[PowExpression] = _check_eval_error_pow
_eval_err_handler[NPV_PowExpression] = _check_eval_error_pow
_eval_err_handler[UnaryFunctionExpression] = _check_eval_error_unary
_eval_err_handler[NPV_UnaryFunctionExpression] = _check_eval_error_unary


class _EvalErrorWalker(StreamBasedExpressionVisitor):
    def __init__(self, config: ConfigDict):
        super().__init__()
        self._warn_list = list()
        self._config = config

    def exitNode(self, node, data):
        """
        callback to be called as the visitor moves from the leaf
        nodes back to the root node.

        Args:
            node: a pyomo expression node
            data: not used in this walker
        """
        if type(node) in _eval_err_handler:
            _eval_err_handler[type(node)](node, self._warn_list, self._config)
        return self._warn_list


# TODO: Rename and redirect once old DegeneracyHunter is removed
@document_kwargs_from_configdict(DHCONFIG)
class DegeneracyHunter2:
    """
    Degeneracy Hunter is a tool for identifying Irreducible Degenerate Sets (IDS) in
    Pyomo models.

    Original implementation by Alex Dowling.

    Args:

        model: model to be diagnosed. The DegeneracyHunter does not support indexed Blocks.

    """

    def __init__(self, model, **kwargs):
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
        self.config = DHCONFIG(kwargs)

        # Get Jacobian and NLP
        self.jacobian, self.nlp = get_jacobian(
            self._model, equality_constraints_only=True
        )

        # Placeholder for solver - deferring construction lets us unit test more easily
        self.solver = None

        # Create placeholders for results
        self.degenerate_set = {}
        self.irreducible_degenerate_sets = []

    def _get_solver(self):
        if self.solver is None:
            self.solver = SolverFactory(self.config.solver)

            options = self.config.solver_options
            if options is None:
                options = {}

            self.solver.options = options

        return self.solver

    def _prepare_candidates_milp(self):
        """
        Prepare MILP to find candidate equations for consider for IDS

        Args:

            None

        Returns:

            m_fc: Pyomo model to find candidates

        """
        _log.info("Building MILP model.")

        # Create Pyomo model for irreducible degenerate set
        m_dh = ConcreteModel()

        # Create index for constraints
        m_dh.C = Set(initialize=range(self.jacobian.shape[0]))
        m_dh.V = Set(initialize=range(self.jacobian.shape[1]))

        # Specify minimum size for nu to be considered non-zero
        M = self.config.M
        m_small = self.config.m_small

        # Add variables
        m_dh.nu = Var(
            m_dh.C,
            bounds=(-M - m_small, M + m_small),
            initialize=1.0,
        )
        m_dh.y_pos = Var(m_dh.C, domain=Binary)
        m_dh.y_neg = Var(m_dh.C, domain=Binary)
        m_dh.abs_nu = Var(m_dh.C, bounds=(0, M + m_small))

        m_dh.pos_xor_neg = Constraint(m_dh.C)

        # Constraint to enforce set is degenerate
        if issparse(self.jacobian):
            J = self.jacobian.tocsc()

            def eq_degenerate(m_dh, v):
                if np.sum(np.abs(J[:, v])) > self.config.trivial_constraint_tolerance:
                    # Find the columns with non-zero entries
                    C_ = find(J[:, v])[0]
                    return sum(J[c, v] * m_dh.nu[c] for c in C_) == 0
                else:
                    # This variable does not appear in any constraint
                    return Constraint.Skip

        else:
            J = self.jacobian

            def eq_degenerate(m_dh, v):
                if np.sum(np.abs(J[:, v])) > self.config.trivial_constraint_tolerance:
                    return sum(J[c, v] * m_dh.nu[c] for c in m_dh.C) == 0
                else:
                    # This variable does not appear in any constraint
                    return Constraint.Skip

        m_dh.degenerate = Constraint(m_dh.V, rule=eq_degenerate)

        # When y_pos = 1, nu >= m_small
        # When y_pos = 0, nu >= - m_small
        def eq_pos_lower(b, c):
            return b.nu[c] >= -m_small + 2 * m_small * b.y_pos[c]

        m_dh.pos_lower = Constraint(m_dh.C, rule=eq_pos_lower)

        # When y_pos = 1, nu <= M + m_small
        # When y_pos = 0, nu <= m_small
        def eq_pos_upper(b, c):
            return b.nu[c] <= M * b.y_pos[c] + m_small

        m_dh.pos_upper = Constraint(m_dh.C, rule=eq_pos_upper)

        # When y_neg = 1, nu <= -m_small
        # When y_neg = 0, nu <= m_small
        def eq_neg_upper(b, c):
            return b.nu[c] <= m_small - 2 * m_small * b.y_neg[c]

        m_dh.neg_upper = Constraint(m_dh.C, rule=eq_neg_upper)

        # When y_neg = 1, nu >= -M - m_small
        # When y_neg = 0, nu >= - m_small
        def eq_neg_lower(b, c):
            return b.nu[c] >= -M * b.y_neg[c] - m_small

        m_dh.neg_lower = Constraint(m_dh.C, rule=eq_neg_lower)

        # Absolute value
        def eq_abs_lower(b, c):
            return -b.abs_nu[c] <= b.nu[c]

        m_dh.abs_lower = Constraint(m_dh.C, rule=eq_abs_lower)

        def eq_abs_upper(b, c):
            return b.nu[c] <= b.abs_nu[c]

        m_dh.abs_upper = Constraint(m_dh.C, rule=eq_abs_upper)

        # At least one constraint must be in the degenerate set
        m_dh.degenerate_set_nonempty = Constraint(
            expr=sum(m_dh.y_pos[c] + m_dh.y_neg[c] for c in m_dh.C) >= 1
        )

        # Minimize the L1-norm of nu
        m_dh.obj = Objective(expr=sum(m_dh.abs_nu[c] for c in m_dh.C))

        self.candidates_milp = m_dh

    def _identify_candidates(self):
        eq_con_list = self.nlp.get_pyomo_equality_constraints()

        for i in self.candidates_milp.C:
            # Check if constraint is included
            if self.candidates_milp.abs_nu[i]() > self.config.m_small * MMULT:
                # If it is, save the value of nu
                if eq_con_list is None:
                    name = i
                else:
                    name = eq_con_list[i]
                self.degenerate_set[name] = self.candidates_milp.nu[i]()

    def _solve_candidates_milp(self, tee: bool = False):
        """Solve MILP to generate set of candidate equations

        Arguments:

            tee: print solver output (default = False)

        """
        _log.info("Solving Candidates MILP model.")

        results = self._get_solver().solve(self.candidates_milp, tee=tee)

        self.degenerate_set = {}

        if check_optimal_termination(results):
            # We found a degenerate set
            self._identify_candidates()
        else:
            _log.debug(
                "Solver did not return an optimal termination condition for "
                "Candidates MILP. This probably indicates the system is full rank."
            )

    def _prepare_ids_milp(self):
        """
        Prepare MILP to compute the irreducible degenerate set

        """
        _log.info("Building MILP model to compute irreducible degenerate set.")

        n_eq = self.jacobian.shape[0]
        n_var = self.jacobian.shape[1]

        # Create Pyomo model for irreducible degenerate set
        m_dh = ConcreteModel()

        # Create index for constraints
        m_dh.C = Set(initialize=range(n_eq))
        m_dh.V = Set(initialize=range(n_var))

        # Specify minimum size for nu to be considered non-zero
        M = self.config.M

        # Add variables
        m_dh.nu = Var(m_dh.C, bounds=(-M, M), initialize=1.0)
        m_dh.y = Var(m_dh.C, domain=Binary)

        # Constraint to enforce set is degenerate
        if issparse(self.jacobian):
            J = self.jacobian.tocsc()

            def eq_degenerate(m_dh, v):
                # Find the columns with non-zero entries
                C = find(J[:, v])[0]
                if len(C) == 0:
                    # Catch for edge-case of trivial constraint 0==0
                    return Constraint.Skip
                return sum(J[c, v] * m_dh.nu[c] for c in C) == 0

        else:
            J = self.jacobian

            def eq_degenerate(m_dh, v):
                return sum(J[c, v] * m_dh.nu[c] for c in m_dh.C) == 0

        m_dh.degenerate = Constraint(m_dh.V, rule=eq_degenerate)

        def eq_lower(m_dh, c):
            return -M * m_dh.y[c] <= m_dh.nu[c]

        m_dh.lower = Constraint(m_dh.C, rule=eq_lower)

        def eq_upper(m_dh, c):
            return m_dh.nu[c] <= M * m_dh.y[c]

        m_dh.upper = Constraint(m_dh.C, rule=eq_upper)

        m_dh.obj = Objective(expr=sum(m_dh.y[c] for c in m_dh.C))

        self.ids_milp = m_dh

    def _get_ids(self):
        # Create empty dictionary
        ids_ = {}

        eq_con_list = self.nlp.get_pyomo_equality_constraints()

        for i in self.ids_milp.C:
            # Check if constraint is included
            if self.ids_milp.y[i]() > YTOL:
                # If it is, save the value of nu
                ids_[eq_con_list[i]] = self.ids_milp.nu[i]()

        return ids_

    def _solve_ids_milp(self, cons: Constraint, tee: bool = False):
        """Solve MILP to check if equation 'cons' is a significant component
        in an irreducible degenerate set

        Args:
            cons: constraint to consider
            tee: Boolean, print solver output (default = False)

        Returns:
            ids: dictionary containing the IDS

        """
        _log.info(f"Solving IDS MILP for constraint {cons.name}.")
        eq_con_list = self.nlp.get_pyomo_equality_constraints()
        cons_idx = eq_con_list.index(cons)

        # Fix weight on candidate equation
        self.ids_milp.nu[cons_idx].fix(1.0)

        # Solve MILP
        results = self._get_solver().solve(self.ids_milp, tee=tee)

        self.ids_milp.nu[cons_idx].unfix()

        if check_optimal_termination(results):
            # We found an irreducible degenerate set
            return self._get_ids()
        else:
            raise ValueError(
                f"Solver did not return an optimal termination condition for "
                f"IDS MILP with constraint {cons.name}."
            )

    def find_irreducible_degenerate_sets(self, tee=False):
        """
        Compute irreducible degenerate sets

        Args:
            tee: Print solver output logs to screen (default=False)

        """

        # Solve to find candidate equations
        _log.info("Searching for Candidate Equations")
        self._prepare_candidates_milp()
        self._solve_candidates_milp(tee=tee)

        # Find irreducible degenerate sets
        # Check if degenerate_set is not empty
        if self.degenerate_set:

            _log.info("Searching for Irreducible Degenerate Sets")
            self._prepare_ids_milp()

            # Loop over candidate equations
            count = 1
            for k in self.degenerate_set:
                print(f"Solving MILP {count} of {len(self.degenerate_set)}.")
                _log.info_high(f"Solving MILP {count} of {len(self.degenerate_set)}.")

                # Check if equation is a major element of an IDS
                ids_ = self._solve_ids_milp(cons=k, tee=tee)

                if ids_ is not None:
                    self.irreducible_degenerate_sets.append(ids_)

                count += 1
        else:
            _log.debug("No candidate equations found")

    def report_irreducible_degenerate_sets(self, stream=None, tee: bool = False):
        """
        Print a report of all the Irreducible Degenerate Sets (IDS) identified in
        model.

        Args:
            stream: I/O object to write report to (default = stdout)
            tee: whether to write solver output logs to screen

        Returns:
            None

        """
        if stream is None:
            stream = sys.stdout

        self.find_irreducible_degenerate_sets(tee=tee)

        stream.write("=" * MAX_STR_LENGTH + "\n")
        stream.write("Irreducible Degenerate Sets\n")

        if self.irreducible_degenerate_sets:
            for i, s in enumerate(self.irreducible_degenerate_sets):
                stream.write(f"\n{TAB}Irreducible Degenerate Set {i}")
                stream.write(f"\n{TAB*2}nu{TAB}Constraint Name")
                for k, v in s.items():
                    value_string = f"{v:.1f}"
                    sep = (2 + len(TAB) - len(value_string)) * " "
                    stream.write(f"\n{TAB*2}{value_string}{sep}{k.name}")
                stream.write("\n")
        else:
            stream.write(
                f"\n{TAB}No candidate equations. The Jacobian is likely full rank.\n"
            )

        stream.write("\n" + "=" * MAX_STR_LENGTH + "\n")


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
            self._eq_con_list = self.nlp.get_pyomo_equality_constraints()
            self._var_list = self.nlp.get_pyomo_variables()

            self.candidate_eqns = None

        elif type(block_or_jac) is np.array:  # pylint: disable=unidiomatic-typecheck
            raise NotImplementedError(
                "Degeneracy Hunter currently only supports analyzing a Pyomo model"
            )

            # # TODO: Need to refactor, document, and test support for Jacobian
            # self.jac_eq = block_or_jac
            # self._eq_con_list = None

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
            print(str(i) + ": " + self._var_list[i].name)
        print("")
        print("Row:    Constraint")
        for i in np.where(abs(self.u[:, n_calc - 1]) > tol)[0]:
            print(str(i) + ": " + self._eq_con_list[i].name)

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
            self.candidates_milp, self.solver, self._eq_con_list, tee
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
                    self.dh_milp, self.solver, c, self._eq_con_list, tee
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


def psweep_runner_validator(val):
    """Domain validator for Parameter Sweep runners

    Args:
        val : value to be checked

    Returns:
        TypeError if val is not a valid callback
    """
    if issubclass(val, ParameterSweepBase):
        return val

    raise ValueError("Workflow runner must be a subclass of ParameterSweepBase.")


CACONFIG = ConfigDict()
CACONFIG.declare(
    "input_specification",
    ConfigValue(
        domain=is_psweepspec,
        doc="ParameterSweepSpecification object defining inputs to be sampled",
    ),
)
CACONFIG.declare(
    "workflow_runner",
    ConfigValue(
        default=SequentialSweepRunner,
        domain=psweep_runner_validator,
        doc="Parameter sweep workflow runner",
    ),
)
CACONFIG.declare(
    "solver_options",
    ConfigValue(
        domain=None,
        description="Options to pass to IPOPT.",
    ),
)
CACONFIG.declare(
    "halt_on_error",
    ConfigValue(
        default=False,
        domain=bool,
        doc="Whether to halt execution of parameter sweep on encountering a solver error (default=False).",
    ),
)


class IpoptConvergenceAnalysis:
    """
    Tool to perform a parameter sweep of model checking for numerical issues and
    convergence characteristics. Users may specify an IDAES ParameterSweep class to
    perform the sweep (default is SequentialSweepRunner).
    """

    CONFIG = CACONFIG()

    def __init__(self, model, solver_obj=None, **kwargs):
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
        self.config = self.CONFIG(kwargs)

        self._model = model

        if solver_obj is None:
            solver_obj = SolverFactory("ipopt")
        self._solver_obj = solver_obj

        if self.config.solver_options is not None:
            solver_obj.options = self.config.solver_options

        self._psweep = self.config.workflow_runner(
            input_specification=self.config.input_specification,
            build_model=self._build_model,
            rebuild_model=True,
            run_model=self._run_model,
            build_outputs=self._build_outputs,
            halt_on_error=self.config.halt_on_error,
            handle_solver_error=self._recourse,
            solver=solver_obj,
        )

    @property
    def results(self):
        """
        Returns the results of the IpoptConvergenceAnalysis run
        """
        return self._psweep.results

    @property
    def samples(self):
        """
        Returns the set of input samples for convergence analysis (pandas DataFrame)
        """
        return self._psweep.get_input_samples()

    def run_convergence_analysis(self):
        """
        Execute convergence analysis sweep by calling execute_parameter_sweep
        method in specified runner.

        Returns:
            dict of results from parameter sweep
        """
        return self._psweep.execute_parameter_sweep()

    def run_convergence_analysis_from_dict(self, input_dict: dict):
        """
        Execute convergence analysis sweep using specification defined in dict.

        Args:
            input_dict: dict to load specification from

        Returns:
            dict of results from parameter sweep
        """
        self.from_dict(input_dict)
        return self.run_convergence_analysis()

    def run_convergence_analysis_from_file(self, filename: str):
        """
        Execute convergence analysis sweep using specification defined in json file.

        Args:
            filename: name of file to load specification from as string

        Returns:
            dict of results from parameter sweep
        """
        self.from_json_file(filename)
        return self.run_convergence_analysis()

    def compare_convergence_to_baseline(
        self, filename: str, rel_tol: float = 0.1, abs_tol: float = 1
    ):
        """
        Run convergence analysis and compare results to those defined in baseline file.

        Args:
            filename: name of baseline file to load specification from as string
            rel_tol: relative tolerance to use for comparing number of iterations
            abs_tol: absolute tolerance to use for comparing number of iterations

        Returns:
            dict containing lists of differences between convergence analysis run and baseline
        """
        with open(filename, "r") as f:
            # Load file manually so we have saved results
            jdict = json.load(f)
        f.close()

        # Run convergence analysis from dict
        self.run_convergence_analysis_from_dict(jdict)

        # Compare results
        return self._compare_results_to_dict(jdict, rel_tol=rel_tol, abs_tol=abs_tol)

    def assert_baseline_comparison(
        self, filename: str, rel_tol: float = 0.1, abs_tol: float = 1
    ):
        """
        Run convergence analysis and assert no differences in results to those defined
        in baseline file.

        Args:
            filename: name of baseline file to load specification from as string
            rel_tol: relative tolerance to use for comparing number of iterations
            abs_tol: absolute tolerance to use for comparing number of iterations

        Raises:
            AssertionError if results of convergence analysis do not match baseline
        """
        diffs = self.compare_convergence_to_baseline(
            filename, rel_tol=rel_tol, abs_tol=abs_tol
        )

        if any(len(v) != 0 for v in diffs.values()):
            raise AssertionError("Convergence analysis does not match baseline")

    def report_convergence_summary(self, stream=None):
        """
        Reports a brief summary of the model convergence run.

        Args:
            stream: Optional output stream to print results to.

        Returns:
            None

        """
        if stream is None:
            stream = sys.stdout

        successes = 0
        failures = 0
        runs_w_restoration = 0
        runs_w_regulariztion = 0
        runs_w_num_iss = 0

        for v in self.results.values():
            # Check for differences
            if v["success"]:
                successes += 1
            else:
                failures += 1

            if v["results"]["iters_in_restoration"] > 0:
                runs_w_restoration += 1
            if v["results"]["iters_w_regularization"] > 0:
                runs_w_regulariztion += 1
            if v["results"]["numerical_issues"] > 0:
                runs_w_num_iss += 1

        stream.write(
            f"Successes: {successes}, Failures {failures} ({100*successes/(successes+failures)}%)\n"
        )
        stream.write(f"Runs with Restoration: {runs_w_restoration}\n")
        stream.write(f"Runs with Regularization: {runs_w_regulariztion}\n")
        stream.write(f"Runs with Numerical Issues: {runs_w_num_iss}\n")

    def to_dict(self):
        """
        Serialize specification and current results to dict form

        Returns:
            dict
        """
        return self._psweep.to_dict()

    def from_dict(self, input_dict):
        """
        Load specification and results from dict.

        Args:
            input_dict: dict to load from

        Returns:
            None
        """
        return self._psweep.from_dict(input_dict)

    def to_json_file(self, filename):
        """
        Write specification and results to json file.

        Args:
            filename: name of file to write to as string

        Returns:
            None
        """
        return self._psweep.to_json_file(filename)

    def from_json_file(self, filename):
        """
        Load specification and results from json file.

        Args:
            filename: name of file to load from as string

        Returns:
            None
        """
        return self._psweep.from_json_file(filename)

    def _build_model(self):
        # Create new instance of model by cloning
        return self._model.clone()

    def _run_model(self, model, solver):
        # Run model using IPOPT and collect stats
        (
            status,
            iters,
            iters_in_restoration,
            iters_w_regularization,
            time,
        ) = self._run_ipopt_with_stats(model, solver)

        run_stats = [
            iters,
            iters_in_restoration,
            iters_w_regularization,
            time,
        ]

        success = check_optimal_termination(status)

        return success, run_stats

    @staticmethod
    def _build_outputs(model, run_stats):
        # Run model diagnostics numerical checks
        dt = DiagnosticsToolbox(model=model)

        warnings = False
        try:
            dt.assert_no_numerical_warnings()
        except AssertionError:
            warnings = True

        # Compile Results
        return {
            "iters": run_stats[0],
            "iters_in_restoration": run_stats[1],
            "iters_w_regularization": run_stats[2],
            "time": run_stats[3],
            "numerical_issues": warnings,
        }

    @staticmethod
    def _recourse(model):
        # Return a default dict indicating no results
        return {
            "iters": -1,
            "iters_in_restoration": -1,
            "iters_w_regularization": -1,
            "time": -1,
            "numerical_issues": -1,
        }

    @staticmethod
    def _parse_ipopt_output(ipopt_file):
        # PArse IPOPT logs and return key metrics
        # ToDO: Check for final iteration with regularization or restoration

        iters = 0
        iters_in_restoration = 0
        iters_w_regularization = 0
        time = 0
        # parse the output file to get the iteration count, solver times, etc.
        with open(ipopt_file, "r") as f:
            parseline = False
            for line in f:
                if line.startswith("iter"):
                    # This marks the start of the iteration logging, set parseline True
                    parseline = True
                elif line.startswith("Number of Iterations....:"):
                    # Marks end of iteration logging, set parseline False
                    parseline = False
                    tokens = line.split()
                    iters = int(tokens[3])
                elif parseline:
                    # Line contains details of an iteration, look for restoration or regularization
                    tokens = line.split()
                    try:
                        if not tokens[6] == "-":
                            # Iteration with regularization
                            iters_w_regularization += 1
                        if tokens[0].endswith("r"):
                            # Iteration in restoration
                            iters_in_restoration += 1
                    except IndexError:
                        # Blank line at end of iteration list, so assume we hit this
                        pass
                elif line.startswith(
                    "Total CPU secs in IPOPT (w/o function evaluations)   ="
                ):
                    tokens = line.split()
                    time += float(tokens[9])
                elif line.startswith(
                    "Total CPU secs in NLP function evaluations           ="
                ):
                    tokens = line.split()
                    time += float(tokens[8])

        return iters, iters_in_restoration, iters_w_regularization, time

    def _run_ipopt_with_stats(self, model, solver, max_iter=500, max_cpu_time=120):
        # Solve model using provided solver (assumed to be IPOPT) and parse logs
        # ToDo: Check that the "solver" is, in fact, IPOPT

        TempfileManager.push()
        tempfile = TempfileManager.create_tempfile(suffix="ipopt_out", text=True)
        opts = {
            "output_file": tempfile,
            "max_iter": max_iter,
            "max_cpu_time": max_cpu_time,
        }

        status_obj = solver.solve(model, options=opts, tee=True)

        (
            iters,
            iters_in_restoration,
            iters_w_regularization,
            time,
        ) = self._parse_ipopt_output(tempfile)

        TempfileManager.pop(remove=True)
        return status_obj, iters, iters_in_restoration, iters_w_regularization, time

    def _compare_results_to_dict(
        self,
        compare_dict: dict,
        rel_tol: float = 0.1,
        abs_tol: float = 1,
    ):
        # Compare results
        success_diff = []
        iters_diff = []
        restore_diff = []
        reg_diff = []
        num_iss_diff = []

        for k, v in self.results.items():
            # Get sample result from compare dict
            try:
                comp = compare_dict["results"][k]
            except KeyError:
                # Reading from json often converts ints to strings
                # Check to see if the index as a string works
                comp = compare_dict["results"][str(k)]

            # Check for differences
            if v["success"] != comp["success"]:
                success_diff.append(k)
            if not isclose(
                v["results"]["iters"],
                comp["results"]["iters"],
                rel_tol=rel_tol,
                abs_tol=abs_tol,
            ):
                iters_diff.append(k)
            if not isclose(
                v["results"]["iters_in_restoration"],
                comp["results"]["iters_in_restoration"],
                rel_tol=rel_tol,
                abs_tol=abs_tol,
            ):
                restore_diff.append(k)
            if not isclose(
                v["results"]["iters_w_regularization"],
                comp["results"]["iters_w_regularization"],
                rel_tol=rel_tol,
                abs_tol=abs_tol,
            ):
                reg_diff.append(k)
            if v["results"]["numerical_issues"] != comp["results"]["numerical_issues"]:
                num_iss_diff.append(k)

        return {
            "success": success_diff,
            "iters": iters_diff,
            "iters_in_restoration": restore_diff,
            "iters_w_regularization": reg_diff,
            "numerical_issues": num_iss_diff,
        }


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


def compute_ill_conditioning_certificate(
    model,
    target_feasibility_tol: float = 1e-06,
    ratio_cutoff: float = 1e-04,
    direction: str = "row",
):
    """
    Finds constraints (rows) or variables (columns) in the model Jacobian that
    may be contributing to ill conditioning.

    This method is based on work published in:

    Klotz, E., Identification, Assessment, and Correction of Ill-Conditioning and
    Numerical Instability in Linear and Integer Programs, Informs 2014, pgs. 54-108
    https://pubsonline.informs.org/doi/epdf/10.1287/educ.2014.0130

    Args:
        model: model to be analysed
        target_feasibility_tol: target tolerance for solving ill conditioning problem
        ratio_cutoff: cut-off for reporting ill conditioning
        direction: 'row' (default, constraints) or 'column' (variables)

    Returns:
        list of strings reporting ill-conditioned variables/constraints and their
        associated y values
    """
    _log.warning(
        "Ill conditioning checks are a beta capability. Please be aware that "
        "the name, location, and API for this may change in future releases."
    )
    # Thanks to B. Knueven for this implementation

    if direction not in ["row", "column"]:
        raise ValueError(
            f"Unrecognised value for direction ({direction}). "
            "Must be 'row' or 'column'."
        )

    jac, nlp = get_jacobian(model)

    inverse_target_kappa = 1e-16 / target_feasibility_tol

    # Set up the components we will analyze, either row or column
    if direction == "row":
        components = nlp.get_pyomo_constraints()
        components_set = RangeSet(0, len(components) - 1)
        results_set = RangeSet(0, nlp.n_primals() - 1)
        jac = jac.transpose().tocsr()
    else:  # direction == "column"
        components = nlp.get_pyomo_variables()
        components_set = RangeSet(0, len(components) - 1)
        results_set = RangeSet(0, nlp.n_constraints() - 1)
        jac = jac.tocsr()

    # Build test problem
    inf_prob = ConcreteModel()

    inf_prob.y_pos = Var(components_set, bounds=(0, None))
    inf_prob.y_neg = Var(components_set, bounds=(0, None))
    inf_prob.y = Expression(components_set, rule=lambda m, i: m.y_pos[i] - m.y_neg[i])

    inf_prob.res_pos = Var(results_set, bounds=(0, None))
    inf_prob.res_neg = Var(results_set, bounds=(0, None))
    inf_prob.res = Expression(
        results_set, rule=lambda m, i: m.res_pos[i] - m.res_neg[i]
    )

    def b_rule(b, i):
        lhs = 0.0

        row = jac.getrow(i)
        for j, val in zip(row.indices, row.data):
            lhs += val * b.y[j]

        return lhs == b.res[i]

    inf_prob.by = Constraint(results_set, rule=b_rule)

    # Normalization of y
    inf_prob.normalize = Constraint(
        expr=1 == sum(inf_prob.y_pos.values()) - sum(inf_prob.y_neg.values())
    )

    inf_prob.y_norm = Var()
    inf_prob.y_norm_constr = Constraint(
        expr=inf_prob.y_norm
        == sum(inf_prob.y_pos.values()) + sum(inf_prob.y_neg.values())
    )

    inf_prob.res_norm = Var()
    inf_prob.res_norm_constr = Constraint(
        expr=inf_prob.res_norm
        == sum(inf_prob.res_pos.values()) + sum(inf_prob.res_neg.values())
    )

    # Objective -- minimize residual
    inf_prob.min_res = Objective(expr=inf_prob.res_norm)

    solver = SolverFactory("cbc")  # TODO: Consider making this an option

    # tighten tolerances  # TODO: If solver is an option, need to allow user options
    solver.options["primalT"] = target_feasibility_tol * 1e-1
    solver.options["dualT"] = target_feasibility_tol * 1e-1

    results = solver.solve(inf_prob, tee=False)
    if not check_optimal_termination(results):
        # TODO: maybe we should tighten tolerances first?
        raise RuntimeError("Ill conditioning diagnostic problem infeasible")

    result_norm = inf_prob.res_norm.value
    if result_norm < 0.0:
        # TODO: try again with tighter tolerances?
        raise RuntimeError(
            "Ill conditioning diagnostic problem has numerically troublesome solution"
        )
    if result_norm >= inverse_target_kappa:
        return []

    # find an equivalent solution which minimizes y_norm
    inf_prob.min_res.deactivate()
    inf_prob.res_norm.fix()

    inf_prob.min_y = Objective(expr=inf_prob.y_norm)

    # if this problem is numerically infeasible, we can still report something to the user
    results = solver.solve(inf_prob, tee=False, load_solutions=False)
    if check_optimal_termination(results):
        inf_prob.solutions.load_from(results)

    ill_cond = []
    slist = sorted(
        inf_prob.y, key=lambda dict_key: abs(value(inf_prob.y[dict_key])), reverse=True
    )
    cutoff = None
    for i in slist:
        if cutoff is None:
            cutoff = abs(value(inf_prob.y[i])) * ratio_cutoff
        val = value(inf_prob.y[i])
        if abs(val) < cutoff:
            break
        ill_cond.append((components[i], val))

    return ill_cond


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


class ConstraintTermAnalysisVisitor(EXPR.StreamBasedExpressionVisitor):
    """
    Expression walker for checking Constraints for problematic terms.

    This walker will walk the expression and look for summation terms
    with mismatched magnitudes or potential cancellations.

    Args:
        term_mismatch_tolerance: tolerance to use when determining mismatched
            terms
        term_cancellation_tolerance: tolerance to use when identifying
            possible cancellation of terms
        term_zero_tolerance: tolerance for considering terms equal to zero
        max_canceling_terms: maximum number of terms to consider when looking
            for canceling combinations (None = consider all possible combinations)
        max_cancellations_per_node: maximum number of cancellations to collect
            for a single node. Collection will terminate when this many cancellations
            have been identified (None = collect all cancellations)

    Returns:
        list of values for top-level summation terms
        list of terms with mismatched magnitudes
        list of terms with potential cancellations
        bool indicating whether expression is a constant
    """

    def __init__(
        self,
        term_mismatch_tolerance: float = 1e6,
        term_cancellation_tolerance: float = 1e-4,
        term_zero_tolerance: float = 1e-10,
        max_canceling_terms: int = 4,
        max_cancellations_per_node: int = 5,
    ):
        super().__init__()

        # Tolerance attributes
        self._log_mm_tol = log10(term_mismatch_tolerance)
        self._sum_tol = term_cancellation_tolerance
        self._zero_tolerance = term_zero_tolerance
        self._max_canceling_terms = max_canceling_terms
        self._max_cancellations_per_node = max_cancellations_per_node

        # Placeholders for collecting results
        self.canceling_terms = ComponentMap()
        self.mismatched_terms = ComponentMap()

        # Flag for if cancellation collection hit limit
        self._cancellation_tripped = False

    def _get_value_for_sum_subexpression(self, child_data):
        # child_data is a tuple, with the 0-th element being the node values
        if isinstance(child_data[0][0], str):
            # Values may be a list containing a string in some cases (e.g. external functions)
            # Return the string in this case
            return child_data[0][0]
        return sum(i for i in child_data[0])

    def _generate_combinations(self, inputs, equality=False):
        # We want to test all combinations of terms for cancellation

        # The number of combinations we check depends on whether this is an (in)equality
        # expression or a sum node deeper in the expression tree.
        # We expect (in)equalities to generally sum to 0 (0 == expr) thus we want to
        # check if any subset of the sum terms sum to zero (i.e. are any terms unnecessary).
        # For other sum nodes, we need to check for any combination of terms.

        # Maximum number of terms to include in combinations
        max_comb = len(inputs)
        if equality:
            # Subtract 1 if (in)equality node
            max_comb += -1
        # We also have a limit on the maximum number of terms to consider
        if self._max_canceling_terms is not None:
            max_comb = min(max_comb, self._max_canceling_terms)

        # Single terms cannot cancel, thus we want all combinations of length 2 to max terms
        # Note the need for +1 due to way range works
        combo_range = range(2, max_comb + 1)

        # Collect combinations of terms in an iterator
        # In order to identify the terms in each cancellation, we will pair each value
        # with its index in the input set as a tuple using enumerate
        for i in chain.from_iterable(
            combinations(enumerate(inputs), r) for r in combo_range
        ):
            # Yield each combination in the set
            yield i

    def _check_sum_cancellations(self, values_list, equality=False):
        # First, strip any terms with value 0 as they do not contribute to cancellation
        # We do this to keep the number of possible combinations as small as possible
        stripped = [i for i in values_list if abs(i) >= self._zero_tolerance]

        cancellations = []

        if len(stripped) == 0:
            # If the stripped list is empty, there are no non-zero terms
            # We can stop here and return False as there are no possible cancellations
            return cancellations

        # For scaling of tolerance, we want to compare to the largest absolute value of
        # the input values
        max_value = abs(max(stripped, key=abs))

        for i in self._generate_combinations(stripped, equality):
            # Generate combinations will return a list of combinations of input terms
            # each element of the list will be a 2-tuple representing a term in the
            # input list with the first value being the position in the input set and
            # the second being the value

            # Check if the sum of values in the combination are below tolerance
            if abs(sum(j[1] for j in i)) <= self._sum_tol * max_value:
                # If so, record combination as canceling
                cancellations.append(i)

            # Terminate loop if we have reached the max cancellations to collect
            if len(cancellations) >= self._max_cancellations_per_node:
                self._cancellation_tripped = True
                break

        return cancellations

    def _perform_checks(self, node, child_data):
        # Perform checks for problematic expressions
        # First, need to check to see if any child data is a list
        # This indicates a sum expression
        const = True

        for d in child_data:
            # We will check for canceling terms here, rather than the sum itself, to handle special cases
            # We want to look for cases where a sum term results in a value much smaller
            # than the terms of the sum
            # Each element of child_data is a tuple where the 0-th element is the node values
            if isinstance(d[0][0], str):
                # Values may be a list containing a string in some cases (e.g. external functions)
                # Skip if this is the case
                pass
            else:
                for c in self._check_sum_cancellations(d[0]):
                    if node in self.canceling_terms.keys():
                        self.canceling_terms[node].append(c)
                    else:
                        self.canceling_terms[node] = [c]

            # Expression is not constant if any child is not constant
            # Element 1 is a bool indicating if the child node is constant
            if not d[1]:
                const = False

        # Return any problematic terms found
        return const

    def _check_base_type(self, node):
        if isinstance(node, VarData):
            const = node.fixed
        else:
            const = True
        return [value(node)], const, False

    def _check_equality_expression(self, node, child_data):
        # (In)equality expressions are a special case of sum expressions
        # child_data has two elements; 0 is the LHS of the (in)equality and 1 is the RHS
        # Each of these then contains three elements; 0 is a list of values for the sum components,
        # 1 is a bool indicating if the child node term is constant, and 2 is a bool indicating if
        # the child ode is a sum expression.

        # First, to check for cancellation we need to negate one side of the expression
        # mdata will contain the new set of child_data with negated values
        mdata = []
        # child_data[0][0] has the values of the LHS of the (in)equality, and we will negate these
        vals = []
        for j in child_data[0][0]:
            vals.append(-j)
        # Append the negated values along with whether the node is constant to mdata
        mdata.append((vals, child_data[0][1]))
        # child_data[1] is the RHS, so we take this as it appears and append to mdata
        mdata.append(child_data[1])

        # Next, call the method to check the sum expression
        vals, const, _ = self._check_sum_expression(node, mdata)

        # Next, we need to check for canceling terms.
        # child_data[x][2] indicates if a node is a sum expression or not
        if not child_data[0][2] and not child_data[1][2]:
            # If both sides are not sum expressions, we have the form a == b
            # Simple lLinking constraints do not need to be checked
            pass
        # Next, we can ignore any term that has already been flagged as mismatched
        elif node in self.mismatched_terms.keys():
            pass
        # We can also ignore any case where one side of the (in)equality is constant
        # I.e. if either child_node[x][1] is True
        elif any(d[1] for d in mdata):
            pass
        else:
            # Check for cancellation
            # First, collect terms from both sides
            # Note: outer loop comes first in list comprehension
            t = [i for d in mdata for i in d[0]]

            # Then check for cancellations
            for c in self._check_sum_cancellations(t, equality=True):
                if node in self.canceling_terms.keys():
                    self.canceling_terms[node].append(c)
                else:
                    self.canceling_terms[node] = [c]

        return vals, const, False

    def _check_product_expr(self, node, child_data):
        # We expect child_data to be a tuple of len 2 (a*b)
        # If both or neither a and b are sums (xor), handle like other expressions
        if not child_data[0][2] ^ child_data[1][2]:
            return self._check_general_expr(node, child_data)
        else:
            # Here we have the case of a sum and a multiplier
            # For this case, we will make the result look like a sum with the
            # multiplier applied
            # This is important for scaling of the sum terms, as cancellation should be
            # Checked on the scaled value

            # First, check if both terms are constants - if so we can just get the value of
            # this node and move on.
            if child_data[0][1] and child_data[1][1]:
                return self._check_general_expr(node, child_data)

            # The length of the values (child_data[i][0]) of the multiplier will be 1
            # We can just iterate over all terms in both value lists and multiply
            # each pair
            vals = [i * j for i in child_data[0][0] for j in child_data[1][0]]

            # Return the list of values, not constant, is a sum expression (apparent)
            return vals, False, True

    def _check_div_expr(self, node, child_data):
        # We expect child_data to be a tuple of len 2 (a/b)
        # If the numerator is not a sum, handle like general expression
        # child_data[0] is numerator, child_data[1] is denominator
        if not child_data[0][2]:
            return self._check_general_expr(node, child_data)
        else:
            # If the numerator is a sum, we will treat this as if the division
            # were applied to each term in the numerator separately
            # This is important for scaling of the sum terms, as cancellation should be
            # Checked on the scaled value

            # First, check if the numerator is constant - if so we can just get the value of
            # this node and move on.
            if child_data[0][1]:
                return self._check_general_expr(node, child_data)

            # Next, we need to get the value of the denominator
            denom = self._get_value_for_sum_subexpression(child_data[1])

            try:
                vals = [i / denom for i in child_data[0][0]]
            except ZeroDivisionError:
                raise ZeroDivisionError(
                    f"Error in ConstraintTermAnalysisVisitor: found division with denominator of 0 "
                    f"({str(node)})."
                )

            # Return the list of values, not constant, is a sum expression (apparent)
            return vals, False, True

    def _check_negation_expr(self, node, child_data):
        # Here we need to defer checking of cancellations too due to different ways
        # these can appear in an expression.
        # We will simply negate all values on the child node values (child_data[0][0])
        # and pass on the rest.
        vals = [-i for i in child_data[0][0]]

        # Return the list of values, not constant, is a sum expression (apparent)
        return vals, child_data[0][1], child_data[0][2]

    def _check_general_expr(self, node, child_data):
        const = self._perform_checks(node, child_data)

        try:
            # pylint: disable-next=protected-access
            val = node._apply_operation(
                list(map(self._get_value_for_sum_subexpression, child_data))
            )
        except ValueError:
            raise ValueError(
                f"Error in ConstraintTermAnalysisVisitor: error evaluating {str(node)}."
            )
        except ZeroDivisionError:
            raise ZeroDivisionError(
                f"Error in ConstraintTermAnalysisVisitor: found division with denominator of 0 "
                f"({str(node)})."
            )
        except Exception as err:
            # Catch and re-raise any other error that occurs
            _log.exception(f"Unexpected {err=}, {type(err)=}")
            raise

        return [val], const, False

    def _check_other_expression(self, node, child_data):
        const = self._perform_checks(node, child_data)

        # First, need to get value of input terms, which may be sub-expressions
        input_mag = []
        for i in child_data:
            input_mag.append(self._get_value_for_sum_subexpression(i))

        # Next, create a copy of the function with expected magnitudes as inputs
        newfunc = node.create_node_with_local_data(input_mag)

        # Evaluate new function and return the value along with check results
        return [value(newfunc)], const, False

    def _check_ranged_expression(self, node, child_data):
        # child_data should have 3 elements, LHS, middle term, and RHS
        lhs_vals, lhs_const, _ = self._check_equality_expression(node, child_data[:2])
        rhs_vals, rhs_const, _ = self._check_equality_expression(node, child_data[1:])

        # Constant is a bit vague in terms ranged expressions.
        # We will call the term constant only if all parts are constant
        const = lhs_const and rhs_const

        # For values, we need to avoid double counting the middle term
        # Also for sign convention, we will negate the outer terms
        vals = lhs_vals + [-rhs_vals[1]]

        return vals, const, False

    def _check_sum_expression(self, node, child_data):
        # Sum expressions need special handling
        # For sums, collect all child values into a list
        vals = []
        # We will check for cancellation in this node at the next level
        # Pyomo is generally good at simplifying compound sums
        const = True
        # Collect data from child nodes
        for d in child_data:
            vals.append(self._get_value_for_sum_subexpression(d))

            # Expression is not constant if any child is not constant
            # Element 1 is a bool indicating if node is constant
            if not d[1]:
                const = False

        # Check for mismatched terms
        if len(vals) > 1:
            absvals = [abs(v) for v in vals if abs(v) >= self._zero_tolerance]
            if len(absvals) > 0:
                vl = max(absvals)
                vs = min(absvals)
                if vl != vs:
                    if vs == 0:
                        # This branch is reachable only if the user
                        # set _zero_tolerance to 0 or a negative number
                        diff = log10(vl)
                    else:
                        diff = log10(vl / vs)

                    if diff >= self._log_mm_tol:
                        self.mismatched_terms[node] = (vl, vs)

        return vals, const, True

    node_type_method_map = {
        EXPR.EqualityExpression: _check_equality_expression,
        EXPR.InequalityExpression: _check_equality_expression,
        EXPR.RangedExpression: _check_ranged_expression,
        EXPR.SumExpression: _check_sum_expression,
        EXPR.NPV_SumExpression: _check_sum_expression,
        EXPR.ProductExpression: _check_product_expr,
        EXPR.MonomialTermExpression: _check_product_expr,
        EXPR.NPV_ProductExpression: _check_product_expr,
        EXPR.DivisionExpression: _check_div_expr,
        EXPR.NPV_DivisionExpression: _check_div_expr,
        EXPR.PowExpression: _check_general_expr,
        EXPR.NPV_PowExpression: _check_general_expr,
        EXPR.NegationExpression: _check_negation_expr,
        EXPR.NPV_NegationExpression: _check_negation_expr,
        EXPR.AbsExpression: _check_general_expr,
        EXPR.NPV_AbsExpression: _check_general_expr,
        EXPR.UnaryFunctionExpression: _check_general_expr,
        EXPR.NPV_UnaryFunctionExpression: _check_general_expr,
        EXPR.Expr_ifExpression: _check_other_expression,
        EXPR.ExternalFunctionExpression: _check_other_expression,
        EXPR.NPV_ExternalFunctionExpression: _check_other_expression,
        EXPR.LinearExpression: _check_sum_expression,
    }

    def exitNode(self, node, data):
        """
        Method to call when exiting node to check for potential issues.
        """
        # Return [node values], constant (bool), sum expression (bool)
        # first check if the node is a leaf
        nodetype = type(node)

        if nodetype in native_types:
            return [node], True, False

        node_func = self.node_type_method_map.get(nodetype, None)
        if node_func is not None:
            return node_func(self, node, data)

        if not node.is_expression_type():
            # this is a leaf, but not a native type
            if nodetype is _PyomoUnit:
                # Unit have no value, so return 1 as a placeholder
                return [1], True, False

            # Var or Param
            return self._check_base_type(node)
            # might want to add other common types here

        # not a leaf - check if it is a named expression
        if (
            hasattr(node, "is_named_expression_type")
            and node.is_named_expression_type()
        ):
            return self._check_other_expression(node, data)

        raise TypeError(
            f"An unhandled expression node type: {str(nodetype)} was encountered while "
            f"analyzing constraint terms {str(node)}"
        )

    def walk_expression(self, expr):
        """
        Main method to call to walk an expression and return analysis.

        Args:
            expr - expression to be analyzed

        Returns:
            list of values of top-level additive terms
            ComponentSet containing any mismatched terms
            ComponentSet containing any canceling terms
            Bool indicating whether expression is a constant
        """
        # Create new holders for collected terms
        self.canceling_terms = ComponentMap()
        self.mismatched_terms = ComponentMap()

        # Call parent walk_expression method
        vals, const, _ = super().walk_expression(expr)

        # Return results
        return (
            vals,
            self.mismatched_terms,
            self.canceling_terms,
            const,
            self._cancellation_tripped,
        )
