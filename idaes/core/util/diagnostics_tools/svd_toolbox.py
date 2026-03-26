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
This module contains the SVD Toolbox.
"""

__author__ = "Alexander Dowling, Douglas Allan, Andrew Lee, Robby Parker, Ben Knueven"

import sys
from inspect import signature

import numpy as np
from scipy.linalg import svd
from scipy.sparse.linalg import svds

from pyomo.core.base.block import BlockData
from pyomo.core.base.var import VarData
from pyomo.core.base.constraint import ConstraintData
from pyomo.common.config import (
    ConfigDict,
    ConfigValue,
    document_kwargs_from_configdict,
    PositiveInt,
    NonNegativeFloat,
)
from pyomo.common.errors import PyomoException

from idaes.core.util.model_statistics import (
    greybox_block_set,
)
from idaes.core.scaling.util import (
    get_jacobian,
)
from idaes.core.util.diagnostics_tools.writer_utils import (
    write_report_section,
    MAX_STR_LENGTH,
    TAB,
)
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


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
        write_report_section(
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
        write_report_section(
            stream=stream,
            lines_list=vars_in_cons,
            title=f"The following variables are involved in {constraint.name}:",
            header="=",
            footer="=",
        )
