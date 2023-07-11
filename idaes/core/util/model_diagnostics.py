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
from pyomo.util.check_units import assert_units_consistent
from pyomo.core.base.units_container import UnitsError
from pyomo.contrib.incidence_analysis import IncidenceGraphInterface
from pyomo.core.expr.visitor import identify_variables
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
from pyomo.contrib.pynumero.asl import AmplInterface

import idaes.core.util.scaling as iscale
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
import idaes.logger as idaeslog

_log = idaeslog.getLogger(__name__)


class DiagnosticsToolbox:
    def __init__(self, model: Block):
        if not isinstance(model, Block):
            raise ValueError("model argument must be an instance of a Pyomo Block.")
        self.model = model

        # Create placeholders for data
        self._activated__block_set = ComponentSet()
        self._deactivated_block_set = ComponentSet()
        self._activated_equalities_set = ComponentSet()
        self._deactivated_equalities_set = ComponentSet()
        self._activated_inequalities_set = ComponentSet()
        self._deactivated_inequalities_set = ComponentSet()
        self._activated_objectives_set = ComponentSet()
        self._deactivated_objectives_set = ComponentSet()
        self._variables_fixed_to_zero_set = ComponentSet()
        self._variables_in_activated_constraints_set = ComponentSet()
        self._fixed_variables_in_activated_constraints_set = ComponentSet()
        self._unfixed_variables_in_activated_constraints_set = ComponentSet()
        self._external_fixed_variables_in_activated_constraints_set = ComponentSet()
        self._external_unfixed_variables_in_activated_constraints_set = ComponentSet()
        self._variables_not_in_activated_constraints_set = ComponentSet()
        self._fixed_variables_not_in_activated_constraints_set = ComponentSet()
        self._degrees_of_freedom = None

        self._constraints_with_inconsistent_units = ComponentSet()

        self._var_dm_partition = None
        self._con_dm_partition = None
        self._uc_var = None
        self._uc_con = None
        self._oc_var = None
        self._oc_con = None

    def collect_model_statistics(self):
        # For now, just use model_statistics tools.
        # In future, we may want to look at reworking these to avoid repeatedly
        # iterating over the model.

        # TODO: Variables with bounds

        # Block Statistics
        self._activated_block_set = activated_blocks_set(self.model)
        self._deactivated_block_set = deactivated_blocks_set(self.model)

        # # Constraint statistics
        self._activated_equalities_set = activated_equalities_set(self.model)
        self._deactivated_equalities_set = deactivated_equalities_set(self.model)
        self._activated_inequalities_set = activated_inequalities_set(self.model)
        self._deactivated_inequalities_set = deactivated_inequalities_set(self.model)

        # # Objective statistics
        self._activated_objectives_set = activated_objectives_set(self.model)
        self._deactivated_objectives_set = deactivated_objectives_set(self.model)

        # Variable statistics
        self._variables_in_activated_constraints_set = (
            variables_in_activated_constraints_set(self.model)
        )
        self._fixed_variables_in_activated_constraints_set = ComponentSet()
        self._unfixed_variables_in_activated_constraints_set = ComponentSet()
        self._external_fixed_variables_in_activated_constraints_set = ComponentSet()
        self._external_unfixed_variables_in_activated_constraints_set = ComponentSet()

        def var_in_block(var, block):
            parent = var.parent_block()
            while parent is not None:
                if parent is block:
                    return True
                parent = parent.parent_block()
            return False

        for v in self._variables_in_activated_constraints_set:
            if v.fixed:
                self._fixed_variables_in_activated_constraints_set.add(v)
                if not var_in_block(v, self.model):
                    self._external_fixed_variables_in_activated_constraints_set.add(v)
            else:
                self._unfixed_variables_in_activated_constraints_set.add(v)
                if not var_in_block(v, self.model):
                    self._external_unfixed_variables_in_activated_constraints_set.add(v)

        # Set of variables fixed to 0
        self._variables_fixed_to_zero_set = ComponentSet()
        for v in self.model.component_data_objects(Var, descend_into=True):
            if v.fixed and value(v) == 0:
                self._variables_fixed_to_zero_set.add(v)

        # TODO: Need to see if this includes inequalities or not
        self._variables_not_in_activated_constraints_set = (
            variables_not_in_activated_constraints_set(self.model)
        )
        # Set of Unused fixed variables
        self._fixed_variables_not_in_activated_constraints_set = ComponentSet()
        for v in self._variables_not_in_activated_constraints_set:
            if v.fixed:
                self._fixed_variables_not_in_activated_constraints_set.add(v)

        # Calculate DoF
        self._degrees_of_freedom = degrees_of_freedom(self.model)

    def check_unit_consistency(self):
        # Check unit consistency of each constraint
        self._constraints_with_inconsistent_units = ComponentSet()
        for c in self.model.component_data_objects(Constraint, descend_into=True):
            try:
                assert_units_consistent(c)
            except UnitsError:
                self._constraints_with_inconsistent_units.add(c)

    def check_dulmage_mendelsohn_partition(self):
        self._var_dm_partition = None
        self._con_dm_partition = None
        self._uc_var = None
        self._uc_con = None
        self._oc_var = None
        self._oc_con = None

        igraph = IncidenceGraphInterface(self.model)
        self._var_dm_partition, self._con_dm_partition = igraph.dulmage_mendelsohn()

        # Collect under- and order-constrained sub-system
        self._uc_var = (
            self._var_dm_partition.unmatched + self._var_dm_partition.underconstrained
        )
        self._uc_con = self._con_dm_partition.underconstrained
        self._oc_var = self._var_dm_partition.overconstrained
        self._oc_con = (
            self._con_dm_partition.overconstrained + self._con_dm_partition.unmatched
        )

    # TODO: Block triangularization analysis
    # Number and size of blocks, polynomial degree of 1x1 blocks, simple pivot check of moderate sized sub-blocks?

    def report_structural_issues(self, rerun_analysis=True, stream=stdout):
        # Potential evaluation errors
        # High Index

        # Run checks unless told not to
        if rerun_analysis:
            self.collect_model_statistics()
            self.check_unit_consistency()
            self.check_dulmage_mendelsohn_partition()

        # Collect warnings
        tab = " " * 4
        warnings = []
        if self._degrees_of_freedom != 0:
            dstring = "Degrees"
            if self._degrees_of_freedom == abs(1):
                dstring = "Degree"
            warnings.append(
                f"\n{tab}WARNING: {self._degrees_of_freedom} {dstring} of Freedom"
            )
        if len(self._constraints_with_inconsistent_units) > 0:
            cstring = "Constraints"
            if len(self._constraints_with_inconsistent_units) == 1:
                cstring = "Constraint"
            warnings.append(
                f"\n{tab}WARNING: {len(self._constraints_with_inconsistent_units)} "
                f"{cstring} with inconsistent units"
            )
        if any(
            len(x) > 0 for x in [self._uc_var, self._uc_con, self._oc_var, self._oc_con]
        ):
            warnings.append(
                f"\n{tab}WARNING: Structural singularity found\n"
                f"{tab*2}Under-Constrained Set: {len(self._uc_var)} "
                f"variables, {len(self._uc_con)} constraints\n"
                f"{tab * 2}Over-Constrained Set: {len(self._oc_var)} "
                f"variables, {len(self._oc_con)} constraints"
            )

        # Collect cautions
        cautions = []
        if len(self._variables_fixed_to_zero_set) > 0:
            vstring = "variables"
            if len(self._variables_fixed_to_zero_set) == 1:
                vstring = "variable"
            cautions.append(
                f"\n{tab}Caution: {len(self._variables_fixed_to_zero_set)} "
                f"{vstring} fixed to 0"
            )
        if len(self._variables_not_in_activated_constraints_set) > 0:
            vstring = "variables"
            if len(self._variables_not_in_activated_constraints_set) == 1:
                vstring = "variable"
            cautions.append(
                f"\n{tab}Caution: {len(self._variables_not_in_activated_constraints_set)} "
                f"unused {vstring} "
                f"({len(self._fixed_variables_not_in_activated_constraints_set)} fixed)"
            )

        # Generate report
        max_str_length = 84
        stream.write("\n" + "=" * max_str_length + "\n")
        stream.write("Model Statistics\n\n")
        stream.write(
            f"{tab}Activated Blocks: {len(self._activated_block_set)} "
            f"(Deactivated: {len(self._deactivated_block_set)})\n"
        )
        stream.write(
            f"{tab}Free Variables in Activated Constraints: "
            f"{len(self._unfixed_variables_in_activated_constraints_set)} "
            f"(External: {len(self._external_unfixed_variables_in_activated_constraints_set)})\n"
        )
        stream.write(
            f"{tab}Fixed Variables in Activated Constraints: "
            f"{len(self._fixed_variables_in_activated_constraints_set)} "
            f"(External: {len(self._external_fixed_variables_in_activated_constraints_set)})\n"
        )
        stream.write(
            f"{tab}Activated Equality Constraints: {len(self._activated_equalities_set)} "
            f"(Deactivated: {len(self._deactivated_equalities_set)})\n"
        )
        stream.write(
            f"{tab}Activated Inequality Constraints: {len(self._activated_inequalities_set)} "
            f"(Deactivated: {len(self._deactivated_inequalities_set)})\n"
        )
        stream.write(
            f"{tab}Activated Objectives: {len(self._activated_objectives_set)} "
            f"(Deactivated: {len(self._deactivated_objectives_set)})\n"
        )

        stream.write("\n" + "-" * max_str_length + "\n")
        if len(warnings) > 0:
            stream.write(f"{len(warnings)} WARNINGS\n")
            for w in warnings:
                stream.write(w)
        else:
            stream.write("No warnings found!\n")

        stream.write("\n\n" + "-" * max_str_length + "\n")
        if len(cautions) > 0:
            stream.write(f"{len(cautions)} Cautions\n")
            for c in cautions:
                stream.write(c)
        else:
            stream.write("No cautions found!\n")

        stream.write("\n\n" + "=" * max_str_length + "\n")


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
            self.jac_eq = iscale.get_jacobian(
                self.block, equality_constraints_only=True
            )[0]

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
                "User wanted constraints and variables associated "
                f"with the {n_calc}-th smallest singular value, "
                f"but only {n_sv} small singular values have been "
                "calculated. Run svd_analysis again and specify "
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
