# -*- coding: utf-8 -*-
#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""Degeneracy Hunter is a collection of utility functions to assist in mathematical
modeling in Pyomo.
"""

__author__ = "Alexander Dowling, Douglas Allan"


from operator import itemgetter

import pyomo.environ as pyo
from pyomo.core.expr.visitor import identify_variables
from pyomo.contrib.pynumero.interfaces.pyomo_nlp import PyomoNLP
import numpy as np
from scipy.linalg import svd
from scipy.sparse.linalg import svds, norm
from scipy.sparse import issparse, find

from idaes.core.util.model_statistics import (
    large_residuals_set,
    variables_near_bounds_set,
)
import idaes.core.util.scaling as iscale


class DegeneracyHunter:
    def __init__(self, block_or_jac, solver=None):
        """Initialize Degeneracy Hunter Object

        Arguments:
            block_or_jac: Pyomo model or Jacobian
            solver: Pyomo SolverFactory

        Notes:
            Passing a Jacobian to Degeneracy Hunter is current untested.

        """

        block_like = False
        try:
            block_like = issubclass(block_or_jac.ctype, pyo.Block)
        except AttributeError:
            pass

        if block_like:

            # Add Pyomo model to the object
            self.block = block_or_jac

            # setup pynumero interface
            self.nlp = PyomoNLP(self.block)

            # Get the scaled Jacobian of equality constraints
            self.jac_eq = iscale.get_jacobian(
                self.block, equality_constraints_only=True
            )[0]

            # Create a list of equality constraint names
            self.eq_con_list = self.nlp.get_pyomo_equality_constraints()
            self.var_list = self.nlp.get_pyomo_variables()

            self.candidate_eqns = None

        elif type(block_or_jac) is np.array:

            raise NotImplementedError(
                "Degeneracy Hunter currently only supports analyzing a Pyomo model"
            )

            # TODO: Need to refactor, document, and test support for Jacobian
            self.jac_eq = block_or_jac

            self.eq_con_list = None

        else:

            raise TypeError("Check the type for 'block_or_jac'")

        # number of equality constraints, variables
        self.n_eq = self.jac_eq.shape[0]
        self.n_var = self.jac_eq.shape[1]

        # Define default candidate equations (enumerate)
        candidate_eqns = range(self.n_eq)

        # Initialize solver
        if solver is None:
            # TODO: Test performance with open solvers such as cbc
            self.solver = pyo.SolverFactory("gurobi")
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
            block : model to be studied
            tol : residual threshold for inclusion in ComponentSet
            print_level: controls to extend of output to the screen
                0: nothing printed
                1: only name of constraint printed
                2: each constraint is pretty printed
                3: pretty print each constraint, then print value for included variable
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
            block : model to be studied
            tol : residual threshold for inclusion in ComponentSet (default = 1e-5)
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

        Argument:
            jac_eq Jacobian of equality constraints [matrix]
            M: largest value for nu

        Returns:
            m_dh: Pyomo model to calculate irreducible degenerate sets
        """

        n_eq = jac_eq.shape[0]
        n_var = jac_eq.shape[1]

        # Create Pyomo model for irreducible degenerate set
        m_dh = pyo.ConcreteModel()

        # Create index for constraints
        m_dh.C = pyo.Set(initialize=range(n_eq))

        m_dh.V = pyo.Set(initialize=range(n_var))

        # Add variables
        m_dh.nu = pyo.Var(m_dh.C, bounds=(-M, M), initialize=1.0)
        m_dh.y = pyo.Var(m_dh.C, domain=pyo.Binary)

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

        m_dh.degenerate = pyo.Constraint(m_dh.V, rule=eq_degenerate)

        def eq_lower(m_dh, c):
            return -M * m_dh.y[c] <= m_dh.nu[c]

        m_dh.lower = pyo.Constraint(m_dh.C, rule=eq_lower)

        def eq_upper(m_dh, c):
            return m_dh.nu[c] <= M * m_dh.y[c]

        m_dh.upper = pyo.Constraint(m_dh.C, rule=eq_upper)

        m_dh.obj = pyo.Objective(expr=sum(m_dh.y[c] for c in m_dh.C))

        return m_dh

    # TODO: Refactor, this should not be a staticmethod
    @staticmethod
    def _prepare_find_candidates_milp(jac_eq, M=1e5, m_small=1e-5):
        """
        Prepare MILP to find candidate equations for consider for IDS

        Argument:
            jac_eq Jacobian of equality constraints [matrix]
            M: maximum value for nu
            m_small: smallest value for nu to be considered non-zero

        Returns:
            m_fc: Pyomo model to find candidates
        """

        n_eq = jac_eq.shape[0]
        n_var = jac_eq.shape[1]

        # Create Pyomo model for irreducible degenerate set
        m_dh = pyo.ConcreteModel()

        # Create index for constraints
        m_dh.C = pyo.Set(initialize=range(n_eq))

        m_dh.V = pyo.Set(initialize=range(n_var))

        # Specify minimum size for nu to be considered non-zero
        m_dh.m_small = m_small

        # Add variables
        m_dh.nu = pyo.Var(m_dh.C, bounds=(-M - m_small, M + m_small), initialize=1.0)
        m_dh.y_pos = pyo.Var(m_dh.C, domain=pyo.Binary)
        m_dh.y_neg = pyo.Var(m_dh.C, domain=pyo.Binary)
        m_dh.abs_nu = pyo.Var(m_dh.C, bounds=(0, M + m_small))

        # Positive exclusive or negative
        def eq_pos_xor_negative(m, c):
            return m.y_pos[c] + m.y_neg[c] <= 1

        m_dh.pos_xor_neg = pyo.Constraint(m_dh.C)

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
                    return pyo.Constraint.Skip

        else:
            m_dh.J = jac_eq

            def eq_degenerate(m_dh, v):
                if np.sum(np.abs(m_dh.J[:, v])) > 1e-6:
                    return sum(m_dh.J[c, v] * m_dh.nu[c] for c in m_dh.C) == 0
                else:
                    # This variable does not appear in any constraint
                    return pyo.Constraint.Skip

        m_dh.pprint()

        m_dh.degenerate = pyo.Constraint(m_dh.V, rule=eq_degenerate)

        # When y_pos = 1, nu >= m_small
        # When y_pos = 0, nu >= - m_small
        def eq_pos_lower(m_dh, c):
            return m_dh.nu[c] >= -m_small + 2 * m_small * m_dh.y_pos[c]

        m_dh.pos_lower = pyo.Constraint(m_dh.C, rule=eq_pos_lower)

        # When y_pos = 1, nu <= M + m_small
        # When y_pos = 0, nu <= m_small
        def eq_pos_upper(m_dh, c):
            return m_dh.nu[c] <= M * m_dh.y_pos[c] + m_small

        m_dh.pos_upper = pyo.Constraint(m_dh.C, rule=eq_pos_upper)

        # When y_neg = 1, nu <= -m_small
        # When y_neg = 0, nu <= m_small
        def eq_neg_upper(m_dh, c):
            return m_dh.nu[c] <= m_small - 2 * m_small * m_dh.y_neg[c]

        m_dh.neg_upper = pyo.Constraint(m_dh.C, rule=eq_neg_upper)

        # When y_neg = 1, nu >= -M - m_small
        # When y_neg = 0, nu >= - m_small
        def eq_neg_lower(m_dh, c):
            return m_dh.nu[c] >= -M * m_dh.y_neg[c] - m_small

        m_dh.neg_lower = pyo.Constraint(m_dh.C, rule=eq_neg_lower)

        # Absolute value
        def eq_abs_lower(m_dh, c):
            return -m_dh.abs_nu[c] <= m_dh.nu[c]

        m_dh.abs_lower = pyo.Constraint(m_dh.C, rule=eq_abs_lower)

        def eq_abs_upper(m_dh, c):
            return m_dh.nu[c] <= m_dh.abs_nu[c]

        m_dh.abs_upper = pyo.Constraint(m_dh.C, rule=eq_abs_upper)

        # At least one constraint must be in the degenerate set
        m_dh.degenerate_set_nonempty = pyo.Constraint(
            expr=sum(m_dh.y_pos[c] + m_dh.y_neg[c] for c in m_dh.C) >= 1
        )

        # Minimize the L1-norm of nu
        m_dh.obj = pyo.Objective(expr=sum(m_dh.abs_nu[c] for c in m_dh.C))

        return m_dh

    # TODO: Refactor, this should not be a staticmethod
    @staticmethod
    def _check_candidate_ids(ids_milp, solver, c, eq_con_list=None, tee=False):
        """Solve MILP to check if equation 'c' is a significant component in an irreducible
        degenerate set

        Arguments:
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

        if pyo.check_optimal_termination(results):
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
            candidate_eqns: either None or list of indicies
            degenerate_set: either None or dictionary containing the degenerate_set
        """

        results = solver.solve(candidates_milp, tee=tee)

        if pyo.check_optimal_termination(results):
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
            Nothing

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

        Returns
        -------
        None.

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
        """Print variable, bounds, and value
        Argument:
            v: variable
        Return:
            nothing
        """
        print(v, "\t\t", v.lb, "\t", v.value, "\t", v.ub)
