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
This module contains an implementation of the Degeneracy Hunter algorithm in Pyomo and IDAES. See this paper for an explanation of how the algorithm works:

Dowling, A. W., & Biegler, L. T. (2015). Degeneracy hunter: An algorithm for determining irreducible sets of degenerate constraints in mathematical programs. In Computer Aided Chemical Engineering (Vol. 37, pp. 809-814). Elsevier.
"""

__author__ = "Alexander Dowling, Douglas Allan, Andrew Lee, Robby Parker, Ben Knueven"

import sys

import numpy as np
from scipy.sparse import issparse, find

from pyomo.environ import (
    Binary,
    check_optimal_termination,
    ConcreteModel,
    Constraint,
    Objective,
    Set,
    SolverFactory,
    Var,
)
from pyomo.core.base.block import BlockData
from pyomo.common.config import (
    ConfigDict,
    ConfigValue,
    document_kwargs_from_configdict,
    NonNegativeFloat,
)

from idaes.core.util.model_statistics import (
    greybox_block_set,
)
from idaes.core.scaling.util import (
    get_jacobian,
)
import idaes.logger as idaeslog

from idaes.core.util.diagnostics_tools.writer_utils import (
    MAX_STR_LENGTH,
    TAB,
)

_log = idaeslog.getLogger(__name__)


# Constants for Degeneracy Hunter
YTOL = 0.9
MMULT = 0.99


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


@document_kwargs_from_configdict(DHCONFIG)
class DegeneracyHunter:
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
