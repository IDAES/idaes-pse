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
Initializer class for implementing Block Triangularization initialization
"""

from pyomo.environ import check_optimal_termination
from pyomo.common.config import Bool, ConfigDict, ConfigValue
from pyomo.contrib.incidence_analysis import (
    IncidenceGraphInterface,
    solve_strongly_connected_components,
)

from idaes.core.initialization.initializer_base import (
    InitializerBase,
    InitializationStatus,
)
from idaes.core.solvers import get_solver
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.model_statistics import number_activated_constraints

__author__ = "Andrew Lee"


class BlockTriangularizationInitializer(InitializerBase):
    """
    Block Triangularization based Initializer object.

    This Initializer should be suitable for most models, but may struggle to initialize
    tightly coupled systems of equations.

    This Initializer uses the common workflow defined in InitializerBase and calls
    the Pyomo solve_strongly_connected_components function to initialize the model.

    """

    # TODO: Block solver is IPOPT for now, as fsolve struggles with VLE
    CONFIG = InitializerBase.CONFIG()
    CONFIG.declare(
        "block_solver",
        ConfigValue(
            default="ipopt_v2",
            description="Solver to use for NxN blocks",
        ),
    )
    CONFIG.declare(
        "block_solver_options",
        ConfigDict(
            implicit=True,
            description="Dict of options to pass to block solver",
            doc="Dict of options to use to set solver.options.",
        ),
    )
    CONFIG.block_solver_options.declare(
        "tol",
        ConfigValue(
            default=1e-8,
            domain=float,
            description="Convergence tolerance for block solver",
        ),
    )
    CONFIG.block_solver_options.declare(
        "max_iter",
        ConfigValue(
            default=200,
            domain=int,
            description="Iteration limit for block solver",
        ),
    )
    CONFIG.declare(
        "block_solver_writer_config",
        ConfigDict(
            implicit=True,
            description="Dict of writer_config arguments to pass to block solver",
        ),
    )
    CONFIG.block_solver_writer_config.declare(
        "linear_presolve",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Whether to use linear presolver with block solver",
        ),
    )
    CONFIG.block_solver_writer_config.declare(
        "scale_model",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Whether to apply model scaling with block solver",
        ),
    )
    CONFIG.declare(
        "block_solver_call_options",
        ConfigDict(
            implicit=True,
            description="Dict of arguments to pass to solver.solve call",
            doc="Dict of arguments to be passed as part of the solver.solve "
            "call, such as tee=True.",
        ),
    )
    CONFIG.declare(
        "calculate_variable_options",
        ConfigDict(
            implicit=True,
            description="Dict of options to pass to 1x1 block solver",
            doc="Dict of options to pass to calc_var_kwds argument in "
            "solve_strongly_connected_components method.",
        ),
    )
    CONFIG.declare(
        "skip_final_solve",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Skip solving the block after block triangularization finishes",
            doc="Skip solving the block after block triangularization finishes."
            "This solve may be necessary to decrease the scaled constraint residuals "
            "below the specified tolerance because, until Pyomo issue #3785 is addressed, "
            "solve_strongly_connected_components does not take scaling into account.",
        ),
    )

    def precheck(self, model):
        """
        Check for perfect matching in model.

        If this fails, it indicates a structural singularity in the model.
        """
        super().precheck(model)

        if model.is_indexed():
            for d in model.values():
                self._check_matching(d)
        else:
            self._check_matching(model)

    def _check_matching(self, block_data):
        """
        Run incidence analysis on given block data and check matching.
        """
        igraph = IncidenceGraphInterface(block_data, include_inequality=False)
        matching = igraph.maximum_matching()
        if len(matching) != len(igraph.variables):
            self._update_summary(
                block_data, "status", InitializationStatus.PrecheckFailed
            )
            raise InitializationError(
                f"Perfect matching not found for {block_data.name}. "
                f"This suggests that the model is structurally singular."
            )

    def initialization_routine(self, model):
        """
        Call Block Triangularization solver on model.
        """
        solver = get_solver(
            self.config.block_solver,
            solver_options=self.config.block_solver_options,
            writer_config=self.config.block_solver_writer_config,
        )

        if model.is_indexed():
            for d in model.values():
                self._solve_block_data(d, solver)
        else:
            self._solve_block_data(model, solver)

    def _solve_block_data(self, block_data, solver):
        """
        Call solve_strongly_connected_components on a given BlockData.
        """
        if number_activated_constraints(block_data) == 0:
            # Nothing to solve. Additionally, the .nl writer
            # throws an exception if a solver is called on a
            # block with no active constraints.
            return

        # TODO: Can we get more diagnostic output from this method?
        results_list = solve_strongly_connected_components(
            block_data,
            solver=solver,
            solve_kwds=self.config.block_solver_call_options,
            calc_var_kwds=self.config.calculate_variable_options,
        )
        for results in results_list:
            if results is not None and not check_optimal_termination(results):
                raise InitializationError(
                    f"Block Triangularization failed with solver status: {results['Solver']}."
                )

        # Until Pyomo issue #3785 is addressed, solve_strongly_connected_components
        # does not take scaling into account. However, the postcheck does, so we
        # perform an extra solve to make sure all constraint residuals are small enough.
        if not self.config.skip_final_solve:
            results = solver.solve(block_data)
            if not check_optimal_termination(results):
                raise InitializationError(
                    f"Could not solve {block_data.name} after block triangularization finished."
                )
