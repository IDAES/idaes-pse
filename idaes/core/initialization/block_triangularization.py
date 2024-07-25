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
Initializer class for implementing Block Triangularization initialization
"""
from pyomo.environ import SolverFactory
from pyomo.common.config import Bool, ConfigDict, ConfigValue
from pyomo.contrib.incidence_analysis import (
    IncidenceGraphInterface,
    solve_strongly_connected_components,
)

from idaes.core.initialization.initializer_base import (
    InitializerBase,
    InitializationStatus,
)
from idaes.core.util.exceptions import InitializationError

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
            default=False,
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
        # TODO: For now, go directly through solver factory as default solver
        # options cause failures. Most of these appear to be due to scaling,
        # so hopefully we can fix these later.
        solver = SolverFactory(
            self.config.block_solver,
            options=self.config.block_solver_options,
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
        # TODO: Can we get diagnostic output from this method?
        solve_strongly_connected_components(
            block_data,
            solver=solver,
            solve_kwds=self.config.block_solver_call_options,
            calc_var_kwds=self.config.calculate_variable_options,
        )
