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
Initializer class for implementing Block Triangularization initialization
"""
from pyomo.environ import SolverFactory
from pyomo.common.config import ConfigDict, ConfigValue
from pyomo.contrib.incidence_analysis import (
    IncidenceGraphInterface,
    solve_strongly_connected_components,
)

from idaes.core.initialization.initializer_base import (
    InitializerBase,
    InitializationStatus,
)
from idaes.core.util.exceptions import InitializationError
from idaes.core.solvers import get_solver

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
            default="ipopt",
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
    CONFIG.declare(
        "block_solver_call_options",
        ConfigDict(
            implicit=True,
            description="Dict of arguments to pass to solver.solve call",
            doc="Dict of arguments to be passed as part of the solver.solve "
            "call, such as tee=True/",
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
        if self.config.block_solver is not None:
            solver = SolverFactory(self.config.block_solver)
            solver.options.update(self.config.block_solver_options)
        else:
            solver = get_solver(options=self.config.block_solver_options)

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
