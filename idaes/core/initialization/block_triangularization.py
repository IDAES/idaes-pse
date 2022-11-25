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
"""
Initializer class for implementing Block Triangularization initialization
"""
from pyomo.environ import check_optimal_termination, Constraint, Var
from pyomo.common.config import ConfigValue
from pyomo.contrib.incidence_analysis.util import solve_strongly_connected_components

from idaes.core.initialization.initializer_base import InitializerBase
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.model_statistics import degrees_of_freedom
from idaes.core.solvers import get_solver

__author__ = "Andrew Lee"


class BlockTriangularizationInitializer(InitializerBase):
    CONFIG = InitializerBase.CONFIG()
    CONFIG.declare(
        "block_solver",
        ConfigValue(
            default=None,  # TODO: Use a square-problem solver here
            description="Solver to use for NxN blocks",
        ),
    )
    CONFIG.declare(
        "block_solver_options",
        ConfigValue(
            default={},
            description="Dict of options to pass to block solver",
        ),
    )
    CONFIG.declare(
        "calculate_variable_options",
        ConfigValue(
            default={},
            description="Dict of options to pass to 1x1 block solver",
            doc="Dict of options to pass to calc_var_kwds argument in "
            "solve_strongly_connected_components method",
        ),
    )

    def precheck(self, model):
        super().precheck(model)

        # TODO: Add addition checks here for perfect matching, etc.

    def initialization_routine(self, model):
        if self.config.block_solver is not None:
            solver = self.config.block_solver
        else:
            solver = get_solver()

        solve_strongly_connected_components(
            model,
            solver=solver,
            solve_kwds=self.config.block_solver_options,
            calc_var_kwds=self.config.calculate_variable_options,
        )
