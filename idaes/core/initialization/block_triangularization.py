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

from idaes.core.initialization.initializer_base import InitializerBase
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.model_statistics import degrees_of_freedom

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

    def initialize(self, model, initial_guesses=None):
        # 1. Get current model state
        init_state = self.get_current_state(model)

        # 2. Load initial guesses
        self.load_initial_guesses(model, initial_guesses)

        # 3. Fix states to make square
        self.fix_input_states(model)

        # 4. Prechecks
        self.precheck(model)

        # 5. try: Call block-triangularization solver
        try:
            self.initialize_model(model)
        # 6. finally: Restore model state
        finally:
            self.restore_model_state(model, init_state)

        # 7. Check convergence
        self.postcheck(model)

    def initialize_model(self, model):
        pass
