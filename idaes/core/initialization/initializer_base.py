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
Base class for initializer objects
"""
from pyomo.environ import check_optimal_termination, Constraint, value, Var
from pyomo.common.config import ConfigBlock, ConfigValue

from idaes.core.util.exceptions import InitializationError
from idaes.core.util.model_statistics import degrees_of_freedom

__author__ = "Andrew Lee"


class InitializerBase:
    CONFIG = ConfigBlock()

    CONFIG.declare(
        "constraint_tolerance",
        ConfigValue(
            default=1e-6,
            domain=float,
            description="Tolerance for checking constraint convergence",
        ),
    )

    def __init__(self):
        self.config = self.CONFIG()

        self.postcheck_summary = {}

    def initialize(self, model, initial_guesses=None):
        raise NotImplementedError()

    def get_current_state(self, model):
        pass

    def load_initial_guesses(self, model, initial_guesses):
        pass

    def fix_input_states(self, model):
        pass

    def precheck(self, model):
        if not degrees_of_freedom(model) == 0:
            raise InitializationError(
                f"Degrees of freedom for {model.name} were not equal to zero during "
                f"initialization (DoF = {degrees_of_freedom(model)})."
            )

    def initialize_model(self, model):
        raise NotImplementedError()

    def restore_model_state(self, model, initial_state):
        pass

    def postcheck(self, model, results_obj=None):
        if results_obj is not None:
            self.postcheck_summary = {"solver_status": check_optimal_termination(model)}
            if not self.postcheck_summary["solver_status"]:
                # Final solver call did not return optimal
                raise InitializationError(
                    f"{model.name} failed to initialize successfully: solver did not return "
                    "optimal termination. Please check the output logs for more information."
                )
        else:
            # Need to manually check initialization
            # First, check that all non-stale Vars have values
            uninit_vars = []
            for v in model.component_data_objects(Var, descend_into=True):
                if not v.stale and v.value is None:
                    uninit_vars.appeand(v)
            # Next check for unconverged equality constraints
            uninit_const = []
            for c in model.component_data_objects(Constraint, descend_into=True):
                if c.upper is not None and c.lower is not None and c.upper == c.lower:
                    if abs(value(c.body - c.lb)) >= self.config.constraint_tolerance:
                        uninit_const.append(c)

            self.postcheck_summary = {
                "uninitialized_vars": uninit_vars,
                "unconverged_constraints": uninit_const,
            }

            if len(uninit_const) > 0 or len(uninit_vars) > 0:
                raise InitializationError(
                    f"{model.name} failed to initialize successfully: uninitialized variables or "
                    "unconverged equality constraints detected. Please check postcheck summary for more information."
                )
