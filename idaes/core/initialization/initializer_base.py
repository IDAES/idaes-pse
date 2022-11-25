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
from pyomo.environ import (
    BooleanVar,
    Block,
    check_optimal_termination,
    Constraint,
    Param,
    value,
    Var,
)
from pyomo.core.base.var import _VarData
from pyomo.common.config import ConfigBlock, ConfigValue

from idaes.core.util import to_json, from_json, StoreSpec
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog

__author__ = "Andrew Lee"

# Set up logger
_log = idaeslog.getLogger(__name__)


StoreState = StoreSpec(
    data_classes={
        Var._ComponentDataClass: (("fixed",), None),
        BooleanVar._ComponentDataClass: (("fixed",), None),
        Block._ComponentDataClass: (("active",), None),
        Constraint._ComponentDataClass: (("active",), None),
    }
)
StoreValues = StoreSpec(
    data_classes={
        Var._ComponentDataClass: (("value",), None),
        BooleanVar._ComponentDataClass: (("value",), None),
    }
)


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
        self.initial_state = {}

    def initialize(self, model, initial_guesses=None, json_file=None):
        # 1. Get current model state
        init_state = self.get_current_state(model)

        # 2. Load initial guesses
        self.load_initial_guesses(
            model, initial_guesses=initial_guesses, json_file=json_file
        )

        # 3. Fix states to make square
        self.fix_initialization_states(model)

        # 4. Prechecks
        self.precheck(model)

        # 5. try: Call block-triangularization solver
        try:
            self.initialization_routine(model)
        # 6. finally: Restore model state
        finally:
            self.restore_model_state(model)

        # 7. Check convergence
        self.postcheck(model)

    def get_current_state(self, model):
        self.initial_state = to_json(model, wts=StoreState, return_dict=True)

        return self.initial_state

    def load_initial_guesses(self, model, initial_guesses=None, json_file=None):
        if initial_guesses is not None and json_file is not None:
            raise ValueError(
                "Cannot provide both a set of initial guesses and a json file to load."
            )

        if initial_guesses is not None:
            self._load_values_from_dict(model, initial_guesses)
        elif json_file is not None:
            from_json(model, fname=json_file, wts=StoreValues)
        else:
            _log.info_high("No initial guesses provided during initialization.")

    def fix_initialization_states(self, model):
        model.fix_initialization_states()

    def precheck(self, model):
        if not degrees_of_freedom(model) == 0:
            raise InitializationError(
                f"Degrees of freedom for {model.name} were not equal to zero during "
                f"initialization (DoF = {degrees_of_freedom(model)})."
            )

    def initialization_routine(self, model):
        raise NotImplementedError()

    def restore_model_state(self, model):
        from_json(model, sd=self.initial_state, wts=StoreState)

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
                if v.value is None:
                    uninit_vars.append(v)
            # Next check for unconverged equality constraints
            uninit_const = []
            for c in model.component_data_objects(Constraint, descend_into=True):
                try:
                    if (
                        c.upper is not None
                        and c.lower is not None
                        and c.upper == c.lower
                    ):
                        if (
                            abs(value(c.body - c.lb))
                            >= self.config.constraint_tolerance
                        ):
                            uninit_const.append(c)
                    elif (
                        c.upper is not None
                        and value(c.body) > c.upper + self.config.constraint_tolerance
                    ):
                        uninit_const.append(c)
                    elif (
                        c.lower is not None
                        and value(c.body) < c.lower + self.config.constraint_tolerance
                    ):
                        uninit_const.append(c)
                except ValueError:
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

    def _load_values_from_dict(self, model, initial_guesses):
        for c, v in initial_guesses.items():
            component = model.find_component(c)

            if not isinstance(component, (Var, _VarData)):
                raise TypeError(
                    f"Component {c} is not a Var. Initial guesses should only contain values for variables."
                )
            else:
                if component.is_indexed():
                    for i in component.values():
                        i.set_value(v)
                else:
                    component.set_value(v)
