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
from enum import Enum

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

from idaes.core.util.model_serializer import to_json, from_json, StoreSpec, _only_fixed
from idaes.core.util.exceptions import InitializationError
from idaes.core.util.model_statistics import degrees_of_freedom
import idaes.logger as idaeslog

__author__ = "Andrew Lee"


class InitializationStatus(Enum):
    """
    Enum of expected outputs from Initialization routines.
    """

    Ok = 1  # Succesfully converged to tolerance
    Failed = 0  # Run, but failed to converge to tolerance
    DoF = -1  # Failed due to Degrees of Freedom issue
    Error = -2  # Exception raised during execution (other than DoF or convergence)


StoreState = StoreSpec(
    data_classes={
        Var._ComponentDataClass: (("fixed", "value"), _only_fixed),
        BooleanVar._ComponentDataClass: (("fixed",), None),
        Block._ComponentDataClass: (("active",), None),
        Constraint._ComponentDataClass: (("active",), None),
    }
)


class InitializerBase:
    """
    Base class for Initializer objects.

    This implements a default workflow nad methods for common tasks.
    Developers should feel free to overload these as necessary.
    """

    CONFIG = ConfigBlock()

    CONFIG.declare(
        "constraint_tolerance",
        ConfigValue(
            default=1e-5,
            domain=float,
            description="Tolerance for checking constraint convergence",
        ),
    )
    CONFIG.declare(
        "output_level",
        ConfigValue(
            default=idaeslog.NOTSET,
            description="Set output level for logging messages",
        ),
    )

    def __init__(self, **kwargs):
        self.config = self.CONFIG(kwargs)

        self.initial_state = None
        self.postcheck_summary = {}
        self.status = None

    def get_logger(self, model):
        return idaeslog.getInitLogger(model.name, self.config.output_level)

    def initialize(
        self, model: Block, initial_guesses: dict = None, json_file: str = None
    ):
        """
        Execute full initialization routine.

        Args:
            model: Pyomo model to be initialized.
            initial_guesses: dict of initial guesses to load.
            json_file: file name of json file to load initial guesses from as str.

        Note - can only provide one of initial_guesses or json_file.

        Returns:
            InitializationStatus Enum
        """
        # 1. Get current model state
        self.get_current_state(model)

        # 2. Load initial guesses
        self.load_initial_guesses(
            model, initial_guesses=initial_guesses, json_file=json_file
        )

        # 3. Fix states to make square
        self.fix_initialization_states(model)

        # 4. Prechecks
        self.precheck(model)

        # 5. try: Call specified initialization routine
        try:
            results = self.initialization_routine(model)
        # 6. finally: Restore model state
        finally:
            self.restore_model_state(model)

        # 7. Check convergence
        return self.postcheck(model, results_obj=results)

    def get_current_state(self, model: Block):
        """
        Get and store current state of variables (fixed/unfixed) and constraints/objectives
        activated/deactivated in model.

        Args:
            model: Pyomo model to get state from.

        Returns:
            dict serializing current model state.
        """
        self.initial_state = to_json(model, wts=StoreState, return_dict=True)

        return self.initial_state

    def load_initial_guesses(
        self, model: Block, initial_guesses: dict = None, json_file: str = None
    ):
        """
        Load initial guesses for variables into model.

        Args:
            model: Pyomo model to be initialized.
            initial_guesses: dict of initial guesses to load.
            json_file: file name of json file to load initial guesses from as str.

        Note - can only provide one of initial_guesses or json_file.

        Returns:
            None

        Raises:
            ValueError if both initial_guesses and json_file are provided.
        """
        if initial_guesses is not None and json_file is not None:
            self.status = InitializationStatus.Error
            raise ValueError(
                "Cannot provide both a set of initial guesses and a json file to load."
            )

        if initial_guesses is not None:
            self._load_values_from_dict(model, initial_guesses)
        elif json_file is not None:
            # TODO: What should we load here?
            # As this is just initialization, for now I am only loading variable values
            from_json(
                model, fname=json_file, wts=StoreSpec().value(only_not_fixed=True)
            )
        else:
            self.get_logger(model).info_high(
                "No initial guesses provided during initialization."
            )

    def fix_initialization_states(self, model: Block):
        """
        Call to model.fix_initialization_states method. Method will pass if
        fix_initialization_states not found.

        Args:
            model: Pyomo Block to fix states on.

        Returns:
            None
        """
        try:
            model.fix_initialization_states()
        except InitializationError:
            pass

    def precheck(self, model: Block):
        """
        Check for satisfied degrees of freedom before running initialization.

        Args:
            model: Pyomo Block to fix states on.

        Returns:
            None

        Raises:
            InitializationError if Degrees of Freedom do not equal 0.
        """
        if not degrees_of_freedom(model) == 0:
            self.status = InitializationStatus.DoF
            raise InitializationError(
                f"Degrees of freedom for {model.name} were not equal to zero during "
                f"initialization (DoF = {degrees_of_freedom(model)})."
            )

    def initialization_routine(self, model: Block):
        """
        Placeholder method to run initialization routine. Derived classes should overload
        this with the desired routine.

        Args:
            model: Pyomo Block to initialize.

        Returns:
            Overloaded method should return a Pyomo solver results object is available,
            otherwise None

        Raises:
            NotImplementedError
        """
        self.status = InitializationStatus.Error
        raise NotImplementedError()

    def restore_model_state(self, model: Block):
        """
        Restore model state to that stored in self.initial_state.

        Args:
            model: Pyomo Block ot restore state on.

        Returns:
            None

        Raises:
            ValueError if no initial state is stored.
        """
        if self.initial_state is not None:
            from_json(model, sd=self.initial_state, wts=StoreState)
        else:
            self.status = InitializationStatus.Error
            raise ValueError("No initial state stored.")

    def postcheck(self, model: Block, results_obj: dict = None):
        """
        Check the model has been converged after initialization.

        If a results_obj is provided, this will be checked using check_optimal_termination,
        otherwise this will walk all constraints in the model and check that they are within
        tolerance (set via the Initializer constraint_tolerance config argument).

        Args:
            model: model to be checked for convergence.
            results_obj: Pyomo solver results dict (if applicable, default=None).

        Returns:
            InitialationStatus Enum
        """
        if results_obj is not None:
            self.postcheck_summary = {
                "solver_status": check_optimal_termination(results_obj)
            }
            if not self.postcheck_summary["solver_status"]:
                # Final solver call did not return optimal
                self.status = InitializationStatus.Failed
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
                self.status = InitializationStatus.Failed
                raise InitializationError(
                    f"{model.name} failed to initialize successfully: uninitialized variables or "
                    "unconverged equality constraints detected. Please check postcheck summary for more information."
                )

        self.status = InitializationStatus.Ok
        return self.status

    def _load_values_from_dict(self, model, initial_guesses):
        """
        Internal method to iterate through items in initial_guesses and set value if Var and not fixed.
        """
        for c, v in initial_guesses.items():
            component = model.find_component(c)

            if not isinstance(component, (Var, _VarData)):
                self.status = InitializationStatus.Error
                raise TypeError(
                    f"Component {c} is not a Var. Initial guesses should only contain values for variables."
                )
            else:
                if component.is_indexed():
                    for i in component.values():
                        if not i.fixed:
                            i.set_value(v)
                elif not component.fixed:
                    component.set_value(v)

    def get_submodel_initializer(self, submodel: Block):
        """
        Lookup Initializer object to use for specified sub-model.

        This starts by checking the local mapping of objects and types, then falls back
        to object.default_initializer, followed by a global default (if defined).

        Args:
            submodel: sub-model to get default initializer for

        Returns:

        """
        # TODO: For MWE, return submodel - this will mean we run submodel.initialize()
        return submodel

        # TODO: Prototype code for getting initializer
        initializer = None

        if hasattr(submodel, "params"):
            # For StateBlocks and ReactionBlocks, look to the associated parameter block
            submodel = submodel.params

        if submodel in self.submodel_initializers:
            # First look for specific model instance
            initializer = self.submodel_initializers[submodel]
        elif type(submodel) in self.submodel_initializers:
            # Then look for types
            initializer = self.submodel_initializers[type(submodel)]
        else:
            # Then try the model's default initializer
            initializer = submodel.default_initializer

        if initializer is None:
            # If initializer is still None, try the master initializer's default
            initializer = self.default_submodel_initializer

        if initializer is None:
            # If we still have no initializer, log a warning and keep going
            self.get_logger(submodel).warning(
                "No Initializer found - attempting to continue."
            )

        return initializer
