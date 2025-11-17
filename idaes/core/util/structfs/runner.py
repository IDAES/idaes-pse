###############################################################################
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
###############################################################################
"""
Run functions in a module in a defined, named, sequence.
"""

# stdlib
import logging
from typing import Callable, Optional, Tuple, Sequence, TypeVar

__author__ = "Dan Gunter (LBNL)"

_log = logging.Logger(__name__)


class Step:
    """Step to run by the `Runner`."""

    def __init__(self, name: str, func: Callable):
        """Constructor

        Args:
            name: Name of the step
            func: Function to call to execute the step
        """
        self.name: str = name
        self.func: Callable = func
        self.substeps: list[Tuple[str, Callable]] = []

    def add_substep(self, name: str, func: Callable):
        """Add a sub-step to this step.
        Substeps are run in the order given.

        Args:
            name: Name of substep
            func: Function to call to execute this substep
        """
        self.substeps.append((name, func))


# Python 3.9-compatible forward reference
ActionType = TypeVar("ActionType", bound="Action")


class Runner:
    """Run a set of defined steps."""

    def __init__(self, steps: Sequence[str]):
        """Constructor.

        Args:
            steps: List of step names
        """
        self._context = {}
        self._actions: dict[str, ActionType] = {}
        self._step_names = list(steps)
        self._steps: dict[str, Step] = {}
        self.reset()

    def __getitem__(self, key):
        """Look for key in `context`"""
        return self._context[key]

    def add_step(self, name: str, func: Callable):
        """Add a step.

        Steps are executed by calling `func(context)`,
        where `context` is a dict (or dict-like) object
        that is used to pass state between steps.

        Args:
            name: Add a step to be executed
            func: Function to execute for the step.

        Raises:
            KeyError: _description_
        """
        step_name = self.normalize_name(name)

        if step_name not in self._step_names:
            raise KeyError(f"Unknown step: {step_name}")
        self._steps[step_name] = Step(step_name, func)

    def add_substep(self, base_name, name, func):
        """Add a substep for a given step.

        Substeps are all executed, in the order added,
        immediately after their base step is executed.

        Args:
            base_name: Step name
            name: Substep name_
            func: Function to execute

        Raises:
            KeyError: Base step or substep is not found
            ValueError: Base step does not have any substeps
        """
        substep_name = self.normalize_name(name)
        base_step_name = self.normalize_name(base_name)
        if base_step_name not in self._step_names:
            raise KeyError(
                f"Unknown base step {base_step_name} for substep {substep_name}"
            )
        try:
            step = self._steps[base_step_name]
        except KeyError:
            raise ValueError(
                f"Empty base step {base_step_name} for substep {substep_name}"
            )
        step.add_substep(substep_name, func)

    def run_step(self, name):
        """Syntactic sugar for calling `run_steps` for a single step."""
        self.run_steps(from_name=name, to_name=name)

    def run_steps(self, from_name: str = "", to_name: str = ""):
        """Run steps from `from_name` to step `to_name`.

        Args:
            from_name: First step to run
            to_name: Last step to run

        Raises:
            KeyError: Unknown or undefined step given
            ValueError: Steps out of order (`from` after `to`)
        """
        if not self._steps:
            return  # nothing to do, no steps defined

        names = (self.normalize_name(from_name), self.normalize_name(to_name))

        step_range = [-1, -1]
        for i, step_name in enumerate(names):
            if step_name == "":
                idx = self._first_step() if i == 0 else self._last_step()
            else:
                try:
                    idx = self._step_names.index(step_name)
                except ValueError:
                    raise KeyError(f"Unknown step: {step_name}")
                if step_name not in self._steps:
                    raise KeyError(f"Empty step: {step_name}")
            step_range[i] = idx

        if step_range[0] > step_range[1]:
            raise ValueError(
                "Steps out of order: {names[0]}={step_range[0]} > {names[1]}={step_range[1]}"
            )

        for action in self._actions.values():
            action.before_run()

        for i in range(step_range[0], step_range[1] + 1):
            step = self._steps.get(self._step_names[i], None)
            if step:
                step.func(self._context)

        for action in self._actions.values():
            action.after_run()

    def reset(self):
        """Reset runner internal state, especially the context."""
        self._context = {}

    def add_action(self, name: str, action_class: type, *args, **kwargs):
        """Add a named action.

        Args:
            name: Arbitrary name for the action, used to get/remove it
            action_class: Subclass of Action to use
            args: Positional arguments passed to `action_class` constructor
            kwargs: Keyword arguments passed to `action_class` constructor
        """
        obj = action_class(self, *args, **kwargs)
        self._actions[name] = obj

    def get_action(self, name: str) -> ActionType:
        """Get an action object.

        Args:
            name: Name of action (as provided to `add_action`)

        Returns:
            ActionType: Action object

        Raises:
            KeyError: If action name does not match any known action
        """
        return self._actions[name]

    def remove_action(self, name: str):
        """Remove an action object.

        Args:
            name: Name of action (as provided to `add_action`)

        Raises:
            KeyError: If action name does not match any known action
        """
        del self._actions[name]

    def _first_step(self):
        for i, name in enumerate(self._step_names):
            if name in self._steps:
                return i
        assert False, "No first step defined"  # should not get here

    def _last_step(self):
        for i in range(len(self._step_names) - 1, -1, -1):
            name = self._step_names[i]
            if name in self._steps:
                return i
        return -1

    @staticmethod
    def normalize_name(s: Optional[str]) -> str:
        """Normalize a step name.
        Args:
            s: Step name

        Returns:
            normalized name
        """
        return "" if s is None else s.lower()

    def _step_begin(self, name: str):
        for action in self._actions.values():
            action.before_step(name)

    def _step_end(self, name: str):
        for action in self._actions.values():
            action.after_step(name)

    def step(self, name: str):
        """Decorator function for creating a new step.

        Args:
            name: Step name

        Returns:
            Decorator function.
        """

        def step_decorator(func):

            def wrapper(*args, **kwargs):
                self._step_begin(name)
                result = func(*args, **kwargs)
                self._step_end(name)
                return result

            self.add_step(name, wrapper)

            return wrapper

        return step_decorator

    def substep(self, base: str, name: str):
        """Decorator function for creating a new substep.

        Substeps are not run directly, and must have an already
        existing base step as their parent.

        Args:
            base: Base step name
            name: Substep name

        Returns:
            Decorator function.
        """

        def step_decorator(func):

            def wrapper(*args, **kwargs):
                self._step_begin(name)
                result = func(*args, **kwargs)
                self._step_end(name)
                return result

            self.add_substep(base, name, wrapper)

            return wrapper

        return step_decorator


class Action:
    """Do something before and/or after each step and/or run performed by a `Runner`."""

    def __init__(self, runner: Runner, log: Optional[logging.Logger] = None):
        """Constructor

        Args:
            runner: Reference to the runner that will trigger this action.
            log: Logger to use when logging informational or error messages
        """
        self._runner = runner
        if log is None:
            log = _log
        self.log = log

    def before_step(self, step_name: str):
        """Perform this action before the named step.

        Args:
            step_name: Name of the step
        """
        return

    def after_step(self, step_name: str):
        """Perform this action after the named step.

        Args:
            step_name: Name of the step
        """
        return

    def before_run(self):
        """Perform this action before a run starts."""
        return

    def after_run(self):
        """Perform this action after a run ends."""
        return
