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
from typing import Callable, Optional, Tuple, Sequence

__author__ = "Dan Gunter (LBNL)"

_log = logging.Logger(__name__)


class Step:
    def __init__(self, name: str, func: Callable):
        self.name: str = name
        self.func: Callable | None = func
        self.substeps: list[Tuple[str, Callable]] = []

    def add_substep(self, name: str, func: Callable):
        self.substeps.append((name, func))


type ActionType = Action


class Runner:
    """Run a set of defined steps."""

    def __init__(self, steps: Sequence[str]):
        """Constructor.

        Args:
            steps: List of step names
        """
        self._actions: dict[str, ActionType] = {}
        self._step_names = list(steps)
        self._steps: dict[str, Step] = {}
        self.reset()

    def __getitem__(self, key):
        return self._context[key]

    def add_step(self, name: str, func: Callable):
        step_name = self._norm_name(name)

        if step_name not in self._step_names:
            raise KeyError(f"Unknown step: {step_name}")
        self._steps[step_name] = Step(step_name, func)

    def add_substep(self, base_name, name, func):
        substep_name = self._norm_name(name)
        base_step_name = self._norm_name(base_name)
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

        names = (self._norm_name(from_name), self._norm_name(to_name))

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
        self._context = {}

    def add_action(self, name: str, action_class: type, *args, **kwargs):
        obj = action_class(self, *args, **kwargs)
        self._actions[name] = obj

    def get_action(self, name: str) -> ActionType:
        return self._actions[name]

    def remove_action(self, name: str):
        del self._actions[name]

    def _step_index(self, name: str):
        return self._step_names.index(name)

    def _first_step(self):
        for i, name in enumerate(self._step_names):
            if name in self._steps:
                return i
        return -1

    def _last_step(self):
        for i in range(len(self._step_names) - 1, -1, -1):
            name = self._step_names[i]
            if name in self._steps:
                return i
        return -1

    @staticmethod
    def _norm_name(s: str | None) -> str:
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
    """Do something before and/or after each step / run.

    ```
    class GreetingAction(Action):
        def before_step(self, name, runner):
            print(f"Hello, step {name}")
        def after_step(self, name, runner):
            print(f"Goodbye, step {name}")

    class DivaAction(Action):
        def after_step(self, name, runner):
            if name == "act":
                print(f"I am now an ac-TORRRR!")

    myrunner = Runner(steps=["plan", "act", "inspect", "revise"])
    myrunner.add_action(GreetingAction)
    myrunner.add_action(DivaAction)
    ```
    """

    def __init__(self, runner: Runner, log: logging.Logger | None = None):
        self._runner = runner
        if log is None:
            log = _log
        self.log = log

    def before_step(self, step_name: str):
        return

    def after_step(self, step_name: str):
        return

    def before_run(self):
        return

    def after_run(self):
        return
