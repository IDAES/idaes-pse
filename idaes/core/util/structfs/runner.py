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
Run functions in a module in a defined, named, sequence.
"""

# stdlib
import logging
from typing import Callable, Optional, Tuple, Sequence, TypeVar

__author__ = "Dan Gunter (LBNL)"

_log = logging.Logger(__name__)


class Step:
    """Step to run by the `Runner`."""

    SEP = "::"  # when printing out step::substep

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
            name: The name of substep
            func: Function to call to execute this substep
        """
        self.substeps.append((name, func))


# Python 3.9-compatible forward reference
ActionType = TypeVar("ActionType", bound="Action")


class Runner:
    """Run a set of defined steps."""

    STEP_ANY = "-"

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
            name: Substep name
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
        self.run_steps(first=name, last=name)

    def run_steps(
        self, first: str = "", last: str = "", after: str = "", before: str = ""
    ):
        """Run steps from `first`/`after` to step `last`/`before`.

           Specify only one of the first/after and last/before pairs.

           Use the special value `STEP_ANY` to mean the first or last defined step.

        Args:
            first: First step to run (include)
            after: Run first defined step after this one (exclude)
            last: Last step to run (include)
            before: Run last defined step before this one (exclude)

        Raises:
            KeyError: Unknown or undefined step given
            ValueError: Steps out of order or both first/after or before/last given
        """
        if first and after:
            raise ValueError("Cannot specify both 'after' and 'first'")
        if last and before:
            raise ValueError("Cannot specify both 'before' and 'last'")
        if not self._steps:
            return  # nothing to do, no steps defined
        args = (
            first or after,
            last or before,
            (bool(first) or not bool(after), bool(last) or not bool(before)),
        )
        return self._run_steps(*args)

    def _run_steps(
        self,
        first: str,
        last: str,
        endpoints: tuple[bool, bool],
    ):
        names = (self.normalize_name(first), self.normalize_name(last))
        print(f"@@ RUN STEPS: first={first} last={last} endpoints={endpoints}")

        # get indexes of first/last step
        step_range = [-1, -1]
        for i, step_name in enumerate(names):
            if step_name == self.STEP_ANY:  # meaning first or last defined
                # this will always find a step as long as there is at least one,
                # which we checked before calling this function
                idx = self._find_step(reverse=(i == 1))
            else:
                try:
                    idx = self._step_names.index(step_name)
                except ValueError:
                    raise KeyError(f"Unknown step: {step_name}")
                if step_name not in self._steps:
                    raise KeyError(f"Empty step: {step_name}")
            step_range[i] = idx

        # check that first comes before last
        if step_range[0] > step_range[1]:
            raise ValueError(
                "Steps out of order: {names[0]}={step_range[0]} > {names[1]}={step_range[1]}"
            )

        # execute overall before-run action
        for action in self._actions.values():
            action.before_run()

        # run each (defined) step
        for i in range(step_range[0], step_range[1] + 1):
            # check whether to skip endpoints in range
            if (i == step_range[0] and not endpoints[0]) or (
                i == step_range[1] and not endpoints[1]
            ):
                continue
            # get the step associated with the index
            step = self._steps.get(self._step_names[i], None)
            # if the step is defined, run it
            if step:
                step.func(self._context)

        # execute overall after-run action
        for action in self._actions.values():
            action.after_run()

    def reset(self):
        """Reset runner internal state, especially the context."""
        self._context = {}

    def list_steps(self, all_steps=False) -> list[str]:
        """Get list of [runnable] steps."""
        result = []
        for n in self._step_names:
            if all_steps or (n in self._steps):
                result.append(n)
        return result

    def add_action(self, name: str, action_class: type, *args, **kwargs) -> object:
        """Add a named action.

        Args:
            name: Arbitrary name for the action, used to get/remove it
            action_class: Subclass of Action to use
            args: Positional arguments passed to `action_class` constructor
            kwargs: Keyword arguments passed to `action_class` constructor
        """
        obj = action_class(self, *args, **kwargs)
        self._actions[name] = obj
        return obj

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

    def _find_step(self, reverse=False):
        start_step, end_step, incr = (
            (0, len(self._step_names), 1),
            (len(self._step_names) - 1, -1, -1),
        )[reverse]
        for i in range(start_step, end_step, incr):
            if self._step_names[i] in self._steps:
                return i
        return -1

    @classmethod
    def normalize_name(cls, s: Optional[str]) -> str:
        """Normalize a step name.
        Args:
            s: Step name

        Returns:
            normalized name
        """
        return cls.STEP_ANY if not s else s.lower()

    def _step_begin(self, name: str):
        for action in self._actions.values():
            action.before_step(name)

    def _substep_begin(self, base: str, name: str):
        for action in self._actions.values():
            action.before_substep(base, name)

    def _step_end(self, name: str):
        for action in self._actions.values():
            action.after_step(name)

    def _substep_end(self, base: str, name: str):
        for action in self._actions.values():
            action.after_substep(base, name)

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
                self._substep_begin(base, name)
                result = func(*args, **kwargs)
                self._substep_end(base, name)
                return result

            self.add_substep(base, name, wrapper)

            return wrapper

        return step_decorator


class Action:
    """The Action class implements a simple framework to run arbitrary
    functions before and/or after each step and/or run performed
    by the `Runner` class.

    To create and use your own Action, inherit from this class
    and then define one or more of the methods:

    * before_step - Called before a given step is executed
    * after_step - Called after a given step is executed
    * before/after_substep - Called before/after a named
      substep is executed (these can have arbitrary names)
    * before_run - Called before the first step is executed
    * after_run - Called after the last step is executed

    Then add the action to the `Runner` class (e.g., `FlowsheetRunner`)
    instance with `add_action()`. Note that you pass the action
    *class*, not instance. Additional settings can be passed to
    the created action instance with arguments to `add_action`.
    Also note that the *name* argument is used to retrieve the
    action instance later, as needed.

    ### Example

    Below is a simple example that prints a message
    before/after every step and prints the total number
    of steps run at the end of the run.

    ```{code}
    :name: hellogoodbye

    from idaes.core.util.structfs.runner import Action
    class HelloGoodbye(Action):
        "Example action, for tutorial purposes."

        def __init__(self, runner, hello="hi", goodbye="bye", **kwargs):
            super().__init__(runner, **kwargs)
            self._hello, self._goodbye = hello, goodbye
            self.step_counter = -1

        def before_run(self):
            self.step_counter = 0

        def before_step(self, name):
            print(f">> {self._hello} from step {name}")

        def before_substep(self, name, subname):
            print(f"  >> {self._hello} from sub-step {subname}")

        def after_step(self, name):
            print(f"<< {self._goodbye} from step {name}")
            self.step_counter += 1

        def after_substep(self, name, subname):
            print(f"  << {self._goodbye} from sub-step {subname}")

        def after_run(self):
            print(f"Ran {self.step_counter} steps")
    ```

    You could add the above example to a Runner subclass,
    here called `my_runner`, like this:

    ```{code}
    my_runner.add_action(
        "hg",
        HelloGoodbye,
        hello="Greetings and salutations",
        goodbye="Smell you later",
    )
    ```

    Then, after running steps, you could print
    the value of the *step_counter* attribute with:

    ```{code}
    print(my_runner.get_action("hg").step_counter)
    ```

    See the pre-defined actions in the
    {py:mod}`runner_actions <idaes.core.util.structfs.runner_actions>`
    module, and their usage in the `FlowsheetRunner` class, for more examples.
    """

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

    def before_substep(self, step_name: str, substep_name: str):
        return

    def after_step(self, step_name: str):
        """Perform this action after the named step.

        Args:
            step_name: Name of the step
        """
        return

    def after_substep(self, step_name: str, substep_name: str):
        return

    def before_run(self):
        """Perform this action before a run starts."""
        return

    def after_run(self):
        """Perform this action after a run ends."""
        return

    def report(self) -> dict:
        return {}
