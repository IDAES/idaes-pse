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
Mock structured flowsheet. This is used in place of the 'real' structured flowsheet
when the `idaes_fi` package is not installed.

Usage::

    from idaes.core.util.structfs import load_flowsheet_runner, Steps

    FlowsheetRunner = load_flowsheet_runner()

Regardless of whether the `FlowsheetRunner` is real or mocked, you can
use `run_steps()` to run the flowsheet programmatically::

    FS = FlowsheetRunner()

    @FS.step(Steps.build)
    def build(ctx):
        # ... create model 'm' ...
        ctx.model = m

    @FS.step(Steps.initialize)
    def initialize(ctx):
        m = ctx.model
        # ..etc..

    # ..etc..

    def main():
        FS.run_steps()

The rules for order of running are the same as in the library:

1. If no arguments are passed to the constructor, the names and order of steps
   will be the sequence in `Steps.index`.
2. If a sequence of strings is passed to the constructor, use this sequence as
   the names and order of steps to run.
3. If an explicit empty argument, e.g. `()` or `[]`, is passed to the
   constructor, run the steps in the order in which they are encountered.

For cases (1) and (2), any step name that is not in the sequence will result in a
KeyError being raised by the decorator function (`FS.step`).
````
"""

from typing import Type
from collections.abc import Sequence
from unittest.mock import Mock

__author__ = "Dan Gunter (LBNL)"


def load_flowsheet_runner(force_mock: bool = False) -> Type:
    """Load real or mock FlowsheetRunner.

    Args:
        force_mock: If True, return the mocked, do-nothing-different, wrapper
                    class even if the `idaes_fi` package is available.

    Returns:
        If the `idaes_fi` (flowsheet inspector) package is available, return the
        FlowsheetRunner class from that package. Otherwise, return a mock
        FlowsheetRunner class that will accept the same API but not do anything
        except call the wrapped functions when `run_steps()` is invoked (i.e., it will
        not provide information for the Flowsheet Inspector UI).
    """
    if force_mock:
        clazz = _MockFlowsheetRunner
    else:
        try:
            from idaes_fi.structfs import FlowsheetRunner  # noqa: F401

            clazz = FlowsheetRunner
        except ImportError:
            clazz = _MockFlowsheetRunner
    return clazz


# ----- end of public interface -----#


class _MockFlowsheetRunner:

    mock = True

    def __init__(self, *args, **kwargs):
        self._step_names = kwargs.get("steps", None)
        if self._step_names is None:
            # use pre-defined steps by default
            self._step_names = Steps.index
        elif not self._step_names:
            # normalize empty to tuple, meaning dynamic
            self._step_names = ()
        self.ctx = Mock()
        self._steps = {}

    def __getattr__(self, name):
        return Mock()

    def step(self, name: str, *args, **kwargs):
        def decorator(func):
            if self._step_names != () and name not in self._step_names:
                steppenlist = ", ".join(self._step_names)
                raise KeyError(f"Unknown step: '{name}' not in: {steppenlist}")
            self._steps[name] = func

            def wrapper(ctx=self.ctx):
                func(ctx)

            return wrapper

        return decorator

    def label(self, *args, **kwargs):
        def decorator(func):
            def wrapper(ctx=self.ctx):
                func(ctx)

            return wrapper

        return decorator

    def run_steps(self, *args, **kwargs):
        if self._step_names is None:
            for name in Steps.index:
                (func := self._steps.get(name, None)) and func(self.ctx)
        elif self._step_names:
            for name in self._step_names:
                (func := self._steps.get(name, None)) and func(self.ctx)
        else:  # dynamic case
            for func in self._steps.values():
                func(self.ctx)


class Steps:
    """Names of steps so that editor autocomplete, etc., will help
    to avoid typos.
    """

    build = "build"
    set_solver = "set_solver"
    initialize = "initialize"
    set_operating_conditions = "set_operating_conditions"
    set_scaling = "set_scaling"
    solve_initial = "solve_initial"
    set_autoscaling = "set_autoscaling"
    add_costing = "add_costing"
    initialize_costing = "initialize_costing"
    setup_optimization = "setup_optimization"
    solve_optimization = "solve_optimization"

    index = (
        build,
        set_solver,
        initialize,
        set_operating_conditions,
        set_scaling,
        solve_initial,
        set_autoscaling,
        add_costing,
        initialize_costing,
        setup_optimization,
        solve_optimization,
    )

    @classmethod
    def __len__(cls):
        return len(cls.index)
