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
Mock structured flowsheet.
This is used in place of the real structured flowsheet when
the idaes_fi package is not installed.

Usage::

    from idaes.core.util.structfs import FlowsheetRunner, Steps

Regardless of whether the `FlowsheetRunner` is real or mocked, you can
use `run_steps()` to run the flowsheet programmatically.
In the case of the mock, it will run each step in the order defined in the file.

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
````
"""

from unittest.mock import Mock

__author__ = "Dan Gunter (LBNL)"


class _MockFlowsheetRunner:
    def __init__(self, *args, **kwargs):
        self._args, self._kwargs = args, kwargs
        self.ctx = Mock()
        self._steps = []

    def __getattr__(self, name):
        return Mock()

    def step(self, *args, **kwargs):
        def decorator(func):
            self._steps.append(func)

            def wrapper(ctx=self.ctx):
                func(ctx)

            return wrapper

        return decorator

    def label(self, *args, **kwargs):
        pass

    def run_steps(self, *args, **kwargs):
        for func in self._steps:
            func(self.ctx)


try:
    from idaes_fi.structfs import FlowsheetRunner  # noqa: F401
    from idaes_fi.structfs.common import Steps  # noqa: F401

    real_flowsheet_runner = True
except ImportError:
    FlowsheetRunner = _MockFlowsheetRunner
    Steps = Mock()
    real_flowsheet_runner = False
