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

import pytest

from idaes.core.util.structfs import FlowsheetRunner, _MockFlowsheetRunner, Steps, real_flowsheet_runner

@pytest.mark.unit
@pytest.mark.parametrize("clazz,mock", [(_MockFlowsheetRunner, True), (FlowsheetRunner, False)])
def test_flowsheet_runner(clazz, mock):
    if not mock and not real_flowsheet_runner:
        return  # do nothing

    print(f"\n** test {('real', 'mock')[mock]} flowsheet runner **")

    runner = clazz()
    calls = []

    @runner.step(Steps.build)
    def a(ctx):
        ctx.model = "model"
        calls.append("a")

    @runner.step(Steps.solve_optimization)
    def b(ctx):
        calls.append("b")

    # Call the decorated function and check that the mock method was called
    runner.run_steps()
    assert calls == ["a", "b"]
    