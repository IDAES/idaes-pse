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

from idaes.core.util.structfs import load_flowsheet_runner, Steps


@pytest.mark.unit
@pytest.mark.parametrize("mock", [True, False])
def test_flowsheet_runner(mock):
    FlowsheetRunner = load_flowsheet_runner(force_mock=mock)
    if not mock and hasattr(FlowsheetRunner, "mock"):
        return  # do nothing; real lib not available

    print(f"\n** test {('real', 'mock')[mock]} flowsheet runner **")

    runner = FlowsheetRunner()
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

@pytest.mark.unit
def test_step_order_user():
    rn = load_flowsheet_runner(force_mock=True)(steps=("a", "b", "c"))
    calls = []

    @rn.step("a")
    def step_a(ctx):
        calls.append("a")

    @rn.step("b")
    def step_b(ctx):
        calls.append("b")

    @rn.step("c")
    def step_c(ctx):
        calls.append("c")

    rn.run_steps(save_report=False)
    assert calls == ["a", "b", "c"]


@pytest.mark.unit
def test_step_order_default():
    rn = load_flowsheet_runner(force_mock=True)()
    calls = []

    @rn.step(Steps.initialize)
    def step_b(ctx):
        calls.append("b")

    @rn.step(Steps.solve_initial)
    def step_c(ctx):
        calls.append("c")

    @rn.step(Steps.build)
    def step_a(ctx):
        calls.append("a")

    rn.run_steps(save_report=False)
    assert calls == ["a", "b", "c"]


@pytest.mark.unit
def test_step_order_dynamic():
    rn = load_flowsheet_runner(force_mock=True)(steps=())
    calls = []

    @rn.step(Steps.initialize)
    def step_b(ctx):
        calls.append("b")

    @rn.step(Steps.solve_initial)
    def step_c(ctx):
        calls.append("c")

    @rn.step(Steps.build)
    def step_a(ctx):
        calls.append("a")

    rn.run_steps(save_report=False)
    assert calls == ["b", "c", "a"]  # dynamic order: as registered, not default


@pytest.mark.unit
def test_step_order_user_extra():
    rn = load_flowsheet_runner(force_mock=True)(steps=("a", "b"))
    calls = []

    @rn.step("a")
    def step_a(ctx):
        calls.append("a")

    @rn.step("b")
    def step_b(ctx):
        calls.append("b")

    with pytest.raises(KeyError):

        @rn.step("c")
        def step_c(ctx):
            calls.append("c")


@pytest.mark.unit
def test_step_order_default_extra():
    rn = load_flowsheet_runner(force_mock=True)()
    calls = []

    @rn.step(Steps.build)
    def step_a(ctx):
        calls.append("a")

    @rn.step(Steps.solve_optimization)
    def step_b(ctx):
        calls.append("b")

    with pytest.raises(KeyError):

        @rn.step("new_step")
        def step_c(ctx):
            calls.append("c")

    rn.run_steps(save_report=False)
    assert calls == ["a", "b"]  # only the registered steps run
