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
import pprint
import time

import pytest
from pytest import approx
from .. import runner
from ..runner_actions import Timer, UnitDofChecker
from . import flash_flowsheet


@pytest.mark.unit
def test_class_timer():
    timer = Timer(runner.Runner([]))
    n, m = 2, 3
    for i in range(n):
        timer.before_run()
        time.sleep(0.1)
        for j in range(m):
            name = f"step{j}"
            timer.before_step(name)
            time.sleep(0.1 * (j + 1))
            timer.after_step(name)
        time.sleep(0.1)
        timer.after_run()

    s = timer.summary()
    # tests/test_runner_actions.py [
    # {'run': 0.8005404472351074,
    # 'steps': [('step0', 0.10010385513305664), ('step1', 0.20009303092956543),
    # ('step2', 0.3000965118408203)], 'inclusive': 0.6002933979034424, 'exclusive': 0.20024704933166504},
    #
    # {'run': 0.8004975318908691, 'steps': [('step0', 0.10008621215820312),
    # ('step1', 0.20008587837219238),
    # ('step2', 0.3000912666320801)], 'inclusive': 0.6002633571624756, 'exclusive': 0.20023417472839355}]
    eps = 0.1  # big variance needed for Windows
    for r in s:
        print(f"Timings: {r}")
        assert r["run"] == approx(0.8, abs=eps)
        assert r["inclusive"] + r["exclusive"] == approx(r["run"])
        for i, (name, t) in enumerate(r["steps"]):
            assert name == f"step{i}"
            assert t == approx(0.1 + 0.1 * i, abs=eps)


@pytest.mark.component
def test_timer_runner():
    rn = runner.Runner(["step1", "step2", "step3"])

    def sleepy(context):
        time.sleep(0.1)

    rn.add_step("step1", sleepy)
    rn.add_step("step2", sleepy)
    rn.add_step("step3", sleepy)

    rn.add_action("timer", Timer)

    rn.run_steps()

    s = rn.get_action("timer").summary()

    eps = 0.1  # big variance needed for Windows
    for r in s:
        print(f"Timings: {r}")
        assert r["run"] == approx(0.3, abs=eps)
        assert r["inclusive"] + r["exclusive"] == approx(r["run"])
        for i, (name, t) in enumerate(r["steps"]):
            assert name == f"step{i + 1}"
            assert t == approx(0.1, abs=eps)


@pytest.mark.unit
def test_unit_dof_action_base():
    rn = flash_flowsheet.FS
    rn.reset()

    def check_step(name, data):
        # print(f"@@ check_step {name} data: {data}")
        assert "fs.flash" in data
        if name == "solve_initial":
            assert data["fs.flash"] == 0

    def check_run(step_dof, model_dof):
        assert model_dof == 0

    rn.add_action(
        "check_dof",
        UnitDofChecker,
        "fs",
        ["build", "solve_initial"],
        check_step,
        check_run,
    )

    rn.run_steps("build", "solve_initial")

    pprint.pprint(rn.get_action("check_dof").as_dataframe())


@pytest.mark.unit
def test_unit_dof_action_getters():
    rn = flash_flowsheet.FS
    rn.reset()

    aname = "check_dof"
    rn.add_action(
        aname,
        UnitDofChecker,
        "fs",
        ["build", "solve_initial"],
    )
    rn.run_steps()

    act = rn.get_action(aname)

    steps = act.steps()
    dofs = []
    for s in steps:
        step_dof = act.get_unit_dof(s)
        assert step_dof
        dofs.append(step_dof)
    assert dofs[0] != dofs[1]

    assert act.steps() == act.steps(only_with_data=True)
