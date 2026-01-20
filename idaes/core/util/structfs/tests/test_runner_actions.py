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

    s = timer.get_history()
    # [ {'run': 0.8005404472351074,
    #    'steps': {
    #       'step0': 0.10010385513305664,
    #       'step2': 0.3000965118408203,
    #       'step1': 0.20009303092956543
    #      },
    #     'inclusive': 0.6002933979034424,
    #     'exclusive': 0.20024704933166504
    #    },
    #    ...
    # ]
    eps = 0.1  # big variance needed for Windows
    for r in s:
        print(f"Timings: {r}")
        assert r["run"] == approx(0.8, abs=eps)
        assert r["inclusive"] + r["exclusive"] == approx(r["run"])
        for name, t in r["steps"].items():
            assert t == approx(0.1 + 0.1 * i, abs=eps)


@pytest.mark.component
def test_timer_runner():
    rn = runner.Runner(["step1", "step2", "step3"])

    @rn.step("step1")
    def sleepy1(context):
        time.sleep(0.1)

    @rn.step("step2")
    def sleepy2(context):
        time.sleep(0.1)

    @rn.step("step3")
    def sleepy3(context):
        time.sleep(0.1)

    rn.add_action("timer", Timer)

    rn.run_steps()

    s = rn.get_action("timer").get_history()

    eps = 0.1  # big variance needed for Windows
    for r in s:
        print(f"Timings: {r}")
        assert r["run"] == approx(0.3, abs=eps)
        assert r["inclusive"] + r["exclusive"] == approx(r["run"])
        for name, t in r["steps"].items():
            # assert name == f"step{i + 1}"
            assert t == approx(0.1, abs=eps)


@pytest.mark.unit
def test_unit_dof_action_base():
    rn = flash_flowsheet.FS
    rn.reset()

    def check_step(name, data):
        print(f"check_step {name} data: {data}")
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

    pprint.pprint(rn.get_action("check_dof").get_dof())


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
        step_dof = act.get_dof()[s]
        assert step_dof
        dofs.append(step_dof)
    assert dofs[0] != dofs[1]

    assert act.steps() == act.steps(only_with_data=True)


@pytest.mark.unit
def test_timer_report():
    rn = flash_flowsheet.FS
    rn.reset()
    rn.add_action("timer", Timer)
    rn.run_steps()
    report = rn.get_action("timer").report()
    # {'build': 0.053082942962646484,
    # 'set_operating_conditions': 0.0004742145538330078,
    # 'initialize': 0.22397446632385254,
    # 'set_solver': 7.581710815429688e-05,
    # 'solve_initial': 0.03623509407043457}
    expect_steps = (
        "build",
        "set_operating_conditions",
        "initialize",
        "set_solver",
        "solve_initial",
    )
    assert report
    for step_name in expect_steps:
        assert step_name in report
        assert report[step_name] < 1


@pytest.mark.unit
def test_dof_report():
    rn = flash_flowsheet.FS
    rn.reset()
    check_steps = (
        "build",
        "set_operating_conditions",
        "initialize",
        "solve_initial",
    )
    rn.add_action("dof", UnitDofChecker, "fs", check_steps)
    rn.run_steps()
    report = rn.get_action("dof").report()
    print(f"@@ REPORT:\n{report}")
    assert report
    report_data = report.model_dump()
    assert report_data["model"] == 0  # model has DOF=0
    for step_name in check_steps:
        assert step_name in report_data["steps"]
        for unit, value in report_data["steps"][step_name].items():
            assert value >= 0  # DOF > 0 in all (step, unit)
