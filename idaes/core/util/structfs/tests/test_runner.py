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
import pytest
from ..runner import Runner

## -- setup --

simple = Runner(("notrun-1", "hello", "hello.dude", "world", "notrun-2"))


@simple.step("hello")
def say_hello(context):
    context["greeting"] = "Hello"
    dude("yo")


@simple.substep("hello", "dude")
def dude(s):
    print(f"{s}! this is called from hello, not directly by runner")


@simple.step("world")
def say_to_world(context):
    msg = f"{context['greeting']}, World!"
    print(msg)
    context["greeting"] = msg


# -- end setup --


@pytest.mark.unit
def test_simple_run_all():
    simple.run_steps()
    assert simple._context["greeting"] == "Hello, World!"


@pytest.mark.unit
def test_runner_actions():

    rn = Runner(("a step",))

    def do_nothing(context):
        print("do nothing")

    rn.add_action("nothing", do_nothing)
    rn.get_action("nothing")
    with pytest.raises(KeyError):
        rn.get_action("something")
    rn.remove_action("nothing")
    with pytest.raises(KeyError):
        rn.get_action("nothing")


@pytest.mark.unit
def test_run_steps_order():
    with pytest.raises(ValueError):
        simple.run_steps("world", "hello")
