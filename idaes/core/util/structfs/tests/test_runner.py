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
from .. import runner_actions

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


empty = Runner(("hi", "bye"))

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


@pytest.mark.unit
def test_run_steps_args():
    simple.run_steps(first="hello")
    simple.run_steps(last="world")


@pytest.mark.unit
def test_run_1step():
    simple.run_step("hello")


@pytest.mark.unit
def test_run_empty_steps():
    empty.run_steps()


@pytest.mark.unit
def test_run_bad_steps():
    with pytest.raises(KeyError):
        simple.run_steps("howdy", "pardner")

    with pytest.raises(KeyError):
        simple.run_step("notrun-1")


@pytest.mark.unit
def test_runner_context():
    simple.run_steps()
    assert simple["greeting"]


@pytest.mark.unit
def test_add_bad_step():
    with pytest.raises(KeyError):

        @simple.step("bad")
        def do_bad(ctx):
            return

    with pytest.raises(KeyError):

        @simple.substep("bad", "sub")
        def do_bad2(ctx):
            return

    # undefined step cannot have a substep

    with pytest.raises(ValueError):

        @simple.substep("notrun-1", "sub")
        def do_bad3(ctx):
            return


@pytest.mark.unit
def test_hellogoodbye():
    simple.add_action(
        "hg",
        runner_actions.HelloGoodbye,
        hello="Greetings and salutations",
        goodbye="Smell you later",
    )
    simple.run_steps(first="-", last="-")
    assert simple.get_action("hg").step_counter == 2
