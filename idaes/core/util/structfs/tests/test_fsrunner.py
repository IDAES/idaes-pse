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
from structured_notebook.fsrunner import FlowsheetRunner

# -- setup --

fsr = FlowsheetRunner()


@fsr.step("build")
def build_it(context):
    print("flowsheet - build")
    assert context.model.fs is not None
    add_units(context.model)


@fsr.substep("build", "add_units")
def add_units(m):
    print("flowsheet - add units (substep)")
    assert m.fs is not None


@fsr.step("add_costing")
def add_costing(context):
    print("flowsheet - costing")
    assert context.model is not None


@fsr.step("solve_optimization")
def solve_opt(context):
    print("flowsheet - solve")
    assert context.model is not None
    context["results"] = 123


# -- end setup --


@pytest.mark.unit
def test_run_all():
    fsr.run_steps()
    assert fsr.results == 123


@pytest.mark.unit
def test_rerun():
    fsr.run_steps()
    first_model = fsr.model
    # model not changed
    fsr.run_steps("solve_optimization")
    assert fsr.model == first_model
    # reset forces new model
    fsr.reset()
    fsr.run_steps("solve_optimization")
    assert fsr.model != first_model
    second_model = fsr.model
    # running from build also creates new model
    fsr.run_steps("build", "add_costing")
    assert fsr.model != second_model
