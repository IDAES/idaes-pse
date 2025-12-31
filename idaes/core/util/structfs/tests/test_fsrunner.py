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
from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock
from ..fsrunner import FlowsheetRunner
from .flash_flowsheet import FS as flash_fs

# -- setup --

fsr = FlowsheetRunner()


@fsr.step("build")
def build_it(context):
    print("flowsheet - build")
    context.model = ConcreteModel()
    print(f"@@ build_it: id(model)={id(context.model)}")
    context.model.fs = FlowsheetBlock(dynamic=False)
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
    print(f"@@ test: id(model)={id(fsr.model)}")

    print("-- rerun --")

    # model not changed
    fsr.run_steps("solve_optimization")
    assert fsr.model == first_model


@pytest.mark.unit
def test_rerun_reset():
    fsr.run_steps()
    print(f"@@ test: id(model)={id(fsr.model)}")
    first_model = fsr.model

    print("-- rerun --")

    # reset forces new model
    fsr.reset()
    fsr.run_steps("solve_optimization")
    assert fsr.model != first_model
    second_model = fsr.model


@pytest.mark.unit
def test_rerun_frombuild():
    fsr.run_steps()
    first_model = fsr.model
    print(f"@@ test: id(model)={id(fsr.model)}")

    print("-- rerun --")

    # running from build also creates new model
    fsr.run_steps("build", "add_costing")
    assert fsr.model != first_model


@pytest.mark.unit
def test_annotation():
    runner = flash_fs
    runner.run_steps("build")
    print(runner.timings.history)

    ann = runner.annotate_var  # alias
    flash = runner.model.fs.flash  # alias
    category = "flash"
    kw = {"input_category": category, "output_category": category}

    ann(
        flash.inlet.flow_mol,
        key="fs.flash.inlet.flow_mol",
        title="Inlet molar flow",
        desc="Flash inlet molar flow rate",
        **kw,
    ).fix(1)
    ann(flash.inlet.temperature, units="Centipedes", **kw).fix(368)
    ann(flash.inlet.pressure, **kw).fix(101325)
    ann(flash.inlet.mole_frac_comp[0, "benzene"], **kw).fix(0.5)
    ann(flash.inlet.mole_frac_comp[0, "toluene"], **kw).fix(0.5)
    ann(flash.heat_duty, **kw).fix(0)
    ann(flash.deltaP, is_input=False, **kw).fix(0)

    ann = runner.annotated_vars
    print("-" * 40)
    print(ann)
    print("-" * 40)
    assert ann["fs.flash.inlet.flow_mol"]["title"] == "Inlet molar flow"
    assert (
        ann["fs.flash.inlet.flow_mol"]["description"] == "Flash inlet molar flow rate"
    )
    assert ann["fs.flash.inlet.flow_mol"]["input_category"] == category
    assert ann["fs.flash.inlet.flow_mol"]["output_category"] == category
    assert runner.model.fs.flash.inlet.flow_mol[0].value == 1
    assert ann["fs.flash._temperature_inlet_ref"]["units"] == "Centipedes"
    assert ann["fs.flash.deltaP"]["is_input"] == False
