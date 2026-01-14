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
'''
The 'flowsheet runner' is an API in the
{py:mod}`structfs` subpackage, and in
particular that package's {py:mod}`runner <structfs.runner>` and
{py:mod}`fsrunner <structfs.fsrunner>` modules.

## Overview

The core idea of the
{py:class}`FlowsheetRunner <structfs.fsrunner.FlowsheetRunner>` class is
that flowsheets should follow a standard set of "steps". By standardizing the
naming and ordering of these steps, it becomes easier to build tools that run
and inspect flowsheets. The Python mechanics of this are to put each step in a
function and wrap that function with decorator. The decorator uses a string to
indicate which standard step the function implements.

Once these functions are defined, the API can be used to execute and inspect a
wrapped flowsheet.

The framework can perform arbitrary actions before and after each run,
and before and after a given set of steps. This is implemented with 
the {py:class}`Actions <structfs.runner.Actions>` class 
and methods `add_action`, `get_action`, and `remove_action` on the base
{py:class}`Runner <structfs.runner.Runner>` class.
More details are given below in the Actions section.

## Step 1: Define flowsheet

It is assumed here that you have Python code to build, configure, and run an
IDAES flowsheet. You will first arrange this code to follow the standard "steps"
of a flowsheet workflow, which are listed in the
{py:class}`BaseFlowsheetRunner <structfs.fsrunner.BaseFlowsheetRunner>`
class' `STEPS` attribute. Not all the steps need to be defined: the API will
skip over steps with no definition when executing a range of steps. To make the
code more structured you can also define internal sub-steps, as described later.

The set of defined steps is:

* build - create the flowsheet
* set_operating_conditions - set initial variable values
* set_scaling - set scaling
* initialize - initialize the flowsheet
* set_solver - choose the solver
* solve_initial - perform an initial (square problem) solve
* add_costing - add costing information (if any)
* check_model_structure - check the model for structural issues
* initialize_costing - initialize costing variables
* solve_optimization - setup and solve the optimization problem
* check_model_numerics - check the model for numerical issues

### Example: Flash flowsheet

This is illustrated below with a before/after of an extremely simple flowsheet
with a single Flash unit model.

#### Before

For now, let's assume this flowsheet uses only four of the standard steps:
"build", "set_operating_conditions", "initialize", and "solve_optimization".
Let's also assume you have four functions defined that correspond to these
steps. Below is a sample flowsheet (for a single Flash unit) that we will use as
an example:

```{code} python
:name: before
from pyomo.environ import ConcreteModel, SolverFactory, SolverStatus
from idaes.core import FlowsheetBlock
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.models.unit_models import Flash

def build_model():
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = BTXParameterBlock(
        valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal", state_vars="FTPz"
    )
    m.fs.flash = Flash(property_package=m.fs.properties)
    return m


def set_operating_conditions(m):
    m.fs.flash.inlet.flow_mol.fix(1)
    m.fs.flash.inlet.temperature.fix(368)
    m.fs.flash.inlet.pressure.fix(101325)
    m.fs.flash.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.flash.inlet.mole_frac_comp[0, "toluene"].fix(0.5)
    m.fs.flash.heat_duty.fix(0)
    m.fs.flash.deltaP.fix(0)


def init_model(m):
    m.fs.flash.initialize()


def solve(m):
    solver = SolverFactory("ipopt")
    return solver.solve(m, tee=True)

```

#### After

In order to make this into a
{py:class}`FlowsheetRunner <structfs.fsrunner.FlowsheetRunner>`-wrapped
flowsheet, we need to do make a few changes. The modified file is shown below,
with changed lines highlighted and descriptions below.

```{code}
:name: after
:linenos:
:emphasize-lines: 7,9, 11, 25, 37, 43, 48, 12, 26, 38, 44, 49, 23, 28, 40, 44, 45, 46

from pyomo.environ import ConcreteModel, SolverFactory
from idaes.core import FlowsheetBlock
import idaes.logger as idaeslog
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE \
import ( BTXParameterBlock, )
from idaes.models.unit_models import Flash

from idaes.core.util.structfs.fsrunner import FlowsheetRunner

FS = FlowsheetRunner()

@FS.step("build") 
def build_model(ctx):
    """Build the model."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.properties = BTXParameterBlock(
        valid_phase=("Liq", "Vap"),
        activity_coeff_model="Ideal",
        state_vars="FTPz"
    )
    m.fs.flash = Flash(property_package=m.fs.properties) 
    # assert degrees_of_freedom(m) == 7
    ctx.model = m

@FS.step("set_operating_conditions")
def set_operating_conditions(ctx):
    """Set operating conditions."""
    m = ctx.model
    m.fs.flash.inlet.flow_mol.fix(1)
    m.fs.flash.inlet.temperature.fix(368)
    m.fs.flash.inlet.pressure.fix(101325)
    m.fs.flash.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.flash.inlet.mole_frac_comp[0, "toluene"].fix(0.5)
    m.fs.flash.heat_duty.fix(0)
    m.fs.flash.deltaP.fix(0)

@FS.step("initialize") 
def init_model(ctx): 
    """ "Initialize the model."""
    m = ctx.model
    m.fs.flash.initialize()

@FS.step("set_solver") 
def set_solver(ctx):
    """Set the solver."""
    ctx.solver = SolverFactory("ipopt")

@FS.step("solve_optimization")
def solve_opt(ctx):
    ctx["results"] = ctx.solver.solve(ctx.model, tee=ctx["tee"])
```

Details on the changes:

* **7**: Import the FlowsheetRunner class.
* **9**: Create a global {py:class}`FlowsheetRunner <structfs.fsrunner.FlowsheetRunner>` object, here called `FS`.
* **11, 25, 37, 43, 48**: Add a `@FS.step()` decorator in front of each function
  with the name of the associated step.
* **12, 26, 38, 44, 49**: Make each function take a single argument which is a {py:class}`fsrunner.Context <structfs.fsrunner.Context>` instance used to
  pass state information between functions (here, that argument is named `ctx`).
* **23**: Assign the model created in the "build" step to `ctx.model`, a
  standard attribute of the context object.
* **28, 40**: Replace the direct passing of the model object (in this case,
  called `m`) with a context object that has a `.model` attribute.
* **44-46**: Add a function for the `set_solver` step, to select the solver
  (here, IPOPT).
* **46**: In the "solve_optimization" step, assign the solver result to
  `ctx["results"]`.

## Step 2: Execute and inspect

Once the flowsheet has been 'wrapped' in the flowsheet runner interface, it can
be run and manipulated via the wrapper object. The basic steps to do this are:
import the flowsheet-runner object, build and execute the flowsheet, and inspect
the flowsheet.

For example to run all the steps and get the status of the solve, you
could do this:

```{code}
FS.run_steps()
assert FS.results.solver.status == SolverStatus.ok
```

Some more examples of using the FlowsheetRunner are shown in the
example notebooks found under the `docs/examples/structfs` directory
([docs link](/examples/structfs/index)).

## Actions  

```{autodoc2-docstring} structfs.runner.Action
```

## Annotation

You can also 'annotate' variables for special 
treatment in display, etc. with the
`annotate_var` function in the 
{py:class}`FlowsheetRunner <structfs.fsrunner.FlowsheetRunner>` class.

```{autodoc2-object} structfs.fsrunner.FlowsheetRunner.annotate_var
```

```{autodoc2-docstring} structfs.fsrunner.FlowsheetRunner.annotate_var
:parser: myst
```

'''
