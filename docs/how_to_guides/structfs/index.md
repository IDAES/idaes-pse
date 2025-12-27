# Flowsheet Runner

The 'flowsheet runner' is an API in the
{py:mod}`idaes.core.util.structfs` subpackage, and in
particular that package's {py:mod}`runner <idaes.core.util.structfs.runner>` and
{py:mod}`fsrunner <idaes.core.util.structfs.fsrunner>` modules.

## Overview

The core idea of the
{py:class}`FlowsheetRunner <idaes.core.util.structfs.fsrunner.FlowsheetRunner>` class is
that flowsheets should follow a standard set of "steps". By standardizing the
naming and ordering of these steps, it becomes easier to build tools that run
and inspect flowsheets. The Python mechanics of this are to put each step in a
function and wrap that function with decorator. The decorator uses a string to
indicate which standard step the function implements.

Once these functions are defined, the API can be used to execute and inspect a
wrapped flowsheet.

## Step 1: Define flowsheet

It is assumed here that you have Python code to build, configure, and run an
IDAES flowsheet. You will first arrange this code to follow the standard "steps"
of a flowsheet workflow, which are listed in the
{py:class}`BaseFlowsheetRunner <idaes.core.util.structfs.fsrunner.BaseFlowsheetRunner>`
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

```{code}
from pyomo.environ import ConcreteModel, SolverFactory from idaes.core
import FlowsheetBlock from
idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import
( BTXParameterBlock, ) from idaes.models.unit_models import Flash

def build_model(): m = ConcreteModel() m.fs = FlowsheetBlock(dynamic=False)
m.fs.properties = BTXParameterBlock( valid_phase=("Liq", "Vap"),
activity_coeff_model="Ideal", state_vars="FTPz" ) m.fs.flash =
Flash(property_package=m.fs.properties) return m

def set_operating_conditions(m): m.fs.flash.inlet.flow_mol.fix(1)
m.fs.flash.inlet.temperature.fix(368) m.fs.flash.inlet.pressure.fix(101325)
m.fs.flash.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
m.fs.flash.inlet.mole_frac_comp[0, "toluene"].fix(0.5)
m.fs.flash.heat_duty.fix(0) m.fs.flash.deltaP.fix(0)

def init_model(m): m.fs.flash.initialize(outlvl=idaeslog.INFO)

def solve(m): solver = SolverFactory("ipopt") return = solver.solve(m,
tee=ctx["tee"])
```

#### After

In order to make this into a
{py:class}`FlowsheetRunner <idaes.core.util.structfs.fsrunner.FlowsheetRunner>`-wrapped
flowsheet, we need to do make a few changes. The modified file is shown below,
with changed lines highlighted and descriptions below.

```{code}
:linenos:

from pyomo.environ import ConcreteModel, SolverFactory
from idaes.core import FlowsheetBlock

# Import idaes logger to set output levels
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
        valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal", state_vars="FTPz"
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
    m.fs.flash.initialize(outlvl=idaeslog.INFO)


@FS.step("set_solver") 
def set_solver(ctx):
    """Set the solver."""
    ctx.solver = SolverFactory("ipopt")


@FS.step("solve_optimization")
def solve_opt(ctx):
    ctx["results"] = ctx.solver.solve(ctx.model, tee=ctx["tee"])
```

Details on the changes:

* **10**: Import the FlowsheetRunner class
* **13**: Create a global `FlowsheetRunner` object, here called `FS`
* **15, 28, 41, 48, 54**: Add a `@FS.step()` decorator in front of each function
  with the name of the associated step
* **16, 29, 42, 39, 55**: Make each function take a single argument which is a
  {py:class}`fsrunner.Context <idaes.core.util.structfs.fsrunner.Context>` instance used to
  pass state information between functions (here, that argument is named `ctx`)
* **25**: Assign the model created in the "build" step to `ctx.model`, a
  standard attribute of the context object
* **31, 44**: Replace the direct passing of the model object (in this case,
  called `m`) with a context object that has a `.model` attribute
* **48-51**: Add a function for the `set_solver` step, to select the solver
  (here, IPOPT)
* **56**: In the "solve_optimization" step, assign the solver result to
  `ctx["results"]`

## Step 2: Execute and inspect

Once the flowsheet has been 'wrapped' in the flowsheet runner interface, it can
be run and manipulated via the wrapper object. The basic steps to do this are:
import the flowsheet-runner object, build and execute the flowsheet, and inspect
the flowsheet.

### Examples

Some examples of doing this are shown in the
[structfs examples notebooks](/examples/structfs/index).

