Flowsheet Runner
=================

.. py:currentmodule:: idaes.core.util.structfs

The 'flowsheet runner' is an API in the :py:mod:`idaes.core.util.structfs` subpackage, and in particular the :py:mod:`.runner` and :py:mod:`.fsrunner` modules.

The core idea of the :py:class:`fsrunner.FlowsheetRunner` class is that the code that does this should follow a standard set of "steps", each of which can be represented with a function.
By standardizing the naming and ordering of these steps, it becomes easier to build tools that run and inspect flowsheets.

There are two stages to using the FlowsheetRunner API:

1. wrapping the code that builds and runs a flowsheet and
2. executing and inspecting a wrapped flowsheet.

Create FlowsheetRunner
-----------------------
It is assumed here that you have Python code to build, configure, and run an IDAES flowsheet.
You will first arrange this code to follow the standard "steps" of a flowsheet workflow, which are
listed in the :py:attr:`fsrunner.BaseFlowsheetRunner.STEPS` class attribute.
Not all the steps need to be defined; the API will skip over steps with no definition when executing a range of steps.
To make the code more structured you can define sub-steps, as described later.


Before
++++++

For now, let's assume a simple flowsheet with four steps: "build", "set_operating_conditions", "initialize", and "solve_optimization". 
Let's also assume you have four functions defined that correspond to these steps.
Below is a sample flowsheet (for a single Flash unit) that we will use as an example:

.. code-block:: python

    from pyomo.environ import ConcreteModel, SolverFactory
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
        m.fs.flash.initialize(outlvl=idaeslog.INFO)

    def solve(m):
        solver = SolverFactory("ipopt")
        return = solver.solve(m, tee=ctx["tee"])

After
+++++

In order to make this into a :py:class:`fsrunner.FlowsheetRunner`-wrapped flowsheet, we need to do make a few changes.
The modified file is shown below, with changed lines highlighted and descriptions below.

.. code-block:: python
    :linenos:
    :emphasize-lines: 10, 13, 15, 16, 25, 28, 29, 31, 31, 42, 44, 48, 49, 51, 54, 55, 56

    from pyomo.environ import ConcreteModel, SolverFactory
    from idaes.core import FlowsheetBlock

    # Import idaes logger to set output levels
    import idaes.logger as idaeslog
    from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
        BTXParameterBlock,
    )
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

Details on the changes:

==============  ===========
Lines           Description
==============  ===========
10              Import the FlowsheetRunner class
13              Create a global ``FlowsheetRunner`` object, here called ``FS``
15,28,41,48,54  Add a ``@FS.step()`` decorator in front of each function with the name of the associated step
16,29,42,39,55  Make each function take a single argument which is a :py:class:`fsrunner.Context` instance used to pass state information between functions (here, that argument is named ``ctx``)
31,44           Replace the direct passing of the model object (in this case, called `m`) with a context object that has a ``.model`` attribute
48-51           Add a function for the ``set_solver`` step, to select the solver (here, IPOPT)
25              Assign the model created in the "build" step to ``ctx.model``, a standard attribute of the context object
31,44           At the top of all other steps, alias the ``ctx.model`` variable to the local variable name used in the original function (in this case, ``m``)
56              In the "solve_optimization" step, assign the solver result to ``ctx["results"]``
==============  ===========


Executing and inspecting FlowsheetRunner
-----------------------------------------
Once the flowsheet has been 'wrapped' in the flowsheet runner interface, it can be run and manipulated via the wrapper object.
The basic steps to do this are: import the flowsheet-runner object, build and execute the flowsheet, and inspect the flowsheet.

Some examples of doing this are shown in the :doc:`structfs examples notebooks </examples/structfs/index>`.
