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
"""
Specialize the generic `Runner` class to running a flowsheet,
in `FlowsheetRunner`.
"""
# stdlib
# none

# third-party
from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock

from idaes_connectivity.base import Connectivity
from idaes_connectivity.jupyter import display_connectivity

# package
from .runner import Runner


class Context(dict):
    """Syntactic sugar for the dictionary for the 'context' passed into each
    step of the `FlowsheetRunner` class.
    """

    @property
    def model(self):
        """The model being run."""
        return self["model"]

    @model.setter
    def model(self, value):
        """The model being run."""
        self["model"] = value

    @property
    def solver(self):
        """The solver used to solve the model."""
        return self["solver"]

    @solver.setter
    def solver(self, value):
        """The solver used to solve the model."""
        self["solver"] = value


class BaseFlowsheetRunner(Runner):
    """Specialize the base `Runner` to handle IDAES flowsheets.

    This class pre-determine the name and order of steps to run

    Attributes:
        STEPS: List of defined step names.
    """

    STEPS = (
        "build",
        "set_operating_conditions",
        "set_scaling",
        "initialize",
        "set_solver",
        "solve_initial",
        "add_costing",
        "check_model_structure",
        "initialize_costing",
        "solve_optimization",
        "check_model_numerics",
    )

    def __init__(self, solver=None, tee=False):
        self.build_step = self.STEPS[0]
        self._solver, self._tee = solver, tee
        super().__init__(self.STEPS)  # needs to be last

    def run_steps(self, first: str = "", last: str = ""):
        """Run the steps.

        Before it calls the superclass to run the steps, checks
        if the step name matches the `build_step` attribute and,
        if so, creates an empty Pyomo ConcreteModel to use as
        the base model for the flowsheet.
        """
        from_step_name = self.normalize_name(first)
        if (
            from_step_name == ""
            or from_step_name == self.build_step
            or self._context.model is None
        ):
            self._context.model = self._create_model()
        super().run_steps(first, last)

    def reset(self):
        self._context = Context(solver=self._solver, tee=self._tee, model=None)

    def _create_model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        return m

    @property
    def model(self):
        """Syntactic sugar to return the model."""
        return self._context.model

    @property
    def results(self):
        """Syntactic sugar to return the `results` in the context."""
        return self._context["results"]

    def mark(
        self,
        obj,
        title=None,
        desc=None,
        units=None,
        rounding=3,
        is_input=True,
        is_output=True,
        input_category="main",
        output_category="main",
    ):
        """Annotate a variable"""
        self._marks[obj] = {
            "fullname": obj.getname(),
            "title": title or obj.getname(),
            "description": desc or obj.getname(fully_qualified=True),
            "units": units or str(obj.getunits()),
            "rounding": rounding,
            "is_input": is_input,
            "is_output": is_output,
            "input_category": input_category,
            "output_category": output_category,
        }
        self._marks[obj].update(kwargs)

    def marks_to_table(self) -> list[list[str]]:
        """Translate 'marks' to a table suitable for loading up the UI (Flowsheet Processor).

        Returns:
            list of rows, each row is a list of strings. First row is the header:
                `name,obj,description,ui_units,display_units,rounding,is_input,input_category,is_output,output_category`
        """
        hdr_items = [
            "name",
            "obj",
            "description",
            "ui_units",
            "display_units",
            "rounding",
            "is_input",
            "input_category",
            "is_output",
            "output_category",
        ]
        # add header row
        tbl = [hdr_items]
        # add one row for each marked variable
        for m in self._marks:
            if "display_units" not in m:
                m["display_units"] = m["units"]
            row = [m["title"], m["fullname"], m["description"], m["units"]]
            for key in hdr_items[4:]:
                row.append(m[key])
            tbl.append(row)
        # return table
        return tbl


class FlowsheetRunner(BaseFlowsheetRunner):

    class DegreesOfFreedom:
        def __init__(self, runner):
            from .runner_actions import UnitDofChecker

            self._a = runner.add_action(
                "dof",
                UnitDofChecker,
                "fs",
                ["build", "solve_initial", "solve_optimization"],
            )
            self._rnr = runner

        def model(self):
            return self._a.get_dof_model()

        def __getattr__(self, name):
            """Naming the step prints a summary of that step."""
            if name not in set(self._rnr.list_steps()):
                raise AttributeError(f"No step named '{name}'")
            self._a.summary(step=name)

        def __str__(self):
            return self._a.summary(stream=None)

        def _ipython_display_(self):
            self._a.summary()

    class Timings:
        def __init__(self, runner):
            from .runner_actions import Timer

            self._a: Timer = runner.add_action("t", Timer)

        @property
        def values(self) -> list[dict]:
            return self._a.get_history()

        @property
        def history(self) -> str:
            h = []
            for i in range(len(self._a)):
                h.append(f"== Run {i + 1} ==")
                h.append("")
                h.append(self._a.summary(run_idx=i))
            return "\n".join(h)

        def __str__(self):
            return self._a.summary()

        def _ipython_display_(self):
            self._a._ipython_display_()

    def __init__(self, **kwargs):
        super().__init__(**kwargs)
        self.dof = self.DegreesOfFreedom(self)
        self.timings = self.Timings(self)

    def build(self):
        """Run just the build step"""
        self.run_step("build")

    def solve_initial(self):
        self.run_steps(last="solve_initial")

    def show_diagram(self):
        return display_connectivity(input_model=self.model)
