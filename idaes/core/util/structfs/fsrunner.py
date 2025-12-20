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
from typing import Union

# third-party
from pyomo.environ import ConcreteModel, is_variable_type
from pyomo.environ import units as pyunits
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
        self._ann = {}
        super().__init__(self.STEPS)  # needs to be last

    def run_steps(self, first: str = "", last: str = "", **kwargs):
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
        super().run_steps(first, last, **kwargs)

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

    def annotate_var(
        self,
        variable: object,
        key: str = None,
        title: str = None,
        desc: str = None,
        units: str = None,
        rounding: int = 3,
        is_input: bool = True,
        is_output: bool = True,
        input_category: str = "main",
        output_category: str = "main",
    ) -> object:
        """Annotate a Pyomo variable.

        Args:
            variable: Pyomo variable being annotated
            key: Key for this block in dict. Defaults to object name.
            title: Name / title of the block. Defaults to object name.
            desc: Description of the block. Defaults to object name.
            units: Units. Defaults to string value of native units.
            rounding: Significant digits
            is_input: Is this variable an input
            is_output: Is this variable an output
            input_category: Name of input grouping to display under
            output_category: Name of output grouping to display under

        Returns:
            Input block (for chaining)

        Raises:
            ValueError: if `is_input` and `is_output` are both False
        """
        if not is_input and not is_output:
            raise ValueError("One of 'is_input', 'is_output' must be True")

        qual_name = variable.name
        key = key or variable.name

        self._ann[key] = {
            "var": variable,
            "fullname": qual_name,
            "title": title or qual_name,
            "description": desc or qual_name,
            "units": units or str(pyunits.get_units(variable)),
            "rounding": rounding,
            "is_input": is_input,
            "is_output": is_output,
            "input_category": input_category,
            "output_category": output_category,
        }

        return variable

    @property
    def annotated_vars(self) -> dict[str,]:
        """Get (a copy of) the annotated variables."""
        return self._ann.copy()


class FlowsheetRunner(BaseFlowsheetRunner):

    class DegreesOfFreedom:
        """Wrapper for the UnitDofChecker action"""

        def __init__(self, runner):
            from .runner_actions import UnitDofChecker

            # check DoF after build, initial solve, and optimization solve
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
        """Wrapper for the Timer action"""

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

    def run_steps(self, after=None, before=None, **kwargs):
        if after is not None:
            if "first" in kwargs:
                raise ValueError("Cannot specify both 'after' and 'first'")
            kwargs["first"] = before
            if "endpoints" in kwargs:
                ep = list(kwargs["endpoints"])
                ep[0] = False
                kwargs["endpoints"] = tuple(ep)
            else:
                kwargs["endpoints"] = (False, True)
        if before is not None:
            if "last" in kwargs:
                raise ValueError("Cannot specify both 'before' and 'last'")
            kwargs["last"] = before
            if "endpoints" in kwargs:
                ep = list(kwargs["endpoints"])
                ep[1] = False
                kwargs["endpoints"] = tuple(ep)
            else:
                kwargs["endpoints"] = (True, False)
        return super().run_steps(**kwargs)

    def show_diagram(self):
        return display_connectivity(input_model=self.model)
