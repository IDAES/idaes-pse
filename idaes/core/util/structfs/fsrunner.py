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


class FlowsheetRunner(Runner):
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

    def run_steps(self, from_name: str = "", to_name: str = ""):
        """Run the steps.

        Before it calls the superclass to run the steps, checks
        if the step name matches the `build_step` attribute and,
        if so, creates an empty Pyomo ConcreteModel to use as
        the base model for the flowsheet.
        """
        from_step_name = self._norm_name(from_name)
        if (
            from_step_name == ""
            or from_step_name == self.build_step
            or self._context.model is None
        ):
            self._context.model = self._create_model()
        super().run_steps(from_name, to_name)

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
