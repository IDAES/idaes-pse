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
# stdlib
from typing import Optional

# third-party
from pyomo.environ import ConcreteModel
from idaes.core import FlowsheetBlock

# package
from .runner import Runner


class Context(dict):
    @property
    def model(self):
        return self["model"]

    @model.setter
    def model(self, value):
        self["model"] = value

    @property
    def solver(self):
        return self["solver"]

    @solver.setter
    def solver(self, value):
        self["solver"] = value


class FlowsheetRunner(Runner):
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
        self._component_dof = {}

    def _create_model(self):
        m = ConcreteModel()
        m.fs = FlowsheetBlock(dynamic=False)
        return m

    @property
    def model(self):
        return self._context.model

    @property
    def results(self):
        return self._context["results"]
