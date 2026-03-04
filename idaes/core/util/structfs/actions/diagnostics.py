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
"""
Model diagnostics
"""

from io import StringIO
import re

import pydantic as pc
from idaes.core.util.model_diagnostics import (
    DiagnosticsToolbox,
    _collect_model_statistics,
)


class StructuralIssueData(pc.BaseModel):
    kind: str = pc.Field(description="Issue kind (warning or caution)")
    category: str = pc.Field(default="other", description="Normalized issue category")
    count: int | None = pc.Field(default=None, description="Count parsed from message")
    message: str = pc.Field(default="", description="Primary issue message")
    details: list[str] = pc.Field(
        default_factory=list, description="Additional detail lines for the issue"
    )


class StructuralData(pc.BaseModel):
    statistics: list[str] = pc.Field(
        default_factory=list, description="Model statistics from structural diagnostics"
    )
    warnings: int = pc.Field(
        default=0, description="Number of structural warnings identified by diagnostics"
    )
    warning_details: list[StructuralIssueData] = pc.Field(
        default_factory=list, description="Structured details for structural warnings"
    )
    cautions: int = pc.Field(
        default=0, description="Number of structural cautions identified by diagnostics"
    )
    caution_details: list[StructuralIssueData] = pc.Field(
        default_factory=list, description="Structured details for structural cautions"
    )
    next_steps: list[str] = pc.Field(
        default_factory=list, description="Suggested structural next-step commands"
    )
    report: str = pc.Field(
        default="", description="Raw structural diagnostics report text"
    )


# class NumericalData(pc.BaseModel):
#     """Numerical issues"""

#     pass


class SolverResultData(pc.BaseModel):
    """Solver result"""

    status: bool = pc.Field(description="General status, True if OK")
    termination_condition: str = pc.Field(
        "Provides the specific reason the solver stopped (e.g., 'optimal', 'infeasible', 'unbounded', 'maxTime')."
    )
    objective_value: float = pc.Field(
        default=0,
        description="The optimal objective function value (if a feasible solution was found",
    )


class ModelDiagnosticsData(pc.BaseModel):
    """Data fetched from the IDAES Diagnostics Toolkit"""

    result: SolverResultData = pc.Field(description="Solver result", default={})
    structural: StructuralData = pc.Field(description="Structural issues")


#   numerical: NumericalData = pc.Field(description="Numerical issues", default={})


class ModelDiagnostics:
    """Run the model diagnostics and parse result."""

    def __init__(self):
        self._model = None
        self.structural_issues: bool = False
        self.numerical_issues: bool = False
        self._data = None  # ModelDiagnosticsData()

    def generate(self, model, results) -> ModelDiagnosticsData:
        """Generate the model diagnostics for current model state.

        Returns:
            The current data from the diagnostics
        """
        self._model = model

        # Extract results from Pyomo results object
        solve_ok = results.solver.status.lower() == "ok"
        solver_result = SolverResultData(
            status=solve_ok,
            termination_condition=results.solver.termination_condition,
        )
        if solve_ok and hasattr(results.problem, "objective"):
            solver_result.objective_value = results.problem.objective.value

        d = DiagnosticsToolbox(model)

        print(f"@@ self.structural_issues={self.structural_issues}")
        # Extract structural issues (uses _functions, thus pylint pragma)
        # pylint: disable=W0212
        if self.structural_issues:
            warnings, next_steps = d._collect_structural_warnings()
            cautions = d._collect_structural_cautions()

            stream = StringIO()
            d.report_structural_issues(stream=stream)
            report = stream.getvalue()

            struct = StructuralData()
            struct.statistics = _collect_model_statistics(model)
            struct.warnings = len(warnings)
            struct.warning_details = [
                self._parse_issue_line(msg, kind="warning") for msg in warnings
            ]
            struct.cautions = len(cautions)
            struct.caution_details = [
                self._parse_issue_line(msg, kind="caution") for msg in cautions
            ]
            struct.next_steps = next_steps
            struct.report = report
        else:
            struct = StructuralData()

        self._data = ModelDiagnosticsData(result=solver_result, structural=struct)
        return self._data

    @property
    def data(self) -> ModelDiagnosticsData:
        """Get the model diagnostics data, as populated by last run."""
        return self._data

    @staticmethod
    def _parse_issue_line(raw: str, kind: str) -> StructuralIssueData:
        lines = [ln.strip() for ln in raw.splitlines() if ln.strip()]
        first = lines[0] if lines else raw.strip()
        body = first.split(":", 1)[1].strip() if ":" in first else first
        match = re.search(r"(\d+)", body)
        count = int(match.group(1)) if match is not None else None

        msg_lower = body.lower()
        category = "other"
        if "degree" in msg_lower and "freedom" in msg_lower:
            category = "degrees_of_freedom"
        elif "inconsistent units" in msg_lower:
            category = "inconsistent_units"
        elif "structural singularity" in msg_lower:
            category = "structural_singularity"
        elif "evaluation errors" in msg_lower:
            category = "evaluation_errors"
        elif "fixed to 0" in msg_lower:
            category = "fixed_to_zero"
        elif "unused variable" in msg_lower or "unused variables" in msg_lower:
            category = "unused_variables"

        return StructuralIssueData(
            kind=kind,
            category=category,
            count=count,
            message=body,
            details=lines[1:],
        )
