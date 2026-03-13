# -*- coding: utf-8 -*-
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
Print diagnostics values in Pydantic data models.
"""

import sys

from idaes.core.util.compute_diagnostics import (
    StructuralCautionsData,
    StructuralIssuesData,
    StructuralWarningsData,
    VariableListData,
)


class Report:
    """Interface class to pretty-print the contents of the associated data model
    to the console.

    Subclasses should implement `get_lines()` to return a title
    and list of lines for the content.
    """

    def print(self, stream=None, indent=4, width=80, border=True):
        # use stdout if no stream is given
        if stream is None:
            stream = sys.stdout
        # setup
        tab = " " * indent

        # write top border
        if border:
            stream.write("=" * width)
            stream.write("\n")

        # get title and content
        title, lines = self.get_lines()

        if lines is None:
            stream.write(title)
            stream.write("\n")
        else:
            # write title and divider
            if title:
                stream.write(title)
                stream.write("\n")
                stream.write("-" * width)
                stream.write("\n")

            # write content
            for line in lines:
                stream.write(f"{tab}{line}\n")

        # write bottom divider
        if border:
            stream.write("=" * width)
            stream.write("\n")

    def get_lines(self) -> tuple[str, list[str]]:
        return "", []  # override in subclasses


def _plural(n, word):
    s = "s" if abs(n) > 1 else ""
    return f"{n} {word}{s}"


class VariableListReport(Report):

    def __init__(self, data: VariableListData):
        self._data = data

    def get_lines(self) -> tuple[str, list[str] | None]:
        d = self._data
        n = len(d.variables)
        if n == 0:
            title = f"No model variables {d.description}"
            return title, None
        title = f"Model variables that {d.description} ({n})"
        lines = []
        for i in range(n):
            items = [d.variables[i]]
            if d.details[i]:
                items.append(d.details[i])
            if d.values[i] is not None:
                items.append(f"value={format(d.values[i], d.value_format)}")
            line = " ".join(items)
            lines.append(line)
        return title, lines


class StructuralWarningsReport(Report):

    def __init__(self, data: StructuralWarningsData):
        self._data = data

    def get_lines(self):
        d = self._data
        lines = []
        if d.dof is not None:
            lines.append(f"WARNING: {_plural(d.dof, 'Degree')} of Freedom")
        if d.inconsistent_units is not None:
            lines.append(
                f"WARNING: {_plural(len(d.inconsistent_units), 'Component')} with inconsistent units"
            )
        if d.underconstrained_set is not None or d.overconstrained_set is not None:
            t = " " * 4
            ucv, ucc = len(d.underconstrained_set.variables), len(
                d.underconstrained_set.constraints
            )
            ocv, occ = len(d.overconstrained_set.variables), len(
                d.overconstrained_set.constraints
            )
            lines.extend(
                [
                    f"WARNING: Structural singularity found",
                    f"{t}Under-Constrained Set: {ucv} variables, {ucc} constraints",
                    f"{t}Over-Constrained Set: {ocv} variables, {occ} constraints",
                ]
            )
        if d.evaluation_errors is not None:
            lines.append(
                f"WARNING: Found {len(d.evaluation_errors)} potential evaluation errors."
            )

        return "Structural warnings", lines


class StructuralCautionsReport(Report):

    def __init__(self, data: StructuralCautionsData):
        self._data = data

    def get_lines(self):
        d = self._data
        lines = []
        if d.zero_vars is not None:
            lines.append(f"Caution: {_plural(len(d.zero_vars), 'variable')} fixed to 0")
        if d.unused_vars_free is not None or d.unused_vars_fixed is not None:
            nfree = 0 if d.unused_vars_free is None else len(d.unused_vars_free)
            nfixed = 0 if d.unused_vars_fixed is None else len(d.unused_vars_fixed)
            lines.append(
                f"Caution: {_plural(nfree + nfixed, 'unused variable')} ({nfixed} fixed)"
            )
        return "Structural cautions", lines


class StructuralIssuesReport(Report):

    def __init__(self, data: StructuralCautionsData):
        self._data = data

    def get_lines(self):
        d = self._data
        wt, wlines = d.warnings.get_lines()
        ct, clines = d.cautions.get_lines()
        return (
            "Structural issues",
            [wt, "-" * len(wt)] + wlines + ["", ct, "-" * len(ct)] + clines,
        )
