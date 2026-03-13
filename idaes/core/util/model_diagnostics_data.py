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
"""Pydantic data structures for model diagnostics reporting."""

import sys
from pydantic import BaseModel, Field


class Printable:
    """Interface class to pretty-print the contents to the console.

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


class VariableListData(BaseModel, Printable):
    tag: str
    description: str
    #: Names of variables
    variables: list[str] = Field(default_factory=list)
    #: Optional values of variables
    values: list[float | None] = Field(default_factory=list)
    value_format: str = ".6E"
    #: Optional descriptive details for variables
    details: list[str] = Field(default_factory=list)

    def get_lines(self) -> tuple[str, list[str] | None]:
        n = len(self.variables)
        if n == 0:
            title = f"No model variables {self.description}"
            return title, None
        title = f"Model variables that {self.description} ({n})"
        lines = []
        for i in range(n):
            items = [self.variables[i]]
            if self.details[i]:
                items.append(self.details[i])
            if self.values[i] is not None:
                items.append(f"value={format(self.values[i], self.value_format)}")
            line = " ".join(items)
            lines.append(line)
        return title, lines


def block_list_names(blocks):
    b_items = []
    for b in blocks:
        if hasattr(b, "name"):
            b_items.append(b.name)
        else:  # indexed
            for i in b:
                b_items.append(i.name)
    return b_items


def _plural(n, word):
    s = "s" if abs(n) > 1 else ""
    return f"{n} {word}{s}"


class VCSet(BaseModel):
    variables: list[str]
    constraints: list[str]

    @classmethod
    def from_blocks(cls, var_blocks, const_blocks) -> "VCSet":
        v_items = block_list_names(var_blocks)
        c_items = block_list_names(const_blocks)
        return VCSet(variables=v_items, constraints=c_items)


class EvalErrorData(BaseModel):
    component_name: str
    message: str


class StructuralWarningsData(BaseModel, Printable):
    dof: int | None = None
    inconsistent_units: list[str] | None = None
    underconstrained_set: VCSet | None = None
    overconstrained_set: VCSet | None = None
    evaluation_errors: list[EvalErrorData] | None = None

    def get_lines(self):
        lines = []
        if self.dof is not None:
            lines.append(f"WARNING: {_plural(self.dof, 'Degree')} of Freedom")
        if self.inconsistent_units is not None:
            lines.append(
                f"WARNING: {_plural(len(self.inconsistent_units), 'Component')} with inconsistent units"
            )
        if (
            self.underconstrained_set is not None
            or self.overconstrained_set is not None
        ):
            t = " " * 4
            ucv, ucc = len(self.underconstrained_set.variables), len(
                self.underconstrained_set.constraints
            )
            ocv, occ = len(self.overconstrained_set.variables), len(
                self.overconstrained_set.constraints
            )
            lines.extend(
                [
                    f"WARNING: Structural singularity found",
                    f"{t}Under-Constrained Set: {ucv} variables, {ucc} constraints",
                    f"{t}Over-Constrained Set: {ocv} variables, {occ} constraints",
                ]
            )
        if self.evaluation_errors is not None:
            lines.append(
                f"WARNING: Found {len(self.evaluation_errors)} potential evaluation errors."
            )

        return "Structural warnings", lines


class StructuralCautionsData(BaseModel, Printable):
    zero_vars: list[str] | None = None
    unused_vars_free: list[str] | None = None
    unused_vars_fixed: list[str] | None = None

    def get_lines(self):
        lines = []
        if self.zero_vars is not None:
            lines.append(
                f"Caution: {_plural(len(self.zero_vars), 'variable')} fixed to 0"
            )
        if self.unused_vars_free is not None or self.unused_vars_fixed is not None:
            nfree = 0 if self.unused_vars_free is None else len(self.unused_vars_free)
            nfixed = (
                0 if self.unused_vars_fixed is None else len(self.unused_vars_fixed)
            )
            lines.append(
                f"Caution: {_plural(nfree + nfixed, 'unused variable')} ({nfixed} fixed)"
            )
        return "Structural cautions", lines


class StructuralIssuesData(BaseModel, Printable):
    """Structural issues (warnings, cautions).

    All possibilities are listed, the value will be None if it is not
    an issue for this model.
    """

    warnings: StructuralWarningsData
    cautions: StructuralCautionsData

    def get_lines(self):
        wt, wlines = self.warnings.get_lines()
        ct, clines = self.cautions.get_lines()
        return (
            "Structural issues",
            [wt, "-" * len(wt)] + wlines + ["", ct, "-" * len(ct)] + clines,
        )


class ReportSectionData(BaseModel):
    """Data for a standard diagnosstructuraltics report section."""

    title: str | None = None
    lines_list: list[str] = Field(default_factory=list)
    line_if_empty: str | None = None
    end_line: str | None = None
    header: str = "-"
    footer: str | None = None

    def display(self, stream=None, max_str_length=84, tab=" " * 4):
        if stream is None:
            stream = sys.stdout

        stream.write(f"{self.header * max_str_length}\n")
        if self.title is not None:
            stream.write(f"{self.title}\n\n")
        if len(self.lines_list) > 0:
            for i in self.lines_list:
                stream.write(f"{tab}{i}\n")
        elif self.line_if_empty is not None:
            stream.write(f"{tab}{self.line_if_empty}\n")
        stream.write("\n")
        if self.end_line is not None:
            stream.write(f"{self.end_line}\n")
        if self.footer is not None:
            stream.write(f"{self.footer * max_str_length}\n")


class TextListData(BaseModel):
    """Simple container for a list of text entries."""

    entries: list[str] = Field(default_factory=list)


class NamedValueData(BaseModel):
    """Name/value diagnostics entry with optional details."""

    name: str
    value: float | int | str
    details: str | None = None
    value_format: str | None = None

    def format_line(self, include_value_label: bool = False):
        if self.value_format is None:
            value_str = str(self.value)
        else:
            value_str = format(self.value, self.value_format)

        if include_value_label:
            line = f"{self.name}: value={value_str}"
        else:
            line = f"{self.name}: {value_str}"

        if self.details is not None:
            line += f" {self.details}"
        return line


class NamedValueListData(BaseModel):
    """Container for diagnostics entries with named values."""

    entries: list[NamedValueData] = Field(default_factory=list)


class NamedPairData(BaseModel):
    """Container for a pair of named entities."""

    first: str
    second: str

    def format_line(self):
        return f"{self.first}, {self.second}"


class NamedPairListData(BaseModel):
    """Container for paired diagnostics entries."""

    entries: list[NamedPairData] = Field(default_factory=list)


class ProblematicConstraintTermsData(BaseModel):
    """Detailed result for problematic term analysis in a constraint."""

    constraint_name: str
    issues: list[str] = Field(default_factory=list)
    end_line: str | None = None


class DulmageMendelsohnBlockData(BaseModel):
    """Variable and constraint names for an independent DM block."""

    block_index: int
    variables: list[str] = Field(default_factory=list)
    constraints: list[str] = Field(default_factory=list)


class DulmageMendelsohnPartitionData(BaseModel):
    """Display data for a Dulmage-Mendelsohn partition report."""

    title: str
    blocks: list[DulmageMendelsohnBlockData] = Field(default_factory=list)
    header: str = "="
    footer: str = "="


# class StructuralIssuesData(BaseModel):
#     """Computed data for the structural issues report."""

#     statistics: list[str] = Field(default_factory=list)
#     warnings: list[str] = Field(default_factory=list)
#     cautions: list[str] = Field(default_factory=list)
#     next_steps: list[str] = Field(default_factory=list)


class NumericalIssuesData(BaseModel):
    """Computed data for the numerical issues report."""

    statistics: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    cautions: list[str] = Field(default_factory=list)
    next_steps: list[str] = Field(default_factory=list)
