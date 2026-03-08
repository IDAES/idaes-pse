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


class ReportSectionData(BaseModel):
    """Data for a standard diagnostics report section."""

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


class StructuralIssuesData(BaseModel):
    """Computed data for the structural issues report."""

    statistics: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    cautions: list[str] = Field(default_factory=list)
    next_steps: list[str] = Field(default_factory=list)


class NumericalIssuesData(BaseModel):
    """Computed data for the numerical issues report."""

    statistics: list[str] = Field(default_factory=list)
    warnings: list[str] = Field(default_factory=list)
    cautions: list[str] = Field(default_factory=list)
    next_steps: list[str] = Field(default_factory=list)
