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
This module contains utility functions for reporting model diagnostics.
"""

__author__ = "Alexander Dowling, Douglas Allan, Andrew Lee, Robby Parker, Ben Knueven"


MAX_STR_LENGTH = 84
TAB = " " * 4


def write_report_section(
    stream,
    lines_list,
    title=None,
    line_if_empty=None,
    end_line=None,
    header="-",
    footer=None,
):
    """
    Writes output in standard format for report and display methods.

    Args:
        stream: stream to write to
        lines_list: list containing lines to be written in body of report
        title: title to be put at top of report
        line_if_empty: line to be written if lines_list is empty
        end_line: line to be written at end of report
        header: character to use to write header separation line
        footer: character to use to write footer separation line

    Returns:
        None

    """
    stream.write(f"{header * MAX_STR_LENGTH}\n")
    if title is not None:
        stream.write(f"{title}\n\n")
    if len(lines_list) > 0:
        for i in lines_list:
            stream.write(f"{TAB}{i}\n")
    elif line_if_empty is not None:
        stream.write(f"{TAB}{line_if_empty}\n")
    stream.write("\n")
    if end_line is not None:
        stream.write(f"{end_line}\n")
    if footer is not None:
        stream.write(f"{footer * MAX_STR_LENGTH}\n")
