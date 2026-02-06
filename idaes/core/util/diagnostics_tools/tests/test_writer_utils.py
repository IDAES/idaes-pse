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
This module contains tests for diagnostics writer utils.
"""
from io import StringIO

import pytest


from idaes.core.util.diagnostics_tools.writer_utils import (
    write_report_section,
)

__author__ = "Alex Dowling, Douglas Allan, Andrew Lee"


@pytest.mark.unit
def test_write_report_section_all():
    stream = StringIO()

    write_report_section(
        stream=stream,
        lines_list=["a", "b", "c"],
        title="foo",
        line_if_empty="bar",
        end_line="baz",
        header="-",
        footer="=",
    )

    expected = """------------------------------------------------------------------------------------
foo

    a
    b
    c

baz
====================================================================================
"""
    assert stream.getvalue() == expected


@pytest.mark.unit
def test_write_report_section_no_lines():
    stream = StringIO()

    write_report_section(
        stream=stream,
        lines_list=[],
        title="foo",
        line_if_empty="bar",
        end_line="baz",
        header="-",
        footer="=",
    )

    expected = """------------------------------------------------------------------------------------
foo

    bar

baz
====================================================================================
"""
    assert stream.getvalue() == expected


@pytest.mark.unit
def test_write_report_section_lines_only():
    stream = StringIO()

    write_report_section(
        stream=stream,
        lines_list=["a", "b", "c"],
    )

    expected = """------------------------------------------------------------------------------------
    a
    b
    c

"""
    assert stream.getvalue() == expected
