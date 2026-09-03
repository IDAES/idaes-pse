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
Tests for idaes.commands.base
"""

import logging

from click.testing import CliRunner
import pytest

from idaes.commands import base


@pytest.fixture
def runner():
    return CliRunner()


@pytest.mark.unit
@pytest.mark.parametrize(
    "verbosity,expected",
    [
        (3, logging.DEBUG),
        (2, logging.INFO),
        (1, logging.WARN),
        (0, logging.ERROR),
        (-1, logging.FATAL),
        (-2, logging.FATAL + 1),
    ],
)
def test_level_from_verbosity(verbosity, expected):
    assert base.level_from_verbosity(verbosity) == expected


@pytest.mark.unit
def test_how_to_report_an_error_with_bars():
    msg = base.how_to_report_an_error()
    assert "github.com/idaes/idaes-pse/issues" in msg
    assert msg.startswith("-")


@pytest.mark.unit
def test_how_to_report_an_error_embedded():
    msg = base.how_to_report_an_error(embed=True)
    assert "github.com/idaes/idaes-pse/issues" in msg
    assert not msg.startswith("-")


@pytest.mark.unit
def test_command_base_verbose(runner):
    result = runner.invoke(base.command_base, ["--verbose", "copyright"])
    assert result.exit_code == 0


@pytest.mark.unit
def test_command_base_quiet(runner):
    result = runner.invoke(base.command_base, ["--quiet", "copyright"])
    assert result.exit_code == 0


@pytest.mark.unit
def test_command_base_verbose_quiet_conflict(runner):
    result = runner.invoke(base.command_base, ["--verbose", "--quiet", "copyright"])
    assert result.exit_code != 0


@pytest.mark.unit
def test_copyright(runner):
    result = runner.invoke(base.command_base, ["copyright"])
    assert result.exit_code == 0
    assert "IDAES" in result.output


@pytest.mark.unit
def test_import_time(runner):
    result = runner.invoke(base.command_base, ["import-time"])
    assert result.exit_code == 0
    assert "Time:" in result.output
