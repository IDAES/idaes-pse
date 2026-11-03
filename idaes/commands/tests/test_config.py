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
Tests for idaes.commands.config

These commands have been removed and now only report that removal; the tests
confirm that reporting behavior.
"""

from click.testing import CliRunner
import pytest

from idaes.commands import config


@pytest.fixture
def runner():
    return CliRunner()


@pytest.mark.unit
@pytest.mark.parametrize(
    "command",
    [
        "config_write",
        "config_file",
        "config_set",
        "config_display",
    ],
)
def test_config_command_removed(runner, command):
    result = runner.invoke(getattr(config, command))
    assert result.exit_code == 0
    assert "removed" in result.output


@pytest.mark.unit
def test_config_command_ignores_extra_args(runner):
    result = runner.invoke(config.config_set, ["some_key", "some_value"])
    assert result.exit_code == 0
    assert "removed" in result.output
