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
Tests for idaes.commands.data_directory
"""

from click.testing import CliRunner
import pytest

import idaes
from idaes.commands import data_directory


@pytest.fixture
def runner():
    return CliRunner()


@pytest.mark.unit
def test_data_directory_show(runner):
    result = runner.invoke(data_directory.data_directory)
    assert result.exit_code == 0
    assert idaes.data_directory in result.output


@pytest.mark.unit
def test_data_directory_exists(runner):
    result = runner.invoke(data_directory.data_directory, ["--exists"])
    assert result.exit_code == 0
    assert result.output.strip() in ("True", "False")


@pytest.mark.unit
def test_data_directory_create(runner):
    result = runner.invoke(data_directory.data_directory, ["--create"])
    assert result.exit_code == 0
    assert "Creating" in result.output


@pytest.mark.unit
def test_bin_directory_show(runner):
    result = runner.invoke(data_directory.bin_directory)
    assert result.exit_code == 0
    assert idaes.bin_directory in result.output


@pytest.mark.unit
def test_bin_directory_exists(runner):
    result = runner.invoke(data_directory.bin_directory, ["--exists"])
    assert result.exit_code == 0
    assert result.output.strip() in ("True", "False")


@pytest.mark.unit
def test_bin_directory_create(runner):
    result = runner.invoke(data_directory.bin_directory, ["--create"])
    assert result.exit_code == 0
    assert "Creating" in result.output
