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
Tests for idaes.commands.env_info
"""

import json
import os
import tempfile

from click.testing import CliRunner
import pytest

from idaes.commands import env_info


@pytest.fixture
def runner():
    return CliRunner()


@pytest.mark.unit
def test_environment_info(runner):
    result = runner.invoke(env_info.environment_info)
    assert result.exit_code == 0
    assert result.output.strip() != ""


@pytest.mark.unit
def test_environment_info_extra_solver(runner):
    result = runner.invoke(env_info.environment_info, ["--solver", "ipopt"])
    assert result.exit_code == 0


@pytest.mark.unit
def test_environment_info_json(runner):
    with tempfile.TemporaryDirectory() as tmpdir:
        json_path = os.path.join(tmpdir, "env.json")
        result = runner.invoke(env_info.environment_info, ["--json", json_path])
        assert result.exit_code == 0
        assert os.path.isfile(json_path)
        with open(json_path, "r") as f:
            data = json.load(f)
        assert isinstance(data, dict)
