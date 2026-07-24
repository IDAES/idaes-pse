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
Tests for idaes.commands.convergence

These commands have been deprecated and now only report that they no longer
work; the tests confirm that reporting behavior.
"""

from click.testing import CliRunner
import pytest

from idaes.commands import convergence


@pytest.fixture
def runner():
    return CliRunner()


@pytest.mark.unit
def test_convergence_sample_deprecated(runner):
    result = runner.invoke(
        convergence.convergence_sample,
        ["-e", "SomeClass", "-s", "sample.json", "-N", "5"],
    )
    assert result.exit_code == 0
    assert "deprecated" in result.output


@pytest.mark.unit
def test_convergence_eval_deprecated(runner):
    result = runner.invoke(convergence.convergence_eval, ["-s", "sample.json"])
    assert result.exit_code == 0
    assert "deprecated" in result.output


@pytest.mark.unit
def test_convergence_search_deprecated(runner):
    result = runner.invoke(convergence.convergence_search)
    assert result.exit_code == 0
    assert "deprecated" in result.output
