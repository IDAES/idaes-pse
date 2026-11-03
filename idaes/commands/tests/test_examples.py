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
Tests for idaes.commands.examples

This command has been removed and now only reports where to find the examples
package; the tests confirm that reporting behavior.
"""

from click.testing import CliRunner
import pytest

from idaes.commands import examples


@pytest.fixture
def runner():
    return CliRunner()


@pytest.mark.unit
def test_get_examples_removed(runner):
    result = runner.invoke(examples.get_examples)
    assert result.exit_code == 0
    assert "idaes-examples" in result.output
