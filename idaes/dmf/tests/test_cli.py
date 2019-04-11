"""
Test the `idaes.dmf.cli` module.
"""
# stdlib
import os

# third-party
from click.testing import CliRunner
import pytest

# package
from idaes.dmf import cli

__author__ = "Dan Gunter"


@pytest.fixture
def runner():
    # trivial, but allows later flexibility
    return CliRunner()

# Note: at some point, we may put some more involved
# tests in here, but for now all the tests are in the
# dmf-cli.rst file in testcode directives.