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


def test_init_args(runner):
    with runner.isolated_filesystem():
        result = runner.invoke(
            cli.init, ["--path", "ws", "--create", "--name", "x", "--desc", "X"]
        )
        assert result.exit_code == 0
        # second time, no --create required
        result = runner.invoke(cli.init, ["--path", "ws"])
        assert result.exit_code == 0


def test_init_prompt(runner):
    with runner.isolated_filesystem():
        result = runner.invoke(cli.init, ["--path", "ws", "--create"], input="x\nX\n")
        assert result.exit_code == 0
        # second time, no --create required
        result = runner.invoke(cli.init, ["--path", "ws"])
        assert result.exit_code == 0
