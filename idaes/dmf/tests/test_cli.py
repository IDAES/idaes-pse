"""
Test the `idaes.dmf.cli` module.
"""
# stdlib
import logging
from urllib.parse import urlparse

# third-party
import click
from click.testing import CliRunner
import pytest

# package
from idaes.dmf import cli

__author__ = "Dan Gunter"


@pytest.fixture
def runner():
    # trivial, but allows later flexibility
    return CliRunner()


def test_verbosity():
    assert cli.level_from_verbosity(10) == logging.DEBUG
    assert cli.level_from_verbosity(2) == logging.INFO
    assert cli.level_from_verbosity(1) == logging.WARNING
    assert cli.level_from_verbosity(0) == logging.ERROR
    assert cli.level_from_verbosity(-1) == logging.FATAL
    assert cli.level_from_verbosity(-2) > logging.FATAL


class ContextFailed(Exception):
    pass


class MockContext:
    def __init__(self, command, parent=None):
        self.command = command
        self.parent = parent

    def fail(self, msg):
        raise ContextFailed()


def test_aliases():
    ag = cli.AliasedGroup(aliases={"clark": "superman", "bruce": "batman",
                                   "carter": "hawkman", "barbara": "batgirl"},
                          commands={"superman": None, "batman": None,
                                    "hawkman": None, "batgirl": None})
    context = MockContext("superman")
    cmd = ag.get_command(context, "clark")
    # match alias and command
    pytest.raises(ContextFailed, ag.get_command, context, "b")
    # match multiple aliases
    pytest.raises(ContextFailed, ag.get_command, context, "c")
    # match multiple commands (superheroes)
    pytest.raises(ContextFailed, ag.get_command, context, "bat")


def test_url():
    ut = cli.URLType()
    u = "http://other.org"
    result = ut.convert(u, None, None)
    assert result == urlparse(u)
    assert ut.convert(urlparse(u), None, None) == urlparse(u)
