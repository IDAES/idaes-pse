##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Test the `idaes.dmf.cli` module.
"""
# stdlib
import logging
from urllib.parse import urlparse
import sys

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


@pytest.mark.unit
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


@pytest.mark.unit
def test_aliases():
    ag = cli.AliasedGroup(
        aliases={
            "clark": "superman",
            "bruce": "batman",
            "carter": "hawkman",
            "barbara": "batgirl",
        },
        commands={"superman": None, "batman": None, "hawkman": None, "batgirl": None},
    )
    context = MockContext("superman")
    cmd = ag.get_command(context, "clark")
    # match alias and command
    pytest.raises(ContextFailed, ag.get_command, context, "b")
    # match multiple aliases
    pytest.raises(ContextFailed, ag.get_command, context, "c")
    # match multiple commands (superheroes)
    pytest.raises(ContextFailed, ag.get_command, context, "bat")


@pytest.mark.unit
def test_url():
    ut = cli.URLType()
    u = "http://other.org"
    result = ut.convert(u, None, None)
    assert result == urlparse(u)
    assert ut.convert(urlparse(u), None, None) == urlparse(u)
