#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Tests for DMF 'workspace' module.
"""
# stdlib
import logging
import os

# third-party
import pytest

# local
from idaes.core.dmf import workspace, errors
from idaes.core.dmf.util import TempDir

# for testing
from .util import init_logging

__author__ = "Dan Gunter"

init_logging()
_log = logging.getLogger(__name__)

BADPATH = os.path.join("should", "never", "exist", "seriously", "no")


@pytest.fixture
def wsdir():
    with TempDir() as d:
        yield d
        print("Removing temporary workspace directory: {}".format(d))


# Tests
# -----


@pytest.mark.unit
def test_ws_init_notfound():
    with pytest.raises(errors.WorkspaceNotFoundError):
        workspace.Workspace(BADPATH)


@pytest.mark.unit
def test_ws_init_noconf(wsdir):
    with pytest.raises(errors.WorkspaceConfNotFoundError):
        workspace.Workspace(wsdir)


@pytest.mark.unit
def test_ws_init_badconf(wsdir):
    conf = open(os.path.join(wsdir, workspace.Workspace.WORKSPACE_CONFIG), "w")
    conf.write("note: this config is wack\n")
    conf.close()
    with pytest.raises(errors.WorkspaceConfMissingField):
        workspace.Workspace(wsdir)


@pytest.mark.unit
def test_ws_init_create_nodefaults(wsdir):
    ws = workspace.Workspace(wsdir, create=True)
    # if defaults were set, then HTML documentation paths
    # will not be empty, so test this value
    assert ws.get_doc_paths() == []


@pytest.mark.unit
def test_ws_init_create_defaults(wsdir):
    ws = workspace.Workspace(wsdir, create=True, add_defaults=True)
    # if defaults were set, then HTML documentation paths
    # will not be empty, so test this value
    assert ws.get_doc_paths() != []


@pytest.mark.unit
def test_ws_accessors(wsdir):
    ws = workspace.Workspace(wsdir, create=True, add_defaults=True)
    assert ws.wsid != ""
    assert os.path.normpath(ws.root) == os.path.normpath(wsdir)
    _ = ws.name
    _ = ws.description
