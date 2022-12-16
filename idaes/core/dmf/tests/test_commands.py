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
Tests for idaes.core.dmf.commands
"""
import logging
import os
import shutil

#
import pytest

#
from idaes.core.dmf import dmfbase, commands, errors
from idaes.core.dmf.util import mkdtemp
from .util import init_logging

__author__ = "Dan Gunter"

init_logging()
_log = logging.getLogger(__name__)


@pytest.fixture(scope="function")
def wspath():
    dirname = mkdtemp()
    yield dirname
    # teardown
    shutil.rmtree(dirname)


@pytest.mark.unit
def test_workspace_init(wspath):
    commands.workspace_init(wspath, {"some": "metadata"})
    # try again. Should work, since it's OK to init twice
    commands.workspace_init(wspath, {"some": "metadata"})


@pytest.mark.unit
def test_workspace_info(wspath):
    commands.workspace_init(wspath, {"some": "metadata"})
    commands.workspace_info(wspath)

    #    subdir = os.path.join(wspath, 'stuff')
    #    os.mkdir(subdir)
    #    commands.workspace_info(subdir)

    notasubdir = os.path.join(wspath, "nope")
    try:
        commands.workspace_info(notasubdir)
        assert False, "Nonexistent subdir workspace info success"
    except errors.CommandError:
        pass


@pytest.mark.unit
def test_find_html_docs(wspath):
    filedir = os.path.dirname(__file__)
    docpath = os.path.join(filedir, "..", "docs", "build", "html")
    if os.path.exists(docpath):
        commands.workspace_init(wspath, {}, html_paths=[docpath])
        dmfobj = dmfbase.DMF(wspath)
        filenames = commands.find_html_docs(dmfobj, dmfobj)
        assert len(filenames) > 0
