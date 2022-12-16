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
Tests for idaes.core.dmf.experiment module
"""
import logging
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Union

import pytest

from idaes.core.dmf import experiment, errors, DMF
from idaes.core.dmf.resource import Predicates
from .util import init_logging

__author__ = "Dan Gunter"

init_logging()
_log = logging.getLogger(__name__)


scratch_dir: Union[str, None] = None
scratch_path: Union[Path, None] = None


def setup_module(module):
    global scratch_dir, scratch_path
    scratch_dir = TemporaryDirectory(prefix="idaes.core.dmf_")  # easier to remove later
    scratch_path = Path(scratch_dir.name)


def teardown_module(module):
    global scratch_dir
    del scratch_dir


@pytest.mark.unit
def test_init():
    tmp_dir = scratch_path / "init"
    dmf = DMF(path=tmp_dir, create=True)
    exp = experiment.Experiment(dmf, name="try1", desc="Nice try")
    assert exp.name == "try1"
    assert exp.id


@pytest.mark.component
def test_create_many():
    tmp_dir = scratch_path / "create_many"
    dmf = DMF(path=tmp_dir, create=True)
    for i in range(100):
        exp = experiment.Experiment(dmf, name="try{i}".format(i=i))
        assert exp.name == "try{i}".format(i=i)
        assert exp.id


@pytest.mark.unit
def test_remove_workflow():
    tmp_dir = scratch_path / "remove_workflow"
    dmf = DMF(path=tmp_dir, create=True)
    # A workflow of copy/remove
    # make an experiment
    e1 = experiment.Experiment(dmf, name="one", version="0.0.1")
    # make a new version of it
    e2 = e1.copy(version="0.0.2")
    # link the two together
    e1.link(e2, predicate=Predicates.version)
    # remove the first one (what happens to the link?)
    e1.remove()
    # check that the first one can't be used any more
    with pytest.raises(errors.BadResourceError):
        e1.update()
    with pytest.raises(errors.BadResourceError):
        e1.link(e2, predicate=Predicates.version)
    # check that the copy can still be modified
    e2.v["desc"] = "This is a copy of e1"
    e2.update()  # this fixes relations in the DB
    # check that the link is gone, i.e. there are no
    # relations in e2 any more
    assert len(e2.v["relations"]) == 0
