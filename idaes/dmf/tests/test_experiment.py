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
Tests for idaes.dmf.experiment module
"""
import logging
import pytest
import sys

from idaes.dmf import experiment, errors, resource
from .util import tmp_dmf  # noqa -- this is a fixture
from .util import init_logging, tmp_dmf

__author__ = "Dan Gunter <dkgunter@lbl.gov>"

if sys.platform.startswith("win"):
    pytest.skip("skipping DMF tests on Windows", allow_module_level=True)

init_logging()
_log = logging.getLogger(__name__)


@pytest.mark.unit
def test_init(tmp_dmf):
    exp = experiment.Experiment(tmp_dmf, name="try1", desc="Nice try")
    assert exp.v["name"] == "try1"
    assert exp.id


@pytest.mark.component
def test_create_many(tmp_dmf):
    for i in range(100):
        exp = experiment.Experiment(tmp_dmf, name="try{i}".format(i=i))
        assert exp.v["name"] == "try{i}".format(i=i)
        assert exp.id


@pytest.mark.unit
def test_remove_workflow(tmp_dmf):
    # A workflow of copy/remove
    # make an experiment
    e1 = experiment.Experiment(tmp_dmf, name="one", version="0.0.1")
    # make a new version of it
    e2 = e1.copy(version="0.0.2")
    # link the two together
    e1.link(e2, predicate=resource.PR_VERSION)
    # remove the first one (what happens to the link?)
    e1.remove()
    # check that the first one can't be used any more
    with pytest.raises(errors.BadResourceError):
        e1.update()
    with pytest.raises(errors.BadResourceError):
        e1.link(e2, predicate=resource.PR_VERSION)
    # check that the copy can still be modified
    e2.v["desc"] = "This is a copy of e1"
    e2.update()  # this fixes relations in the DB
    # check that the link is gone, i.e. there are no
    # relations in e2 any more
    assert len(e2.v["relations"]) == 0

