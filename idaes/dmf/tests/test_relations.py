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
Test the 'relations' connecting Resource instances,
which spans the modules idaes.dmf.{dmf, resource, resourcedb}
"""
# stdlib
import logging
import sys

# third-party
import pytest

# local
from idaes.dmf import experiment, resource

# for testing
from .util import tmp_dmf  # need for fixture
from .util import init_logging

__author__ = "Dan Gunter <dkgunter@lbl.gov>"

if sys.platform.startswith("win"):
    pytest.skip("skipping DMF tests on Windows", allow_module_level=True)

init_logging()
_log = logging.getLogger(__name__)


@pytest.mark.unit
def test_create_relation_in_resource():
    a = resource.Resource()
    b = resource.Resource()
    resource.create_relation_args(a, "contains", b)
    assert len(a.v["relations"]) == 1
    assert len(b.v["relations"]) == 1
    # bad type
    with pytest.raises(TypeError):
        resource.create_relation("foo", "contains", b)


@pytest.mark.unit
def test_relation_in_experiment(tmp_dmf):
    e1 = experiment.Experiment(tmp_dmf, name="1")
    a = resource.Resource(value={"name": "foo"})
    e1.add(a)
    assert len(a.v["relations"]) == 1
    assert len(e1.v["relations"]) == 1


@pytest.mark.unit
def test_relation_with_remove(tmp_dmf):
    e1 = experiment.Experiment(tmp_dmf, name="1")
    n, added = 10, []
    for i in range(n):
        a = resource.Resource({"name": "foo"})
        e1.add(a)
        added.append(a)
    assert len(e1.v["relations"]) == n
    # remove, then update e1
    for a in added:
        tmp_dmf.remove(identifier=a.id)
        e1.update()
        # relation to removed 'a' should be gone
        n -= 1
        assert (len(e1.v["relations"])) == n


@pytest.mark.unit
def test_find_related(tmp_dmf):
    #
    #  r0
    #   | uses
    #   v
    #   r1
    #   | version
    #   v
    #   r2
    #   /\
    #  /  \ derived
    # v   v
    # r3  r4
    #
    r = [resource.Resource({"name": "r{}".format(i)}) for i in range(5)]
    # r3 <-- derived <-- r2 <-- version <-- r1
    cr = resource.create_relation_args  # shortcut
    cr(r[2], resource.PR_DERIVED, r[3])
    cr(r[1], resource.PR_VERSION, r[2])
    # r4 <-- derived <-- r2
    cr(r[2], resource.PR_DERIVED, r[4])
    # r0 -- Uses --> r1
    cr(r[0], resource.PR_USES, r[1])
    # add to dmf
    for i in range(5):
        tmp_dmf.add(r[i])
    # outgoing from r0 should include 1,2,3,4
    names = []
    for d, rr, m in tmp_dmf.find_related(r[0], meta=["name"]):
        names.append(m["name"])
    names.sort()
    assert names == ["r1", "r2", "r3", "r4"]
    # incoming to r4 should include r0, r1, r2
    names = []
    for d, rr, m in tmp_dmf.find_related(r[4], meta=["name"], outgoing=False):
        names.append(m["name"])
    names.sort()
    assert names == ["r0", "r1", "r2"]


@pytest.mark.unit
def test_circular(tmp_dmf):
    #
    # r0 -> derived -> r1 -> derived >- r2 -+
    #  ^                                    |
    #  +------------------------------------+
    #     uses
    r = [resource.Resource({"name": "r{}".format(i)}) for i in range(3)]
    resource.create_relation_args(r[0], resource.PR_DERIVED, r[1])
    resource.create_relation_args(r[1], resource.PR_DERIVED, r[2])
    resource.create_relation_args(r[2], resource.PR_USES, r[0])
    for rr in r:
        tmp_dmf.add(rr)
    # outgoing from r0
    names = []
    for d, rr, m in tmp_dmf.find_related(r[0], meta=["name"]):
        names.append(m["name"])
    names.sort()
    assert names == ["r0", "r1", "r2"]
    # incoming to r1
    names = []
    for d, rr, m in tmp_dmf.find_related(r[0], meta=["name"], outgoing=False):
        names.append(m["name"])
    names.sort()
    assert names == ["r0", "r1", "r2"]
    # reducing depth shortens output
    names = []
    for d, rr, m in tmp_dmf.find_related(r[0], meta=["name"], maxdepth=2):
        names.append(m["name"])
    names.sort()
    assert names == ["r1", "r2"]
    names = []
    for d, rr, m in tmp_dmf.find_related(r[0], meta=["name"], maxdepth=1):
        names.append(m["name"])
    names.sort()
    assert names == ["r1"]
