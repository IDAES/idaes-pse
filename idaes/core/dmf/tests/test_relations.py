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
Test the 'relations' connecting Resource instances,
which spans the modules idaes.core.dmf.{dmf, resource, resourcedb}
"""
# stdlib
import logging
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Union

# third-party
import pytest

# local
from idaes.core.dmf import experiment, resource, DMF
from idaes.core.dmf.resource import Predicates


# for testing
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
def test_create_relation_in_resource():
    a = resource.Resource()
    b = resource.Resource()
    resource.create_relation(a, "contains", b)
    assert len(a.v["relations"]) == 1
    assert len(b.v["relations"]) == 1
    # bad type
    with pytest.raises(TypeError):
        resource.create_relation("foo", "contains", b)


@pytest.mark.unit
def test_relation_in_experiment():
    tmp_dir = scratch_path / "relation_in_experiment"
    dmf = DMF(path=tmp_dir, create=True)
    e1 = experiment.Experiment(dmf, name="1")
    a = resource.Resource(value={"name": "foo"})
    e1.add(a)
    assert len(a.v["relations"]) == 1
    assert len(e1.v["relations"]) == 1


@pytest.mark.unit
def test_relation_with_remove():
    tmp_dir = scratch_path / "relation_with_remove"
    dmf = DMF(path=tmp_dir, create=True)
    e1 = experiment.Experiment(dmf, name="1")
    n, added = 10, []
    for i in range(n):
        a = resource.Resource({"name": "foo"})
        e1.add(a)
        added.append(a)
    assert len(e1.v["relations"]) == n
    # remove, then update e1
    for a in added:
        dmf.remove(identifier=a.id)
        e1.update()
        # relation to removed 'a' should be gone
        n -= 1
        assert (len(e1.v["relations"])) == n


@pytest.mark.unit
def test_find_related():
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
    tmp_dir = scratch_path / "find_related"
    dmf = DMF(path=tmp_dir, create=True)
    r = [resource.Resource({"name": "r{}".format(i)}) for i in range(5)]
    # r3 <-- derived <-- r2 <-- version <-- r1
    cr = resource.create_relation  # shortcut
    cr(r[2], Predicates.derived, r[3])
    cr(r[1], Predicates.version, r[2])
    # r4 <-- derived <-- r2
    cr(r[2], Predicates.derived, r[4])
    # r0 -- Uses --> r1
    cr(r[0], Predicates.uses, r[1])
    # add to dmf
    for i in range(5):
        dmf.add(r[i])
    # outgoing from r0 should include 1,2,3,4
    names = []
    for d, rr, m in dmf.find_related(r[0], meta=["aliases"]):
        names.append(m["aliases"][0])
    names.sort()
    assert names == ["r1", "r2", "r3", "r4"]
    # incoming to r4 should include r0, r1, r2
    names = []
    for d, rr, m in dmf.find_related(r[4], meta=["aliases"], outgoing=False):
        names.append(m["aliases"][0])
    names.sort()
    assert names == ["r0", "r1", "r2"]


@pytest.mark.unit
def test_circular():
    #
    # r0 -> derived -> r1 -> derived >- r2 -+
    #  ^                                    |
    #  +------------------------------------+
    #     uses
    tmp_dir = scratch_path / "circular"
    dmf = DMF(path=tmp_dir, create=True)
    r = [resource.Resource({"name": "r{}".format(i)}) for i in range(3)]
    resource.create_relation(r[0], Predicates.derived, r[1])
    resource.create_relation(r[1], Predicates.derived, r[2])
    resource.create_relation(r[2], Predicates.uses, r[0])
    for rr in r:
        dmf.add(rr)
    # outgoing from r0
    names = []
    for d, rr, m in dmf.find_related(r[0], meta=["aliases"]):
        names.append(m["aliases"][0])
    names.sort()
    assert names == ["r0", "r1", "r2"]
    # incoming to r1
    names = []
    for d, rr, m in dmf.find_related(r[0], meta=["aliases"], outgoing=False):
        names.append(m["aliases"][0])
    names.sort()
    assert names == ["r0", "r1", "r2"]
    # reducing depth shortens output
    names = []
    for d, rr, m in dmf.find_related(r[0], meta=["aliases"], maxdepth=2):
        names.append(m["aliases"][0])
    names.sort()
    assert names == ["r1", "r2"]
    names = []
    for d, rr, m in dmf.find_related(r[0], meta=["aliases"], maxdepth=1):
        names.append(m["aliases"][0])
    names.sort()
    assert names == ["r1"]
