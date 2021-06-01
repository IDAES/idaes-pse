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
Tests for idaes.dmf.propindex module.
"""
# stdlib
import logging
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Union

# third-party
import pytest

# package
import idaes
from idaes.dmf import propindex
from idaes.dmf import DMF
from idaes.dmf import resource

# for testing
from .util import init_logging
from . import for_propindex

__author__ = "Dan Gunter"

init_logging()
_log = logging.getLogger(__name__)

scratch_dir: Union[str, None] = None
scratch_path: Union[Path, None] = None


def setup_module(module):
    global scratch_dir, scratch_path
    scratch_dir = TemporaryDirectory(prefix="idaes_dmf_")  # easier to remove later
    scratch_path = Path(scratch_dir.name)


def teardown_module(module):
    global scratch_dir
    del scratch_dir


@pytest.mark.unit
def test_dmfvisitor():
    with pytest.raises(TypeError):
        propindex.DMFVisitor("o_O", "a")
    propindex.DMFVisitor("o_O")


@pytest.mark.unit
def test_index_property_metadata():
    tmp_dir = scratch_path / "index_property_metadata"
    dmf = DMF(path=tmp_dir, create=True)
    propindex.index_property_metadata(
        dmf, pkg=idaes.dmf, expr=".*IndexMePlease[0-9]", exclude_testdirs=False
    )
    # Check the resource
    for rsrc in dmf.find():
        assert rsrc.v[rsrc.TYPE_FIELD] == resource.ResourceTypes.code
        # print('@@ GOT RESOURCE:\n{}'.format(rsrc.v))


@pytest.mark.unit
def test_index_multiple_versions():
    tmp_dir = scratch_path / "index_multiple_versions"
    dmf = DMF(path=tmp_dir, create=True)

    v1, v2, v3 = "1.0.0", "6.6.6", "9.9.0"
    # index initial version
    propindex.index_property_metadata(
        dmf,
        pkg=idaes.dmf,
        expr=".*IndexMePlease[0-9]",
        exclude_testdirs=False,
        default_version=v1,
    )
    # index again
    propindex.index_property_metadata(
        dmf,
        pkg=idaes.dmf,
        expr=".*IndexMePlease[0-9]",
        exclude_testdirs=False,
        default_version=v2,
    )
    # check that we now have two resources, and
    # a relation between them
    rlist = list(dmf.find({}))
    assert len(rlist) == 2
    rcodes = [r.v["codes"][0] for r in rlist]
    if rcodes[0]["version"][:3] == ("6", "6", "6"):
        first, second = 1, 0
    else:
        first, second = 0, 1

    # Debugging only
    # print('CODES:')
    # print(' - first -')
    # print(rcodes[first])
    # print(rlist[first].v['relations'])
    # print(' - second -')
    # print(rcodes[second])
    # print(rlist[second].v['relations'])

    # Each resource has 1 relation
    assert len(rlist[first].v["relations"]) == 1
    assert len(rlist[second].v["relations"]) == 1
    first_rel = rlist[first].v["relations"][0]
    second_rel = rlist[second].v["relations"][0]
    # First resource is pointed at by second
    assert first_rel[resource.RR_ROLE] == resource.RR_OBJ
    assert first_rel[resource.RR_PRED] == resource.Predicates.version
    assert first_rel[resource.RR_ID] == rlist[second].id
    # Second resource points at first
    assert second_rel[resource.RR_ROLE] == resource.RR_SUBJ
    assert second_rel[resource.RR_PRED] == resource.Predicates.version
    assert second_rel[resource.RR_ID] == rlist[first].id

    # Add the same version
    propindex.index_property_metadata(
        dmf,
        pkg=idaes.dmf,
        expr=".*IndexMePlease[0-9]",
        exclude_testdirs=False,
        default_version=v2,
    )
    # check that we still have two resources
    rlist = list(dmf.find({}))
    assert len(rlist) == 2

    # Now add another version
    propindex.index_property_metadata(
        dmf,
        pkg=idaes.dmf,
        expr=".*IndexMePlease[0-9]",
        exclude_testdirs=False,
        default_version=v3,
    )
    # check that we now have three resources
    rlist = list(dmf.find({}))
    assert len(rlist) == 3
    # check that we have 0 <--> 1 <--> 2
    # first sort by version and save that in the 'indexes' array
    indexes = [(r.v["codes"][0]["version"], i) for i, r in enumerate(rlist)]
    indexes.sort()
    # pull out relations into 'rel' array, in version order
    rel = [rlist[indexes[i][1]].v["relations"] for i in range(3)]
    # check first resource's relations
    assert len(rel[0]) == 1
    # 0 <-- 1
    assert rel[0][0][resource.RR_ID] == rlist[indexes[1][1]].id
    assert rel[0][0][resource.RR_ROLE] == resource.RR_OBJ
    # check second resource's relations
    assert len(rel[1]) == 2
    for j in range(2):
        if rel[1][j][resource.RR_ROLE] == resource.RR_SUBJ:
            # 1 --> 0
            assert rel[1][j][resource.RR_ID] == rlist[indexes[0][1]].id
        else:
            # 1 <-- 2
            assert rel[1][j][resource.RR_ID] == rlist[indexes[2][1]].id
            assert rel[1][j][resource.RR_ROLE] == resource.RR_OBJ
    # check third resource's relations
    # check third resource's relations
    assert len(rel[2]) == 1
    # 2 --> 1
    assert rel[2][0][resource.RR_ID] == rlist[indexes[1][1]].id
    assert rel[2][0][resource.RR_ROLE] == resource.RR_SUBJ

