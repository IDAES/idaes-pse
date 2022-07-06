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
Tests for idaes.core.dmf.resource2 module
"""
# stdlib
from datetime import datetime
import json
import logging
import math
import os
import shutil

# third-party
import pytest

# local
from idaes.core.dmf import resource
from idaes.core.dmf.resource import Predicates
from idaes.core.dmf.util import mkdtemp

# for testing
from .util import init_logging

__author__ = "Dan Gunter"

init_logging()
_log = logging.getLogger(__name__)

test_version = "1.2.3"


@pytest.fixture
def example_resource():
    r = resource.Resource()
    r.v["version_info"]["version"] = test_version
    r.v["collaborators"] = [
        {"name": "Clark Kent", "email": "ckent@dailyplanet.com"},
        {"name": "Superman", "email": "sman@fortress.solitude.org"},
    ]
    r.v["datafiles"].append({"path": "/etc/passwd"})
    r.v["aliases"] = ["test_resource_full"]
    r.v["tags"] = ["test", "resource"]
    # r.relations  -- deal with this separately
    r.data = {"arbitrary": {"values": [1, 2, 3]}}
    return r


@pytest.fixture
def default_resource():
    return resource.Resource()


@pytest.mark.unit
def test_resource_roundtrip(example_resource):
    """Build up a resource with all attributes,
    then make sure it serializes and deserializes.
    """
    r = example_resource  # alias
    # make sure we can serialize it
    r_json = json.dumps(r.v, indent=2)
    # make sure we can deserialize its JSON
    r2_value = json.loads(r_json)
    # reconstruct from the deserialized JSON
    r2 = resource.Resource(value=r2_value)
    # compare them
    assert r2.v["version_info"] == r.v["version_info"]
    assert r2.data["arbitrary"]["values"][0] == 1


_propd = [
    {
        "name": "Viscosity Value",
        "units": "mPa-s",
        "values": [2.6, 6.2],
        "error_type": "absolute",
        "errors": [0.06, 0.004],
        "type": "property",
    },
    {"name": "r", "units": "", "values": [0.2, 1000], "type": "state"},
]

_propm = {
    "datatype": "MEA",
    "info": "J. Chem. Eng. Data, 2009, Vol 54, pg. 3096-30100",
    "notes": "r is MEA weight fraction in aqueous soln.",
    "authors": "Amundsen, T.G., Lars, E.O., Eimer, D.A.",
    "title": "Density and Viscosity of Monoethanolamine + .etc.",
    "date": "1970-01-01",
}


@pytest.mark.unit
def test_validate_default(default_resource):
    default_resource.validate()


@pytest.mark.unit
def test_validate_example(example_resource):
    example_resource.validate()


@pytest.mark.unit
def test_validate_preprocess(default_resource):
    r = default_resource
    r.validate()
    r.v["created"] = "2012-01-01"
    r.v["modified"] = "2012-01-01"
    r.v["version_info"]["created"] = "2011-01-01"
    r.validate()
    with pytest.raises(ValueError):
        r.v["created"] = "None of your business"
        r.validate()


@pytest.fixture
def tmpd():
    d = mkdtemp(prefix="test_resource_", suffix=".idaes")
    yield d
    shutil.rmtree(d)


@pytest.mark.unit
def test_get_datafiles_relative(default_resource, tmpd):
    r = default_resource
    r.v["datafiles_dir"] = tmpd  # paths are now relative to this
    paths = set()
    for i in range(3):
        name = "file{:d}".format(i)
        path = os.path.join(tmpd, name)
        open(path, "w").write("hello {:d}\n".format(i))
        r.v["datafiles"].append({"path": name})
        paths.add(path)
    # check that these, and only these, files are returned
    for f in r.get_datafiles():
        assert str(f) in paths
        paths.remove(str(f))
    assert len(paths) == 0  # all were returned


@pytest.mark.unit
def test_get_datafiles_absolute(default_resource, tmpd):
    r = default_resource
    paths = set()
    for i in range(3):
        name = "file{:d}".format(i)
        path = os.path.join(tmpd, name)
        open(path, "w").write("hello {:d}\n".format(i))
        # note: unlike "..._relative()", add the full path
        r.v["datafiles"].append({"path": path})
        paths.add(path)
    # check that these, and only these, files are returned
    for f in r.get_datafiles():
        assert str(f) in paths
        paths.remove(str(f))
    assert len(paths) == 0  # all were returned


@pytest.mark.unit
def test_create_relation(default_resource, example_resource):
    r1, r2 = default_resource, example_resource
    relation = resource.Triple(r1, Predicates.uses, r2)
    resource.create_relation(relation)
    with pytest.raises(ValueError):  # duplicate will raise ValueError
        resource.create_relation(relation)
    assert len(r1.v["relations"]) == 1
    assert len(r2.v["relations"]) == 1
    assert r1.v["relations"][0]["role"] == "subject"
    assert r2.v["relations"][0]["role"] == "object"
    assert r1.v["relations"][0]["identifier"] == r2.v[r2.ID_FIELD]
    assert r2.v["relations"][0]["identifier"] == r1.v[r2.ID_FIELD]
    # some errors
    with pytest.raises(ValueError):
        resource.create_relation("foo", "bad predicate", "bar")
    # delete relation from subject to test duplicate check for object
    r1.v["relations"] = []
    with pytest.raises(ValueError):  # dup raises ValueError
        resource.create_relation(relation)


@pytest.mark.unit
def test_repr(example_resource):
    txt = example_resource._repr_text_()
    assert len(txt) > 0


@pytest.mark.unit
def test_date_float():
    now = datetime.now()
    now_float = now.timestamp()
    assert resource.date_float(now) == now_float
    # tuples
    now_tuple = tuple(list(now.timetuple()[:6]) + [now.microsecond, now.tzinfo])
    assert resource.date_float(now_tuple) == now_float
    with pytest.raises(ValueError):
        resource.date_float(tuple("garbage"))
    # datetime
    assert resource.date_float(datetime(*now_tuple)) == now_float
    # strings
    with pytest.raises(ValueError):
        resource.date_float("garbage")
    assert resource.date_float(now.isoformat()) == now_float
    # int/float
    assert resource.date_float(now_float) == now_float
    assert resource.date_float(int(now_float)) == math.floor(now_float)
    #    with pytest.raises(ValueError):
    #        resource.date_float(1e12)
    #    with pytest.raises(ValueError):
    #        resource.date_float(1000000000000)
    # none
    with pytest.raises(ValueError):
        resource.date_float(None)


@pytest.mark.unit
def test_version_list():
    f = resource.version_list
    # any tuple prefix is fine
    assert f(1) == f((1,)) == f((1, 0)) == f((1, 0, 0)) == f((1, 0, 0, ""))
    # same thing with a string is fine
    assert f("1") == f("1.0") == f("1.0.0")
    # check all bad possibilities
    with pytest.raises(ValueError):
        f(())
    with pytest.raises(ValueError):
        f(1.0)
    with pytest.raises(ValueError):
        f(("a",))
    with pytest.raises(ValueError):
        f(("a", "b"))
    with pytest.raises(ValueError):
        f((1, 2, 3, None))
    with pytest.raises(ValueError):
        f((1, 2, 3, datetime))
    with pytest.raises(ValueError):
        f("1.2.3.4")  # last bit can't start with '.'
    assert f("1.2.3RC3") == [1, 2, 3, "RC3"]  # extra
    assert f("1.2.3-RC3") == [1, 2, 3, "RC3"]  # stripped "-"


@pytest.mark.unit
def test_format_version():
    assert resource.format_version([1, 2, 3, "RC3"]) == "1.2.3-RC3"


@pytest.mark.unit
def test_identifier_str():
    assert len(resource.identifier_str()) > 1
    assert resource.identifier_str("0" * 32) == "0" * 32
    with pytest.raises(ValueError):
        resource.identifier_str("foobar")
    with pytest.raises(ValueError):
        resource.identifier_str("0" * 31 + "X")


@pytest.mark.unit
def test_dirty_bit():
    dd = resource.Dict({"value": 1})
    assert dd.is_dirty()
    dd.set_clean()
    assert not dd.is_dirty()
    dd["value"] = 2
    assert dd.is_dirty()
    dd.set_clean()
    assert not dd.is_dirty()


@pytest.mark.unit
def test_validate_onlywhendirty(default_resource):
    r = default_resource
    assert r._validations == 0
    r.validate()
    assert r._validations == 1
    r.validate()
    assert r._validations == 1
    r.v["tags"] = ["one tag"]
    r.validate()
    assert r._validations == 2
    r.validate()
    assert r._validations == 2


@pytest.mark.unit
def test_triple_from_resource_relations():
    i = "cookie monster"
    j = "cookies"
    # cookie monster likes cookies
    d = {
        resource.RR_ROLE: resource.RR_SUBJ,
        resource.RR_ID: j,
        resource.RR_PRED: "likes",
    }
    t = resource.triple_from_resource_relations(i, d)
    assert t.subject == i
    assert t.predicate == "likes"
    assert t.object == j
    # cookies are liked by cookie monster
    # result is same relation triple as above
    d = {
        resource.RR_ROLE: resource.RR_OBJ,
        resource.RR_ID: i,
        resource.RR_PRED: "likes",
    }
    t = resource.triple_from_resource_relations(j, d)
    assert t.subject == i
    assert t.predicate == "likes"
    assert t.object == j
