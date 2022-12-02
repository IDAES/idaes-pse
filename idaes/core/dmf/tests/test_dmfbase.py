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
Tests for idaes.core.dmf.dmfbase module

Skip tests that do chmod() except on Linux, as Windows at least leaves
the resulting directories in an un-removable state.
"""
import json
import logging
import os
from pathlib import Path
from tempfile import TemporaryDirectory
from typing import Union

# third-party
import pytest

# package
from idaes.core.dmf import resource
from idaes.core.dmf import errors
from idaes.core.dmf.dmfbase import DMFConfig, DMF
from .util import init_logging
from idaes.core.dmf.util import NamedTemporaryFile

__author__ = "Dan Gunter <dkgunter@lbl.gov>"

init_logging()
_log = logging.getLogger(__name__)

prop_json = [
    {
        "meta": {
            "datatype": "MEA",
            "info": "J. Chem. Eng. Data, 2009, Vol 54, pg. 3096-30100",
            "notes": "r is MEA weight fraction in aqueous soln.",
            "authors": "Amundsen, T.G., Lars, E.O., Eimer, D.A.",
            "title": "Density and Viscosity of Monoethanolamine + .etc.",
            "date": "2009",
        },
        "data": [
            {
                "name": "Viscosity Value",
                "units": "mPa-s",
                "values": [2.6, 6.2],
                "error_type": "absolute",
                "errors": [0.06, 0.004],
                "type": "property",
            },
            {"name": "r", "units": "", "values": [0.2, 1000], "type": "state"},
        ],
    }
]

scratch_dir: Union[str, None] = None
scratch_dmf: Union[DMF, None] = None


def setup_module(module):
    global scratch_dir, scratch_dmf
    scratch_dir = TemporaryDirectory(prefix="idaes.core.dmf_").name
    scratch_dmf = DMF(path=scratch_dir, create=True)


def teardown_module(module):
    global scratch_dmf, scratch_dir
    del scratch_dmf
    del scratch_dir


def add_resources(dmf_obj, num=3, **attrs):
    ids = []
    for i in range(num):
        r = resource.Resource(value=attrs, type_="test")
        r.data = {"i": i}
        dmf_obj.add(r)
        ids.append(r.id)
    return ids


# Tests


@pytest.mark.unit
def test_add_property_data():
    prop = prop_json[0]
    tmp_propdata_file = (Path(scratch_dir) / "add_property_data.json").open("w")
    # Add the resource
    r = resource.Resource(type_="property_data")
    r.set_id()
    r.v["creator"] = {"name": "Dan Gunter"}
    m = prop["meta"]
    work = '{authors}, "{title}". {info}, {date}'.format(**m)
    # Only for type=publication, now
    # r.v["sources"].append({"source": work, "date": m["date"]})
    r.data = {"notes": m["notes"]}
    r.v["tags"].append("MEA")
    r.v["datafiles"].append({"path": tmp_propdata_file.name})
    rid = scratch_dmf.add(r)
    assert rid is not None
    # Retrieve the resource
    r2 = scratch_dmf.fetch_one(rid)
    # Validate the resource
    assert r2.type == "property_data"
    assert "MEA" in r2.v["tags"]
    # Remove the resource
    scratch_dmf.remove(identifier=rid)


@pytest.mark.unit
def test_property_data_file():
    tmp_propdata_file = (Path(scratch_dir) / "property_data_file.json").open("w")
    json.dump(prop_json[0], tmp_propdata_file)
    tmp_propdata_file.close()
    # Add the resource
    r = resource.Resource(type_="property_data")
    r.v["datafiles"].append({"path": tmp_propdata_file.name})
    rid = scratch_dmf.add(r)
    assert rid is not None
    r2 = scratch_dmf.fetch_one(rid)
    path = r2.v["datafiles"][0]["path"]
    f2 = open(path, "r")
    j2 = json.load(f2)
    assert j2 == prop_json[0]


@pytest.mark.unit
def test_find_propertydata():
    # populate DMF with some property data resources
    pj = prop_json[0]
    n, resource_ids = 10, []
    for i in range(n):
        pd = resource.Resource(
            value={"data": pj}, type_=resource.ResourceTypes.property
        )
        resource_ids.append(scratch_dmf.add(pd))
    # get them back again
    filter_ = {"type": resource.ResourceTypes.property}
    pdata = list(scratch_dmf.find(filter_dict=filter_))
    assert len(pdata) == n
    for rid in resource_ids:
        scratch_dmf.remove(identifier=rid)


@pytest.mark.unit
def test_dmf_init_bad_config1():
    tmp_dir = Path(scratch_dir) / "dmf_init_bad_config1"
    tmp_dir.mkdir()
    (tmp_dir / DMF.WORKSPACE_CONFIG).open("w").write("Hello")
    # try to open - should fail
    pytest.raises(errors.WorkspaceError, DMF, path=tmp_dir.name)


@pytest.mark.unit
def test_dmf_init_bad_config2():
    tmp_dir = Path(scratch_dir) / "dmf_init_bad_config2"
    tmp_dir.mkdir()
    (tmp_dir / DMF.WORKSPACE_CONFIG).open("w").write("Hello")
    # try to open - should fail
    pytest.raises(errors.WorkspaceError, DMF, path=tmp_dir.name)


@pytest.mark.unit
def test_dmf_init_logconf():
    tmp_dir = Path(scratch_dir) / "dmf_init_logconf"
    tmp_dir.mkdir()
    (tmp_dir / DMF.WORKSPACE_CONFIG).open("w").write(
        f"""_id: this-is-a-temporary-config
logging:
    idaes.core.dmf.dmfbase:
        level: debug
        output: _stderr_
    root:
        output: _stdout_
    dmf:
        output: _stdout_
    .dmf.experiment:
        output: _stdout_
    # equivalent to previous
    idaes.core.dmf.experiment:
        output: {tmp_dir / "experiment.log"}
    # user
    crazy.little.logger:
        level: error
        output: _stderr_
        """
    )
    DMF(path=str(tmp_dir))


@pytest.mark.unit
def test_dmf_init_logconf_bad():
    tmp_dir = Path(scratch_dir) / "dmf_init_logconf_badlevel"
    tmp_dir.mkdir()
    (tmp_dir / DMF.WORKSPACE_CONFIG).open("w").write(
        """
_id: this-is-a-temporary-config
logging:
    root:
        level: debug
    idaes.core.dmf.util:
        level: "This is not a valid level"
        output: _stderr_
        """
    )
    pytest.raises(errors.DMFError, DMF, path=tmp_dir)
    (tmp_dir / DMF.WORKSPACE_CONFIG).open("w").write(
        """
_id: this-is-a-temporary-config
logging:
    root:
        level: debug
    idaes.core.dmf.util:
        output: {}
        """.format(
            os.path.join(os.path.sep, *map(str, range(10)))
        )
    )
    pytest.raises(errors.DMFError, DMF, path=tmp_dir)


@pytest.mark.unit
def test_dmf_init_workspace_name():
    tmp_dir = Path(scratch_dir) / "dmf_init_workspace_name"
    tmp_dir.mkdir()
    (tmp_dir / DMF.WORKSPACE_CONFIG).open("w").write("_id: this-is-a-temporary-config")
    d = DMF(path=tmp_dir, name="my workspace", desc="It is a space to work")


@pytest.mark.unit
def test_dmf_change_traits():
    tmp_dir = Path(scratch_dir) / "dmf_change_traits"
    tmp_dir.mkdir()
    (tmp_dir / DMF.WORKSPACE_CONFIG).open("w").write("_id: this-is-a-temporary-config")
    d = DMF(path=tmp_dir, name="my workspace", desc="It is a great place to work")
    assert d.db_file
    d.db_file = "newdb.json"
    assert d.db_file == "newdb.json"


@pytest.mark.unit
def test_dmf_add():
    tmp_dir = Path(scratch_dir) / "dmf_add"
    dmf = DMF(path=tmp_dir, create=True)
    r = resource.Resource(value={"desc": "test resource"})
    r.do_copy = True  # copy by default
    # (1) Copy, and don't remove {default behavior}
    tmpf1 = NamedTemporaryFile(delete=False)
    tmpf1.close()
    r.v["datafiles"].append({"path": tmpf1.name})
    # (2) Copy, and remove original
    tmpf2 = NamedTemporaryFile(delete=False)
    tmpf2.close()
    r.v["datafiles"].append({"path": tmpf2.name, "is_tmp": True})
    # (3) Do not copy (or remove)
    tmpf3 = NamedTemporaryFile()
    r.v["datafiles"].append({"path": tmpf3.name, "do_copy": False})

    dmf.add(r)

    os.unlink(tmpf1.name)
    try:
        os.unlink(tmpf2.name)
        assert False, "Expected error"
    except Exception as err:
        pass

    os.unlink(tmpf3.name)
    # This is ignored. It makes no sense to ask the file
    # to be removed, but not copied (just a file delete?!)
    r = resource.Resource(value={"desc": "test resource"})
    r.v["datafiles"].append({"path": "foo", "do_copy": False, "is_tmp": True})
    dmf.add(r)


@pytest.mark.unit
def test_dmf_add_duplicate():
    tmp_dir = Path(scratch_dir) / "dmf_add_duplicate"
    dmf = DMF(path=tmp_dir, create=True)
    r = resource.Resource(value={"desc": "test resource"})
    dmf.add(r)
    pytest.raises(errors.DuplicateResourceError, dmf.add, r)


# @pytest.mark.linux
@pytest.mark.unit
def test_dmf_add_filesystem_err():
    tmp_dir = Path(scratch_dir) / "dmf_add_filesystem_err"
    dmf = DMF(path=tmp_dir, create=True)
    r = resource.Resource(value={"desc": "test resource"})
    # create datafile
    tmpf1 = NamedTemporaryFile(delete=False)
    tmpf1.close()
    r.v["datafiles"].append({"path": tmpf1.name})
    # now, to get an error, move the destination directory
    dest_dir = Path(dmf.root) / dmf.datafile_dir
    moved_dest_dir = Path(str(dest_dir) + "-moved")
    dest_dir.rename(moved_dest_dir)
    # then try to add the resource, which includes copying the file into
    # the (now unwritable) directory
    r.v["datafiles"][0]["do_copy"] = True  # make sure do_copy flag is on
    pytest.raises(errors.DMFError, dmf.add, r)
    # move directory back so we can remove it
    moved_dest_dir.rename(dest_dir)


@pytest.mark.linux
@pytest.mark.unit
def test_dmf_add_tmp_no_copy():
    tmp_dir = Path(scratch_dir) / "dmf_add_tmp_no_copy"
    dmf = DMF(path=tmp_dir, create=True)
    r = resource.Resource(value={"desc": "test resource"})
    # create datafile, with temporary-file flag turned on
    tmp_file = (tmp_dir / "foo").open("w")
    r.v["datafiles"].append({"path": str(tmp_file), "is_tmp": True, "do_copy": True})
    # we want an error trying to COPY this file; to get this,
    # change the permissions of the directory
    os.chmod(tmp_dir, 0o400)
    ok = False
    try:
        dmf.add(r)
    except errors.DMFError:
        ok = True
    if not ok:
        assert False, "DMFError expected"


@pytest.mark.linux
@pytest.mark.unit
def test_dmf_add_tmp_no_unlink():
    tmp_dir = Path(scratch_dir) / "dmf_add_tmp_no_unlink"
    dmf = DMF(path=tmp_dir, create=True)
    r = resource.Resource(value={"desc": "test resource"})
    # create datafile, with temporary-file flag turned on
    tmp_file = (tmp_dir / "foo").open("w")
    r.v["datafiles"].append({"path": str(tmp_file), "is_tmp": True, "do_copy": True})
    # we want an error trying to COPY this file; to get this,
    # change the permissions of the directory
    os.chmod(tmp_dir, 0o500)
    pytest.raises(Exception, dmf.add, r)


@pytest.mark.unit
def test_dmf_update():
    tmp_dir = Path(scratch_dir) / "dmf_update"
    dmf = DMF(path=tmp_dir, create=True)
    ids = add_resources(dmf, 2)
    r1 = dmf.fetch_one(ids[0])
    r1.v[r1.TYPE_FIELD] = "test"
    r1.v["desc"] = "Updated description"
    dmf.update(r1)
    r1b = dmf.fetch_one(ids[0])
    assert r1b.v["desc"] == "Updated description"
    r2 = dmf.fetch_one(ids[1])
    assert r2.v["desc"] != "Updated description"


@pytest.mark.unit
def test_dmf_update_newtype():
    tmp_dir = Path(scratch_dir) / "dmf_update_newtype"
    dmf = DMF(path=tmp_dir, create=True)
    ids = add_resources(dmf, 1)
    r1 = dmf.fetch_one(ids[0])
    r1.v[r1.TYPE_FIELD] = "this type is different"
    try:
        dmf.update(r1)
    except errors.DMFError:
        pass
    else:
        assert False, "DMFError expected for update() with new type"


@pytest.mark.unit
def test_dmf_remove():
    tmp_dir = Path(scratch_dir) / "dmf_remove"
    dmf = DMF(path=tmp_dir, create=True)
    n = 10
    ids = add_resources(dmf, num=n)
    assert dmf.count() == n
    while n > 0:
        n = n - 1
        dmf.remove(ids[n])
        assert dmf.count() == n


@pytest.mark.unit
def test_dmf_remove_filter():
    tmp_dir = Path(scratch_dir) / "dmf_remove_filter"
    dmf = DMF(path=tmp_dir, create=True)
    n = 10
    ids = add_resources(dmf, num=n)
    assert dmf.count() == n
    # remove half of the added resources
    # print("@@ remove half")
    dmf.remove(filter_dict={"data.i": {"$lt": n / 2}})
    n2 = dmf.count()
    assert n2 == n / 2
    # try to remove the same group (should do nothing
    # print("@@ remove more")
    dmf.remove(filter_dict={"data.i": {"$lt": n / 2}})
    n2 = dmf.count()
    assert dmf.count() == n / 2
    # remove the rest
    # print("@@ remove the rest")
    dmf.remove(filter_dict={"data.i": {"$ge": n / 2}})
    assert dmf.count() == 0


@pytest.mark.component
def test_dmf_find():
    tmp_dir = Path(scratch_dir) / "dmf_find"
    dmf = DMF(path=tmp_dir, create=True)
    # populate with batches of records
    # they all have the tag 'all', each batch has 'batch<N>' as well
    # All resources in a batch are given version 1.0.<N>
    # Individual resources will have data of {i: 0..<batchsz-1>}
    batchsz, numbatches = 10, 9
    all_ids = []
    for i in range(numbatches):
        n = batchsz
        batch = "batch{:d}".format(i + 1)
        version = resource.version_list([1, 0, i + 1])
        ids = add_resources(
            dmf, num=n, tags=["all", batch], version_info={"version": version}
        )
        all_ids.extend(ids)
    if _log.isEnabledFor(logging.DEBUG):
        r = dmf.fetch_one(all_ids[0])
        _log.debug("First resource:\n{}".format(r))
    # Find all records, 2 ways
    total_num = batchsz * numbatches
    result = list(dmf.find())
    assert len(result) == total_num
    result = list(dmf.find({"tags": ["all"]}))
    assert len(result) == total_num
    # Find with 'all'
    result = list(dmf.find({"tags!": ["all", "batch1"]}))
    assert len(result) == batchsz


@pytest.mark.unit
def test_dmf_str():
    tmp_dir = Path(scratch_dir) / "dmf_str"
    dmf = DMF(path=tmp_dir, create=True)
    s = str(dmf)
    assert len(s) > 0


#########################
# DMFConfig             #
#########################


@pytest.fixture
def dmfconfig_tmp():
    """Default file is in user's home directory.
    We don't want to actually modify this with a test.
    So switch it out and switch it back when the fixture
    is done.
    """
    default_filename = DMFConfig._filename
    path = Path(scratch_dir) / "config.yaml"
    DMFConfig._filename = str(path)
    yield path.open("w")
    DMFConfig._filename = default_filename


@pytest.fixture
def dmfconfig_none():
    """Default file is in user's home directory.
    Replace it with a nonexistent file.
    """
    default_filename = DMFConfig._filename
    DMFConfig._filename = os.path.join(os.path.sep, "idaes", *map(str, range(10)))
    yield True
    DMFConfig._filename = default_filename


@pytest.mark.unit
def test_dmfconfig_init_defaults_nofile(dmfconfig_none):
    config = DMFConfig()
    assert config.c == DMFConfig.DEFAULTS


@pytest.mark.unit
def test_dmfconfig_init_defaults_emptyfile(dmfconfig_tmp):
    config = DMFConfig()
    assert config.c == DMFConfig.DEFAULTS


@pytest.mark.unit
def test_dmfconfig_init_defaults2(dmfconfig_tmp):
    config = DMFConfig(defaults={"look": "here"})
    assert config.c["look"] == "here"


@pytest.mark.unit
def test_dmfconfig_bad_file(dmfconfig_tmp):
    dmfconfig_tmp.write("{[\n")
    dmfconfig_tmp.flush()
    pytest.raises(ValueError, DMFConfig)


@pytest.mark.unit
def test_dmfconfig_somefile(dmfconfig_tmp):
    dmfconfig_tmp.write("workspace: foobar\n")
    dmfconfig_tmp.flush()
    config = DMFConfig()


@pytest.mark.unit
def test_dmfconfig_save(dmfconfig_tmp):
    config = DMFConfig()
    config.save()


@pytest.mark.unit
def test_dmfconfig_save_nofile(dmfconfig_none):
    config = DMFConfig()
    pytest.raises(IOError, config.save)


@pytest.mark.unit
def test_dmfconfig_attrs(dmfconfig_tmp):
    config = DMFConfig()
    assert config.workspace is not None
