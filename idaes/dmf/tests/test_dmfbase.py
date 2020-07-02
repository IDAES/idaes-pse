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
Tests for idaes.dmf.dmfbase module
"""
import json
import logging
import os
import sys

# third-party
import pytest

# package
from idaes.dmf import resource
from idaes.dmf import errors
from idaes.dmf.dmfbase import DMFConfig, DMF
from idaes.util.system import mkdtemp, NamedTemporaryFile
from .util import init_logging, tmp_dmf, TempDir

__author__ = "Dan Gunter <dkgunter@lbl.gov>"

if sys.platform.startswith("win"):
    pytest.skip("skipping DMF tests on Windows", allow_module_level=True)

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


@pytest.fixture(scope="function")
def tmp_propdata_file(request):
    sdir = getattr(request.module, "scratchdir", "/tmp")
    prop = tmp_propdata()
    tmpf = open(os.path.join(sdir, "resource.json"), "w")
    json.dump(prop, tmpf)
    tmpf.close()
    return tmpf


def tmp_propdata():
    return prop_json[0]


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
def test_add_property_data(tmp_dmf, tmp_propdata_file):
    tmpf, prop = tmp_propdata_file, tmp_propdata()
    # Add the resource
    r = resource.Resource(type_="property_data")
    r.set_id()
    r.v["creator"] = {"name": "Dan Gunter"}
    m = prop["meta"]
    work = '{authors}, "{title}". {info}, {date}'.format(**m)
    r.v["sources"].append({"source": work, "date": m["date"]})
    r.data = {"notes": m["notes"]}
    r.v["tags"].append("MEA")
    r.v["datafiles"].append({"path": tmpf.name})
    rid = tmp_dmf.add(r)
    assert rid is not None
    # Retrieve the resource
    r2 = tmp_dmf.fetch_one(rid)
    # Validate the resource
    assert r2.type == "property_data"
    assert "MEA" in r2.v["tags"]
    # Remove the resource
    tmp_dmf.remove(identifier=rid)


@pytest.mark.unit
def test_property_data_file(tmp_dmf, tmp_propdata_file):
    tmpf, prop = tmp_propdata_file, tmp_propdata()
    # Add the resource
    r = resource.Resource(type_="property_data")
    r.v["datafiles"].append({"path": tmpf.name})
    rid = tmp_dmf.add(r)
    assert rid is not None
    r2 = tmp_dmf.fetch_one(rid)
    path = r2.v["datafiles"][0]["path"]
    f2 = open(path, "r")
    j2 = json.load(f2)
    assert j2 == prop


@pytest.mark.unit
def test_find_propertydata(tmp_dmf):
    # populate DMF with some property data resources
    pj = prop_json[0]
    n = 10
    for i in range(n):
        pd = resource.Resource(value={"data": pj}, type_=resource.TY_PROPERTY)
        tmp_dmf.add(pd)
    # get them back again
    filter_ = {"type": resource.TY_PROPERTY}
    pdata = list(tmp_dmf.find(filter_dict=filter_))
    assert len(pdata) == n


# @pytest.mark.unit
# def test_dmf_init_minimal():
#    pytest.raises(errors.DMFError, DMF)


@pytest.mark.unit
def test_dmf_init_strfile():
    with TempDir() as tmpdir:
        open(os.path.join(tmpdir, DMF.WORKSPACE_CONFIG), "w").write("Hello")
        pytest.raises(errors.WorkspaceError, DMF, path=tmpdir)


@pytest.mark.unit
def test_dmf_init_badfile():
    with TempDir() as tmpdir:
        open(os.path.join(tmpdir, DMF.WORKSPACE_CONFIG), "w").write("Hello: There")
        pytest.raises(errors.WorkspaceError, DMF, path=tmpdir)


@pytest.mark.unit
def test_dmf_init_logconf():
    with TempDir() as tmpdir:
        open(os.path.join(tmpdir, DMF.WORKSPACE_CONFIG), "w").write(
            """
_id: this-is-a-temporary-config
logging:
    idaes.dmf.dmfbase:
        level: debug
        output: _stderr_
    root:
        output: _stdout_
    dmf:
        output: _stdout_
    .dmf.experiment:
        output: _stdout_
    # equivalent to previous
    idaes.dmf.experiment:
        output: /tmp/experiment.log
    # user
    crazy.little.logger:
        level: error
        output: _stderr_
        """
        )
        d = DMF(path=tmpdir)


@pytest.mark.unit
def test_dmf_init_logconf_badlevel():
    with TempDir() as tmpdir:
        open(os.path.join(tmpdir, DMF.WORKSPACE_CONFIG), "w").write(
            """
_id: this-is-a-temporary-config
logging:
    root:
        level: debug
    idaes.dmf.util:
        level: "This is not a valid level"
        output: _stderr_
        """
        )
        pytest.raises(errors.DMFError, DMF, path=tmpdir)


@pytest.mark.unit
def test_dmf_init_logconf_badfile():
    with TempDir() as tmpdir:
        open(os.path.join(tmpdir, DMF.WORKSPACE_CONFIG), "w").write(
            """
_id: this-is-a-temporary-config
logging:
    root:
        level: debug
    idaes.dmf.util:
        output: {}
        """.format(
                os.path.join(os.path.sep, *map(str, range(10)))
            )
        )
        pytest.raises(errors.DMFError, DMF, path=tmpdir)


@pytest.mark.unit
def test_dmf_init_workspace_name():
    with TempDir() as tmpdir:
        open(os.path.join(tmpdir, DMF.WORKSPACE_CONFIG), "w").write(
            "_id: this-is-a-temporary-config"
        )
        d = DMF(path=tmpdir, name="my workspace", desc="It is a great place to work")


@pytest.mark.unit
def test_dmf_change_traits():
    with TempDir() as tmpdir:
        open(os.path.join(tmpdir, DMF.WORKSPACE_CONFIG), "w").write(
            "_id: this-is-a-temporary-config"
        )
        d = DMF(path=tmpdir, name="my workspace", desc="It is a great place to work")
        assert d.db_file
        d.db_file = "newdb.json"
        assert d.db_file == "newdb.json"


@pytest.mark.unit
def test_dmf_add(tmp_dmf):
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

    tmp_dmf.add(r)

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
    tmp_dmf.add(r)


@pytest.mark.unit
def test_dmf_add_duplicate(tmp_dmf):
    r = resource.Resource(value={"desc": "test resource"})
    tmp_dmf.add(r)
    pytest.raises(errors.DuplicateResourceError, tmp_dmf.add, r)


@pytest.mark.unit
def test_dmf_add_filesystem_err(tmp_dmf):
    r = resource.Resource(value={"desc": "test resource"})
    # create datafile
    tmpf1 = NamedTemporaryFile(delete=False)
    tmpf1.close()
    r.v["datafiles"].append({"path": tmpf1.name})
    # now, to get an error, make the DMF datafile path unwritable
    path = os.path.join(tmp_dmf.root, tmp_dmf.datafile_dir)
    os.chmod(path, 000)
    # then try to add the resource, which includes copying the file into
    # the (now unwritable) directory
    pytest.raises(errors.DMFError, tmp_dmf.add, r)
    # make the directory writable again so we can remove it
    os.chmod(path, 0o777)


@pytest.mark.unit
def test_dmf_add_tmp_no_copy(tmp_dmf):
    r = resource.Resource(value={"desc": "test resource"})
    # create datafile, with temporary-file flag turned on
    tmpdir = mkdtemp()
    tmpfile = os.path.join(tmpdir, "foo")
    open(tmpfile, "w")
    r.v["datafiles"].append({"path": tmpfile, "is_tmp": True, "do_copy": True})
    # we want an error trying to COPY this file; to get this,
    # change the permissions of the directory
    os.chmod(tmpdir, 0o400)
    ok = False
    try:
        tmp_dmf.add(r)
    except errors.DMFError:
        ok = True
    finally:
        # change it back and clean up
        os.chmod(tmpdir, 0o700)
        os.unlink(tmpfile)
        os.rmdir(tmpdir)
    if not ok:
        assert False, "DMFError expected"


@pytest.mark.unit
def test_dmf_add_tmp_no_unlink(tmp_dmf):
    r = resource.Resource(value={"desc": "test resource"})
    # create datafile, with temporary-file flag turned on
    tmpdir = mkdtemp()
    tmpfile = os.path.join(tmpdir, "foo")
    open(tmpfile, "w")
    r.v["datafiles"].append({"path": tmpfile, "is_tmp": True, "do_copy": True})
    # we want an error trying to UNLINK this file; to get this,
    # change the permissions of the dir read-only
    os.chmod(tmpdir, 0o500)
    try:
        tmp_dmf.add(r)
    finally:
        # change it back and clean up
        os.chmod(tmpdir, 0o700)
        os.unlink(tmpfile)
        os.rmdir(tmpdir)


@pytest.mark.unit
def test_dmf_update(tmp_dmf):
    ids = add_resources(tmp_dmf, 2)
    r1 = tmp_dmf.fetch_one(ids[0])
    r1.v[r1.TYPE_FIELD] = "test"
    r1.v["desc"] = "Updated description"
    tmp_dmf.update(r1)
    r1b = tmp_dmf.fetch_one(ids[0])
    assert r1b.v["desc"] == "Updated description"
    r2 = tmp_dmf.fetch_one(ids[1])
    assert r2.v["desc"] != "Updated description"


@pytest.mark.unit
def test_dmf_update_newtype(tmp_dmf):
    ids = add_resources(tmp_dmf, 1)
    r1 = tmp_dmf.fetch_one(ids[0])
    r1.v[r1.TYPE_FIELD] = "this type is different"
    try:
        tmp_dmf.update(r1)
    except errors.DMFError:
        pass
    else:
        assert False, "DMFError expected for update() with new type"


@pytest.mark.unit
def test_dmf_remove(tmp_dmf):
    n = 10
    ids = add_resources(tmp_dmf, num=n)
    assert tmp_dmf.count() == n
    while n > 0:
        n = n - 1
        tmp_dmf.remove(ids[n])
        assert tmp_dmf.count() == n


@pytest.mark.unit
def test_dmf_remove_filter(tmp_dmf):
    n = 10
    ids = add_resources(tmp_dmf, num=n)
    assert tmp_dmf.count() == n
    # remove half of the added resources
    print("@@ remove half")
    tmp_dmf.remove(filter_dict={"data.i": {"$lt": n / 2}})
    n2 = tmp_dmf.count()
    assert n2 == n / 2
    # try to remove the same group (should do nothing
    print("@@ remove more")
    tmp_dmf.remove(filter_dict={"data.i": {"$lt": n / 2}})
    n2 = tmp_dmf.count()
    assert tmp_dmf.count() == n / 2
    # remove the rest
    print("@@ remove the rest")
    tmp_dmf.remove(filter_dict={"data.i": {"$ge": n / 2}})
    assert tmp_dmf.count() == 0


@pytest.mark.component
def test_dmf_find(tmp_dmf):
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
            tmp_dmf, num=n, tags=["all", batch], version_info={"version": version}
        )
        all_ids.extend(ids)
    if _log.isEnabledFor(logging.DEBUG):
        r = tmp_dmf.fetch_one(all_ids[0])
        _log.debug("First resource:\n{}".format(r))
    # Find all records, 2 ways
    total_num = batchsz * numbatches
    result = list(tmp_dmf.find())
    assert len(result) == total_num
    result = list(tmp_dmf.find({"tags": ["all"]}))
    assert len(result) == total_num
    # Find with 'all'
    result = list(tmp_dmf.find({"tags!": ["all", "batch1"]}))
    assert len(result) == batchsz


@pytest.mark.unit
def test_dmf_str(tmp_dmf):
    s = str(tmp_dmf)
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
    tmpfile = NamedTemporaryFile()
    DMFConfig._filename = tmpfile.name
    yield tmpfile
    tmpfile.close()
    DMFConfig._filename = default_filename


@pytest.fixture
def dmfconfig_none():
    """Default file is in user's home directory.
    Replace it with a nonexistent file.
    """
    default_filename = DMFConfig._filename
    DMFConfig._filename = os.path.join(os.path.sep, "idaes", *map(str, range(20)))
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
    dmfconfig_tmp.write(b"{[\n")
    dmfconfig_tmp.file.flush()
    pytest.raises(ValueError, DMFConfig)


@pytest.mark.unit
def test_dmfconfig_somefile(dmfconfig_tmp):
    dmfconfig_tmp.write(b"workspace: foobar\n")
    dmfconfig_tmp.file.flush()
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

