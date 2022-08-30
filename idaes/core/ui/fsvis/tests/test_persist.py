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
Tests for the `idaes.core.ui.fsvis.persist` module
"""
# stdlib
import json

# ext
import pytest

# pkg
from idaes.core.ui.fsvis import persist, errors

# === Data ===

bad_data = {"foo": pytest}  # can't serialize as JSON
bad_data_str = "Once upon a time.."
data = {"foo": [{"bar": 123}]}
data_str = json.dumps(data)

# === Tests ===


@pytest.mark.unit
def test_file_data_store(tmp_path):
    p = tmp_path / "test.json"
    store = persist.FileDataStore(p)
    assert p == store.path
    _save_and_load_data(store)


@pytest.mark.unit
def test_memory_data_store():
    store = persist.MemoryDataStore()
    _save_and_load_data(store)


@pytest.mark.unit
def test_compare_and_str(tmp_path):
    p = tmp_path / "test.json"
    fstore = persist.FileDataStore(p)
    mstore = persist.MemoryDataStore()
    assert fstore != mstore
    assert fstore == fstore
    assert mstore == mstore
    str(fstore)
    str(mstore)


@pytest.mark.unit
def test_datastoremanager_create_add(tmp_path):
    p = tmp_path / "test.json"
    fstore = persist.FileDataStore(p)
    mstore = persist.MemoryDataStore()
    dsm = persist.DataStoreManager()
    assert dsm.add("memory", mstore)
    assert not dsm.add("memory", mstore)  # idempotence
    assert dsm.add("file", fstore)
    assert not dsm.add("file", fstore)  # idempotence
    assert dsm.add("file", mstore)  # changed type
    assert dsm.add("file", fstore)  # changed type again
    assert not dsm.add("file", fstore)  # idempotence


@pytest.mark.unit
def test_datastoremanager_save_load(tmp_path):
    dsm = persist.DataStoreManager()
    dsm.add("foo", persist.FileDataStore(tmp_path / "test.json"))
    dsm.add("bar", persist.MemoryDataStore())
    _save_and_load_data_dsm("foo", dsm)
    _save_and_load_data_dsm("bar", dsm)


# === Functions ===


def _save_and_load_data(store):
    with pytest.raises(errors.DatastoreError):
        store.save(bad_data)
    with pytest.raises(errors.DatastoreError):
        store.save(bad_data_str)
    store.save(data)
    assert data == store.load()
    store.save(json.dumps(data))
    assert data == store.load()


def _save_and_load_data_dsm(id_, dsm):
    with pytest.raises(errors.DatastoreError):
        dsm.save(id_, bad_data)
    dsm.save(id_, data)
    result = dsm.load(id_)
    assert data == result
