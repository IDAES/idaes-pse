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
Tests for the `idaes.ui.fsvis.persist` module
"""
import json

import pytest

from idaes.ui.fsvis import persist

# === Data ===

bad_data = {"foo": pytest}  # can't serialize as JSON
data = {"foo": [{"bar": 123}]}
data_str = json.dumps(data)

# === Tests ===

@pytest.mark.unit
def test_file_data_store(tmp_path):
    p = tmp_path / "test.json"
    store = persist.FileDataStore(p)
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
    pytest.raises(ValueError, store.save, bad_data)
    store.save(data)
    result = store.load()
    assert data == result


def _save_and_load_data_dsm(id_, dsm):
    pytest.raises(ValueError, dsm.save, id_, bad_data)
    dsm.save(id_, data)
    result = dsm.load(id_)
    assert data == result
