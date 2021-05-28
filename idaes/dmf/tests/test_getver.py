"""
Tests for idaes.dmf.getver module
"""

import idaes
from idaes import dmf
from idaes.dmf import getver
import uuid

import pytest


@pytest.fixture
def garbage():
    return uuid.uuid4().hex


@pytest.fixture
def make_pip_garbage():
    garbage = uuid.uuid4().hex
    orig, getver.Versioned.PIP = getver.Versioned.PIP, garbage
    yield garbage
    getver.Versioned.PIP = orig


@pytest.mark.unit
def test_clazz():
    v1 = getver.Versioned("idaes")
    v2 = getver.Versioned(idaes)
    assert v1.get_info() == v2.get_info()
    assert v1.name == v2.name


@pytest.mark.unit
def test_func():
    import numpy

    i1 = dmf.get_version_info(idaes)
    i2 = dmf.get_version_info(numpy)
    i3 = dmf.get_version_info("pytest")
    i4 = dmf.get_version_info("traitlets")
    assert i1 != i2 and i2 != i3
    # one type of git hash error
    with pytest.raises(getver.GitHashError):
        _ = getver.Versioned("pytest").git_hash
    # another type of git hash error
    with pytest.raises(getver.GitHashError):
        _ = getver.Versioned("traitlets").git_hash


class fake_module:
    """Solely to exercise a path where the loader has no path attr
    """
    __package__ = "pytest"  # needs to match something real
    class fake_loader:
        pass
    __loader__ = fake_loader


@pytest.mark.unit
def test_more_func():
    import six

    # doesn't have __path__ but has __loader__ with path
    dmf.get_version_info(six)
    # __loader__ without path
    dmf.get_version_info(fake_module)


@pytest.mark.unit
def test_bad_import(garbage):
    with pytest.raises(getver.ModuleImportError):
        getver.Versioned(garbage)


@pytest.mark.unit
def test_bad_pip(make_pip_garbage):
    ver = getver.Versioned("traitlets")
    assert ver is not None


@pytest.mark.unit
def test_repeated_calls():
    import time

    t0 = time.time()
    vv = getver.Versioned("idaes")
    i1 = vv.get_info()
    t1 = time.time()
    for i in range(20):
        i2 = vv.get_info()
        assert i2 == i1
    t2 = time.time()
    # time to get 1st should be much, much longer than subsequent
    assert t2 - t1 < t1 - t0
