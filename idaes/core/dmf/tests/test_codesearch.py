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
Tests for module 'codesearch'.
"""
# stdlib
import logging
import os
import random
import shutil
import sys

# third-party
import pytest

# package
import idaes
from idaes.core.dmf import codesearch
from idaes.core.base import property_meta
from idaes.core.dmf.util import mkdtemp
from .util import init_logging

__author__ = "Dan Gunter"

# if sys.platform.startswith("win"):
#    pytest.skip("skipping DMF tests on Windows", allow_module_level=True)

init_logging()
_log = logging.getLogger(__name__)


# Helper classes and functions


class ListVisitorProperty(codesearch.PropertyMetadataVisitor):
    def __init__(self):
        self.items = []
        super(ListVisitorProperty, self).__init__()

    def visit_metadata(self, obj, meta):
        self.items.append(meta)


class DummyVisitor(codesearch.Visitor):
    def visit(self, obj):
        return True


def _visit(walker):
    lv = ListVisitorProperty()
    walker.walk(lv)
    assert len(lv.items) == 3


# Walker and ModuleClassWalker


@pytest.mark.unit
def test_walker():
    w = codesearch.Walker()
    w.walk(ListVisitorProperty())


# Note: exclude_dirs=contrib is necessary because the 'contrib'
# directories run things on import (!) This should become superfluous
# once the DMF and core are separated from contrib in Gthub.

# Common setup

walker_args = dict(
    class_expr="_TestClass\d?",
    parent_class=codesearch._TestClass,
    suppress_warnings=False,
)

walker_args_noex = dict(parent_class=codesearch._TestClass, suppress_warnings=False)

walker_args_noexclude = walker_args.copy()
walker_args_noexclude.update(
    dict(
        exclude_testdirs=False,
        exclude_tests=False,
        exclude_init=False,
        exclude_setup=False,
    )
)


@pytest.fixture
def package_walker():
    """This walker is also a test. But useful as a base for other things."""
    walker = codesearch.ModuleClassWalker(from_pkg=idaes, **walker_args)
    yield walker


# Test different kinds of walkers


@pytest.mark.component
def test_from_pkg_withexpr():
    walker = codesearch.ModuleClassWalker(from_pkg=idaes, **walker_args)
    _visit(walker)


@pytest.mark.component
def test_from_pkg_noexpr():
    walker = codesearch.ModuleClassWalker(from_pkg=idaes, **walker_args_noex)
    _visit(walker)


@pytest.mark.unit
def test_from_path_only():
    idpath = idaes.__path__[0]
    walker = codesearch.ModuleClassWalker(from_path=idpath, **walker_args)
    _visit(walker)


@pytest.mark.unit
def test_from_path_and_from_pkg():
    idpath = idaes.__path__[0]
    walker = codesearch.ModuleClassWalker(
        from_path=idpath, from_pkg=idaes, **walker_args
    )
    _visit(walker)


@pytest.mark.unit
def test_from_pkg_no_path():
    import os

    try:
        codesearch.ModuleClassWalker(from_pkg=os, **walker_args)
    except IOError:
        pass
    else:
        assert False


@pytest.mark.unit
def test_missing_args():
    pytest.raises(ValueError, codesearch.ModuleClassWalker, class_expr="_TestClass")


@pytest.mark.unit
def test_bad_path():
    try:
        codesearch.ModuleClassWalker(
            from_path="/a/b/c", from_pkg=idaes, class_expr="_TestClass"
        )
    except IOError:
        pass
    else:
        assert False


@pytest.mark.unit
def test_warnings():
    walker = codesearch.ModuleClassWalker(
        from_pkg=idaes,
        class_expr="_TestClass",
        parent_class=codesearch._TestClass,
        suppress_warnings=False,
        exclude_dirs=["contrib"],
    )
    _visit(walker)


@pytest.fixture
def dummy_package():
    random.seed()
    saved_sys_path = sys.path[:]
    # build a bad module in a temporary package
    d = mkdtemp()
    sys.path.append(d)
    foo = os.path.join(d, "foo{}".format(random.randint(1e6, 1e7 - 1)))
    os.mkdir(foo)
    open(os.path.join(foo, "__init__.py"), "w")
    # create some modules
    m = open(os.path.join(foo, "goodmodule.py"), "w")
    m.write("class IndexMe(object):\n    pass\n")
    m.close()
    m = open(os.path.join(foo, "badmodule.py"), "w")
    m.write("This is a bad module.\n")
    m.close()
    yield foo
    # clean up
    shutil.rmtree(d)
    sys.path = saved_sys_path


@pytest.mark.unit
def test_cannot_import(dummy_package):
    w = codesearch.ModuleClassWalker(from_path=dummy_package, suppress_warnings=False)
    w.walk(DummyVisitor())
    assert len(w.get_indexed_classes()) == 1


@pytest.mark.unit
def test_cannot_import_nowarn(dummy_package):
    w = codesearch.ModuleClassWalker(from_path=dummy_package, suppress_warnings=True)
    w.walk(DummyVisitor())
    assert len(w.get_indexed_classes()) == 1


@pytest.mark.unit
def test_get_metadata_unimplemented(dummy_package):
    print("@@ GMU sys.path={}".format(sys.path))
    f = open(os.path.join(dummy_package, "unimpl.py"), "w")
    f.write(
        "class U(object):\n"
        "    @classmethod\n"
        "    def get_metadata(cls):\n"
        "        raise NotImplementedError\n"
    )
    f.close()
    w = codesearch.ModuleClassWalker(from_path=dummy_package, suppress_warnings=True)
    visitor = codesearch.PropertyMetadataVisitor()
    w.walk(visitor)


@pytest.mark.unit
def test_noexclude(dummy_package):
    walker_args_noexclude.update(dict(from_path=dummy_package))
    w = codesearch.ModuleClassWalker(**walker_args_noexclude)
    w.walk(DummyVisitor())


@pytest.mark.unit
def test_get_indexed_classes(package_walker):
    visitor = ListVisitorProperty()
    package_walker.walk(visitor)
    class_list = package_walker.get_indexed_classes()
    assert len(class_list) == len(visitor.items)


# Visitor-related tests


@pytest.mark.unit
def test_visitor():
    v = codesearch.Visitor()
    v.visit("foo")


class MetadataHaver(property_meta.HasPropertyClassMetadata):
    @classmethod
    def define_metadata(cls, m):
        return


@pytest.mark.unit
def test_metadata_visitor():
    mv = codesearch.PropertyMetadataVisitor()
    mv.visit(MetadataHaver())


@pytest.mark.unit
def test_print_metadata_visitor():
    pmv = codesearch.PrintPropertyMetadataVisitor()
    pcm = property_meta.PropertyClassMetadata()
    pcm.add_properties({"foo": {"method": "bogus"}})
    pmv.visit_metadata(MetadataHaver, pcm)


# very low-level tests


@pytest.fixture()
def tmpd():
    d = mkdtemp(suffix="codesearch_get_modules")
    yield d
    shutil.rmtree(d)


@pytest.mark.unit
def test__get_modules(tmpd):
    # make a temporary directory with some files and subdirs
    dirs = [("a",), ("b",), ("c",), ("a", "1"), ("a", "2"), ("a", "2", "i")]
    for d in dirs:
        os.mkdir(os.path.join(tmpd, *d))
    files_skip = [
        ("__init__.py",),
        ("foo.sh",),
        ("a", "bar.py"),
        ("a", "1", "__init__.py"),
        ("a", "2", "i", "__init__.py"),
        ("c", "__init__.py"),
        ("c", "data.dat"),
    ]
    for f in files_skip:
        path = os.path.join(tmpd, *f)
        open(path, "w").write("SKIP\n")
    files_find = [
        ("foo.py",),
        ("a", "bar.py"),
        ("a", "1", "file.py"),
        ("a", "2", "i", "file.py"),
        ("b", "foo.py"),
    ]
    wlk = codesearch.ModuleClassWalker(from_path=tmpd)
    expect_modules = set()
    for f in files_find:
        path = os.path.join(tmpd, *f)
        open(path, "w").write("FIND\n")
        modname = ".".join([wlk._pkg] + list(f))[:-3]  # strip '.py'
        expect_modules.add(modname)
    for mod in wlk._get_modules():
        assert mod in expect_modules
        expect_modules.remove(mod)
    assert not expect_modules
