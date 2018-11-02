"""
Tests for module 'codesearch'.
"""
# stdlib
import os
import shutil
import tempfile
# third-party
import pytest
# package
import idaes
from idaes.dmf import codesearch
from idaes.core import property_meta


# Helper classes and functions


class ListVisitorProperty(codesearch.PropertyMetadataVisitor):
    def __init__(self):
        self.items = []
        super(ListVisitorProperty, self).__init__()

    def visit_metadata(self, obj, meta):
        self.items.append(meta)


def _visit(walker):
    lv = ListVisitorProperty()
    walker.walk(lv)
    assert len(lv.items) == 3


# Walker and ModuleClassWalker


def test_walker():
    w = codesearch.Walker()
    w.walk(ListVisitorProperty())

# Note: exclude_dirs=contrib is necessary because the 'contrib'
# directories run things on import (!) This should become superfluous
# once the DMF and core are separated from contrib in Gthub.


def test_from_pkg():
    walker = codesearch.ModuleClassWalker(from_pkg=idaes,
                                          class_expr='_TestClass',
                                          parent_class=codesearch._TestClass,
                                          suppress_warnings=True,
                                          exclude_dirs=['contrib'])
    _visit(walker)


def test_from_path_only():
    idpath = idaes.__path__[0]
    walker = codesearch.ModuleClassWalker(from_path=idpath,
                                          class_expr='_TestClass',
                                          parent_class=codesearch._TestClass,
                                          suppress_warnings=True,
                                          exclude_dirs=['contrib'])
    _visit(walker)


def test_from_path_and_from_pkg():
    idpath = idaes.__path__[0]
    walker = codesearch.ModuleClassWalker(from_path=idpath,
                                          from_pkg=idaes,
                                          class_expr='_TestClass',
                                          parent_class=codesearch._TestClass,
                                          suppress_warnings=True,
                                          exclude_dirs=['contrib'])
    _visit(walker)


def test_from_pkg_no_path():
    import os
    try:
        codesearch.ModuleClassWalker(from_pkg=os,
                                     class_expr='_TestClass',
                                     parent_class=codesearch._TestClass,
                                     suppress_warnings=True,
                                     exclude_dirs=['contrib'])
    except IOError:
        pass
    else:
        assert False


def test_missing_args():
    pytest.raises(ValueError, codesearch.ModuleClassWalker,
                  class_expr='_TestClass')


def test_bad_path():
    try:
        codesearch.ModuleClassWalker(from_path='/a/b/c', from_pkg=idaes,
                                     class_expr='_TestClass')
    except IOError:
        pass
    else:
        assert False


def test_warnings():
    walker = codesearch.ModuleClassWalker(from_pkg=idaes,
                                          class_expr='_TestClass',
                                          parent_class=codesearch._TestClass,
                                          suppress_warnings=False,
                                          exclude_dirs=['contrib'])
    _visit(walker)


# Visitor-related tests

def test_visitor():
    v = codesearch.Visitor()
    v.visit('foo')


class MetadataHaver(property_meta.HasPropertyClassMetadata):
    @classmethod
    def define_metadata(cls, m):
        return


def test_metadata_visitor():
    mv = codesearch.PropertyMetadataVisitor()
    mv.visit(MetadataHaver())


def test_print_metadata_visitor():
    pmv = codesearch.PrintPropertyMetadataVisitor()
    pcm = property_meta.PropertyClassMetadata()
    pcm.add_properties({'foo': {'method': 'bogus'}})
    pmv.visit_metadata(MetadataHaver, pcm)


# very low-level tests


@pytest.fixture()
def tmpd():
    d = tempfile.mkdtemp(suffix='codesearch_get_modules')
    yield d
    shutil.rmtree(d)


def test__get_modules(tmpd):
    # make a temporary directory with some files and subdirs
    dirs = [('a',), ('b',), ('c',), ('a', '1'), ('a', '2'), ('a', '2', 'i')]
    for d in dirs:
        os.mkdir(os.path.join(tmpd, *d))
    files_skip = [('__init__.py',), ('foo.sh',),
                  ('a', 'bar.py'),
                  ('a', '1', '__init__.py'),
                  ('a', '2', 'i', '__init__.py'),
                  ('c', '__init__.py'), ('c', 'data.dat')]
    for f in files_skip:
        path = os.path.join(tmpd, *f)
        open(path, 'w').write('SKIP\n')
    files_find = [('foo.py',),
                  ('a', 'bar.py'),
                  ('a', '1', 'file.py'),
                  ('a', '2', 'i', 'file.py'),
                  ('b', 'foo.py')]
    wlk = codesearch.ModuleClassWalker(from_path=tmpd)
    expect_modules = set()
    for f in files_find:
        path = os.path.join(tmpd, *f)
        open(path, 'w').write('FIND\n')
        modname = '.'.join([wlk._pkg] + list(f))[:-3]  # strip '.py'
        expect_modules.add(modname)
    for mod in wlk._get_modules():
        assert mod in expect_modules
        expect_modules.remove(mod)
    assert not expect_modules
