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
Test utility functions
"""
# standard library
import importlib
import logging
import os
import shutil
import time
from unittest.mock import MagicMock, patch
import warnings

# third party
import pytest

# local
from idaes.core.dmf import dmfbase
from idaes.core.dmf.util import mkdtemp

__author__ = "Dan Gunter"

scratchdir = None


def init_logging():
    """Init logging for tests."""
    log = logging.getLogger("idaes.core.dmf")
    h = logging.StreamHandler()
    f = logging.Formatter(fmt="%(asctime)s %(name)s [%(levelname)s] %(message)s")
    h.setFormatter(f)
    log.addHandler(h)
    if os.environ.get("TEST_DEBUG", ""):
        log.setLevel(logging.DEBUG)


def eq_(var, val, msg=None):
    """Shorthand for common test assertion."""
    if msg:
        assert var == val, msg
    else:
        assert var == val


def ne_(var, val, msg=None):
    """Shorthand for common test assertion."""
    if msg:
        assert var != val, msg
    else:
        assert var != val


def black_hole():
    """Output destination that swallows up everything."""
    return open("/dev/null", "w")


def patch_modules():
    """Using mock.patch, modify sys.modules to mock up any
    modules not actually importable for tests.

    Returns:
        (mock.mock._patch_dict) Object returned by `mock.patch.dict()`.
    """
    sm = mock_import(
        "core.process_base",
        "idaes_models",
        assign={
            "core.process_base:ProcessBase": "idaes.core.dmf.tests.process_base:ProcessBase"
        },
    )
    # .. can add more to `sm` dict here..
    if sm is not None:
        result = patch.dict("sys.modules", sm)
    else:
        result = None
    return result


def mock_import(name, package, assign=None):
    """Set up mocking for imports of "<package>.<name>".
    Optionally, assign real (importable) modules and/or module members
    to the mocked ones (see `assign`).

    Args:
        name (str): Module path to mock up
        package (str): Package name to mock up
        assign (dict): Mapping where key is mock module and
           value is real module path. Optionally, a ':<member>' can
           be appended to both key and value to assign members of the
           real module to the mock. If this is None, nothing is done.

    Returns:
        (dict) New values for sys.modules dict to be patched.
          If the given `package` imports without error,
          nothing is done and return value is None.
    """
    name, package = name.strip(), package.strip()
    assert name, "module must be non-empty"
    assert package, "package must be non-empty"
    mlist = name.split(".")
    try:
        importlib.import_module(package, mlist[0])
        return None
    except ImportError:
        pass
    pkg = MagicMock()
    # add module path to mock
    modpath, curmod = package, pkg
    sys_modules = {modpath: curmod}
    for m in mlist:
        modpath = modpath + "." + m
        curmod = getattr(curmod, m)
        sys_modules[modpath] = curmod
    # assign real modules/classes to mock ones
    if assign:
        for mock_imp, real_imp in assign.items():
            # extract class, if present
            if ":" in mock_imp:
                mock_imp, mock_class = mock_imp.split(":")
                real_imp, real_class = real_imp.split(":")
            else:
                mock_class, real_class = None, None
            # import real module
            real_pkg = real_imp.split(".")[0]
            real_modules = "." + real_imp[real_imp.index(".") + 1 :]
            real_mod = importlib.import_module(real_modules, real_pkg)
            # assign real module, or class, to mock one
            mock_mod, mock_imp_list = pkg, mock_imp.split(".")
            if mock_class:
                # set mock.module.path.Class = real.module.path.Class
                for m in mock_imp_list:
                    mock_mod = getattr(mock_mod, m)
                real_classobj = getattr(real_mod, real_class)
                setattr(mock_mod, mock_class, real_classobj)
            else:
                # set mock.module.path = real.module.path
                for m in mock_imp_list[:-1]:
                    mock_mod = getattr(mock_mod, m)
                setattr(mock_mod, mock_imp_list[-1], real_mod)
    return sys_modules


class TempDir(object):
    def __init__(self, chdir=False):
        self._d = None
        self._chdir = chdir
        self._origdir = None

    def __enter__(self):
        self._d = mkdtemp(suffix="-dmf")
        if self._chdir:
            self._origdir = os.getcwd()
            os.chdir(self._d)
        return self._d

    def __exit__(self, *args):
        if self._d is not None:
            # remove files in dir
            rmdirs = []
            for dirpath, subdirs, files in os.walk(self._d):
                for f in files:
                    os.unlink(os.path.join(dirpath, f))
                rmdirs.append(dirpath)
            # remove dirs
            while rmdirs:
                path = rmdirs.pop()
                os.rmdir(path)
        self._d = None
        if self._chdir:
            os.chdir(self._origdir)


@pytest.fixture(scope="function")
def tmp_dmf():
    """Test fixture to create a DMF in a temporary
    directory.
    """
    tmpdir = mkdtemp()
    dmf = dmfbase.DMF(path=tmpdir, create=True)
    yield dmf
    removed = False
    for i in range(3):
        try:
            shutil.rmtree(tmpdir)
            removed = True
        except Exception as err:
            time.sleep(1)
    if not removed:
        warnings.warn(f"failed to remove temporary directory: {tmpdir}")
