##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes".
##############################################################################
"""
Tests for DMF Jupyter "magics"
"""
# stdlib
import os
# third-party
import pytest
# local
from idaes.dmf import magics, DMF
from idaes.dmf.resource import Resource
from .util import TempDir


class MockShell(object):
    """Mock object for IPython 'shell'.
    """
    def ev(self, name):
        raise ValueError('Unknown object: {}'.format(name))

###############################################################################
# Fixtures


@pytest.fixture
def magics_impl():
    shell = MockShell()
    return magics.DmfMagicsImpl(shell)


@pytest.fixture
def tmp_dmf():
    with TempDir() as d:
        dmf = DMF(d, name='tempdmf', create=True)
        yield dmf
        del dmf


@pytest.fixture
def tmp_magics():
    with TempDir() as d:
        magics_impl.dmf_init(d, 'create')
        yield magics_impl
        del magics_impl

###############################################################################
# Tests


def test_init_create(magics_impl):
    print('@@ create impl id={}'.format(id(magics_impl)))
    with TempDir() as d:
        magics_impl.dmf_init(d, 'create')
        assert magics_impl.initialized


def test_init_existing(tmp_dmf, magics_impl):
    magics_impl.dmf_init(tmp_dmf.root)


def test_init_extraignored(tmp_dmf, magics_impl):
    magics_impl.dmf_init(tmp_dmf.root, 'foo')  # just generate warning
    magics_impl.dmf_init(tmp_dmf.root, 'foo', 'bar')  # ditto


def test_init_required(magics_impl):
    magics_impl.dmf_info()
    assert magics_impl.last_ok is False
    magics_impl.dmf_help('anything')
    assert magics_impl.last_ok is False


def test_init_badpath(magics_impl):
    nosuchpath = os.path.join(os.path.sep, *map(str, range(10)))
    magics_impl.dmf_init(nosuchpath)
    assert not magics_impl.last_ok
    magics_impl.dmf_init(nosuchpath, 'create')
    assert not magics_impl.last_ok


def test_init_goodpath(magics_impl):
    with TempDir() as goodpath:
        magics_impl.dmf_init(goodpath, 'create')
    assert magics_impl.last_ok


def test_dmf_cmd(magics_impl):
    magics_impl.dmf('not a command')  # unrecognized command
    assert magics_impl.last_ok is False
    magics_impl.dmf('list stuff')  # init required
    assert magics_impl.last_ok is False
    magics_impl.dmf('')  # default is "help", therefore init required
    assert magics_impl.last_ok is False


def test_dmf_workspaces(magics_impl):
    with TempDir() as goodpath:
        # Try with a good path, but no workspace
        # so should return "No valid workspaces found"
        msg = magics_impl.dmf_workspaces(goodpath)
        # make sure first token of response is not a number
        token = msg.split()[0]
        pytest.raises(ValueError, int, token)
        # Create a workspace
        magics_impl.dmf_init(goodpath, 'create')
        # Try again, should return "1 workspace found" (or similar)
        msg = magics_impl.dmf_workspaces(goodpath)
        token = msg.split()[0]
        assert int(token) > 0


def test_dmf_list(magics_impl):
    # Start before setting up DMF; should not work
    magics_impl.dmf('list')
    assert not magics_impl.last_ok
    # Now set up the DMF, add 1 resource, call again. This time,
    # should be OK
    with TempDir() as wspath:
        magics_impl.dmf_init(wspath, 'create')
        rsrc = Resource()
        magics_impl._dmf.add(rsrc)  # XXX: fragile since it accesses _dmf
        magics_impl.dmf('list')
    assert magics_impl.last_ok
