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
Tests for DMF Jupyter "magics"
"""
# stdlib
import os
import sys
import webbrowser

# third-party
import pytest

# local
from idaes.dmf import magics, DMF
from idaes.dmf.magics import DMFMagicError
from idaes.dmf.resource import Resource
from .util import TempDir

if sys.platform.startswith("win"):
    pytest.skip("skipping DMF tests on Windows", allow_module_level=True)

# monkey-patch webbrowser to do nothing
webbrowser.open_new = lambda url: None


class MockShell(object):
    """Mock object for IPython 'shell'.
    """

    def ev(self, name):
        """Mock evaluation.
        For most names, return the DMF class.
        For the special name "FAIL", raise an exception.
        For the special name "NODOCS", return a non-documented Python module.
        """
        if name == "FAIL":
            raise ValueError("Failure")
        elif name == "NODOCS":
            return os
        return DMF


###############################################################################
# Fixtures


@pytest.fixture
def magics_impl():
    shell = MockShell()
    return magics.DmfMagicsImpl(shell)


@pytest.fixture
def tmp_dmf():
    with TempDir() as d:
        dmf = DMF(d, name="tempdmf", create=True)
        yield dmf
        del dmf


@pytest.fixture
def tmp_magics():
    with TempDir() as d:
        magics_impl.dmf_init(d, "create")
        yield magics_impl
        del magics_impl


###############################################################################
# Tests


@pytest.mark.unit
def test_init_create(magics_impl):
    with TempDir() as d:
        magics_impl.dmf_init(d, "create")
        assert magics_impl.initialized


@pytest.mark.unit
def test_init_existing(tmp_dmf, magics_impl):
    magics_impl.dmf_init(tmp_dmf.root)


@pytest.mark.unit
def test_init_extraignored(tmp_dmf, magics_impl):
    magics_impl.dmf_init(tmp_dmf.root, "foo")  # just generate warning
    magics_impl.dmf_init(tmp_dmf.root, "foo", "bar")  # ditto


@pytest.mark.unit
def test_init_required(magics_impl):
    pytest.raises(DMFMagicError, magics_impl.dmf_info)
    pytest.raises(DMFMagicError, magics_impl.dmf_help, "anything")


@pytest.mark.unit
def test_init_badpath(magics_impl):
    nosuchpath = os.path.join(os.path.sep, *map(str, range(10)))
    pytest.raises(DMFMagicError, magics_impl.dmf_init, nosuchpath)
    pytest.raises(DMFMagicError, magics_impl.dmf_init, nosuchpath, "create")


@pytest.mark.unit
def test_init_goodpath(magics_impl):
    with TempDir() as goodpath:
        magics_impl.dmf_init(goodpath, "create")


@pytest.mark.unit
def test_dmf_cmd(magics_impl):
    pytest.raises(DMFMagicError, magics_impl.dmf, "not a command")
    pytest.raises(DMFMagicError, magics_impl.dmf, "list stuff")
    magics_impl.dmf("")  # Empty value is OK


@pytest.mark.unit
def test_dmf_workspaces(magics_impl):
    with TempDir() as goodpath:
        # Try with a good path, but no workspace
        # so should return "No valid workspaces found"
        msg = magics_impl.dmf_workspaces(goodpath)
        # make sure first token of response is not a number
        token = msg.split()[0]
        pytest.raises(ValueError, int, token)
        # Create a workspace
        magics_impl.dmf_init(goodpath, "create")
        # Try again, should return "1 workspace found" (or similar)
        msg = magics_impl.dmf_workspaces(goodpath)
        token = msg.split()[0]
        assert int(token) > 0


@pytest.mark.unit
def test_dmf_list(magics_impl):
    # Start before setting up DMF; should not work
    pytest.raises(DMFMagicError, magics_impl.dmf, "list")
    # Now set up the DMF, add 1 resource, call again. This time,
    # should be OK
    with TempDir() as wspath:
        magics_impl.dmf_init(wspath, "create")
        rsrc = Resource()
        magics_impl._dmf.add(rsrc)  # XXX: fragile since it accesses _dmf
        magics_impl.dmf("list")


@pytest.mark.unit
def test_dmf_info_initrequired(magics_impl):
    # should fail with no DMF
    pytest.raises(DMFMagicError, magics_impl.dmf_info)


@pytest.mark.unit
def test_dmf_info_topics(magics_impl):
    # should fail with DMF, and topics (not implemented)
    with TempDir() as wspath:
        magics_impl.dmf_init(wspath, "create")
        pytest.raises(DMFMagicError, magics_impl.dmf_info, "uptime")


@pytest.mark.unit
def test_dmf_info_notopics(magics_impl):
    # should succeed with DMF, and no topics
    with TempDir() as wspath:
        magics_impl.dmf_init(wspath, "create")
        magics_impl.dmf_info()


@pytest.mark.unit
def test_dmf_info_extrameta(magics_impl):
    # by filling in the metadata, this will test the dmf_info code
    # that prints out list and dict data structures
    with TempDir() as wspath:
        open(os.path.join(wspath, DMF.WORKSPACE_CONFIG), "w").write(
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
        magics_impl.dmf_init(wspath)
        magics_impl.dmf_info()


@pytest.mark.unit
def test_dmf_help(magics_impl):
    with TempDir() as wspath:
        magics_impl.dmf_init(wspath, "create")
        magics_impl.dmf_help()


@pytest.mark.unit
def test_dmf_help_badargs(magics_impl):
    with TempDir() as wspath:
        magics_impl.dmf_init(wspath, "create")
        # This will generate a warning, but is OK. Only the first
        # object is tried (which fails; another warning)
        result = magics_impl.dmf_help("this", "that")
        assert result is None  # None means "OK"


@pytest.mark.unit
def test_dmf_help_obj(magics_impl):
    with TempDir() as wspath:
        magics_impl.dmf_init(wspath, "create")
        # "special" names
        for name in "dmf", "help", "idaes":
            magics_impl.dmf_help(name)
        # object
        magics_impl.dmf_help("idaes.dmf.dmfbase.DMF")
        # failure (still returns None)
        assert magics_impl.dmf_help("FAIL") is None

