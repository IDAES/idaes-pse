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
Tests for DMF Jupyter "magics"
"""
# stdlib
from pathlib import Path
import os
from tempfile import TemporaryDirectory
from typing import Union
import webbrowser

# third-party
import pytest

# local
from idaes.core.dmf import magics, DMF
from idaes.core.dmf.magics import DMFMagicError
from idaes.core.dmf.resource import Resource

# monkey-patch webbrowser to do nothing
webbrowser.open_new = lambda url: None


class MockShell(object):
    """Mock object for IPython 'shell'."""

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


@pytest.fixture
def magics_impl():
    shell = MockShell()
    return magics.DmfMagicsImpl(shell)


scratch_dir: Union[str, None] = None
scratch_path: Union[Path, None] = None


def setup_module(module):
    global scratch_dir, scratch_path
    scratch_dir = TemporaryDirectory(prefix="idaes.core.dmf_")  # easier to remove later
    scratch_path = Path(scratch_dir.name)


def teardown_module(module):
    global scratch_dir
    del scratch_dir


@pytest.mark.unit
def test_init_create(magics_impl):
    tmp_dir = scratch_path / "init_create"
    dmf = DMF(path=tmp_dir, create=True)
    magics_impl.dmf_init(str(tmp_dir), "create")
    assert magics_impl.initialized


@pytest.mark.unit
def test_init_existing(magics_impl):
    tmp_dir = scratch_path / "init_existing"
    dmf = DMF(path=tmp_dir, create=True)
    magics_impl.dmf_init(str(tmp_dir))


@pytest.mark.unit
def test_init_extraignored(magics_impl):
    tmp_dir = scratch_path / "init_extraignored"
    dmf = DMF(path=tmp_dir, create=True)
    magics_impl.dmf_init(str(tmp_dir), "foo")  # just generate warning
    magics_impl.dmf_init(str(tmp_dir), "foo", "bar")  # ditto


@pytest.mark.unit
def test_init_required(magics_impl):
    pytest.raises(DMFMagicError, magics_impl.dmf_status)
    pytest.raises(DMFMagicError, magics_impl.dmf_help, "anything")


@pytest.mark.unit
def test_init_badpath(magics_impl):
    nosuchpath = os.path.join(os.path.sep, *map(str, range(10)))
    pytest.raises(DMFMagicError, magics_impl.dmf_init, nosuchpath)
    pytest.raises(DMFMagicError, magics_impl.dmf_init, nosuchpath, "create")


@pytest.mark.unit
def test_init_goodpath(magics_impl):
    tmp_dir = scratch_path / "init_goodpath"
    magics_impl.dmf_init(str(tmp_dir), "create")


@pytest.mark.unit
def test_dmf_cmd(magics_impl):
    pytest.raises(DMFMagicError, magics_impl.dmf, "not a command")
    pytest.raises(DMFMagicError, magics_impl.dmf, "list stuff")
    magics_impl.dmf("")  # Empty value is OK


@pytest.mark.unit
def test_dmf_workspaces(magics_impl):
    tmp_dir = scratch_path / "dmf_workspaces"
    # Try with a good path, but no workspace
    # so should return "No valid workspaces found"
    msg = magics_impl.dmf_workspaces(tmp_dir)
    # make sure first token of response is not a number
    token = msg.split()[0]
    pytest.raises(ValueError, int, token)
    # Create a workspace
    magics_impl.dmf_init(tmp_dir, "create")
    # Try again, should return "1 workspace found" (or similar)
    msg = magics_impl.dmf_workspaces(tmp_dir)
    token = msg.split()[0]
    assert int(token) > 0


@pytest.mark.unit
def test_dmf_list(magics_impl):
    # Start before setting up DMF; should not work
    pytest.raises(DMFMagicError, magics_impl.dmf, "list")
    # Now set up the DMF, add 1 resource, call again. This time,
    # should be OK
    tmp_dir = scratch_path / "dmf_list"
    magics_impl.dmf_init(tmp_dir, "create")
    rsrc = Resource()
    magics_impl._dmf.add(rsrc)  # XXX: fragile since it accesses _dmf
    magics_impl.dmf("list")


@pytest.mark.unit
def test_dmf_status_initrequired(magics_impl):
    # should fail with no DMF
    pytest.raises(DMFMagicError, magics_impl.dmf_status)


@pytest.mark.unit
def test_dmf_status_topics(magics_impl):
    tmp_dir = scratch_path / "dmf_status_topics"
    # should fail with DMF, and topics (not implemented)
    magics_impl.dmf_init(str(tmp_dir), "create")
    pytest.raises(DMFMagicError, magics_impl.dmf_status, "uptime")
    # should succeed with DMF, and no topics
    magics_impl.dmf_status()


@pytest.mark.unit
def test_dmf_status_extra_meta(magics_impl):
    tmp_dir = scratch_path / "dmf_status_extra_meta"
    tmp_dir.mkdir()
    # by filling in the metadata, this will test the dmf_status code
    # that prints out list and dict data structures
    (tmp_dir / DMF.WORKSPACE_CONFIG).open("w").write(
        f"""
_id: this-is-a-temporary-config
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
        output: {tmp_dir / 'experiment.log'}
    # user
    crazy.little.logger:
        level: error
        output: _stderr_
        """
    )
    magics_impl.dmf_init(str(tmp_dir))
    magics_impl.dmf_status()


@pytest.mark.unit
def test_dmf_help(magics_impl):
    tmp_dir = scratch_path / "dmf_help"
    magics_impl.dmf_init(str(tmp_dir), "create")
    magics_impl.dmf_help()
    # fail with bad args
    result = magics_impl.dmf_help("this", "that")
    assert result is None  # None means "OK"
    # with objects
    # a) "special" names
    for name in "dmf", "help", "idaes":
        magics_impl.dmf_help(name)
    # b) object
    magics_impl.dmf_help("idaes.core.dmf.dmfbase.DMF")
    # c) failure (still returns None)
    assert magics_impl.dmf_help("FAIL") is None
