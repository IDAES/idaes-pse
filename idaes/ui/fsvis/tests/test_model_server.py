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
Tests for model_server module
"""
# stdlib
# ext
import pytest
from pyomo.environ import ConcreteModel

# pkg
from idaes.ui.fsvis import model_server, errors, persist
from idaes.core import FlowsheetBlock
from idaes.generic_models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.generic_models.unit_models import Flash


@pytest.mark.unit
def test_flowsheet_server_class():
    srv = model_server.FlowsheetServer()
    assert srv.port is not None


@pytest.mark.unit
def test_save_flowsheet(flash_model):
    srv = model_server.FlowsheetServer()
    with pytest.raises(errors.ProcessingError):
        srv.save_flowsheet("oscar", {})
    fs = flash_model.fs
    srv.add_flowsheet("oscar", fs, persist.MemoryDataStore())
    with pytest.raises(errors.ProcessingError):
        srv.save_flowsheet("oscar", {"invalid": pytest})


@pytest.mark.unit
def test_update_flowsheet(flash_model):
    srv = model_server.FlowsheetServer()
    with pytest.raises(errors.FlowsheetUnknown):
        srv.update_flowsheet("oscar")
    # add a flowsheet,
    fs = flash_model.fs
    srv.add_flowsheet("oscar", fs, persist.MemoryDataStore())
    # update it with no change
    srv.update_flowsheet("oscar")
    # change and update
    fs.flash.inlet.flow_mol.fix(2)  # orig value = 1
    srv.update_flowsheet("oscar")
    # Put in a bad value, DEPENDS ON PROTECTED ATTR
    srv._flowsheets["oscar"] = None
    with pytest.raises(errors.ProcessingError):
        srv.update_flowsheet("oscar")
    # Delete it, DEPENDS ON PROTECTED ATTR
    del srv._flowsheets["oscar"]
    # Update should fail not-found
    with pytest.raises(errors.FlowsheetNotFoundInMemory):
        srv.update_flowsheet("oscar")


@pytest.fixture(scope="module")
def flash_model():
    """Flash unit model. Use '.fs' attribute to get the flowsheet."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    # Flash properties
    m.fs.properties = BTXParameterBlock(
        default={
            "valid_phase": ("Liq", "Vap"),
            "activity_coeff_model": "Ideal",
            "state_vars": "FTPz",
        }
    )
    # Flash unit
    m.fs.flash = Flash(default={"property_package": m.fs.properties})
    m.fs.flash.inlet.flow_mol.fix(1)
    m.fs.flash.inlet.temperature.fix(368)
    m.fs.flash.inlet.pressure.fix(101325)
    m.fs.flash.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.flash.inlet.mole_frac_comp[0, "toluene"].fix(0.5)
    m.fs.flash.heat_duty.fix(0)
    m.fs.flash.deltaP.fix(0)
    return m


@pytest.mark.component
def test_flowsheet_server_run(flash_model):
    import requests

    srv = model_server.FlowsheetServer()
    srv.start()
    srv.path = "/app"
    resp = requests.get(f"http://localhost:{srv.port}/app")
    assert not resp.ok
    # ok to get /app with bogus id (id is just added to response page)
    resp = requests.get(f"http://localhost:{srv.port}/app?id=1234")
    assert resp.ok
    # not ok to get /fs with bogus id
    resp = requests.get(f"http://localhost:{srv.port}/fs?id=1234")
    assert not resp.ok
    # add the flowsheet
    fs = flash_model.fs
    srv.add_flowsheet("oscar", fs, persist.MemoryDataStore())
    # now /fs should work
    resp = requests.get(f"http://localhost:{srv.port}/fs?id=oscar")
    assert resp.ok
    print("Bogus PUT")
    resp = requests.put(f"http://localhost:{srv.port}/fs")
    assert not resp.ok


class MockRequest:
    """Minimal mock requests.Request object."""

    request_data = ""

    class MockFile:
        def close(self):
            return

        def readline(self, *args):
            return ""

        def read(self, len):
            return bytes(self._data.encode("utf-8"))

    def makefile(self, *args):
        mf = self.MockFile()
        # pass down request data
        mf._data = self.request_data
        return mf

    def sendall(self, b):
        return


class MockServer(model_server.FlowsheetServer):
    pass

class TestFlowsheetServerHandler(model_server.FlowsheetServerHandler):
    _code = None

    def send_error(self, code, *args, **kwargs):
        self._code = code


@pytest.fixture
def put_handler():
    h = TestFlowsheetServerHandler(MockRequest(), "127.0.0.1", MockServer())
    h.path = "/fs?id=identifier"
    h.headers = {"Content-Length": 0}
    h.requestline = ""
    h.request_version = "HTTP/1.1"
    h.command = "PUT"
    return h


@pytest.fixture
def get_handler():
    h = TestFlowsheetServerHandler(MockRequest(), "127.0.0.1", MockServer())
    h.path = "/fs?id=identifier"
    h.headers = {"Content-Length": 0}
    h.requestline = ""
    h.request_version = "HTTP/1.1"
    h.command = "GET"
    return h


@pytest.mark.unit
def test_do_put_base(put_handler):
    put_handler.do_PUT()


@pytest.mark.unit
def test_do_put_invalid_flowsheet(put_handler):
    junk = "hello"
    put_handler.connection.request_data = junk
    put_handler.headers["Content-Length"] = len(junk)
    put_handler.do_PUT()
    assert put_handler._code == 400


@pytest.mark.unit
def test_do_put_invalid_path(put_handler):
    junk = "hello"
    put_handler.path = "/fs"  # missing: "?id=identifier"
    put_handler.connection.request_data = junk
    put_handler.headers["Content-Length"] = len(junk)
    put_handler.do_PUT()
    assert put_handler._code == 400


@pytest.mark.unit
def test_do_put_unknown_error(put_handler):
    def bogus_save(*args):
        raise RuntimeError("totally bogus")
    put_handler.server.save_flowsheet = bogus_save
    put_handler.path = "/fs?id=zzz"
    put_handler.do_PUT()
    assert put_handler._code == 500


@pytest.mark.unit
def test_do_put_malformed_url(put_handler):
    put_handler.path = "/fs?zzz"
    put_handler.do_PUT()
    assert put_handler._code == 500


@pytest.mark.unit
def test_do_get_file(get_handler, tmp_path):
    saved_static = model_server._static_dir
    # modify static path that is used as root for file serving
    model_server._static_dir = tmp_path
    # existing file: OK
    with (tmp_path / "file.txt").open("w") as f:
        f.write("hello")
    get_handler.path = "/file.txt"
    get_handler._code = 1
    get_handler.do_GET()
    assert get_handler._code == 1  # not modified (no error)
    # non-existing file: not OK
    get_handler.path = "/this-file-does-not-exit"
    get_handler.do_GET()
    assert get_handler._code != 200
