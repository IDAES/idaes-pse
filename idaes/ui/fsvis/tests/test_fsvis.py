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
Tests for the IDAES Flowsheet Visualizer (IFV).

These are currently integration tests, because the start/stop the embedded HTTP server.
"""
import glob
import json
import logging
import os
from pathlib import Path
import pytest
import re
import requests
import time

from pyomo.environ import ConcreteModel, SolverFactory, Constraint, value
from idaes.core import FlowsheetBlock
from idaes.generic_models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.generic_models.unit_models import Flash
from idaes.ui.fsvis import fsvis, errors
from idaes.ui.flowsheet import validate_flowsheet


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


@pytest.mark.integration
def test_visualize(flash_model, tmp_path):
    from pathlib import Path

    flowsheet = flash_model.fs
    # Start the visualization server
    result = fsvis.visualize(flowsheet, "Flash", browser=False, save_dir=tmp_path)
    # Get the model
    resp = requests.get(f"http://127.0.0.1:{result.port}/fs?id=Flash")
    data = resp.json()
    # Validate the model
    ok, msg = validate_flowsheet(data)
    assert ok, f"Invalid flowsheet returned: {msg}"
    assert data["model"]["id"] == "Flash"
    assert data["model"]["unit_models"]["flash"]["type"] == "flash"
    assert len(data["cells"]) == 7
    units = [x for x in data["cells"] if x["type"] == "standard.Image"]
    assert len(units) == 4
    unit_images = [Path(x["attrs"]["image"]["xlinkHref"]).name for x in units]
    unit_images.sort()
    assert unit_images == ["feed.svg", "flash.svg", "product.svg", "product.svg"]
    # Modify the model by deleting its one and only component
    flowsheet.del_component("flash")
    # Get the model (again)
    resp = requests.get(f"http://127.0.0.1:{result.port}/fs?id=Flash")
    data = resp.json()
    # Validate the modified model
    expected = {
        "model": {
            "id": "Flash",
            "stream_table": {"columns": ["", "Variable"], "data": [], "index": []},
            "unit_models": {},
            "arcs": {},
        },
        "cells": [],
    }
    assert data == expected


@pytest.mark.integration
def test_save_visualization(flash_model, tmp_path):
    # view logs from the persistence module
    logging.getLogger("idaes.ui.fsvis").setLevel(logging.DEBUG)
    flowsheet = flash_model.fs
    # Start the visualization server, using temporary save location
    save_location = tmp_path / "flash-vis.json"
    fsvis_result = fsvis.visualize(
        flowsheet, "Flash", browser=False, save=save_location, save_dir=tmp_path
    )
    # Check the contents of the saved file are the same as what is returned by the server
    with open(fsvis_result.store.filename) as fp:
        file_data = json.load(fp)
    resp = requests.get(f"http://127.0.0.1:{fsvis_result.port}/fs?id=Flash")
    net_data = resp.json()
    assert file_data == net_data


def _canonicalize(d):
    for cell in d["cells"]:
        if "ports" in cell:
            items = cell["ports"]["items"]
            cell["ports"]["items"] = sorted(items, key=lambda x: x["group"])


@pytest.mark.unit
def test_invoke(flash_model):
    # from inspect import signature -- TODO: use for checking params
    from idaes.ui import fsvis as fsvis_pkg

    functions = {
        "method": getattr(flash_model.fs, "visualize"),
        "package": getattr(fsvis_pkg, "visualize"),
        "module": getattr(fsvis, "visualize"),
    }
    # TODO: check params


@pytest.mark.unit
def test_visualize_fn(flash_model):
    flowsheet = flash_model.fs
    result = fsvis.visualize(flowsheet, browser=False, save=False)
    assert result.store.filename == ""
    #
    for bad_save_as in (1, "/no/such/file/exists.I.hope", flowsheet):
        with pytest.raises(errors.VisualizerError):
            fsvis.visualize(flowsheet, save=bad_save_as, browser=False)


@pytest.mark.unit
def test_flowsheet_name(flash_model, tmp_path):
    raw_name = "Hello World"
    result = fsvis.visualize(
        flash_model.fs, name=raw_name, browser=False, save_dir=tmp_path
    )
    assert re.search(raw_name, result.store.filename)


@pytest.mark.unit
def test_mock_webbrowser(flash_model):
    from idaes.ui.fsvis import fsvis

    wb = fsvis.webbrowser
    for wb_mock in (MockWB(True), MockWB(False)):
        fsvis.webbrowser = wb_mock
        _ = fsvis.visualize(flash_model.fs, save=False)
    fsvis.webbrowser = wb


class MockWB:
    def __init__(self, ok):
        self.ok = ok

    def open(self, *args):
        return self.ok


## Test saving of the status file


@pytest.fixture
def save_files_prefix(tmp_path):
    value = str(tmp_path / "test_visualize")
    # clear out any cruft
    for filename in glob.glob(str(tmp_path / "test_visualize*")):
        os.unlink(filename)
    yield value
    # clear out any cruft (2)
    for filename in glob.glob(str(tmp_path / "test_visualize*")):
        os.unlink(filename)


@pytest.mark.unit
def test_visualize_save_versions(flash_model, save_files_prefix):
    # test versioned file saves
    flowsheet = flash_model.fs
    path = Path(save_files_prefix + "_save")
    work_dir = path.parent
    fs_name = path.name
    for i in range(4):
        save_arg = (True, None)[i % 2]  # try both kinds of 'use default' values
        if i < 3:
            result = fsvis.visualize(
                flowsheet, fs_name, save_dir=work_dir, browser=False, save=save_arg
            )
            if i == 0:
                assert re.search(f"{path.name}.json", result.store.filename)
            else:
                assert re.search(f"{path.name}.*{i}.*\.json", result.store.filename)
        else:
            msv, fsvis.MAX_SAVED_VERSIONS = fsvis.MAX_SAVED_VERSIONS, i - 1
            with pytest.raises(RuntimeError):
                fsvis.visualize(flowsheet, fs_name, save_dir=work_dir, browser=False)
            fsvis.MAX_SAVED_VERSIONS = msv


@pytest.mark.unit
def test_visualize_save_explicit(flash_model, save_files_prefix):
    # test explicit filename
    flowsheet = flash_model.fs
    howdy = Path(save_files_prefix + "_howdy")
    result = fsvis.visualize(flowsheet, "flowsheet", save=howdy, browser=False)
    assert re.search(howdy.name, result.store.filename)
    # overwrite but this time break explicit file into relative name and directory
    result = fsvis.visualize(
        flowsheet,
        "flowsheet",
        save=howdy.name,
        save_dir=howdy.parent,
        browser=False,
        overwrite=True,
    )
    assert re.search(howdy.name, result.store.filename)


@pytest.mark.unit
def test_visualize_save_cannot(flash_model, tmp_path):
    flowsheet = flash_model.fs
    with pytest.raises(errors.VisualizerError):
        fsvis.visualize(flowsheet, "foo", save="foo", save_dir=Path("/a/b/c/d/e/f/g"))


@pytest.mark.unit
def test_visualize_save_overwrite(flash_model, save_files_prefix):
    flowsheet = flash_model.fs
    howdy = Path(save_files_prefix + "_howdy")
    howdy.open("w").write("howdy")
    howdy_stat = os.stat(howdy)
    result = fsvis.visualize(
        flowsheet, "flowsheet", save=howdy, overwrite=True, browser=False
    )
    howdy_stat2 = os.stat(result.store.filename)
    assert (
        howdy_stat2.st_mtime > howdy_stat.st_mtime
    )  # modification time should be later


@pytest.mark.unit
def test_pick_default_save_location():
    from idaes.ui.fsvis.fsvis import _pick_default_save_location as pdsl
    p = pdsl("foo", None)
    assert str(p).endswith("foo.json")
    p = pdsl("foo", Path("/a"))
    assert p == Path("/a") / "foo.json"


@pytest.mark.unit
def test_existing_save_path(tmp_path):
    from idaes.ui.fsvis.fsvis import _handle_existing_save_path as hesp
    name = "foo"
    save_path = tmp_path / (name + ".json")
    # not there
    p = hesp(name, save_path)
    assert p == save_path
    # version 1
    save_path.open("w").write("hello")
    p1 = hesp(name, save_path)
    assert p1 != save_path
    # version 2
    p1.open("w").write("hello")
    p2 = hesp(name, save_path)
    assert str(p2) > str(p1)
    # version too far
    p2.open("w").write("hello")
    with pytest.raises(errors.TooManySavedVersions):
        p3 = hesp(name, save_path, max_versions=2)
    # infinite versions
    p4 = hesp(name, save_path, max_versions=0)
    assert str(p4) > str(p2)
    # overwrite
    p0 = hesp(name, save_path, overwrite=True)
    assert p0 == save_path


@pytest.mark.component
def test_loop_forever():
    from threading import Thread
    for quietness in (True, False):
        thr = Thread(target=fsvis._loop_forever, args=(quietness,))
        thr.setDaemon(True)
        thr.start()
        # wait a while, make sure it's still alive
        time.sleep(3)
        assert thr.is_alive()
    # threads should die when process exits
