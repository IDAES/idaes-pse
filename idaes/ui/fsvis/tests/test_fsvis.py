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
import json
import logging

import pytest
import requests

from pyomo.environ import ConcreteModel, SolverFactory, Constraint, value
from idaes.core import FlowsheetBlock
from idaes.generic_models.properties.activity_coeff_models.BTX_activity_coeff_VLE \
    import BTXParameterBlock
from idaes.generic_models.unit_models import Flash
from idaes.ui.fsvis import fsvis
from idaes.ui.flowsheet import validate_flowsheet


@pytest.fixture(scope="module")
def flash_model():
    """Flash unit model. Use '.fs' attribute to get the flowsheet.
    """
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    # Flash properties
    m.fs.properties = BTXParameterBlock(default={"valid_phase": ('Liq', 'Vap'),
                                                 "activity_coeff_model": "Ideal",
                                                 "state_vars": "FTPz"})
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
def test_visualize(flash_model):
    from pathlib import Path
    flowsheet = flash_model.fs
    # Start the visualization server
    port = fsvis.visualize(flowsheet, "Flash", browser=False, save_as=None)
    # Get the model
    resp = requests.get(f"http://127.0.0.1:{port}/fs?id=Flash")
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
    resp = requests.get(f"http://127.0.0.1:{port}/fs?id=Flash")
    data = resp.json()
    # Validate the modified model
    expected = {'model': {'id': 'Flash', 'unit_models': {}, 'arcs': {}}, 'cells': []}
    assert data == expected


@pytest.mark.integration
def test_save_visualization(flash_model, tmp_path):
    # view logs from the persistence module
    logging.getLogger("idaes.ui.fsvis").setLevel(logging.DEBUG)
    flowsheet = flash_model.fs
    # Start the visualization server, using temporary save location
    save_location = tmp_path / "flash-vis.json"
    port = fsvis.visualize(flowsheet, "Flash", browser=False, save_as=save_location)
    # Check the contents of the saved file are the same as what is returned by the server
    with open(save_location) as fp:
        file_data = json.load(fp)
    resp = requests.get(f"http://127.0.0.1:{port}/fs?id=Flash")
    net_data = resp.json()
    assert file_data == net_data


def _canonicalize(d):
    for cell in d["cells"]:
        if "ports" in cell:
            items = cell["ports"]["items"]
            cell["ports"]["items"] = sorted(items, key=lambda x: x["group"])
