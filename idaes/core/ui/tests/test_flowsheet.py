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
import copy
import json
import numpy as np
from pathlib import Path

import pytest

from idaes.core.ui.flowsheet import (
    FlowsheetSerializer,
    FlowsheetDiff,
    validate_flowsheet,
)
from idaes.models.properties.swco2 import SWCO2ParameterBlock
from idaes.models.unit_models import Heater, PressureChanger, HeatExchanger
from idaes.models.unit_models.pressure_changer import ThermodynamicAssumption
from pyomo.environ import TransformationFactory, ConcreteModel
from pyomo.network import Arc
from idaes.core import FlowsheetBlock
from idaes.models.properties.activity_coeff_models.BTX_activity_coeff_VLE import (
    BTXParameterBlock,
)
from idaes.models.unit_models import Flash, Mixer
from .shared import dict_diff

# === Sample data ===

test_dir = Path(__file__).parent

base_model = {
    "model": {"id": "Model1", "unit_models": {}, "arcs": {}},
    "cells": {},
    "routing_config": {},
}


@pytest.fixture
def models():
    # Build a series of models where each has one more component than
    # the last, and the arcs connect the components in a loop
    models = {}
    unit_types = "mixer", "heater", "stoichiometric_reactor"
    for n in range(1, len(unit_types) + 1):
        model = copy.deepcopy(base_model)
        m = model["model"]
        m["id"] = f"Model{n}"
        m["unit_models"] = {}
        for unit_num in range(n):
            m["unit_models"][f"U{unit_num}"] = {
                "type": unit_types[unit_num],
                "image": unit_types[unit_num] + ".svg",
            }
        m["arcs"] = {}
        if n > 1:
            for arc_num in range(n):
                unit_num = arc_num
                m["arcs"][f"A{arc_num}"] = {
                    "source": f"U{unit_num}",
                    "dest": f"U{(unit_num + 1) % n}",
                    "label": f"stream {arc_num}",
                }
        # add minimal cells for each unit model and arc
        c = model["cells"] = []
        for key, value in m["unit_models"].items():
            c.append(
                {
                    "id": key,
                    "attrs": {
                        "image": {"xlinkHref": "image.svg"},
                        "root": {"title": "TITLE"},
                    },
                }
            )
        for key, value in m["arcs"].items():
            c.append(
                {
                    "id": key,
                    "source": {"id": value["source"]},
                    "target": {"id": value["dest"]},
                    "labels": [{"attrs": {"text": {"text": "LABEL"}}}],
                }
            )
        # done
        models[n] = model
    return models


@pytest.fixture(scope="module")
def demo_flowsheet():
    """Semi-complicated demonstration flowsheet."""
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.BT_props = BTXParameterBlock()
    m.fs.M01 = Mixer(property_package=m.fs.BT_props)
    m.fs.H02 = Heater(property_package=m.fs.BT_props)
    m.fs.F03 = Flash(property_package=m.fs.BT_props)
    m.fs.s01 = Arc(source=m.fs.M01.outlet, destination=m.fs.H02.inlet)
    m.fs.s02 = Arc(source=m.fs.H02.outlet, destination=m.fs.F03.inlet)
    TransformationFactory("network.expand_arcs").apply_to(m.fs)

    m.fs.properties = SWCO2ParameterBlock()
    m.fs.main_compressor = PressureChanger(
        dynamic=False,
        property_package=m.fs.properties,
        compressor=True,
        thermodynamic_assumption=ThermodynamicAssumption.isentropic,
    )

    m.fs.bypass_compressor = PressureChanger(
        dynamic=False,
        property_package=m.fs.properties,
        compressor=True,
        thermodynamic_assumption=ThermodynamicAssumption.isentropic,
    )

    m.fs.turbine = PressureChanger(
        dynamic=False,
        property_package=m.fs.properties,
        compressor=False,
        thermodynamic_assumption=ThermodynamicAssumption.isentropic,
    )
    m.fs.boiler = Heater(
        dynamic=False, property_package=m.fs.properties, has_pressure_change=True
    )
    m.fs.FG_cooler = Heater(
        dynamic=False, property_package=m.fs.properties, has_pressure_change=True
    )
    m.fs.pre_boiler = Heater(
        dynamic=False, property_package=m.fs.properties, has_pressure_change=False
    )
    m.fs.HTR_pseudo_tube = Heater(
        dynamic=False, property_package=m.fs.properties, has_pressure_change=True
    )
    m.fs.LTR_pseudo_tube = Heater(
        dynamic=False, property_package=m.fs.properties, has_pressure_change=True
    )
    return m.fs


@pytest.fixture(scope="module")
def flash_flowsheet():
    # Model and flowsheet
    m = ConcreteModel()
    m.fs = FlowsheetBlock(dynamic=False)
    # Flash properties
    m.fs.properties = BTXParameterBlock(
        valid_phase=("Liq", "Vap"), activity_coeff_model="Ideal", state_vars="FTPz"
    )
    # Flash unit
    m.fs.flash = Flash(property_package=m.fs.properties)
    # TODO: move this to fix(np.NINF, skip_validation=True) once
    # Pyomo#2180 is merged
    m.fs.flash.inlet.flow_mol[:].set_value(np.NINF, True)
    m.fs.flash.inlet.flow_mol.fix()
    m.fs.flash.inlet.temperature.fix(np.inf)
    m.fs.flash.inlet.pressure[:].set_value(np.nan, True)
    m.fs.flash.inlet.pressure.fix()
    m.fs.flash.inlet.mole_frac_comp[0, "benzene"].fix(0.5)
    m.fs.flash.inlet.mole_frac_comp[0, "toluene"].fix(0.5)
    m.fs.flash.heat_duty.fix(0)
    m.fs.flash.deltaP.fix(0)
    return m.fs


@pytest.fixture(scope="module")
def demo_flowsheet_json():
    json_file = test_dir / "demo_flowsheet.json"
    s = json_file.open().read()
    return s


@pytest.fixture(scope="module")
def flash_flowsheet_json():
    json_file = test_dir / "flash_flowsheet.json"
    s = json_file.open().read()
    return s


@pytest.fixture(scope="module")
def serialized_boiler_flowsheet_json():
    json_file = test_dir / "serialized_boiler_flowsheet.json"
    s = json_file.open().read()
    return s


# === Tests ===


@pytest.mark.unit
def test_merge(models):
    """Test the FlowsheetDiff output from the .merge() function."""
    num_models = len(models)

    # With N models, in increasing complexity, test the results of merging
    # each with with the next, including the last with the first.
    for i in range(num_models):
        next_i = (i + 1) % num_models
        old, new = models[i + 1], models[next_i + 1]
        merged = FlowsheetDiff(old, new).merged(do_copy=bool(i % 2))
        assert merged["model"] == new["model"]
        sources, dests, units = [], [], []
        for item in merged["cells"]:
            id_ = item["id"]
            if "source" in item:  # arc
                sources.append(item["source"])
                dests.append(item["target"])
            else:  # unit model
                units.append(id_)
        # Each unit ID will show up exactly once in each of these sets, except
        # when we wrap around to the start where there are no arcs
        expect_unit_ids = sorted([f"U{n}" for n in range(0, next_i + 1)])
        assert expect_unit_ids == sorted(units)
        if next_i == 0:
            assert sources == []
            assert dests == []
        else:
            assert expect_unit_ids == sorted([x["id"] for x in sources])
            assert expect_unit_ids == sorted([x["id"] for x in dests])

    # Test the results of merging each with a changed version of itself
    for i in range(1, num_models + 1):
        old, new = models[i], copy.deepcopy(models[i])
        m = new["model"]
        for key in m["unit_models"]:
            m["unit_models"][key]["image"] = "changed.svg"
        for key in m["arcs"]:
            m["arcs"][key]["label"] = "changed"
        merged = FlowsheetDiff(old, new).merged()
        assert merged["model"] == new["model"]
        for cell in merged["cells"]:
            if "source" in cell:
                # see if label was copied into layout
                assert cell["labels"][0]["attrs"]["text"]["text"] == "changed"
            else:
                assert cell["attrs"]["image"]["xlinkHref"] == "changed.svg"


@pytest.mark.unit
def test_validate_flowsheet(models):
    # these have a type error since they are not iterable at all
    pytest.raises(TypeError, validate_flowsheet, None)
    pytest.raises(TypeError, validate_flowsheet, 123)
    # these are missing the top-level keys (but are sort of iterable, so no type error)
    assert validate_flowsheet("hello")[0] is False
    assert validate_flowsheet([])[0] is False
    # empty one fails
    assert validate_flowsheet({})[0] is False
    # the minimal ones we manually constructed will pass
    for model in models.values():
        assert validate_flowsheet(model)[0]
    # now try tweaks on the minimal ones
    m = models[2]["model"]
    # remove image
    image = m["unit_models"]["U1"]["image"]
    del m["unit_models"]["U1"]["image"]
    assert validate_flowsheet(m)[0] is False
    m["unit_models"]["U1"]["image"] = image  # restore it
    # mess up a unit model ID
    m["unit_models"]["U-FOO"] = m["unit_models"]["U1"]
    del m["unit_models"]["U1"]
    assert validate_flowsheet(m)[0] is False
    m["unit_models"]["U1"] = m["unit_models"]["U-FOO"]
    del m["unit_models"]["U-FOO"]
    # mess up an arc ID
    m["arcs"]["A-FOO"] = m["arcs"]["A1"]
    del m["arcs"]["A1"]
    assert validate_flowsheet(m)[0] is False
    m["arcs"]["A1"] = m["arcs"]["A-FOO"]
    del m["arcs"]["A-FOO"]


def _canonicalize(d):
    for cell in d["cells"]:
        if "ports" in cell:
            items = cell["ports"]["items"]
            cell["ports"]["items"] = sorted(items, key=lambda x: x["id"])
        if "position" in cell:
            cell.pop("position")


@pytest.mark.component
def test_flowsheet_serializer_demo(demo_flowsheet, demo_flowsheet_json):
    """Simple regression test vs. stored data."""
    test_dict = FlowsheetSerializer(demo_flowsheet, "demo").as_dict()
    stored_dict = json.loads(demo_flowsheet_json)
    _canonicalize(test_dict)
    _canonicalize(stored_dict)
    assert json.dumps(test_dict, sort_keys=True) == json.dumps(
        stored_dict, sort_keys=True
    )


@pytest.mark.component
def test_boiler_demo(serialized_boiler_flowsheet_json):
    import idaes.models_extra.power_generation.flowsheets.supercritical_power_plant.boiler_subflowsheet_build as blr

    m, solver = blr.main()
    test_dict = FlowsheetSerializer(m.fs, "boiler").as_dict()
    stored_dict = json.loads(serialized_boiler_flowsheet_json)
    _canonicalize(test_dict)
    _canonicalize(stored_dict)
    test_json = json.dumps(test_dict, sort_keys=True)
    stored_json = json.dumps(stored_dict, sort_keys=True)
    if test_json != stored_json:
        report_failure(test_dict, stored_dict)
        pytest.fail("Serialized flowsheet does not match expected")


@pytest.mark.unit
def test_flowsheet_serializer_flash(flash_flowsheet, flash_flowsheet_json):
    """Simple regression test vs. stored data."""
    test_dict = FlowsheetSerializer(flash_flowsheet, "demo").as_dict()
    stored_dict = json.loads(flash_flowsheet_json)
    _canonicalize(test_dict)
    _canonicalize(stored_dict)
    test_json = json.dumps(test_dict, sort_keys=True)
    stored_json = json.dumps(stored_dict, sort_keys=True)
    if test_json != stored_json:
        report_failure(test_dict, stored_dict)
        pytest.fail("Serialized flowsheet does not match expected")


def report_failure(test_dict, stored_dict):
    test_json, stored_json = (json.dumps(d, indent=2) for d in (test_dict, stored_dict))
    diff = dict_diff(test_dict, stored_dict)
    print("Diff between generated dict and expected dict:")
    print(diff)


# print("---")
# print(f"Generated data (JSON):\n{test_json}")
# print("---")
# print(f"Expected data (JSON):\n{stored_json}")


def _show_json(test=None, stored=None):
    import sys

    print("-" * 60)
    print("TEST VALUE")
    json.dump(test, sys.stdout)
    print()
    print("-" * 60)
    print("STORED VALUE")
    json.dump(stored, sys.stdout)


@pytest.mark.unit
def test_flowsheet_serializer_invalid():
    m = ConcreteModel()
    pytest.raises(ValueError, FlowsheetSerializer, m, "bad")


@pytest.mark.unit
def test_flowsheet_serializer_get_unit_model_type():
    from idaes.core import MaterialBalanceType
    from idaes.models.unit_models.pressure_changer import (
        ThermodynamicAssumption,
    )
    from idaes.models.unit_models.heat_exchanger import (
        delta_temperature_underwood_callback,
    )
    from idaes.models.properties import iapws95
    from pyomo.environ import Set

    # flowsheet
    m = ConcreteModel(name="My Model")
    m.fs = FlowsheetBlock(dynamic=False)
    m.fs.prop_water = iapws95.Iapws95ParameterBlock(
        phase_presentation=iapws95.PhaseType.LG
    )

    # add & test scalar unit model
    m.fs.cond_pump = PressureChanger(
        property_package=m.fs.prop_water,
        material_balance_type=MaterialBalanceType.componentTotal,
        thermodynamic_assumption=ThermodynamicAssumption.pump,
    )
    unit_type = FlowsheetSerializer.get_unit_model_type(m.fs.cond_pump)
    assert unit_type == "pressure_changer"

    # add & test indexed unit model
    m.set_fwh = Set(initialize=[1, 2, 3, 4, 6, 7, 8])
    m.fs.fwh = HeatExchanger(
        m.set_fwh,
        delta_temperature_callback=delta_temperature_underwood_callback,
        hot_side={
            "property_package": m.fs.prop_water,
            "material_balance_type": MaterialBalanceType.componentTotal,
            "has_pressure_change": True,
        },
        cold_side={
            "property_package": m.fs.prop_water,
            "material_balance_type": MaterialBalanceType.componentTotal,
            "has_pressure_change": True,
        },
    )
    unit_type = FlowsheetSerializer.get_unit_model_type(m.fs.fwh)
    assert unit_type == "heat_exchanger"
