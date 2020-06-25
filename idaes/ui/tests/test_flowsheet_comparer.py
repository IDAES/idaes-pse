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
import pytest
from jsonschema.exceptions import ValidationError

from idaes.ui import flowsheet_comparer as fc

from idaes.ui.flowsheet_comparer import Action


@pytest.mark.unit
def test_compare_models():
    model1 = {'model': {
                  'id': 0, 
                  'unit_models': {
                      'M101': {'type': 'mixer', 'image': 'mixer.svg'}, 
                      'H101': {'type': 'heater', 'image': 'heater_2.svg'}, 
                      'R101': {'type': 'stoichiometric_reactor', 'image': 'reactor_s.svg'}, 
                      'F101': {'type': 'flash', 'image': 'flash.svg'}, 
                      'S101': {'type': 'separator', 'image': 'splitter.svg'}, 
                      'C101': {'type': 'pressure_changer', 'image': 'compressor.svg'}, 
                      'F102': {'type': 'flash', 'image': 'flash.svg'}}, 
                  'arcs': {
                      's03': {'source': 'M101', 'dest': 'H101', 'label': "molar flow ('Liq', 'benzene') 0.5\nmolar flow ('Liq', 'toluene') 0.5\nmolar flow ('Liq', 'hydrogen') 0.5\nmolar flow ('Liq', 'methane') 0.5\nmolar flow ('Vap', 'benzene') 0.5\nmolar flow ('Vap', 'toluene') 0.5\nmolar flow ('Vap', 'hydrogen') 0.5\nmolar flow ('Vap', 'methane') 0.5\ntemperature 298.15\npressure 10132"}, 
                      's04': {'source': 'H101', 'dest': 'R101', 'label': "molar flow ('Liq', 'benzene') 0.5\nmolar flow ('Liq', 'toluene') 0.5\nmolar flow ('Liq', 'hydrogen') 0.5\nmolar flow ('Liq', 'methane') 0.5\nmolar flow ('Vap', 'benzene') 0.5\nmolar flow ('Vap', 'toluene') 0.5\nmolar flow ('Vap', 'hydrogen') 0.5\nmolar flow ('Vap', 'methane') 0.5\ntemperature 298.15\npressure 10132"}, 
                      's05': {'source': 'R101', 'dest': 'F101', 'label': "molar flow ('Liq', 'benzene') 0.5\nmolar flow ('Liq', 'toluene') 0.5\nmolar flow ('Liq', 'hydrogen') 0.5\nmolar flow ('Liq', 'methane') 0.5\nmolar flow ('Vap', 'benzene') 0.5\nmolar flow ('Vap', 'toluene') 0.5\nmolar flow ('Vap', 'hydrogen') 0.5\nmolar flow ('Vap', 'methane') 0.5\ntemperature 298.15\npressure 10132"}, 
                      's06': {'source': 'F101', 'dest': 'S101', 'label': "molar flow ('Liq', 'benzene') 0.5\nmolar flow ('Liq', 'toluene') 0.5\nmolar flow ('Liq', 'hydrogen') 0.5\nmolar flow ('Liq', 'methane') 0.5\nmolar flow ('Vap', 'benzene') 0.5\nmolar flow ('Vap', 'toluene') 0.5\nmolar flow ('Vap', 'hydrogen') 0.5\nmolar flow ('Vap', 'methane') 0.5\ntemperature 298.15\npressure 10132"}, 
                      's08': {'source': 'S101', 'dest': 'C101', 'label': "molar flow ('Liq', 'benzene') 0.5\nmolar flow ('Liq', 'toluene') 0.5\nmolar flow ('Liq', 'hydrogen') 0.5\nmolar flow ('Liq', 'methane') 0.5\nmolar flow ('Vap', 'benzene') 0.5\nmolar flow ('Vap', 'toluene') 0.5\nmolar flow ('Vap', 'hydrogen') 0.5\nmolar flow ('Vap', 'methane') 0.5\ntemperature 298.15\npressure 10132"}, 
                      's09': {'source': 'C101', 'dest': 'M101', 'label': "molar flow ('Liq', 'benzene') 0.5\nmolar flow ('Liq', 'toluene') 0.5\nmolar flow ('Liq', 'hydrogen') 0.5\nmolar flow ('Liq', 'methane') 0.5\nmolar flow ('Vap', 'benzene') 0.5\nmolar flow ('Vap', 'toluene') 0.5\nmolar flow ('Vap', 'hydrogen') 0.5\nmolar flow ('Vap', 'methane') 0.5\ntemperature 298.15\npressure 10132"}, 
                      's10': {'source': 'F101', 'dest': 'F102', 'label': "molar flow ('Liq', 'benzene') 0.5\nmolar flow ('Liq', 'toluene') 0.5\nmolar flow ('Liq', 'hydrogen') 0.5\nmolar flow ('Liq', 'methane') 0.5\nmolar flow ('Vap', 'benzene') 0.5\nmolar flow ('Vap', 'toluene') 0.5\nmolar flow ('Vap', 'hydrogen') 0.5\nmolar flow ('Vap', 'methane') 0.5\ntemperature 298.15\npressure 10132"}}}}

    model2 = {'model': {
                'id': 0, 
                'unit_models': {
                    'M101': {'type': 'mixer', 'image': 'mixer.svg'}, 
                    'H101': {'type': 'heater', 'image': 'heater_2.svg'}, 
                    'R101': {'type': 'stoichiometric_reactor', 'image': 'reactor_s.svg'}, 
                    'F101': {'type': 'flash', 'image': 'flash.svg'}, 
                    'S101': {'type': 'separator', 'image': 'splitter.svg'}, 
                    'C101': {'type': 'pressure_changer', 'image': 'compressor.svg'}, 
                    'F102': {'type': 'flash', 'image': 'flash.svg'},
                    'F222': {'type': 'flash', 'image': 'flash.svg'}},
                'arcs': {
                    's03': {'source': 'M101', 'dest': 'H101', 'label': "Hello World!"}, 
                    's04': {'source': 'H101', 'dest': 'R101', 'label': "molar flow ('Liq', 'benzene') 0.5\nmolar flow ('Liq', 'toluene') 0.5\nmolar flow ('Liq', 'hydrogen') 0.5\nmolar flow ('Liq', 'methane') 0.5\nmolar flow ('Vap', 'benzene') 0.5\nmolar flow ('Vap', 'toluene') 0.5\nmolar flow ('Vap', 'hydrogen') 0.5\nmolar flow ('Vap', 'methane') 0.5\ntemperature 298.15\npressure 10132"}, 
                    's05': {'source': 'R101', 'dest': 'F101', 'label': "molar flow ('Liq', 'benzene') 0.5\nmolar flow ('Liq', 'toluene') 0.5\nmolar flow ('Liq', 'hydrogen') 0.5\nmolar flow ('Liq', 'methane') 0.5\nmolar flow ('Vap', 'benzene') 0.5\nmolar flow ('Vap', 'toluene') 0.5\nmolar flow ('Vap', 'hydrogen') 0.5\nmolar flow ('Vap', 'methane') 0.5\ntemperature 298.15\npressure 10132"}, 
                    's06': {'source': 'F101', 'dest': 'S101', 'label': "molar flow ('Liq', 'benzene') 0.5\nmolar flow ('Liq', 'toluene') 0.5\nmolar flow ('Liq', 'hydrogen') 0.5\nmolar flow ('Liq', 'methane') 0.5\nmolar flow ('Vap', 'benzene') 0.5\nmolar flow ('Vap', 'toluene') 0.5\nmolar flow ('Vap', 'hydrogen') 0.5\nmolar flow ('Vap', 'methane') 0.5\ntemperature 298.15\npressure 10132"}, 
                    's08': {'source': 'S101', 'dest': 'C101', 'label': "molar flow ('Liq', 'benzene') 0.5\nmolar flow ('Liq', 'toluene') 0.5\nmolar flow ('Liq', 'hydrogen') 0.5\nmolar flow ('Liq', 'methane') 0.5\nmolar flow ('Vap', 'benzene') 0.5\nmolar flow ('Vap', 'toluene') 0.5\nmolar flow ('Vap', 'hydrogen') 0.5\nmolar flow ('Vap', 'methane') 0.5\ntemperature 298.15\npressure 10132"}, 
                    's09': {'source': 'C101', 'dest': 'M101', 'label': "molar flow ('Liq', 'benzene') 0.5\nmolar flow ('Liq', 'toluene') 0.5\nmolar flow ('Liq', 'hydrogen') 0.5\nmolar flow ('Liq', 'methane') 0.5\nmolar flow ('Vap', 'benzene') 0.5\nmolar flow ('Vap', 'toluene') 0.5\nmolar flow ('Vap', 'hydrogen') 0.5\nmolar flow ('Vap', 'methane') 0.5\ntemperature 298.15\npressure 10132"}, 
                    's10': {'source': 'F101', 'dest': 'F102', 'label': "molar flow ('Liq', 'benzene') 0.5\nmolar flow ('Liq', 'toluene') 0.5\nmolar flow ('Liq', 'hydrogen') 0.5\nmolar flow ('Liq', 'methane') 0.5\nmolar flow ('Vap', 'benzene') 0.5\nmolar flow ('Vap', 'toluene') 0.5\nmolar flow ('Vap', 'hydrogen') 0.5\nmolar flow ('Vap', 'methane') 0.5\ntemperature 298.15\npressure 10132"},
                    's12': {'source': 'M101', 'dest': 'F111', 'label': "molar flow ('Liq', 'benzene') 0.5\nmolar flow ('Liq', 'toluene') 0.5\nmolar flow ('Liq', 'hydrogen') 0.5\nmolar flow ('Liq', 'methane') 0.5\nmolar flow ('Vap', 'benzene') 0.5\nmolar flow ('Vap', 'toluene') 0.5\nmolar flow ('Vap', 'hydrogen') 0.5\nmolar flow ('Vap', 'methane') 0.5\ntemperature 298.15\npressure 10132"}}}}

    diff_model, out_json = fc.compare_models(model1, model2)

    diff_model_truth = {'F222': {'type': 'flash', 'image': 'flash.svg', 'action': Action.ADD.value, "class": "unit model"}, 
                        's03': {'source': 'M101', 'dest': 'H101', 'label': 'Hello World!', 'action': Action.CHANGE.value, "class": "arc"}, 
                        's12': {'source': 'M101', 'dest': 'F111', 'label': "molar flow ('Liq', 'benzene') 0.5\nmolar flow ('Liq', 'toluene') 0.5\nmolar flow ('Liq', 'hydrogen') 0.5\nmolar flow ('Liq', 'methane') 0.5\nmolar flow ('Vap', 'benzene') 0.5\nmolar flow ('Vap', 'toluene') 0.5\nmolar flow ('Vap', 'hydrogen') 0.5\nmolar flow ('Vap', 'methane') 0.5\ntemperature 298.15\npressure 10132", 'action': Action.ADD.value, "class": "arc"}
                        }

    assert diff_model == diff_model_truth

    diff_model, out_json = fc.compare_models(model2, model1)

    diff_model_truth = {'F222': {'type': 'flash', 'image': 'flash.svg', 'action': Action.REMOVE.value, 'class': 'unit model'},
                        's03': {'source': 'M101', 'dest': 'H101', 'label': "molar flow ('Liq', 'benzene') 0.5\nmolar flow ('Liq', 'toluene') 0.5\nmolar flow ('Liq', 'hydrogen') 0.5\nmolar flow ('Liq', 'methane') 0.5\nmolar flow ('Vap', 'benzene') 0.5\nmolar flow ('Vap', 'toluene') 0.5\nmolar flow ('Vap', 'hydrogen') 0.5\nmolar flow ('Vap', 'methane') 0.5\ntemperature 298.15\npressure 10132", 'action': Action.CHANGE.value, "class": "arc"},
                        's12': {'source': 'M101', 'dest': 'F111', 'label': "molar flow ('Liq', 'benzene') 0.5\nmolar flow ('Liq', 'toluene') 0.5\nmolar flow ('Liq', 'hydrogen') 0.5\nmolar flow ('Liq', 'methane') 0.5\nmolar flow ('Vap', 'benzene') 0.5\nmolar flow ('Vap', 'toluene') 0.5\nmolar flow ('Vap', 'hydrogen') 0.5\nmolar flow ('Vap', 'methane') 0.5\ntemperature 298.15\npressure 10132", 'action': Action.REMOVE.value, "class": "arc"}}

    assert diff_model == diff_model_truth


@pytest.mark.unit
def test_compare_models_edge_cases():
    # Both empty models
    existing_model = {"model": {
                        "id": 1,
                        "unit_models": {},
                        "arcs": {}}}

    new_model = {"model": {
                   "id": 2,
                   "unit_models": {},
                   "arcs": {}}}

    diff_model, out_json = fc.compare_models(existing_model, new_model)
    assert diff_model == {}

    # Existing model is empty
    existing_model = {"model": {
                        "id": 1,
                        "unit_models": {},
                        "arcs": {}}}

    new_model = {"model": {
                   'id': 0, 
                   'unit_models': {
                       'M101': {'type': 'mixer', 'image': 'mixer.svg'}, 
                       'H101': {'type': 'heater', 'image': 'heater_2.svg'}}, 
                   'arcs': {
                       's03': {'source': 'M101', 'dest': 'H101', 'label': "hello"}}}}

    diff_model, out_json = fc.compare_models(existing_model, new_model)

    # Since the existing model is empty we return an empty diff_model and  
    # return new_model as the output json
    diff_model_truth = {}
    assert diff_model == diff_model_truth
    assert out_json == new_model

    # New model is empty
    existing_model = {"model": {
                        'id': 0, 
                        'unit_models': {
                            'M101': {'type': 'mixer', 'image': 'mixer.svg'}, 
                            'H101': {'type': 'heater', 'image': 'heater_2.svg'}}, 
                        'arcs': {
                            's03': {'source': 'M101', 'dest': 'H101', 'label': "hello"}}}}

    new_model = {"model": {
                   "id": 1,
                   "unit_models": {},
                   "arcs": {}}}

    diff_model, out_json = fc.compare_models(existing_model, new_model)

    assert diff_model == {'M101': {'type': 'mixer', 'image': 'mixer.svg', 'action': Action.REMOVE.value, 'class': 'unit model'}, 
                          'H101': {'type': 'heater', 'image': 'heater_2.svg', 'action': Action.REMOVE.value, 'class': 'unit model'},
                          's03': {'source': 'M101', 'dest': 'H101', 'label': "hello", 'action': Action.REMOVE.value, 'class': 'arc'}}

    # The models are the same
    model = {"model": {
               'id': 0, 
               'unit_models': {
                   'M101': {'type': 'mixer', 'image': 'mixer.svg'}, 
                   'H101': {'type': 'heater', 'image': 'heater_2.svg'}}, 
               'arcs': {
                   's03': {'source': 'M101', 'dest': 'H101', 'label': "hello"}}}}

    diff_model, out_json = fc.compare_models(model, model)
    assert diff_model == {}


@pytest.mark.unit
def test_compare_models_errors():
    # Both empty
    existing_model = {}
    new_model = {}

    with pytest.raises(KeyError):
        fc.compare_models(existing_model, new_model)

    # New model is missing the id
    existing_model = {"model": {
                        'id': 0, 
                        'unit_models': {
                            'M101': {'type': 'mixer', 'image': 'mixer.svg'}, 
                            'H101': {'type': 'heater', 'image': 'heater_2.svg'}}, 
                        'arcs': {
                            's03': {'source': 'M101', 'dest': 'H101', 'label': "hello"}}}}
    new_model = {"model": {
                   'unit_models': {
                       'M101': {'type': 'mixer', 'image': 'mixer.svg'}, 
                       'H101': {'type': 'heater', 'image': 'heater_2.svg'}}, 
                   'arcs': {
                       's03': {'source': 'M101', 'dest': 'H101', 'label': "hello"}}}}

    with pytest.raises(ValidationError):
        fc.compare_models(existing_model, new_model)

    # Existing model is the wrong format
    existing_model = {"model": {
                        'id': 0, 
                        'arcs': {
                            's03': {'source': 'M101', 'dest': 'H101', 'label': "hello"}}}}
    new_model = {"model": {
                   'id': 0, 
                   'unit_models': {
                       'M101': {'type': 'mixer', 'image': 'mixer.svg'}, 
                       'H101': {'type': 'heater', 'image': 'heater_2.svg'}}, 
                   'arcs': {
                       's03': {'source': 'M101', 'dest': 'H101', 'label': "hello"}}}}

    with pytest.raises(ValidationError):
        fc.compare_models(existing_model, new_model)


@pytest.mark.unit
def test_model_jointjs_conversion():
  # Test unit model addition
  original_jointjs = {"cells": [
                        {"type": "standard.Image", "position": {"x": 100, "y": 100}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "M101", "z": 1, "attrs": {"image": {"xlinkHref": "mixer.svg"}, "label": {"text": "M101"}, "root": {"title": "mixer"}}}, 
                        {"type": "standard.Image", "position": {"x": 200, "y": 200}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "H101", "z": 1, "attrs": {"image": {"xlinkHref": "heater_2.svg"}, "label": {"text": "H101"}, "root": {"title": "heater"}}}, 
                        {"type": "standard.Link", "source": {"anchor": {"name": "right", "args": {"rotate": "false", "padding": 0}}, "id": "H101"}, "target": {"anchor": {"name": "topLeft", "args": {"rotate": "false", "padding": 0}}, "id": "R101"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s04", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "hello", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2}, 
                        {"type": "standard.Link", "source": {"anchor": {"name": "bottomRight", "args": {"rotate": "false", "padding": 0}}, "id": "R101"}, "target": {"anchor": {"name": "left", "args": {"rotate": "false", "padding": 0}}, "id": "F101"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s05", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "world", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2}
                     ]}

  diff_model = {'F222': {'type': 'flash', 'image': 'flash.svg', 'action': Action.ADD.value, "class": "unit model"}}

  new_jointjs = fc.model_jointjs_conversion(diff_model, original_jointjs)
  jointjs_truth = {"cells": [
                        {"type": "standard.Image", "position": {"x": 100, "y": 100}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "M101", "z": 1, "attrs": {"image": {"xlinkHref": "mixer.svg"}, "label": {"text": "M101"}, "root": {"title": "mixer"}}}, 
                        {"type": "standard.Image", "position": {"x": 200, "y": 200}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "H101", "z": 1, "attrs": {"image": {"xlinkHref": "heater_2.svg"}, "label": {"text": "H101"}, "root": {"title": "heater"}}}, 
                        {"type": "standard.Link", "source": {"anchor": {"name": "right", "args": {"rotate": "false", "padding": 0}}, "id": "H101"}, "target": {"anchor": {"name": "topLeft", "args": {"rotate": "false", "padding": 0}}, "id": "R101"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s04", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "hello", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2}, 
                        {"type": "standard.Link", "source": {"anchor": {"name": "bottomRight", "args": {"rotate": "false", "padding": 0}}, "id": "R101"}, "target": {"anchor": {"name": "left", "args": {"rotate": "false", "padding": 0}}, "id": "F101"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s05", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "world", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2},
                        {"type": "standard.Image", "position": {"x": 150, "y": 150}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "F222", "z": 1, "attrs": {"image": {"xlinkHref": "flash.svg"}, "label": {"text": "F222"}, "root": {"title": "flash"}}},
                     ]}

  assert new_jointjs == jointjs_truth

  # Test unit model removal
  original_jointjs = {"cells": [
                        {"type": "standard.Image", "position": {"x": 100, "y": 100}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "M101", "z": 1, "attrs": {"image": {"xlinkHref": "mixer.svg"}, "label": {"text": "M101"}, "root": {"title": "mixer"}}}, 
                        {"type": "standard.Image", "position": {"x": 200, "y": 200}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "H101", "z": 1, "attrs": {"image": {"xlinkHref": "heater_2.svg"}, "label": {"text": "H101"}, "root": {"title": "heater"}}}, 
                        {"type": "standard.Link", "source": {"anchor": {"name": "right", "args": {"rotate": "false", "padding": 0}}, "id": "H101"}, "target": {"anchor": {"name": "topLeft", "args": {"rotate": "false", "padding": 0}}, "id": "R101"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s04", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "hello", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2}, 
                        {"type": "standard.Link", "source": {"anchor": {"name": "bottomRight", "args": {"rotate": "false", "padding": 0}}, "id": "R101"}, "target": {"anchor": {"name": "left", "args": {"rotate": "false", "padding": 0}}, "id": "F101"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s05", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "world", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2}
                     ]}
 
  diff_model = {'M101': {'type': 'mixer', 'image': 'mixer.svg', 'action': Action.REMOVE.value, "class": "unit model"}}

  new_jointjs = fc.model_jointjs_conversion(diff_model, original_jointjs)
  jointjs_truth = {"cells": [
                        {"type": "standard.Image", "position": {"x": 200, "y": 200}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "H101", "z": 1, "attrs": {"image": {"xlinkHref": "heater_2.svg"}, "label": {"text": "H101"}, "root": {"title": "heater"}}}, 
                        {"type": "standard.Link", "source": {"anchor": {"name": "right", "args": {"rotate": "false", "padding": 0}}, "id": "H101"}, "target": {"anchor": {"name": "topLeft", "args": {"rotate": "false", "padding": 0}}, "id": "R101"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s04", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "hello", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2}, 
                        {"type": "standard.Link", "source": {"anchor": {"name": "bottomRight", "args": {"rotate": "false", "padding": 0}}, "id": "R101"}, "target": {"anchor": {"name": "left", "args": {"rotate": "false", "padding": 0}}, "id": "F101"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s05", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "world", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2},
                     ]}

  assert new_jointjs == jointjs_truth

  # Test arc addition
  original_jointjs = {"cells": [
                        {"type": "standard.Image", "position": {"x": 100, "y": 100}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "M101", "z": 1, "attrs": {"image": {"xlinkHref": "mixer.svg"}, "label": {"text": "M101"}, "root": {"title": "mixer"}}}, 
                        {"type": "standard.Image", "position": {"x": 200, "y": 200}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "H101", "z": 1, "attrs": {"image": {"xlinkHref": "heater_2.svg"}, "label": {"text": "H101"}, "root": {"title": "heater"}}}, 
                        {"type": "standard.Link", "source": {"anchor": {"name": "right", "args": {"rotate": "false", "padding": 0}}, "id": "H101"}, "target": {"anchor": {"name": "topLeft", "args": {"rotate": "false", "padding": 0}}, "id": "R101"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s04", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "hello", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2}, 
                        {"type": "standard.Link", "source": {"anchor": {"name": "bottomRight", "args": {"rotate": "false", "padding": 0}}, "id": "R101"}, "target": {"anchor": {"name": "left", "args": {"rotate": "false", "padding": 0}}, "id": "F101"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s05", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "world", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2}
                     ]}

  diff_model = {'s01': {'source': 'M101', 'dest': 'F111', 'label': "foo", 'action': Action.ADD.value, "class": "arc"}}

  new_jointjs = fc.model_jointjs_conversion(diff_model, original_jointjs)
  jointjs_truth = {"cells": [
                        {"type": "standard.Image", "position": {"x": 100, "y": 100}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "M101", "z": 1, "attrs": {"image": {"xlinkHref": "mixer.svg"}, "label": {"text": "M101"}, "root": {"title": "mixer"}}}, 
                        {"type": "standard.Image", "position": {"x": 200, "y": 200}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "H101", "z": 1, "attrs": {"image": {"xlinkHref": "heater_2.svg"}, "label": {"text": "H101"}, "root": {"title": "heater"}}}, 
                        {"type": "standard.Link", "source": {"anchor": {"name": "right", "args": {"rotate": "false", "padding": 0}}, "id": "H101"}, "target": {"anchor": {"name": "topLeft", "args": {"rotate": "false", "padding": 0}}, "id": "R101"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s04", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "hello", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2}, 
                        {"type": "standard.Link", "source": {"anchor": {"name": "bottomRight", "args": {"rotate": "false", "padding": 0}}, "id": "R101"}, "target": {"anchor": {"name": "left", "args": {"rotate": "false", "padding": 0}}, "id": "F101"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s05", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "world", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2},
                        {"type": "standard.Link", "source": {"id": "M101"}, "target": {"id": "F111"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s01", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "foo", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2}
                     ]}

  assert new_jointjs == jointjs_truth

  # Test arc change
  original_jointjs = {"cells": [
                        {"type": "standard.Image", "position": {"x": 100, "y": 100}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "M101", "z": 1, "attrs": {"image": {"xlinkHref": "mixer.svg"}, "label": {"text": "M101"}, "root": {"title": "mixer"}}}, 
                        {"type": "standard.Image", "position": {"x": 200, "y": 200}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "H101", "z": 1, "attrs": {"image": {"xlinkHref": "heater_2.svg"}, "label": {"text": "H101"}, "root": {"title": "heater"}}}, 
                        {"type": "standard.Link", "source": {"anchor": {"name": "right", "args": {"rotate": "false", "padding": 0}}, "id": "H101"}, "target": {"anchor": {"name": "topLeft", "args": {"rotate": "false", "padding": 0}}, "id": "R101"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s04", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "hello", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2}, 
                        {"type": "standard.Link", "source": {"anchor": {"name": "bottomRight", "args": {"rotate": "false", "padding": 0}}, "id": "R101"}, "target": {"anchor": {"name": "left", "args": {"rotate": "false", "padding": 0}}, "id": "F101"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s05", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "world", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2}
                     ]}

  diff_model = {'s04': {'source': 'M101', 'dest': 'F111', 'label': "asdf", 'action': Action.CHANGE.value, "class": "arc"}}

  new_jointjs = fc.model_jointjs_conversion(diff_model, original_jointjs)
  jointjs_truth = {"cells": [
                        {"type": "standard.Image", "position": {"x": 100, "y": 100}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "M101", "z": 1, "attrs": {"image": {"xlinkHref": "mixer.svg"}, "label": {"text": "M101"}, "root": {"title": "mixer"}}}, 
                        {"type": "standard.Image", "position": {"x": 200, "y": 200}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "H101", "z": 1, "attrs": {"image": {"xlinkHref": "heater_2.svg"}, "label": {"text": "H101"}, "root": {"title": "heater"}}},  
                        {"type": "standard.Link", "source": {"anchor": {"name": "bottomRight", "args": {"rotate": "false", "padding": 0}}, "id": "R101"}, "target": {"anchor": {"name": "left", "args": {"rotate": "false", "padding": 0}}, "id": "F101"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s05", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "world", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2},
                        {"type": "standard.Link", "source": {"anchor": {"name": "right", "args": {"rotate": "false", "padding": 0}}, "id": "M101"}, "target": {"anchor": {"name": "topLeft", "args": {"rotate": "false", "padding": 0}}, "id": "F111"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s04", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "asdf", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2}
                     ]}

  assert new_jointjs == jointjs_truth

  # Test arc removal
  original_jointjs = {"cells": [
                        {"type": "standard.Image", "position": {"x": 100, "y": 100}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "M101", "z": 1, "attrs": {"image": {"xlinkHref": "mixer.svg"}, "label": {"text": "M101"}, "root": {"title": "mixer"}}}, 
                        {"type": "standard.Image", "position": {"x": 200, "y": 200}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "H101", "z": 1, "attrs": {"image": {"xlinkHref": "heater_2.svg"}, "label": {"text": "H101"}, "root": {"title": "heater"}}}, 
                        {"type": "standard.Link", "source": {"anchor": {"name": "right", "args": {"rotate": "false", "padding": 0}}, "id": "H101"}, "target": {"anchor": {"name": "topLeft", "args": {"rotate": "false", "padding": 0}}, "id": "R101"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s04", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "hello", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2}, 
                        {"type": "standard.Link", "source": {"anchor": {"name": "bottomRight", "args": {"rotate": "false", "padding": 0}}, "id": "R101"}, "target": {"anchor": {"name": "left", "args": {"rotate": "false", "padding": 0}}, "id": "F101"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s05", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "world", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2}
                     ]}

  diff_model = {'s04': {'source': 'M101', 'dest': 'F111', 'label': "asdf", 'action': Action.REMOVE.value, "class": "arc"}}

  new_jointjs = fc.model_jointjs_conversion(diff_model, original_jointjs)
  jointjs_truth = {"cells": [
                        {"type": "standard.Image", "position": {"x": 100, "y": 100}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "M101", "z": 1, "attrs": {"image": {"xlinkHref": "mixer.svg"}, "label": {"text": "M101"}, "root": {"title": "mixer"}}}, 
                        {"type": "standard.Image", "position": {"x": 200, "y": 200}, "size": {"width": 50, "height": 50}, "angle": 0, "id": "H101", "z": 1, "attrs": {"image": {"xlinkHref": "heater_2.svg"}, "label": {"text": "H101"}, "root": {"title": "heater"}}}, 
                        {"type": "standard.Link", "source": {"anchor": {"name": "bottomRight", "args": {"rotate": "false", "padding": 0}}, "id": "R101"}, "target": {"anchor": {"name": "left", "args": {"rotate": "false", "padding": 0}}, "id": "F101"}, "router": {"name": "orthogonal", "padding": 10}, "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}}, "id": "s05", "labels": [{"attrs": {"rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1}, "text": {"text": "world", "fill": "black", "text-anchor": "left"}}, "position": {"distance": 0.66, "offset": -40}}], "z": 2},
                     ]}

  assert new_jointjs == jointjs_truth