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
import copy
from pprint import pprint
import pytest

from idaes.ui.flowsheet import FlowsheetSerializer, FlowsheetDiff, validate_flowsheet

# === Sample data ===

base_model = {
    "model": {
        "id": "Model1",
        "unit_models": {},
        "arcs": {}
    },
    "cells": {}
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
                "image": unit_types[unit_num] + ".svg"
            }
        m["arcs"] = {}
        if n > 1:
            for arc_num in range(n):
                unit_num = arc_num
                m["arcs"][f"A{arc_num}"] = {
                    "source": f"U{unit_num}",
                    "dest": f"U{(unit_num + 1) % n}",
                    "label": f"stream {arc_num}"
                }
        # add minimal cells for each unit model and arc
        c = model["cells"] = []
        for key, value in m["unit_models"].items():
            c.append({"id": key,
                      "attrs": {"image": {"xlinkHref": "image.svg"}, "root": {"title": "TITLE"}}})
        for key, value in m["arcs"].items():
            c.append({"id": key, "source": {"id": value["source"]}, "target": {"id": value["dest"]},
                      "labels": [{"attrs": {"text": {"text": "LABEL"}}}]})
        # done
        models[n] = model
    return models

# === Tests ===


@pytest.mark.unit
def test_merge(models):
    """Test the FlowsheetDiff output from the .merge() function.
    """
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



# @pytest.mark.unit
# def test_compare_models_edge_cases():
#     # Both empty models
#     existing_model = {"model": {"id": 1, "unit_models": {}, "arcs": {}}}
#
#     new_model = {"model": {"id": 2, "unit_models": {}, "arcs": {}}}
#
#     diff_model, out_json = fc.compare_models(existing_model, new_model)
#     assert diff_model == {}
#
#     # Existing model is empty
#     existing_model = {"model": {"id": 1, "unit_models": {}, "arcs": {}}}
#
#     new_model = {
#         "model": {
#             "id": 0,
#             "unit_models": {
#                 "M101": {"type": "mixer", "image": "mixer.svg"},
#                 "H101": {"type": "heater", "image": "heater_2.svg"},
#             },
#             "arcs": {"s03": {"source": "M101", "dest": "H101", "label": "hello"}},
#         }
#     }
#
#     diff_model, out_json = fc.compare_models(existing_model, new_model)
#
#     # Since the existing model is empty we return an empty diff_model and
#     # return new_model as the output json
#     diff_model_truth = {}
#     assert diff_model == diff_model_truth
#     assert out_json == new_model
#
#     # New model is empty
#     existing_model = {
#         "model": {
#             "id": 0,
#             "unit_models": {
#                 "M101": {"type": "mixer", "image": "mixer.svg"},
#                 "H101": {"type": "heater", "image": "heater_2.svg"},
#             },
#             "arcs": {"s03": {"source": "M101", "dest": "H101", "label": "hello"}},
#         }
#     }
#
#     new_model = {"model": {"id": 1, "unit_models": {}, "arcs": {}}}
#
#     diff_model, out_json = fc.compare_models(existing_model, new_model)
#
#     assert diff_model == {
#         "M101": {
#             "type": "mixer",
#             "image": "mixer.svg",
#             "action": Action.REMOVE.value,
#             "class": "unit model",
#         },
#         "H101": {
#             "type": "heater",
#             "image": "heater_2.svg",
#             "action": Action.REMOVE.value,
#             "class": "unit model",
#         },
#         "s03": {
#             "source": "M101",
#             "dest": "H101",
#             "label": "hello",
#             "action": Action.REMOVE.value,
#             "class": "arc",
#         },
#     }
#
#     # The models are the same
#     model = {
#         "model": {
#             "id": 0,
#             "unit_models": {
#                 "M101": {"type": "mixer", "image": "mixer.svg"},
#                 "H101": {"type": "heater", "image": "heater_2.svg"},
#             },
#             "arcs": {"s03": {"source": "M101", "dest": "H101", "label": "hello"}},
#         }
#     }
#
#     diff_model, out_json = fc.compare_models(model, model)
#     assert diff_model == {}
#
#
# @pytest.mark.unit
# def test_compare_models_errors():
#     # Both empty
#     existing_model = {}
#     new_model = {}
#
#     with pytest.raises(KeyError):
#         fc.compare_models(existing_model, new_model)
#
#     # New model is missing the id
#     existing_model = {
#         "model": {
#             "id": 0,
#             "unit_models": {
#                 "M101": {"type": "mixer", "image": "mixer.svg"},
#                 "H101": {"type": "heater", "image": "heater_2.svg"},
#             },
#             "arcs": {"s03": {"source": "M101", "dest": "H101", "label": "hello"}},
#         }
#     }
#     new_model = {
#         "model": {
#             "unit_models": {
#                 "M101": {"type": "mixer", "image": "mixer.svg"},
#                 "H101": {"type": "heater", "image": "heater_2.svg"},
#             },
#             "arcs": {"s03": {"source": "M101", "dest": "H101", "label": "hello"}},
#         }
#     }
#
#     with pytest.raises(ValidationError):
#         fc.compare_models(existing_model, new_model)
#
#     # Existing model is the wrong format
#     existing_model = {
#         "model": {
#             "id": 0,
#             "arcs": {"s03": {"source": "M101", "dest": "H101", "label": "hello"}},
#         }
#     }
#     new_model = {
#         "model": {
#             "id": 0,
#             "unit_models": {
#                 "M101": {"type": "mixer", "image": "mixer.svg"},
#                 "H101": {"type": "heater", "image": "heater_2.svg"},
#             },
#             "arcs": {"s03": {"source": "M101", "dest": "H101", "label": "hello"}},
#         }
#     }
#
#     with pytest.raises(ValidationError):
#         fc.compare_models(existing_model, new_model)
#
#
# @pytest.mark.unit
# def test_model_jointjs_conversion():
#     # Test unit model addition
#     original_jointjs = {
#         "cells": [
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 100, "y": 100},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "M101",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "mixer.svg"},
#                     "label": {"text": "M101"},
#                     "root": {"title": "mixer"},
#                 },
#             },
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 200, "y": 200},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "H101",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "heater_2.svg"},
#                     "label": {"text": "H101"},
#                     "root": {"title": "heater"},
#                 },
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {
#                     "anchor": {
#                         "name": "right",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "H101",
#                 },
#                 "target": {
#                     "anchor": {
#                         "name": "topLeft",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "R101",
#                 },
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s04",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "hello",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {
#                     "anchor": {
#                         "name": "bottomRight",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "R101",
#                 },
#                 "target": {
#                     "anchor": {
#                         "name": "left",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "F101",
#                 },
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s05",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "world",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#         ]
#     }
#
#     diff_model = {
#         "F222": {
#             "type": "flash",
#             "image": "flash.svg",
#             "action": Action.ADD.value,
#             "class": "unit model",
#         }
#     }
#
#     new_jointjs = fc.model_jointjs_conversion(diff_model, original_jointjs)
#     jointjs_truth = {
#         "cells": [
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 100, "y": 100},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "M101",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "mixer.svg"},
#                     "label": {"text": "M101"},
#                     "root": {"title": "mixer"},
#                 },
#             },
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 200, "y": 200},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "H101",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "heater_2.svg"},
#                     "label": {"text": "H101"},
#                     "root": {"title": "heater"},
#                 },
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {
#                     "anchor": {
#                         "name": "right",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "H101",
#                 },
#                 "target": {
#                     "anchor": {
#                         "name": "topLeft",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "R101",
#                 },
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s04",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "hello",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {
#                     "anchor": {
#                         "name": "bottomRight",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "R101",
#                 },
#                 "target": {
#                     "anchor": {
#                         "name": "left",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "F101",
#                 },
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s05",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "world",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 150, "y": 150},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "F222",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "flash.svg"},
#                     "label": {"text": "F222"},
#                     "root": {"title": "flash"},
#                 },
#             },
#         ]
#     }
#
#     assert new_jointjs == jointjs_truth
#
#     # Test unit model removal
#     original_jointjs = {
#         "cells": [
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 100, "y": 100},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "M101",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "mixer.svg"},
#                     "label": {"text": "M101"},
#                     "root": {"title": "mixer"},
#                 },
#             },
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 200, "y": 200},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "H101",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "heater_2.svg"},
#                     "label": {"text": "H101"},
#                     "root": {"title": "heater"},
#                 },
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {
#                     "anchor": {
#                         "name": "right",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "H101",
#                 },
#                 "target": {
#                     "anchor": {
#                         "name": "topLeft",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "R101",
#                 },
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s04",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "hello",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {
#                     "anchor": {
#                         "name": "bottomRight",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "R101",
#                 },
#                 "target": {
#                     "anchor": {
#                         "name": "left",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "F101",
#                 },
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s05",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "world",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#         ]
#     }
#
#     diff_model = {
#         "M101": {
#             "type": "mixer",
#             "image": "mixer.svg",
#             "action": Action.REMOVE.value,
#             "class": "unit model",
#         }
#     }
#
#     new_jointjs = fc.model_jointjs_conversion(diff_model, original_jointjs)
#     jointjs_truth = {
#         "cells": [
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 200, "y": 200},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "H101",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "heater_2.svg"},
#                     "label": {"text": "H101"},
#                     "root": {"title": "heater"},
#                 },
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {
#                     "anchor": {
#                         "name": "right",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "H101",
#                 },
#                 "target": {
#                     "anchor": {
#                         "name": "topLeft",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "R101",
#                 },
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s04",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "hello",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {
#                     "anchor": {
#                         "name": "bottomRight",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "R101",
#                 },
#                 "target": {
#                     "anchor": {
#                         "name": "left",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "F101",
#                 },
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s05",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "world",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#         ]
#     }
#
#     assert new_jointjs == jointjs_truth
#
#     # Test arc addition
#     original_jointjs = {
#         "cells": [
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 100, "y": 100},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "M101",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "mixer.svg"},
#                     "label": {"text": "M101"},
#                     "root": {"title": "mixer"},
#                 },
#             },
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 200, "y": 200},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "H101",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "heater_2.svg"},
#                     "label": {"text": "H101"},
#                     "root": {"title": "heater"},
#                 },
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {
#                     "anchor": {
#                         "name": "right",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "H101",
#                 },
#                 "target": {
#                     "anchor": {
#                         "name": "topLeft",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "R101",
#                 },
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s04",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "hello",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {
#                     "anchor": {
#                         "name": "bottomRight",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "R101",
#                 },
#                 "target": {
#                     "anchor": {
#                         "name": "left",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "F101",
#                 },
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s05",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "world",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#         ]
#     }
#
#     diff_model = {
#         "s01": {
#             "source": "M101",
#             "dest": "F111",
#             "label": "foo",
#             "action": Action.ADD.value,
#             "class": "arc",
#         }
#     }
#
#     new_jointjs = fc.model_jointjs_conversion(diff_model, original_jointjs)
#     jointjs_truth = {
#         "cells": [
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 100, "y": 100},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "M101",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "mixer.svg"},
#                     "label": {"text": "M101"},
#                     "root": {"title": "mixer"},
#                 },
#             },
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 200, "y": 200},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "H101",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "heater_2.svg"},
#                     "label": {"text": "H101"},
#                     "root": {"title": "heater"},
#                 },
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {
#                     "anchor": {
#                         "name": "right",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "H101",
#                 },
#                 "target": {
#                     "anchor": {
#                         "name": "topLeft",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "R101",
#                 },
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s04",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "hello",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {
#                     "anchor": {
#                         "name": "bottomRight",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "R101",
#                 },
#                 "target": {
#                     "anchor": {
#                         "name": "left",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "F101",
#                 },
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s05",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "world",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {"id": "M101"},
#                 "target": {"id": "F111"},
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s01",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "foo",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#         ]
#     }
#
#     assert new_jointjs == jointjs_truth
#
#     # Test arc change
#     original_jointjs = {
#         "cells": [
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 100, "y": 100},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "M101",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "mixer.svg"},
#                     "label": {"text": "M101"},
#                     "root": {"title": "mixer"},
#                 },
#             },
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 200, "y": 200},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "H101",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "heater_2.svg"},
#                     "label": {"text": "H101"},
#                     "root": {"title": "heater"},
#                 },
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {
#                     "anchor": {
#                         "name": "right",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "H101",
#                 },
#                 "target": {
#                     "anchor": {
#                         "name": "topLeft",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "R101",
#                 },
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s04",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "hello",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {
#                     "anchor": {
#                         "name": "bottomRight",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "R101",
#                 },
#                 "target": {
#                     "anchor": {
#                         "name": "left",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "F101",
#                 },
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s05",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "world",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#         ]
#     }
#
#     diff_model = {
#         "s04": {
#             "source": "M101",
#             "dest": "F111",
#             "label": "asdf",
#             "action": Action.CHANGE.value,
#             "class": "arc",
#         }
#     }
#
#     new_jointjs = fc.model_jointjs_conversion(diff_model, original_jointjs)
#     jointjs_truth = {
#         "cells": [
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 100, "y": 100},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "M101",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "mixer.svg"},
#                     "label": {"text": "M101"},
#                     "root": {"title": "mixer"},
#                 },
#             },
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 200, "y": 200},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "H101",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "heater_2.svg"},
#                     "label": {"text": "H101"},
#                     "root": {"title": "heater"},
#                 },
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {
#                     "anchor": {
#                         "name": "bottomRight",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "R101",
#                 },
#                 "target": {
#                     "anchor": {
#                         "name": "left",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "F101",
#                 },
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s05",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "world",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {
#                     "anchor": {
#                         "name": "right",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "M101",
#                 },
#                 "target": {
#                     "anchor": {
#                         "name": "topLeft",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "F111",
#                 },
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s04",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "asdf",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#         ]
#     }
#
#     assert new_jointjs == jointjs_truth
#
#     # Test arc removal
#     original_jointjs = {
#         "cells": [
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 100, "y": 100},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "M101",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "mixer.svg"},
#                     "label": {"text": "M101"},
#                     "root": {"title": "mixer"},
#                 },
#             },
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 200, "y": 200},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "H101",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "heater_2.svg"},
#                     "label": {"text": "H101"},
#                     "root": {"title": "heater"},
#                 },
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {
#                     "anchor": {
#                         "name": "right",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "H101",
#                 },
#                 "target": {
#                     "anchor": {
#                         "name": "topLeft",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "R101",
#                 },
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s04",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "hello",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {
#                     "anchor": {
#                         "name": "bottomRight",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "R101",
#                 },
#                 "target": {
#                     "anchor": {
#                         "name": "left",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "F101",
#                 },
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s05",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "world",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#         ]
#     }
#
#     diff_model = {
#         "s04": {
#             "source": "M101",
#             "dest": "F111",
#             "label": "asdf",
#             "action": Action.REMOVE.value,
#             "class": "arc",
#         }
#     }
#
#     new_jointjs = fc.model_jointjs_conversion(diff_model, original_jointjs)
#     jointjs_truth = {
#         "cells": [
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 100, "y": 100},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "M101",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "mixer.svg"},
#                     "label": {"text": "M101"},
#                     "root": {"title": "mixer"},
#                 },
#             },
#             {
#                 "type": "standard.Image",
#                 "position": {"x": 200, "y": 200},
#                 "size": {"width": 50, "height": 50},
#                 "angle": 0,
#                 "id": "H101",
#                 "z": 1,
#                 "attrs": {
#                     "image": {"xlinkHref": "heater_2.svg"},
#                     "label": {"text": "H101"},
#                     "root": {"title": "heater"},
#                 },
#             },
#             {
#                 "type": "standard.Link",
#                 "source": {
#                     "anchor": {
#                         "name": "bottomRight",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "R101",
#                 },
#                 "target": {
#                     "anchor": {
#                         "name": "left",
#                         "args": {"rotate": "false", "padding": 0},
#                     },
#                     "id": "F101",
#                 },
#                 "router": {"name": "orthogonal", "padding": 10},
#                 "connector": {
#                     "name": "normal",
#                     "attrs": {"line": {"stroke": "#5c9adb"}},
#                 },
#                 "id": "s05",
#                 "labels": [
#                     {
#                         "attrs": {
#                             "rect": {
#                                 "fill": "#d7dce0",
#                                 "stroke": "#FFFFFF",
#                                 "stroke-width": 1,
#                             },
#                             "text": {
#                                 "text": "world",
#                                 "fill": "black",
#                                 "text-anchor": "left",
#                             },
#                         },
#                         "position": {"distance": 0.66, "offset": -40},
#                     }
#                 ],
#                 "z": 2,
#             },
#         ]
#     }
#
#     assert new_jointjs == jointjs_truth
