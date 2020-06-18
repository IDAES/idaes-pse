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
from enum import Enum
from jsonschema import validate


class UnknownModelDiffException():
    pass


class Action(Enum):
    REMOVE = 1
    ADD = 2
    CHANGE = 3


def compare_models(existing_model, new_model):
    """
    Compares two models that are in this format

    .. code-block:: json

        {
            "model": {
                "id": "id", 
                "unit_models": {
                    "M101": {
                        "image": "mixer.svg", 
                        "type": "mixer"
                    }
                },
                "arcs": {
                    "s03": {
                        "source": "M101", 
                        "dest": "H101", 
                        "label": "molar flow ('Vap', 'hydrogen') 0.5"
                    }
                }
            },
            "cells": [{ "--jointjs code--": "--jointjs code--" }]
        }

    :param existing_model: The current model to compare against
    :param new_model: The new model that has changes
    :return: A diff between the models in this format: 

    .. code-block:: json

        {
            "M111": {
                "type": "flash", 
                "image": "flash.svg", 
                "action": "1", 
                "class": "unit model"
            },
            "s03": {
                "source": "M101", 
                "dest": "F102", 
                "label": "Hello World", 
                "action": "3", 
                "class": "arc"
            }, 
            "s11": {
                "source": "H101", 
                "dest": "F102", 
                "label": "molar flow ('Vap', 'hydrogen') 0.5", 
                "action": "2", 
                "class": "arc"
            }
        }

    """
    diff_model = {}
    model_schema = {
        "type" : "object",
        "properties" : {
            "id" : {"type" : ["number", "string"]},
            "unit_models" : {
                "type" : "object",
            },
            "arcs" : {
                "type" : "object",
            }
        },
        "required" : ["id", "unit_models", "arcs"]
    }

    # Copy the new model into the out_json's model key
    out_json = dict(existing_model)
    try:
        out_json["model"] = new_model["model"]
    except KeyError as error:
        msg = "Unable to find 'model' section of new model json"
        raise KeyError(msg)

    validate(instance=existing_model["model"], schema=model_schema)
    validate(instance=new_model["model"], schema=model_schema)

    try:
        existing_model = existing_model["model"]
    except KeyError as error:
        msg = "Unable to find 'model' section of existing model json"
        raise KeyError(msg)
    # If the existing model is empty then return an empty diff_model 
    # and the full json from the new model
    if existing_model["unit_models"] == {} and \
        existing_model["arcs"] == {}:
        return {}, new_model

    # If the models are the same return an empty diff_model and the full json from
    # the new model. This will happen when the user moves something in the 
    # graph but doesn't change the actual idaes model
    if existing_model == new_model["model"]:
        return {}, new_model

    try:
        new_model = new_model["model"]
    except KeyError as error:
        msg = "Unable to find 'model' section of new model json"
        raise KeyError(msg)

    unit_model_schema = {
        "type" : "object",
        "properties" : {
            "type" : {"type" : "string"},
            "image" : {"type" : "string"}
        },
    }

    # Check for new or changed unit models
    for item, value in new_model["unit_models"].items():
        validate(instance=value, schema=unit_model_schema)
        if item not in existing_model["unit_models"]:
            diff_model[item] = value
            diff_model[item]["action"] = Action.ADD.value
            diff_model[item]["class"] = "unit model"
        elif existing_model["unit_models"][item] != value:
            diff_model[item] = value
            diff_model[item]["action"] = Action.CHANGE.value
            diff_model[item]["class"] = "unit model"
        elif existing_model["unit_models"][item] == value:
            pass
        else:
            msg = ("Unknown diff between new model and existing_model. "
                   "Key: " + str(item) + " new model value: " + str(value) + ", "
                   "existing model value: " + str(existing_model["unit_models"][item]))
            raise UnknownModelDiffException(msg)

    # Check for unit models that have been removed
    for item, value in existing_model["unit_models"].items():
        validate(instance=value, schema=unit_model_schema)
        if item not in new_model["unit_models"]:
            diff_model[item] = value
            diff_model[item]["action"] = Action.REMOVE.value
            diff_model[item]["class"] = "unit model"

    arc_schema = {
        "type" : "object",
        "properties" : {
            "source" : {"type" : "string"},
            "dest" : {"type" : "string"},
            "label" : {"type" : "string"}
        }
    }

    # Check for new or changed arcs
    for item, value in new_model["arcs"].items():
        validate(instance=value, schema=arc_schema)
        if item not in existing_model["arcs"]:
            diff_model[item] = value
            diff_model[item]["action"] = Action.ADD.value
            diff_model[item]["class"] = "arc"
        elif existing_model["arcs"][item] != value:
            diff_model[item] = value
            diff_model[item]["action"] = Action.CHANGE.value
            diff_model[item]["class"] = "arc"
        else:
            pass

    # Check for arcs that have been removed
    for item, value in existing_model["arcs"].items():
        validate(instance=value, schema=arc_schema)
        if item not in new_model["arcs"]:
            diff_model[item] = value
            diff_model[item]["action"] = Action.REMOVE.value
            diff_model[item]["class"] = "arc"


    return diff_model, out_json


def model_jointjs_conversion(diff_model, current_json):
    """
    Converts a model diff (output from compare_models) to jointjs
    :param diff_model: A diff between the models in this format:

    .. code-block:: json

        {
            "M111": {
                "type": "flash", 
                "image": "flash.svg", 
                "action": "1", 
                "class": "unit model"
            }, 
            "s03": {
                "source": "M101", 
                "dest": "F102", 
                "label": "Hello World", 
                "action": "3", 
                "class": "arc"
            }, 
            "s11": {
                "source": "H101", 
                "dest": "F102", 
                "label": "molar flow ('Vap', 'hydrogen') 0.5", 
                "action": "2", 
                "class": "arc"
            }
        }

    :param current_json: The json from jointjs and the current model. In this format:

    .. code-block:: json

        {
            "cells": [{ "--jointjs code--": "--jointjs code--" }],
            "model": {
                "id": "id", 
                "unit_models": {
                    "M101": {
                        "image": "mixer.svg", 
                        "type": "mixer"
                    }
                },
                "arcs": {
                    "s03": {
                        "source": "M101", 
                        "dest": "H101", 
                        "label": "molar flow ('Vap', 'hydrogen') 0.5"
                    }
                }
            }
        }

    """
    new_json = dict(current_json)
    x = 50
    y = 50
    for name, values in diff_model.items():
        x += 100
        y += 100
        found = False
        for item in current_json["cells"]:
            if name == item["id"]:
                # If the action is not removed then update the values
                # If it is removed then just remove it from the new_json
                new_json["cells"].remove(item)

                if values["action"] != Action.REMOVE.value:
                    if values["class"] == "arc":
                        item["id"] = name
                        item["source"]["id"] = values["source"]
                        item["target"]["id"] = values["dest"]
                        item["labels"][0]["attrs"]["text"]["text"] = values["label"]
                        new_json["cells"].append(item)
                    elif values["class"] == "unit model":
                        item["id"] = name
                        item["attrs"]["label"]["text"] = name
                        item["attrs"]["image"]["xlinkHref"] = values["image"]
                        item["attrs"]["root"]["title"] = values["type"]
                        new_json["cells"].append(item)
                    else:
                        # if the class isn't a unit model or an arc then throw an 
                        # exception because we don't know what this is
                        raise UnknownModelDiffException("Unknown model item class")
                found = True
                break

        if not found:
            if values["class"] == "arc":
                new_item = {
                    "type": "standard.Link",
                    "source": {"id": values["source"]},
                    "target": {"id": values["dest"]},
                    "router": {"name": "orthogonal", "padding": 10},
                    "connector": {"name": "normal", 
                                  "attrs": {"line": {"stroke": "#5c9adb"}}},
                    "id": name,
                    "labels": [{
                        "attrs": {
                            "rect": {"fill": "#d7dce0", 
                                     "stroke": "#FFFFFF", 
                                     'stroke-width': 1},
                            "text": {
                                "text": values["label"],
                                "fill": 'black',
                                'text-anchor': 'left',
                            },
                        },
                        "position": {
                            "distance": 0.66,
                            "offset": -40
                        },
                    }],
                    "z": 2
                }
            elif values["class"] == "unit model":
                new_item = {
                    "type": "standard.Image",
                    "id": name,
                    "position": {"x": x, "y": y},
                    "size": {"width": 50, "height": 50},
                    "angle": 0,
                    "z": 1,
                    "attrs": {
                        "image": {"xlinkHref": values["image"]},
                        "label": {"text": name},
                        "root": {"title": values["type"]},
                    }
                }
            else:
                # if the class isn't a unit model or an arc then throw an 
                # exception because we don't know what this is
                raise UnknownModelDiffException("Unknown model item class")

            new_json["cells"].append(new_item)

    return new_json
