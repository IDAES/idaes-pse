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
# stdlib
import copy
import json
from typing import Dict, List, Tuple


class FlowsheetDiff:
    """Compute a flowsheet model 'diff' and use that to compute an updated layout.

    Flowsheets are serialized by the core library. See :func:`validate()` for the required format.

    Example usage::

        diff =  FlowsheetDiff(saved_flowsheet, new_flowsheet)
        merged_flowsheet = {"model": new_flowsheet["model"], "cells": diff.layout}
    """
    def __init__(self, old_flowsheet: Dict, new_flowsheet: Dict, validate=True):
        """Construct with old and new flowsheet, and compute the diff and new layout

        Args:
            old_flowsheet: Old flowsheet value
            new_flowsheet: New flowsheet value
            validate: If True, validate each flowsheet before doing the diff
        Raises:
            ValueError if either flowsheet is invalid
        """
        if validate:
            ok, why = validate_flowsheet(old_flowsheet)
            if not ok:
                raise ValueError(f"invalid old_flowsheet value: {why}")
            ok, why = validate_flowsheet(new_flowsheet)
            if not ok:
                raise ValueError(f"invalid new_flowsheet value: {why}")
        self._old, self._new = old_flowsheet, new_flowsheet
        self._len = None
        self._diff = self._compute_diff()
        self._layout = self._compute_layout() if self._diff else None

    def merged(self, do_copy: bool = False) -> Dict:
        """Return merged flowsheet.

        If the diff is empty, this will be the 'old_flowsheet' (or a copy) passed to the constructor.

        Args:
            do_copy: If True, return a copy so modifying the returned dict won't affect input arguments.

        Returns:
            Merged dict with structure as described in :func:`validate_flowsheet()`.
        """
        if bool(self):
            result = {
                "model": copy.deepcopy(self._new["model"]) if do_copy else self._new["model"],
                "cells": self._layout
            }
        else:
            # diff is empty, return 'old' object
            result = copy.deepcopy(self._old) if do_copy else self._old
        return result

    def __len__(self):
        if self._len is None:
            return 0
        return self._len

    def __bool__(self):
        return len(self) > 0

    def __str__(self):
        return json.dumps(self._diff, indent=2)

    def _compute_diff(self) -> Dict:
        diff = {"add": {}, "remove": {}, "change": {}}
        old_model, new_model = self._old["model"], self._new["model"]
        n = 0
        for cls in "unit_models", "arcs":
            for k in diff.keys():
                diff[k][cls] = {}
            old_data, new_data = old_model[cls], new_model[cls]
            # Add/change
            for key in new_data:
                if key not in old_data:
                    diff["add"][cls][key] = copy.deepcopy(new_data[key])
                    n += 1
                elif old_data[key] != new_data[key]:
                    diff["change"][cls][key] = copy.deepcopy(new_data[key])
                    n += 1
            # Remove
            for key in old_data:
                if key not in new_data:
                    diff["remove"][cls][key] = True
                    n += 1
        self._len = n
        return diff

    def _compute_layout(self) -> List:
        """Based on the old and new flowsheets provided in the constructor, compute the layout
        that should be placed in the "cells" of the new flowsheet.

        Returns:
            New layout to put into the "cells" of the new flowsheet dict
        """
        layout, x, y = [], 50, 50
        # Add
        for cls in self._diff["add"]:
            for id_ in self._diff["add"][cls]:
                values = self._diff["add"][cls][id_]
                if cls == "arcs":
                    new_item = self._new_arc(id_, values)
                else:
                    new_item = self._new_unit_model(id_, values, x, y)
                    x, y = x + 100, y + 100
                layout.append(new_item)
        # Change, remove, and simply copy
        for item in self._old["cells"]:
            id_ = item["id"]
            cls = "arcs" if "source" in item else "unit_models"
            if id_ in self._diff["remove"][cls]:
                continue
            elif id_ in self._diff["change"][cls]:
                values = self._diff["change"][cls][id_]
                if cls == "arcs":
                    new_item = self._update_arc(item, values)
                else:
                    new_item = self._update_unit_model(item, values)
                layout.append(new_item)
            else:
                layout.append(copy.deepcopy(item))
        return layout

    @staticmethod
    def _new_arc(name, values):
        return {
            "type": "standard.Link",
            "source": {"id": values["source"]},
            "target": {"id": values["dest"]},
            "router": {"name": "orthogonal", "padding": 10},
            "connector": {"name": "normal", "attrs": {"line": {"stroke": "#5c9adb"}}},
            "id": name,
            "labels": [
                {
                    "attrs": {
                        "rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", "stroke-width": 1},
                        "text": {
                            "text": values["label"],
                            "fill": "black",
                            "text-anchor": "left",
                        },
                    },
                    "position": {"distance": 0.66, "offset": -40},
                }
            ],
            "z": 2,
        }

    @staticmethod
    def _update_arc(item, values):
        new_item = copy.deepcopy(item)
        new_item["source"]["id"] = values["source"]
        new_item["target"]["id"] = values["dest"]
        new_item["labels"][0]["attrs"]["text"]["text"] = values["label"]
        return new_item

    @staticmethod
    def _new_unit_model(name, values, x, y):
        return {
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
            },
        }

    @staticmethod
    def _update_unit_model(item, values):
        new_item = copy.deepcopy(item)
        new_item["attrs"]["image"]["xlinkHref"] = values["image"]
        new_item["attrs"]["root"]["title"] = values["type"]
        return new_item


def validate_flowsheet(fs: Dict) -> Tuple[bool, str]:
    """Validate a flowsheet.

    Expected format is below.

    .. code-block:: json

        {
            "model": {
                "id": "<model name>",
                "unit_models": {
                    "<component name>": {
                        "image": "<image name>",
                        "type": "<component type name>",
                        ..<more values from model>..
                    },
                    ...
                },
                "arcs": {
                    "<arc name>": {
                        "source": "<component name>",
                        "dest": "<component name>",
                        "label": "<label text>"
                    },
                    ...
                }
            },
            "cells": [
                {"id": "<component_name>", ..<values used by JointJS>..},
                ...
            ]
        }

    Args:
        fs: Flowsheet to validate

    Return:
        Tuple of (True, "") for OK, and (False, "<message>") for failure
    """
    # very quick and dirty validation, but it does make for nice clean error messages
    for key in "model", "cells":
        if key not in fs:
            return False, f"Missing top-level key '{key}'"
    model, component_ids = fs["model"], set()
    for key2 in "id", "unit_models", "arcs":
        if key2 not in model:
            return False, f"The flowsheet model is missing key '{key2}'"
        if key2 == "unit_models":
            for ckey, cval in model[key2].items():
                for key3 in "image", "type":
                    if key3 not in cval:
                        return False, f"Unit model '{ckey}' is missing key '{key3}'"
                component_ids.add(ckey)
        elif key2 == "arcs":
            for akey, aval in model[key2].items():
                for key3 in "source", "dest", "label":
                    if key3 not in aval:
                        return False, f"Arc '{akey}' is missing key '{key3}'"
                component_ids.add(akey)
    cells = fs["cells"]
    for i, cell in enumerate(cells):
        if "id" not in cell:
            return False, f"Cell #{i + 1} is missing key 'id'"
        cell_id = cell["id"]
        if cell_id not in component_ids:
            return False, f"Cell id '{cell_id}' not found in unit models or arcs"
    return True, ""
