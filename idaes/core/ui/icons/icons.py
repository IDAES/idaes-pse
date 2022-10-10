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
#
# Author: Abdelrahman Elbashandy
#################################################################################
import os
import json
from typing import Dict


class UnitModelIcon:
    """Represents icon display information for a given unit model."""

    #: Name of default unit_model to use
    DEFAULT = "default"

    def __init__(self, unit_model: str = None, default: str = DEFAULT):
        """Construct with a given unit model type name.

        Args:
            unit_model: Name of the unit model type. If not given, use value in attribute `DEFAULT`
            default: If the given `unit_model` is not found, use this one; but then if this is falsy raise a KeyError

        Raises:
            KeyError if unit_model name is not found and `default` arg is falsy (e.g. None or "")
        """
        if unit_model is None:
            unit_model = default
        self._model = unit_model

        # Loading the Unit Models mappings
        dir_path = os.path.dirname(os.path.realpath(__file__))
        mappings_file = os.path.join(dir_path, os.pardir, "mappings", "mappings.json")
        with open(mappings_file, "r") as mappings_f:
            self._mapping = json.load(mappings_f)
        self._model_details = self._get_mapping(unit_model, default)
        self._pos = self._build_link_positions()

    def _get_mapping(self, unit_model, default):
        """Find the correct mapping for the given unit_model name."""
        if unit_model in self._mapping:
            return self._mapping[unit_model]

        # Couldn't find unit model and using default model instead
        if not default or default not in self._mapping:
            raise ValueError(
                "Specified unit model doesn't exist, and the default model is not set."
            )
        return self._mapping[default]

    @property
    def icon(self) -> str:
        """Get the name of the icon."""
        # return self._info[0]
        return self._model_details["image"]

    @property
    def link_positions(self) -> Dict:
        """Get the link positions.

        Example result::

            {
                "port_groups": {
                    "in": {
                        "position": {
                            "name": "left",
                            "args": {"x": 15, "y": 0, "dx": 1, "dy": 1},
                        },
                        "attrs": {
                            "rect": {
                                "stroke": "#000000",
                                "stroke-width": 0,
                                "width": 0,
                                "height": 0,
                            }
                        },
                        "markup": "<g><rect/></g>",
                    },
                    "out": {
                        "position": {
                            "name": "left",
                            "args": {"x": 48, "y": 45, "dx": 1, "dy": 1},
                        },
                        "attrs": {
                            "rect": {
                                "stroke": "#000000",
                                "stroke-width": 0,
                                "width": 0,
                                "height": 0,
                            }
                        },
                        "markup": "<g><rect/></g>",
                    },
                },
                "items": []
            }

        Returns:
            The link position (see example result)
        """
        return self._pos

    @property
    def routing_config(self) -> Dict:
        """Get the Unit model routing config to be used to add jointjs vertices
        for layout control within the created graph.

        Example result::

            {
                "in": {
                    "gap": {
                        "direction": "left",
                        "distance": 20
                    }
                }
            }

        Returns:
            The routing configuration (see example result)
        """
        return self._model_details["routing_config"]

    def _build_link_positions(self) -> Dict:
        """Fill in boilerplate based on raw info and place built value in class cache."""
        # build link positions from info
        groups, items = {}, []
        for group_name, group_config in self._model_details["port_groups"].items():
            groups[group_name] = group_config
            groups[group_name].update(
                {
                    "attrs": {
                        "rect": {
                            "stroke": "#000000",
                            "stroke-width": 0,
                            "width": 0,
                            "height": 0,
                        }
                    },
                    "markup": "<g><rect/></g>",
                }
            )

        # set new link positions attr and place in cache
        positions = {"groups": groups, "items": []}
        return positions

    # === Data ===

    # Name is unit name, value is the information for its icon
    # Value is a tuple: icon image, and one or more position tuples: (group [in/out], name [side], (x, y, dx, dy))
    # Notes for updating:
    #  - Use 'cstr' as your template for new entries
    #  - Do not remove in/out entries in existing entries, or arcs won't connect
    # TODO: Move this mapping to its own directory/files.
    _mapping = {}
