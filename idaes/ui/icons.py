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
from typing import Dict


class UnitModelIcon:
    """Represents icon display information for a given unit model.
    """

    _link_positions_map = {}  # cache 'built' link positions

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
            unit_model = self.DEFAULT
        self._model = unit_model
        try:
            self._info = self._mapping[unit_model]
        except KeyError:
            if not default:
                raise
            self._info = self._mapping[self.DEFAULT]
        self._pos = self._build_link_positions()

    @property
    def icon(self) -> str:
        """Get the name of the icon.
        """
        return self._info[0]

    @property
    def link_positions(self) -> Dict:
        """Get the link positions.

        Example result::

            {
                "groups": {
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
                "items": [{"group": "in", "id": "in"}, {"group": "out", "id": "out"}],
            }

        Returns:
            The link position (see example result)
        """
        return self._pos

    def _build_link_positions(self) -> Dict:
        """Fill in boilerplate based on raw info and place built value in class cache.
        Side-effects: set self._pos and add entry to class' _link_positions_map
        """
        # look in cache, return if found
        if self._model in self._link_positions_map:
            return self._link_positions_map[self._model]

        # build link positions from info
        groups, items = {}, []
        for position in self._info[1:]:
            group, name, (x, y, dx, dy) = position
            if group not in groups:
                groups[group] = {}
            groups[group] = {
                "position": {
                    "name": name,
                    "args": {"x": x, "y": y, "dx": dx, "dy": dy},
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
            }
            items.append({"group": group, "id": group})

        # set new link positions attr and place in cache
        positions = {"groups": groups, "items": items}
        self._link_positions_map[self._model] = positions
        return positions

    # === Data ===

    # Name is unit name, value is the information for its icon
    # Value is a tuple: icon image, and one or more position tuples: (group [in/out], name [side], (x, y, dx, dy))
    # Notes for updating:
    #  - Use 'cstr' as your template for new entries
    #  - Do not remove in/out entries in existing entries, or arcs won't connect
    _mapping = {
        "cstr": (
            "reactor_c.svg",
            ("in", "left", (15, 0, 1, 1)),
            ("out", "left", (48, 45, 1, 1)),
        ),
        "flash": (
            "flash.svg",
            ("bottom", "bottom", (25, 50, 1, 1)),
            ("in", "left", (8, 25, 1, 1)),
            ("top", "top", (25, 0, 1, 1)),
        ),
        "gibbs_reactor": (
            "reactor_g.svg",
            ("in", "left", (5, 10, 1, 1)),
            ("out", "left", (45, 45, 1, 1)),
        ),
        "heat_exchanger": (
            "heat_exchanger_1.svg",
            ("in", "left", (2, 25, 1, 1)),
            ("out", "left", (48, 25, 1, 1)),
        ),
        "heater": (
            "heater_2.svg",
            ("in", "left", (6, 25, 1, 1)),
            ("out", "left", (43, 25, 1, 1)),
        ),
        "heat_exchanger_1D": (
            "heat_exchanger_1.svg",
            ("in", "left", (15, 0, 1, 1)),
            ("out", "left", (48, 45, 1, 1)),
        ),
        "mixer": (
            "mixer.svg",
            ("in", "left", (2, 25, 1, 1)),
            ("out", "left", (48, 25, 1, 1)),
        ),
        "plug_flow_reactor": (
            "reactor_pfr.svg",
            ("in", "left", (15, 0, 1, 1)),
            ("out", "left", (48, 45, 1, 1)),
        ),
        "pressure_changer": (
            "compressor.svg",
            ("in", "left", (2, 25, 1, 1)),
            ("out", "left", (48, 25, 1, 1)),
        ),
        "separator": (
            "splitter.svg",
            ("in", "left", (2, 25, 1, 1)),
            ("out", "right", (48, 25, 1, 1)),
        ),
        "stoichiometric_reactor": (
            "reactor_s.svg",
            ("in", "left", (5, 10, 1, 1)),
            ("out", "left", (45, 45, 1, 1)),
        ),
        "equilibrium_reactor": (
            "reactor_e.svg",
            ("in", "left", (5, 10, 1, 1)),
            ("out", "left", (45, 45, 1, 1)),
        ),
        "feed": ("feed.svg", ("out", "left", (48, 25, 1, 1))),
        "product": ("product.svg", ("in", "left", (2, 25, 1, 1))),
        "feed_flash": (
            "feed.svg",
            ("in", "left", (25, 0, 1, 1)),
            ("out", "left", (25, 50, 1, 1)),
        ),
        "statejunction": (
            "NONE",
            ("in", "left", (15, 0, 1, 1)),
            ("out", "left", (48, 45, 1, 1)),
        ),
        "translator": (
            "NONE",
            ("in", "left", (15, 0, 1, 1)),
            ("out", "left", (48, 45, 1, 1)),
        ),
        "packed_column": (
            "packed_column_1.svg",
            ("in", "left", (48, 10, 1, 1)),
            ("out", "left", (48, 40, 1, 1)),
        ),
        "tray_column": (
            "tray_column_1.svg",
            ("in", "left", (48, 10, 1, 1)),
            ("out", "left", (48, 40, 1, 1)),
        ),
        "default": (
            "default.svg",
            ("in", "left", (2, 0, 1, 1)),
            ("out", "left", (48, 50, 1, 1)),
        ),
    }
