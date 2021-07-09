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
Flowsheet-related classes and functions used by the UI.
"""
# stdlib
from collections import defaultdict
import copy
import json
import logging
import re
from typing import Dict, List, Tuple

# third-party
import pandas as pd
import numpy as np
from pyomo.environ import Block, value
from pyomo.network import Arc
from pyomo.network.port import Port

# package
from idaes import logger
from idaes.ui.icons import UnitModelIcon

_log = logger.getLogger(__name__)


class FileBaseNameExistsError(Exception):
    pass


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
                        "...": "more values..."
                    },
                    "...": "etc."
                },
                "arcs": {
                    "<arc name>": {
                        "source": "<component name>",
                        "dest": "<component name>",
                        "label": "<label text>"
                    },
                    "...": "etc."
                }
            },
            "cells": [
                {"id": "<component_name>", "...": "other values used by JointJS.."},
                "..."
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
    # Check if all cell id's are in the model
    cell_ids = set()
    for i, cell in enumerate(cells):
        if "id" not in cell:
            return False, f"Cell #{i + 1} is missing key 'id'"
        cell_id = cell["id"]
        if cell_id not in component_ids:
            return False, f"Cell id '{cell_id}' not found in unit models or arcs"
        cell_ids.add(cell_id)
    # Check if all model id's are in the cells
    if cell_ids != component_ids:
        missing = component_ids - cell_ids
        sfx = "s" if len(missing) > 1 else ""
        return (
            False,
            f"Component id{sfx} {missing} {'are' if sfx else 'is'} not in the layout cells",
        )
    return True, ""


class FlowsheetSerializer:
    """Serializes the flowsheet into one dict with two sections.

    The "model" section contains the id of the flowsheet and the
    unit models and arcs. This will be used to compare the model and convert
    to jointjs. The "cells" section is the jointjs readable code.
    See :py:func:`validate_flowsheet` for details on the format.
    """

    #: Regular expression identifying inlets by last component of ports' name
    INLET_REGEX = re.compile(
        r"^in_|^feed_|^inlet_|^in$|^feed$|^inlet$|_in$|_feed$|_inlet$", re.IGNORECASE
    )
    #: Regular expression identifying outlets by last component of ports' name
    OUTLET_REGEX = re.compile(
        r"^out_|^prod_|^outlet_|^out$|^prod$|^outlet$|_out$|_prod$|_outlet$",
        re.IGNORECASE,
    )

    def __init__(self, flowsheet, name: str, validate: bool = True):
        """Serialize input flowsheet with given name

        Args:
            flowsheet: The flowsheet to serialize
            name: The name of the flowsheet (also called its 'id' in some contexts)
            validate: If True, validate that the flowsheet is a reaonsable IDAES model, first
        Raises:
            ValueError if validation is on and flowsheet is found to be invalid
        """
        # validate, if so directed
        if validate:
            if not hasattr(flowsheet, "component_objects") or not callable(
                flowsheet.component_objects
            ):
                raise ValueError("Flowsheet missing function 'component_objects'")
            try:
                n_obj = 0
                for objtype in Arc, Block:
                    for obj in flowsheet.component_objects(objtype, descend_into=False):
                        n_obj += 1
                        if not hasattr(obj, "getname") or not callable(obj.getname):
                            raise ValueError(
                                "Flowsheet component missing function 'getname'"
                            )
            except ValueError:
                raise
            except Exception as err:
                raise ValueError(
                    f"Flowsheet 'component_objects' cannot be navigated: {err}"
                )
            if n_obj == 0:
                raise ValueError("Flowsheet has no Arcs or unit Blocks")
        # setup
        self.unit_models = {}  # {unit: {"name": unit.getname(), "type": str?}}
        self.arcs = {}  # {Arc.getname(): Arc}
        self.ports = {}  # {Port: parent_unit}
        self.edges = defaultdict(list)  # {name: {"source": unit, "dest": unit}}
        self.orphaned_ports = {}
        self.labels = {}
        self._stream_table_df = None
        self._out_json = {"model": {}}
        self._serialized_contents = defaultdict(dict)
        self._used_ports = set()
        self._known_endpoints = set()
        self._unit_name_used_count = defaultdict(lambda: 0)
        self._sig_figs = (
            5  # Defines the number of significant figures after the decimal place
        )
        self._logger = logger.getLogger(__name__)
        self.name = name
        self.flowsheet = flowsheet
        # serialize
        self._ingest_flowsheet()
        self._construct_output_json()

    def as_dict(self):
        return self._out_json

    def _ingest_flowsheet(self):
        # Stores information on the connectivity and components of the input flowsheet
        self._identify_arcs()

        # Identify and add unit models
        unit_models_map = self._identify_unit_models()
        for model, model_type in unit_models_map.items():
            self._add_unit_model_with_ports(model, model_type)

        self._construct_stream_labels()
        untouched_ports = self._map_edges()
        self._identify_implicit_feeds_and_products(untouched_ports)

    def _identify_arcs(self):
        # Identify the arcs and known endpoints and store them
        for component in self.flowsheet.component_objects(Arc, descend_into=False):
            self.arcs[component.getname()] = component
            self._known_endpoints.add(component.source.parent_block())
            self._known_endpoints.add(component.dest.parent_block())

    def _identify_unit_models(self) -> Dict:
        from idaes.core import UnitModelBlockData  # avoid circular import
        from idaes.core.property_base import PhysicalParameterBlock, StateBlock

        # Create a map of components to their unit type
        components = {}

        # Identify the unit models and ports and store them
        for component in self.flowsheet.component_objects(Block, descend_into=True):
            if isinstance(component, UnitModelBlockData):
                # List of components is the same as the provided one
                components[component] = self.get_unit_model_type(component)
            elif isinstance(component, PhysicalParameterBlock) or isinstance(
                component, StateBlock
            ):
                # skip physical parameter / state blocks
                pass
            else:
                # Find unit models nested within indexed blocks
                type_ = self.get_unit_model_type(component)
                for item in component.parent_component().values():
                    if isinstance(item, UnitModelBlockData):
                        # See if this unit is connected to an arc
                        is_connected = False
                        for arc in self.arcs.values():
                            if (
                                item == arc.source.parent_block()
                                or item == arc.dest.parent_block()
                            ):
                                is_connected = True
                                break
                        # Add to diagram if connected
                        if is_connected:
                            components[item] = type_

        return components

    def _construct_stream_labels(self):
        # Construct the stream labels
        from idaes.core.util.tables import (
            stream_states_dict,
        )  # deferred to avoid circ. import

        # We might have this information from generating self.serialized_components
        # but I (Makayla) don't know how that connects to the stream names so this
        # will be left alone for now
        for stream_name, stream_value in stream_states_dict(self.arcs).items():
            label = ""
            for var, var_value in stream_value.define_display_vars().items():
                var = var.capitalize()

                for k, v in var_value.items():
                    if k is None:
                        label += f"{var} {round(value(v), self._sig_figs)}\n"
                    else:
                        label += f"{var} {k} {round(value(v), self._sig_figs)}\n"
            self.labels[stream_name] = label[:-2]

    def _map_edges(self):
        # Map the arcs to the ports to construct the edges
        used_ports = set()
        for name, arc in self.arcs.items():
            try:  # This is necessary because for internally-nested arcs we may not record ports
                self.edges[name] = {
                    "source": self.ports[arc.source],
                    "dest": self.ports[arc.dest],
                }
                used_ports.add(arc.source)
                used_ports.add(arc.dest)
            except KeyError as error:
                self._logger.error(
                    f"Unable to find port. {name}, {arc.source}, {arc.dest}"
                )

        # Note we're only using the keys from self.ports, {port: parentcomponent}
        return set(self.ports) - used_ports

    def _make_valid_json(self, df):
        """
        Replacing NaN, Infinity, and Negative Infinity with strings to make valid JSON for front end consumption
        Args:
            df: The Pandas dataframe to replace NaN, Infinity, and Negative Infinity with strings
        Return:
            The dataframe that now has valid JSON
        """
        df = df.replace(np.nan, "NaN").replace(np.inf, "Inf").replace(np.NINF, "-Inf")
        return df

    def _add_unit_model_with_ports(self, unit, unit_type):
        unit_name = unit.getname()
        if unit.parent_block() == self.flowsheet:
            # The unit is top-level and therefore should be displayed.
            self.unit_models[unit] = {
                "name": unit_name,
                "type": unit_type,
            }

            self._unit_name_used_count[unit_name] += 1
            for port in unit.component_objects(Port, descend_into=False):
                self.ports[port] = unit

            performance_contents, stream_df = unit.serialize_contents()
            if stream_df is not None and not stream_df.empty:
                # If there is a stream dataframe then we need to reset the index so we can get the variable names
                # and then rename the "index"
                stream_df = stream_df.reset_index().rename(
                    columns={"index": "Variable"}
                )
                stream_df = self._make_valid_json(stream_df)
            self._serialized_contents[unit_name]["stream_contents"] = stream_df

            performance_df = pd.DataFrame()
            if performance_contents:
                # If performance contents is not empty or None then stick it into a dataframe and convert the
                # GeneralVars to actual values
                performance_df = pd.DataFrame(
                    performance_contents["vars"].items(), columns=["Variable", "Value"]
                )
                performance_df["Value"] = performance_df["Value"].map(value)
                performance_df = self._make_valid_json(performance_df)
            self._serialized_contents[unit_name][
                "performance_contents"
            ] = performance_df

        elif unit in self._known_endpoints:
            # Unit is a subcomponent AND it is connected to an Arc. Or maybe it's in
            # an indexed block. Find the top-level parent unit and assign the
            # serialized link to the parent.
            parent_unit = unit.parent_block()
            while not parent_unit == self.flowsheet:
                parent_unit = parent_unit.parent_block()

            self.unit_models[parent_unit] = {
                "name": parent_unit.getname(),
                "type": self.get_unit_model_type(parent_unit),
            }
            for port in unit.component_objects(Port, descend_into=False):
                self.ports[port] = parent_unit
                # Add all of this subcomponent's ports to used_ports; realistically, this is only relevant
                # in a situation where the subcomponent has ports that are not connected to an Arc that we
                # intend to display. Otherwise, the port is marked as "used" when we traverse the Arcs.
                self._used_ports.add(port)
        else:
            # The unit is neither top-level nor connected; do not display this unit, since it is a subcomponent.
            pass

    @staticmethod
    def get_unit_model_type(unit) -> str:
        """Get the 'type' of the unit model."""
        unit_module_type = "unknown"
        if hasattr(unit, "base_class_module"):
            unit_module = unit.base_class_module()
            unit_module_type = unit_module.split(".")[-1]
        else:
            if _log.isEnabledFor(logging.DEBUG):
                _log.debug(f"Could not find type of unit model for: {unit}")
        return unit_module_type

    def _identify_implicit_feeds_and_products(self, untouched_ports):
        """Identify feeds and products not explicitly defined by the user by
        examining the names of ports not connected to arcs.

        This is intended for use on ports of top level unit models.
        It is unclear if it works with nested unit models (probably not).
        """
        for port in sorted(untouched_ports, key=lambda p: str(p)):
            port_name = str(port).split(".")[-1]
            unit_name = self._unique_unit_name(
                port_name
            )  # this becomes the I/O's ID and label
            edge_name = f"s_{unit_name}"
            # identify ports per INLET_REGEX/OUTLET_REGEX
            # then pretend they're unit models
            # then add their edges
            feed_match = self.INLET_REGEX.search(port_name)
            product_match = self.OUTLET_REGEX.search(port_name)
            if feed_match and product_match:
                self._logger.warning(
                    f"Port looks like both a feed and a product: "
                    f"name={port_name} feed-expr={self.INLET_REGEX} "
                    f"product-expr={self.OUTLET_REGEX}"
                )
            elif feed_match or product_match:
                type_ = "feed" if feed_match else "product"
                unit_port = self._PseudoUnit(type_.capitalize(), unit_name)
                self.unit_models[unit_port] = {
                    "name": unit_port.getname(),
                    "type": type_,
                }
                # Add edge. src/dst are reversed for feed vs. product
                src = unit_port if type_ == "feed" else self.ports[port]
                dst = unit_port if type_ != "feed" else self.ports[port]
                self.edges[edge_name] = {"source": src, "dest": dst}
                # Add label
                self.labels[edge_name] = f"{type_} info"
            else:
                self._logger.warning(
                    f"Port is neither a feed nor a product: "
                    f"name={port_name} feed-expr={self.INLET_REGEX} "
                    f"product-expr={self.OUTLET_REGEX}"
                )

    def _unique_unit_name(self, base_name):
        """Prevent name collisions by simply appending a suffix based on how many
        times the name has been used.
        """
        self._unit_name_used_count[base_name] += 1
        return f"{base_name}_{self._unit_name_used_count[base_name]}"

    def _construct_output_json(self):
        self._construct_model_json()
        self._construct_jointjs_json()

    def _construct_model_json(self):
        from idaes.core.util.tables import (
            create_stream_table_dataframe,
        )  # deferred to avoid circular import

        # Get the stream table and add it to the model json
        # Change the index of the pandas dataframe to not be the variables
        self._stream_table_df = (
            create_stream_table_dataframe(self.arcs)
            # Change the index of the pandas dataframe to not be the variables
            .reset_index()
            .rename(columns={"index": "Variable"})
            .reset_index()
            .rename(columns={"index": ""})
            .applymap(
                lambda x: round(x, self._sig_figs) if isinstance(x, (int, float)) else x
            )
        )

        # Change NaNs to None for JSON
        self._stream_table_df = self._stream_table_df.where(
            (pd.notnull(self._stream_table_df)), None
        )

        # Puts df in this format for easier parsing in the javascript table:
        # {'index': ["('Liq', 'benzene')", "('Liq', 'toluene')", "('Liq', 'hydrogen')", "('Liq', 'methane')", "('Vap', 'benzene')", "('Vap', 'toluene')", "('Vap', 'hydrogen')", "('Vap', 'methane')", 'temperature', 'pressure'],
        # 'columns': ['s03', 's04', 's05', 's06', 's08', 's09', 's10'],
        # 'data': [[0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], [0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5], [298.15, 298.15, 298.15, 298.15, 298.15, 298.15, 298.15], [101325.0, 101325.0, 101325.0, 101325.0, 101325.0, 101325.0, 101325.0]]}
        self._out_json["model"]["stream_table"] = self._stream_table_df.to_dict("split")

        self._out_json["model"]["id"] = self.name
        self._out_json["model"]["unit_models"] = {}
        self._out_json["model"]["arcs"] = {}

        for unit_model in self.unit_models.values():
            unit_name = unit_model["name"]
            unit_type = unit_model["type"]
            unit_icon = UnitModelIcon(unit_type)

            unit_contents = {
                "type": unit_type,
                "image": "/images/icons/" + unit_icon.icon,
            }
            if unit_name in self._serialized_contents:
                for pfx in "performance", "stream":
                    content_type = pfx + "_contents"
                    c = (
                        self._serialized_contents[unit_name][content_type]
                        .applymap(
                            lambda x: round(x, self._sig_figs)
                            if isinstance(x, (int, float))
                            else x
                        )
                        .to_dict("index")
                    )
                    # ensure that keys are strings (so it's valid JSON)
                    unit_contents[content_type] = {str(k): v for k, v in c.items()}

            self._out_json["model"]["unit_models"][unit_name] = unit_contents

        for edge, edge_info in self.edges.items():
            self._out_json["model"]["arcs"][edge] = {
                "source": edge_info["source"].getname(),
                "dest": edge_info["dest"].getname(),
                "label": self.labels[edge],
            }

    def _construct_jointjs_json(self):
        self._out_json["cells"] = []

        # Start out in the top left corner until we get a better inital layout
        x_pos = 10
        y_pos = 10
        y_starting_pos = 10

        default_icon = UnitModelIcon()
        for component, unit_attrs in self.unit_models.items():
            unit_icon = UnitModelIcon(unit_attrs["type"])
            try:
                self._create_image_jointjs_json(
                    x_pos,
                    y_pos,
                    unit_attrs["name"],
                    unit_icon.icon,
                    unit_attrs["type"],
                    unit_icon.link_positions,
                )
            except KeyError as e:
                self._logger.info(
                    f'Unable to find icon for {unit_attrs["type"]}. Using default icon'
                )
                self._create_image_jointjs_json(
                    x_pos,
                    y_pos,
                    unit_attrs["name"],
                    default_icon.icon,
                    unit_attrs["type"],
                    default_icon.link_positions,
                )

            # If x_pos it greater than 700 then start another diagonal line
            if x_pos >= 700:
                x_pos = 100
                y_pos = y_starting_pos
                y_starting_pos += 100
            else:
                x_pos += 100
                y_pos += 100

        for name, ports_dict in self.edges.items():
            umst = self.unit_models[ports_dict["source"]]["type"]  # alias
            dest = ports_dict["dest"]
            if hasattr(ports_dict["source"], "vap_outlet"):
                # TODO Figure out how to denote different outlet types. Need to
                # deal with multiple input/output offsets
                for arc in list(self.arcs.values()):
                    if (
                        self.ports[arc.dest] == dest
                        and arc.source == ports_dict["source"].vap_outlet
                    ):
                        source_anchor = "top"
                    else:
                        source_anchor = "bottom"
            else:
                source_anchor = "out"

            # The source_port and dest_port should be replaced by actual names in case there are multiple
            # inlets and outlets between the same two unit models
            source_port = "out"
            dest_port = "in"
            self._create_link_jointjs_json(
                source_port,
                dest_port,
                ports_dict["source"].getname(),
                dest.getname(),
                name,
                self.labels[name],
            )

    def _create_image_jointjs_json(self, x_pos, y_pos, name, image, title, port_groups):
        # Create the jointjs for a given image
        entry = {}
        entry["type"] = "standard.Image"
        entry["position"] = {"x": x_pos, "y": y_pos}
        # The icon width and height default to 50x50 making all icons a square. This will need to be changed
        # when we have more unit models that should not be square. Probaly add it to the icon mapping
        icon_width = 50
        icon_height = 50
        entry["size"] = {"width": icon_width, "height": icon_height}
        # We want the icons to not be at an angle initially
        angle = 0
        entry["angle"] = angle
        entry["id"] = name
        # This defines what layer the icon is on
        z = 1
        entry["z"] = z
        entry["ports"] = port_groups
        entry["attrs"] = {
            "image": {"xlinkHref": "/images/icons/" + image},
            "label": {"text": name},
            "root": {"title": title},
        }
        self._out_json["cells"].append(entry)

    def _create_link_jointjs_json(
        self, source_port, dest_port, source_id, dest_id, name, label
    ):
        # Create the joint js for a link
        # Set the padding to 10. Makayla saw it in a jointjs example
        padding = 10
        # Set the initial offset position for the link labels. Makayla saw these numbers in a jointjs example
        position_distance = 0.66
        position_offset = -40
        z = 2
        entry = {
            "type": "standard.Link",
            "source": {"id": source_id, "port": source_port},
            "target": {"id": dest_id, "port": dest_port},
            "router": {"name": "orthogonal", "padding": padding},
            "connector": {"name": "jumpover", "attrs": {"line": {"stroke": "#5c9adb"}}},
            "id": name,
            "labels": [
                # This label MUST be first or the show/hide will fail
                {
                    "attrs": {
                        # Start with the labels off
                        "rect": {
                            "fill": "#d7dce0",
                            "stroke": "white",
                            "stroke-width": 0,
                            "fill-opacity": 0,
                        },
                        "text": {
                            "text": label,
                            "fill": "black",
                            "text-anchor": "left",
                            "display": "none",
                        },
                    },
                    "position": {
                        "distance": position_distance,
                        "offset": position_offset,
                    },
                },
                {"attrs": {"text": {"text": name}}},
            ],
            "z": z,
        }
        self._out_json["cells"].append(entry)

    class _PseudoUnit:
        """
        Unit-like object in order to emulate missing unit models for implicit feed/product
        eventually may need to actually implement/subclass
        """

        def __init__(self, typename, id):
            self.name = typename
            self.id = id

        def getname(self):
            return self.id


class FlowsheetDiff:
    """Compute a flowsheet model 'diff' and use that to compute an updated layout.

    Flowsheets are serialized by the core library. See :func:`validate_flowsheet` for the required format.

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
                "model": copy.deepcopy(self._new["model"])
                if do_copy
                else self._new["model"],
                "cells": self._layout,
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
                        "rect": {
                            "fill": "#d7dce0",
                            "stroke": "#FFFFFF",
                            "stroke-width": 1,
                        },
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
