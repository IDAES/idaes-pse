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
from collections import defaultdict, deque
import copy
import json
import logging
import re
from typing import Dict, List, Tuple

# third-party
import pandas as pd
import numpy as np
import pint
from pyomo.environ import Block, value
from pyomo.network import Arc
from pyomo.network.port import Port

# package
from idaes import logger
from idaes.core.ui.icons.icons import UnitModelIcon
from idaes.core.ui.icons.positioning import UnitModelsPositioning

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

    # Key that represents feeds
    FEED = "feed"

    # Key that represents products
    PRODUCT = "product"

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
        self.streams = {}  # {Arc.getname(): Arc} or {Port.getname(): Port}
        self.ports = {}  # {Port: parent_unit}
        self.edges = defaultdict(dict)  # {name: {"source": unit, "dest": unit}}
        self.adj_list = defaultdict(set)  # {name: (neighbor1, neighbor2, ...)}
        self.orphaned_ports = {}
        self.labels = {}
        self._stream_table_df = None
        self._ordered_stream_names = deque()
        self._out_json = {"model": {}, "routing_config": {}}
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
        self._positioning_model = None
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

        untouched_ports = self._map_edges()
        self._identify_implicit_feeds_and_products(untouched_ports)

        # Construct Labels (shown in popups boxes in the UI)
        self._construct_stream_labels()

    def _identify_arcs(self):
        # Identify the arcs and known endpoints and store them
        for component in self.flowsheet.component_objects(Arc, descend_into=False):
            self.streams[component.getname()] = component
            self._known_endpoints.add(component.source.parent_block())
            self._known_endpoints.add(component.dest.parent_block())
            self._ordered_stream_names.append(component.getname())

    def _identify_unit_models(self) -> Dict:
        from idaes.core import UnitModelBlockData  # avoid circular import
        from idaes.core.base.property_base import PhysicalParameterBlock, StateBlock

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
                        for stream in self.streams.values():
                            if (
                                item == stream.source.parent_block()
                                or item == stream.dest.parent_block()
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
        for stream_name, stream_value in stream_states_dict(self.streams).items():
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
        for name, stream in self.streams.items():
            try:  # This is necessary because for internally-nested arcs we may not record ports
                src = self.ports[stream.source]
                dst = self.ports[stream.dest]
                self.edges[name] = {
                    "source": src,
                    "dest": dst,
                }
                self.adj_list[src.getname()].add(dst.getname())
                used_ports.add(stream.source)
                used_ports.add(stream.dest)
            except KeyError as error:
                self._logger.error(
                    f"Unable to find port. {name}, {stream.source}, {stream.dest}"
                )

        # Note we're only using the keys from self.ports, {port: parentcomponent}
        untouched_ports = set(self.ports) - used_ports
        return untouched_ports

    def _clean_values(self, df):
        """
        Replacing NaN, Infinity, and Negative Infinity with strings to make valid JSON for front end consumption
        Args:
            df: The Pandas dataframe to replace NaN, Infinity, and Negative Infinity with strings
        Return:
            The dataframe that now has valid JSON
        """
        return df.replace(np.nan, "NaN").replace(np.inf, "Inf").replace(np.NINF, "-Inf")

    def _clean_units(self, df):
        """
        Convert Unit pint object to different string formats in Pandas dataframe
        Args:
            df: The Pandas dataframe to replace pint objects to different string formats
        Return:
            The dataframe that now has valid JSON
        """

        def get_valid_unit_format(pint_unit):
            """Get different formats from the Units' pint object"""
            if isinstance(pint_unit, pint.Unit):
                return {
                    "raw": str(pint_unit),
                    "html": "{:~H}".format(pint_unit),
                    "latex": "{:~L}".format(pint_unit),
                }
            else:
                return {
                    "raw": str(pint_unit),
                    "html": "",
                    "latex": "",
                }

        if "Units" in df.columns:
            df["Units"] = df["Units"].apply(
                lambda pint_unit: get_valid_unit_format(pint_unit)
            )
        return df

    def _make_valid_json(self, df):
        """
        Replacing non-serializable objects to have a final valid JSON object
        Args:
            df: The Pandas dataframe to replace non-serializable objects within
        Return:
            The dataframe that now has valid JSON
        """
        df = self._clean_values(df)
        df = self._clean_units(df)
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
        from idaes.core.util.tables import stream_states_dict

        for port in sorted(untouched_ports, key=lambda p: str(p)):
            port_name = port.getname()
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
                type_ = self.FEED if feed_match else self.PRODUCT
                unit_port = self._PseudoUnit(type_.capitalize(), unit_name)
                self.unit_models[unit_port] = {
                    "name": unit_port.getname(),
                    "type": type_,
                }
                # Add edge. src/dst are reversed for feed vs. product
                src = unit_port if type_ == self.FEED else self.ports[port]
                dst = unit_port if type_ != self.FEED else self.ports[port]
                self.edges[edge_name] = {"source": src, "dest": dst}
                self.adj_list[src.getname()].add(dst.getname())
                # Add label
                self.labels[edge_name] = f"{type_} info"
                # Check if we can add the port to the stream table
                try:
                    # Check if we can extract a state block from this port.
                    # Otherwise this function will raise an error and this
                    # port won't be added to the stream table
                    _ = stream_states_dict({"port": port})

                    self.streams[edge_name] = port

                    # Ordering the feeds vs products in the code to display
                    # the right column order in the stream table
                    if type_ == self.FEED:
                        self._ordered_stream_names.appendleft(edge_name)
                    else:
                        self._ordered_stream_names.append(edge_name)
                except Exception:
                    self._logger.warning(
                        f"Cannot extract state block from Port: "
                        f"name={port_name}. "
                        f"Please add Feed & Product blocks with Arcs to show "
                        f"inlet and outlet stream values in the Stream Table"
                    )
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
        self._positioning_model = UnitModelsPositioning(self.adj_list, self.unit_models)
        self._construct_model_json()
        self._construct_jointjs_json()

    def _construct_model_json(self):
        from idaes.core.util.tables import (
            create_stream_table_ui,
        )  # deferred to avoid circular import

        # Get the stream table and add it to the model json
        # Change the index of the pandas dataframe to not be the variables
        self._stream_table_df = (
            create_stream_table_ui(self.streams)
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

        self._stream_table_df = self._make_valid_json(self._stream_table_df)

        # Order the stream table based on the right order:
        # feed streams -> middle streams -> product streams
        self._ordered_stream_names.appendleft("Units")
        self._ordered_stream_names.appendleft("Variable")
        self._stream_table_df = self._stream_table_df[self._ordered_stream_names]

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

    def _add_port_item(self, cell_index, group, id):
        """Add port item to jointjs element"""
        new_port_item = {"group": group, "id": id}
        if new_port_item not in self._out_json["cells"][cell_index]["ports"]["items"]:
            self._out_json["cells"][cell_index]["ports"]["items"].append(new_port_item)

    def _construct_jointjs_json(self):
        def create_jointjs_image(
            unit_icon: UnitModelIcon, unit_name, unit_type, x_pos, y_pos
        ):
            """Create jointjs element 'standard.Image' type in json format"""
            try:
                return self._create_image_jointjs_json(
                    x_pos,
                    y_pos,
                    unit_name,
                    unit_icon.icon,
                    unit_type,
                    unit_icon.link_positions,
                )
            except KeyError as e:
                self._logger.info(
                    f"Unable to find icon for {unit_type}. Using default icon"
                )
                default_icon = UnitModelIcon()
                return self._create_image_jointjs_json(
                    x_pos,
                    y_pos,
                    unit_name,
                    default_icon.icon,
                    unit_type,
                    default_icon.link_positions,
                )

        self._out_json["cells"] = []

        # Start out in the top left corner until we get a better inital layout
        x_pos = 10
        y_pos = 10

        track_jointjs_elements = {}

        port_index_increment = 0

        # Go through all the edges/links and create the necessary Unit models
        # that are connected to these edges.
        for link_name, ports_dict in self.edges.items():
            src = ports_dict["source"]
            dest = ports_dict["dest"]

            src_unit_name = self.unit_models[src]["name"]
            src_unit_type = self.unit_models[src]["type"]
            src_unit_icon = UnitModelIcon(src_unit_type)

            dest_unit_name = self.unit_models[dest]["name"]
            dest_unit_type = self.unit_models[dest]["type"]
            dest_unit_icon = UnitModelIcon(dest_unit_type)

            if src_unit_name not in track_jointjs_elements:
                x_pos, y_pos = self._positioning_model.get_position(src_unit_name)
                cell_index = create_jointjs_image(
                    src_unit_icon, src_unit_name, src_unit_type, x_pos, y_pos
                )
                track_jointjs_elements[src_unit_name] = cell_index

            if dest_unit_name not in track_jointjs_elements:
                x_pos, y_pos = self._positioning_model.get_position(dest_unit_name)
                cell_index = create_jointjs_image(
                    dest_unit_icon, dest_unit_name, dest_unit_type, x_pos, y_pos
                )
                track_jointjs_elements[dest_unit_name] = cell_index

            # TODO Figure out how to denote different connection direction types.
            # Need to deal with multiple input/output offsets.
            # e.g top, bottom, etc

            # The source_port and dest_port should be replaced by actual names in case there are multiple
            # inlets and outlets between the same two unit models
            src_port = "out"
            dest_port = "in"

            # We create port ids for both ends: Source and Destination elements
            # and pass it to the link jointjs element for accurate port representation
            src_port_id = port_index_increment
            port_index_increment += 1
            dest_port_id = port_index_increment
            port_index_increment += 1

            # Add source port
            self._add_port_item(
                track_jointjs_elements[src_unit_name], src_port, src_port_id
            )
            # Add destination port
            self._add_port_item(
                track_jointjs_elements[dest_unit_name], dest_port, dest_port_id
            )

            link_index = self._create_link_jointjs_json(
                src_port_id,
                dest_port_id,
                src.getname(),
                dest.getname(),
                link_name,
                self.labels[link_name],
            )

            # Add routing config if edge/link has source or destination elements
            # that has routing specifications. e.g. If destination element requires
            # the link to connect horizontally from the left side.
            if link_name not in self._out_json["routing_config"]:
                self._out_json["routing_config"][link_name] = {
                    "cell_index": link_index,
                    "cell_config": {
                        "gap": {
                            "source": {"x": 0, "y": 0},
                            "destination": {"x": 0, "y": 0},
                        }
                    },
                }
            cell_config_gap = self._out_json["routing_config"][link_name][
                "cell_config"
            ]["gap"]
            if (
                src_unit_icon.routing_config
                and src_port in src_unit_icon.routing_config
            ):
                # The port group has to be specified in the routing config
                cell_config_gap["source"] = src_unit_icon.routing_config[src_port][
                    "gap"
                ]

            if (
                dest_unit_icon.routing_config
                and dest_port in dest_unit_icon.routing_config
            ):
                # The port group has to be specified in the routing config
                cell_config_gap["destination"] = dest_unit_icon.routing_config[
                    dest_port
                ]["gap"]

    def _create_image_jointjs_json(self, x_pos, y_pos, name, image, title, port_groups):
        # Create the jointjs for a given image

        # The icon width and height default to 50x50 making all icons a square. This will need to be changed
        # when we have more unit models that should not be square. Probaly add it to the icon mapping
        icon_width = 50
        icon_height = 50
        # We want the icons to not be at an angle initially
        angle = 0
        # This defines what layer the icon is on
        z = 1

        entry = {
            "type": "standard.Image",
            "position": {"x": x_pos, "y": y_pos},
            "size": {"width": icon_width, "height": icon_height},
            "angle": angle,
            "id": name,
            "z": z,
            "ports": port_groups,
            "attrs": {
                "image": {"xlinkHref": "/images/icons/" + image},
                "label": {"text": name},
                "root": {"title": title},
            },
        }
        self._out_json["cells"].append(entry)
        return (
            len(self._out_json["cells"]) - 1
        )  # return the index of the newly added cell

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
            "router": {"name": "manhattan", "padding": padding},
            "connector": {"name": "jumpover", "attrs": {"line": {"stroke": "#5c9adb"}}},
            "id": name,
            "attrs": {"line": {"stroke": "#979797", "stroke-width": 2}},
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
        return (
            len(self._out_json["cells"]) - 1
        )  # return the index of the newly added cell

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
        self._cell_indices = {}
        self._compute_cell_indices()
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
            model = self._new["model"]
            routing_config = self._new["routing_config"]
            if do_copy:
                model = copy.deepcopy(self._new["model"])
                routing_config = copy.deepcopy(self._new["routing_config"])
            result = {
                "model": model,
                "cells": self._layout,
                "routing_config": routing_config,
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

    def _compute_cell_indices(self):
        """Generate a mapping between the cell id and its index in the 'cells' Array.

        Returns:
            dict with structure: cell_id -> index. e.g. {'F101': 2} where 2 is
            the cell index in the 'cells' array.
        """
        for i, cell in enumerate(self._new["cells"]):
            self._cell_indices[cell["id"]] = i

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
        layout = self._new["cells"]

        # Go through each cell in the old model and copy the user changes to
        # the new one
        for cell in self._old["cells"]:
            cell_id = cell["id"]
            if cell_id in self._cell_indices:
                # check if cell is an element (not link) and 'position',
                # 'angle' & 'attrs.label' keys exist in cell
                if "type" in cell and cell["type"] == "standard.Image":
                    if "position" in cell:
                        layout[self._cell_indices[cell_id]]["position"] = cell[
                            "position"
                        ]
                    if "angle" in cell:
                        layout[self._cell_indices[cell_id]]["angle"] = cell["angle"]
                    if "attrs" in cell and "label" in cell["attrs"]:
                        if "attrs" not in layout[self._cell_indices[cell_id]]:
                            layout[self._cell_indices[cell_id]]["attrs"] = {}
                        layout[self._cell_indices[cell_id]]["attrs"]["label"] = cell[
                            "attrs"
                        ]["label"]
                # check if link and if it has 'vertices' & 'labels[1]'
                # keys exist in cell
                elif "type" in cell and cell["type"] == "standard.Link":
                    if "vertices" in cell:
                        layout[self._cell_indices[cell_id]]["vertices"] = cell[
                            "vertices"
                        ]
                    # labels[0] is reserved for the variables info when labels
                    # are enabled in the UI. We are just copying labels[1] as
                    # it has the stream/link name. e.g. s10, s_liq_outlet, etc
                    if "labels" in cell and len(cell["labels"]) >= 2:
                        layout[self._cell_indices[cell_id]]["labels"][1] = cell[
                            "labels"
                        ][1]

        # Update cells' labels or images
        for item in self._old["cells"]:
            id_ = item["id"]
            cls = "arcs" if "source" in item else "unit_models"
            if id_ in self._diff["change"][cls] and id_ in self._cell_indices:
                values = self._diff["change"][cls][id_]
                cell = layout[self._cell_indices[id_]]
                if cls == "arcs":
                    self._update_arc(cell, values)
                else:
                    self._update_unit_model(cell, values)
        return layout

    @staticmethod
    def _update_arc(item, values):
        item["labels"][0]["attrs"]["text"]["text"] = values["label"]

    @staticmethod
    def _update_unit_model(item, values):
        item["attrs"]["image"]["xlinkHref"] = values["image"]
        item["attrs"]["root"]["title"] = values["type"]
