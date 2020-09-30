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
import json
import os
import re

from collections import defaultdict
from pandas import DataFrame

import idaes.logger

from idaes.core import UnitModelBlockData
from idaes.core.util.tables import stream_states_dict
from idaes.ui.link_position_mapping import link_position_mapping
from idaes.ui.icon_mapping import icon_mapping
from pyomo.environ import Block, value
from pyomo.network import Arc
from pyomo.network.port import Port


class FileBaseNameExistsError(Exception):
    pass


class FlowsheetSerializer:
    #: Regular expression identifying inlets by last component of ports' name
    INLET_REGEX = re.compile(r'^in_|^feed_|^inlet_|^in$|^feed$|^inlet$|_in$|_feed$|_inlet$') # TODO case insensitivity
    #: Regular expression identifying outlets by last component of ports' name
    OUTLET_REGEX = re.compile(r'^out_|^prod_|^outlet_|^out$|^prod$|^outlet$|_out$|_prod$|_outlet$')

    def __init__(self):
        self.unit_models = {}  # {unit: {"name": unit.getname(), "type": str?}}
        self.arcs = {}  # {Arc.getname(): Arc}
        self.ports = {}  # {Port: parent_unit}
        self.edges = defaultdict(list)  # {name: {"source": unit, "dest": unit}}
        self.orphaned_ports = {}
        self.labels = {}
        self.out_json = {"model": {}}
        self.serialized_contents = defaultdict(dict)
        self.name = ""
        self.flowsheet = None
        self._used_ports = set()
        self._known_endpoints = set()
        self._used_unit_names = set()
        self._logger = idaes.logger.getLogger(__name__)

    def get_edges(self):
        return self.edges

    def get_ports(self):
        return self.ports

    def get_unit_models(self):
        return self.unit_models

    def serialize(self, flowsheet, name):
        """
        Serializes the flowsheet into one dict with two sections.

        The "model" section contains the id of the flowsheet and the
        unit models and arcs. This will be used to compare the model and convert
        to jointjs

        The "cells" section is the jointjs readable code.

        .. code-block:: json

        {
            "model": {
                "id": "id", 
                "unit_models": {
                    "M101": {
                        "image": "mixer.svg", 
                        "type": "mixer",
                        "performance_contents": {
                            "0": {
                                "Variable": "Heat Duty", 
                                "Value": "0.0"
                            }
                        },
                        "stream_contents": {
                            "0": {
                                "Variable": "temperature", 
                                "Inlet": ".01", 
                                "Outlet": "12"
                            }
                        }
                    }
                },
                "arcs": {
                    "s03": {
                        "source": "M101", 
                        "dest": "H101", 
                        "label": "molar flow ("Vap", "hydrogen") 0.5"
                    }
                }
            },
            "cells": [{ "--jointjs code--": "--jointjs code--" }]
        }

        :param flowsheet: The flowsheet to save. Usually fetched from the model.
        :param name: The name of the flowsheet. This will be used as the model id
        :return: None

        Usage example:
            m = ConcreteModel()
            m.fs = FlowsheetBlock(...)
            ...
            serializer = FlowsheetSerializer()
            serializer.save(m.fs, "output_file")
        """
        self.name = name
        self.flowsheet = flowsheet
        self._ingest_flowsheet()
        self._construct_output_json()
        return self.out_json



    def _ingest_flowsheet(self):
        # Stores information on the connectivity and components of the input flowsheet
        self._identify_arcs()
        self._identify_unit_models()
        self._construct_stream_labels()
        untouched_ports = self._map_edges()
        self._identify_implicit_inlets_and_outlets(untouched_ports)

    def _identify_arcs(self):
        # Identify the arcs and known endpoints and store them
        for component in self.flowsheet.component_objects(Arc, descend_into=False):
            self.arcs[component.getname()] = component
            self._known_endpoints.add(component.source.parent_block())
            self._known_endpoints.add(component.dest.parent_block())

    def _identify_unit_models(self):
        # Identify the unit models and ports and store them
        for component in self.flowsheet.component_objects(Block, descend_into=True):
            # TODO try using component_objects(ctype=X)
            if isinstance(component, UnitModelBlockData):
                self._add_unit_model_with_ports(component)
            else:
                # some unit models are nested within Indexed blocks
                component_object_op = getattr(component, "component_object", None)
                if not callable(component_object_op):
                    for item in component.parent_component().values():
                        if isinstance(item, UnitModelBlockData):
                            if any((item == arc.source.parent_block() or item == arc.dest.parent_block()) for arc in self.arcs.values()):
                                self._add_unit_model_with_ports(item)

    def _construct_stream_labels(self):
        # Construct the stream labels

        # We might have this information from generating self.serialized_components but I (Makayla) don't
        # know how that connects to the stream names so this will be left alone for now
        for stream_name, stream_value in stream_states_dict(self.arcs).items():
            label = ""
            for var, var_value in stream_value.define_display_vars().items():
                var = var.capitalize()

                for k, v in var_value.items():
                    if k is None:
                        label += f"{var} {value(v)}\n"
                    else:
                        label += f"{var} {k} {value(v)}\n"
            self.labels[stream_name] = label[:-2]

    def _map_edges(self):
        # Map the arcs to the ports to construct the edges
        used_ports = set()
        for name, arc in self.arcs.items():
            try:  # this is currently(?) necessary because for internally-nested arcs we may not record ports?
                self.edges[name] = {"source": self.ports[arc.source], 
                                    "dest": self.ports[arc.dest]}
                used_ports.add(arc.source)
                used_ports.add(arc.dest)
            except KeyError as error:
                self._logger.error(f"Unable to find port. {name}, {arc.source}, {arc.dest}")

        # note we're only using the keys from self.ports, {port: parentcomponent}
        return set(self.ports) - used_ports

    def _add_unit_model_with_ports(self, unit):
        unit_name = unit.getname()
        if unit.parent_block() == self.flowsheet:
            # The unit is top-level and therefore should be displayed.
            self.unit_models[unit] = {
                    "name": unit_name,
                    "type": self._get_unit_model_type(unit)
                }

            self._used_unit_names.add(unit_name)
            for port in unit.component_objects(Port, descend_into=False):
                self.ports[port] = unit

            performance_contents, stream_df = unit.serialize_contents()

            if not stream_df.empty:
                # If there is a stream dataframe then we need to reset the index so we can get the variable names
                # and then rename the "index" 
                stream_df = stream_df.reset_index().rename(columns={"index": "Variable"})
            self.serialized_contents[unit_name]["stream_contents"] = stream_df

            performance_df = DataFrame()
            if performance_contents:
                # If performance contents is not empty or None then stick it into a dataframe and convert the
                # GeneralVars to actual values
                performance_df = DataFrame(performance_contents["vars"].items(), columns=["Variable", "Value"])
                performance_df['Value'] = performance_df['Value'].map(lambda v: value(v))
            self.serialized_contents[unit_name]["performance_contents"] = performance_df            

        elif unit in self._known_endpoints:
            # Unit is a subcomponent AND it is connected to an Arc. Or maybe it's in an indexed block TODO CHECK
            # Find the top-level parent unit and assign the serialized link to the parent.
            parent_unit = unit.parent_block()
            while not parent_unit == self.flowsheet:
                parent_unit = parent_unit.parent_block()

            self.unit_models[parent_unit] = {
                "name": parent_unit.getname(),
                "type": self._get_unit_model_type(parent_unit)
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

    def _get_unit_model_type(self, unit):
        # Get the unit models type
        return unit._orig_module.split(".")[-1]  # TODO look for/create equivalent getter, as getname() above

    def _identify_implicit_inlets_and_outlets(self, untouched_ports):
        # Identify feeds and products not explicitly defined by the user by examining the names of ports 
        # not connected to arcs 
        for port in sorted(untouched_ports, key=lambda port: str(port)):  # TODO: sorting required for determinism; better way?
            portname = str(port).split('.')[-1]
            unit_name = self._unique_unit_name(portname)  # this becomes the I/O's ID and label
            edgename = f's_{unit_name}'
            # identify ports per INLET_REGEX/OUTLET_REGEX
            # then pretend they're unit models
            # then add their edges
            inlet_match = self.INLET_REGEX.search(portname)
            if inlet_match:
                # name the feed "unit model" and its connecting edge with the name of the port itself
                feedport = self._PseudoUnit('Feed', unit_name)
                self._used_unit_names.add(unit_name)
                self.unit_models[feedport] = {
                    "name": feedport.getname(),
                    "type": "feed"
                }
                self.edges[edgename] = {
                    "source": feedport,
                    "dest": self.ports[port]  # TODO ensure this works on nested unit models
                }
                self.labels[edgename] = "inlet info"  # TODO - fetch the actual information
                continue

            outlet_match = self.OUTLET_REGEX.search(portname)
            if outlet_match:
                prodport = self._PseudoUnit('Product', unit_name)
                self._used_unit_names.add(unit_name)
                self.unit_models[prodport] = {
                    "name": prodport.getname(),
                    "type": "product"
                }
                self.edges[edgename] = {
                    "source": self.ports[port],  # TODO see above
                    "dest": prodport
                }
                self.labels[edgename] = "outlet info"  # TODO - fetch the actual information
                continue

            # TODO: deal with remaining loose ports here?

    def _unique_unit_name(self, base_name):
        # Prevent name collisions by simply appending a number
        name = base_name
        increment = 0
        while name in self._used_unit_names:
            increment += 1
            name = f'{base_name}_{increment}'
        return name

    def _construct_output_json(self):
        self._construct_model_json()
        self._construct_jointjs_json()

    def _construct_model_json(self):
        self.out_json["model"]["id"] = self.name
        self.out_json["model"]["unit_models"] = {}
        self.out_json["model"]["arcs"] = {}

        for unit_model in self.unit_models.values():
            unit_name = unit_model["name"]
            unit_type = unit_model["type"]
            
            performance_contents = {}
            stream_contents = {}
            if unit_name in self.serialized_contents:
                performance_contents = self.serialized_contents[unit_name]["performance_contents"].to_dict('index')
                stream_contents = self.serialized_contents[unit_name]["stream_contents"].to_dict('index')

            self.out_json["model"]["unit_models"][unit_name] = {
                "type": unit_type,
                "image": "/images/icons/" + icon_mapping(unit_type),
                "performance_contents": performance_contents,
                "stream_contents": stream_contents
            }

        for edge, edge_info in self.edges.items():
            self.out_json["model"]["arcs"][edge] = \
                {"source": edge_info["source"].getname(),
                 "dest": edge_info["dest"].getname(),
                 "label": self.labels[edge]}

    def _construct_jointjs_json(self):
        self.out_json["cells"] = []

        # Start out in the top left corner until we get a better inital layout
        x_pos = 10
        y_pos = 10
        y_starting_pos = 10

        for component, unit_attrs in self.unit_models.items():
            try:
                self._create_image_jointjs_json(
                    self.out_json,
                    x_pos,
                    y_pos,
                    unit_attrs["name"],
                    icon_mapping(unit_attrs["type"]),
                    unit_attrs["type"],
                    link_position_mapping[unit_attrs["type"]]
                )
            except KeyError as e:
                self._logger.info(f'Unable to find icon for {unit_attrs["type"]}. Using default icon')
                self._create_image_jointjs_json(self.out_json, 
                                               x_pos, 
                                               y_pos, 
                                               unit_attrs["name"], 
                                               icon_mapping("default"), 
                                               unit_attrs["type"],
                                               link_position_mapping["default"])

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
                    if (self.ports[arc.dest] == dest and arc.source == ports_dict["source"].vap_outlet):
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
                self.out_json, 
                source_port, 
                dest_port, 
                ports_dict["source"].getname(), 
                dest.getname(), 
                name,
                self.labels[name]
            )

    def _create_image_jointjs_json(self, out_json, x_pos, y_pos, name, image, title, port_groups):
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
        z = (1,)
        entry["z"] = z
        entry["ports"] = port_groups
        entry["attrs"] = {
            "image": {"xlinkHref": "/images/icons/" + image},
            "label": {
                "text": name
            },
            "root": {"title": title},
        }
        out_json["cells"].append(entry)

    def _create_link_jointjs_json(self, out_json, source_port, dest_port, 
                                 source_id, dest_id, name, label):  
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
            "connector": {"name": "normal", 
                          "attrs": {"line": {"stroke": "#5c9adb"}}},
            "id": name,
            "labels": [
                # This label MUST be first or the show/hide will fail
                {"attrs": {
                    # Start with the labels off
                    "rect": {"fill": "#d7dce0", "stroke": "white", "stroke-width": 0, "fill-opacity": 0},
                    "text": {
                        "text": label,
                        "fill": "black",
                        "text-anchor": "left",
                        "display": "none"
                    },
                },
                "position": {
                    "distance": position_distance,
                    "offset": position_offset
                }},
                {"attrs": {
                    "text": {
                        "text": name
                    }
                }}
            ],
            "z": z
        }
        out_json["cells"].append(entry)

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
