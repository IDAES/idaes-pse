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
from collections import defaultdict
import json
import os
import re

from idaes.core import UnitModelBlockData
from idaes.core.util.tables import stream_states_dict
from idaes.ui.link_position_mapping import link_position_mapping
from idaes.ui.icon_mapping import icon_mapping
import idaes.logger

from pyomo.environ import Block, value
from pyomo.network.port import Port
from pyomo.network import Arc
from pyomo.core.base.var import Var
from pyomo.core.base.expression import Expression


class FileBaseNameExistsError(Exception):
    pass


class FlowsheetSerializer:
    #: Regular expression identifying inlets by last component of ports' name
    INLET_REGEX = re.compile(r'^in_|^feed_|^inlet_|^in$|^feed$|^inlet$|_in$|_feed$|_inlet$') # TODO case insensitivity
    r''
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
        self.name = ""
        self._used_unit_names = set()
        self._logger = idaes.logger.getLogger(__name__)

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
                        "type": "mixer"
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
        self.serialize_flowsheet(flowsheet)
        self._construct_output_json()
        return self.out_json

    def serialize_flowsheet(self, flowsheet):
        """Stores information on the connectivity and components of the input flowsheet. """
        for component in flowsheet.component_objects(Arc, descend_into=False):
            self.arcs[component.getname()] = component
            
        for component in flowsheet.component_objects(Block, descend_into=True):
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

        used_ports = set()
        self.edges = {}  # TODO: necessary? do the attributes need clearing?
        for name, arc in self.arcs.items():
            try:  # this is currently(?) necessary because for internally-nested arcs we may not record ports?
                self.edges[name] = {"source": self.ports[arc.source], 
                                    "dest": self.ports[arc.dest]}
                self._logger.info(f'source: {arc.source}; val: {self.ports[arc.source]}')
                used_ports.add(arc.source)
                used_ports.add(arc.dest)
            except KeyError as error:
                print(f"Unable to find port. {name}, {arc.source}, {arc.dest}")  # TODO convert to logging

        # note we're only using the keys from self.ports, {port: parentcomponent}
        untouched_ports = set(self.ports) - used_ports
        self._identify_implicit_inlets_and_outlets(untouched_ports)

    def _identify_implicit_inlets_and_outlets(self, untouched_ports):
        for port in sorted(untouched_ports, key=lambda port: str(port)):  # TODO: sorting required for determinism; better way?
            portname = str(port).split('.')[-1]
            unitname = self._unique_unit_name(portname)  # this becomes the I/O's ID and label
            edgename = f's_{unitname}'
            # identify ports per INLET_REGEX/OUTLET_REGEX
            # then pretend they're unit models
            # then add their edges
            inlet_match = self.INLET_REGEX.search(portname)
            if inlet_match:
                # print(f"  ^ inlet found; parent: {str(self.ports[port])}")
                # name the feed "unit model" and its connecting edge with the name of the port itself
                feedport = self._PseudoUnit('Feed', unitname)
                self._used_unit_names.add(unitname)
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
                # print(f"  ^ outlet found; parent: {str(self.ports[port])}")
                prodport = self._PseudoUnit('Product', unitname)
                self._used_unit_names.add(unitname)
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
        '''Prevent name collisions by simply appending a number'''
        name = base_name
        increment = 0
        while name in self._used_unit_names:
            increment += 1
            name = f'{base_name}_{increment}'
        return name

    def create_image_jointjs_json(self, out_json, x_pos, y_pos, name, image, title, port_groups):
        entry = {}
        entry["type"] = "standard.Image"
        # for now, just tile the positions diagonally
        # TODO Make the default positioning better
        entry["position"] = {"x": x_pos, "y": y_pos}
        # TODO Set the width and height depending on the icon rather than default
        entry["size"] = {"width": 50, "height": 50}
        entry["angle"] = 0
        entry["id"] = name
        entry["z"] = (1,)
        entry["ports"] = port_groups
        entry["attrs"] = {
            "image": {"xlinkHref": "/images/icons/" + image},
            "label": {
                "text": name
            },
            "root": {"title": title},
        }
        out_json["cells"].append(entry)

    def create_link_jointjs_json(self, out_json, source_port, dest_port, 
                                 source_id, dest_id, name, label):      
        entry = {
            "type": "standard.Link",
            "source": {"id": source_id, "port": source_port},
            "target": {"id": dest_id, "port": dest_port},
            "router": {"name": "orthogonal", "padding": 10},
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
                    "distance": 0.66,
                    "offset": -40
                }},
                {"attrs": {
                    "text": {
                        "text": name
                    }
                }}
            ],
            "z": 2
        }
        out_json["cells"].append(entry)

    def get_unit_models(self):
        return self.unit_models

    def get_ports(self):
        return self.ports

    def get_edges(self):
        return self.edges

    def _construct_output_json(self):
        self._construct_model_json()
        self._construct_jointjs_json()

    def _construct_model_json(self):
        self.out_json["model"]["id"] = self.name
        self.out_json["model"]["unit_models"] = {}
        self.out_json["model"]["arcs"] = {}

        for unit_model in self.unit_models.values():
            self.out_json["model"]["unit_models"][unit_model["name"]] = {
                "type": unit_model["type"],
                "image": "/images/icons/" + icon_mapping(unit_model["type"])
            }

        for edge, edge_info in self.edges.items():
            self.out_json["model"]["arcs"][edge] = \
                {"source": edge_info["source"].getname(),
                 "dest": edge_info["dest"].getname(),
                 "label": self.labels[edge]}

    def _construct_jointjs_json(self):
        self.out_json["cells"] = []
        x_pos = 10
        y_pos = 10
        y_starting_pos = 10

        for component, unit_attrs in self.unit_models.items():
            try:
                self.create_image_jointjs_json(
                    self.out_json,
                    x_pos,
                    y_pos,
                    unit_attrs["name"],
                    icon_mapping(unit_attrs["type"]),
                    unit_attrs["type"],
                    link_position_mapping[unit_attrs["type"]]
                )
            except KeyError as e:
                print(f'Unable to find icon for {unit_attrs["type"]}. Using default icon')
                self.create_image_jointjs_json(self.out_json, 
                                               x_pos, 
                                               y_pos, 
                                               unit_attrs["name"], 
                                               icon_mapping("default"), 
                                               unit_attrs["type"],
                                               link_position_mapping["default"])
            if x_pos >= 700:
                x_pos = 100
                y_pos = y_starting_pos
                y_starting_pos += 100
            else:
                x_pos += 100
                y_pos += 100

        id_counter = 0
        for name, ports_dict in self.edges.items():
            umst = self.unit_models[ports_dict["source"]]["type"]  # alias
            dest = ports_dict["dest"]
            if hasattr(ports_dict["source"], "vap_outlet"):
                # TODO Figure out how to denote different outlet types. Need to
                #  deal with multiple input/output offsets
                for arc in list(self.arcs.values()):
                    if (self.ports[arc.dest] == dest and arc.source == ports_dict["source"].vap_outlet):
                        source_anchor = "top"
                    else:
                        source_anchor = "bottom"
            else:
                source_anchor = "out"

            self.create_link_jointjs_json(
                self.out_json, 
                "out", 
                "in", 
                ports_dict["source"].getname(), 
                dest.getname(), 
                name,
                self.labels[name]
            )

    def _add_unit_model_with_ports(self, unit):
        unitname = unit.getname()
        self.unit_models[unit] = {
                "name": unitname,
                "type": unit._orig_module.split(".")[-1]  # TODO look for/create equivalent getter, as getname() above
            }
        self._used_unit_names.add(unitname)
        for subcomponent in unit.component_objects(Port, descend_into=True):
            self.ports[subcomponent] = unit

    # unit-like object in order to emulate missing unit models for implicit inlets/outlets
    # eventually may need to actually implement/subclass
    class _PseudoUnit:
        def __init__(self, typename, id):
            self.name = typename
            self.id = id

        def getname(self):
            return self.id
