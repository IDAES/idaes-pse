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

from idaes.core import UnitModelBlockData
from idaes.core.util.tables import stream_states_dict
from idaes.ui.link_position_mapping import link_position_mapping
from idaes.ui.icon_mapping import icon_mapping

from pyomo.environ import Block, value
from pyomo.network.port import Port
from pyomo.network import Arc
from pyomo.core.base.var import Var
from pyomo.core.base.expression import Expression


class FileBaseNameExistsError(Exception):
    pass


class FlowsheetSerializer:
    def __init__(self):
        self.unit_models = {}
        self.arcs = {}
        self.ports = {}
        self.edges = defaultdict(list)
        self.orphaned_ports = {}
        self.labels = {}
        self.out_json = {"model": {}}
        self.name = ""

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
        for component in flowsheet.component_objects(Arc, descend_into=False):
            self.arcs[component.getname()] = component
            
        for component in flowsheet.component_objects(Block, descend_into=True):
            # TODO try using component_objects(ctype=X)
            if isinstance(component, UnitModelBlockData):
                self.unit_models[component] = {
                    "name": component.getname(), 
                    "type": component._orig_module.split(".")[-1]
                }
                for subcomponent in component.component_objects(Port, descend_into=True):
                    self.ports[subcomponent] = component
            else:
                component_object_op = getattr(component, "component_object", None)
                if not callable(component_object_op):
                    for item in component.parent_component().values():
                        if isinstance(item, UnitModelBlockData):
                            if any((item == arc.source.parent_block() or item == arc.dest.parent_block()) for arc in self.arcs.values()):
                                self.unit_models[item] = {
                                    "name": item.getname(), 
                                    "type": item._orig_module.split(".")[-1]
                                }
                                for subcomponent in item.component_objects(Port, descend_into=True):
                                    self.ports[subcomponent] = item
  
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

        self.edges = {}
        for name, arc in self.arcs.items():
            try:
                self.edges[name] = {"source": self.ports[arc.source], 
                                    "dest": self.ports[arc.dest]}
            except KeyError as error:
                print(f"Unable to find port. {name}, {arc.source}, {arc.dest}")

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

        for edge in self.edges:
            self.out_json["model"]["arcs"][edge] = \
                {"source": self.edges[edge]["source"].getname(),
                 "dest": self.edges[edge]["dest"].getname(),
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
            except KeyError:
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
