##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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
from idaes.dmf.ui.link_position_mapping import link_position_mapping
from idaes.dmf.ui.icon_mapping import icon_mapping

from pyomo.environ import Block
from pyomo.network.port import SimplePort
from pyomo.network import Arc


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

    def serialize(self, flowsheet, file_base_name, overwrite=False):
        """
        Serializes the flowsheet and saves it to a file that can be read by the
        idaes-model-vis  jupyter lab extension.

        :param flowsheet: The flowsheet to save. Usually fetched from the model.
        :param file_base_name: The file prefix to the .idaes.vis file produced.
        The file is created/saved
        in the directory that you ran from Jupyter Lab.
        :param overwrite: Boolean to overwrite an existing file_base_name.idaes.vis.
        If True, the existing file with the same file_base_name will be overwritten.
        This will cause you to lose
        any saved layout. 
        If False and there is an existing file with that file_base_name, you will get
        an error
        message stating that you cannot save a file to the file_base_name
        (and therefore overwriting the saved
        layout). If there is not an existing file with that file_base_name then it
        saves as normal.
        Defaults to False.
        :return: None

        Usage example:
            m = ConcreteModel()
            m.fs = FlowsheetBlock(...)
            ...
            serializer = FlowsheetSerializer()
            serializer.save(m.fs, "output_file")
        """
        vis_file_name = file_base_name + ".idaes.vis"
        if os.path.isfile(vis_file_name) and overwrite is False:
            msg = (
                f"{vis_file_name} already exists. If you wish to overwrite "
                f"this file call save() with overwrite=True. "
                "WARNING: If you overwrite the file, you will lose "
                "your saved layout."
            )
            raise FileBaseNameExistsError(msg)
        else:
            print(f"Creating {vis_file_name}")

        self.serialize_flowsheet(flowsheet)
        self._construct_output_json()

        with open(vis_file_name, "w") as out_file:
            json.dump(self.out_json, out_file)

    def serialize_flowsheet(self, flowsheet):
        for component in flowsheet.component_objects(Block, descend_into=False):
            # TODO try using component_objects(ctype=X)
            if isinstance(component, UnitModelBlockData):
                self.unit_models[component] = {
                    "name": component.getname(), 
                    "type": component._orig_module.split(".")[-1]
                }

                for subcomponent in component.component_objects(descend_into=True):
                    if isinstance(subcomponent, SimplePort):
                        self.ports[subcomponent] = component
  
        for component in flowsheet.component_objects(Arc, descend_into=False):
            self.arcs[component.getname()] = component

        for stream_name, value in stream_states_dict(self.arcs).items():
            label = ""

            for var, var_value in value.define_display_vars().items():
                for stream_type, stream_value in var_value.get_values().items():
                    if stream_type:
                        if var == "flow_mol_phase_comp":
                            var = "Molar Flow"
                        label += f"{var} {stream_type} {stream_value}\n"
                    else:
                        var = var.capitalize()
                        label += f"{var} {stream_value}\n"

            self.labels[stream_name] = label[:-2]

        self.edges = {}
        for name, arc in self.arcs.items():
            self.edges[name] = {"source": self.ports[arc.source], 
                                "dest": self.ports[arc.dest]}

    def create_image_jointjs_json(self, out_json, x_pos, y_pos, name, image, title):
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
        entry["attrs"] = {
            "image": {"xlinkHref": image},
            "label": {"text": name},
            "root": {"title": title},
        }
        out_json["cells"].append(entry)

    def create_link_jointjs_json(self, out_json, source_anchor, dest_anchor, 
                                 source_id, dest_id, name, label):
        entry = {
            "type": "standard.Link",
            "source": {"anchor": source_anchor, "id": source_id},
            "target": {"anchor": dest_anchor, "id": dest_id},
            "router": {"name": "orthogonal", "padding": 10},
            "connector": {"name": "normal", 
                          "attrs": {"line": {"stroke": "#5c9adb"}}},
            "id": name,
            "labels": [{
                "attrs": {
                    "rect": {"fill": "#d7dce0", "stroke": "#FFFFFF", 'stroke-width': 1},
                    "text": {
                        "text": label,
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
        self.out_json["model"]["id"] = 0
        self.out_json["model"]["unit_models"] = {}
        self.out_json["model"]["arcs"] = {}

        for unit_model in self.unit_models.values():
            self.out_json["model"]["unit_models"][unit_model["name"]] = {}
            self.out_json["model"]["unit_models"][unit_model["name"]] = {
                "type": unit_model["type"],
                "image": icon_mapping[unit_model["type"]]
            }

        for edge in self.edges:
            self.out_json["model"]["arcs"][edge] = \
                {"source": self.edges[edge]["source"].getname(),
                 "dest": self.edges[edge]["dest"].getname(),
                 "label": self.labels[edge]}

    def _construct_jointjs_json(self):
        self.out_json["cells"] = []
        x_pos = 100
        y_pos = 100

        for component, unit_attrs in self.unit_models.items():
            try:
                self.create_image_jointjs_json(
                    self.out_json,
                    x_pos,
                    y_pos,
                    unit_attrs["name"],
                    icon_mapping[unit_attrs["type"]],
                    unit_attrs["type"],
                )
            except KeyError:
                self.create_image_jointjs_json(self.out_json, 
                                               x_pos, 
                                               y_pos, 
                                               unit_attrs["name"], 
                                               "default", unit_attrs["type"])

            x_pos += 100
            y_pos += 100

        id_counter = 0
        for name, ports_dict in self.edges.items():
            umst = self.unit_models[ports_dict["source"]]["type"]  # alias
            dest = ports_dict["dest"]

            try:
                if hasattr(ports_dict["source"], "vap_outlet"):
                    # TODO Figure out how to denote different outlet types. Need to
                    #  deal with multiple input/output offsets
                    for arc in list(self.arcs.values()):
                        if (
                            self.ports[arc.dest] == dest
                            and arc.source == ports_dict["source"].vap_outlet
                        ):
                            source_anchor = link_position_mapping[umst][
                                "top_outlet_anchor"
                            ]
                        elif (
                            self.ports[arc.dest] == dest
                            and arc.source == ports_dict["source"].liq_outlet
                        ):
                            source_anchor = link_position_mapping[umst][
                                "bottom_outlet_anchor"
                            ]

                elif "top_outlet_anchor" in link_position_mapping[umst]:
                    source_anchor = \
                        link_position_mapping[umst]["top_outlet_anchor"]
                else:
                    source_anchor = link_position_mapping[umst]["outlet_anchors"]
                    # TODO figure out offsets when mutiple things come
                    #  from/into the same side:
                    # source_anchor["args"]["dy"] = str(100/(len(dest) + 1)) + "%"

            except KeyError:
                source_anchor = link_position_mapping["default"]["outlet_anchors"]
                # TODO figure out offsets when mutiple things come from/into the 
                # same side:
                # source_anchor["args"]["dy"] = str(100/(len(dest) + 1)) + "%"
            try:
                unit_type = self.unit_models[dest]["type"]
                dest_anchor = \
                    link_position_mapping[unit_type]["inlet_anchors"]
            except KeyError:
                dest_anchor = link_position_mapping["default"]["inlet_anchors"]
            self.create_link_jointjs_json(
                self.out_json, 
                source_anchor, 
                dest_anchor, 
                ports_dict["source"].getname(), 
                dest.getname(), 
                name,
                self.labels[name]
            )
            id_counter += 1
