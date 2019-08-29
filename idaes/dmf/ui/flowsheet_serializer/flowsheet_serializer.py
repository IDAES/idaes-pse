from collections import defaultdict
import json
import os

from idaes.core import UnitModelBlockData
from idaes.dmf.ui import icon_mapping, link_position_mapping

from pyomo.environ import Block
from pyomo.network.port import SimplePort
from pyomo.network.arc import SimpleArc

class FileBaseNameExistsError(Exception):
    pass


class FlowsheetSerializer:
    def __init__(self):
        self.unit_models = {}
        self.arcs = []
        self.ports = {}

    def save(self, flowsheet, file_base_name, overwrite=False):
        """
        Serializes the flowsheet and saves it to a file that can be read by the idaes-model-vis 
        jupyter lab extension.
        :param flowsheet: The flowsheet to save. Usually fetched from the model.
        :param file_base_name: The file prefix to the .idaes.vis file produced. The file is created/saved
        in the directory that you ran from Jupyter Lab.
        :param overwrite: Boolean to overwrite an existing file_base_name.idaes.vis.
        If True, the existing file with the same file_base_name will be overwritten. This will cause you to lose
        any saved layout. 
        If False and there is an existing file with that file_base_name, you will get an error 
        message stating that you cannot save a file to the file_base_name (and therefore overwriting the saved 
        layout). If there is not an existing file with that file_base_name then it saves as normal.
        Defaults to False.
        :return: None

        Usage example:
        '''
            m = ConcreteModel()
            m.fs = FlowsheetBlock(...)
            ...
            serializer = FlowsheetSerializer()
            serializer.save(m.fs, 'output_file')
        '''
        """
        for component in flowsheet.component_objects(descend_into=False):
            # TODO try using component_objects(ctype=X)
            if isinstance(component, UnitModelBlockData):
                self.unit_models[component] = {"name": component.getname(), "type": component._orig_module.split('.')[-1]}
                
                for subcomponent in component.component_objects(descend_into=False):
                    if isinstance(subcomponent, SimplePort):
                       self.ports[subcomponent] = component
                        
            elif isinstance(component, SimpleArc): 
               self.arcs.append(component)
            
        edges = defaultdict(list)
        orphaned_ports = set(self.ports.keys())
        for arc in self.arcs:
            edges[self.ports[arc.source]].append(self.ports[arc.dest])
            orphaned_ports.discard(arc.source)
            orphaned_ports.discard(arc.dest)

        vis_file_name = file_base_name + '.idaes.vis'
        if os.path.isfile(vis_file_name) and overwrite == False:
            msg = (f"{vis_file_name} already exists. If you wish to overwrite this file call save() with overwrite=True. "
                   "WARNING: If you overwrite the file, you will lose your saved layout.")
            raise FileBaseNameExistsError(msg)
        else:
            outjson = {}
            print(f"Creating {vis_file_name}")

        outjson['cells'] = []
        x_pos = 50
        y_pos = 50

        for component, unit_attrs in self.unit_models.items():
            print(type(unit_attrs['name']))
            print(unit_attrs['name'])
            self.create_image_json(outjson, x_pos, y_pos, unit_attrs['name'], 
                icon_mapping[unit_attrs['type']], unit_attrs['name'], unit_attrs['type'])
            x_pos += 50
            y_pos += 50

        id_counter = 0
        for source, dests in edges.items():
            for dest in dests:
                if hasattr(source, "vap_outlet"):
                    # TODO Figure out how to denote different outlet types. Need to deal with multiple input/output offsets
                    for arc in self.arcs:
                        if self.ports[arc.dest] == dest and arc.source == source.vap_outlet:
                            source_anchor = link_position_mapping[self.unit_models[source]["type"]]["top_outlet_anchor"]
                        elif self.ports[arc.dest] == dest and arc.source == source.liq_outlet:
                            source_anchor = link_position_mapping[self.unit_models[source]["type"]]["bottom_outlet_anchor"]

                elif "top_outlet_anchor" in link_position_mapping[self.unit_models[source]["type"]]:
                    source_anchor = link_position_mapping[self.unit_models[source]["type"]]["top_outlet_anchor"]
                else:
                    source_anchor = link_position_mapping[self.unit_models[source]["type"]]["outlet_anchors"]
                    # TODO figure out offsets when mutiple things come from/into the same side:
                    # source_anchor["args"]["dy"] = str(100/(len(dests) + 1)) + "%"

                dest_anchor = link_position_mapping[self.unit_models[dest]["type"]]["inlet_anchors"]

                self.create_link_json(outjson, source_anchor, dest_anchor, source.getname(), dest.getname(), id_counter)
                id_counter += 1

        num_open_inlets = 0
        for orphan_port in orphaned_ports:
            unit_model_name = ""
            try:
                if num_open_inlets <= 0:
                    num_open_inlets = len(self.ports[orphan_port].create_inlet_list()) - len(edges[self.ports[orphan_port]])
            except AttributeError:
                num_open_inlets = 0

            if num_open_inlets > 0:
                icon_type = "feed"
            else:
                icon_type = "product"

            self.create_image_json(outjson, x_pos, y_pos, "inlet" + str(id_counter), icon_mapping[icon_type], "", icon_type)
            x_pos += 50
            y_pos += 50

            source_anchor = link_position_mapping[icon_type]["outlet_anchors"]
            dest_anchor = link_position_mapping[self.unit_models[self.ports[orphan_port]]["type"]]["inlet_anchors"]

            if icon_type == "feed":
                source_anchor = link_position_mapping[icon_type]["outlet_anchors"]
                dest_anchor = link_position_mapping[self.unit_models[self.ports[orphan_port]]["type"]]["inlet_anchors"]
                source_id = "inlet" + str(id_counter)
                dest_id = self.ports[orphan_port].getname()
            else:
                source_anchor = link_position_mapping[self.unit_models[self.ports[orphan_port]]["type"]]["inlet_anchors"]
                dest_anchor = link_position_mapping[icon_type]["outlet_anchors"]
                source_id = self.ports[orphan_port].getname()
                dest_id = "inlet" + str(id_counter)

            self.create_link_json(outjson, source_anchor, dest_anchor, source_id, dest_id, id_counter)

            id_counter += 1
            num_open_inlets -= 1

        # TODO handle also reading preexisting file first; is there an elegant way?
        with open(vis_file_name, 'w') as outfile:
            json.dump(outjson, outfile)


    def create_image_json(self, outjson, x_pos, y_pos, id, image, label, title):
        entry = {}
        entry['type'] = 'standard.Image'
        # for now, just tile the positions diagonally
        entry['position'] = {'x': x_pos, 'y': y_pos} # TODO Make the default positioning better
        entry['size'] = {"width": 50, "height": 50} # TODO Set the width and height depending on the icon rather than default
        entry['angle'] = 0
        entry['id'] = id
        entry['z'] = 1,
        entry['attrs'] = {
                'image': {'xlinkHref': image},
                'label': {'text': label},
                'root': {'title': title}
                }
        #outjson['model']['cells'].append(entry)
        outjson['cells'].append(entry)


    def create_link_json(self, outjson, source_anchor, dest_anchor, source_id, dest_id, link_id):
        entry = {
                'type': 'standard.Link',
                'source': {"anchor": source_anchor, 'id': source_id},
                'target': {"anchor": dest_anchor, 'id': dest_id},
                "router": {"name": "orthogonal", "padding": 10},
                "connector": {"name": "jumpover", 'attrs': {'line': {'stroke': '#6FB1E1'}}},
                'id': link_id,
                'z': 2,
                #'labels': [],
                }
        print(entry)
        outjson['cells'].append(entry)