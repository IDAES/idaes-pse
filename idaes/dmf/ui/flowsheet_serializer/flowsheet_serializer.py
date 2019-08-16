from collections import defaultdict
import json
import os

from idaes.core import UnitModelBlockData
from pyomo.network.port import SimplePort
from pyomo.network.arc import SimpleArc
from collections import defaultdict
from idaes.dmf.ui import icon_mapping, link_position_mapping


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
                #self.unit_models[component] = {"name": component.getname(), "type": type(component).__name__}
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
            
            # saved for later
#        for port in orphaned_ports:
#            named_edges["Orphaned"].append(self.ports[port].getname())
#        edges_filename = file_base_name + 'edges.json'
#        with open(edges_filename, 'w') as outfile:  
#            json.dump(named_edges, outfile)
#            
#        named_components = {}
#        for comp in self.unit_models.values():
#            named_components[comp["name"]]= {"type": comp["type"]}
#        nodes_filename = file_base_name + 'nodes.json'
#        with open(nodes_filename, 'w') as outfile:  
#            json.dump(named_components, outfile)


        vis_file_name = file_base_name + '.idaes.vis'
        if os.path.isfile(vis_file_name) and overwrite == False:
            msg = (f"{vis_file_name} already exists. If you wish to overwrite this file call save() with overwrite=True. "
                   "WARNING: If you overwrite the file, you will lose your saved layout.")
            raise FileBaseNameExistsError(msg)
        else:
            outjson = {}
            print(f"Creating {vis_file_name}")
            #outjson['joint_enc'] = []

        #outjson['model']['cells'] = []
        outjson['cells'] = []
        x_pos = 50
        y_pos = 50

        for component, unit_attrs in self.unit_models.items():
            print(unit_attrs['type'])
            print(unit_attrs)
            print(component)
            entry = {}
            entry['type'] = 'standard.Image'
            # for now, just tile the positions diagonally
            entry['position'] = {'x': x_pos, 'y': y_pos} # TODO Make the default positioning better
            x_pos += 100
            y_pos += 50
            entry['size'] = {"width": 100, "height": 100} # TODO Set the width and height depending on the icon rather than default
            entry['angle'] = 0
            entry['id'] = unit_attrs['name']
            entry['z'] = 1,
            entry['attrs'] = {
                    'image': {'xlinkHref': icon_mapping[unit_attrs['type']]},
                    'label': {'text': unit_attrs['name']},
                    'root': {'title': 'joint.shapes.standard.Image'}
                    } # TODO, pending image displaying solution
            #outjson['model']['cells'].append(entry)
            outjson['cells'].append(entry)

        id_counter = 0
        for source, dests in edges.items():
            for dest in dests:
                if "vapor_outlet_anchor" in link_position_mapping[self.unit_models[source]["type"]] and self.unit_models[dest]["type"] != "flash":
                    # TODO Deal with output being on top and on bottom of the flash. This will include figuring out how to denote different outlet types
                    # Need to deal with multiple input/output offsets
                    source_anchor = link_position_mapping[self.unit_models[source]["type"]]["vapor_outlet_anchor"]
                elif "liquid_outlet_anchor" in link_position_mapping[self.unit_models[source]["type"]] and self.unit_models[dest]["type"] == "flash":
                    source_anchor = link_position_mapping[self.unit_models[source]["type"]]["liquid_outlet_anchor"]
                else:
                    source_anchor = link_position_mapping[self.unit_models[source]["type"]]["outlet_anchors"]
                    #if "dy" not in source_anchor["args"]:
                        #source_anchor["args"]["dy"] = link_position_mapping[self.unit_models[source]["type"]]["outlet_anchors"]["args"]["dy"] = str(100/(len(dests) + 1)) + "%"

                print(self.unit_models[source]["type"])
                print(source_anchor)
                print(self.unit_models[dest]["type"])
                print(link_position_mapping[self.unit_models[dest]["type"]]["inlet_anchors"])
                entry = {
                        'type': 'standard.Link',
                        'source': {"anchor": source_anchor, 'id': source.getname()},
                        'target': {"anchor": link_position_mapping[self.unit_models[dest]["type"]]["inlet_anchors"], 'id': dest.getname()},
                        "router": {"name": "orthogonal", "padding": 10},
                        "connector": {"name": "jumpover", 'attrs': {'line': {'stroke': '#6FB1E1'}}},
                        'id': id_counter,# TODO
                        'z': 2,
                        #'labels': [], # sic
                        }
                outjson['cells'].append(entry)
                id_counter += 1

        # TODO handle also reading preexisting file first; is there an elegant way?
        with open(vis_file_name, 'w') as outfile:
            json.dump(outjson, outfile)
