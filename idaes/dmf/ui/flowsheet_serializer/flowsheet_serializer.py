from idaes.core import UnitModelBlockData
from pyomo.network.port import SimplePort
from pyomo.network.arc import SimpleArc
from collections import defaultdict
import json
from idaes.dmf.ui import icon_mapping

class FlowsheetSerializer:
    def __init__(self):
        self.unit_models = {}
        self.arcs = []
        self.ports = {}

    def save(self, flowsheet, file_base_name):
        for component in flowsheet.component_objects(descend_into=False):
            # TODO try using component_objects(ctype=X)
            if isinstance(component, UnitModelBlockData):
                self.unit_models[component] = {"name": component.getname(), "type": type(component).__name__}
                
                for subcomponent in component.component_objects(descend_into=False):
                    if isinstance(subcomponent, SimplePort):
                       self.ports[subcomponent] = component
                        
            elif isinstance(component, SimpleArc): 
               self.arcs.append(component)
            
        edges = {}
        orphaned_ports = set(self.ports.keys())
        for arc in self.arcs:
            edges[(self.ports[arc.source], self.ports[arc.dest])] = arc
            orphaned_ports.discard(arc.source)
            orphaned_ports.discard(arc.dest)


        labeled_edges = defaultdict(list)

        for (source, dest) in edges:
            labeled_edges[source.getname()].append(dest.getname())
            
            # saved for later
#        for port in orphaned_ports:
#            labeled_edges["Orphaned"].append(self.ports[port].getname())
#        edges_filename = file_base_name + 'edges.json'
#        with open(edges_filename, 'w') as outfile:  
#            json.dump(labeled_edges, outfile)
#            
#        named_components = {}
#        for comp in self.unit_models.values():
#            named_components[comp["name"]]= {"type": comp["type"]}
#        nodes_filename = file_base_name + 'nodes.json'
#        with open(nodes_filename, 'w') as outfile:  
#            json.dump(named_components, outfile)


        vis_file_name = file_base_name + '.idaes.vis'
        try:
            with open(vis_file_name, 'r') as oldfile:
                outjson = json.load(oldfile)
                print('loading existing saved vis file')
        except FileNotFoundError:
            outjson = {}
            print('creating new vis file')
            #outjson['joint_enc'] = []

        #outjson['model']['cells'] = []
        outjson['cells'] = []
        x_pos = 50
        y_pos = 50

        for component, unit_attrs in self.unit_models.items():
            entry = {}
            entry['type'] = 'standard.Image'
            # for now, just tile the positions diagonally
            entry['position'] = {'x': x_pos, 'y': y_pos} # TODO
            x_pos += 100
            y_pos += 50
            entry['size'] = {"width": 100, "height": 100} # TODO
            entry['angle'] = 0
            entry['id'] = unit_attrs['name']
            entry['z'] = 1,
            entry['attrs'] = {
                    'image': {'xlinkHref': icon_mapping[unit_attrs['type']]},
                    #'label': {'text':""},
                    'root': {'title': 'joint.shapes.standard.Image'}
                    } # TODO, pending image displaying solution
            #outjson['model']['cells'].append(entry)
            outjson['cells'].append(entry)

        id_counter = 0
        for source, dests in labeled_edges.items():
            for dest in dests:
                entry = {
                        'type': 'standard.Link',
                        'source': {'id': source},
                        'target': {'id': dest},
                        'id': id_counter,# TODO
                        'z': 2,
                        #'labels': [], # sic
                        'attrs': {'line': {'stroke': '#6FB1E1'}}
                        }
                outjson['cells'].append(entry)
                id_counter += 1



        # TODO handle also reading preexisting file first; is there an elegant way?
        with open(vis_file_name, 'w') as outfile:
            json.dump(outjson, outfile)
