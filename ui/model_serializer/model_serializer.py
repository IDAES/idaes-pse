from idaes.core import UnitModelBlockData
from pyomo.network.port import SimplePort
from pyomo.network.arc import SimpleArc
from collections import defaultdict
import json

class ModelSerializer:
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
            
        for port in orphaned_ports:
            labeled_edges["Orphaned"].append(self.ports[port].getname())
        edges_filename = file_base_name + 'edges.json'
        with open(edges_filename, 'w') as outfile:  
            json.dump(labeled_edges, outfile)
            
        named_components = {}
        for comp in self.unit_models.values():
            named_components[comp["name"]]= {"type": comp["type"]}
        nodes_filename = file_base_name + 'nodes.json'
        with open(nodes_filename, 'w') as outfile:  
            json.dump(named_components, outfile)
