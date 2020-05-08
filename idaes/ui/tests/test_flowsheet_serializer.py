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
import pytest
from operator import itemgetter

from idaes.generic_models.flowsheets.demo_flowsheet import (
    build_flowsheet, set_dof, initialize_flowsheet, solve_flowsheet)
from idaes.ui.flowsheet_serializer import FlowsheetSerializer
from idaes.ui.icon_mapping import icon_mapping
from idaes.ui.link_position_mapping import link_position_mapping


def test_serialize_flowsheet():
    # Construct the model from idaes/examples/workshops/Module_2_Flowsheet/Module_2_Flowsheet_Solution.ipynb
    m = build_flowsheet()
  
    fss = FlowsheetSerializer()
    fss.serialize_flowsheet(m.fs)

    unit_models = fss.get_unit_models()
    unit_model_names_types = []
    for unit_model in unit_models:
        unit_model_names_types.append(unit_models[unit_model])

    unit_models_names_type_truth = [{'name': 'M01', 'type': 'mixer'}, 
                                    {'name': 'H02', 'type': 'heater'}, 
                                    {'name': 'F03', 'type': 'flash'}]
  
    set_result = set(tuple(sorted(d.items())) for d in unit_model_names_types)
    set_truth = set(tuple(sorted(d.items())) for d in unit_models_names_type_truth)    
    difference = list(set_truth.symmetric_difference(set_result))

    assert len(difference) == 0
  
    # TODO Figure out how to test ports. Maybe find out if we can find the parent component for the port?
    # ports = fss.get_ports()
    # assert ports == {"<pyomo.network.port.SimplePort object at 0x7fe8d0d79278>": "<idaes.core.process_block._ScalarMixer object at 0x7fe8d0d60360>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0d792e8>": "<idaes.core.process_block._ScalarMixer object at 0x7fe8d0d60360>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0d79358>": "<idaes.core.process_block._ScalarMixer object at 0x7fe8d0d60360>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0d793c8>": "<idaes.core.process_block._ScalarMixer object at 0x7fe8d0d60360>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0d797b8>": "<idaes.core.process_block._ScalarHeater object at 0x7fe8d0db74c8>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0d79828>": "<idaes.core.process_block._ScalarHeater object at 0x7fe8d0db74c8>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0d79a58>": "<idaes.core.process_block._ScalarStoichiometricReactor object at 0x7fe8d0de2ab0>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0d79ac8>": "<idaes.core.process_block._ScalarStoichiometricReactor object at 0x7fe8d0de2ab0>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0d79eb8>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8d0e0fdc8>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0e41128>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8d0e0fdc8>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0e41198>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8d0e0fdc8>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0d79f98>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8d0e0fdc8>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0e41048>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8d0e0fdc8>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0e410b8>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8d0e0fdc8>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0e41278>": "<idaes.core.process_block._ScalarSeparator object at 0x7fe8d0e45708>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0e41588>": "<idaes.core.process_block._ScalarSeparator object at 0x7fe8d0e45708>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0e415f8>": "<idaes.core.process_block._ScalarSeparator object at 0x7fe8d0e45708>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0e41828>": "<idaes.core.process_block._ScalarPressureChanger object at 0x7fe8d0e686c0>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0e41898>": "<idaes.core.process_block._ScalarPressureChanger object at 0x7fe8d0e686c0>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0e41c88>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8e1405cf0>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0e41eb8>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8e1405cf0>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0e41f28>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8e1405cf0>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0e41e48>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8e1405cf0>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0e41dd8>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8e1405cf0>", 
    #                  "<pyomo.network.port.SimplePort object at 0x7fe8d0e41d68>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8e1405cf0>"
    #                  }

    named_edges_results = {}
    edges = fss.get_edges()

    for name, end_points in edges.items():
      named_edges_results[name] = {"source": end_points["source"].getname(), "dest": end_points["dest"].getname()}

    named_edges_truth = {'s01': {'source': 'M01', 'dest': 'H02'}, 
                         's02': {'source': 'H02', 'dest': 'F03'}}

    assert named_edges_results == named_edges_truth


def test_create_image_jointjs_json():
  out_json = {"model": {}}
  out_json["cells"] = []
  x_pos = 0
  y_pos = 0
  component_id = "M101"
  component_type = "mixer"
  image = icon_mapping(component_type)

  fss = FlowsheetSerializer()
  fss.create_image_jointjs_json(out_json, x_pos, y_pos, component_id, image, component_type)
  assert out_json == {'cells': 
                        [{
                            'type': 'standard.Image', 
                            'position': {'x': 0, 'y': 0}, 
                            'size': {'width': 50, 'height': 50}, 
                            'angle': 0, 'id': 'M101', 'z': (1,), 
                            'attrs': {
                                'image': {'xlinkHref': '/images/icons/mixer.svg'}, 
                                'label': {'text': 'M101'}, 
                                'root': {'title': 'mixer'}
                          }
                        }],
                       "model": {}
                      }


def test_create_link_jointjs_json():
    out_json = {"model": {}}
    out_json["cells"] = []
    source_anchor = link_position_mapping["heater"]["outlet_anchors"]
    dest_anchor = link_position_mapping["mixer"]["inlet_anchors"]
    source_id = "M101"
    dest_id = "F101"
    link_id = "s03"
    label = "foo"

    fss = FlowsheetSerializer()
    fss.create_link_jointjs_json(out_json, source_anchor, dest_anchor, source_id, dest_id, link_id, label)
    assert out_json == {'cells': [{
                            'type': 'standard.Link', 
                            'source': {
                              'anchor': {
                                'name': 'right', 
                                'args': {
                                  'rotate': 'false', 
                                  'padding': 0
                                }
                              }, 
                              'id': 'M101'
                            },
                            'target': {
                              'anchor': {
                                'name': 'left', 
                                'args': {
                                  'rotate': 'false', 
                                  'padding': 0
                                }
                              }, 
                              'id': 'F101'
                            }, 
                            'router': {
                              'name': 'orthogonal', 
                              'padding': 10
                            }, 
                            'connector': {
                              'name': 'normal', 
                              'attrs': {
                                'line': {
                                  'stroke': '#5c9adb'
                                }
                              }
                            }, 
                            'id': 's03', 
                            'labels': [{
                              'attrs': {
                                'rect': {
                                  "fill": "#d7dce0", 
                                  "stroke": "#FFFFFF", 
                                  'stroke-width': 1
                                }, 
                                'text': {
                                  'text': 'foo', 
                                  'fill': 'black', 
                                  'text-anchor': 'left'
                                }
                              }, 
                              'position': {
                                'distance': 0.66, 
                                'offset': -40
                              }
                            }], 
                            'z': 2
                          }],
                          "model": {}
                        }
