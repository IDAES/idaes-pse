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

import hda_ideal_VLE as thermo_props
import hda_reaction as reaction_props

from pyomo.environ import (Constraint,
                           Var,
                           ConcreteModel,
                           Expression,
                           Objective,
                           SolverFactory,
                           TransformationFactory,
                           value)
from pyomo.network import Arc, SequentialDecomposition

from idaes.core import FlowsheetBlock

from idaes.dmf.ui.flowsheet_serializer import FlowsheetSerializer
from idaes.dmf.ui.icon_mapping import icon_mapping
from idaes.dmf.ui.link_position_mapping import link_position_mapping

from idaes.generic_models.unit_models import (Flash,
                               PressureChanger,
                               Mixer,
                               Separator as Splitter,
                               Heater,
                               StoichiometricReactor)
from idaes.generic_models.unit_models.pressure_changer import ThermodynamicAssumption


def test_serialize_flowsheet():
    # Construct the model from idaes/examples/workshops/Module_2_Flowsheet/Module_2_Flowsheet_Solution.ipynb
    m = ConcreteModel()
    m.fs = FlowsheetBlock(default={"dynamic": False})
    m.fs.thermo_params = thermo_props.HDAParameterBlock()
    m.fs.reaction_params = reaction_props.HDAReactionParameterBlock(default={"property_package": m.fs.thermo_params})
  
    m.fs.M101 = Mixer(default={"property_package": m.fs.thermo_params,
                               "inlet_list": ["toluene_feed", "hydrogen_feed", "vapor_recycle"]
                               })
  
    m.fs.H101 = Heater(default={"property_package": m.fs.thermo_params,
                                "has_pressure_change": False,
                                "has_phase_equilibrium": True})
    m.fs.R101 = StoichiometricReactor(default={"property_package": m.fs.thermo_params,
                                               "reaction_package": m.fs.reaction_params,
                                               "has_heat_of_reaction": True,
                                               "has_heat_transfer": True,
                                               "has_pressure_change": False
                                               })
    m.fs.F101 = Flash(default={"property_package": m.fs.thermo_params,
                               "has_heat_transfer": True,
                               "has_pressure_change": True
                               })
    m.fs.S101 = Splitter(default={"property_package": m.fs.thermo_params,
                                  "ideal_separation": False,
                                  "outlet_list": ["purge", "recycle"]
                                  })
    m.fs.C101 = PressureChanger(default={"property_package": m.fs.thermo_params,
                                         "compressor": True,
                                         "thermodynamic_assumption": ThermodynamicAssumption.isothermal
                                         })
    m.fs.F102 = Flash(default={"property_package": m.fs.thermo_params,
                               "has_heat_transfer": True,
                               "has_pressure_change": True
                               })
  
    m.fs.s03 = Arc(source=m.fs.M101.outlet, destination=m.fs.H101.inlet)
    m.fs.s04 = Arc(source=m.fs.H101.outlet, destination=m.fs.R101.inlet)
    m.fs.s05 = Arc(source=m.fs.R101.outlet, destination=m.fs.F101.inlet)
    m.fs.s06 = Arc(source=m.fs.F101.vap_outlet, destination=m.fs.S101.inlet)
    m.fs.s08 = Arc(source=m.fs.S101.recycle, destination=m.fs.C101.inlet)
    m.fs.s09 = Arc(source=m.fs.C101.outlet, destination=m.fs.M101.vapor_recycle)
    m.fs.s10 = Arc(source=m.fs.F101.liq_outlet, destination=m.fs.F102.inlet)
  
    fss = FlowsheetSerializer()
    fss.serialize_flowsheet(m.fs)

    unit_models = fss.get_unit_models()
    unit_model_names_types = []
    for unit_model in unit_models:
        unit_model_names_types.append(unit_models[unit_model])

    unit_models_names_type_truth = [{'name': 'M101', 'type': 'mixer'}, 
                                    {'name': 'H101', 'type': 'heater'}, 
                                    {'name': 'R101', 'type': 'stoichiometric_reactor'}, 
                                    {'name': 'F101', 'type': 'flash'}, 
                                    {'name': 'S101', 'type': 'separator'}, 
                                    {'name': 'C101', 'type': 'pressure_changer'}, 
                                    {'name': 'F102', 'type': 'flash'}]
  
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

    print(named_edges_results)

    named_edges_truth = {'s03': {'source': 'M101', 'dest': 'H101'}, 
                         's04': {'source': 'H101', 'dest': 'R101'}, 
                         's05': {'source': 'R101', 'dest': 'F101'}, 
                         's06': {'source': 'F101', 'dest': 'S101'}, 
                         's08': {'source': 'S101', 'dest': 'C101'}, 
                         's09': {'source': 'C101', 'dest': 'M101'}, 
                         's10': {'source': 'F101', 'dest': 'F102'}}

    assert named_edges_results == named_edges_truth


def test_create_image_jointjs_json():
  out_json = {"model": {}}
  out_json["cells"] = []
  x_pos = 0
  y_pos = 0
  component_id = "M101"
  component_type = "mixer"
  image = icon_mapping[component_type]

  fss = FlowsheetSerializer()
  fss.create_image_jointjs_json(out_json, x_pos, y_pos, component_id, image, component_type)
  assert out_json == {'cells': 
                        [{
                            'type': 'standard.Image', 
                            'position': {'x': 0, 'y': 0}, 
                            'size': {'width': 50, 'height': 50}, 
                            'angle': 0, 'id': 'M101', 'z': (1,), 
                            'attrs': {
                                'image': {'xlinkHref': 'mixer.svg'}, 
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
