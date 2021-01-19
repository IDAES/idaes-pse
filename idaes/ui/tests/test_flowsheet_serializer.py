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
import pytest
from operator import itemgetter

import re

from idaes.generic_models.flowsheets.demo_flowsheet import (
    build_flowsheet, set_dof, initialize_flowsheet, solve_flowsheet)
from idaes.generic_models.unit_models.pressure_changer import \
    ThermodynamicAssumption
from idaes.generic_models.properties.swco2 import SWCO2ParameterBlock
from idaes.generic_models.unit_models import Heater, PressureChanger
from idaes.ui.flowsheet_serializer import FlowsheetSerializer
from idaes.ui.icon_mapping import icon_mapping
from idaes.ui.link_position_mapping import link_position_mapping

from pyomo.environ import Expression


@pytest.mark.component
def test_serialize_flowsheet():
    # Construct the model from idaes/examples/workshops/Module_2_Flowsheet/Module_2_Flowsheet_Solution.ipynb
    m = build_flowsheet()
    m.fs.properties = SWCO2ParameterBlock()
    m.fs.main_compressor = PressureChanger(
      default={'dynamic': False,
               'property_package': m.fs.properties,
               'compressor': True,
               'thermodynamic_assumption': ThermodynamicAssumption.isentropic})

    m.fs.bypass_compressor = PressureChanger(
        default={'dynamic': False,
                 'property_package': m.fs.properties,
                 'compressor': True,
                 'thermodynamic_assumption': ThermodynamicAssumption.isentropic})

    m.fs.turbine = PressureChanger(
      default={'dynamic': False,
               'property_package': m.fs.properties,
               'compressor': False,
               'thermodynamic_assumption': ThermodynamicAssumption.isentropic})
    m.fs.boiler = Heater(default={'dynamic': False,
                                  'property_package': m.fs.properties,
                                  'has_pressure_change': True})
    m.fs.FG_cooler = Heater(default={'dynamic': False,
                                     'property_package': m.fs.properties,
                                     'has_pressure_change': True})
    m.fs.pre_boiler = Heater(default={'dynamic': False,
                                      'property_package': m.fs.properties,
                                      'has_pressure_change': False})
    m.fs.HTR_pseudo_tube = Heater(default={'dynamic': False,
                                       'property_package': m.fs.properties,
                                       'has_pressure_change': True})
    m.fs.LTR_pseudo_tube = Heater(default={'dynamic': False,
                                       'property_package': m.fs.properties,
                                       'has_pressure_change': True})

    # _set_numerical_details(m)

    fss = FlowsheetSerializer()
    fss.serialize(m.fs, 'myflowsheet')

    unit_models = fss.get_unit_models()
    unit_model_names_types = []
    for unit_model in unit_models:
        unit_model_names_types.append(unit_models[unit_model])

    unit_models_names_type_truth = [{'name': 'M01', 'type': 'mixer'},
                                    {'name': 'H02', 'type': 'heater'},
                                    {'name': 'F03', 'type': 'flash'},
                                    {'name': 'main_compressor', 'type': 'pressure_changer'},
                                    {'name': 'bypass_compressor', 'type': 'pressure_changer'},
                                    {'name': 'turbine', 'type': 'pressure_changer'},
                                    {'name': 'boiler', 'type': 'heater'},
                                    {'name': 'FG_cooler', 'type': 'heater'},
                                    {'name': 'pre_boiler', 'type': 'heater'},
                                    {'name': 'HTR_pseudo_tube', 'type': 'heater'},
                                    {'name': 'LTR_pseudo_tube', 'type': 'heater'},
                                    ]

    # TODO: Examine whether these are in fact appropriate
    inlet_names = {'inlet', 'inlet_1_1', 'inlet_2_1'}.union({f'inlet_{n}' for n in range(1, 8)})
    outlet_names = {'outlet', 'vap_outlet', 'liq_outlet'}.union({f'outlet_{n}' for n in range(1, 8)})
    feed_and_outlet_names_type_truth = [{'name': name, 'type': 'feed'} for name in inlet_names] + \
                                       [{'name': name, 'type': 'product'} for name in outlet_names]
    unit_models_names_type_truth += feed_and_outlet_names_type_truth

    # canonicalize, because lists of dictionaries are awkward to compare
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

    # TODO: as obove, examine whether this is appropriate treatment for arcs connected to inlets/outlets
    # In particular, note that numbered inlet/outlet names are somewhat arbitrary;
    # they're only deterministic due to sorting upon automatic generation.
    in_out_edges_truth = {
           's_inlet': {'dest': 'FG_cooler', 'source': 'inlet'},
           's_inlet_1_1': {'dest': 'M01', 'source': 'inlet_1_1'},
           's_inlet_1': {'dest': 'HTR_pseudo_tube', 'source': 'inlet_1'},
           's_inlet_2_1': {'dest': 'M01', 'source': 'inlet_2_1'},
           's_inlet_2': {'dest': 'LTR_pseudo_tube', 'source': 'inlet_2'},
           's_inlet_3': {'dest': 'boiler', 'source': 'inlet_3'},
           's_inlet_4': {'dest': 'bypass_compressor', 'source': 'inlet_4'},
           's_inlet_5': {'dest': 'main_compressor', 'source': 'inlet_5'},
           's_inlet_6': {'dest': 'pre_boiler', 'source': 'inlet_6'},
           's_inlet_7': {'dest': 'turbine', 'source': 'inlet_7'},
           's_liq_outlet': {'dest': 'liq_outlet', 'source': 'F03'},
           's_outlet': {'dest': 'outlet', 'source': 'FG_cooler'},
           's_outlet_1': {'dest': 'outlet_1', 'source': 'HTR_pseudo_tube'},
           's_outlet_2': {'dest': 'outlet_2', 'source': 'LTR_pseudo_tube'},
           's_outlet_3': {'dest': 'outlet_3', 'source': 'boiler'},
           's_outlet_4': {'dest': 'outlet_4', 'source': 'bypass_compressor'},
           's_outlet_5': {'dest': 'outlet_5', 'source': 'main_compressor'},
           's_outlet_6': {'dest': 'outlet_6', 'source': 'pre_boiler'},
           's_outlet_7': {'dest': 'outlet_7', 'source': 'turbine'},
           's_vap_outlet': {'dest': 'vap_outlet', 'source': 'F03'}}

    for edge, details in in_out_edges_truth.items():
        named_edges_truth[edge] = details

    assert named_edges_results == named_edges_truth


def _set_numerical_details(m):
    m.fs.turbine.ratioP.fix(1 / 3.68)
    m.fs.turbine.efficiency_isentropic.fix(0.927)
    m.fs.turbine.initialize()
    m.gross_cycle_power_output = \
        Expression(expr=(-m.fs.turbine.work_mechanical[0] -
                         m.fs.main_compressor.work_mechanical[0] -
                         m.fs.bypass_compressor.work_mechanical[0]))
    # account for generator loss = 1.5% of gross power output
    m.net_cycle_power_output = Expression(expr=0.985 * m.gross_cycle_power_output)
    m.total_cycle_power_input = Expression(
        expr=(m.fs.boiler.heat_duty[0] + m.fs.pre_boiler.heat_duty[0] +
              m.fs.FG_cooler.heat_duty[0]))
    m.cycle_efficiency = Expression(
        expr=m.net_cycle_power_output / m.total_cycle_power_input * 100)
    # Expression to compute recovered duty in recuperators
    m.recuperator_duty = Expression(
        expr=(m.fs.HTR_pseudo_tube.heat_duty[0] +
              m.fs.LTR_pseudo_tube.heat_duty[0]))


@pytest.mark.unit
def test_in_out_regex_matching():
    inlet_base_names = ['in', 'feed', 'inlet']
    outlet_base_names = ['out', 'prod', 'outlet']
    inlet_matches = inlet_base_names + \
                    [base + '_toluene' for base in inlet_base_names] + \
                    ['toluene_' + base for base in inlet_base_names] + \
                    ['feed_outlet']  # NOTE: this might be considered pathological input
    outlet_matches = outlet_base_names + \
                    [base + '_hydrogen' for base in outlet_base_names] + \
                    ['hydrogen_' + base for base in outlet_base_names]
    unmatches = [
        '',
        'trout',
        'internal',
        'benzoin_vap',
        'somerandomtext',
        'fresh_produce'
    ]
    pattern = FlowsheetSerializer.INLET_REGEX
    for name in inlet_matches:
        assert pattern.search(name)
    for name in unmatches:
        assert pattern.search(name) is None

    pattern = FlowsheetSerializer.OUTLET_REGEX
    for name in outlet_matches:
        assert pattern.search(name)
    for name in unmatches:
        assert pattern.search(name) is None

@pytest.mark.unit
def test__unique_unit_name():
    fss = FlowsheetSerializer()
    fss._used_unit_names = {'inlet', 'inlet_1', 'outlet', 'toluene_prod'}
    assert fss._unique_unit_name('hydrogen_in') == 'hydrogen_in'
    assert fss._unique_unit_name('toluene') == 'toluene'
    assert fss._unique_unit_name('outlet') == 'outlet_1'
    assert fss._unique_unit_name('inlet') == 'inlet_2'
    assert fss._unique_unit_name('prod') == 'prod'

@pytest.mark.unit
def test_create_image_jointjs_json():
  out_json = {"model": {}}
  out_json["cells"] = []
  x_pos = 0
  y_pos = 0
  component_id = "M101"
  component_type = "mixer"
  image = icon_mapping(component_type)
  port_group = {"groups": {
                    "in":{
                        "position":{
                            "name":"left",
                            "args":{
                               "x":2,
                               "y":25,
                               "dx":1,
                               "dy":1
                            }
                        },
                        "attrs": {
                            "rect": {
                                "stroke": '#000000',
                                'stroke-width': 0,
                                "width": 0,
                                "height": 0
                            }
                        },
                        "markup": '<g><rect/></g>'
                    }
                },
                "items":[
                   {
                      "group":"in",
                      "id":"in"
                   },
                   {
                      "group":"out",
                      "id":"out"
                   }
                ]}

  fss = FlowsheetSerializer()
  fss.create_image_jointjs_json(out_json, x_pos, y_pos, component_id, image, component_type, port_group)
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
                            },
                            'ports': {
                                "groups": {
                                "in":{
                                    "position":{
                                        "name":"left",
                                        "args":{
                                           "x":2,
                                           "y":25,
                                           "dx":1,
                                           "dy":1
                                        }
                                    },
                                    "attrs": {
                                        "rect": {
                                            "stroke": '#000000',
                                            'stroke-width': 0,
                                            "width": 0,
                                            "height": 0
                                        }
                                    },
                                    "markup": '<g><rect/></g>'
                                }
                              },
                              "items":[
                                 {
                                    "group":"in",
                                    "id":"in"
                                 },
                                 {
                                    "group":"out",
                                    "id":"out"
                                 }
                              ]
                            }
                        }],
                       "model": {}
                      }

@pytest.mark.unit
def test_create_link_jointjs_json():
    out_json = {"model": {}}
    out_json["cells"] = []
    source_port = "out"
    dest_port = "in"
    source_id = "M101"
    dest_id = "F101"
    link_id = "s03"
    label = "foo"

    fss = FlowsheetSerializer()
    fss.create_link_jointjs_json(out_json, source_port, dest_port, source_id, dest_id, link_id, label)
    assert out_json == {'cells': [{
                            'type': 'standard.Link',
                            'source': {
                              'port': 'out',
                              'id': 'M101'
                            },
                            'target': {
                              'port': 'in',
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
                              "attrs": {
                                  "rect": {
                                    "fill": '#d7dce0',
                                    "stroke": 'white',
                                    'stroke-width': 0,
                                    "fill-opacity": 0},
                                  "text": {
                                      "text": label,
                                      "fill": 'black',
                                      'text-anchor': 'left',
                                      "display": "none"
                                  },
                              },
                              "position": {
                                  "distance": 0.66,
                                  "offset": -40
                              }},
                              {"attrs": {
                                  "text": {
                                      "text": link_id
                                  }
                              }}],
                            'z': 2
                          }],
                          "model": {}
                        }
