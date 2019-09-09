import pytest

import idaes.dmf.ui.tests.resources.hda_ideal_VLE as thermo_props
import idaes.dmf.ui.tests.resources.hda_reaction as reaction_props

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

from idaes.unit_models import (Flash,
                               PressureChanger,
                               Mixer,
                               Separator as Splitter,
                               Heater,
                               StoichiometricReactor)
from idaes.unit_models.pressure_changer import ThermodynamicAssumption


# We expect this test to fail until we figure out how to test this
@pytest.mark.xfail
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
    unit_models, ports, edges = fss.serialize_flowsheet(m.fs)
  
    # Figure out how to test this because the output of fss.serialize_flowsheet has memory locations
    assert unit_models == {"<idaes.core.process_block._ScalarMixer object at 0x7fe8d0d60360>": {"name": "M101", "type": "mixer"}, 
                           "<idaes.core.process_block._ScalarHeater object at 0x7fe8d0db74c8>": {"name": "H101", "type": "heater"}, 
                           "<idaes.core.process_block._ScalarStoichiometricReactor object at 0x7fe8d0de2ab0>": {"name": "R101", "type": "stoichiometric_reactor"}, 
                           "<idaes.core.process_block._ScalarFlash object at 0x7fe8d0e0fdc8>": {"name": "F101", "type": "flash"}, 
                           "<idaes.core.process_block._ScalarSeparator object at 0x7fe8d0e45708>": {"name": "S101", "type": "separator"}, 
                           "<idaes.core.process_block._ScalarPressureChanger object at 0x7fe8d0e686c0>": {"name": "C101", "type": "pressure_changer"}, 
                           "<idaes.core.process_block._ScalarFlash object at 0x7fe8e1405cf0>": {"name": "F102", "type": "flash"}
                           }
  
    assert ports == {"<pyomo.network.port.SimplePort object at 0x7fe8d0d79278>": "<idaes.core.process_block._ScalarMixer object at 0x7fe8d0d60360>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0d792e8>": "<idaes.core.process_block._ScalarMixer object at 0x7fe8d0d60360>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0d79358>": "<idaes.core.process_block._ScalarMixer object at 0x7fe8d0d60360>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0d793c8>": "<idaes.core.process_block._ScalarMixer object at 0x7fe8d0d60360>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0d797b8>": "<idaes.core.process_block._ScalarHeater object at 0x7fe8d0db74c8>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0d79828>": "<idaes.core.process_block._ScalarHeater object at 0x7fe8d0db74c8>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0d79a58>": "<idaes.core.process_block._ScalarStoichiometricReactor object at 0x7fe8d0de2ab0>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0d79ac8>": "<idaes.core.process_block._ScalarStoichiometricReactor object at 0x7fe8d0de2ab0>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0d79eb8>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8d0e0fdc8>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0e41128>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8d0e0fdc8>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0e41198>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8d0e0fdc8>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0d79f98>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8d0e0fdc8>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0e41048>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8d0e0fdc8>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0e410b8>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8d0e0fdc8>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0e41278>": "<idaes.core.process_block._ScalarSeparator object at 0x7fe8d0e45708>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0e41588>": "<idaes.core.process_block._ScalarSeparator object at 0x7fe8d0e45708>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0e415f8>": "<idaes.core.process_block._ScalarSeparator object at 0x7fe8d0e45708>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0e41828>": "<idaes.core.process_block._ScalarPressureChanger object at 0x7fe8d0e686c0>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0e41898>": "<idaes.core.process_block._ScalarPressureChanger object at 0x7fe8d0e686c0>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0e41c88>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8e1405cf0>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0e41eb8>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8e1405cf0>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0e41f28>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8e1405cf0>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0e41e48>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8e1405cf0>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0e41dd8>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8e1405cf0>", 
                     "<pyomo.network.port.SimplePort object at 0x7fe8d0e41d68>": "<idaes.core.process_block._ScalarFlash object at 0x7fe8e1405cf0>"
                     }
    assert edges == {"<idaes.core.process_block._ScalarMixer object at 0x7fe8d0d60360>": ["<idaes.core.process_block._ScalarHeater object at 0x7fe8d0db74c8>"], 
                     "<idaes.core.process_block._ScalarHeater object at 0x7fe8d0db74c8>": ["<idaes.core.process_block._ScalarStoichiometricReactor object at 0x7fe8d0de2ab0>"], 
                     "<idaes.core.process_block._ScalarStoichiometricReactor object at 0x7fe8d0de2ab0>": ["<idaes.core.process_block._ScalarFlash object at 0x7fe8d0e0fdc8>"], 
                     "<idaes.core.process_block._ScalarFlash object at 0x7fe8d0e0fdc8>": 
                     ["<idaes.core.process_block._ScalarSeparator object at 0x7fe8d0e45708>", 
                     "<idaes.core.process_block._ScalarFlash object at 0x7fe8e1405cf0>"], 
                     "<idaes.core.process_block._ScalarSeparator object at 0x7fe8d0e45708>": ["<idaes.core.process_block._ScalarPressureChanger object at 0x7fe8d0e686c0>"], 
                     "<idaes.core.process_block._ScalarPressureChanger object at 0x7fe8d0e686c0>": ["<idaes.core.process_block._ScalarMixer object at 0x7fe8d0d60360>"]
                     }


def test_create_image_json():
  out_json = {}
  out_json["cells"] = []
  x_pos = 0
  y_pos = 0
  component_id = "M101"
  component_type = "mixer"
  image = icon_mapping[component_type]
  label = "M101"

  fss = FlowsheetSerializer()
  fss.create_image_json(out_json, x_pos, y_pos, component_id, image, label, component_type)
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
                      }]
                     }


def test_create_link_json():
    out_json = {}
    out_json["cells"] = []
    source_anchor = link_position_mapping["heater"]["outlet_anchors"]
    dest_anchor = link_position_mapping["mixer"]["inlet_anchors"]
    source_id = "M101"
    dest_id = "F101"
    link_id = 1

    fss = FlowsheetSerializer()
    fss.create_link_json(out_json, source_anchor, dest_anchor, source_id, dest_id, link_id)
    print(out_json)
    assert out_json == {'cells': 
                        [{
                            'type': 'standard.Link', 
                            'source': {
                                'anchor': {
                                            'name': 'right', 
                                            'args': {'rotate': 'false', 'padding': 0}
                                          }, 
                                'id': 'M101'
                            }, 
                            'target': {
                                'anchor': {
                                            'name': 'left', 
                                            'args': {'rotate': 'false', 'padding': 0}
                                          }, 
                                'id': 'F101'
                            }, 
                            'router': {'name': 'orthogonal', 'padding': 10}, 
                            'connector': {
                                             'name': 'jumpover', 
                                             'attrs': {'line': {'stroke': '#6FB1E1'}}
                                         }, 
                            'id': 1, 'z': 2
                        }]
                       }
