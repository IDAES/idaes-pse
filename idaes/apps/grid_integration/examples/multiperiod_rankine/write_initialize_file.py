import pyomo.environ as pyo
from idaes.apps.rankine.simple_rankine_cycle import create_model, set_inputs, initialize_model, close_flowsheet_loop, add_operating_cost
from idaes.core.util import to_json

m = pyo.ConcreteModel()
m.rankine = create_model(heat_recovery=True)
m.rankine = set_inputs(m.rankine)
m.rankine = initialize_model(m.rankine)
to_json(m.rankine, fname="initialized_rankine_state.json.gz",gz=True, human_read=True)
