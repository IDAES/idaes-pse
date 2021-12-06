import pyomo.environ as pyo
from idaes.apps.grid_integration.examples.multiperiod_rankine.multiperiod_rankine_cycle import create_multiperiod_rankine_model

blk = pyo.Block()
blk.construct()

horizon = 4
mp_rankine = create_multiperiod_rankine_model(n_time_points=horizon)
m = mp_rankine.pyomo_model
blk.rankine = m

blk.HOUR = pyo.Set(initialize = range(horizon))
for t in blk.HOUR:
	print(t)

#Create expression that references underlying power variables
blk.HOUR = pyo.Set(initialize = range(horizon))
blk.P_T = pyo.Expression(blk.HOUR)
blk.tot_cost = pyo.Expression(blk.HOUR)
for (t,b) in enumerate(mp_rankine.get_active_process_blocks()):
    blk.P_T[t] = b.rankine.P_to_grid + b.battery.discharge
    blk.tot_cost[t] = b.rankine.fs.operating_cost

blk.P_T.index_set()
blk.tot_cost.index_set()


