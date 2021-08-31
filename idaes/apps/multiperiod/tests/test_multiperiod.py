import pytest
import pyomo.environ as pyo
from idaes.apps.multiperiod.multiperiod import MultiPeriodModel
import numpy as np

def create_test_model(p):
    m = pyo.ConcreteModel()

    m.x1 = pyo.Var(within=pyo.NonNegativeReals)
    m.y1 = pyo.Var(within=pyo.NonNegativeReals)
    m.x2 = pyo.Var(within=pyo.NonNegativeReals)
    m.y2 = pyo.Var(within=pyo.NonNegativeReals)

    m.c1 = pyo.Constraint(expr = m.x1+m.y1==2)
    m.c2 = pyo.Constraint(expr = m.x1+m.x2<=p)
    m.c3 = pyo.Constraint(expr = m.y1+m.y2>=p)
    return m

def get_link_variables(b1,b2):
    return [(b1.x2,b2.x1),(b1.y2,b2.y1)]

def get_periodic_variables(b1,b2):
    return [(b1.y2,b2.y1),]

@pytest.mark.unit
def test_build_model():
    #data passed to each block time period
    n_time_points = 4
    time_points = np.arange(0,n_time_points)
    time_data = np.linspace(0.5,2,4)
    data_points = [{"p":time_data[i]} for i in range(n_time_points)]
    data_kwargs = dict(zip(time_points,data_points))

    mp = MultiPeriodModel(n_time_points, create_test_model, get_link_variables,
    periodic_variable_func=get_periodic_variables)

    assert mp.current_time == None
    assert mp.n_time_points == 4
    assert mp.create_process_model == create_test_model
    assert mp.get_linking_variable_pairs == get_link_variables
    assert mp.get_periodic_variable_pairs == get_periodic_variables

    #create the multiperiod object
    mp.build_multi_period_model(data_kwargs)
    m = mp.pyomo_model

    assert len(m.TIME) == 4
    assert m.nvariables() == 16
    assert m.nconstraints() == 12 + 6 + 1

    blks = mp.get_active_process_blocks()
    assert len(blks) == 4

@pytest.mark.unit
def test_advance_time():
    n_time_points = 4
    time_points = np.arange(0,n_time_points)
    time_data = np.linspace(0.5,2,4)
    data_points = [{"p":time_data[i]} for i in range(n_time_points)]
    data_kwargs = dict(zip(time_points,data_points))

    mp = MultiPeriodModel(n_time_points, create_test_model, get_link_variables,
    periodic_variable_func=get_periodic_variables)
    mp.build_multi_period_model(data_kwargs)
    m = mp.pyomo_model

    assert(mp.current_time == 0)
    assert len(m.TIME) == 4
    assert(len(m.blocks) == 4)
    mp.advance_time(p=3)
    assert(mp.current_time == 1)
    assert(len(m.blocks) == 5)
    assert(m.blocks[0].process.active == False)

    blks = mp.get_active_process_blocks()
    assert len(blks) == 4
    assert all([blk.active for blk in blks])

    assert m.nvariables() == 16
    assert m.nconstraints() == 12 + 6 + 1
