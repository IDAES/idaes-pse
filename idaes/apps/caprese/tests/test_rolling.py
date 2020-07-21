import pytest
from collections import OrderedDict
from pyomo.environ import *
from idaes.apps.caprese.util import cuid_from_timeslice
from idaes.apps.caprese.rolling import *

def make_model():
    m = ConcreteModel()
    m.time = Set(initialize=range(11))
    m.space = Set(initialize=range(3))
    m.comp = Set(initialize=['a','b'])

    m.v1 = Var(m.time, m.space)
    m.v2 = Var(m.time, m.space, m.comp)

    @m.Block(m.time, m.space)
    def b(b, t, x):
        b.v3 = Var()
        b.v4 = Var(m.comp)
    
    m.v1_refs = [Reference(m.v1[:,x]) for x in m.space]
    m.v2a_refs = [Reference(m.v2[:,x,'a']) for x in m.space]
    m.v2b_refs = [Reference(m.v2[:,x,'b']) for x in m.space]
    m.v3_refs = [Reference(m.b[:,x].v3) for x in m.space]
    m.v4a_refs = [Reference(m.b[:,x].v4['a']) for x in m.space]
    m.v4b_refs = [Reference(m.b[:,x].v4['a']) for x in m.space]

    return m

@pytest.mark.unit
def test_history_init():
    m = make_model()

    # Cases to test:
    #   1. synchronous with good data
    #   2. synchronous with bad (asynchronous) data
    #   3. asynchronous with synchronous data (should be fine)
    #   4. asynchronous with asynchronous data
    
    # 1.
    data = {cuid_from_timeslice(_slice, m.time): [float(i) for i in m.time] 
        for _slice in m.v1_refs}

    history = History(init=data, time=m.time, name='v1_history')

    for cuid in data:
        assert cuid in history
        assert history[cuid] == [float(i) for i in m.time]
        assert history[cuid].time == list(m.time)

    # 4.
    cuid_a = [cuid_from_timeslice(_slice, m.time) for _slice in m.v2a_refs]
    cuid_b = [cuid_from_timeslice(_slice, m.time) for _slice in m.v2b_refs]
    data = {cuid: [(i, float(i)) for i in m.time] for cuid in cuid_a}
    data.update({cuid: [(i, 2*float(i)) for i in m.time if i % 2 == 0]
        for cuid in cuid_b})

    history = History(init=data, time=None, name='v2_history')

    assert history.time is None
    for cuid in cuid_a:
        assert cuid in history
        assert history[cuid] == [float(i) for i in m.time]
        assert history[cuid].time == [i for i in m.time]
    for cuid in cuid_b:
        assert cuid in history
        assert history[cuid] == [2*float(i) for i in m.time if i % 2 == 0]
        assert history[cuid].time == [i for i in m.time if i % 2 == 0]

    # 3.
    data = {cuid_from_timeslice(_slice, m.time): [(i, 3*float(i)) for i in m.time]
            for _slice in m.v3_refs}
    history = History(init=data, time=None, name='v3_history')
    
    assert history.time is None
    for cuid in data:
        assert cuid in history
        assert history[cuid] == [3*float(i) for i in m.time]
        assert history[cuid].time == list(m.time)

    # 2.
    cuid_a = [cuid_from_timeslice(_slice, m.time) for _slice in m.v4a_refs]
    cuid_b = [cuid_from_timeslice(_slice, m.time) for _slice in m.v4b_refs]
    data = {}
    data.update({cuid: [(i, 1.5*float(i)) for i in m.time if i % 3 == 0]})
    data.update({cuid: [(i, 2.5*float(i)) for i in m.time if i % 5 == 0]})

    with pytest.raises(ValueError, match=r'.*must have same length.*'):
        history = History(init=data, time=m.time, name='v4_history')
