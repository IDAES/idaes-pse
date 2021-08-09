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
'''
Tests for dynamic utility functions

Author: Robert Parker
'''
from pyomo.environ import (ConcreteModel, Var, Block, Set,
        TransformationFactory)
from pyomo.dae import (ContinuousSet)

from dyn_utils import *

m = ConcreteModel()

m.time = ContinuousSet(bounds=(0,10))
m.x = ContinuousSet(bounds=(0,10))
m.s = Set(initialize=['a', 'b', 'c'])

m.v = Var(initialize=-1)
m.vt = Var(m.time, initialize=0)
m.vtx = Var(m.time, m.x, initialize=1)
m.vtxs = Var(m.time, m.x, m.s, initialize=2)


disc = TransformationFactory('dae.collocation')
disc.apply_to(m, wrt=m.time, nfe=5, ncp=2, scheme='LAGRANGE-RADAU')
disc.apply_to(m, wrt=m.x, nfe=5, ncp=3, scheme='LAGRANGE-RADAU')

tfe_list = m.time.get_finite_elements()
xfe_list = m.x.get_finite_elements()

m2 = ConcreteModel()
m2.v = Var(initialize=-2)
m2.time = Set(initialize=[0,1,2])

# NOTE: There is another bug somewhere regarding 'equality' of 
# continuous sets. The btx block does not get properly constructed
# when x and t have the same elements.
# 
# May have been fixed by Bethany's most recent dae updates...
m.bt = Block(m.time)
m.btx = Block(m.time, m.x)

for t in m.time:
    m.bt[t].vx = Var(m.x, initialize=1)
    m.bt[t].vxs = Var(m.x, m.s, initialize=2)

for tx in m.btx.index_set():
    m.btx[tx].v = Var(initialize=0)
    m.btx[tx].vs = Var(m.s, initialize=1)
###

def test_get_ik_from_index():
    assert get_ik_from_index(m.time, 4) == (2,0)
    assert get_ik_from_index(m.time, 8.666667) == (4,1)

def test_is_indexed_by():
    assert is_indexed_by(m.vt, m.time) and is_indexed_by(m.vtx, m.time)
    assert not is_indexed_by(m.v, m.time) and not is_indexed_by(m.vt, m.x)
    assert is_indexed_by(m.btx, m.x) and not is_indexed_by(m.bt[0].vx, m.time)

def test_is_implicitly_indexed_by():
    assert is_implicitly_indexed_by(m.vtx, m.x)
    assert is_implicitly_indexed_by(m.bt[tfe_list[2]].vx, m.time)
    assert not is_implicitly_indexed_by(m.bt[tfe_list[2]].vx, m.s)
    assert is_implicitly_indexed_by(m.btx[tfe_list[2], xfe_list[1]].v, m.x)

def test_get_index_set_except():
    info = get_index_set_except(m.vt, m.time)
    set_except = info['set_except']
    index_getter = info['index_getter']
    assert index_getter((), 2.72) == 2.72
    assert None in set_except

    info = get_index_set_except(m.vtx, m.time)
    set_except = info['set_except']
    index_getter = info['index_getter']
    assert (index_getter(xfe_list[2], 1j) == index_getter((xfe_list[2],), 1j) 
            == (1j, 4))
    assert set_except is m.x

    info = get_index_set_except(m.vtxs, m.time)
    set_except = info['set_except']
    index_getter = info['index_getter']
    assert index_getter((xfe_list[1], 'b'), 3.14) == (3.14, xfe_list[1], 'b')
    assert set_except == m.x.cross(m.s)

    info = get_index_set_except(m.btx, m.time)
    set_except = info['set_except']
    index_getter = info['index_getter']
    assert set_except is m.x
    assert index_getter(xfe_list[1], 1) == (1, xfe_list[1])

def test_path_to():
    assert (path_to(m.btx[tfe_list[1], xfe_list[1]].vs['a']) == 
            [('btx', (tfe_list[1], xfe_list[1]))])
    assert (path_to(m.btx[tfe_list[1], xfe_list[1]].vs['a'], include_comp=True) 
            == [('btx', (tfe_list[1], xfe_list[1])), ('vs', 'a')])

def test_path_from_block():
    t, x, c = tfe_list[1], xfe_list[2], 'c'

    m.btx[t, x].blocky = Block()
    blocky = m.btx[t, x].blocky
    blocky.blocko = Block(m.s)
    blocko = blocky.blocko[c]
    blocko.v = Var(initialize=8)
    assert (path_from_block(blocko.v, m.btx[t, x]) 
            == [('blocky', None), ('blocko', c)])
    assert (path_from_block(blocko.v, m, include_comp=True)
            == [('btx', (t, x)), ('blocky', None), ('blocko', c), ('v', None)])

    # Asymmetric structure may cause problems for assumptions made by 
    # copy_values functions. Don't want to introduced this complexity. 
    blocko.del_component(blocko.v)
    blocky.del_component(blocko)
    m.btx[t, x].del_component(blocky)

def test_copy_non_time_indexed_values():
    m.v.set_value(-1)
    m2.v.set_value(-2)
    copy_non_time_indexed_values(m2, m)
    assert m.v.value == m2.v.value

def test_copy_values_at_time():
    t_tgt = tfe_list[1]
    t_src = tfe_list[0]
    m.vt[t_src] = 1
    m.vt[t_tgt] = 2
    for x in m.x:
        m.bt[t_src].vx[x] = 1
        m.bt[t_tgt].vx[x] = 2
    copy_values_at_time(m, m, t_tgt, t_src)
    assert m.vt[t_src] == m.vt[t_tgt] == 1
    assert ([m.bt[t_src].vx[x].value for x in m.x] ==
            [m.bt[t_tgt].vx[x].value for x in m.x])

if __name__ == '__main__':
    test_get_ik_from_index()
    test_is_indexed_by()
    test_is_implicitly_indexed_by()
    test_get_index_set_except()
    test_path_to()
    test_path_from_block()
    test_copy_non_time_indexed_values()
    test_copy_values_at_time()
