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
"""
Tests for dynamic utility methods.
"""

import pytest
from pyomo.environ import (ConcreteModel, Block, Constraint, Var, Set,
        TransformationFactory)
from pyomo.dae import ContinuousSet, DerivativeVar
import idaes.logger as idaeslog
from idaes.core.util.dyn_utils import *

__author__ = "Robert Parker"


# Test explicit/implicit index detection functions
def test_is_indexed_by():
    m = ConcreteModel()
    m.time = ContinuousSet(bounds=(0, 10))
    m.space = ContinuousSet(bounds=(0, 10))
    m.set = Set(initialize=['a', 'b', 'c'])
    m.v = Var()
    m.v1 = Var(m.time)
    m.v2 = Var(m.time, m.space)
    m.v3 = Var(m.set, m.space, m.time)

    @m.Block()
    def b(b):
        b.v = Var()
        b.v1 = Var(m.time)
        b.v2 = Var(m.time, m.space)
        b.v3 = Var(m.set, m.space, m.time)

    @m.Block(m.time)
    def b1(b):
        b.v = Var()
        b.v1 = Var(m.space)
        b.v2 = Var(m.space, m.set)

    @m.Block(m.time, m.space)
    def b2(b):
        b.v = Var()
        b.v1 = Var(m.set)

        @b.Block()
        def b(bl):
            bl.v = Var()
            bl.v1 = Var(m.set)
            bl.v2 = Var(m.time)

    disc = TransformationFactory('dae.collocation')
    disc.apply_to(m, wrt=m.time, nfe=5, ncp=2, scheme='LAGRANGE-RADAU')
    disc.apply_to(m, wrt=m.space, nfe=5, ncp=2, scheme='LAGRANGE-RADAU')

    assert not is_explicitly_indexed_by(m.v, m.time)
    assert is_explicitly_indexed_by(m.b.v2, m.space)

    assert not is_implicitly_indexed_by(m.v1, m.time)
    assert not is_implicitly_indexed_by(m.v2, m.set)
    assert is_implicitly_indexed_by(m.b1[m.time[1]].v2, m.time)

    assert is_implicitly_indexed_by(m.b2[m.time[1], 
        m.space[1]].b.v1, m.time)
    assert (is_implicitly_indexed_by(m.b2[m.time[1], 
        m.space[1]].b.v2, m.time) ==
        is_explicitly_indexed_by(m.b2[m.time[1], 
            m.space[1]].b.v2, m.time))
    assert (not is_implicitly_indexed_by(m.b2[m.time[1], 
        m.space[1]].b.v1, m.set))

    assert (not is_implicitly_indexed_by(m.b2[m.time[1],
        m.space[1]].b.v1, m.space, stop_at=m.b2[m.time[1], m.space[1]]))

# Test get index_set_except and _complete_index
def test_get_index_set_except():
    '''
    Tests:
      For components indexed by 0, 1, 2, 3, 4 sets:
        get_index_set_except one, then two (if any) of those sets
        check two items that should be in set_except
        insert item(s) back into these sets via index_getter
    '''
    m = ConcreteModel()
    m.time = ContinuousSet(bounds=(0, 10))
    m.space = ContinuousSet(bounds=(0, 10))
    m.set1 = Set(initialize=['a', 'b', 'c'])
    m.set2 = Set(initialize=['d', 'e', 'f'])
    m.v = Var()
    m.v1 = Var(m.time)
    m.v2 = Var(m.time, m.space)
    m.v3 = Var(m.time, m.space, m.set1)
    m.v4 = Var(m.time, m.space, m.set1, m.set2)

    # Try a multi-dimensional set
    m.set3 = Set(initialize=[('a', 1), ('b', 2)])
    m.v5 = Var(m.set3)
    m.v6 = Var(m.time, m.space, m.set3)

    disc = TransformationFactory('dae.collocation')
    disc.apply_to(m, wrt=m.time, nfe=5, ncp=2, scheme='LAGRANGE-RADAU')
    disc.apply_to(m, wrt=m.space, nfe=5, ncp=2, scheme='LAGRANGE-RADAU')

    # Want this to give a TypeError
    # info = get_index_set_except(m.v, m.time)

    # Indexed by one set
    info = get_index_set_except(m.v1, m.time)
    set_except = info['set_except']
    index_getter = info['index_getter']
    assert (set_except == [None])
    # Variable is not indexed by anything except time
    # Test that index_getter returns only the new value given,
    # regardless of whether it was part of the set excluded (time):
    assert (index_getter((), -1) == -1)

    # Indexed by two sets
    info = get_index_set_except(m.v2, m.time)
    set_except = info['set_except']
    index_getter = info['index_getter']
    assert (m.space[1] in set_except
                    and m.space.last() in set_except)
    # Here (2,) is the partial index, corresponding to space.
    # Can be provided as a scalar or tuple. 4, the time index,
    # should be inserted before (2,)
    assert (index_getter((2,), 4) == (4, 2))
    assert (index_getter(2, 4) == (4, 2))

    # Case where every set is "omitted," now for multiple sets
    info = get_index_set_except(m.v2, m.space, m.time)
    set_except = info['set_except']
    index_getter = info['index_getter']
    assert (set_except == [None])
    # 5, 7 are the desired index values for space, time 
    # index_getter should put them in the right order for m.v2,
    # even if they are not valid indices for m.v2
    assert (index_getter((), 5, 7) == (7, 5))

    # Indexed by three sets
    info = get_index_set_except(m.v3, m.time)
    # In this case set_except is a product of the two non-time sets
    # indexing v3
    set_except = info['set_except']
    index_getter = info['index_getter']
    assert ((m.space[1], 'b') in set_except
                    and (m.space.last(), 'a') in set_except)
    # index_getter inserts a scalar index into an index of length 2
    assert (index_getter((2, 'b'), 7) == (7, 2, 'b'))

    info = get_index_set_except(m.v3, m.space, m.time)
    # Two sets omitted. Now set_except is just set1
    set_except = info['set_except']
    index_getter = info['index_getter']
    assert ('a' in set_except)
    # index_getter inserts the two new indices in the right order
    assert (index_getter('b', 1.2, 1.1) == (1.1, 1.2, 'b'))

    # Indexed by four sets
    info = get_index_set_except(m.v4, m.set1, m.space)
    # set_except is a product, and there are two indices to insert
    set_except = info['set_except']
    index_getter = info['index_getter']
    assert ((m.time[1], 'd') in set_except)
    assert (index_getter((4, 'f'), 'b', 8) == (4, 8, 'b', 'f'))
    
    # The intended usage of this function looks something like:
    index_set = m.v4.index_set()
    for partial_index in set_except:
        complete_index = index_getter(partial_index, 'a', m.space[2])
        assert (complete_index in index_set)
        # Do something for every index of v4 at 'a' and space[2]

    # Indexed by a multi-dimensional set
    info = get_index_set_except(m.v5, m.set3)
    set_except = info['set_except']
    index_getter = info['index_getter']
    assert (set_except == [None])
    assert (index_getter((), ('a', 1)) == ('a', 1))

    info = get_index_set_except(m.v6, m.set3, m.time)
    set_except = info['set_except']
    index_getter = info['index_getter']
    assert (m.space[1] in set_except)
    assert (index_getter(m.space[1], ('b', 2), m.time[1])
            == (m.time[1], m.space[1], 'b', 2))


def test_fix_and_deactivate():
    m = ConcreteModel()
    m.time = ContinuousSet(bounds=(0, 10))
    m.space = ContinuousSet(bounds=(0, 5))
    m.set1 = Set(initialize=['a', 'b', 'c'])
    m.set2 = Set(initialize=['d', 'e', 'f'])
    m.fs = Block()

    m.fs.v0 = Var(m.space, initialize=1)

    @m.fs.Block()
    def b1(b):
        b.v = Var(m.time, m.space, initialize=1) 
        b.dv = DerivativeVar(b.v, wrt=m.time)

        b.con = Constraint(m.time, m.space, 
                rule=lambda b, t, x: b.dv[t, x] == 7 - b.v[t, x])

        @b.Block(m.time)
        def b2(b, t):
            b.v = Var(initialize=2)

    @m.fs.Block(m.time, m.space)
    def b2(b, t, x):
        b.v = Var(m.set1, initialize=2)

        @b.Block(m.set1)
        def b3(b, c):
            b.v = Var(m.set2, initialize=3)

            @b.Constraint(m.set2)
            def con(b, s):
                return (5*b.v[s] == 
                        m.fs.b2[m.time.first(), m.space.first()].v[c])

    @m.fs.Constraint(m.time)
    def con1(fs, t):
        return fs.b1.v[t, m.space.last()] == 5

    @m.fs.Constraint(m.space)
    def con2(fs, x):
        return fs.b1.v[m.time.first(), x] == fs.v0[x]


    disc = TransformationFactory('dae.collocation') 
    disc.apply_to(m, wrt=m.time, nfe=5, ncp=2, scheme='LAGRANGE-RADAU')
    disc.apply_to(m, wrt=m.space, nfe=5, ncp=2, scheme='LAGRANGE-RADAU')

    for t in m.time:
        m.fs.b1.v[t, m.space.first()].fix()

    active_dict = get_activity_dict(m.fs)
    for comp in m.fs.component_data_objects(Constraint, Block):
        assert active_dict[id(comp)] == True

    deactivate_model_at(m, m.time, m.time[2])
    assert m.fs.con1[m.time[1]].active
    assert not m.fs.con1[m.time[2]].active
    assert m.fs.con2[m.space[1]].active
    assert not m.fs.b1.con[m.time[2], m.space[1]].active
    assert not m.fs.b2[m.time[2], m.space.last()].active
    assert m.fs.b2[m.time[2], m.space.last()].b3['a'].con['e'].active

    deactivate_model_at(m, m.time, [m.time[1], m.time[3]], 
            outlvl=idaeslog.ERROR)
    # Higher outlvl threshold as will encounter warning trying to deactivate
    # disc equations at time.first()
    assert not m.fs.con1[m.time[1]].active
    assert not m.fs.con1[m.time[3]].active
    assert not m.fs.b1.con[m.time[1], m.space[1]].active
    assert not m.fs.b1.con[m.time[3], m.space[1]].active

    init_derivs = get_derivatives_at(m, m.time, m.time.first())[m.time.first()]
    init_deriv_names = [var.name for var in init_derivs]

    assert m.fs.b1.dv[m.time.first(), m.space.first()] in init_derivs
    assert m.fs.b1.dv[m.time[1], m.space[1]].name in init_deriv_names

    deriv_dict = get_derivatives_at(m, m.time, [m.time.first(), m.time.last()])
    deriv_name_dict = {t: [d.name for d in deriv_dict[t]] 
                          for t in deriv_dict.keys()}
    assert m.time.first() in deriv_name_dict.keys()
    assert m.time.last() in deriv_name_dict.keys()
    assert (m.fs.b1.dv[m.time.last(), m.space[1]].name 
            in deriv_name_dict[m.time.last()])
    assert (m.fs.b1.dv[m.time.last(), m.space[1]].name 
            not in deriv_name_dict[m.time.first()])

    vars_unindexed = fix_vars_unindexed_by(m, m.time)
    cons_unindexed = deactivate_constraints_unindexed_by(m, m.time)
    assert m.fs.v0[m.space[1]] in vars_unindexed
    assert m.fs.b1.b2[m.time[1]].v not in vars_unindexed
    assert m.fs.con2[m.space[2]] in cons_unindexed
    assert m.fs.con1[m.time[1]] not in cons_unindexed
    assert not m.fs.con2[m.space[1]].active
    assert m.fs.v0[m.space[1]].fixed

    path = path_from_block(m.fs.b2[m.time[1], m.space[1]].b3['a'].v,
                           m, include_comp=False)
    assert path == [('fs', None), ('b2', (m.time[1], m.space[1])),
                    ('b3', 'a')]
    path = path_from_block(m.fs.b2[m.time[1], m.space[1]].b3['a'].v,
                           m, include_comp=True)
    assert path == [('fs', None), ('b2', (m.time[1], m.space[1])),
                    ('b3', 'a'), ('v', None)]
    path = path_from_block(m.fs.b2[m.time[1], m.space[1]].b3['a'].v['f'],
                           m, include_comp=True)
    assert path == [('fs', None), ('b2', (m.time[1], m.space[1])),
                    ('b3', 'a'), ('v', 'f')]
    path = path_from_block(m.fs.b2[m.time[1], m.space[1]].b3['a'].v['f'],
                           m.fs.b2[m.time[1], m.space[1]], include_comp=True)
    assert path == [('b3', 'a'), ('v', 'f')]
    path = path_from_block(m.fs.b1.con[m.time[1], m.space[1]], m.fs)
    assert path == [('b1', None)]

    m.fs.b1.b2[m.time[1]].v.set_value(-1)
    for x in m.space:
        m.fs.b1.v[m.time[1], x].set_value(-1)
        m.fs.b1.dv[m.time[1], x].set_value(-1)
        for c1 in m.set1:
            m.fs.b2[m.time[1], x].v[c1].set_value(-1)
            for c2 in m.set2:
                m.fs.b2[m.time[1], x].b3[c1].v[c2].set_value(-1)

    copy_values_at_time(m, m, m.time.last(), m.time[1], copy_fixed=False)
    assert m.fs.b1.b2[m.time.last()].v.value == -1
    for x in m.space:
        if x != m.space.first():
            assert m.fs.b1.v[m.time.last(), x].value == -1
        else:
            assert m.fs.b1.v[m.time.last(), x].value == 1
            assert m.fs.b1.v[m.time.last(), x].fixed
        assert m.fs.b1.dv[m.time.last(), x].value == -1
        for c1 in m.set1:
            assert m.fs.b2[m.time.last(), x].v[c1].value == -1
            for c2 in m.set2:
                assert m.fs.b2[m.time[1], x].b3[c1].v[c2].value == -1

def test_copy_non_time_indexed_values():
    m1 = ConcreteModel()
    m1.time = Set(initialize=[1,2,3,4,5])
    m1.v1 = Var(m1.time, initialize=1)
    m1.v2 = Var(initialize=1)
    @m1.Block(['a','b'])
    def b1(b, i):
        b.v3 = Var(initialize=1)

        @b.Block(m1.time)
        def b2(b, t):
            b.v4 = Var(initialize=1)

        @b.Block()
        def b4(b):
            b.v6 = Var(initialize=1)
            
    @m1.Block(m1.time)
    def b3(b, t):
        b.v5 = Var(initialize=1)

    m2 = ConcreteModel()
    m2.time = Set(initialize=[1,2,3,4,5])
    m2.v1 = Var(m2.time, initialize=2)
    m2.v2 = Var(initialize=2)
    @m2.Block(['a','b'])
    def b1(b):
        b.v3 = Var(initialize=2)

        @b.Block(m2.time)
        def b2(b, t):
            b.v4 = Var(initialize=2)

        @b.Block()
        def b4(b):
            b.v6 = Var(initialize=2)

    @m2.Block(m1.time)
    def b3(b, t):
        b.v5 = Var(initialize=2)

    copy_non_time_indexed_values(m1, m2)
    assert m1.v1[1].value != m2.v1[1].value
    assert m1.v2.value == m2.v2.value == 2
    assert m1.b1['a'].v3.value == m2.b1['a'].v3.value == 2
    assert m1.b1['b'].b4.v6.value == m2.b1['b'].b4.v6.value == 2
    assert m1.b3[3].v5.value != m2.b3[3].v5.value


if __name__ == "__main__":
    test_is_indexed_by()
    test_get_index_set_except()
    test_fix_and_deactivate()
    test_copy_non_time_indexed_values()
