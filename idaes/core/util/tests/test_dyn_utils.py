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
from idaes.core.util.dyn_utils import *

#from idaes.core import FlowsheetBlock
#from idaes.core.util.testing import PhysicalParameterTestBlock
#from idaes.core.util.exceptions import ConfigurationError
#from idaes.core.util.initialization import (fix_state_vars,
#                                            revert_state_vars,
#                                            propagate_state,
#                                            solve_indexed_blocks)

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
    self.assertTrue('a' in set_except)
    # index_getter inserts the two new indices in the right order
    self.assertEqual(index_getter('b', 1.2, 1.1), (1.1, 1.2, 'b'))

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

if __name__ == "__main__":
    test_is_indexed_by()
