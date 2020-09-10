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
Tests for Caprese helper utility functions.
"""

import pytest
from pytest import approx
from pyomo.environ import (Block, ConcreteModel,  Constraint, Expression,
                           Set, SolverFactory, Var, value, Objective,
                           TransformationFactory, TerminationCondition,
                           Reference)
from pyomo.network import Arc
from pyomo.common.collections import ComponentSet, ComponentMap
from pyomo.core.expr.visitor import identify_variables
from pyomo.dae.flatten import flatten_dae_components

from idaes.core import (FlowsheetBlock, MaterialBalanceType, EnergyBalanceType,
        MomentumBalanceType)
from idaes.core.util.model_statistics import (degrees_of_freedom, 
        activated_equalities_generator, unfixed_variables_generator)
from idaes.core.util.initialization import initialize_by_time_element
from idaes.core.util.exceptions import ConfigurationError
from idaes.generic_models.unit_models import CSTR, Mixer, MomentumMixingType
from idaes.apps.caprese.util import *
from idaes.apps.caprese.examples.cstr_model import make_model
import idaes.logger as idaeslog

import random

__author__ = "Robert Parker"


# See if ipopt is available and set up solver
solver_available = SolverFactory('ipopt').available()
if solver_available:
    solver = SolverFactory('ipopt')
    solver.options = {'tol': 1e-6,
                      'mu_init': 1e-8,
                      'bound_push': 1e-8,
                      'halt_on_ampl_error': 'yes'}
else:
    solver = None


@pytest.mark.unit
def test_find_comp_in_block():
    m1 = ConcreteModel()

    @m1.Block([1,2,3])
    def b1(b):
        b.v = Var([1,2,3])

    m2 = ConcreteModel()

    @m2.Block([1,2,3])
    def b1(b):
        b.v = Var([1,2,3,4])

    @m2.Block([1,2,3])
    def b2(b):
        b.v = Var([1,2,3])

    v1 = m1.b1[1].v[1]

    assert find_comp_in_block(m2, m1, v1) is m2.b1[1].v[1]

    v2 = m2.b2[1].v[1]
    v3 = m2.b1[3].v[4]

    # These should result in Attribute/KeyErrors
    #find_comp_in_block(m1, m2, v2)
    #find_comp_in_block(m1, m2, v3)
    assert find_comp_in_block(m1, m2, v2, allow_miss=True) is None
    assert find_comp_in_block(m1, m2, v3, allow_miss=True) is None


@pytest.mark.unit
def test_NMPCVarLocator():
    m = ConcreteModel()
    m.time = Set(initialize=[1,2,3])
    m.v = Var(m.time, ['a','b','c'])

    varlist = [Reference(m.v[:,'a']),
               Reference(m.v[:,'b']),
               Reference(m.v[:,'b'])]
    group = NMPCVarGroup(varlist, m.time)

    categ = VariableCategory.DIFFERENTIAL
    locator = NMPCVarLocator(categ, group, 0, is_ic=True)

    assert locator.category == VariableCategory.DIFFERENTIAL
    assert locator.group is group
    assert locator.location == 0
    assert locator.is_ic == True


@pytest.mark.unit
def test_copy_values():
    # Define m1
    m1 = ConcreteModel()
    m1.time = Set(initialize=[1,2,3,4,5])

    m1.v1 = Var(m1.time, initialize=1)
    
    @m1.Block(m1.time)
    def blk(b, t):
        b.v2 = Var(initialize=1)

    # Define m2
    m2 = ConcreteModel()
    m2.time = Set(initialize=[1,2,3,4,5])

    m2.v1 = Var(m2.time, initialize=2)
    
    @m2.Block(m2.time)
    def blk(b, t):
        b.v2 = Var(initialize=2)

    ###

    scalar_vars_1, dae_vars_1 = flatten_dae_components(m1, m1.time, Var)
    scalar_vars_2, dae_vars_2 = flatten_dae_components(m2, m2.time, Var)

    m2.v1[2].set_value(5)
    m2.blk[2].v2.set_value(5)

    copy_values_at_time(dae_vars_1, dae_vars_2, 1, 2)

    for t in m1.time:
        if t != 1:
            assert m1.v1[t].value == 1
            assert m1.blk[t].v2.value == 1
        else:
            assert m1.v1[t].value == 5
            assert m1.blk[t].v2.value == 5


@pytest.mark.unit
def test_find_slices_in_model():
    # Define m1
    m1 = ConcreteModel()
    m1.time = Set(initialize=[1,2,3,4,5])

    m1.v1 = Var(m1.time, initialize=1)
    
    @m1.Block(m1.time)
    def blk(b, t):
        b.v2 = Var(initialize=1)

    # Define m2
    m2 = ConcreteModel()
    m2.time = Set(initialize=[1,2,3,4,5])

    m2.v1 = Var(m2.time, initialize=2)
    
    @m2.Block(m2.time)
    def blk(b, t):
        b.v2 = Var(initialize=2)

    ###

    scalar_vars_1, dae_vars_1 = flatten_dae_components(m1, m1.time, Var)
    scalar_vars_2, dae_vars_2 = flatten_dae_components(m2, m2.time, Var)

    t0_tgt = m1.time.first()
    group = NMPCVarGroup(dae_vars_1, m1.time)
    categ = VariableCategory.ALGEBRAIC
    locator = ComponentMap([(var[t0_tgt], NMPCVarLocator(categ, group, i))
                                for i, var in enumerate(dae_vars_1)])

    tgt_slices = find_slices_in_model(m1, m1.time, m2, m2.time, 
            locator, dae_vars_2)

    dae_var_set_1 = ComponentSet(dae_vars_1)
    assert len(dae_var_set_1) == len(tgt_slices)
    assert len(tgt_slices) == len(dae_vars_2)
    for i, _slice in enumerate(tgt_slices):
        assert dae_vars_2[i].name == _slice.name
        assert _slice in dae_var_set_1


@pytest.mark.skipif(solver is None, reason="Solver not available")
@pytest.mark.component
def test_initialize_by_element_in_range():
    mod = make_model(horizon=2, ntfe=20)
    assert degrees_of_freedom(mod) == 0

    scalar_vars, dae_vars = flatten_dae_components(mod.fs, mod.fs.time, Var)
    diff_vars = [Reference(mod.fs.cstr.control_volume.energy_holdup[:, 'aq']),
                 Reference(mod.fs.cstr.control_volume.material_holdup[:, 'aq', 'S']),
                 Reference(mod.fs.cstr.control_volume.material_holdup[:, 'aq', 'E']),
                 Reference(mod.fs.cstr.control_volume.material_holdup[:, 'aq', 'C']),
                 Reference(mod.fs.cstr.control_volume.material_holdup[:, 'aq', 'P'])]

    initialize_by_element_in_range(mod.fs, mod.fs.time, 0, 1, solver=solver, 
                        dae_vars=dae_vars,
                        time_linking_variables=diff_vars,
                        outlvl=idaeslog.DEBUG,
                        solve_initial_conditions=True)

    assert degrees_of_freedom(mod.fs) == 0

    assert mod.fs.cstr.outlet.conc_mol[1, 'S'].value == approx(10.189, abs=1e-3)
    assert mod.fs.cstr.outlet.conc_mol[1, 'C'].value == approx(0.4275, abs=1e-4)
    assert mod.fs.cstr.outlet.conc_mol[1, 'E'].value == approx(0.0541, abs=1e-4)
    assert mod.fs.cstr.outlet.conc_mol[1, 'P'].value == approx(0.3503, abs=1e-4)

    initialize_by_element_in_range(mod.fs, mod.fs.time, 1, 2, solver=solver, 
                        dae_vars=dae_vars,
                        outlvl=idaeslog.DEBUG)

    assert degrees_of_freedom(mod.fs) == 0
    for con in activated_equalities_generator(mod.fs):
        assert value(con.body) - value(con.upper) < 1e-5

    assert mod.fs.cstr.outlet.conc_mol[2, 'S'].value == approx(11.263, abs=1e-3)
    assert mod.fs.cstr.outlet.conc_mol[2, 'C'].value == approx(0.4809, abs=1e-4)
    assert mod.fs.cstr.outlet.conc_mol[2, 'E'].value == approx(0.0538, abs=1e-4)
    assert mod.fs.cstr.outlet.conc_mol[2, 'P'].value == approx(0.4372, abs=1e-4)


@pytest.mark.component
def test_add_noise_at_time():
    mod = make_model(horizon=2, ntfe=20)
    time = mod.fs.time
    t0 = time.first()
    assert degrees_of_freedom(mod) == 0

    scalar_vars, dae_vars = flatten_dae_components(mod.fs, time, Var)
    diff_vars = [Reference(mod.fs.cstr.control_volume.energy_holdup[:, 'aq']),
                 Reference(mod.fs.cstr.control_volume.material_holdup[:, 'aq', 'S']),
                 Reference(mod.fs.cstr.control_volume.material_holdup[:, 'aq', 'E']),
                 Reference(mod.fs.cstr.control_volume.material_holdup[:, 'aq', 'C']),
                 Reference(mod.fs.cstr.control_volume.material_holdup[:, 'aq', 'P'])]

    for t in time:
        diff_vars[0][t].setlb(290)
        diff_vars[0][t].setub(310)
        for i in range(1,5):
            diff_vars[i][t].setlb(0)
            diff_vars[i][t].setub(1)
            # Pretend this is mole fraction...

    assert diff_vars[0][0].value == 300
    for i in range(1,5):
        assert diff_vars[i][0].value == 0.001

    copy_values_at_time(diff_vars,
                        diff_vars,
                        [t for t in time if t != t0],
                        t0)

    for seed in [4, 8, 15, 16, 23, 42]:
        random.seed(seed)
        weights = [10, 0.001, 0.001, 0.001, 0.001]
        nom_vals = add_noise_at_time(diff_vars, 0, weights=weights)

        assert nom_vals[0][0] == 300
        assert diff_vars[0][0].value != 300
        assert diff_vars[0][0].value == approx(300, abs=2)
        for i in range(1,5):
            assert nom_vals[0][i] == 0.001
            # ^ nom_vals indexed by time, then var-index. This is confusing,
            # might need to change (or only accept one time point at a time)
            assert diff_vars[i][0].value != 0.001
            assert diff_vars[i][0].value == approx(0.001, abs=2e-4)
            # Within four standard deviations should be a safe check

        for i in range(0, 5):
            diff_vars[i][0].set_value(nom_vals[0][i])
        # Reset and try again with new seed

    # Try providing function for uniform random
    rand_fcn = random.uniform

    random_arg_dict = {'range_list': [(295, 305),
                                      (0.001, 0.01),
                                      (0.001, 0.01),
                                      (0.001, 0.01),
                                      (0.001, 0.01)]}

    def args_fcn(i, val, **kwargs):
        # args_fcn expects arguments like this
        range_list = kwargs.pop('range_list', None)
        return range_list[i]

    nom_vals = add_noise_at_time(diff_vars, 0.5, 
                                 random_function=rand_fcn,
                                 args_function=args_fcn,
                                 random_arg_dict=random_arg_dict)

    assert nom_vals[0.5][0] == 300
    assert diff_vars[0][0.5].value != 300
    assert 295 <= diff_vars[0][0.5].value <= 305
    for i in range(1, 5):
        assert nom_vals[0.5][i] == 0.001
        assert diff_vars[i][0.5].value != 0.001
        assert 0.001 <= diff_vars[i][0.5].value <= 0.01

    # Try to get some bound violations
    random_arg_dict = {'range_list': [(295, 305),
                                      (1, 2),
                                      (1, 2),
                                      (1, 2),
                                      (1, 2)]}
    
    nom_vals = add_noise_at_time(diff_vars, 1,
                                 random_function=rand_fcn,
                                 args_function=args_fcn,
                                 random_arg_dict=random_arg_dict,
                                 bound_strategy='push',
                                 bound_push=0.01)
    
    for i in range(1, 5):
        assert diff_vars[i][1].value == 0.99

    random.seed(123)
    with pytest.raises(ValueError) as exc_test:
        # Large weights - one of these lower bounds should fail...
        nom_vals = add_noise_at_time(diff_vars, 1.5, 
                                     bound_strategy='discard',
                                     discard_limit=0,
                                     weights=[1,1,1,1,1],
                                     sig_0=0.05)

    @pytest.mark.unit
    def test_get_violated_bounds_at_time():
        m = ConcreteModel()
        m.time = Set(initialize=[1,2,3])
        m.v = Var(m.time, ['a','b','c'], initialize=5)

        varlist = [Reference(m.v[:,'a']),
                   Reference(m.v[:,'b']),
                   Reference(m.v[:,'c'])]
        group = NMPCVarGroup(varlist, m.time)
        group.set_lb(0, 0)
        group.set_lb(1, 6)
        group.set_lb(2, 0)
        group.set_ub(0, 4)
        group.set_ub(1, 10)
        group.set_ub(2, 10)
        violated = get_violated_bounds_at_time(group, [1,2,3], tolerance=1e-8)
        violated_set = ComponentSet(violated)
        for t in m.time:
            assert m.v[t,'a'] in violated_set
            assert m.v[t,'b'] in violated_set

        violated = get_violated_bounds_at_time(group, 2, tolerance=1e-8)
        violated_set = ComponentSet(violated)
        assert m.v[2,'a'] in violated_set
        assert m.v[2,'b'] in violated_set


if __name__ == '__main__':
    test_find_comp_in_block()
    test_NMPCVarLocator()
    test_copy_values()
    test_find_slices_in_model()
    test_initialize_by_element_in_range()
    test_add_noise_at_time()

