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

from pyomo.environ import SolverFactory, Var, value, Reference
from pyomo.common.collections import ComponentSet, ComponentMap
from pyomo.core.expr.visitor import identify_variables
from pyomo.dae.flatten import flatten_dae_components

from idaes.core.util.model_statistics import (
        degrees_of_freedom, 
        activated_equalities_generator,
        )
from idaes.apps.caprese.util import *
from idaes.apps.caprese.common.config import NoiseBoundOption
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


@pytest.mark.unit
def test_get_violated_bounds():
    bounds = (1., 2.)
    val = 1.5
    assert get_violated_bounds(val, bounds) == (None, 0)

    val = 0.8
    assert get_violated_bounds(val, bounds) == (1., 1)

    val = 2.5
    assert get_violated_bounds(val, bounds) == (2., -1)


@pytest.mark.unit
def test_apply_noise():
    random.seed(1234)
    noise_function = random.gauss
    val_list = [1., 2., 3.]
    params_list = [0.05, 0.05, 0.05]
    new_vals = apply_noise(val_list, params_list, noise_function)

    for val, new_val, p in zip(val_list, new_vals, params_list):
        assert abs(val - new_val) < 3*p

    noise_function = lambda val, rad: random.uniform(val - rad, val + rad)
    new_vals = apply_noise(val_list, params_list, noise_function)
    for val, new_val, p in zip(val_list, new_vals, params_list):
        assert abs(val - new_val) <= p


@pytest.mark.unit
def test_apply_noise_with_bounds():
    NBO = NoiseBoundOption


@pytest.mark.component
def _test_add_noise_at_time():
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
