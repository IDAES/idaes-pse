#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################

"""Basic unit tests for PETSc solver utilities"""
import pytest
import pyomo.environ as pyo
from pyomo.dae.flatten import flatten_dae_components
from pyomo.util.subsystems import (
    TemporarySubsystemManager,
    create_subsystem_block,
)
from idaes.core.solvers import petsc
from idaes.core.solvers.features import dae_with_non_time_indexed_constraint

@pytest.mark.unit
def test_copy_time():
    """test the time copy function.  When this is used, the model is flattened
    and only indexed by time, so testing is pretty simple"""

    m = pyo.ConcreteModel()
    t = [1, 2]
    m.x = pyo.Var(t)
    m.y = pyo.Var(t)

    m.x[1] = 1
    m.x[2] = 2
    m.y[1] = 3
    m.y[2] = 4

    petsc._copy_time([m.x, m.y], 1, 2)

    assert pyo.value(m.x[2]) == 1
    assert pyo.value(m.y[2]) == 3

@pytest.mark.unit
def test_gen_time_disc_eqns():
    m, y1, y2, y3, y4, y5, y6 = dae_with_non_time_indexed_constraint()

    # The model has one time element and drivatives for y[1] to y[5]
    # the final time is 180, so create a list of what the time discretization
    # constraints should be.
    disc_eq = [
        id(m.ydot_disc_eq[180, 1]),
        id(m.ydot_disc_eq[180, 2]),
        id(m.ydot_disc_eq[180, 3]),
        id(m.ydot_disc_eq[180, 4]),
        id(m.ydot_disc_eq[180, 5]),
    ]

    n = 0
    for cs in petsc._generate_time_discretization(m, m.t):
        for c in cs.values():
            assert id(c) in disc_eq
            n += 1

    assert len(disc_eq) == n

@pytest.mark.unit
def test_set_dae_suffix():
    m, y1, y2, y3, y4, y5, y6 = dae_with_non_time_indexed_constraint()
    regular_vars, time_vars = flatten_dae_components(m, m.t, pyo.Var)
    regular_cons, time_cons = flatten_dae_components(m, m.t, pyo.Constraint)
    t = 180
    m.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
    m.scaling_factor[m.ydot[180, 2]] = 10
    constraints = [con[t] for con in time_cons if t in con]
    variables = [var[t] for var in time_vars]
    t_block = create_subsystem_block(constraints, variables)
    petsc._set_dae_suffixes_from_variables(t_block, variables)
    petsc._sub_problem_scaling_suffix(m, t_block)

    assert t_block.dae_suffix[m.ydot[180, 1]] == 2
    assert t_block.dae_suffix[m.ydot[180, 2]] == 2
    assert t_block.dae_suffix[m.ydot[180, 3]] == 2
    assert t_block.dae_suffix[m.ydot[180, 4]] == 2
    assert t_block.dae_suffix[m.ydot[180, 5]] == 2
    assert t_block.scaling_factor[m.ydot[180, 2]] == 10
    assert t_block.dae_suffix[m.y[180, 1]] == 1
    assert t_block.dae_suffix[m.y[180, 2]] == 1
    assert t_block.dae_suffix[m.y[180, 3]] == 1
    assert t_block.dae_suffix[m.y[180, 4]] == 1
    assert t_block.dae_suffix[m.y[180, 5]] == 1
