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
import numpy as np
import json
import os
import pyomo.environ as pyo
import pyomo.dae as pyodae
from pyomo.util.subsystems import (
    TemporarySubsystemManager,
    create_subsystem_block,
)
from idaes.core.solvers import petsc
from idaes.core.solvers.features import dae

def rp_example():
    """This example is done multiple ways to test a few common errors where the
    integrator and fully time-discretized problem differ and alternative working
    formulations.

    The PETSc utilities raise an exception when a differential variable is fixed.
    While this is okay for the fully time-discretized problem, the integrator
    will not correctly link a fixed differential variable (at non-initial time
    points) with a time derivative.
    """
    m = pyo.ConcreteModel()

    m.time = pyodae.ContinuousSet(initialize=(0.0, 10.0))
    m.x = pyo.Var(m.time)
    m.u = pyo.Var(m.time)
    m.dxdt = pyodae.DerivativeVar(m.x, wrt=m.time)

    def diff_eq_rule(m, t):
        return m.dxdt[t] == m.x[t]**2 - m.u[t]
    m.diff_eq = pyo.Constraint(m.time, rule=diff_eq_rule)

    discretizer = pyo.TransformationFactory('dae.finite_difference')
    discretizer.apply_to(m,nfe=1,scheme='BACKWARD')

    for t in m.time:
        m.x[t].fix(2.0*t)

    m.u[0].fix(1.0)

    return m

def rp_example2():
    """ This example is done multiple ways to test a few common errors where the
    integrator and fully time-discretized problem differ, and alternative working
    formulations.

    Rather than fixing a differential variable, we can add an explicit time
    variable for the integrator and fix the time variable to the time index
    for the fully discretized problem. This works as an alternative to the
    fixed differential variable.
    """
    m = pyo.ConcreteModel()

    m.time = pyodae.ContinuousSet(initialize=(0.0, 10.0))
    m.x = pyo.Var(m.time)
    m.u = pyo.Var(m.time)
    m.t = pyo.Var(m.time)
    m.dxdt = pyodae.DerivativeVar(m.x, wrt=m.time)

    def diff_eq_rule(m, t):
        return m.dxdt[t] == m.x[t]**2 - m.u[t]
    m.diff_eq = pyo.Constraint(m.time, rule=diff_eq_rule)

    def x_eq_rule(m, t):
        return m.x[t] == 2.0*m.t[t]
    m.x_eq = pyo.Constraint(m.time, rule=x_eq_rule)

    discretizer = pyo.TransformationFactory('dae.finite_difference')
    discretizer.apply_to(m,nfe=1,scheme='BACKWARD')

    m.u[0].fix(1.0)
    # For fully discretized fix all times at time index
    m.t[0].fix(m.time.first())

    return m

def rp_example3():
    """ This example is done multiple ways to test a few common errors where the
    integrator and fully time-discretized problem differ, and alternative working
    formulations.

    Another way to formulate this problem is to fix the derivative. This doesn't
    work for the integrator because it loses the association between the
    differential variable and its derivative.  Since for users, the result of
    this formulation may be unexpected, the PETSc utilities will raise an
    exception if derivatives are fixed to anything other than 0.
    """
    m = pyo.ConcreteModel()

    m.time = pyodae.ContinuousSet(initialize=(0.0, 10.0))
    m.x = pyo.Var(m.time)
    m.u = pyo.Var(m.time)
    m.dxdt = pyodae.DerivativeVar(m.x, wrt=m.time)

    def diff_eq_rule(m, t):
        return m.dxdt[t] == m.x[t]**2 - m.u[t]
    m.diff_eq = pyo.Constraint(m.time, rule=diff_eq_rule)

    discretizer = pyo.TransformationFactory('dae.finite_difference')
    discretizer.apply_to(m,nfe=1,scheme='BACKWARD')

    m.u[0].fix(1.0)
    m.dxdt[:].fix(2.0)

    return m

def rp_example4():
    """ This example is done multiple ways to test a few common errors where the
    integrator and fully time-discretized problem differ, and alternative working
    formulations.

    Rather than fixing the derivative, we can add a constraint to set the
    derivative.  This should work as intended for both the fully
    time-discretized problem and integrator.
    """
    m = pyo.ConcreteModel()

    m.time = pyodae.ContinuousSet(initialize=(0.0, 10.0))
    m.x = pyo.Var(m.time, initialize=1)
    m.u = pyo.Var(m.time, initialize=1)
    m.dxdt = pyodae.DerivativeVar(m.x, wrt=m.time)

    def diff_eq1_rule(m, t):
        return m.dxdt[t] == m.x[t]**2 - m.u[t]
    m.diff_eq1 = pyo.Constraint(m.time, rule=diff_eq1_rule)

    def diff_eq2_rule(m, t):
        return m.dxdt[t] == 2.0
    m.diff_eq2 = pyo.Constraint(m.time, rule=diff_eq2_rule)

    discretizer = pyo.TransformationFactory('dae.finite_difference')
    discretizer.apply_to(m,nfe=1,scheme='BACKWARD')

    m.u[0].fix(1.0)
    m.x[0].fix(0.0)
    m.diff_eq2[0].deactivate()

    return m

def car_example():
    """This is to test problems where a differential variable doesn't appear in
    a constraint this is based on a Pyomo example here:
    https://github.com/Pyomo/pyomo/blob/main/examples/dae/car_example.py"""
    m = pyo.ConcreteModel()

    m.R = pyo.Param(initialize=0.001) #  Friction factor
    m.L = pyo.Param(initialize=100.0) #  Final position

    m.tau = pyodae.ContinuousSet(bounds=(0,1)) # Unscaled time
    m.time = pyo.Var(m.tau) # Scaled time
    m.tf = pyo.Var()
    m.x = pyo.Var(m.tau,bounds=(0,m.L+50))
    m.v = pyo.Var(m.tau,bounds=(0,None))
    m.a = pyo.Var(m.tau, bounds=(-3.0,1.0),initialize=0)

    m.dtime = pyodae.DerivativeVar(m.time)
    m.dx = pyodae.DerivativeVar(m.x)
    m.dv = pyodae.DerivativeVar(m.v)

    m.obj = pyo.Objective(expr=m.tf)

    def _ode1(m,i):
        if i == 0 :
            return pyo.Constraint.Skip
        return m.dx[i] == m.tf * m.v[i]
    m.ode1 = pyo.Constraint(m.tau, rule=_ode1)

    def _ode2(m,i):
        if i == 0 :
            return pyo.Constraint.Skip
        return m.dv[i] == m.tf*(m.a[i] - m.R*m.v[i]**2)
    m.ode2 = pyo.Constraint(m.tau, rule=_ode2)

    def _ode3(m,i):
        if i == 0:
            return pyo.Constraint.Skip
        return m.dtime[i] == m.tf
    m.ode3 = pyo.Constraint(m.tau, rule=_ode3)

    def _init(m):
        yield m.x[0] == 0
        #yield m.x[1] == m.L
        yield m.v[0] == 0
        yield m.v[1] == 0
        yield m.time[0] == 0
    m.initcon = pyo.ConstraintList(rule=_init)

    discretizer = pyo.TransformationFactory('dae.finite_difference')
    discretizer.apply_to(m,nfe=1,scheme='BACKWARD')
    return m


def dae_with_non_time_indexed_constraint():
    """This provides a DAE model for solver testing. This model contains a non-
    time-indexed variable and constraint and a fixed derivative to test some
    edge cases.

    The problem and expected result are from A test problem from
    https://archimede.dm.uniba.it/~testset/report/chemakzo.pdf.

    Args:
        None

    Returns:
        (tuple): Pyomo ConcreteModel, correct solved value for y[1] to y[6]
    """
    model = pyo.ConcreteModel(name="chemakzo")

    # Set problem parameter values
    model.k = pyo.Param([1,2,3,4], initialize={
        1:18.7,
        2:0.58,
        3:0.09,
        4:0.42})
    model.Ke = pyo.Param(initialize=34.4)
    model.klA = pyo.Param(initialize=3.3)
    model.Ks = pyo.Param(initialize=115.83)
    model.pCO2 = pyo.Param(initialize=0.9)
    # The following parameter H, is best made a parameter, but will use a
    # variable and constraint instead to test non-time-indexed constraints
    #model.H = pyo.Param(initialize=737)

    # Problem variables ydot = dy/dt,
    #    (dy6/dt is not explicitly in the equations, so only 5 ydots)
    model.H = pyo.Var(initialize=737)
    model.t = pyodae.ContinuousSet(bounds=(0,180))
    model.y = pyo.Var(model.t, [1,2,3,4,5,6], initialize=1.0)  #
    model.ydot = pyodae.DerivativeVar(model.y, wrt=model.t) # dy/dt
    model.r = pyo.Var(model.t, [1,2,3,4,5], initialize=1.0)
    model.Fin = pyo.Var(model.t, initialize=1.0)

    # Non-time indexed constraint (just for testing)
    model.H_eqn = pyo.Constraint(expr=model.H==737)

    # Equations
    @model.Constraint(model.t)
    def eq_ydot1(b, t):
        return b.ydot[t, 1] == -2.0*b.r[t, 1] + b.r[t, 2] - b.r[t, 3] - b.r[t, 4]
    @model.Constraint(model.t)
    def eq_ydot2(b, t):
        return b.ydot[t, 2] == -0.5*b.r[t, 1] - b.r[t, 4] - 0.5*b.r[t, 5] + b.Fin[t]
    @model.Constraint(model.t)
    def eq_ydot3(b, t):
        return b.ydot[t, 3] == b.r[t, 1] - b.r[t, 2] + b.r[t, 3]
    @model.Constraint(model.t)
    def eq_ydot4(b, t):
        return b.ydot[t, 4] == -b.r[t, 2] + b.r[t, 3] - 2.0*b.r[t, 4]
    @model.Constraint(model.t)
    def eq_ydot5(b, t):
        return b.ydot[t, 5] == b.r[t, 2] - b.r[t, 3] + b.r[t, 5]
    @model.Constraint(model.t)
    def eq_y6(b, t):
        return b.ydot[t, 6] == b.Ks*b.y[t, 1]*b.y[t, 4] - b.y[t, 6]

    model.ydot[:, 6].fix(0)

    @model.Constraint(model.t)
    def eq_r1(b, t):
        return b.r[t, 1] == b.k[1]*b.y[t, 1]**4*b.y[t, 2]**0.5
    @model.Constraint(model.t)
    def eq_r2(b, t):
        return b.r[t, 2] == b.k[2]*b.y[t, 3]*b.y[t, 4]
    @model.Constraint(model.t)
    def eq_r3(b, t):
        return b.r[t, 3] == b.k[2]/b.Ke*b.y[t, 1]*b.y[t, 5]
    @model.Constraint(model.t)
    def eq_r4(b, t):
        return b.r[t, 4] == b.k[3]*b.y[t, 1]*b.y[t, 4]**2
    @model.Constraint(model.t)
    def eq_r5(b, t):
        return b.r[t, 5] == b.k[4]*b.y[t, 6]**2*b.y[t, 2]**0.5
    @model.Constraint(model.t)
    def eq_Fin(b, t):
        return b.Fin[t] == b.klA*(b.pCO2/b.H - b.y[t, 2])

    # Set initial conditions and solve initial from the values of differential
    # variables (r and y6 well and the derivative vars too).
    y0 = {1:0.444, 2:0.00123, 3:0.0, 4:0.007, 5:0.0} #initial differential vars
    for i in y0:
        model.y[0, i].fix(y0[i])

    model.eq_ydot1[0].deactivate()
    model.eq_ydot2[0].deactivate()
    model.eq_ydot3[0].deactivate()
    model.eq_ydot4[0].deactivate()
    model.eq_ydot5[0].deactivate()

    discretizer = pyo.TransformationFactory('dae.finite_difference')
    discretizer.apply_to(model, nfe=1, scheme='BACKWARD')

    return (
        model,
        0.1150794920661702,
        0.1203831471567715e-2,
        0.1611562887407974,
        0.3656156421249283e-3,
        0.1708010885264404e-1,
        0.4873531310307455e-2,
    )


@pytest.mark.unit
@pytest.mark.skipif(not petsc.petsc_available(), reason="PETSc solver not available")
def test_car():
    m = car_example()
    m.a.fix(1.0)
    m.tf.fix(16.56)

    # solve
    petsc.petsc_dae_by_time_element(
        m,
        time=m.tau,
        ts_options={
            "--ts_type":"cn", # Crank–Nicolson
            "--ts_adapt_type":"basic",
            "--ts_dt":0.01,
        },
    )

    assert pyo.value(m.x[1]) == pytest.approx(131.273, rel=1e-2)


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

    # The model has one time element and derivatives for y[1] to y[5]
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
    for cs in petsc.find_discretization_equations(m, m.t):
        for c in cs.values():
            if c.index()[1] == 6:
                continue
            assert id(c) in disc_eq
            n += 1

    assert len(disc_eq) == n

@pytest.mark.unit
def test_set_dae_suffix():
    m, y1, y2, y3, y4, y5, y6 = dae_with_non_time_indexed_constraint()
    regular_vars, time_vars = pyodae.flatten.flatten_dae_components(m, m.t, pyo.Var)
    regular_cons, time_cons = pyodae.flatten.flatten_dae_components(m, m.t, pyo.Constraint)
    t = 180
    m.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
    m.scaling_factor[m.ydot[180, 2]] = 10
    constraints = [con[t] for con in time_cons if t in con]
    variables = [var[t] for var in time_vars]
    t_block = create_subsystem_block(constraints, variables)
    deriv_diff_map =  petsc._get_derivative_differential_data_map(m, m.t)
    petsc._set_dae_suffixes_from_variables(t_block, variables, deriv_diff_map)
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

    # Make sure deactivating a differential equation makes the variable that
    # would have been differential go algebraic
    m, y1, y2, y3, y4, y5, y6 = dae_with_non_time_indexed_constraint()
    # discretization equations would be deactivated in normal PETSc solve
    for con in petsc.find_discretization_equations(m, m.t):
        con.deactivate()
    # deactivate a differential equation making y4 be algebraic
    m.eq_ydot4[180].deactivate()
    regular_vars, time_vars = pyodae.flatten.flatten_dae_components(m, m.t, pyo.Var)
    regular_cons, time_cons = pyodae.flatten.flatten_dae_components(m, m.t, pyo.Constraint)
    t = 180
    constraints = [con[t] for con in time_cons if t in con]
    variables = [var[t] for var in time_vars]
    t_block = create_subsystem_block(constraints, variables)
    deriv_diff_map = petsc._get_derivative_differential_data_map(m, m.t)
    petsc._set_dae_suffixes_from_variables(t_block, variables, deriv_diff_map)
    petsc._sub_problem_scaling_suffix(m, t_block)

    assert m.y[t, 4] not in t_block.dae_suffix
    assert m.ydot[t, 4] not in t_block.dae_suffix

@pytest.mark.unit
@pytest.mark.skipif(not petsc.petsc_available(), reason="PETSc solver not available")
def test_petsc_read_trajectory():
    """
    Check that the PETSc DAE solver works.
    """
    m, y1, y2, y3, y4, y5, y6 = dae_with_non_time_indexed_constraint()
    m.scaling_factor = pyo.Suffix(direction=pyo.Suffix.EXPORT)
    m.scaling_factor[m.y[180, 1]] = 10 # make sure unscale works
    
    m.y_ref = pyo.Reference(m.y) # make sure references don't get unscaled twice
    petsc.petsc_dae_by_time_element(
        m,
        time=m.t,
        ts_options={
            "--ts_type":"cn", # Crank–Nicolson
            "--ts_adapt_type":"basic",
            "--ts_dt":0.01,
            "--ts_save_trajectory":1,
            "--ts_trajectory_type":"visualization",
        },
        vars_stub="tj_random_junk_123",
        trajectory_save_prefix="tj_random_junk_123",
    )
    assert pytest.approx(y1, rel=1e-3) == pyo.value(m.y[m.t.last(), 1])
    assert pytest.approx(y2, rel=1e-3) == pyo.value(m.y[m.t.last(), 2])
    assert pytest.approx(y3, rel=1e-3) == pyo.value(m.y[m.t.last(), 3])
    assert pytest.approx(y4, rel=1e-3) == pyo.value(m.y[m.t.last(), 4])
    assert pytest.approx(y5, rel=1e-3) == pyo.value(m.y[m.t.last(), 5])
    assert pytest.approx(y6, rel=1e-3) == pyo.value(m.y[m.t.last(), 6])

    tj = petsc.PetscTrajectory(json="tj_random_junk_123_1.json.gz")
    assert tj.get_dt()[0] == pytest.approx(0.01) # if small enough shouldn't be cut
    assert tj.get_vec(m.y[180, 1])[-1] == pytest.approx(y1, rel=1e-3)
    assert tj.get_vec("_time")[-1] == pytest.approx(180)

    times = np.linspace(0, 180, 181)
    vecs = tj.interpolate_vecs(times)
    assert vecs[str(m.y[180, 1])][180] == pytest.approx(y1, rel=1e-3)
    assert vecs["_time"][180] == pytest.approx(180)

    tj.to_json("some_testy_json.json")
    with open("some_testy_json.json", "r") as fp:
        vecs = json.load(fp)
    assert vecs[str(m.y[180, 1])][-1] == pytest.approx(y1, rel=1e-3)
    assert vecs["_time"][-1] == pytest.approx(180)
    os.remove("some_testy_json.json")

    tj.to_json("some_testy_json.json.gz")
    tj2 = petsc.PetscTrajectory(json="some_testy_json.json.gz")
    assert tj2.vecs[str(m.y[180, 1])][-1] == pytest.approx(y1, rel=1e-3)
    assert tj2.vecs["_time"][-1] == pytest.approx(180)
    os.remove("some_testy_json.json.gz")

    tj2 = petsc.PetscTrajectory(vecs=vecs)
    assert tj2.vecs[str(m.y[180, 1])][-1] == pytest.approx(y1, rel=1e-3)
    assert tj2.vecs["_time"][-1] == pytest.approx(180)


@pytest.mark.unit
@pytest.mark.skipif(not petsc.petsc_available(), reason="PETSc solver not available")
def test_rp_example():

    m = rp_example()
    with pytest.raises(RuntimeError):
        petsc.petsc_dae_by_time_element(
            m,
            time=m.time,
        )


@pytest.mark.unit
@pytest.mark.skipif(not petsc.petsc_available(), reason="PETSc solver not available")
def test_rp_example2():

    m = rp_example2()
    petsc.petsc_dae_by_time_element(
        m,
        time=m.time,
        timevar=m.t,
        ts_options={
            "--ts_dt":1,
            "--ts_adapt_type":"none",
        }
    )
    assert pyo.value(m.u[10]) == pytest.approx(398)
    assert pyo.value(m.x[10]) == pytest.approx(20)


@pytest.mark.unit
@pytest.mark.skipif(not petsc.petsc_available(), reason="PETSc solver not available")
def test_rp_example3():

    m = rp_example3()
    with pytest.raises(RuntimeError):
        petsc.petsc_dae_by_time_element(
            m,
            time=m.time,
        )

@pytest.mark.unit
@pytest.mark.skipif(not petsc.petsc_available(), reason="PETSc solver not available")
def test_rp_example4():

    m = rp_example4()
    petsc.petsc_dae_by_time_element(
        m,
        time=m.time,
        ts_options={
            "--ts_dt":1,
            "--ts_adapt_type":"none",
        }
    )
    assert pyo.value(m.u[10]) == pytest.approx(398)
    assert pyo.value(m.x[10]) == pytest.approx(20)
