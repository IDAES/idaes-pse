#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

"""Basic unit tests for PETSc solver utilities"""

import pytest
import re
import numpy as np
import json
import os
from pyomo.common.collections import ComponentSet
from pyomo.environ import (
    assert_optimal_termination,
    Block,
    ConcreteModel,
    Constraint,
    ConstraintList,
    Objective,
    Param,
    Reference,
    Set,
    Suffix,
    TransformationFactory,
    value,
    Var,
)
from pyomo.dae import ContinuousSet, DerivativeVar
from idaes.core.scaling.util import set_scaling_factor
from idaes.core.solvers.petsc_object import (
    PETScIntegrator,
    get_derivative_differential_vardata_map,
)
from idaes.core.solvers.petsc import petsc_available, PetscTrajectory
from idaes.core.util import DiagnosticsToolbox, from_json, StoreSpec


def car_example():
    """This is to test problems where a differential variable doesn't appear in
    a constraint this is based on a Pyomo example here:
    https://github.com/Pyomo/pyomo/blob/main/examples/dae/car_example.py"""
    m = ConcreteModel()

    m.R = Param(initialize=0.001)  #  Friction factor
    m.L = Param(initialize=100.0)  #  Final position

    m.tau = ContinuousSet(bounds=(0, 1))  # Unscaled time
    m.time = Var(m.tau)  # Scaled time
    m.tf = Var()
    m.x = Var(m.tau, bounds=(0, m.L + 50))
    m.v = Var(m.tau, bounds=(0, None))
    m.a = Var(m.tau, bounds=(-3.0, 1.0), initialize=0)

    m.dtime = DerivativeVar(m.time)
    m.dx = DerivativeVar(m.x)
    m.dv = DerivativeVar(m.v)

    m.obj = Objective(expr=m.tf)

    def _ode1(m, i):
        if i == 0:
            return Constraint.Skip
        return m.dx[i] == m.tf * m.v[i]

    m.ode1 = Constraint(m.tau, rule=_ode1)

    def _ode2(m, i):
        if i == 0:
            return Constraint.Skip
        return m.dv[i] == m.tf * (m.a[i] - m.R * m.v[i] ** 2)

    m.ode2 = Constraint(m.tau, rule=_ode2)

    def _ode3(m, i):
        if i == 0:
            return Constraint.Skip
        return m.dtime[i] == m.tf

    m.ode3 = Constraint(m.tau, rule=_ode3)

    def _init(m):
        yield m.x[0] == 0
        # yield m.x[1] == m.L
        yield m.v[0] == 0
        yield m.v[1] == 0
        yield m.time[0] == 0

    m.initcon = ConstraintList(rule=_init)

    discretizer = TransformationFactory("dae.finite_difference")
    discretizer.apply_to(m, nfe=1, scheme="BACKWARD")
    return m


def dae_with_non_time_indexed_constraint(
    nfe=1, ncp=3, transformation_method="dae.finite_difference", scheme="BACKWARD"
):
    """This function provides a DAE model for solver testing. This model contains a non-
    time-indexed variable and constraint and a fixed derivative to test some
    edge cases.

    The problem and expected result are from A test problem from
    https://archimede.dm.uniba.it/~testset/report/chemakzo.pdf.

    Args:
        nfe: Number of finite elements to use in discretization
        ncp: Number of collocation points to use, if using collocation
        transformation_method: Discretization method. Presently,
            options are "dae.finite_difference" and "dae.collocation".
        scheme: Scheme to use in discretization method. Check Pyomo
            DAE documentation for more info

    Returns:
        (tuple): Pyomo ConcreteModel, correct solved value for y[1] to y[6]
    """
    model = ConcreteModel(name="chemakzo")

    # Set problem parameter values
    model.k = Param([1, 2, 3, 4], initialize={1: 18.7, 2: 0.58, 3: 0.09, 4: 0.42})
    model.Ke = Param(initialize=34.4)
    model.klA = Param(initialize=3.3)
    model.Ks = Param(initialize=115.83)
    model.pCO2 = Param(initialize=0.9)
    # The following parameter H, is best made a parameter, but will use a
    # variable and constraint instead to test non-time-indexed constraints
    # model.H = Param(initialize=737)

    # Problem variables ydot = dy/dt,
    #    (dy6/dt is not explicitly in the equations, so only 5 ydots)
    model.H = Var(initialize=737)
    model.t = ContinuousSet(bounds=(0, 180))
    model.y = Var(model.t, [1, 2, 3, 4, 5, 6], initialize=1.0)  #
    model.ydot = DerivativeVar(model.y, wrt=model.t)  # dy/dt
    model.r = Var(model.t, [1, 2, 3, 4, 5], initialize=1.0)
    model.Fin = Var(model.t, initialize=1.0)

    # Non-time indexed constraint (just for testing)
    model.H_eqn = Constraint(expr=model.H == 737)

    # Equations
    @model.Constraint(model.t)
    def eq_ydot1(b, t):
        return b.ydot[t, 1] == -2.0 * b.r[t, 1] + b.r[t, 2] - b.r[t, 3] - b.r[t, 4]

    @model.Constraint(model.t)
    def eq_ydot2(b, t):
        return b.ydot[t, 2] == -0.5 * b.r[t, 1] - b.r[t, 4] - 0.5 * b.r[t, 5] + b.Fin[t]

    @model.Constraint(model.t)
    def eq_ydot3(b, t):
        return b.ydot[t, 3] == b.r[t, 1] - b.r[t, 2] + b.r[t, 3]

    @model.Constraint(model.t)
    def eq_ydot4(b, t):
        return b.ydot[t, 4] == -b.r[t, 2] + b.r[t, 3] - 2.0 * b.r[t, 4]

    @model.Constraint(model.t)
    def eq_ydot5(b, t):
        return b.ydot[t, 5] == b.r[t, 2] - b.r[t, 3] + b.r[t, 5]

    @model.Constraint(model.t)
    def eq_y6(b, t):
        return 0 == b.Ks * b.y[t, 1] * b.y[t, 4] - b.y[t, 6]

    @model.Constraint(model.t)
    def eq_r1(b, t):
        return b.r[t, 1] == b.k[1] * b.y[t, 1] ** 4 * b.y[t, 2] ** 0.5

    @model.Constraint(model.t)
    def eq_r2(b, t):
        return b.r[t, 2] == b.k[2] * b.y[t, 3] * b.y[t, 4]

    @model.Constraint(model.t)
    def eq_r3(b, t):
        return b.r[t, 3] == b.k[2] / b.Ke * b.y[t, 1] * b.y[t, 5]

    @model.Constraint(model.t)
    def eq_r4(b, t):
        return b.r[t, 4] == b.k[3] * b.y[t, 1] * b.y[t, 4] ** 2

    @model.Constraint(model.t)
    def eq_r5(b, t):
        return b.r[t, 5] == b.k[4] * b.y[t, 6] ** 2 * b.y[t, 2] ** 0.5

    @model.Constraint(model.t)
    def eq_Fin(b, t):
        return b.Fin[t] == b.klA * (b.pCO2 / b.H - b.y[t, 2])

    # Set initial conditions and solve initial from the values of differential
    # variables (r and y6 well and the derivative vars too).
    y0 = {1: 0.444, 2: 0.00123, 3: 0.0, 4: 0.007, 5: 0.0}  # initial differential vars
    for i in y0:
        model.y[0, i].fix(y0[i])

    if scheme == "BACKWARD":
        model.eq_ydot1[0].deactivate()
        model.eq_ydot2[0].deactivate()
        model.eq_ydot3[0].deactivate()
        model.eq_ydot4[0].deactivate()
        model.eq_ydot5[0].deactivate()

    discretizer = TransformationFactory(transformation_method)
    if transformation_method == "dae.collocation":
        discretizer.apply_to(model, nfe=nfe, ncp=ncp, scheme=scheme)
    else:
        discretizer.apply_to(model, nfe=nfe, scheme=scheme)
    # y6 has been converted into an algebraic variable. Therefore, we need to
    # deactivate the corresponding discretization equations.
    model.ydot[:, 6].fix(0)
    model.ydot_disc_eq[:, 6].deactivate()

    return (
        model,
        0.1150794920661702,
        0.1203831471567715e-2,
        0.1611562887407974,
        0.3656156421249283e-3,
        0.1708010885264404e-1,
        0.4873531310307455e-2,
    )


class TestGetDerivativeDifferentialVardataMap(object):
    @pytest.mark.unit
    def test_get_derivative_differential_vardata_map_basic_finite_difference(self):
        # Tests basic mapping with a single-index ContinuousSet using Finite Difference
        m = ConcreteModel()
        m.t = ContinuousSet(bounds=(0, 10))
        m.x = Var(m.t)
        m.dxdt = DerivativeVar(m.x, wrt=m.t)

        TransformationFactory("dae.finite_difference").apply_to(
            m, nfe=2, scheme="BACKWARD"
        )

        deriv_map = get_derivative_differential_vardata_map(m, m.t)

        # The derivative at zero won't be mapped because it
        # doesn't have an associated discretization equation
        assert len(deriv_map) == len(m.t) - 1
        for idx in m.t:
            deriv_data = m.dxdt[idx]
            if idx == 0:
                assert deriv_data not in deriv_map
            else:
                assert deriv_data in deriv_map
                assert deriv_map[deriv_data]["diff_var"] is m.x[idx]
                # Finite difference generates discretization equations
                assert deriv_map[deriv_data]["disc_eq"] is m.dxdt_disc_eq[idx]
                assert "cont_eq" not in deriv_map[deriv_data]

    @pytest.mark.unit
    def test_get_derivative_differential_vardata_map_multidimensional_indexing(self):
        # Tests multidimensional indexing:
        # 1. Non-continuous set 'i' is the first index, ContinuousSet 't' is second
        # 2. ContinuousSet 't' is first, non-continuous set 'i' is second
        m = ConcreteModel()
        m.t = ContinuousSet(bounds=(0, 10))
        m.i = Set(initialize=[1, 2])

        m.x = Var(m.i, m.t)
        m.dxdt = DerivativeVar(m.x, wrt=m.t)

        m.y = Var(m.t, m.i)
        m.dydt = DerivativeVar(m.y, wrt=m.t)

        TransformationFactory("dae.finite_difference").apply_to(
            m, nfe=1, scheme="BACKWARD"
        )

        # 1. Map for x (i first, t second)
        deriv_map_x = get_derivative_differential_vardata_map(m, m.t)
        for idx in m.i:
            for t_val in m.t:
                deriv_data = m.dxdt[idx, t_val]
                if t_val == 0:
                    # No discretization equation at t=0 for backward difference
                    assert deriv_data not in deriv_map_x
                else:
                    assert deriv_data in deriv_map_x
                    assert deriv_map_x[deriv_data]["diff_var"] is m.x[idx, t_val]
                    assert (
                        deriv_map_x[deriv_data]["disc_eq"] is m.dxdt_disc_eq[idx, t_val]
                    )

        # 2. Map for y (t first, i second)
        deriv_map_y = get_derivative_differential_vardata_map(m, m.t)
        for t_val in m.t:
            for idx in m.i:
                deriv_data = m.dydt[t_val, idx]
                if t_val == 0:
                    # No discretization equation at t=0 for backward difference
                    assert deriv_data not in deriv_map_y
                else:
                    assert deriv_data in deriv_map_y
                    assert deriv_map_y[deriv_data]["diff_var"] is m.y[t_val, idx]
                    assert (
                        deriv_map_y[deriv_data]["disc_eq"] is m.dydt_disc_eq[t_val, idx]
                    )

    @pytest.mark.unit
    def test_get_derivative_differential_vardata_map_deactivation(self):
        # Tests exclusion when constraint components or specific index data are deactivated
        m = ConcreteModel()
        m.t = ContinuousSet(bounds=(0, 10))
        m.x = Var(m.t)
        m.dxdt = DerivativeVar(m.x, wrt=m.t)

        m.y = Var(m.t)
        m.dydt = DerivativeVar(m.y, wrt=m.t)

        TransformationFactory("dae.finite_difference").apply_to(
            m, nfe=2, scheme="BACKWARD"
        )

        # Deactivate the entire discretization constraint component for x
        m.dxdt_disc_eq.deactivate()

        # Deactivate a single constraint data index for y (at t = 5.0)
        m.dydt_disc_eq[5.0].deactivate()

        deriv_map = get_derivative_differential_vardata_map(m, m.t)

        # x should have no active mappings because the entire constraint component is inactive
        for idx in m.t:
            assert m.dxdt[idx] not in deriv_map

        # y should have an active mapping only at t = 10
        for idx in m.t:
            deriv_data = m.dydt[idx]
            if idx == pytest.approx(10):
                assert deriv_data in deriv_map
                assert deriv_map[deriv_data]["diff_var"] is m.y[idx]
                assert deriv_map[deriv_data]["disc_eq"] is m.dydt_disc_eq[idx]
            else:
                assert deriv_data not in deriv_map

    @pytest.mark.unit
    def test_get_derivative_differential_vardata_map_collocation_lagrange_legendre(
        self,
    ):
        # Tests Lagrange-Legendre collocation, which generates both
        # discretization equations and continuity equations.
        m = ConcreteModel()
        m.t = ContinuousSet(bounds=(0, 10))
        m.x = Var(m.t)
        # For extra fun, put the DerivativeVar on a
        # different block to ensure that both the
        # continuity equation and discretization equation
        # are created on the same block as the DerivativeVar
        m.blk = Block()
        m.blk.dxdt = DerivativeVar(m.x, wrt=m.t)

        TransformationFactory("dae.collocation").apply_to(
            m, nfe=2, ncp=2, scheme="LAGRANGE-LEGENDRE"
        )

        deriv_map = get_derivative_differential_vardata_map(m, m.t)

        # In Lagrange-Legendre collocation:
        # - Collocation points have discretization equations (disc_eq).
        # - Finite element boundaries (except the first point) have continuity equations (cont_eq).
        for idx in m.t:
            deriv_data = m.blk.dxdt[idx]

            has_disc = hasattr(m.blk, "dxdt_disc_eq") and idx in m.blk.dxdt_disc_eq
            has_cont = hasattr(m.blk, "x_t_cont_eq") and idx in m.blk.x_t_cont_eq

            if has_disc or has_cont:
                assert deriv_data in deriv_map
                assert deriv_map[deriv_data]["diff_var"] is m.x[idx]

                if has_disc:
                    assert deriv_map[deriv_data]["disc_eq"] is m.blk.dxdt_disc_eq[idx]
                if has_cont:
                    assert deriv_map[deriv_data]["cont_eq"] is m.blk.x_t_cont_eq[idx]
            else:
                # Point has neither constraint, so it should not be mapped
                assert deriv_data not in deriv_map


@pytest.mark.unit
@pytest.mark.skipif(not petsc_available(), reason="PETSc solver not available")
def test_car():
    m = car_example()
    m.a.fix(1.0)
    m.tf.fix(16.56)

    # solve
    petsc_obj = PETScIntegrator(
        m,
        time_set=m.tau,
        ts_options={
            "--ts_type": "cn",  # Crank–Nicolson
            "--ts_adapt_type": "basic",
            "--ts_dt": 0.01,
        },
        calculate_derivatives=True,
    )
    petsc_obj.dae_by_time_element()

    assert value(m.x[1]) == pytest.approx(131.273, rel=1e-2)
    assert value(m.dtime[1]) == pytest.approx(value(m.time[1] - m.time[0]))
    assert value(m.dx[1]) == pytest.approx(value(m.x[1] - m.x[0]))
    assert value(m.dv[1]) == pytest.approx(value(m.v[1] - m.v[0]))


@pytest.mark.unit
def test_copy_time():
    """test the time copy function.  When this is used, the model is flattened
    and only indexed by time, so testing is pretty simple"""

    m = ConcreteModel()
    t = [1, 2]
    m.x = Var(t)
    m.y = Var(t)

    m.x[1] = 1
    m.x[2] = 2
    m.y[1] = 3
    m.y[2] = 4

    class dummy_integrator:
        time_variables = [m.x, m.y]

    PETScIntegrator._copy_time(dummy_integrator(), 1, 2)

    assert value(m.x[2]) == 1
    assert value(m.y[2]) == 3


@pytest.mark.unit
def test_petsc_initial_condition_problem():
    m, y1, y2, y3, y4, y5, y6 = dae_with_non_time_indexed_constraint()
    set_scaling_factor(m.y[0, 3], 7)  # Make sure that scaling factors are copied over
    set_scaling_factor(m.eq_Fin[0], 11)

    petsc_obj = PETScIntegrator(m, time_set=m.t)
    ic_block = petsc_obj.get_initial_condition_problem(0)
    assert ic_block.scaling_factor[m.y[0, 3]] == 7
    assert ic_block.scaling_factor[m.eq_Fin[0]] == 11

    expected_var_set = ComponentSet([m.H, m.Fin[0]])
    expected_con_set = ComponentSet([m.H_eqn, m.eq_Fin[0], m.eq_y6[0]])
    for i in range(1, 7):
        expected_var_set.add(m.y[0, i])
        # ydot[6] comes along for the ride despite being in zero active constraints
        expected_var_set.add(m.ydot[0, i])
        if not i == 6:
            expected_var_set.add(m.r[0, i])
            expected_con_set.add(m.find_component(f"eq_ydot{i}")[0])
            expected_con_set.add(m.find_component(f"eq_r{i}")[0])

    actual_var_set = ComponentSet(
        vardata for vardata in ic_block.component_data_objects(ctype=Var)
    )
    actual_con_set = ComponentSet(
        condata for condata in ic_block.component_data_objects(ctype=Constraint)
    )
    assert expected_var_set == actual_var_set
    assert expected_con_set == actual_con_set


@pytest.mark.unit
def test_petsc_timestep_problem():
    m, y1, y2, y3, y4, y5, y6 = dae_with_non_time_indexed_constraint()
    set_scaling_factor(m.y[180, 3], 7)  # Make sure that scaling factors are copied over
    set_scaling_factor(m.eq_Fin[180], 11)

    petsc_obj = PETScIntegrator(m, time_set=m.t)
    t_block, diff_vars, fixedness_json = petsc_obj.get_timestep_problem(180)
    assert t_block.scaling_factor[m.y[180, 3]] == 7
    assert t_block.scaling_factor[m.eq_Fin[180]] == 11

    expected_var_set = ComponentSet([m.H, m.Fin[180]])
    expected_con_set = ComponentSet([m.eq_Fin[180], m.eq_y6[180]])
    expected_diff_var_set = ComponentSet()
    for i in range(1, 7):
        expected_var_set.add(m.y[180, i])
        # ydot[6] comes along for the ride despite being in zero active constraints
        expected_var_set.add(m.ydot[180, i])
        if not i == 6:
            expected_diff_var_set.add(m.y[180, i])
            expected_var_set.add(m.r[180, i])
            expected_con_set.add(m.find_component(f"eq_ydot{i}")[180])
            expected_con_set.add(m.find_component(f"eq_r{i}")[180])

    actual_var_set = ComponentSet(
        vardata for vardata in t_block.component_data_objects(ctype=Var)
    )
    actual_con_set = ComponentSet(
        condata for condata in t_block.component_data_objects(ctype=Constraint)
    )
    assert expected_var_set == actual_var_set
    assert expected_con_set == actual_con_set
    assert len(diff_vars) == len(
        ComponentSet(diff_vars)
    )  # Don't want redundant diff vars.
    assert expected_diff_var_set == ComponentSet(diff_vars)

    assert m.H.fixed
    input_var_set = ComponentSet(vardata for vardata in t_block.input_vars.values())
    assert m.H in input_var_set
    assert len(input_var_set) == 1
    from_json(m, sd=fixedness_json, wts=StoreSpec.isfixed())


@pytest.mark.unit
@pytest.mark.skipif(not petsc_available(), reason="PETSc solver not available")
def test_petsc_read_trajectory():
    """
    Check that the PETSc DAE solver works.
    """
    m, y1, y2, y3, y4, y5, y6 = dae_with_non_time_indexed_constraint()
    m.scaling_factor = Suffix(direction=Suffix.EXPORT)
    m.scaling_factor[m.y[180, 1]] = 10  # make sure unscale works

    m.y_ref = Reference(m.y)  # make sure references don't get unscaled twice
    petsc_obj = PETScIntegrator(
        m,
        time_set=m.t,
        ts_options={
            "--ts_type": "cn",  # Crank–Nicolson
            "--ts_adapt_type": "basic",
            "--ts_dt": 0.01,
            "--ts_save_trajectory": 1,
            "--ts_trajectory_type": "visualization",
        },
    )
    res = petsc_obj.dae_by_time_element()
    assert pytest.approx(y1, rel=1e-3) == value(m.y[m.t.last(), 1])
    assert pytest.approx(y2, rel=1e-3) == value(m.y[m.t.last(), 2])
    assert pytest.approx(y3, rel=1e-3) == value(m.y[m.t.last(), 3])
    assert pytest.approx(y4, rel=1e-3) == value(m.y[m.t.last(), 4])
    assert pytest.approx(y5, rel=1e-3) == value(m.y[m.t.last(), 5])
    assert pytest.approx(y6, rel=1e-3) == value(m.y[m.t.last(), 6])

    tj = res.trajectory
    assert tj.get_dt()[0] == pytest.approx(0.01)  # if small enough shouldn't be cut
    assert tj.get_vec(m.y[180, 1])[-1] == pytest.approx(y1, rel=1e-3)
    assert tj.get_vec("_time")[-1] == pytest.approx(180)

    times = np.linspace(0, 180, 181)
    tj2 = tj.interpolate(times)
    assert tj2.get_vec(m.y[180, 1])[180] == pytest.approx(y1, rel=1e-3)
    assert tj2.time[180] == pytest.approx(180)

    tj.to_json("some_testy_json.json")
    with open("some_testy_json.json", "r") as fp:
        vecs = json.load(fp)
    assert vecs[str(m.y[180, 1])][-1] == pytest.approx(y1, rel=1e-3)
    assert vecs["_time"][-1] == pytest.approx(180)
    os.remove("some_testy_json.json")

    tj.to_json("some_testy_json.json.gz")
    tj2 = PetscTrajectory(json="some_testy_json.json.gz")
    assert tj2.vecs[str(m.y[180, 1])][-1] == pytest.approx(y1, rel=1e-3)
    assert tj2.vecs["_time"][-1] == pytest.approx(180)
    os.remove("some_testy_json.json.gz")

    tj2 = PetscTrajectory(vecs=vecs)
    assert tj2.vecs[str(m.y[180, 1])][-1] == pytest.approx(y1, rel=1e-3)
    assert tj2.vecs["_time"][-1] == pytest.approx(180)


@pytest.mark.unit
@pytest.mark.skipif(not petsc_available(), reason="PETSc solver not available")
def test_petsc_read_trajectory_parts():
    """
    Check that the PETSc DAE solver works.
    """
    m, y1, y2, y3, y4, y5, y6 = dae_with_non_time_indexed_constraint(nfe=10)
    m.scaling_factor = Suffix(direction=Suffix.EXPORT)
    m.scaling_factor[m.y[180, 1]] = 10  # make sure unscale works

    m.y_ref = Reference(m.y)  # make sure references don't get unscaled twice

    petsc_obj = PETScIntegrator(
        model=m,
        time_set=m.t,
        ts_options={
            "--ts_type": "cn",  # Crank–Nicolson
            "--ts_adapt_type": "basic",
            "--ts_dt": 0.01,
            "--ts_save_trajectory": 1,
        },
    )

    res = petsc_obj.dae_by_time_element(between=[m.t.first(), m.t.at(4), m.t.last()])
    assert pytest.approx(y1, rel=1e-3) == value(m.y[m.t.last(), 1])
    assert pytest.approx(y2, rel=1e-3) == value(m.y[m.t.last(), 2])
    assert pytest.approx(y3, rel=1e-3) == value(m.y[m.t.last(), 3])
    assert pytest.approx(y4, rel=1e-3) == value(m.y[m.t.last(), 4])
    assert pytest.approx(y5, rel=1e-3) == value(m.y[m.t.last(), 5])
    assert pytest.approx(y6, rel=1e-3) == value(m.y[m.t.last(), 6])

    tj = res.trajectory
    assert tj.get_vec(m.y[180, 1])[-1] == pytest.approx(y1, rel=1e-3)
    assert tj.get_vec("_time")[-1] == pytest.approx(180)
    y1_trj = tj.interpolate_vec(m.t, m.y[180, 1])
    y4_trj = tj.interpolate_vec(m.t, m.y[180, 4])
    for i, t in enumerate(m.t):
        assert y1_trj[i] == pytest.approx(value(m.y[t, 1]))
        assert y4_trj[i] == pytest.approx(value(m.y[t, 4]))


@pytest.mark.skipif(not petsc_available(), reason="PETSc solver not available")
class TestRPExamples:
    """This example is done multiple ways to test a few common errors where the
    integrator and fully time-discretized problem differ, and alternative working
    formulations.
    """

    def make_rp_example(self, time_set=None, nfe=None):
        if time_set is None:
            time_set = [0.0, 10.0]
        if nfe is None:
            nfe = len(time_set) - 1
        m = ConcreteModel()

        m.time = ContinuousSet(initialize=time_set)
        m.x = Var(m.time, initialize=0)
        m.u = Var(m.time, initialize=0)
        m.dxdt = DerivativeVar(m.x, wrt=m.time)

        def diff_eq_rule(m, t):
            return m.dxdt[t] == m.x[t] ** 2 - m.u[t]

        m.diff_eq = Constraint(m.time, rule=diff_eq_rule)

        discretizer = TransformationFactory("dae.finite_difference")
        discretizer.apply_to(m, nfe=nfe, scheme="BACKWARD")
        return m

    @pytest.mark.unit
    def test_fixed_state_variable(self):
        """
        The PETSc utilities raise an exception when a differential variable is fixed.
        While this is okay for the fully time-discretized problem, the integrator
        will not correctly link a fixed differential variable (at non-initial time
        points) with a time derivative.
        """
        m = self.make_rp_example()
        for t in m.time:
            m.x[t].fix(2.0 * t)

        m.u[0].fix(1.0)

        petsc_obj = PETScIntegrator(m, time_set=m.time)
        with pytest.raises(
            RuntimeError,
            match=re.escape(
                "Problem cannot contain a fixed differential variable and "
                "unfixed derivative. Consider either fixing the "
                "corresponding derivative or adding a constraint for the "
                f"differential variable x[10.0] possibly using an "
                "explicit time variable."
            ),
        ):
            _ = petsc_obj.get_timestep_problem(10)

    # TODO Do we want to continue to support this behavior? Because x is given as
    # an explicit function of time, it has effectively become an algebraic equation.
    # TS reports a nonsquare timestep problem, but it solves anyway.
    @pytest.mark.unit
    def test_fixed_time_variable(self):
        """Rather than fixing a differential variable, we can add an explicit time
        variable for the integrator and fix the time variable to the time index
        for the fully discretized problem. This works as an alternative to the
        fixed differential variable.
        """
        m = self.make_rp_example()
        m.t = Var(m.time)

        @m.Constraint(m.time)
        def x_eq(m, t):
            return m.x[t] == 2.0 * m.t[t]

        m.u[0].fix(1.0)
        # For fully discretized fix all times at time index
        m.t[0].fix(m.time.first())

        petsc_obj = PETScIntegrator(
            m,
            m.time,
            time_var=m.t,
            ts_options={
                "--ts_dt": 1,
                "--ts_adapt_type": "none",
            },
        )
        petsc_obj.dae_by_time_element()
        assert value(m.u[10]) == pytest.approx(398)
        assert value(m.x[10]) == pytest.approx(20)

    @pytest.mark.unit
    def test_fixed_derivative(self):
        """Another way to formulate this problem is to fix the derivative. This doesn't
        work for the integrator because it loses the association between the
        differential variable and its derivative.  Since for users, the result of
        this formulation may be unexpected, the PETSc utilities will raise an
        exception if derivatives are fixed to anything other than 0.
        """
        m = self.make_rp_example()
        m.u[0].fix(1.0)
        m.dxdt[:].fix(2.0)

        petsc_obj = PETScIntegrator(m, m.time)
        with pytest.raises(
            RuntimeError,
            match=re.escape(
                "dxdt[10.0] is fixed to a nonzero value 2.0. This is "
                "most likely a modeling error. Instead of fixing the "
                "derivative consider adding a constraint like "
                "dxdt = constant"
            ),
        ):
            petsc_obj.get_timestep_problem(10)

    @pytest.mark.unit
    def test_constrained_derivative(self):
        """Rather than fixing the derivative, we can add a constraint to set the
        derivative.  This should work as intended for both the fully
        time-discretized problem and integrator.
        """
        m = self.make_rp_example()

        @m.Constraint(m.time)
        def ramp_eq(m, t):
            return m.dxdt[t] == 2.0

        m.u[0].fix(1.0)
        m.x[0].fix(0.0)
        m.ramp_eq[0].deactivate()

        petsc_obj = PETScIntegrator(
            m,
            time_set=m.time,
            ts_options={
                "--ts_dt": 1,
                "--ts_adapt_type": "none",
            },
        )
        petsc_obj.dae_by_time_element()
        assert value(m.u[10]) == pytest.approx(398)
        assert value(m.x[10]) == pytest.approx(20)

    @pytest.mark.unit
    def test_ramping_subset(self):
        """Test ramping for a subset of time points"""
        m = self.make_rp_example(nfe=10)
        m.ramp = Var(m.time)

        @m.Constraint(m.time)
        def ramp_eq(m, t):
            return m.dxdt[t] == m.ramp[t]

        m.x[0].fix(0.0)
        for t, v in m.ramp.items():
            if t > 3 and t <= 5:
                v.fix(4)
            else:
                v.fix(0)

        petsc_obj = PETScIntegrator(
            m,
            time_set=m.time,
            ts_options={
                "--ts_dt": 1,
                "--ts_adapt_type": "none",
                "--ts_save_trajectory": 1,
            },
        )
        petsc_obj.dae_by_time_element(between=[0, 3, 5, 10])

        for t in m.time:
            if t < 4:
                assert value(m.ramp[t]) == 0
                assert value(m.x[t]) == pytest.approx(0)
                assert value(m.u[t]) == pytest.approx(0)
            elif t == 4:
                assert value(m.ramp[t]) == 4
                assert value(m.x[t]) == pytest.approx(4)
                assert value(m.u[t]) == pytest.approx(12)
            elif t == 5:
                assert value(m.ramp[t]) == 4
                assert value(m.x[t]) == pytest.approx(8)
                assert value(m.u[t]) == pytest.approx(60)
            else:  # t > 5
                assert value(m.ramp[t]) == 0
                assert value(m.x[t]) == pytest.approx(8)
                assert value(m.u[t]) == pytest.approx(64)

    @pytest.mark.unit
    def test_between_not_in_time_set(self):
        m = self.make_rp_example()
        petsc_obj = PETScIntegrator(
            m,
            time_set=m.time,
            ts_options={
                "--ts_dt": 1,
                "--ts_adapt_type": "none",
                "--ts_save_trajectory": 1,
            },
        )
        with pytest.raises(
            ValueError,
            match=re.escape(
                "Elements of the 'between' argument must be in the time set"
            ),
        ):
            petsc_obj.dae_by_time_element(between=[0, 5, 10])

    @pytest.mark.unit
    def test_calculate_derivatives(self):
        m = self.make_rp_example(time_set=[0, 1, 2], nfe=2)

        m.u.fix(1.0)
        m.x[0].fix(0.0)
        m.diff_eq[0].deactivate()

        petsc_obj = PETScIntegrator(
            m,
            time_set=m.time,
            ts_options={
                "--ts_type": "beuler",
                "--ts_dt": 3e-2,
            },
            calculate_derivatives=True,
        )
        res = petsc_obj.dae_by_time_element(
            between=[0.0, 2.0],
            skip_initial=True,  # With u and x fixed, no variables to solve for at t0
        )
        # No value assigned at 0 because there's no corresponding discretization equation
        assert m.dxdt[0].value is None
        # It should backfill values for interpolated points like t=1
        assert value(m.dxdt[1]) == pytest.approx(-0.7580125427537873, rel=1e-3)
        assert value(m.dxdt[2]) == pytest.approx(-0.2034733131369807, rel=1e-3)

    @pytest.mark.unit
    def test_calculate_derivatives_integrate_first_half(self):
        m = self.make_rp_example(time_set=[0.0, 1.0, 2.0], nfe=2)

        m.u.fix(1.0)
        m.x[0].fix(0.0)
        m.diff_eq[0].deactivate()

        petsc_obj = PETScIntegrator(
            m,
            time_set=m.time,
            ts_options={
                "--ts_type": "beuler",
                "--ts_dt": 3e-2,
            },
            calculate_derivatives=True,
        )
        res = petsc_obj.dae_by_time_element(
            between=[0.0, 1.0],
            skip_initial=True,  # With u and x fixed, no variables to solve for at t0
        )
        # No value assigned at 0 because there's no corresponding discretization equation
        assert m.dxdt[0].value is None
        assert value(m.dxdt[1]) == pytest.approx(-0.7580125427537873, rel=1e-3)
        assert m.dxdt[2].value is None

    @pytest.mark.unit
    def test_calculate_derivatives_integrate_second_half(self):
        m = self.make_rp_example(time_set=[0.0, 1.0, 2.0], nfe=2)

        # Test expects None values
        for t in m.time:
            m.x[t].set_value(None)
            m.u[t].set_value(None)
            m.dxdt[t].set_value(None)

        m.u.fix(1.0)
        m.x[1].fix(-0.7580125427537873)
        m.diff_eq[0].deactivate()

        petsc_obj = PETScIntegrator(
            m,
            time_set=m.time,
            ts_options={
                "--ts_type": "beuler",
                "--ts_dt": 3e-2,
            },
            calculate_derivatives=True,
        )

        res = petsc_obj.dae_by_time_element(
            between=[1.0, 2.0],
            skip_initial=True,  # With u and x fixed, no variables to solve for at t0
        )
        # No value assigned at 0 because there's no corresponding discretization equation
        assert m.dxdt[0].value is None
        assert m.dxdt[1].value is None
        assert value(m.dxdt[2]) == pytest.approx(-0.2034733131369807, rel=1e-3)


@pytest.mark.unit
@pytest.mark.skipif(not petsc_available(), reason="PETSc solver not available")
def test_mixed_derivative_exception():
    # Test exception when mixed space/time derivatives appear in problem
    m = ConcreteModel()

    m.t = ContinuousSet(bounds=(0, 1))
    m.x = ContinuousSet(bounds=(0, 1))

    m.temperature = Var(m.t, m.x)
    m.d2T_dtdx = DerivativeVar(m.temperature, wrt=(m.t, m.x))

    # This problem isn't well-posed, but that shouldn't matter
    @m.Constraint(m.t, m.x)
    def constraint_eqn(b, t, x):
        return b.d2T_dtdx[t, x] == b.temperature[t, x] + 1

    discretizer = TransformationFactory("dae.finite_difference")
    discretizer.apply_to(m, nfe=3, scheme="BACKWARD", wrt=m.t)
    discretizer.apply_to(m, nfe=3, scheme="BACKWARD", wrt=m.x)

    with pytest.raises(
        NotImplementedError,
        match=re.escape(
            "IDAES presently does not support PETSc for second order or higher derivatives "
            "that are differentiated at least once with respect to time. Please reformulate your model so "
            "it does not contain such a derivative (such as by introducing intermediate variables)."
        ),
    ):
        petsc_obj = PETScIntegrator(
            m,
            time_set=m.t,
            ts_options={
                "--ts_dt": 1,
                "--ts_adapt_type": "none",
                "--ts_save_trajectory": 1,
            },
        )


@pytest.mark.unit
@pytest.mark.skipif(not petsc_available(), reason="PETSc solver not available")
def test_petsc_skip_initial_solve():
    """
    Ensure skipping the initial solution works as intended
    """
    m, y1, y2, y3, y4, y5, y6 = dae_with_non_time_indexed_constraint()

    m.y[0, 6].value = 0.35999964

    m.r[0, 1].value = 0.025487429806082887
    m.r[0, 2].value = 0
    m.r[0, 3].value = 0
    m.r[0, 4].value = 1.9580400000000004e-06
    m.r[0, 5].value = 0.0019090002227229192

    m.H.value = 737
    m.Fin[0].value = -2.9149253731343243e-05

    petsc_obj = PETScIntegrator(
        m,
        time_set=m.t,
        ts_options={
            "--ts_type": "cn",  # Crank–Nicolson
            "--ts_adapt_type": "basic",
            "--ts_dt": 0.01,
            "--ts_save_trajectory": 1,
            "--ts_trajectory_type": "visualization",
        },
    )
    res = petsc_obj.dae_by_time_element(skip_initial=True)

    assert pytest.approx(y1, rel=1e-3) == value(m.y[m.t.last(), 1])
    assert pytest.approx(y2, rel=1e-3) == value(m.y[m.t.last(), 2])
    assert pytest.approx(y3, rel=1e-3) == value(m.y[m.t.last(), 3])
    assert pytest.approx(y4, rel=1e-3) == value(m.y[m.t.last(), 4])
    assert pytest.approx(y5, rel=1e-3) == value(m.y[m.t.last(), 5])
    assert pytest.approx(y6, rel=1e-3) == value(m.y[m.t.last(), 6])


@pytest.mark.unit
@pytest.mark.skipif(not petsc_available(), reason="PETSc solver not available")
def test_petsc_traj_previous():
    """
    Make sure trajectory concatenation works as intended
    """
    m, y1, y2, y3, y4, y5, y6 = dae_with_non_time_indexed_constraint(nfe=10)

    # First get a reference trajectory from a single solve
    petsc_obj = PETScIntegrator(
        m,
        time_set=m.t,
        ts_options={
            "--ts_type": "cn",  # Crank–Nicolson
            "--ts_adapt_type": "basic",
            "--ts_dt": 0.01,
            "--ts_save_trajectory": 1,
        },
    )
    res0 = petsc_obj.dae_by_time_element(between=[m.t.first(), m.t.last()])
    tj0 = res0.trajectory

    # Next do it in two steps
    res = petsc_obj.dae_by_time_element(between=[m.t.first(), m.t.at(4)])

    # Fix initial condition for second leg of trajectory
    for j in range(1, 6):
        m.y[m.t.at(4), j].fix()

    res = petsc_obj.dae_by_time_element(
        between=[m.t.at(4), m.t.last()],
        previous_trajectory=res.trajectory,
    )

    assert pytest.approx(y1, rel=1e-3) == value(m.y[m.t.last(), 1])
    assert pytest.approx(y2, rel=1e-3) == value(m.y[m.t.last(), 2])
    assert pytest.approx(y3, rel=1e-3) == value(m.y[m.t.last(), 3])
    assert pytest.approx(y4, rel=1e-3) == value(m.y[m.t.last(), 4])
    assert pytest.approx(y5, rel=1e-3) == value(m.y[m.t.last(), 5])
    assert pytest.approx(y6, rel=1e-3) == value(m.y[m.t.last(), 6])

    tj = res.trajectory
    assert tj.get_vec(m.y[180, 1])[-1] == pytest.approx(y1, rel=1e-3)
    assert tj.get_vec("_time")[-1] == pytest.approx(180)
    y1_trj = tj.interpolate_vec(m.t, m.y[180, 1])
    y4_trj = tj.interpolate_vec(m.t, m.y[180, 4])
    for i, t in enumerate(m.t):
        assert y1_trj[i] == pytest.approx(value(m.y[t, 1]))
        assert y4_trj[i] == pytest.approx(value(m.y[t, 4]))

    t_vec = np.linspace(0, 180, 101)
    y2_tj0 = tj0.interpolate_vec(t_vec, m.y[180, 2])
    y5_tj0 = tj0.interpolate_vec(t_vec, m.y[180, 5])
    y2_tj = tj.interpolate_vec(t_vec, m.y[180, 2])
    y5_tj = tj.interpolate_vec(t_vec, m.y[180, 5])

    for i, t in enumerate(t_vec):
        assert y2_tj[i] == pytest.approx(y2_tj0[i], abs=1e-4)
        assert y5_tj[i] == pytest.approx(y5_tj0[i], abs=1e-4)


@pytest.mark.unit
@pytest.mark.skipif(not petsc_available(), reason="PETSc solver not available")
def test_not_ContinuousSet():
    m = ConcreteModel()

    m.time = Set(initialize=(0.0, 10.0))
    m.x = Var(m.time)
    m.u = Var(m.time)

    with pytest.raises(ValueError, match="ContinuousSet"):
        PETScIntegrator(
            m,
            time_set=m.time,
        )


@pytest.mark.unit
@pytest.mark.skipif(not petsc_available(), reason="PETSc solver not available")
def test_not_discretized():
    m = ConcreteModel()

    m.time = ContinuousSet(initialize=(0.0, 10.0))
    m.x = Var(m.time)
    m.u = Var(m.time)
    m.dxdt = DerivativeVar(m.x, wrt=m.time)

    def diff_eq_rule(m, t):
        return m.dxdt[t] == m.x[t] ** 2 - m.u[t]

    m.diff_eq = Constraint(m.time, rule=diff_eq_rule)

    with pytest.raises(RuntimeError, match="discretized"):
        PETScIntegrator(
            m,
            time_set=m.time,
        )


@pytest.mark.unit
@pytest.mark.skipif(not petsc_available(), reason="PETSc solver not available")
def test_calculate_derivatives_collocation_lr():
    m = ConcreteModel()

    m.time = ContinuousSet(initialize=(0.0, 1.0, 2.0))
    m.x = Var(m.time)
    m.u = Var(m.time)
    m.dxdt = DerivativeVar(m.x, wrt=m.time)

    def diff_eq_rule(m, t):
        return m.dxdt[t] == m.x[t] ** 2 - m.u[t]

    m.diff_eq = Constraint(m.time, rule=diff_eq_rule)

    discretizer = TransformationFactory("dae.collocation")
    discretizer.apply_to(m, nfe=2, scheme="LAGRANGE-RADAU")

    m.u.fix(1.0)
    m.x[0].fix(0.0)

    petsc_obj = PETScIntegrator(
        m,
        time_set=m.time,
        ts_options={
            "--ts_type": "beuler",
            "--ts_dt": 3e-2,
        },
        calculate_derivatives=True,
    )
    res = petsc_obj.dae_by_time_element(
        between=[0.0, 2.0],
        skip_initial=True,  # With u and x fixed, no variables to solve for at t0
    )

    assert m.dxdt[0].value is None
    # It should backfill values for interpolated points like t=1
    assert value(m.dxdt[1]) == pytest.approx(-0.3787483909827891, rel=1e-3)
    assert value(m.dxdt[2]) == pytest.approx(-0.07952012649572815, rel=1e-3)


@pytest.mark.unit
@pytest.mark.skipif(not petsc_available(), reason="PETSc solver not available")
def test_exception_forward_difference():
    m = ConcreteModel()

    m.time = ContinuousSet(initialize=(0.0, 1.0, 2.0))
    m.x = Var(m.time)
    m.u = Var(m.time)
    m.dxdt = DerivativeVar(m.x, wrt=m.time)

    def diff_eq_rule(m, t):
        return m.dxdt[t] == m.x[t] ** 2 - m.u[t]

    m.diff_eq = Constraint(m.time, rule=diff_eq_rule)

    discretizer = TransformationFactory("dae.finite_difference")
    discretizer.apply_to(m, nfe=2, scheme="FORWARD")

    m.u.fix(1.0)
    m.x[0].fix(0.0)

    with pytest.raises(
        RuntimeError,
        match="The PETSc interface supports only the backward difference and Lagrange-Radau "
        "collocation schemes. Found instead the FORWARD Difference scheme.",
    ):
        _ = PETScIntegrator(
            m,
            time_set=m.time,
            ts_options={
                "--ts_type": "beuler",
                "--ts_dt": 3e-2,
            },
            calculate_derivatives=True,
        )


@pytest.mark.unit
@pytest.mark.skipif(not petsc_available(), reason="PETSc solver not available")
def test_exception_central_difference():
    m = ConcreteModel()

    m.time = ContinuousSet(initialize=(0.0, 1.0, 2.0))
    m.x = Var(m.time)
    m.u = Var(m.time)
    m.dxdt = DerivativeVar(m.x, wrt=m.time)

    def diff_eq_rule(m, t):
        return m.dxdt[t] == m.x[t] ** 2 - m.u[t]

    m.diff_eq = Constraint(m.time, rule=diff_eq_rule)

    discretizer = TransformationFactory("dae.finite_difference")
    discretizer.apply_to(m, nfe=2, scheme="CENTRAL")

    m.u.fix(1.0)
    m.x[0].fix(0.0)

    with pytest.raises(
        RuntimeError,
        match=re.escape(
            "The PETSc interface supports only the backward difference and Lagrange-Radau "
            "collocation schemes. Found instead the CENTRAL Difference scheme."
        ),
    ):
        _ = PETScIntegrator(
            m,
            time_set=m.time,
            ts_options={
                "--ts_type": "beuler",
                "--ts_dt": 3e-2,
            },
            calculate_derivatives=True,
        )


@pytest.mark.unit
@pytest.mark.skipif(not petsc_available(), reason="PETSc solver not available")
def test_exception_collocation_ll():
    m = ConcreteModel()

    m.time = ContinuousSet(initialize=(0.0, 1.0, 2.0))
    m.x = Var(m.time)
    m.u = Var(m.time)
    m.dxdt = DerivativeVar(m.x, wrt=m.time)

    def diff_eq_rule(m, t):
        return m.dxdt[t] == m.x[t] ** 2 - m.u[t]

    m.diff_eq = Constraint(m.time, rule=diff_eq_rule)

    discretizer = TransformationFactory("dae.collocation")
    discretizer.apply_to(m, nfe=2, scheme="LAGRANGE-LEGENDRE")

    m.u.fix(1.0)
    m.x[0].fix(0.0)

    with pytest.raises(
        RuntimeError,
        match=re.escape(
            "The PETSc interface supports only the backward difference and Lagrange-Radau "
            "collocation schemes. Found instead the LAGRANGE-LEGENDRE scheme."
        ),
    ):
        _ = PETScIntegrator(
            m,
            time_set=m.time,
            ts_options={
                "--ts_type": "beuler",
                "--ts_dt": 3e-2,
            },
            calculate_derivatives=True,
        )


def make_index_reduction_model(disc_method, nfe, reduce_index):
    m = ConcreteModel()
    m.time = ContinuousSet(initialize=(0.0, 1.0))
    m.comp_set = ["A", "B"]
    m.n = Var(
        m.time,
        m.comp_set,
        initialize=1,
    )
    m.dn_dt = DerivativeVar(m.n, initialize=0)
    m.C_in = Param(m.time, m.comp_set, initialize=1, mutable=True)
    m.F_in = Param(m.time, initialize=1, mutable=True)
    m.V = Param(initialize=1)
    m.V_bar = Param(m.comp_set, initialize={"A": 1, "B": 2})
    m.F_out = Var(m.time, initialize=1)

    @m.Constraint(m.time, m.comp_set)
    def material_balance_eqn(b, t, j):
        return b.dn_dt[t, j] == (
            b.C_in[t, j] * b.F_in[t] - b.n[t, j] / b.V * b.F_out[t]
        )

    @m.Constraint(m.time)
    def volume_eqn(b, t):
        return sum(b.n[t, j] * b.V_bar[j] for j in b.comp_set) == b.V

    discretizer = TransformationFactory(disc_method)
    if disc_method == "dae.collocation":
        discretizer.apply_to(m, nfe=nfe, ncp=2, scheme="LAGRANGE-RADAU")
    else:
        discretizer.apply_to(m, nfe=nfe, scheme="BACKWARD")

    m.n[0, "A"].fix()
    if reduce_index:
        m.dn_dt_disc_eq[:, "B"].deactivate()

        @m.Constraint(m.time)
        def d_volume_dt_eqn(b, t):
            return sum(b.dn_dt[t, j] * b.V_bar[j] for j in b.comp_set) == 0

    else:
        m.n[0, "B"].fix()

    m.n[:, "A"].set_value(1 / 2)
    m.n[:, "B"].set_value(1 / 4)
    m.C_in[:, "A"].set_value(1 / 2)
    m.C_in[:, "B"].set_value(1 / 2)

    return m


@pytest.mark.parametrize("reduce_index", [False, True])
@pytest.mark.parametrize("disc_method", ["dae.finite_difference", "dae.collocation"])
@pytest.mark.parametrize("nfe", [1, 2])
@pytest.mark.unit
@pytest.mark.skipif(not petsc_available(), reason="PETSc solver not available")
def test_find_discretization_equations_deactivated(disc_method, nfe, reduce_index):
    m = make_index_reduction_model(disc_method, nfe, reduce_index)
    out = get_derivative_differential_vardata_map(m, m.time)

    # Make sure that the derivative map gets generated correctly
    for t in m.time:
        if t == 0:
            assert m.dn_dt[t, "A"] not in out
        else:
            assert out[m.dn_dt[t, "A"]]["diff_var"] is m.n[t, "A"]
            assert out[m.dn_dt[t, "A"]]["disc_eq"] is m.dn_dt_disc_eq[t, "A"]
        if reduce_index or t == 0:
            assert m.dn_dt[t, "B"] not in out
        else:
            assert out[m.dn_dt[t, "B"]]["diff_var"] is m.n[t, "B"]
            assert out[m.dn_dt[t, "B"]]["disc_eq"] is m.dn_dt_disc_eq[t, "B"]

    integrator = PETScIntegrator(
        m,
        time_set=m.time,
        ts_options={
            "--ts_type": "beuler",
            "--ts_dt": 3e-2,
        },
        calculate_derivatives=True,
    )

    ic_problem = integrator.get_initial_condition_problem(0)
    dt = DiagnosticsToolbox(ic_problem)
    if reduce_index:
        dt.assert_no_structural_warnings()
    else:
        out = dt.get_dulmage_mendelsohn_partition()
        # 3 for the constraints in the overconstrained set
        # 0 for the first block
        # 0 for the first constraint
        assert out[3][0][0] is m.volume_eqn[0]

    for t in m.time:
        if t == m.time.first():
            continue
        tstep_problem, diff_vars, fixedness = integrator.get_timestep_problem(t)
        diff_vars = ComponentSet(diff_vars)
        dt = DiagnosticsToolbox(tstep_problem)
        for var in diff_vars:
            var.fix()
        assert m.n[t, "A"] in diff_vars
        if reduce_index:
            assert len(diff_vars) == 1
            dt.assert_no_structural_warnings()
        else:
            assert len(diff_vars) == 2
            assert m.n[t, "B"] in diff_vars
            out = dt.get_dulmage_mendelsohn_partition()
            # 3 for the constraints in the overconstrained set
            # 0 for the first block
            # 0 for the first constraint
            assert out[3][0][0] is m.volume_eqn[t]

        from_json(m, sd=fixedness, wts=StoreSpec.isfixed())

    def approx(x):
        return pytest.approx(x, rel=1e-2)

    if reduce_index:
        results = integrator.dae_by_time_element(
            between=[m.time.first(), m.time.last()],
        )
        assert_optimal_termination(results.results[0])
        assert value(m.n[0, "A"]) == 0.5
        assert value(m.n[0, "B"]) == 0.25
        assert value(m.n[1, "A"]) == approx(0.372149)
        assert value(m.n[1, "B"]) == approx(0.313925)
