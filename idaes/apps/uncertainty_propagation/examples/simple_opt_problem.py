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

"""
Consider the following optimization problem.

min f:  p1*x1+ p2*(x2^2) + p1*p2
         s.t  c1: x1 + x2 = p1
              c2: x2 + x3 = p2
              0 <= x1, x2, x3 <= 10
              p1 = 10
              p2 = 5

Variables = (x1, x2, x3)
Parameters (fixed variables) = (p1, p2)
"""

import pyomo.environ as pyo
import numpy as np
import pytest
from idaes.apps.uncertainty_propagation.uncertainties import propagate_uncertainty

if __name__ == "__main__":

    ### Create optimization model
    m = pyo.ConcreteModel()
    m.dual = pyo.Suffix(direction=pyo.Suffix.IMPORT)

    # Define variables
    # m.x1 = pyo.Var(bounds=(0,10))
    # m.x2 = pyo.Var(bounds=(0,10))
    # m.x3 = pyo.Var(bounds=(0,10))

    m.x1 = pyo.Var()
    m.x2 = pyo.Var()
    m.x3 = pyo.Var()

    # Define parameters
    # m.p1 = pyo.Param(initialize=10, mutable=True)
    # m.p2 = pyo.Param(initialize=5, mutable=True)
    # Important Tip: The uncertain parameters need to be defined at Pyomo variables
    m.p1 = pyo.Var(initialize=10)
    m.p2 = pyo.Var(initialize=5)
    # Important Tip: We fix these to first solve with Ipopt. We unfix below
    # before using the uncertainty propagation toolbox.
    m.p1.fix()
    m.p2.fix()

    # Define constraints
    m.con1 = pyo.Constraint(expr=m.x1 + m.x2 - m.p1 == 0)
    m.con2 = pyo.Constraint(expr=m.x2 + m.x3 - m.p2 == 0)

    # Define objective
    m.obj = pyo.Objective(
        expr=m.p1 * m.x1 + m.p2 * (m.x2**2) + m.p1 * m.p2, sense=pyo.minimize
    )

    ### Solve optimization model
    opt = pyo.SolverFactory("ipopt", tee=True)
    opt.solve(m)

    ### Inspect solution
    print("Numeric solution:")
    print("x1 =", m.x1())
    print("x2 =", m.x2())
    print("x3 =", m.x3())
    print(m.dual.display())

    ### Analytic solution
    """
    At the optimal solution, none of the bounds are active. As long as the active set
    does not change (i.e., none of the bounds become active), the
    first order optimality conditions reduce to a simple linear system.
    """

    # dual variables (multipliers)
    v2_ = 0
    v1_ = m.p1()

    # primal variables
    x2_ = (v1_ + v2_) / (2 * m.p2())
    x1_ = m.p1() - x2_
    x3_ = m.p2() - x2_

    print("\nAnalytic solution:")
    print("x1 =", x1_)
    print("x2 =", x2_)
    print("x3 =", x3_)
    print("v1 =", v1_)
    print("v2 =", v2_)

    ### Analytic sensitivity
    """
    Using the analytic solution above, we can compute the sensitivies of x and v to
    perturbations in p1 and p2.


    The matrix dx_dp constains the sensitivities of x to perturbations in p
    """

    # Initialize sensitivity matrix Nx x Np
    # Rows: variables x
    # Columns: parameters p
    dx_dp = np.zeros((3, 2))

    # dx2/dp1 = 1/(2 * p2)
    dx_dp[1, 0] = 1 / (2 * m.p2())

    # dx2/dp2 = -(v1 + v2)/(2 * p2**2)
    dx_dp[1, 1] = -(v1_ + v2_) / (2 * m.p2() ** 2)

    # dx1/dp1 = 1 - dx2/dp1
    dx_dp[0, 0] = 1 - dx_dp[1, 0]

    # dx1/dp2 = 0 - dx2/dp2
    dx_dp[0, 1] = 0 - dx_dp[1, 1]

    # dx3/dp1 = 1 - dx2/dp1
    dx_dp[2, 0] = 0 - dx_dp[1, 0]

    # dx3/dp2 = 0 - dx2/dp2
    dx_dp[2, 1] = 1 - dx_dp[1, 1]

    print("\n\ndx/dp =\n", dx_dp)

    """
    Similarly, we can compute the gradients df_dx, df_dp
    and Jacobians dc_dx, dc_dp
    """

    # Initialize 1 x 3 array to store (\partial f)/(\partial x)
    # Elements: variables x
    df_dx = np.zeros(3)

    # df/dx1 = p1
    df_dx[0] = m.p1()

    # df/dx2 = p2
    df_dx[1] = 2 * m.p2() * x2_

    # df/dx3 = 0

    print("\n\ndf/dx =\n", df_dx)

    # Initialize 1 x 2 array to store (\partial f)/(\partial p)
    # Elements: parameters p
    df_dp = np.zeros(2)

    # df/dxp1 = x1 + p2
    df_dp[0] = x1_ + m.p2()

    # df/dp2 = x2**2 + p1
    df_dp[1] = x2_**2 + m.p1()

    print("\n\ndf/dp =\n", df_dp)

    # Initialize 2 x 3 array to store (\partial c)/(\partial x)
    # Rows: constraints c
    # Columns: variables x
    dc_dx = np.zeros((2, 3))

    # dc1/dx1 = 1
    dc_dx[0, 0] = 1

    # dc1/dx2 = 1
    dc_dx[0, 1] = 1

    # dc2/dx2 = 1
    dc_dx[1, 1] = 1

    # dc2/dx3 = 1
    dc_dx[1, 2] = 1

    # Remaining entries are 0

    print("\n\ndc/dx =\n", dc_dx)

    # Initialize 2 x 2 array to store (\partial c)/(\partial x)
    # Rows: constraints c
    # Columns: variables x
    dc_dp = np.zeros((2, 2))

    # dc1/dp1 = -1
    dc_dp[0, 0] = -1

    # dc2/dp2 = -1
    dc_dp[1, 1] = -1

    # Remaining entries are 0

    print("\n\ndc/dp =\n", dc_dp)

    ### Uncertainty propagation
    """
    Now lets test the uncertainty propagation package. We will assume p has covariance
    sigma_p = [[2, 0], [0, 1]]
    """

    ## Prepare inputs

    # Covariance matrix
    sigma_p = np.array([[2, 0], [0, 1]])

    # Nominal values for uncertain parameters
    theta = {"p1": m.p1(), "p2": m.p2()}

    # Names of uncertain parameters
    theta_names = ["p1", "p2"]

    # Important to unfix the parameters!
    # Otherwise k_aug will complain about too few degrees of freedom
    m.p1.unfix()
    m.p2.unfix()

    ## Run sensitivity toolbox
    results = propagate_uncertainty(m, theta, sigma_p, theta_names)

    ## Compute uncertainty propagation by hand
    # We now compute the expected result with our
    # analytic solution.

    # Take 1: use the formula in the documentation
    tmp_f = df_dp + df_dx @ dx_dp
    sigma_f = tmp_f @ sigma_p @ tmp_f.transpose()
    print("\nsigma_f = ", sigma_f)

    tmp_c = dc_dp + dc_dx @ dx_dp
    sigma_c = tmp_c @ sigma_p @ tmp_c.transpose()
    print("\nsigma_c = \n", sigma_c)

    # Take 2: use an equivalent formula
    print("\n\n")
    grad_f = np.hstack((df_dx, df_dp))
    dsdp = np.vstack((dx_dp, np.eye(2)))
    sigma_f_take_2 = grad_f @ dsdp @ sigma_p @ dsdp.transpose() @ grad_f.transpose()
    # This matches
    print("sigma_f (take 2) = ", sigma_f_take_2)

    dc_ds = np.hstack((dc_dx, dc_dp))
    sigma_c_take_2 = dc_ds @ dsdp @ sigma_p @ dsdp.transpose() @ dc_ds.transpose()
    # This matches
    print("sigma_c (take 2) = ", sigma_c_take_2)

    ## Check results
    # We now check the results for the uncertainty propagation toolbox

    # Check the order from the uncertainty propagation toolbox is as expected
    # TODO: make this test more automated to use this information
    # to reshuffle the analytic solution
    assert results.col == ["x1", "x2", "p1", "p2", "x3"]
    assert results.row == ["con1", "con2", "obj"]

    # Define indices to compare results to analytic solution
    # TODO: make this more robust to use the results above
    var_idx = np.array([True, True, False, False, True])
    theta_idx = np.array([False, False, True, True, False])

    # Check the gradient of the objective w.r.t. x matches
    np.testing.assert_array_almost_equal(results.gradient_f[var_idx], np.array(df_dx))

    # Check the gradient of the objective w.r.t. p (parameters) matches
    np.testing.assert_array_almost_equal(results.gradient_f[theta_idx], np.array(df_dp))

    # Check the Jacobian of the constraints w.r.t. x matches
    np.testing.assert_array_almost_equal(
        results.gradient_c.toarray()[:, var_idx], np.array(dc_dx)
    )

    # Check the Jacobian of the constraints w.r.t. p (parameters) matches
    np.testing.assert_array_almost_equal(
        results.gradient_c.toarray()[:, theta_idx], np.array(dc_dp)
    )

    # Check the NLP sensitivity results for the variables (x) matches
    np.testing.assert_array_almost_equal(
        results.dsdp.toarray()[var_idx, :], np.array(dx_dp)
    )

    # Check the NLP sensitivity results for the parameters (p) matches
    np.testing.assert_array_almost_equal(
        results.dsdp.toarray()[theta_idx, :], np.array([[1, 0], [0, 1]])
    )

    # Check the uncertainty propagation results for the constrains matches
    np.testing.assert_array_almost_equal(results.propagation_c, np.sum(sigma_c))

    # Check the uncertainty propagation results for the objective matches
    assert results.propagation_f == pytest.approx(sigma_f)
