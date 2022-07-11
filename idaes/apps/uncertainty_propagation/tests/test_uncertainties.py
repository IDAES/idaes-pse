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
import sys
import os
from unittest.mock import patch

sys.path.append(os.path.abspath(".."))  # current folder is ~/tests
import numpy as np
import pandas as pd
from scipy import sparse
import pytest
from pytest import approx
from idaes.apps.uncertainty_propagation.uncertainties import (
    quantify_propagate_uncertainty,
    propagate_uncertainty,
    clean_variable_name,
)
from pyomo.opt import SolverFactory
from pyomo.environ import *
import pyomo.contrib.parmest.parmest as parmest

ipopt_available = SolverFactory("ipopt").available()
kaug_available = SolverFactory("k_aug").available()
dotsens_available = SolverFactory("dot_sens").available()


@pytest.mark.skipif(not ipopt_available, reason="The 'ipopt' command is not available")
@pytest.mark.skipif(not kaug_available, reason="The 'k_aug' command is not available")
@pytest.mark.skipif(
    not dotsens_available, reason="The 'dot_sens' command is not available"
)
class TestUncertaintyPropagation:
    @pytest.mark.unit
    def test_quantify_propagate_uncertainty1(self):
        """
        It tests the function quantify_propagate_uncertainty with rooney & biegler's model.
        """
        from idaes.apps.uncertainty_propagation.examples.rooney_biegler import (
            rooney_biegler_model,
            rooney_biegler_model_opt,
        )

        variable_name = ["asymptote", "rate_constant"]
        data = pd.DataFrame(
            data=[[1, 8.3], [2, 10.3], [3, 19.0], [4, 16.0], [5, 15.6], [7, 19.8]],
            columns=["hour", "y"],
        )

        def SSE(model, data):
            expr = sum(
                (data.y[i] - model.response_function[data.hour[i]]) ** 2
                for i in data.index
            )
            return expr

        results = quantify_propagate_uncertainty(
            rooney_biegler_model, rooney_biegler_model_opt, data, variable_name, SSE
        )

        assert results.obj == approx(4.331711213656886)
        np.testing.assert_array_almost_equal(
            results.theta, [19.142575284617866, 0.53109137696521]
        )
        assert list(results.theta.keys()) == ["asymptote", "rate_constant"]
        np.testing.assert_array_almost_equal(results.gradient_f, [0.99506259, 0.945148])
        assert list(results.propagation_c) == []
        np.testing.assert_array_almost_equal(
            results.dsdp.toarray(), [[1.0, 0.0], [0.0, 1.0]]
        )
        np.testing.assert_array_almost_equal(
            results.cov, np.array([[6.30579403, -0.4395341], [-0.4395341, 0.04193591]])
        )
        assert results.propagation_f == pytest.approx(5.45439337747349)

    @pytest.mark.component
    def test_quantify_propagate_uncertainty2(self):
        """
        This is the same test as test_quantify_propagate_uncertainty1,
        but with the second argument of quantify_propagate_uncertainty as Pyomo Concrete Model.
        """
        from idaes.apps.uncertainty_propagation.examples.rooney_biegler import (
            rooney_biegler_model,
        )

        variable_name = ["asymptote", "rate_constant"]
        data = pd.DataFrame(
            data=[[1, 8.3], [2, 10.3], [3, 19.0], [4, 16.0], [5, 15.6], [7, 19.8]],
            columns=["hour", "y"],
        )

        def SSE(model, data):
            expr = sum(
                (data.y[i] - model.response_function[data.hour[i]]) ** 2
                for i in data.index
            )
            return expr

        model_uncertain = ConcreteModel()
        model_uncertain.asymptote = Var(initialize=15)
        model_uncertain.rate_constant = Var(initialize=0.5)
        model_uncertain.obj = Objective(
            expr=model_uncertain.asymptote
            * (1 - exp(-model_uncertain.rate_constant * 10)),
            sense=minimize,
        )

        results = quantify_propagate_uncertainty(
            rooney_biegler_model, model_uncertain, data, variable_name, SSE
        )

        assert results.obj == approx(4.331711213656886)
        np.testing.assert_array_almost_equal(
            results.theta, [19.142575284617866, 0.53109137696521]
        )
        assert list(results.theta.keys()) == ["asymptote", "rate_constant"]
        np.testing.assert_array_almost_equal(results.gradient_f, [0.99506259, 0.945148])
        assert list(results.propagation_c) == []
        np.testing.assert_array_almost_equal(
            results.dsdp.toarray(), [[1.0, 0.0], [0.0, 1.0]]
        )
        np.testing.assert_array_almost_equal(
            results.cov, np.array([[6.30579403, -0.4395341], [-0.4395341, 0.04193591]])
        )
        assert results.propagation_f == pytest.approx(5.45439337747349)

    @pytest.mark.component
    def test_propagate_uncertainty(self):
        """
        It tests the function propagate_uncertainty with rooney & biegler's model.
        """
        from idaes.apps.uncertainty_propagation.examples.rooney_biegler import (
            rooney_biegler_model,
        )

        variable_name = ["asymptote", "rate_constant"]
        data = pd.DataFrame(
            data=[[1, 8.3], [2, 10.3], [3, 19.0], [4, 16.0], [5, 15.6], [7, 19.8]],
            columns=["hour", "y"],
        )

        def SSE(model, data):
            expr = sum(
                (data.y[i] - model.response_function[data.hour[i]]) ** 2
                for i in data.index
            )
            return expr

        parmest_class = parmest.Estimator(
            rooney_biegler_model, data, variable_name, SSE
        )
        obj, theta, cov = parmest_class.theta_est(calc_cov=True, cov_n=len(data.index))
        model_uncertain = ConcreteModel()
        model_uncertain.asymptote = Var(initialize=15)
        model_uncertain.rate_constant = Var(initialize=0.5)
        model_uncertain.obj = Objective(
            expr=model_uncertain.asymptote
            * (1 - exp(-model_uncertain.rate_constant * 10)),
            sense=minimize,
        )

        propagate_results = propagate_uncertainty(
            model_uncertain, theta, cov, variable_name
        )

        np.testing.assert_array_almost_equal(
            propagate_results.gradient_f, [0.9950625870024135, 0.9451480001755206]
        )
        assert list(propagate_results.gradient_c) == []
        np.testing.assert_array_almost_equal(
            propagate_results.dsdp.toarray(), [[1.0, 0.0], [0.0, 1.0]]
        )
        assert list(propagate_results.propagation_c) == []
        assert propagate_results.propagation_f == pytest.approx(5.45439337747349)

    @pytest.mark.component
    def test_propagate_uncertainty1(self):
        """
        It tests the function propagate_uncertainty with
        min f:  p1*x1+ p2*(x2^2) + p1*p2
         s.t  c1: x1 + x2 = p1
              c2: x2 + x3 = p2
              0 <= x1, x2, x3 <= 10
              p1 = 10
              p2 = 5

        Variables = (x1, x2, x3)
        Parameters (fixed variables) = (p1, p2)
        """
        ### Create optimization model
        m = ConcreteModel()
        m.dual = Suffix(direction=Suffix.IMPORT)

        m.x1 = Var()
        m.x2 = Var()
        m.x3 = Var()

        # Define parameters
        m.p1 = Var(initialize=10)
        m.p2 = Var(initialize=5)
        m.p1.fix()
        m.p2.fix()

        # Define constraints
        m.con1 = Constraint(expr=m.x1 + m.x2 - m.p1 == 0)
        m.con2 = Constraint(expr=m.x2 + m.x3 - m.p2 == 0)

        # Define objective
        m.obj = Objective(
            expr=m.p1 * m.x1 + m.p2 * (m.x2**2) + m.p1 * m.p2, sense=minimize
        )

        ### Solve optimization model
        opt = SolverFactory("ipopt", tee=True)
        opt.solve(m)

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
        # Initialize 1 x 2 array to store (\partial f)/(\partial p)
        # Elements: parameters p
        df_dp = np.zeros(2)
        # df/dxp1 = x1 + p2
        df_dp[0] = x1_ + m.p2()
        # df/dp2 = x2**2 + p1
        df_dp[1] = x2_**2 + m.p1()

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
        # Initialize 2 x 2 array to store (\partial c)/(\partial x)
        # Rows: constraints c
        # Columns: variables x
        dc_dp = np.zeros((2, 2))
        # dc1/dp1 = -1
        dc_dp[0, 0] = -1
        # dc2/dp2 = -1
        dc_dp[1, 1] = -1

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

        ## Run package
        results = propagate_uncertainty(m, theta, sigma_p, theta_names)

        ## Check results

        tmp_f = df_dp + df_dx @ dx_dp
        sigma_f = tmp_f @ sigma_p @ tmp_f.transpose()

        tmp_c = dc_dp + dc_dx @ dx_dp
        sigma_c = tmp_c @ sigma_p @ tmp_c.transpose()

        # This currently just checks if the order of the outputs did not change
        # TODO: improve test robustness by using this information to set
        # var_idx and theta_idx. This way the test will still work
        # regardless of the order. In other words, the analytic solution needs to be
        # reordered to match the variable/constraint order from
        # this package. Alternately, the results could be converted into a Pandas dataframe
        assert results.col == ["x1", "x2", "p1", "p2", "x3"]
        assert results.row == ["con1", "con2", "obj"]
        var_idx = np.array([True, True, False, False, True])
        theta_idx = np.array([False, False, True, True, False])

        # Check the gradient of the objective w.r.t. x matches
        np.testing.assert_array_almost_equal(
            results.gradient_f[var_idx], np.array(df_dx)
        )

        # Check the gradient of the objective w.r.t. p (parameters) matches
        np.testing.assert_array_almost_equal(
            results.gradient_f[theta_idx], np.array(df_dp)
        )

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

    @pytest.mark.component
    def test_propagate_uncertainty_error(self):
        """
        It tests a TypeError when the modle_uncertian of function propagate_uncertainty is neither python function nor Pyomo ConcreteModel
        """
        from idaes.apps.uncertainty_propagation.examples.rooney_biegler import (
            rooney_biegler_model,
        )

        variable_name = ["asymptote", "rate_constant"]
        data = pd.DataFrame(
            data=[[1, 8.3], [2, 10.3], [3, 19.0], [4, 16.0], [5, 15.6], [7, 19.8]],
            columns=["hour", "y"],
        )

        def SSE(model, data):
            expr = sum(
                (data.y[i] - model.response_function[data.hour[i]]) ** 2
                for i in data.index
            )
            return expr

        parmest_class = parmest.Estimator(
            rooney_biegler_model, data, variable_name, SSE
        )
        obj, theta, cov = parmest_class.theta_est(calc_cov=True, cov_n=len(data.index))
        model_uncertain = ConcreteModel()
        model_uncertain.asymptote = Var(initialize=15)
        model_uncertain.rate_constant = Var(initialize=0.5)
        model_uncertain.obj = Objective(
            expr=model_uncertain.asymptote
            * (1 - exp(-model_uncertain.rate_constant * 10)),
            sense=minimize,
        )
        with pytest.raises(TypeError):
            propagate_results = propagate_uncertainty(1, theta, cov, variable_name)

    @pytest.mark.integration
    def test_quantify_propagate_uncertainty_NRTL(self):
        """
        It tests the function quantify_propagate_uncertainty with IDAES NRTL model.
        """
        from idaes.apps.uncertainty_propagation.examples.NRTL_model_scripts import (
            NRTL_model,
            NRTL_model_opt,
        )

        variable_name = [
            "fs.properties.tau['benzene','toluene']",
            "fs.properties.tau['toluene','benzene']",
        ]
        current_path = os.path.dirname(os.path.realpath(__file__))
        data = pd.read_csv(os.path.join(current_path, "BT_NRTL_dataset.csv"))

        def SSE(model, data):
            expr = (
                float(data["vap_benzene"])
                - model.fs.flash.vap_outlet.mole_frac_comp[0, "benzene"]
            ) ** 2 + (
                float(data["liq_benzene"])
                - model.fs.flash.liq_outlet.mole_frac_comp[0, "benzene"]
            ) ** 2
            return expr * 1e4

        results = quantify_propagate_uncertainty(
            NRTL_model, NRTL_model_opt, data, variable_name, SSE
        )
        np.testing.assert_array_almost_equal(results.obj, 5.074968578304798)
        np.testing.assert_array_almost_equal(
            np.fromiter(results.theta.values(), dtype=float), [-0.8987624, 1.41048611]
        )
        np.testing.assert_array_almost_equal(results.gradient_f[0], [-0.19649493])
        np.testing.assert_almost_equal(
            results.cov,
            np.array([[0.01194738, -0.02557055], [-0.02557055, 0.05490639]]),
        )
        assert results.propagation_f == pytest.approx(0.0021199499778127204)

    @pytest.mark.integration
    def test_quantify_propagate_uncertainty_NRTL_exception(self):
        """
        It tests an exception error when the ipopt fails for the function quantify_propagate_uncertainty with IDAES NRTL model.
        """
        from idaes.apps.uncertainty_propagation.examples.NRTL_model_scripts import (
            NRTL_model,
            NRTL_model_opt_infeasible,
        )

        variable_name = [
            "fs.properties.tau['benzene','toluene']",
            "fs.properties.tau['toluene','benzene']",
        ]
        current_path = os.path.dirname(os.path.realpath(__file__))
        data = pd.read_csv(os.path.join(current_path, "BT_NRTL_dataset.csv"))

        def SSE(model, data):
            expr = (
                float(data["vap_benzene"])
                - model.fs.flash.vap_outlet.mole_frac_comp[0, "benzene"]
            ) ** 2 + (
                float(data["liq_benzene"])
                - model.fs.flash.liq_outlet.mole_frac_comp[0, "benzene"]
            ) ** 2
            return expr * 1e4

        with pytest.raises(Exception):
            results = quantify_propagate_uncertainty(
                NRTL_model, NRTL_model_opt_infeasible, data, variable_name, SSE
            )

    @pytest.mark.unit
    def test_Exception1(self):
        """
        It tests an ValueError when the tee is not bool for the function quantify_propagate_uncertainty with rooney & biegler's model.
        """
        from idaes.apps.uncertainty_propagation.examples.rooney_biegler import (
            rooney_biegler_model,
            rooney_biegler_model_opt,
        )

        variable_name = ["asymptote", "rate_constant"]
        data = pd.DataFrame(
            data=[[1, 8.3], [2, 10.3], [3, 19.0], [4, 16.0], [5, 15.6], [7, 19.8]],
            columns=["hour", "y"],
        )

        def SSE(model, data):
            expr = sum(
                (data.y[i] - model.response_function[data.hour[i]]) ** 2
                for i in data.index
            )
            return expr

        tee = 1
        with pytest.raises(TypeError):
            results = quantify_propagate_uncertainty(
                rooney_biegler_model,
                rooney_biegler_model_opt,
                data,
                variable_name,
                SSE,
                tee,
            )

    @pytest.mark.unit
    def test_Exception2(self):
        """
        It tests an ValueError when the diagnostic_mode is not bool for the function quantify_propagate_uncertainty with rooney & biegler's model.
        """
        from idaes.apps.uncertainty_propagation.examples.rooney_biegler import (
            rooney_biegler_model,
            rooney_biegler_model_opt,
        )

        variable_name = ["asymptote", "rate_constant"]
        data = pd.DataFrame(
            data=[[1, 8.3], [2, 10.3], [3, 19.0], [4, 16.0], [5, 15.6], [7, 19.8]],
            columns=["hour", "y"],
        )

        def SSE(model, data):
            expr = sum(
                (data.y[i] - model.response_function[data.hour[i]]) ** 2
                for i in data.index
            )
            return expr

        tee = False
        diagnostic_mode = 1
        with pytest.raises(TypeError):
            results = quantify_propagate_uncertainty(
                rooney_biegler_model,
                rooney_biegler_model_opt,
                data,
                variable_name,
                SSE,
                tee,
                diagnostic_mode,
            )

    @pytest.mark.unit
    def test_Exception3(self):
        """
        It tests an ValeError when solver_options is not a dictionary for the function quantify_propagate_uncertainty with rooney & biegler's model.
        """
        from idaes.apps.uncertainty_propagation.examples.rooney_biegler import (
            rooney_biegler_model,
            rooney_biegler_model_opt,
        )

        variable_name = ["asymptote", "rate_constant"]
        data = pd.DataFrame(
            data=[[1, 8.3], [2, 10.3], [3, 19.0], [4, 16.0], [5, 15.6], [7, 19.8]],
            columns=["hour", "y"],
        )

        def SSE(model, data):
            expr = sum(
                (data.y[i] - model.response_function[data.hour[i]]) ** 2
                for i in data.index
            )
            return expr

        tee = False
        diagnostic_mode = False
        solver_options = [1e-8]
        with pytest.raises(TypeError):
            results = quantify_propagate_uncertainty(
                rooney_biegler_model,
                rooney_biegler_model_opt,
                data,
                variable_name,
                SSE,
                tee,
                diagnostic_mode,
                solver_options,
            )

    @pytest.mark.unit
    def test_clean_variable_name1(self):
        """
        It tests the function clean_variable_name when variable names contain ' and spaces.
        """
        theta_names = [
            "fs.properties.tau['benzene', 'toluene']",
            "fs.properties.tau['toluene', 'benzene' ]",
        ]
        theta_names_new, var_dic, clean = clean_variable_name(theta_names)
        theta_names_expected = [
            "fs.properties.tau[benzene,toluene]",
            "fs.properties.tau[toluene,benzene]",
        ]
        assert len(theta_names_expected) == len(theta_names_new)
        assert all([a == b for a, b in zip(theta_names_expected, theta_names_new)])

        assert len(theta_names_expected) == len(var_dic.keys())
        assert all(
            [
                a == b
                for a, b in zip(sorted(theta_names_expected), sorted(var_dic.keys()))
            ]
        )

        assert len(theta_names) == len(var_dic.values())
        assert all(
            [a == b for a, b in zip(sorted(theta_names), sorted(var_dic.values()))]
        )
        assert clean == True

    @pytest.mark.unit
    def test_clean_variable_name2(self):
        """
        It tests the function clean_variable_name when variable names do not contain any ' and spaces.
        """
        theta_names = [
            "fs.properties.tau[benzene,toluene]",
            "fs.properties.tau[toluene,benzene]",
        ]
        theta_names_new, var_dic, clean = clean_variable_name(theta_names)
        assert len(theta_names) == len(theta_names_new)
        assert all([a == b for a, b in zip(theta_names, theta_names_new)])

        assert len(theta_names) == len(var_dic.keys())
        assert all(
            [a == b for a, b in zip(sorted(theta_names), sorted(var_dic.keys()))]
        )

        assert len(theta_names) == len(var_dic.values())
        assert all(
            [a == b for a, b in zip(sorted(theta_names), sorted(var_dic.values()))]
        )
        assert clean == False
