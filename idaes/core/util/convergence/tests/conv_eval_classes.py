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
Test classes for the convergence evaluation testing module

Author: Carl Laird
"""
import pyomo.environ as pe
import idaes.core.util.convergence.convergence_base as cb


class ConvEvalFixedVarMutableParam(cb.ConvergenceEvaluation):
    def __init__(self):
        super(ConvEvalFixedVarMutableParam, self).__init__()

    def get_specification(self):
        s = cb.ConvergenceEvaluationSpecification()

        s.add_sampled_input(
            name="var_a", pyomo_path="var_a", lower=0.1, upper=1.9, mean=1.0, std=0.25
        )

        s.add_sampled_input(
            name="param_b", pyomo_path="param_b", lower=50, upper=150, mean=100, std=10
        )
        return s

    def get_initialized_model(self):
        # create a model with two inputs for the convergence evaluation
        # one that is a fixed variable and one that is a mutable param
        m = pe.ConcreteModel()
        m.var_a = pe.Var(initialize=1.0)
        m.var_a.fix(1.0)
        m.param_b = pe.Param(initialize=100, mutable=True)

        m.x = pe.Var(initialize=2.0)
        m.y = pe.Var(initialize=2.0)

        m.obj = pe.Objective(
            expr=(m.var_a - m.x) ** 2 + m.param_b * (m.y - m.x**2) ** 2
        )

        # return the initialized model
        return m

    def get_solver(self):
        opt = pe.SolverFactory("ipopt")
        return opt


class ConvEvalFixedVarImmutableParam(ConvEvalFixedVarMutableParam):
    def __init__(self):
        super(ConvEvalFixedVarImmutableParam, self).__init__()

    def get_initialized_model(self):
        # create a model with two inputs for the convergence evaluation
        # one that is a fixed variable and one that is a mutable param
        m = pe.ConcreteModel()
        m.var_a = pe.Var(initialize=1.0)
        m.var_a.fix(1.0)
        m.param_b = pe.Param(initialize=100)

        m.x = pe.Var(initialize=2.0)
        m.y = pe.Var(initialize=2.0)

        m.obj = pe.Objective(
            expr=(m.var_a - m.x) ** 2 + m.param_b * (m.y - m.x**2) ** 2
        )

        # return the initialized model
        return m


class ConvEvalUnfixedVarMutableParam(ConvEvalFixedVarMutableParam):
    def __init__(self):
        super(ConvEvalUnfixedVarMutableParam, self).__init__()

    def get_initialized_model(self):
        # create a model with two inputs for the convergence evaluation
        # one that is a fixed variable and one that is a mutable param
        m = pe.ConcreteModel()
        m.var_a = pe.Var(initialize=1.0)
        m.param_b = pe.Param(initialize=100, mutable=True)

        m.x = pe.Var(initialize=2.0)
        m.y = pe.Var(initialize=2.0)

        m.obj = pe.Objective(
            expr=(m.var_a - m.x) ** 2 + m.param_b * (m.y - m.x**2) ** 2
        )

        # return the initialized model
        return m
