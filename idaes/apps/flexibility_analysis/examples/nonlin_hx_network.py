#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
A flexibility analysis example from 

    Grossmann, I. E., & Floudas, C. A. (1987). Active constraint strategy for
    flexibility analysis in chemical processes. Computers & Chemical Engineering,
    11(6), 675-693.
"""
from typing import Tuple, MutableMapping
import pyomo.environ as pe
from pyomo.core.base.block import _BlockData
from pyomo.core.base.param import _ParamData
import idaes.apps.flexibility_analysis as flexibility


def create_model() -> Tuple[
    _BlockData,
    MutableMapping[_ParamData, float],
    MutableMapping[_ParamData, Tuple[float, float]],
]:
    """
    This example is from

    Grossmann, I. E., & Floudas, C. A. (1987). Active constraint strategy for
    flexibility analysis in chemical processes. Computers & Chemical Engineering,
    11(6), 675-693.
    """

    print(
        """This example is based off of \n\n
    Grossmann, I. E., & Floudas, C. A. (1987). Active constraint strategy for 
    flexibility analysis in chemical processes. Computers & Chemical Engineering, 
    11(6), 675-693.\n\n"""
    )

    m = pe.ConcreteModel()

    m.qc = pe.Var()
    m.fh1 = pe.Param(mutable=True, initialize=1.4)

    m.f1 = pe.Constraint(expr=-25 + m.qc * ((1 / m.fh1) - 0.5) + 10 / m.fh1 <= 0)
    m.f2 = pe.Constraint(expr=-190 + 10 / m.fh1 + m.qc / m.fh1 <= 0)
    m.f3 = pe.Constraint(expr=-270 + 250 / m.fh1 + m.qc / m.fh1 <= 0)
    m.f4 = pe.Constraint(expr=260 - 250 / m.fh1 - m.qc / m.fh1 <= 0)

    nominal_values = pe.ComponentMap()
    nominal_values[m.fh1] = 1

    param_bounds = pe.ComponentMap()
    param_bounds[m.fh1] = (1, 1.8)

    return m, nominal_values, param_bounds


def get_var_bounds(m):
    """
    Generate a map with valid variable bounds for
    any possible realization of the uncertain parameters
    """
    res = pe.ComponentMap()
    res[m.qc] = (-1000, 1000)
    return res


def main(method):
    """
    Run the example

    Parameters
    ----------
    method: flexibility.FlexTestMethod
        The method to use for the flexibility test
    """
    m, nominal_values, param_bounds = create_model()
    var_bounds = get_var_bounds(m)
    config = flexibility.FlexTestConfig()
    config.feasibility_tol = 1e-6
    config.terminate_early = False
    config.method = method
    config.minlp_solver = pe.SolverFactory("scip")
    config.sampling_config.solver = pe.SolverFactory("scip")
    config.sampling_config.strategy = flexibility.SamplingStrategy.lhs
    config.sampling_config.num_points = 100
    if method == flexibility.FlexTestMethod.linear_decision_rule:
        config.decision_rule_config = flexibility.LinearDRConfig()
        config.decision_rule_config.solver = pe.SolverFactory("ipopt")
    elif method == flexibility.FlexTestMethod.relu_decision_rule:
        config.decision_rule_config = flexibility.ReluDRConfig()
        config.decision_rule_config.n_layers = 1
        config.decision_rule_config.n_nodes_per_layer = 10
        config.decision_rule_config.epochs = 3000
        config.decision_rule_config.batch_size = 50
        config.decision_rule_config.scale_inputs = True
        config.decision_rule_config.scale_outputs = True
    # results = flexibility.solve_flextest(m=m, uncertain_params=list(nominal_values.keys()),
    #                                      param_nominal_values=nominal_values, param_bounds=param_bounds,
    #                                      controls=[m.qc], valid_var_bounds=var_bounds, config=config)
    # print(results)
    results = flexibility.solve_flex_index(
        m=m,
        uncertain_params=list(nominal_values.keys()),
        param_nominal_values=nominal_values,
        param_bounds=param_bounds,
        controls=[m.qc],
        valid_var_bounds=var_bounds,
        config=config,
    )
    print(results)
    return results


if __name__ == "__main__":
    print("\n\n********************Active Constraint**************************")
    main(flexibility.FlexTestMethod.active_constraint)
    print("\n\n********************Linear Decision Rule**************************")
    main(flexibility.FlexTestMethod.linear_decision_rule)
    print("\n\n********************Vertex Enumeration**************************")
    main(flexibility.FlexTestMethod.vertex_enumeration)
    print("\n\n********************ReLU Decision rule**************************")
    main(flexibility.FlexTestMethod.relu_decision_rule)
