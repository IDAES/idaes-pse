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
from typing import Tuple, Mapping
import random
import numpy as np
import pyomo.environ as pe
from pyomo.core.base.block import _BlockData
from pyomo.contrib.fbbt import interval
import idaes.apps.flexibility_analysis as flexibility


def create_model() -> Tuple[_BlockData, Mapping, Mapping]:
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

    m.uncertain_temps_set = pe.Set(initialize=[1, 3, 5, 8])
    m.uncertain_temps = pe.Param(
        m.uncertain_temps_set, mutable=True, initialize={1: 620, 3: 388, 5: 583, 8: 313}
    )
    nominal_values = pe.ComponentMap()
    for p in m.uncertain_temps.values():
        nominal_values[p] = p.value

    param_bounds = pe.ComponentMap()
    for p in m.uncertain_temps.values():
        param_bounds[p] = (p.value - 10.0, p.value + 10.0)

    m.variable_temps_set = pe.Set(initialize=[2, 4, 6, 7])
    m.variable_temps = pe.Var(m.variable_temps_set, bounds=(0, 1000))
    m.qc = pe.Var(initialize=0)

    m.balances = pe.Constraint([1, 2, 3, 4])
    m.balances[1] = 1.5 * (m.uncertain_temps[1] - m.variable_temps[2]) == 2 * (
        m.variable_temps[4] - m.uncertain_temps[3]
    )
    m.balances[2] = m.uncertain_temps[5] - m.variable_temps[6] == 2 * (
        563 - m.variable_temps[4]
    )
    m.balances[3] = m.variable_temps[6] - m.variable_temps[7] == 3 * (
        393 - m.uncertain_temps[8]
    )
    m.balances[4] = m.qc == 1.5 * (m.variable_temps[2] - 350)

    m.temp_approaches = pe.Constraint([1, 2, 3, 4])
    m.temp_approaches[1] = m.variable_temps[2] >= m.uncertain_temps[3]
    m.temp_approaches[2] = m.variable_temps[6] >= m.variable_temps[4]
    m.temp_approaches[3] = m.variable_temps[7] >= m.uncertain_temps[8]
    m.temp_approaches[4] = m.variable_temps[6] >= 393

    m.performance = pe.Constraint(expr=m.variable_temps[7] <= 323)

    return m, nominal_values, param_bounds


def get_var_bounds(m):
    """
    Generate a map with valid variable bounds for
    any possible realization of the uncertain parameters
    """
    res = pe.ComponentMap()
    for v in m.variable_temps.values():
        res[v] = (100, 1000)
    res[m.qc] = interval.mul(1.5, 1.5, *interval.sub(100, 1000, 350, 350))
    return res


def main(
    flex_index: bool = False,
    method: flexibility.FlexTestMethod = flexibility.FlexTestMethod.active_constraint,
    plot_history=True,
):
    """
    Run the example

    Parameters
    ----------
    flex_index: bool
        If True, the flexibility index will be solved. Otherwise, the flexibility
        test will be solved.
    method: flexibility.FlexTestMethod
        The method to use for the flexibility test
    plot_history: bool
        Only used if method is flexibility.FlexTestMethod.relu_decision_rule;
        Plots the training history for the neural network
    """
    np.random.seed(0)
    random.seed(1)
    m, nominal_values, param_bounds = create_model()
    var_bounds = get_var_bounds(m)
    config = flexibility.FlexTestConfig()
    config.feasibility_tol = 1e-6
    config.terminate_early = False  # TODO: this does not do anything yet
    config.method = method
    config.minlp_solver = pe.SolverFactory("scip")
    config.sampling_config.solver = pe.SolverFactory("appsi_highs")
    config.sampling_config.strategy = "lhs"
    config.sampling_config.num_points = 600
    config.sampling_config.initialization_strategy = "square"
    if method == flexibility.FlexTestMethod.linear_decision_rule:
        config.decision_rule_config = flexibility.LinearDRConfig()
        config.decision_rule_config.solver = pe.SolverFactory("ipopt")
    elif method == flexibility.FlexTestMethod.relu_decision_rule:
        config.decision_rule_config = flexibility.ReluDRConfig()
        config.decision_rule_config.n_layers = 1
        config.decision_rule_config.n_nodes_per_layer = 15
        config.decision_rule_config.epochs = 3000
        config.decision_rule_config.batch_size = 150
        config.decision_rule_config.scale_inputs = True
        config.decision_rule_config.scale_outputs = True
        config.decision_rule_config.plot_history = plot_history
        config.decision_rule_config.tensorflow_seed = 2
    if not flex_index:
        results = flexibility.solve_flextest(
            m=m,
            uncertain_params=list(nominal_values.keys()),
            param_nominal_values=nominal_values,
            param_bounds=param_bounds,
            controls=[m.qc],
            valid_var_bounds=var_bounds,
            config=config,
        )
        print(results)
    else:
        results = flexibility.solve_flex_index(
            m=m,
            uncertain_params=list(nominal_values.keys()),
            param_nominal_values=nominal_values,
            param_bounds=param_bounds,
            controls=[m.qc],
            valid_var_bounds=var_bounds,
            config=config,
            reconstruct_decision_rule=False,
        )
        print(results)
    return results


if __name__ == "__main__":
    print("\n\n********************Active Constraint**************************")
    main(flex_index=True, method=flexibility.FlexTestMethod.active_constraint)
    print("\n\n********************Linear Decision Rule**************************")
    main(flex_index=True, method=flexibility.FlexTestMethod.linear_decision_rule)
    print("\n\n********************Vertex Enumeration**************************")
    main(flex_index=True, method=flexibility.FlexTestMethod.vertex_enumeration)
    print("\n\n********************ReLU Decision rule**************************")
    main(flex_index=True, method=flexibility.FlexTestMethod.relu_decision_rule)
