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
These are tests of the RoundingRegression methodology
Set to solve the Best Subset Selection problem MIP
minimize objective = OLS + complexity_penalty*number_regressors_chosen
Additional hyper-parameters and features will be added over time
"""

from idaes.apps.roundingRegression.RoundingRegression import RoundingRegression
import numpy as np
import pytest


@pytest.mark.unit
def test_simple_regression():
    """Tests the identification of a simple sparse regression"""

    # Create Example Data
    # n = number of data points, p = number of variables, k = true sparsity
    n = 300
    p = 10
    k = 3
    true_variables = np.random.choice([i for i in range(p)], size=k, replace=False)
    true_coefficients = np.zeros(p)
    true_coefficients[true_variables] = 1

    # Predictors
    X = np.random.rand(n, p)

    # Response
    Y = X @ true_coefficients

    # Set Regression penalty (range: 0 - 1)
    complexity_penalty_factor = 0.1

    # Create Instance of Surrogate Modeler
    RR = RoundingRegression(X, Y, complexity_penalty_factor)

    # Build Model
    RR.build_model()

    # Solution
    # MIP Objective
    objective = RR.opt_obj

    # Full Coefficient Array
    coefficients = RR.opt_coeffs

    # Regressors are numbered from 0 to p-1
    chosen_variables = RR.opt_regressors

    # Predict
    Y_predict = X @ coefficients

    # Print
    print("Objective: ", objective)
    print("True Variables: ", true_variables)
    print("Chosen Variables ", chosen_variables)

    # Check if solution makes sense
    assert len(chosen_variables) == np.count_nonzero(
        coefficients
    ), "Number of variables chosen ({}) does not equal number of nonzeros in coefficient vector ({})".format(
        len(chosen_variables), np.count_nonzero(coefficients)
    )

    assert objective >= 0, "Objective ({}) cannot be negative".format(objective)
