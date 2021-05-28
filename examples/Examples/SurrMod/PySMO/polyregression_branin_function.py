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
import pyomo.environ as pyo
from idaes.surrogate.pysmo.polynomial_regression import *
from idaes.surrogate.pysmo import kriging
import pandas as pd
import numpy as np
import os


def main():
    def branin_function(x1, x2):
        pi = 3.1417
        t1 = (x2 - (5.1 * x1 * x1 / (4 * pi * pi)) + (5 * x1 / pi) - 6) ** 2
        t2 = (10 * (1 - 1/(8 * pi))*np.cos(x1))
        y = t1 + t2 + 10
        return y

    # Create 20 random samples for training, 100 for testing
    np.random.seed(100)
    ndata = 25
    nval = 100
    x = np.random.uniform([-5, 0], [10, 15], (ndata, 2))
    xval = np.random.uniform([-5, 0], [10, 15], (nval, 2))
    y = branin_function(x[:, 0], x[:, 1])
    yval = branin_function(xval[:, 0], xval[:, 1])
    xy_data = np.concatenate((x, y.reshape(y.size, 1)), axis=1)

    # Train polynomial model with basis functions similar to ALAMO example: 4th order mononomials, exponents, and first and second degree interaction terms
    train_obj = PolynomialRegression(xy_data, xy_data, maximum_polynomial_order=4, multinomials=1, training_split=0.8, number_of_crossvalidations=5, overwrite=True)
    p = train_obj.get_feature_vector()
    train_obj.set_additional_terms([p[0] * p[0] * p[1],  p[0] * p[1] * p[1], p[0] * p[0] * p[1] * p[1], pyo.exp(p[0]), pyo.exp(p[1])])
    train_obj.training()

    # Evaluate model performance as R2
    y_predict = train_obj.predict_output(xval)
    r2 = kriging.KrigingModel.r2_calculation(yval, y_predict)
    print('\nThe R^2 value for the polynomial over the 100 off-design points is', r2)

    # Print Pyomo expression
    m = pyo.ConcreteModel()
    m.x = pyo.Var([1, 2])
    print(train_obj.generate_expression([m.x[1], m.x[2]]))


if __name__ == "__main__":
    main()
