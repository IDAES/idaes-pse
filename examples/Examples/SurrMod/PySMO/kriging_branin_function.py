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

"""
    Branin function, see Fang paper
        - 81 uniform samples and used to train Kriging model
        - Returned RBF tested on over 5625 off-design points and
        - Parity plot and function plots generated.
"""

from idaes.surrogate.pysmo import kriging as krg
from idaes.surrogate.pysmo import sampling as sp
import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
import pyomo.environ as pyo


def main():
    # Load the necessary data
    current_path = os.path.dirname(os.path.realpath(__file__))
    test_data = pd.read_csv(os.path.join(current_path, 'data_files', 'branin_5625.txt'), sep='\s+', header=None, index_col=None)

    # Select 81 points via uniform sampling on a 9 x 9 grid
    tr_data = sp.UniformSampling(test_data, [9, 9], 'selection')
    training_data = tr_data.sample_points()

    # Fit a Kriging model with regularization to 49 Branin points generated uniformly
    f1 = krg.KrigingModel(training_data, numerical_gradients=False, regularization=True, overwrite=True)
    results_pyomo = f1.training()
    ypred = f1.predict_output(test_data.values[:, :-1])
    r2_pred = f1.r2_calculation(test_data.values[:, -1], ypred)

    # Minimize Kriging surrogate
    m = pyo.ConcreteModel()
    i = pyo.Set(initialize=[1, 2])
    lb = {1: -5, 2: 10}
    ub = {1: 0, 2: 15}

    # bounds for variables
    def f_bounds(m, i):
        return (lb[i], ub[i])

    init_x = {1: 2.5, 2: 7.5}

    # Initial values for x
    def x_init(m, i):
        return init_x[i]

    m.x = pyo.Var(i, bounds=f_bounds, initialize=x_init)
    print('\nThe Kriging expression is: \n', f1.generate_expression([m.x[1], m.x[2]]))
    m.obj = pyo.Objective(expr=f1.generate_expression([m.x[1], m.x[2]]))
    instance = m
    opt = pyo.SolverFactory("ipopt")
    result = opt.solve(instance, tee=True)
    print(pyo.value(instance.obj))
    for p in instance.x:
        print(instance.x[p].value)


if __name__ == "__main__":
    main()




