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
2d-Griewank function, see Griewank (1981) paper
- 310 samples selected using one of the sampling algorithms and used to train RBF
- Returned RBF tested on over 9500 off-design points and
- Samples for training selected randomly.
- Parity plot and function plots generated.
"""

from idaes.surrogate.pysmo.radial_basis_function import *
from idaes.surrogate.pysmo import sampling as sp
import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d


def main():
    # Load the necessary data
    current_path = os.path.dirname(os.path.realpath(__file__))
    data = pd.read_csv(os.path.join(current_path, 'data_files', 'griewank_data.txt'), sep='\s+', header=None, index_col=None)

    # Scale the data and select 310 sample points via Hammersley sampling
    sd = sp.FeatureScaling()
    data_scaled, data_min, data_max = sd.data_scaling_minmax(data)
    no_training_samples = 310
    b = sp.HammersleySampling(data_scaled, no_training_samples, sampling_type='selection')
    training_data = b.sample_points()

    # Fit an RBF model to 310 selected points
    f1 = RadialBasisFunctions(training_data, basis_function='gaussian', solution_method='pyomo', regularization=True, overwrite=True)
    p = f1.get_feature_vector()
    f1.training()

    # Predict values for other points in loaded data not used in training, evaluate R^2
    y_predicted_pyomo = f1.predict_output(data_scaled[:, :-1])
    r2_pyomo = f1.r2_calculation(data_scaled[:, -1], y_predicted_pyomo)
    print(r2_pyomo)

    # Print RBF expression based on headers of input data
    list_vars = []
    for i in p.keys():
        list_vars.append(p[i])
    print('\nThe RBF expression is: \n', f1.generate_expression(list_vars))

    # 3-D plot of errors between prediction and actual values over 10,200 points
    yy = sp.FeatureScaling.data_unscaling_minmax(y_predicted_pyomo, data_min[0, 2], data_max[0, 2])
    x1 = np.linspace(-20, 20, 101)
    x2 = np.linspace(-20, 20, 101)
    X1, X2 = np.meshgrid(x1, x2)
    Y = yy.reshape(101, 101)
    ax = plt.axes(projection='3d')
    ax.plot_surface(X2, X1, Y, cmap='viridis', edgecolor='none')
    ax.set_xlabel('x1')
    ax.set_ylabel('x2')
    ax.set_zlabel('y')
    plt.show()

    plt.plot(data_scaled[:, 2], data_scaled[:, 2], '-', data_scaled[:, 2], y_predicted_pyomo, 'o')
    plt.show()

if __name__ == "__main__":
    main()
