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
    data = pd.read_csv(os.path.join(current_path, 'data_files', 'six_hump_function_data.tab'), sep='\t', header=0, index_col=0)

    # Scale the data and select 200 sampling points via Hammersley sampling
    sd = sp.FeatureScaling()
    data_scaled, data_min, data_max = sd.data_scaling_minmax(data)
    no_training_samples = 150
    b = sp.HammersleySampling(data_scaled, no_training_samples, 'selection')
    training_data = b.sample_points()

    # Fit an RBF model
    f1 = RadialBasisFunctions(training_data, basis_function='gaussian', solution_method='pyomo', regularization=True)
    p = f1.get_feature_vector()
    f1.training()

    # Predict values for other points in loaded data not used in training, evaluate R^2
    y_predicted_pyomo = f1.predict_output(data_scaled[:, :-1])
    r2_pyomo = f1.r2_calculation(data_scaled[:, -1], y_predicted_pyomo)
    print('R2 over 10201 off-design points:', r2_pyomo)

    # Print RBF expression based on headers of input data
    list_vars = []
    for i in p.keys():
        list_vars.append(p[i])
    print('\nThe RBF expression is: \n', f1.generate_expression(list_vars))

    # Contour plots for surrogate and actual function
    fig, (ax1, ax2) = plt.subplots(ncols=2)
    img1 = ax1.contourf(np.linspace(0, 1, 101), np.linspace(0, 1, 101), y_predicted_pyomo.reshape(101, 101), levels=100,
                        cmap='tab20')
    ax1.grid()
    ax1.set_title('Surrogate')
    fig.colorbar(img1, ax=ax1)
    img2 = ax2.contourf(np.linspace(0, 1, 101), np.linspace(0, 1, 101), data_scaled[:, 2].reshape(101, 101), levels=100,
                        cmap='tab20')
    ax2.grid()
    ax2.set_title('Real')
    fig.colorbar(img2, ax=ax2)
    plt.show()


if __name__ == "__main__":
    main()
