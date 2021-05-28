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
from idaes.surrogate.pysmo import kriging as krg
from idaes.surrogate.pysmo import sampling as sp
import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d


def main():
    # Load the necessary data
    current_path = os.path.dirname(os.path.realpath(__file__))
    data = pd.read_csv(os.path.join(current_path, 'data_files', 'six_hump_data_2400.txt'), sep='\s+', header=None, index_col=None)

    # Scale the features only - not compulsory as scaling is done at backend anyway.
    sd = sp.FeatureScaling()
    data_scaled_x, data_min, data_max = sd.data_scaling_minmax(data.values[:, :-1])
    y_r = data.values[:, -1]
    data_scaled = np.concatenate((data_scaled_x, y_r.reshape(y_r.shape[0], 1)), axis=1)

    # Select 100 samples for Kriging training
    no_training_samples = 50
    b = sp.HammersleySampling(data_scaled, no_training_samples, 'selection')
    training_data = b.sample_points()

    # Kriging training
    aa = krg.KrigingModel(training_data, numerical_gradients=False, overwrite=True)
    fv = aa.get_feature_vector()
    ab = aa.training()

    # Print Pyomo expression from input variables
    list_vars = []
    for i in fv.keys():
        list_vars.append(fv[i])
    print('The Kriging expression is: \n eq = ', aa.generate_expression(list_vars))

    # Evaluate Kriging model at points not in the training dataset and calculate R^2
    x_pred = data_scaled[:, :-1]
    y_pred = aa.predict_output(x_pred)
    r2 = aa.r2_calculation(data_scaled[:, -1], y_pred)
    print('The R^2 value is: ', r2)

    # 3D error (deviation) plot
    difference_vector = data_scaled[:, 2] - y_pred[:, 0]
    x1 = np.linspace(-3, 3, 61)
    x2 = np.linspace(-2, 2, 41)
    X1, X2 = np.meshgrid(x1, x2, indexing='ij')  # ij indicates matrix arrangement which is what we have
    Y = difference_vector.reshape(61, 41)
    ax = plt.axes(projection='3d')
    ax.plot_surface(X1, X2, Y, cmap='viridis', edgecolor='none')
    # ax.scatter3D(training_data[:, 0], training_data[:, 1], training_data[:, 2]-y_training_pred[:, 0], c='r', marker='^', s=200, depthshade=False)
    ax.set_xlabel('x1')
    ax.set_ylabel('x2')
    ax.set_zlabel('Error')
    plt.show()


if __name__ == "__main__":
    main()
