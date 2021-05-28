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
Cardinal sine function, see Kuo (2015) dissertation
    - 100 uniformly spaced samples used for training
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
    data = pd.read_csv(os.path.join(current_path, 'data_files', 'cardinal_sine.txt'), sep='\s+', header=None, index_col=None)

    # Fit a multiquadric RBF model to 100 points generated from the cardinal sine function
    data_scaled = data.values
    f1 = RadialBasisFunctions(data_scaled, basis_function='mq', solution_method='pyomo', regularization=True)
    f1.training()

    # Predict values for other 2500 off-design in loaded data, evaluate R^2
    data2 = pd.read_csv(os.path.join(current_path, 'data_files', 'cardinal_sine_2500.txt'), sep='\s+', header=None, index_col=None)
    data2_scaled = data2.values
    y_predicted = f1.predict_output(data2_scaled[:, :-1])
    r2_pyomo = f1.r2_calculation(data2_scaled[:, -1], y_predicted)

    # Plots
    difference_vector = data2_scaled[:, 2] - y_predicted[:, 0]
    x1 = np.linspace(-1, 1, 50)
    x2 = np.linspace(-1, 1, 50)
    X1, X2 = np.meshgrid(x1, x2)
    Y = difference_vector.reshape(50, 50)
    ax = plt.axes(projection='3d')
    ax.plot_surface(X2, X1, Y, cmap='viridis', edgecolor='none')
    ax.set_xlabel('x1')
    ax.set_ylabel('x2')
    ax.set_zlabel('Error')
    plt.show()


if __name__ == "__main__":
    main()
