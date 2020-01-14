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
from idaes.surrogate.pysmo.polynomial_regression import *
from idaes.surrogate.pysmo.sampling import sampling as sp
from pyomo.common.fileutils import PYOMO_ROOT_DIR
import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt

os.path.join(PYOMO_ROOT_DIR, 'contrib', 'surrogates', 'examples', 'data_files')


# 1. Load data
boiler_data = pd.read_csv('boiler_data_full.csv', header=0, index_col=0)
# -------------------------------------------------

# 2. Data for variable estimation extracted and pre-processed

# When FG temp is the goal
x = boiler_data.iloc[:, 0:15]
y = boiler_data.iloc[:, 21]  #31]
h1_data = pd.concat([x, y], axis=1)

# Scaling
h1_data_scaled, data_minimums, data_maximums = sp.FeatureScaling.data_scaling(h1_data)
# -------------------------------------------------

# 3. Data regression - no sample selection!
d = PolynomialRegression(h1_data_scaled, h1_data_scaled, maximum_polynomial_order=3, max_iter=0, number_of_crossvalidations=10, training_split=0.75, multinomials=1, solution_method='pyomo')
Results = d.polynomial_regression_fitting()
# -------------------------------------------------

# 4. Check results - this is how to compare!
xy_data = h1_data_scaled
# arf = d.user_defined_terms([h1_data_scaled[:, 2] ** 3])
x_vector = d.polygeneration(Results.polynomial_order, xy_data[:, :-1]) # , arf)
y_predictions_scaled = np.matmul(x_vector, Results.optimal_weights_array)

# -------------------------------------------------
# 6. Plot and save results

plt.plot(h1_data_scaled[:, -1], h1_data_scaled[:, -1], '-', h1_data_scaled[:, -1], y_predictions_scaled[:, -1], 'o')
plt.show()

# -------------------------------------------------
