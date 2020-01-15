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
from pyomo.common.fileutils import PYOMO_ROOT_DIR
import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d

# data_dir = os.path.join(
    # PYOMO_ROOT_DIR, 'pyomo', 'contrib', 'surrogates', 'examples', 'data_files')


data = pd.read_csv(os.path.join('idaes/surrogate/pysmo/examples/data_files', 'six_hump_data_2400.txt'),
                   sep='\s+', header=None, index_col=None)
sd = sp.FeatureScaling()
data_scaled_x, data_min, data_max = sd.data_scaling_minmax(data.values[:, :-1])
y_r = data.values[:, -1]
data_scaled = np.concatenate(
    (data_scaled_x, y_r.reshape(y_r.shape[0], 1)), axis=1)
no_training_samples = 100
b = sp.HammersleySampling(data_scaled, no_training_samples)
training_data = b.sample_points()

# Kriging training
aa = krg.KrigingModel(training_data)
fv = aa.get_feature_vector()
ab = aa.kriging_training()
print()

list_vars = []
for i in fv.keys():
    list_vars.append(fv[i])
eq = ab.kriging_generate_expression(list_vars)
print('The Kriging expression is: \n eq = ', ab.kriging_generate_expression(list_vars))

# Kriging testing
x_pred = data_scaled[:, :-1]
y_pred = aa.kriging_predict_output(ab, x_pred)
r2 = aa.r2_calculation(data_scaled[:, -1], y_pred)

difference_vector = data_scaled[:, 2] - y_pred[:, 0]
x1 = np.linspace(-3, 3, 61)
x2 = np.linspace(-2, 2, 41)
X1, X2 = np.meshgrid(x1, x2, indexing='ij')  # ij indicates matrix
                                             # arrangement which is what
                                             # we have
Y = difference_vector.reshape(61, 41)
ax = plt.axes(projection='3d')
ax.plot_surface(X1, X2, Y, cmap='viridis', edgecolor='none')
# ax.scatter3D(training_data[:, 0], training_data[:, 1],
#              training_data[:, 2]-y_training_pred[:, 0],
#              c='r', marker='^', s=200, depthshade=False)
ax.set_xlabel('x1')
ax.set_ylabel('x2')
ax.set_zlabel('Error')
plt.show()
