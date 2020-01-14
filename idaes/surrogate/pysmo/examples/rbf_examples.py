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
from idaes.surrogate.pysmo.radial_basis_functions import *
from idaes.surrogate.pysmo import sampling as sp
from pyomo.common.fileutils import PYOMO_ROOT_DIR
import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d

datadir = os.path.join(PYOMO_ROOT_DIR, 'contrib', 'surrogates', 'examples', 'data_files')

# =====================================================================================================================
#
# # 1. Usual functions and solution methods
# """
# Typical test suite used: Three-hump Camel, six-hump Camel, Matyas and Cozad functions
# - Fit functions and generate contour plots
# """
#
data = pd.read_csv(os.path.join(datadir, 'three_humpback_data_v4.csv'), header=0, index_col=0)
# data = pd.read_csv(os.path.join(datadir, 'matyas_function.xls'), header=0, index_col=0)
# data = pd.read_csv(os.path.join(datadir, 'six_hump_function_data.tab'), header=0, index_col=0)
# data = pd.read_csv(os.path.join(datadir, 'cozad_function_data_v2.txt'), header=0, index_col=0)
# data = pd.read_csv(os.path.join(datadir, 'exponential_function_data.csv'), header=0, index_col=0)

sd = sp.FeatureScaling()
data_scaled, data_min, data_max = sd.data_scaling_minmax(data)
no_training_samples = 200
b = sp.HammersleySampling(data_scaled, no_training_samples)
training_data = b.sample_points()

f1 = RadialBasisFunctions(training_data, basis_function='gaussian', solution_method='pyomo', regularization=True)
p = f1.get_feature_vector()
p.pprint()
results_pyomo = f1.rbf_training()

y_predicted_pyomo = f1.rbf_predict_output(results_pyomo, data_scaled[:, :-1])
r2_pyomo = f1.r2_calculation(data_scaled[:, -1], y_predicted_pyomo)

list_vars = []
for i in p.keys():
    list_vars.append(p[i])
print('\nThe RBF expression is: \n',results_pyomo.rbf_generate_expression(list_vars))
# import pyomo.environ as pyo
# m = pyo.ConcreteModel()
# m.x = pyo.Var([1, 2])
# print('\nThe RBF expression is: \n',results_pyomo.rbf_generate_expression([ m.x[1], m.x[2] ]))

# Contour plots for surrogate and actual function
fig, (ax1, ax2) = plt.subplots(ncols=2)
img1 = ax1.contourf(np.linspace(0, 1, 101), np.linspace(0, 1, 101), y_predicted_pyomo.reshape(101, 101), levels=100, cmap='tab20')
ax1.grid()
ax1.set_title('Surrogate')
fig.colorbar(img1, ax=ax1)
img2 = ax2.contourf(np.linspace(0, 1, 101), np.linspace(0, 1, 101), data_scaled[:, 2].reshape(101, 101), levels=100, cmap='tab20')
ax2.grid()
ax2.set_title('Real')
fig.colorbar(img2, ax=ax2)
# plt.show()

# =====================================================================================================================
#
# # 2. 2D-Griewank function with plot
#
# """
# 2d-Griewank function, see Griewank (1981) paper
# - 310 samples selected using one of the sampling algorithms and used to train RBF
# - Returned RBF tested on over 9500 off-design points and
# - Samples for training selected randomly.
# - Parity plot and function plots generated.
# """
#
# data = pd.read_csv(os.path.join(datadir, 'griewank_data.txt'), sep='\s+', header=None, index_col=None)
#
# # First, scale data
# sd = sp.FeatureScaling()
# data_scaled, data_min, data_max = sd.data_scaling_minmax(data)
# no_training_samples = 310
# b = sp.HammersleySampling(data_scaled, no_training_samples)
# training_data = b.sample_points()
# f1 = RadialBasisFunctions(training_data, basis_function='gaussian', solution_method='pyomo', regularization=True)
# p = f1.get_feature_vector()
# p.pprint()
# results_pyomo = f1.rbf_training()
#
# y_predicted_pyomo = f1.rbf_predict_output(results_pyomo, data_scaled.values[:, :-1])
# r2_pyomo = f1.r2_calculation(data_scaled.values[:, -1], y_predicted_pyomo)
#
# list_vars = []
# for i in p.keys():
#     list_vars.append(p[i])
# print('\nThe RBF expression is: \n',results_pyomo.rbf_generate_expression(list_vars))
#
# # Plots
# yy = sp.FeatureScaling.data_unscaling_minmax(y_predicted_pyomo, data_min[0, 2], data_max[0, 2])
# x1 = np.linspace(-20, 20, 101)
# x2 = np.linspace(-20,20, 101)
# X1, X2 = np.meshgrid(x1, x2)
# Y = yy.reshape(101, 101)
# ax = plt.axes(projection='3d')
# ax.plot_surface(X2, X1, Y, cmap='viridis', edgecolor='none')
# ax.set_xlabel('x1')
# ax.set_ylabel('x2')
# ax.set_zlabel('y')
# plt.show()
#
# plt.plot(data_scaled.values[:, 2], data_scaled.values[:, 2], '-', data_scaled.values[:, 2], y_predicted_pyomo, 'o')
# plt.show()
# # =====================================================================================================================

# # 3. Boiler Problem RBF surrogate generation
#
# """
# Boiler problem for IDAES project
# - 200 samples supplied by Miguel for 15 variables and about 15 measured outputs
# - Training done with 150 samples (75%) of data, the remaining used as test data for the RBF fit obtained
# - Samples for training selected randomly.
# """
#
# boiler_data = pd.read_csv('C:/Users/OOAmusat/PycharmProjects/WorkProjects/github_files/data_files/boiler_data_full.csv', header=0, index_col=0)
# x = boiler_data.iloc[:, 0:15]
# y = boiler_data.iloc[:, 21]  #31]
# h1_data = pd.concat([x, y], axis=1)
#
# sd = sp.FeatureScaling()
# h1_data_scaled, data_min, data_max = sd.data_scaling_minmax(h1_data)
# no_training_samples = 150
#
# # np.random.shuffle(h1_data_scaled)
# # training_data = h1_data_scaled[:no_training_samples, :]
# # test_data = h1_data_scaled[no_training_samples:, :]
# h1_data_scaled = h1_data_scaled.sample(frac=1).reset_index(drop=True)
# training_data = h1_data_scaled[:no_training_samples]
# test_data = h1_data_scaled[no_training_samples:]
#
# f1 = RadialBasisFunctions(training_data, basis_function='gaussian', solution_method='pyomo', regularization=True)
# results_pyomo = f1.rbf_training()
# p = f1.get_feature_vector()
# p.pprint()
# y_predicted_pyomo = f1.rbf_predict_output(results_pyomo.weights, test_data.values[:, :-1], results_pyomo.centres, results_pyomo.sigma, results_pyomo.regularization_parameter)
# r2_pyomo = f1.r2_calculation(test_data.values[:, -1], y_predicted_pyomo)
#
# list_vars = []
# for i in p.keys():
#     list_vars.append(p[i])
# print('\nThe RBF expression is: \n',results_pyomo.rbf_generate_expression(list_vars))

# # =====================================================================================================================

# # 4. Cardinal sine function
#
# """
# Cardinal sine function, see Kuo (2015) dissertation
# - 100 uniformly spaced samples used for training
# - Returned RBF tested on over 9500 off-design points and
# - Samples for training selected randomly.
# - Parity plot and function plots generated.
# """
#
# data = pd.read_csv(os.path.join(datadir, 'cardinal_sine.txt'), sep='\s+', header=None, index_col=None)
#
# data_scaled = data.values
# f1 = RadialBasisFunctions(data_scaled, basis_function='mq', solution_method='pyomo', regularization=True)
# results_pyomo = f1.rbf_training()
#
# # Test data
# data = pd.read_csv(os.path.join(datadir, 'cardinal_sine_2500.txt'), sep='\s+', header=None, index_col=None)
# data2_scaled = data2.values
# y_predicted = f1.rbf_predict_output(results_pyomo, data2_scaled[:, :-1])
# r2_pyomo = f1.r2_calculation(data2_scaled[:, -1], y_predicted)
#
# # Plots
# difference_vector = data2_scaled[:, 2] - y_predicted[:, 0]
# x1 = np.linspace(-1, 1, 50)
# x2 = np.linspace(-1, 1, 50)
# X1, X2 = np.meshgrid(x1, x2)
# Y = difference_vector.reshape(50, 50)
# ax = plt.axes(projection='3d')
# ax.plot_surface(X2, X1, Y, cmap='viridis', edgecolor='none')
# ax.set_xlabel('x1')
# ax.set_ylabel('x2')
# ax.set_zlabel('Error')
# plt.show()
