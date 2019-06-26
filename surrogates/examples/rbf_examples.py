from pyomo.contrib import radial_basis_function as rbf
from pyomo.contrib.sampling import sampling as sp
from pyomo.common.fileutils import PYOMO_ROOT_DIR
import pandas as pd
import numpy as np
import os
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d

os.path.join(PYOMO_ROOT_DIR, 'contrib', 'surrogates', 'examples', 'data_files')



# =====================================================================================================================

# # 1. Usual functions and solution methods
# """
# Typical test suite used: Three-hump Camel, six-hump Camel, Matyas and Cozad functions
# - Fit functions and generate contour plots
# """
#
# data = pd.read_csv('three_humpback_data_v4.csv', header=0, index_col=0)
# # # data = pd.read_excel('matyas_function.xls', header=0, index_col=0)
# # # data = pd.read_csv('six_hump_function_data.tab', sep='\t', header=0, index_col=0)
# # # data = pd.read_csv('cozad_function_data_v2.txt', sep='\s+', header=0, index_col=0)
# # # data = pd.read_csv('exponential_function_data.csv', header=0, index_col=0)
#
# sd = sp.FeatureScaling()
# data_scaled, data_min, data_max = sd.data_scaling_minmax(data)
# no_training_samples = 200
# b = sp.HammersleySampling(data_scaled, no_training_samples)
# training_data = b.hs_sample_points()
# f1 = rbf.RadialBasisFunctions(training_data, basis_function='gaussian', solution_method='pyomo', regularization=True)
# results_pyomo = f1.rbf_training()
# y_predicted_pyomo = f1.rbf_predict_output(results_pyomo.weights, data_scaled[:, :-1], results_pyomo.centres, results_pyomo.sigma, results_pyomo.regularization_parameter)
# r2_pyomo = f1.r2_calculation(data_scaled[:, -1], y_predicted_pyomo)
#
# # Contour plots for surrogat eand actual function
# fig, (ax1, ax2) = plt.subplots(ncols=2)
# img1 = ax1.contourf(np.linspace(0, 1, 101), np.linspace(0, 1, 101), y_predicted_pyomo.reshape(101, 101), levels=100, cmap='tab20')
# ax1.grid()
# ax1.set_title('Surrogate')
# fig.colorbar(img1, ax=ax1)
# img2 = ax2.contourf(np.linspace(0, 1, 101), np.linspace(0, 1, 101), data_scaled[:, 2].reshape(101, 101), levels=100, cmap='tab20')
# ax2.grid()
# ax2.set_title('Real')
# fig.colorbar(img2, ax=ax2)
# plt.show()

# =====================================================================================================================

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
# data = pd.read_csv('griewank_data.txt', sep='\s+', header=None, index_col=None)
#
# # First, scale data
# sd = sp.FeatureScaling()
# data_scaled, data_min, data_max = sd.data_scaling_minmax(data)
# no_training_samples = 310
# b = rbf.HammersleySampling(data_scaled, no_training_samples)
# training_data = b.hs_sample_points()
# f1 = rbf.RadialBasisFunctions(training_data, basis_function='gaussian', solution_method='pyomo', regularization=True)
# results_pyomo = f1.rbf_training()
# y_predicted_pyomo = f1.rbf_predict_output(results_pyomo.weights, data_scaled[:, :-1], results_pyomo.centres, results_pyomo.sigma, results_pyomo.regularization_parameter)
# r2_pyomo = f1.r2_calculation(data_scaled[:, -1], y_predicted_pyomo)
#
# # Plots
# yy = rbf.FeatureScaling.data_unscaling_minmax(y_predicted_pyomo, data_min[0, 2], data_max[0, 2])
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
# plt.plot(data_scaled[:, 2], data_scaled[:, 2], '-', data_scaled[:, 2], y_predicted_pyomo, 'o')
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
# boiler_data = pd.read_csv('boiler_data_full.csv', header=0, index_col=0)
# x = boiler_data.iloc[:, 0:15]
# y = boiler_data.iloc[:, 21]  #31]
# h1_data = pd.concat([x, y], axis=1)
#
# sd = sp.FeatureScaling()
# h1_data_scaled, data_min, data_max = sd.data_scaling_minmax(h1_data)
# no_training_samples = 150
#
# np.random.shuffle(h1_data_scaled)
# training_data = h1_data_scaled[:no_training_samples, :]
# test_data = h1_data_scaled[no_training_samples:, :]
#
# f1 = rbf.RadialBasisFunctions(training_data, basis_function='gaussian', solution_method='pyomo', regularization=True)
# results_pyomo = f1.rbf_training()
# y_predicted_pyomo = f1.rbf_predict_output(results_pyomo.weights, test_data[:, :-1], results_pyomo.centres, results_pyomo.sigma, results_pyomo.regularization_parameter)
# r2_pyomo = f1.r2_calculation(test_data[:, -1], y_predicted_pyomo)


# # =====================================================================================================================
#
# # 4. Mass spring problem
#
# """
# Mass-spring problem from Loeven et al. paper (Section 3.1)
# - Data generated from complex combination of sin functions.
# - Range of t= 0-30. Range of omega_n = 0.8 - 1.2
# - Sampling had to be modified - each omega section is sampled separately in order to generate full training data for RBFs
# """
# data = pd.read_csv('mass_spring_data.txt', sep='\s+', header=None, index_col=None)
#
# # First, scale data
# sd = rbf.FeatureScaling()
#
# data_scaled, data_min, data_max = sd.data_scaling_minmax(data)
#
# no_training_samples = 75
#
# p = [0, 1, 2, 3, 4]
# for i in p:
#     q = (i * 60) + i
#     dataset = data_scaled[q:q+61, :]
#     b = sp.LatinHypercubeSampling(dataset, int(no_training_samples/len(p)))
#     td_current = b.lh_sample_points()
#     if i == 0:
#         training_data = td_current
#     else:
#         training_data = np.concatenate((training_data, td_current), axis=0)
#
#
# # # Carry out RBF training
# f1 = rbf.RadialBasisFunctions(training_data, basis_function='gaussian', solution_method='pyomo', regularization=True)
# results_pyomo = f1.rbf_training()
#
# f2 = rbf.RadialBasisFunctions(training_data, basis_function='gaussian', solution_method='bfgs', regularization=True)
# results_bfgs = f2.rbf_training()
#
# f3 = rbf.RadialBasisFunctions(training_data, basis_function='gaussian', solution_method='algebraic', regularization=True)
# results_algebraic = f3.rbf_training()
# #
# y_predicted_pyomo = f1.rbf_predict_output(results_pyomo.weights, data_scaled[:, :-1], results_pyomo.centres, results_pyomo.sigma, results_pyomo.regularization_parameter)
# r2_pyomo = f1.r2_calculation(data_scaled[:, -1], y_predicted_pyomo)
# y_predicted_bfgs = f2.rbf_predict_output(results_bfgs.weights, data_scaled[:, :-1], results_bfgs.centres, results_bfgs.sigma, results_bfgs.regularization_parameter)
# r2_bfgs = f1.r2_calculation(data_scaled[:, -1], y_predicted_bfgs)
# y_predicted_algebraic = f3.rbf_predict_output(results_algebraic.weights, data_scaled[:, :-1], results_algebraic.centres, results_algebraic.sigma, results_algebraic.regularization_parameter)
# r2_algebraic = f1.r2_calculation(data_scaled[:, -1], y_predicted_algebraic)
# print('\nPyomo error:', results_pyomo.rmse, '\nBFGS error:', results_bfgs.rmse, '\nMLE error:', results_algebraic.rmse)
#
#
# # Plot graphs
# tds = training_data[training_data[:, 0].argsort()]
#
# ax1 = plt.subplot(2, 3, 1)
# ax1.plot(tds[0:15, 1], tds[0:15, 2], 'o', data_scaled[0:61, 1], y_predicted_pyomo[0:61, 0], 'r', data_scaled[0:61, 1], data_scaled[0:61, 2], '--')
# ax1.set_title('omega = 0.8')
#
# ax2 = plt.subplot(2, 3, 2)
# ax2.plot(tds[15:30, 1], tds[15:30, 2], 'o', data_scaled[61:122, 1], y_predicted_pyomo[61:122, 0], 'r', data_scaled[61:122, 1], data_scaled[61:122, 2], '--')
# ax2.set_title('omega = 0.9')
#
# ax3 = plt.subplot(2, 3, 3)
# ax3.plot(tds[30:45, 1], tds[30:45, 2], 'o', data_scaled[122:183, 1], y_predicted_pyomo[122:183, 0], 'r', data_scaled[122:183, 1], data_scaled[122:183, 2], '--')
# ax3.set_title('omega = 1')
#
# ax4 = plt.subplot(2, 3, 4)
# ax4.plot(tds[45:60, 1], tds[45:60, 2], 'o', data_scaled[183:244, 1], y_predicted_pyomo[183:244, 0], 'r', data_scaled[183:244, 1], data_scaled[183:244, 2], '--')
# ax4.set_title('omega = 1.1')
#
# ax5 = plt.subplot(2, 3, 5)
# ax5.plot(tds[60:, 1], tds[60:, 2], 'o', data_scaled[244:, 1], y_predicted_pyomo[244:, 0], 'r', data_scaled[244:, 1], data_scaled[244:, 2], '--')
# ax5.set_title('omega = 1.2')
#
# plt.show()

#=====================================================================================================================

# 7. Cardinal sine function

# """
# Cardinal sine function, see Kuo (2015) dissertation
# - 100 uniformly spaced samples used for training
# - Returned RBF tested on over 9500 off-design points and
# - Samples for training selected randomly.
# - Parity plot and function plots generated.
# """
#
# data = pd.read_csv('cardinal_sine.txt', sep='\s+', header=None, index_col=None)
#
# data_scaled = data.values
# f1 = rbf.RadialBasisFunctions(data_scaled, basis_function='mq', solution_method='pyomo', regularization=True)
# results_pyomo = f1.rbf_training()
#
# # Test data
# data2 = pd.read_csv('cardinal_sine_2500.txt', sep='\s+', header=None, index_col=None)
# data2_scaled = data2.values
# y_predicted = f1.rbf_predict_output(results_pyomo.weights, data2_scaled[:, :-1], results_pyomo.centres, results_pyomo.sigma, results_pyomo.regularization_parameter)
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
