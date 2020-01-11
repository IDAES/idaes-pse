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
from pyomo.contrib.surrogates.polynomial_regression import *
from pyomo.contrib.surrogates import sampling as sp
from pyomo.common.fileutils import PYOMO_ROOT_DIR
import pandas as pd
import numpy as np
import os

datadir = os.path.join(
    PYOMO_ROOT_DIR, 'pyomo', 'contrib', 'surrogates', 'examples', 'data_files')

# Load XY data from high fidelity model from tab file using Pandas. Y
# data must be in the last column.
data = pd.read_csv(os.path.join(datadir, 'three_humpback_data_v4.csv'),
                   header=0, index_col=0)

# data = pd.read_excel('matyas_function.xls', header=0, index_col=0)
# data = pd.read_csv('six_hump_function_data.tab', sep='\t', header=0, index_col=0)
# data = pd.read_csv('cozad_function_data_v2.txt', sep='\s+', header=0, index_col=0)
# data = pd.read_csv('mass_spring_data.txt', sep='\s+', header=None, index_col=None)

b = sp.LatinHypercubeSampling(data, 75)
c = b.sample_points()

# Carry out polynomial regression, feeding in the original data and the
# sampled data
d = PolynomialRegression(data, c, maximum_polynomial_order=10, max_iter=20,
                         multinomials=1, solution_method='pyomo')
p = d.get_feature_vector()
p.pprint()
d.set_additional_terms([ pyo.sin(p['X1']), pyo.cos(p['X1']),
                         pyo.sin(p['X2']), pyo.cos(p['X2']) ])
if True:
    results = d.poly_training()
else:
    results = d.polynomial_regression_fitting(
        [ np.sin(c['X1']), np.cos(c['X1']), np.sin(c['X2']), np.cos(c['X2']) ])

m = pyo.ConcreteModel()
m.x = pyo.Var([1,2])

print("")
print(results.generate_expression([m.x[1], m.x[2]]))

