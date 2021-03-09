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
import sys
import os
sys.path.append(os.path.abspath('..')) # current folder is ~/examples
from idaes.apps.uncertainty_propagation.uncertainties import get_sensitivity
import pandas as pd
from idaes.apps.uncertainty_propagation.examples.rooney_biegler import rooney_biegler_model,rooney_biegler_model_opt
from pyomo.environ import *
import pyomo.contrib.parmest.parmest as parmest


variable_name = ['asymptote', 'rate_constant']
data = pd.DataFrame(data=[[1,8.3],
                          [2,10.3],
                          [3,19.0],
                          [4,16.0],
                          [5,15.6],
                          [7,19.8]],
                    columns=['hour', 'y'])

def SSE(model, data):
    expr = sum((data.y[i] - model.response_function[data.hour[i]])**2 for i in data.index)
    return expr


parmest_class = parmest.Estimator(rooney_biegler_model, data,variable_name,SSE)
obj, theta, cov = parmest_class.theta_est(calc_cov=True)

model_uncertain= ConcreteModel()
model_uncertain.asymptote = Var(initialize = 15)
model_uncertain.rate_constant = Var(initialize = 0.5)
model_uncertain.obj = Objective(expr = model_uncertain.asymptote*( 1 - exp(-model_uncertain.rate_constant*10  )  ), sense=minimize)
theta= {'asymptote': 19.142575284617866, 'rate_constant': 0.53109137696521}
for v in variable_name:
    getattr(model_uncertain, v).setlb(theta[v])
    getattr(model_uncertain, v).setub(theta[v])

gradient_f,gradient_f_dic, gradient_c,gradient_c_dic, line_dic =  get_sensitivity(model_uncertain, variable_name)
print(gradient_f_dic)

print(gradient_c_dic)
