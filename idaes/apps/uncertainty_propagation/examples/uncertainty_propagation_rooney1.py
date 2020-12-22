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
from idaes.apps.uncertainty_propagation.uncertainties import quantify_propagate_unucertainty
import pandas as pd
from idaes.apps.uncertainty_propagation.examples.rooney_biegler import rooney_biegler_model
from pyomo.environ import *

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

model_uncertain= ConcreteModel()
model_uncertain.asymptote = Var(initialize = 15)
model_uncertain.rate_constant = Var(initialize = 0.5)
model_uncertain.obj = Objective(expr = model_uncertain.asymptote*( 1 - exp(-model_uncertain.rate_constant*10  )  ), sense=minimize)


obj, theta, cov, propagation_f, propagation_c =  quantify_propagate_unucertainty(rooney_biegler_model,model_uncertain, data, variable_name, SSE)