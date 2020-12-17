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
from idaes.apps.uncertainty_propagation.uncertainties import quantify_propagate_unucertainty
import pandas as pd
from rooney_biegler import rooney_biegler_model,rooney_biegler_model_opt


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


obj, theta, cov, propagation_f, propagation_c =  quantify_propagate_unucertainty(rooney_biegler_model,rooney_biegler_model_opt, data, variable_name, SSE)

print("\nEstiamted theta")
print(theta)
print("\nCovarinace of theta")
print(cov)
print("\nError propagation in the objective function, var(f(theta))")
print(propagation_f)
print("\nError propagation in the constraints, var(c(theta))")
print(propagation_c)


