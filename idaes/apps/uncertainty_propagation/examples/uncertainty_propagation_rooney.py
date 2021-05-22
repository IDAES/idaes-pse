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
from idaes.apps.uncertainty_propagation.uncertainties import quantify_propagate_uncertainty
import pandas as pd
from idaes.apps.uncertainty_propagation.examples.rooney_biegler import rooney_biegler_model,rooney_biegler_model_opt


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


results =  quantify_propagate_uncertainty(rooney_biegler_model,rooney_biegler_model_opt, data, variable_name, SSE)
