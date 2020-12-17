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
Rooney Biegler model, based on Rooney, W. C. and Biegler, L. T. (2001). Design for 
model parameter uncertainty using nonlinear confidence regions. AIChE Journal, 
47(8), 1794-1804.
"""
import pandas as pd
import pyomo.environ as pyo

def rooney_biegler_model(data):
    
    model = pyo.ConcreteModel()

    model.asymptote = pyo.Var(initialize = 15)
    model.rate_constant = pyo.Var(initialize = 0.5)
    
    def response_rule(m, h):
        expr = m.asymptote * (1 - pyo.exp(-m.rate_constant * h))
        return expr
    model.response_function = pyo.Expression(data.hour, rule = response_rule)
    
    def SSE_rule(m):
        return sum((data.y[i] - m.response_function[data.hour[i]])**2 for i in data.index)
    model.SSE = pyo.Objective(rule = SSE_rule, sense=pyo.minimize)
    
    return model


def rooney_biegler_model_opt(theta, theta_names):


    model = pyo.ConcreteModel()

    model.asymptote = pyo.Var(initialize = 15)
    model.rate_constant = pyo.Var(initialize = 0.5)
    
    model.obj = pyo.Objective(expr = model.asymptote*( 1 - pyo.exp(-model.rate_constant*10  )  ), sense=pyo.minimize)
    
    # fix theta with the estimated values    
    for v in theta_names:
        getattr(model, v).setlb(theta[v])
        getattr(model, v).setub(theta[v])

    return model

