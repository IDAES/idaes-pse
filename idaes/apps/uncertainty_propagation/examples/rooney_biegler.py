#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021
# by the software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia University
# Research Corporation, et al.  All rights reserved.
#
# Please see the files COPYRIGHT.md and LICENSE.md for full copyright and
# license information.
#################################################################################
"""
Rooney Biegler model, based on Rooney, W. C. and Biegler, L. T. (2001). Design for 
model parameter uncertainty using nonlinear confidence regions. AIChE Journal, 
47(8), 1794-1804.
"""
import pandas as pd
import pyomo.environ as pyo


def rooney_biegler_model(data):
    """This function generates an instance of the rooney & biegler Pyomo model using 'data' as the input argument

    Parameters
    ----------
    data: pandas DataFrame, list of dictionaries, or list of json file names
        Data that is used to build an instance of the Pyomo model

    Returns
    -------
    m: an instance of the Pyomo model
        for estimating parameters and covariance
    """
    model = pyo.ConcreteModel()

    model.asymptote = pyo.Var(initialize=15)
    model.rate_constant = pyo.Var(initialize=0.5)

    def response_rule(m, h):
        expr = m.asymptote * (1 - pyo.exp(-m.rate_constant * h))
        return expr

    model.response_function = pyo.Expression(data.hour, rule=response_rule)

    return model


def rooney_biegler_model_opt():
    """This function generates an instance of the rooney & biegler Pyomo model

    Returns
    -------
    m: an instance of the Pyomo model
        for uncertainty propagation
    """

    model = pyo.ConcreteModel()

    model.asymptote = pyo.Var(initialize=15)
    model.rate_constant = pyo.Var(initialize=0.5)

    model.obj = pyo.Objective(
        expr=model.asymptote * (1 - pyo.exp(-model.rate_constant * 10)),
        sense=pyo.minimize,
    )
    return model
