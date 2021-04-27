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
# This file applies the ALAMOpython module to the six hump camel problem
# More information on this function can be found at :
#            https://www.sfu.ca/~ssurjano/camel6.html
# This problem utilizes ALAMO's sampling features
from idaes.surrogate import alamopy
import pyomo.environ as pyo

# Import additional python modules for creating the synthetic data
import math
import numpy as np
import pytest


def sixcamel(*x):
    x1, x2 = x
    t1 = np.multiply(
        4.0 - 2.1 * np.power(x1, 2) + np.divide(np.power(x1, 4), 3.0), np.power(x1, 2)
    )
    t2 = np.multiply(4 * np.power(x2, 2) - 4, np.power(x2, 2))
    z = t1 + np.multiply(x1, x2) + t2
    return z

@pytest.mark.unit
def test_multiple_input():

    if not alamopy.has_alamo():
        return False

    # Specify number of points to be used in the training set
    # Validation data can be provided optionally
    ndata = 10
    np.random.seed(2)
    x = np.random.uniform([-2, -1], [2, 1], (ndata, 2))
    z = np.zeros((ndata, 2))
    # specify simulator as examples.sixcamel
    sim = sixcamel
    for i in range(ndata):
        z[i,0] = sim(x[i][0], x[i][1])
        z[i,1] = sim(x[i][0], x[i][1])

    # # Use alamopy's python function wrapper to avoid using ALAMO's I/O format
    almsim = alamopy.wrapwriter(sim)


    # NON GENERIC
    alamo_settings = {'almname': "cam6",
                'monomialpower':(1, 2, 3, 4, 5, 6),
                'multi2power':(1, 2),
                'simulator':almsim,
                'expandoutput':True}

    res = alamopy.doalamo(x, z)

@pytest.mark.unit
def test_gen_pyomo_model():

    if not alamopy.has_alamo():
        return False

    # Specify number of points to be used in the training set
    # Validation data can be provided optionally
    ndata = 20
    np.random.seed(2)
    x = np.random.uniform([-2, -1], [2, 1], (ndata, 2))
    z = np.zeros((ndata, 2))
    # specify simulator as examples.sixcamel
    sim = sixcamel
    for i in range(ndata):
        z[i, 0] = sim(x[i][0], x[i][1])
        z[i, 1] = x[i][0]**2

    # # Use alamopy's python function wrapper to avoid using ALAMO's I/O format
    almsim = alamopy.wrapwriter(sim)

    # NON GENERIC
    alamo_settings = {'almname': "cam6",
                      'monomialpower': (1, 2, 3),
                      'multi2power': (1, 2),
                      'expandoutput': True}

    res = alamopy.doalamo(x, z, **alamo_settings)

    opt = pyo.SolverFactory("ipopt")
    m = pyo.ConcreteModel()
    m.x = pyo.Var(range(1, len(x[0])+1), bounds=(0, 1))
    m.z = pyo.Var()

    pyomo_models = alamopy.generatePyomoExpressions(res, m.x)
    def z_val(model):
        return model.z >= pyomo_models['z2']
    m.valForZ = pyo.Constraint(rule=z_val)

    def min_objective(model):
        return model.z
    m.express_rule = pyo.Objective(rule=min_objective, sense=pyo.minimize)

    pyo_results = opt.solve(m)

    assert pytest.approx(4.4935699883076355e-05) == pyo.value(m.x[1])
    assert pytest.approx(-8.636801831966628e-09) == pyo.value(m.express_rule)

@pytest.mark.unit
def test_single_input_CV():

    if not alamopy.has_alamo():
        return False

    # Specify number of points to be used in the training set
    # Validation data can be provided optionally
    ndata = 10
    x = np.random.uniform([-2, -1], [2, 1], (ndata, 2))
    z = np.zeros((ndata, 1))
    # specify simulator as examples.sixcamel
    sim = sixcamel
    for i in range(ndata):
        z[i,0] = sim(x[i][0], x[i][1])
        # z[i,1] = sim(x[i][0], x[i][1])

    # # Use alamopy's python function wrapper to avoid using ALAMO's I/O format
    almsim = alamopy.wrapwriter(sim)


    # NON GENERIC
    alamo_settings = {'almname': "cam6",
                'monomialpower':(1, 2, 3, 4, 5, 6),
                'multi2power':(1, 2),
                'simulator':almsim,
                'maxiter':20,
                'cvfun':True}

    res = alamopy.doalamo(x, z)
    print(res)
