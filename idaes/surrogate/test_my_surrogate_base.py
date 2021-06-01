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
Tests for surrogate base class and functions.
"""
# standard library
import json
import os
from pathlib import Path
import tempfile

# third-party packages
import pytest
import numpy as np
from idaes.surrogate.main import Pysmo_rbf, Pysmo_kriging, Pysmo_polyregression, \
                                Alamopy, Metrics, GeneralSurrogate
from pyomo.environ import Var, ConcreteModel, Objective

import idaes.surrogate.alamopy as alamopy

@pytest.fixture
def branin_dataset():
    def branin_function(x1, x2):
        pi = 3.1417
        t1 = (x2 - (5.1 * x1 * x1 / (4 * pi * pi)) + (5 * x1 / pi) - 6) ** 2
        t2 = (10 * (1 - 1/(8 * pi))*np.cos(x1))
        y = t1 + t2 + 10
        return y

    # Create 25 random samples for training, 100 for testing
    np.random.seed(100)
    ndata = 25
    nval = 100
    x = np.random.uniform([-5, 0], [10, 15], (ndata, 2))
    y = branin_function(x[:, 0], x[:, 1])

    m = ConcreteModel()
    m.x = Var([1, 2])

    return m, x, y

@pytest.mark.unit
def test_alamopy(branin_dataset):
    m, x, y = branin_dataset

    alamo_settings = {'monomialpower':(1, 2, 3, 4, 5, 6),
                      'multi2power':(1, 2),
                      'expandoutput':True}
    modeler = Alamopy(**alamo_settings)

    modeler.regressed_data(x, y)

    has_alamo_flag = alamopy.multos.has_alamo()
    if has_alamo_flag:

        modeler.build_model()

        print(modeler.get_results())
        m.obj = Objective(expr=modeler._model)
        m.pprint()

        modeler.save_results('results.pickle', overwrite=True)

        check_metrics(modeler.get_results())
        # check_model_returns(modeler.get_model())
    else:
        with pytest.raises(alamopy.AlamoError):
            modeler.build_model()

    return True

@pytest.mark.unit
def test_pysmo_krig(branin_dataset):
    m, x, y = branin_dataset

    pysmo_krg_settings = {'numerical_gradients': True,
                          'regularization': True,
                          'pyomo_vars': [m.x[1], m.x[2]]}
    modeler = Pysmo_kriging(**pysmo_krg_settings)

    modeler.regressed_data(x, y)
    modeler.build_model()

    m.obj = Objective(expr=modeler._model)
    m.pprint()

    modeler.save_results('results.pickle', overwrite=True)

    check_metrics(modeler.get_results())

    return True

@pytest.mark.component
def test_pysmo_rbf(branin_dataset):
    m, x, y = branin_dataset

    pysmo_rbf_settings = {'basis_function': 'gaussian',
                          'regularization': True,
                          'pyomo_vars': [m.x[1], m.x[2]]}
    modeler = Pysmo_rbf(**pysmo_rbf_settings)

    modeler.regressed_data(x, y)
    modeler.build_model()

    m.obj = Objective(expr=modeler._model)
    m.pprint()

    modeler.save_results('results.pickle', overwrite=True)

    check_metrics(modeler.get_results())

    return True

@pytest.mark.unit
def test_pysmo_poly(branin_dataset):
    m, x, y = branin_dataset

    pysmo_pr_settings = {'maximum_polynomial_order':4,
                         'multinomials':1,
                         'pyomo_vars': [m.x[1], m.x[2]],
                         'training_split':0.9,
                         'number_of_crossvalidations': 5,
                         'additional_features_list': ['ft[0] * ft[0] * ft[1] * ft[1]', 'pyo.exp(ft[0])', 'pyo.exp(ft[1])']}
    modeler = Pysmo_polyregression(**pysmo_pr_settings)

    modeler.regressed_data(x, y)
    modeler.build_model()

    m.obj = Objective(expr=modeler._model)
    m.pprint()

    modeler.save_results('results.pickle', overwrite=True)

    check_metrics(modeler.get_results())

    return True


@pytest.mark.integration_test
@pytest.mark.integration
def test_general_interface(branin_dataset):
    m, x, y = branin_dataset

    general_settings = {'alamopy':True, # default
                        'pysmo_polyregression':True, #default
                        'pysmo_kriging':False, #default
                        'pysmo_rbf':False, #default
                        'alamopy_rbf': False,
                        'linear':True,
                        # 'ratio': True,
                        'pyomo_vars': [m.x[1], m.x[2]],
                        'additional_features_list': ['ft[0] * ft[0] * ft[1] * ft[1]', 'pyo.exp(ft[0])',
                                                     'pyo.exp(ft[1])'],
                        'maximum_polynomial_order': 6,
                        'multinomials': True,
                        'regularization': True}

    modeler = GeneralSurrogate(**general_settings)

    modeler.regressed_data(x, y)
    modeler.build_model()

    m.obj = Objective(expr=modeler._model)
    m.pprint()

    modeler.save_results('results.pickle', overwrite=True)

    os.remove('results.pickle')
    os.remove('solution.pickle')

    check_metrics(modeler.get_results())

    return True


# Interface consistency
def check_metrics(model_metrics):

    missing_keys = [Metrics.MSE, Metrics.Order]
    keys = [Metrics.RMSE, Metrics.SSE, Metrics.Time, Metrics.R2]

    for k in keys:
        assert (k in model_metrics), "Results Dictionary is missing %s"%(k)
        model_metrics[k]