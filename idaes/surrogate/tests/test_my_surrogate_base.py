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
import os

# third-party packages
import pytest
import numpy as np
from idaes.surrogate.main import \
    Pysmo_rbf, Pysmo_kriging, Pysmo_polyregression, Metrics, GeneralSurrogate
from pyomo.environ import Block, Var, ConcreteModel, Objective, Set

from idaes.surrogate.my_surrogate_base import SurrogateObject


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
def test_pysmo_krig(branin_dataset):
    m, x, y = branin_dataset

    pysmo_krg_settings = {'numerical_gradients': True,
                          'regularization': True,
                          'pyomo_vars': [m.x[1], m.x[2]]}
    modeler = Pysmo_kriging(**pysmo_krg_settings)

    modeler.regressed_data(x, y)
    modeler.train_surrogate()

    m.obj = Objective(expr=modeler._surrogate)
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
    modeler.train_surrogate()

    m.obj = Objective(expr=modeler._surrogate)
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
    modeler.train_surrogate()

    m.obj = Objective(expr=modeler._surrogate)
    m.pprint()

    modeler.save_results('results.pickle', overwrite=True)

    check_metrics(modeler.get_results())

    return True


@pytest.mark.integration
def test_general_interface(branin_dataset):
    m, x, y = branin_dataset

    general_settings = {'pysmo_polyregression':True, #default
                        'pysmo_kriging':False, #default
                        'pysmo_rbf':False, #default
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
    modeler.train_surrogate()

    m.obj = Objective(expr=modeler._surrogate)
    m.pprint()

    modeler.save_results('results.pickle', overwrite=True)

    os.remove('results.pickle')
    os.remove('solution.pickle')

    check_metrics(modeler.get_results())

    return True


# Interface consistency
# TOOO: What is this method meant for?
def check_metrics(model_metrics):

    missing_keys = [Metrics.MSE, Metrics.Order]
    keys = [Metrics.RMSE, Metrics.SSE, Metrics.Time, Metrics.R2]

    for k in keys:
        assert (k in model_metrics), "Results Dictionary is missing %s"%(k)
        model_metrics[k]


class TestSurrogateObject():
    @pytest.fixture
    def smo(self):
        smo = SurrogateObject("z1 = x1 + x2", ["x1", "x2"], ["z1"])

        return smo

    @pytest.mark.unit
    def test_construct_variables_singleton(self, smo):
        b = Block(concrete=True)

        smo._construct_variables(b)

        assert isinstance(b.x1, Var)
        assert not b.x1.is_indexed()
        assert b.x1.bounds == (None, None)
        assert isinstance(b.x2, Var)
        assert not b.x2.is_indexed()
        assert b.x2.bounds == (None, None)
        assert isinstance(b.z1, Var)
        assert not b.z1.is_indexed()
        assert b.z1.bounds == (None, None)

    @pytest.mark.unit
    def test_construct_variables_singleton_w_bounds(self, smo):
        smo._input_bounds = {"x1": (0, 1), "x2": (10, 20)}

        b = Block(concrete=True)

        smo._construct_variables(b)

        assert isinstance(b.x1, Var)
        assert not b.x1.is_indexed()
        assert b.x1.bounds == (0, 1)
        assert isinstance(b.x2, Var)
        assert not b.x2.is_indexed()
        assert b.x2.bounds == (10, 20)
        assert isinstance(b.z1, Var)
        assert not b.z1.is_indexed()
        assert b.z1.bounds == (None, None)

    @pytest.mark.unit
    def test_construct_variables_indexed(self, smo):
        b = Block(concrete=True)
        b.s = Set(initialize=[1, 2, 3])

        smo._construct_variables(b, index_set=b.s)

        assert isinstance(b.x1, Var)
        assert b.x1.is_indexed()
        assert isinstance(b.x2, Var)
        assert b.x2.is_indexed()
        assert isinstance(b.z1, Var)
        assert b.z1.is_indexed()

        for i in b.s:
            assert b.x1[i].bounds == (None, None)
            assert b.x2[i].bounds == (None, None)
            assert b.z1[i].bounds == (None, None)

    @pytest.mark.unit
    def test_construct_variables_indexed_w_bounds(self, smo):
        smo._input_bounds = {"x1": (0, 1), "x2": (10, 20)}

        b = Block(concrete=True)
        b.s = Set(initialize=[1, 2, 3])

        smo._construct_variables(b, index_set=b.s)

        assert isinstance(b.x1, Var)
        assert b.x1.is_indexed()
        assert isinstance(b.x2, Var)
        assert b.x2.is_indexed()
        assert isinstance(b.z1, Var)
        assert b.z1.is_indexed()

        for i in b.s:
            assert b.x1[i].bounds == (0, 1)
            assert b.x2[i].bounds == (10, 20)
            assert b.z1[i].bounds == (None, None)
