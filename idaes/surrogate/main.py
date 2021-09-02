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
from idaes.surrogate.pysmo import polynomial_regression, radial_basis_function, kriging
from idaes.surrogate.my_surrogate_base import Surrogate
from idaes.surrogate.my_surrogate_base import Metrics

import time
import os.path, pickle
import numpy as np

try:
    import ujson as json
except:
    import json

import pyomo.environ as pyo
from pyomo.common.config import ConfigValue


class GeneralSurrogate(Surrogate):
    CONFIG = Surrogate.CONFIG()

    CONFIG.declare('pyomo_vars', ConfigValue(default=None, domain=list))
    CONFIG.declare('linear', ConfigValue(default=True, domain=bool))
    CONFIG.declare('maximum_polynomial_order', ConfigValue(default=None, domain=int))
    CONFIG.declare('multinomials', ConfigValue(default=True, domain=bool))
    CONFIG.declare('logexp', ConfigValue(default=None, domain=bool))
    CONFIG.declare('sincos', ConfigValue(default=None, domain=bool))
    CONFIG.declare('ratio', ConfigValue(default=False, domain=bool))
    CONFIG.declare('convpen', ConfigValue(default=0, domain=int))
    CONFIG.declare('basis_function', ConfigValue(default=None, domain=str))
    CONFIG.declare('regularization', ConfigValue(default=None, domain=bool))
    CONFIG.declare('metric', ConfigValue(default=Metrics.SSE, domain=str))
    CONFIG.declare('additional_features_list', ConfigValue(default=None, domain=list))
    CONFIG.declare('overwrite', ConfigValue(default=True, domain=bool))
    CONFIG.declare('fname', ConfigValue(default=None, domain=str))

    CONFIG.declare('pysmo_kriging', ConfigValue(default=False, domain=bool))
    CONFIG.declare('pysmo_rbf', ConfigValue(default=False, domain=bool))
    CONFIG.declare('pysmo_polyregression', ConfigValue(default=True, domain=bool))


    def __init__(self, **settings):
        super().__init__(**settings)

        self._pysmo_krg_settings = {}
        self._pysmo_rbf_settings = {}
        self._pysmo_pr_settings = {}

        self._models = []

    def build_model(self):
        super().build_model()

        self.parse_config()
        # models = []

        modeler_krig = Pysmo_kriging(**self._pysmo_krg_settings)
        modeler_krig.regressed_data(self._rdata_in, self._rdata_out)
        if self.config['pysmo_kriging']:
            modeler_krig.build_model()
            self._models.append(modeler_krig)

        modeler_rbf = Pysmo_rbf(**self._pysmo_rbf_settings)
        modeler_rbf.regressed_data(self._rdata_in, self._rdata_out)
        if self.config['pysmo_rbf']:
            modeler_rbf.build_model()
            self._models.append(modeler_rbf)

        modeler_pr = Pysmo_polyregression(**self._pysmo_pr_settings)
        modeler_pr.regressed_data(self._rdata_in, self._rdata_out)
        if self.config['pysmo_polyregression']:
            modeler_pr.build_model()
            self._models.append(modeler_pr)

        best_metric = 10e16
        best_config = {}
        for m in self._models:
            metric = m.get_results()[self.config['metric']]
            # print(m, m.get_results()['SSE'], best_metric, metric)
            if best_metric > float(metric):
                best_metric = metric
                best_config = m.config
                self._results = m.get_results()
                self._model = m.get_model()

        if not self.config['pysmo_rbf']:
            self._models.append(modeler_rbf)

        if not self.config['pysmo_kriging']:
            self._models.append(modeler_krig)

        self.pkl_info['Run settings'] = best_config
        self.pkl_info['Results'] = self._results
        self.pkl_info['Expression'] = self._model

    def parse_config(self):

        self._pysmo_krg_settings = {'numerical_gradients': True,
                                    'regularization': self.config['regularization'],
                                    'pyomo_vars': self.config['pyomo_vars'],
                                    'overwrite': self.config['overwrite']}

        if self.config['basis_function'] is None:
            self.config['basis_function'] = 'gaussian'

        self._pysmo_rbf_settings = {'basis_function': self.config['basis_function'],
                                    'regularization': self.config['regularization'],
                                    'pyomo_vars': self.config['pyomo_vars'],
                                    'overwrite': self.config['overwrite']}

        self._pysmo_pr_settings = {'maximum_polynomial_order': self.config['maximum_polynomial_order'],
                                   'multinomials': self.config['multinomials'],
                                   'pyomo_vars': self.config['pyomo_vars'],
                                   'training_split': 0.9,
                                   'number_of_crossvalidations': 5,
                                   'overwrite': self.config['overwrite']}

        max_order = self.config['maximum_polynomial_order']
        pow_list = np.arange(0.5, max_order + 1, 0.5)

    def get_models(self):
        return self._models


class Pysmo_rbf(Surrogate):
    CONFIG = Surrogate.CONFIG()

    CONFIG.declare('basis_function', ConfigValue(default=None, domain=str))
    CONFIG.declare('solution_method', ConfigValue(default=None, domain=str))
    CONFIG.declare('regularization', ConfigValue(default=None, domain=bool))
    CONFIG.declare('pyomo_vars', ConfigValue(default=None, domain=list))
    CONFIG.declare('overwrite', ConfigValue(default=True, domain=bool))
    CONFIG.declare('fname', ConfigValue(default=None, domain=str))

    def __init__(self, **settings):
        super().__init__(**settings)
        self.pysmo_rbf_results = None

    def build_model(self):
        super().build_model()
        start_time = time.time()
        training_data = np.concatenate((self._rdata_in, self._rdata_out.reshape(self._rdata_out.size, 1)), axis=1)
        self.pyomo_vars = dict(
            (k, self.config[k]) for k in ['pyomo_vars'] if k in self.config)  # Extract variable list
        self.rbf_setup = {k: self.config[k] for k in set(self.config) - set(self.pyomo_vars)}
        prob_init = radial_basis_function.RadialBasisFunctions(training_data, **self.rbf_setup)
        feature_vec = prob_init.get_feature_vector()
        self.pysmo_rbf_results = prob_init.training()
        self._results[Metrics.Time] = time.time() - start_time

        self.handle_results(feature_vec)

        self.pkl_info['Run settings'] = self.config
        self.pkl_info['Results'] = self._results
        self.pkl_info['pysmo result object'] = self.pysmo_rbf_results
        self.pkl_info['Expression'] = self._model
        self.pkl_info['Class initialization'] = prob_init

    def handle_results(self, feature_vec):
        self._results[Metrics.RMSE] = self.pysmo_rbf_results.rmse
        self._results[Metrics.SSE] = self.pysmo_rbf_results.rmse ** 2
        self._results[Metrics.R2] = self.pysmo_rbf_results.R2
        # Generate Pyomo expression
        if self.pyomo_vars:
            self._model = self.pysmo_rbf_results.generate_expression(self.pyomo_vars['pyomo_vars'])
        else:
            list_vars = []
            for i in feature_vec.keys():
                list_vars.append(feature_vec[i])
            # PYLINT-TODO-FIX the method "rbf_generate_expression" does not exist for this class
            # maybe the name of the method should be "generate_expression" like in the self.pyomo_vars case?
            # pylint: disable=no-member
            self._model = self.pysmo_rbf_results.rbf_generate_expression(list_vars)


class Pysmo_kriging(Surrogate):
    CONFIG = Surrogate.CONFIG()

    CONFIG.declare('numerical_gradients', ConfigValue(default=None, domain=bool))
    CONFIG.declare('regularization', ConfigValue(default=None, domain=bool))
    CONFIG.declare('pyomo_vars', ConfigValue(default=None, domain=list))
    CONFIG.declare('overwrite', ConfigValue(default=True, domain=bool))
    CONFIG.declare('fname', ConfigValue(default=None, domain=str))

    def __init__(self, **settings):
        super().__init__(**settings)
        self.pysmo_kriging_results = None

    def build_model(self):
        super().build_model()
        start_time = time.time()
        training_data = np.concatenate((self._rdata_in, self._rdata_out.reshape(self._rdata_out.size, 1)), axis=1)
        self.pyomo_vars = dict(
            (k, self.config[k]) for k in ['pyomo_vars'] if k in self.config)  # Extreact variable list
        self.krg_setup = {k: self.config[k] for k in set(self.config) - set(self.pyomo_vars)}
        prob_init = kriging.KrigingModel(training_data, **self.krg_setup)
        feature_vec = prob_init.get_feature_vector()
        self.pysmo_kriging_results = prob_init.training()
        self._results[Metrics.Time] = time.time() - start_time

        self.handle_results(feature_vec)

        self.pkl_info['pysmo result object'] = self.pysmo_kriging_results
        self.pkl_info['Class initialization'] = prob_init
        self.pkl_info['Run settings'] = self.config
        self.pkl_info['Results'] = self._results
        self.pkl_info['Expression'] = self._model

    def handle_results(self, feature_vec):
        self._results[Metrics.RMSE] = self.pysmo_kriging_results.training_rmse
        self._results[Metrics.SSE] = self.pysmo_kriging_results.training_rmse ** 2
        self._results[Metrics.R2] = self.pysmo_kriging_results.training_R2

        # Generate Pyomo expression
        list_vars = []
        for i in feature_vec.keys():
            list_vars.append(feature_vec[i])
        if self.pyomo_vars:
            self._model = self.pysmo_kriging_results.generate_expression(self.pyomo_vars['pyomo_vars'])
        else:
            list_vars = []
            for i in feature_vec.keys():
                list_vars.append(feature_vec[i])
            # PYLINT-TODO-FIX the method "kriging_generate_expression" does not exist for this class
            # maybe the name of the method should be "generate_expression" like in the self.pyomo_vars case?
            # pylint: disable=no-member
            self._model = self.pysmo_kriging_results.kriging_generate_expression(list_vars)


class Pysmo_polyregression(Surrogate):
    CONFIG = Surrogate.CONFIG()

    CONFIG.declare('number_of_crossvalidations', ConfigValue(default=None, domain=int))
    CONFIG.declare('maximum_polynomial_order', ConfigValue(default=None, domain=int))
    CONFIG.declare('training_split', ConfigValue(default=None, domain=float))
    CONFIG.declare('solution_method', ConfigValue(default=None, domain=str))
    CONFIG.declare('multinomials', ConfigValue(default=None, domain=bool))
    CONFIG.declare('additional_features_list', ConfigValue(default=None, domain=list))
    CONFIG.declare('pyomo_vars', ConfigValue(default=None, domain=list))
    CONFIG.declare('overwrite', ConfigValue(default=True, domain=bool))
    CONFIG.declare('fname', ConfigValue(default=None, domain=str))

    def __init__(self, **settings):
        super().__init__(**settings)
        self.pysmo_kriging_results = None

    def build_model(self):
        super().build_model()
        start_time = time.time()
        training_data = np.concatenate((self._rdata_in, self._rdata_out.reshape(self._rdata_out.size, 1)), axis=1)
        self.pyomo_vars = dict((k, self.config[k]) for k in ['pyomo_vars'] if k in self.config)  # Extract variable list
        self.add_terms = dict((k, self.config[k]) for k in ['additional_features_list'] if k in self.config)
        self.polyreg_setup = {k: self.config[k] for k in set(self.config) - set(self.pyomo_vars) - set(self.add_terms)}
        prob_init = polynomial_regression.PolynomialRegression(training_data, training_data, **self.polyreg_setup)
        feature_vec = prob_init.get_feature_vector()

        if self.config['additional_features_list'] is not None:
            try:
                term_lst = []
                res = [sub.replace('ft', 'feature_vec') for sub in self.add_terms['additional_features_list']]
                for i in res:
                    term_lst.append(eval(i))
                prob_init.set_additional_terms(term_lst)
            except:
                raise ValueError('The additional terms could not be evaluated.')

        self.pysmo_polyregression_results = prob_init.training()
        self._results[Metrics.Time] = time.time() - start_time

        self.handle_results(feature_vec)

        self.pkl_info['pysmo result object'] = self.pysmo_polyregression_results
        self.pkl_info['Class initialization'] = prob_init
        self.pkl_info['Run settings'] = self.config
        self.pkl_info['Results'] = self._results
        self.pkl_info['Expression'] = self._model

    def handle_results(self, feature_vec):
        self._results[Metrics.RMSE] = np.sqrt(self.pysmo_polyregression_results.errors['MSE'])
        self._results[Metrics.SSE] = self.pysmo_polyregression_results.errors['MSE']
        self._results[Metrics.R2] = self.pysmo_polyregression_results.errors['R2']

        # Generate Pyomo expression
        list_vars = []
        for i in feature_vec.keys():
            list_vars.append(feature_vec[i])
        if self.pyomo_vars:
            self._model = self.pysmo_polyregression_results.generate_expression(self.pyomo_vars['pyomo_vars'])
        else:
            list_vars = []
            for i in feature_vec.keys():
                list_vars.append(feature_vec[i])
            self._model = self.pysmo_polyregression_results.generate_expression(list_vars)
