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
from idaes.surrogate import alamopy
from idaes.surrogate.pysmo import polynomial_regression, radial_basis_function, kriging
from idaes.surrogate.my_surrogate_base import Surrogate
from idaes.surrogate.my_surrogate_base import Metrics

from sympy.parsing.sympy_parser import parse_expr
from sympy import sympify
from sympy import symbols, var
import sympy
import time
import os.path, pickle

from pyomo.core.expr.sympy_tools import PyomoSympyBimap, sympy_available, Sympy2PyomoVisitor, sympy2pyomo_expression
import pyomo.environ as pyo
import numpy as np

try:
    import ujson as json
except:
    import json

from pyomo.common.config import ConfigBlock, ConfigValue, ConfigList


class GeneralSurrogate(Surrogate):
    CONFIG = Surrogate.CONFIG()

    CONFIG.declare('pyomo_vars', ConfigValue(default=None, domain=list))
    CONFIG.declare('linear', ConfigValue(default=True, domain=bool))
    CONFIG.declare('maximum_polynomial_order', ConfigValue(default=None, domain=int))
    CONFIG.declare('multinomials', ConfigValue(default=True, domain=bool))
    CONFIG.declare('logexp', ConfigValue(default=None, domain=bool))
    CONFIG.declare('sincos', ConfigValue(default=None, domain=bool))
    CONFIG.declare('ratio', ConfigValue(default=False, domain=bool))
    CONFIG.declare('alamo_modeler', ConfigValue(default=1, domain=int))
    CONFIG.declare('convpen', ConfigValue(default=0, domain=int))
    CONFIG.declare('basis_function', ConfigValue(default=None, domain=str))
    CONFIG.declare('regularization', ConfigValue(default=None, domain=bool))
    CONFIG.declare('metric', ConfigValue(default=Metrics.SSE, domain=str))
    CONFIG.declare('additional_features_list', ConfigValue(default=None, domain=list))
    CONFIG.declare('overwrite', ConfigValue(default=True, domain=bool))
    CONFIG.declare('fname', ConfigValue(default=None, domain=str))

    CONFIG.declare('alamopy', ConfigValue(default=True, domain=bool))
    CONFIG.declare('alamopy_rbf', ConfigValue(default=False, domain=bool))
    CONFIG.declare('pysmo_kriging', ConfigValue(default=False, domain=bool))
    CONFIG.declare('pysmo_rbf', ConfigValue(default=False, domain=bool))
    CONFIG.declare('pysmo_polyregression', ConfigValue(default=True, domain=bool))


    def __init__(self, **settings):
        super().__init__(**settings)

        self._pysmo_krg_settings = {}
        self._pysmo_rbf_settings = {}
        self._pysmo_pr_settings = {}
        self._alamo_settings = {}

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

        modeler_alamo = Alamopy(**self._alamo_settings)
        modeler_alamo.regressed_data(self._rdata_in, self._rdata_out)
        if self.config['alamopy']:
            has_alamo_flag = alamopy.multos.has_alamo()
            if has_alamo_flag:
                modeler_alamo.build_model()
                self._models.append(modeler_alamo)
            else:
                print("ALAMO is not installed, and wasn't used for the General Surrogate.")

        _alamo_settings_rbf = self._alamo_settings
        if self.config['alamopy_rbf']:
            _alamo_settings_rbf['grbfcns'] = self.config['alamopy_rbf']
            modeler_alamo_rbf = Alamopy(**_alamo_settings_rbf)
            modeler_alamo_rbf.regressed_data(self._rdata_in, self._rdata_out)
            if has_alamo_flag:
                modeler_alamo_rbf.build_model()
                self._models.append(modeler_alamo_rbf)

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

        if self.config['linear']:
            self._alamo_settings['linfcns'] = 1
        if self.config['logexp']:
            self._alamo_settings['logfcns'] = 1
            self._alamo_settings['expfcns'] = 1
        if self.config['sincos']:
            self._alamo_settings['sinfcns'] = 1
            self._alamo_settings['cosfcns'] = 1

        max_order = self.config['maximum_polynomial_order']
        pow_list = np.arange(0.5, max_order + 1, 0.5)
        if max_order >= 1:
            self._alamo_settings['monomialpower'] = pow_list
        if self.config['multinomials']:

            self._alamo_settings['multi2power'] = np.arange(0.5, max_order / 2 + 1, 0.5)
            self._alamo_settings['multi3power'] = np.arange(0.5, max_order / 2 + 1, 0.5)
            if self.config['ratio']:
                self._alamo_settings['ratiopower'] = pow_list
        self._alamo_settings['expandoutput'] = True
        self._alamo_settings['modeler'] = self.config['alamo_modeler']
        self._alamo_settings['convpen'] = self.config['convpen']

    def get_models(self):
        return self._models


class Alamopy(Surrogate):
    CONFIG = Surrogate.CONFIG()

    CONFIG.declare('xlabels', ConfigValue(default=None, domain=list))
    CONFIG.declare('zlabels', ConfigValue(default=None, domain=list))
    CONFIG.declare('xval', ConfigValue(default=None, domain=list))
    CONFIG.declare('zval', ConfigValue(default=None, domain=list))
    CONFIG.declare('xmin', ConfigValue(default=None, domain=list))
    CONFIG.declare('xmax', ConfigValue(default=None, domain=list))
    CONFIG.declare('modeler', ConfigValue(default=1, domain=int))
    CONFIG.declare('linfcns', ConfigValue(default=1, domain=int))
    CONFIG.declare('expfcns', ConfigValue(default=0, domain=int))
    CONFIG.declare('logfcns', ConfigValue(default=0, domain=int))
    CONFIG.declare('sinfcns', ConfigValue(default=0, domain=int))
    CONFIG.declare('cosfcns', ConfigValue(default=0, domain=int))
    CONFIG.declare('monomialpower', ConfigValue(default=None, domain=list))
    CONFIG.declare('multi2power', ConfigValue(default=None, domain=list))
    CONFIG.declare('multi3power', ConfigValue(default=None, domain=list))
    CONFIG.declare('ratiopower', ConfigValue(default=None, domain=list))
    CONFIG.declare('grbfcns', ConfigValue(default=0, domain=int))
    CONFIG.declare('convpen', ConfigValue(default=None, domain=float))
    CONFIG.declare('screener', ConfigValue(default=None, domain=int))  # not sure
    CONFIG.declare('almname', ConfigValue(default=None, domain=str))
    CONFIG.declare('showalm', ConfigValue(default=False, domain=bool))
    CONFIG.declare('savescratch', ConfigValue(default=None, domain=str))
    CONFIG.declare('savetrace', ConfigValue(default=None, domain=str))
    CONFIG.declare('expandoutput', ConfigValue(default=True, domain=bool))
    CONFIG.declare('almopt', ConfigValue(default=None, domain=str))
    CONFIG.declare('loo', ConfigValue(default=False, domain=bool))
    CONFIG.declare('lmo', ConfigValue(default=False, domain=bool))  # Check
    CONFIG.declare('maxiter', ConfigValue(default=None, domain=int))
    CONFIG.declare('simulator', ConfigValue(default=None))

    def __init__(self, **settings):
        super().__init__(**settings)
        self.alamopy_results = None

    def build_model(self):
        super().build_model()
        if self._vdata_in is not None and self._vdata_out is not None:
            self.alamopy_results = alamopy.alamo(self._rdata_in, self._rdata_out, xval=self._vdata_in,
                                                 zval=self._vdata_out, **self.config)
        else:
            self.alamopy_results = alamopy.alamo(self._rdata_in, self._rdata_out, **self.config)

        self._res = self.alamopy_results
        self._model = self.alamopy_results['model']

        self.handle_results(self._res)
        self.pkl_info['Run settings'] = self.config.value()
        self.pkl_info['Results'] = self._results
        self.pkl_info['Expression'] = self._model

    def handle_results(self, res):
        # sympy to pyomo converter

        self._results = {}
        if type(self.alamopy_results[Metrics.R2]) is dict:
            label = self.alamopy_results['zlabels'][0]
            self._results[Metrics.RMSE] = self.alamopy_results['rmse'][label]
            self._results[Metrics.SSE] = self.alamopy_results['ssr'][label]
            self._results[Metrics.Time] = self.alamopy_results['totaltime'][label]
            self._results[Metrics.MSE] = float(self.alamopy_results['rmse'][label]) ** 2
            self._results[Metrics.Order] = self.alamopy_results['nbas'][label]
            self._results[Metrics.R2] = self.alamopy_results['R2'][label]
        else:
            self._results[Metrics.RMSE] = self.alamopy_results['rmse']
            self._results[Metrics.SSE] = self.alamopy_results['ssr']
            self._results[Metrics.Time] = self.alamopy_results['totaltime']
            self._results[Metrics.MSE] = float(self.alamopy_results['rmse']) ** 2
            self._results[Metrics.Order] = self.alamopy_results['nbas']
            self._results[Metrics.R2] = self.alamopy_results['R2']

        # Generate pyomo expression
        m = pyo.ConcreteModel()
        m.x = pyo.Var(range(len(self._rdata_in)))

        obj_map = PyomoSympyBimap()
        obj_map.sympy2pyomo = {}
        sympy_locals = {}
        i = 1
        for label in res['xlabels']:
            sympy_locals[label] = sympy.Symbol(label)
            sympy_obj = sympy.Symbol(label)
            obj_map.sympy2pyomo[sympy_obj] = m.x[i]
            i += 1

        model_string = ""
        if type(res['model']) is dict:
            key = list(res['model'].keys())[0]
            model_string = res['model'][key].split('=')[1]
        else:
            model_string = res['model'].split('=')[1]
        model_symp = parse_expr(model_string.replace("^", "**"), local_dict=sympy_locals)
        model_pyomo = sympy2pyomo_expression(model_symp, obj_map)
        self._model = model_pyomo

    def run_test_suite(self, return_all=False, allow_grb=True):

        surrogate_modelers = []
        # all
        alamo_settings1 = {'logfcns': 1,
                           'expfcns': 1,
                           'monomialpower': (0.5, 1, 2, 3, -0.5, -1, -2, -3),
                           'multi2power': (0.5, 1, 2, 3, -0.5, -1, -2, -3),
                           'multi3power': (0.5, 1, 2, 3, -0.5, -1, -2, -3),
                           'ratiopower': (0.5, 1, 2, 3, -0.5, -1, -2, -3),
                           'sinfcns': 1,
                           'cosfcns': 1
                           }
        surrogate_modelers.append(Alamopy(**alamo_settings1))

        # all0
        alamo_settings2 = {'logfcns': 1,
                           'expfcns': 1,
                           'monomialpower': (0.5, 1, 2, 3, -0.5, -1, -2, -3),
                           'multi2power': (0.5, 1, 2, 3, -0.5, -1, -2, -3),
                           'multi3power': (0.5, 1, 2, 3, -0.5, -1, -2, -3),
                           'ratiopower': (0.5, 1, 2, 3, -0.5, -1, -2, -3),
                           }
        surrogate_modelers.append(Alamopy(**alamo_settings2))

        # all1
        if allow_grb:
            alamo_settings3 = {'monomialpower': (0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.3, 1.4, 1.6, 1.8, 2, 2.5,
                                                 3, 3.5, 4, -0.1, -0.2, -0.4, -0.6, -0.8, -1, -1.2, -1.3, -1.4,
                                                 -1.6, -1.8, -2, -2.5, -3, -3.5, -4),
                               'grbfcns': 1
                               }
            surrogate_modelers.append(Alamopy(**alamo_settings3))

        # linear
        alamo_settings4 = {}
        surrogate_modelers.append(Alamopy(**alamo_settings4))

        # logexp
        alamo_settings5 = {'logfcns': 1,
                           'expfcns': 1
                           }
        surrogate_modelers.append(Alamopy(**alamo_settings5))

        # mono
        alamo_settings6 = {'monomialpower': (0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.3, 1.4, 1.6, 1.8, 2)}
        surrogate_modelers.append(Alamopy(**alamo_settings6))

        # mono2
        alamo_settings7 = {'monomialpower': (0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.2, 1.3, 1.4, 1.6, 1.8, 2, 2.5,
                                             3, 3.5, 4, -0.1, -0.2, -0.4, -0.6, -0.8, -1, -1.2, -1.3, -1.4,
                                             -1.6, -1.8, -2, -2.5, -3, -3.5, -4)}
        surrogate_modelers.append(Alamopy(**alamo_settings7))

        # multi2
        alamo_settings8 = {'multi2power': (0.5, 1, 2, 3, -0.5, -1, -2, -3)}
        surrogate_modelers.append(Alamopy(**alamo_settings8))

        # multi3
        alamo_settings9 = {'multi3power': (0.5, 1, 2, 3, -0.5, -1, -2, -3)}
        surrogate_modelers.append(Alamopy(**alamo_settings9))

        # olr
        alamo_settings10 = {'modeler': 6, 'convpen': 0}
        surrogate_modelers.append(Alamopy(**alamo_settings10))

        # pce3
        alamo_settings11 = {'monomialpower': (1, 2, 3),
                            'multi2power': (1, 2)}
        surrogate_modelers.append(Alamopy(**alamo_settings11))

        # quad
        alamo_settings12 = {'monomialpower': (1, 2),
                            'multi2power': (1,)}
        surrogate_modelers.append(Alamopy(**alamo_settings12))

        # ratio
        alamo_settings13 = {'ratiopower': (0.5, 1, 2, 3, -0.5, -1, -2, -3)}
        surrogate_modelers.append(Alamopy(**alamo_settings13))

        # ratio2
        alamo_settings14 = {'ratiopower': (0.1, 0.2, 0.4, 0.6, 0.8, 1, 1.2,
                                           1.3, 1.4, 1.6, 1.8, 2, 2.5, 3, 3.5,
                                           4, -0.1, -0.2, -0.4, -0.6, -0.8, -1,
                                           -1.2, -1.3, -1.4, -1.6, -1.8, -2,
                                           -2.5, -3, -3.5, -4)}
        surrogate_modelers.append(Alamopy(**alamo_settings14))

        # rbf
        if (allow_grb):
            alamo_settings15 = {'grbfcns': 1}
            surrogate_modelers.append(Alamopy(**alamo_settings15))

        # sincos
        alamo_settings16 = {'sinfcns': 1,
                            'cosfcns': 1}
        surrogate_modelers.append(Alamopy(**alamo_settings16))

        best_model = None
        all_metrics = []
        for m in surrogate_modelers:
            m.regressed_data(self._rdata_in, self._rdata_out)
            m.build_model()
            metrics = m.get_results()
            if best_model is None or best_model.get_results()[Metrics.MSE] > metrics[Metrics.MSE]:
                best_model = m

            all_metrics.append(m)

        if return_all:
            return all_metrics

        return best_model


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
