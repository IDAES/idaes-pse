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
Tests for PySMO's family of SurrogateTrainer (PysmoPolyTrainer, PysmoRBFTrainer and PysmoKrigingTrainer)
"""
import pytest
import numpy as np
import pandas as pd
import os
from math import sin, cos, log, exp

from pathlib import Path
from io import StringIO
import re

import pyomo as pyo
from pyomo.environ import ConcreteModel, Var, Constraint
from pyomo.common.tempfiles import TempfileManager

from idaes.core.surrogate.pysmo import (
    polynomial_regression as pr,
    radial_basis_function as rbf,
    kriging as krg,
)

from idaes.core.surrogate.pysmo_surrogate import (
    PysmoPolyTrainer,
    PysmoRBFTrainer,
    PysmoKrigingTrainer,
    PysmoSurrogate,
    PysmoSurrogateTrainingResult,
    PysmoTrainedSurrogate,
)

from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.surrogate.metrics import compute_fit_metrics


dirpath = Path(__file__).parent.resolve()


# String representation of json output for testing


jstring_poly_1 = (
    '{"model_encoding": '
    '{"z1": {"attr": {"regression_data_columns": ["x1", "x2"], '
    '"multinomials": 1, "additional_term_expressions": [], '
    '"optimal_weights_array": [[-75.26111111111476], [-8.815277777775934], [18.81527777777826], [-2.2556956302821618e-13]], '
    '"final_polynomial_order": 1, '
    '"errors": {"MAE": 3.772981926886132e-13, "MSE": 1.5772926701095834e-25, "R2": 1.0, "Adjusted R2": 1.0},'
    ' "extra_terms_feature_vector": ["IndexedParam[x1]", "IndexedParam[x2]"]}, '
    '"map": {"regression_data_columns": "list", "multinomials": "str", '
    '"additional_term_expressions": "list", "optimal_weights_array": "numpy", '
    '"final_polynomial_order": "str", "errors": "str", "extra_terms_feature_vector": "other"}}}, '
    '"input_labels": ["x1", "x2"], '
    '"output_labels": ["z1"], '
    '"input_bounds": {"x1": [0, 5], "x2": [0, 10]}, '
    '"surrogate_type": "poly"}'
)

jstring_poly_2 = (
    '{"model_encoding": '
    '{"z1": {"attr": {"regression_data_columns": ["x1", "x2"], '
    '"multinomials": 1, "additional_term_expressions": [], '
    '"optimal_weights_array": [[-75.26111111111476], [-8.815277777775934], [18.81527777777826], [-2.2556956302821618e-13]], '
    '"final_polynomial_order": 1, '
    '"errors": {"MAE": 3.772981926886132e-13, "MSE": 1.5772926701095834e-25, "R2": 1.0, "Adjusted R2": 1.0},'
    ' "extra_terms_feature_vector": ["IndexedParam[x1]", "IndexedParam[x2]"]}, '
    '"map": {"regression_data_columns": "list", "multinomials": "str", '
    '"additional_term_expressions": "list", "optimal_weights_array": "numpy", '
    '"final_polynomial_order": "str", "errors": "str", "extra_terms_feature_vector": "other"}}, '
    '"z2": {"attr": {"regression_data_columns": ["x1", "x2"], '
    '"multinomials": 1, "additional_term_expressions": [], '
    '"optimal_weights_array": [[-3.0033074724377813], [0.2491731318906352], [1.7508268681094337], [-6.786238238021269e-15]], '
    '"final_polynomial_order": 1, "errors": {"MAE": 1.1901590823981678e-14, "MSE": 1.5225015470765528e-28, "R2": 1.0, "Adjusted R2": 1.0}, '
    '"extra_terms_feature_vector": ["IndexedParam[x1]", "IndexedParam[x2]"]}, '
    '"map": {"regression_data_columns": "list", "multinomials": "str", '
    '"additional_term_expressions": "list", "optimal_weights_array": "numpy", '
    '"final_polynomial_order": "str", "errors": "str", "extra_terms_feature_vector": "other"}}}, '
    '"input_labels": ["x1", "x2"], '
    '"output_labels": ["z1", "z2"], '
    '"input_bounds": null, '
    '"surrogate_type": "poly"}'
)

jstring_poly_3 = (
    '{"model_encoding": '
    '{"z1": {"attr": {"regression_data_columns": ["x1", "x2"], '
    '"multinomials": 0, "additional_term_expressions": ["log(IndexedParam[x1])", "sin(IndexedParam[x2])"], '
    '"optimal_weights_array": [[-14.290243902439855], [6.4274390243899795], [3.572560975609962], [1.9753643165643098e-13], [-4.4048098502003086e-14]], '
    '"final_polynomial_order": 1, '
    '"errors": {"MAE": 1.4210854715202004e-14, "MSE": 2.8188629679897487e-28, "R2": 1.0},'
    ' "extra_terms_feature_vector": ["IndexedParam[x1]", "IndexedParam[x2]"]}, '
    '"map": {"regression_data_columns": "list", "multinomials": "str", '
    '"additional_term_expressions": "list", "optimal_weights_array": "numpy", '
    '"final_polynomial_order": "str", "errors": "str", "extra_terms_feature_vector": "other"}}, '
    '"z2": {"attr": {"regression_data_columns": ["x1", "x2"], '
    '"multinomials": 0, "additional_term_expressions": ["log(IndexedParam[x1])", "sin(IndexedParam[x2])"], '
    '"optimal_weights_array": [[5.704971042443143], [2.4262427606248815], [-0.42624276060821653], [-5.968545102597034e-11], [6.481176706429892e-12]], '
    '"final_polynomial_order": 1, "errors": {"MAE": 3.869645344896829e-12, "MSE": 7.189162598662876e-23, "R2": 1.0}, '
    '"extra_terms_feature_vector": ["IndexedParam[x1]", "IndexedParam[x2]"]}, '
    '"map": {"regression_data_columns": "list", "multinomials": "str", '
    '"additional_term_expressions": "list", "optimal_weights_array": "numpy", '
    '"final_polynomial_order": "str", "errors": "str", "extra_terms_feature_vector": "other"}}}, '
    '"input_labels": ["x1", "x2"], '
    '"output_labels": ["z1", "z2"], '
    '"input_bounds": null, '
    '"surrogate_type": "poly"}'
)

jstring_poly_4 = (
    '{"model_encoding": '
    '{"z1": {"attr": {"regression_data_columns": ["x1", "x2"], '
    '"multinomials": 0, "additional_term_expressions": ["IndexedParam[x1]/IndexedParam[x2]"], '
    '"optimal_weights_array": [[-110.15000000001504], [-17.53750000000189], [27.537500000006148], [-5.3967136315336006e-11]], '
    '"final_polynomial_order": 1, '
    '"errors": {"MAE": 1.0317080523236656e-12, "MSE": 2.126880072091303e-24, "R2": 1.0},'
    ' "extra_terms_feature_vector": ["IndexedParam[x1]", "IndexedParam[x2]"]}, '
    '"map": {"regression_data_columns": "list", "multinomials": "str", '
    '"additional_term_expressions": "other", "optimal_weights_array": "numpy", '
    '"final_polynomial_order": "str", "errors": "str", "extra_terms_feature_vector": "other"}}, '
    '"z2": {"attr": {"regression_data_columns": ["x1", "x2"], '
    '"multinomials": 0, "additional_term_expressions": ["IndexedParam[x1]/IndexedParam[x2]"], '
    '"optimal_weights_array": [[-12.523574144487087], [-2.1308935361219556], [4.1308935361216435], [3.6347869158959156e-12]], '
    '"final_polynomial_order": 1, "errors": {"MAE": 7.762679388179095e-14, "MSE": 6.506051429719772e-27, "R2": 1.0}, '
    '"extra_terms_feature_vector": ["IndexedParam[x1]", "IndexedParam[x2]"]}, '
    '"map": {"regression_data_columns": "list", "multinomials": "str", '
    '"additional_term_expressions": "other", "optimal_weights_array": "numpy", '
    '"final_polynomial_order": "str", "errors": "str", "extra_terms_feature_vector": "other"}}}, '
    '"input_labels": ["x1", "x2"], '
    '"output_labels": ["z1", "z2"], '
    '"input_bounds": null, '
    '"surrogate_type": "poly"}'
)

jstring_rbf = (
    '{"model_encoding": '
    '{"z1": {"attr": {"x_data_columns": ["x1", "x2"], '
    '"x_data": [[0.0, 0.0], [0.25, 0.25], [0.5, 0.5], [0.75, 0.75], [1.0, 1.0]], '
    '"centres": [[0.0, 0.0], [0.25, 0.25], [0.5, 0.5], [0.75, 0.75], [1.0, 1.0]], '
    '"basis_function": "gaussian", '
    '"weights": [[-69.10791015625], [-319807.1317138672], [959336.2551269531], [-959973.7440185547], [320514.66677856445]], '
    '"sigma": 0.05, "regularization_parameter": 0.0, '
    '"rmse": 0.0005986693684275349, "R2": 0.9999971327598984, '
    '"x_data_min": [[1, 5]], "x_data_max": [[5, 9]], "y_data_min": [10], "y_data_max": [50]}, '
    '"map": {"x_data_columns": "list", "x_data": "numpy", "centres": "numpy", '
    '"basis_function": "str", "weights": "numpy", "sigma": "str", "regularization_parameter": "str", '
    '"rmse": "str", "R2": "str", "x_data_min": "numpy", "x_data_max": "numpy", "y_data_min": "numpy", '
    '"y_data_max": "numpy"}}, '
    '"z2": {"attr": {"x_data_columns": ["x1", "x2"], '
    '"x_data": [[0.0, 0.0], [0.25, 0.25], [0.5, 0.5], [0.75, 0.75], [1.0, 1.0]], '
    '"centres": [[0.0, 0.0], [0.25, 0.25], [0.5, 0.5], [0.75, 0.75], [1.0, 1.0]], '
    '"basis_function": "gaussian", "weights": [[-69.10791015625], [-319807.1317138672], [959336.2551269531], [-959973.7440185547], [320514.66677856445]], '
    '"sigma": 0.05, "regularization_parameter": 0.0, '
    '"rmse": 0.0005986693684275349, "R2": 0.9999971327598984, '
    '"x_data_min": [[1, 5]], "x_data_max": [[5, 9]], "y_data_min": [6], "y_data_max": [14]}, '
    '"map": {"x_data_columns": "list", "x_data": "numpy", "centres": "numpy", '
    '"basis_function": "str", "weights": "numpy", "sigma": "str", "regularization_parameter": "str", '
    '"rmse": "str", "R2": "str", "x_data_min": "numpy", "x_data_max": "numpy", "y_data_min": "numpy", '
    '"y_data_max": "numpy"}}}, '
    '"input_labels": ["x1", "x2"], '
    '"output_labels": ["z1", "z2"], '
    '"input_bounds": {"x1": [0, 5], "x2": [0, 10]}, '
    '"surrogate_type": "rbf"}'
)

jstring_krg = (
    '{"model_encoding": '
    '{"z1": {"attr": {"x_data_columns": ["x1", "x2"], '
    '"x_data": [[1, 5], [2, 6], [3, 7], [4, 8], [5, 9]], "x_data_min": [[1, 5]], "x_data_max": [[5, 9]], '
    '"x_data_scaled": [[0.0, 0.0], [0.25, 0.25], [0.5, 0.5], [0.75, 0.75], [1.0, 1.0]], '
    '"optimal_weights": [0.027452451845611077, 0.0010443446337808024], '
    '"optimal_p": 2, "optimal_mean": [[30.00000000077694]], "optimal_variance": [[6503.3113222215325]], '
    '"regularization_parameter": 1.000000000001e-06, '
    '"optimal_covariance_matrix": [[1.000001, 0.9982205353479938, 0.9929011178300284, 0.9840983398813247, 0.971905407660152], '
    "[0.9982205353479938, 1.000001, 0.9982205353479938, 0.9929011178300284, 0.9840983398813247], "
    "[0.9929011178300284, 0.9982205353479938, 1.000001, 0.9982205353479938, 0.9929011178300284], "
    "[0.9840983398813247, 0.9929011178300284, 0.9982205353479938, 1.000001, 0.9982205353479938], "
    "[0.971905407660152, 0.9840983398813247, 0.9929011178300284, 0.9982205353479938, 1.000001]], "
    '"covariance_matrix_inverse": [[108728.9916945844, -240226.85108007095, 82932.18571364644, 121970.72026795016, -73364.51387189297], '
    "[-240226.85108202277, 589985.9891969847, -341158.67300272395, -130592.8567227173, 121970.72027126199], "
    "[82932.18571952915, -341158.67301448685, 516416.75018761755, -341158.6729826693, 82932.18570353556], "
    "[121970.72026201998, -130592.85670691582, -341158.6729945546, 589985.9891699858, -240226.8510697507], "
    "[-73364.51386989365, 121970.72026527137, 82932.18570954115, -240226.85107176506, 108728.99169106234]], "
    '"optimal_y_mu": [[-20.00000000077694], [-10.00000000077694], [-7.769394017032027e-10], [9.99999999922306], [19.99999999922306]], '
    '"training_R2": 0.9999962956016578, "training_rmse": 0.02721910484270722}, '
    '"map": {"x_data_columns": "list", "x_data": "numpy", "x_data_min": "numpy", "x_data_max": "numpy", '
    '"x_data_scaled": "numpy", "optimal_weights": "numpy", "optimal_p": "str", "optimal_mean": "numpy", '
    '"optimal_variance": "numpy", "regularization_parameter": "str", "optimal_covariance_matrix": "numpy", '
    '"covariance_matrix_inverse": "numpy", "optimal_y_mu": "numpy", "training_R2": "str", "training_rmse": "str"}}, '
    '"z2": {"attr": {"x_data_columns": ["x1", "x2"], '
    '"x_data": [[1, 5], [2, 6], [3, 7], [4, 8], [5, 9]], "x_data_min": [[1, 5]], "x_data_max": [[5, 9]], '
    '"x_data_scaled": [[0.0, 0.0], [0.25, 0.25], [0.5, 0.5], [0.75, 0.75], [1.0, 1.0]], '
    '"optimal_weights": [0.02749666901085125, 0.001000000000000049], '
    '"optimal_p": 2, "optimal_mean": [[9.999999999902883]], "optimal_variance": [[260.13320726701056]], '
    '"regularization_parameter": 1e-06, '
    '"optimal_covariance_matrix": [[1.000001, 0.998220543300601, 0.9929011494709431, 0.9840984104422155, 0.9719055315475238], '
    "[0.998220543300601, 1.000001, 0.998220543300601, 0.9929011494709431, 0.9840984104422155], "
    "[0.9929011494709431, 0.998220543300601, 1.000001, 0.998220543300601, 0.9929011494709431], "
    "[0.9840984104422155, 0.9929011494709431, 0.998220543300601, 1.000001, 0.998220543300601], "
    "[0.9719055315475238, 0.9840984104422155, 0.9929011494709431, 0.998220543300601, 1.000001]], "
    '"covariance_matrix_inverse": [[108729.13455237681, -240227.09704128528, 82932.15558036882, 121970.94143487987, -73364.601633614], '
    "[-240227.0970392892, 589986.4681472526, -341158.6596781079, -130593.32427863385, 121970.94144222786], "
    "[82932.15557448889, -341158.6596663887, 516416.7835787105, -341158.659633822, 82932.15555811858], "
    "[121970.94144067129, -130593.32429416949, -341158.6596220617, 589986.4680877628, -240227.09701875152], "
    "[-73364.60163552182, 121970.94144804058, 82932.15555219717, -240227.09701673465, 108729.13454474375]], "
    '"optimal_y_mu": [[-3.999999999902883], [-1.999999999902883], [9.711698112369049e-11], [2.000000000097117], [4.000000000097117]], '
    '"training_R2": 0.9999962956250228, "training_rmse": 0.005443803800474329}, '
    '"map": {"x_data_columns": "list", "x_data": "numpy", "x_data_min": "numpy", "x_data_max": "numpy", '
    '"x_data_scaled": "numpy", "optimal_weights": "numpy", "optimal_p": "str", "optimal_mean": "numpy", '
    '"optimal_variance": "numpy", "regularization_parameter": "str", "optimal_covariance_matrix": "numpy", '
    '"covariance_matrix_inverse": "numpy", "optimal_y_mu": "numpy", "training_R2": "str", "training_rmse": "str"}}}, '
    '"input_labels": ["x1", "x2"], '
    '"output_labels": ["z1", "z2"], '
    '"input_bounds": {"x1": [0, 5], "x2": [0, 10]}, '
    '"surrogate_type": "kriging"}'
)


class TestSurrogateTrainingResult:
    @pytest.fixture
    def pysmo_output_pr(self):
        data = {
            "x1": [1, 2, 3, 4, 5],
            "x2": [5, 6, 7, 8, 9],
            "z1": [10, 20, 30, 40, 50],
        }
        data = pd.DataFrame(data)

        init_pr = pr.PolynomialRegression(
            data, data, maximum_polynomial_order=1, overwrite=True, multinomials=True
        )
        vars = init_pr.get_feature_vector()
        init_pr.training()

        return init_pr, vars

    @pytest.fixture
    def pysmo_output_rbf(self):
        data = {
            "x1": [1, 2, 3, 4, 5],
            "x2": [5, 6, 7, 8, 9],
            "z1": [10, 20, 30, 40, 50],
        }
        data = pd.DataFrame(data)

        init_rbf = rbf.RadialBasisFunctions(
            data, basis_function="linear", overwrite=True
        )
        vars = init_rbf.get_feature_vector()
        init_rbf.training()

        return init_rbf, vars

    @pytest.fixture
    def pysmo_output_krg(self):
        data = {
            "x1": [1, 2, 3, 4, 5],
            "x2": [5, 6, 7, 8, 9],
            "z1": [10, 20, 30, 40, 50],
        }
        data = pd.DataFrame(data)

        init_krg = krg.KrigingModel(data, numerical_gradients=True, overwrite=True)
        vars = init_krg.get_feature_vector()
        init_krg.training()

        return init_krg, vars

    @pytest.mark.unit
    def test_init(self):
        init_func = PysmoSurrogateTrainingResult()
        assert init_func.metrics == {}
        assert init_func._model == None
        assert init_func.expression_str == ""

    @pytest.mark.unit
    def test_model_poly(self, pysmo_output_pr):
        out1, vars = pysmo_output_pr
        init_func_poly = PysmoSurrogateTrainingResult()
        init_func_poly.model = out1
        assert init_func_poly.expression_str == str(
            out1.generate_expression([vars[i] for i in vars.keys()])
        )
        assert init_func_poly._model is not None
        assert isinstance(init_func_poly._model, pr.PolynomialRegression)
        assert init_func_poly._model == out1

    @pytest.mark.unit
    def test_model_rbf(self, pysmo_output_rbf):
        out2, vars = pysmo_output_rbf
        init_func_rbf = PysmoSurrogateTrainingResult()
        init_func_rbf.model = out2
        assert init_func_rbf.expression_str == str(
            out2.generate_expression([vars[i] for i in vars.keys()])
        )
        assert init_func_rbf._model is not None
        assert isinstance(init_func_rbf._model, rbf.RadialBasisFunctions)
        assert init_func_rbf._model == out2

    @pytest.mark.unit
    def test_model_krg(self, pysmo_output_krg):
        out3, vars = pysmo_output_krg
        init_func_krg = PysmoSurrogateTrainingResult()
        init_func_krg.model = out3
        assert init_func_krg.expression_str == str(
            out3.generate_expression([vars[i] for i in vars.keys()])
        )
        assert init_func_krg._model is not None
        assert isinstance(init_func_krg._model, krg.KrigingModel)
        assert init_func_krg._model == out3


class TestTrainedSurrogate:
    @pytest.fixture
    def pysmo_outputs(self):
        data = {
            "x1": [1, 2, 3, 4, 5],
            "x2": [5, 6, 7, 8, 9],
            "z1": [10, 20, 30, 40, 50],
        }
        data = pd.DataFrame(data)

        init_pr = pr.PolynomialRegression(
            data, data, maximum_polynomial_order=1, overwrite=True, multinomials=True
        )
        vars = init_pr.get_feature_vector()
        init_pr.training()

        init_rbf = rbf.RadialBasisFunctions(
            data, basis_function="linear", overwrite=True
        )
        init_rbf.get_feature_vector()
        init_rbf.training()

        init_krg = krg.KrigingModel(data, numerical_gradients=True, overwrite=True)
        init_krg.get_feature_vector()
        init_krg.training()

        return init_pr, init_rbf, init_krg, vars

    @pytest.mark.unit
    def test_init(self):
        init_func = PysmoTrainedSurrogate()
        assert init_func._data == {}
        assert init_func.model_type == ""
        assert init_func.num_outputs == 0
        assert init_func.output_labels == []
        assert init_func.input_labels == None
        assert init_func.input_bounds == None

        init_func1 = PysmoTrainedSurrogate(model_type="poly")
        assert init_func1._data == {}
        assert init_func1.model_type == "poly"
        assert init_func1.num_outputs == 0
        assert init_func1.output_labels == []
        assert init_func1.input_labels == None
        assert init_func1.input_bounds == None

    @pytest.mark.unit
    def test_add_result(self, pysmo_outputs):
        # These need to be tested this way to made sure ``add_result`` builds out model object propoerly.
        out1, out2, out3, vars = pysmo_outputs

        init_func = PysmoTrainedSurrogate()

        outvar = "z1"
        init_func.add_result(outvar, out1)
        assert init_func.output_labels == ["z1"]
        assert init_func._data[outvar] == out1

        outvar = "z2"
        init_func.add_result(outvar, out2)
        assert init_func.output_labels == ["z1", "z2"]
        assert init_func._data[outvar] == out2

        outvar = "z3"
        init_func.add_result(outvar, out3)
        assert init_func.output_labels == ["z1", "z2", "z3"]
        assert init_func._data[outvar] == out3

    @pytest.mark.unit
    def test_get_result(self, pysmo_outputs):
        out1, out2, out3, vars = pysmo_outputs

        init_func = PysmoTrainedSurrogate()

        outvar = "z1"
        init_func.add_result(outvar, out1)
        outvar = "z2"
        init_func.add_result(outvar, out2)
        outvar = "z3"
        init_func.add_result(outvar, out3)

        for i in range(len(init_func.output_labels)):
            assert init_func.get_result(init_func.output_labels[i]) == pysmo_outputs[i]


class TestPysmoPolyTrainer:
    @pytest.fixture
    def pysmo_poly_trainer(self):

        data = {"x1": [1, 2, 3, 4], "x2": [5, 6, 7, 8], "z1": [10, 20, 30, 40]}
        data = pd.DataFrame(data)
        input_labels = ["x1", "x2"]
        output_labels = ["z1"]
        bnds = {"x1": (0, 5), "x2": (0, 10)}
        poly_trainer = PysmoPolyTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=bnds,
            training_dataframe=data,
        )
        return poly_trainer

    @pytest.mark.unit
    def test_defaults(self, pysmo_poly_trainer):
        # Check all defaults
        assert pysmo_poly_trainer.model_type == "poly"
        assert pysmo_poly_trainer.config.maximum_polynomial_order == None
        assert pysmo_poly_trainer.config.multinomials == False
        assert pysmo_poly_trainer.config.number_of_crossvalidations == 3
        assert pysmo_poly_trainer.config.training_split == 0.8
        assert pysmo_poly_trainer.config.solution_method == None
        assert pysmo_poly_trainer.config.extra_features == None

    @pytest.mark.unit
    def test_set_polynomial_order_righttype(self, pysmo_poly_trainer):
        pysmo_poly_trainer.config.maximum_polynomial_order = 3
        assert pysmo_poly_trainer.config.maximum_polynomial_order == 3

    @pytest.mark.unit
    def test_set_polynomial_order_wrongtype(self, pysmo_poly_trainer):
        with pytest.raises(
            ValueError,
            match="invalid value for configuration 'maximum_polynomial_order'",
        ):
            pysmo_poly_trainer.config.maximum_polynomial_order = 3.1

    @pytest.mark.unit
    def test_set_polynomial_order_wrongbounds(self, pysmo_poly_trainer):
        with pytest.raises(
            ValueError,
            match="invalid value for configuration 'maximum_polynomial_order'",
        ):
            pysmo_poly_trainer.config.maximum_polynomial_order = 0

    @pytest.mark.unit
    def test_set_number_of_crossvalidations_righttype(self, pysmo_poly_trainer):
        pysmo_poly_trainer.config.number_of_crossvalidations = 5
        assert pysmo_poly_trainer.config.number_of_crossvalidations == 5

    @pytest.mark.unit
    def test_set_number_of_crossvalidations_wrongtype(self, pysmo_poly_trainer):
        with pytest.raises(
            ValueError,
            match="invalid value for configuration 'number_of_crossvalidations'",
        ):
            pysmo_poly_trainer.config.number_of_crossvalidations = 3.1

    @pytest.mark.unit
    def test_set_number_of_crossvalidations_wrongbounds(self, pysmo_poly_trainer):
        with pytest.raises(
            ValueError,
            match="invalid value for configuration 'number_of_crossvalidations'",
        ):
            pysmo_poly_trainer.config.number_of_crossvalidations = 0

    @pytest.mark.unit
    def test_set_training_split_righttype(self, pysmo_poly_trainer):
        pysmo_poly_trainer.config.training_split = 0.5
        assert pysmo_poly_trainer.config.training_split == 0.5

    @pytest.mark.unit
    def test_set_training_split_wrongbounds(self, pysmo_poly_trainer):
        with pytest.raises(
            ValueError, match="invalid value for configuration 'training_split'"
        ):
            pysmo_poly_trainer.config.training_split = -0.5

    @pytest.mark.unit
    def test_set_solution_method_righttype_1(self, pysmo_poly_trainer):
        pysmo_poly_trainer.config.solution_method = "mle"
        assert pysmo_poly_trainer.config.solution_method == "mle"

    @pytest.mark.unit
    def test_set_solution_method_righttype_2(self, pysmo_poly_trainer):
        pysmo_poly_trainer.config.solution_method = "pyomo"
        assert pysmo_poly_trainer.config.solution_method == "pyomo"

    @pytest.mark.unit
    def test_set_solution_method_righttype_3(self, pysmo_poly_trainer):
        pysmo_poly_trainer.config.solution_method = "bfgs"
        assert pysmo_poly_trainer.config.solution_method == "bfgs"

    @pytest.mark.unit
    def test_set_solution_method_wrongtype(self, pysmo_poly_trainer):
        with pytest.raises(
            ValueError, match="invalid value for configuration 'solution_method'"
        ):
            pysmo_poly_trainer.config.solution_method = "bfgh"

    @pytest.mark.unit
    def test_set_multinomials_righttype_1(self, pysmo_poly_trainer):
        pysmo_poly_trainer.config.multinomials = True
        assert pysmo_poly_trainer.config.multinomials == True

    @pytest.mark.unit
    def test_set_multinomials_righttype_2(self, pysmo_poly_trainer):
        pysmo_poly_trainer.config.multinomials = False
        assert pysmo_poly_trainer.config.multinomials == False

    @pytest.mark.unit
    def test_set_multinomials_righttype_3(self, pysmo_poly_trainer):
        pysmo_poly_trainer.config.multinomials = "False"
        assert pysmo_poly_trainer.config.multinomials == False

    @pytest.mark.unit
    def test_set_multinomials_righttype_4(self, pysmo_poly_trainer):
        pysmo_poly_trainer.config.multinomials = "True"
        assert pysmo_poly_trainer.config.multinomials == True

    @pytest.mark.unit
    def test_set_multinomials_righttype_5(self, pysmo_poly_trainer):
        pysmo_poly_trainer.config.multinomials = 1
        assert pysmo_poly_trainer.config.multinomials == True

    @pytest.mark.unit
    def test_set_multinomials_righttype_6(self, pysmo_poly_trainer):
        pysmo_poly_trainer.config.multinomials = 0
        assert pysmo_poly_trainer.config.multinomials == False

    @pytest.mark.unit
    def test_set_multinomials_wrongtype(self, pysmo_poly_trainer):
        with pytest.raises(ValueError):
            pysmo_poly_trainer.config.multinomials = 2

    @pytest.mark.unit
    def test_set_extra_features_righttype_2(self, pysmo_poly_trainer):
        pysmo_poly_trainer.config.extra_features = ["x1 / x2"]
        assert pysmo_poly_trainer.config.extra_features == ["x1 / x2"]

    @pytest.mark.unit
    def test_set_extra_features_righttype_2(self, pysmo_poly_trainer):
        pysmo_poly_trainer.config.extra_features = ["x1 / x2", "sin(x1)"]
        assert pysmo_poly_trainer.config.extra_features == ["x1 / x2", "sin(x1)"]

    @pytest.mark.unit
    def test_set_extra_features_wrongtype(self, pysmo_poly_trainer):
        with pytest.raises(NameError):
            pysmo_poly_trainer.config.extra_features = x1 / x2

    @pytest.mark.unit
    def test_set_extra_features_wrongtype(self, pysmo_poly_trainer):
        with pytest.raises(ValueError):
            pysmo_poly_trainer.config.extra_features = 10

    @pytest.mark.unit
    def test_create_model_no_extra_features(self, pysmo_poly_trainer):
        pysmo_poly_trainer.config.multinomials = 1
        pysmo_poly_trainer.config.maximum_polynomial_order = 1
        pysmo_poly_trainer.config.solution_method = "mle"
        pysmo_poly_trainer.config.number_of_crossvalidations = 2
        pysmo_poly_trainer.config.training_split = 0.9

        output_label = "z5"
        data = {"x1": [1, 2, 3, 4], "x2": [5, 6, 7, 8], "z1": [10, 20, 30, 40]}
        data = pd.DataFrame(data)

        model = pysmo_poly_trainer._create_model(data, output_label)

        assert (
            model.max_polynomial_order
            == pysmo_poly_trainer.config.maximum_polynomial_order
        )
        assert model.overwrite == True
        assert model.multinomials == pysmo_poly_trainer.config.multinomials
        assert model.solution_method == "mle"
        assert (
            model.number_of_crossvalidations
            == pysmo_poly_trainer.config.number_of_crossvalidations
        )
        assert model.fraction_training == 0.9
        assert model.filename == "solution.pickle"
        assert model.number_of_x_vars == data.shape[1] - 1
        assert model.additional_term_expressions == []
        assert model.extra_terms_feature_vector == None
        np.testing.assert_array_equal(model.original_data, data.values)
        np.testing.assert_array_equal(model.regression_data, data.values)
        assert model.regression_data_columns == data.columns.tolist()[:-1]
        assert list(model.feature_list._data.keys()) == data.columns.tolist()[:-1]

    @pytest.mark.unit
    def test_create_model_with_extra_features(self, pysmo_poly_trainer):
        pysmo_poly_trainer.config.multinomials = 0
        pysmo_poly_trainer.config.maximum_polynomial_order = 2
        pysmo_poly_trainer.config.solution_method = "mle"
        pysmo_poly_trainer.config.number_of_crossvalidations = 2
        pysmo_poly_trainer.config.training_split = 0.9
        pysmo_poly_trainer.config.extra_features = [
            "sin(x1)/cos(x2)",
            "log(x1)*sin(x2)",
            "x1/x2",
        ]

        output_label = "z1"
        data = {
            "x1": [1, 2, 3, 4, 5, 6, 7, 8],
            "x2": [5, 6, 7, 8, 9, 10, 11, 12],
            "z1": [10, 20, 30, 40, 50, 60, 70, 80],
        }
        data = pd.DataFrame(data)

        model = pysmo_poly_trainer._create_model(data, output_label)

        assert model.overwrite == True
        assert model.multinomials == pysmo_poly_trainer.config.multinomials
        assert model.solution_method == "mle"
        assert (
            model.number_of_crossvalidations
            == pysmo_poly_trainer.config.number_of_crossvalidations
        )
        assert (
            model.max_polynomial_order
            == pysmo_poly_trainer.config.maximum_polynomial_order
        )
        assert model.fraction_training == 0.9
        assert model.filename == "solution.pickle"
        assert model.number_of_x_vars == data.shape[1] - 1
        np.testing.assert_array_equal(model.original_data, data.values)
        np.testing.assert_array_equal(model.regression_data, data.values)
        assert model.regression_data_columns == data.columns.tolist()[:-1]
        assert list(model.feature_list._data.keys()) == data.columns.tolist()[:-1]
        assert len(model.additional_term_expressions) == 3
        assert isinstance(model.additional_term_expressions, list)
        assert isinstance(
            model.additional_term_expressions[0],
            pyo.core.expr.numeric_expr.NPV_DivisionExpression,
        )
        assert isinstance(
            model.additional_term_expressions[1],
            pyo.core.expr.numeric_expr.ProductExpression,
        )
        assert isinstance(
            model.additional_term_expressions[2],
            pyo.core.expr.numeric_expr.NPV_DivisionExpression,
        )
        assert (
            str(model.additional_term_expressions[0])
            == "sin(IndexedParam[x1])/cos(IndexedParam[x2])"
        )
        assert (
            str(model.additional_term_expressions[1])
            == "log(IndexedParam[x1])*sin(IndexedParam[x2])"
        )
        assert (
            str(model.additional_term_expressions[2])
            == "IndexedParam[x1]/IndexedParam[x2]"
        )
        assert model.extra_terms_feature_vector == None


class TestPysmoRBFTrainer:
    @pytest.fixture
    def pysmo_rbf_trainer(self):

        data = {"x1": [1, 2, 3, 4], "x2": [5, 6, 7, 8], "z1": [10, 20, 30, 40]}
        data = pd.DataFrame(data)
        input_labels = ["x1", "x2"]
        output_labels = ["z1"]
        bnds = {"x1": (0, 5), "x2": (0, 10)}
        rbf_trainer = PysmoRBFTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=bnds,
            training_dataframe=data,
        )
        return rbf_trainer

    @pytest.mark.unit
    def test_defaults(self, pysmo_rbf_trainer):
        # Check all defaults
        assert pysmo_rbf_trainer.model_type == "None rbf"
        assert pysmo_rbf_trainer.config.basis_function == None
        assert pysmo_rbf_trainer.config.regularization == None
        assert pysmo_rbf_trainer.config.solution_method == None

    @pytest.mark.unit
    def test_set_basis_function_righttype_1(self, pysmo_rbf_trainer):
        pysmo_rbf_trainer.config.basis_function = "linear"
        assert pysmo_rbf_trainer.config.basis_function == "linear"

    @pytest.mark.unit
    def test_set_basis_function_righttype_2(self, pysmo_rbf_trainer):
        pysmo_rbf_trainer.config.basis_function = "cubic"
        assert pysmo_rbf_trainer.config.basis_function == "cubic"

    @pytest.mark.unit
    def test_set_basis_function_righttype_3(self, pysmo_rbf_trainer):
        pysmo_rbf_trainer.config.basis_function = "imq"
        assert pysmo_rbf_trainer.config.basis_function == "imq"

    @pytest.mark.unit
    def test_set_basis_function_righttype_4(self, pysmo_rbf_trainer):
        pysmo_rbf_trainer.config.basis_function = "mq"
        assert pysmo_rbf_trainer.config.basis_function == "mq"

    @pytest.mark.unit
    def test_set_basis_function_righttype_5(self, pysmo_rbf_trainer):
        pysmo_rbf_trainer.config.basis_function = "gaussian"
        assert pysmo_rbf_trainer.config.basis_function == "gaussian"

    @pytest.mark.unit
    def test_set_basis_function_righttype_6(self, pysmo_rbf_trainer):
        pysmo_rbf_trainer.config.basis_function = "spline"
        assert pysmo_rbf_trainer.config.basis_function == "spline"

    @pytest.mark.unit
    def test_set_basis_function_outdomain(self, pysmo_rbf_trainer):
        with pytest.raises(
            ValueError, match="invalid value for configuration 'basis_function'"
        ):
            pysmo_rbf_trainer.config.basis_function = "mqimq"

    @pytest.mark.unit
    def test_set_solution_method_righttype_1(self, pysmo_rbf_trainer):
        pysmo_rbf_trainer.config.solution_method = "algebraic"
        assert pysmo_rbf_trainer.config.solution_method == "algebraic"

    @pytest.mark.unit
    def test_set_solution_method_righttype_2(self, pysmo_rbf_trainer):
        pysmo_rbf_trainer.config.solution_method = "pyomo"
        assert pysmo_rbf_trainer.config.solution_method == "pyomo"

    @pytest.mark.unit
    def test_set_solution_method_righttype_3(self, pysmo_rbf_trainer):
        pysmo_rbf_trainer.config.solution_method = "bfgs"
        assert pysmo_rbf_trainer.config.solution_method == "bfgs"

    @pytest.mark.unit
    def test_set_solution_method_wrongtype(self, pysmo_rbf_trainer):
        with pytest.raises(
            ValueError, match="invalid value for configuration 'solution_method'"
        ):
            pysmo_rbf_trainer.config.solution_method = "mle"

    @pytest.mark.unit
    def test_set_regularization_righttype_1(self, pysmo_rbf_trainer):
        pysmo_rbf_trainer.config.regularization = True
        assert pysmo_rbf_trainer.config.regularization == True

    @pytest.mark.unit
    def test_set_regularization_righttype_2(self, pysmo_rbf_trainer):
        pysmo_rbf_trainer.config.regularization = False
        assert pysmo_rbf_trainer.config.regularization == False

    @pytest.mark.unit
    def test_set_regularization_righttype_3(self, pysmo_rbf_trainer):
        pysmo_rbf_trainer.config.regularization = "False"
        assert pysmo_rbf_trainer.config.regularization == False

    @pytest.mark.unit
    def test_set_regularization_righttype_4(self, pysmo_rbf_trainer):
        pysmo_rbf_trainer.config.regularization = "True"
        assert pysmo_rbf_trainer.config.regularization == True

    @pytest.mark.unit
    def test_set_regularization_righttype_5(self, pysmo_rbf_trainer):
        pysmo_rbf_trainer.config.regularization = 1
        assert pysmo_rbf_trainer.config.regularization == True

    @pytest.mark.unit
    def test_set_regularization_righttype_6(self, pysmo_rbf_trainer):
        pysmo_rbf_trainer.config.regularization = 0
        assert pysmo_rbf_trainer.config.regularization == False

    @pytest.mark.unit
    def test_set_regularization_wrongtype(self, pysmo_rbf_trainer):
        with pytest.raises(ValueError):
            pysmo_rbf_trainer.config.regularization = 2

    @pytest.mark.unit
    def test_create_model_defaults(self, pysmo_rbf_trainer):
        pysmo_rbf_trainer.config.basis_function = None
        pysmo_rbf_trainer.config.regularization = "True"
        pysmo_rbf_trainer.config.solution_method = None

        output_label = "z5"
        data = {"x1": [1, 2, 3, 4], "x2": [5, 6, 7, 8], "z1": [10, 20, 30, 40]}
        data = pd.DataFrame(data)

        model = pysmo_rbf_trainer._create_model(data, output_label)

        assert model.x_data_columns == ["x1", "x2"]
        np.testing.assert_array_equal(model.x_data_unscaled, data.values[:, :-1])
        np.testing.assert_array_equal(model.y_data_unscaled[:, 0], data.values[:, -1])
        assert model.overwrite == True
        assert model.basis_function == "gaussian"
        assert model.regularization == True
        assert model.solution_method == "algebraic"
        # assert model.filename == 'pysmo_Nonerbf_z5.pickle'
        assert list(model.feature_list._data.keys()) == data.columns.tolist()[:-1]

    @pytest.mark.unit
    def test_create_model_cubic(self, pysmo_rbf_trainer):
        pysmo_rbf_trainer.config.basis_function = "cubic"
        pysmo_rbf_trainer.config.regularization = "False"
        pysmo_rbf_trainer.config.solution_method = "pyomo"

        output_label = "z5"
        data = {"x1": [1, 2, 3, 4], "x2": [5, 6, 7, 8], "z1": [10, 20, 30, 40]}
        data = pd.DataFrame(data)

        model = pysmo_rbf_trainer._create_model(data, output_label)

        assert model.x_data_columns == ["x1", "x2"]
        np.testing.assert_array_equal(model.x_data_unscaled, data.values[:, :-1])
        np.testing.assert_array_equal(model.y_data_unscaled[:, 0], data.values[:, -1])
        assert model.overwrite == True
        assert model.basis_function == "cubic"
        assert model.regularization == False
        assert model.solution_method == "pyomo"
        assert model.filename == "solution.pickle"
        assert list(model.feature_list._data.keys()) == data.columns.tolist()[:-1]

    @pytest.mark.unit
    def test_create_model_imq(self, pysmo_rbf_trainer):
        pysmo_rbf_trainer.config.basis_function = "imq"
        pysmo_rbf_trainer.config.regularization = True
        pysmo_rbf_trainer.config.solution_method = "bfgs"

        output_label = "z5"
        data = {"x1": [1, 2, 3, 4], "x2": [5, 6, 7, 8], "z1": [10, 20, 30, 40]}
        data = pd.DataFrame(data)

        model = pysmo_rbf_trainer._create_model(data, output_label)

        assert model.x_data_columns == ["x1", "x2"]
        np.testing.assert_array_equal(model.x_data_unscaled, data.values[:, :-1])
        np.testing.assert_array_equal(model.y_data_unscaled[:, 0], data.values[:, -1])
        assert model.overwrite == True
        assert model.basis_function == "imq"
        assert model.regularization == True
        assert model.solution_method == "bfgs"
        # assert model.filename == 'pysmo_Nonerbf_z5.pickle'
        assert list(model.feature_list._data.keys()) == data.columns.tolist()[:-1]


class TestPysmoKrigingTrainer:
    @pytest.fixture
    def pysmo_krg_trainer(self):

        data = {"x1": [1, 2, 3, 4], "x2": [5, 6, 7, 8], "z1": [10, 20, 30, 40]}
        data = pd.DataFrame(data)
        input_labels = ["x1", "x2"]
        output_labels = ["z1"]
        bnds = {"x1": (0, 5), "x2": (0, 10)}
        krg_trainer = PysmoKrigingTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=bnds,
            training_dataframe=data,
        )
        return krg_trainer

    @pytest.mark.unit
    def test_defaults(self, pysmo_krg_trainer):
        # Check all defaults
        assert pysmo_krg_trainer.model_type == "kriging"
        assert pysmo_krg_trainer.config.numerical_gradients == True
        assert pysmo_krg_trainer.config.regularization == True

    @pytest.mark.unit
    def test_set_regularization_righttype_1(self, pysmo_krg_trainer):
        pysmo_krg_trainer.config.regularization = True
        assert pysmo_krg_trainer.config.regularization == True

    @pytest.mark.unit
    def test_set_regularization_righttype_2(self, pysmo_krg_trainer):
        pysmo_krg_trainer.config.regularization = False
        assert pysmo_krg_trainer.config.regularization == False

    @pytest.mark.unit
    def test_set_regularization_righttype_3(self, pysmo_krg_trainer):
        pysmo_krg_trainer.config.regularization = "False"
        assert pysmo_krg_trainer.config.regularization == False

    @pytest.mark.unit
    def test_set_regularization_righttype_4(self, pysmo_krg_trainer):
        pysmo_krg_trainer.config.regularization = "True"
        assert pysmo_krg_trainer.config.regularization == True

    @pytest.mark.unit
    def test_set_regularization_righttype_5(self, pysmo_krg_trainer):
        pysmo_krg_trainer.config.regularization = 1
        assert pysmo_krg_trainer.config.regularization == True

    @pytest.mark.unit
    def test_set_regularization_righttype_6(self, pysmo_krg_trainer):
        pysmo_krg_trainer.config.regularization = 0
        assert pysmo_krg_trainer.config.regularization == False

    @pytest.mark.unit
    def test_set_regularization_wrongtype(self, pysmo_krg_trainer):
        with pytest.raises(ValueError):
            pysmo_krg_trainer.config.regularization = 2

    @pytest.mark.unit
    def test_set_numerical_gradients_righttype_1(self, pysmo_krg_trainer):
        pysmo_krg_trainer.config.numerical_gradients = True
        assert pysmo_krg_trainer.config.numerical_gradients == True

    @pytest.mark.unit
    def test_set_numerical_gradients_righttype_2(self, pysmo_krg_trainer):
        pysmo_krg_trainer.config.numerical_gradients = False
        assert pysmo_krg_trainer.config.numerical_gradients == False

    @pytest.mark.unit
    def test_set_numerical_gradients_righttype_3(self, pysmo_krg_trainer):
        pysmo_krg_trainer.config.numerical_gradients = "False"
        assert pysmo_krg_trainer.config.numerical_gradients == False

    @pytest.mark.unit
    def test_set_numerical_gradients_righttype_4(self, pysmo_krg_trainer):
        pysmo_krg_trainer.config.numerical_gradients = "True"
        assert pysmo_krg_trainer.config.numerical_gradients == True

    @pytest.mark.unit
    def test_set_numerical_gradients_righttype_5(self, pysmo_krg_trainer):
        pysmo_krg_trainer.config.numerical_gradients = 1
        assert pysmo_krg_trainer.config.numerical_gradients == True

    @pytest.mark.unit
    def test_set_numerical_gradients_righttype_6(self, pysmo_krg_trainer):
        pysmo_krg_trainer.config.numerical_gradients = 0
        assert pysmo_krg_trainer.config.numerical_gradients == False

    @pytest.mark.unit
    def test_set_numerical_gradients_wrongtype(self, pysmo_krg_trainer):
        with pytest.raises(ValueError):
            pysmo_krg_trainer.config.numerical_gradients = 2

    @pytest.mark.unit
    def test_create_model_defaults(self, pysmo_krg_trainer):
        output_label = "z5"
        data = {"x1": [1, 2, 3, 4], "x2": [5, 6, 7, 8], "z1": [10, 20, 30, 40]}
        data = pd.DataFrame(data)

        model = pysmo_krg_trainer._create_model(data, output_label)

        assert model.x_data_columns == ["x1", "x2"]
        np.testing.assert_array_equal(model.x_data, data.values[:, :-1])
        np.testing.assert_array_equal(model.y_data[:, 0], data.values[:, -1])
        assert model.overwrite == True
        assert model.regularization == True
        assert model.num_grads == True
        assert model.filename == "solution.pickle"
        assert model.num_vars == data.shape[1]
        assert list(model.feature_list._data.keys()) == data.columns.tolist()[:-1]

    @pytest.mark.unit
    def test_create_model_no_regularization(self, pysmo_krg_trainer):
        output_label = "z5"
        data = {"x1": [1, 2, 3, 4], "x2": [5, 6, 7, 8], "z1": [10, 20, 30, 40]}
        data = pd.DataFrame(data)

        pysmo_krg_trainer.config.regularization = False
        model = pysmo_krg_trainer._create_model(data, output_label)

        assert model.x_data_columns == ["x1", "x2"]
        np.testing.assert_array_equal(model.x_data, data.values[:, :-1])
        np.testing.assert_array_equal(model.y_data[:, 0], data.values[:, -1])
        assert model.overwrite == True
        assert model.regularization == False
        assert model.num_grads == True
        assert model.filename == "solution.pickle"
        assert model.num_vars == data.shape[1]
        assert list(model.feature_list._data.keys()) == data.columns.tolist()[:-1]

    @pytest.mark.unit
    def test_create_model_no_numerical_grads(self, pysmo_krg_trainer):
        output_label = "z5"
        data = {"x1": [1, 2, 3, 4], "x2": [5, 6, 7, 8], "z1": [10, 20, 30, 40]}
        data = pd.DataFrame(data)

        pysmo_krg_trainer.config.numerical_gradients = "False"
        model = pysmo_krg_trainer._create_model(data, output_label)

        assert model.x_data_columns == ["x1", "x2"]
        np.testing.assert_array_equal(model.x_data, data.values[:, :-1])
        np.testing.assert_array_equal(model.y_data[:, 0], data.values[:, -1])
        assert model.overwrite == True
        assert model.regularization == True
        assert model.num_grads == False
        assert model.filename == "solution.pickle"
        assert model.num_vars == data.shape[1]
        assert list(model.feature_list._data.keys()) == data.columns.tolist()[:-1]


class TestPysmoSurrogate:
    @pytest.fixture
    def pysmo_surr1(self):
        training_data = {
            "x1": [1, 2, 3, 4, 5],
            "x2": [5, 6, 7, 8, 9],
            "z1": [10, 20, 30, 40, 50],
        }  # , 'z2': [6, 8, 10, 12, 14]}
        training_data = pd.DataFrame(training_data)
        validation_data = {
            "x1": [1, 2, 3, 4],
            "x2": [5, 6, 7, 8],
            "z1": [10, 20, 30, 40],
        }  # , 'z2': [6, 8, 10, 12]}#{'x1': [2.5], 'x2': [6.5], 'z1': [25], 'z2': [9]}
        validation_data = pd.DataFrame(validation_data)
        input_labels = ["x1", "x2"]
        output_labels = ["z1"]
        bnds = {"x1": (0, 5), "x2": (0, 10)}

        pysmo_trainer = PysmoPolyTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=bnds,
            training_dataframe=training_data,
            validation_dataframe=validation_data,
            maximum_polynomial_order=1,
            multinomials=True,
            number_of_crossvalidations=3,
        )
        a1 = pysmo_trainer.train_surrogate()
        pysmo_surr1 = PysmoSurrogate(a1, input_labels, output_labels, bnds)

        return pysmo_surr1

    @pytest.fixture
    def pysmo_surr2_poly(self):
        training_data = {
            "x1": [1, 2, 3, 4, 5],
            "x2": [5, 6, 7, 8, 9],
            "z1": [10, 20, 30, 40, 50],
            "z2": [6, 8, 10, 12, 14],
        }
        training_data = pd.DataFrame(training_data)
        validation_data = {
            "x1": [1, 2, 3, 4],
            "x2": [5, 6, 7, 8],
            "z1": [10, 20, 30, 40],
            "z2": [6, 8, 10, 12],
        }  # {'x1': [2.5], 'x2': [6.5], 'z1': [25], 'z2': [9]}
        validation_data = pd.DataFrame(validation_data)
        input_labels = ["x1", "x2"]
        output_labels = ["z1", "z2"]
        bnds = {"x1": (0, 5), "x2": (0, 10)}

        pysmo_trainer = PysmoPolyTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=bnds,
            training_dataframe=training_data,
            validation_dataframe=validation_data,
            maximum_polynomial_order=1,
            multinomials=True,
            number_of_crossvalidations=3,
        )

        a2_poly = pysmo_trainer.train_surrogate()
        pysmo_surr2_poly = PysmoSurrogate(a2_poly, input_labels, output_labels)

        return (a2_poly, pysmo_surr2_poly)

    @pytest.fixture
    def pysmo_surr2_rbf(self):
        training_data = {
            "x1": [1, 2, 3, 4, 5],
            "x2": [5, 6, 7, 8, 9],
            "z1": [10, 20, 30, 40, 50],
            "z2": [6, 8, 10, 12, 14],
        }
        training_data = pd.DataFrame(training_data)
        validation_data = {
            "x1": [1, 2, 3, 4],
            "x2": [5, 6, 7, 8],
            "z1": [10, 20, 30, 40],
            "z2": [6, 8, 10, 12],
        }  # {'x1': [2.5], 'x2': [6.5], 'z1': [25], 'z2': [9]}
        validation_data = pd.DataFrame(validation_data)
        input_labels = ["x1", "x2"]
        output_labels = ["z1", "z2"]
        bnds = {"x1": (0, 5), "x2": (0, 10)}

        pysmo_trainer2 = PysmoRBFTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=bnds,
            training_dataframe=training_data,
            validation_dataframe=validation_data,
            basis_function="gaussian",
            regularization=False,
        )
        a2_rbf = pysmo_trainer2.train_surrogate()
        pysmo_surr2_rbf = PysmoSurrogate(a2_rbf, input_labels, output_labels, bnds)

        return (a2_rbf, pysmo_surr2_rbf)

    @pytest.fixture
    def pysmo_surr2_krg(self):
        training_data = {
            "x1": [1, 2, 3, 4, 5],
            "x2": [7, 6, 9, 5, 8],
            "z1": [10, 20, 30, 40, 50],
            "z2": [6, 8, 10, 12, 14],
        }
        training_data = pd.DataFrame(training_data)
        validation_data = {
            "x1": [1, 2, 3, 4],
            "x2": [5, 6, 7, 8],
            "z1": [10, 20, 30, 40],
            "z2": [6, 8, 10, 12],
        }  # {'x1': [2.5], 'x2': [6.5], 'z1': [25], 'z2': [9]}
        validation_data = pd.DataFrame(validation_data)
        input_labels = ["x1", "x2"]
        output_labels = ["z1", "z2"]
        bnds = {"x1": (0, 5), "x2": (0, 10)}

        np.random.seed(0)
        pysmo_trainer3 = PysmoKrigingTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=bnds,
            training_dataframe=training_data,
            validation_dataframe=validation_data,
            regularization=True,
            numerical_gradients=True,
        )
        a2_krg = pysmo_trainer3.train_surrogate()
        pysmo_surr2_krg = PysmoSurrogate(a2_krg, input_labels, output_labels, bnds)

        return (a2_krg, pysmo_surr2_krg)

    @pytest.fixture
    def pysmo_surr3(self):
        x1 = [1, 2, 3, 4, 5, 6]
        x2 = [5, 6, 7, 8, 9, 10]
        z1 = [
            3.5 * x1[i] + 2.5 * x2[i] - 1.5 * (sin(x1[i]) + cos(x2[i]))
            for i in range(len(x1))
        ]
        z2 = [
            3.5 * x1[i] - 2.5 * x2[i] + 0.5 * (sin(x1[i]) + cos(x2[i]))
            for i in range(len(x1))
        ]
        x = {"x1": x1, "x2": x2, "z1": z1, "z2": z2}
        training_data = pd.DataFrame(x, columns={"x1", "x2", "z1", "z2"})

        # training_data = pd.DataFrame(x, columns={'x1', 'x2', 'z1', 'z2'})
        validation_data = {
            "x1": [1, 2, 3, 4],
            "x2": [5, 6, 7, 8],
            "z1": [10, 20, 30, 40],
            "z2": [6, 8, 10, 12],
        }  # {'x1': [2.5], 'x2': [6.5], 'z1': [25], 'z2': [9]}
        validation_data = pd.DataFrame(validation_data)
        input_labels = ["x1", "x2"]
        output_labels = ["z1", "z2"]
        bnds = {"x1": (0, 10), "x2": (0, 10)}

        pysmo_trainer = PysmoPolyTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=bnds,
            training_dataframe=training_data,
            validation_dataframe=validation_data,
            maximum_polynomial_order=1,
            multinomials=False,
            extra_features=["sin(x1)", "cos(x2)"],
            number_of_crossvalidations=10,
            solution_method="mle",
        )

        a3 = pysmo_trainer.train_surrogate()
        pysmo_surr3 = PysmoSurrogate(a3, input_labels, output_labels)

        return a3, pysmo_surr3

    @pytest.fixture
    def pysmo_surr4(self):
        training_data = {
            "x1": [1, 2, 3, 4, 5],
            "x2": [5, 6, 7, 8, 9],
            "z1": [10, 20, 30, 40, 50],
            "z2": [6, 8, 10, 12, 14],
        }
        training_data = pd.DataFrame(training_data)
        validation_data = {
            "x1": [1, 2, 3, 4],
            "x2": [5, 6, 7, 8],
            "z1": [10, 20, 30, 40],
            "z2": [6, 8, 10, 12],
        }  # {'x1': [2.5], 'x2': [6.5], 'z1': [25], 'z2': [9]}
        validation_data = pd.DataFrame(validation_data)
        input_labels = ["x1", "x2"]
        output_labels = ["z1", "z2"]
        bnds = {"x1": (0, 5), "x2": (0, 10)}

        pysmo_trainer = PysmoPolyTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=bnds,
            training_dataframe=training_data,
            validation_dataframe=validation_data,
            maximum_polynomial_order=1,
            multinomials=False,
            number_of_crossvalidations=3,
            extra_features=["x1 / x2"],
        )

        a2 = pysmo_trainer.train_surrogate()
        pysmo_surr2 = PysmoSurrogate(a2, input_labels, output_labels)

        return pysmo_surr2

    @pytest.fixture
    def pysmo_surr5_rbf(self):
        training_data = {"x1": [1, 2, 3, 4, 5], "z1": [10, 20, 30, 40, 50]}
        training_data = pd.DataFrame(training_data)
        validation_data = {"x1": [1, 2, 3, 4], "z1": [10, 20, 30, 40]}
        validation_data = pd.DataFrame(validation_data)
        input_labels = ["x1"]
        output_labels = ["z1"]
        bnds = {"x1": (0, 5)}

        pysmo_trainer = PysmoRBFTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=bnds,
            training_dataframe=training_data,
            validation_dataframe=validation_data,
            basis_function="cubic",
        )
        a5_rbf = pysmo_trainer.train_surrogate()
        pysmo_surr5_rbf = PysmoSurrogate(a5_rbf, input_labels, output_labels, bnds)

        return a5_rbf, pysmo_surr5_rbf

    @pytest.fixture
    def pysmo_surr5_krg(self):
        training_data = {"x1": [1, 2, 3, 4, 5], "z1": [10, 20, 30, 40, 50]}
        training_data = pd.DataFrame(training_data)
        validation_data = {"x1": [1, 2, 3, 4], "z1": [10, 20, 30, 40]}
        validation_data = pd.DataFrame(validation_data)
        input_labels = ["x1"]
        output_labels = ["z1"]
        bnds = {"x1": (0, 5)}

        pysmo_trainer2 = PysmoKrigingTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=bnds,
            training_dataframe=training_data,
            validation_dataframe=validation_data,
        )
        a5_krg = pysmo_trainer2.train_surrogate()
        pysmo_surr5_krg = PysmoSurrogate(a5_krg, input_labels, output_labels, bnds)

        return a5_krg, pysmo_surr5_krg

    @pytest.fixture
    def pysmo_surr6(self):
        x1 = [1, 2, 3, 4, 5, 6]
        x2 = [5, 6, 7, 8, 9, 10]
        z1 = [
            3.5 * x1[i] + 2.5 * x2[i] - 1.5 * (exp(x1[i] / x2[i]))
            for i in range(len(x1))
        ]
        z2 = [
            3.5 * x1[i] - 2.5 * x2[i] + 0.5 * (exp(x1[i] / x2[i]))
            for i in range(len(x1))
        ]
        x = {"x1": x1, "x2": x2, "z1": z1, "z2": z2}
        training_data = pd.DataFrame(x, columns={"x1", "x2", "z1", "z2"})

        # training_data = pd.DataFrame(x, columns={'x1', 'x2', 'z1', 'z2'})
        validation_data = {
            "x1": [1, 2, 3, 4],
            "x2": [5, 6, 7, 8],
            "z1": [10, 20, 30, 40],
            "z2": [6, 8, 10, 12],
        }  # {'x1': [2.5], 'x2': [6.5], 'z1': [25], 'z2': [9]}
        validation_data = pd.DataFrame(validation_data)
        input_labels = ["x1", "x2"]
        output_labels = ["z1", "z2"]
        bnds = {"x1": (0, 10), "x2": (0, 10)}

        pysmo_trainer = PysmoPolyTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=bnds,
            training_dataframe=training_data,
            validation_dataframe=validation_data,
            maximum_polynomial_order=1,
            multinomials=False,
            extra_features=["exp(x1/x2)"],
            number_of_crossvalidations=10,
            solution_method="mle",
        )

        a6 = pysmo_trainer.train_surrogate()
        pysmo_surr6 = PysmoSurrogate(a6, input_labels, output_labels)

        return a6, pysmo_surr6

    @pytest.mark.unit
    def test_evaluate_unisurrogate_poly(self, pysmo_surr1):
        # Test ``evaluate_surrogate`` for one output with interaction terms
        x = [
            -2,
            -1.8,
            -1.6,
            -1.4,
            -1.2,
            -1.0,
            -0.8,
            -0.6,
            -0.4,
            -0.2,
            0,
            0.2,
            0.4,
            0.6,
            0.8,
            1.0,
            1.2,
            1.4,
            1.6,
            1.8,
            2.0,
        ]

        inputs = np.array([np.tile(x, len(x)), np.repeat(x, len(x))])
        inputs = pd.DataFrame(inputs.transpose(), columns=["x1", "x2"])

        out = pysmo_surr1.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert pytest.approx(out["z1"][i], rel=1e-8) == (
                -75.26111111111476
                - 8.815277777775934 * inputs["x1"][i]
                + 18.81527777777826 * inputs["x2"][i]
                - 2.2556956302821618e-13 * (inputs["x2"][i] * inputs["x1"][i])
            )

    @pytest.mark.unit
    def test_populate_block_unisurrogate_poly(self, pysmo_surr1):
        # Test ``populate_block`` for one output with interaction terms
        blk = SurrogateBlock(concrete=True)

        blk.build_model(pysmo_surr1)
        blk.display()

        assert isinstance(blk.inputs, Var)
        assert blk.inputs["x1"].bounds == (0, 5)
        assert blk.inputs["x2"].bounds == (0, 10)
        assert isinstance(blk.outputs, Var)
        assert blk.outputs["z1"].bounds == (None, None)
        assert isinstance(blk.pysmo_constraint, Constraint)
        assert len(blk.pysmo_constraint) == 1
        assert str(blk.pysmo_constraint["z1"].body) == (
            "outputs[z1] - (-75.26111111111476 - 8.815277777775934*inputs[x1] + 18.81527777777826*inputs[x2] - 2.2556956302821618e-13*(inputs[x2]*inputs[x1]))"
        )

    @pytest.mark.unit
    def test_evaluate_multisurrogate_poly(self, pysmo_surr2_poly):
        # Test ``evaluate_surrogate`` for multiple output polynomials with interaction terms
        x = [
            -2,
            -1.8,
            -1.6,
            -1.4,
            -1.2,
            -1.0,
            -0.8,
            -0.6,
            -0.4,
            -0.2,
            0,
            0.2,
            0.4,
            0.6,
            0.8,
            1.0,
            1.2,
            1.4,
            1.6,
            1.8,
            2.0,
        ]

        inputs = np.array([np.tile(x, len(x)), np.repeat(x, len(x))])
        inputs = pd.DataFrame(inputs.transpose(), columns=["x1", "x2"])

        _, poly_trained = pysmo_surr2_poly
        out = poly_trained.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert pytest.approx(out["z1"][i], rel=1e-8) == (
                -75.26111111111476
                - 8.815277777775934 * inputs["x1"][i]
                + 18.81527777777826 * inputs["x2"][i]
                - 2.2556956302821618e-13 * (inputs["x2"][i] * inputs["x1"][i])
            )
            assert pytest.approx(out["z2"][i], rel=1e-8) == (
                -3.0033074724377813
                + 0.2491731318906352 * inputs["x1"][i]
                + 1.7508268681094337 * inputs["x2"][i]
                - 6.786238238021269e-15 * (inputs["x2"][i] * inputs["x1"][i])
            )

    @pytest.mark.unit
    def test_populate_block_multisurrogate_poly(self, pysmo_surr2_poly):
        # Test ``populate_block`` for multiple output polynomials with interaction terms
        blk = SurrogateBlock(concrete=True)

        (
            _,
            poly_trained,
        ) = pysmo_surr2_poly
        blk.build_model(poly_trained)

        assert isinstance(blk.inputs, Var)
        assert blk.inputs["x1"].bounds == (None, None)
        assert blk.inputs["x2"].bounds == (None, None)
        assert isinstance(blk.outputs, Var)
        assert blk.outputs["z1"].bounds == (None, None)
        assert blk.outputs["z2"].bounds == (None, None)
        assert isinstance(blk.pysmo_constraint, Constraint)
        assert len(blk.pysmo_constraint) == 2
        assert str(blk.pysmo_constraint["z1"].body) == (
            "outputs[z1] - (-75.26111111111476 - 8.815277777775934*inputs[x1] + 18.81527777777826*inputs[x2] - 2.2556956302821618e-13*(inputs[x2]*inputs[x1]))"
        )
        assert str(blk.pysmo_constraint["z2"].body) == (
            "outputs[z2] - (-3.0033074724377813 + 0.2491731318906352*inputs[x1] + 1.7508268681094337*inputs[x2] - 6.786238238021269e-15*(inputs[x2]*inputs[x1]))"
        )

    @pytest.mark.unit
    def test_evaluate_multisurrogate_poly_trigfuncs1(self, pysmo_surr3):
        # Test ``evaluate_surrogate`` for multiple output polynomials with trig terms
        x = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]

        inputs = np.array([np.tile(x, len(x)), np.repeat(x, len(x))])
        inputs = pd.DataFrame(inputs.transpose(), columns=["x1", "x2"])

        sol, poly_trained = pysmo_surr3
        out = poly_trained.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert pytest.approx(out["z1"][i], rel=1e-6) == (
                sol._data["z1"]._model.optimal_weights_array[0, 0]
                + sol._data["z1"]._model.optimal_weights_array[1, 0] * inputs["x1"][i]
                + sol._data["z1"]._model.optimal_weights_array[2, 0] * inputs["x2"][i]
                + sol._data["z1"]._model.optimal_weights_array[3, 0]
                * sin(inputs["x1"][i])
                + sol._data["z1"]._model.optimal_weights_array[4, 0]
                * cos(inputs["x2"][i])
            )
            assert pytest.approx(out["z2"][i], rel=1e-6) == (
                sol._data["z2"]._model.optimal_weights_array[0, 0]
                + sol._data["z2"]._model.optimal_weights_array[1, 0] * inputs["x1"][i]
                + sol._data["z2"]._model.optimal_weights_array[2, 0] * inputs["x2"][i]
                + sol._data["z2"]._model.optimal_weights_array[3, 0]
                * sin(inputs["x1"][i])
                + sol._data["z2"]._model.optimal_weights_array[4, 0]
                * cos(inputs["x2"][i])
            )

    @pytest.mark.unit
    def test_populate_block_trigfuncs1(self, pysmo_surr3):
        blk = SurrogateBlock(concrete=True)

        sol, poly_trained = pysmo_surr3
        blk.build_model(poly_trained)

        assert isinstance(blk.inputs, Var)
        assert isinstance(blk.outputs, Var)
        assert isinstance(blk.pysmo_constraint, Constraint)
        assert len(blk.pysmo_constraint) == 2
        assert str(blk.pysmo_constraint["z1"].body) == (
            "outputs[z1] - (-{} + {}*inputs[x1] + {}*inputs[x2] - {}*sin(inputs[x1]) - {}*cos(inputs[x2]))".format(
                abs(sol._data["z1"]._model.optimal_weights_array[0, 0]),
                abs(sol._data["z1"]._model.optimal_weights_array[1, 0]),
                abs(sol._data["z1"]._model.optimal_weights_array[2, 0]),
                abs(sol._data["z1"]._model.optimal_weights_array[3, 0]),
                abs(sol._data["z1"]._model.optimal_weights_array[4, 0]),
            )
        )
        assert str(blk.pysmo_constraint["z2"].body) == (
            "outputs[z2] - (-{} + {}*inputs[x1] - {}*inputs[x2] + {}*sin(inputs[x1]) + {}*cos(inputs[x2]))"
        ).format(
            abs(sol._data["z2"]._model.optimal_weights_array[0, 0]),
            abs(sol._data["z2"]._model.optimal_weights_array[1, 0]),
            abs(sol._data["z2"]._model.optimal_weights_array[2, 0]),
            abs(sol._data["z2"]._model.optimal_weights_array[3, 0]),
            abs(sol._data["z2"]._model.optimal_weights_array[4, 0]),
        )

    @pytest.mark.unit
    def test_evaluate_multisurrogate_poly_trigfuncs2(self, pysmo_surr6):
        # Test ``evaluate_surrogate`` for multiple output polynomials with log terms
        x = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]

        inputs = np.array([np.tile(x, len(x)), np.repeat(x, len(x))])
        inputs = pd.DataFrame(inputs.transpose(), columns=["x1", "x2"])

        sol, poly_trained = pysmo_surr6
        out = poly_trained.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert pytest.approx(out["z1"][i], rel=1e-6) == (
                sol._data["z1"]._model.optimal_weights_array[0, 0]
                + sol._data["z1"]._model.optimal_weights_array[1, 0] * inputs["x1"][i]
                + sol._data["z1"]._model.optimal_weights_array[2, 0] * inputs["x2"][i]
                + sol._data["z1"]._model.optimal_weights_array[3, 0]
                * exp(inputs["x1"][i] / inputs["x2"][i])
            )
            assert pytest.approx(out["z2"][i], rel=1e-6) == (
                sol._data["z2"]._model.optimal_weights_array[0, 0]
                + sol._data["z2"]._model.optimal_weights_array[1, 0] * inputs["x1"][i]
                + sol._data["z2"]._model.optimal_weights_array[2, 0] * inputs["x2"][i]
                + sol._data["z2"]._model.optimal_weights_array[3, 0]
                * exp(inputs["x1"][i] / inputs["x2"][i])
            )

    @pytest.mark.unit
    def test_populate_block_trigfuncs2(self, pysmo_surr6):
        blk = SurrogateBlock(concrete=True)

        sol, poly_trained = pysmo_surr6
        blk.build_model(poly_trained)

        assert isinstance(blk.inputs, Var)
        assert isinstance(blk.outputs, Var)
        assert isinstance(blk.pysmo_constraint, Constraint)
        assert len(blk.pysmo_constraint) == 2
        assert str(blk.pysmo_constraint["z1"].body) == (
            "outputs[z1] - (-{} + {}*inputs[x1] + {}*inputs[x2] - {}*exp(inputs[x1]/inputs[x2]))".format(
                abs(sol._data["z1"]._model.optimal_weights_array[0, 0]),
                abs(sol._data["z1"]._model.optimal_weights_array[1, 0]),
                abs(sol._data["z1"]._model.optimal_weights_array[2, 0]),
                abs(sol._data["z1"]._model.optimal_weights_array[3, 0]),
            )
        )
        assert str(blk.pysmo_constraint["z2"].body) == (
            "outputs[z2] - (-{} + {}*inputs[x1] - {}*inputs[x2] + {}*exp(inputs[x1]/inputs[x2]))".format(
                abs(sol._data["z2"]._model.optimal_weights_array[0, 0]),
                abs(sol._data["z2"]._model.optimal_weights_array[1, 0]),
                abs(sol._data["z2"]._model.optimal_weights_array[2, 0]),
                abs(sol._data["z2"]._model.optimal_weights_array[3, 0]),
            )
        )

    @pytest.mark.unit
    def test_evaluate_multisurrogate_poly_userdef(self, pysmo_surr4):
        # Test ``evaluate_surrogate`` for multiple output polynomials with ratio-type user-defined features
        x = [
            -2,
            -1.8,
            -1.6,
            -1.4,
            -1.2,
            -1.0,
            -0.8,
            -0.6,
            -0.4,
            -0.2,
            0.2,
            0.4,
            0.6,
            0.8,
            1.0,
            1.2,
            1.4,
            1.6,
            1.8,
            2.0,
            2.2,
        ]

        inputs = np.array([np.tile(x, len(x)), np.repeat(x, len(x))])
        inputs = pd.DataFrame(inputs.transpose(), columns=["x1", "x2"])

        out = pysmo_surr4.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert pytest.approx(out["z1"][i], rel=1e-8) == (
                -110.15000000001504
                - 17.53750000000189 * inputs["x1"][i]
                + 27.537500000006148 * inputs["x2"][i]
                - 5.3967136315336006e-11 * (inputs["x1"][i] / inputs["x2"][i])
            )
            assert pytest.approx(out["z2"][i], rel=1e-8) == (
                -12.523574144487087
                - 2.1308935361219556 * inputs["x1"][i]
                + 4.1308935361216435 * inputs["x2"][i]
                + 3.6347869158959156e-12 * (inputs["x1"][i] / inputs["x2"][i])
            )

    @pytest.mark.unit
    def test_populate_block_multisurrogate_poly_userdef(self, pysmo_surr4):
        # Test ``populate_block`` for multiple output polynomials with interaction terms
        blk = SurrogateBlock(concrete=True)

        blk.build_model(pysmo_surr4)

        assert isinstance(blk.inputs, Var)
        assert blk.inputs["x1"].bounds == (None, None)
        assert blk.inputs["x2"].bounds == (None, None)
        assert isinstance(blk.outputs, Var)
        assert blk.outputs["z1"].bounds == (None, None)
        assert blk.outputs["z2"].bounds == (None, None)
        assert isinstance(blk.pysmo_constraint, Constraint)
        assert len(blk.pysmo_constraint) == 2
        assert str(blk.pysmo_constraint["z1"].body) == (
            "outputs[z1] - (-110.15000000001504 - 17.53750000000189*inputs[x1] + 27.537500000006148*inputs[x2] - 5.3967136315336006e-11*(inputs[x1]/inputs[x2]))"
        )
        assert str(blk.pysmo_constraint["z2"].body) == (
            "outputs[z2] - (-12.523574144487087 - 2.1308935361219556*inputs[x1] + 4.1308935361216435*inputs[x2] + 3.6347869158959156e-12*(inputs[x1]/inputs[x2]))"
        )

    @pytest.mark.parametrize(
        "confidence_dict", [{0.99: 3.2498355440153697}, {0.90: 1.8331129326536335}]
    )
    @pytest.mark.unit
    def test_confint_default(self, confidence_dict):
        # Test that the ``get_confidence_intervals`` function returns the correct upper and lower confidence interval bounds.
        training_data = {
            "x1": [1, 2, 3, 4, 5],
            "x2": [5, 6, 7, 8, 9],
            "z1": [10, 20, 30, 40, 50],
            "z2": [6, 8, 10, 12, 14],
        }
        training_data = pd.DataFrame(training_data)
        validation_data = {
            "x1": [1, 2, 3, 4],
            "x2": [5, 6, 7, 8],
            "z1": [10, 20, 30, 40],
            "z2": [6, 8, 10, 12],
        }
        validation_data = pd.DataFrame(validation_data)
        input_labels = ["x1", "x2"]
        output_labels = ["z1", "z2"]
        bnds = {"x1": (0, 5), "x2": (0, 10)}

        pysmo_trainer = PysmoPolyTrainer(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=bnds,
            training_dataframe=training_data,
            validation_dataframe=validation_data,
            maximum_polynomial_order=1,
            multinomials=True,
            number_of_crossvalidations=3,
        )

        a2_poly = pysmo_trainer.train_surrogate()
        for k in confidence_dict.keys():
            confidence = k
            tval = confidence_dict[k]

        output = pysmo_trainer.get_confidence_intervals(a2_poly, confidence)

        reg_coeffs = {
            "z1": np.array(
                [
                    -75.26111111111476,
                    -8.815277777775934,
                    18.81527777777826,
                    -2.2556956302821618e-13,
                ]
            ),
            "z2": np.array(
                [
                    -3.0033074724377813,
                    0.2491731318906352,
                    1.7508268681094337,
                    6.786238238021269e-15,
                ]
            ),
        }

        # Test that output has the right number of dictionary entries
        assert len(output) == len(output_labels)
        for i in output_labels:
            # Test that the lower confidence bounds are correctly calculated.
            assert pytest.approx(output[i]["Conf. int. lower"].values, abs=1e-9) == (
                reg_coeffs[i] - tval * output[i]["Std. error"].values
            )
            # Test that the upper confidence bounds are correctly calculated.
            assert pytest.approx(output[i]["Conf. int. upper"].values, abs=1e-9) == (
                reg_coeffs[i] + tval * output[i]["Std. error"].values
            )

    @pytest.mark.unit
    def test_evaluate_unisurrogate_rbf(self, pysmo_surr5_rbf):
        # Test ``evaluate_surrogate`` for RBF with one input/output
        x = [1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
        inputs = pd.DataFrame(x, columns=["x1"])

        sol, rbf_trained = pysmo_surr5_rbf
        out = rbf_trained.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert pytest.approx(out["z1"][i], rel=1e-8) == (
                10
                + 40
                * (
                    sol._data["z1"]._model.weights[0, 0]
                    * ((((inputs["x1"][i] - 1) / 4) ** 2) ** 0.5) ** 3
                    + sol._data["z1"]._model.weights[1, 0]
                    * ((((inputs["x1"][i] - 1) / 4 - 0.25) ** 2) ** 0.5) ** 3
                    + sol._data["z1"]._model.weights[2, 0]
                    * ((((inputs["x1"][i] - 1) / 4 - 0.5) ** 2) ** 0.5) ** 3
                    + sol._data["z1"]._model.weights[3, 0]
                    * ((((inputs["x1"][i] - 1) / 4 - 0.75) ** 2) ** 0.5) ** 3
                    + sol._data["z1"]._model.weights[4, 0]
                    * ((((inputs["x1"][i] - 1) / 4 - 1.0) ** 2) ** 0.5) ** 3
                )
            )

    @pytest.mark.unit
    def test_populate_block_unisurrogate_rbf(self, pysmo_surr5_rbf):
        # Test ``populate_block`` for RBF with one input/output
        blk = SurrogateBlock(concrete=True)

        sol, rbf_trained = pysmo_surr5_rbf
        blk.build_model(rbf_trained)
        blk.display()

        assert isinstance(blk.inputs, Var)
        assert blk.inputs["x1"].bounds == (0, 5)
        assert isinstance(blk.outputs, Var)
        assert blk.outputs["z1"].bounds == (None, None)
        assert isinstance(blk.pysmo_constraint, Constraint)
        assert len(blk.pysmo_constraint) == 1
        assert str(blk.pysmo_constraint["z1"].body) == (
            "outputs[z1] - (10 + 40*(-{}*((((inputs[x1] - 1)/4)**2)**0.5)**3 + {}*((((inputs[x1] - 1)/4 - 0.25)**2)**0.5)**3 + {}*((((inputs[x1] - 1)/4 - 0.5)**2)**0.5)**3 + {}*((((inputs[x1] - 1)/4 - 0.75)**2)**0.5)**3 - {}*((((inputs[x1] - 1)/4 - 1.0)**2)**0.5)**3))".format(
                abs(sol._data["z1"]._model.weights[0, 0]),
                abs(sol._data["z1"]._model.weights[1, 0]),
                abs(sol._data["z1"]._model.weights[2, 0]),
                abs(sol._data["z1"]._model.weights[3, 0]),
                abs(sol._data["z1"]._model.weights[4, 0]),
            )
        )

    @pytest.mark.unit
    def test_evaluate_multisurrogate_rbf(self, pysmo_surr2_rbf):
        # Test ``evaluate_surrogate`` for RBF with two inputs/outputs
        x = [
            -2,
            -1.8,
            -1.6,
            -1.4,
            -1.2,
            -1.0,
            -0.8,
            -0.6,
            -0.4,
            -0.2,
            0,
            0.2,
            0.4,
            0.6,
            0.8,
            1.0,
            1.2,
            1.4,
            1.6,
            1.8,
            2.0,
        ]

        inputs = np.array([np.tile(x, len(x)), np.repeat(x, len(x))])
        inputs = pd.DataFrame(inputs.transpose(), columns=["x1", "x2"])

        sol, rbf_trained = pysmo_surr2_rbf
        out = rbf_trained.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert pytest.approx(out["z1"][i], rel=1e-6, abs=1e-8) == (
                (
                    10
                    + 40
                    * (
                        sol._data["z1"]._model.weights[0, 0]
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4) ** 2
                                        + ((inputs["x2"][i] - 5) / 4) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                        + sol._data["z1"]._model.weights[1, 0]
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4 - 0.25) ** 2
                                        + ((inputs["x2"][i] - 5) / 4 - 0.25) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                        + sol._data["z1"]._model.weights[2, 0]
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4 - 0.5) ** 2
                                        + ((inputs["x2"][i] - 5) / 4 - 0.5) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                        + sol._data["z1"]._model.weights[3, 0]
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4 - 0.75) ** 2
                                        + ((inputs["x2"][i] - 5) / 4 - 0.75) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                        + sol._data["z1"]._model.weights[4, 0]
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4 - 1.0) ** 2
                                        + ((inputs["x2"][i] - 5) / 4 - 1.0) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                    )
                )
            )
            assert pytest.approx(out["z2"][i], rel=1e-6, abs=1e-8) == (
                (
                    6
                    + 8
                    * (
                        sol._data["z2"]._model.weights[0, 0]
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4) ** 2
                                        + ((inputs["x2"][i] - 5) / 4) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                        + sol._data["z2"]._model.weights[1, 0]
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4 - 0.25) ** 2
                                        + ((inputs["x2"][i] - 5) / 4 - 0.25) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                        + sol._data["z2"]._model.weights[2, 0]
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4 - 0.5) ** 2
                                        + ((inputs["x2"][i] - 5) / 4 - 0.5) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                        + sol._data["z2"]._model.weights[3, 0]
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4 - 0.75) ** 2
                                        + ((inputs["x2"][i] - 5) / 4 - 0.75) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                        + sol._data["z2"]._model.weights[4, 0]
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4 - 1.0) ** 2
                                        + ((inputs["x2"][i] - 5) / 4 - 1.0) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                    )
                )
            )

    @pytest.mark.unit
    def test_populate_block_multisurrogate_rbf(self, pysmo_surr2_rbf):
        # Test ``populate_block`` for RBF with one input/output
        blk = SurrogateBlock(concrete=True)

        sol, rbf_trained = pysmo_surr2 = pysmo_surr2_rbf
        blk.build_model(rbf_trained)
        blk.display()

        assert isinstance(blk.inputs, Var)
        assert blk.inputs["x1"].bounds == (0, 5)
        assert isinstance(blk.outputs, Var)
        assert blk.outputs["z1"].bounds == (None, None)
        assert isinstance(blk.pysmo_constraint, Constraint)
        assert len(blk.pysmo_constraint) == 2

    @pytest.mark.unit
    def test_evaluate_multisurrogate_kriging(self, pysmo_surr2_krg):
        # Test ``evaluate_surrogate`` for kriging with two inputs/outputs
        x = [1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0]

        inputs = np.array([np.tile(x, len(x)), np.repeat([i + 4 for i in x], len(x))])
        inputs = pd.DataFrame(inputs.transpose(), columns=["x1", "x2"])

        _, krg_trained = pysmo_surr2_krg
        out = krg_trained.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert pytest.approx(out["z1"][i], rel=1e-2) == (
                -7507.921579707475
                * exp(
                    -(
                        0.044124735560275304 * ((inputs["x1"][i] - 1) / 4) ** 2
                        + 0.001 * ((inputs["x2"][i] - 5) / 4 - 0.5) ** 2
                    )
                )
                + 15454.174740136026
                * exp(
                    -(
                        0.044124735560275304 * ((inputs["x1"][i] - 1) / 4 - 0.25) ** 2
                        + 0.001 * ((inputs["x2"][i] - 5) / 4 - 0.25) ** 2
                    )
                )
                - 4531.029221784704
                * exp(
                    -(
                        0.044124735560275304 * ((inputs["x1"][i] - 1) / 4 - 0.5) ** 2
                        + 0.001 * ((inputs["x2"][i] - 5) / 4 - 1.0) ** 2
                    )
                )
                - 9320.651261114403
                * exp(
                    -(
                        0.044124735560275304 * ((inputs["x1"][i] - 1) / 4 - 0.75) ** 2
                        + 0.001 * ((inputs["x2"][i] - 5) / 4) ** 2
                    )
                )
                + 5905.427322470428
                * exp(
                    -(
                        0.044124735560275304 * ((inputs["x1"][i] - 1) / 4 - 1.0) ** 2
                        + 0.001 * ((inputs["x2"][i] - 5) / 4 - 0.75) ** 2
                    )
                )
                + 27.201100480392512
            )
            assert pytest.approx(out["z2"][i], rel=1e-2) == (
                -1501.5841957038108
                * exp(
                    -(
                        0.04412462103209179 * ((inputs["x1"][i] - 1) / 4) ** 2
                        + 0.001 * ((inputs["x2"][i] - 5) / 4 - 0.5) ** 2
                    )
                )
                + 3090.83421508662
                * exp(
                    -(
                        0.04412462103209179 * ((inputs["x1"][i] - 1) / 4 - 0.25) ** 2
                        + 0.001 * ((inputs["x2"][i] - 5) / 4 - 0.25) ** 2
                    )
                )
                - 906.2056548932875
                * exp(
                    -(
                        0.04412462103209179 * ((inputs["x1"][i] - 1) / 4 - 0.5) ** 2
                        + 0.001 * ((inputs["x2"][i] - 5) / 4 - 1.0) ** 2
                    )
                )
                - 1864.1297383033075
                * exp(
                    -(
                        0.04412462103209179 * ((inputs["x1"][i] - 1) / 4 - 0.75) ** 2
                        + 0.001 * ((inputs["x2"][i] - 5) / 4) ** 2
                    )
                )
                + 1181.085373813784
                * exp(
                    -(
                        0.04412462103209179 * ((inputs["x1"][i] - 1) / 4 - 1.0) ** 2
                        + 0.001 * ((inputs["x2"][i] - 5) / 4 - 0.75) ** 2
                    )
                )
                + 9.440220252159808
            )

    @pytest.mark.unit
    def test_populate_block_multisurrogate_kriging(self, pysmo_surr2_krg):
        # Test ``populate_block`` for kriging with two inputs/outputs
        blk = SurrogateBlock(concrete=True)

        _, krg_trained = pysmo_surr2_krg
        blk.build_model(krg_trained)
        blk.display()

        assert isinstance(blk.inputs, Var)
        assert blk.inputs["x1"].bounds == (0, 5)
        assert isinstance(blk.outputs, Var)
        assert blk.outputs["z1"].bounds == (None, None)
        assert isinstance(blk.pysmo_constraint, Constraint)
        assert len(blk.pysmo_constraint) == 2

    @pytest.mark.unit
    def test_save_poly1(self, pysmo_surr1):
        # Test save for polynomial regression with single output with bounds supplied
        stream1 = StringIO()
        pysmo_surr1.save(stream1)
        assert re.sub("errors.*?}", "", jstring_poly_1) == re.sub(
            "errors.*?}", "", stream1.getvalue()
        )

    @pytest.mark.unit
    def test_save_poly2(self, pysmo_surr2_poly):
        # Test save for polynomial regression with multiple outputs - most complicated to save
        _, poly_trained = pysmo_surr2_poly
        stream2a = StringIO()
        poly_trained.save(stream2a)
        assert re.sub("errors.*?}", "", jstring_poly_2) == re.sub(
            "errors.*?}", "", stream2a.getvalue()
        )

    @pytest.mark.unit
    def test_save_poly3(self, pysmo_surr4):
        # Test save for polynomial regression with multiple outputs and custom functions - most complicated to save
        stream4 = StringIO()
        pysmo_surr4.save(stream4)
        assert re.sub("errors.*?}", "", jstring_poly_4) == re.sub(
            "errors.*?}", "", stream4.getvalue()
        )

    @pytest.mark.unit
    def test_load_poly1(self):
        # Try for polynomial with single output and bounds
        m = ConcreteModel()
        m.inputs = Var(["x1", "x2"])

        # Validation data
        x = [
            -2,
            -1.8,
            -1.6,
            -1.4,
            -1.2,
            -1.0,
            -0.8,
            -0.6,
            -0.4,
            -0.2,
            0,
            0.2,
            0.4,
            0.6,
            0.8,
            1.0,
            1.2,
            1.4,
            1.6,
            1.8,
            2.0,
        ]
        inputs = np.array([np.tile(x, len(x)), np.repeat(x, len(x))])
        inputs = pd.DataFrame(inputs.transpose(), columns=["x1", "x2"])

        pysmo_surr_poly = PysmoSurrogate.load(StringIO(jstring_poly_1))
        # Assert labels, outputs and bounds
        assert pysmo_surr_poly._input_labels == ["x1", "x2"]
        assert pysmo_surr_poly._output_labels == ["z1"]
        assert pysmo_surr_poly._input_bounds == {"x1": (0, 5), "x2": (0, 10)}
        assert pysmo_surr_poly._trained.model_type == "poly"
        assert pysmo_surr_poly._trained.output_labels == ["z1"]
        assert len(pysmo_surr_poly._trained._data) == 1
        assert list(pysmo_surr_poly._trained._data) == ["z1"]
        # Assert that correcte xpression string was returned
        assert pysmo_surr_poly._trained._data["z1"].expression_str == (
            "-75.26111111111476 -8.815277777775934*IndexedParam[x1] + 18.81527777777826*IndexedParam[x2]"
            " -2.2556956302821618e-13*(IndexedParam[x2]*IndexedParam[x1])"
        )
        # Assert that correct model is returned with generate_expression()
        assert str(
            pysmo_surr_poly._trained._data["z1"]._model.generate_expression(
                [m.inputs["x1"], m.inputs["x2"]]
            )
        ) == (
            "-75.26111111111476 - 8.815277777775934*inputs[x1] + 18.81527777777826*inputs[x2]"
            " - 2.2556956302821618e-13*(inputs[x2]*inputs[x1])"
        )
        # Assert that returned results will evaluate and return correct results
        out = pysmo_surr_poly.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert pytest.approx(out["z1"][i], rel=1e-8) == (
                -75.26111111111476
                - 8.815277777775934 * inputs["x1"][i]
                + 18.81527777777826 * inputs["x2"][i]
                - 2.2556956302821618e-13 * (inputs["x2"][i] * inputs["x1"][i])
            )

    @pytest.mark.unit
    def test_load_poly2(self):
        # Try for polynomial with multiple outputs and trig functions
        m = ConcreteModel()
        m.inputs = Var(["x1", "x2"])

        # Validation data
        x = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]

        inputs = np.array([np.tile(x, len(x)), np.repeat(x, len(x))])
        inputs = pd.DataFrame(inputs.transpose(), columns=["x1", "x2"])

        pysmo_surr_poly = PysmoSurrogate.load(StringIO(jstring_poly_3))
        # Assert labels, outputs and bounds
        assert pysmo_surr_poly._input_labels == ["x1", "x2"]
        assert pysmo_surr_poly._output_labels == ["z1", "z2"]
        assert pysmo_surr_poly._input_bounds == None
        assert pysmo_surr_poly._trained.model_type == "poly"
        assert pysmo_surr_poly._trained.output_labels == ["z1", "z2"]
        assert len(pysmo_surr_poly._trained._data) == 2
        assert list(pysmo_surr_poly._trained._data) == ["z1", "z2"]
        # Assert that correct expression strings were returned for both surrogates
        assert pysmo_surr_poly._trained._data["z1"].expression_str == (
            "-14.290243902439855 + 6.4274390243899795*IndexedParam[x1] + 3.572560975609962*IndexedParam[x2] "
            "+ 1.9753643165643098e-13*log(IndexedParam[x1]) -4.4048098502003086e-14*sin(IndexedParam[x2])"
        )
        assert pysmo_surr_poly._trained._data["z2"].expression_str == (
            "5.704971042443143 + 2.4262427606248815*IndexedParam[x1] -0.42624276060821653*IndexedParam[x2] "
            "-5.968545102597034e-11*log(IndexedParam[x1]) + 6.481176706429892e-12*sin(IndexedParam[x2])"
        )
        # Assert that correct model is returned with generate_expression()
        assert str(
            pysmo_surr_poly._trained._data["z1"]._model.generate_expression(
                [m.inputs["x1"], m.inputs["x2"]]
            )
        ) == (
            "-14.290243902439855 + 6.4274390243899795*inputs[x1] + 3.572560975609962*inputs[x2] "
            "+ 1.9753643165643098e-13*log(inputs[x1]) - 4.4048098502003086e-14*sin(inputs[x2])"
        )
        assert str(
            pysmo_surr_poly._trained._data["z2"]._model.generate_expression(
                [m.inputs["x1"], m.inputs["x2"]]
            )
        ) == (
            "5.704971042443143 + 2.4262427606248815*inputs[x1] - 0.42624276060821653*inputs[x2] "
            "- 5.968545102597034e-11*log(inputs[x1]) + 6.481176706429892e-12*sin(inputs[x2])"
        )
        # Assert that returned results will evaluate and return correct results
        out = pysmo_surr_poly.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert pytest.approx(out["z1"][i], rel=1e-8) == (
                -14.290243902439855
                + 6.4274390243899795 * inputs["x1"][i]
                + 3.572560975609962 * inputs["x2"][i]
                + 1.9753643165643098e-13 * log(inputs["x1"][i])
                - 4.4048098502003086e-14 * sin(inputs["x2"][i])
            )
            assert pytest.approx(out["z2"][i], rel=1e-8) == (
                5.704971042443143
                + 2.4262427606248815 * inputs["x1"][i]
                - 0.42624276060821653 * inputs["x2"][i]
                - 5.968545102597034e-11 * log(inputs["x1"][i])
                + 6.481176706429892e-12 * sin(inputs["x2"][i])
            )

    @pytest.mark.unit
    def test_load_poly3(self):
        # Try for polynomial with multiple outputs and other user-defined functions
        m = ConcreteModel()
        m.inputs = Var(["x1", "x2"])

        # Validation data
        x = [
            -2,
            -1.8,
            -1.6,
            -1.4,
            -1.2,
            -1.0,
            -0.8,
            -0.6,
            -0.4,
            -0.2,
            0.2,
            0.4,
            0.6,
            0.8,
            1.0,
            1.2,
            1.4,
            1.6,
            1.8,
            2.0,
            2.2,
        ]

        inputs = np.array([np.tile(x, len(x)), np.repeat(x, len(x))])
        inputs = pd.DataFrame(inputs.transpose(), columns=["x1", "x2"])

        pysmo_surr_poly = PysmoSurrogate.load(StringIO(jstring_poly_4))
        # Assert labels, outputs and bounds
        assert pysmo_surr_poly._input_labels == ["x1", "x2"]
        assert pysmo_surr_poly._output_labels == ["z1", "z2"]
        assert pysmo_surr_poly._input_bounds == None
        assert pysmo_surr_poly._trained.model_type == "poly"
        assert pysmo_surr_poly._trained.output_labels == ["z1", "z2"]
        assert len(pysmo_surr_poly._trained._data) == 2
        assert list(pysmo_surr_poly._trained._data) == ["z1", "z2"]
        # Assert that correct expression strings were returned for both surrogates
        assert pysmo_surr_poly._trained._data["z1"].expression_str == (
            "-110.15000000001504 -17.53750000000189*IndexedParam[x1] + 27.537500000006148*IndexedParam[x2] "
            "-5.3967136315336006e-11*(IndexedParam[x1]/IndexedParam[x2])"
        )
        assert pysmo_surr_poly._trained._data["z2"].expression_str == (
            "-12.523574144487087 -2.1308935361219556*IndexedParam[x1] + 4.1308935361216435*IndexedParam[x2]"
            " + 3.6347869158959156e-12*(IndexedParam[x1]/IndexedParam[x2])"
        )
        # Assert that correct model is returned with generate_expression()
        assert str(
            pysmo_surr_poly._trained._data["z1"]._model.generate_expression(
                [m.inputs["x1"], m.inputs["x2"]]
            )
        ) == (
            "-110.15000000001504 - 17.53750000000189*inputs[x1] + 27.537500000006148*inputs[x2] "
            "- 5.3967136315336006e-11*(inputs[x1]/inputs[x2])"
        )
        assert str(
            pysmo_surr_poly._trained._data["z2"]._model.generate_expression(
                [m.inputs["x1"], m.inputs["x2"]]
            )
        ) == (
            "-12.523574144487087 - 2.1308935361219556*inputs[x1] + 4.1308935361216435*inputs[x2] "
            "+ 3.6347869158959156e-12*(inputs[x1]/inputs[x2])"
        )
        # Assert that returned results will evaluate and return correct results
        out = pysmo_surr_poly.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert pytest.approx(out["z1"][i], rel=1e-6) == (
                -110.15000000001504
                - 17.53750000000189 * inputs["x1"][i]
                + 27.537500000006148 * inputs["x2"][i]
                - 5.3967136315336006e-11 * (inputs["x1"][i] / inputs["x2"][i])
            )
            assert pytest.approx(out["z2"][i], rel=1e-6) == (
                -12.523574144487087
                - 2.1308935361219556 * inputs["x1"][i]
                + 4.1308935361216435 * inputs["x2"][i]
                + 3.6347869158959156e-12 * (inputs["x1"][i] / inputs["x2"][i])
            )

    @pytest.mark.unit
    def test_load_rbf(self):
        # Try for polynomial with multiple outputs and other user-defined functions
        m = ConcreteModel()
        m.inputs = Var(["x1", "x2"])

        # Validation data
        x = [
            -2,
            -1.8,
            -1.6,
            -1.4,
            -1.2,
            -1.0,
            -0.8,
            -0.6,
            -0.4,
            -0.2,
            0.2,
            0.4,
            0.6,
            0.8,
            1.0,
            1.2,
            1.4,
            1.6,
            1.8,
            2.0,
            2.2,
        ]

        inputs = np.array([np.tile(x, len(x)), np.repeat(x, len(x))])
        inputs = pd.DataFrame(inputs.transpose(), columns=["x1", "x2"])

        pysmo_surr_rbf = PysmoSurrogate.load(StringIO(jstring_rbf))
        # Assert labels, outputs and bounds
        assert pysmo_surr_rbf._input_labels == ["x1", "x2"]
        assert pysmo_surr_rbf._output_labels == ["z1", "z2"]
        assert pysmo_surr_rbf._input_bounds == {"x1": (0, 5), "x2": (0, 10)}
        assert pysmo_surr_rbf._trained.model_type == "rbf"
        assert pysmo_surr_rbf._trained.output_labels == ["z1", "z2"]
        assert len(pysmo_surr_rbf._trained._data) == 2
        assert list(pysmo_surr_rbf._trained._data) == ["z1", "z2"]
        # Assert that correct expression strings were returned for both surrogates
        assert pysmo_surr_rbf._trained._data["z1"].expression_str == (
            "10 + 40*(-69.10791015625*exp(- (0.05*(((IndexedParam[x1] -1)/4)**2 + ((IndexedParam[x2] -5)/4)**2)**0.5)**2) "
            "-319807.1317138672*exp(- (0.05*(((IndexedParam[x1] -1)/4 -0.25)**2 + ((IndexedParam[x2] -5)/4 -0.25)**2)**0.5)**2) "
            "+ 959336.2551269531*exp(- (0.05*(((IndexedParam[x1] -1)/4 -0.5)**2 + ((IndexedParam[x2] -5)/4 -0.5)**2)**0.5)**2) "
            "-959973.7440185547*exp(- (0.05*(((IndexedParam[x1] -1)/4 -0.75)**2 + ((IndexedParam[x2] -5)/4 -0.75)**2)**0.5)**2) "
            "+ 320514.66677856445*exp(- (0.05*(((IndexedParam[x1] -1)/4 -1.0)**2 + ((IndexedParam[x2] -5)/4 -1.0)**2)**0.5)**2))"
        )
        assert pysmo_surr_rbf._trained._data["z2"].expression_str == (
            "6 + 8*(-69.10791015625*exp(- (0.05*(((IndexedParam[x1] -1)/4)**2 + ((IndexedParam[x2] -5)/4)**2)**0.5)**2) "
            "-319807.1317138672*exp(- (0.05*(((IndexedParam[x1] -1)/4 -0.25)**2 + ((IndexedParam[x2] -5)/4 -0.25)**2)**0.5)**2) "
            "+ 959336.2551269531*exp(- (0.05*(((IndexedParam[x1] -1)/4 -0.5)**2 + ((IndexedParam[x2] -5)/4 -0.5)**2)**0.5)**2) "
            "-959973.7440185547*exp(- (0.05*(((IndexedParam[x1] -1)/4 -0.75)**2 + ((IndexedParam[x2] -5)/4 -0.75)**2)**0.5)**2) "
            "+ 320514.66677856445*exp(- (0.05*(((IndexedParam[x1] -1)/4 -1.0)**2 + ((IndexedParam[x2] -5)/4 -1.0)**2)**0.5)**2))"
        )
        # Assert that correct model is returned with generate_expression()
        assert str(
            pysmo_surr_rbf._trained._data["z1"]._model.generate_expression(
                [m.inputs["x1"], m.inputs["x2"]]
            )
        ) == (
            "10 + 40*(-69.10791015625*exp(- (0.05*(((inputs[x1] - 1)/4)**2 + ((inputs[x2] - 5)/4)**2)**0.5)**2) "
            "- 319807.1317138672*exp(- (0.05*(((inputs[x1] - 1)/4 - 0.25)**2 + ((inputs[x2] - 5)/4 - 0.25)**2)**0.5)**2) "
            "+ 959336.2551269531*exp(- (0.05*(((inputs[x1] - 1)/4 - 0.5)**2 + ((inputs[x2] - 5)/4 - 0.5)**2)**0.5)**2) "
            "- 959973.7440185547*exp(- (0.05*(((inputs[x1] - 1)/4 - 0.75)**2 + ((inputs[x2] - 5)/4 - 0.75)**2)**0.5)**2) "
            "+ 320514.66677856445*exp(- (0.05*(((inputs[x1] - 1)/4 - 1.0)**2 + ((inputs[x2] - 5)/4 - 1.0)**2)**0.5)**2))"
        )
        assert str(
            pysmo_surr_rbf._trained._data["z2"]._model.generate_expression(
                [m.inputs["x1"], m.inputs["x2"]]
            )
        ) == (
            "6 + 8*(-69.10791015625*exp(- (0.05*(((inputs[x1] - 1)/4)**2 + ((inputs[x2] - 5)/4)**2)**0.5)**2) "
            "- 319807.1317138672*exp(- (0.05*(((inputs[x1] - 1)/4 - 0.25)**2 + ((inputs[x2] - 5)/4 - 0.25)**2)**0.5)**2) "
            "+ 959336.2551269531*exp(- (0.05*(((inputs[x1] - 1)/4 - 0.5)**2 + ((inputs[x2] - 5)/4 - 0.5)**2)**0.5)**2) "
            "- 959973.7440185547*exp(- (0.05*(((inputs[x1] - 1)/4 - 0.75)**2 + ((inputs[x2] - 5)/4 - 0.75)**2)**0.5)**2) "
            "+ 320514.66677856445*exp(- (0.05*(((inputs[x1] - 1)/4 - 1.0)**2 + ((inputs[x2] - 5)/4 - 1.0)**2)**0.5)**2))"
        )
        # Assert that returned results will evaluate and return correct results
        out = pysmo_surr_rbf.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert pytest.approx(out["z1"][i], rel=1e-6) == (
                (
                    10
                    + 40
                    * (
                        -69.10791015625
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4) ** 2
                                        + ((inputs["x2"][i] - 5) / 4) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                        - 319807.1317138672
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4 - 0.25) ** 2
                                        + ((inputs["x2"][i] - 5) / 4 - 0.25) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                        + 959336.2551269531
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4 - 0.5) ** 2
                                        + ((inputs["x2"][i] - 5) / 4 - 0.5) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                        - 959973.7440185547
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4 - 0.75) ** 2
                                        + ((inputs["x2"][i] - 5) / 4 - 0.75) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                        + 320514.66677856445
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4 - 1.0) ** 2
                                        + ((inputs["x2"][i] - 5) / 4 - 1.0) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                    )
                )
            )
            assert pytest.approx(out["z2"][i], rel=1e-6) == (
                (
                    6
                    + 8
                    * (
                        -69.10791015625
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4) ** 2
                                        + ((inputs["x2"][i] - 5) / 4) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                        - 319807.1317138672
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4 - 0.25) ** 2
                                        + ((inputs["x2"][i] - 5) / 4 - 0.25) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                        + 959336.2551269531
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4 - 0.5) ** 2
                                        + ((inputs["x2"][i] - 5) / 4 - 0.5) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                        - 959973.7440185547
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4 - 0.75) ** 2
                                        + ((inputs["x2"][i] - 5) / 4 - 0.75) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                        + 320514.66677856445
                        * exp(
                            -(
                                (
                                    0.05
                                    * (
                                        ((inputs["x1"][i] - 1) / 4 - 1.0) ** 2
                                        + ((inputs["x2"][i] - 5) / 4 - 1.0) ** 2
                                    )
                                    ** 0.5
                                )
                                ** 2
                            )
                        )
                    )
                )
            )

    @pytest.mark.unit
    def test_load_krg(self):
        # Try for polynomial with multiple outputs and other user-defined functions
        m = ConcreteModel()
        m.inputs = Var(["x1", "x2"])

        # Validation data
        x = [
            -2,
            -1.8,
            -1.6,
            -1.4,
            -1.2,
            -1.0,
            -0.8,
            -0.6,
            -0.4,
            -0.2,
            0.2,
            0.4,
            0.6,
            0.8,
            1.0,
            1.2,
            1.4,
            1.6,
            1.8,
            2.0,
            2.2,
        ]

        inputs = np.array([np.tile(x, len(x)), np.repeat(x, len(x))])
        inputs = pd.DataFrame(inputs.transpose(), columns=["x1", "x2"])

        pysmo_surr_krg = PysmoSurrogate.load(StringIO(jstring_krg))
        # Assert labels, outputs and bounds
        assert pysmo_surr_krg._input_labels == ["x1", "x2"]
        assert pysmo_surr_krg._output_labels == ["z1", "z2"]
        assert pysmo_surr_krg._input_bounds == {"x1": (0, 5), "x2": (0, 10)}
        assert pysmo_surr_krg._trained.model_type == "kriging"
        assert pysmo_surr_krg._trained.output_labels == ["z1", "z2"]
        assert len(pysmo_surr_krg._trained._data) == 2
        assert list(pysmo_surr_krg._trained._data) == ["z1", "z2"]
        # Assert that correct arrays are loaded for 'z1'
        assert pysmo_surr_krg._trained._data["z1"]._model.optimal_mean == np.array(
            [[30.00000000077694]]
        )
        assert pysmo_surr_krg._trained._data["z1"]._model.optimal_variance == np.array(
            [[6503.3113222215325]]
        )
        assert (
            pysmo_surr_krg._trained._data["z1"]._model.regularization_parameter
            == 1.000000000001e-06
        )
        np.testing.assert_array_equal(
            pysmo_surr_krg._trained._data["z1"]._model.optimal_covariance_matrix,
            np.array(
                [
                    [
                        1.000001,
                        0.9982205353479938,
                        0.9929011178300284,
                        0.9840983398813247,
                        0.971905407660152,
                    ],
                    [
                        0.9982205353479938,
                        1.000001,
                        0.9982205353479938,
                        0.9929011178300284,
                        0.9840983398813247,
                    ],
                    [
                        0.9929011178300284,
                        0.9982205353479938,
                        1.000001,
                        0.9982205353479938,
                        0.9929011178300284,
                    ],
                    [
                        0.9840983398813247,
                        0.9929011178300284,
                        0.9982205353479938,
                        1.000001,
                        0.9982205353479938,
                    ],
                    [
                        0.971905407660152,
                        0.9840983398813247,
                        0.9929011178300284,
                        0.9982205353479938,
                        1.000001,
                    ],
                ]
            ),
        )
        np.testing.assert_array_equal(
            pysmo_surr_krg._trained._data["z1"]._model.covariance_matrix_inverse,
            np.array(
                [
                    [
                        108728.9916945844,
                        -240226.85108007095,
                        82932.18571364644,
                        121970.72026795016,
                        -73364.51387189297,
                    ],
                    [
                        -240226.85108202277,
                        589985.9891969847,
                        -341158.67300272395,
                        -130592.8567227173,
                        121970.72027126199,
                    ],
                    [
                        82932.18571952915,
                        -341158.67301448685,
                        516416.75018761755,
                        -341158.6729826693,
                        82932.18570353556,
                    ],
                    [
                        121970.72026201998,
                        -130592.85670691582,
                        -341158.6729945546,
                        589985.9891699858,
                        -240226.8510697507,
                    ],
                    [
                        -73364.51386989365,
                        121970.72026527137,
                        82932.18570954115,
                        -240226.85107176506,
                        108728.99169106234,
                    ],
                ]
            ),
        )
        np.testing.assert_array_equal(
            pysmo_surr_krg._trained._data["z1"]._model.optimal_y_mu,
            np.array(
                [
                    [-20.00000000077694],
                    [-10.00000000077694],
                    [-7.769394017032027e-10],
                    [9.99999999922306],
                    [19.99999999922306],
                ]
            ),
        )
        np.testing.assert_array_equal(
            pysmo_surr_krg._trained._data["z1"]._model.optimal_weights,
            np.array([0.027452451845611077, 0.0010443446337808024]),
        )

        # Assert that correct arrays are loaded for 'z2'
        assert pysmo_surr_krg._trained._data["z2"]._model.optimal_mean == np.array(
            [[9.999999999902883]]
        )
        assert pysmo_surr_krg._trained._data["z2"]._model.optimal_variance == np.array(
            [[260.13320726701056]]
        )
        assert (
            pysmo_surr_krg._trained._data["z2"]._model.regularization_parameter == 1e-06
        )
        np.testing.assert_array_equal(
            pysmo_surr_krg._trained._data["z2"]._model.optimal_covariance_matrix,
            np.array(
                [
                    [
                        1.000001,
                        0.998220543300601,
                        0.9929011494709431,
                        0.9840984104422155,
                        0.9719055315475238,
                    ],
                    [
                        0.998220543300601,
                        1.000001,
                        0.998220543300601,
                        0.9929011494709431,
                        0.9840984104422155,
                    ],
                    [
                        0.9929011494709431,
                        0.998220543300601,
                        1.000001,
                        0.998220543300601,
                        0.9929011494709431,
                    ],
                    [
                        0.9840984104422155,
                        0.9929011494709431,
                        0.998220543300601,
                        1.000001,
                        0.998220543300601,
                    ],
                    [
                        0.9719055315475238,
                        0.9840984104422155,
                        0.9929011494709431,
                        0.998220543300601,
                        1.000001,
                    ],
                ]
            ),
        )
        np.testing.assert_array_equal(
            pysmo_surr_krg._trained._data["z2"]._model.covariance_matrix_inverse,
            np.array(
                [
                    [
                        108729.13455237681,
                        -240227.09704128528,
                        82932.15558036882,
                        121970.94143487987,
                        -73364.601633614,
                    ],
                    [
                        -240227.0970392892,
                        589986.4681472526,
                        -341158.6596781079,
                        -130593.32427863385,
                        121970.94144222786,
                    ],
                    [
                        82932.15557448889,
                        -341158.6596663887,
                        516416.7835787105,
                        -341158.659633822,
                        82932.15555811858,
                    ],
                    [
                        121970.94144067129,
                        -130593.32429416949,
                        -341158.6596220617,
                        589986.4680877628,
                        -240227.09701875152,
                    ],
                    [
                        -73364.60163552182,
                        121970.94144804058,
                        82932.15555219717,
                        -240227.09701673465,
                        108729.13454474375,
                    ],
                ]
            ),
        )
        np.testing.assert_array_equal(
            pysmo_surr_krg._trained._data["z2"]._model.optimal_y_mu,
            np.array(
                [
                    [-3.999999999902883],
                    [-1.999999999902883],
                    [9.711698112369049e-11],
                    [2.000000000097117],
                    [4.000000000097117],
                ]
            ),
        )
        np.testing.assert_array_equal(
            pysmo_surr_krg._trained._data["z2"]._model.optimal_weights,
            np.array([0.02749666901085125, 0.001000000000000049]),
        )

        # Assert that returned results will evaluate and return correct results
        out = pysmo_surr_krg.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert pytest.approx(out["z1"][i], rel=1e-6) == (
                -19894.397849368
                * exp(
                    -(
                        0.027452451845611077 * ((inputs["x1"][i] - 1) / 4) ** 2
                        + 0.0010443446337808024 * ((inputs["x2"][i] - 5) / 4) ** 2
                    )
                )
                + 38162.96786869278
                * exp(
                    -(
                        0.027452451845611077 * ((inputs["x1"][i] - 1) / 4 - 0.25) ** 2
                        + 0.0010443446337808024
                        * ((inputs["x2"][i] - 5) / 4 - 0.25) ** 2
                    )
                )
                - 1.6681948100955743e-06
                * exp(
                    -(
                        0.027452451845611077 * ((inputs["x1"][i] - 1) / 4 - 0.5) ** 2
                        + 0.0010443446337808024 * ((inputs["x2"][i] - 5) / 4 - 0.5) ** 2
                    )
                )
                - 38162.96786638197
                * exp(
                    -(
                        0.027452451845611077 * ((inputs["x1"][i] - 1) / 4 - 0.75) ** 2
                        + 0.0010443446337808024
                        * ((inputs["x2"][i] - 5) / 4 - 0.75) ** 2
                    )
                )
                + 19894.397848724166
                * exp(
                    -(
                        0.027452451845611077 * ((inputs["x1"][i] - 1) / 4 - 1.0) ** 2
                        + 0.0010443446337808024 * ((inputs["x2"][i] - 5) / 4 - 1.0) ** 2
                    )
                )
                + 30.00000000077694
            )
            assert pytest.approx(out["z2"][i], rel=1e-6) == (
                (
                    -3978.867791629029
                    * exp(
                        -(
                            0.02749666901085125 * ((inputs["x1"][i] - 1) / 4) ** 2
                            + 0.001000000000000049 * ((inputs["x2"][i] - 5) / 4) ** 2
                        )
                    )
                    + 7632.569074293324
                    * exp(
                        -(
                            0.02749666901085125
                            * ((inputs["x1"][i] - 1) / 4 - 0.25) ** 2
                            + 0.001000000000000049
                            * ((inputs["x2"][i] - 5) / 4 - 0.25) ** 2
                        )
                    )
                    - 3.5124027300266805e-07
                    * exp(
                        -(
                            0.02749666901085125 * ((inputs["x1"][i] - 1) / 4 - 0.5) ** 2
                            + 0.001000000000000049
                            * ((inputs["x2"][i] - 5) / 4 - 0.5) ** 2
                        )
                    )
                    - 7632.569073828787
                    * exp(
                        -(
                            0.02749666901085125
                            * ((inputs["x1"][i] - 1) / 4 - 0.75) ** 2
                            + 0.001000000000000049
                            * ((inputs["x2"][i] - 5) / 4 - 0.75) ** 2
                        )
                    )
                    + 3978.8677915156522
                    * exp(
                        -(
                            0.02749666901085125 * ((inputs["x1"][i] - 1) / 4 - 1.0) ** 2
                            + 0.001000000000000049
                            * ((inputs["x2"][i] - 5) / 4 - 1.0) ** 2
                        )
                    )
                    + 9.999999999902883
                )
            )

    @pytest.mark.unit
    def test_save_load(self, pysmo_surr1):
        m = ConcreteModel()
        m.inputs = Var(["x1", "x2"])

        # Save and re-load object
        with TempfileManager as tf:
            fname = tf.create_tempfile(suffix=".json")
            pysmo_surr1.save_to_file(fname, overwrite=True)

            assert os.path.isfile(fname)

            with open(fname, "r") as f:
                js = f.read()
            f.close()

            pysmo_load = PysmoSurrogate.load_from_file(fname)

        # Check loaded object
        assert pysmo_load._input_labels == ["x1", "x2"]
        assert pysmo_load._output_labels == ["z1"]
        assert pysmo_load._input_bounds == {"x1": (0, 5), "x2": (0, 10)}
        assert pysmo_load._trained.model_type == "poly"
        assert pysmo_load._trained.output_labels == ["z1"]
        assert len(pysmo_load._trained._data) == 1
        assert list(pysmo_load._trained._data) == ["z1"]
        # Assert that correct expression string was returned
        assert (
            pysmo_load._trained._data["z1"].expression_str
            == pysmo_surr1._trained._data["z1"].expression_str
        )
        # Assert that correct model is returned with generate_expression()
        assert str(
            pysmo_load._trained._data["z1"]._model.generate_expression(
                [m.inputs["x1"], m.inputs["x2"]]
            )
        ) == str(
            pysmo_surr1._trained._data["z1"]._model.generate_expression(
                [m.inputs["x1"], m.inputs["x2"]]
            )
        )
        # Test that `'evaluate_surrogate`` returns same results pre and post saving
        x = [
            -2,
            -1.8,
            -1.6,
            -1.4,
            -1.2,
            -1.0,
            -0.8,
            -0.6,
            -0.4,
            -0.2,
            0,
            0.2,
            0.4,
            0.6,
            0.8,
            1.0,
            1.2,
            1.4,
            1.6,
            1.8,
            2.0,
        ]
        inputs = np.array([np.tile(x, len(x)), np.repeat(x, len(x))])
        inputs = pd.DataFrame(inputs.transpose(), columns=["x1", "x2"])
        out_presave = pysmo_surr1.evaluate_surrogate(inputs)
        out_postsave = pysmo_load.evaluate_surrogate(inputs)
        for i in range(inputs.shape[0]):
            assert (
                pytest.approx(out_presave["z1"][i], rel=1e-8) == out_postsave["z1"][i]
            )

        # Check for clean up
        assert not os.path.isfile(fname)


@pytest.mark.integration
class TestRegressionWorkflow:
    training_data = pd.DataFrame(
        np.array(
            [
                [0.353837234435, 0.99275270941666, 0.762878272854],
                [0.904978848612, -0.746908518721, 0.387963718723],
                [0.643706630938, -0.617496599522, -0.0205375902284],
                [1.29881420688, 0.305594881575, 2.43011137696],
                [1.35791650867, 0.351045058258, 2.36989368612],
                [0.938369314089, -0.525167416293, 0.829756159423],
                [-1.46593541641, 0.383902178482, 1.14054797964],
                [-0.374378293218, -0.689730440659, -0.219122783909],
                [0.690326213554, 0.569364994374, 0.982068847698],
                [-0.961163301329, 0.499471920546, 0.936855365038],
            ]
        ),
        columns=["x1", "x2", "z1"],
    )

    @pytest.fixture(scope="class")
    def pysmo_trainer_poly(self):
        # Test end-to-end workflow with a simple problem.
        bnds = {"x1": (-1.5, 1.5), "x2": (-1.5, 1.5)}
        pysmo_trainer_poly = PysmoPolyTrainer(
            input_labels=["x1", "x2"],
            output_labels=["z1"],
            input_bounds=bnds,
            training_dataframe=TestRegressionWorkflow.training_data,
        )
        pysmo_trainer_poly.config.maximum_polynomial_order = 6
        pysmo_trainer_poly.config.multinomials = True

        res = pysmo_trainer_poly.train_surrogate()

        return res

    def test_pysmo_trainer_poly_results(self, pysmo_trainer_poly):
        assert pysmo_trainer_poly._data is not None
        assert pysmo_trainer_poly.model_type == "poly"
        assert pysmo_trainer_poly.num_outputs == 1
        assert pysmo_trainer_poly.output_labels == ["z1"]
        assert pysmo_trainer_poly.input_labels == ["x1", "x2"]
        assert len(pysmo_trainer_poly._data) == 1

    def test_pysmo_trainer_poly_object(self, pysmo_trainer_poly):
        bnds = {"x1": (-1.5, 1.5), "x2": (-1.5, 1.5)}
        input_labels = ["x1", "x2"]
        output_labels = ["z1"]
        input_bounds = bnds
        pysmo_object = PysmoSurrogate(
            pysmo_trainer_poly, input_labels, output_labels, bnds
        )
        assert isinstance(pysmo_object, PysmoSurrogate)
        assert pysmo_object._input_labels == ["x1", "x2"]
        assert pysmo_object._output_labels == ["z1"]
        assert pysmo_object._input_bounds == {"x1": (-1.5, 1.5), "x2": (-1.5, 1.5)}

        # Check populating a block to finish workflow
        blk = SurrogateBlock(concrete=True)

        blk.build_model(pysmo_object)

        assert isinstance(blk.inputs, Var)
        assert blk.inputs["x1"].bounds == (-1.5, 1.5)
        assert blk.inputs["x2"].bounds == (-1.5, 1.5)
        assert isinstance(blk.outputs, Var)
        assert blk.outputs["z1"].bounds == (None, None)
        assert isinstance(blk.pysmo_constraint, Constraint)
        assert len(blk.pysmo_constraint) == 1

    def test_metrics(self, pysmo_trainer_poly):
        bnds = {"x1": (-1.5, 1.5), "x2": (-1.5, 1.5)}
        input_labels = ["x1", "x2"]
        output_labels = ["z1"]
        input_bounds = bnds
        pysmo_object = PysmoSurrogate(
            pysmo_trainer_poly, input_labels, output_labels, bnds
        )

        metrics = compute_fit_metrics(
            pysmo_object, TestRegressionWorkflow.training_data
        )

        assert isinstance(metrics, dict)

        assert metrics["z1"]["R2"] == pytest.approx(
            pysmo_trainer_poly._data["z1"].metrics["R2"], rel=1e-8
        )
        assert metrics["z1"]["RMSE"] == pytest.approx(
            pysmo_trainer_poly._data["z1"].metrics["RMSE"], rel=1e-8
        )


@pytest.mark.integration
class TestRBFWorkflow:
    training_data = pd.DataFrame(
        np.array(
            [
                [0.353837234435, 0.99275270941666, 0.762878272854],
                [0.904978848612, -0.746908518721, 0.387963718723],
                [0.643706630938, -0.617496599522, -0.0205375902284],
                [1.29881420688, 0.305594881575, 2.43011137696],
                [1.35791650867, 0.351045058258, 2.36989368612],
                [0.938369314089, -0.525167416293, 0.829756159423],
                [-1.46593541641, 0.383902178482, 1.14054797964],
                [-0.374378293218, -0.689730440659, -0.219122783909],
                [0.690326213554, 0.569364994374, 0.982068847698],
                [-0.961163301329, 0.499471920546, 0.936855365038],
            ]
        ),
        columns=["x1", "x2", "z1"],
    )

    @pytest.fixture(scope="class")
    def pysmo_trainer_rbf(self):
        # Test end-to-end workflow with a simple problem.
        bnds = {"x1": (-1.5, 1.5), "x2": (-1.5, 1.5)}
        pysmo_trainer_rbf = PysmoRBFTrainer(
            input_labels=["x1", "x2"],
            output_labels=["z1"],
            input_bounds=bnds,
            training_dataframe=TestRBFWorkflow.training_data,
        )
        pysmo_trainer_rbf.config.regularization = False
        pysmo_trainer_rbf.config.basis_function = "cubic"

        res = pysmo_trainer_rbf.train_surrogate()

        return res

    def test_pysmo_trainer_rbf_results(self, pysmo_trainer_rbf):
        assert pysmo_trainer_rbf._data is not None
        assert pysmo_trainer_rbf.model_type == "rbf"
        assert pysmo_trainer_rbf.num_outputs == 1
        assert pysmo_trainer_rbf.output_labels == ["z1"]
        assert pysmo_trainer_rbf.input_labels == ["x1", "x2"]
        assert len(pysmo_trainer_rbf._data) == 1

    def test_pysmo_trainer_rbf_object(self, pysmo_trainer_rbf):
        bnds = {"x1": (-1.5, 1.5), "x2": (-1.5, 1.5)}
        input_labels = ["x1", "x2"]
        output_labels = ["z1"]
        input_bounds = bnds
        pysmo_object = PysmoSurrogate(
            pysmo_trainer_rbf, input_labels, output_labels, bnds
        )
        assert isinstance(pysmo_object, PysmoSurrogate)
        assert pysmo_object._input_labels == ["x1", "x2"]
        assert pysmo_object._output_labels == ["z1"]
        assert pysmo_object._input_bounds == {"x1": (-1.5, 1.5), "x2": (-1.5, 1.5)}

        # Check populating a block to finish workflow
        blk = SurrogateBlock(concrete=True)

        blk.build_model(pysmo_object)

        assert isinstance(blk.inputs, Var)
        assert blk.inputs["x1"].bounds == (-1.5, 1.5)
        assert blk.inputs["x2"].bounds == (-1.5, 1.5)
        assert isinstance(blk.outputs, Var)
        assert blk.outputs["z1"].bounds == (None, None)
        assert isinstance(blk.pysmo_constraint, Constraint)
        assert len(blk.pysmo_constraint) == 1

    def test_metrics(self, pysmo_trainer_rbf):
        bnds = {"x1": (-1.5, 1.5), "x2": (-1.5, 1.5)}
        input_labels = ["x1", "x2"]
        output_labels = ["z1"]
        input_bounds = bnds
        pysmo_object = PysmoSurrogate(
            pysmo_trainer_rbf, input_labels, output_labels, bnds
        )

        metrics = compute_fit_metrics(pysmo_object, TestRBFWorkflow.training_data)

        assert isinstance(metrics, dict)

        assert metrics["z1"]["R2"] == pytest.approx(
            pysmo_trainer_rbf._data["z1"].metrics["R2"], rel=1e-8
        )
        assert metrics["z1"]["RMSE"] == pytest.approx(
            pysmo_trainer_rbf._data["z1"].metrics["RMSE"], rel=1e-8
        )


@pytest.mark.integration
class TestKrigingWorkflow:
    training_data = pd.DataFrame(
        np.array(
            [
                [0.353837234435, 0.99275270941666, 0.762878272854],
                [0.904978848612, -0.746908518721, 0.387963718723],
                [0.643706630938, -0.617496599522, -0.0205375902284],
                [1.29881420688, 0.305594881575, 2.43011137696],
                [1.35791650867, 0.351045058258, 2.36989368612],
                [0.938369314089, -0.525167416293, 0.829756159423],
                [-1.46593541641, 0.383902178482, 1.14054797964],
                [-0.374378293218, -0.689730440659, -0.219122783909],
                [0.690326213554, 0.569364994374, 0.982068847698],
                [-0.961163301329, 0.499471920546, 0.936855365038],
            ]
        ),
        columns=["x1", "x2", "z1"],
    )

    @pytest.fixture(scope="class")
    def pysmo_trainer_krg(self):
        # Test end-to-end workflow with a simple problem.
        bnds = {"x1": (-1.5, 1.5), "x2": (-1.5, 1.5)}
        pysmo_trainer_krg = PysmoKrigingTrainer(
            input_labels=["x1", "x2"],
            output_labels=["z1"],
            input_bounds=bnds,
            training_dataframe=TestKrigingWorkflow.training_data,
        )
        pysmo_trainer_krg.config.regularization = False
        pysmo_trainer_krg.config.numerical_gradients = False

        np.random.seed(0)
        res = pysmo_trainer_krg.train_surrogate()

        return res

    def test_pysmo_trainer_krg_results(self, pysmo_trainer_krg):
        assert pysmo_trainer_krg._data is not None
        assert pysmo_trainer_krg.model_type == "kriging"
        assert pysmo_trainer_krg.num_outputs == 1
        assert pysmo_trainer_krg.output_labels == ["z1"]
        assert pysmo_trainer_krg.input_labels == ["x1", "x2"]
        assert len(pysmo_trainer_krg._data) == 1

    def test_pysmo_trainer_krg_object(self, pysmo_trainer_krg):
        bnds = {"x1": (-1.5, 1.5), "x2": (-1.5, 1.5)}
        input_labels = ["x1", "x2"]
        output_labels = ["z1"]
        input_bounds = bnds
        pysmo_object = PysmoSurrogate(
            pysmo_trainer_krg, input_labels, output_labels, bnds
        )
        assert isinstance(pysmo_object, PysmoSurrogate)
        assert pysmo_object._input_labels == ["x1", "x2"]
        assert pysmo_object._output_labels == ["z1"]
        assert pysmo_object._input_bounds == {"x1": (-1.5, 1.5), "x2": (-1.5, 1.5)}

        # Check populating a block to finish workflow
        blk = SurrogateBlock(concrete=True)

        blk.build_model(pysmo_object)

        assert isinstance(blk.inputs, Var)
        assert blk.inputs["x1"].bounds == (-1.5, 1.5)
        assert blk.inputs["x2"].bounds == (-1.5, 1.5)
        assert isinstance(blk.outputs, Var)
        assert blk.outputs["z1"].bounds == (None, None)
        assert isinstance(blk.pysmo_constraint, Constraint)
        assert len(blk.pysmo_constraint) == 1

    def test_metrics(self, pysmo_trainer_krg):
        bnds = {"x1": (-1.5, 1.5), "x2": (-1.5, 1.5)}
        input_labels = ["x1", "x2"]
        output_labels = ["z1"]
        input_bounds = bnds
        pysmo_object = PysmoSurrogate(
            pysmo_trainer_krg, input_labels, output_labels, bnds
        )

        metrics = compute_fit_metrics(pysmo_object, TestRBFWorkflow.training_data)

        assert isinstance(metrics, dict)

        assert metrics["z1"]["R2"] == pytest.approx(
            pysmo_trainer_krg._data["z1"].metrics["R2"], rel=1e-8
        )
        assert metrics["z1"]["RMSE"] == pytest.approx(
            pysmo_trainer_krg._data["z1"].metrics["RMSE"], rel=1e-8
        )
