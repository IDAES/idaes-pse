#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for KerasSurrogate
"""
import pytest

pytest.importorskip("lineartree", reason="lineartree not available")
pytest.importorskip("omlt", reason="omlt not available")

import os.path
import pandas as pd
import numpy as np
from pyomo.common.fileutils import this_file_dir
from pyomo.common.tempfiles import TempfileManager
from pyomo.environ import (
    ConcreteModel,
    Var,
    SolverFactory,
    assert_optimal_termination,
    value,
    Objective,
)
from idaes.core.surrogate.linear_tree_surrogate import LinearTreeSurrogate, load_linear_tree_pickle
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.surrogate.sampling.scaling import OffsetScaler

from lineartree import LinearTreeRegressor
from sklearn.linear_model import LinearRegression


rtol = 1e-4
atol = 1e-4

X_small = np.array(
[
    [-0.68984135],
    [0.91672866],
    [-1.05874972],
    [0.95275351],
    [1.03796615],
    [0.45117668],
    [-0.14704376],
    [1.66043409],
    [-0.73972191],
    [-0.8176603],
    [0.96175973],
    [-1.238874],
    [-0.97492265],
    [1.07121986],
    [-0.95379269],
    [-0.86546252],
    [0.8277057],
    [0.50486757],
    [-1.38435899],
    [1.54092856],
]
)

y_small = np.array(
    [
        [0.04296633],
        [-0.78349216],
        [0.27114188],
        [-0.58516476],
        [-0.15997756],
        [-0.37529212],
        [-1.49249696],
        [1.56412122],
        [0.18697725],
        [0.4035928],
        [-0.53231771],
        [-0.02669967],
        [0.36972983],
        [0.09201347],
        [0.44041505],
        [0.46047019],
        [-1.04855941],
        [-0.586915],
        [0.15472157],
        [1.71225268],
    ]
)

def linear_model_tree(X, y):
    regr = LinearTreeRegressor(LinearRegression(), criterion="mse", max_depth=5)
    regr.fit(X, y)
    return regr

def create_lt_model(return_lt_model_only=True):

    mean_x_small = np.mean(X_small)
    std_x_small = np.std(X_small)
    mean_y_small = np.mean(y_small)
    std_y_small = np.std(y_small)
    scaled_x = (X_small - mean_x_small) / std_x_small
    scaled_y = (y_small - mean_y_small) / std_y_small
    scaled_input_bounds = {"X": (np.min(scaled_x), np.max(scaled_x))}
    unscaled_input_bounds = {"X": (np.min(X_small), np.max(X_small))}

    lt_model = linear_model_tree(scaled_x, scaled_y)
    if return_lt_model_only:
        return lt_model

    inputs_scaler = OffsetScaler(
    expected_columns=["X"],
    offset_series=pd.Series({"X": mean_x_small}),
    factor_series=pd.Series({"X": std_x_small}),
    )

    outputs_scaler = OffsetScaler(
        expected_columns=["Y"],
        offset_series=pd.Series({"Y": mean_y_small}),
        factor_series=pd.Series({"Y": std_y_small}),
    )
    input_labels = ["X"]
    output_labels = ["Y"]

    lt_surrogate = LinearTreeSurrogate(
        lt_model=lt_model,
        input_labels=input_labels,
        output_labels=output_labels,
        input_bounds=unscaled_input_bounds,
        input_scaler=inputs_scaler,
        output_scaler=outputs_scaler,
    )

    return lt_surrogate


@pytest.mark.unit
def test_LinearTreeSurrogate_construction_exceptions():
    lt_model = create_lt_model()
    inputs_scaler = OffsetScaler(
        expected_columns=["X"],
        offset_series=pd.Series({"X": 0.0528}),
        factor_series=pd.Series({"X": 1.00003}),
    )
    outputs_scaler = OffsetScaler(
        expected_columns=["Y"],
        offset_series=pd.Series({"Y": 0.0054}),
        factor_series=pd.Series({"Y": 0.7518}),
    )

    # test exceptions for mismatched labels
    with pytest.raises(ValueError) as excinfo:
        lt_surrogate = LinearTreeSurrogate(
            lt_model=lt_model,
            input_labels=["X_wrong"],
            output_labels=["Y"],
            input_bounds={"X_wrong": (-1.3844, 1.6604)},
            input_scaler=inputs_scaler,
            output_scaler=outputs_scaler,
        )
    assert (
        str(excinfo.value)
        == "OMLTSurrogate created with input_labels that do not match the expected columns in the input_scaler.\n"
        "input_labels=['X_wrong']\n"
        "input_scaler.expected_columns()=['X']"
    )

    with pytest.raises(ValueError) as excinfo:
        lt_surrogate = LinearTreeSurrogate(
            lt_model=lt_model,
            input_labels=["X"],
            output_labels=["Y_wrong"],
            input_bounds={"X": (-1.3844, 1.6604)},
            input_scaler=inputs_scaler,
            output_scaler=outputs_scaler,
        )
    assert (
        str(excinfo.value)
        == "OMLTSurrogate created with output_labels that do not match the expected columns in the output_scaler.\n"
        "output_labels=['Y_wrong']\n"
        "output_scaler.expected_columns()=['Y']"
    )

    with pytest.raises(ValueError) as excinfo:
        lt_surrogate = LinearTreeSurrogate(
            lt_model=lt_model,
            input_labels=["X"],
            output_labels=["Y"],
            input_bounds={"X_wrong": (-1.3844, 1.6604)},
            input_scaler=inputs_scaler,
            output_scaler=outputs_scaler,
        )
    assert (
        str(excinfo.value)
        == "The input_labels did not match the keys in input_bounds.\n"
        "input_bounds.keys(): ['X_wrong']\n"
        "input_labels: ['X']"
    )


@pytest.mark.unit
def test_lt_evaluate():
    x = pd.DataFrame({"X": [0.5]})
    lt_surrogate = create_lt_model(
        return_lt_model_only=False
    )
    lt_mod = create_lt_model(
        return_lt_model_only=True
    )
    y = lt_surrogate.evaluate_surrogate(x)
    expected_y = lt_mod.predict(np.array([(0.5 - np.mean(X_small))/np.std(X_small)]).reshape(-1,1))
    assert expected_y[0] == pytest.approx((y.to_numpy()[0][0] - np.mean(y_small))/np.std(y_small))

@pytest.mark.unit
@pytest.mark.skipif(not SolverFactory("glpk").available(False), reason="no glpk")
def test_keras_surrogate_auto_creating_variables_glpk():
    ###
    # Test Linear Tree
    ###
    lt_surrogate = create_lt_model(
        return_lt_model_only=False
    )
    x_test = pd.DataFrame({"X": [0.5]})
    y_test = lt_surrogate.evaluate_surrogate(x_test)

    # Test lineartree bigm
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=lt_surrogate,
        formulation=LinearTreeSurrogate.Formulation.LINEAR_TREE_GDP_BIGM,
    )
    m.surrogate.inputs["X"].fix(0.5)
    solver = SolverFactory("glpk")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {
            "Y": [value(m.surrogate.outputs["Y"])],
        }
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)

    # Test lineartree hull
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=lt_surrogate,
        formulation=LinearTreeSurrogate.Formulation.LINEAR_TREE_GDP_HULL,
    )
    m.surrogate.inputs["X"].fix(0.5)
    solver = SolverFactory("glpk")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {
            "Y": [value(m.surrogate.outputs["Y"])],
        }
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)

@pytest.mark.unit
@pytest.mark.skipif(not SolverFactory("scip").available(False), reason="no scip")
def test_lt_surrogate_auto_creating_variables_scip():
    ###
    # Test 1->2 relu
    ###
    lt_surrogate = create_lt_model(
        return_lt_model_only=False
    )
    x_test = pd.DataFrame({"X": [0.5]})
    y_test = lt_surrogate.evaluate_surrogate(x_test)

    # Test lineartree bigm
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=lt_surrogate,
        formulation=LinearTreeSurrogate.Formulation.LINEAR_TREE_HYBRID_BIGM,
    )
    m.surrogate.inputs["X"].fix(0.5)
    solver = SolverFactory("scip")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {
            "Y": [value(m.surrogate.outputs["Y"])],
        }
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)


@pytest.mark.unit
@pytest.mark.skipif(not SolverFactory("glpk").available(False), reason="no glpk")
def test_lt_surrogate_with_variables():
    lt_surrogate = create_lt_model(
        return_lt_model_only=False
    )
    x_test = pd.DataFrame({"X": [0.5]})
    y_test = lt_surrogate.evaluate_surrogate(x_test)

    # Test provide scalar inputs, auto create outputs
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.X = Var()
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=lt_surrogate,
        input_vars=[m.X],
        formulation=LinearTreeSurrogate.Formulation.LINEAR_TREE_GDP_BIGM,
    )
    m.X.fix(0.5)
    solver = SolverFactory("glpk")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {
            "Y": [value(m.surrogate.outputs["Y"])],
        }
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)

    # Test provide indexed inputs, auto create outputs
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.inputs = Var(["X"])
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=lt_surrogate,
        input_vars=[m.inputs],
        formulation=LinearTreeSurrogate.Formulation.LINEAR_TREE_GDP_BIGM,
    )
    m.inputs["X"].fix(0.5)
    solver = SolverFactory("glpk")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {
            "Y": [value(m.surrogate.outputs["Y"])],
        }
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)

    # Test auto-create inputs, provide scalar outputs
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.Y = Var()
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=lt_surrogate,
        output_vars=[m.Y],
        formulation=LinearTreeSurrogate.Formulation.LINEAR_TREE_GDP_BIGM,
    )
    m.surrogate.inputs["X"].fix(0.5)
    solver = SolverFactory("glpk")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame({"Y": [value(m.Y)]})
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)

    # Test auto-create inputs, provide indexed outputs
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.outputs = Var(["Y"])
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=lt_surrogate,
        output_vars=[m.outputs],
        formulation=LinearTreeSurrogate.Formulation.LINEAR_TREE_GDP_BIGM,
    )
    m.surrogate.inputs["X"].fix(0.5)
    solver = SolverFactory("glpk")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {"Y": [value(m.outputs["Y"])]}
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)

    # Test provide scalar inputs, provide indexed outputs
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.X = Var()
    m.outputs = Var(["Y"])
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=lt_surrogate,
        input_vars=[m.X],
        output_vars=[m.outputs],
        formulation=LinearTreeSurrogate.Formulation.LINEAR_TREE_GDP_BIGM,
    )
    m.X.fix(0.5)
    solver = SolverFactory("glpk")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {"Y": [value(m.outputs["Y"])]}
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)


@pytest.mark.unit
@pytest.mark.skipif(not SolverFactory("glpk").available(False), reason="no glpk")
def test_save_load():
    lt_surrogate = create_lt_model(
        return_lt_model_only=False
    )
    lt_mod = create_lt_model()
    new_lt_surrogate = None
    dname = None
    with TempfileManager.new_context() as tf:
        dname = tf.mkdtemp()
        lt_surrogate.save_to_folder(dname)
        assert os.path.isdir(dname)
        assert os.path.isfile(os.path.join(dname, "idaes_info.json"))

        new_lt_surrogate = LinearTreeSurrogate.load_from_folder(dname)

    # Check for clean up
    assert not os.path.isdir(dname)

    # check surrogate data members
    assert new_lt_surrogate._input_labels == ["X"]
    assert new_lt_surrogate._output_labels == ["Y"]
    assert sorted(new_lt_surrogate._input_bounds.keys()) == ["X"]
    assert new_lt_surrogate._input_bounds["X"][0] == pytest.approx(-1.38435899)
    assert new_lt_surrogate._input_bounds["X"][1] == pytest.approx(1.66043409)

    # check input scaler
    expected_columns = ["X"]
    offset_series = pd.Series(
        {"X": np.mean(X_small)}
    )
    factor_series = pd.Series({"X": np.std(X_small)})
    scaler = new_lt_surrogate._input_scaler
    assert scaler._expected_columns == expected_columns
    pd.testing.assert_series_equal(scaler._offset, offset_series, rtol=rtol, atol=atol)
    pd.testing.assert_series_equal(scaler._factor, factor_series, rtol=rtol, atol=atol)

    # check output scaler
    expected_columns = ["Y"]
    offset_series = pd.Series({"Y": np.mean(y_small)})
    factor_series = pd.Series({"Y": np.std(y_small)})
    scaler = new_lt_surrogate._output_scaler
    assert scaler._expected_columns == expected_columns
    pd.testing.assert_series_equal(scaler._offset, offset_series, rtol=rtol, atol=atol)
    pd.testing.assert_series_equal(scaler._factor, factor_series, rtol=rtol, atol=atol)

    # check evaluation
    x_test = pd.DataFrame({"X": [0.5]})
    y_test = lt_surrogate.evaluate_surrogate(x_test)
    expected_y = lt_mod.predict(np.array([(0.5 - np.mean(X_small))/np.std(X_small)]).reshape(-1,1))
    assert expected_y[0] == pytest.approx((y_test.to_numpy()[0][0] - np.mean(y_small))/np.std(y_small))

    # check solve with pyomo
    x_test = pd.DataFrame({"X": [0.5]})
    y_test = lt_surrogate.evaluate_surrogate(x_test)

    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=lt_surrogate,
        formulation=LinearTreeSurrogate.Formulation.LINEAR_TREE_GDP_BIGM,
    )
    m.surrogate.inputs["X"].fix(0.5)
    solver = SolverFactory("glpk")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {
            "Y": [value(m.surrogate.outputs["Y"])],
        }
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)


@pytest.mark.unit
@pytest.mark.skipif(not SolverFactory("glpk").available(False), reason="no glpk")
def test_noscalers():
    lt_model = create_lt_model()

    input_labels = ["X"]
    output_labels = ["Y"]
    input_bounds = {"X": (np.min(X_small), np.max(X_small))}

    lt_surrogate = LinearTreeSurrogate(
        lt_model=lt_model,
        input_labels=input_labels,
        output_labels=output_labels,
        input_bounds=input_bounds,
    )
    # check solve with pyomo
    x_test = pd.DataFrame({"X": [0.5]})
    y_test = lt_surrogate.evaluate_surrogate(x_test)

    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=lt_surrogate,
        formulation=LinearTreeSurrogate.Formulation.LINEAR_TREE_GDP_BIGM,
    )
    m.surrogate.inputs["X"].fix(0.5)
    solver = SolverFactory("glpk")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {
            "Y": [value(m.surrogate.outputs["Y"])],
        }
    )
    pd.testing.assert_frame_equal(
        y_test, y_test_pyomo, check_dtype=False, rtol=rtol, atol=atol
    )


@pytest.mark.unit
def test_invalid_formulation():
    lt_surrogate = create_lt_model(
        return_lt_model_only=False
    )
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.surrogate = SurrogateBlock()
    with pytest.raises(ValueError) as excinfo:
        m.surrogate.build_model(surrogate_object=lt_surrogate, formulation="foo")
    assert (
        str(excinfo.value) == 'An unrecognized formulation "foo" was passed '
        "to LinearTreeSurrogate.populate_block. Please pass a valid formulation."
    )
