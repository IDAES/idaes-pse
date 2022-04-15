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
Tests for KerasSurrogate
"""
import pytest

pytest.importorskip("tensorflow.keras", reason="tensorflow.keras not available")
pytest.importorskip("omlt", reason="omlt not available")

import os.path
import pandas as pd
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
from idaes.core.surrogate.keras_surrogate import KerasSurrogate, load_keras_json_hd5
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.surrogate.sampling.scaling import OffsetScaler


rtol = 1e-4
atol = 1e-4


def create_keras_model(name="T_data_1_10_10_2_sigmoid", return_keras_model_only=True):
    #    keras_folder_name = os.path.join(this_file_dir(), 'data', 'keras_models', name)
    #    keras_model = keras.models.load_model(keras_folder_name)
    keras_folder_name = os.path.join(this_file_dir(), "data", "keras_models")
    keras_model = load_keras_json_hd5(keras_folder_name, name)
    if return_keras_model_only:
        return keras_model

    if name.startswith("T_data"):
        inputs_scaler = OffsetScaler(
            expected_columns=["Temperature_K"],
            offset_series=pd.Series({"Temperature_K": 369.873737}),
            factor_series=pd.Series({"Temperature_K": 5.776343}),
        )
        outputs_scaler = OffsetScaler(
            expected_columns=["EnthMol", "VapFrac"],
            offset_series=pd.Series({"EnthMol": 59519.621464, "VapFrac": 0.556033}),
            factor_series=pd.Series({"EnthMol": 14805.370235, "VapFrac": 0.432934}),
        )
        input_labels = ["Temperature_K"]
        output_labels = ["EnthMol", "VapFrac"]
        input_bounds = {"Temperature_K": (360, 380)}
    elif name.startswith("PT_data"):
        inputs_scaler = OffsetScaler(
            expected_columns=["Temperature_K", "Pressure_Pa"],
            offset_series=pd.Series(
                {"Temperature_K": 369.983611, "Pressure_Pa": 111421.319811}
            ),
            factor_series=pd.Series(
                {"Temperature_K": 5.836047, "Pressure_Pa": 5917.954504}
            ),
        )
        outputs_scaler = OffsetScaler(
            expected_columns=["EnthMol", "VapFrac"],
            offset_series=pd.Series({"EnthMol": 54599.629980, "VapFrac": 0.403307}),
            factor_series=pd.Series({"EnthMol": 14654.226615, "VapFrac": 0.430181}),
        )
        input_labels = ["Temperature_K", "Pressure_Pa"]
        output_labels = ["EnthMol", "VapFrac"]
        input_bounds = {
            "Temperature_K": (360.0, 380.0),
            "Pressure_Pa": (101325.0, 1.2 * 101325.0),
        }

    keras_surrogate = KerasSurrogate(
        keras_model=keras_model,
        input_labels=input_labels,
        output_labels=output_labels,
        input_bounds=input_bounds,
        input_scaler=inputs_scaler,
        output_scaler=outputs_scaler,
    )

    return keras_surrogate


@pytest.mark.unit
def test_KerasSurrogate_construction_exceptions():
    keras_model = create_keras_model(name="T_data_1_10_10_2_sigmoid")
    inputs_scaler = OffsetScaler(
        expected_columns=["Temperature_K"],
        offset_series=pd.Series({"Temperature_K": 369.977879}),
        factor_series=pd.Series({"Temperature_K": 5.833396}),
    )
    outputs_scaler = OffsetScaler(
        expected_columns=["EnthMol", "VapFrac"],
        offset_series=pd.Series({"EnthMol": 54582.841444, "VapFrac": 0.402814}),
        factor_series=pd.Series({"EnthMol": 14642.206469, "VapFrac": 0.429829}),
    )

    # test exceptions for mismatched labels
    with pytest.raises(ValueError) as excinfo:
        keras_surrogate = KerasSurrogate(
            keras_model=keras_model,
            input_labels=["Temperature"],
            output_labels=["EnthMol", "VapFrac"],
            input_bounds={"Temperature": (360, 380)},
            input_scaler=inputs_scaler,
            output_scaler=outputs_scaler,
        )
    assert (
        str(excinfo.value)
        == "KerasSurrogate created with input_labels that do not match the expected columns in the input_scaler.\n"
        "input_labels=['Temperature']\n"
        "input_scaler.expected_columns()=['Temperature_K']"
    )

    with pytest.raises(ValueError) as excinfo:
        keras_surrogate = KerasSurrogate(
            keras_model=keras_model,
            input_labels=["Temperature_K"],
            output_labels=["EnthMol", "VapFraction"],
            input_bounds={"Temperature_K": (360, 380)},
            input_scaler=inputs_scaler,
            output_scaler=outputs_scaler,
        )
    assert (
        str(excinfo.value)
        == "KerasSurrogate created with output_labels that do not match the expected columns in the output_scaler.\n"
        "output_labels=['EnthMol', 'VapFraction']\n"
        "output_scaler.expected_columns()=['EnthMol', 'VapFrac']"
    )

    with pytest.raises(ValueError) as excinfo:
        keras_surrogate = KerasSurrogate(
            keras_model=keras_model,
            input_labels=["Temperature_K"],
            output_labels=["EnthMol", "VapFraction"],
            input_bounds={"Temperature": (360, 380)},
            input_scaler=inputs_scaler,
            output_scaler=outputs_scaler,
        )
    assert (
        str(excinfo.value)
        == "The input_labels did not match the keys in input_bounds.\n"
        "input_bounds.keys(): ['Temperature']\n"
        "input_labels: ['Temperature_K']"
    )


@pytest.mark.unit
def test_keras_evaluate():
    x = pd.DataFrame({"Temperature_K": [365, 370, 375]})
    keras_surrogate = create_keras_model(
        name="T_data_1_10_10_2_sigmoid", return_keras_model_only=False
    )
    y = keras_surrogate.evaluate_surrogate(x)
    expected_y = pd.DataFrame(
        {
            "EnthMol": [42686.447275464394, 64437.55629480346, 74746.81533134825],
            "VapFrac": [0.07522893986296653, 0.6880093482449651, 0.9981742834334373],
        }
    )
    pd.testing.assert_frame_equal(y, expected_y, rtol=rtol, atol=atol)

    keras_surrogate = create_keras_model(
        name="T_data_1_10_10_2_relu", return_keras_model_only=False
    )
    y = keras_surrogate.evaluate_surrogate(x)
    expected_y = pd.DataFrame(
        {
            "EnthMol": [40840.90843482576, 63775.05649622617, 74718.54808966955],
            "VapFrac": [-0.001302339421510701, 0.6879003097360135, 0.9980593485100269],
        }
    )
    pd.testing.assert_frame_equal(y, expected_y, rtol=rtol, atol=atol)

    x = pd.DataFrame(
        {
            "Temperature_K": [360, 370, 380],
            "Pressure_Pa": [1.05 * 101325, 1.10 * 101325, 1.15 * 101325],
        }
    )
    keras_surrogate = create_keras_model(
        name="PT_data_2_10_10_2_sigmoid", return_keras_model_only=False
    )
    y = keras_surrogate.evaluate_surrogate(x)
    expected_y = pd.DataFrame(
        {
            "EnthMol": [40194.5586954288, 48660.288218426984, 75178.30324367314],
            "VapFrac": [0.002291496299564877, 0.21942246438431742, 0.9996716243380308],
        }
    )
    pd.testing.assert_frame_equal(y, expected_y, rtol=rtol, atol=atol)

    keras_surrogate = create_keras_model(
        name="PT_data_2_10_10_2_relu", return_keras_model_only=False
    )
    y = keras_surrogate.evaluate_surrogate(x)
    expected_y = pd.DataFrame(
        {
            "EnthMol": [39989.02657637386, 48329.89586985675, 75212.99707375483],
            "VapFrac": [-0.0028634298195242547, 0.2095949409658313, 0.9991636803734303],
        }
    )
    pd.testing.assert_frame_equal(y, expected_y, rtol=rtol, atol=atol)


@pytest.mark.unit
@pytest.mark.skipif(not SolverFactory("ipopt").available(False), reason="no Ipopt")
def test_keras_surrogate_auto_creating_variables():
    ###
    # Test 1->2 sigmoid
    ###
    keras_surrogate = create_keras_model(
        name="T_data_1_10_10_2_sigmoid", return_keras_model_only=False
    )
    x_test = pd.DataFrame({"Temperature_K": [370]})
    y_test = keras_surrogate.evaluate_surrogate(x_test)

    # Test full-space
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=keras_surrogate,
        formulation=KerasSurrogate.Formulation.FULL_SPACE,
    )
    m.surrogate.inputs["Temperature_K"].fix(370)
    solver = SolverFactory("ipopt")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {
            "EnthMol": [value(m.surrogate.outputs["EnthMol"])],
            "VapFrac": [value(m.surrogate.outputs["VapFrac"])],
        }
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)

    # Test reduced-space
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=keras_surrogate,
        formulation=KerasSurrogate.Formulation.REDUCED_SPACE,
    )
    m.surrogate.inputs["Temperature_K"].fix(370)
    solver = SolverFactory("ipopt")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {
            "EnthMol": [value(m.surrogate.outputs["EnthMol"])],
            "VapFrac": [value(m.surrogate.outputs["VapFrac"])],
        }
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)

    ###
    # Test 1->2 relu
    ###
    keras_surrogate = create_keras_model(
        name="T_data_1_10_10_2_relu", return_keras_model_only=False
    )
    x_test = pd.DataFrame({"Temperature_K": [370]})
    y_test = keras_surrogate.evaluate_surrogate(x_test)

    # Test relu complementarity
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=keras_surrogate,
        formulation=KerasSurrogate.Formulation.RELU_COMPLEMENTARITY,
    )
    m.surrogate.inputs["Temperature_K"].fix(370)
    solver = SolverFactory("ipopt")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {
            "EnthMol": [value(m.surrogate.outputs["EnthMol"])],
            "VapFrac": [value(m.surrogate.outputs["VapFrac"])],
        }
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)

    ###
    # Test 2->2 sigmoid
    ###
    keras_surrogate = create_keras_model(
        name="PT_data_2_10_10_2_sigmoid", return_keras_model_only=False
    )
    x_test = pd.DataFrame({"Temperature_K": [370], "Pressure_Pa": [1.1 * 101325]})
    y_test = keras_surrogate.evaluate_surrogate(x_test)

    # Test full-space
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=keras_surrogate,
        formulation=KerasSurrogate.Formulation.FULL_SPACE,
    )
    m.surrogate.inputs["Temperature_K"].fix(370)
    m.surrogate.inputs["Pressure_Pa"].fix(1.1 * 101325)
    solver = SolverFactory("ipopt")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {
            "EnthMol": [value(m.surrogate.outputs["EnthMol"])],
            "VapFrac": [value(m.surrogate.outputs["VapFrac"])],
        }
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)

    # Test reduced-space
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=keras_surrogate,
        formulation=KerasSurrogate.Formulation.REDUCED_SPACE,
    )
    m.surrogate.inputs["Temperature_K"].fix(370)
    m.surrogate.inputs["Pressure_Pa"].fix(1.1 * 101325)
    solver = SolverFactory("ipopt")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {
            "EnthMol": [value(m.surrogate.outputs["EnthMol"])],
            "VapFrac": [value(m.surrogate.outputs["VapFrac"])],
        }
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)

    ###
    # Test 2->2 relu
    ###
    keras_surrogate = create_keras_model(
        name="PT_data_2_10_10_2_relu", return_keras_model_only=False
    )
    x_test = pd.DataFrame({"Temperature_K": [370], "Pressure_Pa": [1.1 * 101325]})
    y_test = keras_surrogate.evaluate_surrogate(x_test)

    # Test relu complementarity
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=keras_surrogate,
        formulation=KerasSurrogate.Formulation.RELU_COMPLEMENTARITY,
    )
    m.surrogate.inputs["Temperature_K"].fix(370)
    m.surrogate.inputs["Pressure_Pa"].fix(1.1 * 101325)
    solver = SolverFactory("ipopt")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {
            "EnthMol": [value(m.surrogate.outputs["EnthMol"])],
            "VapFrac": [value(m.surrogate.outputs["VapFrac"])],
        }
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)


@pytest.mark.unit
@pytest.mark.skipif(not SolverFactory("glpk").available(False), reason="no glpk")
def test_keras_surrogate_auto_creating_variables_glpk():
    ###
    # Test 1->2 relu
    ###
    keras_surrogate = create_keras_model(
        name="T_data_1_10_10_2_relu", return_keras_model_only=False
    )
    x_test = pd.DataFrame({"Temperature_K": [370]})
    y_test = keras_surrogate.evaluate_surrogate(x_test)

    # Test relu bigm
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=keras_surrogate,
        formulation=KerasSurrogate.Formulation.RELU_BIGM,
    )
    m.surrogate.inputs["Temperature_K"].fix(370)
    solver = SolverFactory("glpk")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {
            "EnthMol": [value(m.surrogate.outputs["EnthMol"])],
            "VapFrac": [value(m.surrogate.outputs["VapFrac"])],
        }
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)

    ###
    # Test 2->2 relu
    ###
    keras_surrogate = create_keras_model(
        name="PT_data_2_10_10_2_relu", return_keras_model_only=False
    )
    x_test = pd.DataFrame({"Temperature_K": [370], "Pressure_Pa": [1.1 * 101325]})
    y_test = keras_surrogate.evaluate_surrogate(x_test)

    # Test relu bigm
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=keras_surrogate,
        formulation=KerasSurrogate.Formulation.RELU_BIGM,
    )
    m.surrogate.inputs["Temperature_K"].fix(370)
    m.surrogate.inputs["Pressure_Pa"].fix(1.1 * 101325)
    solver = SolverFactory("glpk")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {
            "EnthMol": [value(m.surrogate.outputs["EnthMol"])],
            "VapFrac": [value(m.surrogate.outputs["VapFrac"])],
        }
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)


@pytest.mark.unit
@pytest.mark.skipif(not SolverFactory("ipopt").available(False), reason="no Ipopt")
def test_keras_surrogate_with_variables():
    keras_surrogate = create_keras_model(
        name="PT_data_2_10_10_2_sigmoid", return_keras_model_only=False
    )
    x_test = pd.DataFrame({"Temperature_K": [370], "Pressure_Pa": [1.1 * 101325]})
    y_test = keras_surrogate.evaluate_surrogate(x_test)

    # Test provide scalar inputs, auto create outputs
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.T = Var()
    m.P = Var()
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=keras_surrogate,
        input_vars=[m.T, m.P],
        formulation=KerasSurrogate.Formulation.FULL_SPACE,
    )
    m.T.fix(370)
    m.P.fix(1.1 * 101325)
    solver = SolverFactory("ipopt")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {
            "EnthMol": [value(m.surrogate.outputs["EnthMol"])],
            "VapFrac": [value(m.surrogate.outputs["VapFrac"])],
        }
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)

    # Test provide indexed inputs, auto create outputs
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.inputs = Var(["T", "P"])
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=keras_surrogate,
        input_vars=[m.inputs],
        formulation=KerasSurrogate.Formulation.FULL_SPACE,
    )
    m.inputs["T"].fix(370)
    m.inputs["P"].fix(1.1 * 101325)
    solver = SolverFactory("ipopt")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {
            "EnthMol": [value(m.surrogate.outputs["EnthMol"])],
            "VapFrac": [value(m.surrogate.outputs["VapFrac"])],
        }
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)

    # Test auto-create inputs, provide scalar outputs
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.H = Var()
    m.V = Var()
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=keras_surrogate,
        output_vars=[m.H, m.V],
        formulation=KerasSurrogate.Formulation.FULL_SPACE,
    )
    m.surrogate.inputs["Temperature_K"].fix(370)
    m.surrogate.inputs["Pressure_Pa"].fix(1.1 * 101325)
    solver = SolverFactory("ipopt")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame({"EnthMol": [value(m.H)], "VapFrac": [value(m.V)]})
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)

    # Test auto-create inputs, provide indexed outputs
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.outputs = Var(["H", "V"])
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=keras_surrogate,
        output_vars=[m.outputs],
        formulation=KerasSurrogate.Formulation.FULL_SPACE,
    )
    m.surrogate.inputs["Temperature_K"].fix(370)
    m.surrogate.inputs["Pressure_Pa"].fix(1.1 * 101325)
    solver = SolverFactory("ipopt")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {"EnthMol": [value(m.outputs["H"])], "VapFrac": [value(m.outputs["V"])]}
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)

    # Test provide scalar inputs, provide indexed outputs
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.T = Var()
    m.P = Var()
    m.outputs = Var(["H", "V"])
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=keras_surrogate,
        input_vars=[m.T, m.P],
        output_vars=[m.outputs],
        formulation=KerasSurrogate.Formulation.FULL_SPACE,
    )
    m.T.fix(370)
    m.P.fix(1.1 * 101325)
    solver = SolverFactory("ipopt")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {"EnthMol": [value(m.outputs["H"])], "VapFrac": [value(m.outputs["V"])]}
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)


@pytest.mark.unit
@pytest.mark.skipif(not SolverFactory("ipopt").available(False), reason="no Ipopt")
def test_save_load():
    keras_surrogate = create_keras_model(
        name="PT_data_2_10_10_2_sigmoid", return_keras_model_only=False
    )

    new_keras_surrogate = None
    dname = None
    with TempfileManager.new_context() as tf:
        dname = tf.mkdtemp()
        keras_surrogate.save_to_folder(dname)
        assert os.path.isdir(dname)
        assert os.path.isfile(os.path.join(dname, "idaes_info.json"))

        new_keras_surrogate = KerasSurrogate.load_from_folder(dname)

    # Check for clean up
    assert not os.path.isdir(dname)

    # check surrogate data members
    assert new_keras_surrogate._input_labels == ["Temperature_K", "Pressure_Pa"]
    assert new_keras_surrogate._output_labels == ["EnthMol", "VapFrac"]
    assert sorted(new_keras_surrogate._input_bounds.keys()) == [
        "Pressure_Pa",
        "Temperature_K",
    ]
    assert new_keras_surrogate._input_bounds["Temperature_K"][0] == pytest.approx(360.0)
    assert new_keras_surrogate._input_bounds["Temperature_K"][1] == pytest.approx(380.0)
    assert new_keras_surrogate._input_bounds["Pressure_Pa"][0] == pytest.approx(
        101325.0
    )
    assert new_keras_surrogate._input_bounds["Pressure_Pa"][1] == pytest.approx(
        1.2 * 101325.0
    )

    # check input scaler
    expected_columns = ["Temperature_K", "Pressure_Pa"]
    offset_series = pd.Series(
        {"Temperature_K": 369.983611, "Pressure_Pa": 111421.319811}
    )
    factor_series = pd.Series({"Temperature_K": 5.836047, "Pressure_Pa": 5917.954504})
    scaler = new_keras_surrogate._input_scaler
    assert scaler._expected_columns == expected_columns
    pd.testing.assert_series_equal(scaler._offset, offset_series, rtol=rtol, atol=atol)
    pd.testing.assert_series_equal(scaler._factor, factor_series, rtol=rtol, atol=atol)

    # check output scaler
    expected_columns = ["EnthMol", "VapFrac"]
    offset_series = pd.Series({"EnthMol": 54599.629980, "VapFrac": 0.403307})
    factor_series = pd.Series({"EnthMol": 14654.226615, "VapFrac": 0.430181})
    scaler = new_keras_surrogate._output_scaler
    assert scaler._expected_columns == expected_columns
    pd.testing.assert_series_equal(scaler._offset, offset_series, rtol=rtol, atol=atol)
    pd.testing.assert_series_equal(scaler._factor, factor_series, rtol=rtol, atol=atol)

    # check evaluation
    x_test = pd.DataFrame(
        {
            "Temperature_K": [360, 370, 380],
            "Pressure_Pa": [1.05 * 101325, 1.10 * 101325, 1.15 * 101325],
        }
    )
    y_test = keras_surrogate.evaluate_surrogate(x_test)
    expected_y = pd.DataFrame(
        {
            "EnthMol": [40194.5586954288, 48660.288218426984, 75178.30324367314],
            "VapFrac": [0.002291496299564877, 0.21942246438431742, 0.9996716243380308],
        }
    )
    pd.testing.assert_frame_equal(y_test, expected_y, rtol=rtol, atol=atol)

    # check solve with pyomo
    x_test = pd.DataFrame({"Temperature_K": [370], "Pressure_Pa": [1.1 * 101325]})
    y_test = keras_surrogate.evaluate_surrogate(x_test)

    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=keras_surrogate,
        formulation=KerasSurrogate.Formulation.FULL_SPACE,
    )
    m.surrogate.inputs["Temperature_K"].fix(370)
    m.surrogate.inputs["Pressure_Pa"].fix(1.1 * 101325)
    solver = SolverFactory("ipopt")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {
            "EnthMol": [value(m.surrogate.outputs["EnthMol"])],
            "VapFrac": [value(m.surrogate.outputs["VapFrac"])],
        }
    )
    pd.testing.assert_frame_equal(y_test, y_test_pyomo, rtol=rtol, atol=atol)


@pytest.mark.unit
@pytest.mark.skipif(not SolverFactory("ipopt").available(False), reason="no Ipopt")
def test_noscalers():
    keras_folder_name = os.path.join(this_file_dir(), "data", "keras_models")
    keras_model = load_keras_json_hd5(keras_folder_name, "PT_data_2_10_10_2_sigmoid")

    input_labels = ["Temperature_K", "Pressure_Pa"]
    output_labels = ["EnthMol", "VapFrac"]
    input_bounds = {"Temperature_K": (-3.0, 3.0), "Pressure_Pa": (-3.0, 3.0)}

    keras_surrogate = KerasSurrogate(
        keras_model=keras_model,
        input_labels=input_labels,
        output_labels=output_labels,
        input_bounds=input_bounds,
    )
    # check solve with pyomo
    x_test = pd.DataFrame({"Temperature_K": [0.5], "Pressure_Pa": [0.5]})
    y_test = keras_surrogate.evaluate_surrogate(x_test)

    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=keras_surrogate,
        formulation=KerasSurrogate.Formulation.FULL_SPACE,
    )
    m.surrogate.inputs["Temperature_K"].fix(0.5)
    m.surrogate.inputs["Pressure_Pa"].fix(0.5)
    solver = SolverFactory("ipopt")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)

    y_test_pyomo = pd.DataFrame(
        {
            "EnthMol": [value(m.surrogate.outputs["EnthMol"])],
            "VapFrac": [value(m.surrogate.outputs["VapFrac"])],
        }
    )
    pd.testing.assert_frame_equal(
        y_test, y_test_pyomo, check_dtype=False, rtol=rtol, atol=atol
    )


@pytest.mark.unit
def test_invalid_formulation():
    keras_surrogate = create_keras_model(
        name="T_data_1_10_10_2_sigmoid", return_keras_model_only=False
    )
    m = ConcreteModel()
    m.obj = Objective(expr=1)
    m.surrogate = SurrogateBlock()
    with pytest.raises(ValueError) as excinfo:
        m.surrogate.build_model(surrogate_object=keras_surrogate, formulation="foo")
    assert (
        str(excinfo.value) == 'An unrecognized formulation "foo" was passed '
        "to KerasSurrogate.populate_block. Please pass a valid formulation."
    )
