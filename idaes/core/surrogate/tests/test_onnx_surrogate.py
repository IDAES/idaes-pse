#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Tests for ONNXSurrogate
"""
import pytest

pytest.importorskip("onnx", reason="onnx not available")
pytest.importorskip("omlt", reason="omlt not available")

import os.path
from pyomo.common.fileutils import this_file_dir
from pyomo.common.tempfiles import TempfileManager
from pyomo.environ import (
    ConcreteModel,
    Var,
    SolverFactory,
    assert_optimal_termination,
    value,
)
from idaes.core.surrogate.onnx_surrogate import ONNXSurrogate
from idaes.core.surrogate.surrogate_block import SurrogateBlock
from idaes.core.surrogate.sampling.scaling import OffsetScaler
import json
import onnx

# onnx, onnx_available = attempt_import("onnx")
rtol = 1e-4
atol = 1e-4


def load_onnx_model_data(
    name="net_st_net_5000_STM_100_s_2000000_60_5_tanh_1e-06_4096_tr_15481_Calcite_ST",
):
    onnx_folder_name = os.path.join(this_file_dir(), "data", "onnx_models")
    onnx_model = onnx.load(os.path.join(onnx_folder_name, "{}.onnx".format(name)))
    with open(os.path.join(onnx_folder_name, "{}_idaes_info.json".format(name))) as fd:
        scaler_info = json.load(fd)

    test_inputs = {
        "vars": [
            "feed_pH",
            "pressure_bar_feed",
            "Na",
            "Cl",
            "Ca",
            "Mg",
            "HCO3",
            "SO4",
            "K",
            "Sr",
            "Ba",
            "HCl",
        ],
        "feed_pH": 9.5,
        "pressure_bar_feed": 1.01325,
        "Na": 0.230858556000748,
        "Cl": 0.106701648328453,
        "Ca": 0.0245274696499109,
        "Mg": 0.0311348703689873,
        "HCO3": 0.430482141673564,
        "SO4": 0.182204065,
        "K": 0.000500561,
        "Sr": 0.000761853,
        "Ba": 2.50e-05,
        "HCl": 10,
    }
    test_outputs = {"vars": ["Calcite_ST"], "Calcite_ST": 40.4270772546529}
    return onnx_model, scaler_info, test_inputs, test_outputs


@pytest.mark.unit
@pytest.mark.skipif(not SolverFactory("ipopt").available(False), reason="no Ipopt")
def test_onnx_surrogate_manual_creation():
    ###
    # Test 1->2 sigmoid
    ###
    onnx_model, scaler_info, test_inputs, test_outputs = load_onnx_model_data()
    input_scaler = None
    for key, items in scaler_info.items():
        print(key, items)
    if scaler_info["input_scaler"] is not None:
        input_scaler = OffsetScaler.from_dict(scaler_info["input_scaler"])

    output_scaler = None
    if scaler_info["output_scaler"] is not None:
        output_scaler = OffsetScaler.from_dict(scaler_info["output_scaler"])
    onnx_surrogate = ONNXSurrogate(
        onnx_model,
        input_labels=scaler_info["input_labels"],
        output_labels=scaler_info["output_labels"],
        input_bounds=scaler_info["input_bounds"],
        input_scaler=input_scaler,
        output_scaler=output_scaler,
    )

    m = ConcreteModel()

    output_vars = ["Calcite_ST"]
    m.inputs = Var(test_inputs["vars"])

    m.outputs = Var(test_outputs["vars"])
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=onnx_surrogate,
        input_vars=[m.inputs[input_var] for input_var in test_inputs["vars"]],
        output_vars=[m.outputs[output_var] for output_var in test_outputs["vars"]],
        formulation=ONNXSurrogate.Formulation.REDUCED_SPACE,
    )
    for key in test_inputs["vars"]:
        m.inputs[key].fix(test_inputs[key])

    solver = SolverFactory("ipopt")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)
    assert pytest.approx(test_outputs["Calcite_ST"], rel=1e-3) == value(
        m.outputs["Calcite_ST"]
    )


@pytest.mark.unit
@pytest.mark.skipif(not SolverFactory("ipopt").available(False), reason="no Ipopt")
def test_onnx_surrogate_load_and_save_from_file():
    ###
    # Test 1->2 sigmoid
    ###
    _, _, test_inputs, test_outputs = load_onnx_model_data()

    onnx_surrogate = ONNXSurrogate.load_onnx_model(
        onnx_model_location=os.path.join(this_file_dir(), "data", "onnx_models"),
        model_name="net_st_net_5000_STM_100_s_2000000_60_5_tanh_1e-06_4096_tr_15481_Calcite_ST",
    )
    with TempfileManager.new_context() as tf:
        dname = tf.mkdtemp()
        onnx_surrogate.save_to_folder(dname, "temp_model")

        loaded_onnx_surrogate = ONNXSurrogate.load_onnx_model(
            onnx_model_location=dname,
            model_name="temp_model",
        )
    assert not os.path.isdir(dname)
    m = ConcreteModel()

    output_vars = ["Calcite_ST"]
    m.inputs = Var(test_inputs["vars"])

    m.outputs = Var(test_outputs["vars"])
    m.surrogate = SurrogateBlock()
    m.surrogate.build_model(
        surrogate_object=loaded_onnx_surrogate,
        input_vars=[m.inputs[input_var] for input_var in test_inputs["vars"]],
        output_vars=[m.outputs[output_var] for output_var in test_outputs["vars"]],
        formulation=ONNXSurrogate.Formulation.REDUCED_SPACE,
    )
    for key in test_inputs["vars"]:
        m.inputs[key].fix(test_inputs[key])

    solver = SolverFactory("ipopt")
    status = solver.solve(m, tee=True)
    assert_optimal_termination(status)
    assert pytest.approx(test_outputs["Calcite_ST"], rel=1e-3) == value(
        m.outputs["Calcite_ST"]
    )
