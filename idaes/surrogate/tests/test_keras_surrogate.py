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
import os.path
import pandas as pd
from pyomo.common.fileutils import this_file_dir
from idaes.surrogate.keras_surrogate import KerasSurrogate
from idaes.surrogate.surrogate_block import SurrogateBlock
from idaes.surrogate.sampling.scaling import StandardScaler
import tensorflow.keras as keras


@pytest.mark.unit
def test_KerasSurrogate_create_new():
    keras_folder_name = os.path.join(this_file_dir(), 'keras_models', 'T_data_1_10_10_2_sigmoid')
    inputs_scaler = StandardScaler(
        expected_columns=['Temperature_K'],
        offset_series = pd.Series({'Temperature_K': 369.977879}),
        factor_series = pd.Series({'Temperature_K':5.833396}))
    outputs_scaler = StandardScaler(
        expected_columns=['EnthMol','VapFrac'],
        offset_series = pd.Series({'EnthMol': 54582.841444, 'VapFrac': 0.402814}),
        factor_series = pd.Series({'EnthMol': 14642.206469, 'VapFrac': 0.429829}))

    keras_model = keras.models.load_model(keras_folder_name)

    # test exceptions for mismatched labels
    with pytest.raises(ValueError) as excinfo:
        keras_surrogate = KerasSurrogate(keras_model=keras_model, input_labels=['Temperature'],
                                         output_labels=['EnthMol', 'VapFrac'],
                                         input_bounds=[(360, 380)],
                                         input_scaler=inputs_scaler,
                                         output_scaler=outputs_scaler)
    assert str(excinfo.value) == "KerasSurrogate created with input_labels that do not match the expected columns in the input_scaler.\n" \
        "input_labels=['Temperature']\n" \
        "input_scaler.expected_columns()=['Temperature_K']"

    with pytest.raises(ValueError) as excinfo:
        keras_surrogate = KerasSurrogate(keras_model=keras_model, input_labels=['Temperature_K'],
                                         output_labels=['EnthMol', 'VapFraction'],
                                         input_bounds=[(360, 380)],
                                         input_scaler=inputs_scaler,
                                         output_scaler=outputs_scaler)
    assert str(excinfo.value) == "KerasSurrogate created with output_labels that do not match the expected columns in the output_scaler.\n" \
        "output_labels=['EnthMol', 'VapFraction']\n" \
        "output_scaler.expected_columns()=['EnthMol', 'VapFrac']"


    keras_surrogate = KerasSurrogate(keras_model=keras_model, input_labels=['Temperature_K'],
                                     output_labels=['EnthMol', 'VapFrac'],
                                     input_bounds=[(360, 380)],
                                     input_scaler=inputs_scaler,
                                     output_scaler=outputs_scaler
                                     )
    # create some data and test

    """
    # test two inputs / two outputs
    inputs_scaler = StandardScaler(
        expected_columns=['Temperature_K', 'Pressure_Pa'],
        offset_series = pd.Series({'Temperature_K': 369.977879, 'Pressure_Pa': 111428.04922}),
        factor_series = pd.Series({'Temperature_K':5.833396, 'Pressure_Pa': 5918.746802}))
    outputs_scaler = StandardScaler(
        expected_columns=['EnthMol','VapFrac'],
        offset_series = pd.Series({'EnthMol': 54582.841444, 'VapFrac': 0.402814}),
        factor_series = pd.Series({'EnthMol': 14642.206469, 'VapFrac': 0.429829}))
    """    
    
