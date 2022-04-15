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
Code to create the keras models used in testing
"""
import tensorflow.keras as keras
from keras.callbacks import ModelCheckpoint
from tensorflow.keras.optimizers import Adam
from tensorflow.keras.layers import Dense
from tensorflow.keras.models import Sequential
import pandas as pd
import matplotlib.pyplot as plt

import idaes.core.surrogate.sampling as sampling
from idaes.core.surrogate.sampling.scaling import OffsetScaler
import idaes.core.surrogate.keras_surrogate as foo

print(dir(foo))
from idaes.core.surrogate.keras_surrogate import save_keras_json_hd5
from pyomo.common.fileutils import this_file_dir
import os
import json

keras.utils.set_random_seed(42)

generate_fit_comparison_files = True


def compare_fit(nn, x_test, y_test, input_scaler, output_scaler):
    scaled_x_test = input_scaler.scale(x_test)
    scaled_y_test_hat = nn.predict(scaled_x_test.to_numpy())
    scaled_y_test_hat = pd.DataFrame(
        data=scaled_y_test_hat, columns=["EnthMol", "VapFrac"]
    )
    y_test_hat = output_scaler.unscale(scaled_y_test_hat)
    y_test_hat.rename(
        columns={"EnthMol": "EnthMol_hat", "VapFrac": " VapFrac_hat"}, inplace=True
    )
    df = pd.concat([x_test, y_test, y_test_hat], axis=1)
    df.plot(
        x="Temperature_K", y=["EnthMol", "EnthMol_hat"], marker=".", linestyle="none"
    )
    plt.savefig("fit_{}.png".format(nn.name))
    df.to_csv("fit_{}.csv".format(nn.name), index=False)


df = pd.read_csv(os.path.join(this_file_dir(), "T_data.csv"))
train_val, testing = sampling.split_dataframe(df, [0.8])
input_labels = ["Temperature_K"]
output_labels = ["EnthMol", "VapFrac"]
inputs = train_val[input_labels]
outputs = train_val[output_labels]

input_scaler = OffsetScaler.create_from_mean_std(inputs)
output_scaler = OffsetScaler.create_from_mean_std(outputs)
print("input_scaler for T_data")
print(input_scaler._offset)
print(input_scaler._factor)
print("output_scaler for T_data")
print(output_scaler._offset)
print(output_scaler._factor)
scaled_inputs = input_scaler.scale(inputs)
scaled_outputs = output_scaler.scale(outputs)

x = scaled_inputs.to_numpy()
y = scaled_outputs.to_numpy()

nn = Sequential(name="T_data_1_10_10_2_sigmoid")
nn.add(Dense(units=10, input_dim=1, activation="sigmoid"))
nn.add(Dense(units=10, activation="sigmoid"))
nn.add(Dense(units=2))
nn.compile(optimizer=Adam(), loss="mse")

mcp_save = ModelCheckpoint(
    ".mdl_wts.hdf5", save_best_only=True, monitor="val_loss", mode="min"
)
history = nn.fit(
    x=x,
    y=y,
    validation_split=0.2,
    batch_size=16,
    verbose=1,
    epochs=1000,
    callbacks=[mcp_save],
)
# nn.save(os.path.join(this_file_dir(), 'keras_models', nn.name))
save_keras_json_hd5(nn, os.path.join(this_file_dir(), "keras_models"), nn.name)

if generate_fit_comparison_files:
    # x_test = pd.DataFrame({'Temperature_K': [365, 370, 375]})
    # compare_fit(nn, x_test, testing[output_labels], input_scaler, output_scaler)
    compare_fit(
        nn, testing[input_labels], testing[output_labels], input_scaler, output_scaler
    )

# RELU Network
nn = Sequential(name="T_data_1_10_10_2_relu")
nn.add(Dense(units=10, input_dim=1, activation="relu"))
nn.add(Dense(units=10, activation="relu"))
nn.add(Dense(units=2))
nn.compile(optimizer=Adam(), loss="mse")

mcp_save = ModelCheckpoint(
    ".mdl_wts.hdf5", save_best_only=True, monitor="val_loss", mode="min"
)
history = nn.fit(
    x=x,
    y=y,
    validation_split=0.2,
    batch_size=16,
    verbose=1,
    epochs=500,
    callbacks=[mcp_save],
)
# nn.save(os.path.join(this_file_dir(), 'keras_models', nn.name))
save_keras_json_hd5(nn, os.path.join(this_file_dir(), "keras_models"), nn.name)

if generate_fit_comparison_files:
    # x_test = pd.DataFrame({'Temperature_K': [365, 370, 375]})
    # compare_fit(nn, xtest, testing[output_labels], input_scaler, output_scaler)
    compare_fit(
        nn, testing[input_labels], testing[output_labels], input_scaler, output_scaler
    )

# PT DATA
df = pd.read_csv(os.path.join(this_file_dir(), "PT_data.csv"))
train_val, testing = sampling.split_dataframe(df, [0.8])
input_labels = ["Temperature_K", "Pressure_Pa"]
output_labels = ["EnthMol", "VapFrac"]
inputs = train_val[input_labels]
outputs = train_val[output_labels]

input_scaler = OffsetScaler.create_from_mean_std(inputs)
output_scaler = OffsetScaler.create_from_mean_std(outputs)
scaled_inputs = input_scaler.scale(inputs)
scaled_outputs = output_scaler.scale(outputs)
print("input_scaler for PT_data")
print(input_scaler._offset)
print(input_scaler._factor)
print("output_scaler for PT_data")
print(output_scaler._offset)
print(output_scaler._factor)

x = scaled_inputs.to_numpy()
y = scaled_outputs.to_numpy()

nn = Sequential(name="PT_data_2_10_10_2_sigmoid")
nn.add(Dense(units=10, input_dim=2, activation="sigmoid"))
nn.add(Dense(units=10, activation="sigmoid"))
nn.add(Dense(units=2))
nn.compile(optimizer=Adam(), loss="mse")

mcp_save = ModelCheckpoint(
    ".mdl_wts.hdf5", save_best_only=True, monitor="val_loss", mode="min"
)
history = nn.fit(
    x=x,
    y=y,
    validation_split=0.2,
    batch_size=16,
    verbose=1,
    epochs=200,
    callbacks=[mcp_save],
)
# nn.save(os.path.join(this_file_dir(), 'keras_models', nn.name))
save_keras_json_hd5(nn, os.path.join(this_file_dir(), "keras_models"), nn.name)

if generate_fit_comparison_files:
    compare_fit(
        nn, testing[input_labels], testing[output_labels], input_scaler, output_scaler
    )

# RELU Network
nn = Sequential(name="PT_data_2_10_10_2_relu")
nn.add(Dense(units=10, input_dim=2, activation="relu"))
nn.add(Dense(units=10, activation="relu"))
nn.add(Dense(units=2))
nn.compile(optimizer=Adam(), loss="mse")

mcp_save = ModelCheckpoint(
    ".mdl_wts.hdf5", save_best_only=True, monitor="val_loss", mode="min"
)
history = nn.fit(
    x=x,
    y=y,
    validation_split=0.2,
    batch_size=16,
    verbose=1,
    epochs=50,
    callbacks=[mcp_save],
)
# nn.save(os.path.join(this_file_dir(), 'keras_models', nn.name))
save_keras_json_hd5(nn, os.path.join(this_file_dir(), "keras_models"), nn.name)

if generate_fit_comparison_files:
    compare_fit(
        nn, testing[input_labels], testing[output_labels], input_scaler, output_scaler
    )
