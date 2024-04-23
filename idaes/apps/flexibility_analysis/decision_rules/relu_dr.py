import numpy as np
import tensorflow as tf
from tensorflow import keras
from tensorflow.keras import layers
import pyomo.environ as pe
from pyomo.core.base.var import _GeneralVarData
from pyomo.core.base.block import _BlockData
from typing import MutableMapping, Sequence
from omlt import OmltBlock, OffsetScaling
from omlt.neuralnet import ReluBigMFormulation
from omlt.io import load_keras_sequential
from idaes.apps.flexibility_analysis.indices import _VarIndex
from .relu_dr_config import ReluDRConfig


def construct_relu_decision_rule(
    input_vals: MutableMapping[_GeneralVarData, Sequence[float]],
    output_vals: MutableMapping[_GeneralVarData, Sequence[float]],
    config: ReluDRConfig,
) -> _BlockData:
    tf.random.set_seed(config.tensorflow_seed)
    inputs = list(input_vals.keys())
    outputs = list(output_vals.keys())
    n_samples = len(input_vals[inputs[0]])

    config: ReluDRConfig = config()
    if config.batch_size > n_samples:
        config.batch_size = n_samples

    training_input = np.empty((n_samples, len(inputs)))
    for ndx, inp in enumerate(inputs):
        training_input[:, ndx] = np.array(input_vals[inp], dtype=np.float64)

    training_output = np.empty((n_samples, len(outputs)))
    for ndx, outp in enumerate(outputs):
        training_output[:, ndx] = np.array(output_vals[outp], dtype=np.float64)

    if config.scale_inputs:
        input_mean = training_input.mean(axis=0)
        input_std = training_input.std(axis=0)
        for ndx in range(len(inputs)):
            training_input[:, ndx] = (
                training_input[:, ndx] - input_mean[ndx]
            ) / input_std[ndx]
    else:
        input_mean = [0] * len(inputs)
        input_std = [1] * len(inputs)

    if config.scale_outputs:
        output_mean = training_output.mean(axis=0)
        output_std = training_output.std(axis=0)
        for ndx in range(len(outputs)):
            training_output[:, ndx] = (
                training_output[:, ndx] - output_mean[ndx]
            ) / output_std[ndx]
    else:
        output_mean = [0] * len(outputs)
        output_std = [1] * len(outputs)

    nn = keras.Sequential()
    nn.add(
        layers.Dense(
            units=config.n_nodes_per_layer, input_dim=len(inputs), activation="relu"
        )
    )
    for layer_ndx in range(config.n_layers - 1):
        nn.add(layers.Dense(config.n_nodes_per_layer, activation="relu"))
    nn.add(layers.Dense(len(outputs)))
    if config.learning_rate is None:
        opt = keras.optimizers.Adam()
    else:
        opt = keras.optimizers.Adam(learning_rate=config.learning_rate)
    nn.compile(optimizer=opt, loss="mse")
    history = nn.fit(
        training_input,
        training_output,
        batch_size=config.batch_size,
        epochs=config.epochs,
        # verbose=0,
    )

    if config.plot_history:
        import matplotlib.pyplot as plt

        plt.scatter(history.epoch, history.history["loss"])
        plt.xlabel("Epoch")
        plt.ylabel("Loss")
        plt.yscale("log")
        plt.show()
        plt.close()

    res = pe.Block(concrete=True)
    res.nn = OmltBlock()

    scaler = OffsetScaling(
        offset_inputs=[float(i) for i in input_mean],
        factor_inputs=[float(i) for i in input_std],
        offset_outputs=[float(i) for i in output_mean],
        factor_outputs=[float(i) for i in output_std],
    )
    input_bounds = {
        ndx: (
            (v.lb - input_mean[ndx]) / input_std[ndx],
            (v.ub - input_mean[ndx]) / input_std[ndx],
        )
        for ndx, v in enumerate(inputs)
    }
    net = load_keras_sequential(nn, scaler, input_bounds)
    formulation = ReluBigMFormulation(net)
    res.nn.build_formulation(formulation)

    res.input_set = pe.Set()
    res.input_links = pe.Constraint(res.input_set)
    for ndx, v in enumerate(inputs):
        key = _VarIndex(v, None)
        res.input_set.add(key)
        res.input_links[key] = v == res.nn.inputs[ndx]

    res.output_set = pe.Set()
    res.output_links = pe.Constraint(res.output_set)
    for ndx, v in enumerate(outputs):
        key = _VarIndex(v, None)
        res.output_set.add(key)
        res.output_links[key] = v == res.nn.outputs[ndx]

    return res
