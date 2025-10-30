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
Interface for importing Keras models into IDAES
"""
# TODO: Missing docstrings
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring

import json
import os.path

import pandas as pd

from pyomo.common.dependencies import attempt_import

from idaes.core.surrogate.sampling.scaling import OffsetScaler

from idaes.core.surrogate.omlt_base_surrogate_class import OMLTSurrogate

# pylint: disable=possibly-used-before-assignment
keras, keras_available = attempt_import("tensorflow.keras")
omlt, omlt_available = attempt_import("omlt")

if omlt_available:
    from omlt.neuralnet.nn_formulation import (
        FullSpaceSmoothNNFormulation,
        ReducedSpaceSmoothNNFormulation,
        ReluBigMFormulation,
        ReluComplementarityFormulation,
    )

    if keras_available:
        from omlt.io import load_keras_sequential


class KerasSurrogate(OMLTSurrogate):
    def __init__(
        self,
        keras_model,
        input_labels,
        output_labels,
        input_bounds,
        input_scaler=None,
        output_scaler=None,
    ):
        """
        Standard SurrogateObject for surrogates based on Keras models.
        Utilizes the OMLT framework for importing Keras models to IDAES.

        Contains methods to both populate a Pyomo Block with constraints
        representing the surrogate and to evaluate the surrogate a set of user
        provided points.

        This constructor should only be used when first creating the surrogate within IDAES.
        Once created, this object can be stored to disk using save_to_folder and loaded
        with load_from_folder

        Args:
           keras_model: Keras Sequential model
              This is the Keras Sequential model that will be loaded. Note that
              specialized layers may not be supported at this time.
           input_labels: list of str
              The ordered list of labels corresponding to the inputs in the keras model
           output_labels: list of str
              The ordered list of labels corresponding to the outputs in the keras model
           input_bounds: None of dict of tuples
              Keys correspond to each of the input labels and values are the tuples of
              bounds (lb, ub)
           input_scaler: None or OffsetScaler
              The scaler to be used for the inputs. If None, then no scaler is used
           output_scaler: None of OffsetScaler
              The scaler to be used for the outputs. If None, then no scaler is used
        """
        super().__init__(
            input_labels=input_labels,
            output_labels=output_labels,
            input_bounds=input_bounds,
            input_scaler=input_scaler,
            output_scaler=output_scaler,
        )
        self._keras_model = keras_model

    def populate_block(self, block, additional_options=None):
        """
        Method to populate a Pyomo Block with the keras model constraints.

        Args:
           block: Pyomo Block component
              The block to be populated with variables and/or constraints.
           additional_options: dict or None
              If not None, then should be a dict with the following keys;
              'formulation': KerasSurrogate.Formulation
              The formulation to use with OMLT. Possible values are FULL_SPACE,
              REDUCED_SPACE, RELU_BIGM, or RELU_COMPLEMENTARITY (default is FULL_SPACE)
        """
        formulation = additional_options.pop(
            "formulation", KerasSurrogate.Formulation.FULL_SPACE
        )
        omlt_scaling, scaled_input_bounds = self.generate_omlt_scaling_objecets()

        # pylint: disable-next=used-before-assignment
        net = load_keras_sequential(
            self._keras_model,
            scaling_object=omlt_scaling,
            scaled_input_bounds=scaled_input_bounds,
        )

        if formulation == KerasSurrogate.Formulation.FULL_SPACE:
            formulation_object = FullSpaceSmoothNNFormulation(net)
        elif formulation == KerasSurrogate.Formulation.REDUCED_SPACE:
            formulation_object = ReducedSpaceSmoothNNFormulation(net)
        elif formulation == KerasSurrogate.Formulation.RELU_BIGM:
            formulation_object = ReluBigMFormulation(net)
        elif formulation == KerasSurrogate.Formulation.RELU_COMPLEMENTARITY:
            formulation_object = ReluComplementarityFormulation(net)
        else:
            raise ValueError(
                'An unrecognized formulation "{}" was passed to '
                "KerasSurrogate.populate_block. Please pass a valid "
                "formulation.".format(formulation)
            )
        self.populate_block_with_net(block, formulation_object)

    def evaluate_surrogate(self, inputs):
        """
        Method to evaluate Keras model at a set of input values.

        Args:
            inputs: numpy array of input values. First dimension of array
                must match the number of input variables.

        Returns:
            outputs: numpy array of values for all outputs evaluated at input
                points.
        """
        x = inputs
        if self._input_scaler is not None:
            x = self._input_scaler.scale(x)
        y = self._keras_model.predict(x.to_numpy())

        # y is a numpy array, make it a dataframe
        y = pd.DataFrame(
            data=y, columns=self.output_labels(), index=inputs.index, dtype="float64"
        )
        if self._output_scaler is not None:
            y = self._output_scaler.unscale(y)
        return y

    def save_to_folder(self, keras_folder_name, keras_model_name="idaes_keras_model"):
        """
        Save the surrogate object to disk by providing the name of the
        folder to contain the keras model and additional IDAES metadata

        Args:
           folder_name: str
              The name of the folder to contain the Keras model and additional
              IDAES metadata
        """
        self._keras_model.save(
            os.path.join(keras_folder_name, keras_model_name + ".keras")
        )
        info = dict()
        info["input_scaler"] = None
        if self._input_scaler is not None:
            info["input_scaler"] = self._input_scaler.to_dict()
        info["output_scaler"] = None
        if self._output_scaler is not None:
            info["output_scaler"] = self._output_scaler.to_dict()

        # serialize information from the base class
        info["input_labels"] = self.input_labels()
        info["output_labels"] = self.output_labels()
        info["input_bounds"] = self.input_bounds()

        with open(os.path.join(keras_folder_name, "idaes_info.json"), "w") as fd:
            json.dump(info, fd)

    @classmethod
    def load_from_folder(cls, keras_folder_name, keras_model_name="idaes_keras_model"):
        """
        Load the surrogate object from disk by providing the name of the
        folder holding the keras model

        Args:
           folder_name: str
              The name of the folder containing the Keras model and additional
              IDAES metadata

        Returns: an instance of KerasSurrogate
        """

        keras_model = keras.models.load_model(
            os.path.join(keras_folder_name, keras_model_name + ".keras")
        )

        with open(os.path.join(keras_folder_name, "idaes_info.json")) as fd:
            info = json.load(fd)

        input_scaler = None
        if info["input_scaler"] is not None:
            input_scaler = OffsetScaler.from_dict(info["input_scaler"])

        output_scaler = None
        if info["output_scaler"] is not None:
            output_scaler = OffsetScaler.from_dict(info["output_scaler"])

        return KerasSurrogate(
            keras_model=keras_model,
            input_labels=info["input_labels"],
            output_labels=info["output_labels"],
            input_bounds=info["input_bounds"],
            input_scaler=input_scaler,
            output_scaler=output_scaler,
        )


def save_keras_json_hd5(nn, path, name):
    json_model = nn.to_json()
    with open(os.path.join(path, "{}.json".format(name)), "w") as json_file:
        json_file.write(json_model)
    nn.save(os.path.join(path, "{}.keras".format(name)))
    nn.save_weights(os.path.join(path, "{}.weights.h5".format(name)))


def load_keras_json_hd5(path, name):
    nn = keras.models.load_model(os.path.join(path, "{}.keras".format(name)))
    nn.load_weights(os.path.join(path, "{}.weights.h5".format(name)))
    return nn
