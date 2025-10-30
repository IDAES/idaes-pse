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
Interface for importing ONNX models into IDAES
"""
# TODO: Missing docstrings
# pylint: disable=missing-class-docstring
# pylint: disable=missing-function-docstring
# TODO: Importing protected  _ACTIVATION_OP_TYPES as not exposed in distributed version
# pylint: disable=W0123
from enum import Enum
import json
import os.path

from pyomo.common.dependencies import attempt_import

from idaes.core.surrogate.sampling.scaling import OffsetScaler

from idaes.core.surrogate.omlt_base_surrogate_class import OMLTSurrogate

# pylint: disable=possibly-used-before-assignment
onnx, onnx_available = attempt_import("onnx")
omlt, omlt_available = attempt_import("omlt")

if omlt_available:
    from omlt.neuralnet import (
        FullSpaceSmoothNNFormulation,
        ReducedSpaceSmoothNNFormulation,
        ReluBigMFormulation,
        ReluComplementarityFormulation,
    )
    import omlt.io as omltio

    if onnx_available:
        from omlt.io import load_onnx_neural_network, write_onnx_model_with_bounds


class ONNXSurrogate(OMLTSurrogate):
    def __init__(
        self,
        onnx_model,
        input_labels,
        output_labels,
        input_bounds,
        input_scaler=None,
        output_scaler=None,
    ):
        """
        Standard SurrogateObject for surrogates based on ONNX models.
        Utilizes the OMLT framework for importing ONNX models to IDAES.

        Contains methods to both populate a Pyomo Block with constraints
        representing the surrogate and to evaluate the surrogate a set of user
        provided points.

        This constructor should only be used when first creating the surrogate within IDAES.
        Once created, this object can be stored to disk using save_to_folder and loaded
        with load_from_folder

        Args:
           onnx_model: Onnx model file to be loaded.
           input_labels: list of str
              The ordered list of labels corresponding to the inputs in the onnx model
           output_labels: list of str
              The ordered list of labels corresponding to the outputs in the onnx model
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

        self._onnx_model = onnx_model

    class Formulation(Enum):
        FULL_SPACE = 1
        REDUCED_SPACE = 2
        RELU_BIGM = 3
        RELU_COMPLEMENTARITY = 4

    def populate_block(self, block, additional_options=None):
        """
        Method to populate a Pyomo Block with the onnx model constraints.

        Args:
           block: Pyomo Block component
              The block to be populated with variables and/or constraints.
           additional_options: dict or None
              If not None, then should be a dict with the following keys;
              'formulation': ONNXSurrogate.Formulation
              The formulation to use with OMLT. Possible values are FULL_SPACE,
              REDUCED_SPACE, RELU_BIGM, or RELU_COMPLEMENTARITY (default is FULL_SPACE)
        """
        formulation = additional_options.pop(
            "formulation", ONNXSurrogate.Formulation.REDUCED_SPACE
        )
        omlt_scaling, scaled_input_bounds = self.generate_omlt_scaling_objecets()

        # omlt takes *scaled* input bounds as a dictionary with int keys
        input_bounds = dict(enumerate(self.input_bounds().values()))
        scaled_input_bounds = omlt_scaling.get_scaled_input_expressions(input_bounds)
        scaled_input_bounds = {i: tuple(bnd) for i, bnd in scaled_input_bounds.items()}

        # TODO: remove this once new OMLT 1.2 is made available and includes tanh support
        # overrides default available activation functions for ONNX, tanh is not listed in 1.1 but is supported

        omltio.onnx_parser._ACTIVATION_OP_TYPES = [  # pylint: disable=protected-access
            "Relu",
            "Sigmoid",
            "LogSoftmax",
            "Tanh",
        ]

        # pylint: disable-next=used-before-assignment
        net = load_onnx_neural_network(
            self._onnx_model,
            scaling_object=omlt_scaling,
            input_bounds=scaled_input_bounds,
        )

        if formulation == ONNXSurrogate.Formulation.FULL_SPACE:
            formulation_object = FullSpaceSmoothNNFormulation(net)
        elif formulation == ONNXSurrogate.Formulation.REDUCED_SPACE:
            formulation_object = ReducedSpaceSmoothNNFormulation(net)
        elif formulation == ONNXSurrogate.Formulation.RELU_BIGM:
            formulation_object = ReluBigMFormulation(net)
        elif formulation == ONNXSurrogate.Formulation.RELU_COMPLEMENTARITY:
            formulation_object = ReluComplementarityFormulation(net)
        else:
            raise ValueError(
                'An unrecognized formulation "{}" was passed to '
                "ONNXSurrogate.populate_block. Please pass a valid "
                "formulation.".format(formulation)
            )
        self.populate_block_with_net(block, formulation_object)

    def evaluate_surrogate(self, inputs):
        """
        Method to evaluate ONNX model at a set of input values.

        Args:
            inputs: numpy array of input values. First dimension of array
                must match the number of input variables.

        Returns:
            outputs: numpy array of values for all outputs evaluated at input
                points.
        """
        raise NotImplementedError

    def save_to_folder(self, save_location, save_name):
        """
        Save the surrogate object to disk by providing the location to store the
                model in and its name, as well as additional IDAES metadata

        Args:
           save_location: str
              The name of the folder to contain the ONNX model and additional
              IDAES metadata
            save_name: str
                The name for the model
        """

        # pylint: disable-next=used-before-assignment
        write_onnx_model_with_bounds(
            os.path.join(save_location, "{}.onnx".format(save_name)),
            onnx_model=self._onnx_model,
            input_bounds=None,
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

        with open(
            os.path.join(save_location, "{}_idaes_info.json".format(save_name)), "w"
        ) as fd:
            json.dump(info, fd)

    @classmethod
    def load_onnx_model(cls, onnx_model_location, model_name):
        """
        Load the surrogate object from disk by providing the name of the
        folder holding the onnx model and its name, including accompanying json file that includes following
        structure:

            "input_scaler":{
                "expected_columns":[list of input_keys],
                "offset":{"input_key":offset_value,etc.},
                "factor":{"input_key":factor_value (e.g. multiplier),etc.}}

            "output_scaler":{
                "expected_columns":[list of output_keys],
                "offset":{"output_key":offset_value,etc.},
                "factor":{"output_key":factor_value (e.g. multiplier),etc.}}

            "input_bounds":{"input_key":[low_bound,high_bound],etc.}
            "input_labels":[list of input_keys]
            "output_labels":[list of output_keys]

        Args:
           folder_name: str
              The name of the folder containing the onnx model and additional
              IDAES metadata
            model_name: str
              The name of the model to load in the folder

        Returns: an instance of ONNXSurrogate
        """
        onnx_model = onnx.load(
            os.path.join(onnx_model_location, "{}.onnx".format(model_name))
        )
        with open(
            os.path.join(onnx_model_location, "{}_idaes_info.json".format(model_name))
        ) as fd:
            info = json.load(fd)

        input_scaler = None
        if info["input_scaler"] is not None:
            input_scaler = OffsetScaler.from_dict(info["input_scaler"])

        output_scaler = None
        if info["output_scaler"] is not None:
            output_scaler = OffsetScaler.from_dict(info["output_scaler"])

        return ONNXSurrogate(
            onnx_model=onnx_model,
            input_labels=info["input_labels"],
            output_labels=info["output_labels"],
            input_bounds=info["input_bounds"],
            input_scaler=input_scaler,
            output_scaler=output_scaler,
        )
