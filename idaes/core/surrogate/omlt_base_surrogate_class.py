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

from enum import Enum
import numpy as np

from pyomo.common.dependencies import attempt_import

from idaes.core.surrogate.base.surrogate_base import SurrogateBase
from idaes.core.surrogate.sampling.scaling import OffsetScaler

# pylint: disable=possibly-used-before-assignment
keras, keras_available = attempt_import("tensorflow.keras")
omlt, omlt_available = attempt_import("omlt")

if omlt_available:
    from omlt import OmltBlock, OffsetScaling


class OMLTSurrogate(SurrogateBase):
    def __init__(
        self,
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
           onnx_model: Onnx model file to be loaded.
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
        )

        # make sure we are using the standard scaler
        if (
            input_scaler is not None
            and not isinstance(input_scaler, OffsetScaler)
            or output_scaler is not None
            and not isinstance(output_scaler, OffsetScaler)
        ):
            raise NotImplementedError("KerasSurrogate only supports the OffsetScaler.")

        # check that the input labels match
        if input_scaler is not None and input_scaler.expected_columns() != input_labels:
            raise ValueError(
                "KerasSurrogate created with input_labels that do not match"
                " the expected columns in the input_scaler.\n"
                "input_labels={}\n"
                "input_scaler.expected_columns()={}".format(
                    input_labels, input_scaler.expected_columns()
                )
            )

        # check that the output labels match
        if (
            output_scaler is not None
            and output_scaler.expected_columns() != output_labels
        ):
            raise ValueError(
                "KerasSurrogate created with output_labels that do not match"
                " the expected columns in the output_scaler.\n"
                "output_labels={}\n"
                "output_scaler.expected_columns()={}".format(
                    output_labels, output_scaler.expected_columns()
                )
            )

        self._input_scaler = input_scaler
        self._output_scaler = output_scaler

    class Formulation(Enum):
        FULL_SPACE = 1
        REDUCED_SPACE = 2
        RELU_BIGM = 3
        RELU_COMPLEMENTARITY = 4

    def generate_omlt_scaling_objecets(self):
        offset_inputs = np.zeros(self.n_inputs())
        factor_inputs = np.ones(self.n_inputs())
        offset_outputs = np.zeros(self.n_outputs())
        factor_outputs = np.ones(self.n_outputs())
        if self._input_scaler:
            offset_inputs = self._input_scaler.offset_series()[
                self.input_labels()
            ].to_numpy()
            factor_inputs = self._input_scaler.factor_series()[
                self.input_labels()
            ].to_numpy()
        if self._output_scaler:
            offset_outputs = self._output_scaler.offset_series()[
                self.output_labels()
            ].to_numpy()
            factor_outputs = self._output_scaler.factor_series()[
                self.output_labels()
            ].to_numpy()

        omlt_scaling = OffsetScaling(
            offset_inputs=offset_inputs,
            factor_inputs=factor_inputs,
            offset_outputs=offset_outputs,
            factor_outputs=factor_outputs,
        )

        # omlt takes *scaled* input bounds as a dictionary with int keys
        input_bounds = dict(enumerate(self.input_bounds().values()))
        scaled_input_bounds = omlt_scaling.get_scaled_input_expressions(input_bounds)
        scaled_input_bounds = {i: tuple(bnd) for i, bnd in scaled_input_bounds.items()}
        return omlt_scaling, scaled_input_bounds

    def populate_block_with_net(self, block, formulation_object):
        """
        Method to populate a Pyomo Block with the omlt model constraints and build its formulation.

        Args:
           block: Pyomo Block component
              The block to be populated with variables and/or constraints.
           formulation_object: omlt loaded network formulation
        """

        block.nn = OmltBlock()
        block.nn.build_formulation(
            formulation_object,
        )

        # input/output variables need to be constrained to be equal
        # auto-created variables that come from OMLT.
        input_idx_by_label = {s: i for i, s in enumerate(self._input_labels)}
        input_vars_as_dict = block.input_vars_as_dict()

        @block.Constraint(self._input_labels)
        def input_surrogate_ties(m, input_label):
            return (
                input_vars_as_dict[input_label]
                == block.nn.inputs[input_idx_by_label[input_label]]
            )

        output_idx_by_label = {s: i for i, s in enumerate(self._output_labels)}
        output_vars_as_dict = block.output_vars_as_dict()

        @block.Constraint(self._output_labels)
        def output_surrogate_ties(m, output_label):
            return (
                output_vars_as_dict[output_label]
                == block.nn.outputs[output_idx_by_label[output_label]]
            )
