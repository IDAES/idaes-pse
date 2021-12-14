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
Interface for importing Keras models into IDAES
"""
import numpy as np
import json
import os.path
import keras
from idaes.surrogate.base.surrogate_base import SurrogateBase
from idaes.surrogate.sampling.scaling import StandardScaler
from omlt import OmltBlock, OffsetScaling
from omlt.neuralnet import FullSpaceContinuousFormulation, ReducedSpaceContinuousFormulation, \
    ReLUBigMFormulation, load_keras_sequential


class KerasSurrogate(SurrogateBase):
    """
    Standard SurrogateObject for surrogates based on Keras models. 
    Utilizes the OMLT framework for importing Keras models to IDAES.

    Contains methods to both populate a Pyomo Block with constraints
    representing the surrogate and to evaluate the surrogate a set of user
    provided points.
    """
    def __init__(self, keras_model, input_labels, output_labels, input_bounds,
                 input_scaler=None, output_scaler=None):
        super().__init__(input_labels=input_labels, output_labels=output_labels, \
                         input_bounds=input_bounds)

        # make sure we are using the standard scaler
        if input_scaler is not None and type(input_scaler) is not StandardScaler or \
           output_scaler is not None and type(output_scaler) is not StandardScaler:
            raise NotImplementedError('KerasSurrogate only supports the StandardScaler.')

        # check that the input labels match
        if input_scaler is not None and input_scaler.expected_columns() != input_labels:
            raise ValueError('KerasSurrogate created with input_labels that do not match'
                             ' the expected columns in the input_scaler.\n'
                             'input_labels={}\n'
                             'input_scaler.expected_columns()={}'.format(
                                 input_labels,
                                 input_scaler.expected_columns()))

        # check that the output labels match
        if output_scaler is not None and output_scaler.expected_columns() != output_labels:
            raise ValueError('KerasSurrogate created with output_labels that do not match'
                             ' the expected columns in the output_scaler.\n'
                             'output_labels={}\n'
                             'output_scaler.expected_columns()={}'.format(
                                 output_labels,
                                 output_scaler.expected_columns()))

        self._input_scaler = input_scaler
        self._output_scaler = output_scaler
        self._keras_model = keras_model
        

    def populate_block(self, block, **kwargs):
        """
        Method to populate a Pyomo Block with the keras model constraints.

        Args:
           block: Pyomo Block component
              The block to be populated with variables and/or constraints.
           formulation: string
              The formulation to use with OMLT. Possible values are 'full-space', 
              'reduced-space', 'relu-bigm', or 'relu-complementarity' (default is 'full-space')
        """
        formulation = kwargs.pop('formulation', 'full-space')
        offset_inputs = np.zeros(self.n_inputs())
        factor_inputs = np.ones(self.n_inputs())
        offset_outputs = np.zeros(self.n_outputs())
        factor_outputs = np.ones(self.n_outputs())
        if self._input_scaler:
            offset_inputs = self._input_scaler.offset_series()[self.input_labels()].to_numpy()
            factor_inputs = self._input_scaler.factor_series()[self.input_labels()].to_numpy()
        if self._output_scaler:
            offset_outputs = self._output_scaler.offset_series()[self.output_labels()].to_numpy()
            factor_outputs = self._output_scaler.factor_series()[self.output_labels()].to_numpy()

        omlt_scaling = OffsetScaling(offset_inputs=scaling_x_offset,
                                      factor_inputs=scaling_x_factor,
                                      offset_outputs=scaling_y_offset,
                                      factor_outputs=scaling_y_factor)
        net = load_keras_sequential(self._keras_model, omlt_scaling, self.input_bounds())

        if formulation == 'full-space':
            formulation_object = FullSpaceContinuousFormulation(net)
        elif formulation == 'reduced-space':
            formulation_object = ReducedSpaceContinuousFormulation(net)
        elif formulation == 'relu-bigm':
            formulation_object = ReLUBigMFormulation(net)
        elif formulation == 'relu-complementarity':
            formulation_object = ReLUComplementarityFormulation(net)

        block.neural_network.build_formulation(formulation_object,
                                               input_vars=block._input_vars_as_list(),
                                               output_vars=block._output_vars_as_list())



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
        y = pd.Dataframe(data=y, columns=self.output_labels(), index=inputs.index)

        if self._output_scaler is not None:
            y = self._output_scaler.unscale(y)

        return y

    @classmethod
    def save_to_file(cls, keras_folder_name):
        """
        Save the surrogate object to disk by providing the name of the
        folder to contain the keras model and additional IDAES metadata

        Args:
           folder_name: str
              The name of the folder to contain the Keras model and additional
              IDAES metadata
        """
        self._keras_model.save(keras_folder_name)
        info = dict()
        info['input_scaler'] = None
        if self._input_scaler is not None:
            info['input_scaler'] = self._input_scaler.to_dict()
        info['output_scaler'] = None
        if self._output_scaler is not None:
            info['output_scaler'] = self._output_scaler.to_dict()
        
        # serialize information from the base class
        info['input_labels'] = self.input_labels()
        info['output_labels'] = self.output_labels()
        info['input_bounds'] = self.input_bounds()

        with open(os.path.join(keras_folder_name, 'idaes_info.json'), 'w') as fd:
            json.dump(info, fd)

    @classmethod
    def load_from_file(cls, keras_folder_name):
        """
        Load the surrogate object from disk by providing the name of the
        folder holding the keras model

        Args:
           folder_name: str
              The name of the folder containing the Keras model and additional
              IDAES metadata

        Returns: an instance of KerasSurrogate
        """
        keras_model = keras.models.load_model(keras_folder_name)
        with open(os.path.join(keras_folder_name, 'idaes_info.json')) as fd:
            info = json.load(fd)

        input_scaler = None
        if info['input_scaler'] is not None:
            input_scaler = StandardScaler.from_dict(info['input_scaler'])

        output_scaler = None
        if info['output_scaler'] is not None:
            output_scaler = StandardScaler.from_dict(info['output_scaler'])

        return KerasSurrogate(keras_model=keras_model,
                              input_labels=info['input_labels'],
                              output_labels=info['output_labels'],
                              input_bounds=info['input_bounds'],
                              input_scaler=input_scaler,
                              output_scaler=output_scaler)
    
        

    
