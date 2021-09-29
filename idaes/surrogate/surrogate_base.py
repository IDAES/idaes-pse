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
Common Surrogate interface for IDAES.
"""
from pyomo.environ import Var
from pyomo.common.config import ConfigBlock, ConfigValue, ConfigList
from pyomo.core.base.global_set import UnindexedComponent_set
import os.path, pickle


class SurrogateTrainer:
    CONFIG = ConfigBlock()

    # TODO: self._surrogate is *not* a pyomo expression
    def __init__(self, input_labels, output_labels, input_bounds=None, **settings):
        """
        TODO : document this method
        Args:
           input_labels: list
              list of labels corresponding to the inputs (in order)
           output_labels: list
              list of labels corresponding to the outputs (in order)
           input_bounds: None, or dict of tuples
              if None, these are set later from the provided data
              if provided, it should be a dictionary where the keys correspond
              to the input label, and the values are tuples of bounds (lower,upper)
           settings: keyword arguments
              configuration options for the derived class
        """

        # Set the config block from passed settings
        self.config = self.CONFIG(settings)

        # Objects known to the base class
        self._training_status = dict()
        self._training_metrics = dict()

        # TODO: make these exceptions
        assert input_labels is not None
        assert output_labels is not None
        self._input_labels = input_labels
        self._output_labels = output_labels
        self._input_bounds = input_bounds

        self._training_data_in = None
        self._training_data_out = None
        self._validation_data_in = None
        self._validation_data_out = None

    def n_inputs(self):
        return len(self._input_labels)

    def n_outputs(self):
        return len(self._output_labels)

    def train_surrogate(self):
        """
        This method should be overridden by the derived classes.

        The ``train_surrogate`` method trains a surrogate model to an input dataset.
        It calls the core method which is called during surrogate generation: ``train_surrogate`` sets up the surrogate problem,
        trains the surrogate, computes the metrics, creates a results object and generates the Pyomo representation of the model.
        It accepts no user input, inheriting the information passed in class initialization.
        """
        raise NotImplementedError('train_surrogate called, but not implemented on the derived class')
    

    def get_surrogate(self):  # SurrogateObject of the appropriate derived class
        """
        The ``get_surrogate`` method returns the result of the surrogate training process as as an IDAES surrogate object of the appropriate derived class.

        This surrogate object can be used to evaluate the surrogate, build an IDAES 
        block, or save the surrogate for use later

        Returns:
            derived from SurrogateBase
        """
        raise NotImplementedError('get_surrogate called, but not implemented on the derived class')

    # TODO: methods to report key metrics

    # Data Handling
    def set_input_labels(self, labels):
        # TODO: argument validation, docs and tests
        self._input_labels = labels

    def set_output_labels(self, labels):
        # TODO: argument validation, docs and tests
        self._output_labels = labels

    def set_input_bounds(self, bounds):
        self._input_bounds = bounds

    def get_training_data(self):
        """
        The ``get_training_data`` method returns the input data used in training the surrogate model.
        Returns:
            Tuple : Tuple of two elements containing input (samples) and output (output values).

        """
        return (self._training_data_in, self._training_data_out)

    def set_training_data(self, input_data, output_data):  # 2D Numparray
        """
        The ``set_training_data`` method initializes the Surrogate class with the data for training the surrogate model.
        Args:
            input_data (NumPy Array) : Two-dimensional NumPy Array containing the samples/features.
            output_data (NumPy Array)  : Two-dimensional NumPy Array containing the output values.
        """
        # TODO: Support Pandas
        # TODO: Testing of data shape based on labels
        self._training_data_in = input_data
        self._training_data_out = output_data

    def get_validation_data(self):
        """
        The ``get_validation_data`` method returns the data supplied for validating the surrogate model.
        Returns:
            Tuple : Tuple of two elements containing validation data samples and output values.

        """
        return (self._validation_data_in, self._validation_data_out)

    def set_validation_data(self, validation_inputs, validation_outputs):  # 2D Numparray
        """
        The ``set_validation_data`` method initializes the Surrogate class with data for validating/testing the surrogate model after generation.
        Args:
            validation_inputs (NumPy Array) : Two-dimensional NumPy Array containing the validation samples/features.
            validation_outputs (NumPy Array) : Two-dimensional NumPy Array containing the output values of the validation samples.
        """
        # TODO: Support Pandas
        # TODO: Testing of data shape based on labels
        self._validation_data_in = validation_inputs
        self._validation_data_out = validation_outputs

    # TODO: Method to get validation metrics

    # TODO: Methods to save and load SurrogateTrainers


class SurrogateBase():
    """
    Base class for standard IDAES Surrogate Object
    """

    def __init__(self, surrogate, input_labels=None, output_labels=None,
                 input_bounds=None):
        self._surrogate = surrogate
        self._input_labels = input_labels
        self._output_labels = output_labels
        self._input_bounds = input_bounds  # dict of bounds for each label

    @property
    def n_inputs(self):
        return len(self._input_labels)

    @property
    def n_outputs(self):
        return len(self._output_labels)

    @property
    def input_labels(self):
        return self._input_labels

    @property
    def output_labels(self):
        return self._output_labels

    @property
    def input_bounds(self):
        return self._input_bounds

    def populate_block(self, block, **kwargs):
        """
        Placeholder method to populate a Pyomo Block with surrogate model
        constraints.

        Args:
            block: Pyomo Block component to be populated with constraints.
            variables: dict mapping surrogate variable labels to existing
                Pyomo Vars (default=None). If no mapping provided,
                construct_variables will be called to create a set of new Vars.
            index_set: (optional) if provided, this will be used to index the
                Constraints created. This must match the indexing Set of the
                Vars provided in the variables argument.

        Returns:
            None
        """
        raise NotImplementedError(
            "SurrogateModel class has not implemented populate_block method.")

    def evaluate_surrogate(self, inputs):
        """
        Placeholder method to evaluate surrogate model at a set of user
        provided values.

        Args:
            inputs: numpy array of input values

        Returns:
            output: numpy array of output values evaluated at inputs
        """
        raise NotImplementedError(
            "SurrogateModel class has not implemented an evaluate_surrogate "
            "method.")

    # todo: this should serialize to a stream
    def save(self, filename):
        """
        Save an instance of this surrogate to be used in a model later
        """
        raise NotImplementedError('"save" should be implemented in the derived'
                                  ' SurrogateObject class')

    # TODO: it is recommended that you add a "load" static method to build
    #       a derived surrogate object from the file on disk
    def load(self, filename):
        """
        Load an instance of this surrogate from a file
        """
        raise NotImplementedError('"load" should be implemented in the derived'
                                  ' SurrogateObject class')
