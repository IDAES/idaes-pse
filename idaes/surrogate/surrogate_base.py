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
from pyomo.common.config import ConfigBlock

from idaes.surrogate.metrics import compute_fit_metrics


class SurrogateTrainer(object):
    CONFIG = ConfigBlock()

    def __init__(self, input_labels, output_labels,
                 training_dataframe, validation_dataframe=None,
                 input_bounds=None, **settings):
        """
        This is the base class for IDAES surrogate training objects.

        Args:
           input_labels: list
              list of labels corresponding to the inputs (in order)
           output_labels: list
              list of labels corresponding to the outputs (in order)
           input_bounds: None, or dict of tuples
              if None, these are set later from the provided data
              if provided, it should be a dictionary where the keys correspond
              to the input label, and the values are tuples of bounds
              (lower, upper)
           training_dataframe: pandas DataFrame
              Pandas DataFrame corresponding to the training data. Columns must
              include all the labels in input_labels and output_labels
           validation_dataframe: pandas DataFrame or None
             Pandas DateFrame corresponding to the validation data. Columns
             must include all the labels in input_labels and output_labels. If
             None is passed, then no validation data will be used. Some
             derived surrogate trainers may require validation data, while
             others may not.
        """
        # Set the config block from passed settings
        self.config = self.CONFIG(settings)

        # We must have at least one input label and one output label
        if input_labels is None or len(input_labels) < 1 \
           or output_labels is None or len(output_labels) < 1:
            raise ValueError(
                'SurrogateTrainer requires a list of input_labels and a list '
                'of output_labels which must both have a length of at '
                'least one')

        self._input_labels = list(input_labels)
        self._output_labels = list(output_labels)

        # check that the input and output labels do not overlap
        all_labels = set(self._input_labels)
        all_labels.update(self._output_labels)
        if len(all_labels) != (
                len(self._input_labels) + len(self._output_labels)):
            raise ValueError(
                'Duplicate label found in input_labels and/or output_labels.')

        # create the data members for training and validation data
        self._training_dataframe = training_dataframe
        self._validation_dataframe = validation_dataframe

        # check that all input labels and output labels are in the dataframes
        if set(self._input_labels) - set(self._training_dataframe.columns):
            raise ValueError('An input label was specified that was not '
                             'found in the training data.')
        if (self._validation_dataframe is not None and
                set(self._input_labels) -
                set(self._validation_dataframe.columns)):
            raise ValueError('An input label was specified that was not '
                             'found in the validation data.')
        if set(self._output_labels) - set(self._training_dataframe.columns):
            raise ValueError('An output label was specified that was not '
                             'found in the training data.')
        if (self._validation_dataframe is not None and
                set(self._output_labels) -
                set(self._validation_dataframe.columns)):
            raise ValueError('An output label was specified that was not '
                             'found in the validation data.')

        if input_bounds is not None:
            self._input_bounds = dict(input_bounds)
        else:
            # get the bounds from the data
            mx = self._training_dataframe.max().to_dict()
            mn = self._training_dataframe.min().to_dict()
            self._input_bounds = {
                k: (mn[k], mx[k]) for k in self._input_labels}

    def n_inputs(self):
        return len(self._input_labels)

    def n_outputs(self):
        return len(self._output_labels)

    def input_labels(self):
        return self._input_labels

    def output_labels(self):
        return self._output_labels

    def input_bounds(self):
        return self._input_bounds

    def train_surrogate(self):
        """
        This method should be overridden by the derived classes.

        The ``train_surrogate`` method is used to train a surrogate model
        using data provided in set_training_data. This method should return an
        instance of a derived surrogate object (from SurrogateBase)

        Returns:
           tuple : (bool, surrogate object, message) where bool indicates
           status of model training, surrogate object is an instance of a
           class derived from SurrogateBase, and message is string containing
           additional information from the trainer.

        """
        raise NotImplementedError(
            'train_surrogate called, but not implemented on the derived class')


class SurrogateBase():
    def __init__(
            self, input_labels=None, output_labels=None, input_bounds=None):
        """
        Base class for standard IDAES Surrogate object. This class is
        responsible for being able to load/save a surrogate model, evaluate
        the model given an input dataframe, and populating a block to provide
        an EO representation of the surrogate for solving in IDAES.

        Args:
           input_labels: list
              list of labels corresponding to the inputs (in order)
           output_labels: list
              list of labels corresponding to the outputs (in order)
           input_bounds: dict of tuples
              A dictionary where the keys correspond to the input label,
              and the values are tuples of bounds (lower,upper). These
              should represent the valid range for the input variables
        """
        self._input_labels = input_labels
        self._output_labels = output_labels
        self._input_bounds = input_bounds

    def n_inputs(self):
        return len(self._input_labels)

    def n_outputs(self):
        return len(self._output_labels)

    def input_labels(self):
        return self._input_labels

    def output_labels(self):
        return self._output_labels

    def input_bounds(self):
        return self._input_bounds

    def populate_block(self, block, **kwargs):
        """
        Method to populate a Pyomo Block with surrogate model
        constraints and variables.

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

    def evaluate_surrogate(self, dataframe):
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

    # TODO: this should serialize to a stream instead of a file
    def save(self, filename):
        """
        Save an instance of this surrogate to be used in a model later
        """
        raise NotImplementedError('"save" should be implemented in the derived'
                                  ' SurrogateObject class')

    @staticmethod
    def load(self, filename):
        """
        Load an instance of this surrogate from a file
        """
        raise NotImplementedError('"load" should be implemented in the derived'
                                  ' SurrogateObject class')

    def compute_fit_metrics(self, data):
        """
        This method computes a variety of metrics regarding the fit of the
        surrogate to the data provided.

        Args:
           data : pandas DataFrame
              pandas DataFrame that includes columns for the inputs and the
              outputs. Metrics for the quality of fit will be computed across
              all the rows in the provided data.

        Returns:
           TrainingMetrics object
        """
        return compute_fit_metrics(self, data)
