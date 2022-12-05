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


class SurrogateTrainer(object):
    CONFIG = ConfigBlock()

    def __init__(
        self,
        input_labels,
        output_labels,
        training_dataframe,
        validation_dataframe=None,
        input_bounds=None,
        **settings
    ):
        """
        This is the base class for IDAES surrogate training objects.

        Args:
           input_labels: list
              list of labels corresponding to the inputs (in order)
           output_labels: list
              list of labels corresponding to the outputs (in order)
           training_dataframe: pandas DataFrame
              Pandas DataFrame corresponding to the training data. Columns must
              include all the labels in input_labels and output_labels
           validation_dataframe: pandas DataFrame or None
             Pandas DateFrame corresponding to the validation data. Columns
             must include all the labels in input_labels and output_labels. If
             None is passed, then no validation data will be used. Some
             derived surrogate trainers may require validation data, while
             others may not.
           input_bounds: None, or dict of tuples
              if None, these are set later from the provided data
              if provided, it should be a dictionary where the keys correspond
              to the input label, and the values are tuples of bounds
              (lower, upper)
           settings: additional keyword arguments
              These are additional keyword arguments that are passed to the CONFIG
              for the derived class.
        """
        # Set the config block from passed settings
        self.config = self.CONFIG(settings)

        # We must have at least one input label and one output label
        if (
            input_labels is None
            or len(input_labels) < 1
            or output_labels is None
            or len(output_labels) < 1
        ):
            raise ValueError(
                "SurrogateTrainer requires a list of input_labels and a list "
                "of output_labels which must both have a length of at "
                "least one"
            )

        self._input_labels = list(input_labels)
        self._output_labels = list(output_labels)

        # check that the input and output labels do not overlap
        all_labels = set(self._input_labels)
        all_labels.update(self._output_labels)
        if len(all_labels) != (len(self._input_labels) + len(self._output_labels)):
            raise ValueError(
                "Duplicate label found in input_labels and/or output_labels."
            )

        # create the data members for training and validation data
        self._training_dataframe = training_dataframe
        self._validation_dataframe = validation_dataframe

        # check that all input labels and output labels are in the dataframes
        diff = set(self._input_labels) - set(self._training_dataframe.columns)
        if diff:
            raise ValueError(
                "The following input labels were not found in "
                "the training data columns: {}.".format(diff)
            )
        if self._validation_dataframe is not None:
            diff = set(self._input_labels) - set(self._validation_dataframe.columns)
            if diff:
                raise ValueError(
                    "The following input labels were not found in "
                    "the validation data columns: {}.".format(diff)
                )

        diff = set(self._output_labels) - set(self._training_dataframe.columns)
        if diff:
            raise ValueError(
                "The following output labels were not found in "
                "the training data columns: {}.".format(diff)
            )

        if self._validation_dataframe is not None:
            diff = set(self._output_labels) - set(self._validation_dataframe.columns)
            if diff:
                raise ValueError(
                    "The following output labels were not found in "
                    "the validation data columns: {}.".format(diff)
                )

        if input_bounds is not None:
            self._input_bounds = dict(input_bounds)
            # check that the labels match the input labels
            if sorted(input_bounds.keys()) != sorted(self._input_labels):
                raise ValueError(
                    "The input_labels did not match the keys in input_bounds.\n"
                    "input_bounds.keys(): {}\n"
                    "input_labels: {}".format(
                        sorted(input_bounds.keys()), sorted(input_labels)
                    )
                )
        else:
            # get the bounds from the data
            mx = self._training_dataframe.max().to_dict()
            mn = self._training_dataframe.min().to_dict()
            self._input_bounds = {k: (mn[k], mx[k]) for k in self._input_labels}

    def n_inputs(self):
        """
        The number of inputs for the surrogate

        Returns: float
        """
        return len(self._input_labels)

    def n_outputs(self):
        """
        The number of outputs for the surrogate

        Returns: float
        """
        return len(self._output_labels)

    def input_labels(self):
        """
        The ordered list of labels for the inputs

        Returns: list of strings
        """
        return list(self._input_labels)

    def output_labels(self):
        """
        The ordered list of labels for the outputs

        Returns: list of strings
        """
        return list(self._output_labels)

    def input_bounds(self):
        """
        The dictionary of input bounds. The keys of the dictionary correspond
        to the labels for the inputs. The values are tuples of (lower_bound, upper_bound)

        Returns: dict
        """
        if self._input_bounds:
            return dict(self._input_bounds)
        return None

    def train_surrogate(self):
        """
        The ``train_surrogate`` method is used to train a surrogate model
        using data provided in set_training_data. This method should return an
        instance of a derived surrogate object (from SurrogateBase)

        This method should be overridden by the derived classes.

        Returns:
           tuple : (bool, surrogate object, message) where bool indicates
           status of model training, surrogate object is an instance of a
           class derived from SurrogateBase, and message is string containing
           additional information from the trainer.

        """
        raise NotImplementedError(
            "train_surrogate called, but not implemented on the derived class"
        )


class SurrogateBase:
    def __init__(self, input_labels=None, output_labels=None, input_bounds=None):
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

        # check that the input and output labels do not overlap
        all_labels = set(self._input_labels)
        all_labels.update(self._output_labels)
        if len(all_labels) != (len(self._input_labels) + len(self._output_labels)):
            raise ValueError(
                "Duplicate label found in input_labels and/or output_labels."
            )

        if input_bounds is not None:
            self._input_bounds = dict(input_bounds)
            # check that the labels match the input labels
            if sorted(input_bounds.keys()) != sorted(self._input_labels):
                raise ValueError(
                    "The input_labels did not match the keys in input_bounds.\n"
                    "input_bounds.keys(): {}\n"
                    "input_labels: {}".format(
                        sorted(input_bounds.keys()), sorted(input_labels)
                    )
                )

    def n_inputs(self):
        """
        The number of inputs for the surrogate

        Returns: float
        """
        return len(self._input_labels)

    def n_outputs(self):
        """
        The number of outputs for the surrogate

        Returns: float
        """
        return len(self._output_labels)

    def input_labels(self):
        """
        The ordered list of labels for the inputs

        Returns: list of strings
        """
        return list(self._input_labels)

    def output_labels(self):
        """
        The ordered list of labels for the outputs

        Returns: list of strings
        """
        return list(self._output_labels)

    def input_bounds(self):
        """
        The dictionary of input bounds. The keys of the dictionary correspond
        to the labels for the inputs. The values are tuples of (lower_bound, upper_bound)

        Returns: dict
        """
        if self._input_bounds:
            return dict(self._input_bounds)
        return None

    def populate_block(self, block, additional_options=None):
        """
        Method to populate a Pyomo Block with surrogate model
        constraints and variables.

        Derived classes must overload this method.

        Args:
            block: Pyomo Block
               Component to be populated with constraints.
            additional_options: dict
               Additional options passed through from SurrogateBlock.build_model

        Returns:
            None
        """
        raise NotImplementedError(
            "SurrogateModel class has not implemented populate_block method."
        )

    def evaluate_surrogate(self, dataframe):
        """
        Method to method to evaluate surrogate model at a set of user
        provided values.

        Derived classes must overload this method

        Args:
           dataframe: pandas DataFrame
              The dataframe of input values to be used in the evaluation. The dataframe
              needs to contain a column corresponding to each of the input labels. Additional
              columns are fine, but are not used.

        Returns:
            output: pandas Dataframe
              Returns a dataframe of the the output values evaluated at the provided inputs.
              The index of the output dataframe should match the index of the provided inputs.
        """
        raise NotImplementedError(
            "SurrogateModel class has not implemented an evaluate_surrogate " "method."
        )

    def save_to_file(self, filename, overwrite=False):
        """
        This method saves an instance of the surrogate to a file so the model
        can be used later.

        Args:
           filename : str
              The path of the filename where the model will be saved
           overwrite : bool
              If True, this method will overwrite the file if it exists. If False
              and the file already exists, this will throw an error.
        """
        arg = "x"
        if overwrite:
            arg = "w"

        with open(filename, arg) as fd:
            self.save(fd)

    def save(self, strm):
        """
        Save an instance of this surrogate to the strm so the model can be used later.
        This method should be overloaded in derived surrogate classes.

        Args:
           strm: IO.TextIO
              This is the python stream like a file object or StringIO that will be used
              to serialize the surrogate object. This methods will often write a string
              of json data to the stream, but hte format for derived classes need not be json.
        """
        raise NotImplementedError(
            '"save" should be implemented in the' " class derived from SurrogateBase"
        )

    @classmethod
    def load_from_file(cls, filename):
        """
        This method creates a new surrogate object by loading the model from the provided file.

        Args:
           filename : str
              The name of the file from which to load the model.

        Returns: an instance of the derived class or None if it failed to load
        """
        with open(filename, "r") as fd:
            return cls.load(fd)

    @classmethod
    def load(cls, strm):
        """
        Create an instance of a surrogate from a stream. This method should
        be overloaded in derived surrogate classes.

        Args:
           strm: IO.TextIO
              This is the python stream like a file object or StringIO containing
              the data required to load the surrogate. This is often, but does not
              need to be json data.

        Returns: an instance of the derived class or None if it failed to load
        """
        raise NotImplementedError(
            '"load" should be implemented in the' " class derived from SurrogateBase"
        )
