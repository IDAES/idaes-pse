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
from pathlib import Path
from typing import Dict
import yaml
from pyomo.common.config import ConfigBlock, ConfigValue, ConfigList
import os.path, pickle


# from mypy_extensions import TypedDict

class Metrics:
    """

    Names for known types of metrics.

    Use these as keys in dictionaries, e.g.:
    m = {Metrics.RMSE: self.rmse}
    When adding attributes to this class, please include a comment with the
    prefix "#:" immediately above it, so Sphinx knows that this is documentation
    for the attribute.

    """

    #: Root mean-squared error
    RMSE = "RMSE"

    #: Mean-squared error
    MSE = "MSE"

    #: Sum of squared error
    SSE = "SSE"

    #: Time
    Time = "Time"

    #: Order
    Order = "Order"

    #: R-squared
    R2 = "R2"


class ConfigurationError(Exception):
    """
    IDAES exception to be used when configuration arguments are incorrect
    or inconsistent.
    """
    pass


# Single Surrogate Modeler
class Surrogate:
    CONFIG = ConfigBlock()

    # Common Declarations
    # CONFIG.declare('overwrite', ConfigValue(default=None, domain=bool))
    # CONFIG.declare('fname', ConfigValue(default=None, domain=str))

    def __init__(self, **settings):
        """
        Initialization fo the Surrogate Class.
        Keyword Args:
            settings            : Dictionary of user-defined configurations and settings for the selected surrogate model tool(s).
        Returns:
            *self** object containing all the input information and run results, including:
                - **self.config** (Dict)                                    : The configuration of the selected tool
                - **self._results** (Tuple)                                 : Performance metrics of the surrogate(s) trained
                - **self._model**   (Pyomo Expression)                      : Pyomo representation of resulting surrogate
                - **self._r_data_in**, **self._r_data_out** (NumPy Array)   : Sample points and output values used in training the surrogate
                - **self._v_data_in**, **self._v_data_out** (NumPy Array)   : Validation/test sample points and their true output values
                - **self.pkl_info** (Python Object)                         : Python object containing relevant surrogate model information.
        """

        # Config
        self.config = self.CONFIG(settings)
        self.modeler = None

        # Results
        self._results = {}
        self._metrics = None
        self._model = None
        self._b_built = False  # flag for regression

        # Data
        self._rdata_in = None
        self._rdata_out = None
        self._vdata_in = None
        self._vdata_out = None

        self.pkl_info = None

    def modify_config(self, **kwargs):
        """
        The ``modify_config`` method allows users to define a new surrogate instance simply by modifying one or more of the
        settings or keywords of a previously-defined surrogate instance.
        The values of the new keywords defined in **settings will directly replace the previous values stored in self.CONFIG.
        All other keywords remain unchanged.
        Args:
            **kwargs (dict)             : Dictionary containing (key, val) entries for the settings to be modified in self.CONFIG
        """
        _b_built = False
        self.config = self.CONFIG(kwargs)

    # Build

    def build_model(self):
        """
        The ``build_model`` method trains a surrogate model to an input dataset.
        It calls the core method which is called during surrogate generation: ``build_model`` sets up the surrogate problem,
        trains the surrogate, computes the metrics, creates a results object and generates the Pyomo representation of the model.
        It accepts no user input, inheriting the information passed in class initialization.
        """
        self._b_built = True

        self.pkl_info = {'In data': self._rdata_in.tolist(),
                         'Out data': self._rdata_out.tolist()}

        pass

    def update_model(self):
        """
        The ``update_model`` trains a new surrogate model based on the updated configuration/set-up defined by
        calling ``modify_config``
        It accepts no user input, inheriting the information provided in ``modify_config``.
        """
        self.build_model()
        pass

    # def generate_expression(self, variable_list): # TODO
    #     pass

    # Get Results

    def get_model(self):  # Pyomo Expression
        """
        The ``get_model`` method returns the result of the surrogate training process as a Pyomo Expression
        Returns:
            Pyomo Expression    : Pyomo expression of surrogate model trained.
        """
        return self._model


    def get_results(self):  # Metrics Object
        """
        The ``get_results`` method returns the performance metrics of the surrogate model(s) generated.
        Returns:
            Metrics of the surrogates based on the training data, including but not limited to the:
                - root mean squared error (RMSE),
                - mean squared error (MSE),
                - :math:`R^{2}`of the surrogate fit, and
                - model run time.
        """
        return self._results

    # Data Handling

    def get_regressed_data(self):
        """
        The ``get_regressed_data`` method returns the input data used in regressing the surrogate model.
        Returns:
            Tuple        : Tuple of two elements containing r_in (samples) and r_out (output values).

        """
        return (self._rdata_in, self._rdata_out)

    def regressed_data(self, r_in, r_out):  # 2D Numparray
        """
        The ``regressed_data`` method initializes the Surrogate class with the data for training the surrogate model.
        Args:
            r_in  (NumPy Array)  : Two-dimensional NumPy Array containing the samples/features.
            r_out (NumPy Array)  : Two-dimensional NumPy Array containing the output values.
        """
        self._rdata_in = r_in
        self._rdata_out = r_out

    def get_validation_data(self):
        """
        The ``get_validation_data`` method returns the data supplied for validating the surrogate model.
        Returns:
            Tuple        : Tuple of two elements containing validation data v_in (samples) and v_out (output values).

        """
        return (self._vdata_in, self._vdata_out)

    def validation_data(self, v_in, v_out):  # 2D Numparray
        """
        The ``validation_data`` method initializes the Surrogate class with data for validating/testing the surrogate model after generation.
        Args:
            v_in  (NumPy Array)  : Two-dimensional NumPy Array containing the validation samples/features.
            v_out (NumPy Array)  : Two-dimensional NumPy Array containing the output values of the validation samples.
        """
        self._vdata_in = v_in
        self._vdata_out = v_out

    # Using regressed model
    # PYLINT-TODO: check if adding self as arg to fix pylint "undefined-variable 'self'" is valid
    def calculate_outputs(self, inputs):  # 2D Numparray, use pyomo expression
        """
        ``calculate_outputs`` evaluates the output predictions from the surrogate for an array of input samples **inputs**
        Args:
            inputs(NumPy Array)     : Two-dimensional NumPy Array containing the sample points to be evaluated.
        Returns:
            outputs(NumPy Array)    : NumPy Array containing the output predictions from the surrogate model.
        """
        outputs = self._model(inputs)
        return outputs

    # Additional Metrics - MUST NOT OVERWRITE MODELER METRICS

    # def get_trained_metrics():  # 2D Nparray, 2D Numparray # TODO
    #     if not _b_built:
    #         print("Warning: No surrogate model regressed.")
    #         return
    #     pass

    # PYLINT-TODO: check if adding as an argument and using self._b_built in the body is valid
    def get_validated_metrics(self, xval, zval):  # 2D, 2D Numparray, use pyomo expression
        """
        ``get_validated_metrics`` evaluates the performance metrics for the surrogate model based on a set of off-design points (xval, zval)
        Args:
            xval (NumPy Array)      : Two-dimensional array containing a set of off-design samples.
            zval (NumPy Array)      : Array containing a true output values at off-design sample points.
        Returns:
            Metrics of the surrogate model evaluated based on the off-design data set (xval, zval), including but not limited to the:
                - root mean squared error (RMSE),
                - mean squared error (MSE), and
                - :math:`R^{2}`of the surrogate fit based on the off-design data points.
        """
        if not self._b_built:
            print("Warning: No surrogate model regressed.")
            return
        pass

    def save_results(self, filename, overwrite=False):
        """
        The ``save_results`` method saves the results of the surrogate run in a pickle object
        Args:
            filename (str)      : The name of the file to be saved. Must be of extension .pickle
            overwrite (bool)    : Boolean controlling whether any existing file with the same name is overwritten or not. Default is False.
        Raises:
            Exception:
                * **filename** is not a string or has the wrong extension type
             Exception:
                * A file with the name **filename** already exists and **overwrite** is False
             Exception:
                * A problem is encountered while trying to save the file.
        """
        # Ensure overwrite option, when entered, is boolean
        if not isinstance(overwrite, bool):
            raise Exception('overwrite must be boolean.')
        # Check if filename is a string with pickle extension
        if not isinstance(filename, str) or os.path.splitext(filename)[-1].lower() != '.pickle':
            raise Exception('filename must be a string with extension ".pickle". Please correct.')
        # If overwite is false, throw up error if the filename already exists in the destination folder
        if os.path.exists(filename) and overwrite is False:
            raise Exception(filename, 'already exists!.\n')
        # Try to save
        try:
            filehandler = open(filename, 'wb')
            pickle.dump(self.pkl_info, filehandler)
            print('\nResults saved in ', str(filename))
        except:
            raise Exception('File could not be saved.')

    def load_results(self, filename):
        """
        ``load_results`` loads the results of a saved run 'filename.obj'.
        Args:
            filename(str)       : Pickle object file containing previous solution to be loaded.
        Returns:
            **self** object containing set-up information and results of the run that created **filename**.
        Raises:
            Exception:
                * **filename** does not exist in the directory.
            Exception:
                * A problem is encountered while trying to load **filename**.
        """
        if os.path.exists(filename) is False:
            raise Exception(filename, 'does not exist in directory.')
        try:
            filehandler = open(filename, 'rb')
            loaded_res = pickle.load(filehandler)
            self.config = loaded_res['Run settings']
            self._results = loaded_res['Results']
            self._model = loaded_res['Expression']
            self._rdata_in = loaded_res['In data']
            self._rdata_out = loaded_res['Out data']
            print(str(filename), 'successfully loaded.')
            return
        except:
            raise Exception('File could not be loaded.')
