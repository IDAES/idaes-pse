"""
Common Surrogate interface for IDAES.
"""
from pathlib import Path
from typing import Dict
import yaml
from pyomo.common.config import ConfigBlock
import os.path, pickle
# from mypy_extensions import TypedDict

class Metrics:
    """Names for known types of metrics.
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
    pass


# Single Surrogate Modeler
class Surrogate:
    CONFIG = ConfigBlock()

    # Common Declarations
    # CONFIG.declare('overwrite', ConfigValue(default=None, domain=bool))
    # CONFIG.declare('fname', ConfigValue(default=None, domain=str))

    def __init__(self, **settings):

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
        _b_built = False
        self.config = self.CONFIG(kwargs)

    # Build

    def build_model(self):
        self._b_built = True

        self.pkl_info = {'In data': self._rdata_in.tolist(),
                         'Out data': self._rdata_out.tolist()}
                         # 'Run settings': self.config}

        pass

    def update_model(self):
        self.build_model()
        pass

    def generate_expression(self, variable_list):
        pass

    # Get Results

    def get_model(self):  # Pyomo Expression
        return self._model

    def get_metrics(self):  # Metrics Object
        return self._metrics

    def get_results(self):  # Pyomo Expression, Metrics Object
        return self._results

    # Data Handling

    def get_regressed_data(self):
        return (self._rdata_in, self._rdata_out)

    def regressed_data(self, r_in, r_out):  # 2D Numparray
        self._rdata_in = r_in
        self._rdata_out = r_out

    def get_validation_data(self):
        return (self._vdata_in, self._vdata.out)

    def validation_data(self, v_in, v_out):  # 2D Numparray
        self._vdata_in = v_in
        self._vdata_out = v_out

    # Using regressed model
    def calculate_outputs(inputs):  # 2D Numparray, use pyomo expression
        outputs = self._model(inputs)
        return

    # Additional Metrics - MUST NOT OVERWRITE MODELER METRICS

    def get_trained_metrics():  # 2D Nparray, 2D Numparray
        if not _b_built:
            print("Warning: No surrogate model regressed.")
            return
        pass

    def get_validated_metrics(xval, zval):  # 2D, 2D Numparray, use pyomo expression
        if not _b_built:
            print("Warning: No surrogate model regressed.")
            return
        pass

    def save_results(self, filename, overwrite=False):
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
