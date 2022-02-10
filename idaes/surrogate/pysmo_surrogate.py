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
from enum import Enum
import subprocess
from io import StringIO
import sys
import os
import numpy as np
import pandas as pd
import json
import jsonpickle

from pyomo.environ import Constraint, value, sin, cos, log, exp, Set, Reals
from pyomo.common.config import ConfigValue, In, Path, ListOf, Bool
from pyomo.common.tee import TeeStream
from pyomo.common.fileutils import Executable
from pyomo.common.tempfiles import TempfileManager
from pyomo.common.config import PositiveInt, PositiveFloat

from idaes.surrogate.base.surrogate_base import SurrogateTrainer, SurrogateBase
from idaes.core.util.exceptions import ConfigurationError
from idaes.surrogate.pysmo import polynomial_regression as pr, radial_basis_function as rbf, kriging as krg
import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)



class PysmoPolyTrainer(SurrogateTrainer):

    CONFIG = SurrogateTrainer.CONFIG()

    CONFIG.declare('maximum_polynomial_order', ConfigValue(
        default=None,
        domain=PositiveInt,
        description="Maximum order of univariate terms. Maximum value is 10."))

    CONFIG.declare('number_of_crossvalidations', ConfigValue(
        default=3,
        domain=PositiveInt,
        description="Number of crossvalidations."))

    CONFIG.declare('training_split', ConfigValue(
        default=0.8,
        domain=PositiveFloat,
        description="Training-testing data split for PySMO."))

    CONFIG.declare('solution_method', ConfigValue(
        default=None,
        domain=In(['pyomo', 'mle', 'bfgs']),
        description="Method for solving regression problem. Must be one of the options ['pyomo', 'mle', 'bfgs']. "))

    CONFIG.declare('multinomials', ConfigValue(
        default=False,
        domain=Bool,
        description="Option for bi-variate pairwise terms in final polynomial"))

    # I/O file options
    CONFIG.declare("fname", ConfigValue(
        default= 'solution.pickle',
        domain=str,
        description="File name for PySMO result - must be a .pickle file. If this "
        "option is not None, then working files will not be deleted."))

    CONFIG.declare("overwrite", ConfigValue(
        default=True,
        domain=Bool,
        description="Flag indicating whether existing files can be overwritten."))



    def __init__(self, **settings):
        super().__init__(**settings)

        self._results = {}


    def train_surrogate(self):
        self._results['model type'] = 'poly'
        self._results['No. outputs'] = len(self._output_labels)
        self._results['models'] = {}
        self._results['metrics'] = {}
        self._results['metrics']['RMSE'] = {}
        self._results['metrics']['R2'] = {}
        self._results['pysmo_results'] = {}

        print(self.config.solution_method)

        for i in range(len(self._output_labels)):
            # Create each dataframe
            pysmo_input = pd.concat([
                self._training_dataframe[self._input_labels], 
                self._training_dataframe[ [self._output_labels[i]] ]
                ], axis=1)
    
            # Train model
            model = pr.PolynomialRegression(pysmo_input, pysmo_input, 
                maximum_polynomial_order = self.config.maximum_polynomial_order,
                training_split = self.config.training_split,
                solution_method = self.config.solution_method,
                multinomials = self.config.multinomials,
                fname = self._output_labels[i] + '_' + self.config.fname,
                overwrite = self.config.overwrite,
                number_of_crossvalidations = self.config.number_of_crossvalidations
                )
            variable_headers = model.get_feature_vector()
            model.training()

            # Extract variables
            vars = []
            for k in variable_headers.keys():
                vars.append(variable_headers[k])

            # Document results
            self._results['pysmo_results'][self._output_labels[i]] = model
            self._results['models'][self._output_labels[i]] = str(model.generate_expression(vars))
            self._results['metrics']['RMSE'][self._output_labels[i]] = model.errors['MSE'] ** 0.5
            self._results['metrics']['R2'][self._output_labels[i]]  = model.errors['R2']
            
            print('\nModel for output', self._output_labels[i], 'trained successfully.\n')

        return self

    def display_pysmo_results(self):
        return self._results['models']



class PysmoRBFTrainer(SurrogateTrainer):

    CONFIG = SurrogateTrainer.CONFIG()

    CONFIG.declare('basis_function', ConfigValue(
        default=None,
        domain=In(['linear', 'cubic', 'gaussian', 'mq', 'imq', 'spline']),
        description="Basis function for RBF."))

    CONFIG.declare('solution_method', ConfigValue(
        default=None,
        domain=In(['pyomo', 'algebraic', 'bfgs']),
        description="Method for solving RBF problem. Must be an instance of 'SolutionMethod (Enum)' "))

    CONFIG.declare('regularization', ConfigValue(
        default=None,
        domain=Bool,
        description="Option for regularization - results in a regression rather than interpolation. "
        "Produces more generalizable models. Useful for noisy data."))

    # I/O file options
    CONFIG.declare("fname", ConfigValue(
        default= 'solution.pickle',
        domain=str,
        description="File name for PySMO result - must be a .pickle file. If this "
        "option is not None, then working files will not be deleted."))

    CONFIG.declare("overwrite", ConfigValue(
        default=True,
        domain=Bool,
        description="Flag indicating whether existing files can be overwritten."))



    def __init__(self, **settings):
        super().__init__(**settings)

        self._results = {}


    def train_surrogate(self):
        self._results['model type'] =  self.config.basis_function + ' rbf'
        self._results['No. outputs'] = len(self._output_labels)
        self._results['models'] = {}
        self._results['metrics'] = {}
        self._results['metrics']['RMSE'] = {}
        self._results['metrics']['R2'] = {}
        self._results['pysmo_results'] = {}

        for i in range(len(self._output_labels)):
            # Create each dataframe
            pysmo_input = pd.concat([
                self._training_dataframe[self._input_labels], 
                self._training_dataframe[ [self._output_labels[i]] ]
                ], axis=1)
    
            # Train model
            model = rbf.RadialBasisFunctions(pysmo_input, 
                basis_function = self.config.basis_function,
                solution_method = self.config.solution_method,
                regularization = self.config.regularization,
                fname = self._output_labels[i] + '_' + self.config.fname,
                overwrite = self.config.overwrite
                )
            variable_headers = model.get_feature_vector()
            model.training()

            # Extract variables
            vars = []
            for k in variable_headers.keys():
                vars.append(variable_headers[k])

            # Document results
            self._results['pysmo_results'][self._output_labels[i]] = model
            self._results['models'][self._output_labels[i]] = str(model.generate_expression(vars))
            self._results['metrics']['RMSE'][self._output_labels[i]] = model.rmse
            self._results['metrics']['R2'][self._output_labels[i]]  = model.R2
            
            print('\nModel for output', self._output_labels[i], 'trained successfully.\n')

        return self

    def display_pysmo_results(self):
        return self._results['models']



class PysmoKrigingTrainer(SurrogateTrainer):

    CONFIG = SurrogateTrainer.CONFIG()

    CONFIG.declare('numerical_gradients', ConfigValue(
        default=True,
        domain=Bool,
        description="Coice of whether numerical gradients are used in kriging model training." 
        "Determines choice of optimization algorithm: Basinhopping (False) or BFGS (True)." 
        "Using the numerical gradient option leads to quicker (but in complex cases possible sub-optimal) convergence"))

    CONFIG.declare('regularization', ConfigValue(
        default=True,
        domain=Bool,
        description="Option for regularization - results in a regression rather than interpolation. "
        "Produces more generalizable models. Useful for noisy data."))

    # I/O file options
    CONFIG.declare("fname", ConfigValue(
        default= 'solution.pickle',
        domain=str,
        description="File name for PySMO result - must be a .pickle file. If this "
        "option is not None, then working files will not be deleted."))

    CONFIG.declare("overwrite", ConfigValue(
        default=True,
        domain=Bool,
        description="Flag indicating whether existing files can be overwritten."))


    def __init__(self, **settings):
        super().__init__(**settings)

        self._results = {}


    def train_surrogate(self):
        self._results['model type'] = 'kriging'
        self._results['No. outputs'] = len(self._output_labels)
        self._results['models'] = {}
        self._results['metrics'] = {}
        self._results['metrics']['RMSE'] = {}
        self._results['metrics']['R2'] = {}
        self._results['pysmo_results'] = {}

        for i in range(len(self._output_labels)):
            # Create each dataframe
            pysmo_input = pd.concat([
                self._training_dataframe[self._input_labels], 
                self._training_dataframe[ [self._output_labels[i]] ]
                ], axis=1)
    
            # Train model
            model = krg.KrigingModel(pysmo_input, 
                numerical_gradients = self.config.numerical_gradients,
                regularization = False,#self.config.regularization,
                fname = self._output_labels[i] + '_' + self.config.fname,
                overwrite = self.config.overwrite
                )
            variable_headers = model.get_feature_vector()
            model.training()

            # Extract variables
            vars = []
            for k in variable_headers.keys():
                vars.append(variable_headers[k])

            # Document results
            self._results['pysmo_results'][self._output_labels[i]] = model
            self._results['models'][self._output_labels[i]] = str(model.generate_expression(vars))
            self._results['metrics']['RMSE'][self._output_labels[i]] = model.training_rmse
            self._results['metrics']['R2'][self._output_labels[i]]  = model.training_R2
            
            print('\nModel for output', self._output_labels[i], 'trained successfully.\n')

        return self

    def display_pysmo_results(self):
        return self._results['models']





class PysmoSurrogate(SurrogateBase):
    def __init__(self, surrogate_expressions, input_labels, output_labels, input_bounds=None):
        super().__init__(input_labels, output_labels, input_bounds)
        self._surrogate_expressions = surrogate_expressions

    def evaluate_surrogate(self, inputs):
        """
        Method to method to evaluate the ALAMO surrogate model at a set of user
        provided values.

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
        inputdata = inputs[self._input_labels].to_numpy()
        outputs = np.zeros(shape=(inputs.shape[0], len(self._output_labels)))

        for i in range(inputdata.shape[0]):
            row_data = inputdata[i, :].reshape(1, len(self._input_labels))
            for o in range(len(self._output_labels)):
                o_name = self._output_labels[o]
                outputs[i, o] = self._surrogate_expressions._results['pysmo_results'][self._output_labels[o]].predict_output(row_data)

        return pd.DataFrame(data=outputs, index=inputs.index, columns=self._output_labels)


    def populate_block(self, block, additional_options=None):
        """
        Method to populate a Pyomo Block with surrogate model constraints.

        Args:
            block: Pyomo Block component to be populated with constraints.
            additional_options: None
               No additional options are required for this surrogate object
        Returns:
            None
        """

        # TODO: do we need to add the index_set stuff back in?
        output_set = Set(initialize=self._output_labels, ordered=True)
        def pysmo_rule(b, o):
            in_vars = block.input_vars_as_dict()
            out_vars = block.output_vars_as_dict()
            return out_vars[o] == self._surrogate_expressions._results['pysmo_results'][o].generate_expression(list(in_vars.values()))

        block.pysmo_constraint = Constraint(output_set, rule=pysmo_rule)


    def save(self, strm):
        """
        Save an instance of this surrogate to the strm so the model can be used later.

        Args:
           strm: IO.TextIO
              This is the python stream like a file object or StringIO that will be used
              to serialize the surrogate object. This method writes a string
              of json data to the stream.
        """
        obj = jsonpickle.encode(self._surrogate_expressions)
        json.dump({"models": obj,
                   "input_labels": self._input_labels,
                   "output_labels": self._output_labels,
                   "input_bounds": self._input_bounds},
                  strm)
        # j = json.dumps({k: vc[k].tolist() if isinstance(vc[k], np.ndarray) else vc[k] for k in p2._surrogate_expressions._results['pysmo_results']['z1'].__dict__})
        # j = json.dumps({k: vc[k].tolist() if isinstance(vc[k], np.ndarray) else to_json(vc[k]) if isinstance(vc[k], (pyomo.core.base.param._ParamData, pyomo.core.base.param.Param)) else vc[k].to_json() if isinstance(vc[k], (pd.Series, pd.DataFrame)) else vc[k] for k in p2._surrogate_expressions._results['pysmo_results']['z1'].__dict__})
    @classmethod
    def load(cls, strm):
        """
        Create an instance of a surrogate from a stream.

        Args:
           strm: stream
              This is the python stream containing the data required to load the surrogate.
              This is often, but does not need to be a string of json data.

        Returns: an instance of the derived class or None if it failed to load
        """
        d = json.load(strm)

        surrogate_expressions = d["surrogate"]
        input_labels = d["input_labels"]
        output_labels = d["output_labels"]

        # # Need to convert list of bounds to tuples
        # input_bounds = {}
        # for k, v in d["input_bounds"].items():
        #     input_bounds[k] = tuple(v)

        # return AlamoSurrogate(surrogate_expressions=surrogate_expressions,
        #                       input_labels=input_labels,
        #                       output_labels=output_labels,
        #                       input_bounds=input_bounds)
