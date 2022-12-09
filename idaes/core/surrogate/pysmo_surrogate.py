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


# stdlib
import io
import json
from json import JSONEncoder, JSONDecodeError
import logging
from typing import Dict, Union

# third-party
import numpy as np
import pandas as pd

# package
import pyomo.core as pc
from pyomo.environ import Constraint, sin, cos, log, exp, Set, Param
from pyomo.common.config import ConfigValue, In, Bool, PositiveInt, PositiveFloat
from idaes.core.surrogate.base.surrogate_base import SurrogateTrainer, SurrogateBase
from idaes.core.surrogate.pysmo import (
    polynomial_regression as pr,
    radial_basis_function as rbf,
    kriging as krg,
)
import idaes.logger as idaeslog
from idaes.core.util import to_json


__author__ = "Oluwamayowa Amusat"
# Logging
# -------
_log = idaeslog.getLogger(__name__)

# Global variables
# ----------------
GLOBAL_FUNCS = {"sin": sin, "cos": cos, "log": log, "exp": exp}


class PysmoSurrogateTrainingResult:
    """
    This is an internal helper class for the PysmoTrainer class.

    It stores the results of a trained PySMO model for a single output.

    The class attributes stored are the model metrics, the PySMO model
    and a string of the surrogate model result.

    """

    def __init__(self):
        self.metrics = {}  # Metrics for the fit
        self._model = None  # Surrogate model object (set/get as property)
        self.expression_str = (
            ""  # Pyomo expression, as a string, for additional variables
        )

    @property
    def model(self):
        return self._model

    @model.setter
    def model(self, value):
        """

        This is the internal helper function that sets the 'expression_str' attribute by
        calling Pysmo's ``generate_expression()`` function.

        It also sets the model attribute.

        """
        # variable_names = list(value.get_feature_vector().values())
        if hasattr(value, "regression_data_columns"):
            extra_terms_feature_vector = list(
                value.feature_list[i] for i in value.regression_data_columns
            )
        else:
            extra_terms_feature_vector = list(
                value.feature_list[i] for i in value.x_data_columns
            )

        # extra_terms_feature_vector = list(value.feature_list[i] for i in value.regression_data_columns)
        self.expression_str = str(
            value.generate_expression(list(extra_terms_feature_vector))
        )
        self._model = value


class PysmoTrainedSurrogate:
    """
    A surrogate that is trained for one or more outputs. Contains internal (helper) functions.

    This object is used for saving and loading the surrogate model along with its metrics and expression.

    """

    def __init__(self, model_type=""):
        """
        This function adds and displays the results for a trained PySMO surrogate.

        Args:
            model_type: string representing the model type ('poly', 'kriging' or 'rbf')

        Returns:
            self object with additional attributes representing the input models, output labels, number of outputs and model type.
        """
        self._data = {}
        self.model_type = model_type
        self.num_outputs = 0
        self.output_labels = []
        self.input_labels, self.input_bounds = None, None

    def add_result(self, output_name, result):
        """
        This function updates the ``self._data`` and ``self.output_labels`` attributes with the trained model object and output name respectively.

        """
        self._data[output_name] = result
        self.output_labels.append(output_name)
        # self.num_outputs += 1

    def get_result(self, output_name) -> PysmoSurrogateTrainingResult:
        """
        This function returns the ``PysmoSurrogateTrainingResult`` result object for a apecified output variable.

        Args:
            output_name: output variable label of interest (string)

        Returns:
            ``PysmoSurrogateTrainingResult`` result object for output variable.

        """
        return self._data[output_name]

    def display_pysmo_results(self):
        """
        This function returns a string representation of the surrogate expressions for all output variables
        """
        for k in self.output_labels:
            print(k, ":", self._data[k].expression_str)


# Surrogate trainer classes
# -------------------------


class PysmoTrainer(SurrogateTrainer):
    """Base class for Pysmo surrogate trainer classes."""

    # Initialize with configuration for base SurrogateTrainer
    CONFIG = SurrogateTrainer.CONFIG()

    # Subclasses must override this with a specific surrogate model type name
    model_type = "base"

    def __init__(self, **settings):
        super().__init__(**settings)
        self._trained = PysmoTrainedSurrogate(model_type=self.model_type)

    def train_surrogate(self) -> PysmoTrainedSurrogate:
        """
        General workflow method for training a PySMO surrogate.

        Takes the existing data set and executes the PySMO workflow to create
        a PysmoTrainedSurrogate object containing the trained model based on the current configuration arguments.

        The PySMOTrainedSurrogate object is the expected input for the PysmoSurrogate class.

        Args:
            None

        Returns:
            Python object : an instance of PysmoTrainedSurrogate representing the trained surrogate

        """
        self._trained.num_outputs = len(self._output_labels)
        self._trained.input_labels = self._input_labels
        if hasattr(self, "_input_bounds"):
            self._trained.input_bounds = self._input_bounds
        self._training_main_loop()
        return self._trained

    def _create_model(
        self, pysmo_input: pd.DataFrame, output_label: str
    ) -> Union[pr.PolynomialRegression, rbf.RadialBasisFunctions, krg.KrigingModel]:
        """Subclasses must override this and make it return a PySMO model."""
        raise NotImplementedError(
            "Sub-class fail to implement overload ``_create_model`` method."
        )

    def _get_metrics(self, model) -> Dict:
        """Subclasses should override this to return a dict of metrics for the model."""
        return {}

    def _training_main_loop(self):
        for i, output_label in enumerate(self._output_labels):
            # Create input dataframe
            pysmo_input = pd.concat(
                [
                    self._training_dataframe[self._input_labels],
                    self._training_dataframe[[self._output_labels[i]]],
                ],
                axis=1,
            )
            # Create and train model
            model = self._create_model(pysmo_input, output_label)
            model.training()
            # Store results
            result = PysmoSurrogateTrainingResult()
            result.model = model
            result.metrics = self._get_metrics(model)
            self._trained.add_result(output_label, result)
            # Log the status
            _log.info(f"Model for output {output_label} trained successfully")


class PysmoPolyTrainer(PysmoTrainer):
    """
    Standard SurrogateTrainer for PySMO's polynomial models.

    This defines a set of configuration options for PySMO, along with
    methods to train the polynomial model and return basic metrics (R2, RMSE).

    """

    model_type = "poly"

    CONFIG = PysmoTrainer.CONFIG()

    CONFIG.declare(
        "maximum_polynomial_order",
        ConfigValue(
            default=None,
            domain=PositiveInt,
            description="Maximum order of univariate terms. Maximum value is 10.",
        ),
    )

    CONFIG.declare(
        "number_of_crossvalidations",
        ConfigValue(
            default=3, domain=PositiveInt, description="Number of crossvalidations."
        ),
    )

    CONFIG.declare(
        "training_split",
        ConfigValue(
            default=0.8,
            domain=PositiveFloat,
            description="Training-testing data split for PySMO.",
        ),
    )

    CONFIG.declare(
        "solution_method",
        ConfigValue(
            default=None,
            domain=In(["pyomo", "mle", "bfgs"]),
            description="Method for solving regression problem. Must be one of the options ['pyomo', 'mle', 'bfgs']. ",
        ),
    )

    CONFIG.declare(
        "multinomials",
        ConfigValue(
            default=False,
            domain=Bool,
            description="Option for bi-variate pairwise terms in final polynomial",
        ),
    )

    CONFIG.declare(
        "extra_features",
        ConfigValue(
            default=None,
            domain=list,
            description="List of extra features to be considered for regression (if any), e.g. ['x1 / x2']. ",
        ),
    )

    def __init__(self, **settings):
        super().__init__(**settings)

    def _create_model(self, pysmo_input, output_label):
        model = pr.PolynomialRegression(
            pysmo_input,
            pysmo_input,
            maximum_polynomial_order=self.config.maximum_polynomial_order,
            training_split=self.config.training_split,
            solution_method=self.config.solution_method,
            multinomials=self.config.multinomials,
            number_of_crossvalidations=self.config.number_of_crossvalidations,
            overwrite=True,
        )
        variable_headers = model.get_feature_vector()
        if self.config.extra_features is not None:
            # create additional terms
            try:
                add_terms = self.config.extra_features
                for j in model.regression_data_columns:
                    add_terms = [
                        add_terms[k].replace(j, "variable_headers['" + str(j) + "']")
                        for k in range(0, len(add_terms))
                    ]
                model.set_additional_terms(
                    [
                        eval(m, GLOBAL_FUNCS, {"variable_headers": variable_headers})
                        for m in add_terms
                    ]
                )
            except:
                raise ValueError("Additional features could not be constructed.")
        return model

    def _get_metrics(self, model):
        return {"RMSE": model.errors["MSE"] ** 0.5, "R2": model.errors["R2"]}

    def get_confidence_intervals(
        self, model: PysmoTrainedSurrogate, confidence: float = 0.95
    ) -> Dict:
        """
        Compute confidence intervals for the regression patamaters.

        Args:
            model           : A PysmoTrainedSurrogate object
            confidence      : Required confidence interval level, default = 0.95 (95%)

        Returns:
            dict(<dict>)    : Dictionary object containing confidence intervals for all regressed parameters.

                              The dictionary keys are the output variables originally supplied during model training.

                              The dictionary values are dataframes containing four columns:

                                - Regression coeff.  : The regression coefficients for the trained model
                                - Std. errors        : The standard error on the estimated coefficient
                                - Conf. int. lower   : Lower confidence bounds for the estimated regression parameters
                                - Conf. int. upper   : Upper confidence bounds for the estimated regression parameters
        """
        confint_dict = {}
        for i in model.output_labels:
            confint_dict[i] = model._data[i].model.confint_regression(confidence)
        return confint_dict


class PysmoRBFTrainer(PysmoTrainer):
    """

    Standard SurrogateTrainer for PySMO's Radial Basis Function (RBF) models.

    This defines a set of configuration options for PySMO, along with
    methods to train the polynomial model and return basic metrics (R2, RMSE).

    """

    base_model_type = "rbf"
    model_type = "rbf"

    CONFIG = SurrogateTrainer.CONFIG()

    CONFIG.declare(
        "basis_function",
        ConfigValue(
            default=None,
            domain=In(["linear", "cubic", "gaussian", "mq", "imq", "spline"]),
            description="Basis function for RBF.",
        ),
    )

    CONFIG.declare(
        "solution_method",
        ConfigValue(
            default=None,
            domain=In(["pyomo", "algebraic", "bfgs"]),
            description="Method for solving RBF problem. Must be an instance of 'SolutionMethod (Enum)' ",
        ),
    )

    CONFIG.declare(
        "regularization",
        ConfigValue(
            default=None,
            domain=Bool,
            description="Option for regularization - results in a regression rather than interpolation. "
            "Produces more generalizable models. Useful for noisy data.",
        ),
    )

    def __init__(self, **settings):
        super().__init__(**settings)
        self.model_type = f"{self.config.basis_function} {self.base_model_type}"

    def _create_model(self, pysmo_input, output_label):
        model = rbf.RadialBasisFunctions(
            pysmo_input,
            basis_function=self.config.basis_function,
            solution_method=self.config.solution_method,
            regularization=self.config.regularization,
            overwrite=True,
        )
        variable_headers = model.get_feature_vector()
        return model

    def _get_metrics(self, model) -> Dict:
        return {"R2": model.R2, "RMSE": model.rmse}


class PysmoKrigingTrainer(PysmoTrainer):
    """

    Standard SurrogateTrainer for PySMO's kriging models.

    This defines a set of configuration options for PySMO, along with
    methods to train the polynomial model and return mbasic metrics (R2, RMSE).

    """

    model_type = "kriging"

    CONFIG = PysmoTrainer.CONFIG()

    CONFIG.declare(
        "numerical_gradients",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Choice of whether numerical gradients are used in Kriging model training."
            "Determines choice of optimization algorithm: Basinhopping (False) or BFGS (True)."
            "Using the numerical gradient option leads to quicker (but in complex cases possible sub-optimal)"
            " convergence",
        ),
    )

    CONFIG.declare(
        "regularization",
        ConfigValue(
            default=True,
            domain=Bool,
            description="Option for regularization - results in a regression rather than interpolation. "
            "Produces more generalizable models. Useful for noisy data.",
        ),
    )

    def __init__(self, **settings):
        super().__init__(**settings)

    def _create_model(self, pysmo_input, output_label):
        model = krg.KrigingModel(
            pysmo_input,
            numerical_gradients=self.config.numerical_gradients,
            regularization=self.config.regularization,
            overwrite=True,
        )
        variable_headers = model.get_feature_vector()
        return model

    def _get_metrics(self, model):
        return {"RMSE": model.training_rmse, "R2": model.training_R2}


class PysmoSurrogate(SurrogateBase):
    """PySMO surrogate model API."""

    def __init__(
        self,
        trained_surrogates: PysmoTrainedSurrogate,
        input_labels,
        output_labels,
        input_bounds=None,
    ):
        """Create PySMO surrogate model.
        Args:
            trained_surrogates: Results of training surrogates.
            input_labels:
            output_labels:
            input_bounds:
        """
        super().__init__(input_labels, output_labels, input_bounds)
        self._trained = trained_surrogates
        self._trained.input_labels, self._trained.input_bounds = (
            input_labels,
            input_bounds,
        )

    def evaluate_surrogate(self, inputs: pd.DataFrame) -> pd.DataFrame:
        """Evaluate the surrogate model at a set of user-provided values.

        Args:
            inputs: The dataframe of input values to be used in the evaluation.
                The dataframe needs to contain a column corresponding to each of the input labels.
                Additional columns are fine, but are not used.

        Returns:
            output: A dataframe of the the output values evaluated at the provided inputs.
                The index of the output dataframe should match the index of the provided inputs.
        """
        inputdata = inputs[self._input_labels].to_numpy()
        outputs = np.zeros(shape=(inputs.shape[0], len(self._output_labels)))

        for i in range(inputdata.shape[0]):
            row_data = inputdata[i, :].reshape(1, len(self._input_labels))
            for j, output_label in enumerate(self._output_labels):
                result = self._trained.get_result(output_label)
                outputs[i, j] = result.model.predict_output(row_data)

        return pd.DataFrame(
            data=outputs, index=inputs.index, columns=self._output_labels
        )

    def populate_block(self, block, additional_options=None):
        """Populate a Pyomo Block with surrogate model constraints.

        Args:
            block: Pyomo Block component to be populated with constraints.
            additional_options: None
                No additional options are required for this surrogate object

        Returns:
            None
        """

        output_set = Set(initialize=self._output_labels, ordered=True)

        def pysmo_rule(b, o):
            in_vars = block.input_vars_as_dict()
            out_vars = block.output_vars_as_dict()
            return out_vars[o] == self._trained.get_result(o).model.generate_expression(
                list(in_vars.values())
            )

        block.pysmo_constraint = Constraint(output_set, rule=pysmo_rule)

    def save(self, stream: io.TextIOBase):
        """Save this surrogate to the provided output stream so the model can be used later.

        Args:
           stream: Output stream for serialized surrogate object.

        Returns:
            None
        """
        # All custom serialization is done by the class provided to json.dump
        json.dump(self._trained, stream, cls=TrainedSurrogateEncoder)

    @classmethod
    def load(cls, stream):
        """Create an instance of a surrogate from a stream.

        Args:
            stream: <_io.StringIO>
                This is the python stream containing the data required to load the surrogate.
                This is often, but does not need to be a string of json data.

        Returns:
            An instance of the derived class or None if it failed to load
        """
        stream.seek(0)
        try:
            json_obj = json.load(
                stream, object_pairs_hook=TrainedSurrogateDecoder.decode_pairs
            )
            obj = PysmoSurrogate(
                trained_surrogates=json_obj[TSEBase.MODEL_KEY],
                input_labels=json_obj[TSEBase.INPUT_KEY],
                output_labels=json_obj[TSEBase.OUTPUT_KEY],
                input_bounds=json_obj[TSEBase.BOUNDS_KEY],
            )
        except JSONDecodeError as err:
            obj = None
            _log.error(f"Error decoding surrogate model (stream={stream}): {err}")
        return obj


# Serialization and deserialization classes
# ------------------------------------------


class TSEBase:
    # top-level keys shared by encoder/decoder classes
    MODEL_KEY = "model_encoding"
    INPUT_KEY = "input_labels"
    OUTPUT_KEY = "output_labels"
    BOUNDS_KEY = "input_bounds"
    TYPE_KEY = "surrogate_type"
    MODEL_ATTR_KEY = "attr"
    MODEL_MAP_KEY = "map"
    MODEL_METRICS_KEY = "errors"


class TrainedSurrogateEncoder(JSONEncoder, TSEBase):
    # Attributes to encode
    attrs = {
        "final_polynomial_order",
        "multinomials",
        "optimal_weights_array",
        "extra_terms_feature_vector",
        "additional_term_expressions",
        "regression_data_columns",
        "errors",
        "centres",
        "x_data_columns",
        "x_data_min",
        "x_data_max",
        "basis_function",
        "self.sigma",
        "y_data_min",
        "y_data_max",
        "weights",
        "sigma",
        "regularization_parameter",
        "R2",
        "rmse",
        "optimal_weights",
        "optimal_p",
        "optimal_mean",
        "optimal_variance",
        "regularization_parameter",
        "optimal_covariance_matrix",
        "covariance_matrix_inverse",
        "optimal_y_mu",
        "training_R2",
        "training_rmse",
        "x_data",
        "x_data_scaled",
    }

    def default(self, obj):
        """
        This function conforms to the interface expected by the json package.
        """
        if isinstance(obj, PysmoTrainedSurrogate):
            return self._encode_surrogate(obj)
        return super().default(obj)

    def _encode_surrogate(self, obj):
        """Encode the trained surrogate for all outputs."""
        models_enc = {
            o: self._encode_model(obj.get_result(o).model) for o in obj.output_labels
        }
        return {
            self.MODEL_KEY: models_enc,
            self.INPUT_KEY: obj.input_labels,
            self.OUTPUT_KEY: obj.output_labels,
            self.BOUNDS_KEY: obj.input_bounds,
            self.TYPE_KEY: obj.model_type,
        }

    def _encode_model(self, model):
        """Encode the surrogate model for a single output."""
        attr_values, attr_types = {}, {}
        # Walk through all attributes in the model that are also in self.attrs
        for name, val in (
            (v, getattr(model, v)) for v in vars(model) if v in self.attrs
        ):
            e = self._encode_attr(val)
            if e is not None:
                attr_values[name], attr_types[name] = self._encode_attr(val)
        return {self.MODEL_ATTR_KEY: attr_values, self.MODEL_MAP_KEY: attr_types}

    @staticmethod
    def _encode_attr(value):
        def is_pyomo_value(v):
            return isinstance(
                v,
                (
                    pc.base.param._ParamData,
                    pc.base.param.Param,
                    pc.expr.numeric_expr.NPV_ProductExpression,
                    pc.expr.numeric_expr.NPV_DivisionExpression,
                ),
            )

        if isinstance(value, np.ndarray):
            return value.tolist(), "numpy"
        elif isinstance(value, (pd.Series, pd.DataFrame)):
            return value.to_json(orient="index"), "pandas"
        elif is_pyomo_value(value):
            return to_json(value, return_dict=True), "pyomo"
        elif isinstance(value, list):
            if len(value) > 0 and is_pyomo_value(value[0]):
                return [str(k) for k in value], "other"  # skip
            return [str(v) for v in value], "list"
        else:
            return value, "str"


class TrainedSurrogateDecoder(TSEBase):
    """Decode a surrogate that was encoded into JSON by :class:`TrainedSurrogateEncoder`."""

    # TODO: Need to test/debug all these methods

    @classmethod
    def decode_pairs(cls, pairs):
        d = {}
        model_json, model_type = None, None
        for key, value in pairs:
            if _log.isEnabledFor(logging.DEBUG):
                _log.debug(f"JSON decode pairs. key={key}, value={value}")
            if key == cls.MODEL_KEY:
                model_json = value
            elif key == cls.TYPE_KEY:
                model_type = value
            elif key == cls.BOUNDS_KEY:
                # convert list values into tuples
                if value is None:
                    d[key] = None
                else:
                    d[key] = {k: tuple(v) for k, v in value.items() if value != None}
            else:
                d[key] = value
        if model_json is None:
            # If this is not present, assume that we are not at the outer level, so stop here
            return d
        if model_type is None:
            raise JSONDecodeError(
                f"Missing key '{cls.TYPE_KEY}' in JSON object, so type of model is unknown",
                doc="",
                pos=0,
            )
        # choose decoder for model type
        base_type = model_type.split(" ")[-1]  # for "<foo> rbf" types, picks out 'rbf'
        try:
            model_decoder = getattr(cls, f"_decode_{base_type}_model")
        except AttributeError:
            raise JSONDecodeError(
                f"No decoder found for model type '{base_type}'",
                doc="",
                pos=0,
            )
        # decode model
        _log.info(f"Decode surrogate. type={model_type}")
        trained = PysmoTrainedSurrogate(model_type)
        for output_label, model_data in model_json.items():
            surr_mod = model_decoder(
                model_type,
                model_data[cls.MODEL_ATTR_KEY],
                model_data[cls.MODEL_MAP_KEY],
            )
            result = PysmoSurrogateTrainingResult()
            result.model = surr_mod
            # result.metrics = model_data[cls.MODEL_ATTR_KEY][cls.MODEL_METRICS_KEY]
            trained.add_result(output_label, result)
        d[cls.MODEL_KEY] = trained
        # return decoded dict
        return d

    @classmethod
    def _decode_poly_model(cls, type_, attr, mapping):
        # Construct model with minimal arguments necessary
        columns = ["x", "y"]  # attr["regression_data_columns"]
        orig_data = pd.DataFrame(
            {c: list(range(10)) for c in columns}
        )  # cls.pd_decode(attr["original_data"])
        regr_data = pd.DataFrame(
            {c: list(range(10)) for c in columns}
        )  # cls.pd_decode(attr["regression_data"])
        max_order = 1  # int(attr["final_polynomial_order"])
        _log.debug(f"Reconstructing PolynomialRegression surrogate.")
        model = pr.PolynomialRegression(orig_data, regr_data, max_order)
        # Set model attributes from saved attributes
        _log.debug("Setting attributes on constructed surrogate model")
        for k, v in attr.items():
            if k not in [
                "feature_list",
                "extra_terms_feature_vector",
                "additional_term_expressions",
            ]:
                type_ = mapping[k]
                decoder = decoders.get(type_, null_decode)
                if _log.isEnabledFor(logging.DEBUG):
                    _log.debug(f"Decode attribute. decoder={decoder}, attr={k}:{type_}")
                decoded_value = decoder(v)
                setattr(model, k, decoded_value)
        p = Param(model.regression_data_columns, mutable=True, initialize=0)
        p.index_set().construct()
        p.construct()
        model.feature_list = p
        model.extra_terms_feature_vector = list(
            model.feature_list[i] for i in model.regression_data_columns
        )
        # Re-create function objects from additional terms
        list_terms = cls._poly_decode_vars(attr["additional_term_expressions"], p)
        model.additional_term_expressions = [
            eval(m, GLOBAL_FUNCS, {"p": p}) for m in list_terms
        ]

        # Remove all attributes created by default
        for jk in list(vars(model).keys()):
            if jk not in attr.keys() and jk != "feature_list":
                delattr(model, jk)

        # Done: return new model
        return model

    @staticmethod
    def _poly_decode_vars(xyz, p):
        """Decode the variables in a polynomial expression."""
        list_idx_vars = [p._data[i].local_name for i in p._data.keys()]
        list_vars = ['p["' + str(i) + '"]' for i in p.keys()]
        pyomo_vars_expr = xyz
        for i in range(0, len(list_idx_vars)):
            pyomo_vars_expr = [
                var_name.replace(list_idx_vars[i], list_vars[i])
                for var_name in pyomo_vars_expr
            ]
        #  return [eval(r, {}, {"p":p}) for r in pyomo_vars_expr]
        return pyomo_vars_expr

    @classmethod
    def _decode_rbf_model(cls, type_, attr, mapping):
        # Construct model with minimal arguments necessary
        columns = ["x", "y"]
        XY_data = pd.DataFrame({c: list(range(10)) for c in columns})
        model = rbf.RadialBasisFunctions(XY_data)  # XXX
        # Set model attributes from saved attributes
        for k, v in attr.items():
            setattr(model, k, decoders.get(mapping[k], null_decode)(v))
        p = Param(model.x_data_columns, mutable=True, initialize=0)
        p.index_set().construct()
        p.construct()
        model.feature_list = p

        # Remove all attributes created by default
        for jk in list(vars(model).keys()):
            if jk not in attr.keys() and jk != "feature_list":
                delattr(model, jk)
        # Done: return new model
        return model

    @classmethod
    def _decode_kriging_model(cls, type_, attr, mapping):
        # Construct model with minimal arguments necessary
        columns = ["x", "y"]
        XY_data = pd.DataFrame({c: list(range(10)) for c in columns})
        model = krg.KrigingModel(XY_data)
        # Set model attributes from saved attributes
        for k, v in attr.items():
            setattr(model, k, decoders.get(mapping[k], null_decode)(v))
        p = Param(model.x_data_columns, mutable=True, initialize=0)
        p.index_set().construct()
        p.construct()
        model.feature_list = p

        # Remove all attributes created by default
        for jk in list(vars(model).keys()):
            if jk not in attr.keys() and jk != "feature_list":
                delattr(model, jk)

        # Done: return new model
        return model


# Helper methods for decoding Pandas and Numpy


def pd_decode(v):
    return pd.read_json(v, orient="index")


def numpy_decode(v):
    return np.array(v)


def null_decode(v):
    return v


decoders = {"numpy": numpy_decode, "pandas": pd_decode}
