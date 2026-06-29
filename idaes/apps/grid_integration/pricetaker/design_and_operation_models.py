#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2026 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################

from pyomo.environ import Binary, Constraint, Expression, NonNegativeReals, Param, Var
from pyomo.common.config import (
    Bool,
    ConfigDict,
    ConfigValue,
    PositiveInt,
    NonNegativeFloat,
)
from idaes.core.base.process_base import declare_process_block_class
from idaes.core.base.process_base import ProcessBlockData
from idaes.core.util.config import ConfigurationError, is_in_range
import idaes.logger as idaeslog

_logger = idaeslog.getLogger(__name__)


def _format_data(coeffs):
    """
    Helper function to correctly format the surrogate model coefficients
    """
    if isinstance(coeffs, (int, float, Param, Var)):
        # Here, correlation is assumed to be y = a_1 * design_var
        return [0, coeffs]

    if isinstance(coeffs, (list, tuple)):
        # Here, correlation is assume to be
        # y = a_0 + a_1 * design_var + a_2 * design_var**2 + ...
        return coeffs

    raise ConfigurationError(
        f"Unrecognized data structure {type(coeffs)} for auxiliary variable coefficients."
    )


def is_valid_variable_design_data(data: dict):
    """Validates the arguments received for the variable design case"""
    # Ensure that design_var and design_var_bounds are present
    if "design_var" not in data:
        raise ConfigurationError("design_var is not specified")

    if "design_var_bounds" not in data:
        raise ConfigurationError("design_var_bounds is not specified")

    new_data = {
        "design_var": data.pop("design_var"),
        "design_var_bounds": data.pop("design_var_bounds"),
        "auxiliary_vars": {},
    }

    # The rest of the variables are assumed to be correlations
    for var_name, coeff in data.items():
        new_data["auxiliary_vars"][var_name] = _format_data(coeff)

    return new_data


def is_valid_polynomial_surrogate_data(data: dict):
    """Validates the arguments received for the variable design case"""
    # Ensure that operation_var is present present
    if "operation_var" not in data:
        raise ConfigurationError("operation_var is not specified")

    new_data = {
        "operation_var": data.pop("operation_var"),
        "auxiliary_vars": {},
    }

    # The rest of the variables are assumed to be correlations
    for var_name, coeff in data.items():
        new_data["auxiliary_vars"][var_name] = _format_data(coeff)

    return new_data


def is_valid_startup_types(data):
    """Validate if the startup_types received is valid"""

    if not isinstance(data, dict):
        raise TypeError("Data must be a dictionary.")

    if len(data) == 0:
        # set to None instead of an empty dict.
        raise ConfigurationError("Received an empty dictionary for startup types")

    if len(data) == 1:
        raise ConfigurationError(
            "At least two startup types must be defined for the unit/process. \n"
            'if startup_types = {"hot": 4, "warm": 8, "cold": 12}, \n'
            "then when downtime is between 4 to 8 -> hot startup; \n"
            "when downtime is between 8 to 12 -> warm startup; \n"
            "when downtime is greater than 12 -> cold startup. \n"
        )

    # check the key names and make sure space is handled correctly.
    # use a method to check if a string is a valid variable name. string.isidentifier() works for now.

    for key, value in data.items():
        if not isinstance(key, str):
            raise TypeError("key must be a valid string.")

        if not key.isidentifier():
            raise ConfigurationError(
                f"Key '{key}' is not a valid Python variable name. "
                "Keys must be valid identifiers."
            )

        if not isinstance(value, int):
            raise TypeError("value must be an int")

    # Values must be unique, as they correspond to different startup types
    if len(data.values()) > len(set(data.values())):
        raise ConfigurationError(
            "Startup time for two or more startup types is the same."
        )

    # Return a dictionary after sorting based on values
    return dict(sorted(data.items(), key=lambda item: item[1]))


# pylint: disable = attribute-defined-outside-init, too-many-ancestors
# pylint: disable = invalid-name, logging-fstring-interpolation
@declare_process_block_class("DesignModel")
class DesignModelData(ProcessBlockData):
    """
    Builds the 'design model' for a unit/process.

    Args:
        model_func: Function that builds the design model
        model_args: Dictionary containing the arguments needed for model_func

    The class defines `install_unit` binary variable that takes the value 1
    if the unit is built/installed, and 0 otherwise.

    Function model_func must declare all the necessary design variables,
    relations among design variables, capital cost correlations, and fixed O&M
    cost correlations. The function must also define attributes `capex` for
    capital cost, and `fom` for fixed O&M cost. If not defined, these attributes
    will be set to zero.

    Example Usage:

    .. code-block:: python

        def my_design_model(m, p_min, p_max, cost):
            m.power = Var()
            m.min_capacity = Constraint(
                expr=p_min * m.install_unit <= m.power
            )
            m.max_capacity = Constraint(
                expr=m.power <= p_max * m.install_unit
            )

            # capex and fom must either be a constant, or Var, or Expression
            m.capex = Expression(expr=cost["capex"] * m.power)
            m.fom = Expression(expr=cost["fom"] * m.power)

        m = ConcreteModel()
        m.unit_1 = DesignModel(
            model_func=my_design_model,
            model_args={
                "p_min": 150, "p_max": 600, "cost": {"capex": 10, "fom": 1},
            },
        )
    """

    CONFIG = ConfigDict()
    CONFIG.declare(
        "model_func",
        ConfigValue(
            doc="Function that builds the design model",
        ),
    )
    CONFIG.declare(
        "model_args",
        ConfigValue(
            default={},
            doc="Dictionary containing arguments needed for model_func",
        ),
    )
    CONFIG.declare(
        "fixed_design_data",
        ConfigValue(
            domain=dict,
            doc="Dictionary containing parameters associated with unit/process design",
        ),
    )
    CONFIG.declare(
        "variable_design_data",
        ConfigValue(
            domain=is_valid_variable_design_data,
            doc="Dictionary containing variables associated with unit/process design",
        ),
    )

    def build(self):
        super().build()

        self.install_unit = Var(
            within=Binary,
            doc="Binary: 1, if the unit is installed, 0 otherwise",
        )

        if self.config.fixed_design_data is not None:
            # The design is fixed, so all the desired quantities are defined as parameters
            for param_name, param_value in self.config.fixed_design_data.items():
                setattr(self, param_name, Param(initialize=param_value, mutable=True))

        elif self.config.variable_design_data is not None:
            # Design is a decision variable, so define desired variables and expressions
            self._build_variable_design_model()

        elif self.config.model_func is not None:
            # User has a custom design model. Call the function that builds the design model
            self.config.model_func(self, **self.config.model_args)

        else:
            # Function that builds the design model is not specified
            _logger.warning(
                "The function that builds the design model is not specified."
                "model_func must declare all the necessary design variables,"
                "relations among design variables, capital cost correlations,"
                "and fixed operating and maintenance cost correlations."
            )
            return

        # Check if capital and fixed O&M costs are defined
        if not hasattr(self, "capex"):
            _logger.warning(
                "'capex' attribute is not set for the design model "
                f"{self.name}. Setting the capital cost of the unit to zero."
            )
            self.capex = 0

        if not hasattr(self, "fom"):
            _logger.warning(
                "'fom' attribute is not set for the design model "
                f"{self.name}. Setting the fixed O&M cost of the unit to zero."
            )
            self.fom = 0

    def _build_variable_design_model(self):
        data = self.config.variable_design_data

        # Define the design variable
        setattr(
            self,
            data["design_var"],
            Var(
                doc="Design variable",
                within=NonNegativeReals,
                bounds=(0, data["design_var_bounds"][1]),
            ),
        )

        # Add bound constraints on the design variable
        design_var = getattr(self, data["design_var"])
        self.design_lb_constraint = Constraint(
            expr=self.install_unit * data["design_var_bounds"][0] <= design_var,
            doc="Ensures that the design is above the lb if the unit is built",
        )
        self.design_ub_constraint = Constraint(
            expr=design_var <= self.install_unit * data["design_var_bounds"][1],
            doc="Ensures that the design is less than ub if the unit is built",
        )

        # Build a polynomial correlation for all auxiliary variables
        for var_name, coeff in data["auxiliary_vars"].items():
            expr = coeff[0] * self.install_unit + sum(
                coeff[i] * design_var**i for i in range(1, len(coeff))
            )
            setattr(self, var_name, Expression(expr=expr))


@declare_process_block_class("OperationModel")
class OperationModelData(ProcessBlockData):
    """
    Builds the 'operation model' for a unit/process.

    Args:
        model_func: Function that builds the operation model
        model_args: Dictionary containing the arguments needed for model_func

    The class defines `op_mode`, `startup`, and `shutdown` binary variables
    to track the operation, startup, and shutdown of the unit/process.

    Function model_func must declare all the necessary operation variables,
    relations among operation variables, and variable O&M cost correlations.

    Example Usage:

    .. code-block:: python

        def my_operation_model(m, design_blk):
            m.power = Var()
            m.fuel_flow = Var()
            ...

        m = ConcreteModel()
        m.unit_1 = DesignModel(
            model_func=my_design_model,
            model_args={
                "p_min": 150, "p_max": 600, "cost": {"capex": 10, "fom": 1},
            },
        )
        m.op_unit_1 = OperationModel(
            model_func=my_operation_model, model_args={"design_blk": m.unit_1},
        )
    """

    CONFIG = ConfigDict()
    CONFIG.declare(
        "model_func",
        ConfigValue(
            doc="Function that builds the design model",
        ),
    )
    CONFIG.declare(
        "model_args",
        ConfigValue(
            default={},
            doc="Dictionary containing arguments needed for model_func",
        ),
    )
    CONFIG.declare(
        "declare_op_vars",
        ConfigValue(
            default=True,
            domain=Bool,
            doc=(
                "Boolean flag to determine if `op_mode`, `startup`, "
                "and `shutdown` vars should be defined"
            ),
        ),
    )
    CONFIG.declare(
        "declare_lmp_param",
        ConfigValue(
            default=True,
            domain=Bool,
            doc=(
                "Boolean flag to determine if LMP data should "
                "automatically be appended to the model"
            ),
        ),
    )
    CONFIG.declare(
        "polynomial_surrogate_data",
        ConfigValue(
            domain=is_valid_polynomial_surrogate_data,
            doc="Dictionary containing polynomial surrogate data",
        ),
    )

    CONFIG.declare(
        "design_block_name",
        ConfigValue(
            domain=str,
            doc="Object containing design model",
        ),
    )
    CONFIG.declare(
        "commodity",
        ConfigValue(
            domain=str,
            doc="Variable on which capacity and ramping limits are needed",
        ),
    )
    CONFIG.declare(
        "capacity",
        ConfigValue(doc="Maximum capacity of the commodity"),
    )
    CONFIG.declare(
        "op_range_lb",
        ConfigValue(
            domain=is_in_range(0, 1),
            doc="Minimum stable operation of the commodity as a fraction of capacity",
        ),
    )
    CONFIG.declare(
        "minimum_up_time",
        ConfigValue(
            domain=PositiveInt,
            doc=(
                "Total uptime (must be >= 1), e.g., 4 time periods. "
                "Uptime must include the minimum uptime and the time "
                "required for shutdown."
            ),
        ),
    )
    CONFIG.declare(
        "minimum_down_time",
        ConfigValue(
            domain=PositiveInt,
            doc=(
                "Total downtime (must be >= 1), e.g., 4 time periods. "
                "Downtime must include the minimum downtime and the time "
                "required for startup"
            ),
        ),
    )
    CONFIG.declare(
        "startup_rate",
        ConfigValue(
            domain=is_in_range(0, 1),
            doc="Startup rate as a fraction of the maximum capacity",
        ),
    )
    CONFIG.declare(
        "shutdown_rate",
        ConfigValue(
            domain=is_in_range(0, 1),
            doc="Shutdown rate as a fraction of the maximum capacity",
        ),
    )
    CONFIG.declare(
        "rampup_rate",
        ConfigValue(
            domain=is_in_range(0, 1),
            doc="Ramp up rate as a fraction of the maximum capacity",
        ),
    )
    CONFIG.declare(
        "rampdown_rate",
        ConfigValue(
            domain=is_in_range(0, 1),
            doc="Ramp down rate as a fraction of the maximum capacity",
        ),
    )

    CONFIG.declare(
        "startup_types",
        ConfigValue(
            domain=is_valid_startup_types,
            doc="Dictionary of startup types and transition times for the unit/process",
        ),
    )

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()

        # Build the model
        if self.config.declare_op_vars:
            self.op_mode = Var(
                within=Binary,
                doc="Binary: 1 if the unit is operating, 0 otherwise",
            )

            self.startup = Var(
                within=Binary,
                doc="Binary: 1 if the startup is initiated, 0 otherwise",
            )

            self.shutdown = Var(
                within=Binary,
                doc="Binary: 1 if the shutdown is initiated, 0 otherwise",
            )

        if self.config.startup_types is not None:
            self.startup_type_vars = Var(
                list(self.config.startup_types.keys()),
                within=Binary,
                doc="Binary: 1 if the startup type is active, 0 otherwise",
            )

        if self.config.declare_lmp_param:
            self.LMP = Param(
                initialize=1,
                mutable=True,
                doc="Time-varying locational marginal prices (LMPs) [in $/MWh]",
            )

        if self.config.polynomial_surrogate_data is not None:
            # Build polynomial-type operation model
            setattr(self, self.config.polynomial_surrogate_data["operation_var"], Var())
            self.build_polynomial_surrogates(
                surrogates=self.config.polynomial_surrogate_data["auxiliary_vars"],
                op_var=getattr(
                    self, self.config.polynomial_surrogate_data["operation_var"]
                ),
                declare_variables=True,
            )

        elif self.config.model_func is not None:
            # User has a custom operation model
            self.config.model_func(self, **self.config.model_args)

        else:
            _logger.warning(
                "The function that builds the operation model is not specified."
                "model_func must declare all the necessary operation variables,"
                "relations among operation variables, and variable"
                "operating and maintenance cost correlations."
            )

    def build_polynomial_surrogates(
        self, surrogates: dict, op_var: Var, declare_variables: bool = True
    ):
        """Builds polynomial-type surrgogate models"""
        surrogate_expressions = {
            expr_name: coeffs[0] * self.op_mode
            + sum(coeffs[i] * op_var**i for i in range(1, len(coeffs)))
            for expr_name, coeffs in surrogates.items()
        }
        self.build_expressions(surrogate_expressions, declare_variables)

    def build_expressions(self, expressions: dict, declare_variables: bool = False):
        """Declares user-defined expressions"""
        for expr_name, expr in expressions.items():
            if declare_variables:
                # Declare an auxiliary variable for each expression
                setattr(self, expr_name, Var())
                setattr(
                    self,
                    "compute_" + expr_name,
                    Constraint(expr=getattr(self, expr_name) == expr),
                )

            else:
                # Declare the expression as a Pyomo Expression
                setattr(self, expr_name, Expression(expr=expr))


def _is_valid_data_type_for_storage_model(obj):
    """Domain validator for config arguments of the storage model"""
    if isinstance(obj, (int, float, Param, Var, Expression)):
        # If the input is an unnamed expression, this check fails.
        # It also fails if the input is an indexed Var, Param, Expression
        return obj

    # Doing the following for unnamed expressions and indexed components
    try:
        # Assume that the input argument is a Pyomo component
        if any(
            (obj.is_expression_type(), obj.is_variable_type(), obj.is_parameter_type())
        ):
            return obj

        raise ConfigurationError(
            "Received an invalid Pyomo component as an argument for StorageModel"
        )

    # pylint: disable = raise-missing-from
    except AttributeError:
        raise ConfigurationError("Received an unsupported data type for StorageModel")


@declare_process_block_class("StorageModel")
class StorageModelData(ProcessBlockData):
    """
    Builds the 'storage model' for a unit/process.
    """

    CONFIG = ConfigDict()
    CONFIG.declare(
        "time_interval",
        ConfigValue(
            default=1,
            domain=NonNegativeFloat,
            doc="Length of each time interval",
        ),
    )
    CONFIG.declare(
        "charge_efficiency",
        ConfigValue(
            default=1,
            domain=is_in_range(0, 1),
            doc="Efficiency associated with charging",
        ),
    )
    CONFIG.declare(
        "discharge_efficiency",
        ConfigValue(
            default=1,
            domain=is_in_range(0, 1),
            doc="Efficiency associated with discharging",
        ),
    )
    CONFIG.declare(
        "min_holdup",
        ConfigValue(
            domain=_is_valid_data_type_for_storage_model,
            doc="Minimum holdup required",
        ),
    )
    CONFIG.declare(
        "max_holdup",
        ConfigValue(
            domain=_is_valid_data_type_for_storage_model,
            doc="Maximum holdup feasible",
        ),
    )
    CONFIG.declare(
        "max_charge_rate",
        ConfigValue(
            domain=_is_valid_data_type_for_storage_model,
            doc="Maximum charge rate allowed",
        ),
    )
    CONFIG.declare(
        "max_discharge_rate",
        ConfigValue(
            domain=_is_valid_data_type_for_storage_model,
            doc="Maximum discharge rate allowed",
        ),
    )

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()

        self.initial_holdup = Var(
            within=NonNegativeReals,
            doc="Holdup/charge at the beginning of the time interval",
        )
        self.final_holdup = Var(
            within=NonNegativeReals, doc="Holdup/charge at the end of the time interval"
        )
        self.charge_rate = Var(within=NonNegativeReals, doc="Charge rate")
        self.discharge_rate = Var(within=NonNegativeReals, doc="Discharge rate")

        def _add_upper_bound(var_name, bound):
            if isinstance(bound, (int, float)) or bound is None:
                # Set it as a bound if the input is either an int or a float.
                # But if a Pyomo component (Param, Var, Expression) is provided, then
                # set it as a constraint.
                getattr(self, var_name).setub(bound)

            else:
                setattr(
                    self,
                    var_name + "_ub_con",
                    Constraint(
                        expr=getattr(self, var_name) <= bound,
                        doc=f"Constrains the maximum value of {var_name}",
                    ),
                )

        _add_upper_bound("charge_rate", self.config.max_charge_rate)
        _add_upper_bound("discharge_rate", self.config.max_discharge_rate)
        _add_upper_bound("final_holdup", self.config.max_holdup)
        _add_upper_bound("initial_holdup", self.config.max_holdup)

        # pylint: disable = no-member
        if (
            isinstance(self.config.min_holdup, (int, float))
            or self.config.min_holdup is None
        ):
            self.final_holdup.setlb(self.config.min_holdup)
            self.initial_holdup.setlb(self.config.min_holdup)

        else:
            # Set a lower bound on holdup
            self.final_holdup_lb_con = Constraint(
                expr=self.final_holdup >= self.config.min_holdup,
                doc="Constrains the minimum value of final_holdup",
            )
            self.initial_holdup_lb_con = Constraint(
                expr=self.initial_holdup >= self.config.min_holdup,
                doc="Constrains the minimum value of initial_holdup",
            )

        # Mass balance/charge balance/tracking holdup
        self.track_holdup_constraint = Constraint(
            expr=(
                self.final_holdup - self.initial_holdup
                == (
                    self.config.charge_efficiency * self.charge_rate
                    - (self.discharge_rate / self.config.discharge_efficiency)
                )
                * self.config.time_interval
            ),
            doc="Models variation in holdup level over the time interval",
        )
