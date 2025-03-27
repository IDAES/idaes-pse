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

from pyomo.environ import Binary, Param, Var
from pyomo.common.config import Bool, ConfigDict, ConfigValue
from idaes.core.base.process_base import declare_process_block_class
from idaes.core.base.process_base import ProcessBlockData
import idaes.logger as idaeslog

_logger = idaeslog.getLogger(__name__)


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

    def build(self):
        super().build()

        self.install_unit = Var(
            within=Binary,
            doc="Binary: 1, if the unit is installed, 0 otherwise",
        )

        if self.config.model_func is None:
            # Function that builds the design model is not specified
            _logger.warning(
                "The function that builds the design model is not specified."
                "model_func must declare all the necessary design variables,"
                "relations among design variables, capital cost correlations,"
                "and fixed operating and maintenance cost correlations."
            )
            return

        # Call the function that builds the design model
        self.config.model_func(self, **self.config.model_args)

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
            doc="Boolean flag to determine if `op_mode`, `startup`, and `shutdown` vars should be defined",
        ),
    )
    CONFIG.declare(
        "declare_lmp_param",
        ConfigValue(
            default=True,
            domain=Bool,
            doc="Boolean flag to determine if LMP data should automatically be appended to the model",
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

        if self.config.declare_lmp_param:
            self.LMP = Param(
                initialize=1,
                mutable=True,
                doc="Time-varying locational marginal prices (LMPs) [in $/MWh]",
            )

        if self.config.model_func is None:
            _logger.warning(
                "The function that builds the operation model is not specified."
                "model_func must declare all the necessary operation variables,"
                "relations among operation variables, and variable"
                "operating and maintenance cost correlations."
            )
            return

        # Call the function that builds the operation model
        self.config.model_func(self, **self.config.model_args)
