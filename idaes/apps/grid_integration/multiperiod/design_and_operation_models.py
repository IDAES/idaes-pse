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

from pyomo.environ import (
    Var,
    Binary,
    Expression,
    NonNegativeReals,
)

from idaes.core.base.process_base import declare_process_block_class
from idaes.models.unit_models import SkeletonUnitModelData
from pyomo.common.config import ConfigValue, In


@declare_process_block_class("DesignModel")
class DesignModelData(SkeletonUnitModelData):
    """
    Class for containing design model ...
    """

    CONFIG = SkeletonUnitModelData.CONFIG()
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
        "declare_cost_expr",
        ConfigValue(
            default=True,
            domain=In([True, False]),
            doc="Should capex and fom attributes be defined?",
        ),
    )

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()

        # Build the model
        self.config.model_func(self, **self.config.model_args)

        cost_expr_list = {
            "capex": "CAPEX of the unit/process",
            "fom": "Annual Fixed O&M cost for the unit/process",
        }

        if self.config.declare_cost_expr:
            for ce in cost_expr_list.keys():
                if not hasattr(self, ce):
                    # TODO: Use IDAES warning
                    print(
                        f"declare_cost_expr is set to True, but {ce} is not defined in {self.config.model_func.__name__}"
                    )
                    print(f"setting {ce} to zero in {self.name}")
                    setattr(self, ce, Expression(expr=0, doc=cost_expr_list[ce]))


@declare_process_block_class("OperationModel")
class OperationModelData(SkeletonUnitModelData):
    """
    Class for containing design model ...
    """

    CONFIG = SkeletonUnitModelData.CONFIG()
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
        "declare_cost_expr",
        ConfigValue(
            default=True,
            domain=In([True, False]),
            doc="Should capex and fom attributes be defined?",
        ),
    )
    CONFIG.declare(
        "declare_op_vars",
        ConfigValue(
            default=True,
            domain=In([True, False]),
            doc="Should op_mode, startup, shutdown vars be defined?",
        ),
    )

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()

        # Build the model
        if self.config.declare_op_vars:
            self.op_mode = Var(
                within=Binary,
                doc="Binary var: 1 if the unit is operating, 0 otherwise",
            )

            self.startup = Var(
                within=Binary,
                doc="Binary var: 1 if the startup is initiated, 0 otherwise",
            )

            self.shutdown = Var(
                within=Binary,
                doc="Binary var: 1 if the shutdown is initiated, 0 otherwise",
            )
            self.aux_startup = Var(
                within=NonNegativeReals,
                doc= "Auxiliary variable for the product of startup and capacity"
            )
            
            self.aux_shutdown = Var(
                within=NonNegativeReals,
                doc= "Auxiliary variable for the product of shutdown and capacity"
            )

            self.aux_op_mode= Var(
                within=NonNegativeReals,
                doc="Auxiliary variable for the product of op_mode and capacity"
            )
            

        self.config.model_func(self, **self.config.model_args)

        cost_expr_list = {
            "non_fuel_vom": "Variable O&M cost excluding fuel/power",
            "fuel_cost": "Cost of fuel/power required for the unit/process",
            "elec_revenue": "Electricity revenue generated by the unit/process",
            "carbon_price": "Carbon price associated with CO2 emissions",
        }

        if self.config.declare_cost_expr:
            for ce in cost_expr_list.keys():
                if not hasattr(self, ce):
                    # TODO: Use IDAES warning
                    print(
                        f"declare_cost_expr is set to True, but {ce} is not defined in {self.config.model_func.__name__}"
                    )
                    print(f"setting {ce} to zero in {self.name}")
                    setattr(self, ce, Expression(expr=0, doc=cost_expr_list[ce]))
        if self.config.capacity_var is not None:
            self._add_capacity_aux_vars()

    def _add_capacity_aux_vars(self):
        v1 = self.config.capacity_var.name.split(".")[-1] + "_op_mode"
        v2 = self.config.capacity_var.name.split(".")[-1] + "_startup"
        v3 = self.config.capacity_var.name.split(".")[-1] + "_shutdown"

        setattr(self,Var(),v1)
        setattr(self,Var(),v2)
        setattr(self,Var(),v3)