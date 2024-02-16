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
    Param,
    Binary,
    Expression,
    NonNegativeReals,
    Constraint,
)
from functools import reduce
from idaes.core.base.process_base import declare_process_block_class
from idaes.models.unit_models import SkeletonUnitModelData
from pyomo.common.config import ConfigValue, In

def deepgetattr(obj, attr):
    return reduce(getattr, attr.split('.'), obj)
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

    CONFIG.declare(
        "declare_build_vars",
        ConfigValue(
            default=True,
            domain=In([True, False]),
            doc="Should build_unit be declared/defined?",
        ),
    )

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()
        
        
        if self.config.declare_build_vars:
            self.build_unit = Var(
                within=Binary,
                doc="Binary var: 1 if we choose to build the unit, 0 otherwise",
            )
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
        "design_blk",
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
    CONFIG.declare(
        "append_lmp_data",
        ConfigValue(
            default=True,
            domain=In([True, False]),
            doc="Should LMP data automatically be appended to the model?"
        )
    )

    CONFIG.declare(
        "capacity_var",
        ConfigValue(
            default=None,
            doc="capacity variable"

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
        
        if self.config.append_lmp_data:
            self.LMP = Param(
                initialize=1, 
                mutable=True,
                doc="Parameter: Will be updated to LMP value at given time"
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
        if isinstance(self.config.capacity_var, str):
            self._add_capacity_aux_vars()

    def _add_capacity_aux_vars(self):
        aux_op_mode = self.config.capacity_var.split(".")[-1] + "_op_mode"
        aux_startup = self.config.capacity_var.split(".")[-1] + "_startup"
        aux_shutdown = self.config.capacity_var.split(".")[-1] + "_shutdown"
        
        setattr(self,aux_op_mode, Var())
        setattr(self,aux_startup,Var())
        setattr(self,aux_shutdown, Var())
        
        aux_var_op_mode = getattr(self,aux_op_mode)
        aux_var_startup = getattr(self,aux_startup)
        aux_var_shutdown = getattr(self,aux_shutdown)
        for blk in self.config.model_args:
            if isinstance(self.config.model_args[blk],DesignModelData):
                design_blk = self.config.model_args[blk]
        
        capacity = deepgetattr(design_blk,self.config.capacity_var)
        capacity_ub = deepgetattr(design_blk,self.config.capacity_var + ".ub")
        capacity_lb = deepgetattr(design_blk,self.config.capacity_var + ".lb")
        
        setattr(self, self.config.capacity_var.split(".")[-1] + "_linearization_con1",
                Constraint(
            expr = (design_blk.build_unit - self.op_mode - self.startup - self.shutdown)*capacity_lb <= 
                    capacity - aux_var_op_mode -  aux_var_startup - aux_var_shutdown
        ))

        setattr(self, self.config.capacity_var.split(".")[-1] + "_linearization_con2",
                Constraint(
            expr = (design_blk.build_unit - self.op_mode - self.startup - self.shutdown)*capacity_ub >= 
                    capacity - aux_var_op_mode -  aux_var_startup - aux_var_shutdown
        ))

        setattr(self,self.config.capacity_var.split(".")[-1] + "_linearization_con3", 
                 Constraint(
            expr = ( self.op_mode*capacity_lb <= aux_var_op_mode  )))
        
        setattr(self, self.config.capacity_var.split(".")[-1] + "_linearization_con4", 
                Constraint(
            expr = ( self.op_mode*capacity_ub >= aux_var_op_mode   )))
        
        setattr(self, self.config.capacity_var.split(".")[-1] + "_linearization_con5", 
                Constraint(
            expr = ( self.startup*capacity_lb <= aux_var_startup   )))
        
        setattr(self,self.config.capacity_var.split(".")[-1] + "_linearization_con6", 
                 Constraint(
            expr = ( self.startup*capacity_ub >= aux_var_startup )))
        
        setattr(self, self.config.capacity_var.split(".")[-1] + "_linearization_con7", 
                Constraint(
            expr = ( self.shutdown*capacity_lb <= aux_var_shutdown)))
        
        setattr(self, self.config.capacity_var.split(".")[-1] + "_linearization_con8", 
                Constraint(
            expr = ( self.shutdown*capacity_ub >= aux_var_shutdown )))



        

