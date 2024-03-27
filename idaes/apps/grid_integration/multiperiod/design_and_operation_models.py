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
    try:
        return reduce(getattr, attr.split("."), obj)
    except:
        raise ValueError(
            f"The attribute `{attr}` does not exist on the {obj.__class__.__name__} instance passed. Make sure that `{attr}` is defined for that object."
        )


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

    # noinspection PyAttributeOutsideInit
    def build(self):
        super().build()

        self.config.model_func(self, **self.config.model_args)


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
            doc="Should LMP data automatically be appended to the model?",
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
                doc="Parameter: Will be updated to LMP value at given time",
            )

        self.config.model_func(self, **self.config.model_args)
