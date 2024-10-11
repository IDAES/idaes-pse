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
"""
Base class for Scalers

Author: Andrew Lee
"""

from pyomo.common.config import (
    Bool,
    ConfigDict,
    ConfigValue,
    String_ConfigFormatter,
)
from pyomo.core.base.constraint import ConstraintData
from pyomo.core.base.var import VarData

from idaes.core.scaling.util import get_scaling_factor, set_scaling_factor
import idaes.logger as idaeslog

# Set up logger
_log = idaeslog.getLogger(__name__)


# Common base ConfigBlock for all Scalers
CONFIG = ConfigDict()
CONFIG.declare(
    "zero_tolerance",
    ConfigValue(
        default=1e-12,
        domain=float,
        description="Value at which a variable will be considered equal to zero for scaling.",
    ),
)
CONFIG.declare(
    "max_variable_scaling_factor",
    ConfigValue(
        default=1e10,
        domain=float,
        description="Maximum value for variable scaling factors.",
    ),
)
CONFIG.declare(
    "min_variable_scaling_factor",
    ConfigValue(
        default=1e-10,
        domain=float,
        description="Minimum value for variable scaling factors.",
    ),
)
CONFIG.declare(
    "max_constraint_scaling_factor",
    ConfigValue(
        default=1e10,
        domain=float,
        description="Maximum value for constraint scaling factors.",
    ),
)
CONFIG.declare(
    "min_constraint_scaling_factor",
    ConfigValue(
        default=1e-10,
        domain=float,
        description="Minimum value for constraint scaling factors.",
    ),
)
CONFIG.declare(
    "overwrite",
    ConfigValue(
        default=False,
        domain=Bool,
        description="Whether to overwrite existing scaling factors.",
    ),
)


class ScalerBase:
    """
    Base class for IDAES Scaler objects

    Contains a number of methods useful for scaling models.

    """

    CONFIG = CONFIG()

    def __init__(self, **kwargs):
        self.config = self.CONFIG(kwargs)

    def __init_subclass__(cls, **kwargs):
        super().__init_subclass__(**kwargs)

        # Handle cases where a class doc string was not set
        if cls.__doc__ is None:
            cls.__doc__ = ""

        cls.__doc__ = cls.__doc__ + cls.CONFIG.generate_documentation(
            format=String_ConfigFormatter(
                block_start="%s\n",
                block_end="",
                item_start="%s\n",
                item_body="%s",
                item_end="\n",
            ),
            indent_spacing=4,
            width=66,
        )

    def get_scaling_factor(self, component):
        """
        Get scaling factor for component.

        Alias for get_scaling_factor utility function.

        Args:
            component: component to get scaling factor for

        Returns:
            float - scaling factor

        Raises:
            TypeError if component is a Block
        """
        return get_scaling_factor(component)

    def set_variable_scaling_factor(
        self, variable, scaling_factor: float, overwrite: bool = None
    ):
        """
        Set scaling factor for variable.

        Scaling factor is limited by min_variable_scaling_factor and max_variable_scaling_factor.

        Args:
            variable: VarData component to set scaling factor for.
            scaling_factor: nominal scaling factor to apply. May be limited by max and min values.
            overwrite: whether to overwrite existing scaling factor (if present).
              Defaults to Scaler config setting.

        Returns:
            None

        Raises:
            TypeError if variable is not an instance of VarData
        """
        if not isinstance(variable, VarData):
            raise TypeError(f"{variable} is not a variable (or is indexed).")
        self._set_scaling_factor(
            component=variable,
            component_type="variable",
            scaling_factor=scaling_factor,
            overwrite=overwrite,
        )

    def set_constraint_scaling_factor(
        self, constraint, scaling_factor: float, overwrite: bool = None
    ):
        """
        Set scaling factor for constraint.

        Scaling factor is limited by min_constraint_scaling_factor and max_constraint_scaling_factor.

        Args:
            constraint: ConstraintData component to set scaling factor for.
            scaling_factor: nominal scaling factor to apply. May be limited by max and min values.
            overwrite: whether to overwrite existing scaling factor (if present).
              Defaults to Scaler config setting.

        Returns:
            None

        Raises:
            TypeError if constraint is not an instance of ConstraintData
        """
        if not isinstance(constraint, ConstraintData):
            raise TypeError(f"{constraint} is not a constraint (or is indexed).")
        self._set_scaling_factor(
            component=constraint,
            component_type="constraint",
            scaling_factor=scaling_factor,
            overwrite=overwrite,
        )

    def _set_scaling_factor(
        self, component, component_type, scaling_factor, overwrite=None
    ):
        """
        PRIVATE METHOD

        The purpose of this method is to apply the correct max and min limits to scaling factors
        before setting them. Which set of limits to apply is determined solely by the component_type
        argument, and no type checking is performed.
        """
        # If overwrite not provided, use scaler config
        # This allows developers to have some more control over when to overwrite or not
        if overwrite is None:
            overwrite = self.config.overwrite

        if component_type == "variable":
            maxsf = self.config.max_variable_scaling_factor
            minsf = self.config.min_variable_scaling_factor
        elif component_type == "constraint":
            maxsf = self.config.max_constraint_scaling_factor
            minsf = self.config.min_constraint_scaling_factor
        else:
            raise ValueError("Invalid value for component_type.")

        if scaling_factor > maxsf:
            _log.debug(
                f"Scaling factor for {component.name} limited by maximum value "
                f"(max_sf: {maxsf} < sf: {scaling_factor})"
            )
            scaling_factor = maxsf
        elif scaling_factor < minsf:
            _log.debug(
                f"Scaling factor for {component.name} limited by minimum value "
                f"(min_sf: {minsf} > sf: {scaling_factor})"
            )
            scaling_factor = minsf

        set_scaling_factor(component, scaling_factor, overwrite=overwrite)
