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
"""
Tool for automatically scaling arc constraints on a Pyomo block.

Author: Douglas Allan
"""

from pyomo.environ import Constraint
from pyomo.network import Arc

import idaes.logger as idaeslog
from idaes.core.scaling.custom_scaler_base import (
    ConstraintScalingScheme,
    CustomScalerBase,
)

_log = idaeslog.getLogger(__name__)


class ArcConstraintScaler(CustomScalerBase):
    """
    Scaler object to automatically scale all arc constraints on a Pyomo block
    (such as a Flowsheet object) using scaling methods from CustomScalerBase.

    Config arguments from ScalerBase and CustomScaler base should be set to,
    for example, make sure that the maximum and minimum scaling thresholds
    are set appropriately.
    """

    def scale_model(self, **kwargs):
        """
        Do not use this method, use one of the scale_arc_constraints
        methods instead.
        """
        raise NotImplementedError(
            "The scale_model method is not implemented for the ArcConstraintScaler. "
            "Use one of the scale_arc_constraints methods instead."
        )

    def scale_arc_constraints_by_nominal_value(
        self,
        blk,
        scheme=ConstraintScalingScheme.inverseMaximum,
        overwrite: bool = False,
        descend_into: bool = False,
    ):
        """
        Iterate over a block to find Arc objects present, then scale
        the corresponding arc constraints using one of the supported
        nominal value constraint scaling schemes.

        Args:
            blk: Pyomo Block to iterate over to find arc constraints
            scheme: ConstraintScalingScheme Enum indicating method to apply for
              determining constraint scaling. The default is inverseMaximum.
            overwrite: Whether to overwrite existing scaling factors.
              The default is False.
            descend_into: Whether to recursively descend into child blocks to
              search for arc constraints there to scale. The default is False.
        Returns:
            None
        """

        for arc in blk.component_data_objects(ctype=Arc, descend_into=descend_into):
            for con in arc.expanded_block.component_data_objects(ctype=Constraint):
                self.scale_constraint_by_nominal_value(
                    con, scheme=scheme, overwrite=overwrite
                )

    def scale_arc_constraints_by_nominal_derivative_norm(
        self, blk, norm: int = 2, overwrite: bool = False, descend_into: bool = False
    ):
        """
        Iterate over a block to find Arc objects present, then scale
        the corresponding arc constraints using the norm of the
        derivative of the constraint.

        Args:
            blk: Pyomo Block to iterate over to find arc constraints
            norm: Type of norm to use for scaling. Must be a positive integer.
              The default is 2.
            overwrite: Whether to overwrite existing scaling factors. The
              default is False.
            descend_into: Whether to recursively descend into child blocks to
              search for arc constraints there to scale. The default is False.
        Returns:
            None
        """
        for arc in blk.component_data_objects(ctype=Arc, descend_into=descend_into):
            for con in arc.expanded_block.component_data_objects(ctype=Constraint):
                self.scale_constraint_by_nominal_derivative_norm(
                    con,
                    norm=norm,
                    overwrite=overwrite,
                )
