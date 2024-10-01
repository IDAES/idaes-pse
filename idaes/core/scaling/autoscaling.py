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
Tools for automatically scaling models based on current state.

Author: Andrew Lee
"""

import scipy as sp

from pyomo.environ import (
    Block,
    Constraint,
    value,
    Var,
)
from pyomo.core.base.block import BlockData
from pyomo.core.base.constraint import ConstraintData
from pyomo.core.base.var import VarData

from idaes.core.scaling.scaling_base import ScalerBase
from idaes.core.scaling.util import get_scaling_factor
from idaes.core.util.scaling import get_jacobian


class AutoScaler(ScalerBase):
    """
    IDAES Autoscaling Toolbox

    Contains a number of methods useful for automatically scaling models based
    on the current model state (i.e., variable values). Users should be aware
    of the limitations of autoscaling however, which only consider the current variable
    values and are thus heavily dependent on model initialization.

    """

    CONFIG = ScalerBase.CONFIG()

    def scale_variables_by_magnitude(self, blk_or_var, descend_into: bool = True):
        """
        Calculate scaling factors for all variables in a model based on their
        current magnitude. Variables with no value are assigned a scaling factor of 1.

        Args:
            blk_or_var: block or variable object to calculate scaling factors for
            descend_into: if blk_or_var is a Block, whether to descend into any sub-Blocks
              (default=True)

        Returns:
            None

        Raises:
            TypeError is blk_or_var is not a Block or Var.
        """
        if isinstance(blk_or_var, BlockData):
            # Scalar Block or element of Indexed Block
            # As scaling factors live with the components parent, do not descend into
            # sub-blocks here
            for v in blk_or_var.component_data_objects(Var, descend_into=False):
                self._vardata_by_magnitude(v)

            # Next, get all child blocks and call scale_variables_by_magnitude recursively
            if descend_into:
                for b in blk_or_var.component_data_objects(Block, descend_into=False):
                    self.scale_variables_by_magnitude(b, descend_into=descend_into)

        elif isinstance(blk_or_var, Block):
            # Indexed Block
            for b in blk_or_var.values():
                for v in b.component_data_objects(Var, descend_into=descend_into):
                    self._vardata_by_magnitude(v)

        elif isinstance(blk_or_var, VarData):
            # Scalar Var or element of Indexed Var
            self._vardata_by_magnitude(blk_or_var)

        elif isinstance(blk_or_var, Var):
            # Indexed Var
            for v in blk_or_var.values():
                self._vardata_by_magnitude(v)
        else:
            raise TypeError(f"{blk_or_var.name} is not a block or variable.")

    def scale_constraints_by_jacobian_norm(
        self, blk_or_cons, norm: int = 2, descend_into: bool = True
    ):
        """
        Calculate scaling factors for all constraints in a model based on the norm of
        the Jacobian matrix, accounting for any variable scaling factors.

        Args:
            blk_or_cons: block or constraint to calculate scaling factors for
            norm: type of norm to use for scaling. Must be a positive integer.
            descend_into: if blk_or_cons is a Block, whether to descend into any sub-Blocks
              (default=True)

        Returns:
            None

        Raises:
            TypeError is blk_or_cons is not a Block or Constraint
            ValueError if norm is not a positive integer
        """
        # Validate norm
        # First cast to int
        norm = int(norm)
        if norm < 1:
            raise ValueError(
                f"Invalid value for norm in scale_constraints_by_jacobian_norm ({norm}). "
                "Value must be a positive integer."
            )

        # We want to avoid generating the Jacobian and NLP more than once, so first we
        # will identify the top-level block and collect all the constraints of interest
        # as a list or iterator

        if isinstance(blk_or_cons, BlockData):
            # Scalar Block or element of Indexed Block
            if descend_into:
                # Scale all constraints, so pass con_list=None
                con_list = None
            else:
                # Otherwise, get an iterator of constraints only on block
                con_list = blk_or_cons.component_data_objects(
                    Constraint, descend_into=False
                )
            jblock = blk_or_cons

        elif isinstance(blk_or_cons, Block):
            # Indexed Block
            # Underlying tools do not work for Indexed blocks, so
            # use parent block instead and collect constraints of interest
            con_list = []
            for b in blk_or_cons.values():
                for c in b.component_data_objects(
                    Constraint, descend_into=descend_into
                ):
                    con_list.append(c)
            jblock = blk_or_cons.parent_block()

        elif isinstance(blk_or_cons, ConstraintData):
            # Scalar Constraint or element of Indexed Constraint
            con_list = [blk_or_cons]
            jblock = blk_or_cons.parent_block()

        elif isinstance(blk_or_cons, Constraint):
            # Indexed Constraint
            con_list = blk_or_cons.values()
            jblock = blk_or_cons.parent_block()

        else:
            raise TypeError(f"{blk_or_cons.name} is not a block or constraint.")

        # Once we have a single target block and list of constraints, call the scaler method
        # once for all the constraints
        self._con_by_norm(jblock, con_list=con_list, norm=norm)

    def scale_model(self, model, norm: int = 2, descend_into: bool = True):
        """
        Apply auto-scaling routine to model.

        Args:
            model: model to be scaled
            norm: type of norm to use for scaling. Must be a positive integer.
            descend_into: if sub-Blocks are present, whether to descend into any sub-Blocks

        Returns:
            None
        """
        self.scale_variables_by_magnitude(blk_or_var=model, descend_into=descend_into)
        self.scale_constraints_by_jacobian_norm(
            blk_or_cons=model, norm=norm, descend_into=descend_into
        )

    def _vardata_by_magnitude(self, vardata):
        if vardata.value is None:
            sf = 1.0
        else:
            val = abs(value(vardata))
            if val <= self.config.zero_tolerance:
                sf = 1.0
            else:
                sf = 1.0 / val

        self._set_scaling_factor(vardata, "variable", sf)

    def _con_by_norm(self, blk, con_list=None, norm=2):
        # Get scaled Jacobian - we want to consider any existing scaling
        # We will account for existing scaling factors later and update them
        jac, nlp = get_jacobian(blk, scaled=True)

        if con_list is None:
            con_list = nlp.get_pyomo_equality_constraints()

        # Use scipy to get all the norms
        # Should be more efficient that iterating in Python
        axis = (
            1  # Could make this an argument to also support variable-based norm scaling
        )
        if jac.format == "csr":
            jac_norms = sp.sparse.linalg.norm(jac, ord=norm, axis=axis)
        else:
            jac_norms = sp.linalg.norm(jac, ord=norm, axis=axis)

        # Iterate over constraints of interest and apply scaling factors
        for c in con_list:
            c_idx = nlp.get_pyomo_equality_constraints().index(c)

            # Get any existing scaling factor for this constraint
            sf_old = get_scaling_factor(c)
            if sf_old is None:
                sf_old = 1.0

            # Get norm for this constraint
            n = jac_norms[c_idx]

            if n <= self.config.zero_tolerance:
                sf = sf_old
            else:
                sf = sf_old / n

            self._set_scaling_factor(c, "constraint", sf)
