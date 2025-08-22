#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Methods for calculating heat of reaction
"""
# TODO: Missing docstrings
# pylint: disable=missing-function-docstring

from pyomo.environ import units as pyunits, Var, value

from idaes.core import MaterialFlowBasis
from idaes.core.scaling import CustomScalerBase, get_scaling_factor
from idaes.core.util.constants import Constants
from idaes.core.util.misc import set_param_from_config


# -----------------------------------------------------------------------------
# Constant dh_rxn
class ConstantEnthalpyRxn(CustomScalerBase):
    """
    Scaler object for the constant_dh_rxn method for calculating
    the heat of reaction
    """

    def variable_scaling_routine(
        self,
        model,
        reaction,
        overwrite: bool = False,
    ):
        rblock = getattr(model.params, "reaction_" + reaction)
        units = model.params.get_metadata().derived_units
        sf_T = get_scaling_factor(model.state_ref.temperature)
        sf_R = value(
            1 / pyunits.convert(Constants.gas_constant, to_units=units["gas_constant"])
        )

        if model.is_property_constructed("dh_rxn"):
            self.set_component_scaling_factor(
                model.dh_rxn[reaction], sf_T * sf_R, overwrite=overwrite
            )

    def constraint_scaling_routine(self, model, reaction, overwrite: bool = False):
        # No constraints generated
        pass


class constant_dh_rxn:
    """Methods for constant heat of reaction."""

    default_scaler = ConstantEnthalpyRxn

    @staticmethod
    def build_parameters(rblock, config):
        units = rblock.parent_block().get_metadata().derived_units

        rbasis = rblock.parent_block().config.reaction_basis
        if rbasis == MaterialFlowBasis.molar:
            basis = "mole"
        elif rbasis == MaterialFlowBasis.mass:
            basis = "mass"

        rblock.dh_rxn_ref = Var(
            doc="Specific heat of reaction at reference state",
            units=units["energy_" + basis],
        )

        set_param_from_config(rblock, param="dh_rxn_ref", config=config)

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        return rblock.dh_rxn_ref

    @staticmethod
    def calculate_scaling_factors(b, rblock):
        v = abs(value(rblock.dh_rxn_ref))

        # Need to make sure dh_rxn is not 0 to avoid division by 0
        if v != 0:
            return 1 / abs(value(rblock.dh_rxn_ref))
        else:
            return 1
