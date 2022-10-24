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
"""
Methods for calculating rate constants
"""
from pyomo.environ import exp, Var, units as pyunits

from idaes.core import MaterialFlowBasis
from idaes.models.properties.modular_properties.base.utility import ConcentrationForm
from idaes.core.util.misc import set_param_from_config
from idaes.core.util.constants import Constants as c
from idaes.core.util.exceptions import BurntToast, ConfigurationError


# -----------------------------------------------------------------------------
# Constant dh_rxn
class arrhenius:
    @staticmethod
    def build_parameters(rblock, config):
        parent = rblock.parent_block()
        units = parent.get_metadata().derived_units

        rbasis = parent.config.reaction_basis
        if rbasis == MaterialFlowBasis.molar:
            r_base = units.AMOUNT
        elif rbasis == MaterialFlowBasis.mass:
            r_base = units.MASS
        else:
            raise BurntToast(
                "{} for unexpected reaction basis {}. This should not happen "
                "so please contact the IDAES developers with this bug.".format(
                    rblock.name, rbasis
                )
            )

        c_form = config.concentration_form
        if c_form is None:
            raise ConfigurationError(
                "{} concentration_form configuration argument was not set. "
                "Please ensure that this argument is included in your "
                "configuration dict.".format(rblock.name)
            )
        elif (
            c_form == ConcentrationForm.moleFraction
            or c_form == ConcentrationForm.massFraction
            or c_form == ConcentrationForm.activity
        ):
            r_units = r_base * units.VOLUME**-1 * units.TIME**-1
        else:
            order = 0
            for p, j in parent.config.property_package._phase_component_set:
                order += -rblock.reaction_order[p, j].value

            if c_form == ConcentrationForm.molarity:
                c_units = units.DENSITY_MOLE
            elif c_form == ConcentrationForm.molality:
                c_units = units.AMOUNT * units.MASS**-1
            elif c_form == ConcentrationForm.partialPressure:
                c_units = units.PRESSURE
            else:
                raise BurntToast(
                    "{} received unrecognised ConcentrationForm ({}). "
                    "This should not happen - please contact the IDAES "
                    "developers with this bug.".format(rblock.name, c_form)
                )

            r_units = r_base * units.LENGTH**-3 * units.TIME**-1 * c_units**order

        rblock.arrhenius_const = Var(
            doc="Arrhenius constant (pre-exponential factor)", units=r_units
        )
        set_param_from_config(rblock, param="arrhenius_const", config=config)

        rblock.energy_activation = Var(doc="Activation energy", units=units.ENERGY_MOLE)

        set_param_from_config(rblock, param="energy_activation", config=config)

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        units = rblock.parent_block().get_metadata().derived_units

        return rblock.arrhenius_const * exp(
            -rblock.energy_activation
            / (pyunits.convert(c.gas_constant, to_units=units.GAS_CONSTANT) * T)
        )
