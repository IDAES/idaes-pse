##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2020, by the
# software owners: The Regents of the University of California, through
# Lawrence Berkeley National Laboratory,  National Technology & Engineering
# Solutions of Sandia, LLC, Carnegie Mellon University, West Virginia
# University Research Corporation, et al. All rights reserved.
#
# Please see the files COPYRIGHT.txt and LICENSE.txt for full copyright and
# license information, respectively. Both files are also available online
# at the URL "https://github.com/IDAES/idaes-pse".
##############################################################################
"""
Methods for calculating rate constants
"""
from pyomo.environ import exp, Var, units as pyunits

from idaes.core import MaterialFlowBasis
from idaes.generic_models.properties.core.generic.generic_reaction import \
    ConcentrationForm
from idaes.core.util.constants import Constants as c
from idaes.core.util.exceptions import BurntToast, ConfigurationError


# -----------------------------------------------------------------------------
# Constant dh_rxn
class arrhenius():
    def build_parameters(rblock, config):
        parent = rblock.parent_block()
        base_units = parent.get_metadata().default_units
        e_units = (base_units["mass"] *
                   base_units["length"]**2 *
                   base_units["time"]**-2 *
                   base_units["amount"]**-1)

        rbasis = parent.config.reaction_basis
        if rbasis == MaterialFlowBasis.molar:
            r_base = base_units["amount"]
        elif rbasis == MaterialFlowBasis.mass:
            r_base = base_units["mass"]
        else:
            raise BurntToast(
                "{} for unexpected reaction basis {}. This should not happen "
                "so please contact the IDAES developers with this bug."
                .format(rblock.name, rbasis))

        c_form = config.concentration_form
        if c_form is None:
            raise ConfigurationError(
                "{} concentration_form configuration argument was not set. "
                "Please ensure that this argument is included in your "
                "configuration dict.".format(rblock.name))
        elif (c_form == ConcentrationForm.moleFraction or
              c_form == ConcentrationForm.massFraction):
            r_units = r_base*base_units["length"]**-3*base_units["time"]**-1
        else:
            order = 0
            for p, j in parent.config.property_package._phase_component_set:
                order += -rblock.reaction_order[p, j].value

            if (c_form == ConcentrationForm.molarity or
                    c_form == ConcentrationForm.activity):
                c_units = base_units["amount"]*base_units["length"]**-3
            elif c_form == ConcentrationForm.molality:
                c_units = base_units["amount"]*base_units["mass"]**-1
            elif c_form == ConcentrationForm.partialPressure:
                c_units = (base_units["mass"] *
                           base_units["length"]**-1 *
                           base_units["time"]**-2)
            else:
                raise BurntToast(
                    "{} get_concentration_term received unrecognised "
                    "ConcentrationForm ({}). This should not happen - please "
                    "contact the IDAES developers with this bug."
                    .format(rblock.name, c_form))

            r_units = (r_base *
                       base_units["length"]**-3 *
                       base_units["time"]**-1 *
                       c_units**order)

        rblock.arrhenius_const = Var(
                initialize=config.parameter_data["arrhenius_const"],
                doc="Arrhenius constant (pre-exponential factor)",
                units=r_units)

        rblock.energy_activation = Var(
                initialize=config.parameter_data["energy_activation"],
                doc="Activation energy",
                units=e_units)

    def return_expression(b, rblock, r_idx, T):
        base_units = rblock.parent_block().get_metadata().default_units
        R_units = (base_units["mass"] *
                   base_units["length"]**2 *
                   base_units["temperature"]**-1 *
                   base_units["amount"]**-1 *
                   base_units["time"]**-2)

        return rblock.arrhenius_const * exp(
            -rblock.energy_activation / (
                pyunits.convert(c.gas_constant, to_units=R_units)*T))
