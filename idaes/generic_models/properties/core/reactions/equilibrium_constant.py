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
Methods for calculating equilibrium constants
"""
from pyomo.environ import exp, Var, units as pyunits

from idaes.generic_models.properties.core.generic.generic_reaction import \
    ConcentrationForm
from idaes.core.util.misc import set_param_from_config
from idaes.core.util.constants import Constants as c
from idaes.core.util.exceptions import BurntToast, ConfigurationError


# -----------------------------------------------------------------------------
# Constant dh_rxn
class van_t_hoff():

    @staticmethod
    def build_parameters(rblock, config):
        parent = rblock.parent_block()
        units = parent.get_metadata().derived_units

        c_form = config.concentration_form
        if c_form is None:
            raise ConfigurationError(
                "{} concentration_form configuration argument was not set. "
                "Please ensure that this argument is included in your "
                "configuration dict.".format(rblock.name))
        elif (c_form == ConcentrationForm.moleFraction or
              c_form == ConcentrationForm.massFraction):
            e_units = None
        else:
            order = 0
            for p, j in parent.config.property_package._phase_component_set:
                order += rblock.reaction_order[p, j].value

            if (c_form == ConcentrationForm.molarity or
                    c_form == ConcentrationForm.activity):
                c_units = units["density_mole"]
            elif c_form == ConcentrationForm.molality:
                c_units = units["amount"]*units["mass"]**-1
            elif c_form == ConcentrationForm.partialPressure:
                c_units = units["pressure"]
            else:
                raise BurntToast(
                    "{} get_concentration_term received unrecognised "
                    "ConcentrationForm ({}). This should not happen - please "
                    "contact the IDAES developers with this bug."
                    .format(rblock.name, c_form))

            e_units = c_units**order

        rblock.k_eq_ref = Var(
                doc="Equilibrium constant at reference state",
                units=e_units)
        set_param_from_config(rblock, param="k_eq_ref", config=config)

        rblock.T_eq_ref = Var(
                doc="Reference temperature for equilibrium constant",
                units=units["temperature"])
        set_param_from_config(rblock, param="T_eq_ref", config=config)

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        units = rblock.parent_block().get_metadata().derived_units

        return rblock.k_eq_ref * exp(
            -(b.dh_rxn[r_idx] /
              pyunits.convert(c.gas_constant,
                              to_units=units["gas_constant"])) *
            (1/T - 1/rblock.T_eq_ref))
