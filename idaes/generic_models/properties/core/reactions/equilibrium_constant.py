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
Methods for calculating equilibrium constants
"""
from pyomo.environ import exp, Var, units as pyunits, value

from idaes.generic_models.properties.core.generic.generic_reaction import \
    ConcentrationForm
from .dh_rxn import constant_dh_rxn
from idaes.core.util.misc import set_param_from_config
from idaes.core.util.constants import Constants as c
from idaes.core.util.exceptions import BurntToast, ConfigurationError


# -----------------------------------------------------------------------------
# Constant Keq
class ConstantKeq():

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

            try:
                # This will work for Reaction Packages
                pc_set = parent.config.property_package._phase_component_set
            except AttributeError:
                # Need to allow for inherent reactions in Property Packages
                if not parent._electrolyte:
                    # In most cases ,should have _phase_component_set
                    pc_set = parent._phase_component_set
                else:
                    # However, for electrolytes need true species set
                    pc_set = parent.true_phase_component_set

            for p, j in pc_set:
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

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        return rblock.k_eq_ref

    @staticmethod
    def calculate_scaling_factors(b, rblock):
        return 1/value(rblock.k_eq_ref)


# -----------------------------------------------------------------------------
# van t'Hoff equation (constant dh_rxn)
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

            try:
                # This will work for Reaction Packages
                pc_set = parent.config.property_package._phase_component_set
            except AttributeError:
                # Need to allow for inherent reactions in Property Packages
                if not parent._electrolyte:
                    # In most cases ,should have _phase_component_set
                    pc_set = parent._phase_component_set
                else:
                    # However, for electrolytes need true species set
                    pc_set = parent.true_phase_component_set

            for p, j in pc_set:
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

    @staticmethod
    def calculate_scaling_factors(b, rblock):
        return 1/value(rblock.k_eq_ref)


# -----------------------------------------------------------------------------
# Constant dh_rxn and ds_rxn
class gibbs_energy():

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
        elif (c_form == ConcentrationForm.molarity or
              c_form == ConcentrationForm.activity):

            order = 0
            try:
                # This will work for Reaction Packages
                pc_set = parent.config.property_package._phase_component_set
            except AttributeError:
                # Need to allow for inherent reactions in Property Packages
                if not parent._electrolyte:
                    # In most cases ,should have _phase_component_set
                    pc_set = parent._phase_component_set
                else:
                    # However, for electrolytes need true species set
                    pc_set = parent.true_phase_component_set

            for p, j in pc_set:
                order += rblock.reaction_order[p, j].value

            rblock._keq_units = (pyunits.convert(1*pyunits.mol/pyunits.L,
                                                 units["density_mole"]))**order
        else:
            raise ConfigurationError(
                "{} calculation of equilibrium constant based on Gibbs energy "
                "is only supported for molarity or activity forms. "
                "Currently selected form: {}".format(rblock.name, c_form))

        # Check that heat of reaction is constant
        if config.heat_of_reaction is not constant_dh_rxn:
            raise ConfigurationError(
                "{} calculating equilibrium constants from Gibbs energy "
                "assumes constant heat of reaction. Please ensure you are "
                "using the constant_dh_rxn method for this reaction"
                .format(rblock.name))

        rblock.ds_rxn_ref = Var(
                doc="Specific molar entropy of reaction at reference state",
                units=units["entropy_mole"])
        set_param_from_config(rblock, param="ds_rxn_ref", config=config)

        rblock.T_eq_ref = Var(
                doc="Reference temperature for equilibrium constant",
                units=units["temperature"])
        set_param_from_config(rblock, param="T_eq_ref", config=config)

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        units = rblock.parent_block().get_metadata().derived_units

        R = pyunits.convert(c.gas_constant, to_units=units["gas_constant"])

        return (exp(
            (-rblock.dh_rxn_ref / (R*T)) +
            (rblock.ds_rxn_ref / R)) * rblock._keq_units)

    @staticmethod
    def calculate_scaling_factors(b, rblock):
        units = rblock.parent_block().get_metadata().derived_units
        R = pyunits.convert(c.gas_constant, to_units=units["gas_constant"])

        keq_val = value(exp(-rblock.dh_rxn_ref/(R*rblock.T_eq_ref) +
                            rblock.ds_rxn_ref/R) * rblock._keq_units)
        return 1/keq_val
