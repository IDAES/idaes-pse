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

from idaes.core.util.constants import Constants as c


# -----------------------------------------------------------------------------
# Constant dh_rxn
class arrhenius():
    def build_parameters(rblock, config):
        base_units = rblock.parent_block().get_metadata().default_units
        e_units = (base_units["mass"] *
                   base_units["length"]**2 *
                   base_units["time"]**-2 *
                   base_units["amount"]**-1)

        rblock.arrhenius_const = Var(
                initialize=config.parameter_data["arrhenius_const"],
                doc="Arrhenius constant (pre-exponential factor)")

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
