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
from pyomo.environ import exp, Var

from idaes.core.util.constants import Constants as c


# -----------------------------------------------------------------------------
# Constant dh_rxn
class arrhenius():
    def build_parameters(rblock, config):
        rblock.arrhenius_const = Var(
                initialize=config.parameter_data["arrhenius_const"],
                doc="Arrhenius constant (pre-exponential factor)")

        rblock.energy_activation = Var(
                initialize=config.parameter_data["energy_activation"],
                doc="Activation energy")

    def return_expression(b, rblock, r_idx, T):
        return rblock.arrhenius_const * exp(
            -rblock.energy_activation / (c.gas_constant*T))
