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

from idaes.core.util.constants import Constants as c


# -----------------------------------------------------------------------------
# Constant dh_rxn
class van_t_hoff():
    def build_parameters(rblock, config):
        base_units = rblock.parent_block().get_metadata().default_units

        rblock.k_eq_ref = Var(
                initialize=config.parameter_data["k_eq_ref"],
                doc="Equilibrium constant at reference state")

        rblock.T_eq_ref = Var(
                initialize=config.parameter_data["T_eq_ref"],
                doc="Reference temperature for equilibrium constant",
                units=base_units["temperature"])

    def return_expression(b, rblock, r_idx, T):
        units = rblock.parent_block().get_metadata().derived_units

        return rblock.k_eq_ref * exp(
            -(b.dh_rxn[r_idx] /
              pyunits.convert(c.gas_constant,
                              to_units=units["gas_constant"])) *
            (1/T - 1/rblock.T_eq_ref))
