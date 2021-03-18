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
Methods for defining equibilibrium reactions
"""
from idaes.core.util.math import safe_log

from idaes.generic_models.properties.core.generic.generic_reaction import \
    get_concentration_term


# Smooth parameter to use in safe_log approximations
# Needs to be small due to small magnitude of many equilibrium constants
EPS = 1e-15


# ----------------------------------------------------------------------------
class power_law_equil():

    @staticmethod
    def build_parameters(rblock, config):
        pass

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        e = None

        if hasattr(b.params, "_electrolyte") and b.params._electrolyte:
            pc_set = b.params.true_phase_component_set
        else:
            pc_set = b.phase_component_set

        # Get reaction orders and construct power law expression
        for p, j in pc_set:
            o = rblock.reaction_order[p, j]

            if e is None and o != 0:
                e = get_concentration_term(b, r_idx)[p, j]**o
            elif e is not None and o != 0:
                e = e*get_concentration_term(b, r_idx)[p, j]**o

        return b.k_eq[r_idx] == e


# ----------------------------------------------------------------------------
class log_power_law_equil():

    @staticmethod
    def build_parameters(rblock, config):
        pass

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        e = None

        if hasattr(b.params, "_electrolyte") and b.params._electrolyte:
            pc_set = b.params.true_phase_component_set
        else:
            pc_set = b.phase_component_set

        # Get reaction orders and construct power law expression
        for p, j in pc_set:
            o = rblock.reaction_order[p, j]

            if e is None and o != 0:
                e = o*safe_log(get_concentration_term(b, r_idx)[p, j], eps=EPS)
            elif e is not None and o != 0:
                e = e + o*safe_log(
                    get_concentration_term(b, r_idx)[p, j], eps=EPS)

        return safe_log(b.k_eq[r_idx], eps=EPS) == e
