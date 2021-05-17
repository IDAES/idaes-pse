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
from pyomo.environ import Param, units as pyunits

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

            if e is None and o.value != 0:
                e = get_concentration_term(b, r_idx)[p, j]**o
            elif e is not None and o.value != 0:
                e = e*get_concentration_term(b, r_idx)[p, j]**o

        return b.k_eq[r_idx] == e

    @staticmethod
    def calculate_scaling_factors(b, sf_keq):
        return sf_keq


# ----------------------------------------------------------------------------
class log_power_law_equil():

    @staticmethod
    def build_parameters(rblock, config):
        rblock.eps = Param(default=EPS,
                           mutable=True,
                           doc="Smoothing factor for safe log function")

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

            if e is None and o.value != 0:
                # Need to strip units from concentration term (if applicable)
                c = get_concentration_term(b, r_idx)[p, j]
                u = pyunits.get_units(c)
                if u is not None:
                    # Has units, so divide conc by units
                    expr = c/u
                else:
                    # Units is None, so just use conc
                    expr = c
                e = o*safe_log(expr, eps=rblock.eps)
            elif e is not None and o.value != 0:
                # Need to strip units from concentration term (if applicable)
                c = get_concentration_term(b, r_idx)[p, j]
                u = pyunits.get_units(c)
                if u is not None:
                    # Has units, so divide conc by units
                    expr = c/u
                else:
                    # Units is None, so just use conc
                    expr = c
                e = e + o*safe_log(expr, eps=rblock.eps)

        # Need to check units on k_eq as well
        u = pyunits.get_units(b.k_eq[r_idx])
        if u is not None:
            # Has units, so divide k_eq by units
            expr = b.k_eq[r_idx]/u
        else:
            # Units is None, so just use k_eq
            expr = b.k_eq[r_idx]

        return safe_log(expr, eps=rblock.eps) == e

    @staticmethod
    def calculate_scaling_factors(b, sf_keq):
        return 1
