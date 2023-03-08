#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Methods for defining reaction rates
"""
# TODO: Missing docstrings
# pylint: disable=missing-function-docstring

from idaes.models.properties.modular_properties.base.utility import (
    get_concentration_term,
)


# -----------------------------------------------------------------------------
class power_law_rate:
    """Methods for power-law rate expressions."""

    @staticmethod
    def build_parameters(rblock, config):
        pass

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        e = None
        # Get reaction orders and construct power law expression
        for p, j in b.phase_component_set:
            o = rblock.reaction_order[p, j]

            if e is None and o.value != 0:
                e = get_concentration_term(b, r_idx)[p, j] ** o
            elif e is not None and o.value != 0:
                e = e * get_concentration_term(b, r_idx)[p, j] ** o

        return b.k_rxn[r_idx] * e
