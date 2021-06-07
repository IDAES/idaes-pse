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
Methods for defining equilibrium reactions
"""
from pyomo.environ import Param, units as pyunits

from idaes.core.util.math import safe_log, smooth_max
from idaes.core.util.exceptions import ConfigurationError

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


# ----------------------------------------------------------------------------
class solubility_product():
    """
    Complementariity formulation for solid precipitation
    Thanks to Larry Biegler for the formulation

    Like any phase equilibrium, solid precipitation is complicated by
    conditional constraint. For precipitation, the solubility product,
    Ksp = f(C) only holds if solids are present in the system. That is, if the
    system is subsaturated and no solids are present, Ksp >= f(C) instead.

    Letting S represent solids in the system and Q = Ksp - f(C), we have the
    following:

        * If S >= 0, Q = 0
        * If S = 0, Q >= 0

    Thus, only one of S and Q can be greater than zero at any time.
    This can be written in the form of a complementarity constraint as:

        * S - MAX(0, S-Q) == 0

    where S is assumed to be the sum of the flowrates any solids formed in the
    reaction. This allows for multiple solid products, and only applies the
    solubility product if all flowrates are non-zero.

    Note that the smooth maximum approximation requires the use of a smoothing
    factor eps, and that the value of this factor needs to be adjusted for the
    specific problem. as an initial guideline, it is suggested that this be at
    least one to two orders of magnitude smaller than the solubility product.
    """

    @staticmethod
    def build_parameters(rblock, config):
        rblock.eps = Param(mutable=True,
                           initialize=1e-4,
                           doc="Smoothing parameter for smooth maximum")

    @staticmethod
    def return_expression(b, rblock, r_idx, T):
        e = None
        s = None

        if hasattr(b.params, "_electrolyte") and b.params._electrolyte:
            pc_set = b.params.true_phase_component_set
        else:
            pc_set = b.phase_component_set

        # Get reaction orders and construct power law expression
        for p, j in pc_set:
            p_obj = b.state_ref.params.get_phase(p)

            # First, build E
            o = rblock.reaction_order[p, j]

            if e is None and o.value != 0:
                e = get_concentration_term(b, r_idx)[p, j]**o
            elif e is not None and o.value != 0:
                e = e*get_concentration_term(b, r_idx)[p, j]**o

            if p_obj.is_solid_phase():
                # If solid phase, identify S
                r_config = b.params.config.equilibrium_reactions[r_idx]

                try:
                    stoic = r_config.stoichiometry[p, j]
                    if s is None and stoic != 0:
                        s = b.state_ref.flow_mol_phase_comp[p, j]
                    elif s is not None and stoic != 0:
                        s += b.state_ref.flow_mol_phase_comp[p, j]
                except KeyError:
                    pass

        if s is None:
            # Catch for not finding a solid phase
            raise ConfigurationError(
                "{} did not find a solid phase component for precipitation "
                "reaction {}. This is likely due to the reaction "
                "configuration.".format(b.name, r_idx))

        Q = b.k_eq[r_idx] - e

        return s - smooth_max(0, s-Q, rblock.eps) == 0

    @staticmethod
    def calculate_scaling_factors(b, sf_keq):
        return sf_keq
