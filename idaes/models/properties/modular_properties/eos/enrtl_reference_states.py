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
Reference state sub-methods for eNRTL activity coefficient method.

Only applicable to liquid/electrolyte phases

Reference state classes should implement the following methods:
    - ref_state : expressions for x and X at the reference state
    - ndIdn : derivative term for long-range contribution

Reference:

Song, Y. and Chen, C.-C., Symmetric Electrolyte Nonrandom Two-Liquid Activity
Coefficient Model, Ind. Eng. Chem. Res., 2009, Vol. 48, pgs. 7788–7797
"""
from pyomo.environ import Expression

from idaes.models.properties.modular_properties.base.utility import (
    get_component_object as cobj,
)
import idaes.logger as idaeslog


# Set up logger
_log = idaeslog.getLogger(__name__)


# Near-zero value to use for unsymmetric reference state
EPS = 1e-20


class Unsymmetric(object):
    """
    Sub-methods for the symmetric (fused-salt) reference state
    """

    @staticmethod
    def ref_state(b, pname):
        def rule_x_ref(b, i):
            if i in b.params.solvent_set or i in b.params.solute_set:
                # Eqn 66
                return b.mole_frac_phase_comp_true[pname, i] / sum(
                    b.mole_frac_phase_comp_true[pname, j]
                    for j in b.params.solvent_set | b.params.solute_set
                )
            else:
                return EPS

        b.add_component(
            pname + "_x_ref", Expression(b.params.true_species_set, rule=rule_x_ref)
        )

    @staticmethod
    def ndIdn(b, pname, i):
        # Eqn 71
        return 0


class Symmetric(object):
    """
    Sub-methods for the symmetric (fused-salt) reference state
    """

    @staticmethod
    def ref_state(b, pname):
        def rule_x_ref(b, i):
            if i in b.params.ion_set:
                # Eqn 66
                return b.mole_frac_phase_comp_true[pname, i] / sum(
                    b.mole_frac_phase_comp_true[pname, j] for j in b.params.ion_set
                )
            else:
                return 0

        b.add_component(
            pname + "_x_ref", Expression(b.params.true_species_set, rule=rule_x_ref)
        )

    @staticmethod
    def ndIdn(b, pname, i):
        # Eqn 75
        return 0.5 * sum(
            cobj(b, j).config.charge ** 2 * ndxdn(b, pname, i, j)
            for j in b.params.ion_set
        )


def ndxdn(b, pname, i, j):
    x0 = getattr(b, pname + "_x_ref")

    # Delta function used in Eqns 73-76 (not defined in paper)
    if i == j:
        delta = 1
    else:
        delta = 0

    # Eqn 76
    return (delta - x0[j]) / sum(
        b.mole_frac_phase_comp_true[pname, k] for k in b.params.ion_set
    )
