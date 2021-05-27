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

from idaes.generic_models.properties.core.generic.utility import (
    get_component_object as cobj)
import idaes.logger as idaeslog


# Set up logger
_log = idaeslog.getLogger(__name__)


class Unsymmetric(object):
    """
    Sub-methods for the symmetric (fused-salt) reference state
    """
    @staticmethod
    def ref_state(b, pname):
        def rule_x_ref(b, i):
            if i in b.params.solvent_set or i in b.params.solute_set:
                # Eqn 66
                return (b.mole_frac_phase_comp_true[pname, i] /
                        sum(b.mole_frac_phase_comp_true[pname, j]
                            for j in b.params.solvent_set|b.params.solute_set))
            else:
                return 0

        b.add_component(pname+"_x_ref",
                        Expression(b.params.true_species_set, rule=rule_x_ref))

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
                return (b.mole_frac_phase_comp_true[pname, i] /
                        sum(b.mole_frac_phase_comp_true[pname, j]
                            for j in b.params.ion_set))
            else:
                return 0

        b.add_component(pname+"_x_ref",
                        Expression(b.params.true_species_set, rule=rule_x_ref))

    @staticmethod
    def ndIdn(b, pname, i):
        # Eqn 75
        return 0.5*sum(cobj(b, j).config.charge**2*ndxdn(b, pname, i, j)
                       for j in b.params.ion_set)


def ndxdn(b, pname, i, j):
    x0 = getattr(b, pname+"_x_ref")

    # Delta function used in Eqns 73-76 (not defined in paper)
    if i == j:
        delta = 1
    else:
        delta = 0

    # Eqn 76
    return ((delta - x0[j]) /
            sum(b.mole_frac_phase_comp_true[pname, k]
                for k in b.params.ion_set))
