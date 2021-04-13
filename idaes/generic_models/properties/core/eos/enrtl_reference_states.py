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

import idaes.logger as idaeslog


# Set up logger
_log = idaeslog.getLogger(__name__)


class Symmetric(object):
    """
    Sub-methods for the symmetric (fused-salt) reference state
    """
    def ref_state(b, pname):
        def rule_x_ref(b, i):
            if i in b.params.ion_set:
                return (b.mole_frac_phase_comp_true[pname, i] /
                        sum(b.mole_frac_phase_comp_true[pname, j]
                            for j in b.params.ion_set))
            else:
                return 0

        b.add_component(pname+"_x_ref",
                        Expression(b.params.true_species_set, rule=rule_x_ref))

    def delta(i, j):
        # Delta function used in Eqns 73-76 (not defined in paper)
        if i == j:
            return 1
        else:
            return 0

    def ndIdn(b, i, j):
        pass
