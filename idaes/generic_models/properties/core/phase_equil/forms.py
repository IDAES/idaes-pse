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
Library of common forms for phase equilibrium constraints
"""


def fugacity(b, phase1, phase2, comp):
    pp = (phase1, phase2)
    return (b.params.get_phase(phase1)
            .config.equation_of_state.fug_phase_comp_eq(b, phase1, comp, pp) ==
            b.params.get_phase(phase2)
            .config.equation_of_state.fug_phase_comp_eq(b, phase2, comp, pp))


def log_fugacity(b, phase1, phase2, comp):
    pp = (phase1, phase2)
    return (b.params.get_phase(phase1).config.equation_of_state
            .log_fug_phase_comp_eq(b, phase1, comp, pp) ==
            b.params.get_phase(phase2).config.equation_of_state
            .log_fug_phase_comp_eq(b, phase2, comp, pp))
