###############################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES), and is copyright (c) 2018-2021 by the
# software owners: The Regents of the University of California, through Lawrence
# Berkeley National Laboratory,  National Technology & Engineering Solutions of
# Sandia, LLC, Carnegie Mellon University, West Virginia University Research
# Corporation, et al.  All rights reserved.
#
# NOTICE.  This Software was developed under funding from the U.S. Department of
# Energy and the U.S. Government consequently retains certain rights. As such, the
# U.S. Government has been granted for itself and others acting on its behalf a
# paid-up, nonexclusive, irrevocable, worldwide license in the Software to
# reproduce, distribute copies to the public, prepare derivative works, and
# perform publicly and display publicly, and to permit other to do so.
#
###############################################################################
"""
Library of common forms for phase equilibrium constraints
"""
import idaes.core.util.scaling as iscale


class fugacity():

    @staticmethod
    def return_expression(b, phase1, phase2, comp):
        pp = (phase1, phase2)
        return (b.params.get_phase(phase1).config.equation_of_state
                .fug_phase_comp_eq(b, phase1, comp, pp) ==
                b.params.get_phase(phase2).config.equation_of_state
                .fug_phase_comp_eq(b, phase2, comp, pp))

    @staticmethod
    def calculate_scaling_factors(b, phase1, phase2, comp):
        with b.lock_attribute_creation_context():
            if hasattr(b, "fug_phase_comp_eq"):
                sf_1 = iscale.get_scaling_factor(
                    b.fug_phase_comp[phase1, comp], default=1, warning=True)
                sf_2 = iscale.get_scaling_factor(
                    b.fug_phase_comp[phase2, comp], default=1, warning=True)

                return min(sf_1, sf_2)
            elif hasattr(b, "fug_phase_comp"):
                sf_1 = iscale.get_scaling_factor(
                    b.fug_phase_comp[phase1, comp], default=1, warning=True)
                sf_2 = iscale.get_scaling_factor(
                    b.fug_phase_comp[phase2, comp], default=1, warning=True)

                return min(sf_1, sf_2)
            else:
                return 1


class log_fugacity():

    @staticmethod
    def return_expression(b, phase1, phase2, comp):
        pp = (phase1, phase2)
        return (b.params.get_phase(phase1).config.equation_of_state
                .log_fug_phase_comp_eq(b, phase1, comp, pp) ==
                b.params.get_phase(phase2).config.equation_of_state
                .log_fug_phase_comp_eq(b, phase2, comp, pp))

    @staticmethod
    def calculate_scaling_factors(b, phase1, phase2, comp):
        # For log fugacity, assume already well scaled
        return 1
