#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Library of common forms for phase equilibrium constraints
"""
# TODO: Missing docstrings
# pylint: disable=missing-function-docstring

import idaes.core.util.scaling as iscale
from idaes.core.scaling import CustomScalerBase


class FugacityScaler(CustomScalerBase):
    """
    Scaling method for the fugacity form of phase equilibrium
    """

    def variable_scaling_routine(self, model, index, overwrite: bool = False):
        # No variables added
        pass

    def constraint_scaling_routine(self, model, index, overwrite: bool = False):
        p1, p2, j = index
        self.scale_constraint_by_nominal_value(
            model.equilibrium_constraint[p1, p2, j], overwrite=overwrite
        )


class fugacity:
    """Phase equilibrium through equating fugacity"""

    default_scaler = FugacityScaler

    @staticmethod
    def return_expression(b, phase1, phase2, comp):
        pp = (phase1, phase2)
        return b.params.get_phase(phase1).config.equation_of_state.fug_phase_comp_eq(
            b, phase1, comp, pp
        ) == b.params.get_phase(phase2).config.equation_of_state.fug_phase_comp_eq(
            b, phase2, comp, pp
        )

    @staticmethod
    def calculate_scaling_factors(b, phase1, phase2, comp):
        sf_xp1 = iscale.get_scaling_factor(
            b.mole_frac_phase_comp[phase1, comp],
            default=1e3,
            warning=True,
        )
        sf_xp2 = iscale.get_scaling_factor(
            b.mole_frac_phase_comp[phase2, comp],
            default=1e3,
            warning=True,
        )
        sf_x = min(sf_xp1, sf_xp2)
        sf_P = iscale.get_scaling_factor(b.pressure, default=1e-5, warning=True)
        return sf_x * sf_P


class LogFugacityScaler(CustomScalerBase):
    """
    Scaling method for the logfugacity form of phase equilibrium
    """

    def variable_scaling_routine(self, model, index, overwrite: bool = False):
        # No variables added
        pass

    def constraint_scaling_routine(self, model, index, overwrite: bool = False):
        p1, p2, j = index
        if (p1, j) in model.phase_component_set and (
            p2,
            j,
        ) in model.phase_component_set:
            self.set_component_scaling_factor(
                model.equilibrium_constraint[p1, p2, j], 1, overwrite=overwrite
            )


class log_fugacity:
    """Phase equilibrium through equating log of fugacity."""

    default_scaler = LogFugacityScaler

    @staticmethod
    def return_expression(b, phase1, phase2, comp):
        pp = (phase1, phase2)
        return b.params.get_phase(
            phase1
        ).config.equation_of_state.log_fug_phase_comp_eq(
            b, phase1, comp, pp
        ) == b.params.get_phase(
            phase2
        ).config.equation_of_state.log_fug_phase_comp_eq(
            b, phase2, comp, pp
        )

    @staticmethod
    def calculate_scaling_factors(b, phase1, phase2, comp):
        # For log fugacity, assume already well scaled
        return 1
