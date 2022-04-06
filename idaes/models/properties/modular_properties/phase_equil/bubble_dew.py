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
from pyomo.environ import Constraint

from idaes.models.properties.modular_properties.base.utility import (
    get_method,
    get_component_object as cobj,
)
import idaes.core.util.scaling as iscale


class IdealBubbleDew:
    # -------------------------------------------------------------------------
    # Bubble temperature methods
    # This approach can only be used when both liquid and vapor phases use
    # Ideal properties
    # Henry's Law components also cause issues due to the need to (potentially)
    # calcualte concentrations at the bubble and dew points
    @staticmethod
    def temperature_bubble(b):
        try:

            def rule_bubble_temp(b, p1, p2):
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    l_only_comps,
                    v_only_comps,
                ) = _valid_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif v_only_comps != []:
                    # Non-condensables present, no bubble point
                    return Constraint.Skip

                return (
                    sum(
                        b.mole_frac_comp[j]
                        * get_method(b, "pressure_sat_comp", j)(
                            b, cobj(b, j), b.temperature_bubble[p1, p2]
                        )
                        for j in vl_comps
                    )
                    + sum(
                        b.mole_frac_comp[j]
                        * b.params.get_component(j)
                        .config.henry_component[l_phase]["method"]
                        .return_expression(b, l_phase, j, b.temperature_bubble[p1, p2])
                        for j in henry_comps
                    )
                    - b.pressure
                ) == 0

            b.eq_temperature_bubble = Constraint(
                b.params._pe_pairs, rule=rule_bubble_temp
            )
        except AttributeError:
            b.del_component(b.eq_temperature_bubble)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_bubble_temp(b, p1, p2, j):
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                l_only_comps,
                v_only_comps,
            ) = _valid_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif v_only_comps != []:
                # Non-condensables present, no bubble point
                return Constraint.Skip

            if j in vl_comps:
                return b._mole_frac_tbub[p1, p2, j] * b.pressure == (
                    b.mole_frac_comp[j]
                    * get_method(b, "pressure_sat_comp", j)(
                        b, cobj(b, j), b.temperature_bubble[p1, p2]
                    )
                )
            elif j in henry_comps:
                return b._mole_frac_tbub[p1, p2, j] * b.pressure == (
                    b.mole_frac_comp[j]
                    * b.params.get_component(j)
                    .config.henry_component[l_phase]["method"]
                    .return_expression(b, l_phase, j, b.temperature_bubble[p1, p2])
                )
            else:
                return b._mole_frac_tbub[p1, p2, j] == 0

        b.eq_mole_frac_tbub = Constraint(
            b.params._pe_pairs, b.component_list, rule=rule_mole_frac_bubble_temp
        )

    @staticmethod
    def scale_temperature_bubble(b, overwrite=True):
        sf_P = iscale.get_scaling_factor(b.pressure, default=1e-5, warning=True)
        sf_mf = iscale.get_scaling_factor(b.mole_frac_comp, default=1e3, warning=True)

        for pp in b.params._pe_pairs:
            for j in b.component_list:
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    l_only_comps,
                    v_only_comps,
                ) = _valid_VL_component_list(b, pp)
                if l_phase is None or v_phase is None:
                    continue
                elif v_only_comps != []:
                    continue

                if j in vl_comps:
                    sf = sf_P * sf_mf
                else:
                    sf = sf_mf

                iscale.constraint_scaling_transform(
                    b.eq_temperature_bubble[pp[0], pp[1]], sf_P, overwrite=overwrite
                )
                iscale.constraint_scaling_transform(
                    b.eq_mole_frac_tbub[pp[0], pp[1], j], sf, overwrite=overwrite
                )

    # -------------------------------------------------------------------------
    # Dew temperature methods
    @staticmethod
    def temperature_dew(b):
        try:

            def rule_dew_temp(b, p1, p2):
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    l_only_comps,
                    v_only_comps,
                ) = _valid_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif l_only_comps != []:
                    # Non-vaporisables present, no dew point
                    return Constraint.Skip

                return (
                    b.pressure
                    * (
                        sum(
                            b.mole_frac_comp[j]
                            / get_method(b, "pressure_sat_comp", j)(
                                b, cobj(b, j), b.temperature_dew[p1, p2]
                            )
                            for j in vl_comps
                        )
                        + sum(
                            b.mole_frac_comp[j]
                            / b.params.get_component(j)
                            .config.henry_component[l_phase]["method"]
                            .return_expression(b, l_phase, j, b.temperature_dew[p1, p2])
                            for j in henry_comps
                        )
                    )
                    - 1
                    == 0
                )

            b.eq_temperature_dew = Constraint(b.params._pe_pairs, rule=rule_dew_temp)
        except AttributeError:
            b.del_component(b.eq_temperature_dew)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_dew_temp(b, p1, p2, j):
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                l_only_comps,
                v_only_comps,
            ) = _valid_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif l_only_comps != []:
                # Non-vaporisables present, no dew point
                return Constraint.Skip

            if j in vl_comps:
                return (
                    b._mole_frac_tdew[p1, p2, j]
                    * get_method(b, "pressure_sat_comp", j)(
                        b, cobj(b, j), b.temperature_dew[p1, p2]
                    )
                    == b.mole_frac_comp[j] * b.pressure
                )
            elif j in henry_comps:
                return (
                    b._mole_frac_tdew[p1, p2, j]
                    * b.params.get_component(j)
                    .config.henry_component[l_phase]["method"]
                    .return_expression(b, l_phase, j, b.temperature_dew[p1, p2])
                    == b.mole_frac_comp[j] * b.pressure
                )
            else:
                return b._mole_frac_tdew[p1, p2, j] == 0

        b.eq_mole_frac_tdew = Constraint(
            b.params._pe_pairs, b.component_list, rule=rule_mole_frac_dew_temp
        )

    @staticmethod
    def scale_temperature_dew(b, overwrite=True):
        sf_P = iscale.get_scaling_factor(b.pressure, default=1e-5, warning=True)
        sf_mf = iscale.get_scaling_factor(b.mole_frac_comp, default=1e3, warning=True)

        for pp in b.params._pe_pairs:
            for j in b.component_list:
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    l_only_comps,
                    v_only_comps,
                ) = _valid_VL_component_list(b, pp)
                if l_phase is None or v_phase is None:
                    continue
                elif v_only_comps != []:
                    continue

                if j in vl_comps:
                    sf = sf_P * sf_mf
                else:
                    sf = sf_mf

                # b.eq_temperature_dew is well-scaled by default
                iscale.constraint_scaling_transform(
                    b.eq_mole_frac_tdew[pp[0], pp[1], j], sf, overwrite=overwrite
                )

    # -------------------------------------------------------------------------
    # Bubble pressure methods
    @staticmethod
    def pressure_bubble(b):
        try:

            def rule_bubble_press(b, p1, p2):
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    l_only_comps,
                    v_only_comps,
                ) = _valid_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif v_only_comps != []:
                    # Non-condensables present, no bubble point
                    return Constraint.Skip

                return b.pressure_bubble[p1, p2] == (
                    sum(b.mole_frac_comp[j] * b.pressure_sat_comp[j] for j in vl_comps)
                    + sum(
                        b.mole_frac_comp[j] * b.henry[l_phase, j] for j in henry_comps
                    )
                )

            b.eq_pressure_bubble = Constraint(
                b.params._pe_pairs, rule=rule_bubble_press
            )
        except AttributeError:
            b.del_component(b.eq_pressure_bubble)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_bubble_press(b, p1, p2, j):
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                l_only_comps,
                v_only_comps,
            ) = _valid_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif v_only_comps != []:
                # Non-condensables present, no bubble point
                return Constraint.Skip

            if j in vl_comps:
                return (
                    b._mole_frac_pbub[p1, p2, j] * b.pressure_bubble[p1, p2]
                    == b.mole_frac_comp[j] * b.pressure_sat_comp[j]
                )
            if j in henry_comps:
                return (
                    b._mole_frac_pbub[p1, p2, j] * b.pressure_bubble[p1, p2]
                    == b.mole_frac_comp[j] * b.henry[l_phase, j]
                )
            else:
                return b._mole_frac_pbub[p1, p2, j] == 0

        b.eq_mole_frac_pbub = Constraint(
            b.params._pe_pairs, b.component_list, rule=rule_mole_frac_bubble_press
        )

    @staticmethod
    def scale_pressure_bubble(b, overwrite=True):
        sf_P = iscale.get_scaling_factor(b.pressure, default=1e-5, warning=True)
        sf_mf = iscale.get_scaling_factor(b.mole_frac_comp, default=1e3, warning=True)

        for pp in b.params._pe_pairs:
            for j in b.component_list:
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    l_only_comps,
                    v_only_comps,
                ) = _valid_VL_component_list(b, pp)
                if l_phase is None or v_phase is None:
                    continue
                elif v_only_comps != []:
                    continue

                if j in vl_comps:
                    sf = sf_P * sf_mf
                else:
                    sf = sf_mf

                iscale.constraint_scaling_transform(
                    b.eq_pressure_bubble[pp[0], pp[1]], sf_P, overwrite=overwrite
                )
                iscale.constraint_scaling_transform(
                    b.eq_mole_frac_pbub[pp[0], pp[1], j], sf, overwrite=overwrite
                )

    # -------------------------------------------------------------------------
    # Dew pressure methods
    @staticmethod
    def pressure_dew(b):
        try:

            def rule_dew_press(b, p1, p2):
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    l_only_comps,
                    v_only_comps,
                ) = _valid_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif l_only_comps != []:
                    # Non-vaporisables present, no dew point
                    return Constraint.Skip

                return 0 == 1 - b.pressure_dew[p1, p2] * (
                    sum(b.mole_frac_comp[j] / b.pressure_sat_comp[j] for j in vl_comps)
                    + sum(
                        b.mole_frac_comp[j] / b.henry[l_phase, j] for j in henry_comps
                    )
                )

            b.eq_pressure_dew = Constraint(b.params._pe_pairs, rule=rule_dew_press)
        except AttributeError:
            b.del_component(b.eq_pressure_dew)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_dew_press(b, p1, p2, j):
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                l_only_comps,
                v_only_comps,
            ) = _valid_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif l_only_comps != []:
                # Non-vaporisables present, no dew point
                return Constraint.Skip

            if j in vl_comps:
                return (
                    b._mole_frac_pdew[p1, p2, j] * b.pressure_sat_comp[j]
                    == b.mole_frac_comp[j] * b.pressure_dew[p1, p2]
                )
            elif j in henry_comps:
                return (
                    b._mole_frac_pdew[p1, p2, j] * b.henry[l_phase, j]
                    == b.mole_frac_comp[j] * b.pressure_dew[p1, p2]
                )
            else:
                return b._mole_frac_pdew[p1, p2, j] == 0

        b.eq_mole_frac_pdew = Constraint(
            b.params._pe_pairs, b.component_list, rule=rule_mole_frac_dew_press
        )

    @staticmethod
    def scale_pressure_dew(b, overwrite=True):
        sf_P = iscale.get_scaling_factor(b.pressure, default=1e-5, warning=True)
        sf_mf = iscale.get_scaling_factor(b.mole_frac_comp, default=1e3, warning=True)

        for pp in b.params._pe_pairs:
            for j in b.component_list:
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    l_only_comps,
                    v_only_comps,
                ) = _valid_VL_component_list(b, pp)
                if l_phase is None or v_phase is None:
                    continue
                elif v_only_comps != []:
                    continue

                if j in vl_comps:
                    sf = sf_P * sf_mf
                else:
                    sf = sf_mf

                # b.eq_pressure_dew is well-scaled by default
                iscale.constraint_scaling_transform(
                    b.eq_mole_frac_pdew[pp[0], pp[1], j], sf, overwrite=overwrite
                )


class LogBubbleDew:
    # -------------------------------------------------------------------------
    # Bubble temperature methods
    @staticmethod
    def temperature_bubble(b):
        try:

            def rule_bubble_temp(b, p1, p2, j):
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    l_only_comps,
                    v_only_comps,
                ) = _valid_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif v_only_comps != []:
                    # Non-condensables present, no bubble point
                    return Constraint.Skip

                l_eos = b.params.get_phase(l_phase).config.equation_of_state
                v_eos = b.params.get_phase(v_phase).config.equation_of_state

                if j in vl_comps or j in henry_comps:
                    return l_eos.log_fug_phase_comp_Tbub(
                        b, l_phase, j, (p1, p2)
                    ) == v_eos.log_fug_phase_comp_Tbub(b, v_phase, j, (p1, p2))
                else:
                    return b._mole_frac_tbub[p1, p2, j] == 0

            b.eq_temperature_bubble = Constraint(
                b.params._pe_pairs, b.component_list, rule=rule_bubble_temp
            )
        except AttributeError:
            b.del_component(b.eq_temperature_bubble)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_bubble_temp(b, p1, p2):
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                l_only_comps,
                v_only_comps,
            ) = _valid_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif v_only_comps != []:
                # Non-condensables present, no bubble point
                return Constraint.Skip

            return 1 == (
                sum(b._mole_frac_tbub[p1, p2, j] for j in vl_comps)
                + sum(b._mole_frac_tbub[p1, p2, j] for j in henry_comps)
            )

        b.eq_mole_frac_tbub = Constraint(
            b.params._pe_pairs, rule=rule_mole_frac_bubble_temp
        )

    @staticmethod
    def scale_temperature_bubble(b, overwrite=True):
        sf_mf = iscale.get_scaling_factor(b.mole_frac_comp, default=1e3, warning=True)

        for pp in b.params._pe_pairs:
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                l_only_comps,
                v_only_comps,
            ) = _valid_VL_component_list(b, pp)
            if l_phase is None or v_phase is None:
                continue
            elif v_only_comps != []:
                continue

            # Assume b.eq_temperature_bubble is well-scaled
            iscale.constraint_scaling_transform(
                b.eq_mole_frac_tbub[pp[0], pp[1]], sf_mf, overwrite=overwrite
            )

    # -------------------------------------------------------------------------
    # Dew temperature methods
    @staticmethod
    def temperature_dew(b):
        try:

            def rule_dew_temp(b, p1, p2, j):
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    l_only_comps,
                    v_only_comps,
                ) = _valid_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif l_only_comps != []:
                    # Non-vapouriszbles present, no dew point
                    return Constraint.Skip

                l_eos = b.params.get_phase(l_phase).config.equation_of_state
                v_eos = b.params.get_phase(v_phase).config.equation_of_state

                if j in vl_comps or j in henry_comps:
                    return l_eos.log_fug_phase_comp_Tdew(
                        b, l_phase, j, (p1, p2)
                    ) == v_eos.log_fug_phase_comp_Tdew(b, v_phase, j, (p1, p2))
                else:
                    return b._mole_frac_tdew[p1, p2, j] == 0

            b.eq_temperature_dew = Constraint(
                b.params._pe_pairs, b.component_list, rule=rule_dew_temp
            )
        except AttributeError:
            b.del_component(b.eq_temperature_dew)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_dew_temp(b, p1, p2):
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                l_only_comps,
                v_only_comps,
            ) = _valid_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif l_only_comps != []:
                # Non-vaporisables present, no dew point
                return Constraint.Skip

            return 1 == (
                sum(b._mole_frac_tdew[p1, p2, j] for j in vl_comps)
                + sum(b._mole_frac_tdew[p1, p2, j] for j in henry_comps)
            )

        b.eq_mole_frac_tdew = Constraint(
            b.params._pe_pairs, rule=rule_mole_frac_dew_temp
        )

    @staticmethod
    def scale_temperature_dew(b, overwrite=True):
        sf_mf = iscale.get_scaling_factor(b.mole_frac_comp, default=1e3, warning=True)

        for pp in b.params._pe_pairs:
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                l_only_comps,
                v_only_comps,
            ) = _valid_VL_component_list(b, pp)
            if l_phase is None or v_phase is None:
                continue
            elif v_only_comps != []:
                continue

            # Assume b.eq_temperature_dew is well-scaled
            iscale.constraint_scaling_transform(
                b.eq_mole_frac_tdew[pp[0], pp[1]], sf_mf, overwrite=overwrite
            )

    # -------------------------------------------------------------------------
    # Bubble pressure methods
    @staticmethod
    def pressure_bubble(b):
        try:

            def rule_bubble_press(b, p1, p2, j):
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    l_only_comps,
                    v_only_comps,
                ) = _valid_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif v_only_comps != []:
                    # Non-condensables present, no bubble point
                    return Constraint.Skip

                l_eos = b.params.get_phase(l_phase).config.equation_of_state
                v_eos = b.params.get_phase(v_phase).config.equation_of_state

                if j in vl_comps or j in henry_comps:
                    return l_eos.log_fug_phase_comp_Pbub(
                        b, l_phase, j, (p1, p2)
                    ) == v_eos.log_fug_phase_comp_Pbub(b, v_phase, j, (p1, p2))
                else:
                    return b._mole_frac_pbub[p1, p2, j] == 0

            b.eq_pressure_bubble = Constraint(
                b.params._pe_pairs, b.component_list, rule=rule_bubble_press
            )
        except AttributeError:
            b.del_component(b.eq_pressure_bubble)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_bubble_press(b, p1, p2):
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                l_only_comps,
                v_only_comps,
            ) = _valid_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif v_only_comps != []:
                # Non-condensables present, no bubble point
                return Constraint.Skip

            return 1 == (
                sum(b._mole_frac_pbub[p1, p2, j] for j in vl_comps)
                + sum(b._mole_frac_pbub[p1, p2, j] for j in henry_comps)
            )

        b.eq_mole_frac_pbub = Constraint(
            b.params._pe_pairs, rule=rule_mole_frac_bubble_press
        )

    @staticmethod
    def scale_pressure_bubble(b, overwrite=True):
        sf_mf = iscale.get_scaling_factor(b.mole_frac_comp, default=1e3, warning=True)

        for pp in b.params._pe_pairs:
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                l_only_comps,
                v_only_comps,
            ) = _valid_VL_component_list(b, pp)
            if l_phase is None or v_phase is None:
                continue
            elif v_only_comps != []:
                continue

            # Assume b.eq_pressure_bubble is well-scaled
            iscale.constraint_scaling_transform(
                b.eq_mole_frac_pbub[pp[0], pp[1]], sf_mf, overwrite=overwrite
            )

    # -------------------------------------------------------------------------
    # Dew pressure methods
    @staticmethod
    def pressure_dew(b):
        try:

            def rule_dew_press(b, p1, p2, j):
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    l_only_comps,
                    v_only_comps,
                ) = _valid_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif v_only_comps != []:
                    # Non-condensables present, no bubble point
                    return Constraint.Skip

                l_eos = b.params.get_phase(l_phase).config.equation_of_state
                v_eos = b.params.get_phase(v_phase).config.equation_of_state

                if j in vl_comps or j in henry_comps:
                    return l_eos.log_fug_phase_comp_Pdew(
                        b, l_phase, j, (p1, p2)
                    ) == v_eos.log_fug_phase_comp_Pdew(b, v_phase, j, (p1, p2))
                else:
                    return b._mole_frac_pdew[p1, p2, j] == 0

            b.eq_pressure_dew = Constraint(
                b.params._pe_pairs, b.component_list, rule=rule_dew_press
            )
        except AttributeError:
            b.del_component(b.eq_pressure_dew)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_dew_press(b, p1, p2):
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                l_only_comps,
                v_only_comps,
            ) = _valid_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif l_only_comps != []:
                # Non-vaporisables present, no dew point
                return Constraint.Skip

            return 1 == (
                sum(b._mole_frac_pdew[p1, p2, j] for j in vl_comps)
                + sum(b._mole_frac_pdew[p1, p2, j] for j in henry_comps)
            )

        b.eq_mole_frac_pdew = Constraint(
            b.params._pe_pairs, rule=rule_mole_frac_dew_press
        )

    @staticmethod
    def scale_pressure_dew(b, overwrite=True):
        sf_mf = iscale.get_scaling_factor(b.mole_frac_comp, default=1e3, warning=True)

        for pp in b.params._pe_pairs:
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                l_only_comps,
                v_only_comps,
            ) = _valid_VL_component_list(b, pp)
            if l_phase is None or v_phase is None:
                continue
            elif v_only_comps != []:
                continue

            # Assume b.eq_pressure_dew is well-scaled
            iscale.constraint_scaling_transform(
                b.eq_mole_frac_pdew[pp[0], pp[1]], sf_mf, overwrite=overwrite
            )


def _valid_VL_component_list(blk, pp):
    vl_comps = []
    henry_comps = []
    l_only_comps = []
    v_only_comps = []

    pparams = blk.params
    l_phase = None
    v_phase = None
    if pparams.get_phase(pp[0]).is_liquid_phase():
        l_phase = pp[0]
    elif pparams.get_phase(pp[0]).is_vapor_phase():
        v_phase = pp[0]

    if pparams.get_phase(pp[1]).is_liquid_phase():
        l_phase = pp[1]
    elif pparams.get_phase(pp[1]).is_vapor_phase():
        v_phase = pp[1]

    # Only need to do this for V-L pairs, so check
    if l_phase is not None and v_phase is not None:
        for j in blk.params.component_list:
            if (l_phase, j) in blk.phase_component_set and (
                v_phase,
                j,
            ) in blk.phase_component_set:
                cobj = pparams.get_component(j)
                if cobj.config.henry_component is not None and (
                    pp[0] in cobj.config.henry_component
                    or pp[1] in cobj.config.henry_component
                ):
                    henry_comps.append(j)
                else:
                    vl_comps.append(j)
            elif (l_phase, j) in blk.phase_component_set:
                l_only_comps.append(j)
            elif (v_phase, j) in blk.phase_component_set:
                v_only_comps.append(j)

    return l_phase, v_phase, vl_comps, henry_comps, l_only_comps, v_only_comps
