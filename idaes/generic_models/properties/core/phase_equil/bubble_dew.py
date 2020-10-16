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
from pyomo.environ import Constraint, log

from idaes.generic_models.properties.core.generic.utility import \
        get_method, get_component_object as cobj


class IdealBubbleDew():
    # -------------------------------------------------------------------------
    # Bubble temperature methods
    @staticmethod
    def temperature_bubble(b):
        try:
            def rule_bubble_temp(b, p1, p2):
                (l_phase,
                 v_phase,
                 vl_comps,
                 l_only_comps,
                 v_only_comps) = _valid_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif v_only_comps != []:
                    # Non-condensables present, no bubble point
                    return Constraint.Skip

                return (sum(b.mole_frac_comp[j] *
                            get_method(b, "pressure_sat_comp", j)(
                                b, cobj(b, j), b.temperature_bubble[p1, p2])
                            for j in vl_comps) -
                        b.pressure) == 0
            b.eq_temperature_bubble = Constraint(b.params._pe_pairs,
                                                 rule=rule_bubble_temp)
        except AttributeError:
            b.del_component(b.eq_temperature_bubble)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_bubble_temp(b, p1, p2, j):
            (l_phase,
             v_phase,
             vl_comps,
             l_only_comps,
             v_only_comps) = _valid_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif v_only_comps != []:
                # Non-condensables present, no bubble point
                return Constraint.Skip

            if j in vl_comps:
                return b._mole_frac_tbub[p1, p2, j]*b.pressure == (
                    b.mole_frac_comp[j] *
                    get_method(b, "pressure_sat_comp", j)(
                        b, cobj(b, j), b.temperature_bubble[p1, p2]))
            else:
                return b._mole_frac_tbub[p1, p2, j] == 0
        b.eq_mole_frac_tbub = Constraint(b.params._pe_pairs,
                                         b.params.component_list,
                                         rule=rule_mole_frac_bubble_temp)

    # -------------------------------------------------------------------------
    # Dew temperature methods
    @staticmethod
    def temperature_dew(b):
        try:
            def rule_dew_temp(b, p1, p2):
                (l_phase,
                 v_phase,
                 vl_comps,
                 l_only_comps,
                 v_only_comps) = _valid_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif l_only_comps != []:
                    # Non-vaporisables present, no dew point
                    return Constraint.Skip

                return (b.pressure*sum(
                            b.mole_frac_comp[j] /
                            get_method(b, "pressure_sat_comp", j)(
                                b, cobj(b, j), b.temperature_dew[p1, p2])
                            for j in vl_comps) - 1 ==
                        0)
            b.eq_temperature_dew = Constraint(b.params._pe_pairs,
                                              rule=rule_dew_temp)
        except AttributeError:
            b.del_component(b.eq_temperature_dew)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_dew_temp(b, p1, p2, j):
            (l_phase,
             v_phase,
             vl_comps,
             l_only_comps,
             v_only_comps) = _valid_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif l_only_comps != []:
                # Non-vaporisables present, no dew point
                return Constraint.Skip

            if j in vl_comps:
                return (b._mole_frac_tdew[p1, p2, j] *
                        get_method(b, "pressure_sat_comp", j)(
                            b, cobj(b, j), b.temperature_dew[p1, p2]) ==
                        b.mole_frac_comp[j]*b.pressure)
            else:
                return b._mole_frac_tdew[p1, p2, j] == 0

        b.eq_mole_frac_tdew = Constraint(b.params._pe_pairs,
                                         b.params.component_list,
                                         rule=rule_mole_frac_dew_temp)

    # -------------------------------------------------------------------------
    # Bubble pressure methods
    @staticmethod
    def pressure_bubble(b):
        try:
            def rule_bubble_press(b, p1, p2):
                (l_phase,
                 v_phase,
                 vl_comps,
                 l_only_comps,
                 v_only_comps) = _valid_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif v_only_comps != []:
                    # Non-condensables present, no bubble point
                    return Constraint.Skip

                return b.pressure_bubble[p1, p2] == sum(
                        b.mole_frac_comp[j] *
                        get_method(b, "pressure_sat_comp", j)(
                            b, cobj(b, j), b.temperature)
                        for j in vl_comps)
            b.eq_pressure_bubble = Constraint(b.params._pe_pairs,
                                              rule=rule_bubble_press)
        except AttributeError:
            b.del_component(b.eq_pressure_bubble)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_bubble_press(b, p1, p2, j):
            (l_phase,
             v_phase,
             vl_comps,
             l_only_comps,
             v_only_comps) = _valid_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif v_only_comps != []:
                # Non-condensables present, no bubble point
                return Constraint.Skip

            if j in vl_comps:
                return (b._mole_frac_pbub[p1, p2, j] *
                        b.pressure_bubble[p1, p2]) == (
                    b.mole_frac_comp[j] *
                    get_method(b, "pressure_sat_comp", j)(
                        b, cobj(b, j), b.temperature))
            else:
                return b._mole_frac_pbub[p1, p2, j] == 0

        b.eq_mole_frac_pbub = Constraint(b.params._pe_pairs,
                                         b.params.component_list,
                                         rule=rule_mole_frac_bubble_press)

    # -------------------------------------------------------------------------
    # Dew pressure methods
    @staticmethod
    def pressure_dew(b):
        try:
            def rule_dew_press(b, p1, p2):
                (l_phase,
                 v_phase,
                 vl_comps,
                 l_only_comps,
                 v_only_comps) = _valid_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif l_only_comps != []:
                    # Non-vaporisables present, no dew point
                    return Constraint.Skip

                return 0 == 1 - b.pressure_dew[p1, p2]*sum(
                        b.mole_frac_comp[j] /
                        get_method(b, "pressure_sat_comp", j)(
                            b, cobj(b, j), b.temperature)
                        for j in vl_comps)
            b.eq_pressure_dew = Constraint(b.params._pe_pairs,
                                           rule=rule_dew_press)
        except AttributeError:
            b.del_component(b.eq_pressure_dew)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_dew_press(b, p1, p2, j):
            (l_phase,
             v_phase,
             vl_comps,
             l_only_comps,
             v_only_comps) = _valid_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif l_only_comps != []:
                # Non-vaporisables present, no dew point
                return Constraint.Skip

            if j in vl_comps:
                return (b._mole_frac_pdew[p1, p2, j] *
                        get_method(b, "pressure_sat_comp", j)(
                            b, cobj(b, j), b.temperature) ==
                        b.mole_frac_comp[j]*b.pressure_dew[p1, p2])
            else:
                return b._mole_frac_pdew[p1, p2, j] == 0

        b.eq_mole_frac_pdew = Constraint(b.params._pe_pairs,
                                         b.params.component_list,
                                         rule=rule_mole_frac_dew_press)


class LogBubbleDew():
    # -------------------------------------------------------------------------
    # Bubble temperature methods
    @staticmethod
    def temperature_bubble(b):
        try:
            def rule_bubble_temp(b, p1, p2, j):
                (l_phase,
                 v_phase,
                 vl_comps,
                 l_only_comps,
                 v_only_comps) = _valid_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif v_only_comps != []:
                    # Non-condensables present, no bubble point
                    return Constraint.Skip

                l_eos = b.params.get_phase(l_phase).config.equation_of_state
                v_eos = b.params.get_phase(v_phase).config.equation_of_state

                if j in vl_comps:
                    return (
                        log(b.mole_frac_comp[j]) +
                        l_eos.log_fug_coeff_phase_comp_Tbub(
                            b, l_phase, j, (p1, p2)) ==
                        log(b._mole_frac_tbub[p1, p2, j]) +
                        v_eos.log_fug_coeff_phase_comp_Tbub(
                            b, v_phase, j, (p1, p2)))
                else:
                    return b._mole_frac_tbub[p1, p2, j] == 0

            b.eq_temperature_bubble = Constraint(b.params._pe_pairs,
                                                 b.params.component_list,
                                                 rule=rule_bubble_temp)
        except AttributeError:
            b.del_component(b.eq_temperature_bubble)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_bubble_temp(b, p1, p2):
            (l_phase,
             v_phase,
             vl_comps,
             l_only_comps,
             v_only_comps) = _valid_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif v_only_comps != []:
                # Non-condensables present, no bubble point
                return Constraint.Skip

            return 1e3 == 1e3*sum(b._mole_frac_tbub[p1, p2, j]
                                  for j in vl_comps)
        b.eq_mole_frac_tbub = Constraint(b.params._pe_pairs,
                                         rule=rule_mole_frac_bubble_temp)

    # -------------------------------------------------------------------------
    # Dew temperature methods
    @staticmethod
    def temperature_dew(b):
        try:
            def rule_dew_temp(b, p1, p2, j):
                (l_phase,
                 v_phase,
                 vl_comps,
                 l_only_comps,
                 v_only_comps) = _valid_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif l_only_comps != []:
                    # Non-vapouriszbles present, no dew point
                    return Constraint.Skip

                l_eos = b.params.get_phase(l_phase).config.equation_of_state
                v_eos = b.params.get_phase(v_phase).config.equation_of_state

                if j in vl_comps:
                    return (
                        log(b._mole_frac_tdew[p1, p2, j]) +
                        l_eos.log_fug_coeff_phase_comp_Tdew(
                            b, l_phase, j, (p1, p2)) ==
                        log(b.mole_frac_comp[j]) +
                        v_eos.log_fug_coeff_phase_comp_Tdew(
                            b, v_phase, j, (p1, p2)))
                else:
                    return b._mole_frac_tdew[p1, p2, j] == 0
            b.eq_temperature_dew = Constraint(b.params._pe_pairs,
                                              b.params.component_list,
                                              rule=rule_dew_temp)
        except AttributeError:
            b.del_component(b.eq_temperature_dew)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_dew_temp(b, p1, p2):
            (l_phase,
             v_phase,
             vl_comps,
             l_only_comps,
             v_only_comps) = _valid_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif l_only_comps != []:
                # Non-vaporisables present, no dew point
                return Constraint.Skip

            return 1e3 == 1e3*sum(b._mole_frac_tdew[p1, p2, j]
                                  for j in vl_comps)
        b.eq_mole_frac_tdew = Constraint(b.params._pe_pairs,
                                         rule=rule_mole_frac_dew_temp)

    # -------------------------------------------------------------------------
    # Bubble pressure methods
    @staticmethod
    def pressure_bubble(b):
        try:
            def rule_bubble_press(b, p1, p2, j):
                (l_phase,
                 v_phase,
                 vl_comps,
                 l_only_comps,
                 v_only_comps) = _valid_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif v_only_comps != []:
                    # Non-condensables present, no bubble point
                    return Constraint.Skip

                l_eos = b.params.get_phase(l_phase).config.equation_of_state
                v_eos = b.params.get_phase(v_phase).config.equation_of_state

                if j in vl_comps:
                    return (
                        log(b.mole_frac_comp[j]) +
                        l_eos.log_fug_coeff_phase_comp_Pbub(
                            b, l_phase, j, (p1, p2)) ==
                        log(b._mole_frac_pbub[p1, p2, j]) +
                        v_eos.log_fug_coeff_phase_comp_Pbub(
                            b, v_phase, j, (p1, p2)))
                else:
                    return b._mole_frac_pbub[p1, p2, j] == 0
            b.eq_pressure_bubble = Constraint(b.params._pe_pairs,
                                              b.params.component_list,
                                              rule=rule_bubble_press)
        except AttributeError:
            b.del_component(b.eq_pressure_bubble)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_bubble_press(b, p1, p2):
            (l_phase,
             v_phase,
             vl_comps,
             l_only_comps,
             v_only_comps) = _valid_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif v_only_comps != []:
                # Non-condensables present, no bubble point
                return Constraint.Skip

            return 1e3 == 1e3*sum(b._mole_frac_pbub[p1, p2, j]
                                  for j in vl_comps)
        b.eq_mole_frac_pbub = Constraint(b.params._pe_pairs,
                                         rule=rule_mole_frac_bubble_press)

    # -------------------------------------------------------------------------
    # Dew pressure methods
    @staticmethod
    def pressure_dew(b):
        try:
            def rule_dew_press(b, p1, p2, j):
                (l_phase,
                 v_phase,
                 vl_comps,
                 l_only_comps,
                 v_only_comps) = _valid_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif v_only_comps != []:
                    # Non-condensables present, no bubble point
                    return Constraint.Skip

                l_eos = b.params.get_phase(l_phase).config.equation_of_state
                v_eos = b.params.get_phase(v_phase).config.equation_of_state

                if j in vl_comps:
                    return (
                        log(b._mole_frac_p_dew[p1, p2, j]) +
                        l_eos.log_fug_coeff_phase_comp_Pdew(
                            b, l_phase, j, (p1, p2)) ==
                        log(b.mole_frac_comp[j]) +
                        v_eos.log_fug_coeff_phase_comp_Pdew(
                            b, v_phase, j, (p1, p2)))
                else:
                    return b._mole_frac_pdew[p1, p2, j] == 0
            b.eq_pressure_dew = Constraint(b.params._pe_pairs,
                                           b.params.component_list,
                                           rule=rule_dew_press)
        except AttributeError:
            b.del_component(b.eq_pressure_dew)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_dew_press(b, p1, p2):
            (l_phase,
             v_phase,
             vl_comps,
             l_only_comps,
             v_only_comps) = _valid_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif l_only_comps != []:
                # Non-vaporisables present, no dew point
                return Constraint.Skip

            return 1e3 == 1e3*sum(b._mole_frac_pdew[p1, p2, j]
                                  for j in vl_comps)
        b.eq_mole_frac_pdew = Constraint(b.params._pe_pairs,
                                         rule=rule_mole_frac_dew_press)


def _valid_VL_component_list(blk, pp):
    vl_comps = []
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
            if ((l_phase, j) in pparams._phase_component_set and
                    (v_phase, j) in pparams._phase_component_set):
                vl_comps.append(j)
            elif (l_phase, j) in pparams._phase_component_set:
                l_only_comps.append(j)
            elif (v_phase, j) in pparams._phase_component_set:
                v_only_comps.append(j)

    return l_phase, v_phase, vl_comps, l_only_comps, v_only_comps
