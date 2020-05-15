##############################################################################
# Institute for the Design of Advanced Energy Systems Process Systems
# Engineering Framework (IDAES PSE Framework) Copyright (c) 2018-2019, by the
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

from idaes.generic_models.properties.core.generic.generic_property import \
        get_method, get_component_object as cobj


def _vle_check(b, p1, p2):
    # Bubble nad dew point only make sense for a vapor-liquid pair
    # Check if this is the case
    p1_obj = b.params.get_phase(p1)
    p2_obj = b.params.get_phase(p2)
    # One phase must be liquid nad the other vapor
    if not (p1_obj.is_liquid_phase() or p2_obj.is_liquid_phase()):
        return False
    elif not (p1_obj.is_vapor_phase() or p2_obj.is_vapor_phase()):
        return False
    return True


def _identify_phases(b, p1, p2):
    v_eos = None
    l_eos = None

    p1_obj = b.params.get_phase(p1)
    p2_obj = b.params.get_phase(p2)
    # Identify liquid and vapor phases
    if p1_obj.is_liquid_phase():
        l_phase = p1
        l_eos = p1_obj.config.equation_of_state
    elif p2_obj.is_liquid_phase():
        l_phase = p2
        l_eos = p2_obj.config.equation_of_state

    if p1_obj.is_vapor_phase():
        v_phase = p1
        v_eos = p1_obj.config.equation_of_state
    elif p2_obj.is_vapor_phase():
        v_phase = p2
        v_eos = p2_obj.config.equation_of_state

    return l_phase, l_eos, v_phase, v_eos


class IdealBubbleDew():
    # -------------------------------------------------------------------------
    # Bubble temperature methods
    def temperature_bubble(b):
        try:
            def rule_bubble_temp(b, p1, p2):
                if not _vle_check(b, p1, p2):
                    return Constraint.Skip

                return (sum(b.mole_frac_comp[j] *
                            get_method(b, "pressure_sat_comp", j)(
                                b, cobj(b, j), b.temperature_bubble[p1, p2])
                            for j in b.params.component_list) -
                        b.pressure) == 0
            b.eq_temperature_bubble = Constraint(b.params._pe_pairs,
                                                 rule=rule_bubble_temp)
        except AttributeError:
            b.del_component(b.eq_temperature_bubble)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_bubble_temp(b, p1, p2, j):
            if not _vle_check(b, p1, p2):
                return Constraint.Skip
            return b._mole_frac_tbub[p1, p2, j]*b.pressure == (
                b.mole_frac_comp[j] *
                get_method(b, "pressure_sat_comp", j)(
                    b, cobj(b, j), b.temperature_bubble[p1, p2]))
        b.eq_mole_frac_tbub = Constraint(b.params._pe_pairs,
                                         b.params.component_list,
                                         rule=rule_mole_frac_bubble_temp)

    # -------------------------------------------------------------------------
    # Dew temperature methods
    def temperature_dew(b):
        try:
            def rule_dew_temp(b, p1, p2):
                if not _vle_check(b, p1, p2):
                    return Constraint.Skip
                return (b.pressure*sum(
                            b.mole_frac_comp[j] /
                            get_method(b, "pressure_sat_comp", j)(
                                b, cobj(b, j), b.temperature_dew[p1, p2])
                            for j in b.params.component_list) - 1 ==
                        0)
            b.eq_temperature_dew = Constraint(b.params._pe_pairs,
                                              rule=rule_dew_temp)
        except AttributeError:
            b.del_component(b.eq_temperature_dew)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_dew_temp(b, p1, p2, j):
            if not _vle_check(b, p1, p2):
                return Constraint.Skip
            return (b._mole_frac_tdew[p1, p2, j] *
                    get_method(b, "pressure_sat_comp", j)(
                        b, cobj(b, j), b.temperature_dew[p1, p2]) ==
                    b.mole_frac_comp[j]*b.pressure)
        b.eq_mole_frac_tdew = Constraint(b.params._pe_pairs,
                                         b.params.component_list,
                                         rule=rule_mole_frac_dew_temp)

    # -------------------------------------------------------------------------
    # Bubble pressure methods
    def pressure_bubble(b):
        try:
            def rule_bubble_press(b, p1, p2):
                if not _vle_check(b, p1, p2):
                    return Constraint.Skip
                return b.pressure_bubble[p1, p2] == sum(
                        b.mole_frac_comp[j] *
                        get_method(b, "pressure_sat_comp", j)(
                            b, cobj(b, j), b.temperature)
                        for j in b.params.component_list)
            b.eq_pressure_bubble = Constraint(b.params._pe_pairs,
                                              rule=rule_bubble_press)
        except AttributeError:
            b.del_component(b.eq_pressure_bubble)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_bubble_press(b, p1, p2, j):
            if not _vle_check(b, p1, p2):
                return Constraint.Skip
            return b._mole_frac_pbub[p1, p2, j]*b.pressure_bubble[p1, p2] == (
                b.mole_frac_comp[j] *
                get_method(b, "pressure_sat_comp", j)(
                    b, cobj(b, j), b.temperature))
        b.eq_mole_frac_pbub = Constraint(b.params._pe_pairs,
                                         b.params.component_list,
                                         rule=rule_mole_frac_bubble_press)

    # -------------------------------------------------------------------------
    # Dew pressure methods
    def pressure_dew(b):
        try:
            def rule_dew_press(b, p1, p2):
                if not _vle_check(b, p1, p2):
                    return Constraint.Skip
                return 0 == 1 - b.pressure_dew[p1, p2]*sum(
                        b.mole_frac_comp[j] /
                        get_method(b, "pressure_sat_comp", j)(
                            b, cobj(b, j), b.temperature)
                        for j in b.params.component_list)
            b.eq_pressure_dew = Constraint(b.params._pe_pairs,
                                           rule=rule_dew_press)
        except AttributeError:
            b.del_component(b.eq_pressure_dew)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_dew_press(b, p1, p2, j):
            if not _vle_check(b, p1, p2):
                return Constraint.Skip
            return (b._mole_frac_pdew[p1, p2, j] *
                    get_method(b, "pressure_sat_comp", j)(
                        b, cobj(b, j), b.temperature) ==
                    b.mole_frac_comp[j]*b.pressure_dew[p1, p2])
        b.eq_mole_frac_pdew = Constraint(b.params._pe_pairs,
                                         b.params.component_list,
                                         rule=rule_mole_frac_dew_press)


class LogBubbleDew():
    # -------------------------------------------------------------------------
    # Bubble temperature methods
    def temperature_bubble(b):
        try:
            def rule_bubble_temp(b, p1, p2, j):
                l_phase, l_eos, v_phase, v_eos = _identify_phases(b, p1, p2)

                # If one or both of v_phase and l_phase is None, this is not
                # a vapor-liquid pair, so sjip constraint
                if v_eos is None and l_eos is None:
                    return Constraint.Skip

                return (
                    log(b.mole_frac_comp[j]) +
                    log(l_eos.fug_phase_comp_Tbub(b, l_phase, j, (p1, p2))) ==
                    log(b._mole_frac_tbub[p1, p2, j]) +
                    log(v_eos.fug_phase_comp_Tbub(b, v_phase, j, (p1, p2))))
            b.eq_temperature_bubble = Constraint(b.params._pe_pairs,
                                                 b.params.component_list,
                                                 rule=rule_bubble_temp)
        except AttributeError:
            b.del_component(b.eq_temperature_bubble)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_bubble_temp(b, p1, p2):
            if not _vle_check(b, p1, p2):
                return Constraint.Skip
            return 1e3 == 1e3*sum(b._mole_frac_tbub[p1, p2, j]
                                  for j in b.params.component_list)
        b.eq_mole_frac_tbub = Constraint(b.params._pe_pairs,
                                         rule=rule_mole_frac_bubble_temp)

    # -------------------------------------------------------------------------
    # Dew temperature methods
    def temperature_dew(b):
        try:
            def rule_dew_temp(b, p1, p2, j):
                l_phase, l_eos, v_phase, v_eos = _identify_phases(b, p1, p2)

                # If one or both of v_phase and l_phase is None, this is not
                # a vapor-liquid pair, so sjip constraint
                if v_eos is None and l_eos is None:
                    return Constraint.Skip

                return (
                    log(b._mole_frac_tdew[p1, p2, j]) +
                    log(l_eos.fug_phase_comp_Tdew(b, l_phase, j, (p1, p2))) ==
                    log(b.mole_frac_comp[j]) +
                    log(v_eos.fug_phase_comp_Tdew(b, v_phase, j, (p1, p2))))
            b.eq_temperature_dew = Constraint(b.params._pe_pairs,
                                              b.params.component_list,
                                              rule=rule_dew_temp)
        except AttributeError:
            b.del_component(b.eq_temperature_dew)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_dew_temp(b, p1, p2):
            l_phase, l_eos, v_phase, v_eos = _identify_phases(b, p1, p2)

            # If one or both of v_phase and l_phase is None, this is not
            # a vapor-liquid pair, so sjip constraint
            if v_eos is None and l_eos is None:
                return Constraint.Skip

            return 1e3 == 1e3*sum(b._mole_frac_tdew[p1, p2, j]
                                  for j in b.params.component_list)
        b.eq_mole_frac_tdew = Constraint(b.params._pe_pairs,
                                         rule=rule_mole_frac_dew_temp)

    # -------------------------------------------------------------------------
    # Bubble pressure methods
    def pressure_bubble(b):
        try:
            def rule_bubble_press(b, p1, p2, j):
                l_phase, l_eos, v_phase, v_eos = _identify_phases(b, p1, p2)

                # If one or both of v_phase and l_phase is None, this is not
                # a vapor-liquid pair, so sjip constraint
                if v_eos is None and l_eos is None:
                    return Constraint.Skip

                return (
                    log(b.mole_frac_comp[j]) +
                    log(l_eos.fug_phase_comp_Pbub(b, l_phase, j, (p1, p2))) ==
                    log(b._mole_frac_pbub[p1, p2, j]) +
                    log(v_eos.fug_phase_comp_Pbub(b, v_phase, j, (p1, p2))))
            b.eq_pressure_bubble = Constraint(b.params._pe_pairs,
                                              b.params.component_list,
                                              rule=rule_bubble_press)
        except AttributeError:
            b.del_component(b.eq_pressure_bubble)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_bubble_press(b, p1, p2):
            l_phase, l_eos, v_phase, v_eos = _identify_phases(b, p1, p2)

            # If one or both of v_phase and l_phase is None, this is not
            # a vapor-liquid pair, so sjip constraint
            if v_eos is None and l_eos is None:
                return Constraint.Skip

            return 1e3 == 1e3*sum(b._mole_frac_pbub[p1, p2, j]
                                  for j in b.params.component_list)
        b.eq_mole_frac_pbub = Constraint(b.params._pe_pairs,
                                         rule=rule_mole_frac_bubble_press)

    # -------------------------------------------------------------------------
    # Dew pressure methods
    def pressure_dew(b):
        try:
            def rule_dew_press(b, p1, p2, j):
                l_phase, l_eos, v_phase, v_eos = _identify_phases(b, p1, p2)

                # If one or both of v_phase and l_phase is None, this is not
                # a vapor-liquid pair, so sjip constraint
                if v_eos is None and l_eos is None:
                    return Constraint.Skip

                return (
                    log(b._mole_frac_p_dew[p1, p2, j]) +
                    log(l_eos.fug_phase_comp_Pdew(b, l_phase, j, (p1, p2))) ==
                    log(b.mole_frac_comp[j]) +
                    log(v_eos.fug_phase_comp_Pdew(b, v_phase, j, (p1, p2))))
            b.eq_pressure_dew = Constraint(b.params._pe_pairs,
                                           b.params.component_list,
                                           rule=rule_dew_press)
        except AttributeError:
            b.del_component(b.eq_pressure_dew)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_dew_press(b, p1, p2):
            l_phase, l_eos, v_phase, v_eos = _identify_phases(b, p1, p2)

            # If one or both of v_phase and l_phase is None, this is not
            # a vapor-liquid pair, so sjip constraint
            if v_eos is None and l_eos is None:
                return Constraint.Skip

            return 1e3 == 1e3*sum(b._mole_frac_pdew[p1, p2, j]
                                  for j in b.params.component_list)
        b.eq_mole_frac_pdew = Constraint(b.params._pe_pairs,
                                         rule=rule_mole_frac_dew_press)
