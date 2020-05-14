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


class IdealBubbleDew():
    # -------------------------------------------------------------------------
    # Bubble temperature methods
    def temperature_bubble(b):
        try:
            def rule_bubble_temp(b, p1, p2):
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
            return (b._mole_frac_pdew[p1, p2, j] *
                    get_method(b, "pressure_sat_comp", j)(
                        b, cobj(b, j), b.temperature) ==
                    b.mole_frac_comp[j]*b.pressure_dew[p1, p2])
        b.eq_mole_frac_pdew = Constraint(b.params._pe_pairs,
                                         b.params.component_list,
                                         rule=rule_mole_frac_dew_press)
