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
from pyomo.environ import Constraint


# -----------------------------------------------------------------------------
# Bubble temperture methods
def bubble_temp_ideal(b):
    def rule_bubble_temp(b):
        return sum(b.mole_frac[j] *
                   b._params.config.pressure_sat(b, b.temperature_bubble, j)
                   for j in b._params.component_list) - \
            b.pressure == 0
    b.eq_temperature_bubble = Constraint(rule=rule_bubble_temp)

    def rule_mole_frac_bubble_temp(b, j):
        return b._mole_frac_tbub[j]*b.pressure == \
            b.mole_frac[j]*b._params.config.pressure_sat(b,
                                                         b.temperature_bubble,
                                                         j)
    b.eq_mole_frac_tbub = Constraint(b._params.component_list,
                                     rule=rule_mole_frac_bubble_temp)


#def bubble_temp_log(b):
#    b._sum_mole_frac_tbub = Constraint(
#            expr=1e3 == 1e3*sum(b._mole_frac_tbub[j]
#                                for j in b._params.component_list))
#
#    def rule_bubble_temp(b, j):
#        return b._params.config.log_fug_liq_tbub(j) == b._params.config.log_fug_vap_tbub(j)
#    b.eq_temperature_bubble = Constraint(b._params.component_list,
#                                         rule=rule_bubble_temp)


# -----------------------------------------------------------------------------
# Dew temperture methods
def dew_temp_ideal(b):
    def rule_dew_temp(b):
        return (b.pressure*sum(b.mole_frac[j] /
                               b._params.config.pressure_sat(b,
                                                             b.temperature_dew,
                                                             j)
                               for j in b._params.component_list) - 1 ==
                0)
    b.eq_temperature_dew = Constraint(rule=rule_dew_temp)

    def rule_mole_frac_dew_temp(b, j):
        return (b._mole_frac_tdew[j] *
                b._params.config.pressure_sat(b, b.temperature_dew, j) ==
                b.mole_frac[j]*b.pressure)
    b.eq_mole_frac_tdew = Constraint(b._params.component_list,
                                     rule=rule_mole_frac_dew_temp)


# -----------------------------------------------------------------------------
# Bubble pressure methods
def bubble_press_ideal(b):
    def rule_bubble_press(b):
        return b.pressure_bubble == sum(
                b.mole_frac[j] *
                b._params.config.pressure_sat(b, b.temperature, j)
                for j in b._params.component_list)
    b.eq_pressure_bubble = Constraint(rule=rule_bubble_press)

    def rule_mole_frac_bubble_press(b, j):
        return b._mole_frac_pbub[j]*b.pressure_bubble == \
            b.mole_frac[j]*b._params.config.pressure_sat(b,
                                                         b.temperature,
                                                         j)
    b.eq_mole_frac_pbub = Constraint(b._params.component_list,
                                     rule=rule_mole_frac_bubble_press)


# -----------------------------------------------------------------------------
# Dew temperture methods
def dew_press_ideal(b):
    def rule_dew_press(b):
        return 0 == 1 - b.pressure_dew*sum(
                            b.mole_frac[j] /
                            b._params.config.pressure_sat(b, b.temperature, j)
                            for j in b._params.component_list)
    b.eq_pressure_dew = Constraint(rule=rule_dew_press)

    def rule_mole_frac_dew_press(b, j):
        return (b._mole_frac_pdew[j] *
                b._params.config.pressure_sat(b, b.temperature, j) ==
                b.mole_frac[j]*b.pressure_dew)
    b.eq_mole_frac_pdew = Constraint(b._params.component_list,
                                     rule=rule_mole_frac_dew_press)
