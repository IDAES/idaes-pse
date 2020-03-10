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
"""
Implementation of the formulation proposed in:

Burgard, A.P., Eason, J.P., Eslick, J.C., Ghouse, J.H., Lee, A., Biegler, L.T.,
Miller, D.C., 2018, A Smooth, Square Flash Formulation for Equation-Oriented
Flowsheet Optimization. Proceedings of the 13th International Symposium on
Process Systems Engineering – PSE 2018, July 1-5, 2018, San Diego.
"""
from pyomo.environ import Constraint, Expression, Param, Var
from idaes.core.util.math import smooth_max, smooth_min


def phase_equil(b):
    # Definition of equilibrium temperature for smooth VLE
    b._teq = Var(initialize=b.temperature.value,
                 doc='Temperature for calculating phase equilibrium')
    b._t1 = Var(initialize=b.temperature.value,
                doc='Intermediate temperature for calculating Teq')

    b.eps_1 = Param(default=0.01,
                    mutable=True,
                    doc='Smoothing parameter for Teq')
    b.eps_2 = Param(default=0.0005,
                    mutable=True,
                    doc='Smoothing parameter for Teq')

    # PSE paper Eqn 13
    def rule_t1(b):
        return b._t1 == smooth_max(b.temperature,
                                   b.temperature_bubble,
                                   b.eps_1)
    b._t1_constraint = Constraint(rule=rule_t1)

    # PSE paper Eqn 14
    # TODO : Add option for supercritical extension
    def rule_teq(b):
        return b._teq == smooth_min(b._t1,
                                    b.temperature_dew,
                                    b.eps_2)
    b._teq_constraint = Constraint(rule=rule_teq)

    def rule_tr_eq(b, i):
        return b._teq / b.params.temperature_crit_comp[i]
    b._tr_eq = Expression(b.params.component_list,
                          rule=rule_tr_eq,
                          doc='Component reduced temperatures [-]')

    def rule_equilibrium(b, j):
        return (b.params.config.equation_of_state["Vap"].fug_phase_comp(
                    b, "Vap", j) ==
                b.params.config.equation_of_state["Liq"].fug_phase_comp(
                    b, "Liq", j))
    b.equilibrium_constraint = \
        Constraint(b.params.component_list, rule=rule_equilibrium)


def phase_equil_initialization(b):
    for c in b.component_objects(Constraint):
        # Activate equilibrium constraints
        if c.local_name in ("total_flow_balance",
                            "component_flow_balances",
                            "equilibrium_constraint",
                            "sum_mole_frac",
                            "_t1_constraint",
                            "_teq_constraint"):
            c.activate()
