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
from pyomo.environ import Constraint, Param, Var, value
from idaes.core.util.math import smooth_max, smooth_min


def phase_equil(b, phase_pair):
    # This method is called via StateBlock.build, thus does not need clean-up
    # try/except statements
    suffix = "_"+phase_pair[0]+"_"+phase_pair[1]

    # Definition of equilibrium temperature for smooth VLE
    b.add_component("_t1"+suffix, Var(
            initialize=b.temperature.value,
            doc='Intermediate temperature for calculating Teq'))
    _t1 = getattr(b, "_t1"+suffix)

    b.add_component("eps_1"+suffix, Param(default=0.01,
                                          mutable=True,
                                          doc='Smoothing parameter for Teq'))
    b.add_component("eps_2"+suffix, Param(default=0.0005,
                                          mutable=True,
                                          doc='Smoothing parameter for Teq'))
    eps_1 = getattr(b, "eps_1"+suffix)
    eps_2 = getattr(b, "eps_1"+suffix)

    # PSE paper Eqn 13
    def rule_t1(b):
        return _t1 == smooth_max(b.temperature,
                                 b.temperature_bubble,
                                 eps_1)
    b.add_component("_t1_constraint"+suffix, Constraint(rule=rule_t1))

    # PSE paper Eqn 14
    # TODO : Add option for supercritical extension
    def rule_teq(b):
        return b._teq[phase_pair] == smooth_min(_t1,
                                                b.temperature_dew,
                                                eps_2)
    b.add_component("_teq_constraint"+suffix, Constraint(rule=rule_teq))


def phase_equil_initialization(b, phase_pair):
    suffix = "_"+phase_pair[0]+"_"+phase_pair[1]

    for c in b.component_objects(Constraint):
        # Activate equilibrium constraints
        if c.local_name in ("_t1_constraint"+suffix,
                            "_teq_constraint"+suffix):
            c.activate()


def calculate_teq(b, phase_pair):
    suffix = "_"+phase_pair[0]+"_"+phase_pair[1]
    _t1 = getattr(b, "_t1"+suffix)

    _t1.value = max(value(b.temperature), b.temperature_bubble.value)
    b._teq[phase_pair].value = min(_t1.value, b.temperature_dew.value)
