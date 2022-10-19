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
"""
Implementation of the formulation proposed in:

Burgard, A.P., Eason, J.P., Eslick, J.C., Ghouse, J.H., Lee, A., Biegler, L.T.,
Miller, D.C., 2018, A Smooth, Square Flash Formulation for Equation-Oriented
Flowsheet Optimization. Proceedings of the 13th International Symposium on
Process Systems Engineering – PSE 2018, July 1-5, 2018, San Diego.
"""
from pyomo.environ import Constraint, Param, Var, value, units as pyunits
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.math import smooth_max, smooth_min
from idaes.models.properties.modular_properties.phase_equil.bubble_dew import (
    _valid_VL_component_list,
)
import idaes.core.util.scaling as iscale


def phase_equil(b, phase_pair):
    # This method is called via StateBlock.build, thus does not need clean-up
    # try/except statements
    suffix = "_" + phase_pair[0] + "_" + phase_pair[1]

    # Smooth VLE assumes a liquid and a vapor phase, so validate this
    (
        l_phase,
        v_phase,
        vl_comps,
        henry_comps,
        l_only_comps,
        v_only_comps,
    ) = _valid_VL_component_list(b, phase_pair)

    if l_phase is None or v_phase is None:
        raise ConfigurationError(
            "{} Generic Property Package phase pair {}-{} was set to use "
            "Smooth VLE formulation, however this is not a vapor-liquid pair.".format(
                b.params.name, phase_pair[0], phase_pair[1]
            )
        )

    # Definition of equilibrium temperature slack variables for smooth VLE
    t_units = b.params.get_metadata().default_units["temperature"]
    if v_only_comps == []:
        b.add_component(
            "s" + suffix,
            Var(b.params.phase_list,
                initialize=0.0,
                bounds=(0, None),
                doc='Slack variable for equilibrium temperature',
                units=t_units),
        )
        s = getattr(b, "s" + suffix)

        # Equilibrium temperature
        def rule_tbar(b):
            return b._tbar[phase_pair] - b.temperature - \
                s['Vap'] + s['Liq'] == 0
        b.add_component("_tbar_constraint" + suffix, Constraint(rule=rule_tbar))

    else:
        b._tbar[phase_pair] = b.temperature

    b.add_component(
        "eps" + suffix,
        Param(default=1e-04,
              mutable=True,
              doc='Smoothing parameter for complementarities',
              units=t_units),
    )
    eps = getattr(b, "eps" + suffix)
    
    b.add_component(
        "gp" + suffix,
        Var(b.params.phase_list,
            initialize=0.0,
            bounds=(0, None),
            doc='Slack variable for cubic second derivative for phase p',
            units=pyunits.dimensionless),
    )
    gp = getattr(b, "gp" + suffix)
    
    b.add_component(
        "gn" + suffix,
        Var(b.params.phase_list,
            initialize=0.0,
            bounds=(0, None),
            doc='Slack variable for cubic second derivative for phase p',
            units=pyunits.dimensionless),
    )
    gn = getattr(b, "gn" + suffix)

    def rule_temperature_slack_complementarity(b, p):
        flow_phase = b.flow_mol_phase[p]
        if b.params.config.supercritical_extension:
            if p == 'Vap':
                return smooth_min(s[p], flow_phase, eps) == 0
            else:
                return Constraint.Skip
        else:
            return smooth_min(s[p], flow_phase, eps) == 0
    b.add_component(
        "temperature_slack_complementarity" + suffix,
        Constraint(
            b.params.phase_list,
            rule=rule_temperature_slack_complementarity,
        ),
    )

    def rule_cubic_root_complementarity(b, p):
        p1, p2 = phase_pair
        return b.cubic_second_derivative[p1, p2, p] == gp[p] - gn[p]
    b.add_component(
        "cubic_root_complementarity" + suffix,
        Constraint(b.params.phase_list, rule=rule_cubic_root_complementarity),
    )

    def rule_cubic_slack_complementarity(b, p):
        flow_phase = b.flow_mol_phase[p]
        if p == 'Vap':
            return smooth_min(gn[p], flow_phase, eps) == 0
        else:
            if b.params.config.supercritical_extension:
                return smooth_min(gp[p] + s[p], flow_phase, eps) == 0
            else:
                return smooth_min(gp[p], flow_phase, eps) == 0
    b.add_component(
        "cubic_slack_complementarity" + suffix,
        Constraint(b.params.phase_list, rule=rule_cubic_slack_complementarity),
    )

    if b.params.config.supercritical_extension:
        b.add_component(
        "pp" + suffix,
        Var(initialize=0.0,
            bounds=(0, None),
            doc="Artificial pressure variable to check if P > Pc",
            units=pyunits.dimensionless)
        )
        pp = getattr(b, "pp" + suffix)
        
        b.add_component(
        "pn" + suffix,
        Var(initialize=0.0,
            bounds=(0, None),
            doc="Artificial pressure variable to check if P < Pc",
            units=pyunits.dimensionless)
        )
        pn = getattr(b, "pn" + suffix)
    
    
        def rule_pbar(b):
            """
            pbar = P - Pc * Pm - eps_4 / 4
            """
            return b.pbar[phase_pair] == b.pressure - pp - eps / 4
        b.add_component("pbar_constraint" + suffix, Constraint(rule=rule_pbar))

        # Rule P+
        def rule_pressure_comparison_1(b):
            """
            P+ = max(0, P - Pc)
            """
            return pp == smooth_max(0, b.pressure - b.pressure_crit_mix, eps)
        b.add_component(
            "eq_pressure_comparison_1" + suffix,
            Constraint(rule=rule_pressure_comparison_1),
        )
        
        # Rule P-
        def rule_pressure_comparison_2(b):
            """
            P- = max(0, Pc - P)
            """
            return pn == smooth_max(0, b.pressure_crit_mix - b.pressure, eps)
        b.add_component(
            "eq_pressure_comparison_2" + suffix,
            Constraint(rule=rule_pressure_comparison_2),
        )
    
        def rule_vapor_flow_complementarity(b):
            """
            0 <= P+ _|_ V >= 0
            or: 
            0 == min(P+, V)
            """
            return smooth_min(pp, b.flow_mol_phase['Vap'], eps) == 0
        b.add_component(
            "vapor_flow_complementarity" + suffix,
            Constraint(rule=rule_vapor_flow_complementarity),
        )


def calculate_scaling_factors(b, phase_pair):
    suffix = "_" + phase_pair[0] + "_" + phase_pair[1]
    sf_T = iscale.get_scaling_factor(b.temperature, default=1, warning=True)

    try:
        tbar_cons = getattr(b, "_tbar_constraint" + suffix)
        iscale.set_scaling_factor(b._tbar[phase_pair], sf_T)
        iscale.constraint_scaling_transform(tbar_cons, sf_T, overwrite=False)
    except AttributeError:
        pass


def phase_equil_initialization(b, phase_pair):
    suffix = "_" + phase_pair[0] + "_" + phase_pair[1]

    for c in b.component_objects(Constraint):
        # Activate equilibrium constraints
        if c.local_name in ("_tbar_constraint" + suffix,):
            c.deactivate()


def calculate_tbar(b, phase_pair):
    suffix = "_" + phase_pair[0] + "_" + phase_pair[1]

    if (hasattr(b, "temperature_bubble") and
        hasattr(b, "temperature_dew")):
        b._tbar[phase_pair].value = \
            0.5 * value(b.temperature_bubble[phase_pair] +
                        b.temperature_dew[phase_pair])

    else:
        b._tbar[phase_pair].value = value(b.temperature_bubble[phase_pair])


def calculate_pbar(b, phase_pair):
    if hasattr(b, "pbar"):
        suffix = "_" + phase_pair[0] + "_" + phase_pair[1]
        
        pp = getattr(b, "pp" + suffix)
        pn = getattr(b, "pn" + suffix)
        eps = getattr(b, "eps" + suffix)
        
        pp.value = max(0, value(b.pressure - b.pressure_crit_mix))
        pn.value = max(0, value(b.pressure_crit_mix - b.pressure))
        b.pbar[phase_pair].value = value(b.pressure - pp - eps / 4)



def fix_tbar(b, phase_pair):
    suffix = "_" + phase_pair[0] + "_" + phase_pair[1]
    b._tbar[phase_pair].fix()


def unfix_tbar(b, phase_pair):
    suffix = "_" + phase_pair[0] + "_" + phase_pair[1]
    b._tbar[phase_pair].unfix()


def calculate_temperature_slacks(b, phase_pair):
    suffix = "_" + phase_pair[0] + "_" + phase_pair[1]
    
    s = getattr(b, "s" + suffix)
    
    if b._tbar[phase_pair].value > b.temperature.value:
        s['Vap'].value = value(b._tbar[phase_pair] - b.temperature)
        s['Liq'].value = 0
    elif b._tbar[phase_pair].value < b.temperature.value:
        s['Vap'].value = 0
        s['Liq'].value = value(b.temperature - b._tbar[phase_pair])
    else:
        s['Vap'].value = 0
        s['Liq'].value = 0


def calculate_ceos_derivative_slacks(b, phase_pair):
    suffix = "_" + phase_pair[0] + "_" + phase_pair[1]
    p1, p2 = phase_pair
    
    gp = getattr(b, "gp" + suffix)
    gn = getattr(b, "gn" + suffix)
    
    if value(b.cubic_second_derivative[p1, p2, 'Liq']) < 0:
        gp['Liq'].value = 0
        gn['Liq'].value = value(-b.cubic_second_derivative[p1, p2, 'Liq'])
    if value(b.cubic_second_derivative[p1, p2, 'Vap']) > 0:
        gp['Vap'].value = value(b.cubic_second_derivative[p1, p2, 'Vap'])
        gn['Vap'].value = 0


def bubble_dew_method(b):
    t_units = b.params.get_metadata().default_units["temperature"]
    p_units = pyunits.Pa
    try:
        b.temperature_bubble = Var(
            b.params._pe_pairs,
            initialize=b.temperature,
            doc="Bubble point temperature of mixture",
            bounds=(b.temperature.lb, b.temperature.ub),
            units=t_units,
        )
        b.temperature_dew = Var(
            b.params._pe_pairs,
            initialize=b.temperature,
            doc="Dew point temperature of mixture",
            bounds=(b.temperature.lb, b.temperature.ub),
            units=t_units,
        )
        
        for pp in b.params._pe_pairs:
            b.temperature_bubble[pp].fix()
            b.temperature_dew[pp].fix()
            
    except AttributeError:
        b.del_component(b.temperature_bubble)
        b.del_component(b.temperature_dew)
        raise


# -----------------------------------------------------------------------------
class SmoothVLE(object):
    # Deprecation: Eventually replace static methods with class methods
    phase_equil = phase_equil
    phase_equil_initialization = phase_equil_initialization
    calculate_tbar = calculate_tbar
    calculate_pbar = calculate_pbar
    fix_tbar = fix_tbar
    unfix_tbar = unfix_tbar
    calculate_temperature_slacks = calculate_temperature_slacks
    calculate_ceos_derivative_slacks = calculate_ceos_derivative_slacks
    calculate_scaling_factors = calculate_scaling_factors
    bubble_dew_method = bubble_dew_method
