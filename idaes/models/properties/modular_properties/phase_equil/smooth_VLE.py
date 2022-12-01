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
        s = Var(b.params.phase_list,
                initialize=0.0,
                bounds=(0, None),
                doc='Slack variable for equilibrium temperature',
                units=t_units)
        b.add_component("s" + suffix, s)

        # Equilibrium temperature
        def rule_tbar(b):
            return b._tbar[phase_pair] - b.temperature - \
                s['Vap'] + s['Liq'] == 0
        b.add_component("_tbar_constraint" + suffix, Constraint(rule=rule_tbar))

    else:
        b._tbar[phase_pair] = b.temperature

    eps = Param(default=1e-04,
                mutable=True,
                doc='Smoothing parameter for complementarities',
                units=t_units)
    b.add_component("eps" + suffix, eps)
    
    gp = Var(b.params.phase_list,
             initialize=0.0,
             bounds=(0, None),
             doc='Slack variable for cubic second derivative for phase p',
             units=pyunits.dimensionless)
    b.add_component("gp" + suffix, gp)
    
    gn = Var(b.params.phase_list,
             initialize=0.0,
             bounds=(0, None),
             doc='Slack variable for cubic second derivative for phase p',
             units=pyunits.dimensionless)
    b.add_component("gn" + suffix, gn)

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
        _pp = Var(initialize=0.0,
                 bounds=(0, None),
                 doc="Artificial pressure variable to check if P > Pc",
                 units=pyunits.dimensionless)
        b.add_component("_pp" + suffix, _pp)
        
        _pn = Var(initialize=0.0,
                 bounds=(0, None),
                 doc="Artificial pressure variable to check if P < Pc",
                 units=pyunits.dimensionless)
        b.add_component("_pn" + suffix, _pn)

        def rule_pbar(b):
            """
            _pbar = P - Pc * Pm - eps / 4
            """
            return b._pbar[phase_pair] == b.pressure - _pp - eps / 4
        b.add_component("_pbar_constraint" + suffix, Constraint(rule=rule_pbar))

        # Rule P+
        def rule_pressure_comparison_1(b):
            """
            P+ = max(0, P - Pc)
            """
            return _pp == smooth_max(0, b.pressure - b.pressure_crit_mix, eps)
        b.add_component(
            "eq_pressure_comparison_1" + suffix,
            Constraint(rule=rule_pressure_comparison_1),
        )
        
        # Rule P-
        def rule_pressure_comparison_2(b):
            """
            P- = max(0, Pc - P)
            """
            return _pn == smooth_max(0, b.pressure_crit_mix - b.pressure, eps)
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
            return smooth_min(_pp, b.flow_mol_phase['Vap'], eps) == 0
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


def calculate_pbar(b, phase_pair):
    if hasattr(b, "_pbar"):
        suffix = "_" + phase_pair[0] + "_" + phase_pair[1]
        
        _pp = getattr(b, "_pp" + suffix)
        _pn = getattr(b, "_pn" + suffix)
        eps = getattr(b, "eps" + suffix)
        
        _pp.value = max(0, value(b.pressure - b.pressure_crit_mix))
        _pn.value = max(0, value(b.pressure_crit_mix - b.pressure))
        b._pbar[phase_pair].value = value(b.pressure - _pp - eps / 4)


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


# -----------------------------------------------------------------------------
class SmoothVLE(object):
    # Deprecation: Eventually replace static methods with class methods
    phase_equil = phase_equil
    phase_equil_initialization = phase_equil_initialization
    calculate_pbar = calculate_pbar
    calculate_temperature_slacks = calculate_temperature_slacks
    calculate_ceos_derivative_slacks = calculate_ceos_derivative_slacks
    calculate_scaling_factors = calculate_scaling_factors
