#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2023 by the software owners: The Regents of the
# University of California, through Lawrence Berkeley National Laboratory,
# National Technology & Engineering Solutions of Sandia, LLC, Carnegie Mellon
# University, West Virginia University Research Corporation, et al.
# All rights reserved.  Please see the files COPYRIGHT.md and LICENSE.md
# for full copyright and license information.
#################################################################################
"""
Implementation of the formulation proposed in:

Dabadghao, V., Ghouse, J., Eslick, J., Lee, A., Burgard, A., Miller, D., Biegler, L., 2022,
A complementarity-based vapor-liquid equilibrium formulation for equation-oriented simulation
and optimization. AIChE Journal, DOI: 10.1002/aic.18029
"""
# TODO: Pylint complains about variables with _x names as they are built by other classes
# pylint: disable=protected-access

from pyomo.environ import Constraint, Param, units as pyunits, Var, value

from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.math import smooth_min
from idaes.models.properties.modular_properties.base.utility import (
    identify_VL_component_list,
)
from idaes.models.properties.modular_properties.base.utility import (
    estimate_Tbub,
    estimate_Tdew,
)
import idaes.core.util.scaling as iscale


# -----------------------------------------------------------------------------
class SmoothVLE2:
    """
    Improved Vapor-Liquid Equilibrium complementarity formulation for Cubic Equations of State.
    """

    @staticmethod
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
        ) = identify_VL_component_list(b, phase_pair)

        if l_phase is None or v_phase is None:
            raise ConfigurationError(
                f"{b.params.name} Generic Property Package phase pair {phase_pair[0]}-{phase_pair[1]} "
                "was set to use Smooth VLE formulation, however this is not a vapor-liquid pair."
            )

        # Definition of equilibrium temperature for smooth VLE
        uom = b.params.get_metadata().default_units
        t_units = uom.TEMPERATURE
        f_units = uom.AMOUNT / uom.TIME

        s = Var(
            b.params.phase_list,
            initialize=0.0,
            bounds=(0, None),
            doc="Slack variable for equilibrium temperature",
            units=pyunits.dimensionless,
        )
        b.add_component("s" + suffix, s)

        # Equilibrium temperature
        def rule_teq(b):
            if b.params.get_phase(phase_pair[0]).is_vapor_phase():
                vapor_phase = phase_pair[0]
                liquid_phase = phase_pair[1]
            else:
                vapor_phase = phase_pair[1]
                liquid_phase = phase_pair[0]
            return (
                b._teq[phase_pair]
                - b.temperature
                - s[vapor_phase] * t_units
                + s[liquid_phase] * t_units
                == 0
            )

        b.add_component("_tbar_constraint" + suffix, Constraint(rule=rule_teq))

        eps = Param(
            default=1e-04,
            mutable=True,
            doc="Smoothing parameter for complementarities",
            units=f_units,
        )
        b.add_component("eps" + suffix, eps)

        gp = Var(
            b.params.phase_list,
            initialize=0.0,
            bounds=(0, None),
            doc="Slack variable for cubic second derivative for phase p",
            units=pyunits.dimensionless,
        )
        b.add_component("gp" + suffix, gp)

        gn = Var(
            b.params.phase_list,
            initialize=0.0,
            bounds=(0, None),
            doc="Slack variable for cubic second derivative for phase p",
            units=pyunits.dimensionless,
        )
        b.add_component("gn" + suffix, gn)

        def rule_temperature_slack_complementarity(b, p):
            flow_phase = b.flow_mol_phase[p]

            return smooth_min(s[p] * f_units, flow_phase, eps) == 0

        b.add_component(
            "temperature_slack_complementarity" + suffix,
            Constraint(
                b.params.phase_list,
                rule=rule_temperature_slack_complementarity,
            ),
        )

        def rule_cubic_root_complementarity(b, p):
            p1, p2 = phase_pair
            pobj = b.params.get_phase(p)
            cname = pobj.config.equation_of_state_options["type"].name
            cubic_second_derivative = getattr(
                b,
                "_" + cname + "_cubic_second_derivative",
            )
            return cubic_second_derivative[p1, p2, p] == gp[p] - gn[p]

        b.add_component(
            "cubic_root_complementarity" + suffix,
            Constraint(b.params.phase_list, rule=rule_cubic_root_complementarity),
        )

        def rule_cubic_slack_complementarity(b, p):
            flow_phase = b.flow_mol_phase[p]
            if b.params.get_phase(p).is_vapor_phase():
                return smooth_min(gn[p] * f_units, flow_phase, eps) == 0
            else:
                return smooth_min(gp[p] * f_units, flow_phase, eps) == 0

        b.add_component(
            "cubic_slack_complementarity" + suffix,
            Constraint(b.params.phase_list, rule=rule_cubic_slack_complementarity),
        )

    @staticmethod
    def calculate_scaling_factors(b, phase_pair):
        suffix = "_" + phase_pair[0] + "_" + phase_pair[1]
        sf_T = iscale.get_scaling_factor(b.temperature, default=1, warning=True)

        try:
            teq_cons = getattr(b, "_teq_constraint" + suffix)
            iscale.set_scaling_factor(b._teq[phase_pair], sf_T)
            iscale.constraint_scaling_transform(teq_cons, sf_T, overwrite=False)
        except AttributeError:
            pass

    @staticmethod
    def calculate_teq(blk, pp):
        # ---------------------------------------------------------------------
        # If present, initialize bubble and dew point calculations, and
        # equilibrium temperature _teq
        T_units = blk.params.get_metadata().default_units.TEMPERATURE

        liquid_phase, _, raoult_comps, henry_comps, _, _ = identify_VL_component_list(
            blk, pp
        )

        Tbub = estimate_Tbub(blk, T_units, raoult_comps, henry_comps, liquid_phase)
        Tdew = estimate_Tbub(blk, T_units, raoult_comps, henry_comps, liquid_phase)

        assert False
        suffix = "_" + phase_pair[0] + "_" + phase_pair[1]

        if hasattr(b, "eq_temperature_bubble"):
            _t1 = getattr(b, "_t1" + suffix)
            _t1.value = max(
                value(b.temperature), b.temperature_bubble[phase_pair].value
            )
        else:
            _t1 = b.temperature

        if hasattr(b, "eq_temperature_dew"):
            b._teq[phase_pair].value = min(
                _t1.value, b.temperature_dew[phase_pair].value
            )
        else:
            b._teq[phase_pair].value = _t1.value

        # ---------------------------------------------------------------------
        # Initialize sV and sL slacks
        _calculate_temperature_slacks(blk, pp)

        # ---------------------------------------------------------------------
        # If flash, initialize g+ and g- slacks
        _calculate_ceos_derivative_slacks(blk, pp)

    @staticmethod
    def phase_equil_initialization(b, phase_pair):
        suffix = "_" + phase_pair[0] + "_" + phase_pair[1]

        for c in b.component_objects(Constraint):
            # Activate equilibrium constraints
            if c.local_name in ("_teq_constraint" + suffix,):
                c.deactivate()


def _calculate_temperature_slacks(b, phase_pair):
    suffix = "_" + phase_pair[0] + "_" + phase_pair[1]

    s = getattr(b, "s" + suffix)

    if b.params.get_phase(phase_pair[0]).is_vapor_phase():
        vapor_phase = phase_pair[0]
        liquid_phase = phase_pair[1]
    else:
        vapor_phase = phase_pair[1]
        liquid_phase = phase_pair[0]

    if b._teq[phase_pair].value > b.temperature.value:
        s[vapor_phase].value = value(b._teq[phase_pair] - b.temperature)
        s[liquid_phase].value = 0
    elif b._teq[phase_pair].value < b.temperature.value:
        s[vapor_phase].value = 0
        s[liquid_phase].value = value(b.temperature - b._teq[phase_pair])
    else:
        s[vapor_phase].value = 0
        s[liquid_phase].value = 0


def _calculate_ceos_derivative_slacks(b, phase_pair):
    suffix = "_" + phase_pair[0] + "_" + phase_pair[1]
    p1, p2 = phase_pair

    gp = getattr(b, "gp" + suffix)
    gn = getattr(b, "gn" + suffix)

    if b.params.get_phase(phase_pair[0]).is_vapor_phase():
        vapor_phase = phase_pair[0]
        liquid_phase = phase_pair[1]
    else:
        vapor_phase = phase_pair[1]
        liquid_phase = phase_pair[0]

    vapobj = b.params.get_phase(vapor_phase)
    liqobj = b.params.get_phase(liquid_phase)
    cname_vap = vapobj.config.equation_of_state_options["type"].name
    cname_liq = liqobj.config.equation_of_state_options["type"].name
    cubic_second_derivative_vap = getattr(
        b, "_" + cname_vap + "_cubic_second_derivative"
    )
    cubic_second_derivative_liq = getattr(
        b, "_" + cname_liq + "_cubic_second_derivative"
    )

    if value(cubic_second_derivative_liq[p1, p2, liquid_phase]) < 0:
        gp[liquid_phase].value = 0
        gn[liquid_phase].value = value(
            -cubic_second_derivative_liq[p1, p2, liquid_phase]
        )
    if value(cubic_second_derivative_vap[p1, p2, vapor_phase]) > 0:
        gp[vapor_phase].value = value(cubic_second_derivative_vap[p1, p2, vapor_phase])
        gn[vapor_phase].value = 0
