#################################################################################
# The Institute for the Design of Advanced Energy Systems Integrated Platform
# Framework (IDAES IP) was produced under the DOE Institute for the
# Design of Advanced Energy Systems (IDAES).
#
# Copyright (c) 2018-2024 by the software owners: The Regents of the
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

from pyomo.environ import (
    Constraint,
    Expression,
    Param,
    Set,
    units as pyunits,
    Var,
    value,
)

from idaes.core.util.exceptions import ConfigurationError
from idaes.core.util.math import smooth_min
from idaes.models.properties.modular_properties.base.utility import (
    identify_VL_component_list,
)
from idaes.models.properties.modular_properties.base.utility import (
    estimate_Tbub,
    estimate_Tdew,
)
from idaes.models.properties.modular_properties.eos.ceos import (
    calculate_equilibrium_cubic_coefficients,
)
import idaes.core.util.scaling as iscale
from idaes.core.scaling import CustomScalerBase


# Small value for initializing slack variables
EPS_INIT = 1e-4


class CubicComplementarityVLEScaler(CustomScalerBase):
    """
    Scaler for CubicComplementarityVLE
    """

    def variable_scaling_routine(self, model, phase_pair, overwrite: bool = False):
        suffix = "_" + phase_pair[0] + "_" + phase_pair[1]
        sf_T = self.get_scaling_factor(model.temperature)
        if model.is_property_constructed("_teq_constraint" + suffix):
            teq = model._teq[phase_pair]  # pylint: disable=protected-access
            self.set_component_scaling_factor(teq, sf_T, overwrite=overwrite)

    def constraint_scaling_routine(self, model, phase_pair, overwrite: bool = False):
        suffix = "_" + phase_pair[0] + "_" + phase_pair[1]
        if model.is_property_constructed("_teq_constraint" + suffix):
            teq = model._teq[phase_pair]  # pylint: disable=protected-access
            sf_T = self.get_scaling_factor(teq)
            teq_cons = getattr(model, "_teq_constraint" + suffix)
            self.set_component_scaling_factor(teq_cons, sf_T, overwrite=overwrite)


# -----------------------------------------------------------------------------
class CubicComplementarityVLE:
    """
    Improved Vapor-Liquid Equilibrium complementarity formulation for Cubic Equations of State.
    """

    default_scaler = CubicComplementarityVLEScaler

    @staticmethod
    def phase_equil(b, phase_pair):
        """
        Method for constructing phase equilibrium variables and constraints
        """
        # This method is called via StateBlock.build, thus does not need clean-up
        # try/except statements
        suffix = "_" + phase_pair[0] + "_" + phase_pair[1]

        # Smooth VLE assumes a liquid and a vapor phase, so validate this
        (
            l_phase,
            v_phase,
            _,
            _,
            _,
            _,
        ) = identify_VL_component_list(b, phase_pair)

        if l_phase is None or v_phase is None:
            raise ConfigurationError(
                f"{b.params.name} - Generic Property Package phase pair {phase_pair[0]}-{phase_pair[1]} "
                "was set to use CubicComplementarityVLE formulation, however this is not a vapor-liquid pair."
            )

        try:
            lobj = b.params.get_phase(l_phase)
            ctype = lobj._cubic_type  # pylint: disable=protected-access
            cname = lobj.config.equation_of_state_options["type"].name
            vobj = b.params.get_phase(v_phase)

            if (
                ctype != vobj._cubic_type  # pylint: disable=protected-access
                or lobj.config.equation_of_state_options
                != vobj.config.equation_of_state_options
            ):
                raise ConfigurationError(
                    f"{b.params.name} - CubicComplementarityVLE formulation requires that both "
                    "phases use the same type of cubic equation of state."
                )
        except AttributeError:
            raise ConfigurationError(
                f"{b.params.name} - CubicComplementarityVLE formulation only supports cubic equations of state."
            )

        # Definition of equilibrium temperature for smooth VLE
        uom = b.params.get_metadata().default_units
        t_units = uom.TEMPERATURE
        f_units = uom.AMOUNT / uom.TIME

        vl_phase_set = Set(initialize=[phase_pair[0], phase_pair[1]])
        b.add_component("_vle_set" + suffix, vl_phase_set)

        s = Var(
            vl_phase_set,
            initialize=EPS_INIT,
            bounds=(0, None),
            doc="Slack variable for equilibrium temperature",
            units=t_units,
        )
        b.add_component("s" + suffix, s)

        # Equilibrium temperature
        def rule_teq(b):
            return (
                b._teq[phase_pair]  # pylint: disable=protected-access
                - b.temperature
                - s[v_phase]
                + s[l_phase]
                == 0
            )

        b.add_component("_teq_constraint" + suffix, Constraint(rule=rule_teq))

        # Epsilon variables will be given units of flow, as this is usually what dominates
        # the complementarity equations.
        eps_t = Param(
            default=1,
            mutable=True,
            doc="Smoothing parameter for temperature complementarity",
            units=f_units,
        )
        b.add_component("eps_t" + suffix, eps_t)

        eps_z = Param(
            default=1,
            mutable=True,
            doc="Smoothing parameter for cubic root complementarities",
            units=f_units,
        )
        b.add_component("eps_z" + suffix, eps_z)

        gp = Var(
            vl_phase_set,
            initialize=EPS_INIT,
            bounds=(0, None),
            doc="Slack variable for cubic second derivative for phase p",
            units=pyunits.dimensionless,
        )
        b.add_component("gp" + suffix, gp)

        gn = Var(
            vl_phase_set,
            initialize=EPS_INIT,
            bounds=(0, None),
            doc="Slack variable for cubic second derivative for phase p",
            units=pyunits.dimensionless,
        )
        b.add_component("gn" + suffix, gn)

        def rule_temperature_slack_complementarity(b, p):
            flow_phase = b.flow_mol_phase[p]

            return smooth_min(s[p] * f_units / t_units, flow_phase, eps_t) == 0

        b.add_component(
            "temperature_slack_complementarity" + suffix,
            Constraint(
                vl_phase_set,
                rule=rule_temperature_slack_complementarity,
            ),
        )

        def rule_cubic_second_derivative(b, p):
            p1, p2 = phase_pair

            _b, _, _ = calculate_equilibrium_cubic_coefficients(
                b, cname, ctype, p1, p2, p
            )
            z = b.compress_fact_phase[p]

            return 6 * z + 2 * _b

        b.add_component(
            "cubic_second_derivative" + suffix,
            Expression(vl_phase_set, rule=rule_cubic_second_derivative),
        )

        def rule_cubic_root_complementarity(b, p):
            der = getattr(b, "cubic_second_derivative" + suffix)
            return der[p] == gp[p] - gn[p]

        b.add_component(
            "cubic_root_complementarity" + suffix,
            Constraint(vl_phase_set, rule=rule_cubic_root_complementarity),
        )

        def rule_cubic_slack_complementarity(b, p):
            flow_phase = b.flow_mol_phase[p]
            if b.params.get_phase(p).is_vapor_phase():
                return smooth_min(gn[p] * f_units, flow_phase, eps_z) == 0
            else:
                return smooth_min(gp[p] * f_units, flow_phase, eps_z) == 0

        b.add_component(
            "cubic_slack_complementarity" + suffix,
            Constraint(vl_phase_set, rule=rule_cubic_slack_complementarity),
        )

    @staticmethod
    def calculate_scaling_factors(b, phase_pair):
        """
        Method to calculate scaling factors for phase equilibrium
        """
        suffix = "_" + phase_pair[0] + "_" + phase_pair[1]
        sf_T = iscale.get_scaling_factor(b.temperature, default=1, warning=True)

        if b.is_property_constructed("_teq_constraint" + suffix):
            teq_cons = getattr(b, "_teq_constraint" + suffix)
            # pylint: disable-next=protected-access
            iscale.set_scaling_factor(b._teq[phase_pair], sf_T)
            iscale.constraint_scaling_transform(teq_cons, sf_T, overwrite=False)
        else:
            pass

    @staticmethod
    def calculate_teq(blk, pp):
        """
        Method to calculate initial guess for equilibrium temperature
        """
        # ---------------------------------------------------------------------
        # If present, initialize bubble and dew point calculations, and
        # equilibrium temperature _teq
        T_units = blk.params.get_metadata().default_units.TEMPERATURE

        (
            liquid_phase,
            vapor_phase,
            raoult_comps,
            henry_comps,
            l_only_comps,
            v_only_comps,
        ) = identify_VL_component_list(blk, pp)

        if len(v_only_comps) == 0:
            if blk.is_property_constructed("temperature_bubble"):
                Tbub = value(blk.temperature_bubble[pp])
            else:
                Tbub = estimate_Tbub(
                    blk, T_units, raoult_comps, henry_comps, liquid_phase
                )
            t1 = max(value(blk.temperature), Tbub)
        else:
            t1 = value(blk.temperature)

        if len(l_only_comps) == 0:
            if blk.is_property_constructed("temperature_dew"):
                Tdew = value(blk.temperature_bubble[pp])
            else:
                Tdew = estimate_Tdew(
                    blk, T_units, raoult_comps, henry_comps, liquid_phase
                )
            t2 = min(t1, Tdew)
        else:
            t2 = t1

        blk._teq[pp].set_value(t2)  # pylint: disable=protected-access

        # ---------------------------------------------------------------------
        # Initialize sV and sL slacks
        _calculate_temperature_slacks(blk, pp, liquid_phase, vapor_phase)

        # ---------------------------------------------------------------------
        # If flash, initialize g+ and g- slacks
        _calculate_ceos_derivative_slacks(blk, pp, liquid_phase, vapor_phase)

    @staticmethod
    def phase_equil_initialization(b, phase_pair):
        """
        Method to initialize phase equilibrium
        """
        suffix = "_" + phase_pair[0] + "_" + phase_pair[1]

        for c in b.component_objects(Constraint):
            # Activate equilibrium constraints
            if c.local_name in (
                "_teq_constraint" + suffix,
                "temperature_slack_complementarity" + suffix,
                "cubic_root_complementarity" + suffix,
                "cubic_slack_complementarity" + suffix,
            ):
                c.activate()


def _calculate_temperature_slacks(b, phase_pair, liquid_phase, vapor_phase):
    suffix = "_" + phase_pair[0] + "_" + phase_pair[1]

    s = getattr(b, "s" + suffix)

    # pylint: disable-next=protected-access
    if value(b._teq[phase_pair]) > value(b.temperature):
        # pylint: disable-next=protected-access
        s[vapor_phase].set_value(value(b._teq[phase_pair] - b.temperature))
        s[liquid_phase].set_value(EPS_INIT)
    # pylint: disable-next=protected-access
    elif value(b._teq[phase_pair]) < value(b.temperature):
        s[vapor_phase].set_value(EPS_INIT)
        # pylint: disable-next=protected-access
        s[liquid_phase].set_value(value(b.temperature - b._teq[phase_pair]))
    else:
        s[vapor_phase].set_value(EPS_INIT)
        s[liquid_phase].set_value(EPS_INIT)


def _calculate_ceos_derivative_slacks(b, phase_pair, liquid_phase, vapor_phase):
    suffix = "_" + phase_pair[0] + "_" + phase_pair[1]

    gp = getattr(b, "gp" + suffix)
    gn = getattr(b, "gn" + suffix)
    der = getattr(b, "cubic_second_derivative" + suffix)

    if value(der[liquid_phase]) < 0:
        gp[liquid_phase].set_value(EPS_INIT)
        gn[liquid_phase].set_value(value(-der[liquid_phase]))
    else:
        gp[liquid_phase].set_value(EPS_INIT)
        gn[liquid_phase].set_value(EPS_INIT)

    if value(der[vapor_phase]) > 0:
        gp[vapor_phase].set_value(value(der[vapor_phase]))
        gn[vapor_phase].set_value(EPS_INIT)
    else:
        gp[vapor_phase].set_value(EPS_INIT)
        gn[vapor_phase].set_value(EPS_INIT)
