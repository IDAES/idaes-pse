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
Modular methods for calculating bubble and dew points

Authors: Andrew Lee, Douglas Allan
"""
# TODO: Pylint complains about variables with _x names as they are built by other classes
# pylint: disable=protected-access

from pyomo.environ import Constraint

from idaes.models.properties.modular_properties.base.utility import (
    get_method,
    get_component_object as cobj,
    identify_VL_component_list,
)
import idaes.core.util.scaling as iscale
from idaes.core.util.exceptions import ConfigurationError
from idaes.core.scaling import (
    ConstraintScalingScheme,
    CustomScalerBase,
)


class IdealBubbleDewScaler(CustomScalerBase):
    """
    Scaling method for the IdealBubbleDew method for bubble/dew point calculations.
    No new variables are created, so only constraints need to be scaled.
    """

    def variable_scaling_routine(self, model, overwrite: bool = False):
        pass

    def constraint_scaling_routine(self, model, overwrite: bool = False):
        sf_P = self.get_scaling_factor(model.pressure)
        sf_mf = {}
        for i, v in model.mole_frac_comp.items():
            sf_mf[i] = self.get_scaling_factor(v)

        for pp in model.params._pe_pairs:
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,  # pylint: disable=W0612
                l_only_comps,
                v_only_comps,
            ) = identify_VL_component_list(model, pp)
            if l_phase is None or v_phase is None:
                continue

            if len(v_only_comps) == 0 and model.is_property_constructed(
                "eq_temperature_bubble"
            ):
                self.set_component_scaling_factor(
                    model.eq_temperature_bubble[pp[0], pp[1]], sf_P, overwrite=overwrite
                )
                for j in model.component_list:
                    # TODO Henry
                    if j in vl_comps:
                        sf = sf_P * sf_mf[j]
                    else:
                        sf = sf_mf[j]
                    self.set_component_scaling_factor(
                        model.eq_mole_frac_tbub[pp[0], pp[1], j],
                        sf,
                        overwrite=overwrite,
                    )

            if len(l_only_comps) == 0 and model.is_property_constructed(
                "eq_temperature_dew"
            ):
                # eq_temperature_dew is well-scaled by default
                self.set_component_scaling_factor(
                    model.eq_temperature_dew[pp[0], pp[1]], 1, overwrite=overwrite
                )
                for j in model.component_list:
                    # TODO Henry
                    if j in vl_comps:
                        sf = sf_P * sf_mf[j]
                    else:
                        sf = sf_mf[j]
                    self.set_component_scaling_factor(
                        model.eq_mole_frac_tdew[pp[0], pp[1], j],
                        sf,
                        overwrite=overwrite,
                    )

            if len(v_only_comps) == 0 and model.is_property_constructed(
                "eq_pressure_bubble"
            ):
                self.set_component_scaling_factor(
                    model.eq_pressure_bubble[pp[0], pp[1]], sf_P, overwrite=overwrite
                )
                for j in model.component_list:
                    # TODO Henry
                    if j in vl_comps:
                        sf = sf_P * sf_mf[j]
                    else:
                        sf = sf_mf[j]
                    self.set_component_scaling_factor(
                        model.eq_mole_frac_pbub[pp[0], pp[1], j],
                        sf,
                        overwrite=overwrite,
                    )
            if len(l_only_comps) == 0 and model.is_property_constructed(
                "eq_pressure_dew"
            ):
                # eq_pressure_dew is well-scaled by default
                self.set_component_scaling_factor(
                    model.eq_pressure_dew[pp[0], pp[1]], 1, overwrite=overwrite
                )
                for j in model.component_list:
                    # TODO Henry
                    if j in vl_comps:
                        sf = sf_P * sf_mf[j]
                    else:
                        sf = sf_mf[j]
                    self.set_component_scaling_factor(
                        model.eq_mole_frac_pdew[pp[0], pp[1], j],
                        sf,
                        overwrite=overwrite,
                    )


class IdealBubbleDew:
    """Bubble and dew point calculations for ideal systems."""

    default_scaler = IdealBubbleDewScaler

    # -------------------------------------------------------------------------
    # Bubble temperature methods
    # This approach can only be used when both liquid and vapor phases use
    # Ideal properties
    # Henry's Law components also cause issues due to the need to (potentially)
    # calculate concentrations at the bubble and dew points
    @staticmethod
    def temperature_bubble(b):
        """
        Rule for calculating bubble temperature
        """
        _non_vle_phase_check(b)
        try:

            def rule_bubble_temp(b, p1, p2):
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    _,
                    v_only_comps,
                ) = identify_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif v_only_comps != []:
                    # Non-condensables present, no bubble point
                    return Constraint.Skip

                return (
                    sum(
                        b.mole_frac_comp[j]
                        * get_method(b, "pressure_sat_comp", j)(
                            b, cobj(b, j), b.temperature_bubble[p1, p2]
                        )
                        for j in vl_comps
                    )
                    + sum(
                        b.mole_frac_comp[j]
                        * b.params.get_component(j)
                        .config.henry_component[l_phase]["method"]
                        .return_expression(b, l_phase, j, b.temperature_bubble[p1, p2])
                        for j in henry_comps
                    )
                    - b.pressure
                ) == 0

            b.eq_temperature_bubble = Constraint(
                b.params._pe_pairs, rule=rule_bubble_temp
            )
        except AttributeError:
            b.del_component(b.eq_temperature_bubble)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_bubble_temp(b, p1, p2, j):
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                _,
                v_only_comps,
            ) = identify_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif v_only_comps != []:
                # Non-condensables present, no bubble point
                return Constraint.Skip

            if j in vl_comps:
                return b._mole_frac_tbub[p1, p2, j] * b.pressure == (
                    b.mole_frac_comp[j]
                    * get_method(b, "pressure_sat_comp", j)(
                        b, cobj(b, j), b.temperature_bubble[p1, p2]
                    )
                )
            elif j in henry_comps:
                return b._mole_frac_tbub[p1, p2, j] * b.pressure == (
                    b.mole_frac_comp[j]
                    * b.params.get_component(j)
                    .config.henry_component[l_phase]["method"]
                    .return_expression(b, l_phase, j, b.temperature_bubble[p1, p2])
                )
            else:
                return b._mole_frac_tbub[p1, p2, j] == 0

        b.eq_mole_frac_tbub = Constraint(
            b.params._pe_pairs, b.component_list, rule=rule_mole_frac_bubble_temp
        )

    @staticmethod
    def scale_temperature_bubble(b, overwrite=True):
        """
        Scaling method for bubble temperature
        """
        sf_P = iscale.get_scaling_factor(b.pressure, default=1e-5, warning=True)
        sf_mf = iscale.min_scaling_factor(
            b.mole_frac_comp.values(), default=1e3, warning=True
        )

        for pp in b.params._pe_pairs:
            for j in b.component_list:
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    _,
                    _,
                    v_only_comps,
                ) = identify_VL_component_list(b, pp)
                if l_phase is None or v_phase is None:
                    continue
                elif v_only_comps != []:
                    continue

                if j in vl_comps:
                    sf = sf_P * sf_mf
                else:
                    sf = sf_mf

                iscale.constraint_scaling_transform(
                    b.eq_temperature_bubble[pp[0], pp[1]], sf_P, overwrite=overwrite
                )
                iscale.constraint_scaling_transform(
                    b.eq_mole_frac_tbub[pp[0], pp[1], j], sf, overwrite=overwrite
                )

    # -------------------------------------------------------------------------
    # Dew temperature methods
    @staticmethod
    def temperature_dew(b):
        """
        Method for calculating dew temperature
        """
        _non_vle_phase_check(b)
        try:

            def rule_dew_temp(b, p1, p2):
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    l_only_comps,
                    _,
                ) = identify_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif l_only_comps != []:
                    # Non-volatiles present, no dew point
                    return Constraint.Skip

                return (
                    b.pressure
                    * (
                        sum(
                            b.mole_frac_comp[j]
                            / get_method(b, "pressure_sat_comp", j)(
                                b, cobj(b, j), b.temperature_dew[p1, p2]
                            )
                            for j in vl_comps
                        )
                        + sum(
                            b.mole_frac_comp[j]
                            / b.params.get_component(j)
                            .config.henry_component[l_phase]["method"]
                            .return_expression(b, l_phase, j, b.temperature_dew[p1, p2])
                            for j in henry_comps
                        )
                    )
                    - 1
                    == 0
                )

            b.eq_temperature_dew = Constraint(b.params._pe_pairs, rule=rule_dew_temp)
        except AttributeError:
            b.del_component(b.eq_temperature_dew)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_dew_temp(b, p1, p2, j):
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                l_only_comps,
                _,
            ) = identify_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif l_only_comps != []:
                # Non-volatiles present, no dew point
                return Constraint.Skip

            if j in vl_comps:
                return (
                    b._mole_frac_tdew[p1, p2, j]
                    * get_method(b, "pressure_sat_comp", j)(
                        b, cobj(b, j), b.temperature_dew[p1, p2]
                    )
                    == b.mole_frac_comp[j] * b.pressure
                )
            elif j in henry_comps:
                return (
                    b._mole_frac_tdew[p1, p2, j]
                    * b.params.get_component(j)
                    .config.henry_component[l_phase]["method"]
                    .return_expression(b, l_phase, j, b.temperature_dew[p1, p2])
                    == b.mole_frac_comp[j] * b.pressure
                )
            else:
                return b._mole_frac_tdew[p1, p2, j] == 0

        b.eq_mole_frac_tdew = Constraint(
            b.params._pe_pairs, b.component_list, rule=rule_mole_frac_dew_temp
        )

    @staticmethod
    def scale_temperature_dew(b, overwrite=True):
        """
        Scaling method for dew temperature
        """
        sf_P = iscale.get_scaling_factor(b.pressure, default=1e-5, warning=True)
        sf_mf = iscale.min_scaling_factor(
            b.mole_frac_comp.values(), default=1e3, warning=True
        )

        for pp in b.params._pe_pairs:
            for j in b.component_list:
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    _,
                    l_only_comps,
                    _,
                ) = identify_VL_component_list(b, pp)
                if l_phase is None or v_phase is None:
                    continue
                elif l_only_comps != []:
                    continue

                if j in vl_comps:
                    sf = sf_P * sf_mf
                else:
                    sf = sf_mf

                # b.eq_temperature_dew is well-scaled by default
                iscale.constraint_scaling_transform(
                    b.eq_mole_frac_tdew[pp[0], pp[1], j], sf, overwrite=overwrite
                )

    # -------------------------------------------------------------------------
    # Bubble pressure methods
    @staticmethod
    def pressure_bubble(b):
        """
        Method for calculating bubble pressure
        """
        _non_vle_phase_check(b)
        try:

            def rule_bubble_press(b, p1, p2):
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    _,
                    v_only_comps,
                ) = identify_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif v_only_comps != []:
                    # Non-condensables present, no bubble point
                    return Constraint.Skip

                return b.pressure_bubble[p1, p2] == (
                    sum(b.mole_frac_comp[j] * b.pressure_sat_comp[j] for j in vl_comps)
                    + sum(
                        b.mole_frac_comp[j] * b.henry[l_phase, j] for j in henry_comps
                    )
                )

            b.eq_pressure_bubble = Constraint(
                b.params._pe_pairs, rule=rule_bubble_press
            )
        except AttributeError:
            b.del_component(b.eq_pressure_bubble)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_bubble_press(b, p1, p2, j):
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                _,
                v_only_comps,
            ) = identify_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif v_only_comps != []:
                # Non-condensables present, no bubble point
                return Constraint.Skip

            if j in vl_comps:
                return (
                    b._mole_frac_pbub[p1, p2, j] * b.pressure_bubble[p1, p2]
                    == b.mole_frac_comp[j] * b.pressure_sat_comp[j]
                )
            if j in henry_comps:
                return (
                    b._mole_frac_pbub[p1, p2, j] * b.pressure_bubble[p1, p2]
                    == b.mole_frac_comp[j] * b.henry[l_phase, j]
                )
            else:
                return b._mole_frac_pbub[p1, p2, j] == 0

        b.eq_mole_frac_pbub = Constraint(
            b.params._pe_pairs, b.component_list, rule=rule_mole_frac_bubble_press
        )

    @staticmethod
    def scale_pressure_bubble(b, overwrite=True):
        """
        Scaling method for bubble pressure
        """
        sf_P = iscale.get_scaling_factor(b.pressure, default=1e-5, warning=True)
        sf_mf = iscale.get_scaling_factor(b.mole_frac_comp, default=1e3, warning=True)

        for pp in b.params._pe_pairs:
            for j in b.component_list:
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    _,
                    _,
                    v_only_comps,
                ) = identify_VL_component_list(b, pp)
                if l_phase is None or v_phase is None:
                    continue
                elif v_only_comps != []:
                    continue

                if j in vl_comps:
                    sf = sf_P * sf_mf
                else:
                    sf = sf_mf

                iscale.constraint_scaling_transform(
                    b.eq_pressure_bubble[pp[0], pp[1]], sf_P, overwrite=overwrite
                )
                iscale.constraint_scaling_transform(
                    b.eq_mole_frac_pbub[pp[0], pp[1], j], sf, overwrite=overwrite
                )

    # -------------------------------------------------------------------------
    # Dew pressure methods
    @staticmethod
    def pressure_dew(b):
        """
        Method for calculating dew pressure
        """
        _non_vle_phase_check(b)
        try:

            def rule_dew_press(b, p1, p2):
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    l_only_comps,
                    _,
                ) = identify_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif l_only_comps != []:
                    # Non-volatiles present, no dew point
                    return Constraint.Skip

                return 0 == 1 - b.pressure_dew[p1, p2] * (
                    sum(b.mole_frac_comp[j] / b.pressure_sat_comp[j] for j in vl_comps)
                    + sum(
                        b.mole_frac_comp[j] / b.henry[l_phase, j] for j in henry_comps
                    )
                )

            b.eq_pressure_dew = Constraint(b.params._pe_pairs, rule=rule_dew_press)
        except AttributeError:
            b.del_component(b.eq_pressure_dew)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_dew_press(b, p1, p2, j):
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                l_only_comps,
                _,
            ) = identify_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif l_only_comps != []:
                # Non-volatiles present, no dew point
                return Constraint.Skip

            if j in vl_comps:
                return (
                    b._mole_frac_pdew[p1, p2, j] * b.pressure_sat_comp[j]
                    == b.mole_frac_comp[j] * b.pressure_dew[p1, p2]
                )
            elif j in henry_comps:
                return (
                    b._mole_frac_pdew[p1, p2, j] * b.henry[l_phase, j]
                    == b.mole_frac_comp[j] * b.pressure_dew[p1, p2]
                )
            else:
                return b._mole_frac_pdew[p1, p2, j] == 0

        b.eq_mole_frac_pdew = Constraint(
            b.params._pe_pairs, b.component_list, rule=rule_mole_frac_dew_press
        )

    @staticmethod
    def scale_pressure_dew(b, overwrite=True):
        """
        Scaling method for dew pressure
        """
        sf_P = iscale.get_scaling_factor(b.pressure, default=1e-5, warning=True)
        sf_mf = iscale.get_scaling_factor(b.mole_frac_comp, default=1e3, warning=True)

        for pp in b.params._pe_pairs:
            for j in b.component_list:
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    _,
                    _,
                    v_only_comps,
                ) = identify_VL_component_list(b, pp)
                if l_phase is None or v_phase is None:
                    continue
                elif v_only_comps != []:
                    continue

                if j in vl_comps:
                    sf = sf_P * sf_mf
                else:
                    sf = sf_mf

                # b.eq_pressure_dew is well-scaled by default
                iscale.constraint_scaling_transform(
                    b.eq_mole_frac_pdew[pp[0], pp[1], j], sf, overwrite=overwrite
                )


class LogBubbleDewScaler(CustomScalerBase):
    """
    Scaling method for the LogBubbleDew scaler
    """

    def variable_scaling_routine(self, model, overwrite: bool = False):
        pass

    def constraint_scaling_routine(self, model, overwrite: bool = False):
        for pp in model.params._pe_pairs:
            (
                l_phase,
                v_phase,
                _,
                _,
                l_only_comps,
                v_only_comps,
            ) = identify_VL_component_list(model, pp)
            if l_phase is None or v_phase is None:
                continue
            if len(v_only_comps) == 0 and model.is_property_constructed(
                "eq_temperature_bubble"
            ):
                self.scale_constraint_by_nominal_value(
                    model.eq_mole_frac_tbub[pp[0], pp[1]],
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )
                for j in model.component_list:
                    # Either a log-form constraint or setting
                    # a mole fraction that is used in no other
                    # equation to zero.
                    self.set_component_scaling_factor(
                        model.eq_temperature_bubble[pp[0], pp[1], j],
                        1,
                        overwrite=overwrite,
                    )
            if len(l_only_comps) == 0 and model.is_property_constructed(
                "eq_temperature_dew"
            ):
                self.scale_constraint_by_nominal_value(
                    model.eq_mole_frac_tdew[pp[0], pp[1]],
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )
                for j in model.component_list:
                    # Either a log-form constraint or setting
                    # a mole fraction that is used in no other
                    # equation to zero.
                    self.set_component_scaling_factor(
                        model.eq_temperature_dew[pp[0], pp[1], j],
                        1,
                        overwrite=overwrite,
                    )
            if len(v_only_comps) == 0 and model.is_property_constructed(
                "eq_pressure_bubble"
            ):
                self.scale_constraint_by_nominal_value(
                    model.eq_mole_frac_pbub[pp[0], pp[1]],
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )
                for j in model.component_list:
                    # Either a log-form constraint or setting
                    # a mole fraction that is used in no other
                    # equation to zero.
                    self.set_component_scaling_factor(
                        model.eq_pressure_bubble[pp[0], pp[1], j],
                        1,
                        overwrite=overwrite,
                    )

            if len(l_only_comps) == 0 and model.is_property_constructed(
                "eq_pressure_dew"
            ):
                self.scale_constraint_by_nominal_value(
                    model.eq_mole_frac_pdew[pp[0], pp[1]],
                    scheme=ConstraintScalingScheme.inverseMaximum,
                    overwrite=overwrite,
                )
                for j in model.component_list:
                    # Either a log-form constraint or setting
                    # a mole fraction that is used in no other
                    # equation to zero.
                    self.set_component_scaling_factor(
                        model.eq_pressure_dew[pp[0], pp[1], j], 1, overwrite=overwrite
                    )


class LogBubbleDew:
    """General bubble and dew point calculations (log formulation)."""

    default_scaler = LogBubbleDewScaler

    # -------------------------------------------------------------------------
    # Bubble temperature methods
    @staticmethod
    def temperature_bubble(b):
        """
        Method for constructing bubble temperature constraint
        """
        try:

            def rule_bubble_temp(b, p1, p2, j):
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    _,
                    v_only_comps,
                ) = identify_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif v_only_comps != []:
                    # Non-condensables present, no bubble point
                    return Constraint.Skip

                l_eos = b.params.get_phase(l_phase).config.equation_of_state
                v_eos = b.params.get_phase(v_phase).config.equation_of_state

                if j in vl_comps or j in henry_comps:
                    return l_eos.log_fug_phase_comp_Tbub(
                        b, l_phase, j, (p1, p2)
                    ) == v_eos.log_fug_phase_comp_Tbub(b, v_phase, j, (p1, p2))
                else:
                    return b._mole_frac_tbub[p1, p2, j] == 0

            b.eq_temperature_bubble = Constraint(
                b.params._pe_pairs, b.component_list, rule=rule_bubble_temp
            )
        except AttributeError:
            b.del_component(b.eq_temperature_bubble)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_bubble_temp(b, p1, p2):
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                _,
                v_only_comps,
            ) = identify_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif v_only_comps != []:
                # Non-condensables present, no bubble point
                return Constraint.Skip

            return 1 == (
                sum(b._mole_frac_tbub[p1, p2, j] for j in vl_comps)
                + sum(b._mole_frac_tbub[p1, p2, j] for j in henry_comps)
            )

        b.eq_mole_frac_tbub = Constraint(
            b.params._pe_pairs, rule=rule_mole_frac_bubble_temp
        )

    @staticmethod
    def scale_temperature_bubble(b, overwrite=True):
        """
        Method for scaling bubble temperature
        """
        sf_mf = iscale.get_scaling_factor(b.mole_frac_comp, default=1e3, warning=True)

        for pp in b.params._pe_pairs:
            (
                l_phase,
                v_phase,
                _,
                _,
                _,
                v_only_comps,
            ) = identify_VL_component_list(b, pp)
            if l_phase is None or v_phase is None:
                continue
            elif v_only_comps != []:
                continue

            # Assume b.eq_temperature_bubble is well-scaled
            iscale.constraint_scaling_transform(
                b.eq_mole_frac_tbub[pp[0], pp[1]], sf_mf, overwrite=overwrite
            )

    # -------------------------------------------------------------------------
    # Dew temperature methods
    @staticmethod
    def temperature_dew(b):
        """
        Method for constructing dew temperature constraint
        """
        try:

            def rule_dew_temp(b, p1, p2, j):
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    l_only_comps,
                    _,
                ) = identify_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif l_only_comps != []:
                    # Non-vapouriszbles present, no dew point
                    return Constraint.Skip

                l_eos = b.params.get_phase(l_phase).config.equation_of_state
                v_eos = b.params.get_phase(v_phase).config.equation_of_state

                if j in vl_comps or j in henry_comps:
                    return l_eos.log_fug_phase_comp_Tdew(
                        b, l_phase, j, (p1, p2)
                    ) == v_eos.log_fug_phase_comp_Tdew(b, v_phase, j, (p1, p2))
                else:
                    return b._mole_frac_tdew[p1, p2, j] == 0

            b.eq_temperature_dew = Constraint(
                b.params._pe_pairs, b.component_list, rule=rule_dew_temp
            )
        except AttributeError:
            b.del_component(b.eq_temperature_dew)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_dew_temp(b, p1, p2):
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                l_only_comps,
                _,
            ) = identify_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif l_only_comps != []:
                # Non-volatiles present, no dew point
                return Constraint.Skip

            return 1 == (
                sum(b._mole_frac_tdew[p1, p2, j] for j in vl_comps)
                + sum(b._mole_frac_tdew[p1, p2, j] for j in henry_comps)
            )

        b.eq_mole_frac_tdew = Constraint(
            b.params._pe_pairs, rule=rule_mole_frac_dew_temp
        )

    @staticmethod
    def scale_temperature_dew(b, overwrite=True):
        """
        Method for scaling dew temperature
        """
        sf_mf = iscale.get_scaling_factor(b.mole_frac_comp, default=1e3, warning=True)

        for pp in b.params._pe_pairs:
            (
                l_phase,
                v_phase,
                _,
                _,
                _,
                v_only_comps,
            ) = identify_VL_component_list(b, pp)
            if l_phase is None or v_phase is None:
                continue
            elif v_only_comps != []:
                continue

            # Assume b.eq_temperature_dew is well-scaled
            iscale.constraint_scaling_transform(
                b.eq_mole_frac_tdew[pp[0], pp[1]], sf_mf, overwrite=overwrite
            )

    # -------------------------------------------------------------------------
    # Bubble pressure methods
    @staticmethod
    def pressure_bubble(b):
        """
        Method for constructing bubble pressure
        """
        try:

            def rule_bubble_press(b, p1, p2, j):
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    _,
                    v_only_comps,
                ) = identify_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif v_only_comps != []:
                    # Non-condensables present, no bubble point
                    return Constraint.Skip

                l_eos = b.params.get_phase(l_phase).config.equation_of_state
                v_eos = b.params.get_phase(v_phase).config.equation_of_state

                if j in vl_comps or j in henry_comps:
                    return l_eos.log_fug_phase_comp_Pbub(
                        b, l_phase, j, (p1, p2)
                    ) == v_eos.log_fug_phase_comp_Pbub(b, v_phase, j, (p1, p2))
                else:
                    return b._mole_frac_pbub[p1, p2, j] == 0

            b.eq_pressure_bubble = Constraint(
                b.params._pe_pairs, b.component_list, rule=rule_bubble_press
            )
        except AttributeError:
            b.del_component(b.eq_pressure_bubble)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_bubble_press(b, p1, p2):
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                _,
                v_only_comps,
            ) = identify_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif v_only_comps != []:
                # Non-condensables present, no bubble point
                return Constraint.Skip

            return 1 == (
                sum(b._mole_frac_pbub[p1, p2, j] for j in vl_comps)
                + sum(b._mole_frac_pbub[p1, p2, j] for j in henry_comps)
            )

        b.eq_mole_frac_pbub = Constraint(
            b.params._pe_pairs, rule=rule_mole_frac_bubble_press
        )

    @staticmethod
    def scale_pressure_bubble(b, overwrite=True):
        """
        Method for scaling bubble pressure
        """
        sf_mf = iscale.get_scaling_factor(b.mole_frac_comp, default=1e3, warning=True)

        for pp in b.params._pe_pairs:
            (
                l_phase,
                v_phase,
                _,
                _,
                _,
                v_only_comps,
            ) = identify_VL_component_list(b, pp)
            if l_phase is None or v_phase is None:
                continue
            elif v_only_comps != []:
                continue

            # Assume b.eq_pressure_bubble is well-scaled
            iscale.constraint_scaling_transform(
                b.eq_mole_frac_pbub[pp[0], pp[1]], sf_mf, overwrite=overwrite
            )

    # -------------------------------------------------------------------------
    # Dew pressure methods
    @staticmethod
    def pressure_dew(b):
        """
        Method constructing dew pressure constraints
        """
        try:

            def rule_dew_press(b, p1, p2, j):
                (
                    l_phase,
                    v_phase,
                    vl_comps,
                    henry_comps,
                    _,
                    v_only_comps,
                ) = identify_VL_component_list(b, (p1, p2))

                if l_phase is None or v_phase is None:
                    # Not a VLE pair
                    return Constraint.Skip
                elif v_only_comps != []:
                    # Non-condensables present, no bubble point
                    return Constraint.Skip

                l_eos = b.params.get_phase(l_phase).config.equation_of_state
                v_eos = b.params.get_phase(v_phase).config.equation_of_state

                if j in vl_comps or j in henry_comps:
                    return l_eos.log_fug_phase_comp_Pdew(
                        b, l_phase, j, (p1, p2)
                    ) == v_eos.log_fug_phase_comp_Pdew(b, v_phase, j, (p1, p2))
                else:
                    return b._mole_frac_pdew[p1, p2, j] == 0

            b.eq_pressure_dew = Constraint(
                b.params._pe_pairs, b.component_list, rule=rule_dew_press
            )
        except AttributeError:
            b.del_component(b.eq_pressure_dew)
            raise

        # Don't need a try/except here, will pass if first constraint did
        def rule_mole_frac_dew_press(b, p1, p2):
            (
                l_phase,
                v_phase,
                vl_comps,
                henry_comps,
                l_only_comps,
                _,
            ) = identify_VL_component_list(b, (p1, p2))

            if l_phase is None or v_phase is None:
                # Not a VLE pair
                return Constraint.Skip
            elif l_only_comps != []:
                # Non-volatiles present, no dew point
                return Constraint.Skip

            return 1 == (
                sum(b._mole_frac_pdew[p1, p2, j] for j in vl_comps)
                + sum(b._mole_frac_pdew[p1, p2, j] for j in henry_comps)
            )

        b.eq_mole_frac_pdew = Constraint(
            b.params._pe_pairs, rule=rule_mole_frac_dew_press
        )

    @staticmethod
    def scale_pressure_dew(b, overwrite=True):
        """
        Method for scaling dew pressure
        """
        sf_mf = iscale.get_scaling_factor(b.mole_frac_comp, default=1e3, warning=True)

        for pp in b.params._pe_pairs:
            (
                l_phase,
                v_phase,
                _,
                _,
                _,
                v_only_comps,
            ) = identify_VL_component_list(b, pp)
            if l_phase is None or v_phase is None:
                continue
            elif v_only_comps != []:
                continue

            # Assume b.eq_pressure_dew is well-scaled
            iscale.constraint_scaling_transform(
                b.eq_mole_frac_pdew[pp[0], pp[1]], sf_mf, overwrite=overwrite
            )


def _non_vle_phase_check(blk):
    if len(blk.phase_list) > 2:
        raise ConfigurationError(
            "Ideal assumption for calculating bubble and/or dew points is only valid "
            "for systems with two phases. Please use LogBubbleDew approach instead."
        )
