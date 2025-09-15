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

Burgard, A.P., Eason, J.P., Eslick, J.C., Ghouse, J.H., Lee, A., Biegler, L.T.,
Miller, D.C., 2018, A Smooth, Square Flash Formulation for Equation-Oriented
Flowsheet Optimization. Proceedings of the 13th International Symposium on
Process Systems Engineering – PSE 2018, July 1-5, 2018, San Diego.
"""
# TODO: Pylint complains about variables with _x names as they are built by other classes
# pylint: disable=protected-access

from pyomo.environ import Constraint, Param, Var, value
from idaes.core.util.math import smooth_max, smooth_min
from idaes.models.properties.modular_properties.base.utility import (
    identify_VL_component_list,
)
import idaes.core.util.scaling as iscale
from idaes.core.scaling import CustomScalerBase


class SmoothVLEScaler(CustomScalerBase):
    """
    Scaling method for the SmoothVLE method for phase equilibrium
    """

    def variable_scaling_routine(self, model, phase_pair, overwrite: bool = False):
        suffix = "_" + phase_pair[0] + "_" + phase_pair[1]
        sf_T = self.get_scaling_factor(model.temperature)

        if model.is_property_constructed("_t1" + suffix):
            t1 = getattr(model, "_t1" + suffix)
            self.set_component_scaling_factor(t1, sf_T, overwrite=overwrite)
        # _teq is scaled in main method

    def constraint_scaling_routine(self, model, phase_pair, overwrite: bool = False):
        suffix = "_" + phase_pair[0] + "_" + phase_pair[1]

        if model.is_property_constructed("_t1_constraint" + suffix):
            t1 = getattr(model, "_t1" + suffix)
            t1_con = getattr(model, "_t1_constraint" + suffix)
            self.scale_constraint_by_component(t1_con, t1, overwrite=overwrite)

        _teq_cons = getattr(model, "_teq_constraint" + suffix)
        self.scale_constraint_by_component(
            _teq_cons, model._teq[phase_pair], overwrite=overwrite
        )


# -----------------------------------------------------------------------------
class SmoothVLE(object):
    """Methods for constructing equations associated with Smooth VLE formulation."""

    default_scaler = SmoothVLEScaler

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
            _,
            _,
            _,
            _,
            l_only_comps,
            v_only_comps,
        ) = identify_VL_component_list(b, phase_pair)

        # Definition of equilibrium temperature for smooth VLE
        t_units = b.params.get_metadata().default_units.TEMPERATURE
        if v_only_comps == []:
            b.add_component(
                "_t1" + suffix,
                Var(
                    initialize=b.temperature.value,
                    doc="Intermediate temperature for calculating Teq",
                    units=t_units,
                ),
            )
            _t1 = getattr(b, "_t1" + suffix)

            b.add_component(
                "eps_1" + suffix,
                Param(
                    default=0.01,
                    mutable=True,
                    doc="Smoothing parameter for Teq",
                    units=t_units,
                ),
            )
            eps_1 = getattr(b, "eps_1" + suffix)

            # PSE paper Eqn 13
            def rule_t1(b):
                return _t1 == smooth_max(
                    b.temperature, b.temperature_bubble[phase_pair], eps_1
                )

            b.add_component("_t1_constraint" + suffix, Constraint(rule=rule_t1))
        else:
            _t1 = b.temperature

        if l_only_comps == []:
            b.add_component(
                "eps_2" + suffix,
                Param(
                    default=0.0005,
                    mutable=True,
                    doc="Smoothing parameter for Teq",
                    units=t_units,
                ),
            )
            eps_2 = getattr(b, "eps_2" + suffix)

            # PSE paper Eqn 14
            # TODO : Add option for supercritical extension
            def rule_teq(b):
                return b._teq[phase_pair] == smooth_min(
                    _t1, b.temperature_dew[phase_pair], eps_2
                )

        elif v_only_comps == []:

            def rule_teq(b):
                return b._teq[phase_pair] == _t1

        else:

            def rule_teq(b):
                return b._teq[phase_pair] == b.temperature

        b.add_component("_teq_constraint" + suffix, Constraint(rule=rule_teq))

    @staticmethod
    def calculate_scaling_factors(b, phase_pair):
        """
        Method to calculate scaling factors for phase equilibrium
        """
        suffix = "_" + phase_pair[0] + "_" + phase_pair[1]
        sf_T = iscale.get_scaling_factor(b.temperature, default=1, warning=True)

        try:
            _t1 = getattr(b, "_t1" + suffix)
            _t1_cons = getattr(b, "_t1_constraint" + suffix)
            iscale.set_scaling_factor(_t1, sf_T)
            iscale.constraint_scaling_transform(_t1_cons, sf_T, overwrite=False)
        except AttributeError:
            pass

        _teq_cons = getattr(b, "_teq_constraint" + suffix)
        iscale.constraint_scaling_transform(_teq_cons, sf_T, overwrite=False)

    @staticmethod
    def phase_equil_initialization(b, phase_pair):
        """
        Method to initialize phase equilibrium
        """
        suffix = "_" + phase_pair[0] + "_" + phase_pair[1]

        for c in b.component_objects(Constraint):
            # Activate equilibrium constraints
            if c.local_name in ("_t1_constraint" + suffix, "_teq_constraint" + suffix):
                c.activate()

    @staticmethod
    def calculate_teq(b, phase_pair):
        """
        Method to calculate initial guess for equilibrium temperature
        """
        suffix = "_" + phase_pair[0] + "_" + phase_pair[1]

        if hasattr(b, "eq_temperature_bubble"):
            _t1 = getattr(b, "_t1" + suffix)
            _t1.set_value(
                max(value(b.temperature), value(b.temperature_bubble[phase_pair]))
            )
        else:
            _t1 = b.temperature

        if hasattr(b, "eq_temperature_dew"):
            b._teq[phase_pair].set_value(
                min(value(_t1), value(b.temperature_dew[phase_pair]))
            )
        else:
            b._teq[phase_pair].set_value(value(_t1))
