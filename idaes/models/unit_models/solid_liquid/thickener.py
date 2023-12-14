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
Thickener unit model.

This model extends the SLSeparator unit model by adding constraints that relate
area and vessel height to the liquid recovery fraction.

"""
# Import Python libraries
import logging

# Import Pyomo libraries
from pyomo.environ import Var, units

# Import IDAES cores
from idaes.core import (
    declare_process_block_class,
)
from idaes.models.unit_models.solid_liquid.sl_separator import (
    SLSeparatorData,
)


__author__ = "Andrew Lee"


# Set up logger
logger = logging.getLogger("idaes.unit_model")


@declare_process_block_class("Thickener0D")
class Thickener0DData(SLSeparatorData):
    """
    Thickener0D Unit Model Class
    """

    CONFIG = SLSeparatorData.CONFIG()
    # TODO: Method to calculate pinch settling velocity

    def build(self):
        """
        Begin building model (pre-DAE transformation).

        Args:
            None

        Returns:
            None
        """
        logger.warning(
            "The Thickener0D model is currently a beta capability and will "
            "likely change in the next release as a more predictive version is "
            "developed."
        )
        # Call super().build to setup dynamics
        super().build()

        # Add additional variables and constraints
        s_metadata = self.solid_state.params.get_metadata()
        # TODO: Add support for molar basis

        self.area = Var(
            initialize=1,
            units=s_metadata.get_derived_units("area"),
            doc="Cross sectional area of thickener",
        )

        self.height = Var(
            initialize=1,
            units=s_metadata.get_derived_units("length"),
            doc="Total depth of thickener",
        )

        self.height_clarification = Var(
            initialize=1,
            units=s_metadata.get_derived_units("length"),
            doc="Depth of clarification zone",
        )

        self.settling_velocity_pinch = Var(
            self.flowsheet().time,
            initialize=1,
            units=s_metadata.get_derived_units("velocity"),
            doc="Settling velocity of suspension at pinch point",
        )

        self.liquid_solid_pinch = Var(
            self.flowsheet().time,
            initialize=0.2,
            units=units.dimensionless,
            doc="Liquid-solid ratio (mass basis) at pinch point",
        )

        self.liquid_solid_underflow = Var(
            self.flowsheet().time,
            initialize=0.2,
            units=units.dimensionless,
            doc="Liquid-solid ratio (mass basis) at underflow",
        )

        self.settling_time = Var(
            self.flowsheet().time,
            initialize=1,
            units=s_metadata.get_derived_units("time"),
            doc="Settling time of suspension",
        )

        @self.Constraint(
            self.flowsheet().time, doc="Constraint to calculate L/S ratio at underflow"
        )
        def underflow_sl_constraint(b, t):
            return b.liquid_solid_underflow[t] * b.solid_state[
                t
            ].flow_mass == units.convert(
                b.split.retained_state[t].flow_mass,
                to_units=s_metadata.get_derived_units("flow_mass"),
            )

        @self.Constraint(
            self.flowsheet().time, doc="Constraint to estimate cross-sectional area"
        )
        def cross_sectional_area_constraint(b, t):
            return b.area * b.liquid_inlet_state[
                t
            ].dens_mass * b.settling_velocity_pinch[t] == b.solid_state[t].flow_mass * (
                b.liquid_solid_pinch[t] - b.liquid_solid_underflow[t]
            )

        @self.Constraint(
            self.flowsheet().time, doc="Constraint to estimate height of thickener"
        )
        def height_constraint(b, t):
            return (b.height - b.height_clarification) * b.area * b.solid_state[
                t
            ].dens_mass == units.convert(
                b.settling_time[t]
                * (
                    b.solid_state[t].flow_mass
                    + units.convert(
                        b.solid_state[t].dens_mass / b.liquid_inlet_state[t].dens_mass,
                        to_units=units.dimensionless,
                    )
                    * (
                        b.liquid_inlet_state[t].flow_mass
                        + b.split.retained_state[t].flow_mass
                    )
                    / 2
                ),
                to_units=s_metadata.get_derived_units("mass"),
            )

    def _get_performance_contents(self, time_point=0):

        return {
            "vars": {
                "Area": self.area,
                "Height": self.height,
                "Liquid Recovery": self.liquid_recovery[time_point],
                "Underflow L/S": self.liquid_solid_underflow[time_point],
                "Pinch L/S": self.liquid_solid_pinch[time_point],
                "Critical Settling Velocity": self.settling_velocity_pinch[time_point],
                "Settling Time": self.settling_time[time_point],
            }
        }
